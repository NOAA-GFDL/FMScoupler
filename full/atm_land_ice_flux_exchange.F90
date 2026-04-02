!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> \file
!> \brief Performs flux calculations and exchange grid operations for atmosphere, land and ice

module atm_land_ice_flux_exchange_mod

  ! atmos_drivers
  use atmos_model_mod, only: &
       atmos_data_type, & ! derived type containing fields needed for flux exchange between components
       land_ice_atmos_boundary_type ! derived type containing quantities going from land and ice to atmos

  ! MOM6/SIS2
  use ice_model_mod, only: &
       ice_data_type, & ! derived type holding ice model data
       land_ice_boundary_type, & ! derived type for flux exchange between land and sea ice 
       ocean_ice_boundary_type ! derived type for flux exchange between ocean and sea ice
  use ice_model_mod, only: &
       atmos_ice_boundary_type,  ! derived type for flux exchange between atmosphere and sea ice
  use ice_model_mod, only: &
       update_ice_atm_deposition_flux ! updates fluxes of type "air_sea_deposition

  ! Land_lad2
  use land_model_mod, only: &
       land_data_type, & ! derived type to pass information from land to atmosphere
       atmos_land_boundary_type ! derived type to pass information from atmosphere to land

  ! FMSCoupler/shared
  use surface_flux_mod, only: &
       surface_flux, & ! subroutine to compute fluxes on exchange grids
       surface_flux_init ! subroutine to initialize surface_flux_mod

  ! FMSCoupler/full
  use atmos_ocean_fluxes_calc_mod, only: &
       atmos_ocean_fluxes_calc ! subroutine to computes gas fluxes for atmosphere and ocean

  ! FMSCoupler/full
  use atmos_ocean_dep_fluxes_calc_mod, only: &
       atmos_ocean_dep_fluxes_calc ! subroutine to compute ocean and atmosphere deposition gas fluxes

  ! land_lad2
#ifndef _USE_LEGACY_LAND_
  use land_model_mod, only: &
       set_default_diag_filter, & ! subroutine to set default tile diagnostic selector
       register_tiled_diag_field ! subroutine to register diag field within the land model
  use land_model_mod, only: &
       send_tile_data, & ! subroutine to save data in buffer within the land model for the registered field
       dump_tile_diag_fields ! subroutine for workaround with Intel compilers and OpenMP
#endif

  ! am5_phys
#ifndef use_AM3_physics
  use atmos_tracer_driver_mod, only: &
       atmos_tracer_has_surf_setl_flux, & ! function returns True if tracer sedimentation flux > 0 at bottom of atm
       get_atmos_tracer_surf_setl_flux ! subroutine to retrieve tracer sedimentation flux at bottom of atm
  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data_down 
  use atmos_cmip_diag_mod, only: &
       register_cmip_diag_field_2d ! function to register CMIP diagnostic data
  use atmos_global_diag_mod, only: &
       register_global_diag_field, & ! function that calls FMS/register_diag_field for globally averaged data
       get_global_diag_field_id, & ! function to retrieve internally-tracked id of the global diag field
       send_global_diag ! function that calls FMS/diag_manager/send_data for global fields
#ifndef _USE_LEGACY_LAND_
  use land_model_mod, only: &
       send_global_land_diag ! function to save land model field on unstructured grid for global integral
#endif
#endif


#ifdef SCM
  ! option to override various surface boundary conditions for SCM
  use scm_forc_mod,  only: do_specified_flux, scm_surface_flux, &
       do_specified_tskin, TSKIN, &
       do_specified_albedo, ALBEDO_OBS, &
       do_specified_rough_leng, ROUGH_MOM, ROUGH_HEAT,  &
       do_specified_land
#endif

  use FMS
  use FMSconstants, only: &
       rdgas, & ! 287.04, gas constant for dry air [J/kg/deg]
       rvgas, & ! 461.50, gas constant for water vapor
       cp_air, & ! RDGAS/KAPPA , specific heat capacity of dry air at constant pressure [J/kg/deg]
       stefan, & ! 5.6734e-8, Stefan-Boltzmann constant [W/m^2/deg^4]
       WTMAIR, & ! 2.896440e+01,  molecular weight of air [amu]
       HLV, & ! 2.500e6, latent heat of evaporation [J/kg]
       HLF, & ! 3.34e5, latent heat of fusion [J/kg]
       Radius, & ! 6371.0e+3, radius of the Earth [m]
       PI, & ! 3.14159265358979323846
       CP_OCEAN, & ! 3989.24495292815, specific heat capacity [J/kg/deg]
       WTMCO2, & ! 44.00995,  molecular weight of carbon dioxide [amu]
       WTMC, & ! 12.00000,  molecular weight of carbon [amu]
       EPSLN, & ! 1.0e-40, a small number to prevent divide by zero exceptions
       GRAV, & ! 9.80, acceleration due to gravity [m/s^2]
       WTMH2O ! WTMAIR*(RDGAS/RVGAS) molecular weight of water [amu]

  implicit none
  include 'netcdf.inc'
  private

  public :: &
       atm_land_ice_flux_exchange_init, &
       sfc_boundary_layer, &
       generate_sfc_xgrid, &
       flux_down_from_atmos, &
       flux_up_to_atmos, &
       flux_atmos_to_ocean, &
       flux_ex_arrays_dealloc,&
       atm_stock_integrate, &
       send_ice_mask_sic
  
  !> coupler version number that is set automatically during compile time  
  character(len=128) :: version = '$Id$'

  !> coupler tag that is set automatically during compile time
  character(len=128) :: tag = '$Name$'

  !> FmsXgridXmap_type that holds the exchange grid between different components
  type(FmsXgridXmap_type), save :: xmap_sfc

  !> total number of exchange grid cells
  integer :: n_xgrid_sfc=0

  !> module name used when registering variable for diag_manager
  character(len=4), parameter :: mod_name = 'flux'

  integer :: &
       !> diag_manager register field id for 'drag coefficient for moisture'
       id_drag_moist, & 
       !> diag_manager register field id for 'drag coefficient for heat'
       id_drag_heat, &
       !> diag_manager register field id for 'drag coefficient for momentum'
       id_drag_mom, &
       !> diag_manager register field id for 'surface roughness for moisture'
       id_rough_moist, &
       !> diag_manager register field id for 'surface roughness for heat'
       id_rough_heat, &
       !> diag_manager register field id for 'surface roughness for momentum'
       id_rough_mom, &
       ! diag_manager register field id for 'fractional amount of sea ice'
       id_land_mask, &
       !> diag_manager register field id for 'fractional amount of land'
       id_ice_mask, &
       !> diag_manager register field id for 'friction velocity'
       id_u_star, &
       !> diag_manager register field id for 'bouyancy scale'
       id_b_star, &
       !> diag_manager register field id for 'moisture scale'
       id_q_star, &
       !> diag_manager register field id for 'zonal wind stress'
       id_u_flux, &
       !> diag_manager register field id for 'meridional wind stress'
       id_v_flux, &
       !> diag_manager register field id for 'surface temperature'
       id_t_surf, &
       !> diag_manager register field id for 'surface temperature from ocean output'
       id_t_ocean, &
       !> diag_manager register field id for 'sensible heat flux'
       id_t_flux, &
       !> diag_manager register field id for 'net (down-up) longwave flux'
       id_r_flux, &
       !> diag_manager register field id for 'evaporation rate'
       id_q_flux, &
       !> diag_manager register field id for 'sea level pressure'
       id_slp, &
       !> diag_manager register field id for 'temperature at lowest atmospheric level'
       id_t_atm, &
       !> diag_manager register field id for 'u wind component at lowest atmospheric level'
       id_u_atm, &
       !> diag_manager register field id for 'v wind component at lowest atmospheric level'
       id_v_atm, &
       !> diag_manager register field id for 'wind speed for flux calculations'
       id_wind, &
       !> diag_manager register field id for 'surface air virtual potential temperature'
       id_thv_atm, &
       !> diag_manager register field id for 'surface virtual potential temperature'
       id_thv_surf, &
       !> diag_manager register field id for 'temperature at z_ref_heat'
       id_t_ref, &
       !> diag_manager register field id for 'relative humidity at z_ref_heat'
       id_rh_ref, &
       !> diag_manager register field id for 'zonal wind component at z_ref_mom'
       id_u_ref, &
       !> diag_manager register field id for 'meridional wind component at z_ref_mom'
       id_v_ref, &
       !> diag_manager register field id for 'absolute value of wind at z_ref_mom'
       id_wind_ref, &
       !> diag_manager register field id for 'ref height interp factor for for heat'
       id_del_h, &
       !> diag_manager register field id for 'ref height interp factor for momentum'
       id_del_m, &
       !> diag_manager register field id for 'ref height interp factor for moisture'
       id_del_q, &
       !> diag_manager register field id for 'topographic scaling fractor for momentum drag'
       id_rough_scale, &
       !> diag_manager register field id for 'canopy air temperature'
       id_t_ca, &
       !> diag_manager register field id for 'height of lowest atmospheric level'
       id_z_atm, &
       !> diag_manager register field id for 'pressure at lowest atmospheric level'
       id_p_atm, &
       !> diag_manager register field id for 'gust scale'
       id_gust, &
       !> diag_manager register field id for 'temperature at z_ref_heat over land'
       id_t_ref_land, &
       !> diag_manager register field id for 'relative humidity at z_ref_heat over land'
       id_rh_ref_land, &
       !> diag_manager register field id for 'zonal wind component at z_ref_mom over land'
       id_u_ref_land, &
       !>diag_manager register field id for 'meridional wind component at z_ref_mom over land'
       id_v_ref_land, &
       !> diag_manager register field id for 'specific humidity at z_ref_heat'
       id_q_ref, &
       !> diag_manager register field id for 'specific humidity at z_ref_heat over land'
       id_q_ref_land, &
       !> diag_manager register field id for 'evaporation rate over land'
       id_q_flux_land, &
       !> diag_manager register field id for 'relative humidity at z_ref_heat'
       id_rh_ref_cmip, &
       !>diag_manager register field id for 'near-surface specific humidity on land use tile'
       id_hussLut_land, &
       !> diag_manager register field id for 'near-surface air temperature at
       !! z_ref_heat above displacement height on land-use tile'
       id_tasLut_land, &
       !> diag_manager register field id for 'sensible heat flux over land'
       id_t_flux_land

  integer :: &
       !> diag_manager register field id for 'co2 dry volume mixing ratio at lowest atmospheric level'
       id_co2_atm_dvmr, &
       !> diag_manager register field id for 'c02 dry volume mixing ratio at surface'
       id_co2_surf_dvmr

  integer :: &
       !> diag_manager register field id for 'concentration of co2 to be passed to land/photosynthesis'
       id_co2_bot, &
       !> diag_manager register field id for 'concentration of co2 to be passed to ocean'
       id_co2_flux_pcair_atm, &
       !> diag_manager register field id for 'concentration of o2 to be passed to to ocean'
       id_o2_flux_pcair_atm

  integer, allocatable :: &
       !> array of diag_manager register field ids for 'tracers at lowest atmospheric level'
       id_tr_atm(:), &
       !> array of diag_manager register field ids for 'tracers at surface'
       id_tr_surf(:), &
       !> array of diag_manager register field ids for 'tracers fluxes'
       id_tr_flux(:), &
       !> array of diag_manager register field ids for 'flux of co2 concentration in [mol/m2*s]'
       id_tr_mol_flux(:), &
       !> array of diag_manager register field ids for 'tracers at z_ref_heat'
       id_tr_ref(:), &
       !> array of diag_manager register field ids for 'tracer flux at z_ref_heat over land'
       id_tr_ref_land(:)


  integer, allocatable ::  &
       !> array of diag_manager register field ids for 'gross flux of tracer concentration over land in [mol/m2*s]'
       id_tr_mol_flux0(:) 

  integer, allocatable :: &
       !> array of diag_manager register field ids for 'flux of tracer concentration over land in [kg/m2*s]'
       id_tr_flux_land(:), &
       !> array of diag_manager register field ids for 'flux of tracer concentration over land in [mol/m2*s]'
       id_tr_mol_flux_land(:)
  integer, allocatable :: &
       !> array of diag_manager register field ids for 'deposition velocity at lowest atmospheric level over land'
       !! Used only when _USE_LEGACY_LAND_ macro is set at compile time
       id_tr_con_atm_land(:), &
       !> array of diag_manager register field id for 'deposition velocity at reference height over land'
       id_tr_con_ref_land(:)
  integer, allocatable :: &
       !> array of diag_manager register field ids for 'deposition velocity at lowest atmospheric level (atm)'.
       !! Used only when _USE_LEGACY_LAND_ macro is set at compile time
       id_tr_con_atm(:), &
       !> array of diag_manager register field ids for 'deposition velocity at reference height (atm)'       
       id_tr_con_ref(:)

  ! id's for cmip specific fields
  integer :: &
       !> diag_manager register field id for 'near-surface air temperature' (cmip)
       id_tas, &
       !> diag_manager register field id for 'eastward near-surface wind' (cmip)       
       id_uas, &
       !> diag_manager register field id for 'northward near-surface wind' (cmip)
       id_vas, &
       !> diag_manager register field id for 'surface temperature' (cmip)
       id_ts, &
       !> diag_manager register field id for 'air pressure at sea level' (cmip)
       id_psl, &
       !> diag_manager register field id for 'near-surface wind speed' (cmip)
       id_sfcWind, &
       !> diag_manager register field id for 'surface downward eastward wind stress' (cmip)
       id_tauu, &
       !> diag_manager register field id for 'surface downward northward wind stress' (cmip)
       id_tauv, &
       !> diag_manager register field id for 'near-surface relative humidty' (cmip)
       id_hurs, &
       !> diag_manager register field id for 'near-surface specific humidity' (cmip)
       id_huss, &
       !> diag_manager register field id for 'water evaporation flux' (cmip)
       id_evspsbl, &
       !> diag_manager register field id for 'surface upward latent heat flux' (cmip)
       id_hfls, &
       !> diag_manager register field id for 'surface upward sensible heat flux' (cmip)
       id_hfss, &
       !> diag_manager register field id for 'near-surface relative humidty' (cmip)
       id_rhs, &
       !> diag_manager register field id for 'fraction of the grid cell occupied by land' (cmip)
       id_sftlf, &
       !> diag_manager register field id for 'sea surface temperature' (cmip)
       id_tos, &
       !> diag_manager register field id for 'sea ice area fraction' (cmip)
       id_sic, &
       !> diag_manager register field id for 'surface temperature on land or sea ice' (cmip)
       id_tslsi, &
       !> diag_manager register field id for 'near surface height' (cmip)
       id_height2m, &
       !> diag_manager register field id for 'near surface height' (cmip)
       id_height10m

  integer :: &
       !> diag_manager register field id for 'global integral of water evaporation flux'
       id_evspsbl_g, &
       !> diag_manager register field id for 'global integral of surface temperature'
       id_ts_g, &
       !> diag_manager register field id for 'global integral of near-surface air temperature'
       id_tas_g, &
       !> diag_manager register field id for 'global integral of near-surface air temperature on land'
       id_tasl_g, &
       !> diag_manager register field id for 'global integral of surface upward sensible heat flux'
       id_hfss_g, &
       !> diag_manager register field id for 'global integral of surface upward latent heat flux'
       id_hfls_g, &
       !> diag_manager register field id for 'global integral of near-surface relative humidty'
       id_rls_g

  !>If true, saves land_mask, sftlf, height2m, and height10m once per file at first call to sf_boundary_layer
  logical :: first_static = .true.

  !> true if atm_land_ice_flux_exchnge_init has been called
  logical :: do_init = .true.

  !> first or second order conservative remapping onto exchange grid
  integer :: remap_method = 1

  !> rdgas/rvgas
  real, parameter :: d622 = rdgas/rvgas
  !> 1.0-d622
  real, parameter :: d378 = 1.0-d622
  !> d378/d622  
  real, parameter :: d608   = d378/d622
  !> freezing point of water at 1 atm [K]
  real, parameter :: tfreeze = 273.15
  !> frac_precip
  real, allocatable, dimension(:,:) :: frac_precip

  !> Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
  real  :: z_ref_heat =  2. 
  !> Reference height (meters) for mementum diagnostics (u_ref, v_ref, del_m)
  real :: z_ref_mom  = 10.

  !> do_forecast
  logical :: do_forecast = .false.

  !> OpenMP number of thread.  Do loops on the exchange grid are parallelized into noblocks
  integer :: nblocks = 1
  
  !> If true, convert liquid precip to snow when t_ref < tfreeze
  !! Used for atm override experiments where liquid and frozen precip are combined
  logical :: partition_fprec_from_lprec = .FALSE.

  !> If true, scale mass of liqud preciptation
  logical :: scale_precip_2d = .false.

  !>  Initializing OpenMP parameter
  integer :: my_nblocks = 1
  integer, allocatable :: &
       !> starting do loop indices for OpenMP thread
       block_start(:), &
       !> ending do loop indices for OpenMP thread
       block_end(:)


  real, allocatable, dimension(:) :: &
       !> surface temperature for radiation calc on exchange grid [K]
       !! Note, T canopy is only differet from t_surf over vegetated land
       ex_t_surf, &
       !> no documentation
       ex_t_surf_miz, &  
       !> near-surface (canopy) air temperature on exchange grid [K]
       ex_t_ca, &
       !> surface pressure on exchange grid on the exchange grid
       ex_p_surf, &
       !> surface pressure on exchange grid
       ex_slp, &
       !> sens heat flux on the exchange grid
       ex_flux_t, &
       !> longwave radiation flux on the exchange grid
       ex_flux_lw, &
       !> d(sens.heat.flux)/d(T canopy) on the exchange grid
       ex_dhdt_surf, &
       !> d(water.vap.flux)/d(T canopy) on the exchange grid
       ex_dedt_surf, &
       !> d(water.vap.flux)/d(q canopy) on the exchange grid
       ex_dqsatdt_surf, &
       !> dt/mass * dedet_surf * gamma on the exchange grid
       ex_e_q_n, &
       !> d(LW flux)/d(T surf) on the exchange grid
       ex_drdt_surf, &
       !> d(sens.heat.flux)/d(T atm) on the exchange grid
       ex_dhdt_atm, &
       !> u stress on atmosphere on the exchange grid
       ex_flux_u, &
       !> v stress on atmosphere on the exchange grid
       ex_flux_v, &
       !> d(stress)/d(u) on the exchange grid
       ex_dtaudu_atm, &
       !> d(stress)/d(v) on the exchange grid
       ex_dtaudv_atm, &
       !> mask array of seaice fractions on the exchange grid.
       !! Takes value of 1 when there is any open water in the OCN grid cell.
       !! Takes value of 0 when there is no open water in the OCN grid cell (i.e,
       !! totally covered with ice or land).  Note, ex_seawater should not be
       !! mistaken with ex_avail where ex_avail is 1 for all OCN grid cells regardless
       !! of sea-ice coverage.
       ex_seawater, &
       !> no documentation
       ex_albedo_vis_dir_fix, &
       !> no documentation
       ex_albedo_nir_dir_fix, &
       !> no documentation
       ex_albedo_vis_dif_fix, &
       !> no documentation
       ex_albedo_nir_dif_fix, &
       !> q drag coefficient on the exchange grid
       ex_drag_q, &
       !> drag coefficient for heat on the exchange grid       
       ex_cd_t, &
       !> drag coefficient for momentum on the exchange grid
       ex_cd_m, &
       !> boyuancy scale on the exchange grid
       ex_b_star, &
       !> friction velocity on exchange grid
       ex_u_star, &
       !> wind speed on exchange grid
       ex_wind, &
       !> height of lowest atmospheric level on exchange grid
       ex_z_atm, &
       !> deposition velocity at lowest atmospheric level on the exchange grid
       ex_con_atm

  
#ifdef SCM
  real, allocatable, dimension(:) :: &
       ex_dhdt_surf_forland, &
       ex_dedt_surf_forland, &
       ex_dedq_surf_forland
#endif

  real, allocatable, dimension(:,:) :: &
       !> surface temperature for radiation calc on exchange grid [K]
       ex_tr_surf, & 
       !> tracer fluxes on the exchange grid
       ex_flux_tr, &
       !> d(tracer flux)/d(surf tracer) on the exchange grid
       ex_dfdtr_surf, &
       !> d(tracer flux)/d(atm tracer) on the exchange grid
       ex_dfdtr_atm, &
       !> coefficient in implicit scheme on the exchange grid
       ex_e_tr_n, &
       !> coefficient in implicit scheme on the exchange grid
       ex_f_tr_delt_n 

  real, allocatable, dimension(:,:) :: &
       !> deposition velocity at reference height on the exchange grid
       ex_tr_con_ref, & 
       !> deposition velocity at atmospheric height on the exchange grid
       ex_tr_con_atm    
  
  logical, allocatable, dimension(:) :: &
       !> true where exchange grid cell is over ocean and/or seaice
       ex_avail,& 
       !> true where exchange grid cell is over land
       ex_land 
  real, allocatable, dimension(:) :: &
       !> no documentation
       ex_e_t_n, &
       !> no documentation
       ex_f_t_delt_n

  !> number of prognostic tracers in the atmos model
  integer :: n_atm_tr  
  !> number of prognostic tracers in the atmos model
  integer :: n_atm_tr_tot  
  !> number of prognostic tracers in the land model
  integer :: n_lnd_tr
  !> number of prognostic tracers in the land model
  integer :: n_lnd_tr_tot 
  !> number of tracers exchanged between models
  integer :: n_exch_tr
  !> number of gex fields exchanged between land and atmosphere
  integer :: n_gex_atm2lnd 
  !> number of gex fields exchanged between atmosphere and land
  integer :: n_gex_lnd2atm 

  !> derived type to hold the index of the tracer in atm, ice, land models
  type :: tracer_ind_type
     integer :: atm, ice, lnd 
  end type tracer_ind_type

  !> table of tracers passed through flux exchange
  type(tracer_ind_type), allocatable :: tr_table(:) 

  !> derived type to hold index of the tracer on the exchange grid, ice, and land models
  type :: tracer_exch_ind_type
     integer :: exch = 0
     integer :: ice = 0
     integer :: lnd = 0 
  end type tracer_exch_ind_type

  !> map atm tracers to exchange, ice and land variables
  type(tracer_exch_ind_type), allocatable :: tr_table_map(:)

  !> specific humidity index.  Initialized as NO_TRACER
  integer :: isphum = NO_TRACER
  !> co2 tracer index.  Initialized as NO_TRACER
  integer :: ico2 = NO_TRACER
  !> nh3 tracer index.  Initialized as NO_TRACER
  integer :: inh3 = NO_TRACER

  !> atm gas fields,  Used as place holder for atmospheric fields
  type(FmsCoupler1dBC_type), pointer :: ex_gas_fields_atm=>NULL() 

  !> ice gas fields.  Used as place holder for ice fields 
  type(FmsCoupler1dBC_type), pointer :: ex_gas_fields_ice=>NULL()

  !> gas flux fields.  Used as place holder for intermediate calculations such as piston velocities
  type(FmsCoupler1dBC_type), pointer :: ex_gas_fluxes=>NULL() 

  interface put_logical_to_real
     module procedure put_logical_to_real_sg
     module procedure put_logical_to_real_ug
  end interface

  integer :: &
       !> number of x gridpoints in the atm compute domain
       ni_atm, &
       !> number of y gridpoints in the atm compute domain
       nj_atm
  
  integer, parameter :: &
       !> flag to set boundary_type%xtype when grids are physically different and data between model components
       !! needs to be exchanged via the exchange grid
       REGRID=1, &
       !> flag to set boundary_type%xtype when grids are physically same, but differ in domain decomposition.
       REDIST=2, &
       !> flag to set boundary_type%xtype when grids and domaisn are identical and data can be
       !! copied directly beteween components
       DIRECT=3
  
  integer :: &
       !> FMS clock id for profiling general processes
       cplClock, &
       !> FMS clock id for profiling sfc_boundary_layer
       sfcClock, &
       !> FMS clock id for profiling flux down from atmosphere
       fluxAtmDnClock, &
       !> FMS clock for profiling exchange grid generation
       regenClock, &
       !> FMS clock for profiling flux up to atmosphere
         fluxAtmUpClock
  
  integer :: &
       !> exchange grid index for xgrid_stock_move.  Set to value of 1
       X1_GRID_ATM, &
       !> exchange grid index for xgrid_stock_move.  Set to value of 2
       X1_GRID_ICE, & 
       !> exchange grid index for xgrid_stock_move.  Set to value of 3
       X1_GRID_LND 

  real :: &
       !> atmospheric timestep [s]
       Dt_atm, &
       !> coupled timestep [s]
       Dt_cpl
  
  integer :: &
       !> number of x gridpoints in ice compute domain
       nxc_ice=0, &
       !> number of y gridpoints in ice compute domain
       nyc_ice=0, &
       !> number of vertical levels in ice
       nk_ice=0
  
  integer :: &
       !> number of x gridpoints in land compute domain
       nxc_lnd=0, &
       !> number of y gridpoints in land compute domain
       nyc_lnd=0

contains

  !> Subroutine atm_land_ice_flux_exchange_init initializes atm_land_ice_flux_exchange_mod by
  !! allocating and seting default values for module level variable; and calling initialization routines
  !! in FMS modules.  This subroutine must be called before calling any other public procedures in this
  !! module
  subroutine atm_land_ice_flux_exchange_init(Time, Atm, Land, Ice, atmos_ice_boundary, land_ice_atmos_boundary, &
       Dt_atm_in, Dt_cpl_in, z_ref_heat_in, z_ref_mom_in, do_area_weighted_flux_in, do_forecast_in, &
       partition_fprec_from_lprec_in, scale_precip_2d_in, nblocks_in, cplClock_in, ex_gas_fields_atm_in, &
       ex_gas_fields_ice_in, ex_gas_fluxes_in)

    !> current model time    
    type(FmsTime_type), intent(in) :: Time
    !> derived data type holding atmosphere boundary data
    type(atmos_data_type), intent(inout) :: Atm
    !> derived data type holding land boundary data
    type(land_data_type), intent(in) :: Land
    !>derived data type holding ice boundary data
    type(ice_data_type), intent(inout) :: Ice 
    !> derived type holding properties and fluxes passed from atmosphere to ice
    type(atmos_ice_boundary_type), intent(inout) :: atmos_ice_boundary
    !> derived type holding properties and fluxes passed from exchange grid to atmosphere, land, and ice
    type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
    !> used to set dt_atm (atmosphere time step [s]) in the module
    real, intent(in) :: Dt_atm_in
    !> used to set dt_cpl (coupled time step [s]) in the module
    real, intent(in) :: Dt_cpl_in
    !> used to set z_ref_heat (reference height for temperature and relative humidity diagnostics [m]) in the module
    real, intent(in) :: z_ref_heat_in,
    !> used to set z_ref_mom (reference height for momentum diagnostics [m]) in the module
    real, intent(in) :: z_ref_mom_in
    !> used to set scale_precip_2d in the module.  if true, rescale Atm%lprec
    logical, intent(in) :: scale_precip_2d_in
    !> used to set do_area_weighted_flux in the module.  if true, divide flux by area
    logical, intent(in) :: do_area_weighted_flux_in
    !> used to set do_forecast in the module
    !! if true, and #ifdef AM3_physics,, put atm%surf_diff%sst_miz on the exchange grid
    logical, intent(in) :: do_forecast_in
    !> used to set partition_fprec_from_lprec in the module
    !! if true, convert liquid precip to snow when t_ref < tfreeze
    logical, intent(in) :: partition_fprec_from_lprec_in
    !> used to set nblocks (number of OpenMP blocks) in the module.
    integer, intent(in) :: nblocks_in
    !> used to set cplClock in the module.
    !! The clock is used to measure processes mainly used for development and debugging
    integer, intent(in)  :: cplClock_in
    
    type(FmsCoupler1dBC_type), intent(in), target :: &
         !> used to set ex_gas_fields_atm in the module.
         !! Contains atm surface variables used for computing atm-ocean gas fluxes and flux-regulating parameters
         ex_gas_fields_atm_in, &
         !> used to set ex_gas_fields_ice in the module.  Contains ice-top and ocean surface variables
         !! used for computing atm-ocean gas fluxes and flux-regulating parameters
         ex_gas_fields_ice_in, &
         !> used to set ex_gas_fluxes in the module that is used to exchange gas/tracer fluxes between atm and ocean
         ex_gas_fluxes_in
    
    character(len=48), parameter :: module_name = 'atm_land_ice_flux_exchange_mod'
    character(len=64), parameter :: sub_name = 'atm_land_ice_flux_init'
    character(len=256), parameter :: note_header = '==>Note from '//trim(module_name)//'('//trim(sub_name)//'):'

    integer :: &
         i, & ! temporary index do loop
         n ! temporary index for counting
    integer :: &
         outunit, & ! returned value from fms_mpp_stdout() 
         logunit ! returned value from fms_mpp_stdlog()
    integer :: &
         is, & ! starting x-index on compute domain
         ie, & ! ending x-index on compute domain
         js, & ! starting y-index on compute domain
         je, & ! ending y-index on compute domain
         kd ! number of levels in the z direction

    character(32) :: tr_name ! dummy variable to hold name of tracers
    logical :: found ! dummy variable to search through tracer index in ex_gas_fluxes

    !> INITIALIZE MODULE-LEVEL VARIABLES
    !{
    Dt_atm = Dt_atm_in
    Dt_cpl = Dt_cpl_in
    z_ref_heat = z_ref_heat_in
    z_ref_mom = z_ref_mom_in
    do_area_weighted_flux = do_area_weighted_flux_in
    do_forecast = do_forecast_in
    partition_fprec_from_lprec = partition_fprec_from_lprec_in
    scale_precip_2d = scale_precip_2d_in
    nblocks = nblocks_in
    cplClock = cplClock_in
    ex_gas_fields_atm => ex_gas_fields_atm_in
    ex_gas_fields_ice => ex_gas_fields_ice_in
    ex_gas_fluxes     => ex_gas_fluxes_in
    !}
    

    !> GET FILE UNIT FOR STDOUT AND STDLOG FOR INTERNAL LOGGING PURPOSES
    !{
    outunit = fms_mpp_stdout()
    logunit = fms_mpp_stdlog()
    !}

    !< ALLOCATE OPENMP BLOCK_START ND BLOCK_END HOLDING DO LOOP INDICES
    !{
    allocate(block_start(nblocks), block_end(nblocks))
    !}

    
    !> FROM THE TRACER TABLE, GET THE TOTAL NUMBER TRACERS,
    !! TOTAL NUMBER OF SPECIFIC HUMIDITY TRACER,AND
    !! TOTAL NUMBER PROGNOSTIC TRACERS IN ATMOSPHERE AND LAND MODELS
    !{
    call fms_tracer_manager_get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, num_prog=n_atm_tr)
    call fms_tracer_manager_get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, num_prog=n_lnd_tr)
    !}

    
    !> CONSTRUCT THE TRACER TABLE (TR_TABLE):
    !! FOR EACH TRACER, RECORD THE TRACER_INDEX IN THE ATM MODEL, ICE MODEL, AND LAND MODEL
    !! SKIP ALL ATMOS TRACERS THAT DO NOT HAVE CORRESPONDING SURFACE TRACERS
    !{
    allocate(tr_table(n_atm_tr))
    allocate(tr_table_map(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, i, tr_name )
       tr_table(n)%atm = i
       tr_table(n)%ice = fms_tracer_manager_get_tracer_index ( MODEL_ICE,  tr_name )
       tr_table_map(i)%ice = tr_table(n)%ice
       tr_table(n)%lnd = fms_tracer_manager_get_tracer_index ( MODEL_LAND, tr_name )
       tr_table_map(i)%lnd = tr_table(n)%lnd
       if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) then
          tr_table_map(i)%exch = n
          n = n + 1
       endif
    enddo
    n_exch_tr = n - 1
    !}

    
    !> GET THE TOTAL NUMBER OF GENERIC EXCHANGE FIELDS BETWEEN ATMOSPHERE AND LAND
    !{
    !generic exchange
    n_gex_atm2lnd = fms_gex_get_n_ex(MODEL_ATMOS,MODEL_LAND)
    if (fms_mpp_root_pe().eq.fms_mpp_pe()) write(*,*) 'atm_land_ice_flux_exchange_init [gex]',n_gex_atm2lnd
    n_gex_lnd2atm = fms_gex_get_n_ex(MODEL_LAND,MODEL_ATMOS)
    if (fms_mpp_root_pe().eq.fms_mpp_pe()) write(*,*) 'atm_land_ice_flux_exchange_init [gex]',n_gex_lnd2atm

    do n = 1, ex_gas_fluxes%num_bcs
       if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then
          found = .false.
          do i = 1, n_exch_tr
             if (ex_gas_fluxes%bc(n)%atm_tr_index .eq. tr_table(i)%atm) then
                found = .true.
                exit
             endif
          enddo
          if (.not. found) then
             n_exch_tr = n_exch_tr + 1
             tr_table(n_exch_tr)%atm = ex_gas_fluxes%bc(n)%atm_tr_index
             tr_table(n_exch_tr)%ice = NO_TRACER ! because ocean-atm gas fluxes are not held in the ice model as tracers
             tr_table(n_exch_tr)%lnd = NO_TRACER ! because this would have been found above
             tr_table_map(n_exch_tr)%exch = n_exch_tr
             tr_table_map(n_exch_tr)%ice = tr_table(n_exch_tr)%ice
             tr_table_map(n_exch_tr)%lnd = tr_table(n_exch_tr)%lnd
          endif
       endif
    enddo
    write(outunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
    write(logunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
    do i = 1,n_exch_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
       write(outunit,*)'Tracer field name :'//trim(tr_name)
       write(logunit,*)'Tracer field name :'//trim(tr_name)
    enddo
    !}

    
    !>  GET THE TRACER INDEX OF SPECIFIC HUMIDITY
    !{
    ! +fix-me-slm+ specific humidity may not be present if we are running with
    ! dry atmosphere. Besides, model may use mixing ratio ('mix_rat') (?). However,
    ! some atmos code also assumes 'sphum' is present, so for now the following
    ! code may be good enough.
    do i = 1,n_exch_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
       if(fms_mpp_lowercase(tr_name)=='sphum') then
          isphum = i
       endif
       ! jgj: find out which exchange tracer is co2
       if(fms_mpp_lowercase(tr_name)=='co2') then
          ico2 = i
          write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',ico2
       endif
       if(fms_mpp_lowercase(tr_name)=='nh3') then
          inh3 = i
          write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',inh3
       endif
    enddo

    if (isphum==NO_TRACER) then
       call fms_error_mesg('atm_land_ice_flux_exchange_mod',&
            'tracer "sphum" must be present in the atmosphere', FATAL )
    endif

    if (ico2==NO_TRACER) then
       call fms_error_mesg('atm_land_ice_flux_exchange_mod',&
            'tracer "co2" not present in the atmosphere', NOTE )
    endif
    !}
    
    !> INITIALIZE FRAC_PRECIP IF SCALE_PRECIP_2D IS TRUE
    !{
    call fms_mpp_domains_get_compute_domain(Atm%domain, is, ie, js, je)

    if (scale_precip_2d) then
       allocate(frac_precip(is:ie,js:je))
       frac_precip=0.0
    endif
    !}

    
    !> SET UP THE EXCHANGE GRID AND SET X1_GRID_ATM = 1, X1_GRID_ICE = 2, AND X1_GRID_LAND = 3
    !! SETS XMAP_SFC(1)%GRIDS FOR ATM, XMAP_SFC(2)%GRIDS FOR ICE, XMAP_SFC(3)%GRIDS FOR LAND
    !{
    call fms_xgrid_init(remap_method)
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_setup_xmap(xmap_sfc, ['ATM', 'OCN', 'LND'], &
         [Atm%Domain, Ice%Domain, Land%Domain], &
         "INPUT/grid_spec.nc", Atm%grid, lnd_ug_domain=Land%ug_domain)
#else
    call fms_xgrid_setup_xmap(xmap_sfc, ['ATM', 'OCN', 'LND'], &
         [Atm%Domain, Ice%Domain, Land%Domain], &
         "INPUT/grid_spec.nc", Atm%grid)
#endif
    ! exchange grid indices
    X1_GRID_ATM = 1
    X1_GRID_ICE = 2
    X1_GRID_LND = 3
    call generate_sfc_xgrid( Land, Ice )
    if (n_xgrid_sfc.eq.1) write (*,'(a,i6,6x,a)') 'PE = ', fms_mpp_pe(), 'Surface exchange size equals one.'
    !}

    
    !>  INITIALIZE SURFACE_FLUX MODULE
    !{
    call surface_flux_init()
    !}

    
    !> INITIALLIZE FMS DIAG_INTEGRAL FIELDS FOR EVAP, T_SURF, T_REF GLOBAL INTEGRAL QUANTITIES
    !{
    !! call diag_integral_field_init ('prec', 'f6.3')
    call fms_diag_integral_field_init ('evap', 'f6.3')
#ifndef use_AM3_physics
    call fms_diag_integral_field_init ('t_surf', 'f10.3')
    call fms_diag_integral_field_init ('t_ref',  'f10.3')
#endif
    !}


    !> REGISTER FMS DIAGNOSTIC FIELDS IN DIAG_MANAGER
    ! all fields will be output on the atmospheric grid
    !{
    call diag_field_init ( Time, Atm%axes(1:2), Land%axes, Land%pe )
    !}

    !> GET THE SIZE OF THE ATM COMPUTE DOMAIN
    !{
    ni_atm = size(Atm%lon_bnd,1)-1 ! to dimension "diag_atm"
    nj_atm = size(Atm%lon_bnd,2)-1 ! in flux_ocean_to_ice
    !}
    

    !> ALLOCATE ATMOS_ICE_BOUNDARY AND SET FIELDS EQUAL TO ZERO
    !{
    call fms_mpp_domains_get_compute_domain( Ice%domain, is, ie, js, je )
    kd = size(Ice%part_size,3)
    allocate( atmos_ice_boundary%u_flux(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%v_flux(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%u_star(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%t_flux(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%q_flux(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%lw_flux(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_flux_vis_dir(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_flux_vis_dif(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_flux_nir_dir(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_flux_nir_dif(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_down_vis_dir(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_down_vis_dif(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_down_nir_dir(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%sw_down_nir_dif(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%lprec(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%fprec(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%dhdt(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%dedt(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%drdt(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%coszen(is:ie,js:je,kd) )
    allocate( atmos_ice_boundary%p(is:ie,js:je,kd) )
    ! initialize boundary values for override experiments
    atmos_ice_boundary%u_flux=0.0
    atmos_ice_boundary%v_flux=0.0
    atmos_ice_boundary%u_star=0.0
    atmos_ice_boundary%t_flux=0.0
    atmos_ice_boundary%q_flux=0.0
    atmos_ice_boundary%lw_flux=0.0
    atmos_ice_boundary%sw_flux_vis_dir=0.0
    atmos_ice_boundary%sw_flux_vis_dif=0.0
    atmos_ice_boundary%sw_flux_nir_dir=0.0
    atmos_ice_boundary%sw_flux_nir_dif=0.0
    atmos_ice_boundary%sw_down_vis_dir=0.0
    atmos_ice_boundary%sw_down_vis_dif=0.0
    atmos_ice_boundary%sw_down_nir_dir=0.0
    atmos_ice_boundary%sw_down_nir_dif=0.0
    atmos_ice_boundary%lprec=0.0
    atmos_ice_boundary%fprec=0.0
    atmos_ice_boundary%dhdt=0.0
    atmos_ice_boundary%dedt=0.0
    atmos_ice_boundary%drdt=0.0
    atmos_ice_boundary%coszen=0.0
    atmos_ice_boundary%p=0.0

    ! allocate fields for extra fluxes
    ! Copying initialized gas fluxes from exchange grid to atmosphere_ice boundary

    call fms_coupler_type_copy(ex_gas_fluxes, atmos_ice_boundary%fluxes, is, ie, js, je, kd, &
         mod_name, Ice%axes, Time, suffix = '_atm_ice')
    
    !--- Ice%ocean_fields and Ice%ocean_fluxes_top will not be passed to ocean, so these two
    !--- coupler_type_copy calls are moved from ice_ocean_flux_init to here.
    if (.not.fms_coupler_type_initialized(Ice%ocean_fields)) &
         call fms_coupler_type_spawn(ex_gas_fields_ice, Ice%ocean_fields, (/is,is,ie,ie/), &
         (/js,js,je,je/), (/1, kd/), suffix = '_ice')
    call fms_coupler_type_set_diags(Ice%ocean_fields, 'ice_flux', Ice%axes, Time)
    !}

    
    !> ALLOCATE LAND_ICE_ATMOS_BOUNDARY AND SET FIELDS EQUAL TO ZERO EXCEPT FOR
    !! T_OCEAN WHICH IS SET TO 200 K, T_REF TO 273 K, AND ROUGHNESS LENGTHS TO 0.01 m
    !{
    call fms_mpp_domains_get_compute_domain( Atm%domain, is, ie, js, je )
    allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%t_ocean(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%u_ref(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%v_ref(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%t_ref(is:ie,js:je) )  ! PBL depth mods
    allocate( land_ice_atmos_boundary%q_ref(is:ie,js:je) )  ! PBL depth mods
    allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%albedo_vis_dir(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%albedo_nir_dir(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%albedo_vis_dif(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%albedo_nir_dif(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%dt_tr(is:ie,js:je,n_atm_tr) )
    allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%dtaudu(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%q_star(is:ie,js:je) )
#ifndef use_AM3_physics
    allocate( land_ice_atmos_boundary%shflx(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%lhflx(is:ie,js:je) )    
#endif
    allocate( land_ice_atmos_boundary%wind(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%thv_atm(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%thv_surf(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%rough_heat(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%frac_open_sea(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%gex_lnd2atm(is:ie,js:je,n_gex_lnd2atm) )

    ! initialize boundary values for override experiments
    land_ice_atmos_boundary%t=273.0
    land_ice_atmos_boundary%t_ocean=200.0
    land_ice_atmos_boundary%u_ref=0.0
    land_ice_atmos_boundary%v_ref=0.0
    land_ice_atmos_boundary%t_ref=273.0   ! PBL depth mods
    land_ice_atmos_boundary%q_ref=0.0     ! PBL depth mods
    land_ice_atmos_boundary%albedo=0.0
    land_ice_atmos_boundary%albedo_vis_dir=0.0
    land_ice_atmos_boundary%albedo_nir_dir=0.0
    land_ice_atmos_boundary%albedo_vis_dif=0.0
    land_ice_atmos_boundary%albedo_nir_dif=0.0
    land_ice_atmos_boundary%land_frac=0.0
    land_ice_atmos_boundary%dt_t=0.0
    land_ice_atmos_boundary%dt_tr=0.0
    land_ice_atmos_boundary%u_flux=0.0
    land_ice_atmos_boundary%v_flux=0.0
    land_ice_atmos_boundary%dtaudu=0.0
    land_ice_atmos_boundary%dtaudv=0.0
    land_ice_atmos_boundary%u_star=0.0
    land_ice_atmos_boundary%b_star=0.0
    land_ice_atmos_boundary%q_star=0.0
    land_ice_atmos_boundary%gex_lnd2atm=0.0

#ifndef use_AM3_physics
    land_ice_atmos_boundary%shflx=0.0
    land_ice_atmos_boundary%lhflx=0.0
#endif
    land_ice_atmos_boundary%wind=0.0
    land_ice_atmos_boundary%thv_atm=0.0
    land_ice_atmos_boundary%thv_surf=0.0
    land_ice_atmos_boundary%rough_mom=0.01
    land_ice_atmos_boundary%rough_heat=0.01
    land_ice_atmos_boundary%frac_open_sea=0.0
    !}
    
    
    !> COPY EX_GAS_FIELDS_ATM TO ATM%FIELDS
    !{
    ! The first call is no longer necessary, the fluxes will be passed by the land module
    ! The 2nd call is useful in the case of a ocean model only simulation
    call fms_coupler_type_copy(ex_gas_fields_atm, Atm%fields, is, ie, js, je, &
         mod_name, Atm%axes(1:2), Time, suffix = '_atm')
    !}
    
    !> GET THE SIZE OF ICE COMPUTE DOMAIN
    !{
    if( Ice%pe) then
       call fms_mpp_domains_get_compute_domain(Ice%domain, xsize=nxc_ice, ysize=nyc_ice)
       nk_ice = size(Ice%part_size,3)
    endif
    !}

    
    !> GET THE SIZE OF LAND COMPUTE DOMAIN
    !{
    if( Land%pe) then
       call fms_mpp_domains_get_compute_domain(Land%domain, xsize=nxc_lnd, ysize=nyc_lnd)
    endif
    !}

    
    !> DECLARE CLOCKS FOR PROFILING
    !{
    sfcClock = fms_mpp_clock_id( 'SFC boundary layer', flags=fms_clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    fluxAtmDnClock = fms_mpp_clock_id( 'Flux DN from atm', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )
    regenClock = fms_mpp_clock_id( 'XGrid generation', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )
    fluxAtmUpClock = fms_mpp_clock_id( 'Flux UP to atm', flags=fms_clock_flag_default, grain=CLOCK_ROUTINE )
    !}

    
    !> SET DO_INIT = .FALSE. IN ORDER TO AVOID RE-INITIALIZATION THE MODULE
    !! IF THIS SUBROUTINE IS CALLED AGAIN
    !{
    do_init = .false.
    !}

    
  end subroutine atm_land_ice_flux_exchange_init

  !> Subroutine sfc_boundary_layer computes the following fluxes and exchanges the fluxes between the model components:
  !!   t_surf_atm: surface temperature used for radiation [K]
  !!   albedo_atm: surface albedo used for radiation  [dimensionless]
  !!   rough_mom_atm: surface roughness for momentum [m]
  !!   land_frac_atm: fractional area of land beneath an atmospheric grid box
  !!   dtaudu_atm, dtaudv_atm: derivatives of wind stress wrt the lowest level wind speed [Pa/(m/s)]
  !!   flux_u_atm: zonal wind stress  [Pa]
  !!   flux_v_atm: meridional wind stress [Pa]
  !!   u_star_atm: friction velocity [m/s]
  !!   b_star_atm: buoyancy scale  [m2/s]
  !!
  !! To exchange fluxes between model components, data on one model component grid is first mapped onto
  !! the exchange grid.  Fluxes are computed on the exchange grid.  Then, data is mapped from the
  !! exchange grid to the receiving model component grid. Note, computed fields and fluxes can be overwritten
  !! with calls to data_override where data will be overwritten only if the tracers are specified in the
  !! tracer_table.
  !! \note `u_star` and `b_star` are defined so that `u_star**2` is the magnitude of surface stress
  !! divided by density of air at the surface, and `u_star*b_star` is the buoyancy flux at the surface.
  subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Land_Ice_Atmos_Boundary )

    !> timestep
    real, intent(in) :: dt 
    !> current model time
    type(FmsTime_type), intent(in) :: Time 
    !> derived type holding atmosphere boundary data
    type(atmos_data_type), intent(inout) :: Atm
    !> derived type holding land boundary data
    type(land_data_type), intent(inout) :: Land
    !> derived type holding ice boundary data
    type(ice_data_type), intent(inout) :: Ice
    !> derived type holding properties and fluxes passed between land and ice to atm
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary

    real, dimension(n_xgrid_sfc) :: &
         ex_albedo, ! albedo on xgrid
         ex_albedo_vis_dir, & ! albedo for light with wavelength in visible region of the solar spectrum
         ex_albedo_nir_dir, & ! albedo for light with wavelength in near-ir region of the solar spectrum
         ex_albedo_vis_dif, & ! albedo for "diffuse" light with wavelength in visible region of the solar spectrum
         ex_albedo_nir_dif, & ! albedo for "diffuse" light with wavelength in near-ir region of the solar spectrum
         ex_land_frac, &  ! fractional area of land on grid cell
         ex_t_atm, & ! air temperature at the lowest atmospheric level
         ex_p_atm, & ! pressure at the lowest atmospheric level
         ex_u_atm, & ! u wind component at the lowest atmospheric level
         ex_v_atm, & ! v wind component at the lowest atmospheric level
         ex_gust, & ! gust scale
         ex_t_surf4, & ! (surface temperature) ** 4
         ex_u_surf, & ! u wind component at Earth's surface
         ex_v_surf, &  ! v wind component at Earth's surface
         ex_rough_mom, &  ! momentum roughness length
         ex_rough_heat, & ! heat roughness length
         ex_rough_moist, & ! moisture roughness length
         ex_rough_scale, & ! scale factor for topographic roughness calculation
         ex_q_star, & ! turbulent moisture scale
         ex_thv_atm, & ! surface area theta_v
         ex_thv_surf, & ! surface theta_v
         ex_cd_q, & ! moisture exchange coefficient
         ex_ref, &! specific humidity at z_ref_heat
         ex_ref_u, & ! zonal wind component at z_ref_mom
         ex_ref_v, & ! meridional wind component at z_ref_mom
         ex_u10, & ! zonal wind speed at 10m above the surface
         ex_ref2, & ! quantity on exchange grid
         ex_t_ref, & ! temperature at z_ref_heat on exchange grid
         ex_qs_ref, & ! quantity on exchange grid
         ex_qs_ref_cmip, & ! < quantity on exchange grid
         ex_del_m, & ! reference height for interpolation factor for momentum
         ex_del_h, & ! reference height interpolation factor for heat
         ex_del_q, & ! reference height interpation factor for moisture
         ex_frac_open_sea ! open-water mask

    real :: rho ! quantity used to compute ex_tr_con_atm
    real, dimension(n_xgrid_sfc,n_exch_tr) :: &
         ex_tr_atm, & ! concentration of tracer at lowest atmospheric level
         ex_tr_ref ! concentration of tracer at reference height

    real, dimension(n_xgrid_sfc) :: ex_co2_atm_dvmr ! added for co2_atm diagnostic
    real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: &
         diag_atm ! temporary array to hold data

#ifndef _USE_LEGACY_LAND_
    real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2)) :: diag_land
    real, dimension(size(Land%t_ca, 1))                   :: diag_land_ug, tile_size_ug
    real, dimension(nxc_lnd,nyc_lnd)                      :: diag_land_sg, tile_size_sg
    logical, dimension(size(Land%t_ca, 1))                :: mask_ug
    logical, dimension(nxc_lnd,nyc_lnd)                   :: mask_sg
    integer :: k
#else
    real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: &
         diag_land ! temporary array to hold data
#endif
    
    real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2),size(Ice%t_surf,3)) :: &
         sea ! temporary array to hold data
    
    real, dimension(size(Ice%albedo,1),size(Ice%albedo,2),size(Ice%albedo,3)) :: &
         tmp_open_sea ! temporary array to hold data

    real :: zrefm ! reference height for computing surface fluxes with Monin-Obukhov similarity theory
    real :: zrefh ! reference height for computing surface fluxes with Monin-Obukhov similarity theory
    logical :: used ! returned value from data_override.  if true, data was overwritten
    character(32) :: tr_name, tr_units ! tracer name and tracer unit
    integer :: tr, n, m ! tracer indices
    integer :: is, ie, isc, iec, jsc, jec ! domain indices
    integer :: i, l, j, n_gex ! do loop indices
    
    real, dimension(n_xgrid_sfc,n_gex_lnd2atm) ::  ex_gex_lnd2atm ! holds generic, non-tracer fields on exchange grid

    !> CHECK MODULE INITIALIZATION
    !{
    ! [1] check that the module was initialized
    if (do_init) then
       call fms_error_mesg ('atm_land_ice_flux_exchange_mod', 'must call atm_land_ice_flux_exchange_init first', FATAL)
    end if
    !}

    
    !> INITIALIZE CLOCKS FOR PROFILING
    !{
    ! [2]
    !Balaji
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(sfcClock)
    !}

    
    !> ALLOCATE ARRAY FOR EXCHANGE FIELDS.  THE ARRAYS ARE DEALLOCATED IN FLUX_UP_TO_ATMOS
    !{
    ! [2] allocate storage for variables that are also used in flux_up_to_atmos
    allocate ( &
         ex_t_surf   (n_xgrid_sfc),  &
         ex_t_surf_miz(n_xgrid_sfc), &
         ex_p_surf   (n_xgrid_sfc),  &
         ex_slp      (n_xgrid_sfc),  &
         ex_t_ca     (n_xgrid_sfc),  &
         ex_dhdt_surf(n_xgrid_sfc),  &
         ex_dedt_surf(n_xgrid_sfc),  &
         ex_dqsatdt_surf(n_xgrid_sfc),  &
         ex_drdt_surf(n_xgrid_sfc),  &
         ex_dhdt_atm (n_xgrid_sfc),  &
         ex_flux_t   (n_xgrid_sfc),  &
         ex_flux_lw  (n_xgrid_sfc),  &
         ex_drag_q   (n_xgrid_sfc),  &
         ex_avail    (n_xgrid_sfc),  &
         ex_f_t_delt_n(n_xgrid_sfc), &

         ex_tr_surf     (n_xgrid_sfc, n_exch_tr), &
         ex_dfdtr_surf  (n_xgrid_sfc, n_exch_tr), &
         ex_dfdtr_atm   (n_xgrid_sfc, n_exch_tr), &
         ex_flux_tr     (n_xgrid_sfc, n_exch_tr), &
         ex_f_tr_delt_n (n_xgrid_sfc, n_exch_tr), &
         ex_e_tr_n      (n_xgrid_sfc, n_exch_tr), &
         ex_con_atm(n_xgrid_sfc), &

         ex_tr_con_ref(n_xgrid_sfc, n_exch_tr), &
         ex_tr_con_atm(n_xgrid_sfc, n_exch_tr), &

         ! MOD these were moved from local ! so they can be passed to flux down
         ex_flux_u(n_xgrid_sfc),    &
         ex_flux_v(n_xgrid_sfc),    &
         ex_dtaudu_atm(n_xgrid_sfc),&
         ex_dtaudv_atm(n_xgrid_sfc),&
         ex_seawater(n_xgrid_sfc),  &

         ! values added for LM3
         ex_cd_t     (n_xgrid_sfc),  &
         ex_cd_m     (n_xgrid_sfc),  &
         ex_b_star   (n_xgrid_sfc),  &
         ex_u_star   (n_xgrid_sfc),  &
         ex_wind     (n_xgrid_sfc),  &
         ex_z_atm    (n_xgrid_sfc),  &

         ex_e_t_n    (n_xgrid_sfc),  &
         ex_e_q_n    (n_xgrid_sfc),  &
         ex_land     (n_xgrid_sfc)   )

#ifdef SCM
    allocate ( &
         ex_dhdt_surf_forland(n_xgrid_sfc), &
         ex_dedt_surf_forland(n_xgrid_sfc), &
         ex_dedq_surf_forland(n_xgrid_sfc)  )
#endif
    !}
    

    !> ALLOCATE EX_GAS_FIELDS_ICE ARRAYS FOR OCEAN_ICE_BOUNDARY EXCHANGE FIELDS
    !{
    ex_p_surf = 1.0
    ! Actual allocation of exchange fields for ocean_ice boundary
    do n = 1, ex_gas_fields_ice%num_bcs
       do m = 1, ex_gas_fields_ice%bc(n)%num_fields
          if (associated(ex_gas_fields_ice%bc(n)%field(m)%values)) then
             call fms_mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_ice already allocated.' )
          endif
          allocate ( ex_gas_fields_ice%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fields_ice%bc(n)%field(m)%values = 0.0
       enddo
    enddo
    !}

    
    !> ALLOCATE EX_GAS_FIELDS_ATM ARRAYS FOR ATMOSPHERE EXCHANGE FIELDS
    !{
    do n = 1, ex_gas_fields_atm%num_bcs
       do m = 1, ex_gas_fields_atm%bc(n)%num_fields
          if (associated(ex_gas_fields_atm%bc(n)%field(m)%values)) then
             call fms_mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_atm already allocated.' )
          endif
          allocate ( ex_gas_fields_atm%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fields_atm%bc(n)%field(m)%values = 0.0
       enddo
    enddo
    !}

    !> ALLOCATE EX_GAS_FLUXES FOR ADDITIONAL EXCHANGE FIELDS
    !{
    do n = 1, ex_gas_fluxes%num_bcs
       do m = 1, ex_gas_fluxes%bc(n)%num_fields
          if (associated(ex_gas_fluxes%bc(n)%field(m)%values)) then
             call fms_mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fluxes already allocated.' )
          endif
          allocate ( ex_gas_fluxes%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fluxes%bc(n)%field(m)%values = 0.0
       enddo
    enddo
    !}

    
    ! Call the atmosphere tracer driver to gather the data needed for extra gas tracers
    ! For ocean only model
    ! call atmos_get_fields_for_flux(Atm)

    
    !> ON THE EXCHANGE GRID, SET INITIAL VALUES FOR ALBEDO, DRAG COEFFICIENTS, AND OPEN WATER MASK
    !{
    ! [3] initialize some values on exchange grid: this is actually a safeguard
    ! against using undefined values
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_t_surf,ex_u_surf, &
    !$OMP                                  ex_v_surf,ex_albedo,ex_albedo_vis_dir,ex_albedo_nir_dir, &
    !$OMP                                  ex_albedo_vis_dif,ex_albedo_nir_dif,ex_cd_t,ex_cd_m,  &
    !$OMP                                  ex_cd_q,ex_frac_open_sea,n_gex_lnd2atm,ex_gex_lnd2atm) &
    !$OMP                          private(is,ie,n_gex)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is,ie
          ex_t_surf(i)   = 200.
          ex_u_surf(i)   =   0.
          ex_v_surf(i)   =   0.
          ex_albedo(i) = 0.
          ex_albedo_vis_dir(i) = 0.
          ex_albedo_nir_dir(i) = 0.
          ex_albedo_vis_dif(i) = 0.
          ex_albedo_nir_dif(i) = 0.

          !---- do not use if relax time /= 0 ----
          ex_cd_t(i) = 0.0
          ex_cd_m(i) = 0.0
          ex_cd_q(i) = 0.0
          ex_frac_open_sea(i) =0.
       end do
       do n_gex=1,n_gex_lnd2atm
         do i = is,ie
               ex_gex_lnd2atm(i,n_gex)   = 0.0
         enddo
       enddo
    enddo
    !}

    
    !> OVERRIDE SUBSET OF ATM FIELDS.  DATA WILL BE OVERWRITTEN ONLY IF FIELD IS SPECIFIED IN DATA_TABLE
    !{
    ! data_override stuff moved from coupler_main
    call fms_data_override ('ATM', 't_bot',  Atm%t_bot , Time)
    call fms_data_override ('ATM', 'z_bot',  Atm%z_bot , Time)
    call fms_data_override ('ATM', 'p_bot',  Atm%p_bot , Time)
    call fms_data_override ('ATM', 'u_bot',  Atm%u_bot , Time)
    call fms_data_override ('ATM', 'v_bot',  Atm%v_bot , Time)
    call fms_data_override ('ATM', 'p_surf', Atm%p_surf, Time)
    call fms_data_override ('ATM', 'slp',    Atm%slp,    Time)
    call fms_data_override ('ATM', 'gust',   Atm%gust,   Time)
    !}

    
    !> CONVERT CO2 TRACER UNITS TO WET_MMR UNITS
    ! jgj: 2008/07/18
    ! FV atm advects tracers in moist mass mixing ratio: kg co2 /(kg air + kg water)
    ! cubed sphere advects moist mass mixing ratio also (per SJ)
    ! data table co2 overrides for ocean (co2_flux_pcair_atm)
    ! and land (co2_bot) should be in dry vmr (mol/mol) units.
    !  ATM: co2_flux_pcair_atm : to override atm_btm layer to send to ocean
    !  ATM: co2_bot            : to override atm_btm layer to send to land

    ! data override for co2 to be passed to land/photosynthesis (co2_bot)
    ! land co2 data override is in dry_vmr units, so convert to wet_mmr for land model.
    ! co2mmr = (wco2/wair) * co2vmr;  wet_mmr = dry_mmr * (1-Q)
    do tr = 1,n_atm_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr, tr_name )
       call fms_data_override('ATM', trim(tr_name)//'_bot', Atm%tr_bot(:,:,tr), Time, override=used)
       ! conversion for land co2 data override from dry vmr to moist mmr
       if (used .and. fms_mpp_lowercase(trim(tr_name)).eq.'co2') then
          ! 2017/08/08 jgj add co2_bot diagnostic in dry_vmr units
          if ( id_co2_bot > 0 ) used = fms_diag_send_data ( id_co2_bot, Atm%tr_bot(:,:,tr), Time )

          isc = lbound(Atm%tr_bot,1); iec = ubound(Atm%tr_bot,1)
          jsc = lbound(Atm%tr_bot,2); jec = ubound(Atm%tr_bot,2)
          !$OMP parallel do default(none) shared(isc,iec,jsc,jec,Atm,tr,isphum)
          do j = jsc, jec
             do i = isc, iec
                Atm%tr_bot(i,j,tr) = Atm%tr_bot(i,j,tr) * (WTMCO2/WTMAIR) *    &
                     (1.0 - Atm%tr_bot(i,j,isphum))
             enddo
          enddo
       end if
    enddo
    !}


    !> OVERRIDE CO2 VALUES WHERE DATA WILL BE OVERWRITTEN ONLY IF FIELD IS SPECIFIED IN DATA_TABLE, AND
    !! SEND DATA TO THE DIAG_MANAGER BUFFER
    !{
    ! data override for co2 to be passed to ocean (co2_flux_pcair_atm)
    ! atmos_co2.F90 already called: converts tr_bot passed to ocean via gas_flux
    ! from moist mmr to dry vmr.
    do n = 1, atm%fields%num_bcs
       do m = 1, atm%fields%bc(n)%num_fields
          call fms_data_override('ATM', atm%fields%bc(n)%field(m)%name, &
               atm%fields%bc(n)%field(m)%values, Time, override = atm%fields%bc(n)%field(m)%override)
          ex_gas_fields_atm%bc(n)%field(m)%override = atm%fields%bc(n)%field(m)%override
          if ( atm%fields%bc(n)%field(m)%override .and. &
               fms_mpp_lowercase(trim(atm%fields%bc(n)%field(m)%name)) .eq. 'co2_flux_pcair_atm') then
             if( id_co2_flux_pcair_atm > 0 ) &
                  used = fms_diag_send_data ( id_co2_flux_pcair_atm, atm%fields%bc(n)%field(m)%values, Time )
          endif
          if ( atm%fields%bc(n)%field(m)%override .and. &
               fms_mpp_lowercase(trim(atm%fields%bc(n)%field(m)%name)) .eq. 'o2_flux_pcair_atm') then
             if( id_o2_flux_pcair_atm > 0 ) &
                used = fms_diag_send_data ( id_o2_flux_pcair_atm, atm%fields%bc(n)%field(m)%values, Time )
          endif
       enddo
    enddo
    do n = 1, atm%fields%num_bcs
       if (atm%fields%bc(n)%use_atm_pressure) then
          if (.not. atm%fields%bc(n)%field(fms_coupler_ind_psurf)%override) then
             atm%fields%bc(n)%field(fms_coupler_ind_psurf)%values = Atm%p_surf
          endif
       endif       
    enddo
    !}

    
    !> OVERRIDE SUBSET OF ICE AND LAND FIELD.  DATA WILL BE OVERWRITTEN ONLY IF THE FIELD IS SPECIFIED IN DATA_TABLE
    !{
    call fms_data_override ('ICE', 't_surf',     Ice%t_surf,      Time)
    call fms_data_override ('ICE', 'rough_mom',  Ice%rough_mom,   Time)
    call fms_data_override ('ICE', 'rough_heat', Ice%rough_heat,  Time)
    call fms_data_override ('ICE', 'rough_moist',Ice%rough_moist, Time)
    call fms_data_override ('ICE', 'albedo',     Ice%albedo,      Time)
    call fms_data_override ('ICE', 'albedo_vis_dir', Ice%albedo_vis_dir, Time)
    call fms_data_override ('ICE', 'albedo_nir_dir', Ice%albedo_nir_dir, Time)
    call fms_data_override ('ICE', 'albedo_vis_dif', Ice%albedo_vis_dif, Time)
    call fms_data_override ('ICE', 'albedo_nir_dif', Ice%albedo_nir_dif, Time)
    call fms_data_override ('ICE', 'u_surf',     Ice%u_surf,      Time)
    call fms_data_override ('ICE', 'v_surf',     Ice%v_surf,      Time)
    call fms_coupler_type_data_override('ICE', Ice%ocean_fields, Time)
    call fms_coupler_type_send_data(Ice%ocean_fields, Time)
#ifndef _USE_LEGACY_LAND_
    call fms_data_override_ug ('LND', 't_surf',     Land%t_surf,     Time)
    call fms_data_override_ug ('LND', 't_ca',       Land%t_ca,       Time)
    call fms_data_override_ug ('LND', 'rough_mom',  Land%rough_mom,  Time)
    call fms_data_override_ug ('LND', 'rough_heat', Land%rough_heat, Time)
    call fms_data_override_ug ('LND', 'albedo', Land%albedo,     Time)
#else
    call fms_data_override ('LND', 't_surf',     Land%t_surf,     Time)
    call fms_data_override ('LND', 't_ca',       Land%t_ca,       Time)
    call fms_data_override ('LND', 'rough_mom',  Land%rough_mom,  Time)
    call fms_data_override ('LND', 'rough_heat', Land%rough_heat, Time)
    call fms_data_override ('LND', 'albedo', Land%albedo,     Time)
#endif

    ! tracer data override
    do tr = 1, n_lnd_tr
       call fms_tracer_manager_get_tracer_names( MODEL_LAND, tr, tr_name )
#ifndef _USE_LEGACY_LAND_
       call fms_data_override_ug('LND', trim(tr_name)//'_surf', Land%tr(:,:,tr), Time)
    enddo
    call fms_data_override_ug ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
    call fms_data_override_ug ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
    call fms_data_override_ug ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
    call fms_data_override_ug ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)
#else
       call fms_data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
    enddo
    call fms_data_override ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
    call fms_data_override ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
    call fms_data_override ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
    call fms_data_override ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)
#endif
    !}

    
    !> MAP ATM FIELDS ONTO THE EXCHANGE GRID
    !{
    ! [4] put all the qantities we need onto exchange grid
    ! [4.1] put atmosphere quantities onto exchange grid
#ifdef use_AM3_physics
    if (do_forecast) then
       call fms_xgrid_put_to_xgrid (Atm%Surf_diff%sst_miz , 'ATM', ex_t_surf_miz, &
       xmap_sfc, remap_method=remap_method, complete=.false.)
    endif
#endif
    ! put atmosphere bottom layer tracer data onto exchange grid
    do tr = 1,n_exch_tr
       call fms_xgrid_put_to_xgrid (Atm%tr_bot(:,:,tr_table(tr)%atm) , 'ATM', ex_tr_atm(:,tr), xmap_sfc, &
            remap_method=remap_method, complete=.false.)
    enddo
    do n = 1, Atm%fields%num_bcs
      if(ex_gas_fields_atm%bc(n)%flux_type  .ne. 'air_sea_deposition') then
       do m = 1, Atm%fields%bc(n)%num_fields
          call fms_xgrid_put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
               ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method, complete=.false.)
       enddo
      endif
    enddo

    call fms_xgrid_put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%slp,    'ATM', ex_slp,    xmap_sfc, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc, remap_method=remap_method, complete=.true.)
    !}

    
    ! slm, Mar 20 2002: changed order in whith the data transferred from ice and land
    ! grids, to fill t_ca first with t_surf over ocean and then with t_ca from
    ! land, where it is different from t_surf. It is mostly to simplify
    ! diagnostic, since surface_flux calculations distinguish between land and
    ! not-land anyway.


    !> INITIALIZE EX_TR_SURF TO BE THE AMOUNT OF TRACERS AT THE BOTTOM-MOST ATMOSPHERE LAYER
    !{
    ! prefill surface values with atmospheric values before putting tracers
    ! from ice or land, so that gradient is 0 if tracers are not filled
    ex_tr_surf = ex_tr_atm
    !}

    
    !> MAP ICE FIELDS ONTO THE EXCHANGE GRID
    !{
    ! [4.2] put ice quantities onto exchange grid
    ! (assume that ocean quantites are stored in no ice partition)
    ! (note: ex_avail is true at ice and ocean points)
    call fms_xgrid_put_to_xgrid (Ice%t_surf,      'OCN', ex_t_surf,      xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%rough_mom,   'OCN', ex_rough_mom,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%rough_heat,  'OCN', ex_rough_heat,  xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%rough_moist, 'OCN', ex_rough_moist, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%albedo,      'OCN', ex_albedo,      xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%albedo_vis_dir, 'OCN', ex_albedo_vis_dir, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%albedo_nir_dir, 'OCN', ex_albedo_nir_dir, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%albedo_vis_dif, 'OCN', ex_albedo_vis_dif, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%albedo_nir_dif, 'OCN', ex_albedo_nir_dif, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%u_surf,      'OCN', ex_u_surf,      xmap_sfc)
    call fms_xgrid_put_to_xgrid (Ice%v_surf,      'OCN', ex_v_surf,      xmap_sfc)

    tmp_open_sea        = 0.
    tmp_open_sea(:,:,1) = 1.
    call fms_xgrid_put_to_xgrid ( tmp_open_sea,  'OCN', ex_frac_open_sea,   xmap_sfc)

    do n = 1, ice%ocean_fields%num_bcs
       do m = 1, ice%ocean_fields%bc(n)%num_fields
          call fms_xgrid_put_to_xgrid (Ice%ocean_fields%bc(n)%field(m)%values, 'OCN',      &
               ex_gas_fields_ice%bc(n)%field(m)%values, xmap_sfc)
       enddo
    enddo
    !}
    

    !> ON THE EXCHANGE GRID, GENERATE DYNAMIC WET MASK ARRAY WITH VALUE OF 1.O FOR OPEN WATER
    !{
    !Generate a wet mask array on the xgrid which is:
    !       1: where there is any open water in the OCN grid cell
    !       0: where there is no  open water in the OCN grid cell, i.e., totally ice covered or land
    !This is a dynamic wet mask and particularly is not the same as a static land mask because seaice fractions change
    ! as the model runs.
    !One way to create such a mask is if 'OCN' puts an array on the xgrid which is
    ! 1 on every OCN grid cell with the 3rd index equal to 1 (ice category 1 corresponds to open water in the grid cell)
    ! 0 otherwise
    !This wet mask will be used to limit the air-sea flux exchange to areas that are not totally covered by seaice.
    sea = 0.0; sea(:,:,1) = 1.0;
    ex_seawater = 0.0
    call fms_xgrid_put_to_xgrid (sea, 'OCN', ex_seawater, xmap_sfc)
    !}
    

    !Question: Why is the above ex_seawater a dynamic mask array?
    !          From its construction it looks like a static array of 1s and 0s !
    !Answer: The xmap_sfc is dynamic and changes as the model steps because it contains updated information about
    !        seaice fractions. The updated array "ex_seawater" after the above "put" call will be 1 where there
    !        is open water even if those grid cells where previously closed by seaice.
    !        Particularly if we restrict xgrid calculations where  ex_seawater==1
    !        only cells with (partially or totally) open water contribute and totally covered cells won't contribute.

    !Not related to the above comments, it seems that the above ex_frac_open_sea could be replaced by ex_seawater
    !for code cleaning.
    !The following test does not print out anything for a fully coupled model. This asserts that
    !1. The two arrays ex_frac_open_sea and  ex_seawater are the same
    !   Their difference is that ex_frac_open_sea is a local array,
    !   but ex_seawater is a module array used outside this subroutine
    !2. The xgrid array "ex_seawater" is either 0 or 1 and nothing else, just like the "OCN" grid array "sea"
    !Test:
    !do l = 1, my_nblocks
    !   is=block_start(l)
    !   ie=block_end(l)
    !   do i = is, ie
    !      if(ex_frac_open_sea(i) /= ex_seawater(i)) &
    !           print*,'ex_frac_open_sea != ex_seawater ',ex_frac_open_sea(i),ex_seawater(i)
    !      if(ex_seawater(i) /= 0.0 .or.ex_seawater(i) /= 1.0 )&
    !           print*,'ex_seawater !1 or 0 ' ,ex_seawater(i)
    !   enddo
    !enddo

    
    !> ON THE EXCHANGE GRID, INITIALIZE CANOPY TEMPERATURE TO BE THE SURFACE TEMPERATURE
    !{
    ex_t_ca = ex_t_surf ! used to define values over the ocean
    !}
    

    !> ON THE EXCHANGE GRID, INITIALIZE CANOPY TEMPERATURE TO BE THE SURFACE TEMPERATURE
    !{
    ! [4.3] put land quantities onto exchange grid ----
    call fms_xgrid_some(xmap_sfc, ex_land, 'LND')
    !}

    
    !>  MAP LAND EXCHANGE FIELDS ONTO THE EXCHANGE GRID
    !{    
#ifndef _USE_LEGACY_LAND_

#ifdef use_AM3_physics
    if (do_forecast) then
       call fms_xgrid_put_to_xgrid_ug (Land%t_surf,     'LND', ex_t_surf_miz,  xmap_sfc)
       ex_t_ca(:) = ex_t_surf_miz(:)
    end if
#endif

    call fms_xgrid_put_to_xgrid_ug (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
    ex_rough_scale = ex_rough_mom
    call fms_xgrid_put_to_xgrid_ug(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)

    do n_gex=1,n_gex_lnd2atm
      call fms_xgrid_put_to_xgrid_ug (Land%gex_lnd2atm(:,:,n_gex),'LND', &
                                   ex_gex_lnd2atm(:,n_gex),xmap_sfc)
    end do

    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
          call fms_xgrid_put_to_xgrid_ug ( Land%tr(:,:,n), 'LND', ex_tr_surf(:,tr), xmap_sfc )
       else
          ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
          ! fluxes will be 0
       endif
    enddo
#else

#ifdef use_AM3_physics
    if (do_forecast) then
       call fms_xgrid_put_to_xgrid (Land%t_surf,     'LND', ex_t_surf_miz,  xmap_sfc)
       ex_t_ca(:) = ex_t_surf_miz(:)
    end if
#endif
    call fms_xgrid_put_to_xgrid (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
    ex_rough_scale = ex_rough_mom
    call fms_xgrid_put_to_xgrid(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)

    do n_gex=1,n_gex_lnd2atm
      call fms_xgrid_put_to_xgrid (Land%gex_lnd2atm(:,:,n_gex),'LND', &
                                   ex_gex_lnd2atm(:,n_gex),xmap_sfc)
    end do

    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
          call fms_xgrid_put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf(:,tr), xmap_sfc )
       else
          ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
          ! fluxes will be 0
       endif
    enddo
#endif

    ex_land_frac = 0.0
    call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)
    !}
    

#ifdef SCM
    if (do_specified_land) then
       if (do_specified_albedo) then
          ex_albedo = ALBEDO_OBS
          ex_albedo_vis_dir = ALBEDO_OBS
          ex_albedo_nir_dir = ALBEDO_OBS
          ex_albedo_vis_dif = ALBEDO_OBS
          ex_albedo_nir_dif = ALBEDO_OBS
       endif
       if (do_specified_tskin) then
          ex_t_surf = TSKIN
          ex_t_ca   = TSKIN
          ex_tr_surf(:,isphum) = 15.e-3
       endif
       if (do_specified_rough_leng) then
          ex_rough_mom   = ROUGH_MOM
          ex_rough_heat  = ROUGH_HEAT
          ex_rough_moist = ROUGH_HEAT
       endif
    endif
#endif

#ifdef use_AM3_physics
    if (do_forecast) then
       ex_t_surf = ex_t_surf_miz
    end if
#endif


    !> ON THE EXCHANGE GRID, COMPUTE EXPLICIT FLUXES AND TENDENCIES BY CALLING SURFACE_FLUX
    !{
    ! [5] compute explicit fluxes and tendencies at all available points ---
    call fms_xgrid_some(xmap_sfc, ex_avail)
    !$OMP parallel do default(none) shared(my_nblocks,ex_t_atm,ex_tr_atm,ex_u_atm,ex_v_atm, &
    !$OMP                                  ex_p_atm,ex_z_atm,ex_p_surf,ex_t_surf,ex_t_ca, &
    !$OMP                                  ex_tr_surf,ex_u_surf,ex_v_surf,ex_rough_mom, &
    !$OMP                                  ex_rough_heat,ex_rough_moist,ex_rough_scale,    &
    !$OMP                                  ex_gust,ex_flux_t,ex_flux_tr,ex_flux_lw, &
    !$OMP                                  ex_flux_u,ex_flux_v,ex_cd_m,ex_cd_t,ex_cd_q, &
    !$OMP                                  ex_wind,ex_u_star,ex_b_star,ex_q_star,       &
    !$OMP                                  ex_thv_atm,ex_thv_surf,                      &
    !$OMP                                  ex_dhdt_surf,ex_dedt_surf,ex_dfdtr_surf,   &
    !$OMP                                  ex_drdt_surf,ex_dhdt_atm,ex_dfdtr_atm,   &
    !$OMP                                  ex_dtaudu_atm, ex_dtaudv_atm,dt,ex_land, &
    !$OMP                                  ex_seawater,ex_avail,block_start,block_end,isphum) &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call surface_flux (&
            ex_t_atm(is:ie), ex_tr_atm(is:ie,isphum), ex_u_atm(is:ie),ex_v_atm(is:ie),ex_p_atm(is:ie),ex_z_atm(is:ie),&
            ex_p_surf(is:ie),ex_t_surf(is:ie), ex_t_ca(is:ie),  ex_tr_surf(is:ie,isphum),                       &
            ex_u_surf(is:ie), ex_v_surf(is:ie),                                           &
            ex_rough_mom(is:ie), ex_rough_heat(is:ie), ex_rough_moist(is:ie), ex_rough_scale(is:ie),    &
            ex_gust(is:ie),                                                        &
            ex_flux_t(is:ie), ex_flux_tr(is:ie,isphum), ex_flux_lw(is:ie), ex_flux_u(is:ie), ex_flux_v(is:ie),         &
            ex_cd_m(is:ie),   ex_cd_t(is:ie), ex_cd_q(is:ie),                                    &
            ex_wind(is:ie),   ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),                     &
            ex_thv_atm(is:ie),   ex_thv_surf(is:ie),                                             &
            ex_dhdt_surf(is:ie), ex_dedt_surf(is:ie), ex_dfdtr_surf(is:ie,isphum),  ex_drdt_surf(is:ie),        &
            ex_dhdt_atm(is:ie),  ex_dfdtr_atm(is:ie,isphum),  ex_dtaudu_atm(is:ie), ex_dtaudv_atm(is:ie),       &
            dt,                                                             &
            ex_land(is:ie), ex_seawater(is:ie) .gt. 0.0,  ex_avail(is:ie)            )
    enddo

#ifdef SCM
    ! Option to override surface fluxes for SCM
    if (do_specified_flux) then

       call scm_surface_flux ( &
            ex_t_atm, ex_tr_atm(:,isphum),  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
            ex_p_surf,ex_t_surf, ex_t_ca,  ex_tr_surf(:,isphum),                       &
            ex_u_surf, ex_v_surf,                                                      &
            ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale,               &
            ex_gust,                                                                   &
            ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw, ex_flux_u, ex_flux_v,         &
            ex_cd_m,   ex_cd_t, ex_cd_q,                                               &
            ex_wind,   ex_u_star, ex_b_star, ex_q_star,                                &
            ex_thv_atm, ex_thv_surf,                                                   &
            ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
            ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,       &
            dt,                                                                        &
            ex_land, ex_seawater .gt. 0.0,  ex_avail,                                    &
            ex_dhdt_surf_forland,  ex_dedt_surf_forland,  ex_dedq_surf_forland  )

    endif
#endif
    !}


    !  call mpp_clock_end(fluxClock)

    
    !> CALL MONIN_OBUKHOV_MO_PROFILE IN FMS
    !! ON THE EXCHANGE GRID, COMPUTE ZONAL AND MERIDIONAL WINDS AT THE BOUNDARY LAYER AND AT THE REFERENCE HEIGHTS
    !{
    zrefm = 10.0
    zrefh = z_ref_heat
    !      ---- optimize calculation ----
    !$OMP parallel do default(shared) private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call fms_monin_obukhov_mo_profile ( zrefm, zrefh, ex_z_atm(is:ie), ex_rough_mom(is:ie), &
            ex_rough_heat(is:ie), ex_rough_moist(is:ie),          &
            ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),        &
            ex_del_m(is:ie), ex_del_h(is:ie), ex_del_q(is:ie), ex_avail(is:ie)  )
       do i = is,ie
          ex_u10(i) = 0.
          if(ex_avail(i)) then
             ex_ref_u(i) = ex_u_surf(i) + (ex_u_atm(i)-ex_u_surf(i)) * ex_del_m(i)
             ex_ref_v(i) = ex_v_surf(i) + (ex_v_atm(i)-ex_v_surf(i)) * ex_del_m(i)
             ex_u10(i) = sqrt(ex_ref_u(i)**2 + ex_ref_v(i)**2)
          endif
       enddo
       do n = 1, ex_gas_fields_atm%num_bcs
          if (atm%fields%bc(n)%use_10m_wind_speed) then
             if (.not. ex_gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%override) then
                do i = is,ie
                   ex_gas_fields_atm%bc(n)%field(fms_coupler_ind_u10)%values(i) = ex_u10(i)
                enddo
             endif
          endif
       enddo
       !}
       
       
       !> ON THE EXCHANGE GRID, CALCULATE ATMOSPHERIC CONDUCTANCE
       !{
       do i=is,ie
          ex_con_atm(i) = ex_wind(i)*ex_cd_q(i)
       end do
       !}

       
       !> ON THE EXCHANGE GRID, COMPUTE DERIVATIVES OF TRACER FLUXES
       !{
       ! fill derivatives for all tracers
       ! F = C0*u*rho*delta_q, C0*u*rho is the same for all tracers, copy from sphum
       do tr = 1,n_exch_tr
          if (tr==isphum) cycle
          do i = is,ie
             ! slm: ex_dfdtr_surf(:,isphum) is manipulated in surface_flux: it is set to
             ! zero over the ocean, so it is not appropriate to use for other tracers.
             ! However, since flux = rho*Cd*|v|*(q_surf-q_atm), we can simply use negative
             ! dfdtr_atm for the dfdtr_surf derivative. This will break if ever the flux
             ! formulation is changed to be not symmetrical w.r.t. q_surf and q_atm, but
             ! then this whole section will have to be changed.
             ex_dfdtr_atm  (i,tr) =  ex_dfdtr_atm  (i,isphum)
             ex_dfdtr_surf (i,tr) = -ex_dfdtr_atm (i,isphum)
             ex_flux_tr    (i,tr) =  ex_dfdtr_surf(i,tr)*(ex_tr_surf(i,tr)-ex_tr_atm(i,tr))
          enddo
       enddo
    enddo
    !}
    

    ! Combine explicit ocean flux and implicit land flux of extra flux fields.


    !> ON THE EXCHANGE GRID, COMPUTE EXPLICIT FLUXES BETWEEN ATM AND OCEAN
    !{
    call atmos_ocean_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater, ex_t_surf)
    
    do n = 1, ex_gas_fluxes%num_bcs 
       if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then
          m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
          if (id_tr_mol_flux0(m) .gt. 0) then
             call fms_xgrid_get_from_xgrid (&
                  diag_atm, 'ATM', ex_gas_fluxes%bc(n)%field(fms_coupler_ind_flux0)%values(:), xmap_sfc)
             used = fms_diag_send_data ( id_tr_mol_flux0(m), diag_atm, Time )
          end if
       end if
    end do

    ! The following statement is a concise version of what's following and worth
    ! looking into in the future.
    ! ex_flux_tr(:,itracer) = ex_gas_fluxes%bc(itracer_ocn)%field(fms_coupler_ind_flux)%values(:)
    ! where(ex_seawater.gt.0) ex_flux_tr(:,itracer) = F_ocn
    !$OMP parallel do default(shared) private(is,ie,m,tr_units,tr_name)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do n = 1, ex_gas_fluxes%num_bcs
          if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then
             m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
             call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, ex_gas_fluxes%bc(n)%atm_tr_index, tr_name, &
                                                       units=tr_units)
             do i = is,ie
                if (ex_land(i)) cycle  ! over land, don't do anything
                ! on ocean or ice cells, flux is explicit therefore we zero derivatives.
                ex_dfdtr_atm(i,m)  = 0.0
                ex_dfdtr_surf(i,m) = 0.0
                if (ex_seawater(i)>0.0) then
                   if (fms_mpp_lowercase(trim(tr_units)).eq."vmr") then
                      ! in mol/m2/s but from land model it should be in vmr * kg/m2/s
                      ! This was converting to dry vmr (as opposed to ambient vmr)
                      !    ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i)  &
                      !         * WTMAIR*1.0e-3 &
                      !         / (1.-ex_tr_atm(i,isphum))
                      !
                      ! vmr * kg/m2/s = mol(X)/[m2 s] [1/mol(air) * mol(air)] * WTM(air)
                      !
                      ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) &
                           * 1.0e-3*WTMAIR*WTMH2O/((1.-ex_tr_atm(i,isphum))*WTMH2O+ex_tr_atm(i,isphum)*WTMAIR)

                      !end if
                   else
                   ! convert to kg co2/m2/sec for atm
                   ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(fms_coupler_ind_flux)%values(i) * &
                                        ex_gas_fluxes%bc(n)%mol_wt * 1.0e-03
                   end if
                else
                   ex_flux_tr(i,m) = 0.0 ! pure ice exchange cell
                endif
             enddo
          endif
       enddo 
    enddo
    !}
    

    !> OVERRIDE LAND AND ICE TRACER FLUXES.  DATA WILL BE OVERWRITTEN ONLY IF FIELD IS SPECIFIED IN DATA_TABLE
    !{
    ! [5.2] override tracer fluxes and derivatives
    do tr = 1,n_exch_tr
       if( tr_table(tr)%atm == NO_TRACER ) cycle ! it should never happen, though

       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       ! [5.2.1] override tracer flux. Note that "sea" and "diag_land" are repeatedly used
       ! as temporary storage for the values we are overriding fluxes and derivative with,
       ! over ocean and land respectively
#ifndef _USE_LEGACY_LAND_
       call fms_data_override_ug ( 'LND', 'ex_flux_'//trim(tr_name), diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_flux_'//trim(tr_name), sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_flux_tr(:,tr), xmap_sfc )
       ! [5.2.2] override derivative of flux wrt surface concentration
       call fms_data_override_ug ( 'LND', 'ex_dfd'//trim(tr_name)//'_surf', diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_dfdtr_surf(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_surf', sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dfdtr_surf(:,tr), xmap_sfc )
       ! [5.2.3] override derivative of flux wrt atmospheric concentration
       call fms_data_override_ug ( 'LND', 'ex_dfd'//trim(tr_name)//'_atm', diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_dfdtr_atm(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_atm', sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dfdtr_atm(:,tr), xmap_sfc )
    enddo

    ! [5.3] override flux and derivatives for sensible heat flux
    ! [5.3.1] override flux
    call fms_data_override_ug ( 'LND', 'ex_flux_t', diag_land, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_flux_t, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_flux_t', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_flux_t, xmap_sfc )
    ! [5.3.2] override derivative of flux wrt near-surface temperature
    call fms_data_override_ug ( 'LND', 'ex_dhdt_surf', diag_land, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_dhdt_surf, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_dhdt_surf', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dhdt_surf, xmap_sfc )
    ! [5.3.3] override derivative of flux wrt atmospheric temperature
    call fms_data_override_ug ( 'LND', 'ex_dhdt_atm', diag_land, Time,override=used )
    if (used) call fms_xgrid_put_to_xgrid_ug ( diag_land, 'LND', ex_dhdt_atm, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_dhdt_atm', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dhdt_atm, xmap_sfc )
#else
       call fms_data_override ( 'LND', 'ex_flux_'//trim(tr_name), diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_flux_'//trim(tr_name), sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_flux_tr(:,tr), xmap_sfc )
       ! [5.2.2] override derivative of flux wrt surface concentration
       call fms_data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_surf', diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_dfdtr_surf(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_surf', sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dfdtr_surf(:,tr), xmap_sfc )
       ! [5.2.3] override derivative of flux wrt atmospheric concentration
       call fms_data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_atm', diag_land, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_dfdtr_atm(:,tr), xmap_sfc )
       call fms_data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_atm', sea, Time, override=used )
       if(used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dfdtr_atm(:,tr), xmap_sfc )
    enddo

    ! [5.3] override flux and derivatives for sensible heat flux
    ! [5.3.1] override flux
    call fms_data_override ( 'LND', 'ex_flux_t', diag_land, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_flux_t, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_flux_t', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_flux_t, xmap_sfc )
    ! [5.3.2] override derivative of flux wrt near-surface temperature
    call fms_data_override ( 'LND', 'ex_dhdt_surf', diag_land, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_dhdt_surf, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_dhdt_surf', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dhdt_surf, xmap_sfc )
    ! [5.3.3] override derivative of flux wrt atmospheric temperature
    call fms_data_override ( 'LND', 'ex_dhdt_atm', diag_land, Time,override=used )
    if (used) call fms_xgrid_put_to_xgrid ( diag_land, 'LND', ex_dhdt_atm, xmap_sfc )
    call fms_data_override ( 'ICE', 'ex_dhdt_atm', sea, Time, override=used )
    if (used) call fms_xgrid_put_to_xgrid ( sea, 'OCN', ex_dhdt_atm, xmap_sfc )
#endif
    !}
    

    ! NB: names of the override fields are constructed using tracer name and certain
    ! prefixes / suffixes. For example, for the tracer named "sphum" (specific humidity) they will be:
    ! "ex_flux_sphum", "ex_dfdsphum_surf", and "ex_dfdsphum_atm".
    !
    ! For sensible heat flux names are "ex_flux_t", "ex_dhdt_surf", and "ex_dhdt_atm";
    ! despite the name those are actually in energy units, W/m2, W/(m2 degK), and
    ! W/(m2 degK) respectively


    !> ON THE EXCHANGE GRID, COMPUTE T_SURF**4
    !! NOTE, TO COMPUTE LONGWAVE RADIATION FLUXES MORE ACCURATELY, T_SURF**4 IS
    !! PUT ON THE EXCHANGE GRID RATHER THAN T_SURF IN ORDER TO ACCOUNT FOR THE NONLINEARITY IN
    !! THE STEFAN-BOLTZMANN LAW WHERE LONGWAVE_FLUX = STEFAN_BOLTZMANN _CONSTANT * T**4.
    !! ON THE EXCHANGE GRID, AS QUANTITIES ARE REMAPPED, FIELDS ARE AREA-WEIGHTED (AVERAGED)
    !! SUCH THAT OUTPUT_TEMPERATURE = SUM(INPUT_TEMPERATURE * (XGRID_AREA)/(INPUT_GRID_AREA))
    !! WHERE THE SUM IS OVER ALL XGRID CELLS THAT OVERLAP WITH THE OUTPUT CELL. 
    !! BECAUSE OF THIS WEIGHTING, THE COMPUTED FLUX WOULD DIFFER FROM USING <T**4> VS <T>**4
    !! (WHERE <> IS USED TO DENOTE AVERAGING)
    !{
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_avail,  &
    !$OMP                                  ex_drag_q,ex_wind,ex_cd_q,ex_t_surf4,ex_t_surf ) &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          if(ex_avail(i)) ex_drag_q(i) = ex_wind(i)*ex_cd_q(i)
          ! [6] get mean quantities on atmosphere grid
          ! [6.1] compute t surf for radiation
          ex_t_surf4(i) = ex_t_surf(i) ** 4
       enddo
    enddo
    !}

    !> MAP QUANTITIES FROM THE EXCHANGE GRID TO THE ATM GRID
    !{
    ! [6.2] put relevant quantities onto atmospheric boundary
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%t_ocean,           'ATM',   &
         ex_t_surf ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%frac_open_sea,'ATM',ex_frac_open_sea, xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir,    'ATM',   &
         ex_albedo_vis_dir   ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir,    'ATM',   &
         ex_albedo_nir_dir   ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif,    'ATM',   &
         ex_albedo_vis_dif   ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif,    'ATM',   &
         ex_albedo_nir_dif   ,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%rough_heat,'ATM',           &
         ex_rough_heat, xmap_sfc, complete=.false.) !kgao
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc, complete=.false.)

    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc, complete=.true.)

    do n_gex=1,n_gex_lnd2atm
      call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%gex_lnd2atm(:,:,n_gex),  'ATM', ex_gex_lnd2atm(:,n_gex), &
                                     xmap_sfc, complete=.false.)
    end do

    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%u_ref, 'ATM', ex_ref_u     , xmap_sfc, complete=.false.) !bqx
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%v_ref, 'ATM', ex_ref_v     , xmap_sfc, complete=.true.) !bqx

! for shield+mom6 coupling; used by shield pbl schemes (am5 with tke-edmf should do the same)
#ifndef use_AM3_physics
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%shflx,'ATM', ex_flux_t    , xmap_sfc, complete=.false.)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%lhflx,'ATM', ex_flux_tr(:,isphum), xmap_sfc, complete=.true.)
#endif
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%wind,      'ATM', ex_wind      , xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%thv_atm,   'ATM', ex_thv_atm   , xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%thv_surf,  'ATM', ex_thv_surf  , xmap_sfc)

#ifdef use_AM3_physics
    if (do_forecast) then
       call fms_xgrid_get_from_xgrid (Ice%t_surf, 'OCN', ex_t_surf,  xmap_sfc)
    end if
#endif
    !}


    !> ON THE ATM GRID, COMPUTE T**0.25
    !{
    call fms_mpp_domains_get_compute_domain( Atm%domain, isc, iec, jsc, jec )
    !$OMP parallel do default(none) shared(isc,iec,jsc,jec,Land_Ice_Atmos_Boundary ) &
    !$OMP                          private(is,ie)
    do j = jsc, jec
       do i = isc, iec
          Land_Ice_Atmos_Boundary%t(i,j) = Land_Ice_Atmos_Boundary%t(i,j) ** 0.25
       enddo
    enddo
    !}

    
    !> DATA OVERRIDE ATMOSPHERIC QUANTITIES.
    !! DATA_OVERRIDE WILL ONLY OVERWRITE IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    ! fms_data_override calls moved here from coupler_main
    call fms_data_override('ATM', 't',         Land_Ice_Atmos_Boundary%t,         Time)
    call fms_data_override('ATM', 'albedo',    Land_Ice_Atmos_Boundary%albedo,    Time)

    call fms_data_override('ATM', 'albedo_vis_dir',    Land_Ice_Atmos_Boundary%albedo_vis_dir,    Time)
    call fms_data_override('ATM', 'albedo_nir_dir',    Land_Ice_Atmos_Boundary%albedo_nir_dir,    Time)
    call fms_data_override('ATM', 'albedo_vis_dif',    Land_Ice_Atmos_Boundary%albedo_vis_dif,    Time)
    call fms_data_override('ATM', 'albedo_nir_dif',    Land_Ice_Atmos_Boundary%albedo_nir_dif,    Time)
    call fms_data_override('ATM', 'land_frac', Land_Ice_Atmos_Boundary%land_frac, Time)
    call fms_data_override('ATM', 'dt_t',      Land_Ice_Atmos_Boundary%dt_t,      Time)
    do tr=1,n_atm_tr
       call fms_tracer_manager_get_tracer_names(MODEL_ATMOS, tr, tr_name)
       call fms_data_override('ATM', 'dt_'//trim(tr_name), Land_Ice_Atmos_Boundary%dt_tr(:,:,tr), Time)
    enddo
    call fms_data_override('ATM', 'u_flux',    Land_Ice_Atmos_Boundary%u_flux,    Time)
    call fms_data_override('ATM', 'v_flux',    Land_Ice_Atmos_Boundary%v_flux,    Time)
    call fms_data_override('ATM', 'dtaudu',    Land_Ice_Atmos_Boundary%dtaudu,    Time)
    call fms_data_override('ATM', 'dtaudv',    Land_Ice_Atmos_Boundary%dtaudv,    Time)
    call fms_data_override('ATM', 'u_star',    Land_Ice_Atmos_Boundary%u_star,    Time)
    call fms_data_override('ATM', 'b_star',    Land_Ice_Atmos_Boundary%b_star,    Time)
    ! call fms_data_override('ATM', 'q_star',    Land_Ice_Atmos_Boundary%q_star,    Time)
    call fms_data_override('ATM', 'rough_mom', Land_Ice_Atmos_Boundary%rough_mom, Time)
    !}

    
    !> ON THE EXCHANGE GRID, INITIALIZE ARRAYS FOR FIXING THE ALBEDO
    !
    !!  STILL NEEDED   ????
    !! IS THIS CORRECT ??
    allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )

    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end, &
    !$OMP ex_albedo_vis_dir_fix,ex_albedo_nir_dir_fix, ex_albedo_vis_dif_fix,ex_albedo_nir_dif_fix ) &
    !$OMP private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_vis_dir_fix(i) = 0.
          ex_albedo_nir_dir_fix(i) = 0.
          ex_albedo_vis_dif_fix(i) = 0.
          ex_albedo_nir_dif_fix(i) = 0.
       enddo
    enddo
    !}

    !> MAP ALBEDO FIELDS TO THE EXCHANGE GRID
    !{
    call fms_xgrid_put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir, 'ATM',  &
         ex_albedo_vis_dir_fix, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir, 'ATM', &
         ex_albedo_nir_dir_fix, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif, 'ATM',   &
         ex_albedo_vis_dif_fix, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif, 'ATM',  &
         ex_albedo_nir_dif_fix, xmap_sfc, complete=.true.)
    !}


    !> ON THE EXCHANGE GRID, COMPUTE THE ALBEDO FIXING FACTORS
    !{
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end, &
    !$OMP                                  ex_albedo,ex_albedo_vis_dir_fix,ex_albedo_vis_dir, &
    !$OMP                                  ex_albedo_nir_dir,ex_albedo_nir_dir_fix,           &
    !$OMP                                  ex_albedo_vis_dif_fix,ex_albedo_vis_dif,           &
    !$OMP                                  ex_albedo_nir_dif_fix,ex_albedo_nir_dif)           &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_vis_dir_fix(i) = (1.0-ex_albedo_vis_dir(i)) / (1.0-ex_albedo_vis_dir_fix(i))
          ex_albedo_nir_dir_fix(i) = (1.0-ex_albedo_nir_dir(i)) / (1.0-ex_albedo_nir_dir_fix(i))
          ex_albedo_vis_dif_fix(i) = (1.0-ex_albedo_vis_dif(i)) / (1.0-ex_albedo_vis_dif_fix(i))
          ex_albedo_nir_dif_fix(i) = (1.0-ex_albedo_nir_dif(i)) / (1.0-ex_albedo_nir_dif_fix(i))
       enddo
    enddo
    !}
    

#ifdef SCM
    if (do_specified_albedo .and. do_specified_land) then
       ex_albedo_vis_dir_fix = 1.
       ex_albedo_vis_dif_fix = 1.
       ex_albedo_nir_dir_fix = 1.
       ex_albedo_nir_dif_fix = 1.
    endif
#endif

    !=======================================================================
    ! [7] diagnostics section

    
    !> SAVE STATIC FIELDS.  THE STATIC FIELDS WILL BE SAVED ONLY THE FIRST TIME THIS SUBROUTINE IS CALLED
    !! IF FIRST_STATIC = .TRUE.
    !{
    !------- save static fields first time only ------
    if (first_static) then

       !------- land fraction ------
       if ( id_land_mask > 0 ) then
          used = fms_diag_send_data ( id_land_mask, Land_Ice_Atmos_Boundary%land_frac, Time )
       endif
       if ( id_sftlf > 0 ) then
          used = fms_diag_send_data ( id_sftlf, Land_Ice_Atmos_Boundary%land_frac, Time )
       endif
       ! near-surface heights
       if ( id_height2m  > 0) used = fms_diag_send_data ( id_height2m, z_ref_heat, Time )
       if ( id_height10m > 0) used = fms_diag_send_data ( id_height10m, z_ref_mom, Time )

       first_static = .false.
    endif
    !}

    
    !> MAP ATMOSPHERIC DATA FROM THE EXCHANGE GRID TO THE ATM GRID AND SEND DATA TO DIAG_MANAGER BUFFER.
    !! NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FOUND IN THE DIAG_TABLE
    !{    
    !------- Atm fields -----------
    do n = 1, Atm%fields%num_bcs
       do m = 1, Atm%fields%bc(n)%num_fields
          if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then
             if (atm%fields%bc(n)%use_10m_wind_speed .and. m .eq. fms_coupler_ind_u10 .and. &
                 .not. Atm%fields%bc(n)%field(m)%override) then
                call fms_xgrid_get_from_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',     &
                     ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc)
             endif
             if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then
                used = fms_diag_send_data(Atm%fields%bc(n)%field(m)%id_diag, Atm%fields%bc(n)%field(m)%values, Time )
             endif
          endif
       enddo
    enddo

    !------- drag coeff moisture -----------
    if ( id_wind > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_wind, xmap_sfc)
       used = fms_diag_send_data ( id_wind, diag_atm, Time )
    endif
    !------- drag coeff moisture -----------
    if ( id_drag_moist > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_cd_q, xmap_sfc)
       used = fms_diag_send_data ( id_drag_moist, diag_atm, Time )
    endif

    !------- drag coeff heat -----------
    if ( id_drag_heat > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_cd_t, xmap_sfc)
       used = fms_diag_send_data ( id_drag_heat, diag_atm, Time )
    endif

    !------- drag coeff momemtum -----------
    if ( id_drag_mom > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_cd_m, xmap_sfc)
       used = fms_diag_send_data ( id_drag_mom, diag_atm, Time )
    endif

    !------- roughness moisture -----------
    if ( id_rough_moist > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_rough_moist, xmap_sfc)
       used = fms_diag_send_data ( id_rough_moist, diag_atm, Time )
    endif

    !------- roughness heat -----------
    if ( id_rough_heat > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_rough_heat, xmap_sfc)
       used = fms_diag_send_data ( id_rough_heat, diag_atm, Time )
    endif

    !------- roughness momemtum -----------
    used = fms_diag_send_data ( id_rough_mom, Land_Ice_Atmos_Boundary%rough_mom, Time )

    !------- friction velocity -----------
    used = fms_diag_send_data ( id_u_star, Land_Ice_Atmos_Boundary%u_star, Time )

    !------- bouyancy -----------
    used = fms_diag_send_data ( id_b_star, Land_Ice_Atmos_Boundary%b_star, Time )

    !------- moisture scale -----------
    used = fms_diag_send_data ( id_q_star, Land_Ice_Atmos_Boundary%q_star, Time )

    !------- surf and atm virtual potential temperature -----------
    used = fms_diag_send_data ( id_thv_atm,  Land_Ice_Atmos_Boundary%thv_atm,  Time )
    used = fms_diag_send_data ( id_thv_surf, Land_Ice_Atmos_Boundary%thv_surf, Time )

    !-----------------------------------------------------------------------
    !------ diagnostics for fields at bottom atmospheric level ------

    if ( id_t_atm > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_t_atm, xmap_sfc)
       used = fms_diag_send_data ( id_t_atm, diag_atm, Time )
    endif

    if ( id_u_atm > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_u_atm, xmap_sfc)
       used = fms_diag_send_data ( id_u_atm, diag_atm, Time )
    endif

    if ( id_v_atm > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_v_atm, xmap_sfc)
       used = fms_diag_send_data ( id_v_atm, diag_atm, Time )
    endif

    do tr = 1,n_exch_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       if ( id_tr_atm(tr) > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_tr_atm(:,tr), xmap_sfc)
          used = fms_diag_send_data ( id_tr_atm(tr), diag_atm, Time )
       endif
       !!jgj: add dryvmr co2_atm
       ! - slm Mar 25 2010: moved to resolve interdependence of diagnostic fields
       if ( id_co2_atm_dvmr > 0 .and. fms_mpp_lowercase(trim(tr_name))=='co2') then
          ex_co2_atm_dvmr = (ex_tr_atm(:,tr) / (1.0 - ex_tr_atm(:,isphum))) * WTMAIR/WTMCO2
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_co2_atm_dvmr, xmap_sfc)
          used = fms_diag_send_data ( id_co2_atm_dvmr, diag_atm, Time )
       endif
    enddo

    ! - slm, Mar 25, 2002
    if ( id_p_atm > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_p_atm, xmap_sfc)
       used = fms_diag_send_data ( id_p_atm, diag_atm, Time )
    endif
    if ( id_z_atm > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_z_atm, xmap_sfc)
       used = fms_diag_send_data ( id_z_atm, diag_atm, Time )
    endif
    if ( id_gust > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_gust, xmap_sfc)
       used = fms_diag_send_data ( id_gust, diag_atm, Time )
    endif

    ! - bw, Sep 17, 2007
    if ( id_slp > 0 .or. id_psl > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_slp, xmap_sfc)
       if ( id_slp > 0 ) used = fms_diag_send_data ( id_slp, diag_atm, Time )
       if ( id_psl > 0 ) used = fms_diag_send_data ( id_psl, diag_atm, Time )
    endif

    if ( id_t_ocean > 0 ) then
       used = fms_diag_send_data ( id_t_ocean, Land_Ice_Atmos_Boundary%t_ocean, Time )
    endif
    !}


    !>COMPUTE DIAGNOSTIC FIELDS AT REFERENCE LEVELS WITH FMS_MONIN_OBUHKOV_MO_PROFILE
    !{
    !--------- diagnostics for fields at reference level ---------
    !cjg
    !  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
    !       id_u_ref > 0 .or. id_v_ref  > 0 .or. id_wind_ref > 0 .or. &
    !       id_q_ref > 0 .or. id_q_ref_land > 0 .or. &
    !       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
    !       id_rh_ref_cmip >0 .or. &
    !       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then

    zrefm = z_ref_mom
    zrefh = z_ref_heat
    !      ---- optimize calculation ----
    !cjg     if ( id_t_ref <= 0 ) zrefh = zrefm

    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,zrefm,zrefh,ex_z_atm, &
    !$OMP                                  ex_rough_mom,ex_rough_heat,ex_rough_moist,ex_u_star,   &
    !$OMP                                  ex_b_star,ex_q_star,ex_del_m,ex_del_h,ex_del_q,        &
    !$OMP                                  ex_tr_ref,n_exch_tr,id_tr_ref,id_tr_ref_land,          &
    !$OMP                                  ex_tr_con_ref, id_tr_con_ref, id_tr_con_ref_land,      &
    !$OMP                                  ex_tr_con_atm, id_tr_con_atm, id_tr_con_atm_land,      &
    !$OMP                                  ex_avail,ex_ref,ex_tr_surf,ex_tr_atm,isphum,           &
    !$OMP                                  ex_flux_tr,ex_t_atm,ex_p_surf)                         &
    !$OMP                          private(is,ie,rho)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call fms_monin_obukhov_mo_profile ( zrefm, zrefh, ex_z_atm(is:ie),   ex_rough_mom(is:ie), &
            ex_rough_heat(is:ie), ex_rough_moist(is:ie),          &
            ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),        &
            ex_del_m(is:ie), ex_del_h(is:ie), ex_del_q(is:ie), ex_avail(is:ie)  )

       !    ------- reference relative humidity -----------
       !cjg     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
       !cjg          id_rh_ref_cmip > 0 .or. &
       !cjg          id_q_ref > 0 .or. id_q_ref_land >0 ) then
       do i = is,ie
          ex_ref(i) = 1.0e-06
          ex_tr_ref(i,:) = 1.e-20
          if (ex_avail(i) .and. ex_rough_moist(i) > 0.) then
               ex_ref(i)   = ex_tr_surf(i,isphum) + (ex_tr_atm(i,isphum)-ex_tr_surf(i,isphum)) * ex_del_q(i)
               rho         = ex_p_surf(i)/(rdgas * ex_t_atm(i)*(1.0+d608*ex_tr_atm(i,isphum)))
               do tr=1,n_exch_tr
                  if (id_tr_ref(tr).gt.0             &
                       .or. id_tr_ref_land(tr).gt.0  &
                       .or. id_tr_con_ref(tr).gt.0   &
                       .or. id_tr_con_ref_land(tr).gt.0) then
                     ex_tr_ref(i,tr) = ex_tr_surf(i,tr) + (ex_tr_atm(i,tr)-ex_tr_surf(i,tr)) * ex_del_q(i)
                     ex_tr_con_ref(i,tr) = -ex_flux_tr(i,tr)/max(ex_tr_ref(i,tr)*rho,epsln)
                  end if
                  if (id_tr_con_atm(tr).gt.0 .or. id_tr_con_atm_land(tr).gt.0) then
                     ex_tr_con_atm(i,tr) = -ex_flux_tr(i,tr)/max(ex_tr_atm(i,tr)*rho,epsln)
                  end if
               end do
            end if
       enddo
    enddo
    !}


    !> MAP LAND AND ADDITIONAL ATMOSPHERIC FIELDS FROM THE EXCHANGE GRID TO THE COMPONENT GRID AND
    !! SEND DATA TO DIAG_MANAGER BUFFER FOR DIAGNOSTIC OUTPUT
    !! NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FOUND IN THE DIAG_TABLE
    !{
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%q_ref, 'ATM', ex_ref,   xmap_sfc)  ! cjg
    if(id_q_ref > 0) then
       used = fms_diag_send_data(id_q_ref,Land_Ice_Atmos_Boundary%q_ref,Time)
    endif
    if(id_huss > 0) then
       used = fms_diag_send_data(id_huss,Land_Ice_Atmos_Boundary%q_ref,Time)
    endif
    if(id_q_ref_land > 0 .or.id_hussLut_land > 0) then
       !duplicate send_tile_data. We may remove id_q_ref_land in the future.
#ifndef _USE_LEGACY_LAND_
       call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_ref, xmap_sfc)
       call send_tile_data (id_q_ref_land, diag_land)
       call send_tile_data (id_hussLut_land, diag_land)
#else
       call fms_xgrid_get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
       used = fms_diag_send_tile_averaged_data(id_q_ref_land, diag_land, &
            Land%tile_size, Time, mask=Land%mask)
#endif
    endif

    do tr=1,n_exch_tr
       if (id_tr_ref(tr)>0) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_tr_ref(:,tr),   xmap_sfc)  ! cjg
          if(id_tr_ref(tr) > 0) then
             used = fms_diag_send_data(id_tr_ref(tr),diag_atm,Time)
          endif
       end if

       if(id_tr_ref_land(tr) > 0) then
#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_tr_ref(:,tr), xmap_sfc)
          !duplicate send_tile_data. We may remove id_q_ref_land in the future.
          call send_tile_data (id_tr_ref_land(tr), diag_land)
#else
          call fms_xgrid_get_from_xgrid (diag_land, 'LND', ex_tr_ref(:,tr), xmap_sfc)
          used = fms_diag_send_tile_averaged_data(id_tr_ref_land(tr), diag_land, &
               Land%tile_size, Time, mask=Land%mask)
#endif
       endif
    end do

    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_t_ref,ex_avail, &
    !$OMP                                  ex_rough_heat,ex_t_ca,ex_t_atm,ex_p_surf,ex_qs_ref,ex_del_h,      &
    !$OMP                                  ex_ref,ex_qs_ref_cmip,ex_ref2 ) &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is,ie
          ex_t_ref(i) = 200.
          if ( ex_avail(i) .and. ex_rough_heat(i) > 0. ) &
               ex_t_ref(i) = ex_t_ca(i) + (ex_t_atm(i)-ex_t_ca(i)) * ex_del_h(i)
       enddo
       call fms_sat_vapor_pres_compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref(is:ie), q = ex_ref(is:ie))
       call fms_sat_vapor_pres_compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref_cmip(is:ie),  &
            q = ex_ref(is:ie), es_over_liq_and_ice = .true.)
       do i = is,ie
          if(ex_avail(i)) then
             ! remove cap on relative humidity -- this mod requested by cjg, ljd
             !RSH    ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)
             ex_ref2(i)   = 100.*ex_ref(i)/ex_qs_ref_cmip(i)
             ex_ref(i)    = 100.*ex_ref(i)/ex_qs_ref(i)
          endif
       enddo
    enddo

    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%t_ref, 'ATM', ex_t_ref, xmap_sfc)  ! cjg

    if ( id_rh_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
       call fms_xgrid_get_from_xgrid_ug (diag_land,'LND', ex_ref, xmap_sfc)
       call send_tile_data (id_rh_ref_land, diag_land)
#else
       call fms_xgrid_get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
       used = fms_diag_send_tile_averaged_data ( id_rh_ref_land, diag_land, &
            Land%tile_size, Time, mask = Land%mask )
#endif
    endif

    if(id_rh_ref > 0) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
       used = fms_diag_send_data ( id_rh_ref, diag_atm, Time )
    endif

    if(id_rh_ref_cmip > 0 .or. id_hurs > 0 .or. id_rhs > 0) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref2, xmap_sfc)
       if (id_rh_ref_cmip > 0) used = fms_diag_send_data ( id_rh_ref_cmip, diag_atm, Time )
       if (id_hurs > 0)        used = fms_diag_send_data ( id_hurs, diag_atm, Time )
       if (id_rhs  > 0)        used = fms_diag_send_data ( id_rhs,  diag_atm, Time )
    endif
    !cjg  endif

    !    ------- reference temp -----------
#ifdef use_AM3_physics
    if ( id_t_ref > 0 .or. id_t_ref_land > 0 .or. id_tasLut_land > 0 ) then
       where (ex_avail) &
            ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
       if (id_t_ref_land > 0.or.id_tasLut_land > 0) then
#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_get_from_xgrid_ug(diag_land, 'LND', ex_ref, xmap_sfc)
          if (id_t_ref_land > 0)  call send_tile_data (id_t_ref_land, diag_land)
          if (id_tasLut_land > 0) call send_tile_data (id_tasLut_land, diag_land)
#else
          call fms_xgrid_get_from_xgrid(diag_land, 'LND', ex_ref, xmap_sfc)
          if (id_t_ref_land > 0) used = fms_diag_send_tile_averaged_data ( id_t_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif

       if ( id_t_ref > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          used = fms_diag_send_data ( id_t_ref, diag_atm, Time )
       endif
    endif
#else
    where (ex_avail) &
         ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
    if (id_t_ref_land > 0 .or. id_tasLut_land > 0 .or. id_tasl_g > 0) then
       ! t_ref diagnostic at land points only
#ifndef _USE_LEGACY_LAND_
       call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_ref, xmap_sfc)
       if (id_t_ref_land > 0)  call send_tile_data (id_t_ref_land, diag_land)
       if (id_tasLut_land > 0) call send_tile_data (id_tasLut_land, diag_land)
       if (id_tasl_g > 0) then
         used = send_global_land_diag ( get_global_diag_field_id(id_tasl_g), &
                                diag_land, Time, Land%tile_size, Land%mask, Land )
       endif
#else
       call fms_xgrid_get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
       if (id_t_ref_land > 0) used = fms_diag_send_tile_averaged_data ( id_t_ref_land, diag_land, &
            Land%tile_size, Time, mask = Land%mask )
#endif
    endif
    ! t_ref diagnostic at all atmos points
    call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
    if ( id_t_ref > 0 ) used = fms_diag_send_data ( id_t_ref, diag_atm, Time )
    if ( id_tas > 0 )   used = fms_diag_send_data ( id_tas, diag_atm, Time )
    call fms_sum_diag_integral_field ('t_ref',  diag_atm)
    if ( id_tas_g > 0 )  used = send_global_diag ( id_tas_g, diag_atm, Time )
#endif

    !    ------- reference u comp -----------
    if ( id_u_ref > 0 .or. id_u_ref_land > 0 .or. id_uas > 0) then
       where (ex_avail) &
            ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
       if ( id_u_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_get_from_xgrid_ug ( diag_land, 'LND', ex_ref, xmap_sfc )
          call send_tile_data ( id_u_ref_land, diag_land )
#else
          call fms_xgrid_get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
          used = fms_diag_send_tile_averaged_data ( id_u_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif
       if ( id_u_ref > 0 .or. id_uas > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          if ( id_u_ref > 0 ) used = fms_diag_send_data ( id_u_ref, diag_atm, Time )
          if ( id_uas > 0 )   used = fms_diag_send_data ( id_uas, diag_atm, Time )
       endif
    endif

    !    ------- reference v comp -----------
    if ( id_v_ref > 0 .or. id_v_ref_land > 0 .or. id_vas > 0 ) then
       where (ex_avail) &
            ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
       if ( id_v_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_get_from_xgrid_ug ( diag_land, 'LND', ex_ref, xmap_sfc )
          call send_tile_data ( id_v_ref_land, diag_land )
#else
          call fms_xgrid_get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
          used = fms_diag_send_tile_averaged_data ( id_v_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif
       if ( id_v_ref > 0 .or. id_vas > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          if ( id_v_ref > 0 ) used = fms_diag_send_data ( id_v_ref, diag_atm, Time )
          if ( id_vas   > 0 ) used = fms_diag_send_data ( id_vas, diag_atm, Time )
       endif
    endif

    !    ------- reference-level absolute wind -----------
    if ( id_wind_ref > 0 .or. id_sfcWind > 0 ) then
       where (ex_avail) &
            ex_ref = sqrt((ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m)**2 &
            +(ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m)**2)
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
       if ( id_wind_ref > 0 ) used = fms_diag_send_data ( id_wind_ref, diag_atm, Time )
       if ( id_sfcWind  > 0 ) used = fms_diag_send_data ( id_sfcWind, diag_atm, Time )
    endif

    !    ------- interp factor for heat ------
    if ( id_del_h > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_del_h, xmap_sfc)
       used = fms_diag_send_data ( id_del_h, diag_atm, Time )
    endif

    !    ------- interp factor for momentum ------
    if ( id_del_m > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_del_m, xmap_sfc)
       used = fms_diag_send_data ( id_del_m, diag_atm, Time )
    endif

    !    ------- interp factor for moisture ------
    if ( id_del_q > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_del_q, xmap_sfc)
       used = fms_diag_send_data ( id_del_q, diag_atm, Time )
    endif

    !cjg  endif
    ! topographic roughness scale
    if(id_rough_scale>0) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM',&
            (log(ex_z_atm/ex_rough_mom+1.0)/log(ex_z_atm/ex_rough_scale+1.0))**2, xmap_sfc)
       used = fms_diag_send_data(id_rough_scale, diag_atm, Time)
    endif
    !}


    !> END CLOCKS FOR PROFILING
    !{
    call fms_mpp_clock_end(sfcClock)
    call fms_mpp_clock_end(cplClock)
    !}

  end subroutine sfc_boundary_layer

  !> Subroutine flux_down_from_atmos corrects for the implicit treatment of atmospheric diffisuve fluxes
  !! in flux exchange from atm to land and ice
  subroutine flux_down_from_atmos (Time, Atm, Land, Ice, Atmos_boundary, Land_boundary, Ice_boundary )

    !> current model time
    type(FmsTime_type), intent(in) :: Time
    !> derived data type holding atmosphere boundary data
    type(atmos_data_type), intent(inout) :: Atm 
    !> derived data type holding land boundary data
    type(land_data_type), intent(in) :: Land
    !> derived data type holding ice boundary dat
    type(ice_data_type), intent(in) :: Ice
    !< derived data type holding properties and fluxes passed from exchange grid to atmosphere land and ice
    type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
    !> derived data type holding properties and fluxes passed from atmosphere to land
    type(atmos_land_boundary_type), intent(inout):: Land_boundary
    !> derived data type holding properties and fluxes passed from atmosphere to ice
    type(atmos_ice_boundary_type), intent(inout):: Ice_boundary

    ! the following variables are on the exchange grid
    real, dimension(n_xgrid_sfc) :: &
         ex_flux_sw, & ! net sfc shortwave radiation flux
         ex_flux_lwd, & ! longwave radiation downward flux at surface
         ex_flux_sw_dir, & ! direct shortwave radiation flux
         ex_flux_sw_dif, & ! downward diffuse visible shortwave radiation flux ar the surface
         ex_flux_sw_down_vis_dir, & ! downward direct visible shortwave radiation flux at the surface
         ex_flux_sw_down_total_dir,  & ! downward direct total shortwave radiation flux at the surface
         ex_flux_sw_down_vis_dif, & ! downward diffuse visible shortwave radiation flux at the surface
         ex_flux_sw_down_total_dif,  & ! downward diffusive total shortwave radiation flux at the surface
         ex_flux_sw_vis, & ! net visible shortwave radiation flux at the surface
         ex_flux_sw_vis_dir, & ! net direct visible shortwave radiation flux at the surface
         ex_flux_sw_vis_dif, & ! net diffuse visible shortwave radiation flux at the surface
         ex_lprec, & ! liquid precipitation
         ex_fprec, & ! frozen precitipation
         ex_tprec, & ! temperature of precipitation, currently equal to atm temperature
         ex_u_star_smooth, & ! friction velocity correction
#ifdef use_AM3_physics
         ex_coszen ! cosine of the zenith angle
#else
         ex_coszen, & ! cosine of the zenith angle
         ex_setl_flux, & ! tracer sedimentation flux from the lowest atmophere layer (positive down)
         ex_dsetl_dtr  ! derivative of setl_flux w.r.t. tracer concentration
#endif

    ! tracer seditation flux from the lowest atmospehere layer on the atm grid
    real :: setl_flux(size(Atm%tr_bot,1),size(Atm%tr_bot,2))

    ! derivative of setl_dtr from the lowest atmosphere layer on the atm grid
    real :: dsetl_dtr(size(Atm%tr_bot,1),size(Atm%tr_bot,2))
    

    ! temporary arrays
    real, dimension(n_xgrid_sfc) :: ex_gamma, ex_dtmass, ex_delta_t, ex_delta_u, ex_delta_v, ex_dflux_t

    ! generic exchange fields between atm and land
    real, dimension(n_xgrid_sfc,n_gex_atm2lnd) ::  ex_gex_atm2lnd

    real, dimension(n_xgrid_sfc,n_exch_tr) :: &
         ex_delta_tr, & ! tracer tendencies
         ex_dflux_tr    ! fracer flux change

    real    :: cp_inv ! inverse heat capacity at constant pressure
    logical :: used, ov ! used in calls to fms
    integer :: ier ! used in calls to fms
    integer :: is_atm, ie_atm, js_atm, je_atm, j ! do loop indices

    character(32) :: tr_name ! name of the tracer
    integer :: tr, n, m ! tracer indices
    integer :: is, ie, l, i, n_gex ! do loop indices


    !> START CLOCKS FOR PROFILING
    !{
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(fluxAtmDnClock)
    !}

    !> INITIALIZE REUSABLE FLAG.  DATA_OVERRIDE WILL RETURN OV=.TRUE. IF DATA WAS OVERWRITTEN
    !{
    ov = .FALSE.
    !}

    
    !> OVERRIDE ATM SHORTWAVE AND LONGWAVE DIRECT AND DOWNWARD DIFFUSIVE FLUXES
    !! NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE XARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    ! fms_data_override calls moved here from coupler_main
    call fms_data_override ('ATM', 'flux_sw',  Atm%flux_sw, Time)
    call fms_data_override ('ATM', 'flux_sw_dir',  Atm%flux_sw_dir, Time)
    call fms_data_override ('ATM', 'flux_sw_dif',  Atm%flux_sw_dif, Time)
    call fms_data_override ('ATM', 'flux_sw_down_vis_dir',  Atm%flux_sw_down_vis_dir, Time)
    call fms_data_override ('ATM', 'flux_sw_down_vis_dif',  Atm%flux_sw_down_vis_dif, Time)
    call fms_data_override ('ATM', 'flux_sw_down_total_dir',  Atm%flux_sw_down_total_dir, Time)
    call fms_data_override ('ATM', 'flux_sw_down_total_dif',  Atm%flux_sw_down_total_dif, Time)
    call fms_data_override ('ATM', 'flux_sw_vis',  Atm%flux_sw_vis, Time)
    call fms_data_override ('ATM', 'flux_sw_vis_dir',  Atm%flux_sw_vis_dir, Time)
    call fms_data_override ('ATM', 'flux_sw_vis_dif',  Atm%flux_sw_vis_dif, Time)
    call fms_data_override ('ATM', 'flux_lw',  Atm%flux_lw, Time)
    call fms_data_override ('ATM', 'lprec',    Atm%lprec,   Time)
    !}


    !> SCALE LIQUID PRECIPITATION BY FRAC_PRECIP IF SCALE_PRECIP_2D IS TRUE
    !! SCALE_PRECIP_2D IS SET DURING MODULE INITIALIZATION CALL TO ATM_LAND_ICE_FLUX_EXCHANGE_INIT
    !! FRAC_PRECIP VALUES ARE SET WITH DATA_OVERRIDE
    !{
    if (scale_precip_2d) then
       call fms_mpp_domains_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       call fms_data_override ('ATM', 'precip_scale2d',    frac_precip,   Time)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             Atm%lprec(i,j) = Atm%lprec(i,j)*frac_precip(i,j)
          enddo
       enddo
    endif
    !}

    
    !> PARTITION PRECIPTATION TO LIQUID PRECIPITATION AND FROZEN PRECIPITATION IF PARTITION_FPREC_FROM_LPREC = .TRUE.
    !! PARTIION_FPREC_FROM_LPREC IS SET AS PART OF MODULE INITIALIZATION CALL IN ATM_LAND_ICE_FLUX_EXCHANGE
    !{
    if (partition_fprec_from_lprec .and. Atm%pe) then
       call fms_mpp_domains_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             if (Atm%t_bot(i,j) < tfreeze) then
                Atm%fprec(i,j) = Atm%lprec(i,j)
                Atm%lprec(i,j) = 0.0
             endif
          enddo
       enddo
    endif
    !}


    !> OVERRIDE ATM FPREC, COZEN, AND SURF_DIFF FIELDS.
    !! NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE ARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    call fms_data_override ('ATM', 'fprec',    Atm%fprec,   Time)
    call fms_data_override ('ATM', 'coszen',   Atm%coszen,  Time)
    call fms_data_override ('ATM', 'dtmass',   Atm%Surf_Diff%dtmass, Time)
    call fms_data_override ('ATM', 'delta_t',  Atm%Surf_Diff%delta_t, Time)
    call fms_data_override ('ATM', 'dflux_t',  Atm%Surf_Diff%dflux_t, Time)
    do tr = 1,n_atm_tr
       call fms_tracer_manager_get_tracer_names(MODEL_ATMOS,tr,tr_name)
       call fms_data_override ('ATM', 'delta_'//trim(tr_name),  Atm%Surf_Diff%delta_tr(:,:,tr), Time)
       call fms_data_override ('ATM', 'dflux_'//trim(tr_name),  Atm%Surf_Diff%dflux_tr(:,:,tr), Time)
    enddo
    !}
    

    !> MAP ATMOSPHERE QUANTITIES ONTO THE EXCHANGE GRID
    !{
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_sw_dir, &
    !$OMP                                  ex_flux_sw_vis_dir,ex_flux_sw_dif,ex_delta_u,    &
    !$OMP                                  ex_flux_sw_vis_dif,ex_flux_lwd,ex_delta_v,       &
    !$OMP                                  ex_gex_atm2lnd,n_gex_atm2lnd)                    &
    !$OMP                          private(is,ie,n_gex)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_flux_sw_dir(i)     = 0.0
          ex_flux_sw_vis_dir(i) = 0.0
          ex_flux_sw_dif(i)     = 0.0
          ex_flux_sw_vis_dif(i) = 0.0
          ex_flux_lwd(i)        = 0.0
          ex_delta_u(i)         = 0.0
          ex_delta_v(i)         = 0.0
       enddo
       do n_gex=1,n_gex_atm2lnd
          do i = is, ie
             ex_gex_atm2lnd(i,n_gex)   = 0.0
          end do
       end do
    enddo
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_dir, 'ATM', ex_flux_sw_dir, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_vis_dir, 'ATM', ex_flux_sw_vis_dir, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_dif, 'ATM', ex_flux_sw_dif, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_vis_dif, 'ATM', ex_flux_sw_vis_dif, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_down_vis_dir, 'ATM', ex_flux_sw_down_vis_dir, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_down_total_dir,'ATM', ex_flux_sw_down_total_dir, xmap_sfc,complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_down_vis_dif,'ATM', ex_flux_sw_down_vis_dif, xmap_sfc,complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%flux_sw_down_total_dif, 'ATM',ex_flux_sw_down_total_dif, xmap_sfc,complete=.false.)

    !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
    !  if (fms_mpp_pe()== fms_mpp_root_pe()) print *,'LPREC', ccc
    
    !!$  if(do_area_weighted_flux) then
    !!$     call put_to_xgrid (Atm%lprec * AREA_ATM_MODEL,   'ATM', ex_lprec, xmap_sfc)
    !!$     call put_to_xgrid (Atm%fprec * AREA_ATM_MODEL,   'ATM', ex_fprec, xmap_sfc)
    !!$  else
    call fms_xgrid_put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc, complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%t_bot,   'ATM', ex_tprec, xmap_sfc, complete=.false.)
    !!$  endif

    do n_gex=1,n_gex_atm2lnd
       call fms_xgrid_put_to_xgrid (Atm%gex_atm2lnd(:,:,n_gex),'ATM',ex_gex_atm2lnd(:,n_gex),xmap_sfc,complete=.false.)
    end do

    call fms_xgrid_put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc, complete=.true.)
    call fms_xgrid_put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc, remap_method=remap_method, complete=.false.)
    !}

    
    !> ON THE EXCHANGE GRID, UPDATE U AND V STRESS
    !{
    ! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
    ! on exchange grid instead of the stresses themselves so that only the
    ! implicit corrections are filtered through the atmospheric grid not the
    ! stresses themselves
    call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%delta_u, 'ATM', ex_delta_u, xmap_sfc, remap_method=remap_method, &
                                 complete=.false.)
    call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%delta_v, 'ATM', ex_delta_v, xmap_sfc, remap_method=remap_method, &
                                 complete=.true.)

    ! MOD update stresses using atmos delta's but derivatives on exchange grid
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_u,ex_delta_u, &
    !$OMP                                  ex_dtaudu_atm,ex_dtaudv_atm,ex_flux_v,ex_delta_v )     &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_flux_u(i) = ex_flux_u(i) + ex_delta_u(i)*ex_dtaudu_atm(i)
          ex_flux_v(i) = ex_flux_v(i) + ex_delta_v(i)*ex_dtaudv_atm(i)
       enddo
    enddo
    !}
    

    !> ON THE EXCHANGE GRID, FIX SHORTWAVE RADIATION FLUX OF VISIBLE LIGHT TO TAKE INTO ACCOUNT FOR ALBEDO VARIATION
    !{
    !---- adjust sw flux for albedo variations on exch grid ----
    !---- adjust 4 categories (vis/nir dir/dif) separately  ----
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_sw_dir, &
    !$OMP                                  ex_flux_sw_vis_dir,ex_albedo_nir_dir_fix,        &
    !$OMP                                  ex_albedo_vis_dir_fix,ex_flux_sw_dif,            &
    !$OMP                                  ex_flux_sw_vis_dif,ex_flux_sw_vis,ex_flux_sw,    &
    !$OMP                                  ex_albedo_nir_dif_fix,ex_albedo_vis_dif_fix   )  &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) - ex_flux_sw_vis_dir(i)     ! temporarily nir/dir
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) * ex_albedo_nir_dir_fix(i)  ! fix nir/dir
          ex_flux_sw_vis_dir(i) = ex_flux_sw_vis_dir(i) * ex_albedo_vis_dir_fix(i) ! fix vis/dir
          ex_flux_sw_dir(i) = ex_flux_sw_dir(i) + ex_flux_sw_vis_dir(i)     ! back to total dir

          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) - ex_flux_sw_vis_dif(i)     ! temporarily nir/dif
          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) * ex_albedo_nir_dif_fix(i)  ! fix nir/dif
          ex_flux_sw_vis_dif(i) = ex_flux_sw_vis_dif(i) * ex_albedo_vis_dif_fix(i) ! fix vis/dif
          ex_flux_sw_dif(i) = ex_flux_sw_dif(i) + ex_flux_sw_vis_dif(i)     ! back to total dif

          ex_flux_sw_vis(i) = ex_flux_sw_vis_dir(i) + ex_flux_sw_vis_dif(i) ! legacy, remove later
          ex_flux_sw(i)     = ex_flux_sw_dir(i)     + ex_flux_sw_dif(i)     ! legacy, remove later
       enddo
    enddo

    !!$  ex_flux_sw_dir = ex_flux_sw_dir - ex_flux_sw_vis_dir            ! temporarily nir/dir
    !!$  ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_nir_dir_fix         ! fix nir/dir
    !!$  ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix ! fix vis/dir
    !!$  ex_flux_sw_dir = ex_flux_sw_dir + ex_flux_sw_vis_dir            ! back to total dir
    !!$
    !!$  ex_flux_sw_dif = ex_flux_sw_dif - ex_flux_sw_vis_dif            ! temporarily nir/dif
    !!$  ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_nir_dif_fix         ! fix nir/dif
    !!$  ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix ! fix vis/dif
    !!$  ex_flux_sw_dif = ex_flux_sw_dif + ex_flux_sw_vis_dif            ! back to total dif
    !!$
    !!$  ex_flux_sw_vis = ex_flux_sw_vis_dir + ex_flux_sw_vis_dif        ! legacy, remove later
    !!$  ex_flux_sw     = ex_flux_sw_dir     + ex_flux_sw_dif            ! legacy, remove later

    deallocate ( ex_albedo_vis_dir_fix )
    deallocate ( ex_albedo_nir_dir_fix )
    deallocate ( ex_albedo_vis_dif_fix )
    deallocate ( ex_albedo_nir_dif_fix )
    !}
    

    !> ON THE EXCHANGE GRID, ADJUST FLUXES FOR IMPLICIT DEPENDENCE
    !{
    do tr = 1,n_exch_tr
       n = tr_table(tr)%atm
       call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%delta_tr(:,:,n), 'ATM', ex_delta_tr(:,tr), xmap_sfc, complete=.false.)
       call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%dflux_tr(:,:,n), 'ATM', ex_dflux_tr(:,tr), xmap_sfc, complete=.false.)
    enddo

    call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc, complete=.false. )
    call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc, complete=.false. )
    call fms_xgrid_put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc, complete=.true. )

#ifndef use_AM3_physics
    ! Get sedimentation flux. Has to be here (instead of sfc_boundary_layer sub)
    ! because of time stepping order: sedimentation fluxes are calculated in
    ! update_atmos_model_down (in atmos_tracer_driver), but sfc_boundary_layer
    ! is called before that.
    do tr = 1,n_exch_tr
       if (atmos_tracer_has_surf_setl_flux(tr_table(tr)%atm)) then
          call get_atmos_tracer_surf_setl_flux (tr_table(tr)%atm, setl_flux, dsetl_dtr)
          call fms_xgrid_put_to_xgrid(setl_flux, 'ATM', ex_setl_flux, xmap_sfc)
          call fms_xgrid_put_to_xgrid(dsetl_dtr, 'ATM', ex_dsetl_dtr, xmap_sfc)
          where (ex_avail)
             ! minus sign is because sedimentation is positive down
             ex_flux_tr(:,tr)   = ex_flux_tr(:,tr)   - ex_setl_flux(:)
             ex_dfdtr_atm(:,tr) = ex_dfdtr_atm(:,tr) - ex_dsetl_dtr(:)
          end where
       endif
    enddo
#endif

    cp_inv = 1.0/cp_air

    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_lw,ex_flux_lwd, &
    !$OMP                                  ex_avail,ex_gamma,ex_dtmass,ex_dflux_t,ex_dhdt_atm,      &
    !$OMP                                  cp_inv,ex_e_t_n,ex_dhdt_surf,ex_f_t_delt_n,ex_delta_t,   &
    !$OMP                                  ex_flux_t,ex_dflux_tr,isphum,ex_dfdtr_atm,ex_e_q_n,      &
    !$OMP                                  ex_dedt_surf,n_exch_tr,ex_e_tr_n,ex_dfdtr_surf,          &
    !$OMP                                  ex_f_tr_delt_n,ex_delta_tr,ex_flux_tr )                  &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          !----- compute net longwave flux (down-up) -----
          ! (note: lw up already in ex_flux_lw)
          ex_flux_lw(i) = ex_flux_lwd(i) - ex_flux_lw(i)
          if (ex_avail(i) ) then

             ! temperature
             ex_gamma(i)      =  1./ (1.0 - ex_dtmass(i)*(ex_dflux_t(i) + ex_dhdt_atm(i)*cp_inv))
             ex_e_t_n(i)      =  ex_dtmass(i)*ex_dhdt_surf(i)*cp_inv*ex_gamma(i)
             ex_f_t_delt_n(i) = (ex_delta_t(i) + ex_dtmass(i) * ex_flux_t(i)*cp_inv) * ex_gamma(i)

             ex_flux_t (i)    =  ex_flux_t(i)        + ex_dhdt_atm(i) * ex_f_t_delt_n(i)
             ex_dhdt_surf(i)  =  ex_dhdt_surf(i)     + ex_dhdt_atm(i) * ex_e_t_n(i)

             ! moisture
             !     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))
             ! here it looks like two derivatives with different units are added together,
             ! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
             ! regions of exchange grid, so that if one of them is not zero the other is, and
             ! vice versa.
             !     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma
             !     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma
             !     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n
             !     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
             !     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
             ! moisture vs. surface temperture, assuming saturation
             ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,isphum) + ex_dfdtr_atm(i,isphum)))
             ex_e_q_n(i)      =  ex_dtmass(i) * ex_dedt_surf(i) * ex_gamma(i)
             ex_dedt_surf(i)  =  ex_dedt_surf(i) + ex_dfdtr_atm(i,isphum) * ex_e_q_n(i)
             do tr = 1,n_exch_tr
                ex_gamma(i)   =  1.0 / (1.0 - ex_dtmass(i)*(ex_dflux_tr(i,tr) + ex_dfdtr_atm(i,tr)))

                ex_e_tr_n(i,tr)      =  ex_dtmass(i)*ex_dfdtr_surf(i,tr)*ex_gamma(i)
                ex_f_tr_delt_n(i,tr) = (ex_delta_tr(i,tr)+ex_dtmass(i)*ex_flux_tr(i,tr))*ex_gamma(i)

                ex_flux_tr(i,tr)     =  ex_flux_tr(i,tr) + ex_dfdtr_atm(i,tr)*ex_f_tr_delt_n(i,tr)
                ex_dfdtr_surf(i,tr)  =  ex_dfdtr_surf(i,tr) + ex_dfdtr_atm(i,tr)*ex_e_tr_n(i,tr)
             enddo
          endif
       enddo
    enddo 
    !}
    

    !> MAP FLUXES FROM THE EXCHANGE GRID TO THE LAND GRID AND OVERRIDE FIELDS
    !! WITH DATA_OVERRIDE WHERE DATA WILL BE OVERWRITTEN IF THE FIELD IS SPECIFIED IN DATA_TABLE
    !{
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,  xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,  xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
#ifdef SCM
    if (do_specified_land .and. do_specified_flux) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%dhdt,  'LND', ex_dhdt_surf_forland, xmap_sfc)
    else
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%dhdt,  'LND', ex_dhdt_surf, xmap_sfc)
    endif
#else
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
#endif
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

    call fms_xgrid_get_from_xgrid_ug (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
    call fms_xgrid_get_from_xgrid_ug (Land_boundary%tprec,   'LND', ex_tprec,     xmap_sfc)
    !!$  if(do_area_weighted_flux) then
    !!$     ! evap goes here???
    !!$     do k = 1, size(Land_boundary%lprec, dim=3)
    !!$        ! Note: we divide by AREA_ATM_MODEL, which should be the same as
    !!$        ! AREA_LND_MODEL (but the latter may not be defined)
    !!$        call divide_by_area(data=Land_boundary%lprec(:,:,k), area=AREA_ATM_MODEL)
    !!$        call divide_by_area(data=Land_boundary%fprec(:,:,k), area=AREA_ATM_MODEL)
    !!$     enddo
    !!$  endif

    if(associated(Land_boundary%drag_q)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
       call fms_data_override_ug('LND', 'drag_q', Land_boundary%drag_q,  Time )
    endif
    if(associated(Land_boundary%lwdn_flux)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
       call fms_data_override_ug('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
    endif
    if(associated(Land_boundary%cd_m)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
       call fms_data_override_ug('LND', 'cd_m', Land_boundary%cd_m, Time )
    endif
    if(associated(Land_boundary%cd_t)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
       call fms_data_override_ug('LND', 'cd_t', Land_boundary%cd_t, Time )
    endif
    if(associated(Land_boundary%bstar)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
       call fms_data_override_ug('LND', 'bstar',  Land_boundary%bstar, Time )
    endif
    if(associated(Land_boundary%ustar)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
       call fms_data_override_ug('LND', 'ustar',  Land_boundary%ustar, Time )
    endif
    if(associated(Land_boundary%wind)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
       call fms_data_override_ug('LND', 'wind',  Land_boundary%wind, Time )
    endif
    if(associated(Land_boundary%z_bot)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
       call fms_data_override_ug('LND', 'z_bot',  Land_boundary%z_bot, Time )
    endif
    if (associated(Land_boundary%con_atm)) then
       call fms_xgrid_get_from_xgrid_ug (Land_boundary%con_atm, 'LND', ex_con_atm, xmap_sfc)
    end if

    if (associated(Land_boundary%gex_atm2lnd)) then
       do n_gex=1,n_gex_atm2lnd
          call fms_xgrid_get_from_xgrid_ug (Land_boundary%gex_atm2lnd(:,:,n_gex),'LND',ex_gex_atm2lnd(:,n_gex),xmap_sfc)
          !add data_override here
       end do
    end if
#else
    call fms_xgrid_get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
#ifdef SCM
    if (do_specified_land .and. do_specified_flux) then
       call fms_xgrid_get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf_forland, xmap_sfc)
    else
       call fms_xgrid_get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf, xmap_sfc)
    endif
#else
    call fms_xgrid_get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
#endif
    call fms_xgrid_get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

    call fms_xgrid_get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
    call fms_xgrid_get_from_xgrid (Land_boundary%tprec,   'LND', ex_tprec,     xmap_sfc)

    if (associated(Land_boundary%gex_atm2lnd)) then
       do n_gex=1,n_gex_atm2lnd
          call fms_xgrid_get_from_xgrid (Land_boundary%gex_atm2lnd(:,:,n_gex), 'LND', ex_gex_atm2lnd(:,n_gex), xmap_sfc)
       end do
    end if

    !!$  if(do_area_weighted_flux) then
    !!$     ! evap goes here???
    !!$     do k = 1, size(Land_boundary%lprec, dim=3)
    !!$        ! Note: we divide by AREA_ATM_MODEL, which should be the same as
    !!$        ! AREA_LND_MODEL (but the latter may not be defined)
    !!$        call divide_by_area(data=Land_boundary%lprec(:,:,k), area=AREA_ATM_MODEL)
    !!$        call divide_by_area(data=Land_boundary%fprec(:,:,k), area=AREA_ATM_MODEL)
    !!$     enddo
    !!$  endif

    if(associated(Land_boundary%drag_q)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
       call fms_data_override('LND', 'drag_q', Land_boundary%drag_q,  Time )
    endif
    if(associated(Land_boundary%lwdn_flux)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
       call fms_data_override('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
    endif
    if(associated(Land_boundary%cd_m)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
       call fms_data_override('LND', 'cd_m', Land_boundary%cd_m, Time )
    endif
    if(associated(Land_boundary%cd_t)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
       call fms_data_override('LND', 'cd_t', Land_boundary%cd_t, Time )
    endif
    if(associated(Land_boundary%bstar)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
       call fms_data_override('LND', 'bstar',  Land_boundary%bstar, Time )
    endif
    if(associated(Land_boundary%ustar)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
       call fms_data_override('LND', 'ustar',  Land_boundary%ustar, Time )
    endif
    if(associated(Land_boundary%wind)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
       call fms_data_override('LND', 'wind',  Land_boundary%wind, Time )
    endif
    if(associated(Land_boundary%z_bot)) then
       call fms_xgrid_get_from_xgrid (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
       call fms_data_override('LND', 'z_bot',  Land_boundary%z_bot, Time )
    endif
#endif
    Land_boundary%tr_flux = 0.0
    Land_boundary%dfdtr = 0.0
    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then

#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_get_from_xgrid_ug (Land_boundary%tr_flux(:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
          call fms_xgrid_get_from_xgrid_ug (Land_boundary%dfdtr(:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
#else
          call fms_xgrid_get_from_xgrid (Land_boundary%tr_flux(:,:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
          call fms_xgrid_get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
#endif
#ifdef SCM
          if (do_specified_land .and. do_specified_flux .and. tr.eq.isphum) then
#ifndef _USE_LEGACY_LAND_
             call fms_xgrid_get_from_xgrid_ug (Land_boundary%dfdtr(:,:,n),   'LND', ex_dedq_surf_forland(:), xmap_sfc)
#else
             call fms_xgrid_get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dedq_surf_forland(:), xmap_sfc)
#endif
          endif
#endif
       endif
    enddo
    !}

    
    !> OVERRIDE LAND FLUXES.  NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE
    !! ARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    !  current time is Time: is that ok? not available in land_data_type
    !Balaji: data_override calls moved here from coupler_main
#ifndef _USE_LEGACY_LAND_
    call fms_data_override_ug('LND', 't_flux',  Land_boundary%t_flux,  Time )
    call fms_data_override_ug('LND', 'lw_flux', Land_boundary%lw_flux, Time )
    call fms_data_override_ug('LND', 'sw_flux', Land_boundary%sw_flux, Time )
    call fms_data_override_ug('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
    call fms_data_override_ug('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
    call fms_data_override_ug('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
    call fms_data_override_ug('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )

    call fms_data_override_ug('LND', 'lprec',   Land_boundary%lprec,   Time )
    call fms_data_override_ug('LND', 'fprec',   Land_boundary%fprec,   Time )
    call fms_data_override_ug('LND', 'dhdt',    Land_boundary%dhdt,    Time )
    call fms_data_override_ug('LND', 'drdt',    Land_boundary%drdt,    Time )
    call fms_data_override_ug('LND', 'p_surf',  Land_boundary%p_surf,  Time )
    do tr = 1,n_lnd_tr
       call fms_tracer_manager_get_tracer_names(MODEL_LAND, tr, tr_name)
       call fms_data_override_ug('LND', trim(tr_name)//'_flux', Land_boundary%tr_flux(:,:,tr), Time)
       call fms_data_override_ug('LND', 'dfd'//trim(tr_name),   Land_boundary%dfdtr  (:,:,tr), Time)
#else
    call fms_data_override('LND', 't_flux',  Land_boundary%t_flux,  Time )
    call fms_data_override('LND', 'lw_flux', Land_boundary%lw_flux, Time )
    call fms_data_override('LND', 'sw_flux', Land_boundary%sw_flux, Time )
    call fms_data_override('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
    call fms_data_override('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
    call fms_data_override('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
    call fms_data_override('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )

    call fms_data_override('LND', 'lprec',   Land_boundary%lprec,   Time )
    call fms_data_override('LND', 'fprec',   Land_boundary%fprec,   Time )
    call fms_data_override('LND', 'dhdt',    Land_boundary%dhdt,    Time )
    call fms_data_override('LND', 'drdt',    Land_boundary%drdt,    Time )
    call fms_data_override('LND', 'p_surf',  Land_boundary%p_surf,  Time )
    do tr = 1,n_lnd_tr
       call fms_tracer_manager_get_tracer_names(MODEL_LAND, tr, tr_name)
       call fms_data_override('LND', trim(tr_name)//'_flux', Land_boundary%tr_flux(:,:,:,tr), Time)
       call fms_data_override('LND', 'dfd'//trim(tr_name),   Land_boundary%dfdtr  (:,:,:,tr), Time)
#endif
    enddo
    !}
    

    !> MAP ICE FIELDS FROM THE EXCHANGE GRID TO THE ICE GRID
    !{
    call fms_xgrid_get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_tr(:,isphum), xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_flux_vis_dir,  'OCN', ex_flux_sw_vis_dir,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_flux_nir_dir,  'OCN', ex_flux_sw_dir,xmap_sfc)
    ! ice & ocean use these 4: dir/dif nir/vis
    Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_flux_nir_dir - Ice_boundary%sw_flux_vis_dir
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_flux_vis_dif,  'OCN', ex_flux_sw_vis_dif,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_flux_nir_dif,  'OCN', ex_flux_sw_dif,xmap_sfc)
    ! ice & ocean use these 4: dir/dif nir/vis
    Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_flux_nir_dif - Ice_boundary%sw_flux_vis_dif

    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_down_vis_dir,  'OCN', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_down_nir_dir,  'OCN', ex_flux_sw_down_total_dir, xmap_sfc)
    ! ice & ocean use these 4: dir/dif nir/vis
    Ice_boundary%sw_down_nir_dir = Ice_boundary%sw_down_nir_dir - Ice_boundary%sw_down_vis_dir

    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_down_vis_dif,  'OCN', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%sw_down_nir_dif,  'OCN', ex_flux_sw_down_total_dif,xmap_sfc)
    ! ice & ocean use these 4: dir/dif nir/vis
    Ice_boundary%sw_down_nir_dif = Ice_boundary%sw_down_nir_dif - Ice_boundary%sw_down_vis_dif

    call fms_xgrid_get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%u_star,   'OCN', ex_u_star,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%p,        'OCN', ex_slp,       xmap_sfc) ! mw mod

    call fms_xgrid_get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
    call fms_xgrid_get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
    !!$  if (do_area_weighted_flux) then
    !!$     where (AREA_ATM_SPHERE /= 0)
    !!$        Ice_boundary%lprec = Ice_boundary%lprec * AREA_ATM_MODEL/AREA_ATM_SPHERE
    !!$        Ice_boundary%fprec = Ice_boundary%fprec * AREA_ATM_MODEL/AREA_ATM_SPHERE
    !!$     end where
    !!$  endif
    !!$  if(do_area_weighted_flux) then
    !!$     do k = 1, size(Ice_boundary%lprec, dim=3)
    !!$        call divide_by_area(data=Ice_boundary%lprec(:,:,k), area=AREA_ATM_SPHERE)
    !!$        call divide_by_area(data=Ice_boundary%fprec(:,:,k), area=AREA_ATM_SPHERE)
    !!$     enddo
    !!$  endif

    ! Extra fluxes
    do n = 1, Ice_boundary%fluxes%num_bcs
      if(ex_gas_fluxes%bc(n)%flux_type  .ne. 'air_sea_deposition') then
       do m = 1, Ice_boundary%fluxes%bc(n)%num_fields
          call fms_xgrid_get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
               ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)
       enddo
      endif
    enddo
    !}


    !> OVERRIDE ICE FIELDS.  NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE
    !! ARRAY IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    ! data_override calls moved here from coupler_main
    call fms_data_override('ICE', 'u_flux', Ice_boundary%u_flux,  Time)
    call fms_data_override('ICE', 'v_flux', Ice_boundary%v_flux,  Time)
    call fms_data_override('ICE', 't_flux', Ice_boundary%t_flux,  Time)
    call fms_data_override('ICE', 'q_flux', Ice_boundary%q_flux,  Time)
    call fms_data_override('ICE', 'lw_flux',Ice_boundary%lw_flux, Time)
    call fms_data_override('ICE', 'lw_flux_dn',Ice_boundary%lw_flux, Time, override=ov)
    if (ov) then
       Ice_boundary%lw_flux = Ice_boundary%lw_flux - stefan*Ice%t_surf**4
    endif
    call fms_data_override('ICE', 'sw_flux_nir_dir',Ice_boundary%sw_flux_nir_dir, Time)
    call fms_data_override('ICE', 'sw_flux_vis_dir',Ice_boundary%sw_flux_vis_dir, Time)
    call fms_data_override('ICE', 'sw_flux_nir_dif',Ice_boundary%sw_flux_nir_dif, Time, override=ov)
    call fms_data_override('ICE', 'sw_flux_vis_dif',Ice_boundary%sw_flux_vis_dif, Time)
    call fms_data_override('ICE', 'sw_flux_vis_dir_dn',Ice_boundary%sw_down_vis_dir, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_vis_dir = Ice_boundary%sw_down_vis_dir*(1.0-Ice%albedo_vis_dir)
    endif
    call fms_data_override('ICE', 'sw_flux_vis_dif_dn',Ice_boundary%sw_down_vis_dif, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_vis_dif = Ice_boundary%sw_down_vis_dif*(1.0-Ice%albedo_vis_dif)
    endif
    call fms_data_override('ICE', 'sw_flux_nir_dir_dn',Ice_boundary%sw_down_nir_dir, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_down_nir_dir*(1.0-Ice%albedo_nir_dir)
    endif
    call fms_data_override('ICE', 'sw_flux_nir_dif_dn',Ice_boundary%sw_down_nir_dif, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_down_nir_dif*(1.0-Ice%albedo_nir_dif)
    endif
    call fms_data_override('ICE', 'lprec',  Ice_boundary%lprec,   Time)
    call fms_data_override('ICE', 'fprec',  Ice_boundary%fprec,   Time)
    call fms_data_override('ICE', 'dhdt',   Ice_boundary%dhdt,    Time)
    call fms_data_override('ICE', 'dedt',   Ice_boundary%dedt,    Time)
    call fms_data_override('ICE', 'drdt',   Ice_boundary%drdt,    Time)
    call fms_data_override('ICE', 'coszen', Ice_boundary%coszen,  Time)
    call fms_data_override('ICE', 'p',      Ice_boundary%p,       Time)

    call fms_coupler_type_data_override('ICE', Ice_boundary%fluxes, Time)
    call fms_coupler_type_send_data(Ice_boundary%fluxes, Time)
    !}
    

    !> COMPUTE STOCK CHANGES BETWEEN COMPONENTS
    !{
    ! Atm -> Lnd (precip)
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_stock_move_ug( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_WATER),  &
         & TO   = fms_stock_constants_lnd_stock(ISTOCK_WATER), &
         & stock_ug_data3d = (Land_boundary%lprec + Land_boundary%fprec), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Lnd) ')

    ! Atm -> Lnd (heat)
    call fms_xgrid_stock_move_ug( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_HEAT),  &
         & TO   = fms_stock_constants_lnd_stock(ISTOCK_HEAT), &
         & stock_ug_data3d = (-Land_boundary%t_flux + Land_boundary%lw_flux +  Land_boundary%sw_flux - &
                              Land_boundary%fprec*HLF), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Lnd) ')
#else
    call fms_xgrid_stock_move( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_WATER),  &
         & TO   = fms_stock_constants_lnd_stock(ISTOCK_WATER), &
         & stock_data3d = (Land_boundary%lprec + Land_boundary%fprec), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Lnd) ')

    ! Atm -> Lnd (heat)
    call fms_xgrid_stock_move( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_HEAT),  &
         & TO   = fms_stock_constants_lnd_stock(ISTOCK_HEAT), &
         & stock_data3d = (-Land_boundary%t_flux + Land_boundary%lw_flux +  Land_boundary%sw_flux - &
                            Land_boundary%fprec*HLF), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Lnd) ')
#endif

    ! Atm -> Ice (precip)
    call fms_xgrid_stock_move( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_WATER), &
         & TO   = fms_stock_constants_ice_stock(ISTOCK_WATER), &
         & stock_data3d = (Ice_boundary%lprec + Ice_boundary%fprec), &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Ice) ')

    ! Atm -> Ice (heat)
    call fms_xgrid_stock_move( &
         & FROM = fms_stock_constants_atm_stock(ISTOCK_HEAT), &
         & TO   = fms_stock_constants_ice_stock(ISTOCK_HEAT), &
         & stock_data3d = (-Ice_boundary%t_flux + Ice_boundary%lw_flux - Ice_boundary%fprec*HLF + &
                            Ice_boundary%sw_flux_vis_dir + &
         Ice_boundary%sw_flux_vis_dif + Ice_boundary%sw_flux_nir_dir + Ice_boundary%sw_flux_nir_dif), &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Ice) ')
    !}

    
    deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)


    !> SEND U_FLUX AND V_FLUX TO THE DIAG_MANAGER BUFFER
    !! NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FOUND IN DIAG_TABLE.YAML
    !{
    !-------------------- diagnostics section ------------------------------

    !------- zonal wind stress -----------
    used = fms_diag_send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
    used = fms_diag_send_data ( id_tauu,  -Atmos_boundary%u_flux, Time )

    !------- meridional wind stress -----------
    used = fms_diag_send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
    used = fms_diag_send_data ( id_tauv,  -Atmos_boundary%v_flux, Time )
    !}


    !> END CLOCK FOR PROFILING
    !{
    call fms_mpp_clock_end(fluxAtmDnClock)
    call fms_mpp_clock_end(cplClock)
    !}

  end subroutine flux_down_from_atmos

  !> Subroutine generate_sfc_xgrid updates the fractional area of the land-ice exchange grid, where
  !! the fractional area measures the portion of the exchange grid cell that correspoonds to land and to ice.
  !! This subroutine reduces the size of the exchange grid by eliminating exchange grid cells that are
  !! pure land or pure ice (i.e., eliminate side 2 tiles with fractional area value of 0.0)
  subroutine generate_sfc_xgrid( Land, Ice )

    !> derived data type to specify land boundary data
    type(land_data_type), intent(in) :: Land 
    !> derived data type to specify ice boundary dat
    type(ice_data_type),  intent(in) :: Ice 

    ! compute domain indices
    integer :: isc, iec, jsc, jec 


    !> INITIALIZE CLOCK FOR PROFILING
    !{
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(regenClock)
    !}
    

    !> GET ICE COMPUTE DOMAIN INDICES
    !{
    call fms_mpp_domains_get_compute_domain(Ice%Domain, isc, iec, jsc, jec)
    !}
    

    !> UPDATE FRACTIONAL AREAS OF THE EXCHANGE GRID THAT ARE ICE AND LAND
    !{
    call fms_xgrid_set_frac_area (Ice%part_size(isc:iec,jsc:jec,:) , 'OCN', xmap_sfc)
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_set_frac_area_ug (Land%tile_size, 'LND', xmap_sfc)
#else
    call fms_xgrid_set_frac_area (Land%tile_size, 'LND', xmap_sfc)
#endif
    !}


    !> UPDATE THE NUMBER OF EXCHANGE GRID CELLS SAVED IN THE MODULE
    !{
    n_xgrid_sfc = max(fms_xgrid_count(xmap_sfc),1)
    if(n_xgrid_sfc .GE. nblocks) then
       my_nblocks = nblocks
       call fms_mpp_domains_compute_extent(1, n_xgrid_sfc, nblocks, block_start, block_end)
    else
       my_nblocks = 1
       block_start(1) = 1
       block_end(1) = n_xgrid_sfc
    endif
    !}
    

    !> END CLOCK FOR PROFILING
    !{
    call fms_mpp_clock_end(regenClock)
    call fms_mpp_clock_end(cplClock)
    !}
    
  end subroutine generate_sfc_xgrid

  !> Subroutine flux_up_to_atmos corrects the fluxes to take into account
  !! the new surface temperatures in land and ice models.
  !! Final increments for temperature and specific humidity in the
  !! lowest atmospheric layer are computed and returned to the atmospheric model
  !! in order to finalize the increments in the rest of the atmosphere.
  !!
  !! The following elements of the land_ice_atmos_boundary_type are computed:
  !!   dt_t = temperature change at the lowest atmospheric level [K]
  !!   dt_q = specific humidity change at the lowest atmospheric level [kg/kg]
  subroutine flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary, Land_boundary, Ice_boundary )
    !> current model time
    type(FmsTime_type), intent(in)  :: Time !< Current time
    !> derived data type holding land boundary data
    type(land_data_type), intent(inout) :: Land
    !> derived data type holding ice boundary data
    type(ice_data_type),  intent(inout) :: Ice 
    !> derived data type holding properties and fluxes passed from exchange grid to the atmosphere, land and ice
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary 
    !> derived data type holding properties and fluxes passed from atmosphere to land
    type(atmos_land_boundary_type), intent(inout) :: Land_boundary
    !> derived data type holding properties and fluxes passed from atmosphere to ice
    type(atmos_ice_boundary_type),  intent(inout) :: Ice_boundary

    ! arrays on exchange grid
    real, dimension(n_xgrid_sfc) ::  &
         ex_t_surf_new, &
         ex_dt_t_surf,  &
         ex_delta_t_n,  &
         ex_t_ca_new,   &
         ex_dt_t_ca,    &
         ex_icetemp,    &
         ex_land_frac,  &
         ex_temp

    real, dimension(n_xgrid_sfc,n_exch_tr) :: &
         ! updated tracer values at the surface on exchange grid
         ex_tr_surf_new, & 
         ! tendency of tracers at the surface on exchange grid
         ex_dt_tr_surf, & 
         ! array on exchange grid
         ex_delta_tr_n

    ! added for co2_surf diagnostic, where co2_surf_dvmr is the updated CO2 tracer values at the surface (dry vmr)
    real, dimension(n_xgrid_sfc) :: ex_co2_surf_dvmr   

    real, dimension(size(Land_Ice_Atmos_Boundary%dt_t,1),size(Land_Ice_Atmos_Boundary%dt_t,2)) :: &
         diag_atm, &
         evap_atm, &
         frac_atm

#ifndef _USE_LEGACY_LAND_
    real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2)) :: data_lnd, diag_land
#else
    real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2), size(Land_boundary%lprec,3)) :: data_lnd,&
                                                                                                              diag_land
#endif
    real, dimension(size(Ice_boundary%lprec,1), size(Ice_boundary%lprec,2), size(Ice_boundary%lprec,3)) :: data_ice
    real, dimension(size(Ice%albedo,1),size(Ice%albedo,2),size(Ice%albedo,3)) ::  icegrid

    logical :: used
    integer :: tr ! tracer index
    character(32) :: tr_name, tr_units ! tracer name
    integer :: n, i, m, ier

    integer :: is, ie, l


    !> START CLOCK FOR PROFILING
    !{
    call fms_mpp_clock_begin(cplClock)
    call fms_mpp_clock_begin(fluxAtmUpClock)
    !}

    
    !> OVERRIDE ICE%T_SURF, LAND%T_CA, LAND%T_SURF AND LAND SURFACE TRACERS
    !! NOTE, DATA_OVERRIDE WILL ONLY OVERWRITE DATA IF THE FIELD IS SPECIFIED IN THE DATA_TABLE
    !{
    ! data_override calls moved here from coupler_main
    call fms_data_override ( 'ICE', 't_surf', Ice%t_surf,  Time)
#ifndef _USE_LEGACY_LAND_
    call fms_data_override_ug ( 'LND', 't_ca',   Land%t_ca,   Time)
    call fms_data_override_ug ( 'LND', 't_surf', Land%t_surf, Time)
#else
    call fms_data_override ( 'LND', 't_ca',   Land%t_ca,   Time)
    call fms_data_override ( 'LND', 't_surf', Land%t_surf, Time)
#endif
    do tr = 1, n_lnd_tr
       call fms_tracer_manager_get_tracer_names( MODEL_LAND, tr, tr_name )
#ifndef _USE_LEGACY_LAND_
       call fms_data_override_ug('LND', trim(tr_name)//'_surf', Land%tr(:,:,tr), Time)
#else
       call fms_data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
#endif
    enddo
    !}

    
    !> INITIALIZE EX_T_SURF_NEW = 200.0
    !{
    ex_t_surf_new = 200.0
    !}


    !> MAP ICE%T_SURF, LAND%T_CA AND LAND%T_SURF ONTO THE EXCHANGE GRID
    !{
    call fms_xgrid_put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
    ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_put_to_xgrid_ug (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
    call fms_xgrid_put_to_xgrid_ug (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)
#else
    call fms_xgrid_put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
    call fms_xgrid_put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)
#endif
    !}

    
    !  call escomp(ex_t_ca_new, ex_q_surf_new)
    !  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new)
    !  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

    
#ifdef SCM
    if (do_specified_flux .and. do_specified_land) then
       ex_t_surf_new = ex_t_surf
       ex_t_ca_new   = ex_t_ca
    endif
#endif

    
    !> ON THE EXCHANGE GRID, COMPUTE CHANGES IN SURFACE TEMPERATURE AND RADIATIVE TEMPERATURE
    !{
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          if(ex_avail(i)) then
             ex_dt_t_ca(i)  = ex_t_ca_new(i)   - ex_t_ca(i)   ! changes in near-surface T
             ex_dt_t_surf(i) = ex_t_surf_new(i) - ex_t_surf(i) ! changes in radiative T
          endif
       enddo

       if (do_forecast) then
          do i = is, ie
             if(ex_avail(i) .and. (.not.ex_land(i))) then
                ex_dt_t_ca  (i) = 0.
                ex_dt_t_surf(i) = 0.
             endif
          enddo
       end if
       !}
       
       
       !> ON THE EXCHANGE GRID,
       !! UPDATE FLUXES AND ATMOSPHERIC INCREMENTS FOR IMPLICIT DEPENDENCE ON SURFACE TEMPERATURE
       !{
       do tr = 1,n_exch_tr
          ! set up updated surface tracer field so that flux to atmos for absent
          ! tracers is zero
          do i = is,ie
             if(.not.ex_avail(i)) cycle
             if (ex_dfdtr_surf(i,tr)/=0.0) then
                ex_dt_tr_surf(i,tr) = -ex_flux_tr(i,tr)/ex_dfdtr_surf(i,tr)
             else
                ex_dt_tr_surf(i,tr) = 0.0
             endif
             ex_tr_surf_new(i,tr) = ex_tr_surf(i,tr)+ex_dt_tr_surf(i,tr)
          enddo
       enddo
    enddo
    ! get all tracers available from land, and calculate changes in near-tracer field
    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
#ifndef _USE_LEGACY_LAND_
          call fms_xgrid_put_to_xgrid_ug ( Land%tr(:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
#else
          call fms_xgrid_put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
#endif
       endif
    enddo
    !}
    

    ! get all tracers available from ocean here

    
    !> ON THE EXCHANGE GRID, UPDATE TRACER TENDENCIES IN THE ATMOSPHERE
    !{
    ! update tracer tendencies in the atmosphere
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do tr = 1,n_exch_tr
          do i = is, ie
             if(ex_avail(i)) then
                ex_dt_tr_surf(i,tr) = ex_tr_surf_new(i,tr) - ex_tr_surf(i,tr)
                ex_delta_tr_n(i,tr) = ex_f_tr_delt_n(i,tr) + ex_dt_tr_surf(i,tr) * ex_e_tr_n(i,tr)
                ex_flux_tr(i,tr)    = ex_flux_tr(i,tr)     + ex_dt_tr_surf(i,tr) * ex_dfdtr_surf(i,tr)
             endif
          enddo
       enddo

       ! re-calculate fluxes of specific humidity over ocean
       do i = is, ie
          if(ex_avail(i) .and. (.not.ex_land(i))) then
             ! note that in this region (over ocean) ex_dt_t_surf == ex_dt_t_ca
             ex_delta_tr_n(i,isphum)  = ex_f_tr_delt_n(i,isphum) + ex_dt_t_surf(i) * ex_e_q_n(i)
             ex_flux_tr(i,isphum)     = ex_flux_tr(i,isphum)     + ex_dt_t_surf(i) * ex_dedt_surf(i)
          endif
       enddo
    enddo    

    do tr=1,n_exch_tr
       ! get updated tracer tendency on the atmospheic grid
       n=tr_table(tr)%atm
       call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%dt_tr(:,:,n), 'ATM', ex_delta_tr_n(:,tr), xmap_sfc)
    enddo

    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_delta_t_n(i) = 0.0
          if(ex_avail(i)) then
             ex_flux_t(i)    = ex_flux_t(i)  + ex_dt_t_ca(i)   * ex_dhdt_surf(i)
             ex_flux_lw(i)   = ex_flux_lw(i) - ex_dt_t_surf(i) * ex_drdt_surf(i)
             ex_delta_t_n(i) = ex_f_t_delt_n(i)  + ex_dt_t_ca(i)*ex_e_t_n(i)
          endif
       enddo
    enddo
    !}
    

    !> MAP DT_T, SHFLX, and LHFLX FIELDS IN LAND_ICE_ATMOS_BOUNDARY FROM THE EXCHANGE GRID TO THE ATMOSPERE GRID
    !{
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)
#ifndef use_AM3_physics
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%shflx,'ATM', ex_flux_t    , xmap_sfc) !miz
    call fms_xgrid_get_from_xgrid (Land_Ice_Atmos_Boundary%lhflx,'ATM', ex_flux_tr(:,isphum), xmap_sfc)!miz
#endif
    !}
    

    !> MAP DATA FROM THE EXCHANGE GRID TO OCN/ATM/LND GRID AND SEND DATA TO THE DIAG_MANAGER BUFFER.
    !! NOTE, DATA WILL ONLY BE OUTPUTTED IF VARIABLE SPECIFICATION IS FOUND IN DIAG_TABLE.YAML
    !{
    !-------------------- diagnostics section ------------------------------

    !------- new surface temperature -----------
#ifdef use_AM3_physics
    if ( id_t_surf > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
       used = fms_diag_send_data ( id_t_surf, diag_atm, Time )
    endif
#else
    call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
    if ( id_t_surf > 0 ) then
       used = fms_diag_send_data ( id_t_surf, diag_atm, Time )
    endif
    if ( id_ts > 0 ) then
       used = fms_diag_send_data ( id_ts, diag_atm, Time )
    endif
    call fms_sum_diag_integral_field ('t_surf', diag_atm)
#ifndef use_AM3_physics
    if ( id_ts_g > 0 ) used = send_global_diag ( id_ts_g, diag_atm, Time )
#endif
    !------- new surface temperature only over open ocean -----------
    if ( id_tos > 0 ) then
       ex_icetemp = 0.0
       icegrid = 0.0; icegrid(:,:,1) = 1.0
       call fms_xgrid_put_to_xgrid ( icegrid, 'OCN', ex_icetemp, xmap_sfc)
       ex_temp = ex_t_surf_new * ex_icetemp
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_temp, xmap_sfc)
       call fms_xgrid_get_from_xgrid (frac_atm, 'ATM', ex_icetemp, xmap_sfc)
       where (frac_atm > 0.0)
          diag_atm = (diag_atm/frac_atm) ! - tfreeze  CMIP6 in degK
          frac_atm = 1.0
       elsewhere
          diag_atm = 0.0
          frac_atm = 0.0
       endwhere
       used = fms_diag_send_data ( id_tos, diag_atm, Time, rmask=frac_atm )
    endif

    !------- new surface temperature only over land and sea-ice -----------
    if ( id_tslsi > 0 ) then
       ex_land_frac = 0.0
       call put_logical_to_real (Land%mask, 'LND', ex_land_frac, xmap_sfc)
       icegrid = 1.0; icegrid(:,:,1) = 0.
       ex_icetemp = 0.
       call fms_xgrid_put_to_xgrid (icegrid, 'OCN', ex_icetemp, xmap_sfc)
       ex_icetemp = ex_icetemp + ex_land_frac
       ex_temp = ex_t_surf_new * ex_icetemp
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_temp, xmap_sfc)
       call fms_xgrid_get_from_xgrid (frac_atm, 'ATM', ex_icetemp, xmap_sfc)
       where (frac_atm > 0.0)
          diag_atm = diag_atm/frac_atm
          frac_atm = 1.0
       elsewhere
          diag_atm = 0.0
          frac_atm = 0.0
       endwhere
       used = fms_diag_send_data ( id_tslsi, diag_atm, Time, rmask=frac_atm )
    endif
#endif


    ! + slm, Mar 27 2002
    ! ------ new canopy temperature --------
    !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
    !   but this will be changed in future version of the land madel
    if ( id_t_ca > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_t_ca_new, xmap_sfc)
       used = fms_diag_send_data ( id_t_ca, diag_atm, Time )
    endif

    !------- updated surface tracer fields ------
    do tr=1,n_exch_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       if ( id_tr_surf(tr) > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_tr_surf_new(:,tr), xmap_sfc)
          used = fms_diag_send_data ( id_tr_surf(tr), diag_atm, Time )
       endif
       ! - slm Mar 25, 2010: moved to resolve interdependence of diagnostic fields
       if ( id_co2_surf_dvmr > 0 .and. fms_mpp_lowercase(trim(tr_name))=='co2') then
          ex_co2_surf_dvmr = (ex_tr_surf_new(:,tr) / (1.0 - ex_tr_surf_new(:,isphum))) * WTMAIR/WTMCO2
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_co2_surf_dvmr, xmap_sfc)
          used = fms_diag_send_data ( id_co2_surf_dvmr, diag_atm, Time )
       endif
    enddo

    !------- sensible heat flux -----------
    if ( id_t_flux > 0 .or. id_hfss > 0 .or. id_hfss_g > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
       if ( id_t_flux > 0 ) used = fms_diag_send_data ( id_t_flux, diag_atm, Time )
       if ( id_hfss   > 0 ) used = fms_diag_send_data ( id_hfss, diag_atm, Time )
#ifndef use_AM3_physics
       if ( id_hfss_g > 0 ) used = send_global_diag ( id_hfss_g, diag_atm, Time )
#endif
    endif

    !------- net longwave flux -----------
    if ( id_r_flux > 0 .or. id_rls_g > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
       if ( id_r_flux > 0 ) used = fms_diag_send_data ( id_r_flux, diag_atm, Time )
#ifndef use_AM3_physics
       if ( id_rls_g  > 0 ) used = send_global_diag ( id_rls_g, diag_atm, Time )
#endif
    endif

    !------- tracer fluxes ------------
    ! tr_mol_flux diagnostic will be correct for co2 tracer only.
    ! will need update code to use correct molar mass for tracers other than co2
    do tr=1,n_exch_tr
       if ( id_tr_flux(tr) > 0 .or. id_tr_mol_flux(tr) > 0 ) then
          call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_flux_tr(:,tr), xmap_sfc)
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name, units=tr_units )
          if (id_tr_flux(tr) > 0 ) &
               used = fms_diag_send_data ( id_tr_flux(tr), diag_atm, Time )
          !     if (id_tr_mol_flux(tr) > 0 ) &
          !          used = fms_diag_end_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMCO2, Time)
          ! 2017/08/08 jgj - replaced 2 lines above by the following
          if (id_tr_mol_flux(tr) > 0 .and. fms_mpp_lowercase(trim(tr_name))=='co2') then
             used = fms_diag_send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMCO2, Time)
             !sometimes in 2018 f1p for vmr tracers
            elseif (id_tr_mol_flux(tr) > 0 .and. fms_mpp_lowercase(trim(tr_units)).eq."vmr") then
               ! if (ocn_atm_flux_vmr_bug) then
               !    call fms_xgrid_get_from_xgrid (diag_atm, 'ATM',  &
               !         ex_flux_tr(:,tr)*(1.-ex_tr_surf_new(:,isphum)), xmap_sfc)
               !    used = fms_diag_send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMAIR, Time)
               ! else
               !flux is in vmr * kg/m2/s. Divide by MW_air
               call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', &
                    ex_flux_tr(:,tr)*((1.-ex_tr_surf_new(:,isphum))*WTMH2O+ex_tr_surf_new(:,isphum)*WTMAIR) &
                    / (1e-3*WTMAIR*WTMH2O) , &
                    xmap_sfc)
               used = fms_diag_send_data ( id_tr_mol_flux(tr), diag_atm, Time)
           endif
        endif
        if ( id_tr_con_atm(tr) > 0 ) then
           call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_tr_con_atm(:,tr), xmap_sfc)
           used = fms_diag_send_data ( id_tr_con_atm(tr), diag_atm, Time )
        end if
        if ( id_tr_con_ref(tr) > 0 ) then
           call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_tr_con_ref(:,tr), xmap_sfc)
           used = fms_diag_send_data ( id_tr_con_ref(tr), diag_atm, Time )
        end if
    enddo

#ifndef _USE_LEGACY_LAND_
    if ( id_t_flux_land > 0 ) then
       call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_flux_t, xmap_sfc)
       call send_tile_data ( id_t_flux_land, diag_land )
    endif
    !------- tracer fluxes for land
    do tr=1,n_exch_tr
       if ( id_tr_flux_land(tr) > 0 .or. id_tr_mol_flux_land(tr) > 0 ) then
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name, units=tr_units )
          call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc)
          if (id_tr_flux_land(tr) > 0 ) &
                 call send_tile_data (id_tr_flux_land(tr), diag_land )
          if (id_tr_mol_flux_land(tr) > 0) then
             if (fms_mpp_lowercase(trim(tr_name))=='co2') then
                call send_tile_data (id_tr_mol_flux_land(tr), diag_land*1000./WTMCO2)
             elseif (fms_mpp_lowercase(trim(tr_units)).eq.'vmr') then
                !flux is in vmr * kg/m2/s. Divide by MW_air
                call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', &
                     ex_flux_tr(:,tr)*((1.-ex_tr_surf_new(:,isphum))*WTMH2O+ex_tr_surf_new(:,isphum)*WTMAIR) &
                     /(1e-3*WTMAIR*WTMH2O) , &
                     xmap_sfc)
                call send_tile_data ( id_tr_mol_flux_land(tr), diag_land)
             endif
          endif
       endif
    enddo

    !-------- tracer deposition velocity
    do tr=1,n_exch_tr
       if ( id_tr_con_atm_land(tr) > 0 ) then
          call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_tr_con_atm(:,tr), xmap_sfc)
          call send_tile_data (id_tr_con_atm_land(tr), diag_land )
       endif
       if ( id_tr_con_ref_land(tr) > 0 ) then
          call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_tr_con_ref(:,tr), xmap_sfc )
          call send_tile_data (id_tr_con_ref_land(tr), diag_land )
       endif
    enddo

#endif

    !-----------------------------------------------------------------------
    !---- accumulate global integral of evaporation (mm/day) -----
    call fms_xgrid_get_from_xgrid (evap_atm, 'ATM', ex_flux_tr(:,isphum), xmap_sfc)
    if( id_q_flux > 0 )  used = fms_diag_send_data ( id_q_flux, evap_atm, Time)
    if( id_evspsbl > 0 ) used = fms_diag_send_data ( id_evspsbl, evap_atm, Time)
    if( id_hfls    > 0 ) used = fms_diag_send_data ( id_hfls, HLV*evap_atm, Time)
#ifndef use_AM3_physics
    if( id_hfls_g  > 0 ) used = send_global_diag ( id_hfls_g, HLV*evap_atm, Time)
#endif

#ifndef _USE_LEGACY_LAND_
    if( id_q_flux_land > 0 ) then
       call fms_xgrid_get_from_xgrid_ug (diag_land, 'LND', ex_flux_tr(:,isphum), xmap_sfc)
       call send_tile_data (id_q_flux_land, diag_land)
#else
    if( id_q_flux_land > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_land, 'LND', ex_flux_tr(:,isphum), xmap_sfc)
       used = fms_diag_send_tile_averaged_data(id_q_flux_land, diag_land, &
            Land%tile_size, Time, mask=Land%mask)
#endif
    endif
    call fms_sum_diag_integral_field ('evap', evap_atm*86400.)
#ifndef use_AM3_physics
    if (id_evspsbl_g > 0) used = send_global_diag ( id_evspsbl_g, evap_atm, Time )
#endif

#ifndef _USE_LEGACY_LAND_
    call send_tile_data (id_q_flux_land, diag_land)
    ! need this to avoid diag issues with tiling changes in update_land_slow
    call dump_tile_diag_fields(Time)
    !}


    !>  COMPUTE STOCK EXCHANGE BETWEEN MODEL COMPONENTS
    !{
    call fms_xgrid_get_from_xgrid_ug(data_lnd, 'LND', ex_flux_tr(:,isphum), xmap_sfc)    
    
    ! Lnd -> Atm (evap)
    call fms_xgrid_stock_move_ug( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_WATER), &
         & FROM = fms_stock_constants_lnd_stock(ISTOCK_WATER), &
         & stock_ug_data3d = data_lnd, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Lnd->ATm) ')

    ! Lnd -> Atm (heat lost through evap)
    call fms_xgrid_stock_move_ug( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_HEAT), &
         & FROM = fms_stock_constants_lnd_stock(ISTOCK_HEAT), &
         & stock_ug_data3d = data_lnd * HLV, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Lnd->ATm) ')
#else
    call fms_xgrid_get_from_xgrid(data_lnd, 'LND', ex_flux_tr(:,isphum), xmap_sfc)

    ! Lnd -> Atm (evap)
    call fms_xgrid_stock_move( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_WATER), &
         & FROM = fms_stock_constants_lnd_stock(ISTOCK_WATER), &
         & stock_data3d = data_lnd, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Lnd->ATm) ')

    ! Lnd -> Atm (heat lost through evap)
    call fms_xgrid_stock_move( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_HEAT), &
         & FROM = fms_stock_constants_lnd_stock(ISTOCK_HEAT), &
         & stock_data3d = data_lnd * HLV, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Lnd->ATm) ')
#endif

    call fms_xgrid_get_from_xgrid(data_ice, 'OCN', ex_flux_tr(:,isphum), xmap_sfc)

    ! Ice -> Atm (evap)
    call fms_xgrid_stock_move( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_WATER), &
         & FROM = fms_stock_constants_ice_stock(ISTOCK_WATER), &
         & stock_data3d = data_ice, &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Ice->ATm) ')

    ! Ice -> Atm (heat lost through evap)
    call fms_xgrid_stock_move( &
         & TO   = fms_stock_constants_atm_stock(ISTOCK_HEAT), &
         & FROM = fms_stock_constants_ice_stock(ISTOCK_HEAT), &
         & stock_data3d = data_ice * HLV, &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Ice->ATm) ')
    !}
   

    !> END CLOCK FOR PROFILING
    !{
    call fms_mpp_clock_end(fluxAtmUpClock)
    call fms_mpp_clock_end(cplClock)
    !}
    
  end subroutine flux_up_to_atmos

  !> Subroutine flux_ex_arrays_dealloc deallocates the model-level ex_* arrays that were
  !! allocated in sfc_boundary_layer
  subroutine flux_ex_arrays_dealloc()

    integer :: m,n

    deallocate ( &
         ex_t_surf   ,  &
         ex_t_surf_miz, &
         ex_p_surf   ,  &
         ex_slp      ,  &
         ex_t_ca     ,  &
         ex_dhdt_surf,  &
         ex_dedt_surf,  &
         ex_dqsatdt_surf,  &
         ex_drdt_surf,  &
         ex_dhdt_atm ,  &
         ex_flux_t   ,  &
         ex_flux_lw  ,  &
         ex_drag_q   ,  &
         ex_avail    ,  &
         ex_f_t_delt_n, &
         ex_tr_surf  ,  &
         ex_tr_con_ref, &
         ex_tr_con_atm, &

         ex_dfdtr_surf  , &
         ex_dfdtr_atm   , &
         ex_flux_tr     , &
         ex_f_tr_delt_n , &
         ex_e_tr_n      , &

         ex_e_t_n    ,  &
         ex_e_q_n    ,  &
                                ! values added for LM3
         ex_cd_t     ,  &
         ex_cd_m     ,  &
         ex_b_star   ,  &
         ex_u_star   ,  &
         ex_wind     ,  &
         ex_z_atm    ,  &
         ex_con_atm  ,  &
         ex_seawater ,  &
         ex_land        )

#ifdef SCM
    deallocate ( &
         ex_dhdt_surf_forland, &
         ex_dedt_surf_forland, &
         ex_dedq_surf_forland  )
#endif

    ! Extra fluxes
    do n = 1, ex_gas_fields_ice%num_bcs
       do m = 1, ex_gas_fields_ice%bc(n)%num_fields
          deallocate ( ex_gas_fields_ice%bc(n)%field(m)%values )
          nullify ( ex_gas_fields_ice%bc(n)%field(m)%values )
       enddo 
    enddo 

    do n = 1, ex_gas_fields_atm%num_bcs
       do m = 1, ex_gas_fields_atm%bc(n)%num_fields
          deallocate ( ex_gas_fields_atm%bc(n)%field(m)%values )
          nullify ( ex_gas_fields_atm%bc(n)%field(m)%values )
       enddo
    enddo
    
    do n = 1, ex_gas_fluxes%num_bcs
       do m = 1, ex_gas_fluxes%bc(n)%num_fields
          deallocate ( ex_gas_fluxes%bc(n)%field(m)%values )
          nullify ( ex_gas_fluxes%bc(n)%field(m)%values )
       enddo 
    enddo 

  end subroutine flux_ex_arrays_dealloc


  !> Subroutine flux_atmos_to_ocean computes deposition gas fluxes between atmosphere and ocean
  !! This subroutine is called only if the do_flux namelist variable is set to .True.
  subroutine flux_atmos_to_ocean(Time, Atm, Ice_boundary, Ice)

    !> current time
    type(FmsTime_type), intent(in) :: Time
    !> derived data type holding atmosphere boundary data
    type(atmos_data_type), intent(inout):: Atm
    !> derived data type holding properties and fluxes passed from atmosphere to ice
    type(atmos_ice_boundary_type), intent(inout):: Ice_boundary
    !> derived type holding ice boundary tdata
    type(ice_data_type),  intent(inout):: Ice
    
    integer :: n, m
    logical :: used

#ifndef use_AM3_physics
    call atmos_tracer_driver_gather_data_down(Atm%fields, Atm%tr_bot)
#endif

    !> MAP ATMOSPHERE FIELDS TO THE EXCHANGE MAP FOR FLUX EXCHANGE WITH OCEAN
    !{
    !air-sea deposition fluxes
    do n = 1, Atm%fields%num_bcs
       !Do the string copies.
       Atm%fields%bc(n)%flux_type = trim(ex_gas_fluxes%bc(n)%flux_type)
       Atm%fields%bc(n)%implementation = trim(ex_gas_fluxes%bc(n)%implementation)
       if(ex_gas_fields_atm%bc(n)%flux_type  .eq. 'air_sea_deposition') then
          do m = 1, Atm%fields%bc(n)%num_fields 
             call fms_xgrid_put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
                  ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method)
          enddo
       endif
    enddo
    !}

    
    !> ON THE EXCHANGE GRID, CALCULATE OCEAN EXPLICIT FLUX BY CALLING ATMOS_OCEAN_DEP_FLUXES_CALC
    !{
    call atmos_ocean_dep_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater)
    !}


    !> MAP AIR_SEA_DEPOSITION FLUX FROM THE EXCHANGE GRID TO THE ICE GRID FOLLOWED BY
    !! CALL DATA_OVERRIDE WHERE DATA WILL BE OVERWRITTEN IF THE FLUX FIELDS
    !! ARE SPECIFIED IN THE DATA_TABLE.  SEND_DATA TO THE DIAG_MANAGER BUFFER
    !{
    do n = 1, Ice_boundary%fluxes%num_bcs
       if(Ice_boundary%fluxes%bc(n)%flux_type  .eq. 'air_sea_deposition') then
          do m = 1, Ice_boundary%fluxes%bc(n)%num_fields
             call fms_xgrid_get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
                  ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)
             
             call fms_data_override('ICE', Ice_boundary%fluxes%bc(n)%field(m)%name,     &
                  Ice_boundary%fluxes%bc(n)%field(m)%values, Time)
             if ( Ice_boundary%fluxes%bc(n)%field(m)%id_diag > 0 ) then
                used = fms_diag_send_data(Ice_boundary%fluxes%bc(n)%field(m)%id_diag, &
                     Ice_boundary%fluxes%bc(n)%field(m)%values, Time )
             endif
          enddo 
       endif
    enddo
    !}
    

    !> UPDATE ICE FIELDS THAT ARE LABELED AS AIR_SEA_DEOOSITION FLUXES BY CALLING UPDATE_ICE_ATM_DEPOSITION_FLUX
    !{
    call update_ice_atm_deposition_flux( Ice_boundary, Ice )
    !}
    
  end subroutine flux_atmos_to_ocean

  !> Subroutine put_logical_to_real_sg maps 2D logical mask arrays to real arrays
  !! where .true. -> 1.0 and .false. -> 0.0.  The real array is then mapped
  !! onto the exchange grid.   This subroutine is used internally to convert Land%mask
  !! on structured grid for example, when #ifndef _USE_LEGACY_LAND_ is false  
  subroutine put_logical_to_real_sg (mask, id, ex_mask, xmap)

    !> land/ice mask
    logical, intent(in) :: mask(:,:,:)
    !> component id
    character(len=3), intent(in) :: id
    !> mapped  mask on exchange grid
    real, intent(inout) :: ex_mask(:)
    type(FmsXgridXmap_type), intent(inout) :: xmap

    real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask

    !> MAP (LOGICAL) MASK TO RMASK
    !{
    where (mask)
       rmask = 1.0
    elsewhere
       rmask = 0.0
    endwhere
    !}

    
    !> MAP RMASK TO THE EXCHANGE GRID
    !{
    call fms_xgrid_put_to_xgrid(rmask, id, ex_mask, xmap)
    !}

  end subroutine put_logical_to_real_sg


  !> Subroutine put_logical_to_real_ug maps 2D logical mask arrays to real arrays
  !! where .true. -> 1.0 and .false. -> 0.0.  The real array is then mapped
  !! onto the exchange grid.   This subroutine is used internally to convert Land%mask
  !! on unstructured grid (when #ifndef _USE_LEGACY_LAND_ is true)
  subroutine put_logical_to_real_ug (mask, id, ex_mask, xmap)

    !> mask on component grid
    logical, intent(in) :: mask(:,:)
    !> component id
    character(len=3), intent(in) :: id
    !> mapped mask on exchange grid
    real, intent(inout) :: ex_mask(:)    
    type(FmsXgridXmap_type), intent(inout) :: xmap

    real, dimension(size(mask,1),size(mask,2)) :: rmask

    !> MAP (LOGICAL) MASK TO RMASK
    !{
    where (mask)
       rmask = 1.0
    elsewhere
       rmask = 0.0
    endwhere
    !}


    !> MAP RMASK TO THE EXCHANGE GRID
#ifndef _USE_LEGACY_LAND_
    call fms_xgrid_put_to_xgrid_ug(rmask, id, ex_mask, xmap)
#else
    call fms_xgrid_put_to_xgrid (rmask, id, ex_mask, xmap)
#endif
    !}
    
  end subroutine put_logical_to_real_ug


  !> Subroutine diag_field_init registers the diagnostic fields in this module to the diag_manager
  !! Note, diagnostic fields must be registered in diag_manager and all diagnostics fields
  !! must be specified in the diag_table in order for the data to be outputted
  !! to a NetCDF file at the end of the model run.  This subroutine is called
  !! during module initialization in subroutine atm_land_ice_flux_exchange_init
  subroutine diag_field_init ( Time, atmos_axes, land_axes, land_pe )

    !> curent model time
    type(FmsTime_type), intent(in) :: Time
    !> array size for atmospheric diagnostic fields
    integer, intent(in) :: atmos_axes(2)
    !> array size for land diagnostic fields
    integer, intent(in) :: land_axes(:)
    !> land pe number
    logical, intent(in) :: land_pe

    integer :: iref
    character(len=6) :: label_zm, label_zh
    real, dimension(2) :: trange = (/  100., 400. /), &
         vrange = (/ -400., 400. /), &
         frange = (/ -0.01, 1.01 /)
    character(len=32)  :: name, units ! name of the tracer
    character(len=128) :: longname    ! long name of the tracer
    integer            :: tr          ! tracer index
    integer            :: area_id

    !> CONVERT DIAGNOSTIC LABELS FROM INTEGERS TO STRINGS
    !{
    !------ labels for diagnostics -------
    !  (z_ref_mom, z_ref_heat are namelist variables)
    iref = int(z_ref_mom+0.5)
    if ( real(iref) == z_ref_mom ) then
       write (label_zm,105) iref
       if (iref < 10) write (label_zm,100) iref
    else
       write (label_zm,110) z_ref_mom
    endif

    iref = int(z_ref_heat+0.5)
    if ( real(iref) == z_ref_heat ) then
       write (label_zh,105) iref
       if (iref < 10) write (label_zh,100) iref
    else
       write (label_zh,110) z_ref_heat
    endif
    
100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')
    !}


    !> CALL FMS_DIAG_REGISTER_DIAG_FIELD
    !{
    !--------- initialize static diagnostic fields --------------------

    id_land_mask = &
         fms_diag_register_static_field ( mod_name, 'land_mask', atmos_axes,  &
         'fractional amount of land', 'none', &
         range=frange, interp_method = "conserve_order1" )

    !--------- initialize diagnostic fields --------------------

    id_ice_mask = &
         fms_diag_register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
         'fractional amount of sea ice', 'none',  &
         range=frange, interp_method = "conserve_order1" )

    id_wind = &
         fms_diag_register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
         'wind speed for flux calculations', 'm/s', &
         range=(/0.,vrange(2)/) )

    id_drag_moist = &
         fms_diag_register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
         'drag coeff for moisture',    'none'     )

    id_drag_heat  = &
         fms_diag_register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
         'drag coeff for heat',    'none'     )

    id_drag_mom   = &
         fms_diag_register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
         'drag coeff for momentum',     'none'     )

    id_rough_moist = &
         fms_diag_register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
         'surface roughness for moisture',  'm'  )

    id_rough_heat = &
         fms_diag_register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
         'surface roughness for heat',  'm'  )

    id_rough_mom  = &
         fms_diag_register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
         'surface roughness for momentum',  'm'  )

    id_u_star     = &
         fms_diag_register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
         'friction velocity',   'm/s'   )

    id_b_star     = &
         fms_diag_register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
         'buoyancy scale',      'm/s2'   )

    id_q_star     = &
         fms_diag_register_diag_field ( mod_name, 'q_star',     atmos_axes, Time, &
         'moisture scale',      'kg water/kg air'   )

    id_thv_atm = &
         fms_diag_register_diag_field ( mod_name, 'thv_atm', atmos_axes, Time, &
         'surface air virtual potential temperature', 'K')

    id_thv_surf = &
         fms_diag_register_diag_field ( mod_name, 'thv_surf', atmos_axes, Time, &
         'surface virtual potential temperature', 'K')

    id_u_flux     = &
         fms_diag_register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
         'zonal wind stress',     'pa'   )

    id_v_flux     = &
         fms_diag_register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
         'meridional wind stress',     'pa'   )

    id_t_ocean     = &
         fms_diag_register_diag_field ( mod_name, 't_ocean',     atmos_axes, Time, &
         'surface temperature from ocean output',    'deg_k', &
         range=trange    )

    id_t_surf     = &
         fms_diag_register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
         'surface temperature',    'deg_k', &
         range=trange    )

    ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
    id_t_ca       = &
         fms_diag_register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
         'canopy air temperature',    'deg_k', &
         range=trange    )

    ! - slm, Mar 25, 2002
    id_z_atm      = &
         fms_diag_register_diag_field ( mod_name, 'z_atm',     atmos_axes, Time, &
         'height of btm level',    'm')

    id_p_atm      = &
         fms_diag_register_diag_field ( mod_name, 'p_atm',     atmos_axes, Time, &
         'pressure at btm level',    'pa')

    ! - bw, Mar 25, 2002 -- added diagnostic slp
    id_slp      = &
         fms_diag_register_diag_field ( mod_name, 'slp',      atmos_axes, Time, &
         'sea level pressure',    'pa')

    id_gust       = &
         fms_diag_register_diag_field ( mod_name, 'gust',     atmos_axes, Time, &
         'gust scale',    'm/s')

    id_t_flux     = &
         fms_diag_register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
         'sensible heat flux',     'w/m2'    )

    id_r_flux     = &
         fms_diag_register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
         'net (down-up) longwave flux',   'w/m2'    )

    id_t_atm      = &
         fms_diag_register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
         'temperature at btm level',    'deg_k', &
         range=trange     )

    id_u_atm      = &
         fms_diag_register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
         'u wind component at btm level',  'm/s', &
         range=vrange    )

    id_v_atm      = &
         fms_diag_register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
         'v wind component at btm level',  'm/s', &
         range=vrange    )

    id_t_ref      = &
         fms_diag_register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
         'temperature at '//label_zh, 'deg_k' , &
         range=trange      )

    id_rh_ref     = &
         fms_diag_register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
         'relative humidity at '//label_zh, 'percent' )

    id_rh_ref_cmip = &
         fms_diag_register_diag_field ( mod_name, 'rh_ref_cmip',     atmos_axes, Time,   &
         'relative humidity at '//label_zh, 'percent' )

    id_u_ref      = &
         fms_diag_register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
         'zonal wind component at '//label_zm,  'm/s', &
         range=vrange )

    id_v_ref      = &
         fms_diag_register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
         'meridional wind component at '//label_zm, 'm/s', &
         range=vrange )

    id_wind_ref = &
         fms_diag_register_diag_field ( mod_name, 'wind_ref',   atmos_axes, Time,     &
         'absolute value of wind at '//label_zm, 'm/s', &
         range=vrange )

    id_del_h      = &
         fms_diag_register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
         'ref height interp factor for heat', 'none' )
    id_del_m      = &
         fms_diag_register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
         'ref height interp factor for momentum','none' )
    id_del_q      = &
         fms_diag_register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
         'ref height interp factor for moisture','none' )

    if( land_pe ) then
       ! set the default filter (for area and subsampling) for consequent calls to
       ! register_tiled_diag_field
#ifndef _USE_LEGACY_LAND_
       call set_default_diag_filter('land')
       id_t_ref_land = &
            register_tiled_diag_field ( 'flux_land', 't_ref', Land_axes, Time, &
            'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
            range=trange, missing_value =  -100.0)
       id_q_ref_land = &
            register_tiled_diag_field ( 'flux_land', 'q_ref', Land_axes, Time, &
            'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
            missing_value=-1.0)
       id_rh_ref_land= &
            register_tiled_diag_field ( 'flux_land', 'rh_ref', Land_axes, Time,   &
            'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
            missing_value=-999.0)
       id_u_ref_land = &
            register_tiled_diag_field ( 'flux_land', 'u_ref',  Land_axes, Time, &
            'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
            range=vrange, missing_value=-999.0 )
       id_v_ref_land = &
            register_tiled_diag_field ( 'flux_land', 'v_ref',  Land_axes, Time,     &
            'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
            range=vrange, missing_value = -999.0 )
       id_q_flux_land = &
            register_tiled_diag_field( 'flux_land', 'evap', Land_axes, Time, &
            'evaporation rate over land', 'kg/m2/s', missing_value=-1.0 )
       id_t_flux_land = &
            register_tiled_diag_field( 'flux_land', 'shflx', Land_axes, Time, &
            'sensible heat flux', 'W/m2', missing_value=-1.0 )
       id_tasLut_land = &
            register_tiled_diag_field( 'cmor_land', 'tasLut', Land_axes, Time, &
            'Near-Surface Air Temperature ('//trim(label_zh)//' Above Displacement Height) on Land Use Tile', &
            units='K', standard_name='air_temperature', missing_value=-1.0 )
       id_hussLut_land = &
            register_tiled_diag_field( 'cmor_land', 'hussLut', Land_axes, Time, &
            'Near-Surface Specific Humidity on Land Use Tile', '1.0', &
            standard_name='specific_humidity', missing_value=-1.0 )
       allocate(id_tr_flux_land(n_exch_tr))
       allocate(id_tr_mol_flux_land(n_exch_tr))
       allocate(id_tr_con_atm_land(n_exch_tr))
       allocate(id_tr_con_ref_land(n_exch_tr))
       allocate(id_tr_ref_land(n_exch_tr))

       do tr = 1, n_exch_tr
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )

          id_tr_con_atm_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_tot_con_atm', &
               Land_axes, Time, 'vd of '//trim(longname), 'm/s', missing_value=-1.0 )
          id_tr_con_ref_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_tot_con_ref', &
               Land_axes, Time, 'vd of '//trim(longname)//' at '//trim(label_zh), 'm/s', missing_value=-1.0 )

          id_tr_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_flux', &
               Land_axes, Time, 'flux of '//trim(longname), trim(units)//' kg air/(m2 s)', missing_value=-1.0 )
          if ( fms_mpp_lowercase(trim(name))=='co2') then
             id_tr_mol_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_mol_flux', &
                  Land_axes,Time, 'flux of '//trim(longname), 'mol CO2/(m2 s)', missing_value=-1.0 )
          else
             id_tr_mol_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_mol_flux', &
                  Land_axes,Time, 'flux of '//trim(longname), 'mol/(m2 s)', missing_value=-1.0 )
          endif
          ! we skip sphum because it is already available as flux_land/q_ref
          if ( tr .ne. isphum ) then
             id_tr_ref_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_ref', &
                  Land_axes, Time, trim(longname)//' at '//trim(label_zh)//' over land', &
                  trim(units),missing_value=-1.0)
          else
             id_tr_ref_land(tr) = -1
          end if
       enddo
#else
       id_t_ref_land = &
            fms_diag_register_diag_field ( 'flux_land', 't_ref', Land_axes, Time, &
            'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
            range=trange, missing_value =  -100.0)
       id_q_ref_land = &
            fms_diag_register_diag_field ( 'flux_land', 'q_ref', Land_axes, Time, &
            'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
            missing_value=-1.0)
       id_rh_ref_land= &
            fms_diag_register_diag_field ( 'flux_land', 'rh_ref', Land_axes, Time,   &
            'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
            missing_value=-999.0)
       id_u_ref_land = &
            fms_diag_register_diag_field ( 'flux_land', 'u_ref',  Land_axes, Time, &
            'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
            range=vrange, missing_value=-999.0 )
       id_v_ref_land = &
            fms_diag_register_diag_field ( 'flux_land', 'v_ref',  Land_axes, Time,     &
            'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
            range=vrange, missing_value = -999.0 )
       id_q_flux_land = &
            fms_diag_register_diag_field( 'flux_land', 'evap', Land_axes, Time, &
            'evaporation rate over land', 'kg/m2/s', missing_value=-1.0 )
       id_t_flux_land = &
            fms_diag_register_diag_field( 'flux_land', 'shflx', Land_axes, Time, &
            'sensible heat flux', 'W/m2', missing_value=-1.0 )
       id_tasLut_land = &
            fms_diag_register_diag_field( 'cmor_land', 'tasLut', Land_axes, Time, &
            'Near-Surface Air Temperature ('//trim(label_zh)//' Above Displacement Height) on Land Use Tile', &
            units='K', standard_name='air_temperature', missing_value=-1.0 )
       id_hussLut_land = &
            fms_diag_register_diag_field( 'cmor_land', 'hussLut', Land_axes, Time, &
            'Near-Surface Specific Humidity on Land Use Tile', '1.0', &
            standard_name='specific_humidity', missing_value=-1.0 )
       allocate(id_tr_flux_land(n_exch_tr))
       allocate(id_tr_mol_flux_land(n_exch_tr))
       allocate(id_tr_con_atm_land(n_exch_tr)); id_tr_con_atm_land(:)=-1
       allocate(id_tr_con_ref_land(n_exch_tr)); id_tr_con_ref_land(:)=-1
       allocate(id_tr_ref_land(n_exch_tr)); id_tr_ref_land(:)=-1

       do tr = 1, n_exch_tr
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
          id_tr_flux_land(tr) = fms_diag_register_diag_field( 'flux_land', trim(name)//'_flux', Land_axes, Time, &
               'flux of '//trim(longname), trim(units)//' kg air/(m2 s)', missing_value=-1.0 )
          if ( fms_mpp_lowercase(trim(name))=='co2') then
             id_tr_mol_flux_land(tr) = fms_diag_register_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, &
                  Time, 'flux of '//trim(longname), 'mol CO2/(m2 s)', missing_value=-1.0 )
          else
             id_tr_mol_flux_land(tr) = fms_diag_register_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, &
                  Time, 'flux of '//trim(longname), 'mol/(m2 s)', missing_value=-1.0 )
          endif
       enddo
#endif
    endif

    id_q_ref = &
         fms_diag_register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
         'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)

    id_rough_scale = &
         fms_diag_register_diag_field ( mod_name, 'rough_scale', atmos_axes, Time, &
         'topographic scaling factor for momentum drag','1' )
    !-----------------------------------------------------------------------

    allocate(id_tr_atm(n_exch_tr))
    allocate(id_tr_surf(n_exch_tr))
    allocate(id_tr_flux(n_exch_tr))
    allocate(id_tr_mol_flux(n_exch_tr))
    allocate(id_tr_mol_flux0(n_exch_tr))
    allocate(id_tr_con_atm(n_exch_tr))
    allocate(id_tr_con_ref(n_exch_tr))
    allocate(id_tr_ref(n_exch_tr))

    do tr = 1, n_exch_tr
       call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )

       id_tr_con_atm(tr) = fms_diag_register_diag_field( mod_name, trim(name)//'_tot_con_atm', atmos_axes, Time, &
            'vd of '//trim(longname), 'm/s', missing_value=-1.0 )
       id_tr_con_ref(tr) = fms_diag_register_diag_field( mod_name, trim(name)//'_tot_con_ref', atmos_axes, Time, &
            'vd of '//trim(longname)//' at '//trim(label_zh), 'm/s', missing_value=-1.0 )

       id_tr_atm(tr) = fms_diag_register_diag_field (mod_name, trim(name)//'_atm', atmos_axes, Time, &
            trim(longname)//' at btm level', trim(units))
       id_tr_surf(tr) = fms_diag_register_diag_field (mod_name, trim(name)//'_surf', atmos_axes, Time, &
            trim(longname)//' at the surface', trim(units))
       id_tr_flux(tr) = fms_diag_register_diag_field(mod_name, trim(name)//'_flux', atmos_axes, Time, &
            'flux of '//trim(longname), trim(units)//' kg air/(m2 s)')
       if ( tr .ne. isphum ) then
          id_tr_ref(tr) = fms_diag_register_diag_field (mod_name, trim(name)//'_ref',  atmos_axes, Time, &
               trim(longname)//' at '//trim(label_zh), trim(units),missing_value=-1.0)
       else
          id_tr_ref(tr) = -1
       end if
       !! add dryvmr co2_surf and co2_atm
       if ( fms_mpp_lowercase(trim(name))=='co2') then
          ! - slm Mar 25, 2010: moved registration of mol_flux inside 'if' to disable
          ! saving incorrect results (mol fluxes for other tracers computed with CO2 molar
          ! mass)
          id_tr_mol_flux(tr) = fms_diag_register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
               'flux of '//trim(longname), 'mol CO2/(m2 s)')
          id_co2_atm_dvmr = fms_diag_register_diag_field (mod_name, trim(name)//'_atm_dvmr', atmos_axes, Time, &
               trim(longname)//' at btm level', 'mol CO2 /mol air')
          id_co2_surf_dvmr = fms_diag_register_diag_field (mod_name, trim(name)//'_surf_dvmr', atmos_axes, Time, &
               trim(longname)//' at the surface', 'mol CO2 /mol air')
       else
!f1p
          id_tr_mol_flux(tr) = fms_diag_register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
               'flux of '//trim(longname), 'mol/(m2 s)')
       endif
!f1p
       id_tr_mol_flux0(tr) = fms_diag_register_diag_field(mod_name, trim(name)//'_mol_flux_atm0', atmos_axes, Time, &
            'gross flux of '//trim(longname), 'mol/(m2 s)')

    enddo

    ! 2017/08/08 jgj add diagnostics for co2 data overrides even if co2 is not a tracer
    ! register data calls not needed here for co2_flux_pcair_atm and o2_flux_pcair_atm as this happens elsewhere
    id_co2_bot = fms_diag_register_diag_field (mod_name, 'co2_bot', atmos_axes, Time, &
           'co2_bot from data_override', 'ppmv')

    ! id_nh3_flux_atm0 = fms_diag_register_diag_field (mod_name, 'nh3_flux_atm0', atmos_axes, Time, &
    !        'nh3 flux out of the ocean assuming not nh3 in the atmosphere', 'mol/m2/s')


    id_q_flux = fms_diag_register_diag_field( mod_name, 'evap',       atmos_axes, Time, &
         'evaporation rate',        'kg/m2/s'  )

    !--------------------------------------------------------------------
    !    retrieve the diag_manager id for the area diagnostic,
    !    needed for cmorizing various diagnostics.
    !--------------------------------------------------------------------
    area_id = fms_diag_get_field_id ('dynamics', 'area')
    if (area_id .eq. DIAG_FIELD_NOT_FOUND) call fms_error_mesg &
         ('diag_field_init in atm_land_ice_flux_exchange_mod', &
         'diagnostic field "dynamics", "area" is not in the diag_table', NOTE)

    !-----------------------------------------------------------------------
    !  register cmip variable names
    !-----------------------------------------------------------------------
    ! NOTE: add extra dimension reference level fields?  height2m, height10m
    !       for now we will handle this with an attribute

    id_height2m = &
        fms_diag_register_static_field ( mod_name, 'height2m', (/null_axis_id/), &
                             'Height', 'm', standard_name = 'height' )
    if ( id_height2m > 0 ) then
       call fms_diag_field_add_attribute( id_height2m, 'axis', 'Z' )
       call fms_diag_field_add_attribute( id_height2m, 'positive', 'up' )
    endif

    id_height10m = &
        fms_diag_register_static_field ( mod_name, 'height10m', (/null_axis_id/), &
                             'Height', 'm', standard_name = 'height' )
    if ( id_height10m > 0 ) then
       call fms_diag_field_add_attribute( id_height10m, 'axis', 'Z' )
       call fms_diag_field_add_attribute( id_height10m, 'positive', 'up' )
    endif

#ifdef use_AM3_physics
    id_tas      = &
         fms_diag_register_diag_field ( mod_name, 'tas', atmos_axes, Time, &
         'Near-Surface Air Temperature', 'K' , &
         standard_name = 'air_temperature', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=trange )
    if ( id_tas > 0 .and. id_height2m > 0) &
       call fms_diag_field_add_attribute( id_tas, 'coordinates', 'height2m' )

    id_uas      = &
         fms_diag_register_diag_field ( mod_name, 'uas', atmos_axes, Time, &
         'Eastward Near-Surface Wind', 'm s-1', &
         standard_name = 'eastward_wind', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_uas > 0 .and. id_height10m > 0) &
       call fms_diag_field_add_attribute( id_uas, 'coordinates', 'height10m' )

    id_vas      = &
         fms_diag_register_diag_field ( mod_name, 'vas', atmos_axes, Time, &
         'Northward Near-Surface Wind', 'm s-1', &
         standard_name = 'northward_wind', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_vas > 0 .and. id_height10m > 0 ) &
       call fms_diag_field_add_attribute( id_vas, 'coordinates', 'height10m' )

    id_sfcWind = &
         fms_diag_register_diag_field ( mod_name, 'sfcWind', atmos_axes, Time, &
         'Near-Surface Wind Speed', 'm s-1', &
         standard_name = 'wind_speed', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_sfcWind > 0 .and. id_height10m > 0 ) &
       call fms_diag_field_add_attribute( id_sfcWind, 'coordinates', 'height10m' )

    id_huss = &
         fms_diag_register_diag_field ( mod_name, 'huss', atmos_axes, Time, &
         'Near-Surface Specific Humidity', '1.0', &
         standard_name = 'specific_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_huss > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_huss, 'coordinates', 'height2m' )

    id_hurs = &
         fms_diag_register_diag_field ( mod_name, 'hurs', atmos_axes, Time, &
         'Near-Surface Relative Humidity', '%', &
         standard_name = 'relative_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_hurs > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_hurs, 'coordinates', 'height2m' )

    id_rhs = &
         fms_diag_register_diag_field ( mod_name, 'rhs', atmos_axes, Time, &
         'Near-Surface Relative Humidity', '%', &
         standard_name = 'relative_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_rhs > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_rhs, 'coordinates', 'height2m' )

    id_ts = &
         fms_diag_register_diag_field ( mod_name, 'ts', atmos_axes, Time, &
         'Surface Temperature', 'K', &
         standard_name = 'surface_temperature', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=trange )

    id_psl = &
         fms_diag_register_diag_field ( mod_name, 'psl', atmos_axes, Time, &
         'Sea Level Pressure', 'Pa', &
         standard_name = 'air_pressure_at_sea_level', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_tauu = &
         fms_diag_register_diag_field ( mod_name, 'tauu', atmos_axes, Time, &
         'Surface Downward Eastward Wind Stress', 'Pa', &
         standard_name = 'surface_downward_eastward_stress', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_tauv = &
         fms_diag_register_diag_field ( mod_name, 'tauv', atmos_axes, Time, &
         'Surface Downward Northward Wind Stress', 'Pa', &
         standard_name = 'surface_downward_northward_stress', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_hfss = &
         fms_diag_register_diag_field ( mod_name, 'hfss', atmos_axes, Time, &
         'Surface Upward Sensible Heat Flux', 'W m-2', &
         standard_name = 'surface_upward_sensible_heat_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_hfls = &
         fms_diag_register_diag_field ( mod_name, 'hfls', atmos_axes, Time, &
         'Surface Upward Latent Heat Flux', 'W m-2', &
         standard_name = 'surface_upward_latent_heat_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_hfls > 0 ) call fms_diag_field_add_attribute( id_hfls, 'comment', 'Lv*evap' )

    id_evspsbl = &
         fms_diag_register_diag_field( mod_name, 'evspsbl', atmos_axes, Time, &
         'Evaporation', 'kg m-2 s-1', &
         standard_name = 'water_evaporation_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_sftlf = &
         fms_diag_register_static_field ( mod_name, 'sftlf', atmos_axes,  &
         'Fraction of the Grid Cell Occupied by Land', '1.0', &
         standard_name = 'land_area_fraction', area=area_id, &
         interp_method = "conserve_order1" )

    id_tslsi = &
         fms_diag_register_diag_field ( mod_name, 'tslsi', atmos_axes, Time,  &
         'Surface Temperature Where Land or Sea Ice', 'K', &
         standard_name = 'surface_temperature', area=area_id, &
         mask_variant=.true., missing_value=CMOR_MISSING_VALUE )

    id_tos = &
         fms_diag_register_diag_field ( mod_name, 'tos', atmos_axes, Time,  &
         'Sea Surface Temperature', 'K', &
         standard_name = 'sea_surface_temperature', area=area_id, &
         mask_variant=.true., missing_value=CMOR_MISSING_VALUE )

    id_sic = &
         fms_diag_register_diag_field ( mod_name, 'sic', atmos_axes, Time,  &
         'Sea Ice Area Fraction', '1.0', &
         standard_name = 'sea_ice_area_fraction', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_sic > 0 ) call fms_diag_field_add_attribute( id_sic, 'comment', &
         'averaged over the ocean portion of grid box' )
#else
    id_tas = register_cmip_diag_field_2d ( mod_name, 'tas', Time, &
                            'Near-Surface Air Temperature', 'K' , &
                             standard_name='air_temperature' )
    if ( id_tas > 0 .and. id_height2m > 0) &
       call fms_diag_field_add_attribute( id_tas, 'coordinates', 'height2m' )

    id_uas = register_cmip_diag_field_2d ( mod_name, 'uas', Time, &
                           'Eastward Near-Surface Wind', 'm s-1', &
                            standard_name='eastward_wind' )
    if ( id_uas > 0 .and. id_height10m > 0) &
       call fms_diag_field_add_attribute( id_uas, 'coordinates', 'height10m' )

    id_vas = register_cmip_diag_field_2d ( mod_name, 'vas', Time, &
                          'Northward Near-Surface Wind', 'm s-1', &
                           standard_name='northward_wind' )
    if ( id_vas > 0 .and. id_height10m > 0 ) &
       call fms_diag_field_add_attribute( id_vas, 'coordinates', 'height10m' )

    id_sfcWind = register_cmip_diag_field_2d ( mod_name, 'sfcWind', Time, &
                                      'Near-Surface Wind Speed', 'm s-1', &
                                       standard_name='wind_speed' )
    if ( id_sfcWind > 0 .and. id_height10m > 0 ) &
       call fms_diag_field_add_attribute( id_sfcWind, 'coordinates', 'height10m' )

    id_huss = register_cmip_diag_field_2d ( mod_name, 'huss', Time, &
                           'Near-Surface Specific Humidity', '1.0', &
                            standard_name='specific_humidity' )
    if ( id_huss > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_huss, 'coordinates', 'height2m' )

    id_hurs = register_cmip_diag_field_2d ( mod_name, 'hurs', Time, &
                             'Near-Surface Relative Humidity', '%', &
                              standard_name='relative_humidity' )
    if ( id_hurs > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_hurs, 'coordinates', 'height2m' )

    id_rhs = register_cmip_diag_field_2d ( mod_name, 'rhs', Time, &
                           'Near-Surface Relative Humidity', '%', &
                            standard_name='relative_humidity' )
    if ( id_rhs > 0 .and. id_height2m > 0 ) &
       call fms_diag_field_add_attribute( id_rhs, 'coordinates', 'height2m' )

    id_ts = register_cmip_diag_field_2d ( mod_name, 'ts', Time, &
                                    'Surface Temperature', 'K', &
                            standard_name='surface_temperature' )

    id_psl = register_cmip_diag_field_2d ( mod_name, 'psl', Time, &
                                      'Sea Level Pressure', 'Pa', &
                        standard_name='air_pressure_at_sea_level' )

    id_tauu = register_cmip_diag_field_2d ( mod_name, 'tauu', Time, &
                     'Surface Downward Eastward Wind Stress', 'Pa', &
                   standard_name='surface_downward_eastward_stress' )

    id_tauv = register_cmip_diag_field_2d ( mod_name, 'tauv', Time, &
                    'Surface Downward Northward Wind Stress', 'Pa', &
                  standard_name='surface_downward_northward_stress' )

    id_hfss = register_cmip_diag_field_2d ( mod_name, 'hfss', Time, &
                      'Surface Upward Sensible Heat Flux', 'W m-2', &
                  standard_name='surface_upward_sensible_heat_flux' )

    id_hfls = register_cmip_diag_field_2d ( mod_name, 'hfls', Time, &
                        'Surface Upward Latent Heat Flux', 'W m-2', &
                    standard_name='surface_upward_latent_heat_flux' )
    if ( id_hfls > 0 ) call fms_diag_field_add_attribute( id_hfls, 'comment', 'Lv*evap' )

    id_evspsbl = register_cmip_diag_field_2d ( mod_name, 'evspsbl', Time, &
                                             'Evaporation', 'kg m-2 s-1', &
                                   standard_name='water_evaporation_flux' )

    id_sftlf = fms_diag_register_static_field ( mod_name, 'sftlf', atmos_axes,  &
                  'Fraction of the Grid Cell Occupied by Land', '1.0', &
                     standard_name='land_area_fraction', area=area_id, &
                     interp_method='conserve_order1' )

    id_tslsi = register_cmip_diag_field_2d ( mod_name, 'tslsi', Time,  &
                     'Surface Temperature Where Land or Sea Ice', 'K', &
                                  standard_name='surface_temperature', &
                                     mask_variant=.true. )

    ! tos,sic are ocean,seaIce fields on the atmos grid
    ! useful for amip-type runs

    id_tos = register_cmip_diag_field_2d ( mod_name, 'tos', Time,  &
                                   'Sea Surface Temperature', 'K', &
                          standard_name='sea_surface_temperature', &
                          mask_variant=.true. )

    id_sic = register_cmip_diag_field_2d ( mod_name, 'sic', Time,  &
                                   'Sea Ice Area Fraction', '1.0', &
                            standard_name='sea_ice_area_fraction' )
    if ( id_sic > 0 ) call fms_diag_field_add_attribute( id_sic, 'comment', &
         'averaged over the ocean portion of grid box' )

    !----- initialize global integrals for netCDF output -----
    id_evspsbl_g = register_global_diag_field ( 'evspsbl', Time, &
                                    'Evaporation', 'mm d-1', &
                          standard_name='water_evaporation_flux' )

    id_ts_g = register_global_diag_field ( 'ts', Time, &
                                     'Surface Temperature', 'K', &
                             standard_name='surface_temperature' )

    id_tas_g = register_global_diag_field ( 'tas', Time, &
                           'Near-Surface Air Temperature', 'K' , &
                                 standard_name='air_temperature' )
    if ( id_tas_g > 0 .and. id_height2m > 0) &
         call fms_diag_field_add_attribute ( get_global_diag_field_id(id_tas_g), 'coordinates', 'height2m' )

    id_tasl_g = register_global_diag_field ( 'tasl', Time, &
                           'Near-Surface Air Temperature (Land Only)', 'K' , &
                                 standard_name='air_temperature' )
#if defined(_USE_LEGACY_LAND_) || defined(use_AM3_physics)
    if(id_tasl_g>0) then
       call fms_mpp_error(WARNING, "diag_field_init: field tasl is registered, but macro "// &
             "_USE_LEGACY_LAND_ or use_AM3_physics is defined, no data will be written out")
    endif
#endif
    if ( id_tasl_g > 0 .and. id_height2m > 0) &
         call fms_diag_field_add_attribute ( get_global_diag_field_id(id_tasl_g), 'coordinates', 'height2m' )

    id_hfss_g = register_global_diag_field ( 'hfss', Time, &
                   'Surface Upward Sensible Heat Flux', 'W m-2', &
               standard_name='surface_upward_sensible_heat_flux' )

    id_hfls_g = register_global_diag_field ( 'hfls', Time, &
                  'Surface Upward Latent Heat Flux', 'W m-2', &
                  standard_name='surface_upward_latent_heat_flux')
    if ( id_hfls_g > 0 ) &
         call fms_diag_field_add_attribute( get_global_diag_field_id(id_hfls_g), 'comment', 'Lv*evap' )

    id_rls_g = register_global_diag_field ( 'rls', Time, &
                   'Net Longwave Surface Radiation', 'W m-2', &
               standard_name='surface_net_longwave_flux' )

#endif
    !}

  end subroutine diag_field_init


  !> Subroutine divide_by_area divides data on a grid by the grid cell area only for
  !! cells with non-zero area.  This subroutine iscurrently not used.  Note, a
  !! similar subroutine also exists in ice_ocean_flux_exchange_mod
  subroutine divide_by_area(data, area)

    !> data to be divided
    real, intent(inout) :: data(:,:)
    !> area used as denominator
    real, intent(in) :: area(:,:)

    !> CHECK TO ENSURE SHAPE OF DATA IS THE SAME AS SHAPE OF AREA
    !! IF SHAPES MISMATCH, RETURN 
    !{
    if(size(data, dim=1) /= size(area, dim=1) .or. size(data, dim=2) /= size(area, dim=2)) then
       ! no op
       return
    endif
    !}
    
    !> DIVIDE DATA BY GRID CELL AREA WHERE AREA /= 0.0
    !{
    where(area /= 0.0)
       data = data / area
    end where
    !}

  end subroutine divide_by_area

  
  !> Subroutine send_ice_mask_sic sends the ice mask to diag_manager.
  !! If the variables ice_mask or sic have been registered with diag_manager,
  !! this subroutine  maps the fractional amount of sea ice
  !! from the OCN grid to the ATM grid and sends the data to the diag_manager buffer.
  ! This was called inside flux_ocean_to_ice.  Why?
  subroutine send_ice_mask_sic(Time)

    !> current model time
    type(FmsTime_type),  intent(in) :: Time

    real, dimension(nxc_ice, nyc_ice, nk_ice) :: ice_frac
    real, dimension(n_xgrid_sfc)               :: ex_ice_frac
    real, dimension(ni_atm, nj_atm)           :: diag_atm, ocean_frac
    logical :: used

    !> IF ID_ICE_MASK > 0 OR ID_SIC > 0    
    if ( id_ice_mask > 0 .or. id_sic > 0) then

       !> INITIALIZE ICE_FRAC
       !{
       ice_frac        = 1.
       ice_frac(:,:,1) = 0.
       ex_ice_frac     = 0.
       !}
       
       !> MAP ICE_MASK FROM THE OCN GRID TO THE EXCHANGE GRID
       call fms_xgrid_put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)

       !> MAP ICE_MASK FROM THE EXCHANGE GRID TO THE ATM GRID
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)

       !> IF ID_ICE_MASK > 0, SEND ICE_MASK TO THE DIAG_MANAGER BUFFER
       if ( id_ice_mask > 0 ) used = fms_diag_send_data ( id_ice_mask, diag_atm, Time )

       !> FOR CMIP, IF ID_SIC > 0, COMPUTE SEA ICE FRACTIONAL AREA FOR ATM GRID CELLS THAT ARE OVER THE OCEAN
       !! AND NORMALIZE AREA BY THE FRACTION OF ATMOS GRID CELL THAT IS OCEAN
       !{
       if ( id_sic > 0) then
          ice_frac = 1.
          ex_ice_frac = 0.
          call fms_xgrid_put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
          call fms_xgrid_get_from_xgrid (ocean_frac, 'ATM', ex_ice_frac, xmap_sfc)
          where (ocean_frac > 0.0)
             diag_atm = min(1., diag_atm/ocean_frac) ! CMIP6 as fraction
             ocean_frac = 1.0
          elsewhere
             diag_atm = 0.0
             ocean_frac = 0.0
          endwhere
          used = fms_diag_send_data ( id_sic, diag_atm, Time, rmask=ocean_frac )
       endif
    endif
    !}

  end subroutine send_ice_mask_sic

  
  !> Subroutine atm_stock_integrate integrates over the total precipitation
  !! (liquid and frozen) in the atmosphere and multiply the integrated value by
  !! the timestep dt.  This subroutine is called in flux_exchange_mod/flux_check_stocks
  subroutine atm_stock_integrate(Atm, res)

    !> derived type holding the atmosphere boundary data
    type(atmos_data_type), intent(in) :: Atm
    !> integrated value
    real, intent(out) :: res

    integer :: ier

    !> CALL FMS_XGRID_STOCK_INTEGRATE    
    call fms_xgrid_stock_integrate_2d(Atm%lprec + Atm%fprec, xmap=xmap_sfc, delta_t=Dt_atm, &
         & radius=Radius, res=res, ier=ier)

  end subroutine atm_stock_integrate

end module atm_land_ice_flux_exchange_mod
