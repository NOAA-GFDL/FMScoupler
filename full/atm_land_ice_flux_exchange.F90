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

module atm_land_ice_flux_exchange_mod

!! Components
  use ocean_model_mod,    only: ocean_model_init_sfc, ocean_model_flux_init, ocean_model_data_get
  use   atmos_model_mod,  only: atmos_data_type, land_ice_atmos_boundary_type
  use   ocean_model_mod,  only: ocean_public_type, ice_ocean_boundary_type
  use   ocean_model_mod,  only: ocean_state_type
  use   ice_model_mod,    only: ice_data_type, land_ice_boundary_type, ocean_ice_boundary_type
  use   ice_model_mod,    only: atmos_ice_boundary_type, Ice_stock_pe
  use   ice_model_mod,    only: update_ice_atm_deposition_flux
  use    land_model_mod,  only: land_data_type, atmos_land_boundary_type
  use  surface_flux_mod,  only: surface_flux, surface_flux_init
  use land_model_mod,          only: Lnd_stock_pe
  use ocean_model_mod,         only: Ocean_stock_pe
  use atmos_model_mod,         only: Atm_stock_pe
  use atmos_ocean_fluxes_mod,  only: atmos_ocean_fluxes_init
  use atmos_ocean_fluxes_calc_mod, only: atmos_ocean_fluxes_calc
  use atmos_ocean_dep_fluxes_calc_mod, only: atmos_ocean_dep_fluxes_calc

!! Conditional Imports
#ifndef _USE_LEGACY_LAND_
  use    land_model_mod,  only: set_default_diag_filter, register_tiled_diag_field
  use    land_model_mod,  only: send_tile_data, dump_tile_diag_fields
#endif

#ifdef use_AM3_physics
  use atmos_tracer_driver_mod, only: atmos_tracer_flux_init
#else
  use atmos_tracer_driver_mod, only: atmos_tracer_flux_init, &
       atmos_tracer_has_surf_setl_flux, get_atmos_tracer_surf_setl_flux
  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data_down
  use atmos_cmip_diag_mod,   only: register_cmip_diag_field_2d
  use atmos_global_diag_mod, only: register_global_diag_field, &
                                   get_global_diag_field_id, &
                                   send_global_diag
#ifndef _USE_LEGACY_LAND_
  use land_model_mod,        only: send_global_land_diag
#endif
#endif


#ifdef SCM
  ! option to override various surface boundary conditions for SCM
  use scm_forc_mod,            only: do_specified_flux, scm_surface_flux,             &
                                     do_specified_tskin, TSKIN,                       &
                                     do_specified_albedo, ALBEDO_OBS,                 &
                                     do_specified_rough_leng, ROUGH_MOM, ROUGH_HEAT,  &
                                     do_specified_land
#endif

!! FMS
use FMS
use FMSconstants, only: rdgas, rvgas, cp_air, stefan, WTMAIR, HLV, HLF, Radius, &
                        PI, CP_OCEAN, WTMCO2, WTMC, EPSLN, GRAV

  implicit none
  include 'netcdf.inc'
  private

  public :: atm_land_ice_flux_exchange_init,   &
            sfc_boundary_layer,   &
            generate_sfc_xgrid,   &
            flux_down_from_atmos, &
            flux_up_to_atmos,     &
            flux_atmos_to_ocean,  &
            flux_ex_arrays_dealloc,&
            atm_stock_integrate,  &
            send_ice_mask_sic

  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id$'
  character(len=128) :: tag = '$Name$'
  !-----------------------------------------------------------------------
  !---- exchange grid maps -----

  type(xmap_type), save :: xmap_sfc

  integer         :: n_xgrid_sfc=0

  !-----------------------------------------------------------------------
  !-------- namelist (for diagnostics) ------

  character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
             id_rough_moist, id_rough_heat, id_rough_mom,    &
             id_land_mask,   id_ice_mask,     &
             id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
             id_t_surf, id_t_flux, id_r_flux, id_q_flux, id_slp,      &
             id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
             id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref, id_wind_ref,  &
             id_del_h,  id_del_m,  id_del_q,  id_rough_scale,         &
             id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm, id_gust, &
             id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land, &
             id_q_ref,  id_q_ref_land, id_q_flux_land, id_rh_ref_cmip, &
             id_hussLut_land, id_tasLut_land, id_t_flux_land
  integer :: id_co2_atm_dvmr, id_co2_surf_dvmr
! 2017/08/15 jgj added
  integer :: id_co2_bot, id_co2_flux_pcair_atm, id_o2_flux_pcair_atm

  integer, allocatable :: id_tr_atm(:), id_tr_surf(:), id_tr_flux(:), id_tr_mol_flux(:)
  integer, allocatable :: id_tr_mol_flux0(:) !f1p
  integer, allocatable :: id_tr_flux_land(:), id_tr_mol_flux_land(:)

  ! id's for cmip specific fields
  integer :: id_tas, id_uas, id_vas, id_ts, id_psl, &
             id_sfcWind, id_tauu, id_tauv, &
             id_hurs, id_huss, id_evspsbl, id_hfls, id_hfss, &
             id_rhs, id_sftlf, id_tos, id_sic, id_tslsi, &
             id_height2m, id_height10m

  ! globally averaged diagnostics
  integer :: id_evspsbl_g, id_ts_g, id_tas_g, id_tasl_g, id_hfss_g, id_hfls_g, id_rls_g

  logical :: first_static = .true.
  logical :: do_init = .true.
  integer :: remap_method = 1

  real, parameter :: bound_tol = 1e-7

  real, parameter :: d622 = rdgas/rvgas
  real, parameter :: d378 = 1.0-d622
  real, parameter :: tfreeze = 273.15
  real, allocatable, dimension(:,:) :: frac_precip

  !--- the following is from flux_exchange_nml
  real    :: z_ref_heat =  2. !< Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
  real    :: z_ref_mom  = 10. !< Reference height (meters) for mementum diagnostics (u_ref, v_ref, del_m)
  logical :: ex_u_star_smooth_bug = .false. !< By default, the global exchange grid \c u_star will not be interpolated
                                            !! from atmospheric grid, this is different from Jakarta behavior and will
                                            !! change answers.  So to preserve Jakarta behavior and reproduce answers
                                            !! explicitly set this namelist variable to .true. in input.nml.
  logical :: sw1way_bug = .false.
  logical :: do_area_weighted_flux = .FALSE.
  logical :: do_forecast = .false.
  integer :: nblocks = 1
  logical :: partition_fprec_from_lprec = .FALSE. !< option for ATM override experiments where liquid+frozen precip are combined
                                                  !! This option will convert liquid precip to snow when t_ref is less than
                                                  !! tfreeze parameter
  logical :: scale_precip_2d = .false.

  integer              :: my_nblocks = 1
  integer, allocatable :: block_start(:), block_end(:)

  ! ---- allocatable module storage --------------------------------------------
  real, allocatable, dimension(:) :: &
                                ! NOTE: T canopy is only differet from t_surf over vegetated land
       ex_t_surf,    &   !< surface temperature for radiation calc, degK
       ex_t_surf_miz,&   !< miz
       ex_t_ca,      &   !< near-surface (canopy) air temperature, degK
       ex_p_surf,    &   !< surface pressure
       ex_slp,       &   !< surface pressure

       ex_flux_t,    &   !< sens heat flux
       ex_flux_lw,   &   !< longwave radiation flux

       ex_dhdt_surf, &   !< d(sens.heat.flux)/d(T canopy)
       ex_dedt_surf, &   !< d(water.vap.flux)/d(T canopy)
       ex_dqsatdt_surf, &   !< d(water.vap.flux)/d(q canopy)
       ex_e_q_n,     &
       ex_drdt_surf, &   !< d(LW flux)/d(T surf)
       ex_dhdt_atm,  &   !< d(sens.heat.flux)/d(T atm)
       ex_flux_u,    &   !< u stress on atmosphere
       ex_flux_v,    &   !< v stress on atmosphere
       ex_dtaudu_atm,&   !< d(stress)/d(u)
       ex_dtaudv_atm,&   !< d(stress)/d(v)
       ex_seawater,  &
       ex_albedo_fix,&
       ex_albedo_vis_dir_fix,&
       ex_albedo_nir_dir_fix,&
       ex_albedo_vis_dif_fix,&
       ex_albedo_nir_dif_fix,&
       ex_old_albedo,&   !< old value of albedo for downward flux calculations
       ex_drag_q,    &   !< q drag.coeff.
       ex_cd_t,      &
       ex_cd_m,      &
       ex_b_star,    &
       ex_u_star,    &
       ex_wind,      &
       ex_z_atm

#ifdef SCM
  real, allocatable, dimension(:) :: &
       ex_dhdt_surf_forland, &
       ex_dedt_surf_forland, &
       ex_dedq_surf_forland
#endif

  real, allocatable, dimension(:,:) :: &
       ex_tr_surf,    & !< near-surface tracer fields
       ex_flux_tr,    & !< tracer fluxes
       ex_dfdtr_surf, & !< d(tracer flux)/d(surf tracer)
       ex_dfdtr_atm,  & !< d(tracer flux)/d(atm tracer)
       ex_e_tr_n,     & !< coefficient in implicit scheme
       ex_f_tr_delt_n   !< coefficient in implicit scheme

  logical, allocatable, dimension(:) :: &
       ex_avail,     &   !< true where data on exchange grid are available
       ex_land           !< true if exchange grid cell is over land
  real, allocatable, dimension(:) :: &
       ex_e_t_n,      &
       ex_f_t_delt_n

  integer :: n_atm_tr  !< number of prognostic tracers in the atmos model
  integer :: n_atm_tr_tot  !< number of prognostic tracers in the atmos model
  integer :: n_lnd_tr  !< number of prognostic tracers in the land model
  integer :: n_lnd_tr_tot  !< number of prognostic tracers in the land model
  integer :: n_exch_tr !< number of tracers exchanged between models

  type :: tracer_ind_type
     integer :: atm, ice, lnd !< indices of the tracer in the respective models
  end type tracer_ind_type
  type(tracer_ind_type), allocatable :: tr_table(:) !< table of tracer indices
  type :: tracer_exch_ind_type
     integer :: exch = 0  !< exchange grid index
     integer :: ice = 0   !< ice model index
     integer :: lnd = 0   !< land model index
  end type tracer_exch_ind_type
  type(tracer_exch_ind_type), allocatable :: tr_table_map(:) !< map atm tracers to exchange, ice and land variables
  integer :: isphum = NO_TRACER       !< index of specific humidity tracer in tracer table
  integer :: ico2   = NO_TRACER       !< index of co2 tracer in tracer table
  integer :: inh3   = NO_TRACER       !< index of nh3 tracer in tracer table
  type(coupler_1d_bc_type), pointer :: ex_gas_fields_atm=>NULL() !< gas fields in atm
                                                                 !< Place holder for various atmospheric fields.
  type(coupler_1d_bc_type), pointer :: ex_gas_fields_ice=>NULL() ! gas fields on ice
  type(coupler_1d_bc_type), pointer :: ex_gas_fluxes=>NULL()     ! gas flux
                                                                 !< Place holder of intermediate calculations, such as
                                                                 !< piston velocities etc.

  interface put_logical_to_real
     module procedure put_logical_to_real_sg
     module procedure put_logical_to_real_ug
  end interface

  integer :: ni_atm, nj_atm !< to do atmos diagnostic from flux_ocean_to_ice
  real, dimension(3) :: ccc !< for conservation checks
  !Balaji, sets boundary_type%xtype
  !  REGRID: grids are physically different, pass via exchange grid
  !  REDIST: same physical grid, different decomposition, must move data around
  !  DIRECT: same physical grid, same domain decomposition, can directly copy data
  integer, parameter :: REGRID=1, REDIST=2, DIRECT=3
  integer :: cplClock, sfcClock, fluxAtmDnClock, regenClock, fluxAtmUpClock

  ! Exchange grid indices
  integer :: X1_GRID_ATM, X1_GRID_ICE, X1_GRID_LND
  real    :: Dt_atm, Dt_cpl
  integer :: nxc_ice=0, nyc_ice=0, nk_ice=0
  integer :: nxc_lnd=0, nyc_lnd=0

contains

  !#######################################################################
  !> \brief Initialization routine.
  !!
  !! Initializes the interpolation routines,diagnostics and boundary data
  !!
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere resolution"
  !!    The atmosphere grid size from file grid_spec.nc is not compatible with the atmosphere
  !!    resolution from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)"
  !!    The longitude from file grid_spec.nc ( from field yba ) is different from the longitude from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)"
  !!    The longitude from file grid_spec.nc ( from field xba ) is different from the longitude from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere latitudes (see grid_spec.nc)"
  !!    The latitude from file grid_spec.nc is different from the latitude from atmosphere model.
  subroutine atm_land_ice_flux_exchange_init(Time, Atm, Land, Ice, atmos_ice_boundary, land_ice_atmos_boundary, &
                                             Dt_atm_in, Dt_cpl_in, z_ref_heat_in, z_ref_mom_in,                 &
                                             ex_u_star_smooth_bug_in, sw1way_bug_in, do_area_weighted_flux_in,  &
                                             do_forecast_in, partition_fprec_from_lprec_in, scale_precip_2d_in, &
                                             nblocks_in, cplClock_in, ex_gas_fields_atm_in, &
                                             ex_gas_fields_ice_in, ex_gas_fluxes_in)
    type(time_type),                   intent(in)    :: Time !< The model's current time
    type(atmos_data_type),             intent(inout) :: Atm  !< A derived data type to specify atmosphere boundary data
    type(land_data_type),              intent(in)    :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),               intent(inout) :: Ice  !< A derived data type to specify ice boundary data
    type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary !< A derived data type to specify properties and fluxes passed from atmosphere to ice
    type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary !< A derived data type to specify properties and fluxes passed from exchange grid to
    !! the atmosphere, land and ice
    real,                 intent(in)    :: Dt_atm_in !< Atmosphere time step in seconds
    real,                 intent(in)    :: Dt_cpl_in !< Coupled time step in seconds
    real,                 intent(in)    :: z_ref_heat_in, z_ref_mom_in
    logical,              intent(in)    :: ex_u_star_smooth_bug_in, scale_precip_2d_in
    logical,              intent(in)    :: sw1way_bug_in, do_area_weighted_flux_in
    logical,              intent(in)    :: do_forecast_in, partition_fprec_from_lprec_in
    integer,              intent(in)    :: nblocks_in
    integer,              intent(in)    :: cplClock_in
    type(coupler_1d_bc_type), intent(in), target :: ex_gas_fields_atm_in, ex_gas_fields_ice_in, ex_gas_fluxes_in

    character(len=48), parameter :: module_name = 'flux_exchange_mod'
    character(len=64), parameter    :: sub_name = 'flux_exchange_init'
    character(len=256), parameter   :: note_header = '==>Note from ' // trim(module_name) //     &
         '(' // trim(sub_name) // '):'
    integer        :: i, n
    integer        :: outunit, logunit
    integer :: is, ie, js, je, kd
    character(32) :: tr_name
    logical       :: found
    character(32)  :: method
    character(512) :: parameters
    real           :: value

    Dt_atm = Dt_atm_in
    Dt_cpl = Dt_cpl_in
    z_ref_heat = z_ref_heat_in
    z_ref_mom = z_ref_mom_in
    ex_u_star_smooth_bug = ex_u_star_smooth_bug_in
    sw1way_bug = sw1way_bug_in
    do_area_weighted_flux = do_area_weighted_flux_in
    do_forecast = do_forecast_in
    partition_fprec_from_lprec = partition_fprec_from_lprec_in
    scale_precip_2d = scale_precip_2d_in
    nblocks = nblocks_in
    cplClock = cplClock_in
    ex_gas_fields_atm => ex_gas_fields_atm_in
    ex_gas_fields_ice => ex_gas_fields_ice_in
    ex_gas_fluxes     => ex_gas_fluxes_in

    outunit = stdout(); logunit = stdlog()

    allocate(block_start(nblocks), block_end(nblocks))

    !----- find out number of atmospheric prognostic tracers and index of specific
    !      humidity in the tracer table
    call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
         num_prog=n_atm_tr)
    call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
         num_prog=n_lnd_tr)

    ! assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    allocate(tr_table_map(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
       call get_tracer_names( MODEL_ATMOS, i, tr_name )
       tr_table(n)%atm = i
       tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
       tr_table_map(i)%ice = tr_table(n)%ice
       tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
       tr_table_map(i)%lnd = tr_table(n)%lnd
       if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) then
          tr_table_map(i)%exch = n
          n = n + 1
       endif
    enddo
    n_exch_tr = n - 1
    !
    !     Set up tracer table entries for ocean-atm gas fluxes where the names of tracers in the
    !     atmosphere and ocean may not be equal
    !
    do n = 1, ex_gas_fluxes%num_bcs  !{
       if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
          found = .false.
          do i = 1, n_exch_tr  !{
             if (ex_gas_fluxes%bc(n)%atm_tr_index .eq. tr_table(i)%atm) then
                found = .true.
                exit
             endif
          enddo  !} i
          if (.not. found) then
             n_exch_tr = n_exch_tr + 1
             tr_table(n_exch_tr)%atm = ex_gas_fluxes%bc(n)%atm_tr_index
             tr_table(n_exch_tr)%ice = NO_TRACER ! because ocean-atm gas fluxes are not held in the ice model as tracers
             tr_table(n_exch_tr)%lnd = NO_TRACER ! because this would have been found above
             tr_table_map(n_exch_tr)%exch = n_exch_tr
             tr_table_map(n_exch_tr)%ice = tr_table(n_exch_tr)%ice
             tr_table_map(n_exch_tr)%lnd = tr_table(n_exch_tr)%lnd
          endif
       endif  !}
    enddo  !} n
    write(outunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
    write(logunit,*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
    do i = 1,n_exch_tr
       call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
       write(outunit,*)'Tracer field name :'//trim(tr_name)
       write(logunit,*)'Tracer field name :'//trim(tr_name)
    enddo

    ! find out which tracer is specific humidity

    ! +fix-me-slm+ specific humidity may not be present if we are running with
    ! dry atmosphere. Besides, model may use mixing ratio ('mix_rat') (?). However,
    ! some atmos code also assumes 'sphum' is present, so for now the following
    ! code may be good enough.

    do i = 1,n_exch_tr
       call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
       if(lowercase(tr_name)=='sphum') then
          isphum = i
       endif
       ! jgj: find out which exchange tracer is co2
       if(lowercase(tr_name)=='co2') then
          ico2 = i
          write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',ico2
       endif
       if(lowercase(tr_name)=='nh3') then
          inh3 = i
          write(outunit,*)'Exchange tracer index for '//trim(tr_name),' : ',inh3
       endif
    enddo

    if (isphum==NO_TRACER) then
       call error_mesg('atm_land_ice_flux_exchange_mod',&
            'tracer "sphum" must be present in the atmosphere', FATAL )
    endif

    if (ico2==NO_TRACER) then
       call error_mesg('atm_land_ice_flux_exchange_mod',&
            'tracer "co2" not present in the atmosphere', NOTE )
    endif

    !--------- read gridspec file ------------------
    !only atmos pelists needs to do it here, ocean model will do it elsewhere


    !
    ! check atmosphere and grid_spec.nc have same atmosphere lat/lon boundaries
    !
    call mpp_get_compute_domain(Atm%domain, is, ie, js, je)

    if (scale_precip_2d) then
       allocate(frac_precip(is:ie,js:je))
       frac_precip=0.0
    endif

    call xgrid_init(remap_method)
#ifndef _USE_LEGACY_LAND_
    call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
         (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
         "INPUT/grid_spec.nc", Atm%grid, lnd_ug_domain=Land%ug_domain)
#else
    call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
         (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
         "INPUT/grid_spec.nc", Atm%grid)
#endif
    ! exchange grid indices
    X1_GRID_ATM = 1; X1_GRID_ICE = 2; X1_GRID_LND = 3;
    call generate_sfc_xgrid( Land, Ice )
    if (n_xgrid_sfc.eq.1) write (*,'(a,i6,6x,a)') 'PE = ', mpp_pe(), 'Surface exchange size equals one.'

    call surface_flux_init()

    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    !----- initialize quantities for global integral package -----

    !! call diag_integral_field_init ('prec', 'f6.3')
    call diag_integral_field_init ('evap', 'f6.3')
#ifndef use_AM3_physics
    call diag_integral_field_init ('t_surf', 'f10.3') !miz
    call diag_integral_field_init ('t_ref',  'f10.3') !miz
#endif

    !-----------------------------------------------------------------------
    !----- initialize diagnostic fields -----
    !----- all fields will be output on the atmospheric grid -----

    call diag_field_init ( Time, Atm%axes(1:2), Land%axes, Land%pe )
    ni_atm = size(Atm%lon_bnd,1)-1 ! to dimension "diag_atm"
    nj_atm = size(Atm%lon_bnd,2)-1 ! in flux_ocean_to_ice

    !Balaji

    !allocate atmos_ice_boundary
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
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
    ! initialize boundary values for override experiments (mjh)
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

    !         allocate fields for extra fluxes
    ! Copying initialized gas fluxes from exchange grid to atmosphere_ice boundary

    call coupler_type_copy(ex_gas_fluxes, atmos_ice_boundary%fluxes, is, ie, js, je, kd,    &
         mod_name, Ice%axes, Time, suffix = '_atm_ice')

    !--- Ice%ocean_fields and Ice%ocean_fluxes_top will not be passed to ocean, so these two
    !--- coupler_type_copy calls are moved from ice_ocean_flux_init to here.
    if (.not.coupler_type_initialized(Ice%ocean_fields)) &
      call coupler_type_spawn(ex_gas_fields_ice, Ice%ocean_fields, (/is,is,ie,ie/), &
                              (/js,js,je,je/), (/1, kd/), suffix = '_ice')
    call coupler_type_set_diags(Ice%ocean_fields, 'ice_flux', Ice%axes, Time)

    !allocate land_ice_atmos_boundary
    call mpp_get_compute_domain( Atm%domain, is, ie, js, je )
    allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%u_ref(is:ie,js:je) )  ! bqx
    allocate( land_ice_atmos_boundary%v_ref(is:ie,js:je) )  ! bqx
    allocate( land_ice_atmos_boundary%t_ref(is:ie,js:je) )  ! cjg: PBL depth mods
    allocate( land_ice_atmos_boundary%q_ref(is:ie,js:je) )  ! cjg: PBL depth mods
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
    allocate( land_ice_atmos_boundary%shflx(is:ie,js:je) )!miz
    allocate( land_ice_atmos_boundary%lhflx(is:ie,js:je) )!miz
#endif
    allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
    allocate( land_ice_atmos_boundary%frac_open_sea(is:ie,js:je) )
    ! initialize boundary values for override experiments (mjh)
    land_ice_atmos_boundary%t=273.0
    land_ice_atmos_boundary%u_ref=0.0   ! bqx
    land_ice_atmos_boundary%v_ref=0.0   ! bqx
    land_ice_atmos_boundary%t_ref=273.0   ! cjg: PBL depth mods
    land_ice_atmos_boundary%q_ref=0.0     ! cjg: PBL depth mods
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
#ifndef use_AM3_physics
    land_ice_atmos_boundary%shflx=0.0
    land_ice_atmos_boundary%lhflx=0.0
#endif
    land_ice_atmos_boundary%rough_mom=0.01
    land_ice_atmos_boundary%frac_open_sea=0.0

    ! allocate fields for extra tracers
    ! The first call is no longer necessary, the fluxes will be passed by the land module
    ! The 2nd call is useful in the case of a ocean model only simulation
    !
    call coupler_type_copy(ex_gas_fields_atm, Atm%fields, is, ie, js, je,                   &
         mod_name, Atm%axes(1:2), Time, suffix = '_atm')

    if( Ice%pe) then
       call mpp_get_compute_domain(Ice%domain, xsize=nxc_ice, ysize=nyc_ice)
       nk_ice = size(Ice%part_size,3)
    endif

    if( Land%pe) then
       call mpp_get_compute_domain(Land%domain, xsize=nxc_lnd, ysize=nyc_lnd)
    endif

    !Balaji: clocks on atm%pe only
    sfcClock = mpp_clock_id( 'SFC boundary layer', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    fluxAtmDnClock = mpp_clock_id( 'Flux DN from atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    regenClock = mpp_clock_id( 'XGrid generation', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxAtmUpClock = mpp_clock_id( 'Flux UP to atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )

    do_init = .false.

  end subroutine atm_land_ice_flux_exchange_init

  !#######################################################################
  !> \brief Computes explicit fluxes as well as derivatives that will be used to compute an implicit flux correction.
  !!
  !!
  !!  The following quantities in the land_ice_atmos_boundary_type are computed:
  !!
  !! <pre>
  !!         t_surf_atm = surface temperature (used for radiation)    (K)
  !!         albedo_atm = surface albedo      (used for radiation)    (nondimensional)
  !!      rough_mom_atm = surface roughness for momentum (m)
  !!      land_frac_atm = fractional area of land beneath an atmospheric
  !!                      grid box
  !!         dtaudu_atm, dtaudv_atm = derivatives of wind stress w.r.t. the
  !!                                  lowest level wind speed  (Pa/(m/s))
  !!         flux_u_atm = zonal wind stress  (Pa)
  !!         flux_v_atm = meridional wind stress (Pa)
  !!         u_star_atm = friction velocity (m/s)
  !!         b_star_atm = buoyancy scale    (m2/s)
  !! </pre>
  !! \note `u_star` and `b_star` are defined so that `u_star**2` is the magnitude
  !!           of surface stress divided by density of air at the surface,
  !!           and `u_star*b_star` is the buoyancy flux at the surface.
  !!
  !! \throw FATAL, "must call atm_land_ice_flux_exchange_init first"
  !!    atm_land_ice_flux_exchange_init has not been called before calling sfc_boundary_layer.
  subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Land_Ice_Atmos_Boundary )
    real,                  intent(in)     :: dt !< Time step
    type(time_type),       intent(in)     :: Time !< Current time
    type(atmos_data_type), intent(inout)  :: Atm !< A derived data type to specify atmosphere boundary data
    type(land_data_type),  intent(inout)  :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),   intent(inout)  :: Ice !< A derived data type to specify ice boundary data
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary !< A derived data type to specify properties and
                                                                                 !! fluxes passed from exchange grid to the atmosphere,
                                                                                 !! land and ice

    ! ---- local vars ----------------------------------------------------------
    real, dimension(n_xgrid_sfc) :: &
         ex_albedo,     &
         ex_albedo_vis_dir,     &
         ex_albedo_nir_dir,     &
         ex_albedo_vis_dif,     &
         ex_albedo_nir_dif,     &
         ex_land_frac,  &
         ex_t_atm,      &
         ex_p_atm,      &
         ex_u_atm, ex_v_atm,    &
         ex_gust,       &
         ex_t_surf4,    &
         ex_u_surf, ex_v_surf,  &
         ex_rough_mom, ex_rough_heat, ex_rough_moist, &
         ex_rough_scale,&
         ex_q_star,     &
         ex_cd_q,       &
         ex_ref, ex_ref_u, ex_ref_v, ex_u10, &
         ex_ref2,       &
         ex_t_ref,      &
         ex_qs_ref,     &
         ex_qs_ref_cmip,     &
         ex_del_m,      &
         ex_del_h,      &
         ex_del_q,      &
         ex_frac_open_sea

    real, dimension(n_xgrid_sfc,n_exch_tr) :: ex_tr_atm
    ! jgj: added for co2_atm diagnostic
    real, dimension(n_xgrid_sfc)           :: ex_co2_atm_dvmr
    real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: diag_atm
#ifndef _USE_LEGACY_LAND_
    real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2)) :: diag_land
    real, dimension(size(Land%t_ca, 1))                   :: diag_land_ug, tile_size_ug
    real, dimension(nxc_lnd,nyc_lnd)                      :: diag_land_sg, tile_size_sg
    logical, dimension(size(Land%t_ca, 1))                :: mask_ug
    logical, dimension(nxc_lnd,nyc_lnd)                   :: mask_sg
    integer :: k
#else
    real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
#endif
    real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2),size(Ice%t_surf,3)) :: sea
    real, dimension(size(Ice%albedo,1),size(Ice%albedo,2),size(Ice%albedo,3)) ::  tmp_open_sea
    real    :: zrefm, zrefh
    logical :: used
    character(32) :: tr_name, tr_units ! tracer name
    integer :: tr, n, m ! tracer indices
    integer :: i
    integer :: is,ie,l,j
    integer :: isc,iec,jsc,jec

    ! [1] check that the module was initialized
    if (do_init) call error_mesg ('atm_land_ice_flux_exchange_mod',  &
         'must call atm_land_ice_flux_exchange_init first', FATAL)
    !Balaji
    call mpp_clock_begin(cplClock)
    call mpp_clock_begin(sfcClock)
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

    ex_p_surf = 1.0
    ! Actual allocation of exchange fields for ocean_ice boundary
    do n = 1, ex_gas_fields_ice%num_bcs  !{
       do m = 1, ex_gas_fields_ice%bc(n)%num_fields  !{
          if (associated(ex_gas_fields_ice%bc(n)%field(m)%values)) then  !{
             call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_ice already allocated.' )
          endif  !}
          allocate ( ex_gas_fields_ice%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fields_ice%bc(n)%field(m)%values = 0.0
       enddo  !} m
    enddo  !} n

    do n = 1, ex_gas_fields_atm%num_bcs  !{
       do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
          if (associated(ex_gas_fields_atm%bc(n)%field(m)%values)) then  !{
             call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fields_atm already allocated.' )
          endif  !}
          allocate ( ex_gas_fields_atm%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fields_atm%bc(n)%field(m)%values = 0.0
       enddo  !} m
    enddo  !} n

    do n = 1, ex_gas_fluxes%num_bcs  !{
       do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
          if (associated(ex_gas_fluxes%bc(n)%field(m)%values)) then  !{
             call mpp_error( FATAL, 'sfc_boundary_layer: ex_gas_fluxes already allocated.' )
          endif  !}
          allocate ( ex_gas_fluxes%bc(n)%field(m)%values(n_xgrid_sfc) )
          ex_gas_fluxes%bc(n)%field(m)%values = 0.0
       enddo  !} m
    enddo  !} n

    !
    !       Call the atmosphere tracer driver to gather the data needed for extra gas tracers
    ! For ocean only model

    !  call atmos_get_fields_for_flux(Atm)

    ! [3] initialize some values on exchange grid: this is actually a safeguard
    ! against using undefined values
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_t_surf,ex_u_surf, &
    !$OMP                                  ex_v_surf,ex_albedo,ex_albedo_vis_dir,ex_albedo_nir_dir, &
    !$OMP                                  ex_albedo_vis_dif,ex_albedo_nir_dif,ex_cd_t,ex_cd_m,  &
    !$OMP                                  ex_cd_q,ex_frac_open_sea)                             &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is,ie
          ex_t_surf(i)   = 200.
          ex_u_surf(i)   =   0.
          ex_v_surf(i)   =   0.
          ex_albedo(i) = 0. ! bw
          ex_albedo_vis_dir(i) = 0.
          ex_albedo_nir_dir(i) = 0.
          ex_albedo_vis_dif(i) = 0.
          ex_albedo_nir_dif(i) = 0.

          !---- do not use if relax time /= 0 ----
          ex_cd_t(i) = 0.0
          ex_cd_m(i) = 0.0
          ex_cd_q(i) = 0.0
          ex_frac_open_sea(i) =0.
       enddo
    enddo
    !-----------------------------------------------------------------------
    !Balaji: data_override stuff moved from coupler_main
    call data_override ('ATM', 't_bot',  Atm%t_bot , Time)
    call data_override ('ATM', 'z_bot',  Atm%z_bot , Time)
    call data_override ('ATM', 'p_bot',  Atm%p_bot , Time)
    call data_override ('ATM', 'u_bot',  Atm%u_bot , Time)
    call data_override ('ATM', 'v_bot',  Atm%v_bot , Time)
    call data_override ('ATM', 'p_surf', Atm%p_surf, Time)
    call data_override ('ATM', 'slp',    Atm%slp,    Time)
    call data_override ('ATM', 'gust',   Atm%gust,   Time)
    !
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
    !
    do tr = 1,n_atm_tr
       call get_tracer_names( MODEL_ATMOS, tr, tr_name )
       call data_override('ATM', trim(tr_name)//'_bot', Atm%tr_bot(:,:,tr), Time, override=used)
       ! conversion for land co2 data override from dry vmr to moist mmr
       if (used .and. lowercase(trim(tr_name)).eq.'co2') then
          ! 2017/08/08 jgj add co2_bot diagnostic in dry_vmr units
          if ( id_co2_bot > 0 ) used = send_data ( id_co2_bot, Atm%tr_bot(:,:,tr), Time )

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

    ! data override for co2 to be passed to ocean (co2_flux_pcair_atm)
    ! atmos_co2.F90 already called: converts tr_bot passed to ocean via gas_flux
    ! from moist mmr to dry vmr.
    do n = 1, atm%fields%num_bcs  !{
       do m = 1, atm%fields%bc(n)%num_fields  !{
          call data_override('ATM', atm%fields%bc(n)%field(m)%name,      &
               atm%fields%bc(n)%field(m)%values, Time, override = atm%fields%bc(n)%field(m)%override)
          ex_gas_fields_atm%bc(n)%field(m)%override = atm%fields%bc(n)%field(m)%override
          ! 2017/08/08 jgj add co2_flux_pcair_atm diagnostic
          if ( atm%fields%bc(n)%field(m)%override .and. lowercase(trim(atm%fields%bc(n)%field(m)%name)) .eq. 'co2_flux_pcair_atm') then
             if( id_co2_flux_pcair_atm > 0 ) used = send_data ( id_co2_flux_pcair_atm, atm%fields%bc(n)%field(m)%values, Time )
          endif
          ! 2017/08/15 jgj add o2_flux_pcair_atm diagnostic
          if ( atm%fields%bc(n)%field(m)%override .and. lowercase(trim(atm%fields%bc(n)%field(m)%name)) .eq. 'o2_flux_pcair_atm') then
             if( id_o2_flux_pcair_atm > 0 ) used = send_data ( id_o2_flux_pcair_atm, atm%fields%bc(n)%field(m)%values, Time )
          endif
       enddo  !} m
    enddo  !} n
    do n = 1, atm%fields%num_bcs  !{
       if (atm%fields%bc(n)%use_atm_pressure) then  !{
          if (.not. atm%fields%bc(n)%field(ind_psurf)%override) then  !{
             atm%fields%bc(n)%field(ind_psurf)%values = Atm%p_surf
          endif  !}
       endif  !}
    enddo  !} n
    call data_override ('ICE', 't_surf',     Ice%t_surf,      Time)
    call data_override ('ICE', 'rough_mom',  Ice%rough_mom,   Time)
    call data_override ('ICE', 'rough_heat', Ice%rough_heat,  Time)
    call data_override ('ICE', 'rough_moist',Ice%rough_moist, Time)
    call data_override ('ICE', 'albedo',     Ice%albedo,      Time)
    call data_override ('ICE', 'albedo_vis_dir', Ice%albedo_vis_dir, Time)
    call data_override ('ICE', 'albedo_nir_dir', Ice%albedo_nir_dir, Time)
    call data_override ('ICE', 'albedo_vis_dif', Ice%albedo_vis_dif, Time)
    call data_override ('ICE', 'albedo_nir_dif', Ice%albedo_nir_dif, Time)
    call data_override ('ICE', 'u_surf',     Ice%u_surf,      Time)
    call data_override ('ICE', 'v_surf',     Ice%v_surf,      Time)
    call coupler_type_data_override('ICE', Ice%ocean_fields, Time)
    call coupler_type_send_data(Ice%ocean_fields, Time)
#ifndef _USE_LEGACY_LAND_
    call data_override_ug ('LND', 't_surf',     Land%t_surf,     Time)
    call data_override_ug ('LND', 't_ca',       Land%t_ca,       Time)
    call data_override_ug ('LND', 'rough_mom',  Land%rough_mom,  Time)
    call data_override_ug ('LND', 'rough_heat', Land%rough_heat, Time)
    call data_override_ug ('LND', 'albedo', Land%albedo,     Time)
#else
    call data_override ('LND', 't_surf',     Land%t_surf,     Time)
    call data_override ('LND', 't_ca',       Land%t_ca,       Time)
    call data_override ('LND', 'rough_mom',  Land%rough_mom,  Time)
    call data_override ('LND', 'rough_heat', Land%rough_heat, Time)
    call data_override ('LND', 'albedo', Land%albedo,     Time)
#endif

    ! tracer data override
    do tr = 1, n_lnd_tr
       call get_tracer_names( MODEL_LAND, tr, tr_name )
#ifndef _USE_LEGACY_LAND_
       call data_override_ug('LND', trim(tr_name)//'_surf', Land%tr(:,:,tr), Time)
    enddo
    call data_override_ug ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
    call data_override_ug ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
    call data_override_ug ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
    call data_override_ug ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)
#else
       call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
    enddo
    call data_override ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
    call data_override ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
    call data_override ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
    call data_override ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)
#endif

    !---- put atmosphere quantities onto exchange grid ----

    ! [4] put all the qantities we need onto exchange grid
    ! [4.1] put atmosphere quantities onto exchange grid
#ifdef use_AM3_physics
    if (do_forecast) then
       call put_to_xgrid (Atm%Surf_diff%sst_miz , 'ATM', ex_t_surf_miz, xmap_sfc, remap_method=remap_method, complete=.false.)
    endif
#endif
    ! put atmosphere bottom layer tracer data onto exchange grid
    do tr = 1,n_exch_tr
       call put_to_xgrid (Atm%tr_bot(:,:,tr_table(tr)%atm) , 'ATM', ex_tr_atm(:,tr), xmap_sfc, &
            remap_method=remap_method, complete=.false.)
    enddo
    do n = 1, Atm%fields%num_bcs  !{
      if(ex_gas_fields_atm%bc(n)%flux_type  .ne. 'air_sea_deposition') then
       do m = 1, Atm%fields%bc(n)%num_fields  !{
          call put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
               ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method, complete=.false.)
       enddo  !} m
      endif
    enddo  !} n

    call put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%slp,    'ATM', ex_slp,    xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc, remap_method=remap_method, complete=.true.)

    ! slm, Mar 20 2002: changed order in whith the data transferred from ice and land
    ! grids, to fill t_ca first with t_surf over ocean and then with t_ca from
    ! land, where it is different from t_surf. It is mostly to simplify
    ! diagnostic, since surface_flux calculations distinguish between land and
    ! not-land anyway.

    ! prefill surface values with atmospheric values before putting tracers
    ! from ice or land, so that gradient is 0 if tracers are not filled
    ex_tr_surf = ex_tr_atm

    ! [4.2] put ice quantities onto exchange grid
    ! (assume that ocean quantites are stored in no ice partition)
    ! (note: ex_avail is true at ice and ocean points)
    call put_to_xgrid (Ice%t_surf,      'OCN', ex_t_surf,      xmap_sfc)
    call put_to_xgrid (Ice%rough_mom,   'OCN', ex_rough_mom,   xmap_sfc)
    call put_to_xgrid (Ice%rough_heat,  'OCN', ex_rough_heat,  xmap_sfc)
    call put_to_xgrid (Ice%rough_moist, 'OCN', ex_rough_moist, xmap_sfc)
    call put_to_xgrid (Ice%albedo,      'OCN', ex_albedo,      xmap_sfc)
    call put_to_xgrid (Ice%albedo_vis_dir, 'OCN', ex_albedo_vis_dir, xmap_sfc)
    call put_to_xgrid (Ice%albedo_nir_dir, 'OCN', ex_albedo_nir_dir, xmap_sfc)
    call put_to_xgrid (Ice%albedo_vis_dif, 'OCN', ex_albedo_vis_dif, xmap_sfc)
    call put_to_xgrid (Ice%albedo_nir_dif, 'OCN', ex_albedo_nir_dif, xmap_sfc)
    call put_to_xgrid (Ice%u_surf,      'OCN', ex_u_surf,      xmap_sfc)
    call put_to_xgrid (Ice%v_surf,      'OCN', ex_v_surf,      xmap_sfc)

    tmp_open_sea        = 0.
    tmp_open_sea(:,:,1) = 1.
    call put_to_xgrid ( tmp_open_sea,  'OCN', ex_frac_open_sea,   xmap_sfc)

    do n = 1, ice%ocean_fields%num_bcs  !{
       do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
          call put_to_xgrid (Ice%ocean_fields%bc(n)%field(m)%values, 'OCN',      &
               ex_gas_fields_ice%bc(n)%field(m)%values, xmap_sfc)
       enddo  !} m
    enddo  !} n
    sea = 0.0; sea(:,:,1) = 1.0;
    ex_seawater = 0.0
    call put_to_xgrid (sea,             'OCN', ex_seawater,    xmap_sfc)
    ex_t_ca = ex_t_surf ! slm, Mar 20 2002 to define values over the ocean

    ! [4.3] put land quantities onto exchange grid ----
    call some(xmap_sfc, ex_land, 'LND')

#ifndef _USE_LEGACY_LAND_

#ifdef use_AM3_physics
    if (do_forecast) then
       call put_to_xgrid_ug (Land%t_surf,     'LND', ex_t_surf_miz,  xmap_sfc)
       ex_t_ca(:) = ex_t_surf_miz(:)
    end if
#endif

    call put_to_xgrid_ug (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
    call put_to_xgrid_ug (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
    call put_to_xgrid_ug (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
    call put_to_xgrid_ug (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
    call put_to_xgrid_ug (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
    call put_to_xgrid_ug (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
    call put_to_xgrid_ug (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
    call put_to_xgrid_ug (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
    call put_to_xgrid_ug (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
    call put_to_xgrid_ug (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
    ex_rough_scale = ex_rough_mom
    call put_to_xgrid_ug(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)

    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
          call put_to_xgrid_ug ( Land%tr(:,:,n), 'LND', ex_tr_surf(:,tr), xmap_sfc )
       else
          ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
          ! fluxes will be 0
       endif
    enddo
#else

#ifdef use_AM3_physics
    if (do_forecast) then
       call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf_miz,  xmap_sfc)
       ex_t_ca(:) = ex_t_surf_miz(:)
    end if
#endif
    call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
    call put_to_xgrid (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
    call put_to_xgrid (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
    call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
    call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
    call put_to_xgrid (Land%albedo,     'LND', ex_albedo,      xmap_sfc)
    call put_to_xgrid (Land%albedo_vis_dir,     'LND', ex_albedo_vis_dir,   xmap_sfc)
    call put_to_xgrid (Land%albedo_nir_dir,     'LND', ex_albedo_nir_dir,   xmap_sfc)
    call put_to_xgrid (Land%albedo_vis_dif,     'LND', ex_albedo_vis_dif,   xmap_sfc)
    call put_to_xgrid (Land%albedo_nir_dif,     'LND', ex_albedo_nir_dif,   xmap_sfc)
    ex_rough_scale = ex_rough_mom
    call put_to_xgrid(Land%rough_scale, 'LND', ex_rough_scale, xmap_sfc)

    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
          call put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf(:,tr), xmap_sfc )
       else
          ! do nothing, since ex_tr_surf is prefilled with ex_tr_atm, and therefore
          ! fluxes will be 0
       endif
    enddo
#endif

    ex_land_frac = 0.0
    call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)

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

    ! [5] compute explicit fluxes and tendencies at all available points ---
    call some(xmap_sfc, ex_avail)
    !$OMP parallel do default(none) shared(my_nblocks,ex_t_atm,ex_tr_atm,ex_u_atm,ex_v_atm, &
    !$OMP                                  ex_p_atm,ex_z_atm,ex_p_surf,ex_t_surf,ex_t_ca, &
    !$OMP                                  ex_tr_surf,ex_u_surf,ex_v_surf,ex_rough_mom, &
    !$OMP                                  ex_rough_heat,ex_rough_moist,ex_rough_scale,    &
    !$OMP                                  ex_gust,ex_flux_t,ex_flux_tr,ex_flux_lw, &
    !$OMP                                  ex_flux_u,ex_flux_v,ex_cd_m,ex_cd_t,ex_cd_q, &
    !$OMP                                  ex_wind,ex_u_star,ex_b_star,ex_q_star,       &
    !$OMP                                  ex_dhdt_surf,ex_dedt_surf,ex_dfdtr_surf,   &
    !$OMP                                  ex_drdt_surf,ex_dhdt_atm,ex_dfdtr_atm,   &
    !$OMP                                  ex_dtaudu_atm, ex_dtaudv_atm,dt,ex_land, &
    !$OMP                                  ex_seawater,ex_avail,block_start,block_end,isphum) &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call surface_flux (&
            ex_t_atm(is:ie), ex_tr_atm(is:ie,isphum),  ex_u_atm(is:ie), ex_v_atm(is:ie),  ex_p_atm(is:ie),  ex_z_atm(is:ie),  &
            ex_p_surf(is:ie),ex_t_surf(is:ie), ex_t_ca(is:ie),  ex_tr_surf(is:ie,isphum),                       &
            ex_u_surf(is:ie), ex_v_surf(is:ie),                                           &
            ex_rough_mom(is:ie), ex_rough_heat(is:ie), ex_rough_moist(is:ie), ex_rough_scale(is:ie),    &
            ex_gust(is:ie),                                                        &
            ex_flux_t(is:ie), ex_flux_tr(is:ie,isphum), ex_flux_lw(is:ie), ex_flux_u(is:ie), ex_flux_v(is:ie),         &
            ex_cd_m(is:ie),   ex_cd_t(is:ie), ex_cd_q(is:ie),                                    &
            ex_wind(is:ie),   ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),                     &
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
            ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
            ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,       &
            dt,                                                                        &
            ex_land, ex_seawater .gt. 0.0,  ex_avail,                                    &
            ex_dhdt_surf_forland,  ex_dedt_surf_forland,  ex_dedq_surf_forland  )

    endif
#endif

    !  call mpp_clock_end(fluxClock)
    zrefm = 10.0
    zrefh = z_ref_heat
    !      ---- optimize calculation ----
    !$OMP parallel do default(shared) private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call mo_profile ( zrefm, zrefh, ex_z_atm(is:ie), ex_rough_mom(is:ie), &
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
       do n = 1, ex_gas_fields_atm%num_bcs  !{
          if (atm%fields%bc(n)%use_10m_wind_speed) then  !{
             if (.not. ex_gas_fields_atm%bc(n)%field(ind_u10)%override) then  !{
                do i = is,ie
                   ex_gas_fields_atm%bc(n)%field(ind_u10)%values(i) = ex_u10(i)
                enddo
             endif  !}
          endif  !}
       enddo  !} n
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
    enddo ! end of block loop

    ! Combine explicit ocean flux and implicit land flux of extra flux fields.

    ! Calculate ocean explicit flux here
    call atmos_ocean_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater, ex_t_surf)

   do n = 1, ex_gas_fluxes%num_bcs  !{
      if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
         m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
         if (id_tr_mol_flux0(m) .gt. 0) then
            call get_from_xgrid (diag_atm, 'ATM', ex_gas_fluxes%bc(n)%field(ind_flux0)%values(:), xmap_sfc)
            used = send_data ( id_tr_mol_flux0(m), diag_atm, Time )
         end if
      end if
   end do


    ! The following statement is a concise version of what's following and worth
    ! looking into in the future.
    ! ex_flux_tr(:,itracer) = ex_gas_fluxes%bc(itracer_ocn)%field(ind_flux)%values(:)
    ! where(ex_seawater.gt.0) ex_flux_tr(:,itracer) = F_ocn
    !$OMP parallel do default(shared) private(is,ie,m,tr_units,tr_name)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do n = 1, ex_gas_fluxes%num_bcs  !{
          if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
             m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
             call get_tracer_names( MODEL_ATMOS, ex_gas_fluxes%bc(n)%atm_tr_index, tr_name, units=tr_units)
             do i = is,ie  !{
                if (ex_land(i)) cycle  ! over land, don't do anything
                ! on ocean or ice cells, flux is explicit therefore we zero derivatives.
                ex_dfdtr_atm(i,m)  = 0.0
                ex_dfdtr_surf(i,m) = 0.0
                if (ex_seawater(i)>0.0) then
                   if (lowercase(trim(tr_units)).eq."vmr") then
                      ! in mol/m2/s but from land model it should be in vmr * kg/m2/s
                      ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(ind_flux)%values(i) * WTMAIR*1.0e-3 &
                           / (1.-ex_tr_atm(i,isphum))
                   else
                   ! jgj: convert to kg co2/m2/sec for atm
                   ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(ind_flux)%values(i) * ex_gas_fluxes%bc(n)%mol_wt * 1.0e-03
                   end if
                else
                   ex_flux_tr(i,m) = 0.0 ! pure ice exchange cell
                endif  !}
             enddo  !} i
          endif  !}
       enddo  !} n
    enddo ! l

    ! [5.2] override tracer fluxes and derivatives
    do tr = 1,n_exch_tr
       if( tr_table(tr)%atm == NO_TRACER ) cycle ! it should never happen, though

       call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       ! [5.2.1] override tracer flux. Note that "sea" and "diag_land" are repeatedly used
       ! as temporary storage for the values we are overriding fluxes and derivative with,
       ! over ocean and land respectively
#ifndef _USE_LEGACY_LAND_
       call data_override_ug ( 'LND', 'ex_flux_'//trim(tr_name), diag_land, Time, override=used )
       if(used) call put_to_xgrid_ug ( diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_flux_'//trim(tr_name), sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_flux_tr(:,tr), xmap_sfc )
       ! [5.2.2] override derivative of flux wrt surface concentration
       call data_override_ug ( 'LND', 'ex_dfd'//trim(tr_name)//'_surf', diag_land, Time, override=used )
       if(used) call put_to_xgrid_ug ( diag_land, 'LND', ex_dfdtr_surf(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_surf', sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_surf(:,tr), xmap_sfc )
       ! [5.2.3] override derivative of flux wrt atmospheric concentration
       call data_override_ug ( 'LND', 'ex_dfd'//trim(tr_name)//'_atm', diag_land, Time, override=used )
       if(used) call put_to_xgrid_ug ( diag_land, 'LND', ex_dfdtr_atm(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_atm', sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_atm(:,tr), xmap_sfc )
    enddo

    ! [5.3] override flux and derivatives for sensible heat flux
    ! [5.3.1] override flux
    call data_override_ug ( 'LND', 'ex_flux_t', diag_land, Time, override=used )
    if (used) call put_to_xgrid_ug ( diag_land, 'LND', ex_flux_t, xmap_sfc )
    call data_override ( 'ICE', 'ex_flux_t', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_flux_t, xmap_sfc )
    ! [5.3.2] override derivative of flux wrt near-surface temperature
    call data_override_ug ( 'LND', 'ex_dhdt_surf', diag_land, Time, override=used )
    if (used) call put_to_xgrid_ug ( diag_land, 'LND', ex_dhdt_surf, xmap_sfc )
    call data_override ( 'ICE', 'ex_dhdt_surf', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_surf, xmap_sfc )
    ! [5.3.3] override derivative of flux wrt atmospheric temperature
    call data_override_ug ( 'LND', 'ex_dhdt_atm', diag_land, Time,override=used )
    if (used) call put_to_xgrid_ug ( diag_land, 'LND', ex_dhdt_atm, xmap_sfc )
    call data_override ( 'ICE', 'ex_dhdt_atm', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_atm, xmap_sfc )
#else
       call data_override ( 'LND', 'ex_flux_'//trim(tr_name), diag_land, Time, override=used )
       if(used) call put_to_xgrid ( diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_flux_'//trim(tr_name), sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_flux_tr(:,tr), xmap_sfc )
       ! [5.2.2] override derivative of flux wrt surface concentration
       call data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_surf', diag_land, Time, override=used )
       if(used) call put_to_xgrid ( diag_land, 'LND', ex_dfdtr_surf(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_surf', sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_surf(:,tr), xmap_sfc )
       ! [5.2.3] override derivative of flux wrt atmospheric concentration
       call data_override ( 'LND', 'ex_dfd'//trim(tr_name)//'_atm', diag_land, Time, override=used )
       if(used) call put_to_xgrid ( diag_land, 'LND', ex_dfdtr_atm(:,tr), xmap_sfc )
       call data_override ( 'ICE', 'ex_dfd'//trim(tr_name)//'_atm', sea, Time, override=used )
       if(used) call put_to_xgrid ( sea, 'OCN', ex_dfdtr_atm(:,tr), xmap_sfc )
    enddo

    ! [5.3] override flux and derivatives for sensible heat flux
    ! [5.3.1] override flux
    call data_override ( 'LND', 'ex_flux_t', diag_land, Time, override=used )
    if (used) call put_to_xgrid ( diag_land, 'LND', ex_flux_t, xmap_sfc )
    call data_override ( 'ICE', 'ex_flux_t', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_flux_t, xmap_sfc )
    ! [5.3.2] override derivative of flux wrt near-surface temperature
    call data_override ( 'LND', 'ex_dhdt_surf', diag_land, Time, override=used )
    if (used) call put_to_xgrid ( diag_land, 'LND', ex_dhdt_surf, xmap_sfc )
    call data_override ( 'ICE', 'ex_dhdt_surf', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_surf, xmap_sfc )
    ! [5.3.3] override derivative of flux wrt atmospheric temperature
    call data_override ( 'LND', 'ex_dhdt_atm', diag_land, Time,override=used )
    if (used) call put_to_xgrid ( diag_land, 'LND', ex_dhdt_atm, xmap_sfc )
    call data_override ( 'ICE', 'ex_dhdt_atm', sea, Time, override=used )
    if (used) call put_to_xgrid ( sea, 'OCN', ex_dhdt_atm, xmap_sfc )
#endif

    ! NB: names of the override fields are constructed using tracer name and certain
    ! prefixes / suffixes. For example, for the tracer named "sphum" (specific humidity) they will be:
    ! "ex_flux_sphum", "ex_dfdsphum_surf", and "ex_dfdsphum_atm".
    !
    ! For sensible heat flux names are "ex_flux_t", "ex_dhdt_surf", and "ex_dhdt_atm";
    ! despite the name those are actually in energy units, W/m2, W/(m2 degK), and
    ! W/(m2 degK) respectively

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

    ! [6.2] put relevant quantities onto atmospheric boundary
    call get_from_xgrid (Land_Ice_Atmos_Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%frac_open_sea,'ATM',ex_frac_open_sea, xmap_sfc)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir,    'ATM',   &
         ex_albedo_vis_dir   ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir,    'ATM',   &
         ex_albedo_nir_dir   ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif,    'ATM',   &
         ex_albedo_vis_dif   ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif,    'ATM',   &
         ex_albedo_nir_dif   ,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc, complete=.false.)

    call get_from_xgrid (Land_Ice_Atmos_Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc, complete=.false.)
    call get_from_xgrid (Land_Ice_Atmos_Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc, complete=.true.)

    call get_from_xgrid (Land_Ice_Atmos_Boundary%u_ref,     'ATM', ex_ref_u     , xmap_sfc, complete=.false.) !bqx
    call get_from_xgrid (Land_Ice_Atmos_Boundary%v_ref,     'ATM', ex_ref_v     , xmap_sfc, complete=.true.) !bqx

#ifdef use_AM3_physics
    if (do_forecast) then
       call get_from_xgrid (Ice%t_surf, 'OCN', ex_t_surf,  xmap_sfc)
    end if
#endif

    call mpp_get_compute_domain( Atm%domain, isc, iec, jsc, jec )
    !$OMP parallel do default(none) shared(isc,iec,jsc,jec,Land_Ice_Atmos_Boundary ) &
    !$OMP                          private(is,ie)
    do j = jsc, jec
       do i = isc, iec
          Land_Ice_Atmos_Boundary%t(i,j) = Land_Ice_Atmos_Boundary%t(i,j) ** 0.25
       enddo
    enddo
    !Balaji: data_override calls moved here from coupler_main
    call data_override('ATM', 't',         Land_Ice_Atmos_Boundary%t,         Time)
    call data_override('ATM', 'albedo',    Land_Ice_Atmos_Boundary%albedo,    Time)

    call data_override('ATM', 'albedo_vis_dir',    Land_Ice_Atmos_Boundary%albedo_vis_dir,    Time)
    call data_override('ATM', 'albedo_nir_dir',    Land_Ice_Atmos_Boundary%albedo_nir_dir,    Time)
    call data_override('ATM', 'albedo_vis_dif',    Land_Ice_Atmos_Boundary%albedo_vis_dif,    Time)
    call data_override('ATM', 'albedo_nir_dif',    Land_Ice_Atmos_Boundary%albedo_nir_dif,    Time)
    call data_override('ATM', 'land_frac', Land_Ice_Atmos_Boundary%land_frac, Time)
    call data_override('ATM', 'dt_t',      Land_Ice_Atmos_Boundary%dt_t,      Time)
    do tr=1,n_atm_tr
       call get_tracer_names(MODEL_ATMOS, tr, tr_name)
       call data_override('ATM', 'dt_'//trim(tr_name), Land_Ice_Atmos_Boundary%dt_tr(:,:,tr), Time)
    enddo
    call data_override('ATM', 'u_flux',    Land_Ice_Atmos_Boundary%u_flux,    Time)
    call data_override('ATM', 'v_flux',    Land_Ice_Atmos_Boundary%v_flux,    Time)
    call data_override('ATM', 'dtaudu',    Land_Ice_Atmos_Boundary%dtaudu,    Time)
    call data_override('ATM', 'dtaudv',    Land_Ice_Atmos_Boundary%dtaudv,    Time)
    call data_override('ATM', 'u_star',    Land_Ice_Atmos_Boundary%u_star,    Time)
    call data_override('ATM', 'b_star',    Land_Ice_Atmos_Boundary%b_star,    Time)
    ! call data_override('ATM', 'q_star',    Land_Ice_Atmos_Boundary%q_star,    Time)
    call data_override('ATM', 'rough_mom', Land_Ice_Atmos_Boundary%rough_mom, Time)

    ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
    ! on exchange grid
    ! allocate ( ex_old_albedo(n_xgrid_sfc)  )
    ! ex_old_albedo = ex_albedo

    !!  STILL NEEDED   ????
    !! IS THIS CORRECT ??
    allocate ( ex_albedo_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
    allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )

    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_albedo_fix, &
    !$OMP                                  ex_albedo_vis_dir_fix,ex_albedo_nir_dir_fix,    &
    !$OMP                                  ex_albedo_vis_dif_fix,ex_albedo_nir_dif_fix )   &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_fix(i) = 0.
          ex_albedo_vis_dir_fix(i) = 0.
          ex_albedo_nir_dir_fix(i) = 0.
          ex_albedo_vis_dif_fix(i) = 0.
          ex_albedo_nir_dif_fix(i) = 0.
       enddo
    enddo


    call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo, 'ATM',  ex_albedo_fix, xmap_sfc, complete=.false.)
    call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir, 'ATM',  &
         ex_albedo_vis_dir_fix, xmap_sfc, complete=.false.)
    call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir, 'ATM', &
         ex_albedo_nir_dir_fix, xmap_sfc, complete=.false.)
    call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif, 'ATM',   &
         ex_albedo_vis_dif_fix, xmap_sfc, complete=.false.)
    call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif, 'ATM',  &
         ex_albedo_nir_dif_fix, xmap_sfc, complete=.true.)
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_albedo_fix,    &
    !$OMP                                  ex_albedo,ex_albedo_vis_dir_fix,ex_albedo_vis_dir, &
    !$OMP                                  ex_albedo_nir_dir,ex_albedo_nir_dir_fix,           &
    !$OMP                                  ex_albedo_vis_dif_fix,ex_albedo_vis_dif,           &
    !$OMP                                  ex_albedo_nir_dif_fix,ex_albedo_nir_dif)           &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is, ie
          ex_albedo_fix(i) = (1.0-ex_albedo(i)) / (1.0-ex_albedo_fix(i))
          ex_albedo_vis_dir_fix(i) = (1.0-ex_albedo_vis_dir(i)) / (1.0-ex_albedo_vis_dir_fix(i))
          ex_albedo_nir_dir_fix(i) = (1.0-ex_albedo_nir_dir(i)) / (1.0-ex_albedo_nir_dir_fix(i))
          ex_albedo_vis_dif_fix(i) = (1.0-ex_albedo_vis_dif(i)) / (1.0-ex_albedo_vis_dif_fix(i))
          ex_albedo_nir_dif_fix(i) = (1.0-ex_albedo_nir_dif(i)) / (1.0-ex_albedo_nir_dif_fix(i))
       enddo
    enddo

#ifdef SCM
    if (do_specified_albedo .and. do_specified_land) then
       ex_albedo_fix = 1.
       ex_albedo_vis_dir_fix = 1.
       ex_albedo_vis_dif_fix = 1.
       ex_albedo_nir_dir_fix = 1.
       ex_albedo_nir_dif_fix = 1.
    endif
#endif

    !=======================================================================
    ! [7] diagnostics section

    !------- save static fields first time only ------
    if (first_static) then

       !------- land fraction ------
       if ( id_land_mask > 0 ) then
          used = send_data ( id_land_mask, Land_Ice_Atmos_Boundary%land_frac, Time )
       endif
       if ( id_sftlf > 0 ) then
          used = send_data ( id_sftlf, Land_Ice_Atmos_Boundary%land_frac, Time )
       endif
       ! near-surface heights
       if ( id_height2m  > 0) used = send_data ( id_height2m, z_ref_heat, Time )
       if ( id_height10m > 0) used = send_data ( id_height10m, z_ref_mom, Time )

       first_static = .false.
    endif

    !------- Atm fields -----------
    do n = 1, Atm%fields%num_bcs  !{
       do m = 1, Atm%fields%bc(n)%num_fields  !{
          if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
             if (atm%fields%bc(n)%use_10m_wind_speed .and. m .eq. ind_u10 .and. .not. Atm%fields%bc(n)%field(m)%override) then  !{
                call get_from_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',     &
                     ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc)
             endif  !}
             if ( Atm%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
                used = send_data(Atm%fields%bc(n)%field(m)%id_diag, Atm%fields%bc(n)%field(m)%values, Time )
             endif  !}
          endif  !}
       enddo  !} m
    enddo  !} n

    !------- drag coeff moisture -----------
    if ( id_wind > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_wind, xmap_sfc)
       used = send_data ( id_wind, diag_atm, Time )
    endif
    !------- drag coeff moisture -----------
    if ( id_drag_moist > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_cd_q, xmap_sfc)
       used = send_data ( id_drag_moist, diag_atm, Time )
    endif

    !------- drag coeff heat -----------
    if ( id_drag_heat > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_cd_t, xmap_sfc)
       used = send_data ( id_drag_heat, diag_atm, Time )
    endif

    !------- drag coeff momemtum -----------
    if ( id_drag_mom > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_cd_m, xmap_sfc)
       used = send_data ( id_drag_mom, diag_atm, Time )
    endif

    !------- roughness moisture -----------
    if ( id_rough_moist > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_rough_moist, xmap_sfc)
       used = send_data ( id_rough_moist, diag_atm, Time )
    endif

    !------- roughness heat -----------
    if ( id_rough_heat > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_rough_heat, xmap_sfc)
       used = send_data ( id_rough_heat, diag_atm, Time )
    endif

    !------- roughness momemtum -----------
    used = send_data ( id_rough_mom, Land_Ice_Atmos_Boundary%rough_mom, Time )

    !------- friction velocity -----------
    used = send_data ( id_u_star, Land_Ice_Atmos_Boundary%u_star, Time )

    !------- bouyancy -----------
    used = send_data ( id_b_star, Land_Ice_Atmos_Boundary%b_star, Time )

    !------- moisture scale -----------
    used = send_data ( id_q_star, Land_Ice_Atmos_Boundary%q_star, Time )

    !-----------------------------------------------------------------------
    !------ diagnostics for fields at bottom atmospheric level ------

    if ( id_t_atm > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_t_atm, xmap_sfc)
       used = send_data ( id_t_atm, diag_atm, Time )
    endif

    if ( id_u_atm > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_u_atm, xmap_sfc)
       used = send_data ( id_u_atm, diag_atm, Time )
    endif

    if ( id_v_atm > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_v_atm, xmap_sfc)
       used = send_data ( id_v_atm, diag_atm, Time )
    endif

    do tr = 1,n_exch_tr
       call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       if ( id_tr_atm(tr) > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_tr_atm(:,tr), xmap_sfc)
          used = send_data ( id_tr_atm(tr), diag_atm, Time )
       endif
       !!jgj: add dryvmr co2_atm
       ! - slm Mar 25 2010: moved to resolve interdependence of diagnostic fields
       if ( id_co2_atm_dvmr > 0 .and. lowercase(trim(tr_name))=='co2') then
          ex_co2_atm_dvmr = (ex_tr_atm(:,tr) / (1.0 - ex_tr_atm(:,isphum))) * WTMAIR/WTMCO2
          call get_from_xgrid (diag_atm, 'ATM', ex_co2_atm_dvmr, xmap_sfc)
          used = send_data ( id_co2_atm_dvmr, diag_atm, Time )
       endif
    enddo

    ! - slm, Mar 25, 2002
    if ( id_p_atm > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_p_atm, xmap_sfc)
       used = send_data ( id_p_atm, diag_atm, Time )
    endif
    if ( id_z_atm > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_z_atm, xmap_sfc)
       used = send_data ( id_z_atm, diag_atm, Time )
    endif
    if ( id_gust > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_gust, xmap_sfc)
       used = send_data ( id_gust, diag_atm, Time )
    endif

    ! - bw, Sep 17, 2007
    if ( id_slp > 0 .or. id_psl > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_slp, xmap_sfc)
       if ( id_slp > 0 ) used = send_data ( id_slp, diag_atm, Time )
       if ( id_psl > 0 ) used = send_data ( id_psl, diag_atm, Time )
    endif

    !-----------------------------------------------------------------------
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
    !$OMP                                  ex_avail,ex_ref,ex_tr_surf,ex_tr_atm,isphum)           &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       call mo_profile ( zrefm, zrefh, ex_z_atm(is:ie),   ex_rough_mom(is:ie), &
            ex_rough_heat(is:ie), ex_rough_moist(is:ie),          &
            ex_u_star(is:ie), ex_b_star(is:ie), ex_q_star(is:ie),        &
            ex_del_m(is:ie), ex_del_h(is:ie), ex_del_q(is:ie), ex_avail(is:ie)  )

       !    ------- reference relative humidity -----------
       !cjg     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
       !cjg          id_rh_ref_cmip > 0 .or. &
       !cjg          id_q_ref > 0 .or. id_q_ref_land >0 ) then
       do i = is,ie
          ex_ref(i) = 1.0e-06
          if (ex_avail(i)) &
               ex_ref(i)   = ex_tr_surf(i,isphum) + (ex_tr_atm(i,isphum)-ex_tr_surf(i,isphum)) * ex_del_q(i)
       enddo
    enddo
    call get_from_xgrid (Land_Ice_Atmos_Boundary%q_ref, 'ATM', ex_ref,   xmap_sfc)  ! cjg
    if(id_q_ref > 0) then
       used = send_data(id_q_ref,Land_Ice_Atmos_Boundary%q_ref,Time)
    endif
    if(id_huss > 0) then
       used = send_data(id_huss,Land_Ice_Atmos_Boundary%q_ref,Time)
    endif
    if(id_q_ref_land > 0 .or.id_hussLut_land > 0) then
!duplicate send_tile_data. We may remove id_q_ref_land in the future.
#ifndef _USE_LEGACY_LAND_
       call get_from_xgrid_ug (diag_land, 'LND', ex_ref, xmap_sfc)
       call send_tile_data (id_q_ref_land, diag_land)
       call send_tile_data (id_hussLut_land, diag_land)
#else
       call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
       used = send_tile_averaged_data(id_q_ref_land, diag_land, &
            Land%tile_size, Time, mask=Land%mask)
#endif
    endif
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_t_ref,ex_avail, &
    !$OMP                                  ex_t_ca,ex_t_atm,ex_p_surf,ex_qs_ref,ex_del_h,      &
    !$OMP                                  ex_ref,ex_qs_ref_cmip,ex_ref2 ) &
    !$OMP                          private(is,ie)
    do l = 1, my_nblocks
       is=block_start(l)
       ie=block_end(l)
       do i = is,ie
          ex_t_ref(i) = 200.
          if(ex_avail(i)) &
               ex_t_ref(i) = ex_t_ca(i) + (ex_t_atm(i)-ex_t_ca(i)) * ex_del_h(i)
       enddo
       call compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref(is:ie), q = ex_ref(is:ie))
       call compute_qs (ex_t_ref(is:ie), ex_p_surf(is:ie), ex_qs_ref_cmip(is:ie),  &
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

    call get_from_xgrid (Land_Ice_Atmos_Boundary%t_ref, 'ATM', ex_t_ref, xmap_sfc)  ! cjg

    if ( id_rh_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
       call get_from_xgrid_ug (diag_land,'LND', ex_ref, xmap_sfc)
       call send_tile_data (id_rh_ref_land, diag_land)
#else
       call get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
       used = send_tile_averaged_data ( id_rh_ref_land, diag_land, &
            Land%tile_size, Time, mask = Land%mask )
#endif
    endif

    if(id_rh_ref > 0) then
       call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
       used = send_data ( id_rh_ref, diag_atm, Time )
    endif

    if(id_rh_ref_cmip > 0 .or. id_hurs > 0 .or. id_rhs > 0) then
       call get_from_xgrid (diag_atm, 'ATM', ex_ref2, xmap_sfc)
       if (id_rh_ref_cmip > 0) used = send_data ( id_rh_ref_cmip, diag_atm, Time )
       if (id_hurs > 0)        used = send_data ( id_hurs, diag_atm, Time )
       if (id_rhs  > 0)        used = send_data ( id_rhs,  diag_atm, Time )
    endif
    !cjg  endif

    !    ------- reference temp -----------
#ifdef use_AM3_physics
    if ( id_t_ref > 0 .or. id_t_ref_land > 0 .or. id_tasLut_land > 0 ) then
       where (ex_avail) &
            ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
       if (id_t_ref_land > 0.or.id_tasLut_land > 0) then
#ifndef _USE_LEGACY_LAND_
          call get_from_xgrid_ug(diag_land, 'LND', ex_ref, xmap_sfc)
          if (id_t_ref_land > 0)  call send_tile_data (id_t_ref_land, diag_land)
          if (id_tasLut_land > 0) call send_tile_data (id_tasLut_land, diag_land)
#else
          call get_from_xgrid(diag_land, 'LND', ex_ref, xmap_sfc)
          if (id_t_ref_land > 0) used = send_tile_averaged_data ( id_t_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif

       if ( id_t_ref > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          used = send_data ( id_t_ref, diag_atm, Time )
       endif
    endif
#else
    where (ex_avail) &
         ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
    if (id_t_ref_land > 0 .or. id_tasLut_land > 0 .or. id_tasl_g > 0) then
       ! t_ref diagnostic at land points only
#ifndef _USE_LEGACY_LAND_
       call get_from_xgrid_ug (diag_land, 'LND', ex_ref, xmap_sfc)
       if (id_t_ref_land > 0)  call send_tile_data (id_t_ref_land, diag_land)
       if (id_tasLut_land > 0) call send_tile_data (id_tasLut_land, diag_land)
       if (id_tasl_g > 0) then
         used = send_global_land_diag ( get_global_diag_field_id(id_tasl_g), &
                                diag_land, Time, Land%tile_size, Land%mask, Land )
       endif
#else
       call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
       if (id_t_ref_land > 0) used = send_tile_averaged_data ( id_t_ref_land, diag_land, &
            Land%tile_size, Time, mask = Land%mask )
#endif
    endif
    ! t_ref diagnostic at all atmos points
    call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
    if ( id_t_ref > 0 ) used = send_data ( id_t_ref, diag_atm, Time )
    if ( id_tas > 0 )   used = send_data ( id_tas, diag_atm, Time )
    call sum_diag_integral_field ('t_ref',  diag_atm)
    if ( id_tas_g > 0 )  used = send_global_diag ( id_tas_g, diag_atm, Time )
#endif

    !    ------- reference u comp -----------
    if ( id_u_ref > 0 .or. id_u_ref_land > 0 .or. id_uas > 0) then
       where (ex_avail) &
            ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
       if ( id_u_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
          call get_from_xgrid_ug ( diag_land, 'LND', ex_ref, xmap_sfc )
          call send_tile_data ( id_u_ref_land, diag_land )
#else
          call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
          used = send_tile_averaged_data ( id_u_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif
       if ( id_u_ref > 0 .or. id_uas > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          if ( id_u_ref > 0 ) used = send_data ( id_u_ref, diag_atm, Time )
          if ( id_uas > 0 )   used = send_data ( id_uas, diag_atm, Time )
       endif
    endif

    !    ------- reference v comp -----------
    if ( id_v_ref > 0 .or. id_v_ref_land > 0 .or. id_vas > 0 ) then
       where (ex_avail) &
            ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
       if ( id_v_ref_land > 0 ) then
#ifndef _USE_LEGACY_LAND_
          call get_from_xgrid_ug ( diag_land, 'LND', ex_ref, xmap_sfc )
          call send_tile_data ( id_v_ref_land, diag_land )
#else
          call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
          used = send_tile_averaged_data ( id_v_ref_land, diag_land, &
               Land%tile_size, Time, mask = Land%mask )
#endif
       endif
       if ( id_v_ref > 0 .or. id_vas > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
          if ( id_v_ref > 0 ) used = send_data ( id_v_ref, diag_atm, Time )
          if ( id_vas   > 0 ) used = send_data ( id_vas, diag_atm, Time )
       endif
    endif

    !    ------- reference-level absolute wind -----------
    if ( id_wind_ref > 0 .or. id_sfcWind > 0 ) then
       where (ex_avail) &
            ex_ref = sqrt((ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m)**2 &
            +(ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m)**2)
       call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
       if ( id_wind_ref > 0 ) used = send_data ( id_wind_ref, diag_atm, Time )
       if ( id_sfcWind  > 0 ) used = send_data ( id_sfcWind, diag_atm, Time )
    endif

    !    ------- interp factor for heat ------
    if ( id_del_h > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_del_h, xmap_sfc)
       used = send_data ( id_del_h, diag_atm, Time )
    endif

    !    ------- interp factor for momentum ------
    if ( id_del_m > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_del_m, xmap_sfc)
       used = send_data ( id_del_m, diag_atm, Time )
    endif

    !    ------- interp factor for moisture ------
    if ( id_del_q > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_del_q, xmap_sfc)
       used = send_data ( id_del_q, diag_atm, Time )
    endif

    !cjg  endif
    ! topographic roughness scale
    if(id_rough_scale>0) then
       call get_from_xgrid (diag_atm, 'ATM',&
            (log(ex_z_atm/ex_rough_mom+1.0)/log(ex_z_atm/ex_rough_scale+1.0))**2, xmap_sfc)
       used = send_data(id_rough_scale, diag_atm, Time)
    endif

    !Balaji
    call mpp_clock_end(sfcClock)
    call mpp_clock_end(cplClock)

    !=======================================================================

  end subroutine sfc_boundary_layer

  !#######################################################################

  !> Returns fluxes and derivatives corrected for the implicit treatment of atmospheric
  !! diffusive fluxes, as well as the increments in the temperature and specific humidity
  !! of the lowest atmospheric layer due to all explicit processes as well as the diffusive
  !! fluxes through the top of this layer.
  !!
  !!
  !! The following elements from Atmos_boundary are used as input:
  !! <pre>
  !!        flux_u_atm = zonal wind stress (Pa)
  !!        flux_v_atm = meridional wind stress (Pa)
  !! </pre>
  !!
  !! The following elements of Land_boundary are output:
  !! <pre>
  !!       flux_t_land = sensible heat flux (W/m2)
  !!       flux_q_land = specific humidity flux (Kg/(m2 s)
  !!      flux_lw_land = net longwave flux (W/m2), uncorrected for
  !!                     changes in surface temperature
  !!      flux_sw_land = net shortwave flux (W/m2)
  !!         dhdt_land = derivative of sensible heat flux w.r.t.
  !!                     surface temperature (on land model grid)  (W/(m2 K)
  !!         dedt_land = derivative of specific humidity flux w.r.t.
  !!                     surface temperature (on land model grid)  (Kg/(m2 s K)
  !!         drdt_land = derivative of upward longwave flux w.r.t.
  !!                     surface temperature (on land model grid) (W/(m2 K)
  !!        lprec_land = liquid precipitation, mass for one time step
  !!                      (Kg/m2)
  !!        fprec_land = frozen precipitation, mass for one time step
  !!                      (Kg/m2)
  !! </pre>
  !!
  !! The following elements of Ice_boundary are output:
  !! <pre>
  !!        flux_u_ice = zonal wind stress (Pa)
  !!        flux_v_ice = meridional wind stress (Pa)
  !!        coszen_ice = cosine of the zenith angle
  !! </pre>
  subroutine flux_down_from_atmos (Time, Atm, Land, Ice, Atmos_boundary, Land_boundary, Ice_boundary )
    type(time_type),       intent(in)    :: Time !< Current time
    type(atmos_data_type), intent(inout) :: Atm  !< A derived data type to specify atmosphere boundary data
    type(land_data_type),  intent(in)    :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),   intent(in)    :: Ice  !< A derived data type to specify ice boundary data
    type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary !< A derived data type to specify properties and fluxes
                                                                    !! passed from exchange grid to the atmosphere, land and ice
    type(atmos_land_boundary_type),    intent(inout):: Land_boundary !< A derived data type to specify properties and fluxes
                                                                     !! passed from atmosphere to land
    type(atmos_ice_boundary_type),     intent(inout):: Ice_boundary !< A derived data type to specify properties and fluxes passed
                                                                    !! from atmosphere to ice

    real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
         ex_flux_sw_dir,  &
         ex_flux_sw_dif,  &
         ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
         ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
         ex_flux_sw_vis, &
         ex_flux_sw_vis_dir, &
         ex_flux_sw_vis_dif, &
         ex_lprec, ex_fprec,      &
         ex_tprec, & ! temperature of precipitation, currently equal to atm T
         ex_u_star_smooth,        &
#ifdef use_AM3_physics
    ex_coszen
#else
    ex_coszen, &
         ex_setl_flux, & ! tracer sedimentation flux from the lowest atm layer (positive down)
         ex_dsetl_dtr    ! and its derivative w.r.t. the tracer concentration
#endif
    real :: setl_flux(size(Atm%tr_bot,1),size(Atm%tr_bot,2))
    real :: dsetl_dtr(size(Atm%tr_bot,1),size(Atm%tr_bot,2))


    real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
         ex_delta_t, ex_delta_u, ex_delta_v, ex_dflux_t

    real, dimension(n_xgrid_sfc,n_exch_tr) :: &
         ex_delta_tr, & ! tracer tendencies
         ex_dflux_tr    ! fracer flux change

    real    :: cp_inv
    logical :: used
    logical :: ov
    integer :: ier
    integer :: is_atm, ie_atm, js_atm, je_atm, j

    character(32) :: tr_name ! name of the tracer
    integer :: tr, n, m ! tracer indices
    integer :: is, ie, l, i

    !Balaji
    call mpp_clock_begin(cplClock)
    call mpp_clock_begin(fluxAtmDnClock)
    ov = .FALSE.
    !-----------------------------------------------------------------------
    !Balaji: data_override calls moved here from coupler_main
    call data_override ('ATM', 'flux_sw',  Atm%flux_sw, Time)
    call data_override ('ATM', 'flux_sw_dir',  Atm%flux_sw_dir, Time)
    call data_override ('ATM', 'flux_sw_dif',  Atm%flux_sw_dif, Time)
    call data_override ('ATM', 'flux_sw_down_vis_dir',  Atm%flux_sw_down_vis_dir, Time)
    call data_override ('ATM', 'flux_sw_down_vis_dif',  Atm%flux_sw_down_vis_dif, Time)
    call data_override ('ATM', 'flux_sw_down_total_dir',  Atm%flux_sw_down_total_dir, Time)
    call data_override ('ATM', 'flux_sw_down_total_dif',  Atm%flux_sw_down_total_dif, Time)
    call data_override ('ATM', 'flux_sw_vis',  Atm%flux_sw_vis, Time)
    call data_override ('ATM', 'flux_sw_vis_dir',  Atm%flux_sw_vis_dir, Time)
    call data_override ('ATM', 'flux_sw_vis_dif',  Atm%flux_sw_vis_dif, Time)
    call data_override ('ATM', 'flux_lw',  Atm%flux_lw, Time)
    call data_override ('ATM', 'lprec',    Atm%lprec,   Time)

    if (scale_precip_2d) then
       call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       call data_override ('ATM', 'precip_scale2d',    frac_precip,   Time)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             Atm%lprec(i,j) = Atm%lprec(i,j)*frac_precip(i,j)
          enddo
       enddo
    endif

    if (partition_fprec_from_lprec .and. Atm%pe) then
       call mpp_get_compute_domain(Atm%Domain, is_atm, ie_atm, js_atm, je_atm)
       do j=js_atm,je_atm
          do i=is_atm, ie_atm
             if (Atm%t_bot(i,j) < tfreeze) then
                Atm%fprec(i,j) = Atm%lprec(i,j)
                Atm%lprec(i,j) = 0.0
             endif
          enddo
       enddo
    endif

    call data_override ('ATM', 'fprec',    Atm%fprec,   Time)
    call data_override ('ATM', 'coszen',   Atm%coszen,  Time)
    call data_override ('ATM', 'dtmass',   Atm%Surf_Diff%dtmass, Time)
    call data_override ('ATM', 'delta_t',  Atm%Surf_Diff%delta_t, Time)
    call data_override ('ATM', 'dflux_t',  Atm%Surf_Diff%dflux_t, Time)
    do tr = 1,n_atm_tr
       call get_tracer_names(MODEL_ATMOS,tr,tr_name)
       call data_override ('ATM', 'delta_'//trim(tr_name),  Atm%Surf_Diff%delta_tr(:,:,tr), Time)
       call data_override ('ATM', 'dflux_'//trim(tr_name),  Atm%Surf_Diff%dflux_tr(:,:,tr), Time)
    enddo

    !---- put atmosphere quantities onto exchange grid ----

    if(sw1way_bug) then
       call put_to_xgrid (Atm%flux_sw, 'ATM', ex_flux_sw, xmap_sfc, complete=.false.)
       call put_to_xgrid (Atm%flux_sw_vis, 'ATM', ex_flux_sw_vis, xmap_sfc, complete=.false.)
    end if
    !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_sw_dir, &
    !$OMP                                  ex_flux_sw_vis_dir,ex_flux_sw_dif,ex_delta_u,    &
    !$OMP                                  ex_flux_sw_vis_dif,ex_flux_lwd,ex_delta_v )      &
    !$OMP                          private(is,ie)
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
    enddo
    call put_to_xgrid (Atm%flux_sw_dir, 'ATM', ex_flux_sw_dir, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_vis_dir, 'ATM', ex_flux_sw_vis_dir, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_dif, 'ATM', ex_flux_sw_dif, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_vis_dif, 'ATM', ex_flux_sw_vis_dif, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_down_vis_dir, 'ATM', ex_flux_sw_down_vis_dir, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_down_total_dir, 'ATM', ex_flux_sw_down_total_dir, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_down_vis_dif, 'ATM', ex_flux_sw_down_vis_dif, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%flux_sw_down_total_dif, 'ATM', ex_flux_sw_down_total_dif, xmap_sfc, complete=.false.)

    !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
    !  if (mpp_pe()== mpp_root_pe()) print *,'LPREC', ccc

!!$  if(do_area_weighted_flux) then
!!$     call put_to_xgrid (Atm%lprec * AREA_ATM_MODEL,   'ATM', ex_lprec, xmap_sfc)
!!$     call put_to_xgrid (Atm%fprec * AREA_ATM_MODEL,   'ATM', ex_fprec, xmap_sfc)
!!$  else
    call put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc, complete=.false.)
    call put_to_xgrid (Atm%t_bot,   'ATM', ex_tprec, xmap_sfc, complete=.false.)
!!$  endif

    call put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc, complete=.true.)

    call put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc, remap_method=remap_method, complete=.false.)
    if(ex_u_star_smooth_bug) then
       call put_to_xgrid (Atmos_boundary%u_star, 'ATM', ex_u_star_smooth, xmap_sfc, remap_method=remap_method, complete=.false.)
       ex_u_star = ex_u_star_smooth
    endif


    ! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
    ! on exchange grid instead of the stresses themselves so that only the
    ! implicit corrections are filtered through the atmospheric grid not the
    ! stresses themselves
    call put_to_xgrid (Atm%Surf_Diff%delta_u, 'ATM', ex_delta_u, xmap_sfc, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%Surf_Diff%delta_v, 'ATM', ex_delta_v, xmap_sfc, remap_method=remap_method, complete=.true.)

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

    !-----------------------------------------------------------------------
    !---- adjust sw flux for albedo variations on exch grid ----
    !---- adjust 4 categories (vis/nir dir/dif) separately  ----
    if( sw1way_bug ) then ! to reproduce old results, may remove in the next major release.
       !-----------------------------------------------------------------------
       !---- adjust sw flux for albedo variations on exch grid ----
       !$OMP parallel do default(none) shared(my_nblocks,block_start,block_end,ex_flux_sw,ex_albedo_fix, &
       !$OMP                                  ex_flux_sw_vis,ex_albedo_vis_dir_fix,ex_flux_sw_dir,       &
       !$OMP                                  ex_albedo_vis_dif_fix,ex_flux_sw_vis_dir,ex_flux_sw_dif,   &
       !$OMP                                  ex_flux_sw_vis_dif)                                        &
       !$OMP                          private(is,ie)
       do l = 1, my_nblocks
          is=block_start(l)
          ie=block_end(l)
          do i = is, ie
             ex_flux_sw(i) = ex_flux_sw(i) * ex_albedo_fix(i)

             ex_flux_sw_vis(i) = ex_flux_sw_vis(i) * ex_albedo_vis_dir_fix(i)
             ex_flux_sw_dir(i) = ex_flux_sw_dir(i) * ex_albedo_vis_dir_fix(i)
             ex_flux_sw_dif(i) = ex_flux_sw_dif(i) * ex_albedo_vis_dif_fix(i)
             ex_flux_sw_vis_dir(i) = ex_flux_sw_vis_dir(i) * ex_albedo_vis_dir_fix(i)
             ex_flux_sw_vis_dif(i) = ex_flux_sw_vis_dif(i) * ex_albedo_vis_dif_fix(i)
          enddo
       enddo
    else
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
    end if

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

    deallocate ( ex_albedo_fix )
    deallocate ( ex_albedo_vis_dir_fix )
    deallocate ( ex_albedo_nir_dir_fix )
    deallocate ( ex_albedo_vis_dif_fix )
    deallocate ( ex_albedo_nir_dif_fix )

    !-----------------------------------------------------------------------
    !----- adjust fluxes for implicit dependence on atmosphere ----

    do tr = 1,n_exch_tr
       n = tr_table(tr)%atm
       call put_to_xgrid (Atm%Surf_Diff%delta_tr(:,:,n), 'ATM', ex_delta_tr(:,tr), xmap_sfc, complete=.false. )
       call put_to_xgrid (Atm%Surf_Diff%dflux_tr(:,:,n), 'ATM', ex_dflux_tr(:,tr), xmap_sfc, complete=.false. )
    enddo

    call put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc, complete=.false. )
    call put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc, complete=.false. )
    call put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc, complete=.true. )

#ifndef use_AM3_physics
    ! Get sedimentation flux. Has to be here (instead of sfc_boundary_layer sub)
    ! because of time stepping order: sedimentation fluxes are calculated in
    ! update_atmos_model_down (in atmos_tracer_driver), but sfc_boundary_layer
    ! is called before that.
    do tr = 1,n_exch_tr
       if (atmos_tracer_has_surf_setl_flux(tr_table(tr)%atm)) then
          call get_atmos_tracer_surf_setl_flux (tr_table(tr)%atm, setl_flux, dsetl_dtr)
          call put_to_xgrid(setl_flux, 'ATM', ex_setl_flux, xmap_sfc)
          call put_to_xgrid(dsetl_dtr, 'ATM', ex_dsetl_dtr, xmap_sfc)
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
       enddo ! i = is, ie
    enddo !  l = 1, my_nblocks
    !-----------------------------------------------------------------------
    !---- output fields on the land grid -------

#ifndef _USE_LEGACY_LAND_
    call get_from_xgrid_ug (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
#ifdef SCM
    if (do_specified_land .and. do_specified_flux) then
       call get_from_xgrid_ug (Land_boundary%dhdt,  'LND', ex_dhdt_surf_forland, xmap_sfc)
    else
       call get_from_xgrid_ug (Land_boundary%dhdt,  'LND', ex_dhdt_surf, xmap_sfc)
    endif
#else
    call get_from_xgrid_ug (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
#endif
    call get_from_xgrid_ug (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

    call get_from_xgrid_ug (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
    call get_from_xgrid_ug (Land_boundary%tprec,   'LND', ex_tprec,     xmap_sfc)
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
       call get_from_xgrid_ug (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
       call data_override_ug('LND', 'drag_q', Land_boundary%drag_q,  Time )
    endif
    if(associated(Land_boundary%lwdn_flux)) then
       call get_from_xgrid_ug (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
       call data_override_ug('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
    endif
    if(associated(Land_boundary%cd_m)) then
       call get_from_xgrid_ug (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
       call data_override_ug('LND', 'cd_m', Land_boundary%cd_m, Time )
    endif
    if(associated(Land_boundary%cd_t)) then
       call get_from_xgrid_ug (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
       call data_override_ug('LND', 'cd_t', Land_boundary%cd_t, Time )
    endif
    if(associated(Land_boundary%bstar)) then
       call get_from_xgrid_ug (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
       call data_override_ug('LND', 'bstar',  Land_boundary%bstar, Time )
    endif
    if(associated(Land_boundary%ustar)) then
       call get_from_xgrid_ug (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
       call data_override_ug('LND', 'ustar',  Land_boundary%ustar, Time )
    endif
    if(associated(Land_boundary%wind)) then
       call get_from_xgrid_ug (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
       call data_override_ug('LND', 'wind',  Land_boundary%wind, Time )
    endif
    if(associated(Land_boundary%z_bot)) then
       call get_from_xgrid_ug (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
       call data_override_ug('LND', 'z_bot',  Land_boundary%z_bot, Time )
    endif
#else
    call get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
    call get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
    call get_from_xgrid (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call get_from_xgrid (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
    call get_from_xgrid (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call get_from_xgrid (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
    call get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
#ifdef SCM
    if (do_specified_land .and. do_specified_flux) then
       call get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf_forland, xmap_sfc)
    else
       call get_from_xgrid (Land_boundary%dhdt,  'LND', ex_dhdt_surf, xmap_sfc)
    endif
#else
    call get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
#endif
    call get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
    call get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

    call get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
    call get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
    call get_from_xgrid (Land_boundary%tprec,   'LND', ex_tprec,     xmap_sfc)
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
       call get_from_xgrid (Land_boundary%drag_q, 'LND', ex_drag_q,    xmap_sfc)
       call data_override('LND', 'drag_q', Land_boundary%drag_q,  Time )
    endif
    if(associated(Land_boundary%lwdn_flux)) then
       call get_from_xgrid (Land_boundary%lwdn_flux, 'LND', ex_flux_lwd, xmap_sfc)
       call data_override('LND', 'lwdn_flux', Land_boundary%lwdn_flux, Time )
    endif
    if(associated(Land_boundary%cd_m)) then
       call get_from_xgrid (Land_boundary%cd_m, 'LND', ex_cd_m, xmap_sfc)
       call data_override('LND', 'cd_m', Land_boundary%cd_m, Time )
    endif
    if(associated(Land_boundary%cd_t)) then
       call get_from_xgrid (Land_boundary%cd_t, 'LND', ex_cd_t, xmap_sfc)
       call data_override('LND', 'cd_t', Land_boundary%cd_t, Time )
    endif
    if(associated(Land_boundary%bstar)) then
       call get_from_xgrid (Land_boundary%bstar, 'LND', ex_b_star, xmap_sfc)
       call data_override('LND', 'bstar',  Land_boundary%bstar, Time )
    endif
    if(associated(Land_boundary%ustar)) then
       call get_from_xgrid (Land_boundary%ustar, 'LND', ex_u_star, xmap_sfc)
       call data_override('LND', 'ustar',  Land_boundary%ustar, Time )
    endif
    if(associated(Land_boundary%wind)) then
       call get_from_xgrid (Land_boundary%wind, 'LND', ex_wind, xmap_sfc)
       call data_override('LND', 'wind',  Land_boundary%wind, Time )
    endif
    if(associated(Land_boundary%z_bot)) then
       call get_from_xgrid (Land_boundary%z_bot, 'LND', ex_z_atm, xmap_sfc)
       call data_override('LND', 'z_bot',  Land_boundary%z_bot, Time )
    endif
#endif
    Land_boundary%tr_flux = 0.0
    Land_boundary%dfdtr = 0.0
    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then

#ifndef _USE_LEGACY_LAND_
          call get_from_xgrid_ug (Land_boundary%tr_flux(:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
          call get_from_xgrid_ug (Land_boundary%dfdtr(:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
#else
          call get_from_xgrid (Land_boundary%tr_flux(:,:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
          call get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
#endif
#ifdef SCM
          if (do_specified_land .and. do_specified_flux .and. tr.eq.isphum) then
#ifndef _USE_LEGACY_LAND_
             call get_from_xgrid_ug (Land_boundary%dfdtr(:,:,n),   'LND', ex_dedq_surf_forland(:), xmap_sfc)
#else
             call get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dedq_surf_forland(:), xmap_sfc)
#endif
          endif
#endif
       endif
    enddo

    !  current time is Time: is that ok? not available in land_data_type
    !Balaji: data_override calls moved here from coupler_main
#ifndef _USE_LEGACY_LAND_
    call data_override_ug('LND', 't_flux',  Land_boundary%t_flux,  Time )
    call data_override_ug('LND', 'lw_flux', Land_boundary%lw_flux, Time )
    call data_override_ug('LND', 'sw_flux', Land_boundary%sw_flux, Time )
    call data_override_ug('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
    call data_override_ug('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
    call data_override_ug('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
    call data_override_ug('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )

    call data_override_ug('LND', 'lprec',   Land_boundary%lprec,   Time )
    call data_override_ug('LND', 'fprec',   Land_boundary%fprec,   Time )
    call data_override_ug('LND', 'dhdt',    Land_boundary%dhdt,    Time )
    call data_override_ug('LND', 'drdt',    Land_boundary%drdt,    Time )
    call data_override_ug('LND', 'p_surf',  Land_boundary%p_surf,  Time )
    do tr = 1,n_lnd_tr
       call get_tracer_names(MODEL_LAND, tr, tr_name)
       call data_override_ug('LND', trim(tr_name)//'_flux', Land_boundary%tr_flux(:,:,tr), Time)
       call data_override_ug('LND', 'dfd'//trim(tr_name),   Land_boundary%dfdtr  (:,:,tr), Time)
#else
    call data_override('LND', 't_flux',  Land_boundary%t_flux,  Time )
    call data_override('LND', 'lw_flux', Land_boundary%lw_flux, Time )
    call data_override('LND', 'sw_flux', Land_boundary%sw_flux, Time )
    call data_override('LND', 'sw_flux_down_vis_dir', Land_boundary%sw_flux_down_vis_dir, Time )
    call data_override('LND', 'sw_flux_down_total_dir', Land_boundary%sw_flux_down_total_dir, Time )
    call data_override('LND', 'sw_flux_down_vis_dif', Land_boundary%sw_flux_down_vis_dif, Time )
    call data_override('LND', 'sw_flux_down_total_dif', Land_boundary%sw_flux_down_total_dif, Time )

    call data_override('LND', 'lprec',   Land_boundary%lprec,   Time )
    call data_override('LND', 'fprec',   Land_boundary%fprec,   Time )
    call data_override('LND', 'dhdt',    Land_boundary%dhdt,    Time )
    call data_override('LND', 'drdt',    Land_boundary%drdt,    Time )
    call data_override('LND', 'p_surf',  Land_boundary%p_surf,  Time )
    do tr = 1,n_lnd_tr
       call get_tracer_names(MODEL_LAND, tr, tr_name)
       call data_override('LND', trim(tr_name)//'_flux', Land_boundary%tr_flux(:,:,:,tr), Time)
       call data_override('LND', 'dfd'//trim(tr_name),   Land_boundary%dfdtr  (:,:,:,tr), Time)
#endif
    enddo

    !-----------------------------------------------------------------------
    !---- output fields on the ice grid -------

    call get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
    call get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_tr(:,isphum), xmap_sfc)
    call get_from_xgrid (Ice_boundary%sw_flux_vis_dir,  'OCN', ex_flux_sw_vis_dir,   xmap_sfc)
    call get_from_xgrid (Ice_boundary%sw_flux_nir_dir,  'OCN', ex_flux_sw_dir,xmap_sfc)
    Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_flux_nir_dir - Ice_boundary%sw_flux_vis_dir ! ice & ocean use these 4: dir/dif nir/vis

    call get_from_xgrid (Ice_boundary%sw_flux_vis_dif,  'OCN', ex_flux_sw_vis_dif,   xmap_sfc)
    call get_from_xgrid (Ice_boundary%sw_flux_nir_dif,  'OCN', ex_flux_sw_dif,xmap_sfc)
    Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_flux_nir_dif - Ice_boundary%sw_flux_vis_dif ! ice & ocean use these 4: dir/dif nir/vis

    call get_from_xgrid (Ice_boundary%sw_down_vis_dir,  'OCN', ex_flux_sw_down_vis_dir,   xmap_sfc)
    call get_from_xgrid (Ice_boundary%sw_down_nir_dir,  'OCN', ex_flux_sw_down_total_dir, xmap_sfc)
    Ice_boundary%sw_down_nir_dir = Ice_boundary%sw_down_nir_dir - Ice_boundary%sw_down_vis_dir ! ice & ocean use these 4: dir/dif nir/vis

    call get_from_xgrid (Ice_boundary%sw_down_vis_dif,  'OCN', ex_flux_sw_down_vis_dif,   xmap_sfc)
    call get_from_xgrid (Ice_boundary%sw_down_nir_dif,  'OCN', ex_flux_sw_down_total_dif,xmap_sfc)
    Ice_boundary%sw_down_nir_dif = Ice_boundary%sw_down_nir_dif - Ice_boundary%sw_down_vis_dif ! ice & ocean use these 4: dir/dif nir/vis

    call get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
    call get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
    call get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
    call get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
    call get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
    call get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
    call get_from_xgrid (Ice_boundary%u_star,   'OCN', ex_u_star,    xmap_sfc)
    call get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)
    call get_from_xgrid (Ice_boundary%p,        'OCN', ex_slp,       xmap_sfc) ! mw mod

    call get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
    call get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
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
    do n = 1, Ice_boundary%fluxes%num_bcs  !{
      if(ex_gas_fluxes%bc(n)%flux_type  .ne. 'air_sea_deposition') then
       do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
          call get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
               ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)
       enddo  !} m
      endif
    enddo  !} n

    !Balaji: data_override calls moved here from coupler_main
    call data_override('ICE', 'u_flux', Ice_boundary%u_flux,  Time)
    call data_override('ICE', 'v_flux', Ice_boundary%v_flux,  Time)
    call data_override('ICE', 't_flux', Ice_boundary%t_flux,  Time)
    call data_override('ICE', 'q_flux', Ice_boundary%q_flux,  Time)
    call data_override('ICE', 'lw_flux',Ice_boundary%lw_flux, Time)
    call data_override('ICE', 'lw_flux_dn',Ice_boundary%lw_flux, Time, override=ov)
    if (ov) then
       Ice_boundary%lw_flux = Ice_boundary%lw_flux - stefan*Ice%t_surf**4
    endif
    call data_override('ICE', 'sw_flux_nir_dir',Ice_boundary%sw_flux_nir_dir, Time)
    call data_override('ICE', 'sw_flux_vis_dir',Ice_boundary%sw_flux_vis_dir, Time)
    call data_override('ICE', 'sw_flux_nir_dif',Ice_boundary%sw_flux_nir_dif, Time, override=ov)
    call data_override('ICE', 'sw_flux_vis_dif',Ice_boundary%sw_flux_vis_dif, Time)
    call data_override('ICE', 'sw_flux_vis_dir_dn',Ice_boundary%sw_down_vis_dir, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_vis_dir = Ice_boundary%sw_down_vis_dir*(1.0-Ice%albedo_vis_dir)
    endif
    call data_override('ICE', 'sw_flux_vis_dif_dn',Ice_boundary%sw_down_vis_dif, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_vis_dif = Ice_boundary%sw_down_vis_dif*(1.0-Ice%albedo_vis_dif)
    endif
    call data_override('ICE', 'sw_flux_nir_dir_dn',Ice_boundary%sw_down_nir_dir, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_nir_dir = Ice_boundary%sw_down_nir_dir*(1.0-Ice%albedo_nir_dir)
    endif
    call data_override('ICE', 'sw_flux_nir_dif_dn',Ice_boundary%sw_down_nir_dif, Time, override=ov)
    if (ov) then
       Ice_boundary%sw_flux_nir_dif = Ice_boundary%sw_down_nir_dif*(1.0-Ice%albedo_nir_dif)
    endif
    call data_override('ICE', 'lprec',  Ice_boundary%lprec,   Time)
    call data_override('ICE', 'fprec',  Ice_boundary%fprec,   Time)
    call data_override('ICE', 'dhdt',   Ice_boundary%dhdt,    Time)
    call data_override('ICE', 'dedt',   Ice_boundary%dedt,    Time)
    call data_override('ICE', 'drdt',   Ice_boundary%drdt,    Time)
    call data_override('ICE', 'coszen', Ice_boundary%coszen,  Time)
    call data_override('ICE', 'p',      Ice_boundary%p,       Time)

    call coupler_type_data_override('ICE', Ice_boundary%fluxes, Time)

    call coupler_type_send_data(Ice_boundary%fluxes, Time)

    ! compute stock changes

    ! Atm -> Lnd (precip)
#ifndef _USE_LEGACY_LAND_
    call stock_move_ug( &
         & FROM = Atm_stock(ISTOCK_WATER),  &
         & TO   = Lnd_stock(ISTOCK_WATER), &
         & DATA = (Land_boundary%lprec + Land_boundary%fprec), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Lnd) ')

    ! Atm -> Lnd (heat)
    call stock_move_ug( &
         & FROM = Atm_stock(ISTOCK_HEAT),  &
         & TO   = Lnd_stock(ISTOCK_HEAT), &
         & DATA = (-Land_boundary%t_flux + Land_boundary%lw_flux +  Land_boundary%sw_flux - Land_boundary%fprec*HLF), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Lnd) ')
#else
    call stock_move( &
         & FROM = Atm_stock(ISTOCK_WATER),  &
         & TO   = Lnd_stock(ISTOCK_WATER), &
         & DATA = (Land_boundary%lprec + Land_boundary%fprec), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Lnd) ')

    ! Atm -> Lnd (heat)
    call stock_move( &
         & FROM = Atm_stock(ISTOCK_HEAT),  &
         & TO   = Lnd_stock(ISTOCK_HEAT), &
         & DATA = (-Land_boundary%t_flux + Land_boundary%lw_flux +  Land_boundary%sw_flux - Land_boundary%fprec*HLF), &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Lnd) ')
#endif

    ! Atm -> Ice (precip)
    call stock_move( &
         & FROM = Atm_stock(ISTOCK_WATER), &
         & TO   = Ice_stock(ISTOCK_WATER), &
         & DATA = (Ice_boundary%lprec + Ice_boundary%fprec), &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move PRECIP (Atm->Ice) ')

    ! Atm -> Ice (heat)
    call stock_move( &
         & FROM = Atm_stock(ISTOCK_HEAT), &
         & TO   = Ice_stock(ISTOCK_HEAT), &
         & DATA = (-Ice_boundary%t_flux + Ice_boundary%lw_flux - Ice_boundary%fprec*HLF + Ice_boundary%sw_flux_vis_dir + &
         Ice_boundary%sw_flux_vis_dif + Ice_boundary%sw_flux_nir_dir + Ice_boundary%sw_flux_nir_dif), &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & from_side=ISTOCK_BOTTOM, to_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move HEAT (Atm->Ice) ')

    deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)

    !=======================================================================
    !-------------------- diagnostics section ------------------------------

    !------- zonal wind stress -----------
    used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
    used = send_data ( id_tauu,  -Atmos_boundary%u_flux, Time )

    !------- meridional wind stress -----------
    used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
    used = send_data ( id_tauv,  -Atmos_boundary%v_flux, Time )

    !Balaji
    call mpp_clock_end(fluxAtmDnClock)
    call mpp_clock_end(cplClock)
    !=======================================================================

  end subroutine flux_down_from_atmos

  !#######################################################################
  !> \brief Optimizes the exchange grids by eliminating land and ice partitions with no data.
  !!
  !! Optimizes the exchange grids by eliminating land and ice partitions with no data.
  subroutine generate_sfc_xgrid( Land, Ice )
    ! subroutine to regenerate exchange grid eliminating side 2 tiles with 0 frac area
    type(land_data_type), intent(in) :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),  intent(in) :: Ice !< A derived data type to specify ice boundary data

    integer :: isc, iec, jsc, jec

    !Balaji
    call mpp_clock_begin(cplClock)
    call mpp_clock_begin(regenClock)

    call mpp_get_compute_domain(Ice%Domain, isc, iec, jsc, jec)

    call set_frac_area (Ice%part_size(isc:iec,jsc:jec,:) , 'OCN', xmap_sfc)
#ifndef _USE_LEGACY_LAND_
    call set_frac_area_ug (Land%tile_size, 'LND', xmap_sfc)
#else
    call set_frac_area (Land%tile_size, 'LND', xmap_sfc)
#endif
    n_xgrid_sfc = max(xgrid_count(xmap_sfc),1)
    if(n_xgrid_sfc .GE. nblocks) then
       my_nblocks = nblocks
       call mpp_compute_extent(1, n_xgrid_sfc, nblocks, block_start, block_end)
    else
       my_nblocks = 1
       block_start(1) = 1
       block_end(1) = n_xgrid_sfc
    endif

    !Balaji
    call mpp_clock_end(regenClock)
    call mpp_clock_end(cplClock)
    return
  end subroutine generate_sfc_xgrid

  !#######################################################################
  !> \brief  Corrects the fluxes for consistency with the new surface temperatures in land
  !!         and ice models.
  !!
  !! Corrects the fluxes for consistency with the new surface temperatures in land
  !! and ice models. Final increments for temperature and specific humidity in the
  !! lowest atmospheric layer are computed and returned to the atmospheric model
  !! so that it can finalize the increments in the rest of the atmosphere.
  !!
  !! The following elements of the land_ice_atmos_boundary_type are computed:
  !! <pre>
  !!        dt_t  = temperature change at the lowest
  !!                 atmospheric level (deg k)
  !!        dt_q  = specific humidity change at the lowest
  !!                 atmospheric level (kg/kg)
  !! </pre>
  subroutine flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary, Land_boundary, Ice_boundary )
    type(time_type),      intent(in)    :: Time !< Current time
    type(land_data_type), intent(inout) :: Land !< A derived data type to specify ice boundary data
    type(ice_data_type),  intent(inout) :: Ice  !< A derived data type to specify ice boundary data
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary !< A derived data type to specify properties and fluxed
                                                                                 !! passed from exchange grid to the atmosphere, land and ice
    type(atmos_land_boundary_type), intent(inout)     :: Land_boundary
    type(atmos_ice_boundary_type),  intent(inout)     :: Ice_boundary

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
         ex_tr_surf_new,    & ! updated tracer values at the surface
         ex_dt_tr_surf,     & ! tendency of tracers at the surface
         ex_delta_tr_n
    ! jgj: added for co2_surf diagnostic
    real, dimension(n_xgrid_sfc) :: &
         ex_co2_surf_dvmr   ! updated CO2 tracer values at the surface (dry vmr)

    real, dimension(size(Land_Ice_Atmos_Boundary%dt_t,1),size(Land_Ice_Atmos_Boundary%dt_t,2)) :: diag_atm, &
         evap_atm, frac_atm
#ifndef _USE_LEGACY_LAND_
    real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2)) :: data_lnd, diag_land
#else
    real, dimension(size(Land_boundary%lprec,1), size(Land_boundary%lprec,2), size(Land_boundary%lprec,3)) :: data_lnd, diag_land
#endif
    real, dimension(size(Ice_boundary%lprec,1), size(Ice_boundary%lprec,2), size(Ice_boundary%lprec,3)) :: data_ice
    real, dimension(size(Ice%albedo,1),size(Ice%albedo,2),size(Ice%albedo,3)) ::  icegrid
    logical :: used

    integer :: tr       ! tracer index
    character(32) :: tr_name, tr_units ! tracer name
    integer :: n, i, m, ier

    integer :: is, ie, l

    !Balaji
    call mpp_clock_begin(cplClock)
    call mpp_clock_begin(fluxAtmUpClock)
    !-----------------------------------------------------------------------
    !Balaji: data_override calls moved here from coupler_main
    call data_override ( 'ICE', 't_surf', Ice%t_surf,  Time)
#ifndef _USE_LEGACY_LAND_
    call data_override_ug ( 'LND', 't_ca',   Land%t_ca,   Time)
    call data_override_ug ( 'LND', 't_surf', Land%t_surf, Time)
#else
    call data_override ( 'LND', 't_ca',   Land%t_ca,   Time)
    call data_override ( 'LND', 't_surf', Land%t_surf, Time)
#endif
    do tr = 1, n_lnd_tr
       call get_tracer_names( MODEL_LAND, tr, tr_name )
#ifndef _USE_LEGACY_LAND_
       call data_override_ug('LND', trim(tr_name)//'_surf', Land%tr(:,:,tr), Time)
#else
       call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
#endif
    enddo

    !----- compute surface temperature change -----

    ex_t_surf_new = 200.0

    call put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
    ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
#ifndef _USE_LEGACY_LAND_
    call put_to_xgrid_ug (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
    call put_to_xgrid_ug (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)
#else
    call put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
    call put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)
#endif
    !  call escomp(ex_t_ca_new, ex_q_surf_new)
    !  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new)
    !  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

#ifdef SCM
    if (do_specified_flux .and. do_specified_land) then
       ex_t_surf_new = ex_t_surf
       ex_t_ca_new   = ex_t_ca
    endif
#endif

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

       !-----------------------------------------------------------------------
       !-----  adjust fluxes and atmospheric increments for
       !-----  implicit dependence on surface temperature -----
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
    enddo !  l = 1, my_nblocks
    ! get all tracers available from land, and calculate changes in near-tracer field
    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
#ifndef _USE_LEGACY_LAND_
          call put_to_xgrid_ug ( Land%tr(:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
#else
          call put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
#endif
       endif
    enddo

    ! get all tracers available from ocean here

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
       call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_tr(:,:,n), 'ATM', ex_delta_tr_n(:,tr), xmap_sfc)
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

    !-----------------------------------------------------------------------
    !---- get mean quantites on atmospheric grid ----

    call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)
#ifndef use_AM3_physics
    call get_from_xgrid (Land_Ice_Atmos_Boundary%shflx,'ATM', ex_flux_t    , xmap_sfc) !miz
    call get_from_xgrid (Land_Ice_Atmos_Boundary%lhflx,'ATM', ex_flux_tr(:,isphum), xmap_sfc)!miz
#endif

    !=======================================================================
    !-------------------- diagnostics section ------------------------------

    !------- new surface temperature -----------
#ifdef use_AM3_physics
    if ( id_t_surf > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
       used = send_data ( id_t_surf, diag_atm, Time )
    endif
#else
    call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
    if ( id_t_surf > 0 ) then
       used = send_data ( id_t_surf, diag_atm, Time )
    endif
    if ( id_ts > 0 ) then
       used = send_data ( id_ts, diag_atm, Time )
    endif
    call sum_diag_integral_field ('t_surf', diag_atm)
#ifndef use_AM3_physics
    if ( id_ts_g > 0 ) used = send_global_diag ( id_ts_g, diag_atm, Time )
#endif
    !------- new surface temperature only over open ocean -----------
    if ( id_tos > 0 ) then
       ex_icetemp = 0.0
       icegrid = 0.0; icegrid(:,:,1) = 1.0
       call put_to_xgrid ( icegrid, 'OCN', ex_icetemp, xmap_sfc)
       ex_temp = ex_t_surf_new * ex_icetemp
       call get_from_xgrid (diag_atm, 'ATM', ex_temp, xmap_sfc)
       call get_from_xgrid (frac_atm, 'ATM', ex_icetemp, xmap_sfc)
       where (frac_atm > 0.0)
          diag_atm = (diag_atm/frac_atm) ! - tfreeze  CMIP6 in degK
          frac_atm = 1.0
       elsewhere
          diag_atm = 0.0
          frac_atm = 0.0
       endwhere
       used = send_data ( id_tos, diag_atm, Time, rmask=frac_atm )
    endif

    !------- new surface temperature only over land and sea-ice -----------
    if ( id_tslsi > 0 ) then
       ex_land_frac = 0.0
       call put_logical_to_real (Land%mask, 'LND', ex_land_frac, xmap_sfc)
       icegrid = 1.0; icegrid(:,:,1) = 0.
       ex_icetemp = 0.
       call put_to_xgrid (icegrid, 'OCN', ex_icetemp, xmap_sfc)
       ex_icetemp = ex_icetemp + ex_land_frac
       ex_temp = ex_t_surf_new * ex_icetemp
       call get_from_xgrid (diag_atm, 'ATM', ex_temp, xmap_sfc)
       call get_from_xgrid (frac_atm, 'ATM', ex_icetemp, xmap_sfc)
       where (frac_atm > 0.0)
          diag_atm = diag_atm/frac_atm
          frac_atm = 1.0
       elsewhere
          diag_atm = 0.0
          frac_atm = 0.0
       endwhere
       used = send_data ( id_tslsi, diag_atm, Time, rmask=frac_atm )
    endif
#endif


    ! + slm, Mar 27 2002
    ! ------ new canopy temperature --------
    !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
    !   but this will be changed in future version of the land madel
    if ( id_t_ca > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_t_ca_new, xmap_sfc)
       used = send_data ( id_t_ca, diag_atm, Time )
    endif

    !------- updated surface tracer fields ------
    do tr=1,n_exch_tr
       call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
       if ( id_tr_surf(tr) > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_tr_surf_new(:,tr), xmap_sfc)
          used = send_data ( id_tr_surf(tr), diag_atm, Time )
       endif
       !!jgj:  add dryvmr co2_surf
       ! - slm Mar 25, 2010: moved to resolve interdependence of diagnostic fields
       if ( id_co2_surf_dvmr > 0 .and. lowercase(trim(tr_name))=='co2') then
          ex_co2_surf_dvmr = (ex_tr_surf_new(:,tr) / (1.0 - ex_tr_surf_new(:,isphum))) * WTMAIR/WTMCO2
          call get_from_xgrid (diag_atm, 'ATM', ex_co2_surf_dvmr, xmap_sfc)
          used = send_data ( id_co2_surf_dvmr, diag_atm, Time )
       endif
    enddo

    !------- sensible heat flux -----------
    if ( id_t_flux > 0 .or. id_hfss > 0 .or. id_hfss_g > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
       if ( id_t_flux > 0 ) used = send_data ( id_t_flux, diag_atm, Time )
       if ( id_hfss   > 0 ) used = send_data ( id_hfss, diag_atm, Time )
#ifndef use_AM3_physics
       if ( id_hfss_g > 0 ) used = send_global_diag ( id_hfss_g, diag_atm, Time )
#endif
    endif

    !------- net longwave flux -----------
    if ( id_r_flux > 0 .or. id_rls_g > 0 ) then
       call get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
       if ( id_r_flux > 0 ) used = send_data ( id_r_flux, diag_atm, Time )
#ifndef use_AM3_physics
       if ( id_rls_g  > 0 ) used = send_global_diag ( id_rls_g, diag_atm, Time )
#endif
    endif

    !------- tracer fluxes ------------
    ! tr_mol_flux diagnostic will be correct for co2 tracer only.
    ! will need update code to use correct molar mass for tracers other than co2
    do tr=1,n_exch_tr
       if ( id_tr_flux(tr) > 0 .or. id_tr_mol_flux(tr) > 0 ) then
          call get_from_xgrid (diag_atm, 'ATM', ex_flux_tr(:,tr), xmap_sfc)
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name, units=tr_units )
          if (id_tr_flux(tr) > 0 ) &
               used = send_data ( id_tr_flux(tr), diag_atm, Time )
    !     if (id_tr_mol_flux(tr) > 0 ) &
    !          used = send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMCO2, Time)
    ! 2017/08/08 jgj - replaced 2 lines above by the following
          if (id_tr_mol_flux(tr) > 0 .and. lowercase(trim(tr_name))=='co2') then
               used = send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMCO2, Time)
    !sometimes in 2018 f1p for vmr tracers
           elseif (id_tr_mol_flux(tr) > 0 .and. lowercase(trim(tr_units)).eq."vmr") then
              call get_from_xgrid (diag_atm, 'ATM', ex_flux_tr(:,tr)*(1.-ex_tr_surf_new(:,isphum)), xmap_sfc)
              used = send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMAIR, Time)
       endif
       endif
    enddo

#ifndef _USE_LEGACY_LAND_
    if ( id_t_flux_land > 0 ) then
       call get_from_xgrid_ug (diag_land, 'LND', ex_flux_t, xmap_sfc)
       call send_tile_data ( id_t_flux_land, diag_land )
    endif
    !------- tracer fluxes for land
    do tr=1,n_exch_tr
       if ( id_tr_flux_land(tr) > 0 .or. id_tr_mol_flux_land(tr) > 0 ) then
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name, units=tr_units )
          call get_from_xgrid_ug (diag_land, 'LND', ex_flux_tr(:,tr), xmap_sfc)
          if (id_tr_flux_land(tr) > 0 ) &
                 call send_tile_data (id_tr_flux_land(tr), diag_land )
          if (id_tr_mol_flux_land(tr) > 0) then
             if (lowercase(trim(tr_name))=='co2') then
                call send_tile_data (id_tr_mol_flux_land(tr), diag_land*1000./WTMCO2)
             elseif (lowercase(trim(tr_units)).eq.'vmr') then
                call get_from_xgrid_ug (diag_land, 'LND', ex_flux_tr(:,tr)*(1.-ex_tr_surf_new(:,isphum)), xmap_sfc)
                call send_tile_data (id_tr_mol_flux_land(tr), diag_atm*1000./WTMAIR )
             endif
          endif
       endif
    enddo
#endif

    !-----------------------------------------------------------------------
    !---- accumulate global integral of evaporation (mm/day) -----
    call get_from_xgrid (evap_atm, 'ATM', ex_flux_tr(:,isphum), xmap_sfc)
    if( id_q_flux > 0 )  used = send_data ( id_q_flux, evap_atm, Time)
    if( id_evspsbl > 0 ) used = send_data ( id_evspsbl, evap_atm, Time)
    if( id_hfls    > 0 ) used = send_data ( id_hfls, HLV*evap_atm, Time)
#ifndef use_AM3_physics
    if( id_hfls_g  > 0 ) used = send_global_diag ( id_hfls_g, HLV*evap_atm, Time)
#endif

#ifndef _USE_LEGACY_LAND_
    if( id_q_flux_land > 0 ) then
       call get_from_xgrid_ug (diag_land, 'LND', ex_flux_tr(:,isphum), xmap_sfc)
       call send_tile_data (id_q_flux_land, diag_land)
#else
    if( id_q_flux_land > 0 ) then
       call get_from_xgrid (diag_land, 'LND', ex_flux_tr(:,isphum), xmap_sfc)
       used = send_tile_averaged_data(id_q_flux_land, diag_land, &
            Land%tile_size, Time, mask=Land%mask)
#endif
    endif
    call sum_diag_integral_field ('evap', evap_atm*86400.)
#ifndef use_AM3_physics
    if (id_evspsbl_g > 0) used = send_global_diag ( id_evspsbl_g, evap_atm, Time )
#endif

#ifndef _USE_LEGACY_LAND_
    call send_tile_data (id_q_flux_land, diag_land)
    ! need this to avoid diag issues with tiling changes in update_land_slow
    call dump_tile_diag_fields(Time)
    call get_from_xgrid_ug(data_lnd, 'LND', ex_flux_tr(:,isphum), xmap_sfc)

    ! compute stock changes

    ! Lnd -> Atm (evap)
    call stock_move_ug( &
         & TO   = Atm_stock(ISTOCK_WATER), &
         & FROM = Lnd_stock(ISTOCK_WATER), &
         & DATA = data_lnd, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Lnd->ATm) ')

    ! Lnd -> Atm (heat lost through evap)
    call stock_move_ug( &
         & TO   = Atm_stock(ISTOCK_HEAT), &
         & FROM = Lnd_stock(ISTOCK_HEAT), &
         & DATA = data_lnd * HLV, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Lnd->ATm) ')
#else
    call get_from_xgrid(data_lnd, 'LND', ex_flux_tr(:,isphum), xmap_sfc)

    ! compute stock changes

    ! Lnd -> Atm (evap)
    call stock_move( &
         & TO   = Atm_stock(ISTOCK_WATER), &
         & FROM = Lnd_stock(ISTOCK_WATER), &
         & DATA = data_lnd, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Lnd->ATm) ')

    ! Lnd -> Atm (heat lost through evap)
    call stock_move( &
         & TO   = Atm_stock(ISTOCK_HEAT), &
         & FROM = Lnd_stock(ISTOCK_HEAT), &
         & DATA = data_lnd * HLV, &
         & grid_index=X1_GRID_LND, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_SIDE, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Lnd->ATm) ')
#endif

    call get_from_xgrid(data_ice, 'OCN', ex_flux_tr(:,isphum), xmap_sfc)

    ! Ice -> Atm (evap)
    call stock_move( &
         & TO   = Atm_stock(ISTOCK_WATER), &
         & FROM = Ice_stock(ISTOCK_WATER), &
         & DATA = data_ice, &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP (Ice->ATm) ')

    ! Ice -> Atm (heat lost through evap)
    call stock_move( &
         & TO   = Atm_stock(ISTOCK_HEAT), &
         & FROM = Ice_stock(ISTOCK_HEAT), &
         & DATA = data_ice * HLV, &
         & grid_index=X1_GRID_ICE, &
         & xmap=xmap_sfc, &
         & delta_t=Dt_atm, &
         & to_side=ISTOCK_TOP, from_side=ISTOCK_TOP, &
         & radius=Radius, ier=ier, verbose='stock move EVAP*HLV (Ice->ATm) ')

    !Balaji
    call mpp_clock_end(fluxAtmUpClock)
    call mpp_clock_end(cplClock)
  end subroutine flux_up_to_atmos

  subroutine flux_ex_arrays_dealloc
    integer :: m,n

    !=======================================================================
    !---- deallocate module storage ----
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
         ex_seawater ,  &
         ex_land        )

#ifdef SCM
    deallocate ( &
         ex_dhdt_surf_forland, &
         ex_dedt_surf_forland, &
         ex_dedq_surf_forland  )
#endif

    ! Extra fluxes
    do n = 1, ex_gas_fields_ice%num_bcs  !{
       do m = 1, ex_gas_fields_ice%bc(n)%num_fields  !{
          deallocate ( ex_gas_fields_ice%bc(n)%field(m)%values )
          nullify ( ex_gas_fields_ice%bc(n)%field(m)%values )
       enddo  !} m
    enddo  !} n

    do n = 1, ex_gas_fields_atm%num_bcs  !{
       do m = 1, ex_gas_fields_atm%bc(n)%num_fields  !{
          deallocate ( ex_gas_fields_atm%bc(n)%field(m)%values )
          nullify ( ex_gas_fields_atm%bc(n)%field(m)%values )
       enddo  !} m
    enddo  !} n

    do n = 1, ex_gas_fluxes%num_bcs  !{
       do m = 1, ex_gas_fluxes%bc(n)%num_fields  !{
          deallocate ( ex_gas_fluxes%bc(n)%field(m)%values )
          nullify ( ex_gas_fluxes%bc(n)%field(m)%values )
       enddo  !} m
    enddo  !} n

  end subroutine flux_ex_arrays_dealloc

  subroutine flux_atmos_to_ocean(Time, Atm, Ice_boundary, Ice)
  type(time_type),               intent(in)   :: Time         !< Current time
  type(atmos_data_type),         intent(inout):: Atm          !< A derived data type to specify atmosphere boundary data
  type(atmos_ice_boundary_type), intent(inout):: Ice_boundary !< A derived data type to specify properties and fluxes passed
                                                                  !! from atmosphere to ice
  type(ice_data_type),           intent(inout):: Ice

  integer :: n,m
  logical :: used

#ifndef use_AM3_physics
  call atmos_tracer_driver_gather_data_down(Atm%fields, Atm%tr_bot)
#endif

  !air-sea deposition fluxes
  do n = 1, Atm%fields%num_bcs  !{
   !Do the string copies.
   Atm%fields%bc(n)%flux_type = trim(ex_gas_fluxes%bc(n)%flux_type)
   Atm%fields%bc(n)%implementation = trim(ex_gas_fluxes%bc(n)%implementation)
   if(ex_gas_fields_atm%bc(n)%flux_type  .eq. 'air_sea_deposition') then
    do m = 1, Atm%fields%bc(n)%num_fields  !{
      call put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
           ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method)
    enddo  !} m
   endif
  enddo  !} n

  ! Calculate ocean explicit flux here

  call atmos_ocean_dep_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater)

  do n = 1, Ice_boundary%fluxes%num_bcs  !{
     if(Ice_boundary%fluxes%bc(n)%flux_type  .eq. 'air_sea_deposition') then
        do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
           call get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
                ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)

           call data_override('ICE', Ice_boundary%fluxes%bc(n)%field(m)%name,     &
              Ice_boundary%fluxes%bc(n)%field(m)%values, Time)
           if ( Ice_boundary%fluxes%bc(n)%field(m)%id_diag > 0 ) then  !{
              used = send_data(Ice_boundary%fluxes%bc(n)%field(m)%id_diag, Ice_boundary%fluxes%bc(n)%field(m)%values, Time )
           endif  !}
        enddo  !} m
     endif
  enddo  !} n

  call update_ice_atm_deposition_flux( Ice_boundary, Ice )

  end subroutine flux_atmos_to_ocean

  !#######################################################################

  !> \brief Puts land or ice model masks (with partitions) onto the
  !! exchange grid as a real array (1.=true, 0.=false)
  subroutine put_logical_to_real_sg (mask, id, ex_mask, xmap)

    logical         , intent(in)    :: mask(:,:,:)
    character(len=3), intent(in)    :: id
    real            , intent(inout) :: ex_mask(:)
    type(xmap_type), intent(inout) :: xmap

    !-----------------------------------------------------------------------
    !    puts land or ice model masks (with partitions) onto the
    !    exchange grid as a real array (1.=true, 0.=false)
    !-----------------------------------------------------------------------

    real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask

    where (mask)
       rmask = 1.0
    elsewhere
       rmask = 0.0
    endwhere

    call put_to_xgrid(rmask, id, ex_mask, xmap)

  end subroutine put_logical_to_real_sg

  !#######################################################################

  !> \brief Puts land or ice model masks (with partitions) onto the
  !! exchange grid as a real array (1.=true, 0.=false)
  subroutine put_logical_to_real_ug (mask, id, ex_mask, xmap)

    logical         , intent(in)    :: mask(:,:)
    character(len=3), intent(in)    :: id
    real            , intent(inout) :: ex_mask(:)
    type(xmap_type), intent(inout) :: xmap

    !-----------------------------------------------------------------------
    !    puts land or ice model masks (with partitions) onto the
    !    exchange grid as a real array (1.=true, 0.=false)
    !-----------------------------------------------------------------------

    real, dimension(size(mask,1),size(mask,2)) :: rmask

    where (mask)
       rmask = 1.0
    elsewhere
       rmask = 0.0
    endwhere

#ifndef _USE_LEGACY_LAND_
    call put_to_xgrid_ug(rmask, id, ex_mask, xmap)
#else
    call put_to_xgrid (rmask, id, ex_mask, xmap)
#endif

  end subroutine put_logical_to_real_ug


  !#######################################################################

  !> \brief Initializes diagnostic fields that may be output from this
  !! module (the ID numbers may be referenced anywhere in this module)
  subroutine diag_field_init ( Time, atmos_axes, land_axes, land_pe )

    type(time_type), intent(in) :: Time
    integer,         intent(in) :: atmos_axes(2)
    integer,         intent(in) :: land_axes(:)
    logical,         intent(in) :: land_pe

    integer :: iref
    character(len=6) :: label_zm, label_zh
    real, dimension(2) :: trange = (/  100., 400. /), &
         vrange = (/ -400., 400. /), &
         frange = (/ -0.01, 1.01 /)
    character(len=32)  :: name, units ! name of the tracer
    character(len=128) :: longname    ! long name of the tracer
    integer            :: tr          ! tracer index
    integer            :: area_id
    !-----------------------------------------------------------------------
    !  initializes diagnostic fields that may be output from this module
    !  (the id numbers may be referenced anywhere in this module)
    !-----------------------------------------------------------------------

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

    !--------- initialize static diagnostic fields --------------------

    id_land_mask = &
         register_static_field ( mod_name, 'land_mask', atmos_axes,  &
         'fractional amount of land', 'none', &
         range=frange, interp_method = "conserve_order1" )

    !--------- initialize diagnostic fields --------------------

    id_ice_mask = &
         register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
         'fractional amount of sea ice', 'none',  &
         range=frange, interp_method = "conserve_order1" )

    id_wind = &
         register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
         'wind speed for flux calculations', 'm/s', &
         range=(/0.,vrange(2)/) )

    id_drag_moist = &
         register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
         'drag coeff for moisture',    'none'     )

    id_drag_heat  = &
         register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
         'drag coeff for heat',    'none'     )

    id_drag_mom   = &
         register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
         'drag coeff for momentum',     'none'     )

    id_rough_moist = &
         register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
         'surface roughness for moisture',  'm'  )

    id_rough_heat = &
         register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
         'surface roughness for heat',  'm'  )

    id_rough_mom  = &
         register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
         'surface roughness for momentum',  'm'  )

    id_u_star     = &
         register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
         'friction velocity',   'm/s'   )

    id_b_star     = &
         register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
         'buoyancy scale',      'm/s2'   )

    id_q_star     = &
         register_diag_field ( mod_name, 'q_star',     atmos_axes, Time, &
         'moisture scale',      'kg water/kg air'   )

    id_u_flux     = &
         register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
         'zonal wind stress',     'pa'   )

    id_v_flux     = &
         register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
         'meridional wind stress',     'pa'   )

    id_t_surf     = &
         register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
         'surface temperature',    'deg_k', &
         range=trange    )

    ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
    id_t_ca       = &
         register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
         'canopy air temperature',    'deg_k', &
         range=trange    )

    ! - slm, Mar 25, 2002
    id_z_atm      = &
         register_diag_field ( mod_name, 'z_atm',     atmos_axes, Time, &
         'height of btm level',    'm')

    id_p_atm      = &
         register_diag_field ( mod_name, 'p_atm',     atmos_axes, Time, &
         'pressure at btm level',    'pa')

    ! - bw, Mar 25, 2002 -- added diagnostic slp
    id_slp      = &
         register_diag_field ( mod_name, 'slp',      atmos_axes, Time, &
         'sea level pressure',    'pa')

    id_gust       = &
         register_diag_field ( mod_name, 'gust',     atmos_axes, Time, &
         'gust scale',    'm/s')

    id_t_flux     = &
         register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
         'sensible heat flux',     'w/m2'    )

    id_r_flux     = &
         register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
         'net (down-up) longwave flux',   'w/m2'    )

    id_t_atm      = &
         register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
         'temperature at btm level',    'deg_k', &
         range=trange     )

    id_u_atm      = &
         register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
         'u wind component at btm level',  'm/s', &
         range=vrange    )

    id_v_atm      = &
         register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
         'v wind component at btm level',  'm/s', &
         range=vrange    )

    id_t_ref      = &
         register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
         'temperature at '//label_zh, 'deg_k' , &
         range=trange      )

    id_rh_ref     = &
         register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
         'relative humidity at '//label_zh, 'percent' )

    id_rh_ref_cmip = &
         register_diag_field ( mod_name, 'rh_ref_cmip',     atmos_axes, Time,   &
         'relative humidity at '//label_zh, 'percent' )

    id_u_ref      = &
         register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
         'zonal wind component at '//label_zm,  'm/s', &
         range=vrange )

    id_v_ref      = &
         register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
         'meridional wind component at '//label_zm, 'm/s', &
         range=vrange )

    id_wind_ref = &
         register_diag_field ( mod_name, 'wind_ref',   atmos_axes, Time,     &
         'absolute value of wind at '//label_zm, 'm/s', &
         range=vrange )

    id_del_h      = &
         register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
         'ref height interp factor for heat', 'none' )
    id_del_m      = &
         register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
         'ref height interp factor for momentum','none' )
    id_del_q      = &
         register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
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
       do tr = 1, n_exch_tr
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
          id_tr_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_flux', Land_axes, Time, &
               'flux of '//trim(longname), trim(units)//' kg air/(m2 s)', missing_value=-1.0 )
          if ( lowercase(trim(name))=='co2') then
             id_tr_mol_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, Time, &
                  'flux of '//trim(longname), 'mol CO2/(m2 s)', missing_value=-1.0 )
          else
             id_tr_mol_flux_land(tr) = register_tiled_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, Time, &
                  'flux of '//trim(longname), 'mol/(m2 s)', missing_value=-1.0 )
          endif
       enddo
#else
       id_t_ref_land = &
            register_diag_field ( 'flux_land', 't_ref', Land_axes, Time, &
            'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
            range=trange, missing_value =  -100.0)
       id_q_ref_land = &
            register_diag_field ( 'flux_land', 'q_ref', Land_axes, Time, &
            'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
            missing_value=-1.0)
       id_rh_ref_land= &
            register_diag_field ( 'flux_land', 'rh_ref', Land_axes, Time,   &
            'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
            missing_value=-999.0)
       id_u_ref_land = &
            register_diag_field ( 'flux_land', 'u_ref',  Land_axes, Time, &
            'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
            range=vrange, missing_value=-999.0 )
       id_v_ref_land = &
            register_diag_field ( 'flux_land', 'v_ref',  Land_axes, Time,     &
            'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
            range=vrange, missing_value = -999.0 )
       id_q_flux_land = &
            register_diag_field( 'flux_land', 'evap', Land_axes, Time, &
            'evaporation rate over land', 'kg/m2/s', missing_value=-1.0 )
       id_t_flux_land = &
            register_diag_field( 'flux_land', 'shflx', Land_axes, Time, &
            'sensible heat flux', 'W/m2', missing_value=-1.0 )
       id_tasLut_land = &
            register_diag_field( 'cmor_land', 'tasLut', Land_axes, Time, &
            'Near-Surface Air Temperature ('//trim(label_zh)//' Above Displacement Height) on Land Use Tile', &
            units='K', standard_name='air_temperature', missing_value=-1.0 )
       id_hussLut_land = &
            register_diag_field( 'cmor_land', 'hussLut', Land_axes, Time, &
            'Near-Surface Specific Humidity on Land Use Tile', '1.0', &
            standard_name='specific_humidity', missing_value=-1.0 )
       allocate(id_tr_flux_land(n_exch_tr))
       allocate(id_tr_mol_flux_land(n_exch_tr))
       do tr = 1, n_exch_tr
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
          id_tr_flux_land(tr) = register_diag_field( 'flux_land', trim(name)//'_flux', Land_axes, Time, &
               'flux of '//trim(longname), trim(units)//' kg air/(m2 s)', missing_value=-1.0 )
          if ( lowercase(trim(name))=='co2') then
             id_tr_mol_flux_land(tr) = register_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, Time, &
                  'flux of '//trim(longname), 'mol CO2/(m2 s)', missing_value=-1.0 )
          else
             id_tr_mol_flux_land(tr) = register_diag_field( 'flux_land', trim(name)//'_mol_flux', Land_axes, Time, &
                  'flux of '//trim(longname), 'mol/(m2 s)', missing_value=-1.0 )
          endif
       enddo
#endif
    endif

    id_q_ref = &
         register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
         'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)

    id_rough_scale = &
         register_diag_field ( mod_name, 'rough_scale', atmos_axes, Time, &
         'topographic scaling factor for momentum drag','1' )
    !-----------------------------------------------------------------------

    allocate(id_tr_atm(n_exch_tr))
    allocate(id_tr_surf(n_exch_tr))
    allocate(id_tr_flux(n_exch_tr))
    allocate(id_tr_mol_flux(n_exch_tr))
    allocate(id_tr_mol_flux0(n_exch_tr))

    do tr = 1, n_exch_tr
       call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
       id_tr_atm(tr) = register_diag_field (mod_name, trim(name)//'_atm', atmos_axes, Time, &
            trim(longname)//' at btm level', trim(units))
       id_tr_surf(tr) = register_diag_field (mod_name, trim(name)//'_surf', atmos_axes, Time, &
            trim(longname)//' at the surface', trim(units))
       id_tr_flux(tr) = register_diag_field(mod_name, trim(name)//'_flux', atmos_axes, Time, &
            'flux of '//trim(longname), trim(units)//' kg air/(m2 s)')
       !! add dryvmr co2_surf and co2_atm
       if ( lowercase(trim(name))=='co2') then
          ! - slm Mar 25, 2010: moved registration of mol_flux inside 'if' to disable
          ! saving incorrect results (mol fluxes for other tracers computed with CO2 molar
          ! mass)
          id_tr_mol_flux(tr) = register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
               'flux of '//trim(longname), 'mol CO2/(m2 s)')
          id_co2_atm_dvmr = register_diag_field (mod_name, trim(name)//'_atm_dvmr', atmos_axes, Time, &
               trim(longname)//' at btm level', 'mol CO2 /mol air')
          id_co2_surf_dvmr = register_diag_field (mod_name, trim(name)//'_surf_dvmr', atmos_axes, Time, &
               trim(longname)//' at the surface', 'mol CO2 /mol air')
       else
!f1p
          id_tr_mol_flux(tr) = register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
               'flux of '//trim(longname), 'mol/(m2 s)')
       endif
!f1p
       id_tr_mol_flux0(tr) = register_diag_field(mod_name, trim(name)//'_mol_flux_atm0', atmos_axes, Time, &
            'gross flux of '//trim(longname), 'mol/(m2 s)')

    enddo

    ! 2017/08/08 jgj add diagnostics for co2 data overrides even if co2 is not a tracer
    ! register data calls not needed here for co2_flux_pcair_atm and o2_flux_pcair_atm as this happens elsewhere
    id_co2_bot = register_diag_field (mod_name, 'co2_bot', atmos_axes, Time, &
           'co2_bot from data_override', 'ppmv')

    ! id_nh3_flux_atm0 = register_diag_field (mod_name, 'nh3_flux_atm0', atmos_axes, Time, &
    !        'nh3 flux out of the ocean assuming not nh3 in the atmosphere', 'mol/m2/s')


    id_q_flux = register_diag_field( mod_name, 'evap',       atmos_axes, Time, &
         'evaporation rate',        'kg/m2/s'  )

    !--------------------------------------------------------------------
    !    retrieve the diag_manager id for the area diagnostic,
    !    needed for cmorizing various diagnostics.
    !--------------------------------------------------------------------
    area_id = get_diag_field_id ('dynamics', 'area')
    if (area_id .eq. DIAG_FIELD_NOT_FOUND) call error_mesg &
         ('diag_field_init in atm_land_ice_flux_exchange_mod', &
         'diagnostic field "dynamics", "area" is not in the diag_table', NOTE)

    !-----------------------------------------------------------------------
    !  register cmip variable names
    !-----------------------------------------------------------------------
    ! NOTE: add extra dimension reference level fields?  height2m, height10m
    !       for now we will handle this with an attribute

    id_height2m = &
        register_static_field ( mod_name, 'height2m', (/null_axis_id/), &
                             'Height', 'm', standard_name = 'height' )
    if ( id_height2m > 0 ) then
       call diag_field_add_attribute( id_height2m, 'axis', 'Z' )
       call diag_field_add_attribute( id_height2m, 'positive', 'up' )
    endif

    id_height10m = &
        register_static_field ( mod_name, 'height10m', (/null_axis_id/), &
                             'Height', 'm', standard_name = 'height' )
    if ( id_height10m > 0 ) then
       call diag_field_add_attribute( id_height10m, 'axis', 'Z' )
       call diag_field_add_attribute( id_height10m, 'positive', 'up' )
    endif

#ifdef use_AM3_physics
    id_tas      = &
         register_diag_field ( mod_name, 'tas', atmos_axes, Time, &
         'Near-Surface Air Temperature', 'K' , &
         standard_name = 'air_temperature', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=trange )
    if ( id_tas > 0 .and. id_height2m > 0) &
       call diag_field_add_attribute( id_tas, 'coordinates', 'height2m' )

    id_uas      = &
         register_diag_field ( mod_name, 'uas', atmos_axes, Time, &
         'Eastward Near-Surface Wind', 'm s-1', &
         standard_name = 'eastward_wind', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_uas > 0 .and. id_height10m > 0) &
       call diag_field_add_attribute( id_uas, 'coordinates', 'height10m' )

    id_vas      = &
         register_diag_field ( mod_name, 'vas', atmos_axes, Time, &
         'Northward Near-Surface Wind', 'm s-1', &
         standard_name = 'northward_wind', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_vas > 0 .and. id_height10m > 0 ) &
       call diag_field_add_attribute( id_vas, 'coordinates', 'height10m' )

    id_sfcWind = &
         register_diag_field ( mod_name, 'sfcWind', atmos_axes, Time, &
         'Near-Surface Wind Speed', 'm s-1', &
         standard_name = 'wind_speed', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=vrange )
    if ( id_sfcWind > 0 .and. id_height10m > 0 ) &
       call diag_field_add_attribute( id_sfcWind, 'coordinates', 'height10m' )

    id_huss = &
         register_diag_field ( mod_name, 'huss', atmos_axes, Time, &
         'Near-Surface Specific Humidity', '1.0', &
         standard_name = 'specific_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_huss > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_huss, 'coordinates', 'height2m' )

    id_hurs = &
         register_diag_field ( mod_name, 'hurs', atmos_axes, Time, &
         'Near-Surface Relative Humidity', '%', &
         standard_name = 'relative_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_hurs > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_hurs, 'coordinates', 'height2m' )

    id_rhs = &
         register_diag_field ( mod_name, 'rhs', atmos_axes, Time, &
         'Near-Surface Relative Humidity', '%', &
         standard_name = 'relative_humidity', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_rhs > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_rhs, 'coordinates', 'height2m' )

    id_ts = &
         register_diag_field ( mod_name, 'ts', atmos_axes, Time, &
         'Surface Temperature', 'K', &
         standard_name = 'surface_temperature', area=area_id, &
         missing_value=CMOR_MISSING_VALUE, range=trange )

    id_psl = &
         register_diag_field ( mod_name, 'psl', atmos_axes, Time, &
         'Sea Level Pressure', 'Pa', &
         standard_name = 'air_pressure_at_sea_level', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_tauu = &
         register_diag_field ( mod_name, 'tauu', atmos_axes, Time, &
         'Surface Downward Eastward Wind Stress', 'Pa', &
         standard_name = 'surface_downward_eastward_stress', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_tauv = &
         register_diag_field ( mod_name, 'tauv', atmos_axes, Time, &
         'Surface Downward Northward Wind Stress', 'Pa', &
         standard_name = 'surface_downward_northward_stress', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_hfss = &
         register_diag_field ( mod_name, 'hfss', atmos_axes, Time, &
         'Surface Upward Sensible Heat Flux', 'W m-2', &
         standard_name = 'surface_upward_sensible_heat_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_hfls = &
         register_diag_field ( mod_name, 'hfls', atmos_axes, Time, &
         'Surface Upward Latent Heat Flux', 'W m-2', &
         standard_name = 'surface_upward_latent_heat_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_hfls > 0 ) call diag_field_add_attribute( id_hfls, 'comment', 'Lv*evap' )

    id_evspsbl = &
         register_diag_field( mod_name, 'evspsbl', atmos_axes, Time, &
         'Evaporation', 'kg m-2 s-1', &
         standard_name = 'water_evaporation_flux', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )

    id_sftlf = &
         register_static_field ( mod_name, 'sftlf', atmos_axes,  &
         'Fraction of the Grid Cell Occupied by Land', '1.0', &
         standard_name = 'land_area_fraction', area=area_id, &
         interp_method = "conserve_order1" )

    id_tslsi = &
         register_diag_field ( mod_name, 'tslsi', atmos_axes, Time,  &
         'Surface Temperature Where Land or Sea Ice', 'K', &
         standard_name = 'surface_temperature', area=area_id, &
         mask_variant=.true., missing_value=CMOR_MISSING_VALUE )

    id_tos = &
         register_diag_field ( mod_name, 'tos', atmos_axes, Time,  &
         'Sea Surface Temperature', 'K', &
         standard_name = 'sea_surface_temperature', area=area_id, &
         mask_variant=.true., missing_value=CMOR_MISSING_VALUE )

    id_sic = &
         register_diag_field ( mod_name, 'sic', atmos_axes, Time,  &
         'Sea Ice Area Fraction', '1.0', &
         standard_name = 'sea_ice_area_fraction', area=area_id, &
         missing_value=CMOR_MISSING_VALUE )
    if ( id_sic > 0 ) call diag_field_add_attribute( id_sic, 'comment', &
         'averaged over the ocean portion of grid box' )
#else
    id_tas = register_cmip_diag_field_2d ( mod_name, 'tas', Time, &
                            'Near-Surface Air Temperature', 'K' , &
                             standard_name='air_temperature' )
    if ( id_tas > 0 .and. id_height2m > 0) &
       call diag_field_add_attribute( id_tas, 'coordinates', 'height2m' )

    id_uas = register_cmip_diag_field_2d ( mod_name, 'uas', Time, &
                           'Eastward Near-Surface Wind', 'm s-1', &
                            standard_name='eastward_wind' )
    if ( id_uas > 0 .and. id_height10m > 0) &
       call diag_field_add_attribute( id_uas, 'coordinates', 'height10m' )

    id_vas = register_cmip_diag_field_2d ( mod_name, 'vas', Time, &
                          'Northward Near-Surface Wind', 'm s-1', &
                           standard_name='northward_wind' )
    if ( id_vas > 0 .and. id_height10m > 0 ) &
       call diag_field_add_attribute( id_vas, 'coordinates', 'height10m' )

    id_sfcWind = register_cmip_diag_field_2d ( mod_name, 'sfcWind', Time, &
                                      'Near-Surface Wind Speed', 'm s-1', &
                                       standard_name='wind_speed' )
    if ( id_sfcWind > 0 .and. id_height10m > 0 ) &
       call diag_field_add_attribute( id_sfcWind, 'coordinates', 'height10m' )

    id_huss = register_cmip_diag_field_2d ( mod_name, 'huss', Time, &
                           'Near-Surface Specific Humidity', '1.0', &
                            standard_name='specific_humidity' )
    if ( id_huss > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_huss, 'coordinates', 'height2m' )

    id_hurs = register_cmip_diag_field_2d ( mod_name, 'hurs', Time, &
                             'Near-Surface Relative Humidity', '%', &
                              standard_name='relative_humidity' )
    if ( id_hurs > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_hurs, 'coordinates', 'height2m' )

    id_rhs = register_cmip_diag_field_2d ( mod_name, 'rhs', Time, &
                           'Near-Surface Relative Humidity', '%', &
                            standard_name='relative_humidity' )
    if ( id_rhs > 0 .and. id_height2m > 0 ) &
       call diag_field_add_attribute( id_rhs, 'coordinates', 'height2m' )

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
    if ( id_hfls > 0 ) call diag_field_add_attribute( id_hfls, 'comment', 'Lv*evap' )

    id_evspsbl = register_cmip_diag_field_2d ( mod_name, 'evspsbl', Time, &
                                             'Evaporation', 'kg m-2 s-1', &
                                   standard_name='water_evaporation_flux' )

    id_sftlf = register_static_field ( mod_name, 'sftlf', atmos_axes,  &
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
    if ( id_sic > 0 ) call diag_field_add_attribute( id_sic, 'comment', &
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
         call diag_field_add_attribute ( get_global_diag_field_id(id_tas_g), 'coordinates', 'height2m' )

    id_tasl_g = register_global_diag_field ( 'tasl', Time, &
                           'Near-Surface Air Temperature (Land Only)', 'K' , &
                                 standard_name='air_temperature' )
#if defined(_USE_LEGACY_LAND_) || defined(use_AM3_physics)
    if(id_tasl_g>0) then
       call mpp_error(WARNING, "diag_field_init: field tasl is registered, but macro "// &
             "_USE_LEGACY_LAND_ or use_AM3_physics is defined, no data will be written out")
    endif
#endif
    if ( id_tasl_g > 0 .and. id_height2m > 0) &
         call diag_field_add_attribute ( get_global_diag_field_id(id_tasl_g), 'coordinates', 'height2m' )

    id_hfss_g = register_global_diag_field ( 'hfss', Time, &
                   'Surface Upward Sensible Heat Flux', 'W m-2', &
               standard_name='surface_upward_sensible_heat_flux' )

    id_hfls_g = register_global_diag_field ( 'hfls', Time, &
                  'Surface Upward Latent Heat Flux', 'W m-2', &
                  standard_name='surface_upward_latent_heat_flux')
    if ( id_hfls_g > 0 ) &
         call diag_field_add_attribute( get_global_diag_field_id(id_hfls_g), 'comment', 'Lv*evap' )

    id_rls_g = register_global_diag_field ( 'rls', Time, &
                   'Net Longwave Surface Radiation', 'W m-2', &
               standard_name='surface_net_longwave_flux' )

#endif
    !-----------------------------------------------------------------------

  end subroutine diag_field_init


  !######################################################################################
  !> \brief Divide data by area while avoiding zero area elements
  subroutine divide_by_area(data, area)
    real, intent(inout) :: data(:,:)
    real, intent(in)    :: area(:,:)

    if(size(data, dim=1) /= size(area, dim=1) .or. size(data, dim=2) /= size(area, dim=2)) then
       ! no op
       return
    endif

    where(area /= 0.0)
       data = data / area
    end where

  end subroutine divide_by_area

  !#######################################################################
  !> \brief Send out the ice_mask and/or sic data.
  !! This was called inside flux_ocean_to_ice. Why?
  subroutine send_ice_mask_sic(Time)
    type(time_type),         intent(in)  :: Time !< Current time

    real, dimension(nxc_ice, nyc_ice, nk_ice) :: ice_frac
    real, dimension(n_xgrid_sfc)              :: ex_ice_frac
    real, dimension(ni_atm, nj_atm)           :: diag_atm, ocean_frac
    logical :: used

    if ( id_ice_mask > 0 .or. id_sic > 0) then
       ice_frac        = 1.
       ice_frac(:,:,1) = 0.
       ex_ice_frac     = 0.
       call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
       call get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)
       if ( id_ice_mask > 0 ) used = send_data ( id_ice_mask, diag_atm, Time )

       ! ice concentration for only the ocean part of the atmos grid box
       ! normalize ice fraction over entire atmos grid box by the
       ! fraction of atmos grid box that is ocean
       if ( id_sic > 0) then
          ice_frac = 1.
          ex_ice_frac = 0.
          call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
          call get_from_xgrid (ocean_frac, 'ATM', ex_ice_frac, xmap_sfc)
          where (ocean_frac > 0.0)
             diag_atm = min(1., diag_atm/ocean_frac) ! CMIP6 as fraction
             ocean_frac = 1.0
          elsewhere
             diag_atm = 0.0
             ocean_frac = 0.0
          endwhere
          used = send_data ( id_sic, diag_atm, Time, rmask=ocean_frac )
       endif
    endif

  end subroutine send_ice_mask_sic

  !#######################################################################

  subroutine atm_stock_integrate(Atm, res)
    type(atmos_data_type), intent(in) :: Atm
    real,                 intent(out) :: res
    integer :: ier

    call stock_integrate_2d(Atm%lprec + Atm%fprec, xmap=xmap_sfc, delta_t=Dt_atm, &
         & radius=Radius, res=res, ier=ier)

  end subroutine atm_stock_integrate

!#########################################################################

end module atm_land_ice_flux_exchange_mod
