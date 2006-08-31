module flux_exchange_mod
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="V.Balaji@noaa.gov"> V. Balaji </CONTACT>
! <CONTACT EMAIL="Sergey.Malyshev@noaa.gov"> Sergey Malyshev </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, &
       mpp_error, stderr, stdout, stdlog, FATAL, NOTE, mpp_set_current_pelist, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
       CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_ROUTINE, lowercase
                    
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_global_sum, mpp_redistribute, operator(.EQ.)
  use mpp_io_mod,      only: mpp_close

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
!will support 3 types of flux_exchange:
!REGRID: physically distinct grids, via xgrid
!REDIST: same grid, transfer in index space only
!DIRECT: same grid, same decomp, direct copy
  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod, only: ocean_data_type, ice_ocean_boundary_type
  use ice_model_mod,   only: ice_data_type, land_ice_boundary_type, &
       ocean_ice_boundary_type, atmos_ice_boundary_type
  use    land_model_mod, only:  land_data_type, atmos_land_boundary_type

  use  surface_flux_mod, only: surface_flux
  use monin_obukhov_mod, only: mo_profile     

  use xgrid_mod, only: xmap_type, setup_xmap, set_frac_area, &
       put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init

  use diag_integral_mod, only:     diag_integral_field_init, &
       sum_diag_integral_field

  use  diag_manager_mod, only: register_diag_field,  &
       register_static_field, send_data, send_tile_averaged_data

  use  time_manager_mod, only: time_type

  use sat_vapor_pres_mod, only: escomp

  use      constants_mod, only: rdgas, rvgas, cp_air, stefan, WTMAIR

!Balaji
!utilities stuff into use fms_mod
  use fms_mod,                    only: clock_flag_default, check_nml_error, error_mesg
  use fms_mod,                    only: open_namelist_file, write_version_number

  use data_override_mod,          only: data_override
  use coupler_types_mod,          only: coupler_1d_bc_type
  use atmos_ocean_fluxes_mod,     only: atmos_ocean_fluxes_init, atmos_ocean_fluxes_calc
  use ocean_model_mod,            only: ocean_model_init_sfc, ocean_model_flux_init
  use coupler_types_mod,          only: coupler_type_copy
  use coupler_types_mod,          only: ind_psurf, ind_u10
  use atmos_tracer_driver_mod,    only: atmos_tracer_flux_init

  use field_manager_mod,          only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE
  use tracer_manager_mod,         only: get_tracer_index
  use tracer_manager_mod,         only: get_tracer_names, get_number_tracers, NO_TRACER

  implicit none
  include 'netcdf.inc'
private

  character(len=48), parameter :: module_name = 'flux_exchange_mod'

  public :: flux_exchange_init,   &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_ice_to_ocean,    &
     flux_ocean_to_ice

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: flux_exchange.F90,v 13.0.2.1.2.4 2006/07/13 12:29:57 fil Exp $'
  character(len=128) :: tag = '$Name: memphis_2006_08 $'
!-----------------------------------------------------------------------
!---- exchange grid maps -----

type(xmap_type), save :: xmap_sfc, xmap_runoff

integer         :: n_xgrid_sfc,  n_xgrid_runoff

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
     id_rough_moist, id_rough_heat, id_rough_mom,    &
     id_land_mask,   id_ice_mask,     &
     id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
     id_t_surf, id_t_flux, id_r_flux, id_q_flux,              &
     id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
     id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
     id_del_h,  id_del_m,  id_del_q,  id_rough_scale,         &
     id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm, id_gust, &
     id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land, &
     id_q_ref,  id_q_ref_land

integer, allocatable :: id_tr_atm(:), id_tr_surf(:), id_tr_flux(:), id_tr_mol_flux(:)

logical :: first_static = .true.
logical :: do_init = .true.
integer :: remap_method = 1

real, parameter :: bound_tol = 1e-7

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622

!--- namelist interface ------------------------------------------------------
! <NAMELIST NAME="flux_exchange_nml">
!   <DATA NAME="z_ref_heat"  TYPE="real"  DEFAULT="2.0">
!    eference height (meters) for temperature and relative humidity 
!    diagnostics (t_ref,rh_ref,del_h,del_q)
!   </DATA>
!   <DATA NAME="z_ref_mom"  TYPE="real"  DEFAULT="10.0">
!    reference height (meters) for momentum diagnostics (u_ref,v_ref,del_m)
!   </DATA>
!   <DATA NAME="ex_u_star_smooth_bug"  TYPE="logical"  DEFAULT="false">
!    By default, the global exchange grid u_star will not be interpolated from 
!    atmospheric grid, this is different from Jakarta behavior and will
!    change answers. So to perserve Jakarta behavior and reproduce answers
!    explicitly set this namelist variable to .true. in input.nml.
!    Talk to mw, ens for details.
!   </DATA>

  real ::  z_ref_heat =  2.,  &
           z_ref_mom  = 10.
  logical :: ex_u_star_smooth_bug = .false.
  logical :: sw1way_bug = .false.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom, ex_u_star_smooth_bug, sw1way_bug
! </NAMELIST>

! ---- allocatable module storage --------------------------------------------
real, allocatable, dimension(:) :: &
     ! NOTE: T canopy is only differet from t_surf over vegetated land
     ex_t_surf,    &   ! surface temperature for radiation calc, degK
     ex_t_ca,      &   ! near-surface (canopy) air temperature, degK
     ex_p_surf,    &   ! surface pressure

     ex_flux_t,    &   ! sens heat flux
     ex_flux_lw,   &   ! longwave radiation flux

     ex_dhdt_surf, &   ! d(sens.heat.flux)/d(T canopy)
     ex_dedt_surf, &   ! d(water.vap.flux)/d(T canopy)
     ex_dqsatdt_surf, &   ! d(water.vap.flux)/d(q canopy)
     ex_e_q_n,     &
     ex_drdt_surf, &   ! d(LW flux)/d(T surf)
     ex_dhdt_atm,  &   ! d(sens.heat.flux)/d(T atm)
     ex_flux_u,    &   ! u stress on atmosphere
     ex_flux_v,    &   ! v stress on atmosphere
     ex_dtaudu_atm,&   ! d(stress)/d(u)
     ex_dtaudv_atm,&   ! d(stress)/d(v)
     ex_albedo_fix,&
     ex_albedo_vis_dir_fix,&
     ex_albedo_nir_dir_fix,&
     ex_albedo_vis_dif_fix,&
     ex_albedo_nir_dif_fix,&
     ex_old_albedo,&   ! old value of albedo for downward flux calculations
     ex_drag_q,    &   ! q drag.coeff.
     ex_cd_t,      &
     ex_cd_m,      &
     ex_b_star,    &
     ex_u_star,    &
     ex_wind,      &
     ex_z_atm

real, allocatable, dimension(:,:) :: &
     ex_tr_surf,    & ! near-surface tracer fields
     ex_flux_tr,    & ! tracer fluxes
     ex_dfdtr_surf, & ! d(tracer flux)/d(surf tracer)
     ex_dfdtr_atm,  & ! d(tracer flux)/d(atm tracer)
     ex_e_tr_n,     & ! coefficient in implicit scheme 
     ex_f_tr_delt_n   ! coefficient in implicit scheme

logical, allocatable, dimension(:) :: &
     ex_avail,     &   ! true where data on exchange grid are available
     ex_land           ! true if exchange grid cell is over land
real, allocatable, dimension(:) :: &
     ex_e_t_n,      &
     ex_f_t_delt_n

integer :: n_atm_tr  ! number of prognostic tracers in the atmos model
integer :: n_atm_tr_tot  ! number of prognostic tracers in the atmos model
integer :: n_lnd_tr  ! number of prognostic tracers in the land model 
integer :: n_lnd_tr_tot  ! number of prognostic tracers in the land model 
integer :: n_exch_tr ! number of tracers exchanged between models

type :: tracer_ind_type
   integer :: atm, ice, lnd ! indices of the tracer in the respective models
end type 
type(tracer_ind_type), allocatable :: tr_table(:) ! table of tracer indices
type :: tracer_exch_ind_type
   integer :: exch = 0  ! exchange grid index
   integer :: ice = 0   ! ice model index
   integer :: lnd = 0   ! land model index
end type tracer_exch_ind_type
type(tracer_exch_ind_type), allocatable :: tr_table_map(:) ! map atm tracers to exchange, ice and land variables
integer :: isphum = NO_TRACER       ! index of specific humidity tracer in tracer table

type(coupler_1d_bc_type), save        :: ex_gas_fields_atm  ! gas fields in atm
                     ! Place holder for various atmospheric fields.
type(coupler_1d_bc_type), save        :: ex_gas_fields_ice  ! gas fields on ice
type(coupler_1d_bc_type), save        :: ex_gas_fluxes      ! gas flux
                     ! Place holder of intermediate calculations, such as
                     ! piston velocities etc.

integer :: ni_atm, nj_atm ! to do atmos diagnostic from flux_ocean_to_ice
real, dimension(3) :: ccc ! for conservation checks
!Balaji, sets boundary_type%xtype
!  REGRID: grids are physically different, pass via exchange grid
!  REDIST: same physical grid, different decomposition, must move data around
!  DIRECT: same physical grid, same domain decomposition, can directly copy data
integer, parameter :: REGRID=1, REDIST=2, DIRECT=3
!Balaji: clocks moved into flux_exchange
  integer :: cplClock, sfcClock, fluxAtmDnClock, fluxLandIceClock, &
             fluxIceOceanClock, fluxOceanIceClock, regenClock, fluxAtmUpClock, &
             cplOcnClock

  logical :: ocn_pe, ice_pe
  integer, allocatable, dimension(:) :: ocn_pelist, ice_pelist

contains
!  coupler_main control loop
!  --------------------
!
!       DO slow time steps (ocean)
!
!           call flux_ocean_to_ice            (flux_exchange)
!
!           call update_ice_model_slow_up     (ice_model)
!
!           DO fast time steps (atmos)
!
!                call sfc_boundary_layer      (flux_exchange)
!
!                call update_atmos_model_down (atmos_model)
!
!                call flux_down_from_atmos    (flux_exchange)
!
!                call update_land_model_fast  (land_model)
!
!                call flux_land_to_ice        (flux_exchange)
!
!                call update_ice_model_fast   (ice_model)
!
!                call flux_up_to_atmos        (flux_exchange)
!
!                call update_atmos_model_up   (atmos_model)
!
!           END DO
!
!           call update_ice_model_slow_dn     (ice_model)
!
!           call flux_ice_to_ocean            (flux_exchange)
!
!           call update_ocean_model           (ocean_model)
!
!      END DO
!
!   LAND_FAST and ICE_FAST must update the surface temperature
!
! There are 7 subroutines defined in the following 5 include files
! that form the core operations of the coupler
!
! flux_exchange_init
! sfc_boundary_layer
! flux_down_from_atmos
! flux_up_to_atmos
! flux_land_to_ice
! flux_ice_to_ocean
! flux_ocean_to_ice
!
include 'flux_exchange_init.inc'
include 'flux_exchange_sbl.inc'
include 'flux_exchange_atmos.inc'
include 'flux_exchange_land_ice.inc'
include 'flux_exchange_ice_ocean.inc'

!#######################################################################
! <SUBROUTINE NAME="generate_sfc_xgrid">
!  <OVERVIEW>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </OVERVIEW>
!  <DESCRIPTION>
!   Optimizes the exchange grids by eliminating land and ice partitions with no data. 
!  </DESCRIPTION>
!  <TEMPLATE>
!   call generate_sfc_xgrid( Land, Ice )
!		
!  </TEMPLATE>
!  <IN NAME=" Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!  A derived data type to specify ice boundary data.
!  </IN>
!
subroutine generate_sfc_xgrid( Land, Ice )
! subroutine to regenerate exchange grid eliminating side 2 tiles with 0 frac area
    type(land_data_type), intent(in) :: Land
    type(ice_data_type),  intent(in) :: Ice

    integer :: isc, iec, jsc, jec

!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(regenClock)

  call mpp_get_compute_domain(Ice%Domain, isc, iec, jsc, jec)

  call set_frac_area (Ice%part_size(isc:iec,jsc:jec,:) , 'OCN', xmap_sfc)
  call set_frac_area (Land%tile_size, 'LND', xmap_sfc)
  n_xgrid_sfc = max(xgrid_count(xmap_sfc),1)

!Balaji
  call mpp_clock_end(regenClock)
  call mpp_clock_end(cplClock)
  return
end subroutine generate_sfc_xgrid
! </SUBROUTINE>

!#######################################################################

subroutine put_logical_to_real (mask, id, ex_mask, xmap)

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

end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes, land_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)
  integer,         intent(in) :: land_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
       vrange = (/ -400., 400. /), &
       frange = (/ -0.01, 1.01 /)
  character(len=32)  :: name, units ! name of the tracer
  character(len=128) :: longname    ! long name of the tracer
  integer            :: tr          ! tracer index
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
       range=frange )
  
  !--------- initialize diagnostic fields --------------------

  id_ice_mask = &
       register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
       'fractional amount of sea ice', 'none',  &
       range=frange )
  
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

  id_u_ref      = &
       register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
       'zonal wind component at '//label_zm,  'm/s', &
       range=vrange )

  id_v_ref      = &
       register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
       'meridional wind component at '//label_zm, 'm/s', &
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

  ! + slm Jun 02, 2002 -- diagnostics of reference values over the land
  id_t_ref_land = &
       register_diag_field ( mod_name, 't_ref_land', Land_axes, Time, &
       'temperature at '//trim(label_zh)//' over land', 'deg_k' , &
       range=trange, missing_value =  -100.0)
  id_rh_ref_land= &
       register_diag_field ( mod_name, 'rh_ref_land', Land_axes, Time,   &
       'relative humidity at '//trim(label_zh)//' over land', 'percent',       &
       missing_value=-999.0)
  id_u_ref_land = &
       register_diag_field ( mod_name, 'u_ref_land',  Land_axes, Time, &
       'zonal wind component at '//trim(label_zm)//' over land',  'm/s', &
       range=vrange, missing_value=-999.0 )
  id_v_ref_land = &
       register_diag_field ( mod_name, 'v_ref_land',  Land_axes, Time,     &
       'meridional wind component at '//trim(label_zm)//' over land', 'm/s', &
       range=vrange, missing_value = -999.0 )
  ! - slm Jun 02, 2002
  id_q_ref = &
       register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
       'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)
  id_q_ref_land = &
       register_diag_field ( mod_name, 'q_ref_land', Land_axes, Time, &
       'specific humidity at '//trim(label_zh)//' over land', 'kg/kg',          &
       missing_value=-1.0)

  id_rough_scale = &
       register_diag_field ( mod_name, 'rough_scale', atmos_axes, Time, &
       'topographic scaling factor for momentum drag','1' )
!-----------------------------------------------------------------------

  allocate(id_tr_atm(n_exch_tr))
  allocate(id_tr_surf(n_exch_tr))
  allocate(id_tr_flux(n_exch_tr))
  allocate(id_tr_mol_flux(n_exch_tr))

  do tr = 1, n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, name, longname, units )
     id_tr_atm(tr) = register_diag_field (mod_name, trim(name)//'_atm', atmos_axes, Time, &
          trim(longname)//' at btm level', trim(units))
     id_tr_surf(tr) = register_diag_field (mod_name, trim(name)//'_surf', atmos_axes, Time, &
          trim(longname)//' at the surface', trim(units))
     id_tr_flux(tr) = register_diag_field(mod_name, trim(name)//'_flux', atmos_axes, Time, &
          'flux of '//trim(longname), trim(units)//' kg/(m2 s)')
     id_tr_mol_flux(tr) = register_diag_field(mod_name, trim(name)//'_mol_flux', atmos_axes, Time, &
          'flux of '//trim(longname), 'mol/(m2 s)')
  enddo
  id_q_flux = register_diag_field( mod_name, 'evap',       atmos_axes, Time, &
         'evaporation rate',        'kg/m2/s'  )

  end subroutine diag_field_init
! <OVERVIEW>
!   The flux_exchange module provides interfaces to couple the following component 
!   models: atmosphere, ocean, land, and ice. All interpolation between physically 
!   distinct model grids is handled by the exchange grid (xgrid_mod) with the 
!   interpolated quantities being conserved.
! </OVERVIEW>

! <DESCRIPTION>
!  <PRE>
!  1.This version of flux_exchange_mod allows the definition of physically independent
!    grids for atmosphere, land and sea ice. Ice and ocean must share the same physical
!    grid (though the domain decomposition on parallel systems may be different). 
!    Grid information is input through the grid_spec file (URL). The masked region of the
!    land grid and ice/ocean grid must "tile" each other. The masked region of the ice grid
!    and ocean grid must be identical. 
!
!         ATMOSPHERE  |----|----|----|----|----|----|----|----|
!
!               LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
!
!                ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!               OCEAN |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!
!              where  |xxx| = masked grid point
!         
!
!    The atmosphere, land, and ice grids exchange information using the exchange grid xmap_sfc.
!
!    The land and ice grids exchange runoff data using the exchange grid xmap_runoff.
!
!    Transfer of data between the ice bottom and ocean does not require an exchange 
!    grid as the grids are physically identical. The flux routines will automatically
!    detect and redistribute data if their domain decompositions are different.
!
!    To get information from the atmosphere to the ocean it must pass through the 
!    ice model, first by interpolating from the atmospheric grid to the ice grid, 
!    and then transferring from the ice grid to the ocean grid.

!  2.Each component model must have a public defined data type containing specific 
!    boundary fields. A list of these quantities is located in the NOTES of this document. 
!
!  3.The surface flux of sensible heat and surface evaporation can be implicit functions
!    of surface temperature. As a consequence, the parts of the land and sea-ice models 
!    that update the surface temperature must be called on the atmospheric time step 
!
!  4.The surface fluxes of all other tracers and of momentum are assumed to be explicit
!    functions of all surface parameters 
!
!  5.While no explicit reference in made within this module to the implicit treatment 
!    of vertical diffusion in the atmosphere and in the land or sea-ice models, the 
!    module is designed to allow for simultaneous implicit time integration on both 
!    sides of the surface interface. 
!
!  6.Due to #5, the diffusion part of the land and ice models must be called on the 
!    atmospheric time step.
  
!7. Any field passed from one component to another may be "faked" to a
!   constant value, or to data acquired from a file, using the
!   data_override feature of FMS. The fields to override are runtime
!   configurable, using the text file <tt>data_table</tt> for input.
!   See the data_override_mod documentation for more details.
!
!   We DO NOT RECOMMEND exercising the data override capabilities of
!   the FMS coupler until the user has acquired considerable
!   sophistication in running FMS.
!
!   Here is a listing of the override capabilities of the flux_exchange
!   module:
!
!   FROM the atmosphere boundary TO the exchange grid (in sfc_boundary_layer):
!  
!        t_bot, q_bot, z_bot, p_bot, u_bot, v_bot, p_surf, gust
!
!   FROM the ice boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, rough_mom, rough_heat, rough_moist, albedo, u_surf, v_surf
!     
!   FROM the land boundary TO the exchange grid (in sfc_boundary_layer):
!
!        t_surf, t_ca, q_ca, rough_mom, rough_heat, albedo
!
!   FROM the exchange grid TO land_ice_atmos_boundary (in
!   sfc_boundary_layer):
!
!        t, albedo, land_frac, dt_t, dt_q, u_flux, v_flux, dtaudu, dtaudv,
!        u_star, b_star, rough_mom
!   
!   FROM the atmosphere boundary TO the exchange grid (in
!    flux_down_from_atmos):
!
!        flux_sw, flux_lw, lprec, fprec, coszen, dtmass, delta_t,
!        delta_q, dflux_t, dflux_q
!        
!   FROM the exchange grid TO the land boundary (in
!    flux_down_from_atmos):
!
!    t_flux, q_flux, lw_flux, sw_flux, lprec, fprec, dhdt, dedt, dedq,
!    drdt, drag_q, p_surf
!    
!   FROM the exchange grid TO the ice boundary (in flux_down_from_atmos):
!
!        u_flux, v_flux, t_flux, q_flux, lw_flux, lw_flux_dn, sw_flux,
!        sw_flux_dn, lprec, fprec, dhdt, dedt, drdt, coszen, p 
!
!   FROM the land boundary TO the ice boundary (in flux_land_to_ice):
!
!        runoff, calving
!
!   FROM the ice boundary TO the ocean boundary (in flux_ice_to_ocean):
! 
!        u_flux, v_flux, t_flux, q_flux, salt_flux, lw_flux, sw_flux,
!        lprec, fprec, runoff, calving, p
!        
!   FROM the ocean boundary TO the ice boundary (in flux_ocean_to_ice):
!
!        u, v, t, s, frazil, sea_level
!
!   FROM the ice boundary TO the atmosphere boundary (in flux_up_to_atmos):
!
!        t_surf
!
!   FROM the land boundary TO the atmosphere boundary (in
!    flux_up_to_atmos):
!  
!        t_ca, t_surf, q_ca
!
!  See NOTES below for an explanation of the field names.
!  </PRE>
! </DESCRIPTION>
! <DIAGFIELDS>
!   <NETCDF NAME="land_mask" UNITS="none">
!     fractional amount of land
!   </NETCDF>
!   <NETCDF NAME="wind" UNITS="m/s">
!     wind speed for flux calculations
!   </NETCDF>
!   <NETCDF NAME="drag_moist" UNITS="none">
!     drag coeff for moisture
!   </NETCDF>
!   <NETCDF NAME="drag_heat" UNITS="none">
!     drag coeff for heat
!   </NETCDF>
!   <NETCDF NAME="drag_mom" UNITS="none">
!     drag coeff for momentum
!   </NETCDF>
!   <NETCDF NAME="rough_moist" UNITS="m">
!     surface roughness for moisture
!   </NETCDF>
!   <NETCDF NAME="rough_heat" UNITS="m">
!     surface roughness for heat
!   </NETCDF>
!   <NETCDF NAME="rough_mom" UNITS="m">
!     surface roughness for momentum
!   </NETCDF>
!   <NETCDF NAME="u_star" UNITS="m/s">
!     friction velocity
!   </NETCDF>
!   <NETCDF NAME="b_star" UNITS="m/s">
!     buoyancy scale
!   </NETCDF>
!   <NETCDF NAME="q_star" UNITS="kg water/kg air">
!     moisture scale
!   </NETCDF>
!   <NETCDF NAME="t_atm" UNITS="deg_k">
!     temperature at btm level
!   </NETCDF>
!   <NETCDF NAME="u_atm" UNITS="m/s">
!     u wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="v_atm" UNITS="m/s">
!     v wind component at btm level
!   </NETCDF>
!   <NETCDF NAME="q_atm" UNITS="kg/kg">
!     specific humidity at btm level
!   </NETCDF>
!   <NETCDF NAME="p_atm" UNITS="pa">
!     pressure at btm level
!   </NETCDF>
!   <NETCDF NAME="z_atm" UNITS="m">
!     height of btm level
!   </NETCDF>
!   <NETCDF NAME="gust" UNITS="m/s">
!     gust scale 
!   </NETCDF>
!   <NETCDF NAME="rh_ref" UNITS="percent">
!     relative humidity at ref height
!   </NETCDF>
!   <NETCDF NAME="t_ref" UNITS="deg_k">
!    temperature at ref height
!   </NETCDF>
!   <NETCDF NAME="u_ref" UNITS="m/s">
!    zonal wind component at ref height
!   </NETCDF>
!   <NETCDF NAME="v_ref" UNITS="m/s">
!    meridional wind component at ref height 
!   </NETCDF>
!   <NETCDF NAME="del_h" UNITS="none">
!    ref height interp factor for heat 
!   </NETCDF>
!   <NETCDF NAME="del_m" UNITS="none">
!    ref height interp factor for momentum 
!   </NETCDF>
!   <NETCDF NAME="del_q" UNITS="none">
!    ref height interp factor for moisture
!   </NETCDF>
!   <NETCDF NAME="tau_x" UNITS="pa">
!    zonal wind stress
!   </NETCDF>
!   <NETCDF NAME="tau_y" UNITS="pa">
!    meridional wind stress
!   </NETCDF>
!   <NETCDF NAME="ice_mask" UNITS="none">
!    fractional amount of sea ice 
!   </NETCDF>
!   <NETCDF NAME="t_surf" UNITS="deg_k">
!     surface temperature
!   </NETCDF>
!   <NETCDF NAME="t_ca" UNITS="deg_k">
!     canopy air temperature
!   </NETCDF>
!   <NETCDF NAME="q_surf" UNITS="kg/kg">
!     surface specific humidity 
!   </NETCDF>
!   <NETCDF NAME="shflx" UNITS="w/m2">
!     sensible heat flux
!   </NETCDF>
!   <NETCDF NAME="evap" UNITS="kg/m2/s">
!     evaporation rate 
!   </NETCDF>
!   <NETCDF NAME="lwflx" UNITS="w/m2">
!    net (down-up) longwave flux 
!   </NETCDF>

! </DIAGFIELDS>

! <INFO>


!   <NOTE>
!   <PRE>
!
!  MAIN PROGRAM EXAMPLE
!  --------------------
!
!       DO slow time steps (ocean)
!
!           call flux_ocean_to_ice
!
!           call ICE_SLOW_UP
!
!
!           DO fast time steps (atmos)
!
!                call sfc_boundary_layer
!
!                call ATMOS_DOWN
!
!                call flux_down_from_atmos
!
!                call LAND_FAST
!
!                call ICE_FAST
!
!                call flux_up_to_atmos
!
!                call ATMOS_UP
!
!           END DO
!
!           call ICE_SLOW_DN
!
!           call flux_ice_to_ocean
!
!           call OCEAN
!
!      END DO
!
!   LAND_FAST and ICE_FAST must update the surface temperature
!
! =======================================================================
!
! REQUIRED VARIABLES IN DEFINED DATA TYPES FOR COMPONENT MODELS
! --------------------------------------------------------------
!
! type (atmos_boundary_data_type) :: Atm
! type (surf_diff_type) :: Atm%Surf_Diff
!
! real, dimension(:)
!
!    Atm%lon_bnd   longitude axis grid box boundaries in radians
!                  must be monotonic
!    Atm%lat_bnd   latitude axis grid box boundaries in radians
!                  must be monotonic
!
! real, dimension(:,:)
!
!    Atm%t_bot     temperature at lowest model level
!    Atm%q_bot     specific humidity at lowest model level
!    Atm%z_bot     height above the surface for the lowest model level (m)
!    Atm%p_bot     pressure at lowest model level (pa)
!    Atm%u_bot     zonal wind component at lowest model level (m/s)
!    Atm%v_bot     meridional wind component at lowest model level (m/s)
!    Atm%p_surf    surface pressure (pa)
!    Atm%gust      gustiness factor (m/s)
!    Atm%flux_sw   net shortwave flux at the surface
!    Atm%flux_lw   downward longwave flux at the surface
!    Atm%lprec     liquid precipitation (kg/m2)
!    Atm%fprec     water equivalent frozen precipitation (kg/m2)
!    Atm%coszen    cosine of the zenith angle
!
!   (the following five fields are gathered into a data type for convenience in passing
!   this information through the different levels of the atmospheric model --
!   these fields are rlated to the simultaneous implicit time steps in the
!   atmosphere and surface models -- they are described more fully in
!   flux_exchange.tech.ps and
!   in the documntation for vert_diff_mod
!
!
!    Atm%Surf_Diff%dtmass   = dt/mass where dt = atmospheric time step ((i+1) = (i-1) for leapfrog) (s)
!                           mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!    Atm%Surf_Diff%delta_t  increment ((i+1) = (i-1) for leapfrog) in temperature of
!                           lowest atmospheric layer  (K)
!    Atm%Surf_Diff%delta_q  increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!                           lowest atmospheric layer (nondimensional -- Kg/Kg)
!    Atm%Surf_Diff%dflux_t  derivative of implicit part of downward temperature flux at top of lowest
!                           atmospheric layer with respect to temperature
!                           of lowest atmospheric layer (Kg/(m2 s))
!    Atm%Surf_Diff%dflux_q  derivative of implicit part of downward moisture flux at top of lowest
!                           atmospheric layer with respect to specific humidity of
!                           of lowest atmospheric layer (Kg/(m2 s))
!
!
! integer, dimension(4)
!
!    Atm%axes      Axis identifiers returned by diag_axis_init for the
!                  atmospheric model axes: X, Y, Z_full, Z_half.
!
! -----------------------------------------------
!
! type (land_boundary_data_type) :: Land
!
! real, dimension(:)
!
!    Land%lon_bnd     longitude axis grid box boundaries in radians
!                     must be monotonic
!    Land%lat_bnd     latitude axis grid box boundaries in radians
!                     must be monotonic
!
! logical, dimension(:,:,:)
!
!    Land%mask        land/sea mask (true for land)
!    Land%glacier     glacier mask  (true for glacier)
!
! real, dimension(:,:,:)
!
!    Land%tile_size   fractional area of each tile (partition)
!
!    Land%t_surf      surface temperature (deg k)
!    Land%albedo      surface albedo (fraction)
!    Land%rough_mom   surface roughness for momentum (m)
!    Land%rough_heat  surface roughness for heat/moisture (m)
!    Land%stomatal    stomatal resistance
!    Land%snow        snow depth (water equivalent) (kg/m2)
!    Land%water       water depth of the uppermost bucket (kg/m2)
!    Land%max_water   maximum water depth allowed in the uppermost bucket (kg/m2)
!
! -----------------------------------------------
!
!
! type (ice_boundary_data_type) :: Ice
!
! real, dimension(:)
!
!    Ice%lon_bnd       longitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lat_bnd       latitude axis grid box boundaries for temperature points
!                      in radians (must be monotonic)
!    Ice%lon_bnd_uv    longitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!    Ice%lat_bnd_uv    latitude axis grid box boundaries for momentum points
!                      in radians (must be monotonic)
!
! logical, dimension(:,:,:)
!
!    Ice%mask          ocean/land mask for temperature points
!                        (true for ocean, with or without ice)
!    Ice%mask_uv       ocean/land mask for momentum points
!                        (true for ocean, with or without ice)
!    Ice%ice_mask      optional ice mask (true for ice)
!
! real, dimension(:,:,:)
!
!    Ice%part_size     fractional area of each partition of a temperature grid box
!    Ice%part_size_uv  fractional area of each partition of a momentum grid box
!
!    the following fields are located on the ice top grid
!
!    Ice%t_surf        surface temperature (deg k)
!    Ice%albedo        surface albedo (fraction)
!    Ice%rough_mom     surface roughness for momentum (m)
!    Ice%rough_heat    surface roughness for heat/moisture (m)
!    Ice%u_surf        zonal (ocean/ice) current at the surface (m/s)
!    Ice%v_surf        meridional (ocean/ice) current at the surface (m/s)
!
!    the following fields are located on the ice bottom grid
!
!    Ice%flux_u        zonal wind stress (Pa)
!    Ice%flux_v        meridional wind stress (Pa)
!    Ice%flux_t        sensible heat flux (w/m2)
!    Ice%flux_q        specific humidity flux (kg/m2/s)
!    Ice%flux_sw       net (down-up) shortwave flux (w/m2)
!    Ice%flux_lw       net (down-up) longwave flux (w/m2)
!    Ice%lprec         mass of liquid precipitation since last time step (Kg/m2)
!    Ice%fprec         mass of frozen precipitation since last time step (Kg/m2)
!    Ice%runoff        mass of runoff water since last time step (Kg/m2)
!
! -----------------------------------------------
!
! type (ocean_boundary_data_type) :: Ocean
!
! real, dimension(:)
!
!    Ocean%Data%lon_bnd      longitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd      latitude axis grid box boundaries for temperature
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lon_bnd_uv   longitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!    Ocean%Data%lat_bnd_uv   latitude axis grid box boundaries for momentum
!                            points on the ocean DATA GRID (radians)
!
!    Ocean%Ocean%lon_bnd     longitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd     latitude axis grid box boundaries for temperature
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lon_bnd_uv  longitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!    Ocean%Ocean%lat_bnd_uv  latitude axis grid box boundaries for momentum
!                            points on the ocean MODEL GRID (radians)
!
!      Note: The data values in all longitude and latitude grid box boundary
!            array must be monotonic.
!
! logical, dimension(:,:)
!
!    Ocean%Data%mask       ocean/land mask for temperature points on the ocean
!                          DATA GRID (true for ocean)
!    Ocean%Data%mask_uv    ocean/land mask for momentum points on the ocean
!                          DATA GRID (true for ocean)
!
!    Ocean%Ocean%mask      ocean/land mask for temperature points on the ocean
!                          MODEL GRID (true for ocean)
!    Ocean%Ocean%mask_uv   ocean/land mask for momentum points on the ocean
!                          MODEL GRID (true for ocean)
!
! real, dimension(:,:)
!
!    Ocean%t_surf_data  surface temperature on the ocean DATA GRID (deg k)
!
!    Ocean%t_surf       surface temperature on the ocean MODEL GRID (deg k)
!    Ocean%u_surf       zonal ocean current at the surface on the ocean
!                       MODEL GRID (m/s)
!    Ocean%v_surf       meridional ocean current at the surface on the
!                       ocean MODEL GRID (m/s)
!    Ocean%frazil       frazil at temperature points on the ocean MODEL GRID
!
!   </PRE>
!   </NOTE>
! </INFO>

end module flux_exchange_mod

