

module flux_exchange_mod

use atmos_coupled_mod, only: atmos_boundary_data_type
use ocean_coupled_mod, only: ocean_boundary_data_type
use   ice_coupled_mod, only:   ice_boundary_data_type
use    land_model_mod, only:  land_boundary_data_type

use surface_flux_mod, only: surface_flux, surface_profile     

use     exchange_mod, only: boundary_map_type,           &
                            exchange_map_type,           &
                            init_boundary_map,           &
                            lon_lat_size, lon_lat_map,   &
                            complete_side1_boundary_map, &
                            complete_side2_boundary_map, &
                            get_exchange_grid_size,      &
                            get_exchange_grid,           &
                            put_exchange_grid,           &
                            set_frac_area

use diag_integral_mod, only:     diag_integral_field_init, &
                             sum_diag_integral_field

use     utilities_mod, only: file_exist, open_file, check_nml_error,  &
                             error_mesg, FATAL, get_my_pe, close_file

use  diag_manager_mod, only: register_diag_field,  &
                             register_static_field, send_data

use  time_manager_mod, only: time_type

use sat_vapor_pres_mod, only: escomp

use      constants_mod, only: rdgas, rvgas, cp

implicit none
private

public :: flux_exchange_init,   &
          flux_calculation,     &
          flux_down_from_atmos, &
          flux_up_to_atmos,     &
          flux_ice_to_ocean,    &
          flux_ocean_to_ice

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id: flux_exchange.F90,v 1.5 2001/07/05 17:42:57 fms Exp $'
character(len=128) :: tag = '$Name: eugene $'
!-----------------------------------------------------------------------
!---- boundary maps and exchange grid maps -----

type (boundary_map_type) :: bd_map_atm, bd_map_land,                 &
                            bd_map_ocean_data, bd_map_ocean_model,   &
                            bd_map_ice_top, bd_map_ice_bot,          &
                            bd_map_ice_bot_uv, bd_map_ocean_data_uv, &
                            bd_map_ocean_model_uv

type(exchange_map_type), target :: ex_map_top, ex_map_bot, ex_map_bot_uv

integer :: ex_num_top, ex_num_bot, ex_num_bot_uv

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
           id_rough_moist, id_rough_heat, id_rough_mom,    &
           id_land_mask,   id_glac_mask,  id_ice_mask,     &
           id_u_star, id_b_star, id_u_flux, id_v_flux, id_t_surf,   &
           id_t_flux, id_q_flux, id_r_flux,                         &
           id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
           id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
           id_del_h,  id_del_m,  id_del_q

logical :: first_static = .true.
logical :: do_init = .true.

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622

!-----------------------------------------------------------------------

real ::  z_ref_heat =  2.,  &
         z_ref_mom  = 10.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom

!-----------------------------------------------------------------------

!---- allocatable module storage ------

  real, allocatable, dimension(:) ::  ex_t_surf,      &
          ex_dhdt_surf, ex_dedt_surf, ex_drdt_surf,   &
          ex_dhdt_atm,  ex_dedq_atm,                  &
          ex_flux_t, ex_flux_q, ex_flux_lw, ex_albedo_fix
  logical, allocatable, dimension(:) ::  ex_avail

  real, allocatable, dimension(:) :: ex_e_t_n,  ex_f_t_delt_n, &
                                     ex_e_q_n,  ex_f_q_delt_n           

contains

!#######################################################################

 subroutine flux_calculation ( dt, Time, Atm, Land, Ice, land_frac_atm,&
                               t_surf_atm, albedo_atm, rough_mom_atm,  &
                               flux_u_atm, flux_v_atm, dtaudv_atm,     &
                               u_star_atm, b_star_atm                  )

 real,                   intent(in)  :: dt
 type       (time_type), intent(in)  :: Time
 type (atmos_boundary_data_type), intent(in)  :: Atm
 type  (land_boundary_data_type), intent(in)  :: Land
 type   (ice_boundary_data_type), intent(in)  :: Ice
 real, dimension(:,:),   intent(out) :: t_surf_atm, albedo_atm,    &
                                        rough_mom_atm,             &
                                        land_frac_atm, dtaudv_atm, &
                                        flux_u_atm, flux_v_atm,    &
                                        u_star_atm, b_star_atm

 real, dimension(ex_num_top) :: ex_albedo, ex_land_frac,  &
       ex_t_atm,  ex_q_atm, ex_z_atm, ex_p_atm, ex_u_atm, ex_v_atm,  &
       ex_p_surf, ex_gust, ex_t_surf4, ex_u_surf, ex_v_surf,         &
       ex_rough_mom, ex_rough_heat, ex_rough_moist,                  &
       ex_stomatal, ex_snow, ex_water, ex_max_water,                 &
       ex_u_star, ex_b_star, ex_q_star, ex_q_surf, ex_glac_frac,     &
       ex_cd_q, ex_cd_t, ex_cd_m, ex_flux_u, ex_flux_v, ex_dtaudv_atm,&
       ex_ref, ex_t_ref, ex_qs_ref, ex_del_m, ex_del_h, ex_del_q,     &
       ex_wind

 logical, dimension(ex_num_top) :: ex_glacier, ex_land
 real, dimension(size(t_surf_atm,1),size(t_surf_atm,2)) :: diag_atm

 real    :: zrefm, zrefh
 logical :: used


!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('flux_exchange_mod',  &
                 'must call flux_exchange_init first', FATAL)

!-----------------------------------------------------------------------
!------ allocate storage also needed in flux_up_to_atmos -----

   allocate ( ex_t_surf   (ex_num_top),  &
              ex_dhdt_surf(ex_num_top),  &
              ex_dedt_surf(ex_num_top),  &
              ex_drdt_surf(ex_num_top),  &
              ex_dhdt_atm (ex_num_top),  &
              ex_dedq_atm (ex_num_top),  &
              ex_flux_t   (ex_num_top),  &
              ex_flux_q   (ex_num_top),  &
              ex_flux_lw  (ex_num_top),  &
              ex_avail    (ex_num_top)   )

   allocate ( ex_f_t_delt_n   (ex_num_top),  &
              ex_f_q_delt_n   (ex_num_top),  &
              ex_e_t_n        (ex_num_top),  &
              ex_e_q_n        (ex_num_top) )

!--- initialize some values ---
   ex_t_surf   = 200.
   ex_u_surf   =   0.
   ex_v_surf   =   0.
   ex_stomatal =   0.
   ex_snow     =   0.

!---- do not use if relax time /= 0 ----
   ex_cd_t = 0.0
   ex_cd_m = 0.0
   ex_cd_q = 0.0
!-----------------------------------------------------------------------
!---- put atmosphere quantities onto exchange grid ----

   call put_exchange_grid (Atm%t_bot , ex_t_atm , bd_map_atm)
   call put_exchange_grid (Atm%q_bot , ex_q_atm , bd_map_atm)
   call put_exchange_grid (Atm%z_bot , ex_z_atm , bd_map_atm)
   call put_exchange_grid (Atm%p_bot , ex_p_atm , bd_map_atm)
   call put_exchange_grid (Atm%u_bot , ex_u_atm , bd_map_atm)
   call put_exchange_grid (Atm%v_bot , ex_v_atm , bd_map_atm)
   call put_exchange_grid (Atm%p_surf, ex_p_surf, bd_map_atm)
   call put_exchange_grid (Atm%gust,   ex_gust,   bd_map_atm)

!-----------------------------------------------------------------------
!---- put land quantities onto exchange grid ----

   call set_frac_area (Land%tile_size, bd_map_land)

   ex_land = .false.
   call put_exchange_grid  &
             (Land%t_surf       , ex_t_surf     , bd_map_land,  &
                                                  avail=ex_land)
   call put_exchange_grid  &
             (Land%rough_mom    , ex_rough_mom  , bd_map_land)
   call put_exchange_grid  &
             (Land%rough_heat   , ex_rough_heat , bd_map_land)
   call put_exchange_grid  &
             (Land%rough_heat   , ex_rough_moist, bd_map_land)
   call put_exchange_grid  &
             (Land%albedo       , ex_albedo     , bd_map_land)
   call put_exchange_grid  &
             (Land%stomatal     , ex_stomatal   , bd_map_land)
   call put_exchange_grid  &
             (Land%snow         , ex_snow       , bd_map_land)
   call put_exchange_grid  &
             (Land%water        , ex_water      , bd_map_land)
   call put_exchange_grid  &
             (Land%max_water    , ex_max_water  , bd_map_land)

   ex_land_frac = 0.0
   call put_logical_to_real      &
             (Land%mask         , ex_land_frac  , bd_map_land)
   ex_glac_frac = 0.0
   call put_logical_to_real      &
             (Land%glacier      , ex_glac_frac  , bd_map_land)

   where (.not.ex_land) ex_glac_frac = 0.0
   ex_glacier = ex_glac_frac > 0.5

!-----------------------------------------------------------------------
!---- put ice quantities onto exchange grid ----
!---- (assume that ocean quantites are stored in no ice partition) ----
!     (note: ex_avail is true at ice and ocean points)

   ex_avail = .false.
   call put_exchange_grid  &
             (Ice%t_surf       , ex_t_surf     , bd_map_ice_top, &
                                                 avail=ex_avail)
   call put_exchange_grid  &
             (Ice%rough_mom    , ex_rough_mom  , bd_map_ice_top)
   call put_exchange_grid  &
             (Ice%rough_heat   , ex_rough_heat , bd_map_ice_top)
   call put_exchange_grid  &
             (Ice%rough_moist  , ex_rough_moist, bd_map_ice_top)
   call put_exchange_grid  &
             (Ice%albedo       , ex_albedo     , bd_map_ice_top)
   call put_exchange_grid  &
             (Ice%u_surf       , ex_u_surf     , bd_map_ice_top)
   call put_exchange_grid  &
             (Ice%v_surf       , ex_v_surf     , bd_map_ice_top)

!=======================================================================
!---- compute explicit fluxes and tendencies at all available points ---

   ex_avail = ex_land .or. ex_avail

   call surface_flux (ex_t_atm, ex_q_atm, ex_u_atm, ex_v_atm,  &
                      ex_p_atm, ex_z_atm,                      &
                  ex_p_surf, ex_t_surf, ex_u_surf, ex_v_surf,  &
                  ex_rough_mom, ex_rough_heat, ex_rough_moist, &
                  ex_gust, ex_stomatal,                        &
                  ex_snow, ex_water,  ex_max_water,            &
             ex_flux_t, ex_flux_q, ex_flux_lw, ex_flux_u, ex_flux_v,  &
             ex_cd_m,   ex_cd_t, ex_cd_q, ex_wind,                    &
             ex_u_star, ex_b_star, ex_q_star, ex_q_surf,              &
             ex_dhdt_surf, ex_dedt_surf,  ex_drdt_surf,               &
             ex_dhdt_atm,  ex_dedq_atm,   ex_dtaudv_atm,              &
             dt,           ex_land,       ex_glacier,      ex_avail   )

!=======================================================================
!---- get mean quantities on atmosphere grid ----
!---- compute t surf for radiation ----

   ex_t_surf4 = ex_t_surf ** 4

   call get_exchange_grid (ex_t_surf4  , t_surf_atm   , bd_map_atm)
   call get_exchange_grid (ex_albedo   , albedo_atm   , bd_map_atm)
   call get_exchange_grid (ex_rough_mom, rough_mom_atm, bd_map_atm)
   call get_exchange_grid (ex_land_frac, land_frac_atm, bd_map_atm)

   call get_exchange_grid (ex_flux_u,     flux_u_atm, bd_map_atm)
   call get_exchange_grid (ex_flux_v,     flux_v_atm, bd_map_atm)
   call get_exchange_grid (ex_dtaudv_atm, dtaudv_atm, bd_map_atm)
   call get_exchange_grid (ex_u_star,     u_star_atm, bd_map_atm)
   call get_exchange_grid (ex_b_star,     b_star_atm, bd_map_atm)

   t_surf_atm = t_surf_atm ** 0.25
   
!---- save atmos albedo fix on exch grid ----

   allocate ( ex_albedo_fix(ex_num_top) )
   call put_exchange_grid (albedo_atm, ex_albedo_fix, bd_map_atm)
   ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)

!=======================================================================
!-------------------- diagnostics section ------------------------------

!------- save static fields first time only ------
   if (first_static) then

!------- land fraction ------
      if ( id_land_mask > 0 ) then
         used = send_data ( id_land_mask, land_frac_atm, Time )
      endif

!------- glacier fraction -----
      if ( id_glac_mask > 0 ) then
         call get_exchange_grid (ex_glac_frac, diag_atm, bd_map_atm)
         used = send_data ( id_glac_mask, diag_atm, Time )
      endif

      first_static = .false.
   endif

!------- drag coeff moisture -----------
   if ( id_wind > 0 ) then
      call get_exchange_grid (ex_wind, diag_atm, bd_map_atm)
      used = send_data ( id_wind, diag_atm, Time )
   endif
!------- drag coeff moisture -----------
   if ( id_drag_moist > 0 ) then
      call get_exchange_grid (ex_cd_q, diag_atm, bd_map_atm)
      used = send_data ( id_drag_moist, diag_atm, Time )
   endif

!------- drag coeff heat -----------
   if ( id_drag_heat > 0 ) then
      call get_exchange_grid (ex_cd_t, diag_atm, bd_map_atm)
      used = send_data ( id_drag_heat, diag_atm, Time )
   endif

!------- drag coeff momemtum -----------
   if ( id_drag_mom > 0 ) then
      call get_exchange_grid (ex_cd_m, diag_atm, bd_map_atm)
      used = send_data ( id_drag_mom, diag_atm, Time )
   endif

!------- roughness moisture -----------
   if ( id_rough_moist > 0 ) then
      call get_exchange_grid (ex_rough_moist, diag_atm, bd_map_atm)
      used = send_data ( id_rough_moist, diag_atm, Time )
   endif

!------- roughness heat -----------
   if ( id_rough_heat > 0 ) then
      call get_exchange_grid (ex_rough_heat, diag_atm, bd_map_atm)
      used = send_data ( id_rough_heat, diag_atm, Time )
   endif

!------- roughness momemtum -----------
   if ( id_rough_mom > 0 ) then
      used = send_data ( id_rough_mom, rough_mom_atm, Time )
   endif

!------- friction velocity -----------
   if ( id_u_star > 0 ) then
      used = send_data ( id_u_star, u_star_atm, Time )
   endif

!------- bouyancy -----------
   if ( id_b_star > 0 ) then
      used = send_data ( id_b_star, b_star_atm, Time )
   endif

!-----------------------------------------------------------------------
!------ diagnostics for fields at bottom atmospheric level ------

   if ( id_t_atm > 0 ) then
      call get_exchange_grid (ex_t_atm, diag_atm, bd_map_atm)
      used = send_data ( id_t_atm, diag_atm, Time )
   endif

   if ( id_u_atm > 0 ) then
      call get_exchange_grid (ex_u_atm, diag_atm, bd_map_atm)
      used = send_data ( id_u_atm, diag_atm, Time )
   endif

   if ( id_v_atm > 0 ) then
      call get_exchange_grid (ex_v_atm, diag_atm, bd_map_atm)
      used = send_data ( id_v_atm, diag_atm, Time )
   endif

!-----------------------------------------------------------------------
!--------- diagnostics for fields at reference level ---------

   if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
        id_u_ref > 0 .or. id_v_ref  > 0 ) then

       zrefm = z_ref_mom
       zrefh = z_ref_heat
!      ---- optimize calculation ----
       if ( id_t_ref <= 0 ) zrefh = zrefm

       call surface_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
                              ex_rough_heat, ex_rough_moist,          &
                              ex_u_star, ex_b_star, ex_q_star,        &
                              ex_del_m, ex_del_h, ex_del_q, ex_avail  )

!    ------- reference relative humidity -----------
       if ( id_rh_ref > 0 ) then
          ex_ref   = ex_q_surf + (ex_q_atm-ex_q_surf) * ex_del_q
          ex_t_ref = ex_t_surf + (ex_t_atm-ex_t_surf) * ex_del_h
          call escomp (ex_t_ref, ex_qs_ref)
          ex_qs_ref = d622*ex_qs_ref/(ex_p_surf-d378*ex_qs_ref)
          ex_ref    = 100.*ex_ref/ex_qs_ref
          call get_exchange_grid (ex_ref, diag_atm, bd_map_atm)
          used = send_data ( id_rh_ref, diag_atm, Time )
       endif

!    ------- reference temp -----------
       if ( id_t_ref > 0 ) then
          ex_ref = ex_t_surf + (ex_t_atm-ex_t_surf) * ex_del_h
          call get_exchange_grid (ex_ref, diag_atm, bd_map_atm)
          used = send_data ( id_t_ref, diag_atm, Time )
       endif

!    ------- reference u comp -----------
       if ( id_u_ref > 0 ) then
          ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
          call get_exchange_grid (ex_ref, diag_atm, bd_map_atm)
          used = send_data ( id_u_ref, diag_atm, Time )
       endif

!    ------- reference v comp -----------
       if ( id_v_ref > 0 ) then
          ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
          call get_exchange_grid (ex_ref, diag_atm, bd_map_atm)
          used = send_data ( id_v_ref, diag_atm, Time )
       endif

!    ------- interp factor for heat ------
       if ( id_del_h > 0 ) then
          call get_exchange_grid (ex_del_h, diag_atm, bd_map_atm)
          used = send_data ( id_del_h, diag_atm, Time )
       endif

!    ------- interp factor for momentum ------
       if ( id_del_m > 0 ) then
          call get_exchange_grid (ex_del_m, diag_atm, bd_map_atm)
          used = send_data ( id_del_m, diag_atm, Time )
       endif

!    ------- interp factor for moisture ------
       if ( id_del_q > 0 ) then
          call get_exchange_grid (ex_del_q, diag_atm, bd_map_atm)
          used = send_data ( id_del_q, diag_atm, Time )
       endif

   endif

!=======================================================================

 end subroutine flux_calculation

!#######################################################################

 subroutine flux_down_from_atmos (Time, Atm, Land, Ice,             &
                                  flux_u_atm, flux_v_atm,           &
                                  flux_t_land, flux_q_land,         &
                                  flux_lw_land, flux_sw_land,       &
                                  dhdt_land, dedt_land, drdt_land,  &
                                  lprec_land, fprec_land,           &
                                  flux_t_ice, flux_q_ice,           &
                                  flux_lw_ice, flux_sw_ice,         &
                                  dhdt_ice, dedt_ice, drdt_ice,     &
                                  lprec_ice , fprec_ice,            &
                                  flux_u_ice, flux_v_ice, coszen_ice)

 type       (time_type), intent(in)  :: Time
 type (atmos_boundary_data_type), intent(in)  :: Atm
 type  (land_boundary_data_type), intent(in)  :: Land
 type   (ice_boundary_data_type), intent(in)  :: Ice
 real, dimension(:,:),   intent(in)  :: flux_u_atm, flux_v_atm
 real, dimension(:,:,:), intent(out) ::                               &
                                    flux_t_land, flux_q_land,         &
                                    flux_lw_land, flux_sw_land,       &
                                    dhdt_land, dedt_land, drdt_land,  &
                                    lprec_land, fprec_land,           &
                                    flux_t_ice, flux_q_ice,           &
                                    flux_lw_ice, flux_sw_ice,         &
                                    dhdt_ice, dedt_ice, drdt_ice,     &
                                    lprec_ice , fprec_ice,            &
                                    flux_u_ice, flux_v_ice, coszen_ice

 real, dimension(ex_num_top) :: ex_flux_sw, ex_flux_lwd, &
                                ex_lprec, ex_fprec,      &
                                ex_flux_u, ex_flux_v,    &
                                ex_coszen, ex_ice_frac

 real, dimension(ex_num_top) :: ex_gamma  , ex_dtmass,  &
                                ex_delta_t, ex_delta_q, &
                                ex_dflux_t, ex_dflux_q

 real :: ice_frac (size(dhdt_ice,1),size(dhdt_ice,2),size(dhdt_ice,3))
 real :: diag_atm (size(flux_u_atm,1),size(flux_u_atm,2))

 real    :: cp_inv
 logical :: used

!-----------------------------------------------------------------------
!---- put atmosphere quantities onto exchange grid ----

   call put_exchange_grid (Atm%flux_sw,  ex_flux_sw,  bd_map_atm)
   call put_exchange_grid (Atm%flux_lw,  ex_flux_lwd, bd_map_atm)

   call put_exchange_grid (Atm%lprec, ex_lprec, bd_map_atm)
   call put_exchange_grid (Atm%fprec, ex_fprec, bd_map_atm)

   call put_exchange_grid (Atm%coszen, ex_coszen, bd_map_atm)

   call put_exchange_grid (flux_u_atm, ex_flux_u, bd_map_atm)
   call put_exchange_grid (flux_v_atm, ex_flux_v, bd_map_atm)

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----
   
   ex_flux_sw = ex_flux_sw * ex_albedo_fix
   deallocate ( ex_albedo_fix )

!----- compute net longwave flux (down-up) -----
!     (note: lw up already in ex_flux_lw)

   ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----

 
   call put_exchange_grid (Atm%Surf_Diff%dtmass   , ex_dtmass   , bd_map_atm)
   call put_exchange_grid (Atm%Surf_Diff%delta_t  , ex_delta_t  , bd_map_atm)
   call put_exchange_grid (Atm%Surf_Diff%delta_q  , ex_delta_q  , bd_map_atm)
   call put_exchange_grid (Atm%Surf_Diff%dflux_t  , ex_dflux_t  , bd_map_atm)
   call put_exchange_grid (Atm%Surf_Diff%dflux_q  , ex_dflux_q  , bd_map_atm)

 cp_inv = 1.0/cp

 where(ex_avail)

! temperature

   ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
   ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
   ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    

   ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
   ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

! moisture

   ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))
   ex_e_q_n      =  ex_dtmass*ex_dedt_surf*ex_gamma
   ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma    

   ex_flux_q     =  ex_flux_q        + ex_dedq_atm * ex_f_q_delt_n 
   ex_dedt_surf  =  ex_dedt_surf     + ex_dedq_atm * ex_e_q_n   

 endwhere


!-----------------------------------------------------------------------
!---- output fields on the land grid -------

    call get_exchange_grid (ex_flux_t,  flux_t_land,  bd_map_land)
    call get_exchange_grid (ex_flux_q,  flux_q_land,  bd_map_land)
    call get_exchange_grid (ex_flux_sw, flux_sw_land, bd_map_land)
    call get_exchange_grid (ex_flux_lw, flux_lw_land, bd_map_land)
    call get_exchange_grid (ex_dhdt_surf, dhdt_land,  bd_map_land)
    call get_exchange_grid (ex_dedt_surf, dedt_land,  bd_map_land)
    call get_exchange_grid (ex_drdt_surf, drdt_land,  bd_map_land)
    call get_exchange_grid (ex_lprec,    lprec_land,  bd_map_land)
    call get_exchange_grid (ex_fprec,    fprec_land,  bd_map_land)

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

   call get_exchange_grid (ex_flux_t,  flux_t_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_flux_q,  flux_q_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_flux_sw, flux_sw_ice, bd_map_ice_top)
   call get_exchange_grid (ex_flux_lw, flux_lw_ice, bd_map_ice_top)
   call get_exchange_grid (ex_dhdt_surf, dhdt_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_dedt_surf, dedt_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_drdt_surf, drdt_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_lprec,    lprec_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_fprec,    fprec_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_flux_u,  flux_u_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_flux_v,  flux_v_ice,  bd_map_ice_top)
   call get_exchange_grid (ex_coszen,  coszen_ice,  bd_map_ice_top)

!=======================================================================
!-------------------- diagnostics section ------------------------------

!------- zonal wind stress -----------
   if ( id_u_flux > 0 ) then
      used = send_data ( id_u_flux, flux_u_atm, Time )
   endif

!------- meridional wind stress -----------
   if ( id_v_flux > 0 ) then
      used = send_data ( id_v_flux, flux_v_atm, Time )
   endif

!------- diagnostics for sea ice fraction -----------
!---- this should probably be done in slow loop -----

   if ( id_ice_mask > 0 ) then
      ice_frac        = 1.
      ice_frac(:,:,1) = 0.
      ex_ice_frac     = 0.
      call put_exchange_grid (ice_frac, ex_ice_frac, bd_map_ice_top)
      call get_exchange_grid (ex_ice_frac, diag_atm, bd_map_atm)
      used = send_data ( id_ice_mask, diag_atm, Time )
   endif

!=======================================================================

 end subroutine flux_down_from_atmos

!#######################################################################

 subroutine flux_ice_to_ocean ( Ice, flux_u_ocean,  flux_v_ocean, &
                                     flux_t_ocean,  flux_q_ocean, &
                                    flux_sw_ocean, flux_lw_ocean, &
                                      lprec_ocean,   fprec_ocean, &
                                     runoff_ocean                 )

  type (ice_boundary_data_type),   intent(in)  :: Ice
  real, dimension(:,:),   intent(out) :: flux_u_ocean,  flux_v_ocean,  &
                                         flux_t_ocean,  flux_q_ocean,  &
                                         flux_sw_ocean, flux_lw_ocean, &
                                         lprec_ocean ,  fprec_ocean,   &
                                         runoff_ocean

   real, dimension(ex_num_bot) :: ex_flux_t, ex_flux_q, &
                                  ex_flux_sw, ex_flux_lw, &
                                  ex_lprec,   ex_fprec, ex_runoff

   real, dimension(ex_num_bot_uv) :: ex_flux_u, ex_flux_v

!-----------------------------------------------------------------------
!-----  put ice grid quantities onto exchange grid -----
!-----  then get ocean grid fields from exchange grid -----

  call put_exchange_grid (Ice%flux_u,  ex_flux_u , bd_map_ice_bot_uv)
  call put_exchange_grid (Ice%flux_v,  ex_flux_v , bd_map_ice_bot_uv)
  call put_exchange_grid (Ice%flux_t,  ex_flux_t , bd_map_ice_bot)
  call put_exchange_grid (Ice%flux_q,  ex_flux_q , bd_map_ice_bot)
  call put_exchange_grid (Ice%flux_sw, ex_flux_sw, bd_map_ice_bot)
  call put_exchange_grid (Ice%flux_lw, ex_flux_lw, bd_map_ice_bot)
  call put_exchange_grid (Ice%lprec,   ex_lprec,   bd_map_ice_bot)
  call put_exchange_grid (Ice%fprec,   ex_fprec,   bd_map_ice_bot)
  call put_exchange_grid (Ice%runoff,  ex_runoff,  bd_map_ice_bot)

  call get_exchange_grid (ex_flux_u ,  flux_u_ocean, bd_map_ocean_model_uv)
  call get_exchange_grid (ex_flux_v ,  flux_v_ocean, bd_map_ocean_model_uv)
  call get_exchange_grid (ex_flux_t ,  flux_t_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_flux_q ,  flux_q_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_flux_sw, flux_sw_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_flux_lw, flux_lw_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_lprec,     lprec_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_fprec,     fprec_ocean, bd_map_ocean_model)
  call get_exchange_grid (ex_runoff,   runoff_ocean, bd_map_ocean_model)

!-----------------------------------------------------------------------

 end subroutine flux_ice_to_ocean

!#######################################################################

 subroutine flux_ocean_to_ice ( Ocean, Ice, t_surf_ice, frazil_ice, &
                                            u_surf_ice, v_surf_ice  )

  type (ocean_boundary_data_type), intent(in)  :: Ocean
  type (ice_boundary_data_type),   intent(in)  :: Ice
  real, dimension(:,:,:), intent(out) :: t_surf_ice,  frazil_ice, &
                                         u_surf_ice,  v_surf_ice

   real, dimension(ex_num_bot)    :: ex_t_surf, ex_frazil
   real, dimension(ex_num_bot_uv) :: ex_u_surf, ex_v_surf

!-----------------------------------------------------------------------
!-----  put ocean grid fields onto exchange grid -----
!-----  then get ice grid fields from exchange grid data -----
!-----  use new fraction ice areas ------

 call set_frac_area (Ice%part_size,    bd_map_ice_top   )
 call set_frac_area (Ice%part_size   , bd_map_ice_bot   )
 call set_frac_area (Ice%part_size_uv, bd_map_ice_bot_uv)

 call put_exchange_grid (Ocean%t_surf_data, ex_t_surf, bd_map_ocean_data)
 ex_frazil = 0.0
 ex_u_surf = 0.0
 ex_v_surf = 0.0

 call put_exchange_grid (Ocean%t_surf, ex_t_surf, bd_map_ocean_model)
 call put_exchange_grid (Ocean%frazil, ex_frazil, bd_map_ocean_model)
 call put_exchange_grid (Ocean%u_surf, ex_u_surf, bd_map_ocean_model_uv)
 call put_exchange_grid (Ocean%v_surf, ex_v_surf, bd_map_ocean_model_uv)

 call get_exchange_grid (ex_t_surf,    t_surf_ice, bd_map_ice_bot)
 call get_exchange_grid (ex_frazil,    frazil_ice, bd_map_ice_bot)
 call get_exchange_grid (ex_u_surf,    u_surf_ice, bd_map_ice_bot_uv)
 call get_exchange_grid (ex_v_surf,    v_surf_ice, bd_map_ice_bot_uv)

!-----------------------------------------------------------------------

 end subroutine flux_ocean_to_ice

!#######################################################################

 subroutine flux_up_to_atmos ( Time, Land, Ice, dt_t_atm, dt_q_atm )

 type       (time_type), intent(in)  :: Time
 type  (land_boundary_data_type), intent(in)  :: Land
 type   (ice_boundary_data_type), intent(in)  :: Ice
 real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm

  real, dimension(ex_num_top) :: ex_t_surf_new, ex_dt_t_surf,  &
                                 ex_dt_t, ex_dt_q, &
                                 ex_delta_t_n, ex_delta_q_n
  real, dimension(size(dt_t_atm,1),size(dt_t_atm,2)) :: diag_atm, &
                                                        evap_atm
  logical :: used
!-----------------------------------------------------------------------
!----- compute surface temperature change -----

   ex_t_surf_new = 200.
   call put_exchange_grid (Land%t_surf, ex_t_surf_new, bd_map_land)
   call put_exchange_grid (Ice%t_surf, ex_t_surf_new, bd_map_ice_top)

   ex_dt_t_surf = ex_t_surf_new - ex_t_surf

!-----------------------------------------------------------------------
!-----  adjust fluxes and atmospheric increments for 
!-----  implicit dependence on surface temperature -----

   ex_delta_t_n = 0.0
   ex_delta_q_n = 0.0

   where(ex_avail)
     ex_flux_t     = ex_flux_t  + ex_dt_t_surf * ex_dhdt_surf
     ex_flux_q     = ex_flux_q  + ex_dt_t_surf * ex_dedt_surf
     ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
     ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_surf*ex_e_t_n
     ex_delta_q_n  = ex_f_q_delt_n  + ex_dt_t_surf*ex_e_q_n
   endwhere 

!-----------------------------------------------------------------------
!---- get mean quantites on atmospheric grid ----

   call get_exchange_grid (ex_delta_t_n, dt_t_atm, bd_map_atm)
   call get_exchange_grid (ex_delta_q_n, dt_q_atm, bd_map_atm)

!  ---- always get evaporation for diagnostic purposes ----

   call get_exchange_grid (ex_flux_q, evap_atm, bd_map_atm)

!=======================================================================
!-------------------- diagnostics section ------------------------------

!------- new surface temperature -----------
   if ( id_t_surf > 0 ) then
      call get_exchange_grid (ex_t_surf_new, diag_atm, bd_map_atm)
      used = send_data ( id_t_surf, diag_atm, Time )
   endif

!------- sensible heat flux -----------
   if ( id_t_flux > 0 ) then
      call get_exchange_grid (ex_flux_t, diag_atm, bd_map_atm)
      used = send_data ( id_t_flux, diag_atm, Time )
   endif

!------- latent heat flux (see below) -----------

!  if ( id_q_flux > 0 ) then
!!!!  latent  = hlv
!!!!  where (ice_or_snow) latent = hlf + hlv
!     call get_exchange_grid (ex_flux_q, diag_atm, bd_map_atm)
!     used = send_data ( id_q_flux, diag_atm, Time )
!  endif

!------- net longwave flux -----------
   if ( id_r_flux > 0 ) then
      call get_exchange_grid (ex_flux_lw, diag_atm, bd_map_atm)
      used = send_data ( id_r_flux, diag_atm, Time )
   endif

!------- evaporation rate -----------

   if ( id_q_flux > 0 ) then
      used = send_data ( id_q_flux, evap_atm, Time )
   endif

!-----------------------------------------------------------------------
!---- accumulate global integral of evaporation (mm/day) -----

    call sum_diag_integral_field ('evap', evap_atm*86400.)

!=======================================================================
!---- deallocate module storage ----

   deallocate ( ex_t_surf, ex_dhdt_surf, ex_dedt_surf, ex_drdt_surf, &
                ex_dhdt_atm,  ex_dedq_atm,                           &
                ex_flux_t, ex_flux_q, ex_flux_lw, ex_avail           )

   deallocate ( ex_f_t_delt_n, ex_f_q_delt_n, ex_e_t_n, ex_e_q_n)


!-----------------------------------------------------------------------

 end subroutine flux_up_to_atmos

!#######################################################################

 subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean )

 type                (time_type), intent(in)  :: Time
 type (atmos_boundary_data_type), intent(in)  :: Atm
 type  (land_boundary_data_type), intent(in)  :: Land
 type   (ice_boundary_data_type), intent(in)  :: Ice
 type (ocean_boundary_data_type), intent(in)  :: Ocean

 logical, dimension(size(Atm%lon_bnd)-1,size(Atm%lat_bnd)-1) :: atm_mask
 integer :: part_land, part_ice
 integer :: unit, ierr, io, iref

!-----------------------------------------------------------------------
!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=flux_exchange_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'flux_exchange_nml')
      enddo
 10   call close_file (unit)
   endif

!--------- write version number and namelist ------------------

   unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) then
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (unit, nml=flux_exchange_nml)
   endif
   call close_file (unit)

!-----------------------------------------------------------------------
!---- exchange map between atmosphere, land, and ice top -----

   atm_mask  = .true.
   part_land = size(Land%mask,3)
   part_ice  = size (Ice%ice_mask,3)

   call init_boundary_map (bd_map_atm    , 1, ex_map_top)
   call init_boundary_map (bd_map_land   , 2, ex_map_top)
   call init_boundary_map (bd_map_ice_top, 2, ex_map_top)

   call lon_lat_size ( Atm%lon_bnd,  Atm%lat_bnd,  &
                      Land%lon_bnd, Land%lat_bnd,  &
                  atm_mask, Land%mask(:,:,1), 1, part_land, bd_map_atm)

   call lon_lat_size (Atm%lon_bnd, Atm%lat_bnd,  &
                      Ice%lon_bnd, Ice%lat_bnd,  &
                    atm_mask, Ice%mask,   1, part_ice, bd_map_atm)

   call lon_lat_map ( Atm%lon_bnd,  Atm%lat_bnd,     &
                     Land%lon_bnd, Land%lat_bnd,     &
                      atm_mask,    Land%mask(:,:,1), &
                         1,           part_land,     &
                      bd_map_atm ,  bd_map_land      )

   call lon_lat_map ( Atm%lon_bnd,  Atm%lat_bnd,     &
                      Ice%lon_bnd,  Ice%lat_bnd,     &
                      atm_mask,     Ice%mask,        &
                         1,           part_ice ,     &
                      bd_map_atm ,  bd_map_ice_top   )

   call complete_side1_boundary_map (bd_map_atm)

   call complete_side2_boundary_map (bd_map_land, &
                     size(Land%mask,1), size(Land%mask,2), part_land )

   call complete_side2_boundary_map (bd_map_ice_top, &
                     size(Ice%mask,1),  size(Ice%mask,2),  part_ice )

!-----------------------------------------------------------------------
!---- exchange map between ocean and ice bottom -----
!---- at mass (temperature) points ----

   call init_boundary_map (bd_map_ice_bot    , 1, ex_map_bot)
   call init_boundary_map (bd_map_ocean_data , 2, ex_map_bot)
   call init_boundary_map (bd_map_ocean_model, 2, ex_map_bot)

   call lon_lat_size (       Ice%lon_bnd,          Ice%lat_bnd,  &
                      Ocean%Data%lon_bnd,   Ocean%Data%lat_bnd,  &
                                Ice%mask,   Ocean%Data%mask,     &
                       part_ice,        1,        bd_map_ice_bot )

   call lon_lat_size (        Ice%lon_bnd,         Ice%lat_bnd,  &
                      Ocean%Model%lon_bnd, Ocean%Model%lat_bnd,  &
                                Ice%mask,  Ocean%Model%mask,     &
                       part_ice,        1,        bd_map_ice_bot )

   call lon_lat_map (       Ice%lon_bnd,          Ice%lat_bnd,  &
                     Ocean%Data%lon_bnd,   Ocean%Data%lat_bnd,  &
                               Ice%mask,   Ocean%Data%mask,     &
                     part_ice, 1,  bd_map_ice_bot, bd_map_ocean_data )

   call lon_lat_map (        Ice%lon_bnd,          Ice%lat_bnd,  &
                     Ocean%Model%lon_bnd,  Ocean%Model%lat_bnd,  &
                               Ice%mask,   Ocean%Model%mask,     &
                     part_ice, 1,  bd_map_ice_bot, bd_map_ocean_model )

   call complete_side1_boundary_map (bd_map_ice_bot)

   call complete_side2_boundary_map (bd_map_ocean_data, &
                  size(Ocean%Data%mask,1), size(Ocean%Data%mask,2), 1 )

   call complete_side2_boundary_map (bd_map_ocean_model, &
                  size(Ocean%Model%mask,1), size(Ocean%Model%mask,2), 1 )

!---- at velocity points ----

   call init_boundary_map (bd_map_ice_bot_uv    , 1, ex_map_bot_uv)
   call init_boundary_map (bd_map_ocean_data_uv , 2, ex_map_bot_uv)
   call init_boundary_map (bd_map_ocean_model_uv, 2, ex_map_bot_uv)

   call lon_lat_size (       Ice%lon_bnd_uv,         Ice%lat_bnd_uv,  &
                      Ocean%Data%lon_bnd_uv,  Ocean%Data%lat_bnd_uv,  &
                                Ice%mask_uv,  Ocean%Data%mask_uv,     &
                         part_ice,        1,        bd_map_ice_bot_uv )

   call lon_lat_size (        Ice%lon_bnd_uv,         Ice%lat_bnd_uv,  &
                      Ocean%Model%lon_bnd_uv, Ocean%Model%lat_bnd_uv,  &
                                Ice%mask_uv,  Ocean%Model%mask_uv,     &
                         part_ice,        1,        bd_map_ice_bot_uv )

   call lon_lat_map (       Ice%lon_bnd_uv,          Ice%lat_bnd_uv,  &
                     Ocean%Data%lon_bnd_uv,   Ocean%Data%lat_bnd_uv,  &
                               Ice%mask_uv,   Ocean%Data%mask_uv,     &
                     part_ice, 1,  bd_map_ice_bot_uv, bd_map_ocean_data_uv )

   call lon_lat_map (        Ice%lon_bnd_uv,          Ice%lat_bnd_uv,  &
                     Ocean%Model%lon_bnd_uv,  Ocean%Model%lat_bnd_uv,  &
                               Ice%mask_uv,   Ocean%Model%mask_uv,     &
                     part_ice, 1,  bd_map_ice_bot_uv, bd_map_ocean_model_uv )

   call complete_side1_boundary_map (bd_map_ice_bot_uv)

   call complete_side2_boundary_map (bd_map_ocean_data_uv, &
          size(Ocean%Data%mask_uv,1), size(Ocean%Data%mask_uv,2), 1 )

   call complete_side2_boundary_map (bd_map_ocean_model_uv, &
          size(Ocean%Model%mask_uv,1), size(Ocean%Model%mask_uv,2), 1 )

!-----------------------------------------------------------------------

   ex_num_top    = get_exchange_grid_size (bd_map_atm)
   ex_num_bot    = get_exchange_grid_size (bd_map_ocean_data)
   ex_num_bot_uv = get_exchange_grid_size (bd_map_ocean_data_uv)

!-----------------------------------------------------------------------
!----- initialize quantities for global integral package -----

!! call diag_integral_field_init ('prec', 'f6.3')
   call diag_integral_field_init ('evap', 'f6.3')

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

   call diag_field_init ( Time, Atm%axes(1:2) )
!---- done ----

   do_init = .false.

!-----------------------------------------------------------------------

 end subroutine flux_exchange_init

!#######################################################################

 subroutine put_logical_to_real (mask, ex_mask, bd_map)

   logical                 , intent(in)    :: mask(:,:,:)
   real                    , intent(inout) :: ex_mask(:)
   type (boundary_map_type), intent(inout) :: bd_map

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

   call put_exchange_grid (rmask, ex_mask, bd_map)

 end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
                        vrange = (/ -400., 400. /), &
                        frange = (/ -0.01, 1.01 /)
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

   id_glac_mask = &
   register_static_field ( mod_name, 'glac_mask', atmos_axes,  &
                          'fractional amount of glacier', 'none', &
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

   id_t_flux     = &
   register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
                        'sensible heat flux',     'w/m2'    )

   id_q_flux     = &
   register_diag_field ( mod_name, 'evap',       atmos_axes, Time, &
                        'evaporation rate',        'kg/m2/s'  )

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

!-----------------------------------------------------------------------

end subroutine diag_field_init

!#######################################################################

end module flux_exchange_mod

