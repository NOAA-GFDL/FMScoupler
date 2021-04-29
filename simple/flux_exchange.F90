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
!
module flux_exchange_mod

!-----------------------------------------------------------------------
!! Components
use   atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
use    land_model_mod, only:  land_data_type, atmos_land_boundary_type
use     ice_model_mod, only:   ice_data_type, atmos_ice_boundary_type
#ifndef use_AM3_physics
use atmos_cmip_diag_mod,   only: register_cmip_diag_field_2d
#endif
use surface_flux_mod, only: surface_flux, surface_flux_init

!! FMS
use FMS
use FMSconstants, only: RDGAS, RVGAS, CP_AIR, HLV, HLF, PI

implicit none
private

public :: flux_exchange_init,   &
          sfc_boundary_layer,   &
          flux_down_from_atmos, &
          flux_up_to_atmos,     &
          flux_exchange_end

!-----------------------------------------------------------------------
character(len=128) :: version = '$Id$'
character(len=128) :: tag = '$Name$'

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=14), parameter :: mod_name = 'flux'

integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,              &
           id_rough_moist, id_rough_heat, id_rough_mom,             &
           id_u_star, id_b_star, id_q_star, id_u_flux, id_v_flux,   &
           id_t_surf, id_t_flux, id_q_flux, id_r_flux,              &
           id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
           id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,  id_q_ref,    &
           id_del_h,  id_del_m,  id_del_q, id_albedo,  id_gust,     &
           id_t_ca,   id_q_surf, id_q_atm, id_z_atm, id_p_atm,      &
           id_land_mask, id_ice_mask, id_rough_scale,               &
           id_albedo_vis_dir, id_albedo_nir_dir,                    &
           id_albedo_vis_dif, id_albedo_nir_dif

! t_ref(tas), u_ref, v_ref, t_surf, id_wind(check)
! v_flux (wind stress? check)
! q_ref (huss), t_flux (hfss)

! Atm%slp can be saved as psl

! probably don't need id_tos, id_tslsi but check data request

  ! lgs id's for cmip specific fields for aquaplanet
  integer :: id_tas, id_uas, id_vas, id_ts, id_psl, &
             id_sfcWind, id_tauu, id_tauv, &
             id_hurs, id_huss, id_evspsbl, id_hfls, id_hfss, &
             !id_sftlf, id_tos, id_tslsi, id_sic, &
             id_height2m, id_height10m

logical :: first_static = .true.
logical :: do_init = .true.
logical :: do_read_nml = .true.

! index in tracer array for water vapor
integer :: isphum, n_atm_tr_tot, n_atm_tr

!-----------------------------------------------------------------------

real ::  z_ref_heat = 2.,  &
         z_ref_mom  = 10.

logical :: use_existing_grid_spec = .false.
logical :: all_ocean = .true.
logical :: all_land = .false. ! note: if both all_ocean = all_land = .true.
                              !       then all_ocean = .true. will be used

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom,  &
                             use_existing_grid_spec, &
                             all_ocean, all_land

!-----------------------------------------------------------------------

!---- grid indices ----

integer :: is, ie, js, je

!---- allocatable module storage ------

  real, allocatable, dimension(:,:) :: t_surf, t_ca, q_surf, p_surf

  real, allocatable, dimension(:,:) :: e_t_n, f_t_delt_n, &
                                       e_q_n, f_q_delt_n

  real, allocatable, dimension(:,:) :: dhdt_surf, dedt_surf, dedq_surf, &
                                       drdt_surf, dhdt_atm,  dedq_atm,  &
                                       flux_t, flux_q, flux_lw

  real, allocatable, dimension(:,:) :: flux_u, flux_v, drag_q, &
                                       dtaudu_atm, dtaudv_atm

  real, allocatable, dimension(:,:) :: cd_t, cd_m, b_star, u_star, wind

!-----------------------------------------------------------------------

 real, parameter :: d622 = RDGAS/RVGAS
 real, parameter :: d378 = 1.0-d622

 logical :: used

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )

 real,                   intent(in)  :: dt  !< Time step
 type       (time_type), intent(in)  :: Time !< Current time
 type (atmos_data_type), intent(in)  :: Atm  !< A derived data type to specify atmospheric boundary data
 type(land_data_type),  intent(inout)  :: Land !< A derived data type to specify land boundary data
 type(ice_data_type),   intent(inout)  :: Ice !< A derived data type to specify ice boundary data
 type(land_ice_atmos_boundary_type), intent(inout) :: Boundary !< A derived data type to specify properties and
 !! fluxes passed from exchange grid to the atmosphere,

real, dimension(is:ie,js:je) :: u_surf, v_surf, rough_heat, rough_moist, &
                                rough_mom, rough_scale, q_star, cd_q,    &
                                albedo, albedo_vis_dir, albedo_nir_dir,  &
                                albedo_vis_dif, albedo_nir_dif,          &
                                del_m, del_h, del_q, land_frac,          &
                                ref, ref2, t_ref, qs_ref, qs_ref_cmip

logical, dimension(is:ie,js:je) :: mask, seawater, avail
real :: zrefm, zrefh


!-----------------------------------------------------------------------

   if (do_init) call error_mesg ('sfc_boundary_layer',  &
                 'must call simple_surface_init first', FATAL)

!-----------------------------------------------------------------------
!------ allocate storage also needed in flux_up_to_atmos -----

   allocate ( e_t_n      (is:ie, js:je), &
              e_q_n      (is:ie, js:je), &
              f_t_delt_n (is:ie, js:je), &
              f_q_delt_n (is:ie, js:je), &
              dhdt_surf  (is:ie, js:je), &
              dedt_surf  (is:ie, js:je), &
              dedq_surf  (is:ie, js:je), &
              drdt_surf  (is:ie, js:je), &
              dhdt_atm   (is:ie, js:je), &
              dedq_atm   (is:ie, js:je), &
              flux_t     (is:ie, js:je), &
              flux_q     (is:ie, js:je), &
              flux_lw    (is:ie, js:je), &
              flux_u     (is:ie, js:je), &
              flux_v     (is:ie, js:je), &
              dtaudu_atm (is:ie, js:je), &
              dtaudv_atm (is:ie, js:je), &
              drag_q     (is:ie, js:je), &
              t_surf     (is:ie, js:je), &
              t_ca       (is:ie, js:je), &
              p_surf     (is:ie, js:je), &
              q_surf     (is:ie, js:je)  )

   allocate ( cd_t       (is:ie, js:je), &
              cd_m       (is:ie, js:je), &
              b_star     (is:ie, js:je), &
              u_star     (is:ie, js:je), &
              wind       (is:ie, js:je)  )


   u_surf     = 0.0
   v_surf     = 0.0

!---- do not use if relax time /= 0 ----
   cd_t = 0.0
   cd_m = 0.0
   cd_q = 0.0

   avail   = .true.

!---- atmosphere quantities ----

   p_surf = Atm%p_surf

!---- ice quantities ----

 where (Ice%mask)
   t_surf = Ice%t_surf
   t_ca   = Ice%t_surf ! to define values over the ice/ocean
   rough_mom   = Ice%rough_mom
   rough_heat  = Ice%rough_heat
   rough_moist = Ice%rough_moist
   rough_scale = rough_mom
   albedo      = Ice%albedo
   albedo_vis_dir = Ice%albedo_vis_dir
   albedo_nir_dir = Ice%albedo_nir_dir
   albedo_vis_dif = Ice%albedo_vis_dif
   albedo_nir_dif = Ice%albedo_nir_dif
   land_frac = 0.0
 endwhere

!---- land quantities ----

 where (Land%mask(:,:,1))
   t_surf = Land%t_surf(:,:,1)
   t_ca   = Land%t_ca  (:,:,1)
   q_surf = Land%tr    (:,:,1,1)
   rough_mom   = Land%rough_mom  (:,:,1)
   rough_heat  = Land%rough_heat (:,:,1)
   rough_moist = Land%rough_heat (:,:,1)
   rough_scale = Land%rough_scale(:,:,1)
   albedo      = Land%albedo     (:,:,1)
   albedo_vis_dir = Land%albedo_vis_dir(:,:,1)
   albedo_nir_dir = Land%albedo_nir_dir(:,:,1)
   albedo_vis_dif = Land%albedo_vis_dif(:,:,1)
   albedo_nir_dif = Land%albedo_nir_dif(:,:,1)
   land_frac = 1.0
 endwhere

!--- compute explicit fluxes and tendencies at all available points ---

   avail = .true.
   call surface_flux_2d (Atm%t_bot, Atm%tr_bot(:,:,isphum), Atm%u_bot, Atm%v_bot,        &
                      Atm%p_bot, Atm%z_bot,                              &
                      p_surf, t_surf, t_ca, q_surf, u_surf, v_surf,      &
                      rough_mom, rough_heat, rough_moist, rough_scale,   &
                      Atm%gust,                                          &
                      flux_t, flux_q, flux_lw, flux_u, flux_v,           &
                      cd_m,   cd_t, cd_q, wind,                          &
                      u_star, b_star, q_star,                            &
                      dhdt_surf, dedt_surf, dedq_surf, drdt_surf,        &
                      dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,     &
                      dt, Land%mask(:,:,1), seawater, avail              )

 ! additional calculation to avoid passing extra
 ! argument out of surface_flux (data duplication)
   drag_q = wind*cd_q

 ! put relevant quantities onto atmospheric boundary

   Boundary%t = t_surf
   Boundary%albedo = albedo
   Boundary%albedo_vis_dir = albedo_vis_dir
   Boundary%albedo_nir_dir = albedo_nir_dir
   Boundary%albedo_vis_dif = albedo_vis_dif
   Boundary%albedo_nir_dif = albedo_nir_dif
   Boundary%rough_mom   = rough_mom
   Boundary%land_frac = land_frac
   Boundary%u_flux    = flux_u
   Boundary%v_flux    = flux_v
   Boundary%dtaudu    = dtaudu_atm
   Boundary%dtaudv    = dtaudv_atm
   Boundary%u_star    = u_star
   Boundary%b_star    = b_star
   Boundary%q_star    = q_star

!=======================================================================
!-------------------- diagnostics section ------------------------------

 if (first_static) then
   if ( id_land_mask   > 0 ) used = send_data ( id_land_mask,   Boundary%land_frac, Time )
   ! near-surface heights
   if ( id_height2m  > 0) used = send_data ( id_height2m, z_ref_heat, Time )
   if ( id_height10m > 0) used = send_data ( id_height10m, z_ref_mom, Time )

   first_static = .false.
 endif
   if ( id_wind        > 0 ) used = send_data ( id_wind,        wind,         Time )
   if ( id_drag_moist  > 0 ) used = send_data ( id_drag_moist,  cd_q,         Time )
   if ( id_drag_heat   > 0 ) used = send_data ( id_drag_heat,   cd_t,         Time )
   if ( id_drag_mom    > 0 ) used = send_data ( id_drag_mom,    cd_m,         Time )
   if ( id_rough_moist > 0 ) used = send_data ( id_rough_moist, rough_moist,  Time )
   if ( id_rough_heat  > 0 ) used = send_data ( id_rough_heat,  rough_heat,   Time )
   if ( id_rough_mom   > 0 ) used = send_data ( id_rough_mom,   rough_mom,    Time )
   if ( id_u_star      > 0 ) used = send_data ( id_u_star,      u_star,       Time )
   if ( id_b_star      > 0 ) used = send_data ( id_b_star,      b_star,       Time )
   if ( id_q_star      > 0 ) used = send_data ( id_q_star,      q_star,       Time )
   if ( id_t_atm       > 0 ) used = send_data ( id_t_atm,       Atm%t_bot,    Time )
   if ( id_u_atm       > 0 ) used = send_data ( id_u_atm,       Atm%u_bot,    Time )
   if ( id_v_atm       > 0 ) used = send_data ( id_v_atm,       Atm%v_bot,    Time )
   if ( id_q_atm       > 0 ) used = send_data ( id_q_atm,       Atm%tr_bot(:,:,isphum),    Time )
   if ( id_p_atm       > 0 ) used = send_data ( id_p_atm,       Atm%p_bot,    Time )
   if ( id_z_atm       > 0 ) used = send_data ( id_z_atm,       Atm%z_bot,    Time )
   if ( id_gust        > 0 ) used = send_data ( id_gust,        Atm%gust,     Time )
   if ( id_u_flux      > 0 ) used = send_data ( id_u_flux,      flux_u,       Time )
   if ( id_v_flux      > 0 ) used = send_data ( id_v_flux,      flux_v,       Time )
   if ( id_albedo      > 0 ) used = send_data ( id_albedo,      albedo,       Time )
   if ( id_albedo_vis_dir > 0 ) used = send_data ( id_albedo_vis_dir, albedo_vis_dir, Time )
   if ( id_albedo_nir_dir > 0 ) used = send_data ( id_albedo_nir_dir, albedo_nir_dir, Time )
   if ( id_albedo_vis_dif > 0 ) used = send_data ( id_albedo_vis_dif, albedo_vis_dif, Time )
   if ( id_albedo_nir_dif > 0 ) used = send_data ( id_albedo_nir_dif, albedo_nir_dif, Time )

!---- ice fraction ----
   if ( id_ice_mask > 0 ) then
       where (Ice%ice_mask)
          ref = 1.0
       elsewhere
          ref = 0.0
       endwhere
       used = send_data ( id_ice_mask, ref, Time )
   endif

 ! diagnostics for fields at reference level

   if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
        id_u_ref > 0 .or. id_v_ref  > 0 .or. &
        id_q_ref > 0 ) then

        zrefm = z_ref_mom
        zrefh = z_ref_heat
       !---- optimize calculation ----
        if ( id_t_ref <= 0 ) zrefh = zrefm

         call mo_profile ( zrefm, zrefh, Atm%z_bot,           &
                          rough_mom, rough_heat, rough_moist, &
                          u_star, b_star, q_star,             &
                          del_m, del_h, del_q                 )

     !---- reference relative humidity ----
      if ( id_rh_ref > 0 .or.  id_q_ref > 0 .or. id_hurs > 0 .or. id_huss > 0) then
         ref   = q_surf + (Atm%tr_bot(:,:,isphum)-q_surf) * del_q
         if (id_q_ref > 0) used = send_data ( id_q_ref, ref, Time )
         if (id_huss  > 0) used = send_data (id_huss,ref,Time)

         t_ref = t_ca + (Atm%t_bot-t_ca) * del_h
         !call escomp (t_ref, qs_ref)
         call compute_qs (t_ref, p_surf, qs_ref, q = ref)
         call compute_qs (t_ref, p_surf, qs_ref_cmip,  &
            q = ref, es_over_liq_and_ice = .true.)
         qs_ref = d622*qs_ref/(p_surf-d378*qs_ref)

         ref    = 100.*ref/qs_ref
         ref2   = 100.*ref/qs_ref_cmip

         if (id_rh_ref > 0) used = send_data ( id_rh_ref, ref, Time )
         if (id_hurs   > 0) used = send_data ( id_hurs, ref2, Time )
      endif

     !---- reference temperature ----
      if ( id_t_ref > 0 ) then
         ref = t_ca + (Atm%t_bot-t_ca) * del_h
         used = send_data ( id_t_ref, ref, Time )
      endif

     !---- reference u comp ----
      if ( id_u_ref > 0 .or. id_uas > 0) then
         ref = u_surf + (Atm%u_bot-u_surf) * del_m
         used = send_data ( id_u_ref, ref, Time )
      endif
      if ( id_uas > 0 ) used = send_data ( id_uas, ref, Time )

     !---- reference v comp ----
      if ( id_v_ref > 0 .or. id_vas > 0) then
         ref = v_surf + (Atm%v_bot-v_surf) * del_m
         used = send_data ( id_v_ref, ref, Time )
      endif
      if ( id_vas > 0 ) used = send_data ( id_vas, ref, Time )

      !    ------- reference-level absolute wind -----------
      if ( id_sfcWind > 0 ) then
              ref = sqrt((u_surf + (Atm%u_bot-u_surf) * del_m)**2 &
              +(v_surf + (Atm%v_bot-v_surf) * del_m)**2)
         if ( id_sfcWind  > 0 ) used = send_data ( id_sfcWind, ref , Time )
      endif


     !---- interp factors ----
      if ( id_del_h > 0 )  used = send_data ( id_del_h, del_h, Time )
      if ( id_del_m > 0 )  used = send_data ( id_del_m, del_m, Time )
      if ( id_del_q > 0 )  used = send_data ( id_del_q, del_q, Time )

   endif

  ! topographic roughness scale
   if (id_rough_scale > 0) then
     ref = (log(Atm%z_bot/rough_mom+1)/log(Atm%z_bot/rough_scale+1))**2
     used = send_data(id_rough_scale, ref, Time)
  endif

! lgs line below is from atm_land_ice_flux_exchange.F90  what should diag_atm be?
!   if ( id_tas         > 0 ) used = send_data ( id_tas, diag_atm, Time )
   if ( id_tas > 0 ) used = send_data ( id_tas, t_ref, Time )
   if ( id_psl > 0 ) used = send_data ( id_psl, Atm%slp , Time )


!=======================================================================

 end subroutine sfc_boundary_layer

!#######################################################################

subroutine flux_down_from_atmos (Time, Atm, Land, Ice,  &
                         Atmos_boundary, Land_boundary, Ice_boundary )

type       (time_type), intent(in)  :: Time
type (atmos_data_type), intent(in)  :: Atm
type  (land_data_type), intent(in)  :: Land
type   (ice_data_type), intent(in)  :: Ice
type(land_ice_atmos_boundary_type),intent(in)   :: Atmos_boundary
type(atmos_land_boundary_type),    intent(inout):: Land_boundary
type(atmos_ice_boundary_type),     intent(inout):: Ice_boundary

!real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm

real, dimension(is:ie,js:je) :: gamma, dtmass, delta_t, delta_q, &
                                dflux_t, dflux_q, flux, deriv, dt_t_surf

 real, parameter :: CP_INV = 1./CP_AIR

  ! update stresses using atmos delta's

  flux_u = flux_u + Atm%Surf_Diff%delta_u * dtaudu_atm
  flux_v = flux_v + Atm%Surf_Diff%delta_v * dtaudv_atm

!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in flux_lw)

   flux_lw = Atm%flux_lw - flux_lw

!----- adjust fluxes for implicit dependence on atmosphere ----

   dtmass  = Atm%Surf_Diff%dtmass
   delta_t = Atm%Surf_Diff%delta_t
   delta_q = Atm%Surf_Diff%delta_tr(:,:,isphum)
   dflux_t = Atm%Surf_Diff%dflux_t
   dflux_q = Atm%Surf_Diff%dflux_tr(:,:,isphum)

 ! temperature

   gamma      =  1./ (1.0 - dtmass*(dflux_t + dhdt_atm*CP_INV))
   e_t_n      =  dtmass*dhdt_surf*CP_INV*gamma
   f_t_delt_n = (delta_t + dtmass * flux_t*CP_INV) * gamma

   flux_t     =  flux_t        + dhdt_atm * f_t_delt_n
   dhdt_surf  =  dhdt_surf     + dhdt_atm * e_t_n

! moisture

   gamma      =  1./ (1.0 - dtmass*(dflux_q + dedq_atm))
   e_q_n      =  dtmass*(dedt_surf+dedq_surf)*gamma
   f_q_delt_n = (delta_q  + dtmass * flux_q) * gamma

   flux_q     =  flux_q        + dedq_atm * f_q_delt_n
   dedt_surf  =  dedt_surf     + dedq_atm * e_q_n
   dedq_surf  =  dedq_surf     + dedq_atm * e_q_n

!-----------------------------------------------------------------------
!---- output fields on the land grid -------

      Land_boundary%lprec  (:,:,1) = 0.0
      Land_boundary%fprec  (:,:,1) = 0.0

   where (Land%mask(:,:,1))
      Land_boundary%t_flux (:,:,1) = flux_t
      Land_boundary%tr_flux(:,:,1,1) = flux_q
      Land_boundary%sw_flux(:,:,1) = Atm%flux_sw
      Land_boundary%sw_flux_down_vis_dir  (:,:,1) = Atm%flux_sw_down_vis_dir
      Land_boundary%sw_flux_down_total_dir(:,:,1) = Atm%flux_sw_down_total_dir
      Land_boundary%sw_flux_down_vis_dif  (:,:,1) = Atm%flux_sw_down_vis_dif
      Land_boundary%sw_flux_down_total_dif(:,:,1) = Atm%flux_sw_down_total_dif
      Land_boundary%lw_flux(:,:,1) = flux_lw
      Land_boundary%dhdt   (:,:,1) = dhdt_surf
!     Land_boundary%dedt   (:,:,1) = dedt_surf
      Land_boundary%dfdtr(:,:,1,1) = dedq_surf
      Land_boundary%drdt   (:,:,1) = drdt_surf
      Land_boundary%lprec  (:,:,1) = Atm%lprec
      Land_boundary%fprec  (:,:,1) = Atm%fprec
    ! Land_boundary%drag_q (:,:,1) = drag_q
      Land_boundary%p_surf (:,:,1) = p_surf

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

   elsewhere
      Ice_boundary%t_flux  = flux_t
      Ice_boundary%q_flux  = flux_q
      Ice_boundary%sw_flux = Atm%flux_sw
      Ice_boundary%lw_flux = flux_lw
      Ice_boundary%dhdt    = dhdt_surf
      Ice_boundary%dedt    = dedt_surf
      Ice_boundary%drdt    = drdt_surf
      Ice_boundary%lprec   = Atm%lprec
      Ice_boundary%fprec   = Atm%fprec
      Ice_boundary%u_star  = Atmos_boundary%u_star
      Ice_boundary%coszen  = Atm%coszen
   endwhere

   if (associated(Land_boundary%drag_q)) then
      where (Land%mask(:,:,1)) Land_boundary%drag_q(:,:,1) = drag_q
   endif
   if (associated(Land_boundary%lwdn_flux)) then
      where (Land%mask(:,:,1)) Land_boundary%lwdn_flux(:,:,1) = Atm%flux_lw
   endif
   if (associated(Land_boundary%cd_m)) then
      where (Land%mask(:,:,1)) Land_boundary%cd_m(:,:,1) = cd_m
   endif
   if (associated(Land_boundary%cd_t)) then
      where (Land%mask(:,:,1)) Land_boundary%cd_t(:,:,1) = cd_t
   endif
   if (associated(Land_boundary%bstar)) then
      where (Land%mask(:,:,1)) Land_boundary%bstar(:,:,1) = b_star
   endif
   if (associated(Land_boundary%ustar)) then
      where (Land%mask(:,:,1)) Land_boundary%ustar(:,:,1) = u_star
   endif
   if (associated(Land_boundary%wind)) then
      where (Land%mask(:,:,1)) Land_boundary%wind(:,:,1) = wind
   endif
   if (associated(Land_boundary%z_bot)) then
      where (Land%mask(:,:,1)) Land_boundary%z_bot(:,:,1) = Atm%z_bot
   endif

   deallocate ( flux_u, flux_v, dtaudu_atm, dtaudv_atm )

!-----------------------------------------------------------------------
!-------------------- diagnostics section ------------------------------

if ( id_u_flux > 0 ) used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
if ( id_tauu > 0 )   used = send_data ( id_tauu,  -Atmos_boundary%u_flux, Time )
if ( id_v_flux > 0 ) used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
if ( id_tauv > 0 )   used = send_data ( id_tauv,  -Atmos_boundary%v_flux, Time )

!-----------------------------------------------------------------------

end subroutine flux_down_from_atmos

!#######################################################################

subroutine flux_up_to_atmos (Time, Land, Ice, Boundary )

 type(time_type),      intent(in)  :: Time
 type(land_data_type), intent(inout)  :: Land
 type(ice_data_type),  intent(inout)  :: Ice
 type(land_ice_atmos_boundary_type), intent(inout) :: Boundary


 real, dimension(is:ie,js:je) :: t_surf_new, dt_t_surf, &
                                 q_surf_new, dt_q_surf, &
                                 delta_t_n,  delta_q_n, &
                                 t_ca_new,   dt_t_ca

 ! compute surface temperature change

  where (Land%mask(:,:,1))
     t_surf_new = Land%t_surf(:,:,1)
     t_ca_new   = Land%t_ca  (:,:,1)
  elsewhere
     t_surf_new = Ice%t_surf
     t_ca_new   = t_surf_new
  endwhere

 !??????? should this be done in land model ??????
  call escomp (t_surf_new, q_surf_new)
  where (Land%mask(:,:,1))
     q_surf_new = Land%tr(:,:,1,1)
  elsewhere
     q_surf_new = d622*q_surf_new/(p_surf-d378*q_surf_new)
  endwhere

  dt_t_ca   = t_ca_new   - t_ca   ! changes in near-surface T
  dt_t_surf = t_surf_new - t_surf ! changes in radiative T
  dt_q_surf = q_surf_new - q_surf ! changes in near-surface q

 ! adjust fluxes and atmospheric increments for
 ! implicit dependence on surface temperature

  flux_t        = flux_t      + dt_t_ca  *dhdt_surf
  flux_lw       = flux_lw     - dt_t_surf*drdt_surf
  Boundary%dt_t = f_t_delt_n  + dt_t_ca  *e_t_n

  where (Land%mask(:,:,1))
     flux_q                     = flux_q      + dt_q_surf*dedq_surf
     Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_q_surf*e_q_n
  elsewhere
     flux_q                     = flux_q      + dt_t_surf*dedt_surf
     Boundary%dt_tr(:,:,isphum) = f_q_delt_n  + dt_t_surf*e_q_n
  endwhere

!print *, 'PE,dt_t(L)(mn,mx)=',mpp_pe(),minval(Boundary%dt_t,mask=Land%mask(:,:,1)),maxval(Boundary%dt_t,mask=Land%mask(:,:,1))
!print *, 'PE,dt_q(L)(mn,mx)=',mpp_pe(),minval(Boundary%dt_q,mask=Land%mask(:,:,1)),maxval(Boundary%dt_q,mask=Land%mask(:,:,1))
!print *, 'PE,dt_t(I)(mn,mx)=',mpp_pe(),minval(Boundary%dt_t,mask=Ice%mask),maxval(Boundary%dt_t,mask=Ice%mask)
!print *, 'PE,dt_q(I)(mn,mx)=',mpp_pe(),minval(Boundary%dt_q,mask=Ice%mask),maxval(Boundary%dt_q,mask=Ice%mask)

!=======================================================================
!-------------------- diagnostics section ------------------------------

   if ( id_t_surf > 0 ) used = send_data ( id_t_surf, t_surf_new, Time )
   if ( id_t_ca   > 0 ) used = send_data ( id_t_ca,   t_ca_new,   Time )
   if ( id_q_surf > 0 ) used = send_data ( id_q_surf, q_surf_new, Time )
   if ( id_t_flux > 0 ) used = send_data ( id_t_flux, flux_t,     Time )
   if ( id_r_flux > 0 ) used = send_data ( id_r_flux, flux_lw,    Time )
   if ( id_q_flux > 0 ) used = send_data ( id_q_flux, flux_q,     Time )
   if ( id_evspsbl> 0 ) used = send_data ( id_evspsbl, flux_q,    Time )
   if ( id_hfls   > 0 ) used = send_data ( id_hfls,   HLV*flux_q, Time )
   if ( id_hfss   > 0 ) used = send_data ( id_hfss,   flux_t,     Time )
   if ( id_ts     > 0 ) used = send_data ( id_ts,     t_surf_new, Time )

   call sum_diag_integral_field ('evap', flux_q*86400.)

!=======================================================================
!---- deallocate module storage ----


   deallocate (f_t_delt_n, f_q_delt_n, e_t_n, e_q_n)
   deallocate (dhdt_surf, dedt_surf, dedq_surf, &
               drdt_surf, dhdt_atm,  dedq_atm,  &
               flux_t, flux_q, flux_lw, drag_q)
   deallocate (t_surf, p_surf, t_ca, q_surf )
   deallocate (cd_t, cd_m, b_star, u_star, wind)

!-----------------------------------------------------------------------

 end subroutine flux_up_to_atmos

!#######################################################################

 subroutine flux_exchange_init ( Time, Atm, Land, Ice,   &
                          !      atmos_land_boundary,    &
                                 atmos_ice_boundary,     &
                                 land_ice_atmos_boundary )

 type       (time_type), intent(in)  :: Time
 type (atmos_data_type), intent(in)  :: Atm
 type  (land_data_type), intent(in)  :: Land
 type   (ice_data_type), intent(in)  :: Ice
!type(atmos_land_boundary_type),    intent(inout) :: atmos_land_boundary
 type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary
 type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary


 integer :: j, kd
 integer :: isc, iec, jsc, jec
 real :: xx, lat

!-----------------------------------------------------------------------
!------ read namelist ------

   if (do_read_nml) call read_namelist

!-----------------------------------------------------------------------
!----- save the grid indices -----

   call mpp_get_compute_domain (Atm%Domain, is, ie, js, je )

!--------- write version number and namelist ------------------

   call write_version_number (version, tag)
   if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=flux_exchange_nml)


   call diag_integral_field_init ('evap', 'f6.3')
   call diag_field_init ( Time, Atm%axes(1:2) )

!----- find out number of atmospheric prognostic tracers and index of specific
!      humidity in the tracer table
   call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
                            num_prog=n_atm_tr)
   isphum = get_tracer_index( MODEL_ATMOS, 'sphum' )
   if (isphum==NO_TRACER) &
      call error_mesg('flux_exchange_mod', 'Cannot find water vapor in ATM tracer table', FATAL)
!-----------------------------------------------------------------------
!------ allocate atmos_land_boundary ------

   call mpp_get_compute_domain ( Land%Domain, isc, iec, jsc, jec )
   if (isc /= is .or. iec /= ie .or. jsc /= js .or. jec /= je ) &
         call error_mesg ('flux_exchange_init', 'land model '// &
                'domain does not match atmosphere domain', FATAL)
   kd = size(Land%mask,3) ! must be 1 (should check)
 ! allocate( atmos_land_boundary%t_flux (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%q_flux (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%lw_flux(is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%sw_flux(is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%sw_flux_down_vis_dir  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%sw_flux_down_total_dir(is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%sw_flux_down_vis_dif  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%sw_flux_down_total_dif(is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%lprec (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%fprec (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%dhdt  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%dedt  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%dedq  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%drdt  (is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%drag_q(is:ie,js:je,kd) )
 ! allocate( atmos_land_boundary%p_surf(is:ie,js:je,kd) )

!------ initialize boundary values ------

 ! atmos_land_boundary%t_flux=0.0
 ! atmos_land_boundary%q_flux=0.0
 ! atmos_land_boundary%lw_flux=0.0
 ! atmos_land_boundary%sw_flux=0.0
 ! atmos_land_boundary%sw_flux_down_vis_dir=0.0
 ! atmos_land_boundary%sw_flux_down_total_dir=0.0
 ! atmos_land_boundary%sw_flux_down_vis_dif=0.0
 ! atmos_land_boundary%sw_flux_down_total_dif=0.0
 ! atmos_land_boundary%lprec=0.0
 ! atmos_land_boundary%fprec=0.0
 ! atmos_land_boundary%dhdt=0.0
 ! atmos_land_boundary%dedt=0.0
 ! atmos_land_boundary%dedq=0.0
 ! atmos_land_boundary%drdt=0.0
 ! atmos_land_boundary%drag_q=0.0
 ! atmos_land_boundary%p_surf=0.0

!-----------------------------------------------------------------------
!------ allocate atmos ice boundary ------

   call mpp_get_compute_domain ( Ice%Domain, isc, iec, jsc, jec )
   if (isc /= is .or. iec /= ie .or. jsc /= js .or. jec /= je ) &
          call error_mesg ('flux_exchange_init', 'ice model '// &
                'domain does not match atmosphere domain', FATAL)

   allocate( atmos_ice_boundary%u_star(is:ie,js:je) )
   allocate( atmos_ice_boundary%t_flux(is:ie,js:je) )
   allocate( atmos_ice_boundary%q_flux(is:ie,js:je) )
   allocate( atmos_ice_boundary%lw_flux(is:ie,js:je) )
   allocate( atmos_ice_boundary%sw_flux(is:ie,js:je) )
   allocate( atmos_ice_boundary%lprec(is:ie,js:je) )
   allocate( atmos_ice_boundary%fprec(is:ie,js:je) )
   allocate( atmos_ice_boundary%dhdt(is:ie,js:je) )
   allocate( atmos_ice_boundary%dedt(is:ie,js:je) )
   allocate( atmos_ice_boundary%drdt(is:ie,js:je) )
   allocate( atmos_ice_boundary%coszen(is:ie,js:je) )

!------ initialize boundary values ------

   atmos_ice_boundary%u_star=0.0
   atmos_ice_boundary%t_flux=0.0
   atmos_ice_boundary%q_flux=0.0
   atmos_ice_boundary%lw_flux=0.0
   atmos_ice_boundary%sw_flux=0.0
   atmos_ice_boundary%lprec=0.0
   atmos_ice_boundary%fprec=0.0
   atmos_ice_boundary%dhdt=0.0
   atmos_ice_boundary%dedt=0.0
   atmos_ice_boundary%drdt=0.0
   atmos_ice_boundary%coszen=0.0

!-----------------------------------------------------------------------
!------ allocate land-ice-atmos boundary

   allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
   allocate( land_ice_atmos_boundary%t_ref(is:ie,js:je) )
   allocate( land_ice_atmos_boundary%q_ref(is:ie,js:je) )
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

!------ initialize boundary values ------

   land_ice_atmos_boundary%t=273.0
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

   call surface_flux_init()
!-----------------------------------------------------------------------

   do_init = .false.

!-----------------------------------------------------------------------

 end subroutine flux_exchange_init

!#######################################################################

 subroutine read_namelist

 integer :: ierr, io

   read (input_nml_file, nml=flux_exchange_nml, iostat=io)
   ierr = check_nml_error(io, 'flux_exchange_nml')

   do_read_nml = .false.

   end subroutine read_namelist

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
                        vrange = (/ -400., 400. /), &
                        frange = (/ -0.01, 1.01 /)
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

  id_t_ca       = &
       register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
       'canopy air temperature',    'deg_k', &
       range=trange    )

  id_q_atm     = &
       register_diag_field ( mod_name, 'q_atm',     atmos_axes, Time, &
       'specific humidity at btm level',    'kg/kg')

  id_q_surf     = &
       register_diag_field ( mod_name, 'q_surf',     atmos_axes, Time, &
       'surface specific humidity',    'kg/kg')

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
   id_q_ref = &
       register_diag_field ( mod_name, 'q_ref', atmos_axes, Time,     &
       'specific humidity at '//trim(label_zh), 'kg/kg', missing_value=-1.0)

   id_del_h      = &
   register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
                        'ref height interp factor for heat', 'none' )
   id_del_m      = &
   register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
                        'ref height interp factor for momentum','none' )
   id_del_q      = &
   register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
                        'ref height interp factor for moisture','none' )
   id_albedo      = &
   register_diag_field ( mod_name, 'albedo',      atmos_axes, Time,     &
                        'surface albedo','none' )
   id_albedo_vis_dir = &
   register_diag_field ( mod_name, 'albedo_vis_dir', atmos_axes, Time,     &
                        'VIS direct surface albedo','none' )
   id_albedo_nir_dir = &
   register_diag_field ( mod_name, 'albedo_nir_dir', atmos_axes, Time,     &
                        'NIR direct surface albedo','none' )
   id_albedo_vis_dif = &
   register_diag_field ( mod_name, 'albedo_vis_dif', atmos_axes, Time,     &
                        'VIS diffuse surface albedo','none' )
   id_albedo_nir_dif = &
   register_diag_field ( mod_name, 'albedo_nir_dif', atmos_axes, Time,     &
                        'NIR diffuse surface albedo','none' )

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
#ifndef use_AM3_physics
    id_tas = register_cmip_diag_field_2d ( mod_name, 'tas', Time, &
                            'Near-Surface Air Temperature', 'K' , &
                             standard_name='air_temperature' )
    if ( id_tas > 0 .and. id_height2m > 0) &
       call diag_field_add_attribute( id_tas, 'coordinates', 'height2m' )
! in diag table include height2m wherever tas is included

    id_evspsbl = register_cmip_diag_field_2d ( mod_name, 'evspsbl', Time, &
                                             'Evaporation', 'kg m-2 s-1', &
                                   standard_name='water_evaporation_flux' )

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

    id_hfls = register_cmip_diag_field_2d ( mod_name, 'hfls', Time, &
                        'Surface Upward Latent Heat Flux', 'W m-2', &
                    standard_name='surface_upward_latent_heat_flux' )
    if ( id_hfls > 0 ) &
       call diag_field_add_attribute( id_hfls, 'comment', 'Lv*evap' )

    id_hfss = register_cmip_diag_field_2d ( mod_name, 'hfss', Time, &
                      'Surface Upward Sensible Heat Flux', 'W m-2', &
                  standard_name='surface_upward_sensible_heat_flux' )
#endif

!-----------------------------------------------------------------------

end subroutine diag_field_init


!########################################################################

subroutine flux_exchange_end (Atm)

type (atmos_data_type), intent(in)  :: Atm
integer :: unit

end subroutine flux_exchange_end

!############################################################################
! copied from surface_flux_mod

subroutine surface_flux_2d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t,    flux_q,     flux_r,    flux_u,    flux_v,              &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt,        land,       seawater,  avail  )


  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:,:) :: land,  seawater, avail
  real, intent(in),  dimension(:,:) :: &
       t_atm,     q_atm_in,   u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist, rough_scale, gust
  real, intent(out), dimension(:,:) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, intent(inout), dimension(:,:) :: q_surf
  real, intent(in) :: dt

  ! ---- local vars -----------------------------------------------------------
  integer :: j

  do j = 1, size(t_atm,2)
     call surface_flux (                                           &
          t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),     v_atm(:,j),     p_atm(:,j),  z_atm(:,j), &
          p_surf(:,j),    t_surf(:,j),     t_ca(:,j),      q_surf(:,j),                             &
          u_surf(:,j),    v_surf(:,j),                                                              &
          rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),           &
          flux_t(:,j),    flux_q(:,j),     flux_r(:,j),    flux_u(:,j),    flux_v(:,j),             &
          cd_m(:,j),      cd_t(:,j),       cd_q(:,j),                                               &
          w_atm(:,j),     u_star(:,j),     b_star(:,j),     q_star(:,j),                            &
          dhdt_surf(:,j), dedt_surf(:,j),  dedq_surf(:,j),  drdt_surf(:,j),                         &
          dhdt_atm(:,j),  dedq_atm(:,j),   dtaudu_atm(:,j), dtaudv_atm(:,j),                        &
          dt,             land(:,j),       seawater(:,j),  avail(:,j)  )
  end do

end subroutine surface_flux_2d


!############################################################################

end module flux_exchange_mod

