
module surface_flux_mod

!-----------------------------------------------------------------------

use       utilities_mod, only: error_mesg, FATAL, open_file,  &
                               close_file, get_my_pe
use   monin_obukhov_mod, only: mo_drag, mo_profile
use  sat_vapor_pres_mod, only: escomp, descomp, tcheck
use       constants_mod, only: cp, hlv, stefan, rdgas, rvgas, grav

implicit none
private

public  surface_flux, surface_profile

!-----------------------------------------------------------------------

   character(len=128) :: version = '$Id: surface_flux.F90,v 1.5 2001/03/06 19:02:19 fms Exp $'
   character(len=128) :: tag = '$Name: eugene $'

   logical :: do_init = .true.

   real, parameter :: d622   = rdgas/rvgas
   real, parameter :: d378   = 1.-d622
   real, parameter :: d608   = d378/d622
   real, parameter :: hlars  = hlv/rvgas
   real, parameter :: gcp    = grav/cp
   real, parameter :: kappa  = rdgas/cp

contains

!#######################################################################

subroutine surface_flux (                                              &
                 t_atm,     q_atm,      u_atm,     v_atm,              &
                 p_atm,     z_atm,                                     &
                 p_surf,    t_surf,     u_surf,    v_surf,             &
                 rough_mom, rough_heat, rough_moist, gust,  stomatal,  &
                 snow_depth, water_depth,  max_water,                  &
                 flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
                 cd_m,      cd_t,       cd_q,      w_atm,              &
                 u_star,    b_star,     q_star,    q_surf,             &
                 dhdt_surf, dedt_surf,  drdt_surf,                     &
                 dhdt_atm,  dedq_atm,   dtaudv_atm,                    &
                 dt,        land,       glacier,      avail            )

logical, intent(in), dimension(:) :: land, glacier, avail
                 
real, intent(in),  dimension(:) ::                                     &
                 t_atm,     q_atm,      u_atm,     v_atm,              &
                 p_atm,     z_atm,                                     &
                 p_surf,    t_surf,     u_surf,    v_surf,             &
                 rough_mom, rough_heat, rough_moist,  gust,  stomatal, &
                 snow_depth, water_depth, max_water

real, intent(out), dimension(:) ::                                     &
                 flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
                 dhdt_surf, dedt_surf,  drdt_surf,                     &
                 dhdt_atm,  dedq_atm,   dtaudv_atm  

real, intent(inout), dimension(:) :: cd_m, cd_t, cd_q
real, intent(out), dimension(:) :: w_atm, u_star, b_star, q_star, q_surf
real, intent(in) :: dt

!-----------------------------------------------------------------------

real, dimension(size(t_atm)) ::                                      &
                 th_atm,   tv_atm,   th_surf,   th_surf1,            &
                 e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
                 t_surf0,  t_surf1,  t_surf2,   u_dif,     v_dif,    &
                 rho,      rho_drag, drag_t,    drag_m,    drag_q,   &
                 beta

logical, dimension(size(t_atm)) :: bone_dry, snow_on_ground

real, parameter :: del_temp = 0.1    
           ! temperature increment for computation of flux derivatives

real :: del_temp_inv = 1.0/del_temp


integer :: i, nbad
!-----------------------------------------------------------------------

if (do_init) call surface_flux_init

!---- use local value of surf temp ----

   where (avail)
      t_surf0 = t_surf
   elsewhere
      t_surf0 = 200.   !  avoids out-of-bounds in es lookup 
   endwhere

!---- check surface temp range for es lookup ----

   call tcheck ( t_surf0, nbad )

!  ---- terminate when bad values ----
   if ( nbad > 0 ) then
        print *, 'pe, nbad, total = ', get_my_pe(), nbad, count(avail)
        print *, 'min,max ts=',minval(t_surf0),maxval(t_surf0)
        print *, 'ts=0, count=',count(t_surf0 < .0001)
        call error_mesg ( 'surface_flux',  &
        'surface temperatures outside the es lookup table range', FATAL)
   endif

!---- saturation vapor pressure lookup ----

t_surf1 = t_surf0 + del_temp

call escomp (t_surf0, e_sat)   ! saturation vapor pressure
call escomp (t_surf1, e_sat1)  ! perturbed  vapor pressure


!---------- compute only where available ----------

          beta = 1.0

where (land .and. .not.glacier) snow_on_ground = (snow_depth > 0.0) 

where (land .and. .not.glacier .and. .not. snow_on_ground) &
          beta = min(water_depth/(0.75*max_water), 1.0)

where (avail)

   q_sat   = d622*e_sat /(p_surf-d378*e_sat )  
   q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)     
            ! surface specific humidity at saturation

   tv_atm  = t_atm  * (1.0 + d608*q_atm)
            !  virtual temperature

   p_ratio = (p_surf/p_atm)**kappa
   th_atm  = tv_atm * p_ratio
            !  virtual potential temperature, using p_surf as reference 


   th_surf  = t_surf0 * (1.0 + d608*beta*q_sat )
   th_surf1 = t_surf1 * (1.0 + d608*beta*q_sat1)
            !  surface virtual (potential) temperature

    u_dif = u_surf - u_atm
    v_dif = v_surf - v_atm
            !  velocity components relative to surface

    w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + gust*gust)
            !  wind speed, adding on the  "gustiness"

endwhere

!-----------------------------------------------------------------------

!  monin-obukhov similarity theory 

call mo_drag (dt, th_atm, th_surf, z_atm,                &
              rough_mom, rough_heat, rough_moist, w_atm, &
              cd_m, cd_t, cd_q, u_star, b_star, avail    )

!call mo_drag (dt, th_atm, th_surf, z_atm, rough_mom, rough_heat, w_atm, &
!              cd_m, cd_t, u_star, b_star, avail)

!-----------------------------------------------------------------------

where (avail)

!---- surface layer drag coefficients ----

     drag_t = cd_t * w_atm
     drag_q = cd_q * w_atm
     drag_m = cd_m * w_atm

!---- density ----

     rho = p_atm / (rdgas * tv_atm)  

!---- sensible heat flux ----

     rho_drag = cp * drag_t * rho
     flux_t = rho_drag * (th_surf - th_atm)   ! flux of sensible heat (W/m**2)

     dhdt_surf =  rho_drag*(th_surf1 - th_surf)*del_temp_inv
                        ! d(sensible heat flux)/d(surface temperature)

     dhdt_atm  = -rho_drag*p_ratio   
                        ! d(sensible heat flux)/d(atmos temperature)

!---- evaporation ----

     drag_q = 1.0/((1.0/drag_q) + stomatal)

     rho_drag  =  beta * drag_q * rho 
     flux_q    =  rho_drag * (q_sat - q_atm)  
                            ! flux of latent heat  (W/m**2)

     dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
                            ! d(latent heat flux)/d(surface temperature)
     dedq_atm  = -rho_drag   
                            ! d(latent heat flux)/d(atmospheric mixing ratio)

     q_star = flux_q / (u_star * rho)             ! moisture scale
     q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity

!---- upward long wave radiation -----

     t_surf2   = t_surf0* t_surf0 
     flux_r    = stefan * t_surf2 * t_surf2           ! (W/m**2)
     drdt_surf = 4. * stefan * t_surf2 * t_surf0
                            ! d(upward longwave)/d(surface temperature)

!---- stresses ----

     rho_drag   = drag_m * rho
     flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
     flux_v     = rho_drag * v_dif   ! meridional component of stress 
     dtaudv_atm = -rho_drag          ! d(stress component)/d(atmos wind)

elsewhere

!------- zero-out un-available data in output only fields -------
     flux_t     = 0.0
     flux_q     = 0.0
     flux_r     = 0.0
     flux_u     = 0.0
     flux_v     = 0.0
     dhdt_surf  = 0.0
     dedt_surf  = 0.0
     drdt_surf  = 0.0
     dhdt_atm   = 0.0
     dedq_atm   = 0.0
     dtaudv_atm = 0.0
     u_star     = 0.0
     b_star     = 0.0
     q_star     = 0.0
     q_surf     = 0.0
     w_atm      = 0.0
endwhere

bone_dry = .false.
where(land .and. .not.glacier) bone_dry = flux_q > water_depth/dt

where(bone_dry) 
  flux_q = water_depth/dt
  dedt_surf = 0.0
  dedq_atm = 0.0
endwhere


end subroutine surface_flux

!#######################################################################

subroutine surface_profile (zref_mom, zref_heat, z_atm,             &
                            rough_mom, rough_heat, rough_moist,     &
                            u_star, b_star, q_star,                 &
                            del_mom, del_heat, del_moist, avail     )

real,    intent(in)                :: zref_mom, zref_heat
real,    intent(in),  dimension(:) :: z_atm, rough_mom, rough_heat, rough_moist,  &
                                      u_star, b_star, q_star
real,    intent(out), dimension(:) :: del_mom, del_heat, del_moist
logical, intent(in),  dimension(:) :: avail

  real, dimension(size(del_mom)) :: del_mom_local


!---- get scaling at zref separately for momentum and heat ----

        call mo_profile ( zref_mom, z_atm, rough_mom, rough_heat, rough_moist,  &
                          u_star, b_star, q_star, del_mom, del_heat, del_moist, avail )
                 
   if ( zref_mom /= zref_heat ) then
        call mo_profile ( zref_heat, z_atm, rough_mom, rough_heat, rough_moist,  &
                          u_star, b_star, q_star, del_mom_local, del_heat, del_moist,  &
                          avail )
   endif

!-----------------------------------------------------------------------

end subroutine surface_profile

!#######################################################################

subroutine surface_flux_init

integer :: unit

!----- write version number -----

   unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) &
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
   call close_file (unit)
   do_init = .false.

end subroutine surface_flux_init

!#######################################################################

end module surface_flux_mod

