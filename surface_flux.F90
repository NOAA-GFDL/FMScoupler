
module surface_flux_mod

!-----------------------------------------------------------------------

use       utilities_mod, only: error_mesg, FATAL, open_file,  &
                               close_file, get_my_pe, file_exist, &
			       check_nml_error
			       
use   monin_obukhov_mod, only: mo_drag, mo_profile

use  sat_vapor_pres_mod, only: escomp, descomp

use       constants_mod, only: cp, hlv, stefan, rdgas, rvgas, grav

implicit none
private

public  surface_flux

interface surface_flux
    module procedure  surface_flux_0d, surface_flux_1d, surface_flux_2d
end interface

!-----------------------------------------------------------------------

   character(len=128) :: version = '$Id: surface_flux.F90,v 1.6 2002/02/22 19:07:23 fms Exp $'
   character(len=128) :: tagname = '$Name: galway $'
   
   logical :: do_init = .true.

   real, parameter :: d622   = rdgas/rvgas
   real, parameter :: d378   = 1.-d622
   real, parameter :: hlars  = hlv/rvgas
   real, parameter :: gcp    = grav/cp
   real, parameter :: kappa  = rdgas/cp
   real            :: d608   = d378/d622
      ! d608 set to zero at initialization if the use of 
      ! virtual temperatures is turned off in namelist
      
      
! namelist with default values

   logical :: no_neg_q         = .false.  ! for backwards compatibility
   logical :: use_virtual_temp = .true. 
   logical :: alt_gustiness    = .false.
   real    :: gust_const       =  1.0

namelist /surface_flux_nml/ no_neg_q,         &
                            use_virtual_temp, &
                            alt_gustiness,    &
			    gust_const
   

contains

!#######################################################################

subroutine surface_flux_1d (                                           &
                 t_atm,     q_atm_in,   u_atm,     v_atm,              &
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
                 t_atm,     q_atm_in,   u_atm,     v_atm,              &
                 p_atm,     z_atm,                                     &
                 p_surf,    t_surf,     u_surf,    v_surf,             &
                 rough_mom, rough_heat, rough_moist,  gust,  stomatal, &
                 snow_depth, water_depth, max_water

real, intent(out), dimension(:) ::                                     &
                 flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
                 dhdt_surf, dedt_surf,  drdt_surf,                     &
                 dhdt_atm,  dedq_atm,   dtaudv_atm,                    &
		 w_atm,     u_star,     b_star,    q_star,  q_surf,    &
                 cd_m,      cd_t,       cd_q
real, intent(in) :: dt

!-----------------------------------------------------------------------

real, dimension(size(t_atm)) ::                                      &
                 th_atm,   tv_atm,   th_surf,   th_surf1,            &
                 e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
                 t_surf0,  t_surf1,  t_surf2,   u_dif,     v_dif,    &
                 rho,      rho_drag, drag_t,    drag_m,    drag_q,   &
                 beta,     q_atm

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

!!---- check surface temp range for es lookup ----
!Temperature check is done in table lookup, "tcheck" not needed or supported -arl
!
!   call tcheck ( t_surf0, nbad )
!
!!  ---- terminate when bad values ----
!   if ( nbad > 0 ) then
!        print *, 'pe, nbad, total = ', get_my_pe(), nbad, count(avail)
!        print *, 'min,max ts=',minval(t_surf0),maxval(t_surf0)
!        print *, 'ts=0, count=',count(t_surf0 < .0001)
!        call error_mesg ( 'surface_flux',  &
!        'surface temperatures outside the es lookup table range', FATAL)
!   endif
!
!---- saturation vapor pressure lookup ----

t_surf1 = t_surf0 + del_temp

call escomp (t_surf0, e_sat)   ! saturation vapor pressure
call escomp (t_surf1, e_sat1)  ! perturbed  vapor pressure

!----- check for negative atmospheric humidities

  where(avail) q_atm = q_atm_in
  if(no_neg_q) then
    where(avail .and. q_atm_in < 0.0) q_atm = 0.0
  endif

!----- evaporation factor

beta = 1.0
where (avail .and. land .and. .not.glacier) snow_on_ground = (snow_depth > 0.0) 
where (avail .and. land .and. .not.glacier .and. .not. snow_on_ground) &
          beta = min(water_depth/(0.75*max_water), 1.0)

!------ generate information deeded by monin_obukhov

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
endwhere

if(alt_gustiness) then
  where(avail) w_atm = max(sqrt(u_dif*u_dif + v_dif*v_dif), gust_const)
else
  where(avail) w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + gust*gust)
endif

!-----------------------------------------------------------------------

!  monin-obukhov similarity theory 

call mo_drag (th_atm, th_surf, z_atm,                    &
              rough_mom, rough_heat, rough_moist, w_atm, &
              cd_m, cd_t, cd_q, u_star, b_star, avail    )

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
                            ! flux of latent heat  (Kg/(m**2 s))

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
where(avail .and. land .and. .not.glacier) bone_dry = flux_q > water_depth/dt

where(bone_dry) 
  flux_q = water_depth/dt
  dedt_surf = 0.0
  dedq_atm = 0.0
endwhere

end subroutine surface_flux_1d

!#######################################################################

subroutine surface_flux_0d (                                               &
       t_atm_0,      q_atm_0,       u_atm_0,       v_atm_0,                &
       p_atm_0,      z_atm_0,                                              &
       p_surf_0,     t_surf_0,      u_surf_0,      v_surf_0,               &
       rough_mom_0,  rough_heat_0,  rough_moist_0, gust_0,    stomatal_0,  &
       snow_depth_0, water_depth_0, max_water_0,                           &
       flux_t_0,     flux_q_0,      flux_r_0,      flux_u_0,  flux_v_0,    &
       cd_m_0,       cd_t_0,        cd_q_0,        w_atm_0,                &
       u_star_0,     b_star_0,      q_star_0,      q_surf_0,               &
       dhdt_surf_0,  dedt_surf_0,   drdt_surf_0,                           &
       dhdt_atm_0,   dedq_atm_0,    dtaudv_atm_0,                          &
       dt,           land_0,        glacier_0                              )

logical, intent(in) :: land_0, glacier_0
                 
real, intent(in) ::                                                      &
      t_atm_0,      q_atm_0,       u_atm_0,        v_atm_0,              &
      p_atm_0,      z_atm_0,                                             &
      p_surf_0,     t_surf_0,      u_surf_0,       v_surf_0,             &
      rough_mom_0,  rough_heat_0,  rough_moist_0,  gust_0,  stomatal_0,  &
      snow_depth_0, water_depth_0, max_water_0

real, intent(out) ::                                                   &
      flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,     &
      dhdt_surf_0, dedt_surf_0,  drdt_surf_0,                          &
      dhdt_atm_0,  dedq_atm_0,   dtaudv_atm_0,                         &  
      w_atm_0,     u_star_0,     b_star_0,    q_star_0,  q_surf_0,     &
      cd_m_0,      cd_t_0,       cd_q_0
      
real, intent(in)    :: dt

logical, dimension(1) :: land, glacier, avail

real, dimension(1) ::                                             &
          t_atm,      q_atm,       u_atm,        v_atm,           &
          p_atm,      z_atm,                                      &
          p_surf,     t_surf,      u_surf,       v_surf,          &
          rough_mom,  rough_heat,  rough_moist,  gust,  stomatal, &
          snow_depth, water_depth, max_water,                     &
	  cd_m,       cd_t,        cd_q,                          &
          flux_t,    flux_q,     flux_r,    flux_u,  flux_v,      &
          dhdt_surf, dedt_surf,  drdt_surf,                       &
          dhdt_atm,  dedq_atm,   dtaudv_atm,                      &
	  w_atm,     u_star,     b_star,    q_star,  q_surf

avail = .true.

t_atm(1)       = t_atm_0
q_atm(1)       = q_atm_0
u_atm(1)       = u_atm_0
v_atm(1)       = v_atm_0
p_atm(1)       = p_atm_0
z_atm(1)       = z_atm_0
p_surf(1)      = p_surf_0
t_surf(1)      = t_surf_0
u_surf(1)      = u_surf_0
v_surf(1)      = v_surf_0
rough_mom(1)   = rough_mom_0
rough_heat(1)  = rough_heat_0
rough_moist(1) = rough_moist_0
gust(1)        = gust_0
stomatal(1)    = stomatal_0
snow_depth(1)  = snow_depth_0
water_depth(1) = water_depth_0
max_water(1)   = max_water_0

call surface_flux_1d (                                                 &
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
                 dt,        land,       glacier,   avail               )
		 
flux_t_0     = flux_t(1)
flux_q_0     = flux_q(1)
flux_r_0     = flux_r(1)
flux_u_0     = flux_u(1)
flux_v_0     = flux_v(1)
dhdt_surf_0  = dhdt_surf(1)
dedt_surf_0  = dedt_surf(1)
drdt_surf_0  = drdt_surf(1)
dhdt_atm_0   = dhdt_atm(1)
dedq_atm_0   = dedq_atm(1)
dtaudv_atm_0 = dtaudv_atm(1)
w_atm_0      = w_atm(1)
u_star_0     = u_star(1)
b_star_0     = b_star(1)
q_star_0     = q_star(1)
q_surf_0     = q_surf(1)
cd_m_0       = cd_m(1)
cd_t_0       = cd_t(1)
cd_q_0       = cd_q(1)


end subroutine surface_flux_0d

!#######################################################################


subroutine surface_flux_2d (                                      &
       t_atm,      q_atm,       u_atm,       v_atm,               &
       p_atm,      z_atm,                                         &
       p_surf,     t_surf,      u_surf,      v_surf,              &
       rough_mom,  rough_heat,  rough_moist, gust,    stomatal,   &
       snow_depth, water_depth, max_water,                        &
       flux_t,     flux_q,      flux_r,      flux_u,  flux_v,     &
       cd_m,       cd_t,        cd_q,        w_atm,               &
       u_star,     b_star,      q_star,      q_surf,              &
       dhdt_surf,  dedt_surf,   drdt_surf,                        &
       dhdt_atm,   dedq_atm,    dtaudv_atm,                       &
       dt,         land,        glacier                           )

logical, intent(in), dimension(:,:) :: land, glacier
                 
real, intent(in), dimension(:,:) ::                              &
      t_atm,      q_atm,       u_atm,        v_atm,              &
      p_atm,      z_atm,                                         &
      p_surf,     t_surf,      u_surf,       v_surf,             &
      rough_mom,  rough_heat,  rough_moist,  gust,  stomatal,    &
      snow_depth, water_depth, max_water

real, intent(out), dimension(:,:) ::                             &
      flux_t,    flux_q,     flux_r,    flux_u,  flux_v,         &
      dhdt_surf, dedt_surf,  drdt_surf,                          &
      dhdt_atm,  dedq_atm,   dtaudv_atm,                         &  
      w_atm,     u_star,     b_star,    q_star,  q_surf,         &
      cd_m,      cd_t,       cd_q
      
real, intent(in)    :: dt

logical, dimension(size(t_atm,1)) :: avail
integer :: j

avail = .true.


do j = 1, size(t_atm,2)
  call surface_flux_1d (                                                              &
     t_atm(:,j),      q_atm(:,j),       u_atm(:,j),       v_atm(:,j),                 &
     p_atm(:,j),      z_atm(:,j),                                                     &
     p_surf(:,j),     t_surf(:,j),      u_surf(:,j),      v_surf(:,j),                &
     rough_mom(:,j),  rough_heat(:,j),  rough_moist(:,j), gust(:,j),   stomatal(:,j), &
     snow_depth(:,j), water_depth(:,j), max_water(:,j),                               &
     flux_t(:,j),     flux_q(:,j),      flux_r(:,j),      flux_u(:,j), flux_v(:,j),   &
     cd_m(:,j),       cd_t(:,j),        cd_q(:,j),        w_atm(:,j),                 &
     u_star(:,j),     b_star(:,j),      q_star(:,j),      q_surf(:,j),                &
     dhdt_surf(:,j),  dedt_surf(:,j),   drdt_surf(:,j),                               &
     dhdt_atm(:,j),   dedq_atm(:,j),    dtaudv_atm(:,j),			      &
     dt,              land(:,j),        glacier(:,j),     avail                       )
end do
end subroutine surface_flux_2d

!#######################################################################

subroutine surface_flux_init

integer :: unit, ierr, io

!------ read namelist ------

   if ( file_exist('input.nml')) then
      unit = open_file ('input.nml', action='read')
      ierr=1; do while (ierr /= 0)
         read  (unit, nml=surface_flux_nml, iostat=io, end=10)
         ierr = check_nml_error(io,'surface_flux_nml')
      enddo
 10   call close_file (unit)
   endif

!----- write version number -----

   unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) &
        write (unit,'(/,80("="),/(a))') trim(version), trim(tagname)
   call close_file (unit)

   if(.not. use_virtual_temp) d608 = 0.0
   
   do_init = .false.

end subroutine surface_flux_init

!#######################################################################

end module surface_flux_mod

