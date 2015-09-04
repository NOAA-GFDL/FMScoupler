!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief Driver program for the calculation of fluxes on the exchange grids
!!
!! \section surface_flux_config Surface Flux Configuration
!!
!! surface_flux_mod is configured via the surface_flux_nml namelist in the `input.nml` file.
!! The following table are the available namelist variables.
!!
!! <table>
!!   <tr>
!!     <th>Variable Name</th>
!!     <th>Type</th>
!!     <th>Default Value</th>
!!     <th>Description</th>
!!   </tr>
!!   <tr>
!!     <td>no_neg_q</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If q_atm_in (specific humidity) is negative (because of numerical
!!       truncation), then override with 0.</td>
!!   </tr>
!!   <tr>
!!     <td>use_virtual_temp</td>
!!     <td>logical</td>
!!     <td>.TRUE.</td>
!!     <td>If true, use virtual potential temp to calculate the stability of the
!!       surface layer. if false, use potential temp.</td>
!!   </tr>
!!   <tr>
!!     <td>alt_gustiness</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>An alternative formulation for gustiness calculation. A minimum bound
!!       on the wind speed used influx calculations, with the bound equal to
!!       gust_const.</td>
!!   </tr>
!!   <tr>
!!     <td>old_dtaudv</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>The derivative of surface wind stress w.r.t. the zonal wind and
!!       meridional wind are approximated by the same tendency.</td>
!!   </tr>
!!   <tr>
!!     <td>use_mixing_ratio</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>An option to provide capability to run the Manabe Climate form of the
!!       surface flux (coded for legacy purposes). </td>
!!   </tr>
!!   <tr>
!!     <td>gust_const</td>
!!     <td>real</td>
!!     <td>1.0</td>
!!     <td>Constant for alternative gustiness calculation.</td>
!!   </tr>
!!   <tr>
!!     <td>gust_min</td>
!!     <td>real</td>
!!     <td>0.0</td>
!!     <td>Minimum gustiness used when alt_gustiness = false.</td>
!!   </tr>
!!   <tr>
!!     <td>ncar_ocean_flux</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Use NCAR climate model turbulent flux calculation described by Large
!!       and Yeager, NCAR Technical Document, 2004.</td>
!!   </tr>
!!   <tr>
!!     <td>ncar_ocean_flux_orig</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Use NCAR climate model turbulent flux calculation described by Large
!!       and Yeager, NCAR Technical Document, 2004, using the original GFDL
!!       implementation, which contains a bug in the specification of the exchange
!!       coefficient for the sensible heat.  This option is available for legacy
!!       purposes, and is not recommended for new experiments.</td>
!!   </tr>
!!   <tr>
!!     <td>raoult_sat_vap</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Reduce saturation vapor pressures to account for seawater
!!     salinity.</td>
!!   </tr>
!!   <tr>
!!     <td>do_simple</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!! </table
module surface_flux_mod

use             fms_mod, only: close_file, mpp_pe, mpp_root_pe, write_version_number
use             fms_mod, only: file_exist, check_nml_error, open_namelist_file, stdlog
use   monin_obukhov_mod, only: mo_drag, mo_profile, monin_obukhov_init
use  sat_vapor_pres_mod, only: escomp, descomp
use       constants_mod, only: cp_air, hlv, stefan, rdgas, rvgas, grav, vonkarm
use             mpp_mod, only: input_nml_file, FATAL, mpp_error

implicit none
private

! ==== public interface ======================================================
public  surface_flux, surface_flux_init
! ==== end of public interface ===============================================

!> \brief For the calculation of fluxes on the exchange grids.
!!
!! For the calculation of fluxes on the exchange grids.
interface surface_flux
!    module procedure surface_flux_0d
    module procedure surface_flux_1d
    module procedure surface_flux_2d
end interface


!-----------------------------------------------------------------------

character(len=*), parameter :: version = '$Id$'
character(len=*), parameter :: tagname = '$Name$'

logical :: module_is_initialized = .false.

real, parameter :: d622   = rdgas/rvgas
real, parameter :: d378   = 1.-d622
real, parameter :: hlars  = hlv/rvgas
real, parameter :: gcp    = grav/cp_air
real, parameter :: kappa  = rdgas/cp_air
real            :: d608   = d378/d622
      ! d608 set to zero at initialization if the use of
      ! virtual temperatures is turned off in namelist


! ---- namelist with default values ------------------------------------------
logical :: no_neg_q              = .false. !< If a_atm_in (specific humidity) is negative (because of numerical truncation),
                                           !! then override with 0.0
logical :: use_virtual_temp      = .true.  !< If .TRUE., use virtual potential temp to calculate the stability of the surface
                                           !! layer.  If .FALSE., use potential temp.
logical :: alt_gustiness         = .false. !< An alternaive formulation for gustiness calculation.  A minimum bound on the wind
                                           !! speed used influx calculations, with the bound equal to gust_const
logical :: old_dtaudv            = .false. !< The derivative of surface wind stress with respect to the zonal wind and meridional
                                           !! wind are approximated by the same tendency
logical :: use_mixing_ratio      = .false. !< An option to provide capability to run the Manabe Climate form of the surface flux
                                           !! (coded for legacy purposes).
real    :: gust_const            =  1.0 !< Constant for alternative gustiness calculation
real    :: gust_min              =  0.0 !< Minimum gustiness used when alt_gustiness is .FALSE.
logical :: ncar_ocean_flux       = .false. !< Use NCAR climate model turbulent flux calculation described by Large and Yeager,
                                           !! NCAR Technical Document, 2004
logical :: ncar_ocean_flux_orig  = .false. !< Use NCAR climate model turbulent flux calculation described by Large and Yeager,
                                           !! NCAR Technical Document, 2004, using the original GFDL implementation, which
                                           !! contains a bug in the specification of the exchange coefficient for the sensible
                                           !! heat.  This option is available for legacy purposes, and is not recommended for
                                           !! new experiments.
logical :: raoult_sat_vap        = .false. !< Reduce saturation vapor pressure to account for seawater
logical :: do_simple             = .false.


namelist /surface_flux_nml/ no_neg_q,             &
                            use_virtual_temp,     &
                            alt_gustiness,        &
                            gust_const,           &
                            gust_min,             &
                            old_dtaudv,           &
                            use_mixing_ratio,     &
                            ncar_ocean_flux,      &
                            ncar_ocean_flux_orig, &
                            raoult_sat_vap,       &
                            do_simple



contains


! ============================================================================
subroutine surface_flux_1d (                                           &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t, flux_q, flux_r, flux_u, flux_v,                           &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt,        land,      seawater,     avail  )
  ! ---- arguments -----------------------------------------------------------
  logical, intent(in), dimension(:) :: land, & !< Indicates where land exists (.TRUE. if exchange cell is on land
                                       seawater, & !< Indicates where liquid ocean water exists (.TRUE. if exchange cell is on liquid ocean water)
                                       avail !< .TRUE. where the exchange cell is active
  real, intent(in),  dimension(:) :: t_atm, & !< Air temp lowest atmospheric level.
                                     q_atm_in, & !< Mixing ratio at lowest atmospheric level (kg/kg).
                                     u_atm, & !< Zonal wind velocity at lowest atmospheric level.
                                     v_atm, & !< Meridional wind velocity at lowest atmospheric level.
                                     p_atm, & !< Pressure lowest atmospheric level.
                                     z_atm, & !< Height lowest atmospheric level.
                                     t_ca, & !< Air temp at the canopy
                                     p_surf, & !< Pressure at the Earth's surface
                                     t_surf, & !< Temp at the Earth's surface
                                     u_surf, & !< Zonal wind velocity at the Earth's surface
                                     v_surf, & !< Meridional wind velocity at the Earth's surface
                                     rough_mom, & !< Momentum roughness length
                                     rough_heat, & !< Heat roughness length
                                     rough_moist, & !< Moisture roughness length
                                     rough_scale, & !< Scale factor used to topographic roughness calculation
                                     gust !< Gustiness factor
  real, intent(out), dimension(:) :: flux_t, & !< Sensible heat flux
                                     flux_q, & !< Evaporative water flux
                                     flux_r, & !< Radiative energy flux
                                     flux_u, & !< Zonal momentum flux
                                     flux_v, & !< Meridional momentum flux
                                     dhdt_surf, & !< Sensible heat flux temperature sensitivity
                                     dedt_surf, & !< Moisture flux temperature sensitivity
                                     dedq_surf, & !< Moisture flux humidity sensitivity
                                     drdt_surf, & !< Radiative energy flux temperature sensitivity
                                     dhdt_atm, & !< Derivative of sensible heat flux over temp at the lowest atmos level
                                     dedq_atm, & !< Derivative of water vapor flux over temp at the lowest atmos level
                                     dtaudu_atm, & !< Derivative of zonal wind stress with respect to the lowest level zonal wind speed of the atmos
                                     dtaudv_atm, & !< Derivative of meridional wind stress with respect to the lowest level meridional wind speed of the atmos
                                     w_atm, & !< Absolute wind at the lowest atmospheric level
                                     u_star, & !< Turbulent velocity scale
                                     b_star, & !< Turbulent buoyant scale
                                     q_star, & !< Turbulent moisture scale
                                     cd_m, & !< Momentum exchange coefficient
                                     cd_t, & ! Heat exchange coefficient
                                     cd_q !< Moisture exchange coefficient
  real, intent(inout), dimension(:) :: q_surf !< Mixing ratio at the Earth's surface (kg/kg)
  real, intent(in) :: dt !< Time step (it is not used presently)

  ! ---- local constants -----------------------------------------------------
  ! temperature increment and its reciprocal value for comp. of derivatives
  real, parameter:: del_temp=0.1, del_temp_inv=1.0/del_temp

  ! ---- local vars ----------------------------------------------------------
  real, dimension(size(t_atm(:))) ::                          &
       thv_atm,  th_atm,   tv_atm,    thv_surf,            &
       e_sat,    e_sat1,   q_sat,     q_sat1,    p_ratio,  &
       t_surf0,  t_surf1,  u_dif,     v_dif,               &
       rho_drag, drag_t,    drag_m,   drag_q,    rho,      &
       q_atm,    q_surf0,  dw_atmdu,  dw_atmdv,  w_gust

  integer :: i, nbad


  if (.not. module_is_initialized) &
     call mpp_error(FATAL, "surface_flux_1d: surface_flux_init is not called")

  !---- use local value of surf temp ----

  t_surf0 = 200.   !  avoids out-of-bounds in es lookup
  where (avail)
     where (land)
        t_surf0 = t_ca
     elsewhere
        t_surf0 = t_surf
     endwhere
  endwhere

  t_surf1 = t_surf0 + del_temp

  call escomp ( t_surf0, e_sat  )  ! saturation vapor pressure
  call escomp ( t_surf1, e_sat1 )  ! perturbed  vapor pressure

  if(use_mixing_ratio) then
    ! surface mixing ratio at saturation
    q_sat   = d622*e_sat /(p_surf-e_sat )
    q_sat1  = d622*e_sat1/(p_surf-e_sat1)
  elseif(do_simple) then                  !rif:(09/02/09)
    q_sat   = d622*e_sat / p_surf
    q_sat1  = d622*e_sat1/ p_surf
  else
    ! surface specific humidity at saturation
    q_sat   = d622*e_sat /(p_surf-d378*e_sat )
    q_sat1  = d622*e_sat1/(p_surf-d378*e_sat1)
  endif

  ! initilaize surface air humidity according to surface type
  where (land)
     q_surf0 = q_surf ! land calculates it
  elsewhere
     q_surf0 = q_sat  ! everything else assumes saturated sfc humidity
  endwhere

  if (raoult_sat_vap) where (seawater) q_surf0 = 0.98 * q_surf0

  ! check for negative atmospheric humidities
  where(avail) q_atm = q_atm_in
  if(no_neg_q) then
     where(avail .and. q_atm_in < 0.0) q_atm = 0.0
  endif

  ! generate information needed by monin_obukhov
  where (avail)
     p_ratio = (p_surf/p_atm)**kappa

     tv_atm  = t_atm  * (1.0 + d608*q_atm)     ! virtual temperature
     th_atm  = t_atm  * p_ratio                ! potential T, using p_surf as refernce
     thv_atm = tv_atm * p_ratio                ! virt. potential T, using p_surf as reference
     thv_surf= t_surf0 * (1.0 + d608*q_surf0 ) ! surface virtual (potential) T
!     thv_surf= t_surf0                        ! surface virtual (potential) T -- just for testing tun off the q_surf

     u_dif = u_surf - u_atm                    ! velocity components relative to surface
     v_dif = v_surf - v_atm
  endwhere

  if(alt_gustiness) then
     do i = 1, size(avail)
        if (.not.avail(i)) cycle
        w_atm(i) = max(sqrt(u_dif(i)**2 + v_dif(i)**2), gust_const)
        ! derivatives of surface wind w.r.t. atm. wind components
        if(w_atm(i) > gust_const) then
           dw_atmdu(i) = u_dif(i)/w_atm(i)
           dw_atmdv(i) = v_dif(i)/w_atm(i)
        else
           dw_atmdu(i) = 0.0
           dw_atmdv(i) = 0.0
        endif
     enddo
  else
     if (gust_min > 0.0) then
       where(avail)
         w_gust = max(gust,gust_min) ! minimum gustiness
       end where
     else
       where(avail)
         w_gust = gust
       end where
     endif

     where(avail)
        w_atm = sqrt(u_dif*u_dif + v_dif*v_dif + w_gust*w_gust)
        ! derivatives of surface wind w.r.t. atm. wind components
        dw_atmdu = u_dif/w_atm
        dw_atmdv = v_dif/w_atm
     endwhere
  endif

  !  monin-obukhov similarity theory
  call mo_drag (thv_atm, thv_surf, z_atm,                  &
       rough_mom, rough_heat, rough_moist, w_atm,          &
       cd_m, cd_t, cd_q, u_star, b_star, avail             )

  ! override with ocean fluxes from NCAR calculation
  if (ncar_ocean_flux .or. ncar_ocean_flux_orig) then
    call  ncar_ocean_fluxes (w_atm, th_atm, t_surf0, q_atm, q_surf0, z_atm, &
                             seawater, cd_m, cd_t, cd_q, u_star, b_star     )
  end if

  where (avail)
     ! scale momentum drag coefficient on orographic roughness
     cd_m = cd_m*(log(z_atm/rough_mom+1)/log(z_atm/rough_scale+1))**2
     ! surface layer drag coefficients
     drag_t = cd_t * w_atm
     drag_q = cd_q * w_atm
     drag_m = cd_m * w_atm

     ! density
     rho = p_atm / (rdgas * tv_atm)

     ! sensible heat flux
     rho_drag = cp_air * drag_t * rho
     flux_t = rho_drag * (t_surf0 - th_atm)  ! flux of sensible heat (W/m**2)
     dhdt_surf =  rho_drag                   ! d(sensible heat flux)/d(surface temperature)
     dhdt_atm  = -rho_drag*p_ratio           ! d(sensible heat flux)/d(atmos temperature)

     ! evaporation
     rho_drag  =  drag_q * rho
     flux_q    =  rho_drag * (q_surf0 - q_atm) ! flux of water vapor  (Kg/(m**2 s))

     where (land)
        dedq_surf = rho_drag
        dedt_surf = 0
     elsewhere
        dedq_surf = 0
        dedt_surf =  rho_drag * (q_sat1 - q_sat) *del_temp_inv
     endwhere

     dedq_atm  = -rho_drag   ! d(latent heat flux)/d(atmospheric mixing ratio)

     q_star = flux_q / (u_star * rho)             ! moisture scale
     ! ask Chris and Steve K if we still want to keep this for diagnostics
     q_surf = q_atm + flux_q / (rho*cd_q*w_atm)   ! surface specific humidity

     ! upward long wave radiation
     flux_r    =   stefan*t_surf**4               ! (W/m**2)
     drdt_surf = 4*stefan*t_surf**3               ! d(upward longwave)/d(surface temperature)

     ! stresses
     rho_drag   = drag_m * rho
     flux_u     = rho_drag * u_dif   ! zonal      component of stress (Nt/m**2)
     flux_v     = rho_drag * v_dif   ! meridional component of stress

  elsewhere
     ! zero-out un-available data in output only fields
     flux_t     = 0.0
     flux_q     = 0.0
     flux_r     = 0.0
     flux_u     = 0.0
     flux_v     = 0.0
     dhdt_surf  = 0.0
     dedt_surf  = 0.0
     dedq_surf  = 0.0
     drdt_surf  = 0.0
     dhdt_atm   = 0.0
     dedq_atm   = 0.0
     u_star     = 0.0
     b_star     = 0.0
     q_star     = 0.0
     q_surf     = 0.0
     w_atm      = 0.0
  endwhere

  ! calculate d(stress component)/d(atmos wind component)
  dtaudu_atm = 0.0
  dtaudv_atm = 0.0
  if (old_dtaudv) then
     where(avail)
        dtaudv_atm = -rho_drag
        dtaudu_atm = -rho_drag
     endwhere
  else
     where(avail)
        dtaudu_atm = -cd_m*rho*(dw_atmdu*u_dif + w_atm)
        dtaudv_atm = -cd_m*rho*(dw_atmdv*v_dif + w_atm)
     endwhere
  endif

end subroutine surface_flux_1d


!#######################################################################
subroutine surface_flux_0d (                                                 &
     t_atm_0,     q_atm_0,      u_atm_0,     v_atm_0,   p_atm_0, z_atm_0,    &
     p_surf_0,    t_surf_0,     t_ca_0,      q_surf_0,                       &
     u_surf_0,    v_surf_0,                                                  &
     rough_mom_0, rough_heat_0, rough_moist_0, rough_scale_0, gust_0,        &
     flux_t_0,    flux_q_0,     flux_r_0,    flux_u_0,  flux_v_0,            &
     cd_m_0,      cd_t_0,       cd_q_0,                                      &
     w_atm_0,     u_star_0,     b_star_0,     q_star_0,                      &
     dhdt_surf_0, dedt_surf_0,  dedq_surf_0,  drdt_surf_0,                   &
     dhdt_atm_0,  dedq_atm_0,   dtaudu_atm_0, dtaudv_atm_0,                  &
     dt,          land_0,       seawater_0,  avail_0  )

  ! ---- arguments -----------------------------------------------------------
  logical, intent(in) :: land_0, & !< Indicates where land exists (.TRUE. if exchange cell is on land
                         seawater_0, & !< Indicates where liquid ocean water exists (.TRUE. if exchange cell is on liquid ocean water)
                         avail_0 !< .TRUE. where the exchange cell is active
  real, intent(in) :: t_atm_0, & !< Air temp lowest atmospheric level.
                      q_atm_0, & !< Mixing ratio at lowest atmospheric level (kg/kg).
                      u_atm_0, & !< Zonal wind velocity at lowest atmospheric level.
                      v_atm_0, & !< Meridional wind velocity at lowest atmospheric level.
                      p_atm_0, & !< Pressure lowest atmospheric level.
                      z_atm_0, & !< Height lowest atmospheric level.
                      t_ca_0, & !< Air temp at the canopy
                      p_surf_0, & !< Pressure at the Earth's surface
                      t_surf_0, & !< Temp at the Earth's surface
                      u_surf_0, & !< Zonal wind velocity at the Earth's surface
                      v_surf_0, & !< Meridional wind velocity at the Earth's surface
                      rough_mom_0, & !< Momentum roughness length
                      rough_heat_0, & !< Heat roughness length
                      rough_moist_0, & !< Moisture roughness length
                      rough_scale_0, & !< Scale factor used to topographic roughness calculation
                      gust_0 !< Gustiness factor
  real, intent(out) :: flux_t_0, & !< Sensible heat flux
                       flux_q_0, & !< Evaporative water flux
                       flux_r_0, & !< Radiative energy flux
                       flux_u_0, & !< Zonal momentum flux
                       flux_v_0, & !< Meridional momentum flux
                       dhdt_surf_0, & !< Sensible heat flux temperature sensitivity
                       dedt_surf_0, & !< Moisture flux temperature sensitivity
                       dedq_surf_0, & !< Moisture flux humidity sensitivity
                       drdt_surf_0, & !< Radiative energy flux temperature sensitivity
                       dhdt_atm_0, & !< Derivative of sensible heat flux over temp at the lowest atmos level
                       dedq_atm_0, & !< Derivative of water vapor flux over temp at the lowest atmos level
                       dtaudu_atm_0, & !< Derivative of zonal wind stress with respect to the lowest level zonal wind speed of the atmos
                       dtaudv_atm_0, & !< Derivative of meridional wind stress with respect to the lowest level meridional wind speed of the atmos
                       w_atm_0, & !< Absolute wind at the lowest atmospheric level
                       u_star_0, & !< Turbulent velocity scale
                       b_star_0, & !< Turbulent buoyant scale
                       q_star_0, & !< Turbulent moisture scale
                       cd_m_0, & !< Momentum exchange coefficient
                       cd_t_0, & ! Heat exchange coefficient
                       cd_q_0 !< Moisture exchange coefficient
  real, intent(inout) :: q_surf_0 !< Mixing ratio at the Earth's surface (kg/kg)
  real, intent(in) :: dt !< Time step (it is not used presently)

  ! ---- local vars ----------------------------------------------------------
  logical, dimension(1) :: land,  seawater, avail
  real, dimension(1) :: &
       t_atm,     q_atm,      u_atm,     v_atm,              &
       p_atm,     z_atm,      t_ca,                          &
       p_surf,    t_surf,     u_surf,    v_surf,             &
       rough_mom, rough_heat, rough_moist,  rough_scale, gust
  real, dimension(1) :: &
       flux_t,    flux_q,     flux_r,    flux_u,  flux_v,    &
       dhdt_surf, dedt_surf,  dedq_surf, drdt_surf,          &
       dhdt_atm,  dedq_atm,   dtaudu_atm,dtaudv_atm,         &
       w_atm,     u_star,     b_star,    q_star,             &
       cd_m,      cd_t,       cd_q
  real, dimension(1) :: q_surf


  avail = .true.

  t_atm(1)       = t_atm_0
  q_atm(1)       = q_atm_0
  u_atm(1)       = u_atm_0
  v_atm(1)       = v_atm_0
  p_atm(1)       = p_atm_0
  z_atm(1)       = z_atm_0
  t_ca(1)        = t_ca_0
  p_surf(1)      = p_surf_0
  t_surf(1)      = t_surf_0
  u_surf(1)      = u_surf_0
  v_surf(1)      = v_surf_0
  rough_mom(1)   = rough_mom_0
  rough_heat(1)  = rough_heat_0
  rough_moist(1) = rough_moist_0
  rough_scale(1) = rough_scale_0
  gust(1)        = gust_0
  q_surf(1)      = q_surf_0
  land(1)        = land_0
  seawater(1)    = seawater_0
  avail(1)       = avail_0

  call surface_flux_1d (                                                 &
       t_atm,     q_atm,      u_atm,     v_atm,     p_atm,     z_atm,    &
       p_surf,    t_surf,     t_ca,      q_surf,                         &
       u_surf,    v_surf,                                                &
       rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
       flux_t, flux_q, flux_r, flux_u, flux_v,                           &
       cd_m,      cd_t,       cd_q,                                      &
       w_atm,     u_star,     b_star,     q_star,                        &
       dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
       dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
       dt,        land,      seawater, avail  )

  flux_t_0     = flux_t(1)
  flux_q_0     = flux_q(1)
  flux_r_0     = flux_r(1)
  flux_u_0     = flux_u(1)
  flux_v_0     = flux_v(1)
  dhdt_surf_0  = dhdt_surf(1)
  dedt_surf_0  = dedt_surf(1)
  dedq_surf_0  = dedq_surf(1)
  drdt_surf_0  = drdt_surf(1)
  dhdt_atm_0   = dhdt_atm(1)
  dedq_atm_0   = dedq_atm(1)
  dtaudu_atm_0 = dtaudu_atm(1)
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
  logical, intent(in), dimension(:,:) :: land, & !< Indicates where land exists (.TRUE. if exchange cell is on land
                                         seawater, & !< Indicates where liquid ocean water exists (.TRUE. if exchange cell is on liquid ocean water)
                                         avail !< .TRUE. where the exchange cell is active
  real, intent(in),  dimension(:,:) :: t_atm, & !< Air temp lowest atmospheric level.
                                       q_atm_in, & !< Mixing ratio at lowest atmospheric level (kg/kg).
                                       u_atm, & !< Zonal wind velocity at lowest atmospheric level.
                                       v_atm, & !< Meridional wind velocity at lowest atmospheric level.
                                       p_atm, & !< Pressure lowest atmospheric level.
                                       z_atm, & !< Height lowest atmospheric level.
                                       t_ca, & !< Air temp at the canopy
                                       p_surf, & !< Pressure at the Earth's surface
                                       t_surf, & !< Temp at the Earth's surface
                                       u_surf, & !< Zonal wind velocity at the Earth's surface
                                       v_surf, & !< Meridional wind velocity at the Earth's surface
                                       rough_mom, & !< Momentum roughness length
                                       rough_heat, & !< Heat roughness length
                                       rough_moist, & !< Moisture roughness length
                                       rough_scale, & !< Scale factor used to topographic roughness calculation
                                       gust !< Gustiness factor
  real, intent(out), dimension(:,:) :: flux_t, & !< Sensible heat flux
                                       flux_q, & !< Evaporative water flux
                                       flux_r, & !< Radiative energy flux
                                       flux_u, & !< Zonal momentum flux
                                       flux_v, & !< Meridional momentum flux
                                       dhdt_surf, & !< Sensible heat flux temperature sensitivity
                                       dedt_surf, & !< Moisture flux temperature sensitivity
                                       dedq_surf, & !< Moisture flux humidity sensitivity
                                       drdt_surf, & !< Radiative energy flux temperature sensitivity
                                       dhdt_atm, & !< Derivative of sensible heat flux over temp at the lowest atmos level
                                       dedq_atm, & !< Derivative of water vapor flux over temp at the lowest atmos level
                                       dtaudu_atm, & !< Derivative of zonal wind stress with respect to the lowest level zonal wind speed of the atmos
                                       dtaudv_atm, & !< Derivative of meridional wind stress with respect to the lowest level meridional wind speed of the atmos
                                       w_atm, & !< Absolute wind at the lowest atmospheric level
                                       u_star, & !< Turbulent velocity scale
                                       b_star, & !< Turbulent buoyant scale
                                       q_star, & !< Turbulent moisture scale
                                       cd_m, & !< Momentum exchange coefficient
                                       cd_t, & ! Heat exchange coefficient
                                       cd_q !< Moisture exchange coefficient
  real, intent(inout), dimension(:,:) :: q_surf !< Mixing ratio at the Earth's surface (kg/kg)
  real, intent(in) :: dt !< Time step (it is not used presently)

  ! ---- local vars -----------------------------------------------------------
  integer :: j

  do j = 1, size(t_atm,2)
     call surface_flux_1d (                                           &
          t_atm(:,j),     q_atm_in(:,j),   u_atm(:,j),     v_atm(:,j),     p_atm(:,j),     z_atm(:,j),    &
          p_surf(:,j),    t_surf(:,j),     t_ca(:,j),      q_surf(:,j),                                   &
          u_surf(:,j),    v_surf(:,j),                                                                    &
          rough_mom(:,j), rough_heat(:,j), rough_moist(:,j), rough_scale(:,j), gust(:,j),                 &
          flux_t(:,j),    flux_q(:,j),     flux_r(:,j),    flux_u(:,j),    flux_v(:,j),                   &
          cd_m(:,j),      cd_t(:,j),       cd_q(:,j),                                                     &
          w_atm(:,j),     u_star(:,j),     b_star(:,j),     q_star(:,j),                                  &
          dhdt_surf(:,j), dedt_surf(:,j),  dedq_surf(:,j),  drdt_surf(:,j),                               &
          dhdt_atm(:,j),  dedq_atm(:,j),   dtaudu_atm(:,j), dtaudv_atm(:,j),                              &
          dt,             land(:,j),       seawater(:,j),  avail(:,j)  )
  end do
end subroutine surface_flux_2d


! ============================================================================
!> \brief Initialization of the surface flux module--reads the nml.
subroutine surface_flux_init

! ---- local vars ----------------------------------------------------------
  integer :: unit, ierr, io

  ! read namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, surface_flux_nml, iostat=io)
      ierr = check_nml_error(io,'surface_flux_nml')
#else
  if ( file_exist('input.nml')) then
     unit = open_namelist_file ()
     ierr=1;
     do while (ierr /= 0)
        read  (unit, nml=surface_flux_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'surface_flux_nml')
     enddo
10   call close_file (unit)
  endif
#endif

  ! write version number
  call write_version_number(version, tagname)

  unit = stdlog()
  if ( mpp_pe() == mpp_root_pe() )  write (unit, nml=surface_flux_nml)

  if(.not. use_virtual_temp) d608 = 0.0

  call monin_obukhov_init()

  module_is_initialized = .true.

end subroutine surface_flux_init



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> \brief Over-ocean fluxes following Large and Yeager (used in NCAR models)           !
!!
!! Original  code: GFDL.Climate.Model.Info@noaa.gov <br \>
!! Update Jul2007: GFDL.Climate.Model.Info@noaa.gov (ch and ce exchange coeff bugfix)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
                              cd, ch, ce, ustar, bstar       )
real   , intent(in)   , dimension(:) :: u_del, t, ts, q, qs, z
logical, intent(in)   , dimension(:) :: avail
real   , intent(inout), dimension(:) :: cd, ch, ce, ustar, bstar

  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd_rt                                ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h            ! stability parameters
  real :: u, u10, tv, tstar, qstar, z0, xx, stab
  integer, parameter :: n_itts = 2
  integer               i, j

  if(ncar_ocean_flux_orig) then

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2;                !     10b (this code is wrong)
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)**2;                !     10c (this code is wrong)
             end do
         end if
      end do

  else

      do i=1,size(u_del(:))
         if (avail(i)) then
             tv = t(i)*(1+0.608*q(i));
             u = max(u_del(i), 0.5);                                 ! 0.5 m/s floor on wind (undocumented NCAR)
             u10 = u;                                                ! first guess 10m wind

             cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                ! L-Y eqn. 6a
             cd_n10_rt = sqrt(cd_n10);
             ce_n10 =                     34.6 *cd_n10_rt/1e3;       ! L-Y eqn. 6b
             stab = 0.5 + sign(0.5,t(i)-ts(i))
             ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;       ! L-Y eqn. 6c

             cd(i) = cd_n10;                                         ! first guess for exchange coeff's at z
             ch(i) = ch_n10;
             ce(i) = ce_n10;
             do j=1,n_itts                                           ! Monin-Obukhov iteration
                cd_rt = sqrt(cd(i));
                ustar(i) = cd_rt*u;                                   ! L-Y eqn. 7a
                tstar    = (ch(i)/cd_rt)*(t(i)-ts(i));                ! L-Y eqn. 7b
                qstar    = (ce(i)/cd_rt)*(q(i)-qs(i));                ! L-Y eqn. 7c
                bstar(i) = grav*(tstar/tv+qstar/(q(i)+1/0.608));
                zeta     = vonkarm*bstar(i)*z(i)/(ustar(i)*ustar(i)); ! L-Y eqn. 8a
                zeta     = sign( min(abs(zeta),10.0), zeta );         ! undocumented NCAR
                x2 = sqrt(abs(1-16*zeta));                            ! L-Y eqn. 8b
                x2 = max(x2, 1.0);                                    ! undocumented NCAR
                x = sqrt(x2);

                if (zeta > 0) then
                    psi_m = -5*zeta;                                    ! L-Y eqn. 8c
                    psi_h = -5*zeta;                                    ! L-Y eqn. 8c
                else
                    psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)); ! L-Y eqn. 8d
                    psi_h = 2*log((1+x2)/2);                                ! L-Y eqn. 8e
                end if

                u10 = u/(1+cd_n10_rt*(log(z(i)/10)-psi_m)/vonkarm);       ! L-Y eqn. 9
                cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
                cd_n10_rt = sqrt(cd_n10);
                ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
                stab = 0.5 + sign(0.5,zeta)
                ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
                z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic

                xx = (log(z(i)/10)-psi_m)/vonkarm;
                cd(i) = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
                xx = (log(z(i)/10)-psi_h)/vonkarm;
                ch(i) = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10b (corrected code)
                ce(i) = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd(i)/cd_n10) ! 10c (corrected code)
             end do
         end if
      end do

  endif

end subroutine ncar_ocean_fluxes


end module surface_flux_mod
