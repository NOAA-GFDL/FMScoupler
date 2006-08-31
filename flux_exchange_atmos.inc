! This include file contains 2 subroutines
! flux_down_from_atmos
! flux_up_to_atmos
!
! Both subroutines make extensive use of exchange grid global variables
! defined in the main flux_exchange module
!
!#######################################################################

! <SUBROUTINE NAME="flux_down_from_atmos">
!  <OVERVIEW>
!   Returns fluxes and derivatives corrected for the implicit treatment of atmospheric 
!   diffusive fluxes, as well as the increments in the temperature and specific humidity 
!   of the lowest atmospheric layer due to all explicit processes as well as the diffusive 
!   fluxes through the top of this layer. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements from Atmos_boundary are used as input: 
!
!        flux_u_atm = zonal wind stress (Pa)  
!        flux_v_atm = meridional wind stress (Pa)
!
!
!    The following elements of Land_boundary are output: 
!
!       flux_t_land = sensible heat flux (W/m2)
!       flux_q_land = specific humidity flux (Kg/(m2 s)
!      flux_lw_land = net longwave flux (W/m2), uncorrected for
!                     changes in surface temperature
!      flux_sw_land = net shortwave flux (W/m2)
!         dhdt_land = derivative of sensible heat flux w.r.t.
!                     surface temperature (on land model grid)  (W/(m2 K)
!         dedt_land = derivative of specific humidity flux w.r.t.
!                     surface temperature (on land model grid)  (Kg/(m2 s K)
!         drdt_land = derivative of upward longwave flux w.r.t.
!                     surface temperature (on land model grid) (W/(m2 K)
!        lprec_land = liquid precipitation, mass for one time step
!                      (Kg/m2)
!        fprec_land = frozen precipitation, mass for one time step
!                      (Kg/m2)
!
!
!    The following elements of Ice_boundary are output: 
!
!        flux_u_ice = zonal wind stress (Pa)
!        flux_v_ice = meridional wind stress (Pa)
!        coszen_ice = cosine of the zenith angle
!
!   </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_down_from_atmos (Time, Atm, Land, Ice, &
!		Atmos_boundary, Land_boundary, Ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </IN>
!  <INOUT NAME="Land_boundary" TYPE="atmos_land_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to land.
!  </INOUT>
!  <INOUT NAME="Ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!
subroutine flux_down_from_atmos (Time, Atm, Land, Ice, &
     Atmos_boundary, Land_boundary, Ice_boundary )

  type(time_type),       intent(in) :: Time
  type(atmos_data_type), intent(inout) :: Atm
  type(land_data_type),  intent(in) :: Land
  type(ice_data_type),   intent(in) :: Ice
  type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
  type(atmos_land_boundary_type),    intent(inout):: Land_boundary
  type(atmos_ice_boundary_type),     intent(inout):: Ice_boundary

  real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
       ex_flux_sw_dir,  &
                                    ex_flux_sw_dif,  &
      ex_flux_sw_down_vis_dir, ex_flux_sw_down_total_dir,  &
      ex_flux_sw_down_vis_dif, ex_flux_sw_down_total_dif,  &
       ex_flux_sw_vis, &
       ex_flux_sw_vis_dir, &
       ex_flux_sw_vis_dif, &
       ex_lprec, ex_fprec,      &
       ex_u_star_smooth,        &
       ex_coszen

  real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
       ex_delta_t, ex_delta_u, ex_delta_v, ex_dflux_t

  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_delta_tr, & ! tracer tendencies
       ex_dflux_tr    ! fracer flux change

  real    :: cp_inv
  logical :: used
  logical :: ov

  character(32) :: tr_name ! name of the tracer
  integer :: tr, n, m ! tracer indices

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
     call put_to_xgrid (Atm%flux_sw, 'ATM', ex_flux_sw, xmap_sfc)
     call put_to_xgrid (Atm%flux_sw_vis, 'ATM', ex_flux_sw_vis, xmap_sfc)
  end if
  call put_to_xgrid (Atm%flux_sw_dir, 'ATM', ex_flux_sw_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_vis_dir, 'ATM', ex_flux_sw_vis_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_dif, 'ATM', ex_flux_sw_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_vis_dif, 'ATM', ex_flux_sw_vis_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_vis_dir, 'ATM', ex_flux_sw_down_vis_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_total_dir, 'ATM', ex_flux_sw_down_total_dir, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_vis_dif, 'ATM', ex_flux_sw_down_vis_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_sw_down_total_dif, 'ATM', ex_flux_sw_down_total_dif, xmap_sfc)
  call put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc, remap_method=remap_method)
  !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
  !  if (mpp_pe()== mpp_root_pe()) print *,'LPREC', ccc

  call put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc)
  call put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc)

  call put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc)

  if(ex_u_star_smooth_bug) then
     call put_to_xgrid (Atmos_boundary%u_star, 'ATM', ex_u_star_smooth, xmap_sfc, remap_method=remap_method)
     ex_u_star = ex_u_star_smooth
  endif


! MOD changed the following two lines to put Atmos%surf_diff%delta_u and v
! on exchange grid instead of the stresses themselves so that only the 
! implicit corrections are filtered through the atmospheric grid not the
! stresses themselves
  call put_to_xgrid (Atm%Surf_Diff%delta_u, 'ATM', ex_delta_u, xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%Surf_Diff%delta_v, 'ATM', ex_delta_v, xmap_sfc, remap_method=remap_method)

  ! MOD update stresses using atmos delta's but derivatives on exchange grid
  ex_flux_u = ex_flux_u + ex_delta_u*ex_dtaudu_atm
  ex_flux_v = ex_flux_v + ex_delta_v*ex_dtaudv_atm

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----
!---- adjust 4 categories (vis/nir dir/dif) separately  ----
  if( sw1way_bug ) then ! to reproduce old results, may remove in the next major release.
!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----

     ex_flux_sw = ex_flux_sw * ex_albedo_fix


     ex_flux_sw_vis = ex_flux_sw_vis * ex_albedo_vis_dir_fix
     ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_vis_dir_fix
     ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_vis_dif_fix
     ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix
     ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix
  else 
     ex_flux_sw_dir = ex_flux_sw_dir - ex_flux_sw_vis_dir     ! temporarily nir/dir
     ex_flux_sw_dir = ex_flux_sw_dir * ex_albedo_nir_dir_fix  ! fix nir/dir
     ex_flux_sw_vis_dir = ex_flux_sw_vis_dir * ex_albedo_vis_dir_fix ! fix vis/dir
     ex_flux_sw_dir = ex_flux_sw_dir + ex_flux_sw_vis_dir     ! back to total dir

     ex_flux_sw_dif = ex_flux_sw_dif - ex_flux_sw_vis_dif     ! temporarily nir/dif
     ex_flux_sw_dif = ex_flux_sw_dif * ex_albedo_nir_dif_fix  ! fix nir/dif
     ex_flux_sw_vis_dif = ex_flux_sw_vis_dif * ex_albedo_vis_dif_fix ! fix vis/dif
     ex_flux_sw_dif = ex_flux_sw_dif + ex_flux_sw_vis_dif     ! back to total dif

     ex_flux_sw_vis = ex_flux_sw_vis_dir + ex_flux_sw_vis_dif ! legacy, remove later
     ex_flux_sw     = ex_flux_sw_dir     + ex_flux_sw_dif     ! legacy, remove later
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
!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in ex_flux_lw)

  ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----


  call put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc )
  call put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc )
  do tr = 1,n_exch_tr
     n = tr_table(tr)%atm
     call put_to_xgrid (Atm%Surf_Diff%delta_tr(:,:,n), 'ATM', ex_delta_tr(:,tr), xmap_sfc )
     call put_to_xgrid (Atm%Surf_Diff%dflux_tr(:,:,n), 'ATM', ex_dflux_tr(:,tr), xmap_sfc )
  enddo

  cp_inv = 1.0/cp_air

  where(ex_avail)

     ! temperature

     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
     ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
     ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    
     
     ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
     ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

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
     ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,isphum) + ex_dfdtr_atm(:,isphum)))
     ex_e_q_n      =  ex_dtmass * ex_dedt_surf * ex_gamma
     ex_dedt_surf  =  ex_dedt_surf + ex_dfdtr_atm(:,isphum) * ex_e_q_n
  endwhere
  do tr = 1,n_exch_tr
     where(ex_avail)
        ex_gamma   =  1.0 / (1.0 - ex_dtmass*(ex_dflux_tr(:,tr) + ex_dfdtr_atm(:,tr)))

        ex_e_tr_n(:,tr)      =  ex_dtmass*ex_dfdtr_surf(:,tr)*ex_gamma
        ex_f_tr_delt_n(:,tr) = (ex_delta_tr(:,tr)+ex_dtmass*ex_flux_tr(:,tr))*ex_gamma    
     
        ex_flux_tr(:,tr)     =  ex_flux_tr(:,tr) + ex_dfdtr_atm(:,tr)*ex_f_tr_delt_n(:,tr) 
        ex_dfdtr_surf(:,tr)  =  ex_dfdtr_surf(:,tr) + ex_dfdtr_atm(:,tr)*ex_e_tr_n(:,tr)
     endwhere
  enddo
!-----------------------------------------------------------------------
!---- output fields on the land grid -------

  call get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dir, 'LND', ex_flux_sw_down_vis_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dir, 'LND', ex_flux_sw_down_total_dir,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_vis_dif, 'LND', ex_flux_sw_down_vis_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux_down_total_dif, 'LND', ex_flux_sw_down_total_dif,   xmap_sfc)
  call get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

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
     call data_override('LND', 'bstar',  Land_boundary%bstar, Time )
  endif

  Land_boundary%tr_flux(:,:,:,:) = 0.0
  Land_boundary%dfdtr(:,:,:,:) = 0.0
  do tr = 1,n_exch_tr
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        call get_from_xgrid (Land_boundary%tr_flux(:,:,:,n), 'LND', ex_flux_tr(:,tr), xmap_sfc)
        call get_from_xgrid (Land_boundary%dfdtr(:,:,:,n),   'LND', ex_dfdtr_surf(:,tr), xmap_sfc)
     endif
  enddo

!  current time is Time: is that ok? not available in land_data_type
!Balaji: data_override calls moved here from coupler_main
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
  enddo

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

  call get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_tr(:,isphum), xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux,  'OCN', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis,  'OCN', ex_flux_sw_vis,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_dir,  'OCN', ex_flux_sw_dir,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis_dir,  'OCN', ex_flux_sw_vis_dir,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_dif,  'OCN', ex_flux_sw_dif,xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux_vis_dif,  'OCN', ex_flux_sw_vis_dif,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_star,   'OCN', ex_u_star,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%p,        'OCN', ex_p_surf,    xmap_sfc) ! mw mod

! Extra fluxes
  do n = 1, Ice_boundary%fluxes%num_bcs  !{
    do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
      call get_from_xgrid (Ice_boundary%fluxes%bc(n)%field(m)%values, 'OCN',  &
           ex_gas_fluxes%bc(n)%field(m)%values, xmap_sfc)
    enddo  !} m
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
  call data_override('ICE', 'sw_flux',Ice_boundary%sw_flux, Time)
  call data_override('ICE', 'sw_flux_vis',Ice_boundary%sw_flux_vis, Time)
  call data_override('ICE', 'sw_flux_vis_dir',Ice_boundary%sw_flux_vis_dir, Time)
  call data_override('ICE', 'sw_flux_dir',Ice_boundary%sw_flux_dir, Time, override=ov)
  call data_override('ICE', 'sw_flux_vis_dif',Ice_boundary%sw_flux_vis_dif, Time)
  call data_override('ICE', 'sw_flux_dif',Ice_boundary%sw_flux_dif, Time, override=ov)
  call data_override('ICE', 'sw_flux_dn',Ice_boundary%sw_flux, Time, override=ov)
!! ANY CHANGE NEEDED HERE FOR sw_flux_vis ??
! this converts net fluxes to downward fluxes
  if (ov) then
    Ice_boundary%sw_flux = Ice_boundary%sw_flux*(1-Ice%albedo)
! ?? NEEDED:    Ice_boundary%sw_flux_vis = Ice_boundary%sw_flux_vis*(1-Ice%albedo_vis_dir or dif, must choose ??)
! ?? NEEDED:    Ice_boundary%sw_flux_vis_dir = Ice_boundary%sw_flux_vis_dir*(1-Ice%albedo_vis_dir)
! ?? NEEDED:    Ice_boundary%sw_flux_vis_dif = Ice_boundary%sw_flux_vis_dif*(1-Ice%albedo_vis_dif)
! ?? NEEDED:    Ice_boundary%sw_flux_dir = Ice_boundary%sw_flux_dir*(1-Ice%albedo_( select vis or nir) dir)
! ?? NEEDED:    Ice_boundary%sw_flux_dif = Ice_boundary%sw_flux_dif*(1-Ice%albedo_ ( select vis or nir ) dif)
  endif
  call data_override('ICE', 'lprec',  Ice_boundary%lprec,   Time)
  call data_override('ICE', 'fprec',  Ice_boundary%fprec,   Time)
  call data_override('ICE', 'dhdt',   Ice_boundary%dhdt,    Time)
  call data_override('ICE', 'dedt',   Ice_boundary%dedt,    Time)
  call data_override('ICE', 'drdt',   Ice_boundary%drdt,    Time)
  call data_override('ICE', 'coszen', Ice_boundary%coszen,  Time)
  call data_override('ICE', 'p',      Ice_boundary%p,       Time)

  do n = 1, Ice_boundary%fluxes%num_bcs  !{
    do m = 1, Ice_boundary%fluxes%bc(n)%num_fields  !{
      call data_override('ICE', Ice_boundary%fluxes%bc(n)%field(m)%name,     &
           Ice_boundary%fluxes%bc(n)%field(m)%values, Time)
      if ( Ice_boundary%fluxes%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ice_boundary%fluxes%bc(n)%field(m)%id_diag, Ice_boundary%fluxes%bc(n)%field(m)%values, Time )
      endif  !}
    enddo  !} m
  enddo  !} n

  deallocate ( ex_flux_u, ex_flux_v, ex_dtaudu_atm, ex_dtaudv_atm)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- zonal wind stress -----------
  if ( id_u_flux > 0 ) then
     used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
  endif

  !------- meridional wind stress -----------
  if ( id_v_flux > 0 ) then
     used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
  endif

!Balaji
  call mpp_clock_end(fluxAtmDnClock)
  call mpp_clock_end(cplClock)
!=======================================================================

  end subroutine flux_down_from_atmos
! </SUBROUTINE>
!#######################################################################
! <SUBROUTINE NAME="flux_up_to_atmos">
!  <OVERVIEW>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Corrects the fluxes for consistency with the new surface temperatures in land 
!   and ice models. Final increments for temperature and specific humidity in the 
!   lowest atmospheric layer are computed and returned to the atmospheric model
!   so that it can finalize the increments in the rest of the atmosphere. 
!  <PRE>
!
!   The following elements of the land_ice_atmos_boundary_type are computed:
!        dt_t  = temperature change at the lowest
!                 atmospheric level (deg k)
!        dt_q  = specific humidity change at the lowest
!                 atmospheric level (kg/kg)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   Current time.
!  </IN>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Land_Ice_Atmos_Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice. 
!  </INOUT>
!
subroutine flux_up_to_atmos ( Time, Land, Ice, Land_Ice_Atmos_Boundary )

  type(time_type),      intent(in)  :: Time
  type(land_data_type), intent(inout)  :: Land
  type(ice_data_type),  intent(inout)  :: Ice
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary

  real, dimension(n_xgrid_sfc) ::  &
       ex_t_surf_new, &
       ex_dt_t_surf,  &
       ex_delta_t_n,  &
       ex_t_ca_new,   &
       ex_dt_t_ca
  real, dimension(n_xgrid_sfc,n_exch_tr) :: &
       ex_tr_surf_new,    & ! updated tracer values at the surface
       ex_dt_tr_surf,     & ! tendency of tracers at the surface
       ex_delta_tr_n

  real, dimension(size(Land_Ice_Atmos_Boundary%dt_t,1),size(Land_Ice_Atmos_Boundary%dt_t,2)) :: diag_atm, &
       evap_atm
  logical :: used

  integer :: tr       ! tracer index
  character(32) :: tr_name ! tracer name
  integer :: n, i, m


  !Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxAtmUpClock)
  !-----------------------------------------------------------------------
  !Balaji: data_override calls moved here from coupler_main
  call data_override ( 'ICE', 't_surf', Ice%t_surf,  Time)
  call data_override ( 'LND', 't_ca',   Land%t_ca,   Time)
  call data_override ( 'LND', 't_surf', Land%t_surf, Time)
  do tr = 1, n_lnd_tr
     call get_tracer_names( MODEL_LAND, tr, tr_name )
     call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
  enddo

  !----- compute surface temperature change -----

  ex_t_surf_new = 200.0

  call put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
  ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
  call put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
  call put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)

  !  call escomp(ex_t_ca_new, ex_q_surf_new)
  !  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new) 
  !  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

  where (ex_avail)
     ex_dt_t_ca   = ex_t_ca_new   - ex_t_ca   ! changes in near-surface T
     ex_dt_t_surf = ex_t_surf_new - ex_t_surf ! changes in radiative T
  endwhere

  !-----------------------------------------------------------------------
  !-----  adjust fluxes and atmospheric increments for 
  !-----  implicit dependence on surface temperature -----
  do tr = 1,n_exch_tr
     ! set up updated surface tracer field so that flux to atmos for absent 
     ! tracers is zero
     do i = 1, size(ex_avail(:))
        if(.not.ex_avail(i)) cycle
        if (ex_dfdtr_surf(i,tr)/=0) then
           ex_dt_tr_surf(i,tr) = -ex_flux_tr(i,tr)/ex_dfdtr_surf(i,tr)
        else
           ex_dt_tr_surf(i,tr) = 0
        endif
        ex_tr_surf_new(i,tr) = ex_tr_surf(i,tr)+ex_dt_tr_surf(i,tr)
     enddo
     ! get all tracers available from land, and calculate changes in near-tracer field
     n = tr_table(tr)%lnd
     if(n /= NO_TRACER ) then
        call put_to_xgrid ( Land%tr(:,:,:,n), 'LND', ex_tr_surf_new(:,tr), xmap_sfc )
     endif

     ! get all tracers available from ocean here 

     ! update tracer tendencies in the atmosphere
     where (ex_avail)
        ex_dt_tr_surf(:,tr) = ex_tr_surf_new(:,tr) - ex_tr_surf(:,tr)
        ex_delta_tr_n(:,tr) = ex_f_tr_delt_n(:,tr) + ex_dt_tr_surf(:,tr) * ex_e_tr_n(:,tr)
        ex_flux_tr(:,tr)    = ex_flux_tr(:,tr)     + ex_dt_tr_surf(:,tr) * ex_dfdtr_surf(:,tr)
     endwhere
  enddo

  ! re-calculate fluxes of specific humidity over ocean
  where (ex_avail.and..not.ex_land) 
     ! note that in this region (over ocean) ex_dt_t_surf == ex_dt_t_ca
     ex_delta_tr_n(:,isphum)  = ex_f_tr_delt_n(:,isphum) + ex_dt_t_surf * ex_e_q_n
     ex_flux_tr(:,isphum)     = ex_flux_tr(:,isphum)     + ex_dt_t_surf * ex_dedt_surf
  endwhere

  do tr=1,n_exch_tr
     ! get updated tracer tendency on the atmospheic grid
     n=tr_table(tr)%atm
     call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_tr(:,:,n), 'ATM', ex_delta_tr_n(:,tr), xmap_sfc)
  enddo

  ex_delta_t_n = 0.0

  where(ex_avail)
     ex_flux_t     = ex_flux_t  + ex_dt_t_ca   * ex_dhdt_surf
     ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
     ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_ca*ex_e_t_n
  endwhere

  !-----------------------------------------------------------------------
  !---- get mean quantites on atmospheric grid ----

  call get_from_xgrid (Land_Ice_Atmos_Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- new surface temperature -----------
  if ( id_t_surf > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
     used = send_data ( id_t_surf, diag_atm, Time )
  endif


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
     if ( id_tr_surf(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_tr_surf_new(:,tr), xmap_sfc)
        used = send_data ( id_tr_surf(tr), diag_atm, Time )
     endif
  enddo

  !------- sensible heat flux -----------
  if ( id_t_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
     used = send_data ( id_t_flux, diag_atm, Time )
  endif

  !------- net longwave flux -----------
  if ( id_r_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
     used = send_data ( id_r_flux, diag_atm, Time )
  endif

  !------- tracer fluxes ------------
  do tr=1,n_exch_tr
     if ( id_tr_flux(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_flux_tr(:,tr), xmap_sfc)
        used = send_data ( id_tr_flux(tr), diag_atm, Time )
        used = send_data ( id_tr_mol_flux(tr), diag_atm*1000./WTMAIR, Time)
     endif
  enddo

  !-----------------------------------------------------------------------
  !---- accumulate global integral of evaporation (mm/day) -----
  call get_from_xgrid (evap_atm, 'ATM', ex_flux_tr(:,isphum), xmap_sfc)
  if( id_q_flux > 0 ) used = send_data ( id_q_flux, evap_atm, Time)
  call sum_diag_integral_field ('evap', evap_atm*86400.)

  !=======================================================================
  !---- deallocate module storage ----
  deallocate ( &
       ex_t_surf   ,  &
       ex_p_surf   ,  &
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
       
  ex_land        )

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

!Balaji
  call mpp_clock_end(fluxAtmUpClock)
  call mpp_clock_end(cplClock)

!-----------------------------------------------------------------------

end subroutine flux_up_to_atmos
! </SUBROUTINE>
