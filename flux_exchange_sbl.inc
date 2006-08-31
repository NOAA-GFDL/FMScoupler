! This include file contains the subroutine
! sfc_boundary_layer
!
! This subroutine makes extensive use of exchange grid global variables
! defined in the main flux_exchange module
!#######################################################################
! <SUBROUTINE NAME="sfc_boundary_layer">
!  <OVERVIEW>
!   Computes explicit fluxes as well as derivatives that will be used to compute an implicit flux correction. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!  The following quantities in the land_ice_atmos_boundary_type are computed:
!
!     
!         t_surf_atm = surface temperature (used for radiation)    (K)
!         albedo_atm = surface albedo      (used for radiation)    (nondimensional)
!      rough_mom_atm = surface roughness for momentum (m)
!      land_frac_atm = fractional area of land beneath an atmospheric
!                      grid box 
!         dtaudu_atm, dtaudv_atm = derivatives of wind stress w.r.t. the
!                      lowest level wind speed  (Pa/(m/s))
!         flux_u_atm = zonal wind stress  (Pa)
!         flux_v_atm = meridional wind stress (Pa)
!         u_star_atm = friction velocity (m/s)
!         b_star_atm = buoyancy scale    (m2/s)
!
!         (u_star and b_star are defined so that u_star**2 = magnitude
!           of surface stress divided by density of air at the surface, 
!           and u_star*b_star = buoyancy flux at the surface)
!
!   </PRE>
!  </DESCRIPTION>

!  <TEMPLATE>
!   call sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" dt" TYPE="real">
!   time step. 
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <INOUT NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </INOUT>
!  <INOUT NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </INOUT>
!  <INOUT NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </INOUT>
!  <INOUT NAME="Boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.   
!  </INOUT>
!
subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Land_Ice_Atmos_Boundary )

  real,                  intent(in)  :: dt
  type(time_type),       intent(in)  :: Time
  type(atmos_data_type), intent(inout)  :: Atm
  type(land_data_type),  intent(inout)  :: Land
  type(ice_data_type),   intent(inout)  :: Ice
  type(land_ice_atmos_boundary_type), intent(inout) :: Land_Ice_Atmos_Boundary

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
       ex_t_ref,      &
       ex_qs_ref,     &
       ex_del_m,      &
       ex_del_h,      &
       ex_del_q,      &
       ex_seawater

  real, dimension(n_xgrid_sfc,n_exch_tr) :: ex_tr_atm
  real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: diag_atm
  real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
  real, dimension(size(Ice%t_surf,1),size(Ice%t_surf,2),size(Ice%t_surf,3)) :: sea
  real    :: zrefm, zrefh
  logical :: used
  character(32) :: tr_name ! tracer name
  integer :: tr, n, m ! tracer indices
  integer :: i, ind_flux = 1

  ! [1] check that the module was initialized
!   <ERROR MSG="must call flux_exchange_init first " STATUS="FATAL">
!      flux_exchange_init has not been called before calling sfc_boundary_layer.
!   </ERROR>
  if (do_init) call error_mesg ('flux_exchange_mod',  &
       'must call flux_exchange_init first', FATAL)
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(sfcClock)
  ! [2] allocate storage for variables that are also used in flux_up_to_atmos
  allocate ( &
       ex_t_surf   (n_xgrid_sfc),  &
       ex_p_surf   (n_xgrid_sfc),  &
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
  ex_t_surf   = 200.
  ex_u_surf   =   0.
  ex_v_surf   =   0.

  !---- do not use if relax time /= 0 ----
  ex_cd_t = 0.0
  ex_cd_m = 0.0
  ex_cd_q = 0.0
!-----------------------------------------------------------------------
!Balaji: data_override stuff moved from coupler_main
  call data_override ('ATM', 't_bot',  Atm%t_bot , Time)
  call data_override ('ATM', 'z_bot',  Atm%z_bot , Time)
  call data_override ('ATM', 'p_bot',  Atm%p_bot , Time)
  call data_override ('ATM', 'u_bot',  Atm%u_bot , Time)
  call data_override ('ATM', 'v_bot',  Atm%v_bot , Time)
  call data_override ('ATM', 'p_surf', Atm%p_surf, Time)
  call data_override ('ATM', 'gust',   Atm%gust,   Time)
  do tr = 1,n_atm_tr
     call get_tracer_names( MODEL_ATMOS, tr, tr_name )
     call data_override('ATM', trim(tr_name)//'_bot', Atm%tr_bot(:,:,tr), Time)
  enddo
  do n = 1, atm%fields%num_bcs  !{
    do m = 1, atm%fields%bc(n)%num_fields  !{
      call data_override('ATM', atm%fields%bc(n)%field(m)%name,      &
           atm%fields%bc(n)%field(m)%values, Time, override = atm%fields%bc(n)%field(m)%override)
      ex_gas_fields_atm%bc(n)%field(m)%override = atm%fields%bc(n)%field(m)%override
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
  call data_override ('LND', 't_surf',     Land%t_surf,     Time)
  call data_override ('LND', 't_ca',       Land%t_ca,       Time)
  call data_override ('LND', 'rough_mom',  Land%rough_mom,  Time)
  call data_override ('LND', 'rough_heat', Land%rough_heat, Time)
  call data_override ('LND', 'albedo', Land%albedo,     Time)

! tracer data override
  do tr = 1, n_lnd_tr
     call get_tracer_names( MODEL_LAND, tr, tr_name )
     call data_override('LND', trim(tr_name)//'_surf', Land%tr(:,:,:,tr), Time)
  enddo
  do n = 1, ice%ocean_fields%num_bcs  !{
    do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
      call data_override('ICE', ice%ocean_fields%bc(n)%field(m)%name, ice%ocean_fields%bc(n)%field(m)%values, Time)
      if ( Ice%ocean_fields%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ice%ocean_fields%bc(n)%field(m)%id_diag, Ice%ocean_fields%bc(n)%field(m)%values, Time )
      endif  !}
    enddo  !} m
  enddo  !} n
  call data_override ('LND', 'albedo_vis_dir', Land%albedo_vis_dir,Time)
  call data_override ('LND', 'albedo_nir_dir', Land%albedo_nir_dir,Time)
  call data_override ('LND', 'albedo_vis_dif', Land%albedo_vis_dif,Time)
  call data_override ('LND', 'albedo_nir_dif', Land%albedo_nir_dif,Time)

!---- put atmosphere quantities onto exchange grid ----

  ! [4] put all the qantities we need onto exchange grid
  ! [4.1] put atmosphere quantities onto exchange grid
  call put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc, remap_method=remap_method)
  call put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc, remap_method=remap_method)
  do tr = 1,n_exch_tr
     call put_to_xgrid (Atm%tr_bot(:,:,tr_table(tr)%atm) , 'ATM', ex_tr_atm(:,tr), xmap_sfc, &
          remap_method=remap_method)
  enddo
  do n = 1, Atm%fields%num_bcs  !{
    do m = 1, Atm%fields%bc(n)%num_fields  !{
      call put_to_xgrid (Atm%fields%bc(n)%field(m)%values, 'ATM',            &
           ex_gas_fields_atm%bc(n)%field(m)%values, xmap_sfc, remap_method=remap_method)
    enddo  !} m
  enddo  !} n
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

  ex_land_frac = 0.0
  call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)

  ! [5] compute explicit fluxes and tendencies at all available points ---
  call some(xmap_sfc, ex_avail)
  call surface_flux (&
       ex_t_atm, ex_tr_atm(:,isphum),  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_tr_surf(:,isphum),                       &
       ex_u_surf, ex_v_surf,                                           &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_rough_scale,    &
       ex_gust,                                                        &
       ex_flux_t, ex_flux_tr(:,isphum), ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                    &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                     &
       ex_dhdt_surf, ex_dedt_surf, ex_dfdtr_surf(:,isphum),  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dfdtr_atm(:,isphum),  ex_dtaudu_atm, ex_dtaudv_atm,       &
       dt,                                                             &
       ex_land, ex_seawater .gt. 0,  ex_avail                          )

  zrefm = 10.0
  zrefh = z_ref_heat
  !      ---- optimize calculation ----
  call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
       ex_rough_heat, ex_rough_moist,          &
       ex_u_star, ex_b_star, ex_q_star,        &
       ex_del_m, ex_del_h, ex_del_q, ex_avail  )
  ex_ref_u = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
  ex_ref_v = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
  ex_u10 = sqrt(ex_ref_u**2 + ex_ref_v**2)
  do n = 1, ex_gas_fields_atm%num_bcs  !{
     if (atm%fields%bc(n)%use_10m_wind_speed) then  !{
        if (.not. ex_gas_fields_atm%bc(n)%field(ind_u10)%override) then  !{
           ex_gas_fields_atm%bc(n)%field(ind_u10)%values = ex_u10
        endif  !}
     endif  !}
  enddo  !} n
  ! fill derivatives for all tracers
  ! F = C0*u*rho*delta_q, C0*u*rho is the same for all tracers, copy from sphum
  do tr = 1,n_exch_tr
     if (tr==isphum) cycle
     ex_dfdtr_atm  (:,tr) = ex_dfdtr_atm  (:,isphum)
     ex_dfdtr_surf (:,tr) = ex_dfdtr_surf (:,isphum)
     ex_flux_tr    (:,tr) = ex_dfdtr_surf(:,tr)*(ex_tr_surf(:,tr)-ex_tr_atm(:,tr))
  enddo

! Combine explicit ocean flux and implicit land flux of extra flux fields.

  ! Calculate ocean explicit flux here

  call atmos_ocean_fluxes_calc(ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes, ex_seawater)

  ! The following statement is a concise version of what's following and worth
  ! looking into in the future.
  ! ex_flux_tr(:,itracer) = ex_gas_fluxes%bc(itracer_ocn)%field(ind_flux)%values(:)
  ! where(ex_seawater.gt.0) ex_flux_tr(:,itracer) = F_ocn
  do n = 1, ex_gas_fluxes%num_bcs  !{
    if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
      m = tr_table_map(ex_gas_fluxes%bc(n)%atm_tr_index)%exch
      do i = 1, size(ex_seawater(:))  !{
        if (ex_seawater(i) <= 0) then  !{
          if (ex_seawater(i) == 0) then  !{
            cycle                                ! over land, don't do anything
          else  !}{
            ex_flux_tr(i,m) = 0.0             ! pure ice exchange cell
          endif  !}
        else  !}{
    !      on ocean cell, flux is explicit therefore we zero derivitives. 
          ex_dfdtr_atm(i,m)  = 0.0
          ex_dfdtr_surf(i,m) = 0.0
          ex_flux_tr(i,m)    = ex_gas_fluxes%bc(n)%field(ind_flux)%values(i) * WTMAIR/1000.
        endif  !}
      enddo  !} i
    endif  !}
  enddo  !} n

  ! [5.2] override tracer fluxes and derivatives
  do tr = 1,n_exch_tr
     if( tr_table(tr)%atm == NO_TRACER ) cycle ! it should never happen, though

     call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
     ! [5.2.1] override tracer flux. Note that "sea" and "diag_land" are repeatedly used 
     ! as temporary storage for the values we are overriding fluxes and derivative with, 
     ! over ocean and land respectively
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

  ! NB: names of the override fields are constructed using tracer name and certain 
  ! prefixes / suffixes. For example, for the tracer named "sphum" (specific humidity) they will be:
  ! "ex_flux_sphum", "ex_dfdsphum_surf", and "ex_dfdsphum_atm".
  ! 
  ! For sensible heat flux names are "ex_flux_t", "ex_dhdt_surf", and "ex_dhdt_atm"; 
  ! despite the name those are actually in energy units, W/m2, W/(m2 degK), and
  ! W/(m2 degK) respectively

  where (ex_avail) ex_drag_q = ex_wind*ex_cd_q
  ! [6] get mean quantities on atmosphere grid
  ! [6.1] compute t surf for radiation
  ex_t_surf4 = ex_t_surf ** 4

  ! [6.2] put relevant quantities onto atmospheric boundary
  call get_from_xgrid (Land_Ice_Atmos_Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir,    'ATM',   &
                       ex_albedo_vis_dir   ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir,    'ATM',   &
                       ex_albedo_nir_dir   ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif,    'ATM',   &
                       ex_albedo_vis_dif   ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif,    'ATM',   &
                       ex_albedo_nir_dif   ,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc)

  call get_from_xgrid (Land_Ice_Atmos_Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudu,    'ATM', ex_dtaudu_atm, xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc)
  call get_from_xgrid (Land_Ice_Atmos_Boundary%q_star,    'ATM', ex_q_star    , xmap_sfc)

  Land_Ice_Atmos_Boundary%t = Land_Ice_Atmos_Boundary%t ** 0.25
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
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo, 'ATM',  ex_albedo_fix, xmap_sfc)
  ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)

  allocate ( ex_albedo_vis_dir_fix(n_xgrid_sfc) )
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dir, 'ATM',  &
           ex_albedo_vis_dir_fix, xmap_sfc)
  ex_albedo_vis_dir_fix = (1.0-ex_albedo_vis_dir) /  &
(1.0-ex_albedo_vis_dir_fix)
  allocate ( ex_albedo_nir_dir_fix(n_xgrid_sfc) )
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dir, 'ATM', &
 ex_albedo_nir_dir_fix, xmap_sfc)
  ex_albedo_nir_dir_fix = (1.0-ex_albedo_nir_dir) /  &
(1.0-ex_albedo_nir_dir_fix)
  allocate ( ex_albedo_vis_dif_fix(n_xgrid_sfc) )
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_vis_dif, 'ATM',   &
        ex_albedo_vis_dif_fix, xmap_sfc)
   ex_albedo_vis_dif_fix = (1.0-ex_albedo_vis_dif) /   &
       (1.0-ex_albedo_vis_dif_fix)
  allocate ( ex_albedo_nir_dif_fix(n_xgrid_sfc) )
  call put_to_xgrid (Land_Ice_Atmos_Boundary%albedo_nir_dif, 'ATM',  &
 ex_albedo_nir_dif_fix, xmap_sfc)
  ex_albedo_nir_dif_fix = (1.0-ex_albedo_nir_dif) /   &
       (1.0-ex_albedo_nir_dif_fix)

  !=======================================================================
  ! [7] diagnostics section

  !------- save static fields first time only ------
  if (first_static) then

     !------- land fraction ------
     if ( id_land_mask > 0 ) then
        used = send_data ( id_land_mask, Land_Ice_Atmos_Boundary%land_frac, Time )
     endif

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
  if ( id_rough_mom > 0 ) then
     used = send_data ( id_rough_mom, Land_Ice_Atmos_Boundary%rough_mom, Time )
  endif
  
  !------- friction velocity -----------
  if ( id_u_star > 0 ) then
     used = send_data ( id_u_star, Land_Ice_Atmos_Boundary%u_star, Time )
  endif
  
  !------- bouyancy -----------
  if ( id_b_star > 0 ) then
     used = send_data ( id_b_star, Land_Ice_Atmos_Boundary%b_star, Time )
  endif

  !------- moisture scale -----------
  if ( id_q_star > 0 ) then
     used = send_data ( id_q_star, Land_Ice_Atmos_Boundary%q_star, Time )
  endif

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
     if ( id_tr_atm(tr) > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_tr_atm(:,tr), xmap_sfc)
        used = send_data ( id_tr_atm(tr), diag_atm, Time )
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

  !-----------------------------------------------------------------------
  !--------- diagnostics for fields at reference level ---------
  
  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
       id_u_ref > 0 .or. id_v_ref  > 0 .or. &
       id_q_ref > 0 .or. id_q_ref_land > 0 .or. &
       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then
     
     zrefm = z_ref_mom
     zrefh = z_ref_heat
     !      ---- optimize calculation ----
     if ( id_t_ref <= 0 ) zrefh = zrefm
     
     call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
          ex_rough_heat, ex_rough_moist,          &
          ex_u_star, ex_b_star, ex_q_star,        &
          ex_del_m, ex_del_h, ex_del_q, ex_avail  )

     !    ------- reference relative humidity -----------
     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 .or. &
          id_q_ref > 0 .or. id_q_ref_land >0 ) then
        ex_ref   = ex_tr_surf(:,isphum) + (ex_tr_atm(:,isphum)-ex_tr_surf(:,isphum)) * ex_del_q
        if(id_q_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data(id_q_ref,diag_atm,Time)
        endif
        if(id_q_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data(id_q_ref_land, diag_land, &
                Land%tile_size, Time, mask=Land%mask)
        endif
        where (ex_avail) &
           ex_t_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        call escomp (ex_t_ref, ex_qs_ref)
        ex_qs_ref = d622*ex_qs_ref/(ex_p_surf-d378*ex_qs_ref)
        ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)

        if ( id_rh_ref_land > 0 ) then
           call get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_rh_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if(id_rh_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_rh_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference temp -----------
     if ( id_t_ref > 0 .or. id_t_ref_land > 0 ) then
        ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        if (id_t_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_tile_averaged_data ( id_t_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_t_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_t_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference u comp -----------
     if ( id_u_ref > 0 .or. id_u_ref_land > 0) then
        ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
        if ( id_u_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_u_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_u_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_u_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference v comp -----------
     if ( id_v_ref > 0 .or. id_v_ref_land > 0 ) then
        ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
        if ( id_v_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_tile_averaged_data ( id_v_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_v_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_v_ref, diag_atm, Time )
        endif
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

  endif
  ! topographic roughness scale
  if(id_rough_scale>0) then
     call get_from_xgrid (diag_atm, 'ATM',&
          (log(ex_z_atm/ex_rough_mom+1)/log(ex_z_atm/ex_rough_scale+1))**2, xmap_sfc)
     used = send_data(id_rough_scale, diag_atm, Time)
  endif

!Balaji
  call mpp_clock_end(sfcClock)
  call mpp_clock_end(cplClock)

!=======================================================================

end subroutine sfc_boundary_layer
! </SUBROUTINE>
