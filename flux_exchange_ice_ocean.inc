! This include file contains 2 subroutines:
! flux_ice_to_ocean
! flux_ocean_to_ice
!
! Both subroutines make extensive use of exchange grid global variables
! defined in flux_exchange module
!#######################################################################
! <SUBROUTINE NAME="flux_ice_to_ocean">
!  <OVERVIEW>
!   Takes the ice model state (fluxes at the bottom of the ice) and interpolates it to the ocean model grid. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!   The following quantities are transferred from the Ice to the ice_ocean_boundary_type: 
!
!       flux_u = zonal wind stress (Pa)
!       flux_v = meridional wind stress (Pa)
!       flux_t = sensible heat flux (W/m2)
!       flux_q = specific humidity flux (Kg/m2/s)
!    flux_salt = salt flux (Kg/m2/s)
!      flux_sw = net (down-up) shortwave flux (W/m2)
!      flux_lw = net (down-up) longwave flux (W/m2)
!        lprec = mass of liquid precipitation since last
!                      time step (Kg/m2)
!        fprec = mass of frozen precipitation since last
!                time step (Kg/m2)
!       runoff = mass (?) of runoff since last time step
!                       (Kg/m2)
!       p_surf = surface pressure (Pa)
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ice_to_ocean ( Time, Ice, Ocean, Ice_Ocean_Boundary )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="Ice_Ocean_Boundary" TYPE="ice_ocean_boundary_type">
!   A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!
subroutine flux_ice_to_ocean ( Time, Ice, Ocean, Ice_Ocean_Boundary )
  type(time_type),        intent(in) :: Time
  type(ice_data_type),   intent(in)  :: Ice
  type(ocean_data_type), intent(in)  :: Ocean
!  real, dimension(:,:),   intent(out) :: flux_u_ocean,  flux_v_ocean,  &
!                                         flux_t_ocean,  flux_q_ocean,  &
!                                         flux_sw_ocean, flux_lw_ocean, &
!                                         lprec_ocean,   fprec_ocean,   &
!                                         runoff_ocean,  calving_ocean, &
!                                         flux_salt_ocean, p_surf_ocean
  type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary

  integer       :: m
  integer       :: n
  logical       :: used

!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxIceOceanClock)

  do n = 1, Ice%ocean_fluxes%num_bcs  !{
    do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
      if ( Ice%ocean_fluxes%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ice%ocean_fluxes%bc(n)%field(m)%id_diag, Ice%ocean_fluxes%bc(n)%field(m)%values, Time)
      endif  !}
    enddo  !} m
  enddo  !} n

  select case (Ice_Ocean_Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Ice_Ocean_Boundary%u_flux   ) )Ice_Ocean_Boundary%u_flux    = Ice%flux_u
     if( ASSOCIATED(Ice_Ocean_Boundary%v_flux   ) )Ice_Ocean_Boundary%v_flux    = Ice%flux_v
     if( ASSOCIATED(Ice_Ocean_Boundary%t_flux   ) )Ice_Ocean_Boundary%t_flux    = Ice%flux_t
     if( ASSOCIATED(Ice_Ocean_Boundary%q_flux   ) )Ice_Ocean_Boundary%q_flux    = Ice%flux_q
     if( ASSOCIATED(Ice_Ocean_Boundary%salt_flux) )Ice_Ocean_Boundary%salt_flux = Ice%flux_salt
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux  ) )Ice_Ocean_Boundary%sw_flux   = Ice%flux_sw
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis  ) )Ice_Ocean_Boundary%sw_flux_vis   = Ice%flux_sw_vis
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux_dir  ) )Ice_Ocean_Boundary%sw_flux_dir   = Ice%flux_sw_dir
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux_dif  ) )Ice_Ocean_Boundary%sw_flux_dif   = Ice%flux_sw_dif
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dir  ) )Ice_Ocean_Boundary%sw_flux_vis_dir   = Ice%flux_sw_vis_dir
     if( ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dif  ) )Ice_Ocean_Boundary%sw_flux_vis_dif   = Ice%flux_sw_vis_dif
     if( ASSOCIATED(Ice_Ocean_Boundary%lw_flux  ) )Ice_Ocean_Boundary%lw_flux   = Ice%flux_lw
     if( ASSOCIATED(Ice_Ocean_Boundary%lprec    ) )Ice_Ocean_Boundary%lprec     = Ice%lprec
     if( ASSOCIATED(Ice_Ocean_Boundary%fprec    ) )Ice_Ocean_Boundary%fprec     = Ice%fprec
     if( ASSOCIATED(Ice_Ocean_Boundary%runoff   ) )Ice_Ocean_Boundary%runoff    = Ice%runoff
     if( ASSOCIATED(Ice_Ocean_Boundary%calving  ) )Ice_Ocean_Boundary%calving   = Ice%calving
     if( ASSOCIATED(Ice_Ocean_Boundary%p        ) )Ice_Ocean_Boundary%p         = Ice%p_surf
! Extra fluxes
     do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
       do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
         if ( associated(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values) ) then  !{
           Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values = Ice%ocean_fluxes%bc(n)%field(m)%values
         endif  !}
       enddo  !} m
     enddo  !} n

  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if (ASSOCIATED(Ice_Ocean_Boundary%u_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_u, Ocean%Domain, Ice_Ocean_Boundary%u_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%v_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_v, Ocean%Domain, Ice_Ocean_Boundary%v_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%t_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_t, Ocean%Domain, Ice_Ocean_Boundary%t_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%q_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_q, Ocean%Domain, Ice_Ocean_Boundary%q_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%salt_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_salt, Ocean%Domain, Ice_Ocean_Boundary%salt_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw, Ocean%Domain, Ice_Ocean_Boundary%sw_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis, Ocean%Domain, Ice_Ocean_Boundary%sw_flux_vis)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux_dir)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_dir, Ocean%Domain, Ice_Ocean_Boundary%sw_flux_dir)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux_dif)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_dif, Ocean%Domain, Ice_Ocean_Boundary%sw_flux_dif)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dir)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis_dir, Ocean%Domain, Ice_Ocean_Boundary%sw_flux_vis_dir)
     if (ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dif)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw_vis_dif, Ocean%Domain, Ice_Ocean_Boundary%sw_flux_vis_dif)
     if (ASSOCIATED(Ice_Ocean_Boundary%lw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_lw, Ocean%Domain, Ice_Ocean_Boundary%lw_flux)
     if (ASSOCIATED(Ice_Ocean_Boundary%lprec)) &
          call mpp_redistribute(Ice%Domain, Ice%lprec, Ocean%Domain, Ice_Ocean_Boundary%lprec)
     if (ASSOCIATED(Ice_Ocean_Boundary%fprec)) &
          call mpp_redistribute(Ice%Domain, Ice%fprec, Ocean%Domain, Ice_Ocean_Boundary%fprec)
     if (ASSOCIATED(Ice_Ocean_Boundary%runoff)) &
          call mpp_redistribute(Ice%Domain, Ice%runoff, Ocean%Domain, Ice_Ocean_Boundary%runoff)
     if (ASSOCIATED(Ice_Ocean_Boundary%calving)) &
          call mpp_redistribute(Ice%Domain, Ice%calving, Ocean%Domain, Ice_Ocean_Boundary%calving)
     if (ASSOCIATED(Ice_Ocean_Boundary%p)) &
          call mpp_redistribute(Ice%Domain, Ice%p_surf, Ocean%Domain, Ice_Ocean_Boundary%p)

! Extra fluxes
     do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
       do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
         if ( associated(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values) ) then  !{
           call mpp_redistribute(Ice%Domain, Ice%ocean_fluxes%bc(n)%field(m)%values,  &
                Ocean%Domain, Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values)
         endif  !}
       enddo  !} m
     enddo  !} n

  case DEFAULT
!   <ERROR MSG="Ice_Ocean_Boundary%xtype must be DIRECT or REDIST." STATUS="FATAL">
!      The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
!   </ERROR>
     call mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Ice_Ocean_Boundary%xtype must be DIRECT or REDIST.' )
  end select
!Balaji: moved data_override calls here from coupler_main
  if( ocn_pe )then
      call mpp_set_current_pelist(ocn_pelist)
      call data_override('OCN', 'u_flux',    Ice_Ocean_Boundary%u_flux   , Time )
      call data_override('OCN', 'v_flux',    Ice_Ocean_Boundary%v_flux   , Time )
      call data_override('OCN', 't_flux',    Ice_Ocean_Boundary%t_flux   , Time )
      call data_override('OCN', 'q_flux',    Ice_Ocean_Boundary%q_flux   , Time )
      call data_override('OCN', 'salt_flux', Ice_Ocean_Boundary%salt_flux, Time )
      call data_override('OCN', 'lw_flux',   Ice_Ocean_Boundary%lw_flux  , Time )
      call data_override('OCN', 'sw_flux',   Ice_Ocean_Boundary%sw_flux  , Time )
      call data_override('OCN', 'sw_flux_vis',   Ice_Ocean_Boundary%sw_flux_vis  , Time )
      call data_override('OCN', 'sw_flux_dir',   Ice_Ocean_Boundary%sw_flux_dir  , Time )
      call data_override('OCN', 'sw_flux_dif',   Ice_Ocean_Boundary%sw_flux_dif  , Time )
      call data_override('OCN', 'sw_flux_vis_dir',   Ice_Ocean_Boundary%sw_flux_vis_dir  , Time )
      call data_override('OCN', 'sw_flux_vis_dif',   Ice_Ocean_Boundary%sw_flux_vis_dif  , Time )
      call data_override('OCN', 'lprec',     Ice_Ocean_Boundary%lprec    , Time )
      call data_override('OCN', 'fprec',     Ice_Ocean_Boundary%fprec    , Time )
      call data_override('OCN', 'runoff',    Ice_Ocean_Boundary%runoff   , Time )
      call data_override('OCN', 'calving',   Ice_Ocean_Boundary%calving  , Time )
      call data_override('OCN', 'p',         Ice_Ocean_Boundary%p        , Time )

! Extra fluxes
      do n = 1, Ice_Ocean_Boundary%fluxes%num_bcs  !{
         do m = 1, Ice_Ocean_Boundary%fluxes%bc(n)%num_fields  !{
             call data_override('OCN', Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%name,   &
                  Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Time)
            if ( Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%id_diag > 0 ) then  !{
              used = send_data(Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%id_diag,        &
                   Ice_Ocean_Boundary%fluxes%bc(n)%field(m)%values, Time )
            endif  !}
         enddo  !} m
      enddo  !} n
   endif
   call mpp_set_current_pelist()

!Balaji
  call mpp_clock_end(fluxIceOceanClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ice_to_ocean
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="flux_ocean_to_ice">
!  <OVERVIEW>
!   Takes the ocean model state and interpolates it onto the bottom of the ice. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following quantities are transferred from the Ocean to the ocean_ice_boundary_type: 
!
!        t_surf = surface temperature (deg K)
!        frazil = frazil (???)
!        u_surf = zonal ocean current/ice motion (m/s)
!        v_surf = meridional ocean current/ice motion (m/s
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_ocean_to_ice ( Time, Ocean, Ice, Ocean_Ice_Boundary)
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME=" Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Ocean_Ice_Boundary" TYPE="ocean_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!
subroutine flux_ocean_to_ice ( Time, Ocean, Ice, Ocean_Ice_Boundary )
  type(time_type),       intent(in)  :: Time
  type(ocean_data_type), intent(in)  :: Ocean
  type(ice_data_type),   intent(in)  :: Ice
!  real, dimension(:,:),   intent(out) :: t_surf_ice, u_surf_ice, v_surf_ice, &
!                                         frazil_ice, s_surf_ice, sea_lev_ice
  type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary
  real, dimension(size(Ocean_Ice_Boundary%t,1),size(Ocean_Ice_Boundary%t,2),size(Ice%part_size,3)) &
       :: ice_frac
  real, dimension(:), allocatable :: ex_ice_frac
  real, dimension(ni_atm, nj_atm) :: diag_atm
  logical :: used
  integer       :: m
  integer       :: n

!Balaji
  call mpp_clock_begin(cplOcnClock)
  call mpp_clock_begin(fluxOceanIceClock)

  do n = 1, Ocean%fields%num_bcs  !{
    do m = 1, Ocean%fields%bc(n)%num_fields  !{
      if ( Ocean%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
        used = send_data(Ocean%fields%bc(n)%field(m)%id_diag, Ocean%fields%bc(n)%field(m)%values, Time)
      endif  !}
    enddo  !} m
  enddo  !} n

  select case (Ocean_Ice_Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Ocean_Ice_Boundary%u) )Ocean_Ice_Boundary%u = Ocean%u_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%v) )Ocean_Ice_Boundary%v = Ocean%v_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%t) )Ocean_Ice_Boundary%t = Ocean%t_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%s) )Ocean_Ice_Boundary%s = Ocean%s_surf
     if( ASSOCIATED(Ocean_Ice_Boundary%frazil) )Ocean_Ice_Boundary%frazil = Ocean%frazil
     if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )Ocean_Ice_Boundary%sea_level = Ocean%sea_lev

! Extra fluxes
     do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
       do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
         if ( associated(Ocean_Ice_Boundary%fields%bc(n)%field(m)%values) ) then  !{
           Ocean_Ice_Boundary%fields%bc(n)%field(m)%values = Ocean%fields%bc(n)%field(m)%values
         endif  !}
       enddo  !} m
     enddo  !} n
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if( ASSOCIATED(Ocean_Ice_Boundary%u) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%u_surf, Ice%Domain, Ocean_Ice_Boundary%u)
     if( ASSOCIATED(Ocean_Ice_Boundary%v) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%v_surf, Ice%Domain, Ocean_Ice_Boundary%v)
     if( ASSOCIATED(Ocean_Ice_Boundary%t) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%t_surf, Ice%Domain, Ocean_Ice_Boundary%t)
     if( ASSOCIATED(Ocean_Ice_Boundary%s) )                     &
          call mpp_redistribute(Ocean%Domain, Ocean%s_surf, Ice%Domain, Ocean_Ice_Boundary%s)
     if( ASSOCIATED(Ocean_Ice_Boundary%frazil) )                &
          call mpp_redistribute(Ocean%Domain, Ocean%frazil, Ice%Domain, Ocean_Ice_Boundary%frazil)
     if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )             &
          call mpp_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%Domain, Ocean_Ice_Boundary%sea_level)

! Extra fluxes
     do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
       do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
         if ( associated(Ocean_Ice_Boundary%fields%bc(n)%field(m)%values) ) then  !{
           call mpp_redistribute(Ocean%Domain, Ocean%fields%bc(n)%field(m)%values,    &
                Ice%Domain, Ocean_Ice_Boundary%fields%bc(n)%field(m)%values)
         endif  !}
       enddo  !} m
     enddo  !} n
  case DEFAULT
!   <ERROR MSG="Ocean_Ice_Boundary%xtype must be DIRECT or REDIST." STATUS="FATAL">
!     The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
!   </ERROR>
     call mpp_error( FATAL, 'FLUX_OCEAN_TO_ICE: Ocean_Ice_Boundary%xtype must be DIRECT or REDIST.' )
  end select
  if( ice_pe )then
      call mpp_set_current_pelist(ice_pelist)

!Balaji: data_override moved here from coupler_main
      call data_override('ICE', 'u',         Ocean_Ice_Boundary%u,         Time)
      call data_override('ICE', 'v',         Ocean_Ice_Boundary%v,         Time)
      call data_override('ICE', 't',         Ocean_Ice_Boundary%t,         Time)
      call data_override('ICE', 's',         Ocean_Ice_Boundary%s,         Time)
      call data_override('ICE', 'frazil',    Ocean_Ice_Boundary%frazil,    Time)
      call data_override('ICE', 'sea_level', Ocean_Ice_Boundary%sea_level, Time)

! Extra fluxes
      do n = 1, Ocean_Ice_Boundary%fields%num_bcs  !{
         do m = 1, Ocean_Ice_Boundary%fields%bc(n)%num_fields  !{
            call data_override('ICE', Ocean_Ice_Boundary%fields%bc(n)%field(m)%name,    &
                 Ocean_Ice_Boundary%fields%bc(n)%field(m)%values, Time)
            if ( Ocean_Ice_Boundary%fields%bc(n)%field(m)%id_diag > 0 ) then  !{
              used = send_data(Ocean_Ice_Boundary%fields%bc(n)%field(m)%id_diag,        &
                   Ocean_Ice_Boundary%fields%bc(n)%field(m)%values, Time )
            endif  !}
         enddo  !} m
      enddo  !} n
   endif
   call mpp_set_current_pelist()
  
  if ( id_ice_mask > 0 ) then
     allocate ( ex_ice_frac(n_xgrid_sfc) )
     ice_frac        = 1.
     ice_frac(:,:,1) = 0.
     ex_ice_frac     = 0.
     call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
     call get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)
     used = send_data ( id_ice_mask, diag_atm, Time )
     deallocate ( ex_ice_frac )
  endif

!Balaji
  call mpp_clock_end(fluxOceanIceClock)
  call mpp_clock_end(cplOcnClock)
!-----------------------------------------------------------------------

  end subroutine flux_ocean_to_ice
! </SUBROUTINE>
