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

module ice_ocean_flux_exchange_mod

!! FMS
  use FMS
  use FMSconstants, only: HLF, HLV, CP_OCEAN
!! Components
  use ice_model_mod,       only: ice_data_type, ocean_ice_boundary_type
  use ocean_model_mod,     only: ocean_public_type, ice_ocean_boundary_type
  use ocean_model_mod,     only: ocean_state_type, ocean_model_data_get
  use ocean_model_mod,     only: ocean_model_init_sfc

  implicit none ; private


  public :: ice_ocean_flux_exchange_init, &
       flux_ice_to_ocean, flux_ice_to_ocean_finish, &
       flux_ocean_to_ice, flux_ocean_to_ice_finish, &
       flux_ice_to_ocean_stocks,&
       flux_ocean_from_ice_stocks

  !Balaji, sets boundary_type%xtype
  !  REGRID: grids are physically different, pass via exchange grid
  !  REDIST: same physical grid, different decomposition, must move data around
  !  DIRECT: same physical grid, same domain decomposition, can directly copy data
  integer, parameter :: REGRID=1, REDIST=2, DIRECT=3

  logical :: debug_stocks = .false.
  logical :: do_area_weighted_flux = .false.

  integer :: cplOcnClock, fluxOceanIceClock, fluxIceOceanClock
  real    :: Dt_cpl
  integer, allocatable :: slow_ice_ocean_pelist(:)

contains

  subroutine ice_ocean_flux_exchange_init(Time, Ice, Ocean, Ocean_state, ice_ocean_boundary, &
                                          ocean_ice_boundary, Dt_cpl_in, debug_stocks_in,    &
                                          do_area_weighted_flux_in, ex_gas_fields_ice, ex_gas_fluxes, &
                                          do_ocean, slow_ice_ocean_pelist_in )

    type(time_type),               intent(in)    :: Time !< The model's current time
    type(ice_data_type),           intent(inout) :: Ice !< A derived data type to specify ice boundary data
    type(ocean_public_type),       intent(inout) :: Ocean !< A derived data type to specify ocean boundary data
    type(ocean_state_type),        pointer       :: Ocean_state
    type(ice_ocean_boundary_type), intent(inout) :: ice_ocean_boundary !< A derived data type to specify properties and fluxes passed from ice to ocean
    type(ocean_ice_boundary_type), intent(inout) :: ocean_ice_boundary !< A derived data type to specify properties and fluxes passed from ocean to ice
    real,                          intent(in)    :: Dt_cpl_in
    logical,                       intent(in)    :: debug_stocks_in
    logical,                       intent(in)    :: do_area_weighted_flux_in
    type(coupler_1d_bc_type),      intent(in)    :: ex_gas_fields_ice, ex_gas_fluxes
    logical,                       intent(in)    :: do_ocean
    integer, dimension(:),         intent(in)    :: slow_ice_ocean_pelist_in
    integer              :: is, ie, js, je

    Dt_cpl = Dt_cpl_in
    debug_stocks = debug_stocks_in
    do_area_weighted_flux = do_area_weighted_flux_in

    !ocean_ice_boundary and ice_ocean_boundary must be done on all PES
    !domain boundaries will assure no space is allocated on non-relevant PEs.
    call mpp_get_compute_domain( Ice%slow_Domain_NH, is, ie, js, je )
    !allocate ocean_ice_boundary
    allocate( ocean_ice_boundary%u(is:ie,js:je) )
    allocate( ocean_ice_boundary%v(is:ie,js:je) )
    allocate( ocean_ice_boundary%t(is:ie,js:je) )
    allocate( ocean_ice_boundary%s(is:ie,js:je) )
    !frazil and sea_level are optional, if not present they should be nullified
    allocate( ocean_ice_boundary%frazil(is:ie,js:je) )
    allocate( ocean_ice_boundary%sea_level(is:ie,js:je) )
    ! initialize boundary fields for override experiments (mjh)
    ocean_ice_boundary%u=0.0
    ocean_ice_boundary%v=0.0
    ocean_ice_boundary%t=273.0
    ocean_ice_boundary%s=0.0
    ocean_ice_boundary%frazil=0.0
    ocean_ice_boundary%sea_level=0.0

    ! allocate fields for extra tracers in ocean_ice_boundary
    if (.not.coupler_type_initialized(ocean_ice_boundary%fields)) &
      call coupler_type_spawn(ex_gas_fields_ice, ocean_ice_boundary%fields, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix='_ocn_ice')
    if (Ice%pe) &
      call coupler_type_set_diags(ocean_ice_boundary%fields, "ice_flux", Ice%axes(1:2), Time)

    !
    ! allocate fields and fluxes for extra tracers for the Ice type
    !
    if (.not.coupler_type_initialized(Ice%ocean_fluxes)) &
      call coupler_type_spawn(ex_gas_fluxes, Ice%ocean_fluxes, (/is,is,ie,ie/), &
                              (/js,js,je,je/),  suffix = '_ice')

    ! This was never being sent, so comment it out for now.
    ! if (Ice%pe) &
    !   call coupler_type_set_diags(Ice%ocean_fluxes, "ice_flux", Ice%axes(1:2), Time)


    !allocate ice_ocean_boundary
    call mpp_get_compute_domain( Ocean%domain, is, ie, js, je )
    !ML ocean only requires t, q, lw, sw, fprec, calving
    !AMIP ocean needs no input fields
    !choice of fields will eventually be done at runtime
    !via field_manager
    allocate( ice_ocean_boundary%u_flux   (is:ie,js:je) ) ;         ice_ocean_boundary%u_flux = 0.0
    allocate( ice_ocean_boundary%v_flux   (is:ie,js:je) ) ;         ice_ocean_boundary%v_flux = 0.0
    allocate( ice_ocean_boundary%t_flux   (is:ie,js:je) ) ;         ice_ocean_boundary%t_flux = 0.0
    allocate( ice_ocean_boundary%q_flux   (is:ie,js:je) ) ;         ice_ocean_boundary%q_flux = 0.0
    allocate( ice_ocean_boundary%salt_flux(is:ie,js:je) ) ;         ice_ocean_boundary%salt_flux = 0.0
    allocate( ice_ocean_boundary%lw_flux  (is:ie,js:je) ) ;         ice_ocean_boundary%lw_flux = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dir  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_vis_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_vis_dif  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_vis_dif = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dir  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_nir_dir = 0.0
    allocate( ice_ocean_boundary%sw_flux_nir_dif  (is:ie,js:je) ) ; ice_ocean_boundary%sw_flux_nir_dif = 0.0
    allocate( ice_ocean_boundary%lprec    (is:ie,js:je) ) ;         ice_ocean_boundary%lprec = 0.0
    allocate( ice_ocean_boundary%fprec    (is:ie,js:je) ) ;         ice_ocean_boundary%fprec = 0.0
    allocate( ice_ocean_boundary%runoff   (is:ie,js:je) ) ;         ice_ocean_boundary%runoff = 0.0
    allocate( ice_ocean_boundary%calving  (is:ie,js:je) ) ;         ice_ocean_boundary%calving = 0.0
    allocate( ice_ocean_boundary%runoff_hflx   (is:ie,js:je) ) ;    ice_ocean_boundary%runoff_hflx = 0.0
    allocate( ice_ocean_boundary%calving_hflx  (is:ie,js:je) ) ;    ice_ocean_boundary%calving_hflx = 0.0
    allocate( ice_ocean_boundary%p        (is:ie,js:je) ) ;         ice_ocean_boundary%p = 0.0
    allocate( ice_ocean_boundary%mi       (is:ie,js:je) ) ;         ice_ocean_boundary%mi = 0.0
    !Allocating iceberg fields, if the corresponding fields are assosiated in the sea ice model(s)
    if (associated(Ice%ustar_berg)) then
      allocate( ice_ocean_boundary%ustar_berg (is:ie,js:je) ) ;     ice_ocean_boundary%ustar_berg = 0.0
    endif
    if (associated(Ice%area_berg)) then
      allocate( ice_ocean_boundary%area_berg  (is:ie,js:je) ) ;     ice_ocean_boundary%area_berg = 0.0
    endif
    if (associated(Ice%mass_berg)) then
      allocate( ice_ocean_boundary%mass_berg  (is:ie,js:je) ) ;     ice_ocean_boundary%mass_berg = 0.0
    endif
    ! Copy the stagger indication variables from the ice processors the ocean
    ! PEs and vice versa.  The defaults are large negative numbers, so the
    ! global max here picks out only values that have been set on active PEs.
    call mpp_max(Ice%flux_uv_stagger)
    call mpp_max(Ocean%stagger)
    ice_ocean_boundary%wind_stagger = Ice%flux_uv_stagger
    if(do_ocean) then
       ocean_ice_boundary%stagger = Ocean%stagger
    else
       ocean_ice_boundary%stagger = AGRID
    endif

    ! allocate fields for extra tracer fluxes in ice_ocean_boundary
    if (.not.coupler_type_initialized(ice_ocean_boundary%fluxes)) &
      call coupler_type_spawn(ex_gas_fluxes, ice_ocean_boundary%fluxes, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix='_ice_ocn')
    if (Ocean%is_ocean_pe) &
      call coupler_type_set_diags(ice_ocean_boundary%fluxes, "ocean_flux", Ocean%axes(1:2), Time)

    ! This typically only occurs on non-ocean PEs.
    if (.not.coupler_type_initialized(Ocean%fields)) &
      call coupler_type_spawn(ex_gas_fields_ice, Ocean%fields, (/is,is,ie,ie/), &
                              (/js,js,je,je/), suffix = '_ocn')

    ! initialize boundary values for override experiments
    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%slow_Domain_NH )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype

    !       initialize the Ocean type for extra fields for surface fluxes
    ! Same allocation of arrays and stuff
    !       (this must be done after the Ocean fields are allocated as the fields on the Ocean%fields
    !       are read in in this subroutine)
    !

    if ( Ocean%is_ocean_pe ) then
       call mpp_set_current_pelist(Ocean%pelist)
       call ocean_model_init_sfc(Ocean_state, Ocean)
    endif
    call mpp_set_current_pelist()

    !z1l check the flux conservation.
    if(debug_stocks) call check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      allocate(slow_ice_ocean_pelist(size(slow_ice_ocean_pelist_in(:))))
      slow_ice_ocean_pelist = slow_ice_ocean_pelist_in
      call mpp_set_current_pelist(slow_ice_ocean_pelist)
      cplOcnClock = mpp_clock_id( 'Ice-ocean coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
      fluxIceOceanClock = mpp_clock_id( 'Flux ice to ocean', flags=clock_flag_default, grain=CLOCK_ROUTINE )
      fluxOceanIceClock = mpp_clock_id( 'Flux ocean to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    endif

  end subroutine ice_ocean_flux_exchange_init


  !#######################################################################
  !! \brief Takes the ice model state (fluxes at the bottom of the ice) and interpolates it to the ocean model grid.
  !!
  !! The following quantities are transferred from the Ice to the ice_ocean_boundary_type:
  !! <pre>
  !!       flux_u = zonal wind stress (Pa)
  !!       flux_v = meridional wind stress (Pa)
  !!       flux_t = sensible heat flux (W/m2)
  !!       flux_q = specific humidity flux (Kg/m2/s)
  !!    flux_salt = salt flux (Kg/m2/s)
  !!      flux_sw = net (down-up) shortwave flux (W/m2)
  !!      flux_lw = net (down-up) longwave flux (W/m2)
  !!        lprec = mass of liquid precipitation since last
  !!                  time step (Kg/m2)
  !!        fprec = mass of frozen precipitation since last
  !!                  time step (Kg/m2)
  !!       runoff = mass of runoff since last time step (Kg/m2)
  !!       runoff = mass of calving since last time step (Kg/m2)
  !!       p_surf = surface pressure (Pa)
  !! </pre>
  subroutine flux_ice_to_ocean ( Time, Ice, Ocean, Ice_Ocean_Boundary )

    type(time_type),                 intent(in)  :: Time !< Current time
    type(ice_data_type),             intent(in)  :: Ice  !< A derived data type to specify ice boundary data
    type(ocean_public_type),         intent(in)  :: Ocean !< A derived data type to specify ocean boundary data
    type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary !< A derived data type to specify properties and fluxes
                                                         !! passed from ice to ocean

    integer       :: m
    integer       :: n
    logical       :: used

    call mpp_clock_begin(cplOcnClock)
    call mpp_clock_begin(fluxIceOceanClock)

    if(ASSOCIATED(Ice_Ocean_Boundary%u_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_u, Ice_Ocean_Boundary%u_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%v_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_v, Ice_Ocean_Boundary%v_flux, Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%p     ) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%p_surf, Ice_Ocean_Boundary%p     , Ice_Ocean_Boundary%xtype, .FALSE. )

    if(ASSOCIATED(Ice_Ocean_Boundary%mi    ) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%mi,     Ice_Ocean_Boundary%mi    , Ice_Ocean_Boundary%xtype, .FALSE. )

    ! Extra fluxes
    if (Ice_Ocean_Boundary%xtype == DIRECT) then
       call coupler_type_copy_data(Ice%ocean_fluxes, Ice_Ocean_Boundary%fluxes)
    else
       call coupler_type_redistribute_data(Ice%ocean_fluxes, Ice%slow_Domain_NH, &
                     Ice_Ocean_Boundary%fluxes, ocean%Domain, complete=.true.)
    endif

    !--- The following variables may require conserved flux exchange from ice to ocean because the
    !--- ice area maybe different from ocean area.
    if(ASSOCIATED(Ice_Ocean_Boundary%t_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_t, Ice_Ocean_Boundary%t_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%salt_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_salt, Ice_Ocean_Boundary%salt_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_nir_dir, Ice_Ocean_Boundary%sw_flux_nir_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_nir_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_nir_dif, Ice_Ocean_Boundary%sw_flux_nir_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dir) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_vis_dir, Ice_Ocean_Boundary%sw_flux_vis_dir, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%sw_flux_vis_dif) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_sw_vis_dif, Ice_Ocean_Boundary%sw_flux_vis_dif, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%lw_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_lw, Ice_Ocean_Boundary%lw_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%lprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%lprec, Ice_Ocean_Boundary%lprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%fprec) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%fprec, Ice_Ocean_Boundary%fprec, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%runoff) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%runoff, Ice_Ocean_Boundary%runoff, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%calving) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%calving, Ice_Ocean_Boundary%calving, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%ustar_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
       Ice%ustar_berg, Ice_Ocean_Boundary%ustar_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%area_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%area_berg, Ice_Ocean_Boundary%area_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%mass_berg) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%mass_berg, Ice_Ocean_Boundary%mass_berg, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%runoff_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%runoff_hflx, Ice_Ocean_Boundary%runoff_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%calving_hflx) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%calving_hflx, Ice_Ocean_Boundary%calving_hflx, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    if(ASSOCIATED(Ice_Ocean_Boundary%q_flux) ) call flux_ice_to_ocean_redistribute( Ice, Ocean, &
         Ice%flux_q, Ice_Ocean_Boundary%q_flux, Ice_Ocean_Boundary%xtype, do_area_weighted_flux )

    call mpp_clock_end(fluxIceOceanClock)
    call mpp_clock_end(cplOcnClock)
    !-----------------------------------------------------------------------

  end subroutine flux_ice_to_ocean

  !> flux_ice_to_ocean_finish carrries out a final set of tasks that should only occur on
  !! the ocean processors, including data override and perhaps saving diagnostics.
  subroutine flux_ice_to_ocean_finish ( Time, Ice_Ocean_Boundary )

    type(time_type),                 intent(in)  :: Time !< Current time
    type(ice_ocean_boundary_type), intent(inout) :: Ice_Ocean_Boundary !< A derived data type to specify properties and fluxes
                                                         !! passed from ice to ocean

    call data_override('OCN', 'u_flux',    Ice_Ocean_Boundary%u_flux   , Time )
    call data_override('OCN', 'v_flux',    Ice_Ocean_Boundary%v_flux   , Time )
    call data_override('OCN', 't_flux',    Ice_Ocean_Boundary%t_flux   , Time )
    call data_override('OCN', 'q_flux',    Ice_Ocean_Boundary%q_flux   , Time )
    call data_override('OCN', 'salt_flux', Ice_Ocean_Boundary%salt_flux, Time )
    call data_override('OCN', 'lw_flux',   Ice_Ocean_Boundary%lw_flux  , Time )
    call data_override('OCN', 'sw_flux_nir_dir', Ice_Ocean_Boundary%sw_flux_nir_dir, Time )
    call data_override('OCN', 'sw_flux_nir_dif', Ice_Ocean_Boundary%sw_flux_nir_dif, Time )
    call data_override('OCN', 'sw_flux_vis_dir', Ice_Ocean_Boundary%sw_flux_vis_dir, Time )
    call data_override('OCN', 'sw_flux_vis_dif', Ice_Ocean_Boundary%sw_flux_vis_dif, Time )
    call data_override('OCN', 'lprec',     Ice_Ocean_Boundary%lprec    , Time )
    call data_override('OCN', 'fprec',     Ice_Ocean_Boundary%fprec    , Time )
    call data_override('OCN', 'runoff',    Ice_Ocean_Boundary%runoff   , Time )
    call data_override('OCN', 'calving',   Ice_Ocean_Boundary%calving  , Time )
    call data_override('OCN', 'runoff_hflx',    Ice_Ocean_Boundary%runoff_hflx   , Time )
    call data_override('OCN', 'calving_hflx',   Ice_Ocean_Boundary%calving_hflx  , Time )
    call data_override('OCN', 'p',         Ice_Ocean_Boundary%p        , Time )
    call data_override('OCN', 'mi',        Ice_Ocean_Boundary%mi       , Time )

   !Are these if statements needed, or does data_override routine check if variable is associated?
    if (ASSOCIATED(Ice_Ocean_Boundary%ustar_berg) ) &
      call data_override('OCN', 'ustar_berg', Ice_Ocean_Boundary%ustar_berg, Time )
    if (ASSOCIATED(Ice_Ocean_Boundary%area_berg)  ) &
      call data_override('OCN', 'area_berg',  Ice_Ocean_Boundary%area_berg , Time )
    if (ASSOCIATED(Ice_Ocean_Boundary%mass_berg)  ) &
      call data_override('OCN', 'mass_berg',  Ice_Ocean_Boundary%mass_berg , Time )

    ! Extra fluxes
    call coupler_type_data_override('OCN', Ice_Ocean_Boundary%fluxes, Time )

    ! The send_data call for any of the Ice_Ocean_Boundary fluxes would go here.
    call coupler_type_send_data(Ice_Ocean_Boundary%fluxes, Time )

  end subroutine flux_ice_to_ocean_finish

  !#######################################################################
  !> \brief Takes the ocean model state and interpolates it onto the bottom of the ice.
  !!
  !! The following quantities are transferred from the Ocean to the ocean_ice_boundary_type:
  !! <pre>
  !!        t_surf = surface temperature (deg K)
  !!        frazil = frazil fluxes since the last coupling step (J/m2)
  !!        u_surf = zonal ocean current/ice motion (m/s)
  !!        v_surf = meridional ocean current/ice motion (m/s)
  !!        v_surf = meridional ocean current/ice motion (m/s)
  !!       sea_lev = sea level used to drive ice accelerations (m)
  !! </pre>
  !!
  !! \throw FATAL, "Ocean_Ice_Boundary%xtype must be DIRECT or REDIST."
  !!    The value of variable xtype of ice_ocean_boundary_type data must be DIRECT or REDIST.
  subroutine flux_ocean_to_ice ( Time, Ocean, Ice, Ocean_Ice_Boundary )

    type(time_type),                 intent(in)  :: Time  !< Current time
    type(ocean_public_type),         intent(in)  :: Ocean !< A derived data type to specify ocean boundary data
    type(ice_data_type),             intent(in)  :: Ice   !< A derived data type to specify ice boundary data
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary !< A derived data type to specify properties and fluxes
                                                          !! passed from ocean to ice
    real, allocatable, dimension(:,:) :: tmp
    integer       :: m
    integer       :: n
    logical       :: used

    call mpp_clock_begin(cplOcnClock)
    call mpp_clock_begin(fluxOceanIceClock)

    select case (Ocean_Ice_Boundary%xtype)
    case(DIRECT)
       !same grid and domain decomp for ocean and ice
       if( ASSOCIATED(Ocean_Ice_Boundary%u) )Ocean_Ice_Boundary%u = Ocean%u_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%v) )Ocean_Ice_Boundary%v = Ocean%v_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%t) )Ocean_Ice_Boundary%t = Ocean%t_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%s) )Ocean_Ice_Boundary%s = Ocean%s_surf
       if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )Ocean_Ice_Boundary%sea_level = Ocean%sea_lev
       if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
          if(do_area_weighted_flux) then
             Ocean_Ice_Boundary%frazil = Ocean%frazil * Ocean%area
             call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
          else
             Ocean_Ice_Boundary%frazil = Ocean%frazil
          endif
       endif

       ! Extra fluxes
       call coupler_type_copy_data(Ocean%fields, Ocean_Ice_Boundary%fields)

    case(REDIST)
       !same grid, different domain decomp for ocean and ice
       if( ASSOCIATED(Ocean_Ice_Boundary%u) )                     &
            call mpp_redistribute(Ocean%Domain, Ocean%u_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%u)
       if( ASSOCIATED(Ocean_Ice_Boundary%v) )                     &
            call mpp_redistribute(Ocean%Domain, Ocean%v_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%v)
       if( ASSOCIATED(Ocean_Ice_Boundary%t) )                     &
            call mpp_redistribute(Ocean%Domain, Ocean%t_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%t)
       if( ASSOCIATED(Ocean_Ice_Boundary%s) )                     &
            call mpp_redistribute(Ocean%Domain, Ocean%s_surf, Ice%slow_Domain_NH, Ocean_Ice_Boundary%s)

       if( ASSOCIATED(Ocean_Ice_Boundary%sea_level) )             &
            call mpp_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%slow_Domain_NH, Ocean_Ice_Boundary%sea_level)

       if( ASSOCIATED(Ocean_Ice_Boundary%frazil) ) then
          if(do_area_weighted_flux) then
             if (Ocean%is_ocean_pe) then
               allocate(tmp(size(Ocean%area,1), size(Ocean%area,2)))
               tmp(:,:) = Ocean%frazil(:,:) * Ocean%area(:,:)
             endif
             call mpp_redistribute( Ocean%Domain, tmp, Ice%slow_Domain_NH, Ocean_Ice_Boundary%frazil)
             if (Ice%slow_ice_pe) &
               call divide_by_area(data=Ocean_Ice_Boundary%frazil, area=Ice%area)
             if (Ocean%is_ocean_pe) deallocate(tmp)
          else
             call mpp_redistribute(Ocean%Domain, Ocean%frazil, Ice%slow_Domain_NH, Ocean_Ice_Boundary%frazil)
          endif
       endif

       ! Extra fluxes
       call coupler_type_redistribute_data(Ocean%fields, Ocean%Domain, &
                     Ocean_Ice_Boundary%fields, Ice%slow_Domain_NH)
    case DEFAULT
       call mpp_error( FATAL, 'flux_ocean_to_ice: Ocean_Ice_Boundary%xtype must be DIRECT or REDIST.' )
    end select

    call mpp_clock_end(fluxOceanIceClock)
    call mpp_clock_end(cplOcnClock)
    !-----------------------------------------------------------------------

  end subroutine flux_ocean_to_ice

  !> flux_ocean_to_ice_finish carrries out a final set of tasks that should only occur on
  !! the slow-ice processors, including data override and perhaps saving diagnostics.
  subroutine flux_ocean_to_ice_finish( Time, Ice, Ocean_Ice_Boundary )

    type(time_type),                 intent(in)  :: Time  !< Current time
    type(ice_data_type),             intent(in)  :: Ice   !< A derived data type to specify ice boundary data
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_Ice_Boundary !< A derived data type to specify properties and fluxes
                                                          !! passed from ocean to ice
    real          :: from_dq

    call data_override('ICE', 'u',         Ocean_Ice_Boundary%u,         Time)
    call data_override('ICE', 'v',         Ocean_Ice_Boundary%v,         Time)
    call data_override('ICE', 't',         Ocean_Ice_Boundary%t,         Time)
    call data_override('ICE', 's',         Ocean_Ice_Boundary%s,         Time)
    call data_override('ICE', 'frazil',    Ocean_Ice_Boundary%frazil,    Time)
    call data_override('ICE', 'sea_level', Ocean_Ice_Boundary%sea_level, Time)
    call coupler_type_data_override('ICE', Ocean_Ice_Boundary%fields, Time)

    !  Perform diagnostic output for the ocean_ice_boundary fields
    call coupler_type_send_data( Ocean_Ice_Boundary%fields, Time)

    ! frazil (already in J/m^2 so no need to multiply by Dt_cpl)
    from_dq = SUM( Ice%area * Ocean_Ice_Boundary%frazil )
    Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

  end subroutine flux_ocean_to_ice_finish


  !#######################################################################
  !> \brief  Updates Ice and Ocean stocks.
  !!
  !!   Integrate the fluxes over the surface and in time.

  subroutine flux_ice_to_ocean_stocks(Ice)

    type(ice_data_type),   intent(in)  :: Ice !< A derived data type to specify ice boundary data

    real           :: from_dq

    ! fluxes from ice -> ocean, integrate over surface and in time

    ! precip - evap
    from_dq = Dt_cpl * SUM( Ice%area * (Ice%lprec+Ice%fprec-Ice%flux_q) )
    Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_TOP   ) + from_dq

    ! river
    from_dq = Dt_cpl * SUM( Ice%area * (Ice%runoff + Ice%calving) )
    Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_WATER)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_WATER)%dq(ISTOCK_SIDE  ) + from_dq

    ! sensible heat + shortwave + longwave + latent heat
    from_dq = Dt_cpl * SUM( Ice%area * ( &
         &   Ice%flux_sw_vis_dir+Ice%flux_sw_vis_dif &
         & + Ice%flux_sw_nir_dir+Ice%flux_sw_nir_dif + Ice%flux_lw &
         & - (Ice%fprec + Ice%calving)*HLF - Ice%flux_t - Ice%flux_q*HLV) )
    Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    ! heat carried by river + pme (assuming reference temperature of 0 degC and river/pme temp = surface temp)
    ! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE
    from_dq = Dt_cpl * SUM( Ice%area * ( &
         & (Ice%lprec+Ice%fprec-Ice%flux_q + Ice%runoff+Ice%calving)*CP_OCEAN*Ice%SST_C(:,:)) )
    Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_HEAT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq(ISTOCK_SIDE  ) + from_dq

    !SALT flux
    from_dq = Dt_cpl* SUM( Ice%area * ( -Ice%flux_salt ))
    Ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) = Ice_stock(ISTOCK_SALT)%dq(ISTOCK_BOTTOM) - from_dq
    Ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_SALT)%dq(ISTOCK_TOP   ) + from_dq


  end subroutine flux_ice_to_ocean_stocks

  !#######################################################################
  !> \brief  Updates Ocean stocks due to input that the Ocean model gets.
  !!
  !! This subroutine updates the stocks of Ocean by the amount of input that the Ocean gets from Ice component.
  !! Unlike subroutine flux_ice_to_ocean_stocks() that uses Ice%fluxes to update the stocks due to the amount of output from Ice
  !! this subroutine uses Ice_Ocean_boundary%fluxes to calculate the amount of input to the Ocean. These fluxes are the ones
  !! that Ocean model uses internally to calculate its budgets. Hence there should be no difference between this input and what
  !! Ocean model internal diagnostics uses.
  !! This bypasses the possible mismatch in cell areas between Ice and Ocean in diagnosing the stocks of Ocean
  !! and should report a conserving Ocean component regardless of the glitches in fluxes.
  !!
  !! The use of this subroutine in conjunction with  subroutine flux_ice_to_ocean_stocks() will also allow to directly
  !! diagnose the amount "stocks lost in exchange" between Ice and Ocean
  subroutine flux_ocean_from_ice_stocks(ocean_state,Ocean,Ice_Ocean_boundary)
    type(ocean_state_type),        pointer    :: ocean_state
    type(ocean_public_type),       intent(in) :: Ocean
    type(ice_ocean_boundary_type), intent(in) :: Ice_Ocean_Boundary
    real    :: from_dq, cp_ocn
    real, dimension(size(Ice_Ocean_Boundary%lprec,1), size(Ice_Ocean_Boundary%lprec,2)) :: &
      ocean_cell_area, wet, t_surf, t_pme, t_calving, t_runoff, btfHeat
    integer :: isc, iec, jsc, jec

    call mpp_get_compute_domain(Ocean%Domain, isc, iec, jsc, jec)
    call ocean_model_data_get(ocean_state,Ocean,'area'  , ocean_cell_area,isc,jsc)
    call ocean_model_data_get(ocean_state,Ocean,'mask', wet,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_surf', t_surf,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_runoff', t_runoff,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_pme', t_pme,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'t_calving', t_calving,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'btfHeat', btfHeat,isc,jsc )
    call ocean_model_data_get(ocean_state,Ocean,'c_p', cp_ocn )


    ! fluxes from ice -> ocean, integrate over surface and in time

    ! precip - evap
    from_dq = SUM( ocean_cell_area * wet * (Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux) )
    Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) = Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_TOP   ) + from_dq * Dt_cpl

    from_dq = SUM( ocean_cell_area * wet * (Ice_Ocean_Boundary%runoff+Ice_Ocean_Boundary%calving) )
    Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_WATER)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    ! sensible heat + shortwave + longwave + latent heat

    from_dq = SUM( ocean_cell_area * wet *( Ice_Ocean_Boundary%sw_flux_vis_dir + Ice_Ocean_Boundary%sw_flux_vis_dif &
         +Ice_Ocean_Boundary%sw_flux_nir_dir + Ice_Ocean_Boundary%sw_flux_nir_dif &
         +Ice_Ocean_Boundary%lw_flux &
         - (Ice_Ocean_Boundary%fprec + Ice_Ocean_Boundary%calving)*HLF &
         - Ice_Ocean_Boundary%t_flux - Ice_Ocean_Boundary%q_flux*HLV ))

    Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) + from_dq * Dt_cpl

    ! heat carried by river + pme (assuming reference temperature of 0 degC and river/pme temp = surface temp)
    ! Note: it does not matter what the ref temperature is but it must be consistent with that in OCN and ICE

    from_dq = SUM( ocean_cell_area * wet * cp_ocn *&
         ((Ice_Ocean_Boundary%lprec+Ice_Ocean_Boundary%fprec-Ice_Ocean_Boundary%q_flux)*t_pme &
         +Ice_Ocean_Boundary%calving * t_calving &
         +Ice_Ocean_Boundary%runoff  * t_runoff  ))

    Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq * Dt_cpl

    !   Bottom heat flux
    from_dq = - SUM( ocean_cell_area * wet * btfHeat)

    Ocn_stock(ISTOCK_HEAT)%dq_IN( ISTOCK_BOTTOM ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_BOTTOM ) + from_dq * Dt_cpl

    !   Frazil heat

    from_dq =  SUM( ocean_cell_area *wet * Ocean%frazil )
    Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE  ) = Ocn_stock(ISTOCK_HEAT)%dq_IN(ISTOCK_SIDE ) + from_dq

    !SALT flux
    from_dq = SUM( ocean_cell_area * wet * ( -Ice_Ocean_Boundary%salt_flux))
    Ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP  ) = Ocn_stock(ISTOCK_SALT)%dq_IN(ISTOCK_TOP   ) + from_dq  * Dt_cpl


  end subroutine flux_ocean_from_ice_stocks

  !#######################################################################
  !> \brief Performs a globally conservative flux redistribution across ICE/OCN.
  !! Assumes that the ice/ocn grids are the same. If ocean is present,
  !! then assume different mpp domans and redistribute
  !!
  !! \note Should be invoked by all PEs
  subroutine flux_ice_to_ocean_redistribute(ice, ocean, ice_data, ocn_bnd_data, type, do_area_weighted )

    ! Performs a globally conservative flux redistribution across ICE/OCN.
    ! Assumes that the ice/ocn grids are the same. If ocean is present,
    ! then assume different mpp domans and redistribute

    ! should be invoked by all PEs

    type(ice_data_type),              intent(in) :: ice
    type(ocean_public_type),          intent(in) :: ocean
    real, dimension(:,:),             intent(in) :: ice_data
    real, dimension(:,:),            intent(out) :: ocn_bnd_data
    integer,                          intent(in) :: type
    logical,                          intent(in) :: do_area_weighted

    real, allocatable, dimension(:,:) :: tmp

    select case(type)
    case (DIRECT)
       if(do_area_weighted) then
          ocn_bnd_data = ice_data * ice%area
          call divide_by_area(data=ocn_bnd_data, area=ocean%area)
       else
          ocn_bnd_data = ice_data
       endif
    case (REDIST)
       if (do_area_weighted) then
         if ( Ice%slow_ice_pe ) then
            allocate(tmp(size(ice%area,1), size(ice%area,2)))
            tmp(:,:) = ice_data(:,:) * ice%area(:,:)
          endif
          call mpp_redistribute(Ice%slow_Domain_NH, tmp, ocean%Domain, ocn_bnd_data)
          if (ocean%is_ocean_pe) call divide_by_area(ocn_bnd_data, area=ocean%area)
          if (Ice%slow_ice_pe) deallocate(tmp)
       else
          call mpp_redistribute(Ice%slow_Domain_NH, ice_data, ocean%Domain, ocn_bnd_data)
       endif
    case DEFAULT
       call mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Ice_Ocean_Boundary%xtype must be DIRECT or REDIST.' )
    end select

  end subroutine flux_ice_to_ocean_redistribute

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

  !> \brief Check flux conservation for routine flux_ice_to_ocean_redistribute
  !! when do_area_weighted_flux = false and true.
  subroutine check_flux_conservation(Ice, Ocean, Ice_Ocean_Boundary)
    type(ice_data_type),               intent(inout)  :: Ice
    type(ocean_public_type),           intent(inout)  :: Ocean
    type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary

    real, allocatable, dimension(:,:) :: ice_data, ocn_data
    real :: ice_sum, area_weighted_sum, non_area_weighted_sum
    integer :: outunit

    outunit = stdout()
    allocate(ice_data(size(Ice%flux_q,1), size(Ice%flux_q,2) ) )
    allocate(ocn_data(size(Ice_Ocean_Boundary%q_flux,1), size(Ice_Ocean_Boundary%q_flux,2) ) )
    call random_number(ice_data)
    ice_sum = sum(ice_data*ice%area)
    call mpp_sum(ice_sum)
    ocn_data = 0.0
    call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .false.)
    non_area_weighted_sum = sum(ocn_data*ocean%area)
    call mpp_sum(non_area_weighted_sum)
    ocn_data = 0.0
    call flux_ice_to_ocean_redistribute( Ice, Ocean, ice_data, ocn_data, Ice_Ocean_Boundary%xtype, .true.)
    area_weighted_sum = sum(ocn_data*ocean%area)
    call mpp_sum(area_weighted_sum)
    write(outunit,*)"NOTE from flux_exchange_mod: check for flux conservation for flux_ice_to_ocean"
    write(outunit,*)"***** The global area sum of random number on ice domain (input data) is ", ice_sum
    write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
         "do_area_weighted_flux = false is ", non_area_weighted_sum, &
         " and the difference from global input area sum = ", ice_sum - non_area_weighted_sum
    write(outunit,*)"***** The global area sum of data after flux_ice_to_ocean_redistribute with "// &
         "do_area_weighted_flux = true is ", area_weighted_sum, &
         " and the difference from global input area sum = ", ice_sum - area_weighted_sum


  end subroutine check_flux_conservation

end module ice_ocean_flux_exchange_mod
