!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
! utilities used by coupler
!
  use constants_mod, only:    stefan, constants_init
  use time_manager_mod, only: time_type, set_calendar_type, set_time,  &
                              set_date, get_date, days_in_month, month_name,  &
                              operator(+), operator(-), operator (<), &
                              operator (>), operator ( /= ), operator ( / ), &
                              operator (*), thirty_day_months, julian, &
                              NOLEAP, no_calendar

  use  utilities_mod, only: open_file, file_exist, check_nml_error,  &
                            error_mesg, FATAL, WARNING,              &
                            print_version_number, get_my_pe,         &
                            utilities_init, utilities_end,           &
                            close_file, check_system_clock

  use  diag_manager_mod, only: diag_manager_init, DIAG_OCEAN, DIAG_OTHER, DIAG_ALL, diag_manager_end, get_base_date
  use  data_override_mod,only: data_override_init, data_override
!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!
  use  atmos_model_mod, only: atmos_model_init, atmos_model_end, &
                              update_atmos_model_down,           &
                              update_atmos_model_up,             &
                              atmos_data_type, &
                              land_ice_atmos_boundary_type
  use   land_model_mod, only: land_model_init, land_model_end, &
                              land_data_type, atmos_land_boundary_type, &
                              update_land_model_fast, update_land_model_slow

  use    ice_model_mod, only: ice_model_init, ice_model_end,  &
                              update_ice_model_slow_up,          &
                              update_ice_model_fast,          &
                              update_ice_model_slow_dn,          &
                              ice_data_type, land_ice_boundary_type, &
                              ocean_ice_boundary_type, atmos_ice_boundary_type

  use  ocean_model_mod, only: update_ocean_model, ocean_model_init,  &
                              ocean_model_end, ocean_data_type, ice_ocean_boundary_type, &
                              read_ice_ocean_boundary, write_ice_ocean_boundary, &
                              init_default_ice_ocean_boundary
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!
  use flux_exchange_mod, only: flux_exchange_init,   &
                               sfc_boundary_layer,   &
                               generate_sfc_xgrid,   &
                               flux_down_from_atmos, &
                               flux_up_to_atmos,     &
                               flux_land_to_ice,     &
                               flux_ice_to_ocean,    &
                               flux_ocean_to_ice
  use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, CLOCK_COMPONENT
  use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, stderr, mpp_set_current_pelist, mpp_declare_pelist, mpp_error
  use mpp_domains_mod, only: mpp_broadcast_domain, mpp_check_field
  use mpp_set_pelist_mod, only : pelist1, pelist2, atm_pelist, ocn_pelist, group1, group2, &
                                 atm1_pelist, atm2_pelist, ocn1_pelist, ocn2_pelist,       &
                                 mpp_set_pelist_init, mpp_set_pelist_end

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 1.7 2003/04/09 21:09:28 fms Exp $'
  character(len=128) :: tag = '$Name: inchon $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  type (ocean_data_type) :: Ocean

  type(atmos_land_boundary_type)     :: Atmos_land_boundary
  type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  type(land_ice_boundary_type)       :: Land_ice_boundary
  type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
  type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, Time_step_ocean, &
                      Time_step_atmos, Time_step_cpld
  type(time_type) :: Time_atmos, Time_ocean
  integer :: num_ocean_calls, num_atmos_calls, no, na
  integer :: num_cpld_calls, nc

! ----- coupled model initial date -----

  logical :: ocean_seg_start
  logical :: ocean_seg_end
  integer :: date_init(6)
  integer :: calendar_type = -99

!-----------------------------------------------------------------------

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: override = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land
  integer :: dt_ocean = 0  ! ocean tracer timestep
  integer :: dt_cpld  = 0  ! fluxes passed between ice & ocean
  integer,dimension (3)           :: locmax, locmin

  integer ::atmos_npes=0, ocean_npes=0, ice_npes=0, land_npes=0
  logical :: do_atmos =.true., do_land =.true., do_ice =.true., do_ocean=.true.
  logical :: do_flux =.true.
  logical :: concurrent=.FALSE.
  integer :: npes1 = 0, npes2 = 0
  integer :: atm1_npes = 0, atm2_npes = 0, ocn1_npes = 0, ocn2_npes = 0
  logical :: check_parallel = .FALSE.
  logical :: debug_memuse = .FALSE.
  namelist /coupler_nml/ current_date, calendar, override,       &
       months, days, hours, minutes, seconds,  &
       dt_cpld, dt_atmos, dt_ocean, &
       do_atmos, do_land, do_ice, do_ocean, do_flux, &
       atmos_npes, ocean_npes, ice_npes, land_npes, &
       concurrent, npes1, npes2, atm1_npes, ocn1_npes, &
       atm2_npes, ocn2_npes, check_parallel, debug_memuse
  integer :: initclock, atmclock, ocnclock, iceclock, landclock, fluxclock, &
       icefluxclock, ilfluxclock, fluxclockdn, fluxclockup, endclock

  integer :: flags=0
  character(len=80) :: text
  logical :: verbose=.TRUE.
  character(len=128) :: mesg = ""

!#######################################################################

  call mpp_init()
  flags = MPP_CLOCK_SYNC        !you could also choose MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED
!these clocks are on the global pelist
  initclock = mpp_clock_id( 'Initialization', flags=flags, grain=CLOCK_COMPONENT )
  endclock  = mpp_clock_id( 'Termination',    flags=flags, grain=CLOCK_COMPONENT )
  call mpp_clock_begin(initclock)
  
  call utilities_init
  call print_memuse( 'utilities_init' )
  call constants_init

  call coupler_init
  icefluxclock = mpp_clock_id( 'Ice ocean flux', flags=flags, grain=CLOCK_COMPONENT )
  call print_memuse( 'coupler_init' )

  if( Ocean%pe )then
      call mpp_set_current_pelist( Ocean%pelist )
      ocnclock     = mpp_clock_id( 'Ocean', flags=flags, grain=CLOCK_COMPONENT )
  end if
  if( Atm%pe )then
      call mpp_set_current_pelist( Atm%pelist )
      atmclock     = mpp_clock_id( 'Atmosphere',      flags=flags, grain=CLOCK_COMPONENT )
      iceclock     = mpp_clock_id( 'Ice',             flags=flags, grain=CLOCK_COMPONENT )
      landclock    = mpp_clock_id( 'Land',            flags=flags, grain=CLOCK_COMPONENT )
      fluxclock    = mpp_clock_id( 'Surface BL',      flags=flags, grain=CLOCK_COMPONENT )
      ilfluxclock  = mpp_clock_id( 'Ice land flux',   flags=flags, grain=CLOCK_COMPONENT )
      fluxclockdn  = mpp_clock_id( 'Atmos flux down', flags=flags, grain=CLOCK_COMPONENT )
      fluxclockup  = mpp_clock_id( 'Atmos flux up',   flags=flags, grain=CLOCK_COMPONENT )
  end if
  call mpp_set_current_pelist()

  call check_system_clock ('END OF INITIALIZATION')
  call mpp_clock_end (initclock)
 
 if(check_parallel) then
     if(group1) call mpp_set_current_pelist(pelist1)
     if(group2) call mpp_set_current_pelist(pelist2)
  endif

!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop ------

  do nc = 1, num_cpld_calls
     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         call generate_sfc_xgrid( Land, Ice )
     end if
     if(check_parallel) then
        if(group1) call mpp_set_current_pelist(pelist1)
        if(group2) call mpp_set_current_pelist(pelist2)
     else
        call mpp_set_current_pelist()
     endif

     write( text,'(a,i4)' )'nc=', nc
     call print_memuse(text)

   ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
   ! points when running concurrently. The calls as placed next to each other in
   ! concurrent mode to avoid multiple synchronizations within the main loop.

     call mpp_clock_begin(icefluxclock)
     call flux_ocean_to_ice( Ocean, Ice, Ocean_ice_boundary )
     call mpp_clock_end(icefluxclock)
!     call print_memuse( 'flux_ocean_to_ice' )

   ! Update Ice_ocean_boundary; first iteration is supplied by restart     
     if( concurrent .AND. nc.NE.1 )then
         call mpp_clock_begin(icefluxclock)
         call flux_ice_to_ocean( Ice, Ocean, Ice_ocean_boundary )
         call mpp_clock_end(icefluxclock)
     end if

     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         call mpp_clock_begin(iceclock)
         call ocn_ice_bnd_from_data ( Ocean_ice_boundary )
         if (do_ice) call update_ice_model_slow_up( Ocean_ice_boundary, Ice )
         call mpp_clock_end(iceclock)
!         call print_memuse( 'update_ice_model_slow_up' )

!-----------------------------------------------------------------------
!   ------ atmos/fast-land/fast-ice integration loop -------


         do na = 1, num_atmos_calls

            write( text,'(a,2i4)' )'nc,na=', nc, na
            call print_memuse(text)

            Time_atmos = Time_atmos + Time_step_atmos

            call mpp_clock_begin(fluxclock)
            call data_override ('ATM', 't_bot',  Atm%t_bot , Time_atmos)
            call data_override ('ATM', 'q_bot',  Atm%q_bot , Time_atmos)
            call data_override ('ATM', 'z_bot',  Atm%z_bot , Time_atmos)
            call data_override ('ATM', 'p_bot',  Atm%p_bot , Time_atmos)
            call data_override ('ATM', 'u_bot',  Atm%u_bot , Time_atmos)
            call data_override ('ATM', 'v_bot',  Atm%v_bot , Time_atmos)
            call data_override ('ATM', 'p_surf', Atm%p_surf, Time_atmos)
            call data_override ('ATM', 'gust',   Atm%gust,   Time_atmos)
            call data_override ('ICE', 't_surf',     Ice%t_surf,      Time_atmos)
            call data_override ('ICE', 'rough_mom',  Ice%rough_mom,   Time_atmos)
            call data_override ('ICE', 'rough_heat', Ice%rough_heat,  Time_atmos)
            call data_override ('ICE', 'rough_moist',Ice%rough_moist, Time_atmos)
            call data_override ('ICE', 'albedo',     Ice%albedo,      Time_atmos)
            call data_override ('ICE', 'u_surf',     Ice%u_surf,      Time_atmos)
            call data_override ('ICE', 'v_surf',     Ice%v_surf,      Time_atmos)
            call data_override ('LND', 't_surf',     Land%t_surf,     Time_atmos)
            call data_override ('LND', 't_ca',       Land%t_ca,       Time_atmos)
            call data_override ('LND', 'q_ca',       Land%q_ca,       Time_atmos)
            call data_override ('LND', 'rough_mom',  Land%rough_mom,  Time_atmos)
            call data_override ('LND', 'rough_heat', Land%rough_heat, Time_atmos)
            call data_override ('LND', 'albedo', Land%albedo,     Time_atmos)

            if (do_flux) then
              call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, &
                                       Atm, Land, Ice, Land_ice_atmos_boundary )
            end if
            call mpp_clock_end(fluxclock)
!            call print_memuse( 'sfc_boundary_layer' )

!      ---- atmosphere down ----

            call mpp_clock_begin(atmclock)
            call lnd_ice_atm_bnd_from_data ( Land_ice_atmos_boundary )
            if (do_atmos) &
              call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
            call mpp_clock_end(atmclock)
!            call print_memuse( 'update_atmos_model_down' )
            
            call mpp_clock_begin(fluxclockdn)
            call data_override ('ATM', 'flux_sw',  Atm%flux_sw, Time_atmos)
            call data_override ('ATM', 'flux_lw',  Atm%flux_lw, Time_atmos)
            call data_override ('ATM', 'lprec',    Atm%lprec,   Time_atmos)
            call data_override ('ATM', 'fprec',    Atm%fprec,   Time_atmos)
            call data_override ('ATM', 'coszen',   Atm%coszen,  Time_atmos)
            call data_override ('ATM', 'dtmass',   Atm%Surf_Diff%dtmass, Time_atmos)
            call data_override ('ATM', 'delta_t',  Atm%Surf_Diff%delta_t, Time_atmos)
            call data_override ('ATM', 'delta_q',  Atm%Surf_Diff%delta_q, Time_atmos)
            call data_override ('ATM', 'dflux_t',  Atm%Surf_Diff%dflux_t, Time_atmos)
            call data_override ('ATM', 'dflux_q',  Atm%Surf_Diff%dflux_q, Time_atmos)
            call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, &
                 Atmos_land_boundary, &
                 Atmos_ice_boundary )
            call mpp_clock_end(fluxclockdn)
!            call print_memuse( 'flux_down_from_atmos' )
            

!      --------------------------------------------------------------

!      ---- land model ----

            call mpp_clock_begin(landclock)
            call atm_lnd_bnd_from_data ( Atmos_land_boundary )
            if (do_land) &
              call update_land_model_fast( Atmos_land_boundary, Land )
            call mpp_clock_end(landclock)
!            call print_memuse( 'update_land_model_fast' )
            
!      ---- ice model ----


            call mpp_clock_begin(iceclock)
            call atm_ice_bnd_from_data ( Atmos_ice_boundary, Ice )
            if (do_ice) &
              call update_ice_model_fast( Atmos_ice_boundary, Ice )
            call mpp_clock_end(iceclock)
!            call print_memuse( 'update_ice_model_fast' )
            
!      --------------------------------------------------------------
!      ---- atmosphere up ----

            call mpp_clock_begin(fluxclockup)

            call data_override ( 'ICE', 't_surf', Ice%t_surf,  Time_atmos)
            call data_override ( 'LND', 't_ca',   Land%t_ca,   Time_atmos)
            call data_override ( 'LND', 't_surf', Land%t_surf, Time_atmos)
            call data_override ( 'LND', 'q_ca',   Land%q_ca,   Time_atmos)

            call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary )
            call mpp_clock_end(fluxclockup)
!            call print_memuse( 'flux_up_to_atmos' )
            
            call mpp_clock_begin(atmclock)
            if (do_atmos) &
              call update_atmos_model_up( Land_ice_atmos_boundary, Atm )
            call mpp_clock_end(atmclock)
!            call print_memuse( 'update_atmos_model_up' )
            
!--------------

         enddo

!   ------ end of atmospheric time step loop -----
         call mpp_clock_begin(landclock)
         if (do_land) call update_land_model_slow(Land)
         call mpp_clock_end(landclock)
!-----------------------------------------------------------------------

!
!     need flux call to put runoff and p_surf on ice grid
!
         call mpp_clock_begin(ilfluxclock)
         call flux_land_to_ice( Land, Ice, Land_ice_boundary )

         if(check_parallel) then
            if(Atm%pe) call mpp_set_current_pelist(atm_pelist)
            mesg = 'check Land_ice_boundary%runoff'
            call mpp_check_field(Land_ice_boundary%runoff, atm1_pelist, atm2_pelist, &
                                 Ice%domain,mesg)
            if( Atm%pe )  call mpp_set_current_pelist(Atm%pelist)
         endif

         Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?
         call mpp_clock_end(ilfluxclock)

!   ------ slow-ice model ------

         call mpp_clock_begin(iceclock)
         call atm_ice_bnd_from_data ( Atmos_ice_boundary, Ice )
         call lnd_ice_bnd_from_data ( Land_ice_boundary )
         if (do_ice) call update_ice_model_slow_dn( Atmos_ice_boundary, &
                                                    Land_ice_boundary, Ice )
         call mpp_clock_end(iceclock)
         Time = Time_atmos
     end if                     !Atm%pe block

     if( .NOT.concurrent .OR. nc.EQ.1 )then
         if(check_parallel) then
            if(group1) call mpp_set_current_pelist(pelist1)
            if(group2) call mpp_set_current_pelist(pelist2)
         else
            call mpp_set_current_pelist()
         endif
         call mpp_clock_begin(icefluxclock)
         call flux_ice_to_ocean( Ice, Ocean, Ice_ocean_boundary )
         call mpp_clock_end(icefluxclock)
     end if
!     call print_memuse( 'flux_ice_to_ocean' )
     
     if( Ocean%pe )then
         call mpp_set_current_pelist(Ocean%pelist)
         call ice_ocn_bnd_from_data ( Ice_ocean_boundary )
         do no = 1,num_ocean_calls

            write( text,'(a,2i4)' )'nc,no=', nc, no
            call print_memuse(text)

            Time_ocean = Time_ocean + Time_step_ocean

            ocean_seg_start = ( no .eq. 1 )               ! could eliminate these by
            ocean_seg_end   = ( no .eq. num_ocean_calls ) ! putting this loop in
                                                          ! update_ocean_model since
                                                          ! fluxes don't change here

            call mpp_clock_begin(ocnclock)
            if (do_ocean) call update_ocean_model( Ice_ocean_boundary, Ocean, &
                               ocean_seg_start, ocean_seg_end, num_ocean_calls)
            call mpp_clock_end(ocnclock)
!            call print_memuse( 'update_ocean_model' )
            
         enddo
!   ------ end of ocean time step loop -----
!-----------------------------------------------------------------------
         Time = Time_ocean
     end if
!--------------
  enddo

! Need final update of Ice_ocean_boundary for concurrent restart
  if(concurrent)then
     if(check_parallel) then
        if(group1) call mpp_set_current_pelist(pelist1)
        if(group2) call mpp_set_current_pelist(pelist2)
     else
        call mpp_set_current_pelist()
     endif

     call mpp_clock_begin(icefluxclock)
     call flux_ice_to_ocean( Ice, Ocean, Ice_ocean_boundary )
     call mpp_clock_end(icefluxclock)
  endif

  call mpp_set_current_pelist()
!-----------------------------------------------------------------------
  call check_system_clock ('END OF TIME LOOP')

  call mpp_clock_begin(endclock)
  call coupler_end

  call diag_manager_end (Time)
  call mpp_clock_end(endclock)
  call utilities_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
    integer :: unit, log_unit, ierr, io, id, jd, kd, m, i
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    logical :: use_namelist
    integer :: pe, npes
    integer :: atmos_pe_start=0, atmos_pe_end=0, &
               ocean_pe_start=0, ocean_pe_end=0, &
               ice_pe_start=0, ice_pe_end=0, &
               land_pe_start=0, land_pe_end=0
    integer :: atm1_pe_start=0, atm1_pe_end, &
               ocn1_pe_start=0, ocn1_pe_end, &
               atm2_pe_start=0, atm2_pe_end, &
               ocn2_pe_start=0, ocn2_pe_end
    integer :: diag_model_subset=DIAG_ALL
!-----------------------------------------------------------------------
!----- read namelist -------

    unit = open_file ('input.nml', action='read')
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call close_file (unit)

!----- write namelist to logfile (close log_unit later) -----

    log_unit = open_file ('logfile.out', action='append')
    if ( get_my_pe() == 0 ) then
        write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (log_unit, nml=coupler_nml)
    endif

!----- read date and calendar type from restart file -----

    if (file_exist('INPUT/coupler.res')) then
        unit = open_file ('INPUT/coupler.res',  &
             form='native', action='read')
        read  (unit) date
        read  (unit) calendar_type
        call close_file (unit)
        use_namelist = .false.
    else
        use_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---

    if ( use_namelist .or. override ) then

!----- override date with namelist values ------

        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date      = current_date
        endif

!----- override calendar type with namelist value -----

        if (calendar(1:6) == 'julian') then
            calendar_type = julian
        else if (calendar(1:6) == 'NOLEAP') then
            calendar_type = NOLEAP
        else if (calendar(1:10) == 'thirty_day') then
            calendar_type = thirty_day_months
        else if (calendar(1:11) == 'no_calendar') then
            calendar_type = no_calendar
        else if (calendar(1:1) /= ' ') then
            call error_mesg ('program coupler',  &
                 'invalid namelist value for calendar', FATAL)
        else
            call error_mesg ('program coupler',  &
                 'no namelist value for calendar', FATAL)
        endif

    endif

    call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----

    if ( get_my_pe() == 0 ) then
        write (log_unit,16) date(1),trim(month_name(date(2))),date(3:6)
    endif
    call close_file (log_unit)

16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!-----------------------------------------------------------------------
!------ initialize concurrent PEset management ------

!pe information
    pe = mpp_pe()
    npes = mpp_npes()
    if( ice_npes.NE.0 ) &
         call mpp_error( WARNING, 'coupler_init: pelists not yet implemented for ice.' )
    if( land_npes.NE.0 ) &
         call mpp_error( WARNING, 'coupler_init: pelists not yet implemented for land.' )

    if(check_parallel)then
!npes1+npes2 must equal npes
       if(npes1 .EQ. 0) npes1 = npes - npes2
       if(npes2 .EQ. 0) npes2 = npes - npes1
!both npes1 and npes2 should be nonzero

       if(npes .ne. npes1+npes2) call mpp_error(FATAL, 'coupler_init: npes1+npes2 must equal npes ')       
        if( npes1.EQ.0 .OR. npes2.EQ.0 ) &
             call mpp_error( FATAL, 'coupler_init: npes1 or npes2 must be specified for parallel checking' )       

       if( concurrent )then
       ! set up the npes for group1
          if(atm1_npes .EQ.0) atm1_npes = npes1 - ocn1_npes
          if(ocn1_npes .EQ.0) ocn1_npes = npes1 - atm1_npes
          if(npes1 .ne. atm1_npes+ocn1_npes) call mpp_error(FATAL,    &
                 'coupler_init: atm1_npes+ocn1_npes must equal npes1 ')       
          if( atm1_npes.EQ.0 .OR. ocn1_npes.EQ.0 ) &
             call mpp_error( FATAL, 'coupler_init: atm1_npes or ocn1_npes must be specified for parallel checking' )   
       ! set up the npes for group2
          if(atm2_npes .EQ.0) atm2_npes = npes2 - ocn2_npes
          if(ocn2_npes .EQ.0) ocn2_npes = npes2 - atm2_npes
          if(npes2 .ne. atm2_npes+ocn2_npes) call mpp_error(FATAL,    &
                 'coupler_init: atm2_npes+ocn2_npes must equal npes2 ')       
          if( atm2_npes.EQ.0 .OR. ocn2_npes.EQ.0 ) &
             call mpp_error( FATAL, 'coupler_init: atm2_npes or ocn2_npes must be specified for parallel checking' )    
          atm1_pe_start = 0                 ; atm1_pe_end = atm1_npes - 1
          ocn1_pe_start = atm1_npes          ; ocn1_pe_end = npes1 - 1
          atm2_pe_start = npes1             ; atm2_pe_end = npes1 + atm2_npes - 1
          ocn2_pe_start = npes1 + atm2_npes ; ocn2_pe_end = npes - 1
       else      ! don't provide disjoing choice when check parallel
          atm1_npes = npes1; ocn1_npes = npes1
          atm2_npes = npes2; ocn2_npes = npes2
          atm1_pe_start = 0     ; atm1_pe_end = npes1 - 1
          ocn1_pe_start = 0     ; ocn1_pe_end = npes1 - 1
          atm2_pe_start = npes1 ; atm2_pe_end = npes - 1
          ocn2_pe_start = npes1 ; ocn2_pe_end = npes - 1
       endif
       call mpp_set_pelist_init(npes1, npes2, atm1_pe_start, atm1_pe_end,  &
                                ocn1_pe_start, ocn1_pe_end, atm2_pe_start, &
                                atm2_pe_end, ocn2_pe_start, ocn2_pe_end )
       if(group1) then
          atmos_pe_start = atm1_pe_start; atmos_pe_end = atm1_pe_end
          ocean_pe_start = ocn1_pe_start; ocean_pe_end = ocn1_pe_end
          atmos_npes       = atm1_npes    ; ocean_npes     = ocn1_npes
       endif
       if(group2) then
          atmos_pe_start = atm2_pe_start; atmos_pe_end = atm2_pe_end
          ocean_pe_start = ocn2_pe_start; ocean_pe_end = ocn2_pe_end 
          atmos_npes       = atm2_npes    ; ocean_npes     = ocn2_npes                 
       endif
    else
       if( concurrent )then
!atmos_npes + ocean_npes must equal npes
          if( atmos_npes.EQ.0 )atmos_npes = npes - ocean_npes
          if( ocean_npes.EQ.0 )ocean_npes = npes - atmos_npes
!both must now be non-zero
          if( atmos_npes.EQ.0 .OR. ocean_npes.EQ.0 ) &
             call mpp_error( FATAL, 'coupler_init: atmos_npes or ocean_npes must be specified for concurrent coupling.' )
          if( atmos_npes+ocean_npes.NE.npes ) &
             call mpp_error( FATAL, 'coupler_init: atmos_npes+ocean_npes must equal npes for concurrent coupling.' )
          atmos_pe_start = 0
          atmos_pe_end = atmos_npes-1
          ocean_pe_start = atmos_npes
          ocean_pe_end = atmos_npes+ocean_npes-1
       else                        !serial timestepping
          if( atmos_npes.EQ.0 )atmos_npes = npes
          if( ocean_npes.EQ.0 )ocean_npes = npes
          if( max(atmos_npes,ocean_npes).EQ.npes )then !overlapping pelists
             atmos_pe_start = 0
             atmos_pe_end = atmos_npes-1
             ocean_pe_start = 0
             ocean_pe_end = ocean_npes-1
          else                    !disjoint pelists
             if( atmos_npes+ocean_npes.NE.npes ) call mpp_error( FATAL,  &
                 'coupler_init: atmos_npes+ocean_npes must equal npes for serial coupling on disjoint pelists.' )
             atmos_pe_start = 0
             atmos_pe_end = atmos_npes-1
             ocean_pe_start = atmos_npes
             ocean_pe_end = atmos_npes+ocean_npes-1
          end if
       end if
    end if
    allocate( Atm%pelist  (atmos_npes) )
    allocate( Ocean%pelist(ocean_npes) )
    Atm%pelist   = (/(i,i=atmos_pe_start,atmos_pe_end)/)
    Ocean%pelist = (/(i,i=ocean_pe_start,ocean_pe_end)/)
    Atm%pe = atmos_pe_start.LE.pe .AND. pe.LE.atmos_pe_end
    Ocean%pe = ocean_pe_start.LE.pe .AND. pe.LE.ocean_pe_end
    call mpp_declare_pelist( Atm%pelist , '_atm')
    call mpp_declare_pelist( Ocean%pelist , '_ocn')
    if( pe.EQ.mpp_root_pe() )write( stderr(),* ) &
         'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
          atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end
    if( pe.EQ.mpp_root_pe() )write( stderr(),* ) &
         'Concurrent coupling=', concurrent

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

!jwd Fork here is somewhat dangerous. It relies on "no side effects" from
!    diag_manager_init. diag_manager_init or this section should be 
!    re-architected to guarantee this or remove this assumption.
!    For instance, what follows assumes that get_base_date has the same
!    time for both Atm and Ocean pes. While this should be the case, the
!    possible error condition needs to be checked

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        if(atmos_npes /= npes)diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
    elseif( Ocean%pe )then  ! Error check above for disjoint pelists should catch any problem
        call mpp_set_current_pelist(Ocean%pelist)
        if(ocean_npes /= npes)diag_model_subset = DIAG_OCEAN  ! change diag_model_subset from DIAG_ALL
    end if
   call diag_manager_init(DIAG_MODEL_SUBSET=diag_model_subset)   ! initialize diag_manager for processor subset output
    call print_memuse( 'diag_manager_init' )
!-----------------------------------------------------------------------
!------ reset pelist to "full group" ------

    if(check_parallel) then                                               
        if(group1) call mpp_set_current_pelist(pelist1)
        if(group2) call mpp_set_current_pelist(pelist2)
    else                                                       
        call mpp_set_current_pelist()
    endif                              
!----- always override initial/base date with diag_manager value -----

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

!----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)
    Run_length = Time_end - Time

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    unit = open_file ('time_stamp.out', action='write')

    month = month_name(date(2))
    if ( get_my_pe() == 0 ) write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    if ( get_my_pe() == 0 ) write (unit,20) date, month(1:3)

    call close_file (unit)

20  format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_cpld  = set_time (dt_cpld ,0)
    Time_step_ocean = set_time (dt_ocean,0)
    Time_step_atmos = set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_ocean_calls = Time_step_cpld  / Time_step_ocean
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * num_ocean_calls * Time_step_ocean /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of ocean time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'atmos time step is not a multiple of the cpld time step', 2)

! ---- make sure cpld time step is a multiple of ocean time step ----

    if ( num_ocean_calls * Time_step_ocean /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the ocean time step', 2)

!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file
    
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos )
        call print_memuse( 'atmos_model_init' )

!---- land ----------
        call land_model_init( Land, Time_init, Time, Time_step_atmos, Time_step_cpld )
        call print_memuse( 'land_model_init' )

!---- ice -----------
        call ice_model_init( Ice, Time_init, Time, Time_step_atmos, Time_step_cpld )
        call print_memuse( 'ice_model_init' )
        call data_override_init(Atm_domain_in = Atm%domain, Ice_domain_in = Ice%domain, Land_domain_in=Land%domain)
    end if
    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
!---- ocean ---------
        call ocean_model_init( Ocean, Time_init, Time, Time_step_ocean )
        call print_memuse( 'ocean_model_init' )
        call data_override_init(Ocean_domain_in = Ocean%domain )
    end if
    if(check_parallel) then
        if(group1) call mpp_set_current_pelist(pelist1)
        if(group2) call mpp_set_current_pelist(pelist2)
    else
        call mpp_set_current_pelist()
    endif

    call mpp_broadcast_domain(Ice%domain)
    call mpp_broadcast_domain(Ocean%domain)
!-----------------------------------------------------------------------
!---- initialize flux exchange module ----
!    if( pe.EQ.mpp_root_pe()) &
!	write(*,*) 'Begin call data_override *&*&*&*&*&*&*&&**&*&*&*&'

    call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
         atmos_land_boundary, atmos_ice_boundary, land_ice_atmos_boundary, &
         land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

    Time_atmos = Time
    Time_ocean = Time

!---- initialize Ice for ice_ocean_boundary update ----
  ! Only ocean PEs read ice_ocean_boundary; data in ice_ocean_boundary
  ! struct is ignored in non-concurrent mode since atmos
  ! model has run one iteration before ocean model is called.
    if( Ocean%pe )then
      call mpp_set_current_pelist(Ocean%pelist)
      if (use_namelist .or. .NOT.concurrent) then
        call init_default_ice_ocean_boundary(ice_ocean_boundary)
      else
        call read_ice_ocean_boundary('INPUT/coupler_fluxes.res.nc', &
                                     ice_ocean_boundary,Ocean)
      endif
    endif

    if(check_parallel) then
        if(group1) call mpp_set_current_pelist(pelist1)
        if(group2) call mpp_set_current_pelist(pelist2)
    else
        call mpp_set_current_pelist()
    endif

!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----

    unit = open_file ('RESTART/coupler.res',  &
         form='native', action='write')
    call close_file (unit, status='delete')

!-----------------------------------------------------------------------

  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

    integer :: unit, date(6)
!-----------------------------------------------------------------------

    call mpp_set_current_pelist()

!----- compute current date ------

    call get_date (Time, date(1), date(2), date(3),  &
         date(4), date(5), date(6))

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------

    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
        if(concurrent)call write_ice_ocean_boundary('RESTART/coupler_fluxes.res.nc', &
                                      ice_ocean_boundary,Ocean)
        call ocean_model_end (Ocean)
    end if
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call atmos_model_end (Atm)
        call  land_model_end (Land)
        call   ice_model_end (Ice)
    end if
    call mpp_set_current_pelist()

!----- write restart file ------

    if ( get_my_pe() /= 0 ) return

    unit = open_file ('RESTART/coupler.res',  &
         form='native', action='write')
    write (unit) date
    write (unit) calendar_type
    call close_file (unit)

!-----------------------------------------------------------------------

  end subroutine coupler_end

!#######################################################################

  subroutine print_memuse( text, unit )
    use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_min, mpp_max, mpp_sum, stderr
    character(len=*), intent(in) :: text
    integer, intent(in), optional :: unit
    real :: m, mmin, mmax, mavg, mstd
    integer :: mu
    integer :: memuse            !external function

    if( .NOT.verbose )return
    mu = stderr(); if( PRESENT(unit) )mu = unit
    m = memuse()*1e-3
    mmin = m; call mpp_min(mmin)
    mmax = m; call mpp_max(mmax)
    mavg = m; call mpp_sum(mavg); mavg = mavg/mpp_npes()
    mstd = (m-mavg)**2; call mpp_sum(mstd); mstd = sqrt( mstd/mpp_npes() )
    if( debug_memuse ) then
    if( mpp_pe().EQ.mpp_root_pe() )write( mu,'(a32,4es11.3)' ) &
         'Memuse(MB) at '//trim(text)//'=', mmin, mmax, mstd, mavg
    endif
    return
  end subroutine print_memuse

!-----------------------------------------------------------------------
subroutine ocn_ice_bnd_from_data(x)
type (ocean_ice_boundary_type) :: x

  call data_override('ICE', 'u',         x%u,         Time)
  call data_override('ICE', 'v',         x%v,         Time)
  call data_override('ICE', 't',         x%t,         Time)
  call data_override('ICE', 's',         x%s,         Time)
  call data_override('ICE', 'frazil',    x%frazil,    Time)
  call data_override('ICE', 'sea_level', x%sea_level, Time)
end subroutine ocn_ice_bnd_from_data

!-----------------------------------------------------------------------
subroutine lnd_ice_atm_bnd_from_data(x)
type (land_ice_atmos_boundary_type) :: x

  call data_override('ATM', 't',         x%t,         Time)
  call data_override('ATM', 'albedo',    x%albedo,    Time)
  call data_override('ATM', 'land_frac', x%land_frac, Time)
  call data_override('ATM', 'dt_t',      x%dt_t,      Time)
  call data_override('ATM', 'dt_q',      x%dt_q,      Time)
  call data_override('ATM', 'u_flux',    x%u_flux,    Time)
  call data_override('ATM', 'v_flux',    x%v_flux,    Time)
  call data_override('ATM', 'dtaudv',    x%dtaudv,    Time)
  call data_override('ATM', 'u_star',    x%u_star,    Time)
  call data_override('ATM', 'b_star',    x%b_star,    Time)
! call data_override('ATM', 'q_star',    x%q_star,    Time)
  call data_override('ATM', 'rough_mom', x%rough_mom, Time)
end subroutine lnd_ice_atm_bnd_from_data

!-----------------------------------------------------------------------
subroutine atm_lnd_bnd_from_data(x)
type (atmos_land_boundary_type) :: x

  call data_override('LND', 't_flux',  x%t_flux, Time)
  call data_override('LND', 'q_flux',  x%q_flux, Time)
  call data_override('LND', 'lw_flux', x%lw_flux,Time)
  call data_override('LND', 'sw_flux', x%sw_flux,Time)
  call data_override('LND', 'lprec',   x%lprec,  Time)
  call data_override('LND', 'fprec',   x%fprec,  Time)
  call data_override('LND', 'dhdt',    x%dhdt,   Time)
  call data_override('LND', 'dedt',    x%dedt,   Time)
  call data_override('LND', 'dedq',    x%dedq,   Time)
  call data_override('LND', 'drdt',    x%drdt,   Time)
  call data_override('LND', 'drag_q',  x%drag_q, Time)
  call data_override('LND', 'p_surf',  x%p_surf, Time)
end subroutine atm_lnd_bnd_from_data

!-----------------------------------------------------------------------
subroutine atm_ice_bnd_from_data(x, Ice)
type (atmos_ice_boundary_type) :: x
type   (ice_data_type)         :: Ice
 logical :: ov
 ov = .false.

  call data_override('ICE', 'u_flux', x%u_flux,  Time)
  call data_override('ICE', 'v_flux', x%v_flux,  Time)
  call data_override('ICE', 't_flux', x%t_flux,  Time)
  call data_override('ICE', 'q_flux', x%q_flux,  Time)
  call data_override('ICE', 'lw_flux',x%lw_flux, Time)
  call data_override('ICE', 'lw_flux_dn',x%lw_flux, Time, override=ov)
  if (ov) then
    x%lw_flux = x%lw_flux - stefan*Ice%t_surf**4
  endif
  call data_override('ICE', 'sw_flux',x%sw_flux, Time)
  call data_override('ICE', 'sw_flux_dn',x%sw_flux, Time, override=ov)
  if (ov) then
    x%sw_flux = x%sw_flux*(1-Ice%albedo)
  endif
  call data_override('ICE', 'lprec',  x%lprec,   Time)
  call data_override('ICE', 'fprec',  x%fprec,   Time)
  call data_override('ICE', 'dhdt',   x%dhdt,    Time)
  call data_override('ICE', 'dedt',   x%dedt,    Time)
  call data_override('ICE', 'drdt',   x%drdt,    Time)
  call data_override('ICE', 'coszen', x%coszen,  Time)
  call data_override('ICE', 'p',      x%p,       Time)
end subroutine atm_ice_bnd_from_data

!-----------------------------------------------------------------------
subroutine lnd_ice_bnd_from_data(x)
type (land_ice_boundary_type) :: x

  call data_override('ICE', 'runoff' , x%runoff , Time)
  call data_override('ICE', 'calving', x%calving, Time)
end subroutine lnd_ice_bnd_from_data

!-----------------------------------------------------------------------
subroutine ice_ocn_bnd_from_data(x)
type (ice_ocean_boundary_type) :: x

  call data_override('OCN', 'u_flux',    x%u_flux   , Time)
  call data_override('OCN', 'v_flux',    x%v_flux   , Time)
  call data_override('OCN', 't_flux',    x%t_flux   , Time)
  call data_override('OCN', 'q_flux',    x%q_flux   , Time)
  call data_override('OCN', 'salt_flux', x%salt_flux, Time)
  call data_override('OCN', 'lw_flux',   x%lw_flux  , Time)
  call data_override('OCN', 'sw_flux',   x%sw_flux  , Time)
  call data_override('OCN', 'lprec',     x%lprec    , Time)
  call data_override('OCN', 'fprec',     x%fprec    , Time)
  call data_override('OCN', 'runoff',    x%runoff   , Time)
  call data_override('OCN', 'calving',   x%calving  , Time)
  call data_override('OCN', 'p',         x%p        , Time)
end subroutine ice_ocn_bnd_from_data


!##################################################################

end program coupler_main

