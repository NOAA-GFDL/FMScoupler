!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="V.Balaji@noaa.gov"> V. Balaji </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!  A main program that couples component models for atmosphere, ocean, land, 
!  and sea ice on independent grids. 
! </OVERVIEW>

! <DESCRIPTION>
!  This version couples model components representing atmosphere, ocean, land 
!  and sea ice on independent grids. Each model component is represented by a 
!  data type giving the instantaneous model state.
!
!  The component models are coupled to allow implicit vertical diffusion of 
!  heat and moisture at the interfaces of the atmosphere, land, and ice models. 
!  As a result, the atmosphere, land, and ice models all use the same time step. 
!  The atmospheric model has been separated into down and up calls that 
!  correspond to the down and up sweeps of the standard tridiagonal elimination.
!
!  The ocean interface uses explicit mixing. Fluxes to and from the ocean must
!  be passed through the ice model. This includes atmospheric fluxes as well as 
!  fluxes from the land to the ocean (runoff).
!
!  This program contains the model's main time loop. Each iteration of the 
!  main time loop is one coupled (slow) time step. Within this slow time step 
!  loop is a fast time step loop, using the atmospheric time step, where the
!  tridiagonal vertical diffusion equations are solved. Exchange between sea 
!  ice and ocean occurs once every slow timestep.
!
! <PRE>
!      MAIN PROGRAM EXAMPLE
!      --------------------
!
!         DO slow time steps (ocean)
!
!              call flux_ocean_to_ice
!
!              call ICE_SLOW_UP
!
!              DO fast time steps (atmos)
!
!                   call flux_calculation
!
!                   call ATMOS_DOWN
!
!                   call flux_down_from_atmos
!
!                   call LAND_FAST
!
!                   call ICE_FAST
!
!                   call flux_up_to_atmos
!
!                   call ATMOS_UP
!
!              END DO
!
!              call ICE_SLOW_DN
!
!              call flux_ice_to_ocean
!
!              call OCEAN
!
!         END DO

!  </PRE>

! </DESCRIPTION>
! <INFO>
!   <NOTE>
!     <PRE>
!   1.If no value is set for current_date, start_date, or calendar (or default value 
!     specified) then the value from restart file "INPUT/coupler.res" will be used. 
!     If neither a namelist value or restart file value exist the program will fail. 
!   2.The actual run length will be the sum of months, days, hours, minutes, and 
!     seconds. A run length of zero is not a valid option. 
!   3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>

!   <ERROR MSG="no namelist value for current_date " STATUS="FATAL">
!     A namelist value for current_date must be given if no restart file for
!     coupler_main (INPUT/coupler.res) is found. 
!   </ERROR>
!   <ERROR MSG="invalid namelist value for calendar" STATUS="FATAL">
!     The value of calendar must be 'julian', 'noleap', or 'thirty_day'. 
!     See the namelist documentation. 
!   </ERROR>
!   <ERROR MSG="no namelist value for calendar" STATUS="FATAL">
!     If no restart file is present, then a namelist value for calendar 
!     must be specified. 
!   </ERROR>
!   <ERROR MSG="initial time is greater than current time" STATUS="FATAL">
!     If a restart file is present, then the namelist value for either 
!     current_date or start_date was incorrectly set. 
!   </ERROR>
!   <ERROR MSG="run length must be multiple of ocean time step " STATUS="FATAL">
!     There must be an even number of ocean time steps for the requested run length. 
!   </ERROR>
!   <ERROR MSG="final time does not match expected ending time " STATUS="WARNING">
!     This error should probably not occur because of checks done at initialization time. 
!   </ERROR>

! </INFO>

  use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout
  use fms_mod,                 only: read_data, write_data

  use fms_io_mod,              only: fms_io_exit

  use diag_manager_mod,        only: diag_manager_init, diag_manager_end
  use diag_manager_mod,        only: DIAG_OCEAN, DIAG_OTHER, DIAG_ALL, get_base_date

  use field_manager_mod,       only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE

  use tracer_manager_mod,      only: tracer_manager_init, get_tracer_index
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, NO_TRACER

  use coupler_types_mod,       only: coupler_types_init

  use data_override_mod,       only: data_override_init

!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow

  use ice_model_mod,           only: ice_model_init, ice_model_end
  use ice_model_mod,           only: update_ice_model_slow_up
  use ice_model_mod,           only: update_ice_model_fast
  use ice_model_mod,           only: update_ice_model_slow_dn
  use ice_model_mod,           only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,           only: ocean_ice_boundary_type, atmos_ice_boundary_type

  use ocean_model_mod,         only: update_ocean_model, ocean_model_init
  use ocean_model_mod,         only: ocean_model_end, ocean_data_type, ice_ocean_boundary_type
  use ocean_model_mod,         only: read_ice_ocean_boundary, write_ice_ocean_boundary

!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod,       only: flux_exchange_init
  use flux_exchange_mod,       only: sfc_boundary_layer
  use flux_exchange_mod,       only: generate_sfc_xgrid
  use flux_exchange_mod,       only: flux_down_from_atmos
  use flux_exchange_mod,       only: flux_up_to_atmos
  use flux_exchange_mod,       only: flux_land_to_ice
  use flux_exchange_mod,       only: flux_ice_to_ocean
  use flux_exchange_mod,       only: flux_ocean_to_ice

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data

  use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_chksum
  use mpp_mod,                 only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, MAXPES
  use mpp_mod,                 only: stderr, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,                 only: mpp_set_current_pelist, mpp_declare_pelist

  use mpp_io_mod,              only: mpp_open, mpp_close
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use mpp_domains_mod,         only: mpp_broadcast_domain

  use memutils_mod,            only: print_memuse_stats

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 13.0.2.1 2006/06/06 01:19:47 ap Exp $'
  character(len=128) :: tag = '$Name: memphis_2006_08 $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  ! allow members of ocean type to be aliased (ap)
  type (ocean_data_type), target :: Ocean

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
!------ namelist interface -------

! <NAMELIST NAME="coupler_nml">
!   <DATA NAME="current_date"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. 
!   </DATA>
!   <DATA NAME="force_date_from_namelist"  TYPE="logical"  DEFAULT=".false.">
!     Flag that determines whether the namelist variable current_date should 
!     override the date in the restart file INPUT/coupler.res. If the restart 
!     file does not exist then force_date_from_namelist has not effect, the value of current_date 
!     will be used.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
!     'no_calendar' can not be used because the time_manager's date  function are used. 
!     All values must be lowercase.
!   </DATA>
!   <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
!     The number of months that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
!     The number of days that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
!     The number of hours that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
!     The number of minutes that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
!     The number of seconds that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="dt_atmos"  TYPE="integer"  DEFAULT="0">
!     Atmospheric model time step in seconds, including the fast coupling with 
!     land and sea ice. 
!   </DATA>
!   <DATA NAME="dt_ocean"  TYPE="integer"  DEFAULT="0">
!     Ocean model time step in seconds. 
!   </DATA>
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds for coupling between ocean and atmospheric models: 
!     must be an integral multiple of dt_atmos and dt_ocean. This is the "slow" timestep.
!   </DATA>
!  <DATA NAME="do_atmos, do_ocean, do_ice, do_land, do_flux" TYPE="logical">
!  If true (default), that particular model component (atmos, etc.) is run.
!  If false, the execution of that component is skipped. This is used when
!  ALL the output fields sent by that component to the coupler have been
!  overridden using the data_override feature. For advanced users only:
!  if you're not sure, you should leave these values at TRUE.
!  </DATA> 
!  <DATA NAME="concurrent" TYPE="logical">
!  If true, the ocean executes concurrently with the atmosphere-land-ocean
!   on a separate set of PEs.
!  If false (default), the execution is serial: call atmos... followed by
!  call ocean...
!  If using concurrent execution, you must set one of
!   atmos_npes and ocean_npes, see below.
!  </DATA> 
!  <DATA NAME="atmos_npes, ocean_npes" TYPE="integer">
!  If concurrent is set to true, we use these to set the list of PEs on which
!   each component runs.
!  At least one of them must be set to a number between 0 and NPES.
!  If exactly one of these two is set non-zero, the other is set to the
!   remainder from NPES.
!  If both are set non-zero they must add up to NPES.
!  </DATA> 
!  <DATA NAME="use_lag_fluxes" TYPE="logical">
!  If true, then mom4 is forced with SBCs from one coupling timestep ago
!  If false, then mom4 is forced with most recent SBCs.
!  For a leapfrog MOM coupling with dt_cpld=dt_ocean, lag fluxes
!  can be shown to be stable and current fluxes to be unconditionally unstable.
!  For dt_cpld>dt_ocean there is probably sufficient damping.
!  use_lag_fluxes is set to TRUE by default.
!  </DATA>
!  <DATA NAME="n_mask" TYPE="integer">
!    number of region to be masked out. Its value should be less than MAX_PES.
!  </DATA>
!  <DATA NAME="mask_list(2,MAXPES)" TYPE="integer, dimension(2,MAX_MASK_REGION)">
!    The position of the region to be masked out. mask_list(1,:) is the x-layout position
!    and mask_list(2,:) is y-layout position.  
!  </DATA>
!  <DATA NAME="layout_mask" TYPE="integer, dimension(2)">
!   Processor domain layout for all the component model. layout_mask need to be set when and only 
!   when n_mask is greater than 0 ( some domain region is masked out ). When this namelist is set,
!   it will overload the layout in each component model. The default value is (0,0).
!   Currently we require all the component model has the same layout and same grid size.
!  </DATA>
!   <NOTE>
!     <PRE>
!     1.If no value is set for current_date, start_date, or calendar (or default value specified) then the value from restart
!       file "INPUT/coupler.res" will be used. If neither a namelist value or restart file value exist the program will fail. 
!     2.The actual run length will be the sum of months, days, hours, minutes, and seconds. A run length of zero is not a
!       valid option. 
!     3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>
! </NAMELIST>


  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: force_date_from_namelist = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land
  integer :: dt_ocean = 0  ! ocean tracer timestep
  integer :: dt_cpld  = 0  ! fluxes passed between ice & ocean


  integer ::atmos_npes=0, ocean_npes=0, ice_npes=0, land_npes=0
  logical :: do_atmos =.true., do_land =.true., do_ice =.true., do_ocean=.true.
  logical :: do_flux =.true.
  logical :: concurrent=.FALSE.
  logical :: use_lag_fluxes=.TRUE.
  logical :: do_chksum=.FALSE.
  integer :: layout_mask(2) = (/0 , 0/)
  integer :: n_mask = 0
  integer :: mask_list(2, MAXPES), n, m 
  integer, parameter :: mp = 2*MAXPES
  data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, months, days, hours,      &
                         minutes, seconds, dt_cpld, dt_atmos, dt_ocean, do_atmos,    &
                         do_land, do_ice, do_ocean, do_flux, atmos_npes, ocean_npes, &
                         ice_npes, land_npes, concurrent, use_lag_fluxes, do_chksum, &
                         n_mask, layout_mask, mask_list

  integer :: initClock, mainClock, termClock

  character(len=80) :: text
  character(len=48), parameter                    :: mod_name = 'coupler_main_mod'
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_main'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!#######################################################################

  call mpp_init()
!these clocks are on the global pelist
  initClock = mpp_clock_id( 'Initialization' )
  mainClock = mpp_clock_id( 'Main loop' )
  termClock = mpp_clock_id( 'Termination' )
  call mpp_clock_begin(initClock)
  
  call fms_init
  call constants_init

  call coupler_init
  if(do_chksum) call coupler_chksum('coupler_init+', 0)

  call mpp_set_current_pelist()

  call mpp_clock_end (initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop

!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop ------

  do nc = 1, num_cpld_calls
     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         call generate_sfc_xgrid( Land, Ice )
     end if
     call mpp_set_current_pelist()

! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
! points when running concurrently. The calls as placed next to each other in
! concurrent mode to avoid multiple synchronizations within the main loop.
! This is only possible in the serial case when use_lag_fluxes.
     call flux_ocean_to_ice( Time, Ocean, Ice, Ocean_ice_boundary )

! Update Ice_ocean_boundary; first iteration is supplied by restart     
     if( use_lag_fluxes )then
         call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
     end if

     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         if (do_ice) call update_ice_model_slow_up( Ocean_ice_boundary, Ice )

!-----------------------------------------------------------------------
!   ------ atmos/fast-land/fast-ice integration loop -------


         do na = 1, num_atmos_calls

            Time_atmos = Time_atmos + Time_step_atmos

            if (do_atmos) then
              call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
            endif

            if (do_flux) then
              !if(do_chksum) call coupler_chksum('sfc-', (nc-1)*num_atmos_calls+na)
              call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, &
                                       Atm, Land, Ice, Land_ice_atmos_boundary )
              !if(do_chksum) call coupler_chksum('sfc+', (nc-1)*num_atmos_calls+na)
            end if

!      ---- atmosphere down ----

            if (do_atmos) &
              call update_atmos_model_down( Land_ice_atmos_boundary, Atm )

            call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, &
                 Atmos_land_boundary, &
                 Atmos_ice_boundary )
            

!      --------------------------------------------------------------

!      ---- land model ----

            if (do_land) &
              call update_land_model_fast( Atmos_land_boundary, Land )
            
!      ---- ice model ----
            if (do_ice) &
              call update_ice_model_fast( Atmos_ice_boundary, Ice )
            
!      --------------------------------------------------------------
!      ---- atmosphere up ----

            call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary )
            
            if (do_atmos) &
              call update_atmos_model_up( Land_ice_atmos_boundary, Atm )
            
!--------------

         enddo

!   ------ end of atmospheric time step loop -----
         if (do_land) call update_land_model_slow(Atmos_land_boundary,Land)
!-----------------------------------------------------------------------

!
!     need flux call to put runoff and p_surf on ice grid
!
         call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )

         Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?

!   ------ slow-ice model ------

         if (do_ice) call update_ice_model_slow_dn( Atmos_ice_boundary, &
                                                    Land_ice_boundary, Ice )
         Time = Time_atmos
     end if                     !Atm%pe block

     if( .NOT.use_lag_fluxes )then !this will serialize
         call mpp_set_current_pelist()
         call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
     end if
     
     if( Ocean%pe )then
         call mpp_set_current_pelist(Ocean%pelist)
         do no = 1,num_ocean_calls

!            if( mpp_pe().EQ.mpp_root_pe() )write( stderr(),'(a,2i4)' )'nc,no=', nc,no
            Time_ocean = Time_ocean + Time_step_ocean

            ocean_seg_start = ( no .eq. 1 )               ! could eliminate these by
            ocean_seg_end   = ( no .eq. num_ocean_calls ) ! putting this loop in
                                                          ! update_ocean_model since
                                                          ! fluxes don't change here

            if (do_ocean) call update_ocean_model( Ice_ocean_boundary, Ocean, &
                               ocean_seg_start, ocean_seg_end, num_ocean_calls)
            
         enddo
!   ------ end of ocean time step loop -----
!-----------------------------------------------------------------------
         Time = Time_ocean
     end if
!--------------
     if(do_chksum) call coupler_chksum('MAIN_LOOP+', nc)
     write( text,'(a,i4)' )'Main loop at coupling timestep=', nc
     call print_memuse_stats(text)

  enddo

! Need final update of Ice_ocean_boundary for concurrent restart
!  if( concurrent )then
!      call mpp_set_current_pelist()
!      call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
!  endif

  call mpp_set_current_pelist()
!-----------------------------------------------------------------------
  call mpp_clock_end(mainClock)
  call mpp_clock_begin(termClock)

  if(do_chksum) call coupler_chksum('coupler_end-', nc)
  call coupler_end

  call diag_manager_end (Time)
  call mpp_clock_end(termClock)

  call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    character(len=64), parameter    :: sub_name = 'coupler_init'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: warn_header =                                &
         '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: unit,  ierr, io,    m, i
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes
    integer :: atmos_pe_start=0, atmos_pe_end=0, &
               ocean_pe_start=0, ocean_pe_end=0
    integer :: n
    integer :: diag_model_subset=DIAG_ALL
    logical :: other_fields_exist
    logical, allocatable :: maskmap(:,:)
!-----------------------------------------------------------------------
!----- read namelist -------

    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call mpp_close(unit)

!----- write namelist to logfile -----
    call write_version_number (version, tag)
    if( mpp_pe() == mpp_root_pe() )write( stdlog(), nml=coupler_nml )

!----- read date and calendar type from restart file -----

    if( file_exist('INPUT/coupler.res') )then
!Balaji: currently written in binary, needs form=MPP_NATIVE
        call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY )
        read( unit,*,err=999 )calendar_type
        read( unit,* )date_init
        read( unit,* )date
        goto 998 !back to fortran-4
!read old-style coupler.res
999     call mpp_close(unit)
        call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY, form=MPP_NATIVE )
        read(unit)calendar_type
        read(unit)date
998     call mpp_close(unit)
    else
        force_date_from_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---

    if ( force_date_from_namelist ) then

        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date      = current_date
        endif

!----- override calendar type with namelist value -----

        select case( uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        case default
            call mpp_error( FATAL, 'COUPLER_MAIN: coupler_nml entry calendar must be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
        end select
!        if ( uppercase(calendar(1:6)) == 'JULIAN') then
!            calendar_type = JULIAN
!        else if ( uppercase(calendar(1:6)) == 'NOLEAP') then
!            calendar_type = NOLEAP
!        else if ( uppercase(calendar(1:10)) == 'THIRTY_DAY') then
!            calendar_type = THIRTY_DAY_MONTHS
!        else if ( uppercase(calendar(1:11)) == 'NO_CALENDAR') then
!            calendar_type = NO_CALENDAR
!        else if (calendar(1:1) /= ' ') then
!            call error_mesg ('program coupler',  &
!                 'invalid namelist value for calendar', FATAL)
!        else
!            call error_mesg ('program coupler',  &
!                 'no namelist value for calendar', FATAL)
!        endif

    endif

    call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe().EQ.mpp_root_pe() ) &
         write( stdlog(), 16 )date(1),trim(month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

    npes = mpp_npes()
!----- check the value of layout and setup the maskmap for domain layout.
    if( n_mask > 0 ) then
       if(do_atmos .OR. do_land) call mpp_error(FATAL, &
            'program coupler: do_atmos and do_land should be false when n_mask > 0')

       if(concurrent) call mpp_error(FATAL, &
            'program coupler: can not run concurrent run when some regions are masked out')
       if( layout_mask(1)*layout_mask(2) - n_mask .NE. npes ) call mpp_error(FATAL, &
            'program coupler: layout(1)*layout(2) - n_mask should equal to npes when n_mask>0')
       call mpp_error(NOTE, 'program coupler: layout_mask and mask_list is set in coupler_nml, ' // &
                            'the value of layout_mask will override the layout specified in each component model')

       allocate(maskmap(layout_mask(1), layout_mask(2)) )
       maskmap = .TRUE.
       do n=1, n_mask
          if (mask_list(1,n) .gt. layout_mask(1) ) &
             call mpp_error( FATAL, 'program coupler: mask_list elements outside layout defines.' )
          if (mask_list(2,n) .gt. layout_mask(2) ) &
             call mpp_error( FATAL, 'program coupler: mask_list elements outside layout defines.' )
          maskmap(mask_list(1,n),mask_list(2,n)) = .false.
       enddo
       !--- copy maskmap value to each model data type
       allocate(Atm%maskmap(layout_mask(1), layout_mask(2)), Land%maskmap(layout_mask(1), layout_mask(2)) )
       allocate(Ice%maskmap(layout_mask(1), layout_mask(2)), Ocean%maskmap(layout_mask(1), layout_mask(2)))
       Atm%maskmap = maskmap;  Land%maskmap = maskmap
       Ice%maskmap = maskmap;  Ocean%maskmap = maskmap
       deallocate(maskmap)
    else
       if( layout_mask(1)*layout_mask(2) .NE. 0 ) call mpp_error(NOTE, &
            'program coupler: when no region is masked out, layout_mask need not be set' )
    end if


!-----------------------------------------------------------------------
!------ initialize concurrent PEset management ------

!pe information
    pe = mpp_pe()
    if( ice_npes.NE.0 ) &
         call mpp_error( WARNING, 'coupler_init: pelists not yet implemented for ice.' )
    if( land_npes.NE.0 ) &
         call mpp_error( WARNING, 'coupler_init: pelists not yet implemented for land.' )

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
        if( .NOT.use_lag_fluxes )call mpp_error( WARNING, &
             'coupler_init: you have set concurrent=TRUE and use_lag_fluxes=FALSE &
            & in coupler_nml. When not using lag fluxes, components &
            & will synchronize at two points, and thus run serially.' )
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
!messages
    write( text,'(a,2i6)' )'Atmos PE range: ', atmos_pe_start, atmos_pe_end
    call mpp_error( NOTE, 'coupler_init: '//trim(text) )
    write( text,'(a,2i6)' )'Ocean PE range: ', ocean_pe_start, ocean_pe_end
    call mpp_error( NOTE, 'coupler_init: '//trim(text) )
    if( concurrent )then
        call mpp_error( NOTE, 'coupler_init: Running with CONCURRENT coupling.' )
    else
        call mpp_error( NOTE, 'coupler_init: Running with SERIAL coupling.' )
    end if
    if( use_lag_fluxes )then
        call mpp_error( NOTE, 'coupler_init: Sending LAG fluxes to ocean.' )
    else
        call mpp_error( NOTE, 'coupler_init: Sending most recent fluxes to ocean.' )
    end if
    allocate( Atm%pelist  (atmos_npes) )
    allocate( Ocean%pelist(ocean_npes) )
    Atm%pelist   = (/(i,i=atmos_pe_start,atmos_pe_end)/)
    Ocean%pelist = (/(i,i=ocean_pe_start,ocean_pe_end)/)
    Atm%pe = atmos_pe_start.LE.pe .AND. pe.LE.atmos_pe_end
    Ocean%pe = ocean_pe_start.LE.pe .AND. pe.LE.ocean_pe_end
    call mpp_declare_pelist( Atm%pelist,   '_atm' )
    call mpp_declare_pelist( Ocean%pelist, '_ocn' )
    if( concurrent .AND. pe.EQ.mpp_root_pe() )then
        write( stdlog(),'(a)' )'Using concurrent coupling...'
        write( stdlog(),'(a,4i4)' ) &
         'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
          atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end
    end if

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
    call print_memuse_stats( 'diag_manager_init' )
!-----------------------------------------------------------------------
!------ reset pelist to "full group" ------

    call mpp_set_current_pelist()
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

    call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )

    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe().EQ.mpp_root_pe() ) write (unit,20) date, month(1:3)

    call mpp_close(unit)

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
         'atmos time step is not a multiple of the cpld time step', FATAL)

! ---- make sure cpld time step is a multiple of ocean time step ----

    if ( num_ocean_calls * Time_step_ocean /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the ocean time step', FATAL)

!
!       Initialize the tracer manager. This needs to be done on all PEs,
!       before the individual models are initialized.
!

    call tracer_manager_init
!
!       Initialize the coupler types
!

    call coupler_types_init

!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file
    
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos )
        call print_memuse_stats( 'atmos_model_init' )

!---- land ----------
        call land_model_init( Atmos_land_boundary, Land, Time_init, Time, &
             Time_step_atmos, Time_step_cpld )
        call print_memuse_stats( 'land_model_init' )

!---- ice -----------
        call ice_model_init( Ice, Time_init, Time, Time_step_atmos, Time_step_cpld )
        call print_memuse_stats( 'ice_model_init' )
        call data_override_init(Atm_domain_in = Atm%domain, Ice_domain_in = Ice%domain, Land_domain_in=Land%domain)
    end if
    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
!---- ocean ---------
        call ocean_model_init( Ocean, Time_init, Time, Time_step_ocean )
        call print_memuse_stats( 'ocean_model_init' )
        call data_override_init(Ocean_domain_in = Ocean%domain )
    end if
    call mpp_set_current_pelist()

    call mpp_broadcast_domain(Ice%domain)
    call mpp_broadcast_domain(Ocean%domain)
!-----------------------------------------------------------------------
!---- initialize flux exchange module ----
    call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
         atmos_ice_boundary, land_ice_atmos_boundary, &
         land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

    Time_atmos = Time
    Time_ocean = Time

!
!       read in extra fields for the air-sea gas fluxes
!

    if ( Atm%pe ) then
      call mpp_set_current_pelist(Atm%pelist)
      do n = 1, Ice%ocean_fluxes%num_bcs  !{
        other_fields_exist = .false.
        do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
          if (field_exist(Ice%ocean_fluxes%bc(n)%ice_file_in,                           &
                          Ice%ocean_fluxes%bc(n)%field(m)%name)) then
            other_fields_exist = .true.
            write (stdout(),*) trim(note_header), ' Reading restart info for ',         &
                 trim(Ice%ocean_fluxes%bc(n)%field(m)%name), ' from ',                  &
                 trim(Ice%ocean_fluxes%bc(n)%ice_file_in)
            call read_data(Ice%ocean_fluxes%bc(n)%ice_file_in,                          &
                 Ice%ocean_fluxes%bc(n)%field(m)%name,                                  &
                 Ice%ocean_fluxes%bc(n)%field(m)%values, Ice%domain)
          elseif (other_fields_exist) then
            call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                 trim(Ice%ocean_fluxes%bc(n)%field(m)%name) // ' in file ' //           &
                 trim(Ice%ocean_fluxes%bc(n)%ice_file_in))
          endif
        enddo  !} m
      enddo  !} n
    endif
    if ( Ocean%pe ) then
      call mpp_set_current_pelist(Ocean%pelist)
      do n = 1, Ocean%fields%num_bcs  !{
        other_fields_exist = .false.
        do m = 1, Ocean%fields%bc(n)%num_fields  !{
          if (field_exist(Ocean%fields%bc(n)%ocean_file_in,                &
                          Ocean%fields%bc(n)%field(m)%name)) then
            other_fields_exist = .true.
            write (stdout(),*) trim(note_header), ' Reading restart info for ',         &
                 trim(Ocean%fields%bc(n)%field(m)%name), ' from ',         &
                 trim(Ocean%fields%bc(n)%ocean_file_in)
            call read_data(Ocean%fields%bc(n)%ocean_file_in,               &
                 Ocean%fields%bc(n)%field(m)%name,                         &
                 Ocean%fields%bc(n)%field(m)%values, Ocean%domain)
          elseif (other_fields_exist) then
            call mpp_error(FATAL, trim(error_header) // ' Couldn''t find field ' //     &
                 trim(Ocean%fields%bc(n)%field(m)%name) // ' in file ' //  &
                 trim(Ocean%fields%bc(n)%ocean_file_in))
          endif
        enddo  !} m
      enddo  !} n
    endif

    call mpp_set_current_pelist()

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    call mpp_open( unit, 'RESTART/file' )
    call mpp_close(unit, MPP_DELETE)

!-----------------------------------------------------------------------
    call print_memuse_stats('coupler_init')
  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

    integer :: unit, date(6)
    integer :: m
    integer :: n
!-----------------------------------------------------------------------

    call mpp_set_current_pelist()

!----- compute current date ------

    call get_date (Time, date(1), date(2), date(3),  &
         date(4), date(5), date(6))

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------
!the call to fms_io_exit has been moved here
!this will work for serial code or concurrent (disjoint pelists)
!but will fail on overlapping but unequal pelists
    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
!        call write_ice_ocean_boundary('RESTART/coupler_fluxes.res.nc', &
!                                      ice_ocean_boundary,Ocean)
        call ocean_model_end (Ocean)
        do n = 1, Ocean%fields%num_bcs  !{
          do m = 1, Ocean%fields%bc(n)%num_fields  !{
            call write_data(Ocean%fields%bc(n)%ocean_file_out,     &
                 Ocean%fields%bc(n)%field(m)%name,                 &
                 Ocean%fields%bc(n)%field(m)%values, Ocean%domain)
          enddo  !} m
        enddo  !} n
    end if
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call atmos_model_end (Atm)
        call  land_model_end (Atmos_land_boundary, Land)
        call   ice_model_end (Ice)
        do n = 1, Ice%ocean_fluxes%num_bcs  !{
          do m = 1, Ice%ocean_fluxes%bc(n)%num_fields  !{
            call write_data(Ice%ocean_fluxes%bc(n)%ice_file_out,                &
                 Ice%ocean_fluxes%bc(n)%field(m)%name,                          &
                 Ice%ocean_fluxes%bc(n)%field(m)%values, Ice%domain)
          enddo  !} m
        enddo  !} n
    end if
    call fms_io_exit
    call mpp_set_current_pelist()

!----- write restart file ------

    call mpp_open( unit, 'RESTART/coupler.res', nohdrs=.TRUE. )
    if ( mpp_pe().EQ.mpp_root_pe() )then
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write( unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write( unit, '(6i6,8x,a)' )date, &
             'Current model time: year, month, day, hour, minute, second'
    end if
    call mpp_close(unit)

!-----------------------------------------------------------------------

  end subroutine coupler_end

  subroutine coupler_chksum(id, timestep)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep

    type :: tracer_ind_type
       integer :: atm, ice, lnd ! indices of the tracer in the respective models
    end type tracer_ind_type
    integer                            :: n_atm_tr, n_lnd_tr, n_exch_tr
    integer                            :: n_atm_tr_tot, n_lnd_tr_tot
    integer                            :: i, tr, n, m
    type(tracer_ind_type), allocatable :: tr_table(:)
    character(32) :: tr_name

    call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
                             num_prog=n_atm_tr)
    call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
                             num_prog=n_lnd_tr)

    ! assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
       call get_tracer_names( MODEL_ATMOS, i, tr_name )
       tr_table(n)%atm = i
       tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
       tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
       if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) &
            n = n+1
    enddo
    n_exch_tr = n-1

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)

    write(stdout(),*) 'BEGIN CHECKSUM(Atm):: ', id, timestep
    write(stdout(),100) 'atm%t_bot', mpp_chksum(atm%t_bot)
    write(stdout(),100) 'atm%z_bot', mpp_chksum(atm%z_bot)
    write(stdout(),100) 'atm%p_bot', mpp_chksum(atm%p_bot)
    write(stdout(),100) 'atm%u_bot', mpp_chksum(atm%u_bot)
    write(stdout(),100) 'atm%v_bot', mpp_chksum(atm%v_bot)
    write(stdout(),100) 'atm%p_surf', mpp_chksum(atm%p_surf)
    write(stdout(),100) 'atm%gust', mpp_chksum(atm%gust)
    do tr = 1,n_exch_tr
       n = tr_table(tr)%atm
       if(n /= NO_TRACER ) then
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
          write(stdout(),100) 'atm%'//trim(tr_name), mpp_chksum(Atm%tr_bot(:,:,n))
       endif
    enddo

    write(stdout(),100) 'land%t_surf', mpp_chksum(land%t_surf)
    write(stdout(),100) 'land%t_ca', mpp_chksum(land%t_ca)
    write(stdout(),100) 'land%rough_mom', mpp_chksum(land%rough_mom)
    write(stdout(),100) 'land%rough_heat', mpp_chksum(land%rough_heat)
    write(stdout(),100) 'land%rough_scale', mpp_chksum(land%rough_scale)
    do tr = 1,n_exch_tr
       n = tr_table(tr)%lnd
       if(n /= NO_TRACER ) then
          call get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
          write(stdout(),100) 'land%'//trim(tr_name), mpp_chksum(Land%tr(:,:,:,n))
       endif
    enddo

    write(stdout(),100) 'ice%t_surf', mpp_chksum(ice%t_surf)
    write(stdout(),100) 'ice%rough_mom', mpp_chksum(ice%rough_mom)
    write(stdout(),100) 'ice%rough_heat', mpp_chksum(ice%rough_heat)
    write(stdout(),100) 'ice%rough_moist', mpp_chksum(ice%rough_moist)
    write(stdout(),*) 'STOP CHECKSUM(Atm):: ', id, timestep

    !endif

    !if( Ocean%pe )then
        !call mpp_set_current_pelist(Ocean%pelist)

    write(stdout(),*) 'BEGIN CHECKSUM(Ice):: ', id, timestep
    do n = 1, ice%ocean_fields%num_bcs  !{
       do m = 1, ice%ocean_fields%bc(n)%num_fields  !{
          !write(stdout(),101) 'ice%', m, n, mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
          write(stdout(),101) 'ice%',trim(ice%ocean_fields%bc(n)%name), &
               trim(ice%ocean_fields%bc(n)%field(m)%name), mpp_chksum(Ice%ocean_fields%bc(n)%field(m)%values)
       enddo  !} m
    enddo  !} n
    write(stdout(),*) 'STOP CHECKSUM(Ice):: ', id, timestep

    endif

    deallocate(tr_table)

    call mpp_set_current_pelist()

  end subroutine coupler_chksum

  !#######################################################################

  end program coupler_main

