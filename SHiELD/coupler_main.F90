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

program coupler_main

!-----------------------------------------------------------------------
!
!   program that couples component models for the atmosphere,
!   ocean (amip), land, and sea-ice using the exchange module
!
!-----------------------------------------------------------------------

use FMS
use FMSconstants,    only: fmsconstants_init
use atmos_model_mod, only: atmos_model_init, atmos_model_end,  &
                           update_atmos_model_dynamics,        &
                           update_atmos_radiation_physics,     &
                           update_atmos_model_state,           &
                           atmos_data_type, atmos_model_restart
!--- FMS old io
use fms_io_mod, only: fms_io_exit!< This can't be removed until fms_io is not used at all

implicit none

!-----------------------------------------------------------------------
  character(len=128) :: version = 'unknown'
  character(len=128) :: tag = 'FMSCoupler_SHiELD'

!-----------------------------------------------------------------------
!---- model defined-types ----

 type (atmos_data_type) :: Atm

!-----------------------------------------------------------------------
! ----- coupled model time -----

   type (time_type) :: Time_atmos, Time_init, Time_end,  &
                       Time_step_atmos, Time_step_ocean, &
                       Time_restart, Time_step_restart,  &
                       Time_start_restart, Time_restart_aux, &
                       Time_step_restart_aux, Time_start_restart_aux, &
                       Time_duration_restart_aux, Time_restart_end_aux

   integer :: num_cpld_calls, num_atmos_calls, nc, na, ret

! ----- coupled model initial date -----

   integer :: date_init(6)
   integer :: calendar_type = -99

! ----- timing flags -----

   integer :: initClock, mainClock, termClock
   integer, parameter :: timing_level = 1

! ----- namelist -----
   integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /) !< The date that the current integration starts with
   character(len=17) :: calendar = '                 '  !< The calendar type used by the current integration.  Valid values are
                                                        !! consistent with the time_manager module: 'gregorian', 'julian',
                                                        !! 'noleap', or 'thirty_day'. The value 'no_calendar' cannot be used
                                                        !! because the time_manager's date !! functions are used.
                                                        !! All values must be lower case.
   logical :: force_date_from_namelist = .false.  !< Flag that determines whether the namelist variable current_date should override
                                                  !! the date in the restart file `INPUT/coupler.res`.  If the restart file does not
                                                  !! exist then force_date_from_namelist has no effect, the value of current_date
                                                  !! will be used.
   integer :: years=0    !< Number of years the current integration will be run
   integer :: months=0   !< Number of months the current integration will be run
   integer :: days=0     !< Number of days the current integration will be run
   integer :: hours=0    !< Number of hours the current integration will be run
   integer :: minutes=0  !< Number of minutes the current integration will be run
   integer :: seconds=0  !< Number of seconds the current integration will be run
   integer :: dt_atmos = 0  !< Atmospheric model time step in seconds
   integer :: dt_ocean = 0  !< Ocean model time step in seconds - NOT USED IN THIS MODEL
   integer :: restart_days = 0  !< Time interval in days to write out intermediate restart files
   integer :: restart_secs = 0  !< Time interval in seconds to write out intermediate restart files
   integer :: restart_start_days = 0  !< Start time in days to write out intermediate restart files
   integer :: restart_start_secs = 0  !< Start time in seconds to write out intermediate restart files
   integer :: restart_days_aux = 0  !< Time interval in days for auxiliary restart files
   integer :: restart_secs_aux = 0  !< Time interval in seconds for auxiliary restart files
   integer :: restart_start_days_aux = 0  !< Start time in days for auxiliary restart files
   integer :: restart_start_secs_aux = 0  !< Start time in days for auxiliary restart files
   integer :: restart_duration_days_aux = 0  !< Duration in days for auxiliary restart files
   integer :: restart_duration_secs_aux = 0  !< Duration in seconds for auxiliary restart files
   integer :: atmos_nthreads = 1  !< Number of OpenMP threads to use in the atmosphere
   logical :: use_hyper_thread = .false.  !< If .TRUE>, affinity placement (if activated) will consider virtual cores
                                          !! in the placement algorithm
   integer :: iau_offset = 0


   namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, &
                          years, months, days, hours, minutes, seconds, &
                          iau_offset, dt_atmos, dt_ocean, atmos_nthreads, &
                          use_hyper_thread, restart_secs, restart_days, &
                          restart_start_secs, restart_start_days, &
                          restart_secs_aux, restart_days_aux, &
                          restart_start_secs_aux, restart_start_days_aux, &
                          restart_duration_secs_aux, restart_duration_days_aux

! ----- local variables -----
   character(len=32) :: timestamp
   logical :: intrm_rst, intrm_rst_1step

!#######################################################################

 call fms_init()
 call sat_vapor_pres_init()
 call fmsconstants_init()

 initClock = mpp_clock_id( 'Initialization' )
 call mpp_clock_begin (initClock) !nesting problem

 call coupler_init
 call print_memuse_stats('after coupler init')

 call mpp_set_current_pelist()
 call mpp_clock_end (initClock) !end initialization
 mainClock = mpp_clock_id( 'Main loop' )
 termClock = mpp_clock_id( 'Termination' )
 call mpp_clock_begin(mainClock) !begin main loop

 do nc = 1, num_cpld_calls

    Time_atmos = Time_atmos + Time_step_atmos

    call update_atmos_model_dynamics (Atm)

    call update_atmos_radiation_physics (Atm)

    call update_atmos_model_state (Atm)

!--- intermediate restart
    if (intrm_rst) then
      if (nc /= num_cpld_calls) then
        if (intrm_rst_1step .and. nc == 1) then
          timestamp = date_to_string (Time_atmos)
          call atmos_model_restart(Atm, timestamp)
          call coupler_restart(timestamp)
        endif
        if (Time_atmos == Time_restart .or. Time_atmos == Time_restart_aux) then
          if (Time_atmos == Time_restart) then
            timestamp = date_to_string (Time_restart)
          else
            timestamp = date_to_string (Time_restart_aux)
          endif
          call atmos_model_restart(Atm, timestamp)
          call coupler_restart(timestamp)
          if (Time_atmos == Time_restart) &
              Time_restart = Time_restart + Time_step_restart
          if ((restart_secs_aux > 0 .or. restart_days_aux > 0) .and. &
              Time_atmos == Time_restart_aux .and. &
              Time_restart_aux < Time_restart_end_aux) then
            Time_restart_aux = Time_restart_aux + Time_step_restart_aux
          endif
        endif
      endif
    endif

    call print_memuse_stats('after full step')

 enddo

!-----------------------------------------------------------------------

 call mpp_set_current_pelist()
 call mpp_clock_end(mainClock)
 call mpp_clock_begin(termClock)

 call coupler_end
 call mpp_set_current_pelist()
 call mpp_clock_end(termClock)

 call fms_end

!-----------------------------------------------------------------------

contains

!#######################################################################

   subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
    integer :: total_days, total_seconds, unit, ierr, io
    integer :: n
    integer :: date(6), flags
    type (time_type) :: Run_length
    character(len=9) :: month

    character(len=:), dimension(:), allocatable :: restart_file !< Restart file saved as a string
    integer :: time_stamp_unit !< Unif of the time_stamp file
    integer :: ascii_unit  !< Unit of a dummy ascii file

!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

!----- read namelist -------
!----- for backwards compatibilty read from file coupler.nml -----
   read(input_nml_file, nml=coupler_nml, iostat=io)
   ierr = check_nml_error(io, 'coupler_nml')

!----- write namelist to logfile -----
   call write_version_number (version, tag)
   if (mpp_pe() == mpp_root_pe()) write(stdlog(),nml=coupler_nml)

!----- allocate and set the pelist (to the global pelist) -----
   allocate( Atm%pelist  (mpp_npes()) )
   call mpp_get_current_pelist(Atm%pelist)

!----- read restart file -----
    if (file_exists('INPUT/coupler.res')) then
       call ascii_read('INPUT/coupler.res', restart_file)
       read(restart_file(1), *) calendar_type
       read(restart_file(2), *) date_init
       read(restart_file(3), *) date
       deallocate(restart_file)
    else
       force_date_from_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---
    if ( force_date_from_namelist ) then
      if ( sum(current_date) <= 0 ) then
         call error_mesg ('program coupler',  &
              'no namelist value for current_date', FATAL)
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
          call mpp_error ( FATAL, 'COUPLER_MAIN: coupler_nml entry calendar must '// &
                                  'be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
      end select

    endif

    call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----
    if ( mpp_pe() == mpp_root_pe() ) then
      write (stdlog(),16) date(1),trim(month_name(date(2))),date(3:6)
    endif
 16 format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')

!------ setting affinity ------
!$  call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads)
!$  call omp_set_num_threads(atmos_nthreads)

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------
    call diag_manager_init (TIME_INIT=date)

!----- always override initial/base date with diag_manager value -----
    call get_base_date ( date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------
    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------
    Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                           date_init(4), date_init(5), date_init(6))

    Time_atmos = set_date (date(1), date(2), date(3),  &
                           date(4), date(5), date(6))

!-----------------------------------------------------------------------
!----- compute the ending time (compute days in each month first) -----
!
!   (NOTE: if run length in months then starting day must be <= 28)
    if ( months > 0 .and. date(3) > 28 )     &
        call error_mesg ('program coupler',  &
       'if run length in months then starting day must be <= 28', FATAL)

    Time_end = Time_atmos
    total_days = 0
    do n = 1, months
       total_days = total_days + days_in_month(Time_end)
       Time_end = Time_atmos + set_time (0,total_days)
    enddo

    total_days    = total_days + days
    total_seconds = hours*3600 + minutes*60 + seconds
    Run_length    = set_time (total_seconds,total_days)
    Time_end      = Time_atmos + Run_length

    !Need to pass Time_end into diag_manager for multiple thread case.
    call diag_manager_set_time_end(Time_end)


!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------
    if ( mpp_pe().EQ.mpp_root_pe() ) open(newunit = time_stamp_unit, file='time_stamp.out', status='replace', form='formatted')

    month = month_name(date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
                             date(4), date(5), date(6))
    month = month_name(date(2))
    if ( mpp_pe() == mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)

    if ( mpp_pe().EQ.mpp_root_pe() ) close (time_stamp_unit)

 20 format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------
    Time_step_atmos = set_time (dt_atmos,0)
    Time_step_ocean = set_time (dt_ocean,0)
    num_cpld_calls  = Run_length / Time_step_ocean
    num_atmos_calls = Time_step_ocean / Time_step_atmos

    Time_step_restart = set_time (restart_secs, restart_days)
    if (restart_start_secs > 0 .or. restart_start_days > 0) then
       Time_start_restart = set_time (restart_start_secs, restart_start_days)
       Time_restart = Time_atmos + Time_start_restart
    else
       Time_restart = Time_atmos + Time_step_restart
    end if
    Time_step_restart_aux = set_time (restart_secs_aux, restart_days_aux)
    Time_duration_restart_aux = set_time (restart_duration_secs_aux, restart_duration_days_aux)
    Time_start_restart_aux = set_time (restart_start_secs_aux, restart_start_days_aux)
    Time_restart_aux = Time_atmos + Time_start_restart_aux
    Time_restart_end_aux = Time_restart_aux + Time_duration_restart_aux

    intrm_rst = .false.
    intrm_rst_1step = .false.
    if (restart_days > 0 .or. restart_secs > 0) intrm_rst = .true.
    if (intrm_rst .and. restart_start_secs == 0 .and. &
        restart_start_days == 0) intrm_rst_1step = .true.

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time_atmos ) call error_mesg ('program coupler',  &
                    'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * Time_step_ocean /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of ocean time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_ocean )  &
         call error_mesg ('program coupler',   &
         'atmos time step is not a multiple of the ocean time step', FATAL)

!------ initialize component models ------
    call  atmos_model_init (Atm,  Time_init, Time_atmos, Time_step_atmos, &
                            iau_offset)

    call print_memuse_stats('after atmos model init')

!------ initialize data_override -----
    if (.NOT.Atm%bounded_domain) call data_override_init (Atm_domain_in  = Atm%domain)

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --
    if ( mpp_pe() == 0) then
       open(newunit = ascii_unit, file='RESTART/file', status='replace', form='formatted')
       close(ascii_unit,status="delete")
    endif
!-----------------------------------------------------------------------

   end subroutine coupler_init

!#######################################################################
   subroutine coupler_restart(time_stamp)
    character(len=32), intent(in), optional :: time_stamp

    integer :: restart_unit, date(6)
    character(len=128)                      :: file_res

!----- compute current date ------

      call get_date (Time_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))

!----- write restart file ------

    ! write restart file_name
    file_res = 'RESTART/coupler.res'
    if (present(time_stamp)) then
      file_res = 'RESTART/'//trim(time_stamp)//'.coupler.res'
    endif

    if ( mpp_pe().EQ.mpp_root_pe()) then
       open(newunit = restart_unit, file=file_res, status='replace', form='formatted')
       write(restart_unit, '(i6,8x,a)' )calendar_type, &
            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

       write(restart_unit, '(6i6,8x,a)' )date_init, &
            'Model start time:   year, month, day, hour, minute, second'
       write(restart_unit, '(6i6,8x,a)' )date, &
            'Current model time: year, month, day, hour, minute, second'
       close(restart_unit)
    endif

   end subroutine coupler_restart

!#######################################################################

   subroutine coupler_end

   integer :: unit, date(6)
!-----------------------------------------------------------------------

      call atmos_model_end (Atm)


      call get_date (Time_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))

!----- check time versus expected ending time ----

      if (Time_atmos /= Time_end) call error_mesg ('program coupler',  &
              'final time does not match expected ending time', WARNING)

!----- write restart file ------
    call coupler_restart()

!----- final output of diagnostic fields ----0
   call diag_manager_end (Time_atmos)

!----- to be removed once fms_io is fully deprecated -----
   call fms_io_exit()

!-----------------------------------------------------------------------

   end subroutine coupler_end

!#######################################################################

end program coupler_main

