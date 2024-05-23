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
!> \file
!! \brief Main driver program for full coupler. Provides the capability to couple component
!! models (atmosphere, land, sea ice, and ocean).
!!
!! Please see the [**main page**](index.html) for additional information.

!> \mainpage
!!
!! \brief FMS Coupler provides the capability to couple component models (atmosphere, land, sea ice, and ocean)
!! on different logically rectangular grids.
!!
!! This repository holds 3 separate directories with driver programs for different usages
!! along with modules with routines for common operations.
!!
!! There are currently 3 `coupler_main` driver programs, each with their own directory:
!!   - the original 'full' coupler_main
!!   - a slimmed down 'simple' version
!!   - a [**SHiELD**](https://www.gfdl.noaa.gov/shield/) version for use with the model
!!
!! Additionally, files in the 'shared' directory holds modules used by multiple drivers.
!! The information below is provided for the full coupler, but there is considerable overlap between the other
!! versions. Documentation on all programs and modules is available through the [**files**](files.html) tab.
!!
!! \author Bruce Wyman <Bruce.Wyman@noaa.gov>
!! \author V. Balaji <V.Balaji@noaa.gov>
!!
!! [**coupler_main.F90**](full_2coupler__main_8_f90.html) couples component models for atmosphere, ocean,
!! land and sea ice on independent grids.  It also controls the time integration.
!!
!! This version couples model components representing atmosphere, ocean, land
!! and sea ice on independent grids. Each model component is represented by a
!! data type giving the instantaneous model state.
!!
!! The component models are coupled to allow implicit vertical diffusion of
!! heat and moisture at the interfaces of the atmosphere, land, and ice models.
!! As a result, the atmosphere, land, and ice models all use the same time step.
!! The atmospheric model has been separated into down and up calls that
!! correspond to the down and up sweeps of the standard tridiagonal elimination.
!!
!! The ocean interface uses explicit mixing. Fluxes to and from the ocean must
!! be passed through the ice model. This includes atmospheric fluxes as well as
!! fluxes from the land to the ocean (runoff).
!!
!! This program contains the model's main time loop. Each iteration of the
!! main time loop is one coupled (slow) time step. Within this slow time step
!! loop is a fast time step loop, using the atmospheric time step, where the
!! tridiagonal vertical diffusion equations are solved. Exchange between sea
!! ice and ocean occurs once every slow timestep.
!!
!! \section coupler_namelists  Namelists
!!
!! The three components of coupler: @ref coupler_main , flux_exchange_mod, and surface_flux_mod
!! are configured through three namelists
!! * \ref coupler_config "coupler_nml"
!! * \ref flux_exchange_conf "flux_exchange_nml"
!! * \ref surface_flux_config "surface_flux_nml"
!!
!!
!! \note
!! -# If no value is set for current_date, start_date, or calendar (or default value
!!     specified) then the value from restart file "INPUT/coupler.res" will be used.
!!     If neither a namelist value or restart file value exist the program will fail.
!! -# The actual run length will be the sum of months, days, hours, minutes, and
!!     seconds. A run length of zero is not a valid option.
!! -# The run length must be an intergal multiple of the coupling timestep dt_cpld.
!!
!! \section Main Program Example
!! Below is some pseudo-code to illustrate the runtime loop of the coupler_main drivers.
!!
!! ~~~~~~~~~~{.f90}
!! DO slow time steps (ocean)
!!    call flux_ocean_to_ice
!!
!!    call set_ice_surface_fields
!!
!!    DO fast time steps (atmos)
!!       call flux_calculation
!!
!!       call ATMOS_DOWN
!!
!!       call flux_down_from_atmos
!!
!!       call LAND_FAST
!!
!!       call ICE_FAST
!!
!!       call flux_up_to_atmos
!!
!!       call ATMOS_UP
!!    ENDDO
!!
!!    call ICE_SLOW
!!
!!    call flux_ice_to_ocean
!!
!!    call OCEAN
!! ENDDO
!! ~~~~~~~~~~

!> \page coupler_config Coupler Configuration
!!
!! coupler_main is configured via the coupler_nml namelist in the `input.nml` file.
!! The following table contains the available namelist variables.
!!
!! <table>
!!   <tr>
!!     <th>Variable Name</th>
!!     <th>Type</th>
!!     <th>Default Value</th>
!!     <th>Description</th>
!!   </tr>
!!   <tr>
!!     <td>current_date</td>
!!     <td>integer, dimension(6)</td>
!!     <td>(/0,0,0,0,0,0/)</td>
!!     <td>The date that the current integration starts with.</td>
!!   </tr>
!!   <tr>
!!     <td>force_date_from_namelist</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Flag that determines whether the namelist variable current_date should
!!       override the date in the restart file INPUT/coupler.res. If the restart
!!       file does not exist then force_date_from_namelist has not effect, the
!!       value of current_date will be used.</td>
!!   </tr>
!!   <tr>
!!     <td>calendar</td>
!!     <td>character(len=17)</td>
!!     <td>''</td>
!!     <td>The calendar type used by the current integration. Valid values are
!!       consistent with the time_manager module: 'gregorian', 'julian', 'noleap', or
!!       'thirty_day'. The value 'no_calendar' can not be used because the
!!       time_manager's date function are used. All values must be
!!       lowercase.</td>
!!   </tr>
!!   <tr>
!!     <td>months</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>The number of months that the current integration will be run for.</td>
!!   </tr>
!!   <tr>
!!     <td>days</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>The number of days that the current integration will be run for.</td>
!!   </tr>
!!   <tr>
!!     <td>hours</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>The number of hours that the current integration will be run for.</td>
!!   </tr>
!!   <tr>
!!     <td>minutes</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>The number of minutes that the current integration will be run for.</td>
!!   </tr>
!!   <tr>
!!     <td>seconds</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>The number of seconds that the current integration will be run for.</td>
!!   </tr>
!!   <tr>
!!     <td>dt_atmos</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>Atmospheric model time step in seconds, including the fast coupling
!!       with land and sea ice.</td>
!!   </tr>
!!   <tr>
!!     <td>dt_cpld</td>
!!     <td>integer</td>
!!     <td>0</td>
!!     <td>Time step in seconds for coupling between ocean and atmospheric models:
!!       must be an integral multiple of dt_atmos and dt_ocean. This is the
!!       "slow" timestep.</td>
!!   </tr>
!!   <tr>
!!     <td>do_atmos, do_ocean, do_ice, do_land, do_flux</td>
!!     <td>logical</td>
!!     <td>.TRUE.</td>
!!     <td>If true (default), that particular model component (atmos, etc.) is
!!       run. If false, the execution of that component is skipped. This is used
!!       when ALL the output fields sent by that component to the coupler have
!!       been overridden using the data_override feature. For advanced users
!!       only: if you're not sure, you should leave these values at TRUE.</td>
!!   </tr>
!!   <tr>
!!     <td>concurrent</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If true, the ocean executes concurrently with the atmosphere-land-ocean
!!       on a separate set of PEs. If false (default), the execution is serial:
!!       call atmos... followed by call ocean... If using concurrent execution,
!!       you must set one of atmos_npes and ocean_npes, see below.</td>
!!   </tr>
!!   <tr>
!!     <td>do_concurrent_radiation</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If true then radiation is done concurrently.</td>
!!   </tr>
!!   <tr>
!!     <td>atmos_npes, ocean_npes</td>
!!     <td>integer</td>
!!     <td>none</td>
!!     <td>If concurrent is set to true, we use these to set the list of PEs on
!!       which each component runs. At least one of them must be set to a number
!!       between 0 and NPES. If exactly one of these two is set non-zero, the
!!       other is set to the remainder from NPES. If both are set non-zero they
!!       must add up to NPES.</td>
!!   </tr>
!!   <tr>
!!     <td>atmos_nthreads, ocean_nthreads</td>
!!     <td>integer</td>
!!     <td>1</td>
!!     <td>We set here the number of OpenMP threads to use separately for each
!!       component.</td>
!!   </tr>
!!   <tr>
!!     <td>radiation_nthreads</td>
!!     <td>integer</td>
!!     <td>1</td>
!!     <td>Number of threads to use for the concurrent radiation
!!       when do_concurrent_radiation = .true., otherwise is equal
!!       to atmos_nthreads </td>
!!   </tr>
!!   <tr>
!!     <td>use_lag_fluxes</td>
!!     <td>logical</td>
!!     <td>.TRUE.</td>
!!     <td> If true, the ocean is forced with SBCs from one coupling timestep ago.
!!       If false, the ocean is forced with most recent SBCs.  For an old leapfrog
!!       MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable
!!       and current fluxes to be unconditionally unstable.  For dt_cpld>dt_ocean there
!!       is probably sufficient damping for MOM4.  For more modern ocean models (such as
!!       MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False.
!!       should be much more stable.</td>
!!   </tr>
!!   <tr>
!!     <td>concurrent_ice</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If true, the slow sea-ice is forced with the fluxes that were used for
!!       the fast ice processes one timestep before.  When used in conjuction with
!!       setting slow_ice_with_ocean=true, this approach allows the atmosphere and
!!       ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can be
!!       shown to ameliorate or eliminate several ice-ocean coupled instabilities.</td>
!!   </tr>
!!   <tr>
!!     <td>slow_ice_with_ocean</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If true, the slow sea-ice is advanced on the ocean processors.  Otherwise
!!       the slow sea-ice processes are on the same PEs as the fast sea-ice.</td>
!!   </tr>
!!   <tr>
!!     <td>restart_interval</td>
!!     <td>integer, dimension(6)</td>
!!     <td>(/0,0,0,0,0,0/)</td>
!!     <td>The time interval that write out intermediate restart file. The format
!!       is (yr,mo,day,hr,min,sec). When restart_interval is all zero, no
!!       intermediate restart file will be written out.</td>
!!   </tr>
!!   <tr>
!!     <td>do_debug</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>If .TRUE. print additional debugging messages.</td>
!!   </tr>
!!   <tr>
!!     <td>do_chksum</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Turns on/off checksum of certain variables.</td>
!!   </tr>
!!   <tr>
!!     <td>do_endpoint_chksum</td>
!!     <td>logical</td>
!!     <td>.TRUE.</td>
!!     <td>Report checksums of the start and end states of certain variables.</td>
!!   </tr>
!! </table>
!!
!! \note
!! -# If no value is set for current_date, start_date, or calendar (or default value specified) then the value from
!!    restart file "INPUT/coupler.res" will be used. If neither a namelist value or restart file value exist the
!!    program will fail.
!! -# The actual run length will be the sum of months, days, hours, minutes, and seconds. A run length of zero is not a
!!    valid option.
!! -# The run length must be an intergal multiple of the coupling timestep dt_cpld.

!> \throw FATAL, "no namelist value for current_date"
!!     A namelist value for current_date must be given if no restart file for
!!     coupler_main (INPUT/coupler.res) is found.
!! \throw FATAL, "invalid namelist value for calendar"
!!     The value of calendar must be 'gregorian', 'julian', 'noleap', or 'thirty_day'.
!!     See the namelist documentation.
!! \throw FATAL, "no namelist value for calendar"
!!     If no restart file is present, then a namelist value for calendar
!!     must be specified.
!! \throw FATAL, "initial time is greater than current time"
!!     If a restart file is present, then the namelist value for either
!!     current_date or start_date was incorrectly set.
!! \throw FATAL, "run length must be multiple of ocean time step"
!!     There must be an even number of ocean time steps for the requested run length.
!! \throw FATAL, "final time does not match expected ending time"
!!     This error should probably not occur because of checks done at initialization time.
program coupler_main

  !--- F90 module for OpenMP
  use omp_lib
  use FMS
  use full_coupler_mod

  use iso_fortran_env
  implicit none

  !> model defined types
  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  ! allow members of ocean type to be aliased (ap)
  type (ocean_public_type), target  :: Ocean
  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  type(atmos_land_boundary_type) :: Atmos_land_boundary
  type(atmos_ice_boundary_type)  :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type)  :: Land_ice_atmos_boundary
  type(land_ice_boundary_type)  :: Land_ice_boundary
  type(ice_ocean_boundary_type) :: Ice_ocean_boundary
  type(ocean_ice_boundary_type) :: Ocean_ice_boundary
  type(ice_ocean_driver_type), pointer :: ice_ocean_driver_CS => NULL()

  type(FmsTime_type) :: Time
  type(FmsTime_type) :: Time_step_atmos, Time_step_cpld
  type(FmsTime_type) :: Time_atmos, Time_ocean
  type(FmsTime_type) :: Time_flux_ice_to_ocean, Time_flux_ocean_to_ice

  integer :: num_atmos_calls, na
  integer :: num_cpld_calls, nc

  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ice_bc_restart => NULL()
  type(FmsNetcdfDomainFile_t), dimension(:), pointer :: Ocn_bc_restart => NULL()

  type(FmsTime_type) :: Time_restart, Time_start, Time_end
  type(FmsTime_type) :: Time_restart_current
  character(len=32) :: timestamp

  integer :: initClock, mainClock, termClock
  integer :: newClock0, newClock1, newClock2, newClock3, newClock4, newClock5, newClock7
  integer :: newClock6f, newClock6s, newClock6e, newClock10f, newClock10s, newClock10e
  integer :: newClock8, newClock9, newClock11, newClock12, newClock13, newClock14, newClocka
  integer :: newClockb, newClockc, newClockd, newClocke, newClockf, newClockg, newClockh, newClocki
  integer :: newClockj, newClockk, newClockl
  integer :: id_atmos_model_init, id_land_model_init, id_ice_model_init
  integer :: id_ocean_model_init, id_flux_exchange_init

  integer :: outunit
  character(len=80) :: text
  integer, allocatable :: ensemble_pelist(:, :)
  integer, allocatable :: slow_ice_ocean_pelist(:)
  integer :: conc_nthreads = 1
  real :: dsec, omp_sec(2)=0.0, imb_sec(2)=0.0

  !> FREDB_ID related variables
  INTEGER :: i, status, arg_count
  CHARACTER(len=256) :: executable_name, arg, fredb_id

#ifdef FREDB_ID
#define xstr(s) str(s)
#define str(s) #s
  fredb_id = xstr(FREDB_ID)
#else
#warning "FREDB_ID not defined. Continuing as normal."
  fredb_id = 'FREDB_ID was not defined (e.g. -DFREDB_ID=...) during preprocessing'
#endif

  arg_count = command_argument_count()
  DO i=0, arg_count
    CALL get_command_argument(i, arg, status=status)
    if (status .ne. 0) then
      write (error_unit,*) 'get_command_argument failed: status = ', status, ' arg = ', i
      stop 1
    end if

    if (i .eq. 0) then
      executable_name = arg
    else if (arg == '--fredb_id') then
      write (output_unit,*) TRIM(fredb_id)
      stop
    end if
  END DO

  if (arg_count .ge. 1) then
    write (error_unit,*) 'Usage: '//TRIM(executable_name)//' [--fredb_id]'
    stop 1
  end if

  call fms_mpp_init()
  !these clocks are on the global pelist
  initClock = fms_mpp_clock_id( 'Initialization' )
  call fms_mpp_clock_begin(initClock)

  call fms_init
  call fmsconstants_init
  call fms_affinity_init


  call coupler_init(Atm, Ocean, Land, Ice, Ocean_state, Atmos_land_boundary, Atmos_ice_boundary, &
    Ocean_ice_boundary, Ice_ocean_boundary, Land_ice_atmos_boundary, Land_ice_boundary, &
    Ice_ocean_driver_CS, Ice_bc_restart, Ocn_bc_restart, ensemble_pelist, slow_ice_ocean_pelist, &
    conc_nthreads, id_atmos_model_init, id_land_model_init, &
    id_ice_model_init, id_ocean_model_init, id_flux_exchange_init, mainClock, termClock, &
    Time_step_cpld, Time_step_atmos, Time_atmos, Time_ocean, num_cpld_calls, num_atmos_calls, &
    Time, Time_start, Time_end, Time_restart, Time_restart_current)

  if (do_chksum) call coupler_chksum('coupler_init+', 0, Atm, Land, Ice)

  call fms_mpp_set_current_pelist()

  call fms_mpp_clock_end (initClock) !end initialization

  call fms_mpp_clock_begin(mainClock) !begin main loop

!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop ------

  if (check_stocks >= 0) then
    call fms_mpp_set_current_pelist()
    call flux_init_stocks(Time, Atm, Land, Ice, Ocean_state)
  endif

  if (Atm%pe) then
    call fms_mpp_set_current_pelist(Atm%pelist)
    newClock1 = fms_mpp_clock_id( 'generate_sfc_xgrid' )
  endif
  if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
     call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
    newClock2 = fms_mpp_clock_id( 'flux_ocean_to_ice' )
    newClock3 = fms_mpp_clock_id( 'flux_ice_to_ocean' )
  endif
  if (Atm%pe) then
    call fms_mpp_set_current_pelist(Atm%pelist)
    newClock5 = fms_mpp_clock_id( 'ATM' )
    newClock7  = fms_mpp_clock_id( ' ATM: atmos loop' )
    newClocka  = fms_mpp_clock_id( '  A-L: atmos_tracer_driver_gather_data' )
    newClockb  = fms_mpp_clock_id( '  A-L: sfc_boundary_layer' )
    newClockl  = fms_mpp_clock_id( '  A-L: update_atmos_model_dynamics')
    if (.not. do_concurrent_radiation) then
      newClockj  = fms_mpp_clock_id( '  A-L: serial radiation' )
    endif
    newClockc  = fms_mpp_clock_id( '  A-L: update_atmos_model_down' )
    newClockd  = fms_mpp_clock_id( '  A-L: flux_down_from_atmos' )
    newClocke  = fms_mpp_clock_id( '  A-L: update_land_model_fast' )
    newClockf  = fms_mpp_clock_id( '  A-L: update_ice_model_fast' )
    newClockg  = fms_mpp_clock_id( '  A-L: flux_up_to_atmos' )
    newClockh  = fms_mpp_clock_id( '  A-L: update_atmos_model_up' )
    if (do_concurrent_radiation) then
      newClockj  = fms_mpp_clock_id( '  A-L: concurrent radiation' )
      newClocki  = fms_mpp_clock_id( '  A-L: concurrent atmos' )
    endif
    newClockk  = fms_mpp_clock_id( '  A-L: update_atmos_model_state')
    newClock8  = fms_mpp_clock_id( ' ATM: update_land_model_slow' )
    newClock9  = fms_mpp_clock_id( ' ATM: flux_land_to_ice' )
  endif
  if (Ice%pe) then
    if (Ice%fast_ice_pe) call fms_mpp_set_current_pelist(Ice%fast_pelist)
    newClock6f = fms_mpp_clock_id( ' Ice: set_ice_surface fast' )
    newClock10f = fms_mpp_clock_id( ' Ice: update_ice_model_slow fast' )

    if (Ice%slow_ice_pe) call fms_mpp_set_current_pelist(Ice%slow_pelist)
    newClock6s = fms_mpp_clock_id( ' Ice: set_ice_surface slow' )
    newClock10s = fms_mpp_clock_id( ' Ice: update_ice_model_slow slow' )
    newClock11 = fms_mpp_clock_id( ' Ice: flux_ice_to_ocean_stocks' )

    call fms_mpp_set_current_pelist(Ice%pelist)
    newClock6e = fms_mpp_clock_id( ' Ice: set_ice_surface exchange' )
    newClock10e = fms_mpp_clock_id( ' Ice: update_ice_model_slow exchange' )
  endif
  if (Ocean%is_ocean_pe) then
    call fms_mpp_set_current_pelist(Ocean%pelist)
    newClock12 = fms_mpp_clock_id( 'OCN' )
  endif
  call fms_mpp_set_current_pelist()
  newClock4 = fms_mpp_clock_id( 'flux_check_stocks' )
  newClock13 = fms_mpp_clock_id( 'intermediate restart' )
  newClock14 = fms_mpp_clock_id( 'final flux_check_stocks' )

  do nc = 1, num_cpld_calls
    if (do_chksum) call coupler_chksum('top_of_coupled_loop+', nc, Atm, Land, Ice)
    call fms_mpp_set_current_pelist()

    if (do_chksum) then
      if (Atm%pe) then
        call fms_mpp_set_current_pelist(Atm%pelist)
        call atmos_ice_land_chksum('MAIN_LOOP-', nc, Atm, Land, Ice, &
                  Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
      endif
      if (Ocean%is_ocean_pe) then
        call fms_mpp_set_current_pelist(Ocean%pelist)
        call ocean_chksum('MAIN_LOOP-', nc, Ocean, Ice_ocean_boundary)
      endif
      call fms_mpp_set_current_pelist()
    endif

    ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
    ! points when running concurrently. The calls are placed next to each other in
    ! concurrent mode to avoid multiple synchronizations within the main loop.
    ! With concurrent_ice, these only occur on the ocean PEs.
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      ! If the slow ice is on a subset of the ocean PEs, use the ocean PElist.
       call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
      call fms_mpp_clock_begin(newClock2)
       !Redistribute quantities from Ocean to Ocean_ice_boundary
       !Ice intent is In.
       !Ice is used only for accessing Ice%area and knowing if we are on an Ice pe
      call flux_ocean_to_ice( Time, Ocean, Ice, Ocean_ice_boundary )
      Time_flux_ocean_to_ice = Time
      call fms_mpp_clock_end(newClock2)

      ! Update Ice_ocean_boundary; the first iteration is supplied by restarts
      if (use_lag_fluxes) then
        call fms_mpp_clock_begin(newClock3)
        call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
        Time_flux_ice_to_ocean = Time
        call fms_mpp_clock_end(newClock3)
      endif
    endif

    if (do_chksum) then
      call coupler_chksum('flux_ocn2ice+', nc, Atm, Land, Ice)
      if (Atm%pe) then
        call fms_mpp_set_current_pelist(Atm%pelist)
        call atmos_ice_land_chksum('fluxocn2ice+', nc, Atm, Land, Ice, &
                  Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
      endif
      if (Ocean%is_ocean_pe) then
        call fms_mpp_set_current_pelist(Ocean%pelist)
        call ocean_public_type_chksum('fluxocn2ice+', nc, Ocean)
      endif
      call fms_mpp_set_current_pelist()
    endif

    ! To print the value of frazil heat flux at the right time the following block
    ! needs to sit here rather than at the end of the coupler loop.
    if (check_stocks > 0) then
      call fms_mpp_clock_begin(newClock4)
      if (check_stocks*((nc-1)/check_stocks) == nc-1 .AND. nc > 1) then
        call fms_mpp_set_current_pelist()
        call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
      endif
      call fms_mpp_clock_end(newClock4)
    endif

    if (do_ice .and. Ice%pe) then
      if (Ice%slow_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call fms_mpp_clock_begin(newClock6s)

        ! This may do data override or diagnostics on Ice_ocean_boundary.
        call flux_ocean_to_ice_finish( Time_flux_ocean_to_ice, Ice, Ocean_Ice_Boundary )

        call unpack_ocean_ice_boundary( Ocean_ice_boundary, Ice )
        if (do_chksum) call slow_ice_chksum('update_ice_slow+', nc, Ice, Ocean_ice_boundary)
        call fms_mpp_clock_end(newClock6s)
      endif

      ! This could be a point where the model is serialized if the fast and
      ! slow ice are on different PEs.
      if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
      call fms_mpp_clock_begin(newClock6e)
      call exchange_slow_to_fast_ice(Ice)
      call fms_mpp_clock_end(newClock6e)

      if (concurrent_ice) then
        ! This call occurs all ice PEs.
        call fms_mpp_clock_begin(newClock10e)
        call exchange_fast_to_slow_ice(Ice)
        call fms_mpp_clock_end(newClock10e)
      endif

      if (Ice%fast_ice_pe) then
        if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%fast_pelist)
        call fms_mpp_clock_begin(newClock6f)
        call set_ice_surface_fields(Ice)
        call fms_mpp_clock_end(newClock6f)
      endif
    endif

    if (Atm%pe) then
      if (.NOT.(do_ice .and. Ice%pe) .OR. (ice_npes .NE. atmos_npes)) &
         call fms_mpp_set_current_pelist(Atm%pelist)

      call fms_mpp_clock_begin(newClock5)
      if (do_chksum) call atmos_ice_land_chksum('set_ice_surface+', nc, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
      call fms_mpp_clock_begin(newClock1)
      call generate_sfc_xgrid( Land, Ice )
      call fms_mpp_clock_end(newClock1)

      call send_ice_mask_sic(Time)

      !-----------------------------------------------------------------------
      !   ------ atmos/fast-land/fast-ice integration loop -------

      call fms_mpp_clock_begin(newClock7)
      do na = 1, num_atmos_calls
        if (do_chksum) call atmos_ice_land_chksum('top_of_atmos_loop-', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)

        Time_atmos = Time_atmos + Time_step_atmos

        if (do_atmos) then
          call fms_mpp_clock_begin(newClocka)
          call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
          call fms_mpp_clock_end(newClocka)
        endif

        if (do_flux) then
          call fms_mpp_clock_begin(newClockb)
          call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, &
               Atm, Land, Ice, Land_ice_atmos_boundary )
          if (do_chksum)  call atmos_ice_land_chksum('sfc+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          call fms_mpp_clock_end(newClockb)
        endif

!$OMP   PARALLEL  &
!$OMP&    NUM_THREADS(conc_nthreads)  &
!$OMP&    DEFAULT(NONE)  &
!$OMP&    PRIVATE(conc_nthreads) &
!$OMP&    SHARED(atmos_nthreads, radiation_nthreads, nc, na, num_atmos_calls, atmos_npes, land_npes, ice_npes) &
!$OMP&    SHARED(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary) &
!$OMP&    SHARED(Ocean_ice_boundary) &
!$OMP&    SHARED(do_debug, do_chksum, do_atmos, do_land, do_ice, do_concurrent_radiation, omp_sec, imb_sec) &
!$OMP&    SHARED(newClockc, newClockd, newClocke, newClockf, newClockg, newClockh, newClocki, newClockj, newClockl)
!$      if (omp_get_thread_num() == 0) then
!$OMP     PARALLEL &
!$OMP&      NUM_THREADS(1) &
!$OMP&      DEFAULT(NONE) &
!$OMP&      PRIVATE(dsec) &
!$OMP&      SHARED(atmos_nthreads, radiation_nthreads, nc, na, num_atmos_calls, atmos_npes, land_npes, ice_npes) &
!$OMP&      SHARED(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary) &
!$OMP&      SHARED(Ocean_ice_boundary) &
!$OMP&      SHARED(do_debug, do_chksum, do_atmos, do_land, do_ice, do_concurrent_radiation, omp_sec, imb_sec) &
!$OMP&      SHARED(newClockc, newClockd, newClocke, newClockf, newClockg, newClockh, newClocki, newClockj, newClockl)
!$        call omp_set_num_threads(atmos_nthreads)
!$        dsec=omp_get_wtime()

          if (do_concurrent_radiation) call fms_mpp_clock_begin(newClocki)

          !      ---- atmosphere dynamics ----
          if (do_atmos) then
            call fms_mpp_clock_begin(newClockl)
            call update_atmos_model_dynamics( Atm )
            call fms_mpp_clock_end(newClockl)
          endif
          if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_dynamics', (nc-1)*num_atmos_calls+na, &
                 Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update dyn')

          !      ---- SERIAL atmosphere radiation ----
          if (.not.do_concurrent_radiation) then
            call fms_mpp_clock_begin(newClockj)
            call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
            call fms_mpp_clock_end(newClockj)
          endif
          if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_radiation(ser)', (nc-1)*num_atmos_calls+na, &
                 Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update serial rad')

          !      ---- atmosphere down ----
          if (do_atmos) then
            call fms_mpp_clock_begin(newClockc)
            call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
            call fms_mpp_clock_end(newClockc)
          endif
          if (do_chksum) call atmos_ice_land_chksum('update_atmos_down+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update down')

          call fms_mpp_clock_begin(newClockd)
          call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                                     Land_ice_atmos_boundary, &
                                     Atmos_land_boundary, &
                                     Atmos_ice_boundary )
          call fms_mpp_clock_end(newClockd)
          if (do_chksum) call atmos_ice_land_chksum('flux_down_from_atmos+', (nc-1)*num_atmos_calls+na, Atm, Land, &
                 Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)

          !      --------------------------------------------------------------
          !      ---- land model ----
          call fms_mpp_clock_begin(newClocke)
          if (do_land .AND. land%pe) then
            if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)
            call update_land_model_fast( Atmos_land_boundary, Land )
          endif
          if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Atm%pelist)
          call fms_mpp_clock_end(newClocke)
          if (do_chksum) call atmos_ice_land_chksum('update_land_fast+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update land')

          !      ---- ice model ----
          call fms_mpp_clock_begin(newClockf)
          if (do_ice .AND. Ice%fast_ice_pe) then
            if (ice_npes .NE. atmos_npes)call fms_mpp_set_current_pelist(Ice%fast_pelist)
            call update_ice_model_fast( Atmos_ice_boundary, Ice )
          endif
          if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Atm%pelist)
          call fms_mpp_clock_end(newClockf)
          if (do_chksum) call atmos_ice_land_chksum('update_ice_fast+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update ice')

          !      --------------------------------------------------------------
          !      ---- atmosphere up ----
          call fms_mpp_clock_begin(newClockg)
          call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary, &
                                 Atmos_land_boundary, Atmos_ice_boundary )
          call fms_mpp_clock_end(newClockg)
          if (do_chksum) call atmos_ice_land_chksum('flux_up2atmos+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)

          call fms_mpp_clock_begin(newClockh)
          if (do_atmos) &
            call update_atmos_model_up( Land_ice_atmos_boundary, Atm)
          call fms_mpp_clock_end(newClockh)
          if (do_chksum) call atmos_ice_land_chksum('update_atmos_up+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          if (do_debug)  call fms_memutils_print_memuse_stats( 'update up')

          call flux_atmos_to_ocean(Time_atmos, Atm, Atmos_ice_boundary, Ice)

          call flux_ex_arrays_dealloc

          !--------------
          if (do_concurrent_radiation) call fms_mpp_clock_end(newClocki)
!$        omp_sec(1) = omp_sec(1) + (omp_get_wtime() - dsec)
!$OMP END PARALLEL
!$      endif
!$      if (omp_get_thread_num() == max(0,omp_get_num_threads()-1)) then
          !      ---- atmosphere radiation ----
          if (do_concurrent_radiation) then
!$OMP PARALLEL &
!$OMP&      NUM_THREADS(1) &
!$OMP&      DEFAULT(NONE) &
!$OMP&      PRIVATE(dsec) &
!$OMP&      SHARED(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Ocean_ice_boundary, Atmos_land_boundary) &
!$OMP&      SHARED(do_chksum, do_debug, omp_sec, num_atmos_calls, na, radiation_nthreads) &
!$OMP&      SHARED(newClockj)
!$          call omp_set_num_threads(radiation_nthreads)
!$          dsec=omp_get_wtime()

            call fms_mpp_clock_begin(newClockj)
            call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
            call fms_mpp_clock_end(newClockj)
!$          omp_sec(2) = omp_sec(2) + (omp_get_wtime() - dsec)
!---CANNOT PUT AN MPP_CHKSUM HERE AS IT REQUIRES THE ABILITY TO HAVE TWO DIFFERENT OPENMP THREADS
!---INSIDE OF MPI AT THE SAME TIME WHICH IS NOT CURRENTLY ALLOWED
!           if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_radiation(conc)', (nc-1)*num_atmos_calls+na, &
!                   Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
            if (do_debug)  call fms_memutils_print_memuse_stats( 'update concurrent rad')
!$OMP END PARALLEL
          endif
!$      endif
!$      imb_sec(omp_get_thread_num()+1) = imb_sec(omp_get_thread_num()+1) - omp_get_wtime()
!$OMP END PARALLEL
!$      imb_sec(1) = imb_sec(1) + omp_get_wtime()
!$      if (do_concurrent_radiation) imb_sec(2) = imb_sec(2) + omp_get_wtime()
!$      call omp_set_num_threads(atmos_nthreads+(conc_nthreads-1)*radiation_nthreads)

        call fms_mpp_clock_begin(newClockk)
        call update_atmos_model_state( Atm )
        if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_state+', (nc-1)*num_atmos_calls+na, Atm, Land, &
                  Ice,Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
        if (do_debug)  call fms_memutils_print_memuse_stats( 'update state')
        call fms_mpp_clock_end(newClockk)

      enddo ! end of na (fast loop)

      call fms_mpp_clock_end(newClock7)

      call fms_mpp_clock_begin(newClock8)
      !   ------ end of atmospheric time step loop -----
      if (do_land .AND. Land%pe) then
        if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)
        call update_land_model_slow(Atmos_land_boundary,Land)
      endif
      if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Atm%pelist)
      !-----------------------------------------------------------------------
      call fms_mpp_clock_end(newClock8)
      if (do_chksum) call atmos_ice_land_chksum('update_land_slow+', nc, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)

      !
      !     need flux call to put runoff and p_surf on ice grid
      !
      call fms_mpp_clock_begin(newClock9)
      call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )
      call fms_mpp_clock_end(newClock9)
      if (do_chksum) call atmos_ice_land_chksum('fluxlnd2ice+', nc, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)

      Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?
      Time = Time_atmos
      call fms_mpp_clock_end(newClock5)
    endif                     !Atm%pe block

    if(Atm%pe) then
     call fms_mpp_clock_begin(newClock5) !Ice is still using ATM pelist and need to be included in ATM clock
                                        !ATM clock is used for load-balancing the coupled models
    endif
    if (do_ice .and. Ice%pe) then

      if (Ice%fast_ice_PE) then
           if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Ice%fast_pelist)
        call fms_mpp_clock_begin(newClock10f)
       ! These two calls occur on whichever PEs handle the fast ice processess.
        call ice_model_fast_cleanup(Ice)

        call unpack_land_ice_boundary(Ice, Land_ice_boundary)
        call fms_mpp_clock_end(newClock10f)
      endif

      if (.not.concurrent_ice) then
        ! This could be a point where the model is serialized.
        if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
        ! This call occurs all ice PEs.
        call fms_mpp_clock_begin(newClock10e)
        call exchange_fast_to_slow_ice(Ice)
        call fms_mpp_clock_end(newClock10e)
      endif

      !   ------ slow-ice model ------

      ! This call occurs on whichever PEs handle the slow ice processess.
      if (Ice%slow_ice_PE .and. .not.combined_ice_and_ocean) then
        if (slow_ice_with_ocean) call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call fms_mpp_clock_begin(newClock10s)
        call update_ice_model_slow(Ice)

        call fms_mpp_clock_begin(newClock11)
        call flux_ice_to_ocean_stocks(Ice)
        call fms_mpp_clock_end(newClock11)
        call fms_mpp_clock_end(newClock10s)
      endif

      if (do_chksum) call slow_ice_chksum('update_ice_slow+', nc, Ice, Ocean_ice_boundary)
     endif  ! End of Ice%pe block

     if(Atm%pe) then
        call fms_mpp_set_current_pelist(Atm%pelist)
        call fms_mpp_clock_end(newClock5)
     endif

    ! Update Ice_ocean_boundary using the newly calculated fluxes.
    if ((concurrent_ice .OR. .NOT.use_lag_fluxes) .and. .not.combined_ice_and_ocean) then
      !this could serialize unless slow_ice_with_ocean is true.
      if ((.not.do_ice) .or. (.not.slow_ice_with_ocean)) &
        call fms_mpp_set_current_pelist()

      if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
        ! If the slow ice is on a subset of the ocean PEs, use the ocean PElist.
        call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
        call fms_mpp_clock_begin(newClock3)
        call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
        Time_flux_ice_to_ocean = Time
        call fms_mpp_clock_end(newClock3)
      endif
    endif

    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      call fms_mpp_clock_begin(newClock12)

      ! This may do data override or diagnostics on Ice_ocean_boundary.
      call flux_ice_to_ocean_finish(Time_flux_ice_to_ocean, Ice_ocean_boundary)

      if (combined_ice_and_ocean) then
        call flux_ice_to_ocean_stocks(Ice)
        call update_slow_ice_and_ocean(ice_ocean_driver_CS, Ice, Ocean_state, Ocean, &
                      Ice_ocean_boundary, Time_ocean, Time_step_cpld )
      else
      if (do_chksum) call ocean_chksum('update_ocean_model-', nc, Ocean, Ice_ocean_boundary)
      ! update_ocean_model since fluxes don't change here

      if (do_ocean) &
        call update_ocean_model( Ice_ocean_boundary, Ocean_state,  Ocean, &
                                 Time_ocean, Time_step_cpld )
      endif

      if (do_chksum) call ocean_chksum('update_ocean_model+', nc, Ocean, Ice_ocean_boundary)
      ! Get stocks from "Ice_ocean_boundary" and add them to Ocean stocks.
      ! This call is just for record keeping of stocks transfer and
      ! does not modify either Ocean or Ice_ocean_boundary
      call flux_ocean_from_ice_stocks(Ocean_state, Ocean, Ice_ocean_boundary)

      call fms_diag_send_complete(Time_step_cpld)
      Time_ocean = Time_ocean +  Time_step_cpld
      Time = Time_ocean

      call fms_mpp_clock_end(newClock12)
    endif

    !--- write out intermediate restart file when needed.
    if (Time >= Time_restart) then
      Time_restart_current = Time
      Time_restart = fms_time_manager_increment_date(Time, restart_interval(1), restart_interval(2), &
           restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
      timestamp = fms_time_manager_date_to_string(time_restart_current)
      outunit= fms_mpp_stdout()
      write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
           trim(timestamp),' is appended as prefix to each restart file name'
      if (Atm%pe) then
        call atmos_model_restart(Atm, timestamp)
        call land_model_restart(timestamp)
        call ice_model_restart(Ice, timestamp)
      endif
      if (Ocean%is_ocean_pe) then
        call ocean_model_restart(Ocean_state, timestamp)
      endif
      call coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                           Time, Time_restart_current, Time_start, Time_end, timestamp)
    endif

    !--------------
    if (do_chksum) call coupler_chksum('MAIN_LOOP+', nc, Atm, Land, Ice)
    write( text,'(a,i6)' )'Main loop at coupling timestep=', nc
    call fms_memutils_print_memuse_stats(text)
    outunit= fms_mpp_stdout()
    if (fms_mpp_pe() == fms_mpp_root_pe() .and. Atm%pe .and. do_concurrent_radiation) then
      write(outunit,102) 'At coupling step ', nc,' of ',num_cpld_calls, &
           ' Atm & Rad (imbalance): ',omp_sec(1),' (',imb_sec(1),')  ',omp_sec(2),' (',imb_sec(2),')'
    endif
    omp_sec(:)=0.
    imb_sec(:)=0.
    call flush(outunit)

  enddo
102 FORMAT(A17,i5,A4,i5,A24,f10.4,A2,f10.4,A3,f10.4,A2,f10.4,A1)

  call fms_mpp_set_current_pelist()
  call fms_mpp_clock_begin(newClock14)
  if (check_stocks >= 0) then
    call fms_mpp_set_current_pelist()
    call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
  endif
  call fms_mpp_clock_end(newClock14)

  call fms_mpp_set_current_pelist()
!-----------------------------------------------------------------------
  call fms_mpp_clock_end(mainClock)
  call fms_mpp_clock_begin(termClock)

  if (do_chksum) call coupler_chksum('coupler_end-', nc, Atm, Land, Ice)
  call coupler_end(Atm, Land, Ice, Ocean, Ocean_state, Land_ice_atmos_boundary, Atmos_ice_boundary,&
    Atmos_land_boundary, Ice_ocean_boundary, Ocean_ice_boundary, Ocn_bc_restart, Ice_bc_restart, &
    Time, Time_start, Time_end, Time_restart_current)


  call fms_mpp_clock_end(termClock)

  call fms_memutils_print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
  call fms_end

!-----------------------------------------------------------------------

end program coupler_main
