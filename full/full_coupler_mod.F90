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
module full_coupler_mod

  use omp_lib !< F90 module for OpenMP

  use FMS !, status_fms=>status
  use FMSconstants, only: fmsconstants_init

#ifdef use_deprecated_io
  use fms_io_mod,              only: fms_io_exit
#endif

! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_dynamics
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type
  use atmos_model_mod,         only: atmos_data_type_chksum
  use atmos_model_mod,         only: lnd_ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: lnd_atm_bnd_type_chksum
  use atmos_model_mod,         only: ice_atm_bnd_type_chksum
  use atmos_model_mod,         only: atmos_model_restart
  use atmos_model_mod,         only: update_atmos_model_radiation
  use atmos_model_mod,         only: update_atmos_model_state

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow
  use land_model_mod,          only: atm_lnd_bnd_type_chksum
  use land_model_mod,          only: land_data_type_chksum
  use land_model_mod,          only: land_model_restart

  use ice_model_mod,           only: ice_model_init, share_ice_domains, ice_model_end, ice_model_restart
  use ice_model_mod,           only: update_ice_model_fast, set_ice_surface_fields
  use ice_model_mod,           only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,           only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod,           only: ice_data_type_chksum, ocn_ice_bnd_type_chksum
  use ice_model_mod,           only: atm_ice_bnd_type_chksum, lnd_ice_bnd_type_chksum
  use ice_model_mod,           only: unpack_ocean_ice_boundary, exchange_slow_to_fast_ice
  use ice_model_mod,           only: ice_model_fast_cleanup, unpack_land_ice_boundary
  use ice_model_mod,           only: exchange_fast_to_slow_ice, update_ice_model_slow

  use ocean_model_mod,         only: update_ocean_model, ocean_model_init,  ocean_model_end
  use ocean_model_mod,         only: ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ocean_model_mod,         only: ocean_model_restart
  use ocean_model_mod,         only: ocean_public_type_chksum, ice_ocn_bnd_type_chksum

  use combined_ice_ocean_driver, only: update_slow_ice_and_ocean, ice_ocean_driver_type
  use combined_ice_ocean_driver, only: ice_ocean_driver_init, ice_ocean_driver_end
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod,       only: flux_exchange_init, gas_exchange_init, sfc_boundary_layer
  use flux_exchange_mod,       only: generate_sfc_xgrid, send_ice_mask_sic
  use flux_exchange_mod,       only: flux_down_from_atmos, flux_up_to_atmos
  use flux_exchange_mod,       only: flux_land_to_ice, flux_ice_to_ocean, flux_ocean_to_ice
  use flux_exchange_mod,       only: flux_ice_to_ocean_finish, flux_ocean_to_ice_finish
  use flux_exchange_mod,       only: flux_check_stocks, flux_init_stocks
  use flux_exchange_mod,       only: flux_ocean_from_ice_stocks, flux_ice_to_ocean_stocks
  use flux_exchange_mod,       only: flux_atmos_to_ocean, flux_ex_arrays_dealloc

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data

  use iso_fortran_env

  implicit none
  private

  public :: atmos_data_type, land_data_type, ice_data_type
  public :: ocean_public_type, ocean_state_type
  public :: atmos_land_boundary_type, atmos_ice_boundary_type, land_ice_atmos_boundary_type
  public :: land_ice_boundary_type, ice_ocean_boundary_type, ocean_ice_boundary_type, ice_ocean_driver_type

  public :: fmsconstants_init

  public :: update_slow_ice_and_ocean
  public :: send_ice_mask_sic
  public :: flux_ice_to_ocean_finish, flux_ice_to_ocean_stocks, flux_ocean_from_ice_stocks

  public :: atmos_model_restart, land_model_restart, ice_model_restart, ocean_model_restart

  public :: atmos_data_type_chksum,   lnd_ice_atm_bnd_type_chksum
  public :: lnd_atm_bnd_type_chksum,  ice_atm_bnd_type_chksum
  public :: atm_lnd_bnd_type_chksum,  land_data_type_chksum
  public :: ice_data_type_chksum,     ocn_ice_bnd_type_chksum
  public :: atm_ice_bnd_type_chksum,  lnd_ice_bnd_type_chksum
  public :: ocean_public_type_chksum, ice_ocn_bnd_type_chksum

  public :: coupler_init, coupler_end, coupler_restart, coupler_intermediate_restart
  public :: coupler_summarize_timestep

  public :: coupler_flux_init_finish_stocks, coupler_flux_check_stocks
  public :: coupler_flux_ocean_to_ice
  public :: coupler_unpack_ocean_ice_boundary, coupler_exchange_slow_to_fast_ice
  public :: coupler_exchange_fast_to_slow_ice, coupler_set_ice_surface_fields

  public :: coupler_generate_sfc_xgrid
  public :: coupler_atmos_tracer_driver_gather_data, coupler_sfc_boundary_layer
  public :: coupler_update_atmos_model_dynamics, coupler_update_atmos_model_down
  public :: coupler_update_atmos_model_radiation, coupler_flux_down_from_atmos
  public :: coupler_update_land_model_fast, coupler_update_ice_model_fast
  public :: coupler_flux_up_to_atmos, coupler_update_atmos_model_up
  public :: coupler_flux_atmos_to_ocean, coupler_update_atmos_model_state

  public :: coupler_update_land_model_slow, coupler_flux_land_to_ice
  public :: coupler_unpack_land_ice_boundary, coupler_flux_ice_to_ocean
  public :: coupler_update_ice_model_slow_and_stocks, coupler_update_ocean_model

  public :: coupler_clock_type, coupler_components_type, coupler_chksum_type

#include <file_version.fh>

  !> namelist interface

  !> The time interval that write out intermediate restart file.
  !! The format is (yr,mo,day,hr,min,sec).  When restart_interval
  !! is all zero, no intermediate restart file will be written out
  integer, dimension(6), public :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)

  !> The date that the current integration starts with.  (See
  !! force_date_from_namelist.)
  integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /)
  !< The calendar type used by the current integration.  Valid values are
  !! consistent with the time_manager module: 'gregorian', 'julian', 'noleap', or 'thirty_day'.
  !! The value 'no_calendar' cannot be used because the time_manager's date
  !! functions are used.  All values must be lower case.

  character(len=17) :: calendar = '                 '

  !> Flag that determines whether the namelist variable current_date should override
  !! the date in the restart file `INPUT/coupler.res`.  If the restart file does not
  !! exist then force_date_from_namelist has no effect, the value of current_date
  !! will be used.
  logical :: force_date_from_namelist = .false.

  integer, public :: months=0  !< Number of months the current integration will be run
  integer, public :: days=0    !< Number of days the current integration will be run
  integer, public :: hours=0   !< Number of hours the current integration will be run
  integer, public :: minutes=0 !< Number of minutes the current integration will be run
  integer, public :: seconds=0 !< Number of seconds the current integration will be run
  integer, public :: dt_atmos = 0 !< Atmospheric model time step in seconds, including the fast
                                  !! coupling with land and sea ice
  integer, public :: dt_cpld  = 0 !< Time step in seconds for coupling between ocean and atmospheric models.  This must
                                  !! be an integral multiple of dt_atmos and dt_ocean.  This is the "slow" timestep.
  integer, public :: atmos_npes=0 !< The number of MPI tasks to use for the atmosphere
  integer, public :: ocean_npes=0 !< The number of MPI tasks to use for the ocean
  integer, public :: ice_npes=0   !< The number of MPI tasks to use for the ice
  integer, public :: land_npes=0  !< The number of MPI tasks to use for the land
  integer, public :: atmos_nthreads=1 !< Number of OpenMP threads to use in the atmosphere
  integer, public :: ocean_nthreads=1 !< Number of OpenMP threads to use in the ocean
  integer, public :: radiation_nthreads=1 !< Number of threads to use for the radiation.

  !> Indicates if this component should be executed.  If .FALSE., then execution is skipped.
  !! This is used when ALL the output fields sent by this component to the coupler have been
  !! overridden  using the data_override feature.  This is for advanced users only.
  logical, public :: do_atmos =.true.
  logical, public :: do_land =.true. !< See do_atmos
  logical, public :: do_ice =.true.  !< See do_atmos
  logical, public :: do_ocean=.true. !< See do_atmos
  logical, public :: do_flux =.true. !< See do_atmos

  !> If .TRUE., the ocean executes concurrently with the atmosphere-land-ice on a separate
  !! set of PEs.  Concurrent should be .TRUE. if concurrent_ice is .TRUE.
  !! If .FALSE., the execution is serial: call atmos... followed by call ocean...
  logical, public :: concurrent=.FALSE.
  logical, public :: do_concurrent_radiation=.FALSE. !< If .TRUE. then radiation is done concurrently

  !> If .TRUE., the ocean is forced with SBCs from one coupling timestep ago.
  !! If .FALSE., the ocean is forced with most recent SBCs.  For an old leapfrog
  !! MOM4 coupling with dt_cpld=dt_ocean, lag fluxes can be shown to be stable
  !! and current fluxes to be unconditionally unstable.  For dt_cpld>dt_ocean there
  !! is probably sufficient damping for MOM4.  For more modern ocean models (such as
  !! MOM5, GOLD or MOM6) that do not use leapfrog timestepping, use_lag_fluxes=.False.
  !! should be much more stable.
  logical, public :: use_lag_fluxes=.TRUE.

  !> If .TRUE., the slow sea-ice is forced with the fluxes that were used for the
  !! fast ice processes one timestep before.  When used in conjuction with setting
  !! slow_ice_with_ocean=.TRUE., this approach allows the atmosphere and
  !! ocean to run concurrently even if use_lag_fluxes=.FALSE., and it can
  !! be shown to ameliorate or eliminate several ice-ocean coupled instabilities.
  logical, public :: concurrent_ice=.FALSE.

  !> If true, the slow sea-ice is advanced on the ocean processors.  Otherwise
  !! the slow sea-ice processes are on the same PEs as the fast sea-ice.
  logical, public :: slow_ice_with_ocean=.FALSE.

  !< If true, there is a single call from the coupler to advance
  !! both the slow sea-ice and the ocean. slow_ice_with_ocean and
  !! concurrent_ice must both be true if combined_ice_and_ocean is true.
  logical, public :: combined_ice_and_ocean=.FALSE.

  logical, public :: do_chksum=.FALSE.         !< If .TRUE., do multiple checksums throughout the execution of the model
  logical, public :: do_endpoint_chksum=.TRUE. !< If .TRUE., do checksums of the initial and final states.
  logical, public :: do_debug=.FALSE.!< If .TRUE. print additional debugging messages.
  integer, public :: check_stocks = 0 !< -1: never 0: at end of run only n>0: every n coupled steps
  logical, public :: use_hyper_thread = .false.

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist,         &
                         months, days, hours, minutes, seconds, dt_cpld, dt_atmos, &
                         do_atmos, do_land, do_ice, do_ocean, do_flux,             &
                         atmos_npes, ocean_npes, ice_npes, land_npes,              &
                         atmos_nthreads, ocean_nthreads, radiation_nthreads,       &
                         concurrent, do_concurrent_radiation, use_lag_fluxes,      &
                         check_stocks, restart_interval, do_debug, do_chksum,      &
                         use_hyper_thread, concurrent_ice, slow_ice_with_ocean,    &
                         do_endpoint_chksum, combined_ice_and_ocean

  !> coupler_clock_type derived type consist of all clock ids that will be set and used
  !! in full coupler_main.
  type coupler_clock_type
    integer :: initialization
    integer :: main
    integer :: generate_sfc_xgrid
    integer :: flux_ocean_to_ice
    integer :: flux_ice_to_ocean
    integer :: atm
    integer :: atmos_loop
    integer :: atmos_tracer_driver_gather_data
    integer :: sfc_boundary_layer
    integer :: update_atmos_model_dynamics
    integer :: update_atmos_model_down
    integer :: flux_down_from_atmos
    integer :: update_land_model_fast
    integer :: update_ice_model_fast
    integer :: flux_up_to_atmos
    integer :: update_atmos_model_up
    integer :: radiation
    integer :: concurrent_atmos
    integer :: update_atmos_model_state
    integer :: update_land_model_slow
    integer :: flux_land_to_ice
    integer :: set_ice_surface_fast
    integer :: update_ice_model_slow_fast
    integer :: set_ice_surface_slow
    integer :: update_ice_model_slow_slow
    integer :: flux_ice_to_ocean_stocks
    integer :: set_ice_surface_exchange
    integer :: update_ice_model_slow_exchange
    integer :: ocean
    integer :: flux_check_stocks
    integer :: intermediate_restart
    integer :: final_flux_check_stocks
    integer :: termination
    integer :: atmos_model_init
    integer :: land_model_init
    integer :: ice_model_init
    integer :: ocean_model_init
    integer :: flux_exchange_init
  end type coupler_clock_type

  type coupler_components_type
    private
    type(atmos_data_type), pointer :: Atm  !< pointer to Atm
    type(land_data_type),  pointer :: Land !< pointer to Land
    type(ice_data_type),   pointer :: Ice  !< pointer to Ice
    type(ocean_public_type), pointer :: Ocean  !< pointer to Ocean
    type(land_ice_atmos_boundary_type), pointer :: Land_ice_atmos_boundary !< pointer to Land_ice_atmos_boundary
    type(atmos_land_boundary_type), pointer :: Atmos_land_boundary !< pointer to Atmos_land_boundary
    type(atmos_ice_boundary_type),  pointer :: Atmos_ice_boundary  !< pointer to Atmos_ice_boundary
    type(land_ice_boundary_type),   pointer :: Land_ice_boundary   !< pointer to Land_ice_boundary
    type(ice_ocean_boundary_type),  pointer :: Ice_ocean_boundary  !< pointer to Ice_ocean_boundary
    type(ocean_ice_boundary_type),  pointer :: Ocean_ice_boundary  !< pointer to Ocean_ice_boundary
  contains
    procedure, public :: initialize_coupler_components_obj
    procedure, public :: get_component  !< subroutine to retrieve the requested component of an object of this type
  end type coupler_components_type

  !> The purpose of objects of coupler_chksum_type is to simplify the list
  !! of arguments required for chksum related subroutines in full_coupler_mod.
  !! The members of this type point to the model components
  type coupler_chksum_type
    private
    type(coupler_components_type), pointer :: components
  contains
    procedure, public :: initialize_coupler_chksum_obj !< associates the pointers above to model components
    procedure, public :: get_components_obj !< subroutine to retrieve the requested component of an object of this type
    procedure, public :: get_atmos_ice_land_ocean_chksums !< subroutine to compute chksums for atmos - ocean
    procedure, public :: get_atmos_ice_land_chksums !< subroutine to compute chksums for atmos_ice_land
    procedure, public :: get_slow_ice_chksums       !< subroutine to compute chskums for slow_ice
    procedure, public :: get_ocean_chksums          !< subroutine to compute chksums for ocean
    procedure, public :: get_coupler_chksums        !< subroutine to compute chksums for select fields
  end type coupler_chksum_type

  character(len=80) :: text
  character(len=48), parameter :: mod_name = 'coupler_main_mod'

  integer :: calendar_type = INVALID_CALENDAR

  !> coupled model initial date
  integer :: date_init(6) = (/ 0, 0, 0, 0, 0, 0 /)

contains

!#######################################################################

!> \brief Initialize all defined exchange grids and all boundary maps
  subroutine coupler_init(Atm, Ocean, Land, Ice, Ocean_state, Atmos_land_boundary, Atmos_ice_boundary, &
      Ocean_ice_boundary, Ice_ocean_boundary, Land_ice_atmos_boundary, Land_ice_boundary,              &
      Ice_ocean_driver_CS, Ice_bc_restart, Ocn_bc_restart, ensemble_pelist, slow_ice_ocean_pelist, conc_nthreads, &
      coupler_clocks, coupler_components_obj, coupler_chksum_obj, Time_step_cpld, Time_step_atmos, Time_atmos, Time_ocean, &
      num_cpld_calls, num_atmos_calls, Time, Time_start, Time_end, Time_restart, Time_restart_current)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm
    type(land_data_type),  intent(inout) :: Land
    type(ice_data_type),   intent(inout) :: Ice
    type(ocean_public_type), intent(inout) :: Ocean
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state
    type(atmos_land_boundary_type),  intent(inout) :: Atmos_land_boundary
    type(atmos_ice_boundary_type),   intent(inout) :: Atmos_ice_boundary
    type(ice_ocean_boundary_type),   intent(inout) :: Ice_ocean_boundary
    type(ocean_ice_boundary_type),   intent(inout) :: Ocean_ice_boundary
    type(land_ice_boundary_type),    intent(inout) :: Land_ice_boundary
    type(ice_ocean_driver_type), pointer, intent(inout) :: Ice_ocean_driver_CS
    type(land_ice_atmos_boundary_type),   intent(inout) :: Land_ice_atmos_boundary
    type(FmsNetcdfDomainFile_t), pointer, dimension(:), intent(inout) :: Ice_bc_restart, Ocn_bc_restart

    integer, intent(inout) :: conc_nthreads
    integer, allocatable, dimension(:,:), intent(inout) :: ensemble_pelist
    integer, allocatable, dimension(:),   intent(inout) :: slow_ice_ocean_pelist

    type(coupler_clock_type), intent(inout)      :: coupler_clocks
    type(coupler_components_type), intent(inout) :: coupler_components_obj
    type(coupler_chksum_type), intent(inout)     :: coupler_chksum_obj

    type(FMSTime_type), intent(inout) :: Time_step_cpld, Time_step_atmos, Time_atmos, Time_ocean
    type(FMSTime_type), intent(inout) :: Time, Time_start, Time_end, Time_restart, Time_restart_current

    integer, intent(inout) :: num_cpld_calls, num_atmos_calls
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    character(len=64), parameter    :: sub_name = 'coupler_init'
    character(len=256), parameter   :: error_header =                               &
         '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
    character(len=256), parameter   :: note_header =                                &
         '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer :: ierr, io, m, i, outunit, logunit, errunit
    integer :: date(6)
    type (FmsTime_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes

    integer :: ens_siz(6), ensemble_size
    integer :: ensemble_id = 1

    integer :: atmos_pe_start=0, atmos_pe_end=0, &
               ocean_pe_start=0, ocean_pe_end=0
    integer :: n
    integer :: diag_model_subset=DIAG_ALL
    logical :: other_fields_exist
    character(len=256) :: err_msg
    integer :: date_restart(6)
    character(len=64)  :: filename, fieldname
    integer :: id_restart, l
    character(len=8)  :: walldate
    character(len=10) :: walltime
    character(len=5)  :: wallzone
    integer           :: wallvalues(8)
    character(len=:), dimension(:), allocatable :: restart_file !< Restart file saved as a string
    integer :: time_stamp_unit !< Unif of the time_stamp file
    integer :: ascii_unit  !< Unit of a dummy ascii file

    type(FmsTime_type) :: Time_init

    type(FmsCoupler1dBC_type), pointer :: &
      gas_fields_atm => NULL(), &  ! A pointer to the type describing the
                                   ! atmospheric fields that will participate in the gas fluxes.
      gas_fields_ocn => NULL(), &  ! A pointer to the type describing the ocean
                                   ! and ice surface fields that will participate in the gas fluxes.
      gas_fluxes => NULL()  ! A pointer to the type describing the
                            ! atmosphere-ocean gas and tracer fluxes.

    integer :: num_ice_bc_restart, num_ocn_bc_restart
!-----------------------------------------------------------------------

    outunit = fms_mpp_stdout()
    errunit = fms_mpp_stderr()
    logunit = fms_mpp_stdlog()

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Entering coupler_init at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

!----- write version to logfile -------
    call fms_write_version_number('FULL_COUPLER_MOD', version)

!----- read namelist -------

    read (fms_mpp_input_nml_file, coupler_nml, iostat=io)
    ierr = check_nml_error (io, 'coupler_nml')

    !----- read date and calendar type from restart file -----
    if (fms2_io_file_exists('INPUT/coupler.res')) then
       call fms2_io_ascii_read('INPUT/coupler.res', restart_file)
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
             'no namelist value for base_date or current_date', FATAL)
      else
        date = current_date
      endif

!----- override calendar type with namelist value -----

      select case( fms_mpp_uppercase(trim(calendar)) )
      case( 'GREGORIAN' )
        calendar_type = GREGORIAN
      case( 'JULIAN' )
        calendar_type = JULIAN
      case( 'NOLEAP' )
        calendar_type = NOLEAP
      case( 'THIRTY_DAY' )
        calendar_type = THIRTY_DAY_MONTHS
      case( 'NO_CALENDAR' )
        calendar_type = NO_CALENDAR
      end select

    endif

    call fms_time_manager_set_calendar_type (calendar_type, err_msg)
    if (err_msg /= '') then
      call fms_mpp_error(FATAL, 'ERROR in coupler_init: '//trim(err_msg))
    endif

    if (concurrent .AND. .NOT.(use_lag_fluxes .OR. concurrent_ice) ) &
      call fms_mpp_error( WARNING, 'coupler_init: you have set concurrent=TRUE, &
            & use_lag_fluxes=FALSE, and concurrent_ice=FALSE &
            & in coupler_nml. When not using lag fluxes, components &
            & will synchronize at two points, and thus run serially.' )
    if (concurrent_ice .AND. .NOT.slow_ice_with_ocean ) call fms_mpp_error(WARNING, &
           'coupler_init: concurrent_ice is true, but slow ice_with_ocean is &
           & false in coupler_nml.  These two flags should both be true to avoid &
           & effectively serializing the run.' )
    if (use_lag_fluxes .AND. concurrent_ice ) call fms_mpp_error(WARNING, &
           'coupler_init: use_lag_fluxes and concurrent_ice are both true. &
           & These two coupling options are intended to be exclusive.' )

    !Check with the ensemble_manager module for the size of ensemble
    !and PE counts for each member of the ensemble.
    !
    !NOTE: ensemble_manager_init renames all the output files (restart and diagnostics)
    !      to show which ensemble member they are coming from.
    !      There also need to be restart files for each member of the ensemble in INPUT.
    !
    !NOTE: if the ensemble_size=1 the input/output files will not be renamed.
    !

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting initializing ensemble_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call fms_ensemble_manager_init() ! init pelists for ensembles
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing ensemble_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    ens_siz = fms_ensemble_manager_get_ensemble_size()
    ensemble_size = ens_siz(1)
    npes = ens_siz(2)

    !Check for the consistency of PE counts
    if (concurrent) then
!atmos_npes + ocean_npes must equal npes
      if (atmos_npes.EQ.0 ) atmos_npes = npes - ocean_npes
      if (ocean_npes.EQ.0 ) ocean_npes = npes - atmos_npes
!both must now be non-zero
      if (atmos_npes.EQ.0 .OR. ocean_npes.EQ.0 ) &
        call fms_mpp_error( FATAL, 'coupler_init: atmos_npes or ocean_npes must be specified for concurrent coupling.' )
      if (atmos_npes+ocean_npes.NE.npes ) &
        call fms_mpp_error( FATAL, 'coupler_init: atmos_npes+ocean_npes must equal npes for concurrent coupling.' )
    else                        !serial timestepping
      if ((atmos_npes.EQ.0) .and. (do_atmos .or. do_land .or. do_ice)) atmos_npes = npes
      if ((ocean_npes.EQ.0) .and. (do_ocean)) ocean_npes = npes
      if (max(atmos_npes,ocean_npes).EQ.npes) then !overlapping pelists
        ! do nothing
      else                    !disjoint pelists
        if (atmos_npes+ocean_npes.NE.npes ) call fms_mpp_error( FATAL,  &
             'coupler_init: atmos_npes+ocean_npes must equal npes for serial coupling on disjoint pelists.' )
      endif
    endif

    if (land_npes == 0 ) land_npes = atmos_npes
    if (land_npes > atmos_npes) call fms_mpp_error(FATAL, 'coupler_init: land_npes > atmos_npes')

    if (ice_npes  == 0 ) ice_npes  = atmos_npes
    if (ice_npes  > atmos_npes) call fms_mpp_error(FATAL, 'coupler_init: ice_npes > atmos_npes')

    allocate( Atm%pelist  (atmos_npes) )
    allocate( Ocean%pelist(ocean_npes) )
    allocate( Land%pelist (land_npes) )
    allocate( Ice%fast_pelist(ice_npes) )

    !Set up and declare all the needed pelists
    call fms_ensemble_manager_ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                               Atm%pelist, Ocean%pelist, Land%pelist, Ice%fast_pelist)

!set up affinities based on threads

    ensemble_id = fms_ensemble_manager_get_ensemble_id()

    if(allocated(ensemble_pelist)) call fms_mpp_error(FATAL, 'ensemble_pelist unexpectedly has already been allocated')
    allocate(ensemble_pelist(1:ensemble_size,1:npes))
    call fms_ensemble_manager_get_ensemble_pelist(ensemble_pelist)

    Atm%pe            = ANY(Atm%pelist   .EQ. fms_mpp_pe())
    Ocean%is_ocean_pe = ANY(Ocean%pelist .EQ. fms_mpp_pe())
    Land%pe           = ANY(Land%pelist  .EQ. fms_mpp_pe())

    Ice%shared_slow_fast_PEs = .not.slow_ice_with_ocean
    ! However, if using a data atmosphere and slow_ice_with_ocean then shared_slow_fast_PEs
    ! will be true. In this case, all procesors do the ocean, slow ice, and fast ice.
    if (slow_ice_with_ocean.and.(.not.do_atmos)) Ice%shared_slow_fast_PEs = .true.
    ! This is where different settings would be applied if the fast and slow
    ! ice occurred on different PEs.
    if (do_atmos) then
      if (Ice%shared_slow_fast_PEs) then
        ! Fast and slow ice processes occur on the same PEs.
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(:) = Ice%fast_pelist(:)
        allocate( Ice%slow_pelist(ice_npes) )
        Ice%slow_pelist(:) = Ice%fast_pelist(:)
        if(concurrent) then
          if(.not.allocated(slow_ice_ocean_pelist)) then
            allocate(slow_ice_ocean_pelist(ocean_npes+ice_npes))
          else
            call fms_mpp_error(FATAL, 'allocation of slow_ice_ocean_pelist unexpectedly has already been allocated')
          end if
          slow_ice_ocean_pelist(1:ice_npes) = Ice%slow_pelist(:)
          slow_ice_ocean_pelist(ice_npes+1:ice_npes+ocean_npes) = Ocean%pelist(:)
        else
          if(ice_npes .GE. ocean_npes) then
             allocate(slow_ice_ocean_pelist(ice_npes))
             slow_ice_ocean_pelist(:) = Ice%slow_pelist(:)
          else
             allocate(slow_ice_ocean_pelist(ocean_npes))
             slow_ice_ocean_pelist(:) = Ocean%pelist(:)
          endif
        endif
      else
        ! Fast ice processes occur a subset of the atmospheric PEs, while
        ! slow ice processes occur on the ocean PEs.
        allocate( Ice%slow_pelist(ocean_npes) )
        Ice%slow_pelist(:) = Ocean%pelist(:)
        allocate( Ice%pelist  (ice_npes+ocean_npes) )
        ! Set Ice%pelist() to be the union of Ice%fast_pelist and Ice%slow_pelist.
        Ice%pelist(1:ice_npes) = Ice%fast_pelist(:)
        Ice%pelist(ice_npes+1:ice_npes+ocean_npes) = Ocean%pelist(:)
        allocate(slow_ice_ocean_pelist(ocean_npes))
        slow_ice_ocean_pelist(:) = Ocean%pelist(:)
      endif
    elseif (.not.do_atmos) then
      ! In the no atmos cases, shared_slow_fast_PEs is not enough to distinguish
      ! the slow and fast ice procesor layout; slow_ice_with_ocean should be used instead.
      if (slow_ice_with_ocean) then
        ! data atmos, using combined ice-ocean driver
        ! Both fast ice and slow ice processes occur on the same PEs,
        ! since the Atmos and Ocean PEs are shared
        allocate( Ice%slow_pelist(ocean_npes) )
        Ice%slow_pelist(:) = Ocean%pelist(:)
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(1:ice_npes) = Ice%fast_pelist(:)
        allocate(slow_ice_ocean_pelist(ocean_npes))
        slow_ice_ocean_pelist(:) = Ocean%pelist(:)
      else
        ! data atmos, not using combined ice-ocean driver
        allocate( Ice%pelist  (ice_npes) )
        Ice%pelist(:) = Ice%fast_pelist(:)
        allocate( Ice%slow_pelist(ice_npes) )
        Ice%slow_pelist(:) = Ice%fast_pelist(:)
        if(ice_npes .GE. ocean_npes) then
           allocate(slow_ice_ocean_pelist(ice_npes))
           slow_ice_ocean_pelist(:) = Ice%slow_pelist(:)
        else
           allocate(slow_ice_ocean_pelist(ocean_npes))
           slow_ice_ocean_pelist(:) = Ocean%pelist(:)
        endif
      endif
    endif
    Ice%fast_ice_pe = ANY(Ice%fast_pelist(:) .EQ. fms_mpp_pe())
    Ice%slow_ice_pe = ANY(Ice%slow_pelist(:) .EQ. fms_mpp_pe())
    Ice%pe = Ice%fast_ice_pe .OR. Ice%slow_ice_pe
    call fms_mpp_declare_pelist(slow_ice_ocean_pelist)
    !--- dynamic threading turned off when affinity placement is in use
!$  call omp_set_dynamic(.FALSE.)
    !--- nested OpenMP enabled for OpenMP concurrent components
!$  call omp_set_max_active_levels(3)

    if (Atm%pe) then
      call fms_mpp_set_current_pelist( Atm%pelist )
!$    if (.not.do_concurrent_radiation) radiation_nthreads=atmos_nthreads
!$    if (do_concurrent_radiation) conc_nthreads=2
      !--- setting affinity
      if (do_concurrent_radiation) then
!$      call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads + radiation_nthreads)
!$      call omp_set_num_threads(atmos_nthreads+radiation_nthreads)
      else
!$      call fms_affinity_set('ATMOS', use_hyper_thread, atmos_nthreads)
!$      call omp_set_num_threads(atmos_nthreads)
      endif
    endif

    !> The pelists need to be set before initializing the clocks
    call coupler_set_clock_ids(coupler_clocks, Atm, Land, Ice, Ocean, ensemble_pelist, &
                               slow_ice_ocean_pelist, ensemble_id)

    !Write out messages on root PEs
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      write( text,'(a,2i6,a,i2.2)' )'Atmos PE range: ', Atm%pelist(1)  , Atm%pelist(atmos_npes)  ,&
           ' ens_', ensemble_id
      call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      if (ocean_npes .gt. 0) then
        write( text,'(a,2i6,a,i2.2)' )'Ocean PE range: ', Ocean%pelist(1), Ocean%pelist(ocean_npes), &
             ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      else
        write( text,'(a,i2.2)' )'Ocean PE range is not set (do_ocean=.false. and concurrent=.false.) for ens_', &
              ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      endif
      write( text,'(a,2i6,a,i2.2)' )'Land PE range: ', Land%pelist(1)  , Land%pelist(land_npes)  ,&
           ' ens_', ensemble_id
      call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      if (.not.concurrent_ice) then
        write( text,'(a,2i6,a,i2.2)' )'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes), &
             ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      elseif (concurrent_ice) then
        if (do_atmos) then
          write( text,'(a,2i6,a,i2.2)' )'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes+ocean_npes), &
               ' ens_', ensemble_id
          call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        elseif ((.not.do_atmos)) then
          write( text,'(a,2i6,a,i2.2)' )'Ice PE range: ', Ice%pelist(1), Ice%pelist(ice_npes), &
               ' ens_', ensemble_id
          call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        endif
        call fms_mpp_error( NOTE, 'coupler_init: Running with CONCURRENT ICE coupling.' )
        write( text,'(a,2i6,a,i2.2)' )'slow Ice PE range: ', Ice%slow_pelist(1), Ice%slow_pelist(ocean_npes), &
             ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
        write( text,'(a,2i6,a,i2.2)' )'fast Ice PE range: ', Ice%fast_pelist(1), Ice%fast_pelist(ice_npes), &
             ' ens_', ensemble_id
        call fms_mpp_error( NOTE, 'coupler_init: '//trim(text) )
      endif

      if (concurrent) then
        call fms_mpp_error( NOTE, 'coupler_init: Running with CONCURRENT coupling.' )

        write( logunit,'(a)' )'Using concurrent coupling...'
        write( logunit,'(a,4i6)' ) &
              'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
              Atm%pelist(1)  , Atm%pelist(atmos_npes), Ocean%pelist(1), Ocean%pelist(ocean_npes)
      else
        call fms_mpp_error( NOTE, 'coupler_init: Running with SERIAL coupling.' )
      endif
      if (use_lag_fluxes) then
        call fms_mpp_error( NOTE, 'coupler_init: Sending LAG fluxes to ocean.' )
      else
        call fms_mpp_error( NOTE, 'coupler_init: Sending most recent fluxes to ocean.' )
      endif
      if (concurrent_ice) call fms_mpp_error( NOTE, &
        'coupler_init: using lagged slow-ice coupling mode.')
      if (combined_ice_and_ocean) call fms_mpp_error( NOTE, &
        'coupler_init: advancing the ocean and slow-ice in a single call.')
      if (combined_ice_and_ocean .and. .not.concurrent_ice) call fms_mpp_error( FATAL, &
        'coupler_init: concurrent_ice must be true if combined_ice_and_ocean is true.')
      if (combined_ice_and_ocean .and. .not.slow_ice_with_ocean) call fms_mpp_error( FATAL, &
        'coupler_init: slow_ice_with_ocean must be true if combined_ice_and_ocean is true.')
    endif

!----- write namelist to logfile -----
    if (fms_mpp_pe() == fms_mpp_root_pe() )write( logunit, nml=coupler_nml )

!----- write current/initial date actually used to logfile file -----

    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) &
      write( logunit, 16 )date(1),trim(fms_time_manager_month_name(date(2))),date(3:6)
16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

!jwd Fork here is somewhat dangerous. It relies on "no side effects" from
!    diag_manager_init. diag_manager_init or this section should be
!    re-architected to guarantee this or remove this assumption.
!    For instance, what follows assumes that get_base_date has the same
!    time for both Atm and Ocean pes. While this should be the case, the
!    possible error condition needs to be checked

    diag_model_subset=DIAG_ALL
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      if (atmos_npes /= npes) diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
    elseif (Ocean%is_ocean_pe) then  ! Error check above for disjoint pelists should catch any problem
      call fms_mpp_set_current_pelist(Ocean%pelist)
      ! The FMS diag manager has a convention that segregates files with "ocean"
      ! in their names from the other files to handle long diag tables.  This
      ! does not work if the ice is on the ocean PEs.
      if ((ocean_npes /= npes) .and. .not.slow_ice_with_ocean) &
        diag_model_subset = DIAG_OCEAN  ! change diag_model_subset from DIAG_ALL
    endif
    if ( fms_mpp_pe() == fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    ! initialize diag_manager for processor subset output
    call fms_diag_init(DIAG_MODEL_SUBSET=diag_model_subset, TIME_INIT=date)
    call fms_memutils_print_memuse_stats( 'diag_manager_init' )
    if ( fms_mpp_pe() == fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing diag_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
!-----------------------------------------------------------------------
!------ reset pelist to "full group" ------

    call fms_mpp_set_current_pelist()
!----- always override initial/base date with diag_manager value -----

    call fms_diag_get_base_date ( date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init = fms_time_manager_set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = fms_time_manager_set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

    Time_start = Time

!----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + fms_time_manager_set_time(0,fms_time_manager_days_in_month(Time_end))
    enddo
    Time_end   = Time_end + fms_time_manager_set_time(hours*3600+minutes*60+seconds, days)
    !Need to pass Time_end into diag_manager for multiple thread case.
    call fms_diag_set_time_end(Time_end)

    Run_length = Time_end - Time

!--- get the time that last intermediate restart file was written out.
    if (fms2_io_file_exists('INPUT/coupler.intermediate.res')) then
       call fms2_io_ascii_read('INPUT/coupler.intermediate.res', restart_file)
       read(restart_file(1), *) date_restart
       deallocate(restart_file)
    else
       date_restart = date
    endif

    Time_restart_current = Time
    if (ALL(restart_interval ==0)) then
       Time_restart = fms_time_manager_increment_date(Time_end, 0, 0, 10, 0, 0, 0)   ! no intermediate restart
    else
       Time_restart = fms_time_manager_set_date(date_restart(1), date_restart(2), date_restart(3),  &
                               date_restart(4), date_restart(5), date_restart(6) )
       Time_restart = fms_time_manager_increment_date(Time_restart, restart_interval(1), restart_interval(2), &
            restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
       if (Time_restart <= Time) call fms_mpp_error(FATAL, &
            '==>Error from program coupler: The first intermediate restart time is no larger than the start time')
    endif

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) &
      open(newunit = time_stamp_unit, file='time_stamp.out', status='replace', form='formatted')

    month = fms_time_manager_month_name(date(2))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)

    call fms_time_manager_get_date (Time_end, date(1), date(2), date(3),  &
                   date(4), date(5), date(6))
    month = fms_time_manager_month_name(date(2))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) write (time_stamp_unit,20) date, month(1:3)

    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) close(time_stamp_unit)

20  format (i6,5i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_cpld  = fms_time_manager_set_time (dt_cpld ,0)
    Time_step_atmos = fms_time_manager_set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * Time_step_cpld  /= Run_length )  &
      call error_mesg ('program coupler',  &
         'run length must be multiple of coupled time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
      call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the atmos time step', FATAL)

!
!       Initialize the tracer manager. This needs to be done on all PEs,
!       before the individual models are initialized.
!

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize tracer_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call fms_tracer_manager_init()
!   Initialize the gas-exchange fluxes so this information can be made
!   available to the individual components.
    call gas_exchange_init(gas_fields_atm, gas_fields_ocn, gas_fluxes)
    call fms_coupler_types_init()
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing tracer_manager at '&
                       //trim(walldate)//' '//trim(walltime)
    endif



!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Beginning to initialize component models at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    if (Atm%pe) then
        call fms_mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
        if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Starting to initialize atmospheric model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif

        call fms_mpp_clock_begin(coupler_clocks%atmos_model_init)
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos, &
                               do_concurrent_radiation)
        call fms_mpp_clock_end(coupler_clocks%atmos_model_init)

        if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
          call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
          write(errunit,*) 'Finished initializing atmospheric model at '&
                           //trim(walldate)//' '//trim(walltime)
        endif
        call fms_memutils_print_memuse_stats( 'atmos_model_init' )
        call fms_data_override_init(Atm_domain_in = Atm%domain)
    endif
!---- land ----------
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize land model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif

      call fms_mpp_clock_begin(coupler_clocks%land_model_init)
      call land_model_init( Atmos_land_boundary, Land, Time_init, Time, &
                            Time_step_atmos, Time_step_cpld )
      call fms_mpp_clock_end(coupler_clocks%land_model_init)

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing land model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'land_model_init' )
      call fms_data_override_init(Land_domain_in = Land%domain)
#ifndef _USE_LEGACY_LAND_
      call fms_data_override_init(Land_domainUG_in = Land%ug_domain)
#endif
    endif
!---- ice -----------
    if (Ice%pe) then  ! This occurs for all fast or slow ice PEs.
      if (Ice%fast_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
      elseif (Ice%slow_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else
        call fms_mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
      endif
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize ice model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif

      call fms_mpp_clock_begin(coupler_clocks%ice_model_init)
      call ice_model_init(Ice, Time_init, Time, Time_step_atmos, &
                           Time_step_cpld, Verona_coupler=.false., &
                          concurrent_ice=concurrent_ice, &
                          gas_fluxes=gas_fluxes, gas_fields_ocn=gas_fields_ocn )
      call fms_mpp_clock_end(coupler_clocks%ice_model_init)

      ! This must be called using the union of the ice PE_lists.
      call fms_mpp_set_current_pelist(Ice%pelist)
      call share_ice_domains(Ice)

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing ice model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'ice_model_init' )
      if (Ice%fast_ice_pe) then
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
        call fms_data_override_init(Ice_domain_in = Ice%domain)
      endif
    endif

!---- ocean ---------
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize ocean model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif

      call fms_mpp_clock_begin(coupler_clocks%ocean_model_init)
      call ocean_model_init( Ocean, Ocean_state, Time_init, Time, &
                             gas_fields_ocn=gas_fields_ocn  )
      call fms_mpp_clock_end(coupler_clocks%ocean_model_init)

      if (concurrent) then
        call fms_mpp_set_current_pelist( Ocean%pelist )
!$      call fms_affinity_set('OCEAN', use_hyper_thread, ocean_nthreads)
!$      call omp_set_num_threads(ocean_nthreads)
      else
        ocean_nthreads = atmos_nthreads
        !--- omp_num_threads has already been set by the Atmos-pes, but set again to ensure
!$      call omp_set_num_threads(ocean_nthreads)
      endif

      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing ocean model at '&
                         //trim(walldate)//' '//trim(walltime)
      endif
      call fms_memutils_print_memuse_stats( 'ocean_model_init' )
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Starting to initialize data_override at '&
                         //trim(walldate)//' '//trim(walltime)
      endif
      call fms_data_override_init(Ocean_domain_in = Ocean%domain )
      if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
        write(errunit,*) 'Finished initializing data_override at '&
                         //trim(walldate)//' '//trim(walltime)
      endif

      if (combined_ice_and_ocean) &
        call ice_ocean_driver_init(ice_ocean_driver_CS, Time_init, Time)

    endif ! end of Ocean%is_ocean_pe

!---------------------------------------------
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finished initializing component models at '&
                       //trim(walldate)//' '//trim(walltime)
    endif
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))

    call fms_mpp_domains_broadcast_domain(Ice%domain)
    call fms_mpp_domains_broadcast_domain(Ice%slow_domain_NH)
    call fms_mpp_domains_broadcast_domain(Ocean%domain)
!-----------------------------------------------------------------------
!---- initialize flux exchange module ----
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Starting to initialize flux_exchange at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

    call fms_mpp_clock_begin(coupler_clocks%flux_exchange_init)
    call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
             atmos_ice_boundary, land_ice_atmos_boundary, &
             land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
         do_ocean, slow_ice_ocean_pelist, dt_atmos=dt_atmos, dt_cpld=dt_cpld)
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    call fms_mpp_clock_end(coupler_clocks%flux_exchange_init)

    call fms_mpp_set_current_pelist()
    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Finsihed initializing flux_exchange at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

    Time_atmos = Time
    Time_ocean = Time

!
!       read in extra fields for the air-sea gas fluxes
!
    if ( Ice%slow_ice_pe ) then
      call fms_mpp_set_current_pelist(Ice%slow_pelist)

      call fms_coupler_type_register_restarts(Ice%ocean_fluxes, Ice_bc_restart, &
             num_ice_bc_restart, Ice%slow_domain_NH, to_read=.true., ocean_restart=.false., directory="INPUT/")

      ! Restore the fields from the restart files
      do l = 1, num_ice_bc_restart
         if(fms2_io_check_if_open(Ice_bc_restart(l))) call fms2_io_read_restart(Ice_bc_restart(l))
      enddo

      ! Check whether the restarts were read successfully.
      call fms_coupler_type_restore_state(Ice%ocean_fluxes, use_fms2_io=.true., &
              test_by_field=.true.)

      do l = 1, num_ice_bc_restart
        if(fms2_io_check_if_open(Ice_bc_restart(l))) call fms2_io_close_file(Ice_bc_restart(l))
      enddo
    endif !< ( Ice%slow_ice_pe )

    if ( Ocean%is_ocean_pe ) then
      call fms_mpp_set_current_pelist(Ocean%pelist)

      call fms_coupler_type_register_restarts(Ocean%fields, Ocn_bc_restart, &
               num_ocn_bc_restart, Ocean%domain, to_read=.true., ocean_restart=.true., directory="INPUT/")

      ! Restore the fields from the restart files
      do l = 1, num_ocn_bc_restart
         if(fms2_io_check_if_open(Ocn_bc_restart(l))) call fms2_io_read_restart(Ocn_bc_restart(l))
      enddo

      ! Check whether the restarts were read successfully.
      call fms_coupler_type_restore_state(Ocean%fields, use_fms2_io=.true., &
              test_by_field=.true.)

      do l = 1, num_ocn_bc_restart
         if(fms2_io_check_if_open(Ocn_bc_restart(l))) call fms2_io_close_file(Ocn_bc_restart(l))
      enddo
    endif !< ( Ocean%is_ocean_pe )

    call fms_mpp_set_current_pelist()

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    if ( fms_mpp_pe().EQ.fms_mpp_root_pe() ) then
       open(newunit = ascii_unit, file='RESTART/file', status='replace', form='formatted')
       close(ascii_unit,status="delete")
    endif

    ! Call to daig_grid_end to free up memory used during regional
    ! output setup
    CALL fms_diag_grid_end()

!-----------------------------------------------------------------------

    !> Initialize coupler_components_obj memebers to point to model components
    call coupler_components_obj%initialize_coupler_components_obj(Atm, Land, Ice, Ocean, Land_ice_atmos_boundary,&
        Atmos_land_boundary, Atmos_ice_boundary, Land_ice_boundary, Ice_ocean_boundary, Ocean_ice_boundary)

    !> Initialize coupler_chksum_obj
    call coupler_chksum_obj%initialize_coupler_chksum_obj(coupler_components_obj)

    if ( do_endpoint_chksum ) then
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('coupler_init+', 0)
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call coupler_chksum_obj%get_slow_ice_chksums('coupler_init+', 0)
      end if
    end if

    call fms_mpp_set_current_pelist()
    call fms_memutils_print_memuse_stats('coupler_init')

    if (fms_mpp_pe().EQ.fms_mpp_root_pe()) then
      call DATE_AND_TIME(walldate, walltime, wallzone, wallvalues)
      write(errunit,*) 'Exiting coupler_init at '&
                       //trim(walldate)//' '//trim(walltime)
    endif

  end subroutine coupler_init

!#######################################################################

  !> This subroutine associates the pointer in an object of coupler_components_type to the model components
  subroutine initialize_coupler_components_obj(this, Atm, Land, Ice, Ocean, Land_ice_atmos_boundary, &
      Atmos_land_boundary, Atmos_ice_boundary, Land_ice_boundary, Ice_ocean_boundary, Ocean_ice_boundary)

    implicit none
    class(coupler_components_type), intent(inout) :: this !< self
    type(atmos_data_type), target, intent(in) :: Atm  !< Atm
    type(land_data_type),  target, intent(in) :: Land !< Land
    type(ice_data_type),   target, intent(in) :: Ice  !< Ice
    type(ocean_public_type), target, intent(in) :: Ocean !< Ocean
    type(land_ice_atmos_boundary_type), target, intent(in) :: Land_ice_atmos_boundary !< Land_ice_atmos_boundary
    type(atmos_land_boundary_type), target, intent(in) :: Atmos_land_boundary !< Atmos_land_boundary
    type(atmos_ice_boundary_type),  target, intent(in) :: Atmos_ice_boundary  !< Atmos_ice_boundary
    type(land_ice_boundary_type),   target, intent(in) :: Land_ice_boundary   !< Land_ice_boundary
    type(ice_ocean_boundary_type),  target, intent(in) :: Ice_ocean_boundary  !< Ice_ocean_boundary
    type(ocean_ice_boundary_type),  target, intent(in) :: Ocean_ice_boundary  !< Ocean_ice_boundary

    this%Atm => Atm
    this%Land => Land
    this%Ice => Ice
    this%Ocean => Ocean
    this%Land_ice_atmos_boundary => Land_ice_atmos_boundary
    this%Atmos_land_boundary => Atmos_land_boundary
    this%Atmos_ice_boundary => Atmos_ice_boundary
    this%Land_ice_boundary => Land_ice_boundary
    this%Ice_ocean_boundary => Ice_ocean_boundary
    this%Ocean_ice_boundary => Ocean_ice_boundary

  end subroutine initialize_coupler_components_obj

  !> Function get_component returns the requested component in the coupler_components_type object
  !! Users are required to provide the component to be retrieved as an input argument.  For example,
  !! coupler_components_obj%get_component(Atm) will return Atm = coupler_components_obj%Atm
  subroutine get_component(this, retrieve_component )

    implicit none
    class(coupler_components_type), intent(in) :: this !< the coupler_components_type object
    class(*), intent(out) :: retrieve_component  !< requested component to be retrieve.
                             !! retrieve_component can be of type atmos_data_type, land_data_type, ice_data_type,
                             !! ocean_public_type, land_ice_atmos_boundary_type, atmos_land_boundary_type,
                             !! atmos_ice_boundary_type, land_ice_boundary_type, ice_ocean_boundary_type,
                             !! ocean_ice_boundary_type

    select type(retrieve_component)
    type is(atmos_data_type) ; retrieve_component = this%Atm
    type is(land_data_type)  ; retrieve_component = this%Land
    type is(ice_data_type)   ; retrieve_component = this%Ice
    type is(ocean_public_type) ; retrieve_component = this%Ocean
    type is(land_ice_atmos_boundary_type) ; retrieve_component = this%Land_ice_atmos_boundary
    type is(atmos_land_boundary_type) ; retrieve_component = this%Atmos_land_boundary
    type is(atmos_ice_boundary_type)  ; retrieve_component = this%Atmos_ice_boundary
    type is(land_ice_boundary_type)   ; retrieve_component = this%Land_ice_boundary
    type is(ice_ocean_boundary_type)  ; retrieve_component = this%Ice_ocean_boundary
    type is(ocean_ice_boundary_type)  ; retrieve_component = this%Ocean_ice_boundary
    class default
      call fms_mpp_error(FATAL, "failure retrieving component in coupler_components_type object, &
                         cannot recognize the type of requested component")
    end select

  end subroutine get_component

  !> This subroutine associates the pointer in an object of coupler_chksum_type to the component models
  subroutine initialize_coupler_chksum_obj(this, components_obj)

    implicit none
    class(coupler_chksum_type), intent(inout) :: this
    type(coupler_components_type), intent(in), target :: components_obj

    this%components => components_obj

  end subroutine initialize_coupler_chksum_obj

  !> This subroutine retrieves coupler_chksum_obj%components_obj
  subroutine get_components_obj(this, components_obj)

    implicit none

    class(coupler_chksum_type), intent(in)     :: this  !< coupler_chksum_type
    type(coupler_components_type), intent(out) :: components_obj !< coupler_components_type to be returned

    components_obj = this%components

  end subroutine get_components_obj

  !> This subroutine finalizes the run including a final call to get_coupler_chksums if do_chksum = .True.
  !! Coupler_restart is called for the final time.
  subroutine coupler_end(Atm, Land, Ice, Ocean, Ocean_state, Land_ice_atmos_boundary, Atmos_ice_boundary,&
                         Atmos_land_boundary, Ice_ocean_boundary, Ocean_ice_boundary, Ocn_bc_restart,    &
                         Ice_bc_restart, current_timestep, Time_current, Time_start, Time_end, Time_restart_current,&
                         coupler_chksum_obj, coupler_clocks)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm  !< Atm
    type(land_data_type),  intent(inout) :: Land !< Land
    type(ice_data_type),   intent(inout) :: Ice  !< Ice
    type(ocean_public_type),          intent(inout) :: Ocean        !< Ocean
    type(ocean_state_type),  pointer, intent(inout) :: Ocean_state  !< Ocean_state
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !<Land_ice_boundary
    type(atmos_ice_boundary_type),  intent(inout) :: Atmos_ice_boundary  !< Atmos_ice_boundary
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary !< Atmos_land_boundary
    type(ice_ocean_boundary_type),  intent(inout) :: Ice_ocean_boundary  !< Ice_ocean_boundary
    type(ocean_ice_boundary_type),  intent(inout) :: Ocean_ice_boundary  !< Ocean_ice_boundary
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ocn_bc_restart !< required for coupler_restart
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ice_bc_restart !< required for coupler_restart
    integer, intent(in) :: current_timestep   !< current_timestep (nc)
    type(coupler_clock_type), intent(in)  :: coupler_clocks     !< coupler_clocks
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj !< required for chksum computations

    type(FmsTime_type), intent(in) :: Time_current   !< Current timestep
    type(FmsTime_type), intent(in) :: Time_start     !< model starting time
    type(FmsTime_type), intent(in) :: Time_end       !< model run time
    type(FmsTime_type), intent(in) :: Time_restart_current !< Time corresponding to last restart time

    call fms_mpp_clock_begin(coupler_clocks%termination)

    if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('coupler_end-', current_timestep)
    if ( do_endpoint_chksum ) then
      call coupler_chksum_obj%get_atmos_ice_land_ocean_chksums('coupler_end', 0)
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
        call coupler_chksum_obj%get_slow_ice_chksums('coupler_end', 0)
      end if
    endif
    call fms_mpp_set_current_pelist()

!----- check time versus expected ending time ----

    if (Time_current /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------
!the call to fms_io_exit has been moved here
!this will work for serial code or concurrent (disjoint pelists)
!but will fail on overlapping but unequal pelists
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      call ocean_model_end (Ocean, Ocean_state, Time_current)
    endif
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      call atmos_model_end ( Atm )
    endif
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      call land_model_end (Atmos_land_boundary, Land)
    endif
    if (Ice%pe) then  ! This happens on all fast or slow ice PEs.
      if (Ice%slow_ice_PE) then
        call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else ! This must be a fast ice PE.
        call fms_mpp_set_current_pelist(Ice%fast_pelist)
      endif
      call ice_model_end (Ice)
    endif

    !----- write restart file ------
    call coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                         Time_current, Time_restart_current, Time_start)

    call fms_diag_end (Time_current)
#ifdef use_deprecated_io
    call fms_io_exit
#endif

    call fms_mpp_set_current_pelist()
    call fms_mpp_clock_end(coupler_clocks%termination)

!-----------------------------------------------------------------------

  end subroutine coupler_end

  !>@brief Register the axis data as a variable in the netcdf file and add some dummy data.
  !! This is needed so the combiner can work correctly when the io_layout is not 1,1
  subroutine add_domain_dimension_data(fileobj)
    type(FmsNetcdfDomainFile_t) :: fileobj !< Fms2io domain decomposed fileobj
    integer, dimension(:), allocatable :: buffer !< Buffer with axis data
    integer :: is, ie !< Starting and Ending indices for data

    call fms2_io_get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call fms2_io_write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)

    call fms2_io_get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call fms2_io_write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)

  end subroutine add_domain_dimension_data


  !> \brief Writing restart file that contains running time and restart file writing time.
  subroutine coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                            Time_current, Time_restart_current, Time_start, time_stamp)

    implicit none

    type(atmos_data_type),   intent(inout) :: Atm  !< Atm
    type(ice_data_type),     intent(inout) :: Ice  !< Ice
    type(ocean_public_type), intent(inout) :: Ocean !< Ocean

    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ocn_bc_restart !< required for restarts
    type(FmsNetcdfDomainFile_t), dimension(:), pointer, intent(inout) :: Ice_bc_restart !< required for restarts
    type(FmsTime_type), intent(in) :: Time_current         !< current model runtime (Time)
    type(FmsTime_type), intent(in) :: Time_restart_current !< current restart time
    type(FmsTime_type), intent(in) :: Time_start           !< model start time
    character(len=*), intent(in),  optional :: time_stamp !< time_stamp for restart

    character(len=128) :: file_run, file_res

    integer :: yr, mon, day, hr, min, sec, date(6), n
    integer ::  num_ice_bc_restart, num_ocn_bc_restart
    integer :: restart_unit !< Unit for the coupler restart file

    call fms_mpp_set_current_pelist()

    ! write restart file
    if (present(time_stamp)) then
      file_run = 'RESTART/'//trim(time_stamp)//'.coupler.res'
      file_res = 'RESTART/'//trim(time_stamp)//'.coupler.intermediate.res'
    else
      file_run = 'RESTART/coupler.res'
      file_res = 'RESTART/coupler.intermediate.res'
    endif

    !----- compute current date ------
    call fms_time_manager_get_date (Time_current, date(1), date(2), date(3), date(4), date(5), date(6))
    if ( fms_mpp_pe().EQ.fms_mpp_root_pe()) then
       open(newunit = restart_unit, file=file_run, status='replace', form='formatted')
       write(restart_unit, '(i6,8x,a)' ) calendar_type, &
            '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

       write(restart_unit, '(6i6,8x,a)' )date_init, 'Model start time:   year, month, day, hour, minute, second'
       write(restart_unit, '(6i6,8x,a)' )date, 'Current model time: year, month, day, hour, minute, second'
       close(restart_unit)
    endif

    if (Time_restart_current > Time_start) then
      if ( fms_mpp_pe().EQ.fms_mpp_root_pe()) then
        open(newunit = restart_unit, file=file_res, status='replace', form='formatted')
        call fms_time_manager_get_date(Time_restart_current, yr,mon,day,hr,min,sec)
        write(restart_unit, '(6i6,8x,a)' )yr,mon,day,hr,min,sec, &
             'Current intermediate restart time: year, month, day, hour, minute, second'
        close(restart_unit)
      endif
    endif

    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      if (associated(Ocn_bc_restart)) deallocate(Ocn_bc_restart)

      call fms_coupler_type_register_restarts(Ocean%fields, Ocn_bc_restart, &
               num_ocn_bc_restart, Ocean%domain, to_read=.false., ocean_restart=.true., directory="RESTART/")
      do n = 1, num_ocn_bc_restart
         if (fms2_io_check_if_open(Ocn_bc_restart(n))) then
             call fms2_io_write_restart(Ocn_bc_restart(n))
             call add_domain_dimension_data(Ocn_bc_restart(n))
             call fms2_io_close_file(Ocn_bc_restart(n))
          endif
       enddo
    endif !< (Ocean%is_ocean_pe)

    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)

      if (associated(Ice_bc_restart)) deallocate(Ice_bc_restart)
      call fms_coupler_type_register_restarts(Ice%ocean_fluxes, Ice_bc_restart, &
           num_ice_bc_restart, Ice%slow_domain_NH, to_read=.false., ocean_restart=.false., directory="RESTART/")
      do n = 1, num_ice_bc_restart
        if (fms2_io_check_if_open(Ice_bc_restart(n))) then
          call fms2_io_write_restart(Ice_bc_restart(n))
          call add_domain_dimension_data(Ice_bc_restart(n))
          call fms2_io_close_file(Ice_bc_restart(n))
        endif
      enddo
    endif !< (Atm%pe)

  end subroutine coupler_restart

!--------------------------------------------------------------------------

!> \brief Print out checksums for several atm, land and ice variables
  subroutine get_coupler_chksums(this, id, timestep)

    implicit none

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id        !< id to label CHECKSUMS in stdout
    integer         , intent(in) :: timestep  !< timestep

    type :: tracer_ind_type
      integer :: atm, ice, lnd ! indices of the tracer in the respective models
    end type tracer_ind_type

    integer :: n_atm_tr, n_lnd_tr, n_exch_tr
    integer :: n_atm_tr_tot, n_lnd_tr_tot
    integer :: i, tr, n, m, outunit
    type(tracer_ind_type), allocatable :: tr_table(:)
    character(32) :: tr_name

    call fms_tracer_manager_get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, num_prog=n_atm_tr)
    call fms_tracer_manager_get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, num_prog=n_lnd_tr)

    ! Assemble the table of tracer number translation by matching names of
    ! prognostic tracers in the atmosphere and surface models; skip all atmos.
    ! tracers that have no corresponding surface tracers.
    allocate(tr_table(n_atm_tr))
    n = 1
    do i = 1,n_atm_tr
      call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, i, tr_name )
      tr_table(n)%atm = i
      tr_table(n)%ice = fms_tracer_manager_get_tracer_index ( MODEL_ICE,  tr_name )
      tr_table(n)%lnd = fms_tracer_manager_get_tracer_index ( MODEL_LAND, tr_name )
      if (tr_table(n)%ice/=NO_TRACER .or. tr_table(n)%lnd/=NO_TRACER) n = n+1
    enddo
    n_exch_tr = n-1

100 FORMAT("CHECKSUM::",A32," = ",Z20)
101 FORMAT("CHECKSUM::",A16,a,'%',a," = ",Z20)

    if (this%components%Atm%pe) then
      call fms_mpp_set_current_pelist(this%components%Atm%pelist)

      outunit = fms_mpp_stdout()
      write(outunit,*) 'BEGIN CHECKSUM(Atm):: ', id, timestep
      write(outunit,100) 'atm%t_bot',  fms_mpp_chksum(this%components%Atm%t_bot)
      write(outunit,100) 'atm%z_bot',  fms_mpp_chksum(this%components%Atm%z_bot)
      write(outunit,100) 'atm%p_bot',  fms_mpp_chksum(this%components%Atm%p_bot)
      write(outunit,100) 'atm%u_bot',  fms_mpp_chksum(this%components%Atm%u_bot)
      write(outunit,100) 'atm%v_bot',  fms_mpp_chksum(this%components%Atm%v_bot)
      write(outunit,100) 'atm%p_surf', fms_mpp_chksum(this%components%Atm%p_surf)
      write(outunit,100) 'atm%gust',   fms_mpp_chksum(this%components%Atm%gust)
      do tr = 1,n_exch_tr
         n = tr_table(tr)%atm
         if (n /= NO_TRACER) then
            call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
            write(outunit,100) 'atm%'//trim(tr_name), fms_mpp_chksum(this%components%Atm%tr_bot(:,:,n))
          endif
       enddo

      write(outunit,100) 'land%t_surf', fms_mpp_chksum(this%components%Land%t_surf)
      write(outunit,100) 'land%t_ca',   fms_mpp_chksum(this%components%Land%t_ca)
      write(outunit,100) 'land%rough_mom',   fms_mpp_chksum(this%components%Land%rough_mom)
      write(outunit,100) 'land%rough_heat',  fms_mpp_chksum(this%components%Land%rough_heat)
      write(outunit,100) 'land%rough_scale', fms_mpp_chksum(this%components%Land%rough_scale)
      do tr = 1,n_exch_tr
        n = tr_table(tr)%lnd
        if (n /= NO_TRACER) then
          call fms_tracer_manager_get_tracer_names( MODEL_ATMOS, tr_table(tr)%atm, tr_name )
#ifndef _USE_LEGACY_LAND_
          write(outunit,100) 'land%'//trim(tr_name), fms_mpp_chksum(this%components%Land%tr(:,:,n))
#else
          write(outunit,100) 'land%'//trim(tr_name), fms_mpp_chksum(this%components%Land%tr(:,:,:,n))
#endif
        endif
      enddo

      write(outunit,100) 'ice%t_surf', fms_mpp_chksum(this%components%Ice%t_surf)
      write(outunit,100) 'ice%rough_mom', fms_mpp_chksum(this%components%Ice%rough_mom)
      write(outunit,100) 'ice%rough_heat', fms_mpp_chksum(this%components%Ice%rough_heat)
      write(outunit,100) 'ice%rough_moist', fms_mpp_chksum(this%components%Ice%rough_moist)
      write(outunit,*) 'STOP CHECKSUM(Atm):: ', id, timestep

    !endif

    !if (Ocean%is_ocean_pe) call mpp_set_current_pelist(Ocean%pelist)

      write(outunit,*) 'BEGIN CHECKSUM(Ice):: ', id, timestep
      call fms_coupler_type_write_chksums(this%components%Ice%ocean_fields, outunit, 'ice%')
      write(outunit,*) 'STOP CHECKSUM(Ice):: ', id, timestep

    endif

    deallocate(tr_table)

    call fms_mpp_set_current_pelist()

  end subroutine get_coupler_chksums

  !#######################################################################

!> \brief This subroutine calls coupler_chksum as well as atmos_ice_land_chksum and ocean_chksum
  subroutine get_atmos_ice_land_ocean_chksums(this, id, timestep)

    implicit none

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !< ID labelling the set of checksums
    integer         , intent(in) :: timestep !< timestep

    if (this%components%Atm%pe) then
      call fms_mpp_set_current_pelist(this%components%Atm%pelist)
      call this%get_atmos_ice_land_chksums(trim(id), timestep)
    endif
    if (this%components%Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(this%components%Ocean%pelist)
      call this%get_ocean_chksums(trim(id), timestep)
    endif

    call fms_mpp_set_current_pelist()

  end subroutine get_atmos_ice_land_ocean_chksums

!> \brief This subroutine calls subroutine that will print out checksums of the elements
!! of the appropriate type.
!! For coupled models typically these types are not defined on all processors.
!! It is assumed that the appropriate pelist has been set before entering this routine.
!! This can be achieved in the following way.
!! ~~~~~~~~~~{.f90}
!! if (Atm%pe) then
!!    call mpp_set_current_pelist(Atm%pelist)
!!    call atmos_ice_land_chksum('MAIN_LOOP-', nc)
!! endif
!! ~~~~~~~~~~
!! If you are on the global pelist before you enter this routine using the above call,
!! you can return to the global pelist by invoking
!! ~~~~~~~~~~{.f90}
!! call mpp_set_current_pelist()
!! ~~~~~~~~~~
!! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine get_atmos_ice_land_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !< id to label CHECKSUMS in stdout
    integer         , intent(in) :: timestep !< timestep

    call atmos_data_type_chksum(     id, timestep, this%components%Atm)
    call lnd_ice_atm_bnd_type_chksum(id, timestep, this%components%Land_ice_atmos_boundary)

    if (this%components%Ice%fast_ice_pe) then
      call fms_mpp_set_current_pelist(this%components%Ice%fast_pelist)
      call ice_data_type_chksum(   id, timestep, this%components%Ice)
      call atm_ice_bnd_type_chksum(id, timestep, this%components%Atmos_ice_boundary)
    endif
    if (this%components%Land%pe) then
      call fms_mpp_set_current_pelist(this%components%Land%pelist)
      call land_data_type_chksum(  id, timestep, this%components%Land)
      call atm_lnd_bnd_type_chksum(id, timestep, this%components%Atmos_land_boundary)
    endif

    call fms_mpp_set_current_pelist(this%components%Atm%pelist)

  end subroutine get_atmos_ice_land_chksums

!> \brief This subroutine calls subroutine that will print out checksums of the elements
!! of the appropriate type.
!! For coupled models typically these types are not defined on all processors.
!! It is assumed that the appropriate pelist has been set before entering this routine.
!! This can be achieved in the following way.
!! ~~~~~~~~~~{.f90}
!! if (Ice%slow_ice_pe) then
!!    call mpp_set_current_pelist(Ice%slow_pelist)
!!    call slow_ice_chksum('MAIN_LOOP-', nc)
!! endif
!! ~~~~~~~~~~
!! If you are on the global pelist before you enter this routine using the above call,
!! you can return to the global pelist by invoking
!! ~~~~~~~~~~{.f90}
!! call mpp_set_current_pelist()
!! ~~~~~~~~~~
!! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine get_slow_ice_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !<id to label CHECKSUMS in stdout
    integer         , intent(in) :: timestep !< timestep

    call ice_data_type_chksum(    id, timestep, this%components%Ice)
    call ocn_ice_bnd_type_chksum( id, timestep, this%components%Ocean_ice_boundary)

  end subroutine get_slow_ice_chksums


!> \brief This subroutine calls subroutine that will print out checksums of the elements
!! of the appropriate type.
!! For coupled models typically these types are not defined on all processors.
!! It is assumed that the appropriate pelist has been set before entering this routine.
!! This can be achieved in the following way.
!! ~~~~~~~~~~{.f90}
!! if (Ocean%is_ocean_pe) then
!!    call mpp_set_current_pelist(Ocean%pelist)
!!    call ocean_chksum('MAIN_LOOP-', nc)
!! endif
!! ~~~~~~~~~~
!! If you are on the global pelist before you enter this routine using the above call,
!! you can return to the global pelist by invoking
!! ~~~~~~~~~~{.f90}
!! call mpp_set_current_pelist()
!! ~~~~~~~~~~
!! after you exit. This is only necessary if you need to return to the global pelist.
  subroutine get_ocean_chksums(this, id, timestep)

    class(coupler_chksum_type), intent(in) :: this !< self
    character(len=*), intent(in) :: id       !< ID labelling the set of CHECKSUMS
    integer         , intent(in) :: timestep !< Timestep

    call ocean_public_type_chksum(id, timestep, this%components%Ocean)
    call ice_ocn_bnd_type_chksum( id, timestep, this%components%Ice_ocean_boundary)

  end subroutine get_ocean_chksums

!> \brief This subroutine sets the ID for clocks used in coupler_main
  subroutine coupler_set_clock_ids(coupler_clocks, Atm, Land, Ice, Ocean, ensemble_pelist,&
                                   slow_ice_ocean_pelist, ensemble_id)

    implicit none

    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks
    type(atmos_data_type),   intent(in) :: Atm   !< Atm, required to retrieve pe information
    type(land_data_type),    intent(in) :: Land  !< Land, required to retrieve pe information
    type(ocean_public_type), intent(in) :: Ocean !< Ocean, required to retrieve pe information
    type(ice_data_type),     intent(in) :: Ice   !< Ice, required to retrieve pe information
    integer, dimension(:),   intent(in) :: slow_ice_ocean_pelist !< slow_ice_oean_pelist
    integer, dimension(:,:), intent(in) :: ensemble_pelist       !< ensemble_pelist
    integer, intent(in) :: ensemble_id !< ensemble_id used as index in ensemble_pelist

    !> initialization clock
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%atmos_model_init = fms_mpp_clock_id( '  Init: atmos_model_init ' )
    endif
    if (Land%pe) then
      call fms_mpp_set_current_pelist(Land%pelist)
      coupler_clocks%land_model_init  = fms_mpp_clock_id( '  Init: land_model_init ' )
    endif
    if (Ice%pe) then
      if (Ice%shared_slow_fast_PEs) then ; call fms_mpp_set_current_pelist(Ice%pelist)
      elseif (Ice%fast_ice_pe)      then ; call fms_mpp_set_current_pelist(Ice%fast_pelist)
      elseif (Ice%slow_ice_pe)      then ; call fms_mpp_set_current_pelist(Ice%slow_pelist)
      else ; call fms_mpp_error(FATAL, "All Ice%pes must be a part of Ice%fast_ice_pe or Ice%slow_ice_pe")
      endif
      coupler_clocks%ice_model_init   = fms_mpp_clock_id( '  Init: ice_model_init ' )
    endif
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      coupler_clocks%ocean_model_init = fms_mpp_clock_id( '  Init: ocean_model_init ' )
    endif
    call fms_mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
    coupler_clocks%flux_exchange_init = fms_mpp_clock_id( '  Init: flux_exchange_init' )

    call fms_mpp_set_current_pelist()
    coupler_clocks%main = fms_mpp_clock_id( 'Main loop' )
    coupler_clocks%termination = fms_mpp_clock_id( 'Termination' )

    If(Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%generate_sfc_xgrid = fms_mpp_clock_id( 'generate_sfc_xgrid' )
    end if
    if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
      coupler_clocks%flux_ocean_to_ice = fms_mpp_clock_id( 'flux_ocean_to_ice' )
      coupler_clocks%flux_ice_to_ocean = fms_mpp_clock_id( 'flux_ice_to_ocean' )
    endif
    if (Atm%pe) then
      call fms_mpp_set_current_pelist(Atm%pelist)
      coupler_clocks%atm         = fms_mpp_clock_id( 'ATM' )
      coupler_clocks%atmos_loop  = fms_mpp_clock_id( ' ATM: atmos loop' )
      coupler_clocks%atmos_tracer_driver_gather_data  &
          = fms_mpp_clock_id( '  A-L: atmos_tracer_driver_gather_data' )
      coupler_clocks%sfc_boundary_layer           = fms_mpp_clock_id( '  A-L: sfc_boundary_layer' )
      coupler_clocks%update_atmos_model_dynamics  = fms_mpp_clock_id( '  A-L: update_atmos_model_dynamics')
      if (.not. do_concurrent_radiation) &
          coupler_clocks%radiation            = fms_mpp_clock_id( '  A-L: serial radiation' )
      coupler_clocks%update_atmos_model_down  = fms_mpp_clock_id( '  A-L: update_atmos_model_down' )
      coupler_clocks%flux_down_from_atmos     = fms_mpp_clock_id( '  A-L: flux_down_from_atmos' )
      coupler_clocks%update_land_model_fast   = fms_mpp_clock_id( '  A-L: update_land_model_fast' )
      coupler_clocks%update_ice_model_fast    = fms_mpp_clock_id( '  A-L: update_ice_model_fast' )
      coupler_clocks%flux_up_to_atmos         = fms_mpp_clock_id( '  A-L: flux_up_to_atmos' )
      coupler_clocks%update_atmos_model_up    = fms_mpp_clock_id( '  A-L: update_atmos_model_up' )
      if (do_concurrent_radiation) then
        coupler_clocks%radiation             = fms_mpp_clock_id( '  A-L: concurrent radiation' )
        coupler_clocks%concurrent_atmos      = fms_mpp_clock_id( '  A-L: concurrent atmos' )
      endif
      coupler_clocks%update_atmos_model_state  = fms_mpp_clock_id( '  A-L: update_atmos_model_state')
      coupler_clocks%update_land_model_slow    = fms_mpp_clock_id( ' ATM: update_land_model_slow' )
      coupler_clocks%flux_land_to_ice          = fms_mpp_clock_id( ' ATM: flux_land_to_ice' )
    endif
    if (Ice%pe) then
      if (Ice%fast_ice_pe) call fms_mpp_set_current_pelist(Ice%fast_pelist)
      coupler_clocks%set_ice_surface_fast       = fms_mpp_clock_id( ' Ice: set_ice_surface fast' )
      coupler_clocks%update_ice_model_slow_fast = fms_mpp_clock_id( ' Ice: update_ice_model_slow fast' )

      if (Ice%slow_ice_pe) call fms_mpp_set_current_pelist(Ice%slow_pelist)
      coupler_clocks%set_ice_surface_slow       = fms_mpp_clock_id( ' Ice: set_ice_surface slow' )
      coupler_clocks%update_ice_model_slow_slow = fms_mpp_clock_id( ' Ice: update_ice_model_slow slow' )
      coupler_clocks%flux_ice_to_ocean_stocks   = fms_mpp_clock_id( ' Ice: flux_ice_to_ocean_stocks' )

      call fms_mpp_set_current_pelist(Ice%pelist)
      coupler_clocks%set_ice_surface_exchange       = fms_mpp_clock_id( ' Ice: set_ice_surface exchange' )
      coupler_clocks%update_ice_model_slow_exchange = fms_mpp_clock_id( ' Ice: update_ice_model_slow exchange' )

    endif
    if (Ocean%is_ocean_pe) then
      call fms_mpp_set_current_pelist(Ocean%pelist)
      coupler_clocks%ocean = fms_mpp_clock_id( 'OCN' )
    endif

    call fms_mpp_set_current_pelist()
    coupler_clocks%flux_check_stocks       = fms_mpp_clock_id( 'flux_check_stocks' )
    coupler_clocks%intermediate_restart    = fms_mpp_clock_id( 'intermediate restart' )
    coupler_clocks%final_flux_check_stocks = fms_mpp_clock_id( 'final flux_check_stocks' )

  end subroutine coupler_set_clock_ids

!> \brief This subroutine calls flux_init_stocks or does the final call to flux_check_stocks
  subroutine coupler_flux_init_finish_stocks(Time, Atm, Land, Ice, Ocean_state, &
                                             coupler_clocks, init_stocks, finish_stocks)

    implicit none

    type(FmsTime_type),    intent(in) :: Time    !< current Time
    type(atmos_data_type), intent(inout) :: Atm  !< Atm
    type(land_data_type),  intent(inout) :: Land !< Land
    type(ice_data_type),   intent(inout) :: Ice  !< Ice
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state    !< Ocean_state
    type(coupler_clock_type), intent(inout)        :: coupler_clocks !< coupler_clocks
    logical, optional, intent(in) :: init_stocks, finish_stocks  !< control flags to either call flux_init_stocks or
                                                                 !! the final flux_check_stocks

    logical :: init, finish !< control flags set to False. by default and takes on the value of init_stocks and
                            !! finish_stocks if these optional arguments are provided.
                            !! If true, either flux_init_stocks or
                            !! final flux_check_stocks will be called.

    init=.False.   ; if(present(init_stocks)) init=init_stocks
    finish=.False. ; if(present(finish_stocks)) finish=finish_stocks

    if(init) then
      call fms_mpp_set_current_pelist()
      call flux_init_stocks(Time, Atm, Land, Ice, Ocean_state)
    else if(finish) then
      call fms_mpp_set_current_pelist()
      call fms_mpp_clock_begin(coupler_clocks%final_flux_check_stocks)
      if (check_stocks >= 0) then
        call fms_mpp_set_current_pelist()
        call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
      endif
      call fms_mpp_clock_end(coupler_clocks%final_flux_check_stocks)
    else
      call fms_mpp_error(FATAL, 'coupler_flux_init_finish_stocks: either init or finish needs to be .True.')
    end if

  end subroutine coupler_flux_init_finish_stocks

  !> \brief This subroutine calls flux_check_stocks.  Clocks and pelists are set before and after
  !! call to flux_check_stocks.
  subroutine coupler_flux_check_stocks(nc, Time, Atm, Land, Ice, Ocean_state, coupler_clocks)

    implicit none

    integer, intent(in) :: nc                       !< current outerloop timestep
    type(FmsTime_type), intent(in) :: Time          !< Time
    type(atmos_data_type), intent(inout) :: Atm     !< Atm
    type(land_data_type), intent(inout)  :: Land    !< Land
    type(ice_data_type), intent(inout)   :: Ice     !< Ice
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state    !< Ocean_state
    type(coupler_clock_type), intent(inout)        :: coupler_clocks !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%flux_check_stocks)
    if (check_stocks*((nc-1)/check_stocks) == nc-1 .AND. nc > 1) then
      call fms_mpp_set_current_pelist()
      call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
    endif
    call fms_mpp_clock_end(coupler_clocks%flux_check_stocks)

  end subroutine coupler_flux_check_stocks

  !> \brief This subroutine calls flux_ocean_to_ice.
  !! Clocks and pelists are set before and after call flux_ocean_to_ice
  subroutine coupler_flux_ocean_to_ice(Ocean, Ice, Ocean_ice_boundary, coupler_clocks, slow_ice_ocean_pelist)

    implicit none

    type(ocean_public_type), intent(inout) :: Ocean  !< Ocean
    type(ice_data_type),     intent(in)    :: Ice    !< Ice
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary !< Ocean_ice_boundary
    type(coupler_clock_type), intent(inout) :: coupler_clocks          !< coupler_clocks
    integer, dimension(:),    intent(in)    :: slow_ice_ocean_pelist   !< slow_ice_ocean_pelist

    !Redistribute quantities from Ocean to Ocean_ice_boundary

    ! If the slow ice is on a subset of the ocean PEs, use the ocean PElist.
    call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)
    call fms_mpp_clock_begin(coupler_clocks%flux_ocean_to_ice)

    !Ice intent is In, used only for accessing Ice%area and knowing if we are on an Ice pe
    call flux_ocean_to_ice(Ocean, Ice, Ocean_ice_boundary)

    call fms_mpp_clock_end(coupler_clocks%flux_ocean_to_ice)

  end subroutine coupler_flux_ocean_to_ice

  !> \brief This subroutine calls flux_ocean_to_ice
  !! Clocks are set before and after call flux_ice_to_ocean. Current pelist is set when optional
  !! arguments are present and set_current_slow_ice_ocean_pelist=.True.
  subroutine coupler_flux_ice_to_ocean(Ice, Ocean, Ice_ocean_boundary, coupler_clocks, &
                                       slow_ice_ocean_pelist, set_current_slow_ice_ocean_pelist)

    implicit none

    type(ice_data_type),     intent(inout)  :: Ice     !< Ice
    type(ocean_public_type), intent(inout)  :: Ocean   !< Ocean
    type(ice_ocean_boundary_type), intent(inout) :: Ice_ocean_boundary !< Ice_ocean_boundary
    type(coupler_clock_type),      intent(inout) :: coupler_clocks     !< coupler_clocks
    integer, dimension(:), optional, intent(in) :: slow_ice_ocean_pelist  !< slow_ice_ocean_pelist
    !> if true, will call mpp_set_current_pelist(slow_ice_ocean_pelist)
    logical,               optional, intent(in) :: set_current_slow_ice_ocean_pelist

    logical :: set_current_slow_ice_ocean_pelist_in !< .F. by default; set to equal set_current_slow_ice_ocean_pelist

    !> mpp_set_current_pelist(slow_ice_ocean_pelist) is not required if coupler_flux_ice_to_ocean is being called after
    !! coupler_flux_ocean_to_ice:  mpp_set_current_pelist(slow_ice_ocean_pelist) is called
    !! in coupler_flux_ocean_to_ice
    set_current_slow_ice_ocean_pelist_in=.False.
    if(present(set_current_slow_ice_ocean_pelist)) &
        set_current_slow_ice_ocean_pelist_in = set_current_slow_ice_ocean_pelist

    ! Update Ice_ocean_boundary; the first iteration is supplied by restarts

    if(set_current_slow_ice_ocean_pelist_in) call fms_mpp_set_current_pelist(slow_ice_ocean_pelist)

    call fms_mpp_clock_begin(coupler_clocks%flux_ice_to_ocean)
    call flux_ice_to_ocean(Ice, Ocean, Ice_ocean_boundary)
    call fms_mpp_clock_end(coupler_clocks%flux_ice_to_ocean)

  end subroutine coupler_flux_ice_to_ocean

  !> \brief This subroutine calls flux_ocean_to_ice_finish and unpack_ocean_ice_boundary.
  !! Clocks and pelists are set before/after the calls.  Checksum is computed if do_chksum=.True.
  subroutine coupler_unpack_ocean_ice_boundary(nc, Time_flux_ocean_to_ice, Ice, Ocean_ice_boundary, coupler_clocks, &
                                               coupler_chksum_obj)

    implicit none

    integer,             intent(in)    :: nc                     !< Current outer loop timestep
    type(FmsTime_type),  intent(inout) :: Time_flux_ocean_to_ice !< Time flux_ocean_to_ice
    type(ice_data_type), intent(inout) :: Ice                    !< Ice
    type(ocean_ice_boundary_type), intent(inout) :: Ocean_ice_boundary  !< Ocean_ice_boundary
    type(coupler_clock_type),      intent(inout) :: coupler_clocks      !< coupler_clocks
    type(coupler_chksum_type),     intent(in)  :: coupler_chksum_obj

    call fms_mpp_set_current_pelist(Ice%slow_pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_slow)

    ! This may do data override or diagnostics on Ice_ocean_boundary.
    call flux_ocean_to_ice_finish( Time_flux_ocean_to_ice, Ice, Ocean_Ice_Boundary )
    call unpack_ocean_ice_boundary( Ocean_ice_boundary, Ice )
    if (do_chksum) call coupler_chksum_obj%get_slow_ice_chksums('update_ice_slow+', nc)

    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_slow)

  end subroutine coupler_unpack_ocean_ice_boundary

  !> This subroutine calls exchange_slow_to_fast_ice
  !! Clocks and pelists are set before/after the calls.
  subroutine coupler_exchange_slow_to_fast_ice(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice                !< Ice
    type(coupler_clock_type), intent(inout) :: coupler_clocks !<coupler_clocks

    ! This could be a point where the model is serialized if the fast and
    ! slow ice are on different PEs.
    if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_exchange)
    call exchange_slow_to_fast_ice(Ice)
    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_exchange)

  end subroutine coupler_exchange_slow_to_fast_ice

  !> \brief This subroutine calls exchange_fast_to_slow_ice.  Clocks are set before and after the call.
  !! The current pelist is set if the optional argument set_ice_current_pelist is set to true.
  subroutine coupler_exchange_fast_to_slow_ice(Ice, coupler_clocks, set_ice_current_pelist)

    implicit none
    type(ice_data_type), intent(inout) :: Ice                 !< Ice
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks
    logical, optional, intent(in) :: set_ice_current_pelist   !< If true, call mpp_set_current_pelist(Ice%pelist)

    logical :: set_ice_current_pelist_in

    set_ice_current_pelist_in = .False.
    if(present(set_ice_current_pelist)) set_ice_current_pelist_in = set_ice_current_pelist

    if(set_ice_current_pelist_in .and. .not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_exchange)
    call exchange_fast_to_slow_ice(Ice)
    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_exchange)

  end subroutine coupler_exchange_fast_to_slow_ice

!> \brief This subroutine calls set_ice_surface_fields.  Clocks and pelist are set before/after the call.
  subroutine coupler_set_ice_surface_fields(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice                 !< Ice
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks

    if (.not.Ice%shared_slow_fast_PEs) call fms_mpp_set_current_pelist(Ice%fast_pelist)
    call fms_mpp_clock_begin(coupler_clocks%set_ice_surface_fast)
    call set_ice_surface_fields(Ice)
    call fms_mpp_clock_end(coupler_clocks%set_ice_surface_fast)

  end subroutine coupler_set_ice_surface_fields

!> \brief This subroutine calls generate_sfc_xgrid.  Clocks are set and before the call
  subroutine coupler_generate_sfc_xgrid(Land, Ice, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land   !< Land
    type(ice_data_type),  intent(inout) :: Ice    !< Ice
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%generate_sfc_xgrid)
    call generate_sfc_xgrid( Land, Ice )
    call fms_mpp_clock_end(coupler_clocks%generate_sfc_xgrid)

  end subroutine coupler_generate_sfc_xgrid

  !> \brief This subroutine calls atmo_tracer_driver_gather_data.
  !! Clocks are set before and after the call.
  subroutine coupler_atmos_tracer_driver_gather_data(Atm, coupler_clocks)

    implicit none

    type(atmos_data_type), intent(inout)    :: Atm !< Atm
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks
    call fms_mpp_clock_begin(coupler_clocks%atmos_tracer_driver_gather_data)
    call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
    call fms_mpp_clock_end(coupler_clocks%atmos_tracer_driver_gather_data)

  end subroutine coupler_atmos_tracer_driver_gather_data

  !> \brief This subroutine calls coupler_sfc_boundary_layer.  Chksums are computed
  !! if do_chksum = .True.  Clocks are set for runtime statistics.
  subroutine coupler_sfc_boundary_layer(Atm, Land, Ice, Land_ice_atmos_boundary, &
                                        Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm  !< Atm
    type(land_data_type), intent(inout)  :: Land !< Land
    type(ice_data_type), intent(inout)   :: Ice  !< Ice
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !< Land_ice_atmos_boundary
    type(FmsTime_type), intent(in) :: Time_atmos           !< Atmos time
    integer, intent(in)            :: current_timestep     !< (nc-1)*num_atmos_cal + na
    type(coupler_chksum_type), intent(in)   :: coupler_chksum_obj
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%sfc_boundary_layer)

    call sfc_boundary_layer( real(dt_atmos), Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary )
    if(do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('sfc+', current_timestep)

    call fms_mpp_clock_end(coupler_clocks%sfc_boundary_layer)

  end subroutine coupler_sfc_boundary_layer

  !> This subroutine calls update_atmos_model_dynamics.  Clocks are set for runtime statistics.  Chksums
  !! and memory usage are computed if do_chksum and do_debug are .True.
  subroutine coupler_update_atmos_model_dynamics(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm !< Atm
    integer,                   intent(in) :: current_timestep   !< Current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj !< coupler_chksum_obj pointing to component types
    type(coupler_clock_type),  intent(inout) :: coupler_clocks  !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_dynamics)
    call update_atmos_model_dynamics(Atm)
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_dynamics)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_dynamics', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update dyn')

  end subroutine coupler_update_atmos_model_dynamics

  !> This subroutine calls update_atmos_model_radiation.  Clocks are set for runtime statistics.
  !! Chksums are computed if do_chksum is .True. and do_concurrent_radiation is .False..  Memory
  !! usage is computed if do_debug is .True.
  subroutine coupler_update_atmos_model_radiation(Atm, Land_ice_atmos_boundary, coupler_clocks, &
                                                  current_timestep, coupler_chksum_obj)

    implicit none

    type(atmos_data_type), intent(inout) :: Atm !< Atm
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !< Land_ice_atmos_boundary
    type(coupler_clock_type),         intent(inout) :: coupler_clocks     !< coupler_clocks
    integer,                   optional, intent(in) :: current_timestep   !< Current timestep
    type(coupler_chksum_type), optional, intent(in) :: coupler_chksum_obj !< points to component types

    character(128) :: memuse_stats_id = 'update serial rad' !< used to label mem usage

    call fms_mpp_clock_begin(coupler_clocks%radiation)
    call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
    call fms_mpp_clock_end(coupler_clocks%radiation)

    if(do_chksum) then
      !> cannot put mpp_chksum for concurrent_radiation as it requires the ability to have two different OpenMP threads
      !! inside of MPI at the same time which is not currently allowed
      if(.not.do_concurrent_radiation) &
          call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_radiation(ser)',current_timestep)
    end if

    if (do_debug) then
      if(do_concurrent_radiation) memuse_stats_id = 'update concurrent rad'
      call fms_memutils_print_memuse_stats(trim(memuse_stats_id))
    end if

  end subroutine coupler_update_atmos_model_radiation

  !> This subroutine calls update_atmos_model_down.  Clocks are set for runtime statistics.  Chksums
  !! and memory usage are computed if do_chksum and do_debug are .True.
  subroutine coupler_update_atmos_model_down(Atm, Land_ice_atmos_boundary, current_timestep, &
                                             coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm  !< Atm
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !<Land ice_atmos_boundary
    integer,                   intent(in) :: current_timestep   !< Current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj !< coupler_chksum_obj pointing to component types
    type(coupler_clock_type),  intent(inout) :: coupler_clocks  !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_down)
    call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_down)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_down+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update down')

  end subroutine coupler_update_atmos_model_down

  !> This subroutine calls flux_down_from_atmos.  Clocks are set for runtime statistics.  Chksums
  !! are computed if do_chksum = .True.
  subroutine coupler_flux_down_from_atmos(Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, &
              Atmos_ice_boundary, Time_atmos, current_timestep, coupler_clocks, coupler_chksum_obj)

    implicit none
    type(atmos_data_type), intent(inout) :: Atm  !< Atm
    type(land_data_type),  intent(inout) :: Land !< Land
    type(ice_data_type),   intent(inout) :: Ice  !< Ice
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !< Land_ice_atmos_boundary
    type(atmos_land_boundary_type),     intent(inout) :: Atmos_land_boundary     !< Atmos_land_boundary
    type(atmos_ice_boundary_type),      intent(inout) :: Atmos_ice_boundary      !< Atmos_ice_boundary
    type(FmsTime_type), intent(in) :: Time_atmos       !<Time_atmos FmsTime_type containing time in seconds
    integer,            intent(in) :: current_timestep !< current_timestep
    type(coupler_clock_type), intent(inout) :: coupler_clocks !<coupler_clocks
    type(coupler_chksum_type), intent(in)   :: coupler_chksum_obj !< used to compute chksum

    call fms_mpp_clock_begin(coupler_clocks%flux_down_from_atmos)
    call flux_down_from_atmos(Time_atmos, Atm, Land, Ice, Land_ice_atmos_boundary, &
                              Atmos_land_boundary, Atmos_ice_boundary )
    call fms_mpp_clock_end(coupler_clocks%flux_down_from_atmos)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('flux_down_from_atmos+', current_timestep)

  end subroutine coupler_flux_down_from_atmos

  !> This subroutine calls update_land_model_fast.  Clocks are set for runtime statistics.  Chksums
  !! and memory usage are computed if do_chksum and do_debug are .True.
  subroutine coupler_update_land_model_fast(Land, Atmos_land_boundary, atm_pelist, current_timestep, &
                                            coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type),           intent(inout) :: Land !< Land
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary !< Atmos_land_boundary
    integer, dimension(:), intent(in) :: atm_pelist !< Atm%pelist to reset the pelist to Atm%pelist
    integer,                   intent(in) :: current_timestep       !< current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj     !< points to component types
    type(coupler_clock_type),  intent(inout) :: coupler_clocks      !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_land_model_fast) !< current pelist=Atm%pelist
    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)

    call update_land_model_fast( Atmos_land_boundary, Land )

    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_land_model_fast)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_land_fast+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update land')

  end subroutine coupler_update_land_model_fast

  !> This subroutine calls update_ice_model_fast.  Clocks are set for runtime statistics.  Chksums
  !! and memory usage are computed if do_chksum and do_debug are .True.
  subroutine coupler_update_ice_model_fast(Ice, Atmos_ice_boundary, atm_pelist, current_timestep, &
                                           coupler_chksum_obj, coupler_clocks)

    implicit none
    type(ice_data_type),           intent(inout) :: Ice                !< Ice
    type(Atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary !< Atmos_ice_boundary
    integer, dimension(:), intent(in) :: atm_pelist !< Atm%pelist to reset the pelist to Atm%pelist
    integer,                     intent(in) :: current_timestep   !< current_timestep
    type(coupler_chksum_type),   intent(in) :: coupler_chksum_obj !< points to component types
    type(coupler_clock_type), intent(inout) :: coupler_clocks     !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_fast)  !< current pelist = Atm%pelist
    if (ice_npes .NE. atmos_npes)call fms_mpp_set_current_pelist(Ice%fast_pelist)

    call update_ice_model_fast( Atmos_ice_boundary, Ice )

    if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_ice_model_fast)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_ice_fast+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update ice')

  end subroutine coupler_update_ice_model_fast

  !> This subroutine calls flux_up_to_atmos.  Clocks are set for runtime statistics.  Chksums
  !! are computed if do_chksum is .True.
  subroutine coupler_flux_up_to_atmos(Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary,&
                                      Time_atmos, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land !< Land
    type(ice_data_type),  intent(inout) :: Ice  !< Ice
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary !< Land_ice_atmos_boundary
    type(atmos_land_boundary_type),     intent(inout) :: Atmos_land_boundary     !< Atmos_land_boundary
    type(atmos_ice_boundary_type),      intent(inout) :: Atmos_ice_boundary      !< Atmos_ice_boundary
    type(FmsTime_type), intent(in) :: Time_atmos         !< Time_atmos, time in seconds
    integer,            intent(in) :: current_timestep   !< current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj !< points to component types
    type(coupler_clock_type),  intent(in) :: coupler_clocks     !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%flux_up_to_atmos)
    call flux_up_to_atmos(Time_atmos, Land, Ice, Land_ice_atmos_boundary, Atmos_land_boundary, Atmos_ice_boundary)
    call fms_mpp_clock_end(coupler_clocks%flux_up_to_atmos)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('flux_up2atmos+', current_timestep)

  end subroutine coupler_flux_up_to_atmos

  !> This subroutine calls update_atmos_model_up.  Clocks are set for runtime statistics.  Chksums
  !! and memory usage are computed if do_chksum and do_debug are .True.
  subroutine coupler_update_atmos_model_up(Atm, Land_ice_atmos_boundary, current_timestep, &
                                           coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type),              intent(inout) :: Atm                      !< Atm
    type(land_ice_atmos_boundary_type), intent(inout) :: Land_ice_atmos_boundary  !< Land_ice_atmos_boundary
    integer,                  intent(in) :: current_timestep   !< current_timestep
    type(coupler_chksum_type),intent(in) :: coupler_chksum_obj !< points to component types
    type(coupler_clock_type), intent(inout) :: coupler_clocks  !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_up)
    call update_atmos_model_up(Land_ice_atmos_boundary, Atm)
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_up)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_up+', current_timestep)
    if (do_debug) call fms_memutils_print_memuse_stats( 'update up')

  end subroutine coupler_update_atmos_model_up

  !> This subroutine calls flux_atmos_to_ocean and calls flux_ex_arrays_dealloc
  subroutine coupler_flux_atmos_to_ocean(Atm, Atmos_ice_boundary, Ice, Time_atmos)

    implicit none
    type(atmos_data_type),         intent(inout) :: Atm  !< Atm
    type(atmos_ice_boundary_type), intent(inout) :: Atmos_ice_boundary !< Atmos_ice_boundary
    type(ice_data_type), intent(inout) :: Ice        !< Ice
    type(FmsTime_type),  intent(in)    :: Time_atmos !< Time in seconds

    call flux_atmos_to_ocean(Time_atmos, Atm, Atmos_ice_boundary, Ice)
    call flux_ex_arrays_dealloc

  end subroutine coupler_flux_atmos_to_ocean

  !> This subroutine calls update_atmos_model_state.  Chksums are mem usage are computed
  !! if do_chksum and do_debug are .True. respectively
  subroutine coupler_update_atmos_model_state(Atm, current_timestep, coupler_chksum_obj, coupler_clocks)

    implicit none
    type(atmos_data_type), intent(inout)  :: Atm               !< Atm
    integer,               intent(in)     :: current_timestep  !< current_timestep
    type(coupler_chksum_type), intent(in)    :: coupler_chksum_obj !< used to compute chksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks     !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_atmos_model_state)
    call update_atmos_model_state( Atm )
    call fms_mpp_clock_end(coupler_clocks%update_atmos_model_state)

    if (do_chksum) &
        call coupler_chksum_obj%get_atmos_ice_land_chksums('update_atmos_model_state+', current_timestep)
    if (do_debug)  call fms_memutils_print_memuse_stats( 'update state')

  end subroutine coupler_update_atmos_model_state

  !> In this subroutine, update_land model_slow is called by the Land%pes.  The atm_pelist are
  !! only required to set the clocks.  Chksums are computed if do_chksum = .True.
  subroutine coupler_update_land_model_slow(Land, Atmos_land_boundary, atm_pelist, current_timestep, &
                                            coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type),           intent(inout) :: Land                 !< Land
    type(atmos_land_boundary_type), intent(inout) :: Atmos_land_boundary  !< Atmos_land_boundary
    integer, dimension(:), intent(in) :: atm_pelist                !< atm_pelist used for clocks
    integer, intent(in) :: current_timestep                        !< current timestep
    type(coupler_chksum_type), intent(in)    :: coupler_chksum_obj !< coupler_chksum_obj for chksum computation
    type(coupler_clock_type),  intent(inout) :: coupler_clocks     !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%update_land_model_slow)

    if (Land%pe) then
      if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Land%pelist)
      call update_land_model_slow(Atmos_land_boundary,Land)
    endif

    if (land_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(atm_pelist)
    call fms_mpp_clock_end(coupler_clocks%update_land_model_slow)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('update_land_slow+', current_timestep)

  end subroutine coupler_update_land_model_slow

  !> This subroutine calls flux_land_to_ice.  Chksums are computed if do_chksum = .True.
  subroutine coupler_flux_land_to_ice(Land, Ice, Land_ice_boundary, Time, current_timestep, &
                                      coupler_chksum_obj, coupler_clocks)

    implicit none
    type(land_data_type), intent(inout) :: Land !< Land
    type(ice_data_type),  intent(inout) :: Ice  !< Ice
    type(land_ice_boundary_type), intent(inout) :: Land_ice_boundary !< Land_ice_boundary
    type(FmsTime_type), intent(in) :: Time         !< Time (in seconds)
    integer, intent(in)       :: current_timestep  !< current timestep
    type(coupler_chksum_type), intent(in)    :: coupler_chksum_obj !< coupler_chksum_obj to compute chksums
    type(coupler_clock_type),  intent(inout) :: coupler_clocks      !< coupler_clocks

    call fms_mpp_clock_begin(coupler_clocks%flux_land_to_ice)
    call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )
    call fms_mpp_clock_end(coupler_clocks%flux_land_to_ice)

    if (do_chksum) call coupler_chksum_obj%get_atmos_ice_land_chksums('fluxlnd2ice+', current_timestep)

  end subroutine coupler_flux_land_to_ice

  !> This subroutine calls ice_model_fast_cleanup and unpack_land_ice_boundary
  subroutine coupler_unpack_land_ice_boundary(Ice, Land_ice_boundary, coupler_clocks)

    implicit none
    type(ice_data_type),          intent(inout) :: Ice               !< Ice
    type(land_ice_boundary_type), intent(inout) :: Land_ice_boundary !< Land_ice_boundary
    type(coupler_clock_type), intent(inout) :: coupler_clocks        !< coupler_clocks

    if (ice_npes .NE. atmos_npes) call fms_mpp_set_current_pelist(Ice%fast_pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_fast)

    !> These two calls occur on whichever PEs handle the fast ice processess.
    call ice_model_fast_cleanup(Ice)
    call unpack_land_ice_boundary(Ice, Land_ice_boundary)

    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_fast)

  end subroutine coupler_unpack_land_ice_boundary

  !> This subroutine calls update_ice_model_slow and flux_ice_to_ocean_stocks
  subroutine coupler_update_ice_model_slow_and_stocks(Ice, coupler_clocks)

    implicit none
    type(ice_data_type), intent(inout) :: Ice    !< Ice
    type(coupler_clock_type), intent(inout) :: coupler_clocks !< coupler_clocks

    if (slow_ice_with_ocean) call fms_mpp_set_current_pelist(Ice%slow_pelist)
    call fms_mpp_clock_begin(coupler_clocks%update_ice_model_slow_slow)

    call update_ice_model_slow(Ice)

    call fms_mpp_clock_begin(coupler_clocks%flux_ice_to_ocean_stocks)
    call flux_ice_to_ocean_stocks(Ice)
    call fms_mpp_clock_end(coupler_clocks%flux_ice_to_ocean_stocks)

    call fms_mpp_clock_end(coupler_clocks%update_ice_model_slow_slow)

  end subroutine coupler_update_ice_model_slow_and_stocks

  !> This subroutine calls update_ocean_model.  Chksums are computed if do_chksum = .True.
  subroutine coupler_update_ocean_model(Ocean, Ocean_state, Ice_ocean_boundary, &
                                        Time_ocean, Time_step_cpld, current_timestep, coupler_chksum_obj)

    implicit none
    type(ocean_public_type),         intent(inout) :: Ocean               !< Ocean
    type(ocean_state_type), pointer, intent(inout) :: Ocean_state         !< Ocean_state
    type(Ice_ocean_boundary_type),   intent(inout) :: Ice_ocean_boundary  !< Ice_ocean_boundary
    type(FmsTime_type), intent(inout) :: Time_ocean   !< Time_ocean
    type(FmsTime_type), intent(in) :: Time_step_cpld  !< total number of timesteps
    integer, intent(in) :: current_timestep           !< current timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj !< used for checksum computation

    call update_ocean_model(Ice_ocean_boundary, Ocean_state,  Ocean, Time_ocean, Time_step_cpld)
    if (do_chksum) call coupler_chksum_obj%get_ocean_chksums('update_ocean_model+', current_timestep)

  end subroutine coupler_update_ocean_model

  !> Thie subroutine calls component restarts and coupler_restart where the intermediate restart files
  !! is produced in the latter calls.  Time_restart is the next timestep where the intermediate restart
  !! file will be written out.  Time_restart_current records the current restart time.
  subroutine coupler_intermediate_restart(Atm, Ice, Ocean, Ocean_state, Ocn_bc_restart, Ice_bc_restart,&
                                          Time_current, Time_restart, Time_restart_current, Time_start)

    implicit none
    type(atmos_data_type),   intent(inout) :: Atm    !< Atm
    type(ice_data_type),     intent(inout) :: Ice    !< Ice
    type(ocean_public_type), intent(inout) ::  Ocean !< Ocean
    type(ocean_state_type),      pointer, intent(inout) :: Ocean_state       !< Ocean_state
    type(FmsNetcdfDomainFile_t), pointer, intent(inout) :: Ocn_bc_restart(:) !< used for coupler type restarts
    type(FmsNetcdfDomainFile_t), pointer, intent(inout) :: Ice_bc_restart(:) !< used for coupler type restarts
    type(FmsTime_type), intent(in) :: Time_current, Time_start  !< current Timestep and model start time
    !> Restart files will be written when Time=>Time_restart.  Time_restart is incremented by restart_interval
    !! Time_restart_current records the current timestep the restart file is being written.
    !! Time_restart_current does not necessary = Time_restart.
    type(FmsTime_type), intent(inout) :: Time_restart, Time_restart_current
    character(len=32) :: timestamp !< Time in string
    integer :: outunit             !< stdout

    Time_restart_current = Time_current

    timestamp = fms_time_manager_date_to_string(Time_restart_current)
    outunit= fms_mpp_stdout()
    write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
                      trim(timestamp),' is appended as prefix to each restart file name'
    if (Atm%pe) then
      call atmos_model_restart(Atm, timestamp)
      call land_model_restart(timestamp)
      call ice_model_restart(Ice, timestamp)
    endif
    if (Ocean%is_ocean_pe) call ocean_model_restart(Ocean_state, timestamp)

    call coupler_restart(Atm, Ice, Ocean, Ocn_bc_restart, Ice_bc_restart, &
                         Time_current, Time_restart_current, Time_start, timestamp)

    Time_restart = fms_time_manager_increment_date(Time_current, restart_interval(1), restart_interval(2), &
                   restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )

  end subroutine coupler_intermediate_restart

  !> This subroutine mainly prints out the current timestep in the stdout.
  !! Chksum is computed if do_chksum = .True.
  subroutine coupler_summarize_timestep(current_timestep, num_cpld_calls, coupler_chksum_obj, &
                                        is_atmos_pe, omp_sec, imb_sec)

    implicit none
    integer, intent(in) :: current_timestep  !< current_timestep, nc
    integer, intent(in) :: num_cpld_calls    !< total number of outerloop timestep
    type(coupler_chksum_type), intent(in) :: coupler_chksum_obj  !< coupler_chksum_obj
    logical, intent(in)               :: is_atmos_pe             !< Atm%pe
    real, dimension(:), intent(inout) :: omp_sec, imb_sec        !< from omp computation

    integer :: outunit        !< stdout
    character(len=80) :: text !< text to be written out to stdout

    if (do_chksum) call coupler_chksum_obj%get_coupler_chksums('MAIN_LOOP+', current_timestep)
    write( text,'(a,i6)' )'Main loop at coupling timestep=', current_timestep
    call fms_memutils_print_memuse_stats(text)
    outunit= fms_mpp_stdout()

    if (fms_mpp_pe() == fms_mpp_root_pe() .and. is_atmos_pe .and. do_concurrent_radiation) &
        write(outunit,102) 'At coupling step ', current_timestep,' of ',num_cpld_calls, ' Atm & Rad (imbalance): ', &
                            omp_sec(1),' (',imb_sec(1),')  ',omp_sec(2),' (',imb_sec(2),')'

    call flush(outunit)

102 format(A17,i5,A4,i5,A24,f10.4,A2,f10.4,A3,f10.4,A2,f10.4,A1)

  end subroutine coupler_summarize_timestep

end module full_coupler_mod
