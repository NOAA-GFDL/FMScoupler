!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
! utilities used by coupler
!
use time_manager_mod, only: time_type, set_calendar_type, set_time,  &
                            set_date, get_date, days_in_month, month_name,  &
                            operator(+), operator(-), operator (<), &
                            operator (>), operator (/=), operator (/), &
                            operator (*), thirty_day_months, julian, &
                            no_leap, no_calendar

use  utilities_mod, only: open_file, file_exist, check_nml_error,  &
                          error_mesg, FATAL, WARNING,              &
                          print_version_number, get_my_pe,         &
                          utilities_init, utilities_end,           &
                          close_file, check_system_clock

use  diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date
!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!
use  atmos_coupled_mod, only: atmos_coupled_init, atmos_coupled_end, &
                            update_atmos_coupled_down,           &
                            update_atmos_coupled_up,             &
                            atmos_boundary_data_type

use   land_model_mod, only: land_model_init, land_model_end, &
                            land_boundary_data_type, &
                            update_land_model_fast, update_land_model_slow

use    ice_model_mod, only: ice_model_init, ice_model_end,  &
                            ice_bottom_to_ice_top,          &
                            update_ice_model_fast,          &
                            update_ice_model_slow,          &
                            ice_data_type

use  ocean_model_mod, only: update_ocean_model, ocean_model_init,  &
                            ocean_model_end, ocean_data_type
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!
use flux_exchange_mod, only: flux_exchange_init,   &
                             flux_calculation,     &
                             flux_down_from_atmos, &
                             flux_up_to_atmos,     &
                             flux_land_to_ice,     &
                             flux_ice_to_ocean,    &
                             flux_ocean_to_ice

implicit none

!-----------------------------------------------------------------------

 character(len=4), parameter :: vers_num = 'v2.0'
character(len=128) :: version = '$Id: coupler_main.F90,v 1.4 2001/10/25 17:51:58 fms Exp $'
character(len=128) :: tag = '$Name: fez $'

!-----------------------------------------------------------------------
!---- model defined-types ----

 type (atmos_boundary_data_type) :: Atm
 type  (land_boundary_data_type) :: Land
 type   (ice_data_type)          :: Ice
 type (ocean_data_type)          :: Ocean
 
!-----------------------------------------------------------------------
!---- storage for fluxes ----

 real, allocatable, dimension(:,:)   ::                         &
    t_surf_atm, albedo_atm, land_frac_atm, dt_t_atm, dt_q_atm,  &
    flux_u_atm, flux_v_atm, dtaudv_atm, u_star_atm, b_star_atm, &
    rough_mom_atm

 real, allocatable, dimension(:,:)   :: flux_u_ocean,  flux_v_ocean, &
                                        flux_t_ocean,  flux_q_ocean, &
                                       flux_lw_ocean, flux_sw_ocean, &
                                         lprec_ocean,   fprec_ocean, &
                                        runoff_ocean, calving_ocean, &
                                     flux_salt_ocean,  p_surf_ocean

 real, allocatable, dimension(:,:,:) ::                                &
    flux_t_land, flux_q_land, flux_lw_land, flux_sw_land,              &
    dhdt_land  , dedt_land  , drdt_land   , lprec_land  , fprec_land

 real, allocatable, dimension(:,:) :: t_surf_ice, s_surf_ice, frazil_ice, &
                                      sea_lev_ice, u_surf_ice, v_surf_ice
 real, allocatable, dimension(:,:,:) ::                                &
    flux_t_ice , flux_q_ice , flux_lw_ice, flux_sw_ice , coszen_ice,   &
    dhdt_ice   , dedt_ice   , drdt_ice   , lprec_ice   , fprec_ice   , &
    flux_u_ice,  flux_v_ice, p_surf_ice

 real, allocatable, dimension(:,:) :: runoff_ice, calving_ice

!-----------------------------------------------------------------------
! ----- coupled model time -----

   type (time_type) :: Time, Time_init, Time_end,  &
                       Time_step_ocean, Time_step_atmos,  &
		       Time_step_cpld
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

      namelist /coupler_nml/ current_date, calendar, override,       &
                             months, days, hours, minutes, seconds,  &
                             dt_cpld, dt_atmos, dt_ocean

!#######################################################################

 call    utilities_init ( )

 call initialize_coupler 

 call check_system_clock ('END OF INITIALIZATION')
!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop -------

 do nc = 1, num_cpld_calls


    call flux_ocean_to_ice (Ocean, Ice,   t_surf_ice,                         &
                                          u_surf_ice,     v_surf_ice,         &
                                          frazil_ice, s_surf_ice, sea_lev_ice )

    call ice_bottom_to_ice_top (Ice, t_surf_ice,                 &
                                     u_surf_ice,     v_surf_ice, &
                                     frazil_ice, s_surf_ice, sea_lev_ice )

!-----------------------------------------------------------------------
!   ------ atmos/fast-land/fast-ice integration loop -------

    do na = 1, num_atmos_calls

       Time = Time + Time_step_atmos

       call flux_calculation  (float(dt_atmos), Time,                  &
                               Atm, Land, Ice, land_frac_atm,          &
                               t_surf_atm, albedo_atm, rough_mom_atm,  &
                               flux_u_atm, flux_v_atm, dtaudv_atm,     &
                               u_star_atm, b_star_atm                  )

!      ---- atmosphere down ----

       call update_atmos_coupled_down  (Atm,    &
                              t_surf_atm, albedo_atm, rough_mom_atm,  &
                              u_star_atm, b_star_atm, land_frac_atm,  &
                              dtaudv_atm, flux_u_atm, flux_v_atm      )


       call flux_down_from_atmos (Time, Atm, Land, Ice,              &
                                  flux_u_atm,   flux_v_atm,          &
                                  flux_t_land,  flux_q_land,         &
                                  flux_lw_land, flux_sw_land,        &
                                  dhdt_land, dedt_land, drdt_land,   &
                                  lprec_land,  fprec_land,           &
                                  flux_t_ice,  flux_q_ice,           &
                                  flux_lw_ice, flux_sw_ice,          &
                                  dhdt_ice, dedt_ice, drdt_ice,      &
                                  lprec_ice ,  fprec_ice,            &
                                  flux_u_ice,  flux_v_ice, coszen_ice)


!      --------------------------------------------------------------

!      ---- land model ----

       call update_land_model_fast (Land, flux_sw_land, flux_lw_land,    &
                                    flux_t_land,  flux_q_land,           &
                                    dhdt_land,   dedt_land,   drdt_land, &
                                    lprec_land, fprec_land               )

!      ---- ice model ----


       call update_ice_model_fast (Ice, flux_u_ice,  flux_v_ice,      &
                                        flux_sw_ice, flux_lw_ice,     &
                                        flux_t_ice,  flux_q_ice,      &
                                  dhdt_ice,   dedt_ice,   drdt_ice,   &
                                  lprec_ice , fprec_ice, coszen_ice   )

!      --------------------------------------------------------------
!      ---- atmosphere up ----

       call flux_up_to_atmos (Time, Land, Ice, dt_t_atm, dt_q_atm)

       call update_atmos_coupled_up   (Atm, land_frac_atm,    &
                                           dt_t_atm, dt_q_atm)

!--------------

    enddo

!   ------ end of atmospheric time step loop -----
      call update_land_model_slow(Land)
!-----------------------------------------------------------------------

!
!     need flux call to put runoff and p_surf on ice grid
!
      call flux_land_to_ice(Land, Ice, runoff_ice, calving_ice)

      p_surf_ice  = 0.0 ! call flux_atmos_to_ice_slow ?

!   ------ slow-ice model ------

      call update_ice_model_slow (Ice, runoff_ice, calving_ice, p_surf_ice)

      call flux_ice_to_ocean ( Ice, flux_u_ocean,  flux_v_ocean, &
                                    flux_t_ocean,  flux_q_ocean, &
                                   flux_sw_ocean, flux_lw_ocean, &
                                     lprec_ocean,   fprec_ocean, &
                                    runoff_ocean, calving_ocean, &
                                 flux_salt_ocean,  p_surf_ocean  )

      do no = 1,num_ocean_calls

        ocean_seg_start = ( no .eq. 1 )               ! could eliminate these by
        ocean_seg_end   = ( no .eq. num_ocean_calls ) ! putting this loop in
                                                      ! update_ocean_model since
                                                      ! fluxes don't change here

        call update_ocean_model (Ocean, flux_u_ocean,   flux_v_ocean,   &
                                        flux_t_ocean,   flux_q_ocean,   &
                                        flux_sw_ocean, flux_lw_ocean,   &
                                        lprec_ocean,     fprec_ocean,   &
                                        runoff_ocean,  calving_ocean,   &
                                     flux_salt_ocean,   p_surf_ocean,   &
                                        ocean_seg_start, ocean_seg_end, &
                                        num_ocean_calls )

      enddo
!--------------




!   ------ end of ocean time step loop -----
!-----------------------------------------------------------------------

!--------------

 enddo

!-----------------------------------------------------------------------
 call check_system_clock ('END OF TIME LOOP')

 call coupler_end

 call diag_manager_end (Time)
 call utilities_end

!-----------------------------------------------------------------------

 stop

contains

!#######################################################################

   subroutine initialize_coupler

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
    integer :: unit, log_unit, ierr, io, id, jd, kd, m, i
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    logical :: use_namelist
!-----------------------------------------------------------------------
!----- read namelist -------

   unit = open_file ('input.nml', action='read')
   ierr=1; do while (ierr /= 0)
          read  (unit, nml=coupler_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'coupler_nml')
   enddo
10 call close_file (unit)

!----- write namelist to logfile (close log_unit later) -----

   log_unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) then
        write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (log_unit, nml=coupler_nml)
   endif

!----- read restart file -----

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
    else if (calendar(1:7) == 'no_leap') then
        calendar_type = no_leap
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

 16 format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

      call diag_manager_init

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
     'ocean time step is not a multiple of the cpld time step', 2)
!-----------------------------------------------------------------------
!------ initialize component models ------
!------ use global grids -------

!---- atmosphere ----

      call  atmos_coupled_init (Atm, Time_init, Time, Time_step_atmos)

!
! write atmos lon/lat bounds for make_xgrids input
!
!     if (get_my_pe()==0) then
!       print *, 'atmosphere lon boundaries:'
!       do i=1,size(Atm%glon_bnd)
!         print *, Atm%glon_bnd(i)*45/atan(1.0)
!       end do
!       print *, 'atmosphere lat boundaries:'
!       do i=1,size(Atm%glat_bnd)
!         print *, Atm%glat_bnd(i)*45/atan(1.0)
!       end do
!     end if

!---- ocean (with same grid as atmosphere) ----
!---- mask is read from restart -----

      call ocean_model_init (Ocean, Time_init, Time, Time_step_ocean, &
                                                         Atm%Domain )

!---- initialize land model ------

      call land_model_init (Land, Time_init, Time, Time_step_atmos, &
                                                   Time_step_cpld, &
                                                   atmos_domain=Atm%Domain)

!---- initialize ice model -----

      call ice_model_init (Ice, Time_init, Time, Time_step_atmos, &
                                                 Time_step_cpld, Ocean%Domain)

!-----------------------------------------------------------------------
!---- setup allocatable storage for fluxes exchanged between models ----
!---- use local grids -----

   id = size(Atm%lon_bnd)-1
   jd = size(Atm%lat_bnd)-1
   allocate  (  t_surf_atm (id,jd),  albedo_atm (id,jd), &
             land_frac_atm (id,jd),  dtaudv_atm (id,jd), &
                  dt_t_atm (id,jd),    dt_q_atm (id,jd), &
                flux_u_atm (id,jd),  flux_v_atm (id,jd), &
                u_star_atm (id,jd),  b_star_atm (id,jd), &
             rough_mom_atm (id,jd)                       )

   id = size(Ocean%t_surf,1)
   jd = size(Ocean%t_surf,2)
   allocate  (  flux_u_ocean (id,jd),  flux_v_ocean (id,jd), &
                flux_t_ocean (id,jd),  flux_q_ocean (id,jd), &
               flux_sw_ocean (id,jd), flux_lw_ocean (id,jd), &
                 lprec_ocean (id,jd),   fprec_ocean (id,jd), &
                runoff_ocean (id,jd), calving_ocean (id,jd), &
              flux_salt_ocean(id,jd),  p_surf_ocean (id,jd)  )

   id = size(Land%lon_bnd)-1
   jd = size(Land%lat_bnd)-1
   kd = size(Land%mask,3)
   allocate  (  flux_t_land (id,jd,kd),     dhdt_land (id,jd,kd), &
                flux_q_land (id,jd,kd),     dedt_land (id,jd,kd), &
                                            drdt_land (id,jd,kd), &
               flux_sw_land (id,jd,kd),  flux_lw_land (id,jd,kd), &
                 lprec_land (id,jd,kd),    fprec_land (id,jd,kd)  )

   id = size(Ice%t_surf,1)
   jd = size(Ice%t_surf,2)
   kd = size(Ice%ice_mask,3)
   allocate  (  t_surf_ice (id,jd), s_surf_ice (id,jd), frazil_ice (id,jd), &
                sea_lev_ice (id,jd), u_surf_ice (id,jd), v_surf_ice (id,jd) )
   allocate  (  flux_t_ice (id,jd,kd),        dhdt_ice (id,jd,kd), &
                flux_q_ice (id,jd,kd),        dedt_ice (id,jd,kd), &
                                              drdt_ice (id,jd,kd), &
               flux_sw_ice (id,jd,kd),     flux_lw_ice (id,jd,kd), &
                 lprec_ice (id,jd,kd),       fprec_ice (id,jd,kd), &
                coszen_ice (id,jd,kd),      flux_u_ice (id,jd,kd), &
                flux_v_ice (id,jd,kd),      p_surf_ice (id,jd,kd)  )
   allocate ( runoff_ice(id,jd), calving_ice(id,jd) )

!-----------------------------------------------------------------------
!---- initialize flux exchange module ----

      call flux_exchange_init ( Time, Atm, Land, Ice, Ocean )

!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----

      unit = open_file ('RESTART/coupler.res',  &
                        form='native', action='write')
      call close_file (unit, status='delete')

!-----------------------------------------------------------------------

   end subroutine initialize_coupler

!#######################################################################

   subroutine coupler_end

   integer :: unit, date(6)
!-----------------------------------------------------------------------

      call atmos_coupled_end (Atm)
      call ocean_model_end (Ocean)
      call  land_model_end (Land)
      call   ice_model_end (Ice)

!----- compute current date ------

      call get_date (Time, date(1), date(2), date(3),  &
                           date(4), date(5), date(6))

!----- check time versus expected ending time ----

      if (Time /= Time_end) call error_mesg ('program coupler',  &
              'final time does not match expected ending time', WARNING)

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

end program coupler_main

