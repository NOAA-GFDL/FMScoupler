! This include file contains the subroutine
! flux_exchange_init
!
! This subroutine makes extensive use of exchange grid global variables
! defined in the main flux_exchange module
!#######################################################################
! <SUBROUTINE NAME="flux_exchange_init">
!  <OVERVIEW>
!   Initialization routine.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initializes the interpolation routines,diagnostics and boundary data
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
!		atmos_ice_boundary, land_ice_atmos_boundary, &
!		land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )
!		
!  </TEMPLATE>
!  <IN NAME=" Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Atm" TYPE="atmos_data_type">
!   A derived data type to specify atmosphere boundary data.
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <IN NAME="Ocean" TYPE="ocean_data_type">
!   A derived data type to specify ocean boundary data.
!  </IN>
!  <INOUT NAME="atmos_ice_boundary" TYPE="atmos_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from atmosphere to ice.
!  </INOUT>
!  <INOUT NAME="land_ice_atmos_boundary" TYPE="land_ice_atmos_boundary_type">
!   A derived data type to specify properties and fluxes passed from exchange grid to
!   the atmosphere, land and ice.
!  </INOUT>
!  <INOUT NAME="land_ice_boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!  <INOUT NAME="ice_ocean_boundary" TYPE="ice_ocean_boundary_type">
!  A derived data type to specify properties and fluxes passed from ice to ocean.
!  </INOUT>
!  <INOUT NAME="ocean_ice_boundary" TYPE="ocean_ice_boundary_type">
!  A derived data type to specify properties and fluxes passed from ocean to ice.
!  </INOUT>
!
subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
       atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

  type(time_type),                   intent(in)  :: Time
  type(atmos_data_type),             intent(inout)  :: Atm
  type(land_data_type),              intent(in)  :: Land
  type(ice_data_type),               intent(inout)  :: Ice
  type(ocean_data_type),             intent(inout)  :: Ocean
! All intent(OUT) derived types with pointer components must be 
! COMPLETELY allocated here and in subroutines called from here;
! NO pointer components should have been allocated before entry if the
! derived type has intent(OUT) otherwise they may be lost.
  type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary
  type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary
  type(land_ice_boundary_type),      intent(inout) :: land_ice_boundary
  type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary
  type(ocean_ice_boundary_type),     intent(inout) :: ocean_ice_boundary

  character(len=64), parameter    :: sub_name = 'flux_exchange_init'
  character(len=256), parameter   :: error_header = '==>Error from ' // trim(module_name) //   &
                                                    '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: warn_header = '==>Warning from ' // trim(module_name) //  &
                                                   '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: note_header = '==>Note from ' // trim(module_name) //     &
                                                   '(' // trim(sub_name) // '):'
  
  integer :: unit, ierr, io,  i
  integer :: rcode, ncid, varid, dims(4), start(4), nread(4)
  integer :: nlon, nlat

  real, dimension(:), allocatable :: atmlonb, atmlatb
  integer :: is, ie, js, je, kd
  integer, dimension(:), allocatable :: pelist  
  character(32) :: tr_name
  logical       :: found

  integer       ::  n

!-----------------------------------------------------------------------

!
!       initialize atmos_ocean_fluxes
! Setting up flux types, allocates the arrays.
!

!
!       ocean_tracer_flux_init is called first since it has the meaningful value to set
!       for the input/output file names for the tracer flux values used in restarts. These
!       values could be set in the field table, and this ordering allows this.
!       atmos_tracer_flux_init is called last since it will use the values set in 
!       ocean_tracer_flux_init with the exception of atm_tr_index, which can only
!       be meaningfully set from the atmospheric model (not from the field table)
!

    call ocean_model_flux_init
    call atmos_tracer_flux_init
    call atmos_ocean_fluxes_init(ex_gas_fluxes, ex_gas_fields_atm, ex_gas_fields_ice)

!-----------------------------------------------------------------------
!----- read namelist -------

    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=flux_exchange_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'flux_exchange_nml')
    enddo
10  call mpp_close(unit)

!----- write namelist to logfile -----
    call write_version_number (version, tag)
    if( mpp_pe() == mpp_root_pe() )write( stdlog(), nml=flux_exchange_nml )

!----- find out number of atmospheric prognostic tracers and index of specific 
!      humidity in the tracer table
  call get_number_tracers (MODEL_ATMOS, num_tracers=n_atm_tr_tot, &
                           num_prog=n_atm_tr)
  call get_number_tracers (MODEL_LAND, num_tracers=n_lnd_tr_tot, &
                           num_prog=n_lnd_tr)

  ! assemble the table of tracer number translation by matching names of
  ! prognostic tracers in the atmosphere and surface models; skip all atmos.
  ! tracers that have no corresponding surface tracers.
  allocate(tr_table(n_atm_tr))
  allocate(tr_table_map(n_atm_tr))
  n = 1
  do i = 1,n_atm_tr
     call get_tracer_names( MODEL_ATMOS, i, tr_name )
     tr_table(n)%atm = i
     tr_table(n)%ice = get_tracer_index ( MODEL_ICE,  tr_name )
     tr_table_map(i)%ice = tr_table(n)%ice
     tr_table(n)%lnd = get_tracer_index ( MODEL_LAND, tr_name )
     tr_table_map(i)%lnd = tr_table(n)%lnd
     if(tr_table(n)%ice/=NO_TRACER.or.tr_table(n)%lnd/=NO_TRACER) then
       tr_table_map(i)%exch = n
       n = n + 1
     endif
  enddo
  n_exch_tr = n - 1
  !
  !     Set up tracer table entries for ocean-atm gas fluxes where the names of tracers in the
  !     atmosphere and ocean may not be equal
  !
  do n = 1, ex_gas_fluxes%num_bcs  !{
    if (ex_gas_fluxes%bc(n)%atm_tr_index .gt. 0) then  !{
      found = .false.
      do i = 1, n_exch_tr  !{
        if (ex_gas_fluxes%bc(n)%atm_tr_index .eq. tr_table(i)%atm) then
          found = .true.
          exit
        endif
      enddo  !} i
      if (.not. found) then
        n_exch_tr = n_exch_tr + 1
        tr_table(n_exch_tr)%atm = ex_gas_fluxes%bc(n)%atm_tr_index
        tr_table(n_exch_tr)%ice = NO_TRACER ! because ocean-atm gas fluxes are not held in the ice model as tracers
        tr_table(n_exch_tr)%lnd = NO_TRACER ! because this would have been found above
        tr_table_map(n_exch_tr)%exch = n_exch_tr
        tr_table_map(n_exch_tr)%ice = tr_table(n_exch_tr)%ice
        tr_table_map(n_exch_tr)%lnd = tr_table(n_exch_tr)%lnd
      endif
    endif  !}
  enddo  !} n
  write(stdout(),*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  write(stdlog(),*) trim(note_header), ' Number of exchanged tracers = ', n_exch_tr
  do i = 1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
     write(stdout(),*)'Tracer field name :'//trim(tr_name)
     write(stdlog(),*)'Tracer field name :'//trim(tr_name)
  enddo

  ! find out which tracer is specific humidity

  ! +fix-me-slm+ specific humidity may not be present if we are running with
  ! dry atmosphere. Besides, model may use mixing ratio ('mix_rat') (?). However,
  ! some atmos code also assumes 'sphum' is present, so for now the following
  ! code may be good enough.

  do i = 1,n_exch_tr
     call get_tracer_names( MODEL_ATMOS, tr_table(i)%atm, tr_name )
     if(lowercase(tr_name)=='sphum') then
        isphum = i
        exit
     endif
  enddo

  if (isphum==NO_TRACER) then
     call error_mesg('flux_exchange_mod',&
          'tracer "sphum" must be present in the atmosphere', FATAL )
  endif

!--------- read gridspec file ------------------
!only atmos pelists needs to do it here, ocean model will do it elsewhere

    ice_pe = Atm%pe
    ocn_pe = Ocean%pe
    allocate( ice_pelist(size(Atm%pelist)) ) !if ice/land become concurrent, this won't be true...
    ice_pelist(:) = Atm%pelist(:)
    allocate( ocn_pelist(size(Ocean%pelist)) )
    ocn_pelist(:) = Ocean%pelist(:)

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        rcode = nf_open('INPUT/grid_spec.nc',0,ncid)
!   <ERROR MSG="cannot open INPUT/grid_spec.nc" STATUS="FATAL">
!      Can not open the file INPUT/grid_spec.nc. The possible reason is 
!      that file grid_spec.nc is not in the directory INPUT.
!   </ERROR>
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot open INPUT/grid_spec.nc', FATAL)

!
! check atmosphere and grid_spec.nc have same atmosphere lat/lon boundaries
!
        rcode = nf_inq_varid(ncid, 'AREA_ATM', varid)
!   <ERROR MSG="cannot find AREA_ATM on INPUT/grid_spec.nc" STATUS="FATAL">
!      File INPUT/grid_spec.nc does not contain the field AREA_ATM.
!   </ERROR>
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot find AREA_ATM on INPUT/grid_spec.nc', &
             FATAL)
        rcode = nf_inq_vardimid(ncid, varid, dims)
        rcode = nf_inq_dimlen(ncid, dims(1), nlon)
        rcode = nf_inq_dimlen(ncid, dims(2), nlat)
        if (nlon+1/=size(Atm%glon_bnd(:)).or.nlat+1/=size(Atm%glat_bnd(:))) then
            if (mpp_pe()==mpp_root_pe()) then
                print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                     'atmosphere has', size(Atm%glon_bnd(:))-1, 'longitudes,', &
                     size(Atm%glat_bnd(:))-1, 'latitudes (see xba.dat and yba.dat)'
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere resolution" STATUS="FATAL">
!      The atmosphere grid size from file grid_spec.nc is not compatible with the atmosphere 
!      resolution from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod',  &
                 'grid_spec.nc incompatible with atmosphere resolution', FATAL)
        end if
        allocate( atmlonb(size(Atm%glon_bnd(:))) )
        allocate( atmlatb(size(Atm%glat_bnd(:))) )
        rcode = nf_inq_varid(ncid, 'xba', varid)
        start = 1; nread = 1; nread(1) = nlon+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlonb)
        rcode = nf_inq_varid(ncid, 'yba', varid)
        start = 1; nread = 1; nread(1) = nlat+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlatb)
        if (maxval(abs(atmlonb-Atm%glon_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == mpp_root_pe()) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY'
                do i=1,size(atmlonb(:))
                   print *,atmlonb(i),Atm%glon_bnd(i)*45/atan(1.0)
                end do
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)" STATUS="FATAL">
!      longitude from file grid_spec.nc ( from field xba ) is different from the longitude from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                 ,FATAL)
        end if
        if (maxval(abs(atmlatb-Atm%glat_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == mpp_root_pe()) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY'
                do i=1,size(atmlatb(:))
                   print *,atmlatb(i),Atm%glat_bnd(i)*45/atan(1.0)
                end do
            end if
!   <ERROR MSG="grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)" STATUS="FATAL">
!      longitude from file grid_spec.nc ( from field yba ) is different from the longitude from atmosphere model.
!   </ERROR>
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                 , FATAL)
        end if

        call xgrid_init(remap_method)

        call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
             (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
             "INPUT/grid_spec.nc")
        call generate_sfc_xgrid( Land, Ice )

        call setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
             (/ Land%Domain, Ice%Domain /),                    &
             "INPUT/grid_spec.nc"             )
        n_xgrid_runoff = max(xgrid_count(xmap_runoff),1)

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!----- initialize quantities for global integral package -----

!! call diag_integral_field_init ('prec', 'f6.3')
        call diag_integral_field_init ('evap', 'f6.3')

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

        call diag_field_init ( Time, Atm%axes(1:2), Land%axes )
        ni_atm = size(Atm%lon_bnd(:))-1 ! to dimension "diag_atm"
        nj_atm = size(Atm%lat_bnd(:))-1 ! in flux_ocean_to_ice

!Balaji
        
!allocate atmos_ice_boundary
        call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
        kd = size(Ice%ice_mask,3)
        allocate( atmos_ice_boundary%u_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%v_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%u_star(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%t_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%q_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dir(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux_vis_dif(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%fprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dhdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dedt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%drdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%coszen(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%p(is:ie,js:je,kd) )
! initialize boundary values for override experiments (mjh)
        atmos_ice_boundary%u_flux=0.0
        atmos_ice_boundary%v_flux=0.0
        atmos_ice_boundary%u_star=0.0
        atmos_ice_boundary%t_flux=0.0
        atmos_ice_boundary%q_flux=0.0
        atmos_ice_boundary%lw_flux=0.0
        atmos_ice_boundary%sw_flux=0.0
        atmos_ice_boundary%sw_flux_vis=0.0
        atmos_ice_boundary%sw_flux_dir=0.0
        atmos_ice_boundary%sw_flux_vis_dir=0.0
        atmos_ice_boundary%sw_flux_dif=0.0
        atmos_ice_boundary%sw_flux_vis_dif=0.0
        atmos_ice_boundary%lprec=0.0
        atmos_ice_boundary%fprec=0.0
        atmos_ice_boundary%dhdt=0.0
        atmos_ice_boundary%dedt=0.0
        atmos_ice_boundary%drdt=0.0
        atmos_ice_boundary%coszen=0.0
        atmos_ice_boundary%p=0.0

!         allocate fields for extra fluxes
! Copying initialized gas fluxes from exchange grid to atmosphere_ice boundary

        call coupler_type_copy(ex_gas_fluxes, atmos_ice_boundary%fluxes, is, ie, js, je, kd,    &
             mod_name, Ice%axes, Time, suffix = '_atm_ice')
  
!allocate land_ice_boundary
        allocate( land_ice_boundary%runoff(is:ie,js:je) )
        allocate( land_ice_boundary%calving(is:ie,js:je) )
! initialize values for override experiments (mjh)
        land_ice_boundary%runoff=0.0
        land_ice_boundary%calving=0.0
!allocate land_ice_atmos_boundary
        call mpp_get_compute_domain( Atm%domain, is, ie, js, je )
        allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dir(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_vis_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo_nir_dif(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_tr(is:ie,js:je,n_atm_tr) )
        allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudu(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%q_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
! initialize boundary values for override experiments (mjh)
        land_ice_atmos_boundary%t=273.0
        land_ice_atmos_boundary%albedo=0.0
        land_ice_atmos_boundary%albedo_vis_dir=0.0
        land_ice_atmos_boundary%albedo_nir_dir=0.0
        land_ice_atmos_boundary%albedo_vis_dif=0.0
        land_ice_atmos_boundary%albedo_nir_dif=0.0
        land_ice_atmos_boundary%land_frac=0.0
        land_ice_atmos_boundary%dt_t=0.0
        land_ice_atmos_boundary%dt_tr=0.0
        land_ice_atmos_boundary%u_flux=0.0
        land_ice_atmos_boundary%v_flux=0.0
        land_ice_atmos_boundary%dtaudu=0.0
        land_ice_atmos_boundary%dtaudv=0.0
        land_ice_atmos_boundary%u_star=0.0
        land_ice_atmos_boundary%b_star=0.0
        land_ice_atmos_boundary%q_star=0.0
        land_ice_atmos_boundary%rough_mom=0.01

! allocate fields for extra tracers
! The first call is no longer necessary, the fluxes will be passed by the land module
! The 2nd call is useful in the case of a ocean model only simulation
!
        call coupler_type_copy(ex_gas_fields_atm, Atm%fields, is, ie, js, je,                   &
             mod_name, Atm%axes(1:2), Time, suffix = '_atm')

!Balaji: clocks on atm%pe only        
    cplClock = mpp_clock_id( 'Land-ice-atm coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    sfcClock = mpp_clock_id( 'SFC boundary layer', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
    fluxAtmDnClock = mpp_clock_id( 'Flux DN from atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxLandIceClock = mpp_clock_id( 'Flux land to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    regenClock = mpp_clock_id( 'XGrid generation', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxAtmUpClock = mpp_clock_id( 'Flux UP to atm', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    end if

    allocate(pelist(size(Atm%pelist(:))+size(Ocean%pelist(:))))
    pelist(1:size(Atm%pelist(:))) = Atm%pelist(:)
    pelist(size(Atm%pelist(:))+1:size(pelist(:))) = Ocean%pelist(:)
    call mpp_set_current_pelist(pelist)
    deallocate(pelist)

!ocean_ice_boundary and ice_ocean_boundary must be done on all PES
!domain boundaries will assure no space is allocated on non-relevant PEs.
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
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

!
! allocate fields for extra tracers
! Copying gas flux fields from ice to ocean_ice boundary

    call coupler_type_copy(ex_gas_fields_ice, ocean_ice_boundary%fields, is, ie, js, je,        &
         mod_name, Ice%axes(1:2), Time, suffix = 'ocn_ice')

!allocate ice_ocean_boundary
    call mpp_get_compute_domain( Ocean%domain, is, ie, js, je )
!ML ocean only requires t, q, lw, sw, fprec, calving
!AMIP ocean needs no input fields
!choice of fields will eventually be done at runtime
!via field_manager
    allocate( ice_ocean_boundary%u_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%v_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%t_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%q_flux   (is:ie,js:je) )
    allocate( ice_ocean_boundary%salt_flux(is:ie,js:je) )
    allocate( ice_ocean_boundary%lw_flux  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis_dir  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_vis_dif  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_dir  (is:ie,js:je) )
    allocate( ice_ocean_boundary%sw_flux_dif  (is:ie,js:je) )
    allocate( ice_ocean_boundary%lprec    (is:ie,js:je) )
    allocate( ice_ocean_boundary%fprec    (is:ie,js:je) )
    allocate( ice_ocean_boundary%runoff   (is:ie,js:je) )
    allocate( ice_ocean_boundary%calving  (is:ie,js:je) )
    allocate( ice_ocean_boundary%p        (is:ie,js:je) )

!
! allocate fields for extra tracers
!

    call coupler_type_copy(ex_gas_fluxes, ice_ocean_boundary%fluxes, is, ie, js, je,    &
         mod_name, Ocean%axes(1:2), Time, suffix = '_ice_ocn')

    call coupler_type_copy(ex_gas_fields_ice, Ocean%fields, is, ie, js, je,                &
         mod_name, Ocean%axes(1:2), Time, suffix = '_ocn')

!pjp Why are the above not initialized to zero?
! initialize boundary values for override experiments
    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%domain )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype


!
! allocate fields amd fluxes for extra tracers for the Ice type
!

    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
    kd = size(Ice%ice_mask,3)
    call coupler_type_copy(ex_gas_fields_ice, Ice%ocean_fields, is, ie, js, je, kd,     &
         mod_name, Ice%axes, Time, suffix = '_ice')

    call coupler_type_copy(ex_gas_fluxes, Ice%ocean_fluxes, is, ie, js, je,             &
         mod_name, Ice%axes(1:2), Time, suffix = '_ice')

    call coupler_type_copy(ex_gas_fluxes, Ice%ocean_fluxes_top, is, ie, js, je, kd,     &
         mod_name, Ice%axes, Time, suffix = '_ice_top')

!       initialize the Ocean type for extra fields for surface fluxes
! Same allocation of arrays and stuff
!       (this must be done after the Ocean fields are allocated as the fields on the Ocean%fields
!       are read in in this subroutine)
!

    if ( Ocean%pe ) then
      call mpp_set_current_pelist(Ocean%pelist)
      call ocean_model_init_sfc(Ocean)
    end if
    call mpp_set_current_pelist()

!BalaCji
    cplOcnClock = mpp_clock_id( 'Ice-ocean coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
    fluxIceOceanClock = mpp_clock_id( 'Flux ice to ocean', flags=clock_flag_default, grain=CLOCK_ROUTINE )
    fluxOceanIceClock = mpp_clock_id( 'Flux ocean to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )
!---- done ----
    do_init = .false.

  end subroutine flux_exchange_init
! </SUBROUTINE>
