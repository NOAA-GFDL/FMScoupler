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
!

module ice_model_mod

!! Components
use   ice_albedo_mod, only:  ice_albedo_init, ice_albedo
use ocean_albedo_mod, only:  compute_ocean_albedo_new
use  ocean_rough_mod, only:  compute_ocean_roughness, fixed_ocean_roughness

!! FMS
use FMS, sst_anom_fms=>sst_anom
use FMSconstants, only: HLV, HLF, TFREEZE, pi

implicit none
private
real :: cmin,cmax

public :: ice_data_type, atmos_ice_boundary_type,               &
          ice_model_init, ice_model_end, update_ice_model_fast, &
          update_ice_model_slow
!----------------------------------------------------------------
!--------- NAMELIST INTERFACE ---------
!
!  use_climo_ice           = use monthly climatological amip ice mask
!  use_annual_ice          = use annual climatology amip ice mask
!  temp_ice_freeze         = temperature at which sea ice melts
!  heat_capacity_ocean     = heat capacity for ocean in simple mixed layer model
!
real    :: diff                     = 2.092
real    :: thickness_min            = 0.10
real    :: specified_ice_thickness  = 2.0
real    :: heat_capacity_ocean      = 1.e07
real    :: temp_ice_freeze          = -1.66
real    :: roughness_ice            = 1.e-4
logical :: mixed_layer_ocean        = .false.
logical :: use_climo_ice            = .false.
logical :: use_annual_ice           = .false.
logical :: use_climo_sst            = .false.
logical :: use_annual_sst           = .false.
character(len=64) :: ice_method = 'prognostic' ! none, uniform, or prognostic
character(len=64) :: sst_method = 'specified'  ! specified, uniform, or mixed_layer
                                               ! Additional sst specifications: 'aqua_planet_#' test cases are derived 
                                               ! from the 2000 paper by Neale and Hoskins, 'A standard test for AGCMs including 
                                               ! their physical parameterizations: I. The proposal, Atmospheric Science Letters'.   
                                               ! The 'aqua_planet_1' testcase corresponds to the 'Control' SST test case and 
                                               ! provides the pattern which is shifted for the subsequent cases.
                                               ! The test cases Control, and aqua_planet_5N-aqua_planet_60N were documented and used 
                                               ! in Burnett et al., 2021, GRL,  https://doi.org/10.1029/2020GL091980
                                               !   aqua_planet_1   = Control profile
                                               !   aqua_planet_2   = Peaked
                                               !   aqua_planet_3   = Flat
                                               !   aqua_planet_4   = Qobs
                                               !   aqua_planet_5   = Control shifted by 5N
                                               !   aqua_planet_6   = 1KEQ
                                               !   aqua_planet_7   = 3KEQ
                                               !   aqua_planet_8   = 3KW1
                                               !   aqua_planet_10N = Control shifted by 10N
                                               !   aqua_planet_15N = Control shifted by 15N
                                               !   aqua_planet_20N = Control shifted by 20N
                                               !   aqua_planet_25N = Control shifted by 25N
                                               !   aqua_planet_30N = Control shifted by 20N
                                               !   aqua_planet_35N = Control shifted by 35N
                                               !   aqua_planet_40N = Control shifted by 30N
                                               !   aqua_planet_45N = Control shifted by 45N
                                               !   aqua_planet_50N = Control shifted by 50N
                                               !   aqua_planet_55N = Control shifted by 55N
                                               !   aqua_planet_60N = Control shifted by 60N
                                               !   aqua_planet_65N = Control shifted by 65N
                                               !   aqua_planet_70N = Control shifted by 70N
                                               !   aqua_planet_75N = Control shifted by 75N
                                               !   aqua_planet_80N = Control shifted by 80N
                                               !   aqua_planet_85N = Control shifted by 85N
                                               !   aqua_planet_90N = Control shifted by 90N
real              :: temp_ice = 270.      ! used when ice_method = 'uniform'
real              :: temp_sst = 280.      ! used when sst_method = 'uniform'
real              :: sst_anom = 0.        ! sst perturbation used for sensitivity experiments
character(len=64) :: interp_method  = "bilinear" ! conservative or bilinear
logical :: do_netcdf_restart = .true.

namelist /ice_model_nml/ diff, thickness_min, specified_ice_thickness,        &
                         heat_capacity_ocean, temp_ice_freeze, roughness_ice, &
                         ice_method, use_climo_ice, use_annual_ice, temp_ice, &
                         sst_method, use_climo_sst, use_annual_sst, temp_sst, &
                         interp_method, do_netcdf_restart, sst_anom

!----------------------------------------------------------------

type ice_data_type
  type(domain2d),pointer                :: Domain

   real,    pointer, dimension(:,:)     :: glon_bnd =>NULL(), &
                                           glat_bnd =>NULL(), &
                                           lon_bnd =>NULL() , &
                                           lat_bnd =>NULL()

   real,    pointer, dimension(:,:)     :: glon =>NULL(), &
                                           glat =>NULL(), &
                                           lon =>NULL(), &
                                           lat =>NULL()

   logical, pointer, dimension(:,:)     :: gmask =>NULL(), &
                                           mask =>NULL(), &
                                           ice_mask =>NULL()

   real,    pointer, dimension(:,:)   :: t_surf =>NULL(), &
                                         albedo =>NULL(), &
                                         albedo_vis_dir =>NULL(), &
                                         albedo_nir_dir =>NULL(), &
                                         albedo_vis_dif =>NULL(), &
                                         albedo_nir_dif =>NULL(), &
                                         rough_mom =>NULL(),&
                                         rough_heat =>NULL(), &
                                         rough_moist =>NULL(), &
                                         thickness =>NULL()

   type (time_type)                   :: Time_Init, Time,  &
                                         Time_step_fast,   &
                                         Time_step_slow
end type ice_data_type

!----------------------------------------------------------------

type :: atmos_ice_boundary_type
  real, dimension(:,:), pointer :: u_star =>NULL(), &
                                   t_flux =>NULL(), &
                                   q_flux =>NULL(), &
                                   lw_flux =>NULL(), &
                                   sw_flux =>NULL(), &
                                   lprec =>NULL(), &
                                   fprec =>NULL()
  real, dimension(:,:), pointer :: dhdt =>NULL(), &
                                   dedt =>NULL(), &
                                   drdt =>NULL(), &
                                   coszen =>NULL(), &
                                   data =>NULL()
  integer :: xtype
end type atmos_ice_boundary_type

!----------------------------------------------------------------

integer :: is, ie, js, je
type(amip_interp_type), save :: Amip_ice, Amip_sst
logical :: module_is_initialized = .false.
character(len=64) :: fname = 'INPUT/ice_model.res.nc'

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

real, parameter :: LATENT = HLV + HLF

contains

!######################################################################

 subroutine update_ice_model_fast( Atmos_boundary, Ice )
 type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
 type (ice_data_type),          intent(inout) :: Ice

 real, dimension(is:ie,js:je) :: ts_new, gamma, flux_i, t_dt_surf, &
                 flux_t_new, flux_q_new, flux_lw_new, flux_sw_new, &
                 flux_u_new, flux_v_new, lprec_new,   fprec_new,   &
                 deriv

 logical, dimension(is:ie,js:je,2) :: mask_ice
 real,    dimension(is:ie,js:je,2) :: thickness_ice, t_surf_ice, albedo
 real,    dimension(is:ie,js:je)   :: albedo_vis_dir, albedo_nir_dir, &
                                      albedo_vis_dif, albedo_nir_dif
 integer :: dt

!-----------------------------------------------------------------------
!
!   updates ice model on the atmospheric (fast) time step
!   averages input quantities to be seen by the ocean
!
!    flux_u  = zonal wind stress
!    flux_v  = meridional wind stress
!    flux_sw = net shortwave radiation (down-up)
!    flux_sw_vis = net visible shortwave radiation (down-up)
!    flux_sw_dir = net direct shortwave radiation (down-up)
!    flux_sw_dif = net diffuse shortwave radiation (down-up)
!    flux_sw_vis_dir = net visible direct shortwave radiation (down-up)
!    flux_sw_vis_dif = net visible diffuse shortwave radiation (down-up)
!    flux_lw = net longwave radiation (down-up)
!    flux_t  = sensible heat flux
!    flux_q  = specific humidity flux
!    lprec   = liquid precipitiation rate (kg/m2/s)
!    fprec   = frozen precipitiation rate (kg/m2/s)
!    coszen  = cosine of the zenith angle
!
!-----------------------------------------------------------------------
!----- implicit update of ice surface temperature -----

if (trim(ice_method) == 'prognostic') then

    where (Ice%ice_mask)
      gamma = diff / max(Ice%thickness,thickness_min)
      flux_i = gamma * (TFREEZE + temp_ice_freeze - Ice%t_surf)

      t_dt_surf = (flux_i + Atmos_boundary%lw_flux  + Atmos_boundary%sw_flux - &
                            Atmos_boundary%t_flux - Atmos_boundary%q_flux*LATENT)  &
                            / (Atmos_boundary%dhdt + Atmos_boundary%dedt*LATENT +  &
                               Atmos_boundary%drdt + gamma)
      ts_new = Ice%t_surf + t_dt_surf
    elsewhere
      ts_new = 0.
    endwhere

  ! update sea ice temperature
  ! new temperature cannot be greater than freezing temperature
    where (Ice%ice_mask .and. ts_new > TFREEZE )
      Ice%t_surf  = TFREEZE
    endwhere

    where (Ice%ice_mask .and. ts_new <= TFREEZE)
      Ice%t_surf  = ts_new
    endwhere

endif


if (trim(ice_method) == 'uniform') then
    where (Ice%ice_mask)
      !t_dt_surf = temp_ice - Ice%t_surf
       Ice%t_surf = temp_ice
    endwhere
endif

!-----------------------------------------------------------------------
!----- simple mixed layer ocean ------

   if (trim(sst_method) == 'mixed_layer') then

       call get_time ( Ice%Time_step_slow, dt )

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux -   &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean
          t_dt_surf = flux_i/(1.0 -deriv)
          ts_new = Ice%t_surf + t_dt_surf
       endwhere

     ! update sea surface temperature
     ! note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif

   if (trim(sst_method) == 'uniform') then
       where (Ice%mask .and. .not.Ice%ice_mask)
          Ice%t_surf = temp_sst
       endwhere
   endif

!-----------------------------------------------------------------------
!------ update ocean/ice surface parameters --------

   !---- over all ocean points (regardless of ice) ----

    call compute_ocean_roughness ( Ice%mask,        &
                        Atmos_boundary%u_star(:,:), &
                              Ice%rough_mom  (:,:), &
                              Ice%rough_heat (:,:), &
                              Ice%rough_moist(:,:) )

    call compute_ocean_albedo_new ( Ice%mask, Atmos_boundary%coszen(:,:),  &
               albedo_vis_dir, albedo_vis_dif, albedo_nir_dir, albedo_nir_dif, Ice%lat )

   !---- over only ice points -----

    where (Ice%ice_mask)
      Ice%rough_mom   = roughness_ice
      Ice%rough_heat  = roughness_ice
      Ice%rough_moist = roughness_ice
    endwhere

    mask_ice     (:,:,2) = Ice%ice_mask
    thickness_ice(:,:,2) = Ice%thickness
    t_surf_ice   (:,:,2) = Ice%t_surf
    albedo       (:,:,2) = Ice%albedo
    albedo       (:,:,1) = Ice%albedo
    call ice_albedo (mask_ice, thickness_ice, t_surf_ice, albedo)

    where (Ice%ice_mask)
       Ice%albedo = albedo(:,:,2)
       Ice%albedo_vis_dir = albedo(:,:,2)
       Ice%albedo_nir_dir = albedo(:,:,2)
       Ice%albedo_vis_dif = albedo(:,:,2)
       Ice%albedo_nir_dif = albedo(:,:,2)
    elsewhere
       Ice%albedo = albedo(:,:,1)
       Ice%albedo_vis_dir = albedo_vis_dir
       Ice%albedo_nir_dir = albedo_vis_dir
       Ice%albedo_vis_dif = albedo_vis_dir
       Ice%albedo_nir_dif = albedo_vis_dir
    endwhere

!-----------------------------------------------------------------------
!--------- advance time -----------------

 Ice%Time = Ice%Time + Ice%Time_step_fast

!-----------------------------------------------------------------------

 end subroutine update_ice_model_fast

!######################################################################

 subroutine update_ice_model_slow( Atmos_boundary, Ice )
 type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
 type(ice_data_type),           intent(inout) :: Ice

 !---- get the specified sea-ice fraction -----

      if (trim(ice_method) == 'prognostic' .or. trim(ice_method) == 'uniform') then

          call prognostic_ice ( Ice )

      endif

 !---- get the specified ocean temperature -----

      if (trim(sst_method) == 'specified') then

          call prognostic_sst ( Ice )

      endif


 end subroutine update_ice_model_slow

!######################################################################

 subroutine prognostic_ice ( Ice )
 type(ice_data_type),           intent(inout) :: Ice

 real, dimension(is:ie, js:je) :: ice_frac

 !---- get the specified sea-ice fraction -----

    call get_amip_ice (Ice%Time, Amip_ice, ice_frac )

  ! determine which grid boxes have ice coverage
    where ( Ice%mask(:,:) .and. ice_frac > 0.5 )
        Ice%thickness(:,:) = specified_ice_thickness
        Ice%ice_mask (:,:) = .true.
        Ice%t_surf   (:,:) = MIN( Ice%t_surf(:,:), TFREEZE )
    elsewhere
        Ice%thickness(:,:) = 0.0
        Ice%ice_mask (:,:) = .false.
        Ice%t_surf   (:,:) = MAX( Ice%t_surf(:,:), TFREEZE )
    endwhere

 end subroutine prognostic_ice

!######################################################################

 subroutine prognostic_sst ( Ice )
 type(ice_data_type),           intent(inout) :: Ice

 real, dimension(is:ie, js:je) :: sea_temp

 !---- get the specified ocean temperature -----

    call get_amip_sst (Ice%Time, Amip_sst, sea_temp )

  ! determine which grid boxes have open ocean coverage ----
    where ( Ice%mask(:,:) .and. .not.Ice%ice_mask(:,:))
        Ice%t_surf   (:,:) = sea_temp
    endwhere

 end subroutine prognostic_sst

!######################################################################

! subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
! type(ocean_ice_boundary_type), intent(in)    :: Ocean_boundary
! type(ice_data_type),           intent(inout) :: Ice

! set temp at open ocean points

! end subroutine update_ice_model_slow_up

!######################################################################

 subroutine ice_model_init ( Ice, Time_Init, Time, &
                             Time_step_fast, Time_step_slow, &
                             glon_bnd, glat_bnd, Atmos_domain )
 type(ice_data_type), intent(inout) :: Ice
 type(time_type)    , intent(in)    :: Time_Init, Time, &
                                       Time_step_fast, Time_step_slow
 real               , intent(in)    :: glon_bnd(:,:), glat_bnd(:,:)
 type(domain2d), intent(in), target :: Atmos_domain

real :: lon0, lond, latd, amp, t_control, dellon, dom_wid, siggy, tempi
 integer :: isg, ieg, jsg, jeg
 integer :: unit, ierr, io, i, j
 integer :: ndim, nvar, natt, ntime, nlon, nlat, mlon, mlat, layout(2)
 logical :: need_ic
 type(FmsNetcdfDomainFile_t) :: land_mask_fileobj !< Land mask domain decomposed fileobj
 type(FmsNetcdfDomainFile_t) :: ice_restart_fileobj !< Ice restart domain decomposed fileobj

 if (module_is_initialized) then
     return
 endif

 !< Read the namelist
 read (input_nml_file, nml=ice_model_nml, iostat=io)
 ierr = check_nml_error(io, 'ice_model_nml')

 do_netcdf_restart = .true. !< Always do netcdf!

 call write_version_number (version, tagname)
 if ( mpp_pe() == mpp_root_pe() ) then
    write (stdlog(), nml=ice_model_nml)
 endif

!---- error checks ----

  if ( trim(ice_method) /= 'none'    .and. &
       trim(ice_method) /= 'uniform' .and. &
       trim(ice_method) /= 'prognostic' ) call error_mesg &
     ('ice_model_init', 'namelist variable ice_method has invalid value', FATAL)

  if ( trim(sst_method) /= 'specified'            .and. &
       trim(sst_method) /= 'uniform'              .and. &
       trim(sst_method) /= 'aqua_planet_1'        .and. &
       trim(sst_method) /= 'aqua_planet_2'        .and. &
       trim(sst_method) /= 'aqua_planet_3'        .and. &
       trim(sst_method) /= 'aqua_planet_4'        .and. &
       trim(sst_method) /= 'aqua_planet_5'        .and. &
       trim(sst_method) /= 'aqua_planet_6'        .and. &
       trim(sst_method) /= 'aqua_planet_7'        .and. &
       trim(sst_method) /= 'aqua_planet_8'        .and. &
       trim(sst_method) /= 'aqua_planet_10N'      .and. &
       trim(sst_method) /= 'aqua_planet_15N'      .and. &
       trim(sst_method) /= 'aqua_planet_20N'      .and. &
       trim(sst_method) /= 'aqua_planet_25N'      .and. &
       trim(sst_method) /= 'aqua_planet_30N'      .and. &
       trim(sst_method) /= 'aqua_planet_35N'      .and. &
       trim(sst_method) /= 'aqua_planet_40N'      .and. &
       trim(sst_method) /= 'aqua_planet_45N'      .and. &
       trim(sst_method) /= 'aqua_planet_50N'      .and. &
       trim(sst_method) /= 'aqua_planet_55N'      .and. &
       trim(sst_method) /= 'aqua_planet_60N'      .and. &
       trim(sst_method) /= 'aqua_planet_65N'      .and. &
       trim(sst_method) /= 'aqua_planet_70N'      .and. &
       trim(sst_method) /= 'aqua_planet_75N'      .and. &
       trim(sst_method) /= 'aqua_planet_80N'      .and. &
       trim(sst_method) /= 'aqua_planet_85N'      .and. &
       trim(sst_method) /= 'aqua_planet_90N'      .and. &
       trim(sst_method) /= 'aqua_walker'          .and. &
       trim(sst_method) /= 'aqua_walker_cos'      .and. &
       trim(sst_method) /= 'aqua_walker_guass'    .and. &
       trim(sst_method) /= 'aqua_walker_guass_b'  .and. &
       trim(sst_method) /= 'aqua_walker_guass_c'  .and. &
       trim(sst_method) /= 'aqua_walker_guass_d'  .and. &
       trim(sst_method) /= 'mixed_layer' ) call error_mesg &
     ('ice_model_init', 'namelist variable sst_method has invalid value', FATAL)

!----------------------------------------------------------

  nlon = size(glon_bnd,1)-1
  nlat = size(glon_bnd,2)-1

!----------------------------------------------------------
!--- set up domain type ---

! if (present(Atmos_domain)) then
!     call mpp_get_layout (Atmos_domain, layout)
! else
!     call mpp_define_layout  ( (/1,nlon,1,nlat/), mpp_npes(), layout )
! endif

! call mpp_define_domains ( (/1,nlon,1,nlat/), layout, Ice%Domain, &
!                           name='ice grid')

  Ice%Domain => Atmos_domain

!----------------------------------------------------------
! get global domain indices
! this assumes that domain2d type has been assigned

  call mpp_get_global_domain ( Ice%Domain, isg, ieg, jsg, jeg )

  allocate ( Ice%glon_bnd (isg:ieg+1,jsg:jeg+1), &
             Ice%glat_bnd (isg:ieg+1,jsg:jeg+1), &
             Ice%glon     (isg:ieg,jsg:jeg), &
             Ice%glat     (isg:ieg,jsg:jeg), &
             Ice%gmask    (isg:ieg,jsg:jeg)  )

!----------------------------------------------------------
!--- read grid info ----

   Ice%glon_bnd = glon_bnd
   Ice%glat_bnd = glat_bnd
  !do j = js, je+1
  !do i = is, ie+1
  !   Ice%glon_bnd(i,j) = glon_bnd(i-isg+1,j-jsg+1)
  !   Ice%glat_bnd(i,j) = glat_bnd(i-isg+1,j-jsg+1)
  !enddo
  !enddo

  ! read the land mask from a file (land=1)
  if (open_file(land_mask_fileobj, 'INPUT/land_mask.nc', 'read', Ice%domain)) then
      call read_data (land_mask_fileobj, 'land_mask', Ice%glon)
      where (Ice%glon > 0.50)
         Ice%gmask = .false.
      elsewhere
         Ice%gmask = .true.
      endwhere
      call close_file(land_mask_fileobj)
  else
      Ice%gmask = .true.  ! aqua-planet
  endif

   ! mid-point of grid box
  !if (is_latlon(glon_bnd,latb_out)) then
  !   do j = js, je
  !   do i = is, ie
  !      Ice%glon(i,j) = (Ice%glon_bnd(i,j  )+Ice%glon_bnd(i+1,j  )+ &
  !                       Ice%glon_bnd(i,j+1)+Ice%glon_bnd(i+1,j+1))*0.25
  !      Ice%glat(i,j) = (Ice%glat_bnd(i,j  )+Ice%glat_bnd(i+1,j  )+ &
  !                       Ice%glat_bnd(i,j+1)+Ice%glat_bnd(i+1,j+1))*0.25
  !   enddo
  !   enddo
  !else
      call get_cell_center (Ice%glon_bnd, Ice%glat_bnd, Ice%glon, Ice%glat)
  !endif

!----------------------------------------------------------
! get compute domain indices

  call mpp_get_compute_domain ( Ice%Domain, is, ie, js, je )

  allocate ( Ice%lon_bnd        (is:ie+1,js:je+1), &
             Ice%lat_bnd        (is:ie+1,js:je+1), &
             Ice%lon            (is:ie, js:je)   , &
             Ice%lat            (is:ie, js:je)   , &
             Ice%ice_mask       (is:ie, js:je)   , &
             Ice%t_surf         (is:ie, js:je)   , &
             Ice%albedo         (is:ie, js:je)   , &
             Ice%albedo_vis_dir (is:ie, js:je)   , &
             Ice%albedo_nir_dir (is:ie, js:je)   , &
             Ice%albedo_vis_dif (is:ie, js:je)   , &
             Ice%albedo_nir_dif (is:ie, js:je)   , &
             Ice%rough_mom      (is:ie, js:je)   , &
             Ice%rough_heat     (is:ie, js:je)   , &
             Ice%rough_moist    (is:ie, js:je)   , &
             Ice%thickness      (is:ie, js:je)   , &
             Ice%mask           (is:ie, js:je)     )

    Ice%lon_bnd    = Ice%glon_bnd(is:ie+1,js:je+1)
    Ice%lat_bnd    = Ice%glat_bnd(is:ie+1,js:je+1)
    Ice%lon        = Ice%glon(is:ie, js:je)
    Ice%lat        = Ice%glat(is:ie, js:je)
    Ice%mask       = Ice%gmask(is:ie, js:je)
    Ice%Time           = Time
    Ice%Time_init      = Time_init
    Ice%Time_step_fast = Time_step_fast
    Ice%Time_step_slow = Time_step_slow

!----------------------------------------------------------
!----------- read restart -------------

need_ic = .false.

if (open_file(ice_restart_fileobj, 'INPUT/ice_model.res.nc', 'read', Ice%domain, is_restart=.true.)) then
   if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
            'Reading NetCDF formatted restart file: INPUT/ice_model.res.nc', NOTE)

   call read_data(ice_restart_fileobj, 'mlon', mlon)
   call read_data(ice_restart_fileobj, 'mlat', mlat)
   if (mlon /= nlon .or. mlat /= nlat )  &
        call error_mesg ('ice_model_init',           &
                        'incorrect resolution on restart', FATAL)

   call ice_register_restart(ice_restart_fileobj, Ice)
   call read_restart(ice_restart_fileobj)
   call close_file(ice_restart_fileobj)
else
  !--- if no restart then no ice ---
      need_ic = .true.
      Ice%t_surf      = TFREEZE   ! + temp_ice_freeze
      Ice%thickness   = 0.0       ! no ice initially
      Ice%albedo      = 0.14
      Ice%albedo_vis_dir = 0.14
      Ice%albedo_nir_dir = 0.14
      Ice%albedo_vis_dif = 0.14
      Ice%albedo_nir_dif = 0.14
      Ice%rough_mom   = 0.0004
      Ice%rough_heat  = 0.0004
      Ice%rough_moist = 0.0004
    ! Ice%ice_mask    = .false. ! no ice initially

    ! fixed roughness
      call fixed_ocean_roughness ( Ice%mask, Ice%rough_mom, &
                                   Ice%rough_heat, Ice%rough_moist )
endif

  ! initialize mask where ice exists
    Ice%ice_mask = Ice%mask .and. Ice%thickness .ge. thickness_min

    call ice_albedo_init (TFREEZE)

  ! analytic distribution with no ice
  ! melt all ice
    if (sst_method == "aqua_planet_1") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**2) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_2") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-min(3.*abs(Ice%lat)/pi,1.)) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_3") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**4) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_4") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = max(min(1.5*Ice%lat,pi*0.5),-pi*0.5) ! use t_surf as work array
            Ice%t_surf = 27.*(1.-0.5*(sin(Ice%t_surf)**2+sin(Ice%t_surf)**4)) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_5") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_6" .or. sst_method == "aqua_planet_7") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        lond = 30.*pi/180.
        latd = 15.*pi/180.
        if (sst_method == "aqua_planet_6") amp = 1.
        if (sst_method == "aqua_planet_7") amp = 3.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
               t_control = 27.*(1.-sin(max(min(1.5*Ice%lat(i,j),pi*0.5),-pi*0.5))**2) + TFREEZE
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = t_control + amp * cos(0.5*pi*min(max(dellon/lond,-1.),1.))**2 * &
                                                   cos(0.5*pi*min(max(Ice%lat(i,j)/latd,-1.),1.))**2
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        !lond = 30.*pi/180.
        lond = pi/36. ! multiply original lond by 1/6
        latd = 15.*pi/180.
        amp = 8.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               !Ice%t_surf(i,j) = 27. + TFREEZE + amp * cos(0.5*pi*min(max(dellon/lond,-1.),1.))**2
               Ice%t_surf(i,j) = 297. + amp * cos(0.5*pi*min(max(dellon/lond,-1.),1.))**2
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker_guass") then ! see Wofsy and Kuang
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        siggy = (pi/8.)*(1./6.) ! sigma in the guassian distribution
        !amp = 1./(siggy*SQRT(2*pi))
        amp = 8.
        do j = js, je
             tempi=real(j)
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = 297. + amp*EXP(-0.5*((dellon)**2.)/(siggy)**2)
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker_guass_b") then ! see Wofsy and Kuang
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        siggy = (pi/12.)*(1./6.) ! sigma in the guassian distribution
        !amp = 1./(siggy*SQRT(2*pi))
        amp = 8.
        do j = js, je
             tempi=real(j)
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = 297. + amp*EXP(-0.5*((dellon)**2.)/(siggy)**2)
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker_guass_c") then ! see Wofsy and Kuang
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        siggy = (pi/16.)*(1./6.) ! sigma in the guassian distribution
        !amp = 1./(siggy*SQRT(2*pi))
        amp = 8.
        do j = js, je
             tempi=real(j)
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = 297. + amp*EXP(-0.5*((dellon)**2.)/(siggy)**2)
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker_guass_d") then ! see Wofsy and Kuang
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        siggy = (pi/20.)*(1./6.) ! sigma in the guassian distribution
        !amp = 1./(siggy*SQRT(2*pi))
        amp = 4.
        do j = js, je
             tempi=real(j)
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = 297. + amp*EXP(-0.5*((dellon)**2.)/(siggy)**2)
            endif
        enddo
        enddo
    else if (sst_method == "aqua_walker_cos") then ! see Wofsy and Kuang
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lond = pi/36. ! multiply original lond by 1/6
        lon0 = 0.
        amp = 8.
        dom_wid=Ice%lon(is,js)-Ice%lon(is,je)
        !dom_wid=real(ie)-real(is)
        do j = js, je
             tempi=real(j)
        do i = is, ie
           if (Ice%mask(i,j)) then
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = 297. - amp * cos(0.5*pi*min(max(dellon/lond,-1.),1.))
            endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_8") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        latd = 30.*pi/180.
        amp = 3.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**2) + TFREEZE + &
                         amp * cos(Ice%lon-lon0) * cos(0.5*pi*min(max(Ice%lat/latd,-1.),1.))**2
        endwhere
    else if (sst_method == "aqua_planet_10N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/18.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/18.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/18.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_15N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/12.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/12.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/12.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_20N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/9.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/9.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/9.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_25N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 5.*pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-5.*pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-5.*pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_30N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/6.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/6.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/6.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_35N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 7.*pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-7.*pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-7.*pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_40N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 2.*pi/9.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-2.*pi/9.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-2.*pi/9.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_45N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/4.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/4.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/4.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_50N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 5.*pi/18.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-5.*pi/18.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-5.*pi/18.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_55N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 11.*pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-11.*pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-11.*pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_60N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/3.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/3.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/3.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_65N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 13.*pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-13.*pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-13.*pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_70N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 7.*pi/18.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-7.*pi/18.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-7.*pi/18.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_75N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 5.*pi/12.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-5.*pi/12.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-5.*pi/12.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_80N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 4.*pi/9.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-4.*pi/9.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-4.*pi/9.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_85N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > 17.*pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-17.*pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-17.*pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_90N") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/2.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/2.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/2.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    endif


!----------------------------------------------------------

  if (trim(ice_method) == 'prognostic' .or. &
      trim(ice_method) == 'uniform') then
      if (trim(interp_method) == "conservative") then
          Amip_ice = amip_interp_new ( Ice%lon_bnd(:,1),     Ice%lat_bnd(1,:),  &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_ice, use_annual=use_annual_ice )
      else if(trim(interp_method) == "bilinear") then
          Amip_ice = amip_interp_new ( Ice%lon,     Ice%lat,          &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_ice, use_annual=use_annual_ice )
      else
          call error_mesg ('ice_model_init', 'interp_method should be '// &
                           'conservative or bilinear', FATAL)
      endif
      ! initialize ice (if needed)
      if (need_ic) then
         call prognostic_ice ( Ice )
      endif
  endif

!----------------------------------------------------------

  if (trim(sst_method) == 'specified') then
      if (trim(interp_method) == "conservative") then
          Amip_sst = amip_interp_new ( Ice%lon_bnd(:,1),     Ice%lat_bnd(1,:),  &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_sst, use_annual=use_annual_sst )
      else if(trim(interp_method) == "bilinear") then
          Amip_sst = amip_interp_new ( Ice%lon,     Ice%lat,          &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_sst, use_annual=use_annual_sst )
      else
          call error_mesg ('ice_model_init', 'interp_method should be '// &
                           'conservative or bilinear', FATAL)
      endif
      ! initialize sst (if needed)
      if (need_ic) then
         call prognostic_sst ( Ice )
      endif
  endif

print *, 'pe,count(ice,all,ocean)=',mpp_pe(),count(Ice%ice_mask),count(Ice%mask),count(Ice%mask .and. .not.Ice%ice_mask)

! add on non-zero sea surface temperature perturbation (namelist option)
! this perturbation may be useful in accessing model sensitivities

  if ( abs(sst_anom) > 0.0001 ) then
    Ice%t_surf(:,:) = Ice%t_surf(:,:) + sst_anom
  endif

!----------------------------------------------------------

  module_is_initialized = .true.

!----------------------------------------------------------

 end subroutine ice_model_init

 subroutine ice_register_restart(fileobj, Ice)

 type(FmsNetcdfDomainFile_t), intent(inout) :: fileobj    !< Ice restart domain decomposed fileobj
 type(ice_data_type), intent(inout)         :: Ice        !< Ice data type
 character(len=8), dimension(3)             ::  dim_names !< Array of dimension names

 dim_names(1) = "xaxis_1"
 dim_names(2) = "yaxis_1"
 dim_names(3) = "Time"

 call register_axis(fileobj, dim_names(1), "x")
 call register_axis(fileobj, dim_names(2), "y")
 call register_axis(fileobj, dim_names(3), unlimited)

 !< Register the domain decomposed dimensions as variables so that the combiner can work
 !! correctly
 call register_field(fileobj, dim_names(1), "double", (/dim_names(1)/))
 call register_field(fileobj, dim_names(2), "double", (/dim_names(2)/))

 call register_restart_field ( fileobj, 't_surf',         Ice%t_surf,         dim_names )
 call register_restart_field ( fileobj, 'thickness',      Ice%thickness,      dim_names )
 call register_restart_field ( fileobj, 'albedo',         Ice%albedo,         dim_names )
 call register_restart_field ( fileobj, 'albedo_vis_dir', Ice%albedo_vis_dir, dim_names )
 call register_restart_field ( fileobj, 'albedo_nir_dir', Ice%albedo_nir_dir, dim_names )
 call register_restart_field ( fileobj, 'albedo_vis_dif', Ice%albedo_vis_dif, dim_names )
 call register_restart_field ( fileobj, 'albedo_nir_dif', Ice%albedo_nir_dif, dim_names )
 call register_restart_field ( fileobj, 'rough_mom',      Ice%rough_mom,      dim_names )
 call register_restart_field ( fileobj, 'rough_heat',     Ice%rough_heat,     dim_names )
 call register_restart_field ( fileobj, 'rough_moist',    Ice%rough_moist,    dim_names )

 end subroutine ice_register_restart

!######################################################################

 subroutine ice_model_end ( Ice )
 type(ice_data_type), intent(inout) :: Ice
 integer :: unit
 character(len=64) :: fname='RESTART/ice_model.res.nc'
 type(FmsNetcdfDomainFile_t) :: ice_restart_fileobj !< Ice restart domain decomposed fileobj

 if (.not.module_is_initialized) return
 if( do_netcdf_restart) then

    if(mpp_pe() == mpp_root_pe() ) then
       call error_mesg ('ice_model_mod', 'Writing NetCDF formatted restart file: RESTART/ice_model.res.nc', NOTE)
    endif

    if (open_file(ice_restart_fileobj, fname, 'overwrite', Ice%domain, is_restart=.true.)) then
        call ice_register_restart(ice_restart_fileobj, Ice)
        call register_field(ice_restart_fileobj, "mlon", "double")
        call register_field(ice_restart_fileobj, "mlat", "double")
        call write_restart(ice_restart_fileobj)
        call write_data(ice_restart_fileobj, 'mlon', size(Ice%gmask,1))
        call write_data(ice_restart_fileobj, 'mlat', size(Ice%gmask,2))
        call add_domain_dimension_data(ice_restart_fileobj)
        call close_file(ice_restart_fileobj)
    endif !< if(open_file)
 endif


  deallocate ( Ice%glon_bnd, Ice%glat_bnd, Ice%glon, Ice%glat, Ice%gmask )
  deallocate ( Ice%lon_bnd, Ice%lat_bnd, Ice%lon, Ice%lat, Ice%ice_mask,       &
               Ice%t_surf, Ice%albedo, Ice%albedo_vis_dir, Ice%albedo_nir_dir, &
               Ice%albedo_vis_dif, Ice%albedo_nir_dif, Ice%rough_mom,          &
               Ice%rough_heat, Ice%rough_moist, Ice%thickness, Ice%mask        )

  module_is_initialized = .false.

 end subroutine ice_model_end

 !< Add_dimension_data: Adds dummy data for the domain decomposed axis
 subroutine add_domain_dimension_data(fileobj)
  type(FmsNetcdfDomainFile_t) :: fileobj !< Fms2io domain decomposed fileobj
  integer, dimension(:), allocatable :: buffer !< Buffer with axis data
  integer :: is, ie !< Starting and Ending indices for data

    call get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)

    call get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)
 end subroutine add_domain_dimension_data

!######################################################################
!               Routines added for computing then
!            mid-point of grid boxes with cubed sphere
!######################################################################

function is_latlon ( lon, lat )
real, intent(in) :: lon(:,:), lat(:,:)

! Determines if the latitude/longitude values form a
! regular latitude/longitude grid (or something else).
!

logical :: is_latlon
integer :: i, j

  is_latlon = .true.

  do j = 2, size(lon,2)
  do i = 1, size(lon,1)
     if (lon(i,j) .ne. lon(i,1)) then
        is_latlon = .false.
        return
     endif
  enddo
  enddo

  do j = 1, size(lat,2)
  do i = 2, size(lat,1)
     if (lat(i,j) .ne. lat(1,j)) then
        is_latlon = .false.
        return
     endif
  enddo
  enddo

end function is_latlon

!########################################################################

  subroutine get_cell_center (lonb, latb, lon, lat)
    !------------------------------------------------------------------!
    ! calculate cell center (lon,lat)                                  !
    ! by averaging cell corner locations (lonb,latb) in Cartesian coor !
    !------------------------------------------------------------------!
    real, dimension(:,:), intent(in)  :: lonb, latb
    real, dimension(:,:), intent(out) :: lon, lat
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real :: sph_coor(2), xyz_corner(3,size(lonb,1),size(latb,2)), xyz_center(3), abs_center
    integer :: i,j, nx,ny

    nx = size(lonb,1)
    ny = size(latb,2)

    do j=1,ny
       do i=1,nx
          sph_coor(1)=lonb(i,j)
          sph_coor(2)=latb(i,j)
          call latlon2xyz(sph_coor, xyz_corner(1,i,j))
       enddo
    enddo
    do j=1,ny-1
       do i=1,nx-1
          xyz_center(:)=0.25*(xyz_corner(:,i  ,j  )+xyz_corner(:,i+1,j  )  &
                             +xyz_corner(:,i  ,j+1)+xyz_corner(:,i+1,j+1))
          abs_center=xyz_center(1)*xyz_center(1)                           &
                    +xyz_center(2)*xyz_center(2)                           &
                    +xyz_center(3)*xyz_center(3)
          xyz_center(:)=xyz_center(:)/abs_center
          call xyz2latlon(xyz_center, sph_coor)
          lon(i,j)=sph_coor(1)
          lat(i,j)=sph_coor(2)
       enddo
    enddo

  end subroutine get_cell_center

!########################################################################

  subroutine latlon2xyz(sph_coor, xyz_coor)
    !------------------------------------------------------------------!
    ! calculate cartesian coordinates from spherical coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! sph_coor [rad]   latlon coordinate                               !
    !                                                                  !
    ! output:                                                          !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !------------------------------------------------------------------!
    real, dimension(2), intent(in)    :: sph_coor
    real, dimension(3), intent(inout) :: xyz_coor

    xyz_coor(1) = cos(sph_coor(2)) * cos(sph_coor(1))
    xyz_coor(2) = cos(sph_coor(2)) * sin(sph_coor(1))
    xyz_coor(3) = sin(sph_coor(2))

  end subroutine latlon2xyz
  !====================================================================!
  subroutine xyz2latlon(xyz_coor, sph_coor)
    !------------------------------------------------------------------!
    ! calculate spherical coordinates from cartesian coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !                                                                  !
    ! output:                                                          !
    ! sph_coor [rad]   latlon coordinate                               !
    !------------------------------------------------------------------!

    real, dimension(3), intent(in)    :: xyz_coor
    real, dimension(2), intent(inout) :: sph_coor
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real :: coslat, radius
    real, parameter :: epsilon=1.e-7
    integer :: i,j


    radius=sqrt(xyz_coor(1)*xyz_coor(1) &
               +xyz_coor(2)*xyz_coor(2) &
               +xyz_coor(3)*xyz_coor(3))

    sph_coor(1)=atan2(xyz_coor(2),xyz_coor(1))
    if (sph_coor(1)<0.) sph_coor(1)=sph_coor(1)+2*PI
    sph_coor(2)=asin(xyz_coor(3)/radius)

  end subroutine xyz2latlon

!######################################################################

end module ice_model_mod

