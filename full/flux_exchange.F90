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

module flux_exchange_mod

!> \author Bruce Wyman <Bruce.Wyman@noaa.gov>
!! \author V. Balaji <V.Balaji@noaa.gov>
!! \author Sergey Malyshev <Sergey.Malyshev@noaa.gov>
!!
!! \brief The flux_exchange module provides interfaces to couple the following component
!!        models: atmosphere, ocean, land, and ice. All interpolation between physically
!!        distinct model grids is handled by the exchange grid (xgrid_mod) with the
!!        interpolated quantities being conserved.
!!
!! -# This version of flux_exchange_mod allows the definition of physically independent
!!    grids for atmosphere, land and sea ice. Ice and ocean must share the same physical
!!    grid (though the domain decomposition on parallel systems may be different).
!!    Grid information is input through the grid_spec file (URL). The masked region of the
!!    land grid and ice/ocean grid must "tile" each other. The masked region of the ice grid
!!    and ocean grid must be identical.
!! <pre>
!!         ATMOSPHERE  |----|----|----|----|----|----|----|----|
!!
!!               LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
!!
!!                ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!!
!!               OCEAN |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
!! </pre>
!!
!!    where \c |xxx| represents a masked grid point
!!
!!    The atmosphere, land, and ice grids exchange information using the exchange grid xmap_sfc.
!!
!!    The land and ice grids exchange runoff data using the exchange grid xmap_runoff.
!!
!!    Transfer of data between the ice bottom and ocean does not require an exchange
!!    grid as the grids are physically identical. The flux routines will automatically
!!    detect and redistribute data if their domain decompositions are different.
!!
!!    To get information from the atmosphere to the ocean it must pass through the
!!    ice model, first by interpolating from the atmospheric grid to the ice grid,
!!    and then transferring from the ice grid to the ocean grid.
!!
!! -# Each component model must have a public defined data type containing specific
!!    boundary fields. A list of these quantities is located in the NOTES of this document.
!!
!! -# The surface flux of sensible heat and surface evaporation can be implicit functions
!!    of surface temperature. As a consequence, the parts of the land and sea-ice models
!!    that update the surface temperature must be called on the atmospheric time step
!!
!! -# The surface fluxes of all other tracers and of momentum are assumed to be explicit
!!    functions of all surface parameters.
!!
!! -# While no explicit reference is made within this module to the implicit treatment
!!    of vertical diffusion in the atmosphere and in the land or sea-ice models, the
!!    module is designed to allow for simultaneous implicit time integration on both
!!    sides of the surface interface.
!!
!! -# Due to #5, the diffusion part of the land and ice models must be called on the
!!    atmospheric time step, although in the case of concurrent-ice coupling, this
!!    version of the sea-ice that is called by the atmosphere may later be replaced
!!    by a version of the ice that is tightly coupled with the ocean.
!!
!! -# The fluxes of additional tracers related to biological quantities or the
!!    air-sea exchange of gases are accomplished by specifying fields that will
!!    be passed between components via the "field_table" and the use of named
!!    fields in the coupler_..._bc_types.
!!
!! -# Any field passed from one component to another may be "faked" to a
!!    constant value, or to data acquired from a file, using the
!!   data_override feature of FMS. The fields to override are runtime
!!   configurable, using the text file <tt>data_table</tt> for input.
!!   See the data_override_mod documentation for more details.
!!
!!   We DO NOT RECOMMEND exercising the data override capabilities of
!!   the FMS coupler until the user has acquired considerable
!!   sophistication in running FMS.
!!
!!   Here is a listing of the override capabilities of the flux_exchange
!!   module:
!!
!!   - FROM the atmosphere boundary TO the exchange grid (in sfc_boundary_layer):
!!
!!        t_bot, q_bot, z_bot, p_bot, u_bot, v_bot, p_surf, slp, gust
!!
!!   - FROM the ice boundary TO the exchange grid (in sfc_boundary_layer):
!!
!!        t_surf, rough_mom, rough_heat, rough_moist, albedo, u_surf, v_surf
!!
!!   - FROM the land boundary TO the exchange grid (in sfc_boundary_layer):
!!
!!        t_surf, t_ca, q_ca, rough_mom, rough_heat, albedo
!!
!!   - FROM the exchange grid TO land_ice_atmos_boundary (in
!!     sfc_boundary_layer):
!!
!!        t, albedo, land_frac, dt_t, dt_q, u_flux, v_flux, dtaudu, dtaudv,
!!        u_star, b_star, rough_mom
!!
!!   - FROM the atmosphere boundary TO the exchange grid (in
!!     flux_down_from_atmos):
!!
!!        flux_sw, flux_lw, lprec, fprec, coszen, dtmass, delta_t,
!!        delta_q, dflux_t, dflux_q
!!
!!   - FROM the exchange grid TO the land boundary (in
!!     flux_down_from_atmos):
!!
!!     t_flux, q_flux, lw_flux, sw_flux, lprec, fprec, dhdt, dedt, dedq,
!!     drdt, drag_q, p_surf
!!
!!   - FROM the exchange grid TO the ice boundary (in flux_down_from_atmos):
!!
!!        u_flux, v_flux, t_flux, q_flux, lw_flux, lw_flux_dn, sw_flux,
!!        sw_flux_dn, lprec, fprec, dhdt, dedt, drdt, coszen, p
!!
!!   - FROM the land boundary TO the ice boundary (in flux_land_to_ice):
!!
!!        runoff, calving
!!
!!   - FROM the ice boundary TO the ocean boundary (in flux_ice_to_ocean):
!!
!!        u_flux, v_flux, t_flux, q_flux, salt_flux, lw_flux, sw_flux,
!!        lprec, fprec, runoff, calving, p, ustar_berg, area_berg, mass_berg
!!
!!   - FROM the ocean boundary TO the ice boundary (in flux_ocean_to_ice):
!!
!!        u, v, t, s, frazil, sea_level
!!
!!   - FROM the ice boundary TO the atmosphere boundary (in flux_up_to_atmos):
!!
!!        t_surf
!!
!!   - FROM the land boundary TO the atmosphere boundary (in
!!     flux_up_to_atmos):
!!
!!        t_ca, t_surf, q_ca
!!
!!   See NOTES below for an explanation of the field names.
!!
!! \section diag_fields Diagnostic Fields
!!
!! The table below contains the available diagnostic fields is the `flux` diagnostic module.
!!
!! Field Name  | Units           | Description
!! ----------- | --------------- | -----------
!! land_mask   | none            | Fractional amount of land
!! wind        | m/s             | Wind speed for flux calculations
!! drag_moist  | none            | Drag coeff for moisture
!! drag_heat   | none            | Drag coeff for heat
!! drag_mom    | none            | Drag coeff for momentum
!! rough_moist | m               | Surface roughness for moisture
!! rough_heat  | m               | Surface roughness for heat
!! rough_mom   | m               | Surface roughness for momentum
!! u_star      | m/s             | Friction velocity
!! b_star      | m/s             | Buoyancy scale
!! q_star      | kg water/kg air | moisture scale
!! t_atm       | deg_k           | temperature at btm level
!! u_atm       | m/s             | u wind component at btm level
!! v_atm       | m/s             | v wind component at btm level
!! q_atm       | kg/kg           | specific humidity at btm level
!! p_atm       | pa              | pressure at btm level
!! z_atm       | m               | height of btm level
!! gust        | m/s             | gust scale
!! rh_ref      | percent         | relative humidity at ref height
!! t_ref       | deg_k           | temperature at ref height
!! u_ref       | m/s             | zonal wind component at ref height
!! v_ref       | m/s             | meridional wind component at ref height
!! del_h       | none            | ref height interp factor for heat
!! del_m       | none            | ref height interp factor for momentum
!! del_q       | none            | ref height interp factor for moisture
!! tau_x       | pa              | zonal wind stress
!! tau_y       | pa              | meridional wind stress
!! ice_mask    | none            | fractional amount of sea ice
!! t_surf      | deg_k           | surface temperature
!! t_ca        | deg_k           | canopy air temperature
!! q_surf      | kg/kg           | surface specific humidity
!! shflx       | w/m2            | sensible heat flux
!! evap        | kg/m2/s         | evaporation rate
!! lwflx       | w/m2            | net (down-up) longwave flux
!!
!! \section flux_exchange_config Flux Exchange Configuration
!!
!! flux_exchange_mod is configured via the flux_exchange_nml namelist in the `input.nml` file.
!! The following table are the available namelist variables.
!!
!! <table>
!!   <tr>
!!     <th>Variable Name</th>
!!     <th>Type</th>
!!     <th>Default Value</th>
!!     <th>Description</th>
!!   </tr>
!!   <tr>
!!     <td>z_ref_heat</td>
!!     <td>real</td>
!!     <td>2.0</td>
!!     <td>Reference height (meters) for temperature and relative humidity
!!       diagnostics (t_ref, rh_ref, del_h, del_q).</td>
!!   </tr>
!!   <tr>
!!     <td>z_ref_mom</td>
!!     <td>real</td>
!!     <td>10.0</td>
!!     <td>Reference height (meters) for mementum diagnostics (u_ref, v_ref,
!!       del_m).</td>
!!   </tr>
!!   <tr>
!!     <td>ex_u_start_smooth_bug</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>By default, the global exchange grid `u_star` will not be interpolated
!!       from atmospheric grid, this is different from Jakarta behavior and will
!!       change answers.  So to preserve Jakarta behavior and reproduce answers
!!       explicitly set this namelist variable to .true. in input.nml.</td>
!!   </tr>
!!   <tr>
!!     <td>sw1way_bug</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!!   <tr>
!!     <td>do_area_weighted_flux</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!!   <tr>
!!     <td>debug_stocks</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!!   <tr>
!!     <td>divert_stocks_report</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!!   <tr>
!!     <td>do_runoff</td>
!!     <td>logical</td>
!!     <td>.TRUE.</td>
!!     <td>Turns on/off the land runoff interpolation to the ocean.</td>
!!   </tr>
!!   <tr>
!!     <td>do_forecast</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td></td>
!!   </tr>
!!   <tr>
!!     <td>nblocks</td>
!!     <td>integer</td>
!!     <td>1</td>
!!     <td>Specify number of blocks that n_xgrid_sfc is divided into. The main
!!       purpose is for Openmp implementation. Normally you may set nblocks to be
!!       coupler_nml atmos_nthreads.</td>
!!   </tr>
!!   <tr>
!!     <td>partition_fprec_from_lprec</td>
!!     <td>logical</td>
!!     <td>.FALSE.</td>
!!     <td>Option for ATM override experiments where liquid+frozen precip are combined.
!!         This option will convert liquid precip to snow when t_ref is less than tfreeze parameter</td>
!!   </tr>
!!   <tr>
!!     <td>scale_precip_2d</td>
!!     <td>logical</td>
!!     <td>.false.</td>
!!     <td>Option to scale the Atm%lprec.
!!         If this varible is set to .true. Atm%lprec will be rescaled by a field read from the data_table</td>
!!   </tr>
!!
!! \section main_example Main Program Example
!!
!! ~~~~~~~~~~{.f90}
!! DO slow time steps (ocean)
!!    call flux_ocean_to_ice
!!
!!    call ICE_SLOW_UP
!!
!!    DO fast time steps (atmos)
!!       call sfc_boundary_layer
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
!!    END DO
!!
!!    call ICE_SLOW_DN
!!
!!    call flux_ice_to_ocean
!!
!!    call OCEAN
!! END DO
!! ~~~~~~~~~~
!!
!! \note  LAND_FAST and ICE_FAST must update the surface temperature
!!
!! \section Required Variables In Defined Data Types For Component Models
!!
!! \subsection Atmosphere
!!
!! ~~~~~~~~~~{.f90}
!! type (atmos_boundary_data_type) :: Atm
!!
!! real, dimension(:) :: Atm%lon_bnd & ! longitude axis grid box boundaries in radians
!!                                     ! must be monotonic
!!                       Atm%lat_bnd   ! latitude axis grid box boundaries in radians
!!                                     ! must be monotonic
!! real, dimension(:,:) :: Atm%t_bot   & ! temperature at lowest model level
!!                         Atm%q_bot   & ! specific humidity at lowest model level
!!                         Atm%z_bot   & !    height above the surface for the lowest model level (m)
!!                         Atm%p_bot   & !    pressure at lowest model level (pa)
!!                         Atm%u_bot   & !    zonal wind component at lowest model level (m/s)
!!                         Atm%v_bot   & !    meridional wind component at lowest model level (m/s)
!!                         Atm%p_surf  & !   surface pressure (pa)
!!                         Atm%slp     & !   sea level pressure (pa)
!!                         Atm%gust    & !   gustiness factor (m/s)
!!                         Atm%flux_sw & !  net shortwave flux at the surface
!!                         Atm%flux_lw & !  downward longwave flux at the surface
!!                         Atm%lprec   & !  liquid precipitation (kg/m2)
!!                         Atm%fprec   & !  water equivalent frozen precipitation (kg/m2)
!!                         Atm%coszen  & !  cosine of the zenith angle
!! integer, dimension(4) :: Atm%axes ! Axis identifiers returned by diag_axis_init for the
!!                                   ! atmospheric model axes: X, Y, Z_full, Z_half.
!! ~~~~~~~~~~
!!
!! The following five fields are gathered into a data type for convenience in passing
!! this information through the different levels of the atmospheric model --
!! these fields are rlated to the simultaneous implicit time steps in the
!! atmosphere and surface models -- they are described more fully in
!! flux_exchange.tech.ps and in the documntation for vert_diff_mod
!!
!! ~~~~~~~~~~{.f90}
!! type (surf_diff_type) :: Atm%Surf_Diff
!!
!! real, dimension(:,:) :: Atm%Surf_Diff%dtmass  & ! dt/mass where dt = atmospheric time step ((i+1) = (i-1) for leapfrog) (s)
!!                                                 ! mass = mass per unit area of lowest atmosphehic layer  (Kg/m2))
!!                         Atm%Surf_Diff%delta_t & ! increment ((i+1) = (i-1) for leapfrog) in temperature of
!!                                                 ! lowest atmospheric layer  (K)
!!                         Atm%Surf_Diff%delta_q & ! increment ((i+1) = (i-1) for leapfrog) in specific humidity of
!!                                                 ! lowest atmospheric layer (nondimensional -- Kg/Kg)
!!                         Atm%Surf_Diff%dflux_t & ! derivative of implicit part of downward temperature flux at top of lowest
!!                                                 ! atmospheric layer with respect to temperature
!!                                                 ! of lowest atmospheric layer (Kg/(m2 s))
!!                         Atm%Surf_Diff%dflux_q   ! derivative of implicit part of downward moisture flux at top of lowest
!!                                                 ! atmospheric layer with respect to specific humidity of
!!                                                 ! of lowest atmospheric layer (Kg/(m2 s))
!! ~~~~~~~~~~
!!
!! \subsection Land
!!
!! ~~~~~~~~~~{.f90}
!! type (land_boundary_data_type) :: Land
!!
!! real, dimension(:) :: Land%lon_bnd & ! longitude axis grid box boundaries in radians
!!                                      ! must be monotonic
!!                       Land%lat_bnd   ! latitude axis grid box boundaries in radians
!!                                      ! must be monotonic
!!
!! logical, dimension(:,:,:) :: Land%mask & ! land/sea mask (true for land)
!!                              Land%glacier ! glacier mask  (true for glacier)
!!
!! real, dimension(:,:,:) :: Land%tile_size  & !  fractional area of each tile (partition)
!!                           Land%t_surf     & ! surface temperature (deg k)
!!                           Land%albedo     & ! surface albedo (fraction)
!!                           Land%rough_mom  & ! surface roughness for momentum (m)
!!                           Land%rough_heat & ! surface roughness for heat/moisture (m)
!!                           Land%stomatal   & ! stomatal resistance
!!                           Land%snow       & ! snow depth (water equivalent) (kg/m2)
!!                           Land%water      & ! water depth of the uppermost bucket (kg/m2)
!!                           Land%max_water    ! maximum water depth allowed in the uppermost bucket (kg/m2)
!! ~~~~~~~~~~
!!
!! \subsection Ice
!!
!! ~~~~~~~~~~
!! type (ice_boundary_data_type) :: Ice
!!
!! real, dimension(:) :: Ice%lon_bnd    & ! longitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd    & ! latitude axis grid box boundaries for temperature points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lon_bnd_uv & ! longitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!                       Ice%lat_bnd_uv   ! latitude axis grid box boundaries for momentum points
!!                                        ! in radians (must be monotonic)
!!
!! logical, dimension(:,:,:) :: Ice%mask    & ! ocean/land mask for temperature points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%mask_uv & ! ocean/land mask for momentum points
!!                                            ! (true for ocean, with or without ice)
!!                              Ice%ice_mask  ! optional ice mask (true for ice)
!!
!! real, dimension(:,:,:) :: Ice%part_size  & ! fractional area of each partition of a temperature grid box
!!                           Ice%part_size_uv ! fractional area of each partition of a momentum grid box
!! ~~~~~~~~~~
!!
!! The following fields are located on the ice top grid
!!
!! ~~~~~~~~~~{.f90}
!! real, dimension(:,:,:) :: Ice%t_surf     & ! surface temperature (deg k)
!!                           Ice%albedo     & ! surface albedo (fraction)
!!                           Ice%rough_mom  & ! surface roughness for momentum (m)
!!                           Ice%rough_heat & ! surface roughness for heat/moisture (m)
!!                           Ice%u_surf     & ! zonal (ocean/ice) current at the surface (m/s)
!!                           Ice%v_surf       ! meridional (ocean/ice) current at the surface (m/s)
!! ~~~~~~~~~~
!!
!! The following fields are located on the ice bottom grid
!!
!! ~~~~~~~~~~{.f90}
!! real, dimension(:,:,:) :: Ice%flux_u  & ! zonal wind stress (Pa)
!!                           Ice%flux_v  & ! meridional wind stress (Pa)
!!                           Ice%flux_t  & ! sensible heat flux (w/m2)
!!                           Ice%flux_q  & ! specific humidity flux (kg/m2/s)
!!                           Ice%flux_sw & ! net (down-up) shortwave flux (w/m2)
!!                           Ice%flux_lw & ! net (down-up) longwave flux (w/m2)
!!                           Ice%lprec   & ! mass of liquid precipitation since last time step (Kg/m2)
!!                           Ice%fprec   & ! mass of frozen precipitation since last time step (Kg/m2)
!!                           Ice%runoff    ! mass of runoff water since last time step (Kg/m2)
!! ~~~~~~~~~~
!!
!! \subsection Ocean
!!
!! ~~~~~~~~~~{.f90}
!! type (ocean_boundary_data_type) :: Ocean
!!
!! real, dimension(:) :: Ocean%Data%lon_bnd     & ! longitude axis grid box boundaries for temperature
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lat_bnd     & ! latitude axis grid box boundaries for temperature
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lon_bnd_uv  & ! longitude axis grid box boundaries for momentum
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Data%lat_bnd_uv  & ! latitude axis grid box boundaries for momentum
!!                                                ! points on the ocean DATA GRID (radians)
!!                       Ocean%Ocean%lon_bnd    & ! longitude axis grid box boundaries for temperature
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lat_bnd    & ! latitude axis grid box boundaries for temperature
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lon_bnd_uv & ! longitude axis grid box boundaries for momentum
!!                                                ! points on the ocean MODEL GRID (radians)
!!                       Ocean%Ocean%lat_bnd_uv & ! latitude axis grid box boundaries for momentum
!!                                                ! points on the ocean MODEL GRID (radians)
!! ~~~~~~~~~~
!!
!! \note The data values in all longitude and latitude grid box boundary
!!       array must be monotonic.
!!
!! ~~~~~~~~~~{.f90}
!! logical, dimension(:,:) :: Ocean%Data%mask    & ! ocean/land mask for temperature points on the ocean
!!                                                 ! DATA GRID (true for ocean)
!!                            Ocean%Data%mask_uv & ! ocean/land mask for momentum points on the ocean
!!                                                 ! DATA GRID (true for ocean)
!!                            Ocean%Ocean%mask   & ! ocean/land mask for temperature points on the ocean
!!                                                 ! MODEL GRID (true for ocean)
!!                            Ocean%Ocean%mask_uv  ! ocean/land mask for momentum points on the ocean
!!                                                 ! MODEL GRID (true for ocean)
!! real, dimension(:,:) :: Ocean%t_surf_data & ! surface temperature on the ocean DATA GRID (deg k)
!!                         Ocean%t_surf      & ! surface temperature on the ocean MODEL GRID (deg k)
!!                         Ocean%u_surf      & ! zonal ocean current at the surface on the ocean
!!                                             ! MODEL GRID (m/s)
!!                         Ocean%v_surf      & ! meridional ocean current at the surface on the
!!                                             ! ocean MODEL GRID (m/s)
!!                         Ocean%frazil        ! frazil at temperature points on the ocean MODEL GRID
!! ~~~~~~~~~~

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
!will support 3 types of flux_exchange:
!REGRID: physically distinct grids, via xgrid
!REDIST: same grid, transfer in index space only
!DIRECT: same grid, same decomp, direct copy

  use FMS
  use FMSconstants, only: rdgas, rvgas, cp_air, stefan, WTMAIR, &
                          HLV, HLF, Radius, PI, CP_OCEAN, WTMCO2, WTMC

!! Components
  use land_model_mod,             only: Lnd_stock_pe
  use ocean_model_mod,            only: Ocean_stock_pe
  use atmos_model_mod,            only: Atm_stock_pe
  use atm_land_ice_flux_exchange_mod, only: atm_land_ice_flux_exchange_init, sfc_boundary_layer
  use atm_land_ice_flux_exchange_mod, only: generate_sfc_xgrid, flux_down_from_atmos
  use atm_land_ice_flux_exchange_mod, only: flux_up_to_atmos, atm_stock_integrate, send_ice_mask_sic
  use atm_land_ice_flux_exchange_mod, only: flux_atmos_to_ocean, flux_ex_arrays_dealloc
  use land_ice_flux_exchange_mod,     only: flux_land_to_ice, land_ice_flux_exchange_init
  use ice_ocean_flux_exchange_mod,    only: ice_ocean_flux_exchange_init
  use ice_ocean_flux_exchange_mod,    only: flux_ocean_to_ice, flux_ocean_to_ice_finish
  use ice_ocean_flux_exchange_mod,    only: flux_ice_to_ocean, flux_ice_to_ocean_finish
  use ice_ocean_flux_exchange_mod,    only: flux_ice_to_ocean_stocks, flux_ocean_from_ice_stocks
  use atmos_model_mod,    only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod,    only: ocean_public_type, ice_ocean_boundary_type
  use ocean_model_mod,    only: ocean_state_type
  use ice_model_mod,      only: ice_data_type, land_ice_boundary_type, &
                                ocean_ice_boundary_type, atmos_ice_boundary_type, Ice_stock_pe
  use land_model_mod,     only: land_data_type, atmos_land_boundary_type
  use atmos_ocean_fluxes_mod,     only: atmos_ocean_fluxes_init, atmos_ocean_type_fluxes_init
  use atmos_ocean_fluxes_calc_mod, only: atmos_ocean_fluxes_calc
  use ocean_model_mod,            only: ocean_model_init_sfc, ocean_model_flux_init
  use atmos_tracer_driver_mod,    only: atmos_tracer_flux_init

  implicit none ; private

  public :: flux_exchange_init, gas_exchange_init, &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_atmos_to_ocean,  &
     flux_ex_arrays_dealloc,&
     flux_ice_to_ocean,    &
     flux_ice_to_ocean_finish, &
     flux_ocean_to_ice,    &
     flux_ocean_to_ice_finish, &
     flux_check_stocks,    &
     flux_init_stocks,     &
     flux_ice_to_ocean_stocks,&
     flux_ocean_from_ice_stocks,&
     send_ice_mask_sic

  !-----------------------------------------------------------------------
  character(len=128) :: version = '$Id$'
  character(len=128) :: tag = '$Name$'

  logical :: do_init = .true.

  real, parameter :: bound_tol = 1e-7

  real, parameter :: d622 = rdgas/rvgas
  real, parameter :: d378 = 1.0-d622

  real :: z_ref_heat =  2. !< Reference height (meters) for temperature and relative humidity diagnostics (t_ref, rh_ref, del_h, del_q)
  real :: z_ref_mom  = 10. !< Reference height (meters) for mementum diagnostics (u_ref, v_ref, del_m)
  logical :: ex_u_star_smooth_bug = .false. !< By default, the global exchange grid \c u_star will not be interpolated
  !! from atmospheric grid, this is different from Jakarta behavior and will
  !! change answers.  So to preserve Jakarta behavior and reproduce answers
  !! explicitly set this namelist variable to .true. in input.nml.
  logical :: sw1way_bug = .false.
  logical :: do_area_weighted_flux = .FALSE.
  logical :: debug_stocks = .FALSE.
  logical :: divert_stocks_report = .FALSE.
  logical :: do_runoff = .TRUE. !< Turns on/off the land runoff interpolation to the ocean
  logical :: do_forecast = .false.
  integer :: nblocks = 1

  logical :: partition_fprec_from_lprec = .FALSE.  !< option for ATM override experiments where liquid+frozen precip are combined
  !! This option will convert liquid precip to snow when t_ref is less than
  !! tfreeze parameter
  real, parameter    :: tfreeze = 273.15
  logical :: scale_precip_2d = .false.

  namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom, ex_u_star_smooth_bug, sw1way_bug,&
       & do_area_weighted_flux, debug_stocks, divert_stocks_report, do_runoff, do_forecast, nblocks,&
       & partition_fprec_from_lprec, scale_precip_2d

  logical :: gas_fluxes_initialized = .false.  ! This is set to true when the following types are initialized.
  type(coupler_1d_bc_type), target :: ex_gas_fields_atm  ! gas fields in atm
      !< Structure containing atmospheric surfacevariables that are used in the
      !! calculation of the atmosphere-ocean gas fluxes, as well as parameters
      !! regulating these fluxes.
  type(coupler_1d_bc_type), target :: ex_gas_fields_ice  ! gas fields atop the ice or ocean
      !< Structure containing ice-top and ocean surface variables that are used
      !! in the calculation of the atmosphere-ocean gas fluxes, as well as parameters
      !! regulating these fluxes.
  type(coupler_1d_bc_type), target :: ex_gas_fluxes      ! gas fluxes between the atm and ocean
      !< A structure for exchanging gas or tracer fluxes between the atmosphere and ocean,
      !! defined by the field table, as well as a place holder of intermediate calculations,
      !! such as piston velocities, and parameters that impact the fluxes.

  integer :: ni_atm, nj_atm !< to do atmos diagnostic from flux_ocean_to_ice
  real, dimension(3) :: ccc !< for conservation checks
  !Balaji: clocks moved into flux_exchange
  integer :: cplClock

  ! Exchange grid indices
  real    :: Dt_atm, Dt_cpl
  real    :: ATM_PRECIP_NEW

contains

  !#######################################################################
  !> \brief Gas and tracer exchange initialization routine.
  !!
  !! This routine causes the field table to be read to determine which fields
  !! will be needed for the exchanges of gasses and tracers between the
  !! atmosphere and ocean.  The metadata for these fields are stored in the
  !! ex_gas_fluxes and ex_gas_fields arrays, although the data is not allocated yet.
  !! This is intended to be called (optionally) prior to flux_exchange_init.
  subroutine gas_exchange_init (gas_fields_atm, gas_fields_ice, gas_fluxes)
    type(coupler_1d_bc_type), optional, pointer :: gas_fields_atm
      !< Pointer to a structure containing atmospheric surface variables that
      !! are used in the calculation of the atmosphere-ocean gas fluxes, as well
      !! as parameters regulating these fluxes.
    type(coupler_1d_bc_type), optional, pointer :: gas_fields_ice
      !< Pointer to a structure containing ice-top and ocean surface variables
      !! that are used in the calculation of the atmosphere-ocean gas fluxes,
      !! as well as parameters regulating these fluxes.
    type(coupler_1d_bc_type), optional, pointer :: gas_fluxes
      !< Pointer to a s structure for exchanging gas or tracer fluxes between the
      !! atmosphere and ocean, defined by the field table, as well as a place holder
      !! of intermediate calculations, such as piston velocities, and parameters
      !! that impact the fluxes.

    if (.not.gas_fluxes_initialized) then
      call atmos_ocean_type_fluxes_init( )
      call ocean_model_flux_init( )
      call atmos_tracer_flux_init( )
      call atmos_ocean_fluxes_init(ex_gas_fluxes, ex_gas_fields_atm, ex_gas_fields_ice)
      gas_fluxes_initialized = .true.
    endif

    if (present(gas_fields_atm)) gas_fields_atm => ex_gas_fields_atm
    if (present(gas_fields_ice)) gas_fields_ice => ex_gas_fields_ice
    if (present(gas_fluxes)) gas_fluxes => ex_gas_fluxes

  end subroutine gas_exchange_init

  !#######################################################################
  !> \brief Initialization routine.
  !!
  !! Initializes the interpolation routines,diagnostics and boundary data
  !!
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere resolution"
  !!    The atmosphere grid size from file grid_spec.nc is not compatible with the atmosphere
  !!    resolution from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)"
  !!    The longitude from file grid_spec.nc ( from field yba ) is different from the longitude from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)"
  !!    The longitude from file grid_spec.nc ( from field xba ) is different from the longitude from atmosphere model.
  !! \throw FATAL, "grid_spec.nc incompatible with atmosphere latitudes (see grid_spec.nc)"
  !!    The latitude from file grid_spec.nc is different from the latitude from atmosphere model.
  subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, Ocean_state,&
       atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary, &
       do_ocean, slow_ice_ocean_pelist, dt_atmos, dt_cpld )

    type(time_type),                   intent(in)     :: Time !< The model's current time
    type(atmos_data_type),             intent(inout)  :: Atm !< A derived data type to specify atmosphere boundary data
    type(land_data_type),              intent(in)     :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),               intent(inout)  :: Ice !< A derived data type to specify ice boundary data
    type(ocean_public_type),           intent(inout)  :: Ocean !< A derived data type to specify ocean boundary data
    type(ocean_state_type),            pointer        :: Ocean_state
    ! All intent(OUT) derived types with pointer components must be
    ! COMPLETELY allocated here and in subroutines called from here;
    ! NO pointer components should have been allocated before entry if the
    ! derived type has intent(OUT) otherwise they may be lost.
    type(atmos_ice_boundary_type),     intent(inout) :: atmos_ice_boundary !< A derived data type to specify properties and fluxes passed from atmosphere to ice
    type(land_ice_atmos_boundary_type),intent(inout) :: land_ice_atmos_boundary !< A derived data type to specify properties and fluxes passed from exchange grid to
    !! the atmosphere, land and ice
    type(land_ice_boundary_type),      intent(inout) :: land_ice_boundary !< A derived data type to specify properties and fluxes passed from land to ice
    type(ice_ocean_boundary_type),     intent(inout) :: ice_ocean_boundary !< A derived data type to specify properties and fluxes passed from ice to ocean
    type(ocean_ice_boundary_type),     intent(inout) :: ocean_ice_boundary !< A derived data type to specify properties and fluxes passed from ocean to ice
    logical,                           intent(in)    :: do_ocean
    integer, dimension(:),             intent(in)    :: slow_ice_ocean_pelist
    integer, optional,                 intent(in)    :: dt_atmos !< Atmosphere time step in seconds
    integer, optional,                 intent(in)    :: dt_cpld !< Coupled time step in seconds

    character(len=64),  parameter   :: grid_file = 'INPUT/grid_spec.nc'
    integer        :: ierr, io
    integer        :: logunit, unit
    character(len=256) :: errmsg
    integer              :: omp_get_num_threads, nthreads

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

    call sat_vapor_pres_init()

    nthreads = 1
    ! assign nblocks to number of threads.
    !$OMP PARALLEL
    !$  nthreads = omp_get_num_threads()
    !$OMP END PARALLEL
    nblocks = nthreads

    !-----------------------------------------------------------------------
    logunit = stdlog()
    !----- read namelist -------

    read (input_nml_file, flux_exchange_nml, iostat=io)
    ierr = check_nml_error (io, 'flux_exchange_nml')

    !----- write namelist to logfile -----
    call write_version_number (version, tag)
    if( mpp_pe() == mpp_root_pe() )write( logunit, nml=flux_exchange_nml )
    if(nblocks<1) call error_mesg ('flux_exchange_mod',  &
         'flux_exchange_nml nblocks must be positive', FATAL)
    if(nblocks .NE. nthreads) then
       write(errmsg, '(a,i3,a,i3)')'flux_exchange_nml nblocks is set to ', nblocks, &
            ' is different from the default value (number of threads) = ', nthreads
       call error_mesg ('flux_exchange_mod', errmsg, NOTE)
    endif

    ! required by stock_move, all fluxes used to update stocks will be zero if dt_atmos,
    ! and dt_cpld are absent
    Dt_atm = 0.0
    Dt_cpl = 0.0
    if(present(dt_atmos)) Dt_atm = real(dt_atmos)
    if(present(dt_cpld )) Dt_cpl = real(dt_cpld)

    call get_ocean_model_area_elements(Ocean%domain, grid_file)

    if( Atm%pe )then
       call mpp_set_current_pelist(Atm%pelist)
       cplClock = mpp_clock_id( 'Land-ice-atm coupler', flags=clock_flag_default, grain=CLOCK_COMPONENT )
       call check_atm_grid(Atm, grid_file)
       call atm_land_ice_flux_exchange_init(Time, Atm, Land, Ice, atmos_ice_boundary, land_ice_atmos_boundary, &
            Dt_atm, Dt_cpl, z_ref_heat, z_ref_mom, ex_u_star_smooth_bug,  &
            sw1way_bug, do_area_weighted_flux, do_forecast,  &
            partition_fprec_from_lprec, scale_precip_2d, nblocks, cplClock, &
            ex_gas_fields_atm, ex_gas_fields_ice, ex_gas_fluxes)

       call land_ice_flux_exchange_init(Land, Ice, land_ice_boundary, Dt_cpl, do_runoff, cplClock)
    end if

    call mpp_set_current_pelist()
    call ice_ocean_flux_exchange_init(Time, Ice, Ocean, Ocean_state,ice_ocean_boundary, ocean_ice_boundary, &
         Dt_cpl, debug_stocks, do_area_weighted_flux, ex_gas_fields_ice, ex_gas_fluxes, do_ocean, slow_ice_ocean_pelist )

    !---- done ----
    do_init = .false.

  end subroutine flux_exchange_init

  !> \brief Check stock values.
  !!
  !! Will print out any difference between the integrated flux (in time
  !! and space) feeding into a component, and the stock stored in that
  !! component.

  subroutine flux_check_stocks(Time, Atm, Lnd, Ice, Ocn_state)

    type(time_type)                           :: Time
    type(atmos_data_type), optional           :: Atm
    type(land_data_type), optional            :: Lnd
    type(ice_data_type), optional             :: Ice
    type(ocean_state_type), optional, pointer :: Ocn_state

    real :: ref_value
    integer :: i


    do i = 1, NELEMS

       if(present(Atm)) then
          ref_value = 0.0
          call Atm_stock_pe(Atm, index=i, value=ref_value)
          if(i==ISTOCK_WATER .and. Atm%pe ) then
             ! decrease the Atm stock by the precip adjustment to reflect the fact that
             ! after an update_atmos_up call, the precip will be that of the future time step.
             ! Thus, the stock call will represent the (explicit ) precip at
             ! the beginning of the preceding time step, and the (implicit) evap at the
             ! end of the preceding time step
             call atm_stock_integrate(Atm, ATM_PRECIP_NEW)
             ref_value = ref_value + ATM_PRECIP_NEW
          endif

          Atm_stock(i)%q_now = ref_value
       endif

       if(present(Lnd)) then
          ref_value = 0.0
          call Lnd_stock_pe(Lnd, index=i, value=ref_value)
          Lnd_stock(i)%q_now = ref_value
       endif

       if(present(Ice)) then
          ref_value = 0.0
          call Ice_stock_pe(Ice, index=i, value=ref_value)
          Ice_stock(i)%q_now = ref_value
       endif

       if(present(Ocn_state)) then
          ref_value = 0.0
          call Ocean_stock_pe(Ocn_state, index=i, value=ref_value)
          Ocn_stock(i)%q_now = ref_value
       endif
    enddo

    call stocks_report(Time)


  end subroutine flux_check_stocks

  !#######################################################################
  !> \brief Initialize stock values.
  !!
  !! This will call the various component stock_pe routines to store the
  !! the initial stock values.

  subroutine flux_init_stocks(Time, Atm, Lnd, Ice, Ocn_state)
    type(time_type) , intent(in) :: Time
    type(atmos_data_type)        :: Atm
    type(land_data_type)         :: Lnd
    type(ice_data_type)          :: Ice
    type(ocean_state_type), pointer :: Ocn_state

    integer :: i

    stocks_file=stdout()
    ! If the divert_stocks_report is set to true, write the stocks to a new file "stocks.out"
    if(mpp_pe()==mpp_root_pe() .and. divert_stocks_report) then
       open(newunit = stocks_file, file='stocks.out', status='replace', form='formatted')
    endif

    ! Initialize stock values
    do i = 1, NELEMS
       call Atm_stock_pe(   Atm , index=i, value=Atm_stock(i)%q_start)

       if(i==ISTOCK_WATER .and. Atm%pe ) then
          call atm_stock_integrate(Atm, ATM_PRECIP_NEW)
          Atm_stock(i)%q_start = Atm_stock(i)%q_start + ATM_PRECIP_NEW
       endif

       call Lnd_stock_pe(   Lnd , index=i, value=Lnd_stock(i)%q_start)
       call Ice_stock_pe(   Ice , index=i, value=Ice_stock(i)%q_start)
       call Ocean_stock_pe( Ocn_state , index=i, value=Ocn_stock(i)%q_start)
    enddo


    call stocks_report_init(Time)


  end subroutine flux_init_stocks

  subroutine check_atm_grid(Atm, grid_file)
    type(atmos_data_type),    intent(in) :: Atm
    character(len=*),         intent(in) :: grid_file

    integer        :: isg, ieg, jsg, jeg
    integer        :: isc, iec, jsc, jec
    integer        :: isd, ied, jsd, jed
    integer        :: isc2, iec2, jsc2, jec2
    integer        :: nxg, nyg, ioff, joff
    integer        :: nlon, nlat, siz(4)
    integer        :: i, j
    type(domain2d) :: domain2
    real, dimension(:,:), allocatable :: tmpx, tmpy
    real, dimension(:),   allocatable :: atmlonb, atmlatb
    character(len=256)              :: atm_mosaic_file, tile_file, buffer

    integer, dimension(:), allocatable :: pes !> Current pelist
    type(FmsNetcdfFile_t) :: grid_file_obj, atm_mosaic_file_obj          !> Fms2io file obj
    type(FmsNetcdfDomainFile_t) :: tile_file_obj          !> Fms2io file obj
    character(len=20) :: dim_names(2) !> Array of dimension names
    integer :: ppos

    call mpp_get_global_domain(Atm%domain, isg, ieg, jsg, jeg, xsize=nxg, ysize=nyg)
    call mpp_get_compute_domain(Atm%domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Atm%domain, isd, ied, jsd, jed)

    !< Open the grid files with pelist argument so that only one pes open/reads the file
    allocate(pes(mpp_npes()))
    call mpp_get_current_pelist(pes)

    if ( .not. open_file(grid_file_obj, grid_file, "read", pelist=pes)) then
         call error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(grid_file), FATAL)
    endif

    if(size(Atm%lon_bnd,1) .NE. iec-isc+2 .OR. size(Atm%lon_bnd,2) .NE. jec-jsc+2) then
       call error_mesg ('atm_land_ice_flux_exchange_mod',  &
            'size of Atm%lon_bnd does not match the Atm computational domain', FATAL)
    endif
    ioff = lbound(Atm%lon_bnd,1) - isc
    joff = lbound(Atm%lon_bnd,2) - jsc

    if(variable_exists(grid_file_obj, "AREA_ATM" ) ) then  ! old grid
       call get_variable_size(grid_file_obj, "AREA_ATM", siz(1:2))
       nlon = siz(1)
       nlat = siz(2)

       if (nlon /= nxg .or. nlat /= nyg) then
          if (mpp_pe()==mpp_root_pe()) then
             print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                  'atmosphere has', nxg, 'longitudes,', &
                  nyg, 'latitudes (see xba.dat and yba.dat)'
          end if
          call error_mesg ('atm_land_ice_flux_exchange_mod',  &
               'grid_spec.nc incompatible with atmosphere resolution', FATAL)
       end if
       allocate( atmlonb(isg:ieg+1) )
       allocate( atmlatb(jsg:jeg+1) )
       call read_data(grid_file_obj, 'xba', atmlonb)
       call read_data(grid_file_obj, 'yba', atmlatb)

       do i=isc, iec+1
          if(abs(atmlonb(i)-Atm%lon_bnd(i+ioff,jsc+joff)*45.0/atan(1.0))>bound_tol) then
             print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ': ', &
                  atmlonb(i),  Atm%lon_bnd(i+ioff,jsc+joff)*45.0/atan(1.0)
             call error_mesg ('atm_land_ice_flux_exchange_mod', &
                  'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                  , FATAL)
          endif
       enddo
       do j=jsc, jec+1
          if(abs(atmlatb(j)-Atm%lat_bnd(isc+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
             print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at j= ',j, ': ', &
                  atmlatb(j),  Atm%lat_bnd(isc+ioff, j+joff)*45.0/atan(1.0)
             call error_mesg ('atm_land_ice_flux_exchange_mod', &
                  'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                  , FATAL)
          endif
       enddo
       deallocate(atmlonb, atmlatb)
    else if(variable_exists(grid_file_obj, "atm_mosaic_file" ) ) then  ! mosaic grid file.
       call read_data(grid_file_obj, 'atm_mosaic_file', atm_mosaic_file)

       if ( .not. open_file(atm_mosaic_file_obj, "INPUT/"//trim(atm_mosaic_file)//"", "read", pelist=pes)) then
           call error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(atm_mosaic_file), FATAL)
       endif

       call read_data(atm_mosaic_file_obj, "gridfiles", buffer, corner=1)

       !< Remove the .tile from the filename to get basename
       ppos = index(trim(buffer),".tile")
       if ( ppos > 0 ) then
          tile_file = buffer(1:ppos-1)//".nc"
       else
          tile_file = buffer
       endif

       call close_file(atm_mosaic_file_obj)

       call mpp_copy_domain(Atm%domain, domain2)
       call mpp_create_super_grid_domain(domain2)
       call mpp_define_io_domain  (domain2, (/1,1/))

       call mpp_get_compute_domain(domain2, isc2, iec2, jsc2, jec2)

       if(isc2 .NE. 2*isc-1 .OR. iec2 .NE. 2*iec+1 .OR. jsc2 .NE. 2*jsc-1 .OR. jec2 .NE. 2*jec+1) then
          call mpp_error(FATAL, 'atm_land_ice_flux_exchange_mod: supergrid domain is not set properly')
       endif

       !< This is will open the correct atm_mosaic_file for the current tile, i.e "C96_grid.tile1.nc"
       if ( .not. open_file(tile_file_obj, "INPUT/"//trim(tile_file)//"", "read", domain2)) then
          call error_mesg ('atm_land_ice_flux_exchange_mod',  &
              & 'Error opening '//trim(tile_file), FATAL)
       endif

       call get_variable_size(tile_file_obj, 'area', siz(1:2))
       nlon = siz(1); nlat = siz(2)
       if( mod(nlon,2) .NE. 0) call mpp_error(FATAL,  &
            'atm_land_ice_flux_exchange_mod: atmos supergrid longitude size can not be divided by 2')
       if( mod(nlat,2) .NE. 0) call mpp_error(FATAL,  &
            'atm_land_ice_flux_exchange_mod: atmos supergrid latitude size can not be divided by 2')
       nlon = nlon/2
       nlat = nlat/2
       if (nlon /= nxg .or. nlat /= nyg) then
          if (mpp_pe()==mpp_root_pe()) then
             print *, 'atmosphere mosaic tile has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                  'atmosphere has', nxg, 'longitudes,', nyg, 'latitudes'
          end if
          call error_mesg ('atm_land_ice_flux_exchange_mod',  &
               'atmosphere mosaic tile grid file incompatible with atmosphere resolution', FATAL)
       end if

       allocate(tmpx(isc2:iec2,jsc2:jec2), tmpy(isc2:iec2,jsc2:jec2) )

       !< Register the dimension of the variables "x" and "y" in the atm_mosaic_file
       call get_variable_dimension_names(tile_file_obj, "x", dim_names)
       call register_axis(tile_file_obj, dim_names(1), "x")
       call register_axis(tile_file_obj, dim_names(2), "y")
       call register_field(tile_file_obj, "x", "double", dim_names)
       call register_field(tile_file_obj, "y", "double", dim_names)

       !< Read the variables "x" and "y" as domain decomposed variables from the atm_moasic_file
       call read_data( tile_file_obj, 'x', tmpx)
       call read_data( tile_file_obj, 'y', tmpy)

       call close_file(tile_file_obj)

       call mpp_deallocate_domain(domain2)

       do j = jsc, jec+1
          do i = isc, iec+1
             if (abs(tmpx(2*i-1,2*j-1)-Atm%lon_bnd(i+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                     tmpx(2*i-1,2*j-1),  Atm%lon_bnd(i+ioff,j+joff)*45.0/atan(1.0)
                call error_mesg ('atm_land_ice_flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere longitudes (see '//trim(tile_file)//')'&
                     ,FATAL)
             end if
             if (abs(tmpy(2*i-1,2*j-1)-Atm%lat_bnd(i+ioff,j+joff)*45.0/atan(1.0))>bound_tol) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY at i= ',i, ', j= ', j, ': ', &
                     tmpy(2*i-1,2*j-1),  Atm%lat_bnd(i+ioff,j+joff)*45.0/atan(1.0)
                call error_mesg ('atm_land_ice_flux_exchange_mod', &
                     'grid_spec.nc incompatible with atmosphere latitudes (see '//trim(tile_file)//')'&
                     ,FATAL)
             end if
          end do
       end do
       deallocate(tmpx, tmpy)
    else
       call mpp_error(FATAL, 'atm_land_ice_flux_exchange_mod: both AREA_ATMxOCN and ocn_mosaic_file does not exist in '//trim(grid_file))
    end if

    call close_file(grid_file_obj)
  end subroutine check_atm_grid

end module flux_exchange_mod
