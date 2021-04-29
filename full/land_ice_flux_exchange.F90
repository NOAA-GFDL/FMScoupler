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

module land_ice_flux_exchange_mod

!! FMS
  use FMS
  use FMSconstants, only: RADIUS
!! Components
  use land_model_mod,      only: land_data_type
  use ice_model_mod,       only: ice_data_type, land_ice_boundary_type

  implicit none
  private


  !---- exchange grid maps -----

  type(xmap_type), save :: xmap_runoff
  integer         :: n_xgrid_runoff=0

  ! Exchange grid indices
  integer :: X2_GRID_LND, X2_GRID_ICE

  public :: flux_land_to_ice, land_ice_flux_exchange_init

  integer :: cplClock, fluxLandIceClock
  logical :: do_runoff
  real    :: Dt_cpl
contains

  subroutine land_ice_flux_exchange_init(Land, Ice, land_ice_boundary, Dt_cpl_in, do_runoff_in, cplClock_in)
    type(land_data_type),         intent(in)    :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),          intent(inout) :: Ice !< A derived data type to specify ice boundary data
    type(land_ice_boundary_type), intent(inout) :: land_ice_boundary !< A derived data type to specify properties
                                                                     !! and fluxes passed from land to ice
    real,                         intent(in)    :: Dt_cpl_in
    logical,                      intent(in)    :: do_runoff_in
    integer,                      intent(in)    :: cplClock_in

    integer :: is, ie, js, je

    do_runoff = do_runoff_in
    cplClock = cplClock_in
    Dt_cpl   = Dt_cpl_in
    fluxLandIceClock = mpp_clock_id( 'Flux land to ice', flags=clock_flag_default, grain=CLOCK_ROUTINE )

    if (do_runoff) then
       call setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
            (/ Land%Domain, Ice%Domain /),                    &
            "INPUT/grid_spec.nc"             )
       ! exchange grid indices
       X2_GRID_LND = 1; X2_GRID_ICE = 2;
       n_xgrid_runoff = max(xgrid_count(xmap_runoff),1)
       if (n_xgrid_runoff.eq.1) write (*,'(a,i6,6x,a)') 'PE = ', mpp_pe(), 'Runoff  exchange size equals one.'
    endif

    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( land_ice_boundary%runoff(is:ie,js:je) )
    allocate( land_ice_boundary%calving(is:ie,js:je) )
    allocate( land_ice_boundary%runoff_hflx(is:ie,js:je) )
    allocate( land_ice_boundary%calving_hflx(is:ie,js:je) )
    ! initialize values for override experiments (mjh)
    land_ice_boundary%runoff=0.0
    land_ice_boundary%calving=0.0
    land_ice_boundary%runoff_hflx=0.0
    land_ice_boundary%calving_hflx=0.0


  end subroutine land_ice_flux_exchange_init

  !#######################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! flux_land_to_ice - translate runoff from land to ice grids                   !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !> \brief Conservative transfer of water and snow discharge from the land model to sea ice/ocean model.
  !!
  !! The following elements are transferred from the Land to the Land_ice_boundary:
  !! <pre>
  !!        discharge --> runoff (kg/m2)
  !!        discharge_snow --> calving (kg/m2)
  !! </pre>
  subroutine flux_land_to_ice( Time, Land, Ice, Land_Ice_Boundary )
    type(time_type),                intent(in) :: Time !< Current time
    type(land_data_type),           intent(in) :: Land !< A derived data type to specify land boundary data
    type(ice_data_type),            intent(in) :: Ice !< A derived data type to specify ice boundary data
    !real, dimension(:,:),         intent(out) :: runoff_ice, calving_ice
    type(land_ice_boundary_type), intent(inout):: Land_Ice_Boundary !< A derived data type to specify properties and fluxes passed
                                                                    !! from land to ice

    integer                         :: ier
    real, dimension(n_xgrid_runoff) :: ex_runoff, ex_calving, ex_runoff_hflx, ex_calving_hflx
    real, dimension(size(Land_Ice_Boundary%runoff,1),size(Land_Ice_Boundary%runoff,2),1) :: ice_buf

    !Balaji
    call mpp_clock_begin(cplClock)
    call mpp_clock_begin(fluxLandIceClock)

    ! ccc = conservation_check(Land%discharge, 'LND', xmap_runoff)
    ! if (mpp_pe()==mpp_root_pe()) print *,'RUNOFF', ccc

    if (do_runoff) then
       call put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
       call put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
       call put_to_xgrid ( Land%discharge_heat,      'LND', ex_runoff_hflx,  xmap_runoff)
       call put_to_xgrid ( Land%discharge_snow_heat, 'LND', ex_calving_hflx, xmap_runoff)
       call get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
       Land_Ice_Boundary%runoff = ice_buf(:,:,1);
       call get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
       Land_Ice_Boundary%calving = ice_buf(:,:,1);
       call get_from_xgrid (ice_buf, 'OCN', ex_runoff_hflx,  xmap_runoff)
       Land_Ice_Boundary%runoff_hflx = ice_buf(:,:,1);
       call get_from_xgrid (ice_buf, 'OCN', ex_calving_hflx, xmap_runoff)
       Land_Ice_Boundary%calving_hflx = ice_buf(:,:,1);
       !Balaji
       call data_override('ICE', 'runoff' , Land_Ice_Boundary%runoff , Time)
       call data_override('ICE', 'calving', Land_Ice_Boundary%calving, Time)
       call data_override('ICE', 'runoff_hflx' , Land_Ice_Boundary%runoff_hflx , Time)
       call data_override('ICE', 'calving_hflx', Land_Ice_Boundary%calving_hflx, Time)

       ! compute stock increment
       ice_buf(:,:,1) = Land_Ice_Boundary%runoff + Land_Ice_Boundary%calving
       call stock_move(from=Lnd_stock(ISTOCK_WATER), to=Ice_stock(ISTOCK_WATER), &
            & grid_index=X2_GRID_ICE, &
            & data=ice_buf, &
            & xmap=xmap_runoff, &
            & delta_t=Dt_cpl, &
            & from_side=ISTOCK_SIDE, to_side=ISTOCK_SIDE, &
            & radius=Radius, ier=ier, verbose='stock move RUNOFF+CALVING (Lnd->Ice) ')
    else
       Land_Ice_Boundary%runoff = 0.0
       Land_Ice_Boundary%calving = 0.0
       Land_Ice_Boundary%runoff_hflx = 0.0
       Land_Ice_Boundary%calving_hflx = 0.0
    endif

    call mpp_clock_end(fluxLandIceClock)
    call mpp_clock_end(cplClock)

  end subroutine flux_land_to_ice


!#######################################################################

end module land_ice_flux_exchange_mod
