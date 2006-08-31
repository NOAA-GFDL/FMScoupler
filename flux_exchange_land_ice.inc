! This include file contains the subroutine
! flux_land_to_ice
!
! This subroutine makes extensive use of exchange grid global variables
! defined in the main flux_exchange module

!#######################################################################
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! flux_land_to_ice - translate runoff from land to ice grids                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! <SUBROUTINE NAME="flux_land_to_ice">
!  <OVERVIEW>
!   Conservative transfer of water and snow discharge from the land model to sea ice/ocean model. 
!  </OVERVIEW>
!  <DESCRIPTION>
!  <PRE>
!    The following elements are transferred from the Land to the Land_ice_boundary: 
!
!        discharge --> runoff (kg/m2)
!        discharge_snow --> calving (kg/m2)
!
!  </PRE>
!  </DESCRIPTION>
!  <TEMPLATE>
!   call flux_land_to_ice(Time, Land, Ice, Land_Ice_Boundary )
!		
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Land" TYPE="land_data_type">
!   A derived data type to specify land boundary data.
!  </IN>
!  <IN NAME="Ice" TYPE="ice_data_type">
!   A derived data type to specify ice boundary data.
!  </IN>
!  <INOUT NAME="Land_Ice_Boundary" TYPE="land_ice_boundary_type">
!   A derived data type to specify properties and fluxes passed from land to ice.
!  </INOUT>
!
subroutine flux_land_to_ice( Time, Land, Ice, Land_Ice_Boundary )
  type(time_type),               intent(in) :: Time
  type(land_data_type),          intent(in) :: Land
  type(ice_data_type),           intent(in) :: Ice
!real, dimension(:,:), intent(out) :: runoff_ice, calving_ice
  type(land_ice_boundary_type),  intent(inout):: Land_Ice_Boundary
  
  real, dimension(n_xgrid_runoff) :: ex_runoff, ex_calving
  real, dimension(size(Land_Ice_Boundary%runoff,1),size(Land_Ice_Boundary%runoff,2),1) :: ice_buf
!Balaji
  call mpp_clock_begin(cplClock)
  call mpp_clock_begin(fluxLandIceClock)

  ! ccc = conservation_check(Land%discharge, 'LND', xmap_runoff)
  ! if (mpp_pe()==mpp_root_pe()) print *,'RUNOFF', ccc

  call put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
  call put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
  call get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
  Land_Ice_Boundary%runoff = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
  Land_Ice_Boundary%calving = ice_buf(:,:,1);
!Balaji
  call data_override('ICE', 'runoff' , Land_Ice_Boundary%runoff , Time)
  call data_override('ICE', 'calving', Land_Ice_Boundary%calving, Time)
  call mpp_clock_end(fluxLandIceClock)
  call mpp_clock_end(cplClock)

end subroutine flux_land_to_ice
! </SUBROUTINE>
