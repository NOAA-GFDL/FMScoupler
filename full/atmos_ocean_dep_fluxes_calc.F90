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
module atmos_ocean_dep_fluxes_calc_mod

  use FMS

  implicit none

  character(len=*), parameter :: mod_name = "aodfc"
contains

  !> \brief atmos_ocean_dep_fluxes_calc
  !
  !! \throw FATAL, "Number of gas fluxes not zero"
  !! \throw FATAL, "atmos_ocean_dep_fluxes_calc: Bad parameter ([gas_fluxes%bc(n)%param(1)]) for air_sea_deposition for [gas_fluxes%bc(n)%name]"
  !! \throw FATAL, "atmos_ocean_dep_fluxes_calc: Unknown implementation ([gas_fluxes%bc(n)%implementation] for [gas_fluxes%bc(n)%name]"
  subroutine atmos_ocean_dep_fluxes_calc(gas_fields_atm, gas_fields_ice, gas_fluxes, seawater)
    type(coupler_1d_bc_type), intent(in)    :: gas_fields_atm !< Structure containing atmospheric surface
                                                              !! variables that are used in the calculation
                                                              !! of the atmosphere-ocean gas fluxes.
    type(coupler_1d_bc_type), intent(in)    :: gas_fields_ice !< Structure containing ice-top and ocean
                                                              !! surface variables that are used in the
                                                              !! calculation of the atmosphere-ocean gas fluxes.
    type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< Structure containing the gas fluxes between
                                                          !! the atmosphere and the ocean and parameters
                                                          !! related to the calculation of these fluxes.
    real, dimension(:),       intent(in)    :: seawater   !< 1 for the open water category, 0 if ice or land.

    character(len=64), parameter    :: sub_name = 'atmos_ocean_dep_fluxes_calc'
    character(len=256), parameter   :: error_header = &
        &'==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

    integer                                 :: n
    integer                                 :: i
    integer                                 :: length
    real, dimension(:), allocatable         :: kw
    real, dimension(:), allocatable         :: cair
    character(len=128)                      :: error_string

    real, parameter :: permeg=1.0e-6

    ! Return if no fluxes to be calculated

    if (gas_fluxes%num_bcs .le. 0) return

    if (.not. associated(gas_fluxes%bc)) then
      if (gas_fluxes%num_bcs .ne. 0) then
        call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
      else
        return
      endif
    endif

    do n = 1, gas_fluxes%num_bcs
      ! only do calculations if the flux has not been overridden
      if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then
        if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then
          if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
            write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
            call mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Bad parameter (' //&
                & trim(error_string) // ') for air_sea_deposition for ' //&
                & trim(gas_fluxes%bc(n)%name))
          endif

          length = size(gas_fluxes%bc(n)%field(1)%values(:))

          if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = &
                    gas_fields_atm%bc(n)%field(ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              endif
            enddo
          elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
            do i = 1, length
              if (seawater(i) == 1.) then
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = &
                    gas_fields_atm%bc(n)%field(ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
              else
                gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              endif
            enddo
          else
            call mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Unknown implementation ('&
                & // trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
          endif
        else
          cycle
        endif
      endif
    enddo
  end subroutine  atmos_ocean_dep_fluxes_calc
end module atmos_ocean_dep_fluxes_calc_mod
