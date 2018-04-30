!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module coupler_flux_tracers_mod
  !> \brief Go as Richard Slater
  !!
  !! \author Richard Slater <Richard.Slater@noaa.gov

  use mpp_mod,           only: stdout, mpp_error, FATAL
  use fms_mod,           only: write_version_number
  use fm_util_mod,       only: fm_util_set_value, fm_util_set_no_overwrite
  use fm_util_mod,       only: fm_util_set_caller, fm_util_reset_no_overwrite
  use fm_util_mod,       only: fm_util_reset_caller
  use field_manager_mod, only: fm_new_list, fm_change_list, fm_dump_list
  use coupler_types_mod, only: ind_alpha, ind_csurf, ind_sc_no
  use coupler_types_mod, only: ind_pcair, ind_u10, ind_psurf
  use coupler_types_mod, only: ind_deposition
  use coupler_types_mod, only: ind_runoff
  use coupler_types_mod, only: ind_flux, ind_deltap, ind_kw


  implicit none

  private

  ! Include variable "version" to be written to log file.
#include<file_version.h>

  public  coupler_flux_tracers_init
  character(len=*), parameter :: MOD_NAME = 'coupler_flux_tracers_mod'

contains

  !#######################################################################
  !> \brief Initialize the coupler flux tracers
  !!
  !! \throw FATAL, "Could not set the \"coupler_mod\" list"
  !! \throw FATAL, "Could not set the \"GOOD\" list"
  !! \throw FATAL, "Could not set the \"/coupler_mod/fluxes\" list"
  !! \throw FATAL, "Could not set the \"/coupler_mod/types\" list"
  !! \throw FATAL, "Could not change to \"/coupler_mod/types\""
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic/implementation\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic/implementation/ocmip2\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic/atm\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic/ice\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux_generic/flux\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/implementation\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/implementation/ocmip2\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/implementation/ocmip2_data\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/implementation/linear\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/atm\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/ice\" list"
  !! \throw FATAL, "Could not set the \"air_sea_gas_flux/flux\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/implementation\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/implementation/dry\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/implementation/wet\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/atm\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/ice\" list"
  !! \throw FATAL, "Could not set the \"air_sea_deposition/flux\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff/implementation\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff/implementation/river\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff/atm\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff/ice\" list"
  !! \throw FATAL, "Could not set the \"land_sea_runoff/flux\" list"
  !! \throw FATAL, "Could not change to \"/\""
  !! \throw FATAL, "Problem dumping /coupler_mod/types tree"
  subroutine coupler_flux_tracers_init()
    character(len=64), parameter    :: SUB_NAME = 'coupler_flux_tracers_init'
    character(len=256), parameter   :: ERROR_HEADER =&
         & '==>Error from ' // trim(MOD_NAME) // '(' // trim(SUB_NAME) // '):'

    integer            :: outunit
    character(len=128) :: error_msg
    logical, save   :: module_is_initialized = .false.

    if (.NOT.module_is_initialized) then
       ! Set other defaults for the fm_util_set_value routines.
       call fm_util_set_no_overwrite(.true.)
       call fm_util_set_caller(sub_name)
       
       ! Be sure that the various lists and fields are defined in the field manager tree.
       if (fm_new_list('/coupler_mod') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "coupler_mod" list')
       endif

       if (fm_new_list('/coupler_mod/GOOD') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "GOOD" list')
       endif
       call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'GOOD', append = .true.)

       if (fm_new_list('/coupler_mod/fluxes') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "/coupler_mod/fluxes" list')
       endif
       call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'fluxes', append = .true.)

       if (fm_new_list('/coupler_mod/types') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "/coupler_mod/types" list')
       endif
       call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'types', append = .true.)

       ! Change to the "/coupler_mod/types" list.
       if (.not. fm_change_list('/coupler_mod/types')) then
          call mpp_error(FATAL, trim(error_header) // ' Could not change to "/coupler_mod/types"')
       endif

       ! Define the air_sea_gas_flux_generic type and add it.
       if (fm_new_list('air_sea_gas_flux_generic') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic" list')
       endif

       ! Add the implementation list.
       if (fm_new_list('air_sea_gas_flux_generic/implementation') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/implementation" list')
       endif

       ! Add the names of the different implementations.
       if (fm_new_list('air_sea_gas_flux_generic/implementation/ocmip2') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/implementation/ocmip2" list')
       endif
       call fm_util_set_value('air_sea_gas_flux_generic/implementation/ocmip2/num_parameters', 2)

       ! Add some scalar quantaties.
       call fm_util_set_value('air_sea_gas_flux_generic/num_flags', 0)
       call fm_util_set_value('air_sea_gas_flux_generic/use_atm_pressure', .true.)
       call fm_util_set_value('air_sea_gas_flux_generic/use_10m_wind_speed', .true.)
       call fm_util_set_value('air_sea_gas_flux_generic/pass_through_ice', .false.)

       ! Add required fields that will come from the atmosphere model.
       if (fm_new_list('air_sea_gas_flux_generic/atm') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/atm" list')
       endif

       call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'pcair',                     index = ind_pcair)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Atmospheric concentration', index = ind_pcair)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'mol/mol',                   index = ind_pcair)
       
       call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'u10',                index = ind_u10)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Wind speed at 10 m', index = ind_u10)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'm/s',                index = ind_u10)
       
       call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'psurf',                        index = ind_psurf)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Surface atmospheric pressure', index = ind_psurf)
       call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'Pa',                           index = ind_psurf)

       ! Add required fields that will come from the ice model.
       if (fm_new_list('air_sea_gas_flux_generic/ice') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/ice" list')
       endif

       call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'alpha',                        index = ind_alpha)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Solubility w.r.t. atmosphere', index = ind_alpha)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'mol/m^3/atm',                  index = ind_alpha)

       call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'csurf',               index = ind_csurf)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Ocean concentration', index = ind_csurf)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'mol/m^3',             index = ind_csurf)

       call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'sc_no',          index = ind_sc_no)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Schmidt number', index = ind_sc_no)
       call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'dimensionless',  index = ind_sc_no)

       ! Add the flux output field(s).
       if (fm_new_list('air_sea_gas_flux_generic/flux') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/flux" list')
       endif

       call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'flux',         index = ind_flux)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Surface flux', index = ind_flux)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'mol/m^2/s',    index = ind_flux)

       call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'deltap',         index = ind_deltap)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Ocean-air delta pressure', index = ind_deltap)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'uatm',    index = ind_deltap)

       call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'kw',         index = ind_kw)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Piston velocity', index = ind_kw)
       call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'm/s',    index = ind_kw)

       ! Define the air_sea_gas_flux type and add it.
       if (fm_new_list('air_sea_gas_flux') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux" list')
       endif

       ! Add the implementation list.
       if (fm_new_list('air_sea_gas_flux/implementation') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation" list')
       endif

       ! Add the names of the different implementations.
       if (fm_new_list('air_sea_gas_flux/implementation/ocmip2') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/ocmip2" list')
       endif
       call fm_util_set_value('air_sea_gas_flux/implementation/ocmip2/num_parameters', 2)
       if (fm_new_list('air_sea_gas_flux/implementation/ocmip2_data') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/ocmip2_data" list')
       endif
       call fm_util_set_value('air_sea_gas_flux/implementation/ocmip2_data/num_parameters', 2)
       if (fm_new_list('air_sea_gas_flux/implementation/linear') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/linear" list')
       endif
       call fm_util_set_value('air_sea_gas_flux/implementation/linear/num_parameters', 3)

       ! Add some scalar quantaties.
       call fm_util_set_value('air_sea_gas_flux/num_flags', 0)
       call fm_util_set_value('air_sea_gas_flux/use_atm_pressure', .true.)
       call fm_util_set_value('air_sea_gas_flux/use_10m_wind_speed', .true.)
       call fm_util_set_value('air_sea_gas_flux/pass_through_ice', .false.)

       ! Add required fields that will come from the atmosphere model.
       if (fm_new_list('air_sea_gas_flux/atm') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/atm" list')
       endif

       call fm_util_set_value('air_sea_gas_flux/atm/name',      'pcair',                     index = ind_pcair)
       call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Atmospheric concentration', index = ind_pcair)
       call fm_util_set_value('air_sea_gas_flux/atm/units',     'mol/mol',                   index = ind_pcair)

       call fm_util_set_value('air_sea_gas_flux/atm/name',      'u10',                index = ind_u10)
       call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Wind speed at 10 m', index = ind_u10)
       call fm_util_set_value('air_sea_gas_flux/atm/units',     'm/s',                index = ind_u10)

       call fm_util_set_value('air_sea_gas_flux/atm/name',      'psurf',                        index = ind_psurf)
       call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Surface atmospheric pressure', index = ind_psurf)
       call fm_util_set_value('air_sea_gas_flux/atm/units',     'Pa',                           index = ind_psurf)

       ! Add required fields that will come from the ice model.
       if (fm_new_list('air_sea_gas_flux/ice') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/ice" list')
       endif

       call fm_util_set_value('air_sea_gas_flux/ice/name',      'alpha',                                                index = ind_alpha)
       call fm_util_set_value('air_sea_gas_flux/ice/long_name', 'Solubility from atmosphere times Schmidt number term', index = ind_alpha)
       call fm_util_set_value('air_sea_gas_flux/ice/units',     'mol/m^3/atm',                                          index = ind_alpha)
       
       call fm_util_set_value('air_sea_gas_flux/ice/name',      'csurf',                                         index = ind_csurf)
       call fm_util_set_value('air_sea_gas_flux/ice/long_name', 'Ocean concentration times Schmidt number term', index = ind_csurf)
       call fm_util_set_value('air_sea_gas_flux/ice/units',     'mol/m^3',                                       index = ind_csurf)

       ! Add the flux output field(s).
       if (fm_new_list('air_sea_gas_flux/flux') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/flux" list')
       endif

       call fm_util_set_value('air_sea_gas_flux/flux/name',      'flux',         index = ind_flux)
       call fm_util_set_value('air_sea_gas_flux/flux/long_name', 'Surface flux', index = ind_flux)
       call fm_util_set_value('air_sea_gas_flux/flux/units',     'mol/m^2/s',    index = ind_flux)

       ! Define the air_sea_deposition type and add it.
       if (fm_new_list('air_sea_deposition') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition" list')
       endif

       ! Add the implementation list.
       if (fm_new_list('air_sea_deposition/implementation') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation" list')
       endif

       ! Add the names of the different implementations.
       if (fm_new_list('air_sea_deposition/implementation/dry') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation/dry" list')
       endif
       call fm_util_set_value('air_sea_deposition/implementation/dry/num_parameters', 1)
       if (fm_new_list('air_sea_deposition/implementation/wet') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation/wet" list')
       endif
       call fm_util_set_value('air_sea_deposition/implementation/wet/num_parameters', 1)

       ! Add some scalar quantaties.
       call fm_util_set_value('air_sea_deposition/num_flags', 0)
       call fm_util_set_value('air_sea_deposition/use_atm_pressure', .false.)
       call fm_util_set_value('air_sea_deposition/use_10m_wind_speed', .false.)
       call fm_util_set_value('air_sea_deposition/pass_through_ice', .true.)

       ! Add required fields that will come from the atmosphere model.
       if (fm_new_list('air_sea_deposition/atm') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/atm" list')
       endif

       call fm_util_set_value('air_sea_deposition/atm/name',      'deposition',             index = ind_deposition)
       call fm_util_set_value('air_sea_deposition/atm/long_name', 'Atmospheric deposition', index = ind_deposition)
       call fm_util_set_value('air_sea_deposition/atm/units',     'kg/m^2/s',               index = ind_deposition)

       ! Add required fields that will come from the ice model.
       if (fm_new_list('air_sea_deposition/ice') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/ice" list')
       endif

       call fm_util_set_value('air_sea_deposition/ice/name',      ' ', index = 0)
       call fm_util_set_value('air_sea_deposition/ice/long_name', ' ', index = 0)
       call fm_util_set_value('air_sea_deposition/ice/units',     ' ', index = 0)

       ! Add the flux output field(s).
       if (fm_new_list('air_sea_deposition/flux') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/flux" list')
       endif

       call fm_util_set_value('air_sea_deposition/flux/name',      'flux',               index = ind_flux)
       call fm_util_set_value('air_sea_deposition/flux/long_name', 'Surface deposition', index = ind_flux)
       call fm_util_set_value('air_sea_deposition/flux/units',     'mol/m^2/s',          index = ind_flux)

       ! Define the land_sea_runoff type and add it.
       if (fm_new_list('land_sea_runoff') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff" list')
       endif

       ! Add the implementation list.
       if (fm_new_list('land_sea_runoff/implementation') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/implementation" list')
       endif

       ! Add the names of the different implementations.
       if (fm_new_list('land_sea_runoff/implementation/river') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/implementation/river" list')
       endif
       call fm_util_set_value('land_sea_runoff/implementation/river/num_parameters', 1)

       ! Add some scalar quantaties.
       call fm_util_set_value('land_sea_runoff/num_flags', 0)
       call fm_util_set_value('land_sea_runoff/use_atm_pressure', .false.)
       call fm_util_set_value('land_sea_runoff/use_10m_wind_speed', .false.)
       call fm_util_set_value('land_sea_runoff/pass_through_ice', .true.)

       ! Add required fields that will come from the land model (the array name is still called "atm").
       if (fm_new_list('land_sea_runoff/atm') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/atm" list')
       endif

       call fm_util_set_value('land_sea_runoff/atm/name',      'runoff',                       index = ind_runoff)
       call fm_util_set_value('land_sea_runoff/atm/long_name', 'Concentration in land runoff', index = ind_runoff)
       call fm_util_set_value('land_sea_runoff/atm/units',     'mol/m^3',                      index = ind_runoff)

       ! Add required fields that will come from the ice model.
       if (fm_new_list('land_sea_runoff/ice') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/ice" list')
       endif

       call fm_util_set_value('land_sea_runoff/ice/name',      ' ', index = 0)
       call fm_util_set_value('land_sea_runoff/ice/long_name', ' ', index = 0)
       call fm_util_set_value('land_sea_runoff/ice/units',     ' ', index = 0)

       ! Add the flux output field(s).
       if (fm_new_list('land_sea_runoff/flux') .le. 0) then
          call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/flux" list')
       endif

       call fm_util_set_value('land_sea_runoff/flux/name',      'flux',                         index = ind_flux)
       call fm_util_set_value('land_sea_runoff/flux/long_name', 'Concentration in land runoff', index = ind_flux)
       call fm_util_set_value('land_sea_runoff/flux/units',     'mol/m^3',                      index = ind_flux)

       ! Change back to root list.
       if (.not. fm_change_list('/')) then
          call mpp_error(FATAL, trim(error_header) // ' Could not change to "/"')
       endif

       ! Reset the defaults for the fm_util_set_value calls.
       call fm_util_reset_no_overwrite
       call fm_util_reset_caller

       module_is_initialized = .true.
       
       ! Dump the coupler_mod types list.
       outunit = stdout()
       write (outunit,*)
       write (outunit,*) 'Dumping coupler_mod/types tree'
       if (.not. fm_dump_list('/coupler_mod/types', recursive = .true.)) then
          call mpp_error(FATAL, trim(error_header) // ' Problem dumping /coupler_mod/types tree')
       endif
    endif
    return
  end subroutine coupler_flux_tracers_init
end module coupler_flux_tracers_mod
