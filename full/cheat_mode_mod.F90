!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @defgroup cheat_mode_mod cheat_mode_mod
!> @ingroup cheat_mode_mod
!> @brief Generate or extract a tarball containing the final contents of a run directory.
!!
!> This module supports "cheat mode", so that the workflow can be tested without the
!! computational cost of an actual model run. A directory is supplied via the
!! `cheat_mode_nml` namelist: if a tarball is found in this directory corresponding
!! to the date range of the current segment, `cheat_mode_invoke` will extract this
!! tarball; otherwise, it will generate it from the content of the run directory.

#ifdef CHEAT_MODE

module cheat_mode_mod
  use FMS
  implicit none
  private

  public :: cheat_mode_tarball_exists, cheat_mode_tarball_path, cheat_mode_init, cheat_mode_invoke

  logical :: cheat_mode_tarball_exists !< Whether to extract (true) or generate (false) a tarball
  character(:), allocatable :: cheat_mode_tarball_path !< Path to the tarball to be extracted or created
  character(:), allocatable :: dir !< Directory containing the tarballs

  namelist /cheat_mode_nml/ dir

!> @addtogroup cheat_mode_mod
!> @{

contains

  !> Initialize cheat mode. Determine the path to the tarball from `cheat_mode_nml` and
  !! from the date range of the current segment.
  subroutine cheat_mode_init(year0, month0, day0, year1, month1, day1)
    integer, intent(in) :: year0, month0, day0 !< Initial year, month, and day of the current segment
    integer, intent(in) :: year1, month1, day1 !< Final year, month, and day of the current segment
    integer :: io_status

    read (fms_input_nml_file, cheat_mode_nml, iostat=io_status)
    io_status = fms_check_nml_error(io_status, "cheat_mode_nml")

    allocate (character(len(dir) + 22) :: cheat_mode_tarball_path)
    write (cheat_mode_tarball_path, '(A,"/",I0.4,I0.2,I0.2,"_",I0.4,I0.2,I0.2,".tar")') &
          dir, year0, month0, day0, year1, month1, day1
    inquire (file=cheat_mode_tarball_path, exist=cheat_mode_tarball_exists)
  end subroutine

  !> Invoke the tar command to either create or extract the tarball. This subroutine
  !! should only be called by the root PE.
  subroutine cheat_mode_invoke
    if (cheat_mode_tarball_exists) then
      call execute_command_line("tar -xf " // cheat_mode_tarball_path // " --touch --overwrite")
    else
      call execute_command_line("tar -cf " // cheat_mode_tarball_path // " `find . -type f -newer input.nml | xargs`")
    endif
  end subroutine cheat_mode_invoke
end module cheat_mode_mod

#endif

!> @}
! close documentation grouping
