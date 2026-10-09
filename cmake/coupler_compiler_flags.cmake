#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
#*
#* FMS Coupler is free software: you can redistribute it and/or modify
#* it under the terms of the GNU Lesser General Public License as
#* published by the Free Software Foundation, either version 3 of the
#* License, or (at your option) any later version.
#*
#* FMS Coupler is distributed in the hope that it will be useful, but
#* WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#* General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS Coupler.
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
################################################################################
# Compiler flags for the FMS coupler, FMS, and all components, taken from the
# mkmf templates (mkmf/templates/linux-*.mk).
#
# CMAKE_BUILD_TYPE maps to the mkmf build targets:
#   Release -> FFLAGS_OPT   (default)
#   Repro   -> FFLAGS_REPRO (REPRO=on)
#   Debug   -> FFLAGS_DEBUG (DEBUG=on)
################################################################################

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE_UPPER STREQUAL "PROD")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Type of build, corresponding to the compiler optimization and options flags when building. Options are 'Release'/'Prod', 'Repro', or 'Debug'."
    FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Repro Debug)
# OPENMP is shared with the FMS build
option(OPENMP          "Build with OpenMP support"                             ON)
option(OVERRIDE_LIMITS "Compile with the compiler's override-limits flag, if any" ON)
option(VERBOSE_FLAGS   "Add the mkmf template's verbose compile/link flags"    OFF)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_flags_GNU_Fortran.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_flags_Intel_Fortran.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_flags_IntelLLVM_Fortran.cmake)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  include(${CMAKE_CURRENT_LIST_DIR}/compiler_flags_NVHPC_Fortran.cmake)
else()
  message(FATAL_ERROR "Unsupported Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'. Supported compilers are GNU, Intel, IntelLLVM, and NVHPC.")
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" _coupler_build_type)
if(NOT _coupler_build_type MATCHES "^(RELEASE|REPRO|DEBUG)$")
  message(FATAL_ERROR "Invalid CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}'. Expected 'Release', 'Repro', or 'Debug'.")
endif()

set(_coupler_fflags "${COUPLER_FFLAGS} ${COUPLER_FFLAGS_${_coupler_build_type}}")
set(_coupler_ldflags "${COUPLER_LDFLAGS}")

if(OPENMP)
  string(APPEND _coupler_fflags " ${COUPLER_FFLAGS_OPENMP}")
  string(APPEND _coupler_ldflags " ${COUPLER_LDFLAGS_OPENMP}")
endif()
if(OVERRIDE_LIMITS)
  string(APPEND _coupler_fflags " ${COUPLER_FFLAGS_OVERRIDE_LIMITS}")
endif()
if(VERBOSE_FLAGS)
  string(APPEND _coupler_fflags " ${COUPLER_FFLAGS_VERBOSE}")
  string(APPEND _coupler_ldflags " ${COUPLER_LDFLAGS_VERBOSE}")
endif()
string(STRIP "${_coupler_fflags}" _coupler_fflags)
string(STRIP "${_coupler_ldflags}" _coupler_ldflags)

# Directory-level options are inherited by subdirectories (FMS) and come after CMAKE_Fortran_FLAGS[_<CONFIG>],
# so they take precedence. SHELL: keeps repeated words like "-check" from being de-duplicated.
if(_coupler_fflags)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:${_coupler_fflags}>")
endif()
if(_coupler_ldflags)
  add_link_options("SHELL:${_coupler_ldflags}")
endif()

# Optimization/debug flags come from the templates, so drop CMake's build-type defaults
# (FMS sets its own for Release/Debug, but ours are applied after them)
foreach(_type RELEASE REPRO DEBUG RELWITHDEBINFO MINSIZEREL)
  set(CMAKE_Fortran_FLAGS_${_type} "")
endforeach()

message(STATUS "Coupler build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Coupler Fortran flags: ${_coupler_fflags}")
message(STATUS "Coupler linker flags: ${_coupler_ldflags}")
