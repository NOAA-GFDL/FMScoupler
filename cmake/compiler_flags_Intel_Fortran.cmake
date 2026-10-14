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
# Intel Fortran flags, from mkmf/templates/linux-intel.mk

# Default is the template's SSE setting
set(INTEL_ISA "-march=core-avx-i -qno-opt-dynamic-align" CACHE STRING "Intel instruction set (ISA) compile option")

set(COUPLER_FFLAGS "-fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -i4 -real-size 64 -nowarn -g -sox -traceback")

set(COUPLER_FFLAGS_RELEASE "-O2 ${INTEL_ISA}")
set(COUPLER_FFLAGS_REPRO "-fp-model source -O2 ${INTEL_ISA}")
set(COUPLER_FFLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -debug variable_locations -fpe0 -ftrapuv ${INTEL_ISA}")

set(COUPLER_FFLAGS_OPENMP          "-qopenmp")
set(COUPLER_FFLAGS_OVERRIDE_LIMITS "-qoverride-limits")
set(COUPLER_FFLAGS_VERBOSE         "-v -V -what -warn all")

set(COUPLER_LDFLAGS         "")
set(COUPLER_LDFLAGS_OPENMP  "-qopenmp")
set(COUPLER_LDFLAGS_VERBOSE "-Wl,-V,--verbose,-cref,-M")
