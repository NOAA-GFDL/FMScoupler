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
# GNU Fortran flags, from mkmf/templates/linux-gnu.mk

set(COUPLER_FFLAGS "-fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -fno-range-check")
# gfortran >= 10 rejects legacy MPI argument mismatches (as in ncrc5-gcc.mk)
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  string(APPEND COUPLER_FFLAGS " -fallow-argument-mismatch")
endif()

set(COUPLER_FFLAGS_RELEASE "-O3")
set(COUPLER_FFLAGS_REPRO "-O2 -fbounds-check")
set(COUPLER_FFLAGS_DEBUG "-O0 -g -W -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow")

set(COUPLER_FFLAGS_OPENMP          "-fopenmp")
set(COUPLER_FFLAGS_OVERRIDE_LIMITS "")
set(COUPLER_FFLAGS_VERBOSE         "-Wall -Wextra")

set(COUPLER_LDFLAGS         "")
set(COUPLER_LDFLAGS_OPENMP  "-fopenmp")
set(COUPLER_LDFLAGS_VERBOSE "")
