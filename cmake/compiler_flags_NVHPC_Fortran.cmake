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
# NVHPC Fortran flags, from mkmf/templates/linux-nvhpc.mk

set(COUPLER_FFLAGS "-i4 -r8 -byteswapio -Mflushz -Mdaz -D_F2000")

set(COUPLER_FFLAGS_RELEASE "-O3 -Mallocatable=95 -ffp-contract=off")
set(COUPLER_FFLAGS_REPRO "-O2 -fno-associative-math -ffp-contract=off")
set(COUPLER_FFLAGS_DEBUG "-O0 -g -Ktrap=fp")

set(COUPLER_FFLAGS_OPENMP          "-mp")
set(COUPLER_FFLAGS_OVERRIDE_LIMITS "")
set(COUPLER_FFLAGS_VERBOSE         "-v -Minform=inform")

set(COUPLER_LDFLAGS         "-byteswapio")
# The template leaves this blank, but nvfortran needs -mp at link time for OpenMP objects
set(COUPLER_LDFLAGS_OPENMP  "-mp")
set(COUPLER_LDFLAGS_VERBOSE "-v")
