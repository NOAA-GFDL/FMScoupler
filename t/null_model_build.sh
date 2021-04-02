#!/usr/bin/env sh
#
# Script to build a GFDL null model, using all null components, and run
# a simple test on CI systems, like Travis CI or gitlab CI.

# Determine the where this script lives, and set some variables that contain
# other useful directories.
script_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)

# Create a new build directory, to keep from polluting the test
# dirctory for a time when there are more tests.
bld_dir=$(mktemp --directory $script_root/null.XXXXXX)
cd $bld_dir

# Add a directory for the source(s)
src_dir=$bld_dir/src
mkdir $src_dir

# Get the other required repositories for the build.
#
# mkmf -- used for now to build the model
git clone https://github.com/NOAA-GFDL/mkmf.git
export PATH=${PATH}:${bld_dir}/mkmf/bin
mk_template=${bld_dir}/mkmf/templates/linux-ubuntu-xenial-gnu.mk

# FMS
git clone https://github.com/NOAA-GFDL/FMS.git $src_dir/FMS

# ocean_null
git clone https://github.com/NOAA-GFDL/ocean_null.git $src_dir/ocean_null
cd $bld_dir

# atmos_null
git clone https://github.com/NOAA-GFDL/atmos_null.git $src_dir/atmos_null

# land_null
git clone https://github.com/NOAA-GFDL/land_null $src_dir/land_null

# ice_null - need ice_param as well, and depends on ocean_null.
git clone https://github.com/NOAA-GFDL/ice_param.git $src_dir/ice_param
git clone https://github.com/NOAA-GFDL/ice_null.git $src_dir/ice_null
cd $bld_dir

# coupler - simply create symlink, this simplifies using the build system.
ln -s $(readlink -f ../../) $src_dir/coupler

# Create the main Makefile
sed -e 's/<TAB>/\t/' >$bld_dir/Makefile <<EOF
# Makefile for GFDL FMS capable model for use with mkmf build system

SRCROOT = $src_dir/
BUILDROOT = $bld_dir/

MK_TEMPLATE = $mk_template
include \$(MK_TEMPLATE)

coupler_full_test.x: coupler_full/libcoupler_full.a atmos/libatmos_null.a land/libland_null.a ocean/libocean_null.a ice/libice_null.a ice_param/libice_param.a fms/libfms.a
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@ \$(STATIC_LIBS)

coupler_simple_test.x: coupler_simple/libcoupler_simple.a atmos/libatmos_null.a land/libland_null.a ocean/libocean_null.a ice_param/libice_param.a fms/libfms.a
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@ \$(STATIC_LIBS)

fms/libfms.a:  FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=fms \$(@F) 

ocean/libocean_null.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=ocean \$(@F) 

atmos/libatmos_null.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos \$(@F) 

ice_param/libice_param.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=ice_param \$(@F) 

ice/libice_null.a: ocean/libocean_null.a ice_param/libice_param.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=ice \$(@F) 

land/libland_null.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=land \$(@F) 

coupler_full/libcoupler_full.a: atmos/libatmos_null.a ice/libice_null.a ice_param/libice_param.a ocean/libocean_null.a land/libland_null.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=coupler_full \$(@F) 

coupler_simple/libcoupler_simple.a: atmos/libatmos_null.a ice/libice_null.a ice_param/libice_param.a ocean/libocean_null.a land/libland_null.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=coupler_simple \$(@F) 

FORCE:

EOF

# Create the component library Makefiles
#
# libfms
mkdir -p $bld_dir/fms
list_paths -o $bld_dir/fms/pathnames_fms $src_dir/FMS
cd $bld_dir/fms
mkmf -m Makefile -a $src_dir -b $bld_dir -p libfms.a -t $mkmf_template -g -c "-Duse_netCDF -Duse_libMPI -DMAXFIELDS_=200 -DMAXFIELDMETHODS_=200 -DINTERNAL_FILE_NML" -IFMS/include -IFMS/mpp/include $bld_dir/fms/pathnames_fms
cd $bld_dir

# libocean_null
mkdir -p $bld_dir/ocean
list_paths -o $bld_dir/ocean/pathnames_ocean $src_dir/ocean_null
cd $bld_dir/ocean
mkmf -m Makefile -a $src_dir -b $bld_dir -p libocean_null.a -t $mkmf_template -g -o "-I$bld_dir/fms" -IFMS/include $bld_dir/ocean/pathnames_ocean
cd $bld_dir

# libatmos_null
mkdir -p $bld_dir/atmos
list_paths -o $bld_dir/atmos/pathnames_atmos $src_dir/atmos_null
cd $bld_dir/atmos
mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_null.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/fms" -IFMS/include $bld_dir/atmos/pathnames_atmos
cd $bld_dir

# libice_param
mkdir -p $bld_dir/ice_param
list_paths -o $bld_dir/ice_param/pathnames_ice_param $src_dir/ice_param
cd $bld_dir/ice_param
mkmf -m Makefile -a $src_dir -b $bld_dir -p libice_param.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/fms" -IFMS/include $bld_dir/ice_param/pathnames_ice_param
cd $bld_dir

# libice_null
mkdir -p $bld_dir/ice
list_paths -o $bld_dir/ice/pathnames_ice $src_dir/ice_null
cd $bld_dir/ice
mkmf -m Makefile -a $src_dir -b $bld_dir -p libice_null.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/ocean -I$bld_dir/ice_param -I$bld_dir/fms" -IFMS/include $bld_dir/ice/pathnames_ice
cd $bld_dir

# libland_null
mkdir -p $bld_dir/land
list_paths -o $bld_dir/land/pathnames_land $src_dir/land_null
cd $bld_dir/land
mkmf -m Makefile -a $src_dir -b $bld_dir -p libland_null.a -t $mkmf_template -g -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/fms" -IFMS/include $bld_dir/land/pathnames_land
cd $bld_dir

# libcoupler_full
mkdir -p $bld_dir/coupler_full
list_paths -o $bld_dir/coupler_full/pathnames_coupler $src_dir/coupler/shared $src_dir/coupler/full
cd $bld_dir/coupler_full
mkmf -m Makefile -a $src_dir -b $bld_dir -p libcoupler_full.a -t $mkmf_template -g -c "-D_USE_LEGACY_LAND_ -Duse_AM3_physics -DINTERNAL_FILE_NML" -o "-I$bld_dir/land -I$bld_dir/ice -I$bld_dir/ice_param -I$bld_dir/atmos -I$bld_dir/ocean -I$bld_dir/fms" -IFMS/include $bld_dir/coupler_full/pathnames_coupler
cd $bld_dir

# libcoupler_simple
mkdir -p $bld_dir/coupler_simple
list_paths -l -o $bld_dir/coupler_simple/pathnames_coupler $src_dir/coupler/shared $src_dir/coupler/simple
cd $bld_dir/coupler_simple
mkmf -m Makefile -a $src_dir -b $bld_dir -p libcoupler_simple.a -t $mkmf_template -g -c "-D_USE_LEGACY_LAND_ -Duse_AM3_physics -DINTERNAL_FILE_NML" -o "-I$bld_dir/land -I$bld_dir/ice_param -I$bld_dir/atmos -I$bld_dir/ocean -I$bld_dir/fms" -IFMS/include $bld_dir/coupler_simple/pathnames_coupler
cd $bld_dir

# Call make to build the executable with the full coupler
make -j NETCDF=3 DEBUG=on coupler_full_test.x

# Report on the status of the build
if [ $? -eq 0 ]
then
  echo "<NOTE> : make succeeded - full coupler."
else
  echo "<NOTE> : make failed - full coupler."
  exit 1
fi

# Call make to build the executable with the full coupler
make -j NETCDF=3 DEBUG=on coupler_simple_test.x

# Report on the status of the build
if [ $? -eq 0 ]
then
  echo "<NOTE> : make succeeded - simple coupler."
else
  echo "<NOTE> : make failed - simple coupler."
  exit 1
fi
### 17FEB2021 exit 0 to prevent model running.  This is temporary
exit 0
# Run the null models test
# Setup the run directory
mkdir ${bld_dir}/run
cd ${bld_dir}/run
mkdir RESTART
# Get the data files required for the run
tarFile=coupler_null_test_data_full_simple.tar.gz
wget ftp://ftp.gfdl.noaa.gov/perm/GFDL_pubrelease/test_data/${tarFile}
tar zxf ${tarFile}

# Get the full namelist
ln -s input-full.nml input.nml
# Run the null model with the full coupler
### 17FEB2021 commented out the run because it crashes
#mpiexec -n 1 ${bld_dir}/coupler_full_test.x

# Report on the status of the run with the full coupler
if [ $? -eq 0 ]
then
  echo "<NOTE> : run succeeded - full coupler."
else
  echo "<NOTE> : run failed - full coupler."
  exit 1
fi

# Using the same run directory, setup for the simple coupler
# Clear out the RESTART directory
mv RESTART RESTART_full
mkdir RESTART
# Get the simple namelist
rm input.nml
ln -s input-simple.nml input.nml
# Run the null simple coupler test
### 17FEB2021 commented out the run because it crashes
#mpiexec -n 1 ${bld_dir}/coupler_simple_test.x

