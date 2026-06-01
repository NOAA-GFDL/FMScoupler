#!/usr/bin/env sh
#
# Script to build a GFDL null model, using all null components, and run
# a simple test on CI systems, like Travis CI or gitlab CI.

script_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
bld_dir=$(mktemp --directory $script_root/null.XXXXXX)
cd $bld_dir
cmake -DCMAKE_BUILD_TYPE=Debug ../..
make -j
make

# Report on the status of the build
if [ $? -eq 0 ]
then
  echo "::note title=Build Succeeded:: null model with simple coupler built successfully."
else
  echo "::error title=Build Failed:: null model with simple coupler failed compilation."
  exit 1
fi

# Run the null models test
# Setup the run directory
mkdir ${bld_dir}/run
cd ${bld_dir}/run
mkdir RESTART
# Get the data files required for the run
tarFile=coupler_null_test_data_full_simple.tar.gz
curl -O ftp://ftp.gfdl.noaa.gov/perm/GFDL_pubrelease/test_data/${tarFile}
tar zxf ${tarFile}

# add an io layout to the full nml
sed -i '22i  io_layout = 1, 1' input-full.nml

# Get the full namelist
ln -s input-full.nml input.nml
# Run the null model with the full coupler
mpiexec -n 1 ${bld_dir}/coupler_full.x

# Report on the status of the run with the full coupler
if [ $? -eq 0 ]
then
  echo "::note title=Run Succeeded - full coupler:: Full coupler null model ran successfully."
else
  echo "::error title=Run Failed - full coupler:: Full coupler null model run failed execution."
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
mpiexec -n 1 ${bld_dir}/coupler_simple.x

# Report on the status of the run with the simple coupler
if [ $? -eq 0 ]
then
  echo "::note title=Run Succeeded - simple coupler:: simple coupler null model ran successfully"
else
  echo "::error title=Run Failed - simple coupler:: simple coupler null model run failed execution."
  exit 1
fi
