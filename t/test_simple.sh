#!/usr/bin/env sh
#
# Run script for a simple test run using the null_model
# This is the same run done in null_model_build.sh, just without any building steps
rundir=$PWD

# Run the null models test
# Setup the run directory
rm -rf run
mkdir ${rundir}/run
cd ${rundir}/run
mkdir RESTART
# Get the data files required for the run
tarFile=coupler_null_test_data_full_simple.tar.gz
curl -O ftp://ftp.gfdl.noaa.gov/perm/GFDL_pubrelease/test_data/${tarFile}
tar zxf ${tarFile}
# Get the simple namelist
ln -s input-simple.nml input.nml
# Run the null simple coupler test
mpiexec -n 1 ${rundir}/coupler_simple.x

# Report on the status of the run with the simple coupler
if [ $? -eq 0 ]
then
  echo "::note title=Run Succeeded - simple coupler:: simple coupler null model ran successfully"
else
  echo "::error title=Run Failed - simple coupler:: simple coupler null model run failed execution."
  exit 1
fi
