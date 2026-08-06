#!/usr/bin/env sh
#
# Run script for a simple test run using the null_model
# This is the same run done in null_model_build.sh, just without any building steps
set -eu
rundir=$(pwd -P)

# Run the null models test
# Setup the run directory
rm -rf "${rundir}/run"
mkdir -p "${rundir}/run/RESTART"
cd "${rundir}/run"
# Get the data files required for the run
tarFile=coupler_null_test_data_full_simple.tar.gz
curl -O ftp://ftp.gfdl.noaa.gov/perm/GFDL_pubrelease/test_data/${tarFile}
tar zxf ${tarFile}
# add an io layout to the full nml
awk 'NR==22{print "  io_layout = 1, 1"}1' input-full.nml > input-full.nml.tmp && mv input-full.nml.tmp input-full.nml
# Get the full namelist
ln -s input-full.nml input.nml
# Run the null model with the full coupler
mpiexec -n 1 ${rundir}/coupler_full.x

# Report on the status of the run with the full coupler
if [ $? -eq 0 ]
then
  echo "::note title=Run Succeeded - full coupler:: Full coupler null model ran successfully."
else
  echo "::error title=Run Failed - full coupler:: Full coupler null model run failed execution."
  exit 1
fi