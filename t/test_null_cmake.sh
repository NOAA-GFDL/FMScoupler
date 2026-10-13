#!/bin/bash
# Standalone null cmake build, creates a new directory in t/ subdirectory
# called null_model, clones source code and performs a build with cmake
# Usage:
# ./null_cmake_build.sh
# 	[COUPLER_TYPE="full|simple"]
#   [BUILD_TYPE="Release/Prod|Debug|Repro"]
#   [GENERATOR="Unix Makefiles|Ninja|etc.."]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workdir="$SCRIPT_DIR/null_model"
mkdir -p "$workdir"
srcdir="$workdir/src"
builddir="$workdir"

build_type="${BUILD_TYPE:-Debug}"
generator="${CMAKE_GENERATOR:-Unix Makefiles}"

# clone null components
if [[ ! -d "$srcdir" ]]; then
  mkdir -p "$srcdir"
  cd "$srcdir"
  # symlink this repos source code so its with the rest
  ln -s ${SCRIPT_DIR}/.. FMScoupler

  git clone https://github.com/NOAA-GFDL/FMS
  git clone https://github.com/NOAA-GFDL/ocean_null
  git clone https://github.com/NOAA-GFDL/atmos_null
  git clone https://github.com/NOAA-GFDL/land_null
  git clone https://github.com/NOAA-GFDL/ice_param
  git clone https://github.com/NOAA-GFDL/ice_null
else
  echo "$srcdir already exists, skipping checkout"
fi

fms_paths=( "${srcdir}/FMS" )
fms_macros="use_netCDF;use_libMPI;MAXFIELDS_=200;MAXFIELDMETHODS_=200;INTERNAL_FILE_NML;HAVE_GETTID"

ocean_paths=( "${srcdir}/ocean_null" )
ocean_flags="-I${srcdir}/FMS/include"
ocean_macros=""

atmos_paths=( "${srcdir}/atmos_null" )
atmos_flags="-I${srcdir}/FMS/include"
atmos_macros="INTERNAL_FILE_NML"

ice_paths=( "${srcdir}/ice_param" "${srcdir}/ice_null" )
ice_flags="-I${srcdir}/FMS/include"
ice_macros="INTERNAL_FILE_NML"

land_paths=( "${srcdir}/land_null" )
land_flags="-I${srcdir}/FMS/include"
land_macros="INTERNAL_FILE_NML"

coupler_macros="_USE_LEGACY_LAND_;use_AM3_physics;INTERNAL_FILE_NML"

coupler_type="${COUPLER_TYPE:-full}"
build_type="${BUILD_TYPE:-Debug}"
generator="${CMAKE_GENERATOR:-Unix Makefiles}"

# cmake lists are semicolon-separated; this makes "${array[*]}" join with semicolons
IFS=";"

cmake -S "$srcdir/FMScoupler" -B "$builddir" \
  -DCOUPLER_TYPE="$coupler_type" \
  -DEXECUTABLE_NAME="null_model.x" \
  -DCMAKE_BUILD_TYPE="$build_type" \
  -G "$generator" \
  -DTEST_NULL_MODEL=on \
  -Dfms_paths="$fms_paths" \
  -Docean_model_paths="${ocean_paths[*]}" \
  -Docean_model_macros="$ocean_macros" \
  -Docean_model_flags="$ocean_flags" \
  -Datmos_model_paths="${atmos_paths[*]}" \
  -Datmos_model_macros="$atmos_macros" \
  -Datmos_model_flags="$atmos_flags" \
  -Dland_model_paths="${land_paths[*]}" \
  -Dland_model_macros="$land_macros" \
  -Dland_model_flags="$land_flags" \
  -Dice_model_paths="${ice_paths[*]}" \
  -Dice_model_macros="$ice_macros" \
  -Dice_model_flags="$ice_flags" \
  -Dcoupler_macros="$coupler_macros"

unset IFS
make -C "$builddir" -j
make -C "$builddir" test
