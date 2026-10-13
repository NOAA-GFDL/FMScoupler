#!/bin/bash
# Standalone ESM4 cmake build using a slightly modified version of the ESM4.5 public release
# creates a new directory in t/ subdirectory
# called esm4, clones source code and performs a build with cmake

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workdir="$SCRIPT_DIR/esm4"
mkdir -p "$workdir"
srcdir="$workdir/src"
builddir="$workdir"

# this is just cloning the same source code as the esm4.5 public release branch on github
if [[ ! -d "$srcdir" ]]; then
  mkdir -p "$srcdir"
  cd "$srcdir"
  # symlink this repos source code so its with the rest
  ln -s ${SCRIPT_DIR}/.. FMScoupler

  git clone https://github.com/NOAA-GFDL/FMS.git FMS
  git -C FMS checkout 342afcc5080423ced4a7d9f255df12ce37386af5

  git clone -b esm4-cmake https://github.com/NOAA-GFDL/FMScoupler.git FMScoupler

  git clone https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere.git GFDL_atmos_cubed_sphere
  git -C GFDL_atmos_cubed_sphere checkout 4bd9f120332d719e143298ca22776e3b6245637d

  git clone https://github.com/NOAA-GFDL/atmos_drivers.git atmos_drivers
  git -C atmos_drivers checkout 8181e97f03c1605a51123aad85b6195df146e3ba

  git clone https://github.com/NOAA-GFDL/atmos_phys.git atmos_phys
  git -C atmos_phys checkout ea9e1d924e7493b3b2ef09d3c89316faedf6c731

  git clone https://github.com/NOAA-GFDL/MOM6.git MOM6
  git -C MOM6 checkout 3d479de4d4c3528633bb450454de4ea7b98e5070
  git -C MOM6 submodule update --init --recursive

  git clone https://github.com/NOAA-GFDL/SIS2.git SIS2
  git -C SIS2 checkout f96e21736008769c758b85ab56709f69c170a8ba

  git clone https://github.com/NOAA-GFDL/icebergs.git icebergs
  git -C icebergs checkout 765b4cab2d82ac2d740b946b61f9d3dc02d733f0

  git clone https://github.com/NOAA-GFDL/ice_param.git ice_param
  git -C ice_param checkout 1d19a57e0be2ba1b668e2f0dd7a7290733b21cfe

  git clone https://github.com/NOAA-GFDL/ocean_BGC.git ocean_BGC
  git -C ocean_BGC checkout 5760648c59983b66aa44591b1a6a0133f852f333
  git -C ocean_BGC submodule update --init --recursive

  git clone https://github.com/NOAA-GFDL/lm4.git lm4p
  git -C lm4p checkout f0fce68dc1c834e086ce7d557e049dd8830c02b8
  # this little edit is only needed when compiling with gfortran
  # which requires all strings in array constructors to be the same size
  sed -i 's/rainfed_Soybean          /rainfed_Soybean         /' lm4p/vegetation/vegn_data.F90
  sed -i 's/second_season/second_season /' lm4p/vegetation/vegn_data.F90

  cd "$workdir"
else
  echo "$srcdir already exists, skipping checkout"
fi

paths_fms="$srcdir/FMS/"
macros_fms="use_libMPI;use_netCDF;use_yaml;MAXFIELDMETHODS_=600;MAXXGRID=1e9"

paths_ocean=( "$srcdir/MOM6/"{config_src/infra/FMS2,config_src/memory/dynamic_nonsymmetric,config_src/drivers/FMS_cap,config_src/external/ODA_hooks,config_src/external/database_comms,config_src/external/stochastic_physics,config_src/external/MARBL,config_src/external/drifters,pkg/GSW-Fortran/{modules,toolbox}/,src/{*,*/*}/} "$srcdir/FMS/"{coupler,include} "$srcdir/"{ocean_BGC/generic_tracers,ocean_BGC/mocsy/src} )
macros_ocean="MAX_FIELDS_=600;NOT_SET_AFFINITY;_USE_MOM6_DIAG;_USE_GENERIC_TRACER;USE_PRECISION=2"
flags_ocean="-I$srcdir/FMS;-I$srcdir/MOM6;-I$srcdir/MOM6/src/framework/;-I$srcdir/ocean_BGC/mocsy/src"

paths_ice=( "$srcdir/SIS2/"{config_src/dynamic,config_src/external/Icepack_interfaces,src} "$srcdir/icebergs/src" "$srcdir/ice_param" )
macros_ice="USE_FMS2_IO"
flags_ice="-I$srcdir/FMS;-I$srcdir/MOM6;-I$srcdir/MOM6/src/framework/"

paths_land=( "$srcdir/lm4p" )
macros_land=""
flags_land="-I$srcdir/FMS;-I$srcdir/FMS/include"

paths_atmos=( "$srcdir/GFDL_atmos_cubed_sphere/"{driver/GFDL,model,GFDL_tools,tools} "$srcdir/atmos_drivers/coupled" )
macros_atmos="CLIMATE_NUDGE;SPMD"
flags_atmos="-I$srcdir/FMS;-I$srcdir/FMS/include;-I$srcdir/GFDL_atmos_cubed_sphere/tools/"

paths_atmos_phys=( "$srcdir/atmos_phys" )
macros_atmos_phys=""
flags_atmos_phys="-I$srcdir/FMS;-I$srcdir/FMS/include;-I$srcdir/atmos_phys/atmos_param/cosp"

build_type="${BUILD_TYPE:-Debug}"
generator="${CMAKE_GENERATOR:-Unix Makefiles}"

IFS=";"

cmake -S "$srcdir/FMScoupler" -B "$builddir" \
  -DCOUPLER_TYPE="full" \
  -DEXECUTABLE_NAME="esm45.x" \
  -DCMAKE_BUILD_TYPE="$build_type" \
  -G "$generator" \
  -DPREPROCESS_LAND=on \
  -Dfms_paths="$paths_fms" \
  -Dfms_macros="${macros_fms}" \
  -Docean_model_paths="${paths_ocean[*]}" \
  -Docean_model_macros="$macros_ocean" \
  -Docean_model_flags="$flags_ocean" \
  -Datmos_model_paths="${paths_atmos[*]}" \
  -Datmos_model_macros="$macros_atmos" \
  -Datmos_model_flags="$flags_atmos" \
  -Dland_model_paths="${paths_land[*]}" \
  -Dland_model_macros="$macros_land" \
  -Dland_model_flags="$flags_land" \
  -Dice_model_paths="${paths_ice[*]}" \
  -Dice_model_macros="$macros_ice" \
  -Dice_model_flags="$flags_ice" \
  -Datmos_physics_library_paths="${paths_atmos_phys[*]}" \
  -Datmos_physics_library_macros="$macros_atmos_phys" \
  -Datmos_physics_library_flags="$flags_atmos_phys"

unset IFS
make -C "$builddir" -j
