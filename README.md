![Build Status](https://github.com/NOAA-GFDL/FMScoupler/actions/workflows/build.yml/badge.svg?branch=main)

# Flexible Modeling System (FMS) Coupler

Instructions for contribuing to this project can be found in the
[Contribution Guidelines](CONTRIBUTING.md).  Contributing to the FMScoupler involves a higher
level of scrutiny and justification.


The Flexible Modeling System (FMS) is a software framework for supporting the
efficient development, construction, execution, and scientific interpretation
of atmospheric, oceanic, and climate system models.

This repository contains the FMS Coupler, meant for use alongside FMS and other component models,
such as the MOM6 ocean model and the AM4 atmosphere model. More technical information for coupling
models is available below.

For additional information on FMS, please see the [FMS github repository.](http://www.github.com/noaa-gfdl/fms)


## Coupling Models with FMS Coupler

FMS Coupler provides the capability to couple component models (atmosphere, land, sea ice, and ocean
) on different logically rectangular grids. Presently, the atmosphere and land models are
constrained to be on standard longitude/latitude grids (not necessarily the same) and the ocean and
sea-ice are constrained to be on the same grid which need only be logically rectangular.
A logically rectangular grid has an array-like set of areas, each of which border one other area in
each of the pseudo-north, -south, -east, and -west directions, except at the poles.
The coupling between the models is designed to conserve fluxes. For coupled models, a grid
specification file is used to initialize the model grids and perform exchanges between the models.
The next sections describe how this file and associated grids are used in the coupler.

### Grid Specification Files
At runtime, the coupled model sets up its grid using a given `grid_spec.nc` file that it reads from
the INPUT subdirectory. This NetCDF file contains grid information for all of the component models
as well as exchange grid information for the coupler to use.

For more information on generating the necessary files for coupled runs, please see the
[FRE-NCtools repository](http://github.com/noaa-gfdl/FRE-NCtools) for information on grid file
generation tools.

In the coupled model, ocean, sea ice, and land models will read their grids from the grid_spec.nc
file. The land model areas must be read from grid_spec.nc rather than calculated because the
grid_spec.nc land areas have been modified in the grid generation process to remove overlaps with
the ocean/sea-ice grid cells. The land mask is set to true where this modified area is positive.
The land model areas are used for conserving runoff on the land grid. The ocean and sea-ice grids
are the same by virtue of initializing from the same fields of grid_spec.nc. Additionally, the sea
ice model uses this grid information to rotate vectors between the longitude/latitude atmosphere
grid and the general ocean/sea-ice grids.

### Exchange Grids

The coupler uses grid_spec.nc to initialize its exchange grids. An exchange grid between two
component model grids is the grid formed with the union of the bounding lines of the component
model grids. The exchange grid is, therefore, the coarsest grid that is a refinement of each of the
component model grids. The coupler uses exchange grids for two purposes:

- conservative interpolation of fields between models uses the exchange grid cell areas as weights and
- the surface flux calculation takes place on the exchange grid thereby using the finest scale data available.

The coupler has two exchange grids. The first is for surface fluxes with the atmosphere on one side
and the land and sea ice on the other. Under FMS the sea ice model serves as the interface to the
ocean model — the atmosphere model never exchanges directly with the ocean model. The second
exchange grid is between the land and the sea ice for runoff. No fluxes are computed on this
exchange grid; it is used solely for conservation.

The coupler’s utility for interfacing to the grid_spec.nc file and performing exchange grid
operations is xgrid_mod (from the FMS repository). Xgrid_mod uses the mpp_domains domain of each of
the models along with information it reads from the grid specification file to determine grid and
processor connectivities. The coupler’s fortran call to initialize the surface exchange grid
(xmap_sfc) is:
```
call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),                 &
                          (/ Atm%Domain, Ice%Domain, Land%Domain /), &
                          "INPUT/grid_spec.nc"                       )
```
Xgrid_mod reads the exchange grids from grid_spec.nc as a sequence of quintuples: the i/j indices of
the intersecting cells of the two participating grids and their areal overlap. The names of the five
fields are generated automatically from the three character ids of the participating grids that
appear in the above initialization call. For example, for atmosphere/sea ice exchange on the
coupler’s surface exchange grid, the following fields are read by xgrid_mod: I_ATM_ATMxOCN, J_ATM_ATMxOCN, I_OCN_ATMxOCN, J_OCN_ATMxOCN, and AREA_ATMxOCN. These fields were placed in grid_spec.nc by the make_xgrids utility.

### Source Code Documentation

A publicly available version of this repository's doxygen-generated source documentation is
available at https://noaa-gfdl.github.io/FMScoupler. It is updated upon every release.

Alternatively, documentation can be generated manually by running `doxygen` in the top level of your
cloned repository.


### Code Structure

This repository contains two executables, drivers for coupled climate models, each with its own `coupler_main.F90`:

- `full/` - the full coupler, for fully coupled models (atmosphere, land, ice, and ocean), including
  the flux exchange modules between all of the components.
- `simple/` - a simplified coupler for simpler models, with its own `flux_exchange.F90`
  and a basic `ice_model.F90`.

Code used by both executables lives in `shared/` (e.g. `surface_flux.F90`). Each executable also has
modules that exist only in its own directory (e.g. `full/flux_exchange.F90` and
`simple/flux_exchange.F90` are separate implementations).

The FMScoupler repository does not contain the component libraries which are required to compile (ie. `ocean_model_mod`, `atmosphere_model_mod`, etc..).
FMS and the atmosphere, land, ice, and ocean components (plus the atmosphere physics library) are provided by separate repositories or
groups of repositories, which vary between models. See the [CMake Build](#cmake-build) section for details.

### Build/Test Scripts

The coupler has been historically compiled using the [mkmf](https://github.com/noaa-gfdl/mkmf) utility in combination with FRE to generate build scripts,
but now a CMake build is also available.

To test compilation and functionality, some standalone build scripts are provided in the `t/` subdirectory:

- `t/test_null_mkmf.sh` clones the null components, builds, and runs a simple initialization test for both full and simple couplers
- `t/test_null_cmake.sh` performs the same tasks but using cmake and ctest
- `t/esm4_build_cmake.sh` builds a slightly modified version of the [ESM4.5 public release](https://github.com/noaa-gfdl/esm4)

The other scripts are used by ctest to do the runtime testing:

- `t/test_simple.sh` performs the same simple coupler run test as `test_null_mkmf.sh`
- `t/test_full.sh` performs the same full coupler run test as `test_null_mkmf.sh`

### CMake Build

When compiling with cmake, it is best to set the environment variables `CC`/`FC` to the compiler executables you intend to use. CMake automagically detects compiler vendor
and sets flags based on what is available, so this can lead to environment issues if the wrong compiler is detected.

Options for compiler flags/FMS:
- `-DCMAKE_BUILD_TYPE="Release|Debug|Repro|Prod"` Determines the compilation flags based on detected compiler (default Release/Prod)
- `-DWITH_YAML` Enables support for yaml IO in FMS

Options for the coupler's code:
- `-DOPENMP`                                  Enable openmp directive support in the coupler and components (TODO...)
- `-DBUILD_FULL_COUPLER`                      build the full coupler executable (default on)
- `-DBUILD_SIMPLE_COUPLER`                    build the simple coupler executable (default off)
- `-DFULL_COUPLER_EXECUTABLE_NAME="string"`   name of the full coupler executable (default `coupler_full.x`)
- `-DSIMPLE_COUPLER_EXECUTABLE_NAME="string"` name of the simple coupler executable (default `coupler_simple.x`)
- `-DPREPROCESS_LAND`                         enables preprocessing the land model source code with `cpp` explicitly. Preprocessed source code are written to `preprocessed_srcs/` in the build directory. This option also excludes `test_*.F90` files from being collected in the land model source.
- `-DOVERRIDE_LIMITS`                         Enable overriding internal compiler limits, only applicable for intel compilers (default off)
- `-DVERBOSE_FLAGS`                           Enable compiler verbosity flags for extra output when compiling
- `-DTEST_NULL_MODEL`                         enables adding null model runscripts for both couplers to ctest. These can be run using `make test` or `ctest`.
- `-Dcoupler_macros`                          Preprocessor macro definitions for the coupler source code

Options for components, set for each dependency (ie. `land_model`, `atmosphere_physics_library`, etc..)
- `-D<component_name>_paths`  Paths for the given component (required)
- `-D<component_name>_macros` Preprocessor definitions when compiling given component
- `-D<component_name>_flags`  Compiler flags when compiling given component

The FMScoupler CMake build has been set up to compile not only the coupler but also component libraries it depends on. It has been modeled after the mkmf tool, which uses a path-based
process for compiling libraries (also called models) that the coupler depends on. This should allow for easily converting existing mkmf build scripts.

For the coupler, the 6 direct dependencies are:
- FMS
- Atmosphere model
- Land model
- Ocean Model
- Ice Model
- Atmosphere physics library (optional when using `atmos_null` as the atmosphere model)

For each of these dependencies, the CMake build sets the following from user input:
- source code paths (which are searched recursively)
- preprocessor macro definitions
- any other compiler flags

The source code paths may be a single directory containing a repository, or a longer list of specific folders and files that are compiled together to form the library.

For the non-FMS components, the list of paths provided are searched recursively, then compiled together as a library and named according to the dependency (eg. `libocean_model.a`).
For FMS, one path is expected, and the whole repository will be built using it's CMake build and any flags/macros provided.
All the build options for the FMS cmake build are inherited and should be used when configuring the coupler if needed.


Here's an example from the ESM4.5 build with the argument values filled in:
```{bash}
cmake -S src/FMScoupler -B esm4 \
    -DBUILD_FULL_COUPLER=on \
    -DFULL_COUPLER_EXECUTABLE_NAME=esm45.x \
    -DCMAKE_BUILD_TYPE=Debug \
    -G "Unix Makefiles" \
    -DPREPROCESS_LAND=on \
    -Dfms_paths=src/FMS/ \
    -Dfms_macros=use_libMPI;use_netCDF;use_yaml;MAXFIELDMETHODS_=600;MAXXGRID=1e9 \
    -Docean_model_paths=src/MOM6/config_src/infra/FMS2;src/MOM6/config_src/memory/dynamic_nonsymmetric;src/MOM6/config_src/drivers/FMS_cap;src/MOM6/config_src/external/ODA_hooks;src/MOM6/config_src/external/database_comms;src/MOM6/config_src/external/stochastic_physics;src/MOM6/config_src/external/MARBL;src/MOM6/config_src/external/drifters;src/MOM6/pkg/GSW-Fortran/modules/;src/MOM6/pkg/GSW-Fortran/toolbox/;src/MOM6/src/ALE/;src/MOM6/src/core/;src/MOM6/src/diagnostics/;src/MOM6/src/equation_of_state/;src/MOM6/src/framework/;src/MOM6/src/ice_shelf/;src/MOM6/src/initialization/;src/MOM6/src/ocean_data_assim/;src/MOM6/src/parameterizations/;src/MOM6/src/tracer/;src/MOM6/src/user/;src/MOM6/src/equation_of_state/TEOS10/;src/MOM6/src/framework/testing/;src/MOM6/src/parameterizations/CVmix/;src/MOM6/src/parameterizations/lateral/;src/MOM6/src/parameterizations/stochastic/;src/MOM6/src/parameterizations/vertical/;src/FMS/coupler;src/FMS/include;src/ocean_BGC/generic_tracers;src/ocean_BGC/mocsy/src \
    -Docean_model_macros=MAX_FIELDS_=600;NOT_SET_AFFINITY;_USE_MOM6_DIAG;_USE_GENERIC_TRACER;USE_PRECISION=2 \
    -Docean_model_flags=-Isrc/FMS;-Isrc/MOM6;-Isrc/MOM6/src/framework/;-Isrc/ocean_BGC/mocsy/src \
    -Datmos_model_paths=src/GFDL_atmos_cubed_sphere/driver/GFDL;src/GFDL_atmos_cubed_sphere/model;src/GFDL_atmos_cubed_sphere/GFDL_tools;src/GFDL_atmos_cubed_sphere/tools;src/atmos_drivers/coupled \
    -Datmos_model_macros=CLIMATE_NUDGE;SPMD \
    -Datmos_model_flags=-Isrc/FMS;-Isrc/FMS/include;-Isrc/GFDL_atmos_cubed_sphere/tools/ \
    -Dland_model_paths=src/lm4p \
    -Dland_model_flags=-Isrc/FMS;-Isrc/FMS/include \
    -Dice_model_paths=src/SIS2/config_src/dynamic;src/SIS2/config_src/external/Icepack_interfaces;src/SIS2/src;src/icebergs/src;src/ice_param \
    -Dice_model_macros=USE_FMS2_IO \
    -Dice_model_flags=-Isrc/FMS;-Isrc/MOM6;-Isrc/MOM6/src/framework/ \
    -Datmos_physics_library_paths=src/atmos_phys \
    -Datmos_physics_library_flags=-Isrc/FMS;-Isrc/FMS/include;-Isrc/atmos_phys/atmos_param/cosp
```

#### Differences with mkmf

The cmake build is intentionally similar to a mkmf compilation but still has some differences in functionality.

For mkmf, you can use any number of libraries, but for the cmake build you are currently restricted to one library per "model"/dependency.

Targets and templates are the other big change. CMake will automagically detect compiler vendors, so templates do not need to be specified.
Instead, the target is determined by the `-DCMAKE_BUILD_TYPE` argument, which follows the same `Prod`,`Repro`, `Debug` naming scheme.
Additional target options are provided as CMake options: `-DOPENMP`, `-DOVERRIDE_LIMITS`, `-DVERBOSE_FLAGS`
Compiler options based off the original mkmf templates are specified in the `cmake` subdirectory.

The extra preprocessing step that mkmf can do via the `--use-cpp` option is available via the `-DPREPROCESS_LAND` option. It only applies to the land model since
that's where its typically used, but this can be changed if needed. This is required because fortran preprocessors don't support the "stringify" `#` operator
in macro functions, so the direct `cpp` command needs to be used prior to compilation.


## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an 'as is' basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

This project code is made available through GitHub but is managed by NOAA-GFDL
at https://gitlab.gfdl.noaa.gov.
