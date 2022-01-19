![Build Status](https://github.com/NOAA-GFDL/FMScoupler/actions/workflows/build.yml/badge.svg?branch=main)

# Flexible Modeling System (FMS) Coupler

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
