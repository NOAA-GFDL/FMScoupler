# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [2021.03] - 2021-08-16
### Fixed
- In the full coupler, corrects a `get_variable_size` call to prevent crashes when running with the SCM

## [2021.02] - 2021-05-20
### Added
- FMS2_IO was implemented to the full coupler:
	- The coupler restart files are now read with fms2_io's ascii_read
	- Ascii writes are now done with fortran's open, close, and write. They are wrapped in an if, so that only the root pe does the io, newunit ensures that the unit number is unique for each file. 
	- The variables named `unit` have been renamed to avoid fortran conflicts. 
	- The coupler type restarts are now written with fms2_io.
	- The grid file is now read with fms2_io in: full/flux_exchange.F90:check_atm_grid
- FMS2_IO was implemented to the simple coupler:
	- Changed ice_model.F90 so that it reads the land_mask files and it reads/writes the ice restart with fms2io.
	- Removed the native formatted restart file code
	- Fms2_io ascii_read is used to read to the coupler_restart
	- Fotran's `open`, `close`, and `write` are used to write the coupler_restart
	- Removed the read_grid_data and get_grid_size subroutines from simple/ice_model.F90. These are never used. 
- Test cases added for varying the latitude of SST maximum in the simple coupler ice model.
### Changed
- Changes all imports from FMS to use the global `FMS` module and the `FMSconstants` module
- Changes to the ice_model restart files:
  - With fms_io, the variables `mlon` and `mlat` are written as: <br>
    		`double mlon(Time, zaxis_1, yaxis_1, xaxis_1) ;`<br>
  		`double mlat(Time, zaxis_1, yaxis_1, xaxis_1) ;`<br>
    where:<br>
	`xaxis_1 = 1 ;`<br>
	`yaxis_1 = 1 ;`<br>
	`zaxis_1 = 1 ;`<br>
	`Time = UNLIMITED ; // (1 currently)`
  - fms_io wrote this as 4d variables as default. In FMS2_io these variables are scalars as they should be. With fms_io, the other non scalar variables where written as:<br>
    	`double t_surf(Time, zaxis_1, yaxis_2, xaxis_2) ;`<br>
    where:<br>
    	`zaxis_1 = 1 ;`
  - In FMS2_io these variables are not a function of the zaxis_1. They are 2d + time as they should be.<br>
    	`double t_surf(Time, yaxis_1, xaxis_1) ;`<br>
    With fms_io, the variables attributes:<br>
    	`long_name = {The same name as the longname}`<br>
	`units = "none"`<br>
    were written by default.
  - FMS2_io does not do this. Users can specify real long_names and units by calling register_variable_attribute.
### Removed
- FMS_io was almost completely removed from FMScoupler and replaced with fms2_io. 
### Tag Commit Hashes
- 2021.02-alpha1 (c1c8044a6c3efb8ddbbd01a3769bbf2610b34937)
- 2021.02-alpha2 (c1c8044a6c3efb8ddbbd01a3769bbf2610b34937)
- 2021.02-beta1 (62415ea3a62145080efcfe078eb889d6adf681a1)
- 2021.02-beta2 (62415ea3a62145080efcfe078eb889d6adf681a1)

## [2021.01] - 2021-03-08

### Added
- SURFACE_FLUX: Adds a new functionality to enable using NCAR surface fluxes in experiments

### Fixed
- SIMPLE_COUPLER: Fixed issue with simpler coupler not calling data_override_init during initialization, will now call if the data_table file exists 

## Tag Commit Hashes
- 2021.01-beta1 (7e7212c6db62aa7916af0f6ada59c5a83355c1b8)
- 2021.01-alpha2 (4c2de8d2210f77cff38b7ecd8d2c06e7333a0d9e)
- 2021.01-alpha1 (36ed0f9cbe6d158f5205c418394f6539e0435237)

## [2020.04] - 2020-12-07

## Tag Commit Hashes
- 2020.04-beta1 (c8d7687653659a06de23c95785038c21857e5792)
- 2020.04-alpha3 (825b948e390fd68fa3bc46200c83623c5b66a293)
- 2020.04-alpha1 (3b4ed9f9327f341d3a690264acd46416293dbe1f)

## [2020.03] - 2020-10-08

### Added
- COUPLER_MAIN: Adds support for gregorian calenders

### Removed
- ATM_LAND_ICE_FLUX_EXCHANGE: Removed deprecated calls for creating and setting up Ice%ocean_fluxes_top

### Tag Commit Hashes
- 2020.03-beta4 (0b99d9ec5601cf907767806e265b252a1720b301)
- 2020.03-beta3 (0b99d9ec5601cf907767806e265b252a1720b301)
- 2020.03-beta2 (0b99d9ec5601cf907767806e265b252a1720b301)

## [2020.02] - 2020-05-01

### Fixed
- ATMOS_OCEAN_FLUXES_CALC: Fixes div_by_zero error in debug when the ocean gas concentration is zero (e.g. over land points). An epsilon value (1e30) is substituted at points that have zero gas concentration.

### Tag Commit Hashes
- 2020.02-beta1 (87e5798ddbb82a5011dfaa0dc0eb3c9231de18b1)
- 2020.02-beta2 (6aa98ccbeda8b254b5ba1ccb46d3ae7379ef7a4c)
- 2020.02-beta3 (6aa98ccbeda8b254b5ba1ccb46d3ae7379ef7a4c)

## [2020.01] - 2020-03-13

### Tag Commit Hashes
- 2020.01-alpha1 (f48e67ba9c343152045182e80df5c14021130d47)
- 2020.01-beta1 (5159c6713c4e2600227041e19135ffdcb5032aff)
- 2020.01-beta2 (5159c6713c4e2600227041e19135ffdcb5032aff)

## [2019.01]

### Tag Commit Hashes
- testing_20190307 (2766a232809143778d77bc6918236e7085044b89)
- testing_20190422 (2766a232809143778d77bc6918236e7085044b89)
- testing_20190705 (14a9be493037a07f058adba947d8ce5af58af5d7)
- testing_20190809 (14578f09a25c8e6101faba18342630af267cdba9)
