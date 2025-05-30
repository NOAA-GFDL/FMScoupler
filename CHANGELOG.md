# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [2025.02] - 2025-04-28
## Added
- FULL: Added calls to the added generic_exchange(gex) module which was added to FMS for this release. These modifications handle dry deposition on the land side, fix errors in ocean flux diagnositics with vmr tracers for ESM4p5, and pass fluxes between the land and atmosphehere. (#151)

## Fixed
- SIMPLE: Added initialization call for generic_exchange module (#155)
- FULL: Add allocations to fix crashes in models that use \_USE\_LEGACY\_LAND\_ (#156)

### Tag Commit Hashes
- 2025.02-beta2  88994b1f49d5ffe40c765ee9dc69e425df3e1f3f
- 2025.02-beta1  57c2538100c20dafa421c6683af8a786eca0c45d
- 2025.02-alpha1 16ec1e4bfc070520cd6eeb44a6712062c122d34d

## [2025.01] - 2025-01-30

## Added
- FULL: Added ocean temperature and rough heat variables to `land_ice_atmos_boundary` to be used and populated on the exchange grid from ocean variables in `sfc_boundary_layer`, before running the atmos physics for shield coupling. (#143)

### Tag Commit Hashes
- 2025.01-beta1  008399dcf1a7a2d8b4733cded744d9a72ca14f10
- 2025.01-alpha1 008399dcf1a7a2d8b4733cded744d9a72ca14f10

## [2024.03.01] - 2024-09-06

## Fixed
- FULL: Adds conditionals to bypass `coupler_generate_sfc_xgrid`, `flux_exchange_init`, and flux exchange related stock routines.

### Tag Commit Hashes
- 2024.03.01 e357d8ffc44c2808628ec7912c827d2ad41978b7

## [2024.03] - 2024-08-22

### Fixed
- FULL: Removes unused FREDB related code that was causing errors with SPEAR (#137)
- SHIELD: Fixed compile error from missing `&` in SHiELD coupler_main (#140)
- FULL: Fixes crashes when the `do_flux` namelist flag is false by adding conditionals for calls to `coupler_flux_down_from_atmos`, `coupler_flux_up_to_atmos` and `coupler_flux_atmos_to_ocean` (#134)
- SHIELD: Fixed a bug where `data_override` was not working correctly due to missing an argument during initialization. (#131)

### Tag Commit Hashes
- 2024.03-beta1 8dc6eff18b67ca9ab83916e2326af797818bf463

## [2024.02] - 2024-07-11

### Known Issues
- SPEAR errors have been reported due to a string length being exceeded during the `fredb` section of the code in the full coupler before initialization. (#136)

### Added
- SIMPLE: Adds option to obtain concurrent drag and heat transfer coefficents via `monin_obukhov` functions. To control the new option, new namelist options have been added to the `surface_flux_nml`. The `rough_scheme_ocean` in `ocean_rough_nml` must be set to `hwrf17` for the options to be effective. (#113)
	- `do_iter_monin_obukhov` controls whether to iterate over ocean with updated roughness length
	- `niter_monin_obukhov` controls the number of times to iterate
	- `use_u10_neutral` if true uses 10m neutral wind rather than the standard 10m wind to obtain `rough_mom`, `rough_heat`, and `rough_moist`

### Changed
- FULL: A number of changes have been made as part of a refactorization effort for the full coupler.
	- A module `full_coupler_mod` has been added and contains routines and imports previously included in the `coupler_main.F90` file. (#104)
	- A new type, `coupler_clock type`, has been added to handle the usage of any clocks used in the full coupler. This will cause changes in the output order of clock statistics (#106)
	- Updates the written version number to use the version in `file_version.h` (#109)
	- Refactors `atmos_ice_land_chksum` and `ocean_chksum` calls into new routines in `full_coupler_mod` (#114)
	- Refactors `flux_ocean_to_ice` and `flux_ice_to_ocean` calls into new routines in `full_coupler_mod` (#117)
	- Refactors exchange updates for ice as well as setting ice surface fields and unpacking ocean ice boundary into new routines (#118)
	- Refactors ice surface exchange grid generation into a new routine (#119)
	- Adds a `coupler_chksum_type` used to compute checksums when `do_chksum` is true and adds wrapper routines for `atmos_tracer_driver_gather_data` and `sfc_boundary_layer` (#123)
	- Refactors `update_atmos_model_dynamics`, `update_atmos_model_radiation`, `update_atmos_model_down`, and `flux_down_from_atmos` into new routines (#124)
	- Refactors `update_land_model_fast`, `update_ice_model_fast`, `update_atmos_model_up`, `flux_up_to_atmos`, and `update_atmos_model_state` into new routines (#125)
	- Refactors `unpack_land_ice_boundary`, `update_ice_model_slow_and_stocks`, and `update_ocean_model` into new routines (#127)
	- Refactors routines for writing intermediate restarts and summarizing timesteps (#128)
- SHIELD: `data_override` has been enabled for nest domains (#102)

### Removed
- FULL: Two variables that enabled bugged behaviour from past releases have been removed, `ex_u_star_smooth_bug` and `sw1way_bug`. (#99)

### Tag Commit Hashes
- 2024.02-alpha1 eda97dd590c7c83f3664fbbacf22a76f3234251d
- 2024.02-alpha2 2b59f6d4e52e1c7c50f6fba09918f4a0b4ad5741
- 2024.02-beta1  2a51893473d84a106e95da33e0f05225ca49a3db
- 2024.02-beta2  0b2c91402236e49d586dd5151c98d329e97f9e57

## [2024.01] - 2024-05-03

### Added
- Adds option to run SHiELD with a gregorian calender (#105)
- Adds a `diag_send_complete` call for FMS diag_manager rewrite(#101)
- Adds documentation for the exchange grid (#98)

### Tag Commit Hashes
- 2024.01-alpha1 6442d387153064644325c96a5e9e2935139d5e3c
- 2024.01-alpha2 6442d387153064644325c96a5e9e2935139d5e3c
- 2024.01-alpha3 f4782c2c033df086eeac29fbbefb4a0bdac1649f
- 2024.01-beta1  f4782c2c033df086eeac29fbbefb4a0bdac1649f
- 2024.01-beta2  f4782c2c033df086eeac29fbbefb4a0bdac1649f
- 2024.01-beta3  d15c35a92ac4f04c57c539eaa56301c8f70d53cf
- 2024.01-alpha4 d15c35a92ac4f04c57c539eaa56301c8f70d53cf
- 2024.01-alpha5 d15c35a92ac4f04c57c539eaa56301c8f70d53cf
- 2024.01-alpha6 4dc9b0f2a85d34b0fdc8477625b91794a77ac747
- 2024.01-beta4  4dc9b0f2a85d34b0fdc8477625b91794a77ac747
- 2024.01-beta5  4dc9b0f2a85d34b0fdc8477625b91794a77ac747

## [2023.04.01] - 2024-06-20
### Fixed
- SIMPLE COUPLER: Fixed issue where the `sst_anom` value used was not the correct value read in from the namelist. (#122)

### Tag Commit Hashes
- 2023.04 6b4f8b5207483eb7a7839a29909a415253a05db8

## [2023.04] - 2023-12-04
### Added
- Adds additional output arguments `thv_atm` amd `thv_surf` to the `surface_flux` interface, as well as calls to xgrid and send data in order to use a new atmosphere boundary layer scheme.

### Changed
- Routines using the `data` argument name explicictly have been updated to match corresponding FMS updates

### Tag Commit Hashes
2023.04-beta1 93ce3642a7951eb11d7d39441911717923dfc768


## [2023.02] - 2023-07-27
### Fixed
- SHARED: Fixes crashes due to uninitialized namelist variables.

### Changed
- Routines/variables from FMS have been updated to include a prefixes containing fms and the subdirectory/module name. This was necessitated by aliases being added to the 'global' libFMS.F90 module in FMS.
- FULL: Adds logic for PE assignment in order to allow for a data atmosphere to be used alongside the combined ice-ocean driver.

### Removed
- Usage of fms_io and mpp_io has been deprecated. If using these modules in a model, you must compile both FMS and FMScoupler with the -Duse_deprecated_io CPP flag.

### Tag Commit Hashes
2023.02-alpha1 9da2d61f74671ccee775553a36439260c9241383
2023.02-alpha2 9da2d61f74671ccee775553a36439260c9241383
2023.02-alpha3 4ca21a7f3dc3649f934ad7b34e1a61b63e589712
2023.02-beta1  78c438457cd49a82f6eaec2d57638ffc5084c688


## [2023.01] - 2023-04-03
### Fixed
- Fixed IO domain related failures coming from ice model for the null model test
### Added
- Added clock optimizations to the SHiELD coupler

### Tag Commit Hashes
2023.01-alpha1 7c47be33b4049a96bbce3d9b4cc165dbb147e751
2023.01-alpha2 7c47be33b4049a96bbce3d9b4cc165dbb147e751
2023.01-alpha3 7c47be33b4049a96bbce3d9b4cc165dbb147e751
2023.01-alpha4 7c47be33b4049a96bbce3d9b4cc165dbb147e751
2023.01-beta1  2571fc016866898255559355b92347cd354082ce
2023.01-beta2  2571fc016866898255559355b92347cd354082ce
2023.01-beta3  2571fc016866898255559355b92347cd354082ce
2023.01-beta4  2571fc016866898255559355b92347cd354082ce
2023.01-beta5  2571fc016866898255559355b92347cd354082ce

## [2022.03] - 2022-08-01
### Added
- Added doxygen comments for the simple and shield couplers, and general layout improvements for the generated site. It is now updated upon releases and hosted at noaa-gfdl.github.io/FMScoupler

### Tag Commit Hashes
2022.03-alpha1 c12af876d4baef20346db8391422b6b6df209c75
2022.03-beta1  6d4d2b3dce8152400c8c15551763835689de8ddb

## [2022.02] - 2022-04-29
### Removed
- Removes grid code and variables from SHiELD/coupler_main and fixes data_override_init
- Removes outdated logic in simple coupler for data_override_init parsing
### Changed
- Changes routine names used for constants in order to compile with recent constants changes to FMS
### Fixed
- FULL: Replaced a deprecated OpenMP routine causing warnings
- SIMPLE: Fixed a missing variable allocation that was causing failures with certain compilers

### Tag Commit Hashes
2022.02-alpha1 de3e3cbca349021a545a500f5ba1af6af22acfae
2022.02-alpha2 c23b6f3ff1f902adf1fa43f8a5c9d2307bd01106
2022.02-beta1  2bb8f35e2f579e738b58c610c35ca9afd7e36358

## [2022.01] - 2022-03-25
### Added
- Added SHiELD main driver program to the repository
- Added some additional information on the coupler to the readme

### Tag Commit Hashes
- 2022.01-alpha1 4707b255c842dd08b3cd65a45b7924e7a9d88720
- 2022.01-beta1  4707b255c842dd08b3cd65a45b7924e7a9d88720
- 2022.01-alpha2 7790c5d8e243eb97f3fa87f15546dadc6d963ed1
- 2022.01-beta2  c9f405a2383451550b9a8aaa8279a2a973d65c90

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
- FMS_io was almost completely removed from FMScoupler and replaced with fms2_io
	- The only remaining usage is `fms_io_exit` calls. These calls are made at the end of each coupler driver and only called if the CPP macro `use_deprecated_io` is set.
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
