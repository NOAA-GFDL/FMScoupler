# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

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
