# `ice_ocean_driver_type`

Defined in `sis2/src/combined_ice_ocean_driver.F90:38`.

`ice_ocean_driver_type` is the control structure for the combined ice-ocean driver module (`combined_ice_ocean_driver`). It is used only when `combine_ice_and_ocean = .true.` in `coupler_nml`, which requires both `concurrent_ice = .true.` and `slow_ice_with_ocean = .true.`. In this mode, the coupler makes a single call to `update_slow_ice_and_ocean` rather than separate calls to `update_ice_model_slow` and `update_ocean_model`.

An instance (conventionally named `ice_ocean_driver_CS`) is held as a pointer in the coupler state and passed opaquely to `update_slow_ice_and_ocean`, `ice_ocean_driver_init`, and `ice_ocean_driver_end`. All members are private to the module.

This module provides a platform for tighter integration between SIS2 and MOM6 and is expected to evolve as coupled ice-ocean numerics develop.

---

## Control flags

| Field | Type | Description |
|---|---|---|
| CS_is_initialized | logical | .true. once ice_ocean_driver_init has completed successfully; guards against use before initialization |
| single_MOM_call | logical | If .true. (default), MOM6 dynamics and thermodynamics are advanced together in a single call to update_ocean_model. If .false., separate calls are made for the two phases |
| intersperse_ice_ocn | logical | If .true., ice and ocean thermodynamic and dynamic updates are interspersed rather than sequential. Requires single_MOM_call=.false. The default is .false. |
| use_intersperse_bug | logical | If .true., retains a bug in the intersperse option where the ocean state was not being passed back to the sea ice between interspersed steps; exists for backward compatibility |

---

## Timestepping

| Field | Type | Description |
|---|---|---|
| dt_coupled_dyn | real | Timestep [s] for coupling ice and ocean dynamics when intersperse_ice_ocn=.true.; set to <0 (default: -1) to use the standard coupled timestep dt_cpld |
