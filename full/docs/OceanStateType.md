# `ocean_state_type`

`ocean_state_type` contains information about the state of the ocean.
Its members are private and are not accessible outside `ocean_model_MOM.F90`.  
The coupler holds only a pointer to it (`type(ocean_state_type), pointer :: Ocean_state`).

---

## PE flag

| Field | Type | Description |
|---|---|---|
| is_ocean_PE | logical | .true. if this PE is part of the ocean pelist; initialised to .false. so non-ocean PEs that hold the pointer do nothing |

---

## Time

| Field | Type | Description |
|---|---|---|
| Time | type(time_type) | The ocean model's master clock; advanced each call to update_ocean_model |
| Time_dyn | type(time_type) | The ocean model's dynamics time; equals Time after a complete timestep but may lag during a split dynamics/thermodynamics step |

---

## Timestepping control

| Field | Type | Description |
|---|---|---|
| nstep | integer | Counter for the number of update_ocean calls that advanced the dynamics |
| nstep_thermo | integer | Counter for the number of update_ocean calls that advanced the thermodynamics |
| single_step_call | logical | If .true. (default), dynamics and thermodynamics are advanced together in a single call. If .false., separate calls are made and dt/dt_therm below are used |
| dt | real | Baroclinic dynamics timestep [T ~> s]; only used when single_step_call=.false. |
| dt_therm | real | Thermodynamics timestep [T ~> s]; only used when single_step_call=.false. |
| thermo_spans_coupling | logical | If .true., thermodynamic and tracer timesteps can span multiple coupled timesteps |
| diabatic_first | logical | If .true., diabatic and thermodynamic processes are applied before the dynamics step |

---

## Restart control

| Field | Type | Description |
|---|---|---|
| Restart_control | integer | Bit-field controlling restart file writing: bit 0 (+1) saves generic restart files; bit 1 (+2) saves time-stamped files. A negative value suppresses restart writing at run end |

---

## Physics switches

| Field | Type | Description |
|---|---|---|
| use_ice_shelf | logical | If .true., the MOM6 ice shelf model (MOM_ice_shelf) is enabled |
| use_waves | logical | If .true., surface wave coupling is active |
| icebergs_alter_ocean | logical | If .true., icebergs can modify ocean dynamics and forcing fluxes |
| calve_ice_shelf_bergs | logical | If .true., icebergs are seeded from ice-shelf flux through the ice front |
| offline_tracer_mode | logical | If .false. (default), full prognostic dynamics and thermodynamics are integrated. If .true., only tracer advection/diffusion is integrated using velocity fields read from a previous run |

---

## Physical constants

| Field | Type | Description |
|---|---|---|
| C_p | real | Specific heat capacity of seawater [J degC⁻¹ kg⁻¹] |
| press_to_z | real | Conversion factor from pressure to ocean depth, typically 1/(ρ₀g) [Z T² R⁻¹ L⁻² ~> m Pa⁻¹] |

---

## Forcing and surface state structures

These nested types carry the mechanical and thermodynamic forcing fields
assembled by the coupler and consumed by `update_ocean_model`.

| Field | Type | Description |
|---|---|---|
| forces | type(mech_forcing) | Mechanical surface forcing: wind stress, sea-level pressure, and wave-related fields |
| fluxes | type(forcing) | Primary thermodynamic ocean forcing: heat flux, freshwater flux, salt flux, shortwave penetration, etc. |
| flux_tmp | type(forcing) | Secondary forcing structure used when multiple coupled timesteps are taken per thermodynamic step; accumulates forcing between thermodynamic updates |
| sfc_state | type(surface) | Ocean surface state fields (SST, SSS, surface currents, boundary layer depth) returned to the coupler after each update |

---

## Grid and vertical coordinate structures

| Field | Type | Description |
|---|---|---|
| grid | type(ocean_grid_type), pointer | MOM6 horizontal grid structure: cell areas, distances, coordinates, metric terms, and land masks |
| GV | type(verticalGrid_type), pointer | MOM6 vertical grid structure: layer thicknesses, target densities (for isopycnal coordinates), and ALE remapping parameters |
| US | type(unit_scale_type), pointer | MOM6 dimensional unit-scaling factors used to convert between external MKS units and MOM6's internal non-dimensionalised units |

---

## Model control structures

| Field | Type | Description |
|---|---|---|
| MOM_CSp | type(MOM_control_struct) | MOM6 master control structure; holds all module-level control structures, parameter settings, and diagnostic handles for the ocean model |
| Ice_shelf_CSp | type(ice_shelf_CS), pointer | Control structure for the MOM6 ice shelf model; null if use_ice_shelf=.false. |
| marine_ice_CSp | type(marine_ice_CS), pointer | Control structure for the marine ice effects module (e.g. iceberg melt parameterisation) |
| Waves | type(wave_parameters_cs), pointer | Control structure for surface wave coupling; null if use_waves=.false. |
| forcing_CSp | type(surface_forcing_CS), pointer | MOM6 surface forcing control structure; handles the translation from coupler boundary conditions to internal MOM6 forcing arrays |
| diag | type(diag_ctrl), pointer | MOM6 diagnostic control structure; manages registration and posting of all ocean diagnostics to the FMS diag_manager |

---

## I/O paths

| Field | Type | Description |
|---|---|---|
| dirs | type(directories) | Structure containing relevant directory paths for input, output, and restart files |
