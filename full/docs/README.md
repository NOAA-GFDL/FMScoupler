# FMSCoupler — Full Coupled Model

## Introduction

`coupler_main` (`coupler_main.F90`) is the top-level driver for the fully coupled Earth-system model that contains
the main time-stepping loops and calls to exchange fluxes between the atmosphere, ocean, land, and sea-ice components. 
Each component can be modeled on an independent grid with a different MPI domain decomposition.

---

## Model Component State Types

The following derived types hold the instantaneous state for each component:

| Variable | Type | Description |
|---|---|---|
| Atm | atmos_data_type | Atmosphere model state at the current timestep |
| Land | land_data_type | Land model (LM4) state at the current timestep |
| Ice | ice_data_type | Sea-ice model (SIS2) state at the current timestep |
| Ocean | ocean_public_type | Public ocean model (MOM6) state; target for the Ocean_state pointer |
| Ocean_state | ocean_state_type (pointer) | Opaque pointer to the full MOM6 interior state |

See [AtmosDataType.md](AtmosDataType.md), 
[LandDataType.md](LandDataType.md), 
[IceDataType.md](IceDataType.md), 
[OceanPublicType.md](OceanPublicType.md), and 
[OceanStateType.md](OceanStateType.md) for detailed field-level documentation.

---

## Boundary Exchange Types

These derived types hold the fields exchanged at each component interface:

| Variable | Type | Description |
|---|---|---|
| Atmos_land_boundary | atmos_land_boundary_type | Fields exchanged between atmosphere and land |
| Atmos_ice_boundary | atmos_ice_boundary_type | Fields exchanged between atmosphere and sea ice |
| Land_Ice_Atmos_Boundary | land_ice_atmos_boundary_type | Aggregated surface state returned from land and ice to the atmosphere |
| Land_ice_boundary | land_ice_boundary_type | Runoff and calving fields passed from land to ice |
| Ice_ocean_boundary | ice_ocean_boundary_type | All fluxes passed from ice to the ocean |
| Ocean_ice_boundary | ocean_ice_boundary_type | Ocean surface state (SST, currents, frazil) passed up to the ice model |
| ice_ocean_driver_CS | ice_ocean_driver_type (pointer) | Control parameters for the combined ice–ocean driver |

See [FLUX.md](FLUX.md) for details on flux exchange.

---

## Time Integration

There are two nested time-integration loops in coupler_main:

- **Slow (coupled) loop** (`nc = 1 … num_cpld_calls`): 
advances the ocean by one coupled timestep `dt_cpld`. 
- **Fast (atmospheric) loop** (`na = 1 … num_atmos_calls`): 
advances the atmosphere, land surface, and fast sea ice by one atmospheric timestep `dt_atmos`. 
The fast loop is nested inside the slow loop.

### Fast Loop

Heat and moisture are exchanged between the atmosphere and the surface (land/ice) 
using a tridiagonal scheme for implicit vertical diffusion.  The main calls are:

1. `coupler_update_atmos_model_down` — forward (downward) sweep from the atmospheric top to the surface.
2. `coupler_update_atmos_model_radiation` - call the radiation driver (see below for do_concurrent_radiation)
3. `coupler_flux_down_from_atmos` — transfer forward-elimination coefficients to the land and ice models.
4. `coupler_update_land_model_fast` and `coupler_update_ice_model_fast` — in addition to updating hydrology
      and canopy physics, update surface temperature.
5. `coupler_flux_up_to_atmos` — apply implicit surface-flux corrections using the updated surface temperatures.
6. `coupler_update_atmos_model_up` — back-substitution (upward sweep), update convection and large-scale condensation.

When `do_concurrent_radiation = .true.`, steps 1 with steps 3-6 will run simultaneously with 
`coupler_update_atmos_model_radiation` using OpenMP threading. 
(see the `atmos_nthreads` and `radiation_nthreads` namelist variables).

### Slow Loop

When `concurrent = .true.`, the ocean model is updated concurrently with the atmosphere and is one coupling step behind.
With MOM6 model, it is recommended to set `use_lag_fluxes = .false.` so that `flux_ice_to_ocean` 
is called after `update_ocean_model`.  For older ocean models, it is recommended to set `use_lag_fluxes = .true.`
to call `flux_ice_to_ocean` before the ocean model update.

** Miscellaneous 
Sea-ice physics can be split into two timescales that can run on different MPI PE sets:

- **Fast ice** (`Ice%fast_ice_pe`): 
thermodynamics and surface-flux coupling at the atmospheric timestep. Fast ice always runs on a subset of the atmosphere PEs.
- **Slow ice** (`Ice%slow_ice_pe`): 
dynamics, freezing/melting, and transport at the coupled (ocean) timestep. The placement of the 
slow-ice PEs depends on the `slow_ice_with_ocean` namelist variable:
  - `slow_ice_with_ocean = .false.` (default): slow and fast ice share the same PEs.
  - `slow_ice_with_ocean = .true.`: slow ice runs on the ocean PEs. In this case, `Ice%pelist` is the union 
                                    of the fast (atmosphere) and slow (ocean) PE sets.

The flag `concurrent_ice = .true.` runs fast- and slow-ice processes concurrently and requires `slow_ice_with_ocean = .true.`.
The flag `combine_ice_and_ocean = .true.` advances the slow ice and ocean together on the ocean PEs. 
Both `concurrent_ice` and `slow_ice_with_ocean` must be `.true.` to use `combine_ice_and_ocean`.

These options are not commonly used.

---

## Pseudocode

For the most common setting where `concurrent = .true.`, `use_lag_forces = .false.`,
`concurrent_ice = .false.`, `do_concurrent_radiation = .false.`, `do_atm = .true.`,
`do_land = .true.`, `do_ice = .true.`, and `do_ocean = .true.`:

```fortran
call initialization routines 

do nc = 1, num_cpld_calls

  ! Redistribute Ocean surface states to Ocean_ice_boundary
  if (Ocean%is_ocean_pe) then
    synchronize ice and ocean pes 
    call flux_ocean_to_ice(...)
  end if

  ! Map Ocean_ice_boundary to Ice
  call unpack_ocean_ice_boundary(...)

  ! Map Ice%sCS (slow control structure) to Ice%fCS (fast constrol structure)
  call exchange_slow_to_fast_ice(...)

  ! Prepare ice surface fields (albedo, T, etc.)
  call set_ice_surface_fields(...)

  ! generate suface exchange grid between land, ice, and atm
  call generate_sfc_xgrid(...)

  do na = 1, num_atmos_calls

    ! Copy Atm%tr_bot -> Atm%fields for gas-exchange tracers
    call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)

    ! Compute surface exchange coefficients and turbulent fluxes on the exchange grid
    call sfc_boundary_layer(...)

    ! Atmosphere dynamical core (FV3)
    call update_atmos_model_dynamics(...)

    ! Update radiation
    call update_atmos_model_radiation(...)

    ! Downward sweep of the implicit tridiagonal diffusion
    call update_atmos_model_down(...)

    ! Apply implicit atmosphere diffusion correction; pass updated surface fluxes
    ! (heat, moisture, momentum) to land and ice boundary types
    call flux_down_from_atmos(...)

    ! Fast land physics (hydrology, canopy, soil temperature)
    call update_land_model_fast(...)

    ! Fast ice thermodynamics (surface energy balance, melt ponds)
    call update_ice_model_fast(...)

    ! Recompute surface fluxes using updated land/ice surface temperatures
    call flux_up_to_atmos(...)

    ! Back-substitution (upward) sweep; convection; large-scale condensation
    call update_atmos_model_up(...)

    ! Remap atmosphere gas/tracer fields onto exchange grid;
    ! compute air-sea deposition fluxes; deallocate exchange-grid arrays
    call flux_atmos_to_ocean(...)
    call flux_ex_arrays_dealloc()

    ! Advance atmosphere diagnostics and tracer state
    call update_atmos_model_state(...)

  end do  ! fast loop

  ! Slow land physics (routing, carbon, DGVM)
  call update_land_model_slow(...)

  ! Interpolate land runoff, calving, and heat fluxes onto ice grid
  call flux_land_to_ice(...)

  ! Reset fast-ice accumulators; copy Land_ice_boundary into Ice
  call ice_model_fast_cleanup(...)
  call unpack_land_ice_boundary(...)

  ! Exchange fast-ice averages to slow-ice side
  call exchange_fast_to_slow_ice(...)

  ! Slow ice physics (dynamics, freezing/melting, transport)
  call update_ice_model_slow(...)

  if (Ocean%is_ocean_pe) then
    
    ! Interpolate ice-bottom fluxes onto the ocean grid -> Ice_ocean_boundary
    call flux_ice_to_ocean(...)
    
    ! Advance ocean state by dt_cpld using Ice_ocean_boundary forcing
    call update_ocean_model(...)
  endif

end do

call ending routines
```

---

## Setting the Model Start Time

The model start time (`Time_start`) is determined in `coupler_init` by the following priority:

1. If `date_init` exists in the `diag_table`, it is used to set the model start time.
2. If `date_init` is not found in the `diag_table` and `INPUT/coupler.res` exists, 
   the start date and calendar type are read from that file. These values can be overridden if 
   `force_date_from_namelist = .true.` and `current_date` with `calendar_type` are defined in `coupler_nml`.
3. If neither source is available, the start date is taken from `current_date` and `calendar` in `coupler_nml`.

---

## MPI Parallelization

The number of MPI processing elements (PEs) for each model component is specified in `coupler_nml`:

```fortran
&coupler_nml
  atmos_npes = 10
  ocean_npes = 10
  ice_npes   = 10
  land_npes  = 10
/
```

Constraints:
- At least one of `atmos_npes` or `ocean_npes` must be specified.
- `land_npes <= atmos_npes`
- `ice_npes <= atmos_npes`
- `atmos_npes + ocean_npes = npes` (total PE count determined by FMS)

When `concurrent = .true.`, 
- The atmosphere and ocean have distinct PE lists.
- `land%pelist` is a subset of `atm%pelist`.
- `ice%pelist = ice%slow_pelist = ice%fast_pelist`, which are a subset of `atm%pelist`.

---

## OpenMP Parallelization

OpenMP thread counts are also configured in `coupler_nml`:

```fortran
&coupler_nml
  do_concurrent_radiation = .true.
  use_hyper_thread        = .true.
  conc_nthreads           = 2
  atmos_nthreads          = 1
  radiation_nthreads      = 1
  ocean_nthreads          = 1
/
```

The model must be compiled with OpenMP for OpenMP threading.
When using the FRE build system, the model will compile with OpenMP
if the target contains the `-openmp` suffix (e.g., `prod-openmp`).

When `do_concurrent_radiation = .true.`, `conc_nthreads` is set to 2: 
thread 0 on the atmosphere PEs runs atmosphere dynamics and physics, while thread 1 runs the radiation code.
When `do_concurrent_radiation = .false.`, radiation runs sequentially after the atmosphere update.
The `atmos_nthreads` variable controls the number of threads used within the atmosphere dynamics itself.

---

## Time Variables

| Variable | Type | Description |
|---|---|---|
| Time_step_atmos | FmsTime_type | timestep for the fast-loop (atmospheric) |
| Time_step_cpld | FmsTime_type | timestep for  the one slow-loop (coupled) |
| Time_atmos | FmsTime_type | Current model time tracked in the fast loop |
| Time_ocean | FmsTime_type | Current model time tracked for the ocean |
| Time_flux_ice_to_ocean | FmsTime_type | Time when flux was exchanged from ice to ocean |
| Time_flux_ocean_to_ice | FmsTime_type | Time when flux was exchange ocean to ice |
| Time_restart | FmsTime_type | Next scheduled time for writing intermediate restarts |
| Time_restart_current | FmsTime_type | Time at which the most recent intermediate restart was written |
| Time_start | FmsTime_type | Model start time |
| Time_end | FmsTime_type | Model end time |

---

## Loop Counters

| Variable | Type | Description |
|---|---|---|
| num_atmos_calls | integer | Number of iterations in the fast loop |
| na | integer | Do-loop counter for the fast loop |
| num_cpld_calls | integer | Number of coupled timesteps in the run |
| nc | integer | Do-loop counter for the slow loop |
| current_timestep | integer | Total number of fast-loop iterations during runtime, 
  equal to `(nc-1)*num_atmos_calls + na` |
