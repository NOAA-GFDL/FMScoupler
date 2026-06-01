# `ice_data_type`

`ice_data_type` is the publicly visible face of the SIS2 sea-ice model. 
It carries the surface state, inter-component fluxes, and PE/domain metadata that the coupler needs.
 All internal SIS2 state is accessed through the private control structures `fCS` (fast ice) and `sCS` (slow ice), 
 which are opaque to the coupler. 

The type supports a split fast-ice / slow-ice architecture. Fast processes (surface thermodynamics, atmosphere–ice flux coupling) run on the atmospheric timestep on atmosphere PEs; slow processes (ice dynamics, freezing/melting, transport) run on the coupled timestep and may run on ocean PEs when `slow_ice_with_ocean = .true.`. See [README-new.md](README-new.md) for a description of how this affects the PE layout.

---

## PE and domain bookkeeping

| Field | Type | Description |
|---|---|---|
| pe | logical | .true. on any PE that participates in ice computation (fast or slow) |
| fast_ice_pe | logical | .true. on PEs in the fast-ice pelist; these PEs handle atmosphere-timestep ice processes |
| slow_ice_pe | logical | .true. on PEs in the slow-ice pelist; these PEs handle coupled-timestep ice dynamics and ice-ocean exchange |
| shared_slow_fast_PEs | logical | .true. when fast and slow ice use the same PE set and domain decomposition (i.e. slow_ice_with_ocean=.false.); .false. when slow ice runs on the ocean PEs |
| pelist | integer(:) | Combined list of all ice PEs (union of fast and slow pelists); used for flux exchange |
| fast_pelist | integer(:) | MPI PE numbers for fast-ice processes |
| slow_pelist | integer(:) | MPI PE numbers for slow-ice processes |
| Domain | type(domain2D) | Copy of the fast-ice FMS domain without halos; used for exchange-grid setup |
| slow_Domain_NH | type(domain2D) | Copy of the slow-ice FMS domain without halos; used for ice-ocean flux redistribution |
| fast_domain | type(domain2D) pointer | Pointer to the fast-ice MPI domain (or an allocated copy on slow-ice PEs) |
| slow_domain | type(domain2D) pointer | Pointer to the slow-ice MPI domain (or an allocated copy on fast-ice PEs) |
| ocean_pt | logical(:,:) | Mask array; .true. at ocean (non-land) points |
| xtype | integer | Transfer mode for ice-ocean flux exchange: DIRECT (3) when ice and ocean share the same decomposition, REDIST (2) when they differ |
| axes | integer(3) | Diag-manager axis IDs for the ice surface grid |
| Time | type(time_type) | The sea-ice model's current clock time |

---

## Surface state fields (fast-ice grid)
These fields are dimensioned `(:, :, n_categories)` and provide per-ice-thickness-category information to the atmosphere each fast timestep. Category 1 is open ocean; the sum of `part_size` over all categories equals 1.

| Field | Type | Description |
|---|---|---|
| part_size | real(:,:,:) | Fractional coverage of the grid cell by each ice thickness category [dimensionless, 0–1]; category 1 is open ocean |
| t_surf | real(:,:,:) | Surface temperature of the ocean (category 1) or each ice thickness category [K] |
| albedo | real(:,:,:) | Broadband surface albedo averaged across all wavelength and orientation bands within each ice category [dimensionless, 0–1] |
| albedo_vis_dir | real(:,:,:) | Surface albedo for direct visible shortwave radiation in each ice category [dimensionless] |
| albedo_nir_dir | real(:,:,:) | Surface albedo for diffuse visible shortwave radiation in each ice category [dimensionless] |
| albedo_vis_dif | real(:,:,:) | Surface albedo for direct near-infrared shortwave radiation in each ice category [dimensionless] |
| albedo_nir_dif | real(:,:,:) | Surface albedo for diffuse near-infrared shortwave radiation in each ice category [dimensionless] |
| rough_mom | real(:,:,:) | Surface roughness length for momentum at the ocean/ice surface, as provided by ocean_rough_mod [m] |
| rough_heat | real(:,:,:) | Surface roughness length for heat at the ocean/ice surface [m] |
| rough_moist | real(:,:,:) | Surface roughness length for moisture at the ocean/ice surface [m] |
| u_surf | real(:,:,:) | Eastward surface velocity of the ocean (category 1, :,:,1) or sea ice [m/s]; used as a lower boundary condition for wind stress |
| v_surf | real(:,:,:) | Northward surface velocity of the ocean (category 1, :,:,1) or sea ice [m/s] |
| flux_uv_stagger | integer | Staggering of u_surf/v_surf relative to tracer points; valid values are AGRID, BGRID_NE, CGRID_NE, BGRID_SW, CGRID_SW (Arakawa notation). Initialized to -999 so that a global max across all PEs propagates the correct value to non-ice PEs |

---

## Ocean surface state (from ocean-to-ice exchange)

| Field | Type | Description |
|---|---|---|
| s_surf | real(:,:) | Ocean surface salinity [g salt / kg seawater]; populated by flux_ocean_to_ice |
| SST_C | real(:,:) | Ocean surface temperature [°C]; used in forecast mode |

---

## Grid and area fields

| Field | Type | Description |
|---|---|---|
| area | real(:,:) | Area of each ocean cell [m²]; land cells have area = 0 and this field can double as a mask |
| mi | real(:,:) | Total ice + snow mass per unit area [kg/m²]; passed to the ocean for pressure loading and to the wave model |

---

## Fluxes passed from ice to ocean (slow-ice side)

These fields are computed by the slow-ice model. 

| Field | Type | Description |
|---|---|---|
| flux_u | real(:,:) | Flux of x-momentum into the ocean (zonal wind/ice stress) [Pa] |
| flux_v | real(:,:) | Flux of y-momentum into the ocean (meridional wind/ice stress) [Pa] |
| flux_t | real(:,:) | Sensible heat flux out of the ocean [W/m²] |
| flux_q | real(:,:) | Evaporative moisture flux out of the ocean [kg/m²/s] |
| flux_lh | real(:,:) | Latent heat flux out of the ocean [W/m²] |
| flux_lw | real(:,:) | Net longwave flux out of the ocean [W/m²] |
| flux_sw_vis_dir | real(:,:) | Direct visible shortwave heat flux into the ocean [W/m²] |
| flux_sw_vis_dif | real(:,:) | Diffuse visible shortwave heat flux into the ocean [W/m²] |
| flux_sw_nir_dir | real(:,:) | Direct near-infrared shortwave heat flux into the ocean [W/m²] |
| flux_sw_nir_dif | real(:,:) | Diffuse near-infrared shortwave heat flux into the ocean [W/m²] |
| flux_salt | real(:,:) | Salt flux out of the ocean (brine rejection / melting) [kg/m²/s] |
| salt_left_behind | real(:,:) | Salt remaining in the ocean during ice growth (not ejected as brine) [kg/m²/s] |
| lprec | real(:,:) | Liquid precipitation flux into the ocean [kg/m²] |
| fprec | real(:,:) | Frozen precipitation flux into the ocean [kg/m²] |
| runoff | real(:,:) | Liquid runoff from land into the ocean [kg/m²] |
| runoff_hflx | real(:,:) | Heat flux associated with liquid runoff, relative to a reference temperature [W/m²] |
| runoff_carbon | real(:,:) | Carbon content of liquid runoff into the ocean [kg/m²] |
| calving | real(:,:) | Calving of ice or frozen freshwater runoff into the ocean [kg/m²] |
| calving_hflx | real(:,:) | Heat flux associated with calving, relative to a reference temperature [W/m²] |
| p_surf | real(:,:) | Pressure at the ocean surface [Pa]; may or may not include atmospheric pressure depending on configuration |
| stress_mag | real(:,:) | Time-mean magnitude of the ice-ocean stress [Pa]; passed to the ocean when pass_stress_mag = .true. in SIS_slow_CS |

---

## Iceberg fields

Present only when the iceberg module is active (`sCS%do_icebergs = .true.`).

| Field | Type | Description |
|---|---|---|
| ustar_berg | real(:,:) | Friction velocity contribution beneath icebergs [m/s]; used for iceberg-ocean drag |
| area_berg | real(:,:) | Fraction of the grid cell covered by icebergs [m²/m²] |
| mass_berg | real(:,:) | Mass of icebergs per unit area [kg/m²] |

---

## Gas and tracer flux boundary conditions

| Field | Type | Description |
|---|---|---|
| ocean_fields | type(coupler_3d_bc_type) | Named surface-state fields shared with the atmosphere (SST, SSS, piston velocities, etc.) for atmosphere-ocean gas exchange; populated by flux_ocean_to_ice |
| ocean_fluxes | type(coupler_2d_bc_type) | Computed gas fluxes from the ice to the ocean for additional tracers |
| ocean_fluxes_top | type(coupler_3d_bc_type) | Gas flux boundary conditions at the top of the ice (atmosphere side); archaic and flagged for eventual removal |

---

## Internal SIS2 control structures

The following fields are private to SIS2 and are not used by other FMS modules. The coupler holds them opaquely and passes them back to SIS2 routines.

| Field | Type | Description |
|---|---|---|
| fCS | type(SIS_fast_CS) pointer | Control structure for the SIS2 fast ice thermodynamics; lives on atmosphere PEs; contains the fast-ice grid (fCS%G), diagnostics, and all fast-timestep state |
| sCS | type(SIS_slow_CS) pointer | Control structure for the SIS2 slow ice dynamics and thermodynamics; may live on ocean PEs when slow_ice_with_ocean=.true.; contains the slow-ice grid (sCS%G), dynamics, tracer registry, and all slow-timestep state |
| icebergs | type(icebergs) pointer | Control structure for the Lagrangian iceberg module; null when do_icebergs=.false. |
| US | type(unit_scale_type) pointer | SIS2 dimensional unit-scaling factors; used to convert between external MKS values and SIS2 internal non-dimensionalized units |
| Ice_restart | type(SIS_restart_CS) pointer | Control structure for writing and reading slow-ice restart files |
| Ice_fast_restart | type(SIS_restart_CS) pointer | Control structure for writing and reading fast-ice restart files |
| OBC | type(ice_OBC_type) pointer | Control structure for ice open boundary conditions; null when OBCs are not configured |
| restart_output_dir | character(240) | Directory path for restart file output; default is './RESTART/' |
