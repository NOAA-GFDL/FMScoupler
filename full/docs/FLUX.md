# Flux Exchange

Authors for initial documentation: Bruce Wyman, V. Balaji, Sergey Malyshev

## Overview

Six modules couple the atmosphere, ocean, land, and sea-ice components through flux exchange:

| Module | Purpose |
|---|---|
| atm_land_ice_flux_exchange | Exchange fluxes between atmosphere, land, and ice |
| atmos_ocean_fluxes_calc | Compute non-deposition gas fluxes between atmosphere and ocean |
| atmos_ocean_dep_fluxes_calc | Compute deposition gas fluxes between atmosphere and ocean |
| ice_ocean_flux_exchange | Exchange fluxes between ice and ocean |
| land_ice_flux_exchange | Exchange fluxes between land and ice |
| flux_exchange | Top-level module that initializes the flux-exchange modules and contains subroutines for stock computation |

## Design Principles

1. The flux exchange supports physically independent atmosphere, land, and sea-ice grids. 
   Ice and ocean must share the same physical grid but their MPI domain decompositions may differ.
2. The masked region of the land grid and the ice-ocean grid must tile each other, i.e., 
   every atmosphere grid cell is covered by either land or ice-ocean, but not both)
3. The masked regions of the ice and ocean grids must be identical.
4. The atmosphere, land, and ice grids exchange information through the surface exchange grid `xmap_sfc` using conservative interpolation (REGRID).
5. The land and ice grids exchange runoff data using the exchange grid `xmap_runoff`.
6. Ice-bottom to ocean transfer does not require an exchange grid because those grids are physically identical. Flux data are moved between PE layouts using `mpp_redistribute` when decompositions differ (REDIST), or copied directly when the decomposition is the same (DIRECT).
7. Atmospheric fluxes reach the ocean through the ice model: first atmosphere to ice via `flux_down_from_atmos`, then ice to ocean via `flux_ice_to_ocean`.
8. Each component model must expose a public data type containing the boundary fields needed by the coupler 
(see [AtmosDataType.md](AtmosDataType.md), [LandDataType.md](LandDataType.md), [IceDataType.md](IceDataType.md), [OceanPublicType.md](OceanPublicType.md)).
9. Sensible heat flux and surface evaporation can depend implicitly on surface temperature.  
   Thus, land and sea-ice temperature updates must run on the atmospheric timestep.
10. Surface fluxes for all other tracers and for momentum are treated as explicit functions of the surface state.
11. The module is designed to support simultaneous implicit time integration on both sides of the surface interface.
12. Because of that implicit coupling, the diffusion part of the land and ice models must also run on the atmospheric timestep.
13. Additional tracer and gas-exchange fluxes are configured through `field_table` and named boundary-condition fields in the coupler boundary types.
14. Any field exchanged between components can be replaced by a constant or file-based value using FMS `data_override` and `data_table`.
15. `update_land_model_fast` and `update_ice_model_fast` must update the surface temperature each atmospheric timestep for the implicit diffusion scheme to be correct.

## Grid Layout

The three component grids tile the sphere. `|xxx|` marks a masked (inactive) grid point:

```
    ATMOSPHERE  |----|----|----|----|----|----|----|----|
          LAND  |---|---|---|---|xxx|xxx|xxx|xxx|xxx|xxx|
           ICE  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
         OCEAN  |xxx|xxx|xxx|xxx|---|---|---|---|---|---|
```

## Transfer Types

There are three transfer modes used internally by the flux exchange:

| Mode | When used |
|---|---|
| REGRID (xtype=1) | Physically distinct grids (atmosphere↔land, atmosphere↔ice); uses the exchange grid with conservative interpolation |
| REDIST (xtype=2) | Same physical grid, different MPI decomposition (ice↔ocean when slow_ice_with_ocean=.true.); uses mpp_redistribute |
| DIRECT (xtype=3) | Same physical grid, same MPI decomposition (ice↔ocean when not using slow_ice_with_ocean); data are copied directly |

---

## Namelist Parameters (`flux_exchange_nml`)

The following parameters can be set in `input.nml` under `&flux_exchange_nml`:

| Parameter | Type | Default | Description |
|---|---|---|---|
| z_ref_heat | real | 2.0 | Reference height [m] for temperature and relative-humidity diagnostics (t_ref, rh_ref, del_h, del_q) |
| z_ref_mom | real | 10.0 | Reference height [m] for momentum diagnostics (u_ref, v_ref, del_m) |
| do_area_weighted_flux | logical | .false. | When .true., fluxes passed to the ocean are multiplied by the ice area fraction before redistribution, so the ocean receives the grid-cell-mean flux rather than the per-ice-area flux |
| debug_stocks | logical | .false. | Enables additional stock-conservation debug output |
| divert_stocks_report | logical | .false. | Redirects stock reporting to stocks.out rather than the standard log |
| do_runoff | logical | .true. | Enables interpolation of land runoff to the ocean; set to .false. to disable |
| do_forecast | logical | .false. | Enables forecast-mode behavior in the flux coupler |
| nblocks | integer | 1 | Number of blocks used to divide n_xgrid_sfc for OpenMP parallelism. In practice this is often set to match coupler_nml%atmos_nthreads; if left at 1 when threading is active, the model will emit a warning and reset it automatically |
| partition_fprec_from_lprec | logical | .false. | For atmosphere-override experiments where liquid and frozen precipitation are combined, converts liquid precipitation to snow when t_ref < tfreeze |
| scale_precip_2d | logical | .false. | Rescales Atm%lprec by a 2-D field read from data_table |

---

## Data Override Capabilities

> **Warning:** The original authors strongly advise against using data override capabilities for "non-experts".

Any field in the lists below can be replaced at runtime with a constant or file-based value by adding a matching entry to `data_table`. A data override is applied only when a matching entry exists; otherwise the model-computed value is used unchanged. See [DATAOVERRIDE.md](DATAOVERRIDE.md) for the complete field listing organized by subroutine.

### `sfc_boundary_layer` — atmosphere boundary to exchange grid

| Field | Description |
|---|---|
| Atm%t_bot | Temperature at the lowest atmospheric level [K] |
| Atm%z_bot | Height of the lowest atmospheric level [m] |
| Atm%p_bot | Pressure at the lowest atmospheric level [Pa] |
| Atm%u_bot | Zonal wind at the lowest atmospheric level [m/s] |
| Atm%v_bot | Meridional wind at the lowest atmospheric level [m/s] |
| Atm%p_surf | Surface pressure [Pa] |
| Atm%slp | Sea-level pressure [Pa] |
| Atm%gust | Gustiness velocity added in quadrature to the resolved wind for surface flux calculations [m/s] |
| atm%fields%bc(n)%field(m)%values | Per-tracer atmospheric surface fields (e.g., tracer concentrations at the lowest model level); iterated over all entries in Atm%fields |

### `sfc_boundary_layer` — ice boundary to exchange grid

| Field | Description |
|---|---|
| Ice%t_surf | Ice/ocean surface skin temperature [K] |
| Ice%rough_mom | Surface roughness length for momentum over ice [m] |
| Ice%rough_heat | Surface roughness length for heat over ice [m] |
| Ice%rough_moist | Surface roughness length for moisture over ice [m] |
| Ice%albedo | Broadband surface albedo over ice [dimensionless] |
| Ice%albedo_vis_dir | Direct-beam visible-band albedo over ice [dimensionless] |
| Ice%albedo_nir_dir | Direct-beam near-infrared albedo over ice [dimensionless] |
| Ice%albedo_vis_dif | Diffuse visible-band albedo over ice [dimensionless] |
| Ice%albedo_nir_dif | Diffuse near-infrared albedo over ice [dimensionless] |
| Ice%u_surf | Zonal surface current/ice velocity [m/s] |
| Ice%v_surf | Meridional surface current/ice velocity [m/s] |
| Ice%ocean_fields | Coupler boundary-condition type holding ocean/ice-top gas and tracer fields used in atmosphere-ocean flux calculations |

### `sfc_boundary_layer` — land boundary to exchange grid

| Field | Description |
|---|---|
| Land%t_surf | Land surface (radiative) temperature [K] |
| Land%t_ca | Canopy air temperature — near-surface air temperature within the plant canopy [K] |
| Land%rough_mom | Surface roughness length for momentum over land [m] |
| Land%rough_heat | Surface roughness length for heat over land [m] |
| Land%albedo | Broadband surface albedo over land [dimensionless] |
| Land%tr | Surface tracer mixing ratio over land; one entry per exchanged tracer |
| Land%albedo_vis_dir | Direct-beam visible-band albedo over land [dimensionless] |
| Land%albedo_nir_dir | Direct-beam near-infrared albedo over land [dimensionless] |
| Land%albedo_vis_dif | Diffuse visible-band albedo over land [dimensionless] |
| Land%albedo_nir_dif | Diffuse near-infrared albedo over land [dimensionless] |

### `sfc_boundary_layer` — exchange grid to `Land_Ice_Atmos_Boundary`

These are the aggregated surface properties returned from land and ice to the atmosphere.

| Field | Description |
|---|---|
| Land_Ice_Atmos_Boundary%t | Surface temperature (area-weighted over land and ice fractions) seen by the atmosphere [K] |
| Land_Ice_Atmos_Boundary%albedo | Broadband surface albedo [dimensionless] |
| Land_Ice_Atmos_Boundary%albedo_vis_dir | Direct-beam visible-band surface albedo [dimensionless] |
| Land_Ice_Atmos_Boundary%albedo_nir_dir | Direct-beam near-infrared surface albedo [dimensionless] |
| Land_Ice_Atmos_Boundary%albedo_vis_dif | Diffuse visible-band surface albedo [dimensionless] |
| Land_Ice_Atmos_Boundary%albedo_nir_dif | Diffuse near-infrared surface albedo [dimensionless] |
| Land_Ice_Atmos_Boundary%land_frac | Fraction of the atmospheric grid cell covered by land [dimensionless] |
| Land_Ice_Atmos_Boundary%dt_t | Implicit temperature correction applied to the atmosphere by the surface flux scheme [K] |
| Land_Ice_Atmos_Boundary%dt_tr | Implicit tracer mixing-ratio correction; one entry per exchanged tracer |
| Land_Ice_Atmos_Boundary%u_flux | Zonal wind stress on the atmosphere [Pa] |
| Land_Ice_Atmos_Boundary%v_flux | Meridional wind stress on the atmosphere [Pa] |
| Land_Ice_Atmos_Boundary%dtaudu | d(wind stress)/d(u) — implicit coupling coefficient for zonal momentum [Pa·s/m] |
| Land_Ice_Atmos_Boundary%dtaudv | d(wind stress)/d(v) — implicit coupling coefficient for meridional momentum [Pa·s/m] |
| Land_Ice_Atmos_Boundary%u_star | Friction velocity (surface turbulent velocity scale) [m/s] |
| Land_Ice_Atmos_Boundary%b_star | Buoyancy scale used in Monin-Obukhov similarity theory [m/s²] |
| Land_Ice_Atmos_Boundary%rough_mom | Roughness length for momentum [m] |

### `flux_down_from_atmos` — atmosphere boundary to exchange grid

| Field | Description |
|---|---|
| Atm%flux_sw | Total net shortwave flux at the surface [W/m²] |
| Atm%flux_sw_dir | Direct-beam component of the net shortwave flux [W/m²] |
| Atm%flux_sw_dif | Diffuse component of the net shortwave flux [W/m²] |
| Atm%flux_sw_down_vis_dir | Downward direct-beam visible shortwave flux [W/m²] |
| Atm%flux_sw_down_vis_dif | Downward diffuse visible shortwave flux [W/m²] |
| Atm%flux_sw_down_total_dir | Downward direct-beam broadband shortwave flux [W/m²] |
| Atm%flux_sw_down_total_dif | Downward diffuse broadband shortwave flux [W/m²] |
| Atm%flux_sw_vis | Net visible-band shortwave flux at the surface [W/m²] |
| Atm%flux_sw_vis_dir | Direct-beam component of the net visible shortwave flux [W/m²] |
| Atm%flux_sw_vis_dif | Diffuse component of the net visible shortwave flux [W/m²] |
| Atm%flux_lw | Net downward longwave flux at the surface [W/m²] |
| Atm%lprec | Liquid precipitation rate [kg/m²/s] |
| frac_precip | 2-D scaling field for liquid precipitation; applied when scale_precip_2d = .true. [dimensionless] |
| Atm%fprec | Frozen (solid) precipitation rate [kg/m²/s] |
| Atm%coszen | Cosine of the solar zenith angle [dimensionless] |
| Atm%Surf_Diff%dtmass | dt/mass — ratio of the timestep to the surface-layer air mass used in the implicit diffusion scheme [s·m²/kg] |
| Atm%Surf_Diff%delta_t | Forward-elimination temperature coefficient from the implicit vertical diffusion scheme [K] |
| Atm%Surf_Diff%dflux_t | d(sensible heat flux)/d(T_surf) — linearization coefficient for the implicit heat flux scheme [W/m²/K] |
| Atm%Surf_Diff%delta_tr | Forward-elimination coefficient for each tracer from the implicit diffusion scheme |
| Atm%Surf_Diff%dflux_tr | d(tracer flux)/d(tracer_surf) — linearization coefficient for the implicit tracer flux scheme |

### `flux_down_from_atmos` — exchange grid to land boundary

| Field | Description |
|---|---|
| Land_boundary%drag_q | Drag coefficient for moisture [dimensionless] |
| Land_boundary%lwdn_flux | Downward longwave radiation flux at the land surface [W/m²] |
| Land_boundary%cd_m | Drag coefficient for momentum over land [dimensionless] |
| Land_boundary%cd_t | Drag coefficient for heat over land [dimensionless] |
| Land_boundary%bstar | Buoyancy scale passed to the land model [m/s²] |
| Land_boundary%ustar | Friction velocity passed to the land model [m/s] |
| Land_boundary%wind | Wind speed at the lowest atmospheric level [m/s] |
| Land_boundary%z_bot | Height of the lowest atmospheric level above the land surface [m] |
| Land_boundary%t_flux | Sensible heat flux into the land surface [W/m²] |
| Land_boundary%lw_flux | Net longwave flux at the land surface [W/m²] |
| Land_boundary%sw_flux | Net shortwave flux at the land surface [W/m²] |
| Land_boundary%sw_flux_down_vis_dir | Downward direct-beam visible shortwave flux over land [W/m²] |
| Land_boundary%sw_flux_down_total_dir | Downward direct-beam broadband shortwave flux over land [W/m²] |
| Land_boundary%sw_flux_down_vis_dif | Downward diffuse visible shortwave flux over land [W/m²] |
| Land_boundary%sw_flux_down_total_dif | Downward diffuse broadband shortwave flux over land [W/m²] |
| Land_boundary%lprec | Liquid precipitation rate over land [kg/m²/s] |
| Land_boundary%fprec | Frozen precipitation rate over land [kg/m²/s] |
| Land_boundary%dhdt | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for land heat flux [W/m²/K] |
| Land_boundary%drdt | d(longwave flux)/d(T_surf) — implicit coupling coefficient for land longwave flux [W/m²/K] |
| Land_boundary%p_surf | Surface pressure over land [Pa] |
| Land_boundary%tr_flux | Flux of each exchanged tracer into the land surface [kg/m²/s] |
| Land_boundary%dfdtr | d(tracer flux)/d(tracer_surf) — implicit coupling coefficient for each land tracer flux |

### `flux_down_from_atmos` — exchange grid to ice boundary

| Field | Description |
|---|---|
| Ice_boundary%u_flux | Zonal wind stress on the ice surface [Pa] |
| Ice_boundary%v_flux | Meridional wind stress on the ice surface [Pa] |
| Ice_boundary%t_flux | Sensible heat flux into the ice surface [W/m²] |
| Ice_boundary%q_flux | Latent heat (moisture) flux into the ice surface [W/m²] |
| Ice_boundary%lw_flux | Net longwave flux at the ice surface [W/m²] |
| Ice_boundary%sw_flux_nir_dir | Direct-beam near-infrared shortwave flux over ice [W/m²] |
| Ice_boundary%sw_flux_vis_dir | Direct-beam visible shortwave flux over ice [W/m²] |
| Ice_boundary%sw_flux_nir_dif | Diffuse near-infrared shortwave flux over ice [W/m²] |
| Ice_boundary%sw_flux_vis_dif | Diffuse visible shortwave flux over ice [W/m²] |
| Ice_boundary%sw_down_vis_dir | Downward direct-beam visible shortwave flux over ice [W/m²] |
| Ice_boundary%sw_down_vis_dif | Downward diffuse visible shortwave flux over ice [W/m²] |
| Ice_boundary%sw_down_nir_dir | Downward direct-beam near-infrared shortwave flux over ice [W/m²] |
| Ice_boundary%sw_down_nir_dif | Downward diffuse near-infrared shortwave flux over ice [W/m²] |
| Ice_boundary%lprec | Liquid precipitation rate over ice [kg/m²/s] |
| Ice_boundary%fprec | Frozen precipitation rate over ice [kg/m²/s] |
| Ice_boundary%dhdt | d(sensible heat flux)/d(T_surf) — implicit coupling coefficient for ice heat flux [W/m²/K] |
| Ice_boundary%dedt | d(latent heat flux)/d(T_surf) — implicit coupling coefficient for ice moisture flux [W/m²/K] |
| Ice_boundary%drdt | d(longwave flux)/d(T_surf) — implicit coupling coefficient for ice longwave flux [W/m²/K] |
| Ice_boundary%coszen | Cosine of the solar zenith angle over ice [dimensionless] |
| Ice_boundary%p | Surface pressure over ice [Pa] |
| Ice_boundary%fluxes | Coupler boundary-condition type holding all per-tracer gas and deposition fluxes from atmosphere to ice |

### `flux_up_to_atmos` — ice boundary to atmosphere

| Field | Description |
|---|---|
| Ice%t_surf | Updated ice surface temperature after the ice model step [K] |

### `flux_up_to_atmos` — land boundary to atmosphere

| Field | Description |
|---|---|
| Land%t_ca | Updated canopy air temperature after the land model step [K] |
| Land%t_surf | Updated land surface temperature after the land model step [K] |
| Land%tr | Updated surface tracer mixing ratio over land; one entry per exchanged tracer |

### `flux_land_to_ice` — land boundary to ice boundary

Land runoff and calving are interpolated from the land grid onto the ice grid using the `xmap_runoff` exchange grid, then passed to the ocean via `flux_ice_to_ocean`. The `do_runoff` namelist flag must be `.true.` for this exchange to occur.

| Field | Description |
|---|---|
| Land_Ice_Boundary%runoff | Liquid runoff (river discharge) from land [kg/m²/s] |
| Land_Ice_Boundary%calving | Solid calving flux (iceberg/glacier discharge) from land [kg/m²/s] |
| Land_Ice_Boundary%runoff_hflx | Sensible heat flux carried by liquid runoff [W/m²] |
| Land_Ice_Boundary%calving_hflx | Sensible heat flux carried by calving discharge [W/m²] |

### `flux_ice_to_ocean` — ice boundary to ocean boundary

All fluxes listed below are passed from the ice model to the ocean. When the ice and ocean PE decompositions differ (`slow_ice_with_ocean = .true.`), these fields are moved using `mpp_redistribute` (REDIST); otherwise they are copied directly (DIRECT).

| Field | Description |
|---|---|
| Ice_Ocean_Boundary%u_flux | Zonal wind/ice stress on the ocean surface [Pa] |
| Ice_Ocean_Boundary%v_flux | Meridional wind/ice stress on the ocean surface [Pa] |
| Ice_Ocean_Boundary%t_flux | Sensible heat flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%q_flux | Latent heat (freshwater) flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%salt_flux | Salt flux from sea ice to the ocean (brine rejection/melting) [kg/m²/s] |
| Ice_Ocean_Boundary%lw_flux | Net longwave flux at the ocean surface [W/m²] |
| Ice_Ocean_Boundary%sw_flux_nir_dir | Direct-beam near-infrared shortwave flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%sw_flux_nir_dif | Diffuse near-infrared shortwave flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%sw_flux_vis_dir | Direct-beam visible shortwave flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%sw_flux_vis_dif | Diffuse visible shortwave flux into the ocean [W/m²] |
| Ice_Ocean_Boundary%lprec | Liquid precipitation reaching the ocean surface [kg/m²/s] |
| Ice_Ocean_Boundary%fprec | Frozen precipitation reaching the ocean surface [kg/m²/s] |
| Ice_Ocean_Boundary%runoff | Liquid runoff from land routed to the ocean [kg/m²/s] |
| Ice_Ocean_Boundary%calving | Solid calving flux routed to the ocean [kg/m²/s] |
| Ice_Ocean_Boundary%runoff_hflx | Sensible heat flux carried by liquid runoff [W/m²] |
| Ice_Ocean_Boundary%calving_hflx | Sensible heat flux carried by calving [W/m²] |
| Ice_Ocean_Boundary%p | Surface atmospheric pressure at the ocean surface [Pa] |
| Ice_Ocean_Boundary%mi | Sea-ice mass per unit area; used for ice-pressure loading on the ocean [kg/m²] |
| Ice_Ocean_Boundary%ustar_berg | Friction velocity beneath icebergs [m/s]; present only when the iceberg module is active |
| Ice_Ocean_Boundary%area_berg | Fractional area of ocean cell covered by icebergs [dimensionless]; present only when the iceberg module is active |
| Ice_Ocean_Boundary%mass_berg | Iceberg mass per unit area [kg/m²]; present only when the iceberg module is active |
| Ice_Ocean_Boundary%fluxes | Coupler boundary-condition type holding all per-tracer gas fluxes from ice/atmosphere to ocean |

### `flux_ocean_to_ice` — ocean boundary to ice boundary

| Field | Description |
|---|---|
| Ocean_Ice_Boundary%u | Zonal ocean surface current velocity [m/s] |
| Ocean_Ice_Boundary%v | Meridional ocean surface current velocity [m/s] |
| Ocean_Ice_Boundary%t | Sea-surface temperature seen by the ice model [K] |
| Ocean_Ice_Boundary%s | Sea-surface salinity seen by the ice model [psu] |
| Ocean_Ice_Boundary%frazil | Frazil ice heat flux from the ocean to the ice [J/m²] |
| Ocean_Ice_Boundary%sea_level | Sea-surface height [m] |
| Ocean_Ice_Boundary%fields | Coupler boundary-condition type holding per-tracer ocean surface fields passed to the ice model |

---

## Diagnostic Fields

All fields below are registered in `atm_land_ice_flux_exchange.F90` inside `diag_field_init`. 
The "Sent in" column identifies the subroutine that calls `send_data` for each field.

### Static fields

| Field | Description | Units | Sent in |
|---|---|---|---|
| land_mask | Fractional amount of land | — | sfc_boundary_layer |
| height2m | Height of the 2 m reference level (scalar axis) | m | sfc_boundary_layer |
| height10m | Height of the 10 m reference level (scalar axis) | m | sfc_boundary_layer |
| sftlf | Fraction of the grid cell occupied by land | 1 | sfc_boundary_layer |

### Atmosphere grid fields

| Field | Description | Units | Sent in |
|---|---|---|---|
| ice_mask | Fractional amount of sea ice | — | diag_sic |
| wind | Wind speed for flux calculations | m/s | sfc_boundary_layer |
| drag_moist | Drag coefficient for moisture | — | sfc_boundary_layer |
| drag_heat | Drag coefficient for heat | — | sfc_boundary_layer |
| drag_mom | Drag coefficient for momentum | — | sfc_boundary_layer |
| rough_moist | Surface roughness length for moisture | m | sfc_boundary_layer |
| rough_heat | Surface roughness length for heat | m | sfc_boundary_layer |
| rough_mom | Surface roughness length for momentum | m | sfc_boundary_layer |
| u_star | Friction velocity | m/s | sfc_boundary_layer |
| b_star | Buoyancy scale | m/s² | sfc_boundary_layer |
| q_star | Moisture scale | kg water/kg air | sfc_boundary_layer |
| thv_atm | Surface air virtual potential temperature | K | sfc_boundary_layer |
| thv_surf | Surface virtual potential temperature | K | sfc_boundary_layer |
| tau_x | Zonal wind stress | Pa | flux_down_from_atmos |
| tau_y | Meridional wind stress | Pa | flux_down_from_atmos |
| t_ocean | Surface temperature from ocean output | K | flux_up_to_atmos |
| t_surf | Surface temperature | K | flux_up_to_atmos |
| t_ca | Canopy air temperature | K | flux_up_to_atmos |
| z_atm | Height of the lowest atmospheric level | m | sfc_boundary_layer |
| p_atm | Pressure at the lowest atmospheric level | Pa | sfc_boundary_layer |
| slp | Sea-level pressure | Pa | sfc_boundary_layer |
| gust | Gustiness scale | m/s | sfc_boundary_layer |
| shflx | Sensible heat flux | W/m² | flux_up_to_atmos |
| lwflx | Net downward longwave flux | W/m² | flux_up_to_atmos |
| t_atm | Temperature at the lowest atmospheric level | K | sfc_boundary_layer |
| u_atm | Zonal wind component at the lowest atmospheric level | m/s | sfc_boundary_layer |
| v_atm | Meridional wind component at the lowest atmospheric level | m/s | sfc_boundary_layer |
| t_ref | Temperature at reference height (z_ref_heat) | K | sfc_boundary_layer |
| rh_ref | Relative humidity at reference height (z_ref_heat) | % | sfc_boundary_layer |
| rh_ref_cmip | Relative humidity at reference height (CMIP name) | % | sfc_boundary_layer |
| u_ref | Zonal wind component at reference height (z_ref_mom) | m/s | sfc_boundary_layer |
| v_ref | Meridional wind component at reference height (z_ref_mom) | m/s | sfc_boundary_layer |
| wind_ref | Wind speed at reference height (z_ref_mom) | m/s | sfc_boundary_layer |
| del_h | Reference-height interpolation factor for heat | — | sfc_boundary_layer |
| del_m | Reference-height interpolation factor for momentum | — | sfc_boundary_layer |
| del_q | Reference-height interpolation factor for moisture | — | sfc_boundary_layer |
| q_ref | Specific humidity at reference height (z_ref_heat) | kg/kg | sfc_boundary_layer |
| rough_scale | Topographic scaling factor for momentum drag | — | sfc_boundary_layer |
| evap | Evaporation rate | kg/m²/s | flux_up_to_atmos |
| co2_bot | CO₂ mixing ratio at the lowest level (from data_override) | ppmv | sfc_boundary_layer |

### Per-tracer atmosphere grid fields

The following fields are registered in a loop over all exchanged tracers (`n_exch_tr`). 
`*` denotes the tracer name from the tracer table.

| Field | Description | Units | Sent in |
|---|---|---|---|
| *_tot_con_atm | Deposition velocity of tracer at the atmosphere level | m/s | flux_up_to_atmos |
| *_tot_con_ref | Deposition velocity of tracer at reference height | m/s | flux_up_to_atmos |
| *_atm | Tracer mixing ratio at the lowest atmospheric level | tracer units | sfc_boundary_layer |
| *_surf | Tracer mixing ratio at the surface | tracer units | flux_up_to_atmos |
| *_flux | Surface flux of tracer | tracer units · kg air/(m²·s) | flux_up_to_atmos |
| *_ref | Tracer mixing ratio at reference height (skipped for sphum) | tracer units | sfc_boundary_layer |
| *_mol_flux | Molar flux of tracer | mol CO₂/(m²·s) or mol/(m²·s) | flux_up_to_atmos |
| *_atm_dvmr | CO₂ mixing ratio at lowest level in dry volume mixing ratio (CO₂ only) | mol CO₂/mol air | sfc_boundary_layer |
| *_surf_dvmr | CO₂ mixing ratio at the surface in dry volume mixing ratio (CO₂ only) | mol CO₂/mol air | flux_up_to_atmos |
| *_mol_flux_atm0 | Gross (one-way) molar flux of tracer | mol/(m²·s) | sfc_boundary_layer |

### CMIP fields

These fields follow CMIP variable naming conventions and are registered via 
`register_cmip_diag_field_2d` (or `fms_diag_register_diag_field` when `use_AM3_physics` is defined).

| Field | Description | Units | Sent in |
|---|---|---|---|
| tas | Near-surface air temperature | K | sfc_boundary_layer |
| uas | Eastward near-surface wind | m/s | sfc_boundary_layer |
| vas | Northward near-surface wind | m/s | sfc_boundary_layer |
| sfcWind | Near-surface wind speed | m/s | sfc_boundary_layer |
| huss | Near-surface specific humidity | 1 | sfc_boundary_layer |
| hurs | Near-surface relative humidity | % | sfc_boundary_layer |
| rhs | Near-surface relative humidity (alternate name) | % | sfc_boundary_layer |
| ts | Surface temperature | K | flux_up_to_atmos |
| psl | Sea-level pressure | Pa | sfc_boundary_layer |
| tauu | Surface downward eastward wind stress | Pa | flux_down_from_atmos |
| tauv | Surface downward northward wind stress | Pa | flux_down_from_atmos |
| hfss | Surface upward sensible heat flux | W/m² | flux_up_to_atmos |
| hfls | Surface upward latent heat flux | W/m² | flux_up_to_atmos |
| evspsbl | Evaporation | kg/m²/s | flux_up_to_atmos |
| tslsi | Surface temperature over land or sea ice | K | flux_up_to_atmos |
| tos | Sea-surface temperature | K | flux_up_to_atmos |
| sic | Sea ice area fraction | 1 | diag_sic |

### Global integral fields

These fields produce globally averaged scalar time-series output and are registered with 
`register_global_diag_field` (only when `use_AM3_physics` is not defined).

| Field | Description | Units | Sent in |
|---|---|---|---|
| evspsbl | Evaporation | mm/day | flux_up_to_atmos |
| ts | Surface temperature | K | flux_up_to_atmos |
| tas | Near-surface air temperature | K | sfc_boundary_layer |
| tasl | Near-surface air temperature (land only) | K | sfc_boundary_layer |
| hfss | Surface upward sensible heat flux | W/m² | flux_up_to_atmos |
| hfls | Surface upward latent heat flux | W/m² | flux_up_to_atmos |
| rls | Net longwave surface radiation | W/m² | flux_up_to_atmos |

### Tiled land fields

These fields are registered on the land axes (only on land PEs) using `register_tiled_diag_field` or, 
in the legacy land path, `fms_diag_register_diag_field`.

| Field | Description | Units | Sent in |
|---|---|---|---|
| t_ref | Temperature at reference height over land | K | sfc_boundary_layer |
| q_ref | Specific humidity at reference height over land | kg/kg | sfc_boundary_layer |
| rh_ref | Relative humidity at reference height over land | % | sfc_boundary_layer |
| u_ref | Zonal wind component at reference height over land | m/s | sfc_boundary_layer |
| v_ref | Meridional wind component at reference height over land | m/s | sfc_boundary_layer |
| evap | Evaporation rate over land | kg/m²/s | flux_up_to_atmos |
| shflx | Sensible heat flux over land | W/m² | flux_up_to_atmos |
| tasLut | Near-surface air temperature on land-use tile (reference height above displacement height) | K | sfc_boundary_layer |
| hussLut | Near-surface specific humidity on land-use tile | 1 | sfc_boundary_layer |
| *_tot_con_atm | Deposition velocity of tracer (new-land path only) | m/s | flux_up_to_atmos |
| *_tot_con_ref | Deposition velocity of tracer at reference height (new-land path only) | m/s | flux_up_to_atmos |
| *_flux | Surface flux of tracer over land | tracer units · kg air/(m²·s) | flux_up_to_atmos |
| *_mol_flux | Molar flux of tracer over land | mol CO₂/(m²·s) or mol/(m²·s) | flux_up_to_atmos |
| *_ref | Tracer mixing ratio at reference height over land (skipped for sphum; new-land path only) | tracer units | sfc_boundary_layer |

---

## Required Variables in Component Data Types

The sections below list the minimum fields each component must provide to the flux exchange. 
See the linked data-type documentation for full field descriptions.

### Atmosphere (`atmos_data_type`)

```fortran
type(atmos_data_type) :: Atm
real, dimension(:)    :: Atm%lon_bnd, Atm%lat_bnd   ! grid-box corner coordinates [radians]
real, dimension(:,:)  :: Atm%t_bot, Atm%q_bot, Atm%z_bot, Atm%p_bot, &
                         Atm%u_bot, Atm%v_bot,      &  ! lowest-level state
                         Atm%p_surf, Atm%slp, Atm%gust, &
                         Atm%flux_sw, Atm%flux_lw, Atm%lprec, Atm%fprec, Atm%coszen
integer, dimension(4) :: Atm%axes                   ! diag_manager axis IDs: x, y, z_full, z_half
```

- `Atm%lon_bnd`, `Atm%lat_bnd`: Grid-box corner coordinates in radians; must be monotonic.
- `Atm%t_bot` … `Atm%v_bot`: State at the lowest atmospheric model level; primary inputs to `sfc_boundary_layer`.
- `Atm%p_surf`, `Atm%slp`, `Atm%gust`: Surface pressure, sea-level pressure, and gustiness factor.
- `Atm%flux_sw`, `Atm%flux_lw`, `Atm%lprec`, `Atm%fprec`, `Atm%coszen`: Surface radiative and precipitation forcing passed to land and ice by `flux_down_from_atmos`.
- `Atm%axes`: Axis IDs returned by `diag_axis_init`; required for diagnostic registration.

The following fields support implicit time-stepping between the atmosphere and the surface models:

```fortran
type(surf_diff_type)  :: Atm%Surf_Diff
real, dimension(:,:)  :: Atm%Surf_Diff%dtmass,  &  ! dt/mass [s·m²/kg]
                         Atm%Surf_Diff%delta_t,  &  ! forward-elim temperature coefficient [K]
                         Atm%Surf_Diff%delta_q,  &  ! forward-elim moisture coefficient [kg/kg]
                         Atm%Surf_Diff%dflux_t,  &  ! d(heat flux)/d(T_surf) [W/m²/K]
                         Atm%Surf_Diff%dflux_q      ! d(moisture flux)/d(q_surf) [kg/m²/s / (kg/kg)]
```

### Land (`land_data_type`)

```fortran
type(land_data_type)        :: Land
real, dimension(:)          :: Land%lon_bnd, Land%lat_bnd
logical, dimension(:,:,:)   :: Land%mask, Land%glacier
real, dimension(:,:,:)      :: Land%tile_size, Land%t_surf, Land%t_ca, Land%q_ca, &
                               Land%albedo, Land%rough_mom, Land%rough_heat,       &
                               Land%stomatal, Land%snow, Land%water, Land%max_water
```

- `Land%lon_bnd`, `Land%lat_bnd`: Grid-box corner coordinates in radians; must be monotonic.
- `Land%mask`: Land-sea mask; `.true.` over land points.
- `Land%glacier`: Glacier mask; `.true.` over glacier points.
- `Land%tile_size`: Fractional area of each land tile within the atmospheric grid cell [0–1].
- `Land%t_surf`, `Land%albedo`, `Land%rough_mom`, `Land%rough_heat`: Surface state used for turbulent flux and radiation calculations.
- `Land%t_ca`, `Land%q_ca`: Canopy air temperature and specific humidity; returned to the atmosphere by `flux_up_to_atmos`.
- `Land%stomatal`, `Land%snow`, `Land%water`, `Land%max_water`: Additional surface properties used in flux parameterizations.

### Ice (`ice_data_type`)

```fortran
type(ice_data_type)         :: Ice
real, dimension(:)          :: Ice%lon_bnd, Ice%lat_bnd, Ice%lon_bnd_uv, Ice%lat_bnd_uv
logical, dimension(:,:,:)   :: Ice%mask, Ice%mask_uv, Ice%ice_mask
real, dimension(:,:,:)      :: Ice%part_size, Ice%part_size_uv
```

- All boundary arrays are in radians and must be monotonic.
- `Ice%mask`, `Ice%mask_uv`: Ocean-land masks for tracer and momentum points, respectively.
- `Ice%ice_mask`: Optional explicit sea-ice mask.
- `Ice%part_size`, `Ice%part_size_uv`: Fractional area of each ice thickness category.

Fields on the ice **top** (atmosphere–ice interface), provided to the atmosphere each fast timestep:

```fortran
real, dimension(:,:,:) :: Ice%t_surf, Ice%albedo, Ice%rough_mom, &
                          Ice%rough_heat, Ice%rough_moist, Ice%u_surf, Ice%v_surf
```

Fields on the ice **bottom** (ice–ocean interface), populated by `flux_down_from_atmos` and `flux_land_to_ice`, then passed to the ocean by `flux_ice_to_ocean`:

```fortran
real, dimension(:,:,:) :: Ice%flux_u,          &  ! zonal wind stress [Pa]
                          Ice%flux_v,          &  ! meridional wind stress [Pa]
                          Ice%flux_t,          &  ! sensible heat flux [W/m²]
                          Ice%flux_q,          &  ! moisture flux [kg/m²/s]
                          Ice%flux_salt,       &  ! salt flux [kg/m²/s]
                          Ice%flux_lw,         &  ! net longwave flux [W/m²]
                          Ice%flux_sw_vis_dir, &  ! direct visible SW [W/m²]
                          Ice%flux_sw_vis_dif, &  ! diffuse visible SW [W/m²]
                          Ice%flux_sw_nir_dir, &  ! direct near-IR SW [W/m²]
                          Ice%flux_sw_nir_dif, &  ! diffuse near-IR SW [W/m²]
                          Ice%lprec,           &  ! liquid precipitation [kg/m²/s]
                          Ice%fprec,           &  ! frozen precipitation [kg/m²/s]
                          Ice%runoff,          &  ! liquid runoff from land [kg/m²/s]
                          Ice%calving,         &  ! solid discharge from land [kg/m²/s]
                          Ice%runoff_hflx,     &  ! heat flux carried by runoff [W/m²]
                          Ice%calving_hflx,    &  ! heat flux carried by calving [W/m²]
                          Ice%p_surf              ! atmospheric surface pressure [Pa]
```

Optional iceberg fields, allocated only when the iceberg module is active:

```fortran
real, dimension(:,:) :: Ice%ustar_berg, &  ! iceberg friction velocity [m/s]
                        Ice%area_berg,  &  ! iceberg area fraction [dimensionless]
                        Ice%mass_berg      ! iceberg mass per unit area [kg/m²]
```

### Ocean (`ocean_public_type`)

```fortran
type(ocean_public_type) :: Ocean
real, dimension(:,:)    :: Ocean%t_surf, Ocean%s_surf, &
                           Ocean%u_surf, Ocean%v_surf, &
                           Ocean%frazil, Ocean%sea_lev
```

- `Ocean%t_surf`, `Ocean%s_surf`: Sea-surface temperature [K] and salinity [ppt]; passed to the ice model as `Ocean_Ice_Boundary%t` and `%s`.
- `Ocean%u_surf`, `Ocean%v_surf`: Surface ocean currents [m/s]; passed to ice as `Ocean_Ice_Boundary%u`, `%v`. Staggering relative to tracer points is indicated by `Ocean%stagger`.
- `Ocean%frazil`: Accumulated heating from frazil ice formation since the last coupling step [J/m²]; passed to ice as `Ocean_Ice_Boundary%frazil`.
- `Ocean%sea_lev`: Sea level corrected for surface pressure [m]; passed to ice as `Ocean_Ice_Boundary%sea_level`.
