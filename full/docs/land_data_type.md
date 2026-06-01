# `land_data_type`

`land_data_type` carries land surface states passed to the couper.
The land model uses an unstructured tile representation where
multiple land-cover types (soil, vegetation,
lakes, etc.) can co-exist within a single atmospheric grid cell.

---

## Surface state fields
dimension (grid_index,tile_number)

| Field | Description |
|---|---|
| Land%tile_size | Fractional coverage of the atmospheric grid cell by this tile [dimensionless, 0–1]; used to area-weight tile quantities back onto the atmosphere grid |
| Land%t_surf | Ground (radiative) surface temperature [K]; used in longwave radiation and sensible heat flux calculations |
| Land%t_ca | Canopy air temperature — near-surface air temperature within the plant canopy layer [K]; differs from t_surf over vegetated tiles |
| Land%albedo | Broadband surface albedo [dimensionless]; legacy field, per-band albedos below are preferred |
| Land%albedo_vis_dir | Surface albedo for direct-beam visible radiation (0.2–0.7 µm) [dimensionless] |
| Land%albedo_nir_dir | Surface albedo for direct-beam near-infrared radiation [dimensionless] |
| Land%albedo_vis_dif | Surface albedo for diffuse visible radiation [dimensionless] |
| Land%albedo_nir_dif | Surface albedo for diffuse near-infrared radiation [dimensionless] |
| Land%rough_mom | Surface roughness length for momentum [m]; used in Monin-Obukhov flux calculations |
| Land%rough_heat | Surface roughness length for heat and tracers [m] |
| Land%rough_scale | Topographic form-drag scaling factor for momentum [dimensionless]; accounts for sub-grid orographic drag |

---

## Tracer fields
dimension (grid_index, tile_number, tracer_index)

| Field | Description |
|---|---|
| Land%tr | Surface tracer mixing ratios on each tile, including canopy air specific humidity as the first tracer; additional tracers (e.g. CO₂) follow the tracer table order |

---

## Discharge fields
dimension(lon, lat)


These fields carry freshwater and heat leaving the land surface and routed
to the ocean/ice via `flux_land_to_ice`.

| Field | Description |
|---|---|
| Land%discharge | Liquid water discharge (river runoff) from land to ocean [kg/m²/s] |
| Land%discharge_heat | Sensible heat carried by liquid discharge, using 0 °C as datum [W/m²] |
| Land%discharge_snow | Solid water (snow/ice) discharge from land to ocean [kg/m²/s] |
| Land%discharge_snow_heat | Sensible heat carried by solid discharge, using 0 °C as datum [W/m²] |

---

## Grid and domain bookkeeping

| Field | Type | Description |
|---|---|---|
| Land%mask | logical(:,:) | .true. where the grid cell contains land; used to gate land-only computations |
| Land%axes(1) | integer | Diag-manager axis ID for the unstructured land grid; used when registering tiled land diagnostics |
| Land%domain | type(domain2D) | FMS structured-grid domain for the land model; used for halo exchanges and exchange-grid setup |
| Land%ug_domain | type(domainUG) | FMS unstructured-grid domain for the land model; carries the tile-based decomposition used by LM4 |
| Land%pelist | integer(:) | List of MPI PE numbers on which the land model is running |
| Land%pe | logical | .true. on PEs that are part of the land pelist; used to gate land stock calculations |

