# `atmos_land_boundary_type`

`atmos_land_boundary_type` carries all data passed from the coupler to the land model (LM4).
All fields are pointers with dimension (grid_idex, tile_number) 

## Radiative and precipitation forcing

| Field | Type | Description |
|---|---|---|
| Atmos_land_boundary%t_flux | real(:,:) | Sensible heat flux into the land surface [W/m²] |
| Atmos_land_boundary%lw_flux | real(:,:) | Net longwave radiation flux at the land surface [W/m²] |
| Atmos_land_boundary%lwdn_flux | real(:,:) | Downward longwave radiation flux at the land surface [W/m²] |
| Atmos_land_boundary%sw_flux | real(:,:) | Net shortwave radiation flux at the land surface [W/m²] |
| Atmos_land_boundary%swdn_flux | real(:,:) | Downward shortwave radiation flux at the land surface [W/m²] |
| Atmos_land_boundary%sw_flux_down_vis_dir | real(:,:) | Downward direct-beam visible shortwave flux [W/m²] |
| Atmos_land_boundary%sw_flux_down_total_dir | real(:,:) | Downward direct-beam total (broadband) shortwave flux [W/m²] |
| Atmos_land_boundary%sw_flux_down_vis_dif | real(:,:) | Downward diffuse visible shortwave flux [W/m²] |
| Atmos_land_boundary%sw_flux_down_total_dif | real(:,:) | Downward diffuse total (broadband) shortwave flux [W/m²] |
| Atmos_land_boundary%lprec | real(:,:) | Liquid precipitation rate [kg/m²/s] |
| Atmos_land_boundary%fprec | real(:,:) | Frozen precipitation rate [kg/m²/s] |
| Atmos_land_boundary%tprec | real(:,:) | Temperature of precipitation [K] |

## Implicit coupling coefficients

These are the derivatives needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and land.

| Field | Type | Description |
|---|---|---|
| Atmos_land_boundary%dhdt | real(:,:) | d(sensible heat flux)/d(T_surf) — derivative of sensible heat flux with respect to surface temperature [W/m²/K] |
| Atmos_land_boundary%dhdq | real(:,:) | d(sensible heat flux)/d(q_surf) — derivative of sensible heat flux with respect to surface specific humidity [W/m²/(kg/kg)] |
| Atmos_land_boundary%drdt | real(:,:) | d(longwave flux)/d(T_surf) — derivative of longwave flux with respect to surface radiative temperature [W/m²/K] |

## Turbulence and boundary-layer state

| Field | Type | Description |
|---|---|---|
| Atmos_land_boundary%cd_m | real(:,:) | Drag coefficient for momentum [dimensionless] |
| Atmos_land_boundary%cd_t | real(:,:) | Drag coefficient for tracers (heat and moisture) [dimensionless] |
| Atmos_land_boundary%ustar | real(:,:) | Turbulent wind scale (friction velocity) [m/s] |
| Atmos_land_boundary%bstar | real(:,:) | Turbulent buoyancy scale [m/s] |
| Atmos_land_boundary%wind | real(:,:) | Absolute wind speed at the bottom of the atmospheric layer [m/s] |
| Atmos_land_boundary%z_bot | real(:,:) | Height of the bottom atmospheric layer above the land surface [m] |
| Atmos_land_boundary%drag_q | real(:,:) | Product of the moisture drag coefficient and wind speed (cd_q × wind); used in land surface moisture flux calculations [m/s] |
| Atmos_land_boundary%p_surf | real(:,:) | Surface pressure [Pa] |

## Tracer fluxes
dimension (grid_index, tile_number, tracer_index)

| Field | Type | Description |
|---|---|---|
| Atmos_land_boundary%tr_flux | real(:,:,:) | Flux of each tracer into the land surface, including water vapor flux; dimensioned (grid_index, tile, tracer) [tracer units · kg air / (m²·s)] |
| Atmos_land_boundary%dfdtr | real(:,:,:) | d(tracer flux)/d(tracer_surf) — derivative of the tracer flux with respect to the surface tracer value, including evaporation over surface specific humidity; dimensioned (grid_index, tile, tracer) |

## Bookkeeping

| Field | Type | Description |
|---|---|---|
| Atmos_land_boundary%xtype | integer | Transfer mode for the atmosphere-to-land exchange: REGRID (1), REDIST (2), or DIRECT (3) |
