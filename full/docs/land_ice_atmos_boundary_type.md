# `land_ice_atmos_boundary_type`

`land_ice_atmos_boundary_type` contains surface quantities going from land and ice to the atmosphere.
All quantities are on the exchange grid. 

## Surface state returned to the atmosphere

| Field | Type | Description |
|---|---|---|
| Land_ice_atmos_boundary%t | real(:,:) | Area-weighted surface temperature seen by the atmosphere for radiation calculations [K]; weighted over land and ice fractions |
| Land_ice_atmos_boundary%t_ocean | real(:,:) | Ocean surface temperature for radiation calculations; sourced from Ice%t_surf through the exchange grid [K] |
| Land_ice_atmos_boundary%albedo | real(:,:) | Broadband surface albedo [dimensionless] |
| Land_ice_atmos_boundary%albedo_vis_dir | real(:,:) | Direct-beam visible-band surface albedo [dimensionless] |
| Land_ice_atmos_boundary%albedo_nir_dir | real(:,:) | Direct-beam near-infrared surface albedo [dimensionless] |
| Land_ice_atmos_boundary%albedo_vis_dif | real(:,:) | Diffuse visible-band surface albedo [dimensionless] |
| Land_ice_atmos_boundary%albedo_nir_dif | real(:,:) | Diffuse near-infrared surface albedo [dimensionless] |
| Land_ice_atmos_boundary%land_frac | real(:,:) | Fraction of the atmospheric grid cell covered by land [dimensionless] |
| Land_ice_atmos_boundary%frac_open_sea | real(:,:) | Non-sea-ice fraction of the grid cell [dimensionless]; complement of the sea-ice concentration |
| Land_ice_atmos_boundary%rough_mom | real(:,:) | Area-weighted surface roughness length for momentum [m] |
| Land_ice_atmos_boundary%rough_heat | real(:,:) | Area-weighted surface roughness length for heat [m] |

## Reference-height diagnostics

These fields are interpolated from the surface-layer profile to the reference heights z_ref_heat and z_ref_mom configured in flux_exchange_nml.

| Field | Type | Description |
|---|---|---|
| Land_ice_atmos_boundary%u_ref | real(:,:) | Zonal wind at the momentum reference height (z_ref_mom) [m/s] |
| Land_ice_atmos_boundary%v_ref | real(:,:) | Meridional wind at the momentum reference height (z_ref_mom) [m/s] |
| Land_ice_atmos_boundary%t_ref | real(:,:) | Air temperature at the heat reference height (z_ref_heat) [K] |
| Land_ice_atmos_boundary%q_ref | real(:,:) | Specific humidity at the heat reference height (z_ref_heat) [kg/kg] |
| Land_ice_atmos_boundary%wind | real(:,:) | Absolute wind speed at the lowest atmospheric model level including gust corrections [m/s] |
| Land_ice_atmos_boundary%thv_atm | real(:,:) | Virtual potential temperature at the lowest atmospheric model level [K] |
| Land_ice_atmos_boundary%thv_surf | real(:,:) | Virtual potential temperature at the surface [K] |

## Implicit flux corrections applied to the atmosphere

These fields are the outputs of the tridiagonal back-substitution step. They correct the atmospheric temperature and tracer tendencies for the implicit surface coupling.

| Field | Type | Description |
|---|---|---|
| Land_ice_atmos_boundary%dt_t | real(:,:) | Temperature tendency correction at the lowest atmospheric level from the implicit surface flux scheme [K/s] |
| Land_ice_atmos_boundary%dt_tr | real(:,:,:) | Tracer mixing-ratio tendency correction at the lowest level; third dimension indexes tracers [tracer units/s] |

## Surface stress and turbulence scales

| Field | Type | Description |
|---|---|---|
| Land_ice_atmos_boundary%u_flux | real(:,:) | Zonal wind stress on the atmosphere [Pa] |
| Land_ice_atmos_boundary%v_flux | real(:,:) | Meridional wind stress on the atmosphere [Pa] |
| Land_ice_atmos_boundary%dtaudu | real(:,:) | d(zonal wind stress)/d(u) — implicit coupling coefficient for zonal momentum [Pa·s/m] |
| Land_ice_atmos_boundary%dtaudv | real(:,:) | d(meridional wind stress)/d(v) — implicit coupling coefficient for meridional momentum [Pa·s/m] |
| Land_ice_atmos_boundary%u_star | real(:,:) | Friction velocity (surface turbulent velocity scale) [m/s] |
| Land_ice_atmos_boundary%b_star | real(:,:) | Buoyancy scale used in Monin-Obukhov similarity theory [m/s²] |
| Land_ice_atmos_boundary%q_star | real(:,:) | Moisture scale used in Monin-Obukhov similarity theory [kg/kg] |
| Land_ice_atmos_boundary%shflx | real(:,:) | Sensible heat flux at the surface [W/m²]; not compiled when use_AM3_physics is defined |
| Land_ice_atmos_boundary%lhflx | real(:,:) | Latent heat flux at the surface [W/m²]; not compiled when use_AM3_physics is defined |

## Generic and bookkeeping fields

| Field | Type | Description |
|---|---|---|
| Land_ice_atmos_boundary%data | real(:,:,:) | Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations |
| Land_ice_atmos_boundary%gex_lnd2atm | real(:,:,:) | Generic exchange fields returned from the land model to the atmosphere (e.g., surface emission fluxes); third dimension indexes the exchange field list |
| Land_ice_atmos_boundary%xtype | integer | Transfer mode for the exchange-grid-to-atmosphere remap: REGRID (1), REDIST (2), or DIRECT (3) |
