# `atmos_ice_boundary_type`

Defined in `sis2/src/ice_boundary_types.F90:51`.

`atmos_ice_boundary_type` holds data passed between atmosphere and the sea-ice model (SIS2).

---

## Wind stress

| Field | Type | Description |
|---|---|---|
| u_flux | real(:,:,:) | True-eastward wind stress from the atmosphere to the ocean or ice in each thickness category, on an A-grid and not rotated to the model grid [Pa] |
| v_flux | real(:,:,:) | True-northward wind stress from the atmosphere to the ocean or ice in each thickness category, on an A-grid and not rotated to the model grid [Pa] |
| u_star | real(:,:,:) | Atmospheric friction velocity on an A-grid [Pa] |

---

## Turbulent heat and moisture fluxes

| Field | Type | Description |
|---|---|---|
| t_flux | real(:,:,:) | Net sensible heat flux from the ocean or ice surface into the atmosphere [W/m²] |
| q_flux | real(:,:,:) | Moisture flux from the ice or ocean to the atmosphere due to evaporation or sublimation [kg/m²/s] |

---

## Implicit coupling coefficients

These are the linearization (derivative) terms needed to close the implicit tridiagonal surface diffusion scheme between the atmosphere and ice.

| Field | Type | Description |
|---|---|---|
| dhdt | real(:,:,:) | d(upward sensible heat flux)/d(T_surf) — derivative of sensible heat flux with respect to surface temperature [W/m²/°C] |
| dedt | real(:,:,:) | d(sublimation+evaporation rate)/d(T_surf) — derivative of the moisture flux with respect to surface temperature [kg/m²/s/°C] |
| drdt | real(:,:,:) | d(net upward longwave flux)/d(T_surf) — derivative of the net upward longwave flux (i.e. -lw_flux) with respect to surface temperature [W/m²/°C] |

---

## Radiative fluxes

| Field | Type | Description |
|---|---|---|
| lw_flux | real(:,:,:) | Net downward longwave radiation flux from the atmosphere into the ice or ocean [W/m²] |
| sw_flux_vis_dir | real(:,:,:) | Net direct visible shortwave radiation flux into the ice or ocean [W/m²] |
| sw_flux_vis_dif | real(:,:,:) | Net diffuse visible shortwave radiation flux into the ice or ocean [W/m²] |
| sw_flux_nir_dir | real(:,:,:) | Net direct near-infrared shortwave radiation flux into the ice or ocean [W/m²] |
| sw_flux_nir_dif | real(:,:,:) | Net diffuse near-infrared shortwave radiation flux into the ice or ocean [W/m²] |
| sw_down_vis_dir | real(:,:,:) | Downward direct visible shortwave radiation flux from the atmosphere [W/m²] |
| sw_down_vis_dif | real(:,:,:) | Downward diffuse visible shortwave radiation flux from the atmosphere [W/m²] |
| sw_down_nir_dir | real(:,:,:) | Downward direct near-infrared shortwave radiation flux from the atmosphere [W/m²] |
| sw_down_nir_dif | real(:,:,:) | Downward diffuse near-infrared shortwave radiation flux from the atmosphere [W/m²] |
| coszen | real(:,:,:) | Cosine of the solar zenith angle averaged over the next radiation timestep (not the timestep used to compute the sw_flux fields) [dimensionless, ≤ 1] |

---

## Precipitation

| Field | Type | Description |
|---|---|---|
| lprec | real(:,:,:) | Liquid precipitation (rain) from the atmosphere onto the ice or ocean in each thickness category [kg/m²/s]; rain falling on snow is currently assumed to drain directly through the ice into the ocean |
| fprec | real(:,:,:) | Frozen precipitation (snowfall, sleet, hail, graupel) from the atmosphere to the ice or ocean [kg/m²/s]; all forms of frozen precipitation are treated as snow in SIS2 |

---

## Atmospheric state

| Field | Type | Description |
|---|---|---|
| p | real(:,:,:) | Atmospheric surface pressure [Pa]; typically ~1×10⁵ Pa |

---

## Bookkeeping

| Field | Type | Description |
|---|---|---|
| xtype | integer | Transfer mode for the atmosphere-to-ice exchange: REGRID (1), REDIST (2), or DIRECT (3) |
| fluxes | type(coupler_3d_bc_type) | Array of additional per-tracer gas and deposition fluxes from the atmosphere to the ice |
