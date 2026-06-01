# `ocean_ice_boundary_type`

`ocean_ice_boundary_type` carries ocean surface state passed from MOM6 to the sea-ice model (SIS2)

---

## Ocean surface velocities

| Field | Type | Description |
|---|---|---|
| u | real(:,:) | x-direction ocean surface velocity at a position determined by stagger [m/s] |
| v | real(:,:) | y-direction ocean surface velocity at a position determined by stagger [m/s] |
| stagger | integer | Spatial staggering of u and v relative to tracer points; default is BGRID_NE |

---

## Thermodynamic state

| Field | Type | Description |
|---|---|---|
| t | real(:,:) | Ocean surface temperature [K] |
| s | real(:,:) | Ocean surface salinity [g salt / kg seawater] |
| frazil | real(:,:) | Frazil heat rejected by the ocean since the last coupling step [J/m²]; delivered to SIS2 so it can account for ocean-side freezing |
| sea_level | real(:,:) | Sea level after adjustment for any surface pressure that the ocean allows to be expressed [m] |

---

## Calving

| Field | Type | Description |
|---|---|---|
| calving | real(:,:) | Mass flux per unit area of ice-shelf flux to be converted to icebergs [kg/m²/s] |
| calving_hflx | real(:,:) | Heat flux associated with calving [W/m²] |

---

## Bookkeeping

| Field | Type | Description |
|---|---|---|
| data | real(:,:,:) | Collective array providing named access to the scalar fields above; used internally for data-override and exchange-grid operations |
| xtype | integer | Transfer mode for the ocean-to-ice exchange: REGRID (1), REDIST (2), or DIRECT (3) |
| fields | type(coupler_2d_bc_type) | Named array of additional per-tracer ocean surface fields (e.g., pCO₂, SSS for gas exchange) passed to the ice model |
