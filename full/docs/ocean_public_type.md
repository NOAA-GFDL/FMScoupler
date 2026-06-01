# `ocean_public_type`

`ocean_public_type` is the publicly visible face of the MOM6 ocean model.
It contains only the surface fields and domain metadata for the coupler.

---

## Domain and PE bookkeeping

| Field | Type | Description |
|---|---|---|
| Domain | type(domain2d) | FMS domain decomposition for the ocean surface fields; defines the MPI tile layout used for coupler remapping |
| is_ocean_pe | logical | .true. on processors that run the ocean model; used throughout the coupler to gate ocean-only code paths |
| pelist | integer(:) | List of MPI PE numbers assigned to the ocean model |
| maskmap | logical(:,:) | Pointer to a mask indicating which logical processors are active for ocean computation; logical processors covering all-land points may not be mapped to physical PEs. Need not be set if all processors are used |
| instance_name | character(32) | Optional name identifying this ocean model instance; used in ensemble runs to disambiguate log messages |

---

## Surface state fields

These fields are set by the ocean model after each update and read by the
coupler to construct the `ocean_ice_boundary_type` passed to the sea-ice
model via `flux_ocean_to_ice`.

| Field | Description |
|---|---|
| t_surf | Sea-surface temperature (SST) on tracer (T) cells [K] |
| s_surf | Sea-surface salinity (SSS) on T cells [ppt] |
| u_surf | Surface current i-velocity at the locations indicated by stagger [m/s] |
| v_surf | Surface current j-velocity at the locations indicated by stagger [m/s] |
| sea_lev | Sea level corrected for surface pressure: dzt(1) + η + p_atm/(ρ₀g) [m]; passed to the ice model as sea_level |
| frazil | Accumulated heating from frazil ice formation in the ocean since the last coupling step [J/m²]; delivered to the ice model so it can account for ocean-side freezing |
| melt_potential | Instantaneous heat available to melt sea ice from below [J/m²]; computed when the ocean boundary layer depth exceeds HFrz |
| OBLD | Ocean boundary layer depth [m]; used to determine the depth over which melt potential is computed |
| area | Grid-cell area of each ocean surface cell [m²]; used for conservative flux remapping |
| calving | Mass per unit area of ice-shelf flux to be converted to icebergs [kg/m²]; passed to the iceberg module |
| calving_hflx | Heat flux associated with calving [W/m²] |

---

## Grid staggering

| Field | Type | Description |
|---|---|---|
| stagger | integer | Arakawa staggering of the surface velocity components (u_surf, v_surf) relative to tracer points. Valid values: AGRID, BGRID_NE, CGRID_NE, BGRID_SW, CGRID_SW. Set to -999 before initialization so a global max can propagate the value to non-ocean PEs |

---

## Tracer gas fields and diagnostics

| Field | Type | Description |
|---|---|---|
| fields | type(coupler_2d_bc_type) | Named arrays of tracer-related ocean surface fields (e.g. pCO₂, O₂ saturation) used in atmosphere-ocean gas flux calculations; populated by atmos_ocean_fluxes_calc |
| avg_kount | integer | Counter tracking the number of contributions accumulated in the running averages stored in this type; used externally by FMSCoupler to manage time averaging |
| axes(2) | integer(2) | Diag-manager axis IDs available for I/O using this surface data |
