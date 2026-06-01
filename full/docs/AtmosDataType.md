# `atmos_data_type`
is the main derived type holding fields and states of the atmosphere.  Below are all the fields in atmos_data_type:

## Grid metadata

| Field | Type | Description |
|---|---|---|
| Atm%domain | type(domain2d) | FMS domain decomposition object for the atmosphere; defines the MPI tile layout and halo widths |
| Atm%axes | integer(4) | Diag-manager axis indices for x, y, pfull, and phalf; used when registering and sending diagnostic fields |
| Atm%lon_bnd | real(:,:) | Longitude of grid-box corners on the local compute domain [radians] |
| Atm%lat_bnd | real(:,:) | Latitude of grid-box corners on the local compute domain [radians] |
| Atm%lon | real(:,:) | Longitude of grid-box centres on the local compute domain [radians] |
| Atm%lat | real(:,:) | Latitude of grid-box centres on the local compute domain [radians] |
| Atm%grid | type(grid_box_type) | Grid geometry needed for second-order conservative remapping on the cubic-sphere exchange grid (see below) |
| Atm%maskmap | logical(:,:) | Pointer to a mask indicating which logical processors are active for ocean code; processors covering all-land points may not be assigned to physical PEs. Dummy field — must be present for compilation but need not be set |

## Atmospheric state at the bottom-most level 

| Field | Type | Description |
|---|---|---|
| Atm%t_bot | real(:,:) | Temperature at the lowest model level [K] |
| Atm%tr_bot | real(:,:,:) | Tracer mixing ratios at the lowest model level; third dimension indexes the tracer table (specific humidity sphum is always present) |
| Atm%z_bot | real(:,:) | Height of the lowest model level above the surface [m] |
| Atm%p_bot | real(:,:) | Pressure at the lowest model level [Pa] |
| Atm%u_bot | real(:,:) | Zonal wind component at the lowest model level [m/s] |
| Atm%v_bot | real(:,:) | Meridional wind component at the lowest model level [m/s] |
| Atm%p_surf | real(:,:) | Surface pressure [Pa] |
| Atm%slp | real(:,:) | Sea-level pressure [Pa] |
| Atm%gust | real(:,:) | Gustiness factor — a minimum wind speed added in quadrature to the resolved wind to account for sub-grid convective gusts in surface flux calculations [m/s] |
| Atm%coszen | real(:,:) | Cosine of the solar zenith angle; used to weight shortwave fluxes and partition direct vs. diffuse radiation [dimensionless] |

### Shortwave and longwave fluxes

| Field | Type | Description |
|---|---|---|
| Atm%flux_sw | real(:,:) | Total net shortwave flux at the surface (absorbed by the surface) [W/m²] |
| Atm%flux_sw_dir | real(:,:) | Direct-beam component of the net shortwave flux [W/m²] |
| Atm%flux_sw_dif | real(:,:) | Diffuse component of the net shortwave flux [W/m²] |
| Atm%flux_sw_down_vis_dir | real(:,:) | Downward direct-beam flux in the visible band (0.2–0.7 µm) [W/m²] |
| Atm%flux_sw_down_vis_dif | real(:,:) | Downward diffuse flux in the visible band [W/m²] |
| Atm%flux_sw_down_total_dir | real(:,:) | Downward direct-beam broadband shortwave flux [W/m²] |
| Atm%flux_sw_down_total_dif | real(:,:) | Downward diffuse broadband shortwave flux [W/m²] |
| Atm%flux_sw_vis | real(:,:) | Net (downward minus reflected) visible-band shortwave flux at the surface [W/m²] |
| Atm%flux_sw_vis_dir | real(:,:) | Direct-beam component of the net visible shortwave flux [W/m²] |
| Atm%flux_sw_vis_dif | real(:,:) | Diffuse component of the net visible shortwave flux [W/m²] |
| Atm%flux_lw | real(:,:) | Net downward longwave flux at the surface [W/m²] |

### Precipitation
| Atm%lprec | real(:,:) | Mass of liquid precipitation accumulated since the last time step [kg/m²]; equivalent to a rate in kg/m²/s when divided by dt_atm |
| Atm%fprec | real(:,:) | Mass of frozen (solid) precipitation accumulated since the last time step [kg/m²] |

## Generic exchange fields

| Field | Type | Description |
|---|---|---|
| Atm%gex_atm2lnd | real(:,:,:) | Generic exchange fields sent from the atmosphere to the land model; the field table defines which quantities are exchanged (e.g. CO₂, aerosol deposition); third dimension indexes the exchange field list |
| Atm%gex_lnd2atm | real(:,:,:) | Generic exchange fields returned from the land model to the atmosphere (e.g. surface emission fluxes); third dimension indexes the exchange field list |
| Atm%fields | type(coupler_2d_bc_type) | Array of additional tracer boundary-condition fields used for atmosphere-ocean gas exchange (CO₂, O₂, CFCs, etc.); registered and populated by atmos_tracer_flux_init |

---

## Implicit surface diffusion coefficients — `Surf_diff`

`Surf_diff` is of type `surf_diff_type` (defined in
`atmos_phys/atmos_param/vert_diff/vert_diff.F90:44`).  It carries the
forward-elimination coefficients from the implicit vertical diffusion
scheme that couples the atmosphere to the surface models. 

| Field | Type | Description |
|---|---|---|
| Atm%Surf_diff%dtmass | real(:,:) | dt/mass — ratio of the atmospheric timestep to the surface-layer air mass [s·m²/kg]; scales flux tendencies to temperature/tracer tendencies |
| Atm%Surf_diff%dflux_t | real(:,:) | d(sensible heat flux)/d(T_surf) — linearisation of the surface heat flux with respect to surface temperature; used to form the implicit coupling term [W/m²/K] |
| Atm%Surf_diff%delta_t | real(:,:) | Forward-elimination coefficient for temperature from the implicit tridiagonal scheme; represents the accumulated atmospheric temperature forcing at the bottom level waiting for the surface response [K] |
| Atm%Surf_diff%delta_u | real(:,:) | Forward-elimination coefficient for zonal wind from the implicit scheme [m/s] |
| Atm%Surf_diff%delta_v | real(:,:) | Forward-elimination coefficient for meridional wind from the implicit scheme [m/s] |
| Atm%Surf_diff%dflux_tr | real(:,:,:) | d(tracer flux)/d(tracer_surf) — linearisation of tracer surface fluxes with respect to surface tracer concentration; third dimension indexes tracers |
| Atm%Surf_diff%delta_tr | real(:,:,:) | Forward-elimination coefficient for each tracer from the implicit scheme; third dimension indexes tracers |
| Atm%Surf_diff%tdt_dyn | real(:,:,:) | Temperature tendency from dynamics (advection, etc.) passed through the diffusion scheme |
| Atm%Surf_diff%qdt_dyn | real(:,:,:) | Moisture tendency from dynamics |
| Atm%Surf_diff%dgz_dyn | real(:,:,:) | Geopotential height tendency from dynamics |
| Atm%Surf_diff%ddp_dyn | real(:,:,:) | Pressure-thickness tendency from dynamics |
| Atm%Surf_diff%tdt_rad | real(:,:,:) | Temperature tendency from radiation; used in the MIZ (marginal ice zone) forecast mode |

---

## Grid geometry — `grid`

`grid` is of type `grid_box_type` (defined in `FMS/exchange/xgrid.F90:287`).
It holds the geometric quantities needed for second-order conservative flux
remapping between the atmosphere and surface component grids on a
cubic-sphere mesh.

| Field | Type | Description |
|---|---|---|
| Atm%grid%dx | real(:,:) | Grid-box width in the x-direction [m] |
| Atm%grid%dy | real(:,:) | Grid-box width in the y-direction [m] |
| Atm%grid%area | real(:,:) | Grid-box area [m²] |
| Atm%grid%edge_w | real(:) | Western edge lengths of grid boxes along the boundary [m] |
| Atm%grid%edge_e | real(:) | Eastern edge lengths [m] |
| Atm%grid%edge_s | real(:) | Southern edge lengths [m] |
| Atm%grid%edge_n | real(:) | Northern edge lengths [m] |
| Atm%grid%en1 | real(:,:,:) | First unit normal vector at grid-box edges; used to project vector fields (winds, stresses) during remapping |
| Atm%grid%en2 | real(:,:,:) | Second unit normal vector at grid-box edges |
| Atm%grid%vlon | real(:,:,:) | Unit vector in the local longitude direction at each grid point; used to rotate between geographic and local coordinate frames during exchange |
| Atm%grid%vlat | real(:,:,:) | Unit vector in the local latitude direction at each grid point |

## Bookkeeping and time

| Field | Type | Description |
|---|---|---|
| Atm%Time | type(time_type) | Current model time; passed to diag_manager send_data calls and to fms_data_override |
| Atm%Time_step | type(time_type) | Atmospheric model timestep duration |
| Atm%Time_init | type(time_type) | Reference (initial) time for the model run |
| Atm%pelist | integer(:) | List of MPI PE numbers on which the atmosphere is running |
| Atm%pe | logical | .true. on PEs that are part of the atmosphere pelist; used to gate atmosphere-only code blocks |
