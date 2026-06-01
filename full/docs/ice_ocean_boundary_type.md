# `ice_ocean_boundary_type`

Defined in two driver-specific files:
- FMS cap: `mom6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90:180`
- NUOPC cap: `mom6/config_src/drivers/nuopc_cap/mom_surface_forcing_nuopc.F90:166`

`ice_ocean_boundary_type` holds all surface forcing passed from the sea-ice model to MOM6 each coupled timestep.

The two cap implementations share a common set of core fields but differ in the additional fields supported. Fields present only in one cap are noted below.

---

## Wind stress

| Field | Type | Description |
|---|---|---|
| u_flux | real(:,:) | i-direction wind/ice stress on the ocean surface [Pa] |
| v_flux | real(:,:) | j-direction wind/ice stress on the ocean surface [Pa] |
| stress_mag | real(:,:) | Time-mean magnitude of the stress on the ocean [Pa]; present when pass_stress_mag=.true. in SIS_slow_CS. FMS cap only. |
| wind_stagger | integer | Spatial discretization of the wind stresses; may be set by the flux-exchange code based on what the sea-ice model provides, otherwise taken from the surface forcing control structure |
| u10_sqr | real(:,:) | Wind speed squared at 10 m height [m²/s²]. NUOPC cap only. |

---

## Heat fluxes

| Field | Type | Description |
|---|---|---|
| t_flux | real(:,:) | Sensible heat flux into the ocean [W/m²] |
| lw_flux | real(:,:) | Net longwave radiation flux into the ocean [W/m²] |
| sw_flux_vis_dir | real(:,:) | Direct visible shortwave radiation into the ocean [W/m²] |
| sw_flux_vis_dif | real(:,:) | Diffuse visible shortwave radiation into the ocean [W/m²] |
| sw_flux_nir_dir | real(:,:) | Direct near-infrared shortwave radiation into the ocean [W/m²] |
| sw_flux_nir_dif | real(:,:) | Diffuse near-infrared shortwave radiation into the ocean [W/m²] |
| seaice_melt_heat | real(:,:) | Heat flux from sea ice and snow melting [W/m²]. NUOPC cap only. |
| swnet_afracr | real(:,:) | Net shortwave radiation multiplied by the atmosphere fraction, positive into the ocean [W/m²]. NUOPC cap only. |
| swpen_ifrac_n | real(:,:,:) | Net shortwave radiation penetrating into ice and ocean, multiplied by ice fraction per thickness category; positive into the ocean [W/m²]; third dimension indexes ice categories. NUOPC cap only. |

---

## Heat content of freshwater fluxes

| Field | Type | Description |
|---|---|---|
| hrofl | real(:,:) | Heat content from liquid runoff [W/m²] |
| hrofi | real(:,:) | Heat content from frozen runoff (calving) [W/m²] |
| hrofl_glc | real(:,:) | Heat content from liquid glacier runoff via the river-routing model [W/m²] |
| hrofi_glc | real(:,:) | Heat content from frozen glacier runoff via the river-routing model [W/m²] |
| hrain | real(:,:) | Heat content from liquid precipitation [W/m²] |
| hsnow | real(:,:) | Heat content from frozen precipitation [W/m²] |
| hevap | real(:,:) | Heat content from evaporation [W/m²] |
| hcond | real(:,:) | Heat content from condensation [W/m²] |

---

## Freshwater and salt fluxes

| Field | Type | Description |
|---|---|---|
| q_flux | real(:,:) | Specific humidity (freshwater) flux into the ocean [kg/m²/s] |
| salt_flux | real(:,:) | Salt flux from sea ice into the ocean (brine rejection / melting) [kg/m²/s] |
| excess_salt | real(:,:) | Salt left behind in the ocean by brine rejection rather than ejected as a salt flux [kg/m²/s]. FMS cap only. |
| seaice_melt | real(:,:) | Water flux due to sea ice and snow melting [kg/m²/s]. NUOPC cap only. |
| lprec | real(:,:) | Mass flux of liquid precipitation into the ocean [kg/m²/s] |
| fprec | real(:,:) | Mass flux of frozen precipitation into the ocean [kg/m²/s] |

---

## Runoff and calving

The FMS cap uses `runoff` / `calving` naming; the NUOPC cap uses `lrunoff` / `frunoff` and additionally carries glacier-sourced fluxes.

| Field | Type | Cap | Description |
|---|---|---|---|
| runoff | real(:,:) | FMS | Mass flux of liquid runoff from land into the ocean [kg/m²/s] |
| runoff_carbon | real(:,:) | FMS | Mass flux of carbon carried by liquid runoff [kg/m²/s] |
| runoff_hflx | real(:,:) | FMS | Heat content of liquid runoff relative to 0 °C [W/m²] |
| calving | real(:,:) | FMS | Mass flux of frozen runoff (calving) into the ocean [kg/m²/s]; offered first to icebergs if active |
| calving_hflx | real(:,:) | FMS | Heat content of frozen runoff relative to 0 °C [W/m²] |
| lrunoff | real(:,:) | NUOPC | Liquid runoff [kg/m²/s] |
| frunoff | real(:,:) | NUOPC | Frozen (ice) runoff [kg/m²/s] |
| lrunoff_glc | real(:,:) | NUOPC | Liquid glacier runoff delivered via the river-routing model [kg/m²/s] |
| frunoff_glc | real(:,:) | NUOPC | Frozen glacier runoff delivered via the river-routing model [kg/m²/s] |

---

## Pressure and ice loading

| Field | Type | Description |
|---|---|---|
| p | real(:,:) | Pressure of overlying ice and atmosphere on the ocean surface [Pa] |
| mi | real(:,:) | Mass of sea ice per unit ocean area [kg/m²]; used for ice-pressure loading |
| ice_rigidity | real(:,:) | Rigidity of sea ice and ice shelves expressed as a divergence-damping coefficient [m³/s]; determined outside the ocean model |
| ice_fraction | real(:,:) | Fractional ice area [dimensionless]. NUOPC cap only. |
| ifrac_n | real(:,:,:) | Ice fraction per ice thickness category [dimensionless]; third dimension indexes categories. NUOPC cap only. |
| ice_ncat | integer | Number of ice categories provided by the coupler; 1 means per-category data is not used. NUOPC cap only. |
| afracr | real(:,:) | Fractional atmosphere coverage relative to the ocean grid cell [dimensionless]. NUOPC cap only. |

---

## Iceberg fields

Present only when the iceberg module is active.

| Field | Type | Description |
|---|---|---|
| ustar_berg | real(:,:) | Frictional velocity beneath icebergs [m/s] |
| area_berg | real(:,:) | Fractional area of the ocean cell covered by icebergs [m²/m²] |
| mass_berg | real(:,:) | Mass of icebergs per unit ocean area [kg/m²] |

---

## Ice-shelf fields

FMS cap only.

| Field | Type | Description |
|---|---|---|
| shelf_sfc_mass_flux | real(:,:) | Mass flux to the surface of the ice sheet [kg/m²/s] |

---

## Biogeochemistry and aerosol deposition

NUOPC cap only. These fields support ocean biogeochemistry modules that require atmospheric deposition forcing.

| Field | Type | Description |
|---|---|---|
| nhx_dep | real(:,:) | Reduced nitrogen (NHx) deposition flux [kg/m²/s] |
| noy_dep | real(:,:) | Oxidized nitrogen (NOy) deposition flux [kg/m²/s] |
| atm_co2_prog | real(:,:) | Prognostic atmospheric CO₂ concentration [ppm] |
| atm_co2_diag | real(:,:) | Diagnostic atmospheric CO₂ concentration [ppm] |
| atm_fine_dust_flux | real(:,:) | Fine dust deposition flux from the atmosphere [kg/m²/s] |
| atm_coarse_dust_flux | real(:,:) | Coarse dust deposition flux from the atmosphere [kg/m²/s] |
| seaice_dust_flux | real(:,:) | Dust flux released from sea ice [kg/m²/s] |
| atm_bc_flux | real(:,:) | Black carbon deposition flux from the atmosphere [kg/m²/s] |
| seaice_bc_flux | real(:,:) | Black carbon flux released from sea ice [kg/m²/s] |

---

## Surface wave coupling

NUOPC cap only. These fields support Langmuir turbulence and wave-driven mixing parameterizations in MOM6.

| Field | Type | Description |
|---|---|---|
| lamult | real(:,:) | Langmuir turbulence enhancement factor [dimensionless] |
| stk_wavenumbers | real(:) | Central wavenumber of each Stokes drift band [rad/m]; dimensioned (num_stk_bands) |
| ustkb | real(:,:,:) | Stokes drift spectrum, zonal component, at u-points [m/s]; third dimension indexes wavenumber bands |
| vstkb | real(:,:,:) | Stokes drift spectrum, meridional component, at v-points [m/s]; third dimension indexes wavenumber bands |
| num_stk_bands | integer | Number of Stokes drift wavenumber bands passed through the coupler |

---

## Bookkeeping

| Field | Type | Description |
|---|---|---|
| xtype | integer | Transfer mode for the ice-to-ocean exchange: REGRID (1), REDIST (2), or DIRECT (3) |
| fluxes | type(coupler_2d_bc_type) | Named array of additional per-tracer passive tracer fluxes from ice/atmosphere to ocean |
