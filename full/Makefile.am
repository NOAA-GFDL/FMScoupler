# This is an automake file for the GFDL AM4 model.
# This builds directory ./full.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
#libland_ice_flux_exchange.la libflux_exchange.la
noinst_LTLIBRARIES = libice_ocean_flux_exchange.la libatm_land_ice_flux_exchange.la
#libcoupler_main.la

# Each convenience library depends on its source.
#libland_ice_flux_exchange_la_SOURCES = land_ice_flux_exchange.F90
#libflux_exchange_la_SOURCES = flux_exchange.F90
libice_ocean_flux_exchange_la_SOURCES = ice_ocean_flux_exchange.F90
libatm_land_ice_flux_exchange_la_SOURCES = atm_land_ice_flux_exchange.F90
#libcoupler_main_la_SOURCES = coupler_main.F90

# Each mod file depends on the .o file.
#land_ice_flux_exchange_mod.mod: land_ice_flux_exchange.$(OBJEXT)
#flux_exchange_mod.mod: flux_exchange.$(OBJEXT)
ice_ocean_flux_exchange_mod.mod: ice_ocean_flux_exchange.$(OBJEXT)
atm_land_ice_flux_exchange_mod.mod: atm_land_ice_flux_exchange.$(OBJEXT)
#coupler_main_mod.mod: coupler_main.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
#flux_exchange.$(OBJEXT): land_ice_flux_exchange_mod.mod ice_ocean_flux_exchange_mod.mod atm_land_ice_flux_exchange_mod.mod
#coupler_main.$(OBJEXT): flux_exchange_mod.mod

# Mod files are built and then installed as headers.
# land_ice_flux_exchange_mod.mod flux_exchange_mod.mod
MODFILES =  ice_ocean_flux_exchange_mod.mod atm_land_ice_flux_exchange_mod.mod
#coupler_main_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = flux_exchange.tech.ps earth_system_model.pdf

CLEANFILES = *.mod
