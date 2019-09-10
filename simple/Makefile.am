# This is an automake file for the GFDL FMScoupler project.
# This builds directory ./simple.
#
# Ed Hartnett 2019-09-10

# Build these uninstalled convenience libraries.
#libice_model.la
noinst_LTLIBRARIES =  libflux_exchange.la libcoupler_main.la

# Each convenience library depends on its source.
#libice_model_la_SOURCES = ice_model.F90
libflux_exchange_la_SOURCES = flux_exchange.F90
libcoupler_main_la_SOURCES = coupler_main.F90

# Each mod file depends on the .o file.
#ice_model_mod.mod: ice_model.$(OBJEXT)
flux_exchange_mod.mod: flux_exchange.$(OBJEXT)
coupler_main_mod.mod: coupler_main.$(OBJEXT)

# Some mods are dependant on other mods in this dir.
coupler_main.$(OBJEXT): ice_model_mod.mod flux_exchange_mod.mod

# Mod files are built and then installed as headers.
MODFILES = ice_model_mod.mod flux_exchange_mod.mod coupler_main_mod.mod
include_HEADERS = $(MODFILES)
BUILT_SOURCES = $(MODFILES)

EXTRA_DIST = 

CLEANFILES = *.mod
