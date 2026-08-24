#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
#*
#* FMS Coupler is free software: you can redistribute it and/or modify
#* it under the terms of the GNU Lesser General Public License as
#* published by the Free Software Foundation, either version 3 of the
#* License, or (at your option) any later version.
#*
#* FMS Coupler is distributed in the hope that it will be useful, but
#* WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#* General Public License for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS Coupler.
#* If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************
# This CMakeLists.txt file is used to build the FMS Coupler, which is the driver for GFDL's coupled models.
# It is designed to be used with the fre-cli build system, which will set the paths to the various components and
# their compile flags. However, it can also be used standalone by setting the paths and flags manually with -D flags
# to cmake.
#
# The "direct" dependencies for the FMS Coupler are:
# - FMS (Flexible Modeling System)
# - Atmosphere Model (e.g. atmos_drivers, usually combined with GFDL_atmos_cubed_sphere)
# - Ocean Model (e.g. MOM6, etc.)
# - Land Model (e.g. LM4, etc.)
# - Ice Model (e.g. SIS2, etc.)
# - Atmosphere Physics Library (e.g. AM4 physics, etc.) - usually required but optional if using atmos_null
#
# These dependencies can be provided by multiple different code repositories, and are usually compiled alongside the 
# libraries they use. The one exception is the atmosphere physics library, which is usually compiled as a separate
# library and linked to the atmosphere model.

cmake_minimum_required(VERSION 3.15)
project(fmscoupler LANGUAGES Fortran)

include(FetchContent)

# define some file locations
set(EXTERNAL_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external)
file(MAKE_DIRECTORY ${EXTERNAL_SOURCE_DIR})

# Option to pick between simple and full coupler builds
set(COUPLER_TYPE "full" CACHE STRING "Determines which coupler executable to build with. Options are 'simple' or 'full'.")

# Set source paths for each component. Model paths are semicolon-separated CMake lists.
# These will typically be set by fre-cli from its compile.yaml file,
# but can also be set manually with -D flags to cmake
set(fms_path                     "" CACHE PATH   "Path to FMS source code")
set(atmos_model_paths            "" CACHE STRING "Paths to Atmosphere Model source code")
set(ocean_model_paths            "" CACHE STRING "Paths to Ocean Model source code")
set(land_model_paths             "" CACHE STRING "Paths to Land Model source code")
set(ice_model_paths              "" CACHE STRING "Paths to Ice Model source code")
set(atmos_physics_library_paths  "" CACHE STRING "Paths to Atmosphere Physics library")

# Macro definitions for each component, these will typically be set by fre-cli from its compile.yaml file,
# but can also be set manually with -D flags to cmake
set(fms_macros                   "" CACHE STRING "Macro definitions for FMS as a semicolon separated list")
set(atmos_model_macros           "" CACHE STRING "Macro definitions for Atmosphere Model as a semicolon separated list")
set(ocean_model_macros           "" CACHE STRING "Macro definitions for Ocean Model as a semicolon separated list")
set(land_macros                  "" CACHE STRING "Macro definitions for Land Model as a semicolon separated list")
set(ice_model_macros             "" CACHE STRING "Macro definitions for Ice Model as a semicolon separated list")
set(atmos_physics_library_macros "" CACHE STRING "Macro definitions for Atmosphere Physics library as a semicolon separated list")

# Other flags set by fre-cli from its compile.yaml file, but can also be set manually with -D flags to cmake
set(fms_flags                   "" CACHE STRING "Compile flags for FMS")
set(atmos_model_flags           "" CACHE STRING "Compile flags for Atmosphere Model")
set(ocean_model_flags           "" CACHE STRING "Compile flags for Ocean Model")
set(land_flags                  "" CACHE STRING "Compile flags for Land Model")
set(ice_model_flags             "" CACHE STRING "Compile flags for Ice Model")
set(atmos_physics_library_flags "" CACHE STRING "Compile flags for Atmosphere Physics library")

# Error checking for required paths.
if(NOT fms_path)
  message(FATAL_ERROR "FMS path not set, use -Dfms_path=PATH")
endif()
if(NOT atmos_model_paths)
  message(FATAL_ERROR "Atmosphere model paths not set, use -Datmos_model_paths=PATH[;PATH...]")
endif()
if(NOT ocean_model_paths)
  message(FATAL_ERROR "Ocean model paths not set, use -Docean_model_paths=PATH[;PATH...]")
endif()
if(NOT land_model_paths)
  message(FATAL_ERROR "Land model paths not set, use -Dland_model_paths=PATH[;PATH...]")
endif()
if(NOT ice_model_paths)
  message(FATAL_ERROR "Ice model paths not set, use -Dice_model_paths=PATH[;PATH...]")
endif()

# fms has its own CMakeLists.txt, so we can just add it as a subdirectory unlike the other components
add_subdirectory(${fms_path} ${CMAKE_CURRENT_BINARY_DIR}/fms_build)
set(FMS_SOURCE_DIR "${fms_path}")

# set variables for component repos, if atmos physics library path is not set, it will be skipped
set(COMPONENTS
  atmos_model
  ocean_model
  land_model
  ice_model
)
if(atmos_physics_library_paths)
  list(APPEND COMPONENTS atmos_physics_library)
endif()

# Check paths and symlink each component repo, this is a bit of a hack to get
# the external components into the coupler source tree.
# TODO look into alternatives
foreach(component_name IN LISTS COMPONENTS)
  file(REMOVE_RECURSE "${EXTERNAL_SOURCE_DIR}/${component_name}")
  foreach(component_path IN LISTS ${component_name}_paths)
    if(NOT EXISTS ${component_path})
      message(FATAL_ERROR "${component_name} path not found: ${component_path}")
    endif()
    get_filename_component(component_dir_name "${component_path}" NAME)
    set(component_link "${EXTERNAL_SOURCE_DIR}/${component_name}/${component_dir_name}")
    if(EXISTS "${component_link}")
      message(FATAL_ERROR "Duplicate ${component_name} source directory name: ${component_dir_name}")
    endif()
    message(STATUS "Creating symlink in ${EXTERNAL_SOURCE_DIR} for ${component_path}")
    file(CREATE_LINK "${component_path}" "${component_link}" SYMBOLIC)
  endforeach()
endforeach()

# Create libraries for each external component
foreach(component_name IN LISTS COMPONENTS)
  file(GLOB_RECURSE CONFIGURE_DEPENDS component_srcs ${EXTERNAL_SOURCE_DIR}/${component_name}/*.F90)
  add_library(${component_name} STATIC ${component_srcs})
  target_include_directories(${component_name} PUBLIC ${EXTERNAL_SOURCE_DIR}/FMS)
  target_link_libraries(${component_name} PUBLIC FMS::fms)
  target_compile_definitions(${component_name} PUBLIC ${${component_name}_macros})
  target_compile_options(${component_name} PUBLIC ${${component_name}_flags})
	set_target_properties(${component_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_mods"
	)
endforeach()


# Select type-specific sources and dependencies. The target definitions below are shared.
set(SHARED_MOD_SRCS ${CMAKE_SOURCE_DIR}/shared/surface_flux.F90)
set(COUPLER_COMMON_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR}/ocean_mods
  ${CMAKE_CURRENT_BINARY_DIR}/atmos_mods
  ${CMAKE_CURRENT_BINARY_DIR}/land_mods
)
set(COUPLER_COMMON_LIBRARIES FMS::fms ocean atmos land)

if(COUPLER_TYPE STREQUAL "simple")
  set(COUPLER_MOD_SRCS
    ${CMAKE_SOURCE_DIR}/simple/flux_exchange.F90
    ${CMAKE_SOURCE_DIR}/simple/ice_model.F90
  )
  set(COUPLER_EXEC_SRC ${CMAKE_SOURCE_DIR}/simple/coupler_main.F90)
  list(APPEND COUPLER_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods)
  list(APPEND COUPLER_SHARED_LIBRARIES ice_param)
  set(COUPLER_NULL_MODEL_MACROS use_AM3_physics)
elseif(COUPLER_TYPE STREQUAL "full")
  set(COUPLER_MOD_SRCS
    ${CMAKE_SOURCE_DIR}/full/atm_land_ice_flux_exchange.F90
    ${CMAKE_SOURCE_DIR}/full/atmos_ocean_dep_fluxes_calc.F90
    ${CMAKE_SOURCE_DIR}/full/atmos_ocean_fluxes_calc.F90
    ${CMAKE_SOURCE_DIR}/full/flux_exchange.F90
    ${CMAKE_SOURCE_DIR}/full/full_coupler_mod.F90
    ${CMAKE_SOURCE_DIR}/full/ice_ocean_flux_exchange.F90
    ${CMAKE_SOURCE_DIR}/full/land_ice_flux_exchange.F90
  )
  set(COUPLER_EXEC_SRC ${CMAKE_SOURCE_DIR}/full/coupler_main.F90)
  list(APPEND COUPLER_INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/full
    ${CMAKE_CURRENT_BINARY_DIR}/ice_mods
    ${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
  )
  list(APPEND COUPLER_SHARED_LIBRARIES ice ice_param)
  list(APPEND COUPLER_EXEC_LIBRARIES ice)
  set(COUPLER_NULL_MODEL_MACROS use_AM3_physics _USE_LEGACY_LAND_)
else()
  message(FATAL_ERROR "Invalid COUPLER_TYPE '${COUPLER_TYPE}'. Expected 'simple' or 'full'.")
endif()

set(COUPLER_SHARED_TARGET coupler_shared_${COUPLER_TYPE})
set(COUPLER_EXECUTABLE_TARGET coupler_${COUPLER_TYPE})
set(COUPLER_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${COUPLER_TYPE}_coupler_mods")

add_library(${COUPLER_SHARED_TARGET} STATIC ${COUPLER_MOD_SRCS} ${SHARED_MOD_SRCS})
set_target_properties(${COUPLER_SHARED_TARGET} PROPERTIES
  COMPILE_FLAGS "-fdefault-real-8"
  Fortran_MODULE_DIRECTORY "${COUPLER_MODULE_DIR}"
)
if(NULL_MODEL_MACROS)
  target_compile_definitions(${COUPLER_SHARED_TARGET} PUBLIC ${COUPLER_NULL_MODEL_MACROS})
endif()
target_include_directories(${COUPLER_SHARED_TARGET} PUBLIC
  ${COUPLER_COMMON_INCLUDE_DIRS}
  ${COUPLER_INCLUDE_DIRS}
)
target_link_libraries(${COUPLER_SHARED_TARGET} PRIVATE
  ${COUPLER_COMMON_LIBRARIES}
  ${COUPLER_SHARED_LIBRARIES}
)

add_executable(${COUPLER_EXECUTABLE_TARGET} ${COUPLER_EXEC_SRC})
set_target_properties(${COUPLER_EXECUTABLE_TARGET} PROPERTIES
  COMPILE_FLAGS "-fdefault-real-8"
  RUNTIME_OUTPUT_NAME "coupler_${COUPLER_TYPE}.x"
)
target_include_directories(${COUPLER_EXECUTABLE_TARGET} PUBLIC
  ${COUPLER_INCLUDE_DIRS}
  ${COUPLER_MODULE_DIR}
  ${COUPLER_COMMON_INCLUDE_DIRS}
)
target_link_libraries(${COUPLER_EXECUTABLE_TARGET} PRIVATE
  ${COUPLER_SHARED_TARGET}
  ${COUPLER_COMMON_LIBRARIES}
  ${COUPLER_EXEC_LIBRARIES}
)

# Test the selected coupler type for null-model builds.
if(TEST_NULL_MODEL)
  enable_testing()
  add_test(NAME full_coupler_null_model_run
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/t/test_full.sh)
  add_test(NAME simple_coupler_null_model_run
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/t/test_simple.sh)
endif()

