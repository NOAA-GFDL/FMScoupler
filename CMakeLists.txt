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

cmake_minimum_required(VERSION 3.15)
project(fmscoupler LANGUAGES Fortran)

include(FetchContent)

# define some file locations 
set(EXTERNAL_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external)
file(MAKE_DIRECTORY ${EXTERNAL_SOURCE_DIR})

# get FMS source, already has cmake so no extra steps needed
FetchContent_Declare(
  FMS
  GIT_REPOSITORY https://github.com/NOAA-GFDL/FMS.git
  GIT_TAG main
  SOURCE_DIR ${EXTERNAL_SOURCE_DIR}/FMS
)
FetchContent_MakeAvailable(FMS)
message(STATUS "FMS clone complete. FMS source is available at ${FMS_SOURCE_DIR}")

# set variables for component repos 
set(COMPONENTS
  ocean_null
  atmos_null
  land_null
  ice_param
  ice_null
)
set(ocean_null_URL https://github.com/NOAA-GFDL/ocean_null.git)
set(atmos_null_URL https://github.com/NOAA-GFDL/atmos_null.git)
set(land_null_URL https://github.com/NOAA-GFDL/land_null.git)
set(ice_param_URL https://github.com/NOAA-GFDL/ice_param.git)
set(ice_null_URL https://github.com/NOAA-GFDL/ice_null.git)

# clone each component repo
foreach(component_name IN LISTS COMPONENTS)
  set(component_url ${${component_name}_URL})
  FetchContent_Declare(
    ${component_name}
    GIT_REPOSITORY ${component_url}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIR}/${component_name}
    GIT_TAG main
  )
  FetchContent_MakeAvailable(${component_name})
endforeach()

# Build libraries for each external component
foreach(component_name IN LISTS COMPONENTS)
  file(GLOB_RECURSE component_srcs 
    ${EXTERNAL_SOURCE_DIR}/${component_name}/*.F90
  )
  add_library(${component_name} STATIC ${component_srcs})
  target_include_directories(${component_name} PUBLIC ${EXTERNAL_SOURCE_DIR}/FMS)
  target_link_libraries(${component_name} PUBLIC FMS::fms)
  target_compile_definitions(${component_name} PUBLIC INTERNAL_FILE_NML)
  set_target_properties(${component_name} PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
	set_target_properties(${component_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_mods"
	)
endforeach()

# ice_null has additional dependencies instead of just fms
target_link_libraries(ice_null PRIVATE FMS::fms ocean_null ice_param)
target_include_directories(ice_null PRIVATE
	"${CMAKE_CURRENT_BINARY_DIR}/ocean_null_mods"
	"${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods"
)

# set some paths
set(SHARED_MOD_SRCS ${CMAKE_SOURCE_DIR}/shared/surface_flux.F90)
set(SIMPLE_MOD_SRCS
  ${CMAKE_SOURCE_DIR}/simple/flux_exchange.F90
  ${CMAKE_SOURCE_DIR}/simple/ice_model.F90
)
set(FULL_MOD_SRCS
  ${CMAKE_SOURCE_DIR}/full/atm_land_ice_flux_exchange.F90
  ${CMAKE_SOURCE_DIR}/full/atmos_ocean_dep_fluxes_calc.F90
  ${CMAKE_SOURCE_DIR}/full/atmos_ocean_fluxes_calc.F90
  ${CMAKE_SOURCE_DIR}/full/flux_exchange.F90
  ${CMAKE_SOURCE_DIR}/full/full_coupler_mod.F90
  ${CMAKE_SOURCE_DIR}/full/ice_ocean_flux_exchange.F90
  ${CMAKE_SOURCE_DIR}/full/land_ice_flux_exchange.F90
)
set(SIMPLE_EXEC_SRC ${CMAKE_SOURCE_DIR}/simple/coupler_main.F90)
set(FULL_EXEC_SRC ${CMAKE_SOURCE_DIR}/full/coupler_main.F90)

# library for simple coupler modules + shared
add_library(coupler_shared_simple STATIC ${SIMPLE_MOD_SRCS} ${SHARED_MOD_SRCS})
set_target_properties(coupler_shared_simple PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/simple_coupler_mods"
)
set_target_properties(coupler_shared_simple PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_compile_definitions(coupler_shared_simple PUBLIC use_AM3_physics _USE_LEGACY_LAND_)
target_include_directories(coupler_shared_simple PUBLIC
	${CMAKE_CURRENT_BINARY_DIR}/ocean_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
)
target_link_libraries(coupler_shared_simple PRIVATE
  FMS::fms
  ocean_null
  atmos_null
  land_null
  ice_param
)

# library for full coupler modules + shared
add_library(coupler_shared_full STATIC ${FULL_MOD_SRCS} ${SHARED_MOD_SRCS})
set_target_properties(coupler_shared_full PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/full_coupler_mods"
)
set_target_properties(coupler_shared_full PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_compile_definitions(coupler_shared_full PUBLIC use_AM3_physics _USE_LEGACY_LAND_)
target_include_directories(coupler_shared_full PUBLIC
	${CMAKE_SOURCE_DIR}/full
	${CMAKE_CURRENT_BINARY_DIR}/ocean_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_null_mods
)
target_link_libraries(coupler_shared_full PRIVATE
  FMS::fms
  ocean_null
  atmos_null
  land_null
  ice_param
  ice_null
)


# simple coupler executable
add_executable(coupler_simple ${COUPLER_SIMPLE_SRCS} ${CMAKE_SOURCE_DIR}/simple/coupler_main.F90)
set_target_properties(coupler_simple PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_include_directories(coupler_simple PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/simple_coupler_mods
	${CMAKE_CURRENT_BINARY_DIR}/ocean_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
)
target_link_libraries(coupler_simple PRIVATE
  coupler_shared_simple
  FMS::fms
  ocean_null
  atmos_null
  land_null
  ice_param
)
set_target_properties(coupler_simple PROPERTIES RUNTIME_OUTPUT_NAME "coupler_simple.x")

# full coupler executable
add_executable(coupler_full ${COUPLER_FULL_SRCS} ${CMAKE_SOURCE_DIR}/full/coupler_main.F90)
set_target_properties(coupler_full PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_include_directories(coupler_full PUBLIC
	${CMAKE_SOURCE_DIR}/full
  ${CMAKE_CURRENT_BINARY_DIR}/full_coupler_mods
	${CMAKE_CURRENT_BINARY_DIR}/ocean_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_null_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_null_mods
)
target_link_libraries(coupler_full PRIVATE
  coupler_shared_full
  FMS::fms
  ocean_null
  atmos_null
  land_null
  ice_param
  ice_null
)
set_target_properties(coupler_full PROPERTIES RUNTIME_OUTPUT_NAME "coupler_full.x")
