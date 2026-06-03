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

# Options for local dependencies
set(fms_path 			  "${CMAKE_SOURCE_DIR}/../FMS"        CACHE STRING "Path to local fms repository")
set(atmos_path			"${CMAKE_SOURCE_DIR}/../atmos_null" CACHE STRING "Path to local atmospere component repository")
set(ocean_path			"${CMAKE_SOURCE_DIR}/../ocean_null" CACHE STRING "Path to local ocean component repository")
set(land_path			  "${CMAKE_SOURCE_DIR}/../land_null"  CACHE STRING "Path to local land component repository")
set(ice_path			  "${CMAKE_SOURCE_DIR}/../ice_null"   CACHE STRING "Path to local ice component repository")
set(ice_param_path  "${CMAKE_SOURCE_DIR}/../ice_param"  CACHE STRING "Path to local ice_param component repository")

# Options for cloning depedencies
option(FETCH_FMS				"Clone FMS repository for building" OFF)
option(FETCH_COMPONENTS "Clone model component repositories for building" OFF)
set(fms_tag 			"main" CACHE STRING "Branch/Tag name for FMS. Only used if -DFETCH_FMS=on")
set(atmos_tag 		"main" CACHE STRING "Branch/Tag name for atmosphere component. Only used if -DFETCH_COMPONENTS=on")
set(ocean_tag 		"main" CACHE STRING "Branch/Tag name for ocean component. Only used if -DFETCH_COMPONENTS=on")
set(land_tag 			"main" CACHE STRING "Branch/Tag name for land component. Only used if -DFETCH_COMPONENTS=on")
set(ice_tag 			"main" CACHE STRING "Branch/Tag name for ice component. Only used if -DFETCH_COMPONENTS=on")
set(ice_param_tag	"main" CACHE STRING "Branch/Tag name for ice_param component. Only used if -DFETCH_COMPONENTS=on")
set(fms_repo "https://github.com/noaa-gfdl/fms.git"
	CACHE STRING "Repository URL for FMS. Only used if -DFETCH_FMS=on")
set(atmos_repo "https://github.com/noaa-gfdl/atmos_null.git"
	CACHE STRING "Repository URL for atmosphere component. Only used if -DFETCH_COMPONENTS=on")
set(ocean_repo "https://github.com/noaa-gfdl/ocean_null.git"
	CACHE STRING "Repository URL for ocean component. Only used if -DFETCH_COMPONENTS=on")
set(land_repo "https://github.com/noaa-gfdl/land_null.git"
	CACHE STRING "Repository URL for land component. Only used if -DFETCH_COMPONENTS=on")
set(ice_repo "https://github.com/noaa-gfdl/ice_null.git"
	CACHE STRING "Repository URL for ice component. Only used if -DFETCH_COMPONENTS=on")
set(ice_param_repo "https://github.com/noaa-gfdl/ice_param.git"
	CACHE STRING "Repository URL for ice_param component. Only used if -DFETCH_COMPONENTS=on")

# TODO some more error checking here to ensure fetching/local options make sense
if(NOT FETCH_FMS)
  message(NOTICE "-DFETCH_FMS=off - expecting local FMS repository at ${fms_path}.")
endif()

if(NOT FETCH_COMPONENTS)
  message(NOTICE
    "-DFETCH_COMPONENTS=off - expecting components locally available in ${CMAKE_SOURCE_DIR}/.. if path not specified"
  )
endif()

# build FMS, cloning if -DFETCH_FMS=on, otherwise expects a local path provided
if(FETCH_FMS)
  FetchContent_Declare(
    FMS
    GIT_REPOSITORY ${fms_repo}
    GIT_TAG ${fms_tag}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIR}/FMS
  )
  FetchContent_MakeAvailable(FMS)
  message(STATUS "FMS cloned to ${FMS_SOURCE_DIR}")
else()
  add_subdirectory(${fms_path} ${CMAKE_CURRENT_BINARY_DIR}/fms_build)
endif()

# set variables for component repos
# ice_param isn't a direct dependency of the coupler, but is almost always needed by
# whatever ice_model you end up using so is included here
set(COMPONENTS
  ocean
  atmos
  land
  ice
  ice_param
)
# clone or symlink each component repo
foreach(component_name IN LISTS COMPONENTS)
  if(FETCH_COMPONENTS)
    set(component_url ${${component_name}_repo})
    message(STATUS "Cloning ${component_url} to ${EXTERNAL_SOURCE_DIR}/${component_name}")
    FetchContent_Declare(
      ${component_name}
      GIT_REPOSITORY ${component_url}
      SOURCE_DIR ${EXTERNAL_SOURCE_DIR}/${component_name}
      GIT_TAG main
    )
    FetchContent_MakeAvailable(${component_name})
  else()
    set(component_path ${${component_name}_path} )
    message(STATUS "Creating symlink in ${EXTERNAL_SOURCE_DIR} for ${component_path}")
    file(CREATE_LINK 
      "${component_path}"
      "${EXTERNAL_SOURCE_DIR}/${component_name}" 
      SYMBOLIC
    )
  endif()
endforeach()

# Build libraries for each external component
foreach(component_name IN LISTS COMPONENTS)
  file(GLOB_RECURSE component_srcs ${EXTERNAL_SOURCE_DIR}/${component_name}/*.F90)
  # TODO if/when cmake is added to component libraries, this should use add_subdirectory instead
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
target_link_libraries(ice PRIVATE FMS::fms ocean ice_param)
target_include_directories(ice PRIVATE
	"${CMAKE_CURRENT_BINARY_DIR}/ocean_mods"
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
	${CMAKE_CURRENT_BINARY_DIR}/ocean_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
)
target_link_libraries(coupler_shared_simple PRIVATE
  FMS::fms
  ocean
  atmos
  land
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
	${CMAKE_CURRENT_BINARY_DIR}/ocean_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_param_mods
)
target_link_libraries(coupler_shared_full PRIVATE
  FMS::fms
  ocean
  atmos
  land
  ice
  ice_param
)


# simple coupler executable
add_executable(coupler_simple ${CMAKE_SOURCE_DIR}/simple/coupler_main.F90)
set_target_properties(coupler_simple PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_include_directories(coupler_simple PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/simple_coupler_mods
	${CMAKE_CURRENT_BINARY_DIR}/ocean_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_mods
)
target_link_libraries(coupler_simple PRIVATE
  coupler_shared_simple
  FMS::fms
  ocean
  atmos
  land
)
set_target_properties(coupler_simple PROPERTIES RUNTIME_OUTPUT_NAME "coupler_simple.x")

# full coupler executable
add_executable(coupler_full ${CMAKE_SOURCE_DIR}/full/coupler_main.F90)
set_target_properties(coupler_full PROPERTIES COMPILE_FLAGS "-fdefault-real-8")
target_include_directories(coupler_full PUBLIC
	${CMAKE_SOURCE_DIR}/full
  ${CMAKE_CURRENT_BINARY_DIR}/full_coupler_mods
	${CMAKE_CURRENT_BINARY_DIR}/ocean_mods
	${CMAKE_CURRENT_BINARY_DIR}/atmos_mods
	${CMAKE_CURRENT_BINARY_DIR}/land_mods
	${CMAKE_CURRENT_BINARY_DIR}/ice_mods
)
target_link_libraries(coupler_full PRIVATE
  coupler_shared_full
  FMS::fms
  ocean
  atmos
  land
  ice
)
set_target_properties(coupler_full PROPERTIES RUNTIME_OUTPUT_NAME "coupler_full.x")

# simple test run, only intended for null_model builds 
enable_testing()

add_test(NAME null_model_run
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/t/test_run.sh)
