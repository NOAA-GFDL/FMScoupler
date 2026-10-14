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
# - Land Model (e.g. LM4, land_lad2, etc.)
# - Ice Model (e.g. SIS2, etc.)
# - Atmosphere Physics Library (e.g. atmos_phys, am5_phys, etc.) - usually required but optional if using atmos_null
#
# These dependencies can be provided by multiple different code repositories, and are usually compiled alongside the
# libraries they use. The one exception is the atmosphere physics library, which is usually compiled as a separate
# library and linked to the atmosphere model.

cmake_minimum_required(VERSION 3.15)
project(fmscoupler LANGUAGES Fortran)

include(FetchContent)

# Options to enable building the full and/or simple coupler executables. Both can be built in the same
# configuration, each with their own set of modules/libraries. The full coupler is built by default.
option(BUILD_FULL_COUPLER "Build the full coupler executable" ON)
option(BUILD_SIMPLE_COUPLER "Build the simple coupler executable" OFF)

if(NOT BUILD_FULL_COUPLER AND NOT BUILD_SIMPLE_COUPLER)
  message(FATAL_ERROR "At least one of BUILD_FULL_COUPLER or BUILD_SIMPLE_COUPLER must be enabled.")
endif()

# Names of the coupler executables, including extensions. Only used for the respective coupler type when enabled
set(FULL_COUPLER_EXECUTABLE_NAME "coupler_full.x" CACHE STRING "Name of the full coupler executable, including extensions")
set(SIMPLE_COUPLER_EXECUTABLE_NAME "coupler_simple.x" CACHE STRING "Name of the simple coupler executable, including extensions")

# this is analagous to using the --use-cpp in mkmf, but will only apply to whatever src is provided for the land_model code
# preprocessed code will end up in the build directory under 'preprocessed_srcs'
option(PREPROCESS_LAND "Uses the GNU cpp command instead of the compiler to preprocess the land model source code. This is needed for lm4p to properly handle strings in macro functions")

set(coupler_macros "" CACHE STRING "CPP macros to use when compiling the source code contained in the FMScoupler repository")

# Set source paths for each component. Model paths should be semicolon-separated lists.
# These will typically be set by fre-cli from its compile.yaml file,
# but can also be set manually with -D flags to cmake
set(fms_paths                    "" CACHE PATH   "Path to FMS source code")
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
set(land_model_macros            "" CACHE STRING "Macro definitions for Land Model as a semicolon separated list")
set(ice_model_macros             "" CACHE STRING "Macro definitions for Ice Model as a semicolon separated list")
set(atmos_physics_library_macros "" CACHE STRING "Macro definitions for Atmosphere Physics library as a semicolon separated list")

# Other flags set by fre-cli from its compile.yaml file, but can also be set manually with -D flags to cmake
set(fms_flags                   "" CACHE STRING "Compile flags for FMS")
set(atmos_model_flags           "" CACHE STRING "Compile flags for Atmosphere Model")
set(ocean_model_flags           "" CACHE STRING "Compile flags for Ocean Model")
set(land_model_flags            "" CACHE STRING "Compile flags for Land Model")
set(ice_model_flags             "" CACHE STRING "Compile flags for Ice Model")
set(atmos_physics_library_flags "" CACHE STRING "Compile flags for Atmosphere Physics library")

# Error checking for required paths.
if(NOT fms_paths)
  message(FATAL_ERROR "FMS path not set, use -Dfms_paths=PATH")
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

# Compiler-specific flags from the mkmf templates are included in cmake/*.cmake. FMS inherits them too.
# This sets the flags using the CMAKE_BUILD_TYPE argument and the detected Fortran compiler
# Build types correspond to mkmf templates, Release/Prod, Repro, and Debug
# CMake builds typically use 'Release' as the default, so Prod is remapped to Release to not mess with cmake's feng shui
# mkmf bonus "targets" are instead handled as options: -DOPENMP, -DOVERRIDE_LIMITS, -DVERBOSE_FLAGS
# if using intel you can also set the ISA instruction set flag via -DINTEL_ISA="..."
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/coupler_compiler_flags.cmake)

# FMS has it's own CMAKE_BUILD_TYPE behavior for Release and Debug, it just sets optimization lvl (and add -g if debug)
# if its set to Prod/Repro for the coupler, we'll reset it to 'Release' for FMS. It'll be overriden by
# whatever optimization the coupler is using.
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Release|Debug)$")
  set(CMAKE_BUILD_TYPE "Release")
endif()
add_subdirectory(${fms_paths} ${CMAKE_CURRENT_BINARY_DIR}/fms_build)
# FMS compiles its sources in OBJECT libraries (fms_f, fms_c, or fms_<kind>_f/_c), which the final fms archive
# only collects, so macros/flags must be applied to those targets or they are silently ignored.
# fms_flags is a space-separated string, so SHELL: is used to split it into separate options.
get_property(_fms_targets DIRECTORY ${fms_paths} PROPERTY BUILDSYSTEM_TARGETS)
list(FILTER _fms_targets INCLUDE REGEX "^fms(_[a-z0-9]+)?(_[cf])?$")
foreach(_fms_tgt IN LISTS _fms_targets)
  if(fms_macros)
    target_compile_definitions(${_fms_tgt} PRIVATE ${fms_macros})
  endif()
  if(fms_flags)
    target_compile_options(${_fms_tgt} PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:${fms_flags}>")
  endif()
endforeach()

unset(CMAKE_BUILD_TYPE)
set(FMS_SOURCE_DIR "${fms_paths}")

# set variables for component repos, if atmos physics library path is not set, it will be skipped
# since its not needed for all builds
set(COMPONENTS
  atmos_model
  ocean_model
  land_model
  ice_model
)
if(atmos_physics_library_paths)
  list(APPEND COMPONENTS atmos_physics_library)
endif()

# Collect source files and include directories from paths for each component
# TODO: add support for embedded CMake builds in every component
foreach(component_name IN LISTS COMPONENTS)
  set(component_srcs "")

  foreach(component_path IN LISTS ${component_name}_paths)
    if(NOT EXISTS "${component_path}")
      message(WARNING "${component_name} path not found: ${component_path}")
      continue()
    endif()
    if(IS_DIRECTORY "${component_path}")
      file(GLOB_RECURSE dir_srcs CONFIGURE_DEPENDS
        "${component_path}/*.F90"
        "${component_path}/*.f90"
        "${component_path}/*.F"
        "${component_path}/*.f"
        "${component_path}/*.c"
        "${component_path}/*.C"
      )
      file(GLOB_RECURSE dir_headers CONFIGURE_DEPENDS
        "${component_path}/*.h"
        "${component_path}/*.H"
        "${component_path}/*.inc"
      )
      foreach(src_file IN LISTS dir_srcs)
        if(EXISTS "${src_file}")
          list(APPEND component_srcs "${src_file}")
        endif()
      endforeach()
      foreach(header_file IN LISTS dir_headers)
        if(EXISTS "${header_file}")
          get_filename_component(header_dir "${header_file}" DIRECTORY)
          list(APPEND component_headers "${header_dir}")
        endif()
      endforeach()

    else()
      if(EXISTS "${component_path}")
        list(APPEND component_srcs "${component_path}")
      endif()
    endif()
  endforeach()

  if(component_srcs)
    set(normalized_component_srcs "")
    foreach(src_file IN LISTS component_srcs)
      get_filename_component(real_src_file "${src_file}" REALPATH)
      list(APPEND normalized_component_srcs "${real_src_file}")
    endforeach()
    set(component_srcs ${normalized_component_srcs})
    list(REMOVE_DUPLICATES component_srcs)
  else()
    message(FATAL_ERROR "No source files found for ${component_name} in ${${component_name}_paths}")
  endif()

  if(component_headers)
    list(REMOVE_DUPLICATES component_headers)
  endif()
  list(APPEND component_headers "${FMS_SOURCE_DIR}/include")
  list(APPEND component_headers "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_mods")

  # compiler-specific flags are set in cmake/coupler_compiler_flags.cmake
  set(compiler_opts "${${component_name}_flags}")
  # cmake will do the same thing as mkmf's --use-cpp option here, preprocess directly with cpp and lable the output with DO_NOT_MODIFY
  # prior to any actual compilation
  if(component_name STREQUAL "land_model" AND PREPROCESS_LAND)
    set(processed_srcs "")
    set(cpp_macros ${${component_name}_macros})
    list(TRANSFORM cpp_macros PREPEND "-D")
    make_directory(${CMAKE_CURRENT_BINARY_DIR}/preprocessed_srcs)
    foreach(src ${component_srcs})
      get_filename_component(fname ${src} NAME_WE)
      # Run standalone cpp on .F90 source files
      set(out_file "${CMAKE_CURRENT_BINARY_DIR}/preprocessed_srcs/${fname}_DO_NOT_MODIFY.f90")
      add_custom_command(
        OUTPUT ${out_file}
        COMMAND cpp -nostdinc -C ${cpp_macros} ${${component_name}_flags} -o "${out_file}" "${src}"
        DEPENDS ${src}
        COMMENT "Preprocessing ${src} with standalone cpp..."
        VERBATIM
      )
      list(APPEND processed_srcs ${out_file})
    endforeach()
    # TODO lm4p has some random test programs that fail to compile, should exclude anything thats not a module from src search
    # for now we'll just exclude files that start with "test_"
    list(FILTER processed_srcs EXCLUDE REGEX "test_.*.f90")
    set(component_srcs ${processed_srcs})
  endif()

  # this writes out our paths for each component, a la list_paths
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_paths.txt" "${component_name} sources: ${component_srcs}")
  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_paths.txt" "\n${component_name} headers: ${component_headers}")

  # Create libraries for each external component
  add_library(${component_name} STATIC ${component_srcs})
  target_include_directories(${component_name} PRIVATE ${component_headers})
  target_link_libraries(${component_name} PUBLIC FMS::fms)
  target_compile_definitions(${component_name} PRIVATE ${${component_name}_macros})
  if(compiler_opts)
    target_compile_options(${component_name} PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:SHELL:${compiler_opts}>")
  endif()
  set_target_properties(${component_name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${component_name}_mods"
  )

endforeach()

# SIS2 (ice_model) uses MOM6 modules, which are only compiled as part of ocean_model.
if(TARGET ice_model AND TARGET ocean_model)
  target_link_libraries(ice_model PUBLIC ocean_model)
endif()

# atmos_model uses modules from the atmosphere physics library, when it's being built.
if(TARGET atmos_model AND TARGET atmos_physics_library)
  target_link_libraries(atmos_model PUBLIC atmos_physics_library)
  target_include_directories(atmos_model PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/atmos_physics_library_mods)
endif()

# Aliases for convenience
if(TARGET ocean_model)
  add_library(ocean ALIAS ocean_model)
endif()
if(TARGET atmos_model)
  add_library(atmos ALIAS atmos_model)
endif()
if(TARGET land_model)
  add_library(land ALIAS land_model)
endif()
if(TARGET ice_model)
  add_library(ice ALIAS ice_model)
endif()

# Select type-specific sources and dependencies. The target definitions below are shared.
set(SHARED_MOD_SRCS ${CMAKE_SOURCE_DIR}/shared/surface_flux.F90)
set(COUPLER_COMMON_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR}/ocean_model_mods
  ${CMAKE_CURRENT_BINARY_DIR}/atmos_model_mods
  ${CMAKE_CURRENT_BINARY_DIR}/land_model_mods
)
set(COUPLER_COMMON_LIBRARIES FMS::fms ocean_model atmos_model land_model)
if(TARGET atmos_physics_library)
  list(APPEND COUPLER_COMMON_LIBRARIES atmos_physics_library)
  list(APPEND COUPLER_COMMON_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/atmos_physics_library_mods)
endif()

set(ENABLED_COUPLER_TYPES "")
if(BUILD_FULL_COUPLER)
  list(APPEND ENABLED_COUPLER_TYPES "full")
endif()
if(BUILD_SIMPLE_COUPLER)
  list(APPEND ENABLED_COUPLER_TYPES "simple")
endif()

# Sets up build for full/simple coupler or both depending on whats requested
foreach(COUPLER_TYPE IN LISTS ENABLED_COUPLER_TYPES)
  # reset per-type variables so they don't carry over between loop iterations
  set(COUPLER_MOD_SRCS "")
  set(COUPLER_EXEC_SRC "")
  set(COUPLER_INCLUDE_DIRS "")
  set(COUPLER_SHARED_LIBRARIES "")
  set(COUPLER_EXEC_LIBRARIES "")

  if(COUPLER_TYPE STREQUAL "simple")
    # The simple coupler provides its own embedded ice model (simple/ice_model.F90, module ice_model_mod),
    # but still needs utility modules (e.g. ice_albedo_mod, ocean_rough_mod) supplied by the ice_model
    # component library. Its own Fortran_MODULE_DIRECTORY is searched before the external ice_model_mods
    # include dir, so the locally-compiled ice_model_mod correctly takes precedence over the external one.
    set(COUPLER_MOD_SRCS
      ${CMAKE_SOURCE_DIR}/simple/flux_exchange.F90
      ${CMAKE_SOURCE_DIR}/simple/ice_model.F90
    )
    set(COUPLER_EXEC_SRC ${CMAKE_SOURCE_DIR}/simple/coupler_main.F90)
    list(APPEND COUPLER_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/ice_model_mods)
    list(APPEND COUPLER_SHARED_LIBRARIES ice_model)
    set(COUPLER_EXECUTABLE_NAME "${SIMPLE_COUPLER_EXECUTABLE_NAME}")
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
      ${CMAKE_CURRENT_BINARY_DIR}/ice_model_mods
    )
    list(APPEND COUPLER_SHARED_LIBRARIES ice_model)
    list(APPEND COUPLER_EXEC_LIBRARIES ice_model)
    set(COUPLER_EXECUTABLE_NAME "${FULL_COUPLER_EXECUTABLE_NAME}")
  else()
    message(FATAL_ERROR "Invalid coupler type '${COUPLER_TYPE}'. Expected 'simple' or 'full'.")
  endif()

  set(COUPLER_SHARED_TARGET coupler_shared_${COUPLER_TYPE})
  set(COUPLER_EXECUTABLE_TARGET coupler_${COUPLER_TYPE})
  set(COUPLER_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${COUPLER_TYPE}_coupler_mods")

  add_library(${COUPLER_SHARED_TARGET} STATIC ${COUPLER_MOD_SRCS} ${SHARED_MOD_SRCS})
  set_target_properties(${COUPLER_SHARED_TARGET} PROPERTIES
    Fortran_MODULE_DIRECTORY "${COUPLER_MODULE_DIR}"
  )
  # The target's own module directory is listed first so that any modules it compiles itself
  # (e.g. simple/ice_model.F90's ice_model_mod) take precedence over any identically-named modules
  # supplied by an external component library (e.g. an ice_model component containing a full ice model).
  target_include_directories(${COUPLER_SHARED_TARGET} PUBLIC
    ${COUPLER_MODULE_DIR}
    ${COUPLER_COMMON_INCLUDE_DIRS}
    ${COUPLER_INCLUDE_DIRS}
  )
  target_link_libraries(${COUPLER_SHARED_TARGET} PRIVATE
    ${COUPLER_COMMON_LIBRARIES}
    ${COUPLER_SHARED_LIBRARIES}
  )

  add_executable(${COUPLER_EXECUTABLE_TARGET} ${COUPLER_EXEC_SRC})
  set_target_properties(${COUPLER_EXECUTABLE_TARGET} PROPERTIES
    RUNTIME_OUTPUT_NAME "${COUPLER_EXECUTABLE_NAME}"
  )
  # COUPLER_MODULE_DIR (holding coupler_shared's own compiled modules) is listed first so that it takes
  # precedence over any identically-named modules in an external component (see note above).
  target_include_directories(${COUPLER_EXECUTABLE_TARGET} PUBLIC
    ${COUPLER_MODULE_DIR}
    ${COUPLER_INCLUDE_DIRS}
    ${COUPLER_COMMON_INCLUDE_DIRS}
  )
  target_link_libraries(${COUPLER_EXECUTABLE_TARGET} PRIVATE
    ${COUPLER_SHARED_TARGET}
    ${COUPLER_COMMON_LIBRARIES}
    ${COUPLER_EXEC_LIBRARIES}
  )

  # add any requested macros for both the executable and coupler modules
  if(coupler_macros)
    target_compile_definitions(${COUPLER_SHARED_TARGET} PUBLIC ${coupler_macros})
    target_compile_definitions(${COUPLER_EXECUTABLE_TARGET} PUBLIC ${coupler_macros})
  endif()
endforeach()

# Testing

option(TEST_NULL_MODEL "Enable a simple runtime test via ctest. Intended for the null component build.")

# Test each enabled coupler type for null-model builds.
if(TEST_NULL_MODEL)
  enable_testing()
  if(BUILD_FULL_COUPLER)
    add_test(NAME full_coupler_null_model_run
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/t/test_full.sh "${FULL_COUPLER_EXECUTABLE_NAME}")
  endif()
  if(BUILD_SIMPLE_COUPLER)
    add_test(NAME simple_coupler_null_model_run
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/t/test_simple.sh "${SIMPLE_COUPLER_EXECUTABLE_NAME}")
  endif()
endif()

