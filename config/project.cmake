#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#

cinch_minimum_required(VERSION v1.0)

#------------------------------------------------------------------------------#
# Set the project name
#------------------------------------------------------------------------------#

project(FleCSI-SP)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

#------------------------------------------------------------------------------#
# Set header suffix regular expression
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

#------------------------------------------------------------------------------#
# If a C++17 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

# We need C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

#------------------------------------------------------------------------------#
# Set some global variables
#------------------------------------------------------------------------------#

# directory information
set(FLECSI_SP_DATA_DIR  "${PROJECT_SOURCE_DIR}/data")
set(FLECSI_SP_TOOLS_DIR "${PROJECT_SOURCE_DIR}/tools")
set(FLECSI_SP_SHARE_DIR "${CMAKE_INSTALL_PREFIX}/share/FleCSI-SP")

# the default test init driver
set(FLECSI_SP_DEFAULT_TEST_INITIALIZATION_DRIVER
  ${FLECSI_SP_TOOLS_DIR}/driver_initialization.cc)

#------------------------------------------------------------------------------#
# FleCSI Library
#
# Ristra libraries come first in case any options depened on what we found.
#------------------------------------------------------------------------------#

file(GLOB _flecsi_contents ${CMAKE_SOURCE_DIR}/flecsi/*)

if ( _flecsi_contents )
  if ( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
    add_subdirectory( ${CMAKE_SOURCE_DIR}/flecsi )
  endif()
  include_directories(${CMAKE_SOURCE_DIR}/flecsi)
  list(APPEND FLECSI_SP_LIBRARIES FleCSI)
  set(FLECSI_SP_RUNTIME_DRIVER
    ${FleCSI_SOURCE_DIR}/flecsi/execution/${FLECSI_RUNTIME_MODEL}/runtime_driver.cc)
else()
  find_package(FleCSI CONFIG REQUIRED)
  include_directories(${FleCSI_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${FleCSI_LIBRARIES})
  set(FLECSI_SP_RUNTIME_DRIVER ${FLECSI_RUNTIME_DRIVER})
endif()

set( FLECSI_SP_RUNTIME_MODEL ${FLECSI_RUNTIME_MODEL} )

if ( FLECSI_SP_RUNTIME_MODEL STREQUAL "mpi" )
  set( ENABLE_MPI ON CACHE BOOL "" FORCE)
  set( FLECSI_SP_UNIT_POLICY MPI )
elseif ( FLECSI_SP_RUNTIME_MODEL STREQUAL "legion" )
  set( FLECSI_SP_UNIT_POLICY LEGION )
else()
  MESSAGE( FATAL_ERROR
    "Unknown FLECSI_SP_RUNTIME_MODEL being used: ${FLECSI_SP_RUNTIME_MODEL}" )
endif()


#------------------------------------------------------------------------------#
# Ristra Library
#------------------------------------------------------------------------------#

file(GLOB _ristra_contents ${CMAKE_SOURCE_DIR}/ristra/*)

if ( _ristra_contents )
  if ( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
    add_subdirectory( ${CMAKE_SOURCE_DIR}/ristra )
  endif()
  include_directories(${RISTRA_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES Ristra)
else()
  find_package(Ristra CONFIG REQUIRED)
  include_directories(${RISTRA_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${RISTRA_LIBRARIES})
endif()


#------------------------------------------------------------------------------#
# Flecsi-sp-specific configuration
#------------------------------------------------------------------------------#

# double or single precision
OPTION (FLECSI_SP_DOUBLE_PRECISION "Use double precision reals"  ON)

if( FLECSI_SP_DOUBLE_PRECISION )
  message(STATUS "Note: Double precision build activated.")
  SET (FLECSI_SP_TEST_TOLERANCE 1.0e-14 CACHE STRING "The testing tolerance")
else()
  message(STATUS "Note: Single precision build activated.")
  SET (FLECSI_SP_TEST_TOLERANCE 1.0e-6 CACHE STRING "The testing tolerance")
endif()

# size of integer ids to use
option( FLECSI_SP_USE_64BIT_IDS "Type of integer to use for ids" ON )

if( FLECSI_SP_USE_64BIT_IDS )
  message(STATUS "Note: using 64 bit integer ids.")
else()
  message(STATUS "Note: using 32 bit integer ids.")
endif()


#------------------------------------------------------------------------------#
# Boost is needed for program options
#------------------------------------------------------------------------------#

find_package(Boost COMPONENTS program_options QUIET)

# this option overrides what will get set in cinch_load_extras()
option(
  ENABLE_BOOST_PROGRAM_OPTIONS
  "Enable Boost program options for command-line flags"
  ${Boost_FOUND}
)

#------------------------------------------------------------------------------#
# Load the cinch extras
#
# Now all options have been set / overriden
#------------------------------------------------------------------------------#

cinch_load_extras()


#------------------------------------------------------------------------------#
# Legion / MPI
#------------------------------------------------------------------------------#

find_package(Legion)

if (Legion_FOUND)
  include_directories(${Legion_INCLUDE_DIRS})
endif()

find_package(MPI)

if (MPI_FOUND)
  set(MPI_LANGUAGE C CACHE STRING "" FORCE)
  include_directories(${MPI_C_INCLUDE_PATH})
endif()


#------------------------------------------------------------------------------#
# Burton Mesh backing filetype
#------------------------------------------------------------------------------#

# Possible options for burton mesh backend
set(FLECSI_SP_BURTON_BACKENDS "EXO;MPAS")

set(FLECSI_SP_BURTON_BACKEND "EXO" CACHE STRING
  "Which underlying mesh file structure to use")

# Set the possible strings to use (only really useful with a GUI/TUI)
# TODO: Should the value be tested? If so, where?
set_property(CACHE FLECSI_SP_BURTON_BACKEND PROPERTY STRINGS
  ${FLECSI_SP_BURTON_BACKENDS})


if(FLECSI_SP_BURTON_BACKEND STREQUAL "EXO")
  find_package(EXODUSII REQUIRED)
  include_directories(${EXODUSII_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${EXODUSII_LIBRARIES})
  add_definitions(-DFLECSI_SP_USE_EXODUS)

elseif(FLECSI_SP_BURTON_BACKEND STREQUAL "MPAS")
  # Note: If the CXX interface is ever needed, add COMPONENTS CXX here.
  find_package(HDF5 REQUIRED)
  include_directories(${HDF5_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${HDF5_LIBRARIES})
  add_definitions(-DFLECSI_SP_USE_MPAS)
else()
  # This will obviously be a problem if there's such a thing as a Burton-less
  # FleCSI-SP build
  message (FATAL_ERROR
    "Unknown or unset FLECSI_SP_BURTON_BACKEND: ${FLECSI_SP_BURTON_BACKEND}")
endif()

# Is this stringification a bad CMAKE practice?
add_definitions("-DFLECSI_SP_USE_${FLECSI_SP_BURTON_BACKEND}")



#------------------------------------------------------------------------------#
# ParMETIS
#------------------------------------------------------------------------------#

# Counter-intuitive variable: set to TRUE to disable test
set(PARMETIS_TEST_RUNS TRUE)
find_package(ParMETIS 4.0)

option(FLECSI_SP_ENABLE_PARMETIS "Enable partitioning with parnetis." ${PARMETIS_FOUND})

if(FLECSI_SP_ENABLE_PARMETIS AND NOT PARMETIS_FOUND)
  message(FATAL_ERROR "Parmetis requested, but not found")
endif()

if(FLECSI_SP_ENABLE_PARMETIS)
  include_directories(${PARMETIS_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${PARMETIS_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
# configure header
#------------------------------------------------------------------------------#

configure_file(${PROJECT_SOURCE_DIR}/config/flecsi-sp-config.h.in
  ${CMAKE_BINARY_DIR}/flecsi-sp-config.h @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/flecsi-sp-config.h DESTINATION include)

include_directories(${CMAKE_BINARY_DIR})

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

cinch_add_library_target(FleCSI-SP flecsi-sp EXPORT_TARGET FleCSI-SPTargets)
target_link_libraries( FleCSI-SP ${FLECSI_SP_LIBRARIES} )

# this has to go here.  Since cinch_add_library_target is a function, it
# cannot propagate anything outside of function scope.
set(FLECSI_SP_BURTON_SPECIALIZATION_INIT
  ${FLECSI_SP_SHARE_DIR}/burton_specialization_init.cc)


#------------------------------------------------------------------------------#
# Extract all project options so they can be exported to the ProjectConfig.cmake
# file.
#------------------------------------------------------------------------------#

get_cmake_property(_variableNames VARIABLES)
string (REGEX MATCHALL "(^|;)FLECSI_SP_[A-Za-z0-9_]*" _matchedVars
  "${_variableNames}")
foreach (_variableName ${_matchedVars})
  set( FLECSI_SP_CONFIG_CODE
    "${FLECSI_SP_CONFIG_CODE}
set(${_variableName} \"${${_variableName}}\")"
  )
endforeach()

#------------------------------------------------------------------------------#
# Prepare variables for ProjectConfig file.
#------------------------------------------------------------------------------#

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY INCLUDE_DIRECTORIES)

foreach(dir ${dirs})
  if(NOT ${dir} MATCHES ${CMAKE_CURRENT_SOURCE_DIR})
    list(APPEND FLECSI_SP_EXTERNAL_INCLUDE_DIRS ${dir})
  endif()
endforeach()

set(FLECSI_SP_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR})
set(FLECSI_SP_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(FLECSI_SP_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSI-SP)

#------------------------------------------------------------------------------#
# Export targets and package.
#------------------------------------------------------------------------------#

export(
  TARGETS FleCSI-SP
  FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSI-SPTargets.cmake
)

export(PACKAGE FleCSI-SP)

#------------------------------------------------------------------------------#
# configure .cmake file (for other projects)
#------------------------------------------------------------------------------#

configure_file(${PROJECT_SOURCE_DIR}/config/FleCSI-SPConfig.cmake.in
  ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSI-SPConfig.cmake @ONLY)

install(
  FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSI-SPConfig.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSI-SP
)

install(
  EXPORT FleCSI-SPTargets
  DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSI-SP
  COMPONENT dev
)


#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
