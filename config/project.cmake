#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#

cinch_minimum_required(VERSION v1.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

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
  set(RISTRA_LIBRARIES Ristra::Ristra CACHE STRING "Ristra library target name (not being set by RistraConfig.cmake")
  list(APPEND FLECSI_SP_LIBRARIES ${RISTRA_LIBRARIES})
endif()


#------------------------------------------------------------------------------#
# Flecsi-sp-specific configuration
#------------------------------------------------------------------------------#

# double or single precision
OPTION (FLECSI_SP_DOUBLE_PRECISION "Use double precision reals"  ON)

if( FLECSI_SP_DOUBLE_PRECISION ) 
  message(STATUS "Note: Double precision build activated.")
  SET (FLECSI_SP_TEST_TOLERANCE 1.0e-13 CACHE STRING "The testing tolerance")
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

find_package(Boost 1.59.0 COMPONENTS program_options QUIET)

# this option overrides what will get set in cinch_load_extras()
option(
  ENABLE_BOOST
  "Enable Boost program options for command-line flags"
  ON
)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND FLECSI_SP_LIBRARIES ${Boost_LIBRARIES} Boost::program_options Boost::boost)

#------------------------------------------------------------------------------#
# Load the cinch extras
#
# Now all options have been set / overriden
#------------------------------------------------------------------------------#

cinch_load_extras()


#------------------------------------------------------------------------------#
# Legion / MPI
#------------------------------------------------------------------------------#
if ( FLECSI_SP_RUNTIME_MODEL STREQUAL "legion" )
  find_package(Legion REQUIRED)

  if (Legion_FOUND)
    include_directories(${Legion_INCLUDE_DIRS})
  endif()
endif()

#TODO: Gate MPI correctly as not all (?) backends will require it
find_package(MPI COMPONENTS C CXX REQUIRED)

if (MPI_FOUND) 
  list(APPEND FLECSI_SP_LIBRARIES MPI::MPI_CXX MPI::MPI_C)
endif()

#------------------------------------------------------------------------------#
# Exodus II
#------------------------------------------------------------------------------#

find_package(EXODUSII QUIET)

option(FLECSI_SP_ENABLE_EXODUS "Enable I/O with exodus." ${EXODUSII_FOUND})

if(FLECSI_SP_ENABLE_EXODUS AND NOT EXODUSII_FOUND)
  message(FATAL_ERROR "Exodus requested, but not found")
endif()

if(FLECSI_SP_ENABLE_EXODUS)
  include_directories(${EXODUSII_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${EXODUSII_LIBRARIES})
endif()


#------------------------------------------------------------------------------#
# Catalyst
#------------------------------------------------------------------------------#

option(FLECSI_SP_ENABLE_CATALYST "Link the sim with Catalyst for in situ" OFF)

if (FLECSI_SP_ENABLE_CATALYST)
  find_package(ParaView REQUIRED)

  if (NOT TARGET ParaView::PythonCatalyst)
    message(FATAL_ERROR
      "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
      "with Catalyst and Python support enabled. Please rebuild ParaView (or "
      "point to a different build of ParaView) with PARAVIEW_USE_PYTHON set "
      "to TRUE")
    #return ()
  endif()

  if (NOT PARAVIEW_USE_MPI)
    message(FATAL_ERROR
      "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
      "with MPI support enabled. Please rebuild ParaView (or point to a "
      "different build of ParaView) with PARAVIEW_USE_MPI set to TRUE")
    #return ()
  endif ()

  message(STATUS "Found Paraview: ${ParaView_DIR}")

  list( APPEND FLECSI_SP_LIBRARIES ParaView::PythonCatalyst VTK::CommonDataModel VTK::ParallelMPI VTK::IOParallelXML)
  add_definitions(-DFLECSI_SP_ENABLE_CATALYST_ON)
  message("Enable catalyst")
endif()


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
  if ( ParMETIS_LIBRARIES AND NOT PARMETIS_LIBRARIES )
    set(PARMETIS_LIBRARIES ${ParMETIS_LIBRARIES})
    set(PARMETIS_INCLUDE_DIRS ${ParMETIS_INCLUDE_DIRS})
  endif()
  include_directories(${PARMETIS_INCLUDE_DIRS})
  list(APPEND FLECSI_SP_LIBRARIES ${PARMETIS_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
# Portage
#------------------------------------------------------------------------------#

find_package(PORTAGE NAMES portage CONFIG QUIET)

option(FLECSI_SP_ENABLE_PORTAGE "Enable Portage Support" ${PORTAGE_FOUND})

if(FLECSI_SP_ENABLE_PORTAGE)
  if(NOT Boost_FOUND)
    message( FATAL_ERROR "Boost is needed for Portage" )
  endif()
  message( STATUS "Portage location: ${PORTAGE_ROOT}" )
  list( APPEND FLECSI_SP_LIBRARIES ${PORTAGE_LIBRARIES} )
endif()

#------------------------------------------------------------------------------#
# HDF5
#------------------------------------------------------------------------------#
find_package(HDF5 REQUIRED)
list( APPEND FLECSI_SP_LIBRARIES ${HDF5_LIBRARIES} )

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

add_subdirectory(apps)

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
