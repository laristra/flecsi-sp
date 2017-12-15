#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Set the minimum Cinch version
#------------------------------------------------------------------------------#

cinch_minimum_required(1.0)

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
# If a C++14 compiler is available, then set the appropriate flags
#------------------------------------------------------------------------------#

# We need C++ 14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)

#------------------------------------------------------------------------------#
# Load the cinch extras
#------------------------------------------------------------------------------#

cinch_load_extras()

#------------------------------------------------------------------------------#
# FleCSI Library
#------------------------------------------------------------------------------#

find_package(FleCSI CONFIG REQUIRED)
list(APPEND FLECSI_SP_LIBRARIES ${FleCSI_LIBRARIES})
include_directories(${FleCSI_INCLUDE_DIRS})

#------------------------------------------------------------------------------#
# Ristra Library
#------------------------------------------------------------------------------#

find_package(Ristra CONFIG REQUIRED)
list(APPEND FLECSI_SP_LIBRARIES ${RISTRA_LIBRARIES})
include_directories(${RISTRA_INCLUDE_DIRS})


#------------------------------------------------------------------------------#
# Set directory information
#------------------------------------------------------------------------------#

set(FLECSI_SP_DATA_DIR "${CMAKE_SOURCE_DIR}/data")
set(FLECSI_SP_LIBRARIES)

#------------------------------------------------------------------------------#
# Some precision setup
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
# Exodus II
#------------------------------------------------------------------------------#

find_package(EXODUSII QUIET)

option(FLECSI_SP_ENABLE_EXODUS "Enable I/O with exodus." ${EXODUSII_FOUND})

if(FLECSI_SP_ENABLE_EXODUS AND NOT EXODUSII_FOUND)
  message(FATAL_ERROR "Exodus requested, but not found")
endif()

if(FLECSI_SP_ENABLE_EXODUS)
  include_directories(${EXODUSII_INCLUDE_DIRS})
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

cinch_add_library_target(FleCSI-SP flecsi-sp)
list(APPEND FLECSI_SP_LIBRARIES FleCSI-SP)

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
