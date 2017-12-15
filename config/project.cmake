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

include(cxx14)

check_for_cxx14_compiler(CXX14_COMPILER)

if(CXX14_COMPILER)
    enable_cxx14()
else()
    message(FATAL_ERROR "C++14 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# Load the cinch extras
#------------------------------------------------------------------------------#

cinch_load_extras()

#------------------------------------------------------------------------------#
# FleCSI Library
#------------------------------------------------------------------------------#

find_package(FleCSI CONFIG REQUIRED)
list(APPEND FleCSI_SP_LIBRARIES ${FleCSI_LIBRARIES})
include_directories(${FleCSI_INCLUDE_DIRS})

#------------------------------------------------------------------------------#
# Ristra Library
#------------------------------------------------------------------------------#

find_package(Ristra CONFIG REQUIRED)
list(APPEND FleCSI_SP_LIBRARIES ${Ristra_LIBRARIES})
include_directories(${Ristra_INCLUDE_DIRS})


#------------------------------------------------------------------------------#
# Set directory information
#------------------------------------------------------------------------------#

set(FleCSI_SP_DATA_DIR "${CMAKE_SOURCE_DIR}/data")
set(FleCSI_SP_LIBRARIES)

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
list(APPEND FleCSI_SP_LIBRARIES FleCSI-SP)

#------------------------------------------------------------------------------#
# configure .cmake file (for other projects)
#------------------------------------------------------------------------------#

export(
  TARGETS FleCSI-SP
  FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FleCSI-SPTargets.cmake
)

export(PACKAGE FleCSI-SP)

set(FLECSI_SP_LIBRARY_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR})
set(FLECSI_SP_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(FLECSI_SP_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR}/cmake/FleCSI-SP)

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
