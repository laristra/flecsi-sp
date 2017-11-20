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
include_directories(${FleCSI_INCLUDE_DIR})

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
  add_definitions( -DFLECSI_SP_DOUBLE_PRECISION )
  SET (TEST_TOLERANCE 1.0e-14 CACHE STRING "The testing tolerance")
else()
  message(STATUS "Note: Single precision build activated.")
  SET (TEST_TOLERANCE 1.0e-6 CACHE STRING "The testing tolerance")
endif()

add_definitions( -DTEST_TOLERANCE=${TEST_TOLERANCE} )

# size of integer ids to use
option( USE_64BIT_IDS "Type of integer to use for ids" ON )

if( USE_64BIT_IDS ) 
  message(STATUS "Note: using 64 bit integer ids.")
  add_definitions( -DUSE_64BIT_IDS )
else()
  message(STATUS "Note: using 32 bit integer ids.")
endif()

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

cinch_add_library_target(FleCSI-SP flecsi-sp)
list(APPEND FleCSI_SP_LIBRARIES FleCSI-SP)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
