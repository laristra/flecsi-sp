#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

project(flecsi-sp)


set( FleCSI_SP_DATA_DIR "${CMAKE_SOURCE_DIR}/data" )
set( FleCSI_SP_LIBRARIES )

#------------------------------------------------------------------------------#
# Some precision setup
#------------------------------------------------------------------------------#

# double or single precision
OPTION (DOUBLE_PRECISION "Use double precision reals"  ON)

if( DOUBLE_PRECISION ) 
  message(STATUS "Note: Double precision build activated.")
  add_definitions( -DDOUBLE_PRECISION )
  SET (TEST_TOLERANCE 1.0e-14 CACHE STRING "The testing tolerance" )
else()
  message(STATUS "Note: Single precision build activated.")
  SET (TEST_TOLERANCE 1.0e-6 CACHE STRING "The testing tolerance" )
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

cinch_add_library_target(flecsi-sp specializations)
list( APPEND FleCSI_SP_LIBRARIES flecsi-sp )

#------------------------------------------------------------------------------#
# Set header suffix regular expression
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
