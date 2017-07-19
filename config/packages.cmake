#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Check for C++14 compiler.
#------------------------------------------------------------------------------#

include(cxx14)

check_for_cxx14_compiler(CXX14_COMPILER)

if(CXX14_COMPILER)
	enable_cxx14()
else()
	message(FATAL_ERROR "C++14 compatible compiler not found")
endif()

#------------------------------------------------------------------------------#
# FleCSI Library
#------------------------------------------------------------------------------#

find_package(FleCSI REQUIRED)
MESSAGE( STATUS ${FleCSI_LIBRARIES} )
list( APPEND FleCSI_SP_LIBRARIES ${FleCSI_LIBRARIES} )
include_directories(${FleCSI_INCLUDE_DIR})

# This is needed to pick up cinch.h.  If a dependency (library,headers) is 
# a cinch project and is installed somewhere else, and it includes cinch.h,
# which cinch.h should get used.  The one in this project, or the one initially
# used to to make the thirdparty library.
include_directories( ${CMAKE_BINARY_DIR} )

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
