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

include_directories(${FleCSI_INCLUDE_DIR})

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
