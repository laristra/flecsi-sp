#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

add_executable( make_mesh
  main.cc
)
target_link_libraries( make_mesh ${EXODUSII_LIBRARIES} ${Boost_LIBRARIES} Boost::boost Boost::program_options ${HDF5_LIBRARIES} MPI::MPI_CXX)
