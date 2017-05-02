/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include "flecsi/io/io.h"
#include "flecsi-sp/burton/burton_mesh.h"



////////////////////////////////////////////////////////////////////////////////
// General Interface
////////////////////////////////////////////////////////////////////////////////

  
auto filter_boundary = []( auto && entities )  
{
  return 
    std::forward<decltype(entities)>(entities).filter( 
      [](auto e) { return e->is_boundary(); } 
    );
};


////////////////////////////////////////////////////////////////////////////////
// Alias mesh types
////////////////////////////////////////////////////////////////////////////////

namespace flecsi {
namespace sp {
namespace burton {

//! \brief The final 2d mesh type
using burton_mesh_2d_t = burton_mesh_t<2>;
//! \brief The final 3d mesh type
using burton_mesh_3d_t = burton_mesh_t<3>;

}
}
}

//! \brief Expose attributes and attachement sites to all namspaces.
//! This is horrible but it has to be done other wise users need to 
//! write stuff like has_attribute_at( flecsi::sp::burton::attributes::persistent,
//! flecsi::sp::burton::attributes::vertex ).
using namespace flecsi::sp::burton::attributes;

////////////////////////////////////////////////////////////////////////////////
// Delayed includes
////////////////////////////////////////////////////////////////////////////////

#include "burton_io_exodus.h"
#include "burton_io_tecplot.h"
#include "burton_io_vtk.h"
#include "burton_io_vtu.h"
#include "burton_io_vtm.h"

////////////////////////////////////////////////////////////////////////////////
// load some things
////////////////////////////////////////////////////////////////////////////////

namespace flecsi {
namespace sp {
namespace burton {

//! \brief bring write/read mesh into the flecsi::sp::burton namespace
using flecsi::io::write_mesh;
//! \brief bring write/read mesh into the flecsi::sp::burton namespace
using flecsi::io::read_mesh;

}
}
}
