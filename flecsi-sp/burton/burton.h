/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes
#include <flecsi-sp/burton/burton_mesh.h>

////////////////////////////////////////////////////////////////////////////////
// General Interface
////////////////////////////////////////////////////////////////////////////////
 
template< typename E >
auto filter_boundary( E && entities )  
{
  return 
    std::forward<decltype(entities)>(entities).filter( 
      [](auto e) { return e->is_boundary(); } 
    );
};


