/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef function_h
#define function_h

#include "flecsi/execution/execution.h"

using namespace flecsi;
using namespace flecsi::execution;

double sum(double x, double y) {
  std::cout << "Executing sum" << std::endl;
  return x+y;
} // sum

register_function(sum);

define_function_type(example_function_t, double, double, double);

void driver(int argc, char ** argv) {

  example_function_t example_sum = function_handle(sum);

  auto result = execute_function(example_sum, 2.0, 2.0);

  std::cout << "Result: " << result << std::endl;

} // driver

#endif // function_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
