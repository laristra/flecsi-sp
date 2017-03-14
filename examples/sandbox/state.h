#ifndef state_h
#define state_h

#include <array>

template<size_t D>
using vector_t = std::array<double, D>;

template<size_t M, size_t N>
using tensor__ = std::array<std::array<double, N>, M>;
using tensor_t = tensor__<3,3>;	

//----------------------------------------------------------------------------//
// Eos functions.
//----------------------------------------------------------------------------//

// FIXME: Add real physics.
double gruneisen(double r, double e) {
	return r+e+1.8;
} // gruneisen

flecsi_register_function(gruneisen);

// FIXME: Add real physics.
double sesame(double r, double e) {
	return r+e+27.0;
} // sesame

flecsi_register_function(sesame);

//----------------------------------------------------------------------------//
// Stress functions.
//----------------------------------------------------------------------------//

// FIXME: Add real physics.
double ptw(double var1, double var2) {
	return var1+var2+3.1415;
} // ptw

flecsi_register_function(ptw);

// FIXME: Add real physics.
double kospall(double var1, double var2) {
	return var1+var2+2.8;
} // kospall

flecsi_register_function(kospall);

//----------------------------------------------------------------------------//
// Data structures.
//----------------------------------------------------------------------------//

flecsi_define_function_type(eos_function_t, double, double, double);
flecsi_define_function_type(stress_function_t, double, double, double);

///
/// Material state structure.
///
template<size_t D>
struct mat_state__ {
	double volume_fraction;
	double density;
	double energy;
	double pressure;
	double sound_speed;

	eos_function_t eos_function;

	double eos(double r, double e) {
		return flecsi_execute_function(eos_function, r, e);
	} // eos

	stress_function_t stress_function;

	double stress(double var1, double var2) {
		return flecsi_execute_function(stress_function, var1, var2);
	} // stress

	vector_t<D> velocity;
}; // struct mat_state__

///
/// Cell state structure.
///
template<size_t D>
struct cell_state__ {
	double mass;
	double density;
	double energy;
	double pressure;
	double sound_speed;

	vector_t<D> velocity;
}; // struct cell_state__

#if 0
struct mat_model_t {
	int eos_type;
	int strength_type;
	int 
}; // struct mat_model
#endif

#ifndef DIMENSION
	using mat_state_t = mat_state__<3>;
	using cell_state_t = cell_state__<3>;
#else
	using mat_state_t = mat_state__<DIMENSION>;
	using cell_state_t = cell_state__<DIMENSION>;
#endif

#endif // state_h
