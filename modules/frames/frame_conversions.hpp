#ifndef FRAME_CONVERSIONS_HPP
#define FRAME_CONVERSIONS_HPP

#include "core/state.hpp"

namespace astrodynamics_lib {

// Converts a Cartesian state from ICRF (ECI) to ITRF (ECEF).
// utc_jd: epoch as UTC Julian Date.
CartState icrf_to_itrf(CartState state, double utc_jd);

// Converts a Cartesian state from ITRF (ECEF) to ICRF (ECI).
// utc_jd: epoch as UTC Julian Date.
CartState itrf_to_icrf(CartState state, double utc_jd);

} // namespace astrodynamics_lib

#endif // FRAME_CONVERSIONS_HPP
