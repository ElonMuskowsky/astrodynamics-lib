#ifndef STATE_HPP
#define STATE_HPP

#include "core/math.hpp"

namespace astrodynamics_lib {

struct CartState {
    Vec3 pos;  // position [km]
    Vec3 vel;  // velocity [km/s]
};

} // namespace astrodynamics_lib

#endif // STATE_HPP
