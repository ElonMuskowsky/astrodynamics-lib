#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

namespace astrodynamics_lib {

constexpr double pi              = 3.14159265358979323846;

constexpr double rad_to_deg      = 180.0 / pi;
constexpr double deg_to_rad      = pi / 180.0;

// Earth
constexpr double mu_earth        = 398600.4415;       // km^3/s^2
constexpr double r_earth         = 6371.0;            // km
constexpr double r_equator_earth = 6378.1363;         // km
constexpr double omega_earth     = 7.292115900231e-5; // rad/s

} // namespace astrodynamics_lib

#endif // CONSTANTS_HPP
