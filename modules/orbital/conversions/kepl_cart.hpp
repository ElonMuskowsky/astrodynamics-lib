#ifndef KEPL_CART_HPP
#define KEPL_CART_HPP

#include <cmath>
#include <stdexcept>

#include "core/math.hpp"
#include "core/state.hpp"
#include "core/constants.hpp"

namespace astrodynamics_lib {

struct KeplElems {
    double sma;   // semi-major axis [km]
    double ecc;   // eccentricity [-]
    double inc;   // inclination [rad]
    double raan;  // right ascension of the ascending node [rad]
    double aop;   // argument of perigee [rad]
    double ta;    // true anomaly [rad]
};

// Reference: Vallado (2013), Algorithm 10 (COE2RV)
inline CartState coe2rv(const KeplElems& e)
{
    if (std::fabs(e.sma * (1.0 - e.ecc)) < 0.001)
        throw std::runtime_error("coe2rv: periapsis < 0.001 km");
    if (std::fabs(1.0 - e.ecc) < 1e-7)
        throw std::runtime_error("coe2rv: parabolic orbit");
    if (e.ecc > 1.0) {
        double lim = pi - std::acos(1.0 / e.ecc);
        if (std::fabs(e.ta) >= lim)
            throw std::runtime_error("coe2rv: true anomaly outside hyperbolic range");
    }

    double cos_ta     = std::cos(e.ta);
    double inf_factor = 1.0 + e.ecc * cos_ta;
    if (inf_factor < 1e-30)
        throw std::runtime_error("coe2rv: near-infinite orbital radius");

    double p = e.sma * (1.0 - e.ecc * e.ecc);
    if (p < 1e-30)
        throw std::runtime_error("coe2rv: semiparameter near zero");

    double r_mag     = p / inf_factor;
    double aol       = e.aop + e.ta;
    double cos_aol   = std::cos(aol),    sin_aol  = std::sin(aol);
    double cos_raan  = std::cos(e.raan), sin_raan = std::sin(e.raan);
    double cos_inc   = std::cos(e.inc),  sin_inc  = std::sin(e.inc);
    double cos_aop   = std::cos(e.aop),  sin_aop  = std::sin(e.aop);

    Vec3 pos = {
        r_mag * (cos_aol*cos_raan - cos_inc*sin_aol*sin_raan),
        r_mag * (cos_aol*sin_raan + cos_inc*sin_aol*cos_raan),
        r_mag * (sin_aol*sin_inc)
    };

    double factor = std::sqrt(mu_earth / p);
    double vp     = -factor * std::sin(e.ta);
    double vq     =  factor * (e.ecc + cos_ta);

    Vec3 vel = {
        ( cos_raan*cos_aop - sin_raan*sin_aop*cos_inc) * vp + (-cos_raan*sin_aop - sin_raan*cos_aop*cos_inc) * vq,
        ( sin_raan*cos_aop + cos_raan*sin_aop*cos_inc) * vp + (-sin_raan*sin_aop + cos_raan*cos_aop*cos_inc) * vq,
        (sin_aop*sin_inc) * vp + (cos_aop*sin_inc) * vq
    };

    return {pos, vel};
}

// Reference: Vallado (2013), Algorithm 9 (RV2COE)
inline KeplElems rv2coe(const Vec3& pos, const Vec3& vel)
{
    double r_mag = norm(pos);
    double v_mag = norm(vel);

    Vec3   h     = cross(pos, vel);
    double h_mag = norm(h);

    Vec3   n     = cross(Vec3{0.0, 0.0, 1.0}, h);
    double n_mag = norm(n);

    Vec3 e_vec = (pos * (v_mag*v_mag - mu_earth/r_mag) - vel * dot(pos, vel)) / mu_earth;

    KeplElems elems{};
    elems.ecc = norm(e_vec);

    double ksi = v_mag*v_mag / 2.0 - mu_earth / r_mag;
    if (std::fabs(1.0 - elems.ecc) <= 1e-29)
        throw std::runtime_error("rv2coe: parabolic orbit, SMA undefined");

    elems.sma = -mu_earth / (2.0 * ksi);
    elems.inc = std::acos(h.z / h_mag);

    if (elems.inc >= 1e-10) {
        elems.raan = std::acos(n.x / n_mag);
        if (n.y < 0.0) elems.raan = 2.0*pi - elems.raan;
    }

    if (elems.ecc < 1e-10) {
        elems.aop = 0.0;
    } else if (elems.inc < 1e-10) {
        elems.aop = std::acos(e_vec.x / elems.ecc);
        if (e_vec.y < 0.0) elems.aop = 2.0*pi - elems.aop;
    } else {
        elems.aop = std::acos(dot(n, e_vec) / (n_mag * elems.ecc));
        if (e_vec.z < 0.0) elems.aop = 2.0*pi - elems.aop;
    }

    if (elems.ecc >= 1e-10) {
        elems.ta = std::acos(dot(e_vec, pos) / (elems.ecc * r_mag));
        if (dot(pos, vel) < 0.0) elems.ta = 2.0*pi - elems.ta;
    } else if (elems.inc >= 1e-10) {
        elems.ta = std::acos(dot(n, pos) / (n_mag * r_mag));
        if (pos.z < 0.0) elems.ta = 2.0*pi - elems.ta;
    } else {
        elems.ta = std::acos(pos.x / r_mag);
        if (vel.y < 0.0) elems.ta = 2.0*pi - elems.ta;
    }

    return elems;
}

} // namespace astrodynamics_lib

#endif // KEPL_CART_HPP
