#include <gtest/gtest.h>
#include <cmath>

#include "orbital/conversions/kepl_cart.hpp"

using namespace astrodynamics_lib;

static constexpr double POS_TOL = 1e-6;   // km  (~1 mm)
static constexpr double VEL_TOL = 1e-9;   // km/s (~1 µm/s)
static constexpr double ANG_TOL = 1e-10;  // rad
static constexpr double DIM_TOL = 1e-10;  // dimensionless

// ============================================================
// coe2rv
// ============================================================

TEST(KeplCart, coe2rv_reference_circular_inclined)
{
    // Reference values provided externally.
    // sma=6678.137 km, ecc=0, inc=28.5 deg, raan=aop=ta=0
    KeplElems e{6678.137, 0.0, 28.5 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    EXPECT_NEAR(s.pos.x, 6678.137,              POS_TOL);
    EXPECT_NEAR(s.pos.y, 0.0,                   POS_TOL);
    EXPECT_NEAR(s.pos.z, 0.0,                   POS_TOL);
    EXPECT_NEAR(s.vel.x, 0.0,                   VEL_TOL);
    EXPECT_NEAR(s.vel.y, 6.7895302977176506,    VEL_TOL);
    EXPECT_NEAR(s.vel.z, 3.6864141730136519,    VEL_TOL);
}

TEST(KeplCart, coe2rv_reference_sun_sync_elliptic)
{
    // Reference values provided externally.
    // Sun-synchronous near-circular orbit: sma=6872.526 km, ecc=0.00135,
    // inc=97.455 deg, raan=330 deg, aop=69.1 deg, ta=290 deg
    KeplElems e{6872.526, 0.00135,
                97.455 * deg_to_rad, 330.0 * deg_to_rad,
                69.1   * deg_to_rad, 290.0 * deg_to_rad};

    CartState s = coe2rv(e);

    EXPECT_NEAR(s.pos.x,  5955.2903210689792104,  POS_TOL);
    EXPECT_NEAR(s.pos.y, -3422.1231087080955149,  POS_TOL);
    EXPECT_NEAR(s.pos.z,  -106.9868679718580040,  POS_TOL);
    EXPECT_NEAR(s.vel.x,    -0.3989591278997762,  VEL_TOL);
    EXPECT_NEAR(s.vel.y,    -0.9110530741229865,  VEL_TOL);
    EXPECT_NEAR(s.vel.z,     7.5540520409130139,  VEL_TOL);
}

TEST(KeplCart, coe2rv_circular_equatorial_invariants)
{
    // For a circular equatorial orbit: |pos|=sma, |vel|=sqrt(mu/sma),
    // pos.z=0, vel.z=0, pos·vel=0.
    KeplElems e{6778.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    double r     = norm(s.pos);
    double v     = norm(s.vel);
    double pdotv = dot(s.pos, s.vel);

    EXPECT_NEAR(r,       6778.0,                         POS_TOL);
    EXPECT_NEAR(v,       std::sqrt(mu_earth / 6778.0),   VEL_TOL);
    EXPECT_NEAR(s.pos.z, 0.0,                            1e-12);
    EXPECT_NEAR(s.vel.z, 0.0,                            1e-12);
    EXPECT_NEAR(pdotv,   0.0,                            1e-9);
}

TEST(KeplCart, coe2rv_circular_polar_invariants)
{
    // inc=90 deg, ta=0, raan=aop=0:
    // pos lies on x-axis, vel points in +z direction.
    KeplElems e{7000.0, 0.0, 90.0 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    double r = norm(s.pos);
    double v = norm(s.vel);

    EXPECT_NEAR(r,       7000.0,                        POS_TOL);
    EXPECT_NEAR(v,       std::sqrt(mu_earth / 7000.0),  VEL_TOL);
    EXPECT_NEAR(s.pos.x, 7000.0,                        POS_TOL);
    EXPECT_NEAR(s.pos.y, 0.0,                           1e-9);
    EXPECT_NEAR(s.pos.z, 0.0,                           1e-9);
    EXPECT_NEAR(s.vel.x, 0.0,                           VEL_TOL);
    EXPECT_NEAR(s.vel.y, 0.0,                           VEL_TOL);
    EXPECT_NEAR(s.vel.z, std::sqrt(mu_earth / 7000.0),  VEL_TOL);
}

TEST(KeplCart, coe2rv_highly_elliptic_invariants)
{
    // Note: ecc=0.95 < 1 — this is highly elliptic, not hyperbolic.
    // Verified via three independent invariants:
    //   1. periapsis radius = sma*(1-ecc)
    //   2. specific energy  = -mu/(2*sma)
    //   3. angular momentum = sqrt(mu*p), p = sma*(1-ecc^2)
    KeplElems e{1.00067e7, 0.95, 28.5 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    double r = norm(s.pos);
    double v = norm(s.vel);

    double r_periapsis = e.sma * (1.0 - e.ecc);
    EXPECT_NEAR(r, r_periapsis, 1e-3);  // 1 m at ta=0

    double energy_actual   = v*v / 2.0 - mu_earth / r;
    double energy_expected = -mu_earth / (2.0 * e.sma);
    EXPECT_NEAR(energy_actual, energy_expected, 1e-9);

    double p          = e.sma * (1.0 - e.ecc * e.ecc);
    double h_expected = std::sqrt(mu_earth * p);
    double h_actual   = norm(cross(s.pos, s.vel));
    EXPECT_NEAR(h_actual, h_expected, 1e-6);
}

// ============================================================
// rv2coe
// ============================================================

TEST(KeplCart, rv2coe_reference_circular_inclined)
{
    // Inverse of coe2rv_reference_circular_inclined.
    Vec3 pos{6678.137, 0.0, 0.0};
    Vec3 vel{0.0, 6.7895302977176506, 3.6864141730136519};

    KeplElems e = rv2coe(pos, vel);

    EXPECT_NEAR(e.sma,  6678.137,           POS_TOL);
    EXPECT_NEAR(e.ecc,  0.0,                1e-6);
    EXPECT_NEAR(e.inc,  28.5 * deg_to_rad,  ANG_TOL);
    EXPECT_NEAR(e.raan, 0.0,                ANG_TOL);
    EXPECT_NEAR(e.aop,  0.0,                ANG_TOL);
    EXPECT_NEAR(e.ta,   0.0,                ANG_TOL);
}

TEST(KeplCart, rv2coe_reference_sun_sync_elliptic)
{
    // Inverse of coe2rv_reference_sun_sync_elliptic.
    Vec3 pos{ 5955.2903210689792104, -3422.1231087080955149, -106.9868679718580040};
    Vec3 vel{   -0.3989591278997762,    -0.9110530741229865,    7.5540520409130139};

    KeplElems e = rv2coe(pos, vel);

    EXPECT_NEAR(e.sma,  6872.526,             POS_TOL);
    EXPECT_NEAR(e.ecc,  0.00135,             1e-7);
    EXPECT_NEAR(e.inc,  97.455 * deg_to_rad, ANG_TOL);
    EXPECT_NEAR(e.raan, 330.0  * deg_to_rad, ANG_TOL);
    EXPECT_NEAR(e.aop,   69.1  * deg_to_rad, ANG_TOL);
    EXPECT_NEAR(e.ta,   290.0  * deg_to_rad, ANG_TOL);
}

TEST(KeplCart, rv2coe_circular_equatorial)
{
    Vec3 pos{6778.0, 0.0, 0.0};
    Vec3 vel{0.0, std::sqrt(mu_earth / 6778.0), 0.0};

    KeplElems e = rv2coe(pos, vel);

    EXPECT_NEAR(e.sma, 6778.0, POS_TOL);
    EXPECT_NEAR(e.ecc, 0.0,    DIM_TOL);
    EXPECT_NEAR(e.inc, 0.0,    ANG_TOL);
}

// ============================================================
// Round-trips: COE -> RV -> COE
// ============================================================

TEST(KeplCart, roundtrip_circular_equatorial)
{
    KeplElems orig{6778.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    KeplElems back = rv2coe(coe2rv(orig).pos, coe2rv(orig).vel);

    EXPECT_NEAR(back.sma, orig.sma, POS_TOL);
    EXPECT_NEAR(back.ecc, 0.0,      DIM_TOL);
    EXPECT_NEAR(back.inc, 0.0,      ANG_TOL);
}

TEST(KeplCart, roundtrip_circular_inclined)
{
    KeplElems orig{6678.137, 0.0, 28.5 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    EXPECT_NEAR(back.sma, orig.sma, POS_TOL);
    EXPECT_NEAR(back.ecc, 0.0,      DIM_TOL);
    EXPECT_NEAR(back.inc, orig.inc, ANG_TOL);
}

TEST(KeplCart, roundtrip_elliptic_inclined)
{
    // Typical orbit: ecc=0.3, inc=51.6 deg (ISS-like inclination, elliptic)
    KeplElems orig{8000.0, 0.3, 51.6 * deg_to_rad,
                   120.0 * deg_to_rad, 45.0 * deg_to_rad, 90.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    EXPECT_NEAR(back.sma,  orig.sma,  POS_TOL);
    EXPECT_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    EXPECT_NEAR(back.inc,  orig.inc,  ANG_TOL);
    EXPECT_NEAR(back.raan, orig.raan, ANG_TOL);
    EXPECT_NEAR(back.aop,  orig.aop,  ANG_TOL);
    EXPECT_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

TEST(KeplCart, roundtrip_highly_elliptic)
{
    // ecc=0.95 (highly elliptic). ta not at periapsis.
    KeplElems orig{1.00067e7, 0.95, 28.5 * deg_to_rad,
                   45.0 * deg_to_rad, 30.0 * deg_to_rad, 60.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    EXPECT_NEAR(back.sma,  orig.sma,  1e-3);
    EXPECT_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    EXPECT_NEAR(back.inc,  orig.inc,  ANG_TOL);
    EXPECT_NEAR(back.raan, orig.raan, ANG_TOL);
    EXPECT_NEAR(back.aop,  orig.aop,  ANG_TOL);
    EXPECT_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

TEST(KeplCart, roundtrip_retrograde)
{
    // Retrograde orbit: inc > 90 deg
    KeplElems orig{7500.0, 0.1, 135.0 * deg_to_rad,
                   60.0 * deg_to_rad, 90.0 * deg_to_rad, 45.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    EXPECT_NEAR(back.sma,  orig.sma,  POS_TOL);
    EXPECT_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    EXPECT_NEAR(back.inc,  orig.inc,  ANG_TOL);
    EXPECT_NEAR(back.raan, orig.raan, ANG_TOL);
    EXPECT_NEAR(back.aop,  orig.aop,  ANG_TOL);
    EXPECT_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

// ============================================================
// Exceptions
// ============================================================

TEST(KeplCart, exceptions_parabolic)
{
    KeplElems e{10000.0, 1.0, 0.5, 0.0, 0.0, 0.0};
    EXPECT_THROW(coe2rv(e), std::runtime_error);
}

TEST(KeplCart, exceptions_periapsis_too_small)
{
    // sma*(1-ecc) = 1.0*(1-0.9999) = 0.0001 km < 0.001 km threshold
    KeplElems e{1.0, 0.9999, 0.5, 0.0, 0.0, 0.0};
    EXPECT_THROW(coe2rv(e), std::runtime_error);
}

TEST(KeplCart, exceptions_hyperbolic_ta_out_of_range)
{
    // ecc=1.5: limit_ta = pi - acos(1/1.5) ≈ 2.30 rad
    // ta=pi exceeds that limit
    KeplElems e{-10000.0, 1.5, 0.5, 0.0, 0.0, pi};
    EXPECT_THROW(coe2rv(e), std::runtime_error);
}
