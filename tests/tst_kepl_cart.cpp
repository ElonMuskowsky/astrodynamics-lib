#include <QtTest>
#include <cmath>

#include "orbital/conversions/kepl_cart.hpp"

using namespace astrodynamics_lib;

// Floating-point near-equality check with descriptive failure message
#define VERIFY_NEAR(actual, expected, tol) \
    QVERIFY2(std::fabs((actual) - (expected)) < (tol), \
             qPrintable(QString("  actual  : %1\n  expected: %2\n  diff    : %3\n  tol     : %4") \
                        .arg(actual, 0, 'g', 15).arg(expected, 0, 'g', 15) \
                        .arg((actual) - (expected), 0, 'g', 6).arg(tol, 0, 'g', 3)))

static constexpr double POS_TOL = 1e-6;   // km  (~1 mm)
static constexpr double VEL_TOL = 1e-9;   // km/s (~1 µm/s)
static constexpr double ANG_TOL = 1e-10;  // rad
static constexpr double DIM_TOL = 1e-10;  // dimensionless

class TstKeplCart : public QObject
{
    Q_OBJECT

private slots:
    // --- coe2rv ---
    void coe2rv_reference_circular_inclined();
    void coe2rv_reference_sun_sync_elliptic();
    void coe2rv_circular_equatorial_invariants();
    void coe2rv_circular_polar_invariants();
    void coe2rv_highly_elliptic_invariants();

    // --- rv2coe ---
    void rv2coe_reference_circular_inclined();
    void rv2coe_reference_sun_sync_elliptic();
    void rv2coe_circular_equatorial();

    // --- round-trips ---
    void roundtrip_circular_equatorial();
    void roundtrip_circular_inclined();
    void roundtrip_elliptic_inclined();
    void roundtrip_highly_elliptic();
    void roundtrip_retrograde();

    // --- exceptions ---
    void exceptions_parabolic();
    void exceptions_periapsis_too_small();
    void exceptions_hyperbolic_ta_out_of_range();
};

// ============================================================
// coe2rv
// ============================================================

void TstKeplCart::coe2rv_reference_circular_inclined()
{
    // Reference values provided externally.
    // sma=6678.137 km, ecc=0, inc=28.5 deg, raan=aop=ta=0
    KeplElems e{6678.137, 0.0, 28.5 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    VERIFY_NEAR(s.pos.x, 6678.137,              POS_TOL);
    VERIFY_NEAR(s.pos.y, 0.0,                   POS_TOL);
    VERIFY_NEAR(s.pos.z, 0.0,                   POS_TOL);
    VERIFY_NEAR(s.vel.x, 0.0,                   VEL_TOL);
    VERIFY_NEAR(s.vel.y, 6.7895302977176506,    VEL_TOL);
    VERIFY_NEAR(s.vel.z, 3.6864141730136519,    VEL_TOL);
}

void TstKeplCart::coe2rv_reference_sun_sync_elliptic()
{
    // Reference values provided externally.
    // Sun-synchronous near-circular orbit: sma=6872.526 km, ecc=0.00135,
    // inc=97.455 deg, raan=69.1 deg, aop=330 deg, ta=290 deg
    KeplElems e{6872.526, 0.00135,
                97.455 * deg_to_rad, 69.1  * deg_to_rad,
                330.0  * deg_to_rad, 290.0 * deg_to_rad};

    CartState s = coe2rv(e);

    VERIFY_NEAR(s.pos.x,  5955.2903210689792104,  POS_TOL);
    VERIFY_NEAR(s.pos.y, -3422.1231087080955149,  POS_TOL);
    VERIFY_NEAR(s.pos.z,  -106.9868679718580040,  POS_TOL);
    VERIFY_NEAR(s.vel.x,    -0.3989591278997762,  VEL_TOL);
    VERIFY_NEAR(s.vel.y,    -0.9110530741229865,  VEL_TOL);
    VERIFY_NEAR(s.vel.z,     7.5540520409130139,  VEL_TOL);
}

void TstKeplCart::coe2rv_circular_equatorial_invariants()
{
    // For a circular equatorial orbit: |pos|=sma, |vel|=sqrt(mu/sma),
    // pos.z=0, vel.z=0, pos·vel=0.
    KeplElems e{6778.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    double r    = norm(s.pos);
    double v    = norm(s.vel);
    double pdotv = dot(s.pos, s.vel);

    VERIFY_NEAR(r,       6778.0,                         POS_TOL);
    VERIFY_NEAR(v,       std::sqrt(mu_earth / 6778.0),   VEL_TOL);
    VERIFY_NEAR(s.pos.z, 0.0,                            1e-12);
    VERIFY_NEAR(s.vel.z, 0.0,                            1e-12);
    VERIFY_NEAR(pdotv,   0.0,                            1e-9);
}

void TstKeplCart::coe2rv_circular_polar_invariants()
{
    // inc=90 deg, ta=0, raan=aop=0:
    // pos lies on x-axis, vel points in +z direction.
    KeplElems e{7000.0, 0.0, 90.0 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s = coe2rv(e);

    double r = norm(s.pos);
    double v = norm(s.vel);

    VERIFY_NEAR(r,       7000.0,                        POS_TOL);
    VERIFY_NEAR(v,       std::sqrt(mu_earth / 7000.0),  VEL_TOL);
    VERIFY_NEAR(s.pos.x, 7000.0,                        POS_TOL);
    VERIFY_NEAR(s.pos.y, 0.0,                           1e-9);
    VERIFY_NEAR(s.pos.z, 0.0,                           1e-9);
    VERIFY_NEAR(s.vel.x, 0.0,                           VEL_TOL);
    VERIFY_NEAR(s.vel.y, 0.0,                           VEL_TOL);
    VERIFY_NEAR(s.vel.z, std::sqrt(mu_earth / 7000.0),  VEL_TOL);
}

void TstKeplCart::coe2rv_highly_elliptic_invariants()
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
    VERIFY_NEAR(r, r_periapsis, 1e-3);  // 1 m at ta=0

    double energy_actual   = v*v / 2.0 - mu_earth / r;
    double energy_expected = -mu_earth / (2.0 * e.sma);
    VERIFY_NEAR(energy_actual, energy_expected, 1e-9);

    double p          = e.sma * (1.0 - e.ecc * e.ecc);
    double h_expected = std::sqrt(mu_earth * p);
    double h_actual   = norm(cross(s.pos, s.vel));
    VERIFY_NEAR(h_actual, h_expected, 1e-6);
}

// ============================================================
// rv2coe
// ============================================================

void TstKeplCart::rv2coe_reference_circular_inclined()
{
    // Inverse of coe2rv_reference_circular_inclined.
    Vec3 pos{6678.137, 0.0, 0.0};
    Vec3 vel{0.0, 6.7895302977176506, 3.6864141730136519};

    KeplElems e = rv2coe(pos, vel);

    VERIFY_NEAR(e.sma, 6678.137,           POS_TOL);
    VERIFY_NEAR(e.ecc, 0.0,                1e-6);   // small but not exact due to fp precision
    VERIFY_NEAR(e.inc, 28.5 * deg_to_rad,  ANG_TOL);
    VERIFY_NEAR(e.raan, 0.0,               ANG_TOL);
    // aop forced to 0 by rv2coe for circular orbit
    VERIFY_NEAR(e.aop, 0.0,                ANG_TOL);
    VERIFY_NEAR(e.ta,  0.0,                ANG_TOL);
}

void TstKeplCart::rv2coe_reference_sun_sync_elliptic()
{
    // Inverse of coe2rv_reference_sun_sync_elliptic.
    Vec3 pos{ 5955.2903210689792104, -3422.1231087080955149, -106.9868679718580040};
    Vec3 vel{   -0.3989591278997762,    -0.9110530741229865,    7.5540520409130139};

    KeplElems e = rv2coe(pos, vel);

    VERIFY_NEAR(e.sma,  6872.526,            POS_TOL);
    VERIFY_NEAR(e.ecc,  0.00135,             1e-7);
    VERIFY_NEAR(e.inc,  97.455 * deg_to_rad, ANG_TOL);
    VERIFY_NEAR(e.raan,  69.1  * deg_to_rad, ANG_TOL);
    VERIFY_NEAR(e.aop,  330.0  * deg_to_rad, ANG_TOL);
    VERIFY_NEAR(e.ta,   290.0  * deg_to_rad, ANG_TOL);
}

void TstKeplCart::rv2coe_circular_equatorial()
{
    Vec3 pos{6778.0, 0.0, 0.0};
    Vec3 vel{0.0, std::sqrt(mu_earth / 6778.0), 0.0};

    KeplElems e = rv2coe(pos, vel);

    VERIFY_NEAR(e.sma, 6778.0, POS_TOL);
    VERIFY_NEAR(e.ecc, 0.0,    DIM_TOL);
    VERIFY_NEAR(e.inc, 0.0,    ANG_TOL);
}

// ============================================================
// Round-trips: COE -> RV -> COE
// ============================================================

void TstKeplCart::roundtrip_circular_equatorial()
{
    KeplElems orig{6778.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    KeplElems back = rv2coe(coe2rv(orig).pos, coe2rv(orig).vel);

    VERIFY_NEAR(back.sma, orig.sma, POS_TOL);
    VERIFY_NEAR(back.ecc, 0.0,      DIM_TOL);
    VERIFY_NEAR(back.inc, 0.0,      ANG_TOL);
}

void TstKeplCart::roundtrip_circular_inclined()
{
    KeplElems orig{6678.137, 0.0, 28.5 * deg_to_rad, 0.0, 0.0, 0.0};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    VERIFY_NEAR(back.sma, orig.sma, POS_TOL);
    VERIFY_NEAR(back.ecc, 0.0,      DIM_TOL);
    VERIFY_NEAR(back.inc, orig.inc, ANG_TOL);
}

void TstKeplCart::roundtrip_elliptic_inclined()
{
    // Typical orbit: ecc=0.3, inc=51.6 deg (ISS-like inclination, elliptic)
    KeplElems orig{8000.0, 0.3, 51.6 * deg_to_rad,
                   120.0 * deg_to_rad, 45.0 * deg_to_rad, 90.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    VERIFY_NEAR(back.sma,  orig.sma,  POS_TOL);
    VERIFY_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    VERIFY_NEAR(back.inc,  orig.inc,  ANG_TOL);
    VERIFY_NEAR(back.raan, orig.raan, ANG_TOL);
    VERIFY_NEAR(back.aop,  orig.aop,  ANG_TOL);
    VERIFY_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

void TstKeplCart::roundtrip_highly_elliptic()
{
    // ecc=0.95 (highly elliptic, user-provided case). ta not at periapsis.
    KeplElems orig{1.00067e7, 0.95, 28.5 * deg_to_rad,
                   45.0 * deg_to_rad, 30.0 * deg_to_rad, 60.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    VERIFY_NEAR(back.sma,  orig.sma,  1e-3);   // larger scale, relax to 1 m
    VERIFY_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    VERIFY_NEAR(back.inc,  orig.inc,  ANG_TOL);
    VERIFY_NEAR(back.raan, orig.raan, ANG_TOL);
    VERIFY_NEAR(back.aop,  orig.aop,  ANG_TOL);
    VERIFY_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

void TstKeplCart::roundtrip_retrograde()
{
    // Retrograde orbit: inc > 90 deg
    KeplElems orig{7500.0, 0.1, 135.0 * deg_to_rad,
                   60.0 * deg_to_rad, 90.0 * deg_to_rad, 45.0 * deg_to_rad};

    CartState s    = coe2rv(orig);
    KeplElems back = rv2coe(s.pos, s.vel);

    VERIFY_NEAR(back.sma,  orig.sma,  POS_TOL);
    VERIFY_NEAR(back.ecc,  orig.ecc,  DIM_TOL);
    VERIFY_NEAR(back.inc,  orig.inc,  ANG_TOL);
    VERIFY_NEAR(back.raan, orig.raan, ANG_TOL);
    VERIFY_NEAR(back.aop,  orig.aop,  ANG_TOL);
    VERIFY_NEAR(back.ta,   orig.ta,   ANG_TOL);
}

// ============================================================
// Exceptions
// ============================================================

void TstKeplCart::exceptions_parabolic()
{
    KeplElems e{10000.0, 1.0, 0.5, 0.0, 0.0, 0.0};
    QVERIFY_THROWS_EXCEPTION(std::runtime_error, coe2rv(e));
}

void TstKeplCart::exceptions_periapsis_too_small()
{
    // sma*(1-ecc) = 1.0*(1-0.9999) = 0.0001 km < 0.001 km threshold
    KeplElems e{1.0, 0.9999, 0.5, 0.0, 0.0, 0.0};
    QVERIFY_THROWS_EXCEPTION(std::runtime_error, coe2rv(e));
}

void TstKeplCart::exceptions_hyperbolic_ta_out_of_range()
{
    // ecc=1.5: limit_ta = pi - acos(1/1.5) ≈ 2.30 rad
    // ta=pi exceeds that limit
    KeplElems e{-10000.0, 1.5, 0.5, 0.0, 0.0, pi};
    QVERIFY_THROWS_EXCEPTION(std::runtime_error, coe2rv(e));
}

QTEST_APPLESS_MAIN(TstKeplCart)
#include "tst_kepl_cart.moc"
