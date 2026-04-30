#ifndef MATH_HPP
#define MATH_HPP

#include <cmath>

namespace astrodynamics_lib {

struct Vec3 { double x, y, z; };
struct Mat3 { double m[3][3]; };

// --- Vec3 operators ---

inline Vec3 operator+(Vec3 a, Vec3 b)   { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b)   { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 operator-(Vec3 a)           { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(Vec3 a, double s) { return {a.x*s, a.y*s, a.z*s}; }
inline Vec3 operator*(double s, Vec3 a) { return {a.x*s, a.y*s, a.z*s}; }
inline Vec3 operator/(Vec3 a, double s) { return {a.x/s, a.y/s, a.z/s}; }

// --- Vec3 functions ---

inline double dot(Vec3 a, Vec3 b)   { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline double norm(Vec3 a)          { return std::sqrt(dot(a, a)); }
inline Vec3   cross(Vec3 a, Vec3 b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

// --- Mat3 operators ---

inline Mat3 operator+(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = a.m[i][j] + b.m[i][j];
    return r;
}

inline Mat3 operator-(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = a.m[i][j] - b.m[i][j];
    return r;
}

inline Mat3 operator*(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                r.m[i][j] += a.m[i][k] * b.m[k][j];
    return r;
}

inline Vec3 operator*(const Mat3& m, Vec3 v) {
    return {
        m.m[0][0]*v.x + m.m[0][1]*v.y + m.m[0][2]*v.z,
        m.m[1][0]*v.x + m.m[1][1]*v.y + m.m[1][2]*v.z,
        m.m[2][0]*v.x + m.m[2][1]*v.y + m.m[2][2]*v.z
    };
}

inline Mat3 operator*(const Mat3& a, double s) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = a.m[i][j] * s;
    return r;
}

inline Mat3 operator*(double s, const Mat3& a) { return a * s; }

inline Mat3 operator/(const Mat3& a, double s) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r.m[i][j] = a.m[i][j] / s;
    return r;
}

} // namespace astrodynamics_lib

#endif // MATH_HPP
