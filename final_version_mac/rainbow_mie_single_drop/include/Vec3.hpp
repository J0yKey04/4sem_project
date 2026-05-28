#pragma once

#include <cmath>

struct Vec3 {
    double x;
    double y;
    double z;

    Vec3 operator+(const Vec3& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }

    Vec3 operator-(const Vec3& v) const {
        return {x - v.x, y - v.y, z - v.z};
    }

    Vec3 operator*(double a) const {
        return {x * a, y * a, z * a};
    }

    Vec3 operator/(double a) const {
        return {x / a, y / a, z / a};
    }
};

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

inline Vec3 normalize(const Vec3& v) {
    double n = norm(v);
    if (n < 1e-12) {
        return {0.0, 0.0, 0.0};
    }
    return v / n;
}

inline double clampDouble(double x, double a, double b) {
    if (x < a) return a;
    if (x > b) return b;
    return x;
}


