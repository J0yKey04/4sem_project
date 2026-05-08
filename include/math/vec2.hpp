#pragma once

#include <cmath>
#include <stdexcept>

class Vec2 {
public:
    double x;
    double y;

    Vec2();
    Vec2(double x, double y);

    double length() const;
    double lengthSquared() const;

    Vec2 normalized() const;

    double dot(const Vec2& other) const;
    double cross(const Vec2& other) const;

    Vec2 operator+(const Vec2& other) const;
    Vec2 operator-(const Vec2& other) const;

    Vec2 operator*(double k) const;
    Vec2 operator/(double k) const;

    Vec2& operator+=(const Vec2& other);
    Vec2& operator-=(const Vec2& other);
    Vec2& operator*=(double k);
    Vec2& operator/=(double k);

    Vec2 operator-() const;

    bool isZero(double eps = 1e-9) const;
};

Vec2 operator*(double k, const Vec2& v);


