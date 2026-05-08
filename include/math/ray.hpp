#pragma once

#include "math/Vec2.hpp"

class Ray {
public:
    Vec2 origin;
    Vec2 direction;

    double wavelengthNm;
    double intensity;

    Ray();

    Ray(
        const Vec2& origin,
        const Vec2& direction,
        double wavelengthNm,
        double intensity
    );

    Vec2 pointAt(double t) const;
};


