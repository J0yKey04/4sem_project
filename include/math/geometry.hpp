#pragma once

#include "math/Vec2.hpp"
#include "math/Ray.hpp"

namespace Geometry {

    bool intersectRayCircle(
        const Ray& ray,
        const Vec2& center,
        double radius,
        double& tHit
    );

    bool intersectRayCircle(
        const Ray& ray,
        const Vec2& center,
        double radius,
        Vec2& hitPoint
    );

    Vec2 circleNormal(
        const Vec2& point,
        const Vec2& center
    );

    double distance(
        const Vec2& a,
        const Vec2& b
    );

    double angleBetween(
        const Vec2& a,
        const Vec2& b
    );

}

