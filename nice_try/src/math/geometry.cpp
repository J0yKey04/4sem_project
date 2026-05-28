#include "math/Geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Geometry {

    bool intersectRayCircle(
        const Ray& ray,
        const Vec2& center,
        double radius,
        double& tHit
    ) {
        if (radius <= 0.0) {
            throw std::runtime_error("Circle radius must be positive");
        }

        const double EPS = 1e-9;

        Vec2 oc = ray.origin - center;

        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;

        double discriminant = b * b - 4.0 * a * c;

        if (discriminant < 0.0) {
            return false;
        }

        double sqrtD = std::sqrt(discriminant);

        double t1 = (-b - sqrtD) / (2.0 * a);
        double t2 = (-b + sqrtD) / (2.0 * a);

        bool t1Valid = t1 > EPS;
        bool t2Valid = t2 > EPS;

        if (t1Valid && t2Valid) {
            tHit = std::min(t1, t2);
            return true;
        }

        if (t1Valid) {
            tHit = t1;
            return true;
        }

        if (t2Valid) {
            tHit = t2;
            return true;
        }

        return false;
    }

    bool intersectRayCircle(
        const Ray& ray,
        const Vec2& center,
        double radius,
        Vec2& hitPoint
    ) {
        double tHit = 0.0;

        if (!intersectRayCircle(ray, center, radius, tHit)) {
            return false;
        }

        hitPoint = ray.pointAt(tHit);
        return true;
    }

    Vec2 circleNormal(
        const Vec2& point,
        const Vec2& center
    ) {
        Vec2 normal = point - center;

        if (normal.isZero()) {
            throw std::runtime_error("Cannot compute circle normal at center");
        }

        return normal.normalized();
    }

    double distance(
        const Vec2& a,
        const Vec2& b
    ) {
        return (a - b).length();
    }

    double angleBetween(
        const Vec2& a,
        const Vec2& b
    ) {
        double lenA = a.length();
        double lenB = b.length();

        if (lenA < 1e-12 || lenB < 1e-12) {
            throw std::runtime_error("Cannot compute angle with zero vector");
        }

        double cosValue = a.dot(b) / (lenA * lenB);

        cosValue = std::clamp(cosValue, -1.0, 1.0);

        return std::acos(cosValue);
    }

}

