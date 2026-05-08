#include "math/Ray.hpp"

Ray::Ray()
    : origin(0.0, 0.0),
      direction(1.0, 0.0),
      wavelengthNm(550.0),
      intensity(1.0)
{
}

Ray::Ray(
    const Vec2& origin,
    const Vec2& direction,
    double wavelengthNm,
    double intensity
)
    : origin(origin),
      direction(direction.normalized()),
      wavelengthNm(wavelengthNm),
      intensity(intensity)
{
}

Vec2 Ray::pointAt(double t) const {
    return origin + direction * t;
}

