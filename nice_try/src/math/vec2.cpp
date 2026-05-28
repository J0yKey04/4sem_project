#include "math/vec2.hpp"

Vec2::Vec2()
    : x(0.0), y(0.0)
{
}

Vec2::Vec2(double x, double y)
    : x(x), y(y)
{
}

double Vec2::length() const {
    return std::sqrt(x * x + y * y);
}

double Vec2::lengthSquared() const {
    return x * x + y * y;
}

Vec2 Vec2::normalized() const {
    double len = length();

    if (len < 1e-12) {
        throw std::runtime_error("Cannot normalize zero vector");
    }

    return Vec2(x / len, y / len);
}

double Vec2::dot(const Vec2& other) const {
    return x * other.x + y * other.y;
}

double Vec2::cross(const Vec2& other) const {
    return x * other.y - y * other.x;
}

Vec2 Vec2::operator+(const Vec2& other) const {
    return Vec2(x + other.x, y + other.y);
}

Vec2 Vec2::operator-(const Vec2& other) const {
    return Vec2(x - other.x, y - other.y);
}

Vec2 Vec2::operator*(double k) const {
    return Vec2(x * k, y * k);
}

Vec2 Vec2::operator/(double k) const {
    if (std::abs(k) < 1e-12) {
        throw std::runtime_error("Division by zero in Vec2");
    }

    return Vec2(x / k, y / k);
}

Vec2& Vec2::operator+=(const Vec2& other) {
    x += other.x;
    y += other.y;
    return *this;
}

Vec2& Vec2::operator-=(const Vec2& other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Vec2& Vec2::operator*=(double k) {
    x *= k;
    y *= k;
    return *this;
}

Vec2& Vec2::operator/=(double k) {
    if (std::abs(k) < 1e-12) {
        throw std::runtime_error("Division by zero in Vec2");
    }

    x /= k;
    y /= k;
    return *this;
}

Vec2 Vec2::operator-() const {
    return Vec2(-x, -y);
}

bool Vec2::isZero(double eps) const {
    return std::abs(x) < eps && std::abs(y) < eps;
}

Vec2 operator*(double k, const Vec2& v) {
    return Vec2(k * v.x, k * v.y);
}

