#pragma once

#include <cmath>
#include <cstdint>

template <typename T>
struct Vec3 {
  T x, y, z;

  Vec3() : x{0}, y{0}, z{0} {}
  Vec3(const T& x, const T& y, const T& z) : x{x}, y{y}, z{z} {}
  Vec3(const T& value) : x{value}, y{value}, z{value} {}
  Vec3(const Vec3<T>& other) : x{other.x}, y{other.y}, z{other.z} {}
  Vec3& operator=(const Vec3<T>& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }
  Vec3 operator-() const { return {-x, -y, -z}; }
  Vec3 operator-(const Vec3<T>& other) const {
    return {x - other.x, y - other.y, z - other.z};
  }
  Vec3 operator+(const Vec3<T>& other) const {
    return {x + other.x, y + other.y, z + other.z};
  }
  Vec3 operator-(int n) const { return {x - n, y - n, z - n}; }
  Vec3 operator+(int n) const { return {x + n, y + n, z + n}; }

  T dot(const Vec3& other) const {
    return x * other.x + y * other.y + z * other.z;
  }
  Vec3 cross(const Vec3& other) const {
    return {y * other.z - x * other.y, z * other.x - x * other.z,
            x * other.y - y * other.x};
  }

  T length() const { return std::sqrt(x * x + y * y + z * z); }

  Vec3 normalize() const {
    T length = this->length();
    return {x / length, y / length, z / length};
  }

  static Vec3 lerp(const Vec3& a, const Vec3& b, float t) {
    return a + t * (b - a);
  }
};

Vec3<float> operator-(float n, const Vec3<float>& v) {
  return {v.x - n, v.y - n, v.z - n};
}
Vec3<float> operator+(float n, const Vec3<float>& v) {
  return {v.x + n, v.y + n, v.z + n};
}
Vec3<float> operator*(float n, const Vec3<float>& v) {
  return {v.x * n, v.y * n, v.z * n};
}
Vec3<float> operator/(float n, const Vec3<float>& v) {
  return {v.x / n, v.y / n, v.z / n};
}
Vec3<float> operator/(const Vec3<float>& v, float n) {
  return {v.x / n, v.y / n, v.z / n};
}
void operator/=(Vec3<float>& v, float n) {
  v = Vec3<float>{v.x / n, v.y / n, v.z / n};
}
