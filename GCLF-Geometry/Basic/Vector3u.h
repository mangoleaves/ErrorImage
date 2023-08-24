#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace GCLF
{
namespace Geometry
{
class Vec3u
{
public:
	typedef Vec3u vector_type;
	typedef uint32_t value_type;
  static constexpr const size_t size_ = 3;
	static constexpr const bool GCLF = true;
public:
	uint32_t _x;
	uint32_t _y;
	uint32_t _z;
public:
	Vec3u() {}
	Vec3u(const uint32_t& v0, const uint32_t& v1, const uint32_t& v2) :_x(v0), _y(v1), _z(v2) {}
	Vec3u(const uint32_t* v) :_x(v[0]), _y(v[1]), _z(v[2]) {}
	Vec3u(const Vec3u& v) :_x(v._x), _y(v._y), _z(v._z) {}
	Vec3u& operator=(const Vec3u& v)
	{
		_x = v._x; _y = v._y; _z = v._z;
		return *this;
	}
	~Vec3u() {}

	inline uint32_t& x(){ return _x; }
	inline uint32_t& y(){ return _y; }
	inline uint32_t& z(){ return _z; }
	inline uint32_t x()const { return _x; }
	inline uint32_t y()const { return _y; }
	inline uint32_t z()const { return _z; }

	inline uint32_t* data() { return reinterpret_cast<uint32_t*>(this); }
	inline const uint32_t* data() const { return reinterpret_cast<const uint32_t*>(this); }

	inline uint32_t& operator[](size_t dim) { return reinterpret_cast<uint32_t*>(this)[dim]; }
	inline const uint32_t& operator[](size_t dim)const { return reinterpret_cast<const uint32_t*>(this)[dim]; }

	inline Vec3u operator-(const Vec3u& rhs) const { return Vec3u(_x - rhs._x, _y - rhs._y, _z - rhs._z); }
	inline Vec3u operator+(const Vec3u& rhs) const { return Vec3u(_x + rhs._x, _y + rhs._y, _z + rhs._z); }
	inline Vec3u operator*(const uint32_t rhs) const { return Vec3u(_x * rhs, _y * rhs, _z * rhs); }
	inline Vec3u operator/(const uint32_t rhs) const { return Vec3u(_x / rhs, _y / rhs, _z / rhs); }
	inline Vec3u& operator-=(const Vec3u& rhs) { _x -= rhs._x; _y -= rhs._y; _z -= rhs._z; return *this; }
	inline Vec3u& operator+=(const Vec3u& rhs) { _x += rhs._x; _y += rhs._y; _z += rhs._z; return *this; }
	inline Vec3u& operator*=(const uint32_t rhs) { _x *= rhs; _y *= rhs; _z *= rhs; return *this; }
	inline Vec3u& operator/=(const uint32_t rhs) { _x /= rhs; _y /= rhs; _z /= rhs; return *this; }

	inline bool operator==(const Vec3u& rhs) const { return _x == rhs._x && _y == rhs._y && _z == rhs._z; }
	inline bool operator!=(const Vec3u& rhs) const { return !(*this == rhs); }
	inline bool all_leq(const uint32_t& rhs)const { return _x <= rhs && _y <= rhs && _z <= rhs; }
	inline bool all_leq(const Vec3u& rhs)const { return _x <= rhs._x && _y <= rhs._y && _z <= rhs._z; }
	inline bool all_geq(const uint32_t& rhs)const { return _x >= rhs && _y >= rhs && _z >= rhs; }
	inline bool all_geq(const Vec3u& rhs)const { return _x >= rhs._x && _y >= rhs._y && _z >= rhs._z; }
	inline bool less_on(size_t dim, const Vec3u& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, uint32_t rhs)const { return operator[](dim) < rhs; }

	inline uint32_t dot(const Vec3u& rhs) const { return _x * rhs._x + _y * rhs._y + _z * rhs._z; }
	inline uint32_t operator|(const Vec3u& rhs)const { return this->dot(rhs); }

	inline uint32_t squaredNorm() const { return _x * _x + _y * _y + _z * _z; }
	inline uint32_t sqrnorm()const { return squaredNorm(); }
	inline double norm() const { return std::sqrt(squaredNorm()); }
	inline double length()const { return norm(); }

	inline void minimize(const Vec3u& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
		if (rhs._z < _z)_z = rhs._z;
	}
	inline void maximize(const Vec3u& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
		if (rhs._z > _z)_z = rhs._z;
	}
	inline Vec3u& vectorize(const uint32_t& s) { _x = s; _y = s;_z = s; return *this; }
};

inline Vec3u operator*(const uint32_t lhs, const Vec3u& rhs) { return rhs * lhs; }

}// namespace Geometry
}// namespace GCLF