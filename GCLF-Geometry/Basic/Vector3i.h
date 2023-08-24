#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace GCLF
{
namespace Geometry
{
class Vec3i
{
public:
	typedef Vec3i vector_type;
	typedef int value_type;
  static constexpr const size_t size_ = 3;
	static constexpr const bool GCLF = true;
public:
	int _x;
	int _y;
	int _z;
public:
	Vec3i() {}
	Vec3i(const int& v0, const int& v1, const int& v2) :_x(v0), _y(v1), _z(v2) {}
	Vec3i(const int* v) :_x(v[0]), _y(v[1]), _z(v[2]) {}
	Vec3i(const Vec3i& v) :_x(v._x), _y(v._y), _z(v._z) {}
	Vec3i& operator=(const Vec3i& v)
	{
		_x = v._x; _y = v._y; _z = v._z;
		return *this;
	}
	~Vec3i() {}

	inline int& x(){ return _x; }
	inline int& y(){ return _y; }
	inline int& z(){ return _z; }
	inline int x()const { return _x; }
	inline int y()const { return _y; }
	inline int z()const { return _z; }

	inline int* data() { return reinterpret_cast<int*>(this); }
	inline const int* data() const { return reinterpret_cast<const int*>(this); }

	inline int& operator[](size_t dim) { return reinterpret_cast<int*>(this)[dim]; }
	inline const int& operator[](size_t dim)const { return reinterpret_cast<const int*>(this)[dim]; }

	inline Vec3i operator-() const { return Vec3i(-_x, -_y, -_z); }
	inline Vec3i operator-(const Vec3i& rhs) const { return Vec3i(_x - rhs._x, _y - rhs._y, _z - rhs._z); }
	inline Vec3i operator+(const Vec3i& rhs) const { return Vec3i(_x + rhs._x, _y + rhs._y, _z + rhs._z); }
	inline Vec3i operator*(const int rhs) const { return Vec3i(_x * rhs, _y * rhs, _z * rhs); }
	inline Vec3i operator/(const int rhs) const { return Vec3i(_x / rhs, _y / rhs, _z / rhs); }
	inline Vec3i& operator-=(const Vec3i& rhs) { _x -= rhs._x; _y -= rhs._y; _z -= rhs._z; return *this; }
	inline Vec3i& operator+=(const Vec3i& rhs) { _x += rhs._x; _y += rhs._y; _z += rhs._z; return *this; }
	inline Vec3i& operator*=(const int rhs) { _x *= rhs; _y *= rhs; _z *= rhs; return *this; }
	inline Vec3i& operator/=(const int rhs) { _x /= rhs; _y /= rhs; _z /= rhs; return *this; }

	inline bool operator==(const Vec3i& rhs) const { return _x == rhs._x && _y == rhs._y && _z == rhs._z; }
	inline bool operator!=(const Vec3i& rhs) const { return !(*this == rhs); }
	inline bool all_leq(const int& rhs)const { return _x <= rhs && _y <= rhs && _z <= rhs; }
	inline bool all_leq(const Vec3i& rhs)const { return _x <= rhs._x && _y <= rhs._y && _z <= rhs._z; }
	inline bool all_geq(const int& rhs)const { return _x >= rhs && _y >= rhs && _z >= rhs; }
	inline bool all_geq(const Vec3i& rhs)const { return _x >= rhs._x && _y >= rhs._y && _z >= rhs._z; }
	inline bool less_on(size_t dim, const Vec3i& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, int rhs)const { return operator[](dim) < rhs; }

	inline int dot(const Vec3i& rhs) const { return _x * rhs._x + _y * rhs._y + _z * rhs._z; }
	inline int operator|(const Vec3i& rhs)const { return this->dot(rhs); }

	inline int squaredNorm() const { return _x * _x + _y * _y + _z * _z; }
	inline int sqrnorm()const { return squaredNorm(); }
	inline double norm() const { return std::sqrt(squaredNorm()); }
	inline double length()const { return norm(); }

	inline void minimize(const Vec3i& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
		if (rhs._z < _z)_z = rhs._z;
	}
	inline void maximize(const Vec3i& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
		if (rhs._z > _z)_z = rhs._z;
	}
	inline Vec3i& vectorize(const int& s) { _x = s; _y = s;_z = s; return *this; }
};

inline Vec3i operator*(const int lhs, const Vec3i& rhs) { return rhs * lhs; }

}// namespace Geometry
}// namespace GCLF