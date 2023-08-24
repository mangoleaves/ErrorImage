#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace GCLF
{
namespace Geometry
{
class Vec2i
{
public:
	typedef Vec2i vector_type;
	typedef int value_type;
  static constexpr const size_t size_ = 2;
	static constexpr const bool GCLF = true;
public:
	int _x;
	int _y;
public:
	Vec2i() {}
	Vec2i(const int& v0, const int& v1) :_x(v0), _y(v1){}
	Vec2i(const int* v) :_x(v[0]), _y(v[1]) {}
	Vec2i(const Vec2i& v) :_x(v._x), _y(v._y) {}

	Vec2i& operator=(const Vec2i& v)
	{
		_x = v._x; _y = v._y;
		return *this;
	}

	~Vec2i() {}

	inline int& x(){ return _x; }
	inline int& y(){ return _y; }
	inline int x()const { return _x; }
	inline int y()const  { return _y; }

	int* data() { return reinterpret_cast<int*>(this); }
	const int* data() const { return reinterpret_cast<const int*>(this); }

	int& operator[](size_t dim) { return reinterpret_cast<int*>(this)[dim]; }
	const int& operator[](size_t dim)const { return reinterpret_cast<const int*>(this)[dim]; }

	inline Vec2i operator-() const { return Vec2i(-_x, -_y); }
	inline Vec2i operator-(const Vec2i& rhs) const { return Vec2i(_x - rhs._x, _y - rhs._y); }
	inline Vec2i operator+(const Vec2i& rhs) const { return Vec2i(_x + rhs._x, _y + rhs._y); }
	inline Vec2i operator*(const int rhs) const { return Vec2i(_x * rhs, _y * rhs); }
	inline Vec2i operator/(const int rhs) const { return Vec2i(_x / rhs, _y / rhs); }
	inline Vec2i& operator-=(const Vec2i& rhs) { _x -= rhs._x; _y -= rhs._y; return *this; }
	inline Vec2i& operator+=(const Vec2i& rhs) { _x += rhs._x; _y += rhs._y; return *this; }
	inline Vec2i& operator*=(const int rhs) { _x *= rhs; _y *= rhs; return *this; }
	inline Vec2i& operator/=(const int rhs) { _x /= rhs; _y /= rhs; return *this; }

	inline bool operator==(const Vec2i& rhs) const { return _x == rhs._x && _y == rhs._y; }
	inline bool operator!=(const Vec2i& rhs) const { return !(*this == rhs); }
	inline bool all_leq(const int& rhs)const { return _x <= rhs && _y <= rhs; }
	inline bool all_leq(const Vec2i& rhs)const { return _x <= rhs._x && _y <= rhs._y; }
	inline bool all_geq(const int& rhs)const { return _x >= rhs && _y >= rhs; }
	inline bool all_geq(const Vec2i& rhs)const { return _x >= rhs._x && _y >= rhs._y; }
	inline bool less_on(size_t dim, const Vec2i& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, int rhs)const { return operator[](dim) < rhs; }

	inline int dot(const Vec2i& rhs) const { return _x * rhs._x + _y * rhs._y; }
	inline int operator|(const Vec2i& rhs)const { return this->dot(rhs); }

	inline int sum() const { return _x + _y; }
	inline int squaredNorm() const { return _x * _x + _y * _y; }
	inline int sqrnorm()const { return squaredNorm(); }
	inline double norm() const { return std::sqrt(squaredNorm()); }
	inline double length()const { return norm(); }

	inline void minimize(const Vec2i& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
	}
	inline void maximize(const Vec2i& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
	}
	inline Vec2i& vectorize(const int& s) { _x = s; _y = s; return *this; }
};

inline Vec2i operator*(const int lhs, const Vec2i& rhs) { return rhs * lhs; }

}// namespace Geometry
}// namespace GCLF