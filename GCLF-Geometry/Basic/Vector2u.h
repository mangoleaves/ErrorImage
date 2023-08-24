#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace GCLF
{
namespace Geometry
{
class Vec2u
{
public:
	typedef Vec2u vector_type;
	typedef uint32_t value_type;
  static constexpr const size_t size_ = 2;
	static constexpr const bool GCLF = true;
public:
	uint32_t _x;
	uint32_t _y;
public:
	Vec2u() {}
	Vec2u(const uint32_t& v0, const uint32_t& v1) :_x(v0), _y(v1){}
	Vec2u(const uint32_t* v) :_x(v[0]), _y(v[1]) {}
	Vec2u(const Vec2u& v) :_x(v._x), _y(v._y) {}

	Vec2u& operator=(const Vec2u& v)
	{
		_x = v._x; _y = v._y;
		return *this;
	}

	~Vec2u() {}

	inline uint32_t& x(){ return _x; }
	inline uint32_t& y(){ return _y; }
	inline uint32_t x()const { return _x; }
	inline uint32_t y()const  { return _y; }

	uint32_t* data() { return reinterpret_cast<uint32_t*>(this); }
	const uint32_t* data() const { return reinterpret_cast<const uint32_t*>(this); }

	uint32_t& operator[](size_t dim) { return reinterpret_cast<uint32_t*>(this)[dim]; }
	const uint32_t& operator[](size_t dim)const { return reinterpret_cast<const uint32_t*>(this)[dim]; }

	inline Vec2u operator-(const Vec2u& rhs) const { return Vec2u(_x - rhs._x, _y - rhs._y); }
	inline Vec2u operator+(const Vec2u& rhs) const { return Vec2u(_x + rhs._x, _y + rhs._y); }
	inline Vec2u operator*(const uint32_t rhs) const { return Vec2u(_x * rhs, _y * rhs); }
	inline Vec2u operator/(const uint32_t rhs) const { return Vec2u(_x / rhs, _y / rhs); }

	inline bool operator==(const Vec2u& rhs) const { return _x == rhs._x && _y == rhs._y; }
	inline bool operator!=(const Vec2u& rhs) const { return !(*this == rhs); }
	inline bool all_leq(const uint32_t& rhs)const { return _x <= rhs && _y <= rhs; }
	inline bool all_leq(const Vec2d& rhs)const { return _x <= rhs._x && _y <= rhs._y; }
	inline bool all_geq(const uint32_t& rhs)const { return _x >= rhs && _y >= rhs; }
	inline bool all_geq(const Vec2d& rhs)const { return _x >= rhs._x && _y >= rhs._y; }
	inline bool less_on(size_t dim, const Vec2u& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, uint32_t rhs)const { return operator[](dim) < rhs; }

	inline uint32_t dot(const Vec2u& rhs) const { return _x * rhs._x + _y * rhs._y; }

	inline uint32_t sum() const { return _x + _y; }
	inline uint32_t squaredNorm() const { return _x * _x + _y * _y; }
	inline uint32_t sqrnorm()const { return squaredNorm(); }
	inline double norm() const { return std::sqrt(squaredNorm()); }
	inline double length()const { return norm(); }

	inline void minimize(const Vec2u& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
	}
	inline void maximize(const Vec2u& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
	}
	inline Vec2u& vectorize(const uint32_t& s) { _x = s; _y = s; return *this; }
};

inline Vec2u operator*(const uint32_t lhs, const Vec2u& rhs) { return rhs * lhs; }

}// namespace Geometry
}// namespace GCLF