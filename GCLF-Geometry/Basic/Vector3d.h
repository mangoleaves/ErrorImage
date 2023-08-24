#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace GCLF
{
namespace Geometry
{
// TODO: put all VecXY types to a new forlder.
class Vec3d
{
public:
	typedef Vec3d vector_type;
	typedef double value_type;
	static constexpr const size_t size_ = 3;
	static constexpr const bool GCLF = true;
public:
	double _x;
	double _y;
	double _z;
public:
	Vec3d() {}
	Vec3d(const double& v0, const double& v1, const double& v2) :_x(v0), _y(v1), _z(v2) {}
	Vec3d(const double* v) :_x(v[0]), _y(v[1]), _z(v[2]) {}
	Vec3d(const Vec3d& v) :_x(v._x), _y(v._y), _z(v._z) {}
	Vec3d& operator=(const Vec3d& v)
	{
		_x = v._x; _y = v._y; _z = v._z;
		return *this;
	}
	~Vec3d() {}

	inline double& x() { return _x; }
	inline double& y() { return _y; }
	inline double& z() { return _z; }
	inline double x()const { return _x; }
	inline double y()const { return _y; }
	inline double z()const { return _z; }

	inline double* data() { return reinterpret_cast<double*>(this); }
	inline const double* data() const { return reinterpret_cast<const double*>(this); }

	inline double& operator[](size_t dim) { return reinterpret_cast<double*>(this)[dim]; }
	inline const double& operator[](size_t dim)const { return reinterpret_cast<const double*>(this)[dim]; }

	inline Vec3d operator-() const { return Vec3d(-_x, -_y, -_z); }
	inline Vec3d operator-(const Vec3d& rhs) const { return Vec3d(_x - rhs._x, _y - rhs._y, _z - rhs._z); }
	inline Vec3d operator+(const Vec3d& rhs) const { return Vec3d(_x + rhs._x, _y + rhs._y, _z + rhs._z); }
	inline Vec3d operator*(const double rhs) const { return Vec3d(_x * rhs, _y * rhs, _z * rhs); }
	inline Vec3d operator/(const double rhs) const { return Vec3d(_x / rhs, _y / rhs, _z / rhs); }
	inline Vec3d& operator-=(const Vec3d& rhs) { _x -= rhs._x; _y -= rhs._y; _z -= rhs._z; return *this; }
	inline Vec3d& operator+=(const Vec3d& rhs) { _x += rhs._x; _y += rhs._y; _z += rhs._z; return *this; }
	inline Vec3d& operator*=(const double rhs) { _x *= rhs; _y *= rhs; _z *= rhs; return *this; }
	inline Vec3d& operator/=(const double rhs) { _x /= rhs; _y /= rhs; _z /= rhs; return *this; }

	inline bool operator==(const Vec3d& rhs) const { return _x == rhs._x && _y == rhs._y && _z == rhs._z; }
	inline bool operator!=(const Vec3d& rhs) const { return !(*this == rhs); }
	inline bool all_leq(const double& rhs)const { return _x <= rhs && _y <= rhs && _z <= rhs; }
	inline bool all_leq(const Vec3d& rhs)const { return _x <= rhs._x && _y <= rhs._y && _z <= rhs._z; }
	inline bool all_geq(const double& rhs)const { return _x >= rhs && _y >= rhs && _z >= rhs; }
	inline bool all_geq(const Vec3d& rhs)const { return _x >= rhs._x && _y >= rhs._y && _z >= rhs._z; }
	inline bool less_on(size_t dim, const Vec3d& rhs)const { return operator[](dim) < rhs[dim]; }
	inline bool less_on(size_t dim, double rhs)const { return operator[](dim) < rhs; }

	inline double dot(const Vec3d& rhs) const { return _x * rhs._x + _y * rhs._y + _z * rhs._z; }
	inline double operator|(const Vec3d& rhs)const { return this->dot(rhs); }

	inline bool operator<(const Vec3d& rhs) const
	{
		if (_x < rhs._x) return true; else if (_x > rhs._x) return false;
		if (_y < rhs._y) return true; else if (_y > rhs._y) return false;
		if (_z < rhs._z) return true; else return false;
	}
	inline bool operator>(const Vec3d& rhs) const { return rhs < *this; }
	inline bool operator<=(const Vec3d& rhs) const { return *this < rhs || *this == rhs; }
	inline bool operator>=(const Vec3d& rhs) const { return rhs <= *this; }

	inline Vec3d cross(const Vec3d& rhs) const
	{
		return Vec3d(
			_y * rhs._z - _z * rhs._y,
			_z * rhs._x - _x * rhs._z,
			_x * rhs._y - _y * rhs._x
		);
	}
	inline Vec3d operator%(const Vec3d& rhs)const { return this->cross(rhs); }

	inline Vec3d& normalize()
	{
		double sn = squaredNorm();
		return sn > 0 ? (*this) /= sqrt(sn) : (*this);
	}
	inline Vec3d normalized() const
	{
		double sn = squaredNorm();
		return sn > 0 ? (*this) / sqrt(sn) : (*this);
	}

	inline double sum() const { return _x + _y + _z; }
	inline double squaredNorm() const { return _x * _x + _y * _y + _z * _z; }
	inline double sqrnorm()const { return squaredNorm(); }
	inline double norm() const { return std::sqrt(squaredNorm()); }
	inline double length()const { return norm(); }
	inline Vec3d pow(double exp) const
	{
		Vec3d result;
		result._x = std::pow(_x, exp);
		result._y = std::pow(_y, exp);
		result._z = std::pow(_z, exp);
		return result;
	}
	inline double pnorm(double exp) const
	{
		double px = std::pow(std::abs(_x), exp);
		double py = std::pow(std::abs(_y), exp);
		double pz = std::pow(std::abs(_z), exp);
		return std::pow(px + py + pz, 1. / exp);
	}

	inline void minimize(const Vec3d& rhs)
	{
		if (rhs._x < _x)_x = rhs._x;
		if (rhs._y < _y)_y = rhs._y;
		if (rhs._z < _z)_z = rhs._z;
	}
	inline void maximize(const Vec3d& rhs)
	{
		if (rhs._x > _x)_x = rhs._x;
		if (rhs._y > _y)_y = rhs._y;
		if (rhs._z > _z)_z = rhs._z;
	}
	inline Vec3d& vectorize(const double& s) { _x = s; _y = s;_z = s; return *this; }

	inline bool isfinite() const { return std::isfinite(_x) && std::isfinite(_y) && std::isfinite(_z); }
	// TODO: extend these functions to all VecXY types.
	Vec3d floor()
	{
		Vec3d dst;
		dst.x() = std::floor(x());
		dst.y() = std::floor(y());
		dst.z() = std::floor(z());
		return dst;
	}
	Vec3d ceil()
	{
		Vec3d dst;
		dst.x() = std::ceil(x());
		dst.y() = std::ceil(y());
		dst.z() = std::ceil(z());
		return dst;
	}
	Vec3d round()
	{
		Vec3d dst;
		dst.x() = std::round(x());
		dst.y() = std::round(y());
		dst.z() = std::round(z());
		return dst;
	}
};

inline Vec3d operator*(const double lhs, const Vec3d& rhs) { return rhs * lhs; }

}// namespace Geometry
}// namespace GCLF