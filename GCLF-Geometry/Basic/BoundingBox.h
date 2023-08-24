#pragma once

#include "Triangle.h"

namespace GCLF
{
namespace Geometry
{
/// <summary>
/// Reference CGAL/Bbox_3.h
/// </summary>
class BoundingBox
{
protected:
	Vec3d minBound, maxBound;
public:
	BoundingBox()
		:minBound(std::numeric_limits<double>::infinity(),
			std::numeric_limits<double>::infinity(),
			std::numeric_limits<double>::infinity()),
		maxBound(-std::numeric_limits<double>::infinity(),
			-std::numeric_limits<double>::infinity(),
			-std::numeric_limits<double>::infinity())
	{}
	BoundingBox(const Vec3d& minB, const Vec3d& maxB)
		:minBound(minB), maxBound(maxB)
	{}
	BoundingBox(const Vec3d& point)
		:minBound(point), maxBound(point)
	{}
	template<template<typename, typename> class STL_Container>
	BoundingBox(const STL_Container<Vec3d, std::allocator<Vec3d>>& points) : BoundingBox()
	{
		for (const auto& p : points)
			*this += p;
	}
	BoundingBox(const Segment& seg)
	{
		minBound = seg.first; minBound.minimize(seg.second);
		maxBound = seg.first; maxBound.maximize(seg.second);
	}
	BoundingBox(const Triangle& tri)
	{
		minBound = tri.ver0;
		maxBound = tri.ver0;
		minBound.minimize(tri.ver1);
		minBound.minimize(tri.ver2);
		maxBound.maximize(tri.ver1);
		maxBound.maximize(tri.ver2);
	}

	inline bool operator==(const BoundingBox& b)const{ return minBound == b.minBound && maxBound == b.maxBound;}

	inline Vec3d& min() { return minBound; }
	inline Vec3d& max() { return maxBound; }
	inline const Vec3d& min() const { return minBound; }
	inline const Vec3d& max() const { return maxBound; }
	inline double min_coord(size_t dim) { return minBound[dim]; }
	inline double max_coord(size_t dim) { return maxBound[dim]; }
	inline void set_min_coord(size_t dim, double val) { minBound[dim] = val; }
	inline void set_max_coord(size_t dim, double val) { maxBound[dim] = val; }

	BoundingBox operator+(const BoundingBox& b) const;
	BoundingBox& operator+=(const BoundingBox& b);
	BoundingBox operator+(const Vec3d& p)const;
	BoundingBox& operator+=(const Vec3d& p);

	void enlarge(const Vec3d& offset);

	size_t longest_axis() const;

	bool do_intersect(const Vec3d& point) const;
	bool do_intersect(const BoundingBox& b) const;
	bool do_intersect(const Segment& segment) const;
	bool do_intersect(const Sphere& sphere) const;
	//bool do_intersect(const ImplicitPlane& plane)const;
	bool do_intersect(const Triangle& tri)const;

	double distance_to_point(const Vec3d& query, Vec3d& dists)const;
private:
	bool do_intersect(const ImplicitPlane& plane, const int swap_dim)const;
};

class DefaultCalcBox
{
public:
	template<typename Primitive>
	BoundingBox operator()(const Primitive& p) const { return BoundingBox(p); }
};

}// namespace Geometry
}// namespace GCLF