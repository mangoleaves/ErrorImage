#pragma once

#include "SimpleTypes.h"

namespace GCLF
{
namespace Geometry
{
// TODO: develope an general bounding box class.
class BoundingBox3i
{
protected:
	Vec3i minBound, maxBound;
public:
	BoundingBox3i()
		:minBound(std::numeric_limits<int>::infinity(),
			std::numeric_limits<int>::infinity(),
			std::numeric_limits<int>::infinity()),
		maxBound(-std::numeric_limits<int>::infinity(),
			-std::numeric_limits<int>::infinity(),
			-std::numeric_limits<int>::infinity())
	{}
	BoundingBox3i(const Vec3i& minB, const Vec3i& maxB)
		:minBound(minB), maxBound(maxB)
	{}
	BoundingBox3i(const Vec3i& point)
		:minBound(point), maxBound(point)
	{}

	inline bool operator==(const BoundingBox3i& b)const{ return minBound == b.minBound && maxBound == b.maxBound;}

	inline Vec3i& min() { return minBound; }
	inline Vec3i& max() { return maxBound; }
	inline const Vec3i& min() const { return minBound; }
	inline const Vec3i& max() const { return maxBound; }
	inline int min_coord(size_t dim) { return minBound[dim]; }
	inline int max_coord(size_t dim) { return maxBound[dim]; }
	inline void set_min_coord(size_t dim, int val) { minBound[dim] = val; }
	inline void set_max_coord(size_t dim, int val) { maxBound[dim] = val; }

	BoundingBox3i operator+(const BoundingBox3i& b) const;
	BoundingBox3i& operator+=(const BoundingBox3i& b);
	BoundingBox3i operator+(const Vec3i& p)const;
	BoundingBox3i& operator+=(const Vec3i& p);

	void enlarge(const Vec3i& offset);

	size_t longest_axis() const;

  bool do_intersect(const BoundingBox3i& box)const;
  void intersect(const BoundingBox3i& box);
};

}// namespace Geometry
}// namespace GCLF