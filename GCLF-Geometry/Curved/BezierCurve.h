#pragma once
#include "BezierCommon.h"

namespace GCLF
{
namespace Geometry
{
using namespace Utils;

class BezierCurve
{
public:

  /* Operations */

  static Vec3d deCasteljau(const Points& ctrl_points, const uint32_t degree, const double t);

  static Vec3d eval_by_Berstain(const Points& ctrl_points, const uint32_t degree, const double t);

  static Points degree_elevation(const Points& ctrl_points, const uint32_t degree) { ASSERT(false, "no implement"); return Points(); }

  static void subdivision(
    const Points& ctrl_points, const uint32_t degree, const double t,
    std::array<Points, 2>& sub_ctrl_points, Vec3d& sub_center_ctrl_pnt)
  {
    ASSERT(false, "no implement");
  }

  static Points gradient_vectors(const Points& ctrl_points, const uint32_t degree) { ASSERT(false, "no implement"); return Points(); }

  /* Visualization */

  static Points linear_approx(const Points& ctrl_points, const uint32_t degree, uint32_t density);

private:
  static const double fac_[];

  static double facT(int n)
  {
    // NOTE: in most cases, subscript n won't overflow. 
    return fac_[n];
  }
};

/// @brief Geometry implementation of Bezier curve
class BezierCurveImpl
{
public:
  BezierCurveImpl() :degree(0) {}
  BezierCurveImpl(const BezierCurveImpl& rhs) noexcept
    :degree(rhs.degree), ctrl_points(rhs.ctrl_points)
  {}
  BezierCurveImpl(BezierCurveImpl&& rhs) noexcept
    :degree(rhs.degree), ctrl_points(std::move(rhs.ctrl_points))
  {}
  BezierCurveImpl& operator=(const BezierCurveImpl& rhs) noexcept
  {
    degree = rhs.degree;
    ctrl_points = rhs.ctrl_points;
    return *this;
  }
  BezierCurveImpl& operator=(BezierCurveImpl&& rhs) noexcept
  {
    degree = rhs.degree;
    ctrl_points = std::move(rhs.ctrl_points);
    return *this;
  }

  /* Properties */
  uint32_t degree;
  Points ctrl_points;
};

class IndexedBezierCurveImpl : public BezierCurveImpl
{
public:
  IndexedBezierCurveImpl() :BezierCurveImpl() {}
  IndexedBezierCurveImpl(const BezierCurveImpl& curve, Index _idx) :
    BezierCurveImpl(curve), idx(_idx)
  {}
  IndexedBezierCurveImpl(BezierCurveImpl&& curve, Index _idx) :
    BezierCurveImpl(std::move(curve)), idx(_idx)
  {}
  IndexedBezierCurveImpl(const IndexedBezierCurveImpl& rhs) noexcept
    :BezierCurveImpl(rhs), idx(rhs.idx)
  {}
  IndexedBezierCurveImpl(IndexedBezierCurveImpl&& rhs) noexcept
    :BezierCurveImpl(std::move(rhs)), idx(rhs.idx)
  {}
  IndexedBezierCurveImpl& operator=(const IndexedBezierCurveImpl& rhs) noexcept
  {
    degree = rhs.degree;
    ctrl_points = rhs.ctrl_points;
    idx = rhs.idx;
    return *this;
  }
  IndexedBezierCurveImpl& operator=(IndexedBezierCurveImpl&& rhs) noexcept
  {
    degree = rhs.degree;
    ctrl_points = std::move(rhs.ctrl_points);
    idx = rhs.idx;
    return *this;
  }

  /* Properties */
  Index idx;
};

}// namespace Geometry
}// namespace GCLF