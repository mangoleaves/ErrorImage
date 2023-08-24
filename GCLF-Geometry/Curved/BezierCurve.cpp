#include "BezierCurve.h"

namespace GCLF
{
namespace Geometry
{

const double BezierCurve::fac_[] =
//0    1    2    3    4     5      6      7       8       9        10
{ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320., 362880., 3628800. };

Vec3d BezierCurve::deCasteljau(const Points& ctrl_points, const uint32_t degree, const double t)
{
  Points intermediate_points = ctrl_points;
  for (uint32_t n = 0;n < degree;n++)
  {
    for (uint32_t i = 0;i < degree - n;i++)
    {
      intermediate_points[i] =
        (1. - t) * intermediate_points[i] +
        t * intermediate_points[i + 1];
    }
  }
  return intermediate_points[0];
}

Vec3d BezierCurve::eval_by_Berstain(const Points& ctrl_points, const uint32_t degree, const double t)
{
  ASSERT(degree <= 10, "too large degree");

  // std::vector<double> coeff;  coeff.reserve(ctrl_points.size());
  // for (uint32_t i = 0;i <= degree;i++)
  //   coeff.push_back((facT(degree) / (facT(i) * facT(degree - i)))
  //     * pow(1 - t, degree - i) * pow(t, i));
  double coeff[11];
  for (uint32_t i = 0;i <= degree;i++)
    coeff[i] = (facT(degree) / (facT(i) * facT(degree - i)))
      * pow(1 - t, degree - i) * pow(t, i);

  Vec3d result(0., 0., 0.);
  for (uint32_t i = 0;i <= degree;i++)
    result += coeff[i] * ctrl_points[i];

  return result;
}

/// @brief Linear approximation of a Bezier curve.
/// @param density approximation density.
/// @return A polyline.
Points BezierCurve::linear_approx(const Points& ctrl_points, const uint32_t degree, uint32_t density)
{
  ASSERT(density >= 1, "insufficiant density.");

  Points polyline;
  double step = 1.0 / density;
  double t = 0.;
  
  for (uint32_t i = 0;i <= density;i++)
  {
    polyline.push_back(eval_by_Berstain(ctrl_points, degree, t));
    t += step;
  }
  return polyline;
}

}// namespace Geometry
}// namespace GCLF