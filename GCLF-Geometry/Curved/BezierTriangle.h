#pragma once
#include "BezierCommon.h"

namespace GCLF
{
namespace Geometry
{
using namespace Utils;

class BezierTriangle
{
public:

  /* Index Map & Conversion */

  static Vec3d map_param_coord_to_bary_coord(const Vec2d& uv);

  // map (i, j) in lower triangle matrix to the index in linear array.
  static inline Index map_mat_idx_to_vec_idx(const Index i, const Index j) { return (i * (i + 1)) / 2 + j; }
  static inline Index map_mat_idx_to_vec_idx(const Index2& ij) { return map_mat_idx_to_vec_idx(ij[0], ij[1]); }

  // map (i, j) in lower triangle matrix to the polar form (i, j, k). n is degree.
  static inline Index3 map_mat_idx_to_polar(const Index i, const Index j, const uint32_t n){ return Index3(n - i, j, i - j); }
  static inline Index3 map_mat_idx_to_polar(const Index2& ij, const uint32_t n){ return map_mat_idx_to_polar(ij[0], ij[1], n); }

  // map polar form to index in lower triangle matrix. n is degree.
  static inline Index2 map_polar_to_mat_idx(const Index3& ijk, const uint32_t n){ return Index2(n - ijk[0], ijk[1]); }

  /*  Traverse  */

  static Index index_of_vert(uint32_t degree, Index vert_idx);

  static Indices indices_on_edge(uint32_t degree, Index edge_idx);

  static Indices indices_inside_edge(uint32_t degree, Index edge_idx);

  static Indices indices_not_on_edge(uint32_t degree, Index edge_idx);

  static Indices indices_corner(uint32_t degree);

  static Indices indices_interior(uint32_t degree);

  static Indices indices_all(uint32_t degree);

  /* Evaluation */

  static Vec3d deCasteljau(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw);

  static const double fac_[];
  static double facT(int n)
  {
    // NOTE: in most cases, subscript n won't overflow. 
    return fac_[n];
  }
  static double Berstein_polynomial(const uint32_t degree, const Index3& ijk, const Vec3d& uvw);
  static Vec3d eval_by_Berstein(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw);

  /* Elevation */

  static Points degree_elevation(const Points& ctrl_points, const uint32_t degree);

  /* Subdivision */

  static void subdivision_1_3(
    const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw,
    std::array<Points, 3>& sub_ctrl_points, Vec3d& sub_center_ctrl_pnt);

  static void subdivision_1_2(
    const Points& ctrl_points, const uint32_t degree, const Index local_edge_idx, const Vec2d& uv,
    std::array<Points, 2>& sub_ctrl_points, Vec3d& sub_split_ctrl_pnt);

  static void subdivision_1_4(
    const Points& ctrl_points, const uint32_t degree,
    std::array<Points, 4>& sub_ctrl_points);

  static Vec3d sub_control_point(
    const Points& ctrl_points,
    Index i, Index j, uint32_t k0, uint32_t k1, uint32_t k2);

  /* Gradient & Jacobian */

  static Points gradient_vectors(
    const Points& ctrl_points, const uint32_t degree, const Index local_edge_idx);

  static Eigen::Matrix2d Jacobian_Psi_2D(const double edge_len);

  static Eigen::Matrix2d Jacobian_Phi_2D(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw);

  static Eigen::Matrix2d Jacobian_tau_2D(const double edge_len, const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw);

  /* Visualization */

  static std::vector<Index3> get_ctrl_triangles(const uint32_t degree);

  static std::vector<Index2> get_junctions_of_ctrl_triangles(const uint32_t degree);

  static void linear_approximation(
    const Points& ctrl_points, const uint32_t degree,
    const uint32_t density, const bool calc_topo,
    std::vector<Vec3d>& points, std::vector<Index3>& triangles);
};

/// @brief Geometry implementation of Bezier triangle
class BezierTriangleImpl
{
public:
  BezierTriangleImpl():degree(0) {}

  /* Properties */

  uint32_t degree;

  // A control points net with degree being 3.
  // We want to put control points into a linear array.
  // v=1 (300)
  //  *(0)
  //  *    *
  //  *    *    *
  //  *(6) *    *    *(9)  u=1 (030)
  // w=1 (003)
  Points ctrl_points;
};

}// namespace Geometry
}// namespace GCLF