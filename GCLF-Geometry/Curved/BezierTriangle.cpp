#include "BezierTriangle.h"
#include <limits>
#include <numeric>


namespace GCLF
{
namespace Geometry
{

/// @brief Map parameter coordinates to barycentric coordinates.
/// @details Use a standard unit 2D triangle as the parameter domain,
/// i.e., \f( u \in [0,1] , v \in [0,1], u+v \in [0,1] \f).
/// When mappint to the barycentric coordinates uvw,
/// we calculate the sub-triangle's area and multiply it with 2.
/// Finally, we have \f(  u <- u, v <- v, w <- 1-u-v \f).
/// @todo optimize this latter.
Vec3d BezierTriangle::map_param_coord_to_bary_coord(const Vec2d& uv)
{
  Vec3d uvw;
  uvw[0] = uv[0];
  uvw[1] = uv[1];
  uvw[2] = 1. - uvw[0] - uvw[1];
  return uvw;
}

/// @brief Find the index of corner point of Bezier triangle.
/// @todo pre-calculate the indices for some degrees.
Index BezierTriangle::index_of_vert(uint32_t degree, Index vert_idx)
{
  switch (vert_idx)
  {
  case 0:return 0;
  case 1:return (degree * (degree + 1)) / 2;
  case 2:return (degree * (degree + 1)) / 2 + degree;
  default:ASSERT(false, "wrong vertex index.");break;
  }
  return -1;
}

/// @brief Find the indices of points on one edge of Bezier triangle.
/// @details The 0-th edge is from the 0-th point to the n(n+1)/2-th point.
/// The 1-st edge: n(n+1)/2-th     ->  n(n+1)/2+n-th point.
/// The 2-st edge: n(n+1)/2+n-th   ->  0-th point.
/// @todo pre-calculate the indices for some degrees.
Indices BezierTriangle::indices_on_edge(uint32_t degree, Index edge_idx)
{
  Indices indices; indices.reserve(degree);
  switch (edge_idx)
  {
  case 0:
  {
    for (Index i = 0;i <= (Index)degree;i++)
    {
      indices.push_back(map_mat_idx_to_vec_idx(i, 0));
    }
  }break;
  case 1:
  {
    for (Index j = 0;j <= (Index)degree;j++)
    {
      indices.push_back(map_mat_idx_to_vec_idx(degree, j));
    }
  }break;
  case 2:
  {
    for (Index i = degree;i >= 0;i--)
    {
      indices.push_back(map_mat_idx_to_vec_idx(i, i));
    }
  }break;
  default:ASSERT(false, "wrong edge index.");break;
  }
  return indices;
}

/// @brief Find the indices of points inside one edge of Bezier triangle.
/// @details The 0-th edge is from the 0-th point to the n(n+1)/2-th point.
/// The 1-st edge: n(n+1)/2-th     ->  n(n+1)/2+n-th point.
/// The 2-st edge: n(n+1)/2+n-th   ->  0-th point.
/// @todo pre-calculate the indices for some degrees.
Indices BezierTriangle::indices_inside_edge(uint32_t degree, Index edge_idx)
{
  Indices indices; indices.reserve(degree);
  switch (edge_idx)
  {
  case 0:
  {
    for (Index i = 1;i < (Index)degree;i++)
    {
      indices.push_back(map_mat_idx_to_vec_idx(i, 0));
    }
  }break;
  case 1:
  {
    for (Index j = 1;j < (Index)degree;j++)
    {
      indices.push_back(map_mat_idx_to_vec_idx(degree, j));
    }
  }break;
  case 2:
  {
    for (Index i = degree - 1;i > 0;i--)
    {
      indices.push_back(map_mat_idx_to_vec_idx(i, i));
    }
  }break;
  default:ASSERT(false, "wrong edge index.");break;
  }
  return indices;
}

Indices BezierTriangle::indices_not_on_edge(uint32_t degree, Index edge_idx)
{
  Indices indices; indices.reserve((degree * (degree + 1)) / 2);
  switch (edge_idx)
  {
  case 0:
  {
    for (Index i = 1;i <= (Index)degree;i++)
    {
      for (Index j = 1;j <= i;j++)
      {
        indices.push_back(map_mat_idx_to_vec_idx(i, j));
      }
    }
  }break;
  case 1:
  {
    for (Index i = 0;i < (Index)degree;i++)
    {
      for (Index j = 0;j <= i;j++)
      {
        indices.push_back(map_mat_idx_to_vec_idx(i, j));
      }
    }
  }break;
  case 2:
  {
    for (Index i = 1;i <= (Index)degree;i++)
    {
      for (Index j = 0;j < i;j++)
      {
        indices.push_back(map_mat_idx_to_vec_idx(i, j));
      }
    }
  }break;
  default:ASSERT(false, "wrong edge index.");break;
  }
  return indices;
}

/// @brief Find the indices of corner points of Bezier triangle.
Indices BezierTriangle::indices_corner(uint32_t degree)
{
  return Indices({
    0,
    (int)((degree * (degree + 1)) / 2),
    (int)((degree * (degree + 1)) / 2 + degree)
    });
}

/// @brief Find the indices of interior points of Bezier triangle.
Indices BezierTriangle::indices_interior(uint32_t degree)
{
  Indices indices; indices.reserve(degree * degree);
  for (Index i = 1;i < (Index)degree;i++)
  {
    Index v = map_mat_idx_to_vec_idx(i, 1);
    for (Index j = 1;j < i;j++, v++)
      indices.push_back(v);
  }
  return indices;
}


/// @brief Find the indices of all points of Bezier triangle.
Indices BezierTriangle::indices_all(uint32_t degree)
{
  Indices indices; indices.resize(((degree + 1) * (degree + 2)) / 2);
  std::iota(indices.begin(), indices.end(), 0);
  return indices;
}

/// @brief Implement the de Casteljau algorithm and evaluate points on the bezier triangle.
/// @param uvw The barycentric coordinate of the point.
/// @return Evaluation result.
Vec3d BezierTriangle::deCasteljau(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw)
{
  // check the barycentric coordinates
  {
    ASSERT(uvw.isfinite() && uvw.all_geq(-bc_eps), "wrong input coordinates for de Casteljau.");
    const double sum = uvw.sum();
    ASSERT(1. - bc_eps <= sum && sum <= 1. + bc_eps, "wrong input coordinates for de Casteljau.");
  }

  // TODO: Use a global variable to store intermediate points.
  // It avoids allocating memory repeatedly.

  Points intermediate_points = ctrl_points;

  for (int n = (int)degree;n > 0;n--)
  {
    for (Index i = 0;i < n;i++)    // i-th row
    {
      Index s = map_mat_idx_to_vec_idx(i, 0);
      Index t = map_mat_idx_to_vec_idx(i + 1, 0);

      for (Index j = 0;j <= i;j++, s++, t++)       // j-th column
      {
        intermediate_points[s] =
          intermediate_points[t + 1] * uvw[0] +
          intermediate_points[s] * uvw[1] +
          intermediate_points[t] * uvw[2];
      }
    }
  }
  return intermediate_points[0];
}

const double BezierTriangle::fac_[] =
//0    1    2    3    4     5      6      7       8       9        10
{ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320., 362880., 3628800. };

/// @brief Calculate the Berstein polynomial, with given degree, polar index and position.
/// @note the relation between polar index and parameter is (u^j * v^i * w^k).
double BezierTriangle::Berstein_polynomial(const uint32_t degree, const Index3& ijk, const Vec3d& uvw)
{
  if (ijk.x() < 0 || ijk.y() < 0 || ijk.z() < 0)
    return 0.0;

  return facT(degree) * pow(uvw.x(), ijk.y()) * pow(uvw.y(), ijk.x()) * pow(uvw.z(), ijk.z())
    / (facT(ijk.x()) * facT(ijk.y()) * facT(ijk.z()));
}

Vec3d BezierTriangle::eval_by_Berstein(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw)
{
  Vec3d result(0, 0, 0);
  for (Index i = 0;i <= (Index)degree;i++)
  {
    for (Index j = 0;j <= i;j++)
    {
      const Vec3d& cp = ctrl_points[map_mat_idx_to_vec_idx(i, j)];
      Index3 ijk = map_mat_idx_to_polar(i, j, degree);
      result += cp * Berstein_polynomial(degree, ijk, uvw);
    }
  }
  return result;
}

Points BezierTriangle::degree_elevation(const Points& ctrl_points, const uint32_t degree)
{
  // we call the degree elevated bezier tri as "up".
  // we call the original bezier tri as "low".
  Points up_ctrl_points;
  const Points& low_ctrl_points = ctrl_points;
  up_ctrl_points.reserve(((degree + 1) * (degree + 2)) / 2);

  for (Index i = 0;i <= (Index)degree + 1;i++)
  {
    for (Index j = 0;j <= i;j++)
    {
      // here, (i, j) is the index of up bezier triangle whose degree is increased by 1.
      Index3 up_polar = map_mat_idx_to_polar(i, j, degree + 1);

      Vec3d up_point(0., 0., 0.);
      // see degree elevation formular in polar form.
      for (Index k = 0;k < 3;k++)
      {
        if (up_polar[k] > 0)
        {
          Index3 low_polar = up_polar;
          low_polar[k] -= 1;
          Index2 low_ij = map_polar_to_mat_idx(low_polar, degree);
          up_point += low_ctrl_points[map_mat_idx_to_vec_idx(low_ij)]
            * (double)up_polar[k];
        }
      }
      up_point /= (double)(degree + 1);
      up_ctrl_points.push_back(up_point);
    }
  }

  return up_ctrl_points;
  // TODO: optimize those Map functions.
}

/// @brief Apply subdivision_1_3 operation on a Bezier triangle (one Bezier triangle is subdivided to three Bezier triangles).
/// @param [in] uvw The barycentric coordinate of subdivided center.
/// @param [out] sub_ctrl_points Three groups of control points of three sub-triangles.
/// @param [out] sub_center_ctrl_pnt Subdivide point.
void BezierTriangle::subdivision_1_3(
  const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw,
  std::array<Points, 3>& sub_ctrl_points, Vec3d& sub_center_ctrl_pnt)
{
  // check the barycentric coordinates
  {
    ASSERT(uvw.isfinite() && uvw.all_geq(-bc_eps), "wrong input coordinates for subdivision.");
    const double sum = uvw.sum();
    ASSERT(1. - bc_eps <= sum && sum <= 1. + bc_eps, "wrong input coordinates for subdivision.");
  }

  // TODO: Use a global variable to store intermediate points.
  // It avoids allocating memory repeatedly.

  Points intermediate_points = ctrl_points;
  for (size_t i = 0;i < 3;i++)
    sub_ctrl_points[i].resize(ctrl_points.size());

  for (uint8_t n = degree;n > 0;n--)
  {
    // put intermediate control points into subdivided triangles.
    Index diff = degree - n;
    for (Index i = 0;i <= n;i++)
    {
      Index src = map_mat_idx_to_vec_idx(i, 0);
      Index dst = map_mat_idx_to_vec_idx(i + diff, diff);
      sub_ctrl_points[0][dst] = intermediate_points[src];
    }
    for (Index j = 0;j <= n;j++)
    {
      Index src_dst = map_mat_idx_to_vec_idx(n, j);
      sub_ctrl_points[1][src_dst] = intermediate_points[src_dst];
    }
    for (Index i = 0;i <= n;i++)
    {
      Index src = map_mat_idx_to_vec_idx(i, i);
      Index dst = map_mat_idx_to_vec_idx(i + diff, i);
      sub_ctrl_points[2][dst] = intermediate_points[src];
    }
    // calculate intermediate control points
    for (Index i = 0;i < n;i++)    // i-th row
    {
      Index s = map_mat_idx_to_vec_idx(i, 0);
      Index t = map_mat_idx_to_vec_idx(i + 1, 0);

      for (Index j = 0;j <= i;j++, s++, t++)       // j-th column
      {
        intermediate_points[s] =
          intermediate_points[t + 1] * uvw[0] +
          intermediate_points[s] * uvw[1] +
          intermediate_points[t] * uvw[2];
      }
    }
  }
  // put the last one control points into three triangles.
  sub_ctrl_points[0][map_mat_idx_to_vec_idx(degree, degree)] = intermediate_points[0];
  sub_ctrl_points[1][0] = intermediate_points[0];
  sub_ctrl_points[2][map_mat_idx_to_vec_idx(degree, 0)] = intermediate_points[0];
  sub_center_ctrl_pnt = intermediate_points[0];
}

/// @brief Apply subdivision_1_2 operation on a Bezier triangle (one Bezier triangle is subdivided to two Bezier triangles).
/// @param [in] local_edge_idx The edge of Bezier triangle to split.
/// @param [in] uv The barycentric coordinate of subdivided center on the edge.
/// @param [out] sub_ctrl_points Three groups of control points of three sub-triangles.
/// @param [out] sub_center_ctrl_pnt Subdivide point.
void BezierTriangle::subdivision_1_2(
  const Points& ctrl_points, const uint32_t degree, const Index local_edge_idx, const Vec2d& uv,
  std::array<Points, 2>& sub_ctrl_points, Vec3d& sub_split_ctrl_pnt)
{
  // check the barycentric coordinates
  {
    ASSERT(uv.isfinite() && uv.all_geq(-bc_eps), "wrong input coordinates for subdivision at edge.");
    const double sum = uv.sum();
    ASSERT(1. - bc_eps <= sum && sum <= 1. + bc_eps, "wrong input coordinates for subdivision at edge.");
  }
  std::array<Points, 3> _sub_ctrl_points;
  Vec3d uvw;
  switch (local_edge_idx)
  {
  case 0:
  {
    uvw = { 0, uv[0], uv[1] };
    subdivision_1_3(ctrl_points, degree, uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
    sub_ctrl_points[0] = _sub_ctrl_points[2];
    sub_ctrl_points[1] = _sub_ctrl_points[1];
  }break;
  case 1:
  {
    uvw = { uv[1], 0, uv[0] };
    subdivision_1_3(ctrl_points, degree, uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
    sub_ctrl_points[0] = _sub_ctrl_points[0];
    sub_ctrl_points[1] = _sub_ctrl_points[2];
  }break;
  case 2:
  {
    uvw = { uv[0], uv[1], 0 };
    subdivision_1_3(ctrl_points, degree, uvw, _sub_ctrl_points, sub_split_ctrl_pnt);
    sub_ctrl_points[0] = _sub_ctrl_points[1];
    sub_ctrl_points[1] = _sub_ctrl_points[0];
  }break;
  }
}

/// @brief Apply subdivision_1_4 operation on a Bezier triangle (one Bezier triangle is subdivided to four Bezier triangles).
/// We subdivide the triangle at three mid-points of edges, so, we don't need barycenter coordinates.
/// @note Ask Guo Jiapeng for a figure that illustrates the algorithm.
void BezierTriangle::subdivision_1_4(
  const Points& ctrl_points, const uint32_t degree,
  std::array<Points, 4>& sub_ctrl_points)
{
  for (size_t i = 0;i < 4;i++)
  {
    sub_ctrl_points[i].clear();
    sub_ctrl_points[i].reserve(ctrl_points.size());
  }

  for (uint32_t i = 0;i <= degree;i++)
  {
    Index2 ij = map_polar_to_mat_idx(Index3(degree - i, 0, 0), degree - i);
    for (uint32_t n = 0;n <= i;n++)
      sub_ctrl_points[0].push_back(sub_control_point(ctrl_points, ij[0], ij[1], i - n, 0, n));

    ij = map_polar_to_mat_idx(Index3(0, 0, degree - i), degree - i);
    for (uint32_t n = 0;n <= i;n++)
      sub_ctrl_points[1].push_back(sub_control_point(ctrl_points, ij[0], ij[1], n, i - n, 0));

    ij = map_polar_to_mat_idx(Index3(0, degree - i, 0), degree - i);
    for (uint32_t n = 0;n <= i;n++)
      sub_ctrl_points[2].push_back(sub_control_point(ctrl_points, ij[0], ij[1], 0, n, i - n));
  }

  for (uint32_t i = 0;i <= degree;i++)
  {
    for (uint32_t j = 0;j <= i;j++)
    {
      // k0,k1,k2 are slightly different from polar index.
      sub_ctrl_points[3].push_back(sub_control_point(ctrl_points, 0, 0, degree - i, i - j, j));
    }
  }
}


/// @brief A recursion formular to compute sub control points for 1-4 subdivision.
/// @param i,j matrix index
/// @param k0,k1,k2 polar index of sub control point (only used in this subdivision algorithm)
/// @return the wanted control point.
Vec3d BezierTriangle::sub_control_point(
  const Points& ctrl_points,
  Index i, Index j, uint32_t k0, uint32_t k1, uint32_t k2)
{
  if (k0 > 0)
  {
    return 0.5 * (sub_control_point(ctrl_points, i, j, k0 - 1, k1, k2) +
      sub_control_point(ctrl_points, i + 1, j, k0 - 1, k1, k2));
  }
  if (k1 > 0)
  {
    return 0.5 * (sub_control_point(ctrl_points, i + 1, j, 0, k1 - 1, k2) +
      sub_control_point(ctrl_points, i + 1, j + 1, 0, k1 - 1, k2));
  }
  if (k2 > 0)
  {
    return 0.5 * (sub_control_point(ctrl_points, i + 1, j + 1, 0, 0, k2 - 1) +
      sub_control_point(ctrl_points, i, j, 0, 0, k2 - 1));
  }
  return ctrl_points[map_mat_idx_to_vec_idx(i, j)];
}


/// @brief Get the gradient vectors in the direction indicated by local_edge_idx.
/// @return gradient vectors
Points BezierTriangle::gradient_vectors(
  const Points& ctrl_points, const uint32_t degree, const Index local_edge_idx)
{
  Points vectors; vectors.reserve((degree * (degree + 1)) / 2);
  switch (local_edge_idx)
  {
  case 0:
  {
    for (Index i = 1;i <= (Index)degree;i++)
    {
      Index s = map_mat_idx_to_vec_idx(i, 0);
      for (Index j = 0;j < i;j++, s++)
        vectors.push_back(ctrl_points[s + 1] - ctrl_points[s]);
    }
  }break;
  case 1:
  {
    for (Index i = 1;i <= (Index)degree;i++)
    {
      Index s = map_mat_idx_to_vec_idx(i, 0);
      Index t = map_mat_idx_to_vec_idx(i - 1, 0);
      for (Index j = 0;j < i;j++, s++, t++)
        vectors.push_back(ctrl_points[t] - ctrl_points[s]);
    }
  }break;
  case 2:
  {
    for (Index i = 1;i <= (Index)degree;i++)
    {
      Index s = map_mat_idx_to_vec_idx(i, 0);
      for (Index j = 0;j < i;j++, s++)
        vectors.push_back(ctrl_points[s] - ctrl_points[s + 1]);
    }
  }break;
  default:ASSERT(false, "wrong local edge index.");
  }
  return vectors;
}

/// @brief The map Psi maps a regular triangle to a right triangle.
/// The regular triangle's edge length is given in parameter.
/// The right triangle has unit length legs (coincident with the coordinate axes u,v).
/// (Refer to BezierGuarding's appendix A.1)
/// @param edge_len the edge length of the regular triangle
/// @return The 2D Jacobian matrix of map Psi
Eigen::Matrix2d BezierTriangle::Jacobian_Psi_2D(const double edge_len)
{
  // a = edge_len
  // J_Psi = inverse ([ a       a / 2    ])
  //                 ([ 0   sqrt(3)a / 2 ])
  Eigen::Matrix2d mat;
  mat << 1., -1. / SQRT_3, 0, 2. / SQRT_3;
  mat /= edge_len;
  return mat;
}

/// @brief The map Phi maps a right triangle to a Bezier triangle.
/// The right triangle has unit length legs (coincident with the coordinate axes u,v).
/// @return The 2D Jacobian matrix of map Phi.
Eigen::Matrix2d BezierTriangle::Jacobian_Phi_2D(const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw)
{
  // J_Phi (u, v) = [ dx / du,  dx / dv ]
  //                [ dy / du,  dy / dv ]

  // calculate the gradient vectors along the direction (u,v) = (1,0) and (u,v) = (0,1)
  Points gv_u = gradient_vectors(ctrl_points, degree, 0);
  Points gv_v = gradient_vectors(ctrl_points, degree, 1);

  Vec3d gradient_u = eval_by_Berstein(gv_u, degree - 1, uvw);
  Vec3d gradient_v = eval_by_Berstein(gv_v, degree - 1, uvw);

  Eigen::Matrix2d mat;
  mat << gradient_u.x(), gradient_v.x(),
    gradient_u.y(), gradient_v.y();
  return mat;
}

/// @brief The map tau maps a regular triangle to a Bezier triangle.
/// Thre regular triangle' edge length is given in parameter.
/// @return The 2D Jacobian matrix of map tau.
Eigen::Matrix2d BezierTriangle::Jacobian_tau_2D(
  const double edge_len, const Points& ctrl_points, const uint32_t degree, const Vec3d& uvw)
{
  Eigen::Matrix2d J_Psi = Jacobian_Psi_2D(edge_len);
  Eigen::Matrix2d J_Phi = Jacobian_Phi_2D(ctrl_points, degree, uvw);
  return J_Phi * J_Psi;
}

std::vector<Index3> BezierTriangle::get_ctrl_triangles(const uint32_t degree)
{
  std::vector<Index3> ctrl_triangles;
  ctrl_triangles.reserve(degree * degree);
  for (Index i = 0;i < (Index)degree;i++)    // i-th row
  {
    Index s = map_mat_idx_to_vec_idx(i, 0);
    Index t = map_mat_idx_to_vec_idx(i + 1, 0);

    for (Index j = 0;j < i;j++, s++, t++)       // j-th column
    {
      ctrl_triangles.emplace_back(Index3(s, t, t + 1));
      ctrl_triangles.emplace_back(Index3(s, t + 1, s + 1));
    }
    ctrl_triangles.emplace_back(Index3(s, t, t + 1)); // for the case when j==i
  }
  return ctrl_triangles;
}

std::vector<Index2> BezierTriangle::get_junctions_of_ctrl_triangles(const uint32_t degree)
{
  std::vector<Index2> junctions;
  junctions.reserve(3 * degree);

  for (Index i = 0;i < (Index)degree;i++)
  {
    Index s = map_mat_idx_to_vec_idx(i, 0);
    Index t = map_mat_idx_to_vec_idx(i + 1, 0);
    junctions.emplace_back(Index2(s, t));
  }

  for (Index j = 0;j < (Index)degree;j++)
  {
    Index s = map_mat_idx_to_vec_idx(degree, j);
    Index t = map_mat_idx_to_vec_idx(degree, j + 1);
    junctions.emplace_back(Index2(s, t));
  }

  for (Index i = (Index)degree;i > 0;i--)
  {
    Index s = map_mat_idx_to_vec_idx(i, i);
    Index t = map_mat_idx_to_vec_idx(i - 1, i - 1);
    junctions.emplace_back(Index2(s, t));
  }
  return junctions;
}

/// @brief Calculate a triangle mesh to represent the linear approximation of a Bezier triangle.
/// @param ctrl_points Control points.
/// @param degree Degree of Bezier triangle.
/// @param density The approximation density. The higher density, the lower approximation error.
/// @param calc_topo If the "triangles" does not change, we do not calculate it repeateddly.
/// @param points Vertices of the mesh.
/// @param triangles Triangle of the mesh. The indices point to the array "points".
void BezierTriangle::linear_approximation(
  const Points& ctrl_points, const uint32_t degree,
  const uint32_t density, const bool calc_topo,
  Points& points, std::vector<Index3>& triangles)
{
  ASSERT(density >= 1, "insufficiant density.");


  double step = 1.0 / density;
  Vec2d param_coord;
  Vec3d bary_coord;

  // calculate points
  points.clear();
  points.reserve(((density + 1) * (density + 2)) / 2);
  for (Index i = 0;i <= (Index)density;i++)
  {
    for (Index j = 0;j <= i;j++)
    {
      param_coord[0] = j * step;                  // param u
      param_coord[1] = (density - i) * step;  // param v
      bary_coord = map_param_coord_to_bary_coord(param_coord);
      points.push_back(deCasteljau(ctrl_points, degree, bary_coord));
    }
  }

  if (!calc_topo)
    return;

  // calculate triangles
  triangles.clear();
  triangles.reserve(density * density);
  for (Index i = 0;i < (Index)density;i++)    // i-th row
  {
    Index s = map_mat_idx_to_vec_idx(i, 0);
    Index t = map_mat_idx_to_vec_idx(i + 1, 0);

    for (Index j = 0;j < i;j++, s++, t++)       // j-th column
    {
      triangles.emplace_back(Index3(s, t, t + 1));
      triangles.emplace_back(Index3(s, t + 1, s + 1));
    }
    triangles.emplace_back(Index3(s, t, t + 1)); // for the case when j==i
  }
}


}// namespace Geometry
}// namespace GCLF