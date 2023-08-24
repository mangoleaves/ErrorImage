#include "CCheckConstraints.h"

namespace GCLF
{
namespace ImageTriSimp
{

/// @brief Check if all faces of mesh are injective.
/// @param mesh The Bezier mesh.
/// @retval true: if all faces are injective.
/// @retval false: if any face is not injective.
bool check_injectivity(const BMeshT* mesh)
{
  BezierTriangleInjectivityChecker checker;
  checker.set_degree(mesh->bezier_degree);

  for (FaceHandle fh : mesh->faces())
  {
    if (mesh->status(fh).deleted())
      continue;
    const BezierFace& bface = mesh->data(fh).bface;
    if (checker.is_definitely_injective_by_determinant(bface.ctrl_pnt) != 1)
      return false;
  }
  return true;
}

static std::vector<Vec3d> uvw_samples;

void initialize_uvw_samples(const uint32_t degree)
{
  uvw_samples.clear();
  // Quadrature for Energy computation
  int quadrature = 36 * (degree - 1);
  quadrature = (std::sqrt(8 * quadrature + 1) - 1.0) / 2.0;
  if (quadrature < 1)
    quadrature = 1;
  for (int i = 0; i < quadrature; i++)
  {
    for (int j = 0; j < quadrature - i; j++)
    {
      if (1 == degree)
      {
        uvw_samples.emplace_back(0.0, 0.0, 1.0);
        i = quadrature + 1;
        break;
      }
      uvw_samples.emplace_back(
        (1.0 * i) / (quadrature - 1.0),
        (1.0 * j) / (quadrature - 1.0),
        1.0 - ((1.0 * i) / (quadrature - 1.0) + (1.0 * j) / (quadrature - 1.0)));
    }
  }
}

double MIPS_energy(const Points& ctrl_pnts, const uint32_t degree)
{
  std::vector<Eigen::Matrix2d> Jacobians;
  // TODO: optimize below codes.
  for (const Vec3d& uvw : uvw_samples)
    Jacobians.push_back(BezierTriangle::Jacobian_tau_2D(1., ctrl_pnts, degree, uvw));

  double MIPS = 0.0;
  for (const Eigen::Matrix2d& J : Jacobians)
  {
    const double& a = J(0, 0);
    const double& b = J(0, 1);
    const double& c = J(1, 0);
    const double& d = J(1, 1);
    double det = (a * d) - (b * c);
    double Frobenius_norm = (a * a) + (b * b) + (c * c) + (d * d);
    MIPS += Frobenius_norm / pow(det * det, 1. / 3.);
    // MIPS += Frobenius_norm / det;
  }
  return MIPS / uvw_samples.size(); // uniform weights for all samples
}

/// @brief Approximate curvature of a quadratic bezier curve.
/// A quadratic bezier curve: from----mid----to.
/// mf = normalize(from - mid), mt = normalize(to - mid).
/// return cos(theta) = mf.dot(mt).
/// cos(theta) is in [-1,1], the closer the value is to -1, the flatter the curve.
double approx_bezier_curvature(const BMeshT* mesh, EdgeHandle eh)
{
  ASSERT(mesh->bezier_degree == 2, "degree is not 2.");
  HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
  VertexHandle from_vh = mesh->from_vertex_handle(hh);
  VertexHandle to_vh = mesh->to_vertex_handle(hh);
  const Vec3d& from_p = mesh->point(from_vh);
  const Vec3d& to_p = mesh->point(to_vh);
  const Vec3d& mid_p = mesh->point(eh);

  Vec3d mf = (from_p - mid_p).normalized();
  Vec3d mt = (to_p - mid_p).normalized();
  return (mf.x() * mt.x() + mf.y() * mt.y());
}

double approx_bezier_sector_angle(const BMeshT* mesh, HalfedgeHandle hh)
{
  ASSERT(mesh->bezier_degree == 2, "degree is not 2.");
  VertexHandle to_vh = mesh->to_vertex_handle(hh);
  HalfedgeHandle next_hh = mesh->next_halfedge_handle(hh);
  const Vec3d& to_p = mesh->point(to_vh);
  const Vec3d& cur_p = mesh->point(mesh->edge_handle(hh));
  const Vec3d& next_p = mesh->point(mesh->edge_handle(next_hh));

  Vec3d ct = (cur_p - to_p).normalized();
  Vec3d nt = (next_p - to_p).normalized();
  return (ct.x() * nt.x() + ct.y() * nt.y());
}

}// namespace ImageTriSimp
}// namespace GCLF