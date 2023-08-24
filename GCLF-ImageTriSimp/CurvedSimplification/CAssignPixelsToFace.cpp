#include "CAssignPixelsToFace.h"
#include "CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

bool CAssignPixelsToFace::operator()(
  BMeshT* _mesh, ImageT* _image,
  Eigen::MatrixXi& assigned_mask,
  const std::vector<Vec2i>& affected_pixels)
{
  initialize(_mesh, _image, assigned_mask, affected_pixels);

  // Calculate the bounding box for mesh.
  BoundingBox3i ibox = box_for_mesh();

  // Emit rays from left to right, starting at the left boundary of the box, one ray for one row of image.
  for (int y = ibox.min().y();y < ibox.max().y();y++)
  {
    // Current ray
    Vec2d ray_start(ibox.min().x() - 0.5, y + 0.5);
    bool succeed_to_assign = false;
    do
    {
      try
      {
        assign_pixels(assigned_mask, ray_start);
        succeed_to_assign = true;
      }
      catch(std::logic_error _)
      {
        UNUSED(_);
        ray_start.y() += 0.02;
        succeed_to_assign = false;
        // Logger::dev_logger->debug("retry with {:.2f}", ray_start.y());
      }
      catch(...) { throw; }
    } while (!succeed_to_assign && ray_start.y() - y < 0.8);
    if (!succeed_to_assign)
    {
      Logger::dev_logger->warn("fail to assign pixels.");
      // throw std::logic_error("fail to assign pixels.");
    }
  }
  return true;
}


void CAssignPixelsToFace::initialize(
  BMeshT* _mesh, ImageT* _image,
  Eigen::MatrixXi& assigned_mask,
  const std::vector<Vec2i>& affected_pixels)
{
  // Mesh is possibly a local mesh when doing local operations.
  mesh = _mesh;
  image = _image;
  // Set some pixels in image to be assigned and
  // reset the assigned pixels in mesh.
  // Then we only assign unassigned pixels to mesh.
  for (const Vec2i& pixel : affected_pixels)
    assigned_mask(pixel.y(), pixel.x()) = 0;

  ASSERT(!mesh->ctrl_pnt.empty() && mesh->n_faces() > 0, "empty mesh.");

  // Clear assigned pixels.
  for (FaceHandle fh : mesh->faces())
    mesh->data(fh).pixels.clear();

  // Build an AABB tree for edges (Bezier curved) of mesh.
  tree = std::make_unique<BezierCurveTree>(*mesh);
}

void CAssignPixelsToFace::assign_pixels(
  Eigen::MatrixXi& assigned_mask,
  const Vec2d& ray_start)
{
  // Use tree to search possible intersected curves.
  auto possible_curves = tree->intersect_curves(ray_start);
  if (possible_curves.empty())
    return;

  // Find the intersection points between rays and edges.
  std::vector<InterPoint> inter_points = find_intersection_points(ray_start, possible_curves);
  if (inter_points.size() < 2)
    return;

  // Sort the intersection points on each ray.
  std::sort(inter_points.begin(), inter_points.end());

  // If some intersection points are close to one vertex,
  //  contract the one ring intersection points of the vertex to one intersection point.
  cluster_and_collapse_inter_points(inter_points);
  if (inter_points.size() < 2)
    return;

  // Enter faces and exist faces from left to right along rays.
  enter_faces_along_the_ray((int)ray_start.y(), inter_points, assigned_mask);
}

/// @brief Caclulate the bounding box for the mesh, and then
/// enlarge it to get an integer bounding box.
BoundingBox3i CAssignPixelsToFace::box_for_mesh()
{
  BoundingBox box(mesh->ctrl_pnt);
  box.min() = box.min().floor();
  box.max() = box.max().ceil();
  BoundingBox3i ibox;
  ibox.min().x() = (int)box.min().x();
  ibox.min().y() = (int)box.min().y();
  ibox.max().x() = (int)box.max().x();
  ibox.max().y() = (int)box.max().y();

  // if (ibox.min().y() < 0 || ibox.max().y() > image->height)
  //  Logger::dev_logger->warn("Mesh size exceeds the image size.");

  ibox.min().y() = std::max(0, ibox.min().y());
  ibox.max().y() = std::min(image->height, ibox.max().y());
  return ibox;
}

/// @brief Find the intersection case between a ray and a Bezier curve.
/// @param ray_start The ray is horizontal to x-axis, i.e., its direction is (1, 0).
/// @param bcurve The Bezier curve.
/// @param inter_type intersection type.
/// @param inter_point intersection point (if it exists).
/// @note: remember that our pixel's width is 1.
/// if you want to reuse this code, adjust the degeneration threshold.
void CAssignPixelsToFace::ray_Bezier_curve_intersect(
  const Vec2d& ray_start, const BezierCurveImpl& bcurve,
  RayBezierCurveInterType& inter_type, std::vector<double>& inter_params)
{
  typedef RayBezierCurveInterType InterType;
  // Symbols:
  // ray_start: Ps. ray_start components: (xs, ys)
  // control points: P0, P1, P2
  // control points components: (x0, y0), (x1, y1), (x2, y2)

  const Vec2d& Ps = ray_start;
  const Vec3d& P0 = bcurve.ctrl_points[0], & P1 = bcurve.ctrl_points[1], & P2 = bcurve.ctrl_points[2];

  // We want to solve the equation:
  // (1-t)^2 * y0 + 2*t*(1-t)*y1 + t^2 * y2 = ys
  // It is equivalent to:
  // t^2 * (y0 - 2y1 + y2) + t * 2 * (y1 - y0) + (y0 - ys) = 0
  // We simplify it to:
  // a * t^2 + b * t + c = 0

  double a = P0.y() - 2 * P1.y() + P2.y();
  double b = 2 * (P1.y() - P0.y());
  double c = P0.y() - Ps.y();

  if (fabs(a) < 1e-5)
  {
    // The quadratic equation degenerates to linear equation.
    if (fabs(b) < 1e-5)
    {
      // The linear equation also degenerates.
      inter_type = fabs(c) < 1e-5 ?
        InterType::NearlyParallel : InterType::NoIntersect;
      return;
    }
    else
    {
      // Solve a linear equation: b * t + c = 0.
      double t = -c / b;
      if (0. <= t && t <= 1.)
      {
        inter_type = InterType::One;
        inter_params.push_back(t);
      }
      else inter_type = InterType::NoIntersect;
      return;
    }
  }
  else
  {
    // We gonna solve the quadratic equation.
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0)
    {
      // There is no solution.
      inter_type = InterType::NoIntersect;
      return;
    }
    double sqrt_disc = std::sqrt(discriminant);
    if (fabs(sqrt_disc / (2 * a)) < 1e-7)
    {
      // The discriminant is near zero, we think there is only one solution.
      double t = -b / (2 * a);
      if (0. <= t && t <= 1.)
      {
        inter_type = InterType::Tangent;
        inter_params.push_back(t);
      }
      else inter_type = InterType::NoIntersect;
      return;
    }
    else
    {
      // There are two solutions.
      double t0 = (-b - sqrt_disc) / (2 * a);
      double t1 = (-b + sqrt_disc) / (2 * a);
      if (t0 > t1) std::swap(t0, t1);
      if (0. <= t0 && t0 <= 1.)
        inter_params.push_back(t0);
      if (0. <= t1 && t1 <= 1.)
        inter_params.push_back(t1);

      inter_type = inter_params.empty() ?
        InterType::NoIntersect :
        (inter_params.size() == 1 ?
          InterType::One : InterType::Two);
      return;
    }
  }
}

/// @brief Find intersection points between ray and curves.
/// @param possible_curves Curves found by AABB tree. They are possibly intersected with ray.
/// @return intersection points
std::vector<CAssignPixelsToFace::InterPoint>
CAssignPixelsToFace::find_intersection_points(
  const Vec2d& ray_start,
  const std::vector<IndexedBezierCurveImpl*>& possible_curves)
{
  typedef RayBezierCurveInterType InterType;
  typedef InterPoint::Type IPType;
  std::vector<InterPoint> inter_points;
  std::vector<double> inter_params; inter_params.reserve(2);
  for (auto curve : possible_curves)
  {
    InterType inter_type;
    inter_params.clear();
    ray_Bezier_curve_intersect(
      ray_start, *static_cast<BezierCurveImpl*>(curve),
      inter_type, inter_params);

    EdgeHandle eh(curve->idx);
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    auto [from_vh, to_vh] = from_to_vh(*mesh, hh);
    if (from_vh.idx() > to_vh.idx())
    {
      hh = mesh->opposite_halfedge_handle(hh);
      std::swap(from_vh, to_vh);
    }
    auto [from_p, to_p] = from_to_p(*mesh, hh);
    // from_p is also control_points[0], to_p is also control_points[2].
    ASSERT(from_p == curve->ctrl_points[0], "edge's control points error.");

    switch (inter_type)
    {
    case InterType::One:
    case InterType::Two:
    {
      for (double param : inter_params)
      {
        if (param == 0.0 || param == 1.0)
        {
          inter_points.emplace_back(
            InterPoint::Type::ThroughVertex,
            param == 0.0 ? from_p : to_p,
            mesh->InvalidEdgeHandle,
            param == 0.0 ? from_vh : to_vh);
        }
        else
        {
          Vec3d p = BezierCurve::eval_by_Berstain(curve->ctrl_points, curve->degree, param);
          VertexHandle closer_vh = (p - from_p).sqrnorm() < (p - to_p).sqrnorm() ? from_vh : to_vh;
          inter_points.emplace_back(IPType::ThroughEdge, p, eh, closer_vh);
        }
      }
    }
    break;
    case InterType::NearlyParallel:
    {
      if (from_p.x() < to_p.x())
        inter_points.emplace_back(IPType::ThroughVertex, from_p, eh, from_vh);
      else
        inter_points.emplace_back(IPType::ThroughVertex, to_p, eh, to_vh);
    }
    break;
    case InterType::Tangent:
    case InterType::NoIntersect: break;
    default: break;
    }
  }
  return inter_points;
}

/// @brief Some intersection points are too close to the same vertex,
/// we find such points and collapse them to the closest vertex,
/// and set the intersection type to ThroughVertex.
void CAssignPixelsToFace::cluster_and_collapse_inter_points(std::vector<InterPoint>& inter_points)
{
  typedef InterPoint::Type IPType;
  std::vector<uint8_t> contract_vertices(mesh->n_vertices(), 0);
  for (const InterPoint& ip : inter_points)
  {
    if (ip.type != IPType::ThroughEdge)
      continue;
    if (contract_vertices[ip.vh.idx()] != 1 &&
      (ip.point - mesh->point(ip.vh)).length() < 1e-4)
      contract_vertices[ip.vh.idx()] = 1;
  }
  for (InterPoint& ip : inter_points)
  {
    if (ip.type != IPType::ThroughEdge)
      continue;
    if (contract_vertices[ip.vh.idx()])
    {
      ip.type = IPType::ThroughVertex;
      ip.point = mesh->point(ip.vh);
      ip.eh.invalidate();
    }
  }
  auto discard_iter = std::unique(inter_points.begin(), inter_points.end());
  inter_points.erase(discard_iter, inter_points.end());
}

/// @brief Enter faces and exist faces from left to right along a ray.
/// Meanwhile assign pixels to entered faces.
/// @param y the ray start.
/// @param inter_points intersection points.
void CAssignPixelsToFace::enter_faces_along_the_ray(
  const int y, const std::vector<InterPoint>& inter_points,
  Eigen::MatrixXi& assigned_mask)
{
  typedef InterPoint::Type IPType;
  FaceHandle current_fh, next_fh;
  current_fh.invalidate(); next_fh.invalidate();

  // Traverse remain faces
  for (size_t i = 1; i < inter_points.size();i++)
  {
    const InterPoint& current_ip = inter_points[i - 1];
    const InterPoint& next_ip = inter_points[i];
    bool enter_boundary = false, exit_boundary = false;

    if (!current_fh.is_valid()) // we gonna enter a face
    {
      enter_boundary = true;
      if (current_ip.type == IPType::ThroughEdge)
      {
        HalfedgeHandle hh = mesh->halfedge_handle(current_ip.eh, 0);
        current_fh = mesh->is_boundary(hh) ?
          mesh->face_handle(mesh->opposite_halfedge_handle(hh)) :
          mesh->face_handle(hh);
      }
      else if (current_ip.type == IPType::ThroughVertex)
        enter_by_vertex(current_ip, current_fh);
    }

    // Assign pixels on the ray to faces.
    int current_x = round(current_ip.point.x());
    int next_x = round(next_ip.point.x());
    // Assign pixels on current face to it
    for (int x = current_x;x < next_x;x++)
    {
      if (assigned_mask(y, x) == 0)
      {
        mesh->data(current_fh).pixels.emplace_back(x, y);
        assigned_mask(y, x) = 1;
      }
    }
    // Enter next face.
    if (next_ip.type == IPType::ThroughEdge)
    {
      if (mesh->is_boundary(next_ip.eh))
        exit_boundary = true;
      else if (!enter_by_edge(next_ip, current_fh, next_fh))
      {
        //Logger::dev_logger->warn("fail to enter by edge.");
        break;
      }
    }
    else if (next_ip.type == IPType::ThroughVertex)
    {
      if (mesh->is_boundary(next_ip.vh))
        exit_boundary = true;
      else if (!enter_by_vertex(next_ip, next_fh))
      {
        //Logger::dev_logger->warn("fail to enter by vertex.");
        break;
      }
    }

    if (exit_boundary)
    {
      current_fh.invalidate();
      next_fh.invalidate();
    }
    else current_fh = next_fh;
  }
}

bool CAssignPixelsToFace::enter_by_edge(
  const InterPoint& inter_point,
  FaceHandle current_fh, FaceHandle& next_fh)
{
  typedef InterPoint::Type IPType;
  EdgeHandle enter_eh = inter_point.eh;
  HalfedgeHandle hh = mesh->halfedge_handle(enter_eh, 0);
  if (mesh->face_handle(hh) == current_fh)
  {
    next_fh = mesh->face_handle(mesh->opposite_halfedge_handle(hh));
    return true;
  }
  if (mesh->face_handle(mesh->opposite_halfedge_handle(hh)) == current_fh)
  {
    next_fh = mesh->face_handle(hh);
    return true;
  }
  //Logger::dev_logger->warn("fail enter edge.");
  throw std::logic_error("fail enter edge.");
  return false;
}

bool CAssignPixelsToFace::enter_by_vertex(
  const InterPoint& inter_point, FaceHandle& next_fh)
{
  VertexHandle vh = inter_point.vh;
  for (HalfedgeHandle voh : mesh->voh_range(vh))
  {
    if (mesh->is_boundary(voh))
      continue;

    HalfedgeHandle prev_hh = mesh->prev_halfedge_handle(voh);
    HalfedgeHandle opp_prev_hh = mesh->opposite_halfedge_handle(prev_hh);
    Vec3d first_tv = mesh->tangent_vec(voh);
    Vec3d second_tv = mesh->tangent_vec(opp_prev_hh);
    if (first_tv.y() < 0.0 && second_tv.y() > 0.0)
    {
      next_fh = mesh->face_handle(voh);
      return true;
    }
    else if (first_tv.y() == 0.0)
    {
      double cur_y = mesh->point(vh).y();
      double next_y = mesh->point(mesh->to_vertex_handle(voh)).y();
      if (next_y <= cur_y)
        next_fh = mesh->face_handle(voh);
      else if (next_y > cur_y)
        next_fh = mesh->face_handle(mesh->opposite_halfedge_handle(voh));
      return true;
    }
  }
  //Logger::dev_logger->warn("fail enter vertex.");
  throw std::logic_error("fail enter vertex.");
  return false;
}

void CAssignPixelsToFace::restart_enter_at_edge(
  const InterPoint& inter_point, FaceHandle& next_fh)
{
  // TODO: restart by using geometry method.
  Logger::dev_logger->warn("Havn't implement restart entering at an edge.");
  Logger::dev_logger->warn("Stop traversing faces along the ray.");
}

ImageT::Color CAssignPixelsToFace::random_color()
{
  return ImageT::Color((double)rand() / RAND_MAX, (double)rand() / RAND_MAX, (double)rand() / RAND_MAX);
}

}// namespace ImageTriSimp
}// namespace GCLF