#include "LLocalOperation.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

LLocalOperation::LLocalOperation(LinearSimplifier* _simplifier)
{
  simplifier = _simplifier;
  mesh = _simplifier->mesh.get();
  image = _simplifier->image.get();
  initialized = false;
}

void LLocalOperation::clear()
{
  initialized = false;
  global_affected_faces.clear();
  affected_pixel_coords.clear();
  local_mesh = nullptr;
  v2v.clear();
  f2f.clear();
  rv2v.clear();
  rf2f.clear();
  local_optimizing_vertex.invalidate();
  vus.clear();
}

BoundingBox LLocalOperation::bounding_box()const
{
  ASSERT(initialized, "Operator is not initialized.");

  BoundingBox box;
  for (VertexHandle vh : local_mesh->vertices())
    box += local_mesh->point(vh);
  return box;
}

Vec3d LLocalOperation::get_variable()const
{
  ASSERT(initialized, "Operator is not initialized.");

  return local_mesh->point(local_optimizing_vertex);
}

std::vector<Vec3d> LLocalOperation::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Vec3d> population;
  BoundingBox bbox = bounding_box();
  double radius = (bbox.max() - bbox.min()).length() * 0.1;
  population.resize(Np, get_variable());
  for (Vec3d& p : population)
  {
    // TODO: consider nearby feature points
    p.x() += (init_dis(init_generator) - 0.5) * radius * 2.0;
    p.y() += (init_dis(init_generator) - 0.5) * radius * 2.0;
  }
  return population;
}

double LLocalOperation::error_before()const
{
  ASSERT(initialized, "Relocater is not initialized.");

  double color_error = 0.0;
  for (FaceHandle fh : global_affected_faces)
  {
    if (!mesh->data(fh).pixels.empty())
      color_error += mesh->data(fh).error;
  }
  return color_error;
}

double LLocalOperation::quality_before()const
{
  double max_quality = -1.0;
  for (FaceHandle fh : global_affected_faces)
  {
    double face_quality = calc_max_cos(mesh, fh);
    if (face_quality > max_quality)
      max_quality = face_quality;
  }
  return max_quality;
}

bool LLocalOperation::update_variable(const Vec3d& new_point, size_t pidx)
{
  ASSERT(initialized, "Collapser is not initialized.");

  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  vt_updated = false;
  // there is no pixel to assign, exit.
  if (affected_pixel_coords.empty())
    return false;

  update_local_mesh(vt_local_mesh, new_point);

  // checks
  if (!check_constraints(vt_local_mesh))
    return false;
  // possibly more checks

  vt_updated = true;
  return true;
}


bool LLocalOperation::check_constraints(const LMeshT& new_local_mesh)const
{
  if (operation_would_cause_flip(new_local_mesh))
    return false;
  return true;
}

double LLocalOperation::error_after(size_t pidx)
{
  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  ASSERT(vt_updated, "havn't updated.");

  auto& behavior = simplifier->behavior;

  // assign pixels to new local mesh.
  try
  {
    auto assigner = behavior.assign_pixels_to_face;
    assigner(&vt_local_mesh, image, vt_pixel_assigned, affected_pixel_coords);
  }
  catch (...)
  {
    return DBL_MAX;
  }

  // calculate color and error on new local mesh.
  for (FaceHandle fh : vt_local_mesh.faces())
  {
    if (vt_local_mesh.status(fh).deleted())
      continue;
    vt_local_mesh.data(fh).color = behavior.calc_color_on_face(
      vt_local_mesh.data(fh).pixels);
    vt_local_mesh.data(fh).error = behavior.calc_error_on_face(
      vt_local_mesh.data(fh).color,
      vt_local_mesh.data(fh).pixels);
  }

  double color_error = std::accumulate(
    vt_local_mesh.faces_begin(), vt_local_mesh.faces_end(), 0.0,
    [&](double v, FaceHandle fh) {return v + vt_local_mesh.data(fh).error;});
  return color_error;
}

double LLocalOperation::quality_after(size_t pidx)
{
  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  ASSERT(vt_updated, "havn't updated.");

  return calc_max_cos(&vt_local_mesh);
}

void LLocalOperation::just_do_it(const Vec3d& new_point)
{
  update_local_mesh(*local_mesh, new_point);
  update_pixel_color_error(*local_mesh);

  update_global_mesh(new_point);

  post_update();
}

void LLocalOperation::update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point)
{
  new_local_mesh.set_point(local_optimizing_vertex, new_point);
}

void LLocalOperation::update_pixel_color_error(LMeshT& local_mesh)
{
  // assign pixels to local mesh.
  auto& behavior = simplifier->behavior;
  behavior.assign_pixels_to_face(
    &local_mesh, image, pixel_assigned, affected_pixel_coords);

  // calculate color and error on local mesh.
  for (FaceHandle fh : local_mesh.faces())
  {
    if (local_mesh.status(fh).deleted())
      continue;
    local_mesh.data(fh).color = behavior.calc_color_on_face(
      local_mesh.data(fh).pixels);
    local_mesh.data(fh).error = behavior.calc_error_on_face(
      local_mesh.data(fh).color,
      local_mesh.data(fh).pixels);
  }
}

bool LLocalOperation::operation_would_cause_flip(const LMeshT& new_local_mesh)const
{
  return check_flip(&new_local_mesh);
}

bool LLocalOperation::is_corner_vertex(VertexHandle vh)
{
  ASSERT(mesh->is_boundary(vh), "vh is not boundary vertex.");
  HalfedgeHandle hh = mesh->halfedge_handle(vh);
  ASSERT(mesh->from_vertex_handle(hh) == vh, "halfedge's from vertex is wrong.");

  const Vec3d& cur_p = mesh->point(vh);
  const Vec3d& next_p = mesh->point(mesh->to_vertex_handle(hh));
  const Vec3d& prev_p = mesh->point(mesh->from_vertex_handle(mesh->prev_halfedge_handle(hh)));
  return !(fabs(cur_p.x() - next_p.x()) < 1e-6 && fabs(cur_p.x() - prev_p.x()) < 1e-6) &&
    !(fabs(cur_p.y() - next_p.y()) < 1e-6 && fabs(cur_p.y() - prev_p.y()) < 1e-6);
}

}// namespace ImageTriSimp
}// namespace GCLF