#include "CLocalOperation.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CLocalOperation::CLocalOperation(CurvedSimplifier* _simplifier)
{
  simplifier = _simplifier;
  mesh = _simplifier->mesh.get();
  image = _simplifier->image.get();
  initialized = false;
}

void CLocalOperation::clear()
{
  initialized = false;
  global_affected_faces.clear();
  affected_pixel_coords.clear();
  local_mesh = nullptr;
  v2v.clear();
  f2f.clear();
  rv2v.clear();
  rf2f.clear();
  local_optimizing_cpi.clear();
  local_affected_edges.clear();
  vus.clear();
}

BoundingBox CLocalOperation::bounding_box()const
{
  ASSERT(initialized, "Operator is not initialized.");

  return BoundingBox(local_mesh->ctrl_pnt);
}

Points CLocalOperation::get_variables()const
{
  ASSERT(initialized, "Operator is not initialized.");

  Points ctrlpnts;  ctrlpnts.reserve(local_optimizing_cpi.size());
  for (Index cpi : local_optimizing_cpi)
    ctrlpnts.push_back(local_mesh->ctrl_pnt[cpi]);
  return ctrlpnts;
}

double CLocalOperation::error_before()const
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

double CLocalOperation::quality_before()const
{
#ifdef MIPS_SUM
  return std::accumulate(
    global_affected_faces.begin(), global_affected_faces.end(), 0.0,
    [&](double v, FaceHandle fh) {return v + MIPS_energy(mesh->data(fh).bface.ctrl_pnt, mesh->data(fh).bface.degree);});
#endif
#ifdef MIPS_MAX
  double max_quality = 0.0;
  for (FaceHandle fh : global_affected_faces)
  {
    double face_quality = MIPS_energy(mesh->data(fh).bface.ctrl_pnt, mesh->data(fh).bface.degree);
    if (face_quality > max_quality)
      max_quality = face_quality;
  }
  return max_quality;
#endif
#ifdef APPROX_CURVATURE_SUM
  return std::accumulate(
    local_affected_edges.begin(), local_affected_edges.end(), 0.0,
    [&](double v, EdgeHandle eh) {return v + approx_bezier_curvature(local_mesh.get(), eh);});
#endif
#ifdef APPROX_CURVATURE_MAX
  double max_quality = -1.0;
  for (EdgeHandle eh : local_affected_edges)
  {
    double edge_quality = approx_bezier_curvature(local_mesh.get(), eh);
    if (edge_quality > max_quality)
      max_quality = edge_quality;
  }
  return max_quality;
#endif
}

bool CLocalOperation::update_variables(const Points& new_ctrlpnts, size_t pidx)
{
  ASSERT(initialized, "Collapser is not initialized.");
  ASSERT(local_optimizing_cpi.size() == new_ctrlpnts.size(), "sizes of variables are not same.");

  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  vt_updated = false;
  // there is no pixel to assign, exit.
  if (affected_pixel_coords.empty())
    return false;

  update_local_ctrlpnts(vt_local_mesh, new_ctrlpnts);

  // checks
  if (!check_constraints(vt_local_mesh))
    return false;
  // possibly more checks

  vt_updated = true;
  return true;
}

bool CLocalOperation::check_constraints(const BMeshT& new_local_mesh)const
{
  if (operation_would_cause_flip(new_local_mesh))
    return false;
  return true;
}

double CLocalOperation::error_after(size_t pidx)
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

double CLocalOperation::quality_after(size_t pidx)
{
  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  ASSERT(vt_updated, "havn't updated.");

#ifdef MIPS_SUM
  return std::accumulate(
    vt_local_mesh.faces_begin(), vt_local_mesh.faces_end(), 0.0,
    [&](double v, FaceHandle fh) {return v + MIPS_energy(vt_local_mesh.data(fh).bface.ctrl_pnt, vt_local_mesh.data(fh).bface.degree);});
#endif
#ifdef MIPS_MAX
  double max_quality = 0.0;
  for (FaceHandle fh : vt_local_mesh.faces())
  {
    if (vt_local_mesh.status(fh).deleted())
      continue;
    double face_quality = MIPS_energy(vt_local_mesh.data(fh).bface.ctrl_pnt, vt_local_mesh.data(fh).bface.degree);
    if (face_quality > max_quality)
      max_quality = face_quality;
  }
  return max_quality;
#endif
#ifdef APPROX_CURVATURE_SUM
  return std::accumulate(
    local_affected_edges.begin(), local_affected_edges.end(), 0.0,
    [&](double v, EdgeHandle eh) {return v + approx_bezier_curvature(&vt_local_mesh, eh);});
#endif
#ifdef APPROX_CURVATURE_MAX
  double max_quality = -1.0;
  for (EdgeHandle eh : local_affected_edges)
  {
    double edge_quality = approx_bezier_curvature(&vt_local_mesh, eh);
    if (edge_quality > max_quality)
      max_quality = edge_quality;
  }
  return max_quality;
#endif
}


bool CLocalOperation::is_quality_bounded(size_t pidx, double curvature_bound, double min_angle_bound)
{
  auto& vt_local_mesh = vus[pidx].vt_local_mesh;
  auto& vt_pixel_assigned = vus[pidx].vt_pixel_assigned;
  auto& vt_updated = vus[pidx].vt_updated;

  ASSERT(vt_updated, "havn't updated.");

  for (EdgeHandle eh : local_affected_edges)
  {
    if (approx_bezier_curvature(&vt_local_mesh, eh) > curvature_bound)
      return false;
  }

  for (EdgeHandle eh : local_affected_edges)
  {
    HalfedgeHandle hh = vt_local_mesh.halfedge_handle(eh, 0);
    if (vt_local_mesh.is_boundary(hh))
      continue;

    if (approx_bezier_sector_angle(&vt_local_mesh, hh) > min_angle_bound)
      return false;
    if (approx_bezier_sector_angle(&vt_local_mesh, vt_local_mesh.next_halfedge_handle(hh)) > min_angle_bound)
      return false;
    if (approx_bezier_sector_angle(&vt_local_mesh, vt_local_mesh.prev_halfedge_handle(hh)) > min_angle_bound)
      return false;
  }
  return true;
}

void CLocalOperation::just_do_it(const Points& new_ctrlpnts)
{
  update_local_ctrlpnts(*local_mesh, new_ctrlpnts);
  update_pixel_color_error(*local_mesh);

  update_global_mesh(new_ctrlpnts);
}

void CLocalOperation::update_local_ctrlpnts(BMeshT& new_local_mesh, const Points& new_ctrlpnts)
{
  // update control points in new local mesh.
  for (size_t i = 0;i < local_optimizing_cpi.size();i++)
    new_local_mesh.ctrl_pnt[local_optimizing_cpi[i]] = new_ctrlpnts[i];
  new_local_mesh.sync_local_ctrlpnt();
  // extra cpi is not stored in VertexData and EdgeData.
}

void CLocalOperation::update_pixel_color_error(BMeshT& local_mesh)
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

bool CLocalOperation::operation_would_cause_flip(const BMeshT& new_local_mesh)const
{
  return !check_injectivity(&new_local_mesh);
}

bool CLocalOperation::is_corner_vertex(VertexHandle vh)
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