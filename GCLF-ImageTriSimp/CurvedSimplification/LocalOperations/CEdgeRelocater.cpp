#include "CEdgeRelocater.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CEdgeRelocater::CEdgeRelocater(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{}

bool CEdgeRelocater::init(EdgeHandle _e, size_t Np)
{
  clear();
  if (mesh->is_boundary(_e))
    return false;

  // Initialize handles
  relocate_eh = _e;
  relocate_hh = mesh->halfedge_handle(relocate_eh, 0);
  relocate_hh_opp = mesh->halfedge_handle(relocate_eh, 1);
  from_vh = mesh->from_vertex_handle(relocate_hh);
  to_vh = mesh->to_vertex_handle(relocate_hh);

  // Initialize neccessary topology data.
  find_affected_faces_pixels();

  // Initialize neccessary geometry data.
  initialize_local_mesh();

  // Initialize virtual local mesh used to virtually update.
  virtual_update_struct tmp_vus;
  tmp_vus.vt_local_mesh = *local_mesh;
  tmp_vus.vt_pixel_assigned = pixel_assigned;
  tmp_vus.vt_updated = false;
  vus.resize(Np, tmp_vus);

  initialized = true;
  return true;
}

void CEdgeRelocater::clear()
{
  CLocalOperation::clear();
  relocate_eh.invalidate();
}

std::vector<Points> CEdgeRelocater::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Points> population;
  population.resize(Np);
  auto first = population.begin();

  // set the control point to the original control point.
  first = std::fill_n(first, Np / 2, get_variables());
  // set the contrl point to the mid-point of two end-points.
  const Vec3d& from_p = mesh->point(from_vh);
  const Vec3d& to_p = mesh->point(to_vh);
  Vec3d center = (from_p + to_p) * 0.5;
  std::fill(first, population.end(), Points({ center }));

  // random disturb control point
  BoundingBox bbox = bounding_box();
  double radius = (bbox.max() - bbox.min()).length() * 0.1;
  for (Points& ps : population)
  {
    for (Vec3d& p : ps)
    {
      // TODO: consider nearby feature points
      p.x() += (init_dis(init_generator) - 0.5) * radius * 2.0;
      p.y() += (init_dis(init_generator) - 0.5) * radius * 2.0;
    }
  }
  return population;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void CEdgeRelocater::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();

  FaceHandle fa = mesh->face_handle(relocate_hh);
  FaceHandle fb = mesh->face_handle(relocate_hh_opp);

  global_affected_faces.insert(fa);
  global_affected_faces.insert(fb);
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void CEdgeRelocater::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Find control points to optimize.
  HalfedgeHandle local_relocate_hh = local_mesh->find_halfedge(v2v.at(from_vh), v2v.at(to_vh));
  EdgeHandle local_relocate_eh = local_mesh->edge_handle(local_relocate_hh);
  local_optimizing_cpi.push_back(local_mesh->data(local_relocate_eh).cpi);
  local_affected_edges.push_back(local_relocate_eh);
}

void CEdgeRelocater::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  HalfedgeHandle local_relocate_hh = local_mesh->find_halfedge(v2v.at(from_vh), v2v.at(to_vh));
  FaceHandle local_fa = local_mesh->face_handle(local_relocate_hh);
  FaceHandle local_fb = local_mesh->opposite_face_handle(local_relocate_hh);
  {// update fa
    FaceHandle global_fh = rf2f[local_fa.idx()];
    // update control points
    mesh->data(global_fh).bface.ctrl_pnt = local_mesh->data(local_fa).bface.ctrl_pnt;
    for (auto [cp_idx, cp] : mesh->data(global_fh).bface.global_ctrlpnt_cs())
      mesh->ctrl_pnt[cp_idx] = cp;
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fa).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fa).color;
    mesh->data(global_fh).error = local_mesh->data(local_fa).error;
  }
  {// update fb
    FaceHandle global_fh = rf2f[local_fb.idx()];
    // update control points
    mesh->data(global_fh).bface.ctrl_pnt = local_mesh->data(local_fb).bface.ctrl_pnt;
    // for (auto [cp_idx, cp] : mesh->data(global_fh).bface.global_ctrlpnt_cs())
    //   mesh->ctrl_pnt[cp_idx] = cp;
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fb).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fb).color;
    mesh->data(global_fh).error = local_mesh->data(local_fb).error;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF