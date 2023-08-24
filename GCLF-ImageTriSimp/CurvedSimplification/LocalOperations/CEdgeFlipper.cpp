#pragma once
#include "CEdgeFlipper.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CEdgeFlipper::CEdgeFlipper(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{}

bool CEdgeFlipper::init(EdgeHandle _e, size_t Np)
{
  clear();
  // can't flip edge
  if (mesh->is_boundary(_e))
    return false;
  if (!mesh->is_flip_ok(_e))
    return false;

  // Initialize handles
  flip_eh = _e;

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

void CEdgeFlipper::clear()
{
  CLocalOperation::clear();
  flip_eh.invalidate();
}

bool CEdgeFlipper::flip_would_decrease_valence()
{
  int target_val_a0 = mesh->is_boundary(ghs.va0) ? 4 : 6;
  int target_val_a1 = mesh->is_boundary(ghs.va1) ? 4 : 6;
  int target_val_b0 = mesh->is_boundary(ghs.vb0) ? 4 : 6;
  int target_val_b1 = mesh->is_boundary(ghs.vb1) ? 4 : 6;

  int deviation_pre =
    std::abs((int)(mesh->valence(ghs.va0)) - target_val_a0) +
    std::abs((int)(mesh->valence(ghs.va1)) - target_val_a1) +
    std::abs((int)(mesh->valence(ghs.vb0)) - target_val_b0) +
    std::abs((int)(mesh->valence(ghs.vb1)) - target_val_b1);
  int deviation_post =
    std::abs((int)(mesh->valence(ghs.va0)) - 1 - target_val_a0) +
    std::abs((int)(mesh->valence(ghs.va1)) + 1 - target_val_a1) +
    std::abs((int)(mesh->valence(ghs.vb0)) - 1 - target_val_b0) +
    std::abs((int)(mesh->valence(ghs.vb1)) + 1 - target_val_b1);

  return deviation_post < deviation_pre;
}

bool CEdgeFlipper::flip_would_exceed_max_valence(uint32_t max_valence)
{
  return mesh->valence(ghs.va1) + 1 > max_valence ||
    mesh->valence(ghs.vb1) + 1 > max_valence;
}

std::vector<Points> CEdgeFlipper::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Points> population;

  // set the control point to the mid-point of two end-points
  const Vec3d& pa1 = mesh->point(ghs.va1);
  const Vec3d& pb1 = mesh->point(ghs.vb1);
  Vec3d center = (pa1 + pb1) * 0.5;
  population.reserve(Np);
  population.resize(Np - 1, { center });

  // random disturb the control point.
  BoundingBox bbox = bounding_box();
  double radius = (bbox.max() - bbox.min()).length() * 0.1;
  for (Points& ps : population)
  {
    for (Vec3d& p : ps)
    {
      p.x() += (init_dis(init_generator) - 0.5) * radius * 2.0;
      p.y() += (init_dis(init_generator) - 0.5) * radius * 2.0;
    }
  }
  population.push_back({ center });

  return population;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void CEdgeFlipper::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();

  ghs.a0 = mesh->halfedge_handle(flip_eh, 0);
  ghs.b0 = mesh->halfedge_handle(flip_eh, 1);
  ghs.a1 = mesh->next_halfedge_handle(ghs.a0);
  ghs.a2 = mesh->next_halfedge_handle(ghs.a1);
  ghs.b1 = mesh->next_halfedge_handle(ghs.b0);
  ghs.b2 = mesh->next_halfedge_handle(ghs.b1);
  ghs.va0 = mesh->to_vertex_handle(ghs.a0);
  ghs.va1 = mesh->to_vertex_handle(ghs.a1);
  ghs.vb0 = mesh->to_vertex_handle(ghs.b0);
  ghs.vb1 = mesh->to_vertex_handle(ghs.b1);
  ghs.fa = mesh->face_handle(ghs.a0);
  ghs.fb = mesh->face_handle(ghs.b0);

  global_affected_faces.insert(ghs.fa);
  global_affected_faces.insert(ghs.fb);
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void CEdgeFlipper::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Intialize local handles
  lhs.a0 = local_mesh->find_halfedge(v2v.at(ghs.vb0), v2v.at(ghs.va0));
  lhs.b0 = local_mesh->opposite_halfedge_handle(lhs.a0);
  lhs.a1 = local_mesh->next_halfedge_handle(lhs.a0);
  lhs.a2 = local_mesh->next_halfedge_handle(lhs.a1);
  lhs.b1 = local_mesh->next_halfedge_handle(lhs.b0);
  lhs.b2 = local_mesh->next_halfedge_handle(lhs.b1);
  lhs.va0 = local_mesh->to_vertex_handle(lhs.a0);
  lhs.va1 = local_mesh->to_vertex_handle(lhs.a1);
  lhs.vb0 = local_mesh->to_vertex_handle(lhs.b0);
  lhs.vb1 = local_mesh->to_vertex_handle(lhs.b1);
  lhs.fa = local_mesh->face_handle(lhs.a0);
  lhs.fb = local_mesh->face_handle(lhs.b0);

  // Flip the local mesh
  EdgeHandle local_flip_eh = local_mesh->edge_handle(lhs.a0);
  local_mesh->flip_2(local_flip_eh);

  // Find control points to optimize.
  local_optimizing_cpi.push_back(local_mesh->data(local_flip_eh).cpi);
  local_affected_edges.push_back(local_flip_eh);
}

void CEdgeFlipper::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  // do real flip operation.
  mesh->flip_2(flip_eh);
  // update control points in global mesh.
  mesh->ctrl_pnt[mesh->data(flip_eh).cpi] = new_ctrlpnts[0];
  for (auto [cpi, cp] : mesh->data(ghs.fa).bface.global_edge_interior_ctrlpnt_ms(ghs.vb1.idx(), ghs.va1.idx()))
    cp = mesh->ctrl_pnt[cpi];
  for (auto [cpi, cp] : mesh->data(ghs.fb).bface.global_edge_interior_ctrlpnt_ms(ghs.va1.idx(), ghs.vb1.idx()))
    cp = mesh->ctrl_pnt[cpi];
  // update pixels, color and error in global mesh.
  mesh->data(ghs.fa).pixels = local_mesh->data(lhs.fa).pixels;
  mesh->data(ghs.fb).pixels = local_mesh->data(lhs.fb).pixels;
  mesh->data(ghs.fa).color = local_mesh->data(lhs.fa).color;
  mesh->data(ghs.fb).color = local_mesh->data(lhs.fb).color;
  mesh->data(ghs.fa).error = local_mesh->data(lhs.fa).error;
  mesh->data(ghs.fb).error = local_mesh->data(lhs.fb).error;
}

}// namespace ImageTriSimp
}// namespace GCLF