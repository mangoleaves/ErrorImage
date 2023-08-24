#include "CEdgeCollapser.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

CEdgeCollapser::CEdgeCollapser(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{}

bool CEdgeCollapser::init(HalfedgeHandle _hh, size_t Np)
{
  clear();
  // Initialize handles
  if (mesh->is_boundary(mesh->edge_handle(_hh)))
    return false;
  collapse_he = _hh;
  collapse_he_opp = mesh->opposite_halfedge_handle(_hh);
  from_vh = mesh->from_vertex_handle(collapse_he);
  to_vh = mesh->to_vertex_handle(collapse_he);
  if (mesh->is_boundary(from_vh) || mesh->is_boundary(to_vh))
    return false;
  if (!mesh->is_collapse_ok(collapse_he))
    return false;

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

void CEdgeCollapser::clear()
{
  CLocalOperation::clear();
  collapse_he.invalidate();
  collapse_he_opp.invalidate();
  from_vh.invalidate();
  to_vh.invalidate();
}

uint32_t CEdgeCollapser::valence_after_collapse()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return local_mesh->valence(v2v.at(to_vh));
}

VertexHandle CEdgeCollapser::collapse_vertex()
{
  ASSERT(initialized, "Collapser it not initialized.");
  return to_vh;
}

Vec3d CEdgeCollapser::find_smooth_target(const Vec3d& original_point)const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  double weight_sum = 0.0;
  for (VertexHandle vv : local_mesh->vv_range(v2v.at(to_vh)))
  {
    double weight = (local_mesh->point(vv) - original_point).length();
    q += local_mesh->point(vv) * weight;
    weight_sum += weight;
  }
  q /= weight_sum;
  return q;
}

std::vector<Points> CEdgeCollapser::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Points> population; population.resize(Np);
  auto first = population.begin();

  // set the control points to the original control points.
  Points variables = get_variables();
  first = std::fill_n(first, Np / 2, variables);  // edge is collapsed to to_vh.

  // set the control points by using Laplacian smoothing.
  VertexHandle local_to_vh = v2v.at(to_vh);
  {
    Points individual;  individual.reserve(variables.size());
    Vec3d smoothed_to_p = find_smooth_target(mesh->point(to_vh));
    individual.push_back(smoothed_to_p);
    for (HalfedgeHandle hh : local_mesh->voh_range(local_to_vh))
      individual.push_back((smoothed_to_p + local_mesh->point(local_mesh->to_vertex_handle(hh))) * 0.5);
    std::fill(first, population.end(), individual);
  }

  // random disturb control points
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
void CEdgeCollapser::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (HalfedgeHandle voh : mesh->voh_range(from_vh))
  {
    FaceHandle fh = mesh->face_handle(voh);
    global_affected_faces.insert(fh);
  }
  for (HalfedgeHandle voh : mesh->voh_range(to_vh))
  {
    FaceHandle fh = mesh->face_handle(voh);
    global_affected_faces.insert(fh);
  }
  for (FaceHandle fh : global_affected_faces)
  {
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
  n_faces_after_collapsing = global_affected_faces.size();
  if (!mesh->is_boundary(collapse_he)) n_faces_after_collapsing--;
  if (!mesh->is_boundary(collapse_he_opp)) n_faces_after_collapsing--;
}

void CEdgeCollapser::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Collapse the local mesh
  HalfedgeHandle local_collapse_he = local_mesh->find_halfedge(v2v.at(from_vh), v2v.at(to_vh));
  VertexHandle local_to_vh = local_mesh->to_vertex_handle(local_collapse_he);
  local_mesh->collapse(local_collapse_he);  // the edge is collapsed to v2v.at(to_vh).

  // Find control points to optimize.
  local_optimizing_cpi.push_back(local_mesh->data(local_to_vh).cpi);
  for (HalfedgeHandle hh : local_mesh->voh_range(local_to_vh))
  {
    EdgeHandle eh = local_mesh->edge_handle(hh);
    local_optimizing_cpi.push_back(local_mesh->data(eh).cpi);
    local_affected_edges.push_back(eh);
  }
}

void CEdgeCollapser::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  // do real collapse operation.
  mesh->collapse(collapse_he);
  // update control points in global mesh.
  VertexHandle local_vh = v2v.at(to_vh);
  for (FaceHandle local_fh : local_mesh->vf_range(local_vh))
  {
    FaceHandle global_fh = rf2f[local_fh.idx()];
    // update control points
    mesh->data(global_fh).bface.ctrl_pnt = local_mesh->data(local_fh).bface.ctrl_pnt;
    for (auto [cp_idx, cp] : mesh->data(global_fh).bface.global_ctrlpnt_cs())
      mesh->ctrl_pnt[cp_idx] = cp;
    // update pixels, color and error
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF