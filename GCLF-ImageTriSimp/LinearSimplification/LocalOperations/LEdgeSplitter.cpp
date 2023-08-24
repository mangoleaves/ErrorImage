#include "LEdgeSplitter.h"
#include "LinearSimplification/LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{
LEdgeSplitter::LEdgeSplitter(LinearSimplifier* _simplifier)
  :LLocalOperation(_simplifier)
{}

/// @brief Initialize neccessary data before splitting.
bool LEdgeSplitter::init(EdgeHandle _e, size_t Np)
{
  clear();
  split_e = _e;

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

void LEdgeSplitter::clear()
{
  LLocalOperation::clear();
  initialized = false;
}

std::vector<Vec3d> LEdgeSplitter::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Vec3d> population; population.reserve(Np);
  Vec3d& from_p = mesh->point(mesh->from_vertex_handle(h0));
  Vec3d& to_p = mesh->point(mesh->to_vertex_handle(h0));
  Vec3d step = (to_p - from_p) / (Np + 1);
  Vec3d start = from_p;
  for (size_t i = 0;i < Np;i++)
  {
    start += step;
    population.push_back(start);
  }
  return population;
}

bool LEdgeSplitter::split_would_cause_over_valence(uint32_t max_valence)
{
  ASSERT(initialized, "splitter not initialized.");

  return mesh->valence(mesh->opposite_vh(h0)) < max_valence &&
    mesh->valence(mesh->opposite_vh(o0)) < max_valence;
}

/// @brief Find affected faces on mesh and affected pixels on image.
void LEdgeSplitter::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  // Initialize handles
  h0 = mesh->halfedge_handle(split_e, 0);
  o0 = mesh->halfedge_handle(split_e, 1);
  // Initialize pixels
  if (!mesh->is_boundary(h0))
  {
    global_affected_faces.insert(mesh->face_handle(h0));
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(mesh->face_handle(h0)).pixels.begin(), mesh->data(mesh->face_handle(h0)).pixels.end());
    // h1 = mesh->next_halfedge_handle(h0); 
    // h2 = mesh->prev_halfedge_handle(h0); 
  }
  if (!mesh->is_boundary(o0))
  {
    global_affected_faces.insert(mesh->face_handle(o0));
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(mesh->face_handle(o0)).pixels.begin(), mesh->data(mesh->face_handle(o0)).pixels.end());
    // o1 = mesh->next_halfedge_handle(o0);
    // o2 = mesh->prev_halfedge_handle(o0);
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void LEdgeSplitter::initialize_local_mesh()
{
  local_mesh = std::make_unique<LMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Split the local mesh
  VertexHandle local_to_v = v2v.at(mesh->to_vertex_handle(h0));
  VertexHandle local_from_v = v2v.at(mesh->to_vertex_handle(o0));
  local_split_e = local_mesh->edge_handle(local_mesh->find_halfedge(local_from_v, local_to_v));
  local_split_v = local_mesh->split(local_split_e, local_mesh->calc_centroid(local_split_e));
  local_optimizing_vertex = local_split_v;

  do_project_x = false;
  do_project_y = false;
  if (mesh->is_boundary(split_e))
  {
    Vec3d& to_p = local_mesh->point(local_to_v);
    Vec3d& from_p = local_mesh->point(local_from_v);

    if (fabs(from_p.y() - to_p.y()) < fabs(from_p.x() - to_p.x()))
    {
      // project to y
      do_project_y = true;
      y_project_to = to_p.y();
    }
    else
    {
      // project to x
      do_project_x = true;
      x_project_to = to_p.x();
    }
  }
}

void LEdgeSplitter::update_local_mesh(LMeshT& new_local_mesh, const Vec3d& new_point)
{
  new_local_mesh.set_point(local_optimizing_vertex, new_point);
  if (do_project_x)
    new_local_mesh.point(local_optimizing_vertex).x() = x_project_to;
  else if (do_project_y)
    new_local_mesh.point(local_optimizing_vertex).y() = y_project_to;
}

void LEdgeSplitter::update_global_mesh(const Vec3d& new_point)
{
  // update global mesh
  // do real split operation.
  split_v = mesh->split(split_e, local_mesh->point(local_split_v));
  // update pixels, color and error in global mesh.
  for (HalfedgeHandle voh : local_mesh->voh_range(local_split_v))
  {
    if (local_mesh->is_boundary(voh))
      continue;
    VertexHandle local_to_vh = local_mesh->to_vertex_handle(voh);
    FaceHandle local_fh = local_mesh->face_handle(voh);
    FaceHandle global_fh = mesh->face_handle(mesh->find_halfedge(split_v, rv2v[local_to_vh.idx()]));
    mesh->data(global_fh).pixels = local_mesh->data(local_fh).pixels;
    mesh->data(global_fh).color = local_mesh->data(local_fh).color;
    mesh->data(global_fh).error = local_mesh->data(local_fh).error;
  }
}

void LEdgeSplitter::post_update()
{
  // update error on edges
#ifdef ONLY_SPLIT_HIGH_ERROR_EDGE
  for (HalfedgeHandle voh : mesh->voh_range(split_v))
  {
    HalfedgeHandle opp_hh = mesh->opposite_halfedge_handle(voh);
    HalfedgeHandle next_hh = mesh->next_halfedge_handle(voh);
    HalfedgeHandle opp_next_hh = mesh->opposite_halfedge_handle(next_hh);
    EdgeHandle out_eh = mesh->edge_handle(voh);
    EdgeHandle next_eh = mesh->edge_handle(next_hh);
    mesh->data(out_eh).error = 0.0;
    mesh->data(next_eh).error = 0.0;
    if (!mesh->is_boundary(voh))
    {
      double e = mesh->data(mesh->face_handle(voh)).error;
      mesh->data(out_eh).error += e;
      mesh->data(next_eh).error += e;
    }
    if (!mesh->is_boundary(opp_hh))
      mesh->data(out_eh).error += mesh->data(mesh->face_handle(opp_hh)).error;
    if (!mesh->is_boundary(opp_next_hh))
      mesh->data(next_eh).error += mesh->data(mesh->face_handle(opp_next_hh)).error;
}
#endif
}

}// namespace ImageTriSimp
}// namespace GCLF