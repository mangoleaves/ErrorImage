#include "CVertexRelocater.h"
#include "CurvedSimplification/CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

#define OPT_ONLY_VERTEX

CVertexRelocater::CVertexRelocater(CurvedSimplifier* _simplifier)
  :CLocalOperation(_simplifier)
{}

bool CVertexRelocater::init(VertexHandle _v, size_t Np)
{
  clear();
  // Initialize handles
  relocate_vh = _v;

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

void CVertexRelocater::clear()
{
  CLocalOperation::clear();
  relocate_vh.invalidate();
}


std::vector<Points> CVertexRelocater::get_init_population(
  std::default_random_engine& init_generator,
  std::uniform_real_distribution<double>& init_dis,
  size_t Np)const
{
  std::vector<Points> population; population.resize(Np);
  auto first = population.begin();

  // set the control points to the original control points.
  Points variables = get_variables();
  std::fill(first, population.end(), variables);  // edge is collapsed to to_vh.

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
void CVertexRelocater::find_affected_faces_pixels()
{
  simplifier->image->reset_mask(pixel_assigned);
  pixel_assigned.setOnes();
  for (HalfedgeHandle voh : mesh->voh_range(relocate_vh))
  {
    FaceHandle fh = mesh->face_handle(voh);
    if (mesh->is_boundary(voh))
      continue;
    global_affected_faces.insert(fh);
    affected_pixel_coords.insert(affected_pixel_coords.end(),
      mesh->data(fh).pixels.begin(), mesh->data(fh).pixels.end());
  }
  for (const Vec2i& pixel : affected_pixel_coords)
    pixel_assigned(pixel.y(), pixel.x()) = 0;
}

void CVertexRelocater::initialize_local_mesh()
{
  local_mesh = std::make_unique<BMeshT>();
  mesh->local_mesh(global_affected_faces, *local_mesh, v2v, f2f, rv2v, rf2f);

  // Find control points to optimize.
  // Find affected local edges.
  VertexHandle local_relocate_vh = v2v.at(relocate_vh);
  local_optimizing_cpi.push_back(local_mesh->data(local_relocate_vh).cpi);
  for (EdgeHandle eh : local_mesh->ve_range(local_relocate_vh))
  {
  #ifndef OPT_ONLY_VERTEX
    local_optimizing_cpi.push_back(local_mesh->data(eh).cpi);
  #endif
    local_affected_edges.push_back(eh);
}
}

void CVertexRelocater::update_global_mesh(const Points& new_ctrlpnts)
{
  // update global mesh
  VertexHandle local_vh = v2v.at(relocate_vh);
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

#undef OPT_ONLY_VERTEX
}// namespace ImageTriSimp
}// namespace GCLF