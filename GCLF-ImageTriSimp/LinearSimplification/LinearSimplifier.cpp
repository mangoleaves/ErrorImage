#include "LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

void LinearSimplifier::initialize(ImageT& input_image, LMeshT& input_mesh, ParamTriangulator& param)
{
  // get configurations (don't set error bound here!)
  set_color_type(param.color_type);
  set_error_type(param.error_type);
  set_quality_bound(param.param_linear.min_angle, param.param_linear.split_min_angle);
  set_max_valence(param.max_valence);
  config.max_simplify_iter = param.param_linear.max_simplify_iter;
  config.convergence_collapse_number = param.param_linear.convergence_collapse_number;
  // initialize DE optimizer.
  optimizer.initialize(this, param);
  // initialize mesh and image.
  image = std::make_unique<ImageT>(input_image);
  mesh = std::make_unique<LMeshT>(input_mesh);
  // scale mesh to cover full image.
  mesh->scale_mesh(*image);
  // initialize color and error
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  // Set the bounding error
  current_error = lp_error();
  // TEST initialize file output
  operation_cnt = 0;
  error_to_file = false;

  Logger::user_logger->info("Initial mesh: v {}, e {}, f {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Initial root_mean_lp_error {}", root_mean_lp_error());
}

void LinearSimplifier::set_color_type(const std::string& color_type)
{
  config.color_type = color_type;
  if (config.color_type == "constant")
  {
    behavior.calc_color_on_face.type =
      LinearSimplifierBehaviors::LCalcColorOnFace::Type::Constant;
  }
  else if (config.color_type == "linear")
  {
    behavior.calc_color_on_face.type =
      LinearSimplifierBehaviors::LCalcColorOnFace::Type::Linear;
  }
  else if (config.color_type == "quadratic")
  {
    behavior.calc_color_on_face.type =
      LinearSimplifierBehaviors::LCalcColorOnFace::Type::Quadratic;
  }
  else
  {
    ASSERT(false, "wrong color type");
  }
}

void LinearSimplifier::set_error_type(const std::string& error_type)
{
  config.error_type = error_type;
  if (config.error_type == "l1")
  {
    behavior.calc_error_on_face.type =
      LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L1Norm;
    behavior.calc_error_on_mesh.p = 1.0;
  }
  else if (config.error_type == "l2")
  {
    behavior.calc_error_on_face.type =
      LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L2Norm;
    behavior.calc_error_on_mesh.p = 2.0;
  }
  else if (config.error_type == "l4")
  {
    behavior.calc_error_on_face.type =
      LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L4Norm;
    behavior.calc_error_on_mesh.p = 4.0;
  }
  else
  {
    ASSERT(false, "wrong error type {}", config.error_type);
  }
}

void LinearSimplifier::set_error_bound(double error_bound)
{
  switch (behavior.calc_error_on_face.type)
  {
  case LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L1Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = error_bound / 255. * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = error_bound * image->width * image->height * 3;
  #endif
    break;
  case LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L2Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = pow(error_bound / 255., 2.) * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = pow(error_bound, 2.) * image->width * image->height * 3;
  #endif
    break;
  case LinearSimplifierBehaviors::LCalcErrorOnFace::Type::L4Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = pow(error_bound / 255., 4.) * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = pow(error_bound, 4.) * image->width * image->height * 3;
  #endif
    break;
  }
  split_to_error_bound();
}

void LinearSimplifier::set_max_valence(uint32_t mv)
{
  config.max_valence = mv;
}

void LinearSimplifier::set_quality_bound(double min_angle, double split_min_angle)
{
  config.quality_bound = std::cos(min_angle / 180. * M_PI);
  config.split_quality_bound = std::cos(split_min_angle / 180. * M_PI);
}

void LinearSimplifier::split_to_error_bound()
{
  size_t max_vertices = image->width * image->height;
  if (current_error <= config.error_bound) return;
  // error_fout = fopen("C:/Everything/MyData/ActiveData/Research/Projects/ImageTriangulationSimplification/Image In Paper/Priority-Queue/data/neither_priority_error.txt", "w");
  // if (!error_fout)
  //  return;

  while (current_error > config.error_bound && mesh->n_vertices() < max_vertices)
  {
    flip(OptStrategy::EOQB);      if (current_error <= config.error_bound) goto SPLIT_END;
    // error_to_file = true;
    for (size_t i = 0;i < 3;i++)
    {
      size_t relocated_verts = relocate_with_priority();
      if (relocated_verts < mesh->n_vertices() * 0.05)
        relocate(OptStrategy::EOQB);
      if (current_error <= config.error_bound) goto SPLIT_END;
    }
    split_edges_with_priority();
    //split_edges(OptStrategy::EOQB);
    //error_to_file = false;
    if (current_error <= config.error_bound) goto SPLIT_END;

    // flip(OptStrategy::EOQB);      if (current_error <= config.error_bound) goto SPLIT_END;
    // relocate_with_priority();     if (current_error <= config.error_bound) goto SPLIT_END;
    // split_faces_with_priority();  if (current_error <= config.error_bound) goto SPLIT_END;
  }

  if (mesh->n_vertices() >= max_vertices)
  {
    // TODO: triangulate image at pixel resolution.
  }

SPLIT_END:
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
}

void LinearSimplifier::simplify()
{
  size_t iter = 0;
  while (iter < config.max_simplify_iter)
  {
    flip(OptStrategy::EOQB);
    for (size_t i = 0;i < 3;i++)
      relocate(OptStrategy::EOQB);
    size_t collapsed_size = collapse_with_priority();

    if (collapsed_size < config.convergence_collapse_number)
    {
      Logger::user_logger->info("convergence.");
      break;
    }
    iter++;
  }

  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
}

/// @brief It's a test function. Emit a ray, find intersection points.
std::vector<Vec3d> LinearSimplifier::intersect_points(const Vec2d& ray_start)
{
  LAssignPixelsToFace assigner;
  Eigen::MatrixXi mask; image->reset_mask(mask);
  assigner.initialize(mesh.get(), image.get(), mask, std::vector<Vec2i>());

  // Use tree to search possible intersected seg.
  auto possible_seg = assigner.tree->intersect_segs(ray_start);

  // Find the intersection points between rays and edges.
  auto inter_points = assigner.find_intersection_points(ray_start, possible_seg);

  // Sort the intersection points on each ray.
  std::sort(inter_points.begin(), inter_points.end());

  // If some intersection points are close to one vertex,
  //  contract the one ring intersection points of the vertex to one intersection point.
  assigner.cluster_and_collapse_inter_points(inter_points);

  // Output intersection points
  std::vector<Vec3d> points;
  points.reserve(inter_points.size());
  for (auto& ip : inter_points)
    points.push_back(ip.point);
  return points;
}

void LinearSimplifier::update_color_error()
{
  // Assign pixels to faces.
  Eigen::MatrixXi mask; image->reset_mask(mask);
  behavior.assign_pixels_to_face(mesh.get(), image.get(), mask, std::vector<Vec2i>());

  // Calculate colors and errors on faces.
  for (FaceHandle fh : mesh->faces())
  {
    mesh->data(fh).color = behavior.calc_color_on_face(fh);
    mesh->data(fh).error = behavior.calc_error_on_face(fh);
  }
}

void LinearSimplifier::initialize_error_on_verts()
{
  for (VertexHandle v : mesh->vertices())
  {
    mesh->data(v).error = 0.;
    for (FaceHandle vf : mesh->vf_range(v))
    {
      mesh->data(v).error += mesh->data(vf).error;
    }
  }
}

void LinearSimplifier::initialize_error_on_edges()
{
  for (EdgeHandle e : mesh->edges())
  {
    mesh->data(e).error = 0.;
    HalfedgeHandle hh0 = mesh->halfedge_handle(e, 0);
    HalfedgeHandle hh1 = mesh->halfedge_handle(e, 1);
    if (!mesh->is_boundary(hh0))
      mesh->data(e).error += mesh->data(mesh->face_handle(hh0)).error;
    if (!mesh->is_boundary(hh1))
      mesh->data(e).error += mesh->data(mesh->face_handle(hh1)).error;
  }
}

size_t LinearSimplifier::relocate(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t relocated_size = 0;
  for (VertexHandle vh : mesh->vertices())
  {
    traversed_size += 1;
    bool do_relocate = false;
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_relocate = relocate_vertex_ebqo(vh);
      break;
    case OptStrategy::EOQB:
      do_relocate = relocate_vertex_eoqb(vh);
      break;
    }
    if (do_relocate)
    {
      relocated_size += 1;
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("relocate {} / {} vertices.", relocated_size, traversed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return relocated_size;
}

size_t LinearSimplifier::collapse(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t collapsed_size = 0;
  for (EdgeHandle eh : mesh->edges())
  {
    if (mesh->status(eh).deleted())
      continue;
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    VertexHandle v0 = mesh->from_vertex_handle(hh), v1 = mesh->to_vertex_handle(hh);
    if (!mesh->is_boundary(eh) && mesh->is_boundary(v0) && mesh->is_boundary(v1))
      continue;

    traversed_size += 1;
    bool do_collapse = false;

    if (!mesh->is_boundary(eh) && !mesh->is_boundary(v0) && !mesh->is_boundary(v1))
    {
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = collapse_edge_ebqo(eh);
        break;
      case OptStrategy::EOQB:
        do_collapse = collapse_edge_eoqb(eh);
        break;
      }
    }
    else if (!mesh->is_boundary(eh))
    {
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = collapse_boundary_vertex_ebqo(eh);
        break;
      case OptStrategy::EOQB:
        do_collapse = collapse_boundary_vertex_eoqb(eh);
        break;
      }
    }
    else
    {
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = collapse_boundary_edge_ebqo(eh);
        break;
      case OptStrategy::EOQB:
        do_collapse = collapse_boundary_edge_eoqb(eh);
        break;
      }
    }

    if (do_collapse)
    {
      collapsed_size += 1;
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
    }
  }
  mesh->garbage_collection();
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("collapse {} / {} edges.", collapsed_size, traversed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return collapsed_size;
}

size_t LinearSimplifier::flip(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t flipped_size = 0;
  for (EdgeHandle eh : mesh->edges())
  {
    traversed_size += 1;
    bool do_flip = false;
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_flip = flip_edge_ebqo(eh);
      break;
    case OptStrategy::EOQB:
      do_flip = flip_edge_eoqb(eh);
      break;
    }
    if (do_flip)
    {
      flipped_size += 1;
      Logger::dev_logger->trace("current error is {:.5f}", current_error);
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("flip {} / {} edges.", flipped_size, traversed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return flipped_size;
}

size_t LinearSimplifier::split_edges(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t splitted_size = 0;
  size_t end_idx = mesh->n_edges();
  for (EdgeHandle eh : mesh->edges())
  {
    if (eh.idx() == end_idx)
      break;
    traversed_size += 1;
    bool do_split = false;
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_split = false;
      break;
    case OptStrategy::EOQB:
    {
      auto splitter = new_edge_splitter();
      // check if flip ok
      if (splitter.init(eh, 1))
      {
        Vec3d opt_pos;
        double error_before, error_after, quality_before, quality_after;
        optimizer.error_optimized_quality_bounded(
          static_cast<LLocalOperation*>(&splitter), config.quality_bound, current_error,
          opt_pos, error_before, error_after, quality_before, quality_after);
        double error_change = error_after - error_before;
        // we must gain larger error change from splitting
        if (error_change < -current_error * 0.0001)
        {
          splitter.just_do_it(opt_pos);
          current_error += error_change;
          Logger::dev_logger->trace("current error is {:.5f}", current_error);
          do_split = true;
        }
        else do_split = false;
      }
      else do_split = false;
    }
    break;
    }
    if (do_split)
    {
      splitted_size += 1;
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
    }
    if (current_error <= config.error_bound)
      break;
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("split {} / {} edges.", splitted_size, traversed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return splitted_size;
}

size_t LinearSimplifier::split_faces(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t splitted_size = 0;
  size_t end_idx = mesh->n_faces();
  for (FaceHandle fh : mesh->faces())
  {
    if (fh.idx() == end_idx)
      break;
    traversed_size += 1;
    bool do_split = false;
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_split = false;
      break;
    case OptStrategy::EOQB:
      do_split = split_face_eoqb(fh);
      break;
    }
    if (do_split)
    {
      splitted_size += 1;
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("split {} / {} edges.", splitted_size, traversed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return splitted_size;
}

size_t LinearSimplifier::relocate_with_priority()
{
  constexpr size_t max_converge_iter = 100; // TODO: move to header configure

#ifdef ONLY_RELOCATE_HIGH_ERROR_VERTEX
  {
    initialize_error_on_verts();

    // calculate verts_error_thres
    std::vector<double> verts_error; verts_error.reserve(mesh->n_vertices());
    for (VertexHandle v : mesh->vertices())
      verts_error.push_back(mesh->data(v).error);
    std::sort(verts_error.begin(), verts_error.end());
    // size_t thres_idx = (size_t)(verts_error.size() * 0.7);
    // verts_error_threshold = verts_error[thres_idx];
    auto first_inf = std::find(verts_error.begin(), verts_error.end(), DBL_MAX);
    verts_error.resize(first_inf - verts_error.begin());
    ASSERT(verts_error.back() < DBL_MAX, "have inf error.");
    double total_error = std::accumulate(verts_error.begin(), verts_error.end(), 0.0);
    total_error *= 0.9;
    double accumulate_error = 0.0;
    for (auto ve_it = verts_error.rbegin(); ve_it != verts_error.rend();ve_it++)
    {
      accumulate_error += *ve_it;
      if (accumulate_error > total_error)
      {
        verts_error_threshold = *ve_it;
        break;
      }
    }
  }
#endif
  init_verts_queue_to_relocate();

  size_t relocated_size = 0;
  size_t converge_iter = 0;
  while (!verts_queue.empty())
  {
    const auto& top = verts_queue.top();
    if (mesh->status(top.handle).deleted()
      || top.timestamp < verts_timestamp[top.handle.idx()]
      || top.timestamp > 15)
    {
      verts_queue.pop();
      continue;
    }
    // try relocate
    auto relocater = new_relocater();
    if (!relocater.init(top.handle, 1/*config.Np*/))
    {
      verts_queue.pop();
      continue;
    }
    // optimize
    Vec3d opt_pos;
    double error_before, error_after, quality_before, quality_after;
    optimizer.error_optimized_quality_bounded(
      static_cast<LLocalOperation*>(&relocater), config.quality_bound, current_error,
      opt_pos, error_before, error_after, quality_before, quality_after);
    // decrease error ?
    if (error_after < error_before)
    {
      // do and update
      relocater.just_do_it(opt_pos);
      double error_change = error_after - error_before;
      current_error += error_change;

      relocated_size += 1;
      // log and output
      Logger::dev_logger->trace("current error is {:.5f}", current_error);
      if (error_to_file)
      {
        operation_cnt += 1;
        fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
      }
      // converge ?
      if (error_change > -current_error * 0.0001)
      {
        if (++converge_iter == max_converge_iter)
          break;
      }
      else converge_iter = 0;
      // update affected vertices
      update_verts_to_relocate(top.handle);
      for (VertexHandle vv : mesh->vv_range(top.handle))
        update_verts_to_relocate(vv);
    }
    verts_timestamp[top.handle.idx()] += 1;
    verts_queue.pop();
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("relocate {} vertices.", relocated_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return relocated_size;
}

void LinearSimplifier::init_verts_queue_to_relocate()
{
  verts_queue = VertexQueue();
  verts_timestamp.clear();
  verts_timestamp.resize(mesh->n_vertices(), 0);

  for (VertexHandle vh : mesh->vertices())
  {
    update_verts_to_relocate(vh);
  }
}

void LinearSimplifier::update_verts_to_relocate(VertexHandle vh)
{
  uint32_t& current_timestamp = verts_timestamp[vh.idx()];
  current_timestamp++;

  if (mesh->status(vh).deleted())
    return;
#ifdef ONLY_RELOCATE_HIGH_ERROR_VERTEX
  if (mesh->data(vh).error < verts_error_threshold)
    return;
#endif

  // simply relocate vertices with larger error.
  verts_queue.emplace(vh, -mesh->data(vh).error, Vec3d(), current_timestamp);
}

size_t LinearSimplifier::collapse_with_priority()
{
  init_edges_queue_to_collapse();
  Logger::user_logger->info("initialize {} edges to collapse.", edges_queue.size());

  size_t collapsed_size = 0;
  while (!edges_queue.empty())
  {
    const auto& top = edges_queue.top();
    if (mesh->status(top.handle).deleted() || top.timestamp < edges_timestamp[top.handle.idx()])
    {
      edges_queue.pop();
      continue;
    }
    if (current_error + top.error_change > config.error_bound)
      break;
    // do collapse
    EdgeHandle eh = top.handle;
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    VertexHandle v0 = mesh->from_vertex_handle(hh), v1 = mesh->to_vertex_handle(hh);
    VertexHandle collapse_vh;
    if (!mesh->is_boundary(eh) && !mesh->is_boundary(v0) && !mesh->is_boundary(v1))
    {
      auto collapser = new_collapser();
      ASSERT(collapser.init(eh, 1/*config.Np*/), "fail to initialize.");  // TODO: accellerate this.
      collapser.just_do_it(top.pos);
      collapse_vh = collapser.collapse_vertex();
    }
    else if (!mesh->is_boundary(eh))
    {
      auto collapser = new_bv_collapser();
      ASSERT(collapser.init(eh, 1/*config.Np*/), "fail to initialize.");
      collapser.just_do_it(top.pos);
      collapse_vh = collapser.collapse_vertex();
    }
    else
    {
      auto collapser = new_be_collapser();
      ASSERT(collapser.init(eh, 1/*config.Np*/), "fail to initialize.");
      collapser.just_do_it(top.pos);
      collapse_vh = collapser.collapse_vertex();
    }
    current_error += top.error_change;
    edges_timestamp[top.handle.idx()] += 1;
    edges_queue.pop();
    // update affected edges
    std::set<EdgeHandle> affected_edges;
    find_edges_affected_by_collapse(collapse_vh, affected_edges);
    for (EdgeHandle eh : affected_edges)
    {
      update_edge_to_collapse(eh);
    }
    collapsed_size += 1;
    Logger::dev_logger->trace("current error is {}", current_error);
    if (collapsed_size % 20 == 0)
      Logger::user_logger->info("collapse {} edges.", collapsed_size);
    if (error_to_file)
    {
      operation_cnt += 1;
      fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
    }
  }
  mesh->garbage_collection();
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("collapse {} edges.", collapsed_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return collapsed_size;
}

void LinearSimplifier::init_edges_queue_to_collapse()
{
  edges_queue = EdgeQueue();
  edges_timestamp.clear();
  edges_timestamp.resize(mesh->n_edges(), 0);

  for (EdgeHandle eh : mesh->edges())
  {
    update_edge_to_collapse(eh);
    if (eh.idx() % 100 == 0)
      Logger::dev_logger->trace("traverse {} edges to collapse", eh.idx());
  }
}

void LinearSimplifier::find_edges_affected_by_collapse(
  VertexHandle collapse_vh, std::set<EdgeHandle>& edges)
{
  edges.clear();
  for (VertexHandle vv : mesh->vv_range(collapse_vh))
    for (EdgeHandle vve : mesh->ve_range(vv))
      edges.insert(vve);
}

void LinearSimplifier::update_edge_to_collapse(EdgeHandle eh)
{
  uint32_t& current_timestamp = edges_timestamp[eh.idx()];
  current_timestamp++;

  if (mesh->status(eh).deleted())
    return;
  HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
  VertexHandle v0 = mesh->from_vertex_handle(hh), v1 = mesh->to_vertex_handle(hh);
  if (!mesh->is_boundary(eh) && mesh->is_boundary(v0) && mesh->is_boundary(v1))
    return;

  Vec3d opt_pos;
  double error_before, error_after, quality_before, quality_after;
  // We use "Error Optimized Quality Bounded".
  if (!mesh->is_boundary(eh) && !mesh->is_boundary(v0) && !mesh->is_boundary(v1))
  {
    auto collapser = new_collapser();
    // check if collapse ok
    if (!collapser.init(eh, 1/*config.Np*/))
      return;
    // check valence
    if (collapser.valence_after_collapse() > config.max_valence)
      return;
    // try to optimize
    optimizer.error_optimized_quality_bounded(
      static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error,
      opt_pos, error_before, error_after, quality_before, quality_after);
    // updpate
    if (error_after < DBL_MAX)
      edges_queue.emplace(eh, error_after - error_before, opt_pos, current_timestamp);
  }
  else if (!mesh->is_boundary(eh))
  {
    auto collapser = new_bv_collapser();
    // check if collapse ok
    if (!collapser.init(eh, 1/*config.Np*/))
      return;
    // check valence
    if (collapser.valence_after_collapse() > config.max_valence)
      return;
    if (collapser.is_collapsed_to_corner())
    {
      // we can't optimize the position
      if (!collapser.check_constraints())
        return;
      // quality bounded ?
      double quality_after = collapser.quality_after();
      if (quality_after >= config.quality_bound)
        return;
      // error change
      double error_before = collapser.error_before();
      double error_after = collapser.error_after();

      if (error_after < DBL_MAX)
        edges_queue.emplace(eh, error_after - error_before, collapser.corner_position(), current_timestamp);
    }
    else
    {
      // try to optimize
      optimizer.error_optimized_quality_bounded(
        static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error,
        opt_pos, error_before, error_after, quality_before, quality_after);
      // updpate
      if (error_after < DBL_MAX)
        edges_queue.emplace(eh, error_after - error_before, opt_pos, current_timestamp);
    }
  }
  else
  {
    auto collapser = new_be_collapser();
    // check if collapse ok
    if (!collapser.init(eh, 1/*config.Np*/))
      return;
    // check valence
    if (collapser.valence_after_collapse() > config.max_valence)
      return;
    if (collapser.is_collapsed_to_corner())
    {
      if (!collapser.check_constraints())
        return;
      // quality bounded ?
      double quality_after = collapser.quality_after();
      if (quality_after > config.quality_bound)
        return;
      // error change
      double error_before = collapser.error_before();
      double error_after = collapser.error_after();

      if (error_after < DBL_MAX)
        edges_queue.emplace(eh, error_after - error_before, collapser.corner_position(), current_timestamp);
    }
    else
    {
      // try to optimize
      optimizer.error_optimized_quality_bounded(
        static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error,
        opt_pos, error_before, error_after, quality_before, quality_after);
      // updpate
      if (error_after < DBL_MAX)
        edges_queue.emplace(eh, error_after - error_before, opt_pos, current_timestamp);
    }
  }
}

size_t LinearSimplifier::split_edges_with_priority()
{
#ifdef ONLY_SPLIT_HIGH_ERROR_EDGE
  {
    initialize_error_on_edges();

    // calculate edges_error_thres
    std::vector<double> edges_error; edges_error.reserve(mesh->n_edges());
    for (EdgeHandle e : mesh->edges())
      edges_error.push_back(mesh->data(e).error);
    std::sort(edges_error.begin(), edges_error.end());
    // size_t thres_idx = (size_t)(edges_error.size() * 0.1);
    // edges_error_threshold = edges_error[thres_idx];
    auto first_inf = std::find(edges_error.begin(), edges_error.end(), DBL_MAX);
    edges_error.resize(first_inf - edges_error.begin());
    ASSERT(edges_error.back() < DBL_MAX, "have inf error.");
    double total_error = std::accumulate(edges_error.begin(), edges_error.end(), 0.0);
    total_error *= 0.9;
    double accumulate_error = 0.0;
    for (auto ee_it = edges_error.rbegin(); ee_it != edges_error.rend();ee_it++)
    {
      accumulate_error += *ee_it;
      if (accumulate_error > total_error)
      {
        edges_error_threshold = *ee_it;
        break;
      }
    }
  }
#endif

  init_edges_queue_to_split();
  Logger::user_logger->info("initialize {} edges to split", edges_queue.size());

  size_t splitted_size = 0;
  std::vector<EdgeHandle> affected_edges;
  while (!edges_queue.empty() && current_error > config.error_bound)
  {
    const auto& top = edges_queue.top();
    if (mesh->status(top.handle).deleted() || top.timestamp < edges_timestamp[top.handle.idx()])
    {
      edges_queue.pop();
      continue;
    }
    if (top.error_change > -current_error * 0.0001)
      break;
    // do split
    auto splitter = new_edge_splitter();
    ASSERT(splitter.init(top.handle, 1/*config.Np*/), "fail to initialize.");  // TODO: accellerate this.
    splitter.just_do_it(top.pos);
    // update
    VertexHandle split_vh = splitter.split_vertex();
    current_error += top.error_change;
    for (uint32_t i = 0;i < mesh->valence(split_vh) - 1;i++)
      edges_timestamp.push_back(0);
    edges_timestamp[top.handle.idx()] += 1;
    edges_queue.pop();
    // update affected edges
    find_edges_affected_by_split(split_vh, affected_edges);
    for (EdgeHandle eh : affected_edges)
    {
      update_edge_to_split(eh);
    }
    // log and output
    splitted_size += 1;
    Logger::dev_logger->trace("current error is {}", current_error);
    if (error_to_file)
    {
      operation_cnt += 1;
      fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("split {} edges.", splitted_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return splitted_size;
}

void LinearSimplifier::init_edges_queue_to_split()
{
  edges_queue = EdgeQueue();
  edges_timestamp.clear();
  edges_timestamp.resize(mesh->n_edges(), 0);

  for (EdgeHandle eh : mesh->edges())
  {
    update_edge_to_split(eh);
  }
}

void LinearSimplifier::find_edges_affected_by_split(VertexHandle split_vh, std::vector<EdgeHandle>& edges)
{
  edges.clear();
  for (HalfedgeHandle voh : mesh->voh_range(split_vh))
  {
    edges.push_back(mesh->edge_handle(voh));
    edges.push_back(mesh->edge_handle(mesh->next_halfedge_handle(voh)));
  }
}

void LinearSimplifier::update_edge_to_split(EdgeHandle eh)
{
  uint32_t& current_timestamp = edges_timestamp[eh.idx()];
  current_timestamp++;

  if (mesh->status(eh).deleted())
    return;
#ifdef ONLY_SPLIT_HIGH_ERROR_EDGE
  if (mesh->data(eh).error < edges_error_threshold)
    return;
#endif

  Vec3d opt_pos;
  double error_before, error_after, quality_before, quality_after;
  // We use "Error Optimized Quality Bounded".
  auto splitter = new_edge_splitter();
  // check if collapse ok
  if (!splitter.init(eh, 1/*config.Np*/))
    return;
  // check valence
  // if (splitter.split_would_cause_over_valence(config.max_valence))
  //  return;
  // try to optimize
  optimizer.error_optimized_quality_bounded(
    static_cast<LLocalOperation*>(&splitter), config.split_quality_bound, current_error,
    opt_pos, error_before, error_after, quality_before, quality_after);
  // updpate
  if (error_after < DBL_MAX)
    edges_queue.emplace(eh, error_after - error_before, opt_pos, current_timestamp);
}

size_t LinearSimplifier::split_faces_with_priority()
{
#ifdef ONLY_SPLIT_HIGH_ERROR_FACE
  {
    // calculate faces_error_thres
    std::vector<double> faces_error; faces_error.reserve(mesh->n_faces());
    for (FaceHandle f : mesh->faces())
      faces_error.push_back(mesh->data(f).error);
    std::sort(faces_error.begin(), faces_error.end());
    auto first_inf = std::find(faces_error.begin(), faces_error.end(), DBL_MAX);
    faces_error.resize(first_inf - faces_error.begin());
    ASSERT(faces_error.back() < DBL_MAX, "have inf error.");
    // size_t thres_idx = (size_t)(faces_error.size() * 0.1);
    // faces_error_threshold = faces_error[thres_idx];
    double total_error = std::accumulate(faces_error.begin(), faces_error.end(), 0.0);
    total_error *= 0.9;
    double accumulate_error = 0.0;
    for (auto fe_it = faces_error.rbegin(); fe_it != faces_error.rend();fe_it++)
    {
      accumulate_error += *fe_it;
      if (accumulate_error > total_error)
      {
        faces_error_threshold = *fe_it;
        break;
      }
    }
  }
#endif
  init_faces_queue_to_split();
  Logger::user_logger->info("initialize {} faces to split.", faces_queue.size());

  size_t splitted_size = 0;
  while (!faces_queue.empty() || current_error > config.error_bound)
  {
    const auto& top = faces_queue.top();
    if (mesh->status(top.handle).deleted() || top.timestamp < faces_timestamp[top.handle.idx()])
    {
      faces_queue.pop();
      continue;
    }
    if (top.error_change > -current_error * 0.0001)
      break;
    // do split
    auto splitter = new_face_splitter();
    ASSERT(splitter.init(top.handle, 1/*config.Np*/), "fail to initialize.");  // TODO: accellerate this.
    splitter.just_do_it(top.pos);
    VertexHandle split_vh = splitter.split_vertex();
    // update
    current_error += top.error_change;
    faces_timestamp.push_back(0);
    faces_timestamp.push_back(0);
    faces_timestamp[top.handle.idx()] += 1;
    faces_queue.pop();
    // update affected faces
    for (HalfedgeHandle voh : mesh->voh_range(split_vh))
    {
      update_face_to_split(mesh->face_handle(voh));
    }
    // log and output
    splitted_size += 1;
    Logger::dev_logger->trace("current error is {}", current_error);
    if (error_to_file)
    {
      operation_cnt += 1;
      fprintf(error_fout, "%zu, %f\n", operation_cnt, std::sqrt(current_error / (image->width * image->height * 3.0)));
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("split {} faces.", splitted_size);
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return splitted_size;
}

void LinearSimplifier::init_faces_queue_to_split()
{
  faces_queue = FaceQueue();
  faces_timestamp.clear();
  faces_timestamp.resize(mesh->n_faces(), 0);

  for (FaceHandle fh : mesh->faces())
  {
    update_face_to_split(fh);
  }
}

void LinearSimplifier::update_face_to_split(FaceHandle fh)
{
  uint32_t& current_timestamp = faces_timestamp[fh.idx()];
  current_timestamp++;

  if (mesh->status(fh).deleted())
    return;

  if (mesh->data(fh).error < faces_error_threshold)
    return;

  Vec3d opt_pos;
  double error_before, error_after, quality_before, quality_after;
  // We use "Error Optimized Quality Bounded".
  auto splitter = new_face_splitter();
  // check if collapse ok
  if (!splitter.init(fh, 1/*config.Np*/))
    return;
  // check valence
  // if (splitter.split_would_cause_over_valence(config.max_valence))
  //  return;
  // try to optimize
  optimizer.error_optimized_quality_bounded(
    static_cast<LLocalOperation*>(&splitter), config.split_quality_bound, current_error,
    opt_pos, error_before, error_after, quality_before, quality_after);
  // updpate
  if (error_after < DBL_MAX)
    faces_queue.emplace(fh, error_after - error_before, opt_pos, current_timestamp);
}

bool LinearSimplifier::relocate_vertex_ebqo(VertexHandle vh)
{
  auto relocater = new_relocater();
  if (!relocater.init(vh, 1 /*config.Np*/))
    return false;

  return optimizer.do_ebqo(
    static_cast<LLocalOperation*>(&relocater), config.error_bound, current_error);
}

bool LinearSimplifier::collapse_edge_ebqo(EdgeHandle eh)
{
  auto collapser = new_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;;
  // try to optimize
  return optimizer.do_ebqo(
    static_cast<LLocalOperation*>(&collapser), config.error_bound, current_error);
}

bool LinearSimplifier::collapse_boundary_edge_ebqo(EdgeHandle eh)
{
  auto collapser = new_be_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;
  if (collapser.is_collapsed_to_corner())
  {
    // we can't optimize the position.
    if (!collapser.check_constraints())
      return false;

    // error bounded ?
    double error_before = collapser.error_before();
    double error_after = collapser.error_after();
    if (current_error - error_before + error_after > config.error_bound)
      return false;

    // quality optimized ?
    double quality_before = collapser.quality_before();
    double quality_after = collapser.quality_after();
    if (quality_after >= quality_before)
      return false;

    collapser.just_do_it(Vec3d());  // the corner point shouldn't be relocated.
    current_error = current_error - error_before + error_after;
    return true;
  }
  else
  {
    // try to optimize
    return optimizer.do_ebqo(
      static_cast<LLocalOperation*>(&collapser), config.error_bound, current_error);
  }
}

bool LinearSimplifier::collapse_boundary_vertex_ebqo(EdgeHandle eh)
{
  auto collapser = new_bv_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;
  if (collapser.is_collapsed_to_corner())
  {
    // we can't optimize the position
    if (!collapser.check_constraints())
      return false;

    // error bounded ?
    double error_before = collapser.error_before();
    double error_after = collapser.error_after();
    if (current_error - error_before + error_after > config.error_bound)
      return false;

    // quality optimized ?
    double quality_before = collapser.quality_before();
    double quality_after = collapser.quality_after();
    if (quality_after >= quality_before)
      return false;

    collapser.just_do_it(Vec3d());  // the corner point shouldn't be relocated.
    current_error = current_error - error_before + error_after;
    return true;
  }
  else
  {
    // try to optimize
    return optimizer.do_ebqo(
      static_cast<LLocalOperation*>(&collapser), config.error_bound, current_error);
  }
}

bool LinearSimplifier::flip_edge_ebqo(EdgeHandle eh)
{
  auto flipper = new_flipper();
  // check if flip ok
  if (!flipper.init(eh))
    return false;
  // check valence
  if (!flipper.flip_would_decrease_valence())
    return false;
  if (!flipper.check_constraints())
    return false;
  // error bounded ?
  double error_before = flipper.error_before();
  double error_after = flipper.error_after();
  if (current_error - error_before + error_after > config.error_bound)
    return false;
  // quality optimized ?
  double quality_before = flipper.quality_before();
  double quality_after = flipper.quality_after();
  if (quality_after >= quality_before)
    return false;

  flipper.just_do_it();
  current_error = current_error - error_before + error_after;
  return true;
}

bool LinearSimplifier::relocate_vertex_eoqb(VertexHandle vh)
{
  auto relocater = new_relocater();
  if (!relocater.init(vh, 1 /*config.Np*/))
    return false;

  return optimizer.do_eoqb(
    static_cast<LLocalOperation*>(&relocater), config.quality_bound, current_error);
}

bool LinearSimplifier::collapse_edge_eoqb(EdgeHandle eh)
{
  auto collapser = new_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;;
  // try to optimize
  return optimizer.do_eoqb(
    static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error);
}

bool LinearSimplifier::collapse_boundary_edge_eoqb(EdgeHandle eh)
{
  auto collapser = new_be_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;
  if (collapser.is_collapsed_to_corner())
  {
    if (!collapser.check_constraints())
      return false;
    // error optimized ?
    double error_before = collapser.error_before();
    double error_after = collapser.error_after();
    if (error_after >= error_before)
      return false;
    // quality bounded ?
    double quality_after = collapser.quality_after();
    if (quality_after > config.quality_bound)
      return false;

    collapser.just_do_it(Vec3d());  // the corner point shouldn't be relocated.
    current_error = current_error - error_before + error_after;
    Logger::dev_logger->trace("current error is {:.5f}", current_error);
    return true;
  }
  else
  {
    // try to optimize
    return optimizer.do_eoqb(
      static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error);
  }
}

bool LinearSimplifier::collapse_boundary_vertex_eoqb(EdgeHandle eh)
{
  auto collapser = new_bv_collapser();
  // check if collapse ok
  if (!collapser.init(eh, 1/*config.Np*/))
    return false;
  // check valence
  if (collapser.valence_after_collapse() > config.max_valence)
    return false;
  if (collapser.is_collapsed_to_corner())
  {
    if (!collapser.check_constraints())
      return false;
    // error optimized ?
    double error_before = collapser.error_before();
    double error_after = collapser.error_after();
    if (error_after >= error_before)
      return false;
    // quality bounded ?
    double quality_after = collapser.quality_after();
    if (quality_after > config.quality_bound)
      return false;

    collapser.just_do_it(Vec3d());  // the corner point shouldn't be relocated.
    current_error = current_error - error_before + error_after;
    Logger::dev_logger->trace("current error is {:.5f}", current_error);
    return true;
  }
  else
  {
    // try to optimize
    return optimizer.do_eoqb(
      static_cast<LLocalOperation*>(&collapser), config.quality_bound, current_error);
  }
}

bool LinearSimplifier::flip_edge_eoqb(EdgeHandle eh)
{
  auto flipper = new_flipper();
  // check if flip ok
  if (!flipper.init(eh))
    return false;
  // check valence
  if (flipper.flip_would_exceed_max_valence(config.max_valence))
    return false;
  if (!flipper.check_constraints())
    return false;
  // quality bounded ?
  double quality_after = flipper.quality_after();
  if (quality_after > config.quality_bound)
    return false;
  // error optimized ?
  double error_before = flipper.error_before();
  double error_after = flipper.error_after();
  if (error_before < error_after)
    return false;

  flipper.just_do_it();
  current_error = current_error - error_before + error_after;
  return true;
}

bool LinearSimplifier::split_edge_eoqb(EdgeHandle eh)
{
  auto splitter = new_edge_splitter();
  // check if flip ok
  if (!splitter.init(eh, 1))
    return false;
  // check valence
  // if (splitter.split_would_cause_over_valence(config.max_valence))
  //  return false;
  return optimizer.do_eoqb(
    static_cast<LLocalOperation*>(&splitter), config.split_quality_bound, current_error);
}

bool LinearSimplifier::split_face_eoqb(FaceHandle fh)
{
  auto splitter = new_face_splitter();
  // check if flip ok
  if (!splitter.init(fh, 1))
    return false;
  // check valence
  // if (splitter.split_would_cause_over_valence(config.max_valence))
  //  return false;
  return optimizer.do_eoqb(
    static_cast<LLocalOperation*>(&splitter), config.split_quality_bound, current_error);
  return true;
}

}// namespace ImageTriSimp
}// namespace GCLF