#include "CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

void CurvedSimplifier::initialize(ImageT& input_image, LMeshT& input_mesh, ParamTriangulator& param)
{
  initialize_uvw_samples(2);
  // set configurations
  set_color_type(param.color_type);
  set_error_type(param.error_type);
  set_min_angle(param.param_curved.min_angle);
  set_max_curvature(param.param_curved.max_curvature);
  set_max_valence(param.max_valence);
  config.max_simplify_iter = param.param_curved.max_simplify_iter;
  config.convergence_collapse_number = param.param_curved.convergence_collapse_number;
  // initialize DE optimizer.
  optimizer.initialize(this, param);
  // initialize mesh and image.
  image = std::make_unique<ImageT>(input_image);
  mesh = std::make_unique<BMeshT>();
  mesh->initialize(&input_mesh);
  // mesh->scale_to(*image);
  // initialize color and error
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  // Set the bounding error
  current_error = lp_error();

  Logger::user_logger->info("Initial mesh: v {}, e {}, f {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Initial root_mean_lp_error {}", root_mean_lp_error());
}

void CurvedSimplifier::initialize(ImageT& input_image, BMeshT& input_mesh, ParamTriangulator& param)
{
  initialize_uvw_samples(2);
  // set configurations
  set_color_type(param.color_type);
  set_error_type(param.error_type);
  set_min_angle(param.param_curved.min_angle);
  set_max_curvature(param.param_curved.max_curvature);
  set_max_valence(param.max_valence);
  config.max_simplify_iter = param.param_curved.max_simplify_iter;
  config.convergence_collapse_number = param.param_curved.convergence_collapse_number;
  // initialize DE optimizer.
  optimizer.initialize(this, param);
  // initialize mesh and image.
  image = std::make_unique<ImageT>(input_image);
  mesh = std::make_unique<BMeshT>(input_mesh);
  // mesh->scale_to(*image);
  // initialize color and error
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  // Set the bounding error
  current_error = lp_error();

  Logger::user_logger->info("Initial mesh: v {}, e {}, f {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Initial root_mean_lp_error {}", root_mean_lp_error());
}

void CurvedSimplifier::set_color_type(const std::string& color_type)
{
  config.color_type = color_type;
  if (config.color_type == "constant")
  {
    behavior.calc_color_on_face.type =
      CurvedSimplifierBehaviors::CCalcColorOnFace::Type::Constant;
  }
  else if (config.color_type == "linear")
  {
    behavior.calc_color_on_face.type =
      CurvedSimplifierBehaviors::CCalcColorOnFace::Type::Linear;
  }
  else if (config.color_type == "quadratic")
  {
    behavior.calc_color_on_face.type =
      CurvedSimplifierBehaviors::CCalcColorOnFace::Type::Quadratic;
  }
  else
  {
    ASSERT(false, "wrong color type");
  }
}

void CurvedSimplifier::set_error_type(const std::string& error_type)
{
  config.error_type = error_type;
  if (config.error_type == "l1")
  {
    behavior.calc_error_on_face.type =
      CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L1Norm;
    behavior.calc_error_on_mesh.p = 1.0;
  }
  else if (config.error_type == "l2")
  {
    behavior.calc_error_on_face.type =
      CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L2Norm;
    behavior.calc_error_on_mesh.p = 2.0;
  }
  else if (config.error_type == "l4")
  {
    behavior.calc_error_on_face.type =
      CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L4Norm;
    behavior.calc_error_on_mesh.p = 4.0;
  }
  else
  {
    ASSERT(false, "wrong error type");
  }
}

void CurvedSimplifier::set_error_bound(double error_bound)
{
  switch (behavior.calc_error_on_face.type)
  {
  case CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L1Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = error_bound / 255. * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = error_bound * image->width * image->height * 3;
  #endif
    break;
  case CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L2Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = pow(error_bound / 255., 2.) * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = pow(error_bound, 2.) * image->width * image->height * 3;
  #endif
    break;
  case CurvedSimplifierBehaviors::CCalcErrorOnFace::Type::L4Norm:
  #if defined(CHANNEL_0_1)
    config.error_bound = pow(error_bound / 255., 4.) * image->width * image->height * 3;
  #elif defined(CHANNEL_0_255)
    config.error_bound = pow(error_bound, 4.) * image->width * image->height * 3;
  #endif
    break;
  }
}

void CurvedSimplifier::set_max_valence(uint32_t mv)
{
  config.max_valence = mv;
}

void CurvedSimplifier::set_min_angle(double min_angle)
{
  config.min_angle = std::cos(min_angle / 180.0 * M_PI);
}

void CurvedSimplifier::set_max_curvature(double max_curvature)
{
  config.max_curvature = std::cos(max_curvature / 180.0 * M_PI);
}

void CurvedSimplifier::simplify()
{
  size_t iter = 0;
  while (iter < config.max_simplify_iter)
  {
    flip(OptStrategy::EOQB);
    for (size_t i = 0;i < 3;i++)
    {
      relocate_edges(OptStrategy::EOQB);
      relocate_vertices(OptStrategy::EOQB);
    }
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
std::vector<Vec3d> CurvedSimplifier::intersect_points(const Vec2d& ray_start)
{
  CAssignPixelsToFace assigner;
  Eigen::MatrixXi mask; image->reset_mask(mask);
  assigner.initialize(mesh.get(), image.get(), mask, std::vector<Vec2i>());

  // Use tree to search possible intersected curves.
  auto possible_curves = assigner.tree->intersect_curves(ray_start);

  // Find the intersection points between rays and edges.
  auto inter_points = assigner.find_intersection_points(ray_start, possible_curves);

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

void CurvedSimplifier::visualize_error()
{
  visualize_image();

  // calculate error for each pixel
  Eigen::MatrixXd error_mat;
  error_mat.resize(image->height, image->width);
  double min_error = DBL_MAX, max_error = -DBL_MAX;
  for (int y = 0;y < image->height;y++)
  {
    for (int x = 0;x < image->width;x++)
    {
      const ImageT::Color& oc = image->pixel_color(x, y);     // original color
      const ImageT::Color& ac = out_image->pixel_color(x, y); // approximate color
      double x_diff = fabs(oc.x() - ac.x());
      double y_diff = fabs(oc.y() - ac.y());
      double z_diff = fabs(oc.z() - ac.z());
      error_mat(y, x) = pow(x_diff, behavior.calc_error_on_mesh.p) + pow(y_diff, behavior.calc_error_on_mesh.p) + pow(z_diff, behavior.calc_error_on_mesh.p);
      if (error_mat(y, x) > max_error)
        max_error = error_mat(y, x);
      if (error_mat(y, x) < min_error)
        min_error = error_mat(y, x);
    }
  }
  // fstream fout;
  // fout.open("C:/Everything/MyData/ActiveData/Research/Projects/ImageTriangulationSimplification/swan/errors.txt", std::fstream::out);
  // if (fout.is_open())
  // {
  //   for (int y = 0;y < image->height;y++)
  //     for (int x = 0;x < image->width;x++)
  //       fout << error_mat(y, x) << std::endl;
  //   fout.close();
  // }
  // normalize erros
  double sum = 0.0, avg = 0.0, std_dev = 0.0;
  // average
  for (int y = 0;y < image->height;y++)
    for (int x = 0;x < image->width;x++)
      sum += error_mat(y, x);
  avg = sum / (image->height * image->width);
  // standard deviation
  sum = 0.0;
  for (int y = 0;y < image->height;y++)
    for (int x = 0;x < image->width;x++)
      sum += pow(error_mat(y, x) - avg, 2.);
  std_dev = sqrt(sum / (image->height * image->width));
  Logger::user_logger->info("error: average {}, standard deviation {}", avg, std_dev);

  // convert to color
  ColorMap color_map(avg + 3 * std_dev, avg - 3 * std_dev);
  for (int y = 0;y < image->height;y++)
  {
    for (int x = 0;x < image->width;x++)
    {
      Vec3d error_color = color_map.MapToColor(error_mat(y, x));
      ImageT::Color& out_color = out_image->pixel_color(x, y);
      out_color.x() = error_color.x() * 255.;
      out_color.y() = error_color.y() * 255.;
      out_color.z() = error_color.z() * 255.;
    }
  }
}

void CurvedSimplifier::visualize_image()
{
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
}

void CurvedSimplifier::update_color_error()
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

size_t CurvedSimplifier::relocate_vertices(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t relocated_size = 0;
  for (VertexHandle vh : mesh->vertices())
  {
    if (mesh->is_boundary(vh))
      continue; // TODO: relocate boundary vertex!
    traversed_size += 1;
    bool do_relocate = false;

    auto relocater = new_v_relocater();
    relocater.init(vh, 1/*config.Np*/);
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_relocate = optimizer.do_ebqo(
        static_cast<CLocalOperation*>(&relocater), config.error_bound, current_error);
      break;
    case OptStrategy::EOQB:
      do_relocate = optimizer.do_eoqb(
        static_cast<CLocalOperation*>(&relocater), config.min_angle, config.max_curvature, current_error);
      break;
    }
    if (do_relocate)
    {
      relocated_size += 1;
      // if (relocated_size % 10 == 0)
      //  Logger::user_logger->info("relocate {} / {} vertices.", relocated_size, traversed_size);
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return relocated_size;
}

size_t CurvedSimplifier::relocate_edges(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t relocated_size = 0;
  for (EdgeHandle eh : mesh->edges())
  {
    if (mesh->is_boundary(eh))
      continue;
    traversed_size += 1;
    bool do_relocate = false;

    auto relocater = new_e_relocater();
    relocater.init(eh, 1/*config.Np*/);
    switch (opt)
    {
    case OptStrategy::EBQO:
      do_relocate = optimizer.do_ebqo(
        static_cast<CLocalOperation*>(&relocater), config.error_bound, current_error);
      break;
    case OptStrategy::EOQB:
      do_relocate = optimizer.do_eoqb(
        static_cast<CLocalOperation*>(&relocater), config.min_angle, config.max_curvature, current_error);
      break;
    }
    if (do_relocate)
    {
      relocated_size += 1;
      // if (relocated_size % 10 == 0)
      //   Logger::user_logger->info("relocate {} / {} edges.", relocated_size, traversed_size);
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return relocated_size;
}

size_t CurvedSimplifier::collapse(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t collapsed_size = 0;

  for (HalfedgeHandle hh : mesh->halfedges())
  {
    EdgeHandle eh = mesh->edge_handle(hh);
    VertexHandle from_vh = mesh->from_vertex_handle(hh), to_vh = mesh->to_vertex_handle(hh);
    if (mesh->status(eh).deleted())
      continue;
    if (!mesh->is_boundary(eh) && mesh->is_boundary(from_vh) && mesh->is_boundary(to_vh))
      continue;

    traversed_size += 1;
    bool do_collapse = false;
    if (!mesh->is_boundary(eh) && !mesh->is_boundary(from_vh) && !mesh->is_boundary(to_vh))
    {
      auto collapser = new_collapser();
      // try collapse to to_vh
      HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
      // check if collapse ok
      if (!collapser.init(hh, 1/*config.Np*/))
        continue;
      // check valence
      if (collapser.valence_after_collapse() > config.max_valence)
        continue;
      // try to optimize
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = optimizer.do_ebqo(
          static_cast<CLocalOperation*>(&collapser), config.error_bound, current_error);
        break;
      case OptStrategy::EOQB:
        do_collapse = optimizer.do_eoqb(
          static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error);
        break;
      }
    }
    else if (!mesh->is_boundary(eh))
    {
      auto collapser = new_bv_collapser();
      // check if collapse ok
      if (!collapser.init(eh, 1/*config.Np*/))
        continue;
      // check valence
      if (collapser.valence_after_collapse() > config.max_valence)
        continue;
      // try to optimize
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = optimizer.do_ebqo(
          static_cast<CLocalOperation*>(&collapser), config.error_bound, current_error);
        break;
      case OptStrategy::EOQB:
        do_collapse = optimizer.do_eoqb(
          static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error);
        break;
      }
    }
    else
    {
      auto collapser = new_be_collapser();
      // check if collapse ok
      if (!collapser.init(eh, 1/*config.Np*/))
        continue;
      // check valence
      if (collapser.valence_after_collapse() > config.max_valence)
        continue;
      // try to optimize
      switch (opt)
      {
      case OptStrategy::EBQO:
        do_collapse = optimizer.do_ebqo(
          static_cast<CLocalOperation*>(&collapser), config.error_bound, current_error);
        break;
      case OptStrategy::EOQB:
        do_collapse = optimizer.do_eoqb(
          static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error);
        break;
      }
    }
    if (do_collapse)
    {
      // succeed to optimize
      collapsed_size += 1;
      if (collapsed_size % 10 == 0)
        Logger::user_logger->info("collapse {} / {} edges.", collapsed_size, traversed_size);
    }
  }
  mesh->garbage_collection();
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return collapsed_size;
}

size_t CurvedSimplifier::flip(OptStrategy opt)
{
  size_t traversed_size = 0;
  size_t flipped_size = 0;
  for (EdgeHandle eh : mesh->edges())
  {
    traversed_size += 1;
    bool do_flip = false;
    auto flipper = new_flipper();
    // check if flip ok
    if (!flipper.init(eh, 1/*config.Np*/))
      continue;
    // try to optimize
    switch (opt)
    {
    case OptStrategy::EBQO:
      if (flipper.flip_would_decrease_valence())
      {
        do_flip = optimizer.do_ebqo(
          static_cast<CLocalOperation*>(&flipper), config.error_bound, current_error);
      }
      break;
    case OptStrategy::EOQB:
      if (!flipper.flip_would_exceed_max_valence(config.max_valence))
      {
        do_flip = optimizer.do_eoqb(
          static_cast<CLocalOperation*>(&flipper), config.min_angle, config.max_curvature, current_error);
      }
      break;
    }
    if (do_flip)
    {
      // succeed to optimize
      flipped_size += 1;
      // if (flipped_size % 10 == 0)
      //   Logger::user_logger->info("flip {} / {} edges.", flipped_size, traversed_size);
    }
  }
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return flipped_size;
}

size_t CurvedSimplifier::collapse_with_priority()
{
  init_edges_queue_to_collapse();

  size_t collapsed_size = 0;
  while (!edges_queue.empty())
  {
    const auto& top = edges_queue.top();
    if (mesh->status(top.eh).deleted() || top.timestamp < edges_timestamp[top.eh.idx()])
    {
      edges_queue.pop();
      continue;
    }
    if (current_error + top.error_change > config.error_bound)
      break;
    // do collapse
    EdgeHandle eh = top.eh;
    HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
    VertexHandle v0 = mesh->from_vertex_handle(hh), v1 = mesh->to_vertex_handle(hh);
    VertexHandle collapse_vh;
    if (!mesh->is_boundary(eh) && !mesh->is_boundary(v0) && !mesh->is_boundary(v1))
    {
      auto collapser = new_collapser();
      ASSERT(collapser.init(top.hh, 1/*config.Np*/), "fail to initialize.");  // TODO: accellerate this.
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
    Logger::dev_logger->trace("error change is {:.5f}, current error is {:.5f}", top.error_change, current_error);
    edges_timestamp[top.eh.idx()] += 1;
    edges_queue.pop();
    // update affected edges
    std::set<EdgeHandle> affected_edges;
    find_edges_affected_by_collapse(collapse_vh, affected_edges);
    for (EdgeHandle eh : affected_edges)
    {
      update_edge_in_queue(eh);
    }
    collapsed_size += 1;
    if (collapsed_size % 10 == 0)
      Logger::user_logger->info("collapse {} edges.", collapsed_size);
  }
  mesh->garbage_collection();
  update_color_error();
  out_image = std::make_unique<ImageT>(behavior.mesh_to_image());
  Logger::user_logger->info("Final mesh: vertices {} edges {} faces {}", mesh->n_vertices(), mesh->n_edges(), mesh->n_faces());
  Logger::user_logger->info("Final root_mean_lp_error {}", root_mean_lp_error());
  return collapsed_size;
}

void CurvedSimplifier::init_edges_queue_to_collapse()
{
  edges_queue = std::priority_queue<EdgePriority>();
  edges_timestamp.resize(mesh->n_edges(), 0);

  for (EdgeHandle eh : mesh->edges())
  {
    update_edge_in_queue(eh);
  }
  Logger::user_logger->info("initialize {} edges to collapse.", edges_queue.size());
}

void CurvedSimplifier::find_edges_affected_by_collapse(
  VertexHandle collapse_vh, std::set<EdgeHandle>& edges)
{
  edges.clear();
  for (VertexHandle vv : mesh->vv_range(collapse_vh))
    for (EdgeHandle vve : mesh->ve_range(vv))
      edges.insert(vve);
}

void CurvedSimplifier::update_edge_in_queue(EdgeHandle eh)
{
  uint32_t& current_timestamp = edges_timestamp[eh.idx()];
  current_timestamp++;

  if (mesh->status(eh).deleted())
    return;
  HalfedgeHandle hh = mesh->halfedge_handle(eh, 0);
  VertexHandle v0 = mesh->from_vertex_handle(hh), v1 = mesh->to_vertex_handle(hh);
  if (!mesh->is_boundary(eh) && mesh->is_boundary(v0) && mesh->is_boundary(v1))
    return;

  Points opt_pos;
  double error_before, error_after;
  // We use "Error Optimized Quality Bounded".
  if (!mesh->is_boundary(eh) && !mesh->is_boundary(v0) && !mesh->is_boundary(v1))
  {
    HalfedgeHandle hh0 = mesh->halfedge_handle(eh, 0);
    HalfedgeHandle hh1 = mesh->halfedge_handle(eh, 1);
    {// try collapse hh0
      auto collapser = new_collapser();
      // check if collapse ok (only need check once)
      if (!collapser.init(hh0, 1/*config.Np*/))
        return;
      // check valence (only need check once)
      if (collapser.valence_after_collapse() > config.max_valence)
        return;
      // try to optimize
      optimizer.error_optimized_quality_bounded(
        static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error,
        opt_pos, error_before, error_after);
    }
    Points opt_pos_;
    double error_before_, error_after_;
    {// try collapse hh1
      auto collapser = new_collapser();
      // check if collapse ok (only need check once)
      ASSERT(collapser.init(hh1, 1/*config.Np*/), "fail to initialize hh1.");
      // try to optimize
      optimizer.error_optimized_quality_bounded(
        static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error,
        opt_pos_, error_before_, error_after_);
    }
    if (error_after < DBL_MAX && error_after <= error_after_)
      edges_queue.emplace(eh, hh0, error_after - error_before, std::move(opt_pos), current_timestamp);
    else if (error_after_ < DBL_MAX && error_after_ < error_after)
      edges_queue.emplace(eh, hh1, error_after_ - error_before_, std::move(opt_pos_), current_timestamp);
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
    // try to optimize
    optimizer.error_optimized_quality_bounded(
      static_cast<CLocalOperation*>(&collapser), config.min_angle, config.max_curvature, current_error,
      opt_pos, error_before, error_after);
    // updpate
    if (error_after < DBL_MAX)
      edges_queue.emplace(eh, mesh->InvalidHalfedgeHandle, error_after - error_before, std::move(opt_pos), current_timestamp);
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
    // try to optimize
    optimizer.error_optimized_quality_bounded(
      static_cast<CLocalOperation*>(&collapser), config.max_curvature, config.min_angle, current_error,
      opt_pos, error_before, error_after);
    // updpate
    if (error_after < DBL_MAX)
      edges_queue.emplace(eh, mesh->InvalidHalfedgeHandle, error_after - error_before, std::move(opt_pos), current_timestamp);
  }
}

}// namespace ImageTriSimp
}// namespace GCLF