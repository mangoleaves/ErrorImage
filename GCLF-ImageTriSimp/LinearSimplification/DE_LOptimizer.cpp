#include "DE_LOptimizer.h"
#include "LinearSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

void DE_LOptimizer::initialize(LinearSimplifier* _simplifier, ParamTriangulator& param)
{
  simplifier = _simplifier;
  Np = param.Np;
  mutation_scale = param.mutation_scale;
  crossover_rate = param.crossover_rate;
  convergence_rate = param.convergence_rate;
  max_consecutive_iter = param.max_consecutive_iter;
  max_DE_iter = param.max_DE_iter;
  init_random();
}

void DE_LOptimizer::init_random()
{
// #define PSEUDO_RANDOM
#ifdef PSEUDO_RANDOM
  init_generator = std::default_random_engine();
  init_dis = std::uniform_real_distribution<double>(0.0, std::nextafter(1.0, DBL_MAX));
  mutate_generator = std::default_random_engine();
  mutate_dis = std::uniform_int_distribution<int>(0, (int)Np - 1);
  cross_generator = std::default_random_engine();
  cross_dis = std::uniform_real_distribution<double>(0.0, std::nextafter(1.0, DBL_MAX));
#else
  std::random_device seed;
  init_generator = std::default_random_engine(seed());
  init_dis = std::uniform_real_distribution<double>(0.0, std::nextafter(1.0, DBL_MAX));
  mutate_generator = std::default_random_engine(seed());
  mutate_dis = std::uniform_int_distribution<int>(0, (int)Np - 1);
  cross_generator = std::default_random_engine(seed());
  cross_dis = std::uniform_real_distribution<double>(0.0, std::nextafter(1.0, DBL_MAX));
#endif
}

void DE_LOptimizer::mutate(
  const BoundingBox& bbox,
  const std::vector<Vec3d>& old_population,
  int best_idx, int old_idx,
  Vec3d& mutation_individual)
{
  int rand1;
  do rand1 = mutate_dis(mutate_generator);
  while (rand1 == old_idx);
  const Vec3d& op1 = old_population[rand1];
  int rand2;
  do rand2 = mutate_dis(mutate_generator);
  while (rand2 == old_idx || rand2 == rand1);
  const Vec3d& op2 = old_population[rand2];
  // int rand3;
  // do rand3 = mutate_dis(mutate_generator);
  // while (rand3 == old_idx || rand3 == rand2 || rand3 == rand1);
  // const Vec3d& op3 = old_population[rand3];

  // mutation_individual = op1 + mutation_scale * (op2 - op3);
  mutation_individual = old_population[old_idx] +
    mutation_scale * (old_population[best_idx] - old_population[old_idx]) +
    mutation_scale * (op1 - op2);
  mutation_individual.maximize(bbox.min());
  mutation_individual.minimize(bbox.max());
}

void DE_LOptimizer::crossover(
  const Vec3d& old_individual,
  const Vec3d& mutation_individual,
  Vec3d& crossover_individual)
{
  crossover_individual.x() = cross_dis(cross_generator) <= crossover_rate ?
    mutation_individual.x() : old_individual.x();
  crossover_individual.y() = cross_dis(cross_generator) <= crossover_rate ?
    mutation_individual.y() : old_individual.y();
  // crossover_individual.z() = cross_dis(cross_generator) <= crossover_rate ?
  //   mutation_individual.z() : old_individual.z();
}

bool DE_LOptimizer::convergence(
  double old_min_error, double new_min_error,
  size_t& consecutive_iter)
{
  if (old_min_error != DBL_MAX)
  {
    double relative_change = (old_min_error - new_min_error) / old_min_error;
    if (relative_change <= convergence_rate)
      consecutive_iter++;
    else
      consecutive_iter = 0;
    if (consecutive_iter == max_consecutive_iter)
    {
      //Logger::dev_logger->trace("loop is convergent.");
      return true;
    }
  }
  return false;
}

void DE_LOptimizer::error_bounded_quality_optimized(
  LLocalOperation* local_op,
  const double error_bound,
  const double current_error,
  Vec3d& opt_pos,
  double& error_before,
  double& error_after,
  double& quality_before,
  double& quality_after)
{
  /****** Initialize DE ******/
  // declare populations
  std::vector<Vec3d> old_population, mutation_population, crossover_population, new_population;
  // declare errors
  error_after = DBL_MAX;
  std::vector<double> errors_after(Np, DBL_MAX);
  // declare qualityies
  quality_after = DBL_MAX;
  std::vector<double> qualities_after(Np, DBL_MAX);
  size_t new_min_quality_idx = 0;
  double new_min_quality = DBL_MAX, old_min_quality = DBL_MAX;
  // define functions
  auto is_error_bounded = [&](double _error_after) -> bool
    {
      return (current_error - error_before + _error_after) <= error_bound;
    };

  // generate initial population
  BoundingBox bbox = local_op->bounding_box();
  Vec3d variable = local_op->get_variable();

  old_population = local_op->get_init_population(init_generator, init_dis, Np);
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
  error_before = local_op->error_before();
  quality_before = local_op->quality_before();

  //#pragma omp parallel for
  for (int i = 0;i < Np;i++)
  {
    if (local_op->update_variable(old_population[i], 0))
    {
      errors_after[i] = local_op->error_after(0);
      qualities_after[i] = is_error_bounded(errors_after[i]) ?
        local_op->quality_after(0) : DBL_MAX;
    }
    else
    {
      errors_after[i] = DBL_MAX;
      qualities_after[i] = DBL_MAX;
    }
  }

  new_min_quality_idx = std::distance(qualities_after.begin(),
    std::min_element(qualities_after.begin(), qualities_after.end()));
  error_after = errors_after[new_min_quality_idx];
  new_min_quality = qualities_after[new_min_quality_idx];
  old_min_quality = new_min_quality;

  // enter loop
  mutation_population = old_population;
  crossover_population = old_population;
  size_t iter = 0, consecutive_iter = 0;
  while (true)
  {
    iter++;
    // generate mutation population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
      mutate(bbox, old_population, new_min_quality_idx, i, mutation_population[i]);
    // generate crossover population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
      crossover(old_population[i], mutation_population[i], crossover_population[i]);
    // generate new population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
    {
      double new_error = local_op->update_variable(crossover_population[i], 0) ? local_op->error_after(0) : DBL_MAX;
      if (is_error_bounded(new_error))
      {
        double new_quality = local_op->quality_after(0);
        if (new_quality < qualities_after[i])
        {
          new_population[i] = crossover_population[i];  // swap faster?
          errors_after[i] = new_error;
          qualities_after[i] = new_quality;
        }
      }
    }
    new_min_quality_idx = std::distance(qualities_after.begin(),
      std::min_element(qualities_after.begin(), qualities_after.end()));
    error_after = errors_after[new_min_quality_idx];
    new_min_quality = qualities_after[new_min_quality_idx];
    // the first end condition : check if loop is convergent
    if (convergence(old_min_quality, new_min_quality, consecutive_iter))
      break;
    // the second end condition : check if loop reaches the maximum iteration number.
    if (iter == max_DE_iter)
    {
      //Logger::dev_logger->trace("maximum iteration number.");
      break;
    }
    // continue to optimize
    old_min_quality = new_min_quality;
    old_population = new_population;
  }

  opt_pos = new_population[new_min_quality_idx];
  quality_after = new_min_quality;
}

bool DE_LOptimizer::do_ebqo(
  LLocalOperation* local_op,
  double error_bound,
  double& current_error)
{
  Vec3d opt_pos;
  double error_before, error_after, quality_before, quality_after;
  error_bounded_quality_optimized(
    local_op, error_bound, current_error,
    opt_pos, error_before, error_after, quality_before, quality_after);

  if ((current_error - error_before + error_after) <= error_bound && quality_after < quality_before)
  {
    local_op->just_do_it(opt_pos);
    current_error = current_error - error_before + error_after;
    Logger::dev_logger->trace("quality before is {:.5f}, after is {:.5f}", quality_before, quality_after);
    Logger::dev_logger->trace("current error is {:.5f}", current_error);
    return true;
  }
  return false;
}

void DE_LOptimizer::error_optimized_quality_bounded(
  LLocalOperation* local_op,
  const double quality_bound,
  const double current_error,
  Vec3d& opt_pos,
  double& error_before,
  double& error_after,
  double& quality_before,
  double& quality_after)
{
  /****** Initialize DE ******/
  // declare populations
  std::vector<Vec3d> old_population, mutation_population, crossover_population, new_population;
  // declare errors
  error_after = DBL_MAX;
  std::vector<double> errors_after(Np, DBL_MAX);
  size_t new_min_error_idx = 0;
  double new_min_error = DBL_MAX, old_min_error = DBL_MAX;
  // declare qualityies
  quality_after = DBL_MAX;
  std::vector<double> qualities_after(Np, DBL_MAX);
  // define functions
  auto is_quality_bounded = [&](double _quality_after) -> bool
    {
      return _quality_after <= quality_bound;
    };
  auto is_quality_improved = [&](double _quality_after) -> bool
    {
      return _quality_after <= quality_before;
    };

  // generate initial population
  BoundingBox bbox = local_op->bounding_box();
  Vec3d variable = local_op->get_variable();

  old_population = local_op->get_init_population(init_generator, init_dis, Np);
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
  error_before = local_op->error_before();
  quality_before = local_op->quality_before();

  //#pragma omp parallel for
  for (int i = 0;i < Np;i++)
  {
    if (local_op->update_variable(old_population[i], 0))
    {
      qualities_after[i] = local_op->quality_after(0);
      errors_after[i] = (is_quality_bounded(qualities_after[i]) || is_quality_improved(qualities_after[i])) ?
        local_op->error_after(0) : DBL_MAX;
    }
    else
    {
      errors_after[i] = DBL_MAX;
      qualities_after[i] = DBL_MAX;
    }
  }

  new_min_error_idx = std::distance(errors_after.begin(),
    std::min_element(errors_after.begin(), errors_after.end()));
  new_min_error = errors_after[new_min_error_idx];
  old_min_error = new_min_error;
  quality_after = qualities_after[new_min_error_idx];

  // enter loop
  mutation_population = old_population;
  crossover_population = old_population;
  size_t iter = 0, consecutive_iter = 0;
  while (true)
  {
    iter++;
    // generate mutation population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
      mutate(bbox, old_population, new_min_error_idx, i, mutation_population[i]);
    // generate crossover population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
      crossover(old_population[i], mutation_population[i], crossover_population[i]);
    // generate new population
  //#pragma omp parallel for
    for (int i = 0;i < Np;i++)
    {
      double new_quality = local_op->update_variable(crossover_population[i], 0) ? local_op->quality_after(0) : DBL_MAX;
      if (is_quality_bounded(new_quality) || is_quality_improved(new_quality))
      {
        double new_error = local_op->error_after(0);
        if (new_error < errors_after[i])
        {
          new_population[i] = crossover_population[i];  // swap faster?
          errors_after[i] = new_error;
          qualities_after[i] = new_quality;
        }
      }
    }
    new_min_error_idx = std::distance(errors_after.begin(),
      std::min_element(errors_after.begin(), errors_after.end()));
    new_min_error = errors_after[new_min_error_idx];
    quality_after = qualities_after[new_min_error_idx];
    // the first end condition : check if loop is convergent
    if (convergence(old_min_error, new_min_error, consecutive_iter))
      break;
    // the second end condition : check if loop reaches the maximum iteration number.
    if (iter == max_DE_iter)
    {
      //Logger::dev_logger->trace("maximum iteration number.");
      break;
    }
    // continue to optimize
    old_min_error = new_min_error;
    old_population = new_population;
  }

  opt_pos = new_population[new_min_error_idx];
  error_after = new_min_error;
}

bool DE_LOptimizer::do_eoqb(
  LLocalOperation* local_op,
  double quality_bound,
  double& current_error)
{
  Vec3d opt_pos;
  double error_before, error_after, quality_before, quality_after;
  error_optimized_quality_bounded(
    local_op, quality_bound, current_error,
    opt_pos, error_before, error_after, quality_before, quality_after);

  if (error_after < error_before)
  {
    local_op->just_do_it(opt_pos);
    current_error = current_error - error_before + error_after;
    Logger::dev_logger->trace("current error is {:.5f}", current_error);
    return true;
  }
  return false;
}
}// namespace ImageTriSimp
}// namespace GCLF