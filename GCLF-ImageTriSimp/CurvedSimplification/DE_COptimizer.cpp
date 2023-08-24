#include "DE_COptimizer.h"
#include "CurvedSimplifier.h"

namespace GCLF
{
namespace ImageTriSimp
{

void DE_COptimizer::initialize(CurvedSimplifier* _simplifier, ParamTriangulator& param)
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

void DE_COptimizer::init_random()
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

void DE_COptimizer::mutate(
  const BoundingBox& bbox,
  const std::vector<Points>& old_population,
  int best_idx, int old_idx,
  Points& mutation_individual)
{
  int rand1;
  do rand1 = mutate_dis(mutate_generator);
  while (rand1 == old_idx);
  const Points& op1 = old_population[rand1];
  int rand2;
  do rand2 = mutate_dis(mutate_generator);
  while (rand2 == old_idx || rand2 == rand1);
  const Points& op2 = old_population[rand2];
  // int rand3;
  // do rand3 = mutate_dis(mutate_generator);
  // while (rand3 == old_idx || rand3 == rand2 || rand3 == rand1);
  // const Points& op3 = old_population[rand3];
  size_t size = op1.size();

  for (size_t i = 0;i < size;i++)
  {
    mutation_individual[i] = old_population[old_idx][i] +
      mutation_scale * (old_population[best_idx][i] - old_population[old_idx][i]) +
      mutation_scale * (op1[i] - op2[i]);
    mutation_individual[i].maximize(bbox.min());
    mutation_individual[i].minimize(bbox.max());
  }
}

void DE_COptimizer::crossover(
  const Points& old_individual,
  const Points& mutation_individual,
  Points& crossover_individual)
{
  for (int i = 0;i < crossover_individual.size();i++)
  {
    crossover_individual[i].x() = cross_dis(cross_generator) <= crossover_rate ?
      mutation_individual[i].x() : old_individual[i].x();
    crossover_individual[i].y() = cross_dis(cross_generator) <= crossover_rate ?
      mutation_individual[i].y() : old_individual[i].y();
    // crossover_individual[i].z() = cross_dis(cross_generator) <= crossover_rate ?
    //   mutation_individual[i].z() : old_individual[i].z();
  }
}

bool DE_COptimizer::convergence(
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

void DE_COptimizer::error_bounded_quality_optimized(
  CLocalOperation* local_op,
  const double error_bound,
  const double current_error,
  Points& opt_pos,
  double& error_before,
  double& error_after,
  double& quality_before,
  double& quality_after)
{
  /****** Initialize DE ******/
  // declare populations
  std::vector<Points> old_population, mutation_population, crossover_population, new_population;
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
  // define boxes
  BoundingBox bbox = local_op->bounding_box();

  // generate initial population
  old_population = local_op->get_init_population(init_generator, init_dis, Np);
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
  error_before = local_op->error_before();
  quality_before = local_op->quality_before();

  //#pragma omp parallel for
  for (int i = 0;i < Np;i++)
  {
    if (local_op->update_variables(old_population[i], 0))
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
      double new_error = local_op->update_variables(crossover_population[i], 0) ?
        local_op->error_after(0) : DBL_MAX;
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

  opt_pos = std::move(new_population[new_min_quality_idx]);
  quality_after = new_min_quality;
}

bool DE_COptimizer::do_ebqo(
  CLocalOperation* local_op,
  double error_bound,
  double& current_error)
{
  Points opt_pos;
  double error_before, error_after, quality_before, quality_after;
  error_bounded_quality_optimized(
    local_op, error_bound, current_error,
    opt_pos, error_before, error_after, quality_before, quality_after
  );
  if ((current_error - error_before + error_after) <= error_bound && quality_after < quality_before)
  {
    local_op->just_do_it(opt_pos);
    current_error = current_error - error_before + error_after;
    Logger::dev_logger->trace("quality after is {:.5f}, current error is {:.5f}", quality_after, current_error);
    return true;
  }
  return false;
}

void DE_COptimizer::error_optimized_quality_bounded(
  CLocalOperation* local_op,
  const double min_angle,
  const double max_curvature,
  const double current_error,
  Points& opt_pos,
  double& error_before,
  double& error_after)
{
  /****** Initialize DE ******/
  // declare populations
  std::vector<Points> old_population, mutation_population, crossover_population, new_population;
  // declare errors
  error_after = DBL_MAX;
  std::vector<double> errors_after(Np, DBL_MAX);
  size_t new_min_error_idx = 0;
  double new_min_error = DBL_MAX, old_min_error = DBL_MAX;
  // define boxes
  BoundingBox bbox = local_op->bounding_box();

  // generate initial population
  old_population = local_op->get_init_population(init_generator, init_dis, Np);
  new_population = old_population;

  /****** Begin DE *******/

  // initialize
  error_before = local_op->error_before();

  //#pragma omp parallel for
  for (int i = 0;i < Np;i++)
  {
    if (local_op->update_variables(old_population[i], 0))
    {
      bool is_quality_bounded = local_op->is_quality_bounded(0, max_curvature, min_angle);
      errors_after[i] = is_quality_bounded ? local_op->error_after(0) : DBL_MAX;
    }
    else
    {
      errors_after[i] = DBL_MAX;
    }
  }

  new_min_error_idx = std::distance(errors_after.begin(),
    std::min_element(errors_after.begin(), errors_after.end()));
  new_min_error = errors_after[new_min_error_idx];
  old_min_error = new_min_error;

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
      bool is_quality_bounded = local_op->update_variables(crossover_population[i], 0) ? local_op->is_quality_bounded(0, max_curvature, min_angle) : false;
      if (is_quality_bounded)
      {
        double new_error = local_op->error_after(0);
        if (new_error < errors_after[i])
        {
          new_population[i] = crossover_population[i];  // swap faster?
          errors_after[i] = new_error;
        }
      }
    }
    new_min_error_idx = std::distance(errors_after.begin(),
      std::min_element(errors_after.begin(), errors_after.end()));
    new_min_error = errors_after[new_min_error_idx];
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

  error_after = errors_after[new_min_error_idx];
  opt_pos = new_population[new_min_error_idx];
}

bool DE_COptimizer::do_eoqb(
  CLocalOperation* local_op,
  double min_angle, double max_curvature,
  double& current_error)
{
  Points opt_pos;
  double error_before, error_after;
  error_optimized_quality_bounded(
    local_op, min_angle, max_curvature, current_error,
    opt_pos, error_before, error_after
  );

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