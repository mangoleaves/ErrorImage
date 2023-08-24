#pragma once
#include <random>
#include "Mesh/LinearMeshDefinition.h"
#include "LocalOperations/LLocalOperation.h"

namespace GCLF
{
namespace ImageTriSimp
{

class LinearSimplifier;
struct ParamTriangulator;

class DE_LOptimizer
{
public:
  void initialize(LinearSimplifier* _simplifier, ParamTriangulator& param);

  void error_bounded_quality_optimized(
    LLocalOperation* local_op,
    const double error_bound,
    const double current_error,
    Vec3d& opt_pos,
    double& error_before,
    double& error_after,
    double& quality_before,
    double& quality_after
  );

  bool do_ebqo(LLocalOperation* local_op,
    double error_bound,
    double& current_error
  );

  void error_optimized_quality_bounded(
    LLocalOperation* local_op,
    const double quality_bound,
    const double current_error,
    Vec3d& opt_pos,
    double& error_before,
    double& error_after,
    double& quality_before,
    double& quality_after
  );

  bool do_eoqb(
    LLocalOperation* local_op,
    double quality_bound,
    double& current_error
  );
private:
  LinearSimplifier* simplifier;

  // configure parameters
  size_t Np;
  double mutation_scale;
  double crossover_rate;
  double convergence_rate;
  size_t max_consecutive_iter;
  size_t max_DE_iter;

  // random generators and distributions
  std::default_random_engine init_generator;
  std::uniform_real_distribution<double> init_dis;
  std::default_random_engine mutate_generator;
  std::uniform_int_distribution<int> mutate_dis;
  std::default_random_engine cross_generator;
  std::uniform_real_distribution<double> cross_dis;

  void init_random();

  void mutate(
    const BoundingBox& bbox,
    const std::vector<Vec3d>& old_population,
    int best_idx, int old_idx,
    Vec3d& mutation_individual
  );
  void crossover(
    const Vec3d& old_individual,
    const Vec3d& mutation_individual,
    Vec3d& crossover_individual
  );
  bool convergence(
    double old_min_error, double new_min_error,
    size_t& consecutive_iter
  );
};

}// namespace ImageTriSimp
}// namespace GCLF