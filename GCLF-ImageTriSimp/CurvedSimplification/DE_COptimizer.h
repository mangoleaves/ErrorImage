#pragma once
#include <random>
#include "Mesh/BezierMeshDefinition.h"
#include "LocalOperations/CLocalOperation.h"

namespace GCLF
{
namespace ImageTriSimp
{

class CurvedSimplifier;
struct ParamTriangulator;

class DE_COptimizer
{
public:
  void initialize(CurvedSimplifier* _simplifier, ParamTriangulator& param);

  void error_bounded_quality_optimized(
    CLocalOperation* local_op,
    const double error_bound,
    const double current_error,
    Points& opt_pos,
    double& error_before,
    double& error_after,
    double& quality_before,
    double& quality_after
  );

  bool do_ebqo(
    CLocalOperation* local_op,
    double error_bound,
    double& current_error
  );

  void error_optimized_quality_bounded(
    CLocalOperation* local_op,
    const double min_angle,
    const double max_curvature,
    const double current_error,
    Points& opt_pos,
    double& error_before,
    double& error_after
  );

  bool do_eoqb(
    CLocalOperation* local_op,
    double min_angle, double max_curvature,
    double& current_error
  );
private:
  CurvedSimplifier* simplifier;

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
    const std::vector<Points>& old_population,
    int best_idx, int old_idx,
    Points& mutation_individual
  );
  void crossover(
    const Points& old_individual,
    const Points& mutation_individual,
    Points& crossover_individual
  );
  bool convergence(
    double old_min_error, double new_min_error,
    size_t& consecutive_iter
  );
};

}// namespace ImageTriSimp
}// namespace GCLF