#include "ParamParser.h"

namespace GCLF
{
namespace ImageTriSimp
{

ParamLinear::ParamLinear()
{
  min_angle = 10.;
  split_min_angle = 5.;

  max_simplify_iter = 20;
  convergence_collapse_number = 10;
}


ParamCurved::ParamCurved()
{
  min_angle = 10.;
  max_curvature = 135.;

  max_simplify_iter = 20;
  convergence_collapse_number = 5;
}


ParamTriangulator::ParamTriangulator()
{
  // color type
  color_type = "constant";
  // error type
  error_type = "l2";
  // DE
  Np = 20;
  mutation_scale = 0.7;
  crossover_rate = 0.9;
  convergence_rate = 1e-2;
  max_consecutive_iter = 5;
  max_DE_iter = 30;
  // Collapse
  max_valence = 10;
}

}// namespace ImageTriSimp
}// namespace GCLF