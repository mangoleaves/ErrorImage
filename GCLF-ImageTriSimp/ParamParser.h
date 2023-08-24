#pragma once
#include <string>
#include "boost/json.hpp"

namespace GCLF
{
namespace ImageTriSimp
{

struct ParamLinear
{
  /* Quality bound */
  double min_angle;         // the minimal angle bound
  double split_min_angle;   // the minimal angle bound used in edge splitting
  /* Initialization */

  /* Simplification */
  size_t max_simplify_iter;
  size_t convergence_collapse_number;

  ParamLinear();

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["min_angle"] = min_angle;
    jo["split_min_angle"] = split_min_angle;
    jo["max_simplify_iter"] = max_simplify_iter;
    jo["convergence_collapse_number"] = convergence_collapse_number;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    min_angle = jo.at("min_angle").as_double();
    split_min_angle = jo.at("split_min_angle").as_double();
    max_simplify_iter = jo.at("max_simplify_iter").as_int64();
    convergence_collapse_number = jo.at("convergence_collapse_number").as_int64();
  }
};

struct ParamCurved
{
  /* Quality bound */
  double min_angle;
  double max_curvature;
  /* Initialization */

  /* Simplification */
  size_t max_simplify_iter;
  size_t convergence_collapse_number;

  ParamCurved();

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["min_angle"] = min_angle;
    jo["max_curvature"] = max_curvature;
    jo["max_simplify_iter"] = max_simplify_iter;
    jo["convergence_collapse_number"] = convergence_collapse_number;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    min_angle = jo.at("min_angle").as_double();
    max_curvature = jo.at("max_curvature").as_double();
    max_simplify_iter = jo.at("max_simplify_iter").as_int64();
    convergence_collapse_number = jo.at("convergence_collapse_number").as_int64();
  }
};

struct ParamTriangulator
{
  /* Color type */
  //  "constant", "linear", "quadratic"
  std::string color_type;
  /* Error type */
  // "l1", "l2", "l4"
  std::string error_type;
  /* Error bound */
  double error_bound;       // in RMSE or other forms

  /* DE */
  size_t Np;
  double mutation_scale;
  double crossover_rate;
  double convergence_rate;
  size_t max_consecutive_iter;
  size_t max_DE_iter;

  /* Collapse */
  uint32_t max_valence;

  ParamLinear param_linear;
  ParamCurved param_curved;

  ParamTriangulator();

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["color_type"] = color_type;
    jo["error_type"] = error_type;
    jo["Np"] = Np;
    jo["mutation_scale"] = mutation_scale;
    jo["crossover_rate"] = crossover_rate;
    jo["convergence_rate"] = convergence_rate;
    jo["max_consecutive_iter"] = max_consecutive_iter;
    jo["max_DE_iter"] = max_DE_iter;
    jo["max_valence"] = max_valence;
    jo["paramLinear"] = param_linear.serialize();
    jo["paramCurved"] = param_curved.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    color_type = jo.at("color_type").as_string();
    error_type = jo.at("error_type").as_string();
    Np = jo.at("Np").as_int64();
    mutation_scale = jo.at("mutation_scale").as_double();
    crossover_rate = jo.at("crossover_rate").as_double();
    convergence_rate = jo.at("convergence_rate").as_double();
    max_consecutive_iter = jo.at("max_consecutive_iter").as_int64();
    max_DE_iter = jo.at("max_DE_iter").as_int64();
    max_valence = jo.at("max_valence").as_int64();
    param_linear.deserialize(jo.at("paramLinear").as_object());
    param_curved.deserialize(jo.at("paramCurved").as_object());
  }
};

}// namespace ImageTriSimp
}// namespace GCLF