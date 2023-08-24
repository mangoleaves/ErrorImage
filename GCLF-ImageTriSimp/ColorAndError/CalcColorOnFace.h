#pragma once
#include "CalcErrorOnFace.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
class CalcColorOnFace
{
public:
  enum class Type
  {
    Constant,
    Linear,
    Quadratic
  }type;
public:
  CalcColorOnFace() = delete;
  CalcColorOnFace(Simplifier* _simplifier) :simplifier(_simplifier) {}

  FaceColor operator()(FaceHandle fh);
  FaceColor operator()(const std::vector<Vec2i>& pixel_coords);
private:
  Simplifier* simplifier;

  int sigma_fun(int iter);
  Eigen::MatrixXd sign_fun(Eigen::MatrixXd x);
  Eigen::MatrixXd get_fun(Eigen::MatrixXd a, Eigen::MatrixXd b);
  Eigen::MatrixXd inverse_fun(Eigen::MatrixXd a);
  Eigen::MatrixXd get_vector_diag(Eigen::MatrixXd a);
  Eigen::MatrixXd Prox_L1(Eigen::MatrixXd x, double lambda);
  void admm(Eigen::MatrixXd& x, Eigen::MatrixXd s, Eigen::MatrixXd lambda, Eigen::MatrixXd A, double normA, Eigen::MatrixXd b, int maxiter, double tol, double sigma);

  FaceColor calc_constant_color_l1norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_linear_color_l1norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_quadratic_color_l1norm(const std::vector<Vec2i>& pixel_coords);

  FaceColor calc_constant_color_l2norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_linear_color_l2norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_quadratic_color_l2norm(const std::vector<Vec2i>& pixel_coords);

  FaceColor calc_constant_color_l4norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_linear_color_l4norm(const std::vector<Vec2i>& pixel_coords);
  FaceColor calc_quadratic_color_l4norm(const std::vector<Vec2i>& pixel_coords);
};

}// namespace ImageTriSimp
}// namespace GCLF

#include "CalcColorOnFace_impl.h"