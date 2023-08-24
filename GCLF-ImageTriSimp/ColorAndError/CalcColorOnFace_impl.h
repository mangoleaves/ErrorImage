#pragma once
#include "CalcColorOnFace.h"
#include "Eigen/QR"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::operator()(FaceHandle fh)
{
  return operator()(simplifier->mesh->data(fh).pixels);
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::operator()(const std::vector<Vec2i>& pixel_coords)
{
  if (pixel_coords.empty())
    return FaceColor();

  typedef CalcErrorOnFace<Simplifier>::Type ErrorType;

  auto& calc_error_on_face = simplifier->behavior.calc_error_on_face;

  switch (type)
  {
  case Type::Constant:
    switch (calc_error_on_face.type)
    {
    case ErrorType::L1Norm:
      return calc_constant_color_l1norm(pixel_coords);
    case ErrorType::L2Norm:
      return calc_constant_color_l2norm(pixel_coords);
    case ErrorType::L4Norm:
      return calc_constant_color_l4norm(pixel_coords);
    }
  case Type::Linear:
    switch (calc_error_on_face.type)
    {
    case ErrorType::L1Norm:
      return calc_linear_color_l1norm(pixel_coords);
    case ErrorType::L2Norm:
      return calc_linear_color_l2norm(pixel_coords);
    case ErrorType::L4Norm:
      return calc_linear_color_l4norm(pixel_coords);
    }
  case Type::Quadratic:
    switch (calc_error_on_face.type)
    {
    case ErrorType::L1Norm:
      return calc_quadratic_color_l1norm(pixel_coords);
    case ErrorType::L2Norm:
      return calc_quadratic_color_l2norm(pixel_coords);
    case ErrorType::L4Norm:
      return calc_quadratic_color_l4norm(pixel_coords);
    }
    break;
  }
  return FaceColor();
}

// TODO: below codes can be highly improved

template<typename Simplifier>
int CalcColorOnFace<Simplifier>::sigma_fun(int iter)
{
  if (iter < 20) return 3;
  else if (iter < 120) return 5;
  else if (iter < 250) return 10;
  else if (iter < 500) return 50;
  else return 100;
}

template<typename Simplifier>
Eigen::MatrixXd CalcColorOnFace<Simplifier>::sign_fun(Eigen::MatrixXd x)
{
  for (int i = 0; i < x.rows(); i++)
  {
    for (int j = 0; j < x.cols(); j++)
    {
      if (x(i, j) > 0)
      {
        x(i, j) = 1;
      }
      else if (x(i, j) < 0)
      {
        x(i, j) = -1;
      }
      else
      {
        x(i, j) = 0;
      }
    }
  }
  return x;
}

template<typename Simplifier>
Eigen::MatrixXd CalcColorOnFace<Simplifier>::get_fun(Eigen::MatrixXd a, Eigen::MatrixXd b)
{
  int nrow = a.rows();
  int ncol = a.cols();
  Eigen::MatrixXd ans = Eigen::MatrixXd::Ones(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      ans(i, j) = a(i, j) * b(i, j);
    }
  }
  return ans;
}

template<typename Simplifier>
Eigen::MatrixXd CalcColorOnFace<Simplifier>::inverse_fun(Eigen::MatrixXd a)
{
  int nrow = a.rows();
  int ncol = a.cols();
  Eigen::MatrixXd ans = Eigen::MatrixXd::Ones(nrow, ncol);
  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      ans(i, j) = 1. / a(i, j);
    }
  }
  return ans;
}


template<typename Simplifier>
Eigen::MatrixXd CalcColorOnFace<Simplifier>::get_vector_diag(Eigen::MatrixXd a)
{
  int nrow = a.rows();
  int ncol = a.cols();
  int n = max(nrow, ncol);
  Eigen::MatrixXd ans = Eigen::MatrixXd::Ones(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      ans(i, j) = a(i);
    }
  }
  return ans;
}

template<typename Simplifier>
Eigen::MatrixXd CalcColorOnFace<Simplifier>::Prox_L1(Eigen::MatrixXd x, double lambda)
{
  int nrow = x.rows();
  int ncol = x.cols();
  Eigen::MatrixXd maty = Eigen::MatrixXd::Ones(nrow, 1);
  maty = maty * lambda;
  Eigen::MatrixXd absx = x.array().abs();
  Eigen::MatrixXd absx_lambda = absx - maty;

  for (int i = 0; i < nrow; i++)
  {
    for (int j = 0; j < ncol; j++)
    {
      absx_lambda(i, j) = max(absx_lambda(i, j), 0.);
    }
  }
  Eigen::MatrixXd signx = sign_fun(x);
  Eigen::MatrixXd xp = get_fun(signx, absx_lambda);//signx.dot(absx - maty);
  return xp;
}

template<typename Simplifier>
void CalcColorOnFace<Simplifier>::admm(Eigen::MatrixXd& x, Eigen::MatrixXd s, Eigen::MatrixXd lambda, Eigen::MatrixXd A, double normA, Eigen::MatrixXd b, int maxiter, double tol, double sigma)
{
  int n = A.rows();
  int p = A.cols();
  double tau = 1.618;
  double prim_win = 0;
  double dual_win = 0;
  Eigen::MatrixXd ATA = A.transpose() * A;
  for (int iter = 0; iter < maxiter; iter++)
  {
    double tau_sig = tau * sigma;
    Eigen::MatrixXd lambda_sig = lambda / sigma;
    double inv_sig = 1 / sigma;
    /*********************** to solve x****************************/

    x = ATA.colPivHouseholderQr().solve((A.transpose() * (s + b - lambda_sig)));

    //cout << x << "\n";
    Eigen::MatrixXd Ax = A * x;

    /********************to solve s * ***************************/

    s = Prox_L1(Ax - b + lambda_sig, inv_sig);

    /*********************update lambda * *****************************/

    Eigen::MatrixXd lambda_old = lambda;
    lambda = lambda_old + tau_sig * (Ax - b - s);
    Eigen::MatrixXd dlambda = lambda - lambda_old;
    /********************** Optimality condition*********************/

    double	pinf = sqrt((dlambda / tau_sig).squaredNorm()) / (0.01 + normA);
    double	dinf = sqrt((A.transpose() * lambda).squaredNorm()) / (0.01 + normA);
    double	ratio = pinf / dinf;

    //double cons_h = sqrt((s - Prox_L1(s + lambda, inv_sig)).squaredNorm());

    double opt_measure = sqrt(pow(pinf, 2) + pow(dinf, 2));
    if (opt_measure <= tol)// || isnan(opt_measure));   
    {
      //cout << opt_measure << " <= " << tol << "\n";
      break;
    }
    /****************** to update the value of sigma****************/
    // so that pinfand dinf can reduce synchronously

    int sigma_update_iter = sigma_fun(iter);
    double sigmascale = 1.1;

    if (ratio < 1)
    {
      prim_win = prim_win + 1;
    }
    else if (ratio > 1)
    {
      dual_win = dual_win + 1;
    }

    if ((iter % sigma_update_iter) == 0)
    {
      double sigmamax = 1e6;
      double sigmamin = 1e-4;
      if (prim_win > max(1., 1.2 * dual_win))
      {
        prim_win = 0;
        sigma = max(sigmamin, sigma / sigmascale);
      }
      else if (dual_win > max(1., 1.2 * prim_win))
      {
        dual_win = 0;
        sigma = min(sigmamax, sigma * sigmascale);
      }
    }

  }
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_constant_color_l1norm(const std::vector<Vec2i>& pixel_coords)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  FaceColor color;
  color.type = FaceColor::Type::Constant;

  vector<double> color_r;
  vector<double> color_g;
  vector<double> color_b;

  vector<double> r_solve;
  vector<double> g_solve;
  vector<double> b_solve;
  r_solve.resize(256);
  g_solve.resize(256);
  b_solve.resize(256);

  for (const auto& pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
  #if defined(CHANNEL_0_1)
    double realr = pixel_color.x() * 255;
    double realg = pixel_color.y() * 255;
    double realb = pixel_color.z() * 255;
  #elif defined(CHANNEL_0_255)
    double& realr = pixel_color.x();
    double& realg = pixel_color.y();
    double& realb = pixel_color.z();
  #endif

    color_r.push_back(realr);
    color_g.push_back(realg);
    color_b.push_back(realb);

    for (int t = 0; t < 256; t++)
    {
      r_solve[t] = r_solve[t] + abs((t - realr));
      g_solve[t] = g_solve[t] + abs((t - realg));
      b_solve[t] = b_solve[t] + abs((t - realb));
    }
  }
#if defined(CHANNEL_0_1)
  double ansR = (min_element(r_solve.begin(), r_solve.end()) - r_solve.begin()) / 255.;
  double ansG = (min_element(g_solve.begin(), g_solve.end()) - g_solve.begin()) / 255.;
  double ansB = (min_element(b_solve.begin(), b_solve.end()) - b_solve.begin()) / 255.;
#elif defined(CHANNEL_0_255)
  double ansR = (min_element(r_solve.begin(), r_solve.end()) - r_solve.begin());
  double ansG = (min_element(g_solve.begin(), g_solve.end()) - g_solve.begin());
  double ansB = (min_element(b_solve.begin(), b_solve.end()) - b_solve.begin());
#endif

  color.cc.x() = ansR;
  color.cc.y() = ansG;
  color.cc.z() = ansB;

  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_linear_color_l1norm(const std::vector<Vec2i>& pixel_coords)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  FaceColor color;
  color.type = FaceColor::Type::Linear;

  int nrow = pixel_coords.size();

  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(nrow, 3);
  Eigen::MatrixXd b1 = Eigen::MatrixXd::Ones(nrow, 1);
  Eigen::MatrixXd b2 = Eigen::MatrixXd::Ones(nrow, 1);
  Eigen::MatrixXd b3 = Eigen::MatrixXd::Ones(nrow, 1);
  int idx = 0;
  for (const Vec2i& pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
    A(idx, 0) = pc.x();
    A(idx, 1) = pc.y();
    A(idx, 2) = 1;
    b1(idx) = pixel_color.x();
    b2(idx) = pixel_color.y();
    b3(idx) = pixel_color.z();
    idx = idx + 1;
  }
  int ncol = A.cols();
  Eigen::MatrixXd x1 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd x2 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd x3 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd lambda0 = Eigen::MatrixXd::Zero(nrow, 1);
  int maxiter = 20000;
  double tol = 1e-5;
  double sigma = 1;
  double normA = sqrt(A.squaredNorm());
  //cout << normA << "\n";
  Eigen::MatrixXd s1 = A * x1 - b1;
  Eigen::MatrixXd s2 = A * x2 - b2;
  Eigen::MatrixXd s3 = A * x3 - b3;
  admm(x1, s1, lambda0, A, normA, b1, maxiter, tol, sigma);
  admm(x2, s2, lambda0, A, normA, b2, maxiter, tol, sigma);
  admm(x3, s3, lambda0, A, normA, b3, maxiter, tol, sigma);

  // set r channel
  color.lc_x.x() = x1(0); color.lc_y.x() = x1(1); color.cc.x() = x1(2);
  // set g channel
  color.lc_x.y() = x2(0); color.lc_y.y() = x2(1); color.cc.y() = x2(2);
  // set b channel
  color.lc_x.z() = x3(0); color.lc_y.z() = x3(1); color.cc.z() = x3(2);

  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_quadratic_color_l1norm(const std::vector<Vec2i>& pixel_coords)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  FaceColor color;
  color.type = FaceColor::Type::Quadratic;
  int nrow = pixel_coords.size();
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(nrow, 6);
  Eigen::MatrixXd b1 = Eigen::MatrixXd::Ones(nrow, 1);
  Eigen::MatrixXd b2 = Eigen::MatrixXd::Ones(nrow, 1);
  Eigen::MatrixXd b3 = Eigen::MatrixXd::Ones(nrow, 1);

  int idx = 0;
  for (const Vec2i& pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
    A(idx, 0) = pc.x() * pc.x();
    A(idx, 1) = pc.x() * pc.y();
    A(idx, 2) = pc.y() * pc.y();

    A(idx, 3) = pc.x();
    A(idx, 4) = pc.y();
    A(idx, 5) = 1;

    b1(idx) = pixel_color.x();
    b2(idx) = pixel_color.y();
    b3(idx) = pixel_color.z();
    idx = idx + 1;
  }
  int ncol = A.cols();
  Eigen::MatrixXd x1 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd x2 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd x3 = Eigen::MatrixXd::Ones(ncol, 1);
  Eigen::MatrixXd lambda0 = Eigen::MatrixXd::Zero(nrow, 1);
  int maxiter = 20000;
  double tol = 1e-5;
  double sigma = 1;
  double normA = sqrt(A.squaredNorm());
  Eigen::MatrixXd s1 = A * x1 - b1;
  Eigen::MatrixXd s2 = A * x2 - b2;
  Eigen::MatrixXd s3 = A * x3 - b3;
  admm(x1, s1, lambda0, A, normA, b1, maxiter, tol, sigma);
  admm(x2, s2, lambda0, A, normA, b2, maxiter, tol, sigma);
  admm(x3, s3, lambda0, A, normA, b3, maxiter, tol, sigma);

  // set r channel
  color.qc_xx.x() = x1(0); color.qc_xy.x() = x1(1); color.qc_yy.x() = x1(2);
  color.lc_x.x() = x1(3); color.lc_y.x() = x1(4); color.cc.x() = x1(5);
  // set g channel
  color.qc_xx.y() = x2(0); color.qc_xy.y() = x2(1); color.qc_yy.y() = x2(2);
  color.lc_x.y() = x2(3); color.lc_y.y() = x2(4); color.cc.y() = x2(5);
  // set b channel
  color.qc_xx.z() = x3(0); color.qc_xy.z() = x3(1); color.qc_yy.z() = x3(2);
  color.lc_x.z() = x3(3); color.lc_y.z() = x3(4); color.cc.z() = x3(5);

  return color;
}

// TODO end 

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_constant_color_l2norm(const std::vector<Vec2i>& pixel_coords)
{
  // Mesh could be linear mesh or Bezier mesh.
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  // Calculate a constant color on a face to represent a set of pixels.
  // The constant color is the average of colors of pixels.
  FaceColor color;
  color.type = FaceColor::Type::Constant;

  // accumulate pixels' color.
  ImageT::Color sum_color(0., 0., 0.);
  for (const Vec2i& pixel_coord : pixel_coords)
    sum_color += image->pixel_color(pixel_coord.x(), pixel_coord.y());
  // calculate average of pixels' color.
  color.cc = sum_color / pixel_coords.size();

  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_linear_color_l2norm(const std::vector<Vec2i>& pixel_coords)
{
  // Mesh could be linear mesh or Bezier mesh.
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  FaceColor color;
  color.type = FaceColor::Type::Linear;
  // Linear color f(x,y) = d*x + e*y + f.
  // For three channel, the coefficients are dr, er, fr, dg, eg, fg, db, eb, fb.
  // Solve three linear systems Ax = b to get coefficients on red, green and blue channels.
  Eigen::Matrix3d A;
  Eigen::Vector3d br, bg, bb;   // red, green, blue

  double xx = 0, xy = 0, yy = 0, x = 0, y = 0, c = pixel_coords.size();
  double brd = 0, bre = 0, brf = 0, bgd = 0, bge = 0, bgf = 0, bbd = 0, bbe = 0, bbf = 0;
  for (const Vec2i& pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
    // accumulate A
    xx += pc.x() * pc.x(); xy += pc.x() * pc.y(); yy += pc.y() * pc.y();
    x += pc.x(); y += pc.y();
    // accumulate b
    brd += pixel_color.x() * pc.x();
    bre += pixel_color.x() * pc.y();
    brf += pixel_color.x();
    bgd += pixel_color.y() * pc.x();
    bge += pixel_color.y() * pc.y();
    bgf += pixel_color.y();
    bbd += pixel_color.z() * pc.x();
    bbe += pixel_color.z() * pc.y();
    bbf += pixel_color.z();
  }
  A << xx, xy, x, xy, yy, y, x, y, c;
  br << brd, bre, brf;
  bg << bgd, bge, bgf;
  bb << bbd, bbe, bbf;
  Eigen::ColPivHouseholderQR<Eigen::Matrix3d> dec(A);
  Eigen::Vector3d res_r = dec.solve(br);   // dr, er, fr
  Eigen::Vector3d res_g = dec.solve(bg);   // dg, eg, fg
  Eigen::Vector3d res_b = dec.solve(bb);   // db, eb, fb

  // set r channel
  color.lc_x.x() = res_r.x(); color.lc_y.x() = res_r.y(); color.cc.x() = res_r.z();
  // set g channel
  color.lc_x.y() = res_g.x(); color.lc_y.y() = res_g.y(); color.cc.y() = res_g.z();
  // set b channel
  color.lc_x.z() = res_b.x(); color.lc_y.z() = res_b.y(); color.cc.z() = res_b.z();
  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_quadratic_color_l2norm(const std::vector<Vec2i>& pixel_coords)
{
  // Mesh could be linear mesh or Bezier mesh.
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  FaceColor color;
  color.type = FaceColor::Type::Quadratic;
  // Linear color f(x,y) = d*x + e*y + f.
  // For three channel, the coefficients are dr, er, fr, dg, eg, fg, db, eb, fb.
  // Solve three linear systems Ax = b to get coefficients on red, green and blue channels.
  Eigen::MatrixXd A;
  A.resize(6, 6); A.setZero();
  Eigen::VectorXd br, bg, bb;   // red, green, blue
  br.resize(6); bg.resize(6); bb.resize(6);
  br.setZero(); bg.setZero(); bb.setZero();

  // double xx = 0, xy = 0, yy = 0, x = 0, y = 0, c = 1.;
  double xnyn[6];

  for (const Vec2i& pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
    // accumulate A
    xnyn[0] = pc.x() * pc.x(); xnyn[1] = pc.x() * pc.y();
    xnyn[2] = pc.y() * pc.y(); xnyn[3] = pc.x();
    xnyn[4] = pc.y();          xnyn[5] = 1.;
    for (int i = 0;i < 6;i++)
      for (int j = 0;j < 6;j++)
        A(i, j) += xnyn[i] * xnyn[j];
    // accumulate b
    for (int i = 0;i < 6;i++)
    {
      br(i) += pixel_color.x() * xnyn[i];
      bg(i) += pixel_color.y() * xnyn[i];
      bb(i) += pixel_color.z() * xnyn[i];
    }
  }
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(A);
  Eigen::VectorXd res_r = dec.solve(br);   // dr, er, fr
  Eigen::VectorXd res_g = dec.solve(bg);   // dg, eg, fg
  Eigen::VectorXd res_b = dec.solve(bb);   // db, eb, fb

  // set r channel
  color.qc_xx.x() = res_r(0); color.qc_xy.x() = res_r(1); color.qc_yy.x() = res_r(2);
  color.lc_x.x() = res_r(3); color.lc_y.x() = res_r(4); color.cc.x() = res_r(5);
  // set g channel
  color.qc_xx.y() = res_g(0); color.qc_xy.y() = res_g(1); color.qc_yy.y() = res_g(2);
  color.lc_x.y() = res_g(3); color.lc_y.y() = res_g(4); color.cc.y() = res_g(5);
  // set b channel
  color.qc_xx.z() = res_b(0); color.qc_xy.z() = res_b(1); color.qc_yy.z() = res_b(2);
  color.lc_x.z() = res_b(3); color.lc_y.z() = res_b(4); color.cc.z() = res_b(5);
  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_constant_color_l4norm(const std::vector<Vec2i>& pixel_coords)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  FaceColor color;
  color.type = FaceColor::Type::Constant;

  vector<double> color_r;
  vector<double> color_g;
  vector<double> color_b;

  vector<double> r_solve;
  vector<double> g_solve;
  vector<double> b_solve;
  r_solve.resize(256);
  g_solve.resize(256);
  b_solve.resize(256);

  for (auto pc : pixel_coords)
  {
    ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
  #if defined(CHANNEL_0_1)
    double realr = pixel_color.x() * 255;
    double realg = pixel_color.y() * 255;
    double realb = pixel_color.z() * 255;
  #elif defined(CHANNEL_0_255)
    double& realr = pixel_color.x();
    double& realg = pixel_color.y();
    double& realb = pixel_color.z();
  #endif

    color_r.push_back(realr);
    color_g.push_back(realg);
    color_b.push_back(realb);

    for (int t = 0; t < 256; t++)
    {
      r_solve[t] = r_solve[t] + ((t - realr) * (t - realr) * (t - realr) * (t - realr));
      g_solve[t] = g_solve[t] + ((t - realg) * (t - realg) * (t - realg) * (t - realg));
      b_solve[t] = b_solve[t] + ((t - realb) * (t - realb) * (t - realb) * (t - realb));
    }
  }
#if defined(CHANNEL_0_1)
  double ansR = (min_element(r_solve.begin(), r_solve.end()) - r_solve.begin()) / 255.;
  double ansG = (min_element(g_solve.begin(), g_solve.end()) - g_solve.begin()) / 255.;
  double ansB = (min_element(b_solve.begin(), b_solve.end()) - b_solve.begin()) / 255.;
#elif defined(CHANNEL_0_255)
  double ansR = (min_element(r_solve.begin(), r_solve.end()) - r_solve.begin());
  double ansG = (min_element(g_solve.begin(), g_solve.end()) - g_solve.begin());
  double ansB = (min_element(b_solve.begin(), b_solve.end()) - b_solve.begin());
#endif

  color.cc.x() = ansR;
  color.cc.y() = ansG;
  color.cc.z() = ansB;

  return color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_linear_color_l4norm(const std::vector<Vec2i>& pixel_coords)
{
  // Mesh could be linear mesh or Bezier mesh.
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  auto& calc_error_on_face = simplifier->behavior.calc_error_on_face;

  // get an initial color for face. it's also the initial guess for optimization.
  FaceColor face_color = calc_linear_color_l2norm(pixel_coords);
  calc_error_on_face.preset_color_and_pixels(&face_color, &pixel_coords);
  double x[9];
  x[0] = face_color.lc_x.x(); x[1] = face_color.lc_y.x(); x[2] = face_color.cc.x();
  x[3] = face_color.lc_x.y(); x[4] = face_color.lc_y.y(); x[5] = face_color.cc.y();
  x[6] = face_color.lc_x.z(); x[7] = face_color.lc_y.z(); x[8] = face_color.cc.z();

  double parameter[20];
  int info[20];
  //initialize
  INIT_HLBFGS(parameter, info);
  parameter[6] = 1e-5;
  info[4] = 2000; // num_iter
  info[6] = 0;    // T
  info[7] = 0;    // with_hessian ? 1 : 0;
  info[10] = 0;
  info[11] = 1;

  auto newiteration = [&](int iter, int call_iter, double* x, double* f, double* g, double* gnorm)
    {
      if (iter % 200 == 0)
        std::cout << iter << ": " << call_iter << " " << *f << " " << *gnorm << std::endl;
    };
  double error_before = calc_error_on_face.calc_error_with_lpnorm(face_color, pixel_coords, 4.);
  HLBFGS_template(/*N*/9, /*M*/3, x, calc_error_on_face, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  double error_after = calc_error_on_face.calc_error_with_lpnorm(face_color, pixel_coords, 4.);
  double changed_percent = (error_after - error_before) / error_before * 100;
#ifdef OUTPUT_HLBFGS
  Logger::dev_logger->debug("error before {:.5f}, error after {:.5f}, change {:.2f}%", error_before, error_after, changed_percent);

  static size_t cnt = 0;
  static double percent_sum = 0.0;
  cnt++;
  percent_sum += changed_percent;
  if (cnt % 10 == 0)
    Logger::dev_logger->warn("change percent sum {:.5f}", percent_sum);
#endif

  return face_color;
}

template<typename Simplifier>
FaceColor CalcColorOnFace<Simplifier>::calc_quadratic_color_l4norm(const std::vector<Vec2i>& pixel_coords)
{
  // Mesh could be linear mesh or Bezier mesh.
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  auto& calc_error_on_face = simplifier->behavior.calc_error_on_face;

  // get an initial color for face. it's also the initial guess for optimization.
  FaceColor face_color = calc_quadratic_color_l2norm(pixel_coords);
  calc_error_on_face.preset_color_and_pixels(&face_color, &pixel_coords);
  double x[18];
  x[0] = face_color.qc_xx.x();  x[1] = face_color.qc_xy.x();  x[2] = face_color.qc_yy.x();
  x[3] = face_color.lc_x.x();  x[4] = face_color.lc_y.x();  x[5] = face_color.cc.x();
  x[6] = face_color.qc_xx.y();  x[7] = face_color.qc_xy.y();  x[8] = face_color.qc_yy.y();
  x[9] = face_color.lc_x.y();  x[10] = face_color.lc_y.y(); x[11] = face_color.cc.y();
  x[12] = face_color.qc_xx.z(); x[13] = face_color.qc_xy.z(); x[14] = face_color.qc_yy.z();
  x[15] = face_color.lc_x.z();  x[16] = face_color.lc_y.z(); x[17] = face_color.cc.z();

  double parameter[20];
  int info[20];
  //initialize
  INIT_HLBFGS(parameter, info);
  parameter[6] = 1e-9;
  info[4] = 2000; // num_iter
  info[6] = 0;    // T
  info[7] = 0;    // with_hessian ? 1 : 0;
  info[10] = 0;
  info[11] = 1;

  auto newiteration = [&](int iter, int call_iter, double* x, double* f, double* g, double* gnorm)
    {
      if (iter % 200 == 0)
        std::cout << iter << ": " << call_iter << " " << *f << " " << *gnorm << std::endl;
    };
  double error_before = calc_error_on_face.calc_error_with_lpnorm(face_color, pixel_coords, 4.);
  HLBFGS_template(/*N*/18, /*M*/3, x, calc_error_on_face, 0, HLBFGS_UPDATE_Hessian, newiteration, parameter, info);
  double error_after = calc_error_on_face.calc_error_with_lpnorm(face_color, pixel_coords, 4.);
  double changed_percent = (error_after - error_before) / error_before * 100;
#ifdef OUTPUT_HLBFGS
  Logger::dev_logger->debug("error before {:.5f}, error after {:.5f}, change {:.2f}%", error_before, error_after, changed_percent);

  static size_t cnt = 0;
  static double percent_sum = 0.0;
  cnt++;
  percent_sum += changed_percent;
  if (cnt % 10 == 0)
    Logger::dev_logger->warn("change percent sum {:.5f}", percent_sum);
#endif

  return face_color;
}

}// namespace ImageTriSimp
}// namespace GCLF