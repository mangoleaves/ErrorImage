#pragma once
#include "CalcErrorOnFace.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
double CalcErrorOnFace<Simplifier>::operator()(FaceHandle fh)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();

  return operator()(mesh->data(fh).color, mesh->data(fh).pixels);
}

template<typename Simplifier>
double CalcErrorOnFace<Simplifier>::operator()(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords)
{
  if (pixel_coords.empty())
    return DBL_MAX;

  switch (type)
  {
  case (Type::L1Norm):
    return calc_error_with_lpnorm(face_color, pixel_coords, 1);
  case (Type::L2Norm):
    return calc_error_with_lpnorm(face_color, pixel_coords, 2);
  case (Type::L4Norm):
    return calc_error_with_lpnorm(face_color, pixel_coords, 4);
  }
  return DBL_MAX;
}

template<typename Simplifier>
double CalcErrorOnFace<Simplifier>::calc_error_with_lpnorm(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords, double p)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  double error = 0.0;

  // NOTE:These branches for face color type could be hidden in color.get_color().
  // But I think explicit branches may be more efficient?
  switch (face_color.type)
  {
  case FaceColor::Type::Constant:
  {
    /// Calculate the error between face color and pixels' color.
    ImageT::Color approx_color = face_color.get_constant_color();
    for (const Vec2i& pixel_coord : pixel_coords)
    {
      ImageT::Color& pixel_color = image->pixel_color(pixel_coord.x(), pixel_coord.y());
      ImageT::Color color_diff = approx_color - pixel_color;
      color_diff._x = abs(color_diff._x); color_diff._y = abs(color_diff._y); color_diff._z = abs(color_diff._z);
      error += color_diff.pow(p).sum();
    }
  }
  break;
  case FaceColor::Type::Linear:
  {
    for (const Vec2i& pixel_coord : pixel_coords)
    {
      ImageT::Color& pixel_color = image->pixel_color(pixel_coord.x(), pixel_coord.y());
      ImageT::Color approx_color = face_color.get_linear_color((double)pixel_coord.x(), (double)pixel_coord.y());
      ImageT::Color color_diff = approx_color - pixel_color;
      color_diff._x = abs(color_diff._x); color_diff._y = abs(color_diff._y); color_diff._z = abs(color_diff._z);
      error += color_diff.pow(p).sum();
    }
  }
  break;
  case FaceColor::Type::Quadratic:
  {
    for (const Vec2i& pixel_coord : pixel_coords)
    {
      ImageT::Color& pixel_color = image->pixel_color(pixel_coord.x(), pixel_coord.y());
      ImageT::Color approx_color = face_color.get_quadratic_color((double)pixel_coord.x(), (double)pixel_coord.y());
      ImageT::Color color_diff = approx_color - pixel_color;
      color_diff._x = abs(color_diff._x); color_diff._y = abs(color_diff._y); color_diff._z = abs(color_diff._z);
      error += color_diff.pow(p).sum();
    }
  }
  break;
  }
  return error;
}

template<typename Simplifier>
std::vector<double> CalcErrorOnFace<Simplifier>::calc_gradient_with_lpnorm(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords, double p)
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  std::vector<double> gradient;

  // 3 color channel, 3 variables for linear f(x,y) = dx + ey + f.
  // d-r, e-r, f-r, d-g, e-g, f-g, d-b, e-b, f-g
  if (face_color.type == FaceColor::Type::Linear)
    gradient.resize(9, 0.0);
  // 3 color channel, 6 variables for quadratic f(x,y) = ax^2 + bxy + cy^2 + dx + ey + f.
  // a-r, b-r, c-r, d-r, e-r, f-r ...
  else if (face_color.type == FaceColor::Type::Quadratic)
    gradient.resize(18, 0.0);

  switch (face_color.type)
  {
  case FaceColor::Type::Linear:
  {
    for (const Vec2i& pixel_coord : pixel_coords)
    {
      ImageT::Color& pixel_color = image->pixel_color(pixel_coord.x(), pixel_coord.y());
      ImageT::Color approx_color = face_color.get_linear_color((double)pixel_coord.x(), (double)pixel_coord.y());
      ImageT::Color color_diff = approx_color - pixel_color;
      ImageT::Color color_diff_p = color_diff.pow(p - 1) * p;

      gradient[0] += color_diff_p.x() * pixel_coord.x();  // partial Error / partial d-r
      gradient[1] += color_diff_p.x() * pixel_coord.y();  // partial Error / partial e-r
      gradient[2] += color_diff_p.x();                    // partial Error / partial f-r
      gradient[3] += color_diff_p.y() * pixel_coord.x();  // d-g
      gradient[4] += color_diff_p.y() * pixel_coord.y();  // e-g
      gradient[5] += color_diff_p.y();                    // f-g
      gradient[6] += color_diff_p.z() * pixel_coord.x();  // d-b
      gradient[7] += color_diff_p.z() * pixel_coord.y();  // e-b
      gradient[8] += color_diff_p.z();                    // f-b
    }
  }
  break;
  case FaceColor::Type::Quadratic:
  {
    for (const Vec2i& pc : pixel_coords)
    {
      ImageT::Color& pixel_color = image->pixel_color(pc.x(), pc.y());
      ImageT::Color approx_color = face_color.get_quadratic_color((double)pc.x(), (double)pc.y());
      ImageT::Color color_diff = approx_color - pixel_color;
      ImageT::Color color_diff_p = color_diff.pow(p - 1) * p;

      double xx = pc.x() * pc.x(), xy = pc.x() * pc.y(), yy = pc.y() * pc.y();

      gradient[0] += color_diff_p.x() * xx;       // partial Error / partial a-r
      gradient[1] += color_diff_p.x() * xy;       // partial Error / partial b-r
      gradient[2] += color_diff_p.x() * yy;       // partial Error / partial c-r
      gradient[3] += color_diff_p.x() * pc.x();   // partial Error / partial d-r
      gradient[4] += color_diff_p.x() * pc.y();   // partial Error / partial e-r
      gradient[5] += color_diff_p.x();            // partial Error / partial f-r
      gradient[6] += color_diff_p.y() * xx;
      gradient[7] += color_diff_p.y() * xy;
      gradient[8] += color_diff_p.y() * yy;
      gradient[9] += color_diff_p.y() * pc.x();
      gradient[10] += color_diff_p.y() * pc.y();
      gradient[11] += color_diff_p.y();
      gradient[12] += color_diff_p.z() * xx;
      gradient[13] += color_diff_p.z() * xy;
      gradient[14] += color_diff_p.z() * yy;
      gradient[15] += color_diff_p.z() * pc.x();
      gradient[16] += color_diff_p.z() * pc.y();
      gradient[17] += color_diff_p.z();
    }
  }
  break;
  }
  return gradient;
}

template<typename Simplifier>
void CalcErrorOnFace<Simplifier>::operator()(int N, double* x, double* prev_x, double* F, double* G)
{
  ASSERT(ps_face_color && ps_pixel_coords, "null color and pixels.");
  ASSERT(type == Type::L4Norm, "This function is only used for L4 norm.");
  if (ps_face_color->type == FaceColor::Type::Linear)
  {
    ASSERT(N == 9, "wrong number of variables.");
  }
  else if (ps_face_color->type == FaceColor::Type::Quadratic)
  {
    ASSERT(N == 18, "wrong number of variables.");
  }

  switch (ps_face_color->type)
  {
  case FaceColor::Type::Linear:
  {
    ps_face_color->lc_x.x() = x[0]; ps_face_color->lc_y.x() = x[1]; ps_face_color->cc.x() = x[2];
    ps_face_color->lc_x.y() = x[3]; ps_face_color->lc_y.y() = x[4]; ps_face_color->cc.y() = x[5];
    ps_face_color->lc_x.z() = x[6]; ps_face_color->lc_y.z() = x[7]; ps_face_color->cc.z() = x[8];

    *F = calc_error_with_lpnorm(*ps_face_color, *ps_pixel_coords, 4.);
    std::vector<double> gradient = calc_gradient_with_lpnorm(*ps_face_color, *ps_pixel_coords, 4.);
    ASSERT(gradient.size() == 9, "wrong number of variables.");
    std::memcpy(G, gradient.data(), 9 * sizeof(double));
  }
  break;
  case FaceColor::Type::Quadratic:
  {
    ps_face_color->qc_xx.x() = x[0]; ps_face_color->qc_xy.x() = x[1]; ps_face_color->qc_yy.x() = x[2];
    ps_face_color->lc_x.x() = x[3]; ps_face_color->lc_y.x() = x[4]; ps_face_color->cc.x() = x[5];
    ps_face_color->qc_xx.y() = x[6]; ps_face_color->qc_xy.y() = x[7]; ps_face_color->qc_yy.y() = x[8];
    ps_face_color->lc_x.y() = x[9]; ps_face_color->lc_y.y() = x[10]; ps_face_color->cc.y() = x[11];
    ps_face_color->qc_xx.z() = x[12]; ps_face_color->qc_xy.z() = x[13]; ps_face_color->qc_yy.z() = x[14];
    ps_face_color->lc_x.z() = x[15]; ps_face_color->lc_y.z() = x[16]; ps_face_color->cc.z() = x[17];

    *F = calc_error_with_lpnorm(*ps_face_color, *ps_pixel_coords, 4.);
    std::vector<double> gradient = calc_gradient_with_lpnorm(*ps_face_color, *ps_pixel_coords, 4.);
    ASSERT(gradient.size() == 18, "wrong number of variables.");
    std::memcpy(G, gradient.data(), 18 * sizeof(double));
  }
  break;
  }
}

}// namespace ImageTriSimp
}// namespace GCLF