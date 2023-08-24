#include "ImageDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{

void ImageT::resize(int w, int h)
{
  width = w;
  height = h;
  pixel_colors.resize(h, w);
}

/// @brief reset mask to zero.
/// @param mask An integer matrix that has the same shape with image.
void ImageT::reset_mask(Eigen::MatrixXi& mask)
{
  if (mask.rows() != height || mask.cols() != width)
    mask.resize(height, width);

  mask.setZero();
}

ImageT::Color FaceColor::get_color(double x, double y)const
{
  switch (type)
  {
  case Type::Constant:
    return cc;
    break;
  case Type::Linear:
    return lc_x * x + lc_y * y + cc;
    break;
  case Type::Quadratic:
    return qc_xx * (x * x) + qc_xy * (x * y) + qc_yy * (y * y)
      + lc_x * x + lc_y * y + cc;
    break;
  default:
    return cc;
    break;
  }
}

ImageT::Color FaceColor::get_constant_color()const
{
  return cc;
}

ImageT::Color FaceColor::get_linear_color(double x, double y)const
{
  return lc_x * x + lc_y * y + cc;
}

ImageT::Color FaceColor::get_quadratic_color(double x, double y)const
{
  return qc_xx * (x * x) + qc_xy * (x * y) + qc_yy * (y * y)
    + lc_x * x + lc_y * y + cc;
}

}// namespace ImageTriSimp
}// namespace GCLF