#pragma once
#include "CalcErrorOnMesh.h"
#include <algorithm>

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
double CalcErrorOnMesh<Simplifier>::lp_error()
{
  auto mesh = simplifier->mesh.get();
  auto image = simplifier->image.get();
  auto out_image = simplifier->out_image.get();
  double error = 0.0;
  // accumulate error
  for (int y = 0;y < image->height;y++)
  {
    for (int x = 0;x < image->width;x++)
    {
      auto& pixel_color = image->pixel_color(x, y);
      auto& approx_color = out_image->pixel_color(x, y);
      auto color_diff = pixel_color - approx_color;
      color_diff.x() = fabs(color_diff.x()); color_diff.y() = fabs(color_diff.y()); color_diff.z() = fabs(color_diff.z());
      error += color_diff.pow(p).sum();
    }
  }
  return error;
}

template<typename Simplifier>
double CalcErrorOnMesh<Simplifier>::mean_lp_error()
{
  auto image = simplifier->image.get();
  return lp_error() / (image->width * image->height * 3);
}

template<typename Simplifier>
double CalcErrorOnMesh<Simplifier>::root_mean_lp_error()
{
  return std::pow(mean_lp_error(), 1. / p);
}

}// namespace ImageTriSimp
}// namespace GCLF