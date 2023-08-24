#pragma once
#include "MeshToImage.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
ImageT MeshToImage<Simplifier>::operator()()
{
  auto* mesh = simplifier->mesh.get();
  ImageT* image = simplifier->image.get();
  ImageT out_img;

  out_img.resize(image->width, image->height);

#if defined(CHANNEL_0_1)
  ImageT::Color max_color(1., 1., 1.), min_color(0., 0., 0.);
#else defined(CHANNEL_0_255)
  ImageT::Color max_color(255., 255., 255.), min_color(0., 0., 0.);
#endif

  for (FaceHandle fh : mesh->faces())
  {
    auto& face_color = mesh->data(fh).color;
    switch (face_color.type)
    {
    case FaceColor::Type::Constant:
    {
      for (const Vec2i& pc : mesh->data(fh).pixels)
      {
        ImageT::Color color = face_color.get_constant_color();
        color.minimize(max_color); color.maximize(min_color);
        out_img.pixel_color(pc.x(), pc.y()) = color;
      }
    }
    break;
    case FaceColor::Type::Linear:
    {
      for (const Vec2i& pc : mesh->data(fh).pixels)
      {
        ImageT::Color color = face_color.get_linear_color((double)pc.x(), (double)pc.y());
        color.minimize(max_color); color.maximize(min_color);
        out_img.pixel_color(pc.x(), pc.y()) = color;
      }
    }
    break;
    case FaceColor::Type::Quadratic:
    {
      for (const Vec2i& pc : mesh->data(fh).pixels)
      {
        ImageT::Color color = face_color.get_quadratic_color((double)pc.x(), (double)pc.y());
        color.minimize(max_color); color.maximize(min_color);
        out_img.pixel_color(pc.x(), pc.y()) = color;
      }
    }
    default:break;
    }
  }
  return out_img;
}

}// namespace ImageTriSimp
}// namespace GCLF