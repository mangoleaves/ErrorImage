#pragma once
#include "Eigen/Core"
#include "Basic/SimpleTypes.h"

namespace GCLF
{
namespace ImageTriSimp
{

/* Pixel:
  A pixel is defines as a square with length being 1.

  Upper left ---- Upper right
     |                 |
     |                 |
     |     Center      |
     |                 |
     |                 |
  Lower Left ---- Lower right

  Abbreviated as UL, UR, CT, LL, LR respectly.
*/

/* Image:
  An image is composed by pixels.
  The order of pixels is from left to right, from lower to upper.
  The coordinate of a pixel is the position that pixel's LL lies in.
  The first pixel lies in the origin (0, 0).

  (0, height) ------------ (width, height)
      |                          |
      |                          |
      |                          |
      |                          |
    (0, 0) ---------------  (width, 0)
*/

/* Matrix:
  When use a matrix to represent an image or something else,
  the matrix has the same shape with the image.
  The first row of matrix is the lowest row (y=0) of the image.
  The first column of matrix is the left most column (x=0) of the image.
  In another words, the matrix is the upside down of the image.
*/

/* Channel:
  float, range [0., 1.] or [0., 255.].
*/

/* Color:
  vector<Channel, 4>, representing red, green, blue, alpha channels.
  TODO: support more color space, e.g., HSV.
*/

// #define CHANNEL_0_1
#define CHANNEL_0_255

class ImageT
{
public:
  typedef double Channel;
  typedef Geometry::Vec3d Color;
public:
  ImageT() = default;

  void resize(int w, int h);

  Color& pixel_color(int x, int y) { return pixel_colors(y, x); }
  const Color& pixel_color(int x, int y) const { return pixel_colors(y, x); }

  Channel& pixel_red(int x, int y) { return pixel_colors(y, x).x(); }
  Channel& pixel_green(int x, int y) { return pixel_colors(y, x).y(); }
  Channel& pixel_blue(int x, int y) { return pixel_colors(y, x).z(); }
  const Channel& pixel_red(int x, int y) const{ return pixel_colors(y, x).x(); }
  const Channel& pixel_green(int x, int y) const{ return pixel_colors(y, x).y(); }
  const Channel& pixel_blue(int x, int y) const{ return pixel_colors(y, x).z(); }

  static Channel& red_channel(Color& c) { return c.x(); }
  static Channel& green_channel(Color& c) { return c.y(); }
  static Channel& blue_channel(Color& c) { return c.z(); }
  static const Channel& red_channel(const Color& c) { return c.x(); }
  static const Channel& green_channel(const Color& c) { return c.y(); }
  static const Channel& blue_channel(const Color& c) { return c.z(); }

  void reset_mask(Eigen::MatrixXi& mask);
public:
  int height, width;
  /// Pixel colors, each integer represents four channels [Alpha:Red:Green:Blue].
  /// Each channel is an 8 bits / 1 byte.
  Eigen::Matrix<Color, -1, -1> pixel_colors;
};

class FaceColor
{
public:
  enum class Type
  {
    Constant,   // f(x,y) = f
    Linear,     // f(x,y) = d*x + e*y + f
    Quadratic   // f(x,y) = a*x^2 + b*xy + c*y^2 + d*x + e*y + f
  }type;

  // constant coefficient
  ImageT::Color cc;
  // linear coefficient for x and y
  ImageT::Color lc_x, lc_y;
  // quadratic coefficient for x^2, xy, y^2
  ImageT::Color qc_xx, qc_xy, qc_yy;
public:
  ImageT::Color get_color(double x, double y)const;
  ImageT::Color get_constant_color()const;
  ImageT::Color get_linear_color(double x, double y)const;
  ImageT::Color get_quadratic_color(double x, double y)const;
};

}// namespace ImageTriSimp
}// namespace GCLF