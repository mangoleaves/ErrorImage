#pragma once
#include "Image/ImageDefinition.h"
#include "Mesh/LinearMeshDefinition.h"
#include "Mesh/BezierMeshDefinition.h"
#include "Mesh/OMUtils.h"
#include "HLBFGS.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
class CalcErrorOnFace
{
public:
  enum class Type
  {
    L1Norm,
    L2Norm,
    L4Norm
  }type;
public:
  CalcErrorOnFace() = delete;
  CalcErrorOnFace(Simplifier* _simplifier) :simplifier(_simplifier) {}

  double operator()(FaceHandle fh);
  double operator()(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords);

  double calc_error_with_lpnorm(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords, double p);
  std::vector<double> calc_gradient_with_lpnorm(const FaceColor& face_color, const std::vector<Vec2i>& pixel_coords, double p);

  void preset_color_and_pixels(FaceColor* _face_color, const std::vector<Vec2i>* _pixel_coords)
  {
    ps_face_color = _face_color;
    ps_pixel_coords = _pixel_coords;
  }
  void unset_color_and_pixels()
  {
    ps_face_color = nullptr;
    ps_pixel_coords = nullptr;
  }
  void operator()(int N, double*, double*, double*, double*);
private:
  Simplifier* simplifier;
  // preset color and pixel coords, used for HLBFGS optimizer.
  FaceColor* ps_face_color;
  const std::vector<Vec2i>* ps_pixel_coords;
};

}// namespace ImageTriSimp
}// namespace GCLF

#include "CalcErrorOnFace_impl.h"