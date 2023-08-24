#pragma once
#include "Basic/SimpleTypes.h"

namespace GCLF
{
using namespace Geometry;
namespace ImageTriSimp
{

class ColorMap
{
private:
  double maxValue, minValue;
public:
  ColorMap();
  ColorMap(double maxV, double minV);

  double GetMaxValue();
  double GetMinValue();
  void SetMaxMinValue(double maxV, double minV);

  Vec3d MapToColor(double value);
};

}// namespace ImageTriSimp
}// namespace GCLF