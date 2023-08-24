#include "ColorMap.h"

namespace GCLF
{
namespace ImageTriSimp
{

ColorMap::ColorMap()
{
  maxValue = 1.0;
  minValue = 0.0;
}

ColorMap::ColorMap(double maxV, double minV)
{
  maxValue = maxV;
  minValue = minV;
}

double ColorMap::GetMaxValue()
{
  return maxValue;
}

double ColorMap::GetMinValue()
{
  return minValue;
}

void ColorMap::SetMaxMinValue(double maxV, double minV)
{
  maxValue = maxV;
  minValue = minV;
}

Vec3d ColorMap::MapToColor(double value)
{
  double mapValue = (value - minValue) / (maxValue - minValue);
  if (mapValue < 0)
  {
    mapValue = 0;
  }
  else if (mapValue > 1)
  {
    mapValue = 1;
  }

  Vec3d result;
  if (mapValue < 0.25)
  {
    result = Vec3d(0, mapValue * 4, 1);   // blue to green
  }
  else if (mapValue < 0.5)
  {
    result = Vec3d(0, 1, 1 - 4 * (mapValue - 0.25));  // blue to green
  }
  else if (mapValue < 0.75)
  {
    result = Vec3d(4 * (mapValue - 0.5), 1, 0);   // green to red
  }
  else
  {
    result = Vec3d(1, 1 - 4 * (mapValue - 0.25), 0);  // green to red
  }
  result.minimize(Vec3d(1., 1., 1.));
  result.maximize(Vec3d(0., 0., 0.));
  return result;
}

}// namespace ImageTriSimp
}// namespace GCLF