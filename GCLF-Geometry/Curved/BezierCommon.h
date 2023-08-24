#pragma once
#include <vector>
#include <iostream>
#include "Logger/Logger.h"
#include "Basic/SimpleTypes.h"
#include "Eigen/Core"

namespace GCLF
{
namespace Geometry
{
constexpr double SQRT_2 = 1.4142135623730950488;
constexpr double SQRT_3 = 1.7320508075688772935;

typedef int Index;
typedef Vec2i Index2;
typedef Vec3i Index3;

typedef std::vector<Index> Indices;
typedef std::vector<Vec3d> Points;

constexpr double bc_eps = 1e-6; // barycentric coordinates error tolerance.

}// namespace Geometry
}// namespace GCLF