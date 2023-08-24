#pragma once
#include "Mesh/BezierMeshDefinition.h"
#include "Curved/BezierTriangleInjectivityChecker.h"

namespace GCLF
{
using namespace Utils;
namespace ImageTriSimp
{

bool check_injectivity(const BMeshT* mesh);

void initialize_uvw_samples(const uint32_t degree);

double MIPS_energy(const Points& ctrl_pnts, const uint32_t degree);

double approx_bezier_curvature(const BMeshT* mesh, EdgeHandle eh);

double approx_bezier_sector_angle(const BMeshT* mesh, HalfedgeHandle hh);

}// namespace ImageTriSimp
}// namespace GCLF