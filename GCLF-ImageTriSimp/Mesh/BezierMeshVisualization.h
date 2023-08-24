#pragma once
#include "BezierMeshDefinition.h"
#include "LinearMeshDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{
namespace OM = OpenMesh;
using namespace Utils;
using namespace Geometry;

using OpenMesh::VertexHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::FaceHandle;

using OpenMesh::SmartVertexHandle;
using OpenMesh::SmartEdgeHandle;
using OpenMesh::SmartHalfedgeHandle;
using OpenMesh::SmartFaceHandle;

void linear_approximation(const BMeshT& bmesh, uint32_t density, LMeshT& lmesh);

void get_ctrl_mesh(const BMeshT& bmesh, LMeshT& cmesh, LMeshT& jmesh);

}// namespace ImageTriSimp
}// namespace GCLF