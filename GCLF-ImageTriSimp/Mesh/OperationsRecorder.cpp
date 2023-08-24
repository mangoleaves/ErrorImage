#pragma once
#include "OperationsRecorder.h"

namespace GCLF
{
namespace ImageTriSimp
{
typedef BezierMeshOperationsRecorder BOR;

BOR::Operation::Operation()
{
  type = OperationType::Undefined;
  vh.invalidate(); eh.invalidate(); hh.invalidate();
}
BOR::Operation::Operation(const Operation& op)noexcept
{
  type = op.type;
  vh = op.vh; eh = op.eh; hh = op.hh;
  points = op.points;
}
BOR::Operation::Operation(Operation&& op)noexcept
{
  type = op.type;
  vh = op.vh; eh = op.eh; hh = op.hh;
  points = std::move(op.points);
}
BOR::Operation& BOR::Operation::operator=(const Operation& op)noexcept
{
  type = op.type;
  vh = op.vh; eh = op.eh; hh = op.hh;
  points = op.points;
  return *this;
}
BOR::Operation& BOR::Operation::operator=(Operation&& op)noexcept
{
  type = op.type;
  vh = op.vh; eh = op.eh; hh = op.hh;
  points = std::move(op.points);
  return *this;
}

void BOR::add_collapse_interior_edge(HalfedgeHandle collapse_hh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::CollapseInteriorEdge;
  op.hh = collapse_hh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_collapse_boundary_edge(EdgeHandle collapse_eh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::CollapseBoundaryEdge;
  op.eh = collapse_eh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_collapse_boundary_vert(EdgeHandle collapse_eh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::CollapseBoundaryVertex;
  op.eh = collapse_eh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_flip(EdgeHandle flip_eh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::Flip;
  op.eh = flip_eh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_relocate_edge(EdgeHandle relocate_eh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::Flip;
  op.eh = relocate_eh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_relocate_vertex(VertexHandle relocate_vh, const Points& new_ctrlpnts)
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::RelocateVertex;
  op.vh = relocate_vh;
  op.points = new_ctrlpnts;
  operations.push_back(std::move(op));
}

void BOR::add_garbage_collection()
{
  if (!recording)
    return;
  Operation op;
  op.type = OperationType::GarbageCollection;
  operations.push_back(std::move(op));
}

void BOR::restart()
{
  *replay_mesh = *start_mesh; op_idx = 0;
}

void BOR::replay_one()
{
  // TODO: how to visualize efficiently?
}

void BOR::write_file(const std::string& mesh_file, const std::string& op_file)
{
  // TODO: how to keep the mesh unchanged during IO?
}

void BOR::read_file(const std::string& mesh_file, const std::string& op_file)
{
  // TODO: how to keep the mesh unchanged during IO?
}

}// namespace ImageTriSimp
}// namespace GCLF