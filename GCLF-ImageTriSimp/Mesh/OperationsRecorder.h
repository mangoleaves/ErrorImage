#pragma once
#include "BezierMeshDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{

class BezierMeshOperationsRecorder
{
public:
  enum class OperationType
  {
    Undefined,
    GarbageCollection,
    RelocateVertex,
    RelocateEdge,
    Flip,
    CollapseInteriorEdge,
    CollapseBoundaryEdge,
    CollapseBoundaryVertex
  };

  struct Operation
  {
    OperationType type;

    VertexHandle vh;
    EdgeHandle eh;
    HalfedgeHandle hh;
    Points points;

    Operation();
    Operation(const Operation& op)noexcept;
    Operation(Operation&& op)noexcept;
    Operation& operator=(const Operation& op)noexcept;
    Operation& operator=(Operation&& op)noexcept;
  };

  std::unique_ptr<BMeshT> start_mesh;
  std::unique_ptr<BMeshT> replay_mesh;
  std::vector<Operation> operations;
public:
  BezierMeshOperationsRecorder() = default;

  void set_start_mesh(const BMeshT& mesh) { start_mesh = std::make_unique<BMeshT>(mesh); operations.clear(); }

  void add_collapse_interior_edge(HalfedgeHandle collapse_hh, const Points& new_ctrlpnts);
  void add_collapse_boundary_edge(EdgeHandle collapse_eh, const Points& new_ctrlpnts);
  void add_collapse_boundary_vert(EdgeHandle collapse_eh, const Points& new_ctrlpnts);
  void add_flip(EdgeHandle flip_eh, const Points& new_ctrlpnts);
  void add_relocate_edge(EdgeHandle relocate_eh, const Points& new_ctrlpnts);
  void add_relocate_vertex(VertexHandle relocate_vh, const Points& new_ctrlpnts);
  void add_garbage_collection();

  void restart();
  void replay_one();

  void write_file(const std::string& mesh_file, const std::string& op_file);
  void read_file(const std::string& mesh_file, const std::string& op_file);
private:
  bool recording = false;
  size_t op_idx;
};

}// namespace ImageTriSimp
}// namespace GCLF