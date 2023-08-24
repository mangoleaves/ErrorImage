#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriConnectivity.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include "Image/ImageDefinition.h"
#include "BezierFaceDefinition.h"
#include "Curved/BezierCurve.h"
#include "OMUtils.h"

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

// piecewise linear mesh.
class LMeshT;

class BMeshT: public OM::TriConnectivity
{
  /****************************/
  /*  Types and Definitions   */
  /****************************/
public:

  class VertexData
  {
  public:
    Index cpi;
  };
  class HalfedgeData
  {
  public:
  };
  class EdgeData
  {
  public:
    Index cpi;
  };
  class FaceData
  {
  public:
    BezierFace bface;

    FaceColor color;
    std::vector<Vec2i> pixels;
    double error;
  public:
    FaceData() = default;
    FaceData(const FaceData& rhs);
    FaceData(FaceData&& rhs)noexcept;
    FaceData& operator=(const FaceData& rhs);
    FaceData& operator=(FaceData&& rhs)noexcept;
  };
  /**********************************/
  /*  Constructors and Destructors  */
  /**********************************/
public:
  BMeshT();
  BMeshT(const BMeshT& rhs);
  BMeshT(BMeshT&& rhs)noexcept;
  BMeshT& operator=(const BMeshT& rhs);
  BMeshT& operator=(BMeshT&& rhs)noexcept;
  /*****************************/
  /*  Bezier Data & Functions  */
  /*****************************/
public:
  // Bezier mesh has a unique degree for all bezier faces.
  uint32_t bezier_degree;
  // Control points of all bezier faces are stored here.
  Points ctrl_pnt;
  std::vector<uint8_t> ctrl_pnt_deleted;

  /**********************/
  /* Global operations  */
  /**********************/

  Index add_ctrl_pnt(const Vec3d& p);
  void delete_ctrl_pnt(Index ctrl_pnt_idx);
  void clear_ctrl_pnt() { ctrl_pnt.clear(); ctrl_pnt_deleted.clear(); }

  void initialize(LMeshT* lmesh);

  void degree_elevation();

  void build_ctrlpnt_idx();

  void sync_local_ctrlpnt();

  void store_extra_cpi_to_ve();

  /**********************/
  /*  Local operations  */
  /**********************/
public:
  VertexHandle split_face(FaceHandle fh, Vec3d uvw);

  VertexHandle split_edge(HalfedgeHandle hh, Vec2d uv);

  bool collapse(HalfedgeHandle hh);

  // Flip:
  // We combine a split and then a collapse operation to do flip operation.

  // Flip only for the second order Bezier mesh.
  bool flip_2(EdgeHandle eh);

  // Relocate:
  // The relocate operation is dependent on the smooth algorithm.
  // The algorithm may influence a region of mesh, i.e., influence several control points.

  // Relocate only for the second order Bezier mesh.
  const Vec3d& point(VertexHandle vh) { return ctrl_pnt[data(vh).cpi]; }
  const Vec3d& point(VertexHandle vh) const { return ctrl_pnt[data(vh).cpi]; }
  const Vec3d& point(EdgeHandle eh) { return ctrl_pnt[data(eh).cpi]; }
  const Vec3d& point(EdgeHandle eh) const { return ctrl_pnt[data(eh).cpi]; }

  void set_point(VertexHandle vh, const Vec3d& p);
  void set_point(EdgeHandle eh, const Vec3d& p);

  Vec3d tangent_vec(HalfedgeHandle hh)const;

  BezierCurveImpl curve_on_edge(EdgeHandle eh)const;

  void local_mesh(const std::set<FaceHandle>& local_faces, BMeshT& local,
    std::map<VertexHandle, VertexHandle>& v2v, std::map<FaceHandle, FaceHandle>& f2f,
    std::vector<VertexHandle>& rv2v, std::vector<FaceHandle>& rf2f)const;

  /*********************/
  /*   Visualization   */
  /*********************/

  std::vector<std::vector<Vec3d>> linear_approx_curves(uint32_t density)const;

  /**********/
  /*   IO   */
  /**********/

  bool read_ascii(std::string file_name);
  bool write_ascii(std::string file_name);

  /******************/
  /*  Process Image */
  /******************/
  void scale_to(ImageT& image);

  /*********************/
  /*  Data management  */
  /*********************/
public:
  void delete_face(FaceHandle fh);

  void garbage_collection();

  void clean();
public:
  VertexData& data(VertexHandle _vh) { return this->property(data_vpph_, _vh); }
  const VertexData& data(VertexHandle _vh) const { return this->property(data_vpph_, _vh); }

  FaceData& data(FaceHandle _fh) { return this->property(data_fpph_, _fh); }
  const FaceData& data(FaceHandle _fh) const { return this->property(data_fpph_, _fh); }

  EdgeData& data(EdgeHandle _eh) { return this->property(data_epph_, _eh); }
  const EdgeData& data(EdgeHandle _eh) const { return this->property(data_epph_, _eh); }

  HalfedgeData& data(HalfedgeHandle _heh) { return this->property(data_hpph_, _heh); }
  const HalfedgeData& data(HalfedgeHandle _heh) const { return this->property(data_hpph_, _heh); }
private:
  typedef OM::VPropHandleT<VertexData>    DataVPropHandle;
  typedef OM::HPropHandleT<HalfedgeData>  DataHPropHandle;
  typedef OM::EPropHandleT<EdgeData>      DataEPropHandle;
  typedef OM::FPropHandleT<FaceData>      DataFPropHandle;

  DataVPropHandle                           data_vpph_;
  DataHPropHandle                           data_hpph_;
  DataEPropHandle                           data_epph_;
  DataFPropHandle                           data_fpph_;

  /***********************/
  /*  Internal functions */
  /***********************/
private:
  void copy_from(const BMeshT& rhs);
  void move_from(BMeshT&& rhs);
};

}// namespace ImageTriSimp
}// namespace GCLF