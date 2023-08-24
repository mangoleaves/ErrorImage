#pragma once
#include <tuple>
#include "Basic/Types.h"
#include "Curved/BezierTriangle.h"
#include "Logger/Logger.h"
#include "Utils/Traversal.h"

namespace GCLF
{
namespace ImageTriSimp
{
using namespace Utils;
using namespace Geometry;

// Forward declaration
class BezierFace;

/// @brief A basic class for iterator.
template<bool IsConstant, bool IsSequential>
class BF_Iterator
{
public:
  using BFacePtr = std::conditional_t<IsConstant, const BezierFace*, BezierFace*>;
  using IndicesPtr = const Indices*;
  using IterT = BF_Iterator<IsConstant, IsSequential>;
protected:
  BFacePtr bface;
  IndicesPtr indices;
  Index i;
public:
  BF_Iterator():indices(nullptr), i(-1) {}
  BF_Iterator(const BF_Iterator& rhs):bface(rhs.bface), indices(rhs.indices), i(rhs.i) {}
  BF_Iterator(BF_Iterator&& rhs):bface(rhs.bface), indices(rhs.indices), i(rhs.i) {}
  BF_Iterator(BFacePtr _bface, IndicesPtr _indices, Index _i):bface(_bface), indices(_indices), i(_i) {}
  ~BF_Iterator() {}
public:
  bool operator!=(const BF_Iterator& rhs) { return i != rhs.i; }
  bool operator==(const BF_Iterator& rhs) { return i == rhs.i; }

  template<typename T = IterT&> std::enable_if_t<IsSequential, T>
  operator++() { i++; return *this; }
  template<typename T = IterT> std::enable_if_t<IsSequential, T>
  operator++(int) { IterT copy_one = *this; i++; return copy_one; }

  template<typename T = IterT&> std::enable_if_t<!IsSequential, T>
  operator++() { i--; return *this; }
  template<typename T = IterT> std::enable_if_t<!IsSequential, T>
  operator++(int) { IterT copy_one = *this; i--; return copy_one; }

  template<typename T = const Index&> std::enable_if_t<IsConstant, T>
  idx() { return bface->ctrl_pnt_idx[(*indices)[i]]; }
  template<typename T = const Vec3d&> std::enable_if_t<IsConstant, T>
  pnt() { return bface->ctrl_pnt[(*indices)[i]]; }
  template<typename T = std::tuple<const Index&,const  Vec3d&>> std::enable_if_t<IsConstant, T>
  operator*() { return std::tuple<const Index&, const Vec3d&>(idx(), pnt()); }

  template<typename T = Index&> std::enable_if_t<!IsConstant, T>
  idx() { return bface->ctrl_pnt_idx[(*indices)[i]]; }
  template<typename T = Vec3d&> std::enable_if_t<!IsConstant, T>
  pnt() { return bface->ctrl_pnt[(*indices)[i]]; }
  template<typename T = std::tuple<Index&, Vec3d&>> std::enable_if_t<!IsConstant, T>
  operator*() { return std::tuple<Index&, Vec3d&>(idx(), pnt()); }
};

template<bool IsConstant, bool IsSequential>
class BF_Range
{
public:
  using BFacePtr = std::conditional_t<IsConstant, const BezierFace*, BezierFace*>;
  using IterT = BF_Iterator<IsConstant, IsSequential>;
protected:
  BFacePtr bface;
  Indices local_indices;
public:
  BF_Range() = delete;
  BF_Range(BFacePtr _bface, Indices&& _local_indices):bface(_bface), local_indices(std::move(_local_indices)) {}
  // TODO: copy constructor.
  // TODO: move constructor.
  size_t size() { return local_indices.size(); }

  template<typename T = IterT> std::enable_if_t<IsSequential, T>
  begin() { return IterT(bface, &local_indices, 0); }
  template<typename T = IterT> std::enable_if_t<IsSequential, T>
  end() { return IterT(bface, &local_indices, local_indices.size()); }

  template<typename T = IterT> std::enable_if_t<!IsSequential, T>
  begin() { return IterT(bface, &local_indices, local_indices.size() - 1); }
  template<typename T = IterT> std::enable_if_t<!IsSequential, T>
  end() { return IterT(bface, &local_indices, -1); }
};

using BF_MutableSequencialRange = BF_Range<false, true>;
using BF_MutableReversedRange = BF_Range<false, false>;
using BF_ConstantSequencialRange = BF_Range<true, true>;
using BF_ConstantReversedRange = BF_Range<true, false>;

// We need to traversal two range at the same time.
// Thus, we need a "zip" class.
template<typename RT_Left, typename RT_Right>
class BF_Zip
{
  typedef typename RT_Left::IterT IterT_Left;
  typedef typename RT_Right::IterT IterT_Right;
private:
  RT_Left range_left;
  RT_Right range_right;
  class IterT
  {
  private:
    IterT_Left iter_left;
    IterT_Right iter_right;
  public:
    IterT(IterT_Left _left, IterT_Right _right):iter_left(_left), iter_right(_right) {}
    bool operator!=(const IterT& rhs) { return iter_left != rhs.iter_left; }
    bool operator==(const IterT& rhs) { return iter_left == rhs.iter_left; }
    IterT& operator++() { iter_left++; iter_right++; return *this; }
    IterT operator++(int) { IterT copy_one = *this; iter_left++; iter_right++; return copy_one; }

    auto operator*() -> std::tuple<decltype(*iter_left), decltype(*iter_right)>
    {
      return std::tuple<decltype(*iter_left), decltype(*iter_right)>(*iter_left, *iter_right);
    }
  };
public:
  BF_Zip() = delete;
  BF_Zip(RT_Left&& _range_left, RT_Right&& _range_right):range_left(std::move(_range_left)), range_right(std::move(_range_right))
  {
    ASSERT(range_left.size() == range_right.size(), "fail to zip two array with different sizes.")
  }
  IterT begin() { return IterT(range_left.begin(), range_right.begin()); }
  IterT end() { return IterT(range_left.end(), range_right.end()); }
};

class BezierFace
{
public:
  uint32_t degree;
  // Indices of three corner vertices. (pointing to TriConnectivity::vertices.)
  // Three corner vertices are 0-th, n(n+1)/2-th and n(n+1)/2+n-th
  // control point in the control points array.
  // The orientation is 0 -> n(n+1)/2 -> n(n+1)/2+n -> 0.
  Index3 vert_idx;
  // Indices of control points. (pointing to BMeshT::control_points.)
  Indices ctrl_pnt_idx;
  // Geometry implementation of Bezier triangle.
  // The control points are stored in the same order with ctrl_pnt_indx.
  Points ctrl_pnt;
public:
  // Constructions & copy & move
  BezierFace():degree(0), vert_idx(-1, -1, -1) {}
  BezierFace(const BezierFace& bface);
  BezierFace(BezierFace&& bface)noexcept;
  BezierFace& operator=(const BezierFace& bface);
  BezierFace& operator=(BezierFace&& bface)noexcept;
public:
  // Traversal functions
  Index find_local_vert_idx(Index v) const;
  Index find_local_edge_idx(Index v_from, Index v_to) const;

  Index local_corner_ctrlpnt_idx(Index v)const
  {
    return BezierTriangle::index_of_vert(degree, find_local_vert_idx(v));
  }
  Indices local_corner_ctrlpnt_idx()const
  {
    return BezierTriangle::indices_corner(degree);
  }
  Indices local_edge_ctrlpnt_idx(Index local_edge_idx)const
  {
    return BezierTriangle::indices_on_edge(degree, local_edge_idx);
  }
  Indices local_edge_ctrlpnt_idx(Index v_from, Index v_to)const
  {
    return BezierTriangle::indices_on_edge(degree, find_local_edge_idx(v_from, v_to));
  }
  Indices local_edge_interior_ctrlpnt_idx(Index local_edge_idx)const
  {
    return BezierTriangle::indices_inside_edge(degree, local_edge_idx);
  }
  Indices local_edge_interior_ctrlpnt_idx(Index v_from, Index v_to)const
  {
    return BezierTriangle::indices_inside_edge(degree, find_local_edge_idx(v_from, v_to));
  }
  Indices local_interior_ctrlpnt_idx()const
  {
    return BezierTriangle::indices_interior(degree);
  }
  Indices local_ctrlpnt_idx()const
  {
    return BezierTriangle::indices_all(degree);
  }
public:
  Vec3d& global_corner_ctrlpnt(Index v) { return ctrl_pnt[local_corner_ctrlpnt_idx(v)]; }
  const Vec3d& global_corner_ctrlpnt(Index v) const { return ctrl_pnt[local_corner_ctrlpnt_idx(v)]; }
  Index& global_corner_ctrlpnt_idx(Index v) { return ctrl_pnt_idx[local_corner_ctrlpnt_idx(v)]; }
  const Index& global_corner_ctrlpnt_idx(Index v) const { return ctrl_pnt_idx[local_corner_ctrlpnt_idx(v)]; }
public:
#define TRAVERSAL(GLOBAL_NAME, LOCAL_NAME)\
  BF_MutableSequencialRange GLOBAL_NAME##_ms()\
  { return BF_MutableSequencialRange(this, LOCAL_NAME()); }\
  BF_MutableReversedRange GLOBAL_NAME##_mr()\
  { return BF_MutableReversedRange(this, LOCAL_NAME()); }\
  BF_ConstantSequencialRange GLOBAL_NAME##_cs()const\
  { return BF_ConstantSequencialRange (this, LOCAL_NAME()); }\
  BF_ConstantReversedRange GLOBAL_NAME##_cr()const\
  { return BF_ConstantReversedRange (this, LOCAL_NAME()); }

  TRAVERSAL(global_corner_ctrlpnt, local_corner_ctrlpnt_idx);
  TRAVERSAL(global_interior_ctrlpnt, local_interior_ctrlpnt_idx);
  TRAVERSAL(global_ctrlpnt, local_ctrlpnt_idx);
#undef TRAVERSAL
#define TRAVERSAL(GLOBAL_NAME, LOCAL_NAME)\
  BF_MutableSequencialRange GLOBAL_NAME##_ms(Index v_from, Index v_to)\
  { return BF_MutableSequencialRange(this, LOCAL_NAME(v_from, v_to)); }\
  BF_MutableReversedRange GLOBAL_NAME##_mr(Index v_from, Index v_to)\
  { return BF_MutableReversedRange(this, LOCAL_NAME(v_from, v_to)); }\
  BF_ConstantSequencialRange GLOBAL_NAME##_cs(Index v_from, Index v_to)const\
  { return BF_ConstantSequencialRange (this, LOCAL_NAME(v_from, v_to)); }\
  BF_ConstantReversedRange GLOBAL_NAME##_cr(Index v_from, Index v_to)const\
  { return BF_ConstantReversedRange (this, LOCAL_NAME(v_from, v_to)); }

  TRAVERSAL(global_edge_ctrlpnt, local_edge_ctrlpnt_idx);
  TRAVERSAL(global_edge_interior_ctrlpnt, local_edge_interior_ctrlpnt_idx);
#undef TRAVERSAL
#define TRAVERSAL(GLOBAL_NAME, LOCAL_NAME)\
BF_MutableSequencialRange GLOBAL_NAME##_ms(Index local_edge_idx)\
  { return BF_MutableSequencialRange(this, LOCAL_NAME(local_edge_idx)); }\
  BF_MutableReversedRange GLOBAL_NAME##_mr(Index local_edge_idx)\
  { return BF_MutableReversedRange(this, LOCAL_NAME(local_edge_idx)); }\
  BF_ConstantSequencialRange GLOBAL_NAME##_cs(Index local_edge_idx)const\
  { return BF_ConstantSequencialRange (this, LOCAL_NAME(local_edge_idx)); }\
  BF_ConstantReversedRange GLOBAL_NAME##_cr(Index local_edge_idx)const\
  { return BF_ConstantReversedRange (this, LOCAL_NAME(local_edge_idx)); }

  TRAVERSAL(global_edge_ctrlpnt, local_edge_ctrlpnt_idx);
  TRAVERSAL(global_edge_interior_ctrlpnt, local_edge_interior_ctrlpnt_idx);
#undef TRAVERSAL
};

}// namespace ImageTriSimp
}// namespace GCLF