#pragma once

#include "AABBNode.h"

namespace GCLF
{
namespace Geometry
{

enum class TraversalSequence
{
  I_DONT_KNOW,
  NONE,
  ONLY_LEFT,
  ONLY_RIGHT,
  LEFT_THEN_RIGHT,
  RIGHT_THEN_LEFT
};

template<typename Kernel>
class AABBTree
{
private:
  typedef Kernel K;
  // Primitive from Kernel
  typedef typename K::Primitive Primitive;
  typedef std::vector<Primitive> Primitives;
  typedef typename Primitives::iterator PrimitiveIter;
  typedef Primitive* PrimitivePtr;
  // CalcBox from Kernel
  typedef typename K::CalcBox CalcBox;
  // SplitPref from Kernel
  typedef typename K::SplitPred SplitPred;

  typedef AABBNode<Primitive> Node;
protected:
  Primitives m_primitives;
  Node* m_p_root_node = nullptr;
public:
  AABBTree() {}

  AABBTree(Primitives&& primitives)
  {
    insert(std::move(primitives));
    build();
  }
  AABBTree(PrimitiveIter first, PrimitiveIter beyond)
  {
    insert(first, beyond);
    build();
  }

  ~AABBTree()
  {
    clear();
  }

  void insert(Primitives&& primitives);
  void insert(PrimitiveIter first, PrimitiveIter beyond);
  void build();

  void clear();
  void clear_nodes();

  inline size_t size() const { return m_primitives.size(); }
  inline bool empty() const { return m_primitives.empty(); }

  template<typename Trait>
  void traversal(Trait& traits);

protected:
  BoundingBox calc_bbox(PrimitiveIter first, PrimitiveIter beyond);

  // build functions
  void split_primitives(PrimitiveIter first, PrimitiveIter beyond, const BoundingBox& box);
  void expand(Node* node, PrimitiveIter first, PrimitiveIter beyond, const size_t range);

  template<typename Trait>
  bool traversal_node(Node* node, Trait& traits, const size_t nb_primitives);
};
}// namespace Geometry
}// namespace GCLF

#include "AABBTree_impl.h"