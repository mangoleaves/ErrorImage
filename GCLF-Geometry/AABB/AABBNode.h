#pragma once

#include "Basic/BoundingBox.h"

namespace GCLF
{
namespace Geometry
{

template<typename Primitive>
class AABBNode
{
public:
  typedef Primitive* PrimitivePtr;
  typedef AABBNode<Primitive> Node;
public:
  BoundingBox m_bbox;
  void* m_p_left_child, * m_p_right_child;
public:
  AABBNode()
    :m_bbox(),
    m_p_left_child(nullptr),
    m_p_right_child(nullptr)
  {}

  ~AABBNode() {}

  inline BoundingBox& bbox() { return m_bbox; }
  inline void*& left_ptr() {return m_p_left_child;}
  inline void*& right_ptr() {return m_p_right_child;}
  inline Node* left_child() { return static_cast<Node*>(m_p_left_child); }
  inline Node* right_child() { return static_cast<Node*>(m_p_right_child); }
  inline Primitive* left_data() { return static_cast<PrimitivePtr>(m_p_left_child); }
  inline Primitive* right_data() { return static_cast<PrimitivePtr>(m_p_right_child); }
};
}// namespace Geometry
}// namespace GCLF