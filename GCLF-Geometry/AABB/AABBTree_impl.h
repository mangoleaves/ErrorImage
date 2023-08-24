#pragma once

#include "AABBTree.h"

namespace GCLF
{
namespace Geometry
{
template<typename Kernel>
void AABBTree<Kernel>::insert(Primitives&& primitives)
{
  clear();
  m_primitives = std::move(primitives);
}

template<typename Kernel>
void AABBTree<Kernel>::insert(PrimitiveIter first, PrimitiveIter beyond)
{
  clear();
  m_primitives.reserve(std::distance(first, beyond));
  while (first != beyond)
  {
    m_primitives.push_back(*first);
    first++;
  }
}

template<typename Kernel>
void AABBTree<Kernel>::build()
{
  clear_nodes();
  if (m_primitives.size() > 1)
  {
    // allocate tree nodes.
    m_p_root_node = new Node[m_primitives.size() - 1]();
    if (m_p_root_node == nullptr)
    {
      clear();
      throw std::bad_alloc();
    }
    // construct AABB tree.
    expand(m_p_root_node, m_primitives.begin(), m_primitives.end(), m_primitives.size());
  }
}

template<typename Kernel>
void AABBTree<Kernel>::clear()
{
  clear_nodes();
  m_primitives.clear();
}

template<typename Kernel>
void AABBTree<Kernel>::clear_nodes()
{
  if (m_p_root_node)
  {
    delete[] m_p_root_node;
  }
  m_p_root_node = nullptr;
}

template<typename Kernel>
template<typename Trait>
void AABBTree<Kernel>::traversal(Trait& traits)
{
  switch (size())
  {
  case 0:
    break;
  case 1:
    traits.intersection(m_primitives[0]);
    break;
  default:
    traversal_node(m_p_root_node, traits, m_primitives.size());
  }
}

template<typename Kernel>
BoundingBox AABBTree<Kernel>::calc_bbox(PrimitiveIter first, PrimitiveIter beyond)
{
  CalcBox cb;
  BoundingBox bbox = cb(*first);
  for (++first; first != beyond; ++first)
  {
    bbox += cb(*first);
  }
  return bbox;
}


template<typename Kernel>
void AABBTree<Kernel>::split_primitives(PrimitiveIter first, PrimitiveIter beyond, const BoundingBox& box)
{
  PrimitiveIter middle = first + (beyond - first) / 2;
  size_t split_dim = box.longest_axis();
  SplitPred pred(split_dim);
  std::nth_element(first, middle, beyond, pred);
}

template<typename Kernel>
void AABBTree<Kernel>::expand(
  Node* node, PrimitiveIter first, PrimitiveIter beyond, const size_t range)
{
  node->bbox() = calc_bbox(first, beyond);
  split_primitives(first, beyond, node->bbox());

  switch (range)
  {
  case 2:
    node->left_ptr() = &(*first);
    node->right_ptr() = &(*(++first));
    break;
  case 3:
    node->left_ptr() = &(*first);
    node->right_ptr() = (void*)(node + 1);
    expand(node->right_child(), first + 1, beyond, 2);
    break;
  default:
    const size_t new_range = range / 2;

    node->left_ptr() = (void*)(node + 1);
    node->right_ptr() = (void*)(node + new_range);
    expand(node->left_child(), first, first + new_range, new_range);
    expand(node->right_child(), first + new_range, beyond, range - new_range);
    break;
  }
}

template<typename Kernel>
template<typename Trait>
bool AABBTree<Kernel>::traversal_node(
  Node* node, Trait& trait, const size_t nb_primitives)
{
  bool go_next = true;
  switch (nb_primitives)
  {
  case 2:
    go_next = trait.intersection(*node->left_data());
    if (go_next)
      go_next = trait.intersection(*node->right_data());
    return go_next;
  case 3:
    go_next = trait.intersection(*node->left_data());
    if (go_next && trait.do_inter(node->right_child()->bbox()))
      go_next = traversal_node(node->right_child(), trait, 2);
    return go_next;
  default:
    TraversalSequence ts = trait.which_traversal_first(node->left_child()->bbox(), node->right_child()->bbox());
    switch (ts)
    {
    case TraversalSequence::I_DONT_KNOW:
    {
      if (trait.do_inter(node->left_child()->bbox()))
      {
        go_next = traversal_node(node->left_child(), trait, nb_primitives / 2);
        if (go_next && trait.do_inter(node->right_child()->bbox()))
          go_next = traversal_node(node->right_child(), trait, nb_primitives - nb_primitives / 2);
        return go_next;
      }
      else if (trait.do_inter(node->right_child()->bbox()))
      {
        return traversal_node(node->right_child(), trait, nb_primitives - nb_primitives / 2);
      }
    }break;
    case TraversalSequence::NONE:break;
    case TraversalSequence::LEFT_THEN_RIGHT:
    {
      go_next = traversal_node(node->left_child(), trait, nb_primitives / 2);
      if (go_next)
        go_next = traversal_node(node->right_child(), trait, nb_primitives - nb_primitives / 2);
      return go_next;
    }break;
    case TraversalSequence::RIGHT_THEN_LEFT:
    {
      go_next = traversal_node(node->right_child(), trait, nb_primitives - nb_primitives / 2);
      if (go_next)
        go_next = traversal_node(node->left_child(), trait, nb_primitives / 2);
      return go_next;
    }break;
    case TraversalSequence::ONLY_LEFT:
    {
      go_next = traversal_node(node->left_child(), trait, nb_primitives / 2);
      return go_next;
    }break;
    case TraversalSequence::ONLY_RIGHT:
    {
      go_next = traversal_node(node->right_child(), trait, nb_primitives - nb_primitives / 2);
      return go_next;
    }break;
    }
    return true;
  }
}
}// namespace Geometry
}// namespace GCLF