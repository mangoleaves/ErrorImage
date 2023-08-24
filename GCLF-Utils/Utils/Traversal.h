#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include "Logger/Logger.h"

namespace GCLF
{
namespace Utils
{

template<
  typename DST_STD_CONTAINER,
  typename IDX_STD_CONTAINER,
  typename SRC_STD_CONTAINER
>
void copy_certain(
  DST_STD_CONTAINER& dst,
  const IDX_STD_CONTAINER& idx,
  const SRC_STD_CONTAINER& src)
{
  ASSERT(idx.size() == src.size(), "Mismatch size.");
  for (size_t i = 0;i < idx.size();i++)
  {
    dst[idx[i]] = src[i];
  }
}

}// namespace Utils
}// namespace GCLF