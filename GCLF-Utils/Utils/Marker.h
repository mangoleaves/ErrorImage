#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>

namespace GCLF
{
namespace Utils
{

class Marker
{
public:
  Marker() = default;
  Marker(size_t primitive_size) { m_mark.resize(primitive_size, 0); }
  Marker(const Marker& rhs) { m_mark = rhs.m_mark; }

  void mark(int idx) { m_mark[idx] = timestamp; }

  inline bool is_marked(int idx) { return m_mark[idx] == timestamp; }

  void unmark_all();
private:
  uint8_t timestamp = 1;
  std::vector<uint8_t> m_mark;
};

}// namespace Utils
}// namespace GCLF