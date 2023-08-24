#include "Marker.h"

namespace GCLF
{
namespace Utils
{

void Marker::unmark_all()
{
  if (timestamp < 255)
  {
    timestamp++;
  }
  else
  {
    timestamp = 1;
    size_t primitive_size = m_mark.size();
    m_mark.clear();
    m_mark.resize(primitive_size, 0);
  }
}

}// namespace Utils
}// namespace GCLF