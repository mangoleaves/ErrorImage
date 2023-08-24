#pragma once
#include "CalcErrorOnFace.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
class CalcErrorOnMesh
{
public:
  CalcErrorOnMesh() = delete;
  CalcErrorOnMesh(Simplifier* _simplifier) :simplifier(_simplifier) {}

  double p = 2.0;

  double lp_error();
  double mean_lp_error();
  double root_mean_lp_error();
private:
  Simplifier* simplifier;
};

}// namespace ImageTriSimp
}// namespace GCLF

#include "CalcErrorOnMesh_impl.h"