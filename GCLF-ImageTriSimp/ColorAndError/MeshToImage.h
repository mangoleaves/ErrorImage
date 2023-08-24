#pragma once
#include "Image/ImageDefinition.h"
#include "Mesh/LinearMeshDefinition.h"
#include "Mesh/BezierMeshDefinition.h"

namespace GCLF
{
namespace ImageTriSimp
{

template<typename Simplifier>
class MeshToImage
{
public:
  MeshToImage() = delete;
  MeshToImage(Simplifier* _simplifier) :simplifier(_simplifier) {}

  ImageT operator()();
private:
  Simplifier* simplifier;
};

}// namespace ImageTriSimp
}// namespace GCLF

#include "MeshToImage_impl.h"