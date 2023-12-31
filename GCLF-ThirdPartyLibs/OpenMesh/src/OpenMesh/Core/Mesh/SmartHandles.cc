/* ========================================================================= *
 *                                                                           *
 *                               OpenMesh                                    *
 *           Copyright (c) 2001-2019, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openmesh.org                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenMesh.                                            *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
 * ========================================================================= */


//== INCLUDES =================================================================

#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>


//== NAMESPACES ===============================================================

namespace OpenMesh {

//TODO: I was not able to leave those in the header. If you find a way to do that, please do.
//GJP: I commented out the code below.
/*
PolyConnectivity::ConstVertexFaceRange SmartVertexHandle::faces() const
{
  assert(mesh() != nullptr);
  return mesh()->vf_range(*this);
}

PolyConnectivity::ConstVertexEdgeRange SmartVertexHandle::edges() const
{
  assert(mesh() != nullptr);
  return mesh()->ve_range(*this);
}

PolyConnectivity::ConstVertexVertexRange SmartVertexHandle::vertices() const
{
  assert(mesh() != nullptr);
  return mesh()->vv_range(*this);
}

PolyConnectivity::ConstVertexIHalfedgeRange SmartVertexHandle::incoming_halfedges() const
{
  assert(mesh() != nullptr);
  return mesh()->vih_range(*this);
}

PolyConnectivity::ConstVertexOHalfedgeRange SmartVertexHandle::outgoing_halfedges() const
{
  assert(mesh() != nullptr);
  return mesh()->voh_range(*this);
}


PolyConnectivity::ConstFaceVertexRange SmartFaceHandle::vertices() const
{
  assert(mesh() != nullptr);
  return mesh()->fv_range(*this);
}

PolyConnectivity::ConstFaceHalfedgeRange SmartFaceHandle::halfedges() const
{
  assert(mesh() != nullptr);
  return mesh()->fh_range(*this);
}

PolyConnectivity::ConstFaceEdgeRange SmartFaceHandle::edges() const
{
  assert(mesh() != nullptr);
  return mesh()->fe_range(*this);
}

PolyConnectivity::ConstFaceFaceRange SmartFaceHandle::faces() const
{
  assert(mesh() != nullptr);
  return mesh()->ff_range(*this);
}*/


//=============================================================================
} // namespace OpenMesh
//=============================================================================

//=============================================================================
