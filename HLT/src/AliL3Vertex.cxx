// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"
#include "AliL3RootTypes.h"
#include "AliL3Logging.h"
#include "AliL3VertexData.h"
#include "AliL3Vertex.h"

/** \class AliL3Vertex
<pre>
//_____________________________________________________________
// AliL3Vertex
//
// Stores the information of the vertex position
//
</pre>
*/

ClassImp(AliL3Vertex)

AliL3Vertex::AliL3Vertex(){
  //
  // default constructor for the AliL3Vertex class. 
  //

  SetZero();  
}

AliL3Vertex::~AliL3Vertex(){
  //
  // destructor
  //
}

void AliL3Vertex::SetZero()
{
  // set vertex to zero
  SetX(0);
  SetY(0);
  SetZ(0);
  SetXErr(1);
  SetYErr(1);
  SetZErr(1);
  fR=0;
  fPhi=0;
  fMWxy = 1.;
}

void AliL3Vertex::Read(const AliL3VertexData *vertex)
{
  // read vertex
  SetX(vertex->fX);
  SetY(vertex->fY);
  SetZ(vertex->fZ);
  SetXErr(vertex->fXErr);
  SetYErr(vertex->fYErr);
  SetZErr(vertex->fZErr);
  fR=0;
  fPhi=0;
  fMWxy = 1.;
}

