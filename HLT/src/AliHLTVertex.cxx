// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTVertexData.h"
#include "AliHLTVertex.h"

/** \class AliHLTVertex
<pre>
//_____________________________________________________________
// AliHLTVertex
//
// Stores the information of the vertex position
//
</pre>
*/

ClassImp(AliHLTVertex)

AliHLTVertex::AliHLTVertex(){
  //
  // default constructor for the AliHLTVertex class. 
  //

  SetZero();  
}

AliHLTVertex::~AliHLTVertex(){
  //
  // destructor
  //
}

void AliHLTVertex::SetZero()
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

void AliHLTVertex::Read(const AliHLTVertexData *vertex)
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

