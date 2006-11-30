// @(#) $Id$
// Original: AliHLTVertex.cxx,v 1.5 2004/07/02 11:41:18 loizides Exp $

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCVertex.h"

/** \class AliHLTTPCVertex
<pre>
//_____________________________________________________________
// AliHLTTPCVertex
//
// Stores the information of the vertex position
//
</pre>
*/

ClassImp(AliHLTTPCVertex)

AliHLTTPCVertex::AliHLTTPCVertex(){
  //
  // default constructor for the AliHLTTPCVertex class. 
  //

  SetZero();  
}

AliHLTTPCVertex::~AliHLTTPCVertex(){
  //
  // destructor
  //
}

void AliHLTTPCVertex::SetZero()
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

void AliHLTTPCVertex::Read(const AliHLTTPCVertexData *vertex)
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

