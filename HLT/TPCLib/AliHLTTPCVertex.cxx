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

AliHLTTPCVertex::AliHLTTPCVertex()
  :
  fX(0.0),
  fY(0.0),  
  fZ(0.0),  
  fPhi(0.0),
  fR(0.0),  
  fXErr(1.0),
  fYErr(1.0),
  fZErr(1.0),
  fMWxy(1.0)
{
  //
  // default constructor for the AliHLTTPCVertex class. 
  //

  SetZero();  
}

AliHLTTPCVertex::AliHLTTPCVertex(const AliHLTTPCVertex&)
  :
  fX(0.0),
  fY(0.0),  
  fZ(0.0),  
  fPhi(0.0),
  fR(0.0),  
  fXErr(1.0),
  fYErr(1.0),
  fZErr(1.0),
  fMWxy(1.0)
{
  //
  // copy constructor not for use
  //
}

AliHLTTPCVertex& AliHLTTPCVertex::operator=(const AliHLTTPCVertex&)
{
  //
  // assignment operator not for use
  //
  fX=0.0;
  fY=0.0;  
  fZ=0.0;  
  fPhi=0.0;
  fR=0.0;  
  fXErr=1.0;
  fYErr=1.0;
  fZErr=1.0;
  fMWxy=1.0;
  return *this;
}

AliHLTTPCVertex::~AliHLTTPCVertex()
{
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
