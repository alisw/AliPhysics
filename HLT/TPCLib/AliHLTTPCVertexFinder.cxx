// @(#) $Id$
// Original: AliHLTVertexFinder.cxx,v 1.9 2005/06/14 10:55:21 cvetan 

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCVertexData.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCVertexArray.h"
#include "AliHLTTPCVertex.h"
#include "AliHLTTPCVertexFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTransform.h"

/** \class AliHLTTPCVertexFinder
<pre>
//_____________________________________________________________
// AliHLTTPCVertexFinder
//
//   Implementation of AliHLTTPCArray 
//   usage:
// 
//   ResetSector();
//   for(n=0;n<NMEMSEC;n++)  
//     Read();
//   FindSectorVertex();
//   SetZ(GetZSector());
//   SetZErr(GetZErrSector());
// 
//
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCVertexFinder)

AliHLTTPCVertexFinder::AliHLTTPCVertexFinder()
{
  //
  // default constructor for the AliHLTTPCVertexFinder class. 
  //

  //Set vertex to zero.
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

AliHLTTPCVertexFinder::~AliHLTTPCVertexFinder()
{
  //
  // destructor
  //
}

void AliHLTTPCVertexFinder::Reset()
{
  //
  // Reset
  //
  ResetSector();
}


void AliHLTTPCVertexFinder::Read(Int_t ncluster, AliHLTTPCSpacePointData* hits )
{
  //
  //  analyze sector 
  // 
  
  const Int_t kseedRow = AliHLTTPCTransform::GetNRows() - 1; 
  const Int_t kfirstRow = kseedRow-32;
  for(Int_t n=0;n<ncluster;n++){
    if(hits[n].fPadRow==kseedRow)
      FillSectorSeed3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy seeds in 3D
    if(hits[n].fPadRow<=kseedRow && hits[n].fPadRow>=kfirstRow)
      FillSector3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy data in 3D
  }
}

void AliHLTTPCVertexFinder::Analyze()
{
  //
  // analyze all
  //
  FindSectorVertex();
  SetZ(GetZSector());
  SetZErr(GetZSectorErr());
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCVertexFinder::Analyze","Result")
  <<AliHLTTPCLog::kDec<<"Vertex: "<<GetZ()<<"  RMS: "<<GetZErr()<<ENDLOG;
}

void AliHLTTPCVertexFinder::Write(AliHLTTPCVertex *vertex) const
{
  //
  // write
  //
  vertex->SetX(GetX());
  vertex->SetY(GetZ());
  vertex->SetZ(GetZ());
  vertex->SetXErr(GetXErr());
  vertex->SetYErr(GetYErr());
  vertex->SetZErr(GetZErr());

  vertex->SetXYWeight(GetXYWeight());
}

void AliHLTTPCVertexFinder::Write(AliHLTTPCVertexData *vertex) const
{
  //
  // write
  //
  vertex->fX=GetX();
  vertex->fY=GetZ();
  vertex->fZ=GetZ();
  vertex->fXErr=GetXErr();
  vertex->fYErr=GetYErr();
  vertex->fZErr=GetZErr();
}
