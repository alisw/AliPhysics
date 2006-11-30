// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTVertexData.h"
#include "AliHLTLogging.h"
#include "AliHLTVertexArray.h"
#include "AliHLTVertex.h"
#include "AliHLTVertexFinder.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTTransform.h"

/** \class AliHLTVertexFinder
<pre>
//_____________________________________________________________
// AliHLTVertexFinder
//
//   Implementation of AliHLTArray 
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

ClassImp(AliHLTVertexFinder)

AliHLTVertexFinder::AliHLTVertexFinder()
{
  //
  // default constructor for the AliHLTVertexFinder class. 
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

AliHLTVertexFinder::~AliHLTVertexFinder()
{
  //
  // destructor
  //
}

void AliHLTVertexFinder::Reset()
{
  //
  // Reset
  //
  ResetSector();
}


void AliHLTVertexFinder::Read(Int_t ncluster, AliHLTSpacePointData* hits )
{
  //
  //  analyze sector 
  // 
  
  const Int_t kseedRow = AliHLTTransform::GetNRows() - 1; 
  const Int_t kfirstRow = kseedRow-32;
  for(Int_t n=0;n<ncluster;n++){
    if(hits[n].fPadRow==kseedRow)
      FillSectorSeed3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy seeds in 3D
    if(hits[n].fPadRow<=kseedRow && hits[n].fPadRow>=kfirstRow)
      FillSector3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy data in 3D
  }
}

void AliHLTVertexFinder::Analyze()
{
  //
  // analyze all
  //
  FindSectorVertex();
  SetZ(GetZSector());
  SetZErr(GetZSectorErr());
  LOG(AliHLTLog::kInformational,"AliHLTVertexFinder::Analyze","Result")
  <<AliHLTLog::kDec<<"Vertex: "<<GetZ()<<"  RMS: "<<GetZErr()<<ENDLOG;
}

void AliHLTVertexFinder::Write(AliHLTVertex *vertex) const
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

void AliHLTVertexFinder::Write(AliHLTVertexData *vertex) const
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
