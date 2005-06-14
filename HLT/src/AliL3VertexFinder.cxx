// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"
#include "AliL3RootTypes.h"
#include "AliL3SpacePointData.h"
#include "AliL3VertexData.h"
#include "AliL3Logging.h"
#include "AliL3VertexArray.h"
#include "AliL3Vertex.h"
#include "AliL3VertexFinder.h"
#include "AliL3SpacePointData.h"
#include "AliL3Transform.h"

/** \class AliL3VertexFinder
<pre>
//_____________________________________________________________
// AliL3VertexFinder
//
//   Implementation of AliL3Array 
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

ClassImp(AliL3VertexFinder)

AliL3VertexFinder::AliL3VertexFinder()
{
  //
  // default constructor for the AliL3VertexFinder class. 
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

AliL3VertexFinder::~AliL3VertexFinder()
{
  //
  // destructor
  //
}

void AliL3VertexFinder::Reset()
{
  //
  // Reset
  //
  ResetSector();
}


void AliL3VertexFinder::Read(Int_t ncluster, AliL3SpacePointData* hits )
{
  //
  //  analyze sector 
  // 
  
  const Int_t kseedRow = AliL3Transform::GetNRows() - 1; 
  const Int_t kfirstRow = kseedRow-32;
  for(Int_t n=0;n<ncluster;n++){
    if(hits[n].fPadRow==kseedRow)
      FillSectorSeed3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy seeds in 3D
    if(hits[n].fPadRow<=kseedRow && hits[n].fPadRow>=kfirstRow)
      FillSector3D(hits[n].fX,hits[n].fY,hits[n].fZ);  //copy data in 3D
  }
}

void AliL3VertexFinder::Analyze()
{
  //
  // analyze all
  //
  FindSectorVertex();
  SetZ(GetZSector());
  SetZErr(GetZSectorErr());
  LOG(AliL3Log::kInformational,"AliL3VertexFinder::Analyze","Result")
  <<AliL3Log::kDec<<"Vertex: "<<GetZ()<<"  RMS: "<<GetZErr()<<ENDLOG;
}

void AliL3VertexFinder::Write(AliL3Vertex *vertex) const
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

void AliL3VertexFinder::Write(AliL3VertexData *vertex) const
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
