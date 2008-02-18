//-------------------------------------------------------------------------
//     PWG2 specific additional information for the AOD Track
//     Stores a per-track information for the AOD track that is not 
//     included in the standard AliAODTrack
//     Author: Adam Kisiel, OSU, Adam.Kisiel@cern.ch
//-------------------------------------------------------------------------

#include <TBits.h>

#include <AliAODTrack.h>
#include "AliPWG2AODTrack.h"

ClassImp(AliPWG2AODTrack)

AliPWG2AODTrack::AliPWG2AODTrack():
  fSharedMap(160),
  fClusterMap(160),
  fAODTrack(NULL)
{
  // Default constructor
  SetTPCNominalEntrancePoint();
  SetTPCNominalExitPoint();
  fSharedMap.ResetAllBits(kFALSE);
  fClusterMap.ResetAllBits(kTRUE);
}

AliPWG2AODTrack::AliPWG2AODTrack(Double_t tpcentr[3],
				 Double_t tpcexit[3],
				 TBits tpcshare,
				 TBits tpcclus,
				 AliAODTrack *track):
  fSharedMap(tpcshare),
  fClusterMap(tpcclus),
  fAODTrack(track)
{
  // Constructor initializing all fields
  SetTPCNominalEntrancePoint(tpcentr);
  SetTPCNominalExitPoint(tpcexit);
}

AliPWG2AODTrack::~AliPWG2AODTrack()
{
}

AliPWG2AODTrack::AliPWG2AODTrack(const AliPWG2AODTrack& trk):
  TObject(),
  fSharedMap(trk.fSharedMap),
  fClusterMap(trk.fClusterMap),
  fAODTrack(trk.fAODTrack)
{
  // Copy constructor
  Double_t tpcp[3];
  trk.GetTPCNominalEntrancePoint(tpcp);
  SetTPCNominalEntrancePoint(tpcp);
  trk.GetTPCNominalExitPoint(tpcp);
  SetTPCNominalExitPoint(tpcp);
}

AliPWG2AODTrack& AliPWG2AODTrack::operator=(const AliPWG2AODTrack& trk)
{
  // Assignment operator
  if(this!=&trk) {
    fSharedMap = trk.fSharedMap;
    fClusterMap = trk.fClusterMap;
    fAODTrack = trk.fAODTrack;

    Double_t tpcp[3];
    trk.GetTPCNominalEntrancePoint(tpcp);
    SetTPCNominalEntrancePoint(tpcp);
    trk.GetTPCNominalExitPoint(tpcp);
    SetTPCNominalExitPoint(tpcp);
  }
  return *this;
}

void AliPWG2AODTrack::GetTPCNominalEntrancePoint(Double_t *tpce) const
{
  // Return TPC entrance point coordinates
  tpce[0] = fTPCNominalEntrancePoint[0];
  tpce[1] = fTPCNominalEntrancePoint[1];
  tpce[2] = fTPCNominalEntrancePoint[2]; 
}

void AliPWG2AODTrack::GetTPCNominalExitPoint(Double_t *tpce) const
{
  // Return TPC exit point coordinates
  tpce[0] = fTPCNominalExitPoint[0];
  tpce[1] = fTPCNominalExitPoint[1];
  tpce[2] = fTPCNominalExitPoint[2]; 
}

void AliPWG2AODTrack::SetTPCNominalEntrancePoint(Double_t *tpce)
{
  // Set TPC entrance point coordinates
  if (tpce) {
    fTPCNominalEntrancePoint[0] = tpce[0];
    fTPCNominalEntrancePoint[1] = tpce[1];
    fTPCNominalEntrancePoint[2] = tpce[2];
  }
  else {
    fTPCNominalEntrancePoint[0] = 0.0;
    fTPCNominalEntrancePoint[1] = 0.0;
    fTPCNominalEntrancePoint[2] = 0.0;
  }
}

void AliPWG2AODTrack::SetTPCNominalExitPoint(Double_t *tpce)
{
  // Set TPC exit point coordinates
  if (tpce) {
    fTPCNominalExitPoint[0] = tpce[0];
    fTPCNominalExitPoint[1] = tpce[1];
    fTPCNominalExitPoint[2] = tpce[2];
  }
  else {
    fTPCNominalExitPoint[0] = 0.0;
    fTPCNominalExitPoint[1] = 0.0;
    fTPCNominalExitPoint[2] = 0.0;
  }
}

const TBits &AliPWG2AODTrack::GetTPCSharedMap() const
{
  return fSharedMap;
}
const TBits &AliPWG2AODTrack::GetTPCClusterMap() const
{
  return fClusterMap;
}

void AliPWG2AODTrack::SetTPCSharedMap(const TBits &bits)
{
  fSharedMap = bits;
}

void AliPWG2AODTrack::SetTPCClusterMap(const TBits &bits)
{
  fClusterMap = bits;
}

void AliPWG2AODTrack::SetAODTrackRef(AliAODTrack *track)
{
  fAODTrack = track;
}

AliAODTrack *AliPWG2AODTrack::GetRefAODTrack()
{
  return (AliAODTrack *) fAODTrack.GetObject();
}

