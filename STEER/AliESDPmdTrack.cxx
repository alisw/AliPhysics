#include "AliESDPmdTrack.h"

ClassImp(AliESDPmdTrack)

AliESDPmdTrack::AliESDPmdTrack (const AliESDPmdTrack& PMDTrack):TObject(PMDTrack)
{
  // Constructor
  fDet    = PMDTrack.fDet;
  fTheta  = PMDTrack.fTheta;
  fPhi    = PMDTrack.fPhi;
  fCluADC = PMDTrack.fCluADC;
  fCluPID = PMDTrack.fCluPID;
}
//--------------------------------------------------------------------------//
AliESDPmdTrack &AliESDPmdTrack::operator=(const AliESDPmdTrack& PMDTrack)
{
  // Copy constructor
  if(&PMDTrack == this) return *this;
  fDet    = PMDTrack.fDet;
  fTheta  = PMDTrack.fTheta;
  fPhi    = PMDTrack.fPhi;
  fCluADC = PMDTrack.fCluADC;
  fCluPID = PMDTrack.fCluPID;
  return *this;
}
