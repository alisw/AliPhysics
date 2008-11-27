// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
// AliDEvent-AliESD interface                                                   //
//=============================================================================

#include <cmath>
#include "AliESDEvent.h"
#include "AliDEventAliceESD.h"

ClassImp(AliDEventAliceESD)

//=============================================================================
AliDEventAliceESD::AliDEventAliceESD() : AliDEvent(), AliESDEvent() 
{
  // constructor 

  printf("%s object created\n",AliDEvent::ClassName());
}
//=============================================================================
AliDEventAliceESD::~AliDEventAliceESD()
{
  // destructor

}
//=============================================================================
Bool_t AliDEventAliceESD::Good() const 
{
  // event cuts

  if (fabs(Zver())>1) return kFALSE;
  if (AliESDEvent::GetPrimaryVertex()->GetZRes()>0.1) return kFALSE;
  return kTRUE;
}
//=============================================================================
Bool_t AliDEventAliceESD::ParticleGood(Int_t i, Int_t pidi) const 
{
  // track cuts and particle id cuts; pidi=0 means take all tracks
  // consider using the standard ESDcut

  // track quality cuts

  AliESDtrack *track = AliESDEvent::GetTrack(i);
  if (!track->IsOn(AliESDtrack::kTPCrefit)) return 0;        // TPC refit
  if (track->GetTPCNcls() < 120) return 0;                   // number of TPC clusters
  const AliExternalTrackParam *tp = track->GetTPCInnerParam();
  if (!tp) return 0;           

  Float_t r,z;
  track->GetImpactParameters(r,z);
  //  if (fabs(z)>0.2) return 0;                          // impact parameter in z
  //  if (fabs(r)>0.1) return 0;                          // impact parameter in xy

  // pid

  if (pidi==0) return 1;
  if (!track->IsOn(AliESDtrack::kTPCpid)) return 0;
  Double_t p[AliPID::kSPECIES];
  track->GetESDpid(p);
  Int_t q = tp->Charge();
  if (pidi == -211) return p[AliPID::kPion]>0.5 && q==-1;
  else if (pidi == 211) return p[AliPID::kPion]>0.5 && q==1;
  else return 0;
}
//=============================================================================
Bool_t AliDEventAliceESD::PairGood(Double_t /*p0*/, Double_t the0, Double_t phi0, 
				   Double_t /*p1*/, Double_t the1, Double_t phi1) const {

  // two-track separation cut

  double r = 85; // TPC entrance radius in cm
  double x0 = r*sin(the0)*cos(phi0);
  double x1 = r*sin(the1)*cos(phi0);
  double y0 = r*sin(the0)*sin(phi0);
  double y1 = r*sin(the1)*sin(phi1);
  double z0 = r*cos(the0);
  double z1 = r*cos(the1);
  double dx = x1-x0;
  double dy = y1-y0;
  double dz = z1-z0;
  double dist2 = dx*dx+dy*dy+dz*dz;
  return dist2>2*2;
}
//=============================================================================
