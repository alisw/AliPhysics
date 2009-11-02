/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2008

//=============================================================================
// AliUnicorEvent-AliESD interface                                                   //
//=============================================================================

#include <cmath>
#include "AliESDEvent.h"
#include "AliUnicorEventAliceESD.h"

ClassImp(AliUnicorEventAliceESD)

//=============================================================================
AliUnicorEventAliceESD::AliUnicorEventAliceESD(AliESDEvent *esd) : AliUnicorEvent(), fESD(esd) 
{
  // constructor 

  //  printf("%s object created\n",ClassName());
  if (!fESD) fESD = new AliESDEvent();
}
//=============================================================================
AliUnicorEventAliceESD::~AliUnicorEventAliceESD()
{
  // destructor

}
//=============================================================================
Bool_t AliUnicorEventAliceESD::Good() const 
{
  // event cuts

  if (fabs(Zver())>1) return kFALSE;
  if (fESD->GetPrimaryVertex()->GetZRes()>0.1) return kFALSE;
  return kTRUE;
}
//=============================================================================
Bool_t AliUnicorEventAliceESD::ParticleGood(Int_t i, Int_t pidi) const 
{
  // track cuts and particle id cuts; pidi=0 means take all species
  // consider using the standard ESDcut

  // track quality cuts

  AliESDtrack *track = fESD->GetTrack(i);
  if (!track->IsOn(AliESDtrack::kTPCrefit)) return 0;        // TPC refit
  if (track->GetTPCNcls() < 100) return 0;                   // number of TPC clusters
  const AliExternalTrackParam *tp = track->GetTPCInnerParam();
  if (!tp) return 0;           

  Float_t r,z;
  track->GetImpactParameters(r,z);
  if (fabs(z)>0.2) return 0;                          // impact parameter in z
  if (fabs(r)>0.1) return 0;                          // impact parameter in xy

  // pid

  if (pidi==0) return 1;
  if (!track->IsOn(AliESDtrack::kTPCpid)) return 0;
  Double_t p[AliPID::kSPECIES];
  track->GetTPCpid(p);
  Int_t q = tp->Charge();
  if (pidi == -211) return p[AliPID::kPion]+p[AliPID::kMuon]>0.5 && q==-1;
  else if (pidi == 211) return p[AliPID::kPion]+p[AliPID::kMuon]>0.5 && q==1;
  else return 0;
}
//=============================================================================
Bool_t AliUnicorEventAliceESD::PairGood(Double_t /*p0*/, Double_t the0, Double_t phi0, 
				Double_t /*p1*/, Double_t the1, Double_t phi1) const {

  // two-track separation cut

  double dthe = the1-the0;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi0);
  return (fabs(dthe)>0.005 || fabs(dphi)>0.030);
}
//=============================================================================
