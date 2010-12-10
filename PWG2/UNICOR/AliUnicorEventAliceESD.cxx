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
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "AliESDVZERO.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliUnicorEventAliceESD.h"

ClassImp(AliUnicorEventAliceESD)

//=============================================================================
AliUnicorEventAliceESD::AliUnicorEventAliceESD(AliESDEvent *esd) : AliUnicorEvent(), fViper(0), fESD(esd)//, fPhysicsSelection(0) 
{
  // constructor 

  //  printf("%s object created\n",ClassName());
  if (!fESD) fESD = new AliESDEvent();
  //  fPhysicsSelection = new AliPhysicsSelection();

  // V0 percentile graph for centrality determination

  //  TFile::Open("$ALICE_ROOT/ANALYSIS/macros/AliCentralityBy1D_137161_v4.root","read");
  TFile::Open("$ALICE_ROOT/ANALYSIS/macros/AliCentralityBy1D_137366_v2.root","read");
  const TH1F *hi = (const TH1F*) gFile->Get("hmultV0_percentile");
  //const TH1D *hi = (const TH1D*) gFile->Get("hNtracks_percentile");
  //const TH1F *hi = (const TH1F*) gFile->Get("hNclusters1_percentile");
  fViper = new TGraph((const TH1*) hi);
  gFile->Close();
}
//=============================================================================
AliUnicorEventAliceESD::~AliUnicorEventAliceESD() {

  // destructor

}
//=============================================================================
Bool_t AliUnicorEventAliceESD::Good() const 
{
  // event cuts

  //  if (!fPhysicsSelection->IsCollisionCandidate(fESD)) return kFALSE;
  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  if (!vtx->GetStatus()) return kFALSE;
  if (fabs(Zver())>1) return kFALSE;
  if (NGoodParticles()<9) return kFALSE;
  return kTRUE;
}
//=============================================================================
Double_t AliUnicorEventAliceESD::Centrality() const {

  // centrality between 0 (central) and 1 (very peripheral)

  // V0 multiplicity, not linearized
  
  AliESDVZERO *v0 = fESD->GetVZEROData();
  float multa = v0->GetMTotV0A();
  float multc = v0->GetMTotV0C();
  double cent = fViper->Eval(multa+multc)/100.0;

  // number of tracks

  // double cent = fViper->Eval(((AliVEvent*) fESD)->GetNumberOfTracks())/100.0;  

  // number of clusters in second layer of SPD

  //  double nspd2 = fESD->GetMultiplicity()->GetNumberOfITSClusters(1);
  //  double zv = fESD->GetPrimaryVertexSPD()->GetZ();
  //  const double pars[] = {8.10030e-01,-2.80364e-03,-7.19504e-04};
  //  zv -= pars[0];
  //  float corr = 1 + zv*(pars[1] + zv*pars[2]);
  //  nspd2 = corr>0? nspd2/corr : -1;
  //  double cent = fViper->Eval(nspd2)/100.0;  

  return cent;
}
//=============================================================================
Bool_t AliUnicorEventAliceESD::ParticleGood(Int_t i, Int_t pidi) const 
{
  // track cuts and particle id cuts; pidi=0 means take all species
  // consider using the standard ESDcut

  // track quality cuts

  AliESDtrack *track = fESD->GetTrack(i);
  if (!track->IsOn(AliESDtrack::kTPCrefit)) return 0;        // TPC refit
  if (!track->IsOn(AliESDtrack::kITSrefit)) return 0;        // ITS refit
  if (track->GetTPCNcls() < 90) return 0;                    // number of TPC clusters
  if (track->GetKinkIndex(0) > 0) return 0;                  // no kink daughters
  const AliExternalTrackParam *tp = GetTrackParam(i);
  if (!tp) return 0;                                         // track param
  if (fabs(tp->Eta())>0.8) return 0;                         // fiducial pseudorapidity

  //  double pi9 = TMath::Pi()/9.0;
  //  double eta = tp->Eta();
  //  double phi = ParticlePhi(i);
  //  if (eta>0 && phi>-8*pi9 && phi<-7*pi9) return 0; // A10
  //  if (eta>0 && phi> 1*pi9 && phi< 2*pi9) return 0; // A01
  //  if (tp->Pt()<0.2) return 0;                      // lower pt cutoff

  Float_t r,z;
  track->GetImpactParametersTPC(r,z);
  if (fabs(z)>1.0) return 0;                          // impact parameter in z
  if (fabs(r)>1.0) return 0;                          // impact parameter in xy
  if (r==0) return 0;

  //TBits shared = track->GetTPCSharedMap();
  //if (shared.CountBits()) return 0;                 // no shared clusters; pragmatic but dangerous

  if (pidi==0) return 1; 

  // pid

  if (!track->IsOn(AliESDtrack::kTPCpid)) return 0;
  Double_t p[AliPID::kSPECIES];
  track->GetTPCpid(p);
  Int_t q = tp->Charge();

  if (pidi == -211) return p[AliPID::kPion]+p[AliPID::kMuon]>0.5 && q==-1;
  else if (pidi == 211) return p[AliPID::kPion]+p[AliPID::kMuon]>0.5 && q==1;
  else if (pidi == -321) return p[AliPID::kKaon]>0.5 && q==-1;
  else if (pidi ==  321) return p[AliPID::kKaon]>0.5 && q==1;
  else if (pidi == -2212) return p[AliPID::kProton]>0.5 && q==-1;
  else if (pidi ==  2212) return p[AliPID::kProton]>0.5 && q==1;
  else if (pidi == -11) return p[AliPID::kElectron]>0.5 && q==1;
  else if (pidi ==  11) return p[AliPID::kElectron]>0.5 && q==-1;
  else return 0;
}
//=============================================================================
Bool_t AliUnicorEventAliceESD::PairGood(double p0, double the0, double phi0, double z0,  
				double p1, double the1, double phi1, double z1) const {

  // two-track separation cut

  double dthe = the1-the0;
  if (fabs(dthe)>0.010) return kTRUE;
  double B = -0.5; // magnetic field in T, needed for helix; should be gotten in the proper way
  double pt0 = p0*sin(the0);
  double pt1 = p1*sin(the1);
  double r = 1.2; // we will calculate the distance between the two tracks at r=1.2 m
  // for (double r=0.8; r<2.5; r+=0.2) {
  double si0 = -0.3*B*z0*r/2/pt0;
  double si1 = -0.3*B*z1*r/2/pt1;
  if (fabs(si0)>=1.0) return kTRUE; // could be done better
  if (fabs(si1)>=1.0) return kTRUE;
  double dphi = phi1 - phi0 + asin(si1) - asin(si0);
  dphi = TVector2::Phi_mpi_pi(dphi);
  if (fabs(dphi)<0.020) return kFALSE;
  // }
  return kTRUE;
}
//=============================================================================
