/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// Class for AOD reconstructed heavy-flavour cascades
//
// Author: X-M. Zhang, zhangxm@iopp.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"

ClassImp(AliAODRecoCascadeHF)
//-----------------------------------------------------------------------------

AliAODRecoCascadeHF::AliAODRecoCascadeHF() :
  AliAODRecoDecayHF2Prong()
{
  //
  // Default Constructor
  //
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
					 Double_t *px, Double_t *py, Double_t *pz,
					 Double_t *d0, Double_t *d0err, Double_t dca) :
  AliAODRecoDecayHF2Prong(vtx2, px, py, pz, d0, d0err, dca)
{
  //
  //  Constructor with AliAODVertex for decay vertex
  //
  SetCharge(charge);
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
					 Double_t *d0, Double_t *d0err, Double_t dca) :
  AliAODRecoDecayHF2Prong(vtx2, d0, d0err, dca)
{
  //
  //  Constructor with decay vertex and without prongs momenta
  //
  SetCharge(charge);
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(const AliAODRecoCascadeHF &source) :
  AliAODRecoDecayHF2Prong(source)
{
  //
  // Copy constructor
  //
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF &AliAODRecoCascadeHF::operator=(const AliAODRecoCascadeHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF2Prong::operator=(source);

  return *this;
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::~AliAODRecoCascadeHF()
{
  //
  // Default Destructor
  //
}
//-----------------------------------------------------------------------------
Double_t AliAODRecoCascadeHF::InvMassDstarKpipi() const 
{
  //
  // 3 prong invariant mass of the D0 daughters and the soft pion
  //
  Double_t e[3];
  if (Charge()>0){
    e[0]=Get2Prong()->EProng(0,211);
    e[1]=Get2Prong()->EProng(1,321);
  }else{
    e[0]=Get2Prong()->EProng(0,321);
    e[1]=Get2Prong()->EProng(1,211);
  }
  e[2]=EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2];
  Double_t minv = TMath::Sqrt(esum*esum-P()*P());

  return minv; 
}
//----------------------------------------------------------------------------
Int_t AliAODRecoCascadeHF::MatchToMC(Int_t pdgabs,Int_t pdgabs2prong,
                                     Int_t *pdgDg,Int_t *pdgDg2prong,
				     TClonesArray *mcArray) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle
  // 

  Int_t ndg=GetNDaughters();
  if(!ndg) {
    AliError("No daughters available");
    return -1;
  }
  
  AliAODRecoDecayHF2Prong *the2Prong = Get2Prong();
  Int_t lab2Prong = the2Prong->MatchToMC(pdgabs2prong,mcArray,2,pdgDg2prong);

  if(lab2Prong<0) return -1;

  Int_t dgLabels[10]={0,0,0,0,0,0,0,0,0,0};

  // loop on daughters and write labels
  for(Int_t i=0; i<ndg; i++) {
    AliVTrack *trk = (AliVTrack*)GetDaughter(i);
    Int_t lab = trk->GetLabel();
    if(lab==-1) { // this daughter is the 2prong
      lab=lab2Prong;
    } else if(lab<-1) {
      printf("daughter with negative label\n");
      continue;
    }
    dgLabels[i] = lab;
  }

  return AliAODRecoDecay::MatchToMC(pdgabs,mcArray,dgLabels,2,2,pdgDg);
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoCascadeHF::SelectDstar(const Double_t *cutsDstar,
					const Double_t *cutsD0,
					Bool_t testD0) const
{
  //
  // cutsDstar[0] = inv. mass half width of D* [GeV]
  // cutsDstar[1] = half width of (M_Kpipi-M_D0) [GeV]
  // cutsDstar[2] = PtMin of pi_s [GeV/c]
  // cutsDstar[3] = PtMax of pi_s [GeV/c]
  // cutsDstar[4] = theta, angle between the pi_s and decay plane of the D0 [rad]
  //
  // cutsD0[0] = inv. mass half width [GeV]   
  // cutsD0[1] = dca [cm]
  // cutsD0[2] = cosThetaStar 
  // cutsD0[3] = pTK [GeV/c]
  // cutsD0[4] = pTPi [GeV/c]
  // cutsD0[5] = d0K [cm]   upper limit!
  // cutsD0[6] = d0Pi [cm]  upper limit!
  // cutsD0[7] = d0d0 [cm^2]
  // cutsD0[8] = cosThetaPoint


  // check that the D0 passes the cuts
  // (if we have a D*+, it has to pass as D0, 
  //  if we have a D*-, it has to pass as D0bar)

  if(testD0) {
    Int_t okD0=0,okD0bar=0;
    Get2Prong()->SelectD0(cutsD0,okD0,okD0bar);
    if((Charge()==+1 && !okD0) || (Charge()==-1 && !okD0bar)) return kFALSE; 
  }
 
  if( (PtProng(0)<cutsDstar[2]) || (PtProng(0)>cutsDstar[3]) ) return kFALSE;

  Double_t mDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t invmDstar = InvMassDstarKpipi();
  if(TMath::Abs(mDstar-invmDstar)>cutsDstar[0]) return kFALSE;

  Double_t mD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if(TMath::Abs((mDstar-mD0)-DeltaInvMass())>cutsDstar[1]) return kFALSE;

  Double_t theta = AngleD0dkpPisoft(); 
  if(theta>cutsDstar[4]) return kFALSE;
  
  return kTRUE;
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoCascadeHF::SelectLctoV0(const Double_t *cutsLctoV0, 
					 Bool_t okLck0sp, Bool_t okLcLpi) const 
{
  // cuts on Lambdac candidates to V0+bachelor
  // (to be passed to AliAODRecoDecayHF3Prong::SelectLctoV0())
  // 0 = inv. mass half width in K0s hypothesis [GeV]   
  // 1 = inv. mass half width in Lambda hypothesis [GeV]   
  // 2 = inv. mass V0 in K0s hypothesis half width [GeV]   
  // 3 = inv. mass V0 in Lambda hypothesis half width [GeV]   
  // 4 = pT min Bachelor track [GeV/c]
  // 5 = pT min V0-Positive track [GeV/c]
  // 6 = pT min V0-Negative track [GeV/c]
  // 7 = dca cut on the V0 (cm)
  // 8 = dca cut on the cascade (cm)

//   if ( !Getv0() || !Getv0PositiveTrack() || !Getv0NegativeTrack() ) 
//     { AliInfo(Form("Not adapted for ESDv0s, return true...")); return false; }

  Double_t mLck0sp,mLcLpi;
  okLck0sp=1; okLcLpi=1;
  
  Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  Double_t mk0sPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  Double_t mLPDG = TDatabasePDG::Instance()->GetParticle(3122)->Mass();

  // k0s + p
  double mk0s = Getv0()->MassK0Short();
  mLck0sp = InvMassLctoK0sP();

  // lambda + pi 
  double mlambda = Getv0()->MassLambda();
  double malambda = Getv0()->MassAntiLambda();
  mLcLpi = InvMassLctoLambdaPi();

  // cut on Lc mass
  //   with k0s p hypothesis
  if(TMath::Abs(mLck0sp-mLcPDG)>cutsLctoV0[0]) okLck0sp = 0;
  //   with Lambda pi hypothesis
  if(TMath::Abs(mLcLpi-mLcPDG)>cutsLctoV0[1]) okLcLpi = 0;
  
  // cuts on the v0 mass
  if(TMath::Abs(mk0s-mk0sPDG)>cutsLctoV0[2]) okLck0sp = 0;
  if( TMath::Abs(mlambda-mLPDG)>cutsLctoV0[3] && 
      TMath::Abs(malambda-mLPDG)>cutsLctoV0[3] ) okLcLpi = 0;
  
  if(!okLck0sp && !okLcLpi) return 0;
  
  // cuts on the minimum pt of the tracks 
  if(TMath::Abs(GetBachelor()->Pt()) < cutsLctoV0[4]) return 0;
  if(TMath::Abs(Getv0PositiveTrack()->Pt()) < cutsLctoV0[5]) return 0;
  if(TMath::Abs(Getv0NegativeTrack()->Pt()) < cutsLctoV0[6]) return 0;
  
  // cut on the v0 dca
  if(TMath::Abs(Getv0()->DcaV0Daughters()) > cutsLctoV0[7]) return 0;
  
  // cut on the cascade dca
  if( TMath::Abs(GetDCA(0))>cutsLctoV0[8] ||
      TMath::Abs(Getv0()->DcaPosToPrimVertex())>cutsLctoV0[8] ||
      TMath::Abs(Getv0()->DcaNegToPrimVertex())>cutsLctoV0[8] ) return 0;
  
  return true; 

}
//-----------------------------------------------------------------------------
Double_t AliAODRecoCascadeHF::AngleD0dkpPisoft() const {
  //
  // Angle of soft pion to D0 decay plane
  // 

  TVector3 p3Trk0(Get2Prong()->PxProng(0),Get2Prong()->PyProng(0),Get2Prong()->PzProng(0)); // from D0
  TVector3 p3Trk1(Get2Prong()->PxProng(1),Get2Prong()->PyProng(1),Get2Prong()->PzProng(1)); // from D0
  TVector3 p3Trk2(PxProng(0),PyProng(0),PzProng(0)); // pi_s

  TVector3 perp = p3Trk0.Cross(p3Trk1);
  Double_t theta = p3Trk2.Angle(perp);
  if(theta>(TMath::Pi()-theta)) theta = TMath::Pi() - theta;
  theta = TMath::Pi()/2. - theta;

  return theta;
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoCascadeHF::TrigonometricalCut() const {
  //  
  // Trigonometrical constraint
  //
  TVector3 p3Trk0(Get2Prong()->PxProng(0),Get2Prong()->PyProng(0),Get2Prong()->PzProng(0)); // from D0
  TVector3 p3Trk1(Get2Prong()->PxProng(1),Get2Prong()->PyProng(1),Get2Prong()->PzProng(1)); // from D0
  TVector3 p3Trk2(PxProng(0),PyProng(0),PzProng(0)); // pi_s

  Double_t alpha = p3Trk0.Angle(p3Trk2);
  Double_t beta = p3Trk1.Angle(p3Trk2);

  Double_t cosphi01 = TMath::Cos(alpha) / TMath::Cos(AngleD0dkpPisoft());
  Double_t cosphi02 = TMath::Cos(beta) / TMath::Cos(AngleD0dkpPisoft());

  Double_t phi01 = TMath::ACos(cosphi01);
  Double_t phi02 = TMath::ACos(cosphi02);
  Double_t phi00 = p3Trk0.Angle(p3Trk1);

  if((phi01>phi00) || (phi02>phi00)) return kFALSE;
  return kTRUE;
}
