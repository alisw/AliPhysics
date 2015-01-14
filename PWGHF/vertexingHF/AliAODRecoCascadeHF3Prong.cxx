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
// Class for AOD reconstructed heavy-flavour cascades 3prong
// Used for Xic->pi Xi pi analysis
//
// Author: Y.S. Watanabe, wyosuke@cns.s.u-tokyo.ac.jp
/////////////////////////////////////////////////////////////

#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TClonesArray.h>
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF3Prong.h"

ClassImp(AliAODRecoCascadeHF3Prong)

//-----------------------------------------------------------------------------

AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong() :
AliAODRecoDecayHF3Prong()
{
  //
  // Default Constructor
  //
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong(AliAODVertex *vtx2, Short_t charge,
						     Double_t *px, Double_t *py, Double_t *pz,
						     Double_t *d0, Double_t *d0err, 
						     Double_t *dca, Double_t sigvert,
						     Double_t dist12,Double_t dist23):
  AliAODRecoDecayHF3Prong(vtx2, px, py, pz, d0, d0err, dca,sigvert,dist12,dist23,charge)
{
  //
  //  Constructor with AliAODVertex for decay vertex
  //
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong(const AliAODRecoCascadeHF3Prong &source) :
  AliAODRecoDecayHF3Prong(source)
{
  //
  // Copy constructor
  //
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong &AliAODRecoCascadeHF3Prong::operator=(const AliAODRecoCascadeHF3Prong &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF3Prong::operator=(source);

  return *this;
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::~AliAODRecoCascadeHF3Prong()
{
  //
  // Default Destructor
  //
}

//----------------------------------------------------------------------------
//Int_t AliAODRecoCascadeHF3Prong::MatchToMC(Int_t pdgabs,Int_t pdgabs3prong,
//                                     Int_t *pdgDg,Int_t *pdgDg3prong,
//				     TClonesArray *mcArray, Bool_t isV0) const
//{
//  //
//  // Check if this candidate is matched to a MC signal
//  // If no, return -1
//  // If yes, return label (>=0) of the AliAODMCParticle
//  // 
//
//  Int_t ndg=GetNDaughters();
//  if(ndg==0) {
//    AliError("No daughters available");
//    return -1;
//  }
//
//  if ( isV0 &&
//       ( (pdgDg[1]==2212 && pdgDg[0]==310) ||
//	 (pdgDg[1]==211 && pdgDg[0]==3122) ) ) {
//    AliWarning(Form("Please, pay attention: first element in AliAODRecoCascadeHF3Prong object must be the bachelor and second one V0. Skipping! (pdgDg[0] = %d, (pdgDg[1] = %d)", pdgDg[0], pdgDg[1]));
//    return -1;
//  }
//
//  Int_t lab3Prong = -1;
//
//  if (!isV0) {
//    AliAODRecoDecayHF2Prong *the2Prong = Get2Prong();
//    lab2Prong = the2Prong->MatchToMC(pdgabs2prong,mcArray,2,pdgDg2prong);
//  } else {
//    AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(Getv0());
//    lab2Prong = theV0->MatchToMC(pdgabs2prong,mcArray,2,pdgDg2prong); // the V0
//  }
//
//  if(lab2Prong<0) return -1;
//
//  Int_t dgLabels[10]={0,0,0,0,0,0,0,0,0,0};
//
//  if (!isV0) {
//    // loop on daughters and write labels
//    for(Int_t i=0; i<ndg; i++) {
//      AliVTrack *trk = dynamic_cast<AliVTrack*>(GetDaughter(i));
//      if(!trk) continue;
//      Int_t lab = trk->GetLabel();
//      if(lab==-1) { // this daughter is the 2prong
//	lab=lab2Prong;
//      } else if(lab<-1) continue;
//      dgLabels[i] = lab;
//    }
//  } else {
//    AliVTrack *trk = dynamic_cast<AliVTrack*>(GetBachelor()); // the bachelor
//    if (!trk) return -1;
//    dgLabels[0] = trk->GetLabel();//TMath::Abs(trk->GetLabel());
//    dgLabels[1] = lab2Prong;
//  }
//
//  Int_t finalLabel = AliAODRecoDecay::MatchToMC(pdgabs,mcArray,dgLabels,2,2,pdgDg);
//
//  if (finalLabel>=0){
//    // debug printouts for Lc->V0 bachelor case
//
//    if ( isV0 && (dgLabels[0]!=-1 && dgLabels[1]!=-1) ) {
//      AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(Getv0());
//      Bool_t onTheFly = theV0->GetOnFlyStatus();
//      if (pdgDg[0]==2212 && pdgDg[1]==310) {
//	AliAODMCParticle*k0s = dynamic_cast<AliAODMCParticle*>(mcArray->At(lab2Prong));
//	if(k0s){
//	  Int_t labK0 = k0s->GetMother();	
//	  AliAODMCParticle*k0bar = dynamic_cast<AliAODMCParticle*>(mcArray->At(labK0));
//	  if(k0bar){
//	    AliDebug(1,Form(" (onTheFly=%1d) LabelV0=%d (%d) -> LabelK0S=%d (%d -> %d %d)",onTheFly,labK0,k0bar->GetPdgCode(),lab2Prong,pdgabs2prong,pdgDg2prong[0],pdgDg2prong[1]));
//	    AliDebug(1,Form(" LabelLc=%d (%d) -> LabelBachelor=%d (%d) LabelV0=%d (%d)",
//			    finalLabel,pdgabs,
//			    dgLabels[0],pdgDg[0],dgLabels[1],pdgDg[1]));
//	  }
//	}
//      } else if (pdgDg[0]==211 && pdgDg[1]==3122) {
//	AliDebug(1,Form(" (onTheFly=%1d) LabelV0=%d (%d -> %d %d)",onTheFly,lab2Prong,pdgabs2prong,pdgDg2prong[0],pdgDg2prong[1]));
//	AliDebug(1,Form(" LabelLc=%d (%d) -> LabelBachelor=%d (%d) LabelV0=%d (%d)",
//			finalLabel,pdgabs,
//		      dgLabels[0],pdgDg[0],dgLabels[1],pdgDg[1]));
//      }
//
//    }
//  }
//
//  return finalLabel;
//
//}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaXiDaughters() const
{
  //
  // DCA between Xi daughters
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaXiDaughters();

}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaV0Daughters() const
{
  //
  // DCA between Cascade-V0 daughters
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaV0Daughters();

}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDecayLength() const
{
  //
  // Decay length of Xi
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();

  if (!casc) 
    return -1.;
  AliAODVertex *vtxPrimary = GetPrimaryVtx();
  Double_t posVtx[3] = {0.,0.,0.};
  vtxPrimary->GetXYZ(posVtx);
  return casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);

}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDecayLengthV0() const
{
  //
  // Decay length of V0 from Xi
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DecayLengthV0();

}
//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascCosPointingAngle() const 
{
  //
  // Xi pointing angle to primary vertex
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -999.;

  AliAODVertex *vtxPrimary = GetPrimaryVtx();
  Double_t posVtx[3] = {0.,0.,0.};
  vtxPrimary->GetXYZ(posVtx);
  return casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
}
//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascCosPointingAngleV0() const 
{
  //
  // Cos pointing angle of V0 to Xi decay vertex
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -999.;
  return casc->CosPointingAngle(casc->GetDecayVertexXi());
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaV0ToPrimVertex() const
{
  //
  // DCA to primary vertex of Cascade-V0 
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaV0ToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaPosToPrimVertex() const
{
  //
  // DCA to primary vertex of Cascade-positive track
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaPosToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaNegToPrimVertex() const
{
  //
  // DCA to primary vertex of Cascade-negative track
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaNegToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaBachToPrimVertex() const
{
  //
  // DCA to primary vertex of Cascade-Bachelor track
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaBachToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassXi() const
{
  //
  // Xi mass
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassXi();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassLambda() const
{
  //
  // Lambda mass of cascade-v0
  //

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassLambda();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassAntiLambda() const
{
  //
  // Anti-Lambda mass of cascade-v0
  //
	
  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassAntiLambda();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::XicCosPointingAngle() const
{
  //
  // Xic pointing angle to primary vertex
  //

  AliAODVertex *vtxPrimary = GetPrimaryVtx();
  AliAODVertex *vtxSecondary = GetSecondaryVtx();

  Double_t dx = vtxSecondary->GetX()-vtxPrimary->GetX();
  Double_t dy = vtxSecondary->GetY()-vtxPrimary->GetY();
  Double_t dl = sqrt(dx*dx+dy*dy);

  Double_t px = Px();
  Double_t py = Py();
  Double_t pt = Pt();

  if(dl>0&&pt>0)
    return (px*dx+py*dy)/pt/dl; 
  else
    return -9999.;
}
