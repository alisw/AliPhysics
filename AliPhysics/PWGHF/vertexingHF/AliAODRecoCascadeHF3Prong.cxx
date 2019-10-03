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

/// \cond CLASSIMP
ClassImp(AliAODRecoCascadeHF3Prong);
/// \endcond

//-----------------------------------------------------------------------------

AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong() :
AliAODRecoDecayHF3Prong()
{
  ///
  /// Default Constructor
  ///
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong(AliAODVertex *vtx2, Short_t charge,
						     Double_t *px, Double_t *py, Double_t *pz,
						     Double_t *d0, Double_t *d0err, 
						     Double_t *dca, Double_t sigvert,
						     Double_t dist12,Double_t dist23):
  AliAODRecoDecayHF3Prong(vtx2, px, py, pz, d0, d0err, dca,sigvert,dist12,dist23,charge)
{
  ///
  ///  Constructor with AliAODVertex for decay vertex
  ///
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::AliAODRecoCascadeHF3Prong(const AliAODRecoCascadeHF3Prong &source) :
  AliAODRecoDecayHF3Prong(source)
{
  ///
  /// Copy constructor
  ///
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong &AliAODRecoCascadeHF3Prong::operator=(const AliAODRecoCascadeHF3Prong &source)
{
  ///
  /// assignment operator
  ///
  if(&source == this) return *this;

  AliAODRecoDecayHF3Prong::operator=(source);

  return *this;
}

//-----------------------------------------------------------------------------
AliAODRecoCascadeHF3Prong::~AliAODRecoCascadeHF3Prong()
{
  ///
  /// Default Destructor
  ///
}


//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaXiDaughters() const
{
  ///
  /// DCA between Xi daughters
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaXiDaughters();

}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaV0Daughters() const
{
  ///
  /// DCA between Cascade-V0 daughters
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaV0Daughters();

}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDecayLength() const
{
  ///
  /// Decay length of Xi
  ///

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
  ///
  /// Decay length of V0 from Xi
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DecayLengthV0();

}
//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascCosPointingAngle() const 
{
  ///
  /// Xi pointing angle to primary vertex
  ///

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
  ///
  /// Cos pointing angle of V0 to Xi decay vertex
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -999.;
  return casc->CosPointingAngle(casc->GetDecayVertexXi());
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaV0ToPrimVertex() const
{
  ///
  /// DCA to primary vertex of Cascade-V0
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaV0ToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaPosToPrimVertex() const
{
  ///
  /// DCA to primary vertex of Cascade-positive track
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaPosToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaNegToPrimVertex() const
{
  ///
  /// DCA to primary vertex of Cascade-negative track
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaNegToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascDcaBachToPrimVertex() const
{
  ///
  /// DCA to primary vertex of Cascade-Bachelor track
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->DcaBachToPrimVertex();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassXi() const
{
  ///
  /// Xi mass
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassXi();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassLambda() const
{
  ///
  /// Lambda mass of cascade-v0
  ///

  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassLambda();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::CascMassAntiLambda() const
{
  ///
  /// Anti-Lambda mass of cascade-v0
  ///
	
  AliAODcascade *casc = (AliAODcascade*)GetCascade();
  if (!casc) 
    return -1.;
  return casc->MassAntiLambda();
}

//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::XicCosPointingAngle() const
{
  ///
  /// Xic pointing angle to primary vertex
  ///

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
//________________________________________________________________________
Double_t AliAODRecoCascadeHF3Prong::BachelorsCosPointingAngle() const
{
  ///
  /// Bachelor pointing angle to primary vertex
  ///

  AliAODVertex *vtxPrimary = GetPrimaryVtx();
  AliAODVertex *vtxSecondary = GetSecondaryVtx();

  Double_t dx = vtxSecondary->GetX()-vtxPrimary->GetX();
  Double_t dy = vtxSecondary->GetY()-vtxPrimary->GetY();
  Double_t dl = sqrt(dx*dx+dy*dy);

  Double_t px = PxProng(0)+PxProng(2);
  Double_t py = PyProng(0)+PyProng(2);
  Double_t pt = sqrt(px*px+py*py);

  if(dl>0&&pt>0)
    return (px*dx+py*dy)/pt/dl; 
  else
    return -9999.;
}
//----------------------------------------------------------------------------
Int_t AliAODRecoCascadeHF3Prong::MatchToMC(Int_t pdgabs,Int_t pdgabscasc,
                                     Int_t *pdgDg,Int_t *pdgDgcasc,Int_t *pdgDgv0
				     ,TClonesArray *mcArray) const
{
  ///
  /// Check if this candidate is matched to a MC signal
  /// If no, return -1
  /// If yes, return label (>=0) of the AliAODMCParticle
  ///

  Int_t ndg=GetNDaughters();
  if(ndg==0) {
    AliError("No daughters available");
    return -1;
  }

  if ( pdgabs!=4232 || pdgDg[0]!=211 || pdgDg[1]!=3312 || pdgDg[2]!=211 ) 
	{
    AliWarning(Form("Please, pay attention: Only pi Xi pi decay is supported now"));
    return -1;
  }

  AliAODcascade *theCascade = dynamic_cast<AliAODcascade*>(GetCascade());
	if(!theCascade) return -1;
  AliAODTrack *trk1 = dynamic_cast<AliAODTrack*>(GetBachelor1()); // the bachelor
  if (!trk1) return -1;
  AliAODTrack *trk2 = dynamic_cast<AliAODTrack*>(GetBachelor2()); // the bachelor
  if (!trk2) return -1;

  Int_t labcasc = MatchToMCCascade(theCascade,pdgabscasc,pdgDgcasc,pdgDgv0,mcArray); // the cascade
  if(labcasc<0) return -1;
	Int_t labtrk1 = trk1->GetLabel();
	if(labtrk1<0) return -1;
	Int_t labtrk2 = trk2->GetLabel();
	if(labtrk2<0) return -1;

  Int_t dgLabels[10]={0,0,0,0,0,0,0,0,0,0};

  dgLabels[0] = labtrk1;
  dgLabels[1] = labcasc;
  dgLabels[2] = labtrk2;

  Int_t finalLabel = MatchToMCXicPlus(pdgabs,mcArray,dgLabels,3,3,pdgDg);

  return finalLabel;

}
//________________________________________________________________________
Int_t AliAODRecoCascadeHF3Prong::MatchToMCCascade(AliAODcascade *theCascade, Int_t pdgabscasc, Int_t *pdgDgcasc, Int_t *pdgDgv0, TClonesArray *mcArray) const // the cascade
{

	AliAODTrack *cptrack = (AliAODTrack*) theCascade->GetDaughter(0);
	if(!cptrack) return -1;
	Int_t label_p = cptrack->GetLabel();
	if(label_p<0) return -1;
	AliAODTrack *cntrack = (AliAODTrack*) theCascade->GetDaughter(1);
	if(!cntrack) return -1;
	Int_t label_n = cntrack->GetLabel();
	if(label_n<0) return -1;
	Int_t labv0 = theCascade->MatchToMC(pdgDgcasc[1],mcArray,2,pdgDgv0);
	if(labv0<0) return -1;
	AliAODMCParticle *mcpartv0= (AliAODMCParticle*) mcArray->At(labv0);

	AliAODTrack *cbtrack = (AliAODTrack*) theCascade->GetDecayVertexXi()->GetDaughter(0);
	if(!cbtrack) return -1;

	Int_t label_b = cbtrack->GetLabel();
	if(label_b<0) return -1;

	AliAODMCParticle *mcpartb= (AliAODMCParticle*) mcArray->At(label_b);
	Int_t pdgb = TMath::Abs(mcpartb->GetPdgCode());
	if(pdgb!=pdgDgcasc[0]) return -1;

	AliAODMCParticle *mcmotherv0=mcpartv0;
	Bool_t isFromXiv0 = kFALSE;
	Int_t labxiv0 = mcmotherv0->GetMother();
	if(labxiv0<0) return -1;
	mcmotherv0 =  (AliAODMCParticle*) mcArray->At(labxiv0);
	if(mcmotherv0){
		Int_t pdg = TMath::Abs(mcmotherv0 ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXiv0 = kTRUE;
		}
	}
	if(!isFromXiv0) return -1;

	AliAODMCParticle *mcmotherb=mcpartb;
	Bool_t isFromXib = kFALSE;
	Int_t labxib = mcmotherb->GetMother();
	if(labxib<0) return -1;
	mcmotherb =  (AliAODMCParticle*) mcArray->At(labxib);
	if(mcmotherb){
		Int_t pdg = TMath::Abs(mcmotherb ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXib = kTRUE;
		}
	}
	if(!isFromXib) return -1;

	if(labxiv0!=labxib) return -1;//Bachelor and V0 should come from the same Xi

	return labxib;
}
//----------------------------------------------------------------------------
Int_t AliAODRecoCascadeHF3Prong::MatchToMCXicPlus(Int_t pdgabs,TClonesArray *mcArray,
				 Int_t dgLabels[10],Int_t ndg,
				 Int_t ndgCk, const Int_t *pdgDg) const
{
  ///
  /// Check if this candidate is matched to a MC signal
  /// If no, return -1
  /// If yes, return label (>=0) of the AliAODMCParticle
  ///
  ///

  Int_t labMom[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t i,j,lab,labMother,pdgMother,pdgPart;
  AliAODMCParticle *part=0;
  AliAODMCParticle *mother=0;
  Bool_t pdgUsed[10]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};

  // loop on daughter labels
  for(i=0; i<ndg; i++) {
    labMom[i]=-1;
    lab = TMath::Abs(dgLabels[i]);
    if(lab<0) {
      printf("daughter with negative label %d\n",lab);
      return -1;
    }
    part = (AliAODMCParticle*)mcArray->At(lab);
    if(!part) { 
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter, if requested
    if(ndgCk>0) {
      pdgPart=TMath::Abs(part->GetPdgCode());
      for(j=0; j<ndg; j++) {
	if(!pdgUsed[j] && pdgPart==pdgDg[j]) {
	  pdgUsed[j]=kTRUE;
	  break;
	}
      }
    }

    mother = part;
    while(mother->GetMother()>=0) {
      labMother=mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(labMother);
      if(!mother) {
	printf("no MC mother particle\n");
	break;
      }
      pdgMother = TMath::Abs(mother->GetPdgCode());
      if(pdgMother==pdgabs) {
	labMom[i]=labMother;
	break;
      } else if(pdgMother>pdgabs || pdgMother<10) {
	break;
      }
    }
    if(labMom[i]==-1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters

  // check if the candidate is signal
  labMother=labMom[0];
  // all labels have to be the same and !=-1
  for(i=0; i<ndg; i++) {
    if(labMom[i]==-1)        return -1;
    if(labMom[i]!=labMother) return -1;
  }

  // check that all daughter PDGs are matched
  if(ndgCk>0) {
    for(i=0; i<ndg; i++) {
      if(pdgUsed[i]==kFALSE) return -1;
    }
  }
 
  return labMother;
}
