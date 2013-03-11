/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and steps
// For Lc->V0+bachelor
// 
// - Based on AliCFVertexingHFCascade -
//
// Contact : A.De Caro - decaro@sa.infn.it
//           Centro 'E.Fermi' - Rome (Italy)
//           INFN and University of Salerno (Italy)
//
//-----------------------------------------------------------------------

#include "TDatabasePDG.h"
#include "TClonesArray.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliCFTaskVertexingHF.h"
#include "AliCFContainer.h"
#include "AliCFVertexingHF.h"
#include "AliCFVertexingHFLctoV0bachelor.h"

#include <Riostream.h>

using std::cout;
using std::endl;

ClassImp(AliCFVertexingHFLctoV0bachelor)

//_________________________________________
AliCFVertexingHFLctoV0bachelor::AliCFVertexingHFLctoV0bachelor():
fGenLcOption(0)
{
  // standard constructor
}

//_____________________________________
AliCFVertexingHFLctoV0bachelor::AliCFVertexingHFLctoV0bachelor(TClonesArray *mcArray, UShort_t originDselection, Int_t lcDecay):
AliCFVertexingHF(mcArray, originDselection),
  fGenLcOption(lcDecay)
{
  // standard constructor

  SetNProngs(3);
  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  for(Int_t iP=0; iP<fProngs; iP++){
    fPtAccCut[iP]=0.1;
    fEtaAccCut[iP]=0.9;
  }

}


//_____________________________________
AliCFVertexingHFLctoV0bachelor& AliCFVertexingHFLctoV0bachelor::operator=(const AliCFVertexingHFLctoV0bachelor& c)
{
  // operator =
 
  if  (this != &c) {
    AliCFVertexingHF::operator=(c);
  }

  return *this;

}

//__________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::SetRecoCandidateParam(AliAODRecoDecayHF *recoCand)
{
  // set the AliAODRecoDecay candidate
  
  Bool_t bSignAssoc = kFALSE;
 
  fRecoCandidate = recoCand;
  if (!fRecoCandidate) {
    AliError("fRecoCandidate not found, problem in assignement\n");
    return bSignAssoc;
  }
  
  if (fRecoCandidate->GetPrimaryVtx()) AliDebug(4,"fReco Candidate has a pointer to PrimVtx\n");
  
  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;
  if ( !(lcV0bachelor->Getv0()) ) {
    AliDebug(1,"It is not a Lc->V0+bachelor candidate");
    return bSignAssoc;
  }

  Int_t pdgCand = 4122;
  Int_t mcLabel = -1;
  Int_t mcLabelK0S = -1;
  Int_t mcLabelLambda = -1;

  // Lc->K0S+p and cc
  Int_t pdgDgLctoV0bachelor[2]={2212,310}; // first bachelor, second V0
  Int_t pdgDgV0toDaughters[2]={211,211};
  mcLabelK0S = lcV0bachelor->MatchToMC(pdgCand,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,fmcArray,kTRUE);

  // Lc->Lambda+pi and cc
  pdgDgLctoV0bachelor[0]=211, pdgDgLctoV0bachelor[1]=3122; // first bachelor, second V0
  pdgDgV0toDaughters[0]=2212,  pdgDgV0toDaughters[1]=211;
  mcLabelLambda = lcV0bachelor->MatchToMC(pdgCand,pdgDgLctoV0bachelor[1],pdgDgLctoV0bachelor,pdgDgV0toDaughters,fmcArray,kTRUE);

  if (mcLabelK0S!=-1 && mcLabelLambda!=-1)
    AliDebug(2,"Strange: current Lc->V0+bachelor candidate has two MC different labels!");

  if (fGenLcOption==kCountAllLctoV0) {
    if (mcLabelK0S!=-1) mcLabel=mcLabelK0S;
    if (mcLabelLambda!=-1) mcLabel=mcLabelLambda;
  }
  else if (fGenLcOption==kCountK0Sp) {
    if (mcLabelK0S!=-1) mcLabel=mcLabelK0S;
    if (mcLabelLambda!=-1) {
      mcLabel=-1;
      fFake = 0;    // fake candidate
      if (fFakeSelection==1) return bSignAssoc;
    }
  }
  else if (fGenLcOption==kCountLambdapi) {
    if (mcLabelLambda!=-1) mcLabel=mcLabelLambda;
    if (mcLabelK0S!=-1) {
      mcLabel=-1;
      fFake = 0;    // fake candidate
      if (fFakeSelection==1) return bSignAssoc;
    }
  }

  if (mcLabel==-1) {
    AliDebug(4,"No mcLabel found for current candidate");
    return bSignAssoc;
  }
  AliDebug(1,Form("Found mcLabel (%d) for current candidate",mcLabel));

  if (fRecoCandidate->NumberOfFakeDaughters()>0){
    fFake = 0;    // fake candidate
    if (fFakeSelection==1) return bSignAssoc;
  }
  if (fRecoCandidate->NumberOfFakeDaughters()==0){
    fFake = 2;    // non-fake candidate
    if (fFakeSelection==2) return bSignAssoc;
  }
  
  SetMCLabel(mcLabel); // fmcLabel=mcLabel
  fmcPartCandidate = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fmcLabel)); 
  if (!fmcPartCandidate){
    AliDebug(3,"No MC object for current candidate");
    return bSignAssoc;
  }
 
  bSignAssoc = kTRUE;
  return bSignAssoc;

}

//______________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::GetGeneratedValuesFromMCParticle(Double_t* vectorMC) 
{
  // 
  // collecting all the necessary info (pt, y, invMassV0, cosPAwrtPVxV0, onTheFlyStatusV0) from MC particle
  // (additional infos: pTbachelor, pTV0pos, pTV0neg, phi, dcaV0, cTV0, cT, cosPA)
  //

  Bool_t bGenValues = kFALSE;

  if (fmcPartCandidate->GetNDaughters()!=2) {
    AliDebug(2,"Lc MC particle doesn't decay in 2 daughters");
    return bGenValues;
  }

  Int_t daughter0lc = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1lc = fmcPartCandidate->GetDaughter(1);
  if (daughter0lc<0 || daughter1lc<0) {
    AliDebug(2,"Lc daughters are not in MC array");
    return bGenValues;
  }

  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0lc));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1lc));
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug(2,"Problems in the MC Daughters\n");
    return bGenValues;
  }

  if ( fGenLcOption==kCountLambdapi &&
       !(TMath::Abs(mcPartDaughter0->GetPdgCode())==3122 &&
	 TMath::Abs(mcPartDaughter1->GetPdgCode())==211) &&
       !(TMath::Abs(mcPartDaughter1->GetPdgCode())==3122 &&
	 TMath::Abs(mcPartDaughter0->GetPdgCode())==211) ) return bGenValues;
  if ( fGenLcOption==kCountK0Sp &&
       !(TMath::Abs(mcPartDaughter0->GetPdgCode())==2212 &&
	 TMath::Abs(mcPartDaughter1->GetPdgCode())==311) &&
       !(TMath::Abs(mcPartDaughter1->GetPdgCode())==2212 &&
	 TMath::Abs(mcPartDaughter0->GetPdgCode())==311) ) return bGenValues;

  if ( (TMath::Abs(mcPartDaughter0->GetPdgCode())==311   &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==2212) ||
       (TMath::Abs(mcPartDaughter0->GetPdgCode())==3122  &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==211) )
    bGenValues = FillVectorFromMCarray(mcPartDaughter1,mcPartDaughter0,vectorMC);
  else if ( (TMath::Abs(mcPartDaughter1->GetPdgCode())==311   &&
	     TMath::Abs(mcPartDaughter0->GetPdgCode())==2212) ||
	    (TMath::Abs(mcPartDaughter1->GetPdgCode())==3122  &&
	     TMath::Abs(mcPartDaughter0->GetPdgCode())==211) )
    bGenValues = FillVectorFromMCarray(mcPartDaughter0,mcPartDaughter1,vectorMC);

  if (!bGenValues)
    AliDebug(2,"There is something wrong in filling MC vector");

  return bGenValues;

}

//____________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{ 
  // read the variables for the container

  Bool_t bFillRecoValues = kFALSE;

  //Get the Lc and the V0 from Lc
  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;

  AliAODTrack* bachelor = (AliAODTrack*)lcV0bachelor->GetBachelor();
  AliAODv0* v0toDaughters = (AliAODv0*)lcV0bachelor->Getv0();
  if (!lcV0bachelor || !bachelor || !v0toDaughters) {
    AliDebug(2,"No V0 or bachelor in this reco candidate, skipping!");
    return bFillRecoValues;
  }

  Bool_t onTheFlyStatus = v0toDaughters->GetOnFlyStatus();
  AliAODTrack* v0positiveTrack = (AliAODTrack*)lcV0bachelor->Getv0PositiveTrack();
  AliAODTrack* v0negativeTrack = (AliAODTrack*)lcV0bachelor->Getv0NegativeTrack();
  if (!v0positiveTrack || !v0negativeTrack) {
    AliDebug(2,"No V0daughters in this reco candidate, skipping!");
    return bFillRecoValues;
  }

  Double_t pt = lcV0bachelor->Pt();
  Double_t rapidity = lcV0bachelor->Y(4122);

  Double_t cosPAwrtPrimVtxV0 = lcV0bachelor->CosV0PointingAngle();

  Double_t pTbachelor = bachelor->Pt();
  Double_t pTV0pos = v0positiveTrack->Pt();
  Double_t pTV0neg = v0negativeTrack->Pt();
  Double_t phi = lcV0bachelor->Phi();
  Double_t dcaV0 = v0toDaughters->GetDCA();
  Double_t cTLc = lcV0bachelor->Ct(4122); // wrt PrimVtx
  //Double_t dcaLc = lcV0bachelor->GetDCA();
  Double_t cosPointingAngleLc = lcV0bachelor->CosPointingAngle();

  Double_t cTV0 = 0.;
  AliAODVertex *vtx0 = (AliAODVertex*)lcV0bachelor->GetPrimaryVtx();
  if (!vtx0) {
    AliDebug(2,"Candidate has not primary vtx");
  } else {
    Double_t primVtxPos[3] = {0.,0.,0.}; vtx0->GetXYZ(primVtxPos);
    if (fGenLcOption==kCountK0Sp) {
      cTV0 = v0toDaughters->Ct(310,primVtxPos);
    } else if (fGenLcOption==kCountLambdapi) {
      cTV0 = v0toDaughters->Ct(3122,primVtxPos);
    }
  }

  Double_t invMassV0 = 0.;
  if (fGenLcOption==kCountLambdapi) {

    Short_t bachelorCharge = bachelor->Charge();
    if (bachelorCharge==1) {
      invMassV0 = v0toDaughters->MassLambda();
    } else if (bachelorCharge==-1) {
      invMassV0 = v0toDaughters->MassAntiLambda();
    }

  } else if (fGenLcOption==kCountK0Sp) {

    invMassV0 = v0toDaughters->MassK0Short();

  }

  vectorReco[0]  = pt;
  vectorReco[1]  = rapidity;
  vectorReco[2]  = phi;
  vectorReco[3]  = cosPAwrtPrimVtxV0;
  vectorReco[4]  = onTheFlyStatus;
  vectorReco[5]  = fCentValue;
  vectorReco[6]  = fFake; // whether the reconstructed candidate was a fake (fFake = 0) or not (fFake = 2) 
  vectorReco[7]  = fMultiplicity;

  if (fConfiguration==AliCFTaskVertexingHF::kSnail) {
    vectorReco[8]  = pTbachelor;
    vectorReco[9]  = pTV0pos;
    vectorReco[10] = pTV0neg;
    vectorReco[11] = invMassV0;
    vectorReco[12] = dcaV0;
    vectorReco[13] = cTV0*1.E4; // in micron
    vectorReco[14] = cTLc*1.E4; // in micron
    vectorReco[15] = cosPointingAngleLc;
  }

  bFillRecoValues = kTRUE;

  return bFillRecoValues;
}

//_____________________________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::CheckMCChannelDecay() const
{ 
  // check the required decay channel

  Bool_t checkCD = kFALSE;
  
  if (fmcPartCandidate->GetNDaughters()!=2) {
    AliDebug(2, Form("The MC particle doesn't decay in 2 particles, skipping!!"));
    return checkCD;
  }

  Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
  if (daughter0<0 || daughter1<0){
    AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
    return checkCD;
  }
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug(2,"Problems in the MC Daughters\n");
    return checkCD;
  }

  // Lc -> Lambda + pion AND cc
  if (fGenLcOption==kCountLambdapi) {

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==3122 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==211  &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==3122)) {
      AliDebug(2, "The Lc MC doesn't decay in Lambda+pion (or cc), skipping!!");
      return checkCD;  
    }

    if (TMath::Abs(mcPartDaughter0->GetPdgCode())==3122) {
      mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); // the bachelor
      mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); // the V0
    }
    if (!mcPartDaughter0 || !mcPartDaughter1) {
      AliDebug(2,"Problems in the MC Daughters\n");
      return checkCD;
    }

    if (mcPartDaughter1->GetNDaughters()!=2) {
      AliDebug(2, "The Lambda MC particle doesn't decay in 2 particles, skipping!!");
      return checkCD;
    }

    Int_t daughter1D0 = mcPartDaughter1->GetDaughter(0);
    Int_t daughter1D1 = mcPartDaughter1->GetDaughter(1);
    if (daughter1D0<0 || daughter1D1<0) {
      AliDebug(2, Form("The Lambda MC particle doesn't have correct daughters, skipping!!"));
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughter1D0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D0));
    AliAODMCParticle* mcPartDaughter1D1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D1));
    if(!mcPartDaughter1D0 || !mcPartDaughter1D1) {
      AliError("The Lambda daughter particle not found in MC array");
      return checkCD;
    }

    if (!(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==211   &&
	  TMath::Abs(mcPartDaughter1D1->GetPdgCode())==2212) &&
	!(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==2212  &&
	  TMath::Abs(mcPartDaughter1D1->GetPdgCode())==211)) {
      AliDebug(2, "The Lambda MC doesn't decay in pi+proton (or cc), skipping!!");
      return checkCD;
    }

  } else if (fGenLcOption==kCountK0Sp) { // Lc -> K0bar + proton AND cc

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==311   &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==2212) &&
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==2212  &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==311)) {
      AliDebug(2, "The Lc MC doesn't decay in K0+proton (or cc), skipping!!");
      return checkCD;  
    }
    
    if (TMath::Abs(mcPartDaughter0->GetPdgCode())==311) {
      mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); // the bachelor
      mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); // the V0
    }

    Int_t daughter = mcPartDaughter1->GetDaughter(0);
    if (daughter<0) {
      AliDebug(2, Form("The K0/K0bar MC particle doesn't have correct daughter, skipping!!"));
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter));
    if(!mcPartDaughter){
      AliError("The K0/K0bar daughter particle not found in MC array");
      return checkCD;
    }

    if (!(TMath::Abs(mcPartDaughter->GetPdgCode())==310)) {
      AliDebug(2, "The K0/K0bar MC doesn't go in K0S, skipping!!");
      return checkCD;
    }

    if (mcPartDaughter->GetNDaughters()!=2) {
      AliDebug(2, "The K0S MC doesn't decay in 2 particles, skipping!!");
      return checkCD;
    }

    Int_t daughterD0 = mcPartDaughter->GetDaughter(0);
    Int_t daughterD1 = mcPartDaughter->GetDaughter(1);
    if (daughterD0<0 || daughterD1<0) {
      AliDebug(2, Form("The K0S MC particle doesn't have correct daughters, skipping!!"));
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughterD0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD0));
    AliAODMCParticle* mcPartDaughterD1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD1));
    if (!mcPartDaughterD0 || !mcPartDaughterD1) {
      AliError("Daughter particle not found in MC array");
      return checkCD;
    }

    if (! ( TMath::Abs(mcPartDaughterD0->GetPdgCode())==211 &&
	    TMath::Abs(mcPartDaughterD1->GetPdgCode())==211 ) ) {
      AliDebug(2, "The K0S MC doesn't decay in pi+ pi-, skipping!!");
      return checkCD;
    }

  }
  
  checkCD = kTRUE;
  return checkCD;
  
}

//_____________________________________________________________
Double_t AliCFVertexingHFLctoV0bachelor::GetEtaProng(Int_t iProng) const 
{
  //
  // getting eta of the prong - overload the mother class method
  //

  Double_t etaProng =-9999;

  if (!fRecoCandidate) {
    AliDebug(2,"No reco candidate selected");
    return etaProng;
  }

  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;
  AliAODTrack* bachelor = (AliAODTrack*)lcV0bachelor->GetBachelor();
  AliAODTrack* v0Pos = (AliAODTrack*)lcV0bachelor->Getv0PositiveTrack();
  AliAODTrack* v0Neg = (AliAODTrack*)lcV0bachelor->Getv0NegativeTrack();
  if (!(lcV0bachelor->Getv0()) || !bachelor || !v0Pos || !v0Neg) {
    AliDebug(2,"No V0 for this reco candidate selected");
    return etaProng;
  }

  if (iProng==0) etaProng = bachelor->Eta();
  else if (iProng==1) etaProng = v0Pos->Eta();
  else if (iProng==2) etaProng = v0Neg->Eta();

  AliDebug(4,Form("Eta value for prong number %1d = %f",iProng,etaProng));

  return etaProng;

}

//_____________________________________________________________

Double_t AliCFVertexingHFLctoV0bachelor::GetPtProng(Int_t iProng) const 
{
  //
  // getting pt of the prong
  //

  Double_t ptProng=-9999.;

  if (!fRecoCandidate) {
    AliDebug(2,"No reco candidate selected");
    return ptProng;
  }

  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;
  AliAODTrack* bachelor = (AliAODTrack*)lcV0bachelor->GetBachelor();
  AliAODTrack* v0Pos = (AliAODTrack*)lcV0bachelor->Getv0PositiveTrack();
  AliAODTrack* v0Neg = (AliAODTrack*)lcV0bachelor->Getv0NegativeTrack();
  if (!(lcV0bachelor->Getv0()) || !bachelor || !v0Pos || !v0Neg) {
    AliDebug(2,"No V0 for this reco candidate selected");
    return ptProng;
  }

  if (iProng==0) ptProng = bachelor->Pt();
  else if (iProng==1) ptProng = v0Pos->Pt();
  else if (iProng==2) ptProng = v0Neg->Pt();
    
  AliDebug(4,Form("Pt value for prong number %1d = %f",iProng,ptProng));

  return ptProng;
  
}

//_____________________________________________________________

Double_t AliCFVertexingHFLctoV0bachelor::Ctau(AliAODMCParticle *mcPartCandidate)
{

  Double_t cTau = 999999.;

  Int_t daughterD0 = mcPartCandidate->GetDaughter(0);
  Int_t daughterD1 = mcPartCandidate->GetDaughter(1);
  if (daughterD0<0 || daughterD1<0) {
    AliDebug(2, Form("The Lc MC particle doesn't have correct daughters, skipping!!"));
    return cTau;
  }

  AliAODMCParticle *mcPartDaughter0 = (AliAODMCParticle*)fmcArray->At(mcPartCandidate->GetDaughter(0));
  AliAODMCParticle *mcPartDaughter1 = (AliAODMCParticle*)fmcArray->At(mcPartCandidate->GetDaughter(1));
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug(2,"The candidate daughter particles not found in MC array");
    return cTau;
  }

  Double_t vtx1[3] = {0,0,0};   // primary vertex		
  Bool_t hasProdVertex = mcPartCandidate->XvYvZv(vtx1);  // cm

  Double_t vtx1daughter[3] = {0,0,0};   // secondary vertex
  Bool_t v0Vertex = mcPartDaughter0->XvYvZv(vtx1daughter);  //cm
  Double_t vtx2daughter[3] = {0,0,0};   // secondary vertex
  Bool_t bachVertex = hasProdVertex && mcPartDaughter1->XvYvZv(vtx2daughter);  //cm

  if (!hasProdVertex || !v0Vertex || !bachVertex) {
    AliDebug(2,"At least one of Prim.vtx, V0vtx, BachelorVtx doesn't exist!");
    return cTau;
  }

  if (TMath::Abs(vtx1daughter[0]-vtx2daughter[0])>1E-5 ||
      TMath::Abs(vtx1daughter[1]-vtx2daughter[1])>1E-5 ||
      TMath::Abs(vtx1daughter[2]-vtx2daughter[2])>1E-5) {
    AliDebug(2,"Bachelor and V0 haven't common vtx!");
    return cTau;
  }

  Double_t decayLength = 0.;
  for (Int_t ii=0; ii<3; ii++) decayLength += (vtx1daughter[ii]-vtx1[ii])*(vtx1daughter[ii]-vtx1[ii]);
  decayLength = TMath::Sqrt(decayLength);

  cTau = decayLength * mcPartCandidate->M()/mcPartCandidate->P();

  AliDebug(2,Form(" cTau(4122)=%f",cTau));

  return cTau;

}

//------------
Bool_t AliCFVertexingHFLctoV0bachelor::SetLabelArray()
{
  //
  // setting the label arrays
  //

  Bool_t checkCD = kFALSE;
  
  if (fmcPartCandidate->GetNDaughters()!=2) {
    AliDebug(2, Form("The MC particle doesn't have 2 daughters, skipping!!"));
    return checkCD;
  }

  Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
  if (daughter0<0 || daughter1<0){
    AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
    return checkCD;
  }

  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug(2,"Problems in the MC Daughters\n");
    return checkCD;
  }


  fLabelArray = new Int_t[fProngs];

  if (fGenLcOption==kCountLambdapi) { // Lc -> Lambda + pion OR cc

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==3122 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==211  &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==3122)) {
      AliDebug(2, "The Lc MC doesn't decay in Lambda+pion (or cc), skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;  
    }

    // it is Lc -> Lambda + pion OR cc
    if (TMath::Abs(mcPartDaughter0->GetPdgCode())==3122) {
      mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); // the bachelor
      mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); // the V0
      Int_t daughterTemp = daughter0;
      daughter0 = daughter1; // the bachelor label
      daughter1 = daughterTemp; // the V0 label
    }

    if (mcPartDaughter1->GetNDaughters()!=2) {
      AliDebug(2, "The Lambda MC particle doesn't decay in 2 particles, skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    Int_t daughter1D0 = mcPartDaughter1->GetDaughter(0);
    Int_t daughter1D1 = mcPartDaughter1->GetDaughter(1);
    if (daughter1D0<0 || daughter1D1<0) {
      AliDebug(2, Form("The Lambda MC particle doesn't have correct daughters, skipping!!"));
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughter1D0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D0));
    AliAODMCParticle* mcPartDaughter1D1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D1));
    if (!mcPartDaughter1D0 || !mcPartDaughter1D1) {
      AliError("The Lambda daughter particles not found in MC array");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    if (!(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==211   &&
	  TMath::Abs(mcPartDaughter1D1->GetPdgCode())==2212) &&
	!(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==2212  &&
	  TMath::Abs(mcPartDaughter1D1->GetPdgCode())==211)) {
      AliDebug(2, "The Lambda MC doesn't decay in pi+proton (or cc), skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    // Lambda -> p+pi OR cc

    fLabelArray[0] = daughter0;//mcPartDaughter0->GetLabel(); // bachelor

    if (fmcPartCandidate->Charge()>0) {

      if (mcPartDaughter1D0->GetPdgCode()==2212) {
	fLabelArray[1] = daughter1D0;//mcPartDaughter1D0->GetLabel(); // proton
	fLabelArray[2] = daughter1D1;//mcPartDaughter1D1->GetLabel(); // pion
      } else if (mcPartDaughter1D1->GetPdgCode()==2212) {
	fLabelArray[1] = daughter1D1;//mcPartDaughter1D1->GetLabel(); // proton
	fLabelArray[2] = daughter1D0;//mcPartDaughter1D0->GetLabel(); // pion
      }

    } else if (fmcPartCandidate->Charge()<0) {

      if (mcPartDaughter1D0->GetPdgCode()==211) {
	fLabelArray[1] = daughter1D0;//mcPartDaughter1D0->GetLabel(); // pion
	fLabelArray[2] = daughter1D1;//mcPartDaughter1D1->GetLabel(); // proton
      } else if (mcPartDaughter1D1->GetPdgCode()==211) {
	fLabelArray[1] = daughter1D1;//mcPartDaughter1D1->GetLabel(); // pion
	fLabelArray[2] = daughter1D0;//mcPartDaughter1D0->GetLabel(); // proton
      }

    }

  } else if (fGenLcOption==kCountK0Sp) { // Lc -> K0bar + proton OR cc

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==311   &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==2212) &&
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==2212  &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==311)) {
      AliDebug(2, "The Lc MC doesn't decay in K0bar+proton (or cc), skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;  
    }

    if (TMath::Abs(mcPartDaughter0->GetPdgCode())==311) {
      mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); // the bachelor
      mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); // the V0
      Int_t daughterTemp = daughter0;
      daughter0 = daughter1; // the bachelor label
      daughter1 = daughterTemp; // the V0 label
    }

    Int_t daughter = mcPartDaughter1->GetDaughter(0);
    if (daughter<0) {
      AliDebug(2, Form("The K0/K0bar MC particle doesn't have correct daughter, skipping!!"));
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter));
    if (!mcPartDaughter) {
      AliError("The K0/K0bar daughter particle not found in MC array");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    if (!(TMath::Abs(mcPartDaughter->GetPdgCode())==310)) {
      AliDebug(2, "The K0/K0bar MC doesn't go in K0S, skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    if (mcPartDaughter->GetNDaughters()!=2) {
      AliDebug(2, "The K0S MC doesn't decay in 2 particles, skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    Int_t daughterD0 = mcPartDaughter->GetDaughter(0);
    Int_t daughterD1 = mcPartDaughter->GetDaughter(1);
    if (daughterD0<0 || daughterD1<0) {
      AliDebug(2, Form("The K0S MC particle doesn't have correct daughters, skipping!!"));
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    AliAODMCParticle* mcPartDaughterD0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD0));
    AliAODMCParticle* mcPartDaughterD1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD1));
    if (!mcPartDaughterD0 || !mcPartDaughterD1) {
      AliError("The K0S daughter particles not found in MC array");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    if (! ( TMath::Abs(mcPartDaughterD0->GetPdgCode())==211 &&
	    TMath::Abs(mcPartDaughterD1->GetPdgCode())==211 ) ) {
      AliDebug(2, "The K0S MC doesn't decay in pi+ pi-, skipping!!");
      delete [] fLabelArray;
      fLabelArray = 0x0;
      return checkCD;
    }

    // K0S -> pi+ pi-

    fLabelArray[0] = daughter0;//mcPartDaughter0->GetLabel(); // bachelor

    if (mcPartDaughterD0->GetPdgCode()==211) {
      fLabelArray[1] = daughterD0;//mcPartDaughterD0->GetLabel(); // pi+
      fLabelArray[2] = daughterD1;//mcPartDaughterD1->GetLabel(); // pi-
      AliDebug(2,Form(" daughter0=%d ------ daughter1=%d ------ dg0->GetLabel()=%d ------ dg1->GetLabel()=%d ",daughterD0,daughterD1,mcPartDaughterD0->GetLabel(),mcPartDaughterD1->GetLabel()));
    } else if (mcPartDaughterD1->GetPdgCode()==211) {
      fLabelArray[1] = daughterD1;//mcPartDaughterD1->GetLabel(); // pi+
      fLabelArray[2] = daughterD0;//mcPartDaughterD0->GetLabel(); // pi-
      AliDebug(2,Form(" daughter0=%d ------ daughter1=%d ------ dg0->GetLabel()=%d ------ dg1->GetLabel()=%d ",daughterD1,daughterD0,mcPartDaughterD1->GetLabel(),mcPartDaughterD0->GetLabel()));
    }
  }

  AliDebug(2,Form(" label0=%d, label1=%d, label2=%d",fLabelArray[0],fLabelArray[1],fLabelArray[2]));
  
  SetAccCut(); // setting the pt and eta acceptance cuts

  checkCD = kTRUE;
  return checkCD;

}
//____________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::FillVectorFromMCarray(AliAODMCParticle *mcPartDaughterBachelor,
							     AliAODMCParticle *mcPartDaughterK0,
							     Double_t *vectorMC)
{
  // fill the vector

  Bool_t bGenValues = kFALSE;

  AliAODMCParticle *mcPartV0DaughterPos = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[1]));
  AliAODMCParticle *mcPartV0DaughterNeg = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[2]));
  AliAODMCParticle *mcPartDaughterV0 = 0x0;

  if(!mcPartV0DaughterPos && !mcPartV0DaughterNeg) return bGenValues;

  if (TMath::Abs(mcPartDaughterK0->GetPdgCode())==311) {
    Int_t daughterK0 = mcPartDaughterK0->GetDaughter(0);
    if (daughterK0<0) {
      AliDebug(2, Form("The K0/K0bar particle doesn't have correct daughter, skipping!!"));
      return bGenValues;
    }
    mcPartDaughterV0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterK0));
    if (!mcPartDaughterV0) {
      AliDebug(2,"The K0/K0bar daughter particle not found in MC array");
      return bGenValues;
    }
    if (TMath::Abs(mcPartDaughterV0->GetPdgCode())!=310) {
      AliDebug(2,"The K0/K0bar daughter particle is not a K0S");
      return bGenValues;
    }
  } else if (TMath::Abs(mcPartDaughterK0->GetPdgCode())==3122) {
    mcPartDaughterV0 = dynamic_cast<AliAODMCParticle*>(mcPartDaughterK0);
    if (!mcPartDaughterV0) {
      AliDebug(2,"The Lambda particle not found in MC array");
      return bGenValues;
    }
  }

  if (!mcPartDaughterV0) {
    AliDebug(2,"V0 particle not found in MC array");
    return bGenValues;
  }

  Double_t cTLc = Ctau(fmcPartCandidate); // by default wrt Primary Vtx
  Double_t pTbach = mcPartDaughterBachelor->Pt(); // get the bachelor pT

  Double_t vtx1[3] = {0,0,0};   // primary vertex		
  Bool_t hasPrimVtx = fmcPartCandidate->XvYvZv(vtx1);  // cm

  // getting vertex from daughters
  Double_t vtx1daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Bool_t hasSecVtx1 = mcPartDaughterBachelor->XvYvZv(vtx1daughter0);  //cm
  Double_t vtx1daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
  Bool_t hasSecVtx2 = mcPartDaughterV0->XvYvZv(vtx1daughter1);  //cm
  if (!hasPrimVtx || !hasSecVtx1 || !hasSecVtx2) {
    AliDebug(2,"At least one of Prim.vtx, V0vtx, BachelorVtx doesn't exist!");
    //return bGenValues;
  }

  if (TMath::Abs(vtx1daughter0[0]-vtx1daughter1[0])>1E-5 ||
      TMath::Abs(vtx1daughter0[1]-vtx1daughter1[1])>1E-5 ||
      TMath::Abs(vtx1daughter0[2]-vtx1daughter1[2])>1E-5) {
    AliError("Daughters have different secondary vertex, skipping the track");
    //return bGenValues;
  }

  // getting the momentum from the daughters
  Double_t px1[2] = {mcPartDaughterBachelor->Px(), mcPartDaughterV0->Px()};
  Double_t py1[2] = {mcPartDaughterBachelor->Py(), mcPartDaughterV0->Py()};
  Double_t pz1[2] = {mcPartDaughterBachelor->Pz(), mcPartDaughterV0->Pz()};

  Int_t nprongs = 2;
  Short_t charge = mcPartDaughterBachelor->Charge();
  Double_t d0[2] = {0.,0.};
  AliAODRecoDecayHF* decayLc = new AliAODRecoDecayHF(vtx1,vtx1daughter0,nprongs,charge,px1,py1,pz1,d0);
  Double_t cosPAwrtPrimVtxLc = decayLc->CosPointingAngle();
  delete decayLc;

  // getting vertex from daughters
  Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Bool_t hasSecVtx3 = mcPartV0DaughterPos->XvYvZv(vtx2daughter0);  //cm
  Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
  Bool_t hasSecVtx4 = mcPartV0DaughterNeg->XvYvZv(vtx2daughter1);  //cm
  if (!hasSecVtx3 || !hasSecVtx4) {
    AliDebug(2,"At least one of V0Posvtx, V0Negtx doesn't exist!");
    //return bGenValues;
  }

  if (TMath::Abs(vtx2daughter0[0]-vtx2daughter1[0])>1E-5 ||
      TMath::Abs(vtx2daughter0[1]-vtx2daughter1[1])>1E-5 ||
      TMath::Abs(vtx2daughter0[2]-vtx2daughter1[2])>1E-5) {
    AliError("Daughters have different secondary vertex, skipping the track");
    //return bGenValues;
  }

  // getting the momentum from the daughters
  Double_t px[2] = {mcPartV0DaughterPos->Px(), mcPartV0DaughterNeg->Px()};		
  Double_t py[2] = {mcPartV0DaughterPos->Py(), mcPartV0DaughterNeg->Py()};		
  Double_t pz[2] = {mcPartV0DaughterPos->Pz(), mcPartV0DaughterNeg->Pz()};

  nprongs = 2;
  charge = 0;
  AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1,vtx2daughter0,nprongs,charge,px,py,pz,d0);
  Double_t cosPAwrtPrimVtxV0 = decay->CosPointingAngle();
  Double_t cTV0 = 0.; //ct
  if (fGenLcOption==kCountK0Sp) {
    cTV0 = decay->Ct(310); // by default wrt Primary Vtx
  } else if (fGenLcOption==kCountLambdapi) {
    cTV0 = decay->Ct(3122); // by default wrt Primary Vtx
  }

  Double_t invMass = 0.; //invMass
  if (fGenLcOption==kCountK0Sp) {
    invMass = decay->InvMass2Prongs(0,1,211,211);
  } else if (fGenLcOption==kCountLambdapi) {
    if (fmcPartCandidate->GetPdgCode() == 4122)
      invMass = decay->InvMass2Prongs(0,1,2212,211);
    else if (fmcPartCandidate->GetPdgCode() ==-4122)
      invMass = decay->InvMass2Prongs(0,1,211,2212);
  }
  delete decay;

  vectorMC[0]  = fmcPartCandidate->Pt();
  vectorMC[1]  = fmcPartCandidate->Y() ;
  vectorMC[2]  = fmcPartCandidate->Phi();
  vectorMC[3]  = cosPAwrtPrimVtxV0;
  vectorMC[4]  = 0; // dummy value x MC, onTheFlyStatus
  vectorMC[5]  = fCentValue; // reconstructed centrality
  vectorMC[6]  = 1; // dummy value x MC, fFake
  vectorMC[7]  = fMultiplicity; // reconstructed multiplicity

  if (fConfiguration==AliCFTaskVertexingHF::kSnail) {
    vectorMC[8]  = pTbach;
    vectorMC[9]  = mcPartV0DaughterPos->Pt();
    vectorMC[10] = mcPartV0DaughterNeg->Pt();
    vectorMC[11] = invMass;
    vectorMC[12] = 0; // dummy value x MC, V0 DCA
    vectorMC[13] = cTV0*1.E4; // in micron
    vectorMC[14] = cTLc*1.E4; // in micron
    vectorMC[15] = cosPAwrtPrimVtxLc;
  }

  bGenValues = kTRUE;
  return bGenValues;

}
