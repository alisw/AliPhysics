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

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "TClonesArray.h"
#include "AliCFVertexingHF.h"
#include "AliESDtrack.h"
#include "TDatabasePDG.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODv0.h"
#include "AliCFVertexingHFLctoV0bachelor.h"
#include "AliCFContainer.h"
#include "AliCFTaskVertexingHF.h"

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
  
  if (fRecoCandidate->GetPrimaryVtx()) AliDebug(3,"fReco Candidate has a pointer to PrimVtx\n");
  
  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;

  if ( !(lcV0bachelor->Getv0()) ) return bSignAssoc;

  Int_t nentries = fmcArray->GetEntriesFast();
  AliDebug(3,Form("nentries = %d\n", nentries));
 
  Int_t pdgCand = 4122;
  Int_t mcLabel = -1;
  Int_t mcLabelK0S = -1;
  Int_t mcLabelLambda = -1;

  // Lc->K0S+p and cc
  Int_t pdgDgLctoV0bachelor[2]={310,2212};
  Int_t pdgDgV0toDaughters[2]={211,211};
  mcLabelK0S = lcV0bachelor->MatchToMC(pdgCand,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,fmcArray,kTRUE);
  // Lc->Lambda+pi and cc
  pdgDgLctoV0bachelor[0]=3122, pdgDgLctoV0bachelor[1]=211;
  pdgDgV0toDaughters[0]=2212,  pdgDgV0toDaughters[1]=211;
  mcLabelLambda = lcV0bachelor->MatchToMC(pdgCand,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,fmcArray,kTRUE);

  if (mcLabelK0S!=-1 && mcLabelLambda!=-1)
    AliInfo("Strange: current Lc->V0+bachelor candidate has two MC different labels!");

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

  if (mcLabel==-1) return bSignAssoc;

  if (fRecoCandidate->NumberOfFakeDaughters()>0){
    fFake = 0;    // fake candidate
    if (fFakeSelection==1) return bSignAssoc;
  }
  if (fRecoCandidate->NumberOfFakeDaughters()==0){
    fFake = 2;    // non-fake candidate
    if (fFakeSelection==2) return bSignAssoc;
  }
  
  SetMCLabel(mcLabel);
  fmcPartCandidate = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fmcLabel)); 

  if (!fmcPartCandidate){
    AliDebug(3,"No part candidate");
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

  if (fmcPartCandidate->GetNDaughters()!=2) return bGenValues;

  Int_t daughter0lc = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1lc = fmcPartCandidate->GetDaughter(1);

  //the V0
  AliAODMCParticle* mcPartDaughterV0=0;
  if (fGenLcOption==kCountLambdapi) {
    mcPartDaughterV0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0lc)); // strange baryon (0)
    if(mcPartDaughterV0){
      AliInfo(Form(" Case Lc->Lambda+pi: (V0) daughter0_Lc=%d (%d)",daughter0lc,mcPartDaughterV0->GetPdgCode()));
    }
  }
  else if (fGenLcOption==kCountK0Sp) {
    AliAODMCParticle* mcPartDaughterK0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1lc)); // strange meson (1)
    if (!mcPartDaughterK0) return bGenValues;
    if (TMath::Abs(mcPartDaughterK0->GetPdgCode())!=311) return bGenValues;
    Int_t daughterK0 = mcPartDaughterK0->GetDaughter(0);
    mcPartDaughterV0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterK0));
    if (!mcPartDaughterV0) return bGenValues;
    if (TMath::Abs(mcPartDaughterV0->GetPdgCode())!=310) return bGenValues;
    AliInfo(Form(" Case Lc->K0S+p: (V0) daughter1_Lc=%d (%d to %d)",daughter1lc,mcPartDaughterK0->GetPdgCode(),
		 mcPartDaughterV0->GetPdgCode()));
  }

  //the bachelor
  AliAODMCParticle* mcPartDaughterBachelor=0;
  if (fGenLcOption==kCountLambdapi) {
    mcPartDaughterBachelor = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1lc)); // meson (1)
    if(mcPartDaughterBachelor){
      AliInfo(Form(" Case Lc->Lambda+pi: (bachelor) daughter1_Lc=%d (%d)",daughter1lc,mcPartDaughterBachelor->GetPdgCode()));
    }
  }
  if (fGenLcOption==kCountK0Sp) {
    mcPartDaughterBachelor = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0lc)); // baryon (0)
    if(mcPartDaughterBachelor){
      AliInfo(Form(" Case Lc->K0S+p: (bachelor) daughter0_Lc=%d (%d)",daughter0lc,mcPartDaughterBachelor->GetPdgCode()));
    }
  }

  if (!mcPartDaughterV0 || !mcPartDaughterBachelor) return bGenValues;

  if (fGenLcOption==kCountLambdapi) {
    if ( !(mcPartDaughterV0->GetPdgCode()==3122 &&
	   mcPartDaughterBachelor->GetPdgCode()==211) )
      AliInfo("It isn't a Lc->Lambda+pion candidate");
  }
  if (fGenLcOption==kCountK0Sp) {
    if ( !(mcPartDaughterV0->GetPdgCode()==310 &&
	   mcPartDaughterBachelor->GetPdgCode()==2212) )
      AliInfo("It isn't a Lc->K0S+proton candidate");
  }

  Double_t vtx1[3] = {0,0,0};   // primary vertex		
  fmcPartCandidate->XvYvZv(vtx1);  // cm

  // getting vertex from daughters
  Double_t vtx1daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Double_t vtx1daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
  mcPartDaughterBachelor->XvYvZv(vtx1daughter0);  //cm
  mcPartDaughterV0->XvYvZv(vtx1daughter1);  //cm
  if (TMath::Abs(vtx1daughter0[0]-vtx1daughter1[0])>1E-5 ||
      TMath::Abs(vtx1daughter0[1]-vtx1daughter1[1])>1E-5 ||
      TMath::Abs(vtx1daughter0[2]-vtx1daughter1[2])>1E-5) {
    AliError("Daughters have different secondary vertex, skipping the track");
    return bGenValues;
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

  //V0 daughters
  Int_t daughter0 = mcPartDaughterV0->GetDaughter(0);
  Int_t daughter1 = mcPartDaughterV0->GetDaughter(1);
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); //V0daughter0
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); //V0daughter1

  if (!mcPartDaughter0 || !mcPartDaughter1) return kFALSE;

  // getting vertex from daughters
  Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
  mcPartDaughter0->XvYvZv(vtx2daughter0);  //cm
  mcPartDaughter1->XvYvZv(vtx2daughter1);  //cm
  if (TMath::Abs(vtx2daughter0[0]-vtx2daughter1[0])>1E-5 ||
      TMath::Abs(vtx2daughter0[1]-vtx2daughter1[1])>1E-5 ||
      TMath::Abs(vtx2daughter0[2]-vtx2daughter1[2])>1E-5) {
    AliError("Daughters have different secondary vertex, skipping the track");
    return bGenValues;
  }
	
  // always instantiate the AliAODRecoDecay with the positive daughter first, the negative second
  AliAODMCParticle* positiveDaugh = mcPartDaughter0;
  AliAODMCParticle* negativeDaugh = mcPartDaughter1;
  if (mcPartDaughter0->GetPdgCode()<0 && mcPartDaughter1->GetPdgCode()>0){
    // inverting in case the positive daughter is the second one
    positiveDaugh = mcPartDaughter1;
    negativeDaugh = mcPartDaughter0;
  }
  // getting the momentum from the daughters
  Double_t px[2] = {positiveDaugh->Px(), negativeDaugh->Px()};		
  Double_t py[2] = {positiveDaugh->Py(), negativeDaugh->Py()};		
  Double_t pz[2] = {positiveDaugh->Pz(), negativeDaugh->Pz()};

  nprongs = 2;
  charge = 0;
  AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1,vtx2daughter0,nprongs,charge,px,py,pz,d0);
  Double_t cosPAwrtPrimVtxV0 = decay->CosPointingAngle();

  //ct
  Double_t cTV0 = 0.;
  if (fGenLcOption==kCountK0Sp) {
    cTV0 = decay->Ct(310); // by default wrt Primary Vtx
  } else if (fGenLcOption==kCountLambdapi) {
    cTV0 = decay->Ct(3122); // by default wrt Primary Vtx
  }

  Double_t cTLc = Ctau(fmcPartCandidate); // by default wrt Primary Vtx

  // get the bachelor pT
  Double_t pTbach = 0.;

  if (TMath::Abs(fmcPartCandidate->GetPdgCode()) == 4122)
    pTbach = mcPartDaughterBachelor->Pt();

  Double_t invMass = 0.;
  if (fGenLcOption==kCountK0Sp) {
    invMass = decay->InvMass2Prongs(0,1,211,211);
  } else if (fGenLcOption==kCountLambdapi) {
    if (fmcPartCandidate->GetPdgCode() == 4122)
      invMass = decay->InvMass2Prongs(0,1,2212,211);
    else if (fmcPartCandidate->GetPdgCode() ==-4122)
      invMass = decay->InvMass2Prongs(0,1,211,2212);
  }

  switch (fConfiguration){
  case AliCFTaskVertexingHF::kSnail:
    vectorMC[0]  = fmcPartCandidate->Pt();
    vectorMC[1]  = fmcPartCandidate->Y() ;
    vectorMC[2]  = fmcPartCandidate->Phi();
    vectorMC[3]  = cosPAwrtPrimVtxV0;
    vectorMC[4]  = 0; // dummy value x MC, onTheFlyStatus
    vectorMC[5]  = fCentValue; // reconstructed centrality
    vectorMC[6]  = 1; // dummy value x MC, fFake
    vectorMC[7]  = fMultiplicity; // reconstructed multiplicity

    vectorMC[8]  = pTbach;
    vectorMC[9]  = positiveDaugh->Pt();
    vectorMC[10] = negativeDaugh->Pt();
    vectorMC[11] = invMass;
    vectorMC[12] = 0; // dummy value x MC, V0 DCA
    vectorMC[13] = cTV0*1.E4; // in micron
    vectorMC[14] = cTLc*1.E4; // in micron
    vectorMC[15] = cosPAwrtPrimVtxLc;
    break;
  case AliCFTaskVertexingHF::kCheetah:
    vectorMC[0]  = fmcPartCandidate->Pt();
    vectorMC[1]  = fmcPartCandidate->Y() ;
    vectorMC[2]  = fmcPartCandidate->Phi();
    vectorMC[3]  = cosPAwrtPrimVtxV0;
    vectorMC[4]  = 0; // dummy value x MC, onTheFlyStatus
    vectorMC[5]  = fCentValue; // reconstructed centrality
    vectorMC[6]  = 1; // dummy value x MC, fFake
    vectorMC[7]  = fMultiplicity; // reconstructed multiplicity
    break;
  }

  delete decay;
  delete decayLc;

  bGenValues = kTRUE;
  return bGenValues;

}
//____________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{ 
  // read the variables for the container

  Bool_t bFillRecoValues = kFALSE;

  //Get the Lc and the V0 from Lc
  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;

  if ( !(lcV0bachelor->Getv0()) ) return bFillRecoValues;

  AliAODVertex *vtx0 = (AliAODVertex*)lcV0bachelor->GetPrimaryVtx();
  if (vtx0) AliDebug(1,"lcV0bachelor has primary vtx");

  AliAODTrack* bachelor = (AliAODTrack*)lcV0bachelor->GetBachelor();
  AliAODv0* v0toDaughters = (AliAODv0*)lcV0bachelor->Getv0();
  Bool_t onTheFlyStatus = v0toDaughters->GetOnFlyStatus();
  AliAODTrack* v0positiveTrack = (AliAODTrack*)lcV0bachelor->Getv0PositiveTrack();
  AliAODTrack* v0negativeTrack = (AliAODTrack*)lcV0bachelor->Getv0NegativeTrack();

  Double_t pt = lcV0bachelor->Pt();
  Double_t rapidity = lcV0bachelor->Y(4122);
  Double_t invMassV0 = 0.;
  Double_t primVtxPos[3]; vtx0->GetXYZ(primVtxPos);
  Double_t cosPAwrtPrimVtxV0 = v0toDaughters->CosPointingAngle(primVtxPos);
  Double_t cTV0 = 0.;

  Double_t pTbachelor = bachelor->Pt();
  Double_t pTV0pos = v0positiveTrack->Pt();
  Double_t pTV0neg = v0negativeTrack->Pt();
  Double_t phi = lcV0bachelor->Phi();
  Double_t dcaV0 = v0toDaughters->GetDCA();
  Double_t cTLc = lcV0bachelor->Ct(4122); // wrt PrimVtx
  //Double_t dcaLc = lcV0bachelor->GetDCA();
  Double_t cosPointingAngleLc = lcV0bachelor->CosPointingAngle();

  if (fGenLcOption==kCountK0Sp) {
    cTV0 = v0toDaughters->Ct(310,primVtxPos);
  } else if (fGenLcOption==kCountLambdapi) {
    cTV0 = v0toDaughters->Ct(3122,primVtxPos);
  }

  Short_t bachelorCharge = bachelor->Charge();
  if (fGenLcOption==kCountLambdapi) {
    if (bachelorCharge==1) {
      invMassV0 = v0toDaughters->MassLambda();
    } else if (bachelorCharge==-1) {
      invMassV0 = v0toDaughters->MassAntiLambda();
    }

  } else if (fGenLcOption==kCountK0Sp) {
    invMassV0 = v0toDaughters->MassK0Short();
  }

  switch (fConfiguration){
  case AliCFTaskVertexingHF::kSnail:
    vectorReco[0]  = pt;
    vectorReco[1]  = rapidity;
    vectorReco[2]  = phi;
    vectorReco[3]  = cosPAwrtPrimVtxV0;
    vectorReco[4]  = onTheFlyStatus;
    vectorReco[5]  = fCentValue;
    vectorReco[6]  = fFake; // whether the reconstructed candidate was a fake (fFake = 0) or not (fFake = 2) 
    vectorReco[7]  = fMultiplicity;

    vectorReco[8]  = pTbachelor;
    vectorReco[9]  = pTV0pos;
    vectorReco[10] = pTV0neg;
    vectorReco[11] = invMassV0;
    vectorReco[12] = dcaV0;
    vectorReco[13] = cTV0*1.E4; // in micron
    vectorReco[14] = cTLc*1.E4; // in micron
    vectorReco[15] = cosPointingAngleLc;

    break;
  case AliCFTaskVertexingHF::kCheetah:
    vectorReco[0]  = pt;
    vectorReco[1]  = rapidity;
    vectorReco[2]  = phi;
    vectorReco[3]  = cosPAwrtPrimVtxV0;
    vectorReco[4]  = onTheFlyStatus;
    vectorReco[5]  = fCentValue;
    vectorReco[6]  = fFake; // whether the reconstructed candidate was a fake (fFake = 0) or not (fFake = 2) 
    vectorReco[7]  = fMultiplicity;
    break;
  }

  bFillRecoValues = kTRUE;

  return bFillRecoValues;
}

//_____________________________________________________________
Bool_t AliCFVertexingHFLctoV0bachelor::CheckMCChannelDecay() const
{ 
  // check the required decay channel

  Bool_t checkCD = kFALSE;
  
  if (fmcPartCandidate->GetNDaughters()!=2) return checkCD;

  Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));

  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug (2,"Problems in the MC Daughters\n");
    return checkCD;
  }

  // Lc -> Lambda + pion AND cc
  if (fGenLcOption==kCountLambdapi) {

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==3122 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==211 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==3122)) {
      AliDebug(2, "The Lc MC doesn't decay in Lambda+pion (or cc), skipping!!");
      return checkCD;  
    } else if (TMath::Abs(mcPartDaughter0->GetPdgCode())==3122) {
      if (mcPartDaughter0->GetNDaughters()!=2) return checkCD;
      Int_t daughter0D0 = mcPartDaughter0->GetDaughter(0);
      AliAODMCParticle* mcPartDaughter0D0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0D0));
      if(!mcPartDaughter0D0){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      Int_t daughter0D1 = mcPartDaughter0->GetDaughter(1);
      AliAODMCParticle* mcPartDaughter0D1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0D1));
      if(!mcPartDaughter0D1){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      if (!(TMath::Abs(mcPartDaughter0D0->GetPdgCode())==211 &&
	    TMath::Abs(mcPartDaughter0D1->GetPdgCode())==2212) &&
	  !(TMath::Abs(mcPartDaughter0D0->GetPdgCode())==2212 &&
	    TMath::Abs(mcPartDaughter0D1->GetPdgCode())==211)) {
	AliDebug(2, "The Lambda MC doesn't decay in pi+proton (or cc), skipping!!");
	return checkCD;
      }
    } else if (TMath::Abs(mcPartDaughter1->GetPdgCode())==3122) {
      if (mcPartDaughter1->GetNDaughters()!=2) return checkCD;
      Int_t daughter1D0 = mcPartDaughter1->GetDaughter(0);
      AliAODMCParticle* mcPartDaughter1D0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D0));
      if(!mcPartDaughter1D0){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      Int_t daughter1D1 = mcPartDaughter1->GetDaughter(1);
      AliAODMCParticle* mcPartDaughter1D1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1D1));
      if(!mcPartDaughter1D1){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      if (!(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==211 &&
	    TMath::Abs(mcPartDaughter1D1->GetPdgCode())==2212) &&
	  !(TMath::Abs(mcPartDaughter1D0->GetPdgCode())==2212 &&
	    TMath::Abs(mcPartDaughter1D1->GetPdgCode())==211)) {
	AliDebug(2, "The Lambda MC doesn't decay in pi+proton (or cc), skipping!!");
	return checkCD;
      }
    }

  } else if (fGenLcOption==kCountK0Sp) {
  // Lc -> K0bar + proton AND cc

    if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==311 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==2212) &&
	!(TMath::Abs(mcPartDaughter0->GetPdgCode())==2212 &&
	  TMath::Abs(mcPartDaughter1->GetPdgCode())==311)) {
      AliDebug(2, "The Lc MC doesn't decay in K0+proton (or cc), skipping!!");
      return checkCD;  
    }
    
    if (TMath::Abs(mcPartDaughter0->GetPdgCode())==311) {
      Int_t daughter = mcPartDaughter0->GetDaughter(0);
      AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter));
      if(!mcPartDaughter){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }

      if (!(TMath::Abs(mcPartDaughter->GetPdgCode())==310)) {
	AliDebug(2, "The K0 (or K0bar) MC doesn't go in K0S, skipping!!");
	return checkCD;
      }
      if (mcPartDaughter->GetNDaughters()!=2) return checkCD;
      Int_t daughterD0 = mcPartDaughter->GetDaughter(0);
      AliAODMCParticle* mcPartDaughterD0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD0));
      if(!mcPartDaughterD0){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      Int_t daughterD1 = mcPartDaughter->GetDaughter(1);
      AliAODMCParticle* mcPartDaughterD1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD1));
      if(!mcPartDaughterD1){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      if (!(TMath::Abs(mcPartDaughterD0->GetPdgCode())==211 &&
	    TMath::Abs(mcPartDaughterD1->GetPdgCode())==211)) {
	AliDebug(2, "The K0S MC doesn't decay in pi+pi, skipping!!");
	return checkCD;
      }

    }

    if (TMath::Abs(mcPartDaughter1->GetPdgCode())==311) {
      Int_t daughter = mcPartDaughter1->GetDaughter(0);
      AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter));
      if (!(TMath::Abs(mcPartDaughter->GetPdgCode())==310)) {
	AliDebug(2, "The K0 (or K0bar) MC doesn't go in K0S, skipping!!");
	return checkCD;
      }
      if (mcPartDaughter->GetNDaughters()!=2) return checkCD;
      Int_t daughterD0 = mcPartDaughter->GetDaughter(0);
      AliAODMCParticle* mcPartDaughterD0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD0));
      Int_t daughterD1 = mcPartDaughter->GetDaughter(1);
      AliAODMCParticle* mcPartDaughterD1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD1));
      if (! ( TMath::Abs(mcPartDaughterD0->GetPdgCode())==211 &&
	      TMath::Abs(mcPartDaughterD1->GetPdgCode())==211 ) ) {
	AliDebug(2, "The K0S MC doesn't decay in pi+pi, skipping!!");
	return checkCD;
      }

    }

  }
  
  checkCD = kTRUE;
  return checkCD;
  
}

//___________________________________________________________

void AliCFVertexingHFLctoV0bachelor::SetPtAccCut(Float_t* ptAccCut)
{
  //
  // setting the pt cut to be used in the Acceptance steps (MC+Reco)
  //

  AliInfo("The 1st element of the pt cut array will correspond to the cut applied to the bachelor - please check that it is correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fPtAccCut[iP]=ptAccCut[iP];
    }
  }
  return;
}

//___________________________________________________________

void AliCFVertexingHFLctoV0bachelor::SetEtaAccCut(Float_t* etaAccCut)
{
  //
  // setting the eta cut to be used in the Acceptance steps (MC+Reco)
  //

  AliInfo("The 1st element of the pt cut array will correspond to the cut applied to the bachelor - please check that it is correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fEtaAccCut[iP]=etaAccCut[iP];
    }
  }
  return;
}

//___________________________________________________________

void AliCFVertexingHFLctoV0bachelor::SetAccCut(Float_t* ptAccCut, Float_t* etaAccCut)
{
  //
  // setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
  //

  AliInfo("The 1st element of the pt cut array will correspond to the cut applied to the bachelor - please check that it is correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fPtAccCut[iP]=ptAccCut[iP];
      fEtaAccCut[iP]=etaAccCut[iP];
    }
  }
  return;
}

//___________________________________________________________

void AliCFVertexingHFLctoV0bachelor::SetAccCut()
{
  //
  // setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
  //

  AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[0])); // bachelor
  if (!mcPartDaughter) return;
  Int_t mother =  mcPartDaughter->GetMother();
  AliAODMCParticle* mcMother = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother)); 
  if (!mcMother) return;
  /*
  if (fGenLcOption==kCountK0Sp) {
    if ( !(TMath::Abs(mcPartDaughter->GetPdgCode())== 2212) )
      AliFatal(Form("Apparently the proton bachelor is not in the first position (%d <- %d), causing a crash!!",
		    mcPartDaughter->GetPdgCode(),mcMother->GetPdgCode()));
  }
  else if (fGenLcOption==kCountLambdapi) {
    if ( !(TMath::Abs(mcPartDaughter->GetPdgCode())== 211) )
      AliFatal(Form("Apparently the pion bachelor is not in the first position (%d <- %d), causing a crash!!",
		    mcPartDaughter->GetPdgCode(),mcMother->GetPdgCode()));
  }
  */
  if (fProngs>0) {
    fPtAccCut[0]=0.3;  // bachelor
    fEtaAccCut[0]=0.9;  // bachelor
    for (Int_t iP=1; iP<fProngs; iP++){
      fPtAccCut[iP]=0.1;
      fEtaAccCut[iP]=0.9;
    }
  }

  return;

}		

//_____________________________________________________________
Double_t AliCFVertexingHFLctoV0bachelor::GetEtaProng(Int_t iProng) const 
{
  //
  // getting eta of the prong - overload the mother class method
  //

  if (fRecoCandidate) {
   
    AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;

    Double_t etaProng =-9999;
    if (iProng==0) etaProng = lcV0bachelor->GetBachelor()->Eta();
    else if (iProng==1) etaProng = lcV0bachelor->Getv0PositiveTrack()->Eta();
    else if (iProng==2) etaProng = lcV0bachelor->Getv0NegativeTrack()->Eta();

    return etaProng;
    
  }
  return 999999;    
}

//_____________________________________________________________

Double_t AliCFVertexingHFLctoV0bachelor::GetPtProng(Int_t iProng) const 
{
  //
  // getting pt of the prong
  //

  Double_t ptProng=-9999.;

  if (!fRecoCandidate) return ptProng;

  AliAODRecoCascadeHF* lcV0bachelor = (AliAODRecoCascadeHF*)fRecoCandidate;

  if (iProng==0) ptProng = lcV0bachelor->GetBachelor()->Pt();
  else if (iProng==1) ptProng = lcV0bachelor->Getv0PositiveTrack()->Pt();
  else if (iProng==2) ptProng = lcV0bachelor->Getv0NegativeTrack()->Pt();
    
  return ptProng;
  
}

//_____________________________________________________________

Double_t AliCFVertexingHFLctoV0bachelor::Ctau(AliAODMCParticle *mcPartCandidate)
{

  Double_t cTau = 999999.;

  Double_t vtx1[3] = {0,0,0};   // primary vertex		
  Bool_t hasProdVertex = mcPartCandidate->XvYvZv(vtx1);  // cm

  AliAODMCParticle *mcPartDaughter0 = (AliAODMCParticle*)fmcArray->At(mcPartCandidate->GetDaughter(0));
  AliAODMCParticle *mcPartDaughter1 = (AliAODMCParticle*)fmcArray->At(mcPartCandidate->GetDaughter(1));
  if (!mcPartDaughter0 && !mcPartDaughter1) return cTau;
  Double_t vtx1daughter[3] = {0,0,0};   // secondary vertex
  if (mcPartDaughter0)
    hasProdVertex = hasProdVertex || mcPartDaughter0->XvYvZv(vtx1daughter);  //cm
  if (mcPartDaughter1)
    hasProdVertex = hasProdVertex || mcPartDaughter1->XvYvZv(vtx1daughter);  //cm

  if (!hasProdVertex) return cTau;

  Double_t decayLength = 0.;
  for (Int_t ii=0; ii<3; ii++) decayLength += (vtx1daughter[ii]-vtx1[ii])*(vtx1daughter[ii]-vtx1[ii]);
  decayLength = TMath::Sqrt(decayLength);

  cTau = decayLength * mcPartCandidate->M()/mcPartCandidate->P();
  AliDebug(2,Form(" cTau(4122)=%f",cTau));
  return cTau;

}
