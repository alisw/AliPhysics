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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and steps
// For D* and other cascades
//
// Author : A.Grelli a.grelli@uu.nl  UTECHT
//-----------------------------------------------------------------------

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "TClonesArray.h"
#include "AliCFVertexingHF.h"
#include "AliESDtrack.h"
#include "TDatabasePDG.h"
#include "AliAODRecoCascadeHF.h"
#include "AliCFVertexingHFCascade.h"
#include "AliCFContainer.h"
#include "AliCFTaskVertexingHF.h"
#include "AliPIDResponse.h"
#include "AliPID.h"

/// \cond CLASSIMP
ClassImp(AliCFVertexingHFCascade);
/// \endcond

//_________________________________________
AliCFVertexingHFCascade::AliCFVertexingHFCascade():
AliCFVertexingHF(),
  fPDGcascade(0),
  fPDGbachelor(0),
  fPDGneutrDaugh(0),
  fPDGneutrDaughForMC(0),
  fPDGneutrDaughPositive(0),
  fPDGneutrDaughNegative(0),
  fPrimVtx(0x0),
  fUseCutsForTMVA(kFALSE),
  fCutOnMomConservation(0.00001)
{
  // default constructor

  SetNProngs(3);
  fPtAccCut = new Float_t[fProngs];
  fEtaAccCut = new Float_t[fProngs];
  // element 0 in the cut arrays corresponds to the soft pion!!!!!!!! Careful when setting the values...
  fPtAccCut[0] = 0.;
  fEtaAccCut[0] = 0.;
  for(Int_t iP=1; iP<fProngs; iP++){
    fPtAccCut[iP] = 0.1;
    fEtaAccCut[iP] = 0.9;
  }

}

//_________________________________________
AliCFVertexingHFCascade::AliCFVertexingHFCascade(TClonesArray *mcArray, UShort_t originDselection):
AliCFVertexingHF(mcArray, originDselection),
  fPDGcascade(0),
  fPDGbachelor(0),
  fPDGneutrDaugh(0),
  fPDGneutrDaughForMC(0),
  fPDGneutrDaughPositive(0),
  fPDGneutrDaughNegative(0),
  fPrimVtx(0x0),
  fUseCutsForTMVA(kFALSE),
  fCutOnMomConservation(0.00001)
{
  // standard constructor

  SetNProngs(3);
  fPtAccCut = new Float_t[fProngs];
  fEtaAccCut = new Float_t[fProngs];
  // element 0 in the cut arrays corresponds to the soft pion!!!!!!!! Careful when setting the values...
  fPtAccCut[0] = 0.;
  fEtaAccCut[0] = 0.;
  for(Int_t iP=1; iP<fProngs; iP++){
    fPtAccCut[iP] = 0.1;
    fEtaAccCut[iP] = 0.9;
  }

}


//_____________________________________
AliCFVertexingHFCascade& AliCFVertexingHFCascade::operator=(const AliCFVertexingHFCascade& c)
{
  // operator =
 
  if  (this != &c) {

    AliCFVertexingHF::operator=(c);

  }
  return *this;
}

//__________________________________________
Bool_t AliCFVertexingHFCascade::SetRecoCandidateParam(AliAODRecoDecayHF *recoCand)
{
  // set the AliAODRecoDecay candidate

  Bool_t bSignAssoc = kFALSE;

  fRecoCandidate = recoCand;
  AliAODRecoCascadeHF* cascade = (AliAODRecoCascadeHF*)recoCand;

  if (!fRecoCandidate) {
    AliError("fRecoCandidate not found, problem in assignement\n");
    return bSignAssoc;
  }

  if ( fRecoCandidate->GetPrimaryVtx()) AliDebug(3,"fReco Candidate has a pointer to PrimVtx\n");

  //Int_t pdgCand = 413;

  Int_t pdgDgCascade[2] = {fPDGneutrDaugh, fPDGbachelor};
  Int_t pdgDgNeutrDaugh[2] = {fPDGneutrDaughPositive, fPDGneutrDaughNegative};

  Int_t nentries = fmcArray->GetEntriesFast();

  AliDebug(3,Form("nentries = %d\n", nentries));
 
  Bool_t isV0 = kFALSE;
  if (fPDGcascade == 4122) {
    isV0 = kTRUE;
    pdgDgCascade[0] = fPDGbachelor;
    pdgDgCascade[1] = fPDGneutrDaugh;
  }
  AliDebug(3, Form("calling MatchToMC with: fPDGcascade = %d, fPDGneutrDaugh = %d, pdgDgCascade[0] = %d, pdgDgCascade[1] = %d, pdgDgNeutrDaugh[0] = %d, pdgDgNeutrDaugh[1] = %d, fmcArray = %p", fPDGcascade, fPDGneutrDaugh, pdgDgCascade[0], pdgDgCascade[1], pdgDgNeutrDaugh[0], pdgDgNeutrDaugh[1], fmcArray));
 Int_t mcLabel = cascade->MatchToMC(fPDGcascade, fPDGneutrDaugh, pdgDgCascade, pdgDgNeutrDaugh, fmcArray, isV0); 
  
  if (mcLabel == -1) return bSignAssoc;

  if (fRecoCandidate->NumberOfFakeDaughters()>0){
    fFake = 0;    // fake candidate
    if (fFakeSelection == 1) return bSignAssoc;
  }
  if (fRecoCandidate->NumberOfFakeDaughters()==0){
    fFake = 2;    // non-fake candidate
    if (fFakeSelection == 2) return bSignAssoc;
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
Bool_t AliCFVertexingHFCascade::GetGeneratedValuesFromMCParticle(Double_t* vectorMC) 
{
  // 
  // collecting all the necessary info (pt, y, cosThetaStar, ptPi, ptKa, cT) from MC particle
  //

  Bool_t bGenValues = kFALSE;

  Int_t daughter0cascade = fmcPartCandidate->GetDaughterLabel(0);
  Int_t daughter1cascade = fmcPartCandidate->GetDaughterLabel(1);

  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0cascade));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1cascade));
  AliAODMCParticle* mcPartDaughterNeutrDaugh = NULL;
  AliAODMCParticle* mcPartDaughterBachelor = NULL;

  // the Neutral Particle (e.g. D0 for D*, K0S for Lc...)
  // for D* the D0 (the neutral) is the first daughter, while for Lc the V0 is the second, so we check the 
  // charge of the daughters to decide which is which
  if (mcPartDaughter0->Charge()/3 == 0){
    mcPartDaughterNeutrDaugh = mcPartDaughter0;
    mcPartDaughterBachelor =  mcPartDaughter1;
  }
  else {
    mcPartDaughterNeutrDaugh = mcPartDaughter1;
    mcPartDaughterBachelor =  mcPartDaughter0;
  }

  if (!mcPartDaughterNeutrDaugh || !mcPartDaughterBachelor) return kFALSE;

  Double_t vtx1[3] = {0,0,0};   // primary vertex
  Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
  Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
  fmcPartCandidate->XvYvZv(vtx1);  // cm

  //Daughters of the neutral particle of the cascade
  Int_t daughter0 = mcPartDaughterNeutrDaugh->GetDaughterLabel(0); // this is the positive
  Int_t daughter1 = mcPartDaughterNeutrDaugh->GetDaughterLabel(1); // this is the negative

  AliAODMCParticle* mcPartNeutrDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartNeutrDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));

  if (!mcPartNeutrDaughter0 || !mcPartNeutrDaughter1) return kFALSE;

  // getting vertex from daughters
  mcPartNeutrDaughter0->XvYvZv(vtx2daughter0);  // cm
  mcPartNeutrDaughter1->XvYvZv(vtx2daughter1);  //cm
  if (TMath::Abs(vtx2daughter0[0] - vtx2daughter1[0])>1E-5 || TMath::Abs(vtx2daughter0[1]- vtx2daughter1[1])>1E-5 || TMath::Abs(vtx2daughter0[2] - vtx2daughter1[2])>1E-5) {
    AliError("Daughters have different secondary vertex, skipping the track");
    return bGenValues;
  }

  Int_t nprongs = 2;
  Short_t charge = 0;
  // always instantiate the AliAODRecoDecay with the positive daughter first, the negative second
  AliAODMCParticle* positiveDaugh = mcPartNeutrDaughter0;
  AliAODMCParticle* negativeDaugh = mcPartNeutrDaughter1;
  if (mcPartNeutrDaughter0->GetPdgCode() < 0 && mcPartNeutrDaughter1->GetPdgCode() > 0){
    // inverting in case the positive daughter is the second one
    positiveDaugh = mcPartNeutrDaughter1;
    negativeDaugh = mcPartNeutrDaughter0;
  }

  // getting the momentum from the daughters
  Double_t px[2] = {positiveDaugh->Px(), negativeDaugh->Px()};		
  Double_t py[2] = {positiveDaugh->Py(), negativeDaugh->Py()};		
  Double_t pz[2] = {positiveDaugh->Pz(), negativeDaugh->Pz()};

  Double_t d0[2] = {0.,0.};

  AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1, vtx2daughter0, nprongs, charge, px, py, pz, d0);
	
  Double_t cosThetaStar = 0.;
  Double_t cosThetaStarNeutrDaugh = 0.;
  Double_t cosThetaStarNeutrDaughBar = 0.;
  cosThetaStarNeutrDaugh = decay->CosThetaStar(1, fPDGneutrDaugh, fPDGneutrDaughPositive, fPDGneutrDaughNegative);
  cosThetaStarNeutrDaughBar = decay->CosThetaStar(0, fPDGneutrDaugh, fPDGneutrDaughNegative, fPDGneutrDaughPositive);
  if (mcPartDaughterNeutrDaugh->GetPdgCode() == fPDGneutrDaughForMC){  // neutral particle
    AliDebug(3, Form("Neutral Daughter, with pdgprong0 = %d, pdgprong1 = %d", mcPartDaughter0->GetPdgCode(), mcPartDaughter1->GetPdgCode()));
    cosThetaStar = cosThetaStarNeutrDaugh;
  }
  else if (mcPartDaughterNeutrDaugh->GetPdgCode() == -fPDGneutrDaughForMC){  // neutral particle bar
    AliDebug(3, Form("Neutral Daughter, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
    cosThetaStar = cosThetaStarNeutrDaughBar;
  }
  else{
    AliWarning(Form("There are problems!! particle was expected to be either with pdg = %d or its antiparticle with pdg = %d, while we have a %d, check...", fPDGneutrDaughForMC, -fPDGneutrDaughForMC, mcPartDaughterNeutrDaugh->GetPdgCode()));
    delete decay;
    return vectorMC;
  }
  if (cosThetaStar < -1 || cosThetaStar > 1) {
    AliWarning(Form("Invalid value for cosine Theta star %f, returning", cosThetaStar));
    delete decay;
    return bGenValues;
  }
	
  Double_t vectorNeutrDaugh[2] = {0.,0.};

  // evaluate the correct cascade
  if (!EvaluateIfCorrectNeutrDaugh(mcPartDaughterNeutrDaugh, vectorNeutrDaugh)) {
    AliDebug(2, "Error! the Neutral Daughter MC doesn't have correct daughters!!");
    delete decay;
    return bGenValues;  
  }	

  //ct
  Double_t cT = decay->Ct(fPDGneutrDaugh);
  // // get the pT of the daughters
  // Double_t pTNeutrDaugh= 0.;
  // Double_t pTBachelor = 0.;
	
  // if (TMath::Abs(fmcPartCandidate->GetPdgCode()) == fPDGcascade) {
  //   pTNeutrDaugh = mcPartDaughterNeutrDaugh->Pt();
  //   pTBachelor = mcPartDaughterBachelor->Pt();
  // }

  AliDebug(3, Form("The candidate has pt = %f, y = %f", fmcPartCandidate->Pt(), fmcPartCandidate->Y()));
  
  Int_t localmult = -1;
  if(fConfiguration == AliCFTaskVertexingHF::kESE) {
    localmult = ComputeLocalMultiplicity(decay->Eta(), decay->Phi(), 0.4);
  }

  switch (fConfiguration){
  case AliCFTaskVertexingHF::kSnail:
    vectorMC[0] = fmcPartCandidate->Pt();
    vectorMC[1] = fmcPartCandidate->Y() ;
    vectorMC[2] = cosThetaStar ;
    vectorMC[3] = vectorNeutrDaugh[0];
    vectorMC[4] = vectorNeutrDaugh[1];
    vectorMC[5] = cT*1.E4 ;  // in micron
    vectorMC[6] = 0.;   // dummy value, meaningless in MC
    vectorMC[7] = -100000.; // dummy value, meaningless in MC, in micron^2
    vectorMC[8] = 1.01;    // dummy value, meaningless in MC
    vectorMC[9] = fmcPartCandidate->Phi(); 
    vectorMC[10] = fzMCVertex;    // z of reconstructed of primary vertex
    vectorMC[11] = fCentValue; // reconstructed centrality
    vectorMC[12] = 1.;           // always filling with 1 at MC level 
    vectorMC[13] = 1.01; // dummy value for cosPointingXY  multiplicity
    vectorMC[14] = 0.; // dummy value for NormalizedDecayLengthXY multiplicity
    vectorMC[15] = fMultiplicity; // reconstructed multiplicity
    break;
  case AliCFTaskVertexingHF::kCheetah:
    vectorMC[0] = fmcPartCandidate->Pt();
    vectorMC[1] = fmcPartCandidate->Y() ;
    vectorMC[2] = cT*1.E4; // in micron
    vectorMC[3] = fmcPartCandidate->Phi();
    vectorMC[4] = fzMCVertex;
    vectorMC[5] = fCentValue;   // dummy value for dca, meaningless in MC
    vectorMC[6] = 1. ;  // fake: always filling with 1 at MC level 
    vectorMC[7] = fMultiplicity;   // dummy value for d0pi, meaningless in MC, in micron
    break;
  case AliCFTaskVertexingHF::kFalcon:
    vectorMC[0] = fmcPartCandidate->Pt();
    vectorMC[1] = fmcPartCandidate->Y() ;
    vectorMC[2] = fCentValue;   // dummy value for dca, meaningless in MC
    vectorMC[3] = fMultiplicity;   // dummy value for d0pi, meaningless in MC, in micron
    break;
  case AliCFTaskVertexingHF::kESE:
    vectorMC[0] = fmcPartCandidate->Pt();
    vectorMC[1] = fmcPartCandidate->Y() ;
    vectorMC[2] = fCentValue;   // centrality
    vectorMC[3] = fMultiplicity;   // multiplicity (diff estimators can be used)
    vectorMC[4] = localmult;   // local multiplicity (Ntracks in R<0.4)
    vectorMC[5] = fq2;   // magnitude of reduced flow vector (computed using TPC tracks)
    break;
  }

  delete decay;
  bGenValues = kTRUE;
  return bGenValues;
}
//____________________________________________
Bool_t AliCFVertexingHFCascade::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{ 
  // read the variables for the container

  Bool_t bFillRecoValues = kFALSE;
  
  //Get the cascade and the neutral particle from it
  AliAODRecoCascadeHF* cascade = (AliAODRecoCascadeHF*)fRecoCandidate;
  AliAODRecoDecay* neutrDaugh = NULL; 
  if (fPDGcascade == 413) neutrDaugh = cascade->Get2Prong();
  else if (fPDGcascade == 4122) neutrDaugh = cascade->Getv0();
  else {
    return kFALSE;
  }

  //if (cascade->GetPrimaryVtx())printf("cascade has primary vtx\n");
  //if (fRecoCandidate->GetPrimaryVtx())printf("fRecoCandidateDstar has primary vtx\n");

  Double_t pt =  cascade->Pt();
  Double_t rapidity =  cascade->Y(fPDGcascade);
  // Double_t invMass = 0.;
  Double_t cosThetaStar = 9999.;
  Double_t pTneutrDaughPos = 0.;
  Double_t pTneutrDaughNeg = 0.;
  Double_t dca = neutrDaugh->GetDCA();
  // Double_t d0neutrDaughPos = 0.;
  // Double_t d0neutrDaughNeg = 0.;
  Double_t d0xd0 = neutrDaugh->Prodd0d0();
  Double_t cosPointingAngle = neutrDaugh->CosPointingAngle(fPrimVtx);
  Double_t phi = cascade->Phi();
  Double_t cosPointingAngleXY = neutrDaugh->CosPointingAngleXY(fPrimVtx);
  Double_t normDecayLengthXY = neutrDaugh->NormalizedDecayLengthXY(fPrimVtx);

  Int_t pdgCode = fmcPartCandidate->GetPdgCode();
 
  // UInt_t pdgDaughCascade[2] = { static_cast<UInt_t>(fPDGbachelor),  static_cast<UInt_t>(fPDGneutrDaugh) };    // bachelor is first daughter of cascade
  // UInt_t pdgDaughBarCascade[2] = { static_cast<UInt_t>(fPDGneutrDaugh),  static_cast<UInt_t>(fPDGbachelor) }; // bachelor is second daughter in case of a cascade-bar

  if (pdgCode > 0){
    cosThetaStar = neutrDaugh->CosThetaStar(1, fPDGneutrDaugh, fPDGneutrDaughPositive, fPDGneutrDaughNegative);
    pTneutrDaughPos = neutrDaugh->PtProng(0);
    pTneutrDaughNeg = neutrDaugh->PtProng(1);
    // d0neutrDaughPos = neutrDaugh->Getd0Prong(0);
    // d0neutrDaughNeg = neutrDaugh->Getd0Prong(1);
    // invMass = neutrDaugh->InvMass(2, pdgDaughCascade);
  }
  else {
    cosThetaStar = neutrDaugh->CosThetaStar(0, fPDGneutrDaugh, fPDGneutrDaughPositive, fPDGneutrDaughNegative);
    pTneutrDaughPos = neutrDaugh->PtProng(1);
    pTneutrDaughNeg = neutrDaugh->PtProng(0);
    // d0neutrDaughPos = neutrDaugh->Getd0Prong(1);
    // d0neutrDaughNeg = neutrDaugh->Getd0Prong(0);
    // invMass = neutrDaugh->InvMass(2, pdgDaughBarCascade);
  }
  
  Double_t cT = neutrDaugh->Ct(fPDGneutrDaugh, fPrimVtx);
  
  Int_t localmult = -1;
  if(fConfiguration == AliCFTaskVertexingHF::kESE) {
    localmult = ComputeLocalMultiplicity(cascade->Eta(), cascade->Phi(), 0.4);
  }

  switch (fConfiguration){
  case AliCFTaskVertexingHF::kSnail:
    vectorReco[0] = pt;
    vectorReco[1] = rapidity;
    vectorReco[2] = cosThetaStar;
    vectorReco[3] = pTneutrDaughPos;
    vectorReco[4] = pTneutrDaughNeg;
    vectorReco[5] = cT*1.E4;  // in micron
    vectorReco[6] = dca*1.E4;  // in micron
    vectorReco[7] = d0xd0*1.E8;  // in micron^2
    vectorReco[8] = cosPointingAngle;  // in micron
    vectorReco[9] = phi;
    vectorReco[10] = fzPrimVertex;    // z of reconstructed of primary vertex
    vectorReco[11] = fCentValue;
    vectorReco[12] = fFake;      // whether the reconstructed candidate was a fake (fFake = 0) or not (fFake = 2) 
    vectorReco[13] = cosPointingAngleXY;
    vectorReco[14] = normDecayLengthXY; // in cm
    vectorReco[15] = fMultiplicity; // reconstructed multiplicity
    break;
  case AliCFTaskVertexingHF::kCheetah:
    vectorReco[0] = pt;
    vectorReco[1] = rapidity;
    vectorReco[2] = cT*1.E4; // in micron
    vectorReco[3] = phi;
    vectorReco[4] = fzPrimVertex;
    vectorReco[5] = fCentValue;
    vectorReco[6] = fFake;
    vectorReco[7] = fMultiplicity; 
    break;
  case AliCFTaskVertexingHF::kFalcon:
    vectorReco[0] = pt;
    vectorReco[1] = rapidity;
    vectorReco[2] = fCentValue;
    vectorReco[3] = fMultiplicity; 
    break;
  case AliCFTaskVertexingHF::kESE:
    vectorReco[0] = pt;
    vectorReco[1] = rapidity;
    vectorReco[2] = fCentValue;   // centrality
    vectorReco[3] = fMultiplicity;   // multiplicity (diff estimators can be used)
    vectorReco[4] = localmult;   // local multiplicity (Ntracks in DeltaEta<0.1 and DeltaPhi<0.1)
    vectorReco[5] = fq2;   // magnitude of reduced flow vector (computed using TPC tracks)
    break;
  }

  bFillRecoValues = kTRUE;

  return bFillRecoValues;
}


//_____________________________________________________________
Bool_t AliCFVertexingHFCascade::CheckMCChannelDecay() const
{ 
  // check the required decay channel

  Bool_t checkCD = kFALSE;
  
  Int_t daughter0 = fmcPartCandidate->GetDaughterLabel(0);
  Int_t daughter1 = fmcPartCandidate->GetDaughterLabel(1);
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));

  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug (2,"Problems in the MC Daughters\n");
    return checkCD;
  }

  if(TMath::Abs(fmcPartCandidate->GetPdgCode()) == 4122 && (daughter1 - daughter0 != 1)) {
    AliDebug(2, Form("The MC particle doesn't have the correct daughters!!"));
    return checkCD;
  }

  if (!(TMath::Abs(mcPartDaughter0->GetPdgCode()) == fPDGneutrDaughForMC &&
	TMath::Abs(mcPartDaughter1->GetPdgCode()) == fPDGbachelor) && 
      !(TMath::Abs(mcPartDaughter0->GetPdgCode()) == fPDGbachelor &&
	TMath::Abs(mcPartDaughter1->GetPdgCode()) == fPDGneutrDaughForMC)) {
    AliDebug(2, Form("The cascade MC doesn't come from a the decay under study, skipping!! (Pdg codes of daughters = %d, %d)", mcPartDaughter0->GetPdgCode(), mcPartDaughter1->GetPdgCode()));
    return checkCD;  
  }
  
  // the Neutral Particle (e.g. D0 for D*, K0S for Lc...)
  AliAODMCParticle* mcPartDaughterNeutrDaugh = NULL;

  // for D* the D0 (the neutral) is the first daughter, while for Lc the V0 is the second, so we check the 
  // charge of the daughters to decide which is which
  AliDebug(3, Form("Charge0 = %d, Charge1 = %d", mcPartDaughter0->Charge()/3, mcPartDaughter1->Charge()/3));
  if (mcPartDaughter0->Charge()/3 != 0){
    mcPartDaughterNeutrDaugh = mcPartDaughter1;
  }
  else {
    mcPartDaughterNeutrDaugh = mcPartDaughter0;
  }

  Double_t vectorNeutrDaugh[2] ={0., 0.};

  // We are looking at a cascade ...evaluate the correct cascade
  if (!EvaluateIfCorrectNeutrDaugh(mcPartDaughterNeutrDaugh, vectorNeutrDaugh)) {
    AliDebug(2, "Error! the NeutrDaugh MC doesn't have correct daughters!!");
    return checkCD;  
  }
   
  checkCD = kTRUE;
  return checkCD;
  
}

//__________________________________________
Bool_t AliCFVertexingHFCascade::EvaluateIfCorrectNeutrDaugh(AliAODMCParticle* neutralDaugh, Double_t* vectorNeutrDaugh)const
{  
  //
  // check wether D0 is decaing into kpi
  //
  
  Bool_t isHadronic = kFALSE;
  AliDebug(2, Form("neutralDaugh = %p, pdg = %d", neutralDaugh, neutralDaugh->GetPdgCode()));

  if (fPDGcascade == 4122) {
    Int_t labelresonanceDaugh = neutralDaugh->GetDaughterLabel(0);
    AliAODMCParticle* resonanceDaugh = dynamic_cast<AliAODMCParticle*>(fmcArray->At(labelresonanceDaugh));
    if (!resonanceDaugh){
      return kFALSE;
    }
    else {
      AliDebug(3, Form("The daughter of the resonant particle is a %d (we are looking for a %d)", resonanceDaugh->GetPdgCode(), fPDGneutrDaugh));
      if (TMath::Abs(resonanceDaugh->GetPdgCode()) != fPDGneutrDaugh){
	return kFALSE;
      }
      else {
	neutralDaugh = resonanceDaugh;
      }
    }
  }

  Int_t daughterNeutrDaugh0 = neutralDaugh->GetDaughterLabel(0);
  Int_t daughterNeutrDaugh1 = neutralDaugh->GetDaughterLabel(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d", daughterNeutrDaugh0, daughterNeutrDaugh1));
  if (daughterNeutrDaugh0 == 0 || daughterNeutrDaugh1 == 0) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return isHadronic;  
  }
  
  Int_t numberOfExpectedDaughters = 2;
  if (TMath::Abs(daughterNeutrDaugh1 - daughterNeutrDaugh0) != numberOfExpectedDaughters-1) { // should be everytime true - see PDGdatabooklet
    AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
    return isHadronic;  
  }
  
  AliAODMCParticle* mcPartDaughterNeutrDaugh0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterNeutrDaugh0));
  AliAODMCParticle* mcPartDaughterNeutrDaugh1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterNeutrDaugh1));
  if (!mcPartDaughterNeutrDaugh0 || !mcPartDaughterNeutrDaugh1) {
    AliWarning("D0 MC analysis: At least one Daughter Particle not found in tree, skipping"); 
    return isHadronic;  
  }
  
  AliDebug(3, Form("Daughter 0 has pdg = %d, daughter 1 has pdg = %d", mcPartDaughterNeutrDaugh0->GetPdgCode(), mcPartDaughterNeutrDaugh1->GetPdgCode()));
  if (!(TMath::Abs(mcPartDaughterNeutrDaugh0->GetPdgCode()) == fPDGneutrDaughPositive &&
	TMath::Abs(mcPartDaughterNeutrDaugh1->GetPdgCode()) == fPDGneutrDaughNegative) && 
      !(TMath::Abs(mcPartDaughterNeutrDaugh0->GetPdgCode()) == fPDGneutrDaughNegative &&
	TMath::Abs(mcPartDaughterNeutrDaugh1->GetPdgCode()) == fPDGneutrDaughPositive)) {
    AliDebug(2, "The neutral particle (MC) doesn't come from the required decay, skipping!!");
    return isHadronic;  
  }
  
  Double_t sumPxDau = mcPartDaughterNeutrDaugh0->Px()+mcPartDaughterNeutrDaugh1->Px();
  Double_t sumPyDau = mcPartDaughterNeutrDaugh0->Py()+mcPartDaughterNeutrDaugh1->Py();
  Double_t sumPzDau = mcPartDaughterNeutrDaugh0->Pz()+mcPartDaughterNeutrDaugh1->Pz();
  Double_t pxMother = neutralDaugh->Px();
  Double_t pyMother = neutralDaugh->Py();
  Double_t pzMother = neutralDaugh->Pz();
  AliDebug(3, Form("pxMother = %f, pyMother = %f, pzMother = %f", pxMother, pyMother, pzMother));
  AliDebug(3, Form("sumPxDau = %f, sumPyDau = %f, sumPzDau = %f", sumPxDau, sumPyDau, sumPzDau));
  if(TMath::Abs(pxMother-sumPxDau)/(TMath::Abs(pxMother)+1.e-13)>fCutOnMomConservation ||
     TMath::Abs(pyMother-sumPyDau)/(TMath::Abs(pyMother)+1.e-13)>fCutOnMomConservation ||
     TMath::Abs(pzMother-sumPzDau)/(TMath::Abs(pzMother)+1.e-13)>fCutOnMomConservation){
    AliDebug(2, "Momentum conservation violated, skipping!!");
    return isHadronic;  
  }

  Double_t pTNeutrDaughPositive = 0;
  Double_t pTNeutrDaughNegative = 0;
  
  if (mcPartDaughterNeutrDaugh0->GetPdgCode() > 0 ) {
    pTNeutrDaughPositive = mcPartDaughterNeutrDaugh0->Pt();
    pTNeutrDaughNegative = mcPartDaughterNeutrDaugh1->Pt();
  }
  else {
    pTNeutrDaughPositive = mcPartDaughterNeutrDaugh1->Pt();
    pTNeutrDaughNegative = mcPartDaughterNeutrDaugh0->Pt();
  }
  
  isHadronic = kTRUE;
  
  vectorNeutrDaugh[0] = pTNeutrDaughPositive;
  vectorNeutrDaugh[1] = pTNeutrDaughNegative;
 
  return isHadronic;

}

//___________________________________________________________

void AliCFVertexingHFCascade::SetPtAccCut(Float_t* ptAccCut)
{
  //
  // setting the pt cut to be used in the Acceptance steps (MC+Reco)
  //

  AliDebug(3, "The 3rd element of the pt cut array will correspond to the cut applied to the soft pion - please check that it is correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fPtAccCut[iP]=ptAccCut[iP];
    }
  }
  return;
}		



//___________________________________________________________

void AliCFVertexingHFCascade::SetEtaAccCut(Float_t* etaAccCut)
{
  //
  // setting the eta cut to be used in the Acceptance steps (MC+Reco)
  //

  AliDebug(3, "The 3rd element of the eta cut array will correspond to the cut applied to the soft pion - please check that it is correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fEtaAccCut[iP] = etaAccCut[iP];
    }
  }
  return;
}
//___________________________________________________________

void AliCFVertexingHFCascade::SetAccCut(Float_t* ptAccCut, Float_t* etaAccCut)
{
  //
  // setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
  //

  AliDebug(3, "The 3rd element of the pt and cut array will correspond to the cut applied to the soft pion - please check that they are correct");
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fPtAccCut[iP]=ptAccCut[iP];
      fEtaAccCut[iP]=etaAccCut[iP];
    }
  }
  return;
}

//___________________________________________________________

void AliCFVertexingHFCascade::SetAccCut()
{
  //
  // setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
  //

  Int_t bachelorPosition = 2;
  if (fPDGcascade == 4122) bachelorPosition = 0;
  AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[bachelorPosition]));  // should be the soft pion...  
  if(!mcPartDaughter) return;
  Int_t mother =  mcPartDaughter->GetMother();
  AliAODMCParticle* mcMother = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother)); 
  if(!mcMother) return;
  
  if (TMath::Abs(mcPartDaughter->GetPdgCode()) != fPDGbachelor || TMath::Abs(mcMother->GetPdgCode()) != fPDGcascade){
    AliDebug(2, Form("Apparently the expected bachelor is not in the third position, causing an error (pdg expected = %d, actual = %d)!!", fPDGbachelor, mcPartDaughter->GetPdgCode()));
    AliDebug(2, "This should be fixed when checking the MC part family in the CF task...");
    return;  
  }	         
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs; iP++){
      fPtAccCut[iP]=0.1;
      fEtaAccCut[iP]=0.9;
    }
    
    if (fPDGcascade != 4122){
      fPtAccCut[2]=0.06;  // soft pion
      fEtaAccCut[2]=0.9;  // soft pion
    }
  }
  return;
}

//_____________________________________________________________
Double_t AliCFVertexingHFCascade::GetEtaProng(Int_t iProng) const 
{
  //
  // getting eta of the prong - overload the mother class method
  //

  if (fRecoCandidate){
   
    AliAODRecoCascadeHF* cascade = (AliAODRecoCascadeHF*)fRecoCandidate;

    Double_t etaProng =-9999;
    AliAODRecoDecay* neutrDaugh=0; 
    Int_t ibachelor = 0;
    if (fPDGcascade == 413) {
      neutrDaugh = cascade->Get2Prong();
    }
    else if (fPDGcascade == 4122) {
      neutrDaugh = cascade->Getv0();
    }
    if (iProng==0) etaProng = neutrDaugh->EtaProng(0);
    if (iProng==1) etaProng = neutrDaugh->EtaProng(1);
    if (iProng==2) etaProng = cascade->EtaProng(ibachelor);
    if (fPDGcascade == 4122){
      if (iProng==2) etaProng = neutrDaugh->EtaProng(0);
      if (iProng==1) etaProng = neutrDaugh->EtaProng(1);
      if (iProng==0) etaProng = cascade->EtaProng(ibachelor);
    }
    
    return etaProng;
    
  }
  return 999999;    
}
//_____________________________________________________________
Double_t AliCFVertexingHFCascade::GetPtProng(Int_t iProng) const 
{
  //
  // getting pt of the prong
  //
  
  if (fRecoCandidate){

    AliAODRecoCascadeHF* cascade = (AliAODRecoCascadeHF*)fRecoCandidate;
    Double_t ptProng= -9999;
    AliAODRecoDecay* neutrDaugh=0;
    Int_t ibachelor = 0;
    if (fPDGcascade == 413) {
      neutrDaugh = cascade->Get2Prong();
    }
    else if (fPDGcascade == 4122) {
      neutrDaugh = cascade->Getv0();
    }
    if (iProng == 0) ptProng = neutrDaugh->PtProng(0);
    if (iProng == 1) ptProng = neutrDaugh->PtProng(1);
    if (iProng == 2) ptProng = cascade->PtProng(ibachelor);
    if (fPDGcascade == 4122) {
      if (iProng == 2) ptProng = neutrDaugh->PtProng(0);
      if (iProng == 1) ptProng = neutrDaugh->PtProng(1);
      if (iProng == 0) ptProng = cascade->PtProng(ibachelor);
    }    
    //	Double_t ptProng = fRecoCandidate->PtProng(iProng);  
    return ptProng;
    
  }
  return 999999;  

}
//_____________________________________________________________
Bool_t AliCFVertexingHFCascade::CheckAdditionalCuts(AliPIDResponse* pidResponse) const {

  // function to check whether the candidate passes the additional cuts defined in the task to get the
  // invariant mass spectra; these cuts are NOT pt-dependent

  if (fPDGcascade == 4122){
    // the method is implemented only in this case so far
    AliAODRecoCascadeHF* cascade = (AliAODRecoCascadeHF*)fRecoCandidate;
    AliAODv0 * v0part = cascade->Getv0();
    AliAODTrack *bachelor = (AliAODTrack*)cascade->GetBachelor();
    Double_t bachelorEta = bachelor->Eta();
    AliAODTrack *v0pos = (AliAODTrack*)v0part->GetDaughter(0);
    AliAODTrack *v0neg = (AliAODTrack*)v0part->GetDaughter(1);
    Double_t v0posEta = v0pos->Eta();
    Double_t v0negEta = v0neg->Eta();

    Bool_t onFlyV0 = v0part->GetOnFlyStatus(); // on-the-flight V0s
    Double_t nSigmaTPCpr=-999.;
    Double_t nSigmaTOFpr=-999.;
    nSigmaTPCpr = pidResponse->NumberOfSigmasTPC(bachelor,(AliPID::kProton));
    nSigmaTOFpr = pidResponse->NumberOfSigmasTOF(bachelor,(AliPID::kProton));
    Double_t ptArm = v0part->PtArmV0();
    Double_t invmassK0s = v0part->MassK0Short();
    Double_t mK0SPDG = TDatabasePDG::Instance()->GetParticle(310)->Mass();

    Bool_t cutsForTMVA = (TMath::Abs(bachelorEta) < 0.8 && TMath::Abs(v0posEta) < 0.8 && TMath::Abs(v0negEta) < 0.8) &&
      ((nSigmaTOFpr < -800) || (TMath::Abs(nSigmaTOFpr) < 3)) &&  
      ((ptArm > 0.01 && ptArm < 0.07) || (ptArm > 0.105)) &&
      ((TMath::Abs(invmassK0s - mK0SPDG)) < 0.01);


    if (!fUseCutsForTMVA) cutsForTMVA = kTRUE;

    Bool_t cutsForInvMassTask = !(onFlyV0) && 
      (cascade->CosV0PointingAngle()>0.99) && 
      (TMath::Abs(nSigmaTPCpr) <= 3) && 
      (v0part->Getd0Prong(0) < 20) && 
      (v0part->Getd0Prong(1) < 20);

    if (cutsForTMVA && cutsForInvMassTask) {	
      // K0Smass cut
      // eta cut
      // TOF PID cut
      // Arm cut
      return kTRUE;
    }
  }

  return kFALSE;

}
