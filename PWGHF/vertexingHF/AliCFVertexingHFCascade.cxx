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

ClassImp(AliCFVertexingHFCascade)


//_________________________________________
  AliCFVertexingHFCascade::AliCFVertexingHFCascade(TClonesArray *mcArray, UShort_t originDselection):
  AliCFVertexingHF(mcArray, originDselection)
{	
  // standard constructor

  SetNProngs(3);
  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  // element 0 in the cut arrays corresponds to the soft pion!!!!!!!! Careful when setting the values...
  fPtAccCut[0]=0.;
  fEtaAccCut[0]=0.;
  for(Int_t iP=1; iP<fProngs; iP++){
	  fPtAccCut[iP]=0.1;
	  fEtaAccCut[iP]=0.9;
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
  AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)recoCand;

  if (!fRecoCandidate) {
    AliError("fRecoCandidate not found, problem in assignement\n");
    return bSignAssoc;
  }
  
  if ( fRecoCandidate->GetPrimaryVtx()) AliDebug(3,"fReco Candidate has a pointer to PrimVtx\n");
  //if (recoCand->GetPrimaryVtx()) printf("Reco Cand has a pointer to PrimVtx\n");
  
  //Int_t pdgCand = 413;

  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};

  Int_t nentries = fmcArray->GetEntriesFast();

  AliDebug(3,Form("nentries = %d\n", nentries));
 
  Int_t mcLabel =  dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,fmcArray); 
  
  if (mcLabel == -1) return bSignAssoc;

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
Bool_t AliCFVertexingHFCascade::GetGeneratedValuesFromMCParticle(Double_t* vectorMC) 
{
  // 
  // collecting all the necessary info (pt, y, cosThetaStar, ptPi, ptKa, cT) from MC particle
  //
	
	Bool_t bGenValues = kFALSE;
	
	
	Int_t daughter0ds = fmcPartCandidate->GetDaughter(0);
	Int_t daughter1ds = fmcPartCandidate->GetDaughter(1);

        //the D0
	AliAODMCParticle* mcPartDaughterD0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0ds));
	AliAODMCParticle* mcPartDaughterPis = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1ds));

	if(!mcPartDaughterD0 || !mcPartDaughterPis) return kFALSE;

	Double_t vtx1[3] = {0,0,0};   // primary vertex		
	Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
	Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
	fmcPartCandidate->XvYvZv(vtx1);  // cm

	//D0 daughters
	Int_t daughter0 = mcPartDaughterD0->GetDaughter(0);
	Int_t daughter1 = mcPartDaughterD0->GetDaughter(1);

	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0)); //D0
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1)); //pis

	if(!mcPartDaughter0 || !mcPartDaughter1) return kFALSE;

	// getting vertex from daughters
	mcPartDaughter0->XvYvZv(vtx2daughter0);  // cm
	mcPartDaughter1->XvYvZv(vtx2daughter1);  //cm
	if (TMath::Abs(vtx2daughter0[0] - vtx2daughter1[0])>1E-5 || TMath::Abs(vtx2daughter0[1]- vtx2daughter1[1])>1E-5 || TMath::Abs(vtx2daughter0[2] - vtx2daughter1[2])>1E-5) {
	  AliError("Daughters have different secondary vertex, skipping the track");
	  return bGenValues;
	}
	
	Int_t nprongs = 2;
	Short_t charge = 0;
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
	
	Double_t d0[2] = {0.,0.};		
	
	AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1,vtx2daughter0,nprongs,charge,px,py,pz,d0);
	
	Double_t cosThetaStar = 0.;
	Double_t cosThetaStarD0 = 0.;
	Double_t cosThetaStarD0bar = 0.;
	cosThetaStarD0 = decay->CosThetaStar(1,421,211,321);
	cosThetaStarD0bar = decay->CosThetaStar(0,421,321,211);
	if (mcPartDaughterD0->GetPdgCode() == 421){  // D0
	  AliDebug(3, Form("D0, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
	  cosThetaStar = cosThetaStarD0;
	}
	else if (mcPartDaughterD0->GetPdgCode() == -421){  // D0bar{
	  AliDebug(3, Form("D0bar, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
	  cosThetaStar = cosThetaStarD0bar;
	}
	else{
	  AliWarning("There are problems!! particle was expected to be either a D0 or a D0bar, check...");
	  return vectorMC;
	}
	if (cosThetaStar < -1 || cosThetaStar > 1) {
	  AliWarning("Invalid value for cosine Theta star, returning");
	  return bGenValues;
	}
	
	Double_t vectorD0[2] ={0.,0.};

	// evaluate the correct cascade
	if (!EvaluateIfD0toKpi(mcPartDaughterD0,vectorD0)) {
	  AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
	  return bGenValues;  
	}	

	//ct
	Double_t cT = decay->Ct(421);
	// get the pT of the daughters
	Double_t pTD0 = 0.;
	Double_t pTpi = 0.;
	
	if (TMath::Abs(fmcPartCandidate->GetPdgCode()) == 413) {
		pTD0 = mcPartDaughterD0->Pt();
		pTpi = mcPartDaughterPis->Pt();
	}

	
	switch (fConfiguration){
	case AliCFTaskVertexingHF::kSnail:
		vectorMC[0] = fmcPartCandidate->Pt();
		vectorMC[1] = fmcPartCandidate->Y() ;
		vectorMC[2] = cosThetaStar ;
		vectorMC[3] = vectorD0[0]; 
		vectorMC[4] = vectorD0[1];
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
	}

	delete decay;
	bGenValues = kTRUE;
	return bGenValues;
}
//____________________________________________
Bool_t AliCFVertexingHFCascade::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{ 
  // read the variables for the container

  Bool_t bFillRecoValues=kFALSE;
  
  //Get the D* and the D0 from D*
  AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)fRecoCandidate;
  AliAODRecoDecayHF2Prong* d0toKpi = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
  

  if (dstarD0pi->GetPrimaryVtx())printf("dstarD0pi has primary vtx\n");
  if (fRecoCandidate->GetPrimaryVtx())printf("fRecoCandidateDstar has primary vtx\n");

  Double_t pt =  dstarD0pi->Pt();
  Double_t rapidity =  dstarD0pi->YDstar();
  Double_t invMass=0.;
  Double_t cosThetaStar = 9999.;
  Double_t pTpi = 0.;
  Double_t pTK = 0.;
  Double_t dca = d0toKpi->GetDCA();
  Double_t d0pi = 0.;
  Double_t d0K = 0.;
  Double_t d0xd0 = d0toKpi->Prodd0d0();
  Double_t cosPointingAngle = d0toKpi->CosPointingAngle();
  Double_t phi = dstarD0pi->Phi();
  Double_t cosPointingAngleXY = d0toKpi->CosPointingAngleXY();
  Double_t normDecayLengthXY = d0toKpi->NormalizedDecayLengthXY();

  Int_t pdgCode = fmcPartCandidate->GetPdgCode();
 
  if (pdgCode > 0){
    cosThetaStar = d0toKpi->CosThetaStarD0();
    pTpi = d0toKpi->PtProng(0);
    pTK = d0toKpi->PtProng(1);
    d0pi = d0toKpi->Getd0Prong(0);
    d0K = d0toKpi->Getd0Prong(1);
    invMass=d0toKpi->InvMassD0();
  }
  else {
    cosThetaStar = d0toKpi->CosThetaStarD0bar();
    pTpi = d0toKpi->PtProng(1);
    pTK = d0toKpi->PtProng(0);
    d0pi = d0toKpi->Getd0Prong(1);
    d0K = d0toKpi->Getd0Prong(0);
    invMass= d0toKpi->InvMassD0bar();
  }
  
  Double_t cT = d0toKpi->CtD0();
  
	switch (fConfiguration){
	case AliCFTaskVertexingHF::kSnail:
		vectorReco[0] = pt;
		vectorReco[1] = rapidity;
		vectorReco[2] = cosThetaStar;
		vectorReco[3] = pTpi;
		vectorReco[4] = pTK;
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
		vectorReco[1] = rapidity ;
		vectorReco[2] = cT*1.E4; // in micron
		vectorReco[3] = phi; 
		vectorReco[4] = fzPrimVertex;
		vectorReco[5] = fCentValue;   
		vectorReco[6] = fFake ; 
		vectorReco[7] = fMultiplicity;  
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

  
  Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
  Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));

  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliDebug (2,"Problems in the MC Daughters\n");
    return checkCD;
  }

  if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==421 &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
      !(TMath::Abs(mcPartDaughter0->GetPdgCode())==211 &&
	TMath::Abs(mcPartDaughter1->GetPdgCode())==421)) {
    AliDebug(2, "The D0 MC doesn't come from a Kpi decay, skipping!!");
    return checkCD;  
  }
  
  Double_t vectorD0[2] ={0.,0.};

  // D* is a cascade ...evaluate the correct cascade
  if (!EvaluateIfD0toKpi(mcPartDaughter0,vectorD0)) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return checkCD;  
  }
   
  checkCD = kTRUE;
  return checkCD;
  
}

//__________________________________________
Bool_t AliCFVertexingHFCascade::EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, Double_t* vectorD0)const
{  
  //
  // chack wether D0 is decaing into kpi
  //
  
  Bool_t isHadronic = kFALSE;
  
  Int_t daughterD00 = neutralDaugh->GetDaughter(0);
  Int_t daughterD01 = neutralDaugh->GetDaughter(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughterD00,daughterD01));
  if (daughterD00 == 0 || daughterD01 == 0) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return isHadronic;  
  }
  
  if (TMath::Abs(daughterD01 - daughterD00) != 1) { // should be everytime true - see PDGdatabooklet
    AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
    return isHadronic;  
  }
  
  AliAODMCParticle* mcPartDaughterD00 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD00));
  AliAODMCParticle* mcPartDaughterD01 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughterD01));
  if (!mcPartDaughterD00 || !mcPartDaughterD01) {
    AliWarning("D0 MC analysis: At least one Daughter Particle not found in tree, skipping"); 
    return isHadronic;  
  }
  
  if (!(TMath::Abs(mcPartDaughterD00->GetPdgCode())==321 &&
	TMath::Abs(mcPartDaughterD01->GetPdgCode())==211) && 
      !(TMath::Abs(mcPartDaughterD00->GetPdgCode())==211 &&
	TMath::Abs(mcPartDaughterD01->GetPdgCode())==321)) {
    AliDebug(2, "The D0 MC doesn't come from a Kpi decay, skipping!!");
    return isHadronic;  
  }
  
  Double_t pTD0pi = 0;
  Double_t pTD0K = 0;
  
  
  if (TMath::Abs(mcPartDaughterD00->GetPdgCode()) == 211) {
    pTD0pi = mcPartDaughterD00->Pt();
    pTD0K = mcPartDaughterD01->Pt();
  }
  else {
    pTD0pi = mcPartDaughterD01->Pt();
    pTD0K  = mcPartDaughterD00->Pt();
  }
  
  isHadronic = kTRUE;
  
  vectorD0[0] = pTD0pi;
  vectorD0[1] = pTD0K;
 
  return isHadronic;

}

//___________________________________________________________

void AliCFVertexingHFCascade::SetPtAccCut(Float_t* ptAccCut)
{
	//
	// setting the pt cut to be used in the Acceptance steps (MC+Reco)
	//

	AliInfo("The 3rd element of the pt cut array will correspond to the cut applied to the soft pion - please check that it is correct");
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

	AliInfo("The 3rd element of the eta cut array will correspond to the cut applied to the soft pion - please check that it is correct");
	if (fProngs>0){
		for (Int_t iP=0; iP<fProngs; iP++){
			fEtaAccCut[iP]=etaAccCut[iP];
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

	AliInfo("The 3rd element of the pt and cut array will correspond to the cut applied to the soft pion - please check that they are correct");
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
  
  AliAODMCParticle* mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[2]));  // should be the soft pion...  
  if(!mcPartDaughter) return;
  Int_t mother =  mcPartDaughter->GetMother();
  AliAODMCParticle* mcMother = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother)); 
  if(!mcMother) return;

  if (TMath::Abs(mcPartDaughter->GetPdgCode())!= 211 || TMath::Abs(mcMother->GetPdgCode())!=413){
    AliFatal("Apparently the soft pion is not in the third position, causing a crash!!");
  }	         
  if (fProngs>0){
    for (Int_t iP=0; iP<fProngs-1; iP++){
      fPtAccCut[iP]=0.1;
      fEtaAccCut[iP]=0.9;
    }
    fPtAccCut[2]=0.06;  // soft pion
    fEtaAccCut[2]=0.9;  // soft pion
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
   
   AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)fRecoCandidate;

   Double_t etaProng =-9999;
    if(iProng==0) etaProng =dstarD0pi->Get2Prong()->EtaProng(0);
    if(iProng==1) etaProng =dstarD0pi->Get2Prong()->EtaProng(1);
    if(iProng==2) etaProng =dstarD0pi->EtaProng(1);
    
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

    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)fRecoCandidate;
    Double_t ptProng= -9999;
    if(iProng==0) ptProng =dstarD0pi->Get2Prong()->PtProng(0);
    if(iProng==1) ptProng =dstarD0pi->Get2Prong()->PtProng(1);
    if(iProng==2) ptProng =dstarD0pi->PtProng(1);
    
    //	Double_t ptProng = fRecoCandidate->PtProng(iProng);  
    return ptProng;
    
  }
  return 999999;  
  
}
