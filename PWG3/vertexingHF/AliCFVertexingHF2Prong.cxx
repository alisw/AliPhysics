/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
// Class for HF corrections as a function of many variables and step 
// Author : C. Zampolli, CERN
// D. Caffarri, Univ & INFN Padova caffarri@pd.infn.it
// Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------

#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "TClonesArray.h"
#include "AliCFVertexingHF.h"
#include "AliESDtrack.h"
#include "TDatabasePDG.h"

#include "AliCFVertexingHF2Prong.h"
#include "AliCFContainer.h"
#include "AliCFTaskVertexingHF.h"

ClassImp(AliCFVertexingHF2Prong)


//_________________________________________
  AliCFVertexingHF2Prong::AliCFVertexingHF2Prong(TClonesArray *mcArray, UShort_t originDselection):
	  AliCFVertexingHF(mcArray, originDselection)
{	
	//
	// constructor
	//

	SetNProngs(2);
	fPtAccCut=new Float_t[fProngs];
	fEtaAccCut=new Float_t[fProngs];
	for(Int_t iP=0; iP<fProngs; iP++){
		fPtAccCut[iP]=0.1;
		fEtaAccCut[iP]=0.9;
	}

}


//_____________________________________
AliCFVertexingHF2Prong& AliCFVertexingHF2Prong::operator=(const AliCFVertexingHF2Prong& c)
{
	//
	// copy constructor	
	//

	if  (this != &c) {		
		AliCFVertexingHF::operator=(c);		
	}
	return *this;
}

//__________________________________________
Bool_t AliCFVertexingHF2Prong::SetRecoCandidateParam(AliAODRecoDecayHF *recoCand)
{  
	//
	// setting the recontructed candidate
	//
	
	Bool_t bSignAssoc = kFALSE;	
	fRecoCandidate = recoCand;
	if (!fRecoCandidate) {
		AliError("fRecoCandidate not found, problem in assignement\n");
		return bSignAssoc;
	}
	
	if (fRecoCandidate->GetPrimaryVtx()) AliDebug(3,"fReco Candidate has a pointer to PrimVtx\n");
	if (recoCand->GetPrimaryVtx()) AliDebug(3,"Reco Cand has a pointer to PrimVtx\n");
	
	Int_t pdgCand = 421;
	Int_t pdgDgD0toKpi[2]={321,211};
	Int_t nentries = fmcArray->GetEntriesFast();	

	AliDebug(3,Form("nentries = %d\n", nentries));
 
	Int_t mcLabel = fRecoCandidate->MatchToMC(pdgCand,fmcArray,2,pdgDgD0toKpi);
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
Bool_t AliCFVertexingHF2Prong::GetGeneratedValuesFromMCParticle(Double_t* vectorMC) 
{
	// 
	// collecting all the necessary info (pt, y, cosThetaStar, ptPi, ptKa, cT) from MC particle
	//
	
	Bool_t bGenValues = kFALSE;
	Double_t vtx1[3] = {0,0,0};   // primary vertex		
	Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
	Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
	fmcPartCandidate->XvYvZv(vtx1);  // cm

	Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
	Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));
	if(!mcPartDaughter0 || !mcPartDaughter1) return bGenValues;

	// getting vertex from daughters
	mcPartDaughter0->XvYvZv(vtx2daughter0);  // cm
	mcPartDaughter1->XvYvZv(vtx2daughter1);  //cm
	if (TMath::Abs(vtx2daughter0[0] - vtx2daughter1[0]) > 1E-5 || TMath::Abs(vtx2daughter0[1] - vtx2daughter1[1]) > 1E-5 || TMath::Abs(vtx2daughter0[2] - vtx2daughter1[2])>1E-5) {
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
	if (fmcPartCandidate->GetPdgCode() == 421){  // D0
		AliDebug(3, Form("D0, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
		cosThetaStar = cosThetaStarD0;
	}
	else if (fmcPartCandidate->GetPdgCode() == -421){  // D0bar{
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
		
	//ct
	Double_t cT = decay->Ct(421);
	// get the pT of the daughters
	Double_t pTpi = 0.;
	Double_t pTK = 0.;
	
	if (TMath::Abs(mcPartDaughter0->GetPdgCode()) == 211) {
		pTpi = mcPartDaughter0->Pt();
		pTK = mcPartDaughter1->Pt();
	}
	else {
		pTpi = mcPartDaughter1->Pt();
		pTK = mcPartDaughter0->Pt();
	}
	
	switch (fConfiguration){
	case AliCFTaskVertexingHF::kSnail:
		vectorMC[0] = fmcPartCandidate->Pt();
		vectorMC[1] = fmcPartCandidate->Y() ;
		vectorMC[2] = cosThetaStar ;
		vectorMC[3] = pTpi ;
		vectorMC[4] = pTK ;
		vectorMC[5] = cT*1.E4 ;  // in micron
		vectorMC[6] = 0.;   // dummy value for dca, meaningless in MC
		vectorMC[7] = -80000.; // dummy value for d0pixd0K, meaningless in MC, in micron^2
		vectorMC[8] = 1.01;    // dummy value for cosPointing, meaningless in MC
		vectorMC[9] = fmcPartCandidate->Phi(); 
		vectorMC[10] = fzMCVertex;    // z of reconstructed of primary vertex
		vectorMC[11] = fCentValue;   //reconstructed centrality 
		vectorMC[12] = 1.;           // fake: always filling with 1 at MC level 
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
Bool_t AliCFVertexingHF2Prong::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{
	//
	// Getting the reconstructed values from the candidate
	// 
	
	Bool_t bFillRecoValues=kFALSE;
	
	AliAODRecoDecayHF2Prong *d0toKpi = (AliAODRecoDecayHF2Prong*)fRecoCandidate;
	
	if (d0toKpi->GetPrimaryVtx())AliDebug(3,"d0toKpi has primary vtx\n");
	if (fRecoCandidate->GetPrimaryVtx())AliDebug(3,"fRecoCandidate has primary vtx\n");
	
	Double_t pt = d0toKpi->Pt();
	Double_t rapidity = d0toKpi->YD0();
	Double_t invMass=0.;
	Double_t cosThetaStar = 9999.;
	Double_t pTpi = 0.;
	Double_t pTK = 0.;
	Double_t dca = d0toKpi->GetDCA();
	Double_t d0pi = 0.;
	Double_t d0K = 0.;
	Double_t d0xd0 = d0toKpi->Prodd0d0();
	Double_t cosPointingAngle = d0toKpi->CosPointingAngle();
	Double_t phi = d0toKpi->Phi();
	Int_t pdgCode = fmcPartCandidate->GetPdgCode();
	Double_t cosPointingAngleXY = d0toKpi->CosPointingAngleXY();
	Double_t normDecayLengthXY = d0toKpi->NormalizedDecayLengthXY();
       
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
		vectorReco[8] = cosPointingAngle; 
		vectorReco[9] = phi;  
		vectorReco[10] = fzPrimVertex;    // z of reconstructed of primary vertex
		vectorReco[11] = fCentValue; //reconstructed centrality 
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
Bool_t AliCFVertexingHF2Prong::CheckMCChannelDecay() const
{ 
	//
	// checking the MC decay channel
	//
	Bool_t checkCD = kFALSE;
	Int_t daughter0 = fmcPartCandidate->GetDaughter(0);
	Int_t daughter1 = fmcPartCandidate->GetDaughter(1);
	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter0));
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter1));
	
	if (!mcPartDaughter0 || !mcPartDaughter1) {
		AliDebug (2,"Problems in the MC Daughters\n");
		return checkCD;
	}
	
	if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==321 &&
	      TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
	    !(TMath::Abs(mcPartDaughter0->GetPdgCode())==211 &&
	      TMath::Abs(mcPartDaughter1->GetPdgCode())==321)) {
		AliDebug(2, "The D0 MC doesn't come from a Kpi decay, skipping!!");
		return checkCD;  
	}
	
	checkCD = kTRUE;
	return checkCD;
	
}

