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

#include "TParticle.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"

#include "AliCFVertexingHF.h"

//___________________________________________________________
AliCFVertexingHF::AliCFVertexingHF() :
	fmcArray(0x0),
	fRecoCandidate(0),
	fmcPartCandidate(0x0),
	fNDaughters(0),
	fNVar(0),
	fzPrimVertex(0),
	fzMCVertex(0),
	fFillFromGenerated(0),
	fOriginDselection(0),
	fKeepDfromB(kFALSE),
	fKeepDfromBOnly(kFALSE),
	fmcLabel(0),
	fProngs(-1)
{
	return;
}



//_____________________________________________________
AliCFVertexingHF::AliCFVertexingHF(TClonesArray *mcArray, UShort_t originDselection) :
	fmcArray(mcArray),
	fRecoCandidate(0),
	fmcPartCandidate(0x0),
	fNDaughters(0),
	fNVar(0),
	fzPrimVertex(0),
	fzMCVertex(0),
	fFillFromGenerated(0),
	fOriginDselection(0),
	fKeepDfromB(kFALSE),
	fKeepDfromBOnly(kFALSE),
	fmcLabel(0),
	fProngs(-1)
{

  
  SetDselection(originDselection);


  return;
}



//_______________________________________________________
AliCFVertexingHF::~AliCFVertexingHF()
{
	if (fmcArray) delete fmcArray;
	if (fRecoCandidate) delete fRecoCandidate;
	if (fmcPartCandidate) delete fmcPartCandidate;
}


//_____________________________________________________
AliCFVertexingHF& AliCFVertexingHF::operator=(const AliCFVertexingHF& c){
	
	if (this!= &c){
		TObject::operator=(c);
		fmcArray = c.fmcArray;
		fRecoCandidate = c.fRecoCandidate;
		fmcPartCandidate = c.fmcPartCandidate;
		fNDaughters = c.fNDaughters;
		fNVar = c.fNVar;
		fzPrimVertex = c.fzPrimVertex;
		fzMCVertex = c.fzMCVertex;
		fFillFromGenerated = c.fFillFromGenerated;
		fOriginDselection = c.fOriginDselection;
		fKeepDfromB = c.fKeepDfromB;
		fKeepDfromBOnly = c.fKeepDfromBOnly;
		fmcLabel = c.fmcLabel;
		
	}
	
	return *this;
}

//____________________________________________________
AliCFVertexingHF::AliCFVertexingHF(const AliCFVertexingHF &c) :
        TObject(c),
	fmcArray(c.fmcArray),
	fRecoCandidate(c.fRecoCandidate),
	fmcPartCandidate(c.fmcPartCandidate),
	fNDaughters(c.fNDaughters),
	fNVar(c.fNVar),
	fzPrimVertex(c.fzPrimVertex),
	fzMCVertex(c.fzMCVertex),
	fFillFromGenerated(c.fFillFromGenerated),
	fOriginDselection (c.fOriginDselection),
	fKeepDfromB (c.fKeepDfromB),
	fKeepDfromBOnly (c.fKeepDfromBOnly),
	fmcLabel(c.fmcLabel),
	fProngs(c.fProngs)
{
   

}

//___________________________________________________________
void AliCFVertexingHF::SetDselection(UShort_t originDselection){

 fOriginDselection = originDselection;
 
 if (fOriginDselection == 0) {
   fKeepDfromB = kFALSE;
   fKeepDfromBOnly = kFALSE;
 }
 
 if (fOriginDselection == 1) {
   fKeepDfromB = kTRUE;
   fKeepDfromBOnly = kTRUE;
 }
 
 if (fOriginDselection == 2) {
   fKeepDfromB = kTRUE;
   fKeepDfromBOnly = kFALSE;
 }
 
 return;

}

//______________________________________________________
void AliCFVertexingHF::SetMCCandidateParam(Int_t label){
	
  fmcPartCandidate = dynamic_cast <AliAODMCParticle*> (fmcArray->At(label));
  fNDaughters = fmcPartCandidate->GetNDaughters();
  return;
}


//____________________________________________________________
Int_t AliCFVertexingHF::MCcquarkCounting(AliAODMCParticle* mcPart) const{
	
 
  Int_t cquarks = 0;
  if (mcPart->GetPdgCode() == 4) cquarks++; 
  if (mcPart->GetPdgCode() == -4) cquarks++; 
  if (!mcPart) {
    AliWarning("Particle not found in tree, skipping\n"); 
    return cquarks;
  } 
  
  return cquarks;
}


//________________________________________________________
Bool_t AliCFVertexingHF::CheckMCPartFamily(AliAODMCParticle */*mcPart*/, TClonesArray */*mcArray*/) const {

  Int_t pdgGranma = CheckOrigin();
  if (pdgGranma == -9999){
    printf ("This particle come from a B decay channel but the we keep only the prompt charm particles\n");	
    return kFALSE;
  }	
  
 if (pdgGranma == -999){
    printf ("This particle come from a prompt charm particles but we want only the ones coming from B\n");	
    return kFALSE;
  }	

  if (!CheckMCDaughters()) return kFALSE;
  if (!CheckMCChannelDecay()) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________________________________
Int_t AliCFVertexingHF::CheckOrigin() const {		
  //
  // checking whether the mother of the particles come from a charm or a bottom quark
  //
  
  Int_t pdgGranma = 0;
  Int_t mother = 0;
  mother = fmcPartCandidate->GetMother();
  Int_t istep = 0;
  while (mother >0 ){
    istep++;
    AliDebug(2,Form("mother at step %d = %d", istep, mother));
    AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother));
    pdgGranma = mcGranma->GetPdgCode();
    AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
    Int_t abspdgGranma = TMath::Abs(pdgGranma);
    if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
      if (!fKeepDfromB) return -9999; //skip particle if come from a B meson.
      
      else{
	return pdgGranma;
      }
    }
    else {
    if (fKeepDfromBOnly) return -999;
    }
    
    mother = mcGranma->GetMother();
  }
  
  return pdgGranma;
}


//___________________________________________
Bool_t AliCFVertexingHF::CheckMCDaughters()const {
	
  AliAODMCParticle *mcPartDaughter;
  Bool_t checkDaughters = kFALSE;

  Int_t label0 = fmcPartCandidate->GetDaughter(0);
  Int_t label1 = fmcPartCandidate->GetDaughter(1);
  if (label1==0 || label0 == 0){
    AliDebug(2, Form("The MC particle doesn't jave correct daughters, skipping!!"));
	return checkDaughters;  
  }

  if (label1 - label0 != fProngs-1){
    AliDebug(2, Form("The MC particle doesn't come from a %d-prong decay, skipping!!", fProngs));
    return checkDaughters;  
  }

  for (Int_t iProng = 0; iProng<fProngs; iProng++){
    mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(label0+iProng));    
    if (!mcPartDaughter) {
      AliWarning("At least one Daughter Particle not found in tree, skipping"); 
      return checkDaughters;  
    }
  }

  checkDaughters = kTRUE;
  return checkDaughters;
}

	

//______________________________________________________
Bool_t AliCFVertexingHF::FillMCContainer(Double_t *containerInputMC)
{
// fill the container for Generator level selection
  Bool_t mcContainerFilled = kFALSE;  

  Double_t* vectorMC = new Double_t[fNVar];
  for (Int_t iVar = 0; iVar<fNVar; iVar++) vectorMC[iVar]= 9999.;
  
  if (GetGeneratedValuesFromMCParticle(&vectorMC[0])){
    for (Int_t iVar = 0; iVar<fNVar; iVar++){
      
      containerInputMC[iVar] = vectorMC[iVar];
    }
    
    mcContainerFilled = kTRUE;		
  }
  delete vectorMC;
  return mcContainerFilled;	
}

//______________________________________________________
Bool_t AliCFVertexingHF::FillRecoContainer(Double_t *containerInput) {  

// fill the container for Reconstrucred level selection

  Bool_t recoContainerFilled = kFALSE;
  Double_t* vectorValues = new Double_t[fNVar];
  Double_t* vectorReco = new Double_t[fNVar];  
  
  for (Int_t iVar = 0; iVar<fNVar; iVar++) {
    vectorValues[iVar]= 9999.;
    vectorReco[iVar]=9999.;
  }

  if(fFillFromGenerated){
    //filled with MC values
    if (GetGeneratedValuesFromMCParticle(&vectorValues[0])){
      for (Int_t iVar = 0; iVar<fNVar; iVar++){
		  containerInput[iVar] = vectorValues[iVar];
		}
		recoContainerFilled = kTRUE;		
		}
	}
	else{
    //filled with Reco values
		
		if (GetRecoValuesFromCandidate(&vectorReco[0])){
			for (Int_t iVar = 0; iVar<fNVar; iVar++){
				containerInput[iVar] = vectorReco[iVar];
			}
			recoContainerFilled = kTRUE;		
		}
	}
  
  delete vectorValues;
  delete vectorReco;
  return recoContainerFilled;	
}

//_____________________________________________________
Bool_t AliCFVertexingHF::MCAcceptanceStep() const
{

  Bool_t bMCAccStep = kFALSE;
  
  AliAODMCParticle *mcPartDaughter;
  Int_t label0 = fmcPartCandidate->GetDaughter(0);
  Int_t label1 = fmcPartCandidate->GetDaughter(1);
  if (label1==0 || label0 == 0){
    AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
    return bMCAccStep;  
  }

  if (label1 - label0 != fProngs-1){
    AliDebug(2, Form("The MC particle doesn't come from a %d-prong decay, skipping!!", fProngs));
	return bMCAccStep;  
  }

  for (Int_t iProng = 0; iProng<fProngs; iProng++){
    mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(label0+iProng));    
    if (!mcPartDaughter) {
      AliWarning("At least one Daughter Particle not found in tree, skipping"); 
      return bMCAccStep;  
    }
    Double_t eta = mcPartDaughter->Eta();
    Double_t pt = mcPartDaughter->Pt();
	  
	  //set values of eta and pt in the constructor.
    if (TMath::Abs(eta) > 0.9 || pt < 0.1){
		AliDebug(3,"At least one daughter has eta>0.9 or pt < 0.1 \n"); 
      return bMCAccStep;
    }

  }

  
  bMCAccStep = kTRUE;
  return bMCAccStep; 
  
}

//_____________________________________________________
Bool_t AliCFVertexingHF::MCRefitStep(AliAODEvent *aodEvent, AliESDtrackCuts *trackCuts) const
{		
  // check on the kTPCrefit and kITSrefit conditions of the daughters
	Bool_t bRefitStep = kFALSE;

	Int_t label0 = fmcPartCandidate->GetDaughter(0);
	Int_t label1 = fmcPartCandidate->GetDaughter(1);

	if (label1==0 || label0 == 0){
		AliDebug(2, Form("The MC particle doesn't jave correct daughters, skipping!!"));
        return bRefitStep;  
	}

	if (label1 - label0 != fProngs-1){
		AliDebug(2, Form("The MC particle doesn't come from a %d-prong decay, skipping!!", fProngs));
		//AliInfo(Form("The MC particle doesn't come from a %d-prong decay, skipping!!", fProngs));
		return bRefitStep;  
	}

	Int_t foundDaughters = 0;

	if (trackCuts->GetRequireTPCRefit() || trackCuts->GetRequireITSRefit()){
    
		for(Int_t iaod =0; iaod<aodEvent->GetNumberOfTracks(); iaod++){
			AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iaod);
			if (track->GetLabel()>= label0 && track->GetLabel()<=label1){
				foundDaughters++;
				printf("daughter %d \n",foundDaughters);
				if(trackCuts->GetRequireTPCRefit()){
					if((track->GetStatus()&AliESDtrack::kTPCrefit)) {
						bRefitStep = kTRUE;
					}
					else {
						AliDebug(3, "Refit cut not passed , missing TPC refit\n");
						AliInfo( "Refit cut not passed , missing TPC refit\n");
						return kFALSE;
					}
				}
				
				if (trackCuts->GetRequireITSRefit()) {
					if((track->GetStatus()&AliESDtrack::kITSrefit)){
						bRefitStep = kTRUE;
					}
					else {
						AliDebug(3, "Refit cut not passed , missing ITS refit\n");
						AliInfo("Refit cut not passed , missing ITS refit\n");
						return kFALSE;
					}
				}
			}      	
			if (foundDaughters == fProngs){
	 
				break;
			}
      
		}
    
	}
  
	return bRefitStep;
}

//____________________________________________________________________________

Bool_t AliCFVertexingHF::RecoStep() 
//check also vertex and ITS Refit and TPC Refit
{ 
  Bool_t bRecoStep = kFALSE;
  Int_t mcLabel = GetMCLabel();
  
  if (mcLabel == -1) {
    AliDebug(2,"No MC particle found");
    return bRecoStep;
  }
  else{
    fmcPartCandidate = (AliAODMCParticle*)fmcArray->At(mcLabel);
    if (!fmcPartCandidate){
      AliWarning("Could not find associated MC in AOD MC tree");
      return bRecoStep;
    }
    
  }
  
  Int_t pdgGranma = CheckOrigin();
  
  if (pdgGranma == -9999){
    printf ("This particle come from a B decay channel but we keep only prompt charm particles\n");
    return bRecoStep;
  }

  if (pdgGranma == -999){
    printf ("This particle come from a  prompt charm particle but we want only the ones coming from B\n");
    return bRecoStep;
  }
  
   
  bRecoStep=kTRUE;
  return bRecoStep;
  
  
}	
//____________________________________________
Double_t AliCFVertexingHF::GetEtaProng(Int_t iProng) const {

  if (fRecoCandidate){
    Double_t etaProng = fRecoCandidate->EtaProng(iProng);  
    return etaProng;
  }
  return 999999;  
}
//______________________________________________________
Double_t AliCFVertexingHF::GetPtProng(Int_t iProng) const {

  if (fRecoCandidate){
    Double_t ptProng = fRecoCandidate->PtProng(iProng);  
    return ptProng;
  }
  return 999999;  
  
}

//____________________________________________________________________

Bool_t AliCFVertexingHF::RecoAcceptStep(AliESDtrackCuts *trackCuts) const
{
  Bool_t bRecoAccStep = kFALSE;
  
  Float_t etaCutMin, ptCutMin, etaCutMax, ptCutMax;
  trackCuts->GetEtaRange(etaCutMin, etaCutMax);
  trackCuts->GetPtRange(ptCutMin, ptCutMax);
  
  Float_t etaProng=0., ptProng=0.; 
  
  for (Int_t iProng =0; iProng<fProngs; iProng++){
    
    etaProng = GetEtaProng(iProng);
    ptProng = GetPtProng(iProng);
    
    Bool_t acceptanceProng = (etaProng>etaCutMin && etaProng<etaCutMax && ptProng>ptCutMin && ptProng<ptCutMax);
    if (!acceptanceProng) {
      printf ("At least one reconstructed prong isn't in the acceptance\n");
      return bRecoAccStep;
    }
  }
  
  bRecoAccStep=kTRUE;
  return bRecoAccStep;
}
//___________________________________________________________

Bool_t AliCFVertexingHF::FillUnfoldingMatrix(Double_t *fill) const{

  fill = new Double_t[4];
  
  if(fmcPartCandidate){
    
    fill[0] = GetPtCand();
    fill[1] = GetYCand();
    
    fill[2] =  fmcPartCandidate->Pt(); 
    fill[3] =  fmcPartCandidate->Y(); 
    
    return kTRUE;
  }

  delete fill;
  return kFALSE;
}


//___________________________________________________________
