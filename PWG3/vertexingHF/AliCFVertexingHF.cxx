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

// $Id$

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
#include "AliAODRecoCascadeHF.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliCFTaskVertexingHF.h"

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
	fProngs(-1),
	fLabelArray(0x0), 
	fCentValue(0.),
	fPtAccCut(0x0),
	fEtaAccCut(0x0),
	fFakeSelection(0),
	fFake(1.), // setting to MC value
	fRejectIfNoQuark(kFALSE),
	fMultiplicity(0.),
	fConfiguration(AliCFTaskVertexingHF::kCheetah) // by default, setting the fast configuration
{
	//
	// constructor
	//


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
	fProngs(-1),
	fLabelArray(0x0),
	fCentValue(0.),
	fPtAccCut(0x0),
	fEtaAccCut(0x0),
	fFakeSelection(0),
	fFake(1.), // setting to MC value
	fRejectIfNoQuark(kFALSE),
	fMultiplicity(0.),
	fConfiguration(AliCFTaskVertexingHF::kCheetah) // by default, setting the fast configuration
{
	//
	// constructor with mcArray
	//
	
	SetDselection(originDselection);
	return;
}

//_______________________________________________________
AliCFVertexingHF::~AliCFVertexingHF()
{
	//
	// destructor
	//

	if (fmcArray) fmcArray = 0x0;
	if (fRecoCandidate) fRecoCandidate = 0x0;
	if (fmcPartCandidate) fmcPartCandidate = 0x0;
	if (fLabelArray){
	  	delete [] fLabelArray;
	  	fLabelArray = 0x0;
	}	
	if (fPtAccCut){
	  	delete [] fPtAccCut;
	  	fPtAccCut = 0x0;
	}	
	if (fEtaAccCut){
	  	delete [] fEtaAccCut;
	  	fEtaAccCut = 0x0;
	}	
}

//_____________________________________________________
AliCFVertexingHF& AliCFVertexingHF::operator=(const AliCFVertexingHF& c)
{	
	//
	// assigment operator
	//

	if (this!= &c){
		TObject::operator=(c);
		fmcArray = new TClonesArray(*(c.fmcArray));
		fRecoCandidate = new AliAODRecoDecayHF(*(c.fRecoCandidate));
		fmcPartCandidate = new AliAODMCParticle(*(c.fmcPartCandidate));
		fNDaughters = c.fNDaughters;
		fNVar = c.fNVar;
		fzPrimVertex = c.fzPrimVertex;
		fzMCVertex = c.fzMCVertex;
		fFillFromGenerated = c.fFillFromGenerated;
		fOriginDselection = c.fOriginDselection;
		fKeepDfromB = c.fKeepDfromB;
		fKeepDfromBOnly = c.fKeepDfromBOnly;
		fmcLabel = c.fmcLabel;
		fProngs=c.fProngs;
		fCentValue=c.fCentValue;
		fFakeSelection=c.fFakeSelection;
		fFake=c.fFake;
		fRejectIfNoQuark=c.fRejectIfNoQuark;
		if (fProngs > 0){
			fLabelArray = new Int_t[fProngs];
                        fPtAccCut = new Float_t[fProngs];
                        fEtaAccCut = new Float_t[fProngs];
			for(Int_t iP=0; iP<fProngs; iP++){
				fLabelArray[iP]=c.fLabelArray[iP];
				fPtAccCut[iP]=c.fPtAccCut[iP];
				fEtaAccCut[iP]=c.fEtaAccCut[iP];
			}
		}
		fMultiplicity=c.fMultiplicity;
		fConfiguration=c.fConfiguration;
	}
	
	return *this;
}

//____________________________________________________
AliCFVertexingHF::AliCFVertexingHF(const AliCFVertexingHF &c) :
        TObject(c),
	fNDaughters(c.fNDaughters),
	fNVar(c.fNVar),
	fzPrimVertex(c.fzPrimVertex),
	fzMCVertex(c.fzMCVertex),
	fFillFromGenerated(c.fFillFromGenerated),
	fOriginDselection (c.fOriginDselection),
	fKeepDfromB (c.fKeepDfromB),
	fKeepDfromBOnly (c.fKeepDfromBOnly),
	fmcLabel(c.fmcLabel),
	fProngs(c.fProngs),
	fLabelArray(0x0),
	fCentValue(c.fCentValue),
	fPtAccCut(0x0),
	fEtaAccCut(0x0),
	fFakeSelection(c.fFakeSelection),
	fFake(c.fFake),
	fRejectIfNoQuark(c.fRejectIfNoQuark),	
	fMultiplicity(c.fMultiplicity),
	fConfiguration(c.fConfiguration)
{  
	//
	//copy constructor
	//
  fmcArray = new TClonesArray(*(c.fmcArray));
  fRecoCandidate = new AliAODRecoDecayHF(*(c.fRecoCandidate));
  fmcPartCandidate = new AliAODMCParticle(*(c.fmcPartCandidate));
	if (fProngs > 0){
		fLabelArray = new Int_t[fProngs];
                fPtAccCut = new Float_t[fProngs];
                fEtaAccCut = new Float_t[fProngs];
		if (c.fLabelArray) memcpy(fLabelArray,c.fLabelArray,fProngs*sizeof(Int_t));
		if (c.fPtAccCut) memcpy(fPtAccCut,c.fPtAccCut,fProngs*sizeof(Int_t));
		if (c.fEtaAccCut) memcpy(fEtaAccCut,c.fEtaAccCut,fProngs*sizeof(Int_t));
	}
}

//___________________________________________________________
void AliCFVertexingHF::SetDselection(UShort_t originDselection)
{
	// setting the way the D0 will be selected
	// 0 --> only from c quarks
	// 1 --> only from b quarks
	// 2 --> from both c quarks and b quarks
		
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
void AliCFVertexingHF::SetMCCandidateParam(Int_t label)
{
	//
	// setting the parameters (candidate and n. daughters)
	//	
	
        fmcPartCandidate = dynamic_cast <AliAODMCParticle*> (fmcArray->At(label));
	if (fmcPartCandidate){
	  fNDaughters = fmcPartCandidate->GetNDaughters();
	}
	else {
	  AliError(Form("Dynamic cast failed, fNdaughters will remain set to %d",fNDaughters));
	}
	return;
}

//____________________________________________________________
Int_t AliCFVertexingHF::MCcquarkCounting(AliAODMCParticle* mcPart) const
{
	//
	// counting the c-quarks
	// 

	Int_t cquarks = 0;
	if (mcPart) {
	  if (mcPart->GetPdgCode() == 4) cquarks++; 
	  if (mcPart->GetPdgCode() == -4) cquarks++; 
	}
	else {
		AliWarning("Particle not found in tree, skipping\n"); 
		return cquarks;
	} 
	
	return cquarks;
}

//________________________________________________________
Bool_t AliCFVertexingHF::CheckMCPartFamily(AliAODMCParticle */*mcPart*/, TClonesArray */*mcArray*/) const 
{
	// 
	//checking the family
	//

	Int_t pdgGranma = CheckOrigin();

	if (pdgGranma == -99999){
		AliDebug(2,"This particle does not have a quark in his genealogy\n");
		return kFALSE;
	}
	if (pdgGranma == -9999){
		AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");	
		return kFALSE;
	}	
	
	if (pdgGranma == -999){
		AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");  
		return kFALSE;
	}	
	
	if (!CheckMCDaughters()) {
	  AliDebug(3, "CheckMCDaughters false");
	  return kFALSE;
	}
	if (!CheckMCChannelDecay()) {
		AliDebug(3,"CheckMCChannelDecay false");
		return kFALSE;
	}
	return kTRUE;
}

//_________________________________________________________________________________________________
Int_t AliCFVertexingHF::CheckOrigin() const 
{		
	//
	// checking whether the mother of the particles come from a charm or a bottom quark
	//
	
	Int_t pdgGranma = 0;
	Int_t mother = 0;
	mother = fmcPartCandidate->GetMother();
	Int_t istep = 0;
	Int_t abspdgGranma =0;
	Bool_t isFromB=kFALSE;
	Bool_t isQuarkFound=kFALSE;
	while (mother >0 ){
		istep++;
		AliDebug(2,Form("mother at step %d = %d", istep, mother));
		AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother));
		if (mcGranma){
			pdgGranma = mcGranma->GetPdgCode();
			AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
			abspdgGranma = TMath::Abs(pdgGranma);
			if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
			  isFromB=kTRUE;
			}
			if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
			mother = mcGranma->GetMother();
		}else{
			AliError("Failed casting the mother particle!");
			break;
		}
	}

	if(fRejectIfNoQuark && !isQuarkFound) return -99999;
	if(isFromB){
	  if (!fKeepDfromB) return -9999; //skip particle if come from a B meson.
	}
	else{
	  if (fKeepDfromBOnly) return -999;
	}
	return pdgGranma;
}

//___________________________________________
Bool_t AliCFVertexingHF::CheckMCDaughters()const 
{
	//
	// checking the daughters
	// at MC level

	AliAODMCParticle *mcPartDaughter;
	Bool_t checkDaughters = kFALSE;
	
	Int_t label0 = fmcPartCandidate->GetDaughter(0);
	Int_t label1 = fmcPartCandidate->GetDaughter(1);
	AliDebug(3,Form("label0 = %d, label1 = %d",label0,label1));
	if (label1==0 || label0 == 0){
		AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
		return checkDaughters;  
	}
	
	if (fLabelArray == 0x0) {
		return checkDaughters;
	}  

	for (Int_t iProng = 0; iProng<fProngs; iProng++){
		mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[iProng]));    
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
	//
	// fill the container for Generator level selection
	//

	Bool_t mcContainerFilled = kFALSE;  
	
	Double_t* vectorMC = new Double_t[fNVar];
	for (Int_t iVar = 0; iVar<fNVar; iVar++) vectorMC[iVar]= 9999.;
	
	if (GetGeneratedValuesFromMCParticle(&vectorMC[0])){
		for (Int_t iVar = 0; iVar<fNVar; iVar++){			
			containerInputMC[iVar] = vectorMC[iVar];
		}		
		mcContainerFilled = kTRUE;		
	}
	delete [] vectorMC;
	vectorMC = 0x0;
	return mcContainerFilled;	
}

//______________________________________________________
Bool_t AliCFVertexingHF::FillRecoContainer(Double_t *containerInput) 
{  
	//	
	// fill the container for Reconstrucred level selection
	//

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
	
	delete [] vectorValues;
	delete [] vectorReco;
	vectorValues = 0x0;
	vectorReco = 0x0;
	return recoContainerFilled;	
}

//_____________________________________________________
Bool_t AliCFVertexingHF::MCAcceptanceStep() const
{
	//
	// checking the MC acceptance step
	//

	Bool_t bMCAccStep = kFALSE;
	
	AliAODMCParticle *mcPartDaughter;
	Int_t label0 = fmcPartCandidate->GetDaughter(0);
	Int_t label1 = fmcPartCandidate->GetDaughter(1);
	if (label1==0 || label0 == 0){
		AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
		return bMCAccStep;  
	}
	
	if (fLabelArray == 0x0) {
		return bMCAccStep;
	}  

	for (Int_t iProng = 0; iProng<fProngs; iProng++){
		mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fLabelArray[iProng]));    
		if (!mcPartDaughter) {
			AliWarning("At least one Daughter Particle not found in tree, skipping"); 
			return bMCAccStep;  
		}
		Double_t eta = mcPartDaughter->Eta();
		Double_t pt = mcPartDaughter->Pt();
		
		//set values of eta and pt in the constructor.
		//		if (TMath::Abs(eta) > 0.9 || pt < 0.1){
		if (TMath::Abs(eta) > fEtaAccCut[iProng] || pt < fPtAccCut[iProng]){
			AliDebug(3,Form("At least one daughter has eta or pt outside the required range (|eta| = %f, pt = %f, should be |eta| < %f, pt > %f \n", TMath::Abs(eta), pt, fEtaAccCut[iProng], fPtAccCut[iProng])); 
			return bMCAccStep;
		}
	}  
	bMCAccStep = kTRUE;
	return bMCAccStep; 
	
}
 //_____________________________________________________
Bool_t AliCFVertexingHF::MCRefitStep(AliAODEvent *aodEvent, AliESDtrackCuts **trackCuts) const
{		
	//
	// check on the kTPCrefit and kITSrefit conditions of the daughters
	//
	Bool_t bRefitStep = kFALSE;
	
	Int_t label0 = fmcPartCandidate->GetDaughter(0);
	Int_t label1 = fmcPartCandidate->GetDaughter(1);
	
	if (label1==0 || label0 == 0){
		AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
		return bRefitStep;  
	}
	
	if (fLabelArray == 0x0) {
		return bRefitStep;
	}  
	
	Int_t foundDaughters = 0;
	Int_t* temp = new Int_t[fProngs];
	for (Int_t ilabel = 0; ilabel<fProngs; ilabel++){
		temp[ilabel] = fLabelArray[ilabel];
	}

	//	if (trackCuts->GetRequireTPCRefit() || trackCuts->GetRequireITSRefit()){
		
	for(Int_t iaod =0; iaod<aodEvent->GetNumberOfTracks(); iaod++){
		AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iaod);
		if(track->GetStatus()&AliESDtrack::kITSpureSA) continue;
		Bool_t foundTrack = kFALSE;
		Int_t prongindex;
		for (Int_t ilabel = 0; ilabel<fProngs; ilabel++){
			AliDebug(3,Form("fLabelArray[%d] = %d, track->GetLabel() = %d",ilabel,fLabelArray[ilabel],TMath::Abs(track->GetLabel())));
			if ((track->GetLabel()<0)&&(fFakeSelection==1)) continue;
			if ((track->GetLabel()>0)&&(fFakeSelection==2)) continue;

			if (TMath::Abs(track->GetLabel()) == temp[ilabel]) {
				foundTrack = kTRUE;
				temp[ilabel] = 0;
				prongindex=ilabel;
				break;
			}
		}
		if (foundTrack){
			foundDaughters++;
			AliDebug(4,Form("daughter %d \n",foundDaughters));
			if(trackCuts[prongindex]->GetRequireTPCRefit()){
				if(track->GetStatus()&AliESDtrack::kTPCrefit) {
					bRefitStep = kTRUE;
				}
				else {
					AliDebug(3, "Refit cut not passed , missing TPC refit\n");
					delete [] temp;
					temp = 0x0;
					return kFALSE;
				}
			}
			
			if (trackCuts[prongindex]->GetRequireITSRefit()) {
				if(track->GetStatus()&AliESDtrack::kITSrefit){
					bRefitStep = kTRUE;
				}
				else {
					AliDebug(3, "Refit cut not passed , missing ITS refit\n");
					delete [] temp;
					temp = 0x0;
					return kFALSE;
				}
			}
		}      	
		if (foundDaughters == fProngs){				
			break;
		}			
	}    
	//} 						
	delete [] temp;
	temp = 0x0;
	if (foundDaughters== fProngs)  return bRefitStep;
	else return kFALSE;
}

//____________________________________________________________________________

Bool_t AliCFVertexingHF::RecoStep() 
{ 
	//
	//check also vertex and ITS Refit and TPC Refit
	//

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
	
	if (pdgGranma == -99999){
		AliDebug(2,"This particle does not have a quark in his genealogy\n");
		return bRecoStep;
	}
	if (pdgGranma == -9999){
		AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only prompt charm particles\n");
		return bRecoStep;
	}

	if (pdgGranma == -999){
		AliDebug(2,"This particle come from a  prompt charm particle but according to the settings of the task, we want only the ones coming from B\n");
		return bRecoStep;
	}
   
	bRecoStep=kTRUE;
	return bRecoStep;  
}	
//____________________________________________
Double_t AliCFVertexingHF::GetEtaProng(Int_t iProng) const 
{
	//
	// getting eta of the prong
	//
	
	if (fRecoCandidate){
		Double_t etaProng = fRecoCandidate->EtaProng(iProng);  
		return etaProng;
	}
	return 999999;  
}
//______________________________________________________
Double_t AliCFVertexingHF::GetPtProng(Int_t iProng) const 
{
	//
	// getting pt of the prong
	//

	if (fRecoCandidate){
		Double_t ptProng = fRecoCandidate->PtProng(iProng);  
		return ptProng;
	}
	return 999999;  
	
}

//____________________________________________________________________

Bool_t AliCFVertexingHF::RecoAcceptStep(AliESDtrackCuts **trackCuts) const
{
	//
	// reco Acceptance step
	//
	
	Bool_t bRecoAccStep = kFALSE;
	
	Float_t etaCutMin, ptCutMin, etaCutMax, ptCutMax;
	
	Float_t etaProng=0., ptProng=0.; 
	
	for (Int_t iProng =0; iProng<fProngs; iProng++){
		
		trackCuts[iProng]->GetEtaRange(etaCutMin, etaCutMax);
		trackCuts[iProng]->GetPtRange(ptCutMin, ptCutMax);
		etaProng = GetEtaProng(iProng);
		ptProng = GetPtProng(iProng);
		
		Bool_t acceptanceProng = (etaProng>etaCutMin && etaProng<etaCutMax && ptProng>ptCutMin && ptProng<ptCutMax);
		if (!acceptanceProng) {
			AliDebug(2,"At least one reconstructed prong isn't in the acceptance\n");
			return bRecoAccStep;
		}
	}
	
	bRecoAccStep=kTRUE;
	return bRecoAccStep;
}
//___________________________________________________________

Bool_t AliCFVertexingHF::FillUnfoldingMatrix(Double_t fill[4]) const
{
	//
	// filling the unfolding matrix
	//
	
	if(fmcPartCandidate){
		
		fill[0] = GetPtCand();
		fill[1] = GetYCand();
		
		fill[2] =  fmcPartCandidate->Pt(); 
		fill[3] =  fmcPartCandidate->Y(); 
		
		return kTRUE;
	}
	
	return kFALSE;
}
//___________________________________________________________

Int_t AliCFVertexingHF::CheckReflexion(Char_t isSign)
{
	//
	// check for reflexion (particle/antiparticle)
	//

	Int_t mcLabel = GetMCLabel();
	
	if (mcLabel == -1) {
		AliDebug(2,"No MC particle found");
		return 0;
	}
	else{
		fmcPartCandidate = (AliAODMCParticle*)fmcArray->At(mcLabel);
		if (!fmcPartCandidate){
			AliWarning("Could not find associated MC in AOD MC tree");
			return 0;
		}    
	}
	
	if(fmcPartCandidate->GetPdgCode()>0) {
		if (isSign == 1){ // I ask for antiparticle only
			AliDebug(2,"candidate is particle, I ask for antiparticle only");
			return 0;
		}
		return 1;  // particle
	}
	else if(fmcPartCandidate->GetPdgCode()<0) {
		if (isSign == 0){ // I ask for particle only
			AliDebug(2,"candidate is antiparticle, I ask for particle only");
			return 0;
		}
		return 2;  // antiparticle
	}
	else return 0;  // ....shouldn't be...

}
//___________________________________________________________

Bool_t AliCFVertexingHF::SetLabelArray()
{
	//
	// setting the label arrays
	//

	Bool_t bLabelArray = kFALSE;

	fLabelArray = new Int_t[fProngs];

	AliAODMCParticle *mcPartDaughter;
	Int_t label0 = fmcPartCandidate->GetDaughter(0);
	Int_t label1 = fmcPartCandidate->GetDaughter(1);
	AliDebug(2,Form("label0 = %d, label1 = %d",label0,label1));
	if (label1==0 || label0 == 0){
		AliDebug(2, Form("The MC particle doesn't have correct daughters, skipping!!"));
		delete [] fLabelArray; 
		fLabelArray = 0x0;  
		return bLabelArray;
	}
	
	if (label1 - label0 == fProngs-1){
		for (Int_t iProng = 0; iProng<fProngs; iProng++){
			mcPartDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(label0+iProng));    
			if (mcPartDaughter){
				fLabelArray[iProng] =  mcPartDaughter->GetLabel();
			}
			else{
				AliError("Failed casting the daughter particle, returning a NULL label array");
				delete [] fLabelArray; 
				fLabelArray = 0x0;  
				return bLabelArray;
			}
		}

	}
	// resonant decay channel
	else if (label1 - label0 == fProngs-2 && fProngs > 2){
		Int_t labelFirstDau = fmcPartCandidate->GetDaughter(0);
		Int_t foundDaughters = 0;
		for(Int_t iDau=0; iDau<fProngs-1; iDau++){
			Int_t iLabelDau = labelFirstDau+iDau;
			AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iLabelDau));
                        if ( ! part ) {
                          AliError("Wrong particle type in fmcArray");
			  delete [] fLabelArray; 
			  fLabelArray = 0x0; 
                          return bLabelArray;
                        }  
			Int_t pdgCode=TMath::Abs(part->GetPdgCode());
			if(pdgCode==211 || pdgCode==321 || pdgCode==2212){
				if (part) {
					fLabelArray[foundDaughters] = part->GetLabel();
					foundDaughters++;
				}
				else{
					AliError("Error while casting particle! returning a NULL array");
					delete [] fLabelArray; 
					fLabelArray = 0x0;  
					return bLabelArray;
				}
			}
			else{
				Int_t nDauRes=part->GetNDaughters();
				if(nDauRes!=2) {
					delete [] fLabelArray; 
					fLabelArray = 0x0;  
					return bLabelArray;
				}
				Int_t labelFirstDauRes = part->GetDaughter(0); 
				for(Int_t iDauRes=0; iDauRes<nDauRes; iDauRes++){
					Int_t iLabelDauRes = labelFirstDauRes+iDauRes;
					AliAODMCParticle* dauRes = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iLabelDauRes));
					if (dauRes){
						fLabelArray[foundDaughters] = dauRes->GetLabel();
						foundDaughters++;
					}
					else{
						AliError("Error while casting resonant daughter! returning a NULL array");
						delete [] fLabelArray; 
						fLabelArray = 0x0;  
						return bLabelArray;
					}
				}
			}
		}
		if (foundDaughters != fProngs){
			delete [] fLabelArray; 
			fLabelArray = 0x0;  
			return bLabelArray;
		}
	}
	// wrong correspondance label <--> prongs
	else{
		delete [] fLabelArray; 
		fLabelArray = 0x0;  
		return bLabelArray;
	}
	SetAccCut();   // setting the pt and eta acceptance cuts
	bLabelArray = kTRUE;
	return bLabelArray;
}

//___________________________________________________________

void AliCFVertexingHF::SetPtAccCut(Float_t* ptAccCut)
{
	//
	// setting the pt cut to be used in the Acceptance steps (MC+Reco)
	//

	if (fProngs>0){
		for (Int_t iP=0; iP<fProngs; iP++){
			fPtAccCut[iP]=ptAccCut[iP];
		}
	}
	return;
}		



//___________________________________________________________

void AliCFVertexingHF::SetEtaAccCut(Float_t* etaAccCut)
{
	//
	// setting the eta cut to be used in the Acceptance steps (MC+Reco)
	//

	if (fProngs>0){
		for (Int_t iP=0; iP<fProngs; iP++){
			fEtaAccCut[iP]=etaAccCut[iP];
		}
	}
	return;
}	
//___________________________________________________________

void AliCFVertexingHF::SetAccCut(Float_t* ptAccCut, Float_t* etaAccCut)
{
	//
	// setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
	//

	if (fProngs>0){
		for (Int_t iP=0; iP<fProngs; iP++){
			fPtAccCut[iP]=ptAccCut[iP];
			fEtaAccCut[iP]=etaAccCut[iP];
		}
	}
	return;
}		

//___________________________________________________________

void AliCFVertexingHF::SetAccCut()
{
	//
	// setting the pt and eta cut to be used in the Acceptance steps (MC+Reco)
	//

	if (fProngs>0){
		for (Int_t iP=0; iP<fProngs; iP++){
			fPtAccCut[iP]=0.1;
			fEtaAccCut[iP]=0.9;
		}
	}
	return;
}		
