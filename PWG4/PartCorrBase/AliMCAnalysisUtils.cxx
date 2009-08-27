/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id: AliMCAnalysisUtils.cxx 21839 2007-10-29 13:49:42Z gustavo $ */

//_________________________________________________________________________
// Class for analysis utils for MC data
// stored in stack or event header.
// Contains:
//  - method to check the origin of a given track/cluster
//  - method to obtain the generated jets
//                
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////
  

// --- ROOT system ---
#include <TMath.h>
#include <TList.h>
#include "TParticle.h"

//---- ANALYSIS system ----
#include "AliMCAnalysisUtils.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"

  ClassImp(AliMCAnalysisUtils)

 //________________________________________________
  AliMCAnalysisUtils::AliMCAnalysisUtils() : 
    TObject(), fCurrentEvent(-1), fDebug(-1), 
    fJetsList(new TList), fMCGenerator("PYTHIA")
{
  //Ctor
}

//____________________________________________________________________________
AliMCAnalysisUtils::AliMCAnalysisUtils(const AliMCAnalysisUtils & mcutils) :   
  TObject(mcutils), fCurrentEvent(mcutils.fCurrentEvent), fDebug(mcutils.fDebug),
  fJetsList(mcutils.fJetsList), fMCGenerator(mcutils.fMCGenerator)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliMCAnalysisUtils & AliMCAnalysisUtils::operator = (const AliMCAnalysisUtils & mcutils)
{
  // assignment operator
  
  if(&mcutils == this) return *this;
  fCurrentEvent = mcutils.fCurrentEvent ;
  fDebug        = mcutils.fDebug;
  fJetsList     = mcutils.fJetsList;
  fMCGenerator  = mcutils.fMCGenerator;
  
  return *this; 
}

//____________________________________________________________________________
AliMCAnalysisUtils::~AliMCAnalysisUtils() 
{
  // Remove all pointers.
  
  if (fJetsList) {
    fJetsList->Clear();
    delete fJetsList ;
  }     
}


//_________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOrigin(const Int_t label, AliCaloTrackReader* reader, const Int_t input = 0) {
	//Play with the montecarlo particles if available
	Int_t tag = 0;
	
	//Select where the information is, ESD-galice stack or AOD mcparticles branch
	if(reader->ReadStack()){
		tag = CheckOriginInStack(label, reader->GetStack());
	}
	else if(reader->ReadAODMCParticles()){
		tag = CheckOriginInAOD(label, reader->GetAODMCParticles(input));
	}
	
	return tag ;
}	

//_________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInStack(const Int_t label, AliStack* stack) {
  // Play with the MC stack if available. Tag particles depending on their origin.
  // Do same things as in CheckOriginInAOD but different input.
	
  if(!stack) {
    printf("AliMCAnalysisUtils::CheckOriginInStack() - Stack is not available, check analysis settings in configuration file, STOP!!\n");
    abort();
  }
  //  printf("label %d, ntrack %d, nprim %d\n",label, stack->GetNtrack(), stack->GetNprimary());
  //   for(Int_t i = 0; i< stack->GetNprimary(); i++){
  //      TParticle *particle =   stack->Particle(i);
  // 			//particle->Print();
  //   }
	

  Int_t tag = 0;
  if(label >= 0 && label <  stack->GetNtrack()){
    //Mother
    TParticle * mom = stack->Particle(label);
    Int_t mPdg = TMath::Abs(mom->GetPdgCode());
    Int_t mStatus =  mom->GetStatusCode() ;
    Int_t iParent =  mom->GetFirstMother() ;
    if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
    
    //GrandParent
    TParticle * parent = new TParticle ;
    Int_t pPdg = -1;
    Int_t pStatus =-1;
    if(iParent > 0){
      parent = stack->Particle(iParent);
      pPdg = TMath::Abs(parent->GetPdgCode());
      pStatus = parent->GetStatusCode();  
    }
    else if(fDebug > 0 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Parent with label %d\n",iParent);

	//Check if mother is converted, if not, get the first non converted mother
    if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && mStatus == 0){
		SetTagBit(tag,kMCConversion);
		//Check if the mother is photon or electron with status not stable
		while ((pPdg == 22 || pPdg == 11) && mStatus != 1) {
			//Mother
			mom = stack->Particle(mom->GetFirstMother());
			mPdg = TMath::Abs(mom->GetPdgCode());
			mStatus =  mom->GetStatusCode() ;
			iParent =  mom->GetFirstMother() ;
			if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
			
			//GrandParent
			if(iParent > 0){
				parent = stack->Particle(iParent);
				pPdg = TMath::Abs(parent->GetPdgCode());
				pStatus = parent->GetStatusCode();  
			}
		}//while	  
    }//mother and parent are electron or photon and have status 0
	else if(mStatus == 0){	
		//Still a conversion but only one electron/photon generated. Just from hadrons
		if(pPdg == 2112 ||  pPdg == 211 ||  
		   pPdg == 321 ||  pPdg == 2212  ||  pPdg == 130  ||  pPdg == 13 ) 
			SetTagBit(tag,kMCConversion);
			mom = stack->Particle(mom->GetFirstMother());
			mPdg = TMath::Abs(mom->GetPdgCode());
		//Comment for the next lines, we do not check the parent of the hadron for the moment.
		//iParent =  mom->GetFirstMother() ;
		//if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
		
		//GrandParent
		//if(iParent >= 0){
		//	parent = stack->Particle(iParent);
		//	pPdg = TMath::Abs(parent->GetPdgCode());
		//}
	}  	  
	// conversion fo electrons/photons checked  	  
	  
    //first check for typical charged particles
    if(mPdg == 13) SetTagBit(tag,kMCMuon);
    else if(mPdg == 211) SetTagBit(tag,kMCPion);
    else if(mPdg == 321) SetTagBit(tag,kMCKaon);
    else if(mPdg == 2212) SetTagBit(tag,kMCProton);
    //check for pi0 and eta (shouldn't happen unless their decays were
    //turned off)
    else if(mPdg == 111)  SetTagBit(tag,kMCPi0);
    else if(mPdg == 221)  SetTagBit(tag,kMCEta);    
	//Photons  
    else if(mPdg == 22){
      SetTagBit(tag,kMCPhoton);
      if(mStatus == 1){ //undecayed particle
		  if(fMCGenerator == "PYTHIA"){
			  if(iParent < 8 && iParent > 5) {//outgoing partons
				  if(pPdg == 22) SetTagBit(tag,kMCPrompt);
				  else SetTagBit(tag,kMCFragmentation);
			  }//Outgoing partons 
			  else  if(iParent <= 5) {
				  SetTagBit(tag, kMCISR); //Initial state radiation
			  }
			  else if(pStatus == 11){//Decay
				  if(pPdg == 111) SetTagBit(tag,kMCPi0Decay);
				  else if (pPdg == 221) SetTagBit(tag, kMCEtaDecay);
				  else SetTagBit(tag,kMCOtherDecay);
			  }//Decay
			  else {
				  printf("AliMCAnalysisUtils::ChecOrigingInAOD() - what is it? Mother mPdg %d, status %d \n    Parent  iParent %d, pPdg %d %s, status %d\n",
						 mPdg, mStatus,iParent, pPdg, parent->GetName(),pStatus);
				  SetTagBit(tag,kMCOtherDecay);//Check
			  }
		  }//PYTHIA

		  else if(fMCGenerator == "HERWIG"){	  
			  if(pStatus < 197){//Not decay
				  while(1){
					  if(parent->GetFirstMother()<=5) break;
					  iParent = parent->GetFirstMother();
					  parent=stack->Particle(iParent);
					  pStatus= parent->GetStatusCode();
					  pPdg = parent->GetPdgCode();
				  }//Look for the parton
	    
				  if(iParent < 8 && iParent > 5) {
					  if(pPdg == 22) SetTagBit(tag,kMCPrompt);
					  else SetTagBit(tag,kMCFragmentation);
				  }
				  else SetTagBit(tag,kMCISR);//Initial state radiation
			  }//Not decay
			  else{//Decay
				  if(pPdg == 111) SetTagBit(tag,kMCPi0Decay); 
				  else if (pPdg == 221) SetTagBit(tag,kMCEtaDecay);
				  else SetTagBit(tag,kMCOtherDecay);
			  }//Decay
		  }//HERWIG
		  
		  else SetTagBit(tag,kMCUnknown);
		  
      }//Status 1 : created by event generator
		
      else if(mStatus == 0){ // geant
		  if(pPdg == 111) SetTagBit(tag,kMCPi0Decay); 
		  else if (pPdg == 221) SetTagBit(tag,kMCEtaDecay);
		  else  SetTagBit(tag,kMCOtherDecay);	
      }//status 0 : geant generated
		
    }//Mother Photon

    //Electron check.  Where did that electron come from?
    else if(mPdg == 11){ //electron
      SetTagBit(tag,kMCElectron);	
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInStack() - Checking ancestors of electrons");
      
      //check first for B and C ancestry, then other possibilities.
      //An electron from a photon parent could have other particles in
      //its history and we would want to know that, right?
      
      if(mStatus == 1) { //electron from event generator
        if (pPdg == 23) SetTagBit(tag,kMCZDecay);         //parent is Z-boson                                                     
        else if (pPdg == 24) SetTagBit(tag,kMCWDecay);    //parent is W-boson                                                     
        else if (pPdg == 111) SetTagBit(tag,kMCPi0Decay); //Pi0 Dalitz decay                                                      
        else if (pPdg == 221) SetTagBit(tag,kMCEtaDecay); //Eta Dalitz decay                                                      
        else { //check the electron's ancestors for B/C contribution                                                     
          Bool_t bAncestor = kFALSE;
          Bool_t cAncestor = kFALSE;
          TParticle * ancestors = stack->Particle(label); //start from electron                                          
          Int_t iAncestors = ancestors->GetFirstMother(); //get its mother index                                         
          if(iAncestors < 0) return tag; //no ancestor, shouldn't happen                                         
          ancestors = stack->Particle(iAncestors); //get electron mother                                                 
          Int_t aPdg = TMath::Abs(ancestors->GetPdgCode()); //get electron mother pdg                                    
          if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOrigin: Scaning the decay chain for bottom/charm generated electron");
          while(1){//check electron's mother's pdg and its mother(s)
		   //pdgs                                                
            if((499 < aPdg && aPdg < 600)||(4999 < aPdg && aPdg < 6000)) bAncestor = kTRUE;
            if((399 < aPdg && aPdg < 500)||(3999 < aPdg && aPdg < 5000)) cAncestor = kTRUE;
            if(bAncestor && cAncestor) break;
            iAncestors = ancestors->GetFirstMother();
            if(iAncestors < 0) break;
            ancestors = stack->Particle(iAncestors);
            aPdg = ancestors->GetPdgCode();
          }//searching for ancestors                                                                                     
          if(bAncestor && cAncestor) SetTagBit(tag,kMCEFromCFromB);//Decay chain has both B and C                                
          else if(bAncestor && !cAncestor) SetTagBit(tag,kMCEFromB);//Decay chain has only B                                     
	  else if(!bAncestor && cAncestor) SetTagBit(tag,kMCEFromC);//Decay chain has only C                                     
        }
        //if it is not from any of the above, where is it from?                                                          
        if(fDebug > 0) printf("Electron from other origin: %s (pdg = %d)\n",parent->GetName(),pPdg);
        SetTagBit(tag,kMCOtherDecay);
      } 
      else if (mStatus == 0) { //electron from GEANT
	
	//Rewind ancestry and check for electron with status == 1
	//if we find one, we'll assume that this object is from an
	//electron but that it may have gone through some showering in
	//material before the detector
	
	//Not a double-counting problem because we are only accessing
	//these histories for MC labels connected to a reco object.
	//If you wanted to use this to sort through the kine stack
	//directly, might it be a problem?
	Bool_t eleFromEvGen = kFALSE;
        Bool_t bAncestor = kFALSE;
        Bool_t cAncestor = kFALSE;

        TParticle * ancestors = stack->Particle(label); //start from electron                                            
        Int_t iAncestors = ancestors->GetFirstMother(); //get its mother index                                           
        if(iAncestors < 0) return tag; //no ancestor, shouldn't happen                                           
        ancestors = stack->Particle(iAncestors); //get electron mother                                                   
        Int_t aStatus = ancestors->GetStatusCode(); //get electron mother's status                                       
        Int_t aPdg = TMath::Abs(ancestors->GetPdgCode()); //get electron mother pdg                                      
        if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOrigin: Scanning the decay chain for bottom/charm generated electron");
        while(1){//check electron's mother's pdg and its mother(s) pdgs
          if(aStatus == 1 && aPdg == 11) eleFromEvGen = kTRUE;
          if(eleFromEvGen && aPdg == 23) SetTagBit(tag,kMCZDecay);
          if(eleFromEvGen && aPdg == 24) SetTagBit(tag,kMCWDecay);
          if(eleFromEvGen && ((499 < aPdg && aPdg < 600)||(4999 < aPdg && aPdg < 6000))) bAncestor = kTRUE;
          if(eleFromEvGen && ((399 < aPdg && aPdg < 500)||(3999 < aPdg && aPdg < 5000))) cAncestor = kTRUE;
          if(bAncestor && cAncestor) break;
          iAncestors = ancestors->GetFirstMother();
          if(iAncestors <= 5) break;
          ancestors = stack->Particle(iAncestors);
          aPdg = ancestors->GetPdgCode();
        }//searching for ancestors                                   
        if(bAncestor && cAncestor) SetTagBit(tag,kMCEFromCFromB);//Decay chain has both B and C                                  
        else if(bAncestor && !cAncestor) SetTagBit(tag,kMCEFromB);//Decay chain has only B                                       
      else if(!bAncestor && cAncestor) SetTagBit(tag,kMCEFromC);//Decay chain has only C                                       
        if(pPdg == 22 || pPdg == 11 || pPdg == 2112 || pPdg == 211 ||
           pPdg == 321 || pPdg == 2212 || pPdg == 130 || pPdg == 13)
          SetTagBit(tag,kMCConversion);
        else if (pPdg == 111)  SetTagBit(tag,kMCPi0Decay);
	else if (pPdg == 221)    SetTagBit(tag,kMCEtaDecay);
        else {
          //if it is not from any of the above, where is it from?                                                        
          if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
          if(fDebug > 0) printf("Electron from other origin: %s (pdg = %d)\n",parent->GetName(),pPdg);
          SetTagBit(tag,kMCOtherDecay);
        }
      } //GEANT check
    }//electron check
    //Cluster was made by something else
    else {
       if(fDebug > 0) printf("\tSetting kMCUnknown for cluster from %s (pdg = %d, Parent pdg = %d)\n",mom->GetName(),mPdg,pPdg);
      SetTagBit(tag,kMCUnknown);
    }
  }//Good label value
  else{
    if(label < 0 ) printf("AliMCAnalysisUtils::CheckOriginInStack() *** bad label or no stack ***:  label %d \n", label);
    if(label >=  stack->GetNtrack()) printf("AliMCAnalysisUtils::CheckOriginInStack() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
    SetTagBit(tag,kMCUnknown);
  }//Bad label
  
  return tag;
  
}


//_________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInAOD(const Int_t label, TClonesArray *mcparticles) {
  // Play with the MCParticles in AOD if available. Tag particles depending on their origin.
  // Do same things as in CheckOriginInStack but different input.
  if(!mcparticles) {
    printf("AliMCAnalysisUtils::CheckOriginInAOD() - AODMCParticles is not available, check analysis settings in configuration file!!\n");
    
  }
  
  //	printf(">>>>>>>>> Second Event, AODMCParticles, nMC %d\n", (GetAODMCParticles(1))->GetEntriesFast());
  //	for(Int_t i = 0; i< (GetAODMCParticles(1))->GetEntriesFast(); i++){
  //			AliAODMCParticle *particle = (AliAODMCParticle*) (GetAODMCParticles(1))->At(i);
  //			printf("** index %d, PDG %d, Flag %d, Mother %d, Daughter %d and %d, E %2.3f, pT %2.3f, phi %2.3f, eta %2.3f\n", 
  //				   i, particle->GetPdgCode(), particle->GetFlag(), particle->GetMother(), particle->GetDaughter(0),particle->GetDaughter(1),
  //				   particle->E(),particle->Pt(), particle->Phi()*TMath::RadToDeg(), particle->Eta());
  //		if(particle->IsPhysicalPrimary()) printf("Is Physical Primary !\n");
  //		if(!particle->IsPrimary()) printf("Not Primary !\n");
  //   }
  
	
  Int_t tag = 0;
  Int_t nprimaries = mcparticles->GetEntriesFast();
  if(label >= 0 && label < nprimaries){
    //Mother
    AliAODMCParticle * mom = (AliAODMCParticle *) mcparticles->At(label);
    Int_t mPdg = TMath::Abs(mom->GetPdgCode());
    Int_t iParent =  mom->GetMother() ;
    if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
    
    //GrandParent
    AliAODMCParticle * parent = new AliAODMCParticle ;
    Int_t pPdg = -1;
    if(iParent >= 0){
      parent = (AliAODMCParticle *) mcparticles->At(iParent);
      pPdg = TMath::Abs(parent->GetPdgCode());
    }
    else if(fDebug > 0 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Parent with label %d\n",iParent);
 
	//Check if mother is converted, if not, get the first non converted mother
	if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && !mom->IsPrimary()){
		  SetTagBit(tag,kMCConversion);
		  //Check if the mother is photon or electron with status not stable
		  while ((pPdg == 22 || pPdg == 11) && !mom->IsPhysicalPrimary()) {
			  //Mother
			  mom = (AliAODMCParticle *) mcparticles->At(mom->GetMother());
			  mPdg = TMath::Abs(mom->GetPdgCode());
			  iParent =  mom->GetMother() ;
			  if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
			  
			  //GrandParent
			  if(iParent >= 0){
				  parent = (AliAODMCParticle *) mcparticles->At(iParent);
				  pPdg = TMath::Abs(parent->GetPdgCode());
			  }
		  }//while	  
	}//mother and parent are electron or photon and have status 0 and parent is photon or electron
	else if(!mom->IsPrimary()){	
		//Still a conversion but only one electron/photon generated. Just from hadrons
		if(pPdg == 2112 ||  pPdg == 211 ||  
		   pPdg == 321 ||  pPdg == 2212  ||  pPdg == 130  ||  pPdg == 13 ) 
			SetTagBit(tag,kMCConversion);
			mom = (AliAODMCParticle *) mcparticles->At(mom->GetMother());
			mPdg = TMath::Abs(mom->GetPdgCode());
			//Comment for next lines, we do not check the parent of the hadron for the moment.
			//iParent =  mom->GetMother() ;
			//if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
		
			//GrandParent
			//if(iParent >= 0){
			//	parent = (AliAODMCParticle *) mcparticles->At(iParent);
			//	pPdg = TMath::Abs(parent->GetPdgCode());
			//}
	}  
	  
	// conversion into electrons/photons checked  
	  
    //first check for typical charged particles                                                                          
    if(mPdg == 13) SetTagBit(tag,kMCMuon);
    else if(mPdg == 211) SetTagBit(tag,kMCPion);
    else if(mPdg == 321) SetTagBit(tag,kMCKaon);
    else if(mPdg == 2212) SetTagBit(tag,kMCProton);
    //check for pi0 and eta (shouldn't happen unless their decays were                                                   
    //turned off)                                                                                                        
    else if(mPdg == 111)  SetTagBit(tag,kMCPi0);
    else if(mPdg == 221)  SetTagBit(tag,kMCEta);    
	//Photons  
    else if(mPdg == 22){
      SetTagBit(tag,kMCPhoton);
      if(mom->IsPhysicalPrimary()){ //undecayed particle
		  if(iParent < 8 && iParent > 5) {//outgoing partons
			  if(pPdg == 22) SetTagBit(tag,kMCPrompt);
			  else SetTagBit(tag,kMCFragmentation);
		  }//Outgoing partons
		  else if(iParent <= 5) {
			  SetTagBit(tag, kMCISR); //Initial state radiation
		  }
		  else if(parent->IsPrimary() && !parent->IsPhysicalPrimary()){//Decay
			  if(pPdg == 111) SetTagBit(tag,kMCPi0Decay);
			  else if (pPdg == 221) SetTagBit(tag, kMCEtaDecay);
			  else SetTagBit(tag,kMCOtherDecay);
		  }//Decay
		  else {
			  printf("AliMCAnalysisUtils::ChecOrigingInAOD() - what is it? Mother mPdg %d, is primary? %d, is physical %d \n    Parent  iParent %d, pPdg %d, is primary? %d, is physical? %d\n",
					  mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary(),iParent, pPdg,parent->IsPrimary(), parent->IsPhysicalPrimary());
			  SetTagBit(tag,kMCOtherDecay);//Check
		  }
		}// Pythia generated
		else if(!mom->IsPrimary()){	//Decays
		  if(pPdg == 111) SetTagBit(tag,kMCPi0Decay); 
		  else if (pPdg == 221) SetTagBit(tag,kMCEtaDecay);
		  else  SetTagBit(tag,kMCOtherDecay);
      }//not primary : geant generated, decays
      else  {
		  //printf("UNKNOWN 1, mom  pdg %d, primary %d, physical primary %d; parent %d, pdg %d, primary %d, physical primary %d \n",
		  //mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary(), iParent, pPdg, parent->IsPrimary(), parent->IsPhysicalPrimary());
		  SetTagBit(tag,kMCUnknown);
      }
    }//Mother Photon
    
    //Electron check.  Where did that electron come from?
    else if(mPdg == 11){ //electron
      SetTagBit(tag,kMCElectron);	
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Checking ancestors of electrons");
      
      //check first for B and C ancestry, then other possibilities.
      //An electron from a photon parent could have other particles in
      //its history and we would want to know that, right?
      Bool_t eleFromEvGen = kFALSE;
      Bool_t bAncestor = kFALSE;
      Bool_t cAncestor = kFALSE;
      
      if(mom->IsPhysicalPrimary()) { //electron from event generator
		  if      (pPdg == 23) SetTagBit(tag,kMCZDecay);  //parent is Z-boson
		  else if (pPdg == 24) SetTagBit(tag,kMCWDecay);  //parent is W-boson
		  else if (pPdg == 111) SetTagBit(tag,kMCPi0Decay); //Pi0 Dalitz decay
		  else if (pPdg == 221) SetTagBit(tag,kMCEtaDecay); //Eta Dalitz decay
		  else { //check the electron's ancestors for B/C contribution
			  bAncestor = kFALSE;
			  cAncestor = kFALSE;
			  AliAODMCParticle * ancestors = (AliAODMCParticle*)mcparticles->At(label); //start from electron
			  Int_t iAncestors = ancestors->GetMother(); //get its mother index
			  if(iAncestors < 0) return tag; //no ancestor, shouldn't happen
			  ancestors = (AliAODMCParticle *) mcparticles->At(iAncestors);
			  Int_t aPdg = TMath::Abs(ancestors->GetPdgCode());
			  if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOrigin: Scanning the decay chain for bottom/charm generated electron");
			  while(1){//check electron's mother's pdg and its mother(s) pdgs
				  if((499 < aPdg && aPdg < 600)||(4999 < aPdg && aPdg < 6000)) bAncestor = kTRUE;
				  if((399 < aPdg && aPdg < 500)||(3999 < aPdg && aPdg < 5000)) cAncestor = kTRUE;
				  if(bAncestor && cAncestor) break;
				  iAncestors = ancestors->GetMother();
				  if(iAncestors < 0) break;
				  ancestors = (AliAODMCParticle *) mcparticles->At(iAncestors);
				  aPdg = ancestors->GetPdgCode();
			  }//searching for ancestors
			  if(bAncestor && cAncestor) SetTagBit(tag,kMCEFromCFromB);//Decay chain has both B and C
			  else if(bAncestor && !cAncestor) SetTagBit(tag,kMCEFromB);//Decay chain has only B
			  else if(!bAncestor && cAncestor) SetTagBit(tag,kMCEFromC);//Decay chain has only C
		  }
		  //if it is not from any of the above, where is it from?
		   if(fDebug > 0) printf("Electron from other origin: %s (pdg = %d)\n",parent->GetName(),pPdg);
		  SetTagBit(tag,kMCOtherDecay);
      } 
      else if (!mom->IsPrimary()) { //electron from GEANT
	
		  //Rewind ancestry and check for electron with status == 1
		  //if we find one, we'll assume that this object is from an
		  //electron but that it may have gone through some showering in
		  //material before the detector
	
		  //Not a double-counting problem because we are only accessing
		  //these histories for MC labels connected to a reco object.
		  //If you wanted to use this to sort through the kine stack
		  //directly, might it be a problem?
		  eleFromEvGen = kFALSE;
		  bAncestor = kFALSE;
		  cAncestor = kFALSE;
	
		  AliAODMCParticle * ancestors = (AliAODMCParticle *) mcparticles->At(label);
		  Int_t iAncestors = ancestors->GetMother();
		  if(iAncestors < 0) return tag; //no ancestor, shouldn't happen
		  ancestors = (AliAODMCParticle *) mcparticles->At(iAncestors); //get electron mother
		  Int_t aPdg = TMath::Abs(ancestors->GetPdgCode()); //get electron mother pdg
		  if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Scanning the decay chain for bottom/charm electrons");

		  while(1){//check electron's mother's pdg and its mother(s) pdgs
			  if(ancestors->IsPhysicalPrimary() && aPdg == 11) eleFromEvGen = kTRUE;
			  if(eleFromEvGen && aPdg == 23) SetTagBit(tag,kMCZDecay);
			  if(eleFromEvGen && aPdg == 24) SetTagBit(tag,kMCWDecay);
			  if(eleFromEvGen && ((499 < aPdg && aPdg < 600)||(4999 < aPdg && aPdg < 6000))) bAncestor = kTRUE;
			  if(eleFromEvGen && ((399 < aPdg && aPdg < 500)||(3999 < aPdg && aPdg < 5000))) cAncestor = kTRUE;
			  if(bAncestor && cAncestor) break;
			  iAncestors = ancestors->GetMother();
			  if(iAncestors <= 5) break; //don't go all the way back
			  ancestors =  (AliAODMCParticle *) mcparticles->At(iAncestors);
			  aPdg = ancestors->GetPdgCode();
		  }//searching for ancestors                                                                                       
		  if(bAncestor && cAncestor) SetTagBit(tag,kMCEFromCFromB);//Decay chain has both B and C
		  else if(bAncestor && !cAncestor) SetTagBit(tag,kMCEFromB);//Decay chain has only B
		  else if(!bAncestor && cAncestor) SetTagBit(tag,kMCEFromC);//Decay chain has only C
		  if(pPdg == 22 || pPdg == 11 || pPdg == 2112 || pPdg == 211 ||
			 pPdg == 321 || pPdg == 2212 || pPdg == 130 || pPdg == 13)
			  SetTagBit(tag,kMCConversion);
		  else if (pPdg == 111)  SetTagBit(tag,kMCPi0Decay);
		  else if (pPdg == 221)  SetTagBit(tag,kMCEtaDecay);
		  else {
          //if it is not from any of the above, where is it from?
			  if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
			   if(fDebug > 0) printf("Electron from other origin: %s (pdg = %d)\n",parent->GetName(),pPdg);
			  SetTagBit(tag,kMCOtherDecay);
		  }
		} //GEANT check
	}//electron check
    //cluster was made by something else
    else {
       if(fDebug > 0) printf("\tSetting kMCUnknown for cluster from %s (pdg = %d, Parent pdg = %d)\n",mom->GetName(),mPdg,pPdg);
      SetTagBit(tag,kMCUnknown);
    }
  }//Good label value
  else{
    if(label < 0 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() *** bad label or no stack ***:  label %d \n", label);
    if(label >=  mcparticles->GetEntriesFast()) printf("AliMCAnalysisUtils::CheckOriginInAOD() *** large label ***:  label %d, n tracks %d \n", label, mcparticles->GetEntriesFast());
    SetTagBit(tag,kMCUnknown);
  }//Bad label
  
  return tag;
  
}



//_________________________________________________________________________
TList * AliMCAnalysisUtils::GetJets(AliCaloTrackReader * reader){
 //Return list of jets (TParticles) and index of most likely parton that originated it.
  AliStack * stack = reader->GetStack();
  Int_t iEvent = reader->GetEventNumber();	
  AliGenEventHeader * geh = reader->GetGenEventHeader();
  if(fCurrentEvent!=iEvent){
    fCurrentEvent = iEvent;
    fJetsList = new TList;
    Int_t nTriggerJets = 0;
    Float_t tmpjet[]={0,0,0,0};
		
    //printf("Event %d %d\n",fCurrentEvent,iEvent);
    //Get outgoing partons
    if(stack->GetNtrack() < 8) return fJetsList;
    TParticle * parton1 =  stack->Particle(6);
    TParticle * parton2 =  stack->Particle(7);
    if(fDebug > 2){
      printf("AliMCAnalysisUtils::GetJets() - parton 6 : %s, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
	     parton1->GetName(),parton1->Pt(),parton1->Energy(),parton1->Phi()*TMath::RadToDeg(),parton1->Eta());
      printf("AliMCAnalysisUtils::GetJets() - parton 7 : %s, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
	     parton2->GetName(),parton2->Pt(),parton2->Energy(),parton2->Phi()*TMath::RadToDeg(),parton2->Eta());
		}
// 		//Trace the jet from the mother parton
// 		Float_t pt  = 0;
// 		Float_t pt1 = 0;
// 		Float_t pt2 = 0;
// 		Float_t e   = 0;
// 		Float_t e1  = 0;
// 		Float_t e2  = 0;
// 		TParticle * tmptmp = new TParticle;
// 		for(Int_t i = 0; i< stack->GetNprimary(); i++){
// 			tmptmp = stack->Particle(i);
		
// 			if(tmptmp->GetStatusCode() == 1){
// 				pt = tmptmp->Pt();
// 				e =  tmptmp->Energy();			
// 				Int_t imom = tmptmp->GetFirstMother();
// 				Int_t imom1 = 0;
// 				//printf("1st imom %d\n",imom);
// 				while(imom > 5){
// 					imom1=imom;
// 					tmptmp = stack->Particle(imom);
// 					imom = tmptmp->GetFirstMother();
// 					//printf("imom %d	\n",imom);
// 				}
// 				//printf("Last imom %d %d\n",imom1, imom);
// 				if(imom1 == 6) {
// 					pt1+=pt;
// 					e1+=e;				
// 				}
// 				else if (imom1 == 7){
// 					pt2+=pt;
// 					e2+=e;					}
// 			}// status 1
				
// 		}// for
		
// 		printf("JET 1, pt %2.2f, e %2.2f; JET 2, pt %2.2f, e %2.2f \n",pt1,e1,pt2,e2);
		
		//Get the jet, different way for different generator
		//PYTHIA
    if(fMCGenerator == "PYTHIA"){
      TParticle * jet =  new TParticle;
      AliGenPythiaEventHeader* pygeh= (AliGenPythiaEventHeader*) geh;
      nTriggerJets =  pygeh->NTriggerJets();
      if(fDebug > 1)
         printf("AliMCAnalysisUtils::GetJets() - PythiaEventHeader: Njets: %d\n",nTriggerJets);
		
      Int_t iparton = -1;
      for(Int_t i = 0; i< nTriggerJets; i++){
	iparton=-1;
	pygeh->TriggerJet(i, tmpjet);
	jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0],tmpjet[1],tmpjet[2],tmpjet[3], 0,0,0,0);
	//Assign an outgoing parton as mother
	Float_t phidiff1 = TMath::Abs(jet->Phi()-parton1->Phi());		
	Float_t phidiff2 = TMath::Abs(jet->Phi()-parton2->Phi());
	if(phidiff1 > phidiff2) jet->SetFirstMother(7);
	else  jet->SetFirstMother(6);
	//jet->Print();
	if(fDebug > 1)
	  printf("AliMCAnalysisUtils::GetJets() - PYTHIA Jet %d: mother %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
		 i, jet->GetFirstMother(),jet->Pt(),jet->Energy(),jet->Phi()*TMath::RadToDeg(),jet->Eta());
	fJetsList->Add(jet);			
      }
    }//Pythia triggered jets
    //HERWIG
    else if (fMCGenerator=="HERWIG"){
      Int_t pdg = -1;		
      //Check parton 1
      TParticle * tmp = parton1;
      if(parton1->GetPdgCode()!=22){
	while(pdg != 94){
	  if(tmp->GetFirstDaughter()==-1) return fJetsList;
	  tmp = stack->Particle(tmp->GetFirstDaughter());
	  pdg = tmp->GetPdgCode();
	}//while
	
	//Add found jet to list
	TParticle *jet1 = new TParticle(*tmp);
	jet1->SetFirstMother(6);
	fJetsList->Add(jet1);
	//printf("jet 1:  first daughter %d, last daughter %d\n", tmp->GetFirstDaughter(), tmp->GetLastDaughter());
	//tmp = stack->Particle(tmp->GetFirstDaughter());
	//tmp->Print();
	//jet1->Print();
	if(fDebug > 1)			
	  printf("AliMCAnalysisUtils::GetJets() - HERWIG Jet 1: mother %d, status %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
		 jet1->GetFirstMother(),jet1->GetStatusCode(),jet1->Pt(),jet1->Energy(),jet1->Phi()*TMath::RadToDeg(),jet1->Eta());
      }//not photon
      
      //Check parton 2
      pdg = -1;
      tmp = parton2;
      Int_t i = -1;
      if(parton2->GetPdgCode()!=22){
	while(pdg != 94){
	  if(tmp->GetFirstDaughter()==-1) return fJetsList;
	  i = tmp->GetFirstDaughter();
	  tmp = stack->Particle(tmp->GetFirstDaughter());
	  pdg = tmp->GetPdgCode();
	}//while
	//Add found jet to list
	TParticle *jet2 = new TParticle(*tmp);
	jet2->SetFirstMother(7);
	fJetsList->Add(jet2);
	//jet2->Print();
	if(fDebug > 1)
	  printf("AliMCAnalysisUtils::GetJets() - HERWIG Jet 2: mother %d, status %d, pt %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
		 jet2->GetFirstMother(),jet2->GetStatusCode(),jet2->Pt(),jet2->Energy(),jet2->Phi()*TMath::RadToDeg(),jet2->Eta());
	//Int_t first =  tmp->GetFirstDaughter();
	//Int_t last  =  tmp->GetLastDaughter();
	//printf("jet 2:  first daughter %d, last daughter %d, pdg %d\n",first, last, tmp->GetPdgCode());
				//	for(Int_t d = first ; d < last+1; d++){
//						tmp = stack->Particle(d);
//						if(i == tmp->GetFirstMother())
//							printf("Daughter n %d, Mother %d, name %s, status %d, pT %2.2f,E %2.2f, phi %2.2f, eta %2.2f \n",
//							d,tmp->GetFirstMother(), tmp->GetName(), tmp->GetStatusCode(),tmp->Pt(),tmp->Energy(),tmp->Phi()*TMath::RadToDeg(),tmp->Eta());			   
//			   }
			  			   //tmp->Print();
      }//not photon
    }//Herwig generated jets
  }
  
  return fJetsList;
}


//________________________________________________________________
void AliMCAnalysisUtils::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
 
 if(! opt)
   return;
 
 printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
 
 printf("Debug level    = %d\n",fDebug);
 printf("MC Generator   = %s\n",fMCGenerator.Data());
 printf(" \n");
 
} 


