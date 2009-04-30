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

//---- ANALYSIS system ----
#include "AliMCAnalysisUtils.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenPythiaEventHeader.h"

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
Int_t AliMCAnalysisUtils::CheckOrigin(const Int_t label, AliStack * stack) const {
  //Play with the MC stack if available
  //Check origin of the candidates, good for PYTHIA
  
  if(!stack) {
    printf("AliMCAnalysisUtils::CheckOrigin() - Stack is not available, check analysis settings in configuration file, STOP!!\n");
    abort();
  }
  //  printf("label %d, ntrack %d, nprim %d\n",label, stack->GetNtrack(), stack->GetNprimary());
  //   for(Int_t i = 0; i< stack->GetNprimary(); i++){
  //      TParticle *particle =   stack->Particle(i);
  // 			//particle->Print();
  //   }
  if(label >= 0 && label <  stack->GetNtrack()){
    //Mother
    TParticle * mom = stack->Particle(label);
    Int_t mPdg = TMath::Abs(mom->GetPdgCode());
    Int_t mStatus =  mom->GetStatusCode() ;
    Int_t iParent =  mom->GetFirstMother() ;
    if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOrigin: Mother is parton %d\n",iParent);
    
    //GrandParent
    TParticle * parent = new TParticle ;
    Int_t pPdg = -1;
    Int_t pStatus =-1;
    if(iParent > 0){
      parent = stack->Particle(iParent);
      pPdg = TMath::Abs(parent->GetPdgCode());
      pStatus = parent->GetStatusCode();  
    }
    else if(fDebug > 0 ) printf("AliMCAnalysisUtils::CheckOrigin: Parent with label %d\n",iParent);
    
    //return tag
    if(mPdg == 22){
      if(mStatus == 1){
	if(fMCGenerator == "PYTHIA"){
	  if(iParent < 8 && iParent > 5) {//outgoing partons
	    if(pPdg == 22) return kMCPrompt;
	    else  return kMCFragmentation;
	  }//Outgoing partons
	  else if(pStatus == 11){//Decay
	    if(pPdg == 111) return kMCPi0Decay ;
	    else if (pPdg == 321)  return kMCEtaDecay ;
	    else  return kMCOtherDecay ;
	  }//Decay
	  else return kMCISR; //Initial state radiation
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
	      if(pPdg == 22) return kMCPrompt;
	      else  return kMCFragmentation;
	    }
	    return kMCISR;//Initial state radiation
	  }//Not decay
	  else{//Decay
	    if(pPdg == 111) return kMCPi0Decay ;
	    else if (pPdg == 321)  return kMCEtaDecay ;
	    else  return kMCOtherDecay ;
	  }//Decay
	}//HERWIG
	else return  kMCUnknown;
      }//Status 1 : Pythia generated
      else if(mStatus == 0){
	if(pPdg ==22 || pPdg ==11|| pPdg == 2112 ||  pPdg == 211 ||  
	   pPdg == 321 ||  pPdg == 2212  ||  pPdg == 130  ||  pPdg == 13 ) 
	  return kMCConversion ;
	if(pPdg == 111) return kMCPi0Decay ;
	else if (pPdg == 221)  return kMCEtaDecay ;
	else  return kMCOtherDecay ;
      }//status 0 : geant generated
    }//Mother Photon
    else if(mPdg == 111)  return kMCPi0 ;
    else if(mPdg == 221)  return kMCEta ;
    else if(mPdg ==11){
//       if(pPdg !=22 &&mStatus == 0) {
// 	printf("Origin electron, pT %f, status %d, parent %s, pstatus %d, vx %f, vy %f, vz %f\n",
// 	       mom->Pt(),mStatus, parent->GetName(),pStatus,mom->Vx(),mom->Vy(), mom->Vz());

//       }

      if(mStatus == 0) {
	if(pPdg ==22) return kMCConversion ;
	else if (pStatus == 0) return kMCConversion;
	else return kMCElectron ;
      }
      else return kMCElectron ;
    }
    else return kMCUnknown;
  }//Good label value
  else{
    if(label < 0 ) printf("AliMCAnalysisUtils::CheckOrigin: *** bad label or no stack ***:  label %d \n", label);
    if(label >=  stack->GetNtrack()) printf("AliMCAnalysisUtils::CheckOrigin: *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
    return kMCUnknown;
  }//Bad label
	
  return kMCUnknown;
  
}

//_________________________________________________________________________
TList * AliMCAnalysisUtils::GetJets(Int_t iEvent, AliStack * stack, AliGenEventHeader * geh) {
 //Return list of jets (TParticles) and index of most likely parton that originated it.
	
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


