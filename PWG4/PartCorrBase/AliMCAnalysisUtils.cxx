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
#include "TDatabasePDG.h"
#include "TVector3.h"

//---- ANALYSIS system ----
#include "AliMCAnalysisUtils.h"
#include "AliCaloTrackReader.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCParticle.h"

ClassImp(AliMCAnalysisUtils)

//________________________________________
AliMCAnalysisUtils::AliMCAnalysisUtils() : 
TObject(), 
fCurrentEvent(-1), 
fDebug(-1), 
fJetsList(new TList), 
fMCGenerator("PYTHIA")
{
  //Ctor
}

//_______________________________________
AliMCAnalysisUtils::~AliMCAnalysisUtils() 
{
  // Remove all pointers.
  
  if (fJetsList) {
    fJetsList->Clear();
    delete fJetsList ;
  }     
}

//_____________________________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckCommonAncestor(const Int_t index1, const Int_t index2, 
                                              const AliCaloTrackReader* reader, 
                                              Int_t & ancPDG, Int_t & ancStatus, 
                                              TLorentzVector & momentum, TVector3 & prodVertex) 
{
  //Check the first common ancestor of 2 clusters, given the most likely labels of the primaries generating such clusters.
  Int_t label1[100];
  Int_t label2[100];
  label1[0]= index1;
  label2[0]= index2;
  Int_t counter1 = 0;
  Int_t counter2 = 0;
  
  if(label1[0]==label2[0]) {
    //printf("AliMCAnalysisUtils::CheckCommonAncestor() - Already the same label: %d\n",label1[0]);
    counter1=1;
    counter2=1;
  }
  else{
    if(reader->ReadAODMCParticles()){
      TClonesArray * mcparticles = reader->GetAODMCParticles(0);
      
      Int_t label=label1[0];
      while(label > -1 && counter1 < 99){
        counter1++;
        AliAODMCParticle * mom = (AliAODMCParticle *) mcparticles->At(label);
        if(mom){
          label  = mom->GetMother() ;
          label1[counter1]=label;
        }
        //printf("\t counter %d, label %d\n", counter1,label);
      }
      //printf("Org label2=%d,\n",label2[0]);
      label=label2[0];
      while(label > -1 && counter2 < 99){
        counter2++;
        AliAODMCParticle * mom = (AliAODMCParticle *) mcparticles->At(label);
        if(mom){
          label  = mom->GetMother() ;
          label2[counter2]=label;
        }
        //printf("\t counter %d, label %d\n", counter2,label);
      }
    }//AOD MC
    else { //Kine stack from ESDs 
      AliStack * stack = reader->GetStack();
      Int_t label=label1[0];
      while(label > -1 && counter1 < 99){
        counter1++;
        TParticle * mom = stack->Particle(label);
        if(mom){
          label  = mom->GetFirstMother() ;
          label1[counter1]=label;
        }
        //printf("\t counter %d, label %d\n", counter1,label);
      }
      //printf("Org label2=%d,\n",label2[0]);
      label=label2[0];
      while(label > -1 && counter2 < 99){
        counter2++;
        TParticle * mom = stack->Particle(label);
        if(mom){
          label  = mom->GetFirstMother() ;
          label2[counter2]=label;
        }
        //printf("\t counter %d, label %d\n", counter2,label);
      }
    }// Kine stack from ESDs
  }//First labels not the same
  
  if((counter1==99 || counter2==99) && fDebug >=0) printf("AliMCAnalysisUtils::CheckCommonAncestor() - Genealogy too large c1: %d, c2= %d\n", counter1, counter2);
  //printf("CheckAncestor:\n");
  Int_t commonparents = 0;
  Int_t ancLabel = -1;
  //printf("counters %d %d \n",counter1, counter2);
  for (Int_t c1 = 0; c1 < counter1; c1++) {
    for (Int_t c2 = 0; c2 < counter2; c2++) {
      if(label1[c1]==label2[c2] && label1[c1]>-1) {
        ancLabel = label1[c1];
        commonparents++;
        if(reader->ReadAODMCParticles()){
          AliAODMCParticle * mom = (AliAODMCParticle *) reader->GetAODMCParticles(0)->At(label1[c1]);
          if (mom) {
            ancPDG    = mom->GetPdgCode();
            ancStatus = mom->GetStatus();
            momentum.SetPxPyPzE(mom->Px(),mom->Py(),mom->Pz(),mom->E());
            prodVertex.SetXYZ(mom->Xv(),mom->Yv(),mom->Zv());
          }
        }
        else {
          TParticle * mom = (reader->GetStack())->Particle(label1[c1]);
          if (mom) {
            ancPDG    = mom->GetPdgCode();
            ancStatus = mom->GetStatusCode();
            mom->Momentum(momentum);
            prodVertex.SetXYZ(mom->Vx(),mom->Vy(),mom->Vz());
          }
        }
        //First ancestor found, end the loops
        counter1=0;
        counter2=0;
      }//Ancestor found
    }//second cluster loop
  }//first cluster loop
  
  return ancLabel;
}

//_____________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOrigin(const Int_t * label, 
                                      const Int_t nlabels, 
                                      const AliCaloTrackReader* reader, 
                                      const Int_t input = 0) 
{
  //Play with the montecarlo particles if available
  Int_t tag = 0;
  
  if(nlabels<=0) {
    printf("AliMCAnalysisUtils::CheckOrigin(nlabel<=0) - No MC labels available, please check!!!\n");
    return kMCBadLabel;
  }
  
  //Select where the information is, ESD-galice stack or AOD mcparticles branch
  if(reader->ReadStack()){
    tag = CheckOriginInStack(label, nlabels, reader->GetStack());
  }
  else if(reader->ReadAODMCParticles()){
    tag = CheckOriginInAOD(label, nlabels, reader->GetAODMCParticles(input));
  }
  
  return tag ;
}

//_____________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOrigin(const Int_t label, 
                                      const AliCaloTrackReader* reader, 
                                      const Int_t input = 0) 
{
  //Play with the montecarlo particles if available
  Int_t tag = 0;
  
  if(label<0) {
    printf("AliMCAnalysisUtils::CheckOrigin(label<0) - No MC labels available, please check!!!\n");
    return kMCBadLabel;
  }
  
  Int_t labels[]={label};
  
  //Select where the information is, ESD-galice stack or AOD mcparticles branch
  if(reader->ReadStack()){
    tag = CheckOriginInStack(labels, 1,reader->GetStack());
  }
  else if(reader->ReadAODMCParticles()){
    tag = CheckOriginInAOD(labels, 1,reader->GetAODMCParticles(input));
  }
  
  return tag ;
}	

//_________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInStack(const Int_t *labels, 
                                             const Int_t nlabels, 
                                             AliStack* stack) 
{
  // Play with the MC stack if available. Tag particles depending on their origin.
  // Do same things as in CheckOriginInAOD but different input.
  
  //generally speaking, label is the MC label of a reconstructed
  //entity (track, cluster, etc) for which we want to know something 
  //about its heritage, but one can also use it directly with stack 
  //particles not connected to reconstructed entities
  
  if(!stack) {
    if (fDebug >=0) 
      printf("AliMCAnalysisUtils::CheckOriginInStack() - Stack is not available, check analysis settings in configuration file, STOP!!\n");
    return -1;
  }
  
  Int_t tag = 0;
  Int_t label=labels[0];//Most significant particle contributing to the cluster
  
  if(label >= 0 && label < stack->GetNtrack()){
    //MC particle of interest is the "mom" of the entity
    TParticle * mom = stack->Particle(label);
    Int_t iMom     = label;
    Int_t mPdgSign = mom->GetPdgCode();
    Int_t mPdg     = TMath::Abs(mPdgSign);
    Int_t mStatus  = mom->GetStatusCode() ;
    Int_t iParent  = mom->GetFirstMother() ;
    if(fDebug > 0 && label < 8 && fMCGenerator!="") printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
    
    //GrandParent of the entity
    TParticle * parent = NULL;
    Int_t pPdg = -1;
    Int_t pStatus =-1;
    if(iParent >= 0){
      parent = stack->Particle(iParent);
      if(parent){
        pPdg = TMath::Abs(parent->GetPdgCode());
        pStatus = parent->GetStatusCode();  
      }
    }
    else if(fDebug > 0 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Parent with label %d\n",iParent);
    
    if(fDebug > 2 ) {
      printf("AliMCAnalysisUtils::CheckOriginInStack() - Cluster most contributing mother and its parent: \n");
      printf("\t Mother label %d, pdg %d, status %d\n",iMom, mPdg, mStatus);
      printf("\t Parent label %d, pdg %d, status %d\n",iParent, pPdg, pStatus);
    }
	  
    //Check if "mother" of entity is converted, if not, get the first non converted mother
    if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && mStatus == 0){
      SetTagBit(tag,kMCConversion);
      //Check if the mother is photon or electron with status not stable
      while ((pPdg == 22 || pPdg == 11) && mStatus != 1) {
        //Mother
        iMom     = mom->GetFirstMother();
        mom      = stack->Particle(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        mStatus  = mom->GetStatusCode() ;
        iParent  = mom->GetFirstMother() ;
        if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
        
        //GrandParent
        if(iParent >= 0){
          parent = stack->Particle(iParent);
          if(parent){
            pPdg = TMath::Abs(parent->GetPdgCode());
            pStatus = parent->GetStatusCode();  
          }
        }
        else {// in case of gun/box simulations
          pPdg    = 0;
          pStatus = 0;
          break;
        }
      }//while	  
      if(fDebug > 2 ) {
        printf("AliMCAnalysisUtils::CheckOriginInStack() - Converted photon/electron: \n");
        printf("\t Mother label %d, pdg %d, status %d\n",iMom, mPdg, mStatus);
        printf("\t Parent label %d, pdg %d, status %d\n",iParent, pPdg, pStatus);
      }
      
    }//mother and parent are electron or photon and have status 0
    else if((mPdg == 22 || mPdg == 11) && mStatus == 0){	
      //Still a conversion but only one electron/photon generated. Just from hadrons but not decays.
      if(pPdg == 2112 ||  pPdg == 211  ||  pPdg == 321 ||
         pPdg == 2212 ||  pPdg == 130  ||  pPdg == 13 ) {
        SetTagBit(tag,kMCConversion);
        iMom     = mom->GetFirstMother();
        mom      = stack->Particle(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        
        if(fDebug > 2 ) {
          printf("AliMCAnalysisUtils::CheckOriginInStack() - Converted hadron: \n");
          printf("\t Mother label %d, pdg %d, status %d\n",iMom, mPdg, mStatus);
        }
      }//hadron converted
      
      //Comment for the next lines, we do not check the parent of the hadron for the moment.
      //iParent =  mom->GetFirstMother() ;
      //if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Mother is parton %d\n",iParent);
      
      //GrandParent
      //if(iParent >= 0){
      //	parent = stack->Particle(iParent);
      //	pPdg = TMath::Abs(parent->GetPdgCode());
      //}
    }  	  
    // conversion into electrons/photons checked  	  
    
    //first check for typical charged particles
    if     (mPdg     ==    13) SetTagBit(tag,kMCMuon);
    else if(mPdg     ==   211) SetTagBit(tag,kMCPion);
    else if(mPdg     ==   321) SetTagBit(tag,kMCKaon);
    else if(mPdgSign ==  2212) SetTagBit(tag,kMCProton);
    else if(mPdgSign == -2212) SetTagBit(tag,kMCAntiProton);
    else if(mPdgSign ==  2112) SetTagBit(tag,kMCNeutron);
    else if(mPdgSign == -2112) SetTagBit(tag,kMCAntiNeutron);
    
    //check for pi0 and eta (shouldn't happen unless their decays were turned off)
    else if(mPdg == 111)  {
      SetTagBit(tag,kMCPi0Decay);
      if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - First mother is directly pi0, not decayed by generator \n");
      CheckOverlapped2GammaDecay(labels,nlabels, iMom, stack, tag); //set to kMCPi0 if 2 gammas in same cluster
    }
    else if(mPdg == 221) {
      SetTagBit(tag,kMCEtaDecay);
      if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - First mother is directly eta, not decayed by generator \n");
      CheckOverlapped2GammaDecay(labels,nlabels, iMom, stack, tag); //set to kMCEta if 2 gammas in same cluster
    }
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
            if(pPdg == 111) {
              SetTagBit(tag,kMCPi0Decay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - PYTHIA pi0 decay photon,  parent pi0 with status 11 \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCPi0 if 2 gammas in same cluster
            }
            else if (pPdg == 221) {
              SetTagBit(tag, kMCEtaDecay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - PYTHIA eta decay photon,  parent pi0 with status 11 \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag);//set to kMCEta if 2 gammas in same cluster
            }
            else SetTagBit(tag,kMCOtherDecay);
          }//Decay
          else {
            if(fDebug > 1 && parent) printf("AliMCAnalysisUtils::CheckOrigingInStack() - what is it in PYTHIA? Wrong generator setting? Mother mPdg %d, status %d \n    Parent  iParent %d, pPdg %d %s, status %d\n",
                                            mPdg, mStatus,iParent, pPdg, parent->GetName(),pStatus);
            if(pPdg == 111) {
              SetTagBit(tag,kMCPi0Decay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - PYTHIA pi0 decay photon,  parent pi0 with status 11 \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCPi0 if 2 gammas in same cluster
            }
            else if (pPdg == 221) {
              SetTagBit(tag, kMCEtaDecay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - PYTHIA eta decay photon,  parent pi0 with status 11 \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag);//set to kMCEta if 2 gammas in same cluster
            }
            else SetTagBit(tag,kMCOtherDecay);
          }
        }//PYTHIA
        
        else if(fMCGenerator == "HERWIG"){	  
          if(pStatus < 197){//Not decay
            while(1){
              if(parent){
                if(parent->GetFirstMother()<=5) break;
                iParent = parent->GetFirstMother();
                parent=stack->Particle(iParent);
                pStatus= parent->GetStatusCode();
                pPdg = TMath::Abs(parent->GetPdgCode());
              } else break;
            }//Look for the parton
            
            if(iParent < 8 && iParent > 5) {
              if(pPdg == 22) SetTagBit(tag,kMCPrompt);
              else SetTagBit(tag,kMCFragmentation);
            }
            else SetTagBit(tag,kMCISR);//Initial state radiation
          }//Not decay
          else{//Decay
            if(pPdg == 111) {
              SetTagBit(tag,kMCPi0Decay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - HERWIG pi0 decay photon \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCPi0 if 2 gammas in same cluster
            }
            else if (pPdg == 221) {
              SetTagBit(tag,kMCEtaDecay);
              if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - HERWIG eta decay photon \n");
              CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCEta if 2 gammas in same cluster
            }
            else SetTagBit(tag,kMCOtherDecay);
          }//Decay
        }//HERWIG
        
        else SetTagBit(tag,kMCUnknown);
        
      }//Status 1 : created by event generator
      
      else if(mStatus == 0){ // geant
        if(pPdg == 111) {
          SetTagBit(tag,kMCPi0Decay);
          if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Transport MC pi0 decay photon \n");
          CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCPi0 if 2 gammas in same cluster
        }
        else if (pPdg == 221) {
          SetTagBit(tag,kMCEtaDecay);
          if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInStack() - Transport MC eta decay photon \n");
          CheckOverlapped2GammaDecay(labels,nlabels, iParent, stack, tag); //set to kMCEta if 2 gammas in same cluster
        }
        else  SetTagBit(tag,kMCOtherDecay);	
      }//status 0 : geant generated
      
    }//Mother Photon
    
    //Electron check.  Where did that electron come from?
    else if(mPdg == 11){ //electron
      if(pPdg == 11 && parent){
        Int_t iGrandma = parent->GetFirstMother();
        if(iGrandma >= 0) {
          TParticle* gma = (TParticle*)stack->Particle(iGrandma); //get mother
          Int_t gPdg = TMath::Abs(gma->GetPdgCode());
          
          if (gPdg == 23) { SetTagBit(tag,kMCZDecay); } //parent is Z-boson
          else if (gPdg == 24) { SetTagBit(tag,kMCWDecay); } //parent is W-boson
        }
      }
      SetTagBit(tag,kMCElectron);	
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInStack() - Checking ancestors of electrons\n");
      if (pPdg == 111) { SetTagBit(tag,kMCPi0Decay); } //Pi0 Dalitz decay
      else if (pPdg == 221) { SetTagBit(tag,kMCEtaDecay); } //Eta Dalitz decay
      else if((499 < pPdg && pPdg < 600)||(4999 < pPdg && pPdg < 6000)) { SetTagBit(tag,kMCEFromB); } //b-->e decay
      else if((399 < pPdg && pPdg < 500)||(3999 < pPdg && pPdg < 5000)) { //check charm decay
        if(parent){
          Int_t iGrandma = parent->GetFirstMother();
          if(iGrandma >= 0) {
            TParticle* gma = (TParticle*)stack->Particle(iGrandma); //get mother of charm
            Int_t gPdg = TMath::Abs(gma->GetPdgCode());
            if((499 < gPdg && gPdg < 600)||(4999 < gPdg && gPdg < 6000)) SetTagBit(tag,kMCEFromCFromB); //b-->c-->e
            else SetTagBit(tag,kMCEFromC); //c-->e 
          } else SetTagBit(tag,kMCEFromC); //c-->e 
        }//parent
      } else {
        //if it is not from any of the above, where is it from?
        if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
        else SetTagBit(tag,kMCOtherDecay);
        if(fDebug > 0 && parent) printf("AliMCAnalysisUtils::CheckOriginInStack() - Status %d Electron from other origin: %s (pPdg = %d) %s (mpdg = %d)\n",mStatus,parent->GetName(),pPdg,mom->GetName(),mPdg);
      }
    }//electron check
    //Cluster was made by something else
    else {
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInStack() - \tSetting kMCUnknown for cluster from %s (pdg = %d, Parent pdg = %d)\n",mom->GetName(),mPdg,pPdg);
      SetTagBit(tag,kMCUnknown);
    }
  }//Good label value
  else{// Bad label 
	  
    if(label < 0 && (fDebug >= 0)) 
      printf("AliMCAnalysisUtils::CheckOriginInStack() *** bad label or no stack ***:  label %d \n", label);
    if(label >=  stack->GetNtrack() &&  (fDebug >= 0)) 
      printf("AliMCAnalysisUtils::CheckOriginInStack() *** large label ***:  label %d, n tracks %d \n", label, stack->GetNtrack());
    SetTagBit(tag,kMCUnknown);
  }//Bad label
  
  return tag;
  
}


//_________________________________________________________________________
Int_t AliMCAnalysisUtils::CheckOriginInAOD(const Int_t *labels, 
                                           const Int_t nlabels, 
                                           const TClonesArray *mcparticles) 
{
  // Play with the MCParticles in AOD if available. Tag particles depending on their origin.
  // Do same things as in CheckOriginInStack but different input.
  if(!mcparticles) {
    if(fDebug >= 0)
      printf("AliMCAnalysisUtils::CheckOriginInAOD() - AODMCParticles is not available, check analysis settings in configuration file!!\n");
    return -1;
  }
	
  Int_t tag = 0;
  Int_t label=labels[0];//Most significant particle contributing to the cluster
  
  Int_t nprimaries = mcparticles->GetEntriesFast();
  if(label >= 0 && label < nprimaries){
    //Mother
    AliAODMCParticle * mom = (AliAODMCParticle *) mcparticles->At(label);
    Int_t iMom     = label;
    Int_t mPdgSign = mom->GetPdgCode();
    Int_t mPdg     = TMath::Abs(mPdgSign);
    Int_t iParent  = mom->GetMother() ;
    if(fDebug > 0 && label < 8 && fMCGenerator!="") printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
    
    //GrandParent
    AliAODMCParticle * parent = NULL ;
    Int_t pPdg = -1;
    if(iParent >= 0){
      parent = (AliAODMCParticle *) mcparticles->At(iParent);
      pPdg = TMath::Abs(parent->GetPdgCode());
    }
    else if(fDebug > 0 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Parent with label %d\n",iParent);
    
    if(fDebug > 2 ) {
      printf("AliMCAnalysisUtils::CheckOriginInAOD() - Cluster most contributing mother and its parent: \n");
      printf("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary());
      if(parent)
        printf("\t Parent label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iParent, pPdg, parent->IsPrimary(), parent->IsPhysicalPrimary());
    }
	  
    //Check if mother is converted, if not, get the first non converted mother
    if((mPdg == 22 || mPdg == 11) && (pPdg == 22 || pPdg == 11) && !mom->IsPrimary()){
      SetTagBit(tag,kMCConversion);
      //Check if the mother is photon or electron with status not stable
      while ((pPdg == 22 || pPdg == 11) && !mom->IsPhysicalPrimary()) {
        //Mother
        iMom     = mom->GetMother();
        mom      = (AliAODMCParticle *) mcparticles->At(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        iParent  = mom->GetMother() ;
        if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
        
        //GrandParent
        if(iParent >= 0 && parent){
          parent = (AliAODMCParticle *) mcparticles->At(iParent);
          pPdg = TMath::Abs(parent->GetPdgCode());
        }
        // printf("\t While Mother label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary());
        // printf("\t While Parent label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iParent, pPdg, parent->IsPrimary(), parent->IsPhysicalPrimary()); 
        
      }//while	
      
      if(fDebug > 2 ) {
        printf("AliMCAnalysisUtils::CheckOriginInAOD() - Converted photon/electron : \n");
        printf("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary());
        if(parent)
          printf("\t Parent label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iParent, pPdg, parent->IsPrimary(), parent->IsPhysicalPrimary());
      }
      
    }//mother and parent are electron or photon and have status 0 and parent is photon or electron
    else if((mPdg == 22 || mPdg == 11) && !mom->IsPrimary()){	
      //Still a conversion but only one electron/photon generated. Just from hadrons
      if(pPdg == 2112 ||  pPdg == 211 ||  pPdg == 321 ||  
         pPdg == 2212 ||  pPdg == 130 ||  pPdg == 13 ) {
        SetTagBit(tag,kMCConversion);
        iMom     = mom->GetMother();
        mom      = (AliAODMCParticle *) mcparticles->At(iMom);
        mPdgSign = mom->GetPdgCode();
        mPdg     = TMath::Abs(mPdgSign);
        
        if(fDebug > 2 ) {
          printf("AliMCAnalysisUtils::CheckOriginInAOD() - Converted hadron : \n");
          printf("\t Mother label %d, pdg %d, Primary? %d, Physical Primary? %d\n",iMom, mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary());
        }
      }//hadron converted
      
      //Comment for next lines, we do not check the parent of the hadron for the moment.
      //iParent =  mom->GetMother() ;
      //if(fDebug > 0 && label < 8 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Mother is parton %d\n",iParent);
      
      //GrandParent
      //if(iParent >= 0){
      //	parent = (AliAODMCParticle *) mcparticles->At(iParent);
      //	pPdg = TMath::Abs(parent->GetPdgCode());
      //}
    }  
    
    //printf("Final mother mPDG %d\n",mPdg);
    
    // conversion into electrons/photons checked  
    
    //first check for typical charged particles
    if     (mPdg     ==    13) SetTagBit(tag,kMCMuon);
    else if(mPdg     ==   211) SetTagBit(tag,kMCPion);
    else if(mPdg     ==   321) SetTagBit(tag,kMCKaon);
    else if(mPdgSign ==  2212) SetTagBit(tag,kMCProton);
    else if(mPdgSign ==  2112) SetTagBit(tag,kMCNeutron);
    else if(mPdgSign == -2212) SetTagBit(tag,kMCAntiProton);
    else if(mPdgSign == -2112) SetTagBit(tag,kMCAntiNeutron);
    
    //check for pi0 and eta (shouldn't happen unless their decays were turned off)
    else if(mPdg == 111)  {
      SetTagBit(tag,kMCPi0Decay);
      if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - First mother is directly pi0, not decayed by generator \n");
      CheckOverlapped2GammaDecay(labels,nlabels, iMom, mcparticles, tag); //set to kMCPi0 if 2 gammas in same cluster
    }
    else if(mPdg == 221)  {
      SetTagBit(tag,kMCEtaDecay);   
      if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - First mother is directly eta, not decayed by generator \n");
      CheckOverlapped2GammaDecay(labels,nlabels, iMom, mcparticles, tag); //set to kMCEta if 2 gammas in same cluster
    }
    //Photons  
    else if(mPdg == 22){
      SetTagBit(tag,kMCPhoton);
      if(mom->IsPhysicalPrimary() && (fMCGenerator=="PYTHIA" || fMCGenerator=="HERWIG")) //undecayed particle
      {
        if(iParent < 8 && iParent > 5 ) {//outgoing partons
          if(pPdg == 22) SetTagBit(tag,kMCPrompt);
          else SetTagBit(tag,kMCFragmentation);
        }//Outgoing partons
        else if(iParent <= 5 && (fMCGenerator=="PYTHIA" || fMCGenerator=="HERWIG")) {
          SetTagBit(tag, kMCISR); //Initial state radiation
        }
        else if(parent && parent->IsPrimary() && !parent->IsPhysicalPrimary()){//Decay
          if(pPdg == 111){
            SetTagBit(tag,kMCPi0Decay);
            if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Generator pi0 decay photon \n");
            CheckOverlapped2GammaDecay(labels,nlabels, iParent, mcparticles, tag); //set to kMCPi0 if 2 gammas in same cluster
          }
          else if (pPdg == 221) {
            SetTagBit(tag, kMCEtaDecay);
            if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Generator eta decay photon \n");
            CheckOverlapped2GammaDecay(labels,nlabels, iParent, mcparticles, tag); //set to kMCEta if 2 gammas in same cluster
          }
          else SetTagBit(tag,kMCOtherDecay);
        }//Decay
        else {
          if(parent)printf("AliMCAnalysisUtils::CheckOriginInAOD() - what is it? Mother mPdg %d, is primary? %d, is physical %d \n    Parent  iParent %d, pPdg %d, is primary? %d, is physical? %d\n",
                           mPdg, mom->IsPrimary(), mom->IsPhysicalPrimary(),iParent, pPdg,parent->IsPrimary(), parent->IsPhysicalPrimary());
          SetTagBit(tag,kMCOtherDecay);//Check
        }
      }//Physical primary
      else if(!mom->IsPrimary()){	//Decays  
        if(pPdg == 111){ 
          SetTagBit(tag,kMCPi0Decay); 
          if(fDebug > 2 ) 
            printf("AliMCAnalysisUtils::CheckOriginInAOD() - Transport MC pi0 decay photon \n");
          CheckOverlapped2GammaDecay(labels,nlabels, iParent, mcparticles, tag); //set to kMCPi0 if 2 gammas in same cluster          
        }
        else if (pPdg == 221) {
          SetTagBit(tag,kMCEtaDecay);
          if(fDebug > 2 ) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Transport MC eta decay photon \n");
          CheckOverlapped2GammaDecay(labels,nlabels, iParent, mcparticles, tag); //set to kMCEta if 2 gammas in same cluster
        }
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
      if(pPdg == 11 && parent){
        Int_t iGrandma = parent->GetMother();
        if(iGrandma >= 0) {
          AliAODMCParticle* gma = (AliAODMCParticle*)mcparticles->At(iGrandma);
          Int_t gPdg = TMath::Abs(gma->GetPdgCode());
          
          if (gPdg == 23) { SetTagBit(tag,kMCZDecay); } //parent is Z-boson
          else if (gPdg == 24) { SetTagBit(tag,kMCWDecay); } //parent is W-boson
        }
      }
      SetTagBit(tag,kMCElectron);	
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Checking ancestors of electrons");
      if (pPdg == 111) { SetTagBit(tag,kMCPi0Decay); } //Pi0 Dalitz decay
      else if (pPdg == 221) { SetTagBit(tag,kMCEtaDecay); } //Eta Dalitz decay
      else if((499 < pPdg && pPdg < 600)||(4999 < pPdg && pPdg < 6000)) { SetTagBit(tag,kMCEFromB);} //b-hadron decay
      else if((399 < pPdg && pPdg < 500)||(3999 < pPdg && pPdg < 5000)) { //c-hadron decay check
        if(parent){
          Int_t iGrandma = parent->GetMother();
          if(iGrandma >= 0) {
            AliAODMCParticle* gma = (AliAODMCParticle*)mcparticles->At(iGrandma); //charm's mother
            Int_t gPdg = TMath::Abs(gma->GetPdgCode());
            if((499 < gPdg && gPdg < 600)||(4999 < gPdg && gPdg < 6000)) SetTagBit(tag,kMCEFromCFromB); //b-->c-->e decay
            else SetTagBit(tag,kMCEFromC); //c-hadron decay
          } else SetTagBit(tag,kMCEFromC); //c-hadron decay
        }//parent
      } else { //prompt or other decay
        TParticlePDG* foo = TDatabasePDG::Instance()->GetParticle(pPdg);
        TParticlePDG* foo1 = TDatabasePDG::Instance()->GetParticle(mPdg);
        if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInAOD() - Electron from other origin: %s (pPdg = %d) %s (mPdg = %d)\n",foo->GetName(), pPdg,foo1->GetName(),mPdg);
        if(pPdg > 10000) SetTagBit(tag,kMCUnknown);
        else SetTagBit(tag,kMCOtherDecay);
      }      
    }//electron check
    //cluster was made by something else
    else {
      if(fDebug > 0) printf("AliMCAnalysisUtils::CheckOriginInAOD() - \tSetting kMCUnknown for cluster with pdg = %d, Parent pdg = %d\n",mPdg,pPdg);
      SetTagBit(tag,kMCUnknown);
    }
  }//Good label value
  else{//Bad label
	  
    if(label < 0 && (fDebug >= 0) ) 
      printf("AliMCAnalysisUtils::CheckOriginInAOD() *** bad label or no mcparticles ***:  label %d \n", label);
    if(label >=  mcparticles->GetEntriesFast() &&  (fDebug >= 0) ) 
      printf("AliMCAnalysisUtils::CheckOriginInAOD() *** large label ***:  label %d, n tracks %d \n", label, mcparticles->GetEntriesFast());
    SetTagBit(tag,kMCUnknown);
    
  }//Bad label
  
  return tag;
  
}

//_________________________________________________________________________
void AliMCAnalysisUtils::CheckOverlapped2GammaDecay(const Int_t *labels, 
                                                    const Int_t nlabels, 
                                                    const Int_t mesonIndex, 
                                                    AliStack *stack, 
                                                    Int_t &tag)
{
  //Check if cluster is formed from the contribution of 2 decay photons from pi0 or eta. Input in stack
  
  if(labels[0] < 0 || labels[0] > stack->GetNtrack() || nlabels <= 1) {
    if(fDebug > 2) printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Exit : label[0] %d, n primaries %d, nlabels %d \n",
                          labels[0],stack->GetNtrack(), nlabels);
    return;
  }
  
  TParticle * meson = stack->Particle(mesonIndex);
  Int_t mesonPdg    = meson->GetPdgCode();
  if(mesonPdg!=111 && mesonPdg!=221){ 
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Wrong pi0/eta PDG : %d \n",mesonPdg);
    return;
  }
  
  if(fDebug > 2) printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - %s, label %d\n",meson->GetName(), mesonIndex);
  
  //Check if meson decayed into 2 daughters or if both were kept.
  if(meson->GetNDaughters() != 2){
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Not overalapped. Number of daughters is %d, not 2 \n",meson->GetNDaughters());
    return;
  }
  
  //Get the daughters
  Int_t iPhoton0 = meson->GetDaughter(0);
  Int_t iPhoton1 = meson->GetDaughter(1);
  TParticle *photon0 = stack->Particle(iPhoton0);
  TParticle *photon1 = stack->Particle(iPhoton1);
  
  //Check if both daughters are photons
  if(photon0->GetPdgCode() != 22 || photon1->GetPdgCode()!=22){
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Not overalapped. PDG:  daughter 1 = %d, of daughter 2 = %d \n",photon0->GetPdgCode(),photon1->GetPdgCode());
    return;
  }
  
  if(fDebug > 2) 
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Daughter labels : photon0 = %d, photon1 = %d \n",iPhoton0,iPhoton1);
  
  //Check if both photons contribute to the cluster
  Bool_t okPhoton0 = kFALSE;
  Bool_t okPhoton1 = kFALSE;
  
  if(fDebug > 3) printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Labels loop:\n");
  
  for(Int_t i = 0; i < nlabels; i++){
    if(fDebug > 3) printf("\t  at begin:label %d/%d: %d, ok? photon1 %d, photon2 %d\n", i+1, nlabels, labels[i], okPhoton0, okPhoton1);
    
    //If we already found both, break the loop
    if(okPhoton0 && okPhoton1) break;
    
    Int_t index = 	labels[i];
    if      (iPhoton0 == index) {
      okPhoton0 = kTRUE;
      continue;
    }
    else if (iPhoton1 == index) {
      okPhoton1 = kTRUE;
      continue;
    }
    
    //Trace back the mother in case it was a conversion
    
    if(index >= stack->GetNtrack()){
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(ESD) Particle index %d larger than size of list %d !!\n",index,stack->GetNtrack());
      continue;
    }
    
    TParticle * daught = stack->Particle(index);
    Int_t tmpindex = daught->GetFirstMother();		
    if(fDebug > 3) printf("\t Conversion? : mother %d\n",tmpindex);
    while(tmpindex>=0){
      //MC particle of interest is the mother
      if(fDebug > 3) printf("\t \t parent index %d\n",tmpindex);
      daught   = stack->Particle(tmpindex);
      if      (iPhoton0 == tmpindex) {
        okPhoton0 = kTRUE;
        break;
      }
      else if (iPhoton1 == tmpindex) {
        okPhoton1 = kTRUE;
        break;
      }
      tmpindex = daught->GetFirstMother();
    }//While to check if pi0/eta daughter was one of these contributors to the cluster
    
    if(i == 0 && (!okPhoton0 && !okPhoton1) && fDebug>=0) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Something happens, first label should be from a photon decay!\n");
    
  }//loop on list of labels
  
  //If both photons contribute tag as the corresponding meson.
  if(okPhoton0 && okPhoton1){
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - %s OVERLAPPED DECAY \n", meson->GetName());
    
    if(mesonPdg == 111) SetTagBit(tag,kMCPi0);
    else  SetTagBit(tag,kMCEta);
  }
  
}	

//__________________________________________________________________________________
void AliMCAnalysisUtils::CheckOverlapped2GammaDecay(const Int_t *labels,
                                                    const Int_t nlabels, 
                                                    const Int_t mesonIndex, 
                                                    const TClonesArray *mcparticles, 
                                                    Int_t & tag )
{
  //Check if cluster is formed from the contribution of 2 decay photons from pi0 or eta. Input in AODMCParticles
  
  if(labels[0] < 0 || labels[0] > mcparticles->GetEntriesFast() || nlabels <= 1) {
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Exit : label[0] %d, n primaries %d, nlabels %d \n",
             labels[0],mcparticles->GetEntriesFast(), nlabels);
    return;
  }
  
  AliAODMCParticle * meson = (AliAODMCParticle *) mcparticles->At(mesonIndex);
  Int_t mesonPdg = meson->GetPdgCode();
  if(mesonPdg != 111 && mesonPdg != 221) {
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Wrong pi0/eta PDG : %d \n",mesonPdg);
    return;
  }
  
  if(fDebug > 2) printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - pdg %d, label %d, ndaughters %d\n", mesonPdg, mesonIndex, meson->GetNDaughters());
  
  
  //Get the daughters
  if(meson->GetNDaughters() != 2){
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(stack) - Not overalapped. Number of daughters is %d, not 2 \n",meson->GetNDaughters());
    return;
  }
  Int_t iPhoton0 = meson->GetDaughter(0);
  Int_t iPhoton1 = meson->GetDaughter(1);
  //if((iPhoton0 == -1) || (iPhoton1 == -1)){
  //	if(fDebug > 2) 
  //		printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Exit : Not overlapped. At least a daughter do not exists : d1 %d, d2 %d \n", iPhoton0, iPhoton1);
  //	return;
  //}	
  AliAODMCParticle *photon0 = (AliAODMCParticle *) mcparticles->At(iPhoton0);
  AliAODMCParticle *photon1 = (AliAODMCParticle *) mcparticles->At(iPhoton1);
  
  //Check if both daughters are photons
  if(photon0->GetPdgCode() != 22 && photon1->GetPdgCode()!=22){
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Not overlapped. PDG:  daughter 1 = %d, of daughter 2 = %d \n",photon0->GetPdgCode(),photon1->GetPdgCode());
    return;
  }
  
  if(fDebug > 2) 
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Daughter labels : photon0 = %d, photon1 = %d \n",iPhoton0,iPhoton1);
  
  //Check if both photons contribute to the cluster
  Bool_t okPhoton0 = kFALSE;
  Bool_t okPhoton1 = kFALSE;
  
  if(fDebug > 3) 
    printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Labels loop:\n");
  
  for(Int_t i = 0; i < nlabels; i++){
    if(fDebug > 3) 
      printf("\t label %d/%d: %d, ok? %d, %d\n", i, nlabels, labels[i], okPhoton0, okPhoton1);
    
    //If we already found both, break the loop
    if(okPhoton0 && okPhoton1) break;
    
    Int_t index = 	labels[i];
    if      (iPhoton0 == index) {
      okPhoton0 = kTRUE;
      continue;
    }
    else if (iPhoton1 == index) {
      okPhoton1 = kTRUE;
      continue;
    }
    
    //Trace back the mother in case it was a conversion
    
    if(index >= mcparticles->GetEntriesFast()){
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) Particle index %d larger than size of list %d !!\n",index,mcparticles->GetEntriesFast());
      continue;
    }
    
    AliAODMCParticle * daught = (AliAODMCParticle*) mcparticles->At(index);
    Int_t tmpindex = daught->GetMother();
    if(fDebug > 3) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Conversion? : mother %d\n",tmpindex);
    
    while(tmpindex>=0){
      
      //MC particle of interest is the mother
      if(fDebug > 3) 
        printf("\t parent index %d\n",tmpindex);
      daught   = (AliAODMCParticle*) mcparticles->At(tmpindex);
      //printf("tmpindex %d\n",tmpindex);
      if      (iPhoton0 == tmpindex) {
        okPhoton0 = kTRUE;
        break;
      }
      else if (iPhoton1 == tmpindex) {
        okPhoton1 = kTRUE;
        break;
      }		
      tmpindex = daught->GetMother();
    }//While to check if pi0/eta daughter was one of these contributors to the cluster
    
    if(i == 0 && (!okPhoton0 && !okPhoton1) && fDebug>=-1 ) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - Something happens, first label should be from a photon decay!\n");
    
  }//loop on list of labels
  
  //If both photons contribute tag as the corresponding meson.
  if(okPhoton0 && okPhoton1){
    if(fDebug > 2) 
      printf("AliMCAnalysisUtils::CheckOverlapped2GammaDecay(AOD) - %s OVERLAPPED DECAY \n",(TDatabasePDG::Instance()->GetParticle(mesonPdg))->GetName());
    
    if(mesonPdg == 111) SetTagBit(tag,kMCPi0);
    else                SetTagBit(tag,kMCEta);
  }	
  
}

//_________________________________________________________________________
TList * AliMCAnalysisUtils::GetJets(const AliCaloTrackReader * reader)
{
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
      TParticle * jet =  0x0;
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


//________________________________________________________
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


