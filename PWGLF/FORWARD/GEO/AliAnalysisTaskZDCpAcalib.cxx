/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//							   //
//	Class to analyze ZDC data			   //
//							   //
/////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskZDCpAcalib.h"

ClassImp(AliAnalysisTaskZDCpAcalib)


//________________________________________________________________________
AliAnalysisTaskZDCpAcalib::AliAnalysisTaskZDCpAcalib():
  AliAnalysisTaskSE(),
    fDebug(0),
    fAnalysisInput("ESD"),
    fIsMCInput(kFALSE),
    fUseSpecialOutput(kFALSE),
    fOutput(0x0),
    fCentralityTree(0x0),
    fTrigClass(""),
    fIsEventSelected(kFALSE),
    fIsV0ATriggered(0),
    fIsV0CTriggered(0)
{   
   // Default constructor

  for(Int_t itow=0; itow<5; itow++){
     fZNCtower[itow]=0.;  
     fZNAtower[itow]=0.;  
     fZNCtowerLG[itow]=0.;
     fZNAtowerLG[itow]=0.;
     fZPCtower[itow]=0.;  
     fZPAtower[itow]=0.;  
     fZPCtowerLG[itow]=0.;
     fZPAtowerLG[itow]=0.;

  }
  
  for(Int_t ihit=0; ihit<4; ihit++){
    fZNCtdc[ihit]=9999;
    fZNAtdc[ihit]=9999;
    fZPCtdc[ihit]=9999;
    fZPAtdc[ihit]=9999;
  }
}   

//________________________________________________________________________
AliAnalysisTaskZDCpAcalib::AliAnalysisTaskZDCpAcalib(const char *name):
  AliAnalysisTaskSE(name),
    fDebug(0),
    fAnalysisInput("ESD"),
    fIsMCInput(kFALSE),
    fUseSpecialOutput(kFALSE),
    fOutput(0x0),
    fCentralityTree(0x0),
    fTrigClass(""),
    fIsEventSelected(kFALSE),
    fIsV0ATriggered(0),
    fIsV0CTriggered(0)
{
  // Default constructor
 
  for(Int_t itow=0; itow<5; itow++){
     fZNCtower[itow]=0.;  
     fZNAtower[itow]=0.;  
     fZNCtowerLG[itow]=0.;
     fZNAtowerLG[itow]=0.;
     fZPCtower[itow]=0.;  
     fZPAtower[itow]=0.;  
     fZPCtowerLG[itow]=0.;
     fZPAtowerLG[itow]=0.;

  }
  
  for(Int_t ihit=0; ihit<4; ihit++){
     fZNCtdc[ihit]=9999;
     fZNAtdc[ihit]=9999;
     fZPCtdc[ihit]=9999;
     fZPAtdc[ihit]=9999;
  }
  
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); 
  
}
 
//________________________________________________________________________
AliAnalysisTaskZDCpAcalib::~AliAnalysisTaskZDCpAcalib()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  } 
  if(fCentralityTree && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fCentralityTree; fCentralityTree=0;
  } 
  
}  

//________________________________________________________________________
void AliAnalysisTaskZDCpAcalib::UserCreateOutputObjects()
{
  // Create the output containers
  if(fDebug>1) printf("AnalysisTaskZDCpAcalib::UserCreateOutputObjects() \n");

  if (fUseSpecialOutput) OpenFile(1);

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();
  //fOutput->SetName("output");

    fCentralityTree = new TTree("fCentralityTree", "Centrality vs. multiplicity tree");
    //
    fCentralityTree->Branch("trigClass",&fTrigClass,"trigClass/C");
    fCentralityTree->Branch("eventSelected",&fIsEventSelected,"eventSelected/O");
    fCentralityTree->Branch("isV0ATriggered", &fIsV0ATriggered,"isV0ATriggered/I");
    fCentralityTree->Branch("isV0CTriggered", &fIsV0CTriggered,"isV0CTriggered/I");
    fCentralityTree->Branch("znctower", fZNCtower, "znctower[5]/F");
    fCentralityTree->Branch("znatower", fZNAtower, "znatower[5]/F");
    fCentralityTree->Branch("znctowerLG", fZNCtowerLG, "znctowerLG[5]/F");
    fCentralityTree->Branch("znatowerLG", fZNAtowerLG, "znatowerLG[5]/F");
    fCentralityTree->Branch("zpctower", fZNCtower, "znctower[5]/F");
    fCentralityTree->Branch("zpatower", fZNAtower, "znatower[5]/F");
    fCentralityTree->Branch("zpctowerLG", fZNCtowerLG, "znctowerLG[5]/F");
    fCentralityTree->Branch("zpatowerLG", fZNAtowerLG, "znatowerLG[5]/F");
    fCentralityTree->Branch("tdcZNC", fZNCtdc, "tdcZNC[4]/I");
    fCentralityTree->Branch("tdcZNA", fZNAtdc, "tdcZNA[4]/I");
    fCentralityTree->Branch("tdcZPC", fZPCtdc, "tdcZPC[4]/I");
    fCentralityTree->Branch("tdcZPA", fZPAtdc, "tdcZPA[4]/I");
      
    fOutput->Add(fCentralityTree);      
    PostData(1, fOutput);
  
}

//________________________________________________________________________
void AliAnalysisTaskZDCpAcalib::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  if(fDebug>1) printf(" **** AliAnalysisTaskZDCpAcalib::UserExec() \n");
  
  if (!InputEvent()) {
    Printf("ERROR: InputEvent not available");
    return;
  }

  if(fAnalysisInput.CompareTo("ESD")==0){
      
      //printf(" \t ***Analizing ESD ev. %d\n",fNev);
      
      AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
      if(!esd) return;
      
      // Select PHYSICS events (type=7, for data)
      if(!fIsMCInput && esd->GetEventType()!=7) return; 
      
      // ***** Trigger selection
      TString triggerClass = esd->GetFiredTriggerClasses();
      sprintf(fTrigClass,"%s",triggerClass.Data());
      
      // use response of AliPhysicsSelection
      fIsEventSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCINT5);       

          
      AliESDVZERO *vzeroAOD = esd->GetVZEROData();
      fIsV0ATriggered = vzeroAOD->GetV0ADecision();
      fIsV0CTriggered = vzeroAOD->GetV0CDecision();
        
      AliESDZDC *esdZDC = esd->GetESDZDC();
       
      const Double_t * towZNC = esdZDC->GetZN1TowerEnergy();
      const Double_t * towZNA = esdZDC->GetZN2TowerEnergy();
      const Double_t * towZPC = esdZDC->GetZP1TowerEnergy();
      const Double_t * towZPA = esdZDC->GetZP2TowerEnergy();
      //
      const Double_t * towZNCLG = esdZDC->GetZN1TowerEnergyLR();
      const Double_t * towZNALG = esdZDC->GetZN2TowerEnergyLR();
      const Double_t * towZPCLG = esdZDC->GetZP1TowerEnergyLR();
      const Double_t * towZPALG = esdZDC->GetZP2TowerEnergyLR();
      //
      for(Int_t it=0; it<5; it++){
         fZNCtower[it] = (Float_t) (towZNC[it]);
         fZNAtower[it] = (Float_t) (towZNA[it]); 
         fZNCtowerLG[it] = (Float_t) (towZNCLG[it]);
         fZNAtowerLG[it] = (Float_t) (towZNALG[it]); 
         fZPCtower[it] = (Float_t) (towZPC[it]);
         fZPAtower[it] = (Float_t) (towZPA[it]); 
         fZPCtowerLG[it] = (Float_t) (towZPCLG[it]);
         fZPAtowerLG[it] = (Float_t) (towZPALG[it]); 
      }

      for(Int_t i=0; i<4; i++){
	   fZNCtdc[i] = esdZDC->GetZDCTDCCorrected(10, i);
	   fZPCtdc[i] = esdZDC->GetZDCTDCCorrected(11, i);
	   fZNAtdc[i] = esdZDC->GetZDCTDCCorrected(12, i);
	   fZPAtdc[i] = esdZDC->GetZDCTDCCorrected(13, i);
      }      
  }   
  
  
  fCentralityTree->Fill();
 
  PostData(1, fOutput);
  
   
}



//________________________________________________________________________
void AliAnalysisTaskZDCpAcalib::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  printf(" **** AliAnalysisTaskZDCpAcalib::Terminate() \n");
}
