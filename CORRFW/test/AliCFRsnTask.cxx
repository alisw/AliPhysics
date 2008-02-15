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

//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------


#include <TROOT.h>
#include <TInterpreter.h>

#include "AliCFRsnTask.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "TChain.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliRsnDaughter.h"
#include "AliCFPair.h"
#include <memory>

//__________________________________________________________________________
AliCFRsnTask::AliCFRsnTask() :
  fRsnPDG(313),
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //Defual ctor
}
//___________________________________________________________________________
AliCFRsnTask::AliCFRsnTask(const Char_t* name) :
  AliAnalysisTask(name,"AliCFRsnTask"),
  fRsnPDG(313),
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFRsnTask","Calling Constructor");
  DefineInput (0,TChain::Class());
  DefineOutput(0,TH1I::Class());
  DefineOutput(1,AliCFContainer::Class());
  //   DefineOutput(2,TList::Class());
}

//___________________________________________________________________________
AliCFRsnTask& AliCFRsnTask::operator=(const AliCFRsnTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTask::operator=(c) ;
    fRsnPDG     = c.fRsnPDG;
    fChain      = c.fChain;
    fESD        = c.fESD;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFRsnTask::AliCFRsnTask(const AliCFRsnTask& c) :
  AliAnalysisTask(c),
  fRsnPDG(c.fRsnPDG),
  fChain(c.fChain),
  fESD(c.fESD),
  fCFManager(c.fCFManager),
  fHistEventsProcessed(c.fHistEventsProcessed)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFRsnTask::~AliCFRsnTask() {
  //
  //destructor
  //
  Info("~AliCFRsnTask","Calling Destructor");
  if (fChain)               delete fChain ;
  if (fESD)                 delete fESD ;
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//___________________________________________________________________________

void AliCFRsnTask::Init()
{

}
//_________________________________________________
void AliCFRsnTask::Exec(Option_t *)
{
  //
  // Main loop function
  //
  Info("Exec","") ;
  // Get the mc truth
  AliMCEventHandler* mcTruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

  if (!mcTruth) Error("Exec","NO MC INFO FOUND... EXITING");

  // transform possible old AliESD into AliESDEvent

  if (fESD->GetAliESDOld()) fESD->CopyFromOldESD(); //transition to new ESD format

  //pass the MC evt handler to the cuts that need it 

  fCFManager->SetEventInfo(mcTruth);

  // Get the MC event 
  AliMCEvent* mcEvent = mcTruth->MCEvent();
  AliStack*   stack   = mcEvent->Stack();

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  Info("Exec","Looping on MC event");
  for (Int_t ipart=0; ipart<stack->GetNprimary(); ipart++) { 
    AliMCParticle *mcPart  = mcEvent->GetTrack(ipart);

    //check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;
    containerInput[0] = mcPart->Pt();
    containerInput[1] = mcPart->Eta() ;
    //fill the container for Gen-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
    
    //check the Acceptance-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;
    //fill the container for Acceptance-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
  }    


  //Now go to rec level
  Info("Exec","Looping on ESD event");

  //SET THE ESD AS EVENT INFO IN RECONSTRUCTION CUTS
  TObjArray* fCutsReco = fCFManager->GetParticleCutsList(AliCFManager::kPartRecCuts);
  TObjArray* fCutsPID = fCFManager->GetParticleCutsList(AliCFManager::kPartSelCuts);
  TObjArrayIter iter1(fCutsReco);
  TObjArrayIter iter2(fCutsPID);
  AliCFCutBase *cut = 0;
  while ( (cut = (AliCFCutBase*)iter1.Next()) ) {
    cut->SetEvtInfo(fESD);
  }
  while ( (cut = (AliCFCutBase*)iter2.Next()) ) {
    cut->SetEvtInfo(fESD);
  }
  
  for (Int_t iTrack1 = 0; iTrack1<fESD->GetNumberOfTracks(); iTrack1++) {
    AliESDtrack* esdTrack1 = fESD->GetTrack(iTrack1);
    //track1 is negative
    if (esdTrack1->Charge()>=0) continue;
    Int_t esdLabel1 = esdTrack1->GetLabel();
    if (esdLabel1<0) continue;

    for (Int_t iTrack2 = 0; iTrack2<fESD->GetNumberOfTracks(); iTrack2++) {
      AliESDtrack* esdTrack2 = fESD->GetTrack(iTrack2);
      //track2 is positive
      if (esdTrack2->Charge()<=0) continue;
      Int_t esdLabel2 = esdTrack2->GetLabel();
      if (esdLabel2<0) continue;
	
      // copy ESDtrack data into RsnDaughter (and make Bayesian PID)
      //if problem with copy constructor, use the following special command to destroy the pointer when exiting scope
      //std::auto_ptr<AliRsnDaughter> track1(AliRsnDaughter::Adopt(esdTrack1,iTrack1));

      AliRsnDaughter* tmp1=AliRsnDaughter::Adopt(esdTrack1,iTrack1);
      AliRsnDaughter* tmp2=AliRsnDaughter::Adopt(esdTrack2,iTrack2);
      AliRsnDaughter track1(*tmp1);
      AliRsnDaughter track2(*tmp2);
      delete tmp1;
      delete tmp2;

      TParticle *part1 = stack->Particle(esdLabel1);
      track1.SetTruePDG(part1->GetPdgCode());
      Int_t mother1 = part1->GetFirstMother();
      track1.SetMother(mother1);
      if (mother1 >= 0) {
	TParticle *mum = stack->Particle(mother1);
	track1.SetMotherPDG(mum->GetPdgCode());
      }
      
      TParticle *part2 = stack->Particle(esdLabel2);
      track2.SetTruePDG(part2->GetPdgCode());
      Int_t mother2 = part2->GetFirstMother();
      track2.SetMother(mother2);
      if (mother2 >= 0) {
	TParticle *mum = stack->Particle(mother2);
	track2.SetMotherPDG(mum->GetPdgCode());
      }
	
      //make a mother resonance from the 2 candidate daughters
      AliRsnDaughter rsn = AliRsnDaughter::Sum(track1,track2) ;
      AliCFPair pair(esdTrack1,esdTrack2);

      //check if true resonance
      if (rsn.GetMotherPDG() != fRsnPDG) continue;
      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,&pair)) continue;

      //check if associated MC resonance passes the cuts
      Int_t motherLabel=rsn.GetLabel();
      if (motherLabel<0) continue ;
      AliMCParticle* mcRsn = mcEvent->GetTrack(motherLabel);
      if (!mcRsn) continue;
      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcRsn)) continue; 
    
      //fill the container
      containerInput[0] = rsn.GetPt() ;
      containerInput[1] = GetRapidity(rsn.GetEnergy(),rsn.GetPz());
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   

      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,&pair)) continue ;
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepSelected);
    }
  }
    
  fHistEventsProcessed->Fill(0);
  PostData(0,fHistEventsProcessed) ;
  PostData(1,fCFManager->GetParticleContainer()) ;
  
  //   TList * list = new TList();
  //   fCFManager->AddQAHistosToList(list);
  //   PostData(2,list) ;
}


//___________________________________________________________________________
void AliCFRsnTask::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  AliAnalysisTask::Terminate();


  Double_t max1 = fCFManager->GetParticleContainer()->ShowProjection(0,0)->GetMaximum();
  Double_t max2 = fCFManager->GetParticleContainer()->ShowProjection(1,0)->GetMaximum();

  fCFManager->GetParticleContainer()->ShowProjection(0,0)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->GetYaxis()->SetRangeUser(0,max1*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->GetYaxis()->SetRangeUser(0,max1*1.2);

  fCFManager->GetParticleContainer()->ShowProjection(1,0)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->GetYaxis()->SetRangeUser(0,max2*1.2);
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->GetYaxis()->SetRangeUser(0,max2*1.2);

  fCFManager->GetParticleContainer()->ShowProjection(0,0)->SetMarkerStyle(23) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->SetMarkerStyle(24) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->SetMarkerStyle(25) ;
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->SetMarkerStyle(26) ;

  fCFManager->GetParticleContainer()->ShowProjection(1,0)->SetMarkerStyle(23) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->SetMarkerStyle(24) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->SetMarkerStyle(25) ;
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->SetMarkerStyle(26) ;

  TCanvas * c =new TCanvas("c","",1400,800);
  c->Divide(4,2);

  c->cd(1);
  fCFManager->GetParticleContainer()->ShowProjection(0,0)->Draw("p");
  c->cd(2);
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->Draw("p");
  c->cd(3);
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->Draw("p");
  c->cd(4);
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->Draw("p");
  c->cd(5);
  fCFManager->GetParticleContainer()->ShowProjection(1,0)->Draw("p");
  c->cd(6);
  fCFManager->GetParticleContainer()->ShowProjection(1,1)->Draw("p");
  c->cd(7);
  fCFManager->GetParticleContainer()->ShowProjection(1,2)->Draw("p");
  c->cd(8);
  fCFManager->GetParticleContainer()->ShowProjection(1,3)->Draw("p");

  c->SaveAs("plots.eps");

  delete fHistEventsProcessed ;
}

//___________________________________________________________________________
void AliCFRsnTask::ConnectInputData(Option_t *) {
  //
  // Initialize branches.
  //
  Info("ConnectInputData","ConnectInputData of task %s\n",GetName());

  fChain = (TChain*)GetInputData(0);
  fChain->SetBranchStatus("*FMD*",0);
  fChain->SetBranchStatus("*CaloClusters*",0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//___________________________________________________________________________
void AliCFRsnTask::CreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("CreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #0
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

//___________________________________________________________________________
Double_t AliCFRsnTask::GetRapidity(Double_t energy, Double_t pz) {
  if (energy == pz || energy == -pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  if (energy < pz) {
    printf("GetRapidity : ERROR : rapidity for 4-vector with E = Pz -- infinite result");
    return 999;
  }
  Double_t y = 0.5 * log((energy + pz) / (energy - pz));
  return y;
} 


