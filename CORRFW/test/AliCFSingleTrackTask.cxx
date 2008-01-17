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
// Example of task running on AliEn (CAF?)
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------


#ifndef ALICFSINGLETRACKTASK_CXX
#define ALICFSINGLETRACKTASK_CXX
#include <TROOT.h>
#include <TInterpreter.h>

#include "AliCFSingleTrackTask.h"
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
ClassImp(AliCFSingleTrackTask)
//__________________________________________________________________________
AliCFSingleTrackTask::AliCFSingleTrackTask() :
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0)
{
//Defual ctor
}
//___________________________________________________________________________
AliCFSingleTrackTask::AliCFSingleTrackTask(const Char_t* name) :
  AliAnalysisTask(name,"AliCFSingleTrackTask"),
  fChain(0x0),
  fESD(0x0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFSingleTrackTask","Calling Constructor");
  DefineInput (0,TChain::Class());
  DefineOutput(0,TH1I::Class());
  DefineOutput(1,AliCFContainer::Class());
  DefineOutput(2,TList::Class());
}

//___________________________________________________________________________
AliCFSingleTrackTask& AliCFSingleTrackTask::operator=(const AliCFSingleTrackTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTask::operator=(c) ;
    fChain      = c.fChain;
    fESD        = c.fESD;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList ;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFSingleTrackTask::AliCFSingleTrackTask(const AliCFSingleTrackTask& c) :
  AliAnalysisTask(c),
  fChain(c.fChain),
  fESD(c.fESD),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fHistEventsProcessed(c.fHistEventsProcessed)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFSingleTrackTask::~AliCFSingleTrackTask() {
  //
  //destructor
  //
  Info("~AliCFSingleTrackTask","Calling Destructor");
  if (fChain)               delete fChain ;
  if (fESD)                 delete fESD ;
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fQAHistList) {fQAHistList->Clear(); delete fQAHistList;}
}

//___________________________________________________________________________

void AliCFSingleTrackTask::Init()
{

}
//_________________________________________________
void AliCFSingleTrackTask::Exec(Option_t *)
{
  //
  // Main loop function
  //
  Info("Exec","") ;
  // Get the mc truth
  AliMCEventHandler* mcTruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());

  if (!mcTruth) Error("Exec","NO MC INFO FOUND... EXITING\n");

  // transform possible old AliESD into AliESDEvent

  if (fESD->GetAliESDOld()) fESD->CopyFromOldESD(); //transition to new ESD format

  //pass the MC evt handler to the cuts that need it 

  fCFManager->SetEventInfo(mcTruth);

  // Get the MC event 
  AliMCEvent* mcEvent = mcTruth->MCEvent();

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  for (Int_t ipart=0; ipart<mcEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = mcEvent->GetTrack(ipart);

    //check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;

    containerInput[0] = (Float_t)mcPart->Pt();
    containerInput[1] = mcPart->Eta() ;
    //fill the container for Gen-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
    
    //check the Acceptance-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;
    //fill the container for Acceptance-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
  }    

  //Now go to rec level
  for (Int_t iTrack = 0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {

    AliESDtrack* track = fESD->GetTrack(iTrack);
    
    fCFManager->FillQABeforeParticleCuts(AliCFManager::kPartRecCuts,track);  
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,track)) continue;
    fCFManager->FillQAAfterParticleCuts(AliCFManager::kPartRecCuts,track);

    // is track associated to particle ?
    if (track->GetLabel()<0) continue;
    AliMCParticle *mcPart  = mcEvent->GetTrack(track->GetLabel());
    
    // check if this track was part of the signal
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue; 
    
    //fill the container
    Double_t mom[3];
    track->GetPxPyPz(mom);
    Double_t pt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    containerInput[0] = pt ;
    containerInput[1] = track->Eta();
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,track)) continue ;
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepSelected);
  }
  
  fHistEventsProcessed->Fill(0);
  PostData(0,fHistEventsProcessed) ;
  PostData(1,fCFManager->GetParticleContainer()) ;
  PostData(2,fQAHistList) ;
}


//___________________________________________________________________________
void AliCFSingleTrackTask::Terminate(Option_t*)
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

//   TCanvas * c1 =new TCanvas("c1","",600,400);
  c->cd(1);
  fCFManager->GetParticleContainer()->ShowProjection(0,0)->Draw("p");
  c->cd(2);
  fCFManager->GetParticleContainer()->ShowProjection(0,1)->Draw("p");
  c->cd(3);
  fCFManager->GetParticleContainer()->ShowProjection(0,2)->Draw("p");
  c->cd(4);
  fCFManager->GetParticleContainer()->ShowProjection(0,3)->Draw("p");

//   TCanvas * c2 =new TCanvas("c2","",600,400);
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
void AliCFSingleTrackTask::ConnectInputData(Option_t *) {
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
void AliCFSingleTrackTask::CreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("CreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #0
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;

  //slot #2
  fQAHistList = new TList();
  fCFManager->InitQAHistos();
  fCFManager->AddQAHistosToList(fQAHistList);
}

#endif
