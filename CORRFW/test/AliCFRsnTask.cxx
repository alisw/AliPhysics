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


#include "AliCFRsnTask.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "TChain.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliRsnDaughter.h"
#include "AliCFPair.h"
#include "AliRsnParticle.h"

//__________________________________________________________________________
AliCFRsnTask::AliCFRsnTask() :
  AliAnalysisTaskSE(),
  fRsnPDG(0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFRsnTask::AliCFRsnTask(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fRsnPDG(0),
  fCFManager(0x0),
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFRsnTask","Calling Constructor");
  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
  //   DefineOutput(3,TList::Class());
}

//___________________________________________________________________________
AliCFRsnTask& AliCFRsnTask::operator=(const AliCFRsnTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
    fRsnPDG     = c.fRsnPDG;
    fCFManager  = c.fCFManager;
    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliCFRsnTask::AliCFRsnTask(const AliCFRsnTask& c) :
  AliAnalysisTaskSE(c),
  fRsnPDG(c.fRsnPDG),
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
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//_________________________________________________
void AliCFRsnTask::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  Info("UserExec","") ;

  AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (!fESD) {
    Error("UserExec","NO ESD FOUND!");
    return;
  }

  if (!fMCEvent) Error("UserExec","NO MC INFO FOUND!");
  fCFManager->SetEventInfo(fMCEvent);

  AliStack*   stack   = fMCEvent->Stack();

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  Info("UserExec","Looping on MC event");
  for (Int_t ipart=0; ipart<stack->GetNprimary(); ipart++) { 
    AliMCParticle *mcPart  = fMCEvent->GetTrack(ipart);

    //check the MC-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;
    containerInput[0] = mcPart->Pt();
    containerInput[1] = mcPart->Y() ;
    //fill the container for Gen-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
    
    //check the Acceptance-level cuts
    if (!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)) continue;
    //fill the container for Acceptance-level selection
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
  }    


  //Now go to rec level
  Info("UserExec","Looping on ESD event");

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

  // Loop on negative tracks
  for (Int_t iTrack1 = 0; iTrack1<fESD->GetNumberOfTracks(); iTrack1++) {
    AliESDtrack* esdTrack1 = fESD->GetTrack(iTrack1);
    //track1 is negative
    if (esdTrack1->Charge()>=0) continue;
    Int_t esdLabel1 = esdTrack1->GetLabel();
    if (esdLabel1<0) continue;

    //Loop on positive tracks
    for (Int_t iTrack2 = 0; iTrack2<fESD->GetNumberOfTracks(); iTrack2++) {
      AliESDtrack* esdTrack2 = fESD->GetTrack(iTrack2);
      //track2 is positive
      if (esdTrack2->Charge()<=0) continue;
      Int_t esdLabel2 = esdTrack2->GetLabel();
      if (esdLabel2<0) continue;
	
      //Create Resonance daughter objects
      AliRsnDaughter* tmp1 = new AliRsnDaughter(esdTrack1);
      AliRsnDaughter* tmp2 = new AliRsnDaughter(esdTrack2);
      AliRsnDaughter track1(*tmp1);
      AliRsnDaughter track2(*tmp2);
      delete tmp1;
      delete tmp2;

      //Set MC information to resonance daughters
      TParticle *part1 = stack->Particle(esdLabel1);
      track1.InitParticle(part1);
      track1.GetParticle()->SetPDG(part1->GetPdgCode());

      Int_t mother1 = part1->GetFirstMother();
      track1.GetParticle()->SetMother(mother1);
      if (mother1 >= 0) {
	TParticle *mum = stack->Particle(mother1);
	track1.GetParticle()->SetMotherPDG(mum->GetPdgCode());
      }
      
      TParticle *part2 = stack->Particle(esdLabel2);
      track2.InitParticle(part2);
      track2.GetParticle()->SetPDG(part2->GetPdgCode());

      Int_t mother2 = part2->GetFirstMother();
      track2.GetParticle()->SetMother(mother2);
      if (mother2 >= 0) {
	TParticle *mum = stack->Particle(mother2);
	track2.GetParticle()->SetMotherPDG(mum->GetPdgCode());
      }
	
      //make a mother resonance from the 2 candidate daughters
      AliRsnDaughter rsn = AliRsnDaughter::Sum(track1,track2) ;
      AliCFPair pair(esdTrack1,esdTrack2); // This object is used for cuts (to be replaced)

      //check if true resonance
      if (rsn.GetParticle()->PDG() != fRsnPDG) continue;
      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartRecCuts,&pair)) continue;

      //check if associated MC resonance passes the cuts
      Int_t motherLabel=rsn.Label();
      if (motherLabel<0) continue ;
      AliMCParticle* mcRsn = fMCEvent->GetTrack(motherLabel);
      if (!mcRsn) continue;
      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcRsn)) continue; 

      //fill the container
      containerInput[0] = rsn.Pt() ;
      containerInput[1] = GetRapidity(rsn.E(),rsn.Pz());
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   

      if (!fCFManager->CheckParticleCuts(AliCFManager::kPartSelCuts,&pair)) continue ;
      fCFManager->GetParticleContainer()->Fill(containerInput,kStepSelected);
    }
  }
    
  fHistEventsProcessed->Fill(0);
  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  
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
void AliCFRsnTask::UserCreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());

  //slot #1
  OpenFile(1);
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


