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
// Task to prepare efficiency and response matrix for unfolding
//-----------------------------------------------------------------------
// Author : R. Vernet, INFN - Catania (it)
//-----------------------------------------------------------------------


#include "TStyle.h"
#include "AliCFTaskForUnfolding.h"
#include "TCanvas.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "THnSparse.h"
#include "TH2D.h"

ClassImp(AliCFTaskForUnfolding)

//__________________________________________________________________________
AliCFTaskForUnfolding::AliCFTaskForUnfolding() :
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliCFTaskForUnfolding::AliCFTaskForUnfolding(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fCFManager(0x0),
  fHistEventsProcessed(0x0),
  fCorrelation(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFTaskForUnfolding","Calling Constructor");

  /*
    DefineInput(0) and DefineOutput(0)
    are taken care of by AliAnalysisTaskSE constructor
  */
  DefineOutput(1,TH1I::Class());
  DefineOutput(2,AliCFContainer::Class());
  DefineOutput(3,THnSparseD::Class());
}

//___________________________________________________________________________
AliCFTaskForUnfolding::~AliCFTaskForUnfolding() {
  //
  //destructor
  //
  Info("~AliCFTaskForUnfolding","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
  if (fHistEventsProcessed) delete fHistEventsProcessed ;
  if (fCorrelation)         delete fCorrelation ;
}

//_________________________________________________
void AliCFTaskForUnfolding::UserExec(Option_t *)
{
  //
  // Main loop function
  //
  AliInfo("") ;

  AliVEvent*    fEvent = fInputEvent ;
  AliVParticle* track ;
  
  if (!fEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  if (!fMCEvent) Error("UserExec","NO MC INFO FOUND");
  
  //pass the MC evt handler to the cuts that need it 
  fCFManager->SetMCEventInfo(fMCEvent);

  // MC-event selection
  Double_t containerInput[2] ;
        
  //loop on the MC event
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);
    
    if (!fCFManager->CheckParticleCuts(0,mcPart)) continue;
    containerInput[0] = (Float_t)mcPart->Pt();
    containerInput[1] = (Float_t)mcPart->Eta();
    fCFManager->GetParticleContainer()->Fill(containerInput,0);
  }    

  //Now go to rec level
  for (Int_t iTrack = 0; iTrack<fEvent->GetNumberOfTracks(); iTrack++) {
    
    track = fEvent->GetTrack(iTrack);
    if (!fCFManager->CheckParticleCuts(1,track)) continue;
    
    Int_t label = track->GetLabel();
    if (label<0) continue;
    AliMCParticle* mcPart = (AliMCParticle*)fMCEvent->GetTrack(label);
    // check if this track was part of the signal
    if (!fCFManager->CheckParticleCuts(0,mcPart)) continue;

    //fill the container
    Double_t mom[3];
    track->PxPyPz(mom);
    Double_t pt=TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    containerInput[0] = pt;
    containerInput[1] = track->Eta();
    fCFManager->GetParticleContainer()->Fill(containerInput,1) ;
    containerInput[0] = mcPart->Pt();
    containerInput[1] = mcPart->Eta();
    fCFManager->GetParticleContainer()->Fill(containerInput,2);
    
    Double_t fill[4]; //fill response matrix
    // dimensions 0&1 : pt,eta (Rec)
    fill[0] = pt ;
    fill[1] = track->Eta();
    // dimensions 2&3 : pt,eta (MC)
    fill[2] = mcPart->Pt();
    fill[3] = mcPart->Eta();
    fCorrelation->Fill(fill);
  }
  
  fHistEventsProcessed->Fill(0);

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fCorrelation) ;
}


//___________________________________________________________________________
void AliCFTaskForUnfolding::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  Info("Terminate","");
  AliAnalysisTaskSE::Terminate();

  gStyle->SetPalette(1);

  //draw some example plots....

  AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
  TH2D* h00 =   cont->ShowProjection(0,1,0) ;
  TH2D* h01 =   cont->ShowProjection(0,1,1) ;
  THnSparseD* hcorr = dynamic_cast<THnSparseD*> (GetOutputData(3));

  TCanvas * c =new TCanvas("c","",800,400);
  c->Divide(2,1);
  c->cd(1);
  h00->Draw("text");
  c->cd(2);
  h01->Draw("text");
  c->SaveAs("spectra.eps");

  TCanvas * c2 =new TCanvas("c2","",800,400);
  c2->Divide(2,1);
  c2->cd(1);
  hcorr->Projection(0,2)->Draw("text");
  c2->cd(2);
  hcorr->Projection(1,3)->Draw("text");
  c2->SaveAs("correlation.eps");
}


//___________________________________________________________________________
void AliCFTaskForUnfolding::UserCreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());

  //slot #1
  OpenFile(1);
  fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;

//   OpenFile(2);
//   OpenFile(3);
}

