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
// Example of task (running locally, on AliEn and CAF),
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------


#ifndef ALIPWG4HighPtSpectra_CXX
#define ALIPWG4HighPtSpectra_CXX

#include "AliPWG4HighPtSpectra.h"

#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliCFContainer.h"
#include "TChain.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisHelperJetTasks.h"

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtSpectra)

//__________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra() : AliAnalysisTask("AliPWG4HighPtSpectra", ""), 
  fReadAODData(0),
  fCFManager(0x0),
  fESD(0),
  fTrackCuts(0),
  fTrigger(0),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0)
  //  fHistEventsProcessed(0x0)
{
  //
  //Default ctor
  //
}
//___________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra(const Char_t* name) :
  AliAnalysisTask(name,""),
  fReadAODData(0),
  fCFManager(0x0),
  fESD(0),
  fTrackCuts(),//new AliESDtrackCuts),
  fTrigger(0),
  fHistList(0),
  fNEventAll(0),
  fNEventSel(0)
  //  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  AliDebug(2,Form("AliPWG4HighPtQAMC","Calling Constructor"));
  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  DefineOutput(0,TList::Class());
  DefineOutput(1,AliCFContainer::Class());
}

//___________________________________________________________________________
AliPWG4HighPtSpectra& AliPWG4HighPtSpectra::operator=(const AliPWG4HighPtSpectra& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTask::operator=(c) ;
    fReadAODData = c.fReadAODData ;
    fCFManager  = c.fCFManager;
    fHistList = c.fHistList;
    fNEventAll = c.fNEventAll;
    fNEventSel = c.fNEventSel;
    //    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//___________________________________________________________________________
AliPWG4HighPtSpectra::AliPWG4HighPtSpectra(const AliPWG4HighPtSpectra& c) :
  AliAnalysisTask(c),
  fReadAODData(c.fReadAODData),
  fCFManager(c.fCFManager),
  fESD(c.fESD),
  fTrackCuts(c.fTrackCuts),
  fTrigger(c.fTrigger),
  fHistList(c.fHistList),
  fNEventAll(c.fNEventAll),
  fNEventSel(c.fNEventSel)
  //  fHistEventsProcessed(c.fHistEventsProcessed)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliPWG4HighPtSpectra::~AliPWG4HighPtSpectra() {
  //
  //destructor
  //
  Info("~AliPWG4HighPtSpectra","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
  if (fNEventAll) delete fNEventAll ;
  if (fNEventSel) delete fNEventSel ;
}
//________________________________________________________________________
void AliPWG4HighPtSpectra::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));
  //  printf(">> AliPWG4HighPtSpectra::ConnectInputData \n");

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliDebug(2,Form("ERROR: Could not read chain from input slot 0"));
  } else {
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
    } else
      fESD = esdH->GetEvent();
  }
 
}
//_________________________________________________
void AliPWG4HighPtSpectra::Exec(Option_t *)//UserExec(Option_t *)
{
  //
  // Main loop function
  //
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::Exec \n"));  

  // All events without selection
  fNEventAll->Fill(0.);

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    return;
  }

  //Trigger selection
  AliAnalysisHelperJetTasks::Trigger trig;
  trig = (AliAnalysisHelperJetTasks::Trigger)fTrigger;
  if (AliAnalysisHelperJetTasks::IsTriggerFired(fESD,trig)){
    AliDebug(2,Form(" Trigger Selection: event ACCEPTED ... "));
  }else{
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    PostData(0,fHistList);
    PostData(1,fCFManager->GetParticleContainer());
    return;
  } 
  //  if(!fESD->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL") || !fESD->IsTriggerClassFired("CSMBB-ABCE-NOPF-ALL")) return;

  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  AliStack* stack = 0x0;
  AliMCEvent* mcEvent = 0x0;

  if(eventHandler) {
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      AliDebug(2,Form("ERROR: Could not retrieve MC event"));
      return;
    }
    
    AliDebug(2,Form("MC particles: %d", mcEvent->GetNumberOfTracks()));
    
    stack = mcEvent->Stack();                //Particles Stack
    
    AliDebug(2,Form("MC particles stack: %d", stack->GetNtrack()));
  }
  
  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  // Need vertex cut
  if (vtx->GetNContributors() < 2){
    PostData(0,fHistList);
    PostData(1,fCFManager->GetParticleContainer());
    return;
  }
  double primVtx[3];
  vtx->GetXYZ(primVtx);
  if(TMath::Abs(primVtx[0]>1. )|| TMath::Abs(primVtx[1]>1.) || TMath::Abs(primVtx[2]>10.)){
    PostData(0,fHistList);
    PostData(1,fCFManager->GetParticleContainer());
    return;
  }
  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));

  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));

  // Selected events for analysis
  fNEventSel->Fill(0.);

  Double_t containerInputRec[1] ;
  Double_t containerInputTPConly[1] ;
  Double_t containerInputMC[1] ;
  //Now go to rec level
  for (Int_t iTrack = 0; iTrack<nTracks; iTrack++) 
    {   
      if(!fESD->GetTrack(iTrack) ) continue;
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if(!(AliExternalTrackParam *)track->GetTPCInnerParam()) continue;
      AliExternalTrackParam *trackTPC = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!track || !trackTPC) continue;

    
      //fill the container
      containerInputRec[0] = track->Pt();
      containerInputTPConly[0] = trackTPC->Pt();

      if (fTrackCuts->AcceptTrack(track)) {
	fCFManager->GetParticleContainer()->Fill(containerInputRec,kStepReconstructed);
	fCFManager->GetParticleContainer()->Fill(containerInputTPConly,kStepReconstructedTPCOnly);
	
	//Only fill the secondary particle container if MC information is available
	if(eventHandler) {
	  Int_t label = TMath::Abs(track->GetLabel());
	  TParticle *particle = stack->Particle(label) ;
	  if(!particle) continue;
	  containerInputMC[0] = particle->Pt();      
	  fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepReconstructedMC);
	  if (!stack->IsPhysicalPrimary(label) ) {
	    fCFManager->GetParticleContainer()->Fill(containerInputRec,kStepSecondaries);
	  }
	}
      }
      
    }

  if(eventHandler) {
    for(int iPart = 1; iPart<(mcEvent->GetNumberOfTracks()); iPart++)//stack->GetNprimary();
      {
	AliMCParticle *mcPart  = (AliMCParticle*)mcEvent->GetTrack(iPart);
	if(!mcPart) continue;
	
	//fill the container
	containerInputMC[0] = mcPart->Pt();
	
	if (!fCFManager->CheckParticleCuts(3,mcPart)) continue ;
	
	int counter;
	
	Float_t trackLengthTPC = mcPart->GetTPCTrackLength(fESD->GetMagneticField(),0.1,counter,3.0);
	
	if(trackLengthTPC>80.) fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepMCtrackable) ;
	
      }
  }

   PostData(0,fHistList);
   PostData(1,fCFManager->GetParticleContainer());

}


//___________________________________________________________________________
void AliPWG4HighPtSpectra::Terminate(Option_t*)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.


}

//___________________________________________________________________________
void AliPWG4HighPtSpectra::CreateOutputObjects() {
  //HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
  //TO BE SET BEFORE THE EXECUTION OF THE TASK
  //
  AliDebug(2,Form("CreateOutputObjects","CreateOutputObjects of task %s", GetName()));

  //slot #1
  OpenFile(0);
  fHistList = new TList();
  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);

}

#endif
