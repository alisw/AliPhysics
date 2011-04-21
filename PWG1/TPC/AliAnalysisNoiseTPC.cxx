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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Class to test TPC warm-up LASER EVENTS rejection.                     //
//                                                                       //
// Author: Alexander Kalweit (GSI)                                       //
// Modified: Jacek Otwinowski (GSI)                                      //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliRecoParam.h"

#include "AliLog.h"

#include "AliAnalysisNoiseTPC.h"


ClassImp(AliAnalysisNoiseTPC)

//________________________________________________________________________
AliAnalysisNoiseTPC::AliAnalysisNoiseTPC() 
: AliAnalysisTaskSE("TaskChargedHadron"), fESD(0), fListHist(0), fESDtrackCuts(0),
  fHistNoiseTracks(0),fTimeBins(0),fTimeStart(0),fTimeEnd(0)
{
  // default Constructor
}


//________________________________________________________________________
AliAnalysisNoiseTPC::AliAnalysisNoiseTPC(const char *name,UInt_t StartTime, UInt_t EndTime, Int_t deltaTime) 
  : AliAnalysisTaskSE(name), fESD(0), fListHist(0), fESDtrackCuts(0),
    fHistNoiseTracks(0), fTimeBins(0), fTimeStart(StartTime), fTimeEnd(EndTime)
{
  //
  // standard constructur which should be used
  //
  Printf("*** CONSTRUCTOR CALLED ****");
  if(deltaTime)
  {
    fTimeBins = TMath::Nint((fTimeEnd-fTimeStart)/deltaTime);  
    printf("fTimeBins %d \n",fTimeBins);
  } else {
    printf("deltaTime is 0 \n"); 
  }

  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());

}




//________________________________________________________________________
void AliAnalysisNoiseTPC::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once

  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(80);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(2);
  //
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kFALSE);
  fESDtrackCuts->SetEtaRange(-0.005,+0.005);
  fESDtrackCuts->SetPtRange(10., 1e10);


  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  /*
  UInt_t StartTime = 1272672000; // 1st of May
  UInt_t EndTime = 1288569600; // 1st of Nov
  Double_t deltaTime = EndTime - StartTime;
  Int_t timeBins = TMath::Nint(deltaTime/(15.*60)); // every 15min
  */

  // event number in file, # noise tracks, isSelectMB, isSelectWarm, has vtx., time
  Int_t    binsHistNoise[6] = { 20000,  100,    2, 2,   2, fTimeBins};
  Double_t xminHistNoise[6] = { -0.5, -0.5, -0.5, -0.5,-0.5, fTimeStart};
  Double_t xmaxHistNoise[6] = {19999.5, 99.5,  1.5,  1.5, 1.5, fTimeEnd};
  fHistNoiseTracks = new THnSparseS("fHistNoiseTracks","noise tracks:  ev in file, # noise tracks, isSelectMB, isSelectWarm, has vtx,time",6,binsHistNoise,xminHistNoise,xmaxHistNoise);
  //
  fListHist->Add(fHistNoiseTracks);
  //
   
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliAnalysisNoiseTPC::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  AliLog::SetGlobalLogLevel(AliLog::kError);
  //
  // Check Monte Carlo information and other access first:
  //
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    //Printf("ERROR: fESD not available");
    return;
  }
  
  if (!fESDtrackCuts) {
    Printf("ERROR: fESDtrackCuts not available");
    return;
  }
  
  //
  // check if event is selected by physics selection class
  //
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
     Printf("ERROR: Could not receive input handler");
     return;
  }

  Bool_t isSelectedMB = kFALSE;
  // check MB
  isSelectedMB = inputHandler->IsEventSelected() & AliVEvent::kMB;

  // get physics selection
  physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
  if(!physicsSelection) return;

  // check Warm Up events
  Bool_t isSelectedWarm = kFALSE;
  triggerAnalysis = physicsSelection->GetTriggerAnalysis();
  if(!triggerAnalysis) return;
  isSelectedWarm = triggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kTPCLaserWarmUp);

  //
  // monitor vertex position
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if(vertex->GetNContributors()<1) {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if(vertex->GetNContributors()<1) vertex = 0x0;
  }  

  Int_t trackCounter = 0;
  for (Int_t i=0; i<fESD->GetNumberOfTracks(); ++i) 
  {
    AliESDtrack *track = fESD->GetTrack(i);
    if (!track) 
      continue;
      
    if (track->GetTPCNcls() < 30) continue;
    //if (track->GetTPCchi2()/track->GetTPCNcls() > 0.3) continue;
    if (TMath::Abs(track->Eta()) > 0.005) continue;
    if (track->Pt() < 4) continue;
    if (track->GetKinkIndex(0) > 0) continue;
    
    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kITSrefit)==1) continue; // explicitly ask for tracks without ITS refit
    if ((status&AliESDtrack::kTPCrefit)==0) continue;
    
    if (track->GetTPCsignal() > 10) continue;          // explicitly ask for tracks without dE/dx
    
    //if (TMath::Abs(track->GetZ()) < 50) continue;
    
    trackCounter++;
  }

  /*
  Float_t dca[2], cov[3];
  // 1st track loop to determine multiplicities
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    AliESDtrack *track = fESD->GetTrack(i);
    if (!track) continue;
    if (fESDtrackCuts->AcceptTrack(track)) {
      UInt_t status = track->GetStatus();
      if ((status&AliESDtrack::kITSrefit)==1) continue; // explicitly ask for tracks without ITS refit
      if (track->GetTPCsignal() > 30) continue;          // explicitly ask for tracks without dE/dx
      //if (track->GetTPCNcls() > 80) continue;
      track->GetImpactParameters(dca, cov);
      if (TMath::Abs(dca[0]) > 10) continue;
      trackCounter++;
      //
    }
  }
  */
 
  // run number, # noise tracks, isSelectedMB, isSelectedWarm, has vtx.
  Bool_t hasVtx = vertex;
  if (trackCounter > 98) trackCounter = 98;
  //Double_t vecNoise[7] = {fESD->GetRunNumber(),fESD->GetEventNumberInFile(), trackCounter, isSelectedMB, isSelectedWarm, hasVtx, fESD->GetTimeStamp()};
  Double_t vecNoise[6] = {fESD->GetEventNumberInFile(), trackCounter, isSelectedMB, isSelectedWarm, hasVtx, fESD->GetTimeStamp()};
  fHistNoiseTracks->Fill(vecNoise);


  if (trackCounter > 0) {
      cout << "NOISE_EVENT_CATEGORY:"<<"\t"<<
	AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()->GetName()<<"\t"<< 
	fESD->GetEventNumberInFile() << "\t" << 
	fESD->GetTimeStamp() <<endl;
  }

  /*
  if (trackCounter > -1) {
    fHistNoiseTracks->Fill(vecNoise);
    if (trackCounter > 25 && trackCounter < 90 && fESD->GetEventSpecie() != AliRecoParam::kCalib) { // dump all availbale info for these events
      cout << "NOISE_EVENT_CATEGORY:"<<"\t"<<
	AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()->GetName()<<"\t"<< 
	fESD->GetEventNumberInFile() << "\t" << 
	fESD->GetTimeStamp() <<endl;
      for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
	AliESDtrack *track = fESD->GetTrack(i);
	if (fESDtrackCuts->AcceptTrack(track)) cout << "NOISE_TRACK:"<<"\t"<<
	  track->GetTPCNcls() <<"\t"<<
	  track->GetTPCsignal() <<"\t"<<
	  track->GetAlpha() <<"\t"<<
	  track->Pt() <<"\t"<<
	  track->GetZ() <<"\t"<<
	  track->Eta() << endl;
      }
    }
    if (fESD->GetEventSpecie() == AliRecoParam::kCalib) {
      cout << "LASER_EVENT_CATEGORY:"<<"\t"<<
	AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree()->GetCurrentFile()->GetName()<<"\t"<< 
	fESD->GetEventNumberInFile() << "\t" << 
	fESD->GetTimeStamp() <<endl;
    }
     
  }
  */
 
  // Post output data  
  PostData(1, fListHist);
  
}      

//________________________________________________________________________
void AliAnalysisNoiseTPC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  Printf("*** CONSTRUCTOR CALLED ****");

}


