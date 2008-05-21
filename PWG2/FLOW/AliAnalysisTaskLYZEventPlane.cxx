/*************************************************************************
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

#include "Riostream.h" //needed as include
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

class AliAnalysisTask;
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliAnalysisTaskLYZEventPlane.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowLYZEventPlane.h"
#include "AliFlowAnalysisWithLYZEventPlane.h"

// AliAnalysisTaskLYZEventPlane:
//
// analysis task for Lee Yang Zeros Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLYZEventPlane)

//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name) : 
  AliAnalysisTask(name, ""), 
  fESD(0), 
  fAOD(0),
  fAnalysisType("ESD"), 
  fLyzEp(0),
  fLyz(0),
  fEventMaker(0),
  fFirstRunFile(0),
  fSecondRunFile(0)
{
  // Constructor
  cout<<"AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineInput(1, TList::Class());
  DefineInput(2, TList::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskLYZEventPlane::ConnectInputData(Option_t *)"<<endl;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    if (fAnalysisType == "MC") {
      // we want to process only MC
      tree->SetBranchStatus("*", kFALSE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else {
	fESD = esdH->GetEvent();
      }
    }
    else if (fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" ) {
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("Tracks.*", kTRUE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else 
	fESD = esdH->GetEvent();
    }
    else if (fAnalysisType == "AOD") {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      }
      else {
	fAOD = aodH->GetEvent();
      }
    }
    else {
      Printf("Wrong analysis type: Only ESD, ESDMC0, ESDMC1, AOD and MC types are allowed!");
      exit(1);

    }
  }  
}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskLYZEventPlane::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD"  || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, , ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }

  //event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  //lee yang zeros event plane
  fLyzEp = new AliFlowLYZEventPlane() ;
  //Analyser
  fLyz = new AliFlowAnalysisWithLYZEventPlane() ;

  //output file
  TString fName = "outputFromLYZEventPlaneAnalysis" ;
  fName += fAnalysisType.Data() ;
  fName += ".root" ;
  fLyz->SetOutFileName( fName.Data() );
     
  // Get data from input slot 1 and 2
  fFirstRunFile = (TFile*)GetInputData(1);
  cerr<<"fFirstRunFile ("<<fFirstRunFile<<")"<<endl;
  if (fFirstRunFile) cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;
  fSecondRunFile = (TFile*)GetInputData(2);
  cerr<<"fSecondRunFile ("<<fSecondRunFile<<")"<<endl;
  if (fSecondRunFile) cerr<<"fSecondRunFile -> IsOpen() = "<<fSecondRunFile -> IsOpen()<<endl;
  
  fLyzEp -> SetFirstRunFile(fFirstRunFile);
  fLyzEp -> SetSecondRunFile(fSecondRunFile);
  fLyz -> SetFirstRunFile(fFirstRunFile);
  fLyz -> SetSecondRunFile(fSecondRunFile);

  fLyzEp-> Init();
  fLyz-> Init();


}

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  if (fAnalysisType == "MC") {
    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event

    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent);
    fLyz->Make(fEvent,fLyzEp);

    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
  }
  else if (fAnalysisType == "ESDMC0") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,mcEvent,0); //0 = kine from ESD, 1 = kine from MC
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
    //delete mcEvent;
  }
  else if (fAnalysisType == "ESDMC1") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }

    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,mcEvent,1); //0 = kine from ESD, 1 = kine from MC
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
    //delete mcEvent;
  }
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
  }
  
  PostData(0,fLyz->GetOutFile());
}      

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  fLyz->Finish();

  delete fLyz;
  delete fLyzEp;
  delete fEventMaker;
}
