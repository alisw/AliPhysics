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

#include "AliAnalysisTaskLeeYangZeros.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"

// AliAnalysisTaskLeeYangZeros:
// analysis task for Lee Yang Zeros method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLeeYangZeros)

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name, Bool_t firstrun) : 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fAOD(0),
  fAnalysisType("ESD"), 
  fLyz(0),
  fEventMaker(0),
  fFirstRunFile(0),
  fFirstRunLYZ(firstrun), //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE)       //set boolean for use sum to initial value
{
  // Constructor
  cout<<"AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  if (!firstrun) DefineInput(1, TList::Class()); //for second loop 
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskLeeYangZeros::ConnectInputData(Option_t *)"<<endl;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
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
    else if (fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1"  ) {
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

    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskLeeYangZeros::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }

  //event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  //Analyser
  fLyz = new AliFlowAnalysisWithLeeYangZeros() ;
   
  fLyz -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fLyz -> SetUseSum(GetUseSumLYZ());       //set use sum true or false

  //output file
  TString fName = "outputFromLeeYangZerosAnalysis" ;
  fName += fAnalysisType.Data() ;
  if (fFirstRunLYZ) {
    fName += "_firstrun.root" ;
  } else {
    fName += "_secondrun.root" ;
  }
  fLyz->SetHistFileName( fName.Data() );
  
  // Get data from input slot 1
  if (GetNinputs() == 2) {                   //if there are two input slots
    fFirstRunFile = (TFile*)GetInputData(1);
    cerr<<"fFirstRunFile ("<<fFirstRunFile<<")"<<endl;
    if (fFirstRunFile) cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;
    
    fLyz -> SetFirstRunFile(fFirstRunFile);
  }
  
  fLyz-> Init();

}

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Exec(Option_t *) 
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
    fLyz->Make(fEvent);
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
    fLyz->Make(fEvent);
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
    fLyz->Make(fEvent);
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
    fLyz->Make(fEvent);
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
    fLyz->Make(fEvent);
    delete fEvent;
  }
  
}      

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Terminate(Option_t *) 
{
  // Called once at the end of the query
  if (GetNinputs() == 2) { 
    cerr<<"fFirstRunFile -> IsOpen() = "<<fFirstRunFile -> IsOpen()<<endl;
  }
  cerr<<"fLyz->GetHistFile() -> IsOpen() = "<<fLyz->GetHistFile() -> IsOpen()<<endl;

  fLyz->Finish();

  PostData(0,fLyz->GetHistFile());

  delete fLyz;
  delete fEventMaker;

  cout<<".....finished"<<endl;
}
