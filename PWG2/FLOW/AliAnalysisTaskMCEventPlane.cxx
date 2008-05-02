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
#include "TFile.h" //needed as include
#include "TList.h"


class AliAnalysisTask;
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"

#include "AliAnalysisTaskMCEventPlane.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowAnalysisWithMCEventPlane.h"

// AliAnalysisTaskMCEventPlane:
//
// analysis task for Monte Carlo Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskMCEventPlane)

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name) : 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fAOD(0),
  fMc(0),
  fEventMaker(0),
  fAnalysisType("MC")
{
  // Constructor
  cout<<"AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskMCEventPlane::ConnectInputData(Option_t *)"<<endl;

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    if (fAnalysisType == "MC") {
      cout<<"!!!!!reading MC kinematics only"<<endl;
      // we want to process only MC
      tree->SetBranchStatus("*", kFALSE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else {
	fESD = esdH->GetEvent();
      }
    }
    else if (fAnalysisType == "ESD") {
      cout<<"!!!!!reading the ESD only"<<endl;
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("Tracks.*", kTRUE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else
	fESD = esdH->GetEvent();
    }
    else if (fAnalysisType == "AOD") {
      cout<<"!!!!!reading the AOD only"<<endl;
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      }
      else {
	fAOD = aodH->GetEvent();
      }
    }
    else {
      Printf("!!!!!Wrong analysis type: Only ESD, AOD and MC types are allowed!");
      exit(1);
      
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskMCEventPlane::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, AOD and MC are allowed."<<endl;
    exit(1);
  }

  //event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  //Analyser
  fMc  = new AliFlowAnalysisWithMCEventPlane() ;

  //output file
  TString fName = "outputFromMCEventPlaneAnalysis" ;
  fName += fAnalysisType.Data() ;
  fName += ".root" ;
  fMc->SetHistFileName( fName.Data() );
    
  fMc-> Init();


}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  
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
  
  AliGenCocktailEventHeader *header = dynamic_cast<AliGenCocktailEventHeader *> (mcEvent-> GenEventHeader()); 
  if (!header) {
    Printf("ERROR: Could not retrieve AliGenCocktailEventHeader");
    return;
  }
  
  TList *lhd = header->GetHeaders();
  if (!lhd) {
    Printf("ERROR: Could not retrieve List of headers");
    return;
  }

  AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(0)); 
  if (!hdh) {
    Printf("ERROR: Could not retrieve AliGenHijingEventHeader");
    return;
  }
    
  Double_t fRP = hdh->ReactionPlaneAngle();
  cout<<"The reactionPlane is "<<hdh->ReactionPlaneAngle()<<endl;
  
  if (fAnalysisType == "MC") {
    // analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent);
    fMc->Make(fEvent,fRP);

    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    // analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    fMc->Make(fEvent,fRP);
    delete fEvent;
  }
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    // analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);
    fMc->Make(fEvent,fRP);
    delete fEvent;
  }

}      

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  fMc->Finish();
  PostData(0,fMc->GetHistFile());

  delete fMc;
  delete fEventMaker;
}
