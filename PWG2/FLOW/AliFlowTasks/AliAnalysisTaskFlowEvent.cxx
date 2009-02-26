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

// ALICE Analysis Framework
class AliAnalysisTask;
#include "AliAnalysisManager.h"

// ESD interface
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

// AOD interface
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

// Monte Carlo Event
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

// ALICE Correction Framework
#include "AliCFManager.h"

// Interface to Event generators to get Reaction Plane Angle
#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEventSimpleMaker.h"

#include "AliAnalysisTaskFlowEvent.h"

// AliAnalysisTaskFlowEvent:
//
// analysis task for filling the flow event
// from MCEvent, ESD, AOD ....
// and put it in an output stream so it can 
// be used by the various flow analysis methods 
// for cuts the correction framework is used
// which also outputs QA histrograms to view
// the effects of the cuts


ClassImp(AliAnalysisTaskFlowEvent)

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name, Bool_t on) : 
  AliAnalysisTask(name, ""), 
//  fOutputFile(NULL),
  fESD(NULL),
  fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(on)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Define here the flow event output
  DefineOutput(0, AliFlowEventSimple::Class());  
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); }  
  // and for testing open an output file
  //  fOutputFile = new TFile("FlowEvents.root","RECREATE");

}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent() : 
  //  fOutputFile(NULL),
  fESD(NULL),
  fAOD(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskFlowEvent::AliAnalysisTaskFlowEvent()"<<endl;
}

//________________________________________________________________________
AliAnalysisTaskFlowEvent::~AliAnalysisTaskFlowEvent()
{
  //
  // Destructor
  //

  // objects in the output list are deleted 
  // by the TSelector dtor (I hope)

}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskFlowEvent::ConnectInputData(Option_t *)"<<endl;

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
    else if (fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1") {
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
      Printf("!!!!!Wrong analysis type: Only ESD, ESDMC0, ESDMC1, AOD and MC types are allowed!");
      exit(1);
      
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::CreateOutputObjects() 
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskFlowEvent::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }

  // Flow Event maker
  fEventMaker = new AliFlowEventSimpleMaker();
}

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  AliFlowEventSimple* fEvent = NULL;
  Double_t fRP = 0.; // the monte carlo reaction plane angle
  AliMCEvent* mcEvent = NULL;
  // See if we can get Monte Carlo Information and if so get the reaction plane

  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandler) {
    mcEvent = eventHandler->MCEvent();
    if (mcEvent) {
      AliGenCocktailEventHeader *header = dynamic_cast<AliGenCocktailEventHeader *> (mcEvent-> GenEventHeader()); 
      if (header) {
	TList *lhd = header->GetHeaders();
	if (lhd) {
	  AliGenHijingEventHeader *hdh = dynamic_cast<AliGenHijingEventHeader *> (lhd->At(0)); 
	  if (hdh) {
	    fRP = hdh->ReactionPlaneAngle();
	    //	    cout<<"The reactionPlane is: "<< fRP <<endl;
	  }
	}
      }
    }
  }

  // set the value of the monte carlo event plane for the flow event
  fEventMaker->SetMCReactionPlaneAngle(fRP);

  // Fill the FlowEventSimple for MC input          
  if (fAnalysisType == "MC") {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }

    // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
    // This handler can return the current MC event
    if (!mcEvent) { Printf("ERROR: Could not retrieve MC event"); return;}

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    // analysis 
    Printf("Number of MC particles: %d", mcEvent->GetNumberOfTracks());
    fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    // here we have the fEvent and want to make it available as an output stream
    // so no delete fEvent;
  }
  // Fill the FlowEventSimple for ESD input  
  else if (fAnalysisType == "ESD") {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }

    if (!fESD) { Printf("ERROR: fESD not available"); return;}
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    // analysis 
    fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
  }
  // Fill the FlowEventSimple for ESD input combined with MC info  
  else if (fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" ) {
    if (!fCFManager1) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fCFManager2) {cout << "ERROR: No pointer to correction framework cuts! " << endl; return; }
    if (!fESD) { Printf("ERROR: fESD not available"); return;}
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    if (!mcEvent) {Printf("ERROR: Could not retrieve MC event"); return;}

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);


    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
  }
  // Fill the FlowEventSimple for AOD input  
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {Printf("ERROR: fAOD not available"); return;}
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    // analysis 
    //For the moment don't use CF //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD,fCFManager1,fCFManager2);
    fEvent = fEventMaker->FillTracks(fAOD);
  }

  //fListHistos->Print();
  //  fOutputFile->WriteObject(fEvent,"myFlowEventSimple");	
  PostData(0,fEvent);
  if (fQA) {
    PostData(1,fQAInt);
    PostData(2,fQADiff); }
} 

//________________________________________________________________________
void AliAnalysisTaskFlowEvent::Terminate(Option_t *) 
{
  // Called once at the end of the query -- do not call in case of CAF

}
