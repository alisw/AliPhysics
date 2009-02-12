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
#include "TProfile.h"
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

#include "AliCFManager.h"

#include "AliGenCocktailEventHeader.h"
#include "AliGenHijingEventHeader.h"

#include "AliAnalysisTaskMCEventPlane.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisWithMCEventPlane.h"

// AliAnalysisTaskMCEventPlane:
//
// analysis task for Monte Carlo Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskMCEventPlane)

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name, Bool_t on) : 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fAOD(0),
  fAnalysisType("MC"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fMc(0),
  fEventMaker(0),
  fListHistos(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(on)
{
  // Constructor
  cout<<"AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class()); 
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); }  
}

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane() : 
  fESD(0),
  fAOD(0),
  fAnalysisType("MC"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fMc(0),
  fEventMaker(0),
  fListHistos(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskMCEventPlane::AliAnalysisTaskMCEventPlane()"<<endl;

}

//________________________________________________________________________
AliAnalysisTaskMCEventPlane::~AliAnalysisTaskMCEventPlane()
{

  //destructor

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
      Printf("!!!!!Wrong analysis type: Only ESD, ESDMC0, ESDMC1, AOD and MC types are allowed!");
      exit(1);
      
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::CreateOutputObjects() 
{
  // Called once
  cout<<"AliAnalysisTaskMCEventPlane::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD"  || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }

  //event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  //Analyser
  fMc  = new AliFlowAnalysisWithMCEventPlane() ;
      
  fMc-> Init();

  if (fMc->GetHistList()) {
    //fMc->GetHistList()->Print();
    fListHistos = fMc->GetHistList();
    //fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }

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

  fCFManager1->SetEventInfo(mcEvent);
  fCFManager2->SetEventInfo(mcEvent);

  Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
  if (mcEvent->GetNumberOfTracks() == -1)
    {
      cout<<"Skipping Event -- No MC information available for this event"<<endl;
      return;
    }
  
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
  //cout<<"The reactionPlane is "<<hdh->ReactionPlaneAngle()<<endl;
  
  if (fAnalysisType == "MC") {
    // analysis 
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
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
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
    fMc->Make(fEvent,fRP);
    delete fEvent;
  }

  else if (fAnalysisType == "ESDMC0") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());

    // analysis 
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,mcEvent,0); //0 = kine from ESD, 1 = kine from MC
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0);
    fMc->Make(fEvent,fRP);
    delete fEvent;
  }

  else if (fAnalysisType == "ESDMC1") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());

    // analysis 
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,mcEvent,1); //0 = kine from ESD, 1 = kine from MC
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1);
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
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);  //no CF yet!
    fMc->Make(fEvent,fRP);
    delete fEvent;
  }

  PostData(0,fListHistos); //here for CAF
  if (fQA) {
    PostData(1,fQAInt);
    PostData(2,fQADiff); }

}      


//________________________________________________________________________
void AliAnalysisTaskMCEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  AliFlowAnalysisWithMCEventPlane* fMcTerm = new AliFlowAnalysisWithMCEventPlane() ;

  //Get output data
  fListHistos = (TList*)GetOutputData(0);
  // cout << "histogram list in Terminate" << endl;
  if (fListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHists = dynamic_cast<AliFlowCommonHist*> 
      (fListHistos->FindObject("AliFlowCommonHistMCEP"));
    AliFlowCommonHistResults *pCommonHistResults = 
      dynamic_cast<AliFlowCommonHistResults*> 
      (fListHistos->FindObject("AliFlowCommonHistResultsMCEP"));

    TProfile *pHistProFlow = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPt_MCEP"));

    TProfile *pHistProIntFlowRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("fHistProIntFlowRP")); 
                               
    TProfile *pHistProDiffFlowPtRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("fHistProDiffFlowPtRP")); 
     
    TProfile *pHistProDiffFlowEtaRP = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("fHistProDiffFlowEtaRP"));
      
    TProfile *pHistProDiffFlowPtPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("fHistProDiffFlowPtPOI")); 
     
    TProfile *pHistProDiffFlowEtaPOI = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("fHistProDiffFlowEtaPOI"));                             

    if (pCommonHists && pCommonHistResults && pHistProFlow && pHistProIntFlowRP && pHistProDiffFlowPtRP && pHistProDiffFlowEtaRP && pHistProDiffFlowPtPOI && pHistProDiffFlowEtaPOI) {
      fMcTerm->SetCommonHists(pCommonHists);
      fMcTerm->SetCommonHistsRes(pCommonHistResults);
      fMcTerm->SetHistProFlow(pHistProFlow);
      fMcTerm->SetHistProIntFlowRP(pHistProIntFlowRP);
      fMcTerm->SetHistProDiffFlowPtRP(pHistProDiffFlowPtRP);      
      fMcTerm->SetHistProDiffFlowEtaRP(pHistProDiffFlowEtaRP);  
      fMcTerm->SetHistProDiffFlowPtPOI(pHistProDiffFlowPtPOI);      
      fMcTerm->SetHistProDiffFlowEtaPOI(pHistProDiffFlowEtaPOI);          
      fMcTerm->Finish();
      PostData(0,fListHistos);
    } else {
      cout<<"WARNING: Histograms needed to run Finish() are not accessable!"<<endl;  }
    
    //fListHistos->Print();
  } else { cout << "histogram list pointer is empty" << endl;}
    
  cout<<"...finished MCEventPlane."<<endl;
}
