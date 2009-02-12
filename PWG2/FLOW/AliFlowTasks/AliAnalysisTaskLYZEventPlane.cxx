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
#include "TProfile.h"
#include "TFile.h"
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

#include "AliAnalysisTaskLYZEventPlane.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZEventPlane.h"
#include "AliFlowAnalysisWithLYZEventPlane.h"

// AliAnalysisTaskLYZEventPlane:
//
// analysis task for Lee Yang Zeros Event Plane
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLYZEventPlane)

//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name, Bool_t on) : 
  AliAnalysisTask(name, ""), 
  fESD(NULL), 
  fAOD(NULL),
  fAnalysisType("ESD"), 
  fLyzEp(NULL),
  fLyz(NULL),
  fEventMaker(NULL),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fListHistos(NULL),
  fSecondRunFile(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(on)
{
  // Constructor
  cout<<"AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineInput(1, TList::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); } 
}

//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane() : 
  fESD(NULL), 
  fAOD(NULL),
  fAnalysisType("ESD"), 
  fLyzEp(NULL),
  fLyz(NULL),
  fEventMaker(NULL),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fListHistos(NULL),
  fSecondRunFile(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fQA(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskLYZEventPlane::AliAnalysisTaskLYZEventPlane()"<<endl;
}


//________________________________________________________________________
AliAnalysisTaskLYZEventPlane::~AliAnalysisTaskLYZEventPlane() 
{
  //destructor

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
     
  // Get data from input slot
  TList* pSecondRunList = (TList*)GetInputData(1);
  if (pSecondRunList) {
    fLyzEp -> SetSecondRunList(pSecondRunList);
    fLyz -> SetSecondRunList(pSecondRunList);
  } else { cout<<"No Second run List!"<<endl; exit(0); }

  fLyzEp-> Init();
  fLyz-> Init();

  if (fLyz->GetHistList()) {
    fListHistos = fLyz->GetHistList();
    fListHistos->Print();
  }
  else { cout<<"ERROR: Could not retrieve histogram list"<<endl;}


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

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
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
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
  }
  else if (fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1") {
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

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    //lee yang zeros analysis 
    AliFlowEventSimple* fEvent = NULL;
    if (fAnalysisType == "ESDMC0") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent,fCFManager1, fCFManager2,0);  //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent,fCFManager1, fCFManager2,0);  //0 = kine from ESD, 1 = kine from MC
    }

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
    //for the moment don't use CF
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD,fCFManager1,fCFManager2);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);
    fLyz->Make(fEvent,fLyzEp);
    delete fEvent;
  }
  
  PostData(0,fListHistos);
  if (fQA) {
    PostData(1,fQAInt);
    PostData(2,fQADiff); }
}      

//________________________________________________________________________
void AliAnalysisTaskLYZEventPlane::Terminate(Option_t *) 
{
  // Called once at the end of the query
  AliFlowAnalysisWithLYZEventPlane* fLyzTerm = new AliFlowAnalysisWithLYZEventPlane() ;
  fListHistos = (TList*)GetOutputData(0);
  //cout << "histogram list in Terminate" << endl;
  if (fListHistos) {
    //Get the common histograms from the output list
    AliFlowCommonHist *pCommonHist = dynamic_cast<AliFlowCommonHist*> 
      (fListHistos->FindObject("AliFlowCommonHistLYZEP"));
    AliFlowCommonHistResults *pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
      (fListHistos->FindObject("AliFlowCommonHistResultsLYZEP"));

    TProfile* pHistProR0theta = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("First_FlowPro_r0theta_LYZ"));

    TProfile* pHistProFlow = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("FlowPro_VPt_LYZEP"));

    TH1F* pHistQsumforChi = dynamic_cast<TH1F*> 
      (fListHistos->FindObject("Flow_QsumforChi_LYZEP"));

    if (pCommonHist && pCommonHistResults && pHistProR0theta &&
	pHistProFlow && pHistQsumforChi ) {
    fLyzTerm->SetCommonHists(pCommonHist);
    fLyzTerm->SetCommonHistsRes(pCommonHistResults);
    fLyzTerm->SetFirstr0theta(pHistProR0theta);
    fLyzTerm->SetHistProFlow(pHistProFlow);
    fLyzTerm->SetHistQsumforChi(pHistQsumforChi);
    fLyzTerm->Finish();
    PostData(0,fListHistos);
    } else { 
      cout<<"WARNING: Histograms needed to run Finish() are not accessable!"<<endl; 
    }

    fListHistos->Print(); 
  } else { cout << "histogram list pointer is empty" << endl;}

  cout<<".....finished LYZ EventPlane"<<endl;  
}


