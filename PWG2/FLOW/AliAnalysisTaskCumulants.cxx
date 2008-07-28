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
#include "TList.h"
#include "TH1.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "AliCumulantsFunctions.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

//class AliAnalysisTask;
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "../../CORRFW/AliCFManager.h"

#include "AliAnalysisTaskCumulants.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowAnalysisWithCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHistResults.h"

#include "AliCumulantsFunctions.h"

// AliAnalysisTaskCumulants:
// analysis task for 
// Cumulant method
// with many authors (N.K. R.S. A.B.)
// who do something

ClassImp(AliAnalysisTaskCumulants)

//________________________________________________________________________
AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name) : 
  AliAnalysisTask(name, ""), 
  fESD(NULL),
  fAOD(NULL),
  fMyCumuAnalysis(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"), 
  fCFManager1(NULL),
  fCFManager2(NULL),
  fListHistos(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
}
//________________________________________________________________________
AliAnalysisTaskCumulants::AliAnalysisTaskCumulants() : 
  fESD(NULL),
  fAOD(NULL), 
  fMyCumuAnalysis(NULL),
  fEventMaker(NULL),
  fAnalysisType("ESD"),
  fCFManager1(NULL),
  fCFManager2(NULL),
  fListHistos(NULL)
{
  // Constructor
  cout<<"AliAnalysisTaskCumulants::AliAnalysisTaskCumulants(const char *name)"<<endl;
  }

//________________________________________________________________________
void AliAnalysisTaskCumulants::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once
  cout<<"AliAnalysisTaskCumulants::ConnectInputData(Option_t *)"<<endl;

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
void AliAnalysisTaskCumulants::CreateOutputObjects() 
{
  // Called at every worker node to initialize
  cout<<"AliAnalysisTaskCumulants::CreateOutputObjects()"<<endl;

  if (!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0"  || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) {
    cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
    exit(1);
  }
 
  //event maker
  fEventMaker = new AliFlowEventSimpleMaker();
  //analyser
  fMyCumuAnalysis = new AliFlowAnalysisWithCumulants() ;
   
  fMyCumuAnalysis->CreateOutputObjects();

  if (fMyCumuAnalysis->GetHistList()) {
    //	fSP->GetHistList()->Print();
    fListHistos = fMyCumuAnalysis->GetHistList();
    //	fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }
 }

//________________________________________________________________________
void AliAnalysisTaskCumulants::Exec(Option_t *) 
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
    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    //cumulant analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    fMyCumuAnalysis->Make(fEvent);
    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    //Cumulant analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    fMyCumuAnalysis->Make(fEvent);
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

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    //Cumulant analysis 
    AliFlowEventSimple* fEvent=NULL;
    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
    fMyCumuAnalysis->Make(fEvent);
    delete fEvent;
    //delete mcEvent;
  }
  
  else if (fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());

    // analysis 
    //For the moment don't use CF //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD,fCFManager1,fCFManager2);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD);
    fMyCumuAnalysis->Make(fEvent);
    delete fEvent;
  }


  PostData(0,fListHistos); 
}

//________________________________________________________________________
void AliAnalysisTaskCumulants::Terminate(Option_t *) 
{  
  //=====================================================================================================
  // Accessing the output file which contains the merged results from all workers;
        fListHistos = (TList*)GetOutputData(0);
        //fListHistos->Print();//printing the list of stored histos
  //=====================================================================================================
  
     if(fListHistos)
     {	
      // Profiles with values of generating functions
      TProfile2D *fIntFlowGenFun=dynamic_cast<TProfile2D*>(fListHistos->FindObject("fIntFlowGenFun")); 
      TProfile3D *fDiffFlowGenFunRe=dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowGenFunRe"));
      TProfile3D *fDiffFlowGenFunIm=dynamic_cast<TProfile3D*>(fListHistos->FindObject("fDiffFlowGenFunIm"));
      
      // Histograms to store final results
      TH1D *fIntFlowResults=dynamic_cast<TH1D*>(fListHistos->FindObject("fIntFlowResults"));
      TH1D *fDiffFlowResults2=dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults2"));
      TH1D *fDiffFlowResults4=dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults4"));
      TH1D *fDiffFlowResults6=dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults6"));
      TH1D *fDiffFlowResults8=dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults8"));
     
      // Avarage multiplicity 
      TProfile *AvMult=dynamic_cast<TProfile*>(fListHistos->FindObject("fHistProAvM"));
      Double_t BvM=AvMult->GetBinContent(1);//avarage multiplicity
   
      // Calling the function which from values of generating functions calculate integrated and differential flow and store the final results into histograms 
      AliCumulantsFunctions FinalResults(fIntFlowGenFun,fDiffFlowGenFunRe,fDiffFlowGenFunIm,fIntFlowResults,fDiffFlowResults2,fDiffFlowResults4,fDiffFlowResults6,fDiffFlowResults8,BvM);    
      FinalResults.Calculate();
    }
    else
    {
     cout<<"histogram list pointer is empty"<<endl;
    }
}





















