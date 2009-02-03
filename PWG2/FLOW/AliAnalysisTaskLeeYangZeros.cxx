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
#include "TProfile.h"

class AliAnalysisTask;
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "../../CORRFW/AliCFManager.h"

#include "AliFlowLYZConstants.h"   
#include "AliAnalysisTaskLeeYangZeros.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowLYZHist1.h"
#include "AliFlowLYZHist2.h"
#include "AliFlowAnalysisWithLeeYangZeros.h"

// AliAnalysisTaskLeeYangZeros:
// analysis task for Lee Yang Zeros method
// Author: Naomi van der Kolk (kolk@nikhef.nl)

ClassImp(AliAnalysisTaskLeeYangZeros)

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name, Bool_t firstrun, Bool_t on) : 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fAOD(0),
  fAnalysisType("ESD"), 
  fCFManager1(NULL),
  fCFManager2(NULL),
  fLyz(0),
  fEventMaker(0),
  fFirstRunFile(0),
  fListHistos(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fFirstRunLYZ(firstrun),  //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE),       //set boolean for use sum to initial value
  fQA(on)
{
  // Constructor
  cout<<"AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros(const char *name)"<<endl;

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  if (!firstrun) DefineInput(1, TList::Class()); //for second loop 
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  
  if(on) {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); }  
} 



//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros() :  
  fESD(0),
  fAOD(0),
  fAnalysisType("ESD"), 
  fCFManager1(NULL),
  fCFManager2(NULL),
  fLyz(0),
  fEventMaker(0),
  fFirstRunFile(0),
  fListHistos(NULL),
  fQAInt(NULL),
  fQADiff(NULL),
  fFirstRunLYZ(kTRUE), //set boolean for firstrun to initial value
  fUseSumLYZ(kTRUE),    //set boolean for use sum to initial value
  fQA(kFALSE)
{
  // Constructor
  cout<<"AliAnalysisTaskLeeYangZeros::AliAnalysisTaskLeeYangZeros()"<<endl;

}

//________________________________________________________________________
AliAnalysisTaskLeeYangZeros::~AliAnalysisTaskLeeYangZeros()
{

  //destructor

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

  // Get data from input slot 1
  if (GetNinputs() == 2) {                   //if there are two input slots
    TList* pFirstRunList = (TList*)GetInputData(1);
    if (pFirstRunList) {
      fLyz -> SetFirstRunList(pFirstRunList);
    } else { cout<<"No first run List!"<<endl; exit(0); }
  }
  
  fLyz -> Init();

  if (fLyz->GetHistList()) {
    fListHistos = fLyz->GetHistList();
    //    fListHistos->Print();
  }
  else {Printf("ERROR: Could not retrieve histogram list"); }
  
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

    fCFManager1->SetEventInfo(mcEvent);
    fCFManager2->SetEventInfo(mcEvent);

    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

    //lee yang zeros analysis 
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
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
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
    fLyz->Make(fEvent);
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
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,mcEvent,0); //0 = kine from ESD, 1 = kine from MC
    AliFlowEventSimple* fEvent=NULL;
    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
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
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fAOD); //no CF yet!
    fLyz->Make(fEvent);
    delete fEvent;
  }
  
  PostData(0,fListHistos); //here for CAF
  if (fQA) {
    PostData(1,fQAInt);
    PostData(2,fQADiff); }
}      

//________________________________________________________________________
void AliAnalysisTaskLeeYangZeros::Terminate(Option_t *) 
{
  // Called once at the end of the query

  const Int_t iNtheta = AliFlowLYZConstants::kTheta;

  AliFlowAnalysisWithLeeYangZeros* fLyzTerm = new AliFlowAnalysisWithLeeYangZeros() ;
  fLyzTerm -> SetFirstRun(GetFirstRunLYZ());   //set first run true or false
  fLyzTerm -> SetUseSum(GetUseSumLYZ());       //set use sum true or false
   
  fListHistos = (TList*)GetOutputData(0);
  //cout << "histogram list in Terminate" << endl;

  if (fListHistos) {

    //define histograms for first and second run
    AliFlowCommonHist *pCommonHist = NULL;
    AliFlowCommonHistResults *pCommonHistResults = NULL;
    TProfile* pHistProVtheta = NULL;
    TProfile* pHistProReDenom = NULL;
    TProfile* pHistProImDenom = NULL;
    TProfile* pHistProReDtheta = NULL;
    TProfile* pHistProImDtheta = NULL;
    TProfile* pHistProVeta = NULL;
    TProfile* pHistProVPt  = NULL;
    AliFlowLYZHist1 *pLYZHist1[iNtheta] = {NULL};      //array of pointers to AliFlowLYZHist1
    AliFlowLYZHist2 *pLYZHist2[iNtheta] = {NULL};      //array of pointers to AliFlowLYZHist2

    if (GetFirstRunLYZ()) { //first run
      //Get the common histograms from the output list
      pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	(fListHistos->FindObject("AliFlowCommonHistLYZ1"));
      pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	(fListHistos->FindObject("AliFlowCommonHistResultsLYZ1"));
    }
    else { //second run
      //Get the common histograms from the output list
      pCommonHist = dynamic_cast<AliFlowCommonHist*> 
	(fListHistos->FindObject("AliFlowCommonHistLYZ2"));
      pCommonHistResults = dynamic_cast<AliFlowCommonHistResults*> 
	(fListHistos->FindObject("AliFlowCommonHistResultsLYZ2"));
    }

    TProfile* pHistProR0theta = dynamic_cast<TProfile*> 
      (fListHistos->FindObject("First_FlowPro_r0theta_LYZ"));

    TH1F* pHistQsumforChi = dynamic_cast<TH1F*> 
      (fListHistos->FindObject("Flow_QsumforChi_LYZ"));

    
    if (GetFirstRunLYZ()) { //for firstrun
      //Get the histograms from the output list
      for(Int_t theta = 0;theta<iNtheta;theta++){
	TString name = "AliFlowLYZHist1_"; 
	name += theta;
	pLYZHist1[theta] = dynamic_cast<AliFlowLYZHist1*> 
	  (fListHistos->FindObject(name));
      }
      pHistProVtheta = dynamic_cast<TProfile*> 
	  (fListHistos->FindObject("First_FlowPro_Vtheta_LYZ"));

      //Set the histogram pointers and call Finish()
      if (pCommonHist && pCommonHistResults && pLYZHist1[0] && 
	  pHistProVtheta && pHistProR0theta && pHistQsumforChi ) {
	fLyzTerm->SetCommonHists(pCommonHist);
	fLyzTerm->SetCommonHistsRes(pCommonHistResults);
	fLyzTerm->SetHist1(pLYZHist1);
	fLyzTerm->SetHistProVtheta(pHistProVtheta);
	fLyzTerm->SetHistProR0theta(pHistProR0theta);
	fLyzTerm->SetHistQsumforChi(pHistQsumforChi);
	fLyzTerm->Finish();
	PostData(0,fListHistos);
      } else { 
	cout<<"WARNING: Histograms needed to run Finish() firstrun are not accessable!"<<endl; 
      }
    } else { //for second run
      //Get the histograms from the output list
      for(Int_t theta = 0;theta<iNtheta;theta++){
	TString name = "AliFlowLYZHist2_"; 
	name += theta;
	pLYZHist2[theta] = dynamic_cast<AliFlowLYZHist2*> 
	  (fListHistos->FindObject(name));
      }
      
      pHistProReDenom = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ReDenom_LYZ"));
      pHistProImDenom = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ImDenom_LYZ"));

      pHistProReDtheta = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ReDtheta_LYZ"));
      pHistProImDtheta = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_ImDtheta_LYZ"));

      pHistProVeta = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_Veta_LYZ"));
      pHistProVPt = dynamic_cast<TProfile*> 
	(fListHistos->FindObject("Second_FlowPro_VPt_LYZ"));

      //Set the histogram pointers and call Finish()
      if (pCommonHist && pCommonHistResults && pLYZHist2[0] && pHistProR0theta && 
	  pHistProReDenom && pHistProImDenom && pHistProVeta && pHistProVPt) {
	fLyzTerm->SetCommonHists(pCommonHist);
	fLyzTerm->SetCommonHistsRes(pCommonHistResults);
	fLyzTerm->SetHist2(pLYZHist2);
	fLyzTerm->SetHistProR0theta(pHistProR0theta);
	fLyzTerm->SetHistProReDenom(pHistProReDenom);
	fLyzTerm->SetHistProImDenom(pHistProImDenom);
	fLyzTerm->SetHistProReDtheta(pHistProReDtheta);
	fLyzTerm->SetHistProImDtheta(pHistProImDtheta);
	fLyzTerm->SetHistProVeta(pHistProVeta);
	fLyzTerm->SetHistProVPt(pHistProVPt);
	fLyzTerm->SetHistQsumforChi(pHistQsumforChi);
	fLyzTerm->Finish();
	PostData(0,fListHistos);
      } else { 
	cout<<"WARNING: Histograms needed to run Finish() secondrun are not accessable!"<<endl; 
      }
    }
          
    //    fListHistos->Print(); 
  }	
  else { cout << "histogram list pointer is empty" << endl;}

  cout<<".....finished LYZ"<<endl;
}
