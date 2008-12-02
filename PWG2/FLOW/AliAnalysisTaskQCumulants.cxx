/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *f
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

/**************************************
 * analysis task for Q-cumulants      * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/
 
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "../../CORRFW/AliCFManager.h"

#include "AliAnalysisTaskQCumulants.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFlowAnalysisWithQCumulants.h"
#include "AliFlowCumuConstants.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHistResults.h"
#include "AliQCumulantsFunctions.h"

ClassImp(AliAnalysisTaskQCumulants)

//================================================================================================================

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name, Bool_t on): 
 AliAnalysisTask(name,""), 
 fESD(NULL),
 fAOD(NULL),
 fQCA(NULL),//Q-cumulant Analysis (QCA) object
 fEventMaker(NULL),
 fAnalysisType("ESD"), 
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(on)
{
 //constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(const char *name)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with a TChain
 DefineInput(0, TChain::Class());
  
 // Output slot #0 writes into a TList container
 DefineOutput(0, TList::Class());  
 if(on) 
 {
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class()); 
 }  
}

AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants(): 
 fESD(NULL),
 fAOD(NULL), 
 fQCA(NULL),//Q-cumulant Analysis (QCA) object
 fEventMaker(NULL),
 fAnalysisType("ESD"),
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),  
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(kFALSE)
{
 //dummy constructor
 cout<<"AliAnalysisTaskQCumulants::AliAnalysisTaskQCumulants()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskQCumulants::ConnectInputData(Option_t *) 
{
 //connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskQCumulants::ConnectInputData(Option_t *)"<<endl;

 TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
 if (!tree) 
 {
  Printf("ERROR: Could not read chain from input slot 0");
 } 
 else 
 {
 //disable all branches and enable only the needed ones
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

//================================================================================================================

void AliAnalysisTaskQCumulants::CreateOutputObjects() 
{
 //called at every worker node to initialize
 cout<<"AliAnalysisTaskQCumulants::CreateOutputObjects()"<<endl;

 
 //OpenFile(0);
 

 if(!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) 
 {
  cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
  exit(1);
 }
 
 //event maker
 fEventMaker = new AliFlowEventSimpleMaker();
 
 //analyser
 fQCA = new AliFlowAnalysisWithQCumulants();
 fQCA->CreateOutputObjects();

 if(fQCA->GetHistList()) 
 {
  fListHistos = fQCA->GetHistList();
  //fListHistos->Print();
 }
 else 
 {
  Printf("ERROR: Could not retrieve histogram list"); 
 }
 
 //PostData(0,fListHistos);
 
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Exec(Option_t *) 
{
 //main loop (called for each event)
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

    //Q-cumulant analysis 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    fQCA->Make(fEvent);
    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    //Q-cumulant analysis    
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    
    fQCA->Make(fEvent); 
    
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

    //Q-cumulant analysis 
    AliFlowEventSimple* fEvent=NULL;
    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
    fQCA->Make(fEvent);
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
    fQCA->Make(fEvent);
    delete fEvent;
  }
  
  PostData(0,fListHistos); 
  if(fQA) 
  {
   PostData(1,fQAInt);
   PostData(2,fQADiff); 
  }
}

//================================================================================================================

void AliAnalysisTaskQCumulants::Terminate(Option_t *) 
{
 //accessing the output list which contains the merged 2D and 3D profiles from all worker nodes
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();
 
 if(fListHistos)
 {	    
  //final results (integrated flow)
  TH1D *intFlowResults = dynamic_cast<TH1D*>(fListHistos->FindObject("fIntFlowResultsQC"));
  
  //final results (differential flow)
  TH1D *diffFlowResults2ndOrder = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults2ndOrderQC"));
  TH1D *diffFlowResults4thOrder = dynamic_cast<TH1D*>(fListHistos->FindObject("fDiffFlowResults4thOrderQC"));
  
  //final results for covariances (1st bin <2*4>-<2>*<4>, 2nd bin <2*6>-<2>*<6>, ...)
  TH1D *covariances = dynamic_cast<TH1D*>(fListHistos->FindObject("fCovariances"));
  
  //common histograms to store the final results for the 2nd order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes2nd = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults2ndOrderQC"));
  
  //common histograms to store the final results for the 4th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes4th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults4thOrderQC"));
  
  //common histograms to store the final results for the 6th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes6th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults6thOrderQC"));
  
  //common histograms to store the final results for the 8th order integrated and differential flow
  AliFlowCommonHistResults *commonHistRes8th = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResults8thOrderQC"));
  
  //average selected multiplicity (for int. flow) 
  TProfile *AvMult = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlowQC"));
  
  //multi-particle correlations calculated from Q-vectors
  TProfile *QCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fQCorrelations"));
  
  //average of products: 1st bin: <2*4>, 2nd bin: <2*6>, ...
  TProfile *QProduct = dynamic_cast<TProfile*>(fListHistos->FindObject("fQProduct"));
  
  //average 2-, 3- and 4-particle correlations per bin 
  TProfile *binned2p_1n1n = dynamic_cast<TProfile*>(fListHistos->FindObject("f2_1n1n"));
  TProfile *binned2p_2n2n = dynamic_cast<TProfile*>(fListHistos->FindObject("f2_2n2n"));
  TProfile *binned3p_2n1n1n = dynamic_cast<TProfile*>(fListHistos->FindObject("f3_2n1n1n"));
  TProfile *binned3p_1n1n2n = dynamic_cast<TProfile*>(fListHistos->FindObject("f3_1n1n2n"));
  TProfile *binned4p_1n1n1n1n = dynamic_cast<TProfile*>(fListHistos->FindObject("f4_1n1n1n1n"));
  
  //average values of Q-vector components (1st bin: <Q_x>, 2nd bin: <Q_y>, 3rd bin: <(Q_x)^2>, 4th bin: <(Q_y)^2>) 
  TProfile *QVectorComponents = dynamic_cast<TProfile*>(fListHistos->FindObject("fQvectorComponents"));
  
  //multi-particle correlations calculated with nested loop 
  TProfile *DirectCorrelations = dynamic_cast<TProfile*>(fListHistos->FindObject("fDirectCorrelations"));
 
 
  fQCA = new AliFlowAnalysisWithQCumulants();  
 
  fQCA->SetIntFlowResults(intFlowResults); 
  fQCA->SetDiffFlowResults2nd(diffFlowResults2ndOrder);
  fQCA->SetDiffFlowResults4th(diffFlowResults4thOrder); 
  fQCA->SetCovariances(covariances); 

  fQCA->SetCommonHistsResults2nd(commonHistRes2nd); 
  fQCA->SetCommonHistsResults4th(commonHistRes4th);
  fQCA->SetCommonHistsResults6th(commonHistRes6th);
  fQCA->SetCommonHistsResults8th(commonHistRes8th);
 
  fQCA->SetAverageMultiplicity(AvMult);
  fQCA->SetQCorrelations(QCorrelations);
  fQCA->SetQProduct(QProduct);
  fQCA->SetQVectorComponents(QVectorComponents);
 
  fQCA->SetTwo_1n1nPerBin(binned2p_1n1n);
  fQCA->SetTwo_2n2nPerBin(binned2p_2n2n);
  fQCA->SetThree_2n1n1nPerBin(binned3p_2n1n1n);
  fQCA->SetThree_1n1n2nPerBin(binned3p_1n1n2n);
  fQCA->SetFour_1n1n1n1nPerBin(binned4p_1n1n1n1n);
 
  fQCA->SetDirectCorrelations(DirectCorrelations);
 
  fQCA->Finish();
 }
 else
 {
  cout<<"histogram list pointer is empty"<<endl;
 }
}





















