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
 *    analysis task for fitting       * 
 *         q-distribution             *
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

#include "AliCFManager.h"

#include "AliAnalysisTaskFittingQDistribution.h"
#include "AliFlowEventSimpleMaker.h"
#include "AliFittingQDistribution.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowCommonHistResults.h"
#include "AliFittingFunctionsForQDistribution.h"

ClassImp(AliAnalysisTaskFittingQDistribution)

//================================================================================================================

AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(const char *name, Bool_t on, Bool_t useWeights): 
 AliAnalysisTask(name,""), 
 fESD(NULL),
 fAOD(NULL),
 fFQDA(NULL),//Fitting Q_Distribution Analysis (FQDA) object
 fEventMaker(NULL),
 fAnalysisType("ESD"), 
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(on),
 fUseWeights(useWeights),
 fUsePhiWeights(kFALSE),
 fListWeights(NULL)
{
 //constructor
 cout<<"AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(const char *name)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with a TChain
 DefineInput(0, TChain::Class());
 
 // Input slot #1 is needed for the weights 
 if(useWeights)
 {
  DefineInput(1, TList::Class());   
 }
  
 // Output slot #0 writes into a TList container
 DefineOutput(0, TList::Class());  
 if(on) 
 {
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class()); 
 }         
}

AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution(): 
 fESD(NULL),
 fAOD(NULL), 
 fFQDA(NULL),//Fitting q-distribution Analysis (FQDA) object
 fEventMaker(NULL),
 fAnalysisType("ESD"),
 fCFManager1(NULL),
 fCFManager2(NULL),
 fListHistos(NULL),  
 fQAInt(NULL),
 fQADiff(NULL),
 fQA(kFALSE),
 fUseWeights(kFALSE),
 fUsePhiWeights(kFALSE),
 fListWeights(NULL)
{
 //dummy constructor
 cout<<"AliAnalysisTaskFittingQDistribution::AliAnalysisTaskFittingQDistribution()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::ConnectInputData(Option_t *) 
{
 //connect ESD or AOD (called once)
 cout<<"AliAnalysisTaskFittingQDistribution::ConnectInputData(Option_t *)"<<endl;

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

void AliAnalysisTaskFittingQDistribution::CreateOutputObjects() 
{
 //called at every worker node to initialize
 cout<<"AliAnalysisTaskFittingQDistribution::CreateOutputObjects()"<<endl;

 
 //OpenFile(0);
 

 if(!(fAnalysisType == "AOD" || fAnalysisType == "ESD" || fAnalysisType == "ESDMC0" || fAnalysisType == "ESDMC1" || fAnalysisType == "MC")) 
 {
  cout<<"WRONG ANALYSIS TYPE! only ESD, ESDMC0, ESDMC1, AOD and MC are allowed."<<endl;
  exit(1);
 }
 
 //event maker
 fEventMaker = new AliFlowEventSimpleMaker();
  
 //analyser
 fFQDA = new AliFittingQDistribution();
 fFQDA->Init();
 
 //weights:
 if(fUseWeights)
 {
  //pass the flags to class:
  if(fUsePhiWeights) fFQDA->SetUsePhiWeights(fUsePhiWeights);
  //get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fListWeights = (TList*)GetInputData(1); 
  }
  //pass the list with weights to class:
  if(fListWeights) fFQDA->SetWeightsList(fListWeights);
 }
 
 if(fFQDA->GetHistList()) 
 {
  fListHistos = fFQDA->GetHistList();
  //fListHistos->Print();
 }
 else 
 {
  Printf("ERROR: Could not retrieve histogram list"); 
 }
 
 //PostData(0,fListHistos);
 
}

//================================================================================================================

void AliAnalysisTaskFittingQDistribution::Exec(Option_t *) 
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

    //fitting q-distribution 
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(mcEvent,fCFManager1,fCFManager2);
    fFQDA->Make(fEvent);
    delete fEvent;
  }
  else if (fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    
    //fitting q-distribution
    AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD,fCFManager1,fCFManager2);
    //AliFlowEventSimple* fEvent = fEventMaker->FillTracks(fESD);
    fFQDA->Make(fEvent);
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

    //fitting q-distribution 
    AliFlowEventSimple* fEvent=NULL;
    if (fAnalysisType == "ESDMC0") { 
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 0); //0 = kine from ESD, 1 = kine from MC
    } else if (fAnalysisType == "ESDMC1") {
      fEvent = fEventMaker->FillTracks(fESD, mcEvent, fCFManager1, fCFManager2, 1); //0 = kine from ESD, 1 = kine from MC
    }
    fFQDA->Make(fEvent);
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
    fFQDA->Make(fEvent);
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

void AliAnalysisTaskFittingQDistribution::Terminate(Option_t *) 
{  
 //accessing the output list
 fListHistos = (TList*)GetOutputData(0);
 //fListHistos->Print();
 
 if(fListHistos)
 {	    
  //final results (integrated flow)
  TH1D *intFlowResults = dynamic_cast<TH1D*>(fListHistos->FindObject("fIntFlowResultsFQD")); 
  
  //sigma^2
  TH1D *sigma2 = dynamic_cast<TH1D*>(fListHistos->FindObject("fSigma2")); 
  
  //common histograms to store the final results for the integrated flow
  AliFlowCommonHistResults *commonHistRes = dynamic_cast<AliFlowCommonHistResults*>(fListHistos->FindObject("AliFlowCommonHistResultsFQD"));
 
  //average selected multiplicity (for int. flow) 
  TProfile *AvMult = dynamic_cast<TProfile*>(fListHistos->FindObject("fAvMultIntFlowFQD"));
  
  //q-distribution
  TH1D *qDist = dynamic_cast<TH1D*>(fListHistos->FindObject("fQDistributionFQD"));
  
  //----------------------------------------------------
  
  fFQDA = new AliFittingQDistribution();
  
  fFQDA->SetIntFlowResults(intFlowResults);
  fFQDA->SetSigma2(sigma2); 
  fFQDA->SetCommonHistsResults(commonHistRes); 
  
  fFQDA->SetAverageMultiplicity(AvMult);
  fFQDA->SetQDistribution(qDist); 
  
  fFQDA->Finish();  
  
  //----------------------------------------------------      
 }
 else
 {
  cout<<"histogram list pointer is empty"<<endl;
 }
}

//================================================================================================================



















