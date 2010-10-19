/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//------------------------------------------------------------------------------
// Implementation of the AliPerformanceTask class. It checks reconstruction performance 
// for the reconstructed vs MC particle tracks under several conditions. For real data 
// the control QA histograms are filled.
//
// The comparison output objects deriving from AliPerformanceObject 
// (e.g. AliPerformanceRes, AliPerformanceEff, AliPerformanceDEdx, AliPerformanceDCA ...) 
// are stored in the output file (details in description of these classes).
// 
// Author: J.Otwinowski 01/04/2009 
// Changes by M.Knichel 15/10/2010
//------------------------------------------------------------------------------

#include "iostream"

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDVertex.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"

#include "AliMCInfo.h"
#include "AliESDRecInfo.h"
#include "AliMCInfoCuts.h"
#include "AliRecInfoCuts.h"
#include "AliComparisonObject.h"
#include "AliPerformanceObject.h"
#include "AliTPCPerformanceSummary.h"
#include "AliPerformanceTPC.h"
#include "AliPerformanceDEdx.h"
#include "AliPerformanceTask.h"


using namespace std;

ClassImp(AliPerformanceTask)

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask() 
  : AliAnalysisTaskSE("Performance")
  , fESD(0)
  , fESDfriend(0)
  , fMC(0)
  , fOutput(0)
  , fOutputSummary(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriend(kFALSE)
  , fUseHLT(kFALSE)
  , fUseTerminate(kTRUE)
{
  // Dummy Constructor
  // should not be used
}

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask(const char *name, const char */*title*/) 
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fESDfriend(0)
  , fMC(0)
  , fOutput(0)
  , fOutputSummary(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriend(kFALSE)
  , fUseHLT(kFALSE)
  , fUseTerminate(kTRUE)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AliPerformanceTask::~AliPerformanceTask()
{
  if (fOutput)     delete fOutput;    fOutput   = 0; 
  if (fOutputSummary) delete fOutputSummary; fOutputSummary = 0;
  if (fCompList)   delete fCompList;  fCompList = 0; 
}

//_____________________________________________________________________________
Bool_t AliPerformanceTask::AddPerformanceObject(AliPerformanceObject *pObj) 
{
  // add comparison object to the list
  if(pObj == 0) {
    Printf("ERROR: Could not add comparison object");
    return kFALSE;
  }

  // add object to the list
  fCompList->AddLast(pObj);
       
return kTRUE;
}

//_____________________________________________________________________________
void AliPerformanceTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // create output list
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();
  
  // create output list
  //fOutputSummary = new TTree;
  
  // add comparison objects to the output
  AliPerformanceObject *pObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AliPerformanceObject *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
    count++;
  }
  Printf("UserCreateOutputObjects(): Number of output comparison objects: %d \n", count);
  
  PostData(1, fOutput);  
  //PostData(2, fOutputSummary);  
}

//_____________________________________________________________________________
void AliPerformanceTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Decide whether to use HLT or Offline ESD
  if(fUseHLT){

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> 
      (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      printf("ERROR: Could not get ESDInputHandler");
      return;
    } 
    fESD = esdH->GetHLTEvent();
  }// end if fUseHLT
  else  
    fESD = (AliESDEvent*) (InputEvent());

  if(fUseESDfriend)
    {
      fESDfriend = static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
      if(!fESDfriend) {
        Printf("ERROR: ESD friends not available");
      }
    }
  
  if(fUseMCInfo) {
      fMC = MCEvent();
  }  

  //
  AliPerformanceObject *pObj=0;

  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }
  
  if (fUseMCInfo && !fMC) {
    Printf("ERROR: MC event not available");
    return;
  }

  if(fUseESDfriend)
  {
    if(!fESDfriend) {
    Printf("ERROR: ESD friends not available");
    }
  }

  // Process comparison
  fPitList->Reset();
  while(( pObj = (AliPerformanceObject *)fPitList->Next()) != NULL) {
    pObj->Exec(fMC,fESD,fESDfriend,fUseMCInfo,fUseESDfriend);
  }

  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AliPerformanceTask::Terminate(Option_t *) 
{
  // Called once at the end 

  if ( !fUseTerminate )
    return;

  // check output data
    fOutputSummary = dynamic_cast<TTree*> (GetOutputData(2));
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutput) {
        Printf("ERROR: AliPerformanceTask::Terminate(): fOutput data not available  ..." );
        return;
   }
    if (fOutputSummary) { delete fOutputSummary; fOutputSummary=0; }      
    AliPerformanceObject* pObj=0;
    AliPerformanceTPC*  pTPC = 0;
    AliPerformanceDEdx* pDEdx = 0;
    TIterator* itOut = fOutput->MakeIterator();
    itOut->Reset();
    while(( pObj = dynamic_cast<AliPerformanceObject*>(itOut->Next())) != NULL) { 
        if (!  pTPC) {  pTPC = dynamic_cast<AliPerformanceTPC*>(pObj); }
        if (! pDEdx) { pDEdx = dynamic_cast<AliPerformanceDEdx*>(pObj); }
    }
    if (! AliCDBManager::Instance()->GetDefaultStorage()) { AliCDBManager::Instance()->SetDefaultStorage("raw://"); }
    TUUID uuid;
    TString tmpFile = gSystem->TempDirectory() + TString("/TPCQASummary.") + uuid.AsString() + TString(".root");
    AliTPCPerformanceSummary::WriteToFile(pTPC, pDEdx, tmpFile.Data());
    TChain* chain = new TChain("tpcQA");
    chain->Add(tmpFile.Data());
    TTree *tree = chain->CopyTree("1");
    if (chain) { delete chain; chain=0; }
    fOutputSummary = tree;
      
     // Post output data.
     PostData(2, fOutputSummary);

}

//_____________________________________________________________________________
void AliPerformanceTask::FinishTaskOutput()
{
    // called once at the end of each job (on the workernode)
    //
    // projects THnSparse to TH1,2,3
    
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutput) {
        Printf("ERROR: AliPerformanceTask::FinishTaskOutput(): fOutput data not available  ..." );
        return;
   }

      AliPerformanceObject* pObj=0;
      TIterator* itOut = fOutput->MakeIterator();  
      itOut->Reset();
      while(( pObj = dynamic_cast<AliPerformanceObject*>(itOut->Next())) != NULL) {
          pObj->SetRunNumber(fCurrentRunNumber);
          pObj->Analyse();
      }
      
     // Post output data.
     PostData(1, fOutput);
}

//_____________________________________________________________________________
Bool_t AliPerformanceTask::Notify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file", count);

  return kTRUE;
}
