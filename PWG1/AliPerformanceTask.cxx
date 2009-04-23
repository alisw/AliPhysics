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
// (e.g. AliPerformanceRes, AliPerformanceEff, AliPerformanceDEdxA, AliPerformanceDCA ...) 
// are stored in the output file (details in description of these classes).
// 
// Author: J.Otwinowski 01/04/2009 
//------------------------------------------------------------------------------

#include "iostream"

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDVertex.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliGeomManager.h"

#include "AliMCInfo.h"
#include "AliESDRecInfo.h"
#include "AliMCInfoCuts.h"
#include "AliRecInfoCuts.h"
#include "AliComparisonRes.h"
#include "AliComparisonEff.h"
#include "AliComparisonDEdx.h"
#include "AliComparisonDCA.h"
#include "AliComparisonObject.h"
#include "AliPerformanceObject.h"
#include "AliPerformanceTask.h"

using namespace std;

ClassImp(AliPerformanceTask)

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask() 
  : AliAnalysisTask("Performance","Detector Performance")
  , fESD(0)
  , fMC(0)
  , fOutput(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
{
  // Dummy Constructor
  // should not be used
}

//_____________________________________________________________________________
AliPerformanceTask::AliPerformanceTask(const char *name, const char *title) 
  : AliAnalysisTask(name, title)
  , fESD(0)
  , fMC(0)
  , fOutput(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AliPerformanceTask::~AliPerformanceTask()
{
  if(fOutput)   delete fOutput;  fOutput =0; 
  if(fCompList)   delete fCompList;  fCompList =0; 
}

//_____________________________________________________________________________
void AliPerformanceTask::ConnectInputData(Option_t *) 
{
  // Connect input data 
  // Called once

  TTree *tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
    return;
  }

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
  } else {
    fESD = esdH->GetEvent();

    // Enable only the needed branches
    //esdH->SetActiveBranches("AliESDHeader Vertex Tracks");
  }

  // use MC information
  if(fUseMCInfo) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
    } else {
      fMC = eventHandler->MCEvent();
    }
  }
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
void AliPerformanceTask::CreateOutputObjects()
{
  // Create histograms
  // Called once

  // create output list
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();

  // add comparison objects to the output
  AliPerformanceObject *pObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AliPerformanceObject *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
    count++;
  }
  Printf("CreateOutputObjects(): Number of output comparison objects: %d \n", count);
}

//_____________________________________________________________________________
void AliPerformanceTask::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliPerformanceObject *pObj=0;

  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }
  
  if (fUseMCInfo && !fMC) {
    Printf("ERROR: MC event not available");
    return;
  }

  // Process comparison
  fPitList->Reset();
  while(( pObj = (AliPerformanceObject *)fPitList->Next()) != NULL) {
    pObj->Exec(fMC,fESD,fUseMCInfo);
  }

  // Post output data.
  PostData(0, fOutput);
}

//_____________________________________________________________________________
void AliPerformanceTask::Terminate(Option_t *) 
{
  // Called one at the end 
  
  // check output data
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: AliPerformanceTask::Terminate(): Output data not avaiable GetOutputData(0)==0x0 ..." );
    return;
  }
}

//_____________________________________________________________________________
Bool_t AliPerformanceTask::Notify()
{
  static Int_t count = 0;
  count++;
  Printf("Processing %d. file", count);

  return kTRUE;
}
