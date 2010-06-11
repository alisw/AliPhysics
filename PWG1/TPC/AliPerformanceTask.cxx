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
#include "AliESDfriend.h"
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
#include "AliComparisonObject.h"
#include "AliPerformanceObject.h"
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
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriend(kFALSE)
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
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriend(kFALSE)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AliPerformanceTask::~AliPerformanceTask()
{
  if(fOutput)     delete fOutput;    fOutput   = 0; 
  if(fCompList)   delete fCompList;  fCompList = 0; 
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
void AliPerformanceTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
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
  
  // check output data
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    Printf("ERROR: AliPerformanceTask::Terminate(): fOutput data not avaiable  ..." );
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
