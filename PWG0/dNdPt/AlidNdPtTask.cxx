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
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliTracker.h"
#include "AliGeomManager.h"

#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "dNdPt/AlidNdPt.h"
#include "dNdPt/AlidNdPtEventCuts.h"
#include "dNdPt/AlidNdPtAcceptanceCuts.h"

#include "dNdPt/AlidNdPtTask.h"

using namespace std;

ClassImp(AlidNdPtTask)

//_____________________________________________________________________________
AlidNdPtTask::AlidNdPtTask(const char *name) 
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fMC(0)
  , fOutput(0)
  , fPitList(0)
  , fCompList(0)
  , fUseMCInfo(kFALSE)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AlidNdPtTask::~AlidNdPtTask()
{
  if(fOutput) delete fOutput;  fOutput =0; 
  if(fCompList) delete fCompList;  fCompList =0; 
}

//____________________________________________________________________________
Bool_t AlidNdPtTask::Notify()
{
  static Int_t count = 0;
  count++;
  //Printf("Processing %d. file: %s", count, ((TTree*) GetInputData(0))->GetCurrentFile()->GetName());
  TTree *chain = (TChain*)GetInputData(0);
  if(chain)
    Printf("Processing %d. file: %s", count, chain->GetCurrentFile()->GetName());

  /*
  TChain *chain = (TChain*)GetInputData(0);
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdH) {
    Printf("ERROR: Could not get ESDInputHandler");
    return kFALSE;
  } else {
    if(chain)
    Printf("chain->GetCurrentFile()->GetName() %s", chain->GetCurrentFile()->GetName());
  }
  */

return kTRUE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtTask::AddAnalysisObject(AlidNdPt *pObj) 
{
  // add analysis object to the list
  if(pObj == 0) {
    Printf("ERROR: Could not add comparison object");
    return kFALSE;
  }

  // add object to the list
  fCompList->AddLast(pObj);
       
return kTRUE;
}

//_____________________________________________________________________________
void AlidNdPtTask::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(1, "RECREATE");

  //
  // create output list
  //
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();

  // add dNdPt analysis objects to the output
  AlidNdPt *pObj=0;
  Int_t count=0;
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AlidNdPt *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
    count++;
  }
  Printf("UserCreateOutputObjects(): Number of output objects: %d \n", count);

  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTask::UserExec(Option_t *) 
{
  //
  // Called for each event
  //

  // ESD event
  fESD = (AliESDEvent*) (InputEvent());
  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }

  // MC event
  if(fUseMCInfo) {
    fMC = MCEvent();
    if (!fMC) {
      Printf("ERROR: MC event not available");
      return;
    }
  }


    AlidNdPt *pObj = 0;
    fPitList->Reset();
    while((pObj = (AlidNdPt *)fPitList->Next()) != NULL) {
      pObj->Process(fESD,fMC);
    }


  // Post output data.
  PostData(1, fOutput);
}

//_____________________________________________________________________________
void AlidNdPtTask::Terminate(Option_t *) 
{
  // Called one at the end 
  
  // check output data
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    Printf("ERROR: AlidNdPtTask::Terminate(): Output data not avaiable GetOutputData(0)==0x0 ..." );
    return;
  }
}
