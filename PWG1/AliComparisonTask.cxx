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
// Implementation of the AliComparisonTask class. It compares properties of the 
// reconstructed and MC particle tracks under several conditions. 
// As the input it requires the tree with AliRecInfo and AliMCInfo branches. Such
// tree can be prepared in advance by runing AliGenInfoMaker and then AliRecInfoMaker
// (details in description of these classes).
//
// The comparison output objects deriving from AliComparisonObject 
// (e.g. AliComparisonRes, AliComparisonEff, AliComparisonDEdxA, AliComparisonDCA ...) 
// are stored in the output file (details in description of these classes).
// 
// Author: J.Otwinowski 04/02/2008 
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
#include "AliESDInputHandler.h"
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
#include "AliComparisonTask.h"

using namespace std;

ClassImp(AliComparisonTask)

Int_t AliComparisonTask::fEvtNumber = 0;

//_____________________________________________________________________________
AliComparisonTask::AliComparisonTask(const char *name) 
  : AliAnalysisTask(name, "")
  , fTree(0)
  , fInfoMC(0)
  , fInfoRC(0)
  , fOutput(0)
  , fPitList(0)
  , fCompList(0)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AliComparisonTask::~AliComparisonTask()
{
  if(fOutput)   delete fOutput;  fOutput =0; 
  if(fCompList)   delete fCompList;  fCompList =0; 
}

//_____________________________________________________________________________
void AliComparisonTask::ConnectInputData(Option_t *) 
{
  // Connect input data 
  // Called once

  fTree = dynamic_cast<TTree*> (GetInputData(0));
  if (!fTree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    fTree->SetBranchStatus("*",1);
  }

  if(fTree->GetBranch("MC") &&  fTree->GetBranch("RC")) {
    fTree->GetBranch("MC")->SetAddress(&fInfoMC);
    fTree->GetBranch("RC")->SetAddress(&fInfoRC);
  } else {
      Printf("ERROR: Could not get MC and RC branches");
  }
}

//_____________________________________________________________________________
Bool_t AliComparisonTask::AddComparisonObject(AliComparisonObject *pObj) 
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
void AliComparisonTask::CreateOutputObjects()
{
  // Create histograms
  // Called once

  // create output list
  fOutput = new TList;
  fOutput->SetOwner();
  fPitList = fOutput->MakeIterator();

  AliComparisonObject *pObj=0;
  Int_t count=0;

  // add comparison objects to the output
  TIterator *pitCompList = fCompList->MakeIterator();
  pitCompList->Reset();
  while(( pObj = (AliComparisonObject *)pitCompList->Next()) != NULL) {
    fOutput->Add(pObj);
	count++;
  }
  Printf("CreateOutputObjects(): Number of output comparison objects: %d \n", count);
}

//_____________________________________________________________________________
Bool_t AliComparisonTask::ReadEntry(Int_t evt) 
{
// Read entry from the tree
  Long64_t centry = fTree->LoadTree(evt);
  if(centry < 0) return kFALSE;

  if(fTree->GetBranch("MC") &&  fTree->GetBranch("RC")) {
    fTree->GetBranch("MC")->SetAddress(&fInfoMC);
    fTree->GetBranch("RC")->SetAddress(&fInfoRC);
  } else {
      Printf("ERROR: Could not get MC and RC branches");
	  return kFALSE;
  }
  fTree->GetEntry(evt);

return kTRUE;
}
//_____________________________________________________________________________
void AliComparisonTask::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliComparisonObject *pObj=0;

  if (!fInfoMC && !fInfoRC) {
    Printf("ERROR: fInfoMC && fInfoRC not available");
    return;
  }

  // Process comparison
  Bool_t status = ReadEntry(fEvtNumber);
  if(status == kTRUE) 
  {
    fPitList->Reset();
    while(( pObj = (AliComparisonObject *)fPitList->Next()) != NULL) {
       pObj->Exec(fInfoMC,fInfoRC);
    }
  }

  if( !( fEvtNumber % 10000) ) { 
    cout << fEvtNumber << endl;
  }
  fEvtNumber++;

  // Post output data.
  PostData(0, fOutput);
}

//_____________________________________________________________________________
void AliComparisonTask::Terminate(Option_t *) 
{
  // Called one at the end 
  
  // check output data
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: AliComparisonTask::Terminate(): Output data not avaiable GetOutputData(0)==0x0 ..." );
    return;
  }
}
