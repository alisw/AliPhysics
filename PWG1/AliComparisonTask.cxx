//------------------------------------------------------------------------------
// Implementation of the AliComparisonTask class. It compares properties of the 
// reconstructed and MC particle tracks under several conditions. 
// As the input it requires the TTree with AliRecInfo and AliMCInfo branches. 
// 
// The comparison output objects deriving from AliComparisonObject 
// (e.g. AliComparisonRes, AliComparisonEff, AliComparisonDEdxA, AliComparisonDCA ...) 
// are stored in the Output.root file.
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
#include "AliMagFMaps.h"
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

Int_t AliComparisonTask::evtNumber = 0;

//_____________________________________________________________________________
AliComparisonTask::AliComparisonTask(const char *name) 
  : AliAnalysisTask(name, "")
  , fTree(0)
  , fInfoMC(0)
  , fInfoRC(0)
  , fOutput(0)
  , fMagField(0)
  , fMagFMap(0)
  , pitList(0)
  , fCompList(0)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  // set default mag. field
  SetMagField();
  
  // create the list for comparison objects
  fCompList = new TList;
}

//_____________________________________________________________________________
AliComparisonTask::~AliComparisonTask()
{
  if(fOutput)   delete fOutput;  fOutput =0; 
  if(fMagFMap)  delete fMagFMap;  fMagFMap =0; 
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
  
  // set mag. field map 
  fMagFMap = new AliMagFMaps("Maps","Maps", 2, 1., 10., fMagField);
  AliTracker::SetFieldMap(fMagFMap,kFALSE);
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
  pitList = fOutput->MakeIterator();

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
  Bool_t status = ReadEntry(evtNumber);
  if(status == kTRUE) 
  {
    pitList->Reset();
    while(( pObj = (AliComparisonObject *)pitList->Next()) != NULL) {
       pObj->Exec(fInfoMC,fInfoRC);
    }
  }

  if( !( evtNumber % 10000) ) { 
    cout << evtNumber << endl;
  }
  evtNumber++;

  // Post output data.
  PostData(0, fOutput);
}

//_____________________________________________________________________________
void AliComparisonTask::Terminate(Option_t *) 
{
  // Called once at the end of the event loop
  TFile *out = new TFile("Output.root","RECREATE");
  out->cd();

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fOutput->Write();
  out->Close();
}
