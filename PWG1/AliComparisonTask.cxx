//------------------------------------------------------------------------------
// Implementation of the AliComparisonTask class. It compares properties of the 
// reconstructed and MC particle tracks under several conditions. 
// As the input it requires the TTree with AliRecInfo and AliMCInfo branches. 
// The comparison output histograms are stored 
// in the comparison objects: AliComparisonRes, AliComparisonEff, 
// AliComparisonDEdx and AliComparisonDCA. Each of these objects also contains 
// selection cuts which were used during filling the histograms.
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
  , fCompRes(0)
  , fCompEff(0)
  , fCompDEdx(0)
  , fCompDCA(0)
  , fOutput(0)
  , fMagField(0)
  , fMagFMap(0)
  , fGeom(0)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());

  // set default mag. field
  SetMagField();
  
  // set default geometry
  SetGeometry();
}

//_____________________________________________________________________________
AliComparisonTask::~AliComparisonTask()
{
  if(fOutput)   delete fOutput;  fOutput =0; 
  if(fMagFMap)  delete fMagFMap;  fMagFMap =0; 
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

  // set geommetry
  AliGeomManager::LoadGeometry(fGeom);
}

//_____________________________________________________________________________
void AliComparisonTask::CreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutput = new TList;
  fOutput->SetOwner();

  if(fCompRes) fOutput->Add(fCompRes);
  else 
     Printf("WARNING: AliComparisonRes is not added to the output");

  if(fCompEff) fOutput->Add(fCompEff);
  else 
    Printf("WARNING: AliComparisonEff is not added to the output");

  if(fCompDEdx) fOutput->Add(fCompDEdx);
  else 
    Printf("WARNING: AliComparisonDEdx is not added to the output");

  if(fCompDCA) fOutput->Add(fCompDCA);
  else 
     Printf("WARNING: AliComparisonDCA is not added to the output");
}

//_____________________________________________________________________________
Bool_t AliComparisonTask::ReadEntry(Int_t evt) 
{
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

  if (!fInfoMC && !fInfoRC) {
    Printf("ERROR: fInfoMC && fInfoRC not available");
    return;
  }

  // Process comparison
  Bool_t status = ReadEntry(evtNumber);
  if(status == kTRUE) 
  {
     if(fCompRes)  fCompRes->Exec(fInfoMC,fInfoRC);
     if(fCompEff)  fCompEff->Exec(fInfoMC,fInfoRC);
     if(fCompDEdx) fCompDEdx->Exec(fInfoMC,fInfoRC);
     if(fCompDCA)  fCompDCA->Exec(fInfoMC,fInfoRC);
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
  cout << "Terminate " << endl;

  TFile *out = new TFile("Output.root","RECREATE");
  out->cd();

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }

  fCompRes = dynamic_cast<AliComparisonRes*> (fOutput->FindObject("AliComparisonRes"));
  if (!fCompRes) {
    Printf("WARNING: AliComparisonRes not available");
    return;
  }

  fCompEff = dynamic_cast<AliComparisonEff*> (fOutput->FindObject("AliComparisonEff"));
  if (!fCompEff) {
    Printf("WARNING: AliComparisonEff not available");
    return;
  }
   
  fCompDEdx = dynamic_cast<AliComparisonDEdx*> (fOutput->FindObject("AliComparisonDEdx"));
  if (!fCompDEdx) {
    Printf("WARNING: AliComparisonDEdx not available");
    return;
  }

  fCompDCA = dynamic_cast<AliComparisonDCA*> (fOutput->FindObject("AliComparisonDCA"));
  if (!fCompDCA) {
    Printf("WARNING: AliComparisonDCA not available");
    return;
  }

  fOutput->Write();
  out->Close();
}
