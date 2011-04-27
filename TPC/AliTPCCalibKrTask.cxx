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

//----------------------------------------------------------------------------
// The AliTPCCalibKrTask class description (TPC Kr calibration).
// The AliTPCCalibKrTask loops over tree of TPC Kr clusters and fills AliTPCCalibKr 
// calibration component. 
// 
// As the input it requires the tree with reconstructed Kr clusters (AliTPCclusterKr objects). 
// The AliTPCCalibKr output calibration component contains an array of TH3F histograms which can be stored 
// in the ouptut file.
//
//  Author: Jacek Otwinowski (J.Otwinowski@gsi.de)

/*
 
// Usage example:
//

// -- Load toolkit
gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
AliXRDPROOFtoolkit tool;

// -- Make chain of files
TChain * chain = tool.MakeChain("KrClusters.txt","Kr","",1000,0);

// -- Create the analysis manager
AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

// -- Calibration component 
AliTPCCalibKr *calibObj = new AliTPCCalibKr;
calibObj->SetIrocHistogram(200,100,6000);
calibObj->SetOrocHistogram(200,100,5500);
calibObj->Init();

// -- Add task
AliTPCCalibKrTask *task = new AliTPCCalibKrTask;
task->SetInputChain(chain);
task->SetTPCCalibKr(calibObj);
mgr->AddTask(task);

// -- Attach input
cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
mgr->ConnectInput(task, 0, cInput);

// -- Attach output
cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer,"outHistFile.root");
mgr->ConnectOutput(task, 0, cOutput);

// -- Run analysis
mgr->InitAnalysis();
mgr->PrintStatus();
mgr->StartAnalysis("local", chain);

*/

// system includes
#include <cstdio>
using namespace std;

//Root includes
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TFile.h>
//AliRoot includes
#include "AliRawReader.h"
#include "AliRawReaderRoot.h"
#include "AliRawReaderDate.h"
#include "AliTPCRawStream.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliMathBase.h"
#include "TTreeStream.h"

//date
#include "event.h"

//header file
#include "AliTPCCalibKr.h"
#include "AliTPCCalibKrTask.h"

Int_t AliTPCCalibKrTask::fEvtNumber = 0;

ClassImp(AliTPCCalibKrTask)

AliTPCCalibKrTask::AliTPCCalibKrTask(const char *name) : 
  AliAnalysisTask(name,""),
  fClustKr(0),
  fTPCCalibKr(0),
  fTree(0),
  fOutput(0)
{
  //
  // default constructor
  //

  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
}

//_____________________________________________________________________
AliTPCCalibKrTask::~AliTPCCalibKrTask() 
{
  //
  // destructor
  //
  if(fOutput) fOutput->Delete();
  delete fOutput; fOutput = 0;
}

//_____________________________________________________________________
void AliTPCCalibKrTask::ConnectInputData(Option_t *)
{
  // Connect input data
  // Called once

  fTree = dynamic_cast<TTree*> (GetInputData(0));

  if(!fTree) { 
   Printf("ERROR: Could not read chain from input");
  }
  else {
   fTree->SetBranchStatus("*",1); 
   fTree->SetBranchStatus("Cl.fCluster",0); 
  }

  // set branch address
  if(!fTree->GetBranch("Cl.")) {
    Printf("ERROR: Could not get Cl. branch from input");
  } else {
   fTree->GetBranch("Cl.")->SetAddress(&fClustKr);
  }
}

//_____________________________________________________________________
void AliTPCCalibKrTask::CreateOutputObjects()
{
  // create object to the output 
  fOutput = new TList;
  fOutput->SetOwner(); // is owner of the fTPCCalibKr objects

  if(fTPCCalibKr) fOutput->Add(fTPCCalibKr);
  //fTPCCalibKr = new AliTPCCalibKr;
  //if(fTPCCalibKr) fOutput->Add(fTPCCalibKr);
  else
     Printf("WARNING: AliTPCCalibKr is not added to the output");
}

//_____________________________________________________________________
Bool_t AliTPCCalibKrTask::ReadEntry(Int_t evt)
{
  // read entry 
  Long64_t centry = fTree->LoadTree(evt);
  if(centry < 0) return kFALSE;

  if(!fTree->GetBranch("Cl.")) 
  {
    Printf("ERROR: Could not get Cl. branch from input");
	return kFALSE;
  } else {
   fTree->GetBranch("Cl.")->SetAddress(&fClustKr);
  }

  fTree->GetEntry(evt);

return kTRUE;
}

//_____________________________________________________________________
void AliTPCCalibKrTask::Exec(Option_t *)
{
  // Main loop
  // Called for each event
  
  // read entry
  if(fClustKr) delete fClustKr; fClustKr=0;
  Bool_t status = ReadEntry(fEvtNumber);
  if(status==kTRUE) 
  {
	  // Process output objects
      if(fClustKr) fTPCCalibKr->Process(fClustKr);
  }
 
  if( !( fEvtNumber % 100000) ) {
    cout << fEvtNumber << endl; }

  fEvtNumber++;

  // Post output data.
  PostData(0, fOutput);
}

//_____________________________________________________________________
void AliTPCCalibKrTask::Terminate(Option_t *) 
{
  // Called once at the end of the event loop
  cout << "Terminate " << endl;

  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {
    Printf("ERROR: fOutput not available");
    return;
  }
  
  fTPCCalibKr = dynamic_cast<AliTPCCalibKr*> (fOutput->FindObject("AliTPCCalibKr"));
  if (!fTPCCalibKr) {
    Printf("WARNING: AliTPCCalibKr not available");
    return;
  }

  for(int i=0; i<72; ++i) {
	 if( fTPCCalibKr->IsCSide(i) == kTRUE )
	   printf("C side chamber: %d, 3D histo entries: %10.f \n",i,fTPCCalibKr->GetHistoKr(i)->GetEntries());

	 if( fTPCCalibKr->IsCSide(i) == kFALSE )
	   printf("A side chamber: %d, 3D histo entries: %10.f \n",i,fTPCCalibKr->GetHistoKr(i)->GetEntries());
  }
}
