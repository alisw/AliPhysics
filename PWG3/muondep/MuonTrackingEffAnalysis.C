// 2008
// Macro for the running of the AliAnalysisTaskMuonTrackingEff
//

// ROOT includes
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>

// PWG3 includes
#include "AliAnalysisTaskMuonTrackingEff.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

// STEER includes
#include "AliESDInputHandler.h"
#include "AliCDBManager.h"

// MUON includes
#include "AliMUONGeometryTransformer.h"



void MuonTrackingEffAnalysis(const Bool_t alien = false, const Int_t run = 100, const char * fileName = "AliESDs.root", const char * geometryFileName = "geometry.root")
{
//Chain construction
    TChain *chain = new TChain("esdTree");
    chain->Add(fileName);

//Load OCD
    AliCDBManager* man = AliCDBManager::Instance();
    TString ocdbPath;
    if(alien)
    {
      ocdbPath = "alien://folder=/alice/data/2008/LHC08a/OCDB";
    }
    else
    {
      ocdbPath = "local://$ALICE_ROOT/OCDB";
    }

    man->SetDefaultStorage(ocdbPath.Data());
    man->SetSpecificStorage("MUON/Calib/Mapping",ocdbPath);
    man->SetRun(run);
 
//Load Geometry
    AliMUONGeometryTransformer* transformer = new AliMUONGeometryTransformer();
    transformer->LoadGeometryData(geometryFileName);


//Make analysis manager:
    AliAnalysisManager* mgr = new AliAnalysisManager("Manager", "Manager");  
    AliAnalysisTaskMuonTrackingEff* ESDTask = new AliAnalysisTaskMuonTrackingEff("ESDTask", transformer);
    AliESDInputHandler* inHandler = new AliESDInputHandler();

    mgr->SetInputEventHandler  (inHandler );
//     mgr->SetDebugLevel(10);
    mgr->AddTask(ESDTask);

//Create containers for input/output
    AliAnalysisDataContainer* cinput1  =
	mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("chistlist1", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, "MuonTrackingChamberEffHistos.root"); 

//Connection
    mgr->ConnectInput (ESDTask, 0, cinput1 );
    mgr->ConnectOutput(ESDTask, 0, coutput1);
   
//Run analysis
    if(mgr->InitAnalysis())
    {
      //mgr->PrintStatus();
      mgr->StartAnalysis("Local", chain);
    }
} 
