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



void MuonTrackingEffAnalysis(const Bool_t alien = false, const Int_t run = 100, const char * fileName = "AliESDs.root", const char * geometryFileName = "geometry.root", Bool_t isCosmicData = kFALSE)
{
//Chain construction
    TChain *chain = new TChain("esdTree");
    chain->Add(fileName);


    // Name of the output file
    const Char_t* outputName = "MuonTrackingChamberEffHistos.root";

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
    AliAnalysisTaskMuonTrackingEffComplete* ESDTask = new AliAnalysisTaskMuonTrackingEffComplete("ESDTask", transformer, isCosmicData);
    AliESDInputHandler* inHandler = new AliESDInputHandler();

    mgr->SetInputEventHandler  (inHandler );
//     mgr->SetDebugLevel(10);
    mgr->AddTask(ESDTask);

//Create containers for input/output
    AliAnalysisDataContainer* cinput0  =	mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer("TracksDetectedPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);
    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer("TotalTracksPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);
    AliAnalysisDataContainer *coutput2 =
	mgr->CreateContainer("EfficiencyPerDE", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);
    AliAnalysisDataContainer *coutput3 =
	mgr->CreateContainer("TracksDetectedPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);
    AliAnalysisDataContainer *coutput4 =
	mgr->CreateContainer("TotalTracksPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);
    AliAnalysisDataContainer *coutput5 =
	mgr->CreateContainer("EfficiencyPerChamber", TClonesArray::Class(),AliAnalysisManager::kOutputContainer, outputName);

//Connection
    mgr->ConnectInput (ESDTask, 0, cinput0 );

    mgr->ConnectOutput(ESDTask, 0, coutput0);
    mgr->ConnectOutput(ESDTask, 1, coutput1);
    mgr->ConnectOutput(ESDTask, 2, coutput2);
    mgr->ConnectOutput(ESDTask, 3, coutput3);
    mgr->ConnectOutput(ESDTask, 4, coutput4);
    mgr->ConnectOutput(ESDTask, 5, coutput5);
  
//Run analysis
    if(mgr->InitAnalysis())
    {
      //mgr->PrintStatus();
      mgr->StartAnalysis("Local", chain);
    }
} 
