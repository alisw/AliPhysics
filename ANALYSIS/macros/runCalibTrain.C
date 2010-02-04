/*
  Template of calibration/filtering  macro using ESD

  .L $ALICE_ROOT/ANALYSIS/macros/runCalibTrain.C
  runCalibTrain(105160);
 
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#endif

void runCalibTrain(Int_t runNumber, const char *inFileName = "AliESDs.root",
		   const char *outFileName = "AliESDs_v1.root") 
{
  
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSIScalib");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWG3muon");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/ConfigCalibTrain.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCopyESD.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskFilterFriend.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskFilterFriendSecond.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskAddObject.C");
    
    // detector tasks

    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AddTaskTPCCalib.C");

    AliLog::SetClassDebugLevel("AliESDEvent",19);
    TChain *chain = new TChain("esdTree");

    // Steering input chain

    chain->Add(inFileName);    
    ConfigCalibTrain(runNumber);

    AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to ESD", "Analysis Manager");
    // mgr->SetDebugLevel(3);

    // Input

    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    inpHandler->SetActiveBranches("ESDfriend*");
    mgr->SetInputEventHandler  (inpHandler);

    // Output

    AliESDHandler* esdHandler   = new AliESDHandler();
    esdHandler->SetOutputFileName(outFileName);
    mgr->SetOutputEventHandler(esdHandler);

    // Steering Tasks 

    AliAnalysisTaskCopyESD *copy = AddTaskCopyESD();
    AliAnalysisTaskFilterFriend* filter = AddTaskFilterFriend();
    AliAnalysisTaskFilterFriendSecond* filter2 = AddTaskFilterFriendSecond();
    AliAnalysisTaskAddObject* add = AddTaskAddObject();
   
    // Detector Tasks
   
    AliAnalysisTask* tTPC = AddTaskTPCCalib(runNumber);

    // Run the analysis

    if (!mgr->InitAnalysis()) {
	    printf("Analysis cannot be started, returning\n");
	    return;
    }

    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
   
    return;
}

