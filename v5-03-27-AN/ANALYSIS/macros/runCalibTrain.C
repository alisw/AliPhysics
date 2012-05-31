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

void runCalibTrain(TString runNumberString, const char *inFileName = "AliESDs.root",
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
    gSystem->Load("libPWGmuon");
    gROOT->LoadMacro("ConfigCalibTrain.C");
    gROOT->LoadMacro("AddTaskCopyESD.C");
    gROOT->LoadMacro("AddTaskFilterFriend.C");
    gROOT->LoadMacro("AddTaskFilterFriendSecond.C");
    gROOT->LoadMacro("AddTaskAddObject.C");

    // detector tasks

    gROOT->LoadMacro("AddTaskTPCCalib.C");

    AliLog::SetClassDebugLevel("AliESDEvent",19);
    TChain *chain = new TChain("esdTree");

    // Steering input chain

    chain->Add(inFileName);
    Int_t runNumber = runNumberString.Atoi();
    printf("runNumber from runCalibTrain = %d\n",runNumber);
    ConfigCalibTrain(runNumber, "raw://");

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
