#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#endif

void runCalibTrain(const char *inFileName = "AliESDs.root",
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

    gROOT->LoadMacro("AddTaskCopyESD.C");
    gROOT->LoadMacro("AddTaskFilterFriend.C");
    gROOT->LoadMacro("AddTaskFilterFriendSecond.C");
    gROOT->LoadMacro("AddTaskAddObject.C");

    AliLog::SetClassDebugLevel("AliESDEvent",19);
    TChain *chain = new TChain("esdTree");
    // Steering input chain
    chain->Add(inFileName);
    AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to ESD", "Analysis Manager");
    //mgr->SetDebugLevel(3);

    // Input
    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    inpHandler->SetActiveBranches("ESDfriend*");
    mgr->SetInputEventHandler  (inpHandler);

    // Output
    AliESDHandler* esdHandler   = new AliESDHandler();
    esdHandler->SetOutputFileName(outFileName);
    mgr->SetOutputEventHandler(esdHandler);

    // Tasks 
    AliAnalysisTaskCopyESD *copy = AddTaskCopyESD();
    
    AliAnalysisTaskFilterFriend* filter = AddTaskFilterFriend();

    AliAnalysisTaskFilterFriendSecond* filter2 = AddTaskFilterFriendSecond();
    
    AliAnalysisTaskAddObject* add = AddTaskAddObject();
    
    //
    // Run the analysis
    //
    if (!mgr->InitAnalysis()) {
	    printf("Analysis cannot be started, returning\n");
	    return;
    }

    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
   
    return;
}

