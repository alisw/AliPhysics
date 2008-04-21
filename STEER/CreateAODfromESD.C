#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TChain.h>
#include <TSystem.h>
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#endif

void CreateAODfromESD(const char *inFileName = "AliESDs.root",
		      const char *outFileName = "AliAOD.root") {
  
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    TChain *chain = new TChain("esdTree");

    // Steering input chain
    chain->Add(inFileName);
    AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to AOD", "Analysis Manager");

    // Input
    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler  (inpHandler);

    // Output
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName(outFileName);
    mgr->SetOutputEventHandler(aodHandler);

    // Task
    AliAnalysisTaskESDfilter *filter = new AliAnalysisTaskESDfilter("Filter");
    mgr->AddTask(filter);

    
    // Pipelining
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain", TChain::Class(),
                                                             AliAnalysisManager::kInputContainer);
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              "default");
  

    mgr->ConnectInput (filter, 0, cinput1 );
    mgr->ConnectOutput(filter, 0, coutput1);

    //
    // Run the analysis
    //
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
}
