#include "AliAnalysisMultPt.h"
void runLocal()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    //Bool_t local = kTRUE;
     //if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    //Bool_t gridTest = kTRUE;
    
    
    // header location
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisMultPt");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // compile the class (locally) with debug symbols
    gInterpreter->LoadMacro("AliAnalysisMultPt.cxx++g");

  
    // load the addtask macro and create the task
    AliAnalysisMultPt *task = reinterpret_cast<AliAnalysisMultPt*>(gInterpreter->ExecuteMacro("macros/AddMultPt.C"));

    if(!task) return 0x0;
    // Add task
    if(!mgr->InitAnalysis()) return;
    //mgr->AddTask(task);
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    //chain->Add("/Users/neginav/AliAODPbMC.root");
    //chain->Add("/Users/neginav/AliAODMC.root");
    chain->Add("/Users/neginav/AliAOD.root");
    //chain->Add("/Users/neginav/AliAOD0001.root");
    //chain->Add("/Users/neginav/AliAOD0003.root");
    //chain->Add("/Users/neginav/AliAOD0004.root");
    
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain); 
    //chain->Draw("aodTree");
}
