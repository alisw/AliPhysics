/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisMultPt:
// Description: Analysis task to get multiplicity
// and pT distributions
// Author: Negin Alizadehvandchali
// (negin.alizadehvandchali@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "AliAnalysisMultPt.h"
void runLocalMultPt()
{
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
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    chain->Add("/Users/neginav/AODs for Pb-MC/AliAOD1.root");

    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
}
