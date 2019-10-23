#ifndef __CLING__
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisNanoAODCuts.h"
#include "AliAnalysisTaskNanoAODFilter.h"
#include "AliAnalysisTaskNanoAODnormalisation.h"
#include "AliAnalysisTaskNanoAODskimming.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliESDtrackCuts.h"
#include "AliNanoAODTrack.h"
#include "AliNanoFilterPID.h"
#include "AliNanoSkimmingPID.h"
#include "AliPhysicsSelectionTask.h"

#include <TChain.h>
#include <TInterpreter.h>
#include <TStopwatch.h>
#endif

void runOnNanoRsn(const char* filePath = "AliAOD.NanoAOD.root") {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
        mgr = new AliAnalysisManager("train");

    AliAODInputHandler* iH = new AliAODInputHandler();
    mgr->SetInputEventHandler(iH);

    AliAnalysisTaskSE* pidRespTask =
        reinterpret_cast<AliAnalysisTaskSE*>(gInterpreter->ExecuteMacro(
            "$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));

    /*
    AliAnalysisTaskSE* task = (AliAnalysisTaskSE*)gInterpreter->ExecuteMacro(
        "$ALICE_PHYSICS/PWGLF/RESONANCES/PostProcessing/Xi1530/"
        "AddTaskXi1530.C");
        */
    gInterpreter->ExecuteMacro(
        "$ALICE_PHYSICS/PWGLF/RESONANCES/PostProcessing/Xi1530/NanoCheck/"
        "AddTaskNanoCheck.C");

    AliAnalysisTaskNanoAODnormalisation::AddTask();

    mgr->InitAnalysis();
    mgr->PrintStatus();

    // NOTE enable this for strict mode in which you get an AliFatal also for
    // fields in the header which are accessed but not available
    // AliNanoAODHeader::SetFatalMode();

    // Create chain of input files
    TChain* chain = new TChain("aodTree");
    chain->Add(filePath);

    TStopwatch watch;
    //AliLog::SetGlobalLogLevel(AliLog::kFatal);
    mgr->StartAnalysis("local", chain, 100000);
    watch.Print();
    gSystem->Exec("mv AnalysisResults.root AnalysisResultsNano.root");
}
