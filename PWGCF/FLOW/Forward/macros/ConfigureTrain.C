#include <string>

#include "TROOT.h"

#include "AliAODHandler.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskValidation.h"
#include "AliOADBPhysicsSelection.h"
#include "AliPhysicsSelectionTask.h"
#include "AliVEvent.h"
#include "AliAnalysisDataSlot.h"
#include "AliForwardFlowRun2Task.h"

enum mode {kRECON, kTRUTH};

void ConfigureTrain() {
  
  // Add mult selection Task
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  gROOT->ProcessLine("AddTaskMultSelection()");

  // PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask* ps = reinterpret_cast<AliPhysicsSelectionTask*>
  // Signature: Bool_t mCAnalysisFlag, Bool_t applyPileupCuts
  //(gROOT->ProcessLine("AddTaskPhysicsSelection(false, true)"));


 // gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/c2/AddTaskValidation.C");
  //AliAnalysisTaskValidation* validation_task =
    //reinterpret_cast<AliAnalysisTaskValidation*>(gROOT->ProcessLine("AddTaskValidation()"));

  //AliAnalysisTaskValidation* ev_sel_task =
    //dynamic_cast<AliAnalysisTaskValidation*>(AliAnalysisManager::GetAnalysisManager()
      //         ->GetExchangeContainers()
        //       ->FindObject("event_selection_xchange"));


  // Add the cumulant tasks last
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/thoresen/AddTaskForwardFlowRun2.C");
  AliForwardFlowRun2Task* mytask =
    reinterpret_cast<AliForwardFlowRun2Task*>
    (gROOT->ProcessLine("AddTaskForwardFlowRun2()"));

  // Add the C2 tasks last
 // mytask->ConnectInput(1, validation_task->GetOutputSlot(2)->GetContainer());
  //t->ConnectInput(1, validation_task->GetOutputSlot(2)->GetContainer());
}