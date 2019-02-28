/**
 * @file   FTAddMyTask.C
 * @author Freja Thoresen <freja.thoresen@cern.ch>
 *
 * @brief  Add Q-cummulant forward task to train
 *
 *
 * @ingroup pwglf_forward_scripts_tasks
 */
/**
 * @defgroup pwglf_forward_flow Flow
 *
 * Code to deal with flow
 *
 * @ingroup pwglf_forward_topical
 */
/**
 * Add Flow task to train
 *
 * @ingroup pwglf_forward_flow
 */
 #include <string>

 #include "TROOT.h"
 #include "TSystem.h"

 #include "AliAODHandler.h"
 #include "AliAnalysisDataContainer.h"
 #include "AliAnalysisManager.h"
 #include "AliOADBPhysicsSelection.h"
 #include "AliPhysicsSelectionTask.h"
 //#include "AddTaskForwardSecondaries.h"
 #include "AliVEvent.h"
 #include "AliAnalysisDataSlot.h"
 #include "AliAnalysisTaskMCSmearing.h"
 #include "AliMultSelectionTask.h"
 #include "AliCentralMCMultiplicityTask.h"
 #include <sstream>
 #include "AliAnalysisTaskMCParticleFilter.h"
 #include "AliForwardMCMultiplicityTask.h"
 #include "AliAODHandler.h"
 #include "AliESDInputHandler.h"
 #include "AliMCEventHandler.h"
 #include "AliVEventHandler.h"
 #include "AliForwardMCMultiplicityTask.h"
 #include "AliForwardCorrectionManager.h"
 #include "TProof.h"
 #include "TFile.h"
 #include "AliPhysicsSelectionTask.h"
 #include "AliCentralitySelectionTask.h"
/**
* Set a parameter on the forward task
*
* @param task Task to set on
* @param call The call
* @param arg  The argument(s) for the call
*/
void SetOnMCFwd(AliAnalysisTask* task,
    const TString& call,
    const TString& arg)
{
  TString cmd;
  cmd.Form("((AliForwardMCMultiplicityTask*)%p)->%s(%s)",
     task, call.Data(), arg.Data());
  gROOT->ProcessLine(cmd);
}


AliAnalysisTaskSE* AddTaskForwardMakeAOD()
{
  std::cout << "______________________________________________________________________________" << std::endl;
  std::cout << "AddTaskForwardMakeAOD" << std::endl;

  // --- Get analysis manager ----------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    Fatal("","No analysis manager to connect to.");

  std::cout << "______________________________________________________________________________" << std::endl;

  // Add mult selection Task
  AliAnalysisManager::SetCommonFileName("forward.root");

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

  AliMultSelectionTask *taskMultSelection = new AliMultSelectionTask("taskMultSelection", "", kFALSE, 1);
  taskMultSelection->SetAlternateOADBforEstimators("LHC15o");

  mgr->AddTask(taskMultSelection);
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  outputFileName += ":MultSelection";
  if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

  Printf("Set OutputFileName : \n %s\n", outputFileName.Data());

  AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cListMultSelection",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );

  //Recommendation: Tree as a single output slot
  mgr->ConnectInput (taskMultSelection, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskMultSelection, 1, coutputList);

    // ______________________________________________________________________________

  Bool_t usePhysicsSelection = kTRUE;
  Bool_t useCentrality = kTRUE;

  if (usePhysicsSelection) {
  // Physics selection task
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
     mgr->RegisterExtraFile("event_stat.root");
     //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kTRUE);

     AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gROOT->ProcessLine("AddTaskPhysicsSelection(kTRUE)"));
     mgr->AddStatisticsTask(AliVEvent::kAny);
   }
  // Centrality (only Pb-Pb)
  if (useCentrality) {
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
     AliCentralitySelectionTask* taskCentrality = reinterpret_cast<AliCentralitySelectionTask*>(gROOT->ProcessLine("AddTaskCentrality()"));
     taskCentrality->SetMCInput();
  }

  std::cout << "______________________________________________________________________________" << std::endl;


  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");

  AliForwardMCMultiplicityTask* mytask = reinterpret_cast<AliForwardMCMultiplicityTask*>(gROOT->ProcessLine("AddTaskForwardMult(kTRUE,0,0,0,0,\"ForwardAODConfig.C\",\"$ALICE_PHYSICS/OADB/PWGLF/FORWARD/CORRECTIONS/data\",0)"));
    //AliForwardMCMultiplicityTask* mytask  = AddTaskForwardMult(kTRUE,0,0,0,0,\"ForwardAODConfig.C\",\"$ALICE_PHYSICS/OADB/PWGLF/FORWARD/CORRECTIONS/data\",0);
  SetOnMCFwd(mytask, "GetTrackDensity().SetMaxConsequtiveStrips", Form("%d",2));


  std::cout << "______________________________________________________________________________" << std::endl;



    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskCentralMult.C");

    AliCentralMCMultiplicityTask* mytask = reinterpret_cast<AliCentralMCMultiplicityTask*>(gROOT->ProcessLine("AddTaskCentralMult(kTRUE,0,0,0,0,\"CentralAODConfig.C\",\"$ALICE_PHYSICS/OADB/PWGLF/FORWARD/CORRECTIONS/data\",0)"));
    

    std::cout << "______________________________________________________________________________" << std::endl;



  AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
  mgr->AddTask(kinefilter);

  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
  coutput1->SetSpecialOutput();

  mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (kinefilter,  0, coutput1);
  AliAnalysisDataContainer *coutputEx = mgr->CreateContainer("cFilterList", TList::Class(),
								   AliAnalysisManager::kOutputContainer,"pyxsec_hists.root");
  coutputEx->SetDataOwned(kTRUE);
  mgr->ConnectOutput (kinefilter,  1,coutputEx);

  std::cout << "______________________________________________________________________________" << std::endl;

  return mytask;
}
/*
 * EOF
 *
 */
