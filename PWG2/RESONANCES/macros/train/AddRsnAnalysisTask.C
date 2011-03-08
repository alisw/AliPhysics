//
// This macro serves to add the RSN analysis task to the steering macro.
//
// Inputs:
//   - dataLabel   = a string with informations about the type of data
//                   which could be needed to be ported to the config macro
//                   to set up some cuts
//   - configMacro = macro which configures the analysis; it has *ALWAYS*
//                   defined inside a function named 'RsnConfigTask()',
//                   whatever the name of the macro itself, whose first two
//                   arguments must have to be the task and the 'dataLabel' argument.
//
Bool_t AddRsnAnalysisTask
(
   Bool_t      isMC,
   Bool_t      isMix,
   const char *options,
   const char *path     = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train",
   const char *taskName = "RSNtask"
)
{
   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create the task and connect with physics selection
   AliRsnAnalysisTask *task = new AliRsnAnalysisTask(taskName);
   task->SetZeroEventPercentWarning(100.0);
   task->SelectCollisionCandidates();
   task->SetMixing(isMix);
   ::Info("AddRsnAnalysisTask.C", "Mixing: %s", (task->IsMixing() ? "YES" : "NO"));
   ::Info("AddRsnAnalysisTask.C", "MC:     %s", (isMC ? "YES" : "NO"));

   // add the task to manager
   mgr->AddTask(task);
   
   // add the event computations with the options to eventually select centrality cut
   gROOT->LoadMacro(Form("%s/AddRsnEventComputations.C", path));
   AddRsnEventComputations(isMC, options);
   
   // load common macro with cuts and axes
   // for cuts and axes, load the support macro
   gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PWG2/RESONANCES");
   gROOT->LoadMacro(Form("%s/CPhiCutsAndAxes.C++", path));
   
   // add all configs for phi
   gROOT->LoadMacro(Form("%s/RsnConfigPhi.C", path));
   RsnConfigPhi(isMC, "tpcpid_tofpid", path, taskName);
   //RsnConfigPhi(isMC, "itspid_tpcpid_tofpid", path, taskName);
   
   // in case of MC, add efficiency tasks
   if (isMC) {
      ::Info("Adding efficiency");
      gROOT->LoadMacro(Form("%s/AddRsnAnalysisTaskEffPhi.C", path));
      AddRsnAnalysisTaskEffPhi(options, "tpcpid_tofpid");
      //AddRsnAnalysisTaskEffPhi(options, "itspid_tpcpid_tofpid");
   }

   // connect input container according to source choice
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

   // create paths for the output in the common file
   Char_t commonPath[500];
   sprintf(commonPath, "%s", AliAnalysisManager::GetCommonFileName());

   // create containers for output
   AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfo", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
   AliAnalysisDataContainer *outputHist = mgr->CreateContainer("RsnHist", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
   mgr->ConnectOutput(task, 1, outputInfo);
   mgr->ConnectOutput(task, 2, outputHist);

   return kTRUE;
}
