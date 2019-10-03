//
// General macro to configure the RSN analysis task.
// It calls all configs desired by the user, by means
// of the boolean switches defined in the first lines.
// ---
// Inputs:
//  1) flag to know if running on MC or data
//  2) path where all configs are stored
// ---
// Returns:
//  kTRUE  --> initialization successful
//  kFALSE --> initialization failed (some config gave errors)
//

Bool_t useEvent       = 1;
Bool_t useTest        = 0;
Bool_t useKaonMonitor = 1;
Bool_t usePhiTPC      = 1;
Bool_t useKStarTPC    = 0;

Bool_t AddAnalysisTaskRsn
(
   Bool_t      isMC,
   Bool_t      useCentrality,
   Bool_t      checkPileUp,
   const char *path = "$(HOME)/code/resonances/alice-rsn-package/PWG2resonances/RESONANCES/macros/test/pulvir"
)
{  
   //
   // -- INITIALIZATION ----------------------------------------------------------------------------
   //
   
   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create the task and connect with physics selection
   AliRsnAnalysisTask *task = new AliRsnAnalysisTask("RSN");
   task->SelectCollisionCandidates();
   mgr->AddTask(task);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   
   // set the check for pileup
   ::Info("AddAnalysisTaskRsn", "Check of pile-up from SPD is %s", (checkPileUp ? "active" : "disabled"));
   cutVertex->SetCheckPileUp(checkPileUp);
      
   // define and fill cut set
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
   
   //
   // -- CONFIGS -----------------------------------------------------------------------------------
   //
   
   // add configs and count how many are added
   // will return kFALSE if none is added
   
   Int_t added = 0;
   
   if (useTest) {
      // NOTE: this is a test and will had its own definition of event cuts
      // ----> will not use those defined above
      gROOT->LoadMacro(Form("%s/RsnConfigTest.C", path));
      if (!RsnConfigTest(task, isMC)) return kFALSE;
      added++; 
   }
   
   if (useEvent) {
      gROOT->LoadMacro(Form("%s/RsnConfigEvent.C", path));
      if (!RsnConfigEvent(task, isMC, useCentrality, eventCuts)) return kFALSE;
      added++;
   }
   
   if (useKaonMonitor) {
      gROOT->LoadMacro(Form("%s/RsnConfigMonitorTPC.C", path));
      if (!RsnConfigMonitorTPC(task, isMC, useCentrality, eventCuts)) return kFALSE;
      added++;
   }
   
   if (usePhiTPC) {
      gROOT->LoadMacro(Form("%s/RsnConfigPhiTPC.C", path));
      if (!RsnConfigPhiTPC(task, isMC, useCentrality, eventCuts)) return kFALSE;
      added++;
   }
   
   if (useKStarTPC) {
      gROOT->LoadMacro(Form("%s/RsnConfigKStarTPC.C", path));
      if (!RsnConfigKStarTPC(task, isMC, useCentrality, eventCuts)) return kFALSE;
      added++;
   }
   
   ::Info("AddAnalysisTaskRsn.C", "Added %d configs", added);
   if (!added) {
      ::Error("AddAnalysisTaskRsn.C", "Cannot process an empty task!");
      return kFALSE;
   }
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   
   const char *file = AliAnalysisManager::GetCommonFileName();
   AliAnalysisDataContainer *output = mgr->CreateContainer("RsnOut", TList::Class(), AliAnalysisManager::kOutputContainer, file);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   return kTRUE;
}
