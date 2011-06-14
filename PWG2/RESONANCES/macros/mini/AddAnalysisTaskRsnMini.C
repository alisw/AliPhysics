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

Bool_t usePhi   = 1;
Bool_t usePhiMC = 1;

AliRsnMiniAnalysisTask * AddAnalysisTaskRsnMini
(
   Bool_t      isMC,
   const char *path,
   Int_t       nmix = 0
)
{  
   //
   // -- INITIALIZATION ----------------------------------------------------------------------------
   //
   
   // retrieve analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

   // create the task and connect with physics selection
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask("RSN");
   mgr->AddTask(task);
   
   // settings
   task->UseMultiplicity("TRACKS");
   
   // set mixing
   task->SetNMix(nmix);
   task->SetMaxDiffVz(2.0);
   task->SetMaxDiffMult(10.0);
   task->SetMaxDiffAngle(1E20);
   
   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   
   // cut on primary vertex:
   // - 2nd argument --> |Vz| range
   // - 3rd argument --> minimum required number of contributors
   // - 4th argument --> tells if TPC stand-alone vertexes must be accepted
   AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
   
   // set the check for pileup
   cutVertex->SetCheckPileUp(kTRUE);
      
   // define and fill cut set
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutVertex);
   eventCuts->SetCutScheme(cutVertex->GetName());
   
   // set cuts in task
   task->SetEventCuts(eventCuts);
   
   //
   // -- CONFIGS -----------------------------------------------------------------------------------
   //
   
   if (usePhi) {
      gROOT->LoadMacro(Form("%s/ConfigPhi.C", path));
      if (!ConfigPhi(task, isMC, "default")) return 0x0;
   }
   
   if (isMC && usePhiMC) {
      gROOT->LoadMacro(Form("%s/ConfigPhiMC.C", path));
      if (!ConfigPhiMC(task, "default")) return 0x0;
   }
   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   
   const char *file = AliAnalysisManager::GetCommonFileName();
   AliAnalysisDataContainer *output = mgr->CreateContainer("RsnOut", TList::Class(), AliAnalysisManager::kOutputContainer, file);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   return task;
}
