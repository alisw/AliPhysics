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

AliRsnMiniAnalysisTask * AddAnalysisTaskPhiRAApPb
(
   Bool_t      isMC,
   Bool_t      isPP,
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
   AliRsnMiniAnalysisTask *task = new AliRsnMiniAnalysisTask("RSN", isMC);
   task->UseESDTriggerMask(AliVEvent::kINT7);
   
   // settings
   if (isPP)
      task->UseMultiplicity("QUALITY");
   else
      task->UseCentrality("V0A");
   
   // set mixing
   task->UseContinuousMix();
   //task->UseBinnedMix();
   task->SetNMix(nmix);
   task->SetMaxDiffVz(5.0);
   task->SetMaxDiffMult(20.0);
   task->SetMaxDiffAngle(1E20);

   mgr->AddTask(task);

   //
   // -- EVENT CUTS (same for all configs) ---------------------------------------------------------
   //
   
   AliRsnCutEventUtils *cutEventUtils = new AliRsnCutEventUtils("cutEventUtils", kTRUE, kTRUE);
   cutEventUtils->SetUseVertexSelection2013pA(kTRUE);
   cutEventUtils->SetMinPlpContribSPD(5);


      
   // define and fill cut set
   AliRsnCutSet *eventCuts = new AliRsnCutSet("eventCuts", AliRsnTarget::kEvent);
   eventCuts->AddCut(cutEventUtils);
   eventCuts->SetCutScheme(cutEventUtils->GetName());
   
   // set cuts in task
   task->SetEventCuts(eventCuts);
   
   //
   // -- EVENT-ONLY COMPUTATIONS -------------------------------------------------------------------
   //
   
   // initialize value computation for multiplicity/centrality
   // second argument tells if the value must be taken from MC
   // (when this can be done)
   // after creating the value, the task returns its ID
   Int_t multID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   
   // create event-related output
   AliRsnMiniOutput *outMult = task->CreateOutput("eventMult", "HIST", "EVENT");
   // set axes, by passing value ID and defining the binning
   if (isPP) 
      outMult->AddAxis(multID, 300, 0.0, 300.0);
   else
      outMult->AddAxis(multID, 100, 0.0, 100.0);
   
   //
   // -- PAIR CUTS (common to all resonances) ------------------------------------------------------
   //
   
   AliRsnCutMiniPair *cutY = new AliRsnCutMiniPair("cutRapidity", AliRsnCutMiniPair::kRapidityRange);
   cutY->SetRangeD(-0.765, -0.165);

   AliRsnCutMiniPair *cutY2 = new AliRsnCutMiniPair("cutRapidity2", AliRsnCutMiniPair::kRapidityRange);
   cutY2->SetRangeD(-0.465, 0.035);
   
   AliRsnCutSet *cutsPair = new AliRsnCutSet("pairCuts", AliRsnTarget::kMother);
   cutsPair->AddCut(cutY);
   cutsPair->SetCutScheme(cutY->GetName());
   
   AliRsnCutSet *cutsPair2 = new AliRsnCutSet("pairCuts2", AliRsnTarget::kMother);
   cutsPair2->AddCut(cutY2);
   cutsPair2->SetCutScheme(cutY2->GetName());

   //
   // -- CONFIGS -----------------------------------------------------------------------------------
   //
   

   gROOT->LoadMacro(Form("%s/ConfigPhiRAApPb.C", path));
   if (!ConfigPhiRAApPb(task, isMC, isPP, "", cutsPair, cutsPair2)) return 0x0;

   
   //
   // -- CONTAINERS --------------------------------------------------------------------------------
   //
   
   const char *file = AliAnalysisManager::GetCommonFileName();
   AliAnalysisDataContainer *output = mgr->CreateContainer("RsnOut", TList::Class(), AliAnalysisManager::kOutputContainer, file);
   mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 1, output);

   return task;
}
