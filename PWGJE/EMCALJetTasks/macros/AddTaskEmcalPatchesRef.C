EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *AddTaskEmcalPatchesRef(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef *task = new EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef("emcalPatchQA");

  // Set Energy thresholds for additional patch selection:
  // These are events with offline patches of a given type where the trigger reached already the plateau
  // These numers are determined as:
  // EMC7: 3.5 GeV
  // EG1:  14 GeV
  // EG2:  8 GeV
  // EJ1:  22 GeV
  // EJ2:  12 GeV
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef::kEPREL0, 5);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef::kEPREG1, 14);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef::kEPREG2, 8);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef::kEPREJ1, 22);
  task->SetOfflineEnergyThreshold(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalPatchesRef::kEPREJ2, 12);

  mgr->AddTask(task);

  TString outfile(mgr->GetCommonFileName());
  outfile += ":PatchQA";

  task->ConnectInput(0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer("PatchResults", TList::Class(), AliAnalysisManager::kOutputContainer, outfile.Data()));

  return task;
}
