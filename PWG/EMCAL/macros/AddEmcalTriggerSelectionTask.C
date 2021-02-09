PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection *AddEmcalTriggerSelectionTask(TString suffix = "") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddEmcalTriggerSelectionTask", "No analysis manager available");
    return NULL;
  }

  TString taskname("EmcalTriggerSelectionTask"), outfilename(TString::Format("%s:EmcalTriggerSelectionTask", mgr->GetCommonFileName()));

  if(suffix.Length()){
    taskname += "_" + suffix;
    outfilename += "_" + suffix;
  }
  PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection *trgseltask = new PWG::EMCAL::AliAnalysisTaskEmcalTriggerSelection(taskname);
  mgr->AddTask(trgseltask);

  mgr->ConnectInput(trgseltask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(trgseltask, 1, mgr->CreateContainer(taskname, AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfilename));

  return trgseltask;
}
