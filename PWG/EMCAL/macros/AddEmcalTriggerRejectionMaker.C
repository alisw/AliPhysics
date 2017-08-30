PWG::EMCAL::AliEmcalTriggerRejectionMaker *AddEmcalTriggerRejectionMaker(const char *dummy = "", const char *suffix = ""){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString taskname(Form("MaxPatchTask%s", suffix));
  PWG::EMCAL::AliEmcalTriggerRejectionMaker *patchtask = new PWG::EMCAL::AliEmcalTriggerRejectionMaker(taskname.Data());
  mgr->AddTask(patchtask);

  TString outputfile(mgr->GetCommonFileName());
  outputfile += TString::Format(":MaxPatch%s", suffix);

  mgr->ConnectInput(patchtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(patchtask, 1, mgr->CreateContainer(Form("MaxPatchData%s", suffix), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data()));

  return patchtask;
}
