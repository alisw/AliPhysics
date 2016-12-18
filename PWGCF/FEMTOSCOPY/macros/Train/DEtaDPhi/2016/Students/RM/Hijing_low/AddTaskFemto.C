

void AddTaskFemto(TString configMacroName, TString containerName="femtolist", TString configMacroParameters="" )
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemto", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  if (TProofMgr::GetListOfManagers()->GetEntries()) {
    gProof->Load(configMacroName);
  }

  AliAnalysisTaskFemtoMJ *taskfemto = new AliAnalysisTaskFemtoMJ("TaskFemto","$ALICE_PHYSICS/"+configMacroName,configMacroParameters,kFALSE);//"$ALICE_PHYSICS/"+configMacroName,configMacroParameters,kTRUE);
  // TFile *filter_file = TFile::Open("alien:///alice/cern.ch/user/r/rmaselek/filters/first_approach/filters.root");
  mgr->AddTask(taskfemto);


  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer(containerName,  TList::Class(),
  							       AliAnalysisManager::kOutputContainer,outputfile);

  mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskfemto, 0, cout_femto);

  return;
}
