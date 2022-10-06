AliAnalysisTaskSpectrumT0* AddAnalysisTaskSpectrumT0(TString name, TString filepathRunTree, TString eventIDs="1011010010 10011010010 1",TString eventIDsStat="1011010010 10011010010 1") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    Error("AddTask", "This task requires an input event handler");
    return NULL;
  }

  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":taskSpectrumT0";
  AliAnalysisTaskSpectrumT0* task = new AliAnalysisTaskSpectrumT0(name.Data());
  bool isRunListReady = task->GetRunNumsFromROOT(filepathRunTree);
  if(isRunListReady==false) {
    Error("AddTask", "Cannot find run list located in tree");
    return NULL;
  }
  const TString stDel = " ";
  bool isEventIDsReady=false;
  TObjArray *arrayEventIDs = eventIDs.Tokenize(stDel);
  for(const auto& en:(*arrayEventIDs)) {
    TObjString *objst = dynamic_cast<TObjString *>(en);
    if(objst==NULL) continue;
    TString st = objst->GetString();
    if(st.IsBin()) {
      task->AddEventID(std::string{st});
      isEventIDsReady=true;
    }
  }
  if(!isEventIDsReady) {
    Error("AddTask", "Cannot add EventIDs");
    return NULL;
  }

  bool isEventIDsStatReady=false;
  TObjArray *arrayEventIDsStat = eventIDsStat.Tokenize(stDel);
  for(const auto& en:(*arrayEventIDsStat)) {
    TObjString *objst = dynamic_cast<TObjString *>(en);
    if(objst==NULL) continue;
    TString st = objst->GetString();
    if(st.IsBin()) {
      task->AddEventIDstat(std::string{st});
      isEventIDsStatReady=true;
    }
  }
  if(!isEventIDsStatReady) {
    Error("AddTask", "Cannot add EventIDs for stat hists");
    return NULL;
  }


  mgr->AddTask(task);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,mgr->CreateContainer("output", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  return task;
}
