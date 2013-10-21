AliFemtoK0Analysis *AddTaskK0Femto(bool FieldPositive = kTRUE, bool OnlineCase = kTRUE, bool MeritCase = kTRUE, float MinDL = 0.0, int MeritCutChoice = 4, float MinSep = 5.0, TString nameSpec = "NoSpec"){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  AliFemtoK0Analysis *K0Task = new AliFemtoK0Analysis("K0Task", FieldPositive, OnlineCase, MeritCase, MinDL, MeritCutChoice, MinSep);
  if(!K0Task) exit(-1);
  mgr->AddTask(K0Task);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGCF.outputK0Analysis_";
  outputFileName += nameSpec;
  if(FieldPositive) outputFileName += "_FieldPos.root";
  else outputFileName += "_FieldNeg.root";
  AliAnalysisDataContainer *coutputK0 = mgr->CreateContainer("MyList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());

  mgr->ConnectInput(K0Task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(K0Task, 1, coutputK0);

  return K0Task;
}


  
