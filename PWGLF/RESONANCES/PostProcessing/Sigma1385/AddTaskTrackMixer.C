AliAnalysisTaskTrackMixer* AddTaskTrackMixer(
    const char* taskname = "TrackMixer",
    const char* option = "",
    int nmix = 10,
    const char* suffix = "") {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }
  TString foption = option;
  AliAnalysisTaskTrackMixer* taskTrackMixer =
      new AliAnalysisTaskTrackMixer(Form("%s%s", taskname, suffix));
  taskTrackMixer->fEventCuts.fCentralityFramework = 1;
  taskTrackMixer->fEventCuts.SetMaxVertexZposition(10);
  taskTrackMixer->fEventCuts.SelectOnlyInelGt0(false);
  std::cout << "AddTaskTrackMixer:: Option: " << option << std::endl;
  if (foption.Contains("HM")) {
    taskTrackMixer->SetHighMult(kTRUE);  // default: kFALSE
    taskTrackMixer->fEventCuts.fTriggerMask =
        AliVEvent::kHighMultV0;  // default: kINT7
  }
  taskTrackMixer->SetnMix(nmix);

  if (!taskTrackMixer)
    return 0x0;
  mgr->AddTask(taskTrackMixer);

  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(taskTrackMixer, 0, cinput);

  AliAnalysisDataContainer* coutputTrackMixer = mgr->CreateContainer(
      Form("%s%s", taskname, suffix), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  mgr->ConnectOutput(taskTrackMixer, 1, coutputTrackMixer);

  return taskTrackMixer;
}
