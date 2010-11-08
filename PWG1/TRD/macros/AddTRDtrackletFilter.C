AliAnalysisTask *AddTRDtrackletFilter()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  //  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/AliTRDtaskTrackletFilter.cxx++g");

  AliTRDonlineTrackletFilter *task = new AliTRDonlineTrackletFilter("tracklet-filter");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();

  AliAnalysisDataContainer *ctracklets =
    mgr->CreateContainer("ctracklets", TTree::Class(),
			 AliAnalysisManager::kOutputContainer,
			 "tracklets.root");

  task->ConnectInput(0, cinput);

  task->ConnectOutput(0, coutput);
  task->ConnectOutput(1, ctracklets);

  return task;
}
