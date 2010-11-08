AliAnalysisTask *AddTRDtrackletQA()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  //  gROOT->LoadMacro("$ALICE_ROOT/PWG1/TRD/AliTRDtrklQA.cxx++g");

  AliTRDonlineTrackletQA *task = new AliTRDonlineTrackletQA("trklQA");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();
  AliAnalysisDataContainer *ctracklets =
    (AliAnalysisDataContainer*) AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("ctracklets");

  AliAnalysisDataContainer *ctrklqa =
    mgr->CreateContainer("ctrklqa", TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 "trkl-qa.root");

  task->ConnectInput(0, cinput);
  task->ConnectInput(1, ctracklets);

  task->ConnectOutput(0, coutput);
  task->ConnectOutput(1, ctrklqa);

  return task;
}
