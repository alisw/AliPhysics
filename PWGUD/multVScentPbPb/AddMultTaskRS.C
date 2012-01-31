AliTrackletTaskUni* AddMultTaskRS(const char* outName="tracklet.root")
{
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("My test train");
  // create our task
  AliTrackletTaskUni *task = new AliTrackletTaskUni("AliTrackletTaskUni");
  // create output container
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(),AliAnalysisManager::kOutputContainer,outName);

  // add our task to the manager
  mgr->AddTask(task);

  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  //
  return task;
}
