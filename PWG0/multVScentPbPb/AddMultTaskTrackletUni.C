
AliTrackletTaskUni* AddMultTaskTrackletUni(const char* outName="tracklet.root", TString nomergeDir="")
{
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("My test train");
  // create our task
  AliTrackletTaskUni *task = new AliTrackletTaskUni("AliTrackletTaskUni");
  task->SetDontMerge(!nomergeDir.IsNull());
  // create output container
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist", TList::Class(),AliAnalysisManager::kOutputContainer,outName);
  if (!nomergeDir.IsNull()) coutput1->SetSpecialOutput();
  // add our task to the manager
  mgr->AddTask(task);

  // finaly connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  if (!nomergeDir.IsNull()) mgr->SetSpecialOutputLocation(nomergeDir.Data()); //"root://alicers01.cern.ch//tmp/myoutput/");
  //
  return task;
}
