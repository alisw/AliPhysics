
void AddTaskFPexample(AliAnalysisManager *mgr,
		      AliAnalysisAlien *plugin,
		      const char *runtype,
		      const bool useRealData,
		      const char *taskname,
		      const Int_t gridRun = -1)
{
  if (!mgr) {
    Printf("ERROR: undefined manager, FPexample won't be added");
    return;
  }
  if (!plugin) {
    Printf("ERROR: undefined alien plugin, FPexample won't be added");
    return;
  }
  // create task
  gROOT->LoadMacro("AliAnalysisTaskFirstPhysics.cxx++g");
  gROOT->LoadMacro("AliAnalysisHistosVertex.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskFPexample.cxx++g");

  AliAnalysisTaskFPexample* task = new AliAnalysisTaskFPexample(taskname);
  task->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
  task->SetCutTrackPt(0.15, 1e10);
  task->SetCutEta(0.8);
  task->SetCutVertexZ(10);
  mgr->AddTask(task);

  // set output root file name for different analysis
  TString outfilename;
  if (runtype == "grid") {
    outfilename = TString::Format("grid_%d_%s.root", gridRun, useRealData ? "data" : "sim");
    plugin->SetDefaultOutputs(kFALSE);
    plugin->SetOutputFiles(outfilename);
  } else {
    outfilename = TString::Format("%s_xx_%s.root", runtype, useRealData ? "data" : "sim");
  }

  // create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.Data());
  // connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
}
