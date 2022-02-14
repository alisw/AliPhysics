AliAnalysisTaskSpectraTPCRun3* AddTaskSpectraTPCRun3(Bool_t mc = kFALSE, Bool_t aod = kTRUE, Bool_t o2pid = kTRUE)
{

  /* check analysis manager */
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTOFTree", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskTOFTree", "cannot get input event handler");
    return NULL;
  }

  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (aod && str.CompareTo("AOD")) {
    Error("AddTaskSpectraTPCRun3", "input data type is not \"AOD\"");
    return NULL;
  }
  if (!aod && str.CompareTo("ESD")) {
    Error("AddTaskSpectraTPCRun3", "input data type is not \"ESD\"");
    return NULL;
  }

  if (mc) {
    /* check MC truth event handler */
    if (!mgr->GetMCtruthEventHandler()) {
      Error("AddTaskSpectraTPCRun3", "cannot get MC truth event handler");
      return NULL;
    }
  }

  /* get common input data container */
  AliAnalysisDataContainer* inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddTaskTOFTree", "cannot get common input container");
    return NULL;
  }

  /* create output data container */
  AliAnalysisDataContainer* outputc = mgr->CreateContainer("SpectraTPCRun3", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  if (!outputc) {
    Error("", "cannot create output container \"Histos\"");
    return NULL;
  }

  /*  create task and connect input/output */
  AliAnalysisTaskSpectraTPCRun3* task = new AliAnalysisTaskSpectraTPCRun3("SpectraTPCRun3", aod, o2pid);
  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, outputc);
  if (mc)
    task->fMCmode = kTRUE;

  /* return task */
  return task;
}
