AliAnalysisTaskParticleEfficiency *
AddAnalysisTaskParticleEfficiency(const Char_t *partName)
{

  /* check analysis manager */
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddAnalysisTaskParticleEfficiency", "cannot get analysis manager");
    return NULL;
  }

  /* check input event handler */
  if (!mgr->GetInputEventHandler()) {
    Error("AddAnalysisTaskParticleEfficiency", "cannot get input event handler");
    return NULL;
  }
  
  /* check input data type */
  TString str = mgr->GetInputEventHandler()->GetDataType();
  if (str.CompareTo("ESD")) {
    Error("AddAnalysisTaskParticleEfficiency", "input data type is not \"ESD\"");
    return NULL;
  }

  /* check MC truth event handler */
  if (!mgr->GetMCtruthEventHandler()) {
    Error("AddAnalysisTaskParticleEfficiency", "cannot get MC truth event handler");
    return NULL;
  }
  
  /* get common input data container */
  AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
  if (!inputc) {
    Error("AddAnalysisTaskParticleEfficiency", "cannot get common input container");
    return NULL;
  }
  
/* create output data container */
  AliAnalysisDataContainer *outputc1 = mgr->CreateContainer(partName, TList::Class(), AliAnalysisManager::kOutputContainer, "ParticleEfficiency.root");
  if (!outputc1) {
    Error("", "cannot create output container \"Histos\"");
    return NULL;
  }

    /*  create task and connect input/output */
  AliAnalysisTaskParticleEfficiency *task = new AliAnalysisTaskParticleEfficiency(partName);
  mgr->ConnectInput(task, 0, inputc);
  mgr->ConnectOutput(task, 1, outputc1);

  /* return task */
  return task;
  
}
