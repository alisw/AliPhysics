TString fileNameBase="AnalysisResults.root";
AliAnalysisTask*  AddAliTwoParticlePIDCorrTask()
{
  TString taskname = "pPbPIDCorr";
  taskname.Append("taskssssssssss");
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPIDCorr", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPIDCorr", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType();

   AliTwoParticlePIDCorr *taskpPbPIDCorr = new AliTwoParticlePIDCorr("TwoParticlePIDCorrTask");
   
   //Trigger - Physics Selection
   taskpPbPIDCorr->SelectCollisionCandidates(AliVEvent::kINT7);
   mgr->AddTask(taskpPbPIDCorr);
   
   AliAnalysisDataContainer *coutFA = mgr->CreateContainer(taskname.Data(), 
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,fileNameBase.Data());

  mgr->ConnectInput(taskpPbPIDCorr, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskpPbPIDCorr, 1, coutFA);
  
  return taskpPbPIDCorr;
}
