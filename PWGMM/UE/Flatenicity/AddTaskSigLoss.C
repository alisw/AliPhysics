
class AliAnalysisDataContainer;

AliAnalysisTaskSigLoss* AddTaskSigLoss(const Char_t *taskname = "SigLoss",
                                       bool doMultSel  = kFALSE,
                                       float lowMultCut = -1,
                                       float highMultCut  = -1,
                                       const char* multEst = "V0M")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSigLoss", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSigLoss", "This task requires an input event handler");
    return nullptr;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type.Contains("AOD")) {
    ::Error("AddTaskSigLoss", "This task requires to run on ESD");
    return nullptr;
  }

  TString kContSuffix("");  
  // Create and configure the task
  AliAnalysisTaskSigLoss* userTask = new AliAnalysisTaskSigLoss("taskSig");
  if (!userTask)
    return 0x0;
  
  if (doMultSel) {
    userTask->SetDoMultSel(doMultSel);
    userTask->SetMultEst(multEst);
    userTask->SetMultiplicityRange(lowMultCut, highMultCut); 
    kContSuffix.Append(Form("_%s-%.1fto%.1f", multEst, (lowMultCut < 0) ?   0 : lowMultCut, (highMultCut  < 0) ? 100 :  highMultCut));
  }
  
  userTask->SetupStandardEventCutsForRun2();
  mgr->AddTask(userTask);

  mgr->ConnectInput(userTask,  0, mgr->GetCommonInputContainer());
   
  TString kContainerName("cListSigLoss");
  kContainerName += kContSuffix;
  mgr->ConnectOutput(userTask, 1, mgr->CreateContainer(kContainerName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

  return userTask;
}
