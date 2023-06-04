///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskFlatenicityMCpred macro                     //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskFlatenicityMCpred *
AddTaskFlatenicityMCpred(const char *taskname = "Flat", bool isV0Malice = false, int nbins=600)

{
  // get the manager via the static access member. since it's static, you don't
  // need an instance of the class to call the function

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // now you create an instance of your task
  AliAnalysisTaskFlatenicityMCpred *taskFlat =
      new AliAnalysisTaskFlatenicityMCpred("taskFlat");
  if (!taskFlat)
    return 0x0;
  taskFlat->SetAnaModeV0Malice(isV0Malice);
taskFlat->SetNbinsNch(nbins);
  mgr->AddTask(taskFlat);

  mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
      taskFlat, 1,
      mgr->CreateContainer(
          Form("cList%s", taskname), TList::Class(),
          AliAnalysisManager::kOutputContainer,
          Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

  return taskFlat;
}
