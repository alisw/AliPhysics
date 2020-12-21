#include "AliMCSpectraWeights.h"
#if defined (__CLING__)
  #include "AliMCSpectraWeightsAnalysisTask.h"
#endif

AliMCSpectraWeightsAnalysisTask* AddTask_phuhn_MCSpectraWeights(const char* collisionSystem, const char* previousTrain)
{
  TString stCollisionSystem(collisionSystem);
  TString stTrainOutputPath(previousTrain);
  printf("Starting Analysis with AliMCSpectraWeights in %s collisions\n", stCollisionSystem.Data());
  if(stTrainOutputPath.Length()>5) //*.root
    printf("including results from train output %s\n", stTrainOutputPath.Data());

  AliMCSpectraWeights* fMCSpectraWeights = new AliMCSpectraWeights(stCollisionSystem.Data(), "fMCSpectraWeights", AliMCSpectraWeights::SysFlag::kNominal);
  if(stTrainOutputPath.Length()>5) fMCSpectraWeights->SetMCSpectraFile(stTrainOutputPath.Data());
  fMCSpectraWeights->Init();

  if(fMCSpectraWeights) printf("AliMCSpectraWeightsAnalysisTask:: obj created\n");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_testing_Unified", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  AliMCSpectraWeightsAnalysisTask *task = new AliMCSpectraWeightsAnalysisTask("AliMCSpectraWeightsAnalysisTask");
  task->SetUseMC(hasMC);
  if(type.Contains("ESD")) task->SetUseESD();
  else task->SetUseAOD();
  task->SetUseMultiplicity(kTRUE);
  task->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kMB ); //kINT7
  //For particle composition
  task->SetMCSpectraWeightObject(fMCSpectraWeights);
  // Debug info
  task->SetDebugLevel(1);
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPt_ParCor",
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      "AnalysisResults.root");

  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return task;

}
