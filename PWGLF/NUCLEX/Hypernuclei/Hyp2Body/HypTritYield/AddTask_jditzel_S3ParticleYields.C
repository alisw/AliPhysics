AliAnalysisTask *AddTask_jditzel_S3ParticleYields() {
  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  TString fileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisTaskS3ParticleYields* task = new AliAnalysisTaskS3ParticleYields("jditzelTaskS3ParticleYields");   
  if(!task) return 0x0;

  task->SelectPIDcheckOnly(kFALSE);
  task->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);
  task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);

    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histogramsheLp", TList::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("treep", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("treehe", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("treeL", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("treeGenheLp", TTree::Class(),AliAnalysisManager::kOutputContainer,mgr->GetCommonFileName());
	mgr->ConnectOutput(task, 1, coutput1);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);
    mgr->ConnectOutput(task, 4, coutput4);
    mgr->ConnectOutput(task, 5, coutput5);
	return task;
}
