AliAnalysisTaskPIDCORR *AddTaskPIDCORR()
{
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		
		return NULL;
	} 
AliAnalysisTaskPIDCORR* task = new AliAnalysisTaskPIDCORR("PIDCORR_pp7TeV");

mgr->AddTask(task);

 TString outputFileName = AliAnalysisManager::GetCommonFileName();
 AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
 AliAnalysisDataContainer *coutput = mgr->CreateContainer("PIDCORR_pp7TeV", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
        

outputFileName.Data()


    // connect input/output
mgr->ConnectInput(task, 0, cinput);
mgr->ConnectOutput(task, 1, coutput);

return task;

}

