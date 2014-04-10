TaskProducer *AddTaskProducer(const char *name)
{
	// pointer to the analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskProducer", "No analysis manager to connect to.");
		return NULL;
	}  
	// create the task
   TaskProducer *task = new TaskProducer(name);
   mgr->AddTask(task);

	// connecting the input/output containers
   TString outfile = mgr->GetCommonFileName();
   // input data feed
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	mgr->ConnectInput (task, 0, cinput0 );
   // producer output
	AliAnalysisDataContainer *coutput1  = mgr->CreateContainer(
                TString::Format("output_%s", name),
                TList::Class(), AliAnalysisManager::kOutputContainer,
                TString::Format("%s:output",outfile.Data()));
	mgr->ConnectOutput(task, 1, coutput1);
   // exchange output
	AliAnalysisDataContainer *coutput2  = mgr->CreateContainer(
                TString::Format("exchange_%s", name),
                TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
   mgr->ConnectOutput(task, 2, coutput2);
   return task;
}
     
   
