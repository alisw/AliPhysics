TaskConsumer *AddTaskConsumer(const char *name, const char *prodname1, const char *prodname2)
{
// Provide as input the name of the consumer task and the name of the 
// producer task
	// pointer to the analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		::Error("AddTaskConsumer", "No analysis manager to connect to.");
		return NULL;
	}  
	// create the task
   TaskConsumer *task = new TaskConsumer(name);
   mgr->AddTask(task);

	// connecting the input/output containers
   TString outfile = mgr->GetCommonFileName();
   // input data feed
	AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
	mgr->ConnectInput (task, 0, cinput0 );
   // producer task
   TaskProducer *prod1 = mgr->GetTask(prodname1);
   TaskProducer *prod2 = mgr->GetTask(prodname2);
   if (!prod1 || !prod2) {
      ::Error("AddTaskConsumer", "Producer task %s or %s not found in the analysis manager", 
              prodname1, prodname2);
      return 0;
   }
   // Connect to exchange container
   AliAnalysisDataContainer *cinput1 = prod1->GetOutputSlot(2)->GetContainer();
   mgr->ConnectInput(task, 1, cinput1);
   AliAnalysisDataContainer *cinput2 = prod2->GetOutputSlot(2)->GetContainer();
   mgr->ConnectInput(task, 2, cinput2);
   
	AliAnalysisDataContainer *coutput1  = mgr->CreateContainer(
                TString::Format("output_%s", name),
                TList::Class(), AliAnalysisManager::kOutputContainer,
                TString::Format("%s:output",outfile.Data()));
	mgr->ConnectOutput(task, 1, coutput1);
   return task;
}
