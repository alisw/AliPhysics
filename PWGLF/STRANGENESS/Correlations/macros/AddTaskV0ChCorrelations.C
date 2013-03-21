// AddTaskV0ChCorrelations.C
//
// AddTask for AliAnalysisTaskV0ChCorrelations task
//
AliAnalysisTaskV0ChCorrelations *AddTaskV0ChCorrelations(const Bool_t bMCtruth=kTRUE,
														 Float_t DcaDToPV = 0.1,
														 Float_t DcaV0D = 1.0
														)
{
  // Creates a V0Ch correlations analysis task and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Error("AddTaskV0ChCorrelations", "No analysis manager to connect to.");
      return NULL;
    }

	// create task
    AliAnalysisTaskV0ChCorrelations* task = new AliAnalysisTaskV0ChCorrelations("V0ChCorrelations_task");
    task->SetAnalysisMC(bMCtruth);
	task->SetDcaDToPV(DcaDToPV);
	task->SetDcaV0D(DcaV0D);
    mgr->AddTask(task);
    
    // Create ONLY the output containers for the data produced by the task.
    // Get and connect other common input/output containers via the manager as below
    //==============================================================================
    TString outputFileName = AliAnalysisManager::GetCommonFileName();

    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
        
    // connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);
        
    // Return task pointer at the end
    return task;
}

