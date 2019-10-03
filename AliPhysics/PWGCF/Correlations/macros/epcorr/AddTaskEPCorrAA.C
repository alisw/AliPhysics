AliAnalysisTaskEPCorrAA* AddTaskEPCorrAA( const char* outputFileName = 0, const char *containerName = "EPCorrAA_TPCOnly_151202", const char* folderName = "PWGCF_EPCorrelation" )
{


    // Get a pointer to the analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        cout<<"No analysis manager found."<<endl;
        return 0x0;
    }


	// Additional tasks will be managed by 'Dependendies' in wagon
/*
    // Centrality
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 
//  taskCentrality->SetPass(2);
    taskCentrality->SelectCollisionCandidates( AliVEvent::kAny );

    // PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *taskPID = AddTaskPIDResponse(0); // 0 (real data), 1 (MC)
    if (!taskPID) { Printf("no PID task"); return; };

    // Physics Selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); // 150608
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(0, kTRUE); // 0 (real), 1 (MC gen or rec.)
//  physSelTask->SelectCollisionCandidates(AliVEvent::kMB);

	// V0 Event Plane
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
    AliVZEROEPSelectionTask *VZEROEPTask = AddTaskVZEROEPSelection();
*/
	// My Task
    AliAnalysisTaskEPCorrAA* EPCorrAATask = new AliAnalysisTaskEPCorrAA(containerName);
//    AliAnalysisTaskSE* EPCorrAATask = new AliAnalysisTaskEPCorrAA(containerName);
//    EPCorrAATask->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented

    // Add the task.
    mgr->AddTask(EPCorrAATask);

    // Data containers.
    AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
    mgr->ConnectInput(EPCorrAATask, 0, cinput);

    if (!outputFileName) {outputFileName = AliAnalysisManager::GetCommonFileName();}

    AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(containerName, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFileName, folderName));

    mgr->ConnectOutput(EPCorrAATask,1,coutput1);

    return EPCorrAATask;

}
