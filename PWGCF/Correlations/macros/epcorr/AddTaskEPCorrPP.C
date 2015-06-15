AliAnalysisTaskEPCorrPP* AddTaskEPCorrPP( const char* outputFileName = 0, const char *containerName = "EPCorrPP", const char* folderName = "PWGCF_EPCorr" )
{

    // Get a pointer to the analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        cout<<"AddTaskDiHadronPID.C -> No analysis manager found."<<endl;
        return 0x0;
    }


	// Additional tasks

/*  // no need in PP
    // Centrality
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 
//  taskCentrality->SetPass(2);
    taskCentrality->SelectCollisionCandidates( AliVEvent::kAny );
*/

    // PID
    gROOT->LoadMacro("$ALICE_ROOT/../src/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *taskPID = AddTaskPIDResponse(0); // 0 (real data), 1 (MC)
    if (!taskPID) { Printf("no phys sel task"); return; };

    // Physics Selection
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); // 150608
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel, kTRUE);
//  physSelTask->SelectCollisionCandidates(AliVEvent::kMB);


	// My Task
    AliAnalysisTaskEPCorrPP* EPCorrPPTask = new AliAnalysisTaskEPCorrPP(containerName);
//    AliAnalysisTaskSE* EPCorrPPTask = new AliAnalysisTaskEPCorrPP(containerName);
    EPCorrPPTask->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented

    // Add the task.
    mgr->AddTask(EPCorrPPTask);

    // Data containers.
    AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
    mgr->ConnectInput(EPCorrPPTask, 0, cinput);

    if (!outputFileName) {outputFileName = AliAnalysisManager::GetCommonFileName();}

    AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(containerName, TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFileName, folderName));

    mgr->ConnectOutput(EPCorrPPTask,1,coutput1);

    return EPCorrPPTask;



}
