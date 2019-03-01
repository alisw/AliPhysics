AliAnalysisTaskMKTest* AddTaskMKTest(TString name = "MKTest")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { return 0; }
    if (!mgr->GetInputEventHandler()) { return 0; }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    //fileName += ":TaskCentBiasMK";      // create a subfolder in the file
    fileName = TString("out_");
    fileName += name;
    fileName += ".root";
    // now we create an instance of your task
    AliAnalysisTaskMKTest* task = new AliAnalysisTaskMKTest(name.Data());   
    if(!task) return 0x0;
    //configure the task  
    /*
    // cuts0 
    AliESDtrackCuts* cuts0 = new AliESDtrackCuts("TPC-ITS with sector edge");
    cuts0->SetRequireTPCRefit(kTRUE);
    //cuts0->SetMinNCrossedRowsTPC(120); 
    cuts0->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    cuts0->SetMaxChi2PerClusterTPC(4);
    cuts0->SetMaxFractionSharedTPCClusters(0.4);
    //cuts0->SetMaxDCAToVertexXY(3.0); remove this cut!
    cuts0->SetRequireITSRefit(kTRUE);
    cuts0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    cuts0->SetMaxChi2PerClusterITS(36.);
    cuts0->SetDCAToVertex2D(kFALSE);
    cuts0->SetRequireSigmaToVertex(kFALSE);
    cuts0->SetMaxDCAToVertexZ(2.0);
    // 7*(0.0026+0.0050/pt^1.01)
    cuts0->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    cuts0->SetAcceptKinkDaughters(kFALSE);
    // tpcc cut
    cuts0->SetMaxChi2TPCConstrainedGlobal(36.);
    // Geometrical-Length Cut
    cuts0->SetCutGeoNcrNcl(3,130,1.5,0.85,0.7);
    cuts0->SetEtaRange(-0.8,0.8);
    
    task->SetESDtrackCuts(cuts0);
     */
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer(name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
