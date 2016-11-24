AliAnalysisTask* AddTaskAnalysisOnlineQA()
{
  
    bool useSparse = 0;
    
    //get the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  

  // check if there is some input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }

  //attach the actual task
    AliPerformanceTask* task = new AliPerformanceTask("PerformanceQA");
    mgr->AddTask(task);

    AliRecInfoCuts *pRecInfoCutsTPC = new AliRecInfoCuts("pRecInfoCutsTPC");
    pRecInfoCutsTPC->SetMaxDCAToVertexXY(2.4);
    pRecInfoCutsTPC->SetMaxDCAToVertexZ(3.0);
    pRecInfoCutsTPC->SetRequireSigmaToVertex(kFALSE);
    pRecInfoCutsTPC->SetRequireTPCRefit(kFALSE);
    pRecInfoCutsTPC->SetAcceptKinkDaughters(kTRUE);
    pRecInfoCutsTPC->SetMinNClustersTPC(70);
    pRecInfoCutsTPC->SetMaxChi2PerClusterTPC(1000000.);
    pRecInfoCutsTPC->SetDCAToVertex2D(kFALSE);
    pRecInfoCutsTPC->SetHistogramsOn(kFALSE);

    AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts("pMCInfoCuts");
    if(pMCInfoCuts) {
        pMCInfoCuts->SetMinTrackLength(70);
    }

    AliPerformanceTPC *pCompTPC0 = new AliPerformanceTPC("AliPerformanceTPC","AliPerformanceTPC",0,kFALSE,-1,0,useSparse);
    pCompTPC0->SetAliMCInfoCuts(pMCInfoCuts);
    pCompTPC0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompTPC0->SetUseTrackVertex(kTRUE);
    task->AddPerformanceObject(pCompTPC0);

    AliPerformanceMatch *pCompMatch0 = new AliPerformanceMatch("AliPerformanceMatchTPCITS","AliPerformanceMatchTPCITS",0,kFALSE,useSparse);
    if(!pCompMatch0) {
        Error("AddTaskPerformanceMatch", "Cannot create AliPerformanceMatchTPCITS");
    }
    pCompMatch0->SetAliRecInfoCuts(pRecInfoCutsTPC);
    pCompMatch0->SetAliMCInfoCuts(pMCInfoCuts);
    task->AddPerformanceObject(pCompMatch0);

    AliPerformanceMatch *pCompConstrain6 = new AliPerformanceMatch("AliPerformanceMatchTPCConstrain","AliPerformanceMatchTPCConstrain",2,kFALSE,useSparse);
    if(!pCompConstrain6) {
        Error("AddTaskPerformanceTPCdEdxQA", "Cannot create AliPerformanceMatchTPCConstrain");  }
    pCompConstrain6->SetAliRecInfoCuts(pRecInfoCutsTPC);
    task->AddPerformanceObject(pCompConstrain6);
    
    delete pRecInfoCutsTPC;
    delete pMCInfoCuts;
    
    
    AliAnalysisDataContainer *input = mgr->GetCommonInputContainer(); //this is standard
    AliAnalysisDataContainer *output1 = mgr->CreateContainer("TPCQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TPC_%s", mgr->GetCommonFileName(), task->GetName()));
    
    mgr->ConnectInput(task,0,input);
    mgr->ConnectOutput(task,1,output1);
    
  return task;
}
