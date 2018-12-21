AliAnalysisTaskSEHFSystPID *AddTaskHFSystPID(int system = 0,
                                            bool readMC = false,
                                            TString trigClass = "",
                                            AliVEvent::EOfflineTriggerTypes trigMask = AliVEvent::kINT7,
                                            TString outputSuffix = "_ppMB_kINT7",
                                            float nsigmafortag = 0.02,
                                            double fracdownsampl = 1.,
                                            double ptmaxdownsampl = 0.) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskSEHFSystPID", "No analysis manager found.");
        return 0;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AliAnalysisTaskSEHFSystPID", "This task requires an input event handler");
        return NULL;
    }

    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("ESD")){
        ::Error("AliAnalysisTaskSEHFSystPID", "This task requires to run on AOD");
        return NULL;
    }

    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetMinNClustersTPC(50);
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetMaxDCAToVertexZ(2.5);
    esdTrackCuts->SetPtRange(0.,1.e10);
  
    //========= Add task for standard analysis to the ANALYSIS manager ====
    AliAnalysisTaskSEHFSystPID *task = new AliAnalysisTaskSEHFSystPID("SystNsigmaPIDtask",system);
    task->SetReadMC(readMC);
    task->SetTriggerInfo(trigClass.Data(),trigMask);
    task->SetESDtrackCuts(esdTrackCuts);
    task->SetNsigmaKaonForTagging(nsigmafortag);
    if(fracdownsampl<1.)task->EnableDownSampling(fracdownsampl,ptmaxdownsampl);
    mgr->AddTask(task);

    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWGHF_D2H_SystNsigmaPID";

    //define input container
    AliAnalysisDataContainer *cinput = mgr->CreateContainer("cinputPID",TChain::Class(),AliAnalysisManager::kInputContainer);
    //define output containers
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("coutputPIDhistos%s",outputSuffix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("coutputPIDtree%s",outputSuffix.Data()),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

    //connect containers
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput (task, 1, coutput1);
    mgr->ConnectOutput (task, 2, coutput2);
    return task;
}
