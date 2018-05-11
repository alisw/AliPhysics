AlidNdPtUnifiedAnalysisTask* AddTask_yelmardb_Unified()
{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    if (!mgr) {
        Error("AddTask_yelmardb_Unified", "No analysis manager found.");
        return 0;
    }

    // Switch off all AliInfo (too much output!!!)
    AliLog::SetGlobalLogLevel(AliLog::kError);
    mgr->SetDebugLevel(0);


    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

    //AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

    AlidNdPtUnifiedAnalysisTask *task = new AlidNdPtUnifiedAnalysisTask("AlidNdPtUnifiedAnalysisTask_yelmardb");
    task->SetUseMC(hasMC);
    if(type.Contains("ESD")) task->SetUseESD();
    else task->SetUseAOD();
    task->SetUseMultiplicity(kTRUE);
    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetTriggerMask(AliVEvent::kINT7 ); //kINT7


    const Int_t multNbins =200;
    Double_t binsMult[multNbins+1];
    for (Int_t ii = 0; ii<=multNbins;ii++){binsMult[ii]=5*ii-0.5;}
    task->SetBinsMultCent(multNbins,binsMult);

    //TFile *trkEffFile = TFile::Open("$UNIFIED/macros/trkEffParametrisation.root");
    //if(!trkEffFile){printf("No trk eff file\n");}
    //else TF1* trkEffParam = (TF1*)trkEffFile->Get("trkEffFun");
    // trkEffParam->SetDirectory(0);
    // trkEffFile->Close();

    // task->SetTrkEffParametrisation(trkEffParam);

    /// Acceptance cuts for tracks

    task->SetMinEta(-1.);
    task->SetMaxEta(1.);
    task->SetMinPt(0.10);

    task->Set2013pA(kFALSE);
    task->Set2015data(kTRUE);

    ///TOF pileup, kTRUE only for Matching efficiency studies
    task->SetTOFbunchCrossing(kFALSE);

    task->SetMeanXYZv(0.0,0.0,0.0);
    task->SetSigmaMeanXYZv(1.0,1.0,10.0);
    task->SetZvtx(30.);
    task->SetEventTriggerRequired(kTRUE);

    // Quality cuts for tracks
    task->SetTPCRefit(kTRUE);
    task->SetITSRefit(kTRUE);
    task->SetKinkDaughters(kFALSE);
    task->SetRatioCrossedRowsOverFindableClustersTPC(0.8);
    task->SetFractionSharedClustersTPC(0.4);
    task->SetMaxchi2perTPCclu(4.);
    task->SetClusterReqITS(kTRUE);
    task->SetMaxchi2perITSclu(36.);
    task->SetDCAtoVertex2D(kFALSE);
    task->SetSigmaToVertex(kFALSE);
    task->SetDCAtoVertexZ(2.0);
    task->SetDCAtoVertexXYPtDep("0.0182+0.0350/pt^1.01");
    task->SetMaxChi2TPCConstrained(36.);
    task->SetMinLenghtInActiveZoneTPC(0);
    task->SetGeometricalCut(kFALSE,2,130,1.5,0.85,0.7); ///if kTRUE comment CrossedRowsTPC cut
    task->SetMinCrossedRowsTPC(120);

    mgr->AddTask(task);

    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPt",
            TList::Class(),
            AliAnalysisManager::kOutputContainer,
            "AnalysisResults.root");



    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);

    return task;

}
