AliAnalysisTask *AddTaskTrackingUncertHMTF(Bool_t readMC = kFALSE,
                                          TString trigClass = "CINT1B",
                                          AliVEvent::EOfflineTriggerTypes trigMask = AliVEvent::kMB,
                                          AliAnalysisTrackingUncertaintiesHMTF::ESpecies_t specie=(AliAnalysisTrackingUncertaintiesHMTF::kSpecPion|AliAnalysisTrackingUncertaintiesHMTF::kSpecKaon),
                                          Double_t MaxDCAxy = 2.4,
                                          Double_t MaxDCAz  = 3.2,
                                          Double_t MaxEta   = 0.8,
                                          Double_t CrossRowsOverFndCltTPC = 0.8) {
    
    //
    // add task of tracking uncertainty
    //
    //
    //get the current analysis manager
    //
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskTrackingUncertHMTF", "No analysis manager found.");
        return 0;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskImpParDistrib", "This task requires an input event handler");
        return NULL;
    }
    
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("AOD")){
        ::Error("AddTaskImpParDistrib", "This task requires to run on ESD");
        return NULL;
    }
    
    //
    //========= Add task for standard analysis to the ANALYSIS manager ====
    //
    AliAnalysisTrackingUncertaintiesHMTF *task    = new AliAnalysisTrackingUncertaintiesHMTF("trackingUncertainty");
    //
    task->SetReadMC(readMC);
    task->SetTriggerClass(trigClass.Data());
    task->SetTriggerMask(trigMask);
    task->SetSpecie(specie);
    task->SetMaxDCAxy(MaxDCAxy);
    task->SetMaxDCAz(MaxDCAz);
    task->SetEtaRange(MaxEta);
    task->SetCrossRowsOverFndCltTPC(CrossRowsOverFndCltTPC);
    
    mgr->AddTask(task);
    ULong64_t SPeciee = task->GetSpecie();
    TString suffix = "_";
    if(SPeciee&AliAnalysisTrackingUncertaintiesHMTF::kSpecElectron) suffix += "e";
    if(SPeciee&AliAnalysisTrackingUncertaintiesHMTF::kSpecPion)     suffix += "Pi";
    if(SPeciee&AliAnalysisTrackingUncertaintiesHMTF::kSpecKaon)     suffix += "K";
    if(SPeciee&AliAnalysisTrackingUncertaintiesHMTF::kSpecProton)   suffix += "p";
    if(SPeciee&AliAnalysisTrackingUncertaintiesHMTF::kAll)          suffix  = "_All";
    //
    //
    //======================================================================
    //              data containers
    //======================================================================
    //            find input container
    //below the trunk version
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    
    //define output containers
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("trackingUncert%s",suffix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    
    //connect containers
    mgr->ConnectInput  (task, 0, cinput );
    mgr->ConnectOutput (task, 1, coutput1);
    //
    //
    //
    return task;
    
}

