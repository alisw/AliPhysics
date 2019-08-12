AliAnalysisTaskHadronPhiCorr *AddTaskHadronPhiCorr(Bool_t isHH = kFALSE, Float_t multLow = 0.0, Float_t multHigh = 100.0, const char* suffix=""){
    //get the current analysis manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskHFE", "No analysis manager found.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskHFE", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
        
    stringstream taskStream;
    taskStream << "hPhiCorr_mult_" << multLow << "_" << multHigh << suffix;

    TString taskName = taskStream.str();

    AliAnalysisTaskHadronPhiCorr *hPhiCorr = new AliAnalysisTaskHadronPhiCorr(taskName.Data(), isHH, multLow, multHigh); 
    hPhiCorr->SelectCollisionCandidates(AliVEvent::kINT7);
    TString filename = "AnalysisResults.root";

    stringstream containerStream;
    if(isHH){
        containerStream << "hhCorr_mult_" << multLow << "_" << multHigh;
    }else{
        containerStream << "phiCorr_mult_" << multLow << "_" << multHigh;
    }
    containerStream << suffix;
    TString containerName = containerStream.str();

    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containerName, TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data());
    mgr->ConnectInput(hPhiCorr, 0, cinput);
    mgr->ConnectOutput(hPhiCorr, 1, coutput1);
    mgr->AddTask(hPhiCorr);

    return hPhiCorr;
}
