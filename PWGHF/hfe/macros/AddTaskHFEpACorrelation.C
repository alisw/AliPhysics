AliAnalysisTaskHFEpACorrelation *AddTaskHFEpACorrelation(
                                                        
                                                        Bool_t  isMC                    = kFALSE,
                                                        Int_t   triggerIndex    = 0,
                                                        Int_t   configIndex     = 0,
                                                        Int_t   centralityIndex = 0,
                                                        Bool_t  isAOD                   = kFALSE,
                                                        Bool_t isEMCal                  = kFALSE,
                                                        char * period                   = "b",
                                                        Int_t EMCalThreshould   = 0,
                                                        Bool_t ispp = kFALSE,
                                                        Int_t HadronPtCut = 0,
                                                        TString ElectronEfficiencyFile = "alien:///alice/cern.ch/user/h/hcorreia/Efficiency/Electron_Tracking.root",
                                                        TString HadronEfficiencyFile = "alien:///alice/cern.ch/user/h/hcorreia/Efficiency/Hadron_Tracking.root"
                                                        )
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
    if (!mgr) {
        ::Error("AddTaskEMCalHFEpA", "No analysis manager to connect to.");
        return NULL;
    }
    
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskEMCalHFEpA", "This task requires an input event handler");
        return NULL;
    }
    
    //_______________________
    //Config Task
    //gROOT->LoadMacro("ConfigEMCalHFEpACorrelation.C");
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pPb/ConfigHFEpACorrelation.C");
    AliAnalysisTaskHFEpACorrelation *task = ConfigHFEpACorrelation(isMC,triggerIndex,configIndex,centralityIndex,isAOD,isEMCal,EMCalThreshould,ispp,HadronPtCut);
    
    //_______________________
    //Trigger
    
    
    if(!isMC && (period=="d" || period=="e" || period=="f"))
    {
        if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
        if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
        if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
    }
    else if(!isMC)
    {
        if(triggerIndex == 0) task->SelectCollisionCandidates(AliVEvent::kINT7);
        if(triggerIndex == 1) task->SelectCollisionCandidates(AliVEvent::kEMC7);
        if(triggerIndex == 2) task->SelectCollisionCandidates(AliVEvent::kEMCEGA);
        //if(triggerIndex == 3) task->SelectCollisionCandidates(AliVEvent::kEMC8);
        //if(triggerIndex == 4) task->SelectCollisionCandidates(AliVEvent::kEMCEJE); //Jet Trigger
    }
    if(period == "pp")
    {
        task->SelectCollisionCandidates(AliVEvent::kMB);
    }
    
    //Seting online efficiency in TH3s
    TFile *fileE = TFile::Open(ElectronEfficiencyFile.Data());
    TH3F* ElecEffHisto = (TH3F*) fileE->Get(Form("Eff_Elec_Config%d_C%d",configIndex,centralityIndex));
    task->SetEfficiencyElectron(ElecEffHisto);

    
    TFile *fileH = TFile::Open(HadronEfficiencyFile.Data());
    TH3F* HadronEffHisto = (TH3F*) fileH->Get(Form("Eff_Hadron_Config%d_C%d",((Int_t)configIndex/1000)*1000,centralityIndex));
    task->SetEfficiencyHadron(HadronEffHisto);
    
    
    mgr->AddTask(task);
    
    TString containerName = mgr->GetCommonFileName();
    containerName += ":HFE_EMCal_pPb_elienos";
    containerName += Form("_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould);
    
    //Create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("chist_eh_%d_%d_%d_%d",triggerIndex,configIndex,centralityIndex,EMCalThreshould), TList::Class(),    AliAnalysisManager::kOutputContainer, containerName.Data());
    
    //Connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);
    
    return task;
}
