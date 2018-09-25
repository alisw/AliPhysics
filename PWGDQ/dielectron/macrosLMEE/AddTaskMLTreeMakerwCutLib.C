AliAnalysisTask *AddTaskMLTreeMakerwCutLib(   
                                                    Int_t trackCut=1,
                                                    Int_t PIDCut=1,
                                                    Int_t evCut=1,
                                                    Double_t centmin=0.,
                                                    Double_t centmax=100.,
                                                    Bool_t SetTPCCorrection=kFALSE,
                                                    Bool_t useAODFilterCuts=kFALSE
                                                    )
{

//get the current analysis manager
AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
if(!mgr){

    Error("AddTask_slehner_TreeMakeWCutLib", "No analysis manager found.");
    return 0;
}



//Base Directory for GRID / LEGO Train  
TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";

TString configLMEECutLib("LMEECutLib_slehner.C");
TString configLMEECutLibPath(configBasePath+configLMEECutLib);

//LOAD CUTLIB
if(gSystem->Exec(Form("ls %s", configLMEECutLibPath.Data()))==0){

    ::Info("AddTask_slehner_TreeMakeWCutLib","loading LMEECutLib: %s",configLMEECutLibPath.Data());
    gROOT->LoadMacro(configLMEECutLibPath.Data());
} 
else{
    ::Info("AddTask_slehner_TreeMakeWCutLib","LMEECutLib not found: %s", configLMEECutLibPath.Data());
    return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
}

//Do we have an MC handler?
Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
::Info("AddTask_slehner_TreeMakeWCutLib","hasMC = %d",hasMC);


LMEECutLib* cutlib = new LMEECutLib();      
AliAnalysisTaskMLTreeMaker *task = new AliAnalysisTaskMLTreeMaker("treemaker");   

if(SetTPCCorrection){
    TH3D mean = cutlib->SetEtaCorrectionTPC(AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE,1);
    TH3D width = cutlib->SetEtaCorrectionTPC(AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE,2);
    task->SetUseCorr(kTRUE);
    task->SetCorrWidthMean((TH3D*)width.Clone(),(TH3D*)mean.Clone());
}
else  task->SetUseCorr(kFALSE);

task->SelectCollisionCandidates(AliVEvent::kINT7);
task->SetupTrackCuts(cutlib->GetTrackCuts(trackCut,PIDCut,useAODFilterCuts));
task->SetupEventCuts(cutlib->GetEventCuts());
task->SetCentralityPercentileRange(centmin,centmax);

mgr->AddTask(task);

TString outputFN = AliAnalysisManager::GetCommonFileName();

AliAnalysisDataContainer *coutESD = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputFN.Data());
mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(task, 1, coutESD);

return task;

}
