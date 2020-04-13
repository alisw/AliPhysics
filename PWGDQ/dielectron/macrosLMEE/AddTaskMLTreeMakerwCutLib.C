AliAnalysisTask *AddTaskMLTreeMakerwCutLib(   
                                                    Int_t trackCut=1,
                                                    Int_t PIDCut=1,
                                                    Int_t evCut=1,
                                                    Double_t centmin=0.,
                                                    Double_t centmax=100.,
                                                    Bool_t SetPIDCorrection=kFALSE,
                                                    Bool_t useAODFilterCuts=kFALSE,
                                                    Bool_t isMC,
                                                    TString TMVAweight,
                                                    Int_t wagonnr=0,
                                                    Bool_t usePileRej=kTRUE,
                                                    Bool_t isUPC=kFALSE
                                                    )
{

//get the current analysis manager
AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
if(!mgr){

    Error("AddTask_slehner_TreeMakerWCutLib", "No analysis manager found.");
    return 0;
}

//Base Directory for GRID / LEGO Train  
TString configBasePath= "./";

TString configLMEECutLib("LMEECutLib_slehner.C");
TString configLMEECutLibPath(configBasePath+configLMEECutLib);
//  gSystem->Exec(TString("alien_cp alien:///alice/cern.ch/user/s/selehner/cutlibs/LMEECutLib_slehner.C ."));
//LOAD CUTLIB
if(gSystem->Exec(Form("ls %s", configLMEECutLibPath.Data()))==0){

    ::Info("AddTask_slehner_TreeMakerWCutLib","loading LMEECutLib: %s",configLMEECutLibPath.Data());
    gROOT->LoadMacro(configLMEECutLibPath.Data());
} 
else{
    ::Info("AddTask_slehner_TreeMakerWCutLib","LMEECutLib not found: %s", configLMEECutLibPath.Data());
    return 0; // if return is not called, the job will fail instead of running wihout this task... (good for local tests, bad for train)
}

//Do we have an MC handler?
Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != 0x0);
::Info("AddTask_slehner_TreeMakerWCutLib","hasMC = %d",hasMC);

TString weightFile=TMVAweight;
::Info("AddTask_slehner_TreeMakerWCutLib","TMVA weights file: alien:///alice/cern.ch/user/s/selehner/TMVAweights/%s",weightFile.Data());

LMEECutLib* cutlib = new LMEECutLib();    
AliAnalysisTaskMLTreeMaker *task = new AliAnalysisTaskMLTreeMaker("treemaker",weightFile);   

if(SetPIDCorrection){
  task->SetUseCorr(kTRUE);
//  for(int det=1; det<4; det++){
  for(int det=1; det<3; det++){
    TH3D mean = cutlib->SetEtaCorrection(det, isMC, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, 1);
    TH3D width = cutlib->SetEtaCorrection(det, isMC, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, 2);
    task->SetCorrWidthMean(det,(TH3D*)width.Clone(),(TH3D*)mean.Clone());
  }
}
else  task->SetUseCorr(kFALSE);

task->isMC(isMC);
if(isUPC)task->SelectCollisionCandidates(AliVEvent::kAny);
else task->SelectCollisionCandidates(AliVEvent::kINT7);
task->SetCentralityPercentileRange(centmin,centmax);
task->SetupTrackCuts(cutlib->GetTrackCuts(trackCut,PIDCut,0,useAODFilterCuts,TMVAweight,isUPC));
task->SetUseTMVA(kTRUE);
task->SetupEventCuts(cutlib->GetEventCuts(centmin, centmax,usePileRej,isUPC));

mgr->AddTask(task);

TString outputFN = AliAnalysisManager::GetCommonFileName();

AliAnalysisDataContainer *coutESD = mgr->CreateContainer(Form("output%d",wagonnr), TList::Class(),AliAnalysisManager::kOutputContainer,outputFN.Data());
mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(task, 1, coutESD);

return task;

}
