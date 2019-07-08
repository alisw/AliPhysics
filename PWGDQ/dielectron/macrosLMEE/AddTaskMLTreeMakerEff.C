AliAnalysisTask *AddTaskMLTreeMakerEff( )
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
  gSystem->Exec(TString("alien_cp alien:///alice/cern.ch/user/s/selehner/cutlibs/LMEECutLib_slehner.C ."));
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

LMEECutLib* cutlib = new LMEECutLib();      
AliAnalysisTaskMLTreeMakerEff *task = new AliAnalysisTaskMLTreeMakerEff("treemaker");   

mgr->AddTask(task);

TString outputFN = AliAnalysisManager::GetCommonFileName();

AliAnalysisDataContainer *coutESD = mgr->CreateContainer(Form("output%d",1), TList::Class(),AliAnalysisManager::kOutputContainer,outputFN.Data());
mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(task, 1, coutESD);

return task;

}
