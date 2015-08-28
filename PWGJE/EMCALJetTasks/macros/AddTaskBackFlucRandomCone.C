AliAnalysisTaskBackFlucRandomCone *AddTaskBackFlucRandomCone(TString ntracks = "PicoTracks"){
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr){
      Error("AddTaskBackFlucRandomCone","No analysis manager found.");
      return 0;
   }
   if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskBackFlucRandomCone", "This task requires an input event handler");
      return NULL;
    }
    TString wagonName = Form("BackFlucRandomCone");
    AliAnalysisTaskBackFlucRandomCone *task = new AliAnalysisTaskBackFlucRandomCone(wagonName);
    AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
    
    //Connnect input
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    //Connect output
    TString contName(wagonName);
    TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,1,coutput1);

    return task;
}
