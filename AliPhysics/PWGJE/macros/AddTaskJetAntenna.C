

AliAnalysisTaskJetAntenna* AddTaskJetAntenna(const char* bRec1,const char* bRec2, UInt_t filterMask = 272 , Float_t ptTrackMin = 0.15, Int_t kTriggerMask=0, Float_t fCutTM=0.15){

   Printf("adding task jet antenna\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetAntenna", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetAntenna", "This task requires an input event handler.");
      return NULL;
   }

     
  

  TString typeRec(bRec1);
  TString typeGen(bRec2);
      
  AliAnalysisTaskJetAntenna *task = new AliAnalysisTaskJetAntenna(Form("JetAntenna_%s_%s_%d_%f",bRec1,bRec2,kTriggerMask,fCutTM));
   


   task->SetBranchNames(bRec1,bRec2);
   task->SetOfflineTrgMask(kTriggerMask);
   task->SetCentMin(0.);
   task->SetCentMax(100.);
   task->SetFilterMask(filterMask); 
   task->SetTMCut(fCutTM);  
  
  

 

   mgr->AddTask(task);


   AliAnalysisDataContainer *coutputJetAntenna = mgr->CreateContainer(Form("pwgjejetantenna_%s_%s_%d_%f",bRec1,bRec2,kTriggerMask,fCutTM), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_jetantenna_%s_%s_%d_%f",AliAnalysisManager::GetCommonFileName(),bRec1,bRec2,kTriggerMask,fCutTM));





   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetAntenna);

   return task;
}
