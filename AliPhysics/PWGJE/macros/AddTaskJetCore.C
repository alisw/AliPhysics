

AliAnalysisTaskJetCore* AddTaskJetCore(const char* bRec1,const char* bRec2, UInt_t filterMask = 272 , Float_t ptTrackMin = 0.15, Int_t kTriggerMask=0, Int_t eventClassMin = 0, Int_t eventClassMax = 4,Int_t kHardest=1,Float_t kTTminr=11,Float_t kTTmaxr=13,Float_t kTTmins=15,Float_t kTTmaxs=19,Int_t kPhiBkg=0){

   Printf("adding task jet response\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetCore", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetCore", "This task requires an input event handler.");
      return NULL;
   }

     
  

  TString typeRec(bRec1);
  TString typeGen(bRec2);
      
  AliAnalysisTaskJetCore *task = new AliAnalysisTaskJetCore(Form("JetCore_%s_%s_%d_%d_%f%f%f%f_%f",bRec1,bRec2,kTriggerMask,kHardest,kTTminr,kTTmaxr,kTTmins,kTTmaxs,kPhiBkg));
   


   task->SetBranchNames(bRec1,bRec2);
   task->SetOfflineTrgMask(kTriggerMask);
   task->SetEvtClassMin(eventClassMin);
   task->SetEvtClassMax(eventClassMax);
   task->SetCentMin(0.);
   task->SetCentMax(100.);
   task->SetFilterMask(filterMask); 
   task->SetFlagHardest(kHardest);
   task->SetTTLowRef(kTTminr);
   task->SetTTUpRef(kTTmaxr);
   task->SetTTLowSig(kTTmins);
   task->SetTTUpSig(kTTmaxs);
   task->SetFlagPhiBkg(kPhiBkg);   
   task->SetJetPtMin(0.);   
   //task->SetAngStructCloseTracks(1);

 

   mgr->AddTask(task);


   AliAnalysisDataContainer *coutputJetCore = mgr->CreateContainer(Form("pwgjejetcore_%s_%s_%d_%d_%f%f%f%f_%f",bRec1,bRec2,kTriggerMask,kHardest,kTTminr,kTTmaxr,kTTmins,kTTmaxs,kPhiBkg), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_jetcore_%s_%s_%d_%d_%f%f%f%f_%f",AliAnalysisManager::GetCommonFileName(),bRec1,bRec2,kTriggerMask,kHardest,kTTminr,kTTmaxr,kTTmins,kTTmaxs,kPhiBkg));





   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetCore);

   return task;
}
