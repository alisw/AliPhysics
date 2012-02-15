

AliAnalysisTaskJetCore* AddTaskJetCore(const char* bRec1,const char* bRec2, UInt_t filterMask = 272 , Float_t ptTrackMin = 0.15,  Int_t eventClassMin = 0, Int_t eventClassMax = 4){

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
      
  AliAnalysisTaskJetCore *task = new AliAnalysisTaskJetCore(Form("JetCore_%s_%s",bRec1,bRec2));
   


   task->SetBranchNames(bRec1,bRec2);
   //task->SetOfflineTrgMask(AliVEvent::kMB);

   task->SetEvtClassMin(eventClassMin);
   task->SetEvtClassMax(eventClassMax);
   task->SetCentMin(0.);
   task->SetCentMax(100.);

 
   
   task->SetJetPtMin(0.);   
   //task->SetAngStructCloseTracks(1);

 

   mgr->AddTask(task);


   AliAnalysisDataContainer *coutputJetCore = mgr->CreateContainer(Form("pwgjejetcore_%s_%s",bRec1,bRec2), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_jetcore_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec1,bRec2));





   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetCore);

   return task;
}
