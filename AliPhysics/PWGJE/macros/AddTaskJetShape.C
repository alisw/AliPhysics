AliAnalysisTaskJetShape* AddTaskJetShape(const char* bRec1, const char* bRec2, const char* bBkg1, const char *bBkg2,
                                         UInt_t filterMask = 272 ,Bool_t kIsPbPb = kFALSE, UInt_t kTriggerMask=0, Int_t eventClassMin = 0, Int_t eventClassMax = 4)
{

   Printf("adding task jet shape\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetShape", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetShape", "This task requires an input event handler.");
      return NULL;
   }

     
  

  TString sRec(bRec1);
  TString sGen(bRec2);
  TString sRecBkg(bBkg1);
  TString sGenBkg(bBkg2);
      
  AliAnalysisTaskJetShape *task = new AliAnalysisTaskJetShape(Form("JetShape_%s_%s_%d",bRec1,bRec2,kTriggerMask));
   
   task->SetBranchNames(sRec,sGen);
   task->SetBackgroundBranch(sRecBkg, sGenBkg);
   task->SetOfflineTrgMask(kTriggerMask);
   task->SetEvtClassMin(eventClassMin);
   task->SetEvtClassMax(eventClassMax);
   task->SetCentMin(0.);
   task->SetCentMax(100.);
   task->SetFilterMask(filterMask); 
   task->SetJetPtCorrMin(20.,20);   
   task->SetPbPb(kIsPbPb);
   mgr->AddTask(task);


   AliAnalysisDataContainer *coutputJetShape = mgr->CreateContainer(Form("pwgjeJetShape_%s_%s_%s_%s_%d",bRec1,bRec2,bBkg1,bBkg2,kTriggerMask), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_JetShape_%s_%s_%d",AliAnalysisManager::GetCommonFileName(),bRec1,bRec2,kTriggerMask));


   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetShape);

   return task;
}
