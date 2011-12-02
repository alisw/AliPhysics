

AliAnalysisTaskJetCore* AddTaskJetCore(const char* bRec1,const char* bRec2, UInt_t filterMask = 272 , Float_t ptTrackMin = 0.15,  Int_t eventClassMin = 0, Int_t eventClassMax = 4){

   Printf("adding task jet core\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetCore", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetCore", "This task requires an input event handler.");
      return NULL;
   }

     
   //jetsin bRec1 branch have larger radius than those in bRec2  

  TString typeRec(bRec1);
  TString typeGen(bRec2);
      
  AliAnalysisTaskJetCore *task = new AliAnalysisTaskJetCore(Form("JetCore_%s_%s",bRec1,bRec2));
   


   task->SetBranchNames(bRec1,bRec2);
   task->SetOfflineTrgMask(AliVEvent::kMB);

   task->SetEvtClassMin(eventClassMin);
   task->SetEvtClassMax(eventClassMax);
   task->SetCentMin(0.);
   task->SetCentMax(100.);

 
   
   task->SetJetPtMin(0.);   
   
   //task->SetRadioFrac(0.2);   //radius of the concentric cone
                                //within which you sum up the pT of
                                //the tracks to compute the core fraction
                                //in method3. It is also the maxium distance
                                //between jet axis in method2. Default is 0.2 

   //task->SetMinDist(0.1);     //minimum distance between jets to be 
                                //concentric in method1. Default is 0.1

 
   task->SetBackgroundBranch("jeteventbackground_clustersAOD_KT04_B0_Filter00272_Cut00150_Skip00"); 

   mgr->AddTask(task);


   AliAnalysisDataContainer *coutputJetCore = mgr->CreateContainer(Form("pwg4jetcore_%s_%s",bRec1,bRec2), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWG4_jetcore_%s_%s",AliAnalysisManager::GetCommonFileName(),bRec1,bRec2));





   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetCore);

   return task;
}
