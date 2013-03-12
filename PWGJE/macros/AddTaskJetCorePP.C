

AliAnalysisTaskJetCorePP* AddTaskJetCorePP(
   const Char_t* branchPrefix="clustersAOD",
   const Char_t* jetAlgo="ANTIKT",
   Float_t jetParameterR = 0.4,  //jet R 
   Int_t   bgMode = 0, 
   UInt_t  trkFilterMask = 272, 
   Float_t trackLowPtCut = 0.15,
   Int_t   skipJet = 0, 
   Int_t   collisionSystem = 0, //pp=0, pPb=1
   Int_t   offlineTriggerMask=AliVEvent::kMB, //MinBias=0 
   Int_t   minContribVtx = 1,
   Float_t vtxZMin = -10.0,
   Float_t vtxZMax = 10.0,
   Float_t centMin = 0.0,
   Float_t centMax = 100.0,
   Float_t triggerEtaCut = 0.9,
   Float_t trackEtaCut = 0.9,
   const Char_t* nonStdFile=""
  ){ 

   Printf("adding task jet response\n");

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskJetCorePP", "No analysis manager to connect to.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskJetCorePP", "This task requires an input event handler.");
      return NULL;
   }

   Float_t jetEtaMin = -0.9 + jetParameterR;
   Float_t jetEtaMax =  0.9 - jetParameterR; 
    
   TString analBranch(branchPrefix);
   TString stJetAlgo(jetAlgo);
   stJetAlgo.ToUpper();
   analBranch = analBranch + "_" + stJetAlgo + Form("%02d",(Int_t) (10*jetParameterR));
   analBranch = analBranch + Form("_B%d",(Int_t) bgMode);
   analBranch = analBranch + Form("_Filter%05d",(UInt_t) trkFilterMask);
   analBranch = analBranch + Form("_Cut%05d",(Int_t) (1000*trackLowPtCut));
   if(analBranch.BeginsWith("clustersAOD"))
      analBranch = analBranch + Form("_Skip%02d",(Int_t) skipJet);
   //clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00   
   //Skip00 none of the most energetic jets is ommited
   //Cut00150  pT min cut on track
   //Filter00272
  
   AliAnalysisTaskJetCorePP *task = new AliAnalysisTaskJetCorePP(Form("JetCorePP_%s_%d",analBranch.Data(),offlineTriggerMask));

   task->SetBranchName(analBranch.Data());
   task->SetNonStdFile(nonStdFile);
   task->SetSystem(collisionSystem); 
   task->SetJetR(jetParameterR);
   task->SetOfflineTrgMask(offlineTriggerMask);
   task->SetMinContribVtx(minContribVtx);
   task->SetVtxZMin(vtxZMin);
   task->SetVtxZMax(vtxZMax);
   task->SetFilterMask(trkFilterMask);
   task->SetCentMin(centMin);
   task->SetCentMax(centMax);
   task->SetJetEtaMin(jetEtaMin);
   task->SetJetEtaMax(jetEtaMax);
   task->SetTriggerEtaCut(triggerEtaCut);
   task->SetTrackEtaCut(trackEtaCut);
   task->SetTrackLowPtCut(trackLowPtCut);
   task->SetDebugLevel(0); //No debug messages 0
   mgr->AddTask(task);

   AliAnalysisDataContainer *coutputJetCorePP = mgr->CreateContainer(
      Form("pwgjejetcorepp_%s_%d",analBranch.Data(),offlineTriggerMask), 
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:PWGJE_jetcorepp_%s_%d",AliAnalysisManager::GetCommonFileName(),analBranch.Data(),offlineTriggerMask)
   );

   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputJetCorePP);

   return task;
}
