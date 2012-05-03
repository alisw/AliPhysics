AliAnalysisTask *AddTaskV0ForRAA(Bool_t anaPP=kFALSE, Int_t cent=0,Int_t centDet=1,Int_t centRange=0, Bool_t mcMode=kFALSE, Bool_t mcTruthMode=kFALSE,Bool_t onFly=kTRUE,Bool_t usePID=kFALSE){
   
  
  
   //--- get the current analysis manager ---//
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_V0ForRAA", "No analysis manager found.");
      return 0;
   }


   // -- check for ESD and MC ---//
   Bool_t hasESD=kFALSE,hasMC=kFALSE;
   AliESDInputHandler *esdH = 
      static_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
   if (esdH) hasESD=kTRUE;
   cout<<"ESD: "<<hasESD<<endl;
   if(!hasESD) return NULL;

   if(mcMode || mcTruthMode){
      AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
	 (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (mcH) hasMC=kTRUE;
      cout<<"MC: "<<hasMC<<endl;
      if(!hasMC) return NULL;
   }
   
   //============= Set Task Name ===================
   TString taskName=("AliAnalysisTaskV0ForRAA.cxx+");
   //===============================================
   //            Load the task
   gROOT->LoadMacro(taskName.Data());
   if (gProof){
      TString taskSO=gSystem->pwd();
      taskSO+="/";
      taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
      gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
   }
  
   //========= Add task to the ANALYSIS manager =====
   TString cutsname = "AliESDtrackCutsV0ForRAA";
   TString taskname = "V0ForRAA";
   TString outname  = "V0ForRAA";

   if(mcMode) {
      cutsname +="_MCreco";
      taskname +="_MCreco";
      outname  +="_MCreco";
   }
   if(mcTruthMode){
      cutsname +="_MCTruth";
      taskname +="_MCTruth";
      outname  +="_MCTruth";
   }
   if(anaPP) {
      cutsname += "_pp";
      taskname += "_pp";
      outname  += "_pp";
   }
   else {
      cutsname +="_cent";
      cutsname += cent;

      taskname +="_cent";
      taskname += cent;

      outname  +="_cent";
      outname  += cent; 
   }

   TString outputname;
   outputname  = outname;
   outputname += ".root";
   
   AliAnalysisTaskV0ForRAA *task = new AliAnalysisTaskV0ForRAA(taskname);

   Double_t minPt=0.0;//15;
   
   //--- esd track cuts V0 daughters ---//
  
   AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts(cutsname);
   //    esdTrackCuts->SetMinNClustersTPC(70);
   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
   esdTrackCuts->SetMinNCrossedRowsTPC(70);
   esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
   
   esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   esdTrackCuts->SetRequireTPCRefit(kTRUE);
   esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
   
   //esdTrackCuts->SetEtaRange(-0.9,+0.9);
   task->SetESDTrackCuts(esdTrackCuts);


   //--- analysis modes ---//
   task->SetAnapp(anaPP);
   task->SetMCMode(mcMode);
   task->SetMCTruthMode(mcTruthMode);
   
   //--- general cuts ---//
   task->SetUseOnthefly(onFly);
   task->SetUsePID(usePID,4.0,2100.0);
   task->SetPrimVertexZCut(10.0,kTRUE);
   task->SetCosOfPointingAngleK(0.99,1000.0);
   task->SetCosOfPointingAngleL(0.99,1000.0);
   task->SetRapidityCutMother(kTRUE,0.6);
   task->SetDoEtaOfMCDaughtersCut(kFALSE,0.9);
   task->SetCtauCut(5.0,3.8,0.5,1.5);
   
   //--- centrality ---//
   task->SetUseCentrality(centDet);     // 0=off, 1=VZERO, 2=SPD
   task->SetUseCentralityBin(cent);     // Bin to be used 0,5,10,20,30,40,50,60,70,80,90,(100=SPDonly)
   task->SetUseCentralityRange(centRange);
   
   task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral );
   
   mgr->AddTask(task);
 
   
   //================================================
   //              data containers
   //================================================
   //            find input container
   
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

   AliAnalysisDataContainer *coutput1 = 
     mgr->CreateContainer(outname, TList::Class(),
			  AliAnalysisManager::kOutputContainer,outputname);
   
   //--- connect containers ---//
   mgr->ConnectInput  (task,  0, cinput );
   mgr->ConnectOutput (task,  1, coutput1);
   
   
   return task;
}
