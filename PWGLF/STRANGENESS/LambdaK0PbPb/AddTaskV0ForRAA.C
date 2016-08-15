AliAnalysisTaskV0ForRAA *AddTaskV0ForRAA(Bool_t anaPP=kFALSE, Bool_t wSDD=kFALSE,Int_t cent=0,Int_t centDet=1,Int_t centRange=0, Bool_t mcMode=kFALSE, Bool_t mcTruthMode=kFALSE,Bool_t usePID=kFALSE,Double_t radCut=0.0,const Char_t * addname=""){
   
  
  
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

   cutsname += addname;
   taskname += addname;
   outname  += addname;
   
   
   AliAnalysisTaskV0ForRAA *task = new AliAnalysisTaskV0ForRAA(taskname);

   Double_t minPt=0.0;
   //task->SetESDTrackCuts(70,4,kTRUE);
   //task->SetESDTrackCutsCharged(70,4,kTRUE);
   //task->SetESDTrackCutsLowPt(70,4,kTRUE);
  

   //Add cuts to task
  

   //--- analysis modes ---//
   task->SetAnapp(anaPP);
   task->SetMCMode(mcMode);
   task->SetMCTruthMode(mcTruthMode);
   task->SelectWithSDD(wSDD);
   //---------- cuts -------------//
   //general cuts
   // task->SetUseOnthefly(kTRUE);
   task->SetUsePID(usePID,3.0,100.0);
   task->SetPrimVertexZCut(10.0,kTRUE);
 
   //rapidity
   task->SetRapidityCutMother(kTRUE,0.5);
   //task->SetDoEtaOfMCDaughtersCut(kFALSE,0.8);
   
   //TPC cuts
   // task->SetCutMoreNclsThanRows(kTRUE);
   // task->SetCutMoreNclsThanFindable(kTRUE);
   task->SetLowPtTPCCutAliESDTrackCut(-1.0);
   //  task->SetRatioFoundOverFindable(0.5);

   //V0 specific cuts
   //task->SetCosOfPointingAngleK(0.99,1000.0);
   //task->SetCosOfPointingAngleL(0.998,1000.0);


   //task->SetArmenterosCutQt(-1.0,6.0,kTRUE,kFALSE);
   
   //task->SetDCAV0ToVertexK0(0.4);
   //task->SetDCAV0ToVertexL(1.2);

   //task->SetDCADaughtersK0(0.23);
   //task->SetDCADaughtersL(0.35);
   //task->SetDCADaughtersAL(0.35);
   
   task->SetDecayRadiusXYMinMax(radCut,1000.0);

   
   //--- centrality ---//
   task->SetUseCentrality(centDet);        // 0=off, 1=VZERO, 2=SPD
   task->SetUseCentralityBin(cent);        // bin to be used 0,5,10,20,30,40,50,60,70,80,90,(100=SPDonly)
   task->SetUseCentralityRange(centRange); // Add centrality bin for increasing original bin range. 
                                           // For cent 60-80%: cent = 60 and centRange = 10
   
   task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kINT7);
   
   mgr->AddTask(task);
 
   
   //================================================
   //              data containers
   //================================================
   //            find input container
   
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

   AliAnalysisDataContainer *coutput1 = 
      mgr->CreateContainer(outname, TList::Class(),
			   AliAnalysisManager::kOutputContainer,Form("%s:simones", AliAnalysisManager::GetCommonFileName()));
   
   //--- connect containers ---//
   mgr->ConnectInput  (task,  0, cinput );
   mgr->ConnectOutput (task,  1, coutput1);
   
   AliLog::SetClassDebugLevel("AliAnalysisTaskV0ForRAA",2);

   return task;
}
