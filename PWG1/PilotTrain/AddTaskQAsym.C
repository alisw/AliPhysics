AliAnalysisTaskQASym * AddTaskQAsym(Int_t runNumber)

{
  // Creates a QA task exploiting simple symmetries phi, eta +/-, charge ...
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskQAsym", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTasQAsym", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
   
 
   //Task for global tracks
   AliAnalysisTaskQASym *task0 = new AliAnalysisTaskQASym("AliAnalysisTaskQASym_Global");
   task0->SetTrackType(0);
   task0->SelectCollisionCandidates();
   //Task for ITS tracks 
   AliAnalysisTaskQASym *task1 = new AliAnalysisTaskQASym("AliAnalysisTaskQASym_ITS");
   task1->SetTrackType(1);
   task1->SetStandAloneTrack(kFALSE);
   task1->SelectCollisionCandidates();
   //Task for ITS tracks SA
   AliAnalysisTaskQASym *task1sa = new AliAnalysisTaskQASym("AliAnalysisTaskQASym_ITS_SA");
   task1sa->SetTrackType(1);
   task1sa->SetStandAloneTrack(kTRUE);
   task1sa->SelectCollisionCandidates();
   //Task for TPC tracks 
   AliAnalysisTaskQASym *task2 = new AliAnalysisTaskQASym("AliAnalysisTaskQASym_TPC");
   task2->SetTrackType(2);
   task2->SelectCollisionCandidates();

   //cuts for global tracks
   AliESDtrackCuts* esdTrackCutsL0 = new AliESDtrackCuts("AliESDtrackCuts0","Global");
   esdTrackCutsL0->SetMinNClustersTPC(70);
   esdTrackCutsL0->SetRequireTPCRefit(kTRUE);
   esdTrackCutsL0->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL0->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL0->SetAcceptKinkDaughters(kFALSE);
   
   //cuts for ITS tracks
   AliESDtrackCuts* esdTrackCutsL1 = new AliESDtrackCuts("AliESDtrackCuts1","ITS");
   esdTrackCutsL1->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL1->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL1->SetAcceptKinkDaughters(kFALSE);
   esdTrackCutsL1->SetRequireITSRefit(kTRUE);
   esdTrackCutsL1->SetRequireITSStandAlone(kTRUE, kFALSE);

   //cuts for ITS tracks SA
   AliESDtrackCuts* esdTrackCutsL1sa = new AliESDtrackCuts("AliESDtrackCuts1","ITS_SA");
   esdTrackCutsL1sa->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL1sa->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL1sa->SetAcceptKinkDaughters(kFALSE);
   esdTrackCutsL1sa->SetRequireITSRefit(kTRUE);
   // esdTrackCutsL1sa->SetRequireITSStandAlone(kTRUE, kTRUE); //cut on SA tracks in AliAnalysisTaskQASym
   
   //cuts for TPC tracks
   AliESDtrackCuts* esdTrackCutsL2 = new AliESDtrackCuts("AliESDtrackCuts2","TPC");
   esdTrackCutsL2->SetRequireTPCRefit(kFALSE);
   esdTrackCutsL2->SetAcceptKinkDaughters(kFALSE);
   //jacek's cuts:
   esdTrackCutsL2->SetMinNClustersTPC(70);
   // cut on max ncl=160 in Task
   esdTrackCutsL2->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL2->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL2->SetMaxChi2PerClusterTPC(3.999);
   //cut minChi=0 in task
   //esdTrackCutsL2->SetPRange(0.15,16); // not needed for QA
   //esdTrackCutsL2->SetEtaRange(-0.8, 0.7999); // not needed for QA
  

   task0->SetCuts(esdTrackCutsL0);
   task1->SetCuts(esdTrackCutsL1);
   task1sa->SetCuts(esdTrackCutsL1sa);
   task2->SetCuts(esdTrackCutsL2);

   mgr->AddTask(task0);
   mgr->AddTask(task1);
   mgr->AddTask(task1sa);
   mgr->AddTask(task2);
  
   AliAnalysisDataContainer *cout0  = 0;
   AliAnalysisDataContainer *cout1  = 0;
   AliAnalysisDataContainer *cout1sa  = 0;
   AliAnalysisDataContainer *cout2  = 0;
   
   if(runNumber>0){ 
    cout0 =  mgr->CreateContainer("QAsymHists_Global",TList::Class(),
				  AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
    cout1 =  mgr->CreateContainer("QAsymHists_ITS",TList::Class(),
				  AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
    cout1sa =  mgr->CreateContainer("QAsymHists_ITS_SA",TList::Class(),
				    AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
    cout2 =  mgr->CreateContainer("QAsymHists_TPC",TList::Class(),
				  AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
   }
   
   else{
      cout0 = mgr->CreateContainer("QAsymHists_Global",TList::Class(),
				 AliAnalysisManager::kOutputContainer, 
				 Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));
      cout1 = mgr->CreateContainer("QAsymHists_ITS",TList::Class(),
				   AliAnalysisManager::kOutputContainer, 
				 Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));
      cout1sa = mgr->CreateContainer("QAsymHists_ITS_SA",TList::Class(),
				     AliAnalysisManager::kOutputContainer, 
				     Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));
      cout2 = mgr->CreateContainer("QAsymHists_TPC",TList::Class(),
				   AliAnalysisManager::kOutputContainer, 
				 Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));
   }


   mgr->ConnectInput  (task0, 0, mgr->GetCommonInputContainer());
   mgr->ConnectInput  (task1, 0, mgr->GetCommonInputContainer());
   mgr->ConnectInput  (task1sa, 0, mgr->GetCommonInputContainer());
   mgr->ConnectInput  (task2, 0, mgr->GetCommonInputContainer());

   mgr->ConnectOutput (task0, 1, cout0);
   mgr->ConnectOutput (task1, 1, cout1);
   mgr->ConnectOutput (task1sa, 1, cout1sa);
   mgr->ConnectOutput (task2, 1, cout2);
  
   return task0;

}


