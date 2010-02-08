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
   //Task for TPC tracks 
   AliAnalysisTaskQASym *task1 = new AliAnalysisTaskQASym("AliAnalysisTaskQASym_TPC");
   task1->SetTrackType(1);
   task1->SelectCollisionCandidates();

   //cuts for global tracks
   AliESDtrackCuts* esdTrackCutsL0 = new AliESDtrackCuts("AliESDtrackCuts0","Global");
   esdTrackCutsL0->SetMinNClustersTPC(70);
   esdTrackCutsL0->SetRequireTPCRefit(kTRUE);
   esdTrackCutsL0->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL0->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL0->SetAcceptKinkDaughters(kFALSE);
   
   //cuts for TPC tracks
   AliESDtrackCuts* esdTrackCutsL1 = new AliESDtrackCuts("AliESDtrackCuts1","TPC");
   esdTrackCutsL1->SetRequireTPCRefit(kFALSE);
   esdTrackCutsL1->SetAcceptKinkDaughters(kFALSE);
   //jacek's cuts:
   esdTrackCutsL1->SetMinNClustersTPC(70);
   // cut on max ncl=160 in Task
   esdTrackCutsL1->SetMaxDCAToVertexXY(3.);
   esdTrackCutsL1->SetMaxDCAToVertexZ(3.);
   esdTrackCutsL1->SetMaxChi2PerClusterTPC(3.999);
   //cut minChi=0 in task
   //esdTrackCutsL1->SetPRange(0.15,16); // not needed for QA
   //esdTrackCutsL1->SetEtaRange(-0.8, 0.7999); // not needed for QA
  

   task0->SetCuts(esdTrackCutsL0);
   task1->SetCuts(esdTrackCutsL1);

   mgr->AddTask(task0);
   mgr->AddTask(task1);
  
   AliAnalysisDataContainer *cout0  = 0;
   AliAnalysisDataContainer *cout1  = 0;
   
   if(runNumber>0){ 
    cout0 =  mgr->CreateContainer("QAsymHists_Global",TList::Class(),
				  AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
    cout1 =  mgr->CreateContainer("QAsymHists_TPC",TList::Class(),
				  AliAnalysisManager::kOutputContainer, Form("run%d.root",runNumber));
   }
   
   else{
      cout0 = mgr->CreateContainer("QAsymHists_Global",TList::Class(),
				 AliAnalysisManager::kOutputContainer, 
				 Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));
      cout1 = mgr->CreateContainer("QAsymHists_TPC",TList::Class(),
				   AliAnalysisManager::kOutputContainer, 
				 Form("%s:PWG1_QAsymHists",AliAnalysisManager::GetCommonFileName()));

   }


   mgr->ConnectInput  (task0, 0, mgr->GetCommonInputContainer());
   mgr->ConnectInput  (task1, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task0, 1, cout0);
   mgr->ConnectOutput (task1, 1, cout1);
  
   return task1;

}


