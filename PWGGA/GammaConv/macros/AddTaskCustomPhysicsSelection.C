void LHC15fPS(AliPhysicsSelection* ps)
{
  // Defaults filling scheme 
  AliOADBFillingScheme * cFS = new AliOADBFillingScheme("Default");
  cFS->SetFillingSchemeName("Default");
  cFS->SetBXIDs("B","");
  cFS->SetBXIDs("A","" );
  cFS->SetBXIDs("AC","" );
  cFS->SetBXIDs("ACE","" );
  cFS->SetBXIDs("C","");
  cFS->SetBXIDs("E",""       );

  // DefaultPP
  AliOADBPhysicsSelection * cPS = new AliOADBPhysicsSelection("cPS");
  Int_t id = 0;
  cPS->AddCollisionTriggerClass(AliVEvent::kMB,"+CINT5-B-NOPF-ALLNOTRD,+C0SMB-B-NOPF-ALLNOTRD","B",      id);
  cPS->AddBGTriggerClass       (AliVEvent::kMB,"+CINT5-A-NOPF-ALLNOTRD,+C0SMB-A-NOPF-ALLNOTRD","A",      id);
  cPS->AddBGTriggerClass       (AliVEvent::kMB,"+CINT5-C-NOPF-ALLNOTRD,+C0SMB-C-NOPF-ALLNOTRD","C",      id);
  cPS->AddBGTriggerClass       (AliVEvent::kMB,"+CINT5-E-NOPF-ALLNOTRD,+C0SMB-E-NOPF-ALLNOTRD","E",      id);
  cPS->AddBGTriggerClass       (AliVEvent::kMB,"+CINT5-ACE-NOPF-ALLNOTRD,+C0SMB-ACE-NOPF-ALLNOTRD","ACE",id);
  cPS->SetHardwareTrigger      (id,"1");
  cPS->SetOfflineTrigger       (id,"1");

  // Trigger analysis defaults
  AliOADBTriggerAnalysis * cTA = new AliOADBTriggerAnalysis("Default");
  cTA->SetZDCCorrParameters(0.5, 0, 4*0.7, 4*0.7);

  ps->SetCustomOADBObjects(cPS, cFS, cTA);
}

AliPhysicsSelectionTask* AddTaskCustomPhysicsSelection(Bool_t mCAnalysisFlag = kFALSE, Bool_t deprecatedFlag = kTRUE, UInt_t computeBG = 0, Bool_t useSpecialOutput=kFALSE) 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }

  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  
  TString inputDataType = inputHandler->GetDataType(); // can be "ESD" or "AOD"

  // Configure analysis
  //===========================================================================
AliPhysicsSelectionTask *task = new AliPhysicsSelectionTask("");
  task->SetUseSpecialOutput(useSpecialOutput); // RS: optionally use special output
  // this makes physics selection to work using AliMultiInputEventHandler
  if (inputHandler && (inputHandler->IsA() == AliMultiInputEventHandler::Class())) {
    AliMultiInputEventHandler *multiInputHandler=(AliMultiInputEventHandler*)inputHandler;
    AliInputEventHandler *ih = multiInputHandler->GetFirstInputEventHandler();
    if (!ih) {
      ::Error("AddTaskPhysicsSelection","ESD or AOD input handler is missing");
      return NULL;
    }
    ih->SetEventSelection(multiInputHandler->GetEventSelection());
    inputDataType = ih->GetDataType(); // can be "ESD" or "AOD"
  }

  mgr->AddTask(task);

  AliPhysicsSelection* physSel = task->GetPhysicsSelection();
  LHC15fPS(physSel);
  if (mCAnalysisFlag)      
    physSel->SetAnalyzeMC();
  if (computeBG)
    physSel->SetComputeBG(computeBG);

  if(!deprecatedFlag) 
    AliFatal("The BG ID flag is deprecated. Please use the OADB to configure the cuts");

  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cstatsout",
                TList::Class(),
                AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
  //		
  if (useSpecialOutput) coutput1->SetSpecialOutput(); // RS: optionally use special output
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);

  return task;
}   
