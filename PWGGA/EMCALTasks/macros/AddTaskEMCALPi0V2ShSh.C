AliAnalysisTaskEMCALPi0V2ShSh *AddTaskEMCALPi0V2ShSh(TString arrayName = "V1_Ecell150_Eseed300_DT0_WT0",
						     AliVEvent::EOfflineTriggerTypes trig = AliVEvent::kCentral + AliVEvent::kSemiCentral + AliVEvent::kMB + AliVEvent::kEMCEGA,
						     Bool_t isPhosCali = kTRUE, Bool_t isCentFlat = kTRUE,
                                                     const Int_t debug = 0
                                                    )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr) {
    ::Error("AddTaskEMCALPi0V2ShSh", "No analysis manager to connect to.");
    return NULL;
  }  
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventplane", "This task requires an input event handler");
    return NULL;
  }		
  // TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  Bool_t ismc = kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
  cout << "AddTaskEMCALPi0V2ShSh - MC config is: " << ismc << endl;  
  if (ismc) return 0;
  
  AliAnalysisTaskEMCALPi0V2ShSh *task = new AliAnalysisTaskEMCALPi0V2ShSh("EMCALPi0V2ShSh");
  task->SetEMCALClusterListName(arrayName);
  task->SelectCollisionCandidates(trig);
  task->IsPHOSCali(isPhosCali);
  task->IsCentFlat(isCentFlat);
  task->SetDebugLevel(debug);
  
  if (!ismc) {		
    mgr->AddTask(task);
    
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("hist", TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s",AliAnalysisManager::GetCommonFileName()));
      
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutput1);
    mgr->SetDebugLevel(debug);
  }
  return task;
}
