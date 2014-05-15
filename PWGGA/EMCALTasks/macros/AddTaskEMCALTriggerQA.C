// $Id$
AliAnalysisTaskEMCALTriggerQA * AddTaskEMCALTriggerQA(Bool_t kSimulation = kFALSE, TString outputFile = "", )
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALTriggerQA", "This task requires an input event handler");
    return NULL;
  }
  
  Bool_t kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  if(kUseKinematics) kSimulation = kTRUE;
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  AliAnalysisTaskEMCALTriggerQA * qatrigger = new AliAnalysisTaskEMCALTriggerQA("QATrigger");  
  if(kSimulation) qatrigger->SwitchOnMCData();
  else            qatrigger->SwitchOffMCData();
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 

  if(inputDataType.Contains("AOD")) qatrigger->SwitchOffV0SignalHistograms();
  else                              qatrigger->SwitchOnV0SignalHistograms();
  
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = 0;
  if(outputFile.Length()==0)
    coutput = mgr->CreateContainer("EMCALQATrigger", TList::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:EMCALQATrigger",outputFile.Data()));
  else
    coutput = mgr->CreateContainer("EMCALQATrigger", TList::Class(), AliAnalysisManager::kOutputContainer,  outputFile.Data());

  mgr->AddTask(qatrigger);
  mgr->ConnectInput  (qatrigger, 0, cinput1);
  mgr->ConnectOutput (qatrigger, 1, coutput);
  
  return qatrigger;
}
