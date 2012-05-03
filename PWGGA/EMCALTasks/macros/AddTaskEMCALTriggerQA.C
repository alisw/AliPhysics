// $Id$
AliAnalysisTaskEMCALTriggerQA * AddTaskEMCALTriggerQA(TString outputFile = "")
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
    
  AliAnalysisTaskEMCALTriggerQA * qatrigger = new AliAnalysisTaskEMCALTriggerQA("QATrigger");  
  
  if(outputFile.Length()==0)outputFile = AliAnalysisManager::GetCommonFileName(); 

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
