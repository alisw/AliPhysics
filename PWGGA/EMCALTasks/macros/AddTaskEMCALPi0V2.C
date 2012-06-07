// $Id: AddTaskEMCALPi0V2.C 56081 2012-05-01 08:57:08Z loizides $

AliAnalysisTaskPi0V2 *AddTaskEMCALPi0V2()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALPi0V2", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskPi0V2* ana = new  AliAnalysisTaskPi0V2("Pi0v2Task");
  ana->SelectCollisionCandidates( AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEMCALP0v2", 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
