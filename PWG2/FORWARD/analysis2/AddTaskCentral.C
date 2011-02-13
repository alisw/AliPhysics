/** 
 * @defgroup pwg2_forward_scripts Scripts used in the analysis
 *
 * @ingroup pwg2_forward
 */
/**
 * @file 
 * @ingroup pwg2_forward_scripts
 * 
 */
/**
 * This is the macro to include the Forward multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_scripts
 */
AliAnalysisTask* AddTaskCentral()
{
  gSystem->Load("libPWG2forward2");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCentral", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliCentralMultiplicityTask* task = new AliCentralMultiplicityTask("Central");
  task->InitManager(1,900,5);
  mgr->AddTask(task);
  
  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Central", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);
  
  return task;
}
