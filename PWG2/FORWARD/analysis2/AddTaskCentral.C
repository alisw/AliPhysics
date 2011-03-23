/**
 * @file   AddTaskCentral.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:13:25 2011
 * 
 * @brief  
 * 
 * 
 * @ingroup pwg2_forward_scripts_tasks
 */

/**
 * This is the macro to include the Central multiplicity in a train.  
 * 
 * @ingroup pwg2_forward_aod
 */
AliAnalysisTask* AddTaskCentral(UShort_t sys=0, UShort_t sNN=0, 
				Short_t field=0)
{
  gSystem->Load("libPWG2forward2");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCentral", "No analysis manager to connect to.");
    return NULL;
  }   

  // --- Make the task and add it to the manager ---------------------
  AliCentralMultiplicityTask* task = new AliCentralMultiplicityTask("Central");
  if(sys>0 && sNN > 0)
    task->GetManager().Init(sys, sNN, field);
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
//
// EOF
//
