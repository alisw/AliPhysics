AliAnalysisTaskPWG4PidDetEx *AddTaskPWG4PidDetEx()
{
// Creates a jet fider task, configures it and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskPWG4PidDetEx", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   
   AliAnalysisTaskPWG4PidDetEx* taskPid = new AliAnalysisTaskPWG4PidDetEx("TaskPID");
   if(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler"))taskPid->SetAnalysisType("AOD");
   else {
     // Either we assume taht w have the cuts already loaded or we just do it again
     // this can be done more elegant when also accepting AOD filtered output
     gROOT->LoadMacro("AddTaskESDfilter.C");
     AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
     trackFilter->AddCuts(CreateCuts(0));
     taskPid->SetTrackFilter(trackFilter);
   }
   mgr->AddTask(taskPid);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_PID = mgr->CreateContainer("histosPID", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4pid.root");

   mgr->ConnectInput  (taskPid, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (taskPid, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (taskPid, 1, coutput1_PID);
   return taskPid;
}

AliAnalysisTaskPWG4PidDetEx *AddTaskPWG4PidDetEx(AliAnalysisManager* mgr,AliAnalysisDataContainer *cinput)
{
  // This is only for running on PROOF with the old root version 5-22-00 
  // and the older version of the AF

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskJets", "No analysis manager to connect to.");
      return NULL;
   }  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskPWG4PidDetEx", "This task requires an input event handler");
      return NULL;
   }

   // Create the task and configure it.
   //===========================================================================
   AliAnalysisTaskPWG4PidDetEx* taskPid = new AliAnalysisTaskPWG4PidDetEx("TaskPID");
   if(mgr->GetInputEventHandler()->InheritsFrom("AliAODInputHandler"))taskPid->SetAnalysisType("AOD");
   else {
     // Either we assume taht w have the cuts already loaded or we just do it again
     // this can be done more elegant when also accepting AOD filtered output
     gROOT->LoadMacro("AddTaskESDfilter.C");
     AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
     trackFilter->AddCuts(CreateCuts(0));
     taskPid->SetTrackFilter(trackFilter);
   }
   mgr->AddTask(taskPid);


   //
   // Create containers for input/output
   AliAnalysisDataContainer *c_aod_pid = mgr->CreateContainer("cAODpid", TTree::Class(),AliAnalysisManager::kExchangeContainer);
   AliAnalysisDataContainer *coutput1_PID = mgr->CreateContainer("histosPID", TList::Class(),AliAnalysisManager::kOutputContainer,"pwg4pid.root");
   mgr->ConnectInput  (taskPid,  0, cinput  );
   mgr->ConnectOutput (taskPid,  0, c_aod_pid);
   mgr->ConnectOutput (taskPid,  1, coutput1_PID);

   return taskPid;

}  
