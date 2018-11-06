AliAnalysisTaskSVtaskMCFilter* AddTaskSVtaskMCFilter(const char* trkcontname   = "tracks", const char* outtrk = "mytracks"){

   //Steering macro for filter of detector level MC tracks. Removes EPOS part of the event

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
     ::Error("AddSVtaskMCFilter", "No analysis manager to connect to.");
     return NULL;
   }

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()){
      ::Error("AddSVtaskMCFilter", "This task requires an input event handler");
      return NULL;
   }



   TString taskFiltername("SVtaskMCFilter");

   AliAnalysisTaskSVtaskMCFilter* task = mgr->GetTask(taskFiltername.Data());
   if(task){
      ::Info("AddTaskSVtaskMCFilter", Form("Task %s already exist, continue",taskFiltername.Data()));
      return task;
   }else{
      ::Info("AddTaskSVtaskMCFilter", "Creating the task");
  
      // create the task
      task = new AliAnalysisTaskSVtaskMCFilter(taskFiltername.Data());
      AliTrackContainer* trkCont = task->AddTrackContainer(trkcontname);//for data, and reco MC
      
      mgr->AddTask(task);
   }

   task->SetInputTracksName(trkcontname);   
   task->SetFilteredTracksName(outtrk);

   // ------ input data ------
   AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
   mgr->ConnectInput(task, 0, cinput);

   ::Info("AddTaskSVtaskMCFilter", "Input and Output connected to the manager");
   return task;
}
