AliAnalysisDeuteronTree *AddTaskDeuteronTree( const TString lMasterJobSessionFlag = "")
{
    ::Info("AddTaskDeuteronTree","Entering AddTask macro");

    
    // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskDeuteronTree", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskDeuteronTree", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

    ::Info("AddTaskDeuteronTree","About to create new task");

    
   // Create and configure the task
   AliAnalysisDeuteronTree *taskDeuteronTree = new AliAnalysisDeuteronTree("taskDeuteronTree");

  ::Info("AddTaskDeuteronTree","About to add task to manager");
 
   mgr->AddTask(taskDeuteronTree);
    
::Info("AddTaskDeuteronTree","Task added to manager");
 
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWGLFDeuteronTree";

   outputFileName += "_PPB";
   
   Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutputList = mgr->CreateContainer("clistDeuteron",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTree",
							     TTree::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
  
   //This one you should merge in file-resident ways...
   coutputTree->SetSpecialOutput();

   mgr->ConnectInput( taskDeuteronTree, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskDeuteronTree, 1, coutputList);
   mgr->ConnectOutput(taskDeuteronTree, 2, coutputTree);
   
   return taskDeuteronTree;
}   
