AliAnalysisTaskTagCreator *AddTaskTagCreation()
{

// Creates tag AOD files

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskTagCreator ", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // Check input/output handlers
   //===============================================================================
   AliAODInputHandler *aod_input_h = (AliAODInputHandler*)mgr->GetInputEventHandler();
   if (!aod_input_h) {
      ::Error("AddTaskTagCreator ", "Input AOD handler does not exist!.");
      return NULL;
   }
      
   // Create the task, add it to the manager and configure it.
   //===========================================================================   

   AliAnalysisTaskTagCreator *tagcreator  = new AliAnalysisTaskTagCreator("TagCreator");
   mgr->AddTask(tagcreator);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
    // Tag container							      
    AliAnalysisDataContainer *cout_tags = mgr->CreateContainer("cTag",TTree::Class(), AliAnalysisManager::kOutputContainer, "AOD.tag.root");
    
    mgr->ConnectInput  (tagcreator ,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput (tagcreator ,  1, cout_tags);

   return tagcreator ;
}   
