AliAnalysisTaskQAMultistrange *AddTaskQAMultistrange(Bool_t isMC = kFALSE) {

   // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQAMultistrange", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskQAMultistrange", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskQAMultistrange *taskcheckcascade = new AliAnalysisTaskQAMultistrange("TaskCheckCascade");
     taskcheckcascade->SetIsMC                       (isMC);
     taskcheckcascade->SetAnalysisType               (type);

   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWGLFStrangeness.outputCheckCascade";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("fListHistMultistrangeQA",
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cfcontCuts",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cfcontMCCuts",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("cfcontMCgen",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   
   mgr->ConnectInput( taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascade, 2, coutput2);  
   mgr->ConnectOutput(taskcheckcascade, 3, coutput3);  
   mgr->ConnectOutput(taskcheckcascade, 4, coutput4);
 
   return taskcheckcascade;
}   
