//=========================================================================//
//                                                                         //
//              A template Analysis AddTask for PMD analysis               //
//            You can copy it and add features to ir as your need          //
//                                                                         //
//                               Satyajit Jena                             //
//                               sjena@cern.ch                             //
//                                13/04/2012                               //
//                                                                         //
//=========================================================================//


void AddAliPMDAnalysisTaskPbPb(const Char_t *taskname="MyTask", Bool_t isMC) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }

    
  AliPMDAnalysisTaskPbPb *task = new AliPMDAnalysisTaskPbPb(taskname);

  if (!task) {
    Error("AliPMDAnalysisTaskPbPb", "Task could not be created.");
    return NULL;
  }
  
  task->SetIsMC(isMC); // Similarly you can add more functions
                       // See one of my task in aliroot
  mgr->AddTask(task);
 
  TString commonname   = Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput 
    = mgr->CreateContainer("PMD", TList::Class(),    
			   AliAnalysisManager::kOutputContainer, commonname);

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  return;
};
