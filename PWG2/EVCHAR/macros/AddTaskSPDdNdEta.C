//=============================================================================
//
// *** AddTaskSPDdNdEta.C ***
//
// This macro initialize a complete AnalysisTask object for SPD dNdEta analysis.
//
//=============================================================================

AliAnalysisTaskSPDdNdEta *AddTaskSPDdNdEta()
{
// Creates an analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskSPDdNdEta", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskSPDdNdEta", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // C. Create the task, add it to manager.
  //==============================================================================

  AliAnalysisTaskSPDdNdEta *taskSPDdNdEta = new AliAnalysisTaskSPDdNdEta("TaskSPDdNdEta");
  mgr->AddTask(taskSPDdNdEta);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================
  taskSPDdNdEta->SetReadMC(kTRUE); // MC
//  taskSPDdNdEta->SetppAnalysis(kTRUE); //pp analysis
  taskSPDdNdEta->SetTrigger(1);  // 0 = notrigger, 1 = MB1 trigger
  taskSPDdNdEta->SetEvtGen(kTRUE);  //to read Pythia data (kFALSE for Phojet)
  
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  AliAnalysisDataContainer *cout_SPDdNdEta= mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, "SPDdNdEtaData.root");  


   mgr->ConnectInput(taskSPDdNdEta, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskSPDdNdEta, 0, cout_SPDdNdEta);

   // Return task pointer at the end
   return taskSPDdNdEta;
}
