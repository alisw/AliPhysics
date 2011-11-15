AliAnalysisTaskSE* AddTaskVZEROPbPb(Int_t runNumber)
{
  // Creates a PbPb QA task for VZERO
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskQAsym", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTasQAsym", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
   // Configure analysis
   //===========================================================================
   
 
 
   AliAnaVZEROPbPb* task = new AliAnaVZEROPbPb("AliAnaVZEROPbPb");
   task->SetEquaMultRange(800,800.);
   task->SetSumEquaMultRange(100,15000.,20000.);

   if(runNumber>166532)
   		task->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVLNWU-,CVLN-");
   else
   		task->SetClassesNames("CTRUE-,CVLN-,CVHN-,CPBI1-,CPBI2-,C0HWU-");
   
	mgr->AddTask(task);
  
   AliAnalysisDataContainer *cout  = mgr->CreateContainer("PbPbVZEROHists",TList::Class(),
							  AliAnalysisManager::kOutputContainer, Form("%s:VZERO_PbPb_Performance", 
												     mgr->GetCommonFileName()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 1, cout);

   return task;
   
  
}


