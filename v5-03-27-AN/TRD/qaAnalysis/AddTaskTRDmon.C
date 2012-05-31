//////////////////////////////////////////////////////////////
// Macro to setup AliAnalysisTaskTRDmon                     //
// for TRD monitoring.                                      //
// ESD handler must be attached to the AliAnalysisManager   //
//                                                          //
// Output:                                                  //
//  TRDmon.root containing a TObjArray of histograms        //
//                                                          //
// 25.02.2010 Ionut Arsene i.c.arsene@gsi.de                //
//////////////////////////////////////////////////////////////

AliAnalysisTaskTRDmon* AddTaskTRDmon(const Char_t* triggerName = "",
				     Bool_t isCollisionTrigger = kTRUE) 
{
  //
  // Configures an AliAnalysisTRDmon task and adds it to the analysis train
  //

  // check the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    Error("AddTaskTRDmon", "The analysis manager is not initialized");
    return 0x0;
  }

  // check the ESD input handler
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if(!type.Contains("ESD")) {
    Error("AddTaskTRDmon", "AliAnalysisTaskTRDmon task needs the manager to have an ESD input handler.");
    return 0x0;
  }

  // configure task
  AliAnalysisTaskTRDmon *task = new AliAnalysisTaskTRDmon("TRDmon");
  task->SetTriggerName(triggerName.Data());
  task->SetIsCollisionEvent(isCollisionTrigger);
  mgr->AddTask(task);

  // connect input
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // connect output
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("TRDmon", TObjArray::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 Form("%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, output);

  return task;
}
