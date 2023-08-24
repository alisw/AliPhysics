//_____________________________________________________________________
 AliAnalysisTask *AddTaskJCorrectionMap(TString taskName = "JCorrectionMapTask"){

     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     //==== Check the analysis type using the event handlers connected to the analysis mgr	
     if (!mgr->GetInputEventHandler() ){
 	    ::Error("AddTaskFFluc", "This task requires an input event handler" );
 	    return NULL;
     }

     //==== JCORRAN TASK
     AliJCorrectionMapTask *FFtask = new AliJCorrectionMapTask(taskName.Data());
     mgr->AddTask((AliAnalysisTask*) FFtask);

     //==== Create containers for input/output
     AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
     //AliAnalysisDataContainer *FFhist = mgr->CreateContainer(Form("%scontainer",FFtask->GetName()),  TList::Class(), AliAnalysisManager::kOutputContainer,  Form("%s:%s",  AliAnalysisManager::GetCommonFileName(), FFtask->GetName()));
     //==== Connect input/output
     mgr->ConnectInput(FFtask, 0, cinput);
     //mgr->ConnectOutput(FFtask, 1, FFhist);

     return FFtask;
}
