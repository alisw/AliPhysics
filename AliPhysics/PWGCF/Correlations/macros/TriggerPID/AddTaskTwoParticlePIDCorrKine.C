
AliTwoParticlePIDCorrKine *AddTaskTwoParticlePIDCorrKine(){
    
 
 AliTwoParticlePIDCorrKine* kine = new  AliTwoParticlePIDCorrKine("");
 
 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
        cout<<"AliAnalysisKine", "No analysis manager to connect to."<<endl;
        return NULL;
  }

// Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddKineBF", "This Kine requires an input event handler");
    return NULL;
  }

TString type = mgr->GetInputEventHandler()->GetDataType();
    
  mgr->AddTask(kine);
    
    
AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

    // connect input/output                                                                                                                                             
mgr->ConnectInput(kine, 0, cinput);
mgr->ConnectOutput(kine, 1, coutput);

return kine;


}




