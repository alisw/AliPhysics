const char* fName ="NetLambdaMCResults.root";

AliAnalysisTaskNetLambdaMCTrad *AddTaskNetLambdaMCTrad( Bool_t lSaveEventTree = kTRUE, Bool_t lSaveV0 = kTRUE, TString lExtraOptions = "", const TString lMasterJobSessionFlag = "")
{
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskNetLambdaMCTrad", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskNetLambdaMCTrad", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    
    // Create the task and configure it.
    //===========================================================================
    AliAnalysisTaskNetLambdaMCTrad* task = new  AliAnalysisTaskNetLambdaMCTrad(lSaveEventTree, lSaveV0,  "task", lExtraOptions);
    mgr->AddTask(task);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
  
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("cList",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 fName );
   
//    if( lSaveEventTree ){
//        AliAnalysisDataContainer *coutputTree = mgr->CreateContainer("cTreeEvent",
//                                                                     TTree::Class(),
//                                                                     AliAnalysisManager::kOutputContainer,
//                                                                     fName );
//        coutputTree->SetSpecialOutput();
//    }
//    if( lSaveV0 ){
//        AliAnalysisDataContainer *coutputTreeV0 = mgr->CreateContainer("cTreeV0",
//                                                                       TTree::Class(),
//                                                                       AliAnalysisManager::kOutputContainer,
//                                                                       fName );
//        coutputTreeV0->SetSpecialOutput();
//    }
    
    
    mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, coutputList);
//    mgr->ConnectOutput(task, 2, coutputTree);
//    mgr->ConnectOutput(task, 3, coutputTreeV0);
    
    return task;

   
}
