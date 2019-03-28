
AliAnalysisTaskStrangeCascadesDiscrete *AddTaskStrangeCascadesDiscrete( Bool_t lSaveEventTree = kTRUE,
                                                                                                     Bool_t lSaveV0 = kTRUE,
                                                                                                     Bool_t lSaveCascade = kTRUE,
                                                                                                     TString lExtraOptions = "",
                                                                                                     const TString lMasterJobSessionFlag = "",
                                                                                                     TString lExtraOutputName = "")
{
    
    
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangeCascadesDiscrete", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangeCascadesDiscrete", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskStrangeCascadesDiscrete *taskAuxiliary =
    new AliAnalysisTaskStrangeCascadesDiscrete(lSaveEventTree,
                                                         lSaveV0,
                                                         
                                                         lSaveCascade,
                                                         "taskAuxiliary",
                                                         
                                                         lExtraOptions);  // const*charname =  Form("taskAuxiliary%s",lExtraOutputName.Data())
    
    mgr->AddTask(taskAuxiliary);
    
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    //output Lists
  /*  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("skorodum_list", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cList.root");
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("skorodum_listK", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListK0Short.root");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("skorodum_listLambda", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListLambda.root");
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("skorodum_listAntiLambda", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListAntiLambda.root");
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("skorodum_listXiMinus", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListXiMinus.root");
    AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("skorodum_listXiPlus", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListXiPlus.root");
    
    AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("skorodum_listomegaminus", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListOmegaMinus.root");
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("skorodum_listomegaplus", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cListOmegaPlus.root");*/
    //output trees
    /*  AliAnalysisDataContainer *coutput8 = 0x0;
     AliAnalysisDataContainer *coutput9 = 0x0;
     AliAnalysisDataContainer *coutput10 = 0x0; */
    
    
    
    // if(lSaveEventTree)
    AliAnalysisDataContainer*  coutput8 = mgr->CreateContainer("skorodum_treeevent", TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,"skorodum_cTreeEvent.root");
    // if(lSaveV0)
    AliAnalysisDataContainer* coutput9 = mgr->CreateContainer("skorodum_treevo", TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,"skorodum_cTreeV0.root");
    
    //  if(lSaveCascade)
    AliAnalysisDataContainer* coutput10 = mgr->CreateContainer("skorodum_treecascade", TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,"skorodum_cTreeCascade.root");
    
    
   // AliAnalysisDataContainer* coutput11 = mgr->CreateContainer("skorodum_numbereventstree", TTree::Class(),AliAnalysisManager::kOutputContainer,"skorodum_cNumberEvents.root");
   // AliAnalysisDataContainer* coutput12 = mgr->CreateContainer("skorodum_treevocascade", TTree::Class(),AliAnalysisManager::kOutputContainer,"skorodum_cV0CascadeTree.root"); //not needed
    
    
    mgr->ConnectInput (taskAuxiliary, 0, cinput);
  /*  mgr->ConnectOutput(taskAuxiliary, 1, coutput0);
    mgr->ConnectOutput(taskAuxiliary, 2, coutput1);
    mgr->ConnectOutput(taskAuxiliary, 3, coutput2);
    mgr->ConnectOutput(taskAuxiliary, 4, coutput3);
    mgr->ConnectOutput(taskAuxiliary, 5, coutput4);
    mgr->ConnectOutput(taskAuxiliary, 6, coutput5);
    mgr->ConnectOutput(taskAuxiliary, 7, coutput6);
    mgr->ConnectOutput(taskAuxiliary, 8, coutput7);*/
    //  if(lSaveEventTree)
    mgr->ConnectOutput(taskAuxiliary, 1, coutput8); //event
    
    //  if(lSaveV0)
    mgr->ConnectOutput(taskAuxiliary, 2, coutput9); //V0
    
    //   if(lSaveCascade)
    mgr->ConnectOutput(taskAuxiliary, 3, coutput10); //cascades
   // mgr->ConnectOutput(taskAuxiliary, 4, coutput11); // number of events, same and less than coutput8
   // mgr->ConnectOutput(taskAuxiliary, 13, coutput12); //combined tree, not needed

    
    
    
    return taskAuxiliary;
    
    
    
    
}


