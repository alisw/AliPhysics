AliAnalysisTaskStrangenessVsMultiplicityRun2 *AddTaskStrangenessVsMultiplicityRun2( Bool_t lSaveEventTree = kTRUE, Bool_t lSaveV0 = kTRUE, Bool_t lSaveCascade = kTRUE, TString lExtraOptions = "", const TString lMasterJobSessionFlag = "", TString lExtraOutputName = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangenessVsMultiplicity", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangenessVsMultiplicity", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskStrangenessVsMultiplicityRun2 *taskAuxiliary = new AliAnalysisTaskStrangenessVsMultiplicityRun2(lSaveEventTree, lSaveV0, lSaveCascade, Form("taskAuxiliary%s",lExtraOutputName.Data()), lExtraOptions);
    mgr->AddTask(taskAuxiliary);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_StrVsMult";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    outputFileName += lExtraOutputName.Data();
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    TString lC1 = "cList";
    TString lC2 = "cListV0";
    TString lC3 = "cListCascade";
    TString lC4 = "cTreeEvent";
    TString lC5 = "cTreeV0";
    TString lC6 = "cTreeCascade";
    
    lC1 += lExtraOutputName.Data();
    lC2 += lExtraOutputName.Data();
    lC3 += lExtraOutputName.Data();
    lC4 += lExtraOutputName.Data();
    lC5 += lExtraOutputName.Data();
    lC6 += lExtraOutputName.Data();
    
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer(lC1.Data(),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    AliAnalysisDataContainer *coutputListV0 = mgr->CreateContainer(lC2.Data(),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    coutputListV0->SetSpecialOutput();
    AliAnalysisDataContainer *coutputListCascade = mgr->CreateContainer(lC3.Data(),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    coutputListCascade->SetSpecialOutput();
    
    AliAnalysisDataContainer *coutputTree = 0x0;
    AliAnalysisDataContainer *coutputTreeV0 = 0x0;
    AliAnalysisDataContainer *coutputTreeCascade = 0x0;
    if( lSaveEventTree ){
        coutputTree = mgr->CreateContainer(lC4.Data(),
                                                                     TTree::Class(),
                                                                     AliAnalysisManager::kOutputContainer,
                                                                     outputFileName );
        coutputTree->SetSpecialOutput();
    }
    if( lSaveV0 ){
        coutputTreeV0 = mgr->CreateContainer(lC5.Data(),
                                                                       TTree::Class(),
                                                                       AliAnalysisManager::kOutputContainer,
                                                                       outputFileName );
        coutputTreeV0->SetSpecialOutput();
    }
    if (lSaveCascade){
        coutputTreeCascade = mgr->CreateContainer(lC6.Data(),
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            outputFileName );
        coutputTreeCascade->SetSpecialOutput();
    }
    //This one you should merge in file-resident ways...
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskAuxiliary, 1, coutputList);
    mgr->ConnectOutput(taskAuxiliary, 2, coutputListV0);
    mgr->ConnectOutput(taskAuxiliary, 3, coutputListCascade);
    
    if ( lSaveEventTree ) mgr->ConnectOutput(taskAuxiliary, 4, coutputTree);
    if ( lSaveV0 )        mgr->ConnectOutput(taskAuxiliary, 5, coutputTreeV0);
    if ( lSaveCascade )   mgr->ConnectOutput(taskAuxiliary, 6, coutputTreeCascade);
    
    return taskAuxiliary;
}   
