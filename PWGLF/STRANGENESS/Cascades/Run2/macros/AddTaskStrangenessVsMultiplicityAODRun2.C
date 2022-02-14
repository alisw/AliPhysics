AliAnalysisTaskStrangenessVsMultiplicityAODRun2 *AddTaskStrangenessVsMultiplicityAODRun2( Bool_t lSaveEventTree = kTRUE, Bool_t lSaveV0 = kTRUE, Bool_t lSaveCascade = kTRUE, TString lExtraOptions = "", const TString lMasterJobSessionFlag = "", TString lExtraOutputName = "")
{
    // Creates, configures and attaches to the train a cascades check task.
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskStrangenessVsMultiplicityAODRun2", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskStrangenessVsMultiplicityAODRun2", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Create and configure the task
    AliAnalysisTaskStrangenessVsMultiplicityAODRun2 *taskAuxiliary = new AliAnalysisTaskStrangenessVsMultiplicityAODRun2(lSaveEventTree, lSaveV0, lSaveCascade, Form("taskAuxiliary%s",lExtraOutputName.Data()), lExtraOptions);
    mgr->AddTask(taskAuxiliary);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_StrVsMult";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    outputFileName += lExtraOutputName.Data();
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    TString lC[11];
    lC[0] = "cList";
    lC[1] = "cListK0Short";
    lC[2] = "cListLambda";
    lC[3] = "cListAntiLambda";
    lC[4] = "cListXiMinus";
    lC[5] = "cListXiPlus";
    lC[6] = "cListOmegaMinus";
    lC[7] = "cListOmegaPlus";
    lC[8] = "cTreeEvent";
    lC[9] = "cTreeV0";
    lC[10] = "cTreeCascade";
    
    for(Int_t iname=0;iname<11;iname++)
    lC[iname] += lExtraOutputName.Data();
    
    AliAnalysisDataContainer *coutputLists[8];
    for(Int_t ilist=0;ilist<8;ilist++){
        coutputLists[ilist] = mgr->CreateContainer(lC[ilist].Data(),
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileName );
    }
    
    AliAnalysisDataContainer *coutputTree = 0x0;
    AliAnalysisDataContainer *coutputTreeV0 = 0x0;
    AliAnalysisDataContainer *coutputTreeCascade = 0x0;
    if( lSaveEventTree ){
        coutputTree = mgr->CreateContainer(lC[8].Data(),
                                                                     TTree::Class(),
                                                                     AliAnalysisManager::kOutputContainer,
                                                                     outputFileName );
        coutputTree->SetSpecialOutput();
    }
    if( lSaveV0 ){
        coutputTreeV0 = mgr->CreateContainer(lC[9].Data(),
                                                                       TTree::Class(),
                                                                       AliAnalysisManager::kOutputContainer,
                                                                       outputFileName );
        coutputTreeV0->SetSpecialOutput();
    }
    if (lSaveCascade){
        coutputTreeCascade = mgr->CreateContainer(lC[10].Data(),
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            outputFileName );
        coutputTreeCascade->SetSpecialOutput();
    }
    //This one you should merge in file-resident ways...
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
    for(Int_t ilist=0;ilist<8;ilist++){
        mgr->ConnectOutput(taskAuxiliary, ilist+1, coutputLists[ilist]);
    }
    if ( lSaveEventTree ) mgr->ConnectOutput(taskAuxiliary, 9, coutputTree);
    if ( lSaveV0 )        mgr->ConnectOutput(taskAuxiliary, 10, coutputTreeV0);
    if ( lSaveCascade )   mgr->ConnectOutput(taskAuxiliary, 11, coutputTreeCascade);
    
    return taskAuxiliary;
}   
