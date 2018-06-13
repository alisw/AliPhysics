AliAnalysisTaskStrangenessVsMultiplicityMCRun2 *AddTaskStrangenessVsMultiplicityMCRun2( Bool_t lSaveEventTree = kTRUE, Bool_t lSaveV0 = kTRUE, Bool_t lSaveCascade = kTRUE, TString lExtraOptions = "", const TString lMasterJobSessionFlag = "", TString lExtraOutputName = "")
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
    AliAnalysisTaskStrangenessVsMultiplicityMCRun2 *taskAuxiliary = new AliAnalysisTaskStrangenessVsMultiplicityMCRun2(lSaveEventTree, lSaveV0, lSaveCascade, Form("taskAuxiliary%s",lExtraOutputName.Data()), lExtraOptions);
    mgr->AddTask(taskAuxiliary);
    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    
    outputFileName += ":PWGLF_StrVsMult";
    if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
    outputFileName += lExtraOutputName.Data();
    
    Printf("Set OutputFileName : \n %s\n", outputFileName.Data() );
    
    TString lC[9];
    lC[0] = "cList";
    lC[1] = "cListV0";
    lC[2] = "cListXiMinus";
    lC[3] = "cListXiPlus";
    lC[4] = "cListOmegaMinus";
    lC[5] = "cListOmegaPlus";
    lC[6] = "cTreeEvent";
    lC[7] = "cTreeV0";
    lC[8] = "cTreeCascade";
    
    for(Int_t iname=0;iname<9;iname++)
        lC[iname] += lExtraOutputName.Data();
    
    AliAnalysisDataContainer *coutputLists[6];
    for(Int_t ilist=0;ilist<6;ilist++){
        coutputLists[ilist] = mgr->CreateContainer(lC[ilist].Data(),
                                                   TList::Class(),
                                                   AliAnalysisManager::kOutputContainer,
                                                   outputFileName );
    }
    
    AliAnalysisDataContainer *coutputTree = 0x0;
    AliAnalysisDataContainer *coutputTreeV0 = 0x0;
    AliAnalysisDataContainer *coutputTreeCascade = 0x0;
    
    if( lSaveEventTree ){
        AliAnalysisDataContainer *coutputTree = mgr->CreateContainer(Form("cTreeEvent%s",lExtraOutputName.Data()),
                                                                     TTree::Class(),
                                                                     AliAnalysisManager::kOutputContainer,
                                                                     outputFileName );
        coutputTree->SetSpecialOutput();
    }
    if( lSaveV0 ){
        AliAnalysisDataContainer *coutputTreeV0 = mgr->CreateContainer(Form("cTreeV0%s",lExtraOutputName.Data()),
                                                                       TTree::Class(),
                                                                       AliAnalysisManager::kOutputContainer,
                                                                       outputFileName );
        coutputTreeV0->SetSpecialOutput();
    }
    if (lSaveCascade){
        AliAnalysisDataContainer *coutputTreeCascade = mgr->CreateContainer(Form("cTreeCascade%s",lExtraOutputName.Data()),
                                                                            TTree::Class(),
                                                                            AliAnalysisManager::kOutputContainer,
                                                                            outputFileName );
        coutputTreeCascade->SetSpecialOutput();
    }
    //This one you should merge in file-resident ways...
    
    //Recommendation: Tree as a single output slot
    mgr->ConnectInput (taskAuxiliary, 0, mgr->GetCommonInputContainer());
    for(Int_t ilist=0;ilist<6;ilist++){
        mgr->ConnectOutput(taskAuxiliary, ilist+1, coutputLists[ilist]);
    
    if ( lSaveEventTree ) mgr->ConnectOutput(taskAuxiliary, 7, coutputTree);
    if ( lSaveV0 )        mgr->ConnectOutput(taskAuxiliary, 8, coutputTreeV0);
    if ( lSaveCascade )   mgr->ConnectOutput(taskAuxiliary, 9, coutputTreeCascade);
    
    return taskAuxiliary;
}   
