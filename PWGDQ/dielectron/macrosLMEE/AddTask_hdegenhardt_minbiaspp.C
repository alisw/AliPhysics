AliAnalysisTask *AddTask_hdegenhardt_minbiaspp(
          char *period = "16d",
          Bool_t getFromAlien=kFALSE,
          TString cFileName = "Config_hdegenhardt_pp.C",
          Char_t* outputFileName = "LMEE_output.root",
          ULong64_t triggerMask = AliVEvent::kINT7,
          Bool_t rejectPileup = kTRUE,
          Int_t pileuprej = AliDielectronEventCuts::kSPDInMultBins
){

    //=== get the current analysis manager ===========================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTask_hdegenhardt_minbiaspp", "No analysis manager found.");
        return 0;
    }

    TString configFunction(cFileName(0,cFileName.Length() - 2));

    if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFunction.Data()))
        gROOT->LoadMacro(cFileName.Data());

    //Do we have an MC handler?
    printf("------------------------------------------------------\n");
    hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    std::cout << "hasMC = " << hasMC << std::endl;
    if(hasMC) kMixing = 0;

    //=== Create the main dielectron task =============================
    AliAnalysisTaskMultiDielectron *task=new AliAnalysisTaskMultiDielectron("MultiDielectron_mb");
    if (!hasMC) task->UsePhysicsSelection();
    task->SetTriggerMask(triggerMask);
    task->SetTriggerOnV0AND(kTRUE); // only for cross-check

    // SPD pile-up rejection in mult. bins
    if (rejectPileup) {
        task->SetRejectPileup(kTRUE);
        task->SetPileupRejTool(pileuprej);
    }

    // randomize daughters
    task->SetRandomizeDaughters(kTRUE);

    //=== Add event filter ============================================

    task->SetEventFilter(GetEventCutsMinBias() );

	
	printf("Number of die configs = %d\n",nDie);
    //add dielectron analysis with different cuts to the task
    Int_t numberConfigs = nDie;//nDie,1
    for (Int_t i=0; i<numberConfigs; ++i){ //nDie defined in config file
        AliDielectron *diele = Config_hdegenhardt_pp(i, kTRUE, period); // second flag - min.bias analysis?
        if(!diele)continue;
        task->AddDielectron(diele);
    }   
    printf("------------------------------------------------------\n\n");

    mgr->AddTask(task);

    //=== create output containers ===========================

    AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("tree_lowmass",
                         TTree::Class(),
                         AliAnalysisManager::kExchangeContainer,
                         outputFileName);

    AliAnalysisDataContainer *cOutputHist1 =
    mgr->CreateContainer("Output_Histos",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

    AliAnalysisDataContainer *cOutputHist2 =
    mgr->CreateContainer("Output_CF",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

    AliAnalysisDataContainer *cOutputHist3 =
    mgr->CreateContainer("Output_EventStat",
                         TH1D::Class(),
                         AliAnalysisManager::kOutputContainer,
                         outputFileName);

    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 0, coutput1 );
    mgr->ConnectOutput(task, 1, cOutputHist1);
    mgr->ConnectOutput(task, 2, cOutputHist2);
    mgr->ConnectOutput(task, 3, cOutputHist3);

    return task;
}
