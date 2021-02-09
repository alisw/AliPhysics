AliAnalysisTask *AddTask_hmurakam_minbiaspp(Bool_t getFromAlien=kFALSE,
					    TString year ="16",
					    Bool_t hasSpline =kFALSE,
					    TString cFileName = "Config_hmurakam_pp.C",
					    Char_t* outputFileName="LMEE.root",
					    ULong64_t triggerMask = AliVEvent::kINT7,
					    Bool_t rejectPileup = kTRUE,
					    Int_t pileuprej = AliDielectronEventCuts::kSPDInMultBins
                                            )
{
    
    //=== get the current analysis manager ===========================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTask_hmurakam_minbiaspp", "No analysis manager found.");
        return 0;
    }
    
    //Base Directory for GRID / LEGO Train
    TString configBasePath= "$ALICE_PHYSICS/PWGDQ/dielectron/macrosLMEE/";
    if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/macrosLMEE/%s .",cFileName.Data()))) ){
        configBasePath=Form("%s/",gSystem->pwd());
    }
    
    TString configFilePath(configBasePath+cFileName);
    
    std::cout << "Configpath: " << configFilePath << std::endl;
    
    TString configFunction(cFileName(0,cFileName.Length() - 2));

    if (!gROOT->GetListOfGlobalFunctions()->FindObject(configFunction.Data()))
        gROOT->LoadMacro(configFilePath.Data());
    
    //Do we have an MC handler?
    Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
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

    task->SetEventFilter(GetEventCutsMinBias());
    
    //add dielectron analysis with different cuts to the task
    for (Int_t i=0; i<nDie; ++i){ //nDie defined in config file
        
        AliDielectron *diele = Config_hmurakam_pp(i,year.Data(),hasSpline);
        if(!diele)continue;
        task->AddDielectron(diele);
    }
    
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

