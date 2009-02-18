void runESDQA_CAF(Char_t *dataset = "/PWG0/COMMON/LHC08c11_10TeV_0.5T", Int_t max_evt=200000)
{
    TProof::Open("lxb6046");

    /*
    gProof->UploadPackage("AF-v4-14.par");
    gProof->EnablePackage("AF-v4-14");
    */
    //    gProof->ClearPackages();
    gProof->UploadPackage("${ALICE_ROOT}/STEERBase.par");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("${ALICE_ROOT}/ESD.par");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("${ALICE_ROOT}/AOD.par");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("${ALICE_ROOT}/ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("${ALICE_ROOT}/ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");
    gProof->UploadPackage("${ALICE_ROOT}/JETAN.par");
    gProof->EnablePackage("JETAN");

    gProof->UploadPackage("${ALICE_ROOT}/PWG4JetTasks.par");
    gProof->EnablePackage("PWG4JetTasks");
    
    // Make the analysis manager

    //
    // Chain from CAF


    AliAnalysisManager *mgr  = new AliAnalysisManager("ESDAODHists");

    mgr->SetDebugLevel(2);
    //AliLog::EnableDebug(kTRUE);
    //AliLog::SetGlobalLogLevel(2);


    // Set of cuts
    // 
    // standard
/*
    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMinNsigmaToVertex(3);
    esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
    //
    // hard
    AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
    esdTrackCutsH->SetMinNClustersTPC(100);
    esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
    esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsH->SetRequireTPCRefit(kTRUE);
    esdTrackCutsH->SetMinNsigmaToVertex(2);
    esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
    //
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsH);
    //
*/  

    // Take it from the library no need to compile directly
    // Standalone does not need ANALYISalice/AOD/JETAN/PWG4JetsTasks
    //    gProof->Load("AliAnaESDSpectraQA.cxx++g");

    AliAnaESDSpectraQA *pwg4QA = new AliAnaESDSpectraQA("ESDSpectraQA");
    mgr->AddTask(pwg4QA);

    AliESDInputHandler *esdHandler = new AliESDInputHandler;
    //esdHandler->SetInactiveBranches("FMD CaloCluster");
    mgr->SetInputEventHandler(esdHandler);

    AliMCEventHandler* mcH = new AliMCEventHandler;
    mgr->SetMCtruthEventHandler(mcH);


    // Create containers for input/output
    // Top ESD container
    AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();

    //AliAnalysisDataContainer *cout_aodex = mgr->CreateContainer("cAodEx", TTree::Class(),
//							      AliAnalysisManager::kExchangeContainer, "default");

    // Histos
    AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("qa_hists", TList::Class(),
    //AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("qa_hists", TObjArray::Class(),
							      AliAnalysisManager::kOutputContainer, "PWG4QAHists.root");

    mgr->ConnectInput (pwg4QA,  0, cin_esd );
    mgr->ConnectOutput (pwg4QA,  0, cout_hist );
   //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    //mgr->StartAnalysis("local",chain);
    mgr->StartAnalysis("proof", dataset, max_evt);
    //delete mgr;
}
