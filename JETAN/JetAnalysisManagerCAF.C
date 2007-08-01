void JetAnalysisManagerCAF()
{
    //
    if (gApplication) gApplication->InitializeGraphics();
    gROOT->LoadMacro("CreateESDChain.C");
    //
    // Connect to proof
    
    TProof::Reset("proof://morsch@lxb6046.cern.ch"); 
    TProof::Open("proof://morsch@lxb6046.cern.ch"); 
    
    //   gProof->SetParallel(1);
//    gProof->ClearPackage("ESD");
//    gProof->ClearPackage("AOD");
    //gProof->ClearPackage("JETAN");
    //gProof->ClearPackage("ANALYSIS");
    
    gProof->ShowEnabledPackages();
    // Enable the ESD Package
    gProof->UploadPackage("ESD.par");
    gProof->EnablePackage("ESD");
     // Enable the AOD Package
    gProof->UploadPackage("AOD.par");
    gProof->EnablePackage("AOD");
     // Enable the Analysis Package
    gProof->UploadPackage("ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");

    // Enable the JETAN Package
    gProof->UploadPackage("JETAN.par");
    gProof->EnablePackage("JETAN");

    //
    gProof->ShowEnabledPackages();
    
    //
    //
    // Create the chain
    //
    TChain* chain = CreateESDChain("test.txt", 200);
    
    //
    // Create the analysis manager
    //
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetEventHandler(aodHandler);
    mgr-> SetDebugLevel(10);

//
//   Jet Finder Task
//

    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);
//
//  ESD Filter Task
//
    //
    // Set of cuts
    // 
    // standard
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
    //
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsH);
    //
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");

    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TH1F::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");

    mgr->ConnectInput (jetana, 0, cinput1);
    mgr->ConnectOutput(jetana, 0, coutput1);
    mgr->ConnectOutput(jetana, 1, coutput2);

    mgr->ConnectInput  (esdfilter,  0, cinput1 );
    mgr->ConnectOutput (esdfilter,  0, coutput1);
    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("proof",chain);
}
