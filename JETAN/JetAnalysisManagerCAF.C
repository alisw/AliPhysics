void JetAnalysisManagerCAF()
{
    //
    if (gApplication) gApplication->InitializeGraphics();
    //
    // Connect to proof
    
//    TProof::Reset("proof://morsch@lxb6046.cern.ch"); 
    TProof::Open("proof://morsch@lxb6046.cern.ch");
//    gProof->SetParallel(1);
//    gProof->ClearPackage("STEERBase");
//    gProof->ClearPackage("ESD");
//    gProof->ClearPackage("AOD");
//    gProof->ClearPackage("JETAN");
//    gProof->ClearPackage("ANALYSIS");
//    gProof->ClearPackage("ANALYSISalice");
    
    gProof->ShowEnabledPackages();

    // Enable the STEERBase Package
    gProof->UploadPackage("STEERBase.par");
    gProof->EnablePackage("STEERBase");
    // Enable the ESD Package
    gProof->UploadPackage("ESD.par");
    gProof->EnablePackage("ESD");
     // Enable the AOD Package
    gProof->UploadPackage("AOD.par");
    gProof->EnablePackage("AOD");
     // Enable the Analysis Package
    gProof->UploadPackage("ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");

    gProof->UploadPackage("ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");

    // Enable the JETAN Package
    gProof->UploadPackage("JETAN.par");
    gProof->EnablePackage("JETAN");

    //
    gProof->ShowEnabledPackages();
    
    //
    //
    // Create the chain
    //
    // TChain* chain = CreateESDChain("test.txt", 200);
 // Input 
    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    //
    // Create the analysis manager
    //
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("jets.root");
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetOutputEventHandler(aodHandler);
    mgr->SetInputEventHandler(inpHandler);
    mgr-> SetDebugLevel(10);


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
    //
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    //
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);

//
//   Jet Finder Task
//

    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    jetana->SetConfigFile("ConfigJetAnalysisAOD.C");
    mgr->AddTask(jetana);
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
    coutput1->SetSpecialOutput();
    
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");

    mgr->ConnectInput  (esdfilter,  0, cinput1 );
    mgr->ConnectOutput (esdfilter,  0, coutput1);

    mgr->ConnectInput (jetana, 0, cinput1);
    mgr->ConnectOutput(jetana, 0, coutput1);
    mgr->ConnectOutput(jetana, 1, coutput2);


    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("proof","/PWG4/arian/jetjetAbove_50_real");
}
