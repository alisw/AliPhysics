void JetAnalysisManagerLoc()
{
      gSystem->Load("libTree.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libJETAN.so");
      gSystem->Load("libPWG0base.so");
     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //
    gROOT->LoadMacro("CreateESDChain.C");
    TChain* chain = new TChain("esdTree");
    chain->Add("/home/morsch/AliRoot/data/data_jets102/AliESDs.root");

    /////////////////////////////////////////////////////////////////////////////////// 
    // Create the analysis manager
    //
    // Input 
    AliESDInputHandler* inpHandler = new AliESDInputHandler();
    // Output
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    // MC Truth
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    

    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetInputEventHandler  (inpHandler);
    mgr->SetOutputEventHandler (aodHandler);
    mgr->SetMCtruthEventHandler(mcHandler);
    mgr-> SetDebugLevel(10);
    /////////////////////////////////////////////////////////////////////////////////// 
    
    
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
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsH);
    //
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);

    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);

    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");

    mgr->ConnectInput  (esdfilter,  0, cinput1  );
    mgr->ConnectOutput (esdfilter,  0, coutput1 );

    mgr->ConnectInput  (jetana,     0, cinput1  );
    mgr->ConnectOutput (jetana,     0, coutput1 );
    mgr->ConnectOutput (jetana,     1, coutput2 );


    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
}
