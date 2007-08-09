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
     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //
    TChain* chain = new TChain("esdTree");
    chain->Add("/home/morsch/AliRoot/data/data_jets102/AliESDs.root");
    chain->Add("/home/morsch/AliRoot/data/data_jets103/AliESDs.root");
    //
    // Create the analysis manager
    //
    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetDebugLevel(10);   
    //
    // Common output service
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    mgr->SetOutputEventHandler(aodHandler);
    //
    // Common MC truth services
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);

    //
    // Jet analysis
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);

    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "aod.root");

    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TH1F::Class(),
							      AliAnalysisManager::kOutputContainer, "histos.root");

    mgr->ConnectInput  (jetana,  0, cinput1  );
    mgr->ConnectOutput (jetana,  0, coutput1 );
    mgr->ConnectOutput (jetana,  1, coutput2 );
    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
}
