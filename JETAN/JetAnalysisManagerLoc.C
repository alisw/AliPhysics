void JetAnalysisManagerLoc()
{
      gSystem->Load("libTree.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libJETAN.so");
     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //
    TChain* chain = new TChain("esdTree");
    chain->Add("/home/morsch/analysis/AliEn/Interactive/esd/001/AliESDs.root");
    chain->Add("/home/morsch/analysis/AliEn/Interactive/esd/002/AliESDs.root");
    //
    // Make the analysis manager
    //
    AliAnalysisManager *mgr     = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr-> SetDebugLevel(10);
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "jets.root");

    mgr->ConnectInput (jetana, 0, cinput1);
    mgr->ConnectOutput(jetana, 0, coutput1);
    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
}
