void JetAnalysisManagerCAF()
{
    //
    if (gApplication) gApplication->InitializeGraphics();
    gROOT->LoadMacro("CreateESDChain.C");
    //
    // Connect to proof

    TProof::Reset("proof://morsch@lxb6046.cern.ch"); 
    TProof::Open("proof://morsch@lxb6046.cern.ch"); 

    //gProof->SetParallel(3);
    gProof->ClearPackage("ESD");
    gProof->ClearPackage("JETAN");
    gProof->ClearPackage("ANALYSIS");
    
    gProof->ShowEnabledPackages();
    // Enable the ESD Package
    gProof->UploadPackage("ESD.par");
    gProof->EnablePackage("ESD");
     // Enable the Analysis Package
    gProof->UploadPackage("ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
    // Enable the JETAN Package
    gProof->UploadPackage("JETAN.par");
    gProof->EnablePackage("JETAN");
    // Load Configuration macro
    //gProof->Load("ConfigJetAnalysis.C");
    //
    gProof->ShowEnabledPackages();
    
    //
    //
    // Create the chain
    //
    TChain* chain = CreateESDChain("test.txt", 200);
    //
    // Make the analysis manager
    //
    AliAnalysisManager *mgr     = new AliAnalysisManager("Jet Manager", "Jet Manager");
    AliAODHandler* aodHandler   = new AliAODHandler();
    mgr->SetEventHandler(aodHandler);
    mgr-> SetDebugLevel(10);
    
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);
    //
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "aod.root");

    mgr->ConnectInput (jetana, 0, cinput1);
    mgr->ConnectOutput(jetana, 0, coutput1);
    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("proof",chain);
}
