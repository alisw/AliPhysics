void JetAnalysisManagerHLTKine()
{
      gSystem->Load("libCGAL.so");
      gSystem->Load("${FASTJET}/lib/libfastjet.so");

      gSystem->Load("libTree.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libANALYSISalice.so");
      gSystem->Load("libJETAN.so");

     //
    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //

    gROOT->LoadMacro("${ALICE_ROOT}/PWG0/CreateESDChain.C");
    TChain* chain = new TChain("TE");
    chain->Add("/home/jthaeder/jet/data/HEAD_2009-03-17/FastGen/kPythia6Jets104_125_14TeV/JET-ETA=-0.2,0.2_JET-ET=50,1000_R=0.4_10ev/galice.root");

    /////////////////////////////////////////////////////////////////////////////////// 
    // Create the analysis manager
    //
    // Output
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    // MC Truth
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);

    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
    mgr->SetOutputEventHandler (aodHandler);
    mgr->SetMCtruthEventHandler(mcHandler);
    mgr->SetDebugLevel(10);

    AliAnalysisTaskKineFilter *kinefilter = new AliAnalysisTaskKineFilter("Kine Filter");
    mgr->AddTask(kinefilter);
    
    
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    //jetana->SetConfigFile("ConfigJetAnalysisFastjetMC.C");
    jetana->SetConfigFile("${ALICE_ROOT}/HLT/JET/macros/tasks/ConfigJetAnalysisHLTMC.C");
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

    mgr->ConnectInput  (kinefilter,     0, cinput1  );
    mgr->ConnectOutput (kinefilter,     0, coutput1 );

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
