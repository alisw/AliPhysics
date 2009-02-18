{
    gSystem->Load("libPhysics.so");
    // Connecting to the PROOF cluster
    TProof::Open("alicecaf");

    // Clear packages if changing ROOT version on CAF or local
    //gProof->ClearPackages();
    
    // Enable the needed packages
    gProof->UploadPackage("STEERBase");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("ESD");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("AOD");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("ANALYSIS");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("ANALYSISalice");
    gProof->EnablePackage("ANALYSISalice");
    gProof->UploadPackage("JETAN");
    gProof->EnablePackage("JETAN");
       
    // Create the analysis manager
    mgr = new AliAnalysisManager("Analysis UE test");

    // Create, add task
    gProof->Load("AliKineTrackCuts.cxx+");
    AliKineTrackCuts* trackCuts = new AliKineTrackCuts("AliKineTrackCuts", "Eta");
    trackCuts->SetEtaRange(-1., 1.);
    
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(trackCuts);
    
    gProof->Load("AliAnalysisTaskKineFilter.cxx+");
    AliAnalysisTaskKineFilter *kinefilter = new AliAnalysisTaskKineFilter("Kine Filter");
    kinefilter->SetTrackFilter(trackFilter);
    mgr->AddTask(kinefilter);
    
    // Create chain of input files
    gROOT->LoadMacro("CreateESDChain.C");
    chain = CreateChain( "TE", "KINE82XX_30K.txt", 200);
 
 
 
    /////////////////////////////////////////////////////////////////////////////////// 
    // Create the analysis manager
    //
    // Input 
    // MC Truth
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);
    
    // Output
    AliAODHandler* aodHandler = new AliAODHandler();
    aodHandler->SetOutputFileName( "aodKine.root" );

    AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
//    mgr->SetInputEventHandler  (inpHandler);
    mgr->SetMCtruthEventHandler(mcHandler);
    mgr->SetOutputEventHandler (aodHandler);
    mgr->SetDebugLevel(10);
    
    /////////////////////////////////////////////////////////////////////////////////// 
    
    //
    // Set of cuts
    // 
    AliKineTrackCuts* trackCuts = new AliKineTrackCuts("AliKineTrackCuts", "Eta");
    trackCuts->SetEtaRange(-1., 1.);
 //   trackCuts->SetPtRange(0.5);
              
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(trackCuts);
    
    
    AliAnalysisTaskKineFilter* kinefilter = new AliAnalysisTaskKineFilter("Kine Filter");
    kinefilter->SetTrackFilter(trackFilter);
    kinefilter->SetDebugLevel(10);
    mgr->AddTask(kinefilter);
    
    //    Analysis Task for Jet
    //  AliAnalysisTaskJets need ConfigJetAnalysis.C macro !!!!
    AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetDebugLevel(10);
    mgr->AddTask(jetana);

    //
    // Create containers for input/output                  
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

    AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer, "histos.root");

    mgr->ConnectInput  (kinefilter,  0, cinput1  );
    mgr->ConnectOutput (kinefilter,  0, coutput1 );

    mgr->ConnectInput  (jetana,     0, cinput1  );
    mgr->ConnectOutput (jetana,     0, coutput1 );
    mgr->ConnectOutput (jetana,     1, coutput2 );


    //
    // Run the analysis
    //    
    if( mgr->InitAnalysis() ) {
      mgr->PrintStatus();
      mgr->StartAnalysis("proof", chain );
    }
}
