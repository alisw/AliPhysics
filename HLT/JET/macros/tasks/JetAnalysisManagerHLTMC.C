void JetAnalysisManagerHLTMC() {

  if ( getenv("FASTJET") ) {
    gSystem->Load("libCGAL");
    gSystem->Load("libfastjet");
    gSystem->Load("libSISConePlugin");
  }

  gSystem->Load("libTree");
  gSystem->Load("libPhysics");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");
 
  gSystem->Load("libHLTbase");
  gSystem->Load("libAliHLTUtil");
  gSystem->Load("libAliHLTJET");

  Int_t debugLevel = 10;
  
  // --------------------------------------------------------------------------------  
    
  if (gApplication) gApplication->InitializeGraphics();

  // --------------------------------------------------------------------------------  
  //
  // Create the chain
  //
  // --------------------------------------------------------------------------------  

  gROOT->LoadMacro("${ALICE_ROOT}/PWG0/CreateESDChain.C");
  TChain* chain = new TChain("TE");

  chain->Add("/home/jthaeder/jet/data/HEAD_2009-10-26/FastGen/kPythia6Jets125_150_14TeV/JET-ETA=-0.2,0.2_JET-ET=50,1000_R=0.4_200ev/galice.root");

  // --------------------------------------------------------------------------------  
  //
  // Create the analysis manager
  //
  // --------------------------------------------------------------------------------  

  // Output
  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("aod.root");

  // MC Truth
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mcHandler->SetReadTR(kFALSE); // for fastgen

  AliAnalysisManager *mgr  = new AliAnalysisManager("Jet Manager", "Jet Manager");
  mgr->SetOutputEventHandler (aodHandler);
  mgr->SetMCtruthEventHandler(mcHandler);
  mgr->SetDebugLevel(debugLevel);

  // --------------------------------------------------------------------------------
  //
  // Set Configfiles
  //
  // --------------------------------------------------------------------------------

  // -- HLT FFSC - MC
#if FFSC
  AliAnalysisTaskJets *taskFFSC = new AliAnalysisTaskJets("JetAnalysisHLTMC_FFSC");
  taskFFSC->SetConfigFile("./tasks/ConfigJetAnalysisHLTMC.C");
  taskFFSC->SetNonStdBranch("jetsHLTMC_FFSC");
  taskFFSC->SetDebugLevel(debugLevel);
  mgr->AddTask(taskFFSC);
#endif

  // -- HLT FastJet - Kt
#if FASTJET_KT
  if ( getenv("FASTJET") ) {
    AliAnalysisTaskJets *taskKt = new AliAnalysisTaskJets("JetAnalysisHLTMC_Kt");
    taskKt->SetConfigFile("./tasks/ConfigJetAnalysisHLTMCKt.C");
    taskKt->SetNonStdBranch("jetsHLTMC_Kt");
    taskKt->SetDebugLevel(debugLevel);
    mgr->AddTask(taskKt);
  }
#endif

#if FASTJET_ANTIKT
  // -- HLT FastJet - AntiKt
  if ( getenv("FASTJET") ) {
    AliAnalysisTaskJets *taskAntiKt = new AliAnalysisTaskJets("JetAnalysisHLTMC_AntiKt");
    taskAntiKt->SetConfigFile("./tasks/ConfigJetAnalysisHLTMCAntiKt.C");
    taskAntiKt->SetNonStdBranch("jetsHLTMC_AntiKt");
    taskAntiKt->SetDebugLevel(debugLevel);
    mgr->AddTask(taskAntiKt);
  }
#endif

  // --------------------------------------------------------------------------------
  //
  // Create containers for input/output
  //
  // --------------------------------------------------------------------------------  

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();

#if FFSC
  AliAnalysisDataContainer *coutputFFSC = mgr->CreateContainer("histos", TList::Class(), 
							       AliAnalysisManager::kOutputContainer, "histos_MC_FFSC.root");
  mgr->ConnectInput  (taskFFSC,     0, cinput1  );
  mgr->ConnectOutput (taskFFSC,     0, coutput1 );
  mgr->ConnectOutput (taskFFSC,     1, coutputFFSC );
#endif
  
#if FASTJET_KT
  if ( getenv("FASTJET") ) {
    AliAnalysisDataContainer *coutputKt = mgr->CreateContainer("histos", TList::Class(), 
							       AliAnalysisManager::kOutputContainer, "histos_MC_KT.root");
    mgr->ConnectInput  (taskKt,     0, cinput1  );
    mgr->ConnectOutput (taskKt,     0, coutput1 );
    mgr->ConnectOutput (taskKt,     1, coutputKt );
  }
#endif

#if FASTJET_ANTIKT
  if ( getenv("FASTJET") ) {
    AliAnalysisDataContainer *coutputAntiKt = mgr->CreateContainer("histos", TList::Class(), 
								   AliAnalysisManager::kOutputContainer, "histos_MC_ANTIKT.root");
    mgr->ConnectInput  (taskAntiKt,     0, cinput1  );
    mgr->ConnectOutput (taskAntiKt,     0, coutput1 );
    mgr->ConnectOutput (taskAntiKt,     1, coutputAntiKt );
  }
#endif

  //  AliAnalysisTaskKineFilter *kinefilter = new AliAnalysisTaskKineFilter("Kine Filter");
  //  mgr->AddTask(kinefilter);
  //  mgr->ConnectInput  (kinefilter,     0, cinput1  );
  //  mgr->ConnectOutput (kinefilter,     0, coutput1 );

  // --------------------------------------------------------------------------------  
  //
  // Run the analysis
  //    
  // --------------------------------------------------------------------------------  

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}
