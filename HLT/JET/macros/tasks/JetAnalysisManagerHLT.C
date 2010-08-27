void JetAnalysisManagerHLT() {

  if ( getenv("FASTJET") ) {
    gSystem->Load("libCGAL.so");
    gSystem->Load("libfastjet.so");
    gSystem->Load("libSISConePlugin.so");
  }

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
  
  gSystem->Load("libHLTbase.so");
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTJET.so");
  
  // --------------------------------------------------------------------------------  
    
  if (gApplication) gApplication->InitializeGraphics();

  // --------------------------------------------------------------------------------  
  //
  // Create the chain
  //
  // --------------------------------------------------------------------------------  

  gROOT->LoadMacro("${ALICE_ROOT}/PWG0/CreateESDChain.C");
  TChain* chain = new TChain("esdTree");
  
  chain->Add("~/jet/data/HEAD_2010-01-08/Gen/kPythia6Jets86_104_14TeV/JET-ETA=-0.2,0.2_JET-ET=10,1000_R=0.4_20ev/AliESDs.root");


  //  chain->Add("~/jet/data/HEAD_2009-06-04/kPythia6Jets104_125_14TeV/JET-ETA=-0.2,0.2_JET-ET=10,1000_R=0.7_100ev/AliESDs.root");
  //chain->Add("~/jet/data/HEAD_2009-06-04/kPythia6Jets104_125_14TeV/JET-ETA=-0.2,0.2_JET-ET=10,1000_R=0.7_1ev/AliESDs.root");

  // --------------------------------------------------------------------------------  
  //
  // Create the analysis manager
  //
  // --------------------------------------------------------------------------------  

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
  mgr->SetDebugLevel(10);

  // --------------------------------------------------------------------------------
  //
  // Set Configfiles
  //
  // --------------------------------------------------------------------------------

  // -- HLT FFSC
  AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysisHLT");
  jetana->SetConfigFile("./tasks/ConfigJetAnalysisHLT.C");
  jetana->SetNonStdBranch("jetsHLT");
  jetana->SetDebugLevel(10);
  mgr->AddTask(jetana);

  // --------------------------------------------------------------------------------
  //
  // Create containers for input/output
  //
  // --------------------------------------------------------------------------------  

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(), AliAnalysisManager::kOutputContainer, "histos.root");

  mgr->ConnectInput  (jetana,     0, cinput1  );
  mgr->ConnectOutput (jetana,     0, coutput1 );
  mgr->ConnectOutput (jetana,     1, coutput2 );


  // --------------------------------------------------------------------------------  
  //
  // Run the analysis
  //    
  // --------------------------------------------------------------------------------  

  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}
