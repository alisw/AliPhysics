void runAODFilterJets()
{

  bool bKineFilter = true;
  bool bJets = true;

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  if(bJets)gSystem->Load("libJETAN.so");
  
  // Create the chain
  //
  TChain *chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
  
  ////////////////////////////////////////////////////////////////////////
  // Create the analysis manager
  //
  // Input 
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  // Output
  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAOD.root");
  // MC Truth
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  AliAnalysisManager *mgr  = new AliAnalysisManager("Filter Manager", "Filter Manager");
  if(bKineFilter){
    mgr->SetMCtruthEventHandler(mcHandler);
  }
  
  mgr->SetInputEventHandler  (inpHandler);
  mgr->SetOutputEventHandler (aodHandler);
  
  mgr->SetDebugLevel(10);
  
  // Filtering of MC particles (decays conversions etc)
  // this task is also needed to set the MCEventHandler
  // to the AODHandler, this will not be needed when
  // AODHandler goes to ANALYSISalice
  AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
  if(bKineFilter)mgr->AddTask(kinefilter);
      
      
  // 
  // soft
  AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
  esdTrackCutsL->SetMinNClustersTPC(50);
  esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsL->SetRequireTPCRefit(kTRUE);
  esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
  esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
  esdTrackCutsL->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
  
  //
  // hard
  AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Har
d");
  esdTrackCutsH->SetMinNClustersTPC(100);
  esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
  esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsH->SetRequireTPCRefit(kTRUE);
  esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
  esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
  esdTrackCutsH->SetRequireSigmaToVertex(kFALSE);
  esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
  //

  AliAnalysisTaskJets *jetana = 0;
  AliAnalysisTaskJets *jetanaMC = 0;

  if(bJets){
    jetana = new AliAnalysisTaskJets("JetAnalysis");
    jetana->SetConfigFile("../ConfigJetAnalysis.C");
    jetana->SetDebugLevel(10);

    jetanaMC = new AliAnalysisTaskJets("JetAnalysisMC");
    jetanaMC->SetDebugLevel(10);
    jetanaMC->SetConfigFile("../ConfigJetAnalysisMC.C");
    jetanaMC->SetNonStdBranch("jetsMC");
    mgr->AddTask(jetanaMC);
    mgr->AddTask(jetana);    
  }

  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsL);
  trackFilter->AddCuts(esdTrackCutsH);
  
  AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
  esdfilter->SetTrackFilter(trackFilter);
  
  mgr->AddTask(esdfilter);
  
    
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
							   AliAnalysisManager::kInputContainer);
      
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
      
  coutput1->SetSpecialOutput();
  
  if(bKineFilter){
    mgr->ConnectInput  (kinefilter,     0, cinput1  );
    mgr->ConnectOutput (kinefilter,     0, coutput1 );
  }

  mgr->ConnectInput  (esdfilter,     0, cinput1  );
  mgr->ConnectOutput (esdfilter,     0, coutput1 );
  
  if(bJets){
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),
                                                              AliAnalysisManager::kOutputContainer, "histosJets.root");

    AliAnalysisDataContainer *coutputMC2 = mgr->CreateContainer("histosMC", TList::Class(),                                                       AliAnalysisManager::kOutputContainer, "histosJetsMC.root");

    mgr->ConnectInput  (jetana,     0, cinput1  );
    mgr->ConnectOutput (jetana,     0, coutput1 );
    mgr->ConnectOutput (jetana,     1, coutput2 );
    mgr->ConnectInput  (jetanaMC,     0, cinput1  );
    mgr->ConnectOutput (jetanaMC,     0, coutput1 );
    mgr->ConnectOutput (jetanaMC,     1, coutputMC2 );
  }

  //
  // Run the analysis
  //    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
      
}
