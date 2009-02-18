void AnalysisTrainMuonCAF(char* fileout = "AliAOD.root", char *datasetname = "myDataSet", Int_t nev=1234567890)
{
// Macro to produce a generic AOD starting from an ESD file. 
// The AOD is filled with two tasks: 
// 1- with the first one (AliAnalysisTaskESDfilter), 
//    all the branches of the AOD are filled apart from the muons. 
// 2- with the second task (AliAnalysisTaskESDMuonFilter) 
//    muons tracks are added to the tracks branch 
// This macro works on the CAF
// R. Arnaldi 4/5/08

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
    
  // Reset user processes if CAF if not responding anymore
  // TProof::Reset("lxb6046");

  // Connect to proof
  TProof::Open("lxb6046"); // may be username@lxb6046 if user not the same as on local

  // Clear packages if changing ROOT version on CAF or local
  // gProof->ClearPackages();
  // Enable proof debugging if needed
  // gProof->SetLogLevel(5);

  // Common packages
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  // Analysis-specific
  // --- Enable the PWG3base Package
  gProof->UploadPackage("PWG3muon.par");
  gProof->EnablePackage("PWG3muon");

  // Chain from files staged on CAF
  // gROOT->LoadMacro("CreateESDChain.C");
  // TChain* chain = CreateESDChain("ESD1503X_v1.txt",3);
  // TChain* chain = CreateESDChain("ESD82XX_30Kshort.txt", 10);
  
  // Chain from datasets
  gROOT->LoadMacro("CreateChainFromDataSet.C");
  ds = gProof->GetDataSet(datasetname)->GetStagedSubset();
  chain = CreateChainFromDataSet(ds, "esdTree");   
  
  // Make the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis train");
  
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("FMD CaloCluster");
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName(fileout);
  //aodHandler->SetOutputFileName("AOD.root");

  mgr->SetInputEventHandler(esdHandler);
  mgr->SetOutputEventHandler(aodHandler);
  
  // Set of cuts plugged into the ESD filter
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
  // hard cuts
  AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
  esdTrackCutsH->SetMinNClustersTPC(100);
  esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
  esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsH->SetRequireTPCRefit(kTRUE);
  esdTrackCutsH->SetMinNsigmaToVertex(2);
  esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
  esdTrackCutsH->SetPRange(0.,2.);
  //
  //  muon cuts
  AliESDMuonTrackCuts* esdMuonTrackCuts = new AliESDMuonTrackCuts("AliESDMuonTrackCuts", "test");
  esdMuonTrackCuts->SetPRange(0.,20.);
  //esdMuonTrackCuts->SetPtRange(0.,0.5);   // example of kinematic cuts that can be applied
  
  // track filter (to reject tracks not surviving the cuts - refers to all particles apart from muons)
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsH);
  
  // muon track filter  (to reject muon tracks not surviving the cuts)
  AliAnalysisFilter* trackMuonFilter = new AliAnalysisFilter("trackMuonFilter");
  trackMuonFilter->AddCuts(esdMuonTrackCuts);

  // ESD filter task putting standard info to output generic AOD 
  AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
  //esdfilter->SetTrackFilter(trackFilter);
  esdfilter->SetDebugLevel(10);
  mgr->AddTask(esdfilter);
  
  // ESD filter task putting muon info to output generic AOD 
  AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
  esdmuonfilter->SetTrackFilter(trackMuonFilter);
  mgr->AddTask(esdmuonfilter);

  // Containers for input/output
  AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();
  // Output AOD container. 
  AliAnalysisDataContainer *cout_aod = mgr->GetCommonOutputContainer();
        						    
  // Connect containers to tasks slots
  mgr->ConnectInput  (esdfilter,  0, cin_esd  );
  mgr->ConnectOutput (esdfilter,  0, cout_aod );

  mgr->ConnectInput  (esdmuonfilter,  0, cin_esd);
  mgr->ConnectOutput (esdmuonfilter,  0, cout_aod );

  //
  // Run the analysis
  //	
  if (mgr->InitAnalysis()) {
      mgr->PrintStatus();
      mgr->StartAnalysis("proof",chain,nev);
  }   
}

