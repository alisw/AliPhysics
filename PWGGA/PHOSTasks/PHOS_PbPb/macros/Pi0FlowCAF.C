void Pi0FlowCAF(const char* dataset="/alice/data/LHC11h_2_*AOD095",
		const Int_t nEvents=-1, 
		const Int_t nEventsSkip=0)
{
  /* $Id$ */

  // Running analysis AliAnalysisTaskPi0Flow in CAF on AOD datasets

  // Enabling authentication with grid certificate
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");

  const TString proofCluster="alice-caf.cern.ch";
  const TString alirootVersion = "VO_ALICE@AliRoot::v5-03-50-AN";

  TProof::Mgr(proofCluster.Data())->SetROOTVersion("VO_ALICE@ROOT::v5-34-01-1"); 
  // Connecting to proof cluster

  TProof::Open(proofCluster.Data()); // all available workers
  // TProof::Open(proofCluster.Data(), "workers=1x"); // one job per worker
  // TProof::Open(proofCluster.Data(), "masteronly"); // master node only
  if (!gProof) {
    Error("Pi0FlowCAF.C", "Connection to proof cluster failed!");
    return;
  }

  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice"));
  // list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice:PWGGAPHOSTasks"));

  gProof->EnablePackage(alirootVersion.Data(),list);

  cout << "Pi0Flow: processing dataset " << dataset << endl;
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");
  gProof->Load("AliCaloPhoton.cxx+g");
  gProof->Load("AliPHOSEPFlattener.cxx+g");
  gProof->Load("AliAnalysisTaskPi0Flow.cxx++g");
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  
  // AOD input handler
  AliAODInputHandler* inputHandler = new AliAODInputHandler();
  mgr->SetInputEventHandler(inputHandler);
  
  // Debug level
  mgr->SetDebugLevel(1);

  // Add pi0 flow task
  AliAnalysisTaskPi0Flow *task1 = new AliAnalysisTaskPi0Flow("Pi0Spectrum");
  // Reduce binning for reduece memory footprint
  const int kNEdges = 2;
  Double_t cbin[kNEdges] = {0., 10.,};
  TArrayD tbin(kNEdges, cbin);
  Int_t    nMixed[kNEdges-1] = {4};
  TArrayI tNMixed(kNEdges-1, nMixed);

  task1->SetCentralityBinning(tbin, tNMixed);
  task1->SelectCollisionCandidates(AliVEvent::kCentral);
  
  mgr->AddTask(task1);

  mgr->ConnectInput(task1, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PHOSPi0Flow",
							    TList::Class(), 
							    AliAnalysisManager::kOutputContainer, 
							    "Pi0Flow.root");
  mgr->ConnectOutput(task1, 1, coutput1);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("proof", dataset, nEvents, nEventsSkip);
  }
  
}
