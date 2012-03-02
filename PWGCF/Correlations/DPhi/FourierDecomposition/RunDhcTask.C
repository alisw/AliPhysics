// Macro to run AliDhcTask - A. Adare, Mar 2012
// Choose dataType = "esd" or "aod".
// Choose runMode  = "local", "proof", "grid", or "terminate".
// Choose proofCluster = "skaf.saske.sk" or "alice-caf.cern.ch".
// When running on proof: use root, not aliroot.

Int_t verbosity      = 0;
Long64_t nEvents     = 123456789;
TString runMode      = "local"; // "proof";
TString dataType     = "aod";
TString localDir     = "/Users/adare/esd/alice/data/2010/LHC10h/000139107/ESDs/pass2/";
TString esdDir       = localDir + TString("10000139107001.120/");
TString aodDir       = localDir + TString("AOD049/0008/");
TString proofDataset = "/alice/data/LHC10h_000137848_p2_AOD049";
TString proofCluster = "aadare@skaf.saske.sk"; // or "aadare@alice-caf.cern.ch";
TString alirootVer   = "VO_ALICE@AliRoot::v4-21-21-AN"; // for PROOF
TString libsBase     = "Tree:Geom:VMC:STEERBase:ESD:AOD";
TString libsExtra    = "ANALYSIS:ANALYSISalice:PWGCFCorrelationsDPhi";
TList* libList = new TList();
TChain* chain = 0;

void RunDhcTask()
{
  // Creates the large file memstat_<pid>.root.
  // Call TMemStat::Show() on the promp to analyze it.
  if (0)
    TMemStat mm("gnubuiltin");

  if (runMode=="local")
    LocalSetup();
  else if (runMode=="proof") {
    ProofSetup();
  }

  AddDhcTask();
  
  /// ---  All set up now...run Analysis!  ---
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr->InitAnalysis()) 
    return;
  if (verbosity > 0)
    mgr->PrintStatus();
  
  if (runMode=="local") {
    mgr->StartAnalysis(runMode.Data(), chain, nEvents);
  }
  else if (runMode=="proof") {
    mgr->StartAnalysis(runMode.Data(), proofDataset.Data(), nEvents);
  }
  return;
}

void LocalSetup()
{
  // Load libraries 
  TObjArray* libs = libsExtra.Tokenize(":");
  for (Int_t i=0; i<libs->GetEntries(); i++) {
    const char* name = libs->At(i)->GetName();
    Int_t ret = gSystem->Load(Form("lib%s", name));
    if (ret > 0)
      Warning("LocalSetup()", "lib%s already loaded", name);
    if (ret < 0)
      Warning("LocalSetup()", "lib%s failed to load", name);
  }
  
  // Add headers
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  
  // Load macros/tasks
  if (dataType == "esd") {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  }
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  // gROOT->LoadMacro("AliPool.cxx+g");  
  // gROOT->LoadMacro("AliDhcTask.cxx+g");  
  
  // Setup chain
  TString chainName = dataType.Contains("aod") ? "aodTree" : "esdTree";
  chain = new TChain(chainName.Data());
  if (dataType.Contains("aod"))
    chain->Add(Form("%sAliAOD.root", aodDir.Data()));
  if (dataType.Contains("esd"))
    chain->Add(Form("%sAliESDs.root", esdDir.Data()));
  return;
}

void ProofSetup()
{
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
      
  // Get PROOF pointer
  TProof *proof = gProof;
  if (proof) {
    TProof::Reset(proofCluster.Data()); // add kTRUE as a 2nd arg for hard reset
    //    proof->Reset("");
    cerr << "Trying to reset proof" << endl;
  } else {
    proof = TProof::Open(proofCluster.Data(), "workers=60");
  }
  if (!proof) {
    cerr << "Connection to " << proofCluster.Data()
	 << " failed!" << endl;
    return;
  }

  if (0) {
    proof->ShowPackages();
    proof->ShowDataSets();
  }
  
  // Load libraries   
  int  ret = 0;
  libList->Add(new TNamed("ALIROOT_EXTRA_LIBS", libsBase.Data()));
  libList->Add(new TNamed("ALIROOT_EXTRA_LIBS", libsExtra.Data()));

  ret = proof->EnablePackage(alirootVer.Data(), libList);
  
  proof->EnablePackage(alirootVer.Data(), 0);
  if (ret) {
    Error("ProofSetup()", "Failed to load all libs.");
    return;
  }

  // Load macros/tasks
  if (dataType == "esd") {
    proof->Load("AddTaskCentrality.C");
  }
  proof->Load("AddTaskPhysicsSelection.C");
  proof->Load("AliPool.cxx+g");  
  proof->Load("AliDhcTask.cxx+g", 0);
  return;
}

void AddDhcTask()
{
  // Need the following macros loaded beforehand:
  // "$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C" (ESD only)
  // "$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C"
  // "AliDhcTask.cxx+g"

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    mgr = new AliAnalysisManager("DhcAnalysis");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  AliAODInputHandler* aodH = new AliAODInputHandler();
  
  if (dataType.Contains("aod"))
    mgr->SetInputEventHandler(aodH);
  if (dataType.Contains("esd"))
    mgr->SetInputEventHandler(esdH);
  
  if (dataType == "esd") {
    AliCentralitySelectionTask *centTask = AddTaskCentrality();
    centTask->SelectCollisionCandidates(AliVEvent::kAny);
    centTask->SetPass(2);
  }
  
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  
  AliDhcTask *dhcTask = new AliDhcTask("DhcTask");
  dhcTask->SelectCollisionCandidates(AliVEvent::kMB);
  dhcTask->SetVerbosity(verbosity);

  // Add task(s)
  mgr->AddTask(dhcTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("dhclist", TList::Class(),   
			 AliAnalysisManager::kOutputContainer, 
			 "dhctask_output.root");
  
  // Connect input/output
  mgr->ConnectInput(dhcTask, 0, cinput);
  mgr->ConnectOutput(dhcTask, 1, coutput);
  
  // Enable debug printouts
  mgr->SetDebugLevel(verbosity);

  return;
}

