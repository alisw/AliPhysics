void runAnaTrigProof(const char *dataset = "LHC10h_000138396_hlt_clustering", Int_t workers=28, Long64_t nentries=10000000, Long64_t firstentry=0)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include ");

  // Connect to Proof
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  Char_t *alienuser = gSystem->Getenv("alien_API_USER");
  cout<<"==> Your AliEn username is: "<<alienuser<<endl;

  TProof *p = TProof::Open(alienuser!=0 ? Form("%s@alice-caf.cern.ch",
					       alienuser) : "alice-caf.cern.ch",
			   workers>0 ? Form("workers=%d",workers) : "");

  //  gProof->GetManager()->SetROOTVersion("VO_ALICE@ROOT::v5-28-00f");
  gProof->EnablePackage("VO_ALICE@AliRoot::v5-02-07-AN");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AliAnaFwdDet");

  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetInactiveBranches("FMD AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend.fTracks Kinks MuonTracks TrdTracks");
  //  esdH->SetReadFriends(kTRUE);
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  // physics and centrality selection
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(kFALSE);
  AliOADBTriggerAnalysis * oadbTrigAnalysis = new
    AliOADBTriggerAnalysis("CustomTA");
  oadbTrigAnalysis->SetZDCCorrParameters(-0.823, -0.653, 4*0.58, 4*0.5);
  physicsSelectionTask->GetPhysicsSelection()->SetCustomOADBObjects(0,0,oadbTrigAnalysis);


  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

  // Create task
  gProof->Load(Form("%s/AliAnaVZEROTrigger.cxx++g",
		    gSystem->pwd()));
  AliAnaVZEROTrigger *task = new AliAnaVZEROTrigger("AliAnaVZEROTrigger");
  task->SetMBTrigName("CMBACS2");
  task->Setup("trigger.txt");
  
  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("coutput", TList::Class(), 
    AliAnalysisManager::kOutputContainer, Form("VZERO.Trigger.%s.root",dataset));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);


  // Enable debug printouts
  //  mgr->SetDebugLevel(3);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("proof", Form("/alice/data/%s",dataset), nentries, firstentry);
}

