void runAnaEPFlatenningProof(Int_t mode = 0, const char *folder = "/alice/data",
			     const char *dataset = "LHC10h_000138396_hlt_clustering",
			     Int_t workers=28,
			     Bool_t usePS = kFALSE,
			     const char *minBias = "CPBI",
			     Int_t firstFile = 0, Int_t lastFile = -1)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libOADB");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include ");

  if (mode==0 || mode==2) {
    // Connect to Proof
    gEnv->SetValue("XSec.GSI.DelegProxy","2");
    Char_t *alienuser = gSystem->Getenv("alien_API_USER");
    cout<<"==> Your AliEn username is: "<<alienuser<<endl;

    //    TProof::Mgr("alice-caf.cern.ch")->SetROOTVersion("VO_ALICE@ROOT::v5-33-02a");

    TProof *p = TProof::Open(alienuser!=0 ? Form("%s@alice-caf.cern.ch",
						 alienuser) : "alice-caf.cern.ch",
			     workers>0 ? Form("workers=%d",workers) : "");

    gProof->EnablePackage("VO_ALICE@AliRoot::v5-03-23-AN");

    if (1) {
      gProof->Exec("TGrid::Connect(\"alien://\")",kTRUE);
      gSystem->Setenv("OADB_PATH","alien:///alice/cern.ch/user/c/cheshkov/OADB");
      gProof->Exec("gSystem->Setenv(\"OADB_PATH\",\"alien:///alice/cern.ch/user/c/cheshkov/OADB\")",kTRUE);
    }
  }

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("EPFlatenningAnalysis");

  Bool_t esdData = kTRUE;
  TString tempStr1 = folder;
  if (tempStr1.Contains("AOD")) esdData = kFALSE;
  TString tempStr2 = dataset;
  if (tempStr2.Contains("AOD")) esdData = kFALSE;
  if (esdData) {
    AliESDInputHandler* esdH = new AliESDInputHandler();
    esdH->SetInactiveBranches("FMD AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend.fTracks Kinks MuonTracks TrdTracks");
    esdH->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);
  }
  else {
    AliAODInputHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  }

  // physics and centrality selection
  if (esdData && usePS) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(kFALSE);
  }

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

  // Create task
  if (mode==0 || mode==2) {
    if (0) {
      gProof->Load(Form("%s/STEER/STEERBase/AliEventplane.cxx++g",
			gSystem->Getenv("ALICE_PHYSICS")));
      gProof->Load(Form("%s/ANALYSIS/AliVZEROEPSelectionTask.cxx++g",
			gSystem->Getenv("ALICE_PHYSICS")));
    }
    gProof->Load(Form("%s/AliAnaVZEROEPFlatenning.cxx++g",
		      gSystem->pwd()));
  }
  else {
    gROOT->LoadMacro(Form("%s/AliAnaVZEROEPFlatenning.cxx++g",
			  gSystem->pwd()));
  }

  gROOT->LoadMacro("$ALICE_PHYSICS/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
  AliVZEROEPSelectionTask *selTask = AddTaskVZEROEPSelection();

  AliAnaVZEROEPFlatenning *task = new AliAnaVZEROEPFlatenning("AliAnaVZEROEPFlatenning");
  task->SetMBTrigName(minBias);
  if (usePS) task->SetUsePhysSel(kTRUE);
  task->SetInput("/home/cheshkov/alice/AliRoot/PWGPP/VZERO/save/VZERO.EPFlatenning.PS.LHC11h_000170162_p1_muon_.root");
  
  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("coutput", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, (!usePS) ? Form("VZERO.EPFlatenning.%s.root",dataset) : Form("VZERO.EPFlatenning.PS.%s.root",dataset));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(3);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if (mode==0)
    mgr->StartAnalysis("proof", Form("%s/%s",folder,dataset));
  else {
    TGrid::Connect("alien://");
    TChain *chain = new TChain("esdTree");
    TGridResult *res = gGrid->Query(folder,"AliESDs.root");
    Int_t nFiles = res->GetEntries();
    if (lastFile < 0) lastFile = nFiles - 1;
    for(Int_t iFile = firstFile; iFile <= lastFile; ++iFile) {
      TString filename = res->GetKey(iFile, "turl");
      if(filename == "") continue;
      chain->AddFile(filename.Data());
    }
    if (mode==2) {
      gProof->Exec("TGrid::Connect(\"alien://\")",kTRUE);
      mgr->StartAnalysis("proof", chain);
    }
    else
      mgr->StartAnalysis("local", chain);
  }
}

