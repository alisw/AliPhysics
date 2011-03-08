void runProofCascadePbPb (
                          TString  proofCluster      = "mnicassi@alice-caf.cern.ch", //skaf.saske.sk",
                          TString  alirootVer        = "VO_ALICE@AliRoot::v4-21-16-AN", 
                          TString  rootVer           = "VO_ALICE@ROOT::v5-27-06d",
                          TString  dataset           = "/alice/sim/LHC11a7_000137161",
//                          TString  dataset           = "/alice/data/LHC10h_000137161_p1_5plus",                        
                          TString  outFileMC         = "CascadePerformance.root",
                          TString  outFileData       = "CascadeAna.root",
                          Bool_t   runperformancetask= kFALSE, 
                          Bool_t   useMC             = kFALSE, 
                          Bool_t   dataonalien       = kFALSE,
                          Float_t  centrlowlim       = 0., 
                          Float_t  centruplim        = 20., 
                          TString  centrest          = "V0M",
                          Short_t  lCollidingSystems = 1,       //0 = pp, 1 = AA
                          Float_t  vtxlim            = 15.,  
                          Bool_t   usecfcontainers   = kTRUE,
                          Bool_t   kextrasel         = kTRUE,
                          Bool_t   acccut            = kFALSE,
                          Int_t    nEvents           = 5000, 
                          Int_t    nEventsSkip       = 0) { //1.0*1e7

  gEnv->SetValue("XSec.GSI.DelegProxy","2");

  TString alirootMode = "";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
  TString extraLibs;
  TList *list = new TList();
  alirootMode="ALIROOT";
  extraLibs+= "ANALYSIS:ANALYSISalice:CORRFW";  

  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  if (dataonalien) list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));

  // REM: same version of AliRoot on client!
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data()); //If not using the default version, do it the first time only
  TProof::Open(proofCluster.Data());
  // enable n workers per machine
//  TProof::Open(proofCluster.Data(),"workers=nx")       
  // enable less workers
//  TProof::Open(proofCluster.Data(),"workers=20"); //For performance reasons, try to avoid it.
  if (!gProof) {
    Error("runProof.C","Connection to AF failed.");
    return;
  }
  gProof->EnablePackage(alirootVer.Data(), list);

  Analysis(dataset.Data(), outFileMC, outFileData, 
           useMC, nEvents, nEventsSkip,
           centrlowlim, centruplim, centrest, 
           vtxlim, usecfcontainers, lCollidingSystems, kextrasel,
           runperformancetask, acccut);

}

//________________________________________________________________________
void Analysis(TString dataset, TString outFileMC, TString outFileData, 
              Bool_t useMC, Int_t nEvents, Int_t nEventsSkip, 
              Float_t centrlowlim, Float_t centruplim, TString centrest,
              Float_t vtxlim, Bool_t usecfcontainers, Short_t  lCollidingSystems,
              Bool_t kextrasel, Bool_t runperformancetask, Bool_t acccut) {

  TString format = GetFormatFromDataSet(dataset);

  // ALICE stuff
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test train");

  InputHandlerSetup(format,useMC);

  // compile analysis task
  if (runperformancetask) gProof->Load("AliAnalysisTaskCheckPerformanceCascadePbPb.cxx++");
  else gProof->Load("AliAnalysisTaskCheckCascadePbPb.cxx++");

  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Cascade analysis");

  // create task
  if (runperformancetask) {
    AliAnalysisTaskCheckPerformanceCascadePbPb *task = new AliAnalysisTaskCheckPerformanceCascadePbPb("AliAnalysisTaskCheckPerformanceCascadePbPb");
    task->SetDebugLevelCascade           (0);
    task->SetUseCFCont                   (usecfcontainers);
  } else {
    AliAnalysisTaskCheckCascadePbPb *task = new AliAnalysisTaskCheckCascadePbPb("AliAnalysisTaskCheckCascadePbPb");
    task->SetUseCFContCascadeCuts       (usecfcontainers);
  }

  task->SetCollidingSystems           (lCollidingSystems); // only for multiplicity binning 
  task->SetAnalysisType               (format);
  task->SetRelaunchV0CascVertexers    (0);                 // used but code is commented out
  task->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
  task->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex 
  task->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
  task->SetQualityCut80TPCcls         (kTRUE);             // rejects tracks that have less than 80 clusters in the TPC
  task->SetExtraSelections            (kextrasel);         // used to add other selection cuts 
  task->SetCentralityLowLim           (centrlowlim);       // setting centrality selection vriables
  task->SetCentralityUpLim            (centruplim);
  task->SetCentralityEst              (centrest);
  task->SetVertexRange                (vtxlim);

  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(useMC,kFALSE);//,kFALSE); // PbPb: use this second flag up to run 133 and for runs with the CMBS2* triggers (from 137135, including 137161) 
//  task->SelectCollisionCandidates();
  
  // for early PbPb runs collected with CMBS1 (up to 137133) manual settings
  if (!useMC) {
    AliPhysicsSelection * physSel = physSelTask->GetPhysicsSelection();
  
    if (dataset=="/alice/data/LHC10h_000137045_p1_plus") {
      physSel->AddCollisionTriggerClass("+CMBS1A-B-NOPF-ALL");   // on the web page
      physSel->AddCollisionTriggerClass("+CMBS1C-B-NOPF-ALL");
      physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
 
      // This are needed only to fill the statistics tables
      physSel->AddBGTriggerClass("+CMBAC-C-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1C-C-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1A-C-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBAC-A-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1C-A-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1A-A-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBAC-E-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1C-E-NOPF-ALL");
      physSel->AddBGTriggerClass("+CMBS1A-E-NOPF-ALL");
    } else {
    // we used these manual settings for run137161 (first and second mult paper)
    physSel->AddCollisionTriggerClass("+CMBS2A-B-NOPF-ALL");
    physSel->AddCollisionTriggerClass("+CMBS2C-B-NOPF-ALL");
    physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
    }

    task->SelectCollisionCandidates(AliVEvent::kUserDefined);
  } else {

    task->SelectCollisionCandidates(AliVEvent::kMB); 
  }

  // centrality selection 
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
  if (useMC) taskCentr->SetMCInput();

  // create output container
  if (runperformancetask) AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, outFileMC);
  else       AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, outFileData);

  // add task to the manager
  mgr->AddTask(task);

  // connect input and output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  // run analysis
  mgr->InitAnalysis();
  // process dataset  
  mgr->StartAnalysis("proof", dataset.Data(), nEvents, nEventsSkip);  // single dataset
//  mgr->StartAnalysis("proof","/alice/sim/LHC11a7_000137161|/alice/sim/LHC11a7_000137162|/alice/sim/LHC11a7_000137430|/alice/sim/LHC11a7_000137431|/alice/sim/LHC11a7_000137432",nEvents, nEventsSkip);  // multiple dataset

}

//________________________________________________________________________
TString GetFormatFromDataSet(TString dataset) {

//  Info("runProof.C","Detecting format from dataset (may take while, depends on network connection)...");
  TString dsTreeName;
  if (dataset.Contains("#")) {
    Info("runProof.C",Form("Detecting format from dataset name '%s' ...",dataset.Data()));
    dsTreeName=dataset(dataset.Last('#'),dataset.Length());
  } else {
    Info("runProof.C",Form("Detecting format from dataset '%s' (may take while, depends on network connection) ...",dataset.Data()));
    TFileCollection *ds = gProof->GetDataSet(dataset.Data());
    if (!ds) {
      Error(Form("Dataset %s doesn't exist on proof cluster!!!!",dataset.Data()));
      return "";
    }
    dsTreeName = ds->GetDefaultTreeName();
  }

  if (dsTreeName.Contains("esdTree")) {
    Info("runProof.C","ESD input format detected ...");
    return "ESD";
  } else if (dsTreeName.Contains("aodTree"))  {
    Info("runProof.C","AOD input format detected ...");
    return "AOD";
  } else {
    Error("runProof.C",Form("Tree %s is not supported !!!",dsTreeName.Data()));
    Error("runProof.C",Form("Maybe set your DS to %s#esdTree or %s#aodTree",dataset.Data(),dataset.Data()));
  }

  return "";
}

//________________________________________________________________________
Bool_t InputHandlerSetup(TString format = "esd", Bool_t useKine = kTRUE) {
  format.ToLower();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisDataContainer *cin = mgr->GetCommonInputContainer();

  if (cin) return;

  if (!format.CompareTo("esd")) {
    AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdInputHandler) {
      Info("CustomAnalysisTaskInputSetup", "Creating esdInputHandler ...");
      esdInputHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdInputHandler);
    }

    if (useKine) {
      AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

      if (!mcInputHandler) {
        Info("CustomAnalysisTaskInputSetup", "Creating mcInputHandler ...");
        AliMCEventHandler* mcInputHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcInputHandler);
      }
    }

  } else if (!format.CompareTo("aod")) {
    AliAODInputHandler *aodInputHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodInputHandler) {
      Info("CustomAnalysisTaskInputSetup", "Creating aodInputHandler ...");
      aodInputHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodInputHandler);
    }
  } else {
    AliWarning("Wrong input format!!! Only ESD and AOD are supported. Skipping Task ...");
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
