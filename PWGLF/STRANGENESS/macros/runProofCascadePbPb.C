void runProofCascadePbPb(
                     TString  proofCluster      = "mnicassi@alice-caf.cern.ch",//kiaf.sdfarm.kr", //skaf.saske.sk"
                     TString  alirootVer        = "VO_ALICE@AliRoot::v5-30-01-AN",
                     TString  rootVer           = "VO_ALICE@ROOT::v5-30-06-1", 
                     TString  dataset           = "/alice/sim/LHC11f5_000139514", 
                     TString  outFileMC         = "CascadePerformance.root",
                     TString  outFileData       = "CascadeAna.root",
                     Bool_t   runperformancetask= kTRUE, 
                     Bool_t   useMC             = kTRUE, 
                     Bool_t   dataonalien       = kFALSE,
                     Float_t  centrlowlim       = 0., 
                     Float_t  centruplim        = 90., 
                     TString  centrest          = "V0M",
                     Float_t  vtxlim            = 10.,  
                     Bool_t   kextrasel         = kFALSE,
                     Bool_t   acccut            = kFALSE,
                     Bool_t   krelaunchvertexers= kFALSE,
                     Int_t    nEvents           = 1.0*1e7, 
                     Int_t    nEventsSkip       = 0) { 

  gEnv->SetValue("XSec.GSI.DelegProxy","2");

  TString alirootMode = "";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
  TString extraLibs;
  TList *list = new TList();
  alirootMode="ALIROOT";
  extraLibs+= "ANALYSIS:OADB:ANALYSISalice:CORRFW";  
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
           vtxlim, kextrasel,
           runperformancetask, acccut, krelaunchvertexers);

}

//________________________________________________________________________
void Analysis(TString dataset, TString outFileMC, TString outFileData, 
              Bool_t useMC, Int_t nEvents, Int_t nEventsSkip, 
              Float_t centrlowlim, Float_t centruplim, TString centrest,
              Float_t vtxlim, 
              Bool_t kextrasel, Bool_t runperformancetask, Bool_t acccut, Bool_t krelaunchvertexers) {


  TString format = GetFormatFromDataSet(dataset);

  // ALICE stuff
  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test train");

  InputHandlerSetup(format,runperformancetask);

  // compile analysis task
  if (runperformancetask) gProof->Load("AliAnalysisTaskCheckPerformanceCascadePbPb.cxx++");
  else gProof->Load("AliAnalysisTaskCheckCascadePbPb.cxx++");

  // physics selection
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(useMC);

  // centrality selection
  cout<<"Format"<<format.Data()<<endl;
  if (!format.CompareTo("ESD")) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentr = AddTaskCentrality();
    if (useMC) {
      taskCentr->SetMCInput();
      taskCentr->DontUseCleaning(); // for injected MC
    }
  }

  
  // add PID response task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(useMC);

  // create task
  if (runperformancetask) {
    AliAnalysisTaskCheckPerformanceCascadePbPb *task = new AliAnalysisTaskCheckPerformanceCascadePbPb("AliAnalysisTaskCheckPerformanceCascadePbPb");
    task->SetApplyAccCut                 (acccut);
    task->SetRejectEventPileUp(kFALSE);

  } else {
    AliAnalysisTaskCheckCascadePbPb *task = new AliAnalysisTaskCheckCascadePbPb("AliAnalysisTaskCheckCascadePbPb");
  }

  task->SetAnalysisType               (format);
  task->SetRelaunchV0CascVertexers    (krelaunchvertexers);                 
  task->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
  task->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex 
  task->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
  task->SetQualityCut80TPCcls         (kTRUE);             // rejects tracks that have less than 80 clusters in the TPC
  task->SetExtraSelections            (kextrasel);         // used to add other selection cuts 
  task->SetCentralityLowLim           (centrlowlim);       // setting centrality selection variables
  task->SetCentralityUpLim            (centruplim);
  task->SetCentralityEst              (centrest);
  task->SetVertexRange                (vtxlim);

  // create output container
  if (runperformancetask) AliAnalysisDataContainer *output = mgr->CreateContainer("clist", TList::Class(), AliAnalysisManager::kOutputContainer, outFileMC);
  else       AliAnalysisDataContainer *output = mgr->CreateContainer("clist", TList::Class(), AliAnalysisManager::kOutputContainer, outFileData);

  // add task to the manager
  mgr->AddTask(task);

  task->SelectCollisionCandidates();

  // connect input and output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  // run analysis
  mgr->InitAnalysis();
  // process dataset  
  mgr->StartAnalysis("proof", dataset.Data(), nEvents, nEventsSkip);  // single dataset
  //mgr->StartAnalysis("proof","/alice/sim/LHC11f5_000139514|/alice/sim/LHC11f5_000139517",nEvents, nEventsSkip);  // multiple dataset

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
