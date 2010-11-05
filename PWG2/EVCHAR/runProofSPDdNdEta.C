void runProofSPDdNdEta (TString proofCluster="mnicassi@alice-caf.cern.ch",//skaf.saske.sk", 
          TString alirootVer="VO_ALICE@AliRoot::v4-21-02-AN", TString rootVer="VO_ALICE@ROOT::v5-27-06a-1", 
          TString dataset="/PWG4/morsch/MC_LHC10f10a",///alice/sim/LHC10g2d_130844",///PWG4/morsch/MC_LHC10f10a",///PWG2/mnicassi/LHC10g1a_130844",//"/alice/sim/LHC10f8a_130844",
          Bool_t useMC=kTRUE, Bool_t kpbpb=kTRUE, Bool_t readtr=kFALSE, Bool_t recotracklets = kFALSE, Bool_t dataonalien = kFALSE,
          Float_t centrlowlim=0., Float_t centruplim=5., TString centrest="",
          Int_t nEvents=1.0*1e7, Int_t nEventsSkip=0) { //1.0*1e7


      gEnv->SetValue("XSec.GSI.DelegProxy","2");

      TString alirootMode="";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
      TString extraLibs;
      TList *list = new TList();
      if (recotracklets) {
        alirootMode="REC";     // $ALICE_ROOT/macros/loadlibsrec.C
        extraLibs= "ITSrec:CDB:Geom:"; // not needed in default aliroot mode
      } else {
        alirootMode="ALIROOT";
      }
      extraLibs+= "ANALYSIS:ANALYSISalice";
      // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
      list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
      list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
      if (recotracklets||dataonalien) list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));

      //REM: same version of AliRoot on client!
      TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data()); //If not using the default version
//      TProof::Open(proofCluster.Data());
      // enable n workers per machine
//      TProof::Open(proofCluster.Data(),"workers=nx")       
      // enable less workers
      TProof::Open(proofCluster.Data(),"workers=20"); //For performance reasons, try to avoid it.
      if (!gProof) {
        Error("runSPDdNdEtaAna.C","Connection to AF failed.");
        return;
      }
      gProof->EnablePackage(alirootVer.Data(), list);

      Analysis(dataset.Data(), useMC, kpbpb, readtr, recotracklets, nEvents, nEventsSkip, centrlowlim, centruplim, centrest);

}

void Analysis(TString dataset, Bool_t useMC, Bool_t kpbpb, Bool_t readtr, Bool_t recotracklets, Int_t nEvents, Int_t nEventsSkip, Float_t centrlowlim, Float_t centruplim, TString centrest) {

  TString format = GetFormatFromDataSet(dataset);

  // ALICE stuff
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Test train");

  InputHandlerSetup(format,useMC,recotracklets);

  // compile the tracklet reconstruction class
  gProof->Load("AliTrackletAlg.cxx++");
  // compile analysis task
  gProof->Load("AliAnalysisTaskSPDdNdEta.cxx++");

  // create manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("SPD analysis");

  // create task
  AliAnalysisTaskSPDdNdEta *task = new AliAnalysisTaskSPDdNdEta("AliAnalysisTaskSPDdNdEta");

  AliTriggerAnalysis::Trigger trigg = AliTriggerAnalysis::kAcceptAll; // to be changed every time
  AliAnalysisTaskSPDdNdEta::MCCentralityBin kmccentrbin = AliAnalysisTaskSPDdNdEta::kall; // idem

  task->SetReadMC(useMC);
  task->SetTrigger(trigg);
  task->SetReadPbPb(kpbpb);
  task->SetReadTrackRefs(readtr);
  task->SetRecoTracklets(recotracklets);
  task->SetMCCentralityBin(kmccentrbin);
  task->SetCentralityLowLim(centrlowlim);
  task->SetCentralityUpLim(centruplim);
  task->SetCentralityEst(centrest);

  // physics selection
//  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
//  physicsSelectionTask = AddTaskPhysicsSelection(useMC);

  // centrality selection
  gProof->Load("$ALICE_ROOT/ANALYSIS/AliCentralitySelectionTask.cxx++g");
  AliCentralitySelectionTask *taskCentr = new AliCentralitySelectionTask("CentralitySelection");
  taskCentr->SetPercentileFile("$ALICE_ROOT/ANALYSIS/macros/test_AliCentralityBy1D.root");
  taskCentr->SetPercentileFile2("./test_AliCentralityByFunction.root");
  mgr->AddTask(taskCentr);
  mgr->ConnectInput (taskCentr,0, mgr->GetCommonInputContainer());


  // create output container
  if (useMC) AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, "SPDdNdEtaCorr.root");
  else AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, "SPDdNdEtaData.root");

  // add task to the manager
  mgr->AddTask(task);

  // connect input and output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output);

  // run analysis
  mgr->InitAnalysis();
  // process dataset  
  mgr->StartAnalysis("proof", dataset.Data(), nEvents, nEventsSkip);

}


TString GetFormatFromDataSet(TString dataset) {

//   Info("runAAF.C","Detecting format from dataset (may take while, depends on network connection)...");
  TString dsTreeName;
  if (dataset.Contains("#")) {
    Info("runSPDdNdEta.C",Form("Detecting format from dataset name '%s' ...",dataset.Data()));
    dsTreeName=dataset(dataset.Last('#'),dataset.Length());
  } else {
    Info("runSPDdNdEta.C",Form("Detecting format from dataset '%s' (may take while, depends on network connection) ...",dataset.Data()));
    TFileCollection *ds = gProof->GetDataSet(dataset.Data());
    if (!ds) {
      Error(Form("Dataset %s doesn't exist on proof cluster!!!!",dataset.Data()));
      return "";
    }
    dsTreeName = ds->GetDefaultTreeName();
  }

  if (dsTreeName.Contains("esdTree")) {
    Info("runSPDdNdEta.C","ESD input format detected ...");
    return "ESD";
  } else if (dsTreeName.Contains("aodTree"))  {
    Info("runSPDdNdEta.C","AOD input format detected ...");
    return "AOD";
  } else {
    Error("runSPDdNdEta.C",Form("Tree %s is not supported !!!",dsTreeName.Data()));
    Error("runSPDdNdEta.C",Form("Maybe set your DS to %s#esdTree or %s#aodTree",dataset.Data(),dataset.Data()));
  }

  return "";
}

Bool_t InputHandlerSetup(TString format = "esd", Bool_t useKine = kTRUE, Bool_t recotracklets)
{
  format.ToLower();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisDataContainer *cin = mgr->GetCommonInputContainer();

  if (cin) return;

  if (!format.CompareTo("esd"))
  {
    AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdInputHandler)
    {
      Info("CustomAnalysisTaskInputSetup", "Creating esdInputHandler ...");
      if(recotracklets) esdInputHandler = new AliESDInputHandlerRP();
      else esdInputHandler = new AliESDInputHandler();
      mgr->SetInputEventHandler(esdInputHandler);
    }

    if (useKine)
    {
      AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

      if (!mcInputHandler)
      {
        Info("CustomAnalysisTaskInputSetup", "Creating mcInputHandler ...");
        AliMCEventHandler* mcInputHandler = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcInputHandler);
      }
    }

  }
  else if (!format.CompareTo("aod"))
  {
    AliAODInputHandler *aodInputHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!aodInputHandler)
    {
      Info("CustomAnalysisTaskInputSetup", "Creating aodInputHandler ...");
      aodInputHandler = new AliAODInputHandler();
      mgr->SetInputEventHandler(aodInputHandler);
    }
  }
  else
  {
    AliWarning("Wrong input format!!! Only ESD and AOD are supported. Skipping Task ...");
    return kFALSE;
  }

  return kTRUE;
}
