void runProofSPDdNdEta_CorrRef4(
                       TString  proofCluster  = "mnicassi@alice-caf.cern.ch", //skaf.saske.sk",
                       TString  alirootVer    = "VO_ALICE@AliRoot::v4-21-04-AN", 
                       TString  rootVer       = "VO_ALICE@ROOT::v5-27-06a-1",
//                     TString  dataset       = "/alice/data/LHC10h_000137161_p1_plusplusplus", 
                       TString  dataset       = "/alice/sim/LHC10h8_000137161",
                       TString  outFileCorr   = "SPDdNdEtaCorr.root",
//                       TString  outFileCorr   = "SPDdNdEtaCorrRot.root",
                       TString  outFileData   = "SPDdNdEtaData.root",
                       TString  percentFile   = "./AliCentralityBy1D_LHC10g2a_100.root",
                       TString  percentFile2  = "./AliCentralityByFunction_LHC10g2a_100.root",
                       Bool_t   useMC         = kTRUE, 
                       Bool_t   readtr        = kFALSE, 
                       Bool_t   recotracklets = kTRUE, 
                       Bool_t   dataonalien   = kFALSE,
                       Float_t  centrlowlim   = 0., 
                       Float_t  centruplim    = 5., 
                       Bool_t   centrest      = kFALSE,
                       Int_t    minClMultLay2 = 4300, 
//                       Int_t    minClMultLay2 = -1,  
                       Int_t    maxClMultLay2 = 1.0*1e5, 
                       Int_t    minV0Mult     = 14674.5,
                       Float_t vtxlim         = 7.,
                       Bool_t partsp          = kTRUE, 
                       Float_t  phiWindow     = 0.8,
                       Float_t  thetaWindow   = 0.025,
                       Float_t  phiShift      = 0.0045,
                       Float_t  removeClFOvl  = kFALSE,
                       Float_t  phiOvlCut     = 0.005,
                       Float_t  zetaOvlCut    = 0.05,
                       Float_t  phiRotAngle   = 0.,  
//                       Float_t  phiRotAngle   = TMath::Pi(),
                       Float_t  phiWindowAna  = 0.08,
                       Int_t    nEvents       = 1.0*1e7, 
                       Int_t    nEventsSkip   = 0) { //1.0*1e7

  gEnv->SetValue("XSec.GSI.DelegProxy","2");

  TString alirootMode = "";    // STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice (default aliroot mode)
  TString extraLibs;
  TList *list = new TList();
  if (recotracklets) {
    alirootMode="REC";     // $ALICE_PHYSICS/macros/loadlibsrec.C
    extraLibs= "ITSrec:CDB:Geom:"; // not needed in default aliroot mode
  } else alirootMode="ALIROOT";
  extraLibs+= "ANALYSIS:ANALYSISalice";

  // sets $ALIROOT_MODE on each worker to let proof to know to run in special mode
  list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  if (recotracklets||dataonalien) list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));

  // REM: same version of AliRoot on client!
//  TProof::Mgr(proofCluster.Data())->SetROOTVersion(rootVer.Data()); //If not using the default version
  TProof::Open(proofCluster.Data());
  // enable n workers per machine
//  TProof::Open(proofCluster.Data(),"workers=nx")       
  // enable less workers
//  TProof::Open(proofCluster.Data(),"workers=20"); //For performance reasons, try to avoid it.
  if (!gProof) {
    Error("runSPDdNdEtaAna.C","Connection to AF failed.");
    return;
  }
  gProof->EnablePackage(alirootVer.Data(), list);

//  gROOT->LoadMacro("AnalysisMacroa.C");
  Analysis(dataset.Data(), outFileCorr, outFileData, percentFile, percentFile2, 
           useMC, readtr, recotracklets, 
           nEvents, nEventsSkip, centrlowlim, centruplim, centrest, minClMultLay2, maxClMultLay2, minV0Mult, vtxlim, partsp,
           phiWindow, thetaWindow, phiShift, removeClFOvl, phiOvlCut, zetaOvlCut, phiRotAngle,phiWindowAna);

}

//________________________________________________________________________
void Analysis(TString dataset, TString outFileCorr, TString outFileData, TString percentFile, TString percentFile2, 
              Bool_t useMC, Bool_t readtr, Bool_t recotracklets, 
              Int_t nEvents, Int_t nEventsSkip, 
              Float_t centrlowlim, Float_t centruplim, Bool_t centrest, Int_t minClMultLay2, Int_t maxClMultLay2, Int_t minV0Mult, Float_t vtxlim, Bool_t partsp,
              Float_t phiWindow, Float_t thetaWindow, Float_t phiShift, Bool_t removeClFOvl, 
              Float_t phiOvlCut, Float_t zetaOvlCut, Float_t phiRotAngle, Float_t phiWindowAna) {

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
  task->SetReadTrackRefs(readtr);
  task->SetRecoTracklets(recotracklets);
  task->SetMCCentralityBin(kmccentrbin);
  task->SetCentralityLowLim(centrlowlim);
  task->SetCentralityUpLim(centruplim);
  task->SetCentralityEst(centrest);
  task->SetMinClusterMultLay2(minClMultLay2);
  task->SetMaxClusterMultLay2(maxClMultLay2);
  task->SetMinV0Mult(minV0Mult); 
  task->SetVertexRange(vtxlim);
  task->SetPartSpecies(partsp);

  task->SetPhiWindow(phiWindow);
  task->SetThetaWindow(thetaWindow); 
  task->SetPhiShift(phiShift);
  task->SetRemoveClustersFromOverlaps(removeClFOvl); 
  task->SetPhiOverlapCut(phiOvlCut);
  task->SetZetaOverlapCut(zetaOvlCut);
  task->SetPhiRotationAngle(phiRotAngle);
  task->SetPhiWindowAna(phiWindowAna);

  // physics selection
  gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask=AddTaskPhysicsSelection(useMC,kFALSE);
  if (!useMC) {
    AliPhysicsSelection * physSel = physSelTask->GetPhysicsSelection();
  
/*    if (dataset=="/alice/data/LHC10h_000137045_p1_plus") {
      physSel->AddCollisionTriggerClass("+CMBS1A-B-NOPF-ALL");
      physSel->AddCollisionTriggerClass("+CMBS1C-B-NOPF-ALL");
      physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");

    } else {*/
      physSel->AddCollisionTriggerClass("+CMBS2A-B-NOPF-ALL");
      physSel->AddCollisionTriggerClass("+CMBS2C-B-NOPF-ALL");
      physSel->AddCollisionTriggerClass("+CMBAC-B-NOPF-ALL");
//    }

    task->SelectCollisionCandidates(AliVEvent::kUserDefined);
  } else {

    task->SelectCollisionCandidates(AliVEvent::kMB);
  }

  // centrality selection
  gProof->Load("./AliCentralitySelectionTask.cxx++g");
  AliCentralitySelectionTask *taskCentr = new AliCentralitySelectionTask("CentralitySelection");
  taskCentr->SetPercentileFile(percentFile);
  taskCentr->SetPercentileFile2(percentFile2);

  mgr->AddTask(taskCentr);
  mgr->ConnectInput (taskCentr,0, mgr->GetCommonInputContainer());

  // create output container
  if (useMC) AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, outFileCorr);
  else       AliAnalysisDataContainer *output = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer, outFileData);

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

//________________________________________________________________________
TString GetFormatFromDataSet(TString dataset) {

//  Info("runAAF.C","Detecting format from dataset (may take while, depends on network connection)...");
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

//________________________________________________________________________
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

//________________________________________________________________________
// MEMO DATASETS ...
/*
/alice/data/LHC10h_000137045_p1           |      60 | /esdTree     | 5.144e+05|    12 GB |  100 %
/alice/data/LHC10h_000137045_p1_lowflux   |      50 | /esdTree     | 4.144e+05|    10 GB |  100 %
/alice/data/LHC10h_000137045_p1_plus      |      53 | /esdTree     | 4.444e+05|    10 GB |  100 %
/alice/data/LHC10h_000137124_p1           |      89 | /esdTree     | 7.152e+04|     3 GB |  100 %
/alice/data/LHC10h_000137125_p1           |      41 | /esdTree     | 1.537e+04|     1 GB |  100 %
/alice/data/LHC10h_000137132_p1           |      87 | /esdTree     | 8.881e+04|     8 GB |  100 %
/alice/data/LHC10h_000137133_p1           |     306 | /esdTree     | 3.521e+05|    29 GB |  100 %
/alice/data/LHC10h_000137161              |    3348 | /esdTree     | 6.279e+05|   407 GB |   97 %

/alice/sim/LHC10h1_000137045              |    1198 | /esdTree     | 1.198e+04|   138 GB |  100 %
/alice/sim/LHC10h1a_000137045             |    1022 | /esdTree     | 1.006e+04|   124 GB |   98 %

*/
//________________________________________________________________________