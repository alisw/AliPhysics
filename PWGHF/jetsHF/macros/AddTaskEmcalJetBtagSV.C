AliAnalysisTaskEmcalJetBtagSV2 *AddTaskEmcalJetBtagSV(const char *trkcontname   = "tracks",
                                                     const char *jetcontname   = "Jets",
                                                     const char *mctrkcontname = "mcparticles",
                                                     const char *mcjetcontname = "MCJets",
                                                     Double_t jetRadius = 0.4,
                                                     const char *type  = "TPC",
                                                     TString fileout   = "standard",
                                                     Bool_t corrMode   =  kFALSE,
                                                     Bool_t doBkgRej   =  kTRUE,
                                                     Bool_t doQAvtx    =  kFALSE,
                                                     Bool_t doFillV0   =  kFALSE,
                                                     Bool_t checkXsec  =  kFALSE,
                                                     Bool_t useWeight  =  kFALSE,
                                                     const char *patt  =  "",
                                                     Int_t debug       = 1,
                                                     Double_t tagRadius = 0.4,
                                                     TString cutflname = "",
                                                     Float_t minC = 0., Float_t maxC = 100.,
                                                     const char *taskname = "HFjetsContainer")
{
  // Mailto: ycorrale@cern.ch

  // Get the AnalysisManager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {

    ::Error("AddTaskEmcalJetBtagSV", "No analysis manager to connect to.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetBtagSV", "This task requires an input event handler");
    return NULL;
  }
  TString name(taskname);
  // Configure analysis task
  AliAnalysisTaskEmcalJetBtagSV2 *hfTask = new AliAnalysisTaskEmcalJetBtagSV2(name);
 
  AliParticleContainer *trkCont  = hfTask->AddTrackContainer(trkcontname);//for data, and reco MC

  TString strType(type);
  AliJetContainer *jetCont = hfTask->AddJetContainer(jetcontname, strType, jetRadius);
  jetCont->ConnectParticleContainer(trkCont);
  if (corrMode) {

    AliParticleContainer *mctrackCont  = hfTask->AddMCParticleContainer(mctrkcontname);//for MC particle level
    
    TString strType(type);
    AliJetContainer *mcjetCont = hfTask->AddJetContainer(mcjetcontname, strType, jetRadius);
  }

  // Set data or corrections mode
  hfTask->SetCorrectionMode(corrMode); // kFALSE for real data
  hfTask->SetDoBkgRejection(doBkgRej);
  hfTask->SetDoFillSecVtxQA(doQAvtx);
  hfTask->SetDoFillV0Trks(doFillV0);
  hfTask->SetCheckMCCrossSection(checkXsec);
  if (useWeight) hfTask->SetUseWeightOn();
  hfTask->SetGenNamePattern(patt);
  //  hfTask->SetDebugLevel(10);
  hfTask->SetDebugLevel(debug);

  hfTask->SetJetContName(jetcontname);
  hfTask->SetTrkContName(trkcontname);

  // choose the method to select on the flavor of the jet and, if needed, select the pt-hard bin via the minimum pt-hard
  // also set mc jet and particle container
  if (corrMode) {

    hfTask->SetMCJetContName(mcjetcontname);
    hfTask->SetMCTrkContName(mctrkcontname);
    hfTask->SetJetTaggingRadius(tagRadius);
  }

  // Define the tagger
  AliHFJetsTaggingVertex *tagger = new AliHFJetsTaggingVertex();

  // Set analysis cuts
  TString strCutFlName(cutflname);
  if (!strCutFlName.IsNull() && !gSystem->AccessPathName(strCutFlName.Data(), kFileExists)) {

    // read cuts from file
    ::Info(Form("Reading cuts from file: %s", strCutFlName.Data()));
    
    TFile *f = TFile::Open(strCutFlName.Data()); 
    AliRDHFJetsCuts *cuts= (AliRDHFCutsD0toKpi *)f->Get("HFJetsCutsVertex");

    // Set centrality
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetTriggerMask(AliVEvent::kAny);
    cuts->SetTriggerClass(trigClass);
    cuts->PrintAll();
    hfTask->SetCuts(cuts);
  } 
  else {

    // define your cuts here
    DefineCutsTask(hfTask, minC, maxC, corrMode);
    DefineCutsTagger(tagger);
  }

  // // Add task to manager
  hfTask->SetTagger(tagger);
  mgr->AddTask(hfTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cInput =   mgr->GetCommonInputContainer();
  mgr->ConnectInput(hfTask, 0, cInput);

  //All containers in a list
  AliAnalysisDataContainer *cOutput = mgr->CreateContainer(name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(hfTask, 1, cOutput);

  delete tagger;
  return hfTask;
}

//------------------------------------------------------
Bool_t DefineCutsTask(AliAnalysisTaskEmcalJetBtagSV2 *task, 
                      Float_t minC, Float_t maxC, Bool_t corrMode) 
{

    // define cuts for task
    AliRDHFJetsCuts *cuts=new AliRDHFJetsCuts();
    // jets
    cuts->SetJetRadius(0.4); // this cut does nothing
    cuts->SetMaxEtaJet(0.5);//0.9-R
    cuts->SetMinPtJet(5.);//
    cuts->SetMaxPtJet(100.);
    // Set centrality
    cuts->SetMinCentrality(minC);
    cuts->SetMaxCentrality(maxC);
    cuts->SetUsePhysicsSelection(kFALSE);
    if (corrMode) cuts->SetTriggerMask(AliVEvent::kAny);
    else {
      cuts->SetOptPileup(1);
      cuts->ConfigurePileupCuts(5, 0.8);
      // cuts->SetUsePhysicsSelection(kFALSE);
      cuts->SetTriggerClass("CINT7");
      cuts->SetTriggerMask(AliVEvent::kINT7);
    } // pPb minbias only
  task->SetCuts(cuts);
  delete cuts;
  return kTRUE;
}

//------------------------------------------------------
Bool_t DefineCutsTagger(AliHFJetsTaggingVertex *tg)
{
  AliRDHFJetsCutsVertex *cuts2 = new AliRDHFJetsCutsVertex("jetCuts");

  AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(1.0,1.e10);

  cuts2->AddTrackCuts(esdTrackCuts);

  // vertexing
  cuts2->SetNprongs(3);
  cuts2->SetIsElec(kFALSE); // kTRUE to select e in jet vertex
  cuts2->SetMinPtHardestTrack(1.0);//default 0.3
  cuts2->SetSecVtxWithKF(kFALSE);//default with StrLinMinDist
  cuts2->SetImpParCut(0.);//default 0
  cuts2->SetDistPrimSec(0.);//default 0
  cuts2->SetCospCut(-1);//default -1

  tg->SetCuts(cuts2);

  delete esdTrackCuts;

  delete cuts2;

  return kTRUE;
}
