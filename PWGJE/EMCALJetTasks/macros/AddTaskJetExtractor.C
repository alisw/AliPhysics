// AddTaskJetExtractor.C

AliAnalysisTaskJetExtractor* AddTaskJetExtractor(
  const char *trackArray         = "tracks",
  const char *jetArray           = "jets",
  const char *rhoObject          = "Rho",
  Double_t    jetRadius          = 0.3,
  Double_t    minJetEta          = 0.6,
  Double_t    minJetPt           = 0.15,
  Double_t    minTrackPt         = 0.15,
  Double_t    minJetAreaPerc     = 0.557,
  const char *suffix             = ""
)
{  
  cout << " ############ MACRO EXECUTION STARTED: AddTaskJetExtractor.C ############\n";
  //==============================================================================
  // Prepare analysis manager, containers, etc.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    ::Error("AddTaskJetExtractor", "No analysis manager to connect to.");
    return NULL;
  }  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetExtractor", "This task requires an input event handler");
    return NULL;
  }
  
  TString name("AliAnalysisTaskJetExtractor");
  if (strcmp(jetArray,"")) {
    name += "_";
    name += jetArray;
  }
  if (strcmp(rhoObject,"")) {
    name += "_";
    name += rhoObject;
  }
  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisDataContainer* contHistos = mgr->CreateContainer(Form("%s_histos", name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJetsHadronCF", AliAnalysisManager::GetCommonFileName()));

  //==============================================================================
  // Adding and configuring tasks

  AliAnalysisTaskJetExtractor* jetTask = new AliAnalysisTaskJetExtractor(name);
  jetTask->SetNeedEmcalGeom(kFALSE);
  jetTask->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont = 0;
  if(!strcmp(trackArray,"mctracks") || !strcmp(trackArray, "mcparticles"))
    trackCont = jetTask->AddMCParticleContainer(trackArray);
  else
    trackCont = jetTask->AddTrackContainer(trackArray);

  trackCont->SetParticlePtCut(minTrackPt);

  AliRDHFJetsCutsVertex* cuts = new AliRDHFJetsCutsVertex("jetCuts");

  // vertexing
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8, 0.8);
  esdTrackCuts->SetPtRange(1.0, 1.e10);

  cuts->AddTrackCuts(esdTrackCuts);
  cuts->SetNprongs(3);
  cuts->SetMinPtHardestTrack(1.0);//default 0.3
  cuts->SetSecVtxWithKF(kFALSE);//default with StrLinMinDist
  cuts->SetImpParCut(0.);//default 0
  cuts->SetDistPrimSec(0.);//default 0
  cuts->SetCospCut(-1);//default -1

  jetTask->SetVertexerCuts(cuts);

  AliJetContainer *jetCont = jetTask->AddJetContainer(jetArray,6,jetRadius);
  if (jetCont) {
    jetCont->SetRhoName(rhoObject);
    jetCont->SetPercAreaCut(minJetAreaPerc);
    jetCont->SetJetPtCut(minJetPt);
    jetCont->SetLeadingHadronType(0);
    jetCont->SetPtBiasJetTrack(minTrackPt);
    jetCont->SetJetEtaLimits(-minJetEta, +minJetEta);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->SetMaxTrackPt(1000);
  }

  mgr->AddTask(jetTask);

  //==============================================================================
  // Finalization

  mgr->ConnectInput  (jetTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (jetTask, 1, contHistos );
 
  cout << " ############ MACRO EXECUTION DONE: AddTaskJetExtractor.C ############\n";
 
  return jetTask;
}
