// AddTaskChargedJetsHadronCF.C

AliAnalysisTaskChargedJetsHadronCF* AddTaskChargedJetsHadronCF(
  UInt_t physSel                 = 0,
  const char *trackArray         = "tracks",
  const char *jetArray           = "jets",
  const char *rhoObject          = "Rho",
  Double_t    jetRadius          = 0.2,
  Double_t    minJetEta          = 0.7,
  Double_t    minLeadingHadronPt = 0.0,
  Double_t    minJetPt           = 0.15,
  Double_t    minTrackPt         = 0.15,
  Double_t    minJetAreaPerc     = 0.557,
  const char *suffix             = ""
)
{  
  cout << " ############ MACRO EXECUTION STARTED: AddTaskChargedJetsHadronCF.C ############\n";
  //==============================================================================
  // Prepare analysis manager, containers, etc.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr)
  {
    ::Error("AddTaskChargedJetsHadronCF", "No analysis manager to connect to.");
    return NULL;
  }  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskChargedJetsHadronCF", "This task requires an input event handler");
    return NULL;
  }
  
  TString name("AliAnalysisTaskChargedJetsHadronCF");
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

  AliAnalysisTaskChargedJetsHadronCF* jetTask = new AliAnalysisTaskChargedJetsHadronCF(name);
  jetTask->SetNeedEmcalGeom(kFALSE);
  jetTask->SetVzRange(-10.,10.);
  jetTask->SetOffTrigger(physSel);
  
  AliTrackContainer *trackCont = jetTask->AddTrackContainer(trackArray);
  trackCont->SetParticlePtCut(minTrackPt);

  AliJetContainer *jetCont = jetTask->AddJetContainer(jetArray,6,jetRadius);
  if (jetCont) {
    jetCont->SetRhoName(rhoObject);
    jetCont->SetPercAreaCut(minJetAreaPerc);
    jetCont->SetJetPtCut(minJetPt);
    jetCont->SetLeadingHadronType(0);
    jetCont->SetPtBiasJetTrack(minLeadingHadronPt);
    jetCont->SetJetEtaLimits(-minJetEta, +minJetEta);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->SetMaxTrackPt(1000);
  }

  mgr->AddTask(jetTask);

  //==============================================================================
  // Finalization

  mgr->ConnectInput  (jetTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (jetTask, 1, contHistos );
 
  cout << " ############ MACRO EXECUTION DONE: AddTaskChargedJetsHadronCF.C ############\n";
 
  return jetTask;
}
