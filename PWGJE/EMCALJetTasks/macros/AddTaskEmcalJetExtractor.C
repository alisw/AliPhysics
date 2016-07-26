// AddTaskEmcalJetExtractor.C

AliAnalysisTaskEmcalJetExtractor* AddTaskEmcalJetExtractor(
  UInt_t physSel                 = 0,
  const char* trackArray         = "tracks",
  const char* jetArray           = "jets",
  const char* rhoObject          = "Rho",
  Double_t    jetRadius          = 0.2,
  const char* jetTreeName,
  Int_t       extractionType,
  Int_t       extractionCriterium,
  Double_t    extractionMinPt,
  Double_t    extractionMaxPt,
  Double_t    extractionPercentage
)
{  
  cout << " ############ MACRO EXECUTION STARTED: AddTaskEmcalJetExtractor.C ############\n";
  //==============================================================================
  // Prepare analysis manager, containers, etc.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetExtractor", "No analysis manager to connect to.");
    return NULL;
  }  
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetExtractor", "This task requires an input event handler");
    return NULL;
  }
  
  // ####### CONTAINER NAME

  TString name(jetTreeName);

  if(extractionType == 0)
    name += "_AliEmcalJet";
  else if(extractionType == 1)
    name += "_AliBasicJet";
  else if(extractionType == 2)
    name += "_AliBasicJetWithConstituents";

  if(extractionCriterium == 0)
    name += "_CriteriumMinBias";
  else if(extractionCriterium == 1)
    name += "_CriteriumSignal";
  else if(extractionCriterium == 2)
    name += "_CriteriumBackground";

  name += "_Pt";
  name += Form("%3.2f-%3.2f", extractionMinPt, extractionMaxPt);

  name += "_Percentage";
  name += Form("%e", extractionPercentage);

  name += "_";
  name += jetArray;

  if (strcmp(rhoObject,"")) {
    name += "_";
    name += rhoObject;
  }

  AliAnalysisDataContainer* contHistos = mgr->CreateContainer(Form("%s", name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:JetExtractor", AliAnalysisManager::GetCommonFileName()));

  //==============================================================================
  // Adding and configuring tasks
  Double_t    minJetEta          = 0.7;
  Double_t    minLeadingHadronPt = 0.0;
  Double_t    minJetPt           = 0.15;
  Double_t    minTrackPt         = 0.15;
  Double_t    minJetAreaPerc     = 0.557;

  AliAnalysisTaskEmcalJetExtractor* jetTask = new AliAnalysisTaskEmcalJetExtractor(name);
  jetTask->SetNeedEmcalGeom(kFALSE);
  jetTask->SetVzRange(-10.,10.);
//  jetTask->SetOffTrigger(physSel);
  jetTask->DefineExtraction(extractionType, extractionCriterium, extractionMinPt, extractionMaxPt, extractionPercentage);

  AliParticleContainer *trackCont = jetTask->AddTrackContainer(trackArray);
  trackCont->SetParticlePtCut(minTrackPt);

  AliJetContainer *jetCont = jetTask->AddJetContainer(jetArray,"USER",jetRadius);
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
 
  cout << " ############ MACRO EXECUTION DONE: AddTaskEmcalJetExtractor.C ############\n";
 
  return jetTask;
}
