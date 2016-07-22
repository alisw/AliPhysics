void InitHistograms(AliDielectron *die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);

TString names=("NANO");
enum { kNANO=0 };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD     = kTRUE;
TString periodLHC = "";

AliDielectron* ConfigLMEE_nano_pp(Int_t cutDefinition, Bool_t hasMC=kFALSE, TString period="", Bool_t useTrackCuts = kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //

  periodLHC = period;
  printf("this is -%s- filtering \n",periodLHC.Data());

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  printf(" Adding %s%s config %s for %s\n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),periodLHC.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("%s",name.Data()));
  die->SetHasMC(hasMC);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,cutDefinition);
  if(useTrackCuts)
    SetupTrackCuts(die,cutDefinition);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // prefilter settings
  die->SetNoPairing();

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // histogram setup
  InitHistograms(die,cutDefinition);

  // cut QA
  die->SetCutQA();

  return die;
}

//______________________________________________________________________________________
void SetupEventCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the event cuts
  //

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","eventCuts");
  if(!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,+10.);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  Bool_t hasMC=die->GetHasMC();

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  //varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 1e30);
  //varCuts->AddCut(AliDielectronVarManager::kEta,         -1.0, 1.0);
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 5.0);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8, 0.8);
 
  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  //  AliDielectronVarCuts *pidVarCuts = new AliDielectronVarCuts("varPIDCuts","varPIDCuts");
  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");

  // TPC inclusion
  //pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-5,5);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3,3);
 
  if(!isESD) cuts->AddCut(trkFilter);
  cuts->AddCut(varCuts);
  cuts->AddCut(trkCuts);
  cuts->AddCut(pidCuts);
  cuts->Print();

}


//______________________________________________________________________________________
void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //
  Bool_t hasMC=die->GetHasMC();

  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());

  //add histograms to event class
  histos->AddClass("Event");
  histos->UserHistogram("Event","","", 150, 0.0, 150.0,   AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","","", 300,-15., +15.0,   AliDielectronVarManager::kZvPrim);
  // candidates monitoring
  histos->UserProfile("Event","","", AliDielectronVarManager::kTracks, 150,0,150., AliDielectronVarManager::kRefMultTPConly);
  histos->AddClass("Track");
  histos->UserHistogram("Track","","", 400,0.2,20.,200,0.,200., AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,kTRUE);

  die->SetHistogramManager(histos);
}
