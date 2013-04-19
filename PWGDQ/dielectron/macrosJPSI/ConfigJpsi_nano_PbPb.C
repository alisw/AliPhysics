void InitHistograms(AliDielectron *die, Int_t cutDefinition);

void SetupEventCuts(AliDielectron *die, Int_t cutDefinition);
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition);
void SetupV0Cuts(   AliDielectron *die, Int_t cutDefinition);
void SetupPairCuts( AliDielectron *die, Int_t cutDefinition);

void ConfigEvtPlane(AliDielectron *die, Int_t cutDefinition);
void SetEtaCorrection(AliDielectron *die);

TString names=("NANO");
enum { kNANO=0 };

TObjArray *arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntries();

Bool_t  isESD = kTRUE;
TString list  = gSystem->Getenv("LIST");

AliDielectron* ConfigJpsi_nano_PbPb(Int_t cutDefinition, Bool_t hasMC=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //

  //ESD handler?
  isESD=(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->IsA()==AliESDInputHandler::Class());

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();
  printf(" Adding %s%s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),list.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("%s",name.Data()));
  die->SetHasMC(hasMC);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  SetupV0Cuts(die,cutDefinition);
  SetupPairCuts(die,cutDefinition);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // PID eta correction
  SetEtaCorrection(die);
  // prefilter settings
  die->SetPreFilterUnlikeOnly();
  //die->SetPreFilterAllSigns();
  //die->SetNoPairing();

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv OUTPUT vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // histogram setup
  //  InitHistograms(die,cutDefinition);

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
  eventCuts->SetCentralityRange(0., 90.);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

}

//______________________________________________________________________________________
void SetupTrackCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the track cuts
  //

  // Quality cuts
  AliDielectronCutGroup* cuts = new AliDielectronCutGroup("cuts","cuts",AliDielectronCutGroup::kCompAND);
  die->GetTrackFilter().AddCuts(cuts);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  //  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  //  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  //  trkFilter->SetMinNCrossedRowsOverFindable(0.6);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  varCuts->AddCut(AliDielectronVarManager::kPt,           0.85, 1e30);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.9,   0.9);
  //  varCuts->AddCut(AliDielectronVarManager::kTPCclsSegments,7.,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 160.0);
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   4.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.0);
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  //  AliDielectronVarCuts *pidVarCuts = new AliDielectronVarCuts("varPIDCuts","varPIDCuts");
  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TOF inclusion
  //  pidVarCuts->AddCut(AliDielectronVarManager::kTOFbeta,      0.2,   0.9, kTRUE);
  pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-3,3.,0.,0.,kFALSE,
		  AliDielectronPID::kIfAvailable);
  // TPC inclusion
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-3.,3.); // when eta correction ON
  //  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,4.0,0.,0.,kTRUE);
  //  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kProton,-100.,3.5,0.,0.,kTRUE);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TENDER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // exclude conversion electrons selected by the tender
  AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","noConv");
  noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);

  // activate the cut sets (order might be CPU timewise important)
  //if(!isESD) cuts->AddCut(trkFilter);
  cuts->AddCut(varCuts);
  cuts->AddCut(trkCuts);
  cuts->AddCut(pidCuts);
  //cuts->AddCut(pidVarCuts);
  //cuts->AddCut(noconv);
  cuts->Print();

}

//______________________________________________________________________________________
void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the V0 cuts
  //

  // Quality cuts (add the gamma filter to the cut group)
  TIter next(die->GetTrackFilter().GetCuts());
  AliAnalysisCuts *cuts;
  while((cuts = (AliAnalysisCuts*)next())) {
    if(cuts->IsA() == AliDielectronCutGroup::Class())  break;
  }

  AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
  gammaV0Cuts->SetPdgCodes(22,11,11);
  gammaV0Cuts->SetDefaultPID(16);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
  gammaV0Cuts->SetExcludeTracks(kTRUE);
  gammaV0Cuts->Print();

 //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
 //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;

  if(cuts)
    ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
  else
    die->GetTrackFilter().AddCuts(gammaV0Cuts);
}

//______________________________________________________________________________________
void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Setup the pair cuts
  //

  // conversion rejection
  Double_t gCut=0.05;
  AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
  gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
  die->GetPairPreFilter().AddCuts(gammaCuts);

  // rapidity selection
  Double_t yCut=0.9;
  AliDielectronVarCuts *rapCut=new AliDielectronVarCuts(Form("|Y|<%.1f",yCut),Form("|Y|<%.1f",yCut));
  rapCut->AddCut(AliDielectronVarManager::kY,-1.*yCut,yCut);
  die->GetPairFilter().AddCuts(rapCut);

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
  histos->UserHistogram("Event","","", 100, 0.0, 100.0,   AliDielectronVarManager::kCentrality);
  die->SetHistogramManager(histos);
}

//______________________________________________________________________________________
void SetEtaCorrection(AliDielectron *die) {

  if (AliDielectronPID::GetCentroidCorrFunction()) return;

  TF2 *fCntrdCorr=0x0;
  TF1 *fWdthCorr=0x0;
  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv DATA vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  if( !die->GetHasMC() ) {
    // 2-dimensional eta correction for the centroid of electron sigmas
    fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
			      0.0, 3000.0, -0.9, +0.9);
    fCntrdCorr->SetParameters(0.723106, 0.23958, -6.31221, -0.687976, 15.912, 0.579609, -11.6901, -0.000354381);
    // 1-dimensional eta correction for the width of electron sigmas
    fWdthCorr = new TF1("fWdthCorr", "pol2", 0.0, 3000.0);
    fWdthCorr->SetParameters(1.06108, 0.000217804,-5.80291e-08);
  }
  else  {
    /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MONTE CARLO vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    // 2-dimensional eta correction for the centroid of electron sigmas
    fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
			      0.0, 3000.0, -0.9, +0.9);
    fCntrdCorr->SetParameters(+0.378611, -0.070831, -3.076778, +0.121977, +8.576097, +0.113009, -5.001368, -0.000181);
    // 1-dimensional eta correction for the width of electron sigmas
    fWdthCorr = new TF1("fWdthCorr", "pol1", 0.0, 3000.0);
    fWdthCorr->SetParameters(+0.881894, +0.000053);
  }

  // apply corrections
  AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kNacc,AliDielectronVarManager::kEta);
  AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kNacc);

}
