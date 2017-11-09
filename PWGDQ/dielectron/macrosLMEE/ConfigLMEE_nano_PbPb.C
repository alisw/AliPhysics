
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

Bool_t  isESD     = kTRUE;
TString periodLHC = "";
TString list  = gSystem->Getenv("LIST");


AliDielectron* ConfigLMEE_nano_PbPb(Int_t cutDefinition, Bool_t hasMC=kFALSE, TString period="")
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
  printf(" Adding %s%s config %s for %s \n",(isESD?"ESD":"AOD"),(hasMC?" MC":""),name.Data(),list.Data());

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("%s",name.Data()));
  die->SetHasMC(hasMC);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  SetupEventCuts(die,cutDefinition);
  SetupTrackCuts(die,cutDefinition);
  // SetupV0Cuts(die,cutDefinition);   // switch off for nanoAODs??
  // SetupPairCuts(die,cutDefinition);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MISC vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // PID eta correction
  //SetEtaCorrection(die); //no eta corrction
  // prefilter settings
  //  if(hasMC)
  die->SetNoPairing();
  //  else die->SetPreFilterUnlikeOnly();
  //die->SetPreFilterAllSigns();

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
  // eventCuts->SetCentralityRange(10., 50., kTRUE);
  eventCuts->Print();
  die->GetEventFilter().AddCuts(eventCuts);

  AliDielectronVarCuts *eventplaneCuts = new AliDielectronVarCuts("eventplaneCuts","eventplaneCuts");
  eventplaneCuts->AddCut(AliDielectronVarManager::kQnTPCrpH2,-999.,kTRUE); // makes sure that the event has an eventplane
  eventplaneCuts->Print();


  die->GetEventFilter().AddCuts(eventplaneCuts);
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

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv FILTER CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqualSPDany); // I think we loose the possibility to use prefilter?

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv TRACK CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *varCuts   = new AliDielectronVarCuts("VarCuts","VarCuts");
  AliDielectronTrackCuts *trkCuts = new AliDielectronTrackCuts("TrkCuts","TrkCuts");
  // specific cuts
  trkCuts->SetITSclusterCut(AliDielectronTrackCuts::kOneOf, 3); // SPD any
  trkCuts->SetRequireITSRefit(kTRUE);
  trkCuts->SetRequireTPCRefit(kTRUE); // not useful when using prefilter

  // standard cuts
  varCuts->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 160.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 100.0);
  varCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,    0.0,   15.0);
  varCuts->AddCut(AliDielectronVarManager::kNclsSITS,     0.0,   3.1); // means 0 and 1 shared Cluster
  varCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,    0.0,   8.0);
  varCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,    80.0, 160.0);

  varCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 8.);
  varCuts->AddCut(AliDielectronVarManager::kEta,         -0.8,   0.8);
  varCuts->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  varCuts->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  varCuts->AddCut(AliDielectronVarManager::kKinkIndex0,   0.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  AliDielectronPID *pidCuts        = new AliDielectronPID("PIDCuts","PIDCuts");
  // TOF
  // pidCuts->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,-4.,4.,0.,0.,kFALSE, AliDielectronPID::kIfAvailable);
  // TPC
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,-4.,4.);
  pidCuts->AddCut(AliDielectronPID::kTPC,AliPID::kPion,-100.,3.5,0.,0.,kTRUE);
  // ITS
  pidCuts->AddCut(AliDielectronPID::kITS,AliPID::kElectron,-4.,4.);

  /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MC PID CUTS vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
  AliDielectronVarCuts *pidCutsMC = new AliDielectronVarCuts("PIDCutsMC","PIDCutsMC");
  pidCutsMC->SetCutType(AliDielectronVarCuts::kAny);
  pidCutsMC->SetCutOnMCtruth(kTRUE);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, -11., -11.);
  pidCutsMC->AddCut(AliDielectronVarManager::kPdgCode, +11., +11.);

  // activate the cut sets (order might be CPU timewise important)
  if(hasMC) {
    //    cuts->AddCut(pidCutsMC);
    if(!isESD) cuts->AddCut(trkFilter);
    cuts->AddCut(varCuts);
    cuts->AddCut(trkCuts);
    cuts->AddCut(pidCuts);
  }
  else {
    if(!isESD) cuts->AddCut(trkFilter);
    cuts->AddCut(varCuts);
    cuts->AddCut(trkCuts);
    cuts->AddCut(pidCuts);
  }
  cuts->Print();
  die->GetTrackFilter().AddCuts(cuts);

}

// //______________________________________________________________________________________
// void SetupV0Cuts(AliDielectron *die, Int_t cutDefinition)
// {
//   //
//   // Setup the V0 cuts
//   //
//
//   // Quality cuts (add the gamma filter to the cut group)
//   TIter next(die->GetTrackFilter().GetCuts());
//   AliAnalysisCuts *cuts;
//   while((cuts = (AliAnalysisCuts*)next())) {
//     if(cuts->IsA() == AliDielectronCutGroup::Class())  break;
//   }
//
//   AliDielectronV0Cuts *gammaV0Cuts = new AliDielectronV0Cuts("IsGamma","IsGamma");
//   gammaV0Cuts->SetPdgCodes(22,11,11);
//   gammaV0Cuts->SetDefaultPID(16);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.0, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.0, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kR,                             3.0,  90.0, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kM,                             0.0,   0.05, kFALSE);
//   //  gammaV0Cuts->AddCut(AliDielectronVarManager::kOpeningAngle,              0.0,   0.1, kFALSE);
//   gammaV0Cuts->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
//   //  gammaV0Cuts->AddCut(AliDielectronVarManager::kArmAlpha,                     -0.35,  0.35, kFALSE); // not sure if it works as expected
//   gammaV0Cuts->SetExcludeTracks(kTRUE);
//   gammaV0Cuts->Print();
//
//  //  const Double_t |cutAlphaG| < 0.35; &&  const Double_t cutQTG < 0.05;
//  //  const Double_t |cutAlphaG2|[2] = {0.6, 0.8}; &&  const Double_t cutQTG2 < 0.04;
//
//   if(cuts)
//     ((AliDielectronCutGroup*)cuts)->AddCut(gammaV0Cuts);
//   else
//     die->GetTrackFilter().AddCuts(gammaV0Cuts);
// }
//
// //______________________________________________________________________________________
// void SetupPairCuts(AliDielectron *die, Int_t cutDefinition)
// {
//   //
//   // Setup the pair cuts
//   //
//
//   // conversion rejection
//   Double_t gCut=0.05;
//   AliDielectronVarCuts *gammaCuts = new AliDielectronVarCuts("GammaCuts","GammaCuts");
//   gammaCuts->AddCut(AliDielectronVarManager::kM,            0.0,   gCut);
//   //  die->GetPairPreFilter().AddCuts(gammaCuts);
//
//   // rapidity selection
//   Double_t yCut=0.9;
//   AliDielectronVarCuts *rapCut=new AliDielectronVarCuts(Form("|Y|<%.1f",yCut),Form("|Y|<%.1f",yCut));
//   rapCut->AddCut(AliDielectronVarManager::kY,-1.*yCut,yCut);
//   die->GetPairFilter().AddCuts(rapCut);
//
// }

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
  histos->UserHistogram("Event","","", 100, 0.0, 100.0,   AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","","", 300,-15., +15.0,   AliDielectronVarManager::kZvPrim);
  histos->AddClass("Track_ev1+");
  histos->UserHistogram("Track_ev1+","","", 1000, 0., 10.,   AliDielectronVarManager::kPt);
  histos->AddClass("Track_ev1-");
  histos->UserHistogram("Track_ev1-","","", 1000, 0., 10.,   AliDielectronVarManager::kPt);
  // candidates monitoring
  histos->UserProfile("Event","","", AliDielectronVarManager::kTracks, 100,0,100., AliDielectronVarManager::kCentralityNew);
  histos->UserProfile("Event","","", AliDielectronVarManager::kPairs,  100,0,100., AliDielectronVarManager::kCentralityNew);

  die->SetHistogramManager(histos);
}
//
// //______________________________________________________________________________________
// void SetEtaCorrection(AliDielectron *die) {
//
//   if (AliDielectronPID::GetCentroidCorrFunction()) return;
//
//   TF2 *fCntrdCorr=0x0;
//   TF1 *fWdthCorr=0x0;
//   /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv DATA vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
//   if( !die->GetHasMC() ) {
//     // 2-dimensional eta correction for the centroid of electron sigmas
//     fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
// 			 //0.0, 3000.0, -0.9, +0.9); // Nacc
// 			 0.0,   90.0, -0.9, +0.9); // centrality
//     //  fCntrdCorr->SetParameters(0.723106, 0.23958, -6.31221, -0.687976, 15.912, 0.579609, -11.6901, -0.000354381); // Nacc
//     fCntrdCorr->SetParameters(+0.149002, +0.214644 , -6.034930, -0.529588, +14.97902, +0.402640, -10.890027, +0.011248); // centrality
//     // 1-dimensional eta correction for the width of electron sigmas
//     // fWdthCorr = new TF1("fWdthCorr", "pol2", 0.0, 3000.0);     // Nacc
//     // fWdthCorr->SetParameters(1.06108, 0.000217804,-5.80291e-08);
//     fWdthCorr = new TF1("fWdthCorr", "pol2", 0.0, 90.0);       // centrality
//     fWdthCorr->SetParameters(+1.290755, -0.005261, +0.000021);
//   }
//   else  {
//     /* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv MONTE CARLO vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
//     // 2-dimensional eta correction for the centroid of electron sigmas
//     fCntrdCorr = new TF2("fCntrdCorr", "[0] + [1]*y + [2]*y*y + [3]*TMath::Power(y,3) + [4]*TMath::Power(y,4) + [5]*TMath::Power(y,5) + [6]*TMath::Power(y,6) + [7]*x",
// 			      0.0, 3000.0, -0.9, +0.9);
//     fCntrdCorr->SetParameters(+0.378611, -0.070831, -3.076778, +0.121977, +8.576097, +0.113009, -5.001368, -0.000181);
//
//     // 1-dimensional eta correction for the width of electron sigmas
//     fWdthCorr = new TF1("fWdthCorr", "pol1", 0.0, 3000.0);
//     fWdthCorr->SetParameters(+0.881894, +0.000053);
//   }
//
//   // apply corrections
//   // AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kNacc,AliDielectronVarManager::kEta);
//   // AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kNacc);
//   AliDielectronPID::SetCentroidCorrFunction(fCntrdCorr,AliDielectronVarManager::kCentralityNew,AliDielectronVarManager::kEta);
//   AliDielectronPID::SetWidthCorrFunction(fWdthCorr,AliDielectronVarManager::kCentralityNew);
//
// }
