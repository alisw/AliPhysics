#include "THn.h"


void      InitHistograms(AliDielectron *die);
void      InitCF(AliDielectron* die);
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx, kMee_fine};

TString names=("cut5_pt200");

// For eta correction
// TString names=("kPbPb2015_pure_electron_pt200;kPbPb2015_pure_pion_pt200");
// TString names=("kPbPb2015_pure_electron_pt200_MC;kPbPb2015_pure_pion_pt200_MC");

// no track cuts
// TString names=("kPbPb2015_noPID_NoTrackQualityCuts_Pt400");

// MC TOF efficiency
// TString names=("kPbPb2015_pure_electron_pt200_woTOF;kPbPb2015_pure_electron_pt200_wTOFhit;kPbPb2015_pure_electron_pt200_wTOFreq;kPbPb2015_pure_electron_pt200_wTOFif;kPbPb2015_PDG_pure_electron_pt200_woTOF;kPbPb2015_PDG_pure_electron_pt200_wTOFhit;kPbPb2015_PDG_pure_electron_pt200_wTOFreq;kPbPb2015_PDG_pure_electron_pt200_wTOFif"); // For evaluating TOF matching efficiency
// Data TOF efficiency
// TString names=("kPbPb2015_pure_electron_pt200_woTOF;kPbPb2015_pure_electron_pt200_wTOFhit;kPbPb2015_pure_electron_pt200_wTOFreq;kPbPb2015_pure_electron_pt200_wTOFif"); // For evaluating TOF matching efficiency

// Contamination Study
// TString names=("kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej;kPbPb2015_pidV0_electron_pt400_woPionRej;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut;kPbPb2015_pidV0_electron_pt400_wITScut;kPbPb2015_pure_pion_pt400_wITScut;kPbPb2015_pure_kaon_pt400_wITScut;kPbPb2015_pure_proton_pt400_wITScut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut;kPbPb2015_pidV0_electron_pt400_wTPCcut;kPbPb2015_pure_pion_pt400_wTPCcut;kPbPb2015_pure_kaon_pt400_wTPCcut;kPbPb2015_pure_proton_pt400_wTPCcut");
// TString names=("kPbPb2015_Pt200_cut5_woTPCelecut;kPbPb2015_pure_electron_pt200_woTPCelecut;kPbPb2015_pure_kaon_pt200_woTPCelecut;kPbPb2015_pure_proton_pt200_woTPCelecut;kPbPb2015_Pt200_cut5_woPionRej;kPbPb2015_pure_pion_pt200_wTPCelecut");


TObjArray *arrNames = names.Tokenize(";");
const Int_t nDie = arrNames->GetEntries();


Bool_t SetTPCCorrection = kFALSE;
Bool_t SetITSCorrection = kFALSE;
Bool_t SetTOFCorrection = kFALSE;
Bool_t SetMixing = kFALSE;
Bool_t SetPairing = kTRUE;
Bool_t SetMCSignal = kFALSE;
Bool_t SetFillPureMCStep = kFALSE;
Bool_t SetAnalysisRun = kFALSE; // Set to true if you only need the basic analysis histograms


// Set standard cuts
void StandardCuts(AnalysisCut* AnaCut){
  AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);
  AnaCut.SetPreFilterType(LMEECutLib::kNoPreFilter);
  AnaCut.SetPIDPre(LMEECutLib::kStandardPre);
  AnaCut.SetTrackSelectionPre(LMEECutLib::kPrefilter_cut1);
  AnaCut.SetPairCutsPre(LMEECutLib::kNoPairCutsPre);

  AnaCut.SetMixing(LMEECutLib::kEventMixing_1);
  AnaCut.SetESDTrackSelection(LMEECutLib::kStandardESD);
  return;
}

AliDielectron* Config_caklein_LMEEPbPb_AOD_2018(TString cutDefinition, Bool_t hasMC = kFALSE, Int_t centrality = 0)
{

  Int_t  Centrality = -1;
  if      (centrality == 0) Centrality = LMEECutLib::kPbPb0090;
  else if (centrality == 1) Centrality = LMEECutLib::kPbPb0010;
  else if (centrality == 2) Centrality = LMEECutLib::kPbPb1050;
  else                      Centrality = LMEECutLib::kPbPb0080;

  //
  // Setup the instance of AliDielectron
  //
  AnalysisCut AnaCut;

  LMEECutLib*  LMcutlib = new LMEECutLib();

  // task name
  TString name = cutDefinition;

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));

  std::cout << "##### Cutsetting #####" << std::endl;
  if(SetTPCCorrection) LMcutlib->SetEtaCorrectionTPC(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE);
  if(SetITSCorrection) LMcutlib->SetEtaCorrectionITS(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE);
  if(SetTOFCorrection) LMcutlib->SetEtaCorrectionTOF(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE);
  // die->SetCutQA();

  cout << "cutDefinition = " << cutDefinition << endl;
  // Setup Analysis Selection

  // ############ ANALYSIS CUTS
  if (cutDefinition == "cut1_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut20_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_20_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_20);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut21_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_21_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_21);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut22_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_22_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_22);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut23_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_23_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_23);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut24_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_24_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_24);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut25_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_25_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_25);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut26_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_26_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_26);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut27_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_27_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_27);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut28_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_28_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_28);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut29_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_29_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_29);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut30_pt200"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_30_pt200);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_30);
    AnaCut.SetCentrality(Centrality);
    AnaCut.SetStandardCut();
  }

  else {
      cout << " =============================== " << endl;
      cout << " ==== INVALID CONFIGURATION ==== " << endl;
      cout << " cutDefinition = " << cutDefinition << endl;
      cout << " =============================== " << endl;
  }



  // __________________________________________________
  // POSSIBLE FURTHER SETTINGS
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  if(SetPairing != kTRUE){
    die->SetNoPairing();
  }

  // Now configure task
  //
  // add centrality selection to event cuts
  die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(AnaCut) );
  // switch off KF Particle
  die->SetUseKF(kFALSE);

  // --------------------------------------------------
  // with Rejection Step (Prefilter)
  // --------------------------------------------------
  if (AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterAllSigns || AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterUnlikeOnly)
  {
    // set initial track filter.
    // the function 'GetPIDCutsPre()' must also call 'GetTrackCutsPre()'!
    die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsPre(AnaCut) );

    // set Prefilter. "remove all tracks from the Single track arrays that pass the cuts in this filter" (comment in AliDielectron.cxx)
    // cuts = REJECTION!!!
    die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(AnaCut) );

    // "apply leg cuts after the pre filter" (comment in AliDielectron.cxx)
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
    die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );

    // set Pairfilter.
    // cuts = SELECTION!!!
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(AnaCut, kFALSE) );

  }
  // --------------------------------------------------
  // without Rejection Step
  // --------------------------------------------------
  else {
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
    die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(AnaCut, kFALSE) );
  }
  // -------------------------------------------------


  if (SetMixing == kTRUE){
    std::cout << "MIXING ACTIVATED" << std::endl;
    AliDielectronMixingHandler* mix = LMcutlib->GetMixingHandler(AnaCut);
    die->SetMixingHandler(mix);
  }

  if (SetMCSignal) SetupMCsignals(die);

  // histogram setup
  InitHistograms(die);

  return die;
}

//______________________________________________________________________________________

void InitHistograms(AliDielectron *die)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
  AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(),die->GetTitle());

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair;Track_Legs;Pre;RejPair;RejTrack");

  //Event class
  histos->AddClass("Event"); // all classes will be stored in 'THashList fHistoList'

  //Track classes
  //to fill also track info from 2nd event loop until 3
  // in AliDielectron.cxx: fgkTrackClassNames[4] = {"ev1+","ev1-","ev2+","ev2-"};
  for (Int_t i=0; i<2; ++i){
    histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
  }

  //Pair classes
  // to fill also mixed event histograms loop until 10
  // fgkPairClassNames[11] = {
  //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
  //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
  //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
  //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
  // };
  if (SetPairing){
    for (Int_t i=0; i<3; ++i){
      histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
      // Legs of final Pairs. Both charges together. No duplicate entries.
      // histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
    }
  }

  //ME and track rot
  if (die->GetMixingHandler()) {
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
  }
  if (die->GetTrackRotator()) {
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(10)));
  }

  //PreFilter Classes
  //to fill also track info from 2nd event loop until 2
  for (Int_t i=0; i<2; ++i){
    // histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
  }

  //Create Classes for Rejected Tracks/Pairs:
  // for (Int_t i=0; i<3; ++i){
  //   histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
  //   // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
  //   histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
  // }

  //add MC signal histograms to pair class
  if(die->GetMCSignals()) {
    for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
      histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      if (SetFillPureMCStep) histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
    }
  }

  if(SetAnalysisRun == kTRUE){
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);

    // histos->UserHistogram("Pair","InvMass",";m_{ee} (GeV/c^{2})",
    //     AliDielectronHelper::MakeLinBinning(100,  0., 100.), AliDielectronVarManager::kM);
    histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event","centrality","",100,0.,100,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",160,0,8.,AliDielectronVarManager::kPt);
    // histos->UserHistogram("Track","Pt_high",";Pt [GeV];#tracks",80,20,100.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","Eta","",100,-1,1,AliDielectronVarManager::kEta);
    die->SetHistogramManager(histos);
    return;
  }
  else {
    // ################# EVENT ##########################
    // ##################################################
    histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event","centrality","",100,0.,100,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","zVertexPrim","",150,-15,15,AliDielectronVarManager::kZvPrim);

    // ################# TRACK ##########################
    // ##################################################
    histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",160,0,8.,AliDielectronVarManager::kPt);
    // histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)", GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);

    // ############## ITS
    // histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)", GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
    AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  -5,5.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    // histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
    //                       GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    // histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
    //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);

    // ############### TPC
    // histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
    //                       AliDielectronHelper::MakeLinBinning(500,  0.,5.), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
    histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    // histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
    // AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  -5,5.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    // histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
    //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    // histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
    //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
    // histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
    //                       GetVector(kEta3D), GetVector(kSigmaEle), AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);

    // ################## TOF
    histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
    // histos->UserHistogram("Track","TOFnSigmaPion_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
    //                       GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
    // ################## 2D-PID
    // histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
    //                       50,0.,2.5, 160,-12.,20., 150,-10.,20.,
    //                       AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
    // histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
    //                       50,0.,2.5, 160,-12.,20., 50,-5.,5.,
    //                       AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);


    // ############# Eta and Phi
    histos->UserHistogram("Track","Eta","",100,-1,1,AliDielectronVarManager::kEta);
    // histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    // histos->UserHistogram("Track","Phi_Pt","",200,0,10.,120,0.,TMath::TwoPi(),AliDielectronVarManager::kPt, AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);
    // histos->UserHistogram("Track","DEta_DPhi","",100,0,1.6, 100,-3.14159,3.14159,AliDielectronVarManager::kDeltaEta,AliDielectronVarManager::kDeltaPhi);

    // ############## DCA
    // histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
    // histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);

    // ############### Quality
    histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
    histos->UserHistogram("Track","NclsSITS",";ITS shared clusters;#tracks",7,-0.5,6.5.,AliDielectronVarManager::kNclsSITS);
    histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
    histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
    histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",101,59.5,161.5,AliDielectronVarManager::kNFclsTPCr);
    histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",50,0.7,1.2,AliDielectronVarManager::kNFclsTPCfCross);
    histos->UserHistogram("Track","Chi2TPCConstrainedVsGlobal",";golden #chi^{2};#tracks",101,-0.5,100.5,AliDielectronVarManager::kChi2TPCConstrainedVsGlobal);
    histos->UserHistogram("Track","TPCsignalN",";TPC cluster with dE/dx;#tracks",101,59.5,161.5,AliDielectronVarManager::kTPCsignalN);
    histos->UserHistogram("Track","TPCnCls_TPCsignalN",";TPC number clusters;TPC cluster with dE/dx;#tracks",101,59.5,161.5,101,59.5,161.5,AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kTPCsignalN);

    // Relevant for shared cluster cut study
    // histos->UserHistogram("Track","clsSITS_pt_centrality",";ITS shared clusters;Pt [GeV];Centrality",7,-0.5,6.5., 200,0,10., 10, 0, 100.,AliDielectronVarManager::kNclsSITS,AliDielectronVarManager::kPt,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Track","ITSnClsSMap","ITSnClsSMap; ITSnClsSMap ;#tracks",64,-0.5,63.5,AliDielectronVarManager::kNclsSMapITS);
    // histos->UserHistogram("Track","SDD1_SDD2","SDD1_SDD2;SDD1;SDD2;#tracks",100,0.,1000.,100,0.,1000., AliDielectronVarManager::kITSsignalSDD1, AliDielectronVarManager::kITSsignalSDD2);
    // histos->UserHistogram("Track","SSD1_SSD2","SSD1_SSD2;SSD1;SSD2;#tracks",100,0.,1000.,100,0.,1000., AliDielectronVarManager::kITSsignalSSD1, AliDielectronVarManager::kITSsignalSSD2);


    histos->UserHistogram("Track","pdgCode",";PDG Code;#tracks",5000+1,-2500.5,2500.5,AliDielectronVarManager::kPdgCode);
    histos->UserHistogram("Track","pdg_Mother",";PDG Code Mother;#tracks",2000+1,-1000.5,1000.5,AliDielectronVarManager::kPdgCodeMother);
    if (SetMCSignal){
      histos->UserHistogram("Track","TrackSource",";track_source;#tracks",64,-0.5, 63.5,AliDielectronVarManager::kMCLegSource);
      histos->UserHistogram("Track","TrackSource_PDGMother",";TrackSource;PDGMother",64,-0.5, 63.5,2001,-1000.5,1000.5,AliDielectronVarManager::kMCLegSource,AliDielectronVarManager::kPdgCodeMother);
    }

    // ################# PAIRS ##########################
    // ##################################################

    //2D and 3D histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee_fine), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    // histos->UserHistogram("Pair","DeltaPhiEP","#Delta#varphi_vs_TPC_EP;#Delta#varphi;#pairs",
    //                       100, -3.14159, 3.14159,
    //                       AliDielectronVarManager::kQnDeltaPhiTPCrpH2);
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);

    // ################# opening angle and PhiV
    // histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
    //                       GetVector(kMee), GetVector(kOpAng), AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);


    // #################Histograms for eta-map creation
    // const Int_t dimensions = 4;
    // Int_t bins[dimensions] = {100, 20, 40, 16};
    // UInt_t value[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kNTrk, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kEta};
    // Double_t xmin[dimensions] = {0., 0., -4, -0.8};
    // Double_t xmax[dimensions] = {10., 20000., 4, 0.8};
    // histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value);
    // //
    // const Int_t dimensions = 4;
    // Int_t bins[dimensions] = {100, 20, 40, 16};
    // UInt_t value[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kNTrk, AliDielectronVarManager::kTOFnSigmaEle, AliDielectronVarManager::kEta};
    // Double_t xmin[dimensions] = {0., 0., -4, -0.8};
    // Double_t xmax[dimensions] = {10., 20000., 4, 0.8};
    // histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value);
    //
    const Int_t dimensions    = 4;
    Int_t bins[dimensions]    = {100, 20,     50, 16  };
    Double_t xmin[dimensions] = {0.,  0.,     -5, -0.8};
    Double_t xmax[dimensions] = {10., 40000.,  5,  0.8};
    UInt_t value_TPCnSigmaEle[dimensions] = {AliDielectronVarManager::kP, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kEta};
    histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value_TPCnSigmaEle);

    die->SetHistogramManager(histos);
  }
}



TVectorD *GetVector(Int_t var)
{
  switch (var)
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(50, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
    case kMee_fine: return AliDielectronHelper::MakeLinBinning( 400, 0, 4);
    case kSigmaEle: return AliDielectronHelper::MakeLinBinning(50,-5.,5.);

    case kSigmaOther: return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);

    case kTPCdEdx: return AliDielectronHelper::MakeLinBinning(140,  0.,140.);

    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                                                   0.10, 0.12, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46,
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86,
                                                                   0.90, 0.94, 0.98, 1.02, 1.06,
                                                                   1.10, 1.20, 1.30, 1.40, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.80,
                                                                   2.90, 3.00, 3.05, 3.10, 3.20, 3.30, 3.50, 4.00, 4.50, 5.00
                                                                   ");

    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46,
                                                                   0.50
                                                                   ");

    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80,
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80,
                                                                   5.00, 6.00, 7.00, 8.00
                                                                   ");

    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45,
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95,
                                                                   2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45,
                                                                   2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95,
                                                                   3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,
                                                                   4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90,
                                                                   5.00, 5.25, 5.50, 5.75,
                                                                   6.00, 6.25, 6.50, 6.75,
                                                                   7.00, 7.5, 8.00
                                                                   ");

    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  }
  //if ( var.EqualTo("p_2D"      , kIgnoreCase) ) return AliDielectronHelper::MakeLinBinning(160,0.,8.);
}

TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //
  //  return vec;
}

void InitCF(AliDielectron* die)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPt, 80, 0.0, 8., kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta, 16, -0.8, 0.8, kTRUE);
  cf->AddVariable(AliDielectronVarManager::kM, 100, 0, 0.14);
  cf->AddVariable(AliDielectronVarManager::kPt, 8, 0., 8.);

  // cf->AddVariable(AliDielectronVarManager::kPIn,80,0.,4.,kTRUE);
  // cf->AddVariable(AliDielectronVarManager::kRefMultTPConly,"0.,100.,250.,500.,1000.,3000.",kTRUE);
  // cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,GetVector(kSigmaEle),kTRUE);
  // cf->AddVariable(AliDielectronVarManager::kEta,GetVector(kEta3D),kTRUE);

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}

void SetupMCsignals(AliDielectron* die){


  bool fFillPureMC = SetFillPureMCStep;

  AliDielectronSignalMC* single_ele_dontcare = new AliDielectronSignalMC("single_ele_dontcare","single_ele_dontcare");
  single_ele_dontcare->SetFillPureMCStep(fFillPureMC);
  single_ele_dontcare->SetLegPDGs(11,0);
  single_ele_dontcare->SetCheckBothChargesLegs(kTRUE,kTRUE);
  single_ele_dontcare->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC* single_ele_primary = new AliDielectronSignalMC("single_ele_primary","single_ele_primary");
  single_ele_primary->SetFillPureMCStep(fFillPureMC);
  single_ele_primary->SetLegPDGs(11,0);
  single_ele_primary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  single_ele_primary->SetLegSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC* single_ele_secondary = new AliDielectronSignalMC("single_ele_secondary","single_ele_secondary");
  single_ele_secondary->SetFillPureMCStep(fFillPureMC);
  single_ele_secondary->SetLegPDGs(11,0);
  single_ele_secondary->SetCheckBothChargesLegs(kTRUE,kTRUE);
  single_ele_secondary->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC* single_ele_finalstate = new AliDielectronSignalMC("single_ele_finalstate","single_ele_finalstate");
  single_ele_finalstate->SetFillPureMCStep(fFillPureMC);
  single_ele_finalstate->SetLegPDGs(11,0);
  single_ele_finalstate->SetCheckBothChargesLegs(kTRUE,kTRUE);
  single_ele_finalstate->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC* single_ele_finalstate_fromBG = new AliDielectronSignalMC("single_ele_finalstate_fromBG","single_ele_finalstate_fromBG");
  single_ele_finalstate_fromBG->SetFillPureMCStep(fFillPureMC);
  single_ele_finalstate_fromBG->SetLegPDGs(11,0);
  single_ele_finalstate_fromBG->SetCheckBothChargesLegs(kTRUE,kTRUE);
  single_ele_finalstate_fromBG->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kDontCare);


  AliDielectronSignalMC* pair_sameMother = new AliDielectronSignalMC("MCpair_sameMother","MCpair_sameMother");
  pair_sameMother->SetFillPureMCStep(fFillPureMC);
  pair_sameMother->SetLegPDGs(11,-11);
  pair_sameMother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // pair_sameMother->SetLegSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);
  pair_sameMother->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  pair_sameMother->SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameMother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. has no effect for final state ele.

  AliDielectronSignalMC* pair_diffMother = new AliDielectronSignalMC("MCpair_diffMother","MCpair_diffMother");
  pair_diffMother->SetFillPureMCStep(fFillPureMC);
  pair_diffMother->SetLegPDGs(11,-11);
  pair_diffMother->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffMother->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //pair_diffM_inclInj->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //mother
  pair_diffMother->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  //pair_diffMother->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. has no effect for final state ele.

  AliDielectronSignalMC* pair_diffM_oneGamma = new AliDielectronSignalMC("MCpair_diffM_oneGamma","MCpair_diffM_oneGamma");
  pair_diffM_oneGamma->SetFillPureMCStep(fFillPureMC);
  pair_diffM_oneGamma->SetLegPDGs(11,-11);
  pair_diffM_oneGamma->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffM_oneGamma->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  //mother
  pair_diffM_oneGamma->SetMothersRelation(AliDielectronSignalMC::kDifferent);

  AliDielectronSignalMC* pair_diffM_oneGamma_LSmm = new AliDielectronSignalMC("MCpair_diffM_oneGamma_LSmm","MCpair_diffM_oneGamma_LSmm");
  pair_diffM_oneGamma_LSmm->SetFillPureMCStep(fFillPureMC);
  pair_diffM_oneGamma_LSmm->SetLegPDGs(11,11);
  pair_diffM_oneGamma_LSmm->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffM_oneGamma_LSmm->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  pair_diffM_oneGamma_LSmm->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffM_oneGamma_LSmm->SetCheckLikeSign(kFALSE,kTRUE);
  pair_diffM_oneGamma_LSmm->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.

  AliDielectronSignalMC* pair_diffM_oneGamma_LSpp = new AliDielectronSignalMC("MCpair_diffM_oneGamma_LSpp","MCpair_diffM_oneGamma_LSpp");
  pair_diffM_oneGamma_LSpp->SetFillPureMCStep(fFillPureMC);
  pair_diffM_oneGamma_LSpp->SetLegPDGs(-11,-11);
  pair_diffM_oneGamma_LSpp->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffM_oneGamma_LSpp->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kSecondary);
  pair_diffM_oneGamma_LSpp->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffM_oneGamma_LSpp->SetCheckLikeSign(kTRUE,kFALSE);
  pair_diffM_oneGamma_LSpp->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.

  AliDielectronSignalMC* pair_sameGamma = new AliDielectronSignalMC("MCpair_sameGamma","MCpair_sameGamma");
  pair_sameGamma->SetFillPureMCStep(fFillPureMC);
  pair_sameGamma->SetLegPDGs(11,-11);
  pair_sameGamma->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_sameGamma->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
  //mother
  pair_sameGamma->SetMothersRelation(AliDielectronSignalMC::kSame);
  pair_sameGamma->SetMotherPDGs(22,22);

  AliDielectronSignalMC* pair_diffGamma = new AliDielectronSignalMC("MCpair_diffGamma","MCpair_diffGamma");
  pair_diffGamma->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma->SetLegPDGs(11,-11);
  pair_diffGamma->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffGamma->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
  //mother
  pair_diffGamma->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma->SetMotherPDGs(22,22);

  AliDielectronSignalMC* pair_diffGamma_LSmm = new AliDielectronSignalMC("MCpair_diffGamma_LSmm","MCpair_diffGamma_LSmm");
  pair_diffGamma_LSmm->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSmm->SetLegPDGs(11,11);
  pair_diffGamma_LSmm->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSmm->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSmm->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSmm->SetCheckLikeSign(kFALSE,kTRUE);
  pair_diffGamma_LSmm->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSmm->SetMotherPDGs(22,22);

  AliDielectronSignalMC* pair_diffGamma_LSpp = new AliDielectronSignalMC("MCpair_diffGamma_LSpp","MCpair_diffGamma_LSpp");
  pair_diffGamma_LSpp->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSpp->SetLegPDGs(-11,-11);
  pair_diffGamma_LSpp->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSpp->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSpp->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSpp->SetCheckLikeSign(kTRUE,kFALSE);
  pair_diffGamma_LSpp->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSpp->SetMotherPDGs(22,22);

  // Diff gamma same grandmother
  AliDielectronSignalMC* pair_diffGamma_sameGrandM = new AliDielectronSignalMC("MCpair_diffGamma_sameGrandM","MCpair_diffGamma_sameGrandM");
  pair_diffGamma_sameGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_sameGrandM->SetLegPDGs(11,-11);
  pair_diffGamma_sameGrandM->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffGamma_sameGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
  //mother
  pair_diffGamma_sameGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_sameGrandM->SetMotherPDGs(22,22);
  //Grandmother
  pair_diffGamma_sameGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kSame);

  AliDielectronSignalMC* pair_diffGamma_LSmm_sameGrandM = new AliDielectronSignalMC("MCpair_diffGamma_LSmm_sameGrandM","MCpair_diffGamma_LSmm_sameGrandM");
  pair_diffGamma_LSmm_sameGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSmm_sameGrandM->SetLegPDGs(11,11);
  pair_diffGamma_LSmm_sameGrandM->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSmm_sameGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSmm_sameGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSmm_sameGrandM->SetCheckLikeSign(kFALSE,kTRUE);
  pair_diffGamma_LSmm_sameGrandM->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSmm_sameGrandM->SetMotherPDGs(22,22);
  pair_diffGamma_LSmm_sameGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kSame);


  AliDielectronSignalMC* pair_diffGamma_LSpp_sameGrandM = new AliDielectronSignalMC("MCpair_diffGamma_LSpp_sameGrandM","MCpair_diffGamma_LSpp_sameGrandM");
  pair_diffGamma_LSpp_sameGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSpp_sameGrandM->SetLegPDGs(-11,-11);
  pair_diffGamma_LSpp_sameGrandM->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSpp_sameGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSpp_sameGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSpp_sameGrandM->SetCheckLikeSign(kTRUE,kFALSE);
  pair_diffGamma_LSpp_sameGrandM->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSpp_sameGrandM->SetMotherPDGs(22,22);
  pair_diffGamma_LSpp_sameGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kSame);

  // Diff gamma diff grandmother
  AliDielectronSignalMC* pair_diffGamma_diffGrandM = new AliDielectronSignalMC("MCpair_diffGamma_diffGrandM","MCpair_diffGamma_diffGrandM");
  pair_diffGamma_diffGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_diffGrandM->SetLegPDGs(11,-11);
  pair_diffGamma_diffGrandM->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pair_diffGamma_diffGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // for e from gamma: kFinalState is empty.
  //mother
  pair_diffGamma_diffGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_diffGrandM->SetMotherPDGs(22,22);
  pair_diffGamma_diffGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);

  AliDielectronSignalMC* pair_diffGamma_LSmm_diffGrandM = new AliDielectronSignalMC("MCpair_diffGamma_LSmm_diffGrandM","MCpair_diffGamma_LSmm_diffGrandM");
  pair_diffGamma_LSmm_diffGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSmm_diffGrandM->SetLegPDGs(11,11);
  pair_diffGamma_LSmm_diffGrandM->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSmm_diffGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSmm_diffGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSmm_diffGrandM->SetCheckLikeSign(kFALSE,kTRUE);
  pair_diffGamma_LSmm_diffGrandM->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSmm_diffGrandM->SetMotherPDGs(22,22);
  pair_diffGamma_LSmm_diffGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);

  AliDielectronSignalMC* pair_diffGamma_LSpp_diffGrandM = new AliDielectronSignalMC("MCpair_diffGamma_LSpp_diffGrandM","MCpair_diffGamma_LSpp_diffGrandM");
  pair_diffGamma_LSpp_diffGrandM->SetFillPureMCStep(fFillPureMC);
  pair_diffGamma_LSpp_diffGrandM->SetLegPDGs(-11,-11);
  pair_diffGamma_LSpp_diffGrandM->SetCheckBothChargesLegs(kFALSE,kFALSE);
  pair_diffGamma_LSpp_diffGrandM->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  pair_diffGamma_LSpp_diffGrandM->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  pair_diffGamma_LSpp_diffGrandM->SetCheckLikeSign(kTRUE,kFALSE);
  pair_diffGamma_LSpp_diffGrandM->SetCheckUnlikeSign(kFALSE); // not needed here but may save a bit of time.
  pair_diffGamma_LSpp_diffGrandM->SetMotherPDGs(22,22);
  pair_diffGamma_LSpp_diffGrandM->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);



  // ##################### "real" pairs from photon conversions in the detector material ##############################
  // AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
  // signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
  // signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
  // signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
  // signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
  // signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // die->AddSignalMC(signalFromResonance_ULS_gammaConv);
  //
  //
  // // #################### D-Mesons
  // AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
  // diEleOpenCharmCharged->SetLegPDGs(11,-11);
  // diEleOpenCharmCharged->SetMotherPDGs(401,401);
  // diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
  // diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
  // diEleOpenCharmCharged->SetPDGforStack(503);
  // diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(diEleOpenCharmCharged);
  //
  // AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
  // diEleOpenCharmNeutral->SetLegPDGs(11,-11);
  // diEleOpenCharmNeutral->SetMotherPDGs(405,405);
  // diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
  // diEleOpenCharmNeutral->SetPDGforStack(503);
  // diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(diEleOpenCharmNeutral);
  //
  // //B meson (3)
  // AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
  // diEleOneOpenB->SetLegPDGs(11,-11);
  // diEleOneOpenB->SetMotherPDGs(401,501);
  // diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOneOpenB->SetGrandMotherPDGs(501,0);
  // diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
  // diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  // die->AddSignalMC(diEleOneOpenB);
  //
  // //B meson (1)(1)
  // AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
  // diEleOpenB->SetLegPDGs(11,-11);
  // diEleOpenB->SetMotherPDGs(501,501);
  // diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(diEleOpenB);
  //
  // //B meson (2)(2)
  // AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
  // diEleOpenBtoD->SetLegPDGs(11,-11);
  // diEleOpenBtoD->SetMotherPDGs(401,401);
  // diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOpenBtoD->SetGrandMotherPDGs(501,501);
  // diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  // die->AddSignalMC(diEleOpenBtoD);
  //
  // //B meson (1)(2)
  // AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
  // diEleOpenBandBtoD->SetLegPDGs        (11,11);
  // diEleOpenBandBtoD->SetMotherPDGs     (401,501);
  // diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
  // diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  // //do i need this?
  // diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
  // die->AddSignalMC(diEleOpenBandBtoD);


    AliDielectronSignalMC* pion = new AliDielectronSignalMC("pion","Pion");
    pion->SetLegPDGs(11,-11);
    pion->SetMotherPDGs(111,-111);
    pion->SetMothersRelation(AliDielectronSignalMC::kSame);
    // pion->SetFillPureMCStep(kTRUE);
    pion->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    pion->SetCheckBothChargesLegs(kTRUE,kTRUE);
    pion->SetCheckBothChargesMothers(kTRUE,kTRUE);


    AliDielectronSignalMC* eta = new AliDielectronSignalMC("eta", "Eta");
    eta->SetLegPDGs(11,-11);
    eta->SetMotherPDGs(221,-221);
    eta->SetMothersRelation(AliDielectronSignalMC::kSame);
    // eta->SetFillPureMCStep(kTRUE);
    eta->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    eta->SetCheckBothChargesLegs(kTRUE,kTRUE);
    eta->SetCheckBothChargesMothers(kTRUE,kTRUE);


    AliDielectronSignalMC* omega = new AliDielectronSignalMC("omega","Omega");
    omega->SetLegPDGs(11,-11);
    omega->SetMotherPDGs(223,-223);
    omega->SetMothersRelation(AliDielectronSignalMC::kSame);
    // omega->SetFillPureMCStep(kTRUE);
    omega->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    omega->SetCheckBothChargesLegs(kTRUE,kTRUE);
    omega->SetCheckBothChargesMothers(kTRUE,kTRUE);


    AliDielectronSignalMC* phi = new AliDielectronSignalMC("phi", "Phi");
    phi->SetLegPDGs(11,-11);
    phi->SetMotherPDGs(333,-333);
    phi->SetMothersRelation(AliDielectronSignalMC::kSame);
    // phi->SetFillPureMCStep(kTRUE);
    phi->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    phi->SetCheckBothChargesLegs(kTRUE,kTRUE);
    phi->SetCheckBothChargesMothers(kTRUE,kTRUE);


    AliDielectronSignalMC* rho = new AliDielectronSignalMC("rho", "Rho");
    rho->SetLegPDGs(11,-11);
    rho->SetMotherPDGs(113,-113);
    rho->SetMothersRelation(AliDielectronSignalMC::kSame);
    // rho->SetFillPureMCStep(kTRUE);
    rho->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    rho->SetCheckBothChargesLegs(kTRUE,kTRUE);
    rho->SetCheckBothChargesMothers(kTRUE,kTRUE);


    AliDielectronSignalMC* etaprime = new AliDielectronSignalMC("etaprime", "Eta_Prime");
    etaprime->SetLegPDGs(11,-11);
    etaprime->SetMotherPDGs(331,-331);
    etaprime->SetMothersRelation(AliDielectronSignalMC::kSame);
    // etaprime->SetFillPureMCStep(kTRUE);
    etaprime->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    etaprime->SetCheckBothChargesLegs(kTRUE,kTRUE);
    etaprime->SetCheckBothChargesMothers(kTRUE,kTRUE);


    // decay channels
    // (1) D -> e X
    // (1) B -> e X
    // (2) B -> D X -> e X Y
    // (3) B -> e D X -> ee X Y always produces ULS pair

    AliDielectronSignalMC* promptJpsiNonRad= new AliDielectronSignalMC("promptJpsiNonRad","Prompt J/psi non-Radiative");   // prompt J/psi (not from beauty decays)
    promptJpsiNonRad->SetLegPDGs(11,-11);
    promptJpsiNonRad->SetMotherPDGs(443,443);
    promptJpsiNonRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
    promptJpsiNonRad->SetMothersRelation(AliDielectronSignalMC::kSame);
    // promptJpsiNonRad->SetFillPureMCStep(kTRUE);
    promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
    promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
    promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
    promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
    //

    AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
    promptJpsiRad->SetLegPDGs(11,-11);
    promptJpsiRad->SetMotherPDGs(443,443);
    promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
    promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
    // promptJpsiRad->SetFillPureMCStep(kTRUE);
    promptJpsiRad->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
    promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
    promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
    promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);





    // die->AddSignalMC(pion);
    // die->AddSignalMC(eta);
    // die->AddSignalMC(omega);
    // die->AddSignalMC(phi);
    // die->AddSignalMC(rho);
    // die->AddSignalMC(etaprime);
    // die->AddSignalMC(promptJpsiNonRad);
    // die->AddSignalMC(promptJpsiRad);

  // die->AddSignalMC(single_ele_dontcare);
  // die->AddSignalMC(single_ele_primary);
  // die->AddSignalMC(single_ele_secondary);
  // die->AddSignalMC(single_ele_finalstate);
  // die->AddSignalMC(single_ele_finalstate_fromBG);

  die->AddSignalMC(pair_sameMother);
  die->AddSignalMC(pair_sameGamma);
  // die->AddSignalMC(pair_diffM_oneGamma);
  // die->AddSignalMC(pair_diffM_oneGamma_LSmm);
  // die->AddSignalMC(pair_diffM_oneGamma_LSpp);
  // die->AddSignalMC(pair_diffGamma);
  // die->AddSignalMC(pair_diffGamma_LSmm);
  // die->AddSignalMC(pair_diffGamma_LSpp);
  //
  // die->AddSignalMC(pair_diffGamma_sameGrandM);
  // die->AddSignalMC(pair_diffGamma_LSmm_sameGrandM);
  // die->AddSignalMC(pair_diffGamma_LSpp_sameGrandM);
  // die->AddSignalMC(pair_diffGamma_diffGrandM);
  // die->AddSignalMC(pair_diffGamma_LSmm_diffGrandM);
  // die->AddSignalMC(pair_diffGamma_LSpp_diffGrandM);

}
