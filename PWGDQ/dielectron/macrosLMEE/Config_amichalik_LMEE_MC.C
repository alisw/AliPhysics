//#include "AliDielectron.h"023023

void      InitHistograms(AliDielectron *die, Int_t cutDefinition);
void      InitCF(AliDielectron* die, Int_t cutDefinition);
void      SetupMCsignals(AliDielectron* die);
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx};


// "ITSTPCTOFif_trkSPDfirst_1_kSemi", "ITSTPCTOFif_trkSPDfirst_1_kSemi_pt300",
//                           "ITSTPCTOFif_trkSPDfirst_kINT7_pt400", "ITSTPCTOFif_trkSPDfirst_kINT7_pt400_woPID",
//                           "kPbPb2015_Pt400_TOFif_TPCele1_5", "kITSTPCTOFif_trkSPDfirst_kINT7_pt100_woPID",
//                           "kPbPb2015_pidV0_pt400", "kPbPb2015_Pt400_ITSSA",
//                           "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq", "kPbPb2015_Pt400_TPCele_symITS_tightTOFif"};

// TString names=("kPion;kOmega;kAll");
// TString names=("kPDGLepton");
// TString names=("kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif");
// TString names=("kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif");
TString names=("kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif;kPDGLepton");
// TString names=("kPDGLepton;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif;");
TObjArray *arrNames = names.Tokenize(";");
const Int_t nDie = arrNames->GetEntries();


Bool_t MCenabled = kTRUE;//needed for LMEEcutlib
Bool_t isQAtask = kFALSE;

Bool_t SetMixing = kFALSE;
Bool_t SetPairing = kTRUE;


AliDielectron* Config_amichalik_LMEEPbPb(Int_t cutDefinition, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE)
{
  //
  // Setup the instance of AliDielectron
  //
  AnalysisCut AnaCut;

  LMEECutLib*  LMcutlib = new LMEECutLib();

  // task name
  TString name=Form("%02d",cutDefinition);
  if (cutDefinition<arrNames->GetEntriesFast())  name=arrNames->At(cutDefinition)->GetName();

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));
  MCenabled=hasMC;

  cout << "cutDefinition = " << cutDefinition << endl;
  // Setup Analysis Selection

  if (cutDefinition==1) {
    AnaCut.SetPIDAna(LMEECutLib::kPDGLepton);
    // AnaCut.SetAdditionalTrackCuts(LMEECutLib::kAll);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);

    AnaCut.SetPreFilterType(LMEECutLib::kNoPreFilter);
    AnaCut.SetPIDPre(LMEECutLib::kStandardPre);
    AnaCut.SetTrackSelectionPre(LMEECutLib::kPrefilter_cut1);
    AnaCut.SetPairCutsPre(LMEECutLib::kNoPairCutsPre);

    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to90);
    AnaCut.SetMixing(LMEECutLib::kEventMixing_1);
    AnaCut.SetESDTrackSelection(LMEECutLib::kStandardESD);
  }

  if (cutDefinition==0) {
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif);
    // AnaCut.SetAdditionalTrackCuts(LMEECutLib::kAll);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetPairCutsAna(LMEECutLib::kNoPairCutsAna);

    AnaCut.SetPreFilterType(LMEECutLib::kNoPreFilter);
    AnaCut.SetPIDPre(LMEECutLib::kStandardPre);
    AnaCut.SetTrackSelectionPre(LMEECutLib::kPrefilter_cut1);
    AnaCut.SetPairCutsPre(LMEECutLib::kNoPairCutsPre);

    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to90);
    AnaCut.SetMixing(LMEECutLib::kEventMixing_1);
    AnaCut.SetESDTrackSelection(LMEECutLib::kStandardESD);
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
  // deactivate pairing to check track cuts or run with loose pid cuts:
  if(SetPairing != kTRUE){
    die->SetNoPairing();
  }
  // apply correct Pre-Filter Scheme, if necessary
  // if (AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterAllSigns)  {
  //   die->SetPreFilterAllSigns();
  // }
  // else if (AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterUnlikeOnly) {
  //   die->SetPreFilterUnlikeOnly();
  // }


  //
  // Now configure task
  //
  // add centrality selection to event cuts
  die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(AnaCut) );
  // switch off KF Particle
  die->SetUseKF(kFALSE);

  AliDielectronCutGroup* trackCuts= new AliDielectronCutGroup("trackCuts","trackCuts",AliDielectronCutGroup::kCompAND);
  // trackCuts->AddCut(LMcutlib->GetAdditionalTrackCuts(AnaCut));
  trackCuts->AddCut(LMcutlib->GetPIDCutsAna(AnaCut));


 // --------------------------------------------------
  // with Rejection Step (Prefilter)
  // --------------------------------------------------
  if (AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterAllSigns || AnaCut.GetPreFilterType() == LMEECutLib::kPreFilterUnlikeOnly)
  {
    if (isESD) {
      die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
      die->GetPairPreFilterLegs().AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) ); // this is redundant!?
    }
    // set initial track filter.
    // the function 'GetPIDCutsPre()' must also call 'GetTrackCutsPre()'!
    die->GetTrackFilter().AddCuts( trackCuts );

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
    if (isESD) {
      die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
    }
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
    die->GetTrackFilter().AddCuts( trackCuts);
    // die->GetTrackFilter().AddCuts( LMcutlib->GetAdditionalTrackCuts(AnaCut) );
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(AnaCut, kFALSE) );
  }
  // -------------------------------------------------

  AliDielectronTrackRotator* rot= 0x0;
  //To save time and as it is not 100% test, rotation switched off
  /*AliDielectronTrackRotator *rot= LMcutlib->GetTrackRotator(selectedPID);
   die->SetTrackRotator(rot);
  */
  if (SetMixing == kTRUE){
    std::cout << "MIXING ACTIVATED" << std::endl;
    AliDielectronMixingHandler* mix = LMcutlib->GetMixingHandler(AnaCut);
    die->SetMixingHandler(mix);
  }

  SetupMCsignals(die);

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  InitHistograms(die,cutDefinition);

  // the last definition uses no cuts and only the QA histograms should be filled!
  //  InitCF(die,cutDefinition);

  return die;
}


void SetupMCsignals(AliDielectron* die){

// ################## "real" final state pairs ############################
  AliDielectronSignalMC* signalFromResonance_ULS = new AliDielectronSignalMC("signalFromResonance_ULS", "signalFromResonance_ULS");
  signalFromResonance_ULS->SetLegPDGs(11,-11);
  signalFromResonance_ULS->SetMothersRelation(AliDielectronSignalMC::kSame);
  signalFromResonance_ULS->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  signalFromResonance_ULS->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_ULS);

// ################# all finalstate pairs #################################
  AliDielectronSignalMC* signalFromResonance_all_ULS = new AliDielectronSignalMC("signalFromResonance_all_ULS", "signalFromResonance_all_ULS");
  signalFromResonance_all_ULS->SetLegPDGs(11,-11);
  // signalFromResonance_diff_ULS->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_all_ULS->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  signalFromResonance_all_ULS->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_all_ULS);


// ############### Like Sign pairs ############################
  AliDielectronSignalMC* signalFromResonance_LS_electron = new AliDielectronSignalMC("signalFromResonance_LS_electron", "signalFromResonance_LS_electron");
  signalFromResonance_LS_electron->SetLegPDGs(11,11);
  signalFromResonance_LS_electron->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_LS_electron->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  // signalFromResonance_LS_electron->SetCheckBothChargesLegs(kTRUE,kTRUE);
  signalFromResonance_LS_electron->SetCheckLikeSign(kFALSE, kTRUE);
  die->AddSignalMC(signalFromResonance_LS_electron);
//###############Like Sign
  AliDielectronSignalMC* signalFromResonance_LS_positron = new AliDielectronSignalMC("signalFromResonance_LS_positron", "signalFromResonance_LS_positron");
  signalFromResonance_LS_positron->SetLegPDGs(-11,-11);
  signalFromResonance_LS_positron->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_LS_positron->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  // signalFromResonance_LS_positron->SetCheckBothChargesLegs(kTRUE,kTRUE);
  signalFromResonance_LS_positron->SetCheckLikeSign(kTRUE, kFALSE);
  die->AddSignalMC(signalFromResonance_LS_positron);





  //###################### with photon conversion ################################

  // ##################### "real" pairs from photon conversions in the detector material ##############################
  AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
  signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
  signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
  signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
  signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
  signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_ULS_gammaConv);
  // #################### all pairs from photon conversions (real and non-real) ############################
  AliDielectronSignalMC* signalFromResonance_diff_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_diff_ULS_gammaConv", "signalFromResonance_diff_ULS_gammaConv");
  signalFromResonance_diff_ULS_gammaConv->SetLegPDGs(11,-11);
  signalFromResonance_diff_ULS_gammaConv->SetMotherPDGs(22,22);
  // signalFromResonance_diff_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_diff_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  signalFromResonance_diff_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_diff_ULS_gammaConv);


// ######################  Like Sign from electron pairs ###############################

  AliDielectronSignalMC* signalFromResonance_diff_LS_positron_gammaConv = new AliDielectronSignalMC("signalFromResonance_diff_LS_positron_gammaConv", "signalFromResonance_diff_LS_positron_gammaConv");
  signalFromResonance_diff_LS_positron_gammaConv->SetLegPDGs(-11,-11);
  signalFromResonance_diff_LS_positron_gammaConv->SetMotherPDGs(22,22);
  // signalFromResonance_diff_LS_positron_gammaConv->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_diff_LS_positron_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  // signalFromResonance_diff_LS_positron_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  signalFromResonance_diff_LS_positron_gammaConv->SetCheckLikeSign(kTRUE, kFALSE);
  die->AddSignalMC(signalFromResonance_diff_LS_positron_gammaConv);

// ##################### Like Sign from positron pairs ####################################
  AliDielectronSignalMC* signalFromResonance_diff_LS_electron_gammaConv = new AliDielectronSignalMC("signalFromResonance_diff_LS_electron_gammaConv", "signalFromResonance_diff_LS_electron_gammaConv");
  signalFromResonance_diff_LS_electron_gammaConv->SetLegPDGs(11,11);
  signalFromResonance_diff_LS_electron_gammaConv->SetMotherPDGs(22,22);
  // signalFromResonance_diff_LS_electron_gammaConv->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  signalFromResonance_diff_LS_electron_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  // signalFromResonance_diff_LS_electron_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  signalFromResonance_diff_LS_electron_gammaConv->SetCheckLikeSign(kFALSE, kTRUE);
  die->AddSignalMC(signalFromResonance_diff_LS_electron_gammaConv);


// ##############################################################
// ################# Pion contamination ###############################
// ##############################################################
  AliDielectronSignalMC* ele_pion_dif = new AliDielectronSignalMC("ele_pion_dif", "Ele_Pion_dif");
  ele_pion_dif->SetLegPDGs(11,211);
  ele_pion_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  ele_pion_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_pion_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(ele_pion_dif);

  AliDielectronSignalMC* ele_pion_same = new AliDielectronSignalMC("ele_pion_same", "Ele_Pion_same");
  ele_pion_same->SetLegPDGs(11,211);
  // ele_pion_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_pion_same->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_pion_same->SetFillPureMCStep(kTRUE);
  ele_pion_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_pion_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(ele_pion_same);

// ################################################
// ############# LS disturbed by pions ######################
// ####################################################
  AliDielectronSignalMC* ele_pion_LS_plus = new AliDielectronSignalMC("ele_pion_LS_plus", "Ele_Pion_same");
  ele_pion_LS_plus->SetLegPDGs(-11,211);
  // ele_pion_LS_plus->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_pion_LS_plus->SetFillPureMCStep(kTRUE);
  ele_pion_LS_plus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_pion_LS_plus->SetCheckLikeSign(kTRUE,kFALSE);
  die->AddSignalMC(ele_pion_LS_plus);

  AliDielectronSignalMC* ele_pion_same_LS_minus = new AliDielectronSignalMC("ele_pion_same_LS_minus", "Ele_Pion_same");
  ele_pion_same_LS_minus->SetLegPDGs(11,-211);
  // ele_pion_same_LS_minus->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_pion_same_LS_minus->SetFillPureMCStep(kTRUE);
  ele_pion_same_LS_minus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_pion_same_LS_minus->SetCheckLikeSign(kFALSE,kTRUE);
  die->AddSignalMC(ele_pion_same_LS_minus);
  // #################################



  // AliDielectronSignalMC* pion_pion_dif = new AliDielectronSignalMC("pion_pion_dif", "Pion_Pion_dif");
  // pion_pion_dif->SetLegPDGs(211,-211);
  // // pion_pion_dif->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // pion_pion_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);  //
  // // pion_pion_dif->SetFillPureMCStep(kTRUE);
  // pion_pion_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // pion_pion_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // die->AddSignalMC(pion_pion_dif);
  //
  //
  // AliDielectronSignalMC* pion_pion_same = new AliDielectronSignalMC("pion_pion_same", "Pion_Pion_same");
  // pion_pion_same->SetLegPDGs(211,-211);
  // // pion_pion_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // pion_pion_same->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // // pion_pion_same->SetFillPureMCStep(kTRUE);
  // pion_pion_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // pion_pion_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // die->AddSignalMC(pion_pion_same);


// ##########################################################################
// ############################ Proton Contamination ########################################
// ##########################################################################################
  AliDielectronSignalMC* ele_proton_dif = new AliDielectronSignalMC("ele_proton_dif", "Ele_Proton_dif");
  ele_proton_dif->SetLegPDGs(11,2212);
  // ele_proton_dif->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_proton_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);  //
  // ele_proton_dif->SetFillPureMCStep(kTRUE);
  ele_proton_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_proton_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // ele_proton_dif->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_proton_dif);

  AliDielectronSignalMC* ele_proton_same = new AliDielectronSignalMC("ele_proton_same", "Ele_Proton_same");
  ele_proton_same->SetLegPDGs(11,2212);
  // ele_proton_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_proton_same->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_proton_same->SetFillPureMCStep(kTRUE);
  ele_proton_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_proton_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  ele_proton_same->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_proton_same);

// #################################################################################
// ######################### LS proton #################################################
// ####################################################################
  AliDielectronSignalMC* ele_proton_LS_plus = new AliDielectronSignalMC("ele_proton_LS_plus", "Ele_Proton_same");
  ele_proton_LS_plus->SetLegPDGs(-11,2212);
  // ele_proton_LS_plus->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // ele_proton_LS_plus->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_proton_LS_plus->SetFillPureMCStep(kTRUE);
  ele_proton_LS_plus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_proton_LS_plus->SetCheckLikeSign(kTRUE,kFALSE);
  ele_proton_LS_plus->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_proton_LS_plus);

  AliDielectronSignalMC* ele_proton_LS_minus = new AliDielectronSignalMC("ele_proton_LS_minus", "Ele_Proton_same");
  ele_proton_LS_minus->SetLegPDGs(11,-2212);
  // ele_proton_LS_minus->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // ele_proton_LS_minus->SetMothersRelation(AliDielectronSignalMC::kSame);
  // ele_proton_LS_minus->SetFillPureMCStep(kTRUE);
  ele_proton_LS_minus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_proton_LS_minus->SetCheckLikeSign(kFALSE,kTRUE);
  ele_proton_LS_minus->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_proton_LS_minus);








  // AliDielectronSignalMC* proton_proton_dif = new AliDielectronSignalMC("proton_proton_dif", "Proton_Proton_dif");
  // proton_proton_dif->SetLegPDGs(2212,-2212);
  // // proton_proton_dif->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // proton_proton_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);  //
  // // proton_proton_dif->SetFillPureMCStep(kTRUE);
  // proton_proton_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // proton_proton_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // // proton_proton_dif->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(proton_proton_dif);
  //
  // AliDielectronSignalMC* proton_proton_same = new AliDielectronSignalMC("proton_proton_same", "Proton_Proton_same");
  // proton_proton_same->SetLegPDGs(2212,-2212);
  // // proton_proton_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // proton_proton_same->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // // proton_proton_same->SetFillPureMCStep(kTRUE);
  // proton_proton_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // proton_proton_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // proton_proton_same->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(proton_proton_same);

// ########################################################################################
// ################################ Kaon Contamination ######################################
// #########################################################################################
  AliDielectronSignalMC* ele_kaon_dif = new AliDielectronSignalMC("ele_kaon_dif", "Ele_Kaon_dif");
  ele_kaon_dif->SetLegPDGs(11,321);
  // ele_kaon_dif->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_kaon_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);  //
  // ele_kaon_dif->SetFillPureMCStep(kTRUE);
  ele_kaon_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_kaon_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // ele_kaon_dif->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_kaon_dif);

  AliDielectronSignalMC* ele_kaon_same = new AliDielectronSignalMC("ele_kaon_same", "Ele_Kaon_same");
  ele_kaon_same->SetLegPDGs(11,321);
  // ele_kaon_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_kaon_same->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // ele_kaon_same->SetFillPureMCStep(kTRUE);
  ele_kaon_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_kaon_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  ele_kaon_same->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_kaon_same);

// ###################################################################################
// ################################ LS Kaon ####################################
// ############################################################################


  AliDielectronSignalMC* ele_kaon_LS_plus = new AliDielectronSignalMC("ele_kaon_LS_plus", "Ele_Kaon_same");
  ele_kaon_LS_plus->SetLegPDGs(-11,321);
  // ele_kaon_LS_plus->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_kaon_LS_plus->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // ele_kaon_LS_plus->SetFillPureMCStep(kTRUE);
  ele_kaon_LS_plus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_kaon_LS_plus->SetCheckLikeSign(kTRUE,kFALSE);
  ele_kaon_LS_plus->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_kaon_LS_plus);

  AliDielectronSignalMC* ele_kaon_Ls_minus = new AliDielectronSignalMC("ele_kaon_Ls_minus", "Ele_Kaon_same");
  ele_kaon_Ls_minus->SetLegPDGs(11,-321);
  // ele_kaon_Ls_minus->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  ele_kaon_Ls_minus->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // ele_kaon_Ls_minus->SetFillPureMCStep(kTRUE);
  ele_kaon_Ls_minus->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  ele_kaon_Ls_minus->SetCheckLikeSign(kFALSE,kTRUE);
  ele_kaon_Ls_minus->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(ele_kaon_Ls_minus);


  // AliDielectronSignalMC* kaon_kaon_dif = new AliDielectronSignalMC("kaon_kaon_dif", "Kaon_Kaon_dif");
  // kaon_kaon_dif->SetLegPDGs(321,-321);
  // // kaon_kaon_dif->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // kaon_kaon_dif->SetMothersRelation(AliDielectronSignalMC::kDifferent);  //
  // // kaon_kaon_dif->SetFillPureMCStep(kTRUE);
  // kaon_kaon_dif->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // kaon_kaon_dif->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // // kaon_kaon_dif->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(kaon_kaon_dif);
  //
  // AliDielectronSignalMC* kaon_kaon_same = new AliDielectronSignalMC("kaon_kaon_same", "Kaon_Kaon_same");
  // kaon_kaon_same->SetLegPDGs(321,-321);
  // // kaon_kaon_same->SetMotherPDGs(pdg_code[any], -pdg_code[any]);
  // kaon_kaon_same->SetMothersRelation(AliDielectronSignalMC::kSame);  //
  // // kaon_kaon_same->SetFillPureMCStep(kTRUE);
  // kaon_kaon_same->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // kaon_kaon_same->SetCheckBothChargesLegs(kTRUE,kTRUE);
  // kaon_kaon_same->SetCheckBothChargesMothers(kTRUE,kTRUE);
  // die->AddSignalMC(kaon_kaon_same);







  AliDielectronSignalMC* pion = new AliDielectronSignalMC("pion","Pion");
  pion->SetLegPDGs(11,-11);
  pion->SetMotherPDGs(111,-111);
  pion->SetMothersRelation(AliDielectronSignalMC::kSame);
  // pion->SetFillPureMCStep(kTRUE);
  pion->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  pion->SetCheckBothChargesLegs(kTRUE,kTRUE);
  pion->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(pion);


  AliDielectronSignalMC* eta = new AliDielectronSignalMC("eta", "Eta");
  eta->SetLegPDGs(11,-11);
  eta->SetMotherPDGs(221,-221);
  eta->SetMothersRelation(AliDielectronSignalMC::kSame);
  // eta->SetFillPureMCStep(kTRUE);
  eta->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  eta->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eta->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(eta);


  AliDielectronSignalMC* omega = new AliDielectronSignalMC("omega","Omega");
  omega->SetLegPDGs(11,-11);
  omega->SetMotherPDGs(223,-223);
  omega->SetMothersRelation(AliDielectronSignalMC::kSame);
  // omega->SetFillPureMCStep(kTRUE);
  omega->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  omega->SetCheckBothChargesLegs(kTRUE,kTRUE);
  omega->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(omega);


  AliDielectronSignalMC* phi = new AliDielectronSignalMC("phi", "Phi");
  phi->SetLegPDGs(11,-11);
  phi->SetMotherPDGs(333,-333);
  phi->SetMothersRelation(AliDielectronSignalMC::kSame);
  // phi->SetFillPureMCStep(kTRUE);
  phi->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  phi->SetCheckBothChargesLegs(kTRUE,kTRUE);
  phi->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(phi);


  AliDielectronSignalMC* rho = new AliDielectronSignalMC("rho", "Rho");
  rho->SetLegPDGs(11,-11);
  rho->SetMotherPDGs(113,-113);
  rho->SetMothersRelation(AliDielectronSignalMC::kSame);
  // rho->SetFillPureMCStep(kTRUE);
  rho->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  rho->SetCheckBothChargesLegs(kTRUE,kTRUE);
  rho->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(rho);


  AliDielectronSignalMC* etaprime = new AliDielectronSignalMC("etaprime", "Eta_Prime");
  etaprime->SetLegPDGs(11,-11);
  etaprime->SetMotherPDGs(331,-331);
  etaprime->SetMothersRelation(AliDielectronSignalMC::kSame);
  // etaprime->SetFillPureMCStep(kTRUE);
  etaprime->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  etaprime->SetCheckBothChargesLegs(kTRUE,kTRUE);
  etaprime->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(etaprime);


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
  promptJpsiNonRad->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  promptJpsiNonRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiNonRad->SetJpsiRadiative(AliDielectronSignalMC::kIsNotRadiative);
  die->AddSignalMC(promptJpsiNonRad);
  //

  AliDielectronSignalMC* promptJpsiRad = new AliDielectronSignalMC("promptJpsiRad","Prompt J/psi Radiative");   // prompt J/psi (not from beauty decays)
  promptJpsiRad->SetLegPDGs(11,-11);
  promptJpsiRad->SetMotherPDGs(443,443);
  promptJpsiRad->SetGrandMotherPDGs(503,503,kTRUE,kTRUE);   // not from beauty hadrons
  promptJpsiRad->SetMothersRelation(AliDielectronSignalMC::kSame);
  // promptJpsiRad->SetFillPureMCStep(kTRUE);
  promptJpsiRad->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  promptJpsiRad->SetCheckBothChargesLegs(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesMothers(kTRUE,kTRUE);
  promptJpsiRad->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  promptJpsiRad->SetJpsiRadiative(AliDielectronSignalMC::kIsRadiative);
  die->AddSignalMC(promptJpsiRad);

  AliDielectronSignalMC* diEleOpenCharmCharged = new AliDielectronSignalMC("DmesonsCharged","di-electrons from open charm D+- mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmCharged->SetLegPDGs(11,-11);
  diEleOpenCharmCharged->SetMotherPDGs(401,401);
  diEleOpenCharmCharged->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmCharged->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  // diEleOpenCharmCharged->SetFillPureMCStep(kTRUE);
  diEleOpenCharmCharged->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmCharged->SetPDGforStack(503);
  diEleOpenCharmCharged->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmCharged->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmCharged);
// //D0 meson (1)(1)
// AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
//   diEleOpenCharmNeutral->SetLegPDGs(11,-11);
//   diEleOpenCharmNeutral->SetMotherPDGs(405,405);
//   diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//   diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
//   diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
//   diEleOpenCharmNeutral->SetPDGforStack(503);
//   diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
//   diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
//   die->AddSignalMC(diEleOpenCharmNeutral);
//B meson (3)
AliDielectronSignalMC* diEleOneOpenB = new AliDielectronSignalMC("B2ee","di-electrons from one B meson");  // dielectrons originating from open charm hadrons
  diEleOneOpenB->SetLegPDGs(11,-11);
  diEleOneOpenB->SetMotherPDGs(401,501);
  diEleOneOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOneOpenB->SetGrandMotherPDGs(501,0);
  diEleOneOpenB->SetCheckMotherGrandmotherRelation(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOneOpenB->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOneOpenB);

//B meson (1)(1)
AliDielectronSignalMC* diEleOpenB = new AliDielectronSignalMC("BMesons","di-electrons from B mesons");  // dielectrons originating from open charm hadrons
  diEleOpenB->SetLegPDGs(11,-11);
  diEleOpenB->SetMotherPDGs(501,501);
  diEleOpenB->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenB->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenB->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenB->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenB);

//B meson (2)(2)
AliDielectronSignalMC* diEleOpenBtoD = new AliDielectronSignalMC("B2D2ee","di-electrons from B->D-> e");  // dielectrons originating from open charm hadrons
  diEleOpenBtoD->SetLegPDGs(11,-11);
  diEleOpenBtoD->SetMotherPDGs(401,401);
  diEleOpenBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBtoD->SetGrandMotherPDGs(501,501);
  diEleOpenBtoD->SetGrandMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenBtoD);

//B meson (1)(2)
AliDielectronSignalMC* diEleOpenBandBtoD = new AliDielectronSignalMC("B2eAndB2D2e","di-electrons from B->e and B->D->e");  // dielectrons originating from open charm hadrons
  diEleOpenBandBtoD->SetLegPDGs        (11,11);
  diEleOpenBandBtoD->SetMotherPDGs     (401,501);
  diEleOpenBandBtoD->SetGrandMotherPDGs(501,0);
  diEleOpenBandBtoD->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenBandBtoD->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesMothers(kTRUE,kTRUE);
  diEleOpenBandBtoD->SetCheckBothChargesGrandMothers(kTRUE,kTRUE);
  //do i need this?
  diEleOpenBandBtoD->SetCheckMotherGrandmotherRelation(kTRUE,kFALSE);
  die->AddSignalMC(diEleOpenBandBtoD);

}


//______________________________________________________________________________________

void InitHistograms(AliDielectron *die, Int_t cutDefinition)
{
  //
  // Initialise the histograms
  //

  //Setup histogram Manager
  AliDielectronHistos *histos=new AliDielectronHistos(die->GetName(),die->GetTitle());

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
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
    // Legs of final Pairs. Both charges together. No duplicate entries.
    histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
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
    histos->AddClass(Form("Pre_%s",AliDielectron::TrackClassName(i)));
  }

  //Create Classes for Rejected Tracks/Pairs:
  for (Int_t i=0; i<3; ++i){
    histos->AddClass(Form("RejPair_%s",AliDielectron::PairClassName(i)));
    // Legs of rejected Pairs. Both charges together. One track can and will make multiple entries.
    histos->AddClass(Form("RejTrack_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistogramsPair(...)'
  }


  //add MC signal histograms to pair class
  if(die->GetMCSignals()) {
    for (Int_t i=0; i<die->GetMCSignals()->GetEntriesFast(); ++i) {
      //histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      histos->AddClass(Form("Pair_%s",die->GetMCSignals()->At(i)->GetName()));
      // histos->AddClass(Form("Track_%s",die->GetMCSignals()->At(i)->GetName()));
      histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      // if (fFillPureMC) histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
    }
  }
  /*
   //track rotation

   histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
   histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(AliDielectron::kEv1PMRot)));
   */

  //TVectorD* vecRunNumbers = AliDielectronHelper::MakeArbitraryBinning("167900, 167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593, 170600");

  //add histograms to event class
  histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  histos->UserHistogram("Event","RunNumber","",AliDielectronHelper::MakeArbitraryBinning("244917,244918"),AliDielectronVarManager::kRunNumber);
  histos->UserHistogram("Event","Centrality","","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","centrality","",100,0.,100,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","centrality_SPD","",100,0.,100,AliDielectronVarManager::kCentralityCL1);
  histos->UserHistogram("Event","nESDTracks","",8000,0,80000,AliDielectronVarManager::kNTrk);
  histos->UserHistogram("Event","Nacc","",8000,0,8000,AliDielectronVarManager::kNacc);
  histos->UserHistogram("Event","RefMultTPConly","",8000,0,8000,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","epTPC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kTPCrpH2uc);
  histos->UserHistogram("Event","epV0AC","",240,-TMath::Pi(),TMath::Pi(),AliDielectronVarManager::kv0ACrpH2);
  histos->UserHistogram("Event","epV0AC_epTPC","",120,-TMath::PiOver2(),TMath::PiOver2(),120,-TMath::PiOver2(),TMath::PiOver2(),AliDielectronVarManager::kTPCrpH2uc,AliDielectronVarManager::kv0ACrpH2);


  //add histograms to Track classes
  // axis labels are set to the values in 'AliDielectronVarManager.cxx' if the histogram title is empty or starts with ';'! [see AliDielectronHistos::AdaptNameTitle(...)]

  histos->UserHistogram("Track","pdg",";pdg;#tracks",2000,-1000., 1000.,AliDielectronVarManager::kPdgCode);
  histos->UserHistogram("Track","pdgMother",";pdg_{Mother};#tracks",2000,-1000., 1000.,AliDielectronVarManager::kPdgCodeMother);
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  histos->UserHistogram("Track","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
  histos->UserHistogram("Track","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
  histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);
  histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",
                        GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);

  // ITS
  histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                        GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
  histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
//histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
//                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
  // TPC
  histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                        GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                        BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                        GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                        GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);

  // if (isQAtask) {
  //   histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
  //                         GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
  //   histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p_{in} (GeV/c)",
  //                         GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
  //                         AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kPIn);
  //   histos->UserHistogram("Track","TPC_dEdx_P_RunNumber",";p_{in} (GeV/c);TPC signal (arb units);run",
  //                         GetVector(kP2D), GetVector(kTPCdEdx), GetVector(kRuns),
  //                         AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
  //   histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",
  //                         GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx),
  //                         AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTPCsignal);
  //
  //   histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
  //                         GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
  //                         AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
  //   histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
  //                         GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
  //                         AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
  //   histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
  //                         GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
  //                         AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
  //   histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
  //                         GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
  //                         AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
  //   histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
  //                         BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
  //                         AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
  //   histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_Nacc",";N_{TPC ref};n#sigma_{ele}^{TPC};N_{acc}",
  //                         BinsToVector(100,0.,5000.), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
  //                         AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
  //
  //   histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  //   histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
  //   histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
  //   histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
  //                         GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
  // }

  // TRD
  // TRD variables need lot of computing time. since a performance update by Julian, they will not be computed if not needed! (status 7.3.14)
  //  histos->UserHistogram("Track","TRDpidPobEle_P","TRD PID probability Electrons;p_{in} (GeV/c);TRD prob Electrons",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobEle);
  //  histos->UserHistogram("Track","TRDpidPobPio_P","TRD PID probability Pions;p_{in} (GeV/c);TRD prob Pions",
  //                        GetVector(kP2D), 100,0.,1.,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTRDprobPio);

  // TOF
  histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
                        GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  // if (isQAtask) {
  //   histos->UserHistogram("Track","TOFnSigmaPio_P",";p_{in} (GeV/c);TOF number of sigmas Pions",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPio);
  //   histos->UserHistogram("Track","TOFnSigmaKao_P",";p_{in} (GeV/c);TOF number of sigmas Kaons",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaKao);
  //   histos->UserHistogram("Track","TOFnSigmaPro_P",";p_{in} (GeV/c);TOF number of sigmas Protons",
  //                         GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaPro);
  // }
  // // 2D-PID
  // if (isQAtask) {
  //   histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
  //                         50,0.,2.5, 160,-12.,20., 150,-10.,20.,
  //                         AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
  //   histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
  //                         50,0.,2.5, 160,-12.,20., 50,-5.,5.,
  //                         AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);
  // }

  // Eta and Phi
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  // DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);

  // Quality
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
  histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);

  // if (isQAtask) {
  //   histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
  //                         160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
  //   histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
  //                         GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
  //   histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
  //                         GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
  //   histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
  //                         100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
  //   histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
  //                         120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
  // }

  if (!isQAtask)
  {
    //add histograms to Pair classes
    histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
    histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
    histos->UserHistogram("Pair","OpeningAngle","",240,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);

    //2D and 3D histograms
    histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                          GetVector(kMee), GetVector(kPtee),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
    histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                          100,-1.,1., 120,0.,TMath::TwoPi(),
                          AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
    histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                          GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                          GetVector(kMee), GetVector(kPtee), GetVector(kOpAng),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                          GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);

    //opening angle and PhiV
    histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                          GetVector(kMee), GetVector(kOpAng),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                          GetVector(kMee), GetVector(kPhiV),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                          GetVector(kPtee), GetVector(kOpAng),
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
    histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                          GetVector(kPtee), GetVector(kPhiV),
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
    histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                          GetVector(kOpAng), GetVector(kPhiV),
                          AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);

    //centrality
    histos->UserHistogram("Pair","InvMass_Centrality",";Inv. Mass [GeV];Centrality;#pairs",
                          GetVector(kMee), BinsToVector(102,-1,101),
                          AliDielectronVarManager::kM, AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Pair","PairPt_Centrality",";Pair Pt [GeV];Centrality;#pairs",
                          GetVector(kPtee), BinsToVector(102,-1,101),
                          AliDielectronVarManager::kPt, AliDielectronVarManager::kCentralityNew);



    //add histograms to Track classes
    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Pre","Px",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);

    histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                          GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                          GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);

    histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    //histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
    //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);

    histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);

  }// end: (!isQAtask)

  die->SetHistogramManager(histos);
}



TVectorD *GetVector(Int_t var)
{
  switch (var)
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);

    case kSigmaEle:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
      else          return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
    case kSigmaOther:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else          return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(140,  0.,140.);
      else          return AliDielectronHelper::MakeLinBinning(140,  0.,140.);

    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46,
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86,
                                                                   0.90, 0.94, 0.98, 1.02, 1.06,
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90,
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00
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
      //2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45,
      //2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95,
      //3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45,

    case kRuns:   return AliDielectronHelper::MakeArbitraryBinning("244917, 244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068, 246390, 246391, 246392");
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



void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);

  cf->AddVariable(AliDielectronVarManager::kCentralityNew,"0.,5.,10.,20.,30.,50.,80.,100.");

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,160,0.,8.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,350,0.,700.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,60,0.,120.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);

  //only in this case write MC truth info
  if (MCenabled) { // more elegant: die->GetHasMC()
    cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}
