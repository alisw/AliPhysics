// #include "AliDielectron.h"
#include "THn.h"


void      InitHistograms(AliDielectron *die);
void      InitCF(AliDielectron* die);
TVectorD* BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD* GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx};

// TString names=("kPbPb2015_Pt400_tightTOFif_cent;kPbPb2015_Pt400_looseTOFif_cent;kPbPb2015_Pt400_tightTOFreq_cent;kPbPb2015_Pt400_tightTOFif_semi;kPbPb2015_Pt400_looseTOFif_semi;kPbPb2015_Pt400_tightTOFreq_semi");
// TString names=("kPbPb2015_Pt400_tightTOFif_MB;kPbPb2015_Pt400_looseTOFif_MB;kPbPb2015_Pt400_tightTOFreq_MB;kPbPb2015_Pt400_symTOFreq_MB");

// TString names=("kPbPb2015_pure_electron_pt400;kPbPb2015_pure_electron_pt400_woTOF;kPbPb2015_Pt400_symTOFreq_MB");
// TString names=("cut0;cut1;cut2;cut3;cut4;cut5;cut6;cut7;cut8;cut9;cut10;cut11;cut12;cut13;cut14;cut15;cut16;cut17;cut18;cut19");
TString names=("cut0;cut1;cut2;cut3;cut4;cut5;cut6;cut7;cut8;cut9");
// TString names=("cut10;cut11;cut12;cut13;cut14;cut15;cut16;cut17;cut18;cut19");

// TString names=("kPbPb2015_Pt400_tightTOFif_semi;kPbPb2015_Pt400_looseTOFif_semi;kPbPb2015_Pt400_tightTOFreq_semi");
// TString names=("kPbPb2015_Pt400_tightTOFif_semi");


// TString names=("kPbPb2015_noPID_NoTrackQualityCuts_Pt400"); // no track cuts
// TString names=("kPbPb2015_TOFif_NoTrackQualityCuts_Pt400");
// TString names=("kPionTOFeff_wTOF;kPionTOFeff_woTOF"); // For evaluating TOF matching efficiency

// Contamination Study
// TString names=("kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej;kPbPb2015_pidV0_electron_pt400_woPionRej;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut;kPbPb2015_pidV0_electron_pt400_wITScut;kPbPb2015_pure_pion_pt400_wITScut;kPbPb2015_pure_kaon_pt400_wITScut;kPbPb2015_pure_proton_pt400_wITScut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut;kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut;kPbPb2015_pidV0_electron_pt400_wTPCcut;kPbPb2015_pure_pion_pt400_wTPCcut;kPbPb2015_pure_kaon_pt400_wTPCcut;kPbPb2015_pure_proton_pt400_wTPCcut");


TObjArray *arrNames = names.Tokenize(";");
const Int_t nDie = arrNames->GetEntries();


Bool_t MCenabled = kFALSE;//needed for LMEEcutlib
Bool_t isQAtask = kFALSE;

Bool_t SetTPCCorrection = kTRUE;
Bool_t SetMixing = kTRUE;
Bool_t SetPairing = kTRUE;
Bool_t SetMCSignal = kFALSE;


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
AliDielectron* Config_caklein_LMEEPbPb_AOD(TString cutDefinition, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE)
{

  //
  // Setup the instance of AliDielectron
  //
  AnalysisCut AnaCut;

  LMEECutLib*  LMcutlib = new LMEECutLib();

  // task name
  TString name = cutDefinition;

  // init AliDielectron
  AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));
  if(SetTPCCorrection) LMcutlib->SetEtaCorrection(die, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, kFALSE);
  // die->SetCutQA();
  //die->SetHasMC(hasMC);
  MCenabled=hasMC;

  cout << "cutDefinition = " << cutDefinition << endl;
  // Setup Analysis Selection

  // ############ ANALYSIS CUTS
  if (cutDefinition == "kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFreq_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_symTOFreq_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_symTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFif_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_looseTOFif_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_looseTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut0"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_0);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_0);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }

// ########## CONTAMINATION STUDY
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetMixing(LMEECutLib::kEventMixing_1);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif_noPionRej);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFreq_noPionRej);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }


  // ############ CLEAN SAMPLES #############

  else if (cutDefinition == "kPbPb2015_pure_electron_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_electron_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kPbPb2015_pure_electron_pt400_woTOF"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_electron_pt400_woTOF);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }

  else if (cutDefinition == "kPbPb2015_pidV0_electron_pt400_wITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pidV0_electron_pt400_wITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pidV0_electron_pt400_woPionRej"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pidV0_electron_pt400_woPionRej);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();

  }
  else if (cutDefinition == "kPbPb2015_pure_pion_pt400_wITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_pion_pt400_wITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_kaon_pt400_wITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_kaon_pt400_wITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_proton_pt400_wITScut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_proton_pt400_wITScut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pidV0_electron_pt400_wTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pidV0_electron_pt400_wTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kV0);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_pion_pt400_wTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_pion_pt400_wTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_kaon_pt400_wTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_kaon_pt400_wTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_proton_pt400_wTPCcut"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_proton_pt400_wTPCcut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_pion_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_pion_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_kaon_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_kaon_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_pure_proton_pt400"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_pure_proton_pt400);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPbSemiCentral);
    AnaCut.SetStandardCut();
  }
  // TOF Efficiency Study
  else if (cutDefinition == "knanoAODTOFeffCut"){
    AnaCut.SetPIDAna(LMEECutLib::knanoAODTOFeffCut);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "knanoAODTOFeffCut_wTOF"){
    AnaCut.SetPIDAna(LMEECutLib::knanoAODTOFeffCut_wTOF);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPionTOFeff_wTOF"){
    AnaCut.SetPIDAna(LMEECutLib::kPionTOFeff_wTOF);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPionTOFeff_woTOF"){
    AnaCut.SetPIDAna(LMEECutLib::kPionTOFeff_woTOF);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }

  // ########## TRACK QUALITY STUDY

  else if (cutDefinition == "kPbPb2015_noPID_NoTrackQualityCuts_Pt400") {
    AnaCut.SetPIDAna(LMEECutLib::kITSTPCTOFif_trkSPDfirst_kINT7_pt400_woPID);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kNone);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_TOFif_NoTrackQualityCuts_Pt400") {
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_TPCele_AsymITS_tightTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kNone);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
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
    if (isESD) {
      die->GetTrackFilter().AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
    }
    // the function 'GetPIDCutsAna()' must also call 'GetTrackCutsAna()'!
    die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
    die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(AnaCut, kFALSE) );
  }
  // -------------------------------------------------

  // AliDielectronTrackRotator* rot= 0x0;
  //To save time and as it is not 100% test, rotation switched off
  /*AliDielectronTrackRotator *rot= LMcutlib->GetTrackRotator(selectedPID);
   die->SetTrackRotator(rot);
  */
  if (SetMixing == kTRUE){
    std::cout << "MIXING ACTIVATED" << std::endl;
    AliDielectronMixingHandler* mix = LMcutlib->GetMixingHandler(AnaCut);
    die->SetMixingHandler(mix);
  }

  if (SetMCSignal) SetupMCsignals(die);

  // histogram setup
  // only if an AliDielectronHistos object is attached to the
  // dielectron framework histograms will be filled
  //
  InitHistograms(die);

  // the last definition uses no cuts and only the QA histograms should be filled!
  //  InitCF(die);

  return die;
}

//______________________________________________________________________________________

void InitHistograms(AliDielectron *die)
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
      // histos->AddClass(Form("Track_%s",die->GetMCSignals()->At(i)->GetName()));
      histos->AddClass(Form("Track_Legs_%s",die->GetMCSignals()->At(i)->GetName()));
      // if (fFillPureMC) histos->AddClass(Form("Pair_%s_MCtruth",die->GetMCSignals()->At(i)->GetName()));
    }
  }

  // ################# EVENT ##########################
  // ##################################################
  histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
  // histos->UserHistogram("Event","Centrality","","-1,0,10,20,30,40,50,60,70,80,90,100,101;#events",AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","centrality","",100,0.,100,AliDielectronVarManager::kCentralityNew);
  histos->UserHistogram("Event","zVertex","",100,-15,15,AliDielectronVarManager::kZv);
  // histos->UserHistogram("Event","V0_vs_RefMult",";# Tracks V0;# RefMult;",AliDielectronHelper::MakeLinBinning(500,  0.,30000.),AliDielectronHelper::MakeLinBinning(500,  0.,15000.),AliDielectronVarManager::kMultV0,AliDielectronVarManager::kRefMult);
  // histos->UserHistogram("Event","V0_vs_NTrk",";# Tracks V0;# kNTrk;",AliDielectronHelper::MakeLinBinning(1000,  0.,30000.),AliDielectronHelper::MakeLinBinning(1000,  0.,15000.),AliDielectronVarManager::kMultV0,AliDielectronVarManager::kNTrk);
  // histos->UserHistogram("Event","V0_vs_TPCRefMult",";# Tracks V0;# RefMultTPConly;",AliDielectronHelper::MakeLinBinning(1000,  0.,30000.),AliDielectronHelper::MakeLinBinning(300,  0.,5000.),AliDielectronVarManager::kMultV0,AliDielectronVarManager::kRefMultTPConly);
  histos->UserHistogram("Event","nEvTPC_eventplaneents",";;ev plane;",AliDielectronHelper::MakeLinBinning(180,  TMath::Pi()/-2.,TMath::Pi()/2.),AliDielectronVarManager::kQnTPCrpH2);

  // ################# TRACK ##########################
  // ##################################################
  histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  // histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)", GetVector(kP2D), GetVector(kP2D), AliDielectronVarManager::kP,AliDielectronVarManager::kPIn);

  // ############## ITS
  // histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)", GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                        AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  -5,5.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  // histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
  //                       AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  0,10.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
  // histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
  //histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
  //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);

  // ############### TPC
  // histos->UserHistogram("Track","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
  //                       AliDielectronHelper::MakeLinBinning(500,  0.,5.), GetVector(kTPCdEdx), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);
  // histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  histos->UserHistogram("Track","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                        AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  -5,5.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
  //                       BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
  //                       BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
  //                       AliDielectronHelper::MakeLinBinning(500,  0.,5.), AliDielectronHelper::MakeLinBinning(100,  0,10.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  // histos->UserHistogram("Track","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
  // histos->UserHistogram("Track","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);
  histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                        GetVector(kEta3D), GetVector(kSigmaEle), AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);

  // histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_P",";p_{in} (GeV/c);Multiplicity_{TPC};n#sigma_{ele}^{TPC}",
  //                       AliDielectronHelper::MakeLinBinning(100,  0.,5.), AliDielectronHelper::MakeLinBinning(60,  0.,6000.), GetVector(kSigmaEle),
  //                       AliDielectronVarManager::kPIn, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTPCnSigmaEle);
  // //
  // histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";p_{in} (GeV/c);#eta;n#sigma_{ele}^{TPC}",
  //                       AliDielectronHelper::MakeLinBinning(100,  0.,5.), GetVector(kEta3D), GetVector(kSigmaEle),
  //                       AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle);
  //
  // histos->UserHistogram("Track","TPCnSigmaEle_Eta_TPCRefMult",";Eta;Multiplicity_{TPC};n#sigma_{ele}^{TPC}",
  //                       GetVector(kEta3D), AliDielectronHelper::MakeLinBinning(60,  0.,6000.), GetVector(kSigmaEle),
  //                       AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Track","TPCnSigmaEle_Phi",";Phi;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
  //                       GetVector(kEta3D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle);


  // ################## TOF
  // histos->UserHistogram("Track","TOFbeta_P",";p_{in} (GeV/c);TOF beta",
  //                       GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFbeta);
  histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);TOF number of sigmas Electrons",
                        GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTOFnSigmaEle);
  // ################## 2D-PID
  // histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
  //                       50,0.,2.5, 160,-12.,20., 150,-10.,20.,
  //                       AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
  // histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
  //                       50,0.,2.5, 160,-12.,20., 50,-5.,5.,
  //                       AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);


  // ############# Eta and Phi
  histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
  histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
  // histos->UserHistogram("Track","Phi_Pt","",200,0,10.,120,0.,TMath::TwoPi(),AliDielectronVarManager::kPt, AliDielectronVarManager::kPhi);
  // histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

  // ############## DCA
  histos->UserHistogram("Track","dXY","",200,-2.,2.,AliDielectronVarManager::kImpactParXY);
  histos->UserHistogram("Track","dZ" ,"",400,-4.,4.,AliDielectronVarManager::kImpactParZ);
  // histos->UserHistogram("Track","dXY_Pt","",200,0,10.,200,-2.,2.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY);
  // histos->UserHistogram("Track","dZ_Pt" ,"",200,0,10.,200,-2.,2.,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParZ);
  // histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);

  // ############### Quality
  histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",7,-0.5,6.5,AliDielectronVarManager::kNclsITS);
  histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kITSchi2Cl);
  histos->UserHistogram("Track","NclsSITS",";ITS shared clusters;#tracks",7,-0.5,6.5.,AliDielectronVarManager::kNclsSITS);
  histos->UserHistogram("Track","NclsSFracITS",";ITS Fraction of shared clusters;#tracks",100,0,1.,AliDielectronVarManager::kNclsSFracITS);
  histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNclsTPC);
  histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
  histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",100,0,1.,AliDielectronVarManager::kNclsSFracTPC);
  histos->UserHistogram("Track","TPCcrossedRows",";TPC crossed rows;#tracks",160,-0.5,159.5,AliDielectronVarManager::kNFclsTPCr);
  histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
  // histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);
  // histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",21,-0.5,20.5,AliDielectronVarManager::kTPCclsDiff);

    // histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
    //                       160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
    // histos->UserHistogram("Track","TPCcrossedRows_Pt",";Pt [GeV];TPC crossed rows",
    //                       GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCr);
    // histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Pt",";Pt [GeV];TPC crossed rows over findable",
    //                       GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kPt,AliDielectronVarManager::kNFclsTPCfCross);
    // histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
    //                       100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
    // histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
    //                       120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);


  // ################# PAIRS ##########################
  // ##################################################

  //2D and 3D histograms
  histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                        GetVector(kMee), GetVector(kPtee),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
  // histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
  //                       100,-1.,1., 120,0.,TMath::TwoPi(),
  //                       AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                        GetVector(kMee), GetVector(kPtee), GetVector(kPhiV),
                        AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
  // histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
  //                       GetVector(kMee), GetVector(kPtee), GetVector(kOpAng),
  //                       AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  // histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
  //                       GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2),
  //                       AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);

  // ################# opening angle and PhiV
  // histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
  //                       GetVector(kMee), GetVector(kOpAng), AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);


  // ################# PREFILTER ##########################
  // ##################################################
  // histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
  // histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
  //                       GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
  // histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
  //                       GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);

  // histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
  //                       GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
  // histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
  // histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
  //histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
  //                      GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);

  // histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
  // histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
  // histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
  // histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
  //                       GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);

  // const Int_t dimensions = 4;
  // Int_t bins[dimensions] = {100, 30, 40, 16};
  // UInt_t value[dimensions] = {AliDielectronVarManager::kPIn, AliDielectronVarManager::kRefMultTPConly, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kEta};
  // Double_t xmin[dimensions] = {0., 0., -4, -0.8};
  // Double_t xmax[dimensions] = {10., 3000., 4, 0.8};
  // histos->UserHistogram("Track", dimensions, bins, xmin, xmax, value);

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

    case kSigmaEle: return AliDielectronHelper::MakeLinBinning(50,-5.,5.);

    case kSigmaOther: return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);

    case kTPCdEdx: return AliDielectronHelper::MakeLinBinning(140,  0.,140.);

    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46,
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86,
                                                                   0.90, 0.94, 0.98, 1.02, 1.06,
                                                                   1.10, 1.20, 1.30, 1.40, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.80,
                                                                   2.90, 3.00, 3.05, 3.10, 3.20, 3.30, 3.50, 4.00, 4.50, 5.00
                                                                   ");
    // case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
    //                                                                0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
    //                                                                1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45,
    //                                                                1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95,
    //                                                                2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45,
    //                                                                2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95,
    //                                                                3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45,
    //                                                                3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95,
    //                                                                4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30, 4.35, 4.40, 4.45,
    //                                                                4.50, 4.55, 4.60, 4.65, 4.70, 4.75, 4.80, 4.85, 4.90, 4.95,
    //                                                                5.00
    //                                                                ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46,
                                                                   0.50
                                                                   ");
                                                                 );
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



void InitCF(AliDielectron* die)
{
  //
  // Setupd the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());

  //pair variables
  // cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  // cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  // cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  //
  // cf->AddVariable(AliDielectronVarManager::kCentralityNew,"0.,5.,10.,20.,30.,50.,80.,100.");

  //leg variables
  cf->AddVariable(AliDielectronVarManager::kPIn,80,0.,4.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kRefMultTPConly,"0.,100.,250.,500.,1000.,3000.",kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCnSigmaEle,GetVector(kSigmaEle),kTRUE);
  cf->AddVariable(AliDielectronVarManager::kEta,GetVector(kEta3D),kTRUE);

  //only in this case write MC truth info
  // if (MCenabled) { // more elegant: die->GetHasMC()
  //   cf->SetStepForMCtruth();
  //   cf->SetStepsForMCtruthOnly();
  //   cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
  //   cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  // }

  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}

void SetupMCsignals(AliDielectron* die){
  // ##################### "real" pairs from photon conversions in the detector material ##############################
  AliDielectronSignalMC* signalFromResonance_ULS_gammaConv = new AliDielectronSignalMC("signalFromResonance_ULS_gammaConv", "signalFromResonance_ULS_gammaConv");
  signalFromResonance_ULS_gammaConv->SetLegPDGs(11,-11);
  signalFromResonance_ULS_gammaConv->SetMotherPDGs(22,22);
  signalFromResonance_ULS_gammaConv->SetMothersRelation(AliDielectronSignalMC::kSame);
  signalFromResonance_ULS_gammaConv->SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary); // kSecondary means decays in the detector
  signalFromResonance_ULS_gammaConv->SetCheckBothChargesLegs(kTRUE,kTRUE);
  die->AddSignalMC(signalFromResonance_ULS_gammaConv);


  // #################### D-Mesons
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

  AliDielectronSignalMC* diEleOpenCharmNeutral = new AliDielectronSignalMC("DmesonsNeutral","di-electrons from open charm D0 mesons no B grandmother");  // dielectrons originating from open charm hadrons
  diEleOpenCharmNeutral->SetLegPDGs(11,-11);
  diEleOpenCharmNeutral->SetMotherPDGs(405,405);
  diEleOpenCharmNeutral->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  diEleOpenCharmNeutral->SetMothersRelation(AliDielectronSignalMC::kDifferent);
  diEleOpenCharmNeutral->SetCheckStackForPDG(kTRUE);
  diEleOpenCharmNeutral->SetPDGforStack(503);
  diEleOpenCharmNeutral->SetCheckBothChargesLegs(kTRUE,kTRUE);
  diEleOpenCharmNeutral->SetCheckBothChargesMothers(kTRUE,kTRUE);
  die->AddSignalMC(diEleOpenCharmNeutral);

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
