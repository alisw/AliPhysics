///
/// \file ConfigMomentumResolutionAnalysis.C
///
/// \brief The configuration macro which sets up identical pion-pion analyses
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if !defined(__CINT__) || defined(__MAKECINT_)

#include "AliFemtoAnalysisPionPion.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"

#include "AliFemtoCorrFctnKStar.h"
#include "AliFemtoModelCorrFctnTrueQ.h"
#include "AliFemtoMomentumResolver.h"

#include <TROOT.h>

#endif

typedef AliFemtoAnalysisPionPion AFAPP;

bool DEFAULT_DO_KT = kFALSE;

// #include "/home/akubera/development/physics/AliPhysics/inst/include/AliFemtoMomentumResolver.h"


struct MacroParams {
  std::vector<int> centrality_ranges;
  std::vector<unsigned char> pair_codes;
  float qinv_bin_size_MeV;
  float qinv_max_GeV;
  bool do_qinv_cf;
  bool do_q3d_cf;
  bool do_deltaeta_deltaphi_cf;
  bool do_avg_sep_cf;
  bool do_kt_q3d;
  bool do_kt_qinv;
  bool do_model_corrfctn;
  bool do_ylm_cf; // not implemented yet
};

void
BuildConfiguration(
  const TString&,
  AliFemtoAnalysisPionPion::AnalysisParams&,
  AliFemtoAnalysisPionPion::CutParams&,
  MacroParams&
);


AliFemtoManager*
// ConfigMomentumResolutionAnalysis(const TString& param_str="")
ConfigFemtoAnalysis(const TString& param_str="")
{
  std::cout << "[ConfigMomentumResolutionAnalysis (PionPion)]\n";

  // gROOT->ProcessLine(".include \"AliFemtoMomentumResolutionMonitor.h\"");


  const double PionMass = 0.13956995;

  // Get the default configurations
  AFAPP::AnalysisParams analysis_config = AFAPP::DefaultConfig();
  AFAPP::CutParams cut_config = AFAPP::DefaultCutConfig();
  MacroParams macro_config;

  // default
  macro_config.do_qinv_cf = true;
  macro_config.do_q3d_cf = true;
  macro_config.do_deltaeta_deltaphi_cf = false;
  macro_config.do_avg_sep_cf = false;
  macro_config.do_kt_q3d = macro_config.do_kt_qinv = DEFAULT_DO_KT;
  macro_config.do_ylm_cf = false;
  macro_config.qinv_bin_size_MeV = 5.0f;
  macro_config.qinv_max_GeV = 1.0f;
  macro_config.do_model_corrfctn = true;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODChain();
    rdr->SetFilterBit(7);
    rdr->SetEPVZERO(kTRUE);
    rdr->SetUseMultiplicity(AliFemtoEventReaderAOD::kCentrality);
    rdr->SetCentralityFlattening(kFALSE);
    rdr->SetReadV0(0);
    // rdr->SetPrimaryVertexCorrectionTPCPoints(kTRUE);
    rdr->SetDCAglobalTrack(kTRUE);
    rdr->SetReadMC(analysis_config.is_mc_analysis);
  manager->SetEventReader(rdr);

  if (macro_config.centrality_ranges.empty()) {
    // macro_config.centrality_ranges.push_back(std::pair<int, int>(0.0, 300));
    macro_config.centrality_ranges.push_back(0);
    macro_config.centrality_ranges.push_back(10.0);
  }

  // default to identical pi+ and pi- analyses
  if (macro_config.pair_codes.empty()) {
    macro_config.pair_codes.push_back(1);
    macro_config.pair_codes.push_back(0);
  }


  AFAPP::PionType PI_PLUS = AliFemtoAnalysisPionPion::kPiPlus,
                 PI_MINUS = AliFemtoAnalysisPionPion::kPiMinus,
                     NOPI = AliFemtoAnalysisPionPion::kNone;

  // loop over centrality ranges
  for (int cent_it = 0; cent_it + 1 < macro_config.centrality_ranges.size(); cent_it += 2) {

    const int mult_low = macro_config.centrality_ranges[cent_it],
             mult_high = macro_config.centrality_ranges[cent_it + 1];

    const TString mult_low_str = TString::Format("%0.2i", mult_low),
                 mult_high_str = TString::Format("%0.2i", mult_high);

    // loop over pair types
    for (int pair_it = 0; pair_it < macro_config.pair_codes.size(); ++pair_it) {

      const int pair_code = macro_config.pair_codes[pair_it];

      AFAPP::PionType type_1 = (pair_code & 0x1) ? PI_PLUS : PI_MINUS,
                      type_2 = NOPI;

      switch (pair_code & 0x6) {
      case 0:
        break;
      case 2:
        type_2 = PI_PLUS;
        break;
      case 4:
        type_2 = PI_MINUS;
        break;
      default:
        cout << "W-ConfigFemtoAnalysis: Invalid pair code " << pair_code << ". Skipping.\n";
        continue;
      }

      const TString pair_type_str =
        TString(type_1 == AliFemtoAnalysisPionPion::kPiPlus ? "pip" : "pim")
        + (type_2 == AliFemtoAnalysisPionPion::kNone ? "" : type_2 == AliFemtoAnalysisPionPion::kPiPlus ? "_pip" : "_pim");

      // build unique analysis name from centrality and pair types
      const TString analysis_name = TString::Format(
        "PiPiAnalysis_%0.2i_%0.2i_%s", mult_low, mult_high, pair_type_str.Data()
      );

      analysis_config.pion_type_1 = type_1;
      analysis_config.pion_type_2 = type_2;

      cut_config.event_CentralityMin = mult_low;
      cut_config.event_CentralityMax = mult_high;

      AliFemtoAnalysisPionPion *analysis = new AliFemtoAnalysisPionPion(analysis_name, analysis_config, cut_config);
      // AliFemtoCutMonitor *mc_resolver = new AliFemtoMomentumResolver()
      // analysis->PairCut()->AddCutMonitor(mc_resolver, new AliFemtoMomentumResolver());

      analysis->AddStanardCutMonitors();

      if (macro_config.do_avg_sep_cf) {
        AliFemtoAvgSepCorrFctn *avgsep_cf = new AliFemtoAvgSepCorrFctn("avg_sep", 100, 0.0, 40.0);
        avgsep_cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
        analysis->AddCorrFctn(avgsep_cf);
      }

      if (macro_config.do_deltaeta_deltaphi_cf) {
        AliFemtoCorrFctnDPhiStarDEta *deta_dphi_cf = new AliFemtoCorrFctnDPhiStarDEta("_", 1.6,
              // 100, 0.0, 1.6,
              // 100, 0.0, 2.0
              // 200, 0.0, 0.5,
              // 200, 0.0, 0.5
              75, -0.1, 0.1,
              75, -0.1, 0.1
        );
        deta_dphi_cf->SetMagneticFieldSign(1);

        analysis->AddCorrFctn(deta_dphi_cf);
      }

      const float QINV_MIN_VAL = 0.0,
                  QINV_MAX_VAL = macro_config.qinv_max_GeV;

      const int QINV_BIN_COUNT = TMath::Abs((QINV_MAX_VAL - QINV_MIN_VAL) * 1000 / macro_config.qinv_bin_size_MeV);

      if (macro_config.do_qinv_cf) {
        TString cf_title = TString("_qinv_") + pair_type_str;
        analysis->AddCorrFctn(
          new AliFemtoQinvCorrFctn(cf_title.Data(), QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL)
        );
      }

      if (macro_config.do_model_corrfctn) {
//         AliFemtoModelCorrFctnQinv *cf = new AliFemtoModelCorrFctnQinv("CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
//         cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);

          AliFemtoModelCorrFctn *model_cf = new AliFemtoModelCorrFctn("_MC_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
          AliFemtoModelManager *model_manager = new AliFemtoModelManager();

          AliFemtoModelWeightGeneratorBasic *weight_gen = new AliFemtoModelWeightGeneratorBasic();
          weight_gen->ShouldPrintEmptyParticleNotification(kFALSE);
          model_manager->AcceptWeightGenerator(weight_gen);

          model_cf->ConnectToManager(model_manager);

        analysis->AddCorrFctn(
          model_cf
        );
      }

      if (macro_config.do_q3d_cf) {
        cf_title = TString("_q3D_") + pair_type_str;
        analysis->AddCorrFctn(new AliFemtoCorrFctn3DLCMSSym(cf_title, 72, 0.25));
      }

      if (macro_config.do_kt_qinv) {
        // create base CF for all Kt-binned cf
        AliFemtoCorrFctn *kt_qinv_cf = macro_config.do_model_corrfctn
                                     //  ? new AliFemtoModelCorrFctnQinv("MC_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL)
                                     ? new AliFemtoModelCorrFctn("_MC_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL)
                                     : new AliFemtoQinvCorrFctn("_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);

        if (macro_config.do_model_corrfctn) {
          AliFemtoModelManager *model_manager = new AliFemtoModelManager();
          model_manager->AcceptWeightGenerator(new AliFemtoModelWeightGeneratorBasic());
          ((AliFemtoModelCorrFctn*)kt_qinv_cf)->ConnectToManager(model_manager);
        }

        AliFemtoKtBinnedCorrFunc *binned = new AliFemtoKtBinnedCorrFunc("KT_Qinv", kt_qinv_cf);
        binned->AddKtRange(0.2, 0.3);
        binned->AddKtRange(0.3, 0.4);
        binned->AddKtRange(0.4, 0.5);
        binned->AddKtRange(0.5, 0.6);
        binned->AddKtRange(0.6, 0.7);
        binned->AddKtRange(0.7, 0.8);
        binned->AddKtRange(0.8, 1.0);
        analysis->AddCorrFctn(binned);

      }

      if (macro_config.do_kt_q3d) {
        AliFemtoCorrFctn *kt_q3d_cf = (false && macro_config.do_model_corrfctn)
                                    ? NULL // new AliFemtoModelCorrFctnQinv("MC_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL)
                                    : new AliFemtoCorrFctn3DLCMSSym(TString("_q3D") + pair_type_str, QINV_BIN_COUNT, 0.5);

        AliFemtoKtBinnedCorrFunc *kt_q3d = new AliFemtoKtBinnedCorrFunc("_KT_Q3D", kt_q3d_cf);
        kt_q3d->AddKtRange(0.2, 0.3);
        kt_q3d->AddKtRange(0.3, 0.4);
        kt_q3d->AddKtRange(0.4, 0.5);
        kt_q3d->AddKtRange(0.5, 0.6);
        kt_q3d->AddKtRange(0.6, 0.7);
        kt_q3d->AddKtRange(0.7, 0.8);
        kt_q3d->AddKtRange(0.8, 1.0);
        analysis->AddCorrFctn(kt_q3d);
      }

      if (macro_config.do_ylm_cf) {
        /*
        analysis->AddCorrFctn(
          new AliFemtoCorrFctnDirectYlm(
            TString::Format("cylm%stpcM%i", chrgs[ichg], imult),
            3,
            nbinssh,
            0.0,
            shqmax,
            runshlcms
          )
        );
        */
      }

      manager->AddAnalysis(analysis);
    }
  }

  return manager;
}


void
BuildConfiguration(const TString &text,
                   AliFemtoAnalysisPionPion::AnalysisParams &a,
                   AliFemtoAnalysisPionPion::CutParams &cut,
                   MacroParams &mac
                   )
{
  std::cout << "I-BuildConfiguration:" << TBase64::Encode(text) << " \n";
  //   std::cout << "   '" << text << "'\n";

  const TString analysis_varname = "a",
                     cut_varname = "cut",
                   macro_varname = "mac";


  TObjArray* lines = text.Tokenize("\n;");

  TIter next_line(lines);
  TObject *line_obj = NULL;

  while (line_obj = next_line()) {

    const TString line = ((TObjString*)line_obj)->String().Strip(TString::kBoth, ' ');

    TString cmd("");

    switch (line[0]) {
    case '$':
      cmd = cut_varname + "." + line(1, line.Length() - 1);
      break;

    case '@':
      cmd = analysis_varname + "." + line(1, line.Length() - 1);
      break;

    case '~':
      cmd = macro_varname + "." + line(1, line.Length() - 1);
      break;

    case '{':
    {
      UInt_t rangeend = text.Index("}");
      if (rangeend == -1) {
        rangeend = line.Length();
      }
      TString centrality_ranges = line(1, rangeend - 1);
      TObjArray *range_groups = centrality_ranges.Tokenize(",");

      TIter next_range_group(range_groups);
      TObjString *range_group = NULL;

      while (range_group = (TObjString*)next_range_group()) {
        TObjArray *subrange = range_group->String().Tokenize(":");
        TIter next_subrange(subrange);
        TObjString *subrange_it = (TObjString *)next_subrange();
        TString prev = TString::Format("%0.2d", subrange_it->String().Atoi());
        while (subrange_it = (TObjString *)next_subrange()) {
          TString next = TString::Format("%0.2d", subrange_it->String().Atoi());

          cmd = macro_varname + ".centrality_ranges.push_back(" + prev + ");";
          gROOT->ProcessLineFast(cmd);

          cmd = macro_varname + ".centrality_ranges.push_back(" + next + ");";
          gROOT->ProcessLineFast(cmd);
          prev = next;
        }
      }
    }
    continue;

    case '+':
      if (line == "+p") {
        cmd = macro_varname + ".pair_codes.push_back(1)";
      }
      else if (line == "+m") {
        cmd = macro_varname + ".pair_codes.push_back(0)";
      }
      else if (line == "+pp") {
        cmd = macro_varname + ".pair_codes.push_back(3)";
      }
      else if (line == "+mm") {
        cmd = macro_varname + ".pair_codes.push_back(4)";
      }
      else if (line == "+pm") {
        cmd = macro_varname + ".pair_codes.push_back(2)";
      }
      else if (line == "+mp") {
        cmd = macro_varname + ".pair_codes.push_back(5)";
      }
      else {
        continue;
      }
      break;

    default:
      continue;
    }

    cmd += ";";

    cout << "I-BuildConfiguration: `" << cmd << "`\n";
    gROOT->ProcessLineFast(cmd);
  }
}
