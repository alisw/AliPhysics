///
/// \file PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/ConfigFemtoAnalysisRun2.C
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

#include <TROOT.h>

#endif

typedef AliFemtoAnalysisPionPion AFAPP;

bool DEFAULT_DO_KT = kFALSE;

struct MacroParams {
  std::vector<int> centrality_ranges;
  std::vector<unsigned char> pair_codes;
  std::vector<float> kt_ranges;

  float qinv_bin_size_MeV;
  float qinv_max_GeV;

  UInt_t delta_phi_bin_count;
  Float_t delta_phi_min;
  Float_t delta_phi_max;

  UInt_t delta_eta_bin_count;
  Float_t delta_eta_min;
  Float_t delta_eta_max;

  UInt_t q3d_bin_count;
  Float_t q3d_maxq;

  UInt_t ylm_max_l;
  UInt_t ylm_ibin;
  Float_t ylm_vmin;
  Float_t ylm_vmax;

  bool ylm_useLCMS;

  bool do_qinv_cf;
  bool do_q3d_cf;
  bool do_deltaeta_deltaphi_cf;
  bool do_avg_sep_cf;
  bool do_trueq_cf;
  bool do_trueq3d_cf;

  bool do_kt_trueq3d_cf;

  bool do_kt_q3d;
  bool do_kt_qinv;
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
ConfigFemtoAnalysis(const TString& param_str="")
{
  std::cout << "[ConfigFemtoAnalysisRun2 (PionPion)]\n";

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
  macro_config.do_kt_trueq3d_cf = DEFAULT_DO_KT;

  macro_config.do_ylm_cf = false;
  macro_config.do_trueq_cf = true;
  macro_config.do_trueq3d_cf = true;

  macro_config.qinv_bin_size_MeV = 5.0f;
  macro_config.qinv_max_GeV = 1.0f;

  macro_config.delta_phi_bin_count = 72;
  macro_config.delta_phi_min = -0.1;
  macro_config.delta_phi_max =  0.1;

  macro_config.delta_eta_bin_count = 72;
  macro_config.delta_eta_min = -0.1;
  macro_config.delta_eta_max =  0.1;

  macro_config.q3d_bin_count = 120;
  macro_config.q3d_maxq = 0.30;

  macro_config.ylm_max_l = 3;
  macro_config.ylm_ibin = 30;
  macro_config.ylm_vmin = 0.0;
  macro_config.ylm_vmax = 0.3;
  macro_config.ylm_useLCMS = false;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

  // Update Configurations
  macro_config.do_kt_trueq3d_cf &= analysis_config.is_mc_analysis;

  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODMultSelection();
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

  AliFemtoModelManager *model_manager = NULL;

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

      analysis->AddStanardCutMonitors();

      if (macro_config.do_avg_sep_cf) {
        AliFemtoAvgSepCorrFctn *avgsep_cf = new AliFemtoAvgSepCorrFctn("avg_sep", 100, 0.0, 40.0);
        avgsep_cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
        analysis->AddCorrFctn(avgsep_cf);
      }

      if (macro_config.do_deltaeta_deltaphi_cf && !analysis_config.is_mc_analysis) {
        AliFemtoCorrFctnDPhiStarDEta *deta_dphi_cf = new AliFemtoCorrFctnDPhiStarDEta("_",
          cut_config.pair_phi_star_radius,
          macro_config.delta_phi_bin_count, macro_config.delta_phi_min, macro_config.delta_phi_max,
      	  macro_config.delta_eta_bin_count, macro_config.delta_eta_min, macro_config.delta_eta_max
        );
        deta_dphi_cf->SetMagneticFieldSign(1);
        analysis->AddCorrFctn(deta_dphi_cf);
      }

      const float QINV_MIN_VAL = 0.0,
                  QINV_MAX_VAL = macro_config.qinv_max_GeV;

      const int QINV_BIN_COUNT = TMath::Abs((QINV_MAX_VAL - QINV_MIN_VAL) * 1000 / macro_config.qinv_bin_size_MeV);

      if (macro_config.do_qinv_cf) {
        // TString cf_title = TString("_qinv_") + pair_type_str;
        TString cf_title("_qinv");
        int bin_count = (int)TMath::Abs(macro_config.qinv_max_GeV * 1000 / macro_config.qinv_bin_size_MeV));
        AliFemtoCorrFctn *cf = new AliFemtoQinvCorrFctn(cf_title.Data(), bin_count, 0.0, macro_config.qinv_max_GeV);
        analysis->AddCorrFctn(cf);
      }

      if (analysis_config.is_mc_analysis) {
          model_manager = new AliFemtoModelManager();
          AliFemtoModelWeightGeneratorBasic *weight_gen = new AliFemtoModelWeightGeneratorBasic();
          weight_gen->ShouldPrintEmptyParticleNotification(kFALSE);
          model_manager->AcceptWeightGenerator(weight_gen);
      }

      if (analysis_config.is_mc_analysis && macro_config.do_deltaeta_deltaphi_cf) {
        AliFemtoModelCorrFctnDEtaDPhiStar *mc_deta_dphistar_cf =
          AliFemtoModelCorrFctnDEtaDPhiStar::Build()
            .GroupOutput(true)
            .Radius(cut_config.pair_phi_star_radius)
            .Phi(macro_config.delta_phi_bin_count, macro_config.delta_phi_min, macro_config.delta_phi_max)
            .Eta(macro_config.delta_eta_bin_count, macro_config.delta_eta_min, macro_config.delta_eta_max)
            .Build();
        mc_deta_dphistar_cf->ConnectToManager(model_manager);

        analysis->AddCorrFctn(mc_deta_dphistar_cf);
      }

      if (macro_config.do_trueq_cf) {
        AliFemtoModelCorrFctn *model_cf = new AliFemtoModelCorrFctnTrueQ("_MC_CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        model_cf->ConnectToManager(model_manager);
        analysis->AddCorrFctn(model_cf);
      }

      if (macro_config.do_q3d_cf) {
        TString cf_title("_q3D");
        analysis->AddCorrFctn(new AliFemtoCorrFctn3DLCMSSym(cf_title, macro_config.q3d_bin_count, macro_config.q3d_maxq));
      }

      if (macro_config.do_trueq3d_cf) {
        // cf_title = TString("Trueq3D_") + pair_type_str;
        TString cf_title("Trueq3D_");
        AliFemtoModelCorrFctnTrueQ3D *m_cf = new AliFemtoModelCorrFctnTrueQ3D(cf_title, macro_config.q3d_bin_count, macro_config.q3d_maxq);
        m_cf->SetManager(model_manager);
        analysis->AddCorrFctn(m_cf);
      }

      if (macro_config.do_kt_qinv) {
        AliFemtoQinvCorrFctn *qinv_cf = new AliFemtoQinvCorrFctn(cf_title.Data(), bin_count, 0.0, macro_config.qinv_max_GeV);
        AliFemtoKtBinnedCorrFunc *kt_qinv = new AliFemtoKtBinnedCorrFunc("KT_Qinv", qinv_cf);
        // loop over (low,high) pairs in the kt_ranges vector
        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_qinv->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_qinv);
      }

      if (macro_config.do_kt_q3d && !macro_config.do_kt_trueq3d_cf) {
        // TString q3d_cf_name = TString("_q3D_") + pair_type_str;
        TString q3d_cf_name("_q3d");
        AliFemtoKtBinnedCorrFunc *kt_q3d = new AliFemtoKtBinnedCorrFunc("KT_Q3D",
          new AliFemtoCorrFctn3DLCMSSym(q3d_cf_name, macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_q3d->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_q3d);
      }

      if (macro_config.do_kt_trueq3d_cf) {
        TString q3d_cf_name("Trueq3D");

        AliFemtoModelCorrFctnTrueQ3D *m_cf = new AliFemtoModelCorrFctnTrueQ3D(q3d_cf_name, macro_config.q3d_bin_count, macro_config.q3d_maxq);
        m_cf->SetManager(model_manager);

        AliFemtoKtBinnedCorrFunc *kt_true_q3d = new AliFemtoKtBinnedCorrFunc("KT_TrueQ3D", m_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_true_q3d->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_true_q3d);
      }

      if (macro_config.do_ylm_cf) {
        TString name = TString::Format("ylm", chrgs[ichg], imult);
        AliFemtoCorrFctnDirectYlm *ylm_cf = new AliFemtoCorrFctnDirectYlm(
            name,
            macro_config.ylm_max_l,
            macro_config.ylm_ibin, // nbinssh,
            macro_config.ylm_vmin,
            macro_config.ylm_vmax,
            macro_config.ylm_useLCMS
          );
        analysis->AddCorrFctn(ylm_cf);
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
  std::cout << "I-BuildConfiguration (BASE64-Encoded): " << TBase64::Encode(text) << " \n";

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

    case '(':
    {
      UInt_t rangeend = text.Index(")");
      if (rangeend == -1) {
        std::cerr << "W-ConfigFemtoAnalysis: " << "Expected closing parens ')' in configuration string. Using rest of line as kT-bins\n";
        rangeend = line.Length();
      }
      TString kt_ranges = line(1, rangeend - 2);
      std::cout << "Loaded kt-ranges: '" << kt_ranges << "'\n";
      TObjArray *range_groups = kt_ranges.Tokenize(",");

      TIter next_range_group(range_groups);
      TObjString *range_group = NULL;

      while (range_group = (TObjString*)next_range_group()) {
        TObjArray *subrange = range_group->String().Tokenize(":");
        TIter next_subrange(subrange);
        TObjString *subrange_it = (TObjString *)next_subrange();
        TString prev = TString::Format("%0.6e", subrange_it->String().Atof());
        while (subrange_it = (TObjString *)next_subrange()) {
          TString next = TString::Format("%0.6e", subrange_it->String().Atof());

          cmd = macro_varname + ".kt_ranges.push_back(" + prev + ");";
          // std::cout << "`" << cmd << "`\n";
          gROOT->ProcessLineFast(cmd);

          cmd = macro_varname + ".kt_ranges.push_back(" + next + ");";
          // std::cout << "`" << cmd << "`\n";
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
