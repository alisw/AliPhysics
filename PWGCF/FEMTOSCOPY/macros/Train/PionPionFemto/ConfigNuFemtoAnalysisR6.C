///
/// \file PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/ConfigNuFemtoAnalysisR6.C
///
/// \brief The configuration macro which sets up identical pion-pion analyses (ROOT6)
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

#if __cplusplus < 201103L
#error "This file requires use of the c++11 standard (ROOT6)"
#endif

#include "AliFemtoAnalysisPionPion.h"

#include "AliFemtoManager.h"
#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoEventReaderESDChainKine.h"
#include "AliFemtoEventReaderAODChain.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoEventReaderAlt.h"

#include "AliFemtoCutMonitorEventMult.h"
#include "AliFemtoCutMonitorParticleYPt.h"

#include "AliFemtoCorrFctnKStar.h"
#include "AliFemtoCorrFctnDEtaDPhiSimple.h"
#include "AliFemtoCorrFctnDEtaDPhiStar.h"

#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoShareQualityCorrFctn.h"
#include "AliFemtoCorrFctnDPhiStarDEta.h"
#include "AliFemtoQinvCorrFctn.h"
#include "AliFemtoCorrFctnQLCMS.h"
#include "AliFemtoModelGausRinvFreezeOutGenerator.h"
#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelWeightGeneratorBasic.h"
#include "AliFemtoModelWeightGeneratorLednicky.h"

#include "AliFemtoModelCorrFctnQinv.h"
#include "AliFemtoModelCorrFctnDEtaDPhiStar.h"
#include "AliFemtoModelCorrFctnDEtaDPhiAK.h"
#include "AliFemtoCorrFctn3DLCMSPosQuad.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoModelCorrFctnTrueQ6D.h"
#include "AliFemtoQinvCorrFctn.h"

#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoCorrFctnQ3D.h"

#include "AliFemtoModelCorrFctnTrueQ3D.h"
#include "AliFemtoModelCorrFctnTrueQ3DByParent.h"
#include "AliFemtoKtBinnedCorrFunc.h"
#include "AliFemtoModelCorrFctnTrueQ.h"
#include "AliFemtoCorrFctnInvMass.h"

#include <TROOT.h>
#include <TBase64.h>
#include <TNamed.h>

#include <iostream>
#include <vector>
#include <set>


using AFAPP = AliFemtoAnalysisPionPion;
using PionAnalysis = AliFemtoAnalysisPionPion;


struct MacroParams : public TNamed {
  MacroParams()
    : TNamed(AFAPP::make_random_string("macro_").Data(), "Macro Parameters")
    , centrality_ranges()
    , pair_codes()
    , kt_ranges()
    , obins()
    , slbins()
    {}

  MacroParams(const TString &name)
    : TNamed(name.Data(), "Macro Parameters")
    , centrality_ranges()
    , pair_codes()
    , kt_ranges()
    , obins()
    , slbins()
    {}

  MacroParams(PionAnalysis::AnalysisParams &a, PionAnalysis::CutParams &c, TString params);

  std::vector<std::pair<double,double>> centrality_ranges;
  std::set<PionAnalysis::PionType> pair_codes;
  std::vector<float> kt_ranges;

  std::vector<double> obins;
  std::vector<double> slbins;

  // monte carlo weight generator
  bool mcwg_lednicky { true };
  bool mcwg_square { false };
  bool mcwg_coul { true };
  bool mcwg_strong { true };
  bool mcwg_3body { true };

  // monte-carlo
  bool mcfg_disable { false };
  bool mcfg_use_1d { true };
  double mcfg_rinv { 5.0 };
  double mcfg_ro { 5.0 };
  double mcfg_rs { 5.0 };
  double mcfg_rl { 5.0 };

  bool eventreader_run1 { false };
  bool eventreader_use_alt { true };
  int eventreader_filter_bit { 7 };
  bool eventreader_multibit { false };
  int eventreader_read_full_mc { false };
  bool eventreader_epvzero { true };
  bool eventreader_vertex_shift { true };
  int eventreader_dca_globaltrack { 1 };
  bool eventreader_centrality_flattening { false };
  int eventreader_use_multiplicity { AliFemtoEventReaderAOD::kCentrality };

  float qinv_bin_size_MeV { 5.0 };
  float qinv_max_GeV { 1.0 };

  UInt_t delta_phi_bin_count { 72 };
  Float_t delta_phi_min { -0.1 };
  Float_t delta_phi_max { 0.1 };

  UInt_t delta_eta_bin_count { 72 };
  Float_t delta_eta_min { -0.1 };
  Float_t delta_eta_max { 0.1 };

  UInt_t q3d_bin_count { 41 };
  Float_t q3d_maxq { 0.1025 };

  UInt_t avgsep_nbins { 100 };
  Float_t avgsep_max { 20.0 };

  bool trueq3d_extra_bins { false };
  bool trueq3d_weighted_denoms { true };

  UInt_t ylm_max_l { 3 };
  UInt_t ylm_ibin { 30 };
  Float_t ylm_vmin { 0.0 };
  Float_t ylm_vmax { 0.3 };
  Bool_t ylm_useLCMS { true };

  bool do_mc_misident_analysis { false };

  bool do_deltaeta_deltaphi_cf { false };
  bool do_avg_sep_cf { false };
  bool do_kt_avg_sep_cf { false };

  bool do_q3d_varbins { false };

  bool cf_qlcms { false };
  bool do_qinv_cf { false };
  bool do_kt_qinv_cf { false };
  bool cf_kt_qlcms { false };
  bool do_q3d_cf { false };
  bool do_kt_q3d_cf { false };
  bool do_kt_lcms_cf { false };
  bool do_kt_prf_cf { false };
  bool do_kt_bf_cf { false };
  bool do_kt_q3dpf_cf { false };
  bool do_pqq3d_cf { false };
  bool do_kt_pqq3d_cf { false };

  bool do_moco6_cf { false };

  bool do_ylm_cf { false };
  bool do_kt_ylm_cf { false };

  bool do_minv_cf { false };
  bool do_kt_minv_cf { false };

  bool do_minvee_cf { false };
  bool do_kt_minvee_cf { false };

  // monte-carlo correlation functions
  bool do_trueq_cf { false };
  bool do_kt_trueq_cf { false };
  bool do_trueqinv_cf { false };
  bool do_kt_trueqinv_cf { false };
  bool do_trueq3d_cf { false };
  bool do_kt_trueq3d_cf { false };
  bool do_trueq3dbp_cf { false };
  bool do_kt_trueq3dbp_cf { false };

  bool do_sharequality_cf { false };

  bool do_truedetadphi_cf { false };
  bool do_kt_truedetadphi_cf { false };

  bool do_detadphi_simple_cf { false };
  bool do_kt_detadphi_simple_cf { false };

  bool do_detadphistar_cf { false };
  bool do_kt_detadphistar_cf { false };
  double phistar_radius { NAN };
  UInt_t detadphistar_nbins { 144 };

  bool RequestsMonteCarloData() const
    {
      return do_trueq_cf ||
             do_kt_trueq_cf ||
             do_trueqinv_cf ||
             do_kt_trueqinv_cf ||
             do_trueq3d_cf ||
             do_kt_trueq3d_cf ||
             do_trueq3dbp_cf ||
             do_kt_trueq3dbp_cf ||
             do_truedetadphi_cf ||
             do_kt_truedetadphi_cf;
    }

  // ClassDef(MacroParams, 1);
};

// AliFemtoCorrFctn* NewKtBinnedCorrFctn(TString name, AliFemtoCorrFctn *cf, std::vector<std::pair<float,float>> kt_ranges)
AliFemtoCorrFctn* NewKtBinnedCorrFctn(TString name, AliFemtoCorrFctn *cf, const std::vector<float> &kt_ranges)
{
  AliFemtoKtBinnedCorrFunc *ktcf = new AliFemtoKtBinnedCorrFunc(name, cf);
  // ktcf->AddKtRanges(kt_ranges);

  // for (auto kt_range : kt_ranges)
  for (size_t kt_idx=0; kt_idx < kt_ranges.size(); kt_idx += 2) {
    const float low = kt_ranges[kt_idx],
                high = kt_ranges[kt_idx+1];
    ktcf->AddKtRange(low, high);
  }

  return ktcf;
}


void AddKtBinnedCorrFctn(AliFemtoAnalysisPionPion &analysis,
                         TString name,
                         AliFemtoCorrFctn *cf,
                         const std::vector<float> &kt_ranges)
{
  auto *ktcf = NewKtBinnedCorrFctn(name, cf, kt_ranges);
  analysis.AddCorrFctn(ktcf);
}

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str="")
{
  std::cout << "[ConfigFemtoAnalysisRun2 Nu (PionPion)]\n";

  const double PionMass = 0.13956995;

  // default parameters for 3D hists with variable binsize (macro config: do_q3d_varbins)
  const std::vector<double>
    DEFAULT_OBINS = {0.        , 0.005     , 0.01028325, 0.01584975, 0.02169951,
                     0.02783251, 0.03424877, 0.04094828, 0.04793103, 0.05519704,
                     0.06274631, 0.07057882, 0.07869458, 0.08709360, 0.09577586,
                     0.10474138, 0.11399015, 0.12352217, 0.13333744, 0.14343596,
                     0.15381773, 0.16448276, 0.17543103, 0.18666256, 0.19817734,
                     0.20997537, 0.22205665, 0.23442118, 0.24706897, 0.26},
    DEFAULT_SBINS = {-0.26      , -0.24519231, -0.23077692, -0.21675385, -0.20312308,
                     -0.18988462, -0.17703846, -0.16458462, -0.15252308, -0.14085385,
                     -0.12957692, -0.11869231, -0.1082    , -0.0981    , -0.08839231,
                     -0.07907692, -0.07015385, -0.06162308, -0.05348462, -0.04573846,
                     -0.03838462, -0.03142308, -0.02485385, -0.01867692, -0.01289231,
                     -0.0075    , -0.0025    ,  0.0025    ,  0.0075    ,  0.01289231,
                      0.01867692,  0.02485385,  0.03142308,  0.03838462,  0.04573846,
                      0.05348462,  0.06162308,  0.07015385,  0.07907692,  0.08839231,
                      0.0981    ,  0.1082    ,  0.11869231,  0.12957692,  0.14085385,
                      0.15252308,  0.16458462,  0.17703846,  0.18988462,  0.20312308,
                      0.21675385,  0.23077692,  0.24519231,  0.26};

  // Get the default configurations
  AFAPP::AnalysisParams analysis_config = AFAPP::DefaultConfig();
  AFAPP::CutParams cut_config = AFAPP::DefaultCutConfig();
  MacroParams macro_config(analysis_config, cut_config, param_str);

  if (cut_config.pion_1_rm_neg_lbl && macro_config.eventreader_read_full_mc) {
    std::cerr << "\nWARNING : Requested reading full-MC AND removing tracks with negative labels\n\n";
  }

  // Begin to build the manager and analyses
  AliFemtoManager *manager = new AliFemtoManager();

  AliFemtoEventReaderAOD *rdr = (macro_config.eventreader_run1)
	                        ? new AliFemtoEventReaderAODChain()
	                      : (macro_config.eventreader_use_alt)
                                ? new AliFemtoEventReaderAlt()
                                : new AliFemtoEventReaderAODMultSelection();

    auto multest = static_cast<AliFemtoEventReaderAOD::EstEventMult>(macro_config.eventreader_use_multiplicity);
    const ULong_t filter_mask = macro_config.eventreader_multibit
                              ? macro_config.eventreader_filter_bit
                              : BIT(macro_config.eventreader_filter_bit);
    rdr->SetFilterMask(filter_mask);
    rdr->SetEPVZERO(macro_config.eventreader_epvzero);
    rdr->SetUseMultiplicity(multest);
    rdr->SetCentralityFlattening(macro_config.eventreader_centrality_flattening);
    rdr->SetReadV0(0);
    rdr->SetPrimaryVertexCorrectionTPCPoints(macro_config.eventreader_vertex_shift);
    rdr->SetDCAglobalTrack(macro_config.eventreader_dca_globaltrack);
    rdr->SetReadMC(analysis_config.is_mc_analysis);
    rdr->SetReadFullMCData(analysis_config.is_mc_analysis && macro_config.eventreader_read_full_mc);
    manager->SetEventReader(rdr);

  if (macro_config.centrality_ranges.empty()) {
    macro_config.centrality_ranges.emplace_back(0, 10.0);
  }

  const AFAPP::PionType
    PI_PLUS = AliFemtoAnalysisPionPion::kPiPlus,
    PI_MINUS = AliFemtoAnalysisPionPion::kPiMinus,
    NOPI = AliFemtoAnalysisPionPion::kNone;

  // default to identical pi+ and pi- analyses
  if (macro_config.pair_codes.empty()) {
    macro_config.pair_codes.insert(PI_PLUS);
    macro_config.pair_codes.insert(PI_MINUS);
  }

  AliFemtoModelManager *model_manager = nullptr;

  // loop over centrality ranges
  for (auto centrality_range : macro_config.centrality_ranges) {

    // loop over pair types
    for (AFAPP::PionType pair_type : macro_config.pair_codes) {

      const char *pair_type_str = pair_type == PI_PLUS ? "pip"
                                : pair_type == PI_MINUS ? "pim"
                                                        : "??";

      // build unique analysis name from centrality and pair types
      const TString analysis_name = Form("PiPiAnalysis_%02g_%02g_%s",
                                         centrality_range.first,
                                         centrality_range.second,
                                         pair_type_str);

      analysis_config.pion_type_1 = pair_type;
      analysis_config.pion_type_2 = NOPI;

      cut_config.event_centrality = centrality_range;

      auto *analysis = new AliFemtoAnalysisPionPion(analysis_name, analysis_config, cut_config);
      analysis->StoreEventReaderConfiguration(*rdr);

      analysis->AddStanardCutMonitors();

      const float QINV_MIN_VAL = 0.0,
                  QINV_MAX_VAL = macro_config.qinv_max_GeV;

      const int QINV_BIN_COUNT = TMath::Abs((QINV_MAX_VAL - QINV_MIN_VAL) * 1000 / macro_config.qinv_bin_size_MeV);

      if (macro_config.do_deltaeta_deltaphi_cf && !analysis_config.is_mc_analysis) {
        auto *deta_dphi_cf = new AliFemtoCorrFctnDPhiStarDEta("",
          cut_config.pair_phi_star_radius,
          macro_config.delta_phi_bin_count, macro_config.delta_phi_min, macro_config.delta_phi_max,
      	  macro_config.delta_eta_bin_count, macro_config.delta_eta_min, macro_config.delta_eta_max
        );
        deta_dphi_cf->SetMagneticFieldSign(1);
        analysis->AddCorrFctn(deta_dphi_cf);
      }

      if (macro_config.do_qinv_cf) {
        AliFemtoCorrFctn *qinv_cf = new AliFemtoQinvCorrFctn("_qinv", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        analysis->AddCorrFctn(qinv_cf);
      }

      if (analysis_config.is_mc_analysis) {
          model_manager = new AliFemtoModelManager();
          AliFemtoModelWeightGenerator *weight_gen = nullptr;

          if (macro_config.mcwg_lednicky) {
            AliFemtoModelWeightGeneratorLednicky *led_weight_gen = new AliFemtoModelWeightGeneratorLednicky();

            if (macro_config.mcwg_square) { led_weight_gen->SetSquare(); } else { led_weight_gen->SetSphere(); }
            if (macro_config.mcwg_coul) { led_weight_gen->SetCoulOn(); } else { led_weight_gen->SetCoulOff(); }
            if (macro_config.mcwg_strong) { led_weight_gen->SetStrongOn(); } else { led_weight_gen->SetStrongOff(); }
            if (macro_config.mcwg_3body) { led_weight_gen->Set3BodyOn(); } else { led_weight_gen->Set3BodyOff(); }

            weight_gen = led_weight_gen;
          }
          else {
            AliFemtoModelWeightGeneratorBasic *basic_weight_gen = new AliFemtoModelWeightGeneratorBasic();
            basic_weight_gen->ShouldPrintEmptyParticleNotification(kFALSE);
            weight_gen = basic_weight_gen;
          }

          AliFemtoModelFreezeOutGenerator *freeze_gen = nullptr;

          if (macro_config.mcfg_disable) {
          } else if (macro_config.mcfg_use_1d) {
            auto *ff = new AliFemtoModelGausRinvFreezeOutGenerator();
            ff->SetSizeInv(macro_config.mcfg_rinv);
            freeze_gen = ff;
          } else {
            auto *ff = new AliFemtoModelGausLCMSFreezeOutGenerator();
            ff->SetSizeOut(macro_config.mcfg_ro);
            ff->SetSizeSide(macro_config.mcfg_rs);
            ff->SetSizeLong(macro_config.mcfg_rl);
            freeze_gen = ff;
          }

          weight_gen->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());
          model_manager->AcceptWeightGenerator(weight_gen);
          model_manager->AcceptFreezeOutGenerator(freeze_gen);
      }

      if (analysis_config.is_mc_analysis && macro_config.do_deltaeta_deltaphi_cf) {
        auto *mc_deta_dphistar_cf =
          AliFemtoModelCorrFctnDEtaDPhiStar::Build()
            .GroupOutput(true)
            .Radius(cut_config.pair_phi_star_radius)
            .Phi(macro_config.delta_phi_bin_count, macro_config.delta_phi_min, macro_config.delta_phi_max)
            .Eta(macro_config.delta_eta_bin_count, macro_config.delta_eta_min, macro_config.delta_eta_max)
            .Build();
        mc_deta_dphistar_cf->ConnectToManager(model_manager);

        analysis->AddCorrFctn(mc_deta_dphistar_cf);
      }

      if (macro_config.do_avg_sep_cf) {
        auto *avgsep_cf = new AliFemtoAvgSepCorrFctn("AvgSep",
                                                     macro_config.avgsep_nbins,
                                                     0.0, macro_config.avgsep_max);
        avgsep_cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);
        analysis->AddCorrFctn(avgsep_cf);
      }

      if (macro_config.do_sharequality_cf) {
        auto *cf = new AliFemtoShareQualityCorrFctn("",
                                                    QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        analysis->AddCorrFctn(cf);
      }

      if (macro_config.do_kt_avg_sep_cf) {
        auto *avgsep_cf = new AliFemtoAvgSepCorrFctn("",
                                                     macro_config.avgsep_nbins,
                                                     0.0, macro_config.avgsep_max);
        avgsep_cf->SetPairType(AliFemtoAvgSepCorrFctn::kTracks);

        auto *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_AvgSep", avgsep_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }


      if (macro_config.do_trueq_cf) {
        AliFemtoModelCorrFctn *trueq_cf = new AliFemtoModelCorrFctnTrueQ("CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        trueq_cf->ConnectToManager(model_manager);
        analysis->AddCorrFctn(trueq_cf);
      }

      if (macro_config.do_kt_trueq_cf) {
        AliFemtoModelCorrFctn *cf = new AliFemtoModelCorrFctnTrueQ("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        cf->ConnectToManager(model_manager);
        AddKtBinnedCorrFctn(*analysis, "KT_Qlcms", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_trueqinv_cf) {
        AliFemtoModelCorrFctn *trueq_cf = new AliFemtoModelCorrFctnQinv("Qinv", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        trueq_cf->ConnectToManager(model_manager);
        analysis->AddCorrFctn(trueq_cf);
      }

      if (macro_config.do_kt_trueqinv_cf) {
        auto *qinv_cf = new AliFemtoModelCorrFctnQinv("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        qinv_cf->ConnectToManager(model_manager);

        auto *kt_cfs = NewKtBinnedCorrFctn("KT_TrueQinv", qinv_cf, macro_config.kt_ranges);
        analysis->AddCorrFctn(kt_cfs);
      }

      if (macro_config.do_trueq3d_cf) {
        AliFemtoModelCorrFctnTrueQ3D *trueq3d_cf = AliFemtoModelCorrFctnTrueQ3D::Build()
                                                          .NamePrefix("TrueQ3D_")
                                                          .BinCount(macro_config.q3d_bin_count)
                                                          .QRange(macro_config.q3d_maxq)
                                                          .Manager(model_manager);
        analysis->AddCorrFctn(trueq3d_cf);
      }

      if (macro_config.do_moco6_cf) {
        auto *moco6_cf = new AliFemtoModelCorrFctnTrueQ6D("MRC6D", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(moco6_cf);
      }

      if (macro_config.do_kt_qinv_cf) {
        AliFemtoQinvCorrFctn *qinv_cf = new AliFemtoQinvCorrFctn("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        auto *kt_qinv_cfs = NewKtBinnedCorrFctn("KT_Qinv", qinv_cf, macro_config.kt_ranges);
        analysis->AddCorrFctn(kt_qinv_cfs);
      }

      if (macro_config.cf_qlcms) {
        auto *cf = new AliFemtoCorrFctnQLCMS("", "Qlcms", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        analysis->AddCorrFctn(cf);
      }

      if (macro_config.cf_kt_qlcms) {
        auto *cf = new AliFemtoCorrFctnQLCMS("", "", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        AddKtBinnedCorrFctn(*analysis, "KT_Qlcms", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_q3d_cf) {
        auto *cf = new AliFemtoCorrFctn3DLCMSSym("Q3D", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(cf);
      }

      if (macro_config.do_kt_q3d_cf && !macro_config.do_kt_trueq3d_cf) {
        auto *cf = new AliFemtoCorrFctn3DLCMSSym("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        AddKtBinnedCorrFctn(*analysis, "KT_Q3D", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_kt_prf_cf && !macro_config.do_kt_trueq3d_cf) {
        auto *cf = new AliFemtoCorrFctnQ3DPF("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        AddKtBinnedCorrFctn(*analysis, "KT_Q3D_PRF", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_kt_bf_cf && !macro_config.do_kt_trueq3d_cf) {
        auto *cf = new AliFemtoCorrFctnQ3DBF("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        AddKtBinnedCorrFctn(*analysis, "KT_Q3D_BF", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_kt_lcms_cf && !macro_config.do_kt_trueq3d_cf) {
        auto *cf = new AliFemtoCorrFctnQ3DLCMS("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        AddKtBinnedCorrFctn(*analysis, "KT_Q3D_LCMS", cf, macro_config.kt_ranges);
      }

      if (macro_config.do_kt_q3dpf_cf && !macro_config.do_kt_trueq3d_cf) {
        auto *cf = new AliFemtoCorrFctnQ3DPF("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        AddKtBinnedCorrFctn(*analysis, "KT_Q3D_PF", cf, macro_config.kt_ranges);
      }

      const std::vector<double>
        obins = macro_config.obins.empty()
              ? DEFAULT_OBINS
              : macro_config.obins,

        slbins = macro_config.slbins.empty()
               ? DEFAULT_SBINS
               : macro_config.slbins;

      if (macro_config.do_pqq3d_cf) {
        auto *pq_cf = macro_config.do_q3d_varbins
                    ? new AliFemtoCorrFctn3DLCMSPosQuad("PQ3D", "", obins, slbins, slbins)
                    : new AliFemtoCorrFctn3DLCMSPosQuad("PQ3D", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(pq_cf);
      }

      if (macro_config.do_kt_pqq3d_cf) {
        auto *pq_cf = macro_config.do_q3d_varbins
                    ? new AliFemtoCorrFctn3DLCMSPosQuad("", "", obins, slbins, slbins)
                    : new AliFemtoCorrFctn3DLCMSPosQuad("", macro_config.q3d_bin_count, macro_config.q3d_maxq);

        AddKtBinnedCorrFctn(*analysis, "KT_PQ3D", pq_cf, macro_config.kt_ranges);
      }

      if (macro_config.do_kt_trueq3d_cf) {
        const std::vector<double>
          obins = macro_config.obins.empty()
                ? DEFAULT_OBINS
                : macro_config.obins,

          slbins = macro_config.slbins.empty()
                 ? DEFAULT_SBINS
                 : macro_config.slbins;

        auto *kt_trueq3d_cf = macro_config.do_q3d_varbins
                            ? new AliFemtoModelCorrFctnTrueQ3D("", obins, slbins, slbins, model_manager)
                            : AliFemtoModelCorrFctnTrueQ3D::Build()
                                  .NamePrefix("")
                                  .PosOut(macro_config.q3d_bin_count, macro_config.q3d_maxq)
                                  .EnableExtraHists(macro_config.trueq3d_extra_bins)
                                  .EnableWeightedDenominators(macro_config.trueq3d_weighted_denoms)
                                  .Manager(model_manager).into_ptr();

        AddKtBinnedCorrFctn(*analysis, "KT_TrueQ3D", kt_trueq3d_cf, macro_config.kt_ranges);
      }

      if (macro_config.do_trueq3dbp_cf) {
        auto *trueq3dbp_cf = new AliFemtoModelCorrFctnTrueQ3DByParent(macro_config.q3d_bin_count,
                                                                      macro_config.q3d_maxq);
        trueq3dbp_cf->SetManager(model_manager);
        analysis->AddCorrFctn(trueq3dbp_cf);
      }

      if (macro_config.do_kt_trueq3dbp_cf) {
        auto *kt_trueq3dbp_cf = new AliFemtoModelCorrFctnTrueQ3DByParent(macro_config.q3d_bin_count,
                                                                         macro_config.q3d_maxq);
        kt_trueq3dbp_cf->SetManager(model_manager);

        AddKtBinnedCorrFctn(*analysis, "KT_TrueQ3DByP", kt_trueq3dbp_cf, macro_config.kt_ranges);
      }

      if (macro_config.do_truedetadphi_cf) {
        auto *truedetadphi_cf = new AliFemtoModelCorrFctnDEtaDPhiAK("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(truedetadphi_cf);
      }

      if (macro_config.do_kt_truedetadphi_cf) {
        auto *kt_truedetadphi_cf = new AliFemtoModelCorrFctnDEtaDPhiAK("",
                                                                       macro_config.q3d_bin_count,
                                                                       macro_config.q3d_maxq);
        kt_truedetadphi_cf->ConnectToManager(model_manager);

        AddKtBinnedCorrFctn(*analysis, "KT_TrueDetaDphi", kt_truedetadphi_cf, macro_config.kt_ranges);
      }

      if (macro_config.do_detadphi_simple_cf) {
        auto *detadphi_simple_cf = new AliFemtoCorrFctnDEtaDPhiSimple("",
                                                                      macro_config.delta_phi_bin_count,
                                                                      macro_config.delta_eta_bin_count);
        detadphi_simple_cf->SetReadHiddenInfo(analysis_config.is_mc_analysis);
        analysis->AddCorrFctn(detadphi_simple_cf);
      }

      if (macro_config.do_kt_detadphi_simple_cf) {
        auto *kt_detadphi_simple_cf = new AliFemtoCorrFctnDEtaDPhiSimple("",
                                                                         macro_config.delta_phi_bin_count,
                                                                         macro_config.delta_eta_bin_count);
        kt_detadphi_simple_cf->SetReadHiddenInfo(analysis_config.is_mc_analysis);

        AddKtBinnedCorrFctn(*analysis, "KT_DetaDphiSimple", kt_detadphi_simple_cf, macro_config.kt_ranges);
      }

      if (macro_config.do_detadphistar_cf) {
        auto *starcf = AliFemtoCorrFctnDEtaDPhiStar::Params()
                        .Suffix("DEtaDPhiStar")
                        .Radius(macro_config.phistar_radius)
                        .NbinsEta(macro_config.detadphistar_nbins)
                        .NbinsPhi(macro_config.detadphistar_nbins)
                        .into_ptr();
        analysis->AddCorrFctn(starcf);
      }

      if (macro_config.do_kt_detadphistar_cf) {
        auto *kt_starcf = AliFemtoCorrFctnDEtaDPhiStar::Params()
                           .Suffix("")
                           .Radius(macro_config.phistar_radius)
                           .NbinsEta(macro_config.detadphistar_nbins)
                           .NbinsPhi(macro_config.detadphistar_nbins)
                           .into_ptr();

        AddKtBinnedCorrFctn(*analysis, "KT_DEtaDPhiStar", kt_starcf, macro_config.kt_ranges);
      }

      if (macro_config.do_minv_cf) {
        auto *minvcf = new AliFemtoCorrFctnInvMass("MINV", 3000, 0.0, 2.0, PionMass, PionMass);
        analysis->AddCorrFctn(minvcf);
      }

      if (macro_config.do_kt_minv_cf) {
        auto *minvcf = new AliFemtoCorrFctnInvMass("", 3000, 0.0, 2.0, PionMass, PionMass);
        AddKtBinnedCorrFctn(*analysis, "KT_MINV", minvcf, macro_config.kt_ranges);
      }

      const double E_MASS = 0.000511;

      if (macro_config.do_minvee_cf) {
        auto *minvcf = new AliFemtoCorrFctnInvMass("MINVEE", 3000, 0.0, 2.0, E_MASS, E_MASS);
        analysis->AddCorrFctn(minvcf);
      }

      if (macro_config.do_kt_minvee_cf) {
        auto *minvcf = new AliFemtoCorrFctnInvMass("", 3000, 0.0, 2.0, E_MASS, E_MASS);
        AddKtBinnedCorrFctn(*analysis, "KT_MINV_EE", minvcf, macro_config.kt_ranges);
      }

      if (macro_config.do_ylm_cf) {
        auto *ylm_cf = new AliFemtoCorrFctnDirectYlm("Ylm",
                                                     macro_config.ylm_max_l,
                                                     macro_config.ylm_ibin, // nbinssh,
                                                     macro_config.ylm_vmin,
                                                     macro_config.ylm_vmax,
                                                     macro_config.ylm_useLCMS);
        analysis->AddCorrFctn(ylm_cf);
      }

      if (macro_config.do_kt_ylm_cf) {
        auto *ylm_cf = new AliFemtoCorrFctnDirectYlm("",
                                                     macro_config.ylm_max_l,
                                                     macro_config.ylm_ibin, // nbinssh,
                                                     macro_config.ylm_vmin,
                                                     macro_config.ylm_vmax,
                                                     macro_config.ylm_useLCMS);

        AddKtBinnedCorrFctn(*analysis, "KT_DirectYlm", ylm_cf, macro_config.kt_ranges);
      }

      manager->AddAnalysis(analysis);

      if (macro_config.do_mc_misident_analysis) {
        if (!analysis_config.is_mc_analysis) {
          throw std::runtime_error("do_mc_misident_analysis requested for non-mc analysis");
        }
        if (!cut_config.cuts_use_attrs) {
          throw std::runtime_error("do_mc_misident_analysis requested for non 'attrs' cut");
        }

        TString misident_name = analysis_name + "_misident";
        AFAPP::AnalysisParams misident_config(analysis_config);
        AFAPP::CutParams misident_cut_config(cut_config);
        misident_cut_config.mc_nonpion_only = true;

        auto *misident_analysis = new AliFemtoAnalysisPionPion(misident_name, misident_config, misident_cut_config);
        misident_analysis->StoreEventReaderConfiguration(*rdr);
        misident_analysis->AddStanardCutMonitors();

        manager->AddAnalysis(misident_analysis);
      }
    }
  }

  return manager;
}


void
BuildConfiguration(const TString &text,
                   AliFemtoAnalysisPionPion::AnalysisParams &a,
                   AliFemtoAnalysisPionPion::CutParams &cut,
                   MacroParams &mac)
{
  std::cout << "I-BuildConfiguration (BASE64-Encoded): " << TBase64::Encode(text) << " \n";

  const TString analysis_varname = a.GetName(),
                    cut_varname = cut.GetName(),
  // const TString analysis_varname = Form("static_cast<AliFemtoAnalysisPionPion::AnalysisParams*>(%s)", a.GetName()),
                    //  cut_varname = Form("static_cast<AliFemtoAnalysisPionPion::CutParams*>(%s)", cut.GetName()),
                   macro_varname = Form("static_cast<MacroParams*>(%s)", mac.GetName());

  gDirectory->Add(&a);
  gDirectory->Add(&cut);
  gDirectory->Add(&mac);

  TObjArray* lines = text.Tokenize("\n;");

  TIter next_line(lines);
  TObject *line_obj = nullptr;

  while ((line_obj = next_line())) {

    const TString line = ((TObjString*)line_obj)->String().Strip(TString::kBoth, ' ');

    TString cmd("");

    switch (line[0]) {
    case '$':
      cmd = cut_varname + "->" + line(1, line.Length() - 1);
      break;

    case '@':
      cmd = analysis_varname + "->" + line(1, line.Length() - 1);
      break;

    case '~':
      cmd = macro_varname + "->" + line(1, line.Length() - 1);
      break;

    case '{':
    {
      UInt_t rangeend = line.Index("}");
      if (rangeend == static_cast<UInt_t>(-1)) {
        rangeend = line.Length();
      }
      TString centrality_ranges = line(1, rangeend - 1);
      std::cout << "Loading centrality-ranges: '" << centrality_ranges << "'\n";
      TObjArray *range_groups = centrality_ranges.Tokenize(",");

      TIter next_range_group(range_groups);
      TObjString *range_group = nullptr;

      while ((range_group = (TObjString*)next_range_group())) {
        TObjArray *subrange = range_group->String().Tokenize(":");
        TIter next_subrange(subrange);
        TObjString *subrange_it = (TObjString *)next_subrange();
        TString prev = TString::Format("%f", subrange_it->String().Atof());
        while ((subrange_it = (TObjString *)next_subrange())) {
          TString next = TString::Format("%f", subrange_it->String().Atof());

          cmd = macro_varname + "->centrality_ranges.emplace_back(" + prev + "," + next + ");";
          gROOT->ProcessLineFast(cmd);

          prev = next;
        }
      }
    }
    continue;

    case '(':
    {
      UInt_t rangeend = line.Index(")");
      if (rangeend == static_cast<UInt_t>(-1)) {
        std::cerr << "W-ConfigFemtoAnalysis: " << "Expected closing parens ')' in configuration string. Using rest of line as kT-bins\n";
        rangeend = line.Length();
      }
      TString kt_ranges = line(1, rangeend - 1);
      std::cout << "Loading kt-ranges: '" << kt_ranges << "'\n";
      TObjArray *range_groups = kt_ranges.Tokenize(",");

      TIter next_range_group(range_groups);
      TObjString *range_group = nullptr;

      while ((range_group = (TObjString*)next_range_group())) {
        TObjArray *subrange = range_group->String().Tokenize(":");
        TIter next_subrange(subrange);
        TObjString *subrange_it = (TObjString *)next_subrange();
        TString prev = TString::Format("%0.6e", subrange_it->String().Atof());
        while ((subrange_it = (TObjString *)next_subrange())) {
          TString next = TString::Format("%0.6e", subrange_it->String().Atof());

          cmd = macro_varname + "->kt_ranges.push_back(" + prev + ");";
          // std::cout << "`" << cmd << "`\n";
          gROOT->ProcessLineFast(cmd);

          cmd = macro_varname + "->kt_ranges.push_back(" + next + ");";
          // std::cout << "`" << cmd << "`\n";
          gROOT->ProcessLineFast(cmd);
          prev = next;
        }
      }
    }
    continue;

    case '+':
      if (line == "+p") {
        cmd = macro_varname + "->pair_codes.insert(AliFemtoAnalysisPionPion::kPiPlus)";
      }
      else if (line == "+m") {
        cmd = macro_varname + "->pair_codes.insert(AliFemtoAnalysisPionPion::kPiMinus)";
      }
      else {
        std::cerr << "\nError: Unknown command: " << line << "\n\n";
        throw std::runtime_error("Unknown command");
      }
      break;

    default:
      continue;
    }

    cmd += ";";

    std::cout << "I-BuildConfiguration: `" << cmd << "`\n";
    Int_t err = 0;
    gROOT->ProcessLineFast(cmd, &err);

    if (err != TInterpreter::EErrorCode::kNoError) {
      throw std::runtime_error(Form("Bad configuration line: `%s`", cmd.Data()));
    }
  }

  gDirectory->Remove(&a);
  gDirectory->Remove(&cut);
  gDirectory->Remove(&mac);

}

inline
MacroParams::MacroParams(PionAnalysis::AnalysisParams &analysis_config,
                         PionAnalysis::CutParams &c,
                         TString params)
  : MacroParams()
{
  // Read parameter string and update configurations
  BuildConfiguration(params, analysis_config, c, *this);


  if (std::isnan(phistar_radius)) {
    phistar_radius = c.pair_phi_star_radius;
  }

  if (RequestsMonteCarloData() && !analysis_config.is_mc_analysis) {
    std::cerr << "\nMonteCarlo correlation function requested when analysis is *not* MC\n\n";
    throw std::runtime_error("Macro misconfiguration : Not MC-CF requested in Non-MC analysis");
  }

}
