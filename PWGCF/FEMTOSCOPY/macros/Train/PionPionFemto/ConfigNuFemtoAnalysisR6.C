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

#include <TROOT.h>
#include <TBase64.h>
#include <TNamed.h>


using AFAPP = AliFemtoAnalysisPionPion;


struct MacroParams : public TNamed {
  MacroParams()
    : TNamed(AFAPP::make_random_string("macro_").Data(), "Macro Parameters")
    {}

  MacroParams(const TString &name)
    : TNamed(name.Data(), "Macro Parameters")
    {}

  std::vector<int> centrality_ranges;
  std::vector<unsigned char> pair_codes;
  std::vector<float> kt_ranges;

  // monte carlo weight generator
  bool mcwg_lednicky { true };
  bool mcwg_square { false };
  bool mcwg_coul { true };
  bool mcwg_strong { true };
  bool mcwg_3body { true };

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

  UInt_t q3d_bin_count { 51 };
  Float_t q3d_maxq { 0.1275 };

  UInt_t avgsep_nbins { 100 };
  Float_t avgsep_max { 20.0 };

  bool trueq3d_extra_bins { false };
  bool trueq3d_weighted_denoms { true };

  UInt_t ylm_max_l { 3 };
  UInt_t ylm_ibin { 30 };
  Float_t ylm_vmin { 0.0 };
  Float_t ylm_vmax { 0.3 };
  Bool_t ylm_useLCMS { true };

  bool do_deltaeta_deltaphi_cf { false };
  bool do_avg_sep_cf { false };
  bool do_kt_avg_sep_cf { false };

  bool do_qinv_cf { false };
  bool do_kt_qinv_cf { false };
  bool do_q3d_cf {false};
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
  double phistar_radius { 1.2 };
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


// ClassImp(MacroParams);

void
BuildConfiguration(
  const TString&,
  AliFemtoAnalysisPionPion::AnalysisParams&,
  AliFemtoAnalysisPionPion::CutParams&,
  MacroParams&
);

/*
AliFemtoManager*
ConfigFemtoAnalysis(const AliFemtoConfigObject& cfg)
{
  AliFemtoConfigObject evreader_cfg = AliFemtoConfigObject::BuildMap()
    ("_class", "AliFemtoEventReaderAODMultSelection")
    ("filter_bit", 7)
    ("epvzero", true)
    ("dca_globaltrack", true)
    ("centrality_flattenint", false)
    ("use_multiplicity", stati_cast<int>(AliFemtoEventReaderAOD::kCentrality));

  AliFemtoConfigObject::Construct<EventReader>(evreader_cfg)
  AliFemtoEventReaderAOD *rdr = new AliFemtoEventReaderAODMultSelection();

}
*/

AliFemtoManager*
ConfigFemtoAnalysis(const TString& param_str="")
{
  std::cout << "[ConfigFemtoAnalysisRun2 Nu (PionPion)]\n";

  const double PionMass = 0.13956995;

  // Get the default configurations
  AFAPP::AnalysisParams analysis_config = AFAPP::DefaultConfig();
  AFAPP::CutParams cut_config = AFAPP::DefaultCutConfig();
  MacroParams macro_config;

  // Read parameter string and update configurations
  BuildConfiguration(param_str, analysis_config, cut_config, macro_config);

  // Update Configurations
  if (macro_config.RequestsMonteCarloData() && !analysis_config.is_mc_analysis) {
    std::cerr << "\nMonteCarlo correlation function requested when analysis is *not* MC\n\n";
    return nullptr;
  }

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
    // macro_config.centrality_ranges.push_back(std::pair<int, int>(0.0, 300));
    macro_config.centrality_ranges.push_back(0);
    macro_config.centrality_ranges.push_back(10.0);
  }

  // default to identical pi+ and pi- analyses
  if (macro_config.pair_codes.empty()) {
    macro_config.pair_codes.push_back(1);
    macro_config.pair_codes.push_back(0);
  }

  AliFemtoModelManager *model_manager = nullptr;

  AFAPP::PionType PI_PLUS = AliFemtoAnalysisPionPion::kPiPlus,
                 PI_MINUS = AliFemtoAnalysisPionPion::kPiMinus,
                     NOPI = AliFemtoAnalysisPionPion::kNone;

  // loop over centrality ranges
  for (int cent_it = 0; cent_it + 1 < macro_config.centrality_ranges.size(); cent_it += 2) {

    const int cent_low = macro_config.centrality_ranges[cent_it],
             cent_high = macro_config.centrality_ranges[cent_it + 1];

    const TString cent_low_str = TString::Format("%0.2i", cent_low),
                 cent_high_str = TString::Format("%0.2i", cent_high);

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
        throw std::runtime_error(Form("Invalid pair code %d", pair_code));
      }

      const TString pair_type_str =
        TString(type_1 == AliFemtoAnalysisPionPion::kPiPlus ? "pip" : "pim")
        + (type_2 == AliFemtoAnalysisPionPion::kNone ? "" : type_2 == AliFemtoAnalysisPionPion::kPiPlus ? "_pip" : "_pim");

      // build unique analysis name from centrality and pair types
      const TString analysis_name = TString::Format(
        "PiPiAnalysis_%0.2i_%0.2i_%s", cent_low, cent_high, pair_type_str.Data()
      );

      analysis_config.pion_type_1 = type_1;
      analysis_config.pion_type_2 = type_2;

      cut_config.event_CentralityMin = cent_low;
      cut_config.event_CentralityMax = cent_high;

      AliFemtoAnalysisPionPion *analysis = new AliFemtoAnalysisPionPion(analysis_name, analysis_config, cut_config);
      analysis->SetTrackFilter(filter_mask);

      analysis->AddStanardCutMonitors();

      const float QINV_MIN_VAL = 0.0,
                  QINV_MAX_VAL = macro_config.qinv_max_GeV;

      const int QINV_BIN_COUNT = TMath::Abs((QINV_MAX_VAL - QINV_MIN_VAL) * 1000 / macro_config.qinv_bin_size_MeV);


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

          weight_gen->SetPairType(AliFemtoModelWeightGenerator::PionPlusPionPlus());
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
        AliFemtoModelCorrFctn *trueq_cf = new AliFemtoModelCorrFctnTrueQ("CF", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        trueq_cf->ConnectToManager(model_manager);
        analysis->AddCorrFctn(trueq_cf);
      }

      if (macro_config.do_kt_trueq_cf) {
        AliFemtoModelCorrFctn *cf = new AliFemtoModelCorrFctnTrueQ("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        cf->ConnectToManager(model_manager);
        auto *kt_qinv_cfs = new AliFemtoKtBinnedCorrFunc("KT_TrueQ", cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_qinv_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_qinv_cfs);
      }

      if (macro_config.do_trueqinv_cf) {
        AliFemtoModelCorrFctn *trueq_cf = new AliFemtoModelCorrFctnQinv("Qinv", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        trueq_cf->ConnectToManager(model_manager);
        analysis->AddCorrFctn(trueq_cf);
      }

      if (macro_config.do_kt_trueqinv_cf) {
        auto *qinv_cf = new AliFemtoModelCorrFctnQinv("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        qinv_cf->ConnectToManager(model_manager);

        auto *kt_qinv_cfs = new AliFemtoKtBinnedCorrFunc("KT_TrueQinv", qinv_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_qinv_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_qinv_cfs);
      }

      if (macro_config.do_q3d_cf) {
        analysis->AddCorrFctn(new AliFemtoCorrFctn3DLCMSSym("Q3D", macro_config.q3d_bin_count, macro_config.q3d_maxq));
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
        AliFemtoModelCorrFctnTrueQ6D *moco6_cf = new AliFemtoModelCorrFctnTrueQ6D("MRC6D", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(moco6_cf);
      }

      if (macro_config.do_kt_qinv_cf) {
        AliFemtoQinvCorrFctn *kt_qinv_cf = new AliFemtoQinvCorrFctn("", QINV_BIN_COUNT, QINV_MIN_VAL, QINV_MAX_VAL);
        AliFemtoKtBinnedCorrFunc *kt_qinv_cfs = new AliFemtoKtBinnedCorrFunc("KT_Qinv", kt_qinv_cf);
        // loop over (low,high) pairs in the kt_ranges vector
        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_qinv_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_qinv_cfs);
      }

      if (macro_config.do_kt_q3d_cf && !macro_config.do_kt_trueq3d_cf) {
        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_Q3D",
          new AliFemtoCorrFctn3DLCMSSym("_q3d", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_kt_prf_cf && !macro_config.do_kt_trueq3d_cf) {
        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_Q3D_PRF",
          new AliFemtoCorrFctnQ3DPF("", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_kt_bf_cf && !macro_config.do_kt_trueq3d_cf) {
        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_Q3D_BF",
          new AliFemtoCorrFctnQ3DBF("", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_kt_lcms_cf && !macro_config.do_kt_trueq3d_cf) {

        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_Q3D_LCMS",
          new AliFemtoCorrFctnQ3DLCMS("", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_kt_q3dpf_cf && !macro_config.do_kt_trueq3d_cf) {
        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_Q3D_PF",
          new AliFemtoCorrFctnQ3DPF("", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_pqq3d_cf) {
        AliFemtoCorrFctn3DLCMSPosQuad *cf = new AliFemtoCorrFctn3DLCMSPosQuad("PQ3D", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(cf);
      }

      if (macro_config.do_kt_pqq3d_cf) {
        AliFemtoKtBinnedCorrFunc *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("KT_PQ3D",
          new AliFemtoCorrFctn3DLCMSPosQuad("", macro_config.q3d_bin_count, macro_config.q3d_maxq));

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_binned_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_binned_cfs);
      }

      if (macro_config.do_kt_trueq3d_cf) {
        AliFemtoModelCorrFctnTrueQ3D *kt_trueq3d_cf = AliFemtoModelCorrFctnTrueQ3D::Build()
                                                          .NamePrefix("")
                                                          .BinCount(macro_config.q3d_bin_count)
                                                          .QRange(macro_config.q3d_maxq)
                                                          .EnableExtraHists(macro_config.trueq3d_extra_bins)
                                                          .EnableWeightedDenominators(macro_config.trueq3d_weighted_denoms)
                                                          .Manager(model_manager);

        AliFemtoKtBinnedCorrFunc *kt_trueq3d_cfs = new AliFemtoKtBinnedCorrFunc("KT_TrueQ3D", kt_trueq3d_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_trueq3d_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_trueq3d_cfs);
      }

      if (macro_config.do_trueq3dbp_cf) {
        AliFemtoModelCorrFctnTrueQ3DByParent *trueq3dbp_cf = new AliFemtoModelCorrFctnTrueQ3DByParent(macro_config.q3d_bin_count,
                                                                                              macro_config.q3d_maxq);
        trueq3dbp_cf->SetManager(model_manager);
        analysis->AddCorrFctn(trueq3dbp_cf);
      }

      if (macro_config.do_kt_trueq3dbp_cf) {
        AliFemtoModelCorrFctnTrueQ3DByParent *kt_trueq3dbp_cf = new AliFemtoModelCorrFctnTrueQ3DByParent(macro_config.q3d_bin_count,
                                                                                              macro_config.q3d_maxq);
        kt_trueq3dbp_cf->SetManager(model_manager);

        AliFemtoKtBinnedCorrFunc *kt_trueq3dbp_cf_container = new AliFemtoKtBinnedCorrFunc("KT_TrueQ3DByP", kt_trueq3dbp_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_trueq3dbp_cf_container->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_trueq3dbp_cf_container);
      }

      if (macro_config.do_truedetadphi_cf) {
        AliFemtoModelCorrFctnDEtaDPhiAK *truedetadphi_cf = new AliFemtoModelCorrFctnDEtaDPhiAK("", macro_config.q3d_bin_count, macro_config.q3d_maxq);
        analysis->AddCorrFctn(truedetadphi_cf);
      }

      if (macro_config.do_kt_truedetadphi_cf) {
        AliFemtoModelCorrFctnDEtaDPhiAK *kt_truedetadphi_cf =
          new AliFemtoModelCorrFctnDEtaDPhiAK("",
                                              macro_config.q3d_bin_count,
                                              macro_config.q3d_maxq);
        kt_truedetadphi_cf->ConnectToManager(model_manager);

        AliFemtoKtBinnedCorrFunc *kt_truedetadphi_cfs = new AliFemtoKtBinnedCorrFunc("KT_TrueDetaDphi", kt_truedetadphi_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_truedetadphi_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_truedetadphi_cfs);
      }

      if (macro_config.do_detadphi_simple_cf) {
        AliFemtoCorrFctnDEtaDPhiSimple *detadphi_simple_cf = new AliFemtoCorrFctnDEtaDPhiSimple("", macro_config.delta_phi_bin_count, macro_config.delta_eta_bin_count);
        detadphi_simple_cf->SetReadHiddenInfo(analysis_config.is_mc_analysis);
        analysis->AddCorrFctn(detadphi_simple_cf);
      }

      if (macro_config.do_kt_detadphi_simple_cf) {
        AliFemtoCorrFctnDEtaDPhiSimple *kt_detadphi_simple_cf = new AliFemtoCorrFctnDEtaDPhiSimple("", macro_config.delta_phi_bin_count, macro_config.delta_eta_bin_count);
        kt_detadphi_simple_cf->SetReadHiddenInfo(analysis_config.is_mc_analysis);

        AliFemtoKtBinnedCorrFunc *kt_detadphi_simple_cfs = new AliFemtoKtBinnedCorrFunc("KT_DetaDphiSimple", kt_detadphi_simple_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_detadphi_simple_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_detadphi_simple_cfs);
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

        AliFemtoKtBinnedCorrFunc *kt_detadphi_star_cfs = new AliFemtoKtBinnedCorrFunc("KT_DEtaDPhiStar", kt_starcf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float low = macro_config.kt_ranges[kt_idx],
                high = macro_config.kt_ranges[kt_idx+1];
          kt_detadphi_star_cfs->AddKtRange(low, high);
        }
        analysis->AddCorrFctn(kt_detadphi_star_cfs);
      }

      if (macro_config.do_ylm_cf) {
        AliFemtoCorrFctnDirectYlm *ylm_cf = new AliFemtoCorrFctnDirectYlm(
            "AliFemtoCorrFctnDirectYlm",
            macro_config.ylm_max_l,
            macro_config.ylm_ibin, // nbinssh,
            macro_config.ylm_vmin,
            macro_config.ylm_vmax,
            macro_config.ylm_useLCMS
          );
        analysis->AddCorrFctn(ylm_cf);
      }

      if (macro_config.do_kt_ylm_cf) {
        AliFemtoCorrFctnDirectYlm *ylm_cf = new AliFemtoCorrFctnDirectYlm(
            "",
            macro_config.ylm_max_l,
            macro_config.ylm_ibin, // nbinssh,
            macro_config.ylm_vmin,
            macro_config.ylm_vmax,
            macro_config.ylm_useLCMS
          );

        AliFemtoKtBinnedCorrFunc *ylm_cfs = new AliFemtoKtBinnedCorrFunc("KT_DirectYlm", ylm_cf);

        for (size_t kt_idx=0; kt_idx < macro_config.kt_ranges.size(); kt_idx += 2) {
          float lo = macro_config.kt_ranges[kt_idx],
                hi = macro_config.kt_ranges[kt_idx+1];
          ylm_cfs->AddKtRange(lo, hi);
        }
        analysis->AddCorrFctn(ylm_cfs);
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
      if (rangeend == -1) {
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
        TString prev = TString::Format("%0.2d", subrange_it->String().Atoi());
        while ((subrange_it = (TObjString *)next_subrange())) {
          TString next = TString::Format("%0.2d", subrange_it->String().Atoi());

          cmd = macro_varname + "->centrality_ranges.push_back(" + prev + ");";
          gROOT->ProcessLineFast(cmd);

          cmd = macro_varname + "->centrality_ranges.push_back(" + next + ");";
          gROOT->ProcessLineFast(cmd);
          prev = next;
        }
      }
    }
    continue;

    case '(':
    {
      UInt_t rangeend = line.Index(")");
      if (rangeend == -1) {
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
        cmd = macro_varname + "->pair_codes.push_back(1)";
      }
      else if (line == "+m") {
        cmd = macro_varname + "->pair_codes.push_back(0)";
      }
      else if (line == "+pp") {
        cmd = macro_varname + "->pair_codes.push_back(3)";
      }
      else if (line == "+mm") {
        cmd = macro_varname + "->pair_codes.push_back(4)";
      }
      else if (line == "+pm") {
        cmd = macro_varname + "->pair_codes.push_back(2)";
      }
      else if (line == "+mp") {
        cmd = macro_varname + "->pair_codes.push_back(5)";
      }
      else {
        continue;
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
