///
/// \file AliFemtoAnalysisPionPion.cxx
///

#include "AliFemtoAnalysisPionPion.h"

#include "AliFemtoCutMonitorPionPion.h"

#include "AliESDtrack.h"

// event readers
#include "AliFemtoEventReaderAODMultSelection.h"

// event cuts
#include "AliFemtoBasicEventCut.h"
#include "AliFemtoEventCutCentrality.h"

// track cuts
#include "AliFemtoV0TrackCut.h"
#include "AliFemtoBasicTrackCut.h"
#include "AliFemtoESDTrackCut.h"

// pair cuts
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoPairCutDetaDphi.h"
#include "AliFemtoPairCutRejectAll.h"

// correlation functions
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoCorrFctnDPhiStarDEta.h"
#include "AliFemtoCorrFctnDirectYlm.h"
#include "AliFemtoCorrFctn3DLCMSSym.h"
#include "AliFemtoModelCorrFctnTrueQ3D.h"

#include "AliFemtoAnalysisPionPionCuts.h"

#include "AliFemtoAnalysisPionPionObjectConstructor.h"

#include <TROOT.h>
#include <TInterpreter.h>

#include <mutex>
#include <sstream>
#include <utility>
#include <cassert>
#include <typeinfo>
#include <random>


static const double PionMass = 0.13956995;
static const int UNKNOWN_CHARGE = -9999;


template<>
struct Configuration<AliFemtoAnalysisPionPion> {

  TString name;
  AliFemtoAnalysisPionPion::AnalysisParams analysis_params;
  AliFemtoAnalysisPionPion::CutParams cut_params;

  Configuration()
    : name("PionFemtoAnalysis")
    , analysis_params()
    , cut_params()
    {
    }

  Configuration(AliFemtoConfigObject cfg)
    : name(cfg.pop_str("name", "PionFemtoAnalysis"))
    , analysis_params()
    , cut_params()
    {
    }

  operator AliFemtoAnalysisPionPion*() const
    { return into_ptr(); }

  AliFemtoAnalysisPionPion* into_ptr() const
    {
      auto *result = new AliFemtoAnalysisPionPion(name.Data(),
                                                  analysis_params,
                                                  cut_params);
      return result;
    }
};


TString AliFemtoAnalysisPionPion::make_random_string(const TString &prefix)
{
  std::random_device rd;
  std::uniform_int_distribution<char> dist('a', 'z');
  TString result = prefix;

  int N = std::uniform_int_distribution<>(6, 12)(rd);
  for (int n = 0; n < N; ++n) {
    result += dist(rd);
  }
  return result;
}


AliFemtoAnalysisPionPion::AnalysisParams::AnalysisParams()
: TNamed(AliFemtoAnalysisPionPion::make_random_string("analysis_").Data(), "Analysis Parameters")
, vertex_bins(1), vertex_min(-10.0), vertex_max(10.0)
, mult_bins(30), mult_min(0.0), mult_max(6480.0)
// , pion_type_1(kNone)
// , pion_type_2(kNone)
, pion_type_1(kPiPlus)
, pion_type_2(kPiPlus)
, num_events_to_mix(6)
, min_coll_size(15)
, verbose(kFALSE)
, enable_pair_monitors(kTRUE)
, group_output_objects(kTRUE)
, output_settings(kFALSE)
, is_mc_analysis(kFALSE)
{
}

void
AliFemtoAnalysisPionPion::AnalysisParams::calc_automult(const AliFemtoAnalysisPionPion::CutParams &cut)
{
  const double
    centlo = cut.event_centrality.first,
    centhi = cut.event_centrality.second,

    bin_width = (mult_max - mult_min) / mult_bins,

    // magic numbers fitting the exp curve
    min_mult = 4293.2 * std::exp(-centhi / 24.01),
    max_mult = 4293.2 * std::exp(-centlo / 24.01),

    mult_min = std::max(0.0, min_mult - 1000);
    mult_max = max_mult + 2000;

  // keep the same bin width
  mult_bins = (mult_max - mult_min) / bin_width;
}

#include <initializer_list>

// typedef Float_t RangeF_t[2];
typedef std::pair<Float_t, Float_t> RangeF_t;

/// Structure for storing all configuration information regarding the
/// analysis event-cut. This is hidden as an implementation detail due
/// the build system not playing well with nested structures; the optimal
/// situation would be to define this structure as a component of the
/// AliFemtoAnalysisPionPion::CutParams structure.
///
///

template <>
AliFemtoConfigObject AliFemtoAnalysisPionPion::GetConfigurationOf<AliFemtoEventCut>(const AliFemtoEventCut &cut)
{
  return AliFemtoConfigObject::From(cut);
}

template <>
AliFemtoConfigObject AliFemtoAnalysisPionPion::GetConfigurationOf<AliFemtoTrackCut>(const AliFemtoTrackCut &cut)
{
  return AliFemtoConfigObject::From(cut);
}


template <>
AliFemtoConfigObject AliFemtoAnalysisPionPion::GetConfigurationOf<AliFemtoPairCut>(const AliFemtoPairCut &cut)
{
  return AliFemtoConfigObject::From(cut);
}

extern template AliFemtoEventCutCentrality* AliFemtoConfigObject::Into<AliFemtoEventCutCentrality>(bool);
extern template AliFemtoBasicEventCut* AliFemtoConfigObject::Into<AliFemtoBasicEventCut>(bool);


/// Configuration for creating a particle cut specifically for usage
/// with pion values.
///
struct CutConfig_Pion {

    RangeF_t pt = {0.15, 2.0},
             eta = {-0.8, 0.8},
             DCA = {0.5, 4.0},
             nSigma = {-3.0, 3.0};

    Int_t charge = 1;
    UInt_t min_tpc_ncls = 80;
    UInt_t min_its_ncls = 0;

    ULong_t status = 0;

    /// momentum to begin using TOF sigma in sigma calculation
    Float_t tpctof_limit = 0.5,
            tof_kaon_reject_pmin = 0.6,
            tof_sigma = NAN,
            tof_kaon_rejection_sigma = 4.0,
            tof_kaon_rejection_pmin = 0.6,

            tof_proton_reject_sigma = 5.0,
            tof_proton_reject_pmin = 0.6,

            electron_rejection_sigma = 3.0;

    Float_t max_impact_xy = .20,
            max_impact_z = .15,
            min_tpc_chi_ndof = 0.0,
            max_tpc_chi_ndof = 2.0,
            max_its_chi_ndof = 5.0,
            sigma = 3.0;

    Bool_t remove_negative_label = kFALSE,
           use_tpctof = kTRUE,
           tof_required = kTRUE,
           remove_kinks = kTRUE;

    Int_t ideal_pid = 211;

    /// default constructor required to use default initialized members
    CutConfig_Pion(){};
};

/// Configuration for the pair
struct CutConfig_Pair {
  Bool_t use_avgsep = { kFALSE };
  Float_t min_delta_eta = { 0.01 },
          min_delta_phi = { 0.033 },
          phi_star_radius = { 1.1 };

  Float_t min_avgsep = { 3.0 };

  Float_t ee_min = { 0.0 };

  Float_t max_share_fraction = { 0.05 },
          max_share_quality = { 1.0 };

  Bool_t remove_same_label { kFALSE },
         TPCOnly { kTRUE };

  Int_t algorithm_code { 2 };

  CutConfig_Pair(){};
};

static const Configuration<AliFemtoBasicEventCut> default_event;
static const CutConfig_Pion default_pion;
static const CutConfig_Pair default_pair;

const AliFemtoAnalysisPionPion::PionType
  default_PionType = AliFemtoAnalysisPionPion::kPiPlus;


static const AliFemtoAnalysisPionPion::AnalysisParams
analysis_params_from_pion_types(AliFemtoAnalysisPionPion::PionType one, AliFemtoAnalysisPionPion::PionType two) {
    auto result = AliFemtoAnalysisPionPion::DefaultConfig();
    result.pion_type_1 = one;
    result.pion_type_2 = two;
    return result;
}


AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion():
  AliFemtoAnalysisPionPion("AliFemtoAnalysisPionPion")
{
}

AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const char *name):
  AliFemtoAnalysisPionPion(name, default_PionType, default_PionType)
{
}


AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const char *name,
                                                   const PionType pion_1,
                                                   const PionType pion_2):
  AliFemtoAnalysisPionPion(name, pion_1, pion_2, DefaultCutConfig())
{
}


AliFemtoAnalysisPionPion
  ::AliFemtoAnalysisPionPion(const char *name,
                             const PionType pion_1,
                             const PionType pion_2,
                             const CutParams &cut_params):
    AliFemtoAnalysisPionPion(name, analysis_params_from_pion_types(pion_1, pion_2), cut_params)
{
}


AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const char *name,
                                                   const AnalysisParams &params,
                                                   const CutParams &cut_params):
  AliFemtoVertexMultAnalysis(params.vertex_bins,
                             params.vertex_min,
                             params.vertex_max,
                             params.mult_bins,
                             params.mult_min,
                             params.mult_max)
  , fAnalysisName(name)
  , fPionType_1(params.pion_type_1)
  , fPionType_2(params.pion_type_2)
  , fGroupOutputObjects(params.group_output_objects)
  , fOutputSettings(params.output_settings)
  , fMCAnalysis(params.is_mc_analysis)
{
  SetVerboseMode(params.verbose);
  SetEnablePairMonitors(params.enable_pair_monitors);
  SetNumEventsToMix(params.num_events_to_mix);
  SetMinSizePartCollection(params.min_coll_size);

  if (fFirstParticleCut == nullptr) {
    SetFirstParticleCut(BuildPionCut1(cut_params));
  }

  if (fPionType_2 != kNone && fSecondParticleCut == nullptr) {
    SetSecondParticleCut(BuildPionCut2(cut_params));
  } else if (fPionType_2 == kNone) {
    fSecondParticleCut = fFirstParticleCut;
  }

  if (fEventCut == nullptr) {
    SetEventCut(BuildEventCut(cut_params));
  }

  if (fPairCut == nullptr) {
    // if (cut_params.cuts_use_attrs) {
    // SetFirstParticleCut(new AliFemtoPairCutPionPionAKAvgSep());
    // } else {
    SetPairCut(BuildPairCut(cut_params));
    // }
  }
}

AliFemtoAnalysisPionPion::AliFemtoAnalysisPionPion(const AliFemtoConfigObject &cfg)
  : AliFemtoVertexMultAnalysis(cfg.get("mix_vertex_z_bins", (ULong64_t)1), cfg.get("mix_vertex_z_min", -10.0), cfg.get("mix_vertex_z_max", 10.0),
                               cfg.get("mix_mult_bins", (ULong64_t)30), cfg.get("mix_mult_min", 0.0), cfg.get("mix_mult_max", 6480.0))
  , fAnalysisName("AliFemtoAnalysisPionPion")
  , fPionType_1(kPiPlus)
  , fPionType_2(kNone)
  , fGroupOutputObjects(true)
  , fOutputSettings(false)
  , fMCAnalysis(false)
{
  std::cout << "Building AliFemtoAnalysisPionPion from configuration object\n" << cfg.Stringify() << "\n";
  cfg.find_and_load("name", fAnalysisName);

  auto defaults = DefaultConfig();

  SetNumEventsToMix(cfg.get("events_to_mix", 6ull));
  SetVerboseMode(cfg.get("verbose", false));
  SetEnablePairMonitors(cfg.get("pair_monitors", true));
  SetMinSizePartCollection(cfg.get("collection_size_min", 10ull));

  CutParams default_params;

  AliFemtoConfigObject subcfg;

  {
    auto *cut = cfg.find_and_load("event_cut", subcfg)
              ? ConstructEventCut(subcfg)
              : BuildEventCut(default_params);
    SetEventCut(cut);
  }

  {
    auto *cut = cfg.find_and_load("track_cut", subcfg)
              ? ConstructParticleCut(subcfg)
              : BuildPionCut1(default_params);
    SetFirstParticleCut(cut);
    SetSecondParticleCut(cut);  // Identical pions
    if (auto *c = dynamic_cast<AliFemtoTrackCutPionPionAK*>(cut)) {
      if (c->charge == -1) {
        fPionType_1 = kPiMinus;
      }
    } else if (auto *c = dynamic_cast<AliFemtoESDTrackCut*>(cut)) {
      if (c->GetCharge() == -1) {
        fPionType_1 = kPiMinus;
      }
    }
  }

  {
    auto *cut = cfg.find_and_load("pair_cut", subcfg)
              ? ConstructPairCut(subcfg)
              : BuildPairCut(default_params);
    SetPairCut(cut);
  }
}

AliFemtoAnalysisPionPion::AnalysisParams
AliFemtoAnalysisPionPion::DefaultConfig()
{
  return AliFemtoAnalysisPionPion::AnalysisParams();
}


AliFemtoAnalysisPionPion::CutParams::CutParams()
  : TNamed(AliFemtoAnalysisPionPion::make_random_string("cut_").Data(), "Cut Params")
  , cuts_use_attrs(true)
  , mc_pion_only(false)
  , mc_nonpion_only(false)
  , event_use_basic(false)
  , event_mult(default_event.multiplicity)
  , event_centrality(default_event.centrality)
  , event_vertex_z(default_event.vertex_z)
  , event_ep_vzero(default_event.ep_psi)
  , event_TriggerSelection(default_event.trigger_selection)
  , event_AcceptBadVertex(default_event.accept_bad_vertex)
  , event_AcceptOnlyPhysics(default_event.accept_bad_vertex)
  , event_zdc_part(0.0)

    // Pion 1
  , pion_1_pt(default_pion.pt)
  , pion_1_eta(default_pion.eta)
  , pion_1_DCA(default_pion.DCA)

  , pion_1_status(default_pion.status)
  // , default_pion.nSigma.first
  // , default_pion.nSigma.second
  , pion_1_tpc_sigma(default_pion.sigma)
  , pion_1_tof_sigma(default_pion.tof_sigma)
  , pion_1_tof_required(default_pion.tof_required)

  , pion_1_tof_kreject_sigma(default_pion.tof_kaon_rejection_sigma)
  , pion_1_tof_kreject_pmin(default_pion.tof_kaon_reject_pmin)

  , pion_1_tof_preject_sigma(default_pion.tof_proton_reject_sigma)
  , pion_1_tof_preject_pmin(default_pion.tof_proton_reject_pmin)

  , pion_1_tpc_ereject_sigma(default_pion.electron_rejection_sigma)

  , pion_1_max_impact_xy(default_pion.max_impact_xy)
  , pion_1_max_impact_z(default_pion.max_impact_z)
  , pion_1_min_tpc_chi_ndof(default_pion.min_tpc_chi_ndof)
  , pion_1_max_tpc_chi_ndof(default_pion.max_tpc_chi_ndof)
  , pion_1_max_its_chi_ndof(default_pion.max_its_chi_ndof)

  , pion_1_min_tpc_ncls(default_pion.min_tpc_ncls)
  , pion_1_min_its_ncls(default_pion.min_its_ncls)

  , pion_1_remove_kinks(default_pion.remove_kinks)
  , pion_1_rm_neg_lbl(default_pion.remove_negative_label)
  , pion_1_use_tpctof(default_pion.use_tpctof)

  , pion_1_ideal_pid(default_pion.ideal_pid)

    // Pion 2
  , pion_2_PtMin(default_pion.pt.first)
  , pion_2_PtMax(default_pion.pt.second)
  , pion_2_EtaMin(default_pion.eta.first)
  , pion_2_EtaMax(default_pion.eta.second)
  , pion_2_DCAMin(default_pion.DCA.first)
  , pion_2_DCAMax(default_pion.DCA.second)

  , pion_2_max_impact_xy(default_pion.max_impact_xy)
  , pion_2_max_impact_z(default_pion.max_impact_z)
  , pion_2_max_tpc_chi_ndof(default_pion.max_tpc_chi_ndof)
  , pion_2_max_its_chi_ndof(default_pion.max_its_chi_ndof)

  , pion_2_min_tpc_ncls(default_pion.min_tpc_ncls)
  , pion_2_remove_kinks(default_pion.remove_kinks)
  , pion_2_set_label(default_pion.remove_negative_label)

    // Pair
  , pair_use_avgsep(default_pair.use_avgsep)
  , pair_TPCOnly(default_pair.TPCOnly)

  , pair_min_avgsep(default_pair.min_avgsep)

  , pair_delta_eta_min(default_pair.min_delta_eta)
  , pair_delta_phi_min(default_pair.min_delta_phi)
  , pair_phi_star_radius(default_pair.phi_star_radius)

  , pair_ee_min(default_pair.ee_min)
  , pair_max_share_quality(default_pair.max_share_quality)
  , pair_max_share_fraction(default_pair.max_share_fraction)
  , pair_remove_same_label(default_pair.remove_same_label)
  , pair_algorithm(default_pair.algorithm_code)
{}


AliFemtoAnalysisPionPion::CutParams
AliFemtoAnalysisPionPion::DefaultCutConfig()
{
  AliFemtoAnalysisPionPion::CutParams params;

  // sanity checks
  assert(params.event_mult == default_event.multiplicity);
  assert(params.pion_1_pt.first == default_pion.pt.first && params.pion_1_pt.second == default_pion.pt.second);
  assert(params.pair_TPCOnly == default_pair.TPCOnly);
  // assert(params.pair_TPCExitSepMin == default_pair_TPCExitSepMin);
  assert(params.pair_delta_eta_min == default_pair.min_delta_eta);
  assert(params.pair_max_share_fraction == default_pair.min_delta_phi);
  assert(params.pair_remove_same_label == default_pair.remove_same_label);

  return params;
}

AliFemtoTrackCut*
AliFemtoAnalysisPionPion::BuildPionCut1(const CutParams &p) const
{
  const int charge = (fPionType_1 == kPiMinus) ? -1 :
                      (fPionType_1 == kPiPlus) ? +1 :
                                                UNKNOWN_CHARGE;
  if (charge == UNKNOWN_CHARGE) {
    std::cerr << "E-AliFemtoAnalysisPionPion::BuildPionCut: Invalid pion type: '" << fPionType_1 << "'\n";
  }

  if (p.cuts_use_attrs) {
    auto *cut = p.mc_nonpion_only ? new AliFemtoTrackCutPionPionMisidentAK()
              :    p.mc_pion_only ? new AliFemtoTrackCutPionPionIdealAK()
                                  : new AliFemtoTrackCutPionPionAK();
    cut->pt_range = p.pion_1_pt;
    cut->eta_range = p.pion_1_eta;
    cut->status = p.pion_1_status;
    cut->ncls_tpc_min = p.pion_1_min_tpc_ncls;
    cut->ncls_its_min = p.pion_1_min_its_ncls;
    cut->charge = charge;
    cut->max_xy = p.pion_1_max_impact_xy;
    cut->max_z = p.pion_1_max_impact_z;

    cut->tpc_sigma_pion = p.pion_1_tpc_sigma;
    cut->tof_sigma_pion = p.pion_1_tof_sigma;
    cut->require_tof = p.pion_1_tof_required;

    cut->tof_kaon_reject_sigma = p.pion_1_tof_kreject_sigma;
    cut->tof_kaon_momentum_limit = p.pion_1_tof_kreject_pmin;
    // cut->tof_sigma_pion = p.pion_1_tof_sigma;

    cut->tof_proton_reject_sigma = p.pion_1_tof_preject_sigma;
    cut->tof_proton_momentum_limit = p.pion_1_tof_preject_pmin;
    cut->electron_tpc_sigma_min = p.pion_1_tpc_ereject_sigma;
    cut->rchi2_tpc_min = p.pion_1_min_tpc_chi_ndof;
    cut->rchi2_tpc_max = p.pion_1_max_tpc_chi_ndof;
    cut->rchi2_its_max = p.pion_1_max_its_chi_ndof;
    cut->remove_neg_label = p.pion_1_rm_neg_lbl;

    if (auto *ideal_cut = dynamic_cast<AliFemtoTrackCutPionPionIdealAK*>(cut)) {
      ideal_cut->ideal_pid = p.pion_1_ideal_pid;
    }

    return cut;
  }


//   AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();
//   cut->SetNSigmaPion(p.pion_1_NSigmaMin, p.pion_1_NSigmaMax);
//   cut->SetDCA(p.pion_1_DCAMin, p.pion_1_DCAMax);

  AliFemtoESDTrackCut *cut = new AliFemtoESDTrackCut();
  cut->SetCharge(charge);
  cut->SetMass(PionMass);
  cut->SetPt(p.pion_1_pt.first, p.pion_1_pt.second);
  cut->SetEta(p.pion_1_eta.first, p.pion_1_eta.second);
  cut->SetRapidity(p.pion_1_eta.first, p.pion_1_eta.second);
  cut->SetMostProbablePion();
  cut->SetNsigma(p.pion_1_tpc_sigma);
  cut->SetNsigmaTPCTOF(p.pion_1_use_tpctof);
//   cut->SetStatus(AliESDtrack::kTPCrefit | AliESDtrack::kITSrefit);

  /// Settings for TPC-Inner Runmode
  cut->SetStatus(AliESDtrack::kTPCin);
  cut->SetminTPCncls(p.pion_1_min_tpc_ncls);
  cut->SetRemoveKinks(p.pion_1_remove_kinks);
  cut->SetLabel(p.pion_1_rm_neg_lbl);
  cut->SetMaxTPCChiNdof(p.pion_1_max_tpc_chi_ndof);
  cut->SetMaxITSChiNdof(p.pion_1_max_its_chi_ndof);
  cut->SetMaxImpactXY(p.pion_1_max_impact_xy);
  cut->SetMaxImpactZ(p.pion_1_max_impact_z);

  return cut;
}

AliFemtoTrackCut*
AliFemtoAnalysisPionPion::BuildPionCut2(const CutParams &p) const
{
  const int charge = (fPionType_2 == kPiMinus) ? -1 :
                      (fPionType_2 == kPiPlus) ? +1 :
                                                UNKNOWN_CHARGE;
  if (charge == UNKNOWN_CHARGE) {
    std::cerr << "E-AliFemtoAnalysisPionPion::BuildPionCut: Invalid pion type: '" << fPionType_2 << "'\n";
  }

//   AliFemtoBasicTrackCut *cut = new AliFemtoBasicTrackCut();
//   cut->SetCharge(charge);
//   cut->SetMass(PionMass);
//   cut->SetNSigmaPion(p.pion_2_NSigmaMin, p.pion_2_NSigmaMax);
//   cut->SetPt(p.pion_2_PtMin, p.pion_2_PtMax);
//   cut->SetRapidity(p.pion_2_EtaMin, p.pion_2_EtaMax);
//   cut->SetDCA(p.pion_2_DCAMin, p.pion_2_DCAMax);
  // AliFemtoConfigObject *cfg;
  // auto cut = cfg->Construct<AliFemtoESDTrackCut>();

  AliFemtoESDTrackCut *cut = new AliFemtoESDTrackCut();
  cut->SetCharge(charge);
  cut->SetMass(PionMass);
  cut->SetEta(p.pion_2_EtaMin, p.pion_2_EtaMax);
  cut->SetMostProbablePion();

  /// Settings for TPC-Inner Runmode
  cut->SetStatus(AliESDtrack::kTPCin);

  cut->SetminTPCncls(p.pion_2_min_tpc_ncls);
  cut->SetRemoveKinks(p.pion_2_remove_kinks);
  cut->SetLabel(p.pion_2_set_label);
  cut->SetMaxTPCChiNdof(p.pion_2_max_tpc_chi_ndof);
  cut->SetMaxImpactXY(p.pion_2_max_impact_xy);
  cut->SetMaxImpactZ(p.pion_2_max_impact_z);

  return cut;
}

AliFemtoEventCut*
AliFemtoAnalysisPionPion::BuildEventCut(const AliFemtoAnalysisPionPion::CutParams& p) const
{
  if (p.cuts_use_attrs) {
    auto *cut = new AliFemtoEventCutPionPionAK();

    cut->cent_range = {p.event_centrality.first, p.event_centrality.second};
    cut->zvert_range = {p.event_vertex_z.first, p.event_vertex_z.second};
    cut->trigger = p.event_TriggerSelection;
    cut->zdc_participants_min = p.event_zdc_part;

    return cut;
  }

  if (p.event_use_basic) {
    AliFemtoBasicEventCut *cut = new AliFemtoBasicEventCut();

    cut->SetEventMult(p.event_mult.first, p.event_mult.second);
    cut->SetVertZPos(p.event_vertex_z.first, p.event_vertex_z.second);
    cut->SetEPVZERO(p.event_ep_vzero.first, p.event_ep_vzero.second);
    cut->SetTriggerSelection(p.event_TriggerSelection);
    cut->SetAcceptBadVertex(p.event_AcceptBadVertex);
    return cut;
  }

  AliFemtoEventCutCentrality *cut = new AliFemtoEventCutCentrality();
  cut->SetCentralityRange(p.event_centrality.first, p.event_centrality.second);
  cut->SetZPosRange(p.event_vertex_z.first, p.event_vertex_z.second);
  cut->SetEPVZERO(p.event_ep_vzero.first, p.event_ep_vzero.second);
  cut->SetTriggerSelection(p.event_TriggerSelection);

//   if (!fMCAnalysis) {
//  cut->SetAcceptOnlyPhysics(p.event_AcceptOnlyPhysics);
//   }

  return cut;
}



AliFemtoPairCut*
AliFemtoAnalysisPionPion::BuildPairCut(const CutParams &p) const
{
  if (p.cuts_use_attrs) {
    if (p.mc_nonpion_only) {
      return new AliFemtoPairCutRejectAll();
    }
    if (p.pair_use_avgsep) {
      auto *pair_cut = new AliFemtoPairCutPionPionAKAvgSep();
      pair_cut->remove_same_label = p.pair_remove_same_label;
      pair_cut->share_quality_max = p.pair_max_share_quality;
      pair_cut->share_fraction_max = p.pair_max_share_fraction;
      pair_cut->ee_minv_min = p.pair_ee_min;

      pair_cut->avgsep_min = p.pair_min_avgsep;

      return pair_cut;

    } else {
      auto *pair_cut = new AliFemtoPairCutPionPionAKDetaDphi();
      pair_cut->remove_same_label = p.pair_remove_same_label;
      pair_cut->share_quality_max = p.pair_max_share_quality;
      pair_cut->share_fraction_max = p.pair_max_share_fraction;
      pair_cut->ee_minv_min = p.pair_ee_min;

      pair_cut->delta_eta_min = p.pair_delta_eta_min;
      pair_cut->delta_phistar_min = p.pair_delta_phi_min;
      pair_cut->phistar_radius = p.pair_phi_star_radius;
      return pair_cut;
    }
  }

  /*
  AliFemtoPairCutAntiGamma *cut = new AliFemtoPairCutAntiGamma();

  // cut->SetTPCOnly(p.pair_TPCOnly);
  // cut->SetTPCExitSepMinimum(p.pair_TPCExitSepMin);
  // cut->SetMinAvgSeparation(0, p.pair_MinAvgSeparationPos);
  // cut->SetMinAvgSeparation(1, p.pair_MinAvgSeparationNeg);
  */


  AliFemtoPairCutDetaDphi *cut = new AliFemtoPairCutDetaDphi(p.pair_delta_eta_min,
                                                             p.pair_delta_phi_min,
                                                             p.pair_phi_star_radius);

  cut->SetCutTechnique(AliFemtoPairCutDetaDphi::Quad);
  cut->SetShareQualityMax(p.pair_max_share_quality);
  cut->SetShareFractionMax(p.pair_max_share_fraction);
  cut->SetRemoveSameLabel(p.pair_remove_same_label);

  return cut;
}

void AliFemtoAnalysisPionPion::AddStanardCutMonitors()
{
  const bool identical = AnalyzeIdenticalParticles();

  if (fEventCut) {
    fEventCut->AddCutMonitor(new AliFemtoCutMonitorPionPion::Event(true, identical, fMCAnalysis, false),
                             new AliFemtoCutMonitorPionPion::Event(false, identical, fMCAnalysis, false));
  } else {
    std::cout << "E-AliFemtoAnalysisPionPion::AddStanardCutMonitors: NO fEventCut!\n";
  }

  const TString p1_type_str = (fPionType_1 == kPiPlus)  ? "Pi+"
                            : (fPionType_1 == kPiMinus) ? "Pi-"
                                                        : "ERROR";

  if (fFirstParticleCut) {

    Long64_t dca_track_method = 1;
    fEventReaderCfg.find_and_load("dca_method", dca_track_method);

    bool use_wide_bins = dca_track_method == 0;

    auto *pass_cut = new AliFemtoCutMonitorPionPion::Pion(true, p1_type_str, fMCAnalysis, false, use_wide_bins),
         *fail_cut = new AliFemtoCutMonitorPionPion::Pion(false, p1_type_str, fMCAnalysis, false, use_wide_bins);

    fail_cut->SetCharge(fPionType_1 == kPiPlus ? 1 : -1);
    fFirstParticleCut->AddCutMonitor(pass_cut, fail_cut);
  }

  if (!identical) {
    const TString p2_type_str = (fPionType_2 == kPiPlus)  ? "Pi+"
                              : (fPionType_2 == kPiMinus) ? "Pi-"
                                                          : "ERROR";

    fSecondParticleCut->AddCutMonitor(new AliFemtoCutMonitorPionPion::Pion(true, p2_type_str, fMCAnalysis),
                                      new AliFemtoCutMonitorPionPion::Pion(false, p2_type_str, fMCAnalysis));

    // Non-Identical Pion Pairs
    if (fPairCut) {
      const TString pair_type_str = p1_type_str + "/" + p2_type_str;
      fPairCut->AddCutMonitor(new AliFemtoCutMonitorPionPion::Pair(true, pair_type_str, fMCAnalysis),
                              new AliFemtoCutMonitorPionPion::Pair(false, pair_type_str, fMCAnalysis));
    }
  }
  // Identical Pion Pairs
  else if (fPairCut) {
    const TString pair_type_str = "(Identical) " + p1_type_str;
    fPairCut->AddCutMonitor(new AliFemtoCutMonitorPionPion::Pair(true, pair_type_str, fMCAnalysis),
                            new AliFemtoCutMonitorPionPion::Pair(false, pair_type_str, fMCAnalysis));
  }

}

static TObjArray* GetPassFailOutputList(const TString &name,
                                        AliFemtoCutMonitorHandler *handler)
{
  AliFemtoCutMonitorCollection *passing_monitors = handler->PassMonitorColl(),
                               *failing_monitors = handler->FailMonitorColl();

  TObjArray *res = new TObjArray(),
            *passout = new TObjArray(),
            *failout = new TObjArray();

  res->SetName(name);
  passout->SetName("pass");
  failout->SetName("fail");

  res->Add(passout);
  res->Add(failout);

  for (const auto &cut_monitor : *passing_monitors) {
    TList *output_list = cut_monitor->GetOutputList();
    passout->AddAll(output_list);
    delete output_list;
  }

  for (const auto &cut_monitor : *failing_monitors) {
    TList *output_list = cut_monitor->GetOutputList();
    failout->AddAll(output_list);
    delete output_list;
  }

  return res;
}

void
AliFemtoAnalysisPionPion::StoreEventReaderConfiguration(const AliFemtoEventReader &evreader)
{
  if (auto *aod = dynamic_cast<const AliFemtoEventReaderAOD*>(&evreader)) {

    auto mult = aod->GetUseMultiplicity();

    const TString
      cent_type = mult == AliFemtoEventReaderAOD::kCentrality ? "V0M"
                : mult == AliFemtoEventReaderAOD::kCentralityV0A ? "V0A"
                : mult == AliFemtoEventReaderAOD::kCentralityV0C ? "V0C"
                : mult == AliFemtoEventReaderAOD::kCentralityZNA ? "ZNA"
                : mult == AliFemtoEventReaderAOD::kCentralityZNC ? "ZNC"
                : mult == AliFemtoEventReaderAOD::kCentralityCL0 ? "CL0"
                : mult == AliFemtoEventReaderAOD::kCentralityCL1 ? "CL1"
                                                                 : "??";

    fEventReaderCfg = AliFemtoConfigObject::BuildMap()
          ("filtermask", aod->GetTrackFilter())
          ("vertex_shift", aod->GetPrimaryVertexCorrectionTPCPoints())
          ("centrality", cent_type)
          ("dca_method", aod->GetDCAglobalTrack());
  } else {
    fEventReaderCfg = AliFemtoConfigObject::BuildMap()
          ("_class", "??");
  }
}

AliFemtoConfigObject AliFemtoAnalysisPionPion::GetConfiguration() const
{
  auto event_cut_cfg = GetConfigurationOf(*fEventCut);
  auto track_cut_cfg = GetConfigurationOf(static_cast<const AliFemtoTrackCut&>(*fFirstParticleCut));
  auto pair_cut_cfg = GetConfigurationOf(*fPairCut);

  return AliFemtoConfigObject::BuildMap()
                  ("_class", "AliFemtoAnalysisPionPion")
                  ("is_mc", fMCAnalysis)
                  ("events_to_mix", NumEventsToMix())
                  ("collection_size_min", fMinSizePartCollection)
                  ("mix_vertex_z_bins", fVertexZBins)
                  ("mix_vertex_z_range", std::make_pair(fVertexZ[0], fVertexZ[1]))
                  ("mix_mult_bins", fMultBins)
                  ("mix_mult_range", std::make_pair(fMult[0], fMult[1]))
                  ("event_reader", fEventReaderCfg)
                  ("event_cut", event_cut_cfg)
                  ("track_cut", track_cut_cfg)
                  ("pair_cut", pair_cut_cfg);
}

TList* AliFemtoAnalysisPionPion::GetOutputList()
{
  TList *outputlist = nullptr;
  TSeqCollection *output = nullptr;

  auto write_configobj = [&] (TCollection *output)
    {
      auto cfg = GetConfiguration();
      output->Add(new AliFemtoConfigObject(std::move(cfg)));
    };

  // get "standard" outputs - that's what we output
  if (!fGroupOutputObjects) {
    outputlist = AliFemtoVertexMultAnalysis::GetOutputList();
    output = outputlist;
    write_configobj(output);
  } else {

    output = new TObjArray();
    output->SetName(fAnalysisName);
    write_configobj(output);

    if (fEventCut)
    output->Add(GetPassFailOutputList("Event", EventCut()));

    if (AnalyzeIdenticalParticles()) {
    if (fFirstParticleCut)
      output->Add(GetPassFailOutputList("Tracks", fFirstParticleCut));
    } else {
      output->Add(GetPassFailOutputList("Track1", fFirstParticleCut));
      output->Add(GetPassFailOutputList("Track2", fSecondParticleCut));
    }

    output->Add(GetPassFailOutputList("Pair", PairCut()));

    for (auto &corr_fctn : *fCorrFctnCollection) {
      TList *cf_output = corr_fctn->GetOutputList();
      output->AddAll(cf_output);
      delete cf_output;
    }

    outputlist = new TList();
    outputlist->Add(output);
  }

  // get the list of settings from each cut and correlation function
  if (fOutputSettings) {
    TString settings("== settings ==\n");

    TList *setting_list = ListSettings();
    TIter next_setting(setting_list);
    while (TObject *setting = next_setting()) {
      TObjString *setting_str = (TObjString*)setting;
      settings += setting_str->String() + "\n";
    }
    delete setting_list;
    output->Add(new TObjString(settings));
  }

  return outputlist;
}

TList* AliFemtoAnalysisPionPion::ListSettings()
{
  TList *setting_list = new TList();

  TString prefix = "AliFemtoAnalysisPionPion";

  setting_list->AddVector(

    new TObjString(prefix + Form(".mc_analysis=%d", fMCAnalysis)),
    new TObjString(prefix + Form(".identical_analysis=%d", AnalyzeIdenticalParticles())),
    new TObjString(prefix + Form(".pion_1_type=%d", fPionType_1)),

  nullptr);

  if (!AnalyzeIdenticalParticles()) {
    setting_list->Add(new TObjString(prefix + Form(".pion_2_type=%d", fPionType_2)));
  }

  TList *parent_list = AliFemtoVertexMultAnalysis::ListSettings();

  setting_list->AddAll(parent_list);
  delete parent_list;

  setting_list->SetOwner(kTRUE);
  return setting_list;
}

void AliFemtoAnalysisPionPion::EventBegin(const AliFemtoEvent* ev)
{
  fEventCut->EventBegin(ev);
  fFirstParticleCut->EventBegin(ev);
  if (!AnalyzeIdenticalParticles()) {
    fSecondParticleCut->EventBegin(ev);
  }
  fPairCut->EventBegin(ev);
  for (auto &corr_fctn : *fCorrFctnCollection) {
    corr_fctn->EventBegin(ev);
  }
}

void AliFemtoAnalysisPionPion::EventEnd(const AliFemtoEvent* ev)
{
  fEventCut->EventEnd(ev);
  fFirstParticleCut->EventEnd(ev);
  if (!AnalyzeIdenticalParticles()) {
    fSecondParticleCut->EventEnd(ev);
  }
  fPairCut->EventEnd(ev);
  for (auto &corr_fctn : *fCorrFctnCollection) {
    corr_fctn->EventEnd(ev);
  }
}


AliFemtoEventReader*
AliFemtoAnalysisPionPion::ConstructEventReader(AliFemtoConfigObject cfg)
{
  return cfg.Into<AliFemtoEventReader>();
}

AliFemtoEventCut*
AliFemtoAnalysisPionPion::ConstructEventCut(AliFemtoConfigObject cfg)
{
  auto result = cfg.Into<AliFemtoEventCut>();
  return result;
}

AliFemtoPairCut*
AliFemtoAnalysisPionPion::ConstructPairCut(AliFemtoConfigObject cfg)
{
  auto result = cfg.Into<AliFemtoPairCut>();
  return result;
}

AliFemtoParticleCut*
AliFemtoAnalysisPionPion::ConstructParticleCut(AliFemtoConfigObject cfg)
{
  auto result = cfg.Into<AliFemtoTrackCut>();
  return result;
}

AliFemtoAnalysis*
AliFemtoAnalysisPionPion::BuildAnalysisFromConfiguration(AliFemtoConfigObject cfg)
{
  std::string classname;
  if (!cfg.pop_and_load("_class", classname)) {
    std::cerr << "[AliFemtoAnalysisPionPion::BuildAnalysisFromConfiguration] "
              << "Could not load string-property 'class' from object:\n"
              << cfg.Stringify(true)
	      << "\n";
    return nullptr;
  }

  #define TRY_CONSTRUCTING_CLASS(__name) (classname == #__name) ? (AliFemtoAnalysis*)(Configuration<__name>(cfg))

  AliFemtoAnalysis *result = TRY_CONSTRUCTING_CLASS(AliFemtoAnalysisPionPion)
                          //  : TRY_CONSTRUCTING_CLASS(AliFemtoEventReaderAODMultSelection)
                           : nullptr;
  #undef TRY_CONSTRUCTING_CLASS

  return result;
}

AliFemtoCorrFctn*
AliFemtoAnalysisPionPion::ConstructCorrelationFunction(AliFemtoConfigObject cfg)
{
  std::string classname;
  if (!cfg.pop_and_load("_class", classname)) {
    TString msg = "Could not load string-property 'class' from object:\n" + cfg.Stringify(true);
    std::cerr << "[AliFemtoAnalysisPionPion::ConstructCorrelationFunction] " << msg;
    return nullptr;
  }

  #define TRY_CONSTRUCTING_CLASS(__name) (classname == #__name) ? (AliFemtoCorrFctn*)(Configuration<__name>(cfg))

  AliFemtoCorrFctn *result = TRY_CONSTRUCTING_CLASS(AliFemtoAvgSepCorrFctn)
                           : TRY_CONSTRUCTING_CLASS(AliFemtoModelCorrFctnTrueQ3D)
                           : TRY_CONSTRUCTING_CLASS(AliFemtoCorrFctnDirectYlm)
                           : TRY_CONSTRUCTING_CLASS(AliFemtoCorrFctn3DLCMSSym)
                           : TRY_CONSTRUCTING_CLASS(AliFemtoCorrFctnDPhiStarDEta)
                           : nullptr;
  #undef TRY_CONSTRUCTING_CLASS

  return result;
}
