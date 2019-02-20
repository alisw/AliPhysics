///
/// \file AliFemtoUser/AliFemtoCutAttrTrack.h
///

#pragma once

#ifndef ALIFEMTOCUTATTRTRACK_H
#define ALIFEMTOCUTATTRTRACK_H

#include "AliFemtoConfigObject.h"
#include "AliFemtoTrack.h"
#include "AliESDtrackCuts.h"

#include <utility>

namespace pwgfemto {

template <typename T1, typename T2>
struct AddTrackCutAttrs : public T1, public T2 {

  bool Pass(const AliFemtoTrack &track)
    {
      return T1::Pass(track) && T2::Pass(track);
    }

  AddTrackCutAttrs()
    : T1()
    , T2()
    {}

  AddTrackCutAttrs(AliFemtoConfigObject &cfg)
    : T1(cfg)
    , T2(cfg)
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddTrackCutAttrs() {}
};

/// Splits cuts into "selection" and "cut"
///
/// Allows you to logically separate cutting on the type of tracks
/// (e.g. GetLabel()) vs physics properties (pT, DCA, ...)
///
template <typename T1, typename T2>
struct TrackSelectionCut : public T1, public T2 {

  bool PassSelection(const AliFemtoTrack &track)
    {
      return T1::Pass(track);
    }

  bool PassCut(const AliFemtoTrack &track)
    {
      return T2::Pass(track);
    }

  bool Pass(const AliFemtoTrack &track)
    {
      return T1::Pass(track) && T2::Pass(track);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  TrackSelectionCut()
    : T1()
    , T2()
    {}

  TrackSelectionCut(AliFemtoConfigObject &cfg)
    : T1(cfg)
    , T2(cfg)
    {}

  virtual ~TrackSelectionCut() {}
};


struct TrackCutAttrStatus {
  ULong_t status;

  bool Pass(const AliFemtoTrack &trk) const
    {
      return (trk.Flags() & status) == status;
    }

  TrackCutAttrStatus()
    : status(0)
    {}

  TrackCutAttrStatus(ULong_t st)
    : status(st)
    {}

  TrackCutAttrStatus(AliFemtoConfigObject &cfg)
    : status(cfg.pop_uint("status", 0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("status", status);
    }

  virtual ~TrackCutAttrStatus() {}
};

struct TrackCutAttrRemoveKinks {

  bool remove_kinks;

  bool Pass(const AliFemtoTrack &track)
    {
      return !remove_kinks
          || !(track.KinkIndex(0)
               || track.KinkIndex(1)
               || track.KinkIndex(2));
    }

  TrackCutAttrRemoveKinks()
    : remove_kinks(false)
    {}

  TrackCutAttrRemoveKinks(AliFemtoConfigObject &cfg)
    : remove_kinks(cfg.pop_bool("remove_kinks", false))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("remove_kinks", remove_kinks);
    }

  virtual ~TrackCutAttrRemoveKinks() {}
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrRemoveFakeITS {

  bool remove_its_fake;

  bool Pass(const AliFemtoTrack &track)
    {
      return !remove_its_fake || track.ITSncls() < 0;
    }

  TrackCutAttrRemoveFakeITS()
    : remove_its_fake(false)
    {}

  TrackCutAttrRemoveFakeITS(AliFemtoConfigObject &cfg)
    : remove_its_fake(cfg.pop_bool("remove_its_fake", false))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("remove_its_fake", remove_its_fake);
    }

  virtual ~TrackCutAttrRemoveFakeITS() {}
};


/// \class TrackCutAttrImpact
/// \brief Cut tracks outside of z and radial
///
struct TrackCutAttrImpact {

  double max_xy;
  double max_z;

  bool Pass(const AliFemtoTrack &track)
    {
      return max_z <= std::abs(track.ImpactZ())
          && max_xy <= track.ImpactD();
    }

  TrackCutAttrImpact()
    // : xy_range(0, 0.25)
    : max_xy(0.25)
    , max_z(0.3)
    {}

  TrackCutAttrImpact(AliFemtoConfigObject &cfg)
    : max_xy(cfg.pop_float("max_xy", 0.25))
    , max_z(cfg.pop_float("max_z", 0.3))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("max_z", max_z);
      cfg.insert("max_xy", max_xy);
    }

  virtual ~TrackCutAttrImpact() {}
};


/// \class TrackCutAttrChi2ITS
/// \brief Cut on reduced chi2 of results in ITS
///
struct TrackCutAttrChi2ITS {

  double max_rchi2_its;

  bool Pass(const AliFemtoTrack &track)
    {
      return (track.ITSncls() > 0)
              && (track.ITSchi2() / track.ITSncls()) > max_rchi2_its;
    }

  TrackCutAttrChi2ITS()
    : max_rchi2_its(3.0)
    {}

  TrackCutAttrChi2ITS(AliFemtoConfigObject &cfg)
    : max_rchi2_its(cfg.pop_float("max_rchi2_its", 3.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("max_rchi2_its", max_rchi2_its);
    }

  virtual ~TrackCutAttrChi2ITS() {}
};


/// \class TrackCutAttrChi2TPC
/// \brief Cut on reduced chi2 of TPC
///
struct TrackCutAttrChi2TPC {

  double max_rchi2_tpc;

  bool Pass(const AliFemtoTrack &track)
    {
      return (track.TPCncls() > 0)
              && (track.TPCchi2() / track.TPCncls()) < max_rchi2_tpc;
    }

  TrackCutAttrChi2TPC()
    : max_rchi2_tpc(3.0)
    {}

  TrackCutAttrChi2TPC(AliFemtoConfigObject &cfg)
    : max_rchi2_tpc(cfg.pop_float("max_rchi2_tpc", 3.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("max_rchi2_tpc", max_rchi2_tpc);
    }

  virtual ~TrackCutAttrChi2TPC() {}
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrRemoveNonTpc {

  bool remove_non_tpc;

  bool Pass(const AliFemtoTrack &track)
    {
      return !remove_non_tpc || track.Label() < 0;
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("remove_non_tpc", remove_non_tpc);
    }

  TrackCutAttrRemoveNonTpc()
    : remove_non_tpc(false)
    {}

  TrackCutAttrRemoveNonTpc(AliFemtoConfigObject &cfg)
    : remove_non_tpc(cfg.pop_bool("remove_non_tpc", false))
    {}

  virtual ~TrackCutAttrRemoveNonTpc() {}
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrCharge {

  static const int DEFAULT;

  int charge;

  bool Pass(const AliFemtoTrack &track)
    {
      return charge != 0 && track.Charge() == charge;
    }

  TrackCutAttrCharge()
    : charge(DEFAULT)
    {}

  TrackCutAttrCharge(AliFemtoConfigObject &cfg)
    : charge(cfg.pop_int("charge", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    { cfg.insert("charge", charge); }

  virtual ~TrackCutAttrCharge() {}
};

/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrMomentum {
  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> momentum_range;

  bool Pass(const AliFemtoTrack &track)
    {
      double p = track.P().Mag();
      return momentum_range.first <= p && p < momentum_range.second;
    }

  TrackCutAttrMomentum()
    : momentum_range(DEFAULT)
    {}

  TrackCutAttrMomentum(AliFemtoConfigObject &cfg)
    : momentum_range(cfg.pop_range("momentum_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    { cfg.insert("momentum_range", momentum_range); }

  virtual ~TrackCutAttrMomentum() {}
};

/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrPt {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const double pt = track.Pt();
      return pt_range.first <= pt && pt < pt_range.second;
    }

  TrackCutAttrPt()
    : pt_range(DEFAULT)
    {}

  TrackCutAttrPt(AliFemtoConfigObject &cfg)
    : pt_range(cfg.pop_range("pt_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    { cfg.insert("pt_range", pt_range); }

  virtual ~TrackCutAttrPt() {}
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrEta {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> eta_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const double eta = track.P().PseudoRapidity();
      return eta_range.first <= eta && eta < eta_range.second;
    }

  TrackCutAttrEta()
    : eta_range(DEFAULT)
    {}

  TrackCutAttrEta(AliFemtoConfigObject &cfg)
    : eta_range(cfg.pop_range("eta_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    { cfg.insert("eta_range", eta_range); }

  virtual ~TrackCutAttrEta() {}
};


/// Rapidity assuming mass
struct TrackCutAttrRapidity {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> rapidity_range;
  double rapidity_mass;

  bool Pass(const AliFemtoTrack &track)
    {
      const auto p = track.P();
      const double
        mass_sqrd = rapidity_mass * rapidity_mass,
        e = ::sqrt(track.P().Mag2() + mass_sqrd),
        rap = p.z() == 0
            ? 0.0
            : 0.5 * std::log((e + p.z()) / (e - p.z()));

      return rapidity_range.first <= rap && rap < rapidity_range.second;
    }

  TrackCutAttrRapidity()
    : rapidity_range(DEFAULT)
    , rapidity_mass(0.0)
    {}

  TrackCutAttrRapidity(double mass)
    : rapidity_range(DEFAULT)
    , rapidity_mass(mass)
    {}

  TrackCutAttrRapidity(AliFemtoConfigObject &cfg)
    : rapidity_range(cfg.pop_range("rapidity_range", DEFAULT))
    , rapidity_mass(cfg.pop_float("rapidity_mass", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("rapidity_range", rapidity_range);
      cfg.insert("rapidity_mass", rapidity_mass);
    }

  virtual ~TrackCutAttrRapidity() {}

protected:
};


/// Join
typedef AddTrackCutAttrs<TrackCutAttrEta, TrackCutAttrRapidity> TrackCutAttrEtaY;


/// Cut on `PidProbPion()` within a range of sigma
struct TrackCutAttrMostProbablePion {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> pidprob_pion_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const float prob = track.PidProbPion();
      return pidprob_pion_range.first <= prob && prob < pidprob_pion_range.second;
    }

  TrackCutAttrMostProbablePion()
    : pidprob_pion_range(-1.0, 2.0)
    {}

  TrackCutAttrMostProbablePion(AliFemtoConfigObject &cfg)
    : pidprob_pion_range(cfg.pop_range("pidprob_pion_range", DEFAULT))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    { cfg.insert("pidprob_pion_range", pidprob_pion_range); }

  virtual ~TrackCutAttrMostProbablePion() {}
};


/// Cut using combination of NSigmaTPCPI && NSigmaTOFPI
///
struct TrackCutAttrSigmaPion {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> nsigma_pion_range;
  bool usetpctof;
  double nsigma;

  bool Pass(const AliFemtoTrack &track)
    {
      double momentum = track.P().Mag(),
             nsgima_tpc = track.NSigmaTPCPi(),
             nsgima_tof = track.NSigmaTOFPi();

      return is_pion_nsigma(momentum, nsgima_tpc, nsgima_tof);
    }

  bool is_pion_nsigma(float mom, float sigtpc, float sigtof)
    {
      if (usetpctof) {
        return mom > 0.5 ? sigtof*sigtof + sigtpc*sigtpc < nsigma * nsigma
                         : TMath::Abs(sigtpc) < nsigma;
      }

      if (mom < 0.65) {
        if (sigtof < -99.0) {
          return (mom < 0.35) ? TMath::Abs(sigtpc) < 3.0
               : (mom < 0.50) ? TMath::Abs(sigtpc) < 3.0
                              : TMath::Abs(sigtpc) < 2.0;
        } else {
          return TMath::Abs(sigtof) < 3.0 && TMath::Abs(sigtpc) < 3.0;
        }
      }
      if (mom < 1.5) {
        return TMath::Abs(sigtof) < 3.0 && TMath::Abs(sigtpc) < 5.0;
      }

      return TMath::Abs(sigtof) < 2.0 && TMath::Abs(sigtpc) < 5.0;
    }

  TrackCutAttrSigmaPion()
    : nsigma_pion_range(DEFAULT)
    , usetpctof(true)
    , nsigma(3.0)
    {}

  TrackCutAttrSigmaPion(AliFemtoConfigObject &cfg)
    : nsigma_pion_range(cfg.pop_range("nsigma_pion_range", DEFAULT))
    , usetpctof(cfg.pop_bool("use_tpctof", true))
    , nsigma(cfg.pop_float("pion_nsigma", 3.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("nsigma_pion_range", nsigma_pion_range);
      cfg.insert("use_tpctof", usetpctof);
      cfg.insert("pion_nsigma", nsigma);
    }

  virtual ~TrackCutAttrSigmaPion() {}
};


/// Cut on the track's nsigma to proton
///
struct TrackCutAttrSigmaProton {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> nsigma_proton_range;
  bool usetpctof;
  double nsigma;

  bool Pass(const AliFemtoTrack &track)
    {
      const double
        momentum = track.P().Mag(),
        nsgima_tpc = track.NSigmaTPCP(),
        nsgima_tof = track.NSigmaTOFP();

      return is_proton_nsigma(momentum, nsgima_tpc, nsgima_tof);
    }

  bool is_proton_nsigma(float mom, float sigtpc, float sigtof)
    {
      if (usetpctof) {
        return mom > 0.5 ? sigtof*sigtof + sigtpc*sigtpc < nsigma*nsigma
                         : TMath::Abs(sigtpc) < nsigma;
      }

      if (mom <= 0.8) {
        return TMath::Abs(sigtpc) < 3.0;
      }
      if (mom < 2.5) {
        return TMath::Abs(sigtof) < 3.0 && TMath::Abs(sigtpc) < 3.0;
      }

      return TMath::Abs(sigtof) < 2.0 && TMath::Abs(sigtpc) < 3.0;
    }

  TrackCutAttrSigmaProton()
    : nsigma_proton_range(-3.0, 3.0)
    , usetpctof(true)
    , nsigma(3)
    {}

  TrackCutAttrSigmaProton(AliFemtoConfigObject &cfg)
    : nsigma_proton_range(cfg.pop_range("nsigma_proton_range", DEFAULT))
    , usetpctof(cfg.pop_bool("proton_usetpctof", true))
    , nsigma(cfg.pop_bool("proton_nsigma", true))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("nsigma_proton_range", nsigma_proton_range);
      cfg.insert("proton_usetpctof", usetpctof);
      cfg.insert("proton_nsigma", nsigma);
    }

  virtual ~TrackCutAttrSigmaProton() {}
};


/// Cut on ITS layer points
struct TrackCutAttrItsCluster {

  AliESDtrackCuts::ITSClusterRequirement cluster_reqs[3];

  bool Pass(const AliFemtoTrack &track)
    {
      for (int i=0; i < 3; ++i) {
        if (!check_ITS_cluster(cluster_reqs[i],
                               track.HasPointOnITSLayer(i * 2),
                               track.HasPointOnITSLayer(i * 2 + 1))) {
          return false;
        }
      }
      return true;
    }

  bool check_ITS_cluster(AliESDtrackCuts::ITSClusterRequirement req,
                         Bool_t clusterL1,
                         Bool_t clusterL2)
    {
      // checks if the cluster requirement is fullfilled (in this case: return kTRUE)

      switch (req) {
        case AliESDtrackCuts::kOff:        return kTRUE;
        case AliESDtrackCuts::kNone:       return !clusterL1 && !clusterL2;
        case AliESDtrackCuts::kAny:        return clusterL1 || clusterL2;
        case AliESDtrackCuts::kFirst:      return clusterL1;
        case AliESDtrackCuts::kOnlyFirst:  return clusterL1 && !clusterL2;
        case AliESDtrackCuts::kSecond:     return clusterL2;
        case AliESDtrackCuts::kOnlySecond: return clusterL2 && !clusterL1;
        case AliESDtrackCuts::kBoth:       return clusterL1 && clusterL2;
      }

      return false;
    }

  TrackCutAttrItsCluster()
    {
      std::fill_n(cluster_reqs, 3, AliESDtrackCuts::kOff);
    }

  TrackCutAttrItsCluster(AliFemtoConfigObject &cfg)
    {
      AliFemtoConfigObject::ArrayValue_t v;
      if (cfg.pop_and_load("cluster_reqs", v)) {
        for (int i=0; i<3; ++i) {
          cluster_reqs[i] = static_cast<AliESDtrackCuts::ITSClusterRequirement>(v[i].as_int());
        }
      } else {
        std::fill_n(cluster_reqs, 3, AliESDtrackCuts::kOff);
      }
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      std::vector<AliFemtoConfigObject> c(cluster_reqs, cluster_reqs+3);
      cfg.insert("cluster_reqs", c);
    }

  virtual ~TrackCutAttrItsCluster() {}
};


/// Rapidity assuming mass
struct TrackCutAttrElectronRejection {

  bool reject_electrons;

  bool Pass(const AliFemtoTrack &track)
    {
      return !reject_electrons
          || !is_electron(track.NSigmaTPCE(),
                          track.NSigmaTPCPi(),
                          track.NSigmaTPCK(),
                          track.NSigmaTPCP());
    }

  static bool is_electron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
    {
      return TMath::Abs(nsigmaTPCE) < 3.0
          && TMath::Abs(nsigmaTPCPi) > 3.0
          && TMath::Abs(nsigmaTPCK) > 3.0
          && TMath::Abs(nsigmaTPCP) > 3.0;
    }

  TrackCutAttrElectronRejection()
    : reject_electrons(true)
    {}

  TrackCutAttrElectronRejection(AliFemtoConfigObject &cfg)
    : reject_electrons(cfg.pop_bool("reject_electrons", true))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("reject_electrons", reject_electrons);
    }

  virtual ~TrackCutAttrElectronRejection() {}
};


} // namespace pwgfemto



#include <AliFemtoTrackCut.h>


/// \class AliFemtoTrackCutAttr
/// \brief Bridge from AliFemtoTrackCut to a metaclass of TrackCut-Attrs
///
template <typename CRTP, typename CutAttrType>
class AliFemtoTrackCutAttr : public AliFemtoTrackCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  AliFemtoTrackCutAttr()
    : AliFemtoTrackCut()
    , CutAttrs()
    {}

  AliFemtoTrackCutAttr(AliFemtoConfigObject &cfg)
    : AliFemtoTrackCut()
    , CutAttrType(cfg)
    {}

  virtual bool Pass(AliFemtoTrack *ev)
    {
      return CutAttrs::Pass(*ev);
    }

  virtual AliFemtoString Report()
    {
      return "AliFemtoTrackCutAttr Report\n";
    }

  virtual TList* ListSettings() const
    {
      TList* list = new TList();
      AppendSettings(*list);
      return list;
    }

  AliFemtoConfigObject GetConfiguration() const
    {
      AliFemtoConfigObject cfg = AliFemtoConfigObject::BuildMap()
                                    ("_class", CRTP::ClassName());
      CutAttrs::FillConfiguration(cfg);
      return cfg;
    }

  virtual void AppendSettings(TCollection &, TString prefix="") const = 0;
  virtual ~AliFemtoTrackCutAttr() = 0;
};



#endif
