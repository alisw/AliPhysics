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
      return status == 0 || (trk.Flags() & status) == status;
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
      if (status == 0) {
        return;
      }
      cfg.insert("status", status);
    }

  virtual ~TrackCutAttrStatus() {}
};

struct TrackCutAttrRemoveKinks {

  bool remove_kinks;

  bool Pass(const AliFemtoTrack &track) const
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

  bool Pass(const AliFemtoTrack &track) const
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

  bool Pass(const AliFemtoTrack &track) const
    {
      return std::abs(track.ImpactZ()) <= max_z
          && track.ImpactD() < max_xy;
    }

  TrackCutAttrImpact()
    // : xy_range(0, 0.25)
    : max_xy(0.25)
    , max_z(0.3)
    {}

  TrackCutAttrImpact(AliFemtoConfigObject &cfg)
    : max_xy(cfg.pop_num("max_xy", 0.25))
    , max_z(cfg.pop_num("max_z", 0.3))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("max_z", max_z);
      cfg.insert("max_xy", max_xy);
    }

  virtual ~TrackCutAttrImpact() {}
};


/// Remove tracks with number of TPC clusters below some threshold
struct TrackCutAttrMinNclsTPC {
  int ncls_tpc_min;

  bool Pass(const AliFemtoTrack &track) const
    {
      return ncls_tpc_min <= track.TPCncls();
    }

  TrackCutAttrMinNclsTPC()
    : ncls_tpc_min(70)
    {}

  TrackCutAttrMinNclsTPC(AliFemtoConfigObject &cfg)
    : ncls_tpc_min(cfg.pop_int("ncls_tpc_min", 70))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("ncls_tpc_min", ncls_tpc_min);
    }

  virtual ~TrackCutAttrMinNclsTPC() {}
};

/// Remove tracks with number of ITS clusters below some threshold
struct TrackCutAttrMinNclsITS {
  int ncls_its_min;

  bool Pass(const AliFemtoTrack &track) const
    {
      return ncls_its_min <= track.ITSncls();
    }

  TrackCutAttrMinNclsITS()
    : ncls_its_min(0)
    {}

  TrackCutAttrMinNclsITS(AliFemtoConfigObject &cfg)
    : ncls_its_min(cfg.pop_int("ncls_its_min", 0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("ncls_its_min", ncls_its_min);
    }

  virtual ~TrackCutAttrMinNclsITS() {}
};


/// \class TrackCutAttrChi2ITS
/// \brief Cut on reduced chi2 of results in ITS
///
/// If rchi2_its_max is negative, this cut is ignored
///
struct TrackCutAttrChi2ITS {

  double rchi2_its_max;
  double rchi2_its_min;

  bool Pass(const AliFemtoTrack &track) const
    {
      if (rchi2_its_max < 0) {
        return true;
      }
      double chi2 = track.ITSchi2perNDF();
      return rchi2_its_min <= chi2 && chi2 < rchi2_its_max;
    }

  TrackCutAttrChi2ITS()
    : rchi2_its_max(-1.0)
    , rchi2_its_min(0.0)
    {}

  TrackCutAttrChi2ITS(AliFemtoConfigObject &cfg)
    : rchi2_its_max(cfg.pop_num("rchi2_its_max", -1.0))
    , rchi2_its_min(cfg.pop_num("rchi2_its_min", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (rchi2_its_max < 0) {
        return;
      }
      cfg.insert("rchi2_its_max", rchi2_its_max);
      if (rchi2_its_min != 0.0) {
        cfg.insert("rchi2_its_min", rchi2_its_min);
      }
    }

  virtual ~TrackCutAttrChi2ITS() {}
};


/// \class TrackCutAttrChi2TPC
/// \brief Track quality cut on reduced chi2 of TPC fit
///
/// This cut is diabled if rchi2_tpc_max is negative.
///
struct TrackCutAttrChi2TPC {

  double rchi2_tpc_max;
  double rchi2_tpc_min;

  bool Pass(const AliFemtoTrack &track) const
    {
      if (rchi2_tpc_max < 0) {
        return true;
      }
      double chi2 = track.TPCchi2perNDF();
      return rchi2_tpc_min <= chi2 && chi2 < rchi2_tpc_max;
    }

  TrackCutAttrChi2TPC()
    : rchi2_tpc_max(-1.0)
    , rchi2_tpc_min(0.0)
    {}

  TrackCutAttrChi2TPC(AliFemtoConfigObject &cfg)
    : rchi2_tpc_max(cfg.pop_num("rchi2_tpc_max", -1.0))
    , rchi2_tpc_min(cfg.pop_num("rchi2_tpc_min", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (rchi2_tpc_max < 0) {
        return;
      }
      cfg.insert("rchi2_tpc_max", rchi2_tpc_max);
      if (rchi2_tpc_min != 0.0) {
        cfg.insert("rchi2_tpc_min", rchi2_tpc_min);
      }
    }

  virtual ~TrackCutAttrChi2TPC() {}
};


/// Remove tracks with negative label
///
///
struct TrackCutAttrRemoveNegLabel {

  bool remove_neg_label;

  bool Pass(const AliFemtoTrack &track) const
    {
      return !remove_neg_label || track.Label() >= 0;
    }

  TrackCutAttrRemoveNegLabel()
    : remove_neg_label(false)
    {}

  TrackCutAttrRemoveNegLabel(AliFemtoConfigObject &cfg)
    : remove_neg_label(cfg.pop_bool("remove_neg_label", false))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (!remove_neg_label) {
        return;
      }
      cfg.insert("remove_neg_label", remove_neg_label);
    }

  virtual ~TrackCutAttrRemoveNegLabel() {}
};


/// Cut based on charge
///
/// If charge is 0, the cut is not applied.
///
struct TrackCutAttrCharge {

  static const int DEFAULT;

  int charge;

  bool Pass(const AliFemtoTrack &track) const
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

/// Remove outside momentum range
struct TrackCutAttrMomentum {
  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> momentum_range;

  bool Pass(const AliFemtoTrack &track) const
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

/// Remove tracks outside a range of transverse momentum
struct TrackCutAttrPt {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoTrack &track) const
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


/// Remove tracks outside a range of pseudorapidity
struct TrackCutAttrEta {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> eta_range;

  bool Pass(const AliFemtoTrack &track) const
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


/// Cut tracks on calculated rapidity with assumed mass
struct TrackCutAttrRapidity {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> rapidity_range;
  double rapidity_mass;

  bool Pass(const AliFemtoTrack &track) const
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
    , rapidity_mass(cfg.pop_num("rapidity_mass", 0.0))
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

  bool Pass(const AliFemtoTrack &track) const
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

/// Equivalent to AliFemtoESDTrackCut with PIDMethodType 'kContour'
/// and `SetMostProbablePion();`
///
struct TrackCutAttrPidContourPion {

  bool tpc_contour_pion;

  bool Pass(const AliFemtoTrack &track) const
    {
      const float
        p = track.P().Mag(),
        tpc_signal = track.TPCsignal();
      return tpc_contour_pion && is_pion_tpc_dedx(p, tpc_signal);
    }

  static bool is_pion_tpc_dedx(float mom, float dedx)
    {
      const double
        a1 = -343.75,
        b1 = 168.125,
        a2 = 0.0,
        b2 = 65.0;

      return dedx < ((mom < 0.32) ? a1 * mom + b1 : a2 * mom + b2);
    }

  TrackCutAttrPidContourPion()
    : tpc_contour_pion(true)
    {}

  TrackCutAttrPidContourPion(AliFemtoConfigObject &cfg)
    : tpc_contour_pion(cfg.pop_bool("tpc_contour_pion", true))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
     cfg.insert("tpc_contour_pion", tpc_contour_pion);
    }

  virtual ~TrackCutAttrPidContourPion() {}
};


/// Cut on the track's TOF pion-time
///
struct TrackCutAttrRejectKaonTofPionTime {

  double tof_kaon_reject_sigma;
  double tof_sigma_pion;
  double tof_momentum_limit;


  bool Pass(const AliFemtoTrack &track) const
    {
      const double
        p = track.P().Mag(),
        tof_signal = track.TOFpionTime();

      if (p < tof_momentum_limit) {
        return true;
      }

      if (!std::isnan(tof_sigma_pion) and tof_sigma_pion <= std::abs(track.NSigmaTOFPi())) {
        return false;
      }

      const double
        decay = 1.50388857 + tof_kaon_reject_sigma * (0.14601506 + tof_kaon_reject_sigma * 0.04899996),
        amplitude = 6688.09178 + tof_kaon_reject_sigma * (5.68606498 + tof_kaon_reject_sigma * 524.181339),

        max_kaon_signal = amplitude * std::exp(-decay * p);

      return tof_signal < max_kaon_signal;
      // return tof_signal < 9150.442842004599 * std::exp(-2.089165525857385 * p);
    }

  TrackCutAttrRejectKaonTofPionTime()
    : tof_kaon_reject_sigma(3.5)
    , tof_sigma_pion(NAN)
    , tof_momentum_limit(0.5)
    {}

  TrackCutAttrRejectKaonTofPionTime(AliFemtoConfigObject &cfg)
    : tof_kaon_reject_sigma(cfg.pop_num("tof_kaon_reject_sigma", 3.5))
    , tof_sigma_pion(cfg.pop_num("tof_sigma_pion", NAN))
    , tof_momentum_limit(cfg.pop_num("tof_momentum_limit", 0.5))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("tof_kaon_reject_sigma", tof_kaon_reject_sigma);
      if (!std::isnan(tof_sigma_pion)) {
        cfg.insert("tof_sigma_pion", tof_sigma_pion);
      }
      cfg.insert("tof_momentum_limit", tof_momentum_limit);
    }

  virtual ~TrackCutAttrRejectKaonTofPionTime()
    { }
};

/// Cut away Kaon TOF sigma
struct TrackCutAttrRejectKaonTofSigma {

  double tof_kaon_reject_sigma;
  double tof_kaon_momentum_limit;

  bool Pass(const AliFemtoTrack &track) const
    {
      if (std::isnan(tof_kaon_reject_sigma)) {
        return true;
      }

      if (track.P().Mag() < tof_kaon_momentum_limit) {
        return true;
      }

      return tof_kaon_reject_sigma < std::abs(track.NSigmaTOFK());
    }

  TrackCutAttrRejectKaonTofSigma()
    : tof_kaon_reject_sigma(2.0)
    , tof_kaon_momentum_limit(0.60)
    {
    }

  TrackCutAttrRejectKaonTofSigma(AliFemtoConfigObject &cfg)
    : tof_kaon_reject_sigma(cfg.pop_num("tof_kaon_reject_sigma", 2.0))
    , tof_kaon_momentum_limit(cfg.pop_num("tof_kaon_momentum_limit", 0.60))
    {
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (std::isnan(tof_kaon_reject_sigma)) {
        return;
      }

      cfg.insert("tof_kaon_reject_sigma", tof_kaon_reject_sigma);
      cfg.insert("tof_kaon_momentum_limit", tof_kaon_momentum_limit);
    }

  virtual ~TrackCutAttrRejectKaonTofSigma()
    { }
};


/// Cut away Kaon TOF sigma
struct TrackCutAttrRejectProtonTofSigma {

  double tof_proton_reject_sigma;
  double tof_proton_momentum_limit;

  bool Pass(const AliFemtoTrack &track) const
    {
      if (track.P().Mag() < tof_proton_momentum_limit) {
        return true;
      }

      return tof_proton_reject_sigma < std::abs(track.NSigmaTOFP());
    }

  TrackCutAttrRejectProtonTofSigma()
    : tof_proton_reject_sigma(2.0)
    , tof_proton_momentum_limit(0.75)
    {
    }

  TrackCutAttrRejectProtonTofSigma(AliFemtoConfigObject &cfg)
    : tof_proton_reject_sigma(cfg.pop_num("tof_proton_reject_sigma", 2.0))
    , tof_proton_momentum_limit(cfg.pop_num("tof_proton_momentum_limit", 0.75))
    {
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("tof_proton_reject_sigma", tof_proton_reject_sigma);
      cfg.insert("tof_proton_momentum_limit", tof_proton_momentum_limit);
    }

  virtual ~TrackCutAttrRejectProtonTofSigma()
    { }
};


///
struct TrackCutAttrTpcSigmaPion {

  double tpc_sigma_pion;

  bool Pass(const AliFemtoTrack &track) const
    {
      return std::abs(track.NSigmaTPCPi()) < tpc_sigma_pion;
    }

  TrackCutAttrTpcSigmaPion()
    : tpc_sigma_pion(2.5)
    { }

  TrackCutAttrTpcSigmaPion(AliFemtoConfigObject &cfg)
    : tpc_sigma_pion(cfg.pop_num("tpc_sigma_pion", 2.5))
    { }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("tpc_sigma_pion", tpc_sigma_pion);
    }

  virtual ~TrackCutAttrTpcSigmaPion()
    { }
};

///
struct TrackCutAttrTofSigmaPion {

  double tof_sigma_pion;
  double tof_pion_momentum_limit;
  bool require_tof;

  bool Pass(const AliFemtoTrack &track) const
    {
      // out of TOF range
      if (track.P().Mag2() < tof_pion_momentum_limit * tof_pion_momentum_limit) {
        return true;
      }

      // Reject TOF "no signal" value
      if (require_tof && track.TOFpionTime() == 100000.0f) {
        return false;
      }

      if (std::isnan(tof_sigma_pion)) {
        return true;
      }

      return std::abs(track.NSigmaTOFPi()) < tof_sigma_pion;
    }

  TrackCutAttrTofSigmaPion()
    : tof_sigma_pion(5.0)
    , tof_pion_momentum_limit(0.5)
    , require_tof(true)
    { }

  TrackCutAttrTofSigmaPion(AliFemtoConfigObject &cfg)
    : tof_sigma_pion(cfg.pop_num("tof_sigma_pion", 5.0))
    , tof_pion_momentum_limit(cfg.pop_num("tof_pion_momentum_limit", 0.5))
    , require_tof(cfg.pop_num("tof_required", true))
    { }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("tof_required", require_tof);
      cfg.insert("tof_pion_momentum_limit", tof_pion_momentum_limit);

      if (!std::isnan(tof_sigma_pion)) {
        cfg.insert("tof_sigma_pion", tof_sigma_pion);
      }
    }

  virtual ~TrackCutAttrTofSigmaPion()
    { }
};


/// Cut using combination of NSigmaTPCPI && NSigmaTOFPI
///
struct TrackCutAttrSigmaPion {

  static const std::pair<double, double> DEFAULT;

  std::pair<double, double> nsigma_pion_range;
  bool usetpctof;
  double nsigma_pion;
  double tpctof_limit;

  bool Pass(const AliFemtoTrack &track) const
    {
      double momentum = track.P().Mag(),
             nsgima_tpc = track.NSigmaTPCPi(),
             nsgima_tof = track.NSigmaTOFPi();

      bool pass = is_pion_nsigma(momentum, nsgima_tpc, nsgima_tof);
      return pass;
    }

  bool is_pion_nsigma(float mom, float sigtpc, float sigtof) const
    {
      if (usetpctof) {
        return mom > tpctof_limit ? sigtof*sigtof + sigtpc*sigtpc < nsigma_pion * nsigma_pion
                                  : TMath::Abs(sigtpc) < nsigma_pion;
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
    , nsigma_pion(3.0)
    , tpctof_limit(0.5)
    {}

  TrackCutAttrSigmaPion(AliFemtoConfigObject &cfg)
    : nsigma_pion_range(cfg.pop_range("nsigma_pion_range", DEFAULT))
    , usetpctof(cfg.pop_bool("use_tpctof", true))
    , nsigma_pion(cfg.pop_num("nsigma_pion", 3.0))
    , tpctof_limit(cfg.pop_num("tpctof_limit", 0.5))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      // cfg.insert("use_tpctof", usetpctof);
      cfg.insert("nsigma_pion", nsigma_pion);
      cfg.insert("tpctof_limit", tpctof_limit);
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

  bool Pass(const AliFemtoTrack &track) const
    {
      const double
        momentum = track.P().Mag(),
        nsgima_tpc = track.NSigmaTPCP(),
        nsgima_tof = track.NSigmaTOFP();

      return is_proton_nsigma(momentum, nsgima_tpc, nsgima_tof);
    }

  bool is_proton_nsigma(float mom, float sigtpc, float sigtof) const
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

  bool Pass(const AliFemtoTrack &track) const
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

  static bool check_ITS_cluster(AliESDtrackCuts::ITSClusterRequirement req,
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


/// Reject tracks close to the electron band
///
struct TrackCutAttrRejectTpcElectron {

  double electron_tpc_sigma_min;

  bool Pass(const AliFemtoTrack &track) const
    {
      return std::abs(track.NSigmaTPCE()) > electron_tpc_sigma_min;
    }

  TrackCutAttrRejectTpcElectron()
    : electron_tpc_sigma_min(3.0)
    { }

  TrackCutAttrRejectTpcElectron(AliFemtoConfigObject &cfg)
    : electron_tpc_sigma_min(cfg.pop_num("electron_tpc_sigma_min", 3.0))
    { }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("electron_tpc_sigma_min", electron_tpc_sigma_min);
    }

  virtual ~TrackCutAttrRejectTpcElectron()
    { }
};

/// Reject tracks within 3 TPC-sigma of electrons and outside 3-sigma
/// of other particles
struct TrackCutAttrElectronRejection {

  bool reject_electrons;

  bool Pass(const AliFemtoTrack &track) const
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
                                  ("_class", static_cast<const CRTP*>(this)->ClassName());
      CutAttrs::FillConfiguration(cfg);
      return cfg;
    }

  virtual void AppendSettings(TCollection &, TString prefix="") const = 0;
  virtual ~AliFemtoTrackCutAttr()
    {}
};


#endif
