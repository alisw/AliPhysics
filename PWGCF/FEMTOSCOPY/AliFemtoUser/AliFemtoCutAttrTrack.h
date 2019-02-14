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

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddTrackCutAttrs() = 0;
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

  virtual ~TrackSelectionCut() = 0;
};


struct TrackCutAttrStatus {
  ULong64_t status;

  TrackCutAttrStatus()
    : status(0)
    {}

  TrackCutAttrStatus(ULong64_t st)
    : status(st)
    {}

  bool Pass(const AliFemtoTrack &track)
    {
      return (track.Flags() & status) == status;
    }

  virtual ~TrackCutAttrStatus() = 0;
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

  virtual ~TrackCutAttrRemoveKinks() = 0;
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrRemoveFakeITS {

  bool remove_its_fake;

  bool Pass(const AliFemtoTrack &track)
    {
      return !remove_its_fake || track.ITSncls() < 0;
    }

  virtual ~TrackCutAttrRemoveFakeITS() = 0;
};


struct TrackCutAttrImpact {

  std::pair<double, double> xy_range;
  double min_z;

  bool Pass(const AliFemtoTrack &track)
    {
      const double
        impact_d = track.ImpactD(),
        impact_z = track.ImpactZ();

      return min_z <= impact_z
             && (xy_range.first <= impact_d) && (impact_d < xy_range.second);
    }

  virtual ~TrackCutAttrImpact() = 0;
};


struct TrackCutAttrChi2ITS {

  double max_rchi2_its;

  bool Pass(const AliFemtoTrack &track)
    {
      return (track.ITSncls() > 0)
              && (track.ITSchi2() / track.ITSncls()) > max_rchi2_its;
    }

  virtual ~TrackCutAttrChi2ITS() = 0;
};

struct TrackCutAttrChi2TPC {

  double max_rchi2_tpc;

  bool Pass(const AliFemtoTrack &track)
    {
      return (track.TPCncls() > 0)
              && (track.TPCchi2() / track.TPCncls()) > max_rchi2_tpc;
    }

  virtual ~TrackCutAttrChi2TPC() = 0;
};

/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrRemoveNonTpc {

  bool remove_non_tpc;

  bool Pass(const AliFemtoTrack &track)
    {
      return !remove_non_tpc || track.Label() < 0;
    }

  virtual ~TrackCutAttrRemoveNonTpc() = 0;
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrCharge {

  int charge;

  bool Pass(const AliFemtoTrack &track)
    {
      return charge != 0 && track.Charge() == charge;
    }

  virtual ~TrackCutAttrCharge() = 0;
};

/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrMomentum {

  int charge;

  bool Pass(const AliFemtoTrack &track)
    {
      return charge != 0 && track.Charge() == charge;
    }
};

/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrPt {

  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const double pt = track.Pt();
      return pt_range.first <= pt && pt < pt_range.second;
    }

  virtual ~TrackCutAttrPt() = 0;
};


/// Remove tracks with "negative" number of ITS clusters
struct TrackCutAttrEta {

  std::pair<double, double> eta_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const double eta = track.P().PseudoRapidity();
      return eta_range.first <= eta && eta < eta_range.second;
    }

  virtual ~TrackCutAttrEta() = 0;
};


/// Rapidity assuming mass
struct TrackCutAttrRapdidity {

  std::pair<double, double> rapidity_range;

  bool Pass(const AliFemtoTrack &track)
    {
      const auto p = track.P();
      const double
        e = ::sqrt(track.P().Mag2() + fMass * fMass),
        rap = p.z() == 0
            ? 0.0
            : 0.5 * std::log((e + p.z()) / (e - p.z()));

      return rapidity_range.first <= rap && rap < rapidity_range.second;
    }

  virtual ~TrackCutAttrRapdidity() = 0;

protected:
  double fMass;
};


/// Join
typedef AddTrackCutAttrs<TrackCutAttrEta, TrackCutAttrRapdidity> TrackCutAttrEtaY;


/// Cut on `PidProbPion()` within a range of sigma
struct TrackCutAttrMostProbablePion {

  std::pair<double, double> pidprob_pion_range;

  TrackCutAttrMostProbablePion()
    : pidprob_pion_range(-1.0, 2.0)
    {}

  bool Pass(const AliFemtoTrack &track)
    {
      const float prob = track.PidProbPion();
      return pidprob_pion_range.first <= prob && prob < pidprob_pion_range.second;
    }

  virtual ~TrackCutAttrMostProbablePion() = 0;
};


/// Cut using combination of NSigmaTPCPI && NSigmaTOFPI
///
struct TrackCutAttrSigmaPion {

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

  virtual ~TrackCutAttrSigmaPion() = 0;
};


struct TrackCutAttrItsCluster {

  AliESDtrackCuts::ITSClusterRequirement cluster_reqs[3];

  TrackCutAttrItsCluster()
    {
      std::fill_n(cluster_reqs, 3, AliESDtrackCuts::kOff);
    }

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


  virtual ~TrackCutAttrItsCluster() = 0;
};


struct TrackCutAttrSigmaProton {

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

  virtual ~TrackCutAttrSigmaProton() = 0;
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


  virtual ~TrackCutAttrElectronRejection() = 0;
};


} // namespace pwgfemto


#endif
