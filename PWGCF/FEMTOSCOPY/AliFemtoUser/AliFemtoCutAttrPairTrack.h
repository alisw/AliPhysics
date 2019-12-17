///
/// \file AliFemtoUser/AliFemtoPairCutTrackAttr.h
///


#pragma once

#ifndef ALIFEMTOPAIRCUTTRACKATTR_H
#define ALIFEMTOPAIRCUTTRACKATTR_H


#include "AliFemtoConfigObject.h"
#include "AliFemtoPair.h"
#include "AliFemtoPairCut.h"

// re-use static Δφ* calculation
#include "AliFemtoPairCutDetaDphi.h"

namespace pwgfemto {

#if __cplusplus < 201103L

// Old C++ not supported
template <typename T>
bool call_pair_passes(T &, const AliFemtoPair &)
{
  return false;
}

#elif true

/// called if the cut expects two tracks
template <typename T>
bool
call_pair_passes(
  typename std::enable_if<
    std::is_same<
      decltype(std::declval<T>().Pass(std::declval<const AliFemtoTrack &>(), std::declval<const AliFemtoTrack &>())),
      bool
      >::value,
    T>::type &obj,
  const AliFemtoPair &pair)
{
  const AliFemtoTrack &track1 = *pair.Track1()->Track(),
                      &track2 = *pair.Track2()->Track();
  return obj.Pass(track1, track2);
}

/// called if the cut expects a pair
template <typename T>
bool
call_pair_passes(
  typename std::enable_if<
    std::is_same<
      decltype(std::declval<T>().Pass(std::declval<const AliFemtoPair &>())),
      bool
      >::value,
    T>::type &obj,
  const AliFemtoPair &pair)
{
  return obj.Pass(pair);
}

#endif

/// \class AddPairCutAttrs
/// \brief Join cut-types together
///
template <typename T1, typename T2>
struct AddPairCutAttrs : public T1, public T2 {

  bool Pass(const AliFemtoPair &pair)
    {
      return call_pair_passes<T1>(static_cast<T1&>(*this), pair)
          && call_pair_passes<T2>(static_cast<T2&>(*this), pair);
    }

  AddPairCutAttrs()
    : T1()
    , T2()
    { }

  AddPairCutAttrs(AliFemtoConfigObject &cfg)
    : T1(cfg)
    , T2(cfg)
    { }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddPairCutAttrs() {}
};


/// \class PairCutTrackAttrAvgSep
/// \brief Range of average separation
///
struct PairCutTrackAttrAvgSep {

  double avgsep_min;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      return calc_avg_sep(track1, track2) > avgsep_min;
    }

  PairCutTrackAttrAvgSep()
    : avgsep_min(0.0)
    {
    }

  PairCutTrackAttrAvgSep(AliFemtoConfigObject &cfg)
    : avgsep_min(cfg.pop_num("avgsep_min", 0.0))
    {
    }

  static double calc_avg_sep(const AliFemtoTrack &track1,
                             const AliFemtoTrack &track2)
    {
      return AliFemtoPair::CalcAvgSepTracks(track1, track2);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("avgsep_min", avgsep_min);
    }

  virtual ~PairCutTrackAttrAvgSep() {}
};


/// Cut on the sum of the pT
struct PairCutTrackAttrPt {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      const double pt_sum = track1.Pt() + track2.Pt();
      return pt_range.first <= pt_sum && pt_sum < pt_range.second;
    }

  PairCutTrackAttrPt()
    : pt_range(DEFAULT)
    {
    }

  PairCutTrackAttrPt(AliFemtoConfigObject &cfg)
    : pt_range(cfg.pop_range("pt_range", DEFAULT))
    {
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("pt_range", pt_range);
    }

  virtual ~PairCutTrackAttrPt() {}
};


/// Cut on share quality and fraction
struct PairCutTrackAttrShareQuality {

  double share_fraction_max;
  double share_quality_max;

  bool Pass(const AliFemtoPair &pair) const
    {
      // quick-return if we don't want to cut
      if (share_fraction_max >= 1.0 && share_quality_max >= 1.0) {
        return true;
      }

      double share_fraction, share_quality;
      pair.CalcTrackShareQualFractions(share_fraction, share_quality);

      return share_fraction <= share_fraction_max
          && share_quality <= share_quality_max;
    }

  PairCutTrackAttrShareQuality()
    : share_fraction_max(1.0)
    , share_quality_max(1.0)
    {}

  PairCutTrackAttrShareQuality(AliFemtoConfigObject &cfg)
    : share_fraction_max(cfg.pop_num("share_fraction_max", 1.0))
    , share_quality_max(cfg.pop_num("share_quality_max", 1.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("share_quality_max", std::min({share_quality_max, 1.0}));
      cfg.insert("share_fraction_max", std::min({share_fraction_max, 1.0}));
    }

  virtual ~PairCutTrackAttrShareQuality() {}
};


/// \class PairCutTrackAttrDetaDphiStar
/// \brief A pair cut which cuts on the Δη Δφ of the pair.
///
/// Pairs pass which have |Δη| > delta_eta_min
/// and √(Δφ² + Δη²) > delta_phi_min
///
/// The difference in phi is calculated by examining the tracks'
/// azimuthal angle at a particular radius, as determined by the
/// magnetic field of the event.
/// Note: fCurrentMagneticField should be set *before* using this cut
/// It is recommended to do this in the EventBegin method of
/// AliFemtoPairCut.
///
/// The default value for this radius is 1.2 (inside the TPC) but
/// this maybe changed via by changing the phistar_radius member.
///
///    \Delta \phi_{min}* = \phi_1 - \phi_2
///                       + arcsin \left( \frac{ z_1 \cdot B_z \cdot R}{2 p_{T1}} \right)
///                       - arcsin \left( \frac{ z_2 \cdot B_z \cdot R}{2 p_{T2}} \right)
///
///
struct PairCutTrackAttrDetaDphiStar {

  float delta_eta_min,
        delta_phistar_min,
        phistar_radius;

  float fCurrentMagneticField;

  PairCutTrackAttrDetaDphiStar()
   : delta_eta_min(0.0)
   , delta_phistar_min(0.0)
   , phistar_radius(1.1)
   , fCurrentMagneticField(0.0)
   {
   }

  PairCutTrackAttrDetaDphiStar(AliFemtoConfigObject &cut)
   : delta_eta_min(cut.pop_num("delta_eta_min", 0.0))
   , delta_phistar_min(cut.pop_num("delta_phistar_min", 0.0))
   , phistar_radius(cut.pop_num("phistar_radius", 1.1))
   , fCurrentMagneticField(0.0)
   {
   }

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      const AliFemtoThreeVector &p1 = track1.P(),
                                &p2 = track2.P();

      const double deta = calc_delta_eta(p1, p2);
      if (std::abs(deta) >= delta_eta_min) {
        return true;
      }

      const double dphi = AliFemtoPairCutDetaDphi::CalculateDPhiStar(
                            p1, track1.Charge(),
                            p2, track2.Charge(),
                            phistar_radius,
                            fCurrentMagneticField);
      if (std::abs(dphi) >= delta_phistar_min) {
        return true;
      }

      const double
        a = deta / delta_eta_min,
        b = dphi / delta_phistar_min;

      // cut within the ellipse
      return a*a + b*b >= 1.0;
    }

  static double calc_delta_eta(const AliFemtoThreeVector &p1,
                               const AliFemtoThreeVector &p2)
    {
      return p1.PseudoRapidity() - p2.PseudoRapidity();
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.Update(AliFemtoConfigObject::BuildMap()
                 ("delta_eta_min", delta_eta_min)
                 ("delta_phistar_min", delta_phistar_min)
                 ("phistar_radius", phistar_radius));
    }

  virtual ~PairCutTrackAttrDetaDphiStar() {}
};


/// \class PairCutTrackAttrDetaDphi
/// \brief A pair cut which cuts on the Δη Δφ of the pair.
///
struct PairCutTrackAttrDetaDphi {

  double min_delta_eta;
  double min_delta_phi;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      const auto
        &p1 = track1.P(),
        &p2 = track2.P();

      const double
        deta = p1.PseudoRapidity() - p1.PseudoRapidity(),
        dphi = p1.Phi() - p2.Phi();

      return min_delta_eta <= std::abs(deta)
          && min_delta_phi <= std::abs(dphi);
    }

  PairCutTrackAttrDetaDphi()
    : min_delta_eta(0.0)
    , min_delta_phi(0.0)
    {}

  PairCutTrackAttrDetaDphi(AliFemtoConfigObject &cfg)
    : min_delta_eta(cfg.pop_num("min_delta_eta", 0.0))
    , min_delta_phi(cfg.pop_num("min_delta_phi", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.Update(AliFemtoConfigObject::BuildMap()
                 ("min_delta_eta", min_delta_eta)
                 ("min_delta_phi", min_delta_phi));
    }

  virtual ~PairCutTrackAttrDetaDphi() = 0;
};


/// Remove pairs of tracks with the same particle
struct PairCutTrackAttrSameLabel {

  bool remove_same_label;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      // either not removing same label, or passes if they are not equal
      return !remove_same_label || abs(track1.Label()) != abs(track2.Label());
    }

  PairCutTrackAttrSameLabel()
    : remove_same_label(false)
    {}

  PairCutTrackAttrSameLabel(AliFemtoConfigObject &cfg)
    : remove_same_label(cfg.pop_bool("remove_same_label", false))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (!remove_same_label) {
        return;
      }
      cfg.insert("remove_same_label", remove_same_label);
    }

  virtual ~PairCutTrackAttrSameLabel() {}
};


/// Cut on Minv of the pair
///
/// This assumes highly relativistic particles and does not need an
/// assumed particle mass.
///
struct PairCutTrackAttrMinv {
protected:
  /// Square of minv range
  std::pair<double, double> fMinvSqrRange;
  double fMass1;
  double fMass2;

public:

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      // const double minv2 = CalcMinvSqrd(track1.P(), track2.P(), fMass1, fMass2);
      const double minv2 = CalcMinvSqrd(track1.P(), track2.P());
      return fMinvSqrRange.first <= minv2 && minv2 < fMinvSqrRange.second;
    }

  /// Minv assuming highly relativistic particles (E >> m)
  ///
  /// Does not require need mass information
  ///
  /// $ M_{inv}^2 = 2 p_{T,1} p_{T,2} (\cosh(\eta_1 - \eta_2) - \cos(\phi_1 - \phi_2)) $
  ///
  static double CalcMinvSqrd(const AliFemtoThreeVector &p1,
                             const AliFemtoThreeVector &p2)
    {
      const double
        pt1 = p1.Perp2(),
        eta1 = p1.PseudoRapidity(),
        phi1 = p1.Phi(),

        pt2 = p2.Perp2(),
        eta2 = p2.PseudoRapidity(),
        phi2 = p2.Phi();

      return 2 * std::sqrt(pt1 * pt2) * (std::cosh(eta1 - eta2) - std::cos(phi1 - phi2));
    }

  /// Minv with assumed masses
  ///
  /// $ M_{inv}^2 = m_1^2 m_2^2 + 2 (E_1 E_2 - \vec{p_1} \cdot \vec{p_2}) $
  ///
  static double CalcMinvSqrd(const AliFemtoThreeVector &p1,
                             const AliFemtoThreeVector &p2,
                             double mass1,
                             double mass2)
    {
      const double
        m1sqrd = mass1*mass1,
        m2sqrd = mass2*mass2,
        e1sqrd = m1sqrd + p1.Mag2(),
        e2sqrd = m2sqrd + p2.Mag2(),
        E1E2 = std::sqrt(e1sqrd*e2sqrd);

      return m1sqrd + m2sqrd + 2*(E1E2 - p1.Dot(p2));
    }

  /// Minv with assumed masses
  ///
  /// $ M_{inv}^2 = m_1^2 m_2^2 + 2 (E_1 E_2 - \vec{p_1} \cdot \vec{p_2}) $
  ///
  static double CalcMinvSqrd(const AliFemtoLorentzVector &p1,
                             const AliFemtoLorentzVector &p2)
    {
      return (p1 + p2).m2();
    }

  void SetMinvRange(double lo, double hi)
    {
      fMinvSqrRange = std::make_pair(lo * lo, hi * hi);
    }

  void SetMinvRange(const std::pair<double, double> &range)
    {
      double lo = range.first * range.first,
             hi = range.second * range.second;
      fMinvSqrRange = std::make_pair(lo, hi);
    }

  std::pair<double, double>
  GetMinvRange() const
    {
      return std::make_pair(std::sqrt(fMinvSqrRange.first),
                            std::sqrt(fMinvSqrRange.second));
    }

  PairCutTrackAttrMinv()
    : fMinvSqrRange(0, 100)
    , fMass1(0.0)
    , fMass2(0.0)
    {}

  PairCutTrackAttrMinv(AliFemtoConfigObject &cfg)
    : fMinvSqrRange(cfg.pop_range("minv_range", std::make_pair(0.0, 100.0)))
    , fMass1(0.0)
    , fMass2(0.0)
    {
      fMinvSqrRange.first *= fMinvSqrRange.first;
      fMinvSqrRange.second *= fMinvSqrRange.second;
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("minv_range", GetMinvRange());
    }

  virtual ~PairCutTrackAttrMinv() {}
};


/// \class PairCutTrackAttrRemoveEE
/// \brief Cut pairs with Minv near electron mass
///
struct PairCutTrackAttrRemoveEE {
  float ee_minv_min;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2) const
    {
      if (ee_minv_min == 0.0) {
        return true;
      }

      const double
        E_MASS = 0.000511,
        minv_sqrd = PairCutTrackAttrMinv::CalcMinvSqrd(track1.P(), track2.P(), E_MASS, E_MASS),
        minv = std::sqrt(minv_sqrd);

      return std::abs(minv - E_MASS) >= ee_minv_min;
    }

  PairCutTrackAttrRemoveEE()
    : ee_minv_min(0.0)
    {}

  PairCutTrackAttrRemoveEE(AliFemtoConfigObject &cfg)
    : ee_minv_min(cfg.pop_num("ee_minv_min", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      if (ee_minv_min != 0.0) {
        cfg.insert("ee_minv_min", ee_minv_min);
      }
    }

  virtual ~PairCutTrackAttrRemoveEE() {}
};

}  // namespace pwgcf



#include "AliFemtoPairCut.h"


/// \class AliFemtoPairCutAttrTracks
/// \brief Bridge from AliFemtoPairCut to a metaclass of PairCut-Attrs
///
/// Note - This expects two tracks, not an AliFemtoPair
///
/// Subclass and implement your method:
///  `const char* GetName() const`
///
///
template <typename CRTP, typename CutAttrType>
class AliFemtoPairCutAttrTracks : public AliFemtoPairCut, public CutAttrType {
public:

  typedef CutAttrType CutAttrs;

  virtual ~AliFemtoPairCutAttrTracks()
    { }

  AliFemtoPairCutAttrTracks()
    : AliFemtoPairCut()
    , CutAttrType()
    {}

  AliFemtoPairCutAttrTracks(AliFemtoConfigObject &cfg)
    : AliFemtoPairCut()
    , CutAttrType(cfg)
    {}

  /// user-written method to return string describing cuts
  virtual AliFemtoString Report()
    { return ""; }

  /// Return a TList of settings
  virtual TList* ListSettings()
    {
      TList* list = new TList();
      AppendSettings(*list);
      return list;
    }

  virtual void AppendSettings(TCollection &) const = 0;

  virtual bool Pass(const AliFemtoPair *pair)
    {
      return CutAttrs::Pass(*pair);
    }

  void StoreConfiguration(AliFemtoConfigObject &cfg) const
    {
      CutAttrs::FillConfiguration(cfg);
      cfg.insert("_class", static_cast<CRTP*>(this)->GetName());
    }

  AliFemtoConfigObject GetConfiguration() const
    {
      AliFemtoConfigObject result = AliFemtoConfigObject::BuildMap()
                                      ("_class", static_cast<const CRTP*>(this)->ClassName());
      CutAttrs::FillConfiguration(result);
      return result;
    }

};


#endif
