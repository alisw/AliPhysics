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


/// \class AddPairCutAttrs
/// \brief Join cut-types together
///
template <typename T1, typename T2>
struct AddPairCutAttrs : public T1, public T2 {

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      return T1::Pass(track1, track2) && T2::Pass(track1, track2);
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

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      return calc_avg_sep(track1, track2) > avgsep_min;
    }

  PairCutTrackAttrAvgSep()
    : avgsep_min(0.0)
    {
    }

  PairCutTrackAttrAvgSep(AliFemtoConfigObject &cfg)
    : avgsep_min(cfg.pop_float("avgsep_min", 0.0))
    {
    }

  static bool is_valid_point(const AliFemtoThreeVector &point)
    {
      // defined in AliFemtoEventReader::GetGlobalPositionAtGlobalRadiiThroughTPC
      return point.x() <= -9000.0
          && point.y() <= -9000.0
          && point.z() <= -9000.0;
    }

  static double calc_avg_sep(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      int i = 0;
      double sep_sum = 0.0;

      for (; i<8; ++i) {
        const auto &p1 = track1.NominalTpcPoint(i),
                   &p2 = track2.NominalTpcPoint(i);

        if (!is_valid_point(p1) || !is_valid_point(p2)) {
          break;
        }
        sep_sum += (p1 - p2).Mag();
      }
      return sep_sum / i;
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("avgsep_min", avgsep_min);
    }

  virtual ~PairCutTrackAttrAvgSep() {}
};


/// Cut on the sum of the PT
struct PairCutTrackAttrPt {

  static const std::pair<double, double> DEFAULT;
  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
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

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      const std::pair<double, double>
        qual_and_frac = calc_share_quality_fraction(track1, track2);

      const double
        share_fraction = qual_and_frac.first,
        share_quality = qual_and_frac.second;

      return share_fraction <= share_fraction_max && share_quality <= share_quality_max;
    }

  static std::pair<double, double>
  calc_share_quality_fraction(const AliFemtoTrack &track1,
                              const AliFemtoTrack &track2)
    {
      const TBits
        &tpc_clusters_1 = track1.TPCclusters(),
        &tpc_sharing_1 = track1.TPCsharing(),
        &tpc_clusters_2 = track2.TPCclusters(),
        &tpc_sharing_2 = track2.TPCsharing(),

        cls1_and_cls2 = tpc_clusters_1 & tpc_clusters_2,
        cls1_xor_cls2 = tpc_clusters_1 ^ tpc_clusters_2,
        shr1_and_shr2 = tpc_sharing_1 & tpc_sharing_2,
        not_sharing = ~const_cast<TBits&>(shr1_and_shr2);

      const double
        cls1_xor_cls2_bits = cls1_xor_cls2.CountBits(),
        nh = cls1_xor_cls2_bits + 2 * cls1_and_cls2.CountBits(),
        ns = 2 * (cls1_and_cls2 & shr1_and_shr2).CountBits(),
        an = cls1_xor_cls2_bits + ns / 2 - (cls1_and_cls2 & not_sharing).CountBits(),

        share_quality = an / nh,
        share_fraction = ns / nh;

      return std::make_pair(share_quality, share_fraction);
    }

  PairCutTrackAttrShareQuality()
    : share_fraction_max(1.0)
    , share_quality_max(1.0)
    {}

  PairCutTrackAttrShareQuality(AliFemtoConfigObject &cfg)
    : share_fraction_max(cfg.pop_float("share_fraction_max", 1.0))
    , share_quality_max(cfg.pop_float("share_quality_max", 1.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("share_quality_max", share_quality_max);
      cfg.insert("share_fraction_max", share_fraction_max);
    }

  virtual ~PairCutTrackAttrShareQuality() {}
};


/// \class PairCutTrackAttrDetaDphiStar
/// \brief A pair cut which cuts on the Δη Δφ of the pair.
///
/// The difference in phi is calculated by examining the tracks' azimuthal angle at a particular radius
/// of all tracks. The default value for this radius is 1.6 (inside the TPC) but this may be changed via
/// the SetR method.
///
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
   , phistar_radius(1.2)
   , fCurrentMagneticField(0.0)
   {
   }

  PairCutTrackAttrDetaDphiStar(AliFemtoConfigObject &cut)
   : delta_eta_min(cut.pop_float("delta_eta_min", 0.0))
   , delta_phistar_min(cut.pop_float("delta_phistar_min", 0.0))
   , phistar_radius(cut.pop_float("phistar_radius", 1.2))
   , fCurrentMagneticField(0.0)
   {
   }

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      const AliFemtoThreeVector &p1 = track1.P(),
                                &p2 = track2.P();

      const double deta = calc_delta_eta(p1, p2);
      if (delta_eta_min < deta) {
        return true;
      }

      const double dphi = AliFemtoPairCutDetaDphi::CalculateDPhiStar(
                            p1, track1.Charge(),
                            p2, track2.Charge(),
                            phistar_radius,
                            fCurrentMagneticField);

      return delta_phistar_min * delta_phistar_min < deta * deta + dphi * dphi;
    }

  static double calc_delta_eta(const AliFemtoThreeVector &p1, const AliFemtoThreeVector &p2)
    {
      return p1.PseudoRapidity() - p2.PseudoRapidity();
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.Update(AliFemtoConfigObject::BuildMap()
                 ("delta_eta_min", delta_eta_min)
                 ("delta_phistar_min", delta_phistar_min)
                 ("phistar_radius", delta_phistar_min));
    }

  virtual ~PairCutTrackAttrDetaDphiStar() {}
};


/// Remove pairs of tracks with the same particle
struct PairCutTrackAttrSameLabel {

  bool remove_same_label;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      return remove_same_label ? abs(track1.Label()) == abs(track2.Label()) : true;
    }

  PairCutTrackAttrSameLabel()
    : remove_same_label(true)
    {}

  PairCutTrackAttrSameLabel(AliFemtoConfigObject &cfg)
    : remove_same_label(cfg.pop_bool("remove_same_label", true))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("remove_same_label", remove_same_label);
    }

  virtual ~PairCutTrackAttrSameLabel() {}
};


struct PairCutTrackAttrMinv {

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

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
      // const double minv2 = CalcMinvSqrd(track1.P(), track2.P(), fMass1, fMass2);
      const double minv2 = CalcMinvSqrd(track1.P(), track2.P());
      return fMinvSqrRange.first <= minv2 && minv2 < fMinvSqrRange.second;
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

protected:
  /// Square of minv range
  std::pair<double, double> fMinvSqrRange;
  double fMass1;
  double fMass2;

  virtual ~PairCutTrackAttrMinv() {}
};


/// \class PairCutTrackAttrRemoveEE
/// \brief Cut pairs with Minv near electron mass
struct PairCutTrackAttrRemoveEE {
  float ee_minv_min;

  bool Pass(const AliFemtoTrack &track1, const AliFemtoTrack &track2)
    {
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
    : ee_minv_min(cfg.pop_float("ee_minv_min", 0.0))
    {}

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("ee_minv_min", ee_minv_min);
    }

  virtual ~PairCutTrackAttrRemoveEE() {}
};

}  // namespace pwgcf



#include "AliFemtoPairCut.h"


/// \class AliFemtoPairCutAttrTracks
/// \brief Bridge from AliFemtoPairCut to a metaclass of PairCut-Attrs
///
/// Note - This expects two tracks
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
      return CutAttrs::Pass(*pair->Track1()->Track(), *pair->Track2()->Track());
    }

  void StoreConfiguration(AliFemtoConfigObject &cfg) const
    {
      CutAttrs::FillConfiguration(cfg);
      cfg.insert("_class", static_cast<CRTP*>(this)->GetName());
    }

  AliFemtoConfigObject GetConfiguration() const
    {
      AliFemtoConfigObject result = AliFemtoConfigObject::BuildMap()
                                      ("_class", CRTP::ClassName());
      CutAttrs::FillConfiguration(result);
      return result;
    }

};




#endif
