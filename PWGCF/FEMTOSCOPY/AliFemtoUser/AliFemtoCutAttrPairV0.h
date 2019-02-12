///
/// \file AliFemtoUser/AliFemtoCutAttrPairV0.h
///

#pragma once

#ifndef ALIFEMTOCUTATTRPAIRV0_H
#define ALIFEMTOCUTATTRPAIRV0_H


#include "AliFemtoConfigObject.h"
#include "AliFemtoV0.h"


namespace pwgfemto {

/// Cut on the sum of the PT
template <typename T1, typename T2>
struct AddPairCutV0Attr : public T1, public T2 {

  bool Pass(const AliFemtoV0 &v0_1, const AliFemtoV0 &v0_2)
    {
      return T1::Pass(v0_1, v0_2) && T2::Pass(v0_1, v0_2);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      T1::FillConfiguration(cfg);
      T2::FillConfiguration(cfg);
    }

  virtual ~AddPairCutV0Attr() = 0;
};


/// Cut on the sum of the PT
///
struct PairCutV0AttrPt {

  std::pair<double, double> pt_range;

  bool Pass(const AliFemtoV0 &trk1, const AliFemtoV0 &trk2)
    {
      const double pt_sum = trk1.PtV0() + trk2.PtV0();
      return pt_range.first <= pt_sum && pt_sum < pt_range.second;
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("pt_range", pt_range);
    }

  virtual ~PairCutV0AttrPt() = 0;
};

/// Simple cut removing pairs sharing a daugther track
struct PairCutV0SharedDaughters {

  bool Pass(const AliFemtoV0 &trk1, const AliFemtoV0 &trk2)
    {
      return trk1.IdNeg() != trk2.IdNeg()
          && trk1.IdPos() != trk2.IdPos();
    }

  void FillConfiguration(AliFemtoConfigObject &) const
    { }

  virtual ~PairCutV0SharedDaughters() = 0;
};


/// \class PairCutV0MinEntranceSep
/// \brief Accept pairs where TPC entrance separation is larger than
///        some parameter
///
struct PairCutV0MinEntranceSep {

  double min_entrance_sep;

  bool Pass(const AliFemtoV0 &trk1, const AliFemtoV0 &trk2)
    {
      const auto
        dpos = trk1.NominalTpcEntrancePointPos() - trk2.NominalTpcEntrancePointPos(),
        dneg = trk1.NominalTpcEntrancePointNeg() - trk2.NominalTpcEntrancePointNeg();

      const double
        min_sep_sqrd = min_entrance_sep * min_entrance_sep;

       return min_sep_sqrd <= dpos.Mag2()
           && min_sep_sqrd <= dneg.Mag2();
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("min_entrance_sep", min_entrance_sep);
    }

  virtual ~PairCutV0MinEntranceSep() = 0;
};

/// \class PairCutV0AvgSep
/// \brief Cut pairs with child-tracks likely split track
///
struct PairCutV0AvgSep {

  double min_avgsep_nn,
         min_avgsep_np,
         min_avgsep_pn,
         min_avgsep_pp;

  PairCutV0AvgSep()
    : min_avgsep_nn(0)
    , min_avgsep_np(0)
    , min_avgsep_pn(0)
    , min_avgsep_pp(0)
    {
    }

  bool Pass(const AliFemtoV0 &trk1, const AliFemtoV0 &trk2)
    {
      return min_avgsep_nn <= avg_separation_child(trk1, false, trk2, false)
          && min_avgsep_np <= avg_separation_child(trk1, false, trk2, true)
          && min_avgsep_pn <= avg_separation_child(trk1, true, trk2, false)
          && min_avgsep_pp <= avg_separation_child(trk1, true, trk2, true);
    }

  void FillConfiguration(AliFemtoConfigObject &cfg) const
    {
      cfg.insert("min_avgsep_nn", min_avgsep_nn);
      cfg.insert("min_avgsep_np", min_avgsep_np);
      cfg.insert("min_avgsep_pn", min_avgsep_pn);
      cfg.insert("min_avgsep_pp", min_avgsep_pp);
    }

  static bool is_unset_vector(const AliFemtoThreeVector& v)
    {
      return v.x() <= -9999.0
          || v.y() <= -9999.0
          || v.z() <= -9999.0;
    }

  static double
  avg_separation_child(const AliFemtoV0& V1, bool pos_child_1,
                       const AliFemtoV0& V2, bool pos_child_2)
    {
      int counter = 0;
      double sep = 0.0;

      for (int i=0; i < 8; i++) {

        const AliFemtoThreeVector
          p1 = pos_child_1
             ? V1.NominalTpcPointPos(counter)
             : V1.NominalTpcPointNeg(counter),

          p2 = pos_child_2
             ? V2.NominalTpcPointPos(counter)
             : V2.NominalTpcPointNeg(counter);

        if (is_unset_vector(p1) || is_unset_vector(p2)) {
          continue;
        }

        sep += (p1 - p2).Mag();
        ++counter;
      }

      return counter==0 ? 0.0 : sep / counter;
    }

  virtual ~PairCutV0AvgSep() = 0;
};


}  // namespace

#endif
