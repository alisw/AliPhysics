///
/// \file AliFemtoKTPairCut.h
///

#pragma once

#ifndef ALIFEMTOKTPAIRCUT_H
#define ALIFEMTOKTPAIRCUT_H

#include "AliFemtoPairCut.h"


/// \class AliFemtoKTPairCut
/// \brief A pair cut which selects pairs based on their
///        transverse momentum kT
///
/// A cut to remove "shared" and "split" pairs
///
/// \author Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
///
class AliFemtoKTPairCut : public AliFemtoPairCut {
public:
  AliFemtoKTPairCut();
  AliFemtoKTPairCut(double lo, double hi);
  AliFemtoKTPairCut(const AliFemtoKTPairCut& c);
  virtual ~AliFemtoKTPairCut();
  AliFemtoKTPairCut& operator=(const AliFemtoKTPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetKTRange(double ktmin, double ktmax);
  void SetPhiRange(double phimin, double phimax);
  void SetPTMin(double ptmin, double ptmax=1000.0);
  virtual bool Pass(const AliFemtoPair* pair);
  virtual bool Pass(const AliFemtoPair* pair, double aRPAngle);

  std::pair<double, double> GetKtRange() const
    { return std::make_pair(fKTMin, fKTMax); }

  std::pair<double, double> GetPhiRange() const
    { return std::make_pair(fPhiMin, fPhiMax); }

  std::pair<double, double> GetPtRange() const
    { return std::make_pair(fPtMin, fPtMax); }

protected:

  Double_t fKTMin;          // Minimum allowed pair transverse momentum
  Double_t fKTMax;          // Maximum allowed pair transverse momentum
  Double_t fPhiMin;         // Minimum angle vs. reaction plane
  Double_t fPhiMax;         // Maximum angle vs. reaction plane
  Double_t fPtMin;          // Minimum per-particle pT
  Double_t fPtMax;          // Maximum per-particle pT

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoKTPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoPairCut* AliFemtoKTPairCut::Clone()
  { return new AliFemtoKTPairCut(*this); }

#endif
