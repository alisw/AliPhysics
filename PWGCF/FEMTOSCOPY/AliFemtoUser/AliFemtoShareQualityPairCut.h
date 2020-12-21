///
/// \file AliFemtoShareQualityPairCut.h
///

#pragma once

#ifndef ALIFEMTOSHAREQUALITYPAIRCUT_H
#define ALIFEMTOSHAREQUALITYPAIRCUT_H


#include "AliFemtoPairCut.h"

/// \class AliFemtoShareQualityPairCut
/// \brief A cut to remove "shared" and "split" pairs
///
/// A pair cut which checks for some pair qualities that
/// attempt to identify split/doubly-reconstructed tracks
///
/// \author Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
///
class AliFemtoShareQualityPairCut : public AliFemtoPairCut{
public:
  /// Build with share & quality limits set to 1.0
  AliFemtoShareQualityPairCut();
  /// Copy limits, pass/fail counts set to zero
  AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut& cut);
  virtual ~AliFemtoShareQualityPairCut();
  /// Copy limits, pass/fail counts reset to zero
  AliFemtoShareQualityPairCut& operator=(const AliFemtoShareQualityPairCut& cut);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoShareQualityPairCut* Clone() const;

  void SetShareQualityMax(Double_t aAliFemtoShareQualityMax);
  Double_t GetAliFemtoShareQualityMax() const;
  Double_t GetShareQualityMax() const;

  void SetShareFractionMax(Double_t aAliFemtoShareFractionMax);
  Double_t GetAliFemtoShareFractionMax() const;
  Double_t GetShareFractionMax() const;

  void SetRemoveSameLabel(Bool_t aRemove);
  Bool_t GetRemoveSameLabel() const;

  /// Putting the equality in sharequality
  bool operator==(const AliFemtoShareQualityPairCut &) const;

 protected:
  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut

 private:
  Double_t fShareQualityMax;   ///< Maximum allowed pair quality
  Double_t fShareFractionMax;  ///< Maximum allowed share fraction
  Bool_t   fRemoveSameLabel;   ///< Pairs with two tracks with the same label will be removed

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityPairCut, 0)
#endif
};

inline AliFemtoShareQualityPairCut::AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut& c) :
  AliFemtoPairCut(c),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(c.fShareQualityMax),
  fShareFractionMax(c.fShareFractionMax),
  fRemoveSameLabel(c.fRemoveSameLabel)
{ /* no-op */ }

inline AliFemtoShareQualityPairCut* AliFemtoShareQualityPairCut::Clone() const
{
  AliFemtoShareQualityPairCut* c = new AliFemtoShareQualityPairCut(*this);
  return c;
}

inline void AliFemtoShareQualityPairCut::SetShareQualityMax(Double_t aShareQualityMax) {
  fShareQualityMax = aShareQualityMax;
}
inline Double_t AliFemtoShareQualityPairCut::GetAliFemtoShareQualityMax() const {
  return fShareQualityMax;
}
inline Double_t AliFemtoShareQualityPairCut::GetShareQualityMax() const
{
  return fShareQualityMax;
}

inline void AliFemtoShareQualityPairCut::SetShareFractionMax(Double_t aShareFractionMax) {
  fShareFractionMax = aShareFractionMax;
}
inline Double_t AliFemtoShareQualityPairCut::GetAliFemtoShareFractionMax() const {
  return fShareFractionMax;
}
inline Double_t AliFemtoShareQualityPairCut::GetShareFractionMax() const
{
  return fShareFractionMax;
}

inline void AliFemtoShareQualityPairCut::SetRemoveSameLabel(Bool_t aRemove) {
  fRemoveSameLabel = aRemove;
}

inline Bool_t AliFemtoShareQualityPairCut::GetRemoveSameLabel() const
{
  return fRemoveSameLabel;
}

#endif
