///
/// \file AliFemtoDeltaPtPairCut.h
/// \author Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
/// \date Oct 18, 2007
///
/// \class AliFemtoDeltaPtPairCut
/// \brief A pair cut which passes pairs solely on the difference of the
///        magnitude of their transverse momentum (pT)
///
/// The goal of this class is to remove "shared" and "split" pairs
///

#ifndef ALIFEMTODELTAPTPAIRCUT_H
#define ALIFEMTODELTAPTPAIRCUT_H

#include "AliFemtoPairCut.h"
#include <TList.h>

class AliFemtoDeltaPtPairCut : public AliFemtoPairCut {
public:
  AliFemtoDeltaPtPairCut();
  AliFemtoDeltaPtPairCut(double lo, double hi);
  AliFemtoDeltaPtPairCut(const AliFemtoDeltaPtPairCut& c);
  virtual ~AliFemtoDeltaPtPairCut();
  AliFemtoDeltaPtPairCut& operator=(const AliFemtoDeltaPtPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetDeltaPtRange(double ktmin, double ktmax);
  virtual bool Pass(const AliFemtoPair* pair);

protected:
  Double_t fDeltaPtMin; ///< Minimum allowed pair transverse momentum
  Double_t fDeltaPtMax; ///< Maximum allowed pair transverse momentum


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoDeltaPtPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoPairCut* AliFemtoDeltaPtPairCut::Clone()
{
  AliFemtoDeltaPtPairCut* c = new AliFemtoDeltaPtPairCut(*this);
  return c;
}

#endif
