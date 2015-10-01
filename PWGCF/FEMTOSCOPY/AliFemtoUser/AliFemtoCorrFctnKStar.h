///
/// \file AliFemtoCorrFctnKStar.h
///

#ifndef ALIFEMTOCORRFCTN_KSTAR_H_
#define ALIFEMTOCORRFCTN_KSTAR_H_

class TH1D;
#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrFctnKStar
/// \brief A simple correlation function which plots the kStar between pairs
///
/// \author Andrew Kubera, Ohio State University, <andrew.kubera@cern.ch>
///
class AliFemtoCorrFctnKStar : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnKStar();
  AliFemtoCorrFctnKStar(const char *title, const int nbins, const float KStarLo, const float KStarHi);

  virtual AliFemtoString Report();
  virtual void Finish();
  virtual TList* GetOutputList();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

protected:
  TH1D *fNumerator;
  TH1D *fDenominator;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnKStar, 0);
  /// \endcond
#endif
};

#endif /* ALIFEMTOCORRFCTN_KSTAR_H_ */
