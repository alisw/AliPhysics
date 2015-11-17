///
/// \file AliFemtoCorrFctnKStar.h
///

#ifndef ALIFEMTOCORRFCTN_KSTAR_H_
#define ALIFEMTOCORRFCTN_KSTAR_H_

class TH1D;
class TH2F;


#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrFctnKStar
/// \brief A simple correlation function which plots the kStar between pairs
///
/// \author Andrew Kubera, Ohio State University, <andrew.kubera@cern.ch>
///
class AliFemtoCorrFctnKStar : public AliFemtoCorrFctn {
public:

  /// Construct default parameters
  AliFemtoCorrFctnKStar();

  /// Construct with histogram parameters
  AliFemtoCorrFctnKStar(const char *title,
                        const int nbins,
                        const float KStarLo,
                        const float KStarHi);

  /// Create runtime report
  virtual AliFemtoString Report();

  /// GetOutputList
  virtual TList* GetOutputList();

  /// Run after analysis is finished
  virtual void Finish();

  /// Add K* of pairs from an event
  virtual void AddRealPair(AliFemtoPair* aPair);

  /// Add K* of pairs from mied events
  virtual void AddMixedPair(AliFemtoPair* aPair);

  /// calculate $m_{T}$ of pair
  static float CalcMt(const AliFemtoPair* aPair);

protected:

  /// K* of pairs in same event
  TH1D *fNumerator;

  /// K* of pairs in different event
  TH1D *fDenominator;

  /// K* vs kT
  TH2F *fNumerator_kT;
  TH2F *fDenominator_kT;

  /// K* vs mT
  TH2F *fNumerator_mT;
  TH2F *fDenominator_mT;

  /// K* vs mT
  TH2F *fNumerator_qq;
  TH2F *fDenominator_qq;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnKStar, 0);
  /// \endcond
#endif
};

#endif /* ALIFEMTOCORRFCTN_KSTAR_H_ */
