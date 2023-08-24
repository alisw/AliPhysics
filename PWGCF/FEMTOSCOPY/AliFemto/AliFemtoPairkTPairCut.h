///
/// \file AliFemtoPairkTPairCut.h
/// \author Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch
///
/// \class AliFemtoPairkTPairCut
/// \brief A pair cut which selects paris based on their betaT value
///
///

#ifndef ALIFEMTOPAIRKTPAIRCUT_H
#define ALIFEMTOPAIRKTPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoPairkTPairCut : public AliFemtoPairCut {
public:
  /// Construct pair cut with default values - masses of both particles are
  /// equivalent to pions, and all betaT values from 0 to 1 million are passed
  ///
  AliFemtoPairkTPairCut();

  /// Construct pair cut with betaT values and expected particle masses
  AliFemtoPairkTPairCut(double minbetat, double maxbetat, double masspart1, double masspart2);

  /// Copy constructor
  AliFemtoPairkTPairCut(const AliFemtoPairkTPairCut& c);

  /// Destructor
  virtual ~AliFemtoPairkTPairCut();
  AliFemtoPairkTPairCut& operator=(const AliFemtoPairkTPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetBetaTRange(double minbetat, double maxbetat);
  void SetParticleMasses(double masspart1, double masspart2);
  virtual bool Pass(const AliFemtoPair* pair);

protected:
  Double_t fBetaTMin;   ///< Minimum allowed BetaT
  Double_t fBetaTMax;   ///< Maximum allowed BetaT
  Double_t fMassPart1;  ///< Mass of the first particle in pair [GeV]
  Double_t fMassPart2;  ///< Mass of the second particle in pair [GeV]

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoPairkTPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoPairCut* AliFemtoPairkTPairCut::Clone() {
  AliFemtoPairkTPairCut* cPairCut = new AliFemtoPairkTPairCut(*this);
  return cPairCut;
}

#endif
