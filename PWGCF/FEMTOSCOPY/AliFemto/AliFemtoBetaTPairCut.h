///
/// \file AliFemtoBetaTPairCut.h
/// \author Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch
///
/// \class AliFemtoBetaTPairCut
/// \brief A pair cut which selects paris based on their betaT value
///
///

#ifndef ALIFEMTOBETATPAIRCUT_H
#define ALIFEMTOBETATPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoBetaTPairCut : public AliFemtoPairCut {
public:
  /// Construct pair cut with default values - masses of both particles are
  /// equivalent to pions, and all betaT values from 0 to 1 million are passed
  ///
  AliFemtoBetaTPairCut();

  /// Construct pair cut with betaT values and expected particle masses
  AliFemtoBetaTPairCut(double minbetat, double maxbetat, double masspart1, double masspart2);

  /// Copy constructor
  AliFemtoBetaTPairCut(const AliFemtoBetaTPairCut& c);

  /// Destructor
  virtual ~AliFemtoBetaTPairCut();
  AliFemtoBetaTPairCut& operator=(const AliFemtoBetaTPairCut& c);

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
  ClassDef(AliFemtoBetaTPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoPairCut* AliFemtoBetaTPairCut::Clone() {
  AliFemtoBetaTPairCut* cPairCut = new AliFemtoBetaTPairCut(*this);
  return cPairCut;
}

#endif
