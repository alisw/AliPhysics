///
/// \file AliFemtokTPairCut.h
/// \author Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch
///
/// \class AliFemtokTPairCut
/// \brief A pair cut which selects paris based on their betaT value
///
///

#ifndef ALIFEMTOKTPAIRCUT_H
#define ALIFEMTOKTPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtokTPairCut : public AliFemtoPairCut {
public:
  /// Construct pair cut with default values - masses of both particles are
  /// equivalent to pions, and all betaT values from 0 to 1 million are passed
  ///
  AliFemtokTPairCut();

  /// Construct pair cut with betaT values and expected particle masses
  AliFemtokTPairCut(double minbetat, double maxbetat, double masspart1, double masspart2);

  /// Copy constructor
  AliFemtokTPairCut(const AliFemtokTPairCut& c);

  /// Destructor
  virtual ~AliFemtokTPairCut();
  AliFemtokTPairCut& operator=(const AliFemtokTPairCut& c);

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
  ClassDef(AliFemtokTPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoPairCut* AliFemtokTPairCut::Clone() {
  AliFemtokTPairCut* cPairCut = new AliFemtokTPairCut(*this);
  return cPairCut;
}

#endif
