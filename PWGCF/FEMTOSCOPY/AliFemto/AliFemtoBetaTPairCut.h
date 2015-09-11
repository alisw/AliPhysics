//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// AliFemtoBetaTPairCut - a pair cut which selects paris based on their     //
// betaT value                                                              //
//                                                                          //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOBETATPAIRCUT_H
#define ALIFEMTOBETATPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoBetaTPairCut : public AliFemtoPairCut{
public:
  AliFemtoBetaTPairCut();
  AliFemtoBetaTPairCut(double minbetat, double maxbetat, double masspart1, double masspart2);
  AliFemtoBetaTPairCut(const AliFemtoBetaTPairCut& c);
  virtual ~AliFemtoBetaTPairCut();
  AliFemtoBetaTPairCut& operator=(const AliFemtoBetaTPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetBetaTRange(double minbetat, double maxbetat);
  void SetParticleMasses(double masspart1, double masspart2);
  virtual bool Pass(const AliFemtoPair* pair);

 protected:
  Double_t fBetaTMin;          // Minimum allowed BetaT
  Double_t fBetaTMax;          // Maximum allowed BetaT
  Double_t fMassPart1;         // Mass of the first particle in pair [GeV]
  Double_t fMassPart2;         // Mass of the second particle in pair [GeV]

#ifdef __ROOT__
  ClassDef(AliFemtoBetaTPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoBetaTPairCut::Clone() { AliFemtoBetaTPairCut* cPairCut = new AliFemtoBetaTPairCut(*this); return cPairCut;}

#endif
