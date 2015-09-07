/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoBetaTPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOBETATPAIRCUT_H
#define ALIFEMTOBETATPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoBetaTPairCut : public AliFemtoPairCut{
public:
  AliFemtoBetaTPairCut();
  AliFemtoBetaTPairCut(double min, double max);
  AliFemtoBetaTPairCut(const AliFemtoBetaTPairCut& c);
  virtual ~AliFemtoBetaTPairCut();
  AliFemtoBetaTPairCut& operator=(const AliFemtoBetaTPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetBetaTRange(double minbetat, double maxbetat);
  virtual bool Pass(const AliFemtoPair* pair);

 protected:
  Double_t fBetaTMin;          // Minimum allowed BetaT
  Double_t fBetaTMax;          // Maximum allowed BetaT

#ifdef __ROOT__
  ClassDef(AliFemtoBetaTPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoBetaTPairCut::Clone() { AliFemtoBetaTPairCut* cPairCut = new AliFemtoBetaTPairCut(*this); return cPairCut;}

#endif
