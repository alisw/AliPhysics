/***************************************************************************
 *
 * $Id$
 *
 * Author: Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/


#ifndef AliFemtoShareQualityPairCut_hh
#define AliFemtoShareQualityPairCut_hh

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "Base/AliFemtoPairCut.h"

class AliFemtoShareQualityPairCut : public AliFemtoPairCut{
public:
  AliFemtoShareQualityPairCut();
  AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut&);
  ~AliFemtoShareQualityPairCut();

  virtual bool Pass(const AliFemtoPair*);
  virtual AliFemtoString Report();
  AliFemtoShareQualityPairCut* Clone();
  void SetShareQualityMax(Double_t aAliFemtoShareQualityMax);
  Double_t GetAliFemtoShareQualityMax();

private:
  long fNPairsPassed;
  long fNPairsFailed;
  Double_t fShareQualityMax;

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityPairCut, 0)
#endif
};

inline AliFemtoShareQualityPairCut::AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut& c) : AliFemtoPairCut(c) {
  fNPairsPassed = 0;
  fNPairsFailed = 0;
  fShareQualityMax = 1.0; // no cut
}

inline AliFemtoShareQualityPairCut* AliFemtoShareQualityPairCut::Clone() { AliFemtoShareQualityPairCut* c = new AliFemtoShareQualityPairCut(*this); return c;}

#endif
