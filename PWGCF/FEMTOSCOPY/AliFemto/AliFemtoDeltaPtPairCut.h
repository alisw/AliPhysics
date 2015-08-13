/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoKTPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoDeltaPtPairCut.h,v 1.1.2.1 2007/10/19 13:28:14 akisiel Exp $
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

#ifndef ALIFEMTODELTAPTPAIRCUT_H
#define ALIFEMTODELTAPTPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoDeltaPtPairCut : public AliFemtoPairCut{
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
  Double_t fDeltaPtMin;          // Minimum allowed pair transverse momentum
  Double_t fDeltaPtMax;          // Maximum allowed pair transverse momentum 


#ifdef __ROOT__
  ClassDef(AliFemtoDeltaPtPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoDeltaPtPairCut::Clone() { AliFemtoDeltaPtPairCut* c = new AliFemtoDeltaPtPairCut(*this); return c;}

#endif
