/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityTPCEntranceSepPairCut - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityTPCEntranceSepPairCut.h,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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


#ifndef ALIFEMTOSHAREQUALITYTPCENTRANCESEPPAIRCUT_H
#define ALIFEMTOSHAREQUALITYTPCENTRANCESEPPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"

class AliFemtoShareQualityTPCEntranceSepPairCut : public AliFemtoShareQualityPairCut{
public:
  AliFemtoShareQualityTPCEntranceSepPairCut();
  AliFemtoShareQualityTPCEntranceSepPairCut(const AliFemtoShareQualityTPCEntranceSepPairCut& c);
  virtual ~AliFemtoShareQualityTPCEntranceSepPairCut();

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetTPCEntranceSepMinimum(double dtpc);
  
 protected:
  Double_t fDTPCMin;          // Minimum allowed pair nominal separation at the entrance to the TPC

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityTPCEntranceSepPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoShareQualityTPCEntranceSepPairCut::Clone() { AliFemtoShareQualityTPCEntranceSepPairCut* c = new AliFemtoShareQualityTPCEntranceSepPairCut(*this); return c;}

#endif
