/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
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


#ifndef ALIFEMTOSHAREQUALITYPAIRCUT_H
#define ALIFEMTOSHAREQUALITYPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoShareQualityPairCut : public AliFemtoPairCut{
public:
  AliFemtoShareQualityPairCut();
  AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut& cut);
  virtual ~AliFemtoShareQualityPairCut();
  
  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoShareQualityPairCut* Clone();
  void SetShareQualityMax(Double_t aAliFemtoShareQualityMax);
  Double_t GetAliFemtoShareQualityMax() const;
  void SetShareFractionMax(Double_t aAliFemtoShareFractionMax);
  Double_t GetAliFemtoShareFractionMax() const;
  void     SetRemoveSameLabel(Bool_t aRemove);
  
 protected:
  long fNPairsPassed;          // Number of pairs consideered that passed the cut 
  long fNPairsFailed;          // Number of pairs consideered that failed the cut

 private:
  Double_t fShareQualityMax;   // Maximum allowed pair quality
  Double_t fShareFractionMax;  // Maximum allowed share fraction
  Bool_t   fRemoveSameLabel;   // If 1 pairs with two tracks with the same label will be removed 


#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityPairCut, 0)
#endif
};

inline AliFemtoShareQualityPairCut::AliFemtoShareQualityPairCut(const AliFemtoShareQualityPairCut& c) : 
  AliFemtoPairCut(c),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0)// no cut
{ /* no-op */ }

inline AliFemtoShareQualityPairCut* AliFemtoShareQualityPairCut::Clone() { AliFemtoShareQualityPairCut* c = new AliFemtoShareQualityPairCut(*this); return c;}

#endif
