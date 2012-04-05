/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityPairCut.h 24360 2008-03-10 09:48:27Z akisiel $
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


#ifndef ALIFEMTOV0TRACKPAIRCUT_H
#define ALIFEMTOV0TRACKPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoV0TrackPairCut : public AliFemtoPairCut{
public:
  AliFemtoV0TrackPairCut();
  AliFemtoV0TrackPairCut(const AliFemtoV0TrackPairCut& cut);
  virtual ~AliFemtoV0TrackPairCut();
  
  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetV0Max(Double_t aAliFemtoV0Max);
  Double_t GetAliFemtoV0Max() const;
  void     SetRemoveSameLabel(Bool_t aRemove);
  void SetTPCOnly(Bool_t tpconly);
  void SetShareQualityMax(Double_t aShareQualityMax);
  void SetShareFractionMax(Double_t aShareFractionMax);
  
 protected:
  long fNPairsPassed;          // Number of pairs consideered that passed the cut 
  long fNPairsFailed;          // Number of pairs consideered that failed the cut

 private:
  Double_t fV0Max;   // Maximum allowed pair quality
  Double_t fShareQualityMax;
  Double_t fShareFractionMax;  // Maximum allowed share fraction
  Bool_t   fRemoveSameLabel;   // If 1 pairs with two tracks with the same label will be removed
  Bool_t fTrackTPCOnly;


#ifdef __ROOT__
  ClassDef(AliFemtoV0TrackPairCut, 0)
#endif
};

inline AliFemtoV0TrackPairCut::AliFemtoV0TrackPairCut(const AliFemtoV0TrackPairCut& c) : 
  AliFemtoPairCut(c),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareQualityMax(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0),
  fTrackTPCOnly(0)
{ /* no-op */ }

inline AliFemtoPairCut* AliFemtoV0TrackPairCut::Clone() { AliFemtoV0TrackPairCut* c = new AliFemtoV0TrackPairCut(*this); return c;}

#endif
