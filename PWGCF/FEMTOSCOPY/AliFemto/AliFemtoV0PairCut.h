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


#ifndef ALIFEMTOV0PAIRCUT_H
#define ALIFEMTOV0PAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoV0PairCut : public AliFemtoPairCut{
public:
  AliFemtoV0PairCut();
  AliFemtoV0PairCut(const AliFemtoV0PairCut& cut);
  virtual ~AliFemtoV0PairCut();
  AliFemtoV0PairCut& operator=(const AliFemtoV0PairCut& cut);
  
  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetV0Max(Double_t aAliFemtoV0Max);
  Double_t GetAliFemtoV0Max() const;
  void     SetRemoveSameLabel(Bool_t aRemove);
  void SetTPCEntranceSepMinimum(double dtpc);
  void SetTPCExitSepMinimum(double dtpc);
  void SetDataType(AliFemtoDataType type);

 protected:
  long fNPairsPassed;          // Number of pairs consideered that passed the cut 
  long fNPairsFailed;          // Number of pairs consideered that failed the cut
  Double_t fV0Max;   // Maximum allowed pair quality
  Double_t fShareFractionMax;  // Maximum allowed share fraction
  Bool_t   fRemoveSameLabel;   // If 1 pairs with two tracks with the same label will be removed 

  AliFemtoDataType fDataType; //Use ESD / AOD / Kinematics.
  Double_t fDTPCMin;          // Minimum allowed pair nominal separation at the entrance to the TPC
  Double_t fDTPCExitMin;          // Minimum allowed pair nominal separation at the exit of the TPC

#ifdef __ROOT__
  ClassDef(AliFemtoV0PairCut, 0)
#endif
};

inline AliFemtoV0PairCut::AliFemtoV0PairCut(const AliFemtoV0PairCut& c) : 
  AliFemtoPairCut(c),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fV0Max(1.0),
  fShareFractionMax(1.0),
  fRemoveSameLabel(0),
  fDataType(kAOD),
  fDTPCMin(0),
  fDTPCExitMin(0)		      

{ /* no-op */ }

inline AliFemtoPairCut* AliFemtoV0PairCut::Clone() { AliFemtoV0PairCut* c = new AliFemtoV0PairCut(*this); return c;}

#endif
