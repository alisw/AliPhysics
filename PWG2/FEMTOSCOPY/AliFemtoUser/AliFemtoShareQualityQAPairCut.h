/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityQAPairCut - a pair cut which checks for some pair     //
// qualities that attempt to identify slit/doubly reconstructed tracks     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityQAPairCut.h 24360 2008-03-10 09:48:27Z akisiel $
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


#ifndef AliFemtoShareQualityQAPairCut_H
#define AliFemtoShareQualityQAPairCut_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoShareQualityQAPairCut : public AliFemtoPairCut{
public:
  AliFemtoShareQualityQAPairCut();
  AliFemtoShareQualityQAPairCut(const AliFemtoShareQualityQAPairCut& cut);
  virtual ~AliFemtoShareQualityQAPairCut();
  
  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetShareQualityMax(Double_t aAliFemtoShareQualityMax);
  void SetShareQualitymin(Double_t aAliFemtoShareQualitymin);
  void SetShareQualityQASwitch(bool aSwitch);
  void SetShareQualityQAExclusionZone(Double_t lo, Double_t hi);
  Double_t GetAliFemtoShareQualityMax() const;
  void SetShareFractionMax(Double_t aAliFemtoShareFractionMax);
  void SetShareFractionmin(Double_t aAliFemtoShareFractionmin);
  void SetShareFractionQASwitch(bool aSwitch);
  void SetShareFractionQAExclusionZone(Double_t lo, Double_t hi);
  Double_t GetAliFemtoShareFractionMax() const;
  void     SetRemoveSameLabel(Bool_t aRemove);
  
 protected:
  long fNPairsPassed;          // Number of pairs consideered that passed the cut 
  long fNPairsFailed;          // Number of pairs consideered that failed the cut

 private:
  Double_t fShareQualityMax;   // Maximum allowed pair quality
  Double_t fShareQualitymin;   // Minimum allowed pair quality
  Double_t fShareFractionMax;  // Maximum allowed share fraction
  Double_t fShareFractionmin;  // Minimum allowed share fraction
  Bool_t   fRemoveSameLabel;   // If 1 pairs with two tracks with the same label will be removed 
  
  bool     fShareQualityQASwitch;             // Turn on QA Exclusion Zone (true=on)
  Double_t fShareQualityQAExclusionZone[2];   // Set Exclusion Zone limits
  bool     fShareFractionQASwitch;             // Turn on QA Exclusion Zone (true=on)
  Double_t fShareFractionQAExclusionZone[2];   // Set Exclusion Zone limits


#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityQAPairCut, 0)
#endif
};

inline AliFemtoShareQualityQAPairCut::AliFemtoShareQualityQAPairCut(const AliFemtoShareQualityQAPairCut& c) : 
  AliFemtoPairCut(c),
  fNPairsPassed(0),
  fNPairsFailed(0),
  fShareQualityMax(1.0),
  fShareQualitymin(-0.5),
  fShareFractionMax(1.0),
  fShareFractionmin(0.0),
  fRemoveSameLabel(0)// no cut
{ 
  fShareQualityQASwitch  = c.fShareQualityQASwitch;
  fShareQualityQAExclusionZone[0]  = c.fShareQualityQAExclusionZone[0];
  fShareQualityQAExclusionZone[1]  = c.fShareQualityQAExclusionZone[1];
  fShareFractionQASwitch = c.fShareFractionQASwitch; 
  fShareFractionQAExclusionZone[0]  = c.fShareFractionQAExclusionZone[0];
  fShareFractionQAExclusionZone[1]  = c.fShareFractionQAExclusionZone[1];
}

inline AliFemtoPairCut* AliFemtoShareQualityQAPairCut::Clone() { AliFemtoShareQualityQAPairCut* c = new AliFemtoShareQualityQAPairCut(*this); return c;}

#endif
