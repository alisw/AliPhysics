/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoShareQualityTPCEntranceSepQAPairCut - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoShareQualityTPCEntranceSepQAPairCut.h,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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


#ifndef AliFemtoShareQualityTPCEntranceSepQAPairCut_H
#define AliFemtoShareQualityTPCEntranceSepQAPairCut_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityQAPairCut.h"

class AliFemtoShareQualityTPCEntranceSepQAPairCut : public AliFemtoShareQualityQAPairCut{
public:
  AliFemtoShareQualityTPCEntranceSepQAPairCut();
  AliFemtoShareQualityTPCEntranceSepQAPairCut(const AliFemtoShareQualityTPCEntranceSepQAPairCut& c);
  virtual ~AliFemtoShareQualityTPCEntranceSepQAPairCut();

  AliFemtoShareQualityTPCEntranceSepQAPairCut& operator=(const AliFemtoShareQualityTPCEntranceSepQAPairCut& aCut);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetTPCEntranceSepMinimum(double dtpc);
  void SetTPCEntranceSepMaximum(double dtpc);
  void SetTPCEntranceSepQASwitch(bool Switch);
  void SetTPCEntranceSepQAExclusionZone(double lo, double hi);
  
 protected:
  Double_t fDTPCMin;                // Minimum allowed pair nominal separation at the entrance to the TPC
  Double_t fDTPCMax;                // Maximum allowed pair nominal separation at the entrance to the TPC
  bool     fDTPCQASwitch;           // Turn on QA Exclusion Zone (true=on)
  Double_t fDTPCQAExclusionZone[2]; // Exclusion Zone for pair nominal separation at the entrance to the TPC

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityTPCEntranceSepQAPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoShareQualityTPCEntranceSepQAPairCut::Clone() { AliFemtoShareQualityTPCEntranceSepQAPairCut* c = new AliFemtoShareQualityTPCEntranceSepQAPairCut(*this); return c;}

#endif
