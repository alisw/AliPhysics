/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoKTPairCut - a pair cut which selects pairs based on their       //
// transverse momentum kT                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoKTPairCut.h,v 1.1.2.1 2007/10/19 13:28:14 akisiel Exp $
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

#ifndef ALIFEMTOKTPAIRCUT_H
#define ALIFEMTOKTPAIRCUT_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"

class AliFemtoKTPairCut : public AliFemtoPairCut{
public:
  AliFemtoKTPairCut();
  AliFemtoKTPairCut(double lo, double hi);
  AliFemtoKTPairCut(const AliFemtoKTPairCut& c);
  virtual ~AliFemtoKTPairCut();
  AliFemtoKTPairCut& operator=(const AliFemtoKTPairCut& c);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetKTRange(double ktmin, double ktmax);
  void SetPhiRange(double phimin, double phimax);
  void SetPTMin(double ptmin, double ptmax=1000.0);
  virtual bool Pass(const AliFemtoPair* pair);
  virtual bool Pass(const AliFemtoPair* pair, double aRPAngle);

 protected:
  Double_t fKTMin;          // Minimum allowed pair transverse momentum
  Double_t fKTMax;          // Maximum allowed pair transverse momentum 
  Double_t fPhiMin;         // Minimum angle vs. reaction plane 
  Double_t fPhiMax;         // Maximum angle vs. reaction plane
  Double_t fPtMin;          // Minimum per-particle pT
  Double_t fPtMax;          // Maximum per-particle pT

#ifdef __ROOT__
  ClassDef(AliFemtoKTPairCut, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoKTPairCut::Clone() { AliFemtoKTPairCut* c = new AliFemtoKTPairCut(*this); return c;}

#endif
