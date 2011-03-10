/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutRadialDistance - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutRadialDistance.h,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
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


#ifndef ALIFEMTOPAIRCUTRADIALDISTANCE_H
#define ALIFEMTOPAIRCUTRADIALDISTANCE_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"

class AliFemtoPairCutRadialDistance : public AliFemtoShareQualityPairCut{
public:
  AliFemtoPairCutRadialDistance();
  AliFemtoPairCutRadialDistance(const AliFemtoPairCutRadialDistance& c);
  virtual ~AliFemtoPairCutRadialDistance();

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetRadialDistanceMinimum(double radius, double dtpc);
  void SetEtaDifferenceMinimum(double etpc);

 protected:
  Double_t fDRadMin;          // Minimum allowed pair separation at the specified radius
  Double_t fRadius;           // Radius at which the separation is calculated
  Double_t fEtaMin;           // Minimum allowed pair separation in eta

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistance, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistance::Clone() { AliFemtoPairCutRadialDistance* c = new AliFemtoPairCutRadialDistance(*this); return c;}

#endif
