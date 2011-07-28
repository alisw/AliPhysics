/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistance - a pair cut which checks                     //
// for some pair qualities that attempt to identify slit/doubly                //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////
/********************************************************************************
 *
 * Authors: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *          Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch
 *          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch
 *
 ********************************************************************************/



#ifndef ALIFEMTOPAIRCUTRADIALDISTANCE_H
#define ALIFEMTOPAIRCUTRADIALDISTANCE_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

class AliFemtoPairCutRadialDistance : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutRadialDistance();
  AliFemtoPairCutRadialDistance(const AliFemtoPairCutRadialDistance& c);
  virtual ~AliFemtoPairCutRadialDistance();

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetPhiStarDifferenceMinimum(double dtpc);
  void SetEtaDifferenceMinimum(double etpc);
  void SetMinimumRadius(double minrad);
  void SetMagneticFieldSign(int magsign);


 protected:
  Double_t fDPhiStarMin;          // Minimum allowed pair separation //at the specified radius
  //Double_t fRadius;           // Radius at which the separation is calculated
  Double_t fEtaMin;           // Minimum allowed pair separation in eta
  Double_t fMinRad;
  Int_t fMagSign;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistance, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistance::Clone() { AliFemtoPairCutRadialDistance* c = new AliFemtoPairCutRadialDistance(*this); return c;}

#endif
