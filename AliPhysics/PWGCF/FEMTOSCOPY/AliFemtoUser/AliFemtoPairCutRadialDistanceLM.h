/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceLM - a pair cut which checks                   //
// for some pair qualities that attempt to identify slit/doubly                //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////
/********************************************************************************
 *
 * Author: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *         Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch
 *         Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch
 *         Jorge Mercado, University of Heidelberg, jmercado@cern.ch
 *
 ********************************************************************************/

#ifndef ALIFEMTOPAIRCUTRADIALDISTANCELM_H
#define ALIFEMTOPAIRCUTRADIALDISTANCELM_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

class AliFemtoPairCutRadialDistanceLM : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutRadialDistanceLM();
  AliFemtoPairCutRadialDistanceLM(const AliFemtoPairCutRadialDistanceLM& c);
  virtual ~AliFemtoPairCutRadialDistanceLM();

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
  ClassDef(AliFemtoPairCutRadialDistanceLM, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistanceLM::Clone() { AliFemtoPairCutRadialDistanceLM* c = new AliFemtoPairCutRadialDistanceLM(*this); return c;}

#endif
