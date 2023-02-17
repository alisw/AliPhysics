/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceDR - a pair cut which checks                     //
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



#ifndef ALIFEMTOPAIRCUTRADIALDISTANCEDR_H
#define ALIFEMTOPAIRCUTRADIALDISTANCEDR_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutRadialDistanceDR : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutRadialDistanceDR();
  AliFemtoPairCutRadialDistanceDR(const AliFemtoPairCutRadialDistanceDR& c);
  virtual ~AliFemtoPairCutRadialDistanceDR();
  AliFemtoPairCutRadialDistanceDR& operator=(const AliFemtoPairCutRadialDistanceDR& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetPhiStarDifferenceMinimum(double dtpc);
  void SetEtaDifferenceMinimum(double etpc);
  void SetMinimumRadius(double minrad);
  void SetMaximumRadius(double maxrad);
  void SetMagneticFieldSign(int magsign);
  void SetMagneticFieldValue(double magval);
  void SetPhiStarMin(Bool_t);


 protected:
  Double_t fDPhiStarMin;          // Minimum allowed pair separation //at the specified radius
  //Double_t fRadius;           // Radius at which the separation is calculated
  Double_t fEtaMin;           // Minimum allowed pair separation in eta
  Double_t fMinRad;
  Double_t fMaxRad;
  Int_t fMagSign;
  Double_t fMagFieldVal; 			// Magnetic field value (default 0.5)
  Bool_t fPhistarmin;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistanceDR, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistanceDR::Clone() { AliFemtoPairCutRadialDistanceDR* c = new AliFemtoPairCutRadialDistanceDR(*this); return c;}

#endif
