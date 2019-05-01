/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceKKdistT - a pair cut which checks                     //
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



#ifndef AliFemtoPairCutRadialDistanceKKdist_H
#define AliFemtoPairCutRadialDistanceKKdist_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutRadialDistanceKKdist : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutRadialDistanceKKdist();
  AliFemtoPairCutRadialDistanceKKdist(const AliFemtoPairCutRadialDistanceKKdist& c);
  virtual ~AliFemtoPairCutRadialDistanceKKdist();
  AliFemtoPairCutRadialDistanceKKdist& operator=(const AliFemtoPairCutRadialDistanceKKdist& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
 // void SetPhiStarDifferenceMinimum(double dtpc);
 // void SetEtaDifferenceMinimum(double etpc);
  void SetAverageSeparation(double minrad);
  //void SetMaximumRadius(double maxrad);
  //void SetMagneticFieldSign(int magsign);
 // void SetPhiStarMin(Bool_t);


 protected:
  //Double_t fDPhiStarMin;          // Minimum allowed pair separation //at the specified radius
  //Double_t fRadius;           // Radius at which the separation is calculated
  //Double_t fEtaMin;           // Minimum allowed pair separation in eta
  Double_t fMinRad;
  //Double_t fMaxRad;
  //Int_t fMagSign;
 // Bool_t fPhistarmin;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistanceKKdist, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistanceKKdist::Clone() { AliFemtoPairCutRadialDistanceKKdist* c = new AliFemtoPairCutRadialDistanceKKdist(*this); return c;}

#endif
