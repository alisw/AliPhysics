/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceV0 - a pair cut which checks                   //
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



#ifndef ALIFEMTOPAIRCUTRADIALDISTANCEV0_H
#define ALIFEMTOPAIRCUTRADIALDISTANCEV0_H

// do I need these lines ?
//#ifndef StMaker_H
//#include "StMaker.h"
//#endif

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutRadialDistanceV0 : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutRadialDistanceV0();
  AliFemtoPairCutRadialDistanceV0(const AliFemtoPairCutRadialDistanceV0& c);
  virtual ~AliFemtoPairCutRadialDistanceV0();
  AliFemtoPairCutRadialDistanceV0& operator=(const AliFemtoPairCutRadialDistanceV0& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();

  void SetPhiStarDifferenceMinimumPrimDaughter1(double dtpc);
  void SetPhiStarDifferenceMinimumPrimDaughter2(double dtpc);
  void SetPhiStarDifferenceMinimumDaughter1Daughter1(double dtpc);
  void SetPhiStarDifferenceMinimumDaughter1Daughter2(double dtpc);
  void SetPhiStarDifferenceMinimumDaughter2Daughter1(double dtpc);
  void SetPhiStarDifferenceMinimumDaughter2Daughter2(double dtpc);

  void SetEtaDifferenceMinimumPrimDaughter1(double etpc);
  void SetEtaDifferenceMinimumPrimDaughter2(double etpc);
  void SetEtaDifferenceMinimumDaughter1Daughter2(double etpc);
  void SetEtaDifferenceMinimumDaughter2Daughter1(double etpc);
  void SetEtaDifferenceMinimumDaughter2Daughter2(double etpc);

  void SetMinimumRadiusPrimDaughter1(double minrad);
  void SetMinimumRadiusPrimDaughter2(double minrad);
  void SetMinimumRadiusDaughter1Daughter1(double minrad);
  void SetMinimumRadiusDaughter1Daughter2(double minrad);
  void SetMinimumRadiusDaughter2Daughter1(double minrad);
  void SetMinimumRadiusDaughter2Daughter2(double minrad);

  void SetMagneticFieldSign(int magsign);


  void SetPhiStarMinPrimDaughter1(Bool_t);
  void SetPhiStarMinPrimDaughter2(Bool_t);
  void SetPhiStarMinDaughter1Daughter1(Bool_t);
  void SetPhiStarMinDaughter1Daughter2(Bool_t);
  void SetPhiStarMinDaughter2Daughter1(Bool_t);
  void SetPhiStarMinDaughter2Daughter2(Bool_t);


 protected:
  Double_t fDPhiStarMinPrimDaughter1;          // Minimum allowed pair separation //at the specified radius
  Double_t fDPhiStarMinPrimDaughter2;          // Minimum allowed pair separation //at the specified radius
  Double_t fDPhiStarMinDaughter1Daughter1;          // Minimum allowed pair separation //at the specified radius
  Double_t fDPhiStarMinDaughter1Daughter2;          // Minimum allowed pair separation //at the specified radius
  Double_t fDPhiStarMinDaughter2Daughter1;          // Minimum allowed pair separation //at the specified radius
  Double_t fDPhiStarMinDaughter2Daughter2;          // Minimum allowed pair separation //at the specified radius
  //Double_t fRadius;           // Radius at which the separation is calculated
  Double_t fEtaMinPrimDaughter1;           // Minimum allowed pair separation in eta
  Double_t fEtaMinPrimDaughter2;           // Minimum allowed pair separation in eta
  Double_t fEtaMinDaughter1Daughter1;           // Minimum allowed pair separation in eta
  Double_t fEtaMinDaughter1Daughter2;           // Minimum allowed pair separation in eta
  Double_t fEtaMinDaughter2Daughter1;           // Minimum allowed pair separation in eta
  Double_t fEtaMinDaughter2Daughter2;           // Minimum allowed pair separation in eta
  Double_t fMinRadPrimDaughter1;
  Double_t fMinRadPrimDaughter2;
  Double_t fMinRadDaughter1Daughter1;
  Double_t fMinRadDaughter1Daughter2;
  Double_t fMinRadDaughter2Daughter1;
  Double_t fMinRadDaughter2Daughter2;
  Int_t fMagSign;

  Bool_t fPhistarminPrimDaughter1;
  Bool_t fPhistarminPrimDaughter2;
  Bool_t fPhistarminDaughter1Daughter1;
  Bool_t fPhistarminDaughter1Daughter2;
  Bool_t fPhistarminDaughter2Daughter1;
  Bool_t fPhistarminDaughter2Daughter2;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistanceV0, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistanceV0::Clone() { AliFemtoPairCutRadialDistanceV0* c = new AliFemtoPairCutRadialDistanceV0(*this); return c;}

#endif
