/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistance - a pair cut which checks                     //
// for some pair qualities that attempt to identify split/doubly               //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC, adjusted to perform asymmetric cut for pair     //
//                                                                             //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOPAIRCUTRADIALDISTANCEASYMMETRIC_H
#define ALIFEMTOPAIRCUTRADIALDISTANCEASYMMETRIC_H

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutRadialDistanceAsymmetric : public AliFemtoPairCutAntiGamma {
 public:
  AliFemtoPairCutRadialDistanceAsymmetric(Double_t aRadiusMin, Double_t aRadiusMax, Bool_t aCalculateRadiusRange, Double_t aDEtaRangeLow, Double_t aDEtaRangeUp, Double_t aDPhiStarRangeLow, Double_t aDPhiStarRangeUp);
  AliFemtoPairCutRadialDistanceAsymmetric(const AliFemtoPairCutRadialDistanceAsymmetric& cPairCut);
  virtual ~AliFemtoPairCutRadialDistanceAsymmetric();
  AliFemtoPairCutRadialDistanceAsymmetric& operator=(const AliFemtoPairCutRadialDistanceAsymmetric& cPairCut);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();

  void SetDPhiStarRangeLow(double minphistar);
  void SetDPhiStarRangeUp(double maxphistar);
  void SetDEtaRangeLow(double mineta);
  void SetDEtaRangeUp(double maxeta);
  void SetMagneticFieldSign(int magsign);
  void SetRadiusMin(double radmin);
  void SetRadiusMax(double radmax);
  void SetCalculateRadiusRange(bool calculaterange);

 protected:
  Double_t fDPhiStarRangeLow;       // Lower range of pair separation cut in DPhiStar [rad]
  Double_t fDPhiStarRangeUp;        // Upper range of pair separation cut in DPhiStar [rad]
  Double_t fDEtaRangeLow;           // Lower range of pair separation cut in DEta
  Double_t fDEtaRangeUp;            // Upper range of pair separation cut in DEta
  Int_t fMagSign;                   // Magnetic field sign (+1/-1)
  Double_t fRadiusMin;              // Minimum radius at which the pair separation is calculated [m]
  Double_t fRadiusMax;              // Maximum radius at which the pair separation is calculated [m]
  Bool_t fCalculateRadiusRange;     // If true: iterate through all radii in range (fRadiusMin, fRadiusMax), if false: perform cut only for fRadiusMin

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutRadialDistanceAsymmetric, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutRadialDistanceAsymmetric::Clone() { AliFemtoPairCutRadialDistanceAsymmetric* cPairCut = new AliFemtoPairCutRadialDistanceAsymmetric(*this); return cPairCut;}

#endif
