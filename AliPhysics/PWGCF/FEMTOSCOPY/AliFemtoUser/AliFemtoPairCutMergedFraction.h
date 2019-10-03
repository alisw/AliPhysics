/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutMergedFraction - a pair cut which checks                     //
// for some pair qualities that attempt to identify split/doubly               //
// reconstructed tracks, calculating the fraction of "merged" points (number   //
// of points in selected range of radius from primary vertex where tracks are  //
// closer than selected distance divided by total number of calculated points) //
//                                                                             //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOPAIRCUTMERGEDFRACTION_H
#define ALIFEMTOPAIRCUTMERGEDFRACTION_H

#include "AliFemtoPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutMergedFraction : public AliFemtoPairCutAntiGamma {
 public:
  AliFemtoPairCutMergedFraction(Double_t aDistanceMax, Double_t aMergedFractionLimit, Double_t aDEtaMax, Double_t aRadiusMin, Double_t aRadiusMax);
  AliFemtoPairCutMergedFraction(const AliFemtoPairCutMergedFraction& cPairCut);
  virtual ~AliFemtoPairCutMergedFraction();
  AliFemtoPairCutMergedFraction& operator=(const AliFemtoPairCutMergedFraction& cPairCut);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCutMergedFraction* Clone() const;

  void SetDistanceMax(double maxdistance);
  void SetMergedFractionLimit(double fractionlimit);
  void SetDEtaMax(double maxeta);
  void SetMagneticFieldSign(int magsign);
  void SetRadiusMin(double radmin);
  void SetRadiusMax(double radmax);
  void SetMagneticFieldValue(double magval);
  void SetMergedFractionDataType(AliFemtoDataType datatype);

 protected:
  Double_t fDistanceMax;            // Maximum distance where "merging" can occur [m]
  Double_t fMergedFractionLimit;    // Fraction of merged points, above which the pair is removed (between 0 and 1)
  Double_t fDEtaMax;                // Maximum value of dEta, where "merging" can occur
  Double_t fRadiusMin;              // Minimum radius at which the pair separation is calculated [m]
  Double_t fRadiusMax;              // Maximum radius at which the pair separation is calculated [m]
  Int_t fMagSign;                   // Magnetic field sign (+1/-1)
  Double_t fMagFieldVal; 			// Magnetic field value (default 0.5)
  AliFemtoDataType fMergedFractionDataType;       // Use ESD/AOD/Kinematics

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutMergedFraction, 0)
#endif
};

inline AliFemtoPairCutMergedFraction* AliFemtoPairCutMergedFraction::Clone() const
{
  return new AliFemtoPairCutMergedFraction(*this);
}

#endif
