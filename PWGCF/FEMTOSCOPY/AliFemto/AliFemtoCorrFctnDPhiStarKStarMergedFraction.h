////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDPhiStarKStarMergedFraction - correlation function for two //
// particle correlations which uses dPhi* and k* as a function variables,     //
// calculates the fraction of "merged" points                                 //
//                                                                            //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDPHISTARKSTARMERGEDFRACTION_H
#define ALIFEMTOCORRFCTNDPHISTARKSTARMERGEDFRACTION_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoCorrFctnDPhiStarKStarMergedFraction : public AliFemtoCorrFctn {
public:

  AliFemtoCorrFctnDPhiStarKStarMergedFraction(const char* title, Double_t aRadiusMin, Double_t aRadiusMax, Double_t aDistanceMax, Double_t aMergedFractionLimit, Double_t aDEtaMax, const int& aKStarBins, Double_t aKStarRangeLow, Double_t aKStarRangeUp, const Int_t& aPhiStarBins, Double_t aPhiStarRangeLow, Double_t aPhiStarRangeUp);
  AliFemtoCorrFctnDPhiStarKStarMergedFraction(const AliFemtoCorrFctnDPhiStarKStarMergedFraction& aCorrFctn);
  virtual ~AliFemtoCorrFctnDPhiStarKStarMergedFraction();

  AliFemtoCorrFctnDPhiStarKStarMergedFraction& operator=(const AliFemtoCorrFctnDPhiStarKStarMergedFraction& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDPhiStarKStarMergedFraction(*this); }


  void SetRadiusMin(double minrad);
  void SetRadiusMax(double maxrad);
  void SetDistanceMax(double maxdist);
  void SetMergedFractionLimit(double frac);
  void SetDEtaMax(double deta);
  void SetMagneticFieldSign(int magsign);

private:

  TH2D *fDPhiStarKStarMergedNumerator;              // Numerator of dPhi* k* function for pairs which are "merged"
  TH2D *fDPhiStarKStarTotalNumerator;               // Numerator of dPhi* k* function for all pairs
  TH2D *fDPhiStarKStarMergedDenominator;            // Denominator of dPhi* k* function for pairs which are "merged"
  TH2D *fDPhiStarKStarTotalDenominator;             // Denominator of dPhi* k* function for all pairs

  double fDPhiStarRangeLow;           // Lower range of dPhi*
  double fDPhiStarRangeUp;            // Upper range of dPhi*

  double fKStarRangeLow;               // Lower range of k*
  double fKStarRangeUp;                // Upper range of k*

  Double_t fDistanceMax;            // Maximum distance where "merging" can occur [m]
  Double_t fMergedFractionLimit;    // Fraction of merged points, above which the pair is removed (between 0 and 1)
  Double_t fDEtaMax;                // Maximum value of dEta, where "merging" can occur

  Double_t fRadiusMin;              // Minimum radius at which the pair separation is calculated [m]
  Double_t fRadiusMax;              // Maximum radius at which the pair separation is calculated [m]

  Int_t fMagSign;                    // Magnetic field sign

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnDPhiStarKStarMergedFraction, 1);
  /// \endcond
#endif
};


#endif

