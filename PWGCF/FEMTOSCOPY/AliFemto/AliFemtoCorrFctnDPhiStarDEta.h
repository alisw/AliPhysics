////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDPhiStarDEta - correlation function for two particle       //
// correlations which uses dPhi* and dEta as a function variables.            //
//                                                                            //
// Authors: Przemyslaw Karczmarczyk przemyslaw.karczmarczyk@cern.ch           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDPHISTARDETA_H
#define ALIFEMTOCORRFCTNDPHISTARDETA_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoCorrFctnDPhiStarDEta : public AliFemtoCorrFctn {
public:

  AliFemtoCorrFctnDPhiStarDEta(const char* title, double radius, const int& aEtaBins, double aEtaRangeLow, double aEtaRangeUp, const int& aPhiStarBins, double aPhiStarRangeLow, double aPhiStarRangeUp);
  AliFemtoCorrFctnDPhiStarDEta(const AliFemtoCorrFctnDPhiStarDEta& aCorrFctn);
  virtual ~AliFemtoCorrFctnDPhiStarDEta();

  AliFemtoCorrFctnDPhiStarDEta& operator=(const AliFemtoCorrFctnDPhiStarDEta& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDPhiStarDEta(*this); }

  void SetMinimumRadius(double minrad);
  void SetMagneticFieldSign(int magsign);

private:

  TH2D *fDPhiStarDEtaNumerator;      // Numerator of dPhiStar dEta function
  TH2D *fDPhiStarDEtaDenominator;    // Denominator of dPhiStar dEta function

  double fEtaRangeLow;               // Lower range of Eta
  double fEtaRangeUp;                // Upper range of PhiStar

  double fPhiStarRangeLow;           // Lower range of PhiStar
  double fPhiStarRangeUp;            // Upper range of PhiStar

  Double_t fMinRad;                  // Set minimum radial distance
  Int_t fMagSign;                    // Magnetic field sign

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnDPhiStarDEta, 1);
  /// \endcond
#endif
};


#endif

