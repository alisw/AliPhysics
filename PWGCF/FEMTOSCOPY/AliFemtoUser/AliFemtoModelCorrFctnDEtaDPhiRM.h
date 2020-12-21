////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelCorrFctnDEtaDPhiRM - A correlation function that analyzes       //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference, contains histograms for MInv distr                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch, rmaselek@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELCORRFCTNDETADPHIRM_H
#define ALIFEMTOMODELCORRFCTNDETADPHIRM_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnDEtaDPhiRM : public AliFemtoModelCorrFctn {
public:
  AliFemtoModelCorrFctnDEtaDPhiRM(const char* title, const int& aPhiBins, const int& aEtaBins, const double m1, const double m2);
  AliFemtoModelCorrFctnDEtaDPhiRM(const AliFemtoModelCorrFctnDEtaDPhiRM& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnDEtaDPhiRM();

  AliFemtoModelCorrFctnDEtaDPhiRM& operator=(const AliFemtoModelCorrFctnDEtaDPhiRM& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:

  TH2D *fDPhiDEtaNumeratorTrue;      // Numerator of dEta dPhi true function
  TH2D *fDPhiDEtaNumeratorFake;      // Numerator of dEta dPhi fake function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function

  TH2D *fDPhiDEtaColNumerator;       // Numerator of colinear dEta dPhi function
  TH2D *fDPhiDEtaColDenominator;     // Denominator of colinear dEta dPhi function

  TH1D *fDPhiNumeratorTrue;          // Numerator of dPhi true correlation
  TH1D *fDPhiNumeratorFake;          // Numerator of dPhi fake correlation
  TH1D *fDPhiDenominator;            // Denominator of dPhi correlation

  TH1D *fDCosNumeratorTrue;           // Numerator of colinearity true correlation
  TH1D *fDCosNumeratorFake;           // Numerator of colinearity fake correlation
  TH1D *fDCosDenominator;            // Denominator of colinearity correlation

  TH2D *fDPhiPtNumerator;            // Numerator of dPhi correlation vs. Pt min
  TH2D *fDPhiPtDenominator;          // Denominator of dPhi correlation vs. Pt min

  TH2D *fDCosPtNumerator;            // Numerator of colinearity correlation vs. Pt min
  TH2D *fDCosPtDenominator;          // Denominator of colinearity correlation vs. Pt min
  //ADDED BY RAFAL MASELEK:
  TH1D *fPtSumDist;
  TH1D *fInvMassDist;

  double fM1; //mass of particle1
  double fM2; //mass of particle2

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnDEtaDPhiRM, 1)
#endif
};


#endif

