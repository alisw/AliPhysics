////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELCORRFCTNDETADPHI_H
#define ALIFEMTOMODELCORRFCTNDETADPHI_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelCorrFctn.h"

class AliFemtoModelCorrFctnDEtaDPhi : public AliFemtoModelCorrFctn {
public:
  AliFemtoModelCorrFctnDEtaDPhi(const char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoModelCorrFctnDEtaDPhi(const AliFemtoModelCorrFctnDEtaDPhi& aCorrFctn);
  virtual ~AliFemtoModelCorrFctnDEtaDPhi();

  AliFemtoModelCorrFctnDEtaDPhi& operator=(const AliFemtoModelCorrFctnDEtaDPhi& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoModelCorrFctn* Clone() const { return new AliFemtoModelCorrFctnDEtaDPhi(*this); }

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

  double fphiL;
  double fphiT;

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctnDEtaDPhi, 1)
#endif
};


#endif

