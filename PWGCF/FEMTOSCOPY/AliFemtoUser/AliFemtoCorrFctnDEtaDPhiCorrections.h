////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhiCorrections - A correlation function that analyzes //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHICORRECTIONS_H
#define ALIFEMTOCORRFCTNDETADPHICORRECTIONS_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TFile.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiCorrections : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kEta=2};
  enum ParticleType {kNoCorrection=0, kPion=1, kKaon=2, kProton=3};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnDEtaDPhiCorrections(char* title, const int& aPhiBins, const int& aEtaBins);
  AliFemtoCorrFctnDEtaDPhiCorrections(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiCorrections();

  AliFemtoCorrFctnDEtaDPhiCorrections& operator=(const AliFemtoCorrFctnDEtaDPhiCorrections& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void SetDoPtAnalysis(int do2d);
  void SetDoCorrections(bool doCorr);
  void SetDoCorrectionsHist(CorrectionType doCorr);
  double CalculateCorrectionWeight(double pT1, double pT2);
  void LoadCorrectionTabFromFile(const char *pTtab, const char *corrTab);
  
  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fDPhiDEtaNumerator;          // Numerator of dEta dPhi function
  TH2D *fDPhiDEtaDenominator;        // Denominator of dEta dPhi function

  TH1D *fDPhiNumerator;              // Numerator of dPhi correlation
  TH1D *fDPhiDenominator;            // Denominator of dPhi correlation

  TH1D *fDCosNumerator;              // Numerator of colinearity correlation
  TH1D *fDCosDenominator;            // Denominator of colinearity correlation

  int   fDoPtAnalysis;               // set to 1 to do 2D Pt analysis

  TH2D *fDPhiPtNumerator;            // Numerator of dPhi correlation vs. Pt min
  TH2D *fDPhiPtDenominator;          // Denominator of dPhi correlation vs. Pt min

  TH2D *fDCosPtNumerator;            // Numerator of colinearity correlation vs. Pt min
  TH2D *fDCosPtDenominator;          // Denominator of colinearity correlation vs. Pt min

  TH1D *fPhi;
  TH1D *fEta;

  TH2D *fYtYtNumerator;
  TH2D *fYtYtDenominator; 

  CorrectionType fIfCorrectionHist;
  bool fIfCorrection;
  THnSparseF *fPtCorrectionsNum;
  THnSparseF *fPtCorrectionsDen;

  THnSparseF *fEtaCorrectionsNum;
  THnSparseF *fEtaCorrectionsDen;

  double* fCorrFactorTab;
  double* fpTab;
  ParticleType fPartType; // particle type for calculations of correction factor

  double fphiL;
  double fphiT;


#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiCorrections, 1)
#endif
};


#endif

