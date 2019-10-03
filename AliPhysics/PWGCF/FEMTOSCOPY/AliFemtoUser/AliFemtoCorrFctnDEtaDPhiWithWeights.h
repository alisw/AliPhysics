////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDETADPHIWITHWEIGHTS_H
#define ALIFEMTOCORRFCTNDETADPHIWITHWEIGHTS_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDEtaDPhiWithWeights : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kEta=2};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnDEtaDPhiWithWeights(const char* title, TH2D *filter1, TH2D *filter2, const int& aPhiBins, const int& aEtaBins);
  AliFemtoCorrFctnDEtaDPhiWithWeights(const AliFemtoCorrFctnDEtaDPhiWithWeights& aCorrFctn);
  virtual ~AliFemtoCorrFctnDEtaDPhiWithWeights();

  AliFemtoCorrFctnDEtaDPhiWithWeights& operator=(const AliFemtoCorrFctnDEtaDPhiWithWeights& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void SetDoPtAnalysis(int do2d);
  void SetDo4DCorrectionHist(CorrectionType doCorr);

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDEtaDPhiWithWeights(*this); }

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
  TH1D *fPtSumDist;

  TH2D *fYtYtNumerator;
  TH2D *fYtYtDenominator;

  TH2D *fYPtWeightsParticle1;
  TH2D *fYPtWeightsParticle2;


  CorrectionType fIfCorrectionHist;
  THnSparseF *fPtCorrectionsNum;
  THnSparseF *fPtCorrectionsDen;

  THnSparseF *fEtaCorrectionsNum;
  THnSparseF *fEtaCorrectionsDen;

  double fphiL;
  double fphiT;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDEtaDPhiWithWeights, 1)
#endif
};


#endif

