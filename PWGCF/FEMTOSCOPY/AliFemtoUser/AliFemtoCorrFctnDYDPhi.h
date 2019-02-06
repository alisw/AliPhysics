////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDYDPhi - A correlation function that analyzes              //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and rapidity (y) difference                                                //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch,                                  //
//          Piotr Modzelewski Piotr.Mateusz.Modzelewski@cern.ch               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNDYDPHI_H
#define ALIFEMTOCORRFCTNDYDPHI_H

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnDYDPhi : public AliFemtoCorrFctn {
public:
  enum CorrectionType {kNone=0, kPt=1, kY=2};
  typedef enum CorrectionType ReadCorrectionType;

  AliFemtoCorrFctnDYDPhi(const char* title, const int& aPhiBins, const int& aYBins, const double& mass);
  AliFemtoCorrFctnDYDPhi(const AliFemtoCorrFctnDYDPhi& aCorrFctn);
  virtual ~AliFemtoCorrFctnDYDPhi();

  AliFemtoCorrFctnDYDPhi& operator=(const AliFemtoCorrFctnDYDPhi& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void SetDoPtAnalysis(int do2d);
  void SetDo4DCorrectionHist(CorrectionType doCorr);

  void WriteHistos();
  virtual TList* GetOutputList();

  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDYDPhi(*this); }

private:

  TH2D *fDPhiDYNumerator;            // Numerator of dY dPhi function
  TH2D *fDPhiDYDenominator;          // Denominator of dY dPhi function

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
  TH1D *fY;
  TH1D *fPtSumDist;

  TH2D *fYtYtNumerator;
  TH2D *fYtYtDenominator;

  CorrectionType fIfCorrectionHist;
  THnSparseF *fPtCorrectionsNum;
  THnSparseF *fPtCorrectionsDen;

  THnSparseF *fYCorrectionsNum;
  THnSparseF *fYCorrectionsDen;

  double fphiL;
  double fphiT;

  double fMass;

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnDYDPhi, 1)
#endif
};


#endif

