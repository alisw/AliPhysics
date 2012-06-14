#ifndef __AliAnaChargedJetResponseMaker_h__
#define __AliAnaChargedJetResponseMaker_h__
#include "Rtypes.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH2.h"
#include "TProfile.h"
#include "THnSparse.h"

class TGraph;
class TGraphErrors;

//
// Measured spectrum defines minimum and maximum pT on the reconstructed axis of the response matrix. To be set with SetMeasuredSpectrum(TH1D *hPtMeasured);
//

class AliAnaChargedJetResponseMaker {
 public:
  AliAnaChargedJetResponseMaker();
  ~AliAnaChargedJetResponseMaker();

  // kParam        = use parametrization of response 
  // kResiduals    = use response as measured, w/o  statistical error propagation
  // kResidualsErr = use response as measured, with statistical error propagation
  enum ResolutionType {kParam,kResiduals,kResidualsErr}; 

  //Setters
  void SetDebugMode(Bool_t b) {fDebug=b;}

  void SetResolutionType(ResolutionType r) {fResolutionType=r;}

  void SetDeltaPtJetsFunc(TF1 *f1)  {fDeltaPt=f1;}
  void SetDeltaPtJetsHist(TH1D *h1) {fhDeltaPt=h1;}
  void SetNDimensions(Int_t dim)    {fDimensions = dim;}
  void SetMeasuredSpectrum(TH1D *hPtMeasured);

  void SetFlatEfficiency(Double_t eff);
  void SetEfficiency(TGraphErrors *grEff);

  void SetPtMinUnfolded(Double_t ptmin)      {fPtMinUnfolded = ptmin;}
  void SetPtMaxUnfolded(Double_t ptmax)      {fPtMaxUnfolded = ptmax;}
  void SetPtMaxUnfoldedHigh(Double_t ptmaxh) {fPtMaxUnfoldedHigh = ptmaxh;}
  void SetBinWidthFactorUnfolded(Int_t fac)  {fBinWidthFactorUnfolded = fac;}
  void SetSkipBinsUnfolded(Int_t skip)       {fSkipBinsUnfolded=skip;}
  void SetExtraBinsUnfolded(Int_t extra)     {fExtraBinsUnfolded=extra;}
  void SetVariableBinning(Bool_t b, double ptmax) {
    fbVariableBinning=b;
    fPtMaxUnfVarBinning=ptmax;
  }
  void SetCalcErrors(Bool_t b) {fbCalcErrors=b;}

  //Setters for merging fine to normal response matrix
  void SetFineFrac(Int_t i = 10)         {fFineFrac = i;}
  void SetRMMergeWeightFunction(TF1 *f1) {f1MergeFunction = f1;}

  //Getters
  TF1       *GetDeltaPtJetsFunc() {return fDeltaPt;}
  TH1D      *GetDeltaPtJetsHist() {return fhDeltaPt;}
  THnSparse *GetMeasuredSpectrum() {return fPtMeasured;}
  THnSparse *GetEfficiency() {return fEfficiency;}
  THnSparse *GetEfficiencyFine() {return fEfficiencyFine;}
  THnSparse *GetResponseMatrix() {return fResponseMatrix;}
  THnSparse *GetResponseMatrixFine() {return fResponseMatrixFine;}

  //Utility functions
  Double_t InterpolateFast(TGraph *gr, Double_t x);
  Double_t InterpolateFast(TH1 *h, Double_t x);

  TH1D *MultiplyResponseGenerated(TH1 *hGen=0, TH2 *hResponse=0,TH1 *hEfficiency=0,Bool_t bDrawSlices=kFALSE);
  TH1D *MultiplyResponseGenerated(TF1 *fGen, TH2 *hResponse,TH1 *hEfficiency);

  void MakeResponseMatrixJetsFineMerged(Int_t skipBins =0, Int_t binWidthFactor = 2, Int_t extraBins = 0, Bool_t bVariableBinning = kFALSE, Double_t ptmin = 0.);

  void InitializeResponseMatrix();
  void InitializeResponseMatrixFine();

  void InitializeEfficiency();
  void InitializeEfficiencyFine();

  void FillResponseMatrixFineAndMerge();

 protected:
  Bool_t      fDebug;
  ResolutionType fResolutionType;
  TF1        *fDeltaPt;
  TH1D       *fhDeltaPt;
  Int_t       fDimensions; //number of dimensions to unfold (class only prepared for 1 dimension -> 2D response matrix)
  Int_t       fDimRec;
  Int_t       fDimGen;
  Double_t    fPtMin;
  Double_t    fPtMax;
  Int_t       fNbins;
  Double_t   *fBinArrayPtRec;
  THnSparse  *fPtMeasured;
  Double_t    fEffFlat;
  THnSparse  *fEfficiency;
  THnSparse  *fEfficiencyFine;
  THnSparse  *fResponseMatrix;
  THnSparse  *fResponseMatrixFine;
  Double_t    fPtMinUnfolded;          //Minimum pt for unfolded spectrum
  Double_t    fPtMaxUnfolded;          //Maximum pt for unfolded spectrum
  Double_t    fPtMaxUnfoldedHigh;      //Extend last bin of unfolded axis up to fPtMaxUnfoldedHigh
  Int_t       fBinWidthFactorUnfolded; //Unfolded bins x times wider than measured
  Int_t       fSkipBinsUnfolded;       //#unfolded bins to be skipped starting from fPtMinUnfolded
  Int_t       fExtraBinsUnfolded;      //Extra unfolded bins for pTUnf>pTMeas
  Bool_t      fbVariableBinning;       //Unfolded bins 2x narrower for pTUnf<pTMaxUnfVarBinning compared to bin width for pTUnf>pTMaxUnfVarBinning
  Double_t    fPtMaxUnfVarBinning;
  TF1        *f1MergeFunction;
  Int_t       fFineFrac;
  Bool_t      fbCalcErrors;

  ClassDef(AliAnaChargedJetResponseMaker,0);
    
};

#endif
