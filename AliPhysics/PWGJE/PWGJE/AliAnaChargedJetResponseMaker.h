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
  AliAnaChargedJetResponseMaker(const AliAnaChargedJetResponseMaker& obj); // copy constructor
  AliAnaChargedJetResponseMaker& operator=(const AliAnaChargedJetResponseMaker& other); // assignment
  virtual ~AliAnaChargedJetResponseMaker() {;}

  // kParam        = use parametrization of response 
  // kResiduals    = use response as measured, w/o  statistical error propagation
  // kResidualsErr = use response as measured, with statistical error propagation
  enum ResolutionType {kParam,kResiduals,kResidualsErr}; 

  //Setters
  virtual void SetDebugMode(Bool_t b) {fDebug=b;}

  virtual void SetResolutionType(ResolutionType r) {fResolutionType=r;}

  virtual void SetDeltaPtJetsFunc(TF1 *f1)  {fDeltaPt=f1;}
  virtual void SetDeltaPtJetsHist(TH1D *h1) {fhDeltaPt=h1;}
  virtual void SetNDimensions(Int_t dim)    {fDimensions = dim;}
  virtual void SetMeasuredSpectrum(TH1D *hPtMeasured);
  virtual void SetMeasuredSpectrumTruncated(TH1D *h1) {fh1MeasuredTruncated=h1;}

  virtual void SetDetectorResponse(TH2 *h2) {fh2DetectorResponse=(TH2D*)h2;}
  virtual void SetDetectorEfficiency(TH1D *h1) {fhEfficiencyDet=h1;}

  virtual void SetFlatEfficiency(Double_t eff);
  virtual void SetEfficiency(TGraphErrors *grEff);

  virtual void SetPtMinUnfolded(Double_t ptmin)      {fPtMinUnfolded = ptmin;}
  virtual void SetPtMaxUnfolded(Double_t ptmax)      {fPtMaxUnfolded = ptmax;}
  virtual void SetPtMaxUnfoldedHigh(Double_t ptmaxh) {fPtMaxUnfoldedHigh = ptmaxh;}
  virtual void SetBinWidthFactorUnfolded(Int_t fac)  {fBinWidthFactorUnfolded = fac;}
  virtual void SetSkipBinsUnfolded(Int_t skip)       {fSkipBinsUnfolded=skip;}
  virtual void SetExtraBinsUnfolded(Int_t extra)     {fExtraBinsUnfolded=extra;}
  virtual void SetVariableBinning(Bool_t b, double ptmax) {
    fbVariableBinning=b;
    fPtMaxUnfVarBinning=ptmax;
  }
  virtual void SetCalcErrors(Bool_t b) {fbCalcErrors=b;}

  //Setters for merging fine to normal response matrix
  virtual void SetFineFrac(Int_t i = 10)         {fFineFrac = i;}
  virtual void SetRMMergeWeightFunction(TF1 *f1) {f1MergeFunction = f1;}

  //Getters
  virtual TF1       *GetDeltaPtJetsFunc() {return fDeltaPt;}
  virtual TH1D      *GetDeltaPtJetsHist() {return fhDeltaPt;}
  virtual THnSparse *GetMeasuredSpectrum() {return fPtMeasured;}
  virtual THnSparse *GetEfficiency() {return fEfficiency;}
  virtual THnSparse *GetEfficiencyFine() {return fEfficiencyFine;}
  virtual THnSparse *GetResponseMatrix() {return fResponseMatrix;}
  virtual THnSparse *GetResponseMatrixFine() {return fResponseMatrixFine;}

  virtual TH2D      *GetDetectorResponseRebin() {return fh2DetectorResponseRebin;}

  virtual Bool_t CheckInputForCombinedResponse();

  virtual TH2D      *GetResponseCombinedFineFull() {return fh2ResponseMatrixCombinedFineFull;}
  virtual TH2D      *GetResponseCombinedFull() {return fh2ResponseMatrixCombinedFull;}

  virtual TH2D      *GetResponseCombined() {return fh2ResponseMatrixCombined;}
  virtual TH1D      *GetEfficiencyCombined() {return fhEfficiencyCombined;}

  static Double_t   GetBetaPerDOFValue(Int_t betaColl = 0, Int_t betaOpt= 0);

  //Utility functions
  virtual Double_t InterpolateFast(TGraph *gr, Double_t x);
  virtual Double_t InterpolateFast(TH1 *h, Double_t x);

  virtual TH1D *MultiplyResponseGenerated(TH1 *hGen=0, TH2 *hResponse=0,TH1 *hEfficiency=0,Bool_t bDrawSlices=kFALSE);
  virtual TH1D *MultiplyResponseGenerated(TF1 *fGen, TH2 *hResponse,TH1 *hEfficiency);

  virtual void MakeResponseMatrixCombined(Int_t skipBins =0, Int_t binWidthFactor = 2, Int_t extraBins = 0, Bool_t bVariableBinning = kFALSE, Double_t ptmin = 0.);

  virtual void MakeResponseMatrixJetsFineMerged(Int_t skipBins =0, Int_t binWidthFactor = 2, Int_t extraBins = 0, Bool_t bVariableBinning = kFALSE, Double_t ptmin = 0.);

  virtual void InitializeResponseMatrix();
  virtual void InitializeResponseMatrixFine();

  virtual void InitializeEfficiency();
  virtual void InitializeEfficiencyFine();

  virtual void FillResponseMatrixFineAndMerge();

  virtual TH2* MakeResponseMatrixRebin(TH2 *hRMFine = 0, TH2 *hRM = 0, Bool_t useFunctionWeight = kFALSE);

  virtual TH2* CreateTruncated2DHisto(TH2 *h2=0, Double_t xmin=-1, Double_t xmax=-1, Double_t ymin=-1, Double_t ymax=-1);
  virtual TH2* TruncateAxisRangeResponseMatrix(TH2 *hRMOrig=0,  Double_t xmin=-1, Double_t xmax=-1, Double_t ymin=-1, Double_t ymax=-1);

  virtual TH2* MultiplityResponseMatrices(TH2 *h2RMDeltaPt, TH2 *h2RMDetector);

  virtual TH2* GetTransposeResponsMatrix(TH2 *h2RM);

  virtual TH2* NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior);

 protected:
  Bool_t      fDebug;
  ResolutionType fResolutionType;
  TF1        *fDeltaPt;
  TH1D       *fhDeltaPt;
  TH1D       *fh1MeasuredTruncated;
  TH2D       *fh2DetectorResponse;
  TH2D       *fh2DetectorResponseRebin;
  TH1D       *fhEfficiencyDet;
  TH2D       *fh2ResponseMatrixCombinedFineFull;
  TH2D       *fh2ResponseMatrixCombinedFull;
  TH2D       *fh2ResponseMatrixCombined;
  TH1D       *fhEfficiencyCombined;
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
