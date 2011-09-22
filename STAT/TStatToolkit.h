#ifndef TSTATTOOLKIT_H
#define TSTATTOOLKIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// some utilities which do net exist in the standard ROOT
//
 
#include "TObject.h"
#include "TVectorD.h"
#include "TMatrixD.h"
//#include "TGraph2D.h"
//#include "TGraph.h"

class TH1F;
class TH1;
class TH3;
class TString;
class TTree;
class TGraph;
class TGraph2D;
 
class TStatToolkit : public TObject
{
 public:
  TStatToolkit();
  virtual ~TStatToolkit();
  //
  //
  //
  static void    EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static void    EvaluateUniExternal(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh, Float_t externalfactor=1);
  static Int_t  Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    
  //
  // HISTOGRAMS TOOLS
  //
  static  void TruncatedMean(const TH1 * his, TVectorD *param, Float_t down=0, Float_t up=1.0, Bool_t verbose=kFALSE);
  static void LTM(TH1F * his, TVectorD *param=0 , Float_t fraction=1,  Bool_t verbose=kFALSE);
  static Double_t  FitGaus(TH1* his, TVectorD *param=0, TMatrixD *matrix=0, Float_t xmin=0, Float_t xmax=0,  Bool_t verbose=kFALSE);
  static Double_t  FitGaus(Float_t *arr, Int_t nBins, Float_t xMin, Float_t xMax, TVectorD *param=0, TMatrixD *matrix=0, Bool_t verbose=kFALSE);
  static Float_t  GetCOG(const Short_t *arr, Int_t nBins, Float_t xMin, Float_t xMax, Float_t *rms=0, Float_t *sum=0);

  static TGraph2D *  MakeStat2D(TH3 * his, Int_t delta0, Int_t delta1, Int_t type);
  static TGraph *  MakeStat1D(TH3 * his, Int_t delta1, Int_t type);
  //
  // Graph tools
  //
  static TGraph * MakeGraphSparse(TTree * tree, const char * expr="Entry", const char * cut="1");
  //
  // Fitting function
  //
  static TString* FitPlane(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000, Bool_t fix0=kFALSE);
  static TString* FitPlaneFixed(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000);
  //
  //Linear fitter helper function
  //
  static TString* FitPlaneConstrain(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000, Double_t constrain=-1);
  static Int_t GetFitIndex(const TString fString, const TString subString);
 static TString FilterFit(const TString &input, const TString filter, TVectorD &vec, TMatrixD &covar);
 static void Update1D(Double_t delta, Double_t sigma, Int_t s1, TMatrixD &param, TMatrixD &covar);
  static void   Constrain1D(const TString &input, const TString filter, TVectorD &param, TMatrixD & covar, Double_t mean, Double_t sigma);
  static TString  MakeFitString(const TString &input, const TVectorD &param, const TMatrixD & covar, Bool_t verbose=kFALSE);

  //
  // TestFunctions:
  //
 static  void TestGausFit(Int_t nhistos=5000);

 ClassDef(TStatToolkit,0) // Various mathematical tools for physics analysis - which are not included in ROOT TMath
 
};
#endif
