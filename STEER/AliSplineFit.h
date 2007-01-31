#ifndef ALISPLINEFIT_H
#define ALISPLINEFIT_H
/* Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TGraph.h"
#include "TGraphSmooth.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TLinearFitter.h"
#include "TDecompSVD.h"
#include "TDecompSparse.h"
#include "TMatrixDSparse.h"
#include "TF1.h"
#include "TH1F.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "TTreeStream.h"

class AliSplineFit : public TObject {
 public:
  AliSplineFit();
  AliSplineFit(const AliSplineFit&);
  ~AliSplineFit();
  AliSplineFit& operator=(const AliSplineFit&);
  Double_t   Eval(Double_t x, Int_t deriv=0) const;
  void       InitKnots(TGraph * graph, Int_t min, Int_t iter, Double_t maxDelta);
  void       MakeKnots0(TGraph * graph, Double_t maxdelta, Int_t minpoints);
  void       SplineFit(Int_t nder);
  void       MakeSmooth(TGraph * graph, Float_t ratio, char * type);
  void       Update(TSpline3 *spline, Int_t nknots);
  Int_t      GetKnots() const {return fN;} 
  Double_t*  GetX() const {return fX;}
  Double_t*  GetY0() const {return fY0;}
  Double_t*  GetY1() const {return fY1;}
  //
  // Test functions
  //
  TGraph * MakeGraph(Double_t xmin, Double_t xmax, Int_t npoints, Int_t deriv=0) const ;
  TGraph * MakeDiff(TGraph * graph) const ;
  TH1F   * MakeDiffHisto(TGraph * graph) const;
  //
  static void Test(Int_t npoints=2000, Int_t ntracks=100, Float_t snoise=0.05);
  //
  static TGraph * GenerGraph(Int_t npoints, Double_t fraction, Double_t s1, Double_t s2, Double_t s3, Int_t der=0);
  static TGraph * GenerNoise(TGraph * graph0, Double_t s0);

 protected:

  //
  // working parameters for spline fit
  //
  Int_t    OptimizeKnots(Int_t nIter);
  Float_t  CheckKnot(Int_t iKnot);
  Bool_t   RefitKnot(Int_t iKnot);
  //
  Bool_t        fBDump;   //  dump debug information flag
  TGraph       *fGraph;   //! initial graph
  Int_t         fNmin;    //  number of points per one knot in iteration 0
  Double_t      fSigma;   //  locally estimated sigma
  Double_t      fMaxDelta;//  maximal deviation of the spline fit 
  Int_t         fN0;      //  number of knots in iteration 0 
  TClonesArray *fParams;  //  object array of parameters in knots
  TClonesArray *fCovars;  //  object array of covariance in knots
  Int_t        *fIndex;   //  [fN0] index of point corresponding to knot
  static TLinearFitter fitterStatic; // static fitter to save processing time
  //
  // 
  //
  Int_t    fN;            //  number of knots after compression 
  Double_t fChi2;         //  chi2 per degree of freedom 
  Double_t *fX;           //  [fN] - xknot value
  Double_t *fY0;          //  [fN] - y value at X
  Double_t *fY1;          //  [fN] - y derivative value at X
  Double_t *fChi2I;       //  [fN] - chi2 on interval
  ClassDef(AliSplineFit, 0);
};
#endif 
