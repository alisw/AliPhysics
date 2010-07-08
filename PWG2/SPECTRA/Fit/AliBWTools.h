// ----------------------------------------------------------------------
//                     AliBWTools
// 
// This class provides some tools which can be useful in the analsis
// of spectra, to fit or transform histograms. See the comments of the
// individual methods for details
//
// Author: M. Floris (CERN)
// ----------------------------------------------------------------------

#ifndef ALIBWTOOLS_H
#define ALIBWTOOLS_H

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TObject.h"
#include "TH1.h"

class TF1;
class TH1D;
class TH1F;
class TGraphErrors;
#endif



class AliBWTools : public TObject {

public:

  AliBWTools();
  ~AliBWTools();

  static TH1 * GetOneOverPtdNdPt(TH1 * hPt) ;
  static TH1 * GetdNdmtFromdNdpt(TH1 * hpt, Double_t mass);
  static TH1 * GetdNdPtFromOneOverPt(TH1 * h1Pt) ;

  static TGraphErrors * ConcatenateGraphs(TGraphErrors * g1,TGraphErrors * g2); 
  static TH1F *         GetHistoFromGraph(TGraphErrors * g, TH1F* hTemplate) ;  
  static TGraphErrors * GetGraphFromHisto(TH1F * h, Bool_t binWidth = kTRUE) ; 
  static TH1F *         CombineHistos(TH1 * h1, TH1 * h2, TH1* htemplate, Float_t renorm1=1.);
  static void GetFromHistoGraphDifferentX(TH1F * h, TF1 * f, TGraphErrors ** gBarycentre, TGraphErrors ** gXlw); 
  static Float_t GetMean(TH1F * h, Float_t min, Float_t max) ; 

  static void GetMean(TF1 * func, Float_t &mean, Float_t &error, Float_t min=0, Float_t max=100);
  static void GetMeanSquare(TF1 * func, Float_t &mean, Float_t &error, Float_t min=0, Float_t max=100) ;

  static Bool_t Fit (TH1 * h, TF1* f, Float_t min, Float_t max) ;

  static Int_t GetLowestNotEmptyBin(TH1*h);
  static Int_t GetHighestNotEmptyBin(TH1*h);
  static Float_t GetLowestNotEmptyBinEdge(TH1*h) { return h->GetBinLowEdge(GetLowestNotEmptyBin(h));}
  static Float_t GetHighestNotEmptyBinEdge(TH1*h) { return h->GetBinLowEdge(GetHighestNotEmptyBin(h)+1);}

  static void GetResiduals(TGraphErrors * gdata, TF1 * func, TH1F ** hres, TGraphErrors ** gres) ;
  static void GetResiduals(TH1F* hdata, TF1 * func, TH1F ** hres, TH1F ** hresVsBin) ;

  static void GetYield(TH1* h, TF1 * f, Double_t &yield, Double_t &yieldError, Float_t min = 0, 
		       Float_t max = 100,  Double_t *partialYields=0, Double_t *partialYieldsErrors=0);

  static TGraphErrors * DivideGraphByFunc (TGraphErrors * g, TF1 * f, Bool_t invert = kFALSE);
  static TGraphErrors * DivideGraphByHisto(TGraphErrors * g, TH1 * h, Bool_t invert = kFALSE);
  static TH1F         * DivideHistoByFunc (TH1F * h, TF1 * f, Bool_t invert = kFALSE);

private:

  static void GetMoment(TString name, TString var, TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max) ;
  static Double_t GetNormalizedFunc(double * x, double* p);

  static TF1 * fFuncForNormalized; // Function used in GetNormalizedFunc

  AliBWTools(const AliBWTools&);            // not implemented
  AliBWTools& operator=(const AliBWTools&); // not implemented

  ClassDef(AliBWTools,1);

};

#endif
