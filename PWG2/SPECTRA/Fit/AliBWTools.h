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

  static TH1 * GetOneOverPtdNdPt(const TH1 * hPt) ;
  static TH1 * GetdNdmtFromdNdpt(const TH1 * hpt, Double_t mass);
  static TH1 * GetdNdPtFromOneOverPt(const TH1 * h1Pt) ;
  static TH1 * GetdNdptFromdNdmt(const TH1 * hmt, Double_t mass) ;

  static TGraphErrors * ConcatenateGraphs(const TGraphErrors * g1,const TGraphErrors * g2); 
  static TH1F *         GetHistoFromGraph(const TGraphErrors * g, const TH1F* hTemplate) ;  
  static TGraphErrors * GetGraphFromHisto(const TH1F * h, Bool_t binWidth = kTRUE) ; 
  static void GetMeanDataAndExtrapolation(const TH1 * hData, TF1 * fExtrapolation, Double_t &mean, Double_t &error, Float_t min=0, Float_t max=100);

  static TH1F *         CombineHistos(const TH1 * h1, TH1 * h2, const TH1* htemplate, Float_t renorm1=1.);
  static TH1F *         Combine3HistosWithErrors(const TH1 * h1,  const TH1 * h2,  const TH1* h3, 
						  TH1 * he1,  TH1 * he2,  TH1 * he3, 
						 const TH1* htemplate, Int_t statFrom = 0, 
						 Float_t renorm1=1., Float_t renorm2=1., Float_t renorm3=1.,
						 TH1 ** hSyst =0, Bool_t errorFromBinContent=kFALSE);
  static void GetFromHistoGraphDifferentX(const TH1F * h, TF1 * f, TGraphErrors ** gBarycentre, TGraphErrors ** gXlw); 
  static Float_t GetMean(TH1F * h, Float_t min, Float_t max, Float_t * error = NULL) ; 

  static void GetMean(TF1 * func, Float_t &mean, Float_t &error, Float_t min=0, Float_t max=100, Int_t normPar = -1);
  static void GetMeanSquare(TF1 * func, Float_t &mean, Float_t &error, Float_t min=0, Float_t max=100, Int_t normPar = -1) ;

  static Bool_t Fit (TH1 * h, TF1* f, Float_t min, Float_t max) ;

  static Int_t GetLowestNotEmptyBin(const TH1*h);
  static Int_t GetHighestNotEmptyBin(const TH1*h);
  static Float_t GetLowestNotEmptyBinEdge(const TH1*h) { return h->GetBinLowEdge(GetLowestNotEmptyBin(h));}
  static Float_t GetHighestNotEmptyBinEdge(const TH1*h) { return h->GetBinLowEdge(GetHighestNotEmptyBin(h)+1);}

  static void GetResiduals(const TGraphErrors * gdata, const TF1 * func, TH1F ** hres, TGraphErrors ** gres) ;
  static void GetResiduals(const TH1F* hdata, const TF1 * func, TH1F ** hres, TH1F ** hresVsBin) ;

  static void GetYield(TH1* h, TF1 * f, Double_t &yield, Double_t &yieldError, Float_t min = 0, 
		       Float_t max = 100,  Double_t *partialYields=0, Double_t *partialYieldsErrors=0);

  static TGraphErrors * DivideGraphByFunc (const TGraphErrors * g, const TF1 * f, Bool_t invert = kFALSE);
  static TGraphErrors * DivideGraphByHisto(const TGraphErrors * g, TH1 * h, Bool_t invert = kFALSE);
  static TH1F         * DivideHistoByFunc (TH1F * h, TF1 * f, Bool_t invert = kFALSE);

  static void WeightedMean(Int_t npoints, const Double_t *x, const Double_t *xerr, Double_t &mean, Double_t &meanerr);

  static void GetValueAndError(TH1 * hdest, const TH1 * hvalue, const TH1 * herror, Bool_t isPercentError) ;  
  static TH1 * GetRelativeError(TH1 * h);
  static void AddHisto(TH1 * hdest, const TH1* hsource, Bool_t getMirrorBins = kFALSE);
  static void GetHistoCombinedErrors(TH1 * hdest, const TH1 * h1) ;
  static TH1F * DivideHistosDifferentBins(const TH1F* h1, const TH1F* h2);
  static Double_t DoIntegral(TH1* h, Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Double_t & error ,
		      Option_t *option, Bool_t doError) ;
  static TGraphErrors * Divide2Graphs(const TGraphErrors * g1, const TGraphErrors * g2);

  static Double_t dMtdptFunction(Double_t *x, Double_t *p) ;
  static Double_t GetdMtdEta(TH1 *hData, TF1 * fExtrapolation, Double_t mass) ;



private:

  AliBWTools(const AliBWTools&);            // not implemented
  AliBWTools& operator=(const AliBWTools&); // not implemented
  static void GetMoment(TString name, TString var, TF1 * func, Float_t &mean, Float_t &error, Float_t min, Float_t max, Int_t normPar = -1) ;

  static TF1 * fdNdptForETCalc;

  ClassDef(AliBWTools,1);

};

#endif
