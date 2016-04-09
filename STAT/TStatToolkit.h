#ifndef TSTATTOOLKIT_H
#define TSTATTOOLKIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//
// some utilities which do net exist in the standard ROOT
//
/// \file TStatToolkit.h
/// \class TStatToolkit
/// \brief Summary of statistics functions
#include "TMath.h"
#include "Riostream.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3.h"
#include "THnBase.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TObjString.h"
#include "TLinearFitter.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TCut.h"
#include "THashList.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
//
// includes neccessary for test functions
//
#include "TSystem.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TTreeStream.h"
#include "AliSysInfo.h"

#include "TObject.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TMatrixD.h"
#include "TMatrixF.h"
#include <float.h>
//#include "TGraph2D.h"
//#include "TGraph.h"
class THashList;

namespace TStatToolkit
{
  enum TStatType {kEntries, kSum, kMean, kRMS, kMedian, kLTM, kLTMRMS}; 
  enum ENormType {kL1, kL2, kLp, kMax, kHamming, kNNormType };   // http://en.wikipedia.org/w/index.php?title=Norm_(mathematics)&oldid=655824636
  //
  //
  void    EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  void    EvaluateUniExternal(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh, Float_t externalfactor=1);
  Int_t  Freq(Int_t n, const Int_t *inlist, Int_t *outlist, Bool_t down);    
  //
  // HISTOGRAMS TOOLS
  //
  template <typename T> 
  void TruncatedMean(const TH1 * his, TVectorT<T> *param, Float_t down=0, Float_t up=1.0, Bool_t verbose=kFALSE);
  void MedianFilter(TH1 * his1D, Int_t nmedian);

  template <typename T> Bool_t  LTMHisto(TH1 * his, TVectorT<T> &param , Float_t fraction=1);
  template <typename T> Int_t*  LTMUnbinned(int np, const T *arr, TVectorT<T> &params , Float_t keep=1.0);

  template <typename T> void Reorder(int np, T *arr, const int *idx);
  //
  template <typename T> 
  void LTM(TH1 * his, TVectorT<T> *param=0 , Float_t fraction=1,  Bool_t verbose=kFALSE);

  template <typename T> 
  Double_t  FitGaus(TH1* his, TVectorT<T> *param=0, TMatrixT<T> *matrix=0, Float_t xmin=0, Float_t xmax=0,  Bool_t verbose=kFALSE);

  template <typename T> 
  Double_t  FitGaus(Float_t *arr, Int_t nBins, Float_t xMin, Float_t xMax, TVectorT<T> *param=0, TMatrixT<T> *matrix=0, Bool_t verbose=kFALSE);
  Float_t  GetCOG(const Short_t *arr, Int_t nBins, Float_t xMin, Float_t xMax, Float_t *rms=0, Float_t *sum=0);

  TGraph2D *  MakeStat2D(TH3 * his, Int_t delta0, Int_t delta1, Int_t type);
  TGraphErrors *  MakeStat1D(TH2 * his, Int_t deltaBin, Double_t fraction, Int_t returnType, Int_t markerStyle, Int_t markerColor);
  //
  // Graph tools
  //
  THashList *AddMetadata(TTree*, const char *vartagName,const char *varTagValue);
  TNamed *GetMetadata(TTree* tree, const char *vartagName);
  TGraph * MakeGraphSparse(TTree * tree, const char * expr="Entry", const char * cut="1",  Int_t mstyle=25, Int_t mcolor=1, Float_t msize=-1, Float_t offset=0.0);
  TGraphErrors * MakeGraphErrors(TTree * tree, const char * expr="Entry", const char * cut="1",  Int_t mstyle=25, Int_t mcolor=1, Float_t msize=-1, Float_t offset=0.0, Int_t entries=-1, Int_t firstEntry=0);

  //
  // Fitting function
  //
  TString* FitPlane(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000, Bool_t fix0=kFALSE);
  TString* FitPlaneFixed(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000);
  //
  //Linear fitter helper function
  //
  TString* FitPlaneConstrain(TTree * tree, const char* drawCommand, const char* formula, const char* cuts, Double_t & chi2, Int_t &npoints,  TVectorD &fitParam, TMatrixD &covMatrix, Float_t frac=-1, Int_t start=0, Int_t stop=10000000, Double_t constrain=-1);
  Int_t GetFitIndex(const TString fString, const TString subString);
 TString FilterFit(const TString &input, const TString filter, TVectorD &vec, TMatrixD &covar);
 void Update1D(Double_t delta, Double_t sigma, Int_t s1, TMatrixD &param, TMatrixD &covar);
  void   Constrain1D(const TString &input, const TString filter, TVectorD &param, TMatrixD & covar, Double_t mean, Double_t sigma);
  TString  MakeFitString(const TString &input, const TVectorD &param, const TMatrixD & covar, Bool_t verbose=kFALSE);
  //
  // TTree function for the trending
  //
  Int_t  MakeStatAlias(TTree * tree, const char * expr, const char * cut, const char * alias);
  Int_t  SetStatusAlias(TTree * tree, const char * expr, const char * cut, const char * alias);
  TMultiGraph*  MakeStatusMultGr(TTree * tree, const char * expr, const char * cut, const char * alias, Int_t igr=0);  
  void  AddStatusPad(TCanvas* c1, Float_t padratio, Float_t bottommargin);
  void  DrawStatusGraphs(TObjArray* oaMultGr);
  TTree*  WriteStatusToTree(TObject* oStatusGr);
  TMultiGraph*  MakeStatusLines(TTree * tree, const char * expr, const char * cut, const char * alias);
  void  MakeSummaryTree(TTree* treeIn, TTreeSRedirector *pcstream, TObjString& sumID, TCut &selection);
  Double_t GetDefaultStat(TTree * tree, const char * var, const char * selection, TStatType statType);
  //
  //
  void MakeDistortionMap(Int_t iter, THnBase * histo, TTreeSRedirector *pcstream, TMatrixD &projectionInfo, Int_t dumpHisto=100,Int_t verbose=kFALSE);
  void MakeDistortionMapFast(THnBase * histo, TTreeSRedirector *pcstream, TMatrixD &projectionInfo, Int_t verbose=0,  Double_t fractionCut=0.1, const char * estimators=0);

  //
  // norm (distance) functions
  //
  void     CombineArray(TTree *tree, TVectorD &values);
  Double_t GetDistance(const TVectorD &values, const ENormType normType,
                              const Bool_t normaliseToEntries=kFALSE, const Double_t pvalue=1.);
  Double_t GetDistance(const Int_t size, const Double_t *values, const ENormType normType,
                              const Bool_t normaliseToEntries=kFALSE, const Double_t pvalue=1.);
  Double_t GetDistance(TTree * tree, const char * var, const char * selection,
                              const ENormType normType, const Bool_t normaliseToEntries=kFALSE, const Double_t pvalue=1.);
  //
  // TTree function for robust draw
  //
  TH1* DrawHistogram(TTree * tree, const char* drawCommand, const char* cuts = "1", const char* hname = "histo", const char* htitle = "histo", Int_t nsigma = 4, Float_t fraction = 0.75, TObjArray *description=0 );
  //
  // TestFunctions:
  //
  void TestGausFit(Int_t nhistos=5000);
  void CheckTreeAliases(TTree * tree, Int_t ncheck);
 //
 // min, max, mean ...
 template <typename T> 
   void GetMinMax(const T* arr, Long64_t n, double &minVal, double &maxVal);
 template <typename T> 
   void GetMinMaxMean(const T* arr, Long64_t n, double &minVal, double &maxVal, double &meanVal);
 
};

//___________________________________________________________
template <typename T>
void TStatToolkit::GetMinMax(const T* arr, Long64_t n, double &minVal, double &maxVal)
{
  // find min, max entries in the array in a single loop
  minVal =  DBL_MAX;
  maxVal = -DBL_MAX;
  for (int i=n;i--;) {
    double val = arr[i];
    if (val<minVal) minVal = val;
    if (val>maxVal) maxVal = val;
  }
}

//___________________________________________________________
template <typename T>
void TStatToolkit::GetMinMaxMean(const T* arr, Long64_t n, double &minVal, double &maxVal, double &meanVal)
{
  // find min, max entries in the array in a single loop
  minVal =  DBL_MAX;
  maxVal = -DBL_MAX;
  meanVal = 0;
  for (int i=n;i--;) {
    double val = arr[i];
    if (val<minVal) minVal = val;
    if (val>maxVal) maxVal = val;
    meanVal += val;
  }
  if (n) meanVal /= n;
}

//___TStatToolkit__________________________________________________________________________
template <typename T> 
void TStatToolkit::TruncatedMean(const TH1 * his, TVectorT<T> *param, Float_t down, Float_t up, Bool_t verbose){
  //
  //
  //
  Int_t nbins    = his->GetNbinsX();
  Float_t nentries = his->GetEntries();
  Float_t sum      =0;
  Float_t mean   = 0;
  Float_t sigma2 = 0;
  Float_t ncumul=0;  
  for (Int_t ibin=1;ibin<nbins; ibin++){
    ncumul+= his->GetBinContent(ibin);
    Float_t fraction = Float_t(ncumul)/Float_t(nentries);
    if (fraction>down && fraction<up){
      sum+=his->GetBinContent(ibin);
      mean+=his->GetBinCenter(ibin)*his->GetBinContent(ibin);
      sigma2+=his->GetBinCenter(ibin)*his->GetBinCenter(ibin)*his->GetBinContent(ibin);      
    }
  }
  mean/=sum;
  sigma2= TMath::Sqrt(TMath::Abs(sigma2/sum-mean*mean));
  if (param){
    (*param)[0] = his->GetMaximum();
    (*param)[1] = mean;
    (*param)[2] = sigma2;
    
  }
  if (verbose)  printf("Mean\t%f\t Sigma2\t%f\n", mean,sigma2);
}

template <typename T> 
Bool_t TStatToolkit::LTMHisto(TH1 *his1D, TVectorT<T> &params , Float_t fraction){
  //
  // LTM : Trimmed mean on histogram - Modified version for binned data
  // 
  // Robust statistic to estimate properties of the distribution
  // To handle binning error special treatment
  // for definition of unbinned data see:
  //     http://en.wikipedia.org/w/index.php?title=Trimmed_estimator&oldid=582847999
  //
  // Function parameters:
  //     his1D   - input histogram
  //     params  - vector with parameters
  //             - 0 - area
  //             - 1 - mean
  //             - 2 - rms 
  //             - 3 - error estimate of mean
  //             - 4 - error estimate of RMS
  //             - 5 - first accepted bin position
  //             - 6 - last accepted  bin position
  //
  Int_t nbins    = his1D->GetNbinsX();
  Int_t nentries = (Int_t)his1D->GetEntries();
  const Double_t kEpsilon=0.0000000001;

  if (nentries<=0) return 0;
  if (fraction>1) fraction=0;
  if (fraction<0) return 0;
  TVectorD vectorX(nbins);
  TVectorD vectorMean(nbins);
  TVectorD vectorRMS(nbins);
  Double_t sumCont=0;
  for (Int_t ibin0=1; ibin0<=nbins; ibin0++) sumCont+=his1D->GetBinContent(ibin0);
  //
  Double_t minRMS=his1D->GetRMS()*10000;
  Int_t maxBin=0;
  //
  for (Int_t ibin0=1; ibin0<nbins; ibin0++){
    Double_t sum0=0, sum1=0, sum2=0;
    Int_t ibin1=ibin0;
    for ( ibin1=ibin0; ibin1<=nbins; ibin1++){
      Double_t cont=his1D->GetBinContent(ibin1);
      Double_t x= his1D->GetBinCenter(ibin1);
      sum0+=cont;
      sum1+=cont*x;
      sum2+=cont*x*x;
      if ( (ibin0!=ibin1) && sum0>=fraction*sumCont) break;
    }
    vectorX[ibin0]=his1D->GetBinCenter(ibin0);
    if (sum0<fraction*sumCont) continue;
    //
    // substract fractions of bin0 and bin1 to keep sum0=fration*sumCont
    //
    Double_t diff = sum0-fraction*sumCont;
    Double_t mean = (sum0>0) ? sum1/sum0:0;
    //
    Double_t x0=his1D->GetBinCenter(ibin0);
    Double_t x1=his1D->GetBinCenter(ibin1);
    Double_t y0=his1D->GetBinContent(ibin0);
    Double_t y1=his1D->GetBinContent(ibin1);
    //
    Double_t d = y0+y1-diff;    //enties to keep 
    Double_t w0=0,w1=0;
    if (y0<=kEpsilon&&y1>kEpsilon){
      w1=d/y1;
    } 
    if (y1<=kEpsilon&&y0>kEpsilon){
      w0=d/y0;
    }
    if (y0>kEpsilon && y1>kEpsilon && x1>x0  ){
      w0 = (d*(x1-mean))/((x1-x0)*y0);
      w1 = (d-y0*w0)/y1;
      //
      if (w0>1) {w1+=(w0-1)*y0/y1; w0=1;}
      if (w1>1) {w0+=(w1-1)*y1/y0; w1=1;}
    }  
    if ( (x1>x0) &&TMath::Abs(y0*w0+y1*w1-d)>kEpsilon*sum0){
      printf(" TStatToolkit::LTMHisto error\n");
    }
    sum0-=y0+y1;
    sum1-=y0*x0;
    sum1-=y1*x1;
    sum2-=y0*x0*x0;
    sum2-=y1*x1*x1;
    //
    Double_t xx0=his1D->GetXaxis()->GetBinUpEdge(ibin0)-0.5*w0*his1D->GetBinWidth(ibin0);
    Double_t xx1=his1D->GetXaxis()->GetBinLowEdge(ibin1)+0.5*w1*his1D->GetBinWidth(ibin1);
    sum0+=y0*w0+y1*w1;
    sum1+=y0*w0*xx0;
    sum1+=y1*w1*xx1;
    sum2+=y0*w0*xx0*xx0;
    sum2+=y1*w1*xx1*xx1;

    //
    // choose the bin with smallest rms
    //
    if (sum0>0){
      vectorMean[ibin0]=sum1/sum0;
      vectorRMS[ibin0]=TMath::Sqrt(TMath::Abs(sum2/sum0-vectorMean[ibin0]*vectorMean[ibin0]));
      if (vectorRMS[ibin0]<minRMS){
	minRMS=vectorRMS[ibin0];
	params[0]=sum0;
	params[1]=vectorMean[ibin0];
	params[2]=vectorRMS[ibin0];
	params[3]=vectorRMS[ibin0]/TMath::Sqrt(sumCont*fraction);
	params[4]=0; // what is the formula for error of RMS???
	params[5]=ibin0;
	params[6]=ibin1;
	params[7]=his1D->GetBinCenter(ibin0);
	params[8]=his1D->GetBinCenter(ibin1);
	maxBin=ibin0;
      }
    }else{
      break;
    }
  }
  return kTRUE;
}

template <typename T> 
Double_t  TStatToolkit::FitGaus(TH1* his, TVectorT<T> *param, TMatrixT<T> */*matrix*/, Float_t xmin, Float_t xmax, Bool_t verbose){
  //
  //  Fit histogram with gaussian function
  //  
  //  Prameters:
  //       return value- chi2 - if negative ( not enough points)
  //       his        -  input histogram
  //       param      -  vector with parameters 
  //       xmin, xmax -  range to fit - if xmin=xmax=0 - the full histogram range used
  //  Fitting:
  //  1. Step - make logarithm
  //  2. Linear  fit (parabola) - more robust - always converge
  //  3. In case of small statistic bins are averaged
  //  
  static TLinearFitter fitter(3,"pol2");
  TVectorD  par(3);
  TVectorD  sigma(3);
  TMatrixD mat(3,3);
  if (his->GetMaximum()<4) return -1;  
  if (his->GetEntries()<12) return -1;  
  if (his->GetRMS()<mat.GetTol()) return -1;
  Float_t maxEstimate   = his->GetEntries()*his->GetBinWidth(1)/TMath::Sqrt((TMath::TwoPi()*his->GetRMS()));
  Int_t dsmooth = TMath::Nint(6./TMath::Sqrt(maxEstimate));

  if (maxEstimate<1) return -1;
  Int_t nbins    = his->GetNbinsX();
  Int_t npoints=0;
  //


  if (xmin>=xmax){
    xmin = his->GetXaxis()->GetXmin();
    xmax = his->GetXaxis()->GetXmax();
  }
  for (Int_t iter=0; iter<2; iter++){
    fitter.ClearPoints();
    npoints=0;
    for (Int_t ibin=1;ibin<nbins+1; ibin++){
      Int_t countB=1;
      Float_t entriesI =  his->GetBinContent(ibin);
      for (Int_t delta = -dsmooth; delta<=dsmooth; delta++){
	if (ibin+delta>1 &&ibin+delta<nbins-1){
	  entriesI +=  his->GetBinContent(ibin+delta);
	  countB++;
	}
      }
      entriesI/=countB;
      Double_t xcenter= his->GetBinCenter(ibin);
      if (xcenter<xmin || xcenter>xmax) continue;
      Double_t error=1./TMath::Sqrt(countB);
      Float_t   cont=2;
      if (iter>0){
	if (par[0]+par[1]*xcenter+par[2]*xcenter*xcenter>20) return 0;
	cont = TMath::Exp(par[0]+par[1]*xcenter+par[2]*xcenter*xcenter);
	if (cont>1.) error = 1./TMath::Sqrt(cont*Float_t(countB));
      }
      if (entriesI>1&&cont>1){
	fitter.AddPoint(&xcenter,TMath::Log(Float_t(entriesI)),error);
	npoints++;
      }
    }  
    if (npoints>3){
      fitter.Eval();
      fitter.GetParameters(par);
    }else{
      break;
    }
  }
  if (npoints<=3){
    return -1;
  }
  fitter.GetParameters(par);
  fitter.GetCovarianceMatrix(mat);
  if (TMath::Abs(par[1])<mat.GetTol()) return -1;
  if (TMath::Abs(par[2])<mat.GetTol()) return -1;
  Double_t chi2 = fitter.GetChisquare()/Float_t(npoints);
  //fitter.GetParameters();
  if (!param)  param  = new TVectorT<T>(3);
  // if (!matrix) matrix = new TMatrixD(3,3); // Covariance matrix to be implemented
  (*param)[1] = par[1]/(-2.*par[2]);
  (*param)[2] = 1./TMath::Sqrt(TMath::Abs(-2.*par[2]));
  (*param)[0] = TMath::Exp(par[0]+ par[1]* (*param)[1] +  par[2]*(*param)[1]*(*param)[1]);
  if (verbose){
    par.Print();
    mat.Print();
    param->Print();
    printf("Chi2=%f\n",chi2);
    TF1 * f1= new TF1("f1","[0]*exp(-(x-[1])^2/(2*[2]*[2]))",his->GetXaxis()->GetXmin(),his->GetXaxis()->GetXmax());
    f1->SetParameter(0, (*param)[0]);
    f1->SetParameter(1, (*param)[1]);
    f1->SetParameter(2, (*param)[2]);    
    f1->Draw("same");
  }
  return chi2;
}

template <typename T> 
void TStatToolkit::LTM(TH1 * his, TVectorT<T> *param , Float_t fraction,  Bool_t verbose){
  //
  // LTM : Trimmed mean on histogram - Modified version for binned data
  //
  // Robust statistic to estimate properties of the distribution
  // See http://en.wikipedia.org/w/index.php?title=Trimmed_estimator&oldid=582847999
  //
  // New faster version is under preparation
  //
  if (!param) return;
  (*param)[0]=0;
  (*param)[1]=0;
  (*param)[2]=0;  
  Int_t nbins    = his->GetNbinsX();
  Int_t nentries = (Int_t)his->GetEntries();
  if (nentries<=0) return;
  Double_t *data  = new Double_t[nentries];
  Int_t npoints=0;
  for (Int_t ibin=1;ibin<nbins; ibin++){
    Double_t entriesI = his->GetBinContent(ibin);
    //Double_t xcenter= his->GetBinCenter(ibin);
    Double_t x0 =  his->GetXaxis()->GetBinLowEdge(ibin);
    Double_t w  =  his->GetXaxis()->GetBinWidth(ibin);
    for (Int_t ic=0; ic<entriesI; ic++){
      if (npoints<nentries){
	data[npoints]= x0+w*Double_t((ic+0.5)/entriesI);
	npoints++;
      }
    }
  }
  Double_t mean, sigma;  
  Int_t npoints2=TMath::Min(Int_t(fraction*Float_t(npoints)),npoints-1);
  npoints2=TMath::Max(Int_t(0.5*Float_t(npoints)),npoints2);
  TStatToolkit::EvaluateUni(npoints, data, mean,sigma,npoints2);
  delete [] data;
  if (verbose)  printf("Mean\t%f\t Sigma2\t%f\n", mean,sigma);if (param){
    (*param)[0] = his->GetMaximum();
    (*param)[1] = mean;
    (*param)[2] = sigma;    
  }
}


template <typename T> 
Double_t  TStatToolkit::FitGaus(Float_t *arr, Int_t nBins, Float_t xMin, Float_t xMax, TVectorT<T> *param, TMatrixT<T> */*matrix*/, Bool_t verbose){
  //
  //  Fit histogram with gaussian function
  //  
  //  Prameters:
  //     nbins: size of the array and number of histogram bins
  //     xMin, xMax: histogram range
  //     param: paramters of the fit (0-Constant, 1-Mean, 2-Sigma)
  //     matrix: covariance matrix -- not implemented yet, pass dummy matrix!!!
  //
  //  Return values:
  //    >0: the chi2 returned by TLinearFitter
  //    -3: only three points have been used for the calculation - no fitter was used
  //    -2: only two points have been used for the calculation - center of gravity was uesed for calculation
  //    -1: only one point has been used for the calculation - center of gravity was uesed for calculation
  //    -4: invalid result!!
  //
  //  Fitting:
  //  1. Step - make logarithm
  //  2. Linear  fit (parabola) - more robust - always converge
  //  
  static TLinearFitter fitter(3,"pol2");
  static TMatrixD mat(3,3);
  static Double_t kTol = mat.GetTol();
  fitter.StoreData(kFALSE);
  fitter.ClearPoints();
  TVectorD  par(3);
  TVectorD  sigma(3);
  TMatrixD matA(3,3);
  TMatrixD b(3,1);
  Float_t rms = TMath::RMS(nBins,arr);
  Float_t max = TMath::MaxElement(nBins,arr);
  Float_t binWidth = (xMax-xMin)/(Float_t)nBins;

  Float_t meanCOG = 0;
  Float_t rms2COG = 0;
  Float_t sumCOG  = 0;

  Float_t entries = 0;
  Int_t nfilled=0;

  for (Int_t i=0; i<nBins; i++){
      entries+=arr[i];
      if (arr[i]>0) nfilled++;
  }

  if (max<4) return -4;
  if (entries<12) return -4;
  if (rms<kTol) return -4;

  Int_t npoints=0;
  //

  //
  for (Int_t ibin=0;ibin<nBins; ibin++){
      Float_t entriesI = arr[ibin];
    if (entriesI>1){
      Double_t xcenter = xMin+(ibin+0.5)*binWidth;
      
      Float_t error    = 1./TMath::Sqrt(entriesI);
      Float_t val = TMath::Log(Float_t(entriesI));
      fitter.AddPoint(&xcenter,val,error);
      if (npoints<3){
	  matA(npoints,0)=1;
	  matA(npoints,1)=xcenter;
	  matA(npoints,2)=xcenter*xcenter;
	  b(npoints,0)=val;
	  meanCOG+=xcenter*entriesI;
	  rms2COG +=xcenter*entriesI*xcenter;
	  sumCOG +=entriesI;
      }
      npoints++;
    }
  }

  
  Double_t chi2 = 0;
  if (npoints>=3){
      if ( npoints == 3 ){
	  //analytic calculation of the parameters for three points
	  matA.Invert();
	  TMatrixD res(1,3);
	  res.Mult(matA,b);
	  par[0]=res(0,0);
	  par[1]=res(0,1);
	  par[2]=res(0,2);
          chi2 = -3.;
      } else {
          // use fitter for more than three points
	  fitter.Eval();
	  fitter.GetParameters(par);
	  fitter.GetCovarianceMatrix(mat);
	  chi2 = fitter.GetChisquare()/Float_t(npoints);
      }
      if (TMath::Abs(par[1])<kTol) return -4;
      if (TMath::Abs(par[2])<kTol) return -4;

      if (!param)  param  = new TVectorT<T>(3);
      //if (!matrix) matrix = new TMatrixT<T>(3,3);  // !!!!might be a memory leek. use dummy matrix pointer to call this function! // Covariance matrix to be implemented

      (*param)[1] = par[1]/(-2.*par[2]);
      (*param)[2] = 1./TMath::Sqrt(TMath::Abs(-2.*par[2]));
      Double_t lnparam0 = par[0]+ par[1]* (*param)[1] +  par[2]*(*param)[1]*(*param)[1];
      if ( lnparam0>307 ) return -4;
      (*param)[0] = TMath::Exp(lnparam0);
      if (verbose){
	  par.Print();
	  mat.Print();
	  param->Print();
	  printf("Chi2=%f\n",chi2);
	  TF1 * f1= new TF1("f1","[0]*exp(-(x-[1])^2/(2*[2]*[2]))",xMin,xMax);
	  f1->SetParameter(0, (*param)[0]);
	  f1->SetParameter(1, (*param)[1]);
	  f1->SetParameter(2, (*param)[2]);
	  f1->Draw("same");
      }
      return chi2;
  }

  if (npoints == 2){
      //use center of gravity for 2 points
      meanCOG/=sumCOG;
      rms2COG /=sumCOG;
      (*param)[0] = max;
      (*param)[1] = meanCOG;
      (*param)[2] = TMath::Sqrt(TMath::Abs(meanCOG*meanCOG-rms2COG));
      chi2=-2.;
  }
  if ( npoints == 1 ){
      meanCOG/=sumCOG;
      (*param)[0] = max;
      (*param)[1] = meanCOG;
      (*param)[2] = binWidth/TMath::Sqrt(12);
      chi2=-1.;
  }
  return chi2;

}

template <typename T> 
Int_t* TStatToolkit::LTMUnbinned(int np, const T *arr, TVectorT<T> &params , Float_t keep)
{
  //
  // LTM : Trimmed mean of unbinned array
  // 
  // Robust statistic to estimate properties of the distribution
  // To handle binning error special treatment
  // for definition of unbinned data see:
  //     http://en.wikipedia.org/w/index.php?title=Trimmed_estimator&oldid=582847999
  //
  // Function parameters:
  //     np      - number of points in the array
  //     arr     - data array (unsorted)
  //     params  - vector with parameters
  //             - 0 - area
  //             - 1 - mean
  //             - 2 - rms 
  //             - 3 - error estimate of mean
  //             - 4 - error estimate of RMS
  //             - 5 - first accepted element (of sorted array)
  //             - 6 - last accepted  element (of sorted array)
  //
  // On success returns index of sorted events 
  //
  static int book = 0;
  static int *index = 0;
  static double* w = 0;
  int keepN = np*keep;
  if (keepN>np) keepN = np;
  if (keepN<2) return 0;
  params[0] = 0.0f;
  if (book<np) {
    delete[] index;
    book = np;
    index = new int[book];
    delete[] w;
    w = new double[book+book];
  }
  //
  double *wx1 = w, *wx2 = wx1+np;
  TMath::Sort(np,arr,index,kFALSE); // sort in increasing order
  // build cumulants
  double sum1=0.0,sum2=0.0;
  for (int i=0;i<np;i++) {
    double x = arr[index[i]];
    wx1[i] = (sum1+=x);
    wx2[i] = (sum2+=x*x);
  }
  double minRMS = sum2+1e6;
  params[0] = keepN;
  int limI = np - keepN+1;
  for (int i=0;i<limI;i++) {
    int limJ = i+keepN-1;
    Double_t sum1 = wx1[limJ] - (i ? wx1[i-1] : 0.0);
    Double_t sum2 = wx2[limJ] - (i ? wx2[i-1] : 0.0);
    double mean = sum1/keepN;
    double rms2 = sum2/keepN - mean*mean;
    //    printf("%d : %d %e %e\n",i,limJ, mean, TMath::Sqrt(rms2));
    if (rms2>minRMS) continue;
    minRMS = rms2;
    params[1] = mean;
    params[2] = rms2;
    params[5] = i;
    params[6] = limJ;
  }
  //
  if (!params[0]) return 0;
  params[2] = TMath::Sqrt(params[2]);
  params[3] = params[2]/TMath::Sqrt(params[0]); // error on mean
  params[4] = params[3]/TMath::Sqrt(2.0); // error on RMS
  return index;
}

template <typename T> 
void TStatToolkit::Reorder(int np, T *arr, const int *idx)
{
  // rearange arr in order given by idx
  T arrc[np];
  memcpy(arrc,arr,np*sizeof(T));
  for (int i=np;i--;) arr[i] = arrc[idx[i]];
}

#endif
