#ifndef ALIMCLOGLFITTER_H
#define ALIMCLOGLFITTER_H


// Code by Martin Völkl:  martin.andreas.volkl@cern.ch
#include "TObject.h"

class TH1;
class TH1D;
class TH2D;
class AliMCLogLFitter;


class AliMCLogLFitter
{
public:
  AliMCLogLFitter(Int_t nMCFunctions, TH1 ** MCDistributions, TH1 * dataHistogram); 
  ~AliMCLogLFitter();
  void Fit(void);
  void ReturnNegLog(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  double ReturnParameter(Int_t i){return fParameter[i];}
  void ChangeFunctions(Int_t nMCFunctions, TH1 ** MCDistributions, TH1 * dataHistogram);
  void SetTemplateWeight(Int_t parameterNr, TH1D * Weights);
  void SetFunctionCorrelation(int function, int CoupleTo, double ParameterRatio);
  void SetFitRange(double lowerBound, double upperBound);
  void SetParameter(int parameterNr, double startValue, double lowerBound=0.01, double upperBound=1.0);
  void IfNoMCParticleIsProbablyA(Int_t parameterNr){fDefaultParameterInLowBins=parameterNr;}
  void FitVeryCarefully(void){fCarefulFitting=true;}
  void SetUseLogL(void){fUseLogL=true;}
  void SetUseChiSq(void){fUseLogL=false;}
  TH1D * ScanVariable(Int_t parameterNr, Double_t min, Double_t max);
  TH2D * ScanVariable2d(Int_t parameterNr, Double_t min, Double_t max, Int_t parameterNr2, Double_t min2, Double_t max2);
  TH1 * ReturnExpectedDiagram(int i){return fFitDiagsMCExpected[i];}
  Double_t GetLogLikelyhood(void);
  void ScanAll(void);
  
private:
  void UpdateExpectationIteratively(Double_t * par, Int_t bin);
  Double_t WeightedIntegralOfSource(Int_t j);
  Double_t fBestLikelyhood; // Internal current local maximum
  Double_t w(int j, int i); // returns weight if available, else 1
  Int_t fNPar; // Number of free parameters (number of templates for fitting)
  Int_t fNCoupledFunctions; // Number of functions with parameters coupled to others, these are always the last in the array
  Int_t fDefaultParameterInLowBins; // Default nonzero parameter for bins with no MC counts
  Double_t * fParameter; // Fit amplitude parameters
  Double_t * fBestParameters; // Variable to save currenct local maximum internally
  Double_t * fStartParameter; // Start parameters for fitting
  Double_t * fUpperLimitParameter; // Upper limit of allowable parameter values
  Double_t * fLowerLimitParameter; // Lower limit of allowable parameter values
  Int_t * fCouplings; // Number of function that is coupled
  Double_t * fCouplingParameter; // Ratios of the Fit parameters
  Int_t flowerFitRange, fupperFitRange; // Fit range given as bin #
  TH1 ** fFitDiagsMC; // The MC templates
  TH1 ** fFitDiagsMCExpected; // Additional free parameters, A_ji, expectation values of the templates
  TH1 * fFitDiagData; // Data input histogram
  Bool_t * WeightsAvailable;
  TH1D ** WeightsHistograms;
  Int_t fiterations; // Iterations within one bin
  Bool_t fCarefulFitting; // Switch for scanning some variables to look for the global maximum
  Bool_t fUseLogL; // Switch to use LogL of Chi2 method
  int fTotalCalls; // Number of times, the LogL was calculated - #iterations of fit
  bool scan2d; // switch to scan two parameters completely
};

namespace SomewhereElseInAliMCLogLFitter
{
  AliMCLogLFitter * fitter;
  void OutsourcedReturnNegLog(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
  {
    fitter->ReturnNegLog(npar, gin, f, par, iflag);
  }
}

#endif
