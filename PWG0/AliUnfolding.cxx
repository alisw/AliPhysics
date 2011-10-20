/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliUnfolding.cxx 31168 2009-02-23 15:18:45Z jgrosseo $ */

// This class allows 1-dimensional unfolding.
// Methods that are implemented are chi2 minimization and bayesian unfolding.
//
//  Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

#include "AliUnfolding.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TVirtualFitter.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TExec.h>
#include "Riostream.h"
#include "TROOT.h"

using namespace std; //required for resolving the 'cout' symbol

TMatrixD* AliUnfolding::fgCorrelationMatrix = 0;
TMatrixD* AliUnfolding::fgCorrelationMatrixSquared = 0;
TMatrixD* AliUnfolding::fgCorrelationCovarianceMatrix = 0;
TVectorD* AliUnfolding::fgCurrentESDVector = 0;
TVectorD* AliUnfolding::fgEntropyAPriori = 0;
TVectorD* AliUnfolding::fgEfficiency = 0;

TAxis* AliUnfolding::fgUnfoldedAxis = 0;
TAxis* AliUnfolding::fgMeasuredAxis = 0;

TF1* AliUnfolding::fgFitFunction = 0;

AliUnfolding::MethodType AliUnfolding::fgMethodType = AliUnfolding::kInvalid;
Int_t AliUnfolding::fgMaxInput  = -1;  // bins in measured histogram
Int_t AliUnfolding::fgMaxParams = -1;  // bins in unfolded histogram = number of fit params
Float_t AliUnfolding::fgOverflowBinLimit = -1;

AliUnfolding::RegularizationType AliUnfolding::fgRegularizationType = AliUnfolding::kPol1;
Float_t AliUnfolding::fgRegularizationWeight = 10000;
Int_t AliUnfolding::fgSkipBinsBegin = 0;
Float_t AliUnfolding::fgMinuitStepSize = 0.1;                 // (usually not needed to be changed) step size in minimization
Float_t AliUnfolding::fgMinuitPrecision = 1e-6;               // minuit precision
Int_t   AliUnfolding::fgMinuitMaxIterations = 5000;           // minuit maximum number of iterations
Bool_t AliUnfolding::fgMinimumInitialValue = kFALSE;          // set all initial values at least to the smallest value among the initial values
Float_t AliUnfolding::fgMinimumInitialValueFix = -1;
Bool_t AliUnfolding::fgNormalizeInput = kFALSE;               // normalize input spectrum
Float_t AliUnfolding::fgNotFoundEvents = 0;
Bool_t AliUnfolding::fgSkipBin0InChi2 = kFALSE;

Float_t AliUnfolding::fgBayesianSmoothing  = 1;           // smoothing parameter (0 = no smoothing)
Int_t   AliUnfolding::fgBayesianIterations = 10;          // number of iterations in Bayesian method

Bool_t AliUnfolding::fgDebug = kFALSE;

Int_t AliUnfolding::fgCallCount = 0;

Int_t AliUnfolding::fgPowern = 5;

Double_t AliUnfolding::fChi2FromFit = 0.;
Double_t AliUnfolding::fPenaltyVal  = 0.;
Double_t AliUnfolding::fAvgResidual = 0.;

Int_t AliUnfolding::fgPrintChi2Details = 0;

TCanvas *AliUnfolding::fgCanvas = 0;
TH1 *AliUnfolding::fghUnfolded = 0;     
TH2 *AliUnfolding::fghCorrelation = 0;  
TH1 *AliUnfolding::fghEfficiency = 0;   
TH1 *AliUnfolding::fghMeasured = 0;     

ClassImp(AliUnfolding)

//____________________________________________________________________
void AliUnfolding::SetUnfoldingMethod(MethodType methodType)
{
  // set unfolding method
  fgMethodType = methodType; 
  
  const char* name = 0;
  switch (methodType)
  {
    case kInvalid: name = "INVALID"; break;
    case kChi2Minimization: name = "Chi2 Minimization"; break;
    case kBayesian: name = "Bayesian unfolding"; break;
    case kFunction: name = "Functional fit"; break;
  }
  Printf("AliUnfolding::SetUnfoldingMethod: %s enabled.", name);
}

//____________________________________________________________________
void AliUnfolding::SetCreateOverflowBin(Float_t overflowBinLimit) 
{ 
  // enable the creation of a overflow bin that includes all statistics below the given limit
  
  fgOverflowBinLimit = overflowBinLimit; 
  
  Printf("AliUnfolding::SetCreateOverflowBin: overflow bin limit set to %f", overflowBinLimit);
}

//____________________________________________________________________
void AliUnfolding::SetSkipBinsBegin(Int_t nBins)
{
  // set number of skipped bins in regularization
  
  fgSkipBinsBegin = nBins;
  
  Printf("AliUnfolding::SetSkipBinsBegin: skipping %d bins at the beginning of the spectrum in the regularization.", fgSkipBinsBegin);
}

//____________________________________________________________________
void AliUnfolding::SetNbins(Int_t nMeasured, Int_t nUnfolded) 
{ 
  // set number of bins in the input (measured) distribution and in the unfolded distribution
  fgMaxInput = nMeasured; 
  fgMaxParams = nUnfolded; 
  
  if (fgCorrelationMatrix)
  {
    delete fgCorrelationMatrix;
    fgCorrelationMatrix = 0;
  }
  if (fgCorrelationMatrixSquared)
  {
    fgCorrelationMatrixSquared = 0;
    delete fgCorrelationMatrixSquared;
  }
  if (fgCorrelationCovarianceMatrix)
  {
    delete fgCorrelationCovarianceMatrix;
    fgCorrelationCovarianceMatrix = 0;
  }
  if (fgCurrentESDVector)
  {
    delete fgCurrentESDVector;
    fgCurrentESDVector = 0;
  }
  if (fgEntropyAPriori)
  {
    delete fgEntropyAPriori;
    fgEntropyAPriori = 0;
  }
  if (fgEfficiency)
  {
    delete fgEfficiency;
    fgEfficiency = 0;
  }
  if (fgUnfoldedAxis)
  {
    delete fgUnfoldedAxis;
    fgUnfoldedAxis = 0;
  }
  if (fgMeasuredAxis)
  {
    delete fgMeasuredAxis;
    fgMeasuredAxis = 0;
  }
  
  Printf("AliUnfolding::SetNbins: Set %d measured bins and %d unfolded bins", nMeasured, nUnfolded);
}

//____________________________________________________________________
void AliUnfolding::SetChi2Regularization(RegularizationType type, Float_t weight)
{
  //
  // sets the parameters for chi2 minimization
  //

  fgRegularizationType = type;
  fgRegularizationWeight = weight;

  Printf("AliUnfolding::SetChi2Regularization --> Regularization set to %d with weight %f", (Int_t) type, weight);
}

//____________________________________________________________________
void AliUnfolding::SetBayesianParameters(Float_t smoothing, Int_t nIterations)
{
  //
  // sets the parameters for Bayesian unfolding
  //

  fgBayesianSmoothing = smoothing;
  fgBayesianIterations = nIterations;

  Printf("AliUnfolding::SetBayesianParameters --> Paramaters set to %d iterations with smoothing %f", fgBayesianIterations, fgBayesianSmoothing);
}

//____________________________________________________________________
void AliUnfolding::SetFunction(TF1* function)
{
  // set function for unfolding with a fit function
  
  fgFitFunction = function;
  
  Printf("AliUnfolding::SetFunction: Set fit function with %d parameters.", function->GetNpar());
}

//____________________________________________________________________
Int_t AliUnfolding::Unfold(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check)
{
  // unfolds with unfolding method fgMethodType
  //
  // parameters:
  //  correlation: response matrix as measured vs. generated
  //  efficiency:  (optional) efficiency that is applied on the unfolded spectrum, i.e. it has to be in unfolded variables. If 0 no efficiency is applied.
  //  measured:    the measured spectrum
  //  initialConditions: (optional) initial conditions for the unfolding. if 0 the measured spectrum is used as initial conditions.
  //  result:      target for the unfolded result
  //  check:       depends on the unfolding method, see comments in specific functions
  //
  //  return code: see UnfoldWithMinuit/UnfoldWithBayesian/UnfoldWithFunction

  if (fgMaxInput == -1)
  {
    Printf("AliUnfolding::Unfold: WARNING. Number of measured bins not set with SetNbins. Using number of bins in measured distribution");
    fgMaxInput = measured->GetNbinsX();
  }
  if (fgMaxParams == -1)
  {
    Printf("AliUnfolding::Unfold: WARNING. Number of unfolded bins not set with SetNbins. Using number of bins in measured distribution");
    fgMaxParams = measured->GetNbinsX();
  }

  if (fgOverflowBinLimit > 0)
    CreateOverflowBin(correlation, measured);
    
  switch (fgMethodType)
  {
    case kInvalid:
    {
      Printf("AliUnfolding::Unfold: ERROR: Unfolding method not set. Use SetUnfoldingMethod. Exiting...");
      return -1;
    }
    case kChi2Minimization:
      return UnfoldWithMinuit(correlation, efficiency, measured, initialConditions, result, check);
    case kBayesian:
      return UnfoldWithBayesian(correlation, efficiency, measured, initialConditions, result);
    case kFunction:
      return UnfoldWithFunction(correlation, efficiency, measured, initialConditions, result);
  }



  return -1;
}

//____________________________________________________________________
void AliUnfolding::SetStaticVariables(TH2* correlation, TH1* measured, TH1* efficiency)
{
  // fill static variables needed for minuit fit

  if (!fgCorrelationMatrix)
    fgCorrelationMatrix = new TMatrixD(fgMaxInput, fgMaxParams);
  if (!fgCorrelationMatrixSquared)
    fgCorrelationMatrixSquared = new TMatrixD(fgMaxInput, fgMaxParams);
  if (!fgCorrelationCovarianceMatrix)
    fgCorrelationCovarianceMatrix = new TMatrixD(fgMaxInput, fgMaxInput);
  if (!fgCurrentESDVector)
    fgCurrentESDVector = new TVectorD(fgMaxInput);
  if (!fgEntropyAPriori)
    fgEntropyAPriori = new TVectorD(fgMaxParams);
  if (!fgEfficiency)
    fgEfficiency = new TVectorD(fgMaxParams);
  if (!fgUnfoldedAxis)
    delete fgUnfoldedAxis;
  fgUnfoldedAxis = new TAxis(*(correlation->GetXaxis()));
  if (!fgMeasuredAxis)
    delete fgMeasuredAxis;
  fgMeasuredAxis = new TAxis(*(correlation->GetYaxis()));    

  fgCorrelationMatrix->Zero();
  fgCorrelationCovarianceMatrix->Zero();
  fgCurrentESDVector->Zero();
  fgEntropyAPriori->Zero();

  // normalize correction for given nPart
  for (Int_t i=1; i<=correlation->GetNbinsX(); ++i)
  {
    Double_t sum = correlation->Integral(i, i, 1, correlation->GetNbinsY());
    if (sum <= 0)
      continue;
    Float_t maxValue = 0;
    Int_t maxBin = -1;
    for (Int_t j=1; j<=correlation->GetNbinsY(); ++j)
    {
      // find most probably value
      if (maxValue < correlation->GetBinContent(i, j))
      {
        maxValue = correlation->GetBinContent(i, j);
        maxBin = j;
      }

      // npart sum to 1
      correlation->SetBinContent(i, j, correlation->GetBinContent(i, j) / sum);// * correlation->GetXaxis()->GetBinWidth(i));
      correlation->SetBinError(i, j, correlation->GetBinError(i, j) / sum);

      if (i <= fgMaxParams && j <= fgMaxInput)
      {
        (*fgCorrelationMatrix)(j-1, i-1) = correlation->GetBinContent(i, j);
        (*fgCorrelationMatrixSquared)(j-1, i-1) = correlation->GetBinContent(i, j) * correlation->GetBinContent(i, j);
      }
    }

    //printf("MPV for Ntrue = %f is %f\n", fCurrentCorrelation->GetXaxis()->GetBinCenter(i), fCurrentCorrelation->GetYaxis()->GetBinCenter(maxBin));
  }
    
  //normalize measured
  Float_t smallestError = 1;
  if (fgNormalizeInput)
  {
    Float_t sumMeasured = measured->Integral();
    measured->Scale(1.0 / sumMeasured);
    smallestError /= sumMeasured;
  }
  
  for (Int_t i=0; i<fgMaxInput; ++i)
  {
    (*fgCurrentESDVector)[i] = measured->GetBinContent(i+1);
    if (measured->GetBinError(i+1) > 0)
    {
      (*fgCorrelationCovarianceMatrix)(i, i) = (Double_t) 1e-6 / measured->GetBinError(i+1) / measured->GetBinError(i+1);
    }
    else // in this case put error of 1, otherwise 0 bins are not added to the chi2...
      (*fgCorrelationCovarianceMatrix)(i, i) = (Double_t) 1e-6 / smallestError / smallestError;

    if ((*fgCorrelationCovarianceMatrix)(i, i) > 1e7)
      (*fgCorrelationCovarianceMatrix)(i, i) = 0;
    //Printf("%d, %e", i, (*fgCorrelationCovarianceMatrix)(i, i));
  }

  // efficiency is expected to match bin width of result
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    (*fgEfficiency)(i) = efficiency->GetBinContent(i+1);
  }

  if (correlation->GetNbinsX() != fgMaxParams || correlation->GetNbinsY() != fgMaxInput)
    cout << "Response histo has incorrect dimensions; expect (" << fgMaxParams << ", " << fgMaxInput << "), got (" << correlation->GetNbinsX() << ", " << correlation->GetNbinsY() << ")" << endl;

}

//____________________________________________________________________
Int_t AliUnfolding::UnfoldWithMinuit(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result, Bool_t check)
{
  //
  // implementation of unfolding (internal function)
  //
  // unfolds <measured> using response from <correlation> and effiency <efficiency>
  // output is in <result>
  // <initialConditions> set the initial values for the minimization, if 0 <measured> is used
  //   negative values in initialConditions mean that the given parameter is fixed to the absolute of the value
  // if <check> is true no unfolding is made, instead only the chi2 without unfolding is printed
  //
  // returns minuit status (0 = success), (-1 when check was set)
  //

  SetStaticVariables(correlation, measured, efficiency);
  
  // Initialize TMinuit via generic fitter interface
  Int_t params = fgMaxParams;
  if (fgNotFoundEvents > 0)
    params++;
    
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, params);
  Double_t arglist[100];
  //  minuit->SetDefaultFitter("Minuit2");

  // disable any output (-1), unfortuantly we do not see warnings anymore then. Have to find another way...
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);

  // however, enable warnings
  //minuit->ExecuteCommand("SET WAR", arglist, 0);

  // set minimization function
  minuit->SetFCN(Chi2Function);

  // set precision
  minuit->SetPrecision(fgMinuitPrecision);

  minuit->SetMaxIterations(fgMinuitMaxIterations);

  for (Int_t i=0; i<fgMaxParams; i++)
    (*fgEntropyAPriori)[i] = 1;

  // set initial conditions as a-priori distribution for MRX regularization
  /*
  for (Int_t i=0; i<fgMaxParams; i++)
    if (initialConditions && initialConditions->GetBinContent(i+1) > 0)
      (*fgEntropyAPriori)[i] = initialConditions->GetBinContent(i+1);
  */

  if (!initialConditions) {
    initialConditions = measured;
  } else {
    Printf("AliUnfolding::UnfoldWithMinuit: Using different initial conditions...");
    //new TCanvas; initialConditions->DrawCopy();
    if (fgNormalizeInput)
      initialConditions->Scale(1.0 / initialConditions->Integral());
  }

  // extract minimum value from initial conditions (if we set a value to 0 it will stay 0)
  Float_t minValue = 1e35;
  if (fgMinimumInitialValueFix < 0)
  {
    for (Int_t i=0; i<fgMaxParams; ++i)
    {
      Int_t bin = initialConditions->GetXaxis()->FindBin(result->GetXaxis()->GetBinCenter(i+1));
      if (initialConditions->GetBinContent(bin) > 0)
	minValue = TMath::Min(minValue, (Float_t) initialConditions->GetBinContent(bin));
    }
  }
  else
    minValue = fgMinimumInitialValueFix;
  
  Double_t* results = new Double_t[fgMaxParams+1];
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    Int_t bin = initialConditions->GetXaxis()->FindBin(result->GetXaxis()->GetBinCenter(i+1));
    results[i] = initialConditions->GetBinContent(bin);

    Bool_t fix = kFALSE;
    if (results[i] < 0)
    {
      fix = kTRUE;
      results[i] = -results[i];
    }
 
    if (!fix && fgMinimumInitialValue && results[i] < minValue)
      results[i] = minValue;
      
    // minuit sees squared values to prevent it from going negative...
    results[i] = TMath::Sqrt(results[i]);

    minuit->SetParameter(i, Form("param%d", i), results[i], (fix) ? 0 : fgMinuitStepSize, 0, 0);
  }
  if (fgNotFoundEvents > 0)
  {
    results[fgMaxParams] = efficiency->GetBinContent(1);
    minuit->SetParameter(fgMaxParams, "vtx0", results[fgMaxParams], fgMinuitStepSize / 100, 0.01, 0.80);
  }
  
  Int_t dummy = 0;
  Double_t chi2 = 0;
  Chi2Function(dummy, 0, chi2, results, 0);
  printf("AliUnfolding::UnfoldWithMinuit: Chi2 of initial parameters is = %f\n", chi2);

  if (check)
  {
    DrawGuess(results);
    delete[] results;
    return -1;
  }

  // first param is number of iterations, second is precision....
  arglist[0] = 1e6;
  //arglist[1] = 1e-5;
  //  minuit->ExecuteCommand("SET PRINT", arglist, 3);
  //  minuit->ExecuteCommand("SCAN", arglist, 0);
  Int_t status = minuit->ExecuteCommand("MIGRAD", arglist, 1);
  Printf("AliUnfolding::UnfoldWithMinuit: MINUIT status is %d", status);
  //printf("!!!!!!!!!!!!!! MIGRAD finished: Starting MINOS !!!!!!!!!!!!!!");
  //minuit->ExecuteCommand("MINOS", arglist, 0);

  if (fgNotFoundEvents > 0)
  {
    results[fgMaxParams] = minuit->GetParameter(fgMaxParams);
    Printf("Efficiency for bin 0 changed from %f to %f", efficiency->GetBinContent(1), results[fgMaxParams]);
    efficiency->SetBinContent(1, results[fgMaxParams]);
  }
  
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    results[i] = minuit->GetParameter(i);
    Double_t value = results[i] * results[i];
   // error is : 2 * (relError on results[i]) * (value) = 2 * (minuit->GetParError(i) / minuit->GetParameter(i)) * (minuit->GetParameter(i) * minuit->GetParameter(i))
    Double_t error = 0;
    if (TMath::IsNaN(minuit->GetParError(i)))
      Printf("WARNING: Parameter %d error is nan", i);
    else 
      error = 2 * minuit->GetParError(i) * results[i];
    
    if (efficiency)
    {	
      //printf("value before efficiency correction: %f\n",value);
      if (efficiency->GetBinContent(i+1) > 0)
      {
	value /= efficiency->GetBinContent(i+1);
	error /= efficiency->GetBinContent(i+1);
      }
      else
      {
	value = 0;
	error = 0;
      }
    }
    //printf("value after efficiency correction: %f +/- %f\n",value,error);
    result->SetBinContent(i+1, value);
    result->SetBinError(i+1, error);
  }

  Int_t tmpCallCount = fgCallCount;
  fgCallCount = 0; // needs to be 0 so that the Chi2Function prints its output
  Chi2Function(dummy, 0, chi2, results, 0);
  
  Printf("AliUnfolding::UnfoldWithMinuit: iterations %d. Chi2 of final parameters is = %f", tmpCallCount, chi2);
  
  delete[] results;

  return status;
}

//____________________________________________________________________
Int_t AliUnfolding::UnfoldWithBayesian(TH2* correlation, TH1* aEfficiency, TH1* measured, TH1* initialConditions, TH1* aResult)
{
  //
  // unfolds a spectrum using the Bayesian method
  //
  
  if (measured->Integral() <= 0)
  {
    Printf("AliUnfolding::UnfoldWithBayesian: ERROR: The measured spectrum is empty");
    return -1;
  }

  const Int_t kStartBin = 0;

  Int_t kMaxM = fgMaxInput;  //<= fCurrentCorrelation->GetNbinsY(); // max measured axis
  Int_t kMaxT = fgMaxParams; //<= fCurrentCorrelation->GetNbinsX(); // max true axis

  // convergence limit: kMaxT * 0.001^2 = kMaxT * 1e-6 (e.g. 250 bins --> 2.5 e-4)
  const Double_t kConvergenceLimit = kMaxT * 1e-6;

  // store information in arrays, to increase processing speed (~ factor 5)
  Double_t* measuredCopy = new Double_t[kMaxM];
  Double_t* measuredError = new Double_t[kMaxM];
  Double_t* prior = new Double_t[kMaxT];
  Double_t* result = new Double_t[kMaxT];
  Double_t* efficiency = new Double_t[kMaxT];
  Double_t* binWidths = new Double_t[kMaxT];

  Double_t** response = new Double_t*[kMaxT];
  Double_t** inverseResponse = new Double_t*[kMaxT];
  for (Int_t i=0; i<kMaxT; i++)
  {
    response[i] = new Double_t[kMaxM];
    inverseResponse[i] = new Double_t[kMaxM];
  }

  // for normalization
  Float_t measuredIntegral = measured->Integral();
  for (Int_t m=0; m<kMaxM; m++)
  {
    measuredCopy[m] = measured->GetBinContent(m+1) / measuredIntegral;
    measuredError[m] = measured->GetBinError(m+1) / measuredIntegral;

    for (Int_t t=0; t<kMaxT; t++)
    {
      response[t][m] = correlation->GetBinContent(t+1, m+1);
      inverseResponse[t][m] = 0;
    }
  }

  for (Int_t t=0; t<kMaxT; t++)
  {
    if (aEfficiency)
    {
      efficiency[t] = aEfficiency->GetBinContent(t+1);
    }
    else
      efficiency[t] = 1;
      
    prior[t] = measuredCopy[t];
    result[t] = 0;
    binWidths[t] = aResult->GetXaxis()->GetBinWidth(t+1);
  }

  // pick prior distribution
  if (initialConditions)
  {
    printf("Using different starting conditions...\n");
    // for normalization
    Float_t inputDistIntegral = initialConditions->Integral();
    for (Int_t i=0; i<kMaxT; i++)
      prior[i] = initialConditions->GetBinContent(i+1) / inputDistIntegral;
  }

  //TH1F* convergence = new TH1F("convergence", "convergence", 200, 0.5, 200.5);
  
  //new TCanvas;
  // unfold...
  for (Int_t i=0; i<fgBayesianIterations || fgBayesianIterations < 0; i++)
  {
    if (fgDebug)
      Printf("AliUnfolding::UnfoldWithBayesian: iteration %i", i);

    // calculate IR from Bayes theorem
    // IR_ji = R_ij * prior_i / sum_k(R_kj * prior_k)

    Double_t chi2Measured = 0;
    for (Int_t m=0; m<kMaxM; m++)
    {
      Float_t norm = 0;
      for (Int_t t = kStartBin; t<kMaxT; t++)
        norm += response[t][m] * prior[t];

      // calc. chi2: (measured - response * prior) / error
      if (measuredError[m] > 0)
      {
        Double_t value = (measuredCopy[m] - norm) / measuredError[m];
        chi2Measured += value * value;
      }

      if (norm > 0)
      {
        for (Int_t t = kStartBin; t<kMaxT; t++)
          inverseResponse[t][m] = response[t][m] * prior[t] / norm;
      }
      else
      {
        for (Int_t t = kStartBin; t<kMaxT; t++)
          inverseResponse[t][m] = 0;
      }
    }
    //Printf("chi2Measured of the last prior is %e", chi2Measured);

    for (Int_t t = kStartBin; t<kMaxT; t++)
    {
      Float_t value = 0;
      for (Int_t m=0; m<kMaxM; m++)
        value += inverseResponse[t][m] * measuredCopy[m];

      if (efficiency[t] > 0)
        result[t] = value / efficiency[t];
      else
        result[t] = 0;
    }
    
    /* 
    // draw intermediate result
    for (Int_t t=0; t<kMaxT; t++)
    {
      aResult->SetBinContent(t+1, result[t]);
    }
    aResult->SetMarkerStyle(24+i);
    aResult->SetMarkerColor(2);
    aResult->DrawCopy((i == 0) ? "P" : "PSAME");
    */
 
    Double_t chi2LastIter = 0;
    // regularization (simple smoothing)
    for (Int_t t=kStartBin; t<kMaxT; t++)
    {
      Float_t newValue = 0;
      
      // 0 bin excluded from smoothing
      if (t > kStartBin+2 && t<kMaxT-1)
      {
        Float_t average = (result[t-1] / binWidths[t-1] + result[t] / binWidths[t] + result[t+1] / binWidths[t+1]) / 3 * binWidths[t];

        // weight the average with the regularization parameter
        newValue = (1 - fgBayesianSmoothing) * result[t] + fgBayesianSmoothing * average;
      }
      else
        newValue = result[t];

      // calculate chi2 (change from last iteration)
      if (prior[t] > 1e-5)
      {
        Double_t diff = (prior[t] - newValue) / prior[t];
        chi2LastIter += diff * diff;
      }

      prior[t] = newValue;
    }
    //printf("Chi2 of %d iteration = %e\n", i, chi2LastIter);
    //convergence->Fill(i+1, chi2LastIter);

    if (fgBayesianIterations < 0 && chi2LastIter < kConvergenceLimit)
    {
      Printf("AliUnfolding::UnfoldWithBayesian: Stopped Bayesian unfolding after %d iterations at chi2(change since last iteration) of %e; chi2Measured of the last prior is %e", i, chi2LastIter, chi2Measured);
      break;
    }
  } // end of iterations

  //new TCanvas; convergence->DrawCopy(); gPad->SetLogy();
  //delete convergence;

  Float_t factor = 1;
  if (!fgNormalizeInput)
    factor = measuredIntegral;
  for (Int_t t=0; t<kMaxT; t++)
    aResult->SetBinContent(t+1, result[t] * factor);

  delete[] measuredCopy;
  delete[] measuredError;
  delete[] prior;
  delete[] result;
  delete[] efficiency;
  delete[] binWidths;

  for (Int_t i=0; i<kMaxT; i++)
  {
    delete[] response[i];
    delete[] inverseResponse[i];
  }
  delete[] response;
  delete[] inverseResponse;
  
  return 0;

  // ********
  // Calculate the covariance matrix, all arguments are taken from NIM,A362,487-498,1995

  /*printf("Calculating covariance matrix. This may take some time...\n");

  // check if this is the right one...
  TH1* sumHist = GetMultiplicityMC(inputRange, eventType)->ProjectionY("sumHist", 1, GetMultiplicityMC(inputRange, eventType)->GetNbinsX());

  Int_t xBins = hInverseResponseBayes->GetNbinsX();
  Int_t yBins = hInverseResponseBayes->GetNbinsY();

  // calculate "unfolding matrix" Mij
  Float_t matrixM[251][251];
  for (Int_t i=1; i<=xBins; i++)
  {
    for (Int_t j=1; j<=yBins; j++)
    {
      if (fCurrentEfficiency->GetBinContent(i) > 0)
        matrixM[i-1][j-1] = hInverseResponseBayes->GetBinContent(i, j) / fCurrentEfficiency->GetBinContent(i);
      else
        matrixM[i-1][j-1] = 0;
    }
  }

  Float_t* vectorn = new Float_t[yBins];
  for (Int_t j=1; j<=yBins; j++)
    vectorn[j-1] = fCurrentESD->GetBinContent(j);

  // first part of covariance matrix, depends on input distribution n(E)
  Float_t cov1[251][251];

  Float_t nEvents = fCurrentESD->Integral(); // N

  xBins = 20;
  yBins = 20;

  for (Int_t k=0; k<xBins; k++)
  {
    printf("In Cov1: %d\n", k);
    for (Int_t l=0; l<yBins; l++)
    {
      cov1[k][l] = 0;

      // sum_j Mkj Mlj n(Ej) * (1 - n(Ej) / N)
      for (Int_t j=0; j<yBins; j++)
        cov1[k][l] += matrixM[k][j] * matrixM[l][j] * vectorn[j]
          * (1.0 - vectorn[j] / nEvents);

      // - sum_i,j (i != j) Mki Mlj n(Ei) n(Ej) / N
      for (Int_t i=0; i<yBins; i++)
        for (Int_t j=0; j<yBins; j++)
        {
          if (i == j)
            continue;
          cov1[k][l] -= matrixM[k][i] * matrixM[l][j] * vectorn[i]
            * vectorn[j] / nEvents;
         }
    }
  }

  printf("Cov1 finished\n");

  TH2F* cov = (TH2F*) hInverseResponseBayes->Clone("cov");
  cov->Reset();

  for (Int_t i=1; i<=xBins; i++)
    for (Int_t j=1; j<=yBins; j++)
      cov->SetBinContent(i, j, cov1[i-1][j-1]);

  new TCanvas;
  cov->Draw("COLZ");

  // second part of covariance matrix, depends on response matrix
  Float_t cov2[251][251];

  // Cov[P(Er|Cu), P(Es|Cu)] term
  Float_t covTerm[100][100][100];
  for (Int_t r=0; r<yBins; r++)
    for (Int_t u=0; u<xBins; u++)
      for (Int_t s=0; s<yBins; s++)
      {
        if (r == s)
          covTerm[r][u][s] = 1.0 / sumHist->GetBinContent(u+1) * hResponse->GetBinContent(u+1, r+1)
            * (1.0 - hResponse->GetBinContent(u+1, r+1));
        else
          covTerm[r][u][s] = - 1.0 / sumHist->GetBinContent(u+1) * hResponse->GetBinContent(u+1, r+1)
            * hResponse->GetBinContent(u+1, s+1);
      }

  for (Int_t k=0; k<xBins; k++)
    for (Int_t l=0; l<yBins; l++)
    {
      cov2[k][l] = 0;
      printf("In Cov2: %d %d\n", k, l);
      for (Int_t i=0; i<yBins; i++)
        for (Int_t j=0; j<yBins; j++)
        {
          //printf("In Cov2: %d %d %d %d\n", k, l, i, j);
          // calculate Cov(Mki, Mlj) = sum{ru},{su} ...
          Float_t tmpCov = 0;
          for (Int_t r=0; r<yBins; r++)
            for (Int_t u=0; u<xBins; u++)
              for (Int_t s=0; s<yBins; s++)
              {
                if (hResponse->GetBinContent(u+1, r+1) == 0 || hResponse->GetBinContent(u+1, s+1) == 0
                  || hResponse->GetBinContent(u+1, i+1) == 0)
                  continue;

                tmpCov += BayesCovarianceDerivate(matrixM, hResponse, fCurrentEfficiency, k, i, r, u)
                  * BayesCovarianceDerivate(matrixM, hResponse, fCurrentEfficiency, l, j, s, u)
                  * covTerm[r][u][s];
              }

          cov2[k][l] += fCurrentESD->GetBinContent(i+1) * fCurrentESD->GetBinContent(j+1) * tmpCov;
        }
    }

  printf("Cov2 finished\n");

  for (Int_t i=1; i<=xBins; i++)
    for (Int_t j=1; j<=yBins; j++)
      cov->SetBinContent(i, j, cov1[i-1][j-1] + cov2[i-1][j-1]);

  new TCanvas;
  cov->Draw("COLZ");*/
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationPol0(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers constant function (pol0)
  //
  // Does not take into account efficiency
  Double_t chi2 = 0;

  for (Int_t i=1+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    Double_t right  = params[i] / fgUnfoldedAxis->GetBinWidth(i+1);
    Double_t left   = params[i-1] / fgUnfoldedAxis->GetBinWidth(i);

    if (left != 0)
    {
      Double_t diff = (right - left);
      chi2 += diff * diff / left / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
    }
  }

  return chi2 / 100.0;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationPol1(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers linear function (pol1)
  //
  // Does not take into account efficiency
  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0)
      continue;

    Double_t right  = params[i] / fgUnfoldedAxis->GetBinWidth(i+1);
    Double_t middle = params[i-1] / fgUnfoldedAxis->GetBinWidth(i);
    Double_t left   = params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1);

    Double_t der1 = (right - middle) / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
    Double_t der2 = (middle - left) / ((fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i-1)) / 2);

    //Double_t diff = (der1 - der2) / middle;
    //chi2 += diff * diff;
    chi2 += (der1 - der2) * (der1 - der2) / middle * fgUnfoldedAxis->GetBinWidth(i);
  }

  return chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationLog(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers logarithmic function (log)
  //
  // Does not take into account efficiency

  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0 || params[i] == 0 || params[i-2] == 0)
     continue;

    Double_t right  = log(params[i] / fgUnfoldedAxis->GetBinWidth(i+1));
    Double_t middle = log(params[i-1] / fgUnfoldedAxis->GetBinWidth(i));
    Double_t left   = log(params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1));
    
    Double_t der1 = (right - middle) / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
    Double_t der2 = (middle - left) / ((fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i-1)) / 2);
    
    //Double_t error = 1. / params[i] + 4. / params[i-1] + 1. / params[i-2];

    //if (fgCallCount == 0)
    //  Printf("%d %f %f", i, (der1 - der2) * (der1 - der2), error);
    chi2 += (der1 - der2) * (der1 - der2);// / error;
  }

  return chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationTotalCurvature(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // minimizes the total curvature (from Unfolding Methods In High-Energy Physics Experiments,
  // V. Blobel (Hamburg U.) . DESY 84/118, Dec 1984. 40pp.
  //
  // Does not take into account efficiency

  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    Double_t right  = params[i];
    Double_t middle = params[i-1];
    Double_t left   = params[i-2];

    Double_t der1 = (right - middle);
    Double_t der2 = (middle - left);

    Double_t diff = (der1 - der2);

    chi2 += diff * diff;
  }

  return chi2 * 1e4;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationEntropy(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // calculates entropy, from
  // The method of reduced cross-entropy (M. Schmelling 1993)
  //
  // Does not take into account efficiency

  Double_t paramSum = 0;
  
  for (Int_t i=fgSkipBinsBegin; i<fgMaxParams; ++i)
    paramSum += params[i];

  Double_t chi2 = 0;
  for (Int_t i=fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    Double_t tmp = params[i] / paramSum;
    //Double_t tmp = params[i];
    if (tmp > 0 && (*fgEntropyAPriori)[i] > 0)
    {
      chi2 += tmp * TMath::Log(tmp / (*fgEntropyAPriori)[i]);
    }
    else
      chi2 += 100;
  }

  return -chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationRatio(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  //
  // Does not take into account efficiency

  Double_t chi2 = 0;

  for (Int_t i=5+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0 || params[i] == 0)
      continue;

    Double_t right  = params[i] / fgUnfoldedAxis->GetBinWidth(i+1);
    Double_t middle = params[i-1] / fgUnfoldedAxis->GetBinWidth(i);
    Double_t left   = params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1);
    Double_t left2   = params[i-3] / fgUnfoldedAxis->GetBinWidth(i-2);
    Double_t left3   = params[i-4] / fgUnfoldedAxis->GetBinWidth(i-3);
    Double_t left4   = params[i-5] / fgUnfoldedAxis->GetBinWidth(i-4);

    //Double_t diff = left / middle - middle / right;
    //Double_t diff = 2 * left / middle - middle / right - left2 / left;
    Double_t diff = 4 * left2 / left - middle / right - left / middle - left3 / left2 - left4 / left3;
    
    chi2 += diff * diff;// / middle;
  }

  return chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationPowerLaw(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers power law with n = -5
  //
  // Does not take into account efficiency

  Double_t chi2 = 0;

  Double_t right = 0.;
  Double_t middle = 0.;
  Double_t left = 0.;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i] < 1e-8 || params[i-1] < 1e-8 || params[i-2] < 1e-8)
      continue;

    if (fgUnfoldedAxis->GetBinWidth(i+1) < 1e-8 || fgUnfoldedAxis->GetBinWidth(i) < 1e-8 || fgUnfoldedAxis->GetBinWidth(i-1) < 1e-8)
      continue;
    
    middle = TMath::Power(params[i-1] / fgUnfoldedAxis->GetBinWidth(i),fgPowern);

    if(middle>0) {
      right  = TMath::Power(params[i] / fgUnfoldedAxis->GetBinWidth(i),fgPowern)/middle;

      left   = TMath::Power(params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1),fgPowern)/middle;
      
      middle = 1.;
      
      Double_t der1 = (right - middle) / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
      Double_t der2 = (middle - left) / ((fgUnfoldedAxis->GetBinWidth(i-1) + fgUnfoldedAxis->GetBinWidth(i-2)) / 2);
    
      chi2 += (der1 - der2) * (der1 - der2)/ ( fgUnfoldedAxis->GetBinWidth(i)/2. + fgUnfoldedAxis->GetBinWidth(i-1) + fgUnfoldedAxis->GetBinWidth(i-2)/2.)/( fgUnfoldedAxis->GetBinWidth(i)/2. + fgUnfoldedAxis->GetBinWidth(i-1) + fgUnfoldedAxis->GetBinWidth(i-2)/2.);// / error;
      //   printf("i: %d chi2 = %f\n",i,chi2);
    }

  }

  return chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationLogLog(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers a powerlaw (linear on a log-log scale)
  //
  // The calculation takes into account the efficiencies

  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0 || params[i] == 0 || params[i-2] == 0)
     continue;
    if ((*fgEfficiency)(i-1) == 0 || (*fgEfficiency)(i) == 0 || (*fgEfficiency)(i-2) == 0)
     continue;


    Double_t right  = log(params[i] / (*fgEfficiency)(i) / fgUnfoldedAxis->GetBinWidth(i));
    Double_t middle = log(params[i-1] / (*fgEfficiency)(i-1) / fgUnfoldedAxis->GetBinWidth(i-1));
    Double_t left   = log(params[i-2] / (*fgEfficiency)(i-2) / fgUnfoldedAxis->GetBinWidth(i-2));
    
    Double_t der1 = (right - middle) / ( log(fgUnfoldedAxis->GetBinCenter(i+1)) - log(fgUnfoldedAxis->GetBinCenter(i)) );
    Double_t der2 = (middle - left) /( log(fgUnfoldedAxis->GetBinCenter(i)) - log(fgUnfoldedAxis->GetBinCenter(i-1)) );

    double tmp = (log(fgUnfoldedAxis->GetBinCenter(i+1)) - log(fgUnfoldedAxis->GetBinCenter(i-1)))/2.;
    Double_t dder = (der1-der2) / tmp;

    chi2 += dder * dder;
  }

  return chi2;
}



//____________________________________________________________________
void AliUnfolding::Chi2Function(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t)
{
  //
  // fit function for minuit
  // does: (m - Ad)W(m - Ad) where m = measured, A correlation matrix, d = guess, W = covariance matrix
  //
  
  // TODO use static members for the variables here to speed up processing (no construction/deconstruction)

  // d = guess
  TVectorD paramsVector(fgMaxParams);
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramsVector[i] = params[i] * params[i];

  // calculate penalty factor
  Double_t penaltyVal = 0;

  switch (fgRegularizationType)
  {
    case kNone:       break;
    case kPol0:       penaltyVal = RegularizationPol0(paramsVector); break;
    case kPol1:       penaltyVal = RegularizationPol1(paramsVector); break;
    case kCurvature:  penaltyVal = RegularizationTotalCurvature(paramsVector); break;
    case kEntropy:    penaltyVal = RegularizationEntropy(paramsVector); break;
    case kLog:        penaltyVal = RegularizationLog(paramsVector); break;
    case kRatio:      penaltyVal = RegularizationRatio(paramsVector); break;
    case kPowerLaw:   penaltyVal = RegularizationPowerLaw(paramsVector); break;
    case kLogLog:     penaltyVal = RegularizationLogLog(paramsVector); break;
  }

  penaltyVal *= fgRegularizationWeight; //beta*PU

  // Ad
  TVectorD measGuessVector(fgMaxInput);
  measGuessVector = (*fgCorrelationMatrix) * paramsVector;

  // Ad - m
  measGuessVector -= (*fgCurrentESDVector);
  
#if 0
  // new error calcuation using error on the guess
  
  // error from guess
  TVectorD errorGuessVector(fgMaxInput);
  //errorGuessVector = (*fgCorrelationMatrixSquared) * paramsVector;
  errorGuessVector = (*fgCorrelationMatrix) * paramsVector;

  Double_t chi2FromFit = 0;
  for (Int_t i=0; i<fgMaxInput; ++i)
  {
    Float_t error = 1;
    if (errorGuessVector(i) > 0)
      error = errorGuessVector(i);
    chi2FromFit += measGuessVector(i) * measGuessVector(i) / error;
  }

#else
  // old error calcuation using the error on the measured
  TVectorD copy(measGuessVector);

  // (Ad - m) W
  // this step can be optimized because currently only the diagonal elements of fgCorrelationCovarianceMatrix are used
  // normal way is like this:
  // measGuessVector *= (*fgCorrelationCovarianceMatrix);
  // optimized way like this:
  for (Int_t i=0; i<fgMaxInput; ++i)
    measGuessVector[i] *= (*fgCorrelationCovarianceMatrix)(i, i);


  if (fgSkipBin0InChi2)
    measGuessVector[0] = 0;

  // (Ad - m) W (Ad - m)
  // the factor 1e6 prevents that a very small number (measGuessVector[i]) is multiplied with a very
  // big number ((*fgCorrelationCovarianceMatrix)(i, i)) (see UnfoldWithMinuit)
  Double_t chi2FromFit = measGuessVector * copy * 1e6;
#endif

  Double_t notFoundEventsConstraint = 0;
  Double_t currentNotFoundEvents = 0;
  Double_t errorNotFoundEvents = 0;
  
  if (fgNotFoundEvents > 0)
  {
    for (Int_t i=0; i<fgMaxParams; ++i)
    {
      Double_t eff = (1.0 / (*fgEfficiency)(i) - 1);
      if (i == 0)
	eff = (1.0 / params[fgMaxParams] - 1);
      currentNotFoundEvents += eff * paramsVector(i);
      errorNotFoundEvents += eff * eff * paramsVector(i); // error due to guess (paramsVector)
      if (i <= 3)
	errorNotFoundEvents += (eff * 0.03) * (eff * 0.03) * paramsVector(i) * paramsVector(i); // error on eff
      //      if ((fgCallCount % 10000) == 0)
	//Printf("%d %f %f %f", i, (*fgEfficiency)(i), paramsVector(i), currentNotFoundEvents);
    }
    errorNotFoundEvents += fgNotFoundEvents;
    // TODO add error on background, background estimate
    
    notFoundEventsConstraint = (currentNotFoundEvents - fgNotFoundEvents) * (currentNotFoundEvents - fgNotFoundEvents) / errorNotFoundEvents;
  }
  
  Float_t avg = 0;
  Float_t sum = 0;
  Float_t currentMult = 0;
  // try to match dn/deta
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    avg += paramsVector(i) * currentMult;
    sum += paramsVector(i);
    currentMult += fgUnfoldedAxis->GetBinWidth(i);
  }
  avg /= sum;
  Float_t chi2Avg = 0; //(avg - 3.73) * (avg - 3.73) * 100;

  chi2 = chi2FromFit + penaltyVal + notFoundEventsConstraint + chi2Avg;

  if ((fgCallCount++ % 1000) == 0)
  {

    Printf("AliUnfolding::Chi2Function: Iteration %d (ev %d %d +- %f) (%f) (%f): %f %f %f %f --> %f", fgCallCount-1, (Int_t) fgNotFoundEvents, (Int_t) currentNotFoundEvents, TMath::Sqrt(errorNotFoundEvents), params[fgMaxParams-1], avg, chi2FromFit, penaltyVal, notFoundEventsConstraint, chi2Avg, chi2);

    //for (Int_t i=0; i<fgMaxInput; ++i)
    //  Printf("%d: %f", i, measGuessVector(i) * copy(i) * 1e6);
  }

  fChi2FromFit = chi2FromFit;
  fPenaltyVal = penaltyVal;
}

//____________________________________________________________________
void AliUnfolding::MakePads() {
    TPad *presult = new TPad("presult","result",0,0.4,1,1);
    presult->SetNumber(1);
    presult->Draw();
    presult->SetLogy();
    TPad *pres = new TPad("pres","residuals",0,0.2,1,0.4);
    pres->SetNumber(2);
    pres->Draw();
    TPad *ppen = new TPad("ppen","penalty",0,0.0,1,0.2);
    ppen->SetNumber(3);
    ppen->Draw();
 
}
//____________________________________________________________________
void AliUnfolding::DrawResults(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TCanvas *canv, Int_t reuseHists,TH1 *unfolded)
{
  // Draw histograms of
  //   - Result folded with response
  //   - Penalty factors
  //   - Chisquare contributions
  // (Useful for debugging/sanity checks and the interactive unfolder)
  //
  // If a canvas pointer is given (canv != 0), it will be used for all
  // plots; 3 pads are made if needed.


  Int_t blankCanvas = 0;
  TVirtualPad *presult = 0;
  TVirtualPad *pres = 0; 
  TVirtualPad *ppen = 0;
  
  if (canv) {
    canv->cd();
    presult = canv->GetPad(1);
    pres = canv->GetPad(2);
    ppen = canv->GetPad(3); 
    if (presult == 0 || pres == 0 || ppen == 0) {
      canv->Clear();
      MakePads();
      blankCanvas = 1;
      presult = canv->GetPad(1);
      pres = canv->GetPad(2);
      ppen = canv->GetPad(3); 
    } 
  }
  else {
    presult = new TCanvas;
    pres = new TCanvas;
    ppen = new TCanvas;
  }


  if (fgMaxInput == -1)
  {
    Printf("AliUnfolding::Unfold: WARNING. Number of measured bins not set with SetNbins. Using number of bins in measured distribution");
    fgMaxInput = measured->GetNbinsX();
  }
  if (fgMaxParams == -1)
  {
    Printf("AliUnfolding::Unfold: WARNING. Number of unfolded bins not set with SetNbins. Using number of bins in measured distribution");
    fgMaxParams = initialConditions->GetNbinsX();
  }

  if (fgOverflowBinLimit > 0)
    CreateOverflowBin(correlation, measured);

  // Uses Minuit methods

  SetStaticVariables(correlation, measured, efficiency);

  Double_t* params = new Double_t[fgMaxParams+1];
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    params[i] = initialConditions->GetBinContent(i+1) * efficiency->GetBinContent(i+1);

    Bool_t fix = kFALSE;
    if (params[i] < 0)
    {
      fix = kTRUE;
      params[i] = -params[i];
    }
    params[i]=TMath::Sqrt(params[i]);

    //cout << "params[" << i << "] " << params[i] << endl;

  } 

  Double_t chi2;
  Int_t dummy;

  //fgPrintChi2Details = kTRUE;
  fgCallCount = 0; // To make sure that Chi2Function prints the components
  Chi2Function(dummy, 0, chi2, params, 0);

  presult->cd();
  TH1 *meas2 = (TH1*)measured->Clone();
  meas2->SetLineColor(1);
  meas2->SetLineWidth(2);
  meas2->SetMarkerColor(meas2->GetLineColor());
  meas2->SetMarkerStyle(20);
  Float_t scale = unfolded->GetXaxis()->GetBinWidth(1)/meas2->GetXaxis()->GetBinWidth(1);
  meas2->Scale(scale);
  if (blankCanvas)
    meas2->DrawCopy();
  else
    meas2->DrawCopy("same");
  //meas2->GetXaxis()->SetLimits(0,fgMaxInput);
  meas2->SetBit(kCannotPick);
  DrawGuess(params, presult, pres, ppen, reuseHists,unfolded);
  delete [] params;
}
//____________________________________________________________________
void AliUnfolding::RedrawInteractive() {
  // 
  // Helper function for interactive unfolding
  //
  DrawResults(fghCorrelation,fghEfficiency,fghMeasured,fghUnfolded,fgCanvas,1,fghUnfolded);
}
//____________________________________________________________________
void AliUnfolding::InteractiveUnfold(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions) 
{
  //
  // Function to do interactive unfolding
  // A canvas is drawn with the unfolding result
  // Change the histogram with your mouse and all histograms 
  // will be updated automatically

  fgCanvas = new TCanvas("UnfoldingCanvas","Interactive unfolding",500,800);
  fgCanvas->cd();

  MakePads();

  if (fghUnfolded) {
    delete fghUnfolded; fghUnfolded = 0;
  }
  
  gROOT->SetEditHistograms(kTRUE);

  fghUnfolded = new TH1F("AluUnfoldingfghUnfolded","Unfolded result (Interactive unfolder",efficiency->GetNbinsX(),efficiency->GetXaxis()->GetXmin(),efficiency->GetXaxis()->GetXmax());
  fghUnfolded->SetLineColor(2);
  fghUnfolded->SetMarkerColor(2);
  fghUnfolded->SetLineWidth(2);

  fghCorrelation = correlation;
  fghEfficiency = efficiency;
  fghMeasured = measured;

  UnfoldWithMinuit(correlation,efficiency,measured,initialConditions,fghUnfolded,kFALSE);

  fgCanvas->cd(1);
  fghUnfolded->Draw();
  DrawResults(correlation,efficiency,measured,fghUnfolded,fgCanvas,kFALSE,fghUnfolded);

  TExec *execRedraw = new TExec("redraw","AliUnfolding::RedrawInteractive()");
  fghUnfolded->GetListOfFunctions()->Add(execRedraw);
}
//____________________________________________________________________
void AliUnfolding::DrawGuess(Double_t *params, TVirtualPad *pfolded, TVirtualPad *pres, TVirtualPad *ppen, Int_t reuseHists,TH1* unfolded)
{
  //
  // draws residuals of solution suggested by params and effect of regularization
  //

  if (pfolded == 0)
    pfolded = new TCanvas;
  if (ppen == 0)
    ppen = new TCanvas;
  if (pres == 0)
    pres = new TCanvas;
  
  // d
  TVectorD paramsVector(fgMaxParams);
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramsVector[i] = params[i] * params[i];

  // Ad
  TVectorD measGuessVector(fgMaxInput);
  measGuessVector = (*fgCorrelationMatrix) * paramsVector;

  TH1* folded = 0;
  if (reuseHists) 
    folded = dynamic_cast<TH1*>(gROOT->FindObject("__hfolded"));
  if (!reuseHists || folded == 0) {
    if (fgMeasuredAxis->GetXbins()->GetArray()) // variable bins
      folded = new TH1F("__hfolded","Folded histo from AliUnfolding",fgMeasuredAxis->GetNbins(),fgMeasuredAxis->GetXbins()->GetArray());
    else
      folded = new TH1F("__hfolded","Folded histo from AliUnfolding",fgMeasuredAxis->GetNbins(),fgMeasuredAxis->GetXmin(),fgMeasuredAxis->GetXmax());
  }

  folded->SetBit(kCannotPick);
  folded->SetLineColor(4);
  folded->SetLineWidth(2);

  for (Int_t ibin =0; ibin < fgMaxInput; ibin++) 
    folded->SetBinContent(ibin+1, measGuessVector[ibin]);

  Float_t scale = unfolded->GetXaxis()->GetBinWidth(1)/folded->GetXaxis()->GetBinWidth(1);
  folded->Scale(scale);

  pfolded->cd();

  folded->Draw("same");

  // Ad - m
  measGuessVector -= (*fgCurrentESDVector);

  TVectorD copy(measGuessVector);

  // (Ad - m) W
  // this step can be optimized because currently only the diagonal elements of fgCorrelationCovarianceMatrix are used
  // normal way is like this:
  // measGuessVector *= (*fgCorrelationCovarianceMatrix);
  // optimized way like this:
  for (Int_t i=0; i<fgMaxInput; ++i)
    measGuessVector[i] *= (*fgCorrelationCovarianceMatrix)(i, i);

  // (Ad - m) W (Ad - m)
  // the factor 1e6 prevents that a very small number (measGuessVector[i]) is multiplied with a very
  // big number ((*fgCorrelationCovarianceMatrix)(i, i)) (see ApplyMinuitFit)
  //Double_t chi2FromFit = measGuessVector * copy * 1e6;

  // draw residuals
  // Double_t pTarray[fgMaxParams+1];
  // for(int i=0; i<fgMaxParams; i++) {
  //   pTarray[i] = fgUnfoldedAxis->GetBinCenter(i)-0.5*fgUnfoldedAxis->GetBinWidth(i);
  // }
  //  pTarray[fgMaxParams] = fgUnfoldedAxis->GetBinCenter(fgMaxParams-1)+0.5*fgUnfoldedAxis->GetBinWidth(fgMaxParams-1);
  //  TH1* residuals = new TH1F("residuals", "residuals", fgMaxParams,pTarray);
  //  TH1* residuals = new TH1F("residuals", "residuals", fgMaxInput, -0.5, fgMaxInput - 0.5);
  // for (Int_t i=0; i<fgMaxInput; ++i)
  //   residuals->SetBinContent(i+1, measGuessVector(i) * copy(i) * 1e6);7
  TH1* residuals = GetResidualsPlot(params);
  residuals->GetXaxis()->SetTitleSize(0.06);
  residuals->GetXaxis()->SetTitleOffset(0.7);
  residuals->GetXaxis()->SetLabelSize(0.07);
  residuals->GetYaxis()->SetTitleSize(0.08);
  residuals->GetYaxis()->SetTitleOffset(0.5);
  residuals->GetYaxis()->SetLabelSize(0.06);
  pres->cd(); residuals->DrawCopy(); gPad->SetLogy();
    

  // draw penalty
  TH1* penalty = GetPenaltyPlot(params);
  penalty->GetXaxis()->SetTitleSize(0.06);
  penalty->GetXaxis()->SetTitleOffset(0.7);
  penalty->GetXaxis()->SetLabelSize(0.07);
  penalty->GetYaxis()->SetTitleSize(0.08);
  penalty->GetYaxis()->SetTitleOffset(0.5);
  penalty->GetYaxis()->SetLabelSize(0.06);
  ppen->cd(); penalty->DrawCopy(); gPad->SetLogy();
}

//____________________________________________________________________
TH1* AliUnfolding::GetResidualsPlot(TH1* corrected)
{
  //
  // MvL: THIS MUST BE INCORRECT. 
  // Need heff to calculate params from TH1 'corrected'
  //

  //
  // fill residuals histogram of solution suggested by params and effect of regularization
  //

  Double_t* params = new Double_t[fgMaxParams];
  for (Int_t i=0; i<TMath::Min(fgMaxParams, corrected->GetNbinsX()); i++)
    params[i] = TMath::Sqrt(TMath::Abs(corrected->GetBinContent(i+1)*(*fgEfficiency)(i)));


  TH1 * plot = GetResidualsPlot(params);
  delete [] params;
  return plot;
}

//____________________________________________________________________
TH1* AliUnfolding::GetResidualsPlot(Double_t* params)
{
  //
  // fill residuals histogram of solution suggested by params and effect of regularization
  //

  // d
  TVectorD paramsVector(fgMaxParams);
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramsVector[i] = params[i] * params[i];

  // Ad
  TVectorD measGuessVector(fgMaxInput);
  measGuessVector = (*fgCorrelationMatrix) * paramsVector;

  // Ad - m
  measGuessVector -= (*fgCurrentESDVector);

  TVectorD copy(measGuessVector);

  // (Ad - m) W
  // this step can be optimized because currently only the diagonal elements of fgCorrelationCovarianceMatrix are used
  // normal way is like this:
  // measGuessVector *= (*fgCorrelationCovarianceMatrix);
  // optimized way like this:
  for (Int_t i=0; i<fgMaxInput; ++i)
    measGuessVector[i] *= (*fgCorrelationCovarianceMatrix)(i, i);

  // (Ad - m) W (Ad - m)
  // the factor 1e6 prevents that a very small number (measGuessVector[i]) is multiplied with a very
  // big number ((*fgCorrelationCovarianceMatrix)(i, i)) (see ApplyMinuitFit)
  //Double_t chi2FromFit = measGuessVector * copy * 1e6;

  // draw residuals
  TH1* residuals = 0;
  if (fgMeasuredAxis->GetXbins()->GetArray()) // variable bins
    residuals = new TH1F("residuals", "residuals;unfolded pos;residual",fgMeasuredAxis->GetNbins(),fgMeasuredAxis->GetXbins()->GetArray());
  else
    residuals = new TH1F("residuals", "residuals;unfolded pos;residual",fgMeasuredAxis->GetNbins(),fgMeasuredAxis->GetXmin(), fgMeasuredAxis->GetXmax());
  //  TH1* residuals = new TH1F("residuals", "residuals", fgMaxInput, -0.5, fgMaxInput - 0.5);

  Double_t sumResiduals = 0.; 
  for (Int_t i=0; i<fgMaxInput; ++i) {
    residuals->SetBinContent(i+1, measGuessVector(i) * copy(i) * 1e6);
    sumResiduals += measGuessVector(i) * copy(i) * 1e6;
  }
  fAvgResidual = sumResiduals/(double)fgMaxInput;
 
  return residuals;
}

//____________________________________________________________________
TH1* AliUnfolding::GetPenaltyPlot(TH1* corrected)
{
  // draws the penalty factors as function of multiplicity of the current selected regularization

  Double_t* params = new Double_t[fgMaxParams];
  for (Int_t i=0; i<TMath::Min(fgMaxParams, corrected->GetNbinsX()); i++)
    params[i] = (*fgEfficiency)(i)*corrected->GetBinContent(i+1);
  
  TH1* penalty = GetPenaltyPlot(params);
  
  delete[] params;
  
  return penalty;
}

//____________________________________________________________________
TH1* AliUnfolding::GetPenaltyPlot(Double_t* params)
{
  // draws the penalty factors as function of multiplicity of the current selected regularization
  
  //TH1* penalty = new TH1F("penalty", ";unfolded multiplicity;penalty factor", fgMaxParams, -0.5, fgMaxParams - 0.5);
  //  TH1* penalty = new TH1F("penalty", ";unfolded pos;penalty factor", fgMaxParams, fgUnfoldedAxis->GetBinCenter(0)-0.5*fgUnfoldedAxis->GetBinWidth(0),fgUnfoldedAxis->GetBinCenter(fgMaxParams)+0.5*fgUnfoldedAxis->GetBinWidth(fgMaxParams) );

  TH1* penalty = 0;
  if (fgUnfoldedAxis->GetXbins()->GetArray())
    penalty = new TH1F("penalty", ";unfolded pos;penalty factor", fgUnfoldedAxis->GetNbins(),fgUnfoldedAxis->GetXbins()->GetArray());
  else
    penalty = new TH1F("penalty", ";unfolded pos;penalty factor", fgUnfoldedAxis->GetNbins(),fgUnfoldedAxis->GetXmin(),fgUnfoldedAxis->GetXmax());

  for (Int_t i=1+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    Double_t diff = 0;
    if (fgRegularizationType == kPol0)
    {
      Double_t right  = params[i] / fgUnfoldedAxis->GetBinWidth(i+1);
      Double_t left   = params[i-1] / fgUnfoldedAxis->GetBinWidth(i);

      if (left != 0)
      {
        Double_t diffTmp = (right - left);
        diff = diffTmp * diffTmp / left / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2) / 100;
      }
    } 
    if (fgRegularizationType == kPol1 && i > 1+fgSkipBinsBegin) 
    {
      if (params[i-1] == 0)
        continue;

      Double_t right  = params[i] / fgUnfoldedAxis->GetBinWidth(i+1);
      Double_t middle = params[i-1] / fgUnfoldedAxis->GetBinWidth(i);
      Double_t left   = params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1);

      Double_t der1 = (right - middle) / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
      Double_t der2 = (middle - left) / ((fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i-1)) / 2);

      diff = (der1 - der2) * (der1 - der2) / middle;
    }

    if (fgRegularizationType == kLog && i > 1+fgSkipBinsBegin) 
    {
      if (params[i-1] == 0)
        continue;
  
      Double_t right  = log(params[i]);
      Double_t middle = log(params[i-1]);
      Double_t left   = log(params[i-2]);
      
      Double_t der1 = (right - middle);
      Double_t der2 = (middle - left);
  
      //Double_t error = 1. / params[i] + 4. / params[i-1] + 1. / params[i-2];
      //Printf("%d %f %f", i, (der1 - der2) * (der1 - der2), error);
       
      diff = (der1 - der2) * (der1 - der2);// / error;
    }
    if (fgRegularizationType == kCurvature && i > 1+fgSkipBinsBegin)
    {
      Double_t right  = params[i];    // params are sqrt
      Double_t middle = params[i-1];
      Double_t left   = params[i-2];
      
      Double_t der1 = (right - middle)/0.5/(fgUnfoldedAxis->GetBinWidth(i-1) + fgUnfoldedAxis->GetBinWidth(i));
      Double_t der2 = (middle - left)/0.5/(fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i+1));
      
      diff = (der1 - der2)/(fgUnfoldedAxis->GetBinWidth(i-1) + fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i-1))*3.0;
      diff = 1e4*diff*diff;
    }
    if (fgRegularizationType == kPowerLaw && i > 1+fgSkipBinsBegin) 
    {

      if (params[i] < 1e-8 || params[i-1] < 1e-8 || params[i-2] < 1e-8)
	continue;
      
      if (fgUnfoldedAxis->GetBinWidth(i+1) < 1e-8 || fgUnfoldedAxis->GetBinWidth(i) < 1e-8 || fgUnfoldedAxis->GetBinWidth(i) < 1e-8)
	continue;
      
      double middle = TMath::Power(params[i-1] / fgUnfoldedAxis->GetBinWidth(i),fgPowern);

      if(middle>0) {
	double right  = TMath::Power(params[i] / fgUnfoldedAxis->GetBinWidth(i+1),fgPowern)/middle;
	
	double left   = TMath::Power(params[i-2] / fgUnfoldedAxis->GetBinWidth(i-1),fgPowern)/middle;
	
	middle = 1.;
	
	Double_t der1 = (right - middle) / ((fgUnfoldedAxis->GetBinWidth(i+1) + fgUnfoldedAxis->GetBinWidth(i)) / 2);
	Double_t der2 = (middle - left) / ((fgUnfoldedAxis->GetBinWidth(i) + fgUnfoldedAxis->GetBinWidth(i-1)) / 2);
	
	diff = (der1 - der2) * (der1 - der2);// / error;
      }
    }

    if (fgRegularizationType == kLogLog && i > 1+fgSkipBinsBegin) 
    {

      if (params[i] < 1e-8 || params[i-1] < 1e-8 || params[i-2] < 1e-8)
	continue;

      Double_t right  = log(params[i] / (*fgEfficiency)(i) / fgUnfoldedAxis->GetBinWidth(i+1));
      Double_t middle = log(params[i-1] / (*fgEfficiency)(i-1) / fgUnfoldedAxis->GetBinWidth(i));
      Double_t left   = log(params[i-2] / (*fgEfficiency)(i-2) / fgUnfoldedAxis->GetBinWidth(i-1));
      
      Double_t der1 = (right - middle) / ( log(fgUnfoldedAxis->GetBinCenter(i+1)) - log(fgUnfoldedAxis->GetBinCenter(i)) );
      Double_t der2 = (middle - left) /( log(fgUnfoldedAxis->GetBinCenter(i)) - log(fgUnfoldedAxis->GetBinCenter(i-1)) );
      
      double tmp = (log(fgUnfoldedAxis->GetBinCenter(i+1)) - log(fgUnfoldedAxis->GetBinCenter(i-1)))/2.;
      Double_t dder = (der1-der2) / tmp;
      
      diff = dder * dder;
    }
    
    penalty->SetBinContent(i, diff*fgRegularizationWeight);
    
    //Printf("%d %f %f %f %f", i-1, left, middle, right, diff);
  }
  
  return penalty;
}

//____________________________________________________________________
void AliUnfolding::TF1Function(Int_t& unused1, Double_t* unused2, Double_t& chi2, Double_t *params, Int_t unused3)
{
  //
  // fit function for minuit
  // uses the TF1 stored in fgFitFunction
  //

  for (Int_t i=0; i<fgFitFunction->GetNpar(); i++)
    fgFitFunction->SetParameter(i, params[i]);

  Double_t* params2 = new Double_t[fgMaxParams];

  for (Int_t i=0; i<fgMaxParams; ++i)
    params2[i] = fgFitFunction->Eval(i);

  Chi2Function(unused1, unused2, chi2, params2, unused3);
  
  delete[] params2;

  if (fgDebug)
    Printf("%f", chi2);
}

//____________________________________________________________________
Int_t AliUnfolding::UnfoldWithFunction(TH2* correlation, TH1* efficiency, TH1* measured, TH1* /* initialConditions */, TH1* aResult)
{
  //
  // correct spectrum using minuit chi2 method applying a functional fit
  //
  
  if (!fgFitFunction)
  {
    Printf("AliUnfolding::UnfoldWithFunction: ERROR fit function not set. Exiting.");
    return -1;
  }    
  
  SetChi2Regularization(kNone, 0);
  
  SetStaticVariables(correlation, measured, efficiency);
  
  // Initialize TMinuit via generic fitter interface
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, fgFitFunction->GetNpar());

  minuit->SetFCN(TF1Function);
  for (Int_t i=0; i<fgFitFunction->GetNpar(); i++)
  {
    Double_t lower, upper;
    fgFitFunction->GetParLimits(i, lower, upper);
    minuit->SetParameter(i, Form("param%d",i), fgFitFunction->GetParameter(i), fgMinuitStepSize, lower, upper);
  }

  Double_t arglist[100];
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  minuit->ExecuteCommand("SCAN", arglist, 0);
  minuit->ExecuteCommand("MIGRAD", arglist, 0);
  //minuit->ExecuteCommand("MINOS", arglist, 0);

  for (Int_t i=0; i<fgFitFunction->GetNpar(); i++)
    fgFitFunction->SetParameter(i, minuit->GetParameter(i));

  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    Double_t value = fgFitFunction->Eval(i);
    if (fgDebug)
      Printf("%d : %f", i, value);
    
    if (efficiency)
    {
      if (efficiency->GetBinContent(i+1) > 0)
      {
        value /= efficiency->GetBinContent(i+1);
      }
      else
        value = 0;
    }
    aResult->SetBinContent(i+1, value);
    aResult->SetBinError(i+1, 0);
  }
  
  return 0;
}

//____________________________________________________________________
void AliUnfolding::CreateOverflowBin(TH2* correlation, TH1* measured)
{
  // Finds the first bin where the content is below fgStatLimit and combines all values for this bin and larger bins
  // The same limit is applied to the correlation  
  
  Int_t lastBin = 0;
  for (Int_t i=1; i<=measured->GetNbinsX(); ++i)
  {
    if (measured->GetBinContent(i) <= fgOverflowBinLimit)
    {
      lastBin = i;
      break;
    }
  }
  
  if (lastBin == 0)
  {
    Printf("AliUnfolding::CreateOverflowBin: WARNING: First bin is already below limit of %f", fgOverflowBinLimit);
    return;
  }
  
  Printf("AliUnfolding::CreateOverflowBin: Bin limit in measured spectrum determined to be %d", lastBin);
  
  TCanvas* canvas = 0;

  if (fgDebug)
  {
    canvas = new TCanvas("StatSolution", "StatSolution", 1000, 800);
    canvas->Divide(2, 2);

    canvas->cd(1);
    measured->SetStats(kFALSE);
    measured->DrawCopy();
    gPad->SetLogy();

    canvas->cd(2);
    correlation->SetStats(kFALSE);
    correlation->DrawCopy("COLZ");
  }

  measured->SetBinContent(lastBin, measured->Integral(lastBin, measured->GetNbinsX()));
  for (Int_t i=lastBin+1; i<=measured->GetNbinsX(); ++i)
  {
    measured->SetBinContent(i, 0);
    measured->SetBinError(i, 0);
  }
  // the error is set to sqrt(N), better solution possible?, sum of relative errors of all contributions???
  measured->SetBinError(lastBin, TMath::Sqrt(measured->GetBinContent(lastBin)));

  Printf("AliUnfolding::CreateOverflowBin: This bin has now %f +- %f entries", measured->GetBinContent(lastBin), measured->GetBinError(lastBin));

  for (Int_t i=1; i<=correlation->GetNbinsX(); ++i)
  {
    correlation->SetBinContent(i, lastBin, correlation->Integral(i, i, lastBin, correlation->GetNbinsY()));
    // the error is set to sqrt(N), better solution possible?, sum of relative errors of all contributions???
    correlation->SetBinError(i, lastBin, TMath::Sqrt(correlation->GetBinContent(i, lastBin)));

    for (Int_t j=lastBin+1; j<=correlation->GetNbinsY(); ++j)
    {
      correlation->SetBinContent(i, j, 0);
      correlation->SetBinError(i, j, 0);
    }
  }

  Printf("AliUnfolding::CreateOverflowBin: Adjusted correlation matrix!");

  if (canvas)
  {
    canvas->cd(3);
    measured->DrawCopy();
    gPad->SetLogy();

    canvas->cd(4);
    correlation->DrawCopy("COLZ");
  }
}

Int_t AliUnfolding::UnfoldGetBias(TH2* correlation, TH1* efficiency, TH1* measured, TH1* initialConditions, TH1* result)
{
  // unfolds and assigns bias as errors with Eq. 19 of Cowan, "a survey of unfolding methods for particle physics"
  // b_i = sum_j dmu_i/dn_j (nu_j - n_j) with nu_j as folded guess, n_j as data
  // dmu_i/dn_j is found numerically by changing the bin content and re-unfolding
  //
  // return code: 0 (success) -1 (error: from Unfold(...) )

  if (Unfold(correlation, efficiency, measured, initialConditions, result) != 0)
    return -1;

  TMatrixD matrix(fgMaxInput, fgMaxParams);
  
  TH1* newResult[4];
  for (Int_t loop=0; loop<4; loop++)
    newResult[loop] = (TH1*) result->Clone(Form("newresult_%d", loop));
    
  // change bin-by-bin and built matrix of effects
  for (Int_t m=0; m<fgMaxInput; m++)
  {
    if (measured->GetBinContent(m+1) < 1)
      continue;
      
    for (Int_t loop=0; loop<4; loop++)
    {
      //Printf("%d %d", i, loop);
      Float_t factor = 1;
      switch (loop)
      {
        case 0: factor = 0.5; break;
        case 1: factor = -0.5; break;
        case 2: factor = 1; break;
        case 3: factor = -1; break;
        default: return -1;
      }
      
      TH1* measuredClone = (TH1*) measured->Clone("measuredClone");
  
      measuredClone->SetBinContent(m+1, measured->GetBinContent(m+1) + factor * TMath::Sqrt(measured->GetBinContent(m+1)));
      //new TCanvas; measuredClone->Draw("COLZ");
    
      newResult[loop]->Reset();
      Unfold(correlation, efficiency, measuredClone, measuredClone, newResult[loop]);
      // WARNING if we leave here, then nothing is calculated
        //return -1;
      
      delete measuredClone;
    }
    
    for (Int_t t=0; t<fgMaxParams; t++)
    {
      // one value
      //matrix(i, j) = (result->GetBinContent(j+1) - newResult->GetBinContent(j+1)) / TMath::Sqrt(changedHist->GetBinContent(1, i+1));
      
      // four values from the derivate (procedure taken from ROOT -- suggestion by Ruben)
      // = 1/2D * [ 8 (f(D/2) - f(-D/2)) - (f(D)-f(-D)) ]/3
      
      /*
      // check formula
      measured->SetBinContent(1, m+1, 100);
      newResult[0]->SetBinContent(t+1, 5 * 0.5 + 10);
      newResult[1]->SetBinContent(t+1, 5 * -0.5 + 10);
      newResult[2]->SetBinContent(t+1, 5 * 1 + 10);
      newResult[3]->SetBinContent(t+1, 5 * -1 + 10);
      */
      
      matrix(m, t) = 0.5 / TMath::Sqrt(measured->GetBinContent(m+1)) * 
        ( 8. * (newResult[0]->GetBinContent(t+1) - newResult[1]->GetBinContent(t+1)) - 
              (newResult[2]->GetBinContent(t+1) - newResult[3]->GetBinContent(t+1))
        ) / 3;
    }
  }
  
  for (Int_t loop=0; loop<4; loop++)
    delete newResult[loop];
  
  // ...calculate folded guess  
  TH1* convoluted = (TH1*) measured->Clone("convoluted");
  convoluted->Reset();
  for (Int_t m=0; m<fgMaxInput; m++)
  {
    Float_t value = 0;
    for (Int_t t = 0; t<fgMaxParams; t++)
    {
      Float_t tmp = correlation->GetBinContent(t+1, m+1) * result->GetBinContent(t+1);
      if (efficiency)
        tmp *= efficiency->GetBinContent(t+1);
      value += tmp;
    }
    convoluted->SetBinContent(m+1, value);
  }
  
  // rescale
  convoluted->Scale(measured->Integral() / convoluted->Integral());
  
  //new TCanvas; convoluted->Draw(); measured->Draw("SAME"); measured->SetLineColor(2);
  //return;
  
  // difference
  convoluted->Add(measured, -1);
  
  // ...and the bias
  for (Int_t t = 0; t<fgMaxParams; t++)
  {
    Double_t error = 0;
    for (Int_t m=0; m<fgMaxInput; m++)
      error += matrix(m, t) * convoluted->GetBinContent(m+1); 
    result->SetBinError(t+1, error);
  }
  
  //new TCanvas; result->Draw(); gPad->SetLogy();
  
  return 0;
}
