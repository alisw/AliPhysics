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

TMatrixD* AliUnfolding::fgCorrelationMatrix = 0;
TMatrixD* AliUnfolding::fgCorrelationCovarianceMatrix = 0;
TVectorD* AliUnfolding::fgCurrentESDVector = 0;
TVectorD* AliUnfolding::fgEntropyAPriori = 0;

TF1* AliUnfolding::fgFitFunction = 0;

AliUnfolding::MethodType AliUnfolding::fgMethodType = AliUnfolding::kInvalid;
Int_t AliUnfolding::fgMaxInput  = -1;  // bins in measured histogram
Int_t AliUnfolding::fgMaxParams = -1;  // bins in unfolded histogram = number of fit params
Float_t AliUnfolding::fgOverflowBinLimit = -1;

AliUnfolding::RegularizationType AliUnfolding::fgRegularizationType = AliUnfolding::kPol1;
Float_t AliUnfolding::fgRegularizationWeight = 10000;
Int_t AliUnfolding::fgSkipBinsBegin = 0;
Float_t AliUnfolding::fgMinuitStepSize = 0.1;                 // (usually not needed to be changed) step size in minimization
Bool_t AliUnfolding::fgNormalizeInput = kFALSE;                  // normalize input spectrum

Float_t AliUnfolding::fgBayesianSmoothing  = 1;           // smoothing parameter (0 = no smoothing)
Int_t   AliUnfolding::fgBayesianIterations = 10;         // number of iterations in Bayesian method

Bool_t AliUnfolding::fgDebug = kFALSE;

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
void AliUnfolding::SetStaticVariables(TH2* correlation, TH1* measured)
{
  // fill static variables needed for minuit fit

  if (!fgCorrelationMatrix)
    fgCorrelationMatrix = new TMatrixD(fgMaxInput, fgMaxParams);
  if (!fgCorrelationCovarianceMatrix)
    fgCorrelationCovarianceMatrix = new TMatrixD(fgMaxInput, fgMaxInput);
  if (!fgCurrentESDVector)
    fgCurrentESDVector = new TVectorD(fgMaxInput);
  if (!fgEntropyAPriori)
    fgEntropyAPriori = new TVectorD(fgMaxParams);

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
      correlation->SetBinContent(i, j, correlation->GetBinContent(i, j) / sum);
      correlation->SetBinError(i, j, correlation->GetBinError(i, j) / sum);

      if (i <= fgMaxParams && j <= fgMaxInput)
        (*fgCorrelationMatrix)(j-1, i-1) = correlation->GetBinContent(i, j);
    }

    //printf("MPV for Ntrue = %f is %f\n", fCurrentCorrelation->GetXaxis()->GetBinCenter(i), fCurrentCorrelation->GetYaxis()->GetBinCenter(maxBin));
  }
    
  //normalize measured
  if (fgNormalizeInput)
    measured->Scale(1.0 / measured->Integral());
  
  for (Int_t i=0; i<fgMaxInput; ++i)
  {
    (*fgCurrentESDVector)[i] = measured->GetBinContent(i+1);
    if (measured->GetBinError(i+1) > 0)
      (*fgCorrelationCovarianceMatrix)(i, i) = (Double_t) 1e-6 / measured->GetBinError(i+1) / measured->GetBinError(i+1);

    if ((*fgCorrelationCovarianceMatrix)(i, i) > 1e7)
      (*fgCorrelationCovarianceMatrix)(i, i) = 0;
    //Printf("%d, %e", i, (*fgCorrelationCovarianceMatrix)(i, i));
  }
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
  // if <check> is true no unfolding is made, instead only the chi2 without unfolding is printed
  //
  // returns minuit status (0 = success), (-1 when check was set)
  //

  SetStaticVariables(correlation, measured);
  
  // Initialize TMinuit via generic fitter interface
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, fgMaxParams);
  Double_t arglist[100];

  // disable any output (-1), unfortuantly we do not see warnings anymore then. Have to find another way...
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);

  // however, enable warnings
  //minuit->ExecuteCommand("SET WAR", arglist, 0);

  // set minimization function
  minuit->SetFCN(Chi2Function);

  for (Int_t i=0; i<fgMaxParams; i++)
    (*fgEntropyAPriori)[i] = 1;

  if (!initialConditions) {
    initialConditions = measured;
  } else {
    Printf("AliUnfolding::UnfoldWithMinuit: Using different initial conditions...");
    //new TCanvas; initialConditions->DrawCopy();
    if (fgNormalizeInput)
      initialConditions->Scale(1.0 / initialConditions->Integral());
  }

  // set initial conditions as a-priori distribution for MRX regularization
  for (Int_t i=0; i<fgMaxParams; i++)
    if (initialConditions->GetBinContent(i+1) > 0)
      (*fgEntropyAPriori)[i] = initialConditions->GetBinContent(i+1);

  Double_t* results = new Double_t[fgMaxParams+1];
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    results[i] = initialConditions->GetBinContent(i+1);

    // minuit sees squared values to prevent it from going negative...
    results[i] = TMath::Sqrt(results[i]);

    minuit->SetParameter(i, Form("param%d", i), results[i], fgMinuitStepSize, 0, 0);
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
  //minuit->ExecuteCommand("SCAN", arglist, 0);
  Int_t status = minuit->ExecuteCommand("MIGRAD", arglist, 1);
  Printf("AliUnfolding::UnfoldWithMinuit: MINUIT status is %d", status);
  //printf("!!!!!!!!!!!!!! MIGRAD finished: Starting MINOS !!!!!!!!!!!!!!");
  //minuit->ExecuteCommand("MINOS", arglist, 0);

  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    results[i] = minuit->GetParameter(i);
    Double_t value = results[i] * results[i];
    // error is : (relError) * (value) = (minuit->GetParError(i) / minuit->GetParameter(i)) * (minuit->GetParameter(i) * minuit->GetParameter(i))
    Double_t error = minuit->GetParError(i) * results[i];
    
    if (efficiency)
    {
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
    
    result->SetBinContent(i+1, value);
    result->SetBinError(i+1, error);
  }
  
  if (fgDebug)
    DrawGuess(results);

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
    if (efficiency)
    {
      efficiency[t] = aEfficiency->GetBinContent(t+1);
    }
    else
      efficiency[t] = 1;
      
    prior[t] = measuredCopy[t];
    result[t] = 0;
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
  
  // unfold...
  for (Int_t i=0; i<fgBayesianIterations || fgBayesianIterations< 0; i++)
  {
    if (fgDebug)
      Printf("AliUnfolding::UnfoldWithBayesian: iteration %i", i);

    // calculate IR from Beyes theorem
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
    
    // draw intermediate result
    /*
    for (Int_t t=0; t<kMaxT; t++)
      aResult->SetBinContent(t+1, result[t]);
    aResult->SetMarkerStyle(20+i);
    aResult->SetMarkerColor(2);
    aResult->DrawCopy("P SAME HIST");
    */

    Double_t chi2LastIter = 0;
    // regularization (simple smoothing)
    for (Int_t t=kStartBin; t<kMaxT; t++)
    {
      Float_t newValue = 0;
      
      // 0 bin excluded from smoothing
      if (t > kStartBin+2 && t<kMaxT-1)
      {
        Float_t average = (result[t-1] + result[t] + result[t+1]) / 3;

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

  for (Int_t t=0; t<kMaxT; t++)
    aResult->SetBinContent(t+1, result[t]);

  delete[] measuredCopy;
  delete[] measuredError;
  delete[] prior;
  delete[] result;
  delete[] efficiency;

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

  Double_t chi2 = 0;

  for (Int_t i=1+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    Double_t right  = params[i];
    Double_t left   = params[i-1];

    if (right != 0)
    {
      Double_t diff = 1 - left / right;
      chi2 += diff * diff;
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

  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0)
      continue;

    Double_t right  = params[i];
    Double_t middle = params[i-1];
    Double_t left   = params[i-2];

    Double_t der1 = (right - middle);
    Double_t der2 = (middle - left);

    Double_t diff = (der1 - der2) / middle;

    chi2 += diff * diff;
  }

  return chi2;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationLog(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers linear function (pol1)

  Double_t chi2 = 0;

  for (Int_t i=2+fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0)
      continue;

    Double_t right  = log(params[i]);
    Double_t middle = log(params[i-1]);
    Double_t left   = log(params[i-2]);

    Double_t der1 = (right - middle);
    Double_t der2 = (middle - left);

    Double_t diff = (der1 - der2) / middle;

    chi2 += diff * diff;
  }

  return chi2 * 100;
}

//____________________________________________________________________
Double_t AliUnfolding::RegularizationTotalCurvature(TVectorD& params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // minimizes the total curvature (from Unfolding Methods In High-Energy Physics Experiments,
  // V. Blobel (Hamburg U.) . DESY 84/118, Dec 1984. 40pp.

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

  Double_t paramSum = 0;
  
  for (Int_t i=fgSkipBinsBegin; i<fgMaxParams; ++i)
    paramSum += params[i];

  Double_t chi2 = 0;
  for (Int_t i=fgSkipBinsBegin; i<fgMaxParams; ++i)
  {
    //Double_t tmp = params[i] / paramSum;
    Double_t tmp = params[i];
    if (tmp > 0 && (*fgEntropyAPriori)[i] > 0)
    {
      chi2 += tmp * TMath::Log(tmp / (*fgEntropyAPriori)[i]);
    }
  }

  return 100.0 + chi2;
}

//____________________________________________________________________
void AliUnfolding::Chi2Function(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t)
{
  //
  // fit function for minuit
  // does: (m - Ad)W(m - Ad) where m = measured, A correlation matrix, d = guess, W = covariance matrix
  //

  // d
  static TVectorD paramsVector(fgMaxParams);
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
  }

  //if (penaltyVal > 1e10)
  //  paramsVector2.Print();

  penaltyVal *= fgRegularizationWeight;

  // Ad
  TVectorD measGuessVector(fgMaxInput);
  measGuessVector = (*fgCorrelationMatrix) * paramsVector;

  // Ad - m
  measGuessVector -= (*fgCurrentESDVector);

  //measGuessVector.Print();

  TVectorD copy(measGuessVector);

  // (Ad - m) W
  // this step can be optimized because currently only the diagonal elements of fgCorrelationCovarianceMatrix are used
  // normal way is like this:
  // measGuessVector *= (*fgCorrelationCovarianceMatrix);
  // optimized way like this:
  for (Int_t i=0; i<fgMaxInput; ++i)
    measGuessVector[i] *= (*fgCorrelationCovarianceMatrix)(i, i);

  //measGuessVector.Print();

  // (Ad - m) W (Ad - m)
  // the factor 1e6 prevents that a very small number (measGuessVector[i]) is multiplied with a very
  // big number ((*fgCorrelationCovarianceMatrix)(i, i)) (see ApplyMinuitFit)
  Double_t chi2FromFit = measGuessVector * copy * 1e6;

  chi2 = chi2FromFit + penaltyVal;

  static Int_t callCount = 0;
  if ((callCount++ % 10000) == 0)
    Printf("AliUnfolding::Chi2Function: Iteration %d: %f %f --> %f", callCount, chi2FromFit, penaltyVal, chi2);
}

//____________________________________________________________________
void AliUnfolding::DrawGuess(Double_t *params)
{
  //
  // draws residuals of solution suggested by params and effect of regularization
  //

  // d
  static TVectorD paramsVector(fgMaxParams);
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramsVector[i] = params[i] * params[i];

  // Ad
  TVectorD measGuessVector(fgMaxInput);
  measGuessVector = (*fgCorrelationMatrix) * paramsVector;

  // Ad - m
  measGuessVector -= (*fgCurrentESDVector);

  TVectorD copy(measGuessVector);
  //copy.Print();

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
  TH1* residuals = new TH1F("residuals", "residuals", fgMaxInput, -0.5, fgMaxInput - 0.5);
  for (Int_t i=0; i<fgMaxInput; ++i)
    residuals->SetBinContent(i+1, measGuessVector(i) * copy(i) * 1e6);
  new TCanvas; residuals->DrawCopy(); gPad->SetLogy();

  // draw penalty
  if (fgRegularizationType != kPol1) {
    Printf("Drawing not supported");
    return;
  }

  TH1* penalty = new TH1F("penalty", "penalty", fgMaxParams, -0.5, fgMaxParams - 0.5);

  for (Int_t i=2+1; i<fgMaxParams; ++i)
  {
    if (params[i-1] == 0)
      continue;

    Double_t right  = params[i];
    Double_t middle = params[i-1];
    Double_t left   = params[i-2];

    Double_t der1 = (right - middle);
    Double_t der2 = (middle - left);

    Double_t diff = (der1 - der2) / middle;

    penalty->SetBinContent(i-1, diff * diff);
    
    //Printf("%d %f %f %f %f", i-1, left, middle, right, diff);
  }
  new TCanvas; penalty->DrawCopy(); gPad->SetLogy();
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
  
  SetStaticVariables(correlation, measured);
  
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
