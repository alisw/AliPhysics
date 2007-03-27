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

/* $Id$ */

// This class is used to store correction maps, raw input and results of the multiplicity
// measurement with the ITS or TPC
// It also contains functions to correct the spectrum using different methods.
//
//  Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

#include "AliMultiplicityCorrection.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TDirectory.h>
#include <TVirtualFitter.h>
#include <TCanvas.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>
#include <TCollection.h>

ClassImp(AliMultiplicityCorrection)

const Int_t AliMultiplicityCorrection::fgMaxParams = 200;
TH1* AliMultiplicityCorrection::fCurrentESD = 0;
TH1* AliMultiplicityCorrection::fCurrentCorrelation = 0;
TH1* AliMultiplicityCorrection::fCurrentEfficiency = 0;
TMatrixF* AliMultiplicityCorrection::fCorrelationMatrix = 0;
TMatrixF* AliMultiplicityCorrection::fCorrelationCovarianceMatrix = 0;
TVectorF* AliMultiplicityCorrection::fCurrentESDVector = 0;
AliMultiplicityCorrection::RegularizationType AliMultiplicityCorrection::fRegularizationType = AliMultiplicityCorrection::kPol1;
Float_t AliMultiplicityCorrection::fRegularizationWeight = 1e4;

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection() :
  TNamed(), fLastChi2MC(0), fLastChi2Residuals(0)
{
  //
  // default constructor
  //

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i] = 0;

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i] = 0;
    fMultiplicityMB[i] = 0;
    fMultiplicityINEL[i] = 0;
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = 0;
    fMultiplicityESDCorrected[i] = 0;
  }
}

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection(const Char_t* name, const Char_t* title) :
  TNamed(name, title), fLastChi2MC(0), fLastChi2Residuals(0)
{
  //
  // named constructor
  //

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  /*Float_t binLimitsVtx[] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
  Float_t binLimitsN[] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5,
                          10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5,
                          20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5,
                          30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5,
                          40.5, 41.5, 42.5, 43.5, 44.5, 45.5, 46.5, 47.5, 48.5, 49.5,
                          50.5, 55.5, 60.5, 65.5, 70.5, 75.5, 80.5, 85.5, 90.5, 95.5,
                          100.5, 105.5, 110.5, 115.5, 120.5, 125.5, 130.5, 135.5, 140.5, 145.5,
                          150.5, 160.5, 170.5, 180.5, 190.5, 200.5, 210.5, 220.5, 230.5, 240.5,
                          250.5, 275.5, 300.5, 325.5, 350.5, 375.5, 400.5, 425.5, 450.5, 475.5,
                          500.5 };
                                //525.5, 550.5, 575.5, 600.5, 625.5, 650.5, 675.5, 700.5, 725.5,
                          //750.5, 775.5, 800.5, 825.5, 850.5, 875.5, 900.5, 925.5, 950.5, 975.5,
                          //1000.5 };

  #define VTXBINNING 10, binLimitsVtx
  #define NBINNING fgMaxParams, binLimitsN*/

  #define NBINNING 251, -0.5, 250.5
  #define VTXBINNING 10, -10, 10

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i] = new TH2F(Form("fMultiplicityESD%d", i), "fMultiplicityESD;vtx-z;Ntracks;Count", VTXBINNING, NBINNING);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i] = dynamic_cast<TH2F*> (fMultiplicityESD[0]->Clone(Form("fMultiplicityVtx%d", i)));
    fMultiplicityVtx[i]->SetTitle("fMultiplicityVtx;vtx-z;Npart");

    fMultiplicityMB[i] = dynamic_cast<TH2F*> (fMultiplicityVtx[0]->Clone(Form("fMultiplicityMB%d", i)));
    fMultiplicityMB[i]->SetTitle("fMultiplicityMB");

    fMultiplicityINEL[i] = dynamic_cast<TH2F*> (fMultiplicityVtx[0]->Clone(Form("fMultiplicityINEL%d", i)));
    fMultiplicityINEL[i]->SetTitle("fMultiplicityINEL");
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = new TH3F(Form("fCorrelation%d", i), "fCorrelation;vtx-z;Npart;Ntracks", VTXBINNING, NBINNING, NBINNING);
    fMultiplicityESDCorrected[i] = new TH1F(Form("fMultiplicityESDCorrected%d", i), "fMultiplicityESDCorrected;Npart;dN/dN", NBINNING);
  }

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
AliMultiplicityCorrection::~AliMultiplicityCorrection()
{
  //
  // Destructor
  //
}

//____________________________________________________________________
Long64_t AliMultiplicityCorrection::Merge(TCollection* list)
{
  // Merge a list of AliMultiplicityCorrection objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections[kESDHists+kMCHists*3+kCorrHists*2];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliMultiplicityCorrection* entry = dynamic_cast<AliMultiplicityCorrection*> (obj);
    if (entry == 0) 
      continue;

    for (Int_t i = 0; i < kESDHists; ++i)
      collections[i].Add(entry->fMultiplicityESD[i]);

    for (Int_t i = 0; i < kMCHists; ++i)
    {
      collections[kESDHists+i].Add(entry->fMultiplicityVtx[i]);
      collections[kESDHists+kMCHists+i].Add(entry->fMultiplicityMB[i]);
      collections[kESDHists+kMCHists*2+i].Add(entry->fMultiplicityINEL[i]);
    }

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists*3+i].Add(entry->fCorrelation[i]);

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists*3+kCorrHists+i].Add(entry->fMultiplicityESDCorrected[i]);

    count++;
  }

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i]->Merge(&collections[i]);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i]->Merge(&collections[kESDHists+i]);
    fMultiplicityMB[i]->Merge(&collections[kESDHists+kMCHists+i]);
    fMultiplicityINEL[i]->Merge(&collections[kESDHists+kMCHists*2+i]);
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
    fCorrelation[i]->Merge(&collections[kESDHists+kMCHists*3+i]);

  for (Int_t i = 0; i < kCorrHists; ++i)
    fMultiplicityESDCorrected[i]->Merge(&collections[kESDHists+kMCHists*3+kCorrHists+i]);

  delete iter;

  return count+1;
}

//____________________________________________________________________
Bool_t AliMultiplicityCorrection::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  // TODO memory leak. old histograms needs to be deleted.

  Bool_t success = kTRUE;

  for (Int_t i = 0; i < kESDHists; ++i)
  {
    fMultiplicityESD[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityESD[i]->GetName()));
    if (!fMultiplicityESD[i])
      success = kFALSE;
  }

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityVtx[i]->GetName()));
    fMultiplicityMB[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityMB[i]->GetName()));
    fMultiplicityINEL[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityINEL[i]->GetName()));
    if (!fMultiplicityVtx[i] || !fMultiplicityMB[i] || !fMultiplicityINEL[i])
      success = kFALSE;
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = dynamic_cast<TH3F*> (gDirectory->Get(fCorrelation[i]->GetName()));
    if (!fCorrelation[i])
      success = kFALSE;
    fMultiplicityESDCorrected[i] = dynamic_cast<TH1F*> (gDirectory->Get(fMultiplicityESDCorrected[i]->GetName()));
    if (!fMultiplicityESDCorrected[i])
      success = kFALSE;
  }

  gDirectory->cd("..");

  return success;
}

//____________________________________________________________________
void AliMultiplicityCorrection::SaveHistograms()
{
  //
  // saves the histograms
  //

  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());

  for (Int_t i = 0; i < kESDHists; ++i)
    if (fMultiplicityESD[i])
      fMultiplicityESD[i]->Write();

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    if (fMultiplicityVtx[i])
      fMultiplicityVtx[i]->Write();
    if (fMultiplicityMB[i])
      fMultiplicityMB[i]->Write();
    if (fMultiplicityINEL[i])
      fMultiplicityINEL[i]->Write();
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    if (fCorrelation[i])
      fCorrelation[i]->Write();
    if (fMultiplicityESDCorrected[i])
      fMultiplicityESDCorrected[i]->Write();
  }

  gDirectory->cd("..");
}

//____________________________________________________________________
void AliMultiplicityCorrection::FillGenerated(Float_t vtx, Bool_t triggered, Bool_t vertex, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll)
{
  //
  // Fills an event from MC
  //

  if (triggered)
  {
    fMultiplicityMB[0]->Fill(vtx, generated05);
    fMultiplicityMB[1]->Fill(vtx, generated10);
    fMultiplicityMB[2]->Fill(vtx, generated15);
    fMultiplicityMB[3]->Fill(vtx, generated20);
    fMultiplicityMB[4]->Fill(vtx, generatedAll);

    if (vertex)
    {
      fMultiplicityVtx[0]->Fill(vtx, generated05);
      fMultiplicityVtx[1]->Fill(vtx, generated10);
      fMultiplicityVtx[2]->Fill(vtx, generated15);
      fMultiplicityVtx[3]->Fill(vtx, generated20);
      fMultiplicityVtx[4]->Fill(vtx, generatedAll);
    }
  }

  fMultiplicityINEL[0]->Fill(vtx, generated05);
  fMultiplicityINEL[1]->Fill(vtx, generated10);
  fMultiplicityINEL[2]->Fill(vtx, generated15);
  fMultiplicityINEL[3]->Fill(vtx, generated20);
  fMultiplicityINEL[4]->Fill(vtx, generatedAll);
}

//____________________________________________________________________
void AliMultiplicityCorrection::FillMeasured(Float_t vtx, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20)
{
  //
  // Fills an event from ESD
  //

  fMultiplicityESD[0]->Fill(vtx, measured05);
  fMultiplicityESD[1]->Fill(vtx, measured10);
  fMultiplicityESD[2]->Fill(vtx, measured15);
  fMultiplicityESD[3]->Fill(vtx, measured20);
}

//____________________________________________________________________
void AliMultiplicityCorrection::FillCorrection(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll, Int_t measured05, Int_t measured10, Int_t measured15, Int_t measured20)
{
  //
  // Fills an event into the correlation map with the information from MC and ESD
  //

  fCorrelation[0]->Fill(vtx, generated05, measured05);
  fCorrelation[1]->Fill(vtx, generated10, measured10);
  fCorrelation[2]->Fill(vtx, generated15, measured15);
  fCorrelation[3]->Fill(vtx, generated20, measured20);

  fCorrelation[4]->Fill(vtx, generatedAll, measured05);
  fCorrelation[5]->Fill(vtx, generatedAll, measured10);
  fCorrelation[6]->Fill(vtx, generatedAll, measured15);
  fCorrelation[7]->Fill(vtx, generatedAll, measured20);
}

//____________________________________________________________________
Double_t AliMultiplicityCorrection::RegularizationPol0(Double_t *params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers constant function (pol0)

  Double_t chi2 = 0;

  for (Int_t i=1; i<fgMaxParams; ++i)
  {
    if (params[i] == 0)
      continue;

    Double_t right  = params[i]; // / fCurrentESD->GetBinWidth(i+1);
    Double_t left   = params[i-1]; // / fCurrentESD->GetBinWidth(i);

    Double_t diff = (right - left) / right;
    chi2 += diff * diff;
  }

  return chi2;
}

//____________________________________________________________________
Double_t AliMultiplicityCorrection::RegularizationPol1(Double_t *params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // prefers linear function (pol1)

  Double_t chi2 = 0;

  for (Int_t i=2; i<fgMaxParams; ++i)
  {
    if (params[i] == 0 || params[i-1] == 0)
      continue;

    Double_t right  = params[i]; // / fCurrentESD->GetBinWidth(i+1);
    Double_t middle = params[i-1]; // / fCurrentESD->GetBinWidth(i);
    Double_t left   = params[i-2]; // / fCurrentESD->GetBinWidth(i-1);

    Double_t der1 = (right - middle); // / fCurrentESD->GetBinWidth(i);
    Double_t der2 = (middle - left); //  / fCurrentESD->GetBinWidth(i-1);

    Double_t diff = (der1 - der2) / middle; /// fCurrentESD->GetBinWidth(i);

    // TODO give additonal weight to big bins
    chi2 += diff * diff; // * fCurrentESD->GetBinWidth(i) * fCurrentESD->GetBinWidth(i-1);

    //printf("%d --> %f\n", i, diff);
  }

  return chi2; // 10000
}

//____________________________________________________________________
Double_t AliMultiplicityCorrection::RegularizationTotalCurvature(Double_t *params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // minimizes the total curvature (from Unfolding Methods In High-Energy Physics Experiments,
  // V. Blobel (Hamburg U.) . DESY 84/118, Dec 1984. 40pp.

  Double_t chi2 = 0;

  for (Int_t i=2; i<fgMaxParams; ++i)
  {
    if (params[i] == 0 || params[i-1] == 0)
      continue;

    Double_t right  = params[i]; // / fCurrentESD->GetBinWidth(i+1);
    Double_t middle = params[i-1]; // / fCurrentESD->GetBinWidth(i);
    Double_t left   = params[i-2]; // / fCurrentESD->GetBinWidth(i-1);

    Double_t der1 = (right - middle); // / fCurrentESD->GetBinWidth(i);
    Double_t der2 = (middle - left); //  / fCurrentESD->GetBinWidth(i-1);

    Double_t secDer = (der1 - der2); // / fCurrentESD->GetBinWidth(i);

    // square and weight with the bin width
    chi2 += secDer * secDer; // * fCurrentESD->GetBinWidth(i) * fCurrentESD->GetBinWidth(i-1);

    //printf("%d --> %f\n", i, secDer);
  }

  return chi2; // 10
}

//____________________________________________________________________
Double_t AliMultiplicityCorrection::RegularizationEntropy(Double_t *params)
{
  // homogenity term for minuit fitting
  // pure function of the parameters
  // calculates entropy, from
  // The method of reduced cross-entropy (M. Schmelling 1993)

  //static Int_t callCount = 0;

  Double_t paramSum = 0;
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramSum += params[i]; // / fCurrentESD->GetBinWidth(i+1);

  Double_t chi2 = 0;
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    Double_t tmp = params[i] / paramSum; // / fCurrentESD->GetBinWidth(i+1);
    if (tmp > 0)
    {
      chi2 += tmp * log(tmp); /* * fCurrentESD->GetBinWidth(i+1);
      //                                   /\
      //                                    ------------------------------------
      // TODO WARNING the absolute fitting values seem to depend on that value |
      // this should be impossible...
      //if ((callCount % 100) == 0)
      //  printf("%f => %f\n", params[i], tmp * log(tmp)); */
    }
  }
  //if ((callCount++ % 100) == 0)
  //  printf("\n");

  return chi2; // 1000
}

//____________________________________________________________________
void AliMultiplicityCorrection::MinuitFitFunction(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t)
{
  //
  // fit function for minuit
  // does: (m - Ad)W(m - Ad) where m = measured, A correlation matrix, d = guess, W = covariance matrix
  //

  static Int_t callCount = 0;

  Double_t chi2FromFit = 0;

  // d
  TVectorF paramsVector(fgMaxParams);
  for (Int_t i=0; i<fgMaxParams; ++i)
    paramsVector[i] = params[i];

  // Ad
  paramsVector = (*fCorrelationMatrix) * paramsVector;

  // Ad - m
  paramsVector -= (*fCurrentESDVector);

  TVectorF copy(paramsVector);

  // (Ad - m) W
  paramsVector *= (*fCorrelationCovarianceMatrix);

  // (Ad - m) W (Ad - m)
  chi2FromFit = paramsVector * copy;

  /*printf("chi2FromFit = %f\n", chi2FromFit);
  chi2FromFit = 0;

  // loop over multiplicity (x axis of fMultiplicityESD)
  for (Int_t i=1; i<=fCurrentESD->GetNbinsX(); ++i)
  {
    if (i > fCurrentCorrelation->GetNbinsY() || i > fgMaxParams)
      break;

    Double_t sum = 0;
    //Double_t error = 0;
    // loop over generated particles (x axis of fCorrelation) that resulted in reconstruction of i tracks
    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsX(); ++j)
    {
      if (j > fgMaxParams)
        break;

      sum += fCurrentCorrelation->GetBinContent(j, i) * params[j-1];

      //if (params[j-1] > 0)
      //  error += fCurrentCorrelation->GetBinError(j, i) * fCurrentCorrelation->GetBinError(j, i) * params[j-1];
    }

    //printf("%f\n", sum);

    Double_t diff = fCurrentESD->GetBinContent(i) - sum;

    // include error
    if (fCurrentESD->GetBinError(i) > 0)
      diff /= fCurrentESD->GetBinError(i);

    //if (fCurrentESD->GetBinContent(i) > 0)
    //  diff /= fCurrentESD->GetBinContent(i);
    //else
    //  diff /= fCurrentESD->Integral();

    // weight with bin width
    //diff *= fCurrentESD->GetBinWidth(i);

    //error = TMath::Sqrt(error) + fCurrentESD->GetBinError(i);
    //if (error <= 0)
    // error = 1;

    //Double_t tmp = diff / error;
    //chi2 += tmp * tmp;
    chi2FromFit += diff * diff;

    //printf("\nExpected sum = %f; Diff for bin %d is %f\n**********************************\n", fCurrentESD->GetBinContent(i), i, diff);
    //printf("Diff for bin %d is %f\n", i, diff);
  }

  printf("chi2FromFit = %f\n", chi2FromFit);*/

  Double_t penaltyVal = 0;

  switch (fRegularizationType)
  {
    case kNone:      break;
    case kPol0:      penaltyVal = RegularizationPol0(params); break;
    case kPol1:      penaltyVal = RegularizationPol1(params); break;
    case kCurvature: penaltyVal = RegularizationTotalCurvature(params); break;
    case kEntropy:   penaltyVal = RegularizationEntropy(params); break;
  }

  penaltyVal *= fRegularizationWeight;

  chi2 = chi2FromFit + penaltyVal;

  if ((callCount++ % 1000) == 0)
    printf("%f %f --> %f\n", chi2FromFit, penaltyVal, chi2);
}

//____________________________________________________________________
void AliMultiplicityCorrection::SetupCurrentHists(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType)
{
  //
  // fills fCurrentESD, fCurrentCorrelation
  // resets fMultiplicityESDCorrected
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  fMultiplicityESDCorrected[correlationID]->Reset();

  fCurrentESD = fMultiplicityESD[inputRange]->ProjectionY();
  fCurrentESD->Sumw2();

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  TH3* hist = fCorrelation[correlationID];
  for (Int_t y=1; y<=hist->GetYaxis()->GetNbins(); ++y)
  {
    for (Int_t z=1; z<=hist->GetZaxis()->GetNbins(); ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  fCurrentCorrelation = hist->Project3D("zy");
  fCurrentCorrelation->Sumw2();

  fCurrentEfficiency = fMultiplicityVtx[inputRange]->ProjectionY("CurrentEfficiency");
  TH2* divisor = 0;
  switch (eventType)
  {
    case kTrVtx : divisor = fMultiplicityVtx[inputRange]; break;
    case kMB: divisor = fMultiplicityMB[inputRange]; break;
    case kINEL: divisor = fMultiplicityINEL[inputRange]; break;
  }
  fCurrentEfficiency->Divide(divisor->ProjectionY());
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyMinuitFit(Int_t inputRange, Bool_t fullPhaseSpace, Bool_t check)
{
  //
  // correct spectrum using minuit chi2 method
  //
  // if check is kTRUE the input MC solution (by definition the right solution) is used
  // no fit is made and just the chi2 is calculated
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace, kTrVtx);

  fCorrelationMatrix = new TMatrixF(fgMaxParams, fgMaxParams);
  fCorrelationCovarianceMatrix = new TMatrixF(fgMaxParams, fgMaxParams);
  fCurrentESDVector = new TVectorF(fgMaxParams);

  // normalize correction for given nPart
  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());
    if (sum <= 0)
      continue;
    Float_t maxValue = 0;
    Int_t maxBin = -1;
    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
    {
      // find most probably value
      if (maxValue < fCurrentCorrelation->GetBinContent(i, j))
      {
        maxValue = fCurrentCorrelation->GetBinContent(i, j);
        maxBin = j;
      }

      // npart sum to 1
      fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
      fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);

      if (i <= fgMaxParams && j <= fgMaxParams)
        (*fCorrelationMatrix)(j-1, i-1) = fCurrentCorrelation->GetBinContent(i, j);
    }

    //printf("MPV for Ntrue = %f is %f\n", fCurrentCorrelation->GetXaxis()->GetBinCenter(i), fCurrentCorrelation->GetYaxis()->GetBinCenter(maxBin));
  }

  // small hack to get around charge conservation for full phase space ;-)
  /*if (fullPhaseSpace)
  {
    for (Int_t i=2; i<=50; i+=2)
      for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
        fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i-1, j));
  }*/

  /*TCanvas* canvas1 = new TCanvas("ApplyMinuitFit", "ApplyMinuitFit", 800, 400);
  canvas1->Divide(2, 1);
  canvas1->cd(1);
  fCurrentESD->DrawCopy();
  canvas1->cd(2);
  fCurrentCorrelation->DrawCopy("COLZ");*/

  // Initialize TMinuit via generic fitter interface
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, fgMaxParams);

  minuit->SetFCN(MinuitFitFunction);

  static Double_t results[fgMaxParams];
  //printf("%x\n", (void*) results);

  TH1* proj = fMultiplicityVtx[mcTarget]->ProjectionY();
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    //results[i] = 1;
    results[i] = fCurrentESD->GetBinContent(i+1);
    if (results[i] < 1e-2)
      results[i] = 1e-2;
    if (check)
      results[i] = proj->GetBinContent(i+1);
    minuit->SetParameter(i, Form("param%d", i), results[i], results[i] / 10, 0.001, fCurrentESD->GetMaximum() * 10);

    (*fCurrentESDVector)[i] = fCurrentESD->GetBinContent(i+1);
    if (fCurrentESD->GetBinError(i+1) > 0)
      (*fCorrelationCovarianceMatrix)(i, i) = 1.0 / (fCurrentESD->GetBinError(i+1) * fCurrentESD->GetBinError(i+1));

    //minuit->SetParameter(i, Form("param%d", i), fCurrentESD->GetBinWidth(i+1), 0.01, 0.001, 50);
  }
  minuit->SetParameter(0, "param0", results[1], ((results[1] > 1) ? (results[1] / 10) : 0), 0.001, fCurrentESD->GetMaximum() * 10);
  //minuit->SetParameter(0, "param0", 1, 0, 1, 1);

  Int_t dummy;
  Double_t chi2;
  MinuitFitFunction(dummy, 0, chi2, results, 0);
  printf("Chi2 of right solution is = %f\n", chi2);

  if (check)
    return;

  Double_t arglist[100];
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  minuit->ExecuteCommand("MIGRAD", arglist, 1);
  //minuit->ExecuteCommand("MINOS", arglist, 0);

  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    results[i] = minuit->GetParameter(i);
    fMultiplicityESDCorrected[correlationID]->SetBinContent(i+1, results[i]);
    fMultiplicityESDCorrected[correlationID]->SetBinError(i+1, 0);
  }

  //printf("Penalty is %f\n", RegularizationTotalCurvature(results));
}

//____________________________________________________________________
void AliMultiplicityCorrection::NormalizeToBinWidth(TH1* hist)
{
  //
  // normalizes a 1-d histogram to its bin width
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
  {
    hist->SetBinContent(i, hist->GetBinContent(i) / hist->GetBinWidth(i));
    hist->SetBinError(i, hist->GetBinError(i) / hist->GetBinWidth(i));
  }
}

//____________________________________________________________________
void AliMultiplicityCorrection::NormalizeToBinWidth(TH2* hist)
{
  //
  // normalizes a 2-d histogram to its bin width (x width * y width)
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); ++i)
    for (Int_t j=1; j<=hist->GetNbinsY(); ++j)
    {
      Double_t factor = hist->GetXaxis()->GetBinWidth(i) * hist->GetYaxis()->GetBinWidth(j);
      hist->SetBinContent(i, j, hist->GetBinContent(i, j) / factor);
      hist->SetBinError(i, j, hist->GetBinError(i, j) / factor);
    }
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyMinuitFitAll()
{
  //
  // fit all eta ranges
  //

  for (Int_t i=0; i<kESDHists; ++i)
  {
    ApplyMinuitFit(i, kFALSE);
    ApplyMinuitFit(i, kTRUE);
  }
}

//____________________________________________________________________
void AliMultiplicityCorrection::DrawHistograms()
{
  printf("ESD:\n");

  TCanvas* canvas1 = new TCanvas("fMultiplicityESD", "fMultiplicityESD", 900, 600);
  canvas1->Divide(3, 2);
  for (Int_t i = 0; i < kESDHists; ++i)
  {
    canvas1->cd(i+1);
    fMultiplicityESD[i]->DrawCopy("COLZ");
    printf("%d --> %f\n", i, (Float_t) fMultiplicityESD[i]->ProjectionY()->GetMean());
  }

  printf("Vtx:\n");

  TCanvas* canvas2 = new TCanvas("fMultiplicityMC", "fMultiplicityMC", 900, 600);
  canvas2->Divide(3, 2);
  for (Int_t i = 0; i < kMCHists; ++i)
  {
    canvas2->cd(i+1);
    fMultiplicityVtx[i]->DrawCopy("COLZ");
    printf("%d --> %f\n", i, fMultiplicityVtx[i]->ProjectionY()->GetMean());
  }

  TCanvas* canvas3 = new TCanvas("fCorrelation", "fCorrelation", 900, 900);
  canvas3->Divide(3, 3);
  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    canvas3->cd(i+1);
    TH3* hist = dynamic_cast<TH3*> (fCorrelation[i]->Clone());
    for (Int_t y=1; y<=hist->GetYaxis()->GetNbins(); ++y)
    {
      for (Int_t z=1; z<=hist->GetZaxis()->GetNbins(); ++z)
      {
        hist->SetBinContent(0, y, z, 0);
        hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
      }
    }
    TH1* proj = hist->Project3D("zy");
    proj->DrawCopy("COLZ");
  }
}

//____________________________________________________________________
void AliMultiplicityCorrection::DrawComparison(const char* name, Int_t inputRange, Bool_t fullPhaseSpace, Bool_t normalizeESD, TH1* mcHist)
{
  Int_t esdCorrId = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);

  TString tmpStr;
  tmpStr.Form("%s_DrawComparison_%d", name, esdCorrId);

  // calculate residual

  // for that we convolute the response matrix with the gathered result
  // if normalizeESD is not set, the histogram is already normalized, this needs to be passed to CalculateMultiplicityESD
  TH2* convoluted = CalculateMultiplicityESD(fMultiplicityESDCorrected[esdCorrId], esdCorrId, !normalizeESD);
  TH1* proj2 = convoluted->ProjectionY("proj2", -1, -1, "e");
  NormalizeToBinWidth(proj2);
  TH1* residual = convoluted->ProjectionY("residual", -1, -1, "e");
  residual->SetTitle("Residuals");

  residual->Add(fCurrentESD, -1);
  residual->Divide(residual, fCurrentESD, 1, 1, "B");

  // TODO fix errors
  for (Int_t i=1; i<residual->GetNbinsX(); ++i)
  {
    proj2->SetBinError(i, 0);
    residual->SetBinError(i, 0);
  }

  TCanvas* canvas1 = new TCanvas(tmpStr, tmpStr, 1200, 800);
  canvas1->Divide(2, 2);

  canvas1->cd(1);
  TH1* proj = (TH1*) mcHist->Clone("proj");
  NormalizeToBinWidth(proj);

  if (normalizeESD)
    NormalizeToBinWidth(fMultiplicityESDCorrected[esdCorrId]);

  proj->GetXaxis()->SetRangeUser(0, 150);
  proj->DrawCopy("HIST");
  gPad->SetLogy();

  //fMultiplicityESDCorrected[esdCorrId]->SetMarkerStyle(3);
  fMultiplicityESDCorrected[esdCorrId]->SetLineColor(2);
  fMultiplicityESDCorrected[esdCorrId]->DrawCopy("SAME HIST");

  Float_t chi2 = 0;
  for (Int_t i=1; i<=100; ++i)
  {
    if (fMultiplicityESDCorrected[esdCorrId]->GetBinContent(i) != 0)
    {
      Float_t value = (proj->GetBinContent(i) - fMultiplicityESDCorrected[esdCorrId]->GetBinContent(i)) / fMultiplicityESDCorrected[esdCorrId]->GetBinContent(i);
      chi2 += value * value;
    }
  }

  printf("Difference (from MC) is %f for bin 1-100\n", chi2);
  fLastChi2MC = chi2;

  canvas1->cd(2);

  NormalizeToBinWidth(fCurrentESD);
  gPad->SetLogy();
  fCurrentESD->GetXaxis()->SetRangeUser(0, 150);
  //fCurrentESD->SetLineColor(2);
  fCurrentESD->DrawCopy("HIST");

  proj2->SetLineColor(2);
  //proj2->SetMarkerColor(2);
  //proj2->SetMarkerStyle(5);
  proj2->DrawCopy("HIST SAME");

  chi2 = 0;
  for (Int_t i=1; i<=100; ++i)
  {
    if (fCurrentESD->GetBinContent(i) != 0)
    {
      Float_t value = (proj2->GetBinContent(i) - fCurrentESD->GetBinContent(i)) / fCurrentESD->GetBinContent(i);
      chi2 += value * value;
    }
  }

  printf("Difference (Residuals) is %f for bin 1-100\n", chi2);
  fLastChi2Residuals = chi2;

  canvas1->cd(3);
  TH1* clone = dynamic_cast<TH1*> (proj->Clone("clone"));
  clone->Divide(fMultiplicityESDCorrected[esdCorrId]);
  clone->SetTitle("Ratio;Npart;MC/ESD");
  clone->GetYaxis()->SetRangeUser(0.8, 1.2);
  clone->GetXaxis()->SetRangeUser(0, 150);
  clone->DrawCopy("HIST");

  /*TLegend* legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(fMultiplicityESDCorrected, "ESD corrected");
  legend->AddEntry(fMultiplicityMC, "MC");
  legend->AddEntry(fMultiplicityESD, "ESD");
  legend->Draw();*/

  canvas1->cd(4);

  residual->GetYaxis()->SetRangeUser(-0.2, 0.2);
  residual->GetXaxis()->SetRangeUser(0, 150);
  residual->DrawCopy();

  canvas1->SaveAs(Form("%s.gif", canvas1->GetName()));
}

//____________________________________________________________________
void AliMultiplicityCorrection::GetComparisonResults(Float_t& mc, Float_t& residuals)
{
  // Returns the chi2 between the MC and the unfolded ESD as well as between the ESD and the folded unfolded ESD
  // These values are computed during DrawComparison, thus this function picks up the
  // last calculation

  mc = fLastChi2MC;
  residuals = fLastChi2Residuals;
}


//____________________________________________________________________
TH2F* AliMultiplicityCorrection::GetMultiplicityMC(Int_t i, EventType eventType)
{
  //
  // returns the corresponding MC spectrum
  //

  switch (eventType)
  {
    case kTrVtx : return fMultiplicityVtx[i]; break;
    case kMB : return fMultiplicityMB[i]; break;
    case kINEL : return fMultiplicityINEL[i]; break;
  }

  return 0;
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyBayesianMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Float_t regPar)
{
  //
  // correct spectrum using bayesian method
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);

  // TODO should be taken from correlation map
  TH1* sumHist = GetMultiplicityMC(inputRange, eventType)->ProjectionY("sumHist", 1, GetMultiplicityMC(inputRange, eventType)->GetNbinsX());

  //new TCanvas;
  //sumHist->Draw();

  // normalize correction for given nPart
  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    //Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());
    Double_t sum = sumHist->GetBinContent(i);
    if (sum <= 0)
      continue;
    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
    {
      // npart sum to 1
      fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
      fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
    }
  }

  //new TCanvas;
  //fCurrentCorrelation->Draw("COLZ");
  //return;

  // normalize correction for given nTrack
  /*for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
  {
    Double_t sum = fCurrentCorrelation->Integral(1, fCurrentCorrelation->GetNbinsX(), j, j);
    if (sum <= 0)
      continue;
    for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsY(); ++i)
    {
      // ntrack sum to 1
      fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
      fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
    }
  }*/

  // FAKE
  fCurrentEfficiency = ((TH2*) fCurrentCorrelation)->ProjectionX("eff");
  //new TCanvas;
  //fCurrentEfficiency->Draw();

  // pick prior distribution
  TH1F* hPrior = (TH1F*)fCurrentESD->Clone("prior");
  Float_t norm = hPrior->Integral();
  for (Int_t t=1; t<=hPrior->GetNbinsX(); t++)
    hPrior->SetBinContent(t, hPrior->GetBinContent(t)/norm);

  // define temp hist
  TH1F* hTemp = (TH1F*)fCurrentESD->Clone("temp");
  hTemp->Reset();

  // just a shortcut
  TH2F* hResponse = (TH2F*) fCurrentCorrelation;

  // histogram to store the inverse response calculated from bayes theorem (from prior and response) IR
  TH2F* hInverseResponseBayes = (TH2F*)hResponse->Clone("pji");
  hInverseResponseBayes->Reset();

  // unfold...
  Int_t iterations = 20;
  for (Int_t i=0; i<iterations; i++)
  {
    //printf(" iteration %i \n", i);

    // calculate IR from Beyes theorem
    // IR_ji = R_ij * prior_i / sum_k(R_kj * prior_k)

    for (Int_t m=1; m<=hResponse->GetNbinsY(); m++)
    {
      Float_t norm = 0;
      for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++)
        norm += hResponse->GetBinContent(t,m) * hPrior->GetBinContent(t);
      for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++)
      {
        if (norm==0)
          hInverseResponseBayes->SetBinContent(t,m,0);
        else
          hInverseResponseBayes->SetBinContent(t,m, hResponse->GetBinContent(t,m) * hPrior->GetBinContent(t)/norm);
      }
    }

    for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++)
    {
      Float_t value = 0;
      for (Int_t m=1; m<=hResponse->GetNbinsY(); m++)
        value += fCurrentESD->GetBinContent(m) * hInverseResponseBayes->GetBinContent(t,m);

      /*if (eventType == kTrVtx)
      {
        hTemp->SetBinContent(t, value);
      }
      else*/
      {
        if (fCurrentEfficiency->GetBinContent(t) > 0)
          hTemp->SetBinContent(t, value / fCurrentEfficiency->GetBinContent(t));
        else
          hTemp->SetBinContent(t, 0);
      }
    }

    // this is the last guess, fill before (!) smoothing
    for (Int_t t=1; t<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); t++)
    {
      fMultiplicityESDCorrected[correlationID]->SetBinContent(t, hTemp->GetBinContent(t));
      fMultiplicityESDCorrected[correlationID]->SetBinError(t, 0.05 * hTemp->GetBinContent(t)); // TODO

      //printf(" bin %d content %f \n", t, fMultiplicityESDCorrected[correlationID]->GetBinContent(t));
    }

    // regularization (simple smoothing)
    TH1F* hTrueSmoothed = (TH1F*) hTemp->Clone("truesmoothed");

    for (Int_t t=2; t<hTrueSmoothed->GetNbinsX(); t++)
    {
      Float_t average = (hTemp->GetBinContent(t-1) / hTemp->GetBinWidth(t-1)
                         + hTemp->GetBinContent(t) / hTemp->GetBinWidth(t)
                         + hTemp->GetBinContent(t+1) / hTemp->GetBinWidth(t+1)) / 3.;
      average *= hTrueSmoothed->GetBinWidth(t);

      // weight the average with the regularization parameter
      hTrueSmoothed->SetBinContent(t, (1-regPar) * hTemp->GetBinContent(t) + regPar * average);
    }

    // calculate chi2 (change from last iteration)
    Double_t chi2 = 0;

    // use smoothed true (normalized) as new prior
    Float_t norm = 1;
    //for (Int_t t=1; t<=hPrior->GetNbinsX(); t++)
    //  norm = norm + hTrueSmoothed->GetBinContent(t);

    for (Int_t t=1; t<hTrueSmoothed->GetNbinsX(); t++)
    {
      Float_t newValue = hTrueSmoothed->GetBinContent(t)/norm;
      Float_t diff = hPrior->GetBinContent(t) - newValue;
      chi2 += (Double_t) diff * diff;

      hPrior->SetBinContent(t, newValue);
    }

    //printf("Chi2 of %d iteration = %.10f\n", i, chi2);

    //if (i % 5 == 0)
    //  DrawComparison(Form("Bayesian_%d", i), mcTarget, correlationID, kTRUE, eventType);

    delete hTrueSmoothed;
  } // end of iterations


  return;

  // ********
  // Calculate the covariance matrix, all arguments are taken from NIM,A362,487-498,1995

  printf("Calculating covariance matrix. This may take some time...\n");

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
  cov->Draw("COLZ");
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyLaszloMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType)
{
  //
  // correct spectrum using bayesian method
  //

  Float_t regPar = 0.05;

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);

  // TODO should be taken from correlation map
  //TH1* sumHist = GetMultiplicityMC(inputRange, eventType)->ProjectionY("sumHist", 1, GetMultiplicityMC(inputRange, eventType)->GetNbinsX());

  // normalize correction for given nPart
  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());
    //Double_t sum = sumHist->GetBinContent(i);
    if (sum <= 0)
      continue;
    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
    {
      // npart sum to 1
      fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
      fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
    }
  }

  new TCanvas;
  fCurrentCorrelation->Draw("COLZ");

  // FAKE
  fCurrentEfficiency = ((TH2*) fCurrentCorrelation)->ProjectionX("eff");

  // pick prior distribution
  TH1F* hPrior = (TH1F*)fCurrentESD->Clone("prior");
  Float_t norm = 1; //hPrior->Integral();
  for (Int_t t=1; t<=hPrior->GetNbinsX(); t++)
    hPrior->SetBinContent(t, hPrior->GetBinContent(t)/norm);

  // zero distribution
  TH1F* zero =  (TH1F*)hPrior->Clone("zero");

  // define temp hist
  TH1F* hTemp = (TH1F*)fCurrentESD->Clone("temp");
  hTemp->Reset();

  // just a shortcut
  TH2F* hResponse = (TH2F*) fCurrentCorrelation;

  // unfold...
  Int_t iterations = 20;
  for (Int_t i=0; i<iterations; i++)
  {
    //printf(" iteration %i \n", i);

    for (Int_t m=1; m<=hResponse->GetNbinsY(); m++)
    {
      Float_t value = 0;
      for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++)
        value += hResponse->GetBinContent(t, m) * hPrior->GetBinContent(t);
      hTemp->SetBinContent(m, value);
      //printf("%d %f %f %f\n", m, zero->GetBinContent(m), hPrior->GetBinContent(m), value);
    }

    // regularization (simple smoothing)
    TH1F* hTrueSmoothed = (TH1F*) hTemp->Clone("truesmoothed");

    for (Int_t t=2; t<hTrueSmoothed->GetNbinsX(); t++)
    {
      Float_t average = (hTemp->GetBinContent(t-1) / hTemp->GetBinWidth(t-1)
                         + hTemp->GetBinContent(t) / hTemp->GetBinWidth(t)
                         + hTemp->GetBinContent(t+1) / hTemp->GetBinWidth(t+1)) / 3.;
      average *= hTrueSmoothed->GetBinWidth(t);

      // weight the average with the regularization parameter
      hTrueSmoothed->SetBinContent(t, (1-regPar) * hTemp->GetBinContent(t) + regPar * average);
    }

    for (Int_t m=1; m<=hResponse->GetNbinsY(); m++)
      hTemp->SetBinContent(m, zero->GetBinContent(m) + hPrior->GetBinContent(m) - hTrueSmoothed->GetBinContent(m));

    // fill guess
    for (Int_t t=1; t<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); t++)
    {
      fMultiplicityESDCorrected[correlationID]->SetBinContent(t, hTemp->GetBinContent(t));
      fMultiplicityESDCorrected[correlationID]->SetBinError(t, 0.05 * hTemp->GetBinContent(t)); // TODO

      //printf(" bin %d content %f \n", t, fMultiplicityESDCorrected[correlationID]->GetBinContent(t));
    }


    // calculate chi2 (change from last iteration)
    Double_t chi2 = 0;

    // use smoothed true (normalized) as new prior
    Float_t norm = 1; //hTrueSmoothed->Integral();

    for (Int_t t=1; t<hTrueSmoothed->GetNbinsX(); t++)
    {
      Float_t newValue = hTrueSmoothed->GetBinContent(t)/norm;
      Float_t diff = hPrior->GetBinContent(t) - newValue;
      chi2 += (Double_t) diff * diff;

      hPrior->SetBinContent(t, newValue);
    }

    printf("Chi2 of %d iteration = %.10f\n", i, chi2);

    if (i % 5 == 0)
      DrawComparison(Form("Laszlo_%d", i), inputRange, fullPhaseSpace, kTRUE, GetMultiplicityMC(mcTarget, eventType)->ProjectionY());

    delete hTrueSmoothed;
  } // end of iterations

  DrawComparison("Laszlo", inputRange, fullPhaseSpace, kTRUE, GetMultiplicityMC(mcTarget, eventType)->ProjectionY());
}

//____________________________________________________________________
Float_t AliMultiplicityCorrection::BayesCovarianceDerivate(Float_t matrixM[251][251], TH2* hResponse,
  TH1* fCurrentEfficiency, Int_t k, Int_t i, Int_t r, Int_t u)
{
  //
  //
  //

  Float_t result = 0;

  if (k == u && r == i)
    result += 1.0 / hResponse->GetBinContent(u+1, r+1);

  if (k == u)
    result -= 1.0 / fCurrentEfficiency->GetBinContent(u+1);

  if (r == i)
    result -= matrixM[u][i] * fCurrentEfficiency->GetBinContent(u+1) / hResponse->GetBinContent(u+1, i+1);

  result *= matrixM[k][i];

  return result;
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyGaussianMethod(Int_t inputRange, Bool_t fullPhaseSpace)
{
  //
  // correct spectrum using a simple Gaussian approach, that is model-dependent
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace, kTrVtx);

  NormalizeToBinWidth((TH2*) fCurrentCorrelation);

  TH1D* correction = dynamic_cast<TH1D*> (fCurrentESD->Clone("GaussianMean"));
  correction->SetTitle("GaussianMean");

  TH1D* correctionWidth = dynamic_cast<TH1D*> (fCurrentESD->Clone("GaussianWidth"));
  correction->SetTitle("GaussianWidth");

  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    TH1D* proj = (dynamic_cast<TH2*> (fCurrentCorrelation))->ProjectionX("_px", i, i+1);
    proj->Fit("gaus", "0Q");
    correction->SetBinContent(i, proj->GetFunction("gaus")->GetParameter(1));
    correctionWidth->SetBinContent(i, proj->GetFunction("gaus")->GetParameter(2));

    continue;

    // draw for debugging
    new TCanvas;
    proj->DrawCopy();
    proj->GetFunction("gaus")->DrawCopy("SAME");
  }

  TH1* target = fMultiplicityESDCorrected[correlationID];

  Int_t nBins = target->GetNbinsX()*10+1;
  Float_t* binning = new Float_t[nBins];
  for (Int_t i=1; i<=target->GetNbinsX(); ++i)
    for (Int_t j=0; j<10; ++j)
      binning[(i-1)*10 + j] = target->GetXaxis()->GetBinLowEdge(i) + target->GetXaxis()->GetBinWidth(i) / 10 * j;

  binning[nBins-1] = target->GetXaxis()->GetBinUpEdge(target->GetNbinsX());

  TH1F* fineBinned = new TH1F("targetFineBinned", "targetFineBinned", nBins-1, binning);

  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    Float_t mean = correction->GetBinContent(i);
    Float_t width = correctionWidth->GetBinContent(i);

    Int_t fillBegin = fineBinned->FindBin(mean - width * 5);
    Int_t fillEnd   = fineBinned->FindBin(mean + width * 5);
    //printf("bin %d mean %f width %f, filling from %d to %d\n", i, mean, width, fillBegin, fillEnd);

    for (Int_t j=fillBegin; j <= fillEnd; ++j)
    {
      fineBinned->AddBinContent(j, TMath::Gaus(fineBinned->GetXaxis()->GetBinCenter(j), mean, width, kTRUE) * fCurrentESD->GetBinContent(i));
    }
  }

  for (Int_t i=1; i<=target->GetNbinsX(); ++i)
  {
    Float_t sum = 0;
    for (Int_t j=1; j<=10; ++j)
      sum += fineBinned->GetBinContent((i-1)*10 + j);
    target->SetBinContent(i, sum / 10);
  }

  delete[] binning;

  DrawComparison("Gaussian", inputRange, fullPhaseSpace, kFALSE, GetMultiplicityMC(mcTarget, kTrVtx)->ProjectionY());
}

//____________________________________________________________________
TH2F* AliMultiplicityCorrection::CalculateMultiplicityESD(TH1* inputMC, Int_t correlationMap, Bool_t normalized)
{
  // runs the distribution given in inputMC through the response matrix identified by
  // correlationMap and produces a measured distribution
  // although it is a TH2F the vertex axis is not used at the moment and all entries are filled in mid-vertex
  // if normalized is set, inputMC is expected to be normalized to the bin width

  if (!inputMC)
    return 0;

  if (correlationMap < 0 || correlationMap >= kCorrHists)
    return 0;

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  TH3* hist = fCorrelation[correlationMap];
  for (Int_t y=1; y<=hist->GetYaxis()->GetNbins(); ++y)
  {
    for (Int_t z=1; z<=hist->GetZaxis()->GetNbins(); ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  TH1* corr = hist->Project3D("zy");
  corr->Sumw2();

  // normalize correction for given nPart
  for (Int_t i=1; i<=corr->GetNbinsX(); ++i)
  {
    Double_t sum = corr->Integral(i, i, 1, corr->GetNbinsY());
    if (sum <= 0)
      continue;

    for (Int_t j=1; j<=corr->GetNbinsY(); ++j)
    {
      // npart sum to 1
      corr->SetBinContent(i, j, corr->GetBinContent(i, j) / sum);
      corr->SetBinError(i, j, corr->GetBinError(i, j) / sum);
    }
  }

  TH2F* target = dynamic_cast<TH2F*> (fMultiplicityESD[0]->Clone(Form("%s_measured", inputMC->GetName())));
  target->Reset();

  for (Int_t meas=1; meas<=corr->GetNbinsY(); ++meas)
  {
    Float_t measured = 0;
    Float_t error = 0;

    for (Int_t gen=1; gen<=target->GetNbinsY(); ++gen)
    {
      Int_t mcGenBin = inputMC->GetXaxis()->FindBin(target->GetYaxis()->GetBinCenter(gen));

      Float_t factor = 1;
      if (normalized)
        factor = inputMC->GetBinWidth(mcGenBin);

      measured += inputMC->GetBinContent(mcGenBin) * factor * corr->GetBinContent(gen, meas);
      error += inputMC->GetBinError(mcGenBin) * factor * corr->GetBinContent(gen, meas);
    }

    // bin width of the <measured> axis has to be taken into account
    //measured /= target->GetYaxis()->GetBinWidth(meas);
    //error /= target->GetYaxis()->GetBinWidth(meas);

    //printf("%f +- %f ; %f +- %f \n", inputMC->GetBinContent(meas), inputMC->GetBinError(meas), measured, error);

    target->SetBinContent(target->GetNbinsX() / 2, meas, measured);
    target->SetBinError(target->GetNbinsX() / 2, meas, error);
  }

  return target;
}

//____________________________________________________________________
void AliMultiplicityCorrection::SetGenMeasFromFunc(TF1* inputMC, Int_t id)
{
  // uses the given function to fill the input MC histogram and generates from that
  // the measured histogram by applying the response matrix
  // this can be used to evaluate if the methods work indepedently of the input
  // distribution
  // WARNING does not respect the vertex distribution, just fills central vertex bin

  if (!inputMC)
    return;

  if (id < 0 || id >= kESDHists)
    return;

  TH2F* mc = fMultiplicityINEL[id];

  mc->Reset();

  for (Int_t i=1; i<=mc->GetNbinsY(); ++i)
  {
    mc->SetBinContent(mc->GetNbinsX() / 2, i, inputMC->Eval(mc->GetYaxis()->GetBinCenter(i)) * mc->GetYaxis()->GetBinWidth(i));
    mc->SetBinError(mc->GetNbinsX() / 2, i, 0);
  }

  //new TCanvas;
  //mc->Draw("COLZ");

  TH1* proj = mc->ProjectionY();

  fMultiplicityESD[id] = CalculateMultiplicityESD(proj, id);
}
