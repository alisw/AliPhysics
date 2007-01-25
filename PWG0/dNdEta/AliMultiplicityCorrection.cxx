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

ClassImp(AliMultiplicityCorrection)

const Int_t AliMultiplicityCorrection::fgMaxParams = 90;
TH1* AliMultiplicityCorrection::fCurrentESD = 0;
TH1* AliMultiplicityCorrection::fCurrentCorrelation = 0;

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection() :
  TNamed()
{
  //
  // default constructor
  //

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i] = 0;

  for (Int_t i = 0; i < kMCHists; ++i)
    fMultiplicityMC[i] = 0;

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = 0;
    fMultiplicityESDCorrected[i] = 0;
  }
}

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection(const Char_t* name, const Char_t* title) :
  TNamed(name, title)
{
  //
  // named constructor
  //

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Float_t binLimitsVtx[] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
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

  #define BINNING fgMaxParams, binLimitsN

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i] = new TH2F(Form("fMultiplicityESD%d", i), "fMultiplicityESD;vtx-z;Ntracks;Count", 10, binLimitsVtx, BINNING);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityMC[i] = dynamic_cast<TH2F*> (fMultiplicityESD[0]->Clone(Form("fMultiplicityMC%d", i)));
    fMultiplicityMC[i]->SetTitle("fMultiplicityMC");
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = new TH3F(Form("fCorrelation%d", i), "fCorrelation;vtx-z;Npart;Ntracks", 10, binLimitsVtx, BINNING, BINNING);
    fMultiplicityESDCorrected[i] = new TH1F(Form("fMultiplicityESDCorrected%d", i), "fMultiplicityESDCorrected;Npart;dN/dN", BINNING);
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
  TList collections[kESDHists+kMCHists+kCorrHists*2];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliMultiplicityCorrection* entry = dynamic_cast<AliMultiplicityCorrection*> (obj);
    if (entry == 0) 
      continue;

    for (Int_t i = 0; i < kESDHists; ++i)
      collections[i].Add(entry->fMultiplicityESD[i]);

    for (Int_t i = 0; i < kMCHists; ++i)
      collections[kESDHists+i].Add(entry->fMultiplicityMC[i]);

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists+i].Add(entry->fCorrelation[i]);

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists+kCorrHists+i].Add(entry->fMultiplicityESDCorrected[i]);

    count++;
  }

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i]->Merge(&collections[i]);

  for (Int_t i = 0; i < kMCHists; ++i)
    fMultiplicityMC[i]->Merge(&collections[kESDHists+i]);

  for (Int_t i = 0; i < kCorrHists; ++i)
    fCorrelation[i]->Merge(&collections[kESDHists+kMCHists+i]);

  for (Int_t i = 0; i < kCorrHists; ++i)
    fMultiplicityESDCorrected[i]->Merge(&collections[kESDHists+kMCHists+kCorrHists+i]);

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
    fMultiplicityMC[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityMC[i]->GetName()));
    if (!fMultiplicityMC[i])
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
    if (fMultiplicityMC[i])
      fMultiplicityMC[i]->Write();

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
void AliMultiplicityCorrection::FillGenerated(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated15, Int_t generated20, Int_t generatedAll)
{
  //
  // Fills an event from MC
  //

  fMultiplicityMC[0]->Fill(vtx, generated05);
  fMultiplicityMC[1]->Fill(vtx, generated10);
  fMultiplicityMC[2]->Fill(vtx, generated15);
  fMultiplicityMC[3]->Fill(vtx, generated20);
  fMultiplicityMC[4]->Fill(vtx, generatedAll);
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

    Double_t right  = params[i] / fCurrentESD->GetBinWidth(i+1);
    Double_t left   = params[i-1] / fCurrentESD->GetBinWidth(i);

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

    Double_t right  = params[i] / fCurrentESD->GetBinWidth(i+1);
    Double_t middle = params[i-1] / fCurrentESD->GetBinWidth(i);
    Double_t left   = params[i-2] / fCurrentESD->GetBinWidth(i-1);

    Double_t der1 = (right - middle) / fCurrentESD->GetBinWidth(i);
    Double_t der2 = (middle - left)  / fCurrentESD->GetBinWidth(i-1);

    Double_t diff = der1 - der2;

    // TODO give additonal weight to big bins
    chi2 += diff * diff * fCurrentESD->GetBinWidth(i) * fCurrentESD->GetBinWidth(i-1);

    //printf("%d --> %f\n", i, diff);
  }

  return chi2 / 1e5 / 2;
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

    Double_t right  = params[i] / fCurrentESD->GetBinWidth(i+1);
    Double_t middle = params[i-1] / fCurrentESD->GetBinWidth(i);
    Double_t left   = params[i-2] / fCurrentESD->GetBinWidth(i-1);

    Double_t der1 = (right - middle) / fCurrentESD->GetBinWidth(i);
    Double_t der2 = (middle - left)  / fCurrentESD->GetBinWidth(i-1);

    Double_t secDer = (der1 - der2) / fCurrentESD->GetBinWidth(i);

    // square and weight with the bin width
    chi2 += secDer * secDer * fCurrentESD->GetBinWidth(i) * fCurrentESD->GetBinWidth(i-1);

    //printf("%d --> %f\n", i, secDer);
  }

  return chi2 / 1e5;
}

//____________________________________________________________________
void AliMultiplicityCorrection::MinuitFitFunction(Int_t&, Double_t*, Double_t& chi2, Double_t *params, Int_t)
{
  //
  // fit function for minuit
  //

  // TODO take errors into account

  static Int_t callCount = 0;

  Double_t chi2FromFit = 0;

  // loop over multiplicity (x axis of fMultiplicityESD)
  for (Int_t i=1; i<=fCurrentESD->GetNbinsX(); ++i)
  {
    if (i > fCurrentCorrelation->GetNbinsY())
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
      //printf("%f  ", sum);
    }

    Double_t diff = fCurrentESD->GetBinContent(i) - sum;
    if (fCurrentESD->GetBinContent(i) > 0)
      diff /= fCurrentESD->GetBinContent(i);
    else
      diff /= fCurrentESD->Integral();

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

  Double_t penaltyVal = RegularizationTotalCurvature(params);

  chi2 = chi2FromFit * chi2FromFit + penaltyVal * penaltyVal;

  if ((callCount++ % 100) == 0)
    printf("%f %f --> %f\n", chi2FromFit, penaltyVal, chi2);
}

//____________________________________________________________________
void AliMultiplicityCorrection::SetupCurrentHists(Int_t inputRange, Bool_t fullPhaseSpace)
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

  SetupCurrentHists(inputRange, fullPhaseSpace);

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
    }

    printf("MPV for Ntrue = %f is %f\n", fCurrentCorrelation->GetXaxis()->GetBinCenter(i), fCurrentCorrelation->GetYaxis()->GetBinCenter(maxBin));
  }

  // small hack to get around charge conservation for full phase space ;-)
  /*if (fullPhaseSpace)
  {
    for (Int_t i=2; i<=50; i+=2)
      for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
        fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i-1, j));
  }*/

  TCanvas* canvas1 = new TCanvas("ApplyMinuitFit", "ApplyMinuitFit", 800, 400);
  canvas1->Divide(2, 1);
  canvas1->cd(1);
  fCurrentESD->DrawCopy();
  canvas1->cd(2);
  fCurrentCorrelation->DrawCopy("COLZ");

  // Initialize TMinuit via generic fitter interface
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, fgMaxParams);

  minuit->SetFCN(MinuitFitFunction);

  Double_t results[fgMaxParams];

  TH1* proj = fMultiplicityMC[mcTarget]->ProjectionY();
  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    //results[i] = 1;
    results[i] = fCurrentESD->GetBinContent(i+1);
    if (check)
      results[i] = proj->GetBinContent(i+1);
    minuit->SetParameter(i, Form("param%d", i), results[i], ((results[i] > 1) ? (results[i] / 10) : 0), 0, fCurrentESD->GetMaximum() * 100);
  }

  Int_t dummy;
  Double_t chi2;
  MinuitFitFunction(dummy, 0, chi2, results, 0);
  printf("Chi2 of right solution is = %f\n", chi2);

  if (check)
    return;

  Double_t arglist[100];
  arglist[0] = 100000;
  //minuit->ExecuteCommand("SET PRINT", arglist, 1);
  minuit->ExecuteCommand("MIGRAD", arglist, 1);
  //minuit->ExecuteCommand("MINOS", arglist, 0);

  for (Int_t i=0; i<fgMaxParams; ++i)
  {
    results[i] = minuit->GetParameter(i);
    fMultiplicityESDCorrected[correlationID]->SetBinContent(i+1, results[i]);
    fMultiplicityESDCorrected[correlationID]->SetBinError(i+1, 0);
  }

  printf("Penalty is %f\n", RegularizationTotalCurvature(results));

  DrawComparison("MinuitChi2", mcTarget, correlationID);
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

  printf("MC:\n");

  TCanvas* canvas2 = new TCanvas("fMultiplicityMC", "fMultiplicityMC", 900, 600);
  canvas2->Divide(3, 2);
  for (Int_t i = 0; i < kMCHists; ++i)
  {
    canvas2->cd(i+1);
    fMultiplicityMC[i]->DrawCopy("COLZ");
    printf("%d --> %f\n", i, fMultiplicityMC[i]->ProjectionY()->GetMean());
  }

  TCanvas* canvas3 = new TCanvas("fCorrelation", "fCorrelation", 900, 900);
  canvas3->Divide(3, 3);
  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    canvas3->cd(i+1);
    TH1* proj = fCorrelation[i]->Project3D("zy");
    proj->DrawCopy("COLZ");
  }
}

//____________________________________________________________________
void AliMultiplicityCorrection::DrawComparison(const char* name, Int_t mcID, Int_t esdCorrId, Bool_t normalizeESD)
{
  TString tmpStr;
  tmpStr.Form("%s_DrawComparison_%d_%d", name, mcID, esdCorrId);

  TCanvas* canvas1 = new TCanvas(tmpStr, tmpStr, 900, 300);
  canvas1->Divide(3, 1);

  canvas1->cd(1);
  TH1* proj = fMultiplicityMC[mcID]->ProjectionY();
  NormalizeToBinWidth(proj);

  if (normalizeESD)
    NormalizeToBinWidth(fMultiplicityESDCorrected[esdCorrId]);

  proj->GetXaxis()->SetRangeUser(0, 200);
  proj->DrawCopy("HIST");
  gPad->SetLogy();

  NormalizeToBinWidth(fCurrentESD);
  fCurrentESD->SetLineColor(2);
  fCurrentESD->DrawCopy("HISTSAME");

  fMultiplicityESDCorrected[esdCorrId]->SetMarkerStyle(3);
  fMultiplicityESDCorrected[esdCorrId]->DrawCopy("SAMEP");

  canvas1->cd(2);
  fMultiplicityESDCorrected[esdCorrId]->GetXaxis()->SetRangeUser(0, 200);
  fMultiplicityESDCorrected[esdCorrId]->DrawCopy("HIST");

  canvas1->cd(3);
  TH1* clone = dynamic_cast<TH1*> (proj->Clone("clone"));
  clone->Divide(fMultiplicityESDCorrected[esdCorrId]);
  clone->SetTitle("Ratio;Npart;MC/ESD");
  clone->GetYaxis()->SetRangeUser(0.8, 1.2);
  clone->GetXaxis()->SetRangeUser(0, 200);
  clone->DrawCopy("HIST");

  /*TLegend* legend = new TLegend(0.6, 0.7, 0.85, 0.85);
  legend->AddEntry(fMultiplicityESDCorrected, "ESD corrected");
  legend->AddEntry(fMultiplicityMC, "MC");
  legend->AddEntry(fMultiplicityESD, "ESD");
  legend->Draw();*/
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyBayesianMethod(Int_t inputRange, Bool_t fullPhaseSpace)
{
  //
  // correct spectrum using bayesian method
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace);

  // normalize correction for given nPart
  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());
    if (sum <= 0)
      continue;
    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
    {
      // npart sum to 1
      fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
      fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
    }
  }

  Float_t regPar = 0.01;

  Float_t norm = 0;
  // pick prior distribution
  TH1F* hPrior = (TH1F*)fCurrentESD->Clone("prior");
  for (Int_t t=1; t<=hPrior->GetNbinsX(); t++)
    norm = norm + hPrior->GetBinContent(t);
  for (Int_t t=1; t<=hPrior->GetNbinsX(); t++) {
    hPrior->SetBinContent(t, hPrior->GetBinContent(t)/norm);

    printf(" bin %d content %f \n", t, hPrior->GetBinContent(t));

  }
  
  // define hist to store guess of true
  TH1F* hTrue  = (TH1F*)fCurrentESD->Clone("prior");
  //  hTrue->Reset();

  // clone response R
  TH2F* hResponse = (TH2F*)fCurrentCorrelation->Clone("pij");

  // histogram to store the inverse response calculated from bayes theorem (from prior and response)
  // IR
  //TAxis* xAxis = hResponse->GetYaxis();
  //TAxis* yAxis = hResponse->GetXaxis();

  TH2F* hInverseResponseBayes = (TH2F*)hResponse->Clone("pji");
  //new TH2F("pji","pji",
  //					 xAxis->GetNbins(),xAxis->GetXbins()->GetArray(),
  //					 yAxis->GetNbins(),yAxis->GetXbins()->GetArray());
  hInverseResponseBayes->Reset();
  
  // unfold...
  Int_t iterations = 50;
  for (Int_t i=0; i<iterations; i++) {
    printf(" iteration %i \n", i);
    
    // calculate IR from Beyes theorem
    // IR_ji = R_ij * prior_i / sum_k(R_kj * prior_k)
    for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++) {
      for (Int_t m=1; m<=hResponse->GetNbinsY(); m++) {

        norm = 0;
        for(Int_t t2 = 1; t2<=hResponse->GetNbinsX(); t2++)
          norm += (hResponse->GetBinContent(t2,m) * hPrior->GetBinContent(t2));

        if (norm==0)
	  hInverseResponseBayes->SetBinContent(t,m,0);
        else
	  hInverseResponseBayes->SetBinContent(t,m, hResponse->GetBinContent(t,m) * hPrior->GetBinContent(t)/norm);
      }
    }
    // calculate true assuming IR is the correct inverse response
    for (Int_t t = 1; t<=hResponse->GetNbinsX(); t++) {
      hTrue ->SetBinContent(t, 0);
      for (Int_t m=1; m<=hResponse->GetNbinsY(); m++)
        hTrue->SetBinContent(t, hTrue->GetBinContent(t) + fCurrentESD->GetBinContent(m)*hInverseResponseBayes->GetBinContent(t,m));
    }

    // regularization (simple smoothing) NB : not good bin width!!!
    TH1F* hTrueSmoothed = (TH1F*)hTrue->Clone("truesmoothed");

    for (Int_t t=2; t<hTrue->GetNbinsX()-1; t++) {
      Float_t average = (hTrue->GetBinContent(t-1) / hTrue->GetBinWidth(t-1)
                         + hTrue->GetBinContent(t) / hTrue->GetBinWidth(t)
                         + hTrue->GetBinContent(t+1) / hTrue->GetBinWidth(t+1)) / 3.;
      average *= hTrueSmoothed->GetBinWidth(t);

      // weight the average with the regularization parameter
      hTrueSmoothed->SetBinContent(t, (1-regPar)*hTrue->GetBinContent(t) + (regPar*average));
    }

    // use smoothed true (normalized) as new prior
    norm = 0;
    for (Int_t t=1; t<=hPrior->GetNbinsX(); t++)
      norm = norm + hTrueSmoothed->GetBinContent(t);
    for (Int_t t=1; t<=hPrior->GetNbinsX(); t++) {
      hPrior->SetBinContent(t, hTrueSmoothed->GetBinContent(t)/norm);
      hTrue->SetBinContent(t, hTrueSmoothed->GetBinContent(t));
    }

    delete hTrueSmoothed;
  } // end of iterations

  for (Int_t t=1; t<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); t++) {
    fMultiplicityESDCorrected[correlationID]->SetBinContent(t, hTrue->GetBinContent(t));
    fMultiplicityESDCorrected[correlationID]->SetBinError(t, 0.05*hTrue->GetBinContent(t));

    printf(" bin %d content %f \n", t, fMultiplicityESDCorrected[correlationID]->GetBinContent(t));
  }
  
  DrawComparison("Bayesian", mcTarget, correlationID);
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyGaussianMethod(Int_t inputRange, Bool_t fullPhaseSpace)
{
  //
  // correct spectrum using a simple Gaussian approach, that is model-dependent
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace);

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

    Int_t fillBegin = fineBinned->FindBin(mean - width * 3);
    Int_t fillEnd   = fineBinned->FindBin(mean + width * 3);
    printf("bin %d mean %f width %f, filling from %d to %d\n", i, mean, width, fillBegin, fillEnd);

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

  DrawComparison("Gaussian", mcTarget, correlationID, kFALSE);
}
