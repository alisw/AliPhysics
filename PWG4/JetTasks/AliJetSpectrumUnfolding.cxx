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
//
// This class is used to store correlation maps, generated and reconstructed data of the jet spectrum
// It also contains functions to correct the spectrum using the bayesian unfolding
//

#include "AliJetSpectrumUnfolding.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TDirectory.h>
#include <TVirtualFitter.h>
#include <TCanvas.h>
#include <TString.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TCollection.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TVectorD.h>

#include <THnSparse.h>

ClassImp(AliJetSpectrumUnfolding)

const Int_t AliJetSpectrumUnfolding::fgkNBINSE  = 50;
const Int_t AliJetSpectrumUnfolding::fgkNBINSZ  = 50;
const Int_t AliJetSpectrumUnfolding::fgkNEVENTS = 500000;
const Double_t AliJetSpectrumUnfolding::fgkaxisLowerLimitE = 0.;
const Double_t AliJetSpectrumUnfolding::fgkaxisLowerLimitZ = 0.;
const Double_t AliJetSpectrumUnfolding::fgkaxisUpperLimitE = 250.;
const Double_t AliJetSpectrumUnfolding::fgkaxisUpperLimitZ = 1.;

Float_t AliJetSpectrumUnfolding::fgBayesianSmoothing  = 1;           // smoothing parameter (0 = no smoothing)
Int_t   AliJetSpectrumUnfolding::fgBayesianIterations = 100;         // number of iterations in Bayesian method

//____________________________________________________________________

AliJetSpectrumUnfolding::AliJetSpectrumUnfolding() :
  TNamed(), fCurrentRec(0), fCurrentCorrelation(0), fRecSpectrum(0),  fGenSpectrum(0),
  fUnfSpectrum(0), fCorrelation(0), fLastChi2MC(0), fLastChi2MCLimit(0), fLastChi2Residuals(0), fRatioAverage(0)
{
  //
  // default constructor
  //

  fGenSpectrum = 0;
  fRecSpectrum = 0;
  fUnfSpectrum = 0;
  fCorrelation = 0;
}

//____________________________________________________________________
AliJetSpectrumUnfolding::AliJetSpectrumUnfolding(const Char_t* name, const Char_t* title) :
  TNamed(name, title), fCurrentRec(0), fCurrentCorrelation(0), fRecSpectrum(0),
  fGenSpectrum(0), fUnfSpectrum(0), fCorrelation(0), fLastChi2MC(0), fLastChi2MCLimit(0), fLastChi2Residuals(0), fRatioAverage(0)
{
  //
  // named constructor
  //

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fRecSpectrum = new TH2F("fRecSpectrum", "Reconstructed Spectrum;E^{jet}_{rec} [GeV];z^{lp}_{rec}",
			  fgkNBINSE, fgkaxisLowerLimitE, fgkaxisUpperLimitE,
			  fgkNBINSZ, fgkaxisLowerLimitZ, fgkaxisUpperLimitZ);
  fGenSpectrum = new TH2F("fGenSpectrum", "Generated Spectrum;E^{jet}_{gen} [GeV];z^{lp}_{gen}", 
			  fgkNBINSE, fgkaxisLowerLimitE, fgkaxisUpperLimitE,
			  fgkNBINSZ, fgkaxisLowerLimitZ, fgkaxisUpperLimitZ);
  fUnfSpectrum = new TH2F("fUnfSpectrum", "Unfolded Spectrum;E^{jet} [GeV];z^{lp}", 
			  fgkNBINSE, fgkaxisLowerLimitE, fgkaxisUpperLimitE,
			  fgkNBINSZ, fgkaxisLowerLimitZ, fgkaxisUpperLimitZ);

  const Int_t nbin[4]={fgkNBINSE, fgkNBINSE, fgkNBINSZ, fgkNBINSZ};
  //arrays for bin limits
  Double_t lowEdge[4] = {fgkaxisLowerLimitE, fgkaxisLowerLimitE, fgkaxisLowerLimitZ, fgkaxisLowerLimitZ};
  Double_t upEdge[4]  = {fgkaxisUpperLimitE, fgkaxisUpperLimitE, fgkaxisUpperLimitZ, fgkaxisUpperLimitZ};

  fCorrelation = new THnSparseF("fCorrelation", "Correlation Function", 4, nbin, lowEdge, upEdge);

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
AliJetSpectrumUnfolding::~AliJetSpectrumUnfolding()
{
  //
  // Destructor
  //

  if (fGenSpectrum)
    delete fGenSpectrum;
  fGenSpectrum = 0;

  if (fRecSpectrum)
    delete fRecSpectrum;
  fRecSpectrum = 0;

  if (fUnfSpectrum)
    delete fUnfSpectrum;
  fUnfSpectrum = 0;
 
  if (fCorrelation)
    delete fCorrelation;
  fCorrelation = 0;

}

//____________________________________________________________________
Long64_t AliJetSpectrumUnfolding::Merge(TCollection* list)
{
  // Merge a list of AliJetSpectrumUnfolding objects with this (needed for
  // PROOF).
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections[4];

  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliJetSpectrumUnfolding* entry = dynamic_cast<AliJetSpectrumUnfolding*> (obj);
    if (entry == 0)
      continue;

    collections[0].Add(entry->fGenSpectrum);
    collections[1].Add(entry->fRecSpectrum);
    collections[2].Add(entry->fUnfSpectrum);
    collections[3].Add(entry->fCorrelation);

    count++;
  }

  fGenSpectrum->Merge(&collections[0]);
  fRecSpectrum->Merge(&collections[1]);
  fUnfSpectrum->Merge(&collections[2]);
  fCorrelation->Merge(&collections[3]);

  delete iter;

  return count+1;
}

//____________________________________________________________________
Bool_t AliJetSpectrumUnfolding::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  Bool_t success = kTRUE;

  // store old histograms to delete them later
  TList oldHistograms;
  oldHistograms.SetOwner(1);

  if (fGenSpectrum)  oldHistograms.Add(fGenSpectrum);
  if (fRecSpectrum)  oldHistograms.Add(fRecSpectrum);
  if (fUnfSpectrum)  oldHistograms.Add(fUnfSpectrum);
  if (fCorrelation)  oldHistograms.Add(fCorrelation);

  // load new histograms
  fGenSpectrum = dynamic_cast<TH2F*> (gDirectory->Get(fGenSpectrum->GetName()));
  if (!fGenSpectrum)
    success = kFALSE;

  fRecSpectrum = dynamic_cast<TH2F*> (gDirectory->Get(fRecSpectrum->GetName()));
  if (!fRecSpectrum)
    success = kFALSE;

  fUnfSpectrum = dynamic_cast<TH2F*> (gDirectory->Get(fUnfSpectrum->GetName()));
  if (!fUnfSpectrum)
    success = kFALSE;

  fCorrelation = dynamic_cast<THnSparseF*> (gDirectory->Get(fCorrelation->GetName()));
  if (!fCorrelation)
    success = kFALSE;

  gDirectory->cd("..");

  // delete old histograms
  oldHistograms.Delete();

  return success;
}

//____________________________________________________________________
void AliJetSpectrumUnfolding::SaveHistograms()
{
  //
  // saves the histograms
  //

  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());

  if (fGenSpectrum)
    fGenSpectrum->Write();

  if (fRecSpectrum)
    fRecSpectrum->Write();

  if (fUnfSpectrum)
    fUnfSpectrum->Write();

  if (fCorrelation)
    fCorrelation->Write();

  gDirectory->cd("..");
}

//____________________________________________________________________
void AliJetSpectrumUnfolding::SetupCurrentHists(Bool_t createBigBin)
{
  //
  // resets fUnfSpectrum
  //

  fUnfSpectrum->Reset();
  fUnfSpectrum->Sumw2();

  fCurrentRec = (TH2F*)fRecSpectrum->Clone("fCurrentRec");
  fCurrentRec->Sumw2();

  fCurrentCorrelation = (THnSparseF*)fCorrelation->Clone("fCurrentCorrelation");  
  fCurrentCorrelation->Sumw2();

  Printf("Correlation Matrix has %ld filled bins", fCurrentCorrelation->GetNbins());

  if (createBigBin)
  {
    Int_t maxBinE = 0, maxBinZ = 0;
    Float_t maxE = 0, maxZ = 0;
    for (Int_t me=1; me<=fCurrentRec->GetNbinsX(); me++)
      for (Int_t mz=1; mz<=fCurrentRec->GetNbinsY(); mz++)
      {
        if (fCurrentRec->GetBinContent(me,mz) <= 5 && me>fgkNBINSE/2 && mz>fgkNBINSZ/2)
        {
          maxBinE = me;
	  maxBinZ = mz;
	  maxE = fCurrentRec->GetXaxis()->GetBinCenter(me);
	  maxZ = fCurrentRec->GetYaxis()->GetBinCenter(mz);
          break;
        }
      }

    if (maxBinE > 0 || maxBinZ > 0)
    {
      printf("Bin limit in measured spectrum is e = %d and z = %d.\n", maxBinE, maxBinZ);
      fCurrentRec->SetBinContent(maxBinE, maxBinZ, fCurrentRec->Integral(maxBinE, fCurrentRec->GetNbinsX(), maxBinZ, fCurrentRec->GetNbinsY()));
      for (Int_t me=maxBinE+1; me<=fCurrentRec->GetNbinsX(); me++)
        for (Int_t mz=maxBinZ+1; mz<=fCurrentRec->GetNbinsY(); mz++)
        {
          fCurrentRec->SetBinContent(me, mz, 0);
          fCurrentRec->SetBinError(me, mz, 0);
        }
      // the error is set to sqrt(N), better solution possible?, sum of relative errors of all contributions???
      fCurrentRec->SetBinError(maxBinE, maxBinZ, TMath::Sqrt(fCurrentRec->GetBinContent(maxBinE, maxBinZ)));

      printf("This bin has now %f +- %f entries\n", fCurrentRec->GetBinContent(maxBinE, maxBinZ), fCurrentRec->GetBinError(maxBinE, maxBinZ));

    /*  for (Int_t te=1; te<=NBINSE; te++)
      {
        for (Int_t tz=1; tz<=NBINSZ; tz++)
	{
          Int_t binMin[4] = {te, maxBinE, tz, maxBinZ};
	  Int_t binMax[4] = {NBINSE, NBINSE, NBINSZ, NBINSZ};
	  Float_t sum=0;
	  for (Int_t ite=te; ite<=NBINSE; ite++)
	    for (Int_t itz=tz; itz<=NBINSZ; itz++)
	      for (Int_t ime=maxBinE; ime<=NBINSE; ime++)
	        for (Int_t imz=maxBinZ; imz<=NBINSZ; imz++)
	        {
	          Int_t bin[4] = {ite, ime, itz, imz};
	          sum += fCurrentCorrelation->GetBinContent(bin);
	        }
	  fCurrentCorrelation->SetBinContent(binMin, sum);
          fCurrentCorrelation->SetBinError(binMin, TMath::Sqrt(fCurrentCorrelation->GetBinContent(binMin)));
          printf("create big bin1, nbins = %d, te  = %d, tz = %d\n", NBINSE, te, tz);
          for (Int_t me=maxBinE; me<=NBINSE; me++)
          {
            for (Int_t mz=maxBinZ; mz<=NBINSZ; mz++)
	    {
	      Int_t bin[4] = {te, me, tz, mz};
              fCurrentCorrelation->SetBinContent(bin, 0.);
              fCurrentCorrelation->SetBinError(bin, 0.);
              printf("create big bin2\n");
	    }
	  }
        }
      }*/
      
      for(Int_t idx = 0; idx<=fCurrentCorrelation->GetNbins(); idx++)
      {
        Int_t bin[4];
        Float_t binContent = fCurrentCorrelation->GetBinContent(idx,bin);
        Float_t binError   = fCurrentCorrelation->GetBinError(idx);
        Int_t binMin[4] = {bin[0], maxBinE, bin[2], maxBinZ};
        if ( (bin[1]>maxBinE && bin[1]<=fgkNBINSE) && (bin[3]>maxBinZ && bin[3]<=fgkNBINSZ) )
        {
	  fCurrentCorrelation->SetBinContent(binMin, binContent + fCurrentCorrelation->GetBinContent(binMin));
          fCurrentCorrelation->SetBinError(binMin, binError + TMath::Sqrt(fCurrentCorrelation->GetBinContent(binMin)));
          fCurrentCorrelation->SetBinContent(bin, 0.);
          fCurrentCorrelation->SetBinError(bin, 0.);         
        } 
        printf("create big bin1, nbins = %d, te  = %d, tz = %d\n", fgkNBINSE, bin[0], bin[1]);
      }

      printf("Adjusted correlation matrix!\n");
    }
  } // end Create Big Bin

}

//____________________________________________________________________
void AliJetSpectrumUnfolding::SetBayesianParameters(Float_t smoothing, Int_t nIterations)
{
  //
  // sets the parameters for Bayesian unfolding
  //

  fgBayesianSmoothing = smoothing;
  fgBayesianIterations = nIterations;

  printf("AliJetSpectrumUnfolding::SetBayesianParameters --> Paramaters set to %d iterations with smoothing %f\n", fgBayesianIterations, fgBayesianSmoothing);
}

//____________________________________________________________________
void AliJetSpectrumUnfolding::NormalizeToBinWidth(TH2* const hist)
{
  //
  // normalizes a 2-d histogram to its bin width (x width * y width)
  //

  for (Int_t i=1; i<=hist->GetNbinsX(); i++)
    for (Int_t j=1; j<=hist->GetNbinsY(); j++)
    {
      Double_t factor = hist->GetXaxis()->GetBinWidth(i) * hist->GetYaxis()->GetBinWidth(j);
      hist->SetBinContent(i, j, hist->GetBinContent(i, j) / factor);
      hist->SetBinError(i, j, hist->GetBinError(i, j) / factor);
    }
}

//____________________________________________________________________
void AliJetSpectrumUnfolding::DrawHistograms()
{
  //
  // draws the histograms of this class
  //

  gStyle->SetPalette(1);

  TCanvas* canvas1 = new TCanvas("fRecSpectrum", "fRecSpectrum", 900, 600);
  gPad->SetLogz();    
  fRecSpectrum->DrawCopy("COLZ");

  TCanvas* canvas2 = new TCanvas("fGenSpectrum", "fGenSpectrum", 900, 600);
  canvas2->cd();
  gPad->SetLogz();    
  fGenSpectrum->DrawCopy("COLZ");

  TCanvas* canvas3 = new TCanvas("fUnfSpectrum", "fUnfSpectrum", 900, 600);
  canvas3->cd();
  gPad->SetLogz();    
  fUnfSpectrum->DrawCopy("COLZ");

  TCanvas* canvas4 = new TCanvas("fCorrelation", "fCorrelation", 500, 500);
  canvas1->Divide(2);  

  canvas4->cd(1);
  gPad->SetLogz();
  TH2D* h0 = fCorrelation->Projection(1,0);
  h0->SetXTitle("E^{jet}_{gen} [GeV]");
  h0->SetYTitle("E^{jet}_{rec} [GeV]");
  h0->SetTitle("Projection: Jet Energy");    
  h0->DrawCopy("colz");

  canvas1->cd(2);
  gPad->SetLogz();  
  TH2D* h1 = fCorrelation->Projection(3,2);
  h1->SetXTitle("z^{lp}_{gen}");
  h1->SetYTitle("z^{lp}_{rec}");
  h1->SetTitle("Projection: Leading Particle Fragmentation");        
  h1->DrawCopy("colz");

}

//____________________________________________________________________
void AliJetSpectrumUnfolding::DrawComparison(const char* name, TH2* const genHist)
{
  // 
  // Draws the copmparison plot (gen,rec and unfolded distributions
  //

  if (fUnfSpectrum->Integral() == 0)
  {
    printf("ERROR. Unfolded histogram is empty\n");
    return;
  }

  //regain measured distribution used for unfolding, because the bins were modified in SetupCurrentHists
  //in create big bin
  fCurrentRec = (TH2F*)fRecSpectrum->Clone();
  fCurrentRec->Sumw2();
  fCurrentRec->Scale(1.0 / fCurrentRec->Integral());

  // normalize unfolded result to 1
  fUnfSpectrum->Scale(1.0 / fUnfSpectrum->Integral());

  // find bin with <= 100 entries. this is used as limit for the chi2 calculation
  Int_t mcBinLimitE = 0, mcBinLimitZ = 0;
  for (Int_t i=0; i<genHist->GetNbinsX(); ++i)
    for (Int_t j=0; j<genHist->GetNbinsY(); ++j)
    {
      if (genHist->GetBinContent(i,j) > 100)
      {
        mcBinLimitE = i;
	mcBinLimitZ = j;
      }
      else
        break;
    }
  Printf("AliJetSpectrumUnfolding::DrawComparison: Gen bin limit is (x,y) = (%d,%d)", mcBinLimitE,mcBinLimitZ);

  // scale to 1 true spectrum
  genHist->Sumw2();
  genHist->Scale(1.0 / genHist->Integral());

  // calculate residual
  // for that we convolute the response matrix with the gathered result
  TH2* tmpRecRecorrected = (TH2*) fUnfSpectrum->Clone("tmpRecRecorrected");
  TH2* convoluted = CalculateRecSpectrum(tmpRecRecorrected);
  if (convoluted->Integral() > 0)
    convoluted->Scale(1.0 / convoluted->Integral());
  else
    printf("ERROR: convoluted is empty. Something went wrong calculating the convoluted histogram.\n");

  TH2* residual = (TH2*) convoluted->Clone("residual");
  residual->SetTitle("(R#otimesUnfolded - Reconstructed)/Reconstructed;E^{jet} [GeV]; z^{lp}");

  fCurrentRec->Scale(1./fCurrentRec->Integral());
  residual->Add(fCurrentRec, -1);
  //residual->Divide(residual, fCurrentRec, 1, 1, "B");

  // draw canvas
  TCanvas* canvas1 = new TCanvas(name, name, 1000, 1000);
  canvas1->Divide(2, 3);
  
  Int_t style = 1;
  const Int_t nRGBs = 5;
  const Int_t nCont = 500;

  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);

  canvas1->cd(1);
  gStyle->SetPalette(style);
  gPad->SetLogz();
  genHist->SetTitle("Generated Spectrum;E^{jet}_{gen} [GeV];z^{lp}");
  genHist->SetStats(0);
  genHist->DrawCopy("colz");

  canvas1->cd(2);
  gStyle->SetPalette(style);
  gPad->SetLogz();
  fUnfSpectrum->SetStats(0);
  fUnfSpectrum->DrawCopy("colz");

  canvas1->cd(3);
  gStyle->SetPalette(style);
  gPad->SetLogz();
  fCurrentRec->SetTitle(fRecSpectrum->GetTitle());
  fCurrentRec->SetStats(0);
  fCurrentRec->DrawCopy("colz");

  canvas1->cd(4);
  gStyle->SetPalette(style);
  gPad->SetLogy();  
  TH1D* projGenX = genHist->ProjectionX();
  projGenX->SetTitle("Projection: Jet Energy; E^{jet} [GeV]");
  TH1D* projUnfX = fUnfSpectrum->ProjectionX();
  TH1D* projRecX = fCurrentRec->ProjectionX();  
  projGenX->SetStats(0);
  projRecX->SetStats(0);
  projUnfX->SetStats(0);  
  projGenX->SetLineColor(8);
  projRecX->SetLineColor(2);  
  projGenX->DrawCopy();
  projUnfX->DrawCopy("same");
  projRecX->DrawCopy("same");  

  TLegend* legend = new TLegend(0.6, 0.85, 0.98, 0.98);
  legend->AddEntry(projGenX, "Generated Spectrum");
  legend->AddEntry(projUnfX, "Unfolded Spectrum");
  legend->AddEntry(projRecX, "Reconstructed Spectrum");  
  //legend->SetFillColor(0);
  legend->Draw("same");

  canvas1->cd(5);
  gPad->SetLogy();
  gStyle->SetPalette(style);
  TH1D* projGenY = genHist->ProjectionY();
  projGenY->SetTitle("Projection: Leading Particle Fragmentation; z^{lp}");
  TH1D* projUnfY = fUnfSpectrum->ProjectionY();
  TH1D* projRecY = fCurrentRec->ProjectionY();  
  projGenY->SetStats(0);
  projRecY->SetStats(0);
  projUnfY->SetStats(0);  
  projGenY->SetLineColor(8);
  projRecY->SetLineColor(2);  
  projGenY->DrawCopy();
  projUnfY->DrawCopy("same");
  projRecY->DrawCopy("same");  

  TLegend* legend1 = new TLegend(0.6, 0.85, 0.98, 0.98);
  legend1->AddEntry(projGenY, "Generated Spectrum");
  legend1->AddEntry(projUnfY, "Unfolded Spectrum");
  legend1->AddEntry(projRecY, "Recontructed Spectrum");  
  //legend1->SetFillColor(0);
  legend1->Draw("same");
  
  // Draw residuals
  canvas1->cd(6);
  gStyle->SetPalette(style);
  gPad->SetLogz();
  residual->SetStats(0);
  residual->DrawCopy("colz");
  
  canvas1->SaveAs(Form("%s.png", canvas1->GetName()));
}


//____________________________________________________________________
void AliJetSpectrumUnfolding::GetComparisonResults(Float_t* const gen, Int_t*  const genLimit, Float_t* const residuals, Float_t* const ratioAverage) const
{
  // Returns the chi2 between the Generated and the unfolded Reconstructed spectrum as well as between the Reconstructed and the folded unfolded
  // These values are computed during DrawComparison, thus this function picks up the
  // last calculation

  if (gen)
    *gen = fLastChi2MC;
  if (genLimit)
    *genLimit = fLastChi2MCLimit;
  if (residuals)
    *residuals = fLastChi2Residuals;
  if (ratioAverage)
    *ratioAverage = fRatioAverage;
}

//____________________________________________________________________
void AliJetSpectrumUnfolding::ApplyBayesianMethod(Float_t regPar, Int_t nIterations, TH2* const initialConditions, Bool_t determineError)
{
  //
  // correct spectrum using bayesian unfolding
  //

  // initialize seed with current time 
  gRandom->SetSeed(0);

  printf("seting up current arrays and histograms...\n");
  SetupCurrentHists(kFALSE); // kFALSE to not create big bin

  // normalize Correlation Map to convert number of events into probabilities
  /*for (Int_t te=1; te<=NBINSE; te++)
    for (Int_t tz=1; tz<=NBINSZ; tz++)
    {
       Int_t bin[4];
       Float_t sum=0.;
       for (Int_t me = 1; me<=NBINSE; me++)
         for (Int_t mz = 1; mz<=NBINSZ; mz++)
         {
           bin[0] = te; bin[1] = me; 
           bin[2] = tz; bin[3] = mz;           
           sum += fCurrentCorrelation->GetBinContent(bin);
         }
       if (sum > 0.)
         for (Int_t me = 1; me<=NBINSE; me++)
           for (Int_t mz = 1; mz<=NBINSZ; mz++)
           {
             bin[0] = te; bin[1] = me; 
             bin[2] = tz; bin[3] = mz;           
             fCurrentCorrelation->SetBinContent(bin, fCurrentCorrelation->GetBinContent(bin)/sum);
	     fCurrentCorrelation->SetBinError(bin, fCurrentCorrelation->GetBinError(bin)/sum);
	   }
    }*/
  Float_t sum[fgkNBINSE+2][fgkNBINSZ+2];
  memset(sum,0,sizeof(Float_t)*(fgkNBINSE+2)*(fgkNBINSZ+2));

  for (Int_t idx=0; idx<=fCurrentCorrelation->GetNbins(); idx++)
  {
    Int_t bin[4];
    Float_t binContent = fCurrentCorrelation->GetBinContent(idx, bin);
    if ( (bin[1]>0 && bin[1]<=fgkNBINSE) && (bin[3]>0 && bin[3]<=fgkNBINSZ) )
      sum[bin[0]][bin[2]] += binContent; 
  }
  
  for (Int_t idx=0; idx<=fCurrentCorrelation->GetNbins(); idx++)
  {
    Int_t bin[4];
    Float_t binContent = fCurrentCorrelation->GetBinContent(idx, bin);
    Float_t binError   = fCurrentCorrelation->GetBinError(bin);
    if (sum[bin[0]][bin[2]]>0 && (bin[1]>0 && bin[1]<=fgkNBINSE) &&
        (bin[3]>0 && bin[3]<=fgkNBINSZ) && (bin[0]>0 && bin[0]<=fgkNBINSE) && (bin[2]>0 && bin[2]<=fgkNBINSZ) )
    {
      fCurrentCorrelation->SetBinContent(bin, binContent/sum[bin[0]][bin[2]]);
      fCurrentCorrelation->SetBinError(bin, binError/sum[bin[0]][bin[2]]);    
    }  
  }
    
  printf("calling UnfoldWithBayesian\n");
  Int_t success = UnfoldWithBayesian(fCurrentCorrelation, fCurrentRec, initialConditions, fUnfSpectrum, regPar, nIterations, kFALSE); 
  
  if ( success != 0)
    return;

  if (!determineError)
  {
    Printf("AliJetSpectrumUnfolding::ApplyBayesianMethod: WARNING: No errors calculated.");
    return;
  }

  // evaluate errors, this is done by randomizing the measured spectrum following Poission statistics
  // this (new) measured spectrum is then unfolded and the different to the result from the "real" measured
  // spectrum calculated. This is performed N times and the maximum difference is taken as the systematic
  // error of the unfolding method itself.

  TH2* maxError = (TH2*) fUnfSpectrum->Clone("maxError");
  maxError->Reset();

  TH2* normalizedResult = (TH2*) fUnfSpectrum->Clone("normalizedResult");
  normalizedResult->Scale(1.0 / normalizedResult->Integral());

  const Int_t kErrorIterations = 20;

  printf("Spectrum unfolded. Determining error (%d iterations)...\n", kErrorIterations);

  TH2* randomized = (TH2*) fCurrentRec->Clone("randomized");
  TH2* result2 = (TH2*) fUnfSpectrum->Clone("result2");
  for (Int_t n=0; n<kErrorIterations; ++n)
  {
    // randomize the content of clone following a poisson with the mean = the value of that bin
    for (Int_t x=1; x<=randomized->GetNbinsX(); x++)
      for (Int_t y=1; y<=randomized->GetNbinsY(); y++)
      {
        Float_t randomValue = fCurrentRec->GetBinContent(x,y);
        TF1* poisson = new TF1("poisson", "TMath::Poisson(x,[0])",randomValue*0.25, randomValue*1.25);
        poisson->SetParameters(randomValue,0.);
        randomValue = poisson->GetRandom();   
        //printf("%e --> %e\n", fCurrentRec->GetBinContent(x,y), (Double_t)randomValue);
        randomized->SetBinContent(x, y, randomValue);
        delete poisson;
      }

    result2->Reset();
    if (UnfoldWithBayesian(fCurrentCorrelation, randomized, initialConditions, result2, regPar, nIterations) != 0)
      return;

    result2->Scale(1.0 / result2->Integral());

    // calculate ratio
    result2->Divide(normalizedResult);

    //new TCanvas; result2->DrawCopy("HIST");

    // find max. deviation
    for (Int_t i=1; i<=result2->GetNbinsX(); i++)
      for (Int_t j=1; j<=result2->GetNbinsY(); j++)
        maxError->SetBinContent(i, j, TMath::Max(maxError->GetBinContent(i,j), TMath::Abs(1 - result2->GetBinContent(i,j))));
  }
  delete randomized;
  delete result2;

  for (Int_t i=1; i<=fUnfSpectrum->GetNbinsX(); i++)
    for (Int_t j=1; j<=fUnfSpectrum->GetNbinsY(); j++)
      fUnfSpectrum->SetBinError(i, j, fUnfSpectrum->GetBinError(i,j) + maxError->GetBinContent(i,j)*fUnfSpectrum->GetBinContent(i,j));

  delete maxError;
  delete normalizedResult;
}

//____________________________________________________________________
Int_t AliJetSpectrumUnfolding::UnfoldWithBayesian(THnSparseF* const correlation, TH2* const measured, TH2* const initialConditions, TH2* const aResult, Float_t regPar, Int_t nIterations, Bool_t calculateErrors)
{
  //
  // unfolds a spectrum
  //

  if (measured->Integral() <= 0)
  {
    Printf("AliJetSpectrumUnfolding::UnfoldWithBayesian: ERROR: The measured spectrum is empty");
    return 1;
  }
  const Int_t nFillesBins = correlation->GetNbins();  
  const Int_t kStartBin = 1;

  const Int_t kMaxTZ = fgkNBINSZ; // max true axis fragmentation function
  const Int_t kMaxMZ = fgkNBINSZ; // max measured axis fragmentation function
  const Int_t kMaxTE = fgkNBINSE; // max true axis energy
  const Int_t kMaxME = fgkNBINSE; // max measured axis energy
  
  printf("NbinsE=%d - NbinsZ=%d\n", fgkNBINSE, fgkNBINSZ);

  // store information in arrays, to increase processing speed 
  Double_t measuredCopy[kMaxME+1][kMaxMZ+1];
  Double_t prior[kMaxTE+1][kMaxTZ+1];
  Double_t errors[kMaxTE+1][kMaxTZ+1];
  Double_t result[kMaxTE+1][kMaxTZ+1];

  THnSparseF *inverseCorrelation;
  inverseCorrelation = (THnSparseF*)correlation->Clone("inverseCorrelation");
  inverseCorrelation->Reset();
  
  Float_t inputDistIntegral = 1;
  if (initialConditions)
  {
    printf("Using different starting conditions...\n");   
    inputDistIntegral = initialConditions->Integral();
  }
  Float_t measuredIntegral = measured->Integral();  
  for (Int_t me=1; me<=kMaxME; me++)
    for (Int_t mz=1; mz<=kMaxMZ; mz++)
    {
      // normalization of the measured spectrum
      measuredCopy[me][mz] = measured->GetBinContent(me,mz) / measuredIntegral;
      errors[me][mz] = measured->GetBinError(me, mz) / measuredIntegral;
      // pick prior distribution and normalize it
      if (initialConditions)
        prior[me][mz] = initialConditions->GetBinContent(me,mz) / inputDistIntegral;
      else
        prior[me][mz] = measured->GetBinContent(me,mz) / measuredIntegral;
    }

  // unfold...
  for (Int_t i=0; i<nIterations; i++)
  {
   // calculate Inverse Correlation Map from Bayes theorem:
   // IR_ji = R_ij * prior_i / sum_k(R_kj * prior_k)
   /*Float_t norm = 0;
   for (Int_t me=1; me<=kMaxME; me++)
      for (Int_t mz=1; mz<=kMaxMZ; mz++)
      {
        norm = 0;
        for (Int_t te=kStartBin; te<=kMaxTE; te++)
	  for (Int_t tz=kStartBin; tz<=kMaxTZ; tz++)
	  {
	    Int_t bin[4] = {te, me, tz, mz};
            norm += correlation->GetBinContent(bin)*prior[te][tz];
          }
        if (norm > 0)
          for (Int_t te = kStartBin; te <= kMaxTE; te++)
	    for (Int_t tz = kStartBin; tz <= kMaxTZ; tz++)
	    {
	      Int_t bin[4] = {te, me, tz, mz};
	      inverseCorrelation->SetBinContent(bin, correlation->GetBinContent(bin)*prior[te][tz]/norm );
            }
        //else
          // inverse response set to '0' wich has been already done in line 2069
      }*/
    inverseCorrelation->Reset();     
    Float_t norm[kMaxTE+2][kMaxTZ+2];
    for (Int_t te=0; te<(kMaxTE+2); te++)
      for (Int_t tz=0; tz<(kMaxTZ+2); tz++)
        norm[te][tz]=0;
    for (Int_t idx=0; idx<=correlation->GetNbins(); idx++)
    {
      Int_t bin[4];
      Float_t binContent = correlation->GetBinContent(idx, bin);
      if (bin[1]>0 && bin[1]<=fgkNBINSE && bin[3]>0 && bin[3]<=fgkNBINSZ &&
          bin[0]>0 && bin[0]<=fgkNBINSE && bin[2]>0 && bin[2]<=fgkNBINSZ)
        norm[bin[1]][bin[3]] += binContent*prior[bin[0]][bin[2]];
    }
    Float_t chi2Measured=0, diff;    
    for (Int_t idx=0; idx<=correlation->GetNbins(); idx++)
    {
      Int_t bin[4];
      Float_t binContent = correlation->GetBinContent(idx, bin);
      if (norm[bin[1]][bin[3]]>0 && bin[1]>0 && bin[1]<=fgkNBINSE && 
          bin[3]>0 && bin[3]<=fgkNBINSZ && bin[0]>0 && bin[2]>0 && bin[0]<=fgkNBINSE && bin[2]<=fgkNBINSZ)
      {    
        inverseCorrelation->SetBinContent(bin, binContent*prior[bin[0]][bin[2]]/norm[bin[1]][bin[3]]);
        if (errors[bin[1]][bin[3]]>0)
        {
          diff = ((measuredCopy[bin[1]][bin[3]]-norm[bin[1]][bin[3]])/(errors[bin[1]][bin[3]]));
          chi2Measured += diff*diff;
        }   
      }
    }
    
    // calculate "generated" spectrum
    for (Int_t te = kStartBin; te<=kMaxTE; te++)
      for (Int_t tz = kStartBin; tz<=kMaxTZ; tz++)
      {
        Float_t value = 0;
        for (Int_t me=1; me<=kMaxME; me++)
	  for (Int_t mz=1; mz<=kMaxMZ; mz++)
	  {
	    Int_t bin[4] = {te, me, tz, mz};
	    value += inverseCorrelation->GetBinContent(bin)*measuredCopy[me][mz];
          }
        result[te][tz] = value;
        //printf("%e\n", result[te][tz]);
      }

    // regularization (simple smoothing)
    Float_t chi2LastIter = 0;
    for (Int_t te=kStartBin; te<=kMaxTE; te++)
      for (Int_t tz=kStartBin; tz<=kMaxTZ; tz++)
      {
        Float_t newValue = 0;
        // 0 bin excluded from smoothing
        if (( te >(kStartBin+1) && te<(kMaxTE-1) ) && ( tz > (kStartBin+1) && tz<(kMaxTZ-1) ))
        {
          Float_t average = ((result[te-1][tz-1] + result[te-1][tz] + result[te-1][tz+1])+(result[te][tz-1] + result[te][tz] + result[te][tz+1])+(result[te+1][tz-1] + result[te+1][tz] + result[te+1][tz+1]))/9.;

          // weight the average with the regularization parameter
          newValue = (1 - regPar) * result[te][tz] + regPar * average;
        }
        else
          newValue = result[te][tz];
        if (prior[te][tz]>1.e-5)
        { 
          diff = ((prior[te][tz]-newValue)/prior[te][tz]); 
          chi2LastIter = diff*diff;
        }  
        prior[te][tz] = newValue;
      }
    //printf(" iteration %d - chi2LastIter = %e - chi2Measured = %e \n", i, chi2LastIter/((Float_t)kMaxTE*(Float_t)kMaxTZ), chi2Measured/((Float_t)kMaxTE*(Float_t)kMaxTZ)); 
    if (chi2LastIter/((Float_t)kMaxTE*(Float_t)kMaxTZ)<5.e-6 && chi2Measured/((Float_t)kMaxTE*(Float_t)kMaxTZ)<5.e-3)
      break;
  } // end of iterations
  
  // propagate errors of the reconstructed distribution through the unfolding
  for (Int_t te = kStartBin; te<=kMaxTE; te++)
      for (Int_t tz = kStartBin; tz<=kMaxTZ; tz++)
      {
        Float_t valueError = 0;
	//        Float_t binError = 0;
        for (Int_t me=1; me<=kMaxME; me++)
	  for (Int_t mz=1; mz<=kMaxMZ; mz++)
	  {
	    Int_t bin[4] = {te, me, tz, mz};
	    valueError += inverseCorrelation->GetBinContent(bin)*inverseCorrelation->GetBinContent(bin)*errors[me][mz]*errors[me][mz];
          }
        //if (errors[te][tz]!=0)printf("errors[%d][%d]=%e\n", te, tz, valueError);
        aResult->SetBinContent(te, tz, prior[te][tz]);
        aResult->SetBinError(te, tz, TMath::Sqrt(valueError));   
      }

  // ***********************************************************************************************************
  // Calculate the covariance matrix, all arguments are taken from G. D'Agostini (p.6-8)
  if (calculateErrors)
  {
    printf("Covariance matrix will be calculated... this will take a lot of time (>1 day) ;)\n");
    
    //Variables and Matrices that will be use along the calculation    
    const Int_t binsV[4] = {fgkNBINSE,fgkNBINSE, fgkNBINSZ, fgkNBINSZ};
    const Double_t lowEdgeV[4] = {fgkaxisLowerLimitE, fgkaxisLowerLimitE, fgkaxisLowerLimitZ, fgkaxisLowerLimitZ};
    const Double_t upEdgeV[4] = {fgkaxisUpperLimitE, fgkaxisUpperLimitE, fgkaxisUpperLimitZ, fgkaxisUpperLimitZ};
    
    const Double_t nTrue = (Double_t)measured->Integral();
    
    THnSparseF *v = new THnSparseF("V","",4, binsV, lowEdgeV, upEdgeV);
    v->Reset();            
    Double_t invCorrContent1, nt;
    Double_t invCorrContent2, v11,v12, v2;        
    // calculate V1 and V2
    for (Int_t idx1=0; idx1<=nFillesBins; idx1++)
    {
      printf("Covariance Matrix calculation: iteration idx1=%d of %d\n", idx1, nFillesBins);
      for (Int_t idx2=0; idx2<=nFillesBins; idx2++)
      {      
        Int_t bin1[4];
        Int_t bin2[4];
        invCorrContent1 = inverseCorrelation->GetBinContent(idx1, bin1);
        invCorrContent2 = inverseCorrelation->GetBinContent(idx2, bin2);        
        v11=0; v12=0; v2=0;
        if(bin1[0]>0 && bin1[0]<=fgkNBINSE && bin1[1]>0 && bin1[1]<=fgkNBINSE && 
           bin1[2]>0 && bin1[2]<=fgkNBINSZ && bin1[3]>0 && bin1[3]<=fgkNBINSZ &&
           bin2[0]>0 && bin2[0]<=fgkNBINSE && bin2[1]>0 && bin2[1]<=fgkNBINSE && 
           bin2[2]>0 && bin2[2]<=fgkNBINSZ && bin2[3]>0 && bin2[3]<=fgkNBINSZ)
        {   
          if (bin1[1]==bin2[1] && bin1[3]==bin2[3])    
            v11 = invCorrContent1*invCorrContent2*measuredCopy[bin1[1]][bin1[3]]
                  *(1. - measuredCopy[bin2[1]][bin2[3]]/nTrue);                       
          else
            v12 = invCorrContent1*invCorrContent2*measuredCopy[bin1[1]][bin1[3]]*
                  measuredCopy[bin2[1]][bin2[3]]/nTrue;
          nt = (Double_t)prior[bin2[0]][bin2[2]];
          v2 = measuredCopy[bin1[1]][bin1[3]]*measuredCopy[bin2[1]][bin2[3]]*
               invCorrContent1*invCorrContent2*
               BayesUncertaintyTerms(inverseCorrelation, correlation, bin1, bin2, nt);   
          Int_t binV[4] = {bin1[0],bin2[0],bin1[2],bin2[2]};         
          v->SetBinContent(binV,v11-v12 + v2);
        }  
      }
    }    

    for(Int_t te = 1; te<=fgkNBINSE; te++)
      for(Int_t tz = 1; tz<=fgkNBINSZ; tz++)
      {
        Int_t binV[4] = {te,te,tz,tz}; 
        aResult->SetBinError( te, tz, v->GetBinContent(binV) );
      }            
      
    TFile* f = new TFile("Covariance_UnfSpectrum.root");
    f->Open("RECREATE");
    v->Write();
    f->Close();    
  }  
  
  return 0;

}

//____________________________________________________________________
Double_t AliJetSpectrumUnfolding::BayesUncertaintyTerms(THnSparseF* const M, THnSparseF* const C, Int_t* const binTM, Int_t* const binTM1, Double_t nt)
{
  //
  // helper function for the covariance matrix of the bayesian method
  //

  Double_t result = 0;
  Float_t term[9];
  Int_t tmpBin[4], tmpBin1[4];
  const Int_t nFilledBins = C->GetNbins();
  if (nt==0)
    return 0;
    
  Float_t corrContent;
  Float_t invCorrContent;

  tmpBin[0] =binTM[0]; tmpBin[1] =binTM[1];  tmpBin[2] =binTM[2]; tmpBin[3] =binTM[3];
  tmpBin1[0]=binTM[0]; tmpBin1[1]=binTM1[1]; tmpBin1[2]=binTM[2]; tmpBin1[3]=binTM1[3];      
  if (C->GetBinContent(tmpBin)!=0 && C->GetBinContent(tmpBin1)!=0)
  {
    if (binTM[0]==binTM1[0] && binTM[2]==binTM1[2])
      term[0] = BayesCov(M, C, tmpBin, tmpBin1)/
                (C->GetBinContent(tmpBin)*C->GetBinContent(tmpBin1));
    term[2] = term[0]*M->GetBinContent(tmpBin1);
  }            
  else
  {
    term[0] = 0;
    term[2] = 0;    
  }
              
  tmpBin[0]=binTM1[0]; tmpBin[1]=binTM[1]; tmpBin[2]=binTM1[2]; tmpBin[3]=binTM[3];
  tmpBin1[0]=binTM1[0]; tmpBin1[1]=binTM1[1]; tmpBin1[2]=binTM1[2]; tmpBin1[3]=binTM1[3];      
  if (C->GetBinContent(tmpBin)!=0 && C->GetBinContent(tmpBin1)!=0)
    term[6] = BayesCov(M, C, tmpBin, tmpBin1)*
              M->GetBinContent(tmpBin)/
              (C->GetBinContent(tmpBin)*C->GetBinContent(tmpBin1));                    
  else 
    term[6] = 0;
  
  for(Int_t idx1=0; idx1<=nFilledBins; idx1++)
  { 
    Int_t bin1[4];
    corrContent    = C->GetBinContent(idx1, bin1); 
    invCorrContent = M->GetBinContent(idx1, bin1); 
    if(bin1[0]>0 && bin1[0]<=fgkNBINSE && bin1[1]>0 && bin1[1]<=fgkNBINSE && 
       bin1[2]>0 && bin1[2]<=fgkNBINSZ && bin1[3]>0 && bin1[3]<=fgkNBINSZ)
    {
      tmpBin[0] =binTM[0]; tmpBin[1] =binTM[1]; tmpBin[2] =binTM[2]; tmpBin[3] =binTM[3];
      tmpBin1[0]=binTM[0]; tmpBin1[1]=bin1[1];  tmpBin1[2]=binTM[2]; tmpBin1[3]=bin1[3];      
      if (C->GetBinContent(tmpBin)!=0 &&
          binTM[0]==binTM1[0] && binTM[2]==binTM1[2])
        term[1] = BayesCov(M, C, tmpBin, tmpBin1)/C->GetBinContent(tmpBin);
      else
        term[1] = 0;

      tmpBin[0] =binTM[0]; tmpBin[1] =bin1[1];   tmpBin[2] =binTM[2]; tmpBin[3] =bin1[3];
      tmpBin1[0]=binTM[0]; tmpBin1[1]=binTM1[1]; tmpBin1[2]=binTM[2]; tmpBin1[3]=binTM1[3];      
      if (C->GetBinContent(tmpBin1)!=0)
      {
        if (binTM[0]==binTM1[0] && binTM[2]==binTM1[2])
          term[3] = BayesCov(M, C, tmpBin, tmpBin1)/
                    C->GetBinContent(tmpBin1);
        term[5] = BayesCov(M, C, tmpBin, tmpBin1)*M->GetBinContent(tmpBin1)/
                  C->GetBinContent(tmpBin1);
      }            
      else
      {
        term[3] = 0;
        term[5] = 0;
      }  
   
      tmpBin[0] =binTM1[0]; tmpBin[1] =binTM[1]; tmpBin[2] =binTM1[2]; tmpBin[3] =binTM[3];
      tmpBin1[0]=binTM1[0]; tmpBin1[1]=bin1[1];  tmpBin1[2]=binTM1[2]; tmpBin1[3]=bin1[3];      
      if (C->GetBinContent(tmpBin)!=0)
        term[7] = BayesCov(M, C, tmpBin, tmpBin1)*M->GetBinContent(tmpBin)/
                  C->GetBinContent(tmpBin);
      else
        term[7] = 0;

      tmpBin[0] =bin1[0]; tmpBin[1] =binTM[1];  tmpBin[2] =bin1[2]; tmpBin[3] =binTM[3];
      tmpBin1[0]=bin1[0]; tmpBin1[1]=binTM1[1]; tmpBin1[2]=bin1[2]; tmpBin1[3]=binTM1[3];      
      if (C->GetBinContent(tmpBin)!=0 && C->GetBinContent(tmpBin1)!=0)
        term[8] = BayesCov(M, C, tmpBin, tmpBin1)*
                  M->GetBinContent(tmpBin)*M->GetBinContent(tmpBin)/
                  (C->GetBinContent(tmpBin)*C->GetBinContent(tmpBin1));
      else 
        term[8] = 0;
                       
      for (Int_t i=0; i<9; i++)
        result += term[i]/nt;                    
    }          
  }
   
  return result;
}

//____________________________________________________________________
Double_t AliJetSpectrumUnfolding::BayesCov(THnSparseF* const M, THnSparseF* const correlation, Int_t* const binTM, Int_t* const bin1)
{

  //
  // get the covariance matrix  
  //


  Double_t result, result1, result2, result3;
  
  if (binTM[0]==bin1[0] && binTM[2]==bin1[2])
  {
    if (correlation->GetBinContent(bin1)!=0) 
      result1 = 1./correlation->GetBinContent(bin1);
    else 
      result1 = 0;
    result2 = 1.;
  }
  else
  {
    result1 = 0;
    result2 = 0;
  }  
    
  if (binTM[1]==bin1[1] && binTM[3]==bin1[3])
  {
    Int_t tmpbin[4] = {bin1[0], binTM[1], bin1[2], binTM[3]};
    if(correlation->GetBinContent(tmpbin)!=0)
      result3 = M->GetBinContent(tmpbin)/correlation->GetBinContent(tmpbin);
    else 
      result3 = 0;
  }
  else
  {
    result1 = 0;
    result3 = 0;
  }
    
  return result = result1 + result2 + result3;
}

//____________________________________________________________________
TH2F* AliJetSpectrumUnfolding::CalculateRecSpectrum(TH2* const inputGen)
{
  // runs the distribution given in inputGen through the correlation histogram identified by
  // fCorrelation and produces a reconstructed spectrum

  if (!inputGen)
    return 0;

  // normalize to convert number of events into probability
  /*for (Int_t te=1; te<=NBINSE; te++)
    for (Int_t tz=1; tz<=NBINSZ; tz++)
    {
       Int_t bin[4];
       Float_t sum=0.;
       for (Int_t me = 1; me<=NBINSE; me++)
         for (Int_t mz = 1; mz<=NBINSZ; mz++)
         {
           bin[0] = te; bin[1] = me; 
           bin[2] = tz; bin[3] = mz;           
           sum += fCorrelation[correlationMap]->GetBinContent(bin);
         }
       if (sum > 0.)
         for (Int_t me = 1; me<=NBINSE; me++)
           for (Int_t mz = 1; mz<=NBINSZ; mz++)
           {
             bin[0] = te; bin[1] = me; 
             bin[2] = tz; bin[3] = mz;           
             fCorrelation[correlationMap]->SetBinContent(bin, fCorrelation[correlationMap]->GetBinContent(bin)/sum);
	     fCorrelation[correlationMap]->SetBinError(bin, fCorrelation[correlationMap]->GetBinError(bin)/sum);
	   }
    }*/  
  // normalize to convert number of events into probability (the following loop is much faster)
  Float_t sum[fgkNBINSE+2][fgkNBINSZ+2];
  memset(sum,0,sizeof(Float_t)*(fgkNBINSE+2)*(fgkNBINSZ+2));

  for (Int_t idx=0; idx<fCorrelation->GetNbins(); idx++)
  {
    Int_t bin[4];
    Float_t binContent = fCorrelation->GetBinContent(idx, bin);
    if (bin[1]>0 && bin[1]<=fgkNBINSE && bin[3]>0 && bin[3]<=fgkNBINSZ){
      sum[bin[0]][bin[2]] += binContent; 
    }
  }
  
  for (Int_t idx=0; idx<fCorrelation->GetNbins(); idx++)
  {
    Int_t bin[4];
    Float_t binContent = fCorrelation->GetBinContent(idx, bin);
    Float_t binError   = fCorrelation->GetBinError(bin);
    if (sum[bin[0]][bin[2]]>0 && bin[1]>0 && bin[1]<=fgkNBINSE && 
        bin[3]>0 && bin[3]<=fgkNBINSZ && bin[0]>0 && bin[2]>0 && bin[0]<=fgkNBINSE && bin[2]<=fgkNBINSZ) 
    {
      fCorrelation->SetBinContent(bin, binContent/sum[bin[0]][bin[2]]);
      fCorrelation->SetBinError(bin, binError/sum[bin[0]][bin[2]]); 
    }  
  }

  TH2F* target = dynamic_cast<TH2F*> (fRecSpectrum->Clone(Form("reconstructed_%s", inputGen->GetName())));
  target->Reset();

  for (Int_t me=1; me<=fgkNBINSE; ++me)
    for (Int_t mz=1; mz<=fgkNBINSZ; ++mz)
    {
      Float_t measured = 0;
      Float_t error = 0;

      for (Int_t te=1; te<=fgkNBINSE; ++te)
        for (Int_t tz=1; tz<=fgkNBINSZ; ++tz)
        {
	  Int_t bin[4] = {te, me, tz, mz};
          measured += inputGen->GetBinContent(te,tz) * fCorrelation->GetBinContent(bin);
          error += inputGen->GetBinError(te,tz) * fCorrelation->GetBinContent(bin);
        }
      target->SetBinContent(me, mz, measured);
      target->SetBinError(me, mz, error);
    }

  return target;
}

//__________________________________________________________________________________________________
void AliJetSpectrumUnfolding::SetGenRecFromFunc(TF2* const inputGen)
{
  // uses the given function to fill the input Generated histogram and generates from that
  // the reconstructed histogram by applying the response histogram
  // this can be used to evaluate if the methods work indepedently of the input
  // distribution

  if (!inputGen)
    return;

  TH2F* histtmp = new TH2F("histtmp", "tmp", 
			  fgkNBINSE, fgkaxisLowerLimitE, fgkaxisUpperLimitE,
			  fgkNBINSZ, fgkaxisLowerLimitZ, fgkaxisUpperLimitZ);

  TH2F* gen  = fGenSpectrum;

  histtmp->Reset();
  gen->Reset();

  histtmp->FillRandom(inputGen->GetName(), fgkNEVENTS);

  for (Int_t i=1; i<=gen->GetNbinsX(); ++i)
    for (Int_t j=1; j<=gen->GetNbinsY(); ++j)
    {
      gen->SetBinContent(i, j, histtmp->GetBinContent(i,j));
      gen->SetBinError(i, j, histtmp->GetBinError(i,j));
    }

  delete histtmp;

  //new TCanvas;
  //gStyle->SetPalette(1);
  //gPad->SetLogz();
  //gen->Draw("COLZ");


  TH2 *recsave = fRecSpectrum;

  fRecSpectrum = CalculateRecSpectrum(gen);
  fRecSpectrum->SetName(recsave->GetName());
  delete recsave;

  return;
}
//________________________________________________________________________________________
