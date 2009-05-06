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
// e.g. chi2 minimization and bayesian unfolding
//
//  Author: Jan.Fiete.Grosse-Oetringhaus@cern.ch

#include "AliMultiplicityCorrection.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>
#include <TCollection.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TProfile2D.h>

ClassImp(AliMultiplicityCorrection)

// These are the areas where the quality of the unfolding results are evaluated
// Default defined here, call SetQualityRegions to change them
// unit is in multiplicity (not in bin!)

// SPD:   peak area - flat area - low stat area
Int_t AliMultiplicityCorrection::fgQualityRegionsB[kQualityRegions] = {1,  20, 70};
Int_t AliMultiplicityCorrection::fgQualityRegionsE[kQualityRegions] = {10, 65, 80};

//____________________________________________________________________
void AliMultiplicityCorrection::SetQualityRegions(Bool_t SPDStudy)
{
  //
  // sets the quality region definition to TPC or SPD
  //

  if (SPDStudy)
  {
    // SPD:   peak area - flat area - low stat area
    fgQualityRegionsB[0] = 1;
    fgQualityRegionsE[0] = 10;

    fgQualityRegionsB[1] = 20;
    fgQualityRegionsE[1] = 65;

    fgQualityRegionsB[2] = 70;
    fgQualityRegionsE[2] = 80;

    Printf("AliMultiplicityCorrection::SetQualityRegions --> Enabled quality regions for SPD");
  }
  else
  {
    // TPC:   peak area - flat area - low stat area
    fgQualityRegionsB[0] = 4;
    fgQualityRegionsE[0] = 12;

    fgQualityRegionsB[1] = 25;
    fgQualityRegionsE[1] = 55;

    fgQualityRegionsB[2] = 88;
    fgQualityRegionsE[2] = 108;

    Printf("AliMultiplicityCorrection::SetQualityRegions --> Enabled quality regions for TPC");
  }
}

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection() :
  TNamed(), fCurrentESD(0), fCurrentCorrelation(0), fCurrentEfficiency(0), fLastBinLimit(0), fLastChi2MC(0), fLastChi2MCLimit(0), fLastChi2Residuals(0), fRatioAverage(0)
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
    fMultiplicityNSD[i] = 0;
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = 0;
    fMultiplicityESDCorrected[i] = 0;
  }

  for (Int_t i = 0; i < kQualityRegions; ++i)
    fQuality[i] = 0;
}

//____________________________________________________________________
AliMultiplicityCorrection* AliMultiplicityCorrection::Open(const char* fileName, const char* folderName)
{
  // opens the given file, reads the multiplicity from the given folder and returns the object
  
  TFile* file = TFile::Open(fileName);
  if (!file)
  {
    Printf("ERROR: Could not open %s", fileName);
    return 0;
  }
  
  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection(folderName, folderName);
  mult->LoadHistograms();
  
  // TODO closing the file does not work here, because the histograms cannot be read anymore. LoadHistograms need to be adapted
  
  return mult;
}

//____________________________________________________________________
AliMultiplicityCorrection::AliMultiplicityCorrection(const Char_t* name, const Char_t* title) :
  TNamed(name, title),
  fCurrentESD(0),
  fCurrentCorrelation(0),
  fCurrentEfficiency(0),
  fLastBinLimit(0),
  fLastChi2MC(0),
  fLastChi2MCLimit(0),
  fLastChi2Residuals(0),
  fRatioAverage(0)
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
  #define NBINNING fgkMaxParams, binLimitsN*/

  #define NBINNING 201, -0.5, 200.5
  
  Double_t vtxRange[] = { 15, 6, 2 };

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i] = new TH2F(Form("fMultiplicityESD%d", i), "fMultiplicityESD;vtx-z;Ntracks;Count", 1, -vtxRange[i], vtxRange[i], NBINNING);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i] = dynamic_cast<TH2F*> (fMultiplicityESD[i%3]->Clone(Form("fMultiplicityVtx%d", i)));
    fMultiplicityVtx[i]->SetTitle("fMultiplicityVtx;vtx-z;Npart");

    fMultiplicityMB[i] = dynamic_cast<TH2F*> (fMultiplicityVtx[i]->Clone(Form("fMultiplicityMB%d", i)));
    fMultiplicityMB[i]->SetTitle("fMultiplicityMB");

    fMultiplicityINEL[i] = dynamic_cast<TH2F*> (fMultiplicityVtx[i]->Clone(Form("fMultiplicityINEL%d", i)));
    fMultiplicityINEL[i]->SetTitle("fMultiplicityINEL");
    
    fMultiplicityNSD[i] = dynamic_cast<TH2F*> (fMultiplicityVtx[i]->Clone(Form("fMultiplicityNSD%d", i)));
    fMultiplicityNSD[i]->SetTitle("fMultiplicityNSD");
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    fCorrelation[i] = new TH3F(Form("fCorrelation%d", i), "fCorrelation;vtx-z;Npart;Ntracks", 1, -vtxRange[i%3], vtxRange[i%3], NBINNING, NBINNING);
    fMultiplicityESDCorrected[i] = new TH1F(Form("fMultiplicityESDCorrected%d", i), "fMultiplicityESDCorrected;Npart;dN/dN", NBINNING);
  }

  TH1::AddDirectory(oldStatus);

  AliUnfolding::SetNbins(120, 120);
  AliUnfolding::SetSkipBinsBegin(1);
  AliUnfolding::SetNormalizeInput(kTRUE);
}

//____________________________________________________________________
AliMultiplicityCorrection::~AliMultiplicityCorrection()
{
  //
  // Destructor
  //

  Printf("AliMultiplicityCorrection::~AliMultiplicityCorrection called");

  for (Int_t i = 0; i < kESDHists; ++i)
  {
    if (fMultiplicityESD[i])
      delete fMultiplicityESD[i];
    fMultiplicityESD[i] = 0;
  }

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    if (fMultiplicityVtx[i])
      delete fMultiplicityVtx[i];
    fMultiplicityVtx[i] = 0;

    if (fMultiplicityMB[i])
      delete fMultiplicityMB[i];
    fMultiplicityMB[i] = 0;

    if (fMultiplicityINEL[i])
      delete fMultiplicityINEL[i];
    fMultiplicityINEL[i] = 0;
  
    if (fMultiplicityNSD[i])
      delete fMultiplicityNSD[i];
    fMultiplicityNSD[i] = 0;
}

  for (Int_t i = 0; i < kCorrHists; ++i)
  {
    if (fCorrelation[i])
      delete fCorrelation[i];
    fCorrelation[i] = 0;

    if (fMultiplicityESDCorrected[i])
      delete fMultiplicityESDCorrected[i];
    fMultiplicityESDCorrected[i] = 0;
  }
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
  TList collections[kESDHists+kMCHists*4+kCorrHists*2];

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
      collections[kESDHists+kMCHists*3+i].Add(entry->fMultiplicityNSD[i]);
    }

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists*4+i].Add(entry->fCorrelation[i]);

    for (Int_t i = 0; i < kCorrHists; ++i)
      collections[kESDHists+kMCHists*4+kCorrHists+i].Add(entry->fMultiplicityESDCorrected[i]);

    count++;
  }

  for (Int_t i = 0; i < kESDHists; ++i)
    fMultiplicityESD[i]->Merge(&collections[i]);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    fMultiplicityVtx[i]->Merge(&collections[kESDHists+i]);
    fMultiplicityMB[i]->Merge(&collections[kESDHists+kMCHists+i]);
    fMultiplicityINEL[i]->Merge(&collections[kESDHists+kMCHists*2+i]);
    fMultiplicityNSD[i]->Merge(&collections[kESDHists+kMCHists*3+i]);
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
    fCorrelation[i]->Merge(&collections[kESDHists+kMCHists*4+i]);

  for (Int_t i = 0; i < kCorrHists; ++i)
    fMultiplicityESDCorrected[i]->Merge(&collections[kESDHists+kMCHists*4+kCorrHists+i]);

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

  // store old hists to delete them later
  TList oldObjects;
  oldObjects.SetOwner(1);
  for (Int_t i = 0; i < kESDHists; ++i)
    if (fMultiplicityESD[i])
      oldObjects.Add(fMultiplicityESD[i]);

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    if (fMultiplicityVtx[i])
      oldObjects.Add(fMultiplicityVtx[i]);
    if (fMultiplicityMB[i])
      oldObjects.Add(fMultiplicityMB[i]);
    if (fMultiplicityINEL[i])
      oldObjects.Add(fMultiplicityINEL[i]);
    if (fMultiplicityNSD[i])
      oldObjects.Add(fMultiplicityNSD[i]);
  }

  for (Int_t i = 0; i < kCorrHists; ++i)
    if (fCorrelation[i])
      oldObjects.Add(fCorrelation[i]);

  // load histograms

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
    fMultiplicityNSD[i] = dynamic_cast<TH2F*> (gDirectory->Get(fMultiplicityNSD[i]->GetName()));
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

  // delete old hists
  oldObjects.Delete();

  return success;
}

//____________________________________________________________________
void AliMultiplicityCorrection::SaveHistograms(const char* dir)
{
  //
  // saves the histograms
  //

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  for (Int_t i = 0; i < kESDHists; ++i)
    if (fMultiplicityESD[i])
    {
      fMultiplicityESD[i]->Write();
      fMultiplicityESD[i]->ProjectionY(Form("%s_px", fMultiplicityESD[i]->GetName()), 1, fMultiplicityESD[i]->GetNbinsX())->Write();
    }

  for (Int_t i = 0; i < kMCHists; ++i)
  {
    if (fMultiplicityVtx[i])
    {
      fMultiplicityVtx[i]->Write();
      fMultiplicityVtx[i]->ProjectionY(Form("%s_px", fMultiplicityVtx[i]->GetName()), 1, fMultiplicityVtx[i]->GetNbinsX())->Write();
    }
    if (fMultiplicityMB[i])
    {
      fMultiplicityMB[i]->Write();
      fMultiplicityMB[i]->ProjectionY(Form("%s_px", fMultiplicityMB[i]->GetName()), 1, fMultiplicityMB[i]->GetNbinsX())->Write();
    }
    if (fMultiplicityINEL[i])
    {
      fMultiplicityINEL[i]->Write();
      fMultiplicityINEL[i]->ProjectionY(Form("%s_px", fMultiplicityINEL[i]->GetName()), 1, fMultiplicityINEL[i]->GetNbinsX())->Write();
    }
    if (fMultiplicityNSD[i])
    {
      fMultiplicityNSD[i]->Write();
      fMultiplicityNSD[i]->ProjectionY(Form("%s_px", fMultiplicityNSD[i]->GetName()), 1, fMultiplicityNSD[i]->GetNbinsX())->Write();
    }
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
void AliMultiplicityCorrection::FillGenerated(Float_t vtx, Bool_t triggered, Bool_t vertex, AliPWG0Helper::MCProcessType processType, Int_t generated05, Int_t generated10, Int_t generated14, Int_t generatedAll)
{
  //
  // Fills an event from MC
  //

  if (triggered)
  {
    fMultiplicityMB[0]->Fill(vtx, generated05);
    fMultiplicityMB[1]->Fill(vtx, generated10);
    fMultiplicityMB[2]->Fill(vtx, generated14);
    fMultiplicityMB[3]->Fill(vtx, generatedAll);

    if (vertex)
    {
      fMultiplicityVtx[0]->Fill(vtx, generated05);
      fMultiplicityVtx[1]->Fill(vtx, generated10);
      fMultiplicityVtx[2]->Fill(vtx, generated14);
      fMultiplicityVtx[3]->Fill(vtx, generatedAll);
    }
  }

  fMultiplicityINEL[0]->Fill(vtx, generated05);
  fMultiplicityINEL[1]->Fill(vtx, generated10);
  fMultiplicityINEL[2]->Fill(vtx, generated14);
  fMultiplicityINEL[3]->Fill(vtx, generatedAll);
  
  if (processType != AliPWG0Helper::kSD)
  {
    fMultiplicityNSD[0]->Fill(vtx, generated05);
    fMultiplicityNSD[1]->Fill(vtx, generated10);
    fMultiplicityNSD[2]->Fill(vtx, generated14);
    fMultiplicityNSD[3]->Fill(vtx, generatedAll);
  }
}

//____________________________________________________________________
void AliMultiplicityCorrection::FillMeasured(Float_t vtx, Int_t measured05, Int_t measured10, Int_t measured14)
{
  //
  // Fills an event from ESD
  //

  fMultiplicityESD[0]->Fill(vtx, measured05);
  fMultiplicityESD[1]->Fill(vtx, measured10);
  fMultiplicityESD[2]->Fill(vtx, measured14);
}

//____________________________________________________________________
void AliMultiplicityCorrection::FillCorrection(Float_t vtx, Int_t generated05, Int_t generated10, Int_t generated14, Int_t generatedAll, Int_t measured05, Int_t measured10, Int_t measured14)
{
  //
  // Fills an event into the correlation map with the information from MC and ESD
  //

  fCorrelation[0]->Fill(vtx, generated05, measured05);
  fCorrelation[1]->Fill(vtx, generated10, measured10);
  fCorrelation[2]->Fill(vtx, generated14, measured14);

  fCorrelation[3]->Fill(vtx, generatedAll, measured05);
  fCorrelation[4]->Fill(vtx, generatedAll, measured10);
  fCorrelation[5]->Fill(vtx, generatedAll, measured14);
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
  fMultiplicityESDCorrected[correlationID]->Sumw2();

  // project without under/overflow bins
  fCurrentESD = fMultiplicityESD[inputRange]->ProjectionY("fCurrentESD", 1, fMultiplicityESD[inputRange]->GetXaxis()->GetNbins());
  fCurrentESD->Sumw2();

  // empty under/overflow bins in x, otherwise Project3D takes them into account
  TH3* hist = fCorrelation[correlationID];
  for (Int_t y=0; y<=hist->GetYaxis()->GetNbins()+1; ++y)
  {
    for (Int_t z=0; z<=hist->GetZaxis()->GetNbins()+1; ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  fCurrentCorrelation = (TH2*) hist->Project3D("zy");
  fCurrentCorrelation->Sumw2();
  
  Printf("AliMultiplicityCorrection::SetupCurrentHists: Statistics information: %.f entries in correlation map; %.f entries in measured spectrum", fCurrentCorrelation->Integral(), fCurrentESD->Integral());

#if 0 // does not help
  // null bins with one entry
  Int_t nNulledBins = 0;
  for (Int_t x=1; x<=fCurrentCorrelation->GetXaxis()->GetNbins(); ++x)
    for (Int_t y=1; y<=fCurrentCorrelation->GetYaxis()->GetNbins(); ++y)
    {
      if (fCurrentCorrelation->GetBinContent(x, y) == 1)
      {
        fCurrentCorrelation->SetBinContent(x, y, 0);
        fCurrentCorrelation->SetBinError(x, y, 0);

        ++nNulledBins;
      }
    }
  Printf("Nulled %d bins", nNulledBins);
#endif

  fCurrentEfficiency = GetEfficiency(inputRange, eventType);
  //fCurrentEfficiency->Rebin(2);
  //fCurrentEfficiency->Scale(0.5);
}

//____________________________________________________________________
TH1* AliMultiplicityCorrection::GetEfficiency(Int_t inputRange, EventType eventType)
{
  //
  // calculates efficiency for given event type
  //
  
  TH1* divisor = 0;
  switch (eventType)
  {
    case kTrVtx : break;
    case kMB: divisor = fMultiplicityMB[inputRange]->ProjectionY("divisor", 1, fMultiplicityMB[inputRange]->GetNbinsX(), "e"); break;
    case kINEL: divisor = fMultiplicityINEL[inputRange]->ProjectionY("divisor", 1, fMultiplicityINEL[inputRange]->GetNbinsX(), "e"); break;
    case kNSD: divisor = fMultiplicityNSD[inputRange]->ProjectionY("divisor", 1, fMultiplicityNSD[inputRange]->GetNbinsX(), "e"); break;
  }
  TH1* eff = fMultiplicityVtx[inputRange]->ProjectionY("CurrentEfficiency", 1, fMultiplicityVtx[inputRange]->GetNbinsX(), "e");
  
  if (eventType == kTrVtx)
  {
    for (Int_t i=0; i<= eff->GetNbinsX()+1; i++)
      eff->SetBinContent(i, 1);
  }
  else
    eff->Divide(eff, divisor, 1, 1, "B");
    
  return eff;
}

//____________________________________________________________________
TH1* AliMultiplicityCorrection::GetTriggerEfficiency(Int_t inputRange)
{
  //
  // calculates efficiency for given event type
  //

  TH1* divisor = fMultiplicityINEL[inputRange]->ProjectionY("divisor", 1, fMultiplicityINEL[inputRange]->GetNbinsX(), "e");
  TH1* eff = fMultiplicityMB[inputRange]->ProjectionY("CurrentEfficiency", 1, fMultiplicityMB[inputRange]->GetNbinsX(), "e");
  eff->Divide(eff, divisor, 1, 1, "B");
  return eff;
}

//____________________________________________________________________
Int_t AliMultiplicityCorrection::ApplyMinuitFit(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Bool_t check, TH1* initialConditions)
{
  //
  // correct spectrum using minuit chi2 method
  //
  // for description of parameters, see AliUnfolding::Unfold
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  
  AliUnfolding::SetCreateOverflowBin(5);
  AliUnfolding::SetUnfoldingMethod(AliUnfolding::kChi2Minimization);
  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);
  
  if (!initialConditions)
  {
    initialConditions = (TH1*) fCurrentESD->Clone("initialConditions");
    initialConditions->Scale(1.0 / initialConditions->Integral());
    if (!check)
    {
      // set minimum value to prevent MINUIT just staying in the small value
      for (Int_t i=1; i<=initialConditions->GetNbinsX(); i++)
        initialConditions->SetBinContent(i, TMath::Max(initialConditions->GetBinContent(i), 1e-3));
    }
  }

  return AliUnfolding::Unfold(fCurrentCorrelation, fCurrentEfficiency, fCurrentESD, initialConditions, fMultiplicityESDCorrected[correlationID], check);
}

//____________________________________________________________________
Int_t AliMultiplicityCorrection::ApplyNBDFit(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType)
{
  //
  // correct spectrum using minuit chi2 method with a NBD function
  //
  // for description of parameters, see AliUnfolding::Unfold
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  
  AliUnfolding::SetUnfoldingMethod(AliUnfolding::kFunction);
  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);
  
	TF1* func = new TF1("nbd", "[0] * TMath::Gamma([2]+x) / TMath::Gamma([2]) / TMath::Gamma(x+1) * pow([1] / ([1]+[2]), x) * pow(1.0 + [1]/[2], -[2])");
	func->SetParNames("scaling", "averagen", "k");
	func->SetParLimits(0, 0, 1000);
	func->SetParLimits(1, 1, 50);
	func->SetParLimits(2, 1, 10);
	func->SetParameters(1, 10, 2);
  AliUnfolding::SetFunction(func);
  
  return AliUnfolding::Unfold(fCurrentCorrelation, fCurrentEfficiency, fCurrentESD, 0, fMultiplicityESDCorrected[correlationID]);
}

//____________________________________________________________________
void AliMultiplicityCorrection::DrawHistograms()
{
  //
  // draws the histograms of this class
  //

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
void AliMultiplicityCorrection::DrawComparison(const char* name, Int_t inputRange, Bool_t fullPhaseSpace, Bool_t /*normalizeESD*/, TH1* mcHist, Bool_t simple)
{
  // draw comparison plots


  //mcHist->Rebin(2);

  Int_t esdCorrId = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);

  TString tmpStr;
  tmpStr.Form("%s_DrawComparison_%d", name, esdCorrId);

  if (fMultiplicityESDCorrected[esdCorrId]->Integral() == 0)
  {
    printf("ERROR. Unfolded histogram is empty\n");
    return;
  }

  //regain measured distribution used for unfolding, because the bins at high mult. were modified in SetupCurrentHists
  fCurrentESD = fMultiplicityESD[esdCorrId]->ProjectionY("fCurrentESD", 1, fMultiplicityESD[inputRange]->GetXaxis()->GetNbins());
  fCurrentESD->Sumw2();
  fCurrentESD->Scale(1.0 / fCurrentESD->Integral());

  // normalize unfolded result to 1
  fMultiplicityESDCorrected[esdCorrId]->Scale(1.0 / fMultiplicityESDCorrected[esdCorrId]->Integral());

  //fCurrentESD->Scale(mcHist->Integral(2, 200));

  //new TCanvas;
  /*TH1* ratio = (TH1*) fMultiplicityESDCorrected[esdCorrId]->Clone("ratio");
  ratio->Divide(mcHist);
  ratio->Draw("HIST");
  ratio->Fit("pol0", "W0", "", 20, 120);
  Float_t scalingFactor = ratio->GetFunction("pol0")->GetParameter(0);
  delete ratio;
  mcHist->Scale(scalingFactor);*/

  // find bin with <= 50 entries. this is used as limit for the chi2 calculation
  Int_t mcBinLimit = 0;
  for (Int_t i=20; i<mcHist->GetNbinsX(); ++i)
  {
    if (mcHist->GetBinContent(i) > 50)
    {
      mcBinLimit = i;
    }
    else
      break;
  }
  Printf("AliMultiplicityCorrection::DrawComparison: MC bin limit is %d", mcBinLimit);
  
  // scale to 1
  mcHist->Sumw2();
  if (mcHist->Integral() > 0)
    mcHist->Scale(1.0 / mcHist->Integral());

  // calculate residual

  // for that we convolute the response matrix with the gathered result
  // if normalizeESD is not set, the histogram is already normalized, this needs to be passed to CalculateMultiplicityESD
  TH1* tmpESDEfficiencyRecorrected = (TH1*) fMultiplicityESDCorrected[esdCorrId]->Clone("tmpESDEfficiencyRecorrected");
  tmpESDEfficiencyRecorrected->Multiply(fCurrentEfficiency);
  TH2* convoluted = CalculateMultiplicityESD(tmpESDEfficiencyRecorrected, esdCorrId);
  TH1* convolutedProj = convoluted->ProjectionY("convolutedProj", -1, -1, "e");
  if (convolutedProj->Integral() > 0)
  {
    convolutedProj->Scale(1.0 / convolutedProj->Integral());
  }
  else
    printf("ERROR: convolutedProj is empty. Something went wrong calculating the convoluted histogram.\n");

  TH1* residual = (TH1*) convolutedProj->Clone("residual");
  residual->SetTitle("Residuals;Ntracks;(folded unfolded measured - measured) / e");

  residual->Add(fCurrentESD, -1);
  //residual->Divide(residual, fCurrentESD, 1, 1, "B");

  TH1* residualHist = new TH1F("residualHist", "residualHist", 51, -5, 5);

  // find bin limit
  Int_t lastBin = 0;
  for (Int_t i=1; i<=fCurrentESD->GetNbinsX(); ++i)
  {
    if (fCurrentESD->GetBinContent(i) <= 5)
    {
      lastBin = i;
      break;
    }
  }
  
  // TODO fix errors
  Float_t chi2 = 0;
  for (Int_t i=1; i<=residual->GetNbinsX(); ++i)
  {
    if (fCurrentESD->GetBinError(i) > 0)
    {
      Float_t value = residual->GetBinContent(i) / fCurrentESD->GetBinError(i);
      if (i > 1 && i <= lastBin)
        chi2 += value * value;
      Printf("%d --> %f (%f)", i, value * value, chi2);
      residual->SetBinContent(i, value);
      residualHist->Fill(value);
    }
    else
    {
      //printf("Residual bin %d set to 0\n", i);
      residual->SetBinContent(i, 0);
    }
    convolutedProj->SetBinError(i, 0);
    residual->SetBinError(i, 0);
  }
  fLastChi2Residuals = chi2;

  new TCanvas; residualHist->DrawCopy();

  //residualHist->Fit("gaus", "N");
  //delete residualHist;

  printf("Difference (Residuals) is %f for bin 2-%d\n", fLastChi2Residuals, lastBin);

  TCanvas* canvas1 = 0;
  if (simple)
  {
    canvas1 = new TCanvas(tmpStr, tmpStr, 1200, 600);
    canvas1->Divide(2, 1);
  }
  else
  {
    canvas1 = new TCanvas(tmpStr, tmpStr, 1200, 1200);
    canvas1->Divide(2, 3);
  }

  canvas1->cd(1);
  canvas1->cd(1)->SetGridx();
  canvas1->cd(1)->SetGridy();
  canvas1->cd(1)->SetTopMargin(0.05);
  canvas1->cd(1)->SetRightMargin(0.05);
  canvas1->cd(1)->SetLeftMargin(0.12);
  canvas1->cd(1)->SetBottomMargin(0.12);
  TH1* proj = (TH1*) mcHist->Clone("proj");

  // normalize without 0 bin
  proj->Scale(1.0 / proj->Integral(2, proj->GetNbinsX()));
  Printf("Normalized without 0 bin!");
  proj->GetXaxis()->SetRangeUser(0, 200);
  proj->GetYaxis()->SetTitleOffset(1.4);
  //proj->SetLabelSize(0.05, "xy");
  //proj->SetTitleSize(0.05, "xy");
  proj->SetTitle(Form(";True multiplicity in |#eta| < %.1f;Entries", (inputRange+1)*0.5));
  proj->SetStats(kFALSE);

  fMultiplicityESDCorrected[esdCorrId]->SetLineColor(2);
  fMultiplicityESDCorrected[esdCorrId]->SetMarkerColor(2);
  //fMultiplicityESDCorrected[esdCorrId]->SetMarkerStyle(5);
  // normalize without 0 bin
  fMultiplicityESDCorrected[esdCorrId]->Scale(1.0 / fMultiplicityESDCorrected[esdCorrId]->Integral(2, fMultiplicityESDCorrected[esdCorrId]->GetNbinsX()));
  Printf("Normalized without 0 bin!");

  if (proj->GetEntries() > 0) {
    proj->DrawCopy("HIST");
    fMultiplicityESDCorrected[esdCorrId]->DrawCopy("SAME HIST E");
  }
  else
    fMultiplicityESDCorrected[esdCorrId]->DrawCopy("HIST E");

  gPad->SetLogy();

  TLegend* legend = new TLegend(0.3, 0.8, 0.93, 0.93);
  legend->AddEntry(proj, "True distribution");
  legend->AddEntry(fMultiplicityESDCorrected[esdCorrId], "Unfolded distribution");
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->Draw();
  // unfortunately does not work. maybe a bug? --> legend->SetTextSizePixels(14);

  canvas1->cd(2);
  canvas1->cd(2)->SetGridx();
  canvas1->cd(2)->SetGridy();
  canvas1->cd(2)->SetTopMargin(0.05);
  canvas1->cd(2)->SetRightMargin(0.05);
  canvas1->cd(2)->SetLeftMargin(0.12);
  canvas1->cd(2)->SetBottomMargin(0.12);

  gPad->SetLogy();
  fCurrentESD->GetXaxis()->SetRangeUser(0, 200);
  //fCurrentESD->SetLineColor(2);
  fCurrentESD->SetTitle(Form(";Measured multiplicity in |#eta| < %.1f;Entries", (inputRange+1)*0.5));
  fCurrentESD->SetStats(kFALSE);
  fCurrentESD->GetYaxis()->SetTitleOffset(1.4);
  //fCurrentESD->SetLabelSize(0.05, "xy");
  //fCurrentESD->SetTitleSize(0.05, "xy");
  fCurrentESD->DrawCopy("HIST E");

  convolutedProj->SetLineColor(2);
  convolutedProj->SetMarkerColor(2);
  convolutedProj->SetMarkerStyle(5);
  //proj2->SetMarkerColor(2);
  //proj2->SetMarkerStyle(5);
  convolutedProj->DrawCopy("HIST SAME P");

  legend = new TLegend(0.3, 0.8, 0.93, 0.93);
  legend->AddEntry(fCurrentESD, "Measured distribution");
  legend->AddEntry(convolutedProj, "R #otimes unfolded distribution", "P");
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->Draw();

  //TH1* diffMCUnfolded = dynamic_cast<TH1*> (proj->Clone("diffMCUnfolded"));
  //diffMCUnfolded->Add(fMultiplicityESDCorrected[esdCorrId], -1);

  /*Float_t chi2 = 0;
  Float_t chi = 0;
  fLastChi2MCLimit = 0;
  Int_t limit = 0;
  for (Int_t i=2; i<=200; ++i)
  {
    if (proj->GetBinContent(i) != 0)
    {
      Float_t value = (proj->GetBinContent(i) - fMultiplicityESDCorrected[esdCorrId]->GetBinContent(i)) / proj->GetBinContent(i);
      chi2 += value * value;
      chi += TMath::Abs(value);

      //printf("%d %f\n", i, chi);

      if (chi2 < 0.2)
        fLastChi2MCLimit = i;

      if (chi < 3)
        limit = i;

    }
  }*/

  /*chi2 = 0;
  Float_t chi = 0;
  Int_t limit = 0;
  for (Int_t i=1; i<=diffMCUnfolded->GetNbinsX(); ++i)
  {
    if (fMultiplicityESDCorrected[esdCorrId]->GetBinError(i) > 0)
    {
      Double_t value = diffMCUnfolded->GetBinContent(i) / fMultiplicityESDCorrected[esdCorrId]->GetBinError(i);
      if (value > 1e8)
        value = 1e8; //prevent arithmetic exception
      else if (value < -1e8)
        value = -1e8;
      if (i > 1)
      {
        chi2 += value * value;
        chi += TMath::Abs(value);
      }
      diffMCUnfolded->SetBinContent(i, value);
    }
    else
    {
      //printf("diffMCUnfolded bin %d set to 0\n", i);
      diffMCUnfolded->SetBinContent(i, 0);
    }
    if (chi2 < 1000)
      fLastChi2MCLimit = i;
    if (chi < 1000)
      limit = i;
    if (i == 150)
      fLastChi2MC = chi2;
  }

  printf("limits %d %d\n", fLastChi2MCLimit, limit);
  fLastChi2MCLimit = limit;

  printf("Difference (from MC) is %f for bin 2-150. Limit is %d.\n", fLastChi2MC, fLastChi2MCLimit);*/

  if (!simple)
  {
    /*canvas1->cd(3);

    diffMCUnfolded->SetTitle("#chi^{2};Npart;(MC - Unfolded) / e(unfolded)");
    //diffMCUnfolded->GetYaxis()->SetRangeUser(-20, 20);
    diffMCUnfolded->GetXaxis()->SetRangeUser(0, 200);
    diffMCUnfolded->DrawCopy("HIST");

    TH1F* fluctuation = new TH1F("fluctuation", "fluctuation", 20, -5, 5);
    for (Int_t i=20; i<=diffMCUnfolded->GetNbinsX(); ++i)
      fluctuation->Fill(diffMCUnfolded->GetBinContent(i));

    //new TCanvas; fluctuation->DrawCopy();
    delete fluctuation;*/

    /*TLegend* legend = new TLegend(0.6, 0.7, 0.85, 0.85);
    legend->AddEntry(fMultiplicityESDCorrected, "ESD corrected");
    legend->AddEntry(fMultiplicityMC, "MC");
    legend->AddEntry(fMultiplicityESD, "ESD");
    legend->Draw();*/

    canvas1->cd(4);
    residual->GetYaxis()->SetRangeUser(-5, 5);
    residual->GetXaxis()->SetRangeUser(0, 200);
    residual->DrawCopy();

    canvas1->cd(5);

    TH1* ratio = (TH1*) fMultiplicityESDCorrected[esdCorrId]->Clone("ratio");
    ratio->Divide(mcHist);
    ratio->SetTitle("Ratio;true multiplicity;Unfolded / MC");
    ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio->GetXaxis()->SetRangeUser(0, 200);
    ratio->SetStats(kFALSE);
    ratio->Draw("HIST");

    Double_t ratioChi2 = 0;
    fRatioAverage = 0;
    fLastChi2MCLimit = 0;
    for (Int_t i=2; i<=150; ++i)
    {
      Float_t value = ratio->GetBinContent(i) - 1;
      if (value > 1e8)
        value = 1e8; //prevent arithmetic exception
      else if (value < -1e8)
        value = -1e8;

      ratioChi2 += value * value;
      fRatioAverage += TMath::Abs(value);

      if (ratioChi2 < 0.1)
        fLastChi2MCLimit = i;
    }
    fRatioAverage /= 149;

    printf("Sum over (ratio-1)^2 (2..150) is %f; average of |ratio-1| is %f\n", ratioChi2, fRatioAverage);
    // TODO FAKE
    fLastChi2MC = ratioChi2;

    // FFT of ratio
    canvas1->cd(6);
    const Int_t kFFT = 128;
    Double_t fftReal[kFFT];
    Double_t fftImag[kFFT];

    for (Int_t i=0; i<kFFT; ++i)
    {
      fftReal[i] = ratio->GetBinContent(i+1+10);
      // test: ;-)
      //fftReal[i] = cos(TMath::Pi() * 5 * 2 * i / 128);
      fftImag[i] = 0;
    }

    FFT(-1, TMath::Nint(TMath::Log(kFFT) / TMath::Log(2)), fftReal, fftImag);

    TH1* fftResult = (TH1*) ratio->Clone("fftResult");
    fftResult->SetTitle("FFT;true multiplicity;coeff. (10...137)");
    TH1* fftResult2 = (TH1*) ratio->Clone("fftResult2");
    TH1* fftResult3 = (TH1*) ratio->Clone("fftResult3");
    fftResult->Reset();
    fftResult2->Reset();
    fftResult3->Reset();
    fftResult->GetYaxis()->UnZoom();
    fftResult2->GetYaxis()->UnZoom();
    fftResult3->GetYaxis()->UnZoom();

    //Printf("AFTER FFT");
    for (Int_t i=0; i<kFFT; ++i)
    {
      //Printf("%d: %f", i, fftReal[i]);
      fftResult->SetBinContent(i+1, fftReal[i]);
      /*if (i != 0 && TMath::Abs(fftReal[i]) > 0.5)
      {
        Printf("Nulled %d", i);
        fftReal[i] = 0;
      }*/
    }

    fftResult->SetLineColor(1);
    fftResult->DrawCopy();


    //Printf("IMAG");
    for (Int_t i=0; i<kFFT; ++i)
    {
      //Printf("%d: %f", i, fftImag[i]);
      fftResult2->SetBinContent(i+1, fftImag[i]);

      fftResult3->SetBinContent(i+1, TMath::Sqrt(fftReal[i] * fftReal[i] + fftImag[i] * fftImag[i]));
    }

    fftResult2->SetLineColor(2);
    fftResult2->DrawCopy("SAME");

    fftResult3->SetLineColor(4);
    fftResult3->DrawCopy("SAME");

    for (Int_t i=1; i<kFFT - 1; ++i)
    {
      if (TMath::Sqrt(fftReal[i] * fftReal[i] + fftImag[i] * fftImag[i]) > 3)
      {
        fftReal[i-1] = 0;
        fftReal[i] = 0;
        fftReal[i+1] = 0;
        fftImag[i-1] = 0;
        fftImag[i] = 0;
        fftImag[i+1] = 0;
        //fftReal[i] = (fftReal[i-1] + fftReal[i+1]) / 2;
        //fftImag[i]  = (fftImag[i-1] + fftImag[i+1]) / 2;
        //Printf("Nulled %d to %f %f", i, fftReal[i], fftImag[i]);
      }
    }

    FFT(1, TMath::Nint(TMath::Log(kFFT) / TMath::Log(2)), fftReal, fftImag);

    TH1* fftResult4 = (TH1*) fftResult3->Clone("fftResult4");
    fftResult4->Reset();

    //Printf("AFTER BACK-TRAFO");
    for (Int_t i=0; i<kFFT; ++i)
    {
      //Printf("%d: %f", i, fftReal[i]);
      fftResult4->SetBinContent(i+1+10, fftReal[i]);
    }

    canvas1->cd(5);
    fftResult4->SetLineColor(4);
    fftResult4->DrawCopy("SAME");

    // plot (MC - Unfolded) / error (MC)
    canvas1->cd(3);

    TH1* diffMCUnfolded2 = dynamic_cast<TH1*> (proj->Clone("diffMCUnfolded2"));
    diffMCUnfolded2->Add(fMultiplicityESDCorrected[esdCorrId], -1);

    Int_t ndfQual[kQualityRegions];
    for (Int_t region=0; region<kQualityRegions; ++region)
    {
      fQuality[region] = 0;
      ndfQual[region] = 0;
    }


    Double_t newChi2 = 0;
    Double_t newChi2Limit150 = 0;
    Int_t ndf = 0;
    for (Int_t i=1; i<=diffMCUnfolded2->GetNbinsX(); ++i)
    {
      Double_t value = 0;
      if (proj->GetBinError(i) > 0)
      {
        value = diffMCUnfolded2->GetBinContent(i) / proj->GetBinError(i);
        newChi2 += value * value;
        if (i > 1 && i <= mcBinLimit)
          newChi2Limit150 += value * value;
        ++ndf;

        for (Int_t region=0; region<kQualityRegions; ++region)
          if (diffMCUnfolded2->GetXaxis()->GetBinCenter(i) >= fgQualityRegionsB[region] - 0.1 && diffMCUnfolded2->GetXaxis()->GetBinCenter(i) <= fgQualityRegionsE[region] + 0.1) // 0.1 to avoid e.g. 3.9999 < 4 problem
          {
            fQuality[region] += TMath::Abs(value);
            ++ndfQual[region];
          }
      }

      diffMCUnfolded2->SetBinContent(i, value);
    }

    // normalize region to the number of entries
    for (Int_t region=0; region<kQualityRegions; ++region)
    {
      if (ndfQual[region] > 0)
        fQuality[region] /= ndfQual[region];
      Printf("Quality parameter %d (%d <= mult <= %d) is %f with %d df", region, fgQualityRegionsB[region], fgQualityRegionsE[region], fQuality[region], ndfQual[region]);
    }

    if (mcBinLimit > 1)
    {
      // TODO another FAKE
      fLastChi2MC = newChi2Limit150 / (mcBinLimit - 1);
      Printf("Chi2 (2..%d) from (MC - Unfolded) / e(MC) is: %.2f ndf is %d --> chi2 / ndf = %.2f", mcBinLimit, newChi2Limit150, mcBinLimit - 1, fLastChi2MC);
    }
    else
      fLastChi2MC = -1;

    Printf("Chi2 (full range) from (MC - Unfolded) / e(MC) is: %.2f ndf is %d --> chi2 / ndf = %.2f", newChi2, ndf, ((ndf > 0) ? newChi2 / ndf : -1));

    diffMCUnfolded2->SetTitle("#chi^{2};Npart;(MC - Unfolded) / e(MC)");
    diffMCUnfolded2->GetYaxis()->SetRangeUser(-5, 5);
    diffMCUnfolded2->GetXaxis()->SetRangeUser(0, 200);
    diffMCUnfolded2->DrawCopy("HIST");
  }

  canvas1->SaveAs(Form("%s.gif", canvas1->GetName()));
}

//____________________________________________________________________
void AliMultiplicityCorrection::FFT(Int_t dir, Int_t m, Double_t *x, Double_t *y)
{
/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/

   Long_t   nn, i, i1, j, k, i2, l, l1, l2;
   Double_t c1, c2, tx, ty, t1, t2, u1, u2, z;

   /* Calculate the number of points */
   nn = 1;
   for (i = 0; i < m; i++)
       nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i= 0; i < nn - 1; i++) {
       if (i < j) {
	   tx = x[i];
	   ty = y[i];
	   x[i] = x[j];
	   y[i] = y[j];
	   x[j] = tx;
	   y[j] = ty;
       }
       k = i2;
       while (k <= j) {
	   j -= k;
	   k >>= 1;
       }
       j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l = 0; l < m; l++) {
       l1 = l2;
       l2 <<= 1;
       u1 = 1.0;
       u2 = 0.0;
       for (j = 0;j < l1; j++) {
	   for (i = j; i < nn; i += l2) {
	       i1 = i + l1;
	       t1 = u1 * x[i1] - u2 * y[i1];
	       t2 = u1 * y[i1] + u2 * x[i1];
	       x[i1] = x[i] - t1;
	       y[i1] = y[i] - t2;
	       x[i] += t1;
	       y[i] += t2;
	   }
	   z =  u1 * c1 - u2 * c2;
	   u2 = u1 * c2 + u2 * c1;
	   u1 = z;
       }
       c2 = TMath::Sqrt((1.0 - c1) / 2.0);
       if (dir == 1)
	   c2 = -c2;
       c1 = TMath::Sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
       for (i=0;i<nn;i++) {
	   x[i] /= (Double_t)nn;
	   y[i] /= (Double_t)nn;
       }
   }
}

//____________________________________________________________________
void AliMultiplicityCorrection::GetComparisonResults(Float_t* mc, Int_t* mcLimit, Float_t* residuals, Float_t* ratioAverage) const
{
  // Returns the chi2 between the MC and the unfolded ESD as well as between the ESD and the folded unfolded ESD
  // These values are computed during DrawComparison, thus this function picks up the
  // last calculation

  if (mc)
    *mc = fLastChi2MC;
  if (mcLimit)
    *mcLimit = fLastChi2MCLimit;
  if (residuals)
    *residuals = fLastChi2Residuals;
  if (ratioAverage)
    *ratioAverage = fRatioAverage;
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
    case kNSD : return fMultiplicityNSD[i]; break;
  }

  return 0;
}

//____________________________________________________________________
void AliMultiplicityCorrection::SetMultiplicityMC(Int_t i, EventType eventType, TH2F* hist)
{
  //
  // returns the corresponding MC spectrum
  //

  switch (eventType)
  {
    case kTrVtx : fMultiplicityVtx[i] = hist; break;
    case kMB : fMultiplicityMB[i] = hist; break;
    case kINEL : fMultiplicityINEL[i] = hist; break;
    case kNSD : fMultiplicityNSD[i] = hist; break;
  }
}

//____________________________________________________________________
TH1* AliMultiplicityCorrection::CalculateStdDev(TH1** results, Int_t max)
{
  // calculate standard deviation of (results[0] - results[k]) k=1...max-1
  // per bin one gets: sigma(r[0] - r[n]) / r[0]

  TH1* standardDeviation = (TH1*) results[0]->Clone("standardDeviation");
  standardDeviation->Reset();

  for (Int_t x=1; x<=results[0]->GetNbinsX(); x++)
  {
    if (results[0]->GetBinContent(x) > 0)
    {
      Double_t average = 0;
      for (Int_t n=1; n<max; ++n)
        average += results[n]->GetBinContent(x) - results[0]->GetBinContent(x);
      average /= max-1;

      Double_t variance = 0;
      for (Int_t n=1; n<max; ++n)
      {
        Double_t value = results[n]->GetBinContent(x) - results[0]->GetBinContent(x) - average;
        variance += value * value;
      }
      variance /= max-1;

      Double_t standardDev = TMath::Sqrt(variance);
      standardDeviation->SetBinContent(x, standardDev / results[0]->GetBinContent(x));
      //Printf("sigma_%d is %f value %f --> error %f", x, standardDev, results[0]->GetBinContent(x), standardDev / results[0]->GetBinContent(x));
    }
  }

  return standardDeviation;
}

//____________________________________________________________________
TH1* AliMultiplicityCorrection::StatisticalUncertainty(AliUnfolding::MethodType methodType, Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Bool_t randomizeMeasured, Bool_t randomizeResponse, TH1* compareTo)
{
  //
  // evaluates the uncertainty that arises from the non-infinite statistics in the response matrix
  // the function unfolds the spectrum using the default response matrix and several modified ones
  // the modified ones are created by randomizing each cell using poisson statistics with the mean = bin value
  // these unfolded results are compared to the first result gained with the default response OR to the histogram given
  // in <compareTo> (optional)
  //
  // returns the error assigned to the measurement
  //

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);

  // initialize seed with current time
  gRandom->SetSeed(0);
  
  if (methodType == AliUnfolding::kChi2Minimization)
    AliUnfolding::SetCreateOverflowBin(5);
  AliUnfolding::SetUnfoldingMethod(methodType);

  const Int_t kErrorIterations = 150;

  TH1* maxError = 0;
  TH1* firstResult = 0;

  TH1** results = new TH1*[kErrorIterations];

  for (Int_t n=0; n<kErrorIterations; ++n)
  {
    Printf("Iteration %d of %d...", n, kErrorIterations);

    SetupCurrentHists(inputRange, fullPhaseSpace, eventType);

    TH1* measured = (TH1*) fCurrentESD->Clone("measured");

    if (n > 0)
    {
      if (randomizeResponse)
      {
        // randomize response matrix
        for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
          for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
            fCurrentCorrelation->SetBinContent(i, j, gRandom->Poisson(fCurrentCorrelation->GetBinContent(i, j)));
      }

      if (randomizeMeasured)
      {
        // randomize measured spectrum
        for (Int_t x=1; x<=measured->GetNbinsX(); x++) // mult. axis
        {
          Int_t randomValue = gRandom->Poisson(fCurrentESD->GetBinContent(x));
          measured->SetBinContent(x, randomValue);
          measured->SetBinError(x, TMath::Sqrt(randomValue));
        }
      }
    }

    // only for bayesian method we have to do it before the call to Unfold...
    if (methodType == AliUnfolding::kBayesian)
    {
      for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
      {
        // with this it is normalized to 1
        Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());

        // with this normalized to the given efficiency
        if (fCurrentEfficiency->GetBinContent(i) > 0)
          sum /= fCurrentEfficiency->GetBinContent(i);
        else
          sum = 0;

        for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
        {
          if (sum > 0)
          {
            fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
            fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
          }
          else
          {
            fCurrentCorrelation->SetBinContent(i, j, 0);
            fCurrentCorrelation->SetBinError(i, j, 0);
          }
        }
      }
    }

    TH1* result = 0;
    if (n == 0 && compareTo)
    {
      // in this case we just store the histogram we want to compare to
      result = (TH1*) compareTo->Clone("compareTo");
      result->Sumw2();
    }
    else
    {
      result = (TH1*) fMultiplicityESDCorrected[correlationID]->Clone(Form("result_%d", n));

      Int_t returnCode = AliUnfolding::Unfold(fCurrentCorrelation, fCurrentEfficiency, measured, 0, result);

      if (returnCode != 0)
        return 0;
    }

    // normalize
    result->Scale(1.0 / result->Integral());

    if (n == 0)
    {
      firstResult = (TH1*) result->Clone("firstResult");

      maxError = (TH1*) result->Clone("maxError");
      maxError->Reset();
    }
    else
    {
      // calculate ratio
      TH1* ratio = (TH1*) firstResult->Clone("ratio");
      ratio->Divide(result);

      // find max. deviation
      for (Int_t x=1; x<=ratio->GetNbinsX(); x++)
        maxError->SetBinContent(x, TMath::Max(maxError->GetBinContent(x), TMath::Abs(1 - ratio->GetBinContent(x))));

      delete ratio;
    }

    results[n] = result;
  }

  // find covariance matrix
  // results[n] is X_x
  // cov. matrix is M_xy = E ( (X_x - E(X_x)) * (X_y - E(X_y))), with E() = expectation value

  Int_t nBins = results[0]->GetNbinsX();
  Float_t lowEdge = results[0]->GetXaxis()->GetBinLowEdge(1);
  Float_t upEdge = results[0]->GetXaxis()->GetBinUpEdge(nBins);

  // find average, E(X)
  TProfile* average = new TProfile("average", "average", nBins, lowEdge, upEdge);
  for (Int_t n=1; n<kErrorIterations; ++n)
    for (Int_t x=1; x<=results[n]->GetNbinsX(); x++)
      average->Fill(results[n]->GetXaxis()->GetBinCenter(x), results[n]->GetBinContent(x));
  //new TCanvas; average->DrawClone();
  
  // find cov. matrix
  TProfile2D* covMatrix = new TProfile2D("covMatrix", "covMatrix", nBins, lowEdge, upEdge, nBins, lowEdge, upEdge);

  for (Int_t n=1; n<kErrorIterations; ++n)
    for (Int_t x=1; x<=results[n]->GetNbinsX(); x++)
      for (Int_t y=1; y<=results[n]->GetNbinsX(); y++)
      {
        // (X_x - E(X_x)) * (X_y - E(X_y)
        Float_t cov = (results[n]->GetBinContent(x) - average->GetBinContent(x)) * (results[n]->GetBinContent(y) - average->GetBinContent(y));
        covMatrix->Fill(results[n]->GetXaxis()->GetBinCenter(x), results[n]->GetXaxis()->GetBinCenter(y), cov);
      }
  TCanvas* c = new TCanvas; c->cd(); covMatrix->DrawCopy("COLZ");

//   // fill 2D histogram that contains deviation from first
//   TH2F* deviations = new TH2F("deviations", "deviations", nBins, lowEdge, upEdge, 1000, -0.01, 0.01);
//   for (Int_t n=1; n<kErrorIterations; ++n)
//     for (Int_t x=1; x<=results[n]->GetNbinsX(); x++)
//       deviations->Fill(results[n]->GetXaxis()->GetBinCenter(x), results[n]->GetBinContent(x) - results[0]->GetBinContent(x));
//   //new TCanvas; deviations->DrawCopy("COLZ");
// 
//   // get standard deviation "by hand"
//   for (Int_t x=1; x<=nBins; x++)
//   {
//     if (results[0]->GetBinContent(x) > 0)
//     {
//       TH1* proj = deviations->ProjectionY("projRMS", x, x, "e");
//       Float_t standardDev = proj->GetRMS(); // this is standard deviation in fact
//       //standardDeviation->SetBinContent(x, standardDev / results[0]->GetBinContent(x));
//       Printf("sigma_%d is %f value %f --> error %f", x, standardDev, results[0]->GetBinContent(x), standardDev / results[0]->GetBinContent(x));
//     }
//   }

  TH1* standardDeviation = CalculateStdDev(results, kErrorIterations);

  // compare maxError to RMS of profile (variable name: average)
  // first: calculate rms in percent of value
  TH1* rmsError = (TH1*) maxError->Clone("rmsError");
  rmsError->Reset();

  // enable error to be standard deviation (see http://root.cern.ch/root/html/TProfile.html#TProfile:SetErrorOption)
  average->SetErrorOption("s");
  for (Int_t x=1; x<=rmsError->GetNbinsX(); x++)
    if (average->GetBinContent(x) > 0)
      rmsError->SetBinContent(x, average->GetBinError(x) / average->GetBinContent(x));

  // find maxError deviation from average (not from "first result"/mc as above)
  TH1* maxError2 = (TH1*) maxError->Clone("maxError2");
  maxError2->Reset();
  for (Int_t n=1; n<kErrorIterations; ++n)
    for (Int_t x=1; x<=results[n]->GetNbinsX(); x++)
      if (average->GetBinContent(x) > 0)
        maxError2->SetBinContent(x, TMath::Max(maxError2->GetBinContent(x), TMath::Abs((results[n]->GetBinContent(x) - average->GetBinContent(x)) / average->GetBinContent(x))));

  //new TCanvas; maxError2->DrawCopy(); rmsError->SetLineColor(2); rmsError->DrawCopy("SAME"); standardDeviation->SetLineColor(3); standardDeviation->DrawCopy("SAME");

  // plot difference between average and MC/first result
  TH1* averageFirstRatio = (TH1*) results[0]->Clone("averageFirstRatio");
  averageFirstRatio->Reset();
  averageFirstRatio->Divide(results[0], average);

  //new TCanvas; results[0]->DrawCopy(); average->SetLineColor(2); average->DrawClone("SAME");
  //new TCanvas; averageFirstRatio->DrawCopy();

  static TH1* temp = 0;
  if (!temp)
  {
    temp = (TH1*) standardDeviation->Clone("temp");
    for (Int_t x=1; x<=results[0]->GetNbinsX(); x++)
      temp->SetBinContent(x, temp->GetBinContent(x) * results[0]->GetBinContent(x));
  }
  else
  {
    // find difference from result[0] as TH2
    TH2F* pulls = new TH2F("pulls", "pulls;multiplicity;difference", nBins, lowEdge, upEdge, 1000, -10, 10);
    for (Int_t n=1; n<kErrorIterations; ++n)
      for (Int_t x=1; x<=results[n]->GetNbinsX(); x++)
        if (temp->GetBinContent(x) > 0)
          pulls->Fill(results[n]->GetXaxis()->GetBinCenter(x), (results[0]->GetBinContent(x) - results[n]->GetBinContent(x)) / temp->GetBinContent(x));
    new TCanvas("pulls", "pulls", 800, 600); pulls->DrawCopy(); pulls->FitSlicesY();
  }

  // clean up
  for (Int_t n=0; n<kErrorIterations; ++n)
    delete results[n];
  delete[] results;

  // fill into result histogram
  for (Int_t i=1; i<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); ++i)
    fMultiplicityESDCorrected[correlationID]->SetBinContent(i, firstResult->GetBinContent(i));

  for (Int_t i=1; i<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); ++i)
    fMultiplicityESDCorrected[correlationID]->SetBinError(i, maxError->GetBinContent(i) * fMultiplicityESDCorrected[correlationID]->GetBinContent(i));

  return standardDeviation;
}

//____________________________________________________________________
void AliMultiplicityCorrection::ApplyBayesianMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType, Float_t regPar, Int_t nIterations, TH1* initialConditions, Bool_t determineError)
{
  //
  // correct spectrum using bayesian method
  //

  // initialize seed with current time
  gRandom->SetSeed(0);

  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);

  // normalize correction for given nPart
  for (Int_t i=1; i<=fCurrentCorrelation->GetNbinsX(); ++i)
  {
    // with this it is normalized to 1
    Double_t sum = fCurrentCorrelation->Integral(i, i, 1, fCurrentCorrelation->GetNbinsY());

    // with this normalized to the given efficiency
    if (fCurrentEfficiency->GetBinContent(i) > 0)
      sum /= fCurrentEfficiency->GetBinContent(i);
    else
      sum = 0;

    for (Int_t j=1; j<=fCurrentCorrelation->GetNbinsY(); ++j)
    {
      if (sum > 0)
      {
        fCurrentCorrelation->SetBinContent(i, j, fCurrentCorrelation->GetBinContent(i, j) / sum);
        fCurrentCorrelation->SetBinError(i, j, fCurrentCorrelation->GetBinError(i, j) / sum);
      }
      else
      {
        fCurrentCorrelation->SetBinContent(i, j, 0);
        fCurrentCorrelation->SetBinError(i, j, 0);
      }
    }
  }

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);

  AliUnfolding::SetBayesianParameters(regPar, nIterations);
  AliUnfolding::SetUnfoldingMethod(AliUnfolding::kBayesian);
  if (AliUnfolding::Unfold(fCurrentCorrelation, fCurrentEfficiency, fCurrentESD, initialConditions, fMultiplicityESDCorrected[correlationID]) != 0)
    return;

  if (!determineError)
  {
    Printf("AliMultiplicityCorrection::ApplyBayesianMethod: WARNING: No errors calculated.");
    return;
  }

  // evaluate errors, this is done by randomizing the measured spectrum following Poission statistics
  // this (new) measured spectrum is then unfolded and the different to the result from the "real" measured
  // spectrum calculated. This is performed N times and the sigma is taken as the statistical
  // error of the unfolding method itself.

  const Int_t kErrorIterations = 20;

  printf("Spectrum unfolded. Determining error (%d iterations)...\n", kErrorIterations);

  TH1* randomized = (TH1*) fCurrentESD->Clone("randomized");
  TH1* resultArray[kErrorIterations+1];
  for (Int_t n=0; n<kErrorIterations; ++n)
  {
    // randomize the content of clone following a poisson with the mean = the value of that bin
    for (Int_t x=1; x<=randomized->GetNbinsX(); x++) // mult. axis
    {
      Int_t randomValue = gRandom->Poisson(fCurrentESD->GetBinContent(x));
      //printf("%d --> %d\n", fCurrentESD->GetBinContent(x), randomValue);
      randomized->SetBinContent(x, randomValue);
      randomized->SetBinError(x, TMath::Sqrt(randomValue));
    }

    TH1* result2 = (TH1*) fMultiplicityESDCorrected[correlationID]->Clone("result2");
    result2->Reset();
    if (AliUnfolding::Unfold(fCurrentCorrelation, fCurrentEfficiency, randomized, initialConditions, result2) != 0)
      return;

    result2->Scale(1.0 / result2->Integral());

    resultArray[n+1] = result2;
  }
  delete randomized;

  resultArray[0] = fMultiplicityESDCorrected[correlationID];
  TH1* error = CalculateStdDev(resultArray, kErrorIterations+1);

  for (Int_t n=0; n<kErrorIterations; ++n)
    delete resultArray[n+1];

  for (Int_t i=1; i<=fMultiplicityESDCorrected[correlationID]->GetNbinsX(); ++i)
    fMultiplicityESDCorrected[correlationID]->SetBinError(i, error->GetBinContent(i) * fMultiplicityESDCorrected[correlationID]->GetBinContent(i));

  delete error;
}

//____________________________________________________________________
Float_t AliMultiplicityCorrection::BayesCovarianceDerivate(Float_t matrixM[251][251], TH2* hResponse, Int_t k, Int_t i, Int_t r, Int_t u)
{
  //
  // helper function for the covariance matrix of the bayesian method
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
void AliMultiplicityCorrection::ApplyLaszloMethod(Int_t inputRange, Bool_t fullPhaseSpace, EventType eventType)
{
  //
  // correct spectrum using bayesian method
  //

  Float_t regPar = 0;

  Int_t correlationID = inputRange + ((fullPhaseSpace == kFALSE) ? 0 : 4);
  Int_t mcTarget = ((fullPhaseSpace == kFALSE) ? inputRange : 4);

  SetupCurrentHists(inputRange, fullPhaseSpace, eventType);
  //normalize ESD
  fCurrentESD->Scale(1.0 / fCurrentESD->Integral());

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
  Int_t iterations = 25;
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
    norm = 1; //hTrueSmoothed->Integral();

    for (Int_t t=1; t<hTrueSmoothed->GetNbinsX(); t++)
    {
      Float_t newValue = hTemp->GetBinContent(t)/norm;
      Float_t diff = hPrior->GetBinContent(t) - newValue;
      chi2 += (Double_t) diff * diff;

      hPrior->SetBinContent(t, newValue);
    }

    printf("Chi2 of %d iteration = %.10f\n", i, chi2);

    //if (i % 5 == 0)
      DrawComparison(Form("Laszlo_%d", i), inputRange, fullPhaseSpace, kTRUE, GetMultiplicityMC(mcTarget, eventType)->ProjectionY());

    delete hTrueSmoothed;
  } // end of iterations

  DrawComparison("Laszlo", inputRange, fullPhaseSpace, kTRUE, GetMultiplicityMC(mcTarget, eventType)->ProjectionY());
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
  //normalize ESD
  fCurrentESD->Scale(1.0 / fCurrentESD->Integral());

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
TH2F* AliMultiplicityCorrection::CalculateMultiplicityESD(TH1* inputMC, Int_t correlationMap)
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
  for (Int_t y=0; y<=hist->GetYaxis()->GetNbins()+1; ++y)
  {
    for (Int_t z=0; z<=hist->GetZaxis()->GetNbins()+1; ++z)
    {
      hist->SetBinContent(0, y, z, 0);
      hist->SetBinContent(hist->GetXaxis()->GetNbins()+1, y, z, 0);
    }
  }

  TH2* corr = (TH2*) hist->Project3D("zy");
  //corr->Rebin2D(2, 1);
  corr->Sumw2();
  Printf("Correction histogram used for convolution has %f entries", corr->Integral());

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

    for (Int_t gen=1; gen<=corr->GetNbinsX(); ++gen)
    {
      Int_t mcGenBin = inputMC->GetXaxis()->FindBin(corr->GetXaxis()->GetBinCenter(gen));

      measured += inputMC->GetBinContent(mcGenBin) * corr->GetBinContent(gen, meas);
      error += inputMC->GetBinError(mcGenBin) * corr->GetBinContent(gen, meas);
    }

    //printf("%f +- %f ; %f +- %f \n", inputMC->GetBinContent(meas), inputMC->GetBinError(meas), measured, error);

    target->SetBinContent(1 + target->GetNbinsX() / 2, meas, measured);
    target->SetBinError(1 + target->GetNbinsX() / 2, meas, error);
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

  // fill histogram used for random generation
  TH1* tmp = fMultiplicityVtx[id]->ProjectionY("tmp");
  tmp->Reset();

  for (Int_t i=1; i<=tmp->GetNbinsX(); ++i)
    tmp->SetBinContent(i, inputMC->Eval(tmp->GetXaxis()->GetBinCenter(i)) * tmp->GetXaxis()->GetBinWidth(i));
    
  TH1* mcRnd = fMultiplicityVtx[id]->ProjectionY("mcRnd");
  mcRnd->Reset();
  mcRnd->FillRandom(tmp, tmp->Integral());
  
  //new TCanvas; tmp->Draw();
  //new TCanvas; mcRnd->Draw();
  
  // and move into 2d histogram
  TH1* mc = fMultiplicityVtx[id];
  mc->Reset();
  for (Int_t i=1; i<=mc->GetNbinsY(); ++i)
  {
    mc->SetBinContent(mc->GetNbinsX() / 2 + 1, i, mcRnd->GetBinContent(i));
    mc->SetBinError(mc->GetNbinsX() / 2 + 1, i, TMath::Sqrt(mcRnd->GetBinContent(i)));
  }
  
  //new TCanvas; mc->Draw("COLZ");

  // now randomize the measured histogram; funcMeasured is used as pilot function to generated the measured entries
  TH1* funcMeasured = CalculateMultiplicityESD(tmp, id)->ProjectionY("funcMeasured");
  
  //new TCanvas; funcMeasured->Draw();
  
  fMultiplicityESD[id]->Reset();
  
  TH1* measRnd = fMultiplicityESD[id]->ProjectionY("measRnd");
  measRnd->FillRandom(funcMeasured, tmp->Integral());
  
  //new TCanvas; measRnd->Draw();
  
  fMultiplicityESD[id]->Reset();
  for (Int_t i=1; i<=fMultiplicityESD[id]->GetNbinsY(); ++i)
  {
    fMultiplicityESD[id]->SetBinContent(fMultiplicityESD[id]->GetNbinsX() / 2 + 1, i, measRnd->GetBinContent(i));
    fMultiplicityESD[id]->SetBinError(fMultiplicityESD[id]->GetNbinsX() / 2 + 1, i, TMath::Sqrt(measRnd->GetBinContent(i)));
  }
}
