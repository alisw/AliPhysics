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

/* $Id: AliUEHist.cxx 20164 2007-08-14 15:31:50Z morsch $ */

//
//
// encapsulate histogram and corrections for one underlying event histogram
//
//
// Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

#include "AliUEHist.h"
#include "AliCFContainer.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TList.h"
#include "TCollection.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliLog.h"
#include "TCanvas.h"
#include "TF1.h"

ClassImp(AliUEHist)

const Int_t AliUEHist::fgkCFSteps = 10;

AliUEHist::AliUEHist(const char* reqHist) : 
  TObject(),
  fkRegions(4),
  fEventHist(0),
  fTrackHistEfficiency(0),
  fEtaMin(0),
  fEtaMax(0),
  fPtMin(0),
  fPtMax(0),
  fContaminationEnhancement(0),
  fCombineMinMax(0),
  fCache(0)
{
  // Constructor
    
  for (Int_t i=0; i<fkRegions; i++)
    fTrackHist[i] = 0;
    
  if (strlen(reqHist) == 0)
    return;
    
  const char* title = "";
    
  // track level
  Int_t nTrackVars = 4; // eta vs pT vs pT,lead (vs delta phi) vs multiplicity
  Int_t iTrackBin[5];
  Double_t* trackBins[5];
  const char* trackAxisTitle[5];
  
  // eta
  iTrackBin[0] = 20;
  Double_t etaBins[20+1];
  for (Int_t i=0; i<=iTrackBin[0]; i++)
    etaBins[i] = -1.0 + 0.1 * i;
  trackBins[0] = etaBins;
  trackAxisTitle[0] = "#eta";
  
  // pT
  iTrackBin[1] = 39;
  Double_t pTBins[] = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
  trackBins[1] = pTBins;
  trackAxisTitle[1] = "p_{T} (GeV/c)";
  
  // pT,lead binning 1
  const Int_t kNLeadingpTBins = 100;
  Double_t leadingpTBins[kNLeadingpTBins+1];
  for (Int_t i=0; i<=kNLeadingpTBins; i++)
    leadingpTBins[i] = 0.5 * i;
  
  // pT,lead binning 2
  const Int_t kNLeadingpTBins2 = 13;
  Double_t leadingpTBins2[] = { 0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0 };
  
  // phi,lead
  const Int_t kNLeadingPhiBins = 40;
  Double_t leadingPhiBins[kNLeadingPhiBins+1];
  for (Int_t i=0; i<=kNLeadingPhiBins; i++)
    leadingPhiBins[i] = -0.5 * TMath::Pi() + 1.0 / 40 * i * TMath::TwoPi();
    
  // multiplicity
  const Int_t kNMultiplicityBins = 15;
  Double_t multiplicityBins[kNMultiplicityBins+1];
  for (Int_t i=0; i<=kNMultiplicityBins; i++)
    multiplicityBins[i] = -0.5 + i;
  multiplicityBins[kNMultiplicityBins] = 200;

  // particle species
  const Int_t kNSpeciesBins = 4; // pi, K, p, rest
  Double_t speciesBins[] = { -0.5, 0.5, 1.5, 2.5, 3.5 };
  
  trackBins[3] = multiplicityBins;
  iTrackBin[3] = kNMultiplicityBins;
  trackAxisTitle[3] = "multiplicity";
  
  // selection depending on requested histogram
  Int_t axis = -1; // 0 = pT,lead, 1 = phi,lead
  if (strcmp(reqHist, "NumberDensitypT") == 0)
  {
    axis = 0;
    title = "d^{2}N_{ch}/d#phid#eta";
  }
  else if (strcmp(reqHist, "NumberDensityPhi") == 0)
  {
    axis = 1;
    title = "d^{2}N_{ch}/d#phid#eta";
  }
  else if (strcmp(reqHist, "SumpT") == 0)
  {
    axis = 0;
    title = "d^{2}#Sigma p_{T}/d#phid#eta";
  }
  else
    AliFatal(Form("Invalid histogram requested: %s", reqHist));
  
  Int_t initRegions = fkRegions;
  
  if (axis == 0)
  {
    trackBins[2] = leadingpTBins;
    iTrackBin[2] = kNLeadingpTBins;
    trackAxisTitle[2] = "leading p_{T} (GeV/c)";
    
  }
  else if (axis == 1)
  {
    nTrackVars = 5;
    initRegions = 1;
  
    iTrackBin[2] = kNLeadingpTBins2;
    trackBins[2] = leadingpTBins2;
    trackAxisTitle[2] = "leading p_{T} (GeV/c)";
    
    iTrackBin[4] = kNLeadingPhiBins;
    trackBins[4] = leadingPhiBins;
    trackAxisTitle[4] = "#phi w.r.t leading track";
  }
    
  for (Int_t i=0; i<initRegions; i++)
  {
    fTrackHist[i] = new AliCFContainer(Form("fTrackHist_%d", i), title, fgkCFSteps, nTrackVars, iTrackBin);
    
    for (Int_t j=0; j<nTrackVars; j++)
    {
      fTrackHist[i]->SetBinLimits(j, trackBins[j]);
      fTrackHist[i]->SetVarTitle(j, trackAxisTitle[j]);
    }
    
    SetStepNames(fTrackHist[i]);
  }
  
  // track 3rd and 4th axis --> event 1st and 2nd axis
  fEventHist = new AliCFContainer("fEventHist", title, fgkCFSteps, 2, iTrackBin+2);
  
  fEventHist->SetBinLimits(0, trackBins[2]);
  fEventHist->SetVarTitle(0, trackAxisTitle[2]);
  
  fEventHist->SetBinLimits(1, trackBins[3]);
  fEventHist->SetVarTitle(1, trackAxisTitle[3]);
  
  SetStepNames(fEventHist);
  
  iTrackBin[2] = kNSpeciesBins;

  fTrackHistEfficiency = new AliCFContainer("fTrackHistEfficiency", "Tracking efficiency", 3, 3, iTrackBin);
  fTrackHistEfficiency->SetBinLimits(0, trackBins[0]);
  fTrackHistEfficiency->SetVarTitle(0, trackAxisTitle[0]);
  fTrackHistEfficiency->SetBinLimits(1, trackBins[1]);
  fTrackHistEfficiency->SetVarTitle(1, trackAxisTitle[1]);
  fTrackHistEfficiency->SetBinLimits(2, speciesBins);
  fTrackHistEfficiency->SetVarTitle(2, "particle species");
}

//_____________________________________________________________________________
AliUEHist::AliUEHist(const AliUEHist &c) :
  TObject(),
  fkRegions(4),
  fEventHist(0),
  fTrackHistEfficiency(0),
  fEtaMin(0),
  fEtaMax(0),
  fPtMin(0),
  fPtMax(0),
  fContaminationEnhancement(0),
  fCombineMinMax(0),
  fCache(0)
{
  //
  // AliUEHist copy constructor
  //

  ((AliUEHist &) c).Copy(*this);
}

//____________________________________________________________________
void AliUEHist::SetStepNames(AliCFContainer* container)
{
  // sets the names of the correction steps
  
  for (Int_t i=0; i<fgkCFSteps; i++)
    container->SetStepTitle(i, GetStepTitle((CFStep) i));
}

//____________________________________________________________________
AliUEHist::~AliUEHist()
{
  // Destructor
  
  for (Int_t i=0; i<fkRegions; i++)
  {
    if (fTrackHist[i])
    {
      delete fTrackHist[i];
      fTrackHist[i] = 0;
    }
  }
     
  if (fEventHist)
  {
    delete fEventHist;
    fEventHist = 0;
  }
  
  if (fTrackHistEfficiency)
  {
    delete fTrackHistEfficiency;
    fTrackHistEfficiency = 0;
  }

  if (fCache)
  {
    delete fCache;
    fCache = 0;
  }
}

//____________________________________________________________________
AliUEHist &AliUEHist::operator=(const AliUEHist &c)
{
  // assigment operator

  if (this != &c)
    ((AliUEHist &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliUEHist::Copy(TObject& c) const
{
  // copy function

  AliUEHist& target = (AliUEHist &) c;

  for (Int_t i=0; i<fkRegions; i++)
    if (fTrackHist[i])
      target.fTrackHist[i] = dynamic_cast<AliCFContainer*> (fTrackHist[i]->Clone());

  if (fEventHist)
    target.fEventHist = dynamic_cast<AliCFContainer*> (fEventHist->Clone());
  
  if (fTrackHistEfficiency)
    target.fTrackHistEfficiency = dynamic_cast<AliCFContainer*> (fTrackHistEfficiency->Clone());
}

//____________________________________________________________________
Long64_t AliUEHist::Merge(TCollection* list)
{
  // Merge a list of AliUEHist objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of objects
  const Int_t kMaxLists = fkRegions+2;
  TList** lists = new TList*[kMaxLists];
  
  for (Int_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliUEHist* entry = dynamic_cast<AliUEHist*> (obj);
    if (entry == 0) 
      continue;

    for (Int_t i=0; i<fkRegions; i++)
      if (entry->fTrackHist[i])
        lists[i]->Add(entry->fTrackHist[i]);
    
    lists[fkRegions]->Add(entry->fEventHist);
    lists[fkRegions+1]->Add(entry->fTrackHistEfficiency);

    count++;
  }
  for (Int_t i=0; i<fkRegions; i++)
    if (fTrackHist[i])
      fTrackHist[i]->Merge(lists[i]);
  
  fEventHist->Merge(lists[fkRegions]);
  fTrackHistEfficiency->Merge(lists[fkRegions+1]);

  for (Int_t i=0; i<kMaxLists; i++)
    delete lists[i];
    
  delete[] lists;

  return count+1;
}

//____________________________________________________________________
void AliUEHist::SetBinLimits(AliCFGridSparse* grid)
{
  // sets the bin limits in eta and pT defined by fEtaMin/Max, fPtMin/Max
  
  if (fEtaMax > fEtaMin)
    grid->SetRangeUser(0, fEtaMin, fEtaMax);
  if (fPtMax > fPtMin)
    grid->SetRangeUser(1, fPtMin, fPtMax);
}  

//____________________________________________________________________
void AliUEHist::ResetBinLimits(AliCFGridSparse* grid)
{
  // resets all bin limits 
  
  for (Int_t i=0; i<grid->GetNVar(); i++)
    if (grid->GetGrid()->GetAxis(i)->TestBit(TAxis::kAxisRange))
      grid->SetRangeUser(i, 0, -1);
}
  
//____________________________________________________________________
void AliUEHist::CountEmptyBins(AliUEHist::CFStep step, Float_t ptLeadMin, Float_t ptLeadMax)
{
  // prints the number of empty bins in the track end event histograms in the given step
  
  Int_t binBegin[4];
  Int_t binEnd[4];
  
  for (Int_t i=0; i<4; i++)
  {
    binBegin[i] = 1;
    binEnd[i]   = fTrackHist[0]->GetNBins(i);
  }
  
  if (fEtaMax > fEtaMin)
  {
    binBegin[0] = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(0)->FindBin(fEtaMin);
    binEnd[0]   = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(0)->FindBin(fEtaMax);
  }
  
  if (fPtMax > fPtMin)
  {
    binBegin[1] = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(fPtMin);
    binEnd[1]   = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(1)->FindBin(fPtMax);
  }
  
  if (ptLeadMax > ptLeadMin)
  {
    binBegin[2] = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(2)->FindBin(ptLeadMin);
    binEnd[2]   = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(2)->FindBin(ptLeadMax);
  }
  
  // start from multiplicity 1
  binBegin[3] = fTrackHist[0]->GetGrid(step)->GetGrid()->GetAxis(3)->FindBin(1);
  
  for (Int_t region=0; region<fkRegions; region++)
  {
    Int_t total = 0;
    Int_t count = 0;
    Int_t vars[4];
    
    for (Int_t i=0; i<4; i++)
      vars[i] = binBegin[i];
      
    AliCFGridSparse* grid = fTrackHist[region]->GetGrid(step);
    while (1)
    {
      if (grid->GetElement(vars) == 0)
      {
        Printf("Empty bin at eta=%.2f pt=%.2f pt_lead=%.2f mult=%.1f", 
          grid->GetBinCenter(0, vars[0]), 
          grid->GetBinCenter(1, vars[1]), 
          grid->GetBinCenter(2, vars[2]), 
          grid->GetBinCenter(3, vars[3])
        );
        count++;
      }
      
      vars[3]++;
      for (Int_t i=3; i>0; i--)
      {
        if (vars[i] == binEnd[i]+1)
        {
          vars[i] = binBegin[i];
          vars[i-1]++;
        }
      }
      
      if (vars[0] == binEnd[0]+1)
        break;
      total++;
    }
  
    Printf("Region %s has %d empty bins (out of %d bins)", GetRegionTitle((Region) region), count, total);
  }
}  

//____________________________________________________________________
TH1D* AliUEHist::GetUEHist(AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax)
{
  // Extracts the UE histogram at the given step and in the given region by projection and dividing tracks by events
  //
  // ptLeadMin, ptLeadMax: Only meaningful for vs. delta phi plot (third axis is ptLead)
  // Histogram has to be deleted by the caller of the function
  
  // unzoom all axes
  ResetBinLimits(fTrackHist[region]->GetGrid(step));
  ResetBinLimits(fEventHist->GetGrid(step));
  
  SetBinLimits(fTrackHist[region]->GetGrid(step));
    
  TH1D* tracks = 0;
  
  if (ptLeadMin < 0)
  {
    tracks = fTrackHist[region]->ShowProjection(2, step);
    tracks->GetYaxis()->SetTitle(fTrackHist[region]->GetTitle());
    if (fCombineMinMax && region == kMin)
    {
      ResetBinLimits(fTrackHist[kMax]->GetGrid(step));
      SetBinLimits(fTrackHist[kMax]->GetGrid(step));
      
      TH1D* tracks2 = fTrackHist[kMax]->ShowProjection(2, step);
      tracks->Add(tracks2);
      
      ResetBinLimits(fTrackHist[kMax]->GetGrid(step));
    }
      
    // normalize to get a density (deta dphi)
    TAxis* axis = fTrackHist[region]->GetGrid(step)->GetAxis(0);
    Float_t phiRegion = TMath::TwoPi() / 3;
    if (!fCombineMinMax && region == kMin)
      phiRegion /= 2;
    Float_t normalization = phiRegion * (axis->GetBinUpEdge(axis->GetLast()) - axis->GetBinLowEdge(axis->GetFirst()));
    //Printf("Normalization: %f", normalization);
    tracks->Scale(1.0 / normalization);
    
    TH1D* events = fEventHist->ShowProjection(0, step);
    tracks->Divide(events);
  }
  else
  {
    fTrackHist[region]->GetGrid(step)->SetRangeUser(2, ptLeadMin, ptLeadMax);
    tracks = fTrackHist[region]->GetGrid(step)->Project(4);
    Printf("Calculated histogram in %.2f <= pT <= %.2f --> %f entries", ptLeadMin, ptLeadMax, tracks->Integral());
    fTrackHist[region]->GetGrid(step)->SetRangeUser(2, 0, -1);
    
    // normalize to get a density (deta dphi)
    TAxis* axis = fTrackHist[region]->GetGrid(step)->GetAxis(0);
    Float_t normalization = fTrackHist[region]->GetGrid(step)->GetAxis(4)->GetBinWidth(1) * (axis->GetBinUpEdge(axis->GetLast()) - axis->GetBinLowEdge(axis->GetFirst()));
    //Printf("Normalization: %f", normalization);
    tracks->Scale(1.0 / normalization);
    
    TH1D* events = fEventHist->ShowProjection(0, step);
    Int_t nEvents = (Int_t) events->Integral(events->FindBin(ptLeadMin), events->FindBin(ptLeadMax));
    if (nEvents > 0)
      tracks->Scale(1.0 / nEvents);
  }

  ResetBinLimits(fTrackHist[region]->GetGrid(step));

  return tracks;
}

//____________________________________________________________________
void AliUEHist::CorrectTracks(CFStep step1, CFStep step2, TH1* trackCorrection, Int_t var1, Int_t var2)
{
  // corrects from step1 to step2 by multiplying the tracks with trackCorrection
  // trackCorrection can be a function of eta (var1 == 0), pT (var1 == 1), leading pT (var1 == 2), multiplicity (var1 == 3), delta phi (var1 == 4)
  // if var2 >= 0 a two dimension correction is assumed in trackCorrection
  //
  // if trackCorrection is 0, just copies content from step1 to step2
  
  for (Int_t region=0; region<fkRegions; region++)
    CorrectTracks(step1, step2, region, trackCorrection, var1, var2);
}

//____________________________________________________________________
void AliUEHist::CorrectTracks(CFStep step1, CFStep step2, Int_t region, TH1* trackCorrection, Int_t var1, Int_t var2)
{
  //
  // see documentation of CorrectTracks above
  //
  
  if (!fTrackHist[region])
    return;
   
  THnSparse* grid = fTrackHist[region]->GetGrid(step1)->GetGrid();
  THnSparse* target = fTrackHist[region]->GetGrid(step2)->GetGrid();
  
  // clear target histogram
  target->Reset();
  
  if (trackCorrection != 0)
  {
    if (grid->GetAxis(var1)->GetNbins() != trackCorrection->GetNbinsX())
      AliFatal(Form("Invalid binning (var1): %d %d", grid->GetAxis(var1)->GetNbins(), trackCorrection->GetNbinsX()));
      
    if (var2 >= 0 && grid->GetAxis(var2)->GetNbins() != trackCorrection->GetNbinsY())
      AliFatal(Form("Invalid binning (var2): %d %d", grid->GetAxis(var2)->GetNbins(), trackCorrection->GetNbinsY()));
  }
  
  // optimized implementation
  for (Int_t binIdx = 0; binIdx < grid->GetNbins(); binIdx++)
  {
    Int_t bins[5];
    Double_t value = grid->GetBinContent(binIdx, bins);
    Double_t error = grid->GetBinError(binIdx);
    
    if (trackCorrection != 0)
    {
      if (var2 < 0)
      {
        value *= trackCorrection->GetBinContent(bins[var1]);
        error *= trackCorrection->GetBinContent(bins[var1]);
      }
      else
      {
        value *= trackCorrection->GetBinContent(bins[var1], bins[var2]);
        error *= trackCorrection->GetBinContent(bins[var1], bins[var2]);
      }
    }
    
    target->SetBinContent(bins, value);
    target->SetBinError(bins, error);
  }
 
  Printf("AliUEHist::CorrectTracks: Corrected from %f to %f entries. Correction histogram: %f entries (integral: %f)", grid->GetEntries(), target->GetEntries(), (trackCorrection) ? trackCorrection->GetEntries() : -1.0, (trackCorrection) ? trackCorrection->Integral() : -1.0); 
}

//____________________________________________________________________
void AliUEHist::CorrectEvents(CFStep step1, CFStep step2, TH1D* eventCorrection, Int_t var)
{
  // corrects from step1 to step2 by multiplying the events with eventCorrection
  // eventCorrection is as function of leading pT (var == 0) or multiplicity (var == 1)
  //
  // if eventCorrection is 0, just copies content from step1 to step2
  
  AliCFGridSparse* grid = fEventHist->GetGrid(step1);
  AliCFGridSparse* target = fEventHist->GetGrid(step2);
  
  // clear target histogram
  target->GetGrid()->Reset();
  
  if (eventCorrection != 0 && grid->GetNBins(var) != eventCorrection->GetNbinsX())
    AliFatal(Form("Invalid binning: %d %d", grid->GetNBins(var), eventCorrection->GetNbinsX()));
  
  Int_t bins[2];
  for (Int_t x = 1; x <= grid->GetNBins(0); x++)
  {
    bins[0] = x;
    for (Int_t y = 1; y <= grid->GetNBins(1); y++)
    {
      bins[1] = y;
      
      Double_t value = grid->GetElement(bins);
      if (value != 0)
      {
        Double_t error = grid->GetElementError(bins);
        
        if (eventCorrection != 0)
        {
          value *= eventCorrection->GetBinContent(bins[var]);
          error *= eventCorrection->GetBinContent(bins[var]);
        }
        
        target->SetElement(bins, value);
        target->SetElementError(bins, error);
      }
    }
  }
  
  Printf("AliUEHist::CorrectEvents: Corrected from %f to %f entries. Correction histogram: %f entries (integral: %f)", grid->GetEntries(), target->GetEntries(), (eventCorrection) ? eventCorrection->GetEntries() : -1.0, (eventCorrection) ? eventCorrection->Integral() : -1.0); 
}

//____________________________________________________________________
void AliUEHist::Correct(AliUEHist* corrections)
{
  // applies the given corrections to extract from the step kCFStepReconstructed all previous steps
  //
  // in this object the data is expected in the step kCFStepReconstructed
  
  // ---- track level
  
  // bias due to migration in leading pT (because the leading particle is not reconstructed, and the subleading is used)
  // extracted as function of leading pT
  for (Int_t region = 0; region < fkRegions; region++)
  {
    if (!fTrackHist[region])
      continue;

    const char* projAxis = "z";
    Int_t secondBin = -1;

    if (fTrackHist[region]->GetNVar() == 5)
    {
      projAxis = "yz";
      secondBin = 4;
    }
    
    #if 0
      TH1* leadingBias = (TH1*) corrections->GetBias(kCFStepReconstructed, kCFStepTracked, region, projAxis); // from MC
      Printf("WARNING: Using MC bias correction");
    #else
      TH1* leadingBias = (TH1*) GetBias(kCFStepBiasStudy, kCFStepReconstructed, region, projAxis);          // from data
    #endif
    CorrectTracks(kCFStepReconstructed, kCFStepTracked, region, leadingBias, 2, secondBin);
    if (region == kMin && fCombineMinMax)
    {
      CorrectTracks(kCFStepReconstructed, kCFStepTracked, kMax, leadingBias, 2, secondBin);
      delete leadingBias;
      break;
    }
    delete leadingBias;
  }
  CorrectEvents(kCFStepReconstructed, kCFStepTracked, 0, 0);
  
  // correct with kCFStepTracked --> kCFStepTrackedOnlyPrim
  TH2D* contamination = corrections->GetTrackingContamination();
  if (corrections->fContaminationEnhancement)
  {
    Printf("Applying contamination enhancement");
    
    for (Int_t x=1; x<=contamination->GetNbinsX(); x++)
      for (Int_t y=1; y<=contamination->GetNbinsX(); y++)
        contamination->SetBinContent(x, y, contamination->GetBinContent(x, y) * corrections->fContaminationEnhancement->GetBinContent(corrections->fContaminationEnhancement->GetXaxis()->FindBin(contamination->GetYaxis()->GetBinCenter(y))));
  }
  CorrectTracks(kCFStepTracked, kCFStepTrackedOnlyPrim, contamination, 0, 1);
  CorrectEvents(kCFStepTracked, kCFStepTrackedOnlyPrim, 0, 0);
  delete contamination;
  
  // correct with kCFStepTrackedOnlyPrim --> kCFStepAnaTopology
  TH2D* efficiencyCorrection = corrections->GetTrackingEfficiencyCorrection();
  CorrectTracks(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, efficiencyCorrection, 0, 1);
  CorrectEvents(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, 0, 0);
  delete efficiencyCorrection;
  
  // copy 
  CorrectTracks(kCFStepAnaTopology, kCFStepVertex, 0, -1);
  CorrectEvents(kCFStepAnaTopology, kCFStepVertex, 0, 0);
  
  // vertex correction on the level of events as function of multiplicity, weighting tracks and events with the same factor
  // practically independent of low pT cut 
  TH1D* vertexCorrection = (TH1D*) corrections->GetEventEfficiency(kCFStepVertex, kCFStepTriggered, 1);
  
  // convert stage from true multiplicity to observed multiplicity by simple conversion factor
  TH1D* vertexCorrectionObs = (TH1D*) vertexCorrection->Clone("vertexCorrection2");
  vertexCorrectionObs->Reset();
  
  TF1* func = new TF1("func", "[1]+[0]/(x-[2])");
  vertexCorrection->Fit(func, "0", "", 0, 3);

  for (Int_t i=1; i<=vertexCorrectionObs->GetNbinsX(); i++)
  {
    Float_t xPos = 1.0 / 0.77 * vertexCorrectionObs->GetXaxis()->GetBinCenter(i);
    if (xPos < 4)
      vertexCorrectionObs->SetBinContent(i, func->Eval(xPos));
    else
      vertexCorrectionObs->SetBinContent(i, vertexCorrection->Interpolate(xPos));
  }
 
  #if 0
  new TCanvas;
  vertexCorrection->DrawCopy();
  vertexCorrectionObs->SetLineColor(2);
  vertexCorrectionObs->DrawCopy("same");
  func->SetRange(0, 4);
  func->DrawClone("same");
  #endif
  
  CorrectTracks(kCFStepVertex, kCFStepTriggered, vertexCorrectionObs, 3);
  CorrectEvents(kCFStepVertex, kCFStepTriggered, vertexCorrectionObs, 1);
  delete vertexCorrectionObs;
  delete vertexCorrection;
  delete func;
  
  // copy 
  CorrectTracks(kCFStepTriggered, kCFStepAll, 0, -1);
  CorrectEvents(kCFStepTriggered, kCFStepAll, 0, 0);
}

//____________________________________________________________________
TH1* AliUEHist::GetTrackEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2, Int_t source)
{
  // creates a track-level efficiency by dividing step2 by step1
  // projected to axis1 and axis2 (optional if >= 0)
  //
  // source: 0 = fTrackHist; 1 = fTrackHistEfficiency
  
  // integrate over regions
  // cache it for efficiency (usually more than one efficiency is requested)
  
  AliCFContainer* sourceContainer = 0;
  
  if (source == 0)
  {
    if (!fCache)
    {
      fCache = (AliCFContainer*) fTrackHist[0]->Clone();
      for (Int_t i = 1; i < fkRegions; i++)
        if (fTrackHist[i])
          fCache->Add(fTrackHist[i]);
    }
    sourceContainer = fCache;
  }
  else if (source == 1)
  {
    sourceContainer = fTrackHistEfficiency;
    // step offset because we start with kCFStepAnaTopology
    step1 = (CFStep) ((Int_t) step1 - (Int_t) kCFStepAnaTopology);
    step2 = (CFStep) ((Int_t) step2 - (Int_t) kCFStepAnaTopology);
  }
  else
    return 0;
        
  // reset all limits and set the right ones except those in axis1 and axis2
  ResetBinLimits(sourceContainer->GetGrid(step1));
  ResetBinLimits(sourceContainer->GetGrid(step2));
  if (fEtaMax > fEtaMin && axis1 != 0 && axis2 != 0)
  {
    sourceContainer->GetGrid(step1)->SetRangeUser(0, fEtaMin, fEtaMax);
    sourceContainer->GetGrid(step2)->SetRangeUser(0, fEtaMin, fEtaMax);
  }
  if (fPtMax > fPtMin && axis1 != 1 && axis2 != 1)
  {
    sourceContainer->GetGrid(step1)->SetRangeUser(1, fPtMin, fPtMax);
    sourceContainer->GetGrid(step2)->SetRangeUser(1, fPtMin, fPtMax);
  }
  
  TH1* measured = 0;
  TH1* generated = 0;
    
  if (axis2 >= 0)
  {
    generated = sourceContainer->Project(axis1, axis2, step1);
    measured = sourceContainer->Project(axis1, axis2, step2);
  }
  else
  {
    generated = sourceContainer->Project(axis1, step1);
    measured = sourceContainer->Project(axis1, step2);
  }
  
  // check for bins with less than 100 entries, print warning
  Int_t binBegin[2];
  Int_t binEnd[2];
  
  binBegin[0] = 1;
  binBegin[1] = 1;
  
  binEnd[0] = generated->GetNbinsX();
  binEnd[1] = generated->GetNbinsY();
  
  if (fEtaMax > fEtaMin)
  {
    if (axis1 == 0)
    {
      binBegin[0] = generated->GetXaxis()->FindBin(fEtaMin);
      binEnd[0]   = generated->GetXaxis()->FindBin(fEtaMax);
    }
    if (axis2 == 0)
    {
      binBegin[1] = generated->GetYaxis()->FindBin(fEtaMin);
      binEnd[1]   = generated->GetYaxis()->FindBin(fEtaMax);
    }
  }
  
  if (fPtMax > fPtMin)
  {
    // TODO this is just checking up to 15 for now
    Float_t ptMax = TMath::Min((Float_t) 15., fPtMax);
    if (axis1 == 1)
    {
      binBegin[0] = generated->GetXaxis()->FindBin(fPtMin);
      binEnd[0]   = generated->GetXaxis()->FindBin(ptMax);
    }
    if (axis2 == 1)
    {
      binBegin[1] = generated->GetYaxis()->FindBin(fPtMin);
      binEnd[1]   = generated->GetYaxis()->FindBin(ptMax);
    }
  }
  
  Int_t total = 0;
  Int_t count = 0;
  Int_t vars[2];
  
  for (Int_t i=0; i<2; i++)
    vars[i] = binBegin[i];
    
  const Int_t limit = 50;
  while (1)
  {
    if (generated->GetDimension() == 1 && generated->GetBinContent(vars[0]) < limit)
    {
      Printf("Empty bin at %s=%.2f (%.2f entries)", generated->GetXaxis()->GetTitle(), generated->GetXaxis()->GetBinCenter(vars[0]), generated->GetBinContent(vars[0]));
      count++;
    } 
    else if (generated->GetDimension() == 2 && generated->GetBinContent(vars[0], vars[1]) < limit)
    {
      Printf("Empty bin at %s=%.2f %s=%.2f (%.2f entries)", 
        generated->GetXaxis()->GetTitle(), generated->GetXaxis()->GetBinCenter(vars[0]),
        generated->GetYaxis()->GetTitle(), generated->GetYaxis()->GetBinCenter(vars[1]),
        generated->GetBinContent(vars[0], vars[1]));
      count++;
    }
    
    vars[1]++;
    if (vars[1] == binEnd[1]+1)
    {
      vars[1] = binBegin[1];
      vars[0]++;
    }
    
    if (vars[0] == binEnd[0]+1)
      break;
    total++;
  }

  Printf("Correction has %d empty bins (out of %d bins)", count, total);
  
  measured->Divide(measured, generated, 1, 1, "B");
  
  delete generated;
  
  ResetBinLimits(sourceContainer->GetGrid(step1));
  ResetBinLimits(sourceContainer->GetGrid(step2));
  
  return measured;
}

//____________________________________________________________________
TH1* AliUEHist::GetEventEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2, Float_t ptLeadMin, Float_t ptLeadMax)
{
  // creates a event-level efficiency by dividing step2 by step1
  // projected to axis1 and axis2 (optional if >= 0)
  
  if (ptLeadMax > ptLeadMin)
  {
    fEventHist->GetGrid(step1)->SetRangeUser(0, ptLeadMin, ptLeadMax);
    fEventHist->GetGrid(step2)->SetRangeUser(0, ptLeadMin, ptLeadMax);
  }
  
  TH1* measured = 0;
  TH1* generated = 0;
    
  if (axis2 >= 0)
  {
    generated = fEventHist->Project(axis1, axis2, step1);
    measured = fEventHist->Project(axis1, axis2, step2);
  }
  else
  {
    generated = fEventHist->Project(axis1, step1);
    measured = fEventHist->Project(axis1, step2);
  }
  
  measured->Divide(measured, generated, 1, 1, "B");
  
  delete generated;
  
  if (ptLeadMax > ptLeadMin)
  {
    fEventHist->GetGrid(step1)->SetRangeUser(0, 0, -1);
    fEventHist->GetGrid(step2)->SetRangeUser(0, 0, -1);
  }
  
  return measured;
}

//____________________________________________________________________
void AliUEHist::WeightHistogram(TH3* hist1, TH1* hist2)
{
  // weights each entry of the 3d histogram hist1 with the 1d histogram hist2 
  // where the matching is done of the z axis of hist1 with the x axis of hist2
  
  if (hist1->GetNbinsZ() != hist2->GetNbinsX())
    AliFatal(Form("Inconsistent binning %d %d", hist1->GetNbinsZ(), hist2->GetNbinsX()));
  
  for (Int_t x=1; x<=hist1->GetNbinsX(); x++)
  {
    for (Int_t y=1; y<=hist1->GetNbinsY(); y++)
    {
      for (Int_t z=1; z<=hist1->GetNbinsZ(); z++)
      {
        if (hist2->GetBinContent(z) > 0)
        {
          hist1->SetBinContent(x, y, z, hist1->GetBinContent(x, y, z) / hist2->GetBinContent(z));
          hist1->SetBinError(x, y, z, hist1->GetBinError(x, y, z) / hist2->GetBinContent(z));
        }
        else
        {
          hist1->SetBinContent(x, y, z, 0);
          hist1->SetBinError(x, y, z, 0);
        }
      }
    }
  }
}  

//____________________________________________________________________
TH1* AliUEHist::GetBias(CFStep step1, CFStep step2, Int_t region, const char* axis, Float_t leadPtMin, Float_t leadPtMax)
{
  // extracts the track-level bias (integrating out the multiplicity) between two steps (dividing step2 by step1)
  // in the given region (sum over all regions is calculated if region == -1)
  // done by weighting the track-level distribution with the number of events as function of leading pT
  // and then calculating the ratio between the distributions
  // projected to axis which is a TH3::Project3D string, e.g. "x", or "yx"
  //   no projection is done if axis == 0
  
  AliCFContainer* tmp = 0;
  
  if (region == -1)
  {
    tmp = (AliCFContainer*) fTrackHist[0]->Clone();
    for (Int_t i = 1; i < fkRegions; i++)
      if (fTrackHist[i])
	tmp->Add(fTrackHist[i]);
  }
  else if (region == kMin && fCombineMinMax)
  {
    tmp = (AliCFContainer*) fTrackHist[kMin]->Clone();
    tmp->Add(fTrackHist[kMax]);
  }
  else
    tmp = fTrackHist[region];
  
  ResetBinLimits(tmp->GetGrid(step1));
  ResetBinLimits(fEventHist->GetGrid(step1));
  SetBinLimits(tmp->GetGrid(step1));
  
  ResetBinLimits(tmp->GetGrid(step2));
  ResetBinLimits(fEventHist->GetGrid(step2));
  SetBinLimits(tmp->GetGrid(step2));
  
  TH1D* events1 = fEventHist->Project(0, step1);
  TH3D* hist1 = tmp->Project(0, tmp->GetNVar()-1, 2, step1);
  WeightHistogram(hist1, events1);
  
  TH1D* events2 = fEventHist->Project(0, step2);
  TH3D* hist2 = tmp->Project(0, tmp->GetNVar()-1, 2, step2);
  WeightHistogram(hist2, events2);
  
  TH1* generated = hist1;
  TH1* measured = hist2;
  
  if (axis)
  {
    if (leadPtMax > leadPtMin)
    {
      hist1->GetZaxis()->SetRangeUser(leadPtMin, leadPtMax);
      hist2->GetZaxis()->SetRangeUser(leadPtMin, leadPtMax);
    }
    
    if (fEtaMax > fEtaMin && !TString(axis).Contains("x"))
    {
      hist1->GetXaxis()->SetRangeUser(fEtaMin, fEtaMax);
      hist2->GetXaxis()->SetRangeUser(fEtaMin, fEtaMax);
    }
   
    generated = hist1->Project3D(axis);
    measured  = hist2->Project3D(axis);
    
    // delete hists here if projection has been done
    delete hist1;
    delete hist2;
  }
  
  measured->Divide(generated);
  
  delete events1;
  delete events2;
  delete generated;
  
  ResetBinLimits(tmp->GetGrid(step1));
  ResetBinLimits(tmp->GetGrid(step2));
  
  if ((region == -1) || (region == kMin && fCombineMinMax))
    delete tmp;
  
  return measured;
}

//____________________________________________________________________
TH2D* AliUEHist::GetTrackingEfficiency()
{
  // extracts the tracking efficiency by calculating the efficiency from step kCFStepAnaTopology to kCFStepTrackedOnlyPrim
  // integrates over the regions and all other variables than pT and eta to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepAnaTopology, kCFStepTrackedOnlyPrim, 0, 1));
}
  
//____________________________________________________________________
TH1D* AliUEHist::GetTrackingEfficiency(Int_t axis)
{
  // extracts the tracking efficiency by calculating the efficiency from step kCFStepAnaTopology to kCFStepTrackedOnlyPrim
  // integrates over the regions and all other variables than pT (axis == 0) and eta (axis == 1) to increase the statistics

  return dynamic_cast<TH1D*> (GetTrackEfficiency(kCFStepAnaTopology, kCFStepTrackedOnlyPrim, axis));
}

//____________________________________________________________________
TH2D* AliUEHist::GetTrackingCorrection()
{
  // extracts the tracking correction by calculating the efficiency from step kCFStepAnaTopology to kCFStepTracked
  // integrates over the regions and all other variables than pT and eta to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepTracked, kCFStepAnaTopology, 0, 1));
}
  
//____________________________________________________________________
TH1D* AliUEHist::GetTrackingCorrection(Int_t axis)
{
  // extracts the tracking correction by calculating the efficiency from step kCFStepAnaTopology to kCFStepTracked
  // integrates over the regions and all other variables than pT (axis == 0) and eta (axis == 1) to increase the statistics

  return dynamic_cast<TH1D*> (GetTrackEfficiency(kCFStepTracked, kCFStepAnaTopology, axis));
}

//____________________________________________________________________
TH2D* AliUEHist::GetTrackingEfficiencyCorrection()
{
  // extracts the tracking correction by calculating the efficiency from step kCFStepAnaTopology to kCFStepTracked
  // integrates over the regions and all other variables than pT and eta to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, 0, 1));
}
  
//____________________________________________________________________
TH1D* AliUEHist::GetTrackingEfficiencyCorrection(Int_t axis)
{
  // extracts the tracking correction by calculating the efficiency from step kCFStepAnaTopology to kCFStepTracked
  // integrates over the regions and all other variables than pT (axis == 0) and eta (axis == 1) to increase the statistics

  return dynamic_cast<TH1D*> (GetTrackEfficiency(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, axis));
}

//____________________________________________________________________
TH2D* AliUEHist::GetTrackingContamination()
{
  // extracts the tracking contamination by secondaries by calculating the efficiency from step kCFStepTrackedOnlyPrim to kCFStepTracked
  // integrates over the regions and all other variables than pT and eta to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepTracked, kCFStepTrackedOnlyPrim, 0, 1));
}
  
//____________________________________________________________________
TH1D* AliUEHist::GetTrackingContamination(Int_t axis)
{
  // extracts the tracking contamination by secondaries by calculating the efficiency from step kCFStepTrackedOnlyPrim to kCFStepTracked
  // integrates over the regions and all other variables than pT (axis == 0) and eta (axis == 1) to increase the statistics

  return dynamic_cast<TH1D*> (GetTrackEfficiency(kCFStepTracked, kCFStepTrackedOnlyPrim, axis));
}

//____________________________________________________________________
const char* AliUEHist::GetRegionTitle(Region region)
{
  // returns the name of the given region
  
  switch (region)
  {
    case kToward:
      return "Towards";
    case kAway:
      return "Away";
    case kMin:
      return (fCombineMinMax) ? "Transverse" : "Min";
    case kMax:
      return "Max";
  }
  
  return 0;
}
  
//____________________________________________________________________
const char* AliUEHist::GetStepTitle(CFStep step)
{
  // returns the name of the given step
  
  switch (step)
  {
    case kCFStepAll:
      return "All events";
    case kCFStepTriggered:
      return "Triggered";
    case kCFStepVertex:
      return "Primary Vertex";
    case kCFStepAnaTopology:
      return "Required analysis topology";
    case kCFStepTrackedOnlyPrim:
      return "Tracked (matched MC, only primaries)";
    case kCFStepTracked:
      return "Tracked (matched MC, all)";
    case kCFStepReconstructed:
      return "Reconstructed";
    case kCFStepRealLeading:
      return "Correct leading particle identified";
    case kCFStepBiasStudy:
      return "Bias study applying tracking efficiency";
    case kCFStepBiasStudy2:
      return "Bias study applying tracking efficiency in two steps";
  }
  
  return 0;
}

//____________________________________________________________________
void AliUEHist::CopyReconstructedData(AliUEHist* from)
{
  // copies those histograms extracted from ESD to this object
  
  // TODO at present only the pointers are copied
  
  for (Int_t region=0; region<4; region++)
  {
    if (!fTrackHist[region])
      continue;
  
    fTrackHist[region]->SetGrid(AliUEHist::kCFStepReconstructed, from->fTrackHist[region]->GetGrid(AliUEHist::kCFStepReconstructed));
    fTrackHist[region]->SetGrid(AliUEHist::kCFStepBiasStudy,     from->fTrackHist[region]->GetGrid(AliUEHist::kCFStepBiasStudy));
  }
    
  fEventHist->SetGrid(AliUEHist::kCFStepReconstructed, from->fEventHist->GetGrid(AliUEHist::kCFStepReconstructed));
  fEventHist->SetGrid(AliUEHist::kCFStepBiasStudy,     from->fEventHist->GetGrid(AliUEHist::kCFStepBiasStudy));
}
