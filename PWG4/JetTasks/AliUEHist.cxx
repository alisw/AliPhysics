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
#include "AliTHn.h"

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
  fCentralityMin(0),
  fCentralityMax(0),
  fZVtxMin(0),
  fZVtxMax(0),
  fContaminationEnhancement(0),
  fCombineMinMax(0),
  fCache(0),
  fHistogramType(reqHist)
{
  // Constructor
    
  for (UInt_t i=0; i<fkRegions; i++)
    fTrackHist[i] = 0;
    
  if (strlen(reqHist) == 0)
    return;
    
  AliLog::SetClassDebugLevel("AliCFContainer", -1);
  AliLog::SetClassDebugLevel("AliCFGridSparse", -3);
    
  const char* title = "";
    
  // track level
  Int_t nTrackVars = 4; // eta vs pT vs pT,lead (vs delta phi) vs multiplicity
  Int_t iTrackBin[6];
  Double_t* trackBins[6];
  const char* trackAxisTitle[6];
  
  // eta
  iTrackBin[0] = 20;
  Double_t etaBins[20+1];
  for (Int_t i=0; i<=iTrackBin[0]; i++)
    etaBins[i] = -1.0 + 0.1 * i;
  trackBins[0] = etaBins;
  trackAxisTitle[0] = "#eta";
  
  // delta eta
  const Int_t kNDeltaEtaBins = 40;
  Double_t deltaEtaBins[kNDeltaEtaBins+1];
  for (Int_t i=0; i<=kNDeltaEtaBins; i++)
    deltaEtaBins[i] = -2.0 + 0.1 * i;
  
  // pT
  iTrackBin[1] = 22;
  Double_t pTBins[] = {0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0};
  trackBins[1] = pTBins;
  trackAxisTitle[1] = "p_{T} (GeV/c)";
  
  // pT,lead binning 1
  const Int_t kNLeadingpTBins = 100;
  Double_t leadingpTBins[kNLeadingpTBins+1];
  for (Int_t i=0; i<=kNLeadingpTBins; i++)
    leadingpTBins[i] = 0.5 * i;
  
  // pT,lead binning 2
  const Int_t kNLeadingpTBins2 = 9;
  Double_t leadingpTBins2[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0 };
  
  // phi,lead; this binning starts at -pi/2 and is modulo 3
  const Int_t kNLeadingPhiBins = 72;
  Double_t leadingPhiBins[kNLeadingPhiBins+1];
  for (Int_t i=0; i<=kNLeadingPhiBins; i++)
    leadingPhiBins[i] = -TMath::Pi() / 2 + 1.0 / kNLeadingPhiBins * i * TMath::TwoPi();
    
  // multiplicity
  const Int_t kNMultiplicityBins = 15;
  Double_t multiplicityBins[kNMultiplicityBins+1];
  for (Int_t i=0; i<=kNMultiplicityBins; i++)
    multiplicityBins[i] = -0.5 + i;
  multiplicityBins[kNMultiplicityBins] = 200;
  
  trackBins[3] = multiplicityBins;
  iTrackBin[3] = kNMultiplicityBins;
  trackAxisTitle[3] = "multiplicity";
  
  // centrality (in %)
  const Int_t kNCentralityBins = 5 + 1 + 9;
  Double_t centralityBins[] = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 };
  
  // particle species
  const Int_t kNSpeciesBins = 4; // pi, K, p, rest
  Double_t speciesBins[] = { -0.5, 0.5, 1.5, 2.5, 3.5 };
  
  // vtx-z axis
  const Int_t kNVertexBins = 7;
  Double_t vertexBins[] = { -7, -5, -3, -1, 1, 3, 5, 7 };
  
  Bool_t useVtxAxis = kFALSE;
  
  // selection depending on requested histogram
  Int_t axis = -1; // 0 = pT,lead, 1 = phi,lead
  if (strcmp(reqHist, "NumberDensitypT") == 0)
  {
    axis = 0;
    title = "d^{2}N_{ch}/d#varphid#eta";
  }
  else if (strcmp(reqHist, "NumberDensityPhi") == 0)
  {
    axis = 1;
    title = "d^{2}N_{ch}/d#varphid#eta";
  }
  else if (strcmp(reqHist, "NumberDensityPhiCentrality") == 0 || strcmp(reqHist, "NumberDensityPhiCentralityVtx") == 0)
  {
    if (strcmp(reqHist, "NumberDensityPhiCentralityVtx") == 0)
    {
      reqHist = "NumberDensityPhiCentrality";
      fHistogramType = reqHist;
      useVtxAxis = kTRUE;
    }
    axis = 2;
    title = "d^{2}N_{ch}/d#varphid#eta";
  }
  else if (strcmp(reqHist, "SumpT") == 0)
  {
    axis = 0;
    title = "d^{2}#Sigma p_{T}/d#varphid#eta";
  }
  else
    AliFatal(Form("Invalid histogram requested: %s", reqHist));
  
  UInt_t initRegions = fkRegions;
  
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
    trackAxisTitle[4] = "#Delta #varphi w.r.t. leading track";
  }
  else if (axis == 2)
  {
    nTrackVars = 5;
    initRegions = 1;
  
    iTrackBin[0] = kNDeltaEtaBins;
    trackBins[0] = deltaEtaBins;
    trackAxisTitle[0] = "#Delta#eta";
  
    iTrackBin[2] = kNLeadingpTBins2;
    trackBins[2] = leadingpTBins2;
    trackAxisTitle[2] = "leading p_{T} (GeV/c)";
    
    trackBins[3] = centralityBins;
    iTrackBin[3] = kNCentralityBins;
    trackAxisTitle[3] = "centrality";
  
    iTrackBin[4] = kNLeadingPhiBins;
    trackBins[4] = leadingPhiBins;
    trackAxisTitle[4] = "#Delta#varphi (rad.)";

    if (useVtxAxis)
    {
      nTrackVars = 6;
      iTrackBin[5] = kNVertexBins;
      trackBins[5] = vertexBins;
      trackAxisTitle[5] = "z-vtx (cm)";
    }
  }
    
  for (UInt_t i=0; i<initRegions; i++)
  {
    if (strcmp(reqHist, "NumberDensityPhiCentrality") == 0)
      fTrackHist[i] = new AliTHn(Form("fTrackHist_%d", i), title, fgkCFSteps, nTrackVars, iTrackBin);
    else
      fTrackHist[i] = new AliCFContainer(Form("fTrackHist_%d", i), title, fgkCFSteps, nTrackVars, iTrackBin);
    
    for (Int_t j=0; j<nTrackVars; j++)
    {
      fTrackHist[i]->SetBinLimits(j, trackBins[j]);
      fTrackHist[i]->SetVarTitle(j, trackAxisTitle[j]);
    }
    
    SetStepNames(fTrackHist[i]);
  }
  
  // event level
  Int_t nEventVars = 2;
  Int_t iEventBin[3];

  // track 3rd and 4th axis --> event 1st and 2nd axis
  iEventBin[0] = iTrackBin[2];
  iEventBin[1] = iTrackBin[3];
  
  // plus track 5th axis (in certain cases)
  if (axis == 2 && useVtxAxis)
  {
    nEventVars = 3;
    iEventBin[2] = iTrackBin[5];
  }
  
  fEventHist = new AliCFContainer("fEventHist", title, fgkCFSteps, nEventVars, iEventBin);
  
  fEventHist->SetBinLimits(0, trackBins[2]);
  fEventHist->SetVarTitle(0, trackAxisTitle[2]);
  
  fEventHist->SetBinLimits(1, trackBins[3]);
  fEventHist->SetVarTitle(1, trackAxisTitle[3]);
  
  if (axis == 2 && useVtxAxis)
  {
    fEventHist->SetBinLimits(2, trackBins[5]);
    fEventHist->SetVarTitle(2, trackAxisTitle[5]);
  }

  SetStepNames(fEventHist);
  
  iTrackBin[2] = kNSpeciesBins;

  fTrackHistEfficiency = new AliCFContainer("fTrackHistEfficiency", "Tracking efficiency", 3, 4, iTrackBin);
  fTrackHistEfficiency->SetBinLimits(0, trackBins[0]);
  fTrackHistEfficiency->SetVarTitle(0, trackAxisTitle[0]);
  fTrackHistEfficiency->SetBinLimits(1, trackBins[1]);
  fTrackHistEfficiency->SetVarTitle(1, trackAxisTitle[1]);
  fTrackHistEfficiency->SetBinLimits(2, speciesBins);
  fTrackHistEfficiency->SetVarTitle(2, "particle species");
  fTrackHistEfficiency->SetBinLimits(3, trackBins[3]);
  fTrackHistEfficiency->SetVarTitle(3, trackAxisTitle[3]);
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
  fCentralityMin(0),
  fCentralityMax(0),
  fZVtxMin(0),
  fZVtxMax(0),
  fContaminationEnhancement(0),
  fCombineMinMax(0),
  fCache(0),
  fHistogramType()
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
  
  for (UInt_t i=0; i<fkRegions; i++)
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

  for (UInt_t i=0; i<fkRegions; i++)
    if (fTrackHist[i])
      target.fTrackHist[i] = dynamic_cast<AliCFContainer*> (fTrackHist[i]->Clone());

  if (fEventHist)
    target.fEventHist = dynamic_cast<AliCFContainer*> (fEventHist->Clone());
  
  if (fTrackHistEfficiency)
    target.fTrackHistEfficiency = dynamic_cast<AliCFContainer*> (fTrackHistEfficiency->Clone());
    
  target.fEtaMin = fEtaMin;
  target.fEtaMax = fEtaMax;
  target.fPtMin = fPtMin;
  target.fPtMax = fPtMax;
  target.fCentralityMin = fCentralityMin;
  target.fCentralityMax = fCentralityMax;
  target.fZVtxMin = fZVtxMin;
  target.fZVtxMax = fZVtxMax;
  
  if (fContaminationEnhancement)
    target.fContaminationEnhancement = dynamic_cast<TH1F*> (fContaminationEnhancement->Clone());
    
  target.fCombineMinMax = fCombineMinMax;
  target.fHistogramType = fHistogramType;
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
  const UInt_t kMaxLists = fkRegions+2;
  TList** lists = new TList*[kMaxLists];
  
  for (UInt_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliUEHist* entry = dynamic_cast<AliUEHist*> (obj);
    if (entry == 0) 
      continue;

    for (UInt_t i=0; i<fkRegions; i++)
      if (entry->fTrackHist[i])
        lists[i]->Add(entry->fTrackHist[i]);
    
    lists[fkRegions]->Add(entry->fEventHist);
    lists[fkRegions+1]->Add(entry->fTrackHistEfficiency);

    count++;
  }
  for (UInt_t i=0; i<fkRegions; i++)
    if (fTrackHist[i])
      fTrackHist[i]->Merge(lists[i]);
  
  fEventHist->Merge(lists[fkRegions]);
  fTrackHistEfficiency->Merge(lists[fkRegions+1]);

  for (UInt_t i=0; i<kMaxLists; i++)
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
  
  for (UInt_t region=0; region<fkRegions; region++)
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
TH1* AliUEHist::GetUEHist(AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Int_t twoD, Bool_t etaNorm, Int_t* normEvents)
{
  // Extracts the UE histogram at the given step and in the given region by projection and dividing tracks by events
  //
  // ptLeadMin, ptLeadMax: Only meaningful for vs. delta phi plot (third axis is ptLead)
  // Histogram has to be deleted by the caller of the function
  //
  // twoD: 0: 1D histogram as function of phi
  //       1: 2D histogram, deltaphi, deltaeta
  //       10: 1D histogram, within |deltaeta| < 0.8
  //       11: 1D histogram, within 0.8 < |deltaeta| < 1.6
  //
  // etaNorm: if kTRUE (default), the distributions are divided by the area in delta eta
  //
  // normEvents: if non-0 the number of events/trigger particles for the normalization is filled
  
  // unzoom all axes
  ResetBinLimits(fTrackHist[region]->GetGrid(step));
  ResetBinLimits(fEventHist->GetGrid(step));
  
  SetBinLimits(fTrackHist[region]->GetGrid(step));
  
  // z vtx
  if (fZVtxMax > fZVtxMin)
  {
    Printf("Using z-vtx range %f --> %f", fZVtxMin, fZVtxMax);
    fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(5)->SetRangeUser(fZVtxMin, fZVtxMax);
    fEventHist->GetGrid(step)->GetGrid()->GetAxis(2)->SetRangeUser(fZVtxMin, fZVtxMax);
  }
    
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
    Float_t normalization = phiRegion;
    if (etaNorm)
      normalization *= (axis->GetBinUpEdge(axis->GetLast()) - axis->GetBinLowEdge(axis->GetFirst()));
    //Printf("Normalization: %f", normalization);
    tracks->Scale(1.0 / normalization);
    
    TH1D* events = fEventHist->ShowProjection(0, step);
    tracks->Divide(events);
  }
  else
  {
    if (multBinEnd >= multBinBegin)
    {
      Printf("Using multiplicity range %d --> %d", multBinBegin, multBinEnd);
      fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(3)->SetRange(multBinBegin, multBinEnd);
      fEventHist->GetGrid(step)->GetGrid()->GetAxis(1)->SetRange(multBinBegin, multBinEnd);
    }
    
    Int_t firstBin = fTrackHist[region]->GetAxis(2, step)->FindBin(ptLeadMin);
    Int_t lastBin = fTrackHist[region]->GetAxis(2, step)->FindBin(ptLeadMax);
    
    Printf("Using leading pT range %d --> %d", firstBin, lastBin);
    
    fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(2)->SetRange(firstBin, lastBin);
    
    if (twoD == 0)
      tracks = (TH1D*) fTrackHist[region]->GetGrid(step)->Project(4);
    else
      tracks = (TH1D*) fTrackHist[region]->GetGrid(step)->Project(4, 0);
      
    Printf("Calculated histogram --> %f tracks", tracks->Integral());
    fTrackHist[region]->GetGrid(step)->SetRangeUser(2, 0, -1);
    
    if (twoD == 10 || twoD == 11)
    {
      Float_t etaLimit = 1.0;
      if (twoD == 10)
      {
        tracks = (TH1D*) ((TH2*) tracks)->ProjectionX("proj", tracks->GetYaxis()->FindBin(-etaLimit + 0.01), tracks->GetYaxis()->FindBin(etaLimit - 0.01))->Clone();
        
        // TODO calculate acc with 2 * (deta - 0.5 * deta*deta / 1.6)
        tracks->Scale(1. / 0.75);
      }
      else if (twoD == 11)
      {
        TH1D* tracksTmp1 = (TH1D*) ((TH2*) tracks)->ProjectionX("proj1", tracks->GetYaxis()->FindBin(-etaLimit * 2 + 0.01), tracks->GetYaxis()->FindBin(-etaLimit - 0.01))->Clone();
        TH1D* tracksTmp2 = ((TH2*) tracks)->ProjectionX("proj2", tracks->GetYaxis()->FindBin(etaLimit + 0.01), tracks->GetYaxis()->FindBin(2 * etaLimit - 0.01));
        
        tracksTmp1->Add(tracksTmp2);
        
        delete tracks;
        tracks = tracksTmp1;
        delete tracksTmp2;
        
        tracks->Scale(1. / 0.25);
      }
    }
    
    // normalize to get a density (deta dphi)
    // TODO normalization may be off for 2d histogram
    Float_t normalization = fTrackHist[region]->GetGrid(step)->GetAxis(4)->GetBinWidth(1);
    
    if (etaNorm)
    {
      TAxis* axis = fTrackHist[region]->GetGrid(step)->GetAxis(0);
      if (strcmp(axis->GetTitle(), "#eta") == 0)
      {
	Printf("Normalizing using eta-axis range");
	normalization *= axis->GetBinUpEdge(axis->GetLast()) - axis->GetBinLowEdge(axis->GetFirst());
      }
      else
      {
	Printf("Normalizing assuming accepted range of |eta| < 0.8");
	normalization *= 0.8 * 2;
      }
    }
    
    //Printf("Normalization: %f", normalization);
    tracks->Scale(1.0 / normalization);
    
    // NOTE fEventHist contains the number of events for the underlying event analysis and the number of trigger particles for the azimuthal correlation analysis. In the latter case the naming is therefore somewhat misleading!
    TH1D* events = fEventHist->ShowProjection(0, step);
    Int_t nEvents = (Int_t) events->Integral(firstBin, lastBin);
    Printf("Calculated histogram --> %d events", nEvents);
    if (normEvents)
      *normEvents = nEvents;
      
    if (nEvents > 0)
      tracks->Scale(1.0 / nEvents);
    
    delete events;
  }

  ResetBinLimits(fTrackHist[region]->GetGrid(step));
  ResetBinLimits(fEventHist->GetGrid(step));

  return tracks;
}

//____________________________________________________________________
TH2* AliUEHist::GetSumOfRatios(AliUEHist* mixed, AliUEHist::CFStep step, AliUEHist::Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Bool_t etaNorm, Bool_t useVertexBins)
{
  // Calls GetUEHist(...) for *each* multiplicity bin and performs a sum of ratios:
  // 1_N [ (same/mixed)_1 + (same/mixed)_2 + (same/mixed)_3 + ... ]
  // where N is the total number of events/trigger particles and the subscript is the multiplicity bin
  //
  // Can only be used for the 2D histogram at present
  //
  // Parameters:
  //   mixed: AliUEHist containing mixed event corresponding to this object
  //   <other parameters> : check documentation of AliUEHist::GetUEHist
  
  Int_t multIter = multBinBegin;
  
  TH2* totalTracks = 0;
  Int_t totalEvents = 0;
  
  Int_t vertexBin = 1;
  TAxis* vertexAxis = fTrackHist[kToward]->GetGrid(0)->GetGrid()->GetAxis(5);
  if (useVertexBins && !vertexAxis)
  {
    Printf("Vertex axis requested but not available");
    return 0;
  }
  
  // vertex bin loop
  while (1)
  {
    if (useVertexBins)
    {
      SetZVtxRange(vertexAxis->GetBinLowEdge(vertexBin) + 0.01, vertexAxis->GetBinUpEdge(vertexBin) - 0.01);
      mixed->SetZVtxRange(vertexAxis->GetBinLowEdge(vertexBin) + 0.01, vertexAxis->GetBinUpEdge(vertexBin) - 0.01);
      vertexBin++;
    }
    
    // multiplicity loop
    while (1)
    {
      Int_t multBinBeginLocal = multBinBegin;
      Int_t multBinEndLocal = multBinEnd;
      
      if (multBinEnd >= multBinBegin)
      {
	multBinBeginLocal = multIter;
	multBinEndLocal = multIter;
	multIter++;
      }
	
      Int_t nEvents = 0;
      TH2* tracks = (TH2*) GetUEHist(step, region, ptLeadMin, ptLeadMax, multBinBeginLocal, multBinEndLocal, 1, etaNorm, &nEvents);
      // undo normalization
      tracks->Scale(nEvents);
      totalEvents += nEvents;
      
      TH2* mixedTwoD = (TH2*) mixed->GetUEHist(step, region, ptLeadMin, ptLeadMax, multBinBeginLocal, multBinEndLocal, 1, etaNorm);
      
      // asssume flat in dphi, gain in statistics
  //     TH1* histMixedproj = mixedTwoD->ProjectionY();
  //     histMixedproj->Scale(1.0 / mixedTwoD->GetNbinsX());
  //     
  //     for (Int_t x=1; x<=mixedTwoD->GetNbinsX(); x++)
  //       for (Int_t y=1; y<=mixedTwoD->GetNbinsY(); y++)
  // 	mixedTwoD->SetBinContent(x, y, histMixedproj->GetBinContent(y));

  //       delete histMixedproj;
  
      // get mixed event normalization by assuming full acceptance at deta of 0 (only works for flat dphi)
/*      Double_t mixedNorm = mixedTwoD->Integral(1, mixedTwoD->GetNbinsX(), mixedTwoD->GetYaxis()->FindBin(-0.01), mixedTwoD->GetYaxis()->FindBin(0.01));
      mixedNorm /= mixedTwoD->GetNbinsX() * (mixedTwoD->GetYaxis()->FindBin(0.01) - mixedTwoD->GetYaxis()->FindBin(-0.01) + 1);
      tracks->Scale(mixedNorm);*/
      
      tracks->Scale(mixedTwoD->Integral() / tracks->Integral());

      tracks->Divide(mixedTwoD);
      
      delete mixedTwoD;
      
      if (!totalTracks)
	totalTracks = tracks;
      else
      {
	totalTracks->Add(tracks);
	delete tracks;
      }

      if (multIter > multBinEnd)
	break;
    }
    
    if (!useVertexBins || vertexBin > vertexAxis->GetNbins())
      break;
  }

  if (useVertexBins)
    totalEvents = vertexAxis->GetNbins();
  Printf("Dividing %f tracks by %d events", totalTracks->Integral(), totalEvents);
  if (totalEvents > 0)
    totalTracks->Scale(1.0 / totalEvents);
  
  return totalTracks;
}

//____________________________________________________________________
TH1* AliUEHist::GetPtHist(CFStep step, Region region, Float_t ptLeadMin, Float_t ptLeadMax, Int_t multBinBegin, Int_t multBinEnd, Float_t phiMin, Float_t phiMax, Float_t etaMin, Float_t etaMax, Bool_t skipPhiNormalization)
{
  // returns a histogram projected to pT,assoc with the provided cuts
  
   // unzoom all axes
  ResetBinLimits(fTrackHist[region]->GetGrid(step));
  ResetBinLimits(fEventHist->GetGrid(step));
  
  TH1D* tracks = 0;
  
    // the efficiency to have find an event depends on leading pT and this is not corrected for because anyway per bin we calculate tracks over events
  // therefore the number density must be calculated here per leading pT bin and then summed

  if (multBinEnd > multBinBegin)
    Printf("Using multiplicity range %d --> %d", multBinBegin, multBinEnd);
  fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(3)->SetRange(multBinBegin, multBinEnd);
  fEventHist->GetGrid(step)->GetGrid()->GetAxis(1)->SetRange(multBinBegin, multBinEnd);
  
  fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(0)->SetRangeUser(etaMin, etaMax);
  fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(4)->SetRangeUser(phiMin, phiMax);
  
  // get real phi cuts due to binning
  Float_t phiMinReal = fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(4)->GetBinLowEdge(fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(4)->GetFirst());
  Float_t phiMaxReal = fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(4)->GetBinUpEdge(fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(4)->GetLast());
  Printf("phi min = %.2f and max = %.2f requested; due to binning range will be from %.2f to %.2f, i.e. a %.0f%% larger window", phiMin, phiMax, phiMinReal, phiMaxReal, (phiMaxReal - phiMinReal) / (phiMax - phiMin) * 100 - 100);
  
  Int_t firstBin = fTrackHist[region]->GetAxis(2, step)->FindBin(ptLeadMin);
  Int_t lastBin = fTrackHist[region]->GetAxis(2, step)->FindBin(ptLeadMax);
  
  TH1D* events = fEventHist->ShowProjection(0, step);
  
  for (Int_t bin=firstBin; bin<=lastBin; bin++)
  {
    fTrackHist[region]->GetGrid(step)->GetGrid()->GetAxis(2)->SetRange(bin, bin);
    
    // project to pT,assoc
    TH1D* tracksTmp = (TH1D*) fTrackHist[region]->GetGrid(step)->Project(1);
    
    Printf("Calculated histogram in bin %d --> %f tracks", bin, tracksTmp->Integral());
    fTrackHist[region]->GetGrid(step)->SetRangeUser(2, 0, -1);
    
    // normalize to get a density (deta dphi)
    Float_t normalization = 1;
    TAxis* axis = fTrackHist[region]->GetGrid(step)->GetAxis(0);
    if (strcmp(axis->GetTitle(), "#eta") == 0)
    {
      Printf("Normalizing using eta-axis range");
      normalization *= axis->GetBinUpEdge(axis->GetLast()) - axis->GetBinLowEdge(axis->GetFirst());
    }
    else
    {
      Printf("Normalizating assuming accepted range of |eta| < 0.8");
      normalization *= 0.8 * 2;
    }
    
    // dphi
    if (!skipPhiNormalization)
      normalization *= phiMaxReal - phiMinReal;
    
    //Printf("Normalization: %f", normalization);
    tracksTmp->Scale(1.0 / normalization);
    
    // and now dpT (bins have different width)
    for (Int_t i=1; i<=tracksTmp->GetNbinsX(); i++)
    {
      tracksTmp->SetBinContent(i, tracksTmp->GetBinContent(i) / tracksTmp->GetXaxis()->GetBinWidth(i));
      tracksTmp->SetBinError  (i, tracksTmp->GetBinError(i) / tracksTmp->GetXaxis()->GetBinWidth(i));
    }
     
    Int_t nEvents = (Int_t) events->GetBinContent(bin);
    if (nEvents > 0)
      tracksTmp->Scale(1.0 / nEvents);
    Printf("Calculated histogram in bin %d --> %d events", bin, nEvents);
    
    if (!tracks)
      tracks = tracksTmp;
    else
    {
      tracks->Add(tracksTmp);
      delete tracksTmp;
    }
  }
  
  delete events;

  ResetBinLimits(fTrackHist[region]->GetGrid(step));
  ResetBinLimits(fEventHist->GetGrid(step));

  return tracks;
}

void AliUEHist::MultiplyHistograms(THnSparse* grid, THnSparse* target, TH1* histogram, Int_t var1, Int_t var2)
{
  // multiplies <grid> with <histogram> and put the result in <target>
  // <grid> has usually more dimensions than histogram. The axis which are used to choose the value 
  // from <histogram> are given in <var1> and <var2>
  //
  // if <histogram> is 0, just copies content from step1 to step2
  
  // clear target histogram
  target->Reset();
  
  if (histogram != 0)
  {
    if (grid->GetAxis(var1)->GetNbins() != histogram->GetNbinsX())
      AliFatal(Form("Invalid binning (var1): %d %d", grid->GetAxis(var1)->GetNbins(), histogram->GetNbinsX()));
      
    if (var2 >= 0 && grid->GetAxis(var2)->GetNbins() != histogram->GetNbinsY())
      AliFatal(Form("Invalid binning (var2): %d %d", grid->GetAxis(var2)->GetNbins(), histogram->GetNbinsY()));
  }

  if (grid->GetNdimensions() > 6)
    AliFatal("Too many dimensions in THnSparse");
  
  Int_t bins[6];
  
  // optimized implementation
  for (Int_t binIdx = 0; binIdx < grid->GetNbins(); binIdx++)
  {
    Double_t value = grid->GetBinContent(binIdx, bins);
    Double_t error = grid->GetBinError(binIdx);
    
    if (histogram != 0)
    {
      if (var2 < 0)
      {
        value *= histogram->GetBinContent(bins[var1]);
        error *= histogram->GetBinContent(bins[var1]);
      }
      else
      {
        value *= histogram->GetBinContent(bins[var1], bins[var2]);
        error *= histogram->GetBinContent(bins[var1], bins[var2]);
      }
    }
    
    target->SetBinContent(bins, value);
    target->SetBinError(bins, error);
  }
}

//____________________________________________________________________
void AliUEHist::CorrectTracks(CFStep step1, CFStep step2, TH1* trackCorrection, Int_t var1, Int_t var2)
{
  // corrects from step1 to step2 by multiplying the tracks with trackCorrection
  // trackCorrection can be a function of eta (var1 == 0), pT (var1 == 1), leading pT (var1 == 2), multiplicity (var1 == 3), delta phi (var1 == 4)
  // if var2 >= 0 a two dimension correction is assumed in trackCorrection
  //
  // if trackCorrection is 0, just copies content from step1 to step2
  
  for (UInt_t region=0; region<fkRegions; region++)
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
  
  MultiplyHistograms(grid, target, trackCorrection, var1, var2);
  
  Printf("AliUEHist::CorrectTracks: Corrected from %f to %f entries. Correction histogram: %f entries (integral: %f)", grid->GetEntries(), target->GetEntries(), (trackCorrection) ? trackCorrection->GetEntries() : -1.0, (trackCorrection) ? trackCorrection->Integral() : -1.0); 
}

//____________________________________________________________________
void AliUEHist::CorrectEvents(CFStep step1, CFStep step2, TH1* eventCorrection, Int_t var1, Int_t var2)
{
  // corrects from step1 to step2 by multiplying the events with eventCorrection
  // eventCorrection is as function of leading pT (var == 0) or multiplicity (var == 1)
  //
  // if eventCorrection is 0, just copies content from step1 to step2
  
  AliCFGridSparse* grid = fEventHist->GetGrid(step1);
  AliCFGridSparse* target = fEventHist->GetGrid(step2);
  
  MultiplyHistograms(grid->GetGrid(), target->GetGrid(), eventCorrection, var1, var2);

  Printf("AliUEHist::CorrectEvents: Corrected from %f to %f entries. Correction histogram: %f entries (integral: %f)", grid->GetEntries(), target->GetEntries(), (eventCorrection) ? eventCorrection->GetEntries() : -1.0, (eventCorrection) ? eventCorrection->Integral() : -1.0); 
}

//____________________________________________________________________
void AliUEHist::Correct(AliUEHist* corrections)
{
  // applies the given corrections to extract from the step kCFStepReconstructed all previous steps
  //
  // in this object the data is expected in the step kCFStepReconstructed
  
  if (fHistogramType.Length() == 0)
  {
    Printf("W-AliUEHist::Correct: fHistogramType not defined. Guessing histogram type...");
    if (fTrackHist[kToward]->GetNVar() < 5)
    {
      if (strcmp(fTrackHist[kToward]->GetTitle(), "d^{2}N_{ch}/d#varphid#eta") == 0)
        fHistogramType = "NumberDensitypT";
      else if (strcmp(fTrackHist[kToward]->GetTitle(), "d^{2}#Sigma p_{T}/d#varphid#eta") == 0)
        fHistogramType = "SumpT";
    }
    else if (fTrackHist[kToward]->GetNVar() == 5)
    {
      if (strcmp(fTrackHist[kToward]->GetGrid(0)->GetGrid()->GetAxis(3)->GetTitle(), "multiplicity") == 0)
        fHistogramType = "NumberDensityPhi";
      else if (strcmp(fTrackHist[kToward]->GetGrid(0)->GetGrid()->GetAxis(3)->GetTitle(), "centrality") == 0)
        fHistogramType = "NumberDensityPhiCentrality";
    }
      
    if (fHistogramType.Length() == 0)
      AliFatal("Cannot figure out histogram type. Try setting it manually...");
  }
  
  Printf("AliUEHist::Correct: Correcting %s...", fHistogramType.Data());
  
  if (strcmp(fHistogramType, "NumberDensitypT") == 0 || strcmp(fHistogramType, "NumberDensityPhi") == 0 || strcmp(fHistogramType, "SumpT") == 0)
  {
    // ---- track level
    
    // bias due to migration in leading pT (because the leading particle is not reconstructed, and the subleading is used)
    // extracted as function of leading pT
    Bool_t biasFromMC = kFALSE;
    const char* projAxis = "z";
    Int_t secondBin = -1;

    if (strcmp(fHistogramType, "NumberDensityPhi") == 0)
    {
      projAxis = "yz";
      secondBin = 4;
    }
    
    for (UInt_t region = 0; region < fkRegions; region++)
    {
      if (!fTrackHist[region])
        continue;
   
      TH1* leadingBiasTracks =  0;
      if (biasFromMC)
      {
        leadingBiasTracks = (TH1*) corrections->GetBias(kCFStepReconstructed, kCFStepTracked, region, projAxis, 0, -1, 1); // from MC
        Printf("WARNING: Using MC bias correction");
      }
      else
        leadingBiasTracks = (TH1*) GetBias(kCFStepBiasStudy, kCFStepReconstructed, region, projAxis, 0, -1, 1);          // from data
        
      CorrectTracks(kCFStepReconstructed, kCFStepTracked, region, leadingBiasTracks, 2, secondBin);
      if (region == kMin && fCombineMinMax)
      {
        CorrectTracks(kCFStepReconstructed, kCFStepTracked, kMax, leadingBiasTracks, 2, secondBin);
        delete leadingBiasTracks;
        break;
      }
      delete leadingBiasTracks;
    }
    
    TH1* leadingBiasEvents = 0;
    if (biasFromMC)
      leadingBiasEvents = (TH1*) corrections->GetBias(kCFStepReconstructed, kCFStepTracked, kToward, projAxis, 0, -1, 2); // from MC
    else
      leadingBiasEvents = (TH1*) GetBias(kCFStepBiasStudy, kCFStepReconstructed, kToward, projAxis, 0, -1, 2);          // from data
      
    //new TCanvas; leadingBiasEvents->DrawCopy();
    CorrectEvents(kCFStepReconstructed, kCFStepTracked, leadingBiasEvents, 0);
    delete leadingBiasEvents;
    
    // --- contamination correction ---
    
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
    
    // --- efficiency correction ---
    // correct single-particle efficiency for associated particles
    // in addition correct for efficiency on trigger particles (tracks AND events)
    
    TH1* efficiencyCorrection = corrections->GetTrackingEfficiencyCorrection();
    // use kCFStepVertex as a temporary step (will be overwritten later anyway)
    CorrectTracks(kCFStepTrackedOnlyPrim, kCFStepVertex, efficiencyCorrection, 0, 1);
    delete efficiencyCorrection;
    
    // correct pT,T
    efficiencyCorrection = corrections->GetTrackEfficiency(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, 1, -1, 2);
    CorrectEvents(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, efficiencyCorrection, 0);
    CorrectTracks(kCFStepVertex, kCFStepAnaTopology, efficiencyCorrection, 2);
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
    // some defaults
    func->SetParameters(0.1, 1, -0.7);
    vertexCorrection->Fit(func, "0I", "", 0, 3);
    for (Int_t i=1; i<=vertexCorrectionObs->GetNbinsX(); i++)
    {
      Float_t xPos = 1.0 / 0.77 * vertexCorrectionObs->GetXaxis()->GetBinCenter(i);
      if (xPos < 3)
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
  else if (strcmp(fHistogramType, "NumberDensityPhiCentrality") == 0)
  {
    // copy 
    CorrectTracks(kCFStepReconstructed, kCFStepTracked, 0, -1);
    CorrectEvents(kCFStepReconstructed, kCFStepTracked, 0, -1);
    
    // Dont use eta in the following, because it is a Delta-eta axis
    
    // contamination correction
    // correct single-particle contamination for associated particles
    
    TH1* contamination = corrections->GetTrackingContamination(1);
    
    if (1)
    {
      Printf("Applying contamination enhancement");
      
      for (Int_t bin = 1; bin <= contamination->GetNbinsX(); bin++)
      {
        printf("%f", contamination->GetBinContent(bin));
        if (contamination->GetBinContent(bin) > 0)
          contamination->SetBinContent(bin, 1.0 + 1.1 * (contamination->GetBinContent(bin) - 1.0));
        printf(" --> %f\n", contamination->GetBinContent(bin));
      }
    }
      
    CorrectTracks(kCFStepTracked, kCFStepTrackedOnlyPrim, contamination, 1);
    delete contamination;    
    
    // correct for additional contamination due to trigger particle around phi ~ 0
    TH2* correlatedContamination = corrections->GetCorrelatedContamination();
    if (1)
    {
      Printf("Applying contamination enhancement");
      
      for (Int_t bin = 1; bin <= correlatedContamination->GetNbinsX(); bin++)
        for (Int_t bin2 = 1; bin2 <= correlatedContamination->GetNbinsY(); bin2++)
        {
          printf("%f", correlatedContamination->GetBinContent(bin, bin2));
          if (correlatedContamination->GetBinContent(bin, bin2) > 0)
            correlatedContamination->SetBinContent(bin, bin2, 1.0 + 1.1 * (correlatedContamination->GetBinContent(bin, bin2) - 1.0));
          printf(" --> %f\n", correlatedContamination->GetBinContent(bin, bin2));
        }
    }
    
    new TCanvas; correlatedContamination->DrawCopy("COLZ");
    CorrectCorrelatedContamination(kCFStepTrackedOnlyPrim, 0, correlatedContamination);
//     Printf("\n\n\nWARNING ---> SKIPPING CorrectCorrelatedContamination\n\n\n");
    
    delete correlatedContamination;
    
    // TODO correct for contamination of trigger particles (for tracks AND events)
    CorrectEvents(kCFStepTracked, kCFStepTrackedOnlyPrim, 0, 0);
    
    // --- efficiency correction ---
    // correct single-particle efficiency for associated particles
    // in addition correct for efficiency on trigger particles (tracks AND events)
    
    // in bins of pT and centrality
    TH1* efficiencyCorrection = corrections->GetTrackingEfficiencyCorrectionCentrality();
    // use kCFStepAnaTopology as a temporary step 
    CorrectTracks(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, efficiencyCorrection, 1, 3);
    delete efficiencyCorrection;
    
    // correct pT,T in bins of pT and centrality
    efficiencyCorrection = corrections->GetTrackEfficiency(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, 1, 3, 2);
    CorrectEvents(kCFStepTrackedOnlyPrim, kCFStepVertex, efficiencyCorrection, 0, 1);
    CorrectTracks(kCFStepAnaTopology, kCFStepVertex, efficiencyCorrection, 2, 3);
    delete efficiencyCorrection;
    
    // no correction for vertex finding efficiency and trigger efficiency needed in PbPb
    // copy 
    CorrectTracks(kCFStepVertex, kCFStepAll, 0, -1);
    CorrectEvents(kCFStepVertex, kCFStepAll, 0, -1);
  }
  else
    AliFatal(Form("Unknown histogram for correction: %s", GetTitle()));
}

//____________________________________________________________________
TH1* AliUEHist::GetTrackEfficiency(CFStep step1, CFStep step2, Int_t axis1, Int_t axis2, Int_t source, Int_t axis3)
{
  // creates a track-level efficiency by dividing step2 by step1
  // projected to axis1 and axis2 (optional if >= 0)
  //
  // source: 0 = fTrackHist; 1 = fTrackHistEfficiency; 2 = fTrackHistEfficiency rebinned for pT,T / pT,lead binning
  
  // integrate over regions
  // cache it for efficiency (usually more than one efficiency is requested)
  
  AliCFContainer* sourceContainer = 0;
  
  if (source == 0)
  {
    if (!fCache)
    {
      fCache = (AliCFContainer*) fTrackHist[0]->Clone();
      for (UInt_t i = 1; i < fkRegions; i++)
        if (fTrackHist[i])
          fCache->Add(fTrackHist[i]);
    }
    sourceContainer = fCache;
  }
  else if (source == 1 || source == 2)
  {
    sourceContainer = fTrackHistEfficiency;
    // step offset because we start with kCFStepAnaTopology
    step1 = (CFStep) ((Int_t) step1 - (Int_t) kCFStepAnaTopology);
    step2 = (CFStep) ((Int_t) step2 - (Int_t) kCFStepAnaTopology);
  }
  else
    return 0;
        
  // reset all limits and set the right ones except those in axis1, axis2 and axis3
  ResetBinLimits(sourceContainer->GetGrid(step1));
  ResetBinLimits(sourceContainer->GetGrid(step2));
  if (fEtaMax > fEtaMin && axis1 != 0 && axis2 != 0 && axis3 != 0)
  {
    Printf("Restricted eta-range to %f %f", fEtaMin, fEtaMax);
    sourceContainer->GetGrid(step1)->SetRangeUser(0, fEtaMin, fEtaMax);
    sourceContainer->GetGrid(step2)->SetRangeUser(0, fEtaMin, fEtaMax);
  }
  if (fPtMax > fPtMin && axis1 != 1 && axis2 != 1 && axis3 != 1)
  {
    sourceContainer->GetGrid(step1)->SetRangeUser(1, fPtMin, fPtMax);
    sourceContainer->GetGrid(step2)->SetRangeUser(1, fPtMin, fPtMax);
  }
  if (fCentralityMax > fCentralityMin && axis1 != 3 && axis2 != 3 && axis3 != 3)
  {
    sourceContainer->GetGrid(step1)->SetRangeUser(3, fCentralityMin, fCentralityMax);
    sourceContainer->GetGrid(step2)->SetRangeUser(3, fCentralityMin, fCentralityMax);
  }
  
  TH1* measured = 0;
  TH1* generated = 0;
    
  if (axis3 >= 0)
  {
    generated = sourceContainer->Project(step1, axis1, axis2, axis3);
    measured = sourceContainer->Project(step2, axis1, axis2, axis3);
  }
  else if (axis2 >= 0)
  {
    generated = sourceContainer->Project(step1, axis1, axis2);
    measured = sourceContainer->Project(step2, axis1, axis2);
  }
  else
  {
    generated = sourceContainer->Project(step1, axis1);
    measured = sourceContainer->Project(step2, axis1);
  }
  
  // check for bins with less than 50 entries, print warning
  Int_t binBegin[3];
  Int_t binEnd[3];
  
  binBegin[0] = 1;
  binBegin[1] = 1;
  binBegin[2] = 1;
  
  binEnd[0] = generated->GetNbinsX();
  binEnd[1] = generated->GetNbinsY();
  binEnd[2] = generated->GetNbinsZ();
  
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
    if (axis3 == 0)
    {
      binBegin[2] = generated->GetZaxis()->FindBin(fEtaMin);
      binEnd[2]   = generated->GetZaxis()->FindBin(fEtaMax);
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
    if (axis3 == 1)
    {
      binBegin[2] = generated->GetZaxis()->FindBin(fPtMin);
      binEnd[2]   = generated->GetZaxis()->FindBin(ptMax);
    }
  }
  
  Int_t total = 0;
  Int_t count = 0;
  Int_t vars[3];
  
  for (Int_t i=0; i<3; i++)
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
    else if (generated->GetDimension() == 3 && generated->GetBinContent(vars[0], vars[1], vars[2]) < limit)
    {
      Printf("Empty bin at %s=%.2f %s=%.2f %s=%.2f (%.2f entries)", 
        generated->GetXaxis()->GetTitle(), generated->GetXaxis()->GetBinCenter(vars[0]),
        generated->GetYaxis()->GetTitle(), generated->GetYaxis()->GetBinCenter(vars[1]),
        generated->GetZaxis()->GetTitle(), generated->GetZaxis()->GetBinCenter(vars[2]),
        generated->GetBinContent(vars[0], vars[1], vars[2]));
      count++;
    }
    
    vars[2]++;
    if (vars[2] == binEnd[2]+1)
    {
      vars[2] = binBegin[2];
      vars[1]++;
    }
    
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
  
  // rebin if required
  if (source == 2)
  {
    TAxis* axis = fEventHist->GetGrid(0)->GetGrid()->GetAxis(0);
    
    if (axis->GetNbins() < measured->GetNbinsX())
    {
      if (axis2 != -1)
      {
        // 2d rebin with variable axis does not exist in root
        
        TH1* tmp = measured;
        measured = new TH2D(Form("%s_rebinned", tmp->GetName()), tmp->GetTitle(), axis->GetNbins(), axis->GetXbins()->GetArray(), tmp->GetNbinsY(), tmp->GetYaxis()->GetXbins()->GetArray());
        for (Int_t x=1; x<=tmp->GetNbinsX(); x++)
          for (Int_t y=1; y<=tmp->GetNbinsY(); y++)
          {
            ((TH2*) measured)->Fill(tmp->GetXaxis()->GetBinCenter(x), tmp->GetYaxis()->GetBinCenter(y), tmp->GetBinContent(x, y));
            measured->SetBinError(x, y, 0); // cannot trust bin error, set to 0
          }
        delete tmp;
        
        tmp = generated;
        generated = new TH2D(Form("%s_rebinned", tmp->GetName()), tmp->GetTitle(), axis->GetNbins(), axis->GetXbins()->GetArray(), tmp->GetNbinsY(), tmp->GetYaxis()->GetXbins()->GetArray());
        for (Int_t x=1; x<=tmp->GetNbinsX(); x++)
          for (Int_t y=1; y<=tmp->GetNbinsY(); y++)
          {
            ((TH2*) generated)->Fill(tmp->GetXaxis()->GetBinCenter(x), tmp->GetYaxis()->GetBinCenter(y), tmp->GetBinContent(x, y));
            generated->SetBinError(x, y, 0); // cannot trust bin error, set to 0
          }
        delete tmp;
      }
      else
      {
        TH1* tmp = measured;
        measured = tmp->Rebin(axis->GetNbins(), Form("%s_rebinned", tmp->GetName()), axis->GetXbins()->GetArray());
        delete tmp;
        
        tmp = generated;
        generated = tmp->Rebin(axis->GetNbins(), Form("%s_rebinned", tmp->GetName()), axis->GetXbins()->GetArray());
        delete tmp;
      }
    }
    else if (axis->GetNbins() > measured->GetNbinsX())
    {
      if (axis2 != -1)
        AliFatal("Rebinning only works for 1d at present");
  
      // this is an unfortunate case where the number of bins has to be increased in principle
      // there is a region where the binning is finner in one histogram and a region where it is the other way round
      // this cannot be resolved in principle, but as we only calculate the ratio the bin in the second region get the same entries
      // only a certain binning is supported here
      if (axis->GetNbins() != 100 || measured->GetNbinsX() != 39)
        AliFatal(Form("Invalid binning --> %d %d", axis->GetNbins(), measured->GetNbinsX()));
      
      Double_t newBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
      
      // reduce binning below 5 GeV/c
      TH1* tmp = measured;
      measured = tmp->Rebin(27, Form("%s_rebinned", tmp->GetName()), newBins);
      delete tmp;
      
      // increase binning above 5 GeV/c
      tmp = measured;
      measured = new TH1F(Form("%s_rebinned2", tmp->GetName()), tmp->GetTitle(), axis->GetNbins(), axis->GetBinLowEdge(1), axis->GetBinUpEdge(axis->GetNbins()));
      for (Int_t bin=1; bin<=measured->GetNbinsX(); bin++)
      {
        measured->SetBinContent(bin, tmp->GetBinContent(tmp->FindBin(measured->GetBinCenter(bin))));
        measured->SetBinError(bin, tmp->GetBinError(tmp->FindBin(measured->GetBinCenter(bin))));
      }
      delete tmp;
      
      // reduce binning below 5 GeV/c
      tmp = generated;
      generated = tmp->Rebin(27, Form("%s_rebinned", tmp->GetName()), newBins);
      delete tmp;
      
      // increase binning above 5 GeV/c
      tmp = generated;
      generated = new TH1F(Form("%s_rebinned2", tmp->GetName()), tmp->GetTitle(), axis->GetNbins(), axis->GetBinLowEdge(1), axis->GetBinUpEdge(axis->GetNbins()));
      for (Int_t bin=1; bin<=generated->GetNbinsX(); bin++)
      {
        generated->SetBinContent(bin, tmp->GetBinContent(tmp->FindBin(generated->GetBinCenter(bin))));
        generated->SetBinError(bin, tmp->GetBinError(tmp->FindBin(generated->GetBinCenter(bin))));
      }
      delete tmp;
    }
  }
  
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
    generated = fEventHist->Project(step1, axis1, axis2);
    measured = fEventHist->Project(step2, axis1, axis2);
  }
  else
  {
    generated = fEventHist->Project(step1, axis1);
    measured = fEventHist->Project(step2, axis1);
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
TH1* AliUEHist::GetBias(CFStep step1, CFStep step2, Int_t region, const char* axis, Float_t leadPtMin, Float_t leadPtMax, Int_t weighting)
{
  // extracts the track-level bias (integrating out the multiplicity) between two steps (dividing step2 by step1)
  // in the given region (sum over all regions is calculated if region == -1)
  // done by weighting the track-level distribution with the number of events as function of leading pT
  // and then calculating the ratio between the distributions
  // projected to axis which is a TH3::Project3D string, e.g. "x", or "yx"
  //   no projection is done if axis == 0
  // weighting: 0 = tracks weighted with events (as discussed above)
  //            1 = only track bias is returned
  //            2 = only event bias is returned
  
  AliCFContainer* tmp = 0;
  
  if (region == -1)
  {
    tmp = (AliCFContainer*) fTrackHist[0]->Clone();
    for (UInt_t i = 1; i < fkRegions; i++)
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
  
  TH1D* events1 = (TH1D*)fEventHist->Project(step1, 0);
  TH3D* hist1 = (TH3D*)tmp->Project(step1, 0, tmp->GetNVar()-1, 2);
  if (weighting == 0)
    WeightHistogram(hist1, events1);
  
  TH1D* events2 = (TH1D*)fEventHist->Project(step2, 0);
  TH3D* hist2 = (TH3D*)tmp->Project(step2, 0, tmp->GetNVar()-1, 2);
  if (weighting == 0)
    WeightHistogram(hist2, events2);
  
  TH1* generated = hist1;
  TH1* measured = hist2;
  
  if (weighting == 0 || weighting == 1)
  {
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
    delete events1;
    delete events2;
  }
  else if (weighting == 2)
  {
    delete hist1;
    delete hist2;
    generated = events1;
    measured = events2;
  }
  
  measured->Divide(generated);
  
  delete generated;
  
  ResetBinLimits(tmp->GetGrid(step1));
  ResetBinLimits(tmp->GetGrid(step2));
  
  if ((region == -1) || (region == kMin && fCombineMinMax))
    delete tmp;
  
  return measured;
}

//____________________________________________________________________
void AliUEHist::CorrectCorrelatedContamination(CFStep step, Int_t region, TH1* trackCorrection)
{
  // corrects for the given factor in a small delta-eta and delta-phi window as function of pT,A and pT,T
  
  if (!fTrackHist[region])
    return;
   
  THnSparse* grid = fTrackHist[region]->GetGrid(step)->GetGrid();
  
  Int_t var1 = 1;
  Int_t var2 = 2;
  
  if (grid->GetAxis(var1)->GetNbins() != trackCorrection->GetNbinsX())
    AliFatal(Form("Invalid binning (var1): %d %d", grid->GetAxis(var1)->GetNbins(), trackCorrection->GetNbinsX()));
    
  if (grid->GetAxis(var2)->GetNbins() != trackCorrection->GetNbinsY())
    AliFatal(Form("Invalid binning (var2): %d %d", grid->GetAxis(var2)->GetNbins(), trackCorrection->GetNbinsY()));
  
  // optimized implementation
  for (Int_t binIdx = 0; binIdx < grid->GetNbins(); binIdx++)
  {
    Int_t bins[6];
    
    Double_t value = grid->GetBinContent(binIdx, bins);
    Double_t error = grid->GetBinError(binIdx);
    
    // check eta and phi axes
    if (TMath::Abs(grid->GetAxis(0)->GetBinCenter(bins[0])) > 0.1)
      continue;
    if (TMath::Abs(grid->GetAxis(4)->GetBinCenter(bins[4])) > 0.1)
      continue;
    
    value *= trackCorrection->GetBinContent(bins[var1], bins[var2]);
    error *= trackCorrection->GetBinContent(bins[var1], bins[var2]);
    
    grid->SetBinContent(bins, value);
    grid->SetBinError(bins, error);
  }
 
  Printf("AliUEHist::CorrectCorrelatedContamination: Corrected.");
}

//____________________________________________________________________
TH2* AliUEHist::GetCorrelatedContamination()
{
  // contamination correlated with the trigger particle is evaluated between step kCFStepTracked and kCFStepTrackedOnlyPrim in the region of delta eta and delta phi < 0.1 (smallest bin!)
  
  Int_t step1 = kCFStepTrackedOnlyPrim;
  Int_t step2 = kCFStepTracked;
  
  fTrackHist[0]->GetGrid(step1)->SetRangeUser(0, -0.01, 0.01); // delta eta
  fTrackHist[0]->GetGrid(step1)->SetRangeUser(4, -0.01, 0.01); // delta phi
  TH2* tracksStep1 = (TH2*) fTrackHist[0]->Project(step1, 1, 2);
  
  fTrackHist[0]->GetGrid(step2)->SetRangeUser(0, -0.01, 0.01); // delta eta
  fTrackHist[0]->GetGrid(step2)->SetRangeUser(4, -0.01, 0.01); // delta phi
  TH2* tracksStep2 = (TH2*) fTrackHist[0]->Project(step2, 1, 2);
  
  tracksStep1->Divide(tracksStep2);
  
  TH1* triggersStep1 = fEventHist->Project(step1, 0);
  TH1* triggersStep2 = fEventHist->Project(step2, 0);
  
  TH1* singleParticle = GetTrackingContamination(1);
  
  for (Int_t x=1; x<=tracksStep1->GetNbinsX(); x++)
    for (Int_t y=1; y<=tracksStep1->GetNbinsY(); y++)
      if (singleParticle->GetBinContent(x) > 0 && triggersStep1->GetBinContent(y) > 0)
        tracksStep1->SetBinContent(x, y, tracksStep1->GetBinContent(x, y) / triggersStep1->GetBinContent(y) * triggersStep2->GetBinContent(y) / singleParticle->GetBinContent(x));
      else
        tracksStep1->SetBinContent(x, y, 0);
        
  delete singleParticle;
  delete tracksStep2;
  delete triggersStep1;
  delete triggersStep2;
        
  return tracksStep1;
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
TH2D* AliUEHist::GetTrackingEfficiencyCentrality()
{
  // extracts the tracking efficiency by calculating the efficiency from step kCFStepAnaTopology to kCFStepTrackedOnlyPrim
  // integrates over the regions and all other variables than pT, centrality to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepAnaTopology, kCFStepTrackedOnlyPrim, 1, 3));
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
TH2D* AliUEHist::GetTrackingEfficiencyCorrectionCentrality()
{
  // extracts the tracking correction by calculating the efficiency from step kCFStepAnaTopology to kCFStepTracked
  // integrates over the regions and all other variables than pT and centrality to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepTrackedOnlyPrim, kCFStepAnaTopology, 1, 3));
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
TH2D* AliUEHist::GetTrackingContaminationCentrality()
{
  // extracts the tracking contamination by secondaries by calculating the efficiency from step kCFStepTrackedOnlyPrim to kCFStepTracked
  // integrates over the regions and all other variables than pT and centrality to increase the statistics
  //
  // returned histogram has to be deleted by the user

  return dynamic_cast<TH2D*> (GetTrackEfficiency(kCFStepTracked, kCFStepTrackedOnlyPrim, 1, 3));
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
    //fTrackHist[region]->SetGrid(AliUEHist::kCFStepTrackedOnlyPrim, from->fTrackHist[region]->GetGrid(AliUEHist::kCFStepTrackedOnlyPrim));
    fTrackHist[region]->SetGrid(AliUEHist::kCFStepBiasStudy,     from->fTrackHist[region]->GetGrid(AliUEHist::kCFStepBiasStudy));
  }
    
  fEventHist->SetGrid(AliUEHist::kCFStepReconstructed, from->fEventHist->GetGrid(AliUEHist::kCFStepReconstructed));
  //fEventHist->SetGrid(AliUEHist::kCFStepTrackedOnlyPrim, from->fEventHist->GetGrid(AliUEHist::kCFStepTrackedOnlyPrim));
  fEventHist->SetGrid(AliUEHist::kCFStepBiasStudy,     from->fEventHist->GetGrid(AliUEHist::kCFStepBiasStudy));
}

//____________________________________________________________________
void AliUEHist::ExtendTrackingEfficiency(Bool_t verbose)
{
  // fits the tracking efficiency at high pT with a constant and fills all bins with this tracking efficiency

  Float_t fitRangeBegin = 5.01;
  Float_t fitRangeEnd = 14.99;
  Float_t extendRangeBegin = 10.01;

  if (fTrackHistEfficiency->GetNVar() == 3)
  {
    TH1* obj = GetTrackingEfficiency(1);
  
    if (verbose)
    {
      new TCanvas; 
      obj->Draw();
    }
    
    obj->Fit("pol0", (verbose) ? "+" : "0+", "SAME", fitRangeBegin, fitRangeEnd);
  
    Float_t trackingEff = obj->GetFunction("pol0")->GetParameter(0);
  
    obj = GetTrackingContamination(1);
  
    if (verbose)
    {
      new TCanvas; 
      obj->Draw();
    }
    
    obj->Fit("pol0", (verbose) ? "+" : "0+", "SAME", fitRangeBegin, fitRangeEnd);
  
    Float_t trackingCont = obj->GetFunction("pol0")->GetParameter(0);
  
    Printf("AliUEHist::ExtendTrackingEfficiency: Fitted efficiency between %f and %f and got %f tracking efficiency and %f tracking contamination correction. Extending from %f onwards (within %f < eta < %f)", fitRangeBegin, fitRangeEnd, trackingEff, trackingCont, extendRangeBegin, fEtaMin, fEtaMax);
  
    // extend for full pT range
    for (Int_t x = fTrackHistEfficiency->GetAxis(0, 0)->FindBin(fEtaMin); x <= fTrackHistEfficiency->GetAxis(0, 0)->FindBin(fEtaMax); x++)
      for (Int_t y = fTrackHistEfficiency->GetAxis(1, 0)->FindBin(extendRangeBegin); y <= fTrackHistEfficiency->GetNBins(1); y++)
        for (Int_t z = 1; z <= fTrackHistEfficiency->GetNBins(2); z++) // particle type axis
        {
          
          Int_t bins[3];
          bins[0] = x;
          bins[1] = y;
          bins[2] = z;
          
          fTrackHistEfficiency->GetGrid(0)->SetElement(bins, 100);
          fTrackHistEfficiency->GetGrid(1)->SetElement(bins, 100.0 * trackingEff);
          fTrackHistEfficiency->GetGrid(2)->SetElement(bins, 100.0 * trackingEff / trackingCont);
        }
  }
  else if (fTrackHistEfficiency->GetNVar() == 4)
  {
    // fit in centrality intervals of 20% for efficiency, one bin for contamination
    Float_t* trackingEff = 0;
    Float_t* trackingCont = 0;
    Float_t centralityBins[] = { 0, 10, 20, 40, 60, 100 };
    Int_t nCentralityBins = 5;
    
    Printf("AliUEHist::ExtendTrackingEfficiency: Fitting efficiencies between %f and %f. Extending from %f onwards (within %f < eta < %f)", fitRangeBegin, fitRangeEnd, extendRangeBegin, fEtaMin, fEtaMax);
    
    // 0 = eff; 1 = cont
    for (Int_t caseNo = 0; caseNo < 2; caseNo++)
    {
      Float_t* target = 0;
      Int_t centralityBinsLocal = nCentralityBins;
      
      if (caseNo == 0)
      {
        trackingEff = new Float_t[centralityBinsLocal];
        target = trackingEff;
      }
      else
      {
        centralityBinsLocal = 1;
        trackingCont = new Float_t[centralityBinsLocal];
        target = trackingCont;
      }
    
      for (Int_t i=0; i<centralityBinsLocal; i++)
      {
	if (centralityBinsLocal == 1)
	  SetCentralityRange(centralityBins[0] + 0.1, centralityBins[nCentralityBins] - 0.1);
	else
	  SetCentralityRange(centralityBins[i] + 0.1, centralityBins[i+1] - 0.1);
        TH1* proj = (caseNo == 0) ? GetTrackingEfficiency(1) : GetTrackingContamination(1);
        if (verbose)
        {
          new TCanvas;
          proj->DrawCopy();
        }
        if ((Int_t) proj->Fit("pol0", (verbose) ? "+" : "Q0+", "SAME", fitRangeBegin, fitRangeEnd) == 0)
          target[i] = proj->GetFunction("pol0")->GetParameter(0);
        else
          target[i] = 0;
        Printf("AliUEHist::ExtendTrackingEfficiency: case %d, centrality %d, eff %f", caseNo, i, target[i]);
      }
    }
  
    // extend for full pT range
    for (Int_t x = fTrackHistEfficiency->GetAxis(0, 0)->FindBin(fEtaMin); x <= fTrackHistEfficiency->GetAxis(0, 0)->FindBin(fEtaMax); x++)
      for (Int_t y = fTrackHistEfficiency->GetAxis(1, 0)->FindBin(extendRangeBegin); y <= fTrackHistEfficiency->GetNBins(1); y++)
        for (Int_t z = 1; z <= fTrackHistEfficiency->GetNBins(2); z++) // particle type axis
        {
          for (Int_t z2 = 1; z2 <= fTrackHistEfficiency->GetNBins(3); z2++) // centrality axis
          {
            
            Int_t bins[4];
            bins[0] = x;
            bins[1] = y;
            bins[2] = z;
            bins[3] = z2;
            
            Int_t z2Bin = 0;
	    while (centralityBins[z2Bin+1] < fTrackHistEfficiency->GetAxis(3, 0)->GetBinCenter(z2))
	      z2Bin++;
	    
            //Printf("%d %d", z2, z2Bin);
            
            fTrackHistEfficiency->GetGrid(0)->SetElement(bins, 100);
            fTrackHistEfficiency->GetGrid(1)->SetElement(bins, 100.0 * trackingEff[z2Bin]);
            if (trackingCont[0] > 0)
              fTrackHistEfficiency->GetGrid(2)->SetElement(bins, 100.0 * trackingEff[z2Bin] / trackingCont[0]);
            else  
              fTrackHistEfficiency->GetGrid(2)->SetElement(bins, 0);
          }
       }
       
     delete[] trackingEff;
     delete[] trackingCont;
   }
   
   SetCentralityRange(0, 0);
}

void AliUEHist::AdditionalDPhiCorrection(Int_t step)
{
  // corrects the dphi distribution with an extra factor close to dphi ~ 0

  Printf("WARNING: In AliUEHist::AdditionalDPhiCorrection.");

  THnSparse* grid = fTrackHist[0]->GetGrid(step)->GetGrid();
  
  // optimized implementation
  for (Int_t binIdx = 0; binIdx < grid->GetNbins(); binIdx++)
  {
    Int_t bins[6];
    Double_t value = grid->GetBinContent(binIdx, bins);
    Double_t error = grid->GetBinError(binIdx);
    
    Float_t binCenter = grid->GetAxis(4)->GetBinCenter(bins[4]);
    if (TMath::Abs(binCenter) < 0.2)
    {
      value *= 0.985;
      error *= 0.985;
    }
    else if (TMath::Abs(binCenter) < 0.3)
    {
      value *= 0.9925;
      error *= 0.9925;
    }
    
    grid->SetBinContent(bins, value);
    grid->SetBinError(bins, error);
  }
}

void AliUEHist::Scale(Double_t factor)
{
  // scales all contained histograms by the given factor
  
  for (Int_t i=0; i<4; i++)
    if (fTrackHist[i])
      fTrackHist[i]->Scale(factor);
  
  fEventHist->Scale(factor);
  fTrackHistEfficiency->Scale(factor);
}

void AliUEHist::Reset()
{
  // resets all contained histograms
  
  for (Int_t i=0; i<4; i++)
    if (fTrackHist[i])
      for (Int_t step=0; step<fTrackHist[i]->GetNStep(); step++)
        fTrackHist[i]->GetGrid(step)->GetGrid()->Reset();
  
  for (Int_t step=0; step<fEventHist->GetNStep(); step++)
    fEventHist->GetGrid(step)->GetGrid()->Reset();
    
  for (Int_t step=0; step<fTrackHistEfficiency->GetNStep(); step++)
    fTrackHistEfficiency->GetGrid(step)->GetGrid()->Reset();
}
