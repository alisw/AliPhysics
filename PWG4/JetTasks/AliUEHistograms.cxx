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

/* $Id: AliUEHistograms.cxx 20164 2007-08-14 15:31:50Z morsch $ */

//
//
// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
//
//
// Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

#include "AliUEHistograms.h"

#include "AliCFContainer.h"
#include "AliVParticle.h"

#include "TList.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"

ClassImp(AliUEHistograms)

AliUEHistograms::AliUEHistograms() : 
  TObject(),
  fNumberDensitypT(0),
  fSumpT(0),
  fNumberDensityPhi(0),
  fCorrelationpT(0),
  fCorrelationEta(0),
  fCorrelationPhi(0),
  fCorrelationR(0),
  fCorrelationLeading2Phi(0),
  fCorrelationMultiplicity(0),
  fEventCount(0),
  fEventCountDifferential(0),
  fVertexContributors(0)
{
  // Constructor
  
  fNumberDensitypT = new AliUEHist("NumberDensitypT");
  fSumpT = new AliUEHist("SumpT");
  fNumberDensityPhi = new AliUEHist("NumberDensityPhi");
  
  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fCorrelationpT  = new TH2F("fCorrelationpT", ";p_{T,lead} (MC);p_{T,lead} (RECO)", 200, 0, 50, 200, 0, 50);
  fCorrelationEta = new TH2F("fCorrelationEta", ";#eta_{T,lead} (MC);#eta_{T,lead} (RECO)", 200, -1, 1, 200, -1, 1);
  fCorrelationPhi = new TH2F("fCorrelationPhi", ";#phi_{T,lead} (MC);#phi_{T,lead} (RECO)", 200, 0, TMath::TwoPi(), 200, 0, TMath::TwoPi());
  fCorrelationR =   new TH2F("fCorrelationR", ";R;p_{T,lead} (MC)", 200, 0, 2, 200, 0, 50);
  fCorrelationLeading2Phi = new TH2F("fCorrelationLeading2Phi", ";#Delta #phi;p_{T,lead} (MC)", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 50);
  fCorrelationMultiplicity = new TH2F("fCorrelationMultiplicity", ";MC tracks;Reco tracks", 100, -0.5, 99.5, 100, -0.5, 99.5);
  
  fEventCount = new TH2F("fEventCount", ";step;event type;count", AliUEHist::fgkCFSteps+2, -2.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
  fEventCount->GetYaxis()->SetBinLabel(1, "ND");
  fEventCount->GetYaxis()->SetBinLabel(2, "SD");
  fEventCount->GetYaxis()->SetBinLabel(3, "DD");
  
  fEventCountDifferential = new TH3F("fEventCountDifferential", ";p_{T,lead};step;event type", 100, 0, 50, AliUEHist::fgkCFSteps, -0.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
  fEventCountDifferential->GetZaxis()->SetBinLabel(1, "ND");
  fEventCountDifferential->GetZaxis()->SetBinLabel(2, "SD");
  fEventCountDifferential->GetZaxis()->SetBinLabel(3, "DD");
  
  fVertexContributors = new TH1F("fVertexContributors", ";contributors;count", 100, -0.5, 99.5);
  
  TH1::AddDirectory(oldStatus);
}

//_____________________________________________________________________________
AliUEHistograms::AliUEHistograms(const AliUEHistograms &c) :
  TObject(),
  fNumberDensitypT(0),
  fSumpT(0),
  fNumberDensityPhi(0),
  fCorrelationpT(0),
  fCorrelationEta(0),
  fCorrelationPhi(0),
  fCorrelationR(0),
  fCorrelationLeading2Phi(0),
  fCorrelationMultiplicity(0),
  fEventCount(0),
  fEventCountDifferential(0),
  fVertexContributors(0)
{
  //
  // AliUEHistograms copy constructor
  //

  ((AliUEHistograms &) c).Copy(*this);
}

//____________________________________________________________________
AliUEHistograms::~AliUEHistograms()
{
  // Destructor
  
  if (fNumberDensitypT)
  {
    delete fNumberDensitypT;
    fNumberDensitypT = 0;
  }
  
  if (fSumpT)
  {
    delete fSumpT;
    fSumpT = 0;
  }
  
  if (fNumberDensityPhi)
  {
    delete fNumberDensityPhi;
    fNumberDensityPhi = 0;
  }
  
  if (fCorrelationpT)
  {
    delete fCorrelationpT;
    fCorrelationpT = 0;
  }
  
  if (fCorrelationEta)
  {
    delete fCorrelationEta;
    fCorrelationEta = 0;
  }
  
  if (fCorrelationPhi)
  {
    delete fCorrelationPhi;
    fCorrelationPhi = 0;
  }
  
  if (fCorrelationR)
  {
    delete fCorrelationR;
    fCorrelationR = 0;
  }

  if (fCorrelationLeading2Phi)
  {
    delete fCorrelationLeading2Phi;
    fCorrelationLeading2Phi = 0;
  }
  
  if (fCorrelationMultiplicity)
  {
    delete fCorrelationMultiplicity;
    fCorrelationMultiplicity = 0;
  }
  
  if (fEventCount)
  {
    delete fEventCount;
    fEventCount = 0;
  }

  if (fEventCountDifferential)
  {
    delete fEventCountDifferential;
    fEventCountDifferential = 0;
  }
  
  if (fVertexContributors)
  {
    delete fVertexContributors;
    fVertexContributors = 0;
  }
}

AliUEHist* AliUEHistograms::GetUEHist(Int_t id)
{
  // returns AliUEHist object, useful for loops
  
  switch (id)
  {
    case 0: return fNumberDensitypT; break;
    case 1: return fSumpT; break;
    case 2: return fNumberDensityPhi; break;
  }
  
  return 0;
}

//____________________________________________________________________
Int_t AliUEHistograms::CountParticles(TList* list, Float_t ptMin)
{
  // counts the number of particles in the list with a pT above ptMin
  // TODO eta cut needed here?
  
  Int_t count = 0;
  for (Int_t j=0; j<list->GetEntries(); j++)
    if (((AliVParticle*) list->At(j))->Pt() > ptMin)
      count++;
      
  return count;
}
  
//____________________________________________________________________
void AliUEHistograms::Fill(Int_t eventType, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max)
{
  // fills the UE event histograms
  //
  // this function needs the leading (track or jet or ...) and four lists of AliVParticles which contain the particles/tracks to be filled in the four regions
  
  // if leading is not set, just fill event statistics
  if (leading)
  {
    Int_t multiplicity = 0;
    
    // TODO configurable?
    Float_t ptMin = 0.15;
    if (leading->Pt() > ptMin)
      multiplicity++;
    
    multiplicity += CountParticles(toward, ptMin);
    multiplicity += CountParticles(away, ptMin);
    multiplicity += CountParticles(min, ptMin);
    multiplicity += CountParticles(max, ptMin);
     
    FillRegion(AliUEHist::kToward, step, leading, toward, multiplicity);
    FillRegion(AliUEHist::kAway,   step, leading, away, multiplicity);
    FillRegion(AliUEHist::kMin,    step, leading, min, multiplicity);
    FillRegion(AliUEHist::kMax,    step, leading, max, multiplicity);
 
    Double_t vars[2];
    vars[0] = leading->Pt();
    vars[1] = multiplicity;
    fNumberDensitypT->GetEventHist()->Fill(vars, step);
    fSumpT->GetEventHist()->Fill(vars, step);
    fNumberDensityPhi->GetEventHist()->Fill(vars, step);
  
    fEventCountDifferential->Fill(leading->Pt(), step, eventType);
  }
  
  FillEvent(eventType, step);
}
  
//____________________________________________________________________
void AliUEHistograms::FillRegion(AliUEHist::Region region, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity)
{
  // loops over AliVParticles in list and fills the given region at the given step
  //
  // See also Fill(...)

  for (Int_t i=0; i<list->GetEntries(); i++)
  {
    AliVParticle* particle = (AliVParticle*) list->At(i);
    
    Double_t vars[5];
    vars[0] = particle->Eta();
    vars[1] = particle->Pt();
    vars[2] = leading->Pt();
    vars[3] = multiplicity;
    vars[4] = leading->Phi() - particle->Phi();
    if (vars[4] > 1.5 * TMath::Pi()) 
      vars[4] -= TMath::TwoPi();
    if (vars[4] < -0.5 * TMath::Pi())
      vars[4] += TMath::TwoPi();

    fNumberDensitypT->GetTrackHist(region)->Fill(vars, step);
    fSumpT->GetTrackHist(region)->Fill(vars, step, particle->Pt());
    
    // fill all in toward region (is anyway as function of delta phi!)
    fNumberDensityPhi->GetTrackHist(AliUEHist::kToward)->Fill(vars, step);
  }
}

//____________________________________________________________________
void AliUEHistograms::Fill(AliVParticle* leadingMC, AliVParticle* leadingReco)
{
  // fills the correlation histograms
  
  if (leadingMC)
  {
    fCorrelationpT->Fill(leadingMC->Pt(), leadingReco->Pt());
    if (leadingMC->Pt() > 0.5)
    {
      fCorrelationEta->Fill(leadingMC->Eta(), leadingReco->Eta());
      fCorrelationPhi->Fill(leadingMC->Phi(), leadingReco->Phi());
    }
    
    Float_t phiDiff = leadingMC->Phi() - leadingReco->Phi();
    if (phiDiff > TMath::Pi())
      phiDiff -= TMath::TwoPi();
    if (phiDiff < -TMath::Pi())
      phiDiff += TMath::TwoPi();
      
    Float_t etaDiff = leadingMC->Eta() - leadingReco->Eta();
    
    fCorrelationR->Fill(TMath::Sqrt(phiDiff * phiDiff + etaDiff * etaDiff), leadingMC->Pt());
    fCorrelationLeading2Phi->Fill(phiDiff, leadingMC->Pt());
  }
  else
  {
    fCorrelationpT->Fill(1.0, leadingReco->Pt());
    if (leadingReco->Pt() > 0.5)
    {
      fCorrelationEta->Fill(0.0, leadingReco->Eta());
      fCorrelationPhi->Fill(0.0, leadingReco->Phi());
    }
  }
}
  
//____________________________________________________________________
void AliUEHistograms::FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, Int_t particleType)
{
  // fills the tracking efficiency objects
  //
  // mc: all primary MC particles
  // recoPrim: reconstructed primaries (again MC particles)
  // recoAll: reconstructed (again MC particles)
  // particleType is: 0 for pion, 1 for kaon, 2 for proton, 3 for others
  
  for (Int_t step=0; step<3; step++)
  {
    TObjArray* list = mc;
    if (step == 1)
      list = recoPrim;
    else if (step == 2)
      list = recoAll;
      
    for (Int_t i=0; i<list->GetEntries(); i++)
    {
      AliVParticle* particle = (AliVParticle*) list->At(i);
      Double_t vars[3];
      vars[0] = particle->Eta();
      vars[1] = particle->Pt();
      vars[2] = particleType;
      
      fNumberDensitypT->GetTrackHistEfficiency()->Fill(vars, step);
      fSumpT->GetTrackHistEfficiency()->Fill(vars, step);
      fNumberDensityPhi->GetTrackHistEfficiency()->Fill(vars, step);
    }
  }
}

//____________________________________________________________________
void AliUEHistograms::FillEvent(Int_t eventType, Int_t step)
{
  // fills the number of events at the given step and the given enty type
  //
  // WARNING: This function is called from Fill, so only call it for steps where Fill is not called
  
  fEventCount->Fill(step, eventType);
}

//____________________________________________________________________
void AliUEHistograms::SetEtaRange(Float_t etaMin, Float_t etaMax)
{
  // sets eta min and max for all contained AliUEHist classes
  
  fNumberDensitypT->SetEtaRange(etaMin, etaMax);
  fSumpT->SetEtaRange(etaMin, etaMax);
  fNumberDensityPhi->SetEtaRange(etaMin, etaMax);
}

//____________________________________________________________________
void AliUEHistograms::SetPtRange(Float_t ptMin, Float_t ptMax)
{
  // sets pT min and max for all contained AliUEHist classes
  
  fNumberDensitypT->SetPtRange(ptMin, ptMax);
  fSumpT->SetPtRange(ptMin, ptMax);
  fNumberDensityPhi->SetPtRange(ptMin, ptMax);
}

//____________________________________________________________________
void AliUEHistograms::SetContaminationEnhancement(TH1F* hist)
{
  // sets the contamination enhancement histogram in all contained AliUEHist classes
  
  fNumberDensitypT->SetContaminationEnhancement(hist);
  fSumpT->SetContaminationEnhancement(hist);
  fNumberDensityPhi->SetContaminationEnhancement(hist);
}  

//____________________________________________________________________
void AliUEHistograms::SetCombineMinMax(Bool_t flag)
{
  // sets pT min and max for all contained AliUEHist classes
  
  fNumberDensitypT->SetCombineMinMax(flag);
  fSumpT->SetCombineMinMax(flag);
  fNumberDensityPhi->SetCombineMinMax(flag);
}

//____________________________________________________________________
void AliUEHistograms::Correct(AliUEHistograms* corrections)
{
  // corrects the contained histograms by calling AliUEHist::Correct
  
  fNumberDensitypT->Correct(corrections->fNumberDensitypT);
  fSumpT->Correct(corrections->fSumpT);
  fNumberDensityPhi->Correct(corrections->fNumberDensityPhi);
}

//____________________________________________________________________
AliUEHistograms &AliUEHistograms::operator=(const AliUEHistograms &c)
{
  // assigment operator

  if (this != &c)
    ((AliUEHistograms &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliUEHistograms::Copy(TObject& c) const
{
  // copy function

  AliUEHistograms& target = (AliUEHistograms &) c;

  if (fNumberDensitypT)
    target.fNumberDensitypT = dynamic_cast<AliUEHist*> (fNumberDensitypT->Clone());

  if (fSumpT)
    target.fSumpT = dynamic_cast<AliUEHist*> (fSumpT->Clone());

  if (fNumberDensityPhi)
    target.fNumberDensityPhi = dynamic_cast<AliUEHist*> (fNumberDensityPhi->Clone());

  if (fCorrelationpT)
    target.fCorrelationpT = dynamic_cast<TH2F*> (fCorrelationpT->Clone());

  if (fCorrelationEta)
    target.fCorrelationEta = dynamic_cast<TH2F*> (fCorrelationEta->Clone());

  if (fCorrelationPhi)
    target.fCorrelationPhi = dynamic_cast<TH2F*> (fCorrelationPhi->Clone());

  if (fCorrelationR)
    target.fCorrelationR = dynamic_cast<TH2F*> (fCorrelationR->Clone());

  if (fCorrelationLeading2Phi)
    target.fCorrelationLeading2Phi = dynamic_cast<TH2F*> (fCorrelationLeading2Phi->Clone());

  if (fCorrelationMultiplicity)
    target.fCorrelationMultiplicity = dynamic_cast<TH2F*> (fCorrelationMultiplicity->Clone());
  
  if (fEventCount)
    target.fEventCount = dynamic_cast<TH2F*> (fEventCount->Clone());

  if (fEventCountDifferential)
    target.fEventCountDifferential = dynamic_cast<TH3F*> (fEventCountDifferential->Clone());
    
  if (fVertexContributors)
    target.fVertexContributors = dynamic_cast<TH1F*> (fVertexContributors->Clone());
}

//____________________________________________________________________
Long64_t AliUEHistograms::Merge(TCollection* list)
{
  // Merge a list of AliUEHistograms objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of objects
  const Int_t kMaxLists = 12;
  TList* lists[kMaxLists];
  
  for (Int_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliUEHistograms* entry = dynamic_cast<AliUEHistograms*> (obj);
    if (entry == 0) 
      continue;

    lists[0]->Add(entry->fNumberDensitypT);
    lists[1]->Add(entry->fSumpT);
    lists[2]->Add(entry->fNumberDensityPhi);
    lists[3]->Add(entry->fCorrelationpT);
    lists[4]->Add(entry->fCorrelationEta);
    lists[5]->Add(entry->fCorrelationPhi);
    lists[6]->Add(entry->fCorrelationR);
    lists[7]->Add(entry->fCorrelationLeading2Phi);
    lists[8]->Add(entry->fCorrelationMultiplicity);
    lists[9]->Add(entry->fEventCount);
    lists[10]->Add(entry->fEventCountDifferential);
    lists[11]->Add(entry->fVertexContributors);

    count++;
  }
    
  fNumberDensitypT->Merge(lists[0]);
  fSumpT->Merge(lists[1]);
  fNumberDensityPhi->Merge(lists[2]);
  fCorrelationpT->Merge(lists[3]);
  fCorrelationEta->Merge(lists[4]);
  fCorrelationPhi->Merge(lists[5]);
  fCorrelationR->Merge(lists[6]);
  fCorrelationLeading2Phi->Merge(lists[7]);
  fCorrelationMultiplicity->Merge(lists[8]);
  fEventCount->Merge(lists[9]);
  fEventCountDifferential->Merge(lists[10]);
  fVertexContributors->Merge(lists[11]);
  
  for (Int_t i=0; i<kMaxLists; i++)
    delete lists[i];

  return count+1;
}

void AliUEHistograms::CopyReconstructedData(AliUEHistograms* from)
{
  // copies those histograms extracted from ESD to this object
  
  fNumberDensitypT->CopyReconstructedData(from->fNumberDensitypT);
  fSumpT->CopyReconstructedData(from->fSumpT);
  fNumberDensityPhi->CopyReconstructedData(from->fNumberDensityPhi);
}
