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
#include "AliAODTrack.h"

#include "TList.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"

ClassImp(AliUEHistograms)

const Int_t AliUEHistograms::fgkUEHists = 3;

AliUEHistograms::AliUEHistograms(const char* name, const char* histograms) : 
  TNamed(name, name),
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
  fVertexContributors(0),
  fCentralityDistribution(0),
  fCentralityCorrelation(0),
  fITSClusterMap(0),
  fSelectCharge(0),
  fRunNumber(0)
{
  // Constructor
  //
  // the string histograms defines which histograms are created:
  //    1 = NumberDensitypT
  //    2 = SumpT
  //    3 = NumberDensityPhi
  //    4 = NumberDensityPhiCentrality (other multiplicity for Pb)
  
  fTwoTrackDistancePt[0] = 0;
  fTwoTrackDistancePt[1] = 0;
  
  TString histogramsStr(histograms);
  
  if (histogramsStr.Contains("1"))
    fNumberDensitypT = new AliUEHist("NumberDensitypT");
  if (histogramsStr.Contains("2"))
    fSumpT = new AliUEHist("SumpT");
  
  if (histogramsStr.Contains("3"))
    fNumberDensityPhi = new AliUEHist("NumberDensityPhi");
  else if (histogramsStr.Contains("4"))
    fNumberDensityPhi = new AliUEHist("NumberDensityPhiCentrality");
  else if (histogramsStr.Contains("5"))
    fNumberDensityPhi = new AliUEHist("NumberDensityPhiCentralityVtx");
  
  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  if (!histogramsStr.Contains("4") && !histogramsStr.Contains("5"))
  {
    fCorrelationpT  = new TH2F("fCorrelationpT", ";p_{T,lead} (MC);p_{T,lead} (RECO)", 200, 0, 50, 200, 0, 50);
    fCorrelationEta = new TH2F("fCorrelationEta", ";#eta_{lead} (MC);#eta_{T,lead} (RECO)", 200, -1, 1, 200, -1, 1);
    fCorrelationPhi = new TH2F("fCorrelationPhi", ";#varphi_{lead} (MC);#varphi_{T,lead} (RECO)", 200, 0, TMath::TwoPi(), 200, 0, TMath::TwoPi());
  }
  else
  {
    fCorrelationpT  = new TH2F("fCorrelationpT", ";Centrality;p_{T} (RECO)", 100, 0, 100.001, 200, 0, 50);
    fCorrelationEta = new TH2F("fCorrelationEta", ";Centrality;#eta (RECO)", 100, 0, 100.001, 200, -1, 1);
    fCorrelationPhi = new TH2F("fCorrelationPhi", ";Centrality;#varphi (RECO)", 100, 0, 100.001, 200, 0, TMath::TwoPi());
  }
  
  fCorrelationR =   new TH2F("fCorrelationR", ";R;p_{T,lead} (MC)", 200, 0, 2, 200, 0, 50);
  fCorrelationLeading2Phi = new TH2F("fCorrelationLeading2Phi", ";#Delta #varphi;p_{T,lead} (MC)", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 50);
  fCorrelationMultiplicity = new TH2F("fCorrelationMultiplicity", ";MC tracks;Reco tracks", 100, -0.5, 99.5, 100, -0.5, 99.5);
  
  if (!histogramsStr.Contains("4") && !histogramsStr.Contains("5"))
  {
    fEventCount = new TH2F("fEventCount", ";step;event type;count", AliUEHist::fgkCFSteps+2, -2.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
    fEventCount->GetYaxis()->SetBinLabel(1, "ND");
    fEventCount->GetYaxis()->SetBinLabel(2, "SD");
    fEventCount->GetYaxis()->SetBinLabel(3, "DD");
  }
  else
  {
    fEventCount = new TH2F("fEventCount", ";step;centrality;count", AliUEHist::fgkCFSteps+2, -2.5, -0.5 + AliUEHist::fgkCFSteps, fNumberDensityPhi->GetEventHist()->GetNBins(1), fNumberDensityPhi->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray());
  }
  
  fEventCountDifferential = new TH3F("fEventCountDifferential", ";p_{T,lead};step;event type", 100, 0, 50, AliUEHist::fgkCFSteps, -0.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
  fEventCountDifferential->GetZaxis()->SetBinLabel(1, "ND");
  fEventCountDifferential->GetZaxis()->SetBinLabel(2, "SD");
  fEventCountDifferential->GetZaxis()->SetBinLabel(3, "DD");
  
  fVertexContributors = new TH1F("fVertexContributors", ";contributors;count", 100, -0.5, 99.5);
  
  if (fNumberDensityPhi)
  {
    fCentralityDistribution = new TH1F("fCentralityDistribution", ";centrality;count", fNumberDensityPhi->GetEventHist()->GetNBins(1), fNumberDensityPhi->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray());
    fCentralityCorrelation = new TH2F("fCentralityCorrelation", ";centrality;multiplicity", fNumberDensityPhi->GetEventHist()->GetNBins(1), fNumberDensityPhi->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray(), 200, 0, 4000);
  }
  
  fITSClusterMap = new TH3F("fITSClusterMap", "; its cluster map; centrality; pT", 256, -0.5, 255.5, 20, 0, 100.001, 100, 0, 20);
  
  TH1::AddDirectory(oldStatus);
}

//_____________________________________________________________________________
AliUEHistograms::AliUEHistograms(const AliUEHistograms &c) :
  TNamed(fName, fTitle),
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
  fVertexContributors(0),
  fCentralityDistribution(0),
  fCentralityCorrelation(0),
  fITSClusterMap(0),
  fSelectCharge(0),
  fRunNumber(0)
{
  //
  // AliUEHistograms copy constructor
  //

  fTwoTrackDistancePt[0] = 0;
  fTwoTrackDistancePt[1] = 0;

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
  
  if (fCentralityDistribution)
  {
    delete fCentralityDistribution;
    fCentralityDistribution = 0;
  }
  
  if (fCentralityCorrelation)
  {
    delete fCentralityCorrelation;
    fCentralityCorrelation = 0;
  }
  
  if (fITSClusterMap)
  {
    delete fITSClusterMap;
    fITSClusterMap = 0;
  }

  for (Int_t i=0; i<2; i++)
    if (fTwoTrackDistancePt[i])
    {
      delete fTwoTrackDistancePt[i];
      fTwoTrackDistancePt[i] = 0;
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
void AliUEHistograms::Fill(Int_t eventType, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max)
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
     
    FillRegion(AliUEHist::kToward, zVtx, step, leading, toward, multiplicity);
    FillRegion(AliUEHist::kAway,   zVtx, step, leading, away, multiplicity);
    FillRegion(AliUEHist::kMin,    zVtx, step, leading, min, multiplicity);
    FillRegion(AliUEHist::kMax,    zVtx, step, leading, max, multiplicity);
 
    Double_t vars[3];
    vars[0] = leading->Pt();
    vars[1] = multiplicity;
    vars[2] = zVtx;
    for (Int_t i=0; i<fgkUEHists; i++)
      if (GetUEHist(i))
        GetUEHist(i)->GetEventHist()->Fill(vars, step);
  
    fEventCountDifferential->Fill(leading->Pt(), step, eventType);
  }
  
  FillEvent(eventType, step);
}
  
//____________________________________________________________________
void AliUEHistograms::FillRegion(AliUEHist::Region region, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity)
{
  // loops over AliVParticles in list and fills the given region at the given step
  //
  // See also Fill(...)

  for (Int_t i=0; i<list->GetEntries(); i++)
  {
    AliVParticle* particle = (AliVParticle*) list->At(i);
    
    Double_t vars[6];
    vars[0] = particle->Eta();
    vars[1] = particle->Pt();
    vars[2] = leading->Pt();
    vars[3] = multiplicity;
    vars[4] = leading->Phi() - particle->Phi();
    if (vars[4] > 1.5 * TMath::Pi()) 
      vars[4] -= TMath::TwoPi();
    if (vars[4] < -0.5 * TMath::Pi())
      vars[4] += TMath::TwoPi();
    vars[5] = zVtx;
    
    if (fNumberDensitypT)
      fNumberDensitypT->GetTrackHist(region)->Fill(vars, step);
      
    if (fSumpT)
      fSumpT->GetTrackHist(region)->Fill(vars, step, particle->Pt());
    
    // fill all in toward region (is anyway as function of delta phi!)
    if (fNumberDensityPhi)
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
void AliUEHistograms::FillCorrelations(Double_t centrality, Float_t zVtx, AliUEHist::CFStep step, TObjArray* particles, TObjArray* mixed, Float_t weight, Bool_t firstTime)
{
  // fills the fNumberDensityPhi histogram
  //
  // this function need a list of AliVParticles which contain the particles/tracks to be filled
  //
  // if mixed is non-0, mixed events are filled, the trigger particle is from particles, the associated from mixed
  // if weight < 0, then the pt of the associated particle is filled as weight
  
  Bool_t fillpT = kFALSE;
  if (weight < 0)
    fillpT = kTRUE;
  
  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* input = (mixed) ? mixed : particles;
  TArrayF eta(input->GetEntriesFast());
  for (Int_t i=0; i<input->GetEntriesFast(); i++)
    eta[i] = ((AliVParticle*) input->At(i))->Eta();
  
  // if particles is not set, just fill event statistics
  if (particles)
  {
    Int_t jMax = particles->GetEntriesFast();
    if (mixed)
      jMax = mixed->GetEntriesFast();
    
    for (Int_t i=0; i<particles->GetEntriesFast(); i++)
    {
      AliVParticle* triggerParticle = (AliVParticle*) particles->At(i);
      
      // some optimization
      Float_t triggerEta = triggerParticle->Eta();
        
      if (!mixed)
      {
        // QA
        fCorrelationpT->Fill(centrality, triggerParticle->Pt());
        fCorrelationEta->Fill(centrality, triggerEta);
        fCorrelationPhi->Fill(centrality, triggerParticle->Phi());
/*        if (dynamic_cast<AliAODTrack*>(triggerParticle))
          fITSClusterMap->Fill(((AliAODTrack*) triggerParticle)->GetITSClusterMap(), centrality, triggerParticle->Pt());*/
      }
        
      for (Int_t j=0; j<jMax; j++)
      {
        if (!mixed && i == j)
          continue;
      
        AliVParticle* particle = 0;
        if (!mixed)
          particle = (AliVParticle*) particles->At(j);
        else
          particle = (AliVParticle*) mixed->At(j);
        
        // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event for cross-checks)
        if (mixed && triggerParticle == particle)
          continue;
        
        if (particle->Pt() > triggerParticle->Pt())
          continue;
          
        if (fSelectCharge > 0)
        {
          // skip like sign
          if (fSelectCharge == 1 && particle->Charge() * triggerParticle->Charge() > 0)
            continue;
            
          // skip unlike sign
          if (fSelectCharge == 2 && particle->Charge() * triggerParticle->Charge() < 0)
            continue;
        }
        
        Double_t vars[6];
        vars[0] = triggerEta - eta[j];
        vars[1] = particle->Pt();
        vars[2] = triggerParticle->Pt();
        vars[3] = centrality;
        vars[4] = triggerParticle->Phi() - particle->Phi();
        if (vars[4] > 1.5 * TMath::Pi()) 
          vars[4] -= TMath::TwoPi();
        if (vars[4] < -0.5 * TMath::Pi())
          vars[4] += TMath::TwoPi();
	vars[5] = zVtx;
	
	if (fillpT)
	  weight = particle->Pt();
    
        // fill all in toward region and do not use the other regions
        fNumberDensityPhi->GetTrackHist(AliUEHist::kToward)->Fill(vars, step, weight);
      }
 
      if (firstTime)
      {
        // once per trigger particle
        Double_t vars[3];
        vars[0] = triggerParticle->Pt();
        vars[1] = centrality;
	vars[2] = zVtx;
        fNumberDensityPhi->GetEventHist()->Fill(vars, step);
      }
    }
  }
  
  fCentralityDistribution->Fill(centrality);
  fCentralityCorrelation->Fill(centrality, particles->GetEntriesFast());
  FillEvent(centrality, step);
}
  
//____________________________________________________________________
void AliUEHistograms::FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, Int_t particleType, Double_t centrality)
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
      
    for (Int_t i=0; i<list->GetEntriesFast(); i++)
    {
      AliVParticle* particle = (AliVParticle*) list->At(i);
      Double_t vars[4];
      vars[0] = particle->Eta();
      vars[1] = particle->Pt();
      vars[2] = particleType;
      vars[3] = centrality;
      
      for (Int_t j=0; j<fgkUEHists; j++)
        if (GetUEHist(j))
          GetUEHist(j)->GetTrackHistEfficiency()->Fill(vars, step);
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
void AliUEHistograms::FillEvent(Double_t centrality, Int_t step)
{
  // fills the number of events at the given step and the given centrality
  //
  // WARNING: This function is called from Fill, so only call it for steps where Fill is not called
  
  fEventCount->Fill(step, centrality);
}

//____________________________________________________________________
void AliUEHistograms::SetEtaRange(Float_t etaMin, Float_t etaMax)
{
  // sets eta min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetEtaRange(etaMin, etaMax);
}

//____________________________________________________________________
void AliUEHistograms::SetPtRange(Float_t ptMin, Float_t ptMax)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetPtRange(ptMin, ptMax);
}

//____________________________________________________________________
void AliUEHistograms::SetZVtxRange(Float_t min, Float_t max)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetZVtxRange(min, max);
}

//____________________________________________________________________
void AliUEHistograms::SetContaminationEnhancement(TH1F* hist)
{
  // sets the contamination enhancement histogram in all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetContaminationEnhancement(hist);
}  

//____________________________________________________________________
void AliUEHistograms::SetCombineMinMax(Bool_t flag)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetCombineMinMax(flag);
}

//____________________________________________________________________
void AliUEHistograms::Correct(AliUEHistograms* corrections)
{
  // corrects the contained histograms by calling AliUEHist::Correct
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Correct(corrections->GetUEHist(i));
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

  if (fCentralityDistribution)
    target.fCentralityDistribution = dynamic_cast<TH1F*> (fCentralityDistribution->Clone());
    
  if (fCentralityCorrelation)
    target.fCentralityCorrelation = dynamic_cast<TH2F*> (fCentralityCorrelation->Clone());

  if (fITSClusterMap)
    target.fITSClusterMap = dynamic_cast<TH3F*> (fITSClusterMap->Clone());
    
  for (Int_t i=0; i<2; i++)
    if (fTwoTrackDistancePt[i])
      target.fTwoTrackDistancePt[i] = dynamic_cast<TH3F*> (fTwoTrackDistancePt[i]->Clone());

  target.fSelectCharge = fSelectCharge;
  target.fRunNumber = fRunNumber;
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
  const Int_t kMaxLists = 17;
  TList* lists[kMaxLists];
  
  for (Int_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliUEHistograms* entry = dynamic_cast<AliUEHistograms*> (obj);
    if (entry == 0) 
      continue;

    if (entry->fNumberDensitypT)
      lists[0]->Add(entry->fNumberDensitypT);
    if (entry->fSumpT)
      lists[1]->Add(entry->fSumpT);
    if (entry->fNumberDensityPhi)
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
    lists[12]->Add(entry->fCentralityDistribution);
    lists[13]->Add(entry->fITSClusterMap);
    if (fTwoTrackDistancePt[0])
      lists[14]->Add(entry->fTwoTrackDistancePt[0]);
    if (fTwoTrackDistancePt[1])
      lists[15]->Add(entry->fTwoTrackDistancePt[1]);
    if (fCentralityCorrelation)
      lists[16]->Add(entry->fCentralityCorrelation);
    count++;
  }
    
  if (fNumberDensitypT)
    fNumberDensitypT->Merge(lists[0]);
  if (fSumpT)
    fSumpT->Merge(lists[1]);
  if (fNumberDensityPhi)
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
  fCentralityDistribution->Merge(lists[12]);
  fITSClusterMap->Merge(lists[13]);
  if (fTwoTrackDistancePt[0])
    fTwoTrackDistancePt[0]->Merge(lists[14]);
  if (fTwoTrackDistancePt[1])
    fTwoTrackDistancePt[1]->Merge(lists[15]);
  if (fCentralityCorrelation)
    fCentralityCorrelation->Merge(lists[16]);
  
  for (Int_t i=0; i<kMaxLists; i++)
    delete lists[i];

  return count+1;
}

void AliUEHistograms::CopyReconstructedData(AliUEHistograms* from)
{
  // copies those histograms extracted from ESD to this object
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->CopyReconstructedData(from->GetUEHist(i));
}

void AliUEHistograms::ExtendTrackingEfficiency(Bool_t verbose)
{
  // delegates to AliUEHists

  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->ExtendTrackingEfficiency(verbose);
}

void AliUEHistograms::Scale(Double_t factor)
{
  // scales all contained histograms by the given factor
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Scale(factor);
      
  TList list;
  list.Add(fCorrelationpT);
  list.Add(fCorrelationEta);
  list.Add(fCorrelationPhi);
  list.Add(fCorrelationR);
  list.Add(fCorrelationLeading2Phi);
  list.Add(fCorrelationMultiplicity);
  list.Add(fEventCount);
  list.Add(fEventCountDifferential);
  list.Add(fVertexContributors);
  list.Add(fCentralityDistribution);
  list.Add(fCentralityCorrelation);
  list.Add(fITSClusterMap);
  list.Add(fTwoTrackDistancePt[0]);
  list.Add(fTwoTrackDistancePt[1]);
  
  for (Int_t i=0; i<list.GetEntries(); i++)
    ((TH1*) list.At(i))->Scale(factor);
}

void AliUEHistograms::Reset()
{
  // delegates to AliUEHists

  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Reset();
}

void AliUEHistograms::TwoTrackEfficiency(TObjArray* tracks, TObjArray* mixed, Float_t bSign)
{
  // takes the input list <tracks> and fills histograms to study two two-track efficiency effects
  // fTwoTrackDistancePt[i] (i = 0 same, i = 1 mixed)
  //
  // the variables have been developed by the HBT group 
  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700
  
  if (!fTwoTrackDistancePt[0])
  {
    fTwoTrackDistancePt[0] = new TH3F("fTwoTrackDistancePt[0]", ";#Delta#eta;#Delta#varphi^{*}_{min};#Delta p_{T}", 100, -0.05, 0.05, 400, -0.2, 0.2, 20, 0, 10);
    fTwoTrackDistancePt[1] = (TH3F*) fTwoTrackDistancePt[0]->Clone("fTwoTrackDistancePt[1]");
  }

  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TArrayF eta1(tracks->GetEntriesFast());
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
    eta1[i] = ((AliVParticle*) tracks->At(i))->Eta();
  
  Int_t jMax = tracks->GetEntriesFast();
  if (mixed)
    jMax = mixed->GetEntriesFast();
  
  TArrayF eta2(jMax);
  if (!mixed)
    eta2 = eta1;
  else
    for (Int_t i=0; i<mixed->GetEntriesFast(); i++)
      eta2[i] = ((AliVParticle*) mixed->At(i))->Eta();

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle1 = (AliVParticle*) tracks->At(i);
    Float_t phi1 = particle1->Phi();
    Float_t pt1 = particle1->Pt();
    Float_t charge1 = particle1->Charge();
    
    for (Int_t j=0; j<jMax; j++)
    {
      if (!mixed && i == j)
	continue;
      
      AliVParticle* particle2 = 0;
      if (mixed)
	particle2 = (AliVParticle*) mixed->At(j);
      else
	particle2 = (AliVParticle*) tracks->At(j);
      Float_t phi2 = particle2->Phi();
      Float_t pt2 = particle2->Pt();
      Float_t charge2 = particle2->Charge();
      
      if (pt2 > pt1)
	continue;
      
//       Double_t dpt = TMath::Abs(pt1 - pt2);
      Float_t deta = eta1[i] - eta2[j];
      Float_t detaabs = TMath::Abs(deta);
      
      // optimization
      if (detaabs > 0.05 && (pt1 < 8 || pt1 > 15))
	continue;
      
      Float_t dphistarmin = 1e5;
      Float_t dphistarminabs = 1e5;

      for (Double_t rad=0.8; rad<2.51; rad+=0.01) 
      {
	Float_t dphistar = phi1 - phi2 - TMath::ASin(charge1 * 0.075 * bSign * rad / pt1) + TMath::ASin(charge2 * 0.075 * bSign * rad / pt2);
	Float_t dphistarabs = TMath::Abs(dphistar);
	
	if (dphistarabs < dphistarminabs)
	{
	  dphistarmin = dphistar;
	  dphistarminabs = dphistarabs;
	}
      }

      Float_t fillPt = pt2;
      
      // analyze region for IAA paper
      if (pt1 < 8 || pt1 > 15)
	fillPt = 0.25;
    
      if (!mixed)
	fTwoTrackDistancePt[0]->Fill(deta, dphistarmin, fillPt);
      else
	fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, fillPt);
    }
  }
}
