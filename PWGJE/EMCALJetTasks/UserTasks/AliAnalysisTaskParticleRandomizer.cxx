/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
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

#include <iostream>
#include <TRandom3.h>
#include <AliLog.h>
#include <TString.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <AliAODTrack.h>
#include <AliPicoTrack.h>
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>
#include "TH1D.h"
#include "TH2D.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskParticleRandomizer.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskParticleRandomizer)
/// \endcond

//_____________________________________________________________________________________________________
AliAnalysisTaskParticleRandomizer::AliAnalysisTaskParticleRandomizer() :
  AliAnalysisTaskEmcal("AliAnalysisTaskParticleRandomizer", kFALSE), fRandomizeInPhi(1), fRandomizeInEta(0), fRandomizeInTheta(0), fRandomizeInPt(0), fMinPhi(0), fMaxPhi(TMath::TwoPi()), fMinEta(-0.9), fMaxEta(+0.9), fMinPt(0), fMaxPt(120), fDistributionV2(0), fDistributionV3(0), fDistributionV4(0), fDistributionV5(0), fInputArrayName(), fOutputArrayName(), fInputArray(0), fOutputArray(0), fJetRemovalRhoObj(), fJetRemovalArrayName(), fJetRemovalArray(0), fJetRemovalPtThreshold(999.), fJetRemovalNLeadingJets(0), fJetEmbeddingArrayName(), fJetEmbeddingArray(0), fRandomPsi3(0), fLeadingJet(0), fSubleadingJet(0), fRandom()
{
// constructor
}

//_____________________________________________________________________________________________________
AliAnalysisTaskParticleRandomizer::~AliAnalysisTaskParticleRandomizer()
{
// destructor
}


//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::UserCreateOutputObjects()
{
  // Check user input
  if(fInputArrayName.IsNull())
    AliWarning(Form("Name of input array not given!"));
  if(fOutputArrayName.IsNull())
    AliFatal(Form("Name of output array not given!"));

  fRandom = new TRandom3(0);
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::ExecOnce()
{
  // Check if arrays are OK
  fInputArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fInputArrayName.Data())));
  if(!fInputArrayName.IsNull() && !fInputArray)
    AliFatal(Form("Input array '%s' not found!", fInputArrayName.Data()));

  // On demand, load also jets
  if(!fJetRemovalArrayName.IsNull())
  {
    fJetRemovalArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetRemovalArrayName.Data())));
    if(!fJetRemovalArray)
      AliError(Form("Jet array '%s' demanded but not found in event!", fJetRemovalArrayName.Data()));
  }

  // On demand, load array for embedding
  if(!fJetEmbeddingArrayName.IsNull())
  {
    fJetEmbeddingArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetEmbeddingArrayName.Data())));
    if(!fJetEmbeddingArray)
      AliError(Form("Embedding array '%s' demanded but not found in event!", fJetEmbeddingArrayName.Data()));
  }

  if((InputEvent()->FindListObject(Form("%s", fOutputArrayName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fInputArrayName.Data()));

  if(fInputArray)
    if(strcmp(fInputArray->GetClass()->GetName(), "AliAODTrack"))
      AliError(Form("Track type %s not yet supported. Use AliAODTrack", fInputArray->GetClass()->GetName()));

  // Copy the input array to the output array
  fOutputArray = new TClonesArray("AliAODTrack");
  fOutputArray->SetName(fOutputArrayName.Data());
  InputEvent()->AddObject(fOutputArray);

  AliAnalysisTaskEmcal::ExecOnce();
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskParticleRandomizer::Run()
{
  fRandomPsi3 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi3
  fRandomPsi4 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi4
  fRandomPsi5 = fRandom->Rndm()*TMath::Pi(); // once per event, create a random value dedicated for Psi5

  // Get leading jets on demand
  if(fJetRemovalNLeadingJets)
    GetLeadingJets(fLeadingJet, fSubleadingJet);

  Int_t accTracks = 0;

  // Add events in input array
  if(fInputArray)
    for(Int_t iPart=0; iPart<fInputArray->GetEntries(); iPart++)
    {
      // Remove particles contained in jet array (on demand)
      if(fJetRemovalArray && IsParticleInJet(iPart))
        continue;

      // Take only particles from the randomization acceptance
      AliAODTrack* inputParticle = static_cast<AliAODTrack*>(fInputArray->At(iPart));
      if(fRandomizeInPhi && (inputParticle->Phi() < fMinPhi  || inputParticle->Phi() >= fMaxPhi) )
        continue;
      if( (fRandomizeInTheta || fRandomizeInEta) && (inputParticle->Eta() < fMinEta  || inputParticle->Eta() >= fMaxEta) )
        continue;

      new ((*fOutputArray)[accTracks]) AliAODTrack(*((AliAODTrack*)fInputArray->At(iPart)));

      // Randomize on demand
      AliAODTrack* particle = static_cast<AliAODTrack*>(fOutputArray->At(accTracks));
      RandomizeTrack(particle);

      accTracks++;
    }
  
  // Add particles for embedding (on demand)
  if(fJetEmbeddingArray)
    for(Int_t iPart=0; iPart<fJetEmbeddingArray->GetEntries(); iPart++)
    {
      // Take only particles from the randomization acceptance
      AliPicoTrack* inputParticle = static_cast<AliPicoTrack*>(fJetEmbeddingArray->At(iPart));
      if(fRandomizeInPhi && (inputParticle->Phi() < fMinPhi  || inputParticle->Phi() >= fMaxPhi) )
        continue;
      if( (fRandomizeInTheta || fRandomizeInEta) && (inputParticle->Eta() < fMinEta  || inputParticle->Eta() >= fMaxEta) )
        continue;

      new ((*fOutputArray)[accTracks]) AliAODTrack(*(GetAODTrack(inputParticle)));

      // Randomize on demand
      AliAODTrack* particle = static_cast<AliAODTrack*>(fOutputArray->At(accTracks));
      RandomizeTrack(particle);

      accTracks++;
    }

 //  std::cout << Form("%i particles from jets removed out of %i tracks. ", fInputArray->GetEntries()-accTracks, fInputArray->GetEntries()) << std::endl;
  return kTRUE;
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::RandomizeTrack(AliAODTrack* particle)
{
  if(fRandomizeInPhi)
    particle->SetPhi(fMinPhi + fRandom->Rndm()*(fMaxPhi-fMinPhi));
  if(fRandomizeInTheta)
  {
    Double_t minTheta = 2.*atan(exp(-fMinEta));
    Double_t maxTheta = 2.*atan(exp(-fMaxEta));
    particle->SetTheta(minTheta  + fRandom->Rndm()*(maxTheta-minTheta));
  }
  if(fRandomizeInEta)
  {
    Double_t randomEta = fMinEta  + fRandom->Rndm()*(fMaxEta-fMinEta);
    Double_t randomTheta = 2.*atan(exp(-randomEta));
    particle->SetTheta(randomTheta);
  }

  if(fRandomizeInPt)
    particle->SetPt(fMinPt  + fRandom->Rndm()*(fMaxPt-fMinPt));

  if(fDistributionV2 || fDistributionV3 || fDistributionV4 || fDistributionV5)
    particle->SetPhi(AddFlow(particle->Phi(), particle->Pt()));
}

//_____________________________________________________________________________________________________
AliAODTrack* AliAnalysisTaskParticleRandomizer::GetAODTrack(AliPicoTrack* track)
{
  AliAODTrack* newTrack = new AliAODTrack();
  newTrack->SetPt(track->Pt());
  newTrack->SetTheta(2.*atan(exp(-track->Eta()))); // there is no setter for eta
  newTrack->SetPhi(track->Phi());
  newTrack->SetCharge(track->Charge());
  newTrack->SetLabel(track->GetLabel());

  // Hybrid tracks (compatible with LHC11h)
  UInt_t filterMap = BIT(8) | BIT(9);
  newTrack->SetIsHybridGlobalConstrainedGlobal();
  newTrack->SetFilterMap(filterMap);
  return newTrack;
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskParticleRandomizer::IsParticleInJet(Int_t part)
{
  // Check if particle w/ index 'part' is contained in one of the removal jets
  // The removal jets can be jets with a pt theshold or the leading jets
  for(Int_t i=0; i<fJetRemovalArray->GetEntries(); i++)
  {
    AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetRemovalArray->At(i));

    // Check if to remove leading jets
    if (fJetRemovalNLeadingJets)
    {
      // leading jet removal mode
      if( (fJetRemovalNLeadingJets == 1) && (tmpJet == fLeadingJet))
      {
        if(tmpJet->ContainsTrack(part)>=0)
          return kTRUE;
        else
          return kFALSE; // we know we are done here: the track is not contained in leading jet
      }
      // leading or leading removal mode
      else if ( (fJetRemovalNLeadingJets == 2) && ((tmpJet == fLeadingJet) || (tmpJet == fSubleadingJet)))
      {
        if(tmpJet->ContainsTrack(part)>=0)
          return kTRUE;
      }
    }
    // Check if to remove jets above threshold
    else
    {
      Double_t tmpPt = tmpJet->Pt() - tmpJet->Area()*GetExternalRho();
      if(tmpPt >= fJetRemovalPtThreshold)
        if(tmpJet->ContainsTrack(part)>=0)
          return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________________________________
Double_t AliAnalysisTaskParticleRandomizer::GetExternalRho()
{
 // Get rho from event.
  AliRhoParameter* rho = 0;
  if (!fJetRemovalRhoObj.IsNull()) 
  {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fJetRemovalRhoObj.Data()));
    if (!rho) {
      AliError(Form("%s: Could not retrieve rho with name %s!", GetName(), fJetRemovalRhoObj.Data())); 
      return 0;
    }
  }
  else
    return 0;

  return (rho->GetVal());
}


//_____________________________________________________________________________________________________
Double_t AliAnalysisTaskParticleRandomizer::AddFlow(Double_t phi, Double_t pt)
{
  // adapted from AliFlowTrackSimple
  Double_t precisionPhi = 1e-10;
  Int_t maxNumberOfIterations  = 200;

  Double_t phi0=phi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;
  Int_t ptBin = 0;

  // Evaluate V2 for track pt/centrality
  Double_t v2 = 0;
  if(fDistributionV2)
  {
    ptBin = fDistributionV2->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV2->GetNbinsX())
      v2 = fDistributionV2->GetBinContent(fDistributionV2->GetNbinsX(), fDistributionV2->GetYaxis()->FindBin(fCent));
    else if(ptBin>0)
      v2 = fDistributionV2->GetBinContent(ptBin, fDistributionV2->GetYaxis()->FindBin(fCent));
  }

  // Evaluate V3 for track pt/centrality
  Double_t v3 = 0;
  if(fDistributionV3)
  {
    ptBin = fDistributionV3->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV3->GetNbinsX())
      v3 = fDistributionV3->GetBinContent(fDistributionV3->GetNbinsX(), fDistributionV3->GetYaxis()->FindBin(fCent));
    else if(ptBin>0)
      v3 = fDistributionV3->GetBinContent(ptBin, fDistributionV3->GetYaxis()->FindBin(fCent));
  }

  // Evaluate V4 for track pt/centrality
  Double_t v4 = 0;
  if(fDistributionV4)
  {
    ptBin = fDistributionV4->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV4->GetNbinsX())
      v4 = fDistributionV4->GetBinContent(fDistributionV4->GetNbinsX(), fDistributionV4->GetYaxis()->FindBin(fCent));
    else if(ptBin>0)
      v4 = fDistributionV4->GetBinContent(ptBin, fDistributionV4->GetYaxis()->FindBin(fCent));
  }

  // Evaluate V5 for track pt/centrality
  Double_t v5 = 0;
  if(fDistributionV5)
  {
    ptBin = fDistributionV5->GetXaxis()->FindBin(pt);
    if(ptBin>fDistributionV5->GetNbinsX())
      v5 = fDistributionV5->GetBinContent(fDistributionV5->GetNbinsX(), fDistributionV5->GetYaxis()->FindBin(fCent));
    else if(ptBin>0)
      v5 = fDistributionV5->GetBinContent(ptBin, fDistributionV5->GetYaxis()->FindBin(fCent));
  }

  // Add all v's
  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=phi; //store last value for comparison
    f =  phi-phi0
        +      v2*TMath::Sin(2.*(phi-(fEPV0+(TMath::Pi()/2.))))
        +2./3.*v3*TMath::Sin(3.*(phi-fRandomPsi3))
        +0.5  *v4*TMath::Sin(4.*(phi-fRandomPsi4))
        +0.4  *v5*TMath::Sin(5.*(phi-fRandomPsi5));

    fp =  1.0+2.0*(
           +v2*TMath::Cos(2.*(phi-(fEPV0+(TMath::Pi()/2.))))
           +v3*TMath::Cos(3.*(phi-fRandomPsi3))
           +v4*TMath::Cos(4.*(phi-fRandomPsi4))
           +v5*TMath::Cos(5.*(phi-fRandomPsi5))); //first derivative

    phi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,phi,precisionPhi)) break;
  }

  return phi;
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::GetLeadingJets(AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading)
{
  // Customized from AliJetContainer::GetLeadingJet()
  // Get the leading+subleading jet; the sorting is according to pt-A*rho

  jetLeading = 0;
  jetSubLeading = 0;

  Double_t     tmpLeadingPt = 0;
  Double_t     tmpSubleadingPt = 0;

  for(Int_t i=0; i<fJetRemovalArray->GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetRemovalArray->At(i));
    if      ( (jet->Pt()-jet->Area()*GetExternalRho()) > tmpLeadingPt )
    {
      jetSubLeading = jetLeading;
      jetLeading = jet;
      tmpSubleadingPt = tmpLeadingPt;
      tmpLeadingPt = jet->Pt()-jet->Area()*GetExternalRho();
    }
    else if ( (jet->Pt()-jet->Area()*GetExternalRho()) > tmpSubleadingPt )
    {
      jetSubLeading = jet;
      tmpSubleadingPt = jet->Pt()-jet->Area()*GetExternalRho();
    }
  }
}
