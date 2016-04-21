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
#include <AliEmcalJet.h>
#include <AliRhoParameter.h>
#include "AliAnalysisTaskParticleRandomizer.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskParticleRandomizer)
/// \endcond

//_____________________________________________________________________________________________________
AliAnalysisTaskParticleRandomizer::AliAnalysisTaskParticleRandomizer() :
  AliAnalysisTaskSE("AliAnalysisTaskParticleRandomizer"), fInitialized(0), fRandomizeInPhi(1), fRandomizeInEta(0), fRandomizeInTheta(0), fRandomizeInPt(0), fMinPhi(0), fMaxPhi(TMath::TwoPi()), fMinEta(-0.9), fMaxEta(+0.9), fMinPt(0), fMaxPt(120), fInputArrayName(), fOutputArrayName(), fInputArray(0), fOutputArray(0), fJetRemovalRhoObj(), fJetRemovalArrayName(), fJetRemovalArray(0), fJetRemovalPtThreshold(999.), fRandom()
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
  if(fInputArrayName == "")
    AliFatal(Form("Name of input array not given!"));
  if(fOutputArrayName == "")
    AliFatal(Form("Name of output array not given!"));

  fRandom = new TRandom3(0);
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::ExecOnce()
{
  // Check if arrays are OK
  fInputArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fInputArrayName.Data())));
  if(!fInputArray)
    AliFatal(Form("Input array '%s' not found!", fInputArrayName.Data()));

  // On demand, load also jets
  if(!fJetRemovalArrayName.IsNull())
  {
    fJetRemovalArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetRemovalArrayName.Data())));
    if(!fJetRemovalArray)
      AliError(Form("Jet array '%s' demanded but not found in event!", fJetRemovalArrayName.Data()));
  }

  if((InputEvent()->FindListObject(Form("%s", fOutputArrayName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fInputArrayName.Data()));

  if(strcmp(fInputArray->GetClass()->GetName(), "AliAODTrack"))
    AliError(Form("Track type %s not yet supported. Use AliAODTrack", fInputArray->GetClass()->GetName()));

  // Copy the input array to the output array
  fOutputArray = new TClonesArray(fInputArray->GetClass()->GetName());
  fOutputArray->SetName(fOutputArrayName.Data());
  InputEvent()->AddObject(fOutputArray);
  fInitialized = kTRUE;
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::UserExec(Option_t *)
{
  // Run once the exec once (must be here to have the input event ready)
  if(!fInitialized)
    ExecOnce();

  Int_t accTracks = 0;
  for(Int_t iPart=0; iPart<fInputArray->GetEntries(); iPart++)
  {
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

    accTracks++;
  }
//  std::cout << Form("%i particles from jets removed out of %i tracks. ", fInputArray->GetEntries()-accTracks, fInputArray->GetEntries()) << std::endl;
}

//_____________________________________________________________________________________________________
Bool_t AliAnalysisTaskParticleRandomizer::IsParticleInJet(Int_t part)
{
  for(Int_t i=0; i<fJetRemovalArray->GetEntries(); i++)
  {
    AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetRemovalArray->At(i));
    Double_t tmpPt = tmpJet->Pt() - tmpJet->Area()*GetExternalRho();

    if(tmpPt >= fJetRemovalPtThreshold)
      if(tmpJet->ContainsTrack(part)>=0)
        return kTRUE;
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
