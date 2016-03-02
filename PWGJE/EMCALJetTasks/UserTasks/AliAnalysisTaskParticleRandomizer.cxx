/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>

#include <TRandom3.h>
#include <AliLog.h>
#include <TString.h>
#include <TMath.h>
#include <TClonesArray.h>
#include "AliAnalysisTaskParticleRandomizer.h"
#include <AliAODTrack.h>


ClassImp(AliAnalysisTaskParticleRandomizer)

//_____________________________________________________________________________________________________
AliAnalysisTaskParticleRandomizer::AliAnalysisTaskParticleRandomizer() :
  AliAnalysisTaskSE("AliAnalysisTaskParticleRandomizer"), fInitialized(0), fRandomizeInPhi(1), fRandomizeInEta(0), fRandomizeInPt(0), fMinPhi(0), fMaxPhi(TMath::TwoPi()), fMinEta(-0.9), fMaxEta(+0.9), fMinPt(0), fMaxPt(120), fInputArrayName(), fOutputArrayName(), fInputArray(0), fOutputArray(0), fRandom()
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

  fRandom = new TRandom3(0); //TODO: Use different seed
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskParticleRandomizer::ExecOnce()
{
  // Check if arrays are OK
  fInputArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fInputArrayName.Data())));
  if(!fInputArray)
    AliFatal(Form("Input array '%s' not found!", fInputArrayName.Data()));

  if((InputEvent()->FindListObject(Form("%s", fOutputArrayName.Data()))))
    AliFatal(Form("Output array '%s' already exists in the event! Rename it.", fInputArrayName.Data()));

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

  for(Int_t iPart=0; iPart<fInputArray->GetEntries(); iPart++)
  {
    // Create the particle according to the used type and write it into the output clones array
    if(!strcmp(fOutputArray->GetClass()->GetName(), "AliAODTrack"))
    {
      new ((*fOutputArray)[iPart]) AliAODTrack(*((AliAODTrack*)fInputArray->At(iPart)));

      // Randomize on demand
      AliAODTrack* particle = static_cast<AliAODTrack*>(fOutputArray->At(iPart));

      if(fRandomizeInPhi)
        particle->SetPhi(fMinPhi + fRandom->Rndm()*(fMaxPhi-fMinPhi));
      if(fRandomizeInEta)
      {
        Double_t minTheta = 2.*atan(exp(-fMinEta));
        Double_t maxTheta = 2.*atan(exp(-fMaxEta));;
        particle->SetTheta(minTheta  + fRandom->Rndm()*(maxTheta-minTheta));
      }
      if(fRandomizeInPt)
        particle->SetPt(fMinPt  + fRandom->Rndm()*(fMaxPt-fMinPt));
    }
    else
      AliFatal(Form("Track type %s not yet supported.", fOutputArray->GetClass()->GetName()));
  }
}
