/**************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include "AliEmcalMCTrackSelector.h"

#include <iostream>

#include <TClonesArray.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliNamedArrayI.h"

#include "AliLog.h"

ClassImp(AliEmcalMCTrackSelector)

AliEmcalMCTrackSelector::AliEmcalMCTrackSelector() : 
  AliAnalysisTaskSE("AliEmcalMCTrackSelector"),
  fParticlesOutName("MCParticlesSelected"),
  fOnlyPhysPrim(kTRUE),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fOnlyHIJING(kFALSE),
  fRejectPhotonMothers(false),
  fEtaMax(1),
  fParticlesMapName(""),
  fInit(kFALSE),
  fParticlesIn(0),
  fParticlesOut(0),
  fParticlesMap(0),
  fEvent(0),
  fMC(0),
  fIsESD(kFALSE),
  fDisabled(kFALSE)
{
}

AliEmcalMCTrackSelector::AliEmcalMCTrackSelector(const char *name) : 
  AliAnalysisTaskSE(name),
  fParticlesOutName("MCParticlesSelected"),
  fOnlyPhysPrim(kTRUE),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fOnlyHIJING(kFALSE),
  fRejectPhotonMothers(false),
  fEtaMax(1),
  fParticlesMapName(""),
  fInit(kFALSE),
  fParticlesIn(0),
  fParticlesOut(0),
  fParticlesMap(0),
  fEvent(0),
  fMC(0),
  fIsESD(kFALSE),
  fDisabled(kFALSE)
{
}

void AliEmcalMCTrackSelector::UserExec(Option_t *) 
{
  if (fDisabled) return;

  if (!fInit) {
    fEvent = InputEvent();
    if (!fEvent) {
      AliErrorStream() << "Could not retrieve event! Returning" << std::endl;
      return;
    }

    if (fEvent->InheritsFrom("AliESDEvent")) fIsESD = kTRUE;
    else fIsESD = kFALSE;

    TObject *obj = fEvent->FindListObject(fParticlesOutName);
    if (obj) { // the output array is already present in the array!
      AliErrorStream() << "The output array "  << fParticlesOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }
    else {  // copy the array from the standard ESD/AOD collections, and filter if requested      

      fParticlesOut = new TClonesArray("AliAODMCParticle");  // the output will always be of AliAODMCParticle, regardless of the input
      fParticlesOut->SetName(fParticlesOutName);
      fEvent->AddObject(fParticlesOut);

      fParticlesMapName = fParticlesOutName;
      fParticlesMapName += "_Map";

      if (fEvent->FindListObject(fParticlesMapName)) {
        AliErrorStream() << "The output array map " << fParticlesMapName << " is already present in the event! Task will be disabled." << std::endl;
        fDisabled = kTRUE;
        return;
      }
      else {
        fParticlesMap = new AliNamedArrayI(fParticlesMapName, 99999);
        fEvent->AddObject(fParticlesMap);
      }

      if (!fIsESD) {
        fParticlesIn = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!fParticlesIn) {
          AliErrorStream() << "Could not retrieve AOD MC particles! Task will be disabled." << std::endl;
          fDisabled = kTRUE;
          return;
        }
        TClass *cl = fParticlesIn->GetClass();
        if (!cl->GetBaseClass("AliAODMCParticle")) {
          AliErrorStream() << GetName() << ": Collection  " << AliAODMCParticle::StdBranchName() << " %s does not contain AliAODMCParticle! Task will be disabled." << std::endl;
          fDisabled = kTRUE;
          fParticlesIn = 0;
          return;
        }
      }
    }

    fMC = MCEvent();
    if (!fMC) {
      AliErrorStream() << "Could not retrieve MC event! Returning" << std::endl;
      fDisabled = kTRUE;
      return;
    }

    fInit = kTRUE;
  }

  if (fIsESD) ConvertMCParticles(fMC, fParticlesOut, fParticlesMap);
  else CopyMCParticles(fParticlesIn, fParticlesOut, fParticlesMap);
}

void AliEmcalMCTrackSelector::ConvertMCParticles(AliMCEvent* mcEvent, TClonesArray* partOut, AliNamedArrayI* partMap)
{

  // clear container (normally a null operation as the event should clean it already)
  partOut->Delete();

  const Int_t Nparticles = mcEvent->GetNumberOfTracks();

  if (partMap) {
    // clear particles map
    partMap->Clear();
    if (partMap->GetSize() <= Nparticles) partMap->Set(Nparticles*2);
  }

  const Int_t nprim = mcEvent->GetNumberOfPrimaries();

  AliDebugStream(3) << "Number of particles: " << Nparticles << std::endl;

  // loop over particles
  for (Int_t iPart = 0, nacc = 0; iPart < Nparticles; iPart++) {

    if (partMap) partMap->AddAt(-1, iPart);

    AliMCParticle* part = static_cast<AliMCParticle*>(mcEvent->GetTrack(iPart));

    if (!part) continue;

    Bool_t isPhysPrim = mcEvent->IsPhysicalPrimary(iPart);

    Int_t flag = 0;
    if (iPart < nprim) flag |= AliAODMCParticle::kPrimary;
    if (isPhysPrim) flag |= AliAODMCParticle::kPhysicalPrim;
    if (mcEvent->IsSecondaryFromWeakDecay(iPart)) flag |= AliAODMCParticle::kSecondaryFromWeakDecay;
    if (mcEvent->IsSecondaryFromMaterial(iPart)) flag |= AliAODMCParticle::kSecondaryFromMaterial;

    // Reject particle if is a photon and mother of another photon
    if (fRejectPhotonMothers) {
      if(TMath::Abs(part->PdgCode()) == 22){
        bool hasPhotonDaughter = false;
        // check if particle has photon daughter
        if(part->GetDaughterFirst() > -1 && part->GetDaughterLast() > -1) {
          for(int idaughter = part->GetDaughterFirst(); idaughter <= part->GetDaughterLast(); idaughter++) {
            AliMCParticle *daughter = static_cast<AliMCParticle *>(mcEvent->GetTrack(idaughter));
            if(TMath::Abs(daughter->PdgCode()) == 22) {
              hasPhotonDaughter = true;
              break;
            }
          }
        }
        if(hasPhotonDaughter){
          AliDebugStream(2) << "Found photon which is the mother of another photon, not selecting ..." << std::endl; 
          continue;
        } 
      }
    }

    AliDebugStream(3) << "Particle " << iPart << ": pt = " << part->Pt() << ", PDG = " << part->PdgCode() <<
        ", mother " << part->GetMother() <<
        ", kPrimary? " <<  Bool_t((flag & AliAODMCParticle::kPrimary) != 0) <<
        ", kPhysicalPrim? " <<  Bool_t((flag & AliAODMCParticle::kPhysicalPrim) != 0) <<
        ", kSecondaryFromWeakDecay? " <<  Bool_t((flag & AliAODMCParticle::kSecondaryFromWeakDecay) != 0) <<
        ", kSecondaryFromMaterial? " << Bool_t((flag & AliAODMCParticle::kSecondaryFromMaterial) != 0) <<
        ", nacc = " << nacc <<
        ", iPart = " << iPart <<
        std::endl;

    // Do not put rejected particles on the output container
    AliAODMCParticle parttmp(part, iPart, flag);
    parttmp.SetGeneratorIndex(part->GetGeneratorIndex());    
    parttmp.SetStatus(part->Particle()->GetStatusCode());
    parttmp.SetMCProcessCode(part->Particle()->GetUniqueID());

    if (!AcceptParticle(&parttmp)) continue;    

    new ((*partOut)[nacc]) AliAODMCParticle(parttmp);
    if (partMap) partMap->AddAt(nacc, iPart);
    nacc++;
  }
}

void AliEmcalMCTrackSelector::CopyMCParticles(TClonesArray* partIn, TClonesArray* partOut, AliNamedArrayI* partMap)
{

  if (!partIn) return;
  AliDebugStream(5) << "Particles in: " << partIn->GetName() << ", out: " << partOut->GetName() << std::endl;

  // clear container (normally a null operation as the event should clean it already)
  partOut->Delete();

  const Int_t Nparticles = partIn->GetEntriesFast();

  if (partMap) {
    // clear particles map
    partMap->Clear();
    if (partMap->GetSize() <= Nparticles) partMap->Set(Nparticles*2);
  }

  AliDebugStream(2) << "Total number of particles = " << Nparticles << std::endl;

  Int_t nacc = 0;

  // loop over particles
  for (Int_t iPart = 0; iPart < Nparticles; iPart++) {
    if (partMap) partMap->AddAt(-1, iPart);

    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(partIn->At(iPart));

    if (!AcceptParticle(part)) continue;

    // Reject particle if is a photon and mother of another photon
    if (fRejectPhotonMothers) {
      if(TMath::Abs(part->PdgCode()) == 22){
        bool hasPhotonDaughter = false;
        // check if particle has photon daughter
        if(part->GetDaughterFirst() > -1 && part->GetDaughterLast() > -1) {
          for(int idaughter = part->GetDaughterFirst(); idaughter <= part->GetDaughterLast(); idaughter++) {
            AliAODMCParticle *daughter = static_cast<AliAODMCParticle *>(partIn->At(idaughter));
            if(TMath::Abs(daughter->PdgCode()) == 22) {
              hasPhotonDaughter = true;
              break;
            }
          }
        }
        if(hasPhotonDaughter){
          AliDebugStream(2) << "Found photon which is the mother of another photon, not selecting ..." << std::endl; 
          continue;
        } 
      }
    }

    if (partMap) partMap->AddAt(nacc, iPart);

    AliAODMCParticle *newPart = new ((*partOut)[nacc]) AliAODMCParticle(*part);
    newPart->SetGeneratorIndex(part->GetGeneratorIndex());

    nacc++;
  }
  AliDebugStream(5) << "Particles in: " << partIn->GetEntries() << ", out: " << partOut->GetEntries() << std::endl;
}

Bool_t AliEmcalMCTrackSelector::AcceptParticle(AliAODMCParticle* part) const
{
  if (!part) return kFALSE;

  Int_t partPdgCode = TMath::Abs(part->PdgCode());

  if (fOnlyHIJING && (part->GetGeneratorIndex() != 0)) return kFALSE;

  if (fEtaMax > 0. && TMath::Abs(part->Eta()) > fEtaMax) return kFALSE;

  if (fRejectNK && (partPdgCode == 130 || partPdgCode == 2112)) return kFALSE;

  if (fChargedMC && part->Charge() == 0) return kFALSE;

  if (fOnlyPhysPrim && !part->IsPhysicalPrimary()) return kFALSE;

  return kTRUE;
}

AliEmcalMCTrackSelector* AliEmcalMCTrackSelector::AddTaskMCTrackSelector(TString outname, Bool_t nk, Bool_t ch, Double_t etamax, Bool_t physPrim)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMCTrackSelector", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMCTrackSelector", "This task requires an input event handler");
    return NULL;
  }

  // Init the task and do settings
  TString name("AliEmcalMCTrackSelector_");
  name += outname;
  AliEmcalMCTrackSelector *eTask = new AliEmcalMCTrackSelector(name);
  eTask->SetParticlesOutName(outname);
  eTask->SetRejectNK(nk);
  eTask->SetChargedMC(ch);
  eTask->SetEtaMax(etamax);
  eTask->SetOnlyPhysPrim(physPrim);

  // Final settings, pass to manager and set the containers
  mgr->AddTask(eTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  return eTask;
}
