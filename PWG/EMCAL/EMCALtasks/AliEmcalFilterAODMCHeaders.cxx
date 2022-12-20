/**************************************************************************************
 * Copyright (C) 2022, Copyright Holders of the ALICE Collaboration                   *
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

////////////////////////////////////////////////////////////////////////////////////////
// Note: Must be run BEFORE track selector task! This task uses information from the  //
// stack, which is not avalable after running the track selector task.                //
////////////////////////////////////////////////////////////////////////////////////////

#include "AliEmcalFilterAODMCHeaders.h"

#include <iostream>

#include <TClonesArray.h>
#include <TPDGCode.h>

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliNamedArrayI.h"
#include "AliAODMCHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODEvent.h"
#include "TObjString.h"


#include "AliLog.h"

ClassImp(AliEmcalFilterAODMCHeaders)

AliEmcalFilterAODMCHeaders::AliEmcalFilterAODMCHeaders() :
  AliAnalysisTaskSE("AliEmcalFilterAODMCHeaders"),
  fParticlesOutName("MCParticlesNotRejected"),
  fParticlesMapName(""),
  fTracksOutName("MCTracksNotRejected"),
  fTracksInName("tracks"),
  fClustersOutName("MCClustersNotRejected"),
  fClustersInName("caloClusters"),
  fInit(kFALSE),
  fParticlesIn(nullptr),
  fParticlesOut(nullptr),
  fParticlesMap(nullptr),
  fTracksIn(nullptr),
  fTracksOut(nullptr),
  fClustersIn(nullptr),
  fClustersOut(nullptr),
  fEvent(nullptr),
  fMC(nullptr),
  fDisabled(kFALSE),
  fDebugLevel(0),
  fAODMCTrackArray(NULL),
  fMCHeader(0)
{
}

AliEmcalFilterAODMCHeaders::AliEmcalFilterAODMCHeaders(const char *name) :
  AliAnalysisTaskSE(name),
  fParticlesOutName("MCParticlesNotRejected"),
  fParticlesMapName(""),
  fTracksOutName("MCTracksNotRejected"),
  fTracksInName("tracks"),
  fClustersOutName("MCClustersNotRejected"),
  fClustersInName("caloClusters"),
  fInit(kFALSE),
  fParticlesIn(nullptr),
  fParticlesOut(nullptr),
  fParticlesMap(nullptr),
  fTracksIn(nullptr),
  fTracksOut(nullptr),
  fClustersIn(nullptr),
  fClustersOut(nullptr),
  fEvent(nullptr),
  fMC(nullptr),
  fDisabled(kFALSE),
  fDebugLevel(0),
  fAODMCTrackArray(NULL),
  fMCHeader(0)
{
}

AliEmcalFilterAODMCHeaders::~AliEmcalFilterAODMCHeaders() {
  if(fAODMCTrackArray){
    delete[] fAODMCTrackArray;
    fAODMCTrackArray = 0x0;
  }
}

//______________________________________________________________________________
void AliEmcalFilterAODMCHeaders::UserExec(Option_t *)
{
  if(fDisabled) return;

  if(!fInit){
    fEvent = InputEvent();
    if(!fEvent){
      AliErrorStream() << "Could not retrieve event! Returning" << std::endl;
      return;
    }
    if(!fEvent->InheritsFrom("AliAODEvent")){
      AliErrorStream() << "Event type AOD not found! Returning..." << std::endl;
      return;
    }

    fMC = MCEvent();
    if(!fMC){
      AliErrorStream() << "Could not retrieve MC event! Returning" << std::endl;
      fDisabled = kTRUE;
      return;
    }

    fMCHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!fMCHeader) {
       AliError(Form("MC header not found! Task '%s' will be disabled!", GetName()));
       return;
    }

    if(fEvent->FindListObject(fParticlesOutName)){ // the output array is already present in the array!
      AliErrorStream() << "The output array "  << fParticlesOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else if(fEvent->FindListObject(fTracksOutName)){
      AliErrorStream() << "The output array "  << fTracksOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else if(fEvent->FindListObject(fClustersOutName)){
      AliErrorStream() << "The output array "  << fClustersOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else{  // copy the array from the standard AOD collections and filter if requested
      fParticlesMapName = fParticlesOutName;
      fParticlesMapName += "_Map";

      if(fEvent->FindListObject(fParticlesMapName)){
        AliErrorStream() << "The output array map " << fParticlesMapName << " is already present in the event! Task will be disabled." << std::endl;
        fDisabled = kTRUE;
        return;
      }else{
        fParticlesMap = new AliNamedArrayI(fParticlesMapName, 99999);
        fEvent->AddObject(fParticlesMap);
      }

      fParticlesOut = new TClonesArray("AliAODMCParticle");  // the output will always be of AliAODMCParticle, regardless of the input
      fParticlesOut->SetName(fParticlesOutName);
      fEvent->AddObject(fParticlesOut);

      fTracksOut = new TClonesArray("AliAODTrack");
      fTracksOut->SetName(fTracksOutName);
      fEvent->AddObject(fTracksOut);

      fClustersOut = new TClonesArray("AliAODCaloCluster");
      fClustersOut->SetName(fClustersOutName);
      fEvent->AddObject(fClustersOut);

      fParticlesIn = static_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if(!fParticlesIn){
        AliErrorStream() << "Could not retrieve AOD MC particles! Task will be disabled." << std::endl;
        fDisabled = kTRUE;
        return;
      }
      TClass *cl = fParticlesIn->GetClass();
      if(!cl->GetBaseClass("AliAODMCParticle")){
        AliErrorStream() << GetName() << ": Collection  " << AliAODMCParticle::StdBranchName() << " %s does not contain AliAODMCParticle! Task will be disabled." << std::endl;
        fDisabled = kTRUE;
        fParticlesIn = 0;
        return;
      }

      fTracksIn = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTracksInName));
      if (!fTracksIn) {
        AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data()));
        return;
      }

      fClustersIn = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fClustersInName));
      if (!fClustersIn) {
        AliError(Form("Could not retrieve clusters %s!", fClustersInName.Data()));
        return;
      }
    }
    fInit = kTRUE;
  }

  ProcessParticles();
  LinkMothers();
  ProcessTracks();
  ProcessClusters();
}

//______________________________________________________________________________
void AliEmcalFilterAODMCHeaders::LinkMothers(){
  for (Int_t iPart = 0; iPart < fParticlesOut->GetEntries(); iPart++) {
    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(fParticlesOut->At(iPart));
    if(part->GetNDaughters() <= 0) continue;
    if(part->GetDaughterFirst() > -1 && part->GetDaughterLast() > -1) {
      ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->SetDaughter(0, fParticlesMap->At(part->GetDaughterFirst()));
      ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->SetDaughter(1, fParticlesMap->At(part->GetDaughterLast()));
      ((AliAODMCParticle*)(fParticlesOut->At( ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->GetDaughterFirst() )))->SetMother(iPart);
      ((AliAODMCParticle*)(fParticlesOut->At( ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->GetDaughterLast() )))->SetMother(iPart);
    }
  }
}

//______________________________________________________________________________
void AliEmcalFilterAODMCHeaders::ProcessClusters(){
  // clear container
  fClustersOut->Delete();

  Int_t ncl = fClustersIn->GetEntriesFast();
  if(ncl == 0){
    return;
  }
  AliDebugStream(3) << "NClusters = " << ncl << std::endl;

  for(Long_t iCluster = 0, nacc = 0; iCluster < ncl; ++iCluster){
    AliVCluster* clus = static_cast<AliAODCaloCluster*>(fClustersIn->At(iCluster));
    if (!clus) continue;

    if(!IsParticleFromBGEvent(clus->GetLabelAt(0))){
      AliDebugStream(3) << "Parent particle for cluster " << iCluster << " is from bg event. Skipping cluster..." << std::endl;
      continue;
    }

    new ((*fClustersOut)[nacc++]) AliAODCaloCluster(*(AliAODCaloCluster*)clus);
  }
  AliDebugStream(3) << "Clusters in:  " << fClustersIn->GetEntries() << ", out: " << fClustersOut->GetEntries() << std::endl;
}

//______________________________________________________________________________
void AliEmcalFilterAODMCHeaders::ProcessTracks(){
  // clear container
  fTracksOut->Delete();

  int ntr = fTracksIn->GetEntriesFast();
  if(ntr == 0){
    return;
  }
  AliDebugStream(3) << "NTracks = " << ntr << std::endl;

  for (Int_t iTrack=0, nacc=0; iTrack < ntr; ++iTrack) {
    AliVTrack* track = static_cast<AliAODTrack*>(fTracksIn->At(iTrack));
    if (!track) continue;

    if(!IsParticleFromBGEvent(TMath::Abs(track->GetLabel()))){
      AliDebugStream(3) << "Parent particle for track " << iTrack << " is from bg event. Skipping track..." << std::endl;
      continue;
    }

    new ((*fTracksOut)[nacc++]) AliAODTrack(*(AliAODTrack*)track);
  }
  AliDebugStream(3) << "Tracks in   : " << fTracksIn->GetEntries() << ", out: " << fTracksOut->GetEntries() << std::endl;
}

//______________________________________________________________________________
void AliEmcalFilterAODMCHeaders::ProcessParticles()
{
  if (!fParticlesIn) return;
  AliDebugStream(5) << "Particles in: " << fParticlesIn->GetName() << ", out: " << fParticlesOut->GetName() << std::endl;
  // clear container (normally a null operation as the event should clean it already)
  fParticlesOut->Delete();

  Int_t Nparticles = fParticlesIn->GetEntriesFast();
  Int_t nacc = 0;

  AliDebugStream(2) << "Total number of particles = " << Nparticles << std::endl;
  if (fParticlesMap) {
    // clear particles map
    fParticlesMap->Clear();
    if (fParticlesMap->GetSize() <= Nparticles) fParticlesMap->Set(Nparticles*2);
  }

  AliDebugStream(2) << "Total number of particles = " << Nparticles << std::endl;

  // loop over particles
  for (Int_t iPart = 0; iPart < Nparticles; iPart++) {
    if (fParticlesMap) fParticlesMap->AddAt(-1, iPart);

    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(fParticlesIn->At(iPart));
    if (!part) continue;

    // rejection criteria
    if(!IsParticleFromBGEvent(iPart)){
      AliDebugStream(3) << "Particle " << iPart << " is from bg event. Skipping particle..." << std::endl;
      continue;
    }

    if (fParticlesMap) fParticlesMap->AddAt(nacc, iPart);

    AliAODMCParticle *newPart = new ((*fParticlesOut)[nacc]) AliAODMCParticle(*part);
    newPart->SetGeneratorIndex(part->GetGeneratorIndex());
    newPart->SetFlag(part->GetFlag());

    nacc++;
  }
  AliDebugStream(3) << "Particles in: " << fParticlesIn->GetEntries() << ", out: " << fParticlesOut->GetEntries() << std::endl;
}

//______________________________________________________________________________
Bool_t AliEmcalFilterAODMCHeaders::IsParticleFromBGEvent(Int_t index){
  if(index < 0) return false; // No Particle

  bool accepted = false;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
    if(!aodMCParticle) return false; // no particle

    if(!aodMCParticle->IsPrimary()){
      if( aodMCParticle->GetMother() < 0) return false;// material particle, return 0
      return IsParticleFromBGEvent(aodMCParticle->GetMother());
    }

    TString nameGen = GetNameGenerator(index,fMCHeader);
    if(!(nameGen.IsWhitespace() || nameGen.Contains("EPOS"))) accepted = true;
    if (fDebugLevel > 1 && !accepted) AliDebugStream(1) << "rejected:" << index << std::endl;
  }

  return accepted;
}

//______________________________________________________________________________
TString AliEmcalFilterAODMCHeaders::GetNameGenerator(Int_t label, AliAODMCHeader* header){
  /// get the name of the generator that produced a given particle
  Int_t nsumpart = 0;
  TList *lh = header->GetCocktailHeaders();
  Int_t nh = lh->GetEntries();
  AliGenEventHeader *gh = 0;
  for(Int_t i=0;i<nh;i++){
    gh = (AliGenEventHeader*)lh->At(i);
    TString genname = gh->GetName();
    Int_t npart = gh->NProduced();
    if(label>=nsumpart && label<(nsumpart+npart)) return genname;
    nsumpart+=npart;
  }
  TString empty = "";
  return empty;
}

//______________________________________________________________________________
AliEmcalFilterAODMCHeaders* AliEmcalFilterAODMCHeaders::AddTaskFilterAODMCHeaders(const TString nParticlesOut, const TString nTracksOut, const TString nClustersOut, const Int_t debug)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFilterAODMCHeaders", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFilterAODMCHeaders", "This task requires an input event handler");
    return NULL;
  }

  // Init the task and do settings
  TString name("AliEmcalFilterAODMCHeaders_");
  name += nParticlesOut;
  AliEmcalFilterAODMCHeaders *eTask = new AliEmcalFilterAODMCHeaders(name);
  eTask->SetParticlesOutName(nParticlesOut);
  eTask->SetTracksOutName(nTracksOut);
  eTask->SetClustersOutName(nClustersOut);
  eTask->SetDebugLevel(debug);

  // Final settings, pass to manager and set the containers
  mgr->AddTask(eTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  return eTask;
}
