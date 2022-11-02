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

#include "AliEmcalRejectMCBackground.h"

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
#include "AliESDEvent.h"
#include "TObjString.h"


#include "AliLog.h"

ClassImp(AliEmcalRejectMCBackground)

AliEmcalRejectMCBackground::AliEmcalRejectMCBackground() :
  AliAnalysisTaskSE("AliEmcalRejectMCBackground"),
  fHeaderList(nullptr),
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
  fEsdEvent(nullptr),
  fIsESD(kFALSE),
  fDisabled(kFALSE),
  fDebugLevel(0),
  fnHeaders(0),
  fAODMCTrackArray(NULL),
  fNotRejectedStart(0),
  fNotRejectedEnd(0),
  fGeneratorNames(0),
  fAddedSignalPDGCode(0),
  fSignalRejection(0)
{
}

AliEmcalRejectMCBackground::AliEmcalRejectMCBackground(const char *name) :
  AliAnalysisTaskSE(name),
  fHeaderList(nullptr),
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
  fEsdEvent(nullptr),
  fIsESD(kFALSE),
  fDisabled(kFALSE),
  fDebugLevel(0),
  fnHeaders(0),
  fAODMCTrackArray(NULL),
  fNotRejectedStart(0),
  fNotRejectedEnd(0),
  fGeneratorNames(0),
  fAddedSignalPDGCode(0),
  fSignalRejection(0)
{
}

AliEmcalRejectMCBackground::~AliEmcalRejectMCBackground() {
  if(fAODMCTrackArray){
    delete[] fAODMCTrackArray;
    fAODMCTrackArray = 0x0;
  }
}

void AliEmcalRejectMCBackground::UserExec(Option_t *)
{
  if(fDisabled) return;

  if(!fInit){
    fEvent = InputEvent();
    if(!fEvent){
      AliErrorStream() << "Could not retrieve event! Returning" << std::endl;
      return;
    }
    if(fEvent->InheritsFrom("AliESDEvent")) fIsESD = kTRUE;
    else if(fEvent->InheritsFrom("AliAODEvent")) fIsESD = kFALSE;
    else{
      AliErrorStream() << "Event type ESD or AOD not found! Returning..." << std::endl;
      fIsESD = kFALSE;
      return;
    }

    fMC = MCEvent();
    if(!fMC){
      AliErrorStream() << "Could not retrieve MC event! Returning" << std::endl;
      fDisabled = kTRUE;
      return;
    }

    if(fIsESD){
      fEsdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
      if (!fEsdEvent) {
        AliErrorStream() << "Could not retrieve ESD event! Returning" << std::endl;
        return;
      }
    }

    TObject *objParticle = fEvent->FindListObject(fParticlesOutName);
    TObject *objTrack    = fEvent->FindListObject(fTracksOutName);
    TObject *objCluster  = fEvent->FindListObject(fClustersOutName);
    if(objParticle){ // the output array is already present in the array!
      AliErrorStream() << "The output array "  << fParticlesOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else if(objTrack){
      AliErrorStream() << "The output array "  << fTracksOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else if(objCluster){
      AliErrorStream() << "The output array "  << fClustersOutName <<  " is already present in the event! Task will be disabled." << std::endl;
      fDisabled = kTRUE;
      return;
    }else{  // copy the array from the standard ESD/AOD collections, and filter if requested
      fParticlesOut = new TClonesArray("AliAODMCParticle");  // the output will always be of AliAODMCParticle, regardless of the input
      fParticlesOut->SetName(fParticlesOutName);
      fEvent->AddObject(fParticlesOut);
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

      if(fIsESD) fTracksOut = new TClonesArray("AliESDtrack");
      else fTracksOut = new TClonesArray("AliAODTrack");
      fTracksOut->SetName(fTracksOutName);
      fEvent->AddObject(fTracksOut);

      if(fIsESD) fClustersOut = new TClonesArray("AliESDCaloCluster");
      else fClustersOut = new TClonesArray("AliAODCaloCluster");
      fClustersOut->SetName(fClustersOutName);
      fEvent->AddObject(fClustersOut);

      if(!fIsESD){
        fParticlesIn = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
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
        //if (!fTracksIn->GetClass()->GetBaseClass("AliVParticle")) {
        //  AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTracksInName.Data()));
        //  return;
        //}

        fClustersIn = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fClustersInName));
        if (!fClustersIn) {
          AliError(Form("Could not retrieve clusters %s!", fClustersInName.Data()));
          return;
        }
        //if (!fClustersIn->GetClass()->GetBaseClass("AliVParticle")) {
        //  AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fClustersInName.Data()));
        //  return;
        //}

      }
    }
    fInit = kTRUE;
  }

  if(fSignalRejection != 0){
    if(fIsESD) GetNotRejectedParticles(fSignalRejection, GetAcceptedHeader(), fMC);
    else GetNotRejectedParticles(fSignalRejection, GetAcceptedHeader(), fEvent);
  }

  if (fIsESD) CreateParticleMap(fEvent, fMC, fParticlesOut, fParticlesMap);
  else CreateParticleMapAOD(fEvent, fMC, fParticlesIn, fParticlesOut, fParticlesMap);

  LinkMothers();

  ProcessTracks();
  ProcessClusters();
}

void AliEmcalRejectMCBackground::LinkMothers(){
  for (Int_t iPart = 0; iPart < fParticlesOut->GetEntries(); iPart++) {
    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(fParticlesOut->At(iPart));
    if(part->GetNDaughters() <= 0) continue;
    if(part->GetDaughterFirst() > -1 && part->GetDaughterLast() > -1) {
      ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->SetDaughter(0, fParticlesMap->At(part->GetDaughterFirst()));
      ((AliAODMCParticle*)(fParticlesOut->At(iPart)))->SetDaughter(1, fParticlesMap->At(part->GetDaughterLast()));
    }
  }
}

void AliEmcalRejectMCBackground::ProcessClusters(){
  // clear container
  fClustersOut->Delete();

  Int_t ncl;
  if(fIsESD) ncl = fEvent->GetNumberOfCaloClusters();
  else ncl = fClustersIn->GetEntriesFast();
  if(ncl == 0){
    return;
  }
  AliDebugStream(3) << "NClusters = " << ncl << std::endl;

  for(Long_t iCluster = 0, nacc = 0; iCluster < ncl; ++iCluster){
    AliVCluster* clus = NULL;
    if(fIsESD) clus = static_cast<AliESDCaloCluster*>(fEvent->GetCaloCluster(iCluster));
    else clus = static_cast<AliAODCaloCluster*>(fClustersIn->At(iCluster));
    if (!clus) continue;

    if(IsParticleFromBGEvent(clus->GetLabelAt(0), fMC, fEvent, fDebugLevel) < 1){
      AliDebugStream(3) << "Parent particle for cluster " << iCluster << " is from bg event. Skipping cluster..." << std::endl;
      continue;
    }

    if(fIsESD) new ((*fClustersOut)[nacc++]) AliESDCaloCluster(*(AliESDCaloCluster*)clus);
    else new ((*fClustersOut)[nacc++]) AliAODCaloCluster(*(AliAODCaloCluster*)clus);
  }
  AliDebugStream(3) << "Clusters in:  " << fClustersIn->GetEntries() << ", out: " << fClustersOut->GetEntries() << std::endl;
}

void AliEmcalRejectMCBackground::ProcessTracks(){
  // clear container
  fTracksOut->Delete();

  int ntr;
  if(fIsESD) ntr = fEsdEvent->GetNumberOfTracks();
  else ntr = fTracksIn->GetEntriesFast();
  if(ntr == 0){
    return;
  }
  AliDebugStream(3) << "NTracks = " << ntr << std::endl;

  for (Int_t iTrack=0, nacc=0; iTrack < ntr; ++iTrack) {
    AliVTrack* track = NULL;
    if(fIsESD) track = static_cast<AliESDtrack*>(fEsdEvent->GetTrack(iTrack));
    else track = static_cast<AliAODTrack*>(fTracksIn->At(iTrack));
    if (!track) continue;

    if(IsParticleFromBGEvent(TMath::Abs(track->GetLabel()), fMC, fEvent, fDebugLevel) < 1){
      AliDebugStream(3) << "Parent particle for track " << iTrack << " is from bg event. Skipping track..." << std::endl;
      continue;
    }

    if(fIsESD) new ((*fTracksOut)[nacc++]) AliESDtrack(*(AliESDtrack*)track);
    else new ((*fTracksOut)[nacc++]) AliAODTrack(*(AliAODTrack*)track);
  }
  AliDebugStream(3) << "Tracks in   : " << fTracksIn->GetEntries() << ", out: " << fTracksOut->GetEntries() << std::endl;
}

void AliEmcalRejectMCBackground::CreateParticleMap(AliVEvent *event, AliMCEvent* mcEvent, TClonesArray* partOut, AliNamedArrayI* partMap)
{
  // clear container (normally a null operation as the event should clean it already)
  partOut->Delete();

  Int_t Nparticles = mcEvent->GetNumberOfTracks();
  Int_t nprim = mcEvent->GetNumberOfPrimaries();
  Int_t nacc = 0;

  if (partMap) {
    // clear particles map
    partMap->Clear();
    if (partMap->GetSize() <= Nparticles) partMap->Set(Nparticles*2);
  }

  AliDebugStream(3) << "Number of particles: " << Nparticles << std::endl;

  // loop over particles
  for (Int_t iPart = 0; iPart < Nparticles; iPart++) {
    if (partMap) partMap->AddAt(-1, iPart);

    AliMCParticle* part = static_cast<AliMCParticle*>(mcEvent->GetTrack(iPart));
    if (!part) continue;

    // Set flags
    Bool_t isPhysPrim = mcEvent->IsPhysicalPrimary(iPart);
    Int_t flag = 0;
    if (iPart < nprim) flag |= AliAODMCParticle::kPrimary;
    if (isPhysPrim) flag |= AliAODMCParticle::kPhysicalPrim;
    if (mcEvent->IsSecondaryFromWeakDecay(iPart)) flag |= AliAODMCParticle::kSecondaryFromWeakDecay;
    if (mcEvent->IsSecondaryFromMaterial(iPart)) flag |= AliAODMCParticle::kSecondaryFromMaterial;
    AliDebugStream(3) << "Particle " << iPart << ": pt = " << part->Pt() << ", PDG = " << part->PdgCode() <<
        ", mother " << part->GetMother() <<
        ", kPrimary? " <<  Bool_t((flag & AliAODMCParticle::kPrimary) != 0) <<
        ", kPhysicalPrim? " <<  Bool_t((flag & AliAODMCParticle::kPhysicalPrim) != 0) <<
        ", kSecondaryFromWeakDecay? " <<  Bool_t((flag & AliAODMCParticle::kSecondaryFromWeakDecay) != 0) <<
        ", kSecondaryFromMaterial? " << Bool_t((flag & AliAODMCParticle::kSecondaryFromMaterial) != 0) <<
        ", nacc = " << nacc <<
        ", iPart = " << iPart <<
        std::endl;

    // rejection criteriaI
    if(IsParticleFromBGEvent(iPart, mcEvent, event, fDebugLevel) < 1) continue;

    if(partMap) partMap->AddAt(nacc, iPart);

    AliAODMCParticle parttmp(part, iPart, flag);
    parttmp.SetGeneratorIndex(part->GetGeneratorIndex());
    parttmp.SetStatus(part->Particle()->GetStatusCode());
    parttmp.SetMCProcessCode(part->Particle()->GetUniqueID());
    new ((*partOut)[nacc]) AliAODMCParticle(parttmp);

    nacc++;
  }
}

void AliEmcalRejectMCBackground::CreateParticleMapAOD(AliVEvent *event, AliMCEvent *mcEvent, TClonesArray* partIn, TClonesArray* partOut, AliNamedArrayI* partMap)
{
  if (!partIn) return;
  AliDebugStream(5) << "Particles in: " << partIn->GetName() << ", out: " << partOut->GetName() << std::endl;
  // clear container (normally a null operation as the event should clean it already)
  partOut->Delete();

  Int_t Nparticles = partIn->GetEntriesFast();
  Int_t nacc = 0;

  AliDebugStream(2) << "Total number of particles = " << Nparticles << std::endl;
  if (partMap) {
    // clear particles map
    partMap->Clear();
    if (partMap->GetSize() <= Nparticles) partMap->Set(Nparticles*2);
  }

  AliDebugStream(2) << "Total number of particles = " << Nparticles << std::endl;

  // loop over particles
  for (Int_t iPart = 0; iPart < Nparticles; iPart++) {
    if (partMap) partMap->AddAt(-1, iPart);

    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(partIn->At(iPart));
    if (!part) continue;

    // rejection criteria
    if(IsParticleFromBGEvent(iPart, mcEvent, event, fDebugLevel) < 1){
      AliDebugStream(3) << "Particle " << iPart << " is from bg event. Skipping particle..." << std::endl;
      continue;
    }

    if (partMap) partMap->AddAt(nacc, iPart);

    AliAODMCParticle *newPart = new ((*partOut)[nacc]) AliAODMCParticle(*part);
    newPart->SetGeneratorIndex(part->GetGeneratorIndex());
    newPart->SetFlag(part->GetFlag());

    nacc++;
  }
  AliDebugStream(3) << "Particles in: " << partIn->GetEntries() << ", out: " << partOut->GetEntries() << std::endl;
}

void AliEmcalRejectMCBackground::GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *event){

  // rejection == 0: No rejection
  // rejection == 1: Min bias header
  // rejection == 2: User string array of accepted headers

  if(rejection==0) return;

  AliGenCocktailEventHeader *cHeader = 0x0;
  AliAODMCHeader *cHeaderAOD         = 0x0;
  Bool_t headerFound                 = kFALSE;
  AliMCEvent *fMCEvent               = 0x0;
  TClonesArray *fMCEventAOD          = 0x0;

  if(event->IsA()==AliMCEvent::Class()){
    if(dynamic_cast<AliMCEvent*>(event)){
      cHeader                  = dynamic_cast<AliGenCocktailEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
      fMCEvent                 = dynamic_cast<AliMCEvent*>(event);
      if(cHeader) headerFound  = kTRUE;
    }
  }
  if(event->IsA()==AliAODEvent::Class()){ // event is a AODEvent in case of AOD
    cHeaderAOD                 = dynamic_cast<AliAODMCHeader*>(event->FindListObject(AliAODMCHeader::StdBranchName()));
    fMCEventAOD                = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if(cHeaderAOD) headerFound = kTRUE;
  }

  if(fDebugLevel > 0) AliDebugStream(3) << "event starts here" << std::endl;

  if(headerFound){
    TList *genHeaders      = 0x0;
    if(cHeader) genHeaders = cHeader->GetHeaders();
    if(cHeaderAOD){
      genHeaders           = cHeaderAOD->GetCocktailHeaders();
      if(genHeaders->GetEntries()==1){
        fSignalRejection = 0;
        return;
      }
    }

    AliGenEventHeader* gh  = 0;
    fnHeaders              = 0;
    Int_t firstindexA      = 0;
    Int_t lastindexA       = -1;

    if(rejection == 1) fnHeaders = 1; // MinBiasHeader

    // first loop over headers in event and accepted headers from HeaderList.
    // purpose: find out number of headers that are both in the event and on the HeaderList (fnHeaders)
    if(rejection == 2){ // TList of Headers Names
      for(Int_t i = 0; i<genHeaders->GetEntries();i++){
        gh                    = (AliGenEventHeader*)genHeaders->At(i);
        TString GeneratorName = gh->GetName();
        lastindexA            = lastindexA + gh->NProduced();
        if (fDebugLevel > 0 ) AliDebugStream(1) << i << "\t" << GeneratorName.Data() << std::endl;
        for(Int_t j = 0; j<HeaderList->GetEntries();j++){
          TString GeneratorInList   = ((TObjString*)HeaderList->At(j))->GetString();
          if(fDebugLevel > 0 )  AliDebugStream(1) << GeneratorInList.Data() << std::endl;
          if(GeneratorInList.Contains(GeneratorName) ){
            if(fDebugLevel > 0 ) AliDebugStream(1) << "accepted" << std::endl;
            if(GeneratorInList.BeginsWith("PARAM") || GeneratorInList.CompareTo("BOX") == 0 ){
              if(fMCEvent){
                if (fDebugLevel > 0 ) AliDebugStream(1) << "cond 2: " << fnHeaders << std::endl;
                fnHeaders++;
                continue;
              }
              if (fMCEventAOD){
                if (fDebugLevel > 0 ) AliDebugStream(1) << "cond 2: " << fnHeaders << std::endl;
                fnHeaders++;
                continue;
              }
              continue;
            }
            if(GeneratorName.CompareTo(GeneratorInList) == 0 ){
              if (fDebugLevel > 0 ) AliDebugStream(1) << "cond 3: "<< fnHeaders << std::endl;
              fnHeaders++;
              continue;
            }
          }
          if(GeneratorName.Contains(GeneratorInList) ){
            if(GeneratorInList.Contains("Pythia") || GeneratorInList.Contains("pythia") || GeneratorInList.Contains("PYTHIA")){
              if (fDebugLevel > 0 ) AliDebugStream(1) << "Pythia header" << std::endl;
              fnHeaders++;
              continue;
            }
          }
        }
        firstindexA       = firstindexA + gh->NProduced();
      }
    }
    if (fDebugLevel > 0 ) AliDebugStream(1) << "number of headers: " <<fnHeaders << std::endl;

    fNotRejectedStart.resize(fnHeaders);
    fNotRejectedEnd.resize(fnHeaders);
    fGeneratorNames.resize(fnHeaders);

    if(rejection == 1){
      fNotRejectedStart[0]    = 0;
      fNotRejectedEnd[0]      = ((AliGenEventHeader*)genHeaders->At(0))->NProduced()-1;
      fGeneratorNames[0]      = ((AliGenEventHeader*)genHeaders->At(0))->GetName();

      if (fDebugLevel > 0 ){
        AliDebugStream(1) << 0 << "\t" <<fGeneratorNames[0] << "\t" << fNotRejectedStart[0] << "\t" <<fNotRejectedEnd[0] << std::endl;
      }
      return;
    }

    // only get here for rejection==2: second loop over event headers and accepted headers from HeaderList.
    // purpose: feed indeces corresponding to tracks from accepted headers into fNotRejectedStart and fNotRejectedEnd
    Int_t firstindex        = 0;
    Int_t lastindex         =  -1;
    Int_t number            = 0;

    for(Int_t i = 0; i<genHeaders->GetEntries();i++){
      gh = (AliGenEventHeader*)genHeaders->At(i);
      TString GeneratorName = gh->GetName();
      lastindex             = lastindex + gh->NProduced();
      if(fDebugLevel > 0 ) AliDebugStream(1) << i << "\t" << GeneratorName.Data() << std::endl;
      for(Int_t j = 0; j<HeaderList->GetEntries();j++){
        TString GeneratorInList = ((TObjString*)HeaderList->At(j))->GetString();
        if(GeneratorInList.Contains(GeneratorName) ){
          if (GeneratorInList.Contains("PARAM") || GeneratorInList.CompareTo("BOX") == 0 ){
            if(fMCEvent){
              fNotRejectedStart[number] = firstindex;
              fNotRejectedEnd[number] = lastindex;
              fGeneratorNames[number] = GeneratorName;
              number++;
              continue;
            }
            if ( fMCEventAOD){
              fNotRejectedStart[number] = firstindex;
              fNotRejectedEnd[number] = lastindex;
              fGeneratorNames[number] = GeneratorName;
              number++;
              continue;
            }
            continue;
          } else if(GeneratorName.CompareTo(GeneratorInList) == 0 ){
            fNotRejectedStart[number] = firstindex;
            fNotRejectedEnd[number] = lastindex;
            fGeneratorNames[number] = GeneratorName;
            if (fDebugLevel > 0 )  AliDebugStream(1) << "Number of particles produced for: " << i << "\t" << GeneratorName.Data() << "\t" << lastindex-firstindex+1 << std::endl;
            number++;
            continue;
          }
        }
        if(GeneratorName.Contains(GeneratorInList) ){
          if(GeneratorInList.Contains("Pythia") || GeneratorInList.Contains("pythia") || GeneratorInList.Contains("PYTHIA")){
            fNotRejectedStart[number] = firstindex;
            fNotRejectedEnd[number] = lastindex;
            fGeneratorNames[number] = GeneratorName;
            number++;
            continue;
          }
        }
      }
      firstindex           = firstindex + gh->NProduced();
    }
    if (fDebugLevel > 0 ) {
      for (Int_t i = 0; i < number; i++){
        AliDebugStream(1) << i << "\t" <<fGeneratorNames[i] << "\t" << fNotRejectedStart[i] << "\t" <<fNotRejectedEnd[i] << std::endl;
      }
    }
  } else { // No Cocktail Header Found
    fNotRejectedStart.resize(1);
    fNotRejectedEnd.resize(1);

    fnHeaders             = 1;
    fNotRejectedStart[0]       = 0;
    fNotRejectedEnd[0]         = static_cast<AliMCEvent*>(event)->GetNumberOfPrimaries()-1;
    if (rejection > 1){
      fNotRejectedStart[0]     = -1;
      fNotRejectedEnd[0]       = -1;
    }

    fGeneratorNames.resize(1);
    fGeneratorNames[0]         = "NoCocktailGeneratorFound";
  }

}

Int_t AliEmcalRejectMCBackground::IsParticleFromBGEvent(Int_t index, AliMCEvent *mcEvent, AliVEvent *InputEvent, Int_t debug ){

  //   if (debug > 2 ) AliDebugStream(1) << index << std::endl;
  if(index < 0) return 0; // No Particle

  Int_t accepted = 0;
  if(!InputEvent || InputEvent->IsA()==AliESDEvent::Class()){
    if(!mcEvent) return 0; // no mcEvent available, return 0
    if(index >= mcEvent->GetNumberOfPrimaries()){ // initial particle is secondary particle
      if( ((AliMCParticle*) mcEvent->GetTrack(index))->GetMother() < 0) return 0; // material particle, return 0
      return IsParticleFromBGEvent(((AliMCParticle*) mcEvent->GetTrack(index))->GetMother(),mcEvent,InputEvent, debug);
    }
    for(Int_t i = 0;i<fnHeaders;i++){
      //       if (debug > 2 ) AliDebugStream(1) << "header " << fGeneratorNames[i].Data() << ":"<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << std::endl;
      if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
        if (debug > 1 ) AliDebugStream(1) << "accepted:" << index << "\t header " << fGeneratorNames[i].Data()  << ": "<< fNotRejectedStart[i] << "\t" << fNotRejectedEnd[i] << std::endl;
        accepted = 1;
        if(i == 0) accepted = 2; // MB Header
      }
    }
    if (debug > 1 && !accepted) AliDebugStream(1) << "rejected:" << index << std::endl;
  }
  else if(InputEvent->IsA()==AliAODEvent::Class()){
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(InputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray){
      AliAODMCParticle *aodMCParticle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
      if(!aodMCParticle) return 0; // no particle

      if(!aodMCParticle->IsPrimary()){
        if( aodMCParticle->GetMother() < 0) return 0;// material particle, return 0
        return IsParticleFromBGEvent(aodMCParticle->GetMother(),mcEvent,InputEvent, debug);
      }
      index = TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index))->GetLabel());

      for(Int_t i = 0;i<fnHeaders;i++){
        if(index >= fNotRejectedStart[i] && index <= fNotRejectedEnd[i]){
          accepted = 1;
          if(i == 0) accepted = 2; // MB Header
        }
      }
    }
  }

  return accepted;
}

AliEmcalRejectMCBackground* AliEmcalRejectMCBackground::AddTaskRejectMCBackground(const TString nParticlesOut, const TString nTracksOut, const TString nClustersOut, const Int_t signalRejection, const Int_t debug)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskRejectMCBackground", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskRejectMCBackground", "This task requires an input event handler");
    return NULL;
  }

  // Specify the accpeted headers.
  TList *HeaderList = new TList();

  TObjString *Header1 = new TObjString("AliGenPythiaEventHeader");
  TObjString *Header2 = new TObjString("Pythia8Jets");

  HeaderList->Add(Header1);
  HeaderList->Add(Header2);

  // Init the task and do settings
  TString name("AliEmcalRejectMCBackground_");
  name += nParticlesOut;
  AliEmcalRejectMCBackground *eTask = new AliEmcalRejectMCBackground(name);
  eTask->SetParticlesOutName(nParticlesOut);
  eTask->SetTracksOutName(nTracksOut);
  eTask->SetClustersOutName(nClustersOut);
  eTask->SetSignalRejection(signalRejection);
  eTask->SetDebugLevel(debug);
  eTask->SetAcceptedHeader(HeaderList);

  // Final settings, pass to manager and set the containers
  mgr->AddTask(eTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  return eTask;
}
