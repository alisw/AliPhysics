/*
 * AliFemtoDreamAnalysis.cxx
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */
#include <vector>
#include "AliLog.h"
#include "AliFemtoDreamAnalysis.h"
#include "TClonesArray.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include <iostream>
ClassImp(AliFemtoDreamAnalysis)
AliFemtoDreamAnalysis::AliFemtoDreamAnalysis()
    : fMVPileUp(false),
      fEvtCutQA(false),
      fQA(),
      fFemtoTrack(nullptr),
      fFemtov0(nullptr),
      fFemtoCasc(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fTrackCuts(nullptr),
      fAntiTrackCuts(nullptr),
      fv0Cuts(nullptr),
      fAntiv0Cuts(nullptr),
      fCascCuts(nullptr),
      fAntiCascCuts(nullptr),
      fPairCleaner(nullptr),
      fControlSample(nullptr),
      fTrackBufferSize(0),
      fGTI(nullptr),
      fConfig(nullptr),
      fPartColl(nullptr),
      fIsMC(false) {

}

AliFemtoDreamAnalysis::AliFemtoDreamAnalysis(
    const AliFemtoDreamAnalysis& analysis)
    : fMVPileUp(analysis.fMVPileUp),
      fEvtCutQA(analysis.fEvtCutQA),
      fQA(nullptr),
      fFemtoTrack(nullptr),
      fFemtov0(nullptr),
      fFemtoCasc(nullptr),
      fEvent(nullptr),
      fEvtCuts(analysis.fEvtCuts),
      fTrackCuts(analysis.fTrackCuts),
      fAntiTrackCuts(analysis.fAntiTrackCuts),
      fv0Cuts(analysis.fv0Cuts),
      fAntiv0Cuts(analysis.fAntiv0Cuts),
      fCascCuts(analysis.fCascCuts),
      fAntiCascCuts(analysis.fAntiCascCuts),
      fPairCleaner(nullptr),
      fControlSample(nullptr),
      fTrackBufferSize(analysis.fTrackBufferSize),
      fGTI(nullptr),
      fConfig(analysis.fConfig),
      fPartColl(nullptr),
      fIsMC(analysis.fIsMC) {
  Init(fIsMC, AliVEvent::kINT7);
}

AliFemtoDreamAnalysis& AliFemtoDreamAnalysis::operator=(
    const AliFemtoDreamAnalysis& analysis) {
  if (this != &analysis) {
    this->fMVPileUp = analysis.fMVPileUp;
    this->fEvtCutQA = analysis.fEvtCutQA;
    this->fQA = nullptr;
    this->fFemtoTrack = nullptr;
    this->fFemtov0 = nullptr;
    this->fFemtoCasc = nullptr;
    this->fEvent = nullptr;
    this->fEvtCuts = analysis.fEvtCuts;
    this->fTrackCuts = analysis.fTrackCuts;
    this->fAntiTrackCuts = analysis.fAntiTrackCuts;
    this->fv0Cuts = analysis.fv0Cuts;
    this->fAntiv0Cuts = analysis.fAntiv0Cuts;
    this->fCascCuts = analysis.fCascCuts;
    this->fAntiCascCuts = analysis.fAntiCascCuts;
    this->fPairCleaner = nullptr;
    this->fControlSample = nullptr;
    this->fTrackBufferSize = analysis.fTrackBufferSize;
    this->fGTI = nullptr;
    this->fConfig = analysis.fConfig;
    this->fPartColl = nullptr;
    this->fIsMC = analysis.fIsMC;
    this->Init(this->fIsMC, AliVEvent::kINT7);
  }
  return *this;
}

AliFemtoDreamAnalysis::~AliFemtoDreamAnalysis() {
  if (fEvent) {
    delete fEvent;
  }
  if (fFemtoTrack) {
    delete fFemtoTrack;
  }
  if (fFemtov0) {
    delete fFemtov0;
  }
  if (fFemtoCasc) {
    delete fFemtoCasc;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fControlSample) {
    delete fControlSample;
  }
}

void AliFemtoDreamAnalysis::Init(bool isMonteCarlo, UInt_t trigger) {
  fIsMC = isMonteCarlo;
  fFemtoTrack = new AliFemtoDreamTrack();
  fFemtoTrack->SetUseMCInfo(isMonteCarlo);

  fFemtov0 = new AliFemtoDreamv0();
  fFemtov0->SetPDGCode(fv0Cuts->GetPDGv0());
  fFemtov0->SetUseMCInfo(isMonteCarlo);
  fFemtov0->SetPDGDaughterPos(fv0Cuts->GetPDGPosDaug());  //order +sign doesnt play a role
  fFemtov0->GetPosDaughter()->SetUseMCInfo(isMonteCarlo);
  fFemtov0->SetPDGDaughterNeg(fv0Cuts->GetPDGNegDaug());  //only used for MC Matching
  fFemtov0->GetNegDaughter()->SetUseMCInfo(isMonteCarlo);

  fFemtoCasc = new AliFemtoDreamCascade();
  fFemtoCasc->SetUseMCInfo(isMonteCarlo);
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fFemtoCasc->SetPDGCode(fCascCuts->GetPDGCodeCasc());
  fFemtoCasc->SetPDGDaugPos(fCascCuts->GetPDGCodePosDaug());
  fFemtoCasc->GetPosDaug()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->SetPDGDaugNeg(fCascCuts->GetPDGCodeNegDaug());
  fFemtoCasc->GetNegDaug()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->SetPDGDaugBach(fCascCuts->GetPDGCodeBach());
  fFemtoCasc->GetBach()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->Setv0PDGCode(fCascCuts->GetPDGv0());
  fEvtCuts->InitQA();
  fTrackCuts->Init();
  fAntiTrackCuts->Init();
  fv0Cuts->Init();
  fAntiv0Cuts->Init();
  fCascCuts->Init();
  fAntiCascCuts->Init();
  fGTI = new AliAODTrack*[fTrackBufferSize];
  fEvent = new AliFemtoDreamEvent(fMVPileUp, fEvtCutQA, trigger,
                                  fEvtCuts->GetUseAliEventCuts());
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  bool MinBooking = !((!fConfig->GetMinimalBookingME())
      || (!fConfig->GetMinimalBookingSample()));
  fPairCleaner = new AliFemtoDreamPairCleaner(4, 4, MinBooking);

  if (!MinBooking) {
    fQA = new TList();
    fQA->SetOwner();
    fQA->SetName("QA");
    fQA->Add(fPairCleaner->GetHistList());
    if (fEvtCutQA) {
      fQA->Add(fEvent->GetEvtCutList());
    }
  }
  if (fConfig->GetUseEventMixing()) {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample = new AliFemtoDreamControlSample(fConfig);
  }
  return;
}

void AliFemtoDreamAnalysis::ResetGlobalTrackReference() {
  //This method was inherited form H. Beck analysis

  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}
void AliFemtoDreamAnalysis::StoreGlobalTrackReference(AliAODTrack *track) {
  //This method was inherited form H. Beck analysis

  //bhohlweg@cern.ch: We ask for the Unique Track ID that points back to the
  //ESD. Seems like global tracks have a positive ID, Tracks with Filterbit
  //128 only have negative ID, this is used to match the Tracks later to their
  //global counterparts

  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see they don't have
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  //     //    "not adding it's info! "
  //     "ID: %d FilterMap: %d\n"
  //     ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }

  // Check that the id is positive
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  // Check id is not too big for buffer
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  // Warn if we overwrite a track
  if (fGTI[trackID]) {
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    // Imagine the other way around, the zero map zero clusters track
    // is stored and the good one wants to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }  // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  //     ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[trackID]) = track;
}

void AliFemtoDreamAnalysis::Make(AliAODEvent *evt) {
  if (!evt) {
    AliFatal("No Input Event");
  }
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) {
    return;
  }
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Particles;
  std::vector<AliFemtoDreamBasePart> AntiParticles;
  fFemtoTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    fFemtoTrack->SetTrack(track, fEvent->GetMultiplicity());
    if (fTrackCuts->isSelected(fFemtoTrack)) {
      Particles.push_back(*fFemtoTrack);
    }
    if (fAntiTrackCuts->isSelected(fFemtoTrack)) {
      AntiParticles.push_back(*fFemtoTrack);
    }
  }
  std::vector<AliFemtoDreamBasePart> Decays;
  std::vector<AliFemtoDreamBasePart> AntiDecays;
  //  Look for the lambda, store it in an event
  //  Get a V0 from the event:
  TClonesArray *v01 = static_cast<TClonesArray*>(evt->GetV0s());
  //number of V0s:

  fFemtov0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  int entriesV0 = v01->GetEntriesFast();
  for (int iv0 = 0; iv0 < entriesV0; iv0++) {
    AliAODv0 *v0 = evt->GetV0(iv0);
    fFemtov0->Setv0(evt, v0, fEvent->GetMultiplicity());
    if (fv0Cuts->isSelected(fFemtov0)) {
      Decays.push_back(*fFemtov0);
    }
    if (fAntiv0Cuts->isSelected(fFemtov0)) {
      AntiDecays.push_back(*fFemtov0);
    }
  }

  std::vector<AliFemtoDreamBasePart> XiDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiDecays;
  int numcascades = evt->GetNumberOfCascades();
  for (int iXi = 0; iXi < numcascades; ++iXi) {
    AliAODcascade *xi = evt->GetCascade(iXi);
    if (!xi)
      continue;
    fFemtoCasc->SetCascade(evt, xi);
    if (fCascCuts->isSelected(fFemtoCasc)) {
      XiDecays.push_back(*fFemtoCasc);
    }
    if (fAntiCascCuts->isSelected(fFemtoCasc)) {
      AntiXiDecays.push_back(*fFemtoCasc);
    }
  }
  //loop once over the MC stack to calculate Efficiency/Purity
  if (fIsMC) {
    AliAODInputHandler *eventHandler =
        dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler());
    AliMCEvent* fMC = eventHandler->MCEvent();

    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fTrackCuts->GetPDGCode()) {
          fTrackCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiTrackCuts->GetPDGCode()) {
          fAntiTrackCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fv0Cuts->GetPDGv0()) {
          fv0Cuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiv0Cuts->GetPDGv0()) {
          fAntiv0Cuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fCascCuts->GetPDGCodeCasc()) {
          fCascCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiCascCuts->GetPDGCodeCasc()) {
          fAntiCascCuts->FillGenerated(mcPart->Pt());
        }
      }
    }
  }
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 0);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);

  fPairCleaner->CleanDecay(&Decays, 0);
  fPairCleaner->CleanDecay(&AntiDecays, 1);
  fPairCleaner->CleanDecay(&XiDecays, 2);
  fPairCleaner->CleanDecay(&AntiXiDecays, 3);

  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(Decays);
  fPairCleaner->StoreParticle(AntiDecays);
  fPairCleaner->StoreParticle(XiDecays);
  fPairCleaner->StoreParticle(AntiXiDecays);
  if (fConfig->GetUseEventMixing()) {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                        fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
  }
}

void AliFemtoDreamAnalysis::Make(AliESDEvent *evt, AliMCEvent *mcEvent) {
  std::cout << "New Event \n";
  if (!evt) {
    AliFatal("No Input Event");
  }
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) {
    return;
  }
  std::vector<AliFemtoDreamBasePart> Particles;
  std::vector<AliFemtoDreamBasePart> AntiParticles;
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *track = static_cast<AliESDtrack *>(evt->GetTrack(iTrack));
    fFemtoTrack->SetTrack(track, mcEvent, fEvent->GetMultiplicity());
    if (fTrackCuts->isSelected(fFemtoTrack)) {
      Particles.push_back(*fFemtoTrack);
    }
    if (fAntiTrackCuts->isSelected(fFemtoTrack)) {
      AntiParticles.push_back(*fFemtoTrack);
    }
  }
  std::vector<AliFemtoDreamBasePart> Decays;
  std::vector<AliFemtoDreamBasePart> AntiDecays;

  for (int iv0 = 0; iv0 < evt->GetNumberOfV0s(); ++iv0) {
    AliESDv0 *v0 = evt->GetV0(iv0);
    fFemtov0->Setv0(evt, mcEvent, v0, fEvent->GetMultiplicity());
    if (fv0Cuts->isSelected(fFemtov0)) {
      Decays.push_back(*fFemtov0);
    }
    if (fAntiv0Cuts->isSelected(fFemtov0)) {
      AntiDecays.push_back(*fFemtov0);
    }
  }

  std::vector<AliFemtoDreamBasePart> XiDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiDecays;
  for (Int_t nCascade = 0; nCascade < evt->GetNumberOfCascades(); ++nCascade) {
    AliESDcascade *esdCascade = evt->GetCascade(nCascade);
    fFemtoCasc->SetCascade(evt, mcEvent, esdCascade);
    if (fCascCuts->isSelected(fFemtoCasc)) {
      XiDecays.push_back(*fFemtoCasc);
    }
    if (fAntiCascCuts->isSelected(fFemtoCasc)) {
      AntiXiDecays.push_back(*fFemtoCasc);
    }
  }
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 0);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);

  fPairCleaner->CleanDecay(&Decays, 0);
  fPairCleaner->CleanDecay(&AntiDecays, 1);
  fPairCleaner->CleanDecay(&XiDecays, 2);
  fPairCleaner->CleanDecay(&AntiXiDecays, 3);

  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(Decays);
  fPairCleaner->StoreParticle(AntiDecays);
  fPairCleaner->StoreParticle(XiDecays);
  fPairCleaner->StoreParticle(AntiXiDecays);

  if (fConfig->GetUseEventMixing()) {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                        fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample->SetEvent(fPairCleaner->GetCleanParticles(),fEvent);
  }
}

