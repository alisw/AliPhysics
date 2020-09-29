#include "AliAnalysisTaskNanoPPCoalescence.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskNanoPPCoalescence)
AliAnalysisTaskNanoPPCoalescence::AliAnalysisTaskNanoPPCoalescence()
  : AliAnalysisTaskSE(),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fConfig(nullptr),
    fIsMC(false),
    fIsMCTruth(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fGTI(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fTrackBufferSize(2500) {
}

AliAnalysisTaskNanoPPCoalescence::AliAnalysisTaskNanoPPCoalescence(const char *name, const bool isMC)
  : AliAnalysisTaskSE(name),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fConfig(nullptr),
    fIsMC(isMC),
    fIsMCTruth(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fGTI(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fTrackBufferSize(2500) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Results
  DefineOutput(5, TList::Class());  //Output for the Results QA
  if (fIsMC) {
    DefineOutput(6, TList::Class());  //Output for the Proton MC
    DefineOutput(7, TList::Class());  //Output for the AntiProton MC
  }
}

//---------------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPPCoalescence::~AliAnalysisTaskNanoPPCoalescence() {
  delete fEvent;
  delete fTrack;
  delete fProtonTrack;
  delete fAntiProtonTrack;
  delete fPairCleaner;
  delete fPartColl;
}

void AliAnalysisTaskNanoPPCoalescence::UserCreateOutputObjects() {

  fGTI = new AliVTrack*[fTrackBufferSize];

  if (!fEvtCuts) {
    AliError("No Event cuts \n");
  } else {
    fEvtCuts->InitQA();
  }

  if (!fProtonTrack) {
    AliError("No Proton cuts \n");
  } else {
    fProtonTrack->Init();
    fProtonList = fProtonTrack->GetQAHists();
    if (fIsMC) {
      fProtonMCList = fProtonTrack->GetMCQAHists();
    }
  }

  if (!fAntiProtonTrack) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProtonTrack->Init();
    fAntiProtonList = fAntiProtonTrack->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fAntiProtonTrack->GetMCQAHists();
    }
  }

  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,fConfig->GetMinimalBookingME());
    std::cout<<fConfig->GetMinimalBookingME()<<endl;
    fPairCleaner = new AliFemtoDreamPairCleaner(0,0,fConfig->GetMinimalBookingME());
  }

  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates(), false);
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  if (!fEvtCuts->GetMinimalBooking()) {
    fEvtList = fEvtCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fResults);
  PostData(5, fResultsQA);


  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(6, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(7, fAntiProtonMCList);
  }

}

void AliAnalysisTaskNanoPPCoalescence::UserExec(Option_t *option) {
  AliVEvent *fInputEvent = InputEvent();

  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }
  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }
  if (!fEvtCuts) {
    AliError("Event Cuts missing");
    return;
  }
  if (!fProtonTrack || !fAntiProtonTrack) {
    AliError("Proton Cuts missing");
    return;
  }

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEvtCuts->isSelected(fEvent))
    return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard NanoAOD");
      return;
    }

    StoreGlobalTrackReference(track);
  }

  static std::vector<AliFemtoDreamBasePart> Proton;
  Proton.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProton;
  AntiProton.clear();
  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
    if (!track)
      continue;
    fTrack->SetTrack(track, fInputEvent, multiplicity);

    if (fIsMCTruth && fIsMC) {
      int mcpdg;
      mcpdg = fTrack->GetMCPDGCode();
      if ((mcpdg == 2212) && (fProtonTrack->isSelected(fTrack))) {
        Proton.push_back(*fTrack);
      }
      if ((mcpdg == -2212) && (fAntiProtonTrack->isSelected(fTrack))) {
        AntiProton.push_back(*fTrack);
      }
    } else {
      if (fProtonTrack->isSelected(fTrack)) {
        Proton.push_back(*fTrack);
      }
      if (fAntiProtonTrack->isSelected(fTrack)) {
        AntiProton.push_back(*fTrack);
      }
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
        if (mcPart->GetPdgCode() == fProtonTrack->GetPDGCode()) {
          fProtonTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiProtonTrack->GetPDGCode()) {
          fAntiProtonTrack->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Proton, &AntiProton, 0);
  fPairCleaner->StoreParticle(Proton);
  fPairCleaner->StoreParticle(AntiProton);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fResults);
  PostData(5, fResultsQA);

  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(6, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(7, fAntiProtonMCList);
  }

}

void AliAnalysisTaskNanoPPCoalescence::ResetGlobalTrackReference() {
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskNanoPPCoalescence::StoreGlobalTrackReference(AliVTrack *track) {
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }
  if (fGTI[trackID]) {

    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

