#include "AliAnalysisTaskAODSigma0Femto.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskAODSigma0Femto)

    //____________________________________________________________________________________________________
    AliAnalysisTaskAODSigma0Femto::AliAnalysisTaskAODSigma0Femto()
    : AliAnalysisTaskSE("AliAnalysisTaskAODSigma0Femto"),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fRandom(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fLambda(nullptr),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fTrackBufferSize(2500),
      fGTI(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fLambdaHistList(nullptr),
      fLambdaHistMCList(nullptr),
      fAntiLambdaHistList(nullptr),
      fAntiLambdaHistMCList(nullptr),
      fPhotonHistList(nullptr),
      fSigmaHistList(nullptr),
      fAntiSigmaHistList(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr) {
  fRandom = new TRandom3();
}

//____________________________________________________________________________________________________
AliAnalysisTaskAODSigma0Femto::AliAnalysisTaskAODSigma0Femto(const char *name,
                                                             const bool isMC)
    : AliAnalysisTaskSE(name),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fRandom(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fLambda(nullptr),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(isMC),
      fIsLightweight(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fTrackBufferSize(2500),
      fGTI(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fLambdaHistList(nullptr),
      fLambdaHistMCList(nullptr),
      fAntiLambdaHistList(nullptr),
      fAntiLambdaHistMCList(nullptr),
      fPhotonHistList(nullptr),
      fSigmaHistList(nullptr),
      fAntiSigmaHistList(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());
  DefineOutput(9, TList::Class());
  DefineOutput(10, TList::Class());
  DefineOutput(11, TList::Class());

  fRandom = new TRandom3();
  if (fIsMC) {
    DefineOutput(12, TList::Class());
    DefineOutput(13, TList::Class());
    DefineOutput(14, TList::Class());
    DefineOutput(15, TList::Class());
  }
}

//____________________________________________________________________________________________________
AliAnalysisTaskAODSigma0Femto::~AliAnalysisTaskAODSigma0Femto() {
  delete fPartColl;
  delete fPairCleaner;
  delete fProtonTrack;
  delete fLambda;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODSigma0Femto::UserExec(Option_t * /*option*/) {
  AliVEvent *fInputEvent = InputEvent();
  if (fIsMC) fMCEvent = MCEvent();

  // PREAMBLE - CHECK EVERYTHING IS THERE
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

  if (!fTrackCutsPartProton || !fTrackCutsPartAntiProton) {
    AliError("Proton Cuts missing");
    return;
  }

  if (!fV0Cuts || !fAntiV0Cuts) {
    AliError("V0 Cuts missing");
    return;
  }

  if (!fSigmaCuts || !fAntiSigmaCuts) {
    AliError("Sigma0 Cuts missing");
    return;
  }

  // EVENT SELECTION
  AliAODEvent *evt = static_cast<AliAODEvent *>(fInputEvent);
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Particles;
  std::vector<AliFemtoDreamBasePart> AntiParticles;
  const int multiplicity = fEvent->GetMultiplicity();
  fProtonTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    fProtonTrack->SetTrack(track, multiplicity);
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      Particles.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      AntiParticles.push_back(*fProtonTrack);
    }
  }

  // LAMBDA SELECTION
  std::vector<AliFemtoDreamBasePart> Decays;
  std::vector<AliFemtoDreamBasePart> AntiDecays;
  fLambda->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
       iv0 < static_cast<TClonesArray *>(evt->GetV0s())->GetEntriesFast();
       iv0++) {
    AliAODv0 *v0 = evt->GetV0(iv0);
    fLambda->Setv0(evt, v0, multiplicity);
    if (fV0Cuts->isSelected(fLambda)) {
      Decays.push_back(*fLambda);
    }
    if (fAntiV0Cuts->isSelected(fLambda)) {
      AntiDecays.push_back(*fLambda);
    }
  }

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  std::vector<AliFemtoDreamBasePart> Gammas;
  fPhotonCuts->PhotonCuts(evt, fMCEvent, fGammaArray, Gammas);

  // Sigma0 selection
  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, Gammas, Decays);

  // Sigma0 selection
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, Gammas, AntiDecays);

  std::vector<AliFemtoDreamBasePart> sigma0particles, sigma0sidebandUp,
      sigma0sidebandLow, antiSigma0particles, antiSigma0sidebandUp,
      antiSigma0sidebandLow;

  CastToVector(sigma0particles, fSigmaCuts->GetSigma());
  CastToVector(sigma0sidebandUp, fSigmaCuts->GetSidebandUp());
  CastToVector(sigma0sidebandLow, fSigmaCuts->GetSidebandDown());

  CastToVector(antiSigma0particles, fAntiSigmaCuts->GetSigma());
  CastToVector(antiSigma0sidebandUp, fAntiSigmaCuts->GetSidebandUp());
  CastToVector(antiSigma0sidebandLow, fAntiSigmaCuts->GetSidebandDown());

  fPairCleaner->CleanTrackAndDecay(&Particles, &sigma0particles, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &antiSigma0particles, 1);
  fPairCleaner->CleanTrackAndDecay(&Particles, &sigma0sidebandUp, 2);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &antiSigma0sidebandUp, 3);
  fPairCleaner->CleanTrackAndDecay(&Particles, &sigma0sidebandLow, 4);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &antiSigma0sidebandLow, 5);

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(sigma0particles);
  fPairCleaner->StoreParticle(antiSigma0particles);
  fPairCleaner->StoreParticle(sigma0sidebandUp);
  fPairCleaner->StoreParticle(antiSigma0sidebandUp);
  fPairCleaner->StoreParticle(sigma0sidebandLow);
  fPairCleaner->StoreParticle(antiSigma0sidebandLow);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());

  // flush the data
  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fLambdaHistList);
  PostData(6, fAntiLambdaHistList);
  PostData(7, fPhotonHistList);
  PostData(8, fSigmaHistList);
  PostData(9, fAntiSigmaHistList);
  PostData(10, fResultList);
  PostData(11, fResultQAList);
  if (fIsMC) {
    PostData(12, fTrackCutHistMCList);
    PostData(13, fAntiTrackCutHistMCList);
    PostData(14, fLambdaHistMCList);
    PostData(15, fAntiLambdaHistMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODSigma0Femto::CastToVector(
    std::vector<AliFemtoDreamBasePart> &particlesOut,
    std::vector<AliFemtoDreamBasePart> &particlesIn) {
  particlesOut.clear();
  // Randomly pick one of the particles in the container
  if (particlesIn.size() > 0) {
    particlesOut.push_back(particlesIn[fRandom->Rndm() * particlesIn.size()]);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODSigma0Femto::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODSigma0Femto::StoreGlobalTrackReference(
    AliAODTrack *track) {
  // see AliFemtoDreamAnalysis for details
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
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODSigma0Femto::UserCreateOutputObjects() {
  fGTI = new AliAODTrack *[fTrackBufferSize];

  fEvent = new AliFemtoDreamEvent(true, !fIsLightweight, fTrigger);

  fProtonTrack = new AliFemtoDreamTrack();
  fProtonTrack->SetUseMCInfo(fIsMC);

  fLambda = new AliFemtoDreamv0();
  fLambda->SetUseMCInfo(fIsMC);
  fLambda->SetPDGCode(fV0Cuts->GetPDGv0());
  fLambda->SetPDGDaughterPos(fV0Cuts->GetPDGPosDaug());
  fLambda->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fLambda->SetPDGDaughterNeg(fV0Cuts->GetPDGNegDaug());
  fLambda->GetNegDaughter()->SetUseMCInfo(fIsMC);

  const int nPairs = 6;
  fPairCleaner =
      new AliFemtoDreamPairCleaner(nPairs, 0, fConfig->GetMinimalBookingME());
  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  fQA = new TList();
  fQA->SetName("QA");
  fQA->SetOwner(kTRUE);

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (!fIsLightweight) {
    if (fV0Reader->GetConversionCuts() &&
        fV0Reader->GetConversionCuts()->GetCutHistograms()) {
      fQA->Add(fV0Reader->GetConversionCuts()->GetCutHistograms());
    }
  }

  if (fEvtCuts) {
    fEvtCuts->InitQA();
    fQA->Add(fEvent->GetEvtCutList());
    if (fEvtCuts->GetHistList() && !fIsLightweight) {
      fEvtHistList = fEvtCuts->GetHistList();
    } else {
      fEvtHistList = new TList();
      fEvtHistList->SetName("EvtCuts");
      fEvtHistList->SetOwner(true);
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }

  if (!fConfig->GetMinimalBookingME() && fPairCleaner &&
      fPairCleaner->GetHistList()) {
    fQA->Add(fPairCleaner->GetHistList());
  }

  fTrackCutsPartProton->Init("Proton");
  // necessary for the non-min booking case
  fTrackCutsPartProton->SetName("Proton");

  if (fTrackCutsPartProton && fTrackCutsPartProton->GetQAHists()) {
    fTrackCutHistList = fTrackCutsPartProton->GetQAHists();
    if (fIsMC && fTrackCutsPartProton->GetMCQAHists() &&
        fTrackCutsPartProton->GetIsMonteCarlo()) {
      fTrackCutHistMCList = fTrackCutsPartProton->GetMCQAHists();
    }
  }

  fTrackCutsPartAntiProton->Init("Anti-proton");
  // necessary for the non-min booking case
  fTrackCutsPartAntiProton->SetName("Anti-proton");

  if (fTrackCutsPartAntiProton && fTrackCutsPartAntiProton->GetQAHists()) {
    fAntiTrackCutHistList = fTrackCutsPartAntiProton->GetQAHists();
    if (fIsMC && fTrackCutsPartAntiProton->GetMCQAHists() &&
        fTrackCutsPartAntiProton->GetIsMonteCarlo()) {
      fAntiTrackCutHistMCList = fTrackCutsPartAntiProton->GetMCQAHists();
    }
  }

  fV0Cuts->Init();
  fAntiV0Cuts->Init();

  if (fV0Cuts->GetQAHists()) {
    fLambdaHistList = fV0Cuts->GetQAHists();
  } else {
    fLambdaHistList = new TList();
    fLambdaHistList->SetName("v0Cuts");
    fLambdaHistList->SetOwner();
  }
  if (fAntiV0Cuts->GetQAHists()) {
    fAntiLambdaHistList = fAntiV0Cuts->GetQAHists();
  } else {
    fAntiLambdaHistList = new TList();
    fAntiLambdaHistList->SetName("Antiv0Cuts");
    fAntiLambdaHistList->SetOwner();
  }
  if (!fV0Cuts->GetMinimalBooking()) {
    if (fV0Cuts->GetIsMonteCarlo()) {
      fLambdaHistMCList = fV0Cuts->GetMCQAHists();
    }
  } else {
    fLambdaHistMCList = new TList();
    fLambdaHistMCList->SetName("MCv0Cuts");
    fLambdaHistMCList->SetOwner();
  }
  if (!fAntiV0Cuts->GetMinimalBooking()) {
    if (fAntiV0Cuts->GetIsMonteCarlo()) {
      fAntiLambdaHistMCList = fAntiV0Cuts->GetMCQAHists();
    }
  } else {
    fAntiLambdaHistMCList = new TList();
    fAntiLambdaHistMCList->SetName("MCAntiv0Cuts");
    fAntiLambdaHistMCList->SetOwner();
  }

  if (fSigmaCuts) fSigmaCuts->InitCutHistograms(TString("Sigma0"));
  if (fAntiSigmaCuts) fAntiSigmaCuts->InitCutHistograms(TString("AntiSigma0"));

  if (fPhotonCuts) {
    fPhotonCuts->InitCutHistograms(TString("Photon"));
    fPhotonHistList = fPhotonCuts->GetCutHistograms();
  } else {
    fPhotonHistList = new TList();
    fPhotonHistList->SetName("V0_Photon");
    fPhotonHistList->SetOwner(true);
  }

  if (fSigmaCuts && fSigmaCuts->GetCutHistograms()) {
    fSigmaHistList = fSigmaCuts->GetCutHistograms();
  }

  if (fAntiSigmaCuts && fAntiSigmaCuts->GetCutHistograms()) {
    fAntiSigmaHistList = fAntiSigmaCuts->GetCutHistograms();
  }

  if (fPartColl && fPartColl->GetHistList()) {
    fResultList = fPartColl->GetHistList();
  }
  if (!fConfig->GetMinimalBookingME() && fPartColl && fPartColl->GetQAList()) {
    fResultQAList = fPartColl->GetQAList();
  } else {
    fResultQAList = new TList();
    fResultQAList->SetName("ResultsQA");
    fResultQAList->SetOwner(true);
  }

  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fLambdaHistList);
  PostData(6, fAntiLambdaHistList);
  PostData(7, fPhotonHistList);
  PostData(8, fSigmaHistList);
  PostData(9, fAntiSigmaHistList);
  PostData(10, fResultList);
  PostData(11, fResultQAList);
  if (fIsMC) {
    PostData(12, fTrackCutHistMCList);
    PostData(13, fAntiTrackCutHistMCList);
    PostData(14, fLambdaHistMCList);
    PostData(15, fAntiLambdaHistMCList);
  }
}
