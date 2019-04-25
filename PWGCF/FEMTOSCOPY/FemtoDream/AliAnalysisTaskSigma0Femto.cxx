#include "AliAnalysisTaskSigma0Femto.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskSigma0Femto)

    //____________________________________________________________________________________________________
    AliAnalysisTaskSigma0Femto::AliAnalysisTaskSigma0Femto()
    : AliAnalysisTaskSE("AliAnalysisTaskSigma0Femto"),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonQA(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fRandom(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fCheckDaughterCF(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fLambdaHistList(nullptr),
      fAntiLambdaHistList(nullptr),
      fPhotonHistList(nullptr),
      fSigmaHistList(nullptr),
      fAntiSigmaHistList(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr) {
  fRandom = new TRandom3();
}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Femto::AliAnalysisTaskSigma0Femto(const char *name,
                                                       const bool isMC)
    : AliAnalysisTaskSE(name),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonQA(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fRandom(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(isMC),
      fIsLightweight(false),
      fCheckDaughterCF(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fQA(nullptr),
      fEvtHistList(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fLambdaHistList(nullptr),
      fAntiLambdaHistList(nullptr),
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
  }
}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Femto::~AliAnalysisTaskSigma0Femto() {
  delete fPartColl;
  delete fPairCleaner;
  delete fProtonTrack;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::UserExec(Option_t * /*option*/) {
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

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (!fSigmaCuts || !fAntiSigmaCuts) {
    AliError("Sigma0 Cuts missing");
    return;
  }

  // EVENT SELECTION
  AliESDEvent *evt = static_cast<AliESDEvent *>(fInputEvent);
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) return;

  // LAMBDA SELECTION
  fV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // LAMBDA SELECTION
  fAntiV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  if (!fIsLightweight) {
    fPhotonQA->PhotonQA(fInputEvent, fMCEvent, fGammaArray);
  }
  std::vector<AliSigma0ParticleV0> gammaConvContainer;
  CastToVector(gammaConvContainer, fInputEvent);

  // Sigma0 selection
  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                 fV0Cuts->GetV0s());

  // Sigma0 selection
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                     fAntiV0Cuts->GetV0s());

  // Convert the Sigma0 into Femto particles
  static std::vector<AliFemtoDreamBasePart> sigma0particles;
  static std::vector<AliFemtoDreamBasePart> antiSigma0particles;
  static std::vector<AliFemtoDreamBasePart> sigma0sidebandUp;
  static std::vector<AliFemtoDreamBasePart> antiSigma0sidebandUp;
  static std::vector<AliFemtoDreamBasePart> sigma0sidebandLow;
  static std::vector<AliFemtoDreamBasePart> antiSigma0sidebandLow;

  static std::vector<AliFemtoDreamBasePart> sigma0lambda;
  static std::vector<AliFemtoDreamBasePart> sigma0photon;
  static std::vector<AliFemtoDreamBasePart> antiSigma0lambda;
  static std::vector<AliFemtoDreamBasePart> antiSigma0photon;

  CastToVector(fSigmaCuts->GetSigma(), sigma0particles, fMCEvent);
  CastToVector(fAntiSigmaCuts->GetSigma(), antiSigma0particles, fMCEvent);
  CastToVector(fSigmaCuts->GetSidebandUp(), sigma0sidebandUp, fMCEvent);
  CastToVector(fAntiSigmaCuts->GetSidebandUp(), antiSigma0sidebandUp, fMCEvent);
  CastToVector(fSigmaCuts->GetSidebandDown(), sigma0sidebandLow, fMCEvent);
  CastToVector(fAntiSigmaCuts->GetSidebandDown(), antiSigma0sidebandLow,
               fMCEvent);

  // Get the Sigma0 daughters
  if (fCheckDaughterCF) {
    static std::vector<AliSigma0ParticleV0> lambdaSigma;
    static std::vector<AliSigma0ParticleV0> photonSigma;
    static std::vector<AliSigma0ParticleV0> lambdaAntiSigma;
    static std::vector<AliSigma0ParticleV0> photonAntiSigma;

    fSigmaCuts->GetLambda(lambdaSigma);
    fSigmaCuts->GetPhoton(photonSigma);
    fAntiSigmaCuts->GetLambda(lambdaAntiSigma);
    fAntiSigmaCuts->GetPhoton(photonAntiSigma);

    CastToVector(lambdaSigma, sigma0lambda, fMCEvent);
    CastToVector(photonSigma, sigma0photon, fMCEvent);
    CastToVector(lambdaAntiSigma, antiSigma0lambda, fMCEvent);
    CastToVector(photonAntiSigma, antiSigma0photon, fMCEvent);
  }

  // PROTON SELECTION
  const int multiplicity = fEvent->GetMultiplicity();
  UInt_t filterBitProton = fTrackCutsPartProton->GetFilterBit();
  bool useTPConlyTrack = (filterBitProton == 128);
  static std::vector<AliFemtoDreamBasePart> particles;
  particles.clear();
  static std::vector<AliFemtoDreamBasePart> antiParticles;
  antiParticles.clear();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(evt->GetTrack(iTrack));
    fProtonTrack->SetTrack(esdTrack, fMCEvent, multiplicity, useTPConlyTrack);
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      particles.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      antiParticles.push_back(*fProtonTrack);
    }
  }

  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0particles, 0);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0particles, 1);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandUp, 2);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandUp, 3);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandLow, 4);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandLow, 5);
  if (fCheckDaughterCF) {
     fPairCleaner->CleanTrackAndDecay(&particles, &sigma0lambda, 6);
     fPairCleaner->CleanTrackAndDecay(&particles, &sigma0photon, 7);
     fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0lambda, 8);
     fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0photon, 9);
   }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(particles);
  fPairCleaner->StoreParticle(antiParticles);
  fPairCleaner->StoreParticle(sigma0particles);
  fPairCleaner->StoreParticle(antiSigma0particles);
  fPairCleaner->StoreParticle(sigma0sidebandUp);
  fPairCleaner->StoreParticle(antiSigma0sidebandUp);
  fPairCleaner->StoreParticle(sigma0sidebandLow);
  fPairCleaner->StoreParticle(antiSigma0sidebandLow);
  if (fCheckDaughterCF) {
     fPairCleaner->StoreParticle(sigma0lambda);
     fPairCleaner->StoreParticle(sigma0photon);
     fPairCleaner->StoreParticle(antiSigma0lambda);
     fPairCleaner->StoreParticle(antiSigma0photon);
   }

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
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticleV0> &container, const AliVEvent *inputEvent) {
  for (int iGamma = 0; iGamma < fGammaArray->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fGammaArray->At(iGamma));
    if (!PhotonCandidate) continue;

    auto pos =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel1());
    auto neg =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel2());
    if (!pos || !neg) continue;

    // pile up check
    if (fPhotonLegPileUpCut) {
      bool posTrackITS =
          (pos->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
           pos->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
      bool negTrackITS =
          (neg->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
           neg->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
      bool posTrackTOF = pos->GetTOFBunchCrossing() == 0;
      bool negTrackTOF = neg->GetTOFBunchCrossing() == 0;

      bool posTrackCombined = (posTrackITS || posTrackTOF);
      bool negTrackCombined = (negTrackITS || negTrackTOF);

      if (!posTrackCombined || !negTrackCombined) continue;
    }

    AliSigma0ParticleV0 phot(PhotonCandidate, pos, neg, inputEvent);
    if (fIsMC) {
      phot.MatchToMC(fMCEvent, 22, {{11, -11}});
    }
    container.push_back(phot);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticlePhotonMother> &sigmaContainer,
    std::vector<AliFemtoDreamBasePart> &particles, const AliMCEvent *mcEvent) {
  particles.clear();

  // Randomly pick one of the particles in the container
  if (sigmaContainer.size() > 0) {
    particles.push_back(
        {sigmaContainer[fRandom->Rndm() * sigmaContainer.size()], mcEvent});
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticleV0> &container,
    std::vector<AliFemtoDreamBasePart> &particles, const AliMCEvent *mcEvent) {
  particles.clear();
  for (const auto &part : container) {
    particles.push_back({part, mcEvent});
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::UserCreateOutputObjects() {
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

  fEvent = new AliFemtoDreamEvent(true, !fIsLightweight, fTrigger);

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

  fProtonTrack = new AliFemtoDreamTrack();
  fProtonTrack->SetUseMCInfo(fIsMC);

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

  if (fV0Cuts) fV0Cuts->InitCutHistograms(TString("Lambda"));
  if (fAntiV0Cuts) fAntiV0Cuts->InitCutHistograms(TString("AntiLambda"));
  if (fSigmaCuts) fSigmaCuts->InitCutHistograms(TString("Sigma0"));
  if (fAntiSigmaCuts) fAntiSigmaCuts->InitCutHistograms(TString("AntiSigma0"));

  if (fV0Cuts && fV0Cuts->GetCutHistograms()) {
    fLambdaHistList = fV0Cuts->GetCutHistograms();
  }

  if (fAntiV0Cuts && fAntiV0Cuts->GetCutHistograms()) {
    fAntiLambdaHistList = fAntiV0Cuts->GetCutHistograms();
  }

  if (!fIsLightweight) {
    fPhotonQA = new AliSigma0V0Cuts();
    fPhotonQA->SetIsMC(fIsMC);
    fPhotonQA->SetLightweight(false);
    fPhotonQA->SetPID(22);
    fPhotonQA->SetPosPID(AliPID::kElectron, 11);
    fPhotonQA->SetNegPID(AliPID::kElectron, -11);
    fPhotonQA->InitCutHistograms(TString("Photon"));
    fPhotonHistList = fPhotonQA->GetCutHistograms();
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

  const int nPairs = (fCheckDaughterCF) ? 10 : 6;
  fPairCleaner =
      new AliFemtoDreamPairCleaner(nPairs, 0, fConfig->GetMinimalBookingME());
  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  if (!fConfig->GetMinimalBookingME() && fPairCleaner &&
      fPairCleaner->GetHistList()) {
    fQA->Add(fPairCleaner->GetHistList());
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
  }
}
