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

  // PROTON SELECTION
  const int multiplicity = fEvent->GetMultiplicity();
  UInt_t filterBitProton = fTrackCutsPartProton->GetFilterBit();
  bool useTPConlyTrack = (filterBitProton == 128);
  static std::vector<AliFemtoDreamBasePart> Particles;
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles;
  AntiParticles.clear();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(evt->GetTrack(iTrack));
    fProtonTrack->SetTrack(esdTrack, fMCEvent, multiplicity, useTPConlyTrack);
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      Particles.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      AntiParticles.push_back(*fProtonTrack);
    }
  }

  // LAMBDA SELECTION
  fV0Cuts->SelectV0(fInputEvent, fMCEvent);
  fAntiV0Cuts->SelectV0(fInputEvent, fMCEvent);

  static std::vector<AliFemtoDreamBasePart> Decays;
  static std::vector<AliFemtoDreamBasePart> AntiDecays;
  CastToVector(fV0Cuts->GetV0s(), Decays, fMCEvent);
  CastToVector(fAntiV0Cuts->GetV0s(), AntiDecays, fMCEvent);

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  if (!fIsLightweight) {
    fPhotonQA->PhotonQA(fInputEvent, fMCEvent, fGammaArray);
  }
  static std::vector<AliFemtoDreamBasePart> Gammas;
  Gammas.clear();
  CastToVector(Gammas, fInputEvent);

  fPairCleaner->CleanDecay(&Decays, 0);
  fPairCleaner->CleanDecay(&AntiDecays, 1);
  fPairCleaner->CleanDecay(&Gammas, 2);
  fPairCleaner->CleanDecayAndDecay(&Decays, &AntiDecays, 3);
  fPairCleaner->CleanDecayAndDecay(&Decays, &Gammas, 4);
  fPairCleaner->CleanDecayAndDecay(&AntiDecays, &Gammas, 5);

  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, Gammas, Decays);
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, Gammas, AntiDecays);

  std::vector<AliFemtoDreamBasePart> sigma0particles, sigma0sidebandUp,
      sigma0sidebandLow, antiSigma0particles, antiSigma0sidebandUp,
      antiSigma0sidebandLow, sigma0lambda, antiSigma0lambda, sigma0photon, antiSigma0photon;

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
  if (fCheckDaughterCF) {
    fPairCleaner->CleanTrackAndDecay(&Particles, &sigma0lambda, 6);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &antiSigma0lambda, 7);
    fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 8);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 9);
    fPairCleaner->CleanTrackAndDecay(&Particles, &sigma0photon, 10);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &antiSigma0photon, 11);
    fPairCleaner->CleanTrackAndDecay(&Particles, &Gammas, 12);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &Gammas, 13);
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(sigma0particles);
  fPairCleaner->StoreParticle(antiSigma0particles);
  fPairCleaner->StoreParticle(sigma0sidebandUp);
  fPairCleaner->StoreParticle(antiSigma0sidebandUp);
  fPairCleaner->StoreParticle(sigma0sidebandLow);
  fPairCleaner->StoreParticle(antiSigma0sidebandLow);
  if (fCheckDaughterCF) {
    fPairCleaner->StoreParticle(sigma0lambda);
    fPairCleaner->StoreParticle(antiSigma0lambda);
    fPairCleaner->StoreParticle(Decays);
    fPairCleaner->StoreParticle(AntiDecays);
    fPairCleaner->StoreParticle(sigma0photon);
    fPairCleaner->StoreParticle(antiSigma0photon);
    fPairCleaner->StoreParticle(Gammas);
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
    std::vector<AliFemtoDreamBasePart> &container, const AliVEvent *inputEvent) {
  for (int iGamma = 0; iGamma < fGammaArray->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fGammaArray->At(iGamma));
    if (!PhotonCandidate) continue;

    auto pos =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel1());
    auto neg =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel2());
    if (!pos || !neg) continue;

    // cut on the DCA to the PV in transverse plane
    PhotonCandidate->CalculateDistanceOfClossetApproachToPrimVtx(inputEvent->GetPrimaryVertex());
    const float DCAr = PhotonCandidate->GetDCArToPrimVtx();
    if (DCAr > 0.75) {
      continue;
    }

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
    container.push_back( { PhotonCandidate, pos, neg, inputEvent });
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliFemtoDreamBasePart> &particlesOut,
    std::vector<AliFemtoDreamBasePart> &particlesIn) {
  particlesOut.clear();
  // Randomly pick one of the particles in the container
  if (particlesIn.size() > 0) {
    particlesOut.push_back(particlesIn[fRandom->Rndm() * particlesIn.size()]);
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

  const int nPairs = (fCheckDaughterCF) ? 14 : 6;
  fPairCleaner =
      new AliFemtoDreamPairCleaner(nPairs, 6, fConfig->GetMinimalBookingME());
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
