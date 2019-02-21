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
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fIsRun1(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fOutputFemto(nullptr),
      fHistCorrelationPSigmaPLambda(),
      fHistCorrelationPSigmaPGamma(),
      fHistCorrelationPLambdaPGamma(),
      fHistCorrelationAntiPAntiSigmaAntiPAntiLambda(),
      fHistCorrelationAntiPAntiSigmaAntiPAntiGamma(),
      fHistCorrelationAntiPAntiLambdaAntiPAntiGamma() {}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Femto::AliAnalysisTaskSigma0Femto(const char *name)
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
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fIsRun1(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fOutputFemto(nullptr),
      fHistCorrelationPSigmaPLambda(),
      fHistCorrelationPSigmaPGamma(),
      fHistCorrelationPLambdaPGamma(),
      fHistCorrelationAntiPAntiSigmaAntiPAntiLambda(),
      fHistCorrelationAntiPAntiSigmaAntiPAntiGamma(),
      fHistCorrelationAntiPAntiLambdaAntiPAntiGamma() {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
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
  fEvent->SetEvent(static_cast<AliESDEvent *>(fInputEvent));
  if (!fEvtCuts->isSelected(fEvent)) return;

  // PROTON SELECTION
  UInt_t filterBitProton = fTrackCutsPartProton->GetFilterBit();
  bool useTPConlyTrack = (filterBitProton == 128);

  static std::vector<AliFemtoDreamBasePart> particles;
  particles.clear();
  static std::vector<AliFemtoDreamBasePart> antiParticles;
  antiParticles.clear();
  const AliESDEvent *evt = static_cast<AliESDEvent *>(fInputEvent);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    fProtonTrack->SetTrack(static_cast<AliESDtrack *>(evt->GetTrack(iTrack)),
                           fMCEvent, fEvent->GetMultiplicity(),
                           useTPConlyTrack);
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      particles.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      antiParticles.push_back(*fProtonTrack);
    }
  }

  // LAMBDA SELECTION
  fV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // LAMBDA SELECTION
  fAntiV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  if (!fIsLightweight) {
    fPhotonQA->PhotonQA(fInputEvent, fGammaArray);
  }
  std::vector<AliSigma0ParticleV0> gammaConvContainer;
  CastToVector(gammaConvContainer, fInputEvent);

  // Sigma0 selection
  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                 fV0Cuts->GetV0s());

  // Sigma0 selection
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                     fAntiV0Cuts->GetV0s());

  // Get the Sigma0 daughters
  static std::vector<AliSigma0ParticleV0> lambdaSigma;
  static std::vector<AliSigma0ParticleV0> photonSigma;
  static std::vector<AliSigma0ParticleV0> lambdaAntiSigma;
  static std::vector<AliSigma0ParticleV0> photonAntiSigma;
  if (!fIsLightweight) {
    fSigmaCuts->GetLambda(lambdaSigma);
    fSigmaCuts->GetPhoton(photonSigma);
    fAntiSigmaCuts->GetLambda(lambdaAntiSigma);
    fAntiSigmaCuts->GetPhoton(photonAntiSigma);
  }

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

  if (!fIsLightweight) {
    CastToVector(lambdaSigma, sigma0lambda, fMCEvent);
    CastToVector(photonSigma, sigma0photon, fMCEvent);
    CastToVector(lambdaAntiSigma, antiSigma0lambda, fMCEvent);
    CastToVector(photonAntiSigma, antiSigma0photon, fMCEvent);
  }

  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0particles, 0);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0particles, 1);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandUp, 2);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandUp, 3);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandLow, 4);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandLow, 5);
  if (!fIsLightweight) {
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
  if (!fIsLightweight) {
    fPairCleaner->StoreParticle(sigma0lambda);
    fPairCleaner->StoreParticle(sigma0photon);
    fPairCleaner->StoreParticle(antiSigma0lambda);
    fPairCleaner->StoreParticle(antiSigma0photon);
  }

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());

  if (!fIsLightweight) {
    FillCorrelationCorrelator(particles, sigma0particles,
                              fSigmaCuts->GetSigma(), false);
    FillCorrelationCorrelator(antiParticles, antiSigma0particles,
                              fAntiSigmaCuts->GetSigma(), true);
  }

  // flush the data
  PostData(1, fOutputContainer);
  PostData(2, fOutputFemto);
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
      const int label = phot.MatchToMC(fMCEvent, 22, {{11, -11}});
    }
    container.push_back(phot);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticlePhotonMother> &sigmaContainer,
    std::vector<AliFemtoDreamBasePart> &particles, const AliMCEvent *mcEvent) {
  particles.clear();
  for (const auto &sigma : sigmaContainer) {
    particles.push_back({sigma, mcEvent});
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
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
    fOutputContainer = nullptr;
  }
  if (fOutputContainer == nullptr) {
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  if (fOutputFemto != nullptr) {
    delete fOutputFemto;
    fOutputFemto = nullptr;
  }
  if (fOutputFemto == nullptr) {
    fOutputFemto = new TList();
    fOutputFemto->SetOwner(kTRUE);
  }

  fQA = new TList();
  fQA->SetName("EventCuts");
  fQA->SetOwner(true);

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (!fIsLightweight) {
    if (fV0Reader->GetEventCuts() &&
        fV0Reader->GetEventCuts()->GetCutHistograms()) {
      fQA->Add(fV0Reader->GetEventCuts()->GetCutHistograms());
    }
    if (fV0Reader->GetConversionCuts() &&
        fV0Reader->GetConversionCuts()->GetCutHistograms()) {
      fOutputContainer->Add(fV0Reader->GetConversionCuts()->GetCutHistograms());
    }
  }

  fEvent = new AliFemtoDreamEvent(true, !fIsLightweight, fTrigger);

  if (fEvtCuts) {
    fEvtCuts->InitQA();
    if (fEvtCuts->GetHistList() && !fIsLightweight) {
      fQA->Add(fEvtCuts->GetHistList());
    }
    if (fEvent->GetEvtCutList() && !fIsLightweight) {
      fQA->Add(fEvent->GetEvtCutList());
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }
  fOutputContainer->Add(fQA);

  fProtonTrack = new AliFemtoDreamTrack();
  fProtonTrack->SetUseMCInfo(fIsMC);

  fTrackCutsPartProton->Init("Proton");
  // necessary for the non-min booking case
  fTrackCutsPartProton->SetName("Proton");

  if (fTrackCutsPartProton && fTrackCutsPartProton->GetQAHists()) {
    fOutputFemto->Add(fTrackCutsPartProton->GetQAHists());
    if (fIsMC && !fTrackCutsPartProton->GetMinimalBooking() &&
        fTrackCutsPartProton->GetMCQAHists()) {
      fTrackCutsPartProton->SetMCName("MC_Proton");
      fOutputFemto->Add(fTrackCutsPartProton->GetMCQAHists());
    }
  }

  fTrackCutsPartAntiProton->Init("Anti-proton");
  // necessary for the non-min booking case
  fTrackCutsPartAntiProton->SetName("Anti-proton");

  if (fTrackCutsPartAntiProton && fTrackCutsPartAntiProton->GetQAHists()) {
    fOutputFemto->Add(fTrackCutsPartAntiProton->GetQAHists());
    if (fIsMC && !fTrackCutsPartAntiProton->GetMinimalBooking() &&
        fTrackCutsPartAntiProton->GetMCQAHists()) {
      fTrackCutsPartAntiProton->SetMCName("MC_Anti-proton");
      fOutputFemto->Add(fTrackCutsPartAntiProton->GetMCQAHists());
    }
  }

  if (fV0Cuts) fV0Cuts->InitCutHistograms(TString("Lambda"));
  if (fAntiV0Cuts) fAntiV0Cuts->InitCutHistograms(TString("AntiLambda"));
  if (fSigmaCuts) fSigmaCuts->InitCutHistograms(TString("Sigma0"));
  if (fAntiSigmaCuts) fAntiSigmaCuts->InitCutHistograms(TString("AntiSigma0"));

  if (fV0Cuts && fV0Cuts->GetCutHistograms()) {
    fOutputContainer->Add(fV0Cuts->GetCutHistograms());
  }

  if (fAntiV0Cuts && fAntiV0Cuts->GetCutHistograms()) {
    fOutputContainer->Add(fAntiV0Cuts->GetCutHistograms());
  }

  if (!fIsLightweight) {
    fPhotonQA = new AliSigma0V0Cuts();
    fPhotonQA->SetLightweight(false);
    fPhotonQA->SetPID(22);
    fPhotonQA->SetPosPID(AliPID::kElectron, 11);
    fPhotonQA->SetNegPID(AliPID::kElectron, -11);
    fPhotonQA->InitCutHistograms(TString("Photon"));
    fOutputContainer->Add(fPhotonQA->GetCutHistograms());
  }

  if (fSigmaCuts && fSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fSigmaCuts->GetCutHistograms());
  }

  if (fAntiSigmaCuts && fAntiSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fAntiSigmaCuts->GetCutHistograms());
  }

  const int nPairs = (fIsLightweight) ? 6 : 10;
  fPairCleaner =
      new AliFemtoDreamPairCleaner(nPairs, 0, fConfig->GetMinimalBookingME());
  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  if (!fConfig->GetMinimalBookingME() && fPairCleaner &&
      fPairCleaner->GetHistList()) {
    fOutputFemto->Add(fPairCleaner->GetHistList());
  }

  if (fPartColl && fPartColl->GetHistList()) {
    fOutputFemto->Add(fPartColl->GetHistList());
  }
  if (!fConfig->GetMinimalBookingME() && fPartColl && fPartColl->GetQAList()) {
    fOutputFemto->Add(fPartColl->GetQAList());
  }

  if (!fIsLightweight) {
    fHistCorrelationPSigmaPLambda[0] =
        new TH2F("fHistCorrelationPSigmaPLambda",
                 "; k* (p-#Sigma^{0}) (GeV/#it{c});  k* "
                 "(p-#Lambda_{#Sigma^{0}}) (GeV/#it{c})",
                 750, 0, 3, 750, 0, 3);
    fHistCorrelationPSigmaPGamma[0] =
        new TH2F("fHistCorrelationPSigmaPGamma",
                 "; k* (p-#Sigma^{0}) (GeV/#it{c});  k* "
                 "(p-#gamma_{#Sigma^{0}}) (GeV/#it{c})",
                 750, 0, 3, 750, 0, 3);
    fHistCorrelationPLambdaPGamma[0] =
        new TH2F("fHistCorrelationPLambdaPGamma",
                 "; k* (p-#Lambda_{#Sigma^{0}}) (GeV/#it{c});  k* "
                 "(p-#gamma_{#Sigma^{0}}) (GeV/#it{c})",
                 750, 0, 3, 750, 0, 3);
    fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[0] =
        new TH2F("fHistCorrelationAntiPAntiSigmaAntiPAntiLambda",
                 "; k* (#bar{p}-#bar{#Sigma}^{0}) (GeV/#it{c});  k* "
                 "(#bar{p}-#bar{#Lambda}_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
                 750, 0, 3, 750, 0, 3);
    fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[0] =
        new TH2F("fHistCorrelationAntiPAntiSigmaAntiPAntiGamma",
                 "; k* (#bar{p}-#bar{#Sigma}^{0}) (GeV/#it{c});  k* "
                 "(#bar{p}-#gamma_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
                 750, 0, 3, 750, 0, 3);
    fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[0] = new TH2F(
        "fHistCorrelationAntiPAntiLambdaAntiPAntiGamma",
        "; k* (#bar{p}-#bar{#Lambda}_{#bar{#Sigma}^{0}}) (GeV/#it{c});  k* "
        "(#bar{p}-#gamma_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
        750, 0, 3, 750, 0, 3);
    fOutputFemto->Add(fHistCorrelationPSigmaPLambda[0]);
    fOutputFemto->Add(fHistCorrelationPSigmaPGamma[0]);
    fOutputFemto->Add(fHistCorrelationPLambdaPGamma[0]);
    fOutputFemto->Add(fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[0]);
    fOutputFemto->Add(fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[0]);
    fOutputFemto->Add(fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[0]);

    if (fIsMC) {
      for (unsigned int i = 1; i < 3; ++i) {
        const char *appendix = (i == 1) ? "SigmaPrim" : "Background";
        fHistCorrelationPSigmaPLambda[i] =
            new TH2F(Form("fHistCorrelationPSigmaPLambda_%s", appendix),
                     "; k* (p-#Sigma^{0}) (GeV/#it{c});  k* "
                     "(p-#Lambda_{#Sigma^{0}}) (GeV/#it{c})",
                     750, 0, 3, 750, 0, 3);
        fHistCorrelationPSigmaPGamma[i] =
            new TH2F(Form("fHistCorrelationPSigmaPGamma_%s", appendix),
                     "; k* (p-#Sigma^{0}) (GeV/#it{c});  k* "
                     "(p-#gamma_{#Sigma^{0}}) (GeV/#it{c})",
                     750, 0, 3, 750, 0, 3);
        fHistCorrelationPLambdaPGamma[i] =
            new TH2F(Form("fHistCorrelationPLambdaPGamma_%s", appendix),
                     "; k* (p-#Lambda_{#Sigma^{0}}) (GeV/#it{c});  k* "
                     "(p-#gamma_{#Sigma^{0}}) (GeV/#it{c})",
                     750, 0, 3, 750, 0, 3);
        fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[i] = new TH2F(
            Form("fHistCorrelationAntiPAntiSigmaAntiPAntiLambda_%s", appendix),
            "; k* (#bar{p}-#bar{#Sigma}^{0}) (GeV/#it{c});  k* "
            "(#bar{p}-#bar{#Lambda}_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
            750, 0, 3, 750, 0, 3);
        fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[i] = new TH2F(
            Form("fHistCorrelationAntiPAntiSigmaAntiPAntiGamma_%s", appendix),
            "; k* (#bar{p}-#bar{#Sigma}^{0}) (GeV/#it{c});  k* "
            "(#bar{p}-#gamma_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
            750, 0, 3, 750, 0, 3);
        fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[i] = new TH2F(
            Form("fHistCorrelationAntiPAntiLambdaAntiPAntiGamma_%s", appendix),
            "; k* (#bar{p}-#bar{#Lambda}_{#bar{#Sigma}^{0}}) "
            "(GeV/#it{c});  k* "
            "(#bar{p}-#gamma_{#bar{#Sigma}^{0}}) (GeV/#it{c})",
            750, 0, 3, 750, 0, 3);
        fOutputFemto->Add(fHistCorrelationPSigmaPLambda[i]);
        fOutputFemto->Add(fHistCorrelationPSigmaPGamma[i]);
        fOutputFemto->Add(fHistCorrelationPLambdaPGamma[i]);
        fOutputFemto->Add(fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[i]);
        fOutputFemto->Add(fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[i]);
        fOutputFemto->Add(fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[i]);
      }
    }
  }

  PostData(1, fOutputContainer);
  PostData(2, fOutputFemto);
}

void AliAnalysisTaskSigma0Femto::FillCorrelationCorrelator(
    const std::vector<AliFemtoDreamBasePart> &particles,
    const std::vector<AliFemtoDreamBasePart> &sigmasFemto,
    const std::vector<AliSigma0ParticlePhotonMother> &sigmas,
    const bool isAnti) const {
  for (const auto &proton : particles) {
    const TVector3 protonVec = proton.GetMomentum();
    if (!proton.UseParticle()) continue;

    unsigned int sigmaCounter = 0;
    for (const auto &sigma : sigmas) {
      const auto &SigmaParticleFemto = sigmasFemto[sigmaCounter++];
      if (!SigmaParticleFemto.UseParticle()) continue;
      const TVector3 sigmaVec =
          TVector3(sigma.GetPx(), sigma.GetPy(), sigma.GetPz());
      const auto lambda = sigma.GetV0();
      const TVector3 lambdaVec =
          TVector3(lambda.GetPx(), lambda.GetPy(), lambda.GetPz());
      const auto photon = sigma.GetPhoton();
      const TVector3 photonVec =
          TVector3(photon.GetPx(), photon.GetPy(), photon.GetPz());

      const float kRelpSigma =
          ComputeRelk(protonVec, 2212, sigmaVec, (isAnti) ? -3212 : 3212);
      const float kRelpLambda =
          ComputeRelk(protonVec, 2212, lambdaVec, (isAnti) ? -3122 : 3122);
      const float kRelpPhoton = ComputeRelk(protonVec, 2212, photonVec, 22);

      if (!isAnti) {
        fHistCorrelationPSigmaPLambda[0]->Fill(kRelpSigma, kRelpLambda);
        fHistCorrelationPSigmaPGamma[0]->Fill(kRelpSigma, kRelpPhoton);
        fHistCorrelationPLambdaPGamma[0]->Fill(kRelpLambda, kRelpPhoton);
      } else {
        fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[0]->Fill(kRelpSigma,
                                                               kRelpLambda);
        fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[0]->Fill(kRelpSigma,
                                                              kRelpPhoton);
        fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[0]->Fill(kRelpLambda,
                                                               kRelpPhoton);
      }

      if (fIsMC) {
        if (SigmaParticleFemto.GetParticleOrigin() ==
            AliFemtoDreamBasePart::kPhysPrimary) {
          if (!isAnti) {
            fHistCorrelationPSigmaPLambda[1]->Fill(kRelpSigma, kRelpLambda);
            fHistCorrelationPSigmaPGamma[1]->Fill(kRelpSigma, kRelpPhoton);
            fHistCorrelationPLambdaPGamma[1]->Fill(kRelpLambda, kRelpPhoton);
          } else {
            fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[1]->Fill(kRelpSigma,
                                                                   kRelpLambda);
            fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[1]->Fill(kRelpSigma,
                                                                  kRelpPhoton);
            fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[1]->Fill(kRelpLambda,
                                                                   kRelpPhoton);
          }
        } else if (SigmaParticleFemto.GetParticleOrigin() ==
                   AliFemtoDreamBasePart::kFake) {
          if (!isAnti) {
            fHistCorrelationPSigmaPLambda[2]->Fill(kRelpSigma, kRelpLambda);
            fHistCorrelationPSigmaPGamma[2]->Fill(kRelpSigma, kRelpPhoton);
            fHistCorrelationPLambdaPGamma[2]->Fill(kRelpLambda, kRelpPhoton);
          } else {
            fHistCorrelationAntiPAntiSigmaAntiPAntiLambda[2]->Fill(kRelpSigma,
                                                                   kRelpLambda);
            fHistCorrelationAntiPAntiSigmaAntiPAntiGamma[2]->Fill(kRelpSigma,
                                                                  kRelpPhoton);
            fHistCorrelationAntiPAntiLambdaAntiPAntiGamma[2]->Fill(kRelpLambda,
                                                                   kRelpPhoton);
          }
        }
      }
    }
  }
}

float AliAnalysisTaskSigma0Femto::ComputeRelk(const TVector3 &Part1Momentum,
                                              const int PDGPart1,
                                              const TVector3 &Part2Momentum,
                                              const int PDGPart2) const {
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  SPtrackCMS = SPtrack;
  TPProngCMS = TPProng;

  SPtrackCMS.Boost(-betax, -betay, -betaz);
  TPProngCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5 * trackRelK.P();
  return results;
}
