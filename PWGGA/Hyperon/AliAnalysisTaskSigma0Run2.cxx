#include "AliAnalysisTaskSigma0Run2.h"

#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliMCEventHandler.h"
#include "AliMagF.h"
#include "AliStack.h"

#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"

ClassImp(AliAnalysisTaskSigma0Run2)

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2()
    : AliAnalysisTaskSE("AliAnalysisTaskSigma0Run2"),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMCStack(nullptr),
      fESDpid(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fQA(nullptr),
      fEventCuts(nullptr),
      fSingleParticleCuts(nullptr),
      fSingleParticleCutsProton(nullptr),
      fV0Cuts(nullptr),
      fV0LambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fPhotonMotherCuts(nullptr),
      fEventContainer(nullptr),
      fElectron(),
      fPositron(),
      fProton(),
      fAntiProton(),
      fLambda(),
      fAntiLambda(),
      fLambdaFemto(),
      fAntiLambdaFemto(),
      fSigma(),
      fAntiSigma(),
      fPhoton(),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsQA(false),
      fOutputContainer() {}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2(const char *name)
    : AliAnalysisTaskSE(name),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMCStack(nullptr),
      fESDpid(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fQA(nullptr),
      fEventCuts(nullptr),
      fSingleParticleCuts(nullptr),
      fSingleParticleCutsProton(nullptr),
      fV0Cuts(nullptr),
      fV0LambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fPhotonMotherCuts(nullptr),
      fEventContainer(nullptr),
      fElectron(),
      fPositron(),
      fProton(),
      fAntiProton(),
      fLambda(),
      fAntiLambda(),
      fLambdaFemto(),
      fAntiLambdaFemto(),
      fSigma(),
      fAntiSigma(),
      fPhoton(),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsQA(false),
      fOutputContainer() {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::UserExec(Option_t * /*option*/) {
  AliVEvent *fInputEvent = InputEvent();
  if (fIsMC) fMCEvent = MCEvent();

  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }

  if (!fEventCuts) {
    AliWarning("No Event Cuts");
    return;
  }
  if (!fEventCuts->EventIsSelected(fInputEvent, fMCEvent)) {
    PostData(1, fOutputContainer);
    return;
  }

  // Electron Selection
  if (fSingleParticleCuts)
    fSingleParticleCuts->SelectSingleParticles(fInputEvent, fMCEvent, fPositron,
                                               fElectron);


  // Proton Selection
  if (fSingleParticleCutsProton)
    fSingleParticleCutsProton->SelectSingleParticles(fInputEvent, fMCEvent, fProton,
                                               fAntiProton);

  // V0 selection for Sigma0
  if (fV0Cuts)
    fV0Cuts->SelectV0s(fInputEvent, fMCEvent, fLambda, fAntiLambda,
                       AliPID::kProton, AliPID::kPion);

  // V0 selection for Lambda
  std::vector<AliSigma0ParticleV0> lambda;
  std::vector<AliSigma0ParticleV0> antiLambda;
  if (fV0LambdaCuts)
    fV0LambdaCuts->SelectV0s(fInputEvent, fMCEvent, lambda, antiLambda,
                             AliPID::kProton, AliPID::kPion);

  // Photon selection
  if (fPhotonCuts) fPhotonCuts->SelectPhotons(fInputEvent, fMCEvent, fPhoton);

  // Photon mother selection
  if (fPhotonMotherCuts)
    fPhotonMotherCuts->SelectPhotonMother(fInputEvent, fMCEvent, fPhoton,
                                          fLambda, fAntiLambda, fSigma,
                                          fAntiSigma);

  // Event container
  if (fEventContainer) {
    fEventContainer->TrackCleaner(fProton, fAntiProton, fLambda, fAntiLambda);
    fEventContainer->ProcessEvent(fInputEvent, fMCEvent, fProton, fAntiProton, fLambda, fAntiLambda, fSigma, fAntiSigma, fPhoton);
  }

  // flush the data
  PostData(1, fOutputContainer);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::UserCreateOutputObjects() {
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
    fOutputContainer = nullptr;
  }
  if (fOutputContainer == nullptr) {
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  // Event selection
  if (fEventCuts && fEventCuts->GetCutHistograms()) {
    fOutputContainer->Add(fEventCuts->GetCutHistograms());
  }

  // Single particle selection
  if (fSingleParticleCuts && fSingleParticleCuts->GetCutHistograms()) {
    fOutputContainer->Add(fSingleParticleCuts->GetCutHistograms());
  }

  // Single particle selection
  if (fSingleParticleCutsProton && fSingleParticleCutsProton->GetCutHistograms()) {
    fOutputContainer->Add(fSingleParticleCutsProton->GetCutHistograms());
  }

  // V0 selection for Sigma0
  if (fV0Cuts && fV0Cuts->GetCutHistograms()) {
    fOutputContainer->Add(fV0Cuts->GetCutHistograms());
  }

  // V0 selection for Lambda
  if (fV0LambdaCuts && fV0LambdaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fV0LambdaCuts->GetCutHistograms());
  }

  // Photon selection
  if (fPhotonCuts && fPhotonCuts->GetCutHistograms()) {
    fOutputContainer->Add(fPhotonCuts->GetCutHistograms());
  }

  // Photon mother selection
  if (fPhotonMotherCuts && fPhotonMotherCuts->GetCutHistograms()) {
    fOutputContainer->Add(fPhotonMotherCuts->GetCutHistograms());
  }

  // Event container
  if (fEventContainer && fEventContainer->GetCutHistograms()) {
    fOutputContainer->Add(fEventContainer->GetCutHistograms());
  }

  PostData(1, fOutputContainer);
}
