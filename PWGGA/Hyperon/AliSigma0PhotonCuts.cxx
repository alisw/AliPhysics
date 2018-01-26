#include "AliSigma0PhotonCuts.h"

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"

#include <iostream>

ClassImp(AliSigma0PhotonCuts)

//____________________________________________________________________________________________________
AliSigma0PhotonCuts::AliSigma0PhotonCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fReaderGammas(nullptr),
      fIsMC(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMCStack(nullptr),
      fDCAr(999.f),
      fDCAz(999.f),
      fCosPA(0.f),
      fPmax(999.f),
      fPtmax(999.f),
      fPtElemax(999.f),
      fPIDResponse(nullptr),
      fHistCuts(nullptr),
      fHistPhotonCuts(nullptr),
      fHistNPhoton(nullptr),
      fHistPhotonPt(nullptr),
      fHistPhotonP(nullptr),
      fHistPhotonInvMass(nullptr),
      fHistPhotonInvMassPt(nullptr),
      fHistPhotonInvMassEta(nullptr),
      fHistPhotonCPAPt(nullptr),
      fHistPhotonR(nullptr),
      fHistPhotonArm(nullptr),
      fHistPhotonDCAz(nullptr),
      fHistPhotonDCAr(nullptr),
      fHistPhotonEtaPhi(nullptr),
      fHistPhotonConvPointX(nullptr),
      fHistPhotonConvPointY(nullptr),
      fHistPhotonConvPointZ(nullptr),
      fHistPhotonEleP(nullptr),
      fHistPhotonElePt(nullptr),
      fHistPhotonEleNsigmaTPC(nullptr),
      fHistPhotonEleNsigmaTPCPion(nullptr),
      fHistPhotonEleTPCsignal(nullptr),
      fHistTwoPhotonPt(nullptr),
      fHistTwoPhotonP(nullptr),
      fHistTwoPhotonInvMass(nullptr),
      fHistTwoPhotonInvMassPt(nullptr),
      fHistTwoPhotonInvMassEta(nullptr),
      fHistMCRecSigma0PhotonPt(nullptr),
      fHistMCRecSigma0PhotonP(nullptr),
      fHistMCRecSigma0PhotonInvMass(nullptr),
      fHistMCRecSigma0PhotonInvMassPt(nullptr),
      fHistMCRecSigma0PhotonInvMassEta(nullptr),
      fHistMCRecSigma0PhotonR(nullptr),
      fHistMCRecSigma0PhotonArm(nullptr),
      fHistMCRecSigma0PhotonDCAz(nullptr),
      fHistMCRecSigma0PhotonDCAr(nullptr),
      fHistMCRecSigma0PhotonConvPointX(nullptr),
      fHistMCRecSigma0PhotonConvPointY(nullptr),
      fHistMCRecSigma0PhotonConvPointZ(nullptr),
      fHistMCRecSigma0PhotonEleP(nullptr),
      fHistMCRecSigma0PhotonElePt(nullptr),
      fHistMCRecSigma0CPAPt(nullptr),
      fHistMCRecSigma0PhotonDCAzPt(nullptr),
      fHistMCRecSigma0PhotonDCArPt(nullptr),
      fHistMCRecSigma0PhotonPsiPairPt(nullptr),
      fHistMCRecSigma0PhotonChi2Pt(nullptr),
      fHistMCFakeSigma0PhotonPt(nullptr),
      fHistMCFakeSigma0PhotonP(nullptr),
      fHistMCFakeSigma0PhotonInvMass(nullptr),
      fHistMCFakeSigma0PhotonInvMassPt(nullptr),
      fHistMCFakeSigma0PhotonInvMassEta(nullptr),
      fHistMCFakeSigma0PhotonR(nullptr),
      fHistMCFakeSigma0PhotonArm(nullptr),
      fHistMCFakeSigma0PhotonDCAz(nullptr),
      fHistMCFakeSigma0PhotonDCAr(nullptr),
      fHistMCFakeSigma0PhotonConvPointX(nullptr),
      fHistMCFakeSigma0PhotonConvPointY(nullptr),
      fHistMCFakeSigma0PhotonConvPointZ(nullptr),
      fHistMCFakeSigma0PhotonEleP(nullptr),
      fHistMCFakeSigma0PhotonElePt(nullptr),
      fHistMCFakeSigma0CPAPt(nullptr),
      fHistMCFakeSigma0PhotonDCAzPt(nullptr),
      fHistMCFakeSigma0PhotonDCArPt(nullptr),
      fHistMCFakeSigma0PhotonPsiPairPt(nullptr),
      fHistMCFakeSigma0PhotonChi2Pt(nullptr),
      fHistMCFakePhotonPt(nullptr),
      fHistMCFakePhotonP(nullptr),
      fHistMCFakePhotonInvMass(nullptr),
      fHistMCFakePhotonInvMassPt(nullptr),
      fHistMCFakePhotonInvMassEta(nullptr),
      fHistMCFakePhotonR(nullptr),
      fHistMCFakePhotonArm(nullptr),
      fHistMCFakePhotonDCAz(nullptr),
      fHistMCFakePhotonDCAr(nullptr),
      fHistMCFakePhotonConvPointX(nullptr),
      fHistMCFakePhotonConvPointY(nullptr),
      fHistMCFakePhotonConvPointZ(nullptr),
      fHistMCFakePhotonEleP(nullptr),
      fHistMCFakePhotonElePt(nullptr),
      fHistMCFakePhotonCPAPt(nullptr),
      fHistMCFakePhotonDCAzPt(nullptr),
      fHistMCFakePhotonDCArPt(nullptr),
      fHistMCFakePhotonPsiPairPt(nullptr),
      fHistMCFakePhotonChi2Pt(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts::AliSigma0PhotonCuts(const AliSigma0PhotonCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fReaderGammas(nullptr),
      fIsMC(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fMCStack(nullptr),
      fDCAr(999.f),
      fDCAz(999.f),
      fCosPA(0.f),
      fPmax(999.f),
      fPtmax(999.f),
      fPtElemax(999.f),
      fPIDResponse(nullptr),
      fHistCuts(nullptr),
      fHistPhotonCuts(nullptr),
      fHistNPhoton(nullptr),
      fHistPhotonPt(nullptr),
      fHistPhotonP(nullptr),
      fHistPhotonInvMass(nullptr),
      fHistPhotonInvMassPt(nullptr),
      fHistPhotonInvMassEta(nullptr),
      fHistPhotonCPAPt(nullptr),
      fHistPhotonR(nullptr),
      fHistPhotonArm(nullptr),
      fHistPhotonDCAz(nullptr),
      fHistPhotonDCAr(nullptr),
      fHistPhotonEtaPhi(nullptr),
      fHistPhotonConvPointX(nullptr),
      fHistPhotonConvPointY(nullptr),
      fHistPhotonConvPointZ(nullptr),
      fHistPhotonEleP(nullptr),
      fHistPhotonElePt(nullptr),
      fHistPhotonEleNsigmaTPC(nullptr),
      fHistPhotonEleNsigmaTPCPion(nullptr),
      fHistPhotonEleTPCsignal(nullptr),
      fHistTwoPhotonPt(nullptr),
      fHistTwoPhotonP(nullptr),
      fHistTwoPhotonInvMass(nullptr),
      fHistTwoPhotonInvMassPt(nullptr),
      fHistTwoPhotonInvMassEta(nullptr),
      fHistMCRecSigma0PhotonPt(nullptr),
      fHistMCRecSigma0PhotonP(nullptr),
      fHistMCRecSigma0PhotonInvMass(nullptr),
      fHistMCRecSigma0PhotonInvMassPt(nullptr),
      fHistMCRecSigma0PhotonInvMassEta(nullptr),
      fHistMCRecSigma0PhotonR(nullptr),
      fHistMCRecSigma0PhotonArm(nullptr),
      fHistMCRecSigma0PhotonDCAz(nullptr),
      fHistMCRecSigma0PhotonDCAr(nullptr),
      fHistMCRecSigma0PhotonConvPointX(nullptr),
      fHistMCRecSigma0PhotonConvPointY(nullptr),
      fHistMCRecSigma0PhotonConvPointZ(nullptr),
      fHistMCRecSigma0PhotonEleP(nullptr),
      fHistMCRecSigma0PhotonElePt(nullptr),
      fHistMCRecSigma0CPAPt(nullptr),
      fHistMCRecSigma0PhotonDCAzPt(nullptr),
      fHistMCRecSigma0PhotonDCArPt(nullptr),
      fHistMCRecSigma0PhotonPsiPairPt(nullptr),
      fHistMCRecSigma0PhotonChi2Pt(nullptr),
      fHistMCFakeSigma0PhotonPt(nullptr),
      fHistMCFakeSigma0PhotonP(nullptr),
      fHistMCFakeSigma0PhotonInvMass(nullptr),
      fHistMCFakeSigma0PhotonInvMassPt(nullptr),
      fHistMCFakeSigma0PhotonInvMassEta(nullptr),
      fHistMCFakeSigma0PhotonR(nullptr),
      fHistMCFakeSigma0PhotonArm(nullptr),
      fHistMCFakeSigma0PhotonDCAz(nullptr),
      fHistMCFakeSigma0PhotonDCAr(nullptr),
      fHistMCFakeSigma0PhotonConvPointX(nullptr),
      fHistMCFakeSigma0PhotonConvPointY(nullptr),
      fHistMCFakeSigma0PhotonConvPointZ(nullptr),
      fHistMCFakeSigma0PhotonEleP(nullptr),
      fHistMCFakeSigma0PhotonElePt(nullptr),
      fHistMCFakeSigma0CPAPt(nullptr),
      fHistMCFakeSigma0PhotonDCAzPt(nullptr),
      fHistMCFakeSigma0PhotonDCArPt(nullptr),
      fHistMCFakeSigma0PhotonPsiPairPt(nullptr),
      fHistMCFakeSigma0PhotonChi2Pt(nullptr),
      fHistMCFakePhotonPt(nullptr),
      fHistMCFakePhotonP(nullptr),
      fHistMCFakePhotonInvMass(nullptr),
      fHistMCFakePhotonInvMassPt(nullptr),
      fHistMCFakePhotonInvMassEta(nullptr),
      fHistMCFakePhotonR(nullptr),
      fHistMCFakePhotonArm(nullptr),
      fHistMCFakePhotonDCAz(nullptr),
      fHistMCFakePhotonDCAr(nullptr),
      fHistMCFakePhotonConvPointX(nullptr),
      fHistMCFakePhotonConvPointY(nullptr),
      fHistMCFakePhotonConvPointZ(nullptr),
      fHistMCFakePhotonEleP(nullptr),
      fHistMCFakePhotonElePt(nullptr),
      fHistMCFakePhotonCPAPt(nullptr),
      fHistMCFakePhotonDCAzPt(nullptr),
      fHistMCFakePhotonDCArPt(nullptr),
      fHistMCFakePhotonPsiPairPt(nullptr),
      fHistMCFakePhotonChi2Pt(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts &AliSigma0PhotonCuts::operator=(
    const AliSigma0PhotonCuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts::~AliSigma0PhotonCuts() {
  if(fPIDResponse) delete fPIDResponse;
  if(fV0Reader) delete fV0Reader;
  if(fReaderGammas) delete fReaderGammas;

}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts *AliSigma0PhotonCuts::DefaultCuts() {
  AliSigma0PhotonCuts *photonCuts = new AliSigma0PhotonCuts();
  photonCuts->SetDCArMax(1E30);
  photonCuts->SetDCAzMax(1E30);
  photonCuts->SetCosPA(-10.f);
  photonCuts->SetPmax(1E30);
  photonCuts->SetPtmax(1E30);
  photonCuts->SetPtElemax(1E30);
  return photonCuts;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonCuts::SelectPhotons(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliAODConversionPhoton> &photon) {
  photon.clear();

  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  // should no longer be used, but the corresponding functions have not been
  // changed in the Photon classes yet
  if (fIsMC) fMCStack = fMCEvent->Stack();

  fInputEvent = inputEvent;

  fV0Reader = (AliV0ReaderV1 *)man->GetTask(fV0ReaderName.Data());
  if (!fV0Reader) {
    std::cout << "Error: No V0 Reader " << fV0ReaderName.Data() << "\n";
    return;
  }

  fReaderGammas =
      fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  fHistNPhoton->Fill(fReaderGammas->GetEntriesFast());

  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();

  for (Int_t iGamma = 0; iGamma < fReaderGammas->GetEntriesFast(); ++iGamma) {
    AliAODConversionPhoton *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fReaderGammas->At(iGamma));
    if (!PhotonCandidate) continue;
    fHistPhotonP->Fill(PhotonCandidate->GetPhotonP());
    fHistPhotonPt->Fill(PhotonCandidate->GetPhotonPt());
    fHistPhotonInvMass->Fill(PhotonCandidate->GetInvMassPair());
    fHistPhotonInvMassPt->Fill(PhotonCandidate->GetPhotonPt(),
                               PhotonCandidate->GetInvMassPair());
    fHistPhotonInvMassEta->Fill(PhotonCandidate->Eta(),
                                PhotonCandidate->GetInvMassPair());
    fHistPhotonR->Fill(PhotonCandidate->GetPhotonPt(),
                       PhotonCandidate->GetConversionRadius());
    fHistPhotonArm->Fill(PhotonCandidate->GetArmenterosAlpha(),
                         PhotonCandidate->GetArmenterosQt());
    PhotonCandidate->CalculateDistanceOfClossetApproachToPrimVtx(vertex);
    fHistPhotonDCAz->Fill(PhotonCandidate->GetDCAzToPrimVtx());
    fHistPhotonDCAr->Fill(PhotonCandidate->GetDCArToPrimVtx());
    fHistPhotonCPAPt->Fill(
        PhotonCandidate->GetPhotonPt(),
        getCosineOfPointingAngle(PhotonCandidate, fInputEvent));

    Double_t convpoint[3];
    PhotonCandidate->GetConversionPoint(convpoint);

    fHistPhotonConvPointX->Fill(convpoint[0]);
    fHistPhotonConvPointY->Fill(convpoint[1]);
    fHistPhotonConvPointZ->Fill(convpoint[2]);

    AliVParticle *posEle =
        fInputEvent->GetTrack(PhotonCandidate->GetTrackLabelPositive());
    AliVParticle *negEle =
        fInputEvent->GetTrack(PhotonCandidate->GetTrackLabelNegative());

    AliVEvent *event = static_cast<AliVEvent *>(fInputEvent);
    AliVTrack *posEleTrack = static_cast<AliVTrack *>(
        static_cast<AliConversionPhotonCuts *>(fV0Reader->GetConversionCuts())
            ->GetTrack(event, PhotonCandidate->GetTrackLabelPositive()));
    AliVTrack *negEleTrack = static_cast<AliVTrack *>(
        static_cast<AliConversionPhotonCuts *>(fV0Reader->GetConversionCuts())
            ->GetTrack(event, PhotonCandidate->GetTrackLabelNegative()));

    fHistPhotonEleP->Fill(posEle->P());
    fHistPhotonEleP->Fill(negEle->P());
    fHistPhotonElePt->Fill(posEle->Pt());
    fHistPhotonElePt->Fill(negEle->Pt());

    fHistPhotonEleNsigmaTPC->Fill(
        posEleTrack->P(),
        fPIDResponse->NumberOfSigmasTPC(posEleTrack, AliPID::kElectron));
    fHistPhotonEleNsigmaTPC->Fill(
        negEleTrack->P(),
        fPIDResponse->NumberOfSigmasTPC(negEleTrack, AliPID::kElectron));
    fHistPhotonEleNsigmaTPCPion->Fill(
        posEleTrack->P(),
        fPIDResponse->NumberOfSigmasTPC(posEleTrack, AliPID::kPion));
    fHistPhotonEleNsigmaTPCPion->Fill(
        negEleTrack->P(),
        fPIDResponse->NumberOfSigmasTPC(negEleTrack, AliPID::kPion));
    fHistPhotonEleTPCsignal->Fill(posEleTrack->P(),
                                  posEleTrack->GetTPCsignal());
    fHistPhotonEleTPCsignal->Fill(negEleTrack->P(),
                                  negEleTrack->GetTPCsignal());

    if (fIsMC) {
      AliMCParticle *photonElePos = static_cast<AliMCParticle *>(
          fMCEvent->GetTrack(PhotonCandidate->GetMCLabelPositive()));
      AliMCParticle *photonEleNeg = static_cast<AliMCParticle *>(
          fMCEvent->GetTrack(PhotonCandidate->GetMCLabelNegative()));
      if (photonElePos->GetMother() != photonEleNeg->GetMother()) continue;
      AliMCParticle *photonMC = static_cast<AliMCParticle *>(
          fMCEvent->GetTrack(photonEleNeg->GetMother()));
      if (!photonMC) continue;
      AliMCParticle *photonMother = static_cast<AliMCParticle *>(
          fMCEvent->GetTrack(photonMC->GetMother()));
      if (!photonMother) continue;
      // take only true photons
      if (IsPhotonCandidate(PhotonCandidate)) {
        // mother is Sigma0
        if (std::abs(photonMother->PdgCode()) == 3212) {
          fHistMCRecSigma0PhotonP->Fill(PhotonCandidate->GetPhotonP());
          fHistMCRecSigma0PhotonPt->Fill(PhotonCandidate->GetPhotonPt());
          fHistMCRecSigma0PhotonInvMass->Fill(
              PhotonCandidate->GetInvMassPair());
          fHistMCRecSigma0PhotonInvMassPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetInvMassPair());
          fHistMCRecSigma0PhotonInvMassEta->Fill(
              PhotonCandidate->Eta(), PhotonCandidate->GetInvMassPair());
          fHistMCRecSigma0PhotonR->Fill(PhotonCandidate->GetPhotonPt(),
                                        PhotonCandidate->GetConversionRadius());
          fHistMCRecSigma0PhotonArm->Fill(PhotonCandidate->GetArmenterosAlpha(),
                                          PhotonCandidate->GetArmenterosQt());
          fHistMCRecSigma0PhotonDCAz->Fill(PhotonCandidate->GetDCAzToPrimVtx());
          fHistMCRecSigma0PhotonDCAr->Fill(PhotonCandidate->GetDCArToPrimVtx());
          fHistMCRecSigma0PhotonConvPointX->Fill(convpoint[0]);
          fHistMCRecSigma0PhotonConvPointY->Fill(convpoint[1]);
          fHistMCRecSigma0PhotonConvPointZ->Fill(convpoint[2]);
          fHistMCRecSigma0PhotonEleP->Fill(posEle->P());
          fHistMCRecSigma0PhotonEleP->Fill(negEle->P());
          fHistMCRecSigma0PhotonElePt->Fill(posEle->Pt());
          fHistMCRecSigma0PhotonElePt->Fill(negEle->Pt());
          fHistMCRecSigma0CPAPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              getCosineOfPointingAngle(PhotonCandidate, fInputEvent));
          fHistMCRecSigma0PhotonDCAzPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetDCAzToPrimVtx());
          fHistMCRecSigma0PhotonDCArPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetDCArToPrimVtx());
          fHistMCRecSigma0PhotonPsiPairPt->Fill(PhotonCandidate->GetPhotonPt(),
                                                PhotonCandidate->GetPsiPair());
          fHistMCRecSigma0PhotonChi2Pt->Fill(PhotonCandidate->GetPhotonPt(),
                                             PhotonCandidate->GetChi2perNDF());
        } else {
          fHistMCFakeSigma0PhotonP->Fill(PhotonCandidate->GetPhotonP());
          fHistMCFakeSigma0PhotonPt->Fill(PhotonCandidate->GetPhotonPt());
          fHistMCFakeSigma0PhotonInvMass->Fill(
              PhotonCandidate->GetInvMassPair());
          fHistMCFakeSigma0PhotonInvMassPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetInvMassPair());
          fHistMCFakeSigma0PhotonInvMassEta->Fill(
              PhotonCandidate->Eta(), PhotonCandidate->GetInvMassPair());
          fHistMCFakeSigma0PhotonR->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetConversionRadius());
          fHistMCFakeSigma0PhotonArm->Fill(
              PhotonCandidate->GetArmenterosAlpha(),
              PhotonCandidate->GetArmenterosQt());
          fHistMCFakeSigma0PhotonDCAz->Fill(
              PhotonCandidate->GetDCAzToPrimVtx());
          fHistMCFakeSigma0PhotonDCAr->Fill(
              PhotonCandidate->GetDCArToPrimVtx());
          fHistMCFakeSigma0PhotonConvPointX->Fill(convpoint[0]);
          fHistMCFakeSigma0PhotonConvPointY->Fill(convpoint[1]);
          fHistMCFakeSigma0PhotonConvPointZ->Fill(convpoint[2]);
          fHistMCFakeSigma0PhotonEleP->Fill(posEle->P());
          fHistMCFakeSigma0PhotonEleP->Fill(negEle->P());
          fHistMCFakeSigma0PhotonElePt->Fill(posEle->Pt());
          fHistMCFakeSigma0PhotonElePt->Fill(negEle->Pt());
          fHistMCFakeSigma0CPAPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              getCosineOfPointingAngle(PhotonCandidate, fInputEvent));
          fHistMCFakeSigma0PhotonDCAzPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetDCAzToPrimVtx());
          fHistMCFakeSigma0PhotonDCArPt->Fill(
              PhotonCandidate->GetPhotonPt(),
              PhotonCandidate->GetDCArToPrimVtx());
          fHistMCFakeSigma0PhotonPsiPairPt->Fill(PhotonCandidate->GetPhotonPt(),
                                                 PhotonCandidate->GetPsiPair());
          fHistMCFakeSigma0PhotonChi2Pt->Fill(PhotonCandidate->GetPhotonPt(),
                                              PhotonCandidate->GetChi2perNDF());
        }
      }
      // check what fake photons look like
      else {
        fHistMCFakePhotonP->Fill(PhotonCandidate->GetPhotonP());
        fHistMCFakePhotonPt->Fill(PhotonCandidate->GetPhotonPt());
        fHistMCFakePhotonInvMass->Fill(PhotonCandidate->GetInvMassPair());
        fHistMCFakePhotonInvMassPt->Fill(PhotonCandidate->GetPhotonPt(),
                                         PhotonCandidate->GetInvMassPair());
        fHistMCFakePhotonInvMassEta->Fill(PhotonCandidate->Eta(),
                                          PhotonCandidate->GetInvMassPair());
        fHistMCFakePhotonR->Fill(PhotonCandidate->GetPhotonPt(),
                                 PhotonCandidate->GetConversionRadius());
        fHistMCFakePhotonArm->Fill(PhotonCandidate->GetArmenterosAlpha(),
                                   PhotonCandidate->GetArmenterosQt());
        fHistMCFakePhotonDCAz->Fill(PhotonCandidate->GetDCAzToPrimVtx());
        fHistMCFakePhotonDCAr->Fill(PhotonCandidate->GetDCArToPrimVtx());
        fHistMCFakePhotonConvPointX->Fill(convpoint[0]);
        fHistMCFakePhotonConvPointY->Fill(convpoint[1]);
        fHistMCFakePhotonConvPointZ->Fill(convpoint[2]);
        fHistMCFakePhotonEleP->Fill(posEle->P());
        fHistMCFakePhotonEleP->Fill(negEle->P());
        fHistMCFakePhotonElePt->Fill(posEle->Pt());
        fHistMCFakePhotonElePt->Fill(negEle->Pt());
        fHistMCFakePhotonCPAPt->Fill(
            PhotonCandidate->GetPhotonPt(),
            getCosineOfPointingAngle(PhotonCandidate, fInputEvent));
        fHistMCFakePhotonDCAzPt->Fill(PhotonCandidate->GetPhotonPt(),
                                      PhotonCandidate->GetDCAzToPrimVtx());
        fHistMCFakePhotonDCArPt->Fill(PhotonCandidate->GetPhotonPt(),
                                      PhotonCandidate->GetDCArToPrimVtx());
        fHistMCFakePhotonPsiPairPt->Fill(PhotonCandidate->GetPhotonPt(),
                                         PhotonCandidate->GetPsiPair());
        fHistMCFakePhotonChi2Pt->Fill(PhotonCandidate->GetPhotonPt(),
                                      PhotonCandidate->GetChi2perNDF());
      }
    }
  }

  for (int i = 0; i < fReaderGammas->GetEntriesFast(); ++i) {
    AliAODConversionPhoton *photonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fReaderGammas->At(i));
    if (!photonCandidate) continue;
    fHistPhotonCuts->Fill(0);

    // cuts
    photonCandidate->CalculateDistanceOfClossetApproachToPrimVtx(vertex);
    if (photonCandidate->GetDCAzToPrimVtx() > fDCAz) continue;
    fHistPhotonCuts->Fill(1);
    if (photonCandidate->GetDCArToPrimVtx() > fDCAr) continue;
    fHistPhotonCuts->Fill(2);
    if (getCosineOfPointingAngle(photonCandidate, fInputEvent) < fCosPA)
      continue;
    fHistPhotonCuts->Fill(3);
    if (photonCandidate->GetPhotonP() > fPmax) continue;
    fHistPhotonCuts->Fill(4);
    if (photonCandidate->GetPhotonPt() > fPtmax) continue;
    fHistPhotonCuts->Fill(5);
    AliVEvent *aodEvent = static_cast<AliVEvent *>(fInputEvent);
    AliVTrack *posEleTrack = static_cast<AliVTrack *>(
        static_cast<AliConversionPhotonCuts *>(fV0Reader->GetConversionCuts())
            ->GetTrack(aodEvent, photonCandidate->GetTrackLabelPositive()));
    AliVTrack *negEleTrack = static_cast<AliVTrack *>(
        static_cast<AliConversionPhotonCuts *>(fV0Reader->GetConversionCuts())
            ->GetTrack(aodEvent, photonCandidate->GetTrackLabelNegative()));
    if (negEleTrack->Pt() > fPtElemax || posEleTrack->Pt() > fPtElemax)
      continue;
    fHistPhotonCuts->Fill(6);
    fHistPhotonEtaPhi->Fill(photonCandidate->Eta(), photonCandidate->Phi());
    photon.push_back(*photonCandidate);
  }

  // Two photon combined histos
  for (auto photon1 = photon.begin(); photon1 != photon.end(); ++photon1) {
    TLorentzVector photonVector1;
    photonVector1.SetXYZM(photon1->Px(), photon1->Py(), photon1->Pz(),
                          photon1->GetInvMassPair());
    for (auto photon2 = photon1 + 1; photon2 != photon.end(); ++photon2) {
      TLorentzVector photonVector2, photonSum;
      photonVector2.SetXYZM(photon2->Px(), photon2->Py(), photon2->Pz(),
                            photon2->GetInvMassPair());
      photonSum = photonVector1 + photonVector2;
      fHistTwoPhotonPt->Fill(photonSum.Pt());
      fHistTwoPhotonP->Fill(photonSum.P());
      fHistTwoPhotonInvMass->Fill(photonSum.M());
      fHistTwoPhotonInvMassPt->Fill(photonSum.Pt(), photonSum.M());
      fHistTwoPhotonInvMassEta->Fill(photonSum.Pt(), photonSum.Eta());
    }
  }
}

//____________________________________________________________________________________________________
double AliSigma0PhotonCuts::getCosineOfPointingAngle(
    const AliConversionPhotonBase *photon, AliVEvent *event) const {
  // calculates the pointing angle of the recalculated V0

  Double_t momV0[3] = {0, 0, 0};
  if (event->IsA() == AliESDEvent::Class()) {
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(event);
    if (!esdEvent) return -999;
    AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());
    if (!v0) return -999;
    v0->GetPxPyPz(momV0[0], momV0[1], momV0[2]);
  }
  if (event->IsA() == AliAODEvent::Class()) {
    momV0[0] = photon->GetPx();
    momV0[1] = photon->GetPy();
    momV0[2] = photon->GetPz();
  }

  // Recalculated V0 Position vector
  double PosV0[3] = {
      photon->GetConversionX() - event->GetPrimaryVertex()->GetX(),
      photon->GetConversionY() - event->GetPrimaryVertex()->GetY(),
      photon->GetConversionZ() - event->GetPrimaryVertex()->GetZ()};

  double momV02 =
      momV0[0] * momV0[0] + momV0[1] * momV0[1] + momV0[2] * momV0[2];
  double PosV02 =
      PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1] + PosV0[2] * PosV0[2];

  double cosinePointingAngle = -999;
  if (momV02 * PosV02 > 0.0) {
    cosinePointingAngle =
        (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2]) /
        std::sqrt(momV02 * PosV02);
  }

  return cosinePointingAngle;
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonCuts::IsPhotonCandidate(
    AliAODConversionPhoton *TruePhotonCandidate) const {
  AliMCParticle *posDaughter = static_cast<AliMCParticle *>(
      fMCEvent->GetTrack(TruePhotonCandidate->GetMCLabelPositive()));
  AliMCParticle *negDaughter = static_cast<AliMCParticle *>(
      fMCEvent->GetTrack(TruePhotonCandidate->GetMCLabelNegative()));

  if (posDaughter == nullptr || negDaughter == nullptr)
    return false;  // One particle does not exist
  int pdgCode[2] = {std::abs(posDaughter->PdgCode()),
                    std::abs(negDaughter->PdgCode())};

  if (posDaughter->GetMother() != negDaughter->GetMother())
    return false;

  else if (posDaughter->GetMother() == -1)
    return false;

  if (pdgCode[0] != 11 || pdgCode[1] != 11)
    return false;  // One Particle is not a electron

  if (posDaughter->PdgCode() == negDaughter->PdgCode())
    return false;  // Same Charge

  AliMCParticle *Photon = static_cast<AliMCParticle *>(
      fMCEvent->GetTrack(posDaughter->GetMother()));
  if (Photon->PdgCode() != 22) return false;  // Mother is no Photon

  if (((posDaughter->GetUniqueID())) != 5 ||
      ((negDaughter->GetUniqueID())) != 5)
    return false;  // check if the daughters come from a conversion

  return true;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonCuts::InitCutHistograms() {
  std::cout << "============================\n"
            << " PHOTON CUT CONFIGURATION \n"
            << " CPA min    " << fCosPA << "\n"
            << " DCA r max  " << fDCAr << "\n"
            << " DCA z max  " << fDCAz << "\n"
            << " P max      " << fPmax << "\n"
            << " Pt max     " << fPtmax << "\n"
            << " Pt max ele " << fPtElemax << "\n"
            << "============================\n";
  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName("PhotonCut_QA");
  }

  // V0Reader QA
  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    std::cout << "Error: No V0 Reader \n";
    return;
  }

  if (fV0Reader) {
    if ((AliConversionPhotonCuts *)fV0Reader->GetConversionCuts() &&
        ((AliConversionPhotonCuts *)fV0Reader->GetConversionCuts())
            ->GetCutHistograms()) {
      fHistograms->Add(
          ((AliConversionPhotonCuts *)fV0Reader->GetConversionCuts())
              ->GetCutHistograms());
    }
    if ((AliConvEventCuts *)fV0Reader->GetEventCuts() &&
        ((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetCutHistograms()) {
      fHistograms->Add(
          ((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetCutHistograms());
    }
    if (fV0Reader->GetProduceV0FindingEfficiency() &&
        fV0Reader->GetV0FindingEfficiencyHistograms()) {
      fHistograms->Add(fV0Reader->GetV0FindingEfficiencyHistograms());
    }
    if (fV0Reader->GetProduceImpactParamHistograms()) {
      fHistograms->Add(fV0Reader->GetImpactParamHistograms());
    }
  }

  fHistCuts = new TProfile("fHistCuts", ";;Cut value", 10, 0, 10);
  fHistCuts->GetXaxis()->SetBinLabel(1, "DCAr max");
  fHistCuts->GetXaxis()->SetBinLabel(2, "DCAz max");
  fHistCuts->GetXaxis()->SetBinLabel(3, "cos#alpha min");
  fHistCuts->GetXaxis()->SetBinLabel(4, "#it{p} max");
  fHistCuts->GetXaxis()->SetBinLabel(5, "#it{p}_{T, #gamma} max");
  fHistCuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T, e^{-}} max");
  fHistograms->Add(fHistCuts);

  fHistCuts->Fill(0.f, fDCAr);
  fHistCuts->Fill(1.f, fDCAz);
  fHistCuts->Fill(2.f, fCosPA);
  fHistCuts->Fill(3.f, fPmax);
  fHistCuts->Fill(4.f, fPtmax);
  fHistCuts->Fill(5.f, fPtElemax);

  fHistPhotonCuts = new TH1F("fHistPhotonCuts", ";;Entries", 10, 0, 10);
  fHistPhotonCuts->GetXaxis()->SetBinLabel(1, "#gamma");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(2, "DCAr");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(3, "DCAz");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(4, "cos#alpha");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(5, "#it{p}");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T, #gamma} max");
  fHistPhotonCuts->GetXaxis()->SetBinLabel(7, "#it{p}_{T, e^{-}} max");
  fHistograms->Add(fHistPhotonCuts);

  fHistNPhoton =
      new TH1F("fHistNPhoton", "; # #gamma per event; Entries", 50, 0, 50);
  fHistPhotonPt = new TH1F("fHistPhotonPt",
                           "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistPhotonP =
      new TH1F("fHistPhotonP", "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistPhotonInvMass =
      new TH1F("fHistPhotonInvMass",
               "; Invariant mass [GeV/#it{c}^{2}]; Entries", 500, 0, 1);
  fHistPhotonInvMassPt =
      new TH2F("fHistPhotonInvMassPt",
               "; #it{p}_{T} [GeV/#it{c}]; Invariant mass [GeV/#it{c}^{2}]",
               500, 0, 10, 500, 0, 1);
  fHistPhotonInvMassEta =
      new TH2F("fHistPhotonInvMassEta",
               "; #eta; Invariant mass [GeV/#it{c}^{2}]", 500, 0, 2, 500, 0, 1);
  fHistPhotonCPAPt =
      new TH2F("fHistPhotonCPAPt", "; #it{p}_{T} [GeV/#it{c}]; cos#alpha", 500,
               0, 10, 500, 0.5, 1);
  fHistPhotonR =
      new TH2F("fHistPhotonR", "; #it{p}_{T} [GeV/#it{c}]; Conversion radius",
               500, 0, 10, 500, 0, 200);
  fHistPhotonArm =
      new TH2F("fHistPhotonArm", "; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1,
               1, 500, 0, 1);
  fHistPhotonDCAz =
      new TH1F("fHistPhotonDCAz", "; DCA z [cm]; Entries", 500, -25, 25);
  fHistPhotonDCAr =
      new TH1F("fHistPhotonDCAr", "; DCA r [cm]; Entries", 500, 0, 25);
  fHistPhotonEtaPhi = new TH2F("fHistPhotonEtaPhi", "; #eta; #phi", 200, -1, 1,
                               200, 0, 2 * TMath::Pi());
  fHistPhotonConvPointX = new TH1F(
      "fHistPhotonConvPointX", "; conv. point x [cm]; Entries", 500, -200, 200);
  fHistPhotonConvPointY = new TH1F(
      "fHistPhotonConvPointY", "; conv. point x [cm]; Entries", 500, -200, 200);
  fHistPhotonConvPointZ = new TH1F(
      "fHistPhotonConvPointZ", "; conv. point x [cm]; Entries", 500, -200, 200);
  fHistPhotonEleP =
      new TH1F("fHistPhotonEleP", "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistPhotonElePt = new TH1F("fHistPhotonElePt",
                              "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);

  fHistPhotonEleNsigmaTPC = new TH2F("fHistPhotonEleNsigmaTPC",
                                     "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC",
                                     500, 0, 10, 500, -10, 10);
  fHistPhotonEleNsigmaTPCPion = new TH2F(
      "fHistPhotonEleNsigmaTPCPion", "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC",
      500, 0, 10, 500, -10, 10);
  fHistPhotonEleTPCsignal = new TH2F("fHistPhotonEleTPCsignal",
                                     "; p [GeV/#it{c}];  TPC d#it{E}/d#it{x}",
                                     500, 0, 5, 500, 0, 200);

  fHistTwoPhotonPt =
      new TH1F("fHistTwoPhotonPt",
               "; #it{p}_{T} #gamma#gamma [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistTwoPhotonP =
      new TH1F("fHistTwoPhotonP", "; #it{p} #gamma#gamma [GeV/#it{c}]; Entries",
               500, 0, 10);
  fHistTwoPhotonInvMass = new TH1F(
      "fHistTwoPhotonInvMass",
      "; Invariant mass #gamma#gamma [GeV/#it{c}^{2}]; Entries", 500, 0, 1);
  fHistTwoPhotonInvMassPt = new TH2F(
      "fHistTwoPhotonInvMassPt",
      "; #it{p}_{T} #gamma#gamma [GeV/#it{c}]; Invariant mass [GeV/#it{c}^{2}]",
      500, 0, 10, 500, 0, 1);
  fHistTwoPhotonInvMassEta =
      new TH2F("fHistTwoPhotonInvMassEta",
               "; #eta; Invariant mass #gamma#gamma [GeV/#it{c}^{2}]", 500, 0,
               2, 500, 0, 1);

  fHistograms->Add(fHistNPhoton);
  fHistograms->Add(fHistPhotonPt);
  fHistograms->Add(fHistPhotonP);
  fHistograms->Add(fHistPhotonInvMass);
  fHistograms->Add(fHistPhotonInvMassPt);
  fHistograms->Add(fHistPhotonInvMassEta);
  fHistograms->Add(fHistPhotonCPAPt);
  fHistograms->Add(fHistPhotonR);
  fHistograms->Add(fHistPhotonArm);
  fHistograms->Add(fHistPhotonDCAz);
  fHistograms->Add(fHistPhotonDCAr);
  fHistograms->Add(fHistPhotonEleP);
  fHistograms->Add(fHistPhotonElePt);
  fHistograms->Add(fHistPhotonEtaPhi);
  fHistograms->Add(fHistPhotonConvPointX);
  fHistograms->Add(fHistPhotonConvPointY);
  fHistograms->Add(fHistPhotonConvPointZ);
  fHistograms->Add(fHistPhotonEleNsigmaTPC);
  fHistograms->Add(fHistPhotonEleNsigmaTPCPion);
  fHistograms->Add(fHistPhotonEleTPCsignal);
  fHistograms->Add(fHistTwoPhotonPt);
  fHistograms->Add(fHistTwoPhotonP);
  fHistograms->Add(fHistTwoPhotonInvMass);
  fHistograms->Add(fHistTwoPhotonInvMassPt);
  fHistograms->Add(fHistTwoPhotonInvMassEta);

  if (fIsMC) {
    if (fHistogramsMC != nullptr) {
      delete fHistogramsMC;
      fHistogramsMC = nullptr;
    }
    if (fHistogramsMC == nullptr) {
      fHistogramsMC = new TList();
      fHistogramsMC->SetOwner(kTRUE);
      fHistogramsMC->SetName("PhotonCut_MC");
    }

    fHistMCRecSigma0PhotonPt =
        new TH1F("fHistMCRecSigma0PhotonPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCRecSigma0PhotonP =
        new TH1F("fHistMCRecSigma0PhotonP", "; #it{p} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCRecSigma0PhotonInvMass =
        new TH1F("fHistMCRecSigma0PhotonInvMass",
                 "; Invariant mass [GeV/#it{c}^{2}]; Entries", 500, 0, 1);
    fHistMCRecSigma0PhotonInvMassPt =
        new TH2F("fHistMCRecSigma0PhotonInvMassPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Invariant mass [GeV/#it{c}^{2}]",
                 500, 0, 10, 500, 0, 1);
    fHistMCRecSigma0PhotonInvMassEta = new TH2F(
        "fHistMCRecSigma0PhotonInvMassEta",
        "; #eta; Invariant mass [GeV/#it{c}^{2}]", 500, 0, 2, 500, 0, 1);
    fHistMCRecSigma0PhotonR =
        new TH2F("fHistMCRecSigma0PhotonR",
                 "; #it{p}_{T} [GeV/#it{c}]; Conversion radius", 500, 0, 10,
                 500, 0, 200);
    fHistMCRecSigma0PhotonArm =
        new TH2F("fHistMCRecSigma0PhotonArm",
                 "; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistMCRecSigma0PhotonDCAz = new TH1F(
        "fHistMCRecSigma0PhotonDCAz", "; DCA z [cm]; Entries", 500, -25, 25);
    fHistMCRecSigma0PhotonDCAr = new TH1F("fHistMCRecSigma0PhotonDCAr",
                                          "; DCA r [cm]; Entries", 500, 0, 25);
    fHistMCRecSigma0PhotonConvPointX =
        new TH1F("fHistMCRecSigma0PhotonConvPointX",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCRecSigma0PhotonConvPointY =
        new TH1F("fHistMCRecSigma0PhotonConvPointY",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCRecSigma0PhotonConvPointZ =
        new TH1F("fHistMCRecSigma0PhotonConvPointZ",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCRecSigma0PhotonEleP =
        new TH1F("fHistMCRecSigma0PhotonEleP", "; #it{p} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCRecSigma0PhotonElePt =
        new TH1F("fHistMCRecSigma0PhotonElePt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCRecSigma0CPAPt = new TH2F("fHistMCRecSigma0CPAPt",
                                     "; #it{p}_{T} [GeV/#it{c}]; cos#alpha",
                                     500, 0, 10, 500, 0.5, 1);
    fHistMCRecSigma0PhotonDCAzPt = new TH2F(
        "fHistMCRecSigma0PhotonDCAzPt", "; #it{p}_{T} [GeV/#it{c}]; DCA z [cm]",
        500, 0, 10, 500, -25, 25);
    fHistMCRecSigma0PhotonDCArPt = new TH2F(
        "fHistMCRecSigma0PhotonDCArPt", "; #it{p}_{T} [GeV/#it{c}]; DCA r [cm]",
        500, 0, 10, 500, 0, 25);
    fHistMCRecSigma0PhotonPsiPairPt = new TH2F(
        "fHistMCRecSigma0PhotonPsiPairPt",
        "; #it{p}_{T} [GeV/#it{c}]; #Psi_{pair}", 500, 0, 10, 500, 0, 5);
    fHistMCRecSigma0PhotonChi2Pt = new TH2F(
        "fHistMCRecSigma0PhotonChi2Pt", "; #it{p}_{T} [GeV/#it{c}]; #chi^{2}",
        500, 0, 10, 500, 0, 100);

    fHistogramsMC->Add(fHistMCRecSigma0PhotonPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonP);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonInvMass);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonInvMassPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonInvMassEta);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonR);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonArm);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonDCAz);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonDCAr);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonConvPointX);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonConvPointY);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonConvPointZ);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonEleP);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonElePt);
    fHistogramsMC->Add(fHistMCRecSigma0CPAPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonDCAzPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonDCArPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonPsiPairPt);
    fHistogramsMC->Add(fHistMCRecSigma0PhotonChi2Pt);

    fHistMCFakeSigma0PhotonPt =
        new TH1F("fHistMCFakeSigma0PhotonPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCFakeSigma0PhotonP =
        new TH1F("fHistMCFakeSigma0PhotonP", "; #it{p} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCFakeSigma0PhotonInvMass =
        new TH1F("fHistMCFakeSigma0PhotonInvMass",
                 "; Invariant mass [GeV/#it{c}^{2}]; Entries", 5000, 0, 1);
    fHistMCFakeSigma0PhotonInvMassPt =
        new TH2F("fHistMCFakeSigma0PhotonInvMassPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Invariant mass [GeV/#it{c}^{2}]",
                 500, 0, 5, 500, 0, 1);
    fHistMCFakeSigma0PhotonInvMassEta = new TH2F(
        "fHistMCFakeSigma0PhotonInvMassEta",
        "; #eta; Invariant mass [GeV/#it{c}^{2}]", 500, 0, 2, 500, 0, 1);
    fHistMCFakeSigma0PhotonR = new TH2F(
        "fHistMCFakeSigma0PhotonR",
        "; #it{p}_{T} [GeV/#it{c}]; Conversion radius", 500, 0, 5, 500, 0, 200);
    fHistMCFakeSigma0PhotonArm =
        new TH2F("fHistMCFakeSigma0PhotonArm",
                 "; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistMCFakeSigma0PhotonDCAz = new TH1F(
        "fHistMCFakeSigma0PhotonDCAz", "; DCA z [cm]; Entries", 500, -25, 25);
    fHistMCFakeSigma0PhotonDCAr = new TH1F("fHistMCFakeSigma0PhotonDCAr",
                                           "; DCA r [cm]; Entries", 500, 0, 25);
    fHistMCFakeSigma0PhotonConvPointX =
        new TH1F("fHistMCFakeSigma0PhotonConvPointX",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCFakeSigma0PhotonConvPointY =
        new TH1F("fHistMCFakeSigma0PhotonConvPointY",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCFakeSigma0PhotonConvPointZ =
        new TH1F("fHistMCFakeSigma0PhotonConvPointZ",
                 "; conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCFakeSigma0PhotonEleP =
        new TH1F("fHistMCFakeSigma0PhotonEleP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCFakeSigma0PhotonElePt =
        new TH1F("fHistMCFakeSigma0PhotonElePt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCFakeSigma0CPAPt = new TH2F("fHistMCFakeSigma0CPAPt",
                                      "; #it{p}_{T} [GeV/#it{c}]; cos#alpha",
                                      500, 0, 10, 500, 0.5, 1);
    fHistMCFakeSigma0PhotonDCAzPt = new TH2F(
        "fHistMCFakeSigma0PhotonDCAzPt",
        "; #it{p}_{T} [GeV/#it{c}]; DCA z [cm]", 500, 0, 10, 500, -25, 25);
    fHistMCFakeSigma0PhotonDCArPt = new TH2F(
        "fHistMCFakeSigma0PhotonDCArPt",
        "; #it{p}_{T} [GeV/#it{c}]; DCA r [cm]", 500, 0, 10, 500, 0, 25);
    fHistMCFakeSigma0PhotonPsiPairPt = new TH2F(
        "fHistMCFakeSigma0PhotonPsiPairPt",
        "; #it{p}_{T} [GeV/#it{c}]; #Psi_{pair}", 500, 0, 10, 500, 0, 5);
    fHistMCFakeSigma0PhotonChi2Pt = new TH2F(
        "fHistMCFakeSigma0PhotonChi2Pt", "; #it{p}_{T} [GeV/#it{c}]; #chi^{2}",
        500, 0, 10, 500, 0, 100);

    fHistogramsMC->Add(fHistMCFakeSigma0PhotonPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonP);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonInvMass);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonInvMassPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonInvMassEta);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonR);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonArm);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonDCAz);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonDCAr);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonConvPointX);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonConvPointY);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonConvPointZ);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonEleP);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonElePt);
    fHistogramsMC->Add(fHistMCFakeSigma0CPAPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonDCAzPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonDCArPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonPsiPairPt);
    fHistogramsMC->Add(fHistMCFakeSigma0PhotonChi2Pt);

    fHistMCFakePhotonPt =
        new TH1F("fHistMCFakePhotonPt", "; #it{p}_{T} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCFakePhotonP = new TH1F("fHistMCFakePhotonP",
                                  "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCFakePhotonInvMass =
        new TH1F("fHistMCFakePhotonInvMass",
                 "; Invariant mass [GeV/#it{c}^{2}]; Entries", 5000, 0, 1);
    fHistMCFakePhotonInvMassPt =
        new TH2F("fHistMCFakePhotonInvMassPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Invariant mass [GeV/#it{c}^{2}]",
                 500, 0, 5, 500, 0, 1);
    fHistMCFakePhotonInvMassEta = new TH2F(
        "fHistMCFakePhotonInvMassEta",
        "; #eta; Invariant mass [GeV/#it{c}^{2}]", 500, 0, 2, 500, 0, 1);
    fHistMCFakePhotonR = new TH2F(
        "fHistMCFakePhotonR", "; #it{p}_{T} [GeV/#it{c}]; Conversion radius",
        500, 0, 5, 500, 0, 200);
    fHistMCFakePhotonArm =
        new TH2F("fHistMCFakePhotonArm", "; #alpha; #it{q}_{T} [GeV/#it{c}]",
                 500, -1, 1, 500, 0, 1);
    fHistMCFakePhotonDCAz = new TH1F("fHistMCFakePhotonDCAz",
                                     "; DCA z [cm]; Entries", 500, -25, 25);
    fHistMCFakePhotonDCAr =
        new TH1F("fHistMCFakePhotonDCAr", "; DCA r [cm]; Entries", 500, 0, 25);
    fHistMCFakePhotonConvPointX =
        new TH1F("fHistMCFakePhotonConvPointX", "; conv. point x [cm]; Entries",
                 500, -200, 200);
    fHistMCFakePhotonConvPointY =
        new TH1F("fHistMCFakePhotonConvPointY", "; conv. point x [cm]; Entries",
                 500, -200, 200);
    fHistMCFakePhotonConvPointZ =
        new TH1F("fHistMCFakePhotonConvPointZ", "; conv. point x [cm]; Entries",
                 500, -200, 200);
    fHistMCFakePhotonEleP = new TH1F(
        "fHistMCFakePhotonEleP", "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCFakePhotonElePt =
        new TH1F("fHistMCFakePhotonElePt", "; #it{p}_{T} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCFakePhotonCPAPt = new TH2F("fHistMCFakePhotonCPAPt",
                                      "; #it{p}_{T} [GeV/#it{c}]; cos#alpha",
                                      500, 0, 10, 500, 0.5, 1);
    fHistMCFakePhotonDCAzPt = new TH2F("fHistMCFakePhotonDCAzPt",
                                       "; #it{p}_{T} [GeV/#it{c}]; DCA z [cm]",
                                       500, 0, 10, 500, -25, 25);
    fHistMCFakePhotonDCArPt = new TH2F("fHistMCFakePhotonDCArPt",
                                       "; #it{p}_{T} [GeV/#it{c}]; DCA r [cm]",
                                       500, 0, 10, 500, 0, 25);
    fHistMCFakePhotonPsiPairPt = new TH2F(
        "fHistMCFakePhotonPsiPairPt", "; #it{p}_{T} [GeV/#it{c}]; #Psi_{pair}",
        500, 0, 10, 500, 0, 5);
    fHistMCFakePhotonChi2Pt = new TH2F("fHistMCFakePhotonChi2Pt",
                                       "; #it{p}_{T} [GeV/#it{c}]; #chi^{2}",
                                       500, 0, 10, 500, 0, 100);
    fHistogramsMC->Add(fHistMCFakePhotonPt);
    fHistogramsMC->Add(fHistMCFakePhotonP);
    fHistogramsMC->Add(fHistMCFakePhotonInvMass);
    fHistogramsMC->Add(fHistMCFakePhotonInvMassPt);
    fHistogramsMC->Add(fHistMCFakePhotonInvMassEta);
    fHistogramsMC->Add(fHistMCFakePhotonR);
    fHistogramsMC->Add(fHistMCFakePhotonArm);
    fHistogramsMC->Add(fHistMCFakePhotonDCAz);
    fHistogramsMC->Add(fHistMCFakePhotonDCAr);
    fHistogramsMC->Add(fHistMCFakePhotonConvPointX);
    fHistogramsMC->Add(fHistMCFakePhotonConvPointY);
    fHistogramsMC->Add(fHistMCFakePhotonConvPointZ);
    fHistogramsMC->Add(fHistMCFakePhotonEleP);
    fHistogramsMC->Add(fHistMCFakePhotonElePt);
    fHistogramsMC->Add(fHistMCFakePhotonCPAPt);
    fHistogramsMC->Add(fHistMCFakePhotonDCAzPt);
    fHistogramsMC->Add(fHistMCFakePhotonDCArPt);
    fHistogramsMC->Add(fHistMCFakePhotonPsiPairPt);
    fHistogramsMC->Add(fHistMCFakePhotonChi2Pt);

    fHistograms->Add(fHistogramsMC);
  }
}
