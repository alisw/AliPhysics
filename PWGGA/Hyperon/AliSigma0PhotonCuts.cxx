#include "AliSigma0PhotonCuts.h"
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODv0.h"
#include "TMath.h"
#include "AliV0ReaderV1.h"
#include "AliNanoAODTrack.h"

ClassImp(AliSigma0PhotonCuts)

//____________________________________________________________________________________________________
AliSigma0PhotonCuts::AliSigma0PhotonCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsBefore(nullptr),
      fHistogramsAfter(nullptr),
      fHistogramsPos(nullptr),
      fHistogramsNeg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fIsLightweight(false),
      fV0ReaderName("NoInit"),
      fIsMC(false),
      f2DArmenterosCut(false),
      fQtMax(0.),
      fPhotonPtMin(0.),
      fPhotonEtaMax(999),
      f2DPsiPairCut(false),
      fPsiPairMax(999),
      fChi2MaxFor2DPsiPair(999),
      fCPAMin(0.),
      fDCAzMax(999.f),
      fDCArMax(999.f),
      fElectronPtMin(0.),
      fElectronEtaMax(999),
      fElectronRatioFindable(0.),
      fElectronNSigmaTPCMax(999),
      fElectronNSigmaTPCMin(-999),
      fDoTransvRadRejection(false),
      fTransvRadRejectionLow(0.),
      fTransvRadRejectionUp(0.),
      fDoPhotonQualityCut(false),
      fPhotonQuality(-1),
      fDoPhotonPileupCut(false),
      fPhotonPileupCut(false),
      fDoExtraPiZeroRejection(false),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistPiZeroMass(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0MassPt(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
      fHistMCV0Pt(nullptr),
      fHistV0Mother(nullptr),
      fHistPsiPairBefore(nullptr),
      fHistPsiPairAfter(nullptr),
      fHistDecayVertexXBefore(nullptr),
      fHistDecayVertexYBefore(nullptr),
      fHistDecayVertexZBefore(nullptr),
      fHistDecayVertexXAfter(nullptr),
      fHistDecayVertexYAfter(nullptr),
      fHistDecayVertexZAfter(nullptr),
      fHistTransverseRadiusBefore(nullptr),
      fHistTransverseRadiusAfter(nullptr),
      fHistCosPABefore(nullptr),
      fHistCosPAAfter(nullptr),
      fHistDCArBefore(nullptr),
      fHistDCArAfter(nullptr),
      fHistDCAzBefore(nullptr),
      fHistDCAzAfterOthersBefore(nullptr),
      fHistDCAzAfterOthersBeforeQuality(),
      fHistDCAzAfter(nullptr),
      fHistDCA(nullptr),
      fHistDecayLength(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistQualityBefore(nullptr),
      fHistQualityAfter(nullptr),
      fHistTomography(nullptr),
      fHistMCTruthPhotonPt(nullptr),
      fHistMCTruthPhotonSigmaPt(nullptr),
      fHistMCPhotonPt(nullptr),
      fHistMCPhotonSigmaPt(nullptr),
      fHistMCTruthPhotonP(nullptr),
      fHistMCTruthPhotonSigmaP(nullptr),
      fHistMCPhotonP(nullptr),
      fHistMCPhotonSigmaP(nullptr),
      fHistMCPhotonSource(nullptr),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleChi2Before(),
      fHistSingleParticleChi2After(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleDCAtoPVBefore(),
      fHistSingleParticleDCAtoPVAfter(),
      fHistSingleParticlePileUp(),
      fHistSingleParticlePID(),
      fPIDResponse(nullptr) {
}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts::AliSigma0PhotonCuts(const AliSigma0PhotonCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsBefore(nullptr),
      fHistogramsAfter(nullptr),
      fHistogramsPos(nullptr),
      fHistogramsNeg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fIsLightweight(false),
      fV0ReaderName("NoInit"),
      fIsMC(false),
      f2DArmenterosCut(false),
      fQtMax(0.),
      fPhotonPtMin(0.),
      fPhotonEtaMax(999),
      f2DPsiPairCut(false),
      fPsiPairMax(999),
      fChi2MaxFor2DPsiPair(999),
      fCPAMin(0.),
      fDCAzMax(999.f),
      fDCArMax(999.f),
      fElectronPtMin(0.),
      fElectronEtaMax(999),
      fElectronRatioFindable(0.),
      fElectronNSigmaTPCMax(999),
      fElectronNSigmaTPCMin(-999),
      fDoTransvRadRejection(false),
      fTransvRadRejectionLow(0.),
      fTransvRadRejectionUp(0.),
      fDoPhotonQualityCut(false),
      fPhotonQuality(-1),
      fDoPhotonPileupCut(false),
      fPhotonPileupCut(false),
      fDoExtraPiZeroRejection(false),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistPiZeroMass(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0MassPt(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
      fHistMCV0Pt(nullptr),
      fHistV0Mother(nullptr),
      fHistPsiPairBefore(nullptr),
      fHistPsiPairAfter(nullptr),
      fHistDecayVertexXBefore(nullptr),
      fHistDecayVertexYBefore(nullptr),
      fHistDecayVertexZBefore(nullptr),
      fHistDecayVertexXAfter(nullptr),
      fHistDecayVertexYAfter(nullptr),
      fHistDecayVertexZAfter(nullptr),
      fHistTransverseRadiusBefore(nullptr),
      fHistTransverseRadiusAfter(nullptr),
      fHistCosPABefore(nullptr),
      fHistCosPAAfter(nullptr),
      fHistDCArBefore(nullptr),
      fHistDCArAfter(nullptr),
      fHistDCAzBefore(nullptr),
      fHistDCAzAfterOthersBefore(nullptr),
      fHistDCAzAfterOthersBeforeQuality(),
      fHistDCAzAfter(nullptr),
      fHistDCA(nullptr),
      fHistDecayLength(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistQualityBefore(nullptr),
      fHistQualityAfter(nullptr),
      fHistTomography(nullptr),
      fHistMCTruthPhotonPt(nullptr),
      fHistMCTruthPhotonSigmaPt(nullptr),
      fHistMCPhotonPt(nullptr),
      fHistMCPhotonSigmaPt(nullptr),
      fHistMCTruthPhotonP(nullptr),
      fHistMCTruthPhotonSigmaP(nullptr),
      fHistMCPhotonP(nullptr),
      fHistMCPhotonSigmaP(nullptr),
      fHistMCPhotonSource(nullptr),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleChi2Before(),
      fHistSingleParticleChi2After(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleDCAtoPVBefore(),
      fHistSingleParticleDCAtoPVAfter(),
      fHistSingleParticlePileUp(),
      fHistSingleParticlePID(),
      fPIDResponse(nullptr) {
}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts &AliSigma0PhotonCuts::operator=(
    const AliSigma0PhotonCuts &ref) {
  // Assignment operator
  if (this == &ref)
    return *this;
  return (*this);
}

AliSigma0PhotonCuts::~AliSigma0PhotonCuts() {
  delete fPIDResponse;
}

//____________________________________________________________________________________________________
AliSigma0PhotonCuts *AliSigma0PhotonCuts::PhotonCuts() {
  AliSigma0PhotonCuts *v0Cuts = new AliSigma0PhotonCuts();
  v0Cuts->SetArmenteros2DCut(true);
  v0Cuts->SetArmenterosQtMax(0.06);
  v0Cuts->SetPhotonPtMin(0.02);
  v0Cuts->SetPhotonEtaMax(999);
  v0Cuts->SetPsiPair2DCut(true);
  v0Cuts->SetPsiPairMax(0.2);
  v0Cuts->SetChi2ForPsiMax(30);
  v0Cuts->SetCPAMin(0.999);
  v0Cuts->SetElectronPtMin(0.05);
  v0Cuts->SetElectronEtaMax(0.9);
  v0Cuts->SetElectronRatioFindable(0.35);
  v0Cuts->SetElectronNSigmaTPCMax(7);
  v0Cuts->SetElectronNSigmaTPCMin(-6);
  v0Cuts->SetDCAzMax(0.5);
  v0Cuts->SetDCArMax(0.75);
  return v0Cuts;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonCuts::PhotonCuts(
    AliAODEvent *inputEvent, AliMCEvent *mcEvent, const TClonesArray *photons,
    std::vector<AliFemtoDreamBasePart> &container) {
  container.clear();

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = static_cast<AliInputEventHandler *>(man
      ->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  fHistNV0->Fill(photons->GetEntriesFast());

  auto event = static_cast<AliVEvent*>(inputEvent);

  // Look at the MC Truth
  if (fIsMC) {
    for (int iPart = 1; iPart < (mcEvent->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) mcEvent->GetTrack(iPart);
      if (!mcPart->IsPhysicalPrimary())
        continue;
      int pdg = mcPart->GetPdgCode();
      if (pdg == 22) {
        fHistMCTruthPhotonPt->Fill(mcPart->Pt());
        fHistMCTruthPhotonP->Fill(mcPart->P());
        const int pdgMother = ((AliAODMCParticle*) (AliAODMCParticle*) mcEvent
            ->GetTrack(mcPart->GetMother()))->GetPdgCode();
        fHistMCPhotonSource->Fill(mcPart->P(), TMath::Abs(pdgMother));

        if (TMath::Abs(pdgMother) == 3212) {
          fHistMCTruthPhotonSigmaPt->Fill(mcPart->Pt());
          fHistMCTruthPhotonSigmaP->Fill(mcPart->P());
        }
      }
    }
  }

  for (int iGamma = 0; iGamma < photons->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate = dynamic_cast<AliAODConversionPhoton *>(photons->At(
        iGamma));
    fHistCuts->Fill(0.f);

    if (!PhotonCandidate)
      continue;

    AliAODv0 *v0 = inputEvent->GetV0(PhotonCandidate->GetV0Index());
    fHistCuts->Fill(1.f);

    auto pos = GetTrack(dynamic_cast<AliVEvent*>(inputEvent),
                        PhotonCandidate->GetTrackLabelPositive());
    auto neg = GetTrack(dynamic_cast<AliVEvent*>(inputEvent),
                        PhotonCandidate->GetTrackLabelNegative());
    if (!pos || !neg)
      continue;
    fHistCuts->Fill(2.f);

    auto nanoPos = dynamic_cast<AliNanoAODTrack *>(pos);
    if(nanoPos) v0 = nullptr; // for NanoAODs we dont have a v0 matching!

    if (fDoExtraPiZeroRejection) {
      if (!PiZeroRejection(PhotonCandidate, photons, iGamma)) {
        continue;
      }
    }

    if (ProcessPhoton(event, mcEvent, PhotonCandidate, v0, pos, neg)) {
      AliFemtoDreamBasePart photon(PhotonCandidate, pos, neg, inputEvent);
      if (fIsMC) {
        TClonesArray *mcarray = dynamic_cast<TClonesArray *>(event->FindListObject(
            "mcparticles"));
        if (!mcarray) {
          AliError("PhotonCuts: MC Array not found");
        }

        RelabelAODPhotonCandidates(PhotonCandidate, event);
        const int labelPos = PhotonCandidate->GetMCLabelPositive();
        const int labelNeg = PhotonCandidate->GetMCLabelNegative();
        AliAODMCParticle * mcPartPos = (AliAODMCParticle*) mcarray->At(labelPos);
        AliAODMCParticle * mcPartNeg = (AliAODMCParticle*) mcarray->At(labelNeg);

        if (mcPartPos && mcPartNeg) {
          if (mcPartPos->GetMother() > -1
              && (mcPartNeg->GetMother() == mcPartPos->GetMother())) {
            AliAODMCParticle *mcParticle = (AliAODMCParticle*) mcarray->At(
                mcPartPos->GetMother());
            photon.SetID(mcPartPos->GetMother());

            if (mcParticle) {
              photon.SetMCParticle(mcParticle, mcEvent);
              photon.SetMCPDGCode(mcParticle->GetPdgCode());
              photon.SetMotherID(mcParticle->GetMother());  // otherwise the Sigma0 is not set properly as the mother of the photon

              if (mcParticle->PdgCode() == 22) {
                fHistMCV0Pt->Fill(photon.GetPt());
                const int pdgMother = ((AliAODMCParticle*) mcarray->At(
                    mcParticle->GetMother()))->GetPdgCode();
                fHistV0Mother->Fill(photon.GetPt(), TMath::Abs(pdgMother));

                fHistMCPhotonPt->Fill(photon.GetPt());
                fHistMCPhotonP->Fill(photon.GetP());
                if (TMath::Abs(pdgMother) == 3212) {
                  fHistMCPhotonSigmaPt->Fill(photon.GetPt());
                  fHistMCPhotonSigmaP->Fill(photon.GetP());
                }
              }
            }
          }
        }
      }
      container.push_back(photon);
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonCuts::ProcessPhoton(AliVEvent* event, AliMCEvent *mcEvent,
                                        AliAODConversionPhoton *PhotonCandidate,
                                        AliAODv0 *v0, AliVTrack *pos,
                                        AliVTrack *neg) {
  AliNanoAODTrack *nanoPos = dynamic_cast<AliNanoAODTrack *>(pos);
  AliNanoAODTrack *nanoNeg = dynamic_cast<AliNanoAODTrack *>(neg);

  // compute the relevant variables
  const float pt = PhotonCandidate->GetPhotonPt();
  const float invMass = PhotonCandidate->GetPhotonMass();
  double momV0[3] = { 0, 0, 0 };
  momV0[0] = PhotonCandidate->Px();
  momV0[1] = PhotonCandidate->Py();
  momV0[2] = PhotonCandidate->Pz();
  double conv[3] { 0, 0, 0 };
  conv[0] = PhotonCandidate->GetConversionX();
  conv[1] = PhotonCandidate->GetConversionY();
  conv[2] = PhotonCandidate->GetConversionZ();
  const AliVVertex *vertex = event->GetPrimaryVertex();
  const float xPV = vertex->GetX();
  const float yPV = vertex->GetY();
  const float zPV = vertex->GetZ();

  PhotonCandidate->CalculateDistanceOfClossetApproachToPrimVtx(vertex);
  const float DCAz = PhotonCandidate->GetDCAzToPrimVtx();
  const float DCAr = PhotonCandidate->GetDCArToPrimVtx();
  const float DCA = std::sqrt(DCAz * DCAz + DCAr * DCAr);

  const double PosV0[3] = { conv[0] - xPV, conv[1] - yPV, conv[2] - zPV };
  // Recalculated V0 Position vector

  const double momV02 = momV0[0] * momV0[0] + momV0[1] * momV0[1]
      + momV0[2] * momV0[2];
  const double PosV02 = PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1]
      + PosV0[2] * PosV0[2];

  const float cosinePointingAngle =
      (momV02 * PosV02 > 0.0) ?
          (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] + PosV0[2] * momV0[2])
              / std::sqrt(momV02 * PosV02) :
          -999.f;

  const float transRadius = PhotonCandidate->GetConversionRadius();
  const float decayLength = std::sqrt(
      conv[0] * conv[0] + conv[1] * conv[1] + conv[2] * conv[2]);

  const float massLambda = ComputeInvMass(v0, pos, neg, 2212, 211);
  const float massK0 = ComputeInvMass(v0, pos, neg, 211, 211);
  const float massAntiLambda = ComputeInvMass(v0, pos, neg, 211, 2212);

  const float armAlpha = PhotonCandidate->GetArmenterosAlpha();
  const float armQt = PhotonCandidate->GetArmenterosQt();

  const float eta = PhotonCandidate->GetPhotonEta();
  const float phi = PhotonCandidate->GetPhotonPhi();
  const float psiPair = PhotonCandidate->GetPsiPair();

  const float posPt = pos->Pt();
  const float negPt = neg->Pt();
  const float posEta = pos->Eta();
  const float negEta = neg->Eta();
  const float dcaDaughterToPVPos = (v0) ? v0->DcaPosToPrimVertex() : 0;
  const short nClsTPCPos = pos->GetTPCNcls();
  const float nCrossedRowsPos =
      (nanoPos) ? nanoPos->GetTPCNCrossedRows() : pos->GetTPCClusterInfo(2, 1);
  const short nFindablePos = pos->GetTPCNclsF();
  const float ratioFindablePos = nCrossedRowsPos
      / (static_cast<float>(nFindablePos) + 1.e-3);
  float chi2Pos =
      (nClsTPCPos > 5) ? (pos->GetTPCchi2() / float(nClsTPCPos - 5)) : -1.;
  if(nanoPos) chi2Pos = nanoPos->Chi2perNDF();

  const float dcaDaughterToPVNeg = (v0) ? v0->DcaNegToPrimVertex() : 0;
  const short nClsTPCNeg = neg->GetTPCNcls();
  const float nCrossedRowsNeg =
      (nanoNeg) ? nanoNeg->GetTPCNCrossedRows() : neg->GetTPCClusterInfo(2, 1);
  const short nFindableNeg = neg->GetTPCNclsF();
  const float ratioFindableNeg = nCrossedRowsNeg
      / (static_cast<float>(nFindableNeg) + 1.e-3);
  float chi2Neg =
      (nClsTPCNeg > 5) ? (neg->GetTPCchi2() / float(nClsTPCNeg - 5)) : -1.;
  if(nanoNeg) chi2Neg = nanoNeg->Chi2perNDF();

  float pidPos = -999.f;
  float pidNeg = -999.f;
  if (fPIDResponse && !(nanoPos && nanoNeg)) {
    pidPos = fPIDResponse->NumberOfSigmasTPC(pos, AliPID::kElectron);
    pidNeg = fPIDResponse->NumberOfSigmasTPC(neg, AliPID::kElectron);
  } else if(nanoPos && nanoNeg) {
    pidPos = nanoPos->GetVar(AliNanoAODTrack::GetPIDIndex(
                        AliNanoAODTrack::kSigmaTPC, AliPID::kElectron));
    pidNeg = nanoNeg->GetVar(AliNanoAODTrack::GetPIDIndex(
                        AliNanoAODTrack::kSigmaTPC, AliPID::kElectron));
  }

  int nClusterITSneg = 0;
  int nClusterITSpos = 0;
  for (int i = 0; i < 6; ++i) {
    nClusterITSneg += neg->HasPointOnITSLayer(i);
    nClusterITSpos += pos->HasPointOnITSLayer(i);
  }

  int photonQuality = 0;
  if (nClusterITSneg > 1 && nClusterITSpos > 1) {
    photonQuality = 3;
  } else if (nClusterITSneg > 1 || nClusterITSpos > 1) {
    photonQuality = 2;
  } else {
    photonQuality = 1;
  }
  bool pileUpPhoton = (photonQuality == 1) ? true : false;

  // BEFORE THE CUTS
  if (!fIsLightweight) {
    fHistDecayVertexXBefore->Fill(pt, TMath::Abs(conv[0]));
    fHistDecayVertexYBefore->Fill(pt, TMath::Abs(conv[1]));
    fHistDecayVertexZBefore->Fill(pt, TMath::Abs(conv[2]));
    fHistTransverseRadiusBefore->Fill(pt, transRadius);
    fHistCosPABefore->Fill(pt, cosinePointingAngle);
    fHistArmenterosBefore->Fill(armAlpha, armQt);
    fHistDCArBefore->Fill(pt, DCAr);
    fHistDCAzBefore->Fill(pt, DCAz);
    fHistPsiPairBefore->Fill(pt, psiPair);
    fHistQualityBefore->Fill(pt, photonQuality);

    fHistSingleParticleEtaBefore[0]->Fill(posPt, pos->Eta());
    fHistSingleParticleChi2Before[0]->Fill(posPt, chi2Pos);
    fHistSingleParticleNclsTPCBefore[0]->Fill(posPt, nClsTPCPos);
    fHistSingleParticleNclsTPCFindableBefore[0]->Fill(posPt, nFindablePos);
    fHistSingleParticleNclsTPCRatioFindableBefore[0]->Fill(posPt,
                                                           ratioFindablePos);
    fHistSingleParticleNcrossedTPCBefore[0]->Fill(posPt, nCrossedRowsPos);
    fHistSingleParticleDCAtoPVBefore[0]->Fill(posPt, dcaDaughterToPVPos);

    fHistSingleParticleEtaBefore[1]->Fill(negPt, neg->Eta());
    fHistSingleParticleChi2Before[1]->Fill(negPt, chi2Neg);
    fHistSingleParticleNclsTPCBefore[1]->Fill(negPt, nClsTPCNeg);
    fHistSingleParticleNclsTPCFindableBefore[1]->Fill(negPt, nFindableNeg);
    fHistSingleParticleNclsTPCRatioFindableBefore[1]->Fill(negPt,
                                                           ratioFindableNeg);
    fHistSingleParticleNcrossedTPCBefore[1]->Fill(negPt, nCrossedRowsNeg);
    fHistSingleParticleDCAtoPVBefore[1]->Fill(negPt, dcaDaughterToPVNeg);
  }

  // DO THE CUTS
  if (f2DArmenterosCut) {
    const float maxPhotonAsymmetry = 0.95;
    if (!(TMath::Power(armAlpha / maxPhotonAsymmetry, 2)
        + TMath::Power(armQt / fQtMax, 2) < 1)) {
      return false;
    }
  } else {
    if (armQt > fQtMax) {
      return false;
    }
  }
  fHistCuts->Fill(3.f);

  if (pt < fPhotonPtMin) {
    return false;
  }
  fHistCuts->Fill(4.f);

  if (TMath::Abs(eta) > fPhotonEtaMax) {
    return false;
  }
  fHistCuts->Fill(5.f);

  if (f2DPsiPairCut) {
    if (!(TMath::Abs(psiPair)
        < -fPsiPairMax / fChi2MaxFor2DPsiPair * PhotonCandidate->GetChi2perNDF()
            + fPsiPairMax)) {
      return false;
    }
  } else {
    if (TMath::Abs(psiPair) > fPsiPairMax) {
      return false;
    }
  }
  fHistCuts->Fill(6.f);

  if (cosinePointingAngle < fCPAMin) {
    return false;
  }
  fHistCuts->Fill(7.f);

  if(DCAr > fDCArMax) {
    return false;
  }
  fHistCuts->Fill(9.f);

  if (posPt < fElectronPtMin || negPt < fElectronPtMin) {
    return false;
  }
  fHistCuts->Fill(10.f);

  if (TMath::Abs(posEta) > fElectronEtaMax
      || TMath::Abs(negEta) > fElectronEtaMax) {
    return false;
  }
  fHistCuts->Fill(11.f);

  if (ratioFindableNeg < fElectronRatioFindable
      || ratioFindablePos < fElectronRatioFindable) {
    return false;
  }
  fHistCuts->Fill(12.f);

  if (pidPos > fElectronNSigmaTPCMax || pidNeg > fElectronNSigmaTPCMax
      || pidPos < fElectronNSigmaTPCMin || pidNeg < fElectronNSigmaTPCMin) {
    return false;
  }
  fHistCuts->Fill(13.f);

  if (fDoTransvRadRejection && transRadius > fTransvRadRejectionLow
      && transRadius < fTransvRadRejectionUp) {
    return false;
  }
  fHistCuts->Fill(14.f);

  if (fDoPhotonQualityCut && photonQuality != fPhotonQuality) {
    return false;
  }
  if (fDoPhotonPileupCut && pileUpPhoton != fPhotonPileupCut) {
    return false;
  }
  fHistCuts->Fill(15.f);

  fHistDCAzAfterOthersBefore->Fill(pt, DCAz);
  fHistDCAzAfterOthersBeforeQuality[photonQuality]->Fill(pt, DCAz);

  if (TMath::Abs(DCAz) > fDCAzMax) {
    return false;
  }
  fHistCuts->Fill(8.f);

  // AFTER THE CUTS
  fHistV0MassPt->Fill(pt, invMass);

  if (!fIsLightweight) {
    fHistK0Mass->Fill(massK0);
    fHistLambdaMass->Fill(massLambda);
    fHistAntiLambdaMass->Fill(massAntiLambda);
    fHistV0Pt->Fill(pt);
    fHistV0Mass->Fill(invMass);
    fHistCosPA->Fill(pt, cosinePointingAngle);
    fHistEtaPhi->Fill(eta, phi);
    fHistPsiPairAfter->Fill(pt, psiPair);
    fHistDecayVertexXAfter->Fill(pt, TMath::Abs(conv[0]));
    fHistDecayVertexYAfter->Fill(pt, TMath::Abs(conv[1]));
    fHistDecayVertexZAfter->Fill(pt, TMath::Abs(conv[2]));
    fHistTransverseRadiusAfter->Fill(pt, transRadius);
    fHistCosPAAfter->Fill(pt, cosinePointingAngle);
    fHistArmenterosAfter->Fill(armAlpha, armQt);
    fHistQualityAfter->Fill(pt, photonQuality);

    fHistDCArAfter->Fill(pt, DCAr);
    fHistDCAzAfter->Fill(pt, DCAz);
    fHistDCA->Fill(pt, DCA);
    fHistDecayLength->Fill(pt, decayLength);

    bool posTrackITS = (pos->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1)
        || pos->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
    bool negTrackITS = (neg->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1)
        || neg->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
    bool posTrackTOF = pos->GetTOFBunchCrossing() == 0;
    bool negTrackTOF = neg->GetTOFBunchCrossing() == 0;

    bool posTrackCombined = (posTrackITS || posTrackTOF);
    bool negTrackCombined = (negTrackITS || negTrackTOF);

    if (!fIsLightweight) {
      if (posTrackITS)
        fHistSingleParticlePileUp[0]->Fill(0.f, posPt);
      if (posTrackTOF)
        fHistSingleParticlePileUp[0]->Fill(1, posPt);
      if (posTrackCombined)
        fHistSingleParticlePileUp[0]->Fill(2, posPt);
      if (!posTrackITS && !posTrackTOF)
        fHistSingleParticlePileUp[0]->Fill(3, posPt);
      if (negTrackITS)
        fHistSingleParticlePileUp[1]->Fill(0.f, negPt);
      if (negTrackTOF)
        fHistSingleParticlePileUp[1]->Fill(1, negPt);
      if (negTrackCombined)
        fHistSingleParticlePileUp[1]->Fill(2, negPt);
      if (!negTrackITS && !negTrackTOF)
        fHistSingleParticlePileUp[1]->Fill(3, negPt);
    }

    fHistSingleParticlePt[0]->Fill(posPt);
    fHistSingleParticleEtaAfter[0]->Fill(posPt, pos->Eta());
    fHistSingleParticleChi2After[0]->Fill(posPt, chi2Pos);
    fHistSingleParticleNclsTPCAfter[0]->Fill(posPt, nClsTPCPos);
    fHistSingleParticleNclsTPCFindableAfter[0]->Fill(posPt, nFindablePos);
    fHistSingleParticleNclsTPCRatioFindableAfter[0]->Fill(posPt,
                                                          ratioFindablePos);
    fHistSingleParticleNcrossedTPCAfter[0]->Fill(posPt, nCrossedRowsPos);
    fHistSingleParticleDCAtoPVAfter[0]->Fill(posPt, dcaDaughterToPVPos);
    fHistSingleParticlePID[0]->Fill(posPt, pidPos);

    fHistSingleParticlePt[1]->Fill(negPt);
    fHistSingleParticleEtaAfter[1]->Fill(negPt, neg->Eta());
    fHistSingleParticleChi2After[1]->Fill(negPt, chi2Neg);
    fHistSingleParticleNclsTPCAfter[1]->Fill(negPt, nClsTPCNeg);
    fHistSingleParticleNclsTPCFindableAfter[1]->Fill(negPt, nFindableNeg);
    fHistSingleParticleNclsTPCRatioFindableAfter[1]->Fill(negPt,
                                                          ratioFindableNeg);
    fHistSingleParticleNcrossedTPCAfter[1]->Fill(negPt, nCrossedRowsNeg);
    fHistSingleParticleDCAtoPVAfter[1]->Fill(negPt, dcaDaughterToPVNeg);
    fHistSingleParticlePID[1]->Fill(negPt, pidNeg);
    fHistTomography->Fill(conv[0], conv[1]);
  }
  return true;
}


bool AliSigma0PhotonCuts::PiZeroRejection(AliAODConversionPhoton* photon,
                                          const TClonesArray *photons,
                                          int iPhoton) {
  TLorentzVector track1, track2;
  track1.SetXYZM(photon->Px(), photon->Py(), photon->Pz(),
                 photon->GetPhotonMass());

  for (int iGamma = 0; iGamma < photons->GetEntriesFast(); ++iGamma) {
    if (iGamma == iPhoton) {
      continue;
    }
    auto *PhotonCandidate = dynamic_cast<AliAODConversionPhoton *>(photons->At(
        iGamma));
    track2.SetXYZM(PhotonCandidate->Px(), PhotonCandidate->Py(),
                   PhotonCandidate->Pz(), PhotonCandidate->GetPhotonMass());
    TLorentzVector trackSum = track1 + track2;
    const double mass = trackSum.M();
    fHistPiZeroMass->Fill(mass);

    if (mass > 0.125 && mass < 0.15) {
      return false;
    }
  }
  return true;
}

///________________________________________________________________________
AliVTrack *AliSigma0PhotonCuts::GetTrack(AliVEvent * event, int label) {
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  if (label == -999999) {
    return nullptr;  // if AOD relabelling goes wrong, immediately return NULL
  }
  AliVTrack *track;
  if (AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data())
      && ((AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data()))->AreAODsRelabeled()) {
    if (event->GetTrack(label)) {
      track = dynamic_cast<AliVTrack*>(event->GetTrack(label));
      return track;
    }
  } else {
    for (Int_t ii = 0; ii < event->GetNumberOfTracks(); ii++) {
      track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
      if (track && track->GetID() == label) {
        return track;
      }
    }
  }
  return nullptr;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonCuts::RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate, AliVEvent *event){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  Bool_t AODLabelPos = kFALSE;
  Bool_t AODLabelNeg = kFALSE;

  for(Int_t i = 0; i<event->GetNumberOfTracks();i++){
    AliVTrack *tempDaughter = static_cast<AliVTrack*>(event->GetTrack(i));
    AliNanoAODTrack *nanoTempDaughter = dynamic_cast<AliNanoAODTrack *>(tempDaughter);
    if(!AODLabelPos){
      if( nanoTempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
      }
    }
    if(!AODLabelNeg){
      if( nanoTempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
      }
    }
    if(AODLabelNeg && AODLabelPos){
      return;
    }
  }
  if(!AODLabelPos || !AODLabelNeg){
    AliError(Form("NO AOD Daughters Found Pos: %i %i Neg: %i %i, setting all labels to -999999",AODLabelPos,PhotonCandidate->GetTrackLabelPositive(),AODLabelNeg,PhotonCandidate->GetTrackLabelNegative()));
    if(!AODLabelNeg){
      PhotonCandidate->SetMCLabelNegative(-999999);
      PhotonCandidate->SetLabelNegative(-999999);
    }
    if(!AODLabelPos){
      PhotonCandidate->SetMCLabelPositive(-999999);
      PhotonCandidate->SetLabelPositive(-999999);
    }
  }

}

//____________________________________________________________________________________________________
float AliSigma0PhotonCuts::ComputeInvMass(AliAODv0 *v0, AliVTrack *pos,
                                          AliVTrack *neg, int pdgPos,
                                          int pgdNeg) const {
  double posP[3], negP[3];
  pos->GetPxPyPz(posP);
  neg->GetPxPyPz(negP);
  float massDP = TDatabasePDG::Instance()->GetParticle(pdgPos)->Mass();
  float massDN = TDatabasePDG::Instance()->GetParticle(pgdNeg)->Mass();
  TLorentzVector trackPos, trackNeg;
  trackPos.SetXYZM(posP[0], posP[1], posP[2], massDP);
  trackNeg.SetXYZM(negP[0], negP[1], negP[2], massDN);
  TLorentzVector trackSum = trackPos + trackNeg;
  return trackSum.M();
}

//____________________________________________________________________________________________________
void AliSigma0PhotonCuts::InitCutHistograms(TString appendix) {
  const float pi = TMath::Pi();

  TH1::AddDirectory(kFALSE);
  TString name;
  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    name = "V0_" + appendix;
    fHistograms->SetName(name);
  }

  if (!fIsLightweight) {
    if (fHistogramsPos != nullptr) {
      delete fHistogramsPos;
      fHistogramsPos = nullptr;
    }
    if (fHistogramsPos == nullptr) {
      fHistogramsPos = new TList();
      fHistogramsPos->SetOwner(kTRUE);
      fHistogramsPos->SetName("V0_PosDaughter");
    }

    if (fHistogramsNeg != nullptr) {
      delete fHistogramsNeg;
      fHistogramsNeg = nullptr;
    }
    if (fHistogramsNeg == nullptr) {
      fHistogramsNeg = new TList();
      fHistogramsNeg->SetOwner(kTRUE);
      fHistogramsNeg->SetName("V0_NegDaughter");
    }
  }

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 17, 0, 17);
  fHistCutBooking->GetXaxis()->SetBinLabel(1, "2-D Armenteros");
  fHistCutBooking->GetXaxis()->SetBinLabel(2, "#it{q}_{T} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(3, "#it{p}_{T} min");
  fHistCutBooking->GetXaxis()->SetBinLabel(4, "#eta max");
  fHistCutBooking->GetXaxis()->SetBinLabel(5, "2-D #Psi_{pair} - #chi^{2}");
  fHistCutBooking->GetXaxis()->SetBinLabel(6, "#Psi_{pair} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(7, "#chi^{2} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(8, "Electron #it{p}_{T} min");
  fHistCutBooking->GetXaxis()->SetBinLabel(9, "Electron #eta max");
  fHistCutBooking->GetXaxis()->SetBinLabel(10, "Electron ratio findable");
  fHistCutBooking->GetXaxis()->SetBinLabel(11, "n#sigma TPC max");
  fHistCutBooking->GetXaxis()->SetBinLabel(12, "n#sigma TPC min");
  fHistCutBooking->GetXaxis()->SetBinLabel(13, "DCA_{z} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(14, "DCA_{r} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(15, "Transv. radius rej. low");
  fHistCutBooking->GetXaxis()->SetBinLabel(16, "Transv. radius rej. up");
  fHistCutBooking->GetXaxis()->SetBinLabel(17, "Photon quality");
  fHistograms->Add(fHistCutBooking);

  fHistCutBooking->Fill(0.f, static_cast<double>(f2DArmenterosCut));
  fHistCutBooking->Fill(1.f, fQtMax);
  fHistCutBooking->Fill(2.f, fPhotonPtMin);
  fHistCutBooking->Fill(3.f, fPhotonEtaMax);
  fHistCutBooking->Fill(4.f, static_cast<double>(f2DPsiPairCut));
  fHistCutBooking->Fill(5.f, fPsiPairMax);
  fHistCutBooking->Fill(6.f, fChi2MaxFor2DPsiPair);
  fHistCutBooking->Fill(7.f, fElectronPtMin);
  fHistCutBooking->Fill(8.f, fElectronEtaMax);
  fHistCutBooking->Fill(9.f, fElectronRatioFindable);
  fHistCutBooking->Fill(10.f, fElectronNSigmaTPCMax);
  fHistCutBooking->Fill(11.f, fElectronNSigmaTPCMin);
  fHistCutBooking->Fill(12.f, fDCAzMax);
  fHistCutBooking->Fill(13.f, fDCArMax);
  fHistCutBooking->Fill(14.f, fTransvRadRejectionLow);
  fHistCutBooking->Fill(15.f, fTransvRadRejectionUp);
  fHistCutBooking->Fill(16.f, fPhotonQuality);

  fHistV0MassPt = new TH2F(
      "InvMassPt", "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
      100, 0, 10, 100, 0., .1);
  fHistograms->Add(fHistV0MassPt);

  fHistNV0 = new TH1F("fHistNV0", ";Number of V0 candidates; Entries", 16, 0,
                      16);
  fHistograms->Add(fHistNV0);

  fHistCuts = new TH1F("fHistCuts", ";;Entries", 15, 0, 15);
  fHistCuts->GetXaxis()->SetBinLabel(1, "Photon candidates");
  fHistCuts->GetXaxis()->SetBinLabel(2, "V_{0} Matching");
  fHistCuts->GetXaxis()->SetBinLabel(3, "Track Matching");
  fHistCuts->GetXaxis()->SetBinLabel(4, "Armenteros");
  fHistCuts->GetXaxis()->SetBinLabel(5, "#it{p}_{T} min");
  fHistCuts->GetXaxis()->SetBinLabel(6, "#eta max");
  fHistCuts->GetXaxis()->SetBinLabel(7, "#Psi_{pair}");
  fHistCuts->GetXaxis()->SetBinLabel(8, "cos#alpha");
  fHistCuts->GetXaxis()->SetBinLabel(9, "DCA_{z}");
  fHistCuts->GetXaxis()->SetBinLabel(10, "DCA_{r}");
  fHistCuts->GetXaxis()->SetBinLabel(11, "Electron #it{p}_{T} min");
  fHistCuts->GetXaxis()->SetBinLabel(12, "Electron #eta max");
  fHistCuts->GetXaxis()->SetBinLabel(13, "Electron ratio findable");
  fHistCuts->GetXaxis()->SetBinLabel(14, "n#sigma TPC");
  fHistCuts->GetXaxis()->SetBinLabel(15, "Transverse radius rejection");
  fHistCuts->GetXaxis()->SetBinLabel(16, "Photon quality");
  fHistograms->Add(fHistCuts);

  if (!fIsLightweight) {
    fHistLambdaMass = new TH1F(
        "fHistLambdaMass",
        "; Invariant mass p#pi^{-} hypothesis (GeV/#it{c}^{2}); Entries", 125,
        1., 1.25);
    fHistograms->Add(fHistLambdaMass);

    fHistAntiLambdaMass = new TH1F(
        "fHistAntiLambdaMass",
        "; Invariant mass #bar{p}#pi hypothesis (GeV/#it{c}^{2}); Entries", 125,
        1., 1.25);
    fHistograms->Add(fHistAntiLambdaMass);

    fHistK0Mass = new TH1F(
        "fHistK0Mass",
        "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries", 125,
        0.4, 0.65);
    fHistograms->Add(fHistK0Mass);

    fHistV0Pt = new TH1F("fHistV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries", 100,
                         0, 10);
    fHistograms->Add(fHistV0Pt);

    fHistV0Mass = new TH1F("fHistV0Mass",
                           "; Invariant mass (GeV/#it{c}^{2}); Entries", 200,
                           0., 0.2);
    fHistograms->Add(fHistV0Mass);

    fHistCosPA = new TH2F("fHistCosPA",
                          "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 100, 0, 10,
                          100, 0.95, 1);
    fHistograms->Add(fHistCosPA);

    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 100, -1, 1, 100, 0,
                           2 * pi);
    fHistograms->Add(fHistEtaPhi);

    if(fDoExtraPiZeroRejection) {
      fHistPiZeroMass = new TH1F("fHistPiZeroMass", "; #it{M}_{#gamma#gamma} (GeV/#it{c}; Entries", 1000, 0, 0.25);
      fHistograms->Add(fHistPiZeroMass);
    }

    if (fIsMC) {
      fHistMCV0Pt = new TH1F("fHistMCV0Pt",
                             "; #it{p}_{T} (GeV/#it{c}); Entries", 100, 0, 10);
      fHistV0Mother = new TH2F("fHistV0Mother",
                               "; #it{p}_{T} (GeV/#it{c}); PDG code mother",
                               100, 0, 10, 4000, 0, 4000);
      fHistograms->Add(fHistMCV0Pt);
      fHistograms->Add(fHistV0Mother);

      fHistMCTruthPhotonPt = new TH1F("fHistMCTruthPhotonPt",
                                      "; #it{p}_{T} (GeV/#it{c}); Entries", 200,
                                      0, 10);
      fHistMCTruthPhotonSigmaPt = new TH1F("fHistMCTruthPhotonSigmaPt",
                                           "; #it{p}_{T} (GeV/#it{c}); Entries",
                                           200, 0, 10);
      fHistMCPhotonPt = new TH1F("fHistMCPhotonPt",
                                 "; #it{p}_{T} (GeV/#it{c}); Entries", 200, 0,
                                 10);
      fHistMCPhotonSigmaPt = new TH1F("fHistMCPhotonSigmaPt",
                                      "; #it{p}_{T} (GeV/#it{c}); Entries", 200,
                                      0, 10);
      fHistograms->Add(fHistMCTruthPhotonPt);
      fHistograms->Add(fHistMCTruthPhotonSigmaPt);
      fHistograms->Add(fHistMCPhotonPt);
      fHistograms->Add(fHistMCPhotonSigmaPt);

      fHistMCTruthPhotonP = new TH1F("fHistMCTruthPhotonP",
                                      "; #it{p} (GeV/#it{c}); Entries", 200,
                                      0, 10);
      fHistMCTruthPhotonSigmaP = new TH1F("fHistMCTruthPhotonSigmaP",
                                           "; #it{p} (GeV/#it{c}); Entries",
                                           200, 0, 10);
      fHistMCPhotonP = new TH1F("fHistMCPhotonP",
                                 "; #it{p} (GeV/#it{c}); Entries", 200, 0,
                                 10);
      fHistMCPhotonSigmaP = new TH1F("fHistMCPhotonSigmaP",
                                      "; #it{p} (GeV/#it{c}); Entries", 200,
                                      0, 10);
      fHistograms->Add(fHistMCTruthPhotonP);
      fHistograms->Add(fHistMCTruthPhotonSigmaP);
      fHistograms->Add(fHistMCPhotonP);
      fHistograms->Add(fHistMCPhotonSigmaP);

      fHistMCPhotonSource = new TH2F(
          "fHistMCPhotonSource", "; #it{p} (GeV/#it{c}); PDG code mother", 200, 0,
          10, 4000, 0, 4000);
      fHistograms->Add(fHistMCPhotonSource);
    }

    fHistSingleParticleCuts[0] = new TH1F("fHistSingleParticleCuts_pos",
                                          ";;Entries", 11, 0, 11);
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "#chi^{2}");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(4, "nCls TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(5, "nCrossed Rows TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(6,
                                                        "Ratio Findable TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(7, "nCls Findable");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(9,
                                                        "Daughter DCA to PV");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(10, "TPC refit");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(11, "Kink");
    fHistogramsPos->Add(fHistSingleParticleCuts[0]);

    fHistSingleParticleCuts[1] = new TH1F("fHistSingleParticleCuts_neg",
                                          ";;Entries", 11, 0, 11);
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(2, "#eta");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(3, "#chi^{2}");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(4, "nCls TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(5, "nCrossed Rows TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(6,
                                                        "Ratio Findable TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(7, "nCls Findable");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(9,
                                                        "Daughter DCA to PV");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(10, "TPC refit");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(11, "Kink");
    fHistogramsNeg->Add(fHistSingleParticleCuts[1]);

    if (fHistogramsBefore != nullptr) {
      delete fHistogramsBefore;
      fHistogramsBefore = nullptr;
    }
    if (fHistogramsBefore == nullptr) {
      fHistogramsBefore = new TList();
      fHistogramsBefore->SetOwner(kTRUE);
      fHistogramsBefore->SetName("Before");
    }

    if (fHistogramsAfter != nullptr) {
      delete fHistogramsAfter;
      fHistogramsAfter = nullptr;
    }
    if (fHistogramsAfter == nullptr) {
      fHistogramsAfter = new TList();
      fHistogramsAfter->SetOwner(kTRUE);
      fHistogramsAfter->SetName("After");
    }

    fHistDecayVertexXBefore = new TH2F(
        "fHistDecayVertexXBefore",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsBefore->Add(fHistDecayVertexXBefore);

    fHistDecayVertexYBefore = new TH2F(
        "fHistDecayVertexYBefore",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsBefore->Add(fHistDecayVertexYBefore);

    fHistDecayVertexZBefore = new TH2F(
        "fHistDecayVertexZBefore",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsBefore->Add(fHistDecayVertexZBefore);

    fHistDecayVertexXAfter = new TH2F(
        "fHistDecayVertexXAfter",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsAfter->Add(fHistDecayVertexXAfter);

    fHistDecayVertexYAfter = new TH2F(
        "fHistDecayVertexYAfter",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsAfter->Add(fHistDecayVertexYAfter);

    fHistDecayVertexZAfter = new TH2F(
        "fHistDecayVertexZAfter",
        "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsAfter->Add(fHistDecayVertexZAfter);

    fHistTransverseRadiusBefore = new TH2F(
        "fHistTransverseRadiusBefore",
        "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsBefore->Add(fHistTransverseRadiusBefore);

    fHistTransverseRadiusAfter = new TH2F(
        "fHistTransverseRadiusAfter",
        "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 50, 0, 10, 200, 0,
        200);
    fHistogramsAfter->Add(fHistTransverseRadiusAfter);

    fHistCosPABefore = new TH2F("fHistCosPABefore",
                                "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 100,
                                0, 10, 100, 0.95, 1);
    fHistogramsBefore->Add(fHistCosPABefore);

    fHistCosPAAfter = new TH2F("fHistCosPAAfter",
                               "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 50, 0,
                               10, 100, 0.95, 1);
    fHistogramsAfter->Add(fHistCosPAAfter);

    fHistDCArBefore = new TH2F(
        "fHistDCArBefore",
        "; #it{p}_{T} (GeV/#it{c}); DCA_{r} (cm)", 50, 0,
        10, 100, 0, 5);
    fHistogramsBefore->Add(fHistDCArBefore);

    fHistDCArAfter = new TH2F(
        "fHistDCArAfter",
        "; #it{p}_{T} (GeV/#it{c}); DCA_{r} (cm)", 50, 0,
        10, 100, 0, 5);
    fHistogramsAfter->Add(fHistDCArAfter);

    fHistDCAzBefore = new TH2F(
        "fHistDCAzBefore",
        "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", 50, 0,
        10, 200, -5, 5);
    fHistogramsBefore->Add(fHistDCAzBefore);

    fHistDCAzAfterOthersBefore = new TH2F(
        "fHistDCAzAfterOthersBefore",
        "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", 50, 0,
        10, 200, -5, 5);
    fHistogramsAfter->Add(fHistDCAzAfterOthersBefore);

    for (int i = 0; i < 4; ++i) {
      fHistDCAzAfterOthersBeforeQuality[i] = new TH2F(
          Form("fHistDCAzAfterOthersBefore_%i", i),
          "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", 50, 0, 10, 200, -5, 5);
      fHistogramsAfter->Add(fHistDCAzAfterOthersBeforeQuality[i]);
    }

    fHistDCAzAfter = new TH2F(
        "fHistDCAzAfter",
        "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", 50, 0,
        10, 200, -5, 5);
    fHistogramsAfter->Add(fHistDCAzAfter);

    fHistDCA = new TH2F("fHistDCA", "; #it{p}_{T} (GeV/#it{c}); DCA to PV (cm)",
                        50, 0, 10, 100, 0, 5);
    fHistogramsAfter->Add(fHistDCA);

    fHistDecayLength = new TH2F("fHistDecayLength",
                                "; #it{p}_{T} (GeV/#it{c}); Decay length (cm)",
                                50, 0, 10, 200, 0, 200);
    fHistogramsAfter->Add(fHistDecayLength);

    fHistArmenterosBefore = new TH2F("fHistArmenterosBefore",
                                     " ; #alpha; #it{q}_{T} (GeV/#it{c})", 200,
                                     -1, 1, 100, 0, 0.25);
    fHistogramsBefore->Add(fHistArmenterosBefore);

    fHistArmenterosAfter = new TH2F("fHistArmenterosAfter",
                                    " ; #alpha; #it{q}_{T} (GeV/#it{c})", 200,
                                    -1, 1, 100, 0, 0.25);
    fHistogramsAfter->Add(fHistArmenterosAfter);

    fHistPsiPairBefore = new TH2F("fHistPsiPairBefore",
                                  "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}", 100,
                                  0, 10, 200, -pi / 2.f, pi / 2.f);
    fHistogramsBefore->Add(fHistPsiPairBefore);

    fHistPsiPairAfter = new TH2F("fHistPsiPairAfter",
                                 "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}", 100,
                                 0, 10, 200, -pi / 2.f, pi / 2.f);
    fHistogramsAfter->Add(fHistPsiPairAfter);

    fHistQualityBefore = new TH2F("fHistQualityBefore",
                                  "; #it{p}_{T} (GeV/#it{c}); Photon quality",
                                  50, 0, 10, 4, 0, 4);
    fHistogramsBefore->Add(fHistQualityBefore);
    fHistQualityAfter = new TH2F("fHistQualityAfter",
                                 "; #it{p}_{T} (GeV/#it{c}); Photon quality",
                                 50, 0, 10, 4, 0, 4);
    fHistogramsAfter->Add(fHistQualityAfter);

    fHistTomography = new TH2F("fHistTomography", "; #it{x} (cm); #it{y} (cm)",
                               2400, -120, 120, 2400, -120, 120);
    fHistogramsAfter->Add(fHistTomography);

    fHistograms->Add(fHistogramsBefore);
    fHistograms->Add(fHistogramsAfter);

    fHistSingleParticlePt[0] = new TH1F("fHistSingleParticlePt_pos",
                                        ";#it{p}_{T}; Entries", 500, 0, 10);
    fHistSingleParticlePt[1] = new TH1F("fHistSingleParticlePt_neg",
                                        ";#it{p}_{T}; Entries", 500, 0, 10);
    fHistSingleParticleEtaBefore[0] = new TH2F(
        "fHistSingleParticleEtaBefore_pos", "; #it{p}_{T} (GeV/#it{c}); #eta",
        50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaBefore[1] = new TH2F(
        "fHistSingleParticleEtaBefore_neg", "; #it{p}_{T} (GeV/#it{c}); #eta",
        50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaAfter[0] = new TH2F("fHistSingleParticleEtaAfter_pos",
                                              "; #it{p}_{T} (GeV/#it{c}); #eta",
                                              50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaAfter[1] = new TH2F("fHistSingleParticleEtaAfter_neg",
                                              "; #it{p}_{T} (GeV/#it{c}); #eta",
                                              50, 0, 10, 100, -2, 2);
    fHistSingleParticleChi2Before[0] = new TH2F(
        "fHistSingleParticleChi2Before_pos",
        "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2Before[1] = new TH2F(
        "fHistSingleParticleChi2Before_neg",
        "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2After[0] = new TH2F(
        "fHistSingleParticleChi2After_pos",
        "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2After[1] = new TH2F(
        "fHistSingleParticleChi2After_neg",
        "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleNclsTPCBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCFindableBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCFindableBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10, 160, 0,
        160);
    fHistSingleParticleNclsTPCFindableBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCFindableBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10, 160, 0,
        160);
    fHistSingleParticleNclsTPCFindableAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCFindableAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10, 160, 0,
        160);
    fHistSingleParticleNclsTPCFindableAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCFindableAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10, 160, 0,
        160);
    fHistSingleParticleNclsTPCRatioFindableBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNcrossedTPCBefore[0] = new TH2F(
        "fHistSingleParticleNcrossedTPCBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCBefore[1] = new TH2F(
        "fHistSingleParticleNcrossedTPCBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCAfter[0] = new TH2F(
        "fHistSingleParticleNcrossedTPCAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCAfter[1] = new TH2F(
        "fHistSingleParticleNcrossedTPCAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleDCAtoPVBefore[0] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVBefore[1] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVAfter[0] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVAfter[1] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticlePileUp[0] = new TH2F(
        "fHistSingleParticlePileUp_pos",
        "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 50, 0, 10);
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(4, "None");
    fHistSingleParticlePileUp[1] = new TH2F(
        "fHistSingleParticlePileUp_neg",
        "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 50, 0, 10);
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(4, "None");

    fHistSingleParticlePID[0] = new TH2F(
        "fHistSingleParticlePID_pos", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        50, 0, 10, 250, -10, 10);
    fHistSingleParticlePID[1] = new TH2F(
        "fHistSingleParticlePID_neg", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        50, 0, 10, 250, -10, 10);

    fHistogramsPos->Add(fHistSingleParticlePt[0]);
    fHistogramsNeg->Add(fHistSingleParticlePt[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleChi2Before[0]);
    fHistogramsNeg->Add(fHistSingleParticleChi2Before[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleChi2After[0]);
    fHistogramsNeg->Add(fHistSingleParticleChi2After[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVAfter[1]);
    fHistogramsPos->Add(fHistSingleParticlePileUp[0]);
    fHistogramsNeg->Add(fHistSingleParticlePileUp[1]);
    fHistogramsPos->Add(fHistSingleParticlePID[0]);
    fHistogramsNeg->Add(fHistSingleParticlePID[1]);
  }
  fHistograms->Add(fHistogramsPos);
  fHistograms->Add(fHistogramsNeg);
}
