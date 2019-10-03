#include "AliAnalysisTaskHypertriton3ML.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVVertex.h"
#include "AliVertexerHyperTriton3Body.h"

#include <Riostream.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom3.h>

#include <array>
#include <map>
#include <stdio.h>
#include <unordered_map>
#include <vector>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHypertriton3ML);

namespace {

bool IsHyperTriton3(const AliVParticle *vPart, AliMCEvent *mcEvent) {
  int nDaughters = 0;

  int vPartPDG   = vPart->PdgCode();
  int vPartLabel = vPart->GetLabel();

  if (!mcEvent->IsPhysicalPrimary(vPartLabel) || (std::abs(vPartPDG) != 1010010030)) return false;

  for (int iD = vPart->GetDaughterFirst(); iD <= vPart->GetDaughterLast(); iD++) {
    AliVParticle *dPart = mcEvent->GetTrack(iD);

    int dPartPDG = dPart->PdgCode();
    if (std::abs(dPartPDG) != 11) nDaughters++;
  }
  if (nDaughters == 3) return true;
  return false;
}

int IsTrueHyperTriton3Candidate(AliESDtrack *t1, AliESDtrack *t2, AliESDtrack *t3, AliMCEvent *mcEvent) {
  if (!mcEvent) return 0;

  int lab1 = std::abs(t1->GetLabel());
  int lab2 = std::abs(t2->GetLabel());
  int lab3 = std::abs(t3->GetLabel());

  if (mcEvent->IsPhysicalPrimary(lab1)) return -1;
  if (mcEvent->IsPhysicalPrimary(lab2)) return -1;
  if (mcEvent->IsPhysicalPrimary(lab3)) return -1;

  AliVParticle *part1 = mcEvent->GetTrack(lab1);
  AliVParticle *part2 = mcEvent->GetTrack(lab2);
  AliVParticle *part3 = mcEvent->GetTrack(lab3);

  if (!part1) return -1;
  if (!part2) return -1;
  if (!part3) return -1;

  int mom1 = part1->GetMother();
  int mom2 = part2->GetMother();
  int mom3 = part3->GetMother();

  if (mom1 != mom2 || mom1 != mom3 || mom2 != mom3) return -1;

  AliVParticle *mom = mcEvent->GetTrack(mom1);
  if (!mom) return -1;

  return (IsHyperTriton3(mom, mcEvent)) ? mom1 : -1;
}

bool HasTOF(AliVTrack *track) {
  const bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const float len       = track->GetIntegratedLength();
  return hasTOFout && hasTOFtime && (len > 350.);
}

} // namespace

AliAnalysisTaskHypertriton3ML::AliAnalysisTaskHypertriton3ML(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fVertexer{}, fListHist{nullptr}, fTreeHyp3{nullptr},
      fInputHandler{nullptr}, fPIDResponse{nullptr}, fMC{mc}, fOnlyTrueCandidates{false}, fDownscaling{false},
      fHistNSigmaDeu{nullptr}, fHistNSigmaP{nullptr}, fHistNSigmaPi{nullptr}, fHistInvMass{nullptr},
      fDownscalingFactorByEvent{1.}, fDownscalingFactorByCandidate{1.}, fMinCanidatePtToSave{0.1},
      fMaxCanidatePtToSave{100.}, fMinITSNcluster{0}, fMinTPCNcluster{70}, fMaxNSigmaTPCDeu{5.}, fMaxNSigmaTPCP{5.},
      fMaxNSigmaTPCPi{5.}, fMaxNSigmaTOFDeu{5.}, fMaxNSigmaTOFP{5.}, fMaxNSigmaTOFPi{5.},
      fVertexerToleranceGuessCompatibility{0}, fVertexerMaxDistanceInit{100.}, fMinCosPA{0.993},
      fMinDCA2PrimaryVtxDeu{0.025}, fMinDCA2PrimaryVtxP{0.025}, fMinDCA2PrimaryVtxPi{0.05}, fMaxPtPion{1.},
      fSHypertriton{}, fRHypertriton{}, fREvent{}, fDeuVector{}, fPVector{}, fPiVector{} {

  // Settings for the custom vertexer
  fVertexer.SetToleranceGuessCompatibility(fVertexerToleranceGuessCompatibility);
  fVertexer.SetMaxDinstanceInit(fVertexerMaxDistanceInit);

  // Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); // Basic Histograms
  DefineOutput(2, TTree::Class()); // Hypertriton Candidates Tree output
}

AliAnalysisTaskHypertriton3ML::~AliAnalysisTaskHypertriton3ML() {
  if (fListHist) {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3) {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }
}

void AliAnalysisTaskHypertriton3ML::UserCreateOutputObjects() {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler           = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse            = fInputHandler->GetPIDResponse();

  fInputHandler->SetNeedField();

  fTreeHyp3 = new TTree("fHypertritonTree", "Hypertriton3 Candidates");
  fTreeHyp3->Branch("REvent", &fREvent);
  fTreeHyp3->Branch("RHypertriton", &fRHypertriton);
  if (man->GetMCtruthEventHandler()) fTreeHyp3->Branch("SHypertriton", &fSHypertriton);

  fListHist = new TList();
  fListHist->SetOwner();
  fEventCuts.AddQAplotsToList(fListHist);

  fHistNSigmaDeu = new TH2D("fHistNSigmaDeu", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Deuteron; Counts", 100, 0., 10.,
                            80, -5.0, 5.0);
  fHistNSigmaP =
      new TH2D("fHistNSigmaP", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Proton; Counts", 100, 0., 10., 80, -5.0, 5.0);
  fHistNSigmaPi =
      new TH2D("fHistNSigmaPi", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Pion; Counts", 100, 0., 10., 80, -5.0, 5.0);

  fHistInvMass = new TH2D("fHistInvMass", ";m_{dp#pi}(GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c}); Counts", 30, 2.96,
                          3.05, 100, 0, 10);

  fListHist->Add(fHistNSigmaDeu);
  fListHist->Add(fHistNSigmaP);
  fListHist->Add(fHistNSigmaPi);

  fListHist->Add(fHistInvMass);

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();
} // end UserCreateOutputObjects

void AliAnalysisTaskHypertriton3ML::UserExec(Option_t *) {
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskHypertriton3ML::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC) {
    ::Fatal("AliAnalysisTaskHypertriton3ML::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  fREvent.fCent = fEventCuts.GetCentrality();

  double primaryVtxPos[3];
  fEventCuts.GetPrimaryVertex()->GetXYZ(primaryVtxPos);

  fREvent.fX = primaryVtxPos[0];
  fREvent.fY = primaryVtxPos[1];
  fREvent.fZ = primaryVtxPos[2];

  unsigned char tgr = 0x0;

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) tgr = 1;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) tgr = 2;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) tgr = 4;

  fREvent.fTrigger = tgr;

  std::unordered_map<int, int> mcMap;

  if (fMC) {
    fSHypertriton.clear();

    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliVParticle *part = mcEvent->GetTrack(iTrack);

      if (!part) {
        ::Warning("AliAnalysisTaskHypertriton3ML::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }

      if (std::abs(part->Y()) > 1.) continue;
      if (!IsHyperTriton3(part, mcEvent)) continue;

      double decayVtx[3]{0.0, 0.0, 0.0};

      AliVParticle *deu{nullptr}, *p{nullptr}, *pi{nullptr};

      for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD) {
        AliVParticle *daughter = mcEvent->GetTrack(iD);

        if (mcEvent->IsSecondaryFromWeakDecay(iD) && daughter) {
          decayVtx[0] = daughter->Xv();
          decayVtx[1] = daughter->Yv();
          decayVtx[2] = daughter->Zv();

          if (std::abs(daughter->PdgCode()) == AliPID::ParticleCode(AliPID::kDeuteron)) deu = daughter;
          if (std::abs(daughter->PdgCode()) == AliPID::ParticleCode(AliPID::kProton)) p = daughter;
          if (std::abs(daughter->PdgCode()) == AliPID::ParticleCode(AliPID::kPion)) pi = daughter;
        }
      }
      if (deu == nullptr) continue;

      SHypertriton3 hyp3s;

      hyp3s.fPdgCode = part->PdgCode();

      hyp3s.fDecayVtxX = decayVtx[0] - part->Xv();
      hyp3s.fDecayVtxY = decayVtx[1] - part->Yv();
      hyp3s.fDecayVtxZ = decayVtx[2] - part->Zv();

      hyp3s.fPxDeu = deu->Px();
      hyp3s.fPyDeu = deu->Py();
      hyp3s.fPzDeu = deu->Pz();
      hyp3s.fPxP   = p->Px();
      hyp3s.fPyP   = p->Py();
      hyp3s.fPzP   = p->Pz();
      hyp3s.fPxPi  = pi->Px();
      hyp3s.fPyPi  = pi->Py();
      hyp3s.fPzPi  = pi->Pz();

      hyp3s.fRecoIndex = -1;

      mcMap[iTrack] = fSHypertriton.size();
      fSHypertriton.push_back(hyp3s);
    }
  }

  fRHypertriton.clear();

  fDeuVector.clear();
  fPVector.clear();
  fPiVector.clear();

  double b = esdEvent->GetMagneticField();

  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if (!track) continue;

    if (((track->GetStatus() & AliVTrack::kTPCrefit) == 0 && (track->GetStatus() & AliVTrack::kITSrefit) == 0) ||
        (track->GetKinkIndex(0) > 0) || (track->GetITSNcls() < fMinITSNcluster) ||
        (track->GetTPCNcls() < fMinTPCNcluster))
      continue;

    float dca[2];
    track->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    if (fMC && fOnlyTrueCandidates) {
      int lab = std::abs(track->GetLabel());
      if (!mcEvent->IsSecondaryFromWeakDecay(lab)) continue;
      AliVParticle *part = mcEvent->GetTrack(lab);
      AliVParticle *moth = mcEvent->GetTrack(part->GetMother());
      if (std::abs(moth->PdgCode()) != 1010010030) continue;
    }

    float nSigmaTPCDeu = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    if (std::abs(nSigmaTPCDeu) < fMaxNSigmaTPCDeu && dcaNorm > fMinDCA2PrimaryVtxDeu) {
      if (HasTOF(track)) {
        if (std::abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron)) < fMaxNSigmaTOFDeu)
          fDeuVector.push_back(track);
      } else {
        fDeuVector.push_back(track);
      }
    }

    float nSigmaTPCP = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    if (std::abs(nSigmaTPCP) < fMaxNSigmaTPCP && dcaNorm > fMinDCA2PrimaryVtxP) {
      if (HasTOF(track)) {
        if (std::abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) < fMaxNSigmaTOFP)
          fPVector.push_back(track);
      } else {
        fPVector.push_back(track);
      }
    }

    float nSigmaTPCPi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    if (std::abs(nSigmaTPCPi) < fMaxNSigmaTPCPi && dcaNorm > fMinDCA2PrimaryVtxPi && track->Pt() < fMaxPtPion) {
      if (HasTOF(track)) {
        if (std::abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion)) < fMaxNSigmaTOFPi)
          fPiVector.push_back(track);
      } else {
        fPiVector.push_back(track);
      }
    }
  }

  // downscaling the output tree saving only a fDownscalingFactorByEvent fraction of the events
  if (!fMC && fDownscaling && ((int)fDeuVector.size() > 0)) {
    if (gRandom->Rndm() > fDownscalingFactorByEvent) return;
  }

  for (const auto &deu : fDeuVector) {
    float nSigmaDeu = fPIDResponse->NumberOfSigmasTPC(deu, AliPID::kDeuteron);

    for (const auto &p : fPVector) {
      if (deu == p) continue;
      if (p->Charge() * deu->Charge() < 0) continue;

      float nSigmaP = fPIDResponse->NumberOfSigmasTPC(p, AliPID::kProton);

      for (const auto &pi : fPiVector) {
        if (p == pi || deu == pi) continue;
        if (pi->Charge() * p->Charge() > 0) continue;

        float nSigmaPi = fPIDResponse->NumberOfSigmasTPC(pi, AliPID::kPion);

        int momLab = IsTrueHyperTriton3Candidate(deu, p, pi, mcEvent);
        if ((momLab == -1) && fOnlyTrueCandidates) continue;

        bool recoVertex = fVertexer.FindDecayVertex(deu, p, pi, b);
        if (!recoVertex) continue;

        AliESDVertex *decayVtx = static_cast<AliESDVertex *>(fVertexer.GetCurrentVertex());

        double decayLenght[3] = {0., 0., 0.};
        decayVtx->GetXYZ(decayLenght);

        for (int i = 0; i < 3; i++) {
          decayLenght[i] -= primaryVtxPos[i];
        }

        double dcaDecayDeu[2], dcaDecayP[2], dcaDecayPi[2];

        deu->PropagateToDCA(decayVtx, b, 1000., dcaDecayDeu);
        p->PropagateToDCA(decayVtx, b, 1000., dcaDecayP);
        pi->PropagateToDCA(decayVtx, b, 1000., dcaDecayPi);

        LVector_t deu4Vector, p4Vector, pi4Vector, hyp4Vector;

        deu4Vector.SetCoordinates(deu->Px(), deu->Py(), deu->Pz(), AliPID::ParticleMass(AliPID::kDeuteron));
        p4Vector.SetCoordinates(p->Px(), p->Py(), p->Pz(), AliPID::ParticleMass(AliPID::kProton));
        pi4Vector.SetCoordinates(pi->Px(), pi->Py(), pi->Pz(), AliPID::ParticleMass(AliPID::kPion));

        hyp4Vector = deu4Vector + p4Vector + pi4Vector;

        float hypPt = hyp4Vector.Pt();
        float hypM  = hyp4Vector.M();

        if ((hypPt < fMinCanidatePtToSave) || (fMaxCanidatePtToSave < hypPt)) continue;
        if (hypM < 2.9 || hypM > 3.2) continue;

        double dTotHyper = std::sqrt(decayLenght[0] * decayLenght[0] + decayLenght[1] * decayLenght[1] +
                                     decayLenght[2] * decayLenght[2]);

        double cosPA =
            hyp4Vector.Px() * decayLenght[0] + hyp4Vector.Py() * decayLenght[1] + hyp4Vector.Pz() * decayLenght[2];

        cosPA /= (dTotHyper * hyp4Vector.P());
        if (cosPA < fMinCosPA) continue;

        if (fMC) {
          fSHypertriton[mcMap[momLab]].fRecoIndex = (fRHypertriton.size());
        }

        // downscaling the output tree saving only a fDownscalingFactorByCandidate fraction of the candidates
        if (!fMC && fDownscaling) {
          if (gRandom->Rndm() > fDownscalingFactorByCandidate) continue;
        }

        RHypertriton3 hyp3r;

        hyp3r.fDecayVtxX = decayVtx->GetX();
        hyp3r.fDecayVtxY = decayVtx->GetY();
        hyp3r.fDecayVtxZ = decayVtx->GetZ();

        hyp3r.fPxDeu = deu->Px();
        hyp3r.fPyDeu = deu->Py();
        hyp3r.fPzDeu = deu->Pz();
        hyp3r.fPxP   = p->Px();
        hyp3r.fPyP   = p->Py();
        hyp3r.fPzP   = p->Pz();
        hyp3r.fPxPi  = pi->Px();
        hyp3r.fPyPi  = pi->Py();
        hyp3r.fPzPi  = pi->Pz();

        hyp3r.fPosXDeu = deu->GetX();
        hyp3r.fPosYDeu = deu->GetY();
        hyp3r.fPosZDeu = deu->GetZ();
        hyp3r.fPosXP   = p->GetX();
        hyp3r.fPosYP   = p->GetY();
        hyp3r.fPosZP   = p->GetZ();
        hyp3r.fPosXPi  = pi->GetX();
        hyp3r.fPosYPi  = pi->GetY();
        hyp3r.fPosZPi  = pi->GetZ();

        hyp3r.fDCAxyDeu = std::abs(dcaDecayDeu[0]);
        hyp3r.fDCAzDeu  = std::abs(dcaDecayDeu[1]);
        hyp3r.fDCAxyP   = std::abs(dcaDecayP[0]);
        hyp3r.fDCAzP    = std::abs(dcaDecayP[1]);
        hyp3r.fDCAxyPi  = std::abs(dcaDecayPi[0]);
        hyp3r.fDCAzPi   = std::abs(dcaDecayPi[1]);

        hyp3r.fNClusterTPCDeu = deu->GetTPCNcls();
        hyp3r.fNClusterTPCP   = p->GetTPCNcls();
        hyp3r.fNClusterTPCPi  = pi->GetTPCNcls();

        hyp3r.fITSClusterMapDeu = deu->GetITSClusterMap();
        hyp3r.fITSClusterMapP   = p->GetITSClusterMap();
        hyp3r.fITSClusterMapPi  = pi->GetITSClusterMap();

        hyp3r.fNSigmaTPCDeu = nSigmaDeu;
        hyp3r.fNSigmaTPCP   = nSigmaP;
        hyp3r.fNSigmaTPCPi  = nSigmaPi;

        hyp3r.fHasTOFDeu = HasTOF(deu);
        hyp3r.fHasTOFP   = HasTOF(p);
        hyp3r.fHasTOFPi  = HasTOF(pi);

        hyp3r.fNSigmaTOFDeu = fPIDResponse->NumberOfSigmasTOF(deu, AliPID::kDeuteron);
        hyp3r.fNSigmaTOFP   = fPIDResponse->NumberOfSigmasTOF(p, AliPID::kProton);
        hyp3r.fNSigmaTOFPi  = fPIDResponse->NumberOfSigmasTOF(pi, AliPID::kPion);

        hyp3r.fTrackChi2Deu       = deu->GetTPCchi2() / (deu->GetTPCNcls() + 1.e-16);
        hyp3r.fTrackChi2P         = p->GetTPCchi2() / (p->GetTPCNcls() + 1.e-16);
        hyp3r.fTrackChi2Pi        = pi->GetTPCchi2() / (pi->GetTPCNcls() + 1.e-16);
        hyp3r.fDecayVertexChi2NDF = decayVtx->GetChi2perNDF();

        hyp3r.fIsMatter = deu->Charge() > 0;

        fRHypertriton.push_back(hyp3r);

        fHistNSigmaDeu->Fill(deu->Pt(), nSigmaDeu);
        fHistNSigmaP->Fill(p->Pt(), nSigmaP);
        fHistNSigmaPi->Fill(pi->Pt(), nSigmaPi);

        fHistInvMass->Fill(hypM, hypPt);
      }
    }
  }

  if (fMC) {
    fTreeHyp3->Fill();

    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
  }
  if (!fMC && (int)fRHypertriton.size() > 0) {
    fTreeHyp3->Fill();

    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
  }
}

void AliAnalysisTaskHypertriton3ML::Terminate(Option_t *) {}

AliAnalysisTaskHypertriton3ML *AliAnalysisTaskHypertriton3ML::AddTask(bool isMC, TString suffix) {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHyperTriton2BodyML", "No analysis manager found.");
    return nullptr;
  }
  mgr->SetDebugLevel(2);

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHyperTriton2BodyML", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHypertriton3ML";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHypertriton3ML *task = new AliAnalysisTaskHypertriton3ML(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("HyperTritonTree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
