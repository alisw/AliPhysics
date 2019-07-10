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
      fInputHandler{nullptr}, fPIDResponse{nullptr}, fMC{mc}, fOnlyTrueCandidates{false}, fHistNSigmaDeu{nullptr},
      fHistNSigmaP{nullptr}, fHistNSigmaPi{nullptr}, fHistInvMass{nullptr}, fHistTPCdEdx{nullptr},
      fMinCanidatePtToSave{0.1}, fMaxCanidatePtToSave{100.}, fMinITSNcluster{1}, fMinTPCNcluster{70},
      fMaxNSigmaTPCDeu{5.0}, fMaxNSigmaTPCP{5.0}, fMaxNSigmaTPCPi{5.0}, fMinCosPA{0.9}, fMinDCA2PrimaryVtxDeu{0.005},
      fMinDCA2PrimaryVtxP{0.005}, fMinDCA2PrimaryVtxPi{0.01}, fCurrentFileName{""}, fSHypertriton{},
      fRHypertriton{}, fREvent{}, fDeuVector{}, fPVector{}, fPiVector{} {

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

  fHistInvMass = new TH2D("fHistInvMass", ";#it{p}_{T} (GeV/#it{c});Invariant Mass(GeV/#it{c^2}); Counts", 100, 0, 10,
                          30, 2.96, 3.05);

  fHistTPCdEdx[0] = new TH2D("fHistTPCdEdxM", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);
  fHistTPCdEdx[1] = new TH2D("fHistTPCdEdxA", ";#it{p} (GeV/#it{c}); dE/dx; Counts", 256, 0, 10.24, 4096, 0, 2048);

  fListHist->Add(fHistNSigmaDeu);
  fListHist->Add(fHistNSigmaP);
  fListHist->Add(fHistNSigmaPi);

  fListHist->Add(fHistInvMass);

  fListHist->Add(fHistTPCdEdx[0]);
  fListHist->Add(fHistTPCdEdx[1]);

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
  }

  fREvent.fCent = fEventCuts.GetCentrality();

  const AliVVertex *primaryVtx = static_cast<const AliVVertex *>(fEventCuts.GetPrimaryVertex());

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

    float nSigmaDeu = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    if (std::abs(nSigmaDeu) < fMaxNSigmaTPCDeu) {
      double dcaPrimaryDeu[2];
      track->PropagateToDCA(primaryVtx, b, 1000., dcaPrimaryDeu);
      double dcaPrimaryDeuNorm = std::sqrt(dcaPrimaryDeu[0] * dcaPrimaryDeu[0] + dcaPrimaryDeu[1] * dcaPrimaryDeu[1]);
      if (dcaPrimaryDeuNorm < fMinDCA2PrimaryVtxDeu) {
        fDeuVector.push_back(track);
      }
    }

    float nSigmaP = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    if (std::abs(nSigmaP) < fMaxNSigmaTPCP) {
      double dcaPrimaryP[2];
      track->PropagateToDCA(primaryVtx, b, 1000., dcaPrimaryP);
      double dcaPrimaryPNorm = std::sqrt(dcaPrimaryP[0] * dcaPrimaryP[0] + dcaPrimaryP[1] * dcaPrimaryP[1]);
      if (dcaPrimaryPNorm < fMinDCA2PrimaryVtxP) {
        fPVector.push_back(track);
      }
    }

    float nSigmaPi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    if (std::abs(nSigmaPi) < fMaxNSigmaTPCPi) {
      double dcaPrimaryPi[2];
      track->PropagateToDCA(primaryVtx, b, 1000., dcaPrimaryPi);
      double dcaPrimaryPiNorm = std::sqrt(dcaPrimaryPi[0] * dcaPrimaryPi[0] + dcaPrimaryPi[1] * dcaPrimaryPi[1]);
      if (dcaPrimaryPiNorm < fMinDCA2PrimaryVtxPi) {
        fPiVector.push_back(track);
      }
    }
  }

  for (const auto &deu : fDeuVector) {
    float nSigmaDeu = fPIDResponse->NumberOfSigmasTPC(deu, AliPID::kDeuteron);

    for (const auto &p : fPVector) {
      if (p->Charge() * deu->Charge() < 0) continue;
      float nSigmaP = fPIDResponse->NumberOfSigmasTPC(p, AliPID::kProton);

      for (const auto &pi : fPiVector) {
        if (pi->Charge() * p->Charge() > 0) continue;
        float nSigmaPi = fPIDResponse->NumberOfSigmasTPC(pi, AliPID::kPion);

        int momLab = IsTrueHyperTriton3Candidate(deu, p, pi, mcEvent);
        if (momLab != -1) cout << momLab << endl;
        if ((momLab == -1) && fOnlyTrueCandidates) continue;

        bool recoVertex = fVertexer.FindDecayVertex(deu, p, pi, b);
        if (!recoVertex) continue;

        AliESDVertex *decayVtx = static_cast<AliESDVertex *>(fVertexer.GetCurrentVertex());

        double decayLenght[3] = {0., 0., 0.};
        decayVtx->GetXYZ(decayLenght);

        for (int i = 0; i < 3; i++) {
          decayLenght[i] -= primaryVtxPos[i];
        }

        double pDeu[3], pP[3], pPi[3], pHyper[3], dcaDecayDeu[2], dcaDecayP[2], dcaDecayPi[2];

        deu->PropagateToDCA(decayVtx, b, 100., dcaDecayDeu);
        p->PropagateToDCA(decayVtx, b, 100., dcaDecayP);
        pi->PropagateToDCA(decayVtx, b, 100., dcaDecayPi);

        deu->GetPxPyPz(pDeu);
        p->GetPxPyPz(pP);
        pi->GetPxPyPz(pPi);

        for (int i = 0; i < 3; i++) {
          pHyper[i] = pDeu[i] + pP[i] + pPi[i];
        }

        double dTotHyper = std::sqrt(decayLenght[0] * decayLenght[0] + decayLenght[1] * decayLenght[1] +
                                     decayLenght[2] * decayLenght[2]);
        double pTotHyper = std::sqrt(pHyper[0] * pHyper[0] + pHyper[1] * pHyper[1] + pHyper[2] * pHyper[2]);

        double cosPA = 0.;

        for (int i = 0; i < 3; i++) {
          cosPA += pHyper[i] * decayLenght[i];
        }

        cosPA /= (dTotHyper * pTotHyper);
        if (cosPA < fMinCosPA) continue;

        if (fMC) {
          fSHypertriton[mcMap[momLab]].fRecoIndex = (fRHypertriton.size());
        }

        RHypertriton3 hyp3r;

        hyp3r.fDecayVtxX = decayVtx->GetX();
        hyp3r.fDecayVtxY = decayVtx->GetY();
        hyp3r.fDecayVtxZ = decayVtx->GetZ();

        hyp3r.fPxDeu = pDeu[0];
        hyp3r.fPyDeu = pDeu[1];
        hyp3r.fPzDeu = pDeu[2];
        hyp3r.fPxP   = pP[0];
        hyp3r.fPyP   = pP[1];
        hyp3r.fPzP   = pP[2];
        hyp3r.fPxPi  = pPi[0];
        hyp3r.fPyPi  = pPi[1];
        hyp3r.fPzPi  = pPi[2];

        hyp3r.fPosXDeu = deu->GetX();
        hyp3r.fPosYDeu = p->GetY();
        hyp3r.fPosZDeu = pi->GetZ();
        hyp3r.fPosXP   = deu->GetX();
        hyp3r.fPosYP   = p->GetY();
        hyp3r.fPosZP   = pi->GetZ();
        hyp3r.fPosXPi  = deu->GetX();
        hyp3r.fPosYPi  = p->GetY();
        hyp3r.fPosZPi  = pi->GetZ();

        hyp3r.fDCAxyDeu = dcaDecayDeu[0];
        hyp3r.fDCAzDeu  = dcaDecayDeu[1];
        hyp3r.fDCAxyP   = dcaDecayP[0];
        hyp3r.fDCAzP    = dcaDecayP[1];
        hyp3r.fDCAxyPi  = dcaDecayPi[0];
        hyp3r.fDCAzPi   = dcaDecayPi[1];

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

        hyp3r.fTrackChi2Deu       = deu->GetTPCchi2();
        hyp3r.fTrackChi2P         = p->GetTPCchi2();
        hyp3r.fTrackChi2Pi        = pi->GetTPCchi2();
        hyp3r.fDecayVertexChi2NDF = decayVtx->GetChi2perNDF();

        hyp3r.fIsMatter = deu->Charge() > 0;

        fRHypertriton.push_back(hyp3r);
      }
    }
  }

  fTreeHyp3->Fill();

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
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
