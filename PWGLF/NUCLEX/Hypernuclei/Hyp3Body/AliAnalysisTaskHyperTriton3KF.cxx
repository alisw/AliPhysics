#include "AliAnalysisTaskHyperTriton3KF.h"

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVVertex.h"

#include <Riostream.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>

#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"

#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

ClassImp(AliAnalysisTaskHyperTriton3KF);

namespace {

struct HelperParticle {
  AliESDtrack* track = nullptr;
  float        nSigmaTPC = -1.f;
  float        nSigmaTOF = -1.f;
  KFParticle   particle = KFParticle();
};

constexpr float kHyperTritonMass{2.99131};

constexpr float kDeuMass{1.87561};
constexpr float kPMass{0.938272};
constexpr float kPiMass{0.13957};
constexpr float kMasses[3]{kDeuMass,kPMass,kPiMass};
constexpr AliPID::EParticleType kAliPID[3]{AliPID::kDeuteron,AliPID::kProton,AliPID::kPion};
const int kPDGs[3]{AliPID::ParticleCode(kAliPID[0]),AliPID::ParticleCode(kAliPID[1]),AliPID::ParticleCode(kAliPID[2])};

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

  if (!part1 || !part2 || !part3) 
    return -1;

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

/// helper functions
template <typename T> double Sq(T a) { return a * a; }
template <typename F> double Hypot(F a, F b, F c) { return std::sqrt(Sq(a) + Sq(b) + Sq(c)); }
template <typename F> double Hypot(F a, F b, F c, F d) { return std::sqrt(Sq(a) + Sq(b) + Sq(c) + Sq(d)); }

}    // namespace

AliAnalysisTaskHyperTriton3KF::AliAnalysisTaskHyperTriton3KF(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fMC{mc}, fREvent{}, fGenHyp{}, fRecHyp{} {
  fTrackCuts.SetMinNClustersTPC(0);
  fTrackCuts.SetEtaRange(-0.9,0.9);
  /// Settings for the custom vertexer

  /// Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());    // Basic Histograms
  DefineOutput(2, TTree::Class());    // Hypertriton Candidates Tree output
}

AliAnalysisTaskHyperTriton3KF::~AliAnalysisTaskHyperTriton3KF() {
  if (fListHist) {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3) {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }

}

void AliAnalysisTaskHyperTriton3KF::UserCreateOutputObjects() {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler           = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse            = fInputHandler->GetPIDResponse();

  fInputHandler->SetNeedField();

  fListHist = new TList();
  fListHist->SetOwner(true);
  fEventCuts.AddQAplotsToList(fListHist);

  fHistNSigmaDeu = new TH2D("fHistNSigmaDeu", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Deuteron; Counts", 100, 0., 10.,
                            80, -5.0, 5.0);
  fHistNSigmaP =
      new TH2D("fHistNSigmaP", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Proton; Counts", 100, 0., 10., 80, -5.0, 5.0);
  fHistNSigmaPi =
      new TH2D("fHistNSigmaPi", ";#it{p}_{T} (GeV/#it{c});n_{#sigma} TPC Pion; Counts", 100, 0., 10., 80, -5.0, 5.0);

  fHistInvMass =
      new TH2D("fHistInvMass", ";m_{dp#pi}(GeV/#it{c^2}); #it{p}_{T} (GeV/#it{c}); Counts", 30, 2.96, 3.05, 100, 0, 10);

  fListHist->Add(fHistNSigmaDeu);
  fListHist->Add(fHistNSigmaP);
  fListHist->Add(fHistNSigmaPi);

  fListHist->Add(fHistInvMass);

  OpenFile(2);
  fTreeHyp3 = new TTree("Hyp3KF","Hypetriton 3 Body with the KFParticle");
  fTreeHyp3->Branch("RCollision", &fREvent);
  fTreeHyp3->Branch("RHyperTriton", &fRecHyp);

  if (man->GetMCtruthEventHandler())
    fTreeHyp3->Branch("SHyperTriton", &fGenHyp);

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();

}    /// end UserCreateOutputObjects

void AliAnalysisTaskHyperTriton3KF::UserExec(Option_t *) {
    // set Magnetic field for KF
  KFParticle::SetField(fInputEvent->GetMagneticField());
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskHyperTriton3KF::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC) {
    ::Fatal("AliAnalysisTaskHyperTriton3KF::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  fGenHyp.clear();
  fRecHyp.clear();

  double pvPos[3], pvCov[6];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pvPos);
  fEventCuts.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);
  fREvent.fX = pvPos[0];
  fREvent.fY = pvPos[1];
  fREvent.fZ = pvPos[2];

  KFPVertex kfPVertex;
  kfPVertex.SetXYZ(pvPos[0],pvPos[1],pvPos[2]);
  kfPVertex.SetCovarianceMatrix(pvCov[0],pvCov[1],pvCov[2],pvCov[3],pvCov[4],pvCov[5]);
  kfPVertex.SetChi2(fEventCuts.GetPrimaryVertex()->GetChi2());
  kfPVertex.SetNDF(fEventCuts.GetPrimaryVertex()->GetNDF());
  kfPVertex.SetNContributors(fEventCuts.GetPrimaryVertex()->GetNContributors());

  KFParticle prodVertex{kfPVertex};

  fREvent.fTrigger = 0u;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) fREvent.fTrigger |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) fREvent.fTrigger |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) fREvent.fTrigger |= kSemiCentral;
  fREvent.fTrigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::unordered_map<int, int> mcMap;
  if (fMC) {
    double mcVtx[3];
    mcEvent->GetPrimaryVertex()->GetXYZ(mcVtx);
    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliVParticle *part = mcEvent->GetTrack(iTrack);

      if (!part) {
        ::Warning("AliAnalysisTaskHyperTriton3KF::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }

      if (std::abs(part->Y()) > 1.) continue;
      if (!IsHyperTriton3(part, mcEvent)) continue;

      double decayVtx[3]{0.0, 0.0, 0.0};

      for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD) {
        AliVParticle *daughter = mcEvent->GetTrack(iD);

        if (mcEvent->IsSecondaryFromWeakDecay(iD) && daughter && std::abs(daughter->PdgCode()) != 11) {
          decayVtx[0] = daughter->Xv();
          decayVtx[1] = daughter->Yv();
          decayVtx[2] = daughter->Zv();
          break;
        }
      }
      SHyperTriton3KF genHyp;
      genHyp.l = Hypot(mcVtx[0] - decayVtx[0], mcVtx[1] - decayVtx[1], mcVtx[2] - decayVtx[2]);
      genHyp.px = part->Px();
      genHyp.py = part->Py();
      genHyp.pz = part->Pz();
      mcMap[iTrack] = fGenHyp.size();
      fGenHyp.emplace_back(genHyp);
    }
  }

  std::vector<HelperParticle> helpers[3];

  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);
    if (!track) continue;

    if (!fTrackCuts.AcceptTrack(track)) continue;

    if (fMC && fOnlyTrueCandidates) {
      int lab = std::abs(track->GetLabel());
      if (!mcEvent->IsSecondaryFromWeakDecay(lab)) continue;
      AliVParticle *part = mcEvent->GetTrack(lab);
      AliVParticle *moth = mcEvent->GetTrack(part->GetMother());
      if (std::abs(moth->PdgCode()) != 1010010030) continue;
    }

    bool candidate[3]{false,false,false};
    float nSigmasTPC[3]{-1.,-1.,-1.}, nSigmasTOF[3]{-1.,-1.,-1.};
    bool hasTOF{HasTOF(track)};
    float dca[2];
    track->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    for (int iT{0}; iT < 3; ++iT) {
      nSigmasTPC[iT] = fPIDResponse->NumberOfSigmasTPC(track, kAliPID[iT]);
      nSigmasTOF[iT] = fPIDResponse->NumberOfSigmasTOF(track, kAliPID[iT]);
      if (std::abs(nSigmasTPC[iT]) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] && 
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        candidate[iT] = (std::abs(nSigmasTOF[iT]) < fTOFsigmas[iT]) || (!hasTOF && !fRequireTOFpid[kDeuteron]);
    }
  
    if (candidate[0] || candidate[1] || candidate[2]) {
      HelperParticle helper;
      helper.track = track;
      double posmom[6],cov[21];
      track->GetXYZ(posmom);
      track->GetPxPyPz(posmom+3);
      track->GetCovarianceXYZPxPyPz(cov);
      for (int iT{0}; iT < 3; ++iT) {
        if (candidate[iT]) {
          helper.particle.Create(posmom,cov,track->Charge(),kMasses[iT]);
          helper.particle.Chi2() = track->GetTPCchi2();
          helper.particle.NDF() = track->GetNumberOfTPCClusters() * 2;
          helper.nSigmaTPC = nSigmasTPC[iT];
          helper.nSigmaTOF = nSigmasTOF[iT];
          helpers[iT].push_back(helper);
        }
      }
    }
  }

  /// downscaling the output tree saving only a fDownscalingFactorByEvent fraction of the events
  if (!fMC && fDownscaling && ((int)helpers[kDeuteron].size() > 0)) {
    if (gRandom->Rndm() > fDownscalingFactorByEvent) return;
  }

  /// if event mixing is enabled takes deuteron from the event mixing pool
  // if (fEnableEventMixing && fApplyML) {
  //   deuterons = GetEventMixingTracks(fREvent.fCent, fREvent.fZ);
  // }

  RHyperTriton3KF recHyp;
  for (const auto &deu : helpers[kDeuteron]) {
    KFParticle oneCandidate;
    oneCandidate.Q() = deu.particle.GetQ();
    oneCandidate.AddDaughter(deu.particle);

    for (const auto &p : helpers[kProton]) {
      if (deu.track == p.track || p.particle.GetQ() * deu.particle.GetQ() < 0)
        continue;
      KFParticle twoCandidate{oneCandidate};
      twoCandidate.AddDaughter(p.particle);

      recHyp.chi2_deuprot = twoCandidate.GetChi2() / twoCandidate.GetNDF();
      if (recHyp.chi2_deuprot > fMaxKFchi2[0])
        continue;

      for (const auto &pi : helpers[kPion]) {
        if (p.track == pi.track || deu.track == pi.track || pi.particle.GetQ() * p.particle.GetQ() > 0) 
          continue;
        KFParticle hyperTriton{twoCandidate};
        hyperTriton.AddDaughter(pi.particle);
        recHyp.chi2_3prongs = hyperTriton.GetChi2() / hyperTriton.GetNDF();
        if (recHyp.chi2_3prongs > fMaxKFchi2[1])
          continue;
        ROOT::Math::XYZVectorF decayVtx{hyperTriton.X() - prodVertex.X(), hyperTriton.Y() - prodVertex.Y(), hyperTriton.Z() - prodVertex.Z()};
        ROOT::Math::XYZVectorF mom{hyperTriton.Px(), hyperTriton.Py(), hyperTriton.Pz()};
        recHyp.cosPA = mom.Dot(decayVtx) / std::sqrt(decayVtx.Mag2() * mom.Mag2());
        hyperTriton.SetProductionVertex(prodVertex);
        recHyp.chi2_topology = hyperTriton.GetChi2() / hyperTriton.GetNDF();
        if (recHyp.chi2_topology > fMaxKFchi2[2])
          continue;
        recHyp.l = hyperTriton.GetDecayLength();
        recHyp.r = hyperTriton.GetDecayLengthXY();
        recHyp.px = hyperTriton.GetPx();
        recHyp.py = hyperTriton.GetPy();
        recHyp.pz = hyperTriton.GetPz();
        recHyp.m = hyperTriton.GetMass();
        bool record{!fMC || !fOnlyTrueCandidates};
        if (fMC) {
          int momId = IsTrueHyperTriton3Candidate(deu.track, p.track, pi.track, mcEvent);
          if (momId >= 0) {
            fGenHyp[mcMap[momId]].reco = fRecHyp.size();
            record = true;
          }
        }
        if (record)
          fRecHyp.emplace_back(recHyp);
      }
    }
  }

  /// if event mixing is enabled fill the event mixing pool with deuterons
  // if (fEnableEventMixing && fApplyML) {
  //   FillEventMixingPool(fREvent.fCent, fREvent.fZ, fDeuVector);
  // }

  fTreeHyp3->Fill();

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
}

void AliAnalysisTaskHyperTriton3KF::Terminate(Option_t *) {}

int AliAnalysisTaskHyperTriton3KF::FindEventMixingCentBin(const float centrality) {
  if (centrality > 90) return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskHyperTriton3KF::FindEventMixingZBin(const float zvtx) {
  if (zvtx > 10. || zvtx < -10.) return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

void AliAnalysisTaskHyperTriton3KF::FillEventMixingPool(const float centrality, const float zvtx,
                                                        std::vector<AliESDtrack *> tracks) {
  int centBin = FindEventMixingCentBin(centrality);
  int zBin    = FindEventMixingZBin(zvtx);

  auto &trackVector = fEventMixingPool[centBin][zBin];

  for (auto &t : tracks) {
    trackVector.emplace_back(AliESDtrack{*t});
  }

  if (trackVector.size() - fEventMixingPoolDepth > 0) trackVector.pop_front();

  return;
}

std::vector<AliESDtrack *> AliAnalysisTaskHyperTriton3KF::GetEventMixingTracks(const float centrality,
                                                                               const float zvtx) {
  int centBin = FindEventMixingCentBin(centrality);
  int zBin    = FindEventMixingZBin(zvtx);

  std::vector<AliESDtrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin]) {
    tmpVector.emplace_back(&v);
  }

  return tmpVector;
}

AliAnalysisTaskHyperTriton3KF *AliAnalysisTaskHyperTriton3KF::AddTask(bool isMC, TString suffix) {
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
    ::Error("AddTaskHyperTriton3KF", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHyperTriton3KF";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHyperTriton3KF *task = new AliAnalysisTaskHyperTriton3KF(isMC, tskname.Data());

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
