#include "AliAnalysisTaskHypertritonO2.h"

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
#include "Math/LorentzVector.h"

#include "AliDataFile.h"
#include <TFile.h>
#include <TSpline.h>

#include "Track.h"

ClassImp(AliAnalysisTaskHypertritonO2);

namespace {

using lVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;


struct HelperParticle {
  o2::track::TrackParCov* track = nullptr;
  float        nSigmaTPC = -1.f;
  float        nSigmaTOF = -1.f;
};

constexpr float kHyperTritonMass{2.99131};

constexpr float kDeuMass{1.87561};
constexpr float kPMass{0.938272};
constexpr float kPiMass{0.13957};
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
  return hasTOFout && hasTOFtime;
}

/// helper functions
template <typename T> double Sq(T a) { return a * a; }
template <typename F> double Hypot(F a, F b, F c) { return std::sqrt(Sq(a) + Sq(b) + Sq(c)); }
template <typename F> double Hypot(F a, F b, F c, F d) { return std::sqrt(Sq(a) + Sq(b) + Sq(c) + Sq(d)); }

}    // namespace

AliAnalysisTaskHypertritonO2::AliAnalysisTaskHypertritonO2(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fMC{mc}, fGenHyp{}, fRecHyp{} {
  fTrackCuts.SetMinNClustersTPC(0);
  fTrackCuts.SetEtaRange(-0.9,0.9);
  /// Settings for the custom vertexer

  /// Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());    // Basic Histograms
  DefineOutput(2, TTree::Class());    // Hypertriton Candidates Tree output
}

AliAnalysisTaskHypertritonO2::~AliAnalysisTaskHypertritonO2() {
  if (fListHist) {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3) {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }

  if (fCosPAsplineFile)
    delete fCosPAsplineFile;

}

void AliAnalysisTaskHypertritonO2::UserCreateOutputObjects() {
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
  fTreeHyp3->Branch("RCentrality", &fCent);
  fTreeHyp3->Branch("RTrigger", &fTrigger);
  fTreeHyp3->Branch("RHyperTriton", &fRecHyp);

  if (man->GetMCtruthEventHandler()) {
    fTreeHyp3->Branch("SHyperTriton", &fGenHyp);
    fTreeHyp3->Branch("SGenRecMap", &fGenRecMap);
    fTreeHyp3->Branch("SGenRecDeutMom", &fGenRecDeutMom);
    fTreeHyp3->Branch("SGenRecProtMom", &fGenRecProtMom);
    fTreeHyp3->Branch("SGenRecPiMom", &fGenRecPiMom);
  }

  fCosPAsplineFile = TFile::Open(AliDataFile::GetFileName(fCosPAsplineName).data());
  if (fCosPAsplineFile) {
    fCosPAspline = (TSpline3*)fCosPAsplineFile->Get("cutSpline");
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();

}    /// end UserCreateOutputObjects

void AliAnalysisTaskHypertritonO2::UserExec(Option_t *) {
    // set Magnetic field for KF
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent) {
    ::Fatal("AliAnalysisTaskHypertritonO2::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC) {
    ::Fatal("AliAnalysisTaskHypertritonO2::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent)) {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  if (!fMC && fDownscaling) {
    if (gRandom->Rndm() > fDownscalingFactorByEvent) return;
  }

  fGenHyp.clear();
  fRecHyp.clear();
  fGenRecMap.clear();
  fGenRecDeutMom.clear();
  fGenRecProtMom.clear();
  fGenRecPiMom.clear();

  double pvPos[3], pvCov[6];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pvPos);
  fEventCuts.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);
  fCent = fEventCuts.GetCentrality();

  fTrigger = 0u;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) fTrigger |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) fTrigger |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) fTrigger |= kSemiCentral;
  fTrigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::unordered_map<int, int> mcMap;
  if (fMC) {
    double mcVtx[3];
    mcEvent->GetPrimaryVertex()->GetXYZ(mcVtx);
    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
      AliVParticle *part = mcEvent->GetTrack(iTrack);

      if (!part) {
        ::Warning("AliAnalysisTaskHypertritonO2::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }

      if (std::abs(part->Y()) > 1.) continue;
      if (!IsHyperTriton3(part, mcEvent)) continue;

      double decayVtx[4]{0.0, 0.0, 0.0, 0.0};

      for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD) {
        AliVParticle *daughter = mcEvent->GetTrack(iD);

        if (mcEvent->IsSecondaryFromWeakDecay(iD) && daughter && std::abs(daughter->PdgCode()) != 11) {
          decayVtx[0] = daughter->Xv();
          decayVtx[1] = daughter->Yv();
          decayVtx[2] = daughter->Zv();
          decayVtx[3] = daughter->Tv() - part->Tv();
          break;
        }
      }
      SHyperTritonO2 genHyp;
      genHyp.ct = Hypot(mcVtx[0] - decayVtx[0], mcVtx[1] - decayVtx[1], mcVtx[2] - decayVtx[2]) * kHyperTritonMass / part->P();
      genHyp.pt = part->Pt();
      genHyp.phi = std::atan2(part->Py(),part->Px());
      genHyp.pz = part->Pz();
      genHyp.t = decayVtx[3];
      genHyp.positive = part->PdgCode() > 0;
      mcMap[iTrack] = fGenHyp.size();
      fGenHyp.emplace_back(genHyp);
    }
  }

  std::vector<HelperParticle> helpers[3][2];

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

    if (!fVertexer.getUseAbsDCA()) {
      float cyy = track->GetSigmaY2(), czz = track->GetSigmaZ2(), cyz = track->GetSigmaZY();
      float detYZ = cyy * czz - cyz * cyz;
      if (detYZ < 0.)
        continue;
    }

    for (int iT{0}; iT < 3; ++iT) {
      nSigmasTPC[iT] = fPIDResponse->NumberOfSigmasTPC(track, kAliPID[iT]);
      nSigmasTOF[iT] = fPIDResponse->NumberOfSigmasTOF(track, kAliPID[iT]);
      bool requireTOFpid = track->P() > fRequireTOFpid[iT];
      if (std::abs(nSigmasTPC[iT]) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] && 
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        candidate[iT] = (std::abs(nSigmasTOF[iT]) < fTOFsigmas[iT]) || (!hasTOF && !requireTOFpid);
    }
  
    if (candidate[0] || candidate[1] || candidate[2]) {
      HelperParticle helper;
      helper.track = static_cast<o2::track::TrackParCov*>((AliExternalTrackParam*)track);
      for (int iT{0}; iT < 3; ++iT) {
        if (candidate[iT]) {
          int chargeIndex = (fSwapSign && iT == 0) ? track->GetSigned1Pt() < 0 : track->GetSigned1Pt() > 0;
          helper.nSigmaTPC = nSigmasTPC[iT];
          helper.nSigmaTOF = nSigmasTOF[iT];
          helpers[iT][chargeIndex].push_back(helper);
        }
      }
    }
  }


  /// if event mixing is enabled takes deuteron from the event mixing pool
  // if (fEnableEventMixing && fApplyML) {
  //   deuterons = GetEventMixingTracks(fREvent.fCent, fREvent.fZ);
  // }

  RHyperTritonO2 recHyp;
  int indices[2][3]{{1,1,0},{0,0,1}};
  for (int idx{0}; idx < 2; ++idx) {
    for (const auto &deu : helpers[kDeuteron][indices[idx][0]]) {
      for (const auto &p : helpers[kProton][indices[idx][1]]) {
        if (deu.track == p.track)
          continue;
        for (const auto &pi : helpers[kPion][indices[idx][2]]) {
          if (p.track == pi.track || deu.track == pi.track) 
            continue;
          
          int nVert{0};
          try {
            nVert = fVertexer.process(*deu.track, *p.track, *pi.track);
          } catch (std::runtime_error& e) {  }
          if (nVert) {
            auto vert = fVertexer.getPCACandidate();

            auto& deuTrack = fVertexer.getTrack(0);
            auto& prTrack = fVertexer.getTrack(1);
            auto& piTrack = fVertexer.getTrack(2);

            lVector ldeu{deuTrack.Pt(), deuTrack.Eta(), deuTrack.Phi(), kDeuMass};
            lVector lpro{prTrack.Pt(), prTrack.Eta(), prTrack.Phi(), kPMass};
            lVector lpi{piTrack.Pt(), piTrack.Eta(), piTrack.Phi(), kPiMass};
            lVector hypertriton{ldeu + lpro + lpi};

            const float mass = hypertriton.mass();
            if (mass < fMassWindow[0] || mass > fMassWindow[1])
              continue;
            ROOT::Math::XYZVectorF decayVtx{(float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2])};

            const float totalMom = hypertriton.P();
            const float len = std::sqrt(decayVtx.Mag2());
            recHyp.cosPA = hypertriton.Vect().Dot(decayVtx) / (totalMom * len);
            const float cosPA = fUseAbsCosPAcut ? std::abs(recHyp.cosPA) : recHyp.cosPA;
            recHyp.ct = len * kHyperTritonMass / totalMom; 
            if (recHyp.ct < fCandidateCtRange[0] || recHyp.ct > fCandidateCtRange[1])
              continue;
            if (fCosPAspline) {
              if (cosPA < fCosPAspline->Eval(recHyp.ct))
                continue;
            } else if (cosPA < fMinCosPA) {
              continue;
            }
            recHyp.candidates = nVert;
            recHyp.r = decayVtx.Rho();
            float hSign = deu.track->Charge() > 0 ? 1. : -1;
            recHyp.pt = hSign * hypertriton.pt();
            recHyp.phi = hypertriton.phi();
            recHyp.pz = hypertriton.pz();
            recHyp.m = mass;

            float dca[2], bCov[3];
            deu.track->GetImpactParameters(dca, bCov);
            recHyp.dca_de = std::hypot(dca[0],dca[1]);
            p.track->GetImpactParameters(dca, bCov);
            recHyp.dca_pr = std::hypot(dca[0],dca[1]);
            pi.track->GetImpactParameters(dca, bCov);
            recHyp.dca_pi = std::hypot(dca[0],dca[1]);

            recHyp.hasTOF_de = HasTOF(deu.track);
            recHyp.hasTOF_pr = HasTOF(p.track);
            recHyp.hasTOF_pi = HasTOF(pi.track);

            recHyp.tofNsig_de = deu.nSigmaTOF;
            recHyp.tofNsig_pr = p.nSigmaTOF;
            recHyp.tofNsig_pi = pi.nSigmaTOF;

            recHyp.tpcNsig_de = deu.nSigmaTPC;
            recHyp.tpcNsig_pr = p.nSigmaTPC;
            recHyp.tpcNsig_pi = pi.nSigmaTPC;

            auto& deuPos = fVertexer.getTrackPos(0); 
            auto& proPos = fVertexer.getTrackPos(1); 
            auto& piPos = fVertexer.getTrackPos(2); 

            recHyp.dca_de_pr = Hypot(deuPos[0] - proPos[0], deuPos[1] - proPos[1], deuPos[2] - proPos[2]);
            recHyp.dca_de_pi = Hypot(deuPos[0] - piPos[0], deuPos[1] - piPos[1], deuPos[2] - piPos[2]);
            recHyp.dca_pr_pi = Hypot(proPos[0] - piPos[0], proPos[1] - piPos[1], proPos[2] - piPos[2]);

            recHyp.dca_de_sv = Hypot(deuPos[0] - vert[0], deuPos[1] - vert[1], deuPos[2] - vert[2]);
            recHyp.dca_pr_sv = Hypot(proPos[0] - vert[0], proPos[1] - vert[1], proPos[2] - vert[2]);
            recHyp.dca_pi_sv = Hypot(piPos[0] - vert[0], piPos[1] - vert[1], piPos[2] - vert[2]);

            recHyp.chi2 = fVertexer.getChi2AtPCACandidate();

            recHyp.tpcClus_de = deu.track->GetTPCsignalN();
            recHyp.tpcClus_pr = p.track->GetTPCsignalN();
            recHyp.tpcClus_pi = pi.track->GetTPCsignalN();

            bool record{!fMC || !fOnlyTrueCandidates};
            if (fMC) {
              int momId = IsTrueHyperTriton3Candidate((AliESDtrack*)deu.track, (AliESDtrack*)p.track, (AliESDtrack*)pi.track, mcEvent);
              record = record || momId >=0;
              if (record) {
                fGenRecMap.push_back(mcMap[momId]);
                fGenRecDeutMom.push_back(deu.track->P());
                fGenRecProtMom.push_back(p.track->P());
                fGenRecPiMom.push_back(pi.track->P());
              }
            }
            if (record)
              fRecHyp.emplace_back(recHyp);
          }
        }
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

void AliAnalysisTaskHypertritonO2::Terminate(Option_t *) {}

int AliAnalysisTaskHypertritonO2::FindEventMixingCentBin(const float centrality) {
  if (centrality > 90) return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskHypertritonO2::FindEventMixingZBin(const float zvtx) {
  if (zvtx > 10. || zvtx < -10.) return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

void AliAnalysisTaskHypertritonO2::FillEventMixingPool(const float centrality, const float zvtx,
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

std::vector<AliESDtrack *> AliAnalysisTaskHypertritonO2::GetEventMixingTracks(const float centrality,
                                                                               const float zvtx) {
  int centBin = FindEventMixingCentBin(centrality);
  int zBin    = FindEventMixingZBin(zvtx);

  std::vector<AliESDtrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin]) {
    tmpVector.emplace_back(&v);
  }

  return tmpVector;
}

AliAnalysisTaskHypertritonO2 *AliAnalysisTaskHypertritonO2::AddTask(bool isMC, TString suffix) {
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
    ::Error("AddTaskHypertritonO2", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHypertritonO2";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHypertritonO2 *task = new AliAnalysisTaskHypertritonO2(isMC, tskname.Data());

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
