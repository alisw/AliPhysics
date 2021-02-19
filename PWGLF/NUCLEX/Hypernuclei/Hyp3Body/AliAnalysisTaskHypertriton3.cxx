#include "AliAnalysisTaskHypertriton3.h"

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

#include "Math/GenVector/Boost.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"

#include "AliDataFile.h"
#include <TFile.h>
#include <TSpline.h>

#include "Track.h"
#include <memory>

#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

ClassImp(AliAnalysisTaskHypertriton3);

namespace
{

  using lVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;

  struct HelperParticle
  {
    o2::track::TrackParCov *track = nullptr;
    int index = -1;
    float nSigmaTPC = -1.f;
    float nSigmaTOF = -1.f;
    KFParticle particle = KFParticle();
  };

  constexpr float kDeuMass{1.87561};
  constexpr float kPMass{0.938272};
  constexpr float kPiMass{0.13957};
  constexpr float kMasses[3]{kDeuMass, kPMass, kPiMass};
  constexpr AliPID::EParticleType kAliPID[3]{AliPID::kDeuteron, AliPID::kProton, AliPID::kPion};
  const int kPDGs[3]{AliPID::ParticleCode(kAliPID[0]), AliPID::ParticleCode(kAliPID[1]), AliPID::ParticleCode(kAliPID[2])};

  bool IsHyperTriton3(const AliVParticle *vPart, AliMCEvent *mcEvent)
  {
    int nDaughters = 0;

    int vPartPDG = vPart->PdgCode();
    int vPartLabel = vPart->GetLabel();

    if (!mcEvent->IsPhysicalPrimary(vPartLabel) || (std::abs(vPartPDG) != 1010010030))
      return false;

    for (int iD = vPart->GetDaughterFirst(); iD <= vPart->GetDaughterLast(); iD++)
    {
      AliVParticle *dPart = mcEvent->GetTrack(iD);

      int dPartPDG = dPart->PdgCode();
      if (std::abs(dPartPDG) != 11)
        nDaughters++;
    }
    if (nDaughters == 3)
      return true;
    return false;
  }

  int IsTrueHyperTriton3Candidate(AliESDtrack *t1, AliESDtrack *t2, AliESDtrack *t3, AliMCEvent *mcEvent)
  {
    if (!mcEvent)
      return 0;

    int lab1 = std::abs(t1->GetLabel());
    int lab2 = std::abs(t2->GetLabel());
    int lab3 = std::abs(t3->GetLabel());

    if (mcEvent->IsPhysicalPrimary(lab1))
      return -1;
    if (mcEvent->IsPhysicalPrimary(lab2))
      return -1;
    if (mcEvent->IsPhysicalPrimary(lab3))
      return -1;

    AliVParticle *part1 = mcEvent->GetTrack(lab1);
    AliVParticle *part2 = mcEvent->GetTrack(lab2);
    AliVParticle *part3 = mcEvent->GetTrack(lab3);

    if (!part1 || !part2 || !part3)
      return -1;

    int mom1 = part1->GetMother();
    int mom2 = part2->GetMother();
    int mom3 = part3->GetMother();

    if (mom1 != mom2 || mom1 != mom3 || mom2 != mom3)
      return -1;

    AliVParticle *mom = mcEvent->GetTrack(mom1);
    if (!mom)
      return -1;

    return (IsHyperTriton3(mom, mcEvent)) ? mom1 : -1;
  }

  bool HasTOF(AliVTrack *track)
  {
    const bool hasTOFout = track->GetStatus() & AliVTrack::kTOFout;
    const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
    return hasTOFout && hasTOFtime;
  }

} // namespace

AliAnalysisTaskHypertriton3::AliAnalysisTaskHypertriton3(bool mc, std::string name)
    : AliAnalysisTaskSE(name.data()), fEventCuts{}, fVertexer{}, fVertexerLambda{}, fMC{mc}
{
  fTrackCuts.SetMinNClustersTPC(0);
  fTrackCuts.SetEtaRange(-0.9, 0.9);
  /// Settings for the custom vertexer

  /// Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); // Basic Histograms
  DefineOutput(2, TTree::Class()); // Hypertriton Candidates Tree output
}

AliAnalysisTaskHypertriton3::~AliAnalysisTaskHypertriton3()
{
  if (fListHist)
  {
    delete fListHist;
    fListHist = nullptr;
  }

  if (fTreeHyp3)
  {
    delete fTreeHyp3;
    fTreeHyp3 = nullptr;
  }

  if (fCosPAsplineFile)
    delete fCosPAsplineFile;

  if (fGenHypKF || fGenHypO2)
  {
    if (fGenHypKF)
      delete fGenHypKF;
    if (fGenHypO2)
      delete fGenHypO2;
  }
  else if (fRecHyp)
    delete fRecHyp;
}

void AliAnalysisTaskHypertriton3::UserCreateOutputObjects()
{
  fCounter = 0;

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  fInputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse = fInputHandler->GetPIDResponse();

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

  fHistDecVertexRes =
      new TH1D("fHistDecVertexRes", "; Resoultion(cm); Counts", 40, -1, 1);
  fListHist->Add(fHistNSigmaDeu);
  fListHist->Add(fHistNSigmaP);
  fListHist->Add(fHistNSigmaPi);
  fListHist->Add(fHistInvMass);
  fListHist->Add(fHistDecVertexRes);

  OpenFile(2);
  fTreeHyp3 = new TTree("Hyp3O2", "Hypetriton 3 Body with the O2 Vertexer");

  if (fMC && man->GetMCtruthEventHandler())
  {
    if (fKF)
    {
      fGenHypKF = new SHyperTriton3KF;
      fRecHyp = (RHyperTriton *)fGenHypKF;
      fTreeHyp3->Branch("SHyperTriton", fGenHypKF);
    }
    else
    {
      fGenHypO2 = new SHyperTriton3O2;
      fRecHyp = (RHyperTriton *)fGenHypO2;
      fTreeHyp3->Branch("SHyperTriton", fGenHypO2);
    }
  }
  else
  {
    if (fKF)
    {
      fRecHyp = new RHyperTriton3KF;
      fTreeHyp3->Branch("RHyperTriton", static_cast<RHyperTriton3KF *>(fRecHyp));
    }
    else
    {
      fRecHyp = new RHyperTriton3O2;
      fTreeHyp3->Branch("RHyperTriton", static_cast<RHyperTriton3O2 *>(fRecHyp));
    }
  }
  fCosPAsplineFile = TFile::Open(AliDataFile::GetFileName(fCosPAsplineName).data());
  if (fCosPAsplineFile)
  {
    fCosPAspline = (TSpline3 *)fCosPAsplineFile->Get("cutSpline");
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);

  AliPDG::AddParticlesToPdgDataBase();

} /// end UserCreateOutputObjects

void AliAnalysisTaskHypertriton3::UserExec(Option_t *)
{

  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent)
  {
    ::Fatal("AliAnalysisTaskHypertriton3::UserExec", "AliESDEvent not found.");
    return;
  }

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent && fMC)
  {
    ::Fatal("AliAnalysisTaskHypertriton3::UserExec", "Could not retrieve MC event");
    return;
  }

  if (!fEventCuts.AcceptEvent(esdEvent))
  {
    PostData(1, fListHist);
    PostData(2, fTreeHyp3);
    return;
  }

  if (!fMC && fDownscaling)
  {
    if (gRandom->Rndm() > fDownscalingFactorByEvent)
      return;
  }

  double pvPos[3], pvCov[6];
  fEventCuts.GetPrimaryVertex()->GetXYZ(pvPos);
  fEventCuts.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);
  fRecHyp->centrality = fEventCuts.GetCentrality();

  fRecHyp->trigger = 0u;
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
    fRecHyp->trigger |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
    fRecHyp->trigger |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    fRecHyp->trigger |= kSemiCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
    fRecHyp->trigger |= kHighMultV0;
  fRecHyp->trigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::vector<HelperParticle> helpers[3][2];
  std::vector<AliESDtrack *> deuPiTracks[2][2];
  std::vector<AliESDtrack *> prPiTracks[2][2];
  std::vector<EventMixingTrack> deuteronsForMixing;

  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
  {
    AliESDtrack *track = esdEvent->GetTrack(iTrack);

    if (!track)
      continue;

    if (!fTrackCuts.AcceptTrack(track))
      continue;

    if (fMC && fOnlyTrueCandidates)
    {
      int lab = std::abs(track->GetLabel());
      if (!mcEvent->IsSecondaryFromWeakDecay(lab))
        continue;
      AliVParticle *part = mcEvent->GetTrack(lab);
      AliVParticle *moth = mcEvent->GetTrack(part->GetMother());
      if (std::abs(moth->PdgCode()) != 1010010030)
        continue;
    }

    bool candidate[3]{false, false, false};
    float nSigmasTPC[3]{-1., -1., -1.}, nSigmasTOF[3]{-1., -1., -1.};
    bool hasTOF{HasTOF(track)};
    float dca[2];
    track->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    if (fUseCovarianceCut)
    {
      float cyy = track->GetSigmaY2(), czz = track->GetSigmaZ2(), cyz = track->GetSigmaZY();
      float detYZ = cyy * czz - cyz * cyz;
      if (detYZ < 0.)
        continue;
    }

    for (int iT{0}; iT < 3; ++iT)
    {
      nSigmasTPC[iT] = fPIDResponse->NumberOfSigmasTPC(track, kAliPID[iT]);
      nSigmasTOF[iT] = fPIDResponse->NumberOfSigmasTOF(track, kAliPID[iT]);
      bool requireTOFpid = track->P() > fRequireTOFpid[iT];
      if (std::abs(nSigmasTPC[iT]) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] &&
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        candidate[iT] = (std::abs(nSigmasTOF[iT]) < fTOFsigmas[iT]) || (!hasTOF && !requireTOFpid);
    }

    if (candidate[0] || candidate[1] || candidate[2])
    {
      HelperParticle helper;
      helper.track = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)track);
      for (int iT{0}; iT < 3; ++iT)
      {
        if (candidate[iT])
        {
          int chargeIndex = (fSwapSign && iT == fMixingTrack) ? track->GetSigned1Pt() < 0 : track->GetSigned1Pt() > 0;
          helper.nSigmaTPC = nSigmasTPC[iT];
          helper.nSigmaTOF = nSigmasTOF[iT];
          if (fKF)
          {
            double posmom[6], cov[21];
            track->GetXYZ(posmom);
            track->GetPxPyPz(posmom + 3);
            track->GetCovarianceXYZPxPyPz(cov);
            helper.particle.Create(posmom, cov, track->Charge(), kMasses[iT]);
            helper.particle.Chi2() = track->GetTPCchi2();
            helper.particle.NDF() = track->GetNumberOfTPCClusters() * 2;
          }
          if (iT == fMixingTrack && fEnableEventMixing)
            deuteronsForMixing.emplace_back(track, nSigmasTPC[iT], nSigmasTOF[iT], 0);
          else
          {
            helpers[iT][chargeIndex].push_back(helper);
            if (fUseDoubleV0s)
            {

              if (iT == 0)
                deuPiTracks[chargeIndex][0].push_back(track);
              if (iT == 1)
                prPiTracks[chargeIndex][0].push_back(track);
              if (iT == 2)
              {
                deuPiTracks[chargeIndex][1].push_back(track);
                prPiTracks[chargeIndex][1].push_back(track);
              }
            }
          }
        }
      }
    }
  }

  if (fEnableEventMixing)
  {
    auto mixingDeuterons = GetEventMixingTracks(fEventCuts.GetCentrality(), pvPos[2]);
    for (auto mixTrack : mixingDeuterons)
    {
      HelperParticle helper;
      AliESDtrack *track = &(mixTrack->track);
      helper.track = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)track);
      int chargeIndex = track->GetSigned1Pt() > 0;
      helper.nSigmaTPC = mixTrack->nSigmaTPC;
      helper.nSigmaTOF = mixTrack->nSigmaTOF;
      if (fKF)
      {
        double posmom[6], cov[21];
        track->GetXYZ(posmom);
        track->GetPxPyPz(posmom + 3);
        track->GetCovarianceXYZPxPyPz(cov);
        helper.particle.Create(posmom, cov, track->Charge(), kMasses[0]);
        helper.particle.Chi2() = track->GetTPCchi2();
        helper.particle.NDF() = track->GetNumberOfTPCClusters() * 2;
      }
      helpers[fMixingTrack][chargeIndex].push_back(helper);
      //      fakeEvent.AddTrack(track);
    }
  }

  lVector hypertriton;
  ROOT::Math::XYZVectorF decayVtx, decayVtxLambda;
  lVector lproL, lpiL;
  std::unordered_map<int, int> mcMap;
  auto fillTreeInfo = [&](std::array<AliESDtrack *, 3> tracks, std::array<float, 3> nSigmaTPC, std::array<float, 3> nSigmaTOF) {
    const float mass = hypertriton.mass();
    if (mass < fMassWindow[0] || mass > fMassWindow[1])
      return false;

    const float totalMom = hypertriton.P();
    const float len = std::sqrt(decayVtx.Mag2());
    fRecHyp->cosPA = hypertriton.Vect().Dot(decayVtx) / (totalMom * len);
    const float cosPA = fUseAbsCosPAcut ? std::abs(fRecHyp->cosPA) : fRecHyp->cosPA;
    fRecHyp->ct = len * kHyperTritonMass / totalMom;
    if (fRecHyp->ct < fCandidateCtRange[0] || fRecHyp->ct > fCandidateCtRange[1])
      return false;
    if (fCosPAspline)
    {
      if (cosPA < fCosPAspline->Eval(fRecHyp->ct))
        return false;
    }
    else if (cosPA < fMinCosPA)
    {
      return false;
    }
    fRecHyp->r = decayVtx.Rho();
    fRecHyp->positive = tracks[0]->Charge() > 0;
    fRecHyp->pt = hypertriton.pt();
    fRecHyp->phi = hypertriton.phi();
    fRecHyp->pz = hypertriton.pz();
    fRecHyp->m = mass;

    float dca[2], bCov[3];
    tracks[0]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_de = std::hypot(dca[0], dca[1]);
    tracks[1]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_pr = std::hypot(dca[0], dca[1]);
    tracks[2]->GetImpactParameters(dca, bCov);
    fRecHyp->dca_pi = std::hypot(dca[0], dca[1]);

    fRecHyp->hasTOF_de = HasTOF(tracks[0]);
    fRecHyp->hasTOF_pr = HasTOF(tracks[1]);
    fRecHyp->hasTOF_pi = HasTOF(tracks[2]);

    fRecHyp->tofNsig_de = nSigmaTOF[0];
    fRecHyp->tofNsig_pr = nSigmaTOF[1];
    fRecHyp->tofNsig_pi = nSigmaTOF[2];

    fRecHyp->tpcNsig_de = nSigmaTPC[0];
    fRecHyp->tpcNsig_pr = nSigmaTPC[1];
    fRecHyp->tpcNsig_pi = nSigmaTPC[2];

    fRecHyp->tpcClus_de = tracks[0]->GetTPCsignalN();
    fRecHyp->tpcClus_pr = tracks[1]->GetTPCsignalN();
    fRecHyp->tpcClus_pi = tracks[2]->GetTPCsignalN();

    fRecHyp->its_clusmap_de = tracks[0]->GetITSClusterMap();
    fRecHyp->its_clusmap_pr = tracks[1]->GetITSClusterMap();
    fRecHyp->its_clusmap_pi = tracks[2]->GetITSClusterMap();

    fRecHyp->is_ITSrefit_de = tracks[0]->GetStatus() & AliVTrack::kITSrefit;
    fRecHyp->is_ITSrefit_pr = tracks[1]->GetStatus() & AliVTrack::kITSrefit;
    fRecHyp->is_ITSrefit_pi = tracks[2]->GetStatus() & AliVTrack::kITSrefit;

    if (fLambdaCheck)
    {
      int nVertLambda{!fUseDoubleV0s};
      if (!fUseDoubleV0s)
      {
        try
        {
          o2::track::TrackParCov *prTrack = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)tracks[1]);
          o2::track::TrackParCov *piTrack = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)tracks[2]);
          nVertLambda = fVertexerLambda.process(*prTrack, *piTrack);
        }
        catch (std::runtime_error &e)
        {
        }
        if (nVertLambda)
        {
          auto vertLambda = fVertexerLambda.getPCACandidate();
          fVertexerLambda.propagateTracksToVertex();
          auto &prTrackL = fVertexerLambda.getTrack(0);
          auto &piTrackL = fVertexerLambda.getTrack(1);
          decayVtxLambda.SetCoordinates((float)(vertLambda[0] - pvPos[0]), (float)(vertLambda[1] - pvPos[1]), (float)(vertLambda[2] - pvPos[2]));
          lproL.SetCoordinates((float)prTrackL.Pt(), (float)prTrackL.Eta(), (float)prTrackL.Phi(), kPMass);
          lpiL.SetCoordinates((float)piTrackL.Pt(), (float)piTrackL.Eta(), (float)piTrackL.Phi(), kPiMass);
        }
      }

      if (nVertLambda)
      {
        auto vertLambda = fVertexerLambda.getPCACandidate();
        fVertexerLambda.propagateTracksToVertex();
        auto &prTrackL = fVertexerLambda.getTrack(0);
        auto &piTrackL = fVertexerLambda.getTrack(1);
        ROOT::Math::XYZVectorF decayVtxLambda{(float)(vertLambda[0] - pvPos[0]), (float)(vertLambda[1] - pvPos[1]), (float)(vertLambda[2] - pvPos[2])};
        lVector lproL{(float)prTrackL.Pt(), (float)prTrackL.Eta(), (float)prTrackL.Phi(), kPMass};
        lVector lpiL{(float)piTrackL.Pt(), (float)piTrackL.Eta(), (float)piTrackL.Phi(), kPiMass};
        lVector lambda{lproL + lpiL};
        fRecHyp->mppi_vert = lambda.mass();
        const float lambdaLen = std::sqrt(decayVtxLambda.Mag2());
        fRecHyp->cosPA_Lambda = lambda.Vect().Dot(decayVtxLambda) / (lambda.P() * lambdaLen);
        fRecHyp->dca_lambda_hyper = Hypot(decayVtx.x() - vertLambda[0], decayVtx.y() - vertLambda[1], decayVtx.z() - vertLambda[2]);
      }
    }
    return true;
  };

  if (fUseDoubleV0s)
  {

    RHyperTriton3O2 &doubleV0sRecHyp = *(RHyperTriton3O2 *)fRecHyp;

    auto acceptTracks = [&](std::array<AliESDtrack *, 3> tracks) -> bool {
      int candidate = 0;
      for (int iT{0}; iT < 3; ++iT)
      {
        float dca[2];
        tracks[iT]->GetImpactParameters(dca[0], dca[1]);
        double dcaNorm = std::hypot(dca[0], dca[1]);
        float nSigmasTPC = fPIDResponse->NumberOfSigmasTPC(tracks[iT], kAliPID[iT]);
        float nSigmasTOF = fPIDResponse->NumberOfSigmasTOF(tracks[iT], kAliPID[iT]);
        bool requireTOFpid = tracks[iT]->P() > fRequireTOFpid[iT];
        if (std::abs(nSigmasTPC) < fTPCsigmas[iT] && dcaNorm > fMinTrackDCA[iT] && tracks[iT]->Pt() < fTrackPtRange[iT][1] &&
            tracks[iT]->Pt() > fTrackPtRange[iT][0] && tracks[iT]->GetTPCsignalN() >= fMinTPCpidClusters[iT])
        {
          candidate += (std::abs(nSigmasTOF) < fTOFsigmas[iT]) || (!HasTOF(tracks[iT]) && !requireTOFpid);
        }
      }
      int accept = candidate < 3 ? 0 : 1;
      return accept;
    };

    auto deuPiV0s = fV0Vertexer.Tracks2V0vertices3Body(deuPiTracks, esdEvent->GetPrimaryVertex(), esdEvent->GetMagneticField(), fPIDResponse);
    auto prPiV0s = fV0Vertexer.Tracks2V0vertices3Body(prPiTracks, esdEvent->GetPrimaryVertex(), esdEvent->GetMagneticField(), fPIDResponse);

    for (auto &deuPiV0 : deuPiV0s)
    {

      int isPiPositive = CheckPionCharge(deuPiTracks, deuPiV0);
      if (isPiPositive < 0)
        continue;

      int piIndex = isPiPositive ? deuPiV0.GetPindex() : deuPiV0.GetNindex();
      int deuIndex = isPiPositive ? deuPiV0.GetNindex() : deuPiV0.GetPindex();

      auto piTrack = deuPiTracks[isPiPositive][1][piIndex];
      auto deuTrack = deuPiTracks[1 - isPiPositive][0][deuIndex];

      for (auto &prPiV0 : prPiV0s)
      {

        int isSecondPiPositive = CheckPionCharge(prPiTracks, prPiV0);
        if (isSecondPiPositive < 0)
          continue;

        auto secondPiTrack = isSecondPiPositive ? prPiTracks[1][1][prPiV0.GetPindex()] : prPiTracks[0][1][prPiV0.GetNindex()];
        if (piTrack != secondPiTrack) // require common pion
          continue;

        int prIndex = isPiPositive ? prPiV0.GetNindex() : prPiV0.GetPindex();
        auto prTrack = prPiTracks[1 - isPiPositive][0][prIndex];

        std::array<AliESDtrack *, 3> tracks{deuTrack, prTrack, piTrack};

        if (tracks[0] == tracks[1] || tracks[1] == tracks[2] || tracks[0] == tracks[2])
          continue;

        if (deuPiTracks[1 - isPiPositive][0].size() < deuIndex)
          continue;

        if (prPiTracks[1 - isPiPositive][0].size() < prIndex)
          continue;

        if (!acceptTracks(tracks))
          continue;

        double firstXYZ[3], secondXYZ[3], firstCov[6], secondCov[6];
        deuPiV0.GetXYZ(firstXYZ[0], firstXYZ[1], firstXYZ[2]);
        prPiV0.GetXYZ(secondXYZ[0], secondXYZ[1], secondXYZ[2]);
        deuPiV0.GetVertex().GetCovarianceMatrix(firstCov);
        prPiV0.GetVertex().GetCovarianceMatrix(secondCov);
        double vert[3], sigmaVert[3];
        constexpr int covIndex[3]{0, 2, 5};
        for (int i{0}; i < 3; ++i)
        {
          vert[i] = (firstXYZ[i] / firstCov[covIndex[i]] + secondXYZ[i] / secondCov[covIndex[i]]) / (1. / firstCov[covIndex[i]] + 1. / secondCov[covIndex[i]]);
          sigmaVert[i] = TMath::Sqrt(1. / (1. / (firstCov[covIndex[i]] * firstCov[covIndex[i]]) + 1. / (secondCov[covIndex[i]] * secondCov[covIndex[i]])));
        }
        double sumCov[6];
        for (int i{0}; i < 6; ++i)
        {
          sumCov[i] = firstCov[i] + secondCov[i];
        }

        std::array<float, 3> nSigmasTPC{fPIDResponse->NumberOfSigmasTPC(tracks[0], kAliPID[0]), fPIDResponse->NumberOfSigmasTPC(tracks[1], kAliPID[1]), fPIDResponse->NumberOfSigmasTPC(tracks[2], kAliPID[2])};
        std::array<float, 3> nSigmasTOF{fPIDResponse->NumberOfSigmasTOF(tracks[0], kAliPID[0]), fPIDResponse->NumberOfSigmasTOF(tracks[1], kAliPID[1]), fPIDResponse->NumberOfSigmasTOF(tracks[2], kAliPID[2])};

        lVector ldeu{(float)tracks[0]->Pt(), (float)tracks[0]->Eta(), (float)tracks[0]->Phi(), kDeuMass};
        lVector lpro{(float)tracks[1]->Pt(), (float)tracks[1]->Eta(), (float)tracks[1]->Phi(), kPMass};
        lVector lpi{(float)tracks[2]->Pt(), (float)tracks[2]->Eta(), (float)tracks[2]->Phi(), kPiMass};
        hypertriton = ldeu + lpro + lpi;
        doubleV0sRecHyp.mppi = (lpro + lpi).mass2();
        doubleV0sRecHyp.mdpi = (ldeu + lpi).mass2();
        ROOT::Math::Boost boostHyper{hypertriton.BoostToCM()};
        auto d{boostHyper(ldeu).Vect()};
        auto lambda{boostHyper(lpro + lpi).Vect()};
        auto p{boostHyper(lpro).Vect()};
        auto pi{boostHyper(lpi).Vect()};
        doubleV0sRecHyp.momDstar = std::sqrt(d.Mag2());
        doubleV0sRecHyp.cosThetaStar = d.Dot(hypertriton.Vect()) / (doubleV0sRecHyp.momDstar * hypertriton.P());
        doubleV0sRecHyp.cosTheta_ProtonPiH = p.Dot(pi) / std::sqrt(p.Mag2() * pi.Mag2());
        decayVtx.SetCoordinates((float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2]));
        doubleV0sRecHyp.candidates = 3;

        double deuPos[3], proPos[3], piPos[3];

        // AliESDVertex
        AliESDVertex *decESDVtx = new AliESDVertex(vert, sigmaVert, "vert");

        for (int i = 0; i < 3; i++)
        {
          tracks[i]->PropagateToDCA(decESDVtx, esdEvent->GetMagneticField(), 25.);
        }

        tracks[0]->GetXYZ(deuPos);
        tracks[1]->GetXYZ(proPos);
        tracks[2]->GetXYZ(piPos);

        if (fMC && IsTrueHyperTriton3Candidate(tracks[0], tracks[1], tracks[2], mcEvent) >= 0)
        {
          int lab = std::abs(tracks[0]->GetLabel());
          if (!mcEvent->IsSecondaryFromWeakDecay(lab))
            continue;
          AliVParticle *part = mcEvent->GetTrack(lab);
          fHistDecVertexRes->Fill(vert[0] - part->Xv() + vert[1] - part->Yv() + vert[2] - part->Zv());
        }

        doubleV0sRecHyp.dca_de_pr = Hypot(deuPos[0] - proPos[0], deuPos[1] - proPos[1], deuPos[2] - proPos[2]);
        if (doubleV0sRecHyp.dca_de_pr > fMaxTrack2TrackDCA[0])
          continue;
        doubleV0sRecHyp.dca_de_pi = Hypot(deuPos[0] - piPos[0], deuPos[1] - piPos[1], deuPos[2] - piPos[2]);
        if (doubleV0sRecHyp.dca_de_pi > fMaxTrack2TrackDCA[1])
          continue;

        doubleV0sRecHyp.dca_pr_pi = Hypot(proPos[0] - piPos[0], proPos[1] - piPos[1], proPos[2] - piPos[2]);
        if (doubleV0sRecHyp.dca_pr_pi > fMaxTrack2TrackDCA[2])
          continue;

        doubleV0sRecHyp.dca_pr_sv = Hypot(proPos[0] - vert[0], proPos[1] - vert[1], proPos[2] - vert[2]);
        if (doubleV0sRecHyp.dca_pr_sv > fMaxTrack2SVDCA[1])
          continue;

        doubleV0sRecHyp.dca_pi_sv = Hypot(piPos[0] - vert[0], piPos[1] - vert[1], piPos[2] - vert[2]);
        if (doubleV0sRecHyp.dca_pi_sv > fMaxTrack2SVDCA[2])
          continue;

        doubleV0sRecHyp.dca_de_sv = Hypot(deuPos[0] - vert[0], deuPos[1] - vert[1], deuPos[2] - vert[2]);
        if (doubleV0sRecHyp.dca_de_sv > fMaxTrack2SVDCA[0])
          continue;
        
        doubleV0sRecHyp.cosPA_Lambda = prPiV0.GetV0CosineOfPointingAngle();
        doubleV0sRecHyp.chi2 = prPiV0.GetChi2V0() + deuPiV0.GetChi2V0();

        for (int i = 0; i < 3; i++)
        {
          tracks[i]->PropagateToDCA(esdEvent->GetPrimaryVertex(), esdEvent->GetMagneticField(), 25.);
        }

        if (!fillTreeInfo(tracks, nSigmasTPC, nSigmasTOF))
          continue;

        bool record{!fMC || !fOnlyTrueCandidates};
        if (fMC)
        {
          int momId = IsTrueHyperTriton3Candidate(tracks[0], tracks[1], tracks[2], mcEvent);
          record = record || momId >= 0;
          if (record)
          {
            if (fKF)
              FillGenHypertriton(fGenHypKF, momId, true, mcEvent);
            else
              FillGenHypertriton(fGenHypO2, momId, true, mcEvent);
            mcMap[momId] = 1;
          }
        }
        if (record)
          fTreeHyp3->Fill();
      }
    }
  }
  else
  {

    fVertexer.setBz(esdEvent->GetMagneticField());
    fVertexerLambda.setBz(esdEvent->GetMagneticField());
    int indices[2][3]{{1, 1, 0}, {0, 0, 1}};

    KFPVertex kfPVertex;
    kfPVertex.SetXYZ(pvPos[0], pvPos[1], pvPos[2]);
    kfPVertex.SetCovarianceMatrix(pvCov[0], pvCov[1], pvCov[2], pvCov[3], pvCov[4], pvCov[5]);
    kfPVertex.SetChi2(fEventCuts.GetPrimaryVertex()->GetChi2());
    kfPVertex.SetNDF(fEventCuts.GetPrimaryVertex()->GetNDF());
    kfPVertex.SetNContributors(fEventCuts.GetPrimaryVertex()->GetNContributors());

    KFParticle prodVertex{kfPVertex};

    RHyperTriton3KF &kfRecHyp = *(RHyperTriton3KF *)fRecHyp;
    RHyperTriton3O2 &o2RecHyp = *(RHyperTriton3O2 *)fRecHyp;

    for (int idx{0}; idx < 2; ++idx)
    {
      for (const auto &deu : helpers[kDeuteron][indices[idx][0]])
      {
        KFParticle oneCandidate;
        if (fKF)
        {
          oneCandidate.Q() = deu.particle.GetQ();
          oneCandidate.AddDaughter(deu.particle);
        }
        for (const auto &p : helpers[kProton][indices[idx][1]])
        {
          if (deu.track == p.track)
            continue;

          KFParticle twoCandidate{oneCandidate};
          if (fKF)
          {
            twoCandidate.AddDaughter(p.particle);
            kfRecHyp.chi2_deuprot = twoCandidate.GetChi2() / twoCandidate.GetNDF();
            if (kfRecHyp.chi2_deuprot > fMaxKFchi2[0] || kfRecHyp.chi2_deuprot < 0.)
              continue;
          }
          for (const auto &pi : helpers[kPion][indices[idx][2]])
          {
            if (p.track == pi.track || deu.track == pi.track || deu.track == p.track)
              continue;

            ROOT::Math::SVector<double, 3U> vert;
            if (!fKF)
            {
              int nVert{0};
              try
              {
                nVert = fVertexer.process(*deu.track, *p.track, *pi.track);
              }
              catch (std::runtime_error &e)
              {
              }
              if (!nVert)
                continue;

              fVertexer.propagateTracksToVertex();
              auto &deuTrack = fVertexer.getTrack(0);
              auto &prTrack = fVertexer.getTrack(1);
              auto &piTrack = fVertexer.getTrack(2);
              lVector ldeu{(float)deuTrack.Pt(), (float)deuTrack.Eta(), (float)deuTrack.Phi(), kDeuMass};
              lVector lpro{(float)prTrack.Pt(), (float)prTrack.Eta(), (float)prTrack.Phi(), kPMass};
              lVector lpi{(float)piTrack.Pt(), (float)piTrack.Eta(), (float)piTrack.Phi(), kPiMass};

              hypertriton = ldeu + lpro + lpi;

              o2RecHyp.mppi = (lpro + lpi).mass2();
              o2RecHyp.mdpi = (ldeu + lpi).mass2();
              ROOT::Math::Boost boostHyper{hypertriton.BoostToCM()};
              auto d{boostHyper(ldeu).Vect()};
              auto lambda{boostHyper(lpro + lpi).Vect()};
              auto p{boostHyper(lpro).Vect()};
              auto pi{boostHyper(lpi).Vect()};
              o2RecHyp.momDstar = std::sqrt(d.Mag2());
              o2RecHyp.cosThetaStar = d.Dot(hypertriton.Vect()) / (o2RecHyp.momDstar * hypertriton.P());
              o2RecHyp.cosTheta_ProtonPiH = p.Dot(pi) / std::sqrt(p.Mag2() * pi.Mag2());
              vert = fVertexer.getPCACandidate();
              decayVtx.SetCoordinates((float)(vert[0] - pvPos[0]), (float)(vert[1] - pvPos[1]), (float)(vert[2] - pvPos[2]));
              o2RecHyp.candidates = nVert;

              double deuPos[3], proPos[3], piPos[3];
              deuTrack.GetXYZ(deuPos);
              prTrack.GetXYZ(proPos);
              piTrack.GetXYZ(piPos);

              o2RecHyp.dca_de_pr = Hypot(deuPos[0] - proPos[0], deuPos[1] - proPos[1], deuPos[2] - proPos[2]);
              if (o2RecHyp.dca_de_pr > fMaxTrack2TrackDCA[0])
                continue;
              o2RecHyp.dca_de_pi = Hypot(deuPos[0] - piPos[0], deuPos[1] - piPos[1], deuPos[2] - piPos[2]);
              if (o2RecHyp.dca_de_pi > fMaxTrack2TrackDCA[1])
                continue;

              o2RecHyp.dca_pr_pi = Hypot(proPos[0] - piPos[0], proPos[1] - piPos[1], proPos[2] - piPos[2]);
              if (o2RecHyp.dca_pr_pi > fMaxTrack2TrackDCA[2])
                continue;

              o2RecHyp.dca_de_sv = Hypot(deuPos[0] - vert[0], deuPos[1] - vert[1], deuPos[2] - vert[2]);
              if (o2RecHyp.dca_de_sv > fMaxTrack2SVDCA[0])
                continue;
              o2RecHyp.dca_pr_sv = Hypot(proPos[0] - vert[0], proPos[1] - vert[1], proPos[2] - vert[2]);
              if (o2RecHyp.dca_pr_sv > fMaxTrack2SVDCA[1])
                continue;
              o2RecHyp.dca_pi_sv = Hypot(piPos[0] - vert[0], piPos[1] - vert[1], piPos[2] - vert[2]);
              if (o2RecHyp.dca_pi_sv > fMaxTrack2SVDCA[2])
                continue;

              o2RecHyp.chi2 = fVertexer.getChi2AtPCACandidate();
            }
            else
            {
              KFParticle kfHyperTriton{twoCandidate};
              kfHyperTriton.AddDaughter(pi.particle);
              kfRecHyp.chi2_3prongs = kfHyperTriton.GetChi2() / kfHyperTriton.GetNDF();
              if (kfRecHyp.chi2_3prongs > fMaxKFchi2[1] || kfRecHyp.chi2_3prongs < 0.)
                continue;
              double mass = kfHyperTriton.GetMass();
              if (mass < fMassWindow[0] || mass > fMassWindow[1])
                continue;
              vert[0] = kfHyperTriton.X();
              vert[1] = kfHyperTriton.Y();
              vert[2] = kfHyperTriton.Z();
              decayVtx.SetCoordinates(kfHyperTriton.X() - prodVertex.X(), kfHyperTriton.Y() - prodVertex.Y(), kfHyperTriton.Z() - prodVertex.Z());
              ROOT::Math::XYZVectorF mom{kfHyperTriton.Px(), kfHyperTriton.Py(), kfHyperTriton.Pz()};
              hypertriton.SetCoordinates(kfHyperTriton.GetPt(), kfHyperTriton.GetEta(), kfHyperTriton.GetPhi(), kHyperTritonMass);

              kfHyperTriton.SetProductionVertex(prodVertex);
              kfRecHyp.chi2_topology = kfHyperTriton.GetChi2() / kfHyperTriton.GetNDF();
              if (kfRecHyp.chi2_topology > fMaxKFchi2[2] || kfRecHyp.chi2_topology < 0.)
                continue;
            }

            std::array<AliESDtrack *, 3> tracks{(AliESDtrack *)deu.track, (AliESDtrack *)p.track, (AliESDtrack *)pi.track};
            std::array<float, 3> nSigmaTPC{deu.nSigmaTPC, p.nSigmaTPC, pi.nSigmaTPC};
            std::array<float, 3> nSigmaTOF{deu.nSigmaTOF, p.nSigmaTOF, pi.nSigmaTOF};
            if (!fillTreeInfo(tracks, nSigmaTPC, nSigmaTOF))
              continue;

            bool record{!fMC || !fOnlyTrueCandidates};
            if (fMC)
            {
              int momId = IsTrueHyperTriton3Candidate((AliESDtrack *)deu.track, (AliESDtrack *)p.track, (AliESDtrack *)pi.track, mcEvent);
              record = record || momId >= 0;
              if (record)
              {
                if (fKF)
                  FillGenHypertriton(fGenHypKF, momId, true, mcEvent);
                else
                  FillGenHypertriton(fGenHypO2, momId, true, mcEvent);
                mcMap[momId] = 1;
              }
            }
            if (record)
              fTreeHyp3->Fill();
          }
        }
      }
    }
    if (fEnableEventMixing)
      FillEventMixingPool(fEventCuts.GetCentrality(), pvPos[2], deuteronsForMixing);
  }

  if (fMC)
  {
    RHyperTriton rec;
    rec.centrality = fRecHyp->centrality;
    rec.trigger = fRecHyp->trigger;
    *fRecHyp = rec;
    for (int iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++)
    {
      AliVParticle *part = mcEvent->GetTrack(iTrack);
      if (!part)
      {
        ::Warning("AliAnalysisTaskHypertriton3::UserExec",
                  "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", iTrack);
        continue;
      }
      if (std::abs(part->Y()) > 1.)
        continue;
      if (!IsHyperTriton3(part, mcEvent))
        continue;
      if (mcMap.find(iTrack) != mcMap.end())
        continue;
      if (fKF)
        FillGenHypertriton(fGenHypKF, iTrack, false, mcEvent);
      else
        FillGenHypertriton(fGenHypO2, iTrack, false, mcEvent);
      fTreeHyp3->Fill();
    }
  }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
}

void AliAnalysisTaskHypertriton3::Terminate(Option_t *) {}

int AliAnalysisTaskHypertriton3::FindEventMixingCentBin(const float centrality)
{
  if (centrality >= 100.)
    return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskHypertriton3::FindEventMixingZBin(const float zvtx)
{
  if (zvtx > 10. || zvtx < -10.)
    return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

int AliAnalysisTaskHypertriton3::CheckPionCharge(std::vector<AliESDtrack *> tracks[2][2], AliESDv0 v0)
{

  double pP[3], nP[3];
  v0.GetPPxPyPz(pP[0], pP[1], pP[2]);
  v0.GetNPxPyPz(nP[0], nP[1], nP[2]);
  int isPiPositive = -1;
  if (tracks[1][1].size() <= int(v0.GetPindex()) && tracks[0][1].size() <= int(v0.GetNindex()))
    return isPiPositive;
  isPiPositive = (tracks[1][1].size() <= int(v0.GetPindex())) ? 0 : -1;
  isPiPositive = (tracks[0][1].size() <= int(v0.GetNindex())) ? 1 : -1;

  if (isPiPositive == -1)
  {
    double posDiff = std::abs(tracks[1][1][v0.GetPindex()]->Px() - pP[0]);
    double negDiff = std::abs(tracks[0][1][v0.GetNindex()]->Px() - nP[0]);
    isPiPositive = posDiff < negDiff ? 1 : 0;
  }
  return isPiPositive;
}

void AliAnalysisTaskHypertriton3::FillEventMixingPool(const float centrality, const float zvtx,
                                                      const std::vector<EventMixingTrack> &tracks)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  auto &trackVector = fEventMixingPool[centBin][zBin];

  for (auto &t : tracks)
    trackVector.emplace_back(t);

  while (trackVector.size() > fEventMixingPoolDepth)
    trackVector.pop_front();

  return;
}

std::vector<EventMixingTrack *> AliAnalysisTaskHypertriton3::GetEventMixingTracks(const float centrality,
                                                                                  const float zvtx)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  std::vector<EventMixingTrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin])
  {
    if (v.used >= fEventMixingPoolMaxReuse)
      continue;
    tmpVector.emplace_back(&(v));
    v.used++;
  }

  return tmpVector;
}

AliAnalysisTaskHypertriton3 *AliAnalysisTaskHypertriton3::AddTask(bool isMC, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHyperTriton2BodyML", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHypertritonO2", "This task requires an input event handler");
    return nullptr;
  }

  TString tskname = "AliAnalysisTaskHypertriton3";
  tskname.Append(suffix.Data());
  AliAnalysisTaskHypertriton3 *task = new AliAnalysisTaskHypertriton3(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("HyperTritonTree%s", suffix.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree3.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
