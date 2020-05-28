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
    float nSigmaTPC = -1.f;
    float nSigmaTOF = -1.f;
    KFParticle particle;
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

  if (fGenHypKF || fGenHypO2) {
    if (fGenHypKF) delete fGenHypKF;
    if (fGenHypO2) delete fGenHypO2;
  } else if (fRecHyp)
    delete fRecHyp;
}

void AliAnalysisTaskHypertriton3::UserCreateOutputObjects()
{
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

  fListHist->Add(fHistNSigmaDeu);
  fListHist->Add(fHistNSigmaP);
  fListHist->Add(fHistNSigmaPi);

  fListHist->Add(fHistInvMass);

  OpenFile(2);
  fTreeHyp3 = new TTree("Hyp3O2", "Hypetriton 3 Body with the O2 Vertexer");

  if (fMC && man->GetMCtruthEventHandler())
  {
    if (fKF) {
      fGenHypKF = new SHyperTriton3KF;
      fRecHyp = (RHyperTriton*)fGenHypKF;
      fTreeHyp3->Branch("SHyperTriton", fGenHypKF);
    } else {
      fGenHypO2 = new SHyperTriton3O2;
      fRecHyp = (RHyperTriton*)fGenHypO2;
      fTreeHyp3->Branch("SHyperTriton", fGenHypO2);
    }
  }
  else {
    if (fKF)
      fRecHyp = new RHyperTriton3KF;
    else
      fRecHyp = new RHyperTriton3O2;
    fTreeHyp3->Branch("RHyperTriton", fRecHyp);
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
  // set Magnetic field for KF
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
  fRecHyp->trigger |= esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  std::vector<HelperParticle> helpers[3][2];
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
          int chargeIndex = (fSwapSign && iT == 0) ? track->GetSigned1Pt() < 0 : track->GetSigned1Pt() > 0;
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
          helpers[iT][chargeIndex].push_back(helper);
        }
      }
    }
  }

  /// if event mixing is enabled takes deuteron from the event mixing pool
  // if (fEnableEventMixing && fApplyML) {
  //   deuterons = GetEventMixingTracks(fREvent.fCent, fREvent.fZ);
  // }

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

  std::unordered_map<int, int> mcMap;
  RHyperTriton3KF& kfRecHyp = *(RHyperTriton3KF*)fRecHyp;
  RHyperTriton3O2& o2RecHyp = *(RHyperTriton3O2*)fRecHyp;;
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
          if (p.track == pi.track || deu.track == pi.track)
            continue;

          lVector hypertriton;
          ROOT::Math::SVector<double, 3U> vert;
          ROOT::Math::XYZVectorF decayVtx;
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

          const float mass = hypertriton.mass();
          if (mass < fMassWindow[0] || mass > fMassWindow[1])
            continue;

          const float totalMom = hypertriton.P();
          const float len = std::sqrt(decayVtx.Mag2());
          fRecHyp->cosPA = hypertriton.Vect().Dot(decayVtx) / (totalMom * len);
          const float cosPA = fUseAbsCosPAcut ? std::abs(fRecHyp->cosPA) : fRecHyp->cosPA;
          fRecHyp->ct = len * kHyperTritonMass / totalMom;
          if (fRecHyp->ct < fCandidateCtRange[0] || fRecHyp->ct > fCandidateCtRange[1])
            continue;
          if (fCosPAspline)
          {
            if (cosPA < fCosPAspline->Eval(fRecHyp->ct))
              continue;
          }
          else if (cosPA < fMinCosPA)
          {
            continue;
          }
          fRecHyp->r = decayVtx.Rho();
          fRecHyp->positive = deu.track->Charge() > 0;
          fRecHyp->pt = hypertriton.pt();
          fRecHyp->phi = hypertriton.phi();
          fRecHyp->pz = hypertriton.pz();
          fRecHyp->m = mass;

          float dca[2], bCov[3];
          deu.track->GetImpactParameters(dca, bCov);
          fRecHyp->dca_de = std::hypot(dca[0], dca[1]);
          p.track->GetImpactParameters(dca, bCov);
          fRecHyp->dca_pr = std::hypot(dca[0], dca[1]);
          pi.track->GetImpactParameters(dca, bCov);
          fRecHyp->dca_pi = std::hypot(dca[0], dca[1]);

          fRecHyp->hasTOF_de = HasTOF(deu.track);
          fRecHyp->hasTOF_pr = HasTOF(p.track);
          fRecHyp->hasTOF_pi = HasTOF(pi.track);

          fRecHyp->tofNsig_de = deu.nSigmaTOF;
          fRecHyp->tofNsig_pr = p.nSigmaTOF;
          fRecHyp->tofNsig_pi = pi.nSigmaTOF;

          fRecHyp->tpcNsig_de = deu.nSigmaTPC;
          fRecHyp->tpcNsig_pr = p.nSigmaTPC;
          fRecHyp->tpcNsig_pi = pi.nSigmaTPC;

          fRecHyp->tpcClus_de = deu.track->GetTPCsignalN();
          fRecHyp->tpcClus_pr = p.track->GetTPCsignalN();
          fRecHyp->tpcClus_pi = pi.track->GetTPCsignalN();

          if (fLambdaCheck)
          {
            int nVertLambda{0};
            try
            {
              nVertLambda = fVertexerLambda.process(*p.track, *pi.track);
            }
            catch (std::runtime_error &e) {}

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
              fRecHyp->dca_lambda_hyper = Hypot(vert[0] - vertLambda[0], vert[1] - vertLambda[1], vert[2] - vertLambda[2]);
            }
          }

          bool record{!fMC || !fOnlyTrueCandidates};
          if (fMC)
          {
            int momId = IsTrueHyperTriton3Candidate((AliESDtrack *)deu.track, (AliESDtrack *)p.track, (AliESDtrack *)pi.track, mcEvent);
            record = record || momId >= 0;
            if (record)
            {
              if (fKF) FillGenHypertriton(fGenHypKF, momId, true, mcEvent);
              else FillGenHypertriton(fGenHypO2, momId, true, mcEvent);
              mcMap[momId] = 1;
            }
          }
          if (record)
            fTreeHyp3->Fill();
        }
      }
    }
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
      if (fKF) FillGenHypertriton(fGenHypKF, iTrack, false, mcEvent);
      else FillGenHypertriton(fGenHypO2, iTrack, false, mcEvent);
      fTreeHyp3->Fill();
    }
  }

  /// if event mixing is enabled fill the event mixing pool with deuterons
  // if (fEnableEventMixing && fApplyML) {
  //   FillEventMixingPool(fREvent.fCent, fREvent.fZ, fDeuVector);
  // }

  PostData(1, fListHist);
  PostData(2, fTreeHyp3);
}

void AliAnalysisTaskHypertriton3::Terminate(Option_t *) {}

int AliAnalysisTaskHypertriton3::FindEventMixingCentBin(const float centrality)
{
  if (centrality > 90)
    return -999;
  return static_cast<int>(centrality / 10);
}

int AliAnalysisTaskHypertriton3::FindEventMixingZBin(const float zvtx)
{
  if (zvtx > 10. || zvtx < -10.)
    return -999.;
  return static_cast<int>((zvtx + 10.) / 2);
}

void AliAnalysisTaskHypertriton3::FillEventMixingPool(const float centrality, const float zvtx,
                                                       std::vector<AliESDtrack *> tracks)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  auto &trackVector = fEventMixingPool[centBin][zBin];

  for (auto &t : tracks)
  {
    trackVector.emplace_back(AliESDtrack{*t});
  }

  if (trackVector.size() - fEventMixingPoolDepth > 0)
    trackVector.pop_front();

  return;
}

std::vector<AliESDtrack *> AliAnalysisTaskHypertriton3::GetEventMixingTracks(const float centrality,
                                                                              const float zvtx)
{
  int centBin = FindEventMixingCentBin(centrality);
  int zBin = FindEventMixingZBin(zvtx);

  std::vector<AliESDtrack *> tmpVector;

  for (auto &v : fEventMixingPool[centBin][zBin])
  {
    tmpVector.emplace_back(&v);
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
  mgr->SetDebugLevel(2);

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
                           AliAnalysisManager::kOutputContainer, Form("HyperTritonTree.root:%s", suffix.Data()));
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}
