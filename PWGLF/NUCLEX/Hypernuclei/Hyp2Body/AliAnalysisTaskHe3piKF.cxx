#include "AliAnalysisTaskHe3piKF.h"

#include "AliAnalysisDataContainer.h"
#include "AliLog.h"

#include <algorithm>
#include <cmath>
#include <string>
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>
#include <TRandom3.h>
#include <TTree.h>
#include "Math/GenVector/Boost.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliESDEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
#include "AliAnalysisUtils.h"
#include "Math/Vector4D.h"
#define HomogeneousField
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskHe3piKF);
///\endcond

namespace
{

  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

  double Sq(double x)
  {
    return x * x;
  }

  constexpr int kHyperPdg{1010010030};
  constexpr double kHyperMass{2.99131}; /// from AliMC.cxx

}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
struct HelperParticle
{
  AliESDtrack *track{nullptr};
  float nSigmaTPC = -1.f;
  KFParticle particle = KFParticle();
};

constexpr float kHe3Mass{2.80839161};
constexpr float kPiMass{0.13957};
constexpr float kMasses[3]{kHe3Mass, kPiMass};
constexpr AliPID::EParticleType kAliPID[2]{AliPID::kHe3, AliPID::kPion};
const int kPDGs[2]{AliPID::ParticleCode(kAliPID[0]), AliPID::ParticleCode(kAliPID[1])};

AliAnalysisTaskHe3piKF::AliAnalysisTaskHe3piKF(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                              fEventCut{false},
                                                                              fMC{isMC}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskHe3piKF::~AliAnalysisTaskHe3piKF()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;
  if (fList)
    delete fList;
  if (fTree)
    delete fTree;
  if (!fMC)
  {
    delete fRecHyper;
  }
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskHe3piKF::UserCreateOutputObjects()
{
  fList = new TList();
  fList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fList);
  fRecHyper = fMC ? &fGenHyper : new MiniHyperKF;

  OpenFile(2);
  fTree = new TTree("HyperTree", "HyperTree");
  if (fMC)
  {
    fTree->Branch("HyperMC", &fGenHyper);
    fMCEvent = MCEvent();
  }
  else
  {
    fTree->Branch("MiniHyper", fRecHyper);
  }

  PostAllData();
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskHe3piKF::UserExec(Option_t *)
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

  if (!fEventCut.AcceptEvent(esdEvent))
  {
    PostAllData();
    return;
  }

  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3], pvCov[6];
  ;
  pvObj->GetXYZ(pv);
  fEventCut.GetPrimaryVertex()->GetCovarianceMatrix(pvCov);

  fRecHyper->fZ = pv[2];

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  fRecHyper->centrality = fEventCut.GetCentrality();

  unsigned char tgr{0u};

  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7)
    tgr |= kINT7;
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral)
    tgr |= kCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    tgr |= kSemiCentral;
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
    tgr |= kHighMultV0;
  int magField = esdEvent->GetMagneticField() > 0 ? kPositiveB : 0;

  KFParticle::SetField(magField);

  fRecHyper->trigger = tgr + magField;

  std::vector<int> checkedHyperLabel;
  fGenHyper.isReconstructed = true;

  std::vector<HelperParticle> helpers[2][2];
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

    bool candidate[2]{false};
    float nSigmasTPC[2]{-1., -1.};
    float dca[2];
    track->GetImpactParameters(dca[0], dca[1]);
    double dcaNorm = std::hypot(dca[0], dca[1]);

    for (int iT{0}; iT < 2; ++iT)
    {
      nSigmasTPC[iT] = (fMC || iT == kPion) ? fPID->NumberOfSigmasTPC(track, kAliPID[iT]) : customNsigma(track->GetTPCmomentum(), track->GetTPCsignal()); // use custom nSigma for He3 from data
      if (std::abs(nSigmasTPC[iT]) < 5 && dcaNorm > fMinTrackDCA[iT] && track->Pt() < fTrackPtRange[iT][1] &&
          track->Pt() > fTrackPtRange[iT][0] && track->GetTPCsignalN() >= 50 && std::abs(track->Eta()) < 0.8 && (track->GetStatus() & AliVTrack::kTPCrefit))
        candidate[iT] = true;
    }

    if (candidate[0] || candidate[1])
    {
      int chargeIndex = track->GetSigned1Pt() > 0;
      HelperParticle helper;
      for (int iT{0}; iT < 2; ++iT)
      {
        if (candidate[iT])
        {
          helper.nSigmaTPC = nSigmasTPC[iT];
          double posmom[6], cov[21];
          track->GetXYZ(posmom);
          track->GetPxPyPz(posmom + 3);
          int charge = track->Charge();
          track->GetCovarianceXYZPxPyPz(cov);
          if (iT == kHe3)
          {
            posmom[3] *= 2;
            posmom[4] *= 2;
            posmom[5] *= 2;
            charge *= 2;
            for (int i = 6; i < 21; i++)
            {
              cov[i] = cov[i] * 2; /// scale mom space entries of cov matrix by 2
              if (i == 9 || i == 13 || i == 14 || i == 18 || i == 19 || i == 20)
              {
                cov[i] = cov[i] * 2; /// scale mom mom entries of cov matrix by 4
              }
            }
          }
          
          helper.particle.Create(posmom, cov, charge, kMasses[iT]);
          helper.particle.Chi2() = track->GetTPCchi2();
          helper.particle.NDF() = track->GetNumberOfTPCClusters() * 2;
          helper.track = track;
          helpers[iT][chargeIndex].push_back(helper);
        }
      }
    }
  }

  KFPVertex kfPVertex;
  kfPVertex.SetXYZ(pv[0], pv[1], pv[2]);
  kfPVertex.SetCovarianceMatrix(pvCov[0], pvCov[1], pvCov[2], pvCov[3], pvCov[4], pvCov[5]);
  kfPVertex.SetChi2(fEventCut.GetPrimaryVertex()->GetChi2());
  kfPVertex.SetNDF(fEventCut.GetPrimaryVertex()->GetNDF());
  kfPVertex.SetNContributors(fEventCut.GetPrimaryVertex()->GetNContributors());
  KFParticle prodVertex{kfPVertex};
  for (int idx{0}; idx < 2; ++idx)
  {

    for (const auto &he3 : helpers[kHe3][1 - idx])
    {

      KFParticle he3Candidate{he3.particle};

      for (const auto &pi : helpers[kPion][idx])
      {
        int hyperLabel{-1};
        if (fMC)
        {
          fGenHyper.pdg = 0;
          auto he3Part = fMCEvent->GetTrack(std::abs(he3.track->GetLabel()));
          auto piPart = fMCEvent->GetTrack(std::abs(pi.track->GetLabel()));

          if ((abs(he3Part->PdgCode()) == AliPID::ParticleCode(AliPID::kHe3) &&
               abs(piPart->PdgCode()) == AliPID::ParticleCode(AliPID::kPion)))
          {
            // Check hyper
            int labMothHe3 = he3Part->GetMother();
            int labMothPi = piPart->GetMother();

            if (labMothPi >= 0 && labMothPi == labMothHe3)
            {
              auto hyper = fMCEvent->GetTrack(labMothPi);
              if (hyper && std::abs(hyper->PdgCode()) == kHyperPdg)
              {
                hyperLabel = labMothPi;
                fGenHyper.pdg = hyper->PdgCode();
                fGenHyper.ptMC = hyper->Pt();
                fGenHyper.etaMC = hyper->Eta();
                fGenHyper.yMC = hyper->Y();
                double ov[3], dv[3];
                hyper->XvYvZv(ov);
                he3Part->XvYvZv(dv);
                fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * hyper->M() / hyper->P();
                // std::cout << "-------------------------" << std::endl;
                // std::cout << "Gen, Prim VTx x: " << ov[0] << ", y: " << ov[1] << ", z: " << ov[2] << std::endl;
                // std::cout << "Gen, Decay VTx x: " << dv[0] << ", y: " << dv[1] << ", z: " << dv[2] << std::endl;
              }
            }
          }
          if (fOnlyTrueCandidates && hyperLabel < 0)
            continue;
        }

        KFParticle kfHyperTriton;
        if (fMassConstrainedFit)
        {
          kfHyperTriton.SetConstructMethod(2);
        }
        kfHyperTriton.AddDaughter(he3Candidate);
        kfHyperTriton.AddDaughter(pi.particle);

        double recoMass = kfHyperTriton.GetMass();
        if (recoMass > fMassRange[1] || recoMass < fMassRange[0])
          continue;

        ROOT::Math::XYZVectorF decayVtx;
        double deltaPos[3]{kfHyperTriton.X() - prodVertex.X(), kfHyperTriton.Y() - prodVertex.Y(), kfHyperTriton.Z() - prodVertex.Z()};
        decayVtx.SetCoordinates(kfHyperTriton.X() - prodVertex.X(), kfHyperTriton.Y() - prodVertex.Y(), kfHyperTriton.Z() - prodVertex.Z());

        // std::cout << "Rec, Prim VTx x: " << prodVertex.X() << ", y: " << prodVertex.Y() << ", z: " << prodVertex.Z() << std::endl;
        // std::cout << "Rec, Dec VTx x: " << kfHyperTriton.X() << ", y: " << kfHyperTriton.Y() << ", z: " << kfHyperTriton.Z() << std::endl;
        // std::cout << "Charge: " << int(kfHyperTriton.GetQ()) << std::endl;

        fRecHyper->fChi2 = kfHyperTriton.GetChi2() / kfHyperTriton.GetNDF();
        if (fRecHyper->fChi2 > fMaxKFchi2 || fRecHyper->fChi2 < 0.)
          continue;

        fRecHyper->Matter = idx == 0 ? true : false;
        auto &he3track = he3.track;
        auto &pitrack = pi.track;
        LVector_t hyperVector;
        hyperVector.SetCoordinates(kfHyperTriton.Px(), kfHyperTriton.Py(), kfHyperTriton.Pz(), kfHyperTriton.GetMass());

        double prongsDCA = he3Candidate.GetDistanceFromParticle(pi.particle);
        double cpa = (deltaPos[0] * hyperVector.px() +
                      deltaPos[1] * hyperVector.py() +
                      deltaPos[2] * hyperVector.pz()) /
                     std::sqrt(hyperVector.P2() * (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2])));

        // std::cout << "-----------------------" << std::endl;
        // std::cout  << std::setprecision(9) << "Den: " << std::sqrt(hyperVector.P2() * (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2]))) << std::endl;
        // std::cout  << std::setprecision(9) << "Num: " << (deltaPos[0] * hyperVector.px() +
        //               deltaPos[1] * hyperVector.py() +
        //               deltaPos[2] * hyperVector.pz())   << std::endl;
        // std::cout << std::setprecision(9) <<"CosPA: " << cpa << std::endl;
        kfHyperTriton.SetProductionVertex(prodVertex);

        fRecHyper->Lrec = sqrt(decayVtx.Mag2());
        fRecHyper->ct = fRecHyper->Lrec * kHyperMass / hyperVector.P();
        if (hyperVector.pt() < fCandidatePtRange[0] || hyperVector.pt() > fCandidatePtRange[1] || fRecHyper->ct < fCandidateCtRange[0] || fRecHyper->ct > fCandidateCtRange[1])
          continue;

        if (cpa < fMinCosPA || prongsDCA > fMaxProngDCA)
          continue;

        fRecHyper->pt = hyperVector.pt();
        fRecHyper->m = recoMass;
        fRecHyper->V0CosPA = cpa;
        fRecHyper->Rapidity = Eta2y(fRecHyper->pt, kHyperMass, hyperVector.eta());
        fRecHyper->V0radius = decayVtx.Rho();
        float z, xy;
        he3track->GetImpactParameters(xy, z);
        fRecHyper->He3ProngPvDCAXY = xy;
        fRecHyper->He3ProngPvDCA = sqrt(xy * xy + z * z);
        pitrack->GetImpactParameters(xy, z);
        fRecHyper->PiProngPvDCAXY = xy;
        fRecHyper->PiProngPvDCA = sqrt(xy * xy + z * z);

        fRecHyper->ProngsDCA = he3Candidate.GetDistanceFromParticle(pi.particle);
        fRecHyper->TPCmomHe3 = he3track->GetTPCmomentum();
        fRecHyper->TPCsignalHe3 = he3track->GetTPCsignal();
        fRecHyper->NitsClustersHe3 = he3track->GetITSNcls();
        fRecHyper->TPCnSigmaPi = pi.nSigmaTPC;
        fRecHyper->TPCnSigmaHe3 = he3.nSigmaTPC;
        fRecHyper->NpidClustersPion = pitrack->GetTPCsignalN();
        fRecHyper->NpidClustersHe3 = he3track->GetTPCsignalN();

        if (hyperLabel != -1)
        {
          if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), hyperLabel) != checkedHyperLabel.end())
          {
            fGenHyper.isDuplicated = true;
          }
          else
          {
            fGenHyper.isDuplicated = false;
            checkedHyperLabel.push_back(hyperLabel);
          }
        }
        fTree->Fill();
      }
    }
  }
  if (fMC)
  {
    fGenHyper.isReconstructed = false;
    // loop on generated
    for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
    {
      auto track = (AliVParticle *)fMCEvent->GetTrack(iT);
      int pdg = std::abs(track->PdgCode());
      if (pdg != kHyperPdg)
      {
        continue;
      }
      if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), iT) != checkedHyperLabel.end())
      {
        continue;
      }

      // if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
      // {
      //   continue;
      // }

      fGenHyper.ptMC = track->Pt();
      fGenHyper.etaMC = track->Eta();
      fGenHyper.yMC = track->Y();
      fGenHyper.pdg = track->PdgCode();
      double ov[3], dv[3];
      track->XvYvZv(ov);
      bool otherDecayChannel{true};
      for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
      {
        auto daugh = (AliVParticle *)fMCEvent->GetTrack(iD);
        if (!daugh)
        {
          continue;
        }
        if (std::abs(daugh->PdgCode()) == AliPID::ParticleCode(AliPID::kHe3))
        {
          otherDecayChannel = false;
          daugh->XvYvZv(dv);
          break;
        }
      }
      if (otherDecayChannel)
        continue;
      fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * track->M() / track->P();
      fTree->Fill();
    }
  }

  PostAllData();
}

AliAnalysisTaskHe3piKF *AliAnalysisTaskHe3piKF::AddTask(bool isMC, TString tskname, TString suffix)
{
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskStrangenessRatios", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskStrangenessRatios", "This task requires an input event handler");
    return nullptr;
  }

  tskname.Append(suffix.Data());
  AliAnalysisTaskHe3piKF *task = new AliAnalysisTaskHe3piKF(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("%s_tree", tskname.Data()), TTree::Class(),
                                                            AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);

  return task;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskHe3piKF::Eta2y(float pt, float m, float eta) const
{
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskHe3piKF::PostAllData()
{
  PostData(1, fList);
  PostData(2, fTree);
}

void AliAnalysisTaskHe3piKF::SetCustomBetheBloch(float res, const float *bethe)
{
  fCustomResolution = res;
  std::copy(bethe, bethe + 5, fCustomBethe);
}

double AliAnalysisTaskHe3piKF::customNsigma(double mom, double sig)
{
  const float bg = mom / AliPID::ParticleMass(AliPID::kHe3);
  const float *p = fCustomBethe;
  const float expS = AliExternalTrackParam::BetheBlochAleph(bg, p[0], p[1], p[2], p[3], p[4]);
  return (sig - expS) / (fCustomResolution * expS);
}