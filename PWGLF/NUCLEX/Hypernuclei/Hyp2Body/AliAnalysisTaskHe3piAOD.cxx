#include "AliAnalysisTaskHe3piAOD.h"

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

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisUtils.h"

#include "Math/Vector4D.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskHe3piAOD);
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
AliAnalysisTaskHe3piAOD::AliAnalysisTaskHe3piAOD(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                fEventCut{false},
                                                                                fMC{isMC}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskHe3piAOD::~AliAnalysisTaskHe3piAOD()
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
void AliAnalysisTaskHe3piAOD::UserCreateOutputObjects()
{

  fList = new TList();
  fList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fList);
  fRecHyper = fMC ? &fGenHyper : new MiniHyper;

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
void AliAnalysisTaskHe3piAOD::UserExec(Option_t *)
{
  AliAODEvent *ev = (AliAODEvent *)InputEvent();
  if (!fEventCut.AcceptEvent(ev))
  {
    PostAllData();
    return;
  }

  AliAODMCHeader *header{nullptr};
  TClonesArray *MCTrackArray{nullptr};
  if (fMC) {
    //OOB pileup
    header = static_cast<AliAODMCHeader *>(ev->FindListObject(AliAODMCHeader::StdBranchName()));
    if (!header)
    {
      AliWarning("No header found.");
      PostAllData();
      return;
    }
    MCTrackArray = dynamic_cast<TClonesArray *>(ev->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!MCTrackArray)
    {
      AliWarning("No MC track array found.");
      PostAllData();
      return;
    }
  }

  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3];
  pvObj->GetXYZ(pv);
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
  int magField = ev->GetMagneticField() > 0 ? kPositiveB : 0;

  fRecHyper->trigger = tgr + magField;

  std::vector<int> checkedHyperLabel;
  fGenHyper.isReconstructed = true;

  for (int iV0{0}; iV0 < ev->GetNumberOfV0s(); ++iV0)
  {
    AliAODv0 *v0{ev->GetV0(iV0)};
    if (!v0)
      continue;
    if (v0->GetOnFlyStatus() != fUseOnTheFly)
      continue;
    //get daughter tracks (positive, negative and bachelor)
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
    if (!pTrack || !nTrack)
    {
      AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks of the lambdas ...\n");
      continue;
    }

    if (!(pTrack->GetStatus() & AliVTrack::kTPCrefit) || !(nTrack->GetStatus() & AliVTrack::kTPCrefit) ||
        pTrack->GetTPCsignalN() < 50 || nTrack->GetTPCsignalN() < 50 ||
        std::abs(pTrack->Eta()) > 0.8 || std::abs(nTrack->Eta()) > 0.8)
    {
      continue;
    }

    int hyperLabel{-1};
    if (fMC)
    {
      fGenHyper.pdg = 0;
      auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(pTrack->GetLabel()));
      auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(std::abs(nTrack->GetLabel()));
      if ((posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kHe3) &&
           negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kPion)) ||
          (negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kHe3) &&
           posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kPion)))
      {
        // Check hyper
        int labMothPos = posPart->GetMother();
        int labMothNeg = negPart->GetMother();
        if (labMothNeg >= 0 && labMothNeg == labMothPos && !AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(labMothNeg, header, MCTrackArray))
        {
          auto hyper = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
          if (hyper && std::abs(hyper->GetPdgCode()) == kHyperPdg)
          {
            hyperLabel = labMothNeg;
            fGenHyper.pdg = hyper->GetPdgCode();
            fGenHyper.ptMC = hyper->Pt();
            fGenHyper.etaMC = hyper->Eta();
            fGenHyper.yMC = hyper->Y();
            double ov[3], dv[3];
            hyper->XvYvZv(ov);
            posPart->XvYvZv(dv);
            fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) * hyper->M() / hyper->P();
          }
        }
      }

      if (fOnlyTrueCandidates && hyperLabel < 0)
        continue;
    }

    double pNsigma{fMC ? fPID->NumberOfSigmasTPC(pTrack, AliPID::kHe3) : customNsigma(pTrack->GetTPCmomentum(), pTrack->GetTPCsignal())};
    double nNsigma{fMC ? fPID->NumberOfSigmasTPC(nTrack, AliPID::kHe3) : customNsigma(nTrack->GetTPCmomentum(), nTrack->GetTPCsignal())};
    if (std::abs(pNsigma) > 5 && std::abs(nNsigma) > 5)
    {
      continue;
    }
    if (std::abs(pNsigma) < 5 && std::abs(nNsigma) < 5)
    {
      continue;
    }

    fRecHyper->Matter = std::abs(pNsigma) < 5;
    auto he3 = fRecHyper->Matter ? pTrack : nTrack;
    auto pion = fRecHyper->Matter ? nTrack : pTrack;

    double sv[3]{v0->GetSecVtxX(), v0->GetSecVtxY(), v0->GetSecVtxZ()};
    double deltaPos[3]{sv[0] - pv[0], sv[1] - pv[1], sv[2] - pv[2]};

    LVector_t he3Vector, piVector, hyperVector;
    double he3P[3]{fRecHyper->Matter ? v0->MomPosX() : v0->MomNegX(), fRecHyper->Matter ? v0->MomPosY() : v0->MomNegY(), fRecHyper->Matter ? v0->MomPosZ() : v0->MomNegZ()};
    double piP[3]{!fRecHyper->Matter ? v0->MomPosX() : v0->MomNegX(), !fRecHyper->Matter ? v0->MomPosY() : v0->MomNegY(), !fRecHyper->Matter ? v0->MomPosZ() : v0->MomNegZ()};

    he3Vector.SetCoordinates(he3P[0] * 2, he3P[1] * 2, he3P[2] * 2, AliPID::ParticleMass(AliPID::kHe3));
    piVector.SetCoordinates(piP[0], piP[1], piP[2], AliPID::ParticleMass(AliPID::kPion));
    hyperVector = piVector + he3Vector;
    if (hyperVector.mass() > fMassRange[1] || hyperVector.mass() < fMassRange[0])
    {
      continue;
    }

    double cpa = (deltaPos[0] * hyperVector.px() +
                  deltaPos[1] * hyperVector.py() +
                  deltaPos[2] * hyperVector.pz()) /
                 std::sqrt(hyperVector.P2() * (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2])));

    fRecHyper->pt = hyperVector.pt();
    fRecHyper->m = hyperVector.mass();
    fRecHyper->V0CosPA = cpa;
    fRecHyper->Rapidity = Eta2y(fRecHyper->pt, kHyperMass, hyperVector.eta());

    fRecHyper->V0radius = v0->RadiusSecVtx();
    fRecHyper->Lrec = v0->DecayLengthV0(pv);
    fRecHyper->ct = fRecHyper->Lrec * kHyperMass / hyperVector.P();

    fRecHyper->He3ProngPvDCA = fRecHyper->Matter ? v0->DcaPosToPrimVertex() : v0->DcaNegToPrimVertex();
    fRecHyper->PiProngPvDCA = fRecHyper->Matter ? v0->DcaNegToPrimVertex() : v0->DcaPosToPrimVertex();
    float _dummy, xy;
    he3->GetImpactParameters(xy, _dummy);
    fRecHyper->He3ProngPvDCAXY = xy;
    pion->GetImpactParameters(xy, _dummy);
    fRecHyper->PiProngPvDCAXY = xy;
    fRecHyper->ProngsDCA = v0->DcaV0Daughters();
    fRecHyper->TPCmomHe3 = he3->GetTPCmomentum();
    fRecHyper->TPCsignalHe3 = he3->GetTPCsignal();
    fRecHyper->NitsClustersHe3 = he3->GetITSNcls();
    fRecHyper->TPCnSigmaPi = fPID->NumberOfSigmasTPC(pion, AliPID::kPion);
    fRecHyper->TPCnSigmaHe3 = fRecHyper->Matter ? pNsigma : nNsigma;
    fRecHyper->NpidClustersPion = pion->GetTPCsignalN();
    fRecHyper->NpidClustersHe3 = he3->GetTPCsignalN();

    if (hyperLabel != -1)
    {
      if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), hyperLabel) != checkedHyperLabel.end())
      {
        fGenHyper.isDuplicated = true;
      } else {
        fGenHyper.isDuplicated = false;
        checkedHyperLabel.push_back(hyperLabel);
      }
    }
    fTree->Fill();
  }

  if (fMC)
  {
    fGenHyper.isReconstructed = false;
    //loop on generated
    for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT)
    {
      auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
      int pdg = std::abs(track->GetPdgCode());
      if (pdg != kHyperPdg)
      {
        continue;
      }
      if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), iT) != checkedHyperLabel.end())
      {
        continue;
      }

      if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iT, header, MCTrackArray))
      {
        continue;
      }
      fGenHyper.ptMC = track->Pt();
      fGenHyper.etaMC = track->Eta();
      fGenHyper.yMC = track->Y();
      fGenHyper.pdg = track->GetPdgCode();
      double ov[3], dv[3];
      track->XvYvZv(ov);
      bool otherDecayChannel{true};
      for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++)
      {
        auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
        if (!daugh)
        {
          continue;
        }
        if (std::abs(daugh->GetPdgCode()) == AliPID::ParticleCode(AliPID::kHe3))
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

AliAnalysisTaskHe3piAOD *AliAnalysisTaskHe3piAOD::AddTask(bool isMC, TString tskname, TString suffix)
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
  AliAnalysisTaskHe3piAOD *task = new AliAnalysisTaskHe3piAOD(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("%s_treeCascades", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer(Form("%s_treeLambda", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput3->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  return task;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskHe3piAOD::Eta2y(float pt, float m, float eta) const
{
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskHe3piAOD::PostAllData()
{
  PostData(1, fList);
  PostData(2, fTree);
}

void AliAnalysisTaskHe3piAOD::SetCustomBetheBloch(float res, const float *bethe)
{
  fCustomResolution = res;
  std::copy(bethe, bethe + 5, fCustomBethe);
}

double AliAnalysisTaskHe3piAOD::customNsigma(double mom, double sig)
{
  const float bg = mom / AliPID::ParticleMass(AliPID::kHe3);
  const float *p = fCustomBethe;
  const float expS = AliExternalTrackParam::BetheBlochAleph(bg, p[0], p[1], p[2], p[3], p[4]);
  return (sig - expS) / (fCustomResolution * expS);
}