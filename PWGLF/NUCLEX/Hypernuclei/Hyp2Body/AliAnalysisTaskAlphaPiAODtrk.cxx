/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
// AliAnalysisTaskAlphaPiAODtrk class
// analysis task for the Hyper Hydrogen 4 Analysis
// Based on the track-track loop
// Copied from AliAnalysisTaskAlphaPiAODtrk class
//
// Author:
// Bong-Hwi Lim bong-hwi.lim@cern.ch
///////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskAlphaPiAODtrk.h"
#include "AliAnalysisTaskNucleiYield.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "AliAnalysisDataContainer.h"
#include "AliLog.h"
using std::string;

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH2F.h>
#include <TList.h>
#include <TRandom3.h>
#include <TTree.h>

// ALIROOT includes
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliVEventHandler.h"
#include "AliVTrack.h"
#include "Math/Vector4D.h"
#include "THistManager.h"

///\cond CLASSIMP
ClassImp(AliAnalysisTaskAlphaPiAODtrk);
///\endcond

namespace {

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

double Sq(double x) { return x * x; }

constexpr int kHyperPdg{1010010040};
constexpr double kHyperMass{3.931}; /// from AliPDG.cxx

} // namespace

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskAlphaPiAODtrk::AliAnalysisTaskAlphaPiAODtrk(bool isMC,
                                                           TString taskname)
    : AliAnalysisTaskSE(taskname.Data()), fEventCut{false}, fHistos{nullptr},
      fMC{isMC} {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskAlphaPiAODtrk::~AliAnalysisTaskAlphaPiAODtrk() {
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;
  if (fHistos)
    delete fHistos;
  if (fTree)
    delete fTree;
  if (!fMC) {
    delete fRecHyper;
  }
}

/// This function creates all the histograms and all the objects in general used
/// during the analysis \return void
///
void AliAnalysisTaskAlphaPiAODtrk::UserCreateOutputObjects() {
  fHistos = new THistManager("hyperhists");
  fEventCut.AddQAplotsToList(
      fHistos->GetListOfHistograms()); // EventCuts QA Histograms

  // QA Histograms
  fHistos->CreateTH2("QA/hTPCPIDAllTracksP",
                     "TPC PID all tracks;#it{p}_{T} (GeV/#it{c});TPC Signal",
                     400, -20, 20, 1000, 0, 1000);
  fHistos->CreateTH2("QA/hTPCPIDAllTracksN",
                     "TPC PID all tracks;#it{p}_{T} (GeV/#it{c});TPC Signal",
                     400, -20, 20, 1000, 0, 1000);
  fHistos->CreateTH2("QA/hTPCPIDAlpha",
                     "TPC PID #alpha;#it{p}_{T} (GeV/#it{c});TPC Signal #alpha",
                     400, -20, 20, 1000, 0, 1000);
  fHistos->CreateTH2(
      "QA/hTPCPIDAntiAlpha",
      "TPC PID #bar{#alpha};#it{p}_{T} (GeV/#it{c});TPC Signal #bar{#alpha}",
      400, -20, 20, 1000, 0, 1000);
  fHistos->CreateTH2(
      "QA/hTPCPIDnSigmaAlpha",
      "TPC PID #alpha;#it{p}_{T} (GeV/#it{c}); n_{#sigma} #alpha", 200, 0, 20,
      100, -5, 5);
  fHistos->CreateTH2(
      "QA/hTPCPIDnSigmaAntiAlpha",
      "TPC PID #bar{#alpha};#it{p}_{T} (GeV/#it{c}); n_{#sigma} #bar{#alpha}",
      200, -20, 0, 100, -5, 5);

  // MC QA Histograms
  // if (fMC) {
  //     fHistos->CreateTH2("QA_MC/hTPCPIDAlpha", "TPC PID #alpha;#it{p}_{T}
  //     (GeV/#it{c});n_{#sigma} #alpha", 200, 0, 20, 200, 0, 200);
  //     fHistos->CreateTH2("QA_MC/hTPCPIDAntiAlpha", "TPC PID
  //     #bar{#alpha};#it{p}_{T} (GeV/#it{c});n_{#sigma} #bar{#alpha}", 200, 0,
  //     20, 200, 0, 200);
  // }

  fRecHyper = fMC ? &fGenHyper : new StructHypertrk;
  OpenFile(2);
  fTree = new TTree("HyperTree", "HyperTree");
  if (fMC) {
    fTree->Branch("HyperMC", &fGenHyper);
    fMCEvent = MCEvent();
  } else {
    fTree->Branch("StructHypertrk", fRecHyper);
  }

  PostAllData();
}

/// This is the function that is evaluated for each event. The analysis code
/// stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskAlphaPiAODtrk::UserExec(Option_t *) {
  AliAODEvent *ev = (AliAODEvent *)InputEvent();
  if (!fEventCut.AcceptEvent(ev)) {
    PostAllData();
    return;
  }

  AliAODMCHeader *header{nullptr};
  TClonesArray *MCTrackArray{nullptr};
  if (fMC) {
    // OOB pileup
    header = static_cast<AliAODMCHeader *>(
        ev->FindListObject(AliAODMCHeader::StdBranchName()));
    if (!header) {
      AliWarning("No header found.");
      PostAllData();
      return;
    }
    MCTrackArray = dynamic_cast<TClonesArray *>(
        ev->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!MCTrackArray) {
      AliWarning("No MC track array found.");
      PostAllData();
      return;
    }
  }

  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3];
  pvObj->GetXYZ(pv);
  fRecHyper->fZ = pv[2];

  /// To perform the majority of the analysis - and also this one - the standard
  /// PID handler is required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl =
      (AliInputEventHandler *)mgr->GetInputEventHandler();
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
  // Trackloop method (for trial)
  // iTrack loop -> alpha
  AliESDtrack *pionESDTrack = 0x0;
  AliESDtrack *alphaESDTrack = 0x0;
  for (int iTrack{0}; iTrack < ev->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *alphaTrack = dynamic_cast<AliAODTrack *>(ev->GetTrack(iTrack));
    if (!alphaTrack) {
      AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter tracks "
                 "of the lambdas ...\n");
      continue;
    }
    if (!((AliAODTrack *)alphaTrack)->TestFilterBit(fFilterBit))
      continue;
    if (std::abs(alphaTrack->Eta()) > 0.8)
      continue;

    if (alphaTrack->GetSign() > 0)
      fHistos->FillTH2("QA/hTPCPIDAllTracksP",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       alphaTrack->GetTPCsignal());
    else
      fHistos->FillTH2("QA/hTPCPIDAllTracksN",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       alphaTrack->GetTPCsignal());
    double pNsigma{fPID->NumberOfSigmasTPC(alphaTrack, fNucleusPID)};
    if (pNsigma < fPIDrange[0] || pNsigma > fPIDrange[1]) {
      continue;
    }
    if (alphaTrack->GetSign() > 0) {
      fHistos->FillTH2("QA/hTPCPIDAlpha",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       alphaTrack->GetTPCsignal());
      fHistos->FillTH2("QA/hTPCPIDnSigmaAlpha",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       pNsigma);
    } else {
      fHistos->FillTH2("QA/hTPCPIDAntiAlpha",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       alphaTrack->GetTPCsignal());
      fHistos->FillTH2("QA/hTPCPIDnSigmaAntiAlpha",
                       alphaTrack->GetTPCmomentum() / alphaTrack->Charge(),
                       pNsigma);
    }

    // jTrack loop - pion
    for (int jTrack{0}; jTrack < ev->GetNumberOfTracks(); ++jTrack) {
      if (jTrack == iTrack)
        continue;
      AliAODTrack *pionTrack =
          dynamic_cast<AliAODTrack *>(ev->GetTrack(jTrack));
      if (!pionTrack) {
        AliWarning("ERROR: Could not retrieve one of the 2 AOD daughter "
                   "tracks of the lambdas ...\n");
        continue;
      }
      if (!((AliAODTrack *)pionTrack)->TestFilterBit(fFilterBit))
        continue;
      if (std::abs(pionTrack->Eta()) > 0.8)
        continue;

      double nNsigma{fPID->NumberOfSigmasTPC(pionTrack, AliPID::kPion)};
      if (std::abs(nNsigma) > 5) {
        continue;
      }

      int hyperLabel{-1};
      if (fMC) {
        fGenHyper.pdg = 0;
        auto posPart = (AliAODMCParticle *)fMCEvent->GetTrack(
            std::abs(alphaTrack->GetLabel()));
        auto negPart = (AliAODMCParticle *)fMCEvent->GetTrack(
            std::abs(pionTrack->GetLabel()));
        if ((posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kAlpha) &&
             negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kPion)) ||
            (negPart->GetPdgCode() == -AliPID::ParticleCode(AliPID::kAlpha) &&
             posPart->GetPdgCode() == AliPID::ParticleCode(AliPID::kPion))) {
          // Check hyper
          int labMothPos = posPart->GetMother();
          int labMothNeg = negPart->GetMother();
          if (labMothNeg >= 0 && labMothNeg == labMothPos &&
              !AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(
                  labMothNeg, header, MCTrackArray)) {
            auto hyper = (AliAODMCParticle *)fMCEvent->GetTrack(labMothNeg);
            if (hyper && std::abs(hyper->GetPdgCode()) == kHyperPdg) {
              hyperLabel = labMothNeg;
              fGenHyper.pdg = hyper->GetPdgCode();
              fGenHyper.ptMC = hyper->Pt();
              fGenHyper.etaMC = hyper->Eta();
              fGenHyper.yMC = hyper->Y();
              double ov[3], dv[3];
              hyper->XvYvZv(ov);
              posPart->XvYvZv(dv);
              fGenHyper.ctMC = std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) +
                                         Sq(ov[2] - dv[2])) *
                               hyper->M() / hyper->P();
            }
          }
        }
        if (fOnlyTrueCandidates && hyperLabel < 0)
          continue;
      }

      fRecHyper->Matter = alphaTrack->GetSign() > 0;

      // Vertexing
      AliVTrack *aAlphaTrack = (AliVTrack *)ev->GetTrack(iTrack);
      AliVTrack *aPionTrack = (AliVTrack *)ev->GetTrack(jTrack);
      alphaESDTrack = new AliESDtrack(aAlphaTrack);
      pionESDTrack = new AliESDtrack(aPionTrack);
      AliExternalTrackParam pionTrk(*pionESDTrack);
      AliExternalTrackParam alphaTrk(*alphaESDTrack);

      Double_t xn, xp;
      Double_t dca = 0.;
      dca = pionESDTrack->GetDCA(alphaESDTrack, magField, xn,
                                 xp); //! dca (Neg to Pos)

      // CORRECTION from AliV0Vertex
      Bool_t corrected = kFALSE;
      if ((pionTrk.GetX() > 3.) && (xn < 3.)) {
        // correct for the beam pipe material
        corrected = kTRUE;
      }
      if ((alphaTrk.GetX() > 3.) && (xp < 3.)) {
        // correct for the beam pipe material
        corrected = kTRUE;
      }
      if (corrected) {
        dca = pionTrk.GetDCA(&alphaTrk, magField, xn, xp);
      }

      pionTrk.PropagateTo(xn, magField);
      alphaTrk.PropagateTo(xp, magField);

      AliESDv0 vertex(pionTrk, jTrack, alphaTrk, iTrack);
      // Float_t CosPointingAngle = std::abs(vertex.GetV0CosineOfPointingAngle(
      //     pv[0], pv[1], pv[2])); // PointingAngle, disable to use better way

      vertex.SetDcaV0Daughters(dca);
      // vertex.SetV0CosineOfPointingAngle(CosPointingAngle); // will be updated
      // below

      double sv[3];
      sv[0] = vertex.Xv();
      sv[1] = vertex.Yv();
      sv[2] = vertex.Zv();

      // Reconstruction
      LVector_t alphaVector, piVector, hyperVector;
      double alphaP[3];
      alphaTrack->GetPxPyPz(alphaP);
      double piP[3];
      pionTrack->GetPxPyPz(piP);

      alphaVector.SetCoordinates(alphaP[0] * 2, alphaP[1] * 2, alphaP[2] * 2,
                                 AliPID::ParticleMass(AliPID::kAlpha));
      piVector.SetCoordinates(piP[0], piP[1], piP[2],
                              AliPID::ParticleMass(AliPID::kPion));
      hyperVector = piVector + alphaVector;
      if (hyperVector.mass() > fMassRange[1] ||
          hyperVector.mass() < fMassRange[0]) {
        continue;
      }
      double deltaPos[3]{sv[0] - pv[0], sv[1] - pv[1], sv[2] - pv[2]};
      double CosPointingAngle =
          (deltaPos[0] * hyperVector.px() + deltaPos[1] * hyperVector.py() +
           deltaPos[2] * hyperVector.pz()) /
          std::sqrt(hyperVector.P2() *
                    (Sq(deltaPos[0]) + Sq(deltaPos[1]) + Sq(deltaPos[2])));

      vertex.SetV0CosineOfPointingAngle(CosPointingAngle);

      // Done
      const float beta = AliAnalysisTaskNucleiYield::HasTOF(alphaTrack, fPID);
      const int hasTOF = beta > 1.e-24 ? 1 : 0;
      double nsigmaTOFkAlpha =
          (hasTOF) ? fPID->NumberOfSigmasTOF(alphaTrack, AliPID::kAlpha) : -999;
      fRecHyper->pt = hyperVector.pt();
      fRecHyper->m = hyperVector.mass();
      fRecHyper->V0CosPA = CosPointingAngle;
      fRecHyper->Rapidity = Eta2y(fRecHyper->pt, kHyperMass, hyperVector.eta());

      fRecHyper->V0radius = TMath::Hypot(sv[0] - pv[0], sv[1] - pv[1]);
      fRecHyper->Lrec = TMath::Sqrt(TMath::Power(sv[0] - pv[0], 2) +
                                    TMath::Power(sv[1] - pv[1], 2) +
                                    TMath::Power(sv[2] - pv[2], 2));
      fRecHyper->ct = fRecHyper->Lrec * kHyperMass / hyperVector.P();
      Float_t b[2], bCov[3]; // Float due to the function input
      alphaTrack->GetImpactParameters(b, bCov);
      fRecHyper->alphaProngPvDCAXY = b[0];
      fRecHyper->alphaProngPvDCA = b[1];
      pionTrack->GetImpactParameters(b, bCov);
      fRecHyper->PiProngPvDCAXY = b[0];
      fRecHyper->PiProngPvDCA = b[1];
      fRecHyper->ProngsDCA = dca;
      fRecHyper->TPCmomalpha = alphaTrack->GetTPCmomentum();
      fRecHyper->TPCsignalalpha = alphaTrack->GetTPCsignal();
      fRecHyper->TOFnSigmaalpha = nsigmaTOFkAlpha;
      fRecHyper->NitsClustersalpha = alphaTrack->GetITSNcls();
      fRecHyper->TPCnSigmaPi =
          fPID->NumberOfSigmasTPC(pionTrack, AliPID::kPion);
      fRecHyper->TPCnSigmaalpha =
          fPID->NumberOfSigmasTPC(alphaTrack, AliPID::kAlpha);
      fRecHyper->NpidClustersPion = pionTrack->GetTPCsignalN();
      fRecHyper->NpidClustersalpha = alphaTrack->GetTPCsignalN();

      if (hyperLabel != -1) {
        if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(),
                      hyperLabel) != checkedHyperLabel.end()) {
          fGenHyper.isDuplicated = true;
        } else {
          fGenHyper.isDuplicated = false;
          checkedHyperLabel.push_back(hyperLabel);
        }
      }
      fTree->Fill();
    }
  }
  if (fMC) {
    fGenHyper.isReconstructed = false;
    // loop on generated
    for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT) {
      auto track = (AliAODMCParticle *)fMCEvent->GetTrack(iT);
      int pdg = std::abs(track->GetPdgCode());
      if (pdg != kHyperPdg) {
        continue;
      }
      if (std::find(checkedHyperLabel.begin(), checkedHyperLabel.end(), iT) !=
          checkedHyperLabel.end()) {
        continue;
      }

      if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(
              iT, header, MCTrackArray)) {
        continue;
      }
      fGenHyper.ptMC = track->Pt();
      fGenHyper.etaMC = track->Eta();
      fGenHyper.yMC = track->Y();
      fGenHyper.pdg = track->GetPdgCode();
      double ov[3], dv[3];
      track->XvYvZv(ov);
      bool otherDecayChannel{true};
      for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast();
           iD++) {
        auto daugh = (AliAODMCParticle *)fMCEvent->GetTrack(iD);
        if (!daugh) {
          continue;
        }
        if (std::abs(daugh->GetPdgCode()) ==
            AliPID::ParticleCode(AliPID::kAlpha)) {
          otherDecayChannel = false;
          daugh->XvYvZv(dv);
          break;
        }
      }
      if (otherDecayChannel)
        continue;
      fGenHyper.ctMC =
          std::sqrt(Sq(ov[0] - dv[0]) + Sq(ov[1] - dv[1]) + Sq(ov[2] - dv[2])) *
          track->M() / track->P();
      fTree->Fill();
    }
  }

  PostAllData();
}

AliAnalysisTaskAlphaPiAODtrk *
AliAnalysisTaskAlphaPiAODtrk::AddTask(bool isMC, TString tskname,
                                      TString suffix) {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHe4piAOD", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHe4piAOD", "This task requires an input event handler");
    return nullptr;
  }

  tskname.Append(suffix.Data());
  AliAnalysisTaskAlphaPiAODtrk *task =
      new AliAnalysisTaskAlphaPiAODtrk(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_output", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(
      Form("%s_treeHyper", tskname.Data()), TTree::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  return task;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskAlphaPiAODtrk::Eta2y(float pt, float m, float eta) const {
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}

void AliAnalysisTaskAlphaPiAODtrk::PostAllData() {
  PostData(1, fHistos->GetListOfHistograms());
  PostData(2, fTree);
}

void AliAnalysisTaskAlphaPiAODtrk::SetCustomBetheBloch(float res,
                                                       const float *bethe) {
  fCustomResolution = res;
  std::copy(bethe, bethe + 5, fCustomBethe);
}

double AliAnalysisTaskAlphaPiAODtrk::customNsigma(double mom, double sig) {
  // const float bg = mom / AliPID::ParticleMass(AliPID::kAlpha);
  // const float *p = fCustomBethe;
  // const float expS =
  //     AliExternalTrackParam::BetheBlochAleph(bg, p[0], p[1], p[2], p[3],
  //     p[4]);
  // return (sig - expS) / (fCustomResolution * expS);
  const int chargeAlpha = 2;
  const float *p = fCustomBethe;
  Double_t expected =
      chargeAlpha * chargeAlpha *
      AliExternalTrackParam::BetheBlochAleph(
          chargeAlpha * mom / AliPID::ParticleMass(AliPID::kAlpha), p[0], p[1],
          p[2], p[3], p[4]);
  Double_t sigma = expected * fCustomResolution;
  if (TMath::IsNaN(expected))
    return -999;
  return (sig - expected) / sigma;
}
