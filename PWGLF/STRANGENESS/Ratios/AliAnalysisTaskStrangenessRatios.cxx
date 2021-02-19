#include "AliAnalysisTaskStrangenessRatios.h"

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

///\cond CLASSIMP
ClassImp(AliAnalysisTaskStrangenessRatios);
///\endcond

namespace {

double Sq(double x) {
  return x * x;
}

constexpr int kLambdaPdg{3122};
constexpr double kLambdaMass{1.115683};
constexpr int kXiPdg{3312};
constexpr double kXiMass{1.32171};
constexpr int kOmegaPdg{3334};
constexpr double kOmegaMass{1.67245};
constexpr double kcTauXi{4.91359839};
constexpr double kcTauOmega{2.46129608};
constexpr double kcTau[2]{kcTauXi, kcTauOmega};

}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskStrangenessRatios::AliAnalysisTaskStrangenessRatios(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                                  fEventCut{false},
                                                                                                  fMC{isMC}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

/// Standard destructor
///
AliAnalysisTaskStrangenessRatios::~AliAnalysisTaskStrangenessRatios()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    return;
  if (fList)
    delete fList;
  if (fTreeXi)
    delete fTreeXi;
  if (fTreeOmega)
    delete fTreeOmega;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskStrangenessRatios::UserCreateOutputObjects()
{

  fList = new TList();
  fList->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(fList);

  PostData(1, fList);

  fRecCascade = fMC ? &fGenCascade : new MiniCascade;

  OpenFile(2);
  fTreeXi = new TTree("XiTree", "Xi Tree");
  OpenFile(3);
  fTreeOmega = new TTree("OmegaTree", "Omega Tree");

  if (fMC)
  {
    fTreeXi->Branch("MiniCascadeMC", &fGenCascade);
    fTreeOmega->Branch("MiniCascadeMC", &fGenCascade);
    fMCEvent = MCEvent();
  }
  else
  {
    fTreeXi->Branch("MiniCascade", fRecCascade);
    fTreeOmega->Branch("MiniCascade", fRecCascade);
  }

  PostData(2, fTreeXi);
  PostData(3, fTreeOmega);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskStrangenessRatios::UserExec(Option_t *)
{
  AliAODEvent *ev = (AliAODEvent*)InputEvent();
  if (!fEventCut.AcceptEvent(ev))
  {
    PostData(1, fList);
    PostData(2, fTreeXi);
    PostData(3, fTreeOmega);
    return;
  }

  double bField{ev->GetMagneticField()};
  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3];
  pvObj->GetXYZ(pv);

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  fRecCascade->centrality = fEventCut.GetCentrality();

  std::vector<int> checkedLabel;
  fGenCascade.isReconstructed = true;
  for (int iCasc = 0; iCasc < ev->GetNumberOfCascades(); iCasc++)
  {
    AliAODcascade *casc = ev->GetCascade(iCasc);
    if (!casc)
      continue;

    //cascade and V0 2D radii
    double vtxCasc[3]{casc->DecayVertexXiX(), casc->DecayVertexXiY(), casc->DecayVertexXiZ()};
    fRecCascade->radius = std::hypot(vtxCasc[0], vtxCasc[1]);
    fRecCascade->radiusV0 = casc->RadiusSecVtx();

    //get daughter tracks (positive, negative and bachelor)
    AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(0));
    AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDaughter(1));
    AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack *>(casc->GetDecayVertexXi()->GetDaughter(0));
    if (!pTrackCasc || !nTrackCasc || !bTrackCasc)
    {
      AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n");
      continue;
    }

    if (!(pTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(nTrackCasc->GetStatus() & AliVTrack::kTPCrefit) || !(bTrackCasc->GetStatus() & AliVTrack::kTPCrefit)) {
      continue;
    }

    if (pTrackCasc->GetTPCsignalN() < 50 || nTrackCasc->GetTPCsignalN() < 50 || bTrackCasc->GetTPCsignalN() < 50) {
      continue;
    }

    if (std::abs(pTrackCasc->Eta()) > 0.9 || std::abs(nTrackCasc->Eta()) > 0.9 || std::abs(bTrackCasc->Eta()) > 0.9) {
      continue;
    }

    if (fMC) {
      fGenCascade.pdg = 0;
      auto posPart = (AliAODMCParticle*)fMCEvent->GetTrack(std::abs(pTrackCasc->GetLabel()));
      auto negPart = (AliAODMCParticle*)fMCEvent->GetTrack(std::abs(nTrackCasc->GetLabel()));
      auto bacPart = (AliAODMCParticle*)fMCEvent->GetTrack(std::abs(bTrackCasc->GetLabel()));
      // Check lambda
      int labMothPos = posPart->GetMother();
      int labMothNeg = negPart->GetMother();
      int labMothBac = bacPart->GetMother();
      auto lambda = (AliAODMCParticle*)fMCEvent->GetTrack(labMothNeg);
      if (lambda && labMothNeg == labMothPos && std::abs(lambda->GetPdgCode()) == kLambdaPdg) {
        int labMothLam = lambda->GetMother();
        auto cascade = (AliAODMCParticle*)fMCEvent->GetTrack(labMothBac);
        if (!cascade) {
          continue;
        }
        int pdgCascade = std::abs(cascade->GetPdgCode());
        if (labMothLam == labMothBac && (pdgCascade == kXiPdg || pdgCascade == kOmegaPdg)) {
          fGenCascade.pdg = cascade->GetPdgCode();
          fGenCascade.ptMC = cascade->Pt();
          fGenCascade.etaMC = cascade->Eta();
          fGenCascade.yMC = cascade->Y();
          double pv[3], sv[3];
          cascade->XvYvZv(pv);
          bacPart->XvYvZv(sv);
          fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * cascade->M() / cascade->P();
          checkedLabel.push_back(labMothBac);
        }
      }
      if (fOnlyTrueCandidates && fGenCascade.pdg == 0)
        continue;
    }

    fRecCascade->matter = casc->AlphaV0() > 0;

    fRecCascade->tpcNsigmaV0Pi = fPID->NumberOfSigmasTPC(fRecCascade->matter ? nTrackCasc : pTrackCasc, AliPID::kPion);
    fRecCascade->tpcNsigmaV0Pr = fPID->NumberOfSigmasTPC(fRecCascade->matter ? pTrackCasc : nTrackCasc, AliPID::kProton);

    fRecCascade->tpcClBach = bTrackCasc->GetTPCsignalN();
    fRecCascade->tpcClV0Pi = (fRecCascade->matter ? nTrackCasc : pTrackCasc)->GetTPCsignalN();
    fRecCascade->tpcClV0Pr = (fRecCascade->matter ? pTrackCasc : nTrackCasc)->GetTPCsignalN();

    //DCA info
    fRecCascade->dcaBachV0 = casc->DcaXiDaughters();
    fRecCascade->dcaBachPV = casc->DcaBachToPrimVertex();
    fRecCascade->dcaV0prPV = fRecCascade->matter ? casc->DcaPosToPrimVertex() : casc->DcaNegToPrimVertex();
    fRecCascade->dcaV0piPV = fRecCascade->matter ? casc->DcaNegToPrimVertex() : casc->DcaPosToPrimVertex();
    fRecCascade->dcaV0tracks = casc->DcaV0Daughters();
    fRecCascade->dcaV0PV = casc->DcaV0ToPrimVertex();

    //cascade and V0 cosine of pointing angle
    fRecCascade->cosPA = casc->CosPointingAngleXi(pv[0], pv[1], pv[2]);
    fRecCascade->cosPAV0 = casc->CosPointingAngle(pv);

    //TOF matching
    fRecCascade->hasTOFhit = !pTrackCasc->GetTOFBunchCrossing(bField) || !nTrackCasc->GetTOFBunchCrossing(bField) || !bTrackCasc->GetTOFBunchCrossing(bField);

    //track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
    fRecCascade->hasITSrefit = (nTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (pTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (bTrackCasc->GetStatus() & AliVTrack::kITSrefit);

    fRecCascade->V0invMassDelta = ((fRecCascade->matter) ? casc->MassLambda() : casc->MassAntiLambda()) - kLambdaMass;

    //transverse momentum and eta
    fRecCascade->pt = std::sqrt(casc->Pt2Xi());

    //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
    fRecCascade->bachBarCosPA = casc->BachBaryonCosPA();
    double p = std::sqrt(casc->Ptot2Xi());
    fRecCascade->eta = 0.5 * std::log((p + casc->MomXiZ()) / (p - casc->MomXiZ() + 1.e-16));

    //distance over total momentum
    double lOverP = std::sqrt((Sq(vtxCasc[0] - pv[0]) + Sq(vtxCasc[1] - pv[1]) + Sq(vtxCasc[2] - pv[2])) / (casc->Ptot2Xi() + 1e-10));
    double ctLambda = std::sqrt(Sq(vtxCasc[0] - casc->GetSecVtxX()) + Sq(vtxCasc[1] - casc->GetSecVtxY()) + Sq(vtxCasc[2] - casc->GetSecVtxZ())) / (casc->P() + 1e-10);

    if (std::abs(casc->MassOmega() - kOmegaMass) * 1000 < 30) {
      fRecCascade->mass = casc->MassOmega();
      fRecCascade->ct = lOverP * kOmegaMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
      fRecCascade->competingMass = std::abs(casc->MassXi() - kXiMass);
      if (IsTopolSelected(true)) {
        fTreeOmega->Fill();
      }
    }
    if (std::abs(casc->MassXi() - kXiMass) * 1000 < 30) {
      fRecCascade->mass = casc->MassXi();
      fRecCascade->ct = lOverP * kXiMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      fRecCascade->competingMass = std::abs(casc->MassOmega() - kOmegaMass);
      if (IsTopolSelected(false)) {
        fTreeXi->Fill();
      }
    }
  }

  if (fMC) {
    fGenCascade.isReconstructed = false;
    for (int iT{0}; iT < fMCEvent->GetNumberOfTracks(); ++iT) {
      auto track = (AliAODMCParticle*)fMCEvent->GetTrack(iT);
      int pdg = std::abs(track->GetPdgCode());
      if (pdg != kXiPdg && pdg != kOmegaPdg) {
        continue;
      }
      if (std::find(checkedLabel.begin(), checkedLabel.end(), iT) != checkedLabel.end()) {
        continue;
      }
      fGenCascade.ptMC = track->Pt();
      fGenCascade.etaMC = track->Eta();
      fGenCascade.yMC = track->Y();
      fGenCascade.pdg = track->GetPdgCode();
      double pv[3], sv[3];
      track->XvYvZv(pv);
      for (int iD = track->GetDaughterFirst(); iD <= track->GetDaughterLast(); iD++) {
        auto daugh = (AliAODMCParticle*)fMCEvent->GetTrack(iD);
        if (!daugh) {
          continue;
        }
        if (std::abs(daugh->GetPdgCode()) == kLambdaPdg) {
          daugh->XvYvZv(sv);
          break;
        }
      }
      fGenCascade.ctMC = std::sqrt(Sq(pv[0] - sv[0]) + Sq(pv[1] - sv[1]) + Sq(pv[2] - sv[2])) * track->M() / track->P();
      (pdg == kXiPdg ? fTreeXi : fTreeOmega)->Fill();
    }
  }

  //  Post output data.
  PostData(1, fList);
  PostData(2, fTreeXi);
  PostData(3, fTreeOmega);
}

AliAnalysisTaskStrangenessRatios *AliAnalysisTaskStrangenessRatios::AddTask(bool isMC, TString tskname, TString suffix)
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
  AliAnalysisTaskStrangenessRatios *task = new AliAnalysisTaskStrangenessRatios(isMC, tskname.Data());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(
      Form("%s_summary", tskname.Data()), TList::Class(),
      AliAnalysisManager::kOutputContainer, "AnalysisResults.root");

  AliAnalysisDataContainer *coutput2 =
      mgr->CreateContainer(Form("%s_treeXi", tskname.Data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "AnalysisResults.root");
  coutput2->SetSpecialOutput();

  AliAnalysisDataContainer *coutput3 =
      mgr->CreateContainer(Form("%s_treeOmega", tskname.Data()), TTree::Class(),
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
bool AliAnalysisTaskStrangenessRatios::IsTopolSelected(bool isOmega)
{
  return fRecCascade->radius > fCutRadius[isOmega] &&
      fRecCascade->radiusV0 > fCutRadiusV0 &&
      fRecCascade->dcaBachPV > fCutDCABachToPV &&
      fRecCascade->dcaV0PV > fCutDCAV0toPV &&
      fRecCascade->dcaV0piPV > fCutDCAV0piToPV &&
      fRecCascade->dcaV0prPV > fCutDCAV0prToPV &&
      fRecCascade->dcaV0tracks < fCutDCAV0tracks &&
      fRecCascade->dcaBachV0 < fCutDCABachToV0[isOmega] &&
      fRecCascade->cosPA > fCutCosPA &&
      fRecCascade->cosPAV0 > fCutCosPAV0 &&
      fRecCascade->dcaV0prPV > fCutDCAV0prToPV &&
      std::abs(Eta2y(fRecCascade->pt, isOmega ? kOmegaMass : kXiMass, fRecCascade->eta)) < fCutY &&
      fRecCascade->tpcNsigmaBach < fCutNsigmaTPC &&
      fRecCascade->tpcNsigmaV0Pr < fCutNsigmaTPC &&
      fRecCascade->tpcNsigmaV0Pi < fCutNsigmaTPC &&
      fRecCascade->ct < fCutCt * kcTau[isOmega] &&
      fRecCascade->competingMass > fCutCompetingMass &&
      fRecCascade->tpcClBach > fCutTPCclu &&
      fRecCascade->tpcClV0Pi > fCutTPCclu &&
      fRecCascade->tpcClV0Pr > fCutTPCclu;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskStrangenessRatios::Eta2y(float pt, float m, float eta) const
{
  return std::asinh(pt / std::hypot(m, pt) * std::sinh(eta));
}
