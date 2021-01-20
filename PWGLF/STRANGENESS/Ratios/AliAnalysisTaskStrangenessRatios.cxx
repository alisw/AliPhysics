#include "AliAnalysisTaskStrangenessRatios.h"

#include "AliAnalysisDataContainer.h"

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

constexpr double kLambdaMass{1.115683};
constexpr double kXiMass{1.32171};
constexpr double kOmegaMass{1.67245};

}

/// Standard and default constructor of the class.
///
/// \param taskname Name of the task
/// \param partname Name of the analysed particle
///
AliAnalysisTaskStrangenessRatios::AliAnalysisTaskStrangenessRatios(bool isMC, TString taskname) : AliAnalysisTaskSE(taskname.Data()),
                                                                                                  fEventCut{false},
                                                                                                  fList{nullptr},
                                                                                                  fTreeXi{nullptr},
                                                                                                  fTreeOmega{nullptr},
                                                                                                  fPID{nullptr},
                                                                                                  fMC{isMC},
                                                                                                  fCutRadiusXi{1.2},
                                                                                                  fCutRadiusOmega{1.0},
                                                                                                  fCutRadiusV0{3.0},
                                                                                                  fCutDCABachToPV{0.1},
                                                                                                  fCutDCAV0toPV{0.1},
                                                                                                  fCutDCAV0piToPV{0.2},
                                                                                                  fCutDCAV0prToPV{0.2},
                                                                                                  fCutDCAV0tracks{1.0},
                                                                                                  fCutDCABachToV0Xi{1.0},
                                                                                                  fCutDCABachToV0Omega{0.6},
                                                                                                  fCutCosPA{0.95},
                                                                                                  fCutCosPAV0{0.95},
                                                                                                  fCutV0MassWindow{0.005},
                                                                                                  fCutY{0.5},
                                                                                                  fCutYDaught{0.8},
                                                                                                  fCutNsigmaTPC{4.0},
                                                                                                  fCutCtXi{15},
                                                                                                  fCutCtOmega{12},
                                                                                                  fCutCtV0{30},
                                                                                                  fCutCompetingMass{0.008},
                                                                                                  fCutTPCclu{70}
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
    return;
  }

  auto pvObj = fEventCut.GetPrimaryVertex();
  double pv[3];
  pvObj->GetXYZ(pv);

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *handl = (AliInputEventHandler *)mgr->GetInputEventHandler();
  fPID = handl->GetPIDResponse();

  fRecCascade->centrality = fEventCut.GetCentrality();

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
    // fCasc_PosTOFBunchCrossing = pTrackCasc->GetTOFBunchCrossing(lMagField);
    // fCasc_NegTOFBunchCrossing = nTrackCasc->GetTOFBunchCrossing(lMagField);
    // fCasc_BacTOFBunchCrossing = bTrackCasc->GetTOFBunchCrossing(lMagField);

    //track status: ( fCasc_NegTrackStatus & AliESDtrack::kITSrefit ) is the codition to check kITSrefit
    fRecCascade->hasITSrefit = (nTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (pTrackCasc->GetStatus() & AliVTrack::kITSrefit) || (bTrackCasc->GetStatus() & AliVTrack::kITSrefit);

    fRecCascade->V0invMassDelta = ((fRecCascade->matter) ? casc->MassLambda() : casc->MassAntiLambda()) - kLambdaMass;

    //transverse momentum and eta
    fRecCascade->pt = std::sqrt(casc->Pt2Xi());
    fRecCascade->pt = std::sqrt(casc->Eta());

    //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
    fRecCascade->bachBarCosPA = casc->BachBaryonCosPA();
    double p = std::sqrt(casc->Ptot2Xi());
    fRecCascade->eta = 0.5 * std::log((p + casc->MomXiZ()) / (p - casc->MomXiZ() + 1.e-16));

    //distance over total momentum
    double lOverP = std::sqrt((Sq(vtxCasc[0] - pv[0]) + Sq(vtxCasc[1] - pv[1]) + Sq(vtxCasc[2] - pv[2])) / (casc->Ptot2Xi() + 1e-10));

    bool isTopolPassed = false;
    if (std::abs(casc->MassOmega() - kOmegaMass) * 1000 < 30) {
      fRecCascade->mass = casc->MassOmega();
      fRecCascade->ct = lOverP * kOmegaMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
      fRecCascade->competingMass = std::abs(casc->MassXi() - kXiMass);
      isTopolPassed = IsTopolSelected(0);
      if(!fMC && !isTopolPassed) continue;
      else if(fMC && isTopolPassed) fGenCascade.isReconstructed = true;
      fTreeOmega->Fill();
    } else if (std::abs(casc->MassXi() - kXiMass) * 1000 < 30) {
      fRecCascade->mass = casc->MassXi();
      fRecCascade->ct = lOverP * kXiMass;
      fRecCascade->tpcNsigmaBach = fPID->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      fRecCascade->competingMass = std::abs(casc->MassOmega() - kOmegaMass);
      isTopolPassed = IsTopolSelected(1);
      if(!fMC && !isTopolPassed) continue;
      else if(fMC && isTopolPassed) fGenCascade.isReconstructed = true;
      fTreeXi->Fill();
    }

    //MC association
    if (fMC)
    {
      fGenCascade.isReconstructed = true;
      // pdgPosDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->PdgCode();
      // pdgNegDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->PdgCode();
      // pdgBachelor = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->PdgCode();
      // int labMothPosDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->GetMother();
      // int labMothNegDaught = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->GetMother();
      // int labMothV0 = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->GetMother();
      // int labMothBach = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->GetMother();
      // pdgV0 = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))->PdgCode();
      // pdgCasc = ((AliAODMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothBach)))->PdgCode();
      // physprim = ((AliMCParticle *)lMCev->GetTrack((int)TMath::Abs(labMothBach)))->IsPhysicalPrimary();
      // if (fisMCassoc && (pdgPosDaught != 2212 || pdgNegDaught != -211 || pdgBachelor != -211 || labMothPosDaught != labMothNegDaught || pdgV0 != 3122 || labMothV0 != labMothBach || pdgCasc != 3312))
      //   assFlag[kxim] = kFALSE;
      // if (fisMCassoc && (pdgPosDaught != 211 || pdgNegDaught != -2212 || pdgBachelor != 211 || labMothPosDaught != labMothNegDaught || pdgV0 != -3122 || labMothV0 != labMothBach || pdgCasc != -3312))
      //   assFlag[kxip] = kFALSE;
      // if (fisMCassoc && (pdgPosDaught != 2212 || pdgNegDaught != -211 || pdgBachelor != -321 || labMothPosDaught != labMothNegDaught || pdgV0 != 3122 || labMothV0 != labMothBach || pdgCasc != 3334))
      //   assFlag[komm] = kFALSE;
      // if (fisMCassoc && (pdgPosDaught != 211 || pdgNegDaught != -2212 || pdgBachelor != 321 || labMothPosDaught != labMothNegDaught || pdgV0 != -3122 || labMothV0 != labMothBach || pdgCasc != -3334))
      //   assFlag[komp] = kFALSE;
    }
  }

  //  Post output data.
  PostData(1, fList);
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
bool AliAnalysisTaskStrangenessRatios::IsTopolSelected(bool isXi)
{
  if( fRecCascade->radius > (isXi) ? fCutRadiusXi : fCutRadiusOmega &&
      fRecCascade->radiusV0 > fCutRadiusV0 &&
      fRecCascade->dcaBachPV > fCutDCABachToPV &&
      fRecCascade->dcaV0PV > fCutDCAV0toPV &&
      fRecCascade->dcaV0piPV > fCutDCAV0piToPV &&
      fRecCascade->dcaV0prPV > fCutDCAV0prToPV &&
      fRecCascade->dcaV0tracks > fCutDCAV0tracks &&
      fRecCascade->dcaBachV0 > (isXi) ? fCutDCABachToV0Xi : fCutDCABachToV0Omega &&
      fRecCascade->cosPA > fCutCosPA &&
      fRecCascade->cosPAV0 > fCutCosPAV0 &&
      fRecCascade->dcaV0prPV > fCutDCAV0prToPV &&
      std::abs(Eta2y(fRecCascade->pt, fRecCascade->mass, fRecCascade->eta)) < fCutY &&
      fRecCascade->tpcNsigmaBach < fCutNsigmaTPC &&
      fRecCascade->tpcNsigmaV0Pr < fCutNsigmaTPC &&
      fRecCascade->tpcNsigmaV0Pi < fCutNsigmaTPC &&
      fRecCascade->ct < (isXi) ? fCutCtXi : fCutCtOmega &&
      fRecCascade->competingMass > fCutCompetingMass &&
      fRecCascade->tpcClBach > fCutTPCclu &&
      fRecCascade->tpcClPi > fCutTPCclu &&
      fRecCascade->tpcClPr > fCutTPCclu &&
  )
  {
    return true;
  }
  return false;
}

//
//____________________________________________________________________________________________
float AliAnalysisTaskStrangenessRatios::Eta2y(float pt, float m, float eta) const
{
  // convert eta to y
  double mt = std::sqrt(m * m + pt * pt);
  return std::asinh(pt / mt * TMath::SinH(eta));
}
