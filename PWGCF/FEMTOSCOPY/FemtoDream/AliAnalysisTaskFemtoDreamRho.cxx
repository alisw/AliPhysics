/*
 * AliAnalysisTaskFemtoDreamRho.cxx
 *
 *  Created on: 19 Jul 2023
 *      Author: M. Korwieser
 */

#include "AliAnalysisTaskFemtoDreamRho.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisManager.h"
#include "TVector3.h"
#include "TLorentzVector.h"

ClassImp(AliAnalysisTaskFemtoDreamRho)
    AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho()
    : AliAnalysisTaskSE(),
      fTrigger(AliVEvent::kINT7),
      fIsMC(false),
      fDoMcTruth(false),
      fDoCleaning(false),
      fDoAncestors(false),
      fDoProjections(false),
      frhoPtThreshold(0.),
      fIsSameCharge(false),
      fIsMCTrueRhoCombBkrg(false),
      fIsMCcheckedCombs(false),
      fOutput(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fTrackneg(nullptr),
      fRhoParticle(nullptr),
      fEventCuts(nullptr),
      fPosPionCuts(nullptr),
      fNegPionCuts(nullptr),
      fPosPionMinvCuts(nullptr),
      fNegPionMinvCuts(nullptr),
      fRhoCuts(nullptr),
      fPosProtonCuts(nullptr),
      fNegProtonCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fArmenterosRhoTrue(nullptr),
      fArmenterosRhoTrue_Reconstr(nullptr),
      fArmenterosNoCommonMother_Pos(nullptr),
      fArmenterosNoCommonMother_Neg(nullptr),
      fArmenterosNoCommonMother_qtDaughBoth(nullptr),
      fArmenterosNoCommonMother_alphaDaughBoth(nullptr),
      fArmenterosNoRhoTrue_Reconstr_Pos(nullptr),
      fArmenterosNoRhoTrue_Reconstr_Neg(nullptr),
      fArmenterosNoRhoTrue_Reconstr_qtDaughBoth(nullptr),
      fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth(nullptr),
      fArmenterosRhoTrue_Reconstr_qtDaughBoth(nullptr),
      fArmenterosRhoTrue_Reconstr_alphaDaughBoth(nullptr),
      fHist2D_massVSpt_RhoTrue(nullptr),
      fHist2D_massVSpt_RhoCandidateCommon(nullptr),
      fHist2D_massVSpt_RhoCandidateUncommon(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateUncommonFullInvM(nullptr),
      fHist1D_pt_RhoTrue(nullptr),
      fHist2D_pt1VSpt2_RhoTrue(nullptr),
      fHist2D_pTvsmT_noPions(nullptr),
      fHist2D_pTvsmT_noPrims(nullptr),
      fHist2D_pTvsmT_noCommonMother(nullptr),
      fHist2D_pTvsmT_noRho(nullptr),
      fHist2D_pTvsmT_noRho_MC(nullptr),
      fHist2D_PDGvsmT_noRho_MC(nullptr),
      fHist2D_pTvsmT_isRho(nullptr),
      fHist2D_pTvsmT_isRho_MC(nullptr),
      fHist2D_PDGvsMInv_CommonAncestorResonances(nullptr),
      fGTI(0),
      fTrackBufferSize(0)
{
}

AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho(const char *name,
                                                           bool isMC, bool doMcTruth, bool doCleaning, bool doAncestors, bool doProjector, float rhoPtThreshold, bool isSameCharge, bool isMCTrueRhoCombBkrg, bool isMCcheckedCombs)
    : AliAnalysisTaskSE(name),
      fTrigger(AliVEvent::kINT7),
      fIsMC(isMC),
      fDoMcTruth(doMcTruth),
      fDoCleaning(doCleaning),
      fDoAncestors(doAncestors),
      fDoProjections(doProjector),
      frhoPtThreshold(rhoPtThreshold),
      fIsSameCharge(isSameCharge),
      fIsMCTrueRhoCombBkrg(isMCTrueRhoCombBkrg),
      fIsMCcheckedCombs(isMCcheckedCombs),
      fOutput(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fTrackneg(nullptr),
      fRhoParticle(nullptr),
      fEventCuts(nullptr),
      fPosPionCuts(nullptr),
      fNegPionCuts(nullptr),
      fPosPionMinvCuts(nullptr),
      fNegPionMinvCuts(nullptr),
      fRhoCuts(nullptr),
      fPosProtonCuts(nullptr),
      fNegProtonCuts(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fArmenterosRhoTrue(nullptr),
      fArmenterosRhoTrue_Reconstr(nullptr),
      fArmenterosNoCommonMother_Pos(nullptr),
      fArmenterosNoCommonMother_Neg(nullptr),
      fArmenterosNoCommonMother_qtDaughBoth(nullptr),
      fArmenterosNoCommonMother_alphaDaughBoth(nullptr),
      fArmenterosNoRhoTrue_Reconstr_Pos(nullptr),
      fArmenterosNoRhoTrue_Reconstr_Neg(nullptr),
      fArmenterosNoRhoTrue_Reconstr_qtDaughBoth(nullptr),
      fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth(nullptr),
      fArmenterosRhoTrue_Reconstr_qtDaughBoth(nullptr),
      fArmenterosRhoTrue_Reconstr_alphaDaughBoth(nullptr),
      fHist2D_massVSpt_RhoTrue(nullptr),
      fHist2D_massVSpt_RhoCandidateCommon(nullptr),
      fHist2D_massVSpt_RhoCandidateUncommon(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances(nullptr),
      fHist2D_massVSpt_RhoCandidateUncommonFullInvM(nullptr),
      fHist1D_pt_RhoTrue(nullptr),
      fHist2D_pt1VSpt2_RhoTrue(nullptr),
      fHist2D_pTvsmT_noPions(nullptr),
      fHist2D_pTvsmT_noPrims(nullptr),
      fHist2D_pTvsmT_noCommonMother(nullptr),
      fHist2D_pTvsmT_noRho(nullptr),
      fHist2D_pTvsmT_noRho_MC(nullptr),
      fHist2D_PDGvsmT_noRho_MC(nullptr),
      fHist2D_pTvsmT_isRho(nullptr),
      fHist2D_pTvsmT_isRho_MC(nullptr),
      fHist2D_PDGvsMInv_CommonAncestorResonances(nullptr),
      fGTI(nullptr),
      fTrackBufferSize(2000)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskFemtoDreamRho::~AliAnalysisTaskFemtoDreamRho() {}

void AliAnalysisTaskFemtoDreamRho::UserCreateOutputObjects()
{
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, fTrigger); // Here the default pT used for the S_T selection can be changed
  fOutput->Add(fEvent->GetEvtCutList());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);
  // fTrackneg = new AliFemtoDreamTrack();
  // fTrackneg->SetUseMCInfo(fIsMC);

  fRhoParticle = new AliFemtoDreamv0();
  fRhoParticle->SetPDGCode(fRhoCuts->GetPDGv0());
  fRhoParticle->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterPos(
      fRhoCuts->GetPDGPosDaug()); // order +sign doesnt play a role
  fRhoParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterNeg(
      fRhoCuts->GetPDGNegDaug()); // only used for MC Matching
  fRhoParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack *[fTrackBufferSize];

  if (!fEventCuts)
  {
    AliFatal("Event Cuts not set!");
  }
  fEventCuts->InitQA();
  fOutput->Add(fEventCuts->GetHistList());

  if (!fPosProtonCuts)
  {
    AliFatal("Track Cuts for Proton not set!");
  }
  fPosProtonCuts->Init();
  fPosProtonCuts->SetName("Proton");
  fOutput->Add(fPosProtonCuts->GetQAHists());
  if (fPosProtonCuts->GetIsMonteCarlo())
  {
    fPosProtonCuts->SetMCName("MCProton");
    fOutput->Add(fPosProtonCuts->GetMCQAHists());
  }

  if (!fNegProtonCuts)
  {
    AliFatal("Track Cuts for AntiProton not set!");
  }
  fNegProtonCuts->Init();
  fNegProtonCuts->SetName("AntiProton");
  fOutput->Add(fNegProtonCuts->GetQAHists());
  if (fNegProtonCuts->GetIsMonteCarlo())
  {
    fNegProtonCuts->SetMCName("MCAntiProton");
    fOutput->Add(fNegProtonCuts->GetMCQAHists());
  }

  if (!fPosPionCuts)
  {
    AliFatal("Track Cuts for positive pion not set!");
  }
  fPosPionCuts->Init();
  fPosPionCuts->SetName("PosPion");
  fOutput->Add(fPosPionCuts->GetQAHists());
  if (fPosPionCuts->GetIsMonteCarlo())
  {
    fPosPionCuts->SetMCName("MCPosPion");
    fOutput->Add(fPosPionCuts->GetMCQAHists());
  }

  if (!fNegPionCuts)
  {
    AliFatal("Track Cuts for negative pion not set!");
  }
  fNegPionCuts->Init();
  fNegPionCuts->SetName("NegPion");
  fOutput->Add(fNegPionCuts->GetQAHists());
  if (fNegPionCuts->GetIsMonteCarlo())
  {
    fNegPionCuts->SetMCName("MCNegPion");
    fOutput->Add(fNegPionCuts->GetMCQAHists());
  }

  if (fDoProjections)
  {
    if (!fPosPionMinvCuts)
    {
      AliFatal("Track Cuts for positive pion Minv not set!");
    }
    fPosPionMinvCuts->Init();
    fPosPionMinvCuts->SetName("PosPionMinv");
    fOutput->Add(fPosPionMinvCuts->GetQAHists());
    if (fPosPionMinvCuts->GetIsMonteCarlo())
    {
      fPosPionMinvCuts->SetMCName("MCPosPionMinv");
      fOutput->Add(fPosPionMinvCuts->GetMCQAHists());
    }

    if (!fNegPionMinvCuts)
    {
      AliFatal("Track Cuts for negative pion Minv not set!");
    }
    fNegPionMinvCuts->Init();
    fNegPionMinvCuts->SetName("NegPionMinv");
    fOutput->Add(fNegPionMinvCuts->GetQAHists());
    if (fNegPionMinvCuts->GetIsMonteCarlo())
    {
      fNegPionMinvCuts->SetMCName("MCNegPionMinv");
      fOutput->Add(fNegPionMinvCuts->GetMCQAHists());
    }
  }

  if (!fRhoCuts)
  {
    AliFatal("Cuts for the Rho not set!");
  }
  fRhoCuts->Init();
  fRhoCuts->SetName("RhoCandidates");
  fOutput->Add(fRhoCuts->GetQAHists()); // add here a dummy cut object for the MC Truth (could also just be the acceptance selections)
  if (fRhoCuts->GetIsMonteCarlo())
  {
    fRhoCuts->SetMCName("MCRhoCandidates");
    fOutput->Add(fRhoCuts->GetMCQAHists());
  }

  fPairCleaner =
      new AliFemtoDreamPairCleaner(0, 0, fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());
  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());
  if (fDoAncestors && fIsMC)
  {
    // out << "Will not be accessed" << std::endl;

    // Create histograms for the Ancestor investigation of the RhoCandidates
    auto *fHistListRhoCandidatesMCAncestors = new TList();
    fHistListRhoCandidatesMCAncestors->SetName("RhoCandidatesMCAncestor");
    fHistListRhoCandidatesMCAncestors->SetOwner();
    // String
    TString NameIngaAnc = "RhoCandidatesMC";

    TString massPtHistName_Common = TString::Format("histInvariantMassPt_Common%s", NameIngaAnc.Data());
    TString massPtHistTitle_Common = TString::Format("%s Invariant Mass vs. pT Common", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommon = new TH2F(massPtHistName_Common, massPtHistTitle_Common, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommon->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommon->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_Uncommon = TString::Format("histInvariantMassPt_Uncommon%s", NameIngaAnc.Data());
    TString massPtHistTitle_Uncommon = TString::Format("%s Invariant Mass vs. pT Uncommon", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateUncommon = new TH2F(massPtHistName_Uncommon, massPtHistTitle_Uncommon, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateUncommon->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateUncommon->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM = TString::Format("histInvariantMassPt_Common_FullMinv%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM = TString::Format("%s Invariant Mass vs. pT Common", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM = new TH2F(massPtHistName_CommonFullInvM, massPtHistTitle_CommonFullInvM, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_NoResonances = TString::Format("histInvariantMassPt_Common_FullMinv_noResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_NoResonances = TString::Format("%s Invariant Mass vs. pT Common No Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances = new TH2F(massPtHistName_CommonFullInvM_NoResonances, massPtHistTitle_CommonFullInvM_NoResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_kShortResonances = TString::Format("histInvariantMassPt_Common_FullMinv_kShortResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_kShortResonances = TString::Format("%s Invariant Mass vs. pT Common KShort Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances = new TH2F(massPtHistName_CommonFullInvM_kShortResonances, massPtHistTitle_CommonFullInvM_kShortResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_rhoResonances = TString::Format("histInvariantMassPt_Common_FullMinv_rhoResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_rhoResonances = TString::Format("%s Invariant Mass vs. pT Common Rho Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances = new TH2F(massPtHistName_CommonFullInvM_rhoResonances, massPtHistTitle_CommonFullInvM_rhoResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_omegaResonances = TString::Format("histInvariantMassPt_Common_FullMinv_omegaResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_omegaResonances = TString::Format("%s Invariant Mass vs. pT Common omega Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances = new TH2F(massPtHistName_CommonFullInvM_omegaResonances, massPtHistTitle_CommonFullInvM_omegaResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_fzeroResonances = TString::Format("histInvariantMassPt_Common_FullMinv_fzeroResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_fzeroResonances = TString::Format("%s Invariant Mass vs. pT Common f0 Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances = new TH2F(massPtHistName_CommonFullInvM_fzeroResonances, massPtHistTitle_CommonFullInvM_fzeroResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_ftwoResonances = TString::Format("histInvariantMassPt_Common_FullMinv_ftwoResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_ftwoResonances = TString::Format("%s Invariant Mass vs. pT Common f2 Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances = new TH2F(massPtHistName_CommonFullInvM_ftwoResonances, massPtHistTitle_CommonFullInvM_ftwoResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_CommonFullInvM_otherResonances = TString::Format("histInvariantMassPt_Common_FullMinv_otherResonances%s", NameIngaAnc.Data());
    TString massPtHistTitle_CommonFullInvM_otherResonances = TString::Format("%s Invariant Mass vs. pT Common other Resonances", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances = new TH2F(massPtHistName_CommonFullInvM_otherResonances, massPtHistTitle_CommonFullInvM_otherResonances, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    TString massPtHistName_UncommonFullInvM = TString::Format("histInvariantMassPt_Uncommon_FullMinv%s", NameIngaAnc.Data());
    TString massPtHistTitle_UncommonFullInvM = TString::Format("%s Invariant Mass vs. pT Uncommon", NameIngaAnc.Data());
    fHist2D_massVSpt_RhoCandidateUncommonFullInvM = new TH2F(massPtHistName_UncommonFullInvM, massPtHistTitle_UncommonFullInvM, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoCandidateUncommonFullInvM->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoCandidateUncommonFullInvM->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    // Ancestor PDG resonance tracking
    fHist2D_PDGvsMInv_CommonAncestorResonances = new TH2F("fHist2D_PDGvsMInv_CommonAncestorResonances", "fHist2D_PDGvsMInv_CommonAncestorResonances", 5001, -2500, 2500, 500, 0.0, 5.0);
    fHist2D_PDGvsMInv_CommonAncestorResonances->GetXaxis()->SetTitle("PDG");
    fHist2D_PDGvsMInv_CommonAncestorResonances->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommon);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateUncommon);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_massVSpt_RhoCandidateUncommonFullInvM);
    fHistListRhoCandidatesMCAncestors->Add(fHist2D_PDGvsMInv_CommonAncestorResonances);

    fOutput->Add(fHistListRhoCandidatesMCAncestors); // Add the list to your output list
  }
  if (fIsMC)
  {
    // Create histograms for RhoMCTrue particles
    auto *fHistListRhoMCTrue = new TList();
    fHistListRhoMCTrue->SetName("RhoMCTrue");
    fHistListRhoMCTrue->SetOwner();
    // String
    TString NameIng = "RhoMCTrue";

    // Create TH1F for pT spectrum
    TString ptHistName = TString::Format("histPtSpectrum_%s", NameIng.Data());
    TString ptHistTitle = TString::Format("%s pT Spectrum", NameIng.Data());
    fHist1D_pt_RhoTrue = new TH1F(ptHistName, ptHistTitle, 1000, 0.0, 5.0);
    fHist1D_pt_RhoTrue->GetXaxis()->SetTitle("pT (GeV/c)");

    // Track pTdaugther1 vs pTdaugther2
    TString pt1pt2HistName = TString::Format("histPtSpectrum_%s", NameIng.Data());
    TString pt1pt2HistTitle = TString::Format("%s pT Spectrum", NameIng.Data());
    fHist2D_pt1VSpt2_RhoTrue = new TH2F(pt1pt2HistName, pt1pt2HistTitle, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pt1VSpt2_RhoTrue->GetXaxis()->SetTitle("pT1 (GeV/c)");
    fHist2D_pt1VSpt2_RhoTrue->GetYaxis()->SetTitle("pT2 (GeV/c)");

    // Create TH2F for invariant mass vs pT
    TString massPtHistName = TString::Format("histInvariantMassPt_%s", NameIng.Data());
    TString massPtHistTitle = TString::Format("%s Invariant Mass vs. pT", NameIng.Data());
    fHist2D_massVSpt_RhoTrue = new TH2F(massPtHistName, massPtHistTitle, 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_massVSpt_RhoTrue->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_massVSpt_RhoTrue->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    // Create histogram to track the armenteros values of MC true rhos
    fArmenterosRhoTrue = new TH2F("fArmenterosRhoTrue", "fArmenterosRhoTrue", 200, -1, 1, 400, 0, 0.8);
    fArmenterosRhoTrue->GetXaxis()->SetTitle("#alpha");
    fArmenterosRhoTrue->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos reconstructed
    fArmenterosRhoTrue_Reconstr = new TH2F("fArmenterosRhoTrue_Reconstr", "fArmenterosRhoTrue_Reconstr", 200, -1, 1, 400, 0, 0.8);
    fArmenterosRhoTrue_Reconstr->GetXaxis()->SetTitle("#alpha");
    fArmenterosRhoTrue_Reconstr->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC false rhos with posDaughter
    fArmenterosNoCommonMother_Pos = new TH2F("fArmenterosNoCommonMother_Pos", "fArmenterosNoCommonMother_Pos", 200, -1, 1, 400, 0, 0.8);
    fArmenterosNoCommonMother_Pos->GetXaxis()->SetTitle("#alpha");
    fArmenterosNoCommonMother_Pos->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoCommonMother_Neg = new TH2F("fArmenterosNoCommonMother_Neg", "fArmenterosNoCommonMother_Neg", 200, -1, 1, 400, 0, 0.8);
    fArmenterosNoCommonMother_Neg->GetXaxis()->SetTitle("#alpha");
    fArmenterosNoCommonMother_Neg->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoCommonMother_qtDaughBoth = new TH2F("fArmenterosNoCommonMother_qtDaughBoth", "fArmenterosNoCommonMother_qtDaughBoth", 100, 0, 0.8, 400, 0, 0.8);
    fArmenterosNoCommonMother_qtDaughBoth->GetXaxis()->SetTitle("#it{q}_{T}_pos (GeV/#it{c})");
    fArmenterosNoCommonMother_qtDaughBoth->GetYaxis()->SetTitle("#it{q}_{T}_neg (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoCommonMother_alphaDaughBoth = new TH2F("fArmenterosNoCommonMother_alphaDaughBoth", "fArmenterosNoCommonMother_alphaDaughBoth", 200, -1, 1, 200, -1, 1);
    fArmenterosNoCommonMother_alphaDaughBoth->GetXaxis()->SetTitle("#alpha_pos");
    fArmenterosNoCommonMother_alphaDaughBoth->GetYaxis()->SetTitle("#alpha_neg");

    // Create histogram to track the armenteros values of MC false rhos with posDaughter
    fArmenterosNoRhoTrue_Reconstr_Pos = new TH2F("fArmenterosNoRhoTrue_Reconstr_Pos", "fArmenterosNoRhoTrue_Reconstr_Pos", 200, -1, 1, 400, 0, 0.8);
    fArmenterosNoRhoTrue_Reconstr_Pos->GetXaxis()->SetTitle("#alpha");
    fArmenterosNoRhoTrue_Reconstr_Pos->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoRhoTrue_Reconstr_Neg = new TH2F("fArmenterosNoRhoTrue_Reconstr_Neg", "fArmenterosNoRhoTrue_Reconstr_Neg", 200, -1, 1, 400, 0, 0.8);
    fArmenterosNoRhoTrue_Reconstr_Neg->GetXaxis()->SetTitle("#alpha");
    fArmenterosNoRhoTrue_Reconstr_Neg->GetYaxis()->SetTitle("#it{q}_{T} (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoRhoTrue_Reconstr_qtDaughBoth = new TH2F("fArmenterosNoRhoTrue_Reconstr_qtDaughBoth", "fArmenterosNoRhoTrue_Reconstr_qtDaughBoth", 400, 0, 0.8, 400, 0, 0.8);
    fArmenterosNoRhoTrue_Reconstr_qtDaughBoth->GetXaxis()->SetTitle("#it{q}_{T}_pos (GeV/#it{c})");
    fArmenterosNoRhoTrue_Reconstr_qtDaughBoth->GetYaxis()->SetTitle("#it{q}_{T}_neg (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth = new TH2F("fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth", "fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth", 200, -1, 1, 200, -1, 1);
    fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth->GetXaxis()->SetTitle("#alpha_pos");
    fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth->GetYaxis()->SetTitle("#alpha_neg");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosRhoTrue_Reconstr_qtDaughBoth = new TH2F("fArmenterosRhoTrue_Reconstr_qtDaughBoth", "fArmenterosRhoTrue_Reconstr_qtDaughBoth", 400, 0, 0.8, 400, 0, 0.8);
    fArmenterosRhoTrue_Reconstr_qtDaughBoth->GetXaxis()->SetTitle("#it{q}_{T}_pos (GeV/#it{c})");
    fArmenterosRhoTrue_Reconstr_qtDaughBoth->GetYaxis()->SetTitle("#it{q}_{T}_neg (GeV/#it{c})");

    // Create histogram to track the armenteros values of MC true rhos  with negDaughter
    fArmenterosRhoTrue_Reconstr_alphaDaughBoth = new TH2F("fArmenterosRhoTrue_Reconstr_alphaDaughBoth", "fArmenterosRhoTrue_Reconstr_alphaDaughBoth", 200, -1, 1, 200, -1, 1);
    fArmenterosRhoTrue_Reconstr_alphaDaughBoth->GetXaxis()->SetTitle("#alpha_pos");
    fArmenterosRhoTrue_Reconstr_alphaDaughBoth->GetYaxis()->SetTitle("#alpha_neg");

    fHistListRhoMCTrue->Add(fArmenterosRhoTrue); // Add the histogram to your output list
    fHistListRhoMCTrue->Add(fArmenterosRhoTrue_Reconstr);
    fHistListRhoMCTrue->Add(fArmenterosNoCommonMother_Pos);
    fHistListRhoMCTrue->Add(fArmenterosNoCommonMother_Neg);
    fHistListRhoMCTrue->Add(fArmenterosNoCommonMother_qtDaughBoth);
    fHistListRhoMCTrue->Add(fArmenterosNoCommonMother_alphaDaughBoth);
    fHistListRhoMCTrue->Add(fArmenterosNoRhoTrue_Reconstr_Pos);
    fHistListRhoMCTrue->Add(fArmenterosNoRhoTrue_Reconstr_Neg);
    fHistListRhoMCTrue->Add(fArmenterosNoRhoTrue_Reconstr_qtDaughBoth);
    fHistListRhoMCTrue->Add(fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth);
    fHistListRhoMCTrue->Add(fArmenterosRhoTrue_Reconstr_qtDaughBoth);
    fHistListRhoMCTrue->Add(fArmenterosRhoTrue_Reconstr_alphaDaughBoth);
    fHistListRhoMCTrue->Add(fHist1D_pt_RhoTrue); // TH1F for pT spectrum
    fHistListRhoMCTrue->Add(fHist2D_massVSpt_RhoTrue);
    fHistListRhoMCTrue->Add(fHist2D_pt1VSpt2_RhoTrue);

    // Histograms for Checking the pT vs mT for the different components
    fHist2D_pTvsmT_noPions = new TH2F("fHist2D_pTvsmT_noPions", "fHist2D_pTvsmT_noPions", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_noPions->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_noPions->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_noPrims = new TH2F("fHist2D_pTvsmT_noPrims", "fHist2D_pTvsmT_noPrims", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_noPrims->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_noPrims->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_noCommonMother = new TH2F("fHist2D_pTvsmT_noCommonMother", "fHist2D_pTvsmT_noCommonMother", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_noCommonMother->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_noCommonMother->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_noRho = new TH2F("fHist2D_pTvsmT_noRho", "fHist2D_pTvsmT_noRho", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_noRho->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_noRho->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_noRho_MC = new TH2F("fHist2D_pTvsmT_noRho_MC", "fHist2D_pTvsmT_noRho_MC", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_noRho_MC->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_noRho_MC->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_PDGvsmT_noRho_MC = new TH2F("fHist2D_PDGvsmT_noRho_MC", "fHist2D_PDGvsmT_noRho_MC", 5001, -2500, 2500, 500, 0.0, 5.0);
    fHist2D_PDGvsmT_noRho_MC->GetXaxis()->SetTitle("PDG");
    fHist2D_PDGvsmT_noRho_MC->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_isRho = new TH2F("fHist2D_pTvsmT_isRho", "fHist2D_pTvsmT_isRho", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_isRho->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_isRho->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");
    fHist2D_pTvsmT_isRho_MC = new TH2F("fHist2D_pTvsmT_isRho_MC", "fHist2D_pTvsmT_isRho_MC", 500, 0.0, 5.0, 500, 0.0, 5.0);
    fHist2D_pTvsmT_isRho_MC->GetXaxis()->SetTitle("pT (GeV/c)");
    fHist2D_pTvsmT_isRho_MC->GetYaxis()->SetTitle("Invariant Mass (GeV/c^2)");

    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_noPions);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_noPrims);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_noCommonMother);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_noRho);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_noRho_MC);
    fHistListRhoMCTrue->Add(fHist2D_PDGvsmT_noRho_MC);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_isRho);
    fHistListRhoMCTrue->Add(fHist2D_pTvsmT_isRho_MC);

    fOutput->Add(fHistListRhoMCTrue); // Add the histogram to your output list
  }

  PostData(1, fOutput);
}

// Helper functions

void AliAnalysisTaskFemtoDreamRho::UserExec(Option_t *)
{
  AliAODEvent *Event = static_cast<AliAODEvent *>(fInputEvent);
  if (!Event)
  {
    AliWarning("No Input Event");
  }

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent))
    return;

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track)
      continue;
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  // fTrackneg->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // First we want to combine all charged pions with each other in the SE in order to find Rhos
  static std::vector<AliFemtoDreamBasePart> Particles; // pi+ candidates
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles; // pi- candidates
  AntiParticles.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles; // Rho candidates
  V0Particles.clear();
  static std::vector<AliFemtoDreamBasePart> Protons; // proton candidates
  Protons.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtons; // anti-proton candidates
  AntiProtons.clear();
  static std::vector<AliFemtoDreamBasePart> RhoMcTruePart; // Rhos verified via MC information
  RhoMcTruePart.clear();
  static std::vector<AliFemtoDreamBasePart> ProtonMcTruePart; // protons verified via MC information
  ProtonMcTruePart.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtonMcTruePart; // anti-protons verified via MC information
  AntiProtonMcTruePart.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles_SameCharge; // "fake" Rho candidates
  V0Particles_SameCharge.clear();

  // for mc matching and QA plots
  static std::vector<int> Particles_Combinations; // pi+ candidates combinations
  Particles_Combinations.clear();
  static std::vector<int> AntiParticles_Combinations; // pi- candidates combinations
  AntiParticles_Combinations.clear();
  static std::vector<TLorentzVector> Particles_Combinations_LV; // pi+ candidates combinations
  Particles_Combinations_LV.clear();
  static std::vector<TLorentzVector> AntiParticles_Combinations_LV; // pi- candidates combinations
  AntiParticles_Combinations_LV.clear();
  static std::vector<int> Ancestor_Combinations; // Same or different Ancestor
  Ancestor_Combinations.clear();
  // for selecting the MC checked rho0 candidates
  static std::vector<AliFemtoDreamBasePart> V0Particles_MC_verified; // rho0 candidates
  V0Particles_MC_verified.clear();

  // for providing the data input of the projection method
  static std::vector<AliFemtoDreamBasePart> Particles_Minv; // pi+ candidates in Minv selection window of M(pipi)
  Particles_Minv.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles_Minv; // pi- candidates in Minv selection window of M(pipi)
  AntiParticles_Minv.clear();
  // for QA plots of tracks with Minv selection
  static std::vector<int> Tracks_Particles_Minv; // Same or different Ancestor
  Tracks_Particles_Minv.clear();
  static std::vector<int> Tracks_AntiParticles_Minv; // Same or different Ancestor
  Tracks_AntiParticles_Minv.clear();

  static float massChargedPion =
      TDatabasePDG::Instance()->GetParticle(fPosPionCuts->GetPDGCode())->Mass(); // as usual to minimize uncert.
  fRhoParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // Loop to identify all charged pions & protons
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack)
  {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track)
      continue;
    fTrack->SetTrack(track);
    if (fPosPionCuts->isSelected(fTrack))
    {
      fTrack->SetInvMass(massChargedPion); // Since we combine these later we set the inv. mass to the PDG value
      Particles.push_back(*fTrack);
      Tracks_Particles_Minv.push_back(iTrack);
    }
    if (fNegPionCuts->isSelected(fTrack))
    {
      fTrack->SetInvMass(massChargedPion); // Since we combine these later we set the inv. mass to the PDG value
      AntiParticles.push_back(*fTrack);
      Tracks_AntiParticles_Minv.push_back(iTrack);
    }
    if (fPosProtonCuts->isSelected(fTrack))
    {
      Protons.push_back(*fTrack);
    }
    if (fNegProtonCuts->isSelected(fTrack))
    {
      AntiProtons.push_back(*fTrack);
    }
  }

  if (fIsSameCharge)
  {
    // Don't pair the same particles.
    for (size_t i = 0; i < Particles.size(); ++i)
    {
      const auto &posPion1 = Particles[i]; // First particle

      for (size_t j = 0; j < Particles.size(); ++j)
      {
        if (i == j)
        {
          continue; // Skip pairing the same particle
        }
        const auto &posPion2 = Particles[j]; // Second particle

        fRhoParticle->Setv0SameCharge(posPion1, posPion2, Event, false, false, true, fIsSameCharge);
      }

      const float pT_rho_candidate = fRhoParticle->GetPt();
      // At a pT > 1.8 GeV we start to see the rho in the M_inv (for now hard-coded can be optimized)
      if (pT_rho_candidate < frhoPtThreshold - 0.0001) //
      {
        continue;
      }

      if (fRhoCuts->isSelected(fRhoParticle)) // Check for proper Rho candidates, just Minv cut and kaon reject.
      {
        V0Particles_SameCharge.push_back(*fRhoParticle);
      }
    }
  }

  // Construct the V0 for the Rho decay, just simple combinatorics for now
  int counter = 0;
  int counter_tracks_Part = -1;
  int counter_tracks_antiPart = -1;
  for (const auto &posPion : Particles)
  { // Build charged pion pairs!
    counter_tracks_Part++;
    for (const auto &negPion : AntiParticles)
    {
      counter_tracks_antiPart++;
      fRhoParticle->Setv0(posPion, negPion, Event, false, false, true);

      const float pT_rho_candidate = fRhoParticle->GetPt();
      // At a pT > 1.8 GeV we start to see the rho in the M_inv (for now hard-coded can be optimized)
      if (pT_rho_candidate < frhoPtThreshold - 0.0001) //
      {
        continue;
      }

      if (fDoAncestors && fIsMC && AncestorIsSelected(fRhoParticle, fRhoCuts)) // Select everything as the RhoCandidate except the mass
      {
        bool isCommon = CommonAncestors(posPion, negPion, Event, true);
        Ancestor_Combinations.push_back(isCommon);

        // prepare plots for the pT vs minv
        if (isCommon) // isCommon
        {
          int pdg_resonance = -99999;
          bool isResonance = CommonResonance(posPion, negPion, pdg_resonance, Event, true);
          FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM);
          if (isResonance)
          {
            FillAncestorHist2D_PDGvsMinv(posPion, negPion, fHist2D_PDGvsMInv_CommonAncestorResonances, pdg_resonance);

            // Use a switch statement to handle different values of pdg_resonance
            switch (pdg_resonance)
            {
            case 311: // K0
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
              break;
            case 310: // K0short
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
              break;
            case 130: // K0long
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
              break;
            case 113: // RhoMeson
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances);
              break;
            case 223: // OmegaMeson
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances);
              break;
            case 9000221: // F0Meson
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances);
              break;
            case 9000223: // F2Meson
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances);
              break;
            default:
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances);
              break;
            }
          }
          else
          {
            FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances);
          }
        }
        else
        {
          FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateUncommonFullInvM);
        }
      }

      if (fDoProjections)
      { // Also include a selection on the rho candidates pT
        Particles_Minv.push_back(posPion);
        AntiParticles_Minv.push_back(negPion); // negPion  TEST ONLY
        // Get QA of the pion tracks used for the projection, careful here tracks will be counted multiple times
        AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(Tracks_Particles_Minv[counter_tracks_Part]));
        fTrack->SetTrack(track);
        fPosPionMinvCuts->isSelected(fTrack);

        // AliAODTrack *trackneg = static_cast<AliAODTrack *>(Event->GetTrack(Tracks_AntiParticles_Minv[counter_tracks_antiPart]));
        // fTrackneg->SetTrack(trackneg);
        // fNegPionMinvCuts->isSelected(fTrackneg);

        // // Check the kStar
        // float kStar_check = RelativePairMomentum_check(negPion, 211, posPion, 211);

        // //  Debug
        // const float invMassRho = fRhoParticle->Getv0Mass();
        // float posPdeb[3], negPdeb[3];
        // posPion.GetMomentum().GetXYZ(posPdeb);
        // negPion.GetMomentum().GetXYZ(negPdeb);
        // TLorentzVector trackPosDeb, trackNegDeb;
        // const float invPiPlus = posPion.GetInvMass();
        // const float invPiMinus = negPion.GetInvMass();
        // trackPosDeb.SetXYZM(posPdeb[0], posPdeb[1], posPdeb[2], invPiPlus);
        // trackNegDeb.SetXYZM(negPdeb[0], negPdeb[1], negPdeb[2], invPiMinus);
        // TLorentzVector trackSumdeb = trackPosDeb + trackNegDeb;
        // const float invMassPions = trackSumdeb.M();
        // if (invMassRho - invMassPions > 0.0001)
        // {
        //   printf("+++++++++++++++++++++++++++++++++++CAUTION!!+++++++++++++++++++++++++++++++++++++\n");
        // }
        // printf("Check initial minv assignment: %.4f(invPiPlus), %.4f(invPiMinus)\n", invPiPlus, invPiMinus);
        // printf("Check values: %.4f(invMassRho), %.4f(invMassPions)\n", invMassRho, invMassPions);
        // printf("Check kStar calculation: %.4f\n", kStar_check);
      }

      if (fRhoCuts->isSelected(fRhoParticle)) // Check for proper Rho candidates, just Minv cut and kaon reject.
      {
        // Also include a selection on the rho pT (for better control of what goes in the Cf)
        V0Particles.push_back(*fRhoParticle);
        if (fIsMC)
        { // store the combinations for the MC matching
          // also store kinematic distributions needed later
          Particles_Combinations.push_back(posPion.GetID());
          AntiParticles_Combinations.push_back(negPion.GetID());
          float posP[3], negP[3];
          posPion.GetMomentum().GetXYZ(posP);
          negPion.GetMomentum().GetXYZ(negP);
          TLorentzVector trackPos, trackNeg;
          trackPos.SetXYZM(posP[0], posP[1], posP[2], posPion.GetInvMass());
          trackNeg.SetXYZM(negP[0], negP[1], negP[2], negPion.GetInvMass());
          Particles_Combinations_LV.push_back(trackPos);
          AntiParticles_Combinations_LV.push_back(trackNeg);
          // temp move this here in order to check the resonances
          /*if (fDoAncestors && fIsMC && AncestorIsSelected(fRhoParticle, fRhoCuts)) // Select everything as the RhoCandidate except the mass
          {
            bool isCommon = CommonAncestors(posPion, negPion, Event, true);
            Ancestor_Combinations.push_back(isCommon);
            // prepare plots for the pT vs minv
            if (isCommon) // isCommon
            {
              int pdg_resonance = -99999;
              bool isResonance = CommonResonance(posPion, negPion, pdg_resonance, Event, true);
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM);
              if (isResonance)
              {
                FillAncestorHist2D_PDGvsMinv(posPion, negPion, fHist2D_PDGvsMInv_CommonAncestorResonances, pdg_resonance);

                // std::cout << "pdg_resonance: " << pdg_resonance << std::endl;

                // Use a switch statement to handle different values of pdg_resonance
                switch (pdg_resonance)
                {
                case 310: // K0short
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
                  break;
                case 130: // K0Long
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
                  break;
                case 311: // K0
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances);
                  break;
                case 113: // RhoMeson
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances);
                  break;
                case 223: // OmegaMeson
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances);
                  break;
                case 9000221: // F0Meson
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances);
                  break;
                case 9000223: // F2Meson
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances);
                  break;
                default:
                  FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances);
                  break;
                }
              }
              else
              {
                FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances);
              }
            }
            else
            {
              FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateUncommonFullInvM);
            }
          }*/
        }
        if (fDoAncestors && fIsMC)
        {
          bool isCommon = CommonAncestors(posPion, negPion, Event, true);
          //  prepare plots for the pT vs minv
          if (isCommon)
          {
            FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateCommon);
          }
          else
          {
            FillAncestorHist2D_pTvsMinv(posPion, negPion, fHist2D_massVSpt_RhoCandidateUncommon);
          }
        }
      }
    }
  }

  // Implement here the matching of the self-reconstructed rhos to the MC truth
  // Logic:
  /*- Check if the MC is available
    - Obtain the daughter IDs and match them to their MC instance
    - Check if the MC daughters are there
    - Check the PDG code of the daughers
    - Check if the daughters are PhysicalPrimary (should be true as this is the ALICE definition)
    - Check if they share their mother
    - Check if the mother is a rho (check if the particle was reconstructed?, this may not be needed as the rho is never measured anyways)
    - Check if the rho is physicalprimary
    - Check the mother of the rho (is most likely -1)
    - Record the kinematic distributions of the rho (true) and also the reconstructed kinematics from the pions in order to also have them with detector effects
    - log also all the of the other possibilies to see if there is some discrimination power (pTvsMassInv, Armenteros, etc. pp.)
    - Record also for these cases the un/common ancestors*/

  if (fIsMC) // fIsMC
  {
    // counter for the labels
    int counterLabels = 0;
    // PDG codes
    const float pdgIdealDaughters = 211.; // charged pions 211
    const float pdgIdealMother = 113.;    // rho 113

    for (const auto &V0Candidates : V0Particles)
    {

      TLorentzVector trackPos = Particles_Combinations_LV[counterLabels];
      TLorentzVector trackNeg = AntiParticles_Combinations_LV[counterLabels];
      TLorentzVector trackSum = trackPos + trackNeg;
      const float daughterInvM = trackSum.M();
      const float motherpT = V0Candidates.GetPt();

      TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(Event->FindListObject(
          AliAODMCParticle::StdBranchName()));
      if (!fArrayMCAOD)
      {
        AliFatal("No MC Array found\n");
      }
      int noPart = fArrayMCAOD->GetEntriesFast();

      // const int posID = V0.GetIDTracks().at(1); // this seems to be wrongly assigned need still to figure out why the class does not save daughter info Update: The track object does not save this would require some changes for now use the construction below
      // const int negID = V0.GetIDTracks().at(0); // this seems to be wrongly assigned need still to figure out why the class does not save daughter info

      const int posID = (int)Particles_Combinations[counterLabels];     // this is a fix for the issue above
      const int negID = (int)AntiParticles_Combinations[counterLabels]; // this is a fix for the issue above

      // Do this only if we explicitly request the ancestors
      if (fDoAncestors)
      {
        // only check those which have the same ancestor
        bool sameAncestor = Ancestor_Combinations[counterLabels];
        counterLabels++; // increase the counter for the combination

        if (!sameAncestor)
        {
          // in this case we have background
          if (fIsMCTrueRhoCombBkrg)
          {
            V0Particles_MC_verified.push_back(V0Candidates);
          }
          continue;
        }
      }

      if (posID > noPart || negID > noPart)
      {
        AliFatal("ID of daughter track has no MCparticle entity\n");
      }

      AliAODMCParticle *mcPartPos = (AliAODMCParticle *)fArrayMCAOD->At(posID);
      AliAODMCParticle *mcPartNeg = (AliAODMCParticle *)fArrayMCAOD->At(negID);

      if (!mcPartPos || !mcPartNeg)
      {
        AliFatal("MC particle for daughters are NULL\n");
      }

      // obtain kinematic variables

      const float posDauPDG = mcPartPos->GetPdgCode();
      const float negDauPDG = mcPartNeg->GetPdgCode();

      // Seperate here the common template according to all contributions?

      if (pdgIdealDaughters != std::abs(posDauPDG) || pdgIdealDaughters != std::abs(negDauPDG))
      {
        // track the PDG for the misidenified PDG1 vs PDG2
        // track their pTvsMass
        fHist2D_pTvsmT_noPions->Fill(motherpT, daughterInvM);

        // in this case we have background
        if (fIsMCTrueRhoCombBkrg)
        {
          V0Particles_MC_verified.push_back(V0Candidates);
        }

        continue; // jump to next combination
      }

      if (!mcPartPos->IsPhysicalPrimary() || !mcPartNeg->IsPhysicalPrimary())
      {
        // track the PDG for the misidenified PDG1 vs PDG2
        // track their pTvsMass
        fHist2D_pTvsmT_noPrims->Fill(motherpT, daughterInvM);

        // in this case we have background
        if (fIsMCTrueRhoCombBkrg)
        {
          V0Particles_MC_verified.push_back(V0Candidates);
        }

        continue; // jump to next combination
      }

      const int motherIDposDau = mcPartPos->GetMother();
      const int motherIDnegDau = mcPartNeg->GetMother();

      if (motherIDposDau > noPart || motherIDnegDau > noPart)
      {
        AliFatal("ID of mother track has no MCparticle entity\n");
      }

      AliAODMCParticle *mcPartPosMother = (AliAODMCParticle *)fArrayMCAOD->At(motherIDposDau);
      AliAODMCParticle *mcPartNegMother = (AliAODMCParticle *)fArrayMCAOD->At(motherIDnegDau);

      if (!mcPartPosMother || !mcPartNegMother)
      {
        AliFatal("MC particle for mother is NULL\n");
      }

      if (motherIDposDau != motherIDnegDau)
      {
        // track their pTvsMass
        // Fill histograms with some data (replace this with your data)
        fHist2D_pTvsmT_noCommonMother->Fill(motherpT, daughterInvM);
        // track Armenteros plot
        float alpha_NoCommonMother_pos = 0;
        float qT_NoCommonMother_pos = 0;
        CalculateAlphaAndQT(V0Candidates, trackPos, trackNeg, alpha_NoCommonMother_pos, qT_NoCommonMother_pos);
        fArmenterosNoCommonMother_Pos->Fill(alpha_NoCommonMother_pos, qT_NoCommonMother_pos);
        float alpha_NoCommonMother_neg = 0;
        float qT_NoCommonMother_neg = 0;
        CalculateAlphaAndQT(V0Candidates, trackNeg, trackPos, alpha_NoCommonMother_neg, qT_NoCommonMother_neg); // switched pos and neg track
        fArmenterosNoCommonMother_Neg->Fill(alpha_NoCommonMother_neg, qT_NoCommonMother_neg);
        fArmenterosNoCommonMother_qtDaughBoth->Fill(qT_NoCommonMother_pos, qT_NoCommonMother_neg);
        fArmenterosNoCommonMother_alphaDaughBoth->Fill(alpha_NoCommonMother_pos, alpha_NoCommonMother_neg);

        // in this case we have background
        // this is the same check as for the common ancestor above, but will not lead to double counting due to the use of continue and enables different flag combs.
        if (fIsMCTrueRhoCombBkrg)
        {
          V0Particles_MC_verified.push_back(V0Candidates);
        }

        continue;
      }

      const float negMotherPDG = mcPartNegMother->GetPdgCode();
      const float posMotherPDG = mcPartPosMother->GetPdgCode();

      if (posMotherPDG != pdgIdealMother)
      { // This will contain everything mini-jet+Resonances
        // track their pTvsMass
        // Fill histograms with some data (replace this with your data)
        fHist2D_pTvsmT_noRho->Fill(motherpT, daughterInvM);
        fHist2D_pTvsmT_noRho_MC->Fill(mcPartPosMother->Pt(), mcPartPosMother->GetCalcMass()); // here the pos mother was chosen could also take the neg one
        fHist2D_PDGvsmT_noRho_MC->Fill(posMotherPDG, mcPartPosMother->GetCalcMass());
        fHist2D_PDGvsmT_noRho_MC->Fill(-negMotherPDG, mcPartNegMother->GetCalcMass());
        // track Armenteros plot
        float alpha_pos = 0;
        float qT_pos = 0;
        CalculateAlphaAndQT(V0Candidates, trackPos, trackNeg, alpha_pos, qT_pos);
        fArmenterosNoRhoTrue_Reconstr_Pos->Fill(alpha_pos, qT_pos);
        float alpha_neg = 0;
        float qT_neg = 0;
        CalculateAlphaAndQT(V0Candidates, trackNeg, trackPos, alpha_neg, qT_neg); // switched pos and neg track
        fArmenterosNoRhoTrue_Reconstr_Neg->Fill(alpha_neg, qT_neg);
        fArmenterosNoRhoTrue_Reconstr_qtDaughBoth->Fill(qT_pos, qT_neg);
        fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth->Fill(alpha_pos, alpha_neg);

        // in this case we have background
        if (fIsMCTrueRhoCombBkrg)
        {
          V0Particles_MC_verified.push_back(V0Candidates);
        }

        continue;
      }
      //  At this point we have a proper rho candidate,
      //  now we can record the kinemtic distribution from the rhoMC instance reconstructed and MCtruth

      fHist2D_pTvsmT_isRho->Fill(motherpT, daughterInvM);
      fHist2D_pTvsmT_isRho_MC->Fill(mcPartPosMother->Pt(), mcPartPosMother->GetCalcMass());
      float alphapos = 0;
      float qTpos = 0;
      CalculateAlphaAndQT(V0Candidates, trackPos, trackNeg, alphapos, qTpos);
      float alphaneg = 0;
      float qTneg = 0;
      CalculateAlphaAndQT(V0Candidates, trackNeg, trackPos, alphaneg, qTneg); // switched pos and neg track
      fArmenterosRhoTrue_Reconstr->Fill(alphapos, qTpos);
      fArmenterosRhoTrue_Reconstr_qtDaughBoth->Fill(qTpos, qTneg);
      fArmenterosRhoTrue_Reconstr_alphaDaughBoth->Fill(alphapos, alphaneg);

      if (!fIsMCTrueRhoCombBkrg)
      {
        V0Particles_MC_verified.push_back(V0Candidates);
      }

      //  track their pTvsMass (reconstr)
      //  track Armenteros plot (reconstr)
      //  track their pTvsMass (mcTruth)
      //  track Armenteros plot (mcTruth)
    }

    // Use the following to study the MC true distributions (inlcudes only acceptancy no momentum resolution)
    // Very simple construction without checking motherIDs etc. or properties of the daughers as the MC vector breaks the particle collection later on.
    ///________________________________________________________________________________________________________
    if (fIsMC && fDoMcTruth) // fDoMcTruth
    {
      AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      if (!eventHandler)
      {
        AliWarning("No eventHandler for kine Dist");
        return;
      }
      AliMCEvent *fMC = eventHandler->MCEvent();
      if (!fMC)
      {
        AliWarning("No fMC for kine Dist");
        return;
      }

      const float bfield = fMC->GetMagneticField();

      TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(
          Event->FindListObject(AliAODMCParticle::StdBranchName()));
      int noPart = fArrayMCAOD->GetEntriesFast();
      int mcpdg;
      AliFemtoDreamBasePart part;
      for (int iPart = 1; iPart < noPart; iPart++)
      {
        AliAODMCParticle *mcPart = (AliAODMCParticle *)fArrayMCAOD->At(iPart);
        if (!(mcPart))
        {
          continue;
        }
        if (mcPart->GetLabel() < 0)
        {
          continue;
        }
        mcpdg = mcPart->GetPdgCode();
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();

        // this is only for debug for now always off, since else we could select only those particles which are not rho but primary and decay into pi+pi- i.e. we miss the combinatorial background. Only a few omega mesons survive this selection as it should be.
        if (fIsMCTrueRhoCombBkrg && false) // Only needed for a fast cross-check but better to comment this part
        {                                  // If we only want the comb. background made out of (pipi)_{combinatorial}p triplets then we simply skip all true rhos,
          // but as we have purity of the pi and p we will still require them to be true pi and p

          if (mcpdg != 113)
          {

            if (pt < frhoPtThreshold - 0.0001) // frhoPtThreshold
            {
              continue;
            }
            int firstDaughter = mcPart->GetDaughterFirst();
            int lastDaughter = mcPart->GetDaughterLast();
            if (firstDaughter > noPart || lastDaughter > noPart)
            {
              continue; // sanity check
            }
            AliAODMCParticle *mcDaughterOne =
                (AliAODMCParticle *)fArrayMCAOD->At(firstDaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?
            AliAODMCParticle *mcDaughterTwo =
                (AliAODMCParticle *)fArrayMCAOD->At(lastDaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?

            if (!mcDaughterOne || !mcDaughterTwo) // MC info of the daughters cannot be accessed
            {
              continue;
            }
            int dpdgOne = mcDaughterOne->GetPdgCode();
            int dpdgTwo = mcDaughterTwo->GetPdgCode();
            const float pdgIdealDaughters = 211.;                                                     // Hardcoded since B.R. into pions is nearly 100%, rejects decays with gammas
            if ((std::abs(dpdgOne) != pdgIdealDaughters) || (std::abs(dpdgTwo) != pdgIdealDaughters)) // MC info of the daughters cannot be accessed
            {
              continue;
            }
            double dptOne = mcDaughterOne->Pt();
            double detaOne = mcDaughterOne->Eta();
            if ((dptOne > 4.0 || dptOne < 0.14) || (detaOne < -0.8 || detaOne > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
            {
              continue;
            }
            double dptTwo = mcDaughterTwo->Pt();
            double detaTwo = mcDaughterTwo->Eta();
            if ((dptTwo > 4.0 || dptTwo < 0.14) || (detaTwo < -0.8 || detaTwo > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
            {
              continue;
            }

            part.SetMCParticleRePart(mcPart);
            part.SetID(mcPart->GetLabel());
            part.SetMCParticle(mcPart, fMC);
            part.SetIDTracks(mcPart->GetLabel());
            // Needed as else the QA flags throw an error
            SetPhiAtRadiusMCTruth(part, bfield);
            // Correct the value of the internal data member
            part.SetInvMass(mcPart->GetCalcMass());
            // here we could select on the Minv
            RhoMcTruePart.push_back(part);

            // QA the distributions of the daughters!
            // Armenteros Podolanski plot for the MC True particles.
            float alpha = 0;
            float qT = 0;
            CalculateAlphaAndQT(mcPart, mcDaughterOne, mcDaughterTwo, alpha, qT);
            // Make a histogram of the pT of the rho vs the m_inv
            // Record the kinematic distribution of the rhos

            // fpTCorrerrorHistogram->Fill(dpt, dpt2); // record correlation betweent the daughter momenta
            // Fill all the hists
            fHist1D_pt_RhoTrue->Fill(part.GetPt());                              // this is with acceptance cuts
            fHist2D_massVSpt_RhoTrue->Fill(part.GetPt(), mcPart->GetCalcMass()); // this is with acceptance cuts
            fArmenterosRhoTrue->Fill(alpha, qT);                                 // this is with acceptance cuts
            fHist2D_pt1VSpt2_RhoTrue->Fill(dptOne, dptTwo);                      // this is with acceptance cuts
            // track momentum difference (this is tracked in the Armenteros-Plot)
          }
        }
        else // Select all the MC true rhos
        {
          if (mcpdg == 113)
          {

            if (pt < frhoPtThreshold - 0.0001) //
            {
              continue;
            }
            int firstDaughter = mcPart->GetDaughterFirst();
            int lastDaughter = mcPart->GetDaughterLast();
            if (firstDaughter > noPart || lastDaughter > noPart)
            {
              continue; // sanity check
            }
            AliAODMCParticle *mcDaughterOne =
                (AliAODMCParticle *)fArrayMCAOD->At(firstDaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?
            AliAODMCParticle *mcDaughterTwo =
                (AliAODMCParticle *)fArrayMCAOD->At(lastDaughter); // also get the second one as we want to see both momenta, maybe even correalted or in amenteros plot?

            if (!mcDaughterOne || !mcDaughterTwo) // MC info of the daughters cannot be accessed
            {
              continue;
            }
            int dpdgOne = mcDaughterOne->GetPdgCode();
            int dpdgTwo = mcDaughterTwo->GetPdgCode();
            const float pdgIdealDaughters = 211.;                                                     // Hardcoded since B.R. into pions is nearly 100%, rejects decays with gammas
            if ((std::abs(dpdgOne) != pdgIdealDaughters) || (std::abs(dpdgTwo) != pdgIdealDaughters)) // MC info of the daughters cannot be accessed
            {
              continue;
            }
            double dptOne = mcDaughterOne->Pt();
            double detaOne = mcDaughterOne->Eta();
            if ((dptOne > 4.0 || dptOne < 0.14) || (detaOne < -0.8 || detaOne > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
            {
              continue;
            }
            double dptTwo = mcDaughterTwo->Pt();
            double detaTwo = mcDaughterTwo->Eta();
            if ((dptTwo > 4.0 || dptTwo < 0.14) || (detaTwo < -0.8 || detaTwo > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
            {
              continue;
            }

            part.SetMCParticleRePart(mcPart);
            part.SetID(mcPart->GetLabel());
            part.SetMCParticle(mcPart, fMC);
            part.SetIDTracks(mcPart->GetLabel());
            // Needed as else the QA flags throw an error
            SetPhiAtRadiusMCTruth(part, bfield);
            // Correct the value of the internal data member
            part.SetInvMass(mcPart->GetCalcMass());
            // here we could select on the Minv
            RhoMcTruePart.push_back(part);

            // QA the distributions of the daughters!
            // Armenteros Podolanski plot for the MC True particles.
            float alpha = 0;
            float qT = 0;
            CalculateAlphaAndQT(mcPart, mcDaughterOne, mcDaughterTwo, alpha, qT);
            // Make a histogram of the pT of the rho vs the m_inv
            // Record the kinematic distribution of the rhos

            // fpTCorrerrorHistogram->Fill(dpt, dpt2); // record correlation betweent the daughter momenta
            // Fill all the hists
            fHist1D_pt_RhoTrue->Fill(part.GetPt());                              // this is with acceptance cuts
            fHist2D_massVSpt_RhoTrue->Fill(part.GetPt(), mcPart->GetCalcMass()); // this is with acceptance cuts
            fArmenterosRhoTrue->Fill(alpha, qT);                                 // this is with acceptance cuts
            fHist2D_pt1VSpt2_RhoTrue->Fill(dptOne, dptTwo);                      // this is with acceptance cuts
            // track momentum difference (this is tracked in the Armenteros-Plot)
          }
        }
        if (mcpdg == 2212) // Select all the MC true protons within acceptance
        {
          if ((pt > 4.05 || pt < 0.5) || (eta < -0.8 || eta > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
          {
            continue;
          }
          part.SetMCParticleRePart(mcPart); // Needed as else the QA flags throw an error
          SetPhiAtRadiusMCTruth(part, bfield);
          ProtonMcTruePart.push_back(part);
        }
        if (mcpdg == -2212) // Select all the MC true anti-protons within acceptance
        {
          if ((pt > 4.05 || pt < 0.5) || (eta < -0.8 || eta > 0.8)) // within acceptance (refine to have the acceptance cuts on the daughters)
          {
            continue;
          }
          part.SetMCParticleRePart(mcPart); // Needed as else the QA flags throw an error
          SetPhiAtRadiusMCTruth(part, bfield);
          AntiProtonMcTruePart.push_back(part);
        }
      }
    }
  }

  ///________________________________________________________________________________________________________

  // For now no default pair cleaning but may change to CleanTrackandTrack of same charge pion and proton
  if (fDoCleaning)
  {
    // fPairCleaner->CleanDecay(&V0Particles, 0); //Shouldn't be needed at the tracks have a opposite charge
    fPairCleaner->CleanTrackAndDecay(&Particles, &Protons, 0);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiProtons, 1);
    fPairCleaner->CleanTrackAndDecay(&Protons, &V0Particles, 2);     // I may need to take there the pion daughters...
    fPairCleaner->CleanTrackAndDecay(&AntiProtons, &V0Particles, 3); // I may need to take there the pion daughters...
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  if (fIsMC && fDoMcTruth) // These will never require cleaning!
  {
    if (fIsMCcheckedCombs)
    {
      // Since we don't care about efficiency effects we take there the true proton sample in order to enchange the statistics by a factor of 5. The purity is high in either case due to the optimized selection criteria.
      fPairCleaner->StoreParticle(V0Particles);          // V0Particles_MC_verified
      fPairCleaner->StoreParticle(ProtonMcTruePart);     // ProtonMcTruePart Could check what happens with the reconstructed protons instead of MC True sample
      fPairCleaner->StoreParticle(AntiProtonMcTruePart); // AntiProtonMcTruePart
    }
    else
    {
      fPairCleaner->StoreParticle(RhoMcTruePart);
      fPairCleaner->StoreParticle(ProtonMcTruePart);
      fPairCleaner->StoreParticle(AntiProtonMcTruePart);
    }
  }
  else
  {
    if (fIsSameCharge)
    {
      fPairCleaner->StoreParticle(V0Particles_SameCharge);
    }
    else
    {
      fPairCleaner->StoreParticle(V0Particles);
    }
    fPairCleaner->StoreParticle(Protons);
    fPairCleaner->StoreParticle(AntiProtons);
  }
  if (fDoProjections)
  {
    fPairCleaner->StoreParticle(Particles_Minv);
    fPairCleaner->StoreParticle(AntiParticles_Minv);
  }

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamRho::ResetGlobalTrackReference()
{
  for (UShort_t i = 0; i < fTrackBufferSize; i++)
  {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamRho::StoreGlobalTrackReference(
    AliAODTrack *track)
{
  // for documentation see AliFemtoDreamAnalysis

  const int trackID = track->GetID();
  if (trackID < 0)
  {
    return;
  }
  if (trackID >= fTrackBufferSize)
  {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }
  if (fGTI[trackID])
  {
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls()))
    {
      return;
    }
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls())
    {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

void AliAnalysisTaskFemtoDreamRho::SetPhiAtRadiusMCTruth(
    AliFemtoDreamBasePart partMC, const float bfield)
{
  // Calculate PhiAtRadius, taken from AliFemtoDreamTrack
  float TPCradii[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};
  float phi0 = partMC.GetPhi().at(0);
  float pt = partMC.GetPt();
  float chg = partMC.GetCharge().at(0);
  std::vector<float> phiatRadius;
  for (int radius = 0; radius < 9; radius++)
  {
    // 20-Feb-2022
    // Avoid NAN in asin for low momentum particle (particularly for pions)
    if (TMath::Abs(0.1 * chg * bfield * 0.3 * TPCradii[radius] * 0.01 / (2. * pt)) < 1.)
    {
      phiatRadius.push_back(phi0 - TMath::ASin(0.1 * chg * bfield * 0.3 * TPCradii[radius] * 0.01 / (2. * pt)));
    } // safety check for asin
  }
  partMC.SetPhiAtRadius(phiatRadius);
}

bool AliAnalysisTaskFemtoDreamRho::CommonAncestors(const AliFemtoDreamBasePart &part1,
                                                   const AliFemtoDreamBasePart &part2, AliAODEvent *Event, bool verbose)
{
  bool IsCommon = false;

  if (part1.GetMotherID() == part2.GetMotherID())
  {
    IsCommon = true;
  }
  else if (part1.GetMotherID() != part2.GetMotherID())
  {
    IsCommon = false;
  }

  return IsCommon;
}

void AliAnalysisTaskFemtoDreamRho::FillAncestorHist2D_pTvsMinv(const AliFemtoDreamBasePart &posDaughter,
                                                               const AliFemtoDreamBasePart &negDaughter,
                                                               TH2F *hist2D)
{
  float posP[3], negP[3];
  posDaughter.GetMomentum().GetXYZ(posP);
  negDaughter.GetMomentum().GetXYZ(negP);
  TLorentzVector trackPos, trackNeg;
  trackPos.SetXYZM(posP[0], posP[1], posP[2], posDaughter.GetInvMass());
  trackNeg.SetXYZM(negP[0], negP[1], negP[2], negDaughter.GetInvMass());
  TLorentzVector trackSum = trackPos + trackNeg;

  float V0pt = trackSum.Pt();
  float V0mInv = trackSum.Mag();

  hist2D->Fill(V0pt, V0mInv);
}

void AliAnalysisTaskFemtoDreamRho::FillAncestorHist2D_PDGvsMinv(const AliFemtoDreamBasePart &posDaughter,
                                                                const AliFemtoDreamBasePart &negDaughter,
                                                                TH2F *hist2D, int &PDG)
{
  float posP[3], negP[3];
  posDaughter.GetMomentum().GetXYZ(posP);
  negDaughter.GetMomentum().GetXYZ(negP);
  TLorentzVector trackPos, trackNeg;
  trackPos.SetXYZM(posP[0], posP[1], posP[2], posDaughter.GetInvMass());
  trackNeg.SetXYZM(negP[0], negP[1], negP[2], negDaughter.GetInvMass());
  TLorentzVector trackSum = trackPos + trackNeg;

  float V0mInv = trackSum.Mag();

  hist2D->Fill(PDG, V0mInv);
}

bool AliAnalysisTaskFemtoDreamRho::AncestorIsSelected(AliFemtoDreamv0 *v0, AliFemtoDreamv0Cuts *V0Selections)
// Use this carefully as for now this is fine as no Selection on the V0Candidate is used in the RhoSelections
// besides the invMass window which must not be used here: see CPAandMassCuts
{
  bool pass = true;
  if (!v0->IsSet())
  {
    pass = false;
  }
  if (pass)
  {
    if (V0Selections->GetCutDaughters())
    {
      if (!V0Selections->AccessDaughtersPassCuts(v0))
      {
        pass = false;
      }
    }
  }
  if (pass)
  {
    if (!V0Selections->AccessMotherPassCuts(v0))
    {
      pass = false;
    }
  }
  if (pass)
  {
    if (V0Selections->GetCutArmenteros() && !V0Selections->AccessArmenterosSelection(v0))
    {
      pass = false;
    }
  }
  if (pass)
  {
    if (V0Selections->GetKaonRejection() && !V0Selections->AccessRejectAsKaon(v0))
    {
      pass = false;
    }
  }
  // if (pass)
  //{
  //   if (!CPAandMassCuts(v0)) // Dont check the CPA and mass for now
  //   {
  //     pass = false;
  //   }
  // }

  return pass;
}

void AliAnalysisTaskFemtoDreamRho::CalculateAlphaAndQT(const AliAODMCParticle *mcPart, const AliAODMCParticle *mcDaughterOne, const AliAODMCParticle *mcDaughterTwo, float &alpha, float &qT)
{
  // Extract momentum components for mcPart
  Double_t pxmcPart = mcPart->Px();
  Double_t pymcPart = mcPart->Py();
  Double_t pzmcPart = mcPart->Pz();

  // Extract momentum components for mcDaughterOne
  Double_t pxmcDaughterOne = mcDaughterOne->Px();
  Double_t pymcDaughterOne = mcDaughterOne->Py();
  Double_t pzmcDaughterOne = mcDaughterOne->Pz();

  // Extract momentum components for mcDaughterTwo
  Double_t pxmcDaughterTwo = mcDaughterTwo->Px();
  Double_t pymcDaughterTwo = mcDaughterTwo->Py();
  Double_t pzmcDaughterTwo = mcDaughterTwo->Pz();

  // Create TVector3 objects
  TVector3 v0P(pxmcPart, pymcPart, pzmcPart);
  TVector3 posP(pxmcDaughterOne, pymcDaughterOne, pzmcDaughterOne);
  TVector3 negP(pxmcDaughterTwo, pymcDaughterTwo, pzmcDaughterTwo);

  // Calculate alpha
  alpha = 1. - 2. / (1. + posP.Dot(v0P) / negP.Dot(v0P));

  // Calculate qT
  qT = posP.Perp(v0P);
}

void AliAnalysisTaskFemtoDreamRho::CalculateAlphaAndQT(const AliFemtoDreamBasePart Part, const TLorentzVector DaughterOne, const TLorentzVector DaughterTwo, float &alpha, float &qT)
{
  // Extract momentum components for mcPart
  Double_t pxPart = Part.GetPx();
  Double_t pyPart = Part.GetPy();
  Double_t pzPart = Part.GetPz();

  // Extract momentum components for mcDaughterOne
  Double_t pxDaughterOne = DaughterOne.Px();
  Double_t pyDaughterOne = DaughterOne.Py();
  Double_t pzDaughterOne = DaughterOne.Pz();

  // Extract momentum components for mcDaughterTwo
  Double_t pxDaughterTwo = DaughterTwo.Px();
  Double_t pyDaughterTwo = DaughterTwo.Py();
  Double_t pzDaughterTwo = DaughterTwo.Pz();

  // Create TVector3 objects
  TVector3 v0P(pxPart, pyPart, pzPart);
  TVector3 posP(pxDaughterOne, pyDaughterOne, pzDaughterOne);
  TVector3 negP(pxDaughterTwo, pyDaughterTwo, pzDaughterTwo);

  // Calculate alpha
  alpha = 1. - 2. / (1. + posP.Dot(v0P) / negP.Dot(v0P));

  // Calculate qT
  qT = posP.Perp(v0P);
}

void AliAnalysisTaskFemtoDreamRho::CalculateAlphaAndQT(const AliFemtoDreamBasePart *Part, const AliFemtoDreamBasePart *DaughterOne, const AliFemtoDreamBasePart *DaughterTwo, float &alpha, float &qT)
{
  // Extract momentum components for mcPart
  Double_t pxPart = Part->GetPx();
  Double_t pyPart = Part->GetPy();
  Double_t pzPart = Part->GetPz();

  // Extract momentum components for mcDaughterOne
  Double_t pxDaughterOne = DaughterOne->GetPx();
  Double_t pyDaughterOne = DaughterOne->GetPy();
  Double_t pzDaughterOne = DaughterOne->GetPz();

  // Extract momentum components for mcDaughterTwo
  Double_t pxDaughterTwo = DaughterTwo->GetPx();
  Double_t pyDaughterTwo = DaughterTwo->GetPy();
  Double_t pzDaughterTwo = DaughterTwo->GetPz();

  // Create TVector3 objects
  TVector3 v0P(pxPart, pyPart, pzPart);
  TVector3 posP(pxDaughterOne, pyDaughterOne, pzDaughterOne);
  TVector3 negP(pxDaughterTwo, pyDaughterTwo, pzDaughterTwo);

  // Calculate alpha
  alpha = 1. - 2. / (1. + posP.Dot(v0P) / negP.Dot(v0P));

  // Calculate qT
  qT = posP.Perp(v0P);
}

bool AliAnalysisTaskFemtoDreamRho::CommonResonance(const AliFemtoDreamBasePart &part1, const AliFemtoDreamBasePart &part2, int &pdg_resonance, AliAODEvent *Event, bool verbose)
{
  bool IsResonance = true;

  if (part1.GetMotherID() != part2.GetMotherID())
  {
    AliFatal("AliAnalysisTaskFemtoDreamRho::CommonResonance: The two particle should have a common mother");
  }

  // std::cout << "part1.GetMotherID(): " << part1.GetMotherID() << "    part2.GetMotherID(): " << part2.GetMotherID() << std::endl;
  //  do the mc matching and figure out if the PDG codes of the mothers are the same as the GetMotherPDG
  TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(Event->FindListObject(
      AliAODMCParticle::StdBranchName()));
  if (!fArrayMCAOD)
  {
    AliFatal("No MC Array found\n");
  }
  // here we need to take the first MC mother of the particles.
  AliAODMCParticle *mcPartPos = (AliAODMCParticle *)fArrayMCAOD->At(part1.GetID());
  AliAODMCParticle *mcPartNeg = (AliAODMCParticle *)fArrayMCAOD->At(part2.GetID());
  if (!mcPartPos || !mcPartNeg)
  {
    AliFatal("MC particle for daughters are NULL\n");
  }
  const float posDauPDG = mcPartPos->GetPdgCode();
  const float negDauPDG = mcPartNeg->GetPdgCode();
  if (211 != std::abs(posDauPDG) || 211 != std::abs(negDauPDG)) // must be pions
  {
    return false;
  }
  // if (!mcPartPos->IsPhysicalPrimary() || !mcPartNeg->IsPhysicalPrimary()) // has to be primary
  // {
  //    return false;
  // }
  if (mcPartPos->GetMother() != mcPartNeg->GetMother()) // make sure that we have the same mother resonance
  {
    return false;
  }

  AliAODMCParticle *mcPartPosMother = (AliAODMCParticle *)fArrayMCAOD->At(mcPartPos->GetMother());
  AliAODMCParticle *mcPartNegMother = (AliAODMCParticle *)fArrayMCAOD->At(mcPartNeg->GetMother());
  const float posMotherPDG = mcPartPosMother->GetPdgCode();
  const float negMotherPDG = mcPartNegMother->GetPdgCode();

  // std::cout << "mcPartPos->GetPdgCode(): " << posDauPDG << "    mcPartNeg->GetPdgCode(): " << negDauPDG << std::endl;
  // std::cout << "part1.GetMotherPDG(): " << part1.GetMotherPDG() << "    part2.GetMotherPDG(): " << part2.GetMotherPDG() << std::endl;

  if (posMotherPDG != negMotherPDG)
  { // the ID is the same, but the PDG different -> Two tracks from same hard scattering but different resonances.
    return false;
  }

  pdg_resonance = posMotherPDG; // keep track of the resonances

  return IsResonance;
}

float AliAnalysisTaskFemtoDreamRho::RelativePairMomentum_check(
    const AliFemtoDreamBasePart &part1, const int pdg1, const AliFemtoDreamBasePart &part2,
    const int pdg2)
{
  TLorentzVector PartOne, PartTwo;
  PartOne.SetXYZM(part1.GetMomentum().X(), part1.GetMomentum().Y(),
                  part1.GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2.GetMomentum().X(), part2.GetMomentum().Y(),
                  part2.GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());
  return RelativePairMomentum_check(PartOne, PartTwo);
}

float AliAnalysisTaskFemtoDreamRho::RelativePairMomentum_check(
    TLorentzVector &PartOne, TLorentzVector &PartTwo)
{
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}