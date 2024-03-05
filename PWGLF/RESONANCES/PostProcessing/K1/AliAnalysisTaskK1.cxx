/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnalysisTaskK1
 *
 *  Test code for the reconstructing K1(1270)
 *  Output could be saved to Tree by using SetFillnTree(kTRUE)
 *    -> can be used for ML input
 *
 *  Author: Bong-Hwi Lim
 *
 */

#include <TDatabasePDG.h>
#include <math.h>
#include <tuple>

#include <iostream>

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliMultSelectionTask.h"
#include "AliPIDResponse.h"
#include "AliLog.h"
#include "TChain.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODv0.h"
#include "AliAnalysisTaskK1.h"
#include "AliAnalysisTaskTrackMixer.h"
#include "AliESDv0.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t pi = TMath::Pi();
const Double_t kPionMass = AliPID::ParticleMass(AliPID::kPion);
const Double_t kKaonMass = AliPID::ParticleMass(AliPID::kKaon);
const Double_t kK892Mass = TDatabasePDG::Instance()->GetParticle(313)->Mass();
const Double_t kRhoMass = TDatabasePDG::Instance()->GetParticle(113)->Mass();

enum
{
  kPrimaryPionPass = BIT(0),
  kSecondaryPionPass = BIT(1),
  kKaonPass = BIT(2),
};
enum
{
  kNormal = 1,
  kAnti = 2
};
enum
{                   // PDG Code
  kK1PCode = 10323, // K1(1270)+
  kK1NCode = 10313, // K1(1270)0 // will not be found
  kPionCode = 211,  // Pion+
  kKaonCode = 321,  // Kaon-
  kK892Code = 313,  // K892
  kRhoCode = 113    // Rho0
};
enum
{
  kK1P = 1,
  kK1N,
  kK1P_MIX,
  kK1N_MIX,
  kK1P_NOT,
  kK1N_NOT,
  kAllType
};
enum
{
  kK1P_GEN = 1, // 1
  kK1N_GEN,
  kK1P_GEN_INEL10,
  kK1N_GEN_INEL10,
  kK1P_GEN_INEL10_IGZ,
  kK1N_GEN_INEL10_IGZ,
  kK1P_GEN_TRIG,
  kK1N_GEN_TRIG,
  kK1P_REC,
  kK1N_REC
};
enum
{
  kAll = 1, // 1
  kINEL10,
  kINEL_trig,
  kINEL_trig_vtx,
  kINEL_trig_vtx10,
  kINELg0, // 6
  kINELg010,
  kINELg0_trig,
  kINELg0_trig_vtx,
  kINELg0_trig_vtx10,
  kSelected // 11
};

class AliAnalysisTaskK1;

ClassImp(AliAnalysisTaskK1) AliAnalysisTaskK1::AliAnalysisTaskK1()
    : AliAnalysisTaskSE(),
      fTrackCuts{nullptr},
      fPIDResponse{nullptr},
      fAODMCHeader{nullptr},
      fList{nullptr},
      fMCArray{nullptr},
      fVertex{nullptr},
      fNanoTree{nullptr},
      fNanoMCTree{nullptr},
      fRecK1daughter{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fSimK1part{0, 0, 0, 0, 0, 0, 0, 0, 0},
      fHn5DK1Data{nullptr},
      fHn5DK1MC{nullptr},
      fHn2DEvtNorm{nullptr},
      hMultiplicity{nullptr},
      hEtaTrack_before{nullptr},
      hDCAPVTrack_before{nullptr},
      hDCArPVTrack_before{nullptr},
      hPtTrack_before{nullptr},
      hEtaTrack_ppion{nullptr},
      hDCAPVTrack_ppion{nullptr},
      hDCArPVTrack_ppion{nullptr},
      hPtTrack_ppion{nullptr},
      hTPCPIDTrack_ppion{nullptr},
      hTPCPIDTrackNsigVspT_ppion{nullptr},
      hEtaTrack_spion{nullptr},
      hDCAPVTrack_spion{nullptr},
      hDCArPVTrack_spion{nullptr},
      hPtTrack_spion{nullptr},
      hTPCPIDTrack_spion{nullptr},
      hTPCPIDTrackNsigVspT_spion{nullptr},
      hEtaTrack_kaon{nullptr},
      hDCAPVTrack_kaon{nullptr},
      hDCArPVTrack_kaon{nullptr},
      hPtTrack_kaon{nullptr},
      hTPCPIDTrack_kaon{nullptr},
      hTPCPIDTrackNsigVspT_kaon{nullptr},
      hK1OA{nullptr},
      hK1PairAsymm{nullptr},
      hInvMass_k892_rho{nullptr},
      hInvMass_Secondary_PiKa{nullptr},
      hK1SecondarypT{nullptr},
      hK1OA_cut{nullptr},
      hK1PairAsymm_cut{nullptr},
      hInvMass_k892_rho_cut{nullptr},
      hInvMass_Secondary_PiKa_cut{nullptr},
      hK1SecondarypT_cut{nullptr},
      hK1OA_MCTrue{nullptr},
      hK1PairAsymm_MCTrue{nullptr},
      hInvMass_k892_rho_MCTrue{nullptr},
      hInvMass_Secondary_PiKa_MCTrue{nullptr},
      hK1SecondarypT_MCTrue{nullptr},
      fIsAOD{false},
      fIsNano{false},
      fSetMixing{false},
      fFillQAPlot{true},
      fIsMC{false},
      fIsPrimaryMC{true},
      fFillTree{false},
      fIsINEL{false},
      fIsHM{false},
      fSkipFillingHistogram{false},
      fSetModeK892orRho{true},
      fEMpool{},
      fBinCent{},
      fBinZ{},
      fPosPV{},
      fMagField{0},
      fCent{-1},
      fCustomEventID{0},
      fnMix{10},
      fCentBin{-1},
      fZbin{-1},
      fFilterBit{32},
      fTPCNsigPrimaryPionCut{3},
      fTOFNsigPrimaryPionCut{-999},
      fPrimaryPionEtaCut{0.8},
      fPrimaryPionZVertexCut{2.0},
      fPrimaryPionXYVertexSigmaCut{3},
      fTPCNsigSecondaryPionCut{3},
      fTOFNsigSecondaryPionCut{-999},
      fSecondaryPionEtaCut{0.8},
      fSecondaryPionZVertexCut{2.0},
      fSecondaryPionXYVertexSigmaCut{3},
      fTPCNsigKaonCut{3},
      fTOFNsigKaonCut{3},
      fKaonEtaCut{0.8},
      fKaonZVertexCut{2.0},
      fKaonXYVertexSigmaCut{3},
      fSecondaryMassWindowCut{0.1},
      fSecondaryRapCut{0.5},
      fK1YCutHigh{0.5},
      fK1YCutLow{-0.5},
      fMinK1OA{0},
      fMaxK1OA{0.87},
      fMinPairAsym{-0.1},
      fMaxPairAsym{1},
      fMinK1MassCutAnotherScenario{0},
      fMaxK1MassCutAnotherScenario{999},
      fMinK1PiKa{0},
      fMaxK1PiKa{999},
      fMinSecondarypTCut{0.0},
      fMaxSecondarypTCut{999},
      fTracks{0}
{
  /// Default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskK1::AliAnalysisTaskK1(const char *name, Bool_t MCcase)
    : AliAnalysisTaskSE(name),
      fTrackCuts{nullptr},
      fPIDResponse{nullptr},
      fAODMCHeader{nullptr},
      fList{nullptr},
      fMCArray{nullptr},
      fVertex{nullptr},
      fNanoTree{nullptr},
      fNanoMCTree{nullptr},
      fRecK1daughter{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      fSimK1part{0, 0, 0, 0, 0, 0, 0, 0, 0},
      fHn5DK1Data{nullptr},
      fHn5DK1MC{nullptr},
      fHn2DEvtNorm{nullptr},
      hMultiplicity{nullptr},
      hEtaTrack_before{nullptr},
      hDCAPVTrack_before{nullptr},
      hDCArPVTrack_before{nullptr},
      hPtTrack_before{nullptr},
      hEtaTrack_ppion{nullptr},
      hDCAPVTrack_ppion{nullptr},
      hDCArPVTrack_ppion{nullptr},
      hPtTrack_ppion{nullptr},
      hTPCPIDTrack_ppion{nullptr},
      hTPCPIDTrackNsigVspT_ppion{nullptr},
      hEtaTrack_spion{nullptr},
      hDCAPVTrack_spion{nullptr},
      hDCArPVTrack_spion{nullptr},
      hPtTrack_spion{nullptr},
      hTPCPIDTrack_spion{nullptr},
      hTPCPIDTrackNsigVspT_spion{nullptr},
      hEtaTrack_kaon{nullptr},
      hDCAPVTrack_kaon{nullptr},
      hDCArPVTrack_kaon{nullptr},
      hPtTrack_kaon{nullptr},
      hTPCPIDTrack_kaon{nullptr},
      hTPCPIDTrackNsigVspT_kaon{nullptr},
      hK1OA{nullptr},
      hK1PairAsymm{nullptr},
      hInvMass_k892_rho{nullptr},
      hInvMass_Secondary_PiKa{nullptr},
      hK1SecondarypT{nullptr},
      hK1OA_cut{nullptr},
      hK1PairAsymm_cut{nullptr},
      hInvMass_k892_rho_cut{nullptr},
      hInvMass_Secondary_PiKa_cut{nullptr},
      hK1SecondarypT_cut{nullptr},
      hK1OA_MCTrue{nullptr},
      hK1PairAsymm_MCTrue{nullptr},
      hInvMass_k892_rho_MCTrue{nullptr},
      hInvMass_Secondary_PiKa_MCTrue{nullptr},
      hK1SecondarypT_MCTrue{nullptr},
      fIsAOD{false},
      fIsNano{false},
      fSetMixing{false},
      fFillQAPlot{true},
      fIsMC{MCcase},
      fIsPrimaryMC{true},
      fFillTree{false},
      fIsINEL{false},
      fIsHM{false},
      fSkipFillingHistogram{false},
      fSetModeK892orRho{true},
      fEMpool{},
      fBinCent{},
      fBinZ{},
      fPosPV{},
      fMagField{0},
      fCent{-1},
      fCustomEventID{0},
      fnMix{10},
      fCentBin{-1},
      fZbin{-1},
      fFilterBit{32},
      fTPCNsigPrimaryPionCut{3},
      fTOFNsigPrimaryPionCut{-999},
      fPrimaryPionEtaCut{0.8},
      fPrimaryPionZVertexCut{2.0},
      fPrimaryPionXYVertexSigmaCut{7.0},
      fTPCNsigSecondaryPionCut{3},
      fTOFNsigSecondaryPionCut{-999},
      fSecondaryPionEtaCut{0.8},
      fSecondaryPionZVertexCut{2.0},
      fSecondaryPionXYVertexSigmaCut{7.0},
      fTPCNsigKaonCut{3},
      fTOFNsigKaonCut{3},
      fKaonEtaCut{0.8},
      fKaonZVertexCut{2.0},
      fKaonXYVertexSigmaCut{7.0},
      fSecondaryMassWindowCut{0.1},
      fSecondaryRapCut{0.5},
      fK1YCutHigh{0.5},
      fK1YCutLow{-0.5},
      fMinK1OA{0},
      fMaxK1OA{3.14},
      fMinPairAsym{-1},
      fMaxPairAsym{1},
      fMinK1MassCutAnotherScenario{0.0},
      fMaxK1MassCutAnotherScenario{999},
      fMinK1PiKa{0},
      fMaxK1PiKa{999},
      fMinSecondarypTCut{0.0},
      fMaxSecondarypTCut{999},
      fTracks{0}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskK1::~AliAnalysisTaskK1()
{
  delete fNanoTree;
  delete fNanoMCTree;
}
//___________________________________________________________________
void AliAnalysisTaskK1::SetCutOpen()
{
  // Set the cut values for the K1 analysis to the loosest values
  SetMaxTPCnSigPrimaryPion(5);
  SetMaxTOFnSigPrimaryPion(-999);
  SetMaxEtaPrimaryPion(0.8);
  SetMaxVertexZPrimaryPion(3);
  SetMaxVertexXYsigPrimaryPion(10);

  SetMaxTPCnSigSecondaryPion(5);
  SetMaxTOFnSigTOFSecondaryPion(-999);
  SetMaxEtaSecondaryPion(0.8);
  SetMaxVertexZSecondaryPion(3);
  SetMaxVertexXYsigSecondaryPion(10);

  SetMaxTPCnSigKaon(5);
  SetMaxTOFnSigKaon(-999);
  SetMaxEtaKaon(0.8);
  SetMaxVertexZKaon(3);
  SetMaxVertexXYsigKaon(10);

  SetSecondaryMassWindow(1);
  SetSecondaryRapidityCut(1);

  SetK1RapidityCutHigh(1);
  SetK1RapidityCutLow(-1);
  SetK1OAMin(0);
  SetK1OAMax(999);
  SetK1PairAssymMin(-1);
  SetK1PairAssymMax(1);
  SetK1AnotherSecnarioMassCutMin(0);
  SetK1AnotherSecnarioMassCutMax(999);
  SetK1PiKaMassCutMin(0);
  SetK1PiKaMassCutMax(999);
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::UserCreateOutputObjects()
{
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  fList = new TList();
  fList->SetOwner(kTRUE);

  // Create THnSparse
  auto binAnti = AxisStr("AType", {"Normal", "Anti"});
  auto binType = (!fIsMC) ? AxisStr("Type", {"K1P", "K1N", "K1P_mix", "K1N_mix"})
                          : AxisStr("Type", {"K1P", "K1N", "K1P_mix", "K1N_mix", "K1P_no", "K1N_no"});
  auto binTypeMC = AxisStr("Type", {"K1P_gen",
                                    "K1N_gen",
                                    "K1P_gen_inel10",
                                    "K1N_gen_inel10",
                                    "K1P_gen_inel10_igz",
                                    "K1N_gen_inel10_igz",
                                    "K1P_gen_trig",
                                    "K1N_gen_trig",
                                    "K1P_rec",
                                    "K1N_rec"});
  std::vector<double> centaxisbin;
  (fIsHM) ? centaxisbin = {0, 0.001, 0.01, 0.05, 0.1}
          : centaxisbin = {-1, 0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100}; // can be use from pp to PbPb
  if (fIsMC)
    centaxisbin = {-1, 0, 0.001, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100}; // for general MC
  fBinCent = AxisVar("Cent", centaxisbin);
  auto binPt = AxisFix("Pt", 150, 0, 15);
  auto binMass = AxisFix("Mass", 700, 0.6, 2.0);
  auto binMassMC = AxisFix("Mass", 450, 0.6, 1.5);
  // fBinZ = AxisFix("Z", 20, -10, 10); // 1cm diff
  fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10}); // moderate diff

  fHn5DK1Data = CreateTHnSparse("K1_data", "K1_data", 5, {binAnti, binType, fBinCent, binPt, binMass}, "s");
  if (fIsMC)
  {
    auto binTypeMCNorm = AxisStr("Type", {"kAll", "kINEL10", "kINEL_trig", "kINEL_trig_vtx",
                                          "kINEL_trig_vtx10", "kINELg0", "kINELg010", "kINELg0_trig",
                                          "kINELg0_trig_vtx", "kINELg0_trig_vtx10", "kSelected"});
    fHn5DK1MC = CreateTHnSparse("K1_mc", "K1_mc", 5, {binAnti, binTypeMC, fBinCent, binPt, binMassMC}, "s");
    fHn2DEvtNorm = CreateTHnSparse("Normalisation", "", 2, {binTypeMCNorm, fBinCent}, "s");
  }
  if (fFillQAPlot)
  {
    // Make a new Folder named "QA" under the fList
    TList *QAList = new TList();
    QAList->SetOwner(kTRUE);
    QAList->SetName("QA");
    fList->Add(QAList);
    // All tracks
    hEtaTrack_before = new TH1D("hEtaTrack_before", "Eta distribution of tracks;#eta;counts", 20, -1.0, 1.0);
    hDCAPVTrack_before = new TH1D("hDCAPVTrack_before", "DCA to PV distribution of tracks;DCA (cm);counts", 30, 0, 3);
    hDCArPVTrack_before = new TH1D("hDCArPVTrack_before", "DCAr to PV distribution of tracks;DCAr (cm);counts", 50, 0, 0.5);
    hPtTrack_before = new TH1D("hPtTrack_before", "#it{p}_{T} distribution of tracks;#it{p}_{T} (GeV/c);counts", 150, 0, 15);
    QAList->Add(hEtaTrack_before);
    QAList->Add(hDCAPVTrack_before);
    QAList->Add(hDCArPVTrack_before);
    QAList->Add(hPtTrack_before);
    // Primary pions
    hEtaTrack_ppion = new TH1D("hEtaTrack_ppion", "Eta distribution of primary pions;#eta;counts", 20, -1.0, 1.0);
    hDCAPVTrack_ppion = new TH1D("hDCAPVTrack_ppion", "DCA to PV distribution of primary pions;DCA (cm);counts", 30, 0, 3);
    hDCArPVTrack_ppion = new TH1D("hDCArPVTrack_ppion", "DCAr to PV distribution of primary pions;DCAr (cm);counts", 50, 0, 0.5);
    hPtTrack_ppion = new TH1D("hPtTrack_ppion", "#it{p}_{T} distribution of primary pions;#it{p}_{T} (GeV/c);counts", 150, 0, 15);
    hTPCPIDTrack_ppion = new TH2D("hTPCPIDTrack_ppion", "TPC PID of primary pions;P (GeV/c);TPC signal", 200, 0, 20, 200, 0, 200);
    hTPCPIDTrackNsigVspT_ppion = new TH2D("hTPCPIDTrackNsigVspT_ppion", "TPC PID of primary pions;P (GeV/c);TPC n#sigma", 100, 0, 10, 200, -10, 10);
    QAList->Add(hEtaTrack_ppion);
    QAList->Add(hDCAPVTrack_ppion);
    QAList->Add(hDCArPVTrack_ppion);
    QAList->Add(hPtTrack_ppion);
    QAList->Add(hTPCPIDTrack_ppion);
    QAList->Add(hTPCPIDTrackNsigVspT_ppion);
    // Secondary pions
    hEtaTrack_spion = new TH1D("hEtaTrack_spion", "Eta distribution of secondary pions;#eta;counts", 20, -1.0, 1.0);
    hDCAPVTrack_spion = new TH1D("hDCAPVTrack_spion", "DCA to PV distribution of secondary pions;DCA (cm);counts", 30, 0, 3);
    hDCArPVTrack_spion = new TH1D("hDCArPVTrack_spion", "DCAr to PV distribution of secondary pions;DCAr (cm);counts", 50, 0, 0.5);
    hPtTrack_spion = new TH1D("hPtTrack_spion", "#it{p}_{T} distribution of secondary pions;#it{p}_{T} (GeV/c);counts", 150, 0, 15);
    hTPCPIDTrack_spion = new TH2D("hTPCPIDTrack_spion", "TPC PID of secondary pions;P (GeV/c);TPC signal", 200, 0, 20, 200, 0, 200);
    hTPCPIDTrackNsigVspT_spion = new TH2D("hTPCPIDTrackNsigVspT_spion", "TPC PID of secondary pions;P (GeV/c);TPC n#sigma", 100, 0, 10, 200, -10, 10);
    QAList->Add(hEtaTrack_spion);
    QAList->Add(hDCAPVTrack_spion);
    QAList->Add(hDCArPVTrack_spion);
    QAList->Add(hPtTrack_spion);
    QAList->Add(hTPCPIDTrack_spion);
    QAList->Add(hTPCPIDTrackNsigVspT_spion);
    // Kaons
    hEtaTrack_kaon = new TH1D("hEtaTrack_kaon", "Eta distribution of kaons;#eta;counts", 20, -1.0, 1.0);
    hDCAPVTrack_kaon = new TH1D("hDCAPVTrack_kaon", "DCA to PV distribution of kaons;DCA (cm);counts", 30, 0, 3);
    hDCArPVTrack_kaon = new TH1D("hDCArPVTrack_kaon", "DCAr to PV distribution of kaons;DCAr (cm);counts", 50, 0, 0.5);
    hPtTrack_kaon = new TH1D("hPtTrack_kaon", "#it{p}_{T} distribution of kaons;#it{p}_{T} (GeV/c);counts", 150, 0, 15);
    hTPCPIDTrack_kaon = new TH2D("hTPCPIDTrack_kaon", "TPC PID of kaons;P (GeV/c);TPC signal", 200, 0, 20, 200, 0, 200);
    hTPCPIDTrackNsigVspT_kaon = new TH2D("hTPCPIDTrackNsigVspT_kaon", "TPC PID of kaons;P (GeV/c);TPC n#sigma", 100, 0, 10, 200, -10, 10);
    QAList->Add(hEtaTrack_kaon);
    QAList->Add(hDCAPVTrack_kaon);
    QAList->Add(hDCArPVTrack_kaon);
    QAList->Add(hPtTrack_kaon);
    QAList->Add(hTPCPIDTrack_kaon);
    QAList->Add(hTPCPIDTrackNsigVspT_kaon);
    // K1
    hK1OA = new TH1D("hK1OA", "Opening angle distribution of K1;Opening angle (rad);counts", 100, 0, 3.14);
    hK1PairAsymm = new TH1D("hK1PairAsymm", "Pair asymmetry distribution of K1;Asymmetry;counts", 100, -1, 1);
    hInvMass_k892_rho = new TH2D("hInvMass_k892_rho", "Invariant mass distribution k892 invmass vs rho inv mass;Mass_{K*} (GeV/c^{2});Mass_{#rho} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
    hInvMass_Secondary_PiKa = new TH2D("hInvMass_Secondary_PiKa", "Invariant mass distribution secondary invmass vs PiKa inv mass;Mass_{secondary} (GeV/c^{2});Mass_{K#pi} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
    hK1SecondarypT = new TH1D("hK1SecondarypT", "Secondary pT distribution of K1;#it{p}_{T} (GeV/c);counts", 100, 0, 10);

    hK1OA_cut = new TH1D("hK1OA_cut", "Opening angle distribution of K1;Opening angle (rad);counts", 100, 0, 3.14);
    hK1PairAsymm_cut = new TH1D("hK1PairAsymm_cut", "Pair asymmetry distribution of K1;Asymmetry;counts", 100, -1, 1);
    hInvMass_k892_rho_cut = new TH2D("hInvMass_k892_rho_cut", "Invariant mass distribution k892 invmass vs rho inv mass;Mass_{K*} (GeV/c^{2});Mass_{#rho} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
    hInvMass_Secondary_PiKa_cut = new TH2D("hInvMass_Secondary_PiKa_cut", "Invariant mass distribution k892 invmass vs rho inv mass;Mass_{secondary} (GeV/c^{2});Mass_{K#pi} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
    hK1SecondarypT_cut = new TH1D("hK1SecondarypT_cut", "Secondary pT distribution of K1;#it{p}_{T} (GeV/c);counts", 100, 0, 10);
    QAList->Add(hK1OA);
    QAList->Add(hK1PairAsymm);
    QAList->Add(hInvMass_k892_rho);
    QAList->Add(hInvMass_Secondary_PiKa);
    QAList->Add(hK1SecondarypT);
    QAList->Add(hK1OA_cut);
    QAList->Add(hK1PairAsymm_cut);
    QAList->Add(hInvMass_k892_rho_cut);
    QAList->Add(hInvMass_Secondary_PiKa_cut);
    QAList->Add(hK1SecondarypT_cut);
    if (fIsMC)
    {
      hK1OA_MCTrue = new TH1D("hK1OA_MCTrue", "Opening angle distribution of K1;Opening angle (rad);counts", 100, 0, 3.14);
      hK1PairAsymm_MCTrue = new TH1D("hK1PairAsymm_MCTrue", "Pair asymmetry distribution of K1;Asymmetry;counts", 100, -1, 1);
      hInvMass_k892_rho_MCTrue = new TH2D("hInvMass_k892_rho_MCTrue", "Invariant mass distribution k892 invmass vs rho inv mass;Mass_{K*} (GeV/c^{2});Mass_{#rho} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
      hInvMass_Secondary_PiKa_MCTrue = new TH2D("hInvMass_Secondary_PiKa_MCTrue", "Invariant mass distribution k892 invmass vs rho inv mass;Mass_{secondary} (GeV/c^{2});Mass_{K#pi} (GeV/c^{2})", 100, 0.6, 2.0, 200, 0, 2);
      hK1SecondarypT_MCTrue = new TH1D("hK1SecondarypT_MCTrue", "Secondary pT distribution of K1;#it{p}_{T} (GeV/c);counts", 100, 0, 10);
      QAList->Add(hK1OA_MCTrue);
      QAList->Add(hK1PairAsymm_MCTrue);
      QAList->Add(hInvMass_k892_rho_MCTrue);
      QAList->Add(hInvMass_Secondary_PiKa_MCTrue);
      QAList->Add(hK1SecondarypT_MCTrue);
    }
  }
  if (fIsHM)
    hMultiplicity = new TH1F("hMultiplicity", "", 100, 0, 0.1);
  else
    hMultiplicity = new TH1F("hMultiplicity", "", 101, -1, 100);

  fEMpool.resize(fBinCent.GetNbins() + 1, std::vector<eventpool>(fBinZ.GetNbins() + 1));
  if (!fSkipFillingHistogram)
  {
    fList->Add(fHn5DK1Data);
    fList->Add(fHn5DK1MC);
    fList->Add(fHn2DEvtNorm);
    fList->Add(hMultiplicity);
    fEventCut.AddQAplotsToList(fList); // QA histograms from AliEventCuts
    PostData(1, fList);
  }
  if (fFillTree)
  {
    OpenFile(1);
    fNanoTree = new TTree("RTree", "K1 daughter tree");
    fNanoTree->Branch("RK1Daughter", &fRecK1daughter);
    PostData(2, fNanoTree);

    if (fIsMC)
    {
      fNanoMCTree = new TTree("STree", "K1 MC tree");
      fNanoMCTree->Branch("SK1Daughter", &fSimK1part);
      PostData(3, fNanoMCTree);
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (!event) // if there is no event, return
  {
    PostAllData();
    AliInfo("Could not retrieve event");
    return;
  }

  // ESD? AOD?
  event->IsA() == AliESDEvent::Class() ? fEvt = dynamic_cast<AliESDEvent *>(event)  // ESD Case
                                       : fEvt = dynamic_cast<AliAODEvent *>(event); // AOD Case

  if (!fIsAOD && (event->IsA() != AliESDEvent::Class()))
    fIsAOD = true; // check if it is AOD or ESD
  if (!fEvt)       // if there is no event, return
  {
    PostAllData();
    return;
  }
  AliInputEventHandler *inputHandler = (AliInputEventHandler *)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  bool IsEvtSelected{false}, IsINEL0True{false}, IsSelectedTrig{false}, IsVtxInZCut{false}, IsGoodVertex{false};

  AliNanoAODHeader *nanoHeader = dynamic_cast<AliNanoAODHeader *>(fInputEvent->GetHeader());
  if (!nanoHeader) // if it's not a nanoAOD, use AliEventCuts
  {
    IsEvtSelected = fEventCut.AcceptEvent(event);
    if (fIsMC)
    {
      if (fIsAOD)
      {
        fMCArray = (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
        fAODMCHeader = (AliAODMCHeader *)fEvt->FindListObject(AliAODMCHeader::StdBranchName());
      }
      fMCEvent = MCEvent();
      IsINEL0True = fEventCut.IsTrueINELgtZero(fEvt, true);
    }
    fCent = AliMultSelectionTask::IsINELgtZERO(event) ? fEventCut.GetCentrality() : -0.5; // for INEL case, we save -0.5
    fPIDResponse = (AliPIDResponse *)inputHandler->GetPIDResponse();
    if (!fPIDResponse)
      AliInfo("No PIDd");
    IsVtxInZCut = fEventCut.PassedCut(AliEventCuts::kVertexPosition);
    IsSelectedTrig = fEventCut.PassedCut(AliEventCuts::kTrigger);
    IsGoodVertex = fEventCut.PassedCut(AliEventCuts::kVertexQuality);
  }
  else // if it's a nanoAOD, use NanoAOD header (no event selection needed)
  {
    if (!fIsNano)
      fIsNano = kTRUE;
    if (fIsMC)
    {
      if (fIsAOD)
        fMCArray = (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
      fMCEvent = MCEvent();
      IsINEL0True = true;
    }
    IsEvtSelected = true; // no event selection needed
    fCent = nanoHeader->GetCentr("V0M");
    static int inel_index = -1;
    if (inel_index < 0)
      inel_index = nanoHeader->GetVarIndex("cstINELgt0");
    if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
      fCent = -0.5;
  }

  if (fIsMC && !fSkipFillingHistogram)
  {
    // Fill Normalisation histogram
    if (IsVtxInZCut)
      FillMCinput(fMCEvent, 1);
    if (IsINEL0True && IsVtxInZCut)
      FillMCinput(fMCEvent, 2);
    if (IsSelectedTrig)
      FillMCinput(fMCEvent, 3);
    FillTHnSparse(fHn2DEvtNorm, {(int)kAll, (double)fCent});
    if (IsINEL0True)
    {
      FillTHnSparse(fHn2DEvtNorm, {(int)kINELg0, (double)fCent});
      if (IsVtxInZCut)
      {
        FillTHnSparse(fHn2DEvtNorm, {(int)kINELg010, (double)fCent});
      }
      if (IsSelectedTrig)
      {
        FillTHnSparse(fHn2DEvtNorm, {(int)kINELg0_trig, (double)fCent});
        if (IsGoodVertex)
        {
          FillTHnSparse(fHn2DEvtNorm, {(int)kINELg0_trig_vtx, (double)fCent});
          if (IsVtxInZCut)
          {
            FillTHnSparse(fHn2DEvtNorm, {(int)kINELg0_trig_vtx10, (double)fCent});
          }
        }
      }
    }
    if (IsVtxInZCut)
    {
      FillTHnSparse(fHn2DEvtNorm, {(int)kINEL10, (double)fCent});
    }
    if (IsSelectedTrig)
    {
      FillTHnSparse(fHn2DEvtNorm, {(int)kINEL_trig, (double)fCent});
      if (IsGoodVertex)
      {
        FillTHnSparse(fHn2DEvtNorm, {(int)kINEL_trig_vtx, (double)fCent});
        if (IsVtxInZCut)
        {
          FillTHnSparse(fHn2DEvtNorm, {(int)kINEL_trig_vtx10, (double)fCent});
        }
      }
    }
  }

  if (!IsEvtSelected)
  {
    PostAllData();
    return;
  }

  if (!fSkipFillingHistogram)
    hMultiplicity->Fill(fCent); // multiplicity distribution for basic event QA
  fCustomEventID = (fIsNano) ? nanoHeader->GetRunNumberIndex() + fEntry : ((unsigned long)(event->GetBunchCrossNumber()) << 32) + event->GetTimeStamp();

  if (fIsMC)
  {
    FillMCinput(fMCEvent, 0);
    FillTHnSparse(fHn2DEvtNorm, {(int)kSelected, (double)fCent});
  }
  if (fIsAOD)
    fVertex = ((AliAODEvent *)fEvt)->GetPrimaryVertex();
  const AliVVertex *pVtx = fEvt->GetPrimaryVertex();
  fPosPV[0] = pVtx->GetX();
  fPosPV[1] = pVtx->GetY();
  fPosPV[2] = pVtx->GetZ();
  fMagField = fEvt->GetMagneticField();

  // Event Mixing pool -----------------------------------------------------
  fZbin = fBinZ.FindBin(fPosPV[2]) - 1;   // Event mixing z-bin
  fCentBin = fBinCent.FindBin(fCent) - 1; // Event mixing cent bin
  if (fIsINEL)
    fCentBin = 0; // for INEL case

  bool checkTrackPools = FillTrackPools();

  if (checkTrackPools)
  {
    if (!fSkipFillingHistogram)
      FillHistograms(); // Fill the histogram
  }
  if (fSetMixing && fGoodPrimaryPionArray.size())
  {
    FillTrackToEventPool(); // use only pion track pool.
  }

  PostAllData();
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::Terminate(Option_t *)
{
  std::cout << "Total processed tracks: " << fTracks << std::endl;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskK1::FillTrackPools()
{
  // Fill the track pool: primary pion, secondary pion, kaon at the same time
  const UInt_t nTracks = fEvt->GetNumberOfTracks();
  Int_t ntrack = 0;
  fGoodPrimaryPionArray.clear();
  fGoodSecondaryPionArray.clear();
  fGoodKaonArray.clear();
  fGoodTracksArray.clear();
  fMCTrueTrack.clear();
  AliVTrack *track = nullptr;
  Float_t bCov[3];
  Double_t nTPCNSigPion{0}, nTPCNSigKaon{0}, nTOFNSigPion{0}, nTOFNSigKaon{0}, lDCAz{0}, lpT{0}, lsigmaDCAr{0}, lDCAr{0}, lEta{0}, lEnergy{0}, lStatus{0};
  Bool_t isPassPrimaryPionSelection{false}, isPassSecondaryPionSelection{false}, isPassKaonSelection{false};
  Bool_t isAcceptedTrack{false};
  Int_t lMotherID{-999}, lPDG{-999}, trackSelResult{-999};
  fTracks += nTracks;
  for (UInt_t it = 0; it < nTracks; it++)
  {
    isPassPrimaryPionSelection = true;
    isPassSecondaryPionSelection = true;
    isPassKaonSelection = true;
    track = (AliVTrack *)fEvt->GetTrack(it);
    if (!track)
      continue;

    if (!IsTrackAccepted(track))
      continue;
    auto [lDCAr, lDCAz, bCov] = GetImpactParametersWrapper(track);

    nTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
    nTPCNSigKaon = GetTPCnSigma(track, AliPID::kKaon);
    nTOFNSigPion = GetTOFnSigma(track, AliPID::kPion);
    nTOFNSigKaon = GetTOFnSigma(track, AliPID::kKaon);

    // if the particle is not pion or kaon in 5 sigma, skip the track
    if (TMath::Abs(nTPCNSigPion) > 5 && TMath::Abs(nTPCNSigKaon) > 5)
      continue;

    lpT = track->Pt();
    if (lpT < 0.15) // minimum pT cut
      continue;
    lsigmaDCAr = (0.0026 + 0.0050 / lpT);
    lEta = track->Eta();

    if (fFillQAPlot)
    {
      // Fill the QA plot for the whole track. we don't know which track is it yet.
      hEtaTrack_before->Fill(lEta);
      hDCAPVTrack_before->Fill(lDCAz);
      hDCArPVTrack_before->Fill(lDCAr);
      hPtTrack_before->Fill(lpT);
    }

    trackSelResult = TrackSelection(track, nTPCNSigPion, nTOFNSigPion, nTPCNSigKaon, nTOFNSigKaon, lDCAz, lDCAr, lsigmaDCAr, lEta);
    if (trackSelResult == 0)
      continue;

    isPassPrimaryPionSelection = (trackSelResult & kPrimaryPionPass);
    isPassSecondaryPionSelection = (trackSelResult & kSecondaryPionPass);
    isPassKaonSelection = (trackSelResult & kKaonPass);

    if (fFillQAPlot)
    {
      if (isPassPrimaryPionSelection)
      {
        hEtaTrack_ppion->Fill(lEta);
        hDCAPVTrack_ppion->Fill(lDCAz);
        hDCArPVTrack_ppion->Fill(lDCAr);
        hPtTrack_ppion->Fill(lpT);
        hTPCPIDTrack_ppion->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
        hTPCPIDTrackNsigVspT_ppion->Fill(lpT, nTPCNSigPion);
      }
      if (isPassSecondaryPionSelection)
      {
        hEtaTrack_spion->Fill(lEta);
        hDCAPVTrack_spion->Fill(lDCAz);
        hDCArPVTrack_spion->Fill(lDCAr);
        hPtTrack_spion->Fill(lpT);
        hTPCPIDTrack_spion->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
        hTPCPIDTrackNsigVspT_spion->Fill(lpT, nTPCNSigPion);
      }
      if (isPassKaonSelection)
      {
        hEtaTrack_kaon->Fill(lEta);
        hDCAPVTrack_kaon->Fill(lDCAz);
        hDCArPVTrack_kaon->Fill(lDCAr);
        hPtTrack_kaon->Fill(lpT);
        hTPCPIDTrack_kaon->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
        hTPCPIDTrackNsigVspT_kaon->Fill(lpT, nTPCNSigKaon);
      }
    }

    if (isPassPrimaryPionSelection)
      fGoodPrimaryPionArray.push_back(it);
    if (isPassSecondaryPionSelection)
      fGoodSecondaryPionArray.push_back(it);
    if (isPassKaonSelection)
      fGoodKaonArray.push_back(it);

    if (fFillTree)
    {
      if (isPassPrimaryPionSelection || isPassSecondaryPionSelection || isPassKaonSelection)
      {
        fRecK1daughter.eventID = fCustomEventID;
        fRecK1daughter.id = track->GetID();
        fRecK1daughter.zVertex = fPosPV[2];
        fRecK1daughter.px = track->Px();
        fRecK1daughter.py = track->Py();
        fRecK1daughter.pz = track->Pz();
        fRecK1daughter.eta = track->Eta();
        fRecK1daughter.dcaxy = lDCAr;
        fRecK1daughter.dcaz = lDCAz;
        fRecK1daughter.tpcNsigmaPi = nTPCNSigPion;
        fRecK1daughter.tofNsigmaPi = nTOFNSigPion;
        fRecK1daughter.tpcNsigmaKa = nTPCNSigKaon;
        fRecK1daughter.tofNsigmaKa = nTOFNSigKaon;
        fRecK1daughter.centrality = fCent;
        fRecK1daughter.flag = 0;
        fRecK1daughter.flag |= (track->Charge() > 0) ? RK1Daughter::kPositive : 0;
        fRecK1daughter.flag |= (nTOFNSigPion > -998) ? RK1Daughter::kHasTOF : 0;
        if (fIsMC)
        {
          if (!fIsAOD)
          {
            fRecK1daughter.motherID = static_cast<AliMCParticle *>(fMCEvent->GetTrack(track->GetLabel()))->GetMother();
            fRecK1daughter.pdg = static_cast<AliMCParticle *>(fMCEvent->GetTrack(track->GetLabel()))->PdgCode();
          }
          else
          {
            fRecK1daughter.motherID = static_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(track->GetLabel())))->GetMother();
            fRecK1daughter.pdg = static_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(track->GetLabel())))->GetPdgCode();
          }
        }
        fNanoTree->Fill();
      }
    }
  }
  if (fGoodPrimaryPionArray.size() > 0 && fGoodSecondaryPionArray.size() > 0 && fGoodKaonArray.size() > 0)
    return true;
  else
    return false;
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::FillHistograms()
{
  AliVTrack *track_primiary_pion, *track_secondary_pion, *track_kaon, *track_primiary_pion_mix;
  Bool_t isAnti, isPirmaryPionPlus;
  Bool_t SkipMixing = kFALSE;
  Int_t primiaryID{0}, secondaryID{0}, kaonID{0}, primiaryMixID{0};
  int sign = kAllType;
  int binAnti = 0;
  Double_t lMassK892{0}, lMassRho{0}, lMassPiKa{0}, lpTK892{0}, lpTRho{0}, lSecondaryRap{0}, lSecondaryMass{0}, lSecondarypT{0};
  Double_t lK1Mass{0}, lK1Rapidity{0}, lK1Pt{0}, lK1Angle{0}, lK1PairAsym{0};
  Double_t lSecondaryMassRef = (fSetModeK892orRho) ? kK892Mass : kRhoMass;

  TLorentzVector vecPrimaryPion, vecSecondaryPion, vecKaon, vecPrimaryPionMix;
  TLorentzVector vecK892, vecK1, tempRho, tempPiKa;
  const UInt_t nPrimaryPions = fGoodPrimaryPionArray.size();
  const UInt_t nSecondaryPions = fGoodSecondaryPionArray.size();
  const UInt_t nKaons = fGoodKaonArray.size();

  // Prepare the track pool for mixing
  tracklist trackpool;
  if (fSetMixing)
  {
    eventpool &ep = fEMpool[fCentBin][fZbin];
    if ((int)ep.size() < (int)fnMix)
      SkipMixing = kTRUE;
    if (!SkipMixing)
    {
      for (auto pool : ep)
      {
        for (auto track : pool)
          trackpool.push_back((AliVTrack *)track);
      }
    }
  }

  // Reconstruction loop (MAIN LOOP)
  for (UInt_t iSecondaryPion = 0; iSecondaryPion < nSecondaryPions; iSecondaryPion++) // secondary pion which is used for K(892)
  {
    track_secondary_pion = (AliVTrack *)fEvt->GetTrack(fGoodSecondaryPionArray[iSecondaryPion]);
    if (!track_secondary_pion)
      continue;
    secondaryID = track_secondary_pion->GetID();
    vecSecondaryPion.SetXYZM(track_secondary_pion->Px(), track_secondary_pion->Py(), track_secondary_pion->Pz(), kPionMass);

    for (UInt_t iKaon = 0; iKaon < nKaons; iKaon++) // kaon which is used for K(892)
    {
      track_kaon = (AliVTrack *)fEvt->GetTrack(fGoodKaonArray[iKaon]);
      if (!track_kaon)
        continue;
      kaonID = track_kaon->GetID();
      if (kaonID == secondaryID)
        continue;

      vecKaon.SetXYZM(track_kaon->Px(), track_kaon->Py(), track_kaon->Pz(), kKaonMass);

      // K(892) reconstruction
      vecK892 = vecSecondaryPion + vecKaon;
      if (track_secondary_pion->Charge() > 0)
        isAnti = true;
      else
        isAnti = false;

      for (UInt_t iPrimaryPion = 0; iPrimaryPion < nPrimaryPions; iPrimaryPion++) // primary pion which is used for K1
      {
        track_primiary_pion = (AliVTrack *)fEvt->GetTrack(fGoodPrimaryPionArray[iPrimaryPion]);
        if (!track_primiary_pion)
          continue;
        primiaryID = track_primiary_pion->GetID();
        if (primiaryID == secondaryID || primiaryID == kaonID)
          continue;

        vecPrimaryPion.SetXYZM(track_primiary_pion->Px(), track_primiary_pion->Py(), track_primiary_pion->Pz(), kPionMass);

        vecK1 = vecK892 + vecPrimaryPion;
        // Y cut of K1 in any case
        lK1Rapidity = vecK1.Rapidity();
        if ((lK1Rapidity > fK1YCutHigh) || (lK1Rapidity < fK1YCutLow))
          continue;

        // Robust selection for both K892 and Rho daughter cases.
        // A: Pion1 (primary pion), B: Pion2 (secondary pion), C: Kaon (kaon)
        // If we want to apply K892 scenario, K892 will be B+C, and A will be primary pion
        // If we want to apply Rho scenario, Rho will be A+B, and C will be primary kaon
        // In anycase, A+C will be a remain combination. (we called it PiKa here.)
        tempRho = vecPrimaryPion + vecSecondaryPion; // Rho scenario
        tempPiKa = vecPrimaryPion + vecKaon;         // Remain combination
        lMassK892 = vecK892.M();
        lMassRho = tempRho.M();
        lMassPiKa = tempPiKa.M();
        lpTK892 = vecK892.Pt();
        lpTRho = tempRho.Pt();
        lSecondaryRap = (fSetModeK892orRho) ? vecK892.Rapidity() : tempRho.Rapidity();
        lSecondaryMass = (fSetModeK892orRho) ? lMassK892 : lMassRho;
        lSecondarypT = (fSetModeK892orRho) ? lpTK892 : lpTRho;

        // OA Cut between K1 daughters
        lK1Angle = (fSetModeK892orRho) ? vecK892.Angle(vecPrimaryPion.Vect()) : tempRho.Angle(vecKaon.Vect());
        // Pair Assym
        lK1PairAsym = (fSetModeK892orRho) ? (vecK892.E() - vecPrimaryPion.E()) / (vecK892.E() + vecPrimaryPion.E()) : (tempRho.E() - vecKaon.E()) / (tempRho.E() + vecKaon.E());

        if (fFillQAPlot)
        {
          hK1OA->Fill(lK1Angle);
          hK1PairAsymm->Fill(lK1PairAsym);
          hInvMass_k892_rho->Fill(lMassK892, lMassRho);
          hInvMass_Secondary_PiKa->Fill(lSecondaryMass, lMassPiKa);
          hK1SecondarypT->Fill(lSecondarypT);
        }
        // Apply Secondary selection cut here.
        // Y cut
        if (abs(lSecondaryRap) > fSecondaryRapCut)
          continue;
        // mass window
        if (abs(lSecondaryMass - lSecondaryMassRef) > fSecondaryMassWindowCut)
          continue;
        // Selection of unselected track pairs (another scenario)
        if ((lMassPiKa > fMaxK1MassCutAnotherScenario) || (lMassPiKa < fMinK1MassCutAnotherScenario))
          continue;
        // Selection of PiKa pair
        if ((lMassPiKa > fMaxK1PiKa) || (lMassPiKa < fMinK1PiKa))
          continue;
        // K1 OA cut
        if ((lK1Angle > fMaxK1OA) || (lK1Angle < fMinK1OA))
          continue;
        // K1 Pair Asymmetry cut
        if ((lK1PairAsym > fMaxPairAsym) || (lK1PairAsym < fMinPairAsym))
          continue;
        // Secondary pT cut
        if ((lSecondarypT < fMinSecondarypTCut) || (lSecondarypT > fMaxSecondarypTCut))
          continue;
        if (fFillQAPlot)
        {
          hK1OA_cut->Fill(lK1Angle);
          hK1PairAsymm_cut->Fill(lK1PairAsym);
          hInvMass_k892_rho_cut->Fill(lMassK892, lMassRho);
          hInvMass_Secondary_PiKa_cut->Fill(lSecondaryMass, lMassPiKa);
          hK1SecondarypT_cut->Fill(lSecondarypT);
        }
        if (track_primiary_pion->Charge() > 0)
          isPirmaryPionPlus = true;
        else
          isPirmaryPionPlus = false;

        (isAnti) ? binAnti = kAnti : binAnti = kNormal;
        (isPirmaryPionPlus) ? sign = kK1P : sign = kK1N;
        lK1Mass = vecK1.M();
        lK1Pt = vecK1.Pt();

        FillTHnSparse(fHn5DK1Data, {(double)binAnti, (double)sign, (double)fCent, lK1Pt, lK1Mass});

        if (fIsMC)
        {
          if (IsTrueK1(primiaryID, secondaryID, kaonID))
          {
            (isPirmaryPionPlus) ? sign = kK1P_REC : sign = kK1N_REC;
            FillTHnSparse(fHn5DK1MC, {(double)binAnti, (double)sign, (double)fCent, lK1Pt, lK1Mass});
            if (fFillQAPlot)
            {
              hK1OA_MCTrue->Fill(lK1Angle);
              hK1PairAsymm_MCTrue->Fill(lK1PairAsym);
              hInvMass_k892_rho_MCTrue->Fill(lMassK892, lMassRho);
              hInvMass_Secondary_PiKa_MCTrue->Fill(lSecondaryMass, lMassPiKa);
              hK1SecondarypT_MCTrue->Fill(lSecondarypT);
            }
          }
          else
          {
            // MC not true bkg
            (isPirmaryPionPlus) ? sign = kK1P_NOT : sign = kK1N_NOT;
            FillTHnSparse(fHn5DK1MC, {(double)binAnti, (double)sign, (double)fCent, lK1Pt, lK1Mass});
          }
        }
      } // primary pion loop

      if (fSetMixing && !SkipMixing && (fCentBin >= 0) && (fZbin >= 0))
      {
        for (UInt_t iPrimaryPionMix = 0; iPrimaryPionMix < trackpool.size(); iPrimaryPionMix++)
        {
          track_primiary_pion_mix = trackpool.at(iPrimaryPionMix);
          primiaryMixID = track_primiary_pion_mix->GetID();
          if (primiaryMixID == secondaryID || primiaryMixID == kaonID)
            continue;
          vecPrimaryPionMix.SetXYZM(track_primiary_pion_mix->Px(), track_primiary_pion_mix->Py(), track_primiary_pion_mix->Pz(), kPionMass);
          vecK1 = vecK892 + vecPrimaryPionMix;
          // Y cut
          lK1Rapidity = vecK1.Rapidity();
          if ((lK1Rapidity > fK1YCutHigh) || (lK1Rapidity < fK1YCutLow))
            continue;
          // Robust selection for both K892 and Rho daughter cases.
          tempRho = vecPrimaryPionMix + vecSecondaryPion; // Rho scenario
          tempPiKa = vecPrimaryPionMix + vecKaon;         // Remain combination
          lMassK892 = vecK892.M();
          lMassRho = tempRho.M();
          lMassPiKa = tempPiKa.M();
          lpTK892 = vecK892.Pt();
          lpTRho = tempRho.Pt();
          lSecondaryRap = (fSetModeK892orRho) ? vecK892.Rapidity() : tempRho.Rapidity();
          lSecondaryMass = (fSetModeK892orRho) ? lMassK892 : lMassRho;
          lSecondarypT = (fSetModeK892orRho) ? lpTK892 : lpTRho;

          // OA Cut between K1 daughters
          lK1Angle = (fSetModeK892orRho) ? vecK892.Angle(vecPrimaryPionMix.Vect()) : tempRho.Angle(vecKaon.Vect());
          // Pair Assym
          lK1PairAsym = (fSetModeK892orRho) ? (vecK892.E() - vecPrimaryPionMix.E()) / (vecK892.E() + vecPrimaryPionMix.E()) : (tempRho.E() - vecKaon.E()) / (tempRho.E() + vecKaon.E());

          // Apply Secondary selection cut here.
          // Y cut
          if (abs(lSecondaryRap) > fSecondaryRapCut)
            continue;
          // mass window
          if (abs(lSecondaryMass - lSecondaryMassRef) > fSecondaryMassWindowCut)
            continue;
          // Selection of unselected track pairs (another scenario)
          if ((lMassPiKa > fMaxK1MassCutAnotherScenario) || (lMassPiKa < fMinK1MassCutAnotherScenario))
            continue;
          // Selection of PiKa pair
          if ((lMassPiKa > fMaxK1PiKa) || (lMassPiKa < fMinK1PiKa))
            continue;
          // K1 OA cut
          if ((lK1Angle > fMaxK1OA) || (lK1Angle < fMinK1OA))
            continue;
          // K1 Pair Asymmetry cut
          if ((lK1PairAsym > fMaxPairAsym) || (lK1PairAsym < fMinPairAsym))
            continue;
          // Secondary pT cut
          if ((lSecondarypT < fMinSecondarypTCut) || (lSecondarypT > fMaxSecondarypTCut))
            continue;

          if (track_primiary_pion_mix->Charge() > 0)
            isPirmaryPionPlus = true;
          else
            isPirmaryPionPlus = false;

          (isAnti) ? binAnti = kAnti : binAnti = kNormal;
          (isPirmaryPionPlus) ? sign = kK1P_MIX : sign = kK1N_MIX;

          FillTHnSparse(fHn5DK1Data, {(double)binAnti, (double)sign, (double)fCent, vecK1.Pt(), vecK1.M()});
        } // mixing track loop
      }   // mixing

    } // kaon loop
  }   // secondary pion loop
}
void AliAnalysisTaskK1::FillMCinput(AliMCEvent *fMCEvent, int Fillbin)
{
  int sign = kAllType;
  int binAnti = 0;
  bool isSecondaryDaughter{false}, isK1{false};
  int kDaugtherCode = (fSetModeK892orRho) ? kK892Code : kRhoCode;

  if (!fIsAOD)
  {
    for (Int_t it = 0; it < fMCEvent->GetNumberOfPrimaries(); it++)
    {
      TParticle *mcInputTrack = (TParticle *)fMCEvent->GetTrack(it)->Particle();
      if (!mcInputTrack)
      {
        Error("UserExec", "Could not receive MC track %d", it);
        continue;
      }
      isSecondaryDaughter = false;
      isK1 = false;

      Int_t k1PdgCode = mcInputTrack->GetPdgCode();
      if (TMath::Abs(k1PdgCode) == kK1PCode)
        isK1 = true;
      if (TMath::Abs(k1PdgCode) == kDaugtherCode)
        isSecondaryDaughter = true;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;

      if (isK1) // Fill histogram for K1
      {
        // Y cut
        if ((mcInputTrack->Y() > fK1YCutHigh) || (mcInputTrack->Y() < fK1YCutLow))
          continue;

        (k1PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
        if (k1PdgCode > 0)
          sign = kK1P_GEN + (int)Fillbin * 2;
        else
          sign = kK1N_GEN + (int)Fillbin * 2;
        FillTHnSparse(fHn5DK1MC, {(double)binAnti, (double)sign, (double)fCent, mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
      }
      if (Fillbin > 0)
        continue;
      if (!isSecondaryDaughter && !isK1)
        continue;
      if (fFillTree) // Fill both K892 and K1
      {
        fSimK1part.eventID = fCustomEventID;
        fSimK1part.zVertex = fPosPV[2];
        fSimK1part.centrality = fCent;
        fSimK1part.id = it;
        fSimK1part.motherID = mcInputTrack->GetMother(0);
        fSimK1part.daughter1 = mcInputTrack->GetDaughter(0);
        fSimK1part.daughter2 = mcInputTrack->GetDaughter(1);
        fSimK1part.pt = mcInputTrack->Pt();
        fSimK1part.eta = mcInputTrack->Eta();
        fSimK1part.pdg = mcInputTrack->GetPdgCode();
        if (mcInputTrack->IsPrimary())
          fSimK1part.flag = SK1Daughter::kPrimary;
        else
          fSimK1part.flag = SK1Daughter::kSecondary;
        fNanoMCTree->Fill();
      }
    }
  }
  else
  {
    for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++)
    {
      AliAODMCParticle *mcInputTrack = (AliAODMCParticle *)fMCArray->At(it);
      if (!mcInputTrack)
      {
        Error("UserExec", "Could not receive MC track %d", it);
        continue;
      }
      isSecondaryDaughter = false;
      isK1 = false;

      Int_t k1PdgCode = mcInputTrack->GetPdgCode();

      if (TMath::Abs(k1PdgCode) == kK1PCode)
        isK1 = true;
      if (TMath::Abs(k1PdgCode) == kDaugtherCode)
        isSecondaryDaughter = true;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;
      if (isK1)
      {
        // Y cut
        if ((mcInputTrack->Y() > fK1YCutHigh) || (mcInputTrack->Y() < fK1YCutLow))
          continue;
        (k1PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
        if (k1PdgCode > 0)
          sign = kK1P_GEN + (int)Fillbin * 2;
        else
          sign = kK1N_GEN + (int)Fillbin * 2;
        FillTHnSparse(fHn5DK1MC, {(double)binAnti, (double)sign, (double)fCent, mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
      }
      if (Fillbin > 0)
        continue;
      if (!isSecondaryDaughter && !isK1)
        continue;
      if (fFillTree)
      {
        fSimK1part.eventID = fCustomEventID;
        fSimK1part.zVertex = fPosPV[2];
        fSimK1part.centrality = fCent;
        fSimK1part.id = it;
        fSimK1part.motherID = mcInputTrack->GetMother();
        fSimK1part.daughter1 = mcInputTrack->GetDaughterLabel(0);
        fSimK1part.daughter2 = mcInputTrack->GetDaughterLabel(1);
        fSimK1part.pt = mcInputTrack->Pt();
        fSimK1part.eta = mcInputTrack->Eta();
        fSimK1part.pdg = mcInputTrack->GetPdgCode();
        if (mcInputTrack->IsPrimary())
          fSimK1part.flag = SK1Daughter::kPrimary;
        else
          fSimK1part.flag = SK1Daughter::kSecondary;
        fNanoMCTree->Fill();
      }
    }
  }
}
Bool_t AliAnalysisTaskK1::IsTrueK1(UInt_t primiaryID, UInt_t secondaryID, UInt_t kaonID)
{
  AliVTrack *track_primiary_pion, *track_secondary_pion, *track_kaon;
  track_primiary_pion = (AliVTrack *)fEvt->GetTrack(primiaryID);
  track_secondary_pion = (AliVTrack *)fEvt->GetTrack(secondaryID);
  track_kaon = (AliVTrack *)fEvt->GetTrack(kaonID);
  if (!track_primiary_pion || !track_secondary_pion || !track_kaon)
    return kFALSE;

  if (!fIsAOD)
  {
    TParticle *MCPrimaryPion = (TParticle *)fMCEvent->GetTrack(TMath::Abs(track_primiary_pion->GetLabel()))->Particle();
    TParticle *MCSecondaryPion = (TParticle *)fMCEvent->GetTrack(TMath::Abs(track_secondary_pion->GetLabel()))->Particle();
    TParticle *MCKaon = (TParticle *)fMCEvent->GetTrack(TMath::Abs(track_kaon->GetLabel()))->Particle();
    if (fSetModeK892orRho)
    {
      // K892 daughter (secondary pion, kaon) check
      if ((TMath::Abs(MCSecondaryPion->GetPdgCode()) != kPionCode || TMath::Abs(MCKaon->GetPdgCode()) != kKaonCode))
        return kFALSE;

      // K892 mother check
      if (MCSecondaryPion->GetMother(0) != MCKaon->GetMother(0))
        return kFALSE;

      // K892 check
      TParticle *MCK892 = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCSecondaryPion->GetMother(0)))->Particle();
      if (TMath::Abs(MCK892->GetPdgCode()) != kK892Code)
        return kFALSE;

      // Pimary pion check
      if (MCPrimaryPion->GetPdgCode() != kPionCode)
        return kFALSE;

      // Same mother check
      if (MCPrimaryPion->GetMother(0) != MCK892->GetMother(0))
        return kFALSE;

      // K1 check
      TParticle *MCK1 = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCPrimaryPion->GetMother(0)))->Particle();
      if (TMath::Abs(MCK1->GetPdgCode()) != kK1PCode)
        return kFALSE;

      return kTRUE;
    }
    else
    {
      // Rho daughter (primary pion, secondary pion) check
      if ((TMath::Abs(MCPrimaryPion->GetPdgCode()) != kPionCode || TMath::Abs(MCSecondaryPion->GetPdgCode()) != kPionCode))
        return kFALSE;

      // Rho mother check
      if (MCPrimaryPion->GetMother(0) != MCSecondaryPion->GetMother(0))
        return kFALSE;

      // Rho check
      TParticle *MCRho = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCPrimaryPion->GetMother(0)))->Particle();
      if (TMath::Abs(MCRho->GetPdgCode()) != kRhoCode)
        return kFALSE;

      // Kaon check
      if (MCKaon->GetPdgCode() != kKaonCode)
        return kFALSE;

      // Same mother check
      if (MCKaon->GetMother(0) != MCRho->GetMother(0))
        return kFALSE;

      // K1 check
      TParticle *MCK1 = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCRho->GetMother(0)))->Particle();
      if (TMath::Abs(MCK1->GetPdgCode()) != kK1PCode)
        return kFALSE;

      return kTRUE;
    }
  }
  else
  {
    AliAODMCParticle *MCPrimaryPion = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_primiary_pion->GetLabel()));
    AliAODMCParticle *MCSecondaryPion = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_secondary_pion->GetLabel()));
    AliAODMCParticle *MCKaon = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_kaon->GetLabel()));

    if (fSetModeK892orRho)
    {
      // K892 daughter (secondary pion, kaon) check
      if ((TMath::Abs(MCSecondaryPion->GetPdgCode()) != kPionCode || TMath::Abs(MCKaon->GetPdgCode()) != kKaonCode))
        return kFALSE;

      // K892 mother check
      if (MCSecondaryPion->GetMother() != MCKaon->GetMother())
        return kFALSE;

      // K892 check
      AliAODMCParticle *MCK892 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCSecondaryPion->GetMother()));
      if (TMath::Abs(MCK892->GetPdgCode()) != kK892Code)
        return kFALSE;

      // Pimary pion check
      if (MCPrimaryPion->GetPdgCode() != kPionCode)
        return kFALSE;

      // Same mother check
      if (MCPrimaryPion->GetMother() != MCK892->GetMother())
        return kFALSE;

      // K1 check
      AliAODMCParticle *MCK1 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCPrimaryPion->GetMother()));
      if (TMath::Abs(MCK1->GetPdgCode()) != kK1PCode)
        return kFALSE;

      return kTRUE;
    }
    else
    {
      // Rho daughter (primary pion, secondary pion) check
      if ((TMath::Abs(MCPrimaryPion->GetPdgCode()) != kPionCode || TMath::Abs(MCSecondaryPion->GetPdgCode()) != kPionCode))
        return kFALSE;

      // Rho mother check
      if (MCPrimaryPion->GetMother() != MCSecondaryPion->GetMother())
        return kFALSE;

      // Rho check
      AliAODMCParticle *MCRho = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCPrimaryPion->GetMother()));
      if (TMath::Abs(MCRho->GetPdgCode()) != kRhoCode)
        return kFALSE;

      // Kaon check
      if (MCKaon->GetPdgCode() != kKaonCode)
        return kFALSE;

      // Same mother check
      if (MCKaon->GetMother() != MCRho->GetMother())
        return kFALSE;

      // K1 check
      AliAODMCParticle *MCK1 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCRho->GetMother()));
      if (TMath::Abs(MCK1->GetPdgCode()) != kK1PCode)
        return kFALSE;

      return kTRUE;
    }
  }
}
THnSparseD *AliAnalysisTaskK1::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt)
{
  TArrayD xmin(ndim), xmax(ndim);
  TArrayI nbins(ndim);
  std::vector<TArrayD> binnings;

  for (int idim = 0; idim < ndim; ++idim)
  {
    const TAxis &axis = bins[idim];
    nbins[idim] = axis.GetNbins();
    xmin[idim] = axis.GetXmin();
    xmax[idim] = axis.GetXmax();
    TArrayD binEdges(nbins[idim] + 1);
    for (int i = 0; i <= nbins[idim]; ++i)
    {
      binEdges[i] = axis.GetBinLowEdge(i + 1);
    }
    binnings.push_back(binEdges);
  }
  THnSparseD *hsparse = new THnSparseD(name.Data(), title.Data(), ndim, nbins.GetArray(), xmin.GetArray(), xmax.GetArray());

  // Set the bin edges for each axis
  for (int id = 0; id < ndim; ++id)
  {
    hsparse->GetAxis(id)->Set(nbins[id], binnings[id].GetArray());
    TString axisName = bins[id].GetName();
    hsparse->GetAxis(id)->SetName(axisName.Data());
  }
  TString optionstring(opt);
  optionstring.ToLower();
  if (optionstring.Contains("s"))
  {
    hsparse->Sumw2();
  }

  return hsparse;
}
Long64_t AliAnalysisTaskK1::FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w)
{
  if (fSkipFillingHistogram)
    return 0;
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  if (int(x.size()) != h->GetNdimensions())
  {
    std::cout << "ERROR : wrong sized of array while Fill " << h->GetName() << std::endl;
    exit(1);
  }
  return h->Fill(&x.front(), w);
}
TAxis AliAnalysisTaskK1::AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(nbin, xmin, xmax);
  axis.SetName(name);
  return axis;
}
TAxis AliAnalysisTaskK1::AxisStr(TString name, std::vector<TString> bin)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
  UInt_t i = 1;
  for (auto blabel : bin)
    ax.SetBinLabel(i++, blabel);
  return ax;
}
TAxis AliAnalysisTaskK1::AxisVar(TString name, std::vector<Double_t> bin)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(bin.size() - 1, &bin.front());
  axis.SetName(name);
  return axis;
}
Double_t AliAnalysisTaskK1::GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type)
{
  AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
  if (nanoT)
  {
    static bool used = false;
    if (!used)
    {
      AliNanoAODTrack::InitPIDIndex();
      used = true;
    }
    return nanoT->GetVar(
        AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, type));
  }
  else
    return fPIDResponse->NumberOfSigmasTPC(track, type);
}
Double_t AliAnalysisTaskK1::GetTOFnSigma(AliVTrack *track, AliPID::EParticleType type)
{
  AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
  if (nanoT)
  {
    if (!nanoT->HasTOFpid())
      return -999.;
    static bool used = false;
    if (!used)
    {
      AliNanoAODTrack::InitPIDIndex();
      used = true;
    }
    return nanoT->GetVar(
        AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTOF, type));
  }
  else
    return fPIDResponse->NumberOfSigmasTOF(track, type);
}
void AliAnalysisTaskK1::FillTrackToEventPool()
{
  // Fill Selected tracks to event mixing pool
  if ((fCentBin < 0) || (fZbin < 0))
    return;
  AliVTrack *goodtrack;

  tracklist *etl;
  eventpool *ep;
  // Event mixing pool

  ep = &fEMpool[fCentBin][fZbin];
  ep->push_back(tracklist());
  etl = &(ep->back());
  // Fill selected tracks
  for (UInt_t i = 0; i < fGoodPrimaryPionArray.size(); i++)
  {
    goodtrack = (AliVTrack *)fEvt->GetTrack(fGoodPrimaryPionArray[i]);
    if (!goodtrack)
      continue;
    etl->push_back((AliVTrack *)goodtrack->Clone());
  }
  if (!fGoodPrimaryPionArray.size())
    ep->pop_back();
  if ((int)ep->size() > (int)fnMix)
  {
    for (auto it : ep->front())
      delete it;
    ep->pop_front();
  }
}
void AliAnalysisTaskK1::PostAllData()
{
  if (!fSkipFillingHistogram)
    PostData(1, fList);
  if (fFillTree)
  {
    PostData(2, fNanoTree);
    if (fIsMC)
      PostData(3, fNanoMCTree);
  }
}
void AliAnalysisTaskK1::GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3])
{
  AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
  if (nanoT)
    nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
  else
    track->GetImpactParameters(p, cov);
}
std::tuple<Float_t, Float_t, std::array<Float_t, 3>> AliAnalysisTaskK1::GetImpactParametersWrapper(AliVTrack *track)
{
  Float_t p[2] = {0, 0};
  Float_t cov[3] = {0, 0, 0};
  GetImpactParam(track, p, cov); // Assuming this is accessible within the context
  return {p[0], p[1], {cov[0], cov[1], cov[2]}};
}
Int_t AliAnalysisTaskK1::TrackSelection(AliVTrack *track, Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t nTPCNSigKaon, Float_t nTOFNSigKaon, Float_t lpT, Float_t lDCAz, Float_t lDCAr, Float_t lEta)
{
  // Check for valid track pointer
  if (!track)
    return 0;

  Int_t result = 0;
  Double_t lsigmaDCAr = (0.0026 + 0.0050 / lpT);

  // Evaluate each track type selection
  if (IsPassPrimaryPionSelection(nTPCNSigPion, nTOFNSigPion, lEta, lDCAz, lDCAr, lsigmaDCAr))
    result |= kPrimaryPionPass;

  if (IsPassSecondaryPionSelection(nTPCNSigPion, nTOFNSigPion, lEta, lDCAz, lDCAr, lsigmaDCAr))
    result |= kSecondaryPionPass;

  if (IsPassKaonSelection(nTPCNSigKaon, nTOFNSigKaon, lEta, lDCAz, lDCAr, lsigmaDCAr))
    result |= kKaonPass;

  return result;
}
bool AliAnalysisTaskK1::IsPassPrimaryPionSelection(Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr)
{
  return TMath::Abs(nTPCNSigPion) <= fTPCNsigPrimaryPionCut &&
         (fTOFNsigPrimaryPionCut <= 0 || TMath::Abs(nTOFNSigPion) <= fTOFNsigPrimaryPionCut) &&
         TMath::Abs(lEta) <= fPrimaryPionEtaCut &&
         lDCAz <= fPrimaryPionZVertexCut &&
         lDCAr <= lsigmaDCAr * fPrimaryPionXYVertexSigmaCut;
}

bool AliAnalysisTaskK1::IsPassSecondaryPionSelection(Float_t nTPCNSigPion, Float_t nTOFNSigPion, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr)
{
  return TMath::Abs(nTPCNSigPion) <= fTPCNsigSecondaryPionCut &&
         (fTOFNsigSecondaryPionCut <= 0 || TMath::Abs(nTOFNSigPion) <= fTOFNsigSecondaryPionCut) &&
         TMath::Abs(lEta) <= fSecondaryPionEtaCut &&
         lDCAz <= fSecondaryPionZVertexCut &&
         lDCAr <= lsigmaDCAr * fSecondaryPionXYVertexSigmaCut;
}

bool AliAnalysisTaskK1::IsPassKaonSelection(Float_t nTPCNSigKaon, Float_t nTOFNSigKaon, Float_t lEta, Float_t lDCAz, Float_t lDCAr, Double_t lsigmaDCAr)
{
  return TMath::Abs(nTPCNSigKaon) <= fTPCNsigKaonCut &&
         (fTOFNsigKaonCut <= 0 || TMath::Abs(nTOFNSigKaon) <= fTOFNsigKaonCut) &&
         TMath::Abs(lEta) <= fKaonEtaCut &&
         lDCAz <= fKaonZVertexCut &&
         lDCAr <= lsigmaDCAr * fKaonXYVertexSigmaCut;
}
bool AliAnalysisTaskK1::IsTrackAccepted(AliVTrack *track)
{
  if (!fIsAOD)
  {
    return fTrackCuts->AcceptTrack(static_cast<AliESDtrack *>(track));
  }
  else if (!fIsNano)
  {
    auto aodTrack = static_cast<AliAODTrack *>(track);
    return aodTrack->TestFilterBit(fFilterBit);
  }
  else
  {
    return static_cast<AliNanoAODTrack *>(track)->TestFilterBit(fFilterBit);
  }
}
