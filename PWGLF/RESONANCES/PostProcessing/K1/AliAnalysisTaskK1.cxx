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

enum
{
  kPrimaryPionPass = BIT(0),
  kSecondaryPionPass = BIT(1),
  kKaonPass = BIT(2)
};
enum
{
  kNormal = 1,
  kAnti = 2
};
enum
{                  // PDG Code
  kK1PCode = 3224, // K1(1270)+
  kK1NCode = 3114, // K1(1270)-
  kPionCode = 211, // Pion+
  kKaonCode = 321, // Kaon-
  kK892Code = 313, // K892
  kRhoCode = 113   // Rho0
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
      fHistos{nullptr},
      fVertex{nullptr},
      fRTree{nullptr},
      fSTree{nullptr},
      fRecK1{},
      fSimK1{},
      fMCArray{nullptr},
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
      fEMpool{},
      fBinCent{},
      fBinZ{},
      fPosPV{},
      fMagField{0},
      fCent{-1},
      fnMix{10},
      fCentBin{-1},
      fZbin{-1},
      fFilterBit{32},
      fTPCNsigPrimaryPionCut{3},
      fPrimaryPionEtaCut{0.8},
      fPrimaryPionZVertexCut{2.0},
      fPrimaryPionXYVertexSigmaCut{3},
      fTPCNsigSecondaryPionCut{3},
      fSecondaryPionEtaCut{0.8},
      fSecondaryPionZVertexCut{2.0},
      fSecondaryPionXYVertexSigmaCut{3},
      fTPCNsigKaonCut{3},
      fKaonEtaCut{0.8},
      fKaonZVertexCut{2.0},
      fKaonXYVertexSigmaCut{3},
      fK892MassWindowCut{0.1},
      fK892RapCut{0.5},
      fK1YCutHigh{0.5},
      fK1YCutLow{-0.5},
      fMinK1OA{0},
      fMaxK1OA{0.87},
      fMinPairAsym{-0.1},
      fMaxPairAsym{1},
      fMinK1PiPi{0},
      fMaxK1PiPi{1},
      fMinK1PiKa{0},
      fMaxK1PiKa{999}
{
  /// Default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskK1::AliAnalysisTaskK1(const char *name, Bool_t MCcase)
    : AliAnalysisTaskSE(name),
      fTrackCuts{nullptr}, 
      fPIDResponse{nullptr},
      fHistos{nullptr},
      fVertex{nullptr},
      fRTree{nullptr},
      fSTree{nullptr},
      fRecK1{},
      fSimK1{},
      fMCArray{nullptr},
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
      fEMpool{},
      fBinCent{},
      fBinZ{},
      fPosPV{},
      fMagField{0},
      fCent{-1},
      fnMix{10},
      fCentBin{-1},
      fZbin{-1},
      fFilterBit{32},
      fTPCNsigPrimaryPionCut{3},
      fPrimaryPionEtaCut{0.8},
      fPrimaryPionZVertexCut{2.0},
      fPrimaryPionXYVertexSigmaCut{7.0},
      fTPCNsigSecondaryPionCut{3},
      fSecondaryPionEtaCut{0.8},
      fSecondaryPionZVertexCut{2.0},
      fSecondaryPionXYVertexSigmaCut{7.0},
      fTPCNsigKaonCut{3},
      fKaonEtaCut{0.8},
      fKaonZVertexCut{2.0},
      fKaonXYVertexSigmaCut{7.0},
      fK892MassWindowCut{0.1},
      fK892RapCut{0.5},
      fK1YCutHigh{0.5},
      fK1YCutLow{-0.5},
      fMinK1OA{0},
      fMaxK1OA{0.87},
      fMinPairAsym{-0.1},
      fMaxPairAsym{1},
      fMinK1PiPi{0},
      fMaxK1PiPi{1},
      fMinK1PiKa{0},
      fMaxK1PiKa{999}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskK1::~AliAnalysisTaskK1()
{
  if (fRTree)
    delete fRTree;
  if (fSTree)
    delete fSTree;
}
//___________________________________________________________________
void AliAnalysisTaskK1::SetCutOpen()
{
  // Set cut open
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::UserCreateOutputObjects()
{
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fHistos = new THistManager("K1hists");
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
          : centaxisbin = {
                -1, 0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100}; // can be use from pp to PbPb
  if (fIsMC)
    centaxisbin = {
        -1, 0, 0.001, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100}; // for general MC
  fBinCent = AxisVar("Cent", centaxisbin);

  auto binPt = AxisFix("Pt", 150, 0, 15);
  auto binMass = AxisFix("Mass", 700, 0.6, 2.0);
  auto binMassMC = AxisFix("Mass", 450, 0.6, 1.5);
  // fBinZ = AxisFix("Z", 20, -10, 10); // 1cm diff
  fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10}); // moderate diff

  CreateTHnSparse("K1_data", "K1_data", 5, {binAnti, binType, fBinCent, binPt, binMass}, "s");
  if (fIsMC)
  {
    auto binTypeMCNorm = AxisStr("Type", {"kAll", "kINEL10", "kINEL_trig", "kINEL_trig_vtx",
                                          "kINEL_trig_vtx10", "kINELg0", "kINELg010", "kINELg0_trig",
                                          "kINELg0_trig_vtx", "kINELg0_trig_vtx10", "kSelected"});
    CreateTHnSparse("K1_mc", "K1_mc", 5, {binAnti, binTypeMC, fBinCent, binPt, binMassMC}, "s");
    CreateTHnSparse("Normalisation", "", 2, {binTypeMCNorm, fBinCent}, "s");
  }
  if (fFillQAPlot)
  {
    fHistos->CreateTH2("QA/hTPCPIDTrack", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDTrackNsigVspT", "", 100, 0, 10, 200, -10, 10);
    fHistos->CreateTH1("QA/hEtaTrack", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QA/hDCAPVTrack", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCArPVTrack", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hPtTrack", "", 150, 0, 15);

    fHistos->CreateTH2("QAcut/hTPCPIDPrimaryPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDPrimaryPionNsigVspT", "", 100, 0, 10, 200, -10, 10);
    fHistos->CreateTH1("QAcut/hEtaPrimaryPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVPrimaryPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVPrimaryPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtPrimaryPion", "", 150, 0, 15);

    fHistos->CreateTH2("QAcut/hTPCPIDSecondaryPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDSecondaryPionNsigVspT", "", 100, 0, 10, 200, -10, 10);
    fHistos->CreateTH1("QAcut/hEtaSecondaryPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVSecondaryPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVSecondaryPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtSecondaryPion", "", 150, 0, 15);

    fHistos->CreateTH2("QAcut/hTPCPIDKaon", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDKaonNsigVspT", "", 100, 0, 10, 200, -10, 10);
    fHistos->CreateTH1("QAcut/hEtaKaon", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVKaon", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVKaon", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtKaon", "", 150, 0, 15);

    fHistos->CreateTH1("QAReso/hK1OA", "", 100, 0, 3.14);
    fHistos->CreateTH1("QAReso/K1PairAsymm", "", 100, -1, 1);
    fHistos->CreateTH2("QAReso/InvMass_piK_pipi", "", 60, 0.6, 1.2, 200, 0, 2.0);
    fHistos->CreateTH2("QAReso/InvMass_piK_pika", "", 60, 0.6, 1.2, 200, 0, 2.0);

    fHistos->CreateTH1("QAResoCut/hK1OA", "", 100, 0, 3.14);
    fHistos->CreateTH1("QAResoCut/K1PairAsymm", "", 100, -1, 1);
    fHistos->CreateTH2("QAResoCut/InvMass_piK_pipi", "", 60, 0.6, 1.2, 200, 0, 2.0);
    fHistos->CreateTH2("QAResoCut/InvMass_piK_pika", "", 60, 0.6, 1.2, 200, 0, 2.0);
  }
  fEventCut.AddQAplotsToList(fHistos->GetListOfHistograms()); // QA histograms from AliEventCuts
  if (fIsHM) fHistos->CreateTH1("hMultiplicity", "", 100, 0, 0.1, "s");
  else fHistos->CreateTH1("hMultiplicity", "", 101, -1, 100, "s");
  fEMpool.resize(fBinCent.GetNbins() + 1, std::vector<eventpool>(fBinZ.GetNbins() + 1));

  PostData(1, fHistos->GetListOfHistograms());
  if (fFillTree)
  {
    OpenFile(1);
    fRTree = new TTree("RTree", "Reconstructed resonance tree");
    fRTree->Branch("RK1Resonance", &fRecK1);
    PostData(2, fRTree);

    if (fIsMC)
    {
      fSTree = new TTree("STree", "Simulated resonance tree");
      fSTree->Branch("SK1Resonance", &fSimK1);
      PostData(3, fSTree);
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (!event) // if there is no event, return
  {
    PostData(1, fHistos->GetListOfHistograms());
    if (fFillTree)
    {
      PostData(2, fRTree);
      if (fIsMC)
      {
        PostData(3, fSTree);
      }
    }
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
    PostData(1, fHistos->GetListOfHistograms());
    if (fFillTree)
    {
      PostData(2, fRTree);
      if (fIsMC)
      {
        PostData(3, fSTree);
      }
    }
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
        fMCArray =
            (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
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

  if (fIsMC)
  {
    // Fill Normalisation histogram
    if (IsVtxInZCut)
      FillMCinput(fMCEvent, 1);
    if (IsINEL0True && IsVtxInZCut)
      FillMCinput(fMCEvent, 2);
    if (IsSelectedTrig)
      FillMCinput(fMCEvent, 3);
    FillTHnSparse("Normalisation", {(int)kAll, (double)fCent});
    if (IsINEL0True)
    {
      FillTHnSparse("Normalisation", {(int)kINELg0, (double)fCent});
      if (IsVtxInZCut)
      {
        FillTHnSparse("Normalisation", {(int)kINELg010, (double)fCent});
      }
      if (IsSelectedTrig)
      {
        FillTHnSparse("Normalisation", {(int)kINELg0_trig, (double)fCent});
        if (IsGoodVertex)
        {
          FillTHnSparse("Normalisation", {(int)kINELg0_trig_vtx, (double)fCent});
          if (IsVtxInZCut)
          {
            FillTHnSparse("Normalisation", {(int)kINELg0_trig_vtx10, (double)fCent});
          }
        }
      }
    }
    if (IsVtxInZCut)
    {
      FillTHnSparse("Normalisation", {(int)kINEL10, (double)fCent});
    }
    if (IsSelectedTrig)
    {
      FillTHnSparse("Normalisation", {(int)kINEL_trig, (double)fCent});
      if (IsGoodVertex)
      {
        FillTHnSparse("Normalisation", {(int)kINEL_trig_vtx, (double)fCent});
        if (IsVtxInZCut)
        {
          FillTHnSparse("Normalisation", {(int)kINEL_trig_vtx10, (double)fCent});
        }
      }
    }
  }

  if (!IsEvtSelected)
  {
    PostData(1, fHistos->GetListOfHistograms());
    if (fFillTree)
    {
      PostData(2, fRTree);
      if (fIsMC)
      {
        PostData(3, fSTree);
      }
    }
    return;
  }

  fHistos->FillTH1("hMultiplicity", (double)fCent); // multiplicity distribution for basic event QA

  if (fIsMC)
  {
    FillMCinput(fMCEvent);
    FillTHnSparse("Normalisation", {(int)kSelected, (double)fCent});
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
    if (!fSkipFillingHistogram) FillHistograms(); // Fill the histogram
    if (fFillTree) FillTree();
  }
  if (fSetMixing && fGoodPrimaryPionArray.size())
  {
    FillTrackToEventPool(); // use only pion track pool.
  }

  PostData(1, fHistos->GetListOfHistograms());
  if (fFillTree)
  {
    PostData(2, fRTree);
    if (fIsMC)
    {
      PostData(3, fSTree);
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskK1::Terminate(Option_t *) {}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskK1::FillTrackPools()
{
  // Fill the track pool: primary pion, secondary pion, kaon at the same time
  const UInt_t nTracks = fEvt->GetNumberOfTracks();
  fGoodPrimaryPionArray.clear();
  fGoodSecondaryPionArray.clear();
  fGoodKaonArray.clear();
  AliVTrack *track = nullptr;
  Float_t b[2];
  Float_t bCov[3];
  Double_t nTPCNSigPion{0}, nTPCNSigKaon{0}, lDCAz{0}, lpT{0}, lsigmaDCAr{0}, lDCAr{0}, lEta{0};
  Bool_t isPassPrimaryPionSelection{false}, isPassSecondaryPionSelection{false}, isPassKaonSelection{false};
  Bool_t isAcceptedTrack{false};

  for (UInt_t it = 0; it < nTracks; it++)
  {
    isPassPrimaryPionSelection = true;
    isPassSecondaryPionSelection = true;
    isPassKaonSelection = true;
    track = (AliVTrack *)fEvt->GetTrack(it);
    if (!track) continue;

    // ---------- Track selection begin ----------
    isAcceptedTrack = false;
    if (!fIsAOD) {
      // ESD Case
      isAcceptedTrack = fTrackCuts->AcceptTrack(static_cast<AliESDtrack*>(track));
    } else if (!fIsNano) {
      // AOD Case, not Nano
      isAcceptedTrack = static_cast<AliAODTrack*>(track)->TestFilterBit(fFilterBit);
    } else {
      // Nano AOD Case
      isAcceptedTrack = static_cast<AliNanoAODTrack*>(track)->TestFilterBit(fFilterBit);
    }
    if (!isAcceptedTrack) continue;
    GetImpactParam(track, b, bCov);

    lDCAz = b[1];
    nTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
    nTPCNSigKaon = GetTPCnSigma(track, AliPID::kKaon);

    // boost up the speed
    // if the particle is not pion or kaon in 5 sigma, skip the track
    if (TMath::Abs(nTPCNSigPion) > 5 && TMath::Abs(nTPCNSigKaon) > 5) continue;

    lpT = track->Pt();
    lsigmaDCAr = (0.0026 + 0.0050 / lpT);
    lDCAr = b[0];
    lEta = track->Eta();

    if (fFillQAPlot) 
    {
      // Fill the QA plot for the whole track. we don't know which track is it yet.
      fHistos->FillTH1("QA/hDCAPVTrack", lDCAz);
      fHistos->FillTH1("QA/hDCArPVTrack", lDCAr);
      fHistos->FillTH1("QA/hEtaTrack", lEta);
      fHistos->FillTH1("QA/hPtTrack", lpT);
      fHistos->FillTH2("QA/hTPCPIDTrack", track->GetTPCmomentum(), track->GetTPCsignal());
      fHistos->FillTH1("QA/hTPCPIDTrackNsigVspT", lpT, nTPCNSigPion);
    }

    if (lpT < 0.15) // minimum pT cut
      continue;

    // Selection for primary pion
    if (TMath::Abs(nTPCNSigPion) > fTPCNsigPrimaryPionCut)
      isPassPrimaryPionSelection = false;
    if (TMath::Abs(lEta) > fPrimaryPionEtaCut)
      isPassPrimaryPionSelection = false;
    if (lDCAz > fPrimaryPionZVertexCut)
      isPassPrimaryPionSelection = false;
    if (lDCAr > lsigmaDCAr * fPrimaryPionXYVertexSigmaCut)
      isPassPrimaryPionSelection = false;

    // Selection for secondary pion
    if (TMath::Abs(nTPCNSigPion) > fTPCNsigSecondaryPionCut)
      isPassSecondaryPionSelection = false;
    if (TMath::Abs(lEta) > fSecondaryPionEtaCut)
      isPassSecondaryPionSelection = false;
    if (lDCAz > fSecondaryPionZVertexCut)
      isPassSecondaryPionSelection = false;
    if (lDCAr > lsigmaDCAr * fSecondaryPionXYVertexSigmaCut)
      isPassSecondaryPionSelection = false;
    
    // Selection for kaon
    if (TMath::Abs(nTPCNSigKaon) > fTPCNsigKaonCut)
      isPassKaonSelection = false;
    if (TMath::Abs(lEta) > fKaonEtaCut)
      isPassKaonSelection = false;
    if (lDCAz > fKaonZVertexCut)
      isPassKaonSelection = false;
    if (lDCAr > lsigmaDCAr * fKaonXYVertexSigmaCut)
      isPassKaonSelection = false;

    if (fFillQAPlot)
    {
      if (isPassPrimaryPionSelection)
      {
        fHistos->FillTH1("QAcut/hDCAPVPrimaryPion", lDCAz);
        fHistos->FillTH1("QAcut/hDCArPVPrimaryPion", lDCAr);
        fHistos->FillTH1("QAcut/hEtaPrimaryPion", lEta);
        fHistos->FillTH1("QAcut/hPtPrimaryPion", lpT);
        fHistos->FillTH2("QAcut/hTPCPIDPrimaryPion", track->GetTPCmomentum(), track->GetTPCsignal());
        fHistos->FillTH1("QAcut/hTPCPIDPrimaryPionNsigVspT", lpT, nTPCNSigPion);
      }
      if (isPassSecondaryPionSelection)
      {
        fHistos->FillTH1("QAcut/hDCAPVSecondaryPion", lDCAz);
        fHistos->FillTH1("QAcut/hDCArPVSecondaryPion", lDCAr);
        fHistos->FillTH1("QAcut/hEtaSecondaryPion", lEta);
        fHistos->FillTH1("QAcut/hPtSecondaryPion", lpT);
        fHistos->FillTH2("QAcut/hTPCPIDSecondaryPion", track->GetTPCmomentum(), track->GetTPCsignal());
        fHistos->FillTH1("QAcut/hTPCPIDSecondaryPionNsigVspT", lpT, nTPCNSigPion);
      }
      if (isPassKaonSelection)
      {
        fHistos->FillTH1("QAcut/hDCAPVKaon", lDCAz);
        fHistos->FillTH1("QAcut/hDCArPVKaon", lDCAr);
        fHistos->FillTH1("QAcut/hEtaKaon", lEta);
        fHistos->FillTH1("QAcut/hPtKaon", lpT);
        fHistos->FillTH2("QAcut/hTPCPIDKaon", track->GetTPCmomentum(), track->GetTPCsignal());
        fHistos->FillTH1("QAcut/hTPCPIDKaonNsigVspT", lpT, nTPCNSigKaon);
      }
    }

    if (isPassPrimaryPionSelection)
      fGoodPrimaryPionArray.push_back(it);
    if (isPassSecondaryPionSelection)
      fGoodSecondaryPionArray.push_back(it);
    if (isPassKaonSelection)
      fGoodKaonArray.push_back(it);
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
  Double_t lK892Mass{0}, lK1Rapidity{0}, lK1Angle{0}, lK1PairAsym{0}, lMassPiPi{0}, lMassPiKa{0};

  TLorentzVector vecPrimaryPion, vecSecondaryPion, vecKaon, vecPrimaryPionMix;
  TLorentzVector vecK892, vecK1, tempPiPi, tempPiKa;
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
      
      // Apply K(892) cut here.

      // Y cut
      if (abs(vecK892.Rapidity()) > fK892RapCut)
        continue;
      // mass window
      lK892Mass = vecK892.M();
      if (abs(lK892Mass - kK892Mass) > fK892MassWindowCut)
        continue;

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
        // Y cut
        lK1Rapidity = vecK1.Rapidity();
        if ((lK1Rapidity > fK1YCutHigh) || (lK1Rapidity < fK1YCutLow))
          continue;
        // OA Cut
        lK1Angle = vecK892.Angle(vecPrimaryPion.Vect());
        // Pair Assym
        lK1PairAsym = (vecK892.E() - vecPrimaryPion.E()) / (vecK892.E() + vecPrimaryPion.E());
        // PiPi, PiKa mass range cut
        tempPiPi = vecPrimaryPion + vecSecondaryPion;
        tempPiKa = vecPrimaryPion + vecKaon;
        lMassPiPi = tempPiPi.M();
        lMassPiKa = tempPiKa.M();
        if (fFillQAPlot)
        {
          fHistos->FillTH1("QAReso/hK1OA", lK1Angle);
          fHistos->FillTH1("QAReso/K1PairAsymm", lK1PairAsym);
          fHistos->FillTH2("QAReso/InvMass_piK_pipi", lK892Mass, lMassPiPi);
          fHistos->FillTH2("QAReso/InvMass_piK_pika", lK892Mass, lMassPiKa);
        }
        if ((lMassPiPi > fMaxK1PiPi)||(lMassPiPi < fMinK1PiPi)) 
          continue;
        if ((lMassPiKa > fMaxK1PiKa)||(lMassPiKa < fMinK1PiKa)) 
          continue;
        if ((lK1Angle > fMaxK1OA)||(lK1Angle < fMinK1OA)) 
          continue;
        if ((lK1PairAsym > fMaxPairAsym)||(lK1PairAsym < fMinPairAsym)) 
          continue;
        if (fFillQAPlot)
        {
          fHistos->FillTH1("QAResoCut/hK1OA", lK1Angle);
          fHistos->FillTH1("QAResoCut/K1PairAsymm", lK1PairAsym);
          fHistos->FillTH2("QAResoCut/InvMass_piK_pipi", lK892Mass, lMassPiPi);
          fHistos->FillTH2("QAResoCut/InvMass_piK_pika", lK892Mass, lMassPiKa);
        }

        if (track_primiary_pion->Charge() > 0)
          isPirmaryPionPlus = true;
        else
          isPirmaryPionPlus = false;

        (isAnti) ? binAnti = kAnti : binAnti = kNormal;
        (isPirmaryPionPlus) ? sign = kK1P : sign = kK1N;

        FillTHnSparse("K1_data", {(double)binAnti, (double)sign, (double)fCent, vecK1.Pt(), vecK1.M()});

        if (fIsMC)
        {
          if (IsTrueK1(primiaryID, secondaryID, kaonID))
          {
            (isPirmaryPionPlus) ? sign = kK1P_REC : sign = kK1N_REC;
            FillTHnSparse("K1_mc", {(double)binAnti, (double)sign, (double)fCent, vecK1.Pt(), vecK1.M()});
          }
          else
          {
            // MC not true bkg
            (isPirmaryPionPlus) ? sign = kK1P_NOT : sign = kK1N_NOT;
            FillTHnSparse("K1_mc", {(double)binAnti, (double)sign, (double)fCent, vecK1.Pt(), vecK1.M()});
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
          if ((vecK1.Rapidity() > fK1YCutHigh) || (vecK1.Rapidity() < fK1YCutLow))
            continue;
          // OA Cut
          lK1Angle = vecK892.Angle(vecPrimaryPionMix.Vect());
          // Pair Assym
          lK1PairAsym = (vecK892.E() - vecPrimaryPionMix.E()) / (vecK892.E() + vecPrimaryPionMix.E());
          // PiPi, PiKa mass range cut
          tempPiPi = vecPrimaryPionMix + vecSecondaryPion;
          tempPiKa = vecPrimaryPionMix + vecKaon;
          lMassPiPi = tempPiPi.M();
          lMassPiKa = tempPiKa.M();
          // PiPi, PiKa mass range cut
          if ((lMassPiPi > fMaxK1PiPi)||(lMassPiPi < fMinK1PiPi)) 
            continue;
          if ((lMassPiKa > fMaxK1PiKa)||(lMassPiKa < fMinK1PiKa)) 
            continue;
          if ((lK1Angle > fMaxK1OA)||(lK1Angle < fMinK1OA)) 
            continue;
          if ((lK1PairAsym > fMaxPairAsym)||(lK1PairAsym < fMinPairAsym)) 
            continue;
          
          if (track_primiary_pion_mix->Charge() > 0)
            isPirmaryPionPlus = true;
          else
            isPirmaryPionPlus = false;

          (isAnti) ? binAnti = kAnti : binAnti = kNormal;
          (isPirmaryPionPlus) ? sign = kK1P_MIX : sign = kK1N_MIX;

          FillTHnSparse("K1_data", {(double)binAnti, (double)sign, (double)fCent, vecK1.Pt(), vecK1.M()});
        } // mixing track loop
      } // mixing

    } // kaon loop
  } // secondary pion loop
}
void AliAnalysisTaskK1::FillTree()
{
  
}
void AliAnalysisTaskK1::FillMCinput(AliMCEvent *fMCEvent, int Fillbin)
{
  int sign = kAllType;
  int binAnti = 0;
  TLorentzVector vecPart1, vecPart2;
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
      Int_t k1PdgCode = mcInputTrack->GetPdgCode();
      if ((TMath::Abs(k1PdgCode) != kK1PCode) &&
          (TMath::Abs(k1PdgCode) != kK1NCode))
        continue;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;
      // Y cut
      if ((mcInputTrack->Y() > fK1YCutHigh) || (mcInputTrack->Y() < fK1YCutLow))
        continue;
      
      (k1PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
      if (TMath::Abs(k1PdgCode) == kK1PCode)
        sign = kK1P_GEN + (int)Fillbin * 2;
      if (TMath::Abs(k1PdgCode) == kK1NCode)
        sign = kK1N_GEN + (int)Fillbin * 2;

      FillTHnSparse("K1_mc", {(double)binAnti, (double)sign, (double)fCent, mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
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

      Int_t k1PdgCode = mcInputTrack->GetPdgCode();

      if ((TMath::Abs(k1PdgCode) != kK1PCode) &&
          (TMath::Abs(k1PdgCode) != kK1NCode))
        continue;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;

      // Y cut
      if ((mcInputTrack->Y() > fK1YCutHigh) || (mcInputTrack->Y() < fK1YCutLow))
        continue;
      
      (k1PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
      if (TMath::Abs(k1PdgCode) == kK1PCode)
        sign = kK1P_GEN + (int)Fillbin * 2;
      if (TMath::Abs(k1PdgCode) == kK1NCode)
        sign = kK1N_GEN + (int)Fillbin * 2;

      FillTHnSparse("K1_mc", {(double)binAnti, (double)sign, (double)fCent, mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
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

    // K892 daughter (secondary pion, kaon) check
    if ((TMath::Abs(MCSecondaryPion->GetPdgCode()) == kPionCode && TMath::Abs(MCKaon->GetPdgCode()) != kKaonCode))
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
    if ((TMath::Abs(MCK1->GetPdgCode()) != kK1PCode) &&
        (TMath::Abs(MCK1->GetPdgCode()) != kK1NCode))
      return kFALSE;

    return kTRUE;
  }
  else
  {
    AliAODMCParticle *MCPrimaryPion = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_primiary_pion->GetLabel()));
    AliAODMCParticle *MCSecondaryPion = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_secondary_pion->GetLabel()));
    AliAODMCParticle *MCKaon = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track_kaon->GetLabel()));

    // K892 daughter (secondary pion, kaon) check
    if ((TMath::Abs(MCSecondaryPion->GetPdgCode()) == kPionCode && TMath::Abs(MCKaon->GetPdgCode()) != kKaonCode))
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
    if ((TMath::Abs(MCK1->GetPdgCode()) != kK1PCode) &&
        (TMath::Abs(MCK1->GetPdgCode()) != kK1NCode))
      return kFALSE;

    return kTRUE;
  }
}

THnSparse *AliAnalysisTaskK1::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  const TAxis *axises[bins.size()];
  for (UInt_t i = 0; i < bins.size(); i++)
    axises[i] = &bins[i];
  THnSparse *h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
  return h;
}
Long64_t AliAnalysisTaskK1::FillTHnSparse(TString name, std::vector<Double_t> x, Double_t w)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  auto hsparse = dynamic_cast<THnSparse *>(fHistos->FindObject(name));
  if (!hsparse)
  {
    std::cout << "ERROR : no " << name << std::endl;
    exit(1);
  }
  return FillTHnSparse(hsparse, x, w);
}
Long64_t AliAnalysisTaskK1::FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  if (int(x.size()) != h->GetNdimensions())
  {
    std::cout << "ERROR : wrong sized of array while Fill " << h->GetName()
              << std::endl;
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
double AliAnalysisTaskK1::GetTPCnSigma(AliVTrack *track, AliPID::EParticleType type)
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
void AliAnalysisTaskK1::GetImpactParam(AliVTrack *track, Float_t p[2], Float_t cov[3])
{
  AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
  if (nanoT)
    nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
  else
    track->GetImpactParameters(p, cov);
}
Int_t AliAnalysisTaskK1::trackSelection(AliVTrack* track, Float_t& nTPCNSigPion, Float_t& nTPCNSigKaon, Float_t lpT, Float_t lDCAz, Float_t lDCAr, Float_t lEta) {
  // We have 3 cases for track selection: primary pion, secondary pion, kaon
  // the return value will be the BIT combination of the track type
  // BIT(1) : kPrimaryPionPass
  // BIT(2) : kSecondaryPionPass
  // BIT(3) : kKaonPass
  
  // Initialize the output value for all pass at the beggining
  Int_t result = 0;

  if (!track) return result;

  if (TMath::Abs(nTPCNSigPion) < fTPCNsigPrimaryPionCut)
    result |= kPrimaryPionPass;
  if (TMath::Abs(lEta) > fPrimaryPionEtaCut)
    result |= kPrimaryPionPass;
  if (lDCAz > fPrimaryPionZVertexCut)
    result |= kPrimaryPionPass;
  Double_t lsigmaDCAr = (0.0026 + 0.0050 / lpT);
  if (lDCAr > lsigmaDCAr * fPrimaryPionXYVertexSigmaCut)
    result |= kPrimaryPionPass;


}