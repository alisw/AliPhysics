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

/* AliAnalysisTaskSigma1385PM
 *
 *  Test code for the reconstructing Sigma(1385)^{+-}
 *  Output could be saved to nTuple by using SetFillnTuple(kTRUE)
 *    -> can be used for TMVA input
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
#include "AliPIDResponse.h"
#include "TChain.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODv0.h"
#include "AliAnalysisTaskSigma1385PM.h"
#include "AliESDv0.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t pi = TMath::Pi();
const Double_t pionMass = AliPID::ParticleMass(AliPID::kPion);
const Double_t v0Mass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

enum { kNormal = 1, kAnti = 2 };
enum {
  kSigmaStarPCode = 3224,  // Sigma(1385)+
  kSigmaStarNCode = 3114,  // Sigma(1385)-
  kXiCode = 3312,          // Xi-
  kXiZeroCode = 3322,      // Xi0
  kLambdaCode = 3122,      // Lambda
  kProtonCode = 2212,      // Proton+
  kPionCode = 211,         // Pion+
  kXiStarCode = 3324,      // Xi(1530)0
  kXiStarMCode = 3314,     // Xi(1530)-
  kLambdaStarCode = 3124,  // Lambda1520
  kSigmaZeroCode = 3212    // Sigma0
};
enum {
  kSigmaStarP = 1,
  kSigmaStarN,
  kSigmaStarP_MIX,
  kSigmaStarN_MIX,
  kSigmaStarP_NOT,
  kSigmaStarN_NOT,
  kAllType
};
enum {
  kSigmaStarP_GEN = 1,  // 1
  kSigmaStarN_GEN,
  kSigmaStarP_GEN_INEL10,
  kSigmaStarN_GEN_INEL10,
  kSigmaStarP_GEN_INEL10_IGZ,
  kSigmaStarN_GEN_INEL10_IGZ,
  kSigmaStarP_GEN_TRIG,
  kSigmaStarN_GEN_TRIG,
  kSigmaStarP_REC,
  kSigmaStarN_REC,
  kSigmaStarPBkg_type1,
  kSigmaStarNBkg_type1,
  kSigmaStarPBkg_type2,
  kSigmaStarNBkg_type2,
  kSigmaStarPBkg_type3,
  kSigmaStarNBkg_type3,
  kSigmaStarPBkg_type4,
  kSigmaStarNBkg_type4,
  kSigmaStarPBkg_type5,
  kSigmaStarNBkg_type5
};
enum {
  kAll = 1,  // 1
  kINEL10,
  kINEL_trig,
  kINEL_trig_vtx,
  kINEL_trig_vtx10,
  kINELg0,  // 6
  kINELg010,
  kINELg0_trig,
  kINELg0_trig_vtx,
  kINELg0_trig_vtx10,
  kSelected  // 11
};

class AliAnalysisTaskSigma1385PM;

ClassImp(AliAnalysisTaskSigma1385PM)
    AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM()
    : AliAnalysisTaskSE(),
      fTrackCuts(nullptr),
      fEventCuts(),
      fPIDResponse(),
      fEvt(nullptr),
      fMCEvent(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fNtupleSigma1385(nullptr),
      fMCArray(nullptr),
      fIsAOD(kFALSE),
      fIsNano(kFALSE),
      fSetMixing(kFALSE),
      fFillQAPlot(kTRUE),
      fIsMC(kFALSE),
      fIsPrimaryMC(kFALSE),
      fFillnTuple(kFALSE),
      fIsINEL(kFALSE),
      fIsHM(kFALSE),
      fUseAsymmCut(kFALSE),
      fEMpool(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fFilterBit(32.0),
      fTPCNsigSigmaStarPionCut(3.0),
      fSigmaStarPionEtaCut(0.8),
      fSigmaStarPionZVertexCut(2.0),
      fSigmaStarPionXYVertexSigmaCut(7.0),
      fTPCNsigLambdaProtonCut(3.0),
      fTPCNsigLambdaPionCut(3.0),
      fDCADistLambdaDaughtersCut(1.6),
      fDCArDistLambdaPVCut(0.3),
      fDCAPositiveTrack(0.05),
      fDCANegativeTrack(0.05),
      fV0CosineOfPointingAngleCut(0.99),
      fMaxLambdaRapidity(0.8),
      fLambdaLowRadius(1.4),
      fLambdaHighRadius(100.0),
      fLambdaLifetime(30.0),
      fV0MassWindowCut(0.01),
      fSigmaStarYCutHigh(0.5),
      fSigmaStarYCutLow(-0.5),
      fSigmaStarAsymmCutHigh(0.0),
      fSigmaStarAsymmCutLow(1.0),
      fGoodTrackArray(),
      fGoodV0Array() {
  /// Default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskSigma1385PM::AliAnalysisTaskSigma1385PM(const char* name,
                                                       Bool_t MCcase)
    : AliAnalysisTaskSE(name),
      fTrackCuts(nullptr),
      fEventCuts(),
      fPIDResponse(),
      fEvt(nullptr),
      fMCEvent(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fNtupleSigma1385(nullptr),
      fMCArray(nullptr),
      fIsAOD(kFALSE),
      fIsNano(kFALSE),
      fSetMixing(kFALSE),
      fFillQAPlot(kTRUE),
      fIsMC(MCcase),
      fIsPrimaryMC(kFALSE),
      fFillnTuple(kFALSE),
      fIsINEL(kFALSE),
      fIsHM(kFALSE),
      fUseAsymmCut(kFALSE),
      fEMpool(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fFilterBit(32),
      fTPCNsigSigmaStarPionCut(3.0),
      fSigmaStarPionEtaCut(0.8),
      fSigmaStarPionZVertexCut(2.0),
      fSigmaStarPionXYVertexSigmaCut(7.0),
      fTPCNsigLambdaProtonCut(3.0),
      fTPCNsigLambdaPionCut(3.0),
      fDCADistLambdaDaughtersCut(1.6),
      fDCArDistLambdaPVCut(0.3),
      fDCAPositiveTrack(0.05),
      fDCANegativeTrack(0.05),
      fV0CosineOfPointingAngleCut(0.99),
      fMaxLambdaRapidity(0.8),
      fLambdaLowRadius(1.4),
      fLambdaHighRadius(100.0),
      fLambdaLifetime(30.0),
      fV0MassWindowCut(0.01),
      fSigmaStarYCutHigh(0.5),
      fSigmaStarYCutLow(-0.5),
      fSigmaStarAsymmCutHigh(0.0),
      fSigmaStarAsymmCutLow(1.0),
      fGoodTrackArray(),
      fGoodV0Array() {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TNtupleD::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskSigma1385PM::~AliAnalysisTaskSigma1385PM() {}
//___________________________________________________________________
void AliAnalysisTaskSigma1385PM::SetCutOpen() {
  // Pion cuts
  SetFilterbitSigmaStarPion(1);
  SetMaxNsigSigmaStarPion(5);
  SetMaxEtaSigmaStarPion(0.8);
  SetMaxVertexZSigmaStarPion(99);
  SetMaxVertexXYsigSigmaStarPion(99);

  // Lambda cuts
  SetMaxNsigV0Proton(5);
  SetMaxNsigV0Pion(5);
  SetMaxDCAPVV0PosDaughter(99);
  SetMaxDCAPVV0NegDaughter(99);
  SetMaxDCAV0daughters(999);
  SetMaxDCAPVV0(999);
  SetMinCPAV0(0.9);
  SetMaxRapidityV0(0.8);
  SetMaxMassWindowV0(999);
  SetLowRadiusV0(0);
  SetHighRadiusV0(999);
  SetLifetimeV0(999);

  // Sigma Star cut
  SetSigmaStarRapidityCutHigh(1);
  SetSigmaStarRapidityCutLow(-1);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::UserCreateOutputObjects() {
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fHistos = new THistManager("Sigma1385hists");
  auto binAnti = AxisStr("AType", {"Normal", "Anti"});
  auto binType = (!fIsMC)
                     ? AxisStr("Type", {"SigmaStarP", "SigmaStarN",
                                        "SigmaStarP_mix", "SigmaStarN_mix"})
                     : AxisStr("Type", {"SigmaStarP", "SigmaStarN",
                                        "SigmaStarP_mix", "SigmaStarN_mix",
                                        "SigmaStarP_no", "SigmaStarN_no"});
  auto binTypeMC = AxisStr("Type", {"SigmaStarP_gen",
                                    "SigmaStarN_gen",
                                    "SigmaStarP_gen_inel10",
                                    "SigmaStarN_gen_inel10",
                                    "SigmaStarP_gen_inel10_igz",
                                    "SigmaStarN_gen_inel10_igz",
                                    "SigmaStarP_gen_trig",
                                    "SigmaStarN_gen_trig",
                                    "SigmaStarP_rec",
                                    "SigmaStarN_rec",
                                    "kSigmaStarPBkg_type1",
                                    "kSigmaStarNBkg_type1",
                                    "kSigmaStarPBkg_type2",
                                    "kSigmaStarNBkg_type2",
                                    "kSigmaStarPBkg_type3",
                                    "kSigmaStarNBkg_type3",
                                    "kSigmaStarPBkg_type4",
                                    "kSigmaStarNBkg_type4",
                                    "kSigmaStarPBkg_type5",
                                    "kSigmaStarNBkg_type5"});

  std::vector<double> centaxisbin;
  (fIsHM)
      ? centaxisbin = {0, 0.001, 0.01, 0.05, 0.1}
      : centaxisbin = {0,  1,  5,  10, 15, 20, 30, 40,
                       50, 60, 70, 80, 90, 100};  // can be use from pp to PbPb
  if (fIsMC)
    centaxisbin = {0,  0.001, 0.01, 0.05, 0.1, 1,  5,  10, 15, 20,
                   30, 40,    50,   60,   70,  80, 90, 100};  // for general MC
  fBinCent = AxisVar("Cent", centaxisbin);
  auto binPt = AxisFix("Pt", 200, 0, 20);
  auto binMass = AxisFix("Mass", 1800, 1.2, 3.0);
  auto binMassMC = AxisFix("Mass", 800, 1.2, 2.0);
  // fBinZ = AxisFix("Z", 20, -10, 10); // 1cm diff
  fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});  // moderate diff

  CreateTHnSparse("Sigma1385_data", "Sigma1385_data", 5,
                  {binAnti, binType, fBinCent, binPt, binMass}, "s");
  if (fIsMC) {
    auto binTypeMCNorm = AxisStr(
        "Type", {"kAll", "kINEL10", "kINEL_trig", "kINEL_trig_vtx",
                 "kINEL_trig_vtx10", "kINELg0", "kINELg010", "kINELg0_trig",
                 "kINELg0_trig_vtx", "kINELg0_trig_vtx10", "kSelected"});
    CreateTHnSparse("Sigma1385_mc", "Sigma1385_mc", 5,
                    {binAnti, binTypeMC, fBinCent, binPt, binMassMC}, "s");
    CreateTHnSparse("Normalisation", "", 2, {binTypeMCNorm, fBinCent}, "s");
  }
  fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
  if (fIsHM)
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 0.1, "s");
  else
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 100, "s");
  if (fFillQAPlot) {
    fHistos->CreateTH2("QA/hTPCPIDPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QA/hEtaPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QA/hDCAPVPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCArPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hPtPion", "", 200, 0, 20);
    fHistos->CreateTH2("QA/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 200, 0,
                       200);
    fHistos->CreateTH2("QA/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QA/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCAlambdaPV", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hDCALambdaPVProton", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hDCALambdaPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hCosPAlambda", "", 50, 0.95, 1.0, "s");
    fHistos->CreateTH1("QA/hLifetimelambda", "", 50, 0, 50, "s");
    fHistos->CreateTH1("QA/hYlambda", "", 140, -0.7, 0.7, "s");
    fHistos->CreateTH1("QA/hMassLambda", "", 80, 1.08, 1.16, "s");
    fHistos->CreateTH1("QA/hLambdaRxy", "", 200, 0, 200);
    fHistos->CreateTH2("QA/hSigmaStarAsymm", "", 100, 0, 1, 100, 0, 100);
    if (fIsMC){
      fHistos->CreateTH2("QA/hSigmaStarAsymm_true", "", 100, 0, 1, 100, 0, 100);
      fHistos->CreateTH2("QA/hSigmaStarAsymm_true_selected", "", 100, 0, 1, 100, 0, 100);
    }
    fHistos->CreateTH1("QA/hLambdaAntiCheck", "", 4, -0.5, 3.5);

    fHistos->CreateTH2("QAcut/hTPCPIDPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QAcut/hEtaPion", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVPion", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtPion", "", 200, 0, 20);
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0,
                       200);
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 200, 0,
                       200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 200, 0,
                       200);
    fHistos->CreateTH1("QAcut/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCAlambdaPV", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hDCALambdaPVProton", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hDCALambdaPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hCosPAlambda", "", 50, 0.95, 1.0, "s");
    fHistos->CreateTH1("QAcut/hLifetimelambda", "", 50, 0, 50, "s");
    fHistos->CreateTH1("QAcut/hYlambda", "", 140, -0.7, 0.7, "s");
    fHistos->CreateTH1("QAcut/hMassLambda", "", 80, 1.08, 1.16, "s");
    fHistos->CreateTH1("QAcut/hLambdaRxy", "", 200, 0, 200);
    fHistos->CreateTH2("QAcut/hSigmaStarAsymm", "", 100, 0, 1, 100, 0, 100);
    if (fIsMC){
      fHistos->CreateTH2("QAcut/hSigmaStarAsymm_true", "", 100, 0, 1, 100, 0, 100);
      fHistos->CreateTH2("QAcut/hSigmaStarAsymm_true_selected", "", 100, 0, 1, 100, 0, 100);
    }
    fHistos->CreateTH1("QAcut/hLambdaAntiCheck", "", 4, -0.5, 3.5);
  }
  fEMpool.resize(fBinCent.GetNbins() + 1,
                 std::vector<eventpool>(fBinZ.GetNbins() + 1));

  fNtupleSigma1385 = new TNtupleD(
      "fNtupleSigma1385", "Sigma1385",
      "PIDSigmaStarPion:DCASigmaStarPionToPrimVertexZ:"
      "DCASigmaStarPionToPrimVertexR:EtaSigmaStarPion:PhiSigmaStarPion:"
      "PIDV0pTrackProton:PIDV0pTrackPion:PIDV0nTrackProton:PIDV0pTrackPion:"
      "DCAV0Daughters:DCAV0ToPrimVertex:CosPointingAngleV0:V0Mass:EtaV0:"
      "PhiV0:MCflag:Antiflag");

  PostData(1, fHistos->GetListOfHistograms());
  PostData(2, fNtupleSigma1385);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::UserExec(Option_t*) {
  AliVEvent* event = InputEvent();
  if (!event) {
    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fNtupleSigma1385);
    AliInfo("Could not retrieve event");
    return;
  }
  AliNanoAODHeader* nanoHeader =
      dynamic_cast<AliNanoAODHeader*>(fInputEvent->GetHeader());

  event->IsA() == AliESDEvent::Class()
      ? fEvt = dynamic_cast<AliESDEvent*>(event)
      : fEvt = dynamic_cast<AliAODEvent*>(event);
  if (!fIsAOD && (event->IsA() != AliESDEvent::Class()))
    fIsAOD = true;
  if (!fEvt) {
    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fNtupleSigma1385);
    return;
  }
  AliInputEventHandler* inputHandler =
      (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
          ->GetInputEventHandler();

  bool IsEvtSelected{false}, IsINEL0True{false}, IsSelectedTrig{false},
      IsVtxInZCut{false}, IsGoodVertex{false};
  if (!nanoHeader) {
    IsEvtSelected = fEventCuts.AcceptEvent(event);
    if (fIsMC) {
      if (fIsAOD)
        fMCArray =
            (TClonesArray*)fEvt->FindListObject("mcparticles");  // AOD Case
      fMCEvent = MCEvent();
      IsINEL0True = fEventCuts.IsTrueINELgtZero(fEvt, true);
    }
    fCent = fEventCuts.GetCentrality(0);
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if (!fPIDResponse)
      AliInfo("No PIDd");
    IsVtxInZCut = fEventCuts.PassedCut(AliEventCuts::kVertexPosition);
    IsSelectedTrig = fEventCuts.PassedCut(AliEventCuts::kTrigger);
    IsGoodVertex = fEventCuts.PassedCut(AliEventCuts::kVertexQuality);
  } else {
    if (!fIsNano)
      fIsNano = kTRUE;
    IsEvtSelected = true;
    fCent = nanoHeader->GetCentr("V0M");
  }

  if (fIsMC) {
    if (IsVtxInZCut)
      FillMCinput(fMCEvent, 1);
    if (IsINEL0True && IsVtxInZCut)
      FillMCinput(fMCEvent, 2);
    if (IsSelectedTrig)
      FillMCinput(fMCEvent, 3);
    FillTHnSparse("Normalisation", {(int)kAll, (double)fCent});
    if (IsINEL0True) {
      FillTHnSparse("Normalisation", {(int)kINELg0, (double)fCent});
      if (IsVtxInZCut) {
        FillTHnSparse("Normalisation", {(int)kINELg010, (double)fCent});
      }
      if (IsSelectedTrig) {
        FillTHnSparse("Normalisation", {(int)kINELg0_trig, (double)fCent});
        if (IsGoodVertex) {
          FillTHnSparse("Normalisation",
                        {(int)kINELg0_trig_vtx, (double)fCent});
          if (IsVtxInZCut) {
            FillTHnSparse("Normalisation",
                          {(int)kINELg0_trig_vtx10, (double)fCent});
          }
        }
      }
    }
    if (IsVtxInZCut) {
      FillTHnSparse("Normalisation", {(int)kINEL10, (double)fCent});
    }
    if (IsSelectedTrig) {
      FillTHnSparse("Normalisation", {(int)kINEL_trig, (double)fCent});
      if (IsGoodVertex) {
        FillTHnSparse("Normalisation", {(int)kINEL_trig_vtx, (double)fCent});
        if (IsVtxInZCut) {
          FillTHnSparse("Normalisation",
                        {(int)kINEL_trig_vtx10, (double)fCent});
        }
      }
    }
  }

  if (!IsEvtSelected) {
    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fNtupleSigma1385);
    return;  // event cut
  }

  fHistos->FillTH1("hMultiplicity", (double)fCent);

  if (fIsMC) {
    FillMCinput(fMCEvent);
    FillTHnSparse("Normalisation", {(int)kSelected, (double)fCent});
  }
  if (fIsAOD)
    fVertex = ((AliAODEvent*)fEvt)->GetPrimaryVertex();
  const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
  fPosPV[0] = pVtx->GetX();
  fPosPV[1] = pVtx->GetY();
  fPosPV[2] = pVtx->GetZ();

  // Event Mixing pool -----------------------------------------------------
  fZbin = fBinZ.FindBin(fPosPV[2]) - 1;    // Event mixing z-bin
  fCentBin = fBinCent.FindBin(fCent) - 1;  // Event mixing cent bin
  if (fIsINEL)
    fCentBin = 0;  // for INEL case

  bool checkPion = GoodTracksSelection();
  bool checkV0 = GoodV0Selection();

  if (checkPion && checkV0) {
    FillTracks();  // Fill the histogram
    if (fFillnTuple)
      FillNtuples();
  }
  if (fSetMixing && fGoodTrackArray.size()) {
    FillTrackToEventPool();  // use only pion track pool.
  }
  PostData(1, fHistos->GetListOfHistograms());
  PostData(2, fNtupleSigma1385);
}
//_____________________________________________________________________________
void AliAnalysisTaskSigma1385PM::Terminate(Option_t*) {}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSigma1385PM::GoodTracksSelection() {
  const UInt_t nTracks = fEvt->GetNumberOfTracks();
  fGoodTrackArray.clear();
  AliVTrack* track;
  Float_t b[2];
  Float_t bCov[3];
  Double_t nTPCNSigPion, pionZ, pionPt, pionSigmaDCA_r, pionDCA_r, lEta;

  for (UInt_t it = 0; it < nTracks; it++) {
    track = (AliVTrack*)fEvt->GetTrack(it);
    if (!track)
      continue;

    // ---------- Track selection begin ----------
    if (!fIsAOD) {
      if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
        continue;
    }  // ESD Case
    else {
      if (!fIsNano && !((AliAODTrack*)track)->TestFilterBit(fFilterBit))
        continue;
    }  // AOD Case

    GetImpactParam(track, b, bCov);

    pionZ = b[1];
    nTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
    pionPt = track->Pt();
    pionSigmaDCA_r = (0.0026 + 0.0050 / pionPt);
    pionDCA_r = b[0];
    lEta = track->Eta();
    if (fFillQAPlot) {
      fHistos->FillTH1("QA/hDCAPVPion", pionZ);
      fHistos->FillTH1("QA/hDCArPVPion", pionDCA_r);
      fHistos->FillTH1("QA/hEtaPion", lEta);
      fHistos->FillTH1("QA/hPtPion", pionPt);
      fHistos->FillTH2("QA/hTPCPIDPion", track->GetTPCmomentum(),
                       track->GetTPCsignal());
    }

    if (TMath::Abs(nTPCNSigPion) > fTPCNsigSigmaStarPionCut)
      continue;
    if (TMath::Abs(lEta) > fSigmaStarPionEtaCut)
      continue;
    if (pionPt < 0.15)
      continue;
    if (pionZ > fSigmaStarPionZVertexCut)
      continue;
    if (pionDCA_r > pionSigmaDCA_r * fSigmaStarPionXYVertexSigmaCut)
      continue;

    if (fFillQAPlot) {
      fHistos->FillTH1("QAcut/hDCAPVPion", pionZ);
      fHistos->FillTH1("QAcut/hDCArPVPion", pionDCA_r);
      fHistos->FillTH1("QAcut/hEtaPion", lEta);
      fHistos->FillTH1("QAcut/hPtPion", pionPt);
      fHistos->FillTH2("QAcut/hTPCPIDPion", track->GetTPCmomentum(),
                       track->GetTPCsignal());
    }

    fGoodTrackArray.push_back(it);
  }
  return fGoodTrackArray.size();
}
Bool_t AliAnalysisTaskSigma1385PM::GoodV0Selection() {
  fGoodV0Array.clear();
  const UInt_t nV0 = fEvt->GetNumberOfV0s();

  AliESDv0* v0ESD;
  AliAODv0* v0AOD;
  Bool_t lPIDLambda, lPIDAntiLambda;
  Double_t lDCADist_LambdaProton_PV, lDCADist_LambdaPion_PV;
  Double_t lDCADistLambda, lDCADistLambda_PV, lLambdaCPA;
  Double_t nTPCNSigProton, nTPCNSigAntiProton, nTPCNSigPion, nTPCNSigAntiPion;
  Double_t v0Position[3];
  Float_t impParPos[2], impParNeg[2];
  Float_t covMatPos[3], covMatNeg[3];
  Double_t fMassV0;
  Double_t radius, lV0TotalMomentum, lLength, lLifetime;
  UInt_t isAnti = 0;
  UInt_t isAntiCheck = 0;

  Bool_t AcceptedV0 = kTRUE;
  if (!fIsAOD) {  // ESD case
    for (UInt_t it = 0; it < nV0; it++) {
      lPIDLambda = kFALSE;
      lPIDAntiLambda = kFALSE;
      AcceptedV0 = kTRUE;
      isAntiCheck = 0;
      v0ESD = ((AliESDEvent*)fEvt)->GetV0(it);
      if (!v0ESD)
        continue;

      if (TMath::Abs(v0ESD->GetPindex()) == TMath::Abs(v0ESD->GetNindex()))
        continue;

      AliESDtrack* pTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
      AliESDtrack* nTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));

      // filter like-sign V0
      if (TMath::Abs(((pTrackV0->GetSign()) - (nTrackV0->GetSign()))) < 0.1)
        AcceptedV0 = kFALSE;

      // PID cuts
      nTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
      nTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
      nTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
      nTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);
      
      if ((TMath::Abs(nTPCNSigProton) > fTPCNsigLambdaProtonCut) &&
              (TMath::Abs(nTPCNSigPion) > fTPCNsigLambdaPionCut) &&
          (TMath::Abs(nTPCNSigAntiProton) > fTPCNsigLambdaProtonCut) &&
              (TMath::Abs(nTPCNSigAntiPion) > fTPCNsigLambdaPionCut))
        AcceptedV0 = kFALSE;

      if (nTrackV0->GetSign() < 0) {
        lPIDLambda = kTRUE;
        v0ESD->ChangeMassHypothesis(kLambda0);
        isAnti = 0;
        if (fFillQAPlot) {
          fHistos->FillTH2("QA/hTPCPIDLambdaProton", pTrackV0->GetTPCmomentum(),
                           pTrackV0->GetTPCsignal());
          fHistos->FillTH2("QA/hTPCPIDLambdaPion", nTrackV0->GetTPCmomentum(),
                           nTrackV0->GetTPCsignal());
        }
      } else if (nTrackV0->GetSign() > 0) {
        v0ESD->ChangeMassHypothesis(kLambda0Bar);
        lPIDAntiLambda = kTRUE;
        isAnti = 1;
        if (fFillQAPlot) {
          fHistos->FillTH2("QA/hTPCPIDAntiLambdaProton",
                           nTrackV0->GetTPCmomentum(),
                           nTrackV0->GetTPCsignal());
          fHistos->FillTH2("QA/hTPCPIDAntiLambdaPion",
                           pTrackV0->GetTPCmomentum(),
                           pTrackV0->GetTPCsignal());
        }
      } else
        AcceptedV0 = kFALSE;

      if (nTrackV0->GetSign() < 0)
        isAntiCheck += 1;
      if (nTrackV0->GetSign() > 0)
        isAntiCheck += 2;
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLambdaAntiCheck", isAntiCheck);

      // DCA cut
      // DCA between Dautgher particles
      lDCADistLambda = TMath::Abs(v0ESD->GetDcaV0Daughters());
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCA_lambdaDaughters", lDCADistLambda);
      if (lDCADistLambda > fDCADistLambdaDaughtersCut)
        AcceptedV0 = kFALSE;  // DCA proton-pion

      // DCA to PV V0
      lDCADistLambda_PV =
          TMath::Abs(v0ESD->GetD(fPosPV[0], fPosPV[1], fPosPV[2]));
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCAlambdaPV", lDCADistLambda_PV);
      if (lDCADistLambda_PV > fDCArDistLambdaPVCut)
        AcceptedV0 = kFALSE;

      // DCA to PV daughter
      // Adopted from AliRsnCutV0 fCustomTrackDCACuts
      GetImpactParam(pTrackV0, impParPos, covMatPos);
      GetImpactParam(nTrackV0, impParNeg, covMatNeg);
      if (covMatPos[0] <= 0 || covMatPos[2] <= 0)
        covMatPos[0] = 0;
      covMatPos[2] = 0;
      if (covMatNeg[0] <= 0 || covMatNeg[2] <= 0)
        covMatNeg[0] = 0;
      covMatNeg[2] = 0;
      if (!isAnti) {
        lDCADist_LambdaProton_PV = impParPos[0];
        lDCADist_LambdaPion_PV = impParNeg[0];
      } else {
        lDCADist_LambdaPion_PV = impParPos[0];
        lDCADist_LambdaProton_PV = impParNeg[0];
      }
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCALambdaPVProton", lDCADist_LambdaProton_PV);
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCALambdaPVPion", lDCADist_LambdaPion_PV);
      if (TMath::Abs(lDCADist_LambdaProton_PV) < fDCAPositiveTrack)
        AcceptedV0 = kFALSE;
      if (TMath::Abs(lDCADist_LambdaPion_PV) < fDCANegativeTrack)
        AcceptedV0 = kFALSE;

      // CPA cut
      lLambdaCPA = TMath::Abs(
          v0ESD->GetV0CosineOfPointingAngle(fPosPV[0], fPosPV[1], fPosPV[2]));
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hCosPAlambda", lLambdaCPA);
      if (lLambdaCPA < fV0CosineOfPointingAngleCut || lLambdaCPA >= 1)
        AcceptedV0 = kFALSE;

      // Rapidity cut
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hYlambda", v0ESD->RapLambda());
      if (TMath::Abs(v0ESD->RapLambda()) > fMaxLambdaRapidity)
        AcceptedV0 = kFALSE;

      // Radius cut
      v0ESD->GetXYZ(v0Position[0], v0Position[1], v0Position[2]);
      radius = TMath::Hypot(v0Position[0], v0Position[1]);
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLambdaRxy", radius);
      if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
        AcceptedV0 = kFALSE;

      // Life time cut
      lV0TotalMomentum = v0ESD->P();
      lLength = TMath::Sqrt(TMath::Power(v0Position[0] - fPosPV[0], 2) +
                            TMath::Power(v0Position[1] - fPosPV[1], 2) +
                            TMath::Power(v0Position[2] - fPosPV[2], 2));
      lLifetime = TMath::Abs(v0Mass * lLength / lV0TotalMomentum);
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLifetimelambda", lLifetime);
      if (lLifetime > fLambdaLifetime)
        AcceptedV0 = kFALSE;

      // Mass window cut
      fMassV0 = v0ESD->GetEffMass();
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hMassLambda", fMassV0);
      if (TMath::Abs(fMassV0 - v0Mass) > fV0MassWindowCut)
        AcceptedV0 = kFALSE;

      // After selection above
      if (AcceptedV0) {
        fGoodV0Array.push_back({it, isAnti});  // for standard V0
        if (fFillQAPlot) {
          if (lPIDLambda) {
            fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                             pTrackV0->GetTPCmomentum(),
                             pTrackV0->GetTPCsignal());
            fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                             nTrackV0->GetTPCmomentum(),
                             nTrackV0->GetTPCsignal());
          }
          if (lPIDAntiLambda) {
            fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaProton",
                             nTrackV0->GetTPCmomentum(),
                             nTrackV0->GetTPCsignal());
            fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaPion",
                             pTrackV0->GetTPCmomentum(),
                             pTrackV0->GetTPCsignal());
          }
          fHistos->FillTH1("QAcut/hDCA_lambdaDaughters", lDCADistLambda);
          fHistos->FillTH1("QAcut/hDCAlambdaPV", lDCADistLambda_PV);
          fHistos->FillTH1("QAcut/hDCALambdaPVProton",
                           lDCADist_LambdaProton_PV);
          fHistos->FillTH1("QAcut/hDCALambdaPVPion", lDCADist_LambdaPion_PV);
          fHistos->FillTH1("QAcut/hCosPAlambda", lLambdaCPA);
          fHistos->FillTH1("QAcut/hLifetimelambda", lLifetime);
          fHistos->FillTH1("QAcut/hYlambda", v0ESD->RapLambda());
          fHistos->FillTH1("QAcut/hLambdaRxy", radius);
          fHistos->FillTH1("QAcut/hMassLambda", fMassV0);
          fHistos->FillTH1("QAcut/hLambdaAntiCheck", isAntiCheck);
        }
      }
    }  // All V0 loop
  }    // ESD case
  else {
    for (UInt_t it = 0; it < nV0; it++) {
      lPIDLambda = kFALSE;
      lPIDAntiLambda = kFALSE;
      AcceptedV0 = kTRUE;
      isAntiCheck = 0;
      v0AOD = ((AliAODEvent*)fEvt)->GetV0(it);
      if (!v0AOD)
        continue;

      if (TMath::Abs(v0AOD->GetPosID()) == TMath::Abs(v0AOD->GetNegID()))
        continue;

      AliAODTrack* pTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack* nTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

      // filter like-sign V0
      if (TMath::Abs(((pTrackV0->GetSign()) - (nTrackV0->GetSign()))) < 0.1)
        AcceptedV0 = kFALSE;

      // PID cuts
      nTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
      nTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
      nTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
      nTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

      if ((TMath::Abs(nTPCNSigProton) > fTPCNsigLambdaProtonCut) &&
              (TMath::Abs(nTPCNSigPion) > fTPCNsigLambdaPionCut) &&
          (TMath::Abs(nTPCNSigAntiProton) > fTPCNsigLambdaProtonCut) &&
              (TMath::Abs(nTPCNSigAntiPion) > fTPCNsigLambdaPionCut))
        AcceptedV0 = kFALSE;

      if (nTrackV0->GetSign() < 0) {
        lPIDLambda = kTRUE;
        isAnti = 0;
        if (fFillQAPlot) {
          fHistos->FillTH2("QA/hTPCPIDLambdaProton", pTrackV0->GetTPCmomentum(),
                           pTrackV0->GetTPCsignal());
          fHistos->FillTH2("QA/hTPCPIDLambdaPion", nTrackV0->GetTPCmomentum(),
                           nTrackV0->GetTPCsignal());
        }
      } else if (nTrackV0->GetSign() > 0) {
        lPIDAntiLambda = kTRUE;
        isAnti = 1;
        if (fFillQAPlot) {
          fHistos->FillTH2("QA/hTPCPIDAntiLambdaProton",
                           nTrackV0->GetTPCmomentum(),
                           nTrackV0->GetTPCsignal());
          fHistos->FillTH2("QA/hTPCPIDAntiLambdaPion",
                           pTrackV0->GetTPCmomentum(),
                           pTrackV0->GetTPCsignal());
        }
      } else
        AcceptedV0 = kFALSE;
      if (nTrackV0->GetSign() < 0)
        isAntiCheck += 1;
      if (nTrackV0->GetSign() > 0)
        isAntiCheck += 2;
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLambdaAntiCheck", isAntiCheck);

      // DCA cut
      lDCADistLambda = TMath::Abs(v0AOD->DcaV0Daughters());
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCA_lambdaDaughters", lDCADistLambda);
      if (lDCADistLambda > fDCADistLambdaDaughtersCut)
        AcceptedV0 = kFALSE;  // DCA proton-pion

      // DCA to PV daughter
      // Adopted from AliRsnCutV0 fCustomTrackDCACuts
      if (!isAnti) {
        lDCADist_LambdaProton_PV = v0AOD->DcaPosToPrimVertex();
        lDCADist_LambdaPion_PV = v0AOD->DcaNegToPrimVertex();
      } else {
        lDCADist_LambdaPion_PV = v0AOD->DcaPosToPrimVertex();
        lDCADist_LambdaProton_PV = v0AOD->DcaNegToPrimVertex();
      }
      if (TMath::Abs(lDCADist_LambdaProton_PV) < fDCAPositiveTrack)
        AcceptedV0 = kFALSE;
      if (TMath::Abs(lDCADist_LambdaPion_PV) < fDCANegativeTrack)
        AcceptedV0 = kFALSE;
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCALambdaPVProton", lDCADist_LambdaProton_PV);
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCALambdaPVPion", lDCADist_LambdaPion_PV);

      // DCA to PV
      lDCADistLambda_PV = TMath::Abs(v0AOD->DcaV0ToPrimVertex());
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hDCAlambdaPV", lDCADistLambda_PV);
      if (lDCADistLambda_PV > fDCArDistLambdaPVCut)
        AcceptedV0 = kFALSE;

      // CPA cut
      lLambdaCPA = TMath::Abs(v0AOD->CosPointingAngle(fVertex));
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hCosPAlambda", lLambdaCPA);
      if (lLambdaCPA < fV0CosineOfPointingAngleCut || lLambdaCPA >= 1)
        AcceptedV0 = kFALSE;

      // Rapidity cut
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hYlambda", v0AOD->RapLambda());
      if (TMath::Abs(v0AOD->RapLambda()) > fMaxLambdaRapidity)
        AcceptedV0 = kFALSE;

      // Radius cut
      radius = v0AOD->RadiusV0();
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLambdaRxy", radius);
      if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
        AcceptedV0 = kFALSE;

      // Life time cut
      lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
      lLength = v0AOD->DecayLength(fPosPV);
      lLifetime = TMath::Abs(v0Mass * lLength / lV0TotalMomentum);
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hLifetimelambda", lLifetime);
      if (lLifetime > fLambdaLifetime)
        AcceptedV0 = kFALSE;

      // Mass window cut
      fMassV0 = -999;
      if (lPIDLambda)
        fMassV0 = v0AOD->MassLambda();
      else if (lPIDAntiLambda)
        fMassV0 = v0AOD->MassAntiLambda();
      if (fFillQAPlot)
        fHistos->FillTH1("QA/hMassLambda", fMassV0);
      if (TMath::Abs(fMassV0 - v0Mass) > fV0MassWindowCut)
        AcceptedV0 = kFALSE;

      // After selection above
      if (AcceptedV0) {
        fGoodV0Array.push_back({it, isAnti});
        if (fFillQAPlot) {
          if (lPIDLambda) {
            fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                             pTrackV0->GetTPCmomentum(),
                             pTrackV0->GetTPCsignal());
            fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                             nTrackV0->GetTPCmomentum(),
                             nTrackV0->GetTPCsignal());
          }
          if (lPIDAntiLambda) {
            fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaProton",
                             nTrackV0->GetTPCmomentum(),
                             nTrackV0->GetTPCsignal());
            fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaPion",
                             pTrackV0->GetTPCmomentum(),
                             pTrackV0->GetTPCsignal());
          }
          fHistos->FillTH1("QAcut/hDCA_lambdaDaughters", lDCADistLambda);
          fHistos->FillTH1("QAcut/hDCAlambdaPV", lDCADistLambda_PV);
          fHistos->FillTH1("QAcut/hDCALambdaPVProton",
                           lDCADist_LambdaProton_PV);
          fHistos->FillTH1("QAcut/hDCALambdaPVPion", lDCADist_LambdaPion_PV);
          fHistos->FillTH1("QAcut/hCosPAlambda", lLambdaCPA);
          fHistos->FillTH1("QAcut/hLifetimelambda", lLifetime);
          fHistos->FillTH1("QAcut/hYlambda", v0AOD->RapLambda());
          fHistos->FillTH1("QAcut/hLambdaRxy", radius);
          fHistos->FillTH1("QAcut/hMassLambda", fMassV0);
          fHistos->FillTH1("QAcut/hLambdaAntiCheck", isAntiCheck);
        }
      }
    }  // All v0 loop
  }    // AOD case

  return fGoodV0Array.size();
}
void AliAnalysisTaskSigma1385PM::FillTracks() {
  AliVTrack* track1;
  AliVTrack* track_mix;
  AliESDv0* v0ESD;
  AliAODv0* v0AOD;
  Bool_t isAnti, isPionPlus;
  Bool_t SkipMixing = kFALSE;
  Int_t pID, nID;
  int sign = kAllType;
  int binAnti = 0;
  double asym = 0.;

  TLorentzVector vecLambda, vecPion, vecPionMix;
  TLorentzVector vecSigmaStar;
  const UInt_t nV0 = fGoodV0Array.size();
  const UInt_t nTracks = fGoodTrackArray.size();

  tracklist trackpool;
  if (fSetMixing) {
    eventpool& ep = fEMpool[fCentBin][fZbin];
    if ((int)ep.size() < (int)fnMix)
      SkipMixing = kTRUE;
    if (!SkipMixing) {
      for (auto pool : ep) {
        for (auto track : pool)
          trackpool.push_back((AliVTrack*)track);
      }
    }
  }
  for (UInt_t i = 0; i < nV0; i++) {
    if (!fIsAOD) {
      v0ESD = ((AliESDEvent*)fEvt)->GetV0(fGoodV0Array[i][0]);
      if (!v0ESD)
        continue;
      AliESDtrack* pTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
      AliESDtrack* nTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
      pID = pTrackV0->GetID();
      nID = nTrackV0->GetID();

      if (fGoodV0Array[i][1] > 0)
        isAnti = true;
      else
        isAnti = false;

      if (!isAnti)
        v0ESD->ChangeMassHypothesis(kLambda0);
      else
        v0ESD->ChangeMassHypothesis(kLambda0Bar);

      vecLambda.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                        v0ESD->GetEffMass());
    } else {
      v0AOD = ((AliAODEvent*)fEvt)->GetV0(fGoodV0Array[i][0]);
      if (!v0AOD)
        continue;
      AliAODTrack* pTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack* nTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
      pID = pTrackV0->GetID();
      nID = nTrackV0->GetID();

      if (fGoodV0Array[i][1] > 0)
        isAnti = true;
      else
        isAnti = false;

      if (!isAnti)
        vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                          v0AOD->MassLambda());
      else
        vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                          v0AOD->MassAntiLambda());
    }

    for (UInt_t j = 0; j < nTracks; j++) {
      track1 = (AliVTrack*)fEvt->GetTrack(fGoodTrackArray[j]);
      if (!track1)
        continue;

      if (track1->GetID() == pID || track1->GetID() == nID)
        continue;

      vecPion.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

      vecSigmaStar = vecLambda + vecPion;  // vecLambda = cascade, vecPion=pion
      // Y cut
      if ((vecSigmaStar.Rapidity() > fSigmaStarYCutHigh) ||
          (vecSigmaStar.Rapidity() < fSigmaStarYCutLow))
        continue;

      // AsymmCut
      if (fUseAsymmCut) {
        auto P1 = (TVector3)vecLambda.Vect();
        auto P2 = (TVector3)vecPion.Vect();
        asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if (fFillQAPlot)
          fHistos->FillTH2("QA/hSigmaStarAsymm", asym, fCent);
        if ((asym < fSigmaStarAsymmCutLow) || (asym > fSigmaStarAsymmCutHigh))
          continue;
        if (fFillQAPlot)
          fHistos->FillTH2("QAcut/hSigmaStarAsymm", asym, fCent);
      }

      if (track1->Charge() > 0)
        isPionPlus = true;
      else
        isPionPlus = false;

      (isAnti) ? binAnti = kAnti : binAnti = kNormal;
      (isPionPlus) ? sign = kSigmaStarP : sign = kSigmaStarN;

      FillTHnSparse("Sigma1385_data",
                    {(double)binAnti, (double)sign, (double)fCent,
                     vecSigmaStar.Pt(), vecSigmaStar.M()});

      if (fIsMC) {
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 0)) {
          (isPionPlus) ? sign = kSigmaStarP_REC : sign = kSigmaStarN_REC;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
          if (fUseAsymmCut && fFillQAPlot)
            fHistos->FillTH2("QAcut/hSigmaStarAsymm_true_selected", asym, fCent);
        } else {
          // MC not true bkg
          (isPionPlus) ? sign = kSigmaStarP_NOT : sign = kSigmaStarN_NOT;
          FillTHnSparse("Sigma1385_data",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 1)) {
          (isPionPlus) ? sign = kSigmaStarPBkg_type1
                       : sign = kSigmaStarNBkg_type1;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 2)) {
          (isPionPlus) ? sign = kSigmaStarPBkg_type2
                       : sign = kSigmaStarNBkg_type2;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 3)) {
          (isPionPlus) ? sign = kSigmaStarPBkg_type3
                       : sign = kSigmaStarNBkg_type3;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 4)) {
          (isPionPlus) ? sign = kSigmaStarPBkg_type4
                       : sign = kSigmaStarNBkg_type4;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j], 5)) {
          (isPionPlus) ? sign = kSigmaStarPBkg_type5
                       : sign = kSigmaStarNBkg_type5;
          FillTHnSparse("Sigma1385_mc",
                        {(double)binAnti, (double)sign, (double)fCent,
                         vecSigmaStar.Pt(), vecSigmaStar.M()});
        }
      }
    }  // pion loop

    if (fSetMixing && !SkipMixing && (fCentBin >= 0) && (fZbin >= 0)) {
      for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
        track_mix = trackpool.at(jt);
        if (track_mix->GetID() == pID || track_mix->GetID() == nID)
          continue;
        vecPionMix.SetXYZM(track_mix->Px(), track_mix->Py(), track_mix->Pz(),
                           pionMass);
        vecSigmaStar = vecLambda + vecPionMix;
        // Y cut
        if ((vecSigmaStar.Rapidity() > fSigmaStarYCutHigh) ||
            (vecSigmaStar.Rapidity() < fSigmaStarYCutLow))
          continue;
        // AsymmCut
        if (fUseAsymmCut) {
          auto P1 = (TVector3)vecLambda.Vect();
          auto P2 = (TVector3)vecPionMix.Vect();
          auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
          if ((asym < fSigmaStarAsymmCutLow) || (asym > fSigmaStarAsymmCutHigh))
            continue;
        }

        if (track_mix->Charge() > 0)
          isPionPlus = true;
        else
          isPionPlus = false;

        if (fGoodV0Array[i][1] > 0)
          isAnti = true;
        else
          isAnti = false;

        (isAnti) ? binAnti = kAnti : binAnti = kNormal;
        (isPionPlus) ? sign = kSigmaStarP_MIX : sign = kSigmaStarN_MIX;

        FillTHnSparse("Sigma1385_data",
                      {(double)binAnti, (double)sign, (double)fCent,
                       vecSigmaStar.Pt(), vecSigmaStar.M()});
      }
    }
  }
}
void AliAnalysisTaskSigma1385PM::FillNtuples() {
  AliVTrack* track1 = nullptr;
  AliESDv0* v0ESD = nullptr;
  AliAODv0* v0AOD = nullptr;
  Double_t nTPCNSigProton, nTPCNSigAntiProton, nTPCNSigPion, nTPCNSigAntiPion;
  Bool_t isAnti, isPionPlus;
  Int_t pID, nID;
  Double_t tmp[17];
  int sign = kAllType;
  int binAnti = 0;

  for (UInt_t i = 0; i < 16; i++)
    tmp[i] = -999;  // initial value

  TLorentzVector vecLambda, vecPion;
  TLorentzVector vecSigmaStar;
  const UInt_t nV0 = fGoodV0Array.size();
  const UInt_t nTracks = fGoodTrackArray.size();
  for (UInt_t i = 0; i < nV0; i++) {
    if (!fIsAOD) {
      v0ESD = ((AliESDEvent*)fEvt)->GetV0(fGoodV0Array[i][0]);
      if (!v0ESD)
        continue;
      AliESDtrack* pTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
      AliESDtrack* nTrackV0 =
          ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));
      pID = pTrackV0->GetID();
      nID = nTrackV0->GetID();

      if (fGoodV0Array[i][1] > 0)
        isAnti = true;
      else
        isAnti = false;

      if (!isAnti)
        v0ESD->ChangeMassHypothesis(kLambda0);
      else
        v0ESD->ChangeMassHypothesis(kLambda0Bar);

      vecLambda.SetXYZM(v0ESD->Px(), v0ESD->Py(), v0ESD->Pz(),
                        v0ESD->GetEffMass());

      nTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
      nTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
      nTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
      nTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

      // nTuple
      tmp[5] = nTPCNSigProton;                          // PIDV0pTrackProton
      tmp[6] = nTPCNSigAntiProton;                      // PIDV0pTrackPion
      tmp[7] = nTPCNSigPion;                            // PIDV0nTrackProton
      tmp[8] = nTPCNSigAntiPion;                        // PIDV0nTrackPion
      tmp[9] = TMath::Abs(v0ESD->GetDcaV0Daughters());  // DCAV0Daughters
      tmp[10] = TMath::Abs(
          v0ESD->GetD(fPosPV[0], fPosPV[1], fPosPV[2]));  // DCAV0ToPrimVertex
      tmp[11] = v0ESD->GetV0CosineOfPointingAngle(
          fPosPV[0], fPosPV[1], fPosPV[2]);  // CosPointingAngleV0
      tmp[12] = v0ESD->GetEffMass();         // V0Mass
      tmp[13] = v0ESD->Eta();                // EtaV0
      tmp[14] = v0ESD->Phi();                // PhiV0
    } else {
      v0AOD = ((AliAODEvent*)fEvt)->GetV0(fGoodV0Array[i][0]);
      if (!v0AOD)
        continue;
      AliAODTrack* pTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
      AliAODTrack* nTrackV0 =
          (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
      pID = pTrackV0->GetID();
      nID = nTrackV0->GetID();

      if (fGoodV0Array[i][1] > 0)
        isAnti = true;
      else
        isAnti = false;

      if (!isAnti)
        vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                          v0AOD->MassLambda());
      else
        vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(),
                          v0AOD->MassAntiLambda());

      nTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
      nTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
      nTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
      nTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

      // nTuple
      tmp[5] = nTPCNSigProton;                           // PIDV0pTrackProton
      tmp[6] = nTPCNSigAntiProton;                       // PIDV0pTrackPion
      tmp[7] = nTPCNSigPion;                             // PIDV0nTrackProton
      tmp[8] = nTPCNSigAntiPion;                         // PIDV0nTrackPion
      tmp[9] = TMath::Abs(v0AOD->DcaV0Daughters());      // DCAV0Daughters
      tmp[10] = TMath::Abs(v0AOD->DcaV0ToPrimVertex());  // DCAV0ToPrimVertex
      tmp[11] = v0AOD->CosPointingAngle(fVertex);        // CosPointingAngleV0
      if (!isAnti)
        tmp[12] = v0AOD->MassLambda();  // V0Mass
      else
        tmp[12] = v0AOD->MassAntiLambda();  // V0Mass
      tmp[13] = v0AOD->Eta();               // EtaV0
      tmp[14] = v0AOD->Phi();               // PhiV0
    }

    for (UInt_t j = 0; j < nTracks; j++) {
      track1 = (AliVTrack*)fEvt->GetTrack(fGoodTrackArray[j]);
      if (!track1)
        continue;

      if (track1->GetID() == pID || track1->GetID() == nID)
        continue;

      vecPion.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

      vecSigmaStar = vecLambda + vecPion;  // vecLambda = cascade, vecPion=pion
      // Y cut
      if ((vecSigmaStar.Rapidity() > fSigmaStarYCutHigh) ||
          (vecSigmaStar.Rapidity() < fSigmaStarYCutLow))
        continue;
      // AsymmCut
      if (fUseAsymmCut) {
        auto P1 = (TVector3)vecLambda.Vect();
        auto P2 = (TVector3)vecPion.Vect();
        auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if ((asym < fSigmaStarAsymmCutLow) || (asym > fSigmaStarAsymmCutHigh))
          continue;
      }

      if (track1->Charge() > 0)
        isPionPlus = true;
      else
        isPionPlus = false;

      if (fGoodV0Array[i][1] > 0)
        isAnti = true;
      else
        isAnti = false;

      (isAnti) ? binAnti = kAnti : binAnti = kNormal;
      (isPionPlus) ? sign = kSigmaStarP : sign = kSigmaStarN;

      tmp[0] = GetTPCnSigma(track1, AliPID::kPion);  // PIDSigmaStarPion
      tmp[1] = TMath::Abs(track1->GetZ() -
                          fPosPV[2]);  // DCASigmaStarPionToPrimVertexZ
      tmp[2] = sqrt(pow(track1->GetX(), 2) +
                    pow(track1->GetY(), 2));  // DCASigmaStarPionToPrimVertexR
      tmp[3] = track1->Eta();                 // EtaSigmaStarPion
      tmp[4] = track1->Phi();                 // PhiSigmaStarPion
      // nTuple -> See above
      // tmp[5] = nTPCNSigProton;//PIDV0pTrackProton
      // tmp[6] = nTPCNSigAntiProton;//PIDV0pTrackPion
      // tmp[7] = nTPCNSigPion;//PIDV0nTrackProton
      // tmp[8] = nTPCNSigAntiPion;//PIDV0nTrackPion
      // tmp[9] = TMath::Abs(v0ESD->GetDcaV0Daughters()); //DCAV0Daughters
      // tmp[10] = TMath::Abs(v0ESD->GetD(fPosPV[0], fPosPV[1], fPosPV[2]));
      // //DCAV0ToPrimVertex tmp[11] =
      // v0ESD->GetV0CosineOfPointingAngle(fPosPV[0], fPosPV[1],
      // fPosPV[2]);//CosPointingAngleV0 tmp[12] = v0ESD->GetEffMass(); //V0Mass
      // tmp[13] = v0ESD->Eta(); //EtaV0
      // tmp[14] = v0ESD->Phi(); //PhiV0

      if (fIsMC) {
        if (IsTrueSigmaStar(fGoodV0Array[i][0], fGoodTrackArray[j]))
          tmp[15] = (double)sign;  // MCflag
        else
          tmp[15] = 3;  // MCflag -> not true
      } else
        tmp[15] = 0;  // MCflag -> data
      tmp[16] = binAnti;
      fNtupleSigma1385->Fill(tmp);
    }  // pion loop
  }
}
void AliAnalysisTaskSigma1385PM::FillMCinput(AliMCEvent* fMCEvent,
                                             int Fillbin) {
  int sign = kAllType;
  int binAnti = 0;
  TLorentzVector vecPart1, vecPart2;
  if (!fIsAOD) {
    for (Int_t it = 0; it < fMCEvent->GetNumberOfPrimaries(); it++) {
      TParticle* mcInputTrack = (TParticle*)fMCEvent->GetTrack(it)->Particle();
      if (!mcInputTrack) {
        Error("UserExec", "Could not receive MC track %d", it);
        continue;
      }
      Int_t v0PdgCode = mcInputTrack->GetPdgCode();
      if ((TMath::Abs(v0PdgCode) != kSigmaStarPCode) &&
          (TMath::Abs(v0PdgCode) != kSigmaStarNCode))
        continue;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;
      // Y cut
      if ((mcInputTrack->Y() > fSigmaStarYCutHigh) ||
          (mcInputTrack->Y() < fSigmaStarYCutLow))
        continue;
      // AsymmCut
      if (fUseAsymmCut) {
        auto mcPart1 =
            (TParticle*)fMCEvent->GetTrack(abs(mcInputTrack->GetDaughter(0)))
                ->Particle();
        auto mcPart2 =
            (TParticle*)fMCEvent->GetTrack(abs(mcInputTrack->GetDaughter(1)))
                ->Particle();
        mcPart1->Momentum(vecPart1);
        mcPart2->Momentum(vecPart2);

        auto P1 = (TVector3)vecPart1.Vect();
        auto P2 = (TVector3)vecPart2.Vect();
        auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if (fFillQAPlot)
          fHistos->FillTH2("QA/hSigmaStarAsymm_true", asym, fCent);
        if ((asym < fSigmaStarAsymmCutLow) || (asym > fSigmaStarAsymmCutHigh))
          continue;
        if (fFillQAPlot)
          fHistos->FillTH2("QAcut/hSigmaStarAsymm_true", asym, fCent);
      }
      (v0PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
      if (TMath::Abs(v0PdgCode) == kSigmaStarPCode)
        sign = kSigmaStarP_GEN + (int)Fillbin * 2;
      if (TMath::Abs(v0PdgCode) == kSigmaStarNCode)
        sign = kSigmaStarN_GEN + (int)Fillbin * 2;

      FillTHnSparse("Sigma1385_mc",
                    {(double)binAnti, (double)sign, (double)fCent,
                     mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
    }
  } else {
    for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++) {
      AliAODMCParticle* mcInputTrack = (AliAODMCParticle*)fMCArray->At(it);
      if (!mcInputTrack) {
        Error("UserExec", "Could not receive MC track %d", it);
        continue;
      }

      Int_t v0PdgCode = mcInputTrack->GetPdgCode();

      if ((TMath::Abs(v0PdgCode) != kSigmaStarPCode) &&
          (TMath::Abs(v0PdgCode) != kSigmaStarNCode))
        continue;
      if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
        continue;

      // Y cut
      if ((mcInputTrack->Y() > fSigmaStarYCutHigh) ||
          (mcInputTrack->Y() < fSigmaStarYCutLow))
        continue;
      // AsymmCut
      if (fUseAsymmCut) {
        auto mcPart1 = (AliAODMCParticle*)fMCArray->At(
            abs(mcInputTrack->GetDaughterFirst()));
        auto mcPart2 = (AliAODMCParticle*)fMCArray->At(
            abs(mcInputTrack->GetDaughterLast()));
        mcPart1->Momentum(vecPart1);
        mcPart2->Momentum(vecPart2);

        auto P1 = (TVector3)vecPart1.Vect();
        auto P2 = (TVector3)vecPart2.Vect();
        auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if (fFillQAPlot)
          fHistos->FillTH2("QA/hSigmaStarAsymm_true", asym, fCent);
        if ((asym < fSigmaStarAsymmCutLow) || (asym > fSigmaStarAsymmCutHigh))
          continue;
        if (fFillQAPlot)
          fHistos->FillTH2("QAcut/hSigmaStarAsymm_true", asym, fCent);
      }

      (v0PdgCode < 0) ? binAnti = kAnti : binAnti = kNormal;
      if (TMath::Abs(v0PdgCode) == kSigmaStarPCode)
        sign = kSigmaStarP_GEN + (int)Fillbin * 2;
      if (TMath::Abs(v0PdgCode) == kSigmaStarNCode)
        sign = kSigmaStarN_GEN + (int)Fillbin * 2;

      FillTHnSparse("Sigma1385_mc",
                    {(double)binAnti, (double)sign, (double)fCent,
                     mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
    }
  }
}
Bool_t AliAnalysisTaskSigma1385PM::IsTrueSigmaStar(UInt_t v0Index,
                                                   UInt_t pionIndex,
                                                   UInt_t BkgCheck) {
  AliVTrack* track1;
  AliESDv0* v0ESD;
  AliAODv0* v0AOD;

  track1 = (AliVTrack*)fEvt->GetTrack(pionIndex);

  if (!fIsAOD) {
    v0ESD = ((AliESDEvent*)fEvt)->GetV0(v0Index);
    if (!v0ESD)
      return kFALSE;
    AliESDtrack* pTrackV0 =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetPindex()));
    AliESDtrack* nTrackV0 =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(v0ESD->GetNindex()));

    TParticle* MCLamD1 =
        (TParticle*)fMCEvent->GetTrack(TMath::Abs(pTrackV0->GetLabel()))
            ->Particle();
    TParticle* MCLamD2 =
        (TParticle*)fMCEvent->GetTrack(TMath::Abs(nTrackV0->GetLabel()))
            ->Particle();
    TParticle* MCLam;
    TParticle* MCPi;
    TParticle* MCSigmaStar;

    // Aditional bgk estimation
    TParticle* MCSigmaZero;
    TParticle* MCXiStar;
    TParticle* MCXi;

    // Lambda daughter (pion, proton) check
    if ((TMath::Abs(MCLamD1->GetPdgCode()) == kProtonCode &&
         TMath::Abs(MCLamD2->GetPdgCode()) != kPionCode) ||
        (TMath::Abs(MCLamD1->GetPdgCode()) == kPionCode &&
         TMath::Abs(MCLamD2->GetPdgCode()) != kProtonCode))
      return kFALSE;
    // Lambda duather's mother check
    if (MCLamD1->GetMother(0) != MCLamD2->GetMother(0))
      return kFALSE;
    // Lambda check
    MCLam = (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLamD1->GetMother(0)))
                ->Particle();
    if (TMath::Abs(MCLam->GetPdgCode()) != kLambdaCode)
      return kFALSE;
    // pion check
    MCPi = (TParticle*)fMCEvent->GetTrack(TMath::Abs(track1->GetLabel()))
               ->Particle();
    if (TMath::Abs(MCPi->GetPdgCode()) != kPionCode)
      return kFALSE;

    switch (BkgCheck) {
      case 0:  // Normal Sigma1385 case
        // pion-Lambda mother check
        if (MCPi->GetMother(0) != MCLam->GetMother(0))
          return kFALSE;
        MCSigmaStar =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                ->Particle();
        if ((TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarPCode) &&
            (TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarNCode))
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 1:  // Bkg type A: from LambdaStar
        if (MCPi->GetMother(0) != MCLam->GetMother(0))
          return kFALSE;
        MCSigmaStar =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                ->Particle();
        if (TMath::Abs(MCSigmaStar->GetPdgCode()) != kLambdaStarCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 2:  // Bkg type B: from Sigma1385 through Sigma0
        // Sigma0 check
        MCSigmaZero =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                ->Particle();
        if (TMath::Abs(MCSigmaZero->GetPdgCode()) != kSigmaZeroCode)
          return kFALSE;
        // Sigma0-pion mother check
        if (MCPi->GetMother(0) != MCSigmaZero->GetMother(0))
          return kFALSE;
        // Sigma1385 check
        MCSigmaStar = (TParticle*)fMCEvent
                          ->GetTrack(TMath::Abs(MCSigmaZero->GetMother(0)))
                          ->Particle();
        if ((TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarPCode) &&
            (TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarNCode))
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 3:  // Bkg type C: Xi(1530)0
        // Xi- check
        MCXi = (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                   ->Particle();
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiCode)
          return kFALSE;
        // Xi-pion mother check
        if (MCPi->GetMother(0) != MCXi->GetMother(0))
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCXi->GetMother(0)))
                ->Particle();
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 4:  // Bkg type C: Xi(1530)- through Xi0
        // Xi- check
        MCXi = (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                   ->Particle();
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiZeroCode)
          return kFALSE;
        // Xi-pion mother check
        if (MCPi->GetMother(0) != MCXi->GetMother(0))
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCXi->GetMother(0)))
                ->Particle();
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarMCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 5:  // Bkg type C: Xi(1530)- through Xi-
               // Not sure that if it is needed since it will make normal Xi-
               // peak.
        // Xi- check
        MCXi = (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCLam->GetMother(0)))
                   ->Particle();
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiCode)
          return kFALSE;
        // Labmda-pion mother check to Xi-
        if (MCLam->GetMother(0) != MCXi->GetMother(0))
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(MCXi->GetMother(0)))
                ->Particle();
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarMCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
        return kFALSE;
    }
  } else {
    v0AOD = ((AliAODEvent*)fEvt)->GetV0(v0Index);
    if (!v0AOD)
      return kFALSE;
    AliAODTrack* pTrackV0 =
        (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
    AliAODTrack* nTrackV0 =
        (AliAODTrack*)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

    AliAODMCParticle* MCLamD1 =
        (AliAODMCParticle*)fMCArray->At(TMath::Abs(pTrackV0->GetLabel()));
    AliAODMCParticle* MCLamD2 =
        (AliAODMCParticle*)fMCArray->At(TMath::Abs(nTrackV0->GetLabel()));
    AliAODMCParticle* MCLam;
    AliAODMCParticle* MCPi;
    AliAODMCParticle* MCSigmaStar;

    // Aditional bgk estimation
    AliAODMCParticle* MCSigmaZero;
    AliAODMCParticle* MCXiStar;
    AliAODMCParticle* MCXi;

    // Lambda daughter (pion, proton) check
    if ((TMath::Abs(MCLamD1->GetPdgCode()) == kProtonCode &&
         TMath::Abs(MCLamD2->GetPdgCode()) != kPionCode) ||
        (TMath::Abs(MCLamD1->GetPdgCode()) == kPionCode &&
         TMath::Abs(MCLamD2->GetPdgCode()) != kProtonCode))
      return kFALSE;
    // Lambda duather's mother check
    if (MCLamD1->GetMother() != MCLamD2->GetMother())
      return kFALSE;
    // Lambda check
    MCLam = (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLamD1->GetMother()));
    if (TMath::Abs(MCLam->GetPdgCode()) != kLambdaCode)
      return kFALSE;
    // pion check
    MCPi = (AliAODMCParticle*)fMCArray->At(TMath::Abs(track1->GetLabel()));
    if (TMath::Abs(MCPi->GetPdgCode()) != kPionCode)
      return kFALSE;

    switch (BkgCheck) {
      case 0:  // Normal Sigma1385 case
        // pion-Lambda mother check
        if (MCPi->GetMother() != MCLam->GetMother())
          return kFALSE;
        MCSigmaStar =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if ((TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarPCode) &&
            (TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarNCode))
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 1:  // Bkg type A: from LambdaStar
        if (MCPi->GetMother() != MCLam->GetMother())
          return kFALSE;
        MCSigmaStar =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if (TMath::Abs(MCSigmaStar->GetPdgCode()) != kLambdaStarCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 2:  // Bkg type B: from Sigma1385 through Sigma0
        // Sigma0 check
        MCSigmaZero =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if (TMath::Abs(MCSigmaZero->GetPdgCode()) != kSigmaZeroCode)
          return kFALSE;
        // Sigma0-pion mother check
        if (MCPi->GetMother() != MCSigmaZero->GetMother())
          return kFALSE;
        // Sigma1385 check
        MCSigmaStar = (AliAODMCParticle*)fMCArray->At(
            TMath::Abs(MCSigmaZero->GetMother()));
        if ((TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarPCode) &&
            (TMath::Abs(MCSigmaStar->GetPdgCode()) != kSigmaStarNCode))
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 3:  // Bkg type C: Xi(1530)0
        // Xi- check
        MCXi = (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiCode)
          return kFALSE;
        // Xi-pion mother check
        if (MCPi->GetMother() != MCXi->GetMother())
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCXi->GetMother()));
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 4:  // Bkg type C: Xi(1530)- through Xi0
        // Xi0 check
        MCXi = (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiZeroCode)
          return kFALSE;
        // Xi-pion mother check
        if (MCPi->GetMother() != MCXi->GetMother())
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCXi->GetMother()));
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarMCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
      case 5:  // Bkg type C: Xi(1530)- through Xi-
               // Not sure that if it is needed since it will make normal Xi-
               // peak.
        // Xi- check
        MCXi = (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCLam->GetMother()));
        if (TMath::Abs(MCXi->GetPdgCode()) != kXiCode)
          return kFALSE;
        // Labmda-pion mother check to Xi-
        if (MCLam->GetMother() != MCXi->GetMother())
          return kFALSE;
        // Xi(1530)0 check
        MCXiStar =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(MCXi->GetMother()));
        if (TMath::Abs(MCXiStar->GetPdgCode()) != kXiStarMCode)
          return kFALSE;
        if (fIsPrimaryMC) {
          if (MCSigmaStar->IsPrimary())
            return kTRUE;
          else
            return kFALSE;
        } else
          return kTRUE;
        break;
        return kFALSE;
    }
  }
}

THnSparse* AliAnalysisTaskSigma1385PM::CreateTHnSparse(TString name,
                                                       TString title,
                                                       Int_t ndim,
                                                       std::vector<TAxis> bins,
                                                       Option_t* opt) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  const TAxis* axises[bins.size()];
  for (UInt_t i = 0; i < bins.size(); i++)
    axises[i] = &bins[i];
  THnSparse* h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
  return h;
}
Long64_t AliAnalysisTaskSigma1385PM::FillTHnSparse(TString name,
                                                   std::vector<Double_t> x,
                                                   Double_t w) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  auto hsparse = dynamic_cast<THnSparse*>(fHistos->FindObject(name));
  if (!hsparse) {
    std::cout << "ERROR : no " << name << std::endl;
    exit(1);
  }
  return FillTHnSparse(hsparse, x, w);
}

Long64_t AliAnalysisTaskSigma1385PM::FillTHnSparse(THnSparse* h,
                                                   std::vector<Double_t> x,
                                                   Double_t w) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  if (int(x.size()) != h->GetNdimensions()) {
    std::cout << "ERROR : wrong sized of array while Fill " << h->GetName()
              << std::endl;
    exit(1);
  }
  return h->Fill(&x.front(), w);
}
TAxis AliAnalysisTaskSigma1385PM::AxisFix(TString name,
                                          int nbin,
                                          Double_t xmin,
                                          Double_t xmax) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(nbin, xmin, xmax);
  axis.SetName(name);
  return axis;
}
TAxis AliAnalysisTaskSigma1385PM::AxisStr(TString name,
                                          std::vector<TString> bin) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
  UInt_t i = 1;
  for (auto blabel : bin)
    ax.SetBinLabel(i++, blabel);
  return ax;
}

TAxis AliAnalysisTaskSigma1385PM::AxisVar(TString name,
                                          std::vector<Double_t> bin) {
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(bin.size() - 1, &bin.front());
  axis.SetName(name);
  return axis;
}
double AliAnalysisTaskSigma1385PM::GetTPCnSigma(AliVTrack* track,
                                                AliPID::EParticleType type) {
  AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
  if (nanoT) {
    static bool used = false;
    if (!used) {
      AliNanoAODTrack::InitPIDIndex();
      used = true;
    }
    return nanoT->GetVar(
        AliNanoAODTrack::GetPIDIndex(AliNanoAODTrack::kSigmaTPC, type));
  } else
    return fPIDResponse->NumberOfSigmasTPC(track, type);
}
void AliAnalysisTaskSigma1385PM::FillTrackToEventPool() {
  // Fill Selected tracks to event mixing pool
  if ((fCentBin < 0) || (fZbin < 0))
    return;
  AliVTrack* goodtrack;

  tracklist* etl;
  eventpool* ep;
  // Event mixing pool

  ep = &fEMpool[fCentBin][fZbin];
  ep->push_back(tracklist());
  etl = &(ep->back());
  // Fill selected tracks
  for (UInt_t i = 0; i < fGoodTrackArray.size(); i++) {
    goodtrack = (AliVTrack*)fEvt->GetTrack(fGoodTrackArray[i]);
    if (!goodtrack)
      continue;
    etl->push_back((AliVTrack*)goodtrack->Clone());
  }
  if (!fGoodTrackArray.size())
    ep->pop_back();
  if ((int)ep->size() > (int)fnMix) {
    for (auto it : ep->front())
      delete it;
    ep->pop_front();
  }
}
void AliAnalysisTaskSigma1385PM::GetImpactParam(AliVTrack* track,
                                                Float_t p[2],
                                                Float_t cov[3]) {
  AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
  if (nanoT)
    nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
  else
    track->GetImpactParameters(p, cov);
}