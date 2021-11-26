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

////////////////////////////////////////////////////////////////////////////
//
//  This class is used to reconstruct the neutral Xi(1530) resonance for
//  Run2 data.
//  This class essentially combines charged Xi candidates from the Xi Vert-
//  exer with primary charged pions.
//
//  author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//        , Beomkyu  KIM (kimb@cern.ch)
//
//  Last Modified Date: 2019/10/05
//
////////////////////////////////////////////////////////////////////////////

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODcascade.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliESDcascade.h"
#include "AliEventCuts.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliVertex.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TParticlePDG.h"
#include "TSystem.h"

// from header
#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "THistManager.h"
// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAnalysisTaskXi1530.h"

// Some constants
const Double_t pi = TMath::Pi();
const Double_t pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t Ximass = 1.32171;
const Double_t massXi1530 = 1.532;
enum {
    kData = 1,
    kLS,
    kMixing,
    kMCReco,
    kMCTrue,
    kMCTruePS,  // 6
    kINEL10,
    kINELg010,
    kAllType
};                                                  // for Physicsl Results
enum { kIsSelected = 1, kPS, kIsMulti, kAllNone };  // for V0M signal QA plot
enum {
    kTrueINELg0 = 1,
    kTrig_TrueINELg0,
    kGoodVtx_TrueINELg0,
    kVzCutted_TrueINELg0,
    kSelected_TrueINELg0,
    kTrue,
    kTrig,
    kGoodVtx,
    kVzCutted,
    kSelected,
};  // for Trigger Efficiency, Vertext Correction
enum {
    kDefaultOption = 1,
    kTPCNsigmaXi1530PionLoose,
    kTPCNsigmaXi1530PionTight,
    kTPCNsigmaLambdaProtonLoose,
    kTPCNsigmaLambdaProtonTight,
    kTPCNsigmaLambdaPionLoose,
    kTPCNsigmaLambdaPionTight,
    kTPCNsigmaBachelorPionLoose,
    kTPCNsigmaBachelorPionTight,
    kXi1530PionZVertexLoose,
    kXi1530PionZVertexTight,
    kDCADistLambdaDaughtersLoose,
    kDCADistLambdaDaughtersTight,
    kDCADistLambdaPVLoose,
    kDCADistLambdaPVTight,
    kV0CosineOfPointingAngleLoose,
    kV0CosineOfPointingAngleTight,
    kCascadeCosineOfPointingAngleLoose,
    kCascadeCosineOfPointingAngleTight,
    kXiMassWindowLoose,
    kXiMassWindowTight,
    kXiTrackCut
};  // for Systematic study
enum {
    kALL = 1,
    kTrigger,
    kInCompleteDAQ,
    kNoBG,
    kNoPileUp,
    kGoodVertex,
    kZvtx10,
    kIENLgtZERO,
    kAliMultSelection
};  // for Event Check
//___________________________________________________________________
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530()
    : AliAnalysisTaskSE("AliAnalysisTaskXi1530"),
      fOption(),
      goodtrackindices(),
      goodcascadeindices(),
      fEMpool() {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//___________________________________________________________________
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530(const char* name,
                                             const char* option)
    : AliAnalysisTaskSE(name),
      fOption(option),
      goodtrackindices(),
      goodcascadeindices(),
      fEMpool() {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
AliAnalysisTaskXi1530::AliAnalysisTaskXi1530(const AliAnalysisTaskXi1530& ap)
    : fOption(ap.fOption),
      goodtrackindices(ap.goodtrackindices),
      goodcascadeindices(ap.goodcascadeindices),
      fEMpool(ap.fEMpool) {}
//___________________________________________________________________
AliAnalysisTaskXi1530& AliAnalysisTaskXi1530::operator=(
    const AliAnalysisTaskXi1530& ap) {
    // assignment operator

    this->~AliAnalysisTaskXi1530();
    new (this) AliAnalysisTaskXi1530(ap);
    return *this;
}
//___________________________________________________________________
AliAnalysisTaskXi1530::~AliAnalysisTaskXi1530() {
    delete fTrackCuts;
    delete fPIDResponse;
}
//________________________________________________________________________
void AliAnalysisTaskXi1530::UserCreateOutputObjects() {
    // TrackCuts for Xi1530--------------------------------------------------
    // Primary pion cut(Xi1530pion)
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    //fTrackCuts->SetPtRange(0.15, 1e20);
    //fTrackCuts->SetEtaRange(-fXi1530PionEtaCut, fXi1530PionEtaCut);
    //fTrackCuts->SetMaxDCAToVertexZ(fXi1530PionZVertexCut_loose);

    // ----------------------------------------------------------------------

    fHistos = new THistManager("Xi1530hists");

    auto binType = AxisStr("Type", {"DATA", "LS", "Mixing", "MCReco", "MCTrue",
                                    "kMCTruePS", "INEL10", "INELg010"});
    std::vector<double> centaxisbin;
    if (!IsMC) {
        if (IsAA && !IsHighMult)
            centaxisbin = {-1, 0,  10, 20, 30, 40, 50,
                           60, 70, 80, 90, 100};  // for AA study
        else if (!IsHighMult)
            centaxisbin = {-1, 0,  1,  5,  10, 15, 20,
                           30, 40, 50, 60, 70, 80, 90, 100};  // for kINT7 study
        else
            centaxisbin = {0, 0.01, 0.03, 0.05, 0.07, 0.1};  // for HM study
    } else
        centaxisbin = {
            -1, 0,  0.01, 0.03, 0.05, 0.07, 0.1, 1, 5,
            10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};   // for kINT7 study

    binCent = AxisVar("Cent", centaxisbin);  // for kINT7 study
    auto binPt = AxisFix("Pt", 200, 0, 20);
    auto binMass = AxisFix("Mass", 1800, 1.2, 3.0);
    binZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});
    auto binType_V0M =
        AxisStr("Type", {"isSelected", "isSelectedPS", "isSelectedMult"});
    // Axis for the systematic study
    SysCheck = {"DefaultOption",
                "TPCNsigmaXi1530PionLoose",
                "TPCNsigmaXi1530PionTight",
                "TPCNsigmaXiLoose",
                "TPCNsigmaXiTight",
                "Xi1530PionZVertexLoose",
                "Xi1530PionZVertexTight",
                "DCADistLambdaDaughtersLoose",
                "DCADistLambdaDaughtersTight",
                "DCADistXiDaughtersLoose",
                "DCADistXiDaughtersTight",
                "DCADistLambdaPVLoose",
                "DCADistLambdaPVTight",
                "V0CosineOfPointingAngleLoose",
                "V0CosineOfPointingAngleTight",
                "CascadeCosineOfPointingAngleLoose",
                "CascadeCosineOfPointingAngleTight",
                "XiMassWindowLoose",
                "XiMassWindowTight",
                "XiTrackCut"};
    binSystematics = AxisStr("Sys", SysCheck);

    CreateTHnSparse("hInvMass_dXi", "InvMass", 4,
                    {binType, binCent, binPt, binMass},
                    "s");  // inv mass distribution of Xi
    CreateTHnSparse("hInvMass", "InvMass", 5,
                    {binSystematics, binType, binCent, binPt, binMass},
                    "s");  // Normal inv mass distribution of Xi1530
    if (IsQAEvent)
        CreateTHnSparse(
            "EventQA/hV0MSignal", "V0MSignal", 4,
            {binType_V0M, binCent, AxisFix("V0MSig", 2500, 0, 25000),
             AxisFix("SPDNtrk", 4000, 0, 4000)},
            "s");
    if (fExoticFinder2)
        CreateTHnSparse("hInvMass_hf", "InvMass", 4,
                        {binType, binCent, binPt, binMass},
                        "s");
    auto binTrklet =
        AxisVar("nTrklet", {0, 5, 10, 15, 20, 25, 30, 35, 40, 100});

    if (IsMC) {
        // To get Trigger efficiency in each trk/V0M Multiplicity region
        auto MCType = AxisStr(
            "Type", {"kTrueINELg0", "kTrig_TrueINELg0", "kGoodVtx_TrueINELg0",
                     "kVzCutted_TrueINELg0", "kSelected_TrueINELg0", "kTrue",
                     "kTrig", "kGoodVtx", "kVzCutted", "kSelected"});
        CreateTHnSparse("htriggered_CINT7", "", 3, {MCType, binCent, binTrklet},
                        "s");
    }
    fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
    if (!IsHighMult)
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7);
    else
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0);
    if (IsMC)
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7);
    // QA Histograms--------------------------------------------------
    //
    if (fQA) {
        if (IsQAEvent) {
            if (IsHighMult) {
                fHistos->CreateTH1("EventQA/hMult_QA", "", 100, 0, 0.1, "s");
                fHistos->CreateTH1("EventQA/hMult_QA_onlyMult", "", 100, 0, 0.1, "s");
                fHistos->CreateTH1("EventQA/hMult_SkippedDataQA", "", 100, 0, 0.1, "s");
                fHistos->CreateTH1("EventQA/hMult_ProcessedDataQA", "", 100, 0, 0.1, "s");
            } else {
                fHistos->CreateTH1("EventQA/hMult_QA", "", 1010, -1, 100, "s");
                fHistos->CreateTH1("EventQA/hMult_QA_onlyMult", "", 101, -1, 100, "s");
                fHistos->CreateTH1("EventQA/hMult_SkippedDataQA", "", 101, -1, 100, "s");
                fHistos->CreateTH1("EventQA/hMult_ProcessedDataQA", "", 101, -1, 100, "s");
            }
        }
        fHistos->CreateTH2("hPhiEta", "", 180, 0, 2 * pi, 40, -2, 2);
        if (IsQAPID) {
            // T P C   P I D
            //// before
            // dEdX
            fHistos->CreateTH2("hTPCPIDLambdaProton", "", 200, 0, 20, 2000, 0, 200);
            fHistos->CreateTH2("hTPCPIDLambdaPion", "", 200, 0, 20, 2000, 0, 200);
            fHistos->CreateTH2("hTPCPIDBachelorPion", "", 200, 0, 20, 2000, 0, 200);
            fHistos->CreateTH2("hTPCPIDXi1530Pion", "", 200, 0, 20, 2000, 0, 200);
            // Signal
            fHistos->CreateTH1("hTPCPIDsignalLambdaProton", "", 100, -5, 5, "s");
            fHistos->CreateTH1("hTPCPIDsignalLambdaPion", "", 100, -5, 5, "s");
            fHistos->CreateTH1("hTPCPIDsignalBachelorPion", "", 100, -5, 5, "s");
            fHistos->CreateTH1("hTPCPIDsignalXi1530Pion", "", 100, -5, 5, "s");
            //
            //// after
            // dEdX
            fHistos->CreateTH2("hTPCPIDLambdaProton_cut", "", 200, 0, 20, 2000, 0,
                            200);
            fHistos->CreateTH2("hTPCPIDLambdaPion_cut", "", 200, 0, 20, 2000, 0,
                            200);
            fHistos->CreateTH2("hTPCPIDBachelorPion_cut", "", 200, 0, 20, 2000, 0,
                            200);
            fHistos->CreateTH2("hTPCPIDXi1530Pion_cut", "", 200, 0, 20, 2000, 0,
                            200);
            // Signal
            fHistos->CreateTH1("hTPCPIDsignalLambdaProton_cut", "", 100, -5, 5,
                            "s");
            fHistos->CreateTH1("hTPCPIDsignalLambdaPion_cut", "", 100, -5, 5, "s");
            fHistos->CreateTH1("hTPCPIDsignalBachelorPion_cut", "", 100, -5, 5,
                            "s");
            fHistos->CreateTH1("hTPCPIDsignalXi1530Pion_cut", "", 100, -5, 5, "s");
            //// Systemtatics
            //
            fHistos->CreateTH1("hTPCPIDsignalLambdaProton_loose", "", 100, -5,
                               5, "s");
            fHistos->CreateTH1("hTPCPIDsignalLambdaPion_loose", "", 100, -5, 5,
                               "s");
            fHistos->CreateTH1("hTPCPIDsignalBachelorPion_loose", "", 100, -5,
                               5, "s");
            fHistos->CreateTH1("hTPCPIDsignalXi1530Pion_loose", "", 100, -5, 5,
                               "s");

            fHistos->CreateTH1("hTPCPIDsignalLambdaProton_tight", "", 100, -5,
                               5, "s");
            fHistos->CreateTH1("hTPCPIDsignalLambdaPion_tight", "", 100, -5, 5,
                               "s");
            fHistos->CreateTH1("hTPCPIDsignalBachelorPion_tight", "", 100, -5,
                               5, "s");
            fHistos->CreateTH1("hTPCPIDsignalXi1530Pion_tight", "", 100, -5, 5,
                               "s");
        }
        

        // D C A
        // between daughters
        // before
        fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters", "", 300, 0, 3, "s");
        fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters", "", 300, 0, 3, "s");
        // after
        fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_cut", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_cut", "", 300, 0, 3, "s");
        // systematics
        fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_loose", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_loose", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_tight", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_tight", "", 300, 0, 3,
                           "s");
        // to PV
        // before
        fHistos->CreateTH1("hDCADist_Xi1530pion_to_PV", "", 300, 0, 3, "s");
        fHistos->CreateTH1("hDCADist_lambda_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_Xi_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_LambdaProton_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_LambdaPion_to_PV", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_BachelorPion_to_PV", "", 500, 0, 5, "s");
        // after
        fHistos->CreateTH1("hDCADist_Xi1530pion_to_PV_cut", "", 300, 0, 3, "s");
        fHistos->CreateTH1("hDCADist_lambda_to_PV_cut", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_Xi_to_PV_cut", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_LambdaProton_to_PV_cut", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("hDCADist_LambdaPion_to_PV_cut", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_BachelorPion_to_PV_cut", "", 500, 0, 5,
                           "s");

        fHistos->CreateTH1("hDCADist_Xi1530pion_to_PV_loose", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_Xi1530pion_to_PV_tight", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_lambda_to_PV_loose", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_lambda_to_PV_tight", "", 500, 0, 5, "s");
        // C P A
        // before
        fHistos->CreateTH1("hCosPA_lambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("hCosPA_Xi", "", 150, 0.85, 1.0, "s");
        // after
        fHistos->CreateTH1("hCosPA_lambda_cut", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("hCosPA_Xi_cut", "", 150, 0.85, 1.0, "s");

        fHistos->CreateTH1("hCosPA_lambda_loose", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("hCosPA_Xi_loose", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("hCosPA_lambda_tight", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("hCosPA_Xi_tight", "", 150, 0.85, 1.0, "s");

        // M a s s   W i n d o w
        fHistos->CreateTH1("hMass_Xi", "", 200, 1.2, 1.4, "s");      // before
        fHistos->CreateTH1("hMass_Xi_cut", "", 200, 1.2, 1.4, "s");  // after

        fHistos->CreateTH1("hMass_Xi_loose", "", 200, 1.2, 1.4, "s");
        fHistos->CreateTH1("hMass_Xi_tight", "", 200, 1.2, 1.4, "s");

        // E t a
        fHistos->CreateTH2("hPhiEta_Xi", "", 180, 0, 2 * pi, 40, -2,
                           2);  // before
        fHistos->CreateTH2("hPhiEta_Xi_cut", "", 180, 0, 2 * pi, 40, -2,
                           2);  // after

        // Radius X - Y
        fHistos->CreateTH2("hLambda_Rxy", "", 400, -200, 200, 400, -200,
                           200);  // before
        fHistos->CreateTH2("hLambda_Rxy_cut", "", 400, -200, 200, 400, -200,
                           200);  // after
        fHistos->CreateTH2("hXi_Rxy", "", 400, -200, 200, 400, -200,
                           200);  // before
        fHistos->CreateTH2("hXi_Rxy_cut", "", 400, -200, 200, 400, -200,
                           200);  // after
        if (fExoticFinder)
            fHistos->CreateTH1("hExoOpenAngle", "", 180, 0, 2 * pi, "s");
        if (IsMC) {
            // For MC True stduy purpose
            fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_TrueMC", "", 300,
                               0, 3, "s");
            fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_TrueMC", "", 300, 0,
                               3, "s");
            fHistos->CreateTH1("hDCADist_LambdaProton_to_PV_TrueMC", "", 500, 0,
                               5, "s");
            fHistos->CreateTH1("hDCADist_LambdaPion_to_PV_TrueMC", "", 500, 0,
                               5, "s");
            fHistos->CreateTH1("hDCADist_BachelorPion_to_PV_TrueMC", "", 500, 0,
                               5, "s");
            fHistos->CreateTH1("hDCADist_lambda_to_PV_TrueMC", "", 500, 0, 5,
                               "s");
            fHistos->CreateTH1("hDCADist_Xi_to_PV_TrueMC", "", 500, 0, 5, "s");
            fHistos->CreateTH2("hPhiEta_Xi_TrueMC", "", 180, 0, 2 * pi, 40, -2,
                               2);
            fHistos->CreateTH2("hLambda_Rxy_TrueMC", "", 400, -200, 200, 400,
                               -200, 200);
            fHistos->CreateTH2("hXi_Rxy_TrueMC", "", 400, -200, 200, 400, -200,
                               200);
            fHistos->CreateTH1("hMC_generated_Y", "", 400, -2, 2, "s");
            fHistos->CreateTH1("hMC_reconstructed_Y", "", 400, -2, 2, "s");
        }
    }
    if (IsQAInvMass) {
        // Invmass Check
        fHistos->CreateTH1("hTotalInvMass_data", "", 1300, 1.2, 2.5, "s");
        fHistos->CreateTH1("hTotalInvMass_LS", "", 1300, 1.2, 2.5, "s");
        fHistos->CreateTH1("hTotalInvMass_Mix", "", 1300, 1.2, 2.5, "s");
        if (fExoticFinder2){
            fHistos->CreateTH1("hTotalInvMass_HFpp", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFnp", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFpn", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFnn", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFppMix", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFnpMix", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFpnMix", "", 2000, 1.5, 3.5, "s");
            fHistos->CreateTH1("hTotalInvMass_HFnnMix", "", 2000, 1.5, 3.5, "s");
        }
    }
    fEMpool.resize(binCent.GetNbins() + 1,
                   std::vector<eventpool>(binZ.GetNbins() + 1));
    PostData(1, fHistos->GetListOfHistograms());
}

//________________________________________________________________________
void AliAnalysisTaskXi1530::UserExec(Option_t*) {
    // Pointer to a event----------------------------------------------------
    AliVEvent* event = InputEvent();
    if (!event) {
        AliInfo("Could not retrieve event");
        return;
    }
    // ----------------------------------------------------------------------

    // NanoAOD --------------------------------------------------------------
    AliNanoAODHeader* nanoHeader =
        dynamic_cast<AliNanoAODHeader*>(event->GetHeader());
    if ((!IsNano) && (nanoHeader))
        IsNano = true;
    // ----------------------------------------------------------------------

    // Connect to ESD tree --------------------------------------------------
    event->IsA() == AliESDEvent::Class()
        ? fEvt = dynamic_cast<AliESDEvent*>(event)
        : fEvt = dynamic_cast<AliAODEvent*>(event);
    if (!fEvt)
        return;
    if (fEvt->IsA() != AliESDEvent::Class())
        IsAOD = kTRUE;
    // ----------------------------------------------------------------------

    // Load InputHandler for each event--------------------------------------
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();
    // ----------------------------------------------------------------------

    // Global Variables -----------------------------------------------------
    // Vertex Check
    const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
    // const AliVVertex* trackVtx  = fEvt->GetPrimaryVertexTracks() ;
    lPosPV[0] = pVtx->GetX();
    lPosPV[1] = pVtx->GetY();
    lPosPV[2] = pVtx->GetZ();

    // Initialize
    bool IsMultSelcted{false}, IsSelectedTrig{false}, IsEvtSelected{false};
    Double_t intensity = 0.;
    bField = fEvt->GetMagneticField();  // bField for track DCA

    if (!IsNano) {
        IsEvtSelected = fEventCuts.AcceptEvent(event);
        // Preparation for MC
        // ---------------------------------------------------
        if (IsMC) {
            if (IsAOD)
                fMCArray = (TClonesArray*)fEvt->FindListObject(
                    "mcparticles");  // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = fEventCuts.IsTrueINELgtZero(fEvt, true);
        }
        // Trigger Check
        // ----------------------------------------------------------
        IsSelectedTrig = fEventCuts.PassedCut(AliEventCuts::kTrigger);

        // Multiplicity(centrality)
        // ---------------------------------------------
        if (!IsAOD) {
            ftrackmult = fEvt->GetMultiplicity()->GetNumberOfTracklets();
        } else {
            ftrackmult =
                ((AliAODEvent*)fEvt)->GetTracklets()->GetNumberOfTracklets();
        }

        // fCent = GetMultiplicty(fEvt);  // Centrality(AA), Multiplicity(pp)
        fCent = AliMultSelectionTask::IsINELgtZERO(event) 
                    ? fEventCuts.GetCentrality() 
                    : -0.5;

        // PID response
        // ----------------------------------------------------------
        fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
        if (!fPIDResponse)
            AliInfo("No PIDd");

        // Vertex Cuts
        Bool_t IsGoodVertex =
            fEventCuts.PassedCut(AliEventCuts::kVertexQuality);
        Bool_t IsVtxInZCut =
            fEventCuts.PassedCut(AliEventCuts::kVertexPosition);

        Bool_t IsINELg0 = fEventCuts.PassedCut(AliEventCuts::kINELgt0);
        AliMultSelection* MultSelection =
            (AliMultSelection*)fEvt->FindListObject("MultSelection");
        IsMultSelcted = MultSelection->IsEventSelected();
        // Pile up rejection
        // -----------------------------------------------------
        Bool_t IsNotPileUp = fEventCuts.PassedCut(AliEventCuts::kPileUp);

        // In Complete DAQ Event Cut
        // --------------------------------------------------
        Bool_t IncompleteDAQ =
            fEventCuts.PassedCut(AliEventCuts::kDAQincomplete);

        IsPS = IsSelectedTrig     // CINT7 Trigger selected
               && IncompleteDAQ  // No IncompleteDAQ
               && IsNotPileUp;    // PileUp rejection

        IsINEL0Rec = IsPS && IsGoodVertex && IsVtxInZCut && IsINELg0;

        if (IsSelectedTrig && IsMultSelcted && fQA && IsQAEvent)
            fHistos->FillTH1("EventQA/hMult_QA_onlyMult", (double)fCent);

        //  Missing Vetex and Trriger Efficiency
        //  ---------------------------------
        if (IsMC) {
            // bin
            if (IsINEL0True) {
                FillTHnSparse("htriggered_CINT7",
                              {(double)kTrueINELg0, (double)fCent, ftrackmult});
                if (IsSelectedTrig)
                    FillTHnSparse(
                        "htriggered_CINT7",
                        {(double)kTrig_TrueINELg0, (double)fCent, ftrackmult});
                if (IsPS && IsGoodVertex)
                    FillTHnSparse("htriggered_CINT7",
                                  {(double)kGoodVtx_TrueINELg0, (double)fCent,
                                   ftrackmult});
                if (IsPS && IsGoodVertex && IsVtxInZCut)
                    FillTHnSparse("htriggered_CINT7",
                                  {(double)kVzCutted_TrueINELg0, (double)fCent,
                                   ftrackmult});
                if (IsINEL0Rec && IsMultSelcted)
                    FillTHnSparse("htriggered_CINT7",
                                  {(double)kSelected_TrueINELg0, (double)fCent,
                                   ftrackmult});
            }
            FillTHnSparse("htriggered_CINT7",
                          {(double)kTrue, (double)fCent, ftrackmult});
            if (IsSelectedTrig)
                FillTHnSparse("htriggered_CINT7",
                              {(double)kTrig, (double)fCent, ftrackmult});
            if (IsPS && IsGoodVertex)
                FillTHnSparse("htriggered_CINT7",
                              {(double)kGoodVtx, (double)fCent, ftrackmult});
            if (IsPS && IsGoodVertex && IsVtxInZCut)
                FillTHnSparse("htriggered_CINT7",
                              {(double)kVzCutted, (double)fCent, ftrackmult});
            if (IsINEL0Rec && IsMultSelcted)
                FillTHnSparse("htriggered_CINT7",
                              {(double)kSelected, (double)fCent, ftrackmult});

            // Signal Loss Correction
            if (!IsAOD) {
                if (IsINEL0True && IsVtxInZCut) {  // INEL>0|10
                    FillMCinput(fMCEvent, 1);
                    FillMCinputdXi(fMCEvent, 1);
                }
                if (IsVtxInZCut) {  // INEL10
                    FillMCinput(fMCEvent, 2);
                    FillMCinputdXi(fMCEvent, 2);
                }
                if (IsSelectedTrig) {  // INEL>0 +|Vz| < 10 cm
                    FillMCinput(fMCEvent, 3);
                    FillMCinputdXi(fMCEvent, 3);
                }
            } else {
                if (IsINEL0True && IsVtxInZCut) {  // INEL>0|10
                    FillMCinputAOD(fMCEvent, 1);
                    FillMCinputdXiAOD(fMCEvent, 1);
                }
                if (IsVtxInZCut) {  // INEL10
                    FillMCinputAOD(fMCEvent, 2);
                    FillMCinputdXiAOD(fMCEvent, 2);
                }
                if (IsSelectedTrig) {  // INEL>0 +|Vz| < 10 cm
                    FillMCinputAOD(fMCEvent, 3);
                    FillMCinputdXiAOD(fMCEvent, 3);
                }
            }
        }
        // V0M Signal QA
        // ---------------------------------------------------------
        // From BeomKyu Kim
        AliVVZERO* lV0 = fEvt->GetVZEROData();
        if (inputHandler->IsEventSelected()) {
            for (int i = 0; i < 64; i++)
                intensity += lV0->GetMultiplicity(i);
        }

        // ----------------------------------------------------------------------
    } else {
        // Preparation for MC
        // ---------------------------------------------------
        if (IsMC) {
            if (IsAOD)
                fMCArray = (TClonesArray*)fEvt->FindListObject(
                    "mcparticles");  // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = true;
        }
        fCent = nanoHeader->GetCentr("V0M");
        static int inel_index = -1;
        if (inel_index < 0) inel_index = nanoHeader->GetVarIndex("cstINELgt0");
        if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
            fCent = -0.5;
        static int v0mValueIndex =
            nanoHeader->GetVarIndex("MultSelection.V0M.Value");
        static int trkValueIndex =
            nanoHeader->GetVarIndex("MultSelection.SPDTracklets.Value");
        intensity = nanoHeader->GetVar(v0mValueIndex);
        ftrackmult = nanoHeader->GetVar(trkValueIndex);

        if (IsQAEvent)
            fHistos->FillTH1("EventQA/hMult_QA_onlyMult", (double)fCent);

        IsSelectedTrig = true;
        IsMultSelcted = true;
        IsEvtSelected = true;
    }

    if (IsSelectedTrig && IsQAEvent) {
        FillTHnSparse("EventQA/hV0MSignal",
                      {kIsSelected, (double)fCent, (double)intensity,
                       (double)ftrackmult});
        if (IsMultSelcted)
            FillTHnSparse("EventQA/hV0MSignal",
                          {kIsMulti, (double)fCent, (double)intensity,
                           (double)ftrackmult});
    }

    // Event Mixing pool -----------------------------------------------------
    zbin = binZ.FindBin(lPosPV[2]) - 1;           // Event mixing z-bin
    centbin = binCent.FindBin(fCent) - 1;  // Event mixing cent bin
    if (isINEL)
        centbin = 0;                       // for INEL case
    // -----------------------------------------------------------------------

    // Check tracks and casade, Fill histo************************************
    if (IsEvtSelected) {  // In AliEventCuts
        // Draw Multiplicity QA plot in only selected event.
        if (fQA && IsQAEvent) {
            fHistos->FillTH1("EventQA/hMult_QA", (double)fCent);

            // V0M signal QA
            FillTHnSparse(
                "EventQA/hV0MSignal",
                {kPS, (double)fCent, (double)intensity, (double)ftrackmult});
        }
        if (IsMC) {  // After All Event cut!
            if (!IsAOD) {
                FillMCinput(fMCEvent, 4);
                FillMCinputdXi(fMCEvent, 4);
            } else {
                FillMCinputAOD(fMCEvent, 4);
                FillMCinputdXiAOD(fMCEvent, 4);
            }
        }
        bool checkPionTrack = this->GoodTracksSelection();
        bool checkCascade = this->GoodCascadeSelection();
        if (checkPionTrack      // If Good track
            && checkCascade) {  // and Good cascade is in
                                // this event,
            if (!IsAOD)
                this->FillTracks();  // Fill the histogram
            else
                this->FillTracksAOD();  // Fill the histogram(AOD)
        }
    }
    // ***********************************************************************

    PostData(1, fHistos->GetListOfHistograms());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskXi1530::GoodTracksSelection() {
    // Choose Good Tracks from AliESDtracks, AliAODTracks
    //    and Save the label of them,
    //    and Save them for event mixing
    //    for the systematic study, this will be done in "loose cut option"
    //
    // it includes several cuts:
    // - TPC PID cut for pion
    // - Eta cut
    // - Z-vertex cut
    // - DCA cut for primary particle
    // - pion mass window cut <-- not using now
    //
    const UInt_t ntracks = fEvt->GetNumberOfTracks();
    goodtrackindices.clear();
    goodtrackfullindices.clear();
    AliVTrack* track;
    tracklist* etl;
    eventpool* ep;
    // Event mixing pool
    if ( (centbin >= 0) && (zbin >= 0) && fsetmixing) {
        ep = &fEMpool[centbin][zbin];
        ep->push_back(tracklist());
        etl = &(ep->back());
    }

    Float_t b[2];
    Float_t bCov[3];
    Double_t pionZ = -999;
    Double_t pionPt, fEta, fTPCNSigPion;

    for (UInt_t it = 0; it < ntracks; it++) {
        track = (AliVTrack*)fEvt->GetTrack(it);
        if (!track){
            AliInfo("no track");
            continue;
        }
        GetImpactParam(track, b, bCov);

        // Track cuts
        if (!IsAOD) {  // ESD Case
            if (!fTrackCuts->AcceptTrack((AliESDtrack*)track)) {
                AliInfo("can't pass the track cut");
                continue;
            }
        }
        else {
            if (!IsNano) {
                if (!((AliAODTrack*)track)->TestFilterBit(fFilterBit))
                    continue;
            } else {
                if (!(static_cast<AliNanoAODTrack*>(track)->TestFilterBit(
                        fFilterBit)))
                    continue;
            }
        }

        pionPt = track->Pt();
        pionZ = b[1];
        fEta = TMath::Abs(track->Eta());
        fTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);

        if (fEta > fXi1530PionEtaCut) {
            AliInfo(Form("Eta cut failed: track eta: %f, cut: %f",
                            track->Eta(), fXi1530PionEtaCut));
            continue;
        }
        if (pionPt < 0.15) {
            AliInfo(Form("pT cut failed: track pT: %f, cut: 0.15", pionPt));
            continue;
        }
        if (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut_loose) {
            AliInfo(Form("TPC nSigma cut failed: track TPC PID: %f, cut: %f",
                         fTPCNSigPion, fTPCNsigXi1530PionCut_loose));
            continue;
        }

        goodtrackfullindices.push_back(it); // General pion track
        
        if (pionZ > fXi1530PionZVertexCut_loose) {
            AliInfo(Form("zVertex cut failed: track DCAz: %f, cut: %f",
                            pionZ, fXi1530PionZVertexCut_loose));
            continue;
        }
        
        // if (TMath::Abs(track->M() - pionmass) > 0.007) continue;
        if (fQA) {
            fHistos->FillTH2("hPhiEta", track->Phi(), track->Eta());
            if(IsQAPID){ 
                fHistos->FillTH2("hTPCPIDXi1530Pion", track->GetTPCmomentum(),
                             track->GetTPCsignal());
                fHistos->FillTH1("hTPCPIDsignalXi1530Pion", fTPCNSigPion);
            }
            fHistos->FillTH1("hDCADist_Xi1530pion_to_PV", pionZ);
        }  // After default cut

        goodtrackindices.push_back(it);
        // Event mixing pool
        if ((centbin >= 0) && (zbin >= 0) && fsetmixing) {
            AliVTrack* track_mix = (AliVTrack*)fEvt->GetTrack(it);
            etl->push_back((AliVTrack*)track_mix->Clone());
        }
    }

    if ((centbin >= 0) && (zbin >= 0) && fsetmixing) {
        if (!goodtrackindices.size())
            ep->pop_back();
            Int_t epsize = ep->size();
            if (epsize > fnMix) {
                for (auto it : ep->front())
                    delete it;
                ep->pop_front();
        }
    }

    return goodtrackindices.size();
}

Bool_t AliAnalysisTaskXi1530::GoodCascadeSelection() {
    // Choose Good Cascade from AliESDcascade, AliAODcascade
    //   and Save the label of them
    //    for the systematic study, this will be done in "loose cut option"
    //
    // it includes several cuts:
    // - daughter particle standard track cut
    // - daughter particle PID cut
    // - DCA cuts between Lambda daughters and Xi daughters
    // - PV DCA(Impact parameter) cut for Xi/Lambda/all daughters (partially
    // applied)
    // - Cosine Pointing Angle cut for Xi and Lamnbda
    // - Mass window cut for Xi
    // - Eta cut
    // - XY Raidus cut(not applied), just for check

    goodcascadeindices.clear();
    const UInt_t ncascade = fEvt->GetNumberOfCascades();

    const AliESDcascade* Xicandidate;
    const AliAODcascade* Xicandidate_aod;
    Double_t LambdaX, LambdaY, LambdaZ;
    Double_t fTPCNSigProton, fTPCNSigLambdaPion, fTPCNSigBachelorPion;
    Double_t fDCADist_LambdaProton_PV, fDCADist_LambdaPion_PV;
    Double_t fDCADist_Lambda_PV, fDCADist_Xi_PV, fDCADist_BachelorPion_PV;
    Double_t fDCADist_Lambda, fDCADist_Xi;
    Double_t cX, cY, cZ;
    Double_t fLambdaCPA, fXiCPA;
    Double_t fMass_Xi;
    Double_t lPosXi[3];
    Double_t lMomXi[3];
    Double_t Xi_momsum, Xi_eta, Xi_phi;
    Float_t b[2];
    Float_t bCov[3];

    Bool_t StandardXi = kTRUE;
    if (!IsAOD) {  // ESD case
        for (UInt_t it = 0; it < ncascade; it++) {
            StandardXi = kTRUE;
            Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(it);
            if (!Xicandidate) {
                AliInfo("No Cascade!");
                continue;
            }

            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetNindex())) {
                AliInfo("Same index p-n");
                continue;
            }
            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetBindex())) {
                AliInfo("Same index p-b");
                continue;
            }
            if (TMath::Abs(Xicandidate->GetNindex()) ==
                TMath::Abs(Xicandidate->GetBindex())) {
                AliInfo("Same index n-b");
                continue;
            }

            AliESDtrack* pTrackXi =
                (AliESDtrack*)fEvt->GetTrack(TMath::Abs(Xicandidate->GetPindex()));
            AliESDtrack* nTrackXi =
                (AliESDtrack*)fEvt->GetTrack(TMath::Abs(Xicandidate->GetNindex()));
            AliESDtrack* bTrackXi =
                (AliESDtrack*)fEvt->GetTrack(TMath::Abs(Xicandidate->GetBindex()));

            // PID cuts for Xi daughters
            if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
                if (fQA && IsQAPID) {
                    fHistos->FillTH2("hTPCPIDLambdaProton",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                }
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
                if (fQA && IsQAPID) {
                    fHistos->FillTH2("hTPCPIDLambdaProton",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                }
            }
            fTPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion
            if (fQA && IsQAPID) {
                fHistos->FillTH2("hTPCPIDBachelorPion",
                                 bTrackXi->GetTPCmomentum(),
                                 bTrackXi->GetTPCsignal());
                fHistos->FillTH1("hTPCPIDsignalLambdaProton", fTPCNSigProton);
                fHistos->FillTH1("hTPCPIDsignalLambdaPion", fTPCNSigLambdaPion);
                fHistos->FillTH1("hTPCPIDsignalBachelorPion",
                                 fTPCNSigBachelorPion);
            }
            if (TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut_loose){
                AliInfo(Form("proton PID cut failed -  value: %f, cut: %f",
                             fTPCNSigProton, fTPCNsigLambdaProtonCut_loose));
                StandardXi = kFALSE;  // PID for proton
            }
            if (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut_loose) {
                AliInfo(Form("pion PID cut failed -  value: %f, cut: %f",
                             fTPCNSigLambdaPion, fTPCNsigLambdaPionCut_loose));
                StandardXi = kFALSE;  // PID for 1st pion
            }
            if (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut_loose) {
                AliInfo(Form("bPion PID cut failed -  value: %f, cut: %f",
                             fTPCNSigBachelorPion,
                             fTPCNsigBachelorPionCut_loose));
                StandardXi = kFALSE;  // PID for 2nd pion
            }

            // DCA cut
            // DCA between Dautgher particles
            fDCADist_Lambda = TMath::Abs(Xicandidate->GetDcaV0Daughters());
            fDCADist_Xi = TMath::Abs(Xicandidate->GetDcaXiDaughters());
            if (fQA) {
                fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters",
                                 fDCADist_Lambda);
                fHistos->FillTH1("hDCADist_Xi_BTW_Daughters", fDCADist_Xi);
            }

            if (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_loose) {
                AliInfo(Form(
                    "DCA Lambda daughters cut failed -  value: %f, cut: %f",
                    fDCADist_Lambda, fDCADist_LambdaDaughtersCut_loose));
                StandardXi = kFALSE;  // DCA proton-pion
            }
            if (fDCADist_Xi > fDCADist_XiDaughtersCut_loose) {
                AliInfo(
                    Form("DCA Xi daughters cut failed -  value: %f, cut: %f",
                         fDCADist_Xi, fDCADist_XiDaughtersCut_loose));
                StandardXi = kFALSE;  // DCA Lambda-pion
            }

            // DCA to PV
            fDCADist_Lambda_PV =
                TMath::Abs(Xicandidate->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
            fDCADist_Xi_PV =
                TMath::Abs(Xicandidate->GetDcascade(lPosPV[0], lPosPV[1], lPosPV[2]));
            if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                fDCADist_LambdaProton_PV = b[0];
                GetImpactParam(nTrackXi, b, bCov);
                fDCADist_LambdaPion_PV = b[0];
            } else {
                GetImpactParam(nTrackXi, b, bCov);
                fDCADist_LambdaProton_PV = b[0];
                GetImpactParam(pTrackXi, b, bCov);
                fDCADist_LambdaPion_PV = b[0];
            }
            GetImpactParam(bTrackXi, b, bCov);
            fDCADist_BachelorPion_PV = b[0];
            if (fQA) {
                fHistos->FillTH1("hDCADist_lambda_to_PV", fDCADist_Lambda_PV);
                fHistos->FillTH1("hDCADist_Xi_to_PV", fDCADist_Xi_PV);
                fHistos->FillTH1("hDCADist_LambdaProton_to_PV",
                                 fDCADist_LambdaProton_PV);
                fHistos->FillTH1("hDCADist_LambdaPion_to_PV",
                                 fDCADist_LambdaPion_PV);
                fHistos->FillTH1("hDCADist_BachelorPion_to_PV",
                                 fDCADist_BachelorPion_PV);
            }

            if (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_loose) {
                AliInfo(Form("DCA Lambda pV cut failed -  value: %f, cut: %f",
                             fDCADist_Lambda_PV, fDCADist_Lambda_PVCut_loose));
                StandardXi = kFALSE;  // DCA Lambda-vertex
            }

            // CPA cut
            Xicandidate->GetXYZcascade(cX, cY, cZ);
            if (!fCPAstudy)
                fLambdaCPA =
                    Xicandidate->GetV0CosineOfPointingAngle(cX, cY, cZ);
            else
                fLambdaCPA =
                    Xicandidate->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);
            fXiCPA =
                Xicandidate->GetCascadeCosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);
            if (fQA) {
                fHistos->FillTH1("hCosPA_lambda", fLambdaCPA);
                fHistos->FillTH1("hCosPA_Xi", fXiCPA);
            }

            if (fLambdaCPA < fV0CosineOfPointingAngleCut_loose) {
                AliInfo(Form("CPA Lambda cut failed -  value: %f, cut: %f",
                             fLambdaCPA, fV0CosineOfPointingAngleCut_loose));
                StandardXi = kFALSE;
            }
            if (fXiCPA < fCascadeCosineOfPointingAngleCut_loose) {
                AliInfo(Form("CPA Xi cut failed -  value: %f, cut: %f", fXiCPA,
                             fCascadeCosineOfPointingAngleCut_loose));
                StandardXi = kFALSE;
            }

            // Mass window cut
            fMass_Xi = Xicandidate->GetEffMassXi();
            if (fQA)
                fHistos->FillTH1("hMass_Xi", fMass_Xi);
            /*
            // Disable for Xi study
            if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_loose)
                StandardXi = kFALSE;
            */
            // Eta cut
            if (TMath::Abs(Xicandidate->Eta()) > fXiEtaCut) {
                AliInfo(Form("Eta cut failed -  value: %f, cut: %f",
                             Xicandidate->Eta(), fXiEtaCut));
                StandardXi = kFALSE;
            }
            if (fQA)
                fHistos->FillTH2("hPhiEta_Xi", Xicandidate->Phi(),
                                 Xicandidate->Eta());

            // XY Raidus cut(experiemntal)
            Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);
            if (fQA)
                fHistos->FillTH2("hLambda_Rxy", LambdaX, LambdaY);
            // if(sqrt( pow(LambdaX,2) + pow(LambdaY,2) ) > 100)
            // StandardXi=kFALSE; // NOT USING

            if (fQA)
                fHistos->FillTH2("hXi_Rxy", cX, cY);
            if ((sqrt(pow(cX, 2) + pow(cY, 2)) > 12) && fExoticFinder) {
                AliInfo(Form("Exotic xy cut failed -  value: %f, cut: 12",
                             sqrt(pow(cX, 2) + pow(cY, 2))));
                StandardXi = kFALSE;  // NOT USING in normal mode
            }
            

            // After selection above
            if (StandardXi) {  // Save only the Xi is good candidate
                if ((Xicandidate->RapXi() > fXi1530RapidityCut_high) ||
                    (Xicandidate->RapXi() < fXi1530RapidityCut_low)) {
                    FillTHnSparse("hInvMass_dXi",
                                  {(int)kData, (double)fCent, Xicandidate->Pt(),
                                   Xicandidate->GetEffMassXi()});
                    if (IsMC)
                        if (IsTrueXi(it)) {
                            FillTHnSparse(
                                "hInvMass_dXi",
                                {(int)kMCReco, (double)fCent, Xicandidate->Pt(),
                                 Xicandidate->GetEffMassXi()});
                        }
                }
                goodcascadeindices.push_back(it);
            }  // for standard Xi
        }      // All Xi loop
    }          // ESD case
    else {
        for (UInt_t it = 0; it < ncascade; it++) {
            StandardXi = kTRUE;
            Xicandidate_aod = ((AliAODEvent*)fEvt)->GetCascade(it);
            if (!Xicandidate_aod) {
                AliInfo("No Cascade!");
                continue;
            }

            if (TMath::Abs(Xicandidate_aod->GetPosID()) ==
                TMath::Abs(Xicandidate_aod->GetNegID())) {
                AliInfo("Same index p-n");
                continue;
            }
            if (TMath::Abs(Xicandidate_aod->GetPosID()) ==
                TMath::Abs(Xicandidate_aod->GetBachID())) {
                AliInfo("Same index p-b");
                continue;
            }
            if (TMath::Abs(Xicandidate_aod->GetNegID()) ==
                TMath::Abs(Xicandidate_aod->GetBachID())) {
                AliInfo("Same index n-b");
                continue;
            }

            AliAODTrack* pTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDaughter(0));
            AliAODTrack* nTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDaughter(1));
            AliAODTrack* bTrackXi =
                (AliAODTrack*)(Xicandidate_aod->GetDecayVertexXi()->GetDaughter(
                    0));

            // PID cuts for Xi daughters
            if (Xicandidate_aod->ChargeXi() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
                if (fQA && IsQAPID) {
                    fHistos->FillTH2("hTPCPIDLambdaProton",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                }
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
                if (fQA && IsQAPID) {
                    fHistos->FillTH2("hTPCPIDLambdaProton",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                }
            }
            fTPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion
            if (fQA && IsQAPID) {
                fHistos->FillTH2("hTPCPIDBachelorPion",
                                 bTrackXi->GetTPCmomentum(),
                                 bTrackXi->GetTPCsignal());
                fHistos->FillTH1("hTPCPIDsignalLambdaProton", fTPCNSigProton);
                fHistos->FillTH1("hTPCPIDsignalLambdaPion", fTPCNSigLambdaPion);
                fHistos->FillTH1("hTPCPIDsignalBachelorPion",
                                 fTPCNSigBachelorPion);
            }
            if (TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut_loose) {
                AliInfo(Form("proton PID cut failed -  value: %f, cut: %f",
                             fTPCNSigProton, fTPCNsigLambdaProtonCut_loose));
                StandardXi = kFALSE;  // PID for proton
            }
            if (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut_loose) {
                AliInfo(Form("pion PID cut failed -  value: %f, cut: %f",
                             fTPCNSigLambdaPion, fTPCNsigLambdaPionCut_loose));
                StandardXi = kFALSE;  // PID for 1st pion
            }
            if (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut_loose) {
                AliInfo(Form("bPion PID cut failed -  value: %f, cut: %f",
                             fTPCNSigBachelorPion,
                             fTPCNsigBachelorPionCut_loose));
                StandardXi = kFALSE;  // PID for 2nd pion
            }

            // DCA cut
            // DCA between Dautgher particles
            fDCADist_Lambda = TMath::Abs(Xicandidate_aod->DcaV0Daughters());
            fDCADist_Xi = TMath::Abs(Xicandidate_aod->DcaXiDaughters());
            if (fQA) {
                fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters",
                                 fDCADist_Lambda);
                fHistos->FillTH1("hDCADist_Xi_BTW_Daughters", fDCADist_Xi);
            }

            if (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_loose) {
                AliInfo(Form(
                    "DCA Lambda daughters cut failed -  value: %f, cut: %f",
                    fDCADist_Lambda, fDCADist_LambdaDaughtersCut_loose));
                StandardXi = kFALSE;  // DCA proton-pion
            }
            if (fDCADist_Xi > fDCADist_XiDaughtersCut_loose) {
                AliInfo(
                    Form("DCA Xi daughters cut failed -  value: %f, cut: %f",
                         fDCADist_Xi, fDCADist_XiDaughtersCut_loose));
                StandardXi = kFALSE;  // DCA Lambda-pion
            }

            // DCA to PV
            fDCADist_Lambda_PV =
                TMath::Abs(Xicandidate_aod->DcaV0ToPrimVertex());
            fDCADist_Xi_PV =
                TMath::Abs(Xicandidate_aod->DcaXiToPrimVertex(lPosPV[0], lPosPV[1], lPosPV[2]));
            if (Xicandidate_aod->ChargeXi() == -1) {  // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                fDCADist_LambdaProton_PV = b[0];
                GetImpactParam(nTrackXi, b, bCov);
                fDCADist_LambdaPion_PV = b[0];
                /*
                fDCADist_LambdaProton_PV =
                    TMath::Abs(Xicandidate_aod->DcaPosToPrimVertex());
                fDCADist_LambdaPion_PV =
                    TMath::Abs(Xicandidate_aod->DcaNegToPrimVertex());
                */
            } else {
                GetImpactParam(nTrackXi, b, bCov);
                fDCADist_LambdaProton_PV = b[0];
                GetImpactParam(pTrackXi, b, bCov);
                fDCADist_LambdaPion_PV = b[0];
                /*
                fDCADist_LambdaProton_PV =
                    TMath::Abs(Xicandidate_aod->DcaNegToPrimVertex());
                fDCADist_LambdaPion_PV =
                    TMath::Abs(Xicandidate_aod->DcaPosToPrimVertex());
                */
            }
            GetImpactParam(bTrackXi, b, bCov);
            fDCADist_BachelorPion_PV = b[0];
            /*
            Double_t fDCADist_BachelorPion_PV =
                TMath::Abs(Xicandidate_aod->DcaBachToPrimVertex());
            */
            if (fQA) {
                fHistos->FillTH1("hDCADist_lambda_to_PV", fDCADist_Lambda_PV);
                fHistos->FillTH1("hDCADist_Xi_to_PV", fDCADist_Xi_PV);
                fHistos->FillTH1("hDCADist_LambdaProton_to_PV",
                                 fDCADist_LambdaProton_PV);
                fHistos->FillTH1("hDCADist_LambdaPion_to_PV",
                                 fDCADist_LambdaPion_PV);
                fHistos->FillTH1("hDCADist_BachelorPion_to_PV",
                                 fDCADist_BachelorPion_PV);
            }

            if (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_loose) {
                AliInfo(Form("DCA Lambda pV cut failed -  value: %f, cut: %f",
                             fDCADist_Lambda_PV, fDCADist_Lambda_PVCut_loose));
                StandardXi = kFALSE;  // DCA Lambda-vertex
            }

            // CPA cut
            lPosXi[0] = Xicandidate_aod->DecayVertexXiX();
            lPosXi[1] = Xicandidate_aod->DecayVertexXiY();
            lPosXi[2] = Xicandidate_aod->DecayVertexXiZ();
            if (!fCPAstudy)
                fLambdaCPA = Xicandidate_aod->CosPointingAngle(lPosXi);
            else
                fLambdaCPA = Xicandidate_aod->CosPointingAngle(lPosPV);
            fXiCPA =
                Xicandidate_aod->CosPointingAngleXi(lPosPV[0], lPosPV[1], lPosPV[2]);
            if (fQA) {
                fHistos->FillTH1("hCosPA_lambda", fLambdaCPA);
                fHistos->FillTH1("hCosPA_Xi", fXiCPA);
            }

            if (fLambdaCPA < fV0CosineOfPointingAngleCut_loose) {
                AliInfo(Form("CPA Lambda cut failed -  value: %f, cut: %f",
                             fLambdaCPA, fV0CosineOfPointingAngleCut_loose));
                StandardXi = kFALSE;
            }
            if (fXiCPA < fCascadeCosineOfPointingAngleCut_loose) {
                AliInfo(Form("CPA Xi cut failed -  value: %f, cut: %f", fXiCPA,
                             fCascadeCosineOfPointingAngleCut_loose));
                StandardXi = kFALSE;
            }

            // Mass window cut
            fMass_Xi = Xicandidate_aod->MassXi();
            if (fQA)
                fHistos->FillTH1("hMass_Xi", fMass_Xi);
            /*
            // Disable for Xi study
            if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_loose)
                StandardXi = kFALSE;
            */

            // Eta cut
            // Eta: 0.5*TMath::Log((P()+Pz())/(P()-Pz()+1.e-13))
            // Phi: TMath::Pi()+TMath::ATan2(-Py(),-Px())
            lMomXi[0] = Xicandidate_aod->MomXiX();
            lMomXi[1] = Xicandidate_aod->MomXiY();
            lMomXi[2] = Xicandidate_aod->MomXiZ();
            Xi_momsum =
                TMath::Sqrt(lMomXi[0] * lMomXi[0] + lMomXi[1] * lMomXi[1] +
                            lMomXi[2] * lMomXi[2]);
            Xi_eta =
                0.5 * TMath::Log((Xi_momsum + lMomXi[2]) /
                                 (Xi_momsum - lMomXi[2] + 1.e-13));
            Xi_phi =
                TMath::Pi() + TMath::ATan2(-lMomXi[1], -lMomXi[0]);
            if (TMath::Abs(Xi_eta) > fXiEtaCut) {
                AliInfo(Form("Eta cut failed -  value: %f, cut: %f", Xi_eta,
                             fXiEtaCut));
                StandardXi = kFALSE;
            }
            if (fQA)
                fHistos->FillTH2("hPhiEta_Xi", Xi_phi, Xi_eta);

            // XY Raidus cut(experiemntal)
            // Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);
            if (fQA)
                fHistos->FillTH2("hLambda_Rxy",
                                 Xicandidate_aod->DecayVertexV0X(),
                                 Xicandidate_aod->DecayVertexV0Y());
            // if(sqrt( Xicandidate->RadiusV0() ) > 100)
            // StandardXi=kFALSE; // NOT USING

            if (fQA)
                fHistos->FillTH2("hXi_Rxy", lPosXi[0], lPosXi[1]);
            if ((sqrt(TMath::Sqrt(TMath::Power(lPosXi[0], 2) +
                                  TMath::Power(lPosXi[1], 2))) > 12) &&
                fExoticFinder) {
                AliInfo("Exotic xy cut failed");
                StandardXi = kFALSE;  // NOT USING in normal mode
            }
            // After selection above
            if (StandardXi) {  // Save only the Xi is good candidate
                if ((Xicandidate_aod->RapXi() > fXi1530RapidityCut_high) ||
                    (Xicandidate_aod->RapXi() < fXi1530RapidityCut_low)) {
                    FillTHnSparse("hInvMass_dXi", {(int)kData, (double)fCent,
                                                Xicandidate_aod->Pt2Xi(),
                                                Xicandidate_aod->MassXi()});
                    if (IsMC)
                        if (IsTrueXi(it)) {
                            FillTHnSparse("hInvMass_dXi",
                                        {(int)kMCReco, (double)fCent,
                                        Xicandidate_aod->Pt2Xi(),
                                        Xicandidate_aod->MassXi()});
                        }
                }
                goodcascadeindices.push_back(it);
            }  // for standard Xi
        }      // All Xi loop
    }          // AOD case

    return goodcascadeindices.size();
}

void AliAnalysisTaskXi1530::FillTracks() {
    AliVTrack* track1;           // charged track, pion
    AliVTrack* track2;           // charged track, pion
    AliESDcascade* Xicandidate;  // Cascade
    tracklist trackpool;

    TLorentzVector temp1, temp2, temp3;
    TLorentzVector vecsum;  // Xi1530 candidate
    TLorentzVector vecsum2;  // for
    Double_t LambdaX, LambdaY, LambdaZ;
    Double_t fTPCNSigProton, fTPCNSigLambdaPion, fTPCNSigBachelorPion,
        fDCADist_Xi_PV, fDCADist_Lambda_PV;

    // for DCA value
    Float_t b[2];
    Float_t bCov[3];
    Double_t pionZ = -999;

    const UInt_t ncascade = goodcascadeindices.size();
    const UInt_t ntracks = goodtrackindices.size();
    const UInt_t ntracks_full = goodtrackfullindices.size();

    for (UInt_t sys = 0; sys < (UInt_t)binSystematics.GetNbins(); sys++) {
        // Systematic study loop.
        // sys = 0 -> Default cut option
        // for more details, please check "SysCheck" in header file.
        AliInfo(Form("Sys check! %s", (const char*)SysCheck.at(sys)));
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate) {
                AliInfo(Form("No Xi! %s", (const char*)SysCheck.at(sys)));
                continue;
            }

            AliESDtrack* pTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetPindex()));
            AliESDtrack* nTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetNindex()));
            AliESDtrack* bTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetBindex()));

            if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }
            fTPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion

            temp1.SetXYZM(Xicandidate->Px(), Xicandidate->Py(),
                          Xicandidate->Pz(), Xicandidate->GetEffMassXi());

            for (UInt_t j = 0; j < ntracks; j++) {
                track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[j]);
                if (!track1) {
                    AliInfo(
                        Form("No track! %s", (const char*)SysCheck.at(sys)));
                    continue;
                }

                if (track1->GetID() == pTrackXi->GetID() ||
                    track1->GetID() == nTrackXi->GetID() ||
                    track1->GetID() == bTrackXi->GetID()) {
                    AliInfo(
                        Form("same track! %s", (const char*)SysCheck.at(sys)));
                    continue;
                }

                // PID Cut Systematic check
                // -------------------------------------------------
                Double_t fTPCNSigPion = GetTPCnSigma(track1, AliPID::kPion);

                // Xi1530Pion PID
                if ((SysCheck.at(sys) != "TPCNsigmaXi1530PionLoose") &&
                    (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut)) {
                    AliInfo(Form("pion PID! %f %s", fTPCNSigPion,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }

                if ((SysCheck.at(sys) == "TPCNsigmaXi1530PionTight") &&
                    (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut_tight))
                    continue;
                // Xi PID
                if (SysCheck.at(sys) != "TPCNsigmaXiLoose") {
                    if ((TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut) ||
                        (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut) ||
                        (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut)) {
                        AliInfo(Form("Xi PID! %f %f %f %s", fTPCNSigProton,
                                     fTPCNSigLambdaPion, fTPCNSigBachelorPion,
                                     (const char*)SysCheck.at(sys)));
                        continue;
                    }
                }
                if (SysCheck.at(sys) == "TPCNsigmaXiTight") {
                    if ((TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut_tight) ||
                        (TMath::Abs(fTPCNSigLambdaPion) >
                         fTPCNsigLambdaPionCut_tight) ||
                        (TMath::Abs(fTPCNSigBachelorPion) >
                         fTPCNsigBachelorPionCut_tight))
                        continue;
                }

                // Xi1530Pion DCA zVetex Check
                GetImpactParam(track1, b, bCov);
                pionZ = b[1];
                if ((SysCheck.at(sys) != "Xi1530PionZVertexLoose") &&
                    (pionZ > fXi1530PionZVertexCut)) {
                    AliInfo(Form("pionZ! %f %s", pionZ,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "Xi1530PionZVertexTight") &&
                    (pionZ > fXi1530PionZVertexCut_tight))
                    continue;

                // DCA between daughters Check
                Double_t fDCADist_Lambda = TMath::Abs(Xicandidate->GetDcaV0Daughters());
                Double_t fDCADist_Xi = TMath::Abs(Xicandidate->GetDcaXiDaughters());
                if ((SysCheck.at(sys) != "DCADistLambdaDaughtersLoose") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut)) {
                    AliInfo(Form("DCADistLambdaDaughters! %f %s",
                                 fDCADist_Lambda,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }

                if ((SysCheck.at(sys) == "DCADistLambdaDaughtersTight") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "DCADistXiDaughtersLoose") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut)) {
                    AliInfo(Form("DCADistXiDaughters! %f %s", fDCADist_Xi,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }

                if ((SysCheck.at(sys) == "DCADistXiDaughtersTight") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut_tight))
                    continue;

                // DCA Lambda to PV Check
                fDCADist_Lambda_PV =
                    TMath::Abs(Xicandidate->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
                fDCADist_Xi_PV =
                    TMath::Abs(Xicandidate->GetDcascade(lPosPV[0], lPosPV[1], lPosPV[2]));
                if ((SysCheck.at(sys) != "DCADistLambdaPVLoose") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut)) {
                    AliInfo(Form("DCADistLambdaPV! %f %s", fDCADist_Lambda_PV,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }

                if ((SysCheck.at(sys) == "DCADistLambdaPVTight") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_tight))
                    continue;

                // CPA Check
                Double_t cX, cY, cZ;
                Xicandidate->GetXYZcascade(cX, cY, cZ);
                Double_t fLambdaCPA;
                if (!fCPAstudy)
                    fLambdaCPA =
                        Xicandidate->GetV0CosineOfPointingAngle(cX, cY, cZ);
                else
                    fLambdaCPA =
                        Xicandidate->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);
                Double_t fXiCPA =
                    Xicandidate->GetCascadeCosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);

                if ((SysCheck.at(sys) != "V0CosineOfPointingAngleLoose") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut)) {
                    AliInfo(Form("V0CosineOfPointingAngle! %f %s", fLambdaCPA,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "V0CosineOfPointingAngleTight") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "CascadeCosineOfPointingAngleLoose") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut)) {
                    AliInfo(Form("CascadeCosineOfPointingAngle! %f %s", fXiCPA,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut_tight))
                    continue;

                // Xi Mass Window Check
                Double_t fMass_Xi = Xicandidate->GetEffMassXi();
                if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_loose) {
                    continue;
                }
                if ((SysCheck.at(sys) != "XiMassWindowLoose") &&
                    (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut)) {
                    AliInfo(Form("XiMassWindow! %f %s", fMass_Xi,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "XiMassWindowTight") &&
                    (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_tight))
                    continue;

                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionmass);

                vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
                // Y cut
                if ((vecsum.Rapidity() > fXi1530RapidityCut_high) ||
                    (vecsum.Rapidity() < fXi1530RapidityCut_low))
                    continue;

                // Opening Angle - Not using in normal mode
                if (fExoticFinder) {
                    Double_t angle = temp1.Angle(temp2.Vect());
                    fHistos->FillTH1("hExoOpenAngle", angle);
                    if (TMath::Abs(angle) < 0.0785398)  // 4.5 degree
                        continue;
                }

                auto sign = kAllType;
                if ((Xicandidate->Charge() == -1 && track1->Charge() == +1) ||
                    (Xicandidate->Charge() == +1 && track1->Charge() == -1))
                    sign = kData;  // Unlike sign -> Data
                else
                    sign = kLS;  // like sign bg

                if (IsMC) {
                    if (IsTrueXi1530(Xicandidate,
                                     track1)) {  // MC Association, if it
                                                 // comes from True Xi1530

                        // True Xi1530 signals
                        FillTHnSparse("hInvMass",
                                      {(double)sys, (double)kMCReco,
                                       (double)fCent, vecsum.Pt(), vecsum.M()});
                        Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);

                        if (fQA) {
                            fHistos->FillTH1("hMC_reconstructed_Y",
                                             vecsum.Rapidity());
                            // For cut study
                            fHistos->FillTH1(
                                "hDCADist_Lambda_BTW_Daughters_TrueMC",
                                TMath::Abs(Xicandidate->GetDcaV0Daughters()));
                            fHistos->FillTH1(
                                "hDCADist_Xi_BTW_Daughters_TrueMC",
                                TMath::Abs(Xicandidate->GetDcaXiDaughters()));
                            if (Xicandidate->Charge() ==
                                -1) {  // Xi- has +proton, -pion
                                GetImpactParam(pTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC", b[0]);
                                GetImpactParam(nTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC", b[0]);
                            } else {
                                GetImpactParam(nTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC", b[0]);
                                GetImpactParam(pTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC", b[0]);
                            }
                            GetImpactParam(bTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_BachelorPion_to_PV_TrueMC", b[0]);

                            fHistos->FillTH1(
                                "hDCADist_lambda_to_PV_TrueMC",
                                TMath::Abs(Xicandidate->GetD(lPosPV[0], lPosPV[1], lPosPV[2])));
                            fHistos->FillTH1(
                                "hDCADist_Xi_to_PV_TrueMC",
                                TMath::Abs(Xicandidate->GetDcascade(lPosPV[0], lPosPV[1], lPosPV[2])));

                            fHistos->FillTH2("hPhiEta_Xi_TrueMC",
                                             Xicandidate->Phi(),
                                             Xicandidate->Eta());
                            fHistos->FillTH2("hLambda_Rxy_TrueMC", LambdaX,
                                             LambdaY);

                            fHistos->FillTH2("hXi_Rxy_TrueMC", cX, cY);
                        }

                    }  // Xi1530 check
                }      // MC
                FillTHnSparse("hInvMass",
                              {(double)sys, (double)sign, (double)fCent,
                               vecsum.Pt(), vecsum.M()});
                if (IsQAInvMass && sys == 0) {
                    if ((int)sign == (int)kData)
                        fHistos->FillTH1("hTotalInvMass_data", vecsum.M());
                    if ((int)sign == (int)kLS)
                        fHistos->FillTH1("hTotalInvMass_LS", vecsum.M());
                }
                if (fExoticFinder2) {
                    if (SysCheck.at(sys) == "DefaultOption"){
                        for (UInt_t k = 0; k < ntracks_full; k++) {
                            track2 =
                                (AliVTrack*)fEvt->GetTrack(goodtrackfullindices[k]);
                            if (!track2)
                                continue;
                            if (track1->GetID() == track2->GetID())
                                continue;
                            temp3.SetXYZM(track2->Px(), track2->Py(), track2->Pz(),
                                        pionmass);

                            vecsum2 =
                                vecsum + temp3;  // vecsum = Xi1530, temp3=pion
                            // Y cut
                            if ((vecsum2.Rapidity() > fXi1530RapidityCut_high) ||
                                (vecsum2.Rapidity() < fXi1530RapidityCut_low))
                                continue;
                            
                            // Mass window
                            double mXi1530 = vecsum.M();
                            if (TMath::Abs(mXi1530 - massXi1530) > 0.05)
                                continue;

                            int sign2 = kData;
                            if (track2->Charge() > 0){
                                if (track1->Charge() > 0)
                                    sign2 = kData;
                                else
                                    sign2 = kLS;
                            }
                            else {
                                if (track1->Charge() > 0)
                                    sign2 = kMixing;
                                else
                                    sign2 = kMCReco;
                            }
                            FillTHnSparse("hInvMass_hf",
                                        {(double)sign2, (double)fCent,
                                        vecsum2.Pt(), vecsum2.M()});
                            if (IsQAInvMass) {
                                if (track2->Charge() > 0) {
                                    if (track1->Charge() > 0)
                                        fHistos->FillTH1("hTotalInvMass_HFpp",
                                                         vecsum2.M());
                                    else
                                        fHistos->FillTH1("hTotalInvMass_HFnp",
                                                         vecsum2.M());
                                } else {
                                    if (track1->Charge() > 0)
                                        fHistos->FillTH1("hTotalInvMass_HFpn",
                                                         vecsum2.M());
                                    else
                                        fHistos->FillTH1("hTotalInvMass_HFnn",
                                                         vecsum2.M());
                                }
                            }
                        }
                    }
                }

                // Fill the QA Histos
                if (fQA) {
                    if (SysCheck.at(sys) == "DefaultOption") {
                        if (IsQAPID)
                            fHistos->FillTH2("hTPCPIDXi1530Pion_cut",
                                             track1->GetTPCmomentum(),
                                             track1->GetTPCsignal());
                        if (Xicandidate->Charge() ==
                            -1) {  // Xi- has +proton, -pion
                            if (IsQAPID) {
                                fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                                 pTrackXi->GetTPCmomentum(),
                                                 pTrackXi->GetTPCsignal());
                                fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                                 nTrackXi->GetTPCmomentum(),
                                                 nTrackXi->GetTPCsignal());
                            }
                            GetImpactParam(pTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_LambdaProton_to_PV_cut", b[0]);
                            GetImpactParam(nTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_LambdaPion_to_PV_cut", b[0]);
                        } else {  // Xi+ has -proton, +pion
                            if (IsQAPID) {
                                fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                                 nTrackXi->GetTPCmomentum(),
                                                 nTrackXi->GetTPCsignal());
                                fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                                 pTrackXi->GetTPCmomentum(),
                                                 pTrackXi->GetTPCsignal());
                            }
                            GetImpactParam(nTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_LambdaProton_to_PV_cut", b[0]);
                            GetImpactParam(pTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_LambdaPion_to_PV_cut", b[0]);
                        }
                        if (IsQAPID) {
                            fHistos->FillTH2("hTPCPIDBachelorPion_cut",
                                             bTrackXi->GetTPCmomentum(),
                                             bTrackXi->GetTPCsignal());

                            // TPC PID Signal
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_cut",
                                             fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_cut",
                                             fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_cut",
                                             fTPCNSigBachelorPion);
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_cut",
                                             fTPCNSigPion);
                        }
                        // DCA QA
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_cut",
                                         fDCADist_Lambda);
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_cut",
                                         fDCADist_Xi);
                        fHistos->FillTH1("hDCADist_lambda_to_PV_cut",
                                         fDCADist_Lambda_PV);
                        fHistos->FillTH1(
                            "hDCADist_Xi_to_PV_cut",
                            TMath::Abs(Xicandidate->GetDcascade(lPosPV[0], lPosPV[1], lPosPV[2])));
                        GetImpactParam(bTrackXi, b, bCov);
                        fHistos->FillTH1(
                            "hDCADist_BachelorPion_to_PV_cut", b[0]);
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_cut",
                                         pionZ);
                        // CPA QA
                        fHistos->FillTH1("hCosPA_lambda_cut", fLambdaCPA);
                        fHistos->FillTH1("hCosPA_Xi_cut", fXiCPA);

                        // Mass window QA
                        fHistos->FillTH1("hMass_Xi_cut", fMass_Xi);

                        // Eta
                        fHistos->FillTH2("hPhiEta_Xi_cut", Xicandidate->Phi(),
                                         Xicandidate->Eta());

                        // XY Radius
                        Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);
                        fHistos->FillTH2("hLambda_Rxy_cut", LambdaX, LambdaY);
                        Double_t cX, cY, cZ;
                        Xicandidate->GetXYZcascade(cX, cY, cZ);
                        fHistos->FillTH2("hXi_Rxy_cut", cX, cY);
                    }
                    // PID
                    if (IsQAPID) {
                        if (SysCheck.at(sys) == "TPCNsigmaXi1530PionLoose")
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_loose",
                                            fTPCNSigPion);
                        if (SysCheck.at(sys) == "TPCNsigmaXi1530PionTight")
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_tight",
                                            fTPCNSigPion);
                        if (SysCheck.at(sys) == "TPCNsigmaXiLoose") {
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_loose",
                                            fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_loose",
                                            fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_loose",
                                            fTPCNSigBachelorPion);
                        }
                        if (SysCheck.at(sys) == "TPCNsigmaXiTight") {
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_tight",
                                            fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_tight",
                                            fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_tight",
                                            fTPCNSigBachelorPion);
                        }
                    }
                    // Xi1530Pion DCA zVetex Check
                    if (SysCheck.at(sys) == "Xi1530PionZVertexLoose")
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_loose",
                                         pionZ);
                    if (SysCheck.at(sys) == "Xi1530PionZVertexTight")
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_tight",
                                         pionZ);

                    // DCA between daughters Check
                    if (SysCheck.at(sys) == "DCADistLambdaDaughtersLoose")
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_loose",
                                         fDCADist_Lambda);
                    if (SysCheck.at(sys) == "DCADistLambdaDaughtersTight")
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_tight",
                                         fDCADist_Lambda);
                    if (SysCheck.at(sys) == "DCADistXiDaughtersLoose")
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_loose",
                                         fDCADist_Xi);
                    if (SysCheck.at(sys) == "DCADistXiDaughtersTight")
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_tight",
                                         fDCADist_Xi);
                    // Lambda DCA zVetex Check
                    if (SysCheck.at(sys) == "DCADistLambdaPVLoose")
                        fHistos->FillTH1("hDCADist_lambda_to_PV_loose",
                                         fDCADist_Lambda_PV);
                    if (SysCheck.at(sys) == "DCADistLambdaPVTight")
                        fHistos->FillTH1("hDCADist_lambda_to_PV_tight",
                                         fDCADist_Lambda_PV);
                    // CPA Check
                    if (SysCheck.at(sys) == "V0CosineOfPointingAngleLoose")
                        fHistos->FillTH1("hCosPA_lambda_loose", fLambdaCPA);
                    if (SysCheck.at(sys) == "V0CosineOfPointingAngleTight")
                        fHistos->FillTH1("hCosPA_lambda_tight", fLambdaCPA);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleLoose")
                        fHistos->FillTH1("hCosPA_Xi_loose", fXiCPA);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight")
                        fHistos->FillTH1("hCosPA_Xi_tight", fXiCPA);

                    // Xi Mass window
                    if (SysCheck.at(sys) == "XiMassWindowLoose")
                        fHistos->FillTH1("hMass_Xi_loose", fMass_Xi);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight")
                        fHistos->FillTH1("hMass_Xi_tight", fMass_Xi);
                }
            }
        }
        AliInfo(Form("Sys check! %.u", sys));
        if ((!fsetsystematics) && (sys == 0))
            break;
    }

    // Event Mixing
    if ((centbin >= 0) && (zbin >= 0) && fsetmixing) {
        eventpool& ep = fEMpool[centbin][zbin];
        Int_t epsize = ep.size();
        if (IsQAEvent) {
            if (epsize < fnMix) {
                fHistos->FillTH1("EventQA/hMult_SkippedDataQA", (double)fCent);
                return;
            }
            fHistos->FillTH1("EventQA/hMult_ProcessedDataQA", (double)fCent);
        }

        Int_t nForSkipSameEvent = 0;
        for (auto pool : ep) {
            if (nForSkipSameEvent == (epsize - 1))
                continue; // same event
            for (auto track : pool)
                trackpool.push_back((AliVTrack*)track);
            nForSkipSameEvent++;
        }
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate)
                continue;
            temp1.SetXYZM(Xicandidate->Px(), Xicandidate->Py(),
                          Xicandidate->Pz(), Xicandidate->GetEffMassXi());

            AliESDtrack* pTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetPindex()));
            AliESDtrack* nTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetNindex()));
            AliESDtrack* bTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetBindex()));

            for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
                track1 = trackpool.at(jt);
                if (track1->GetID() == pTrackXi->GetID() ||
                    track1->GetID() == nTrackXi->GetID() ||
                    track1->GetID() == bTrackXi->GetID())
                    continue;
                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionmass);
                vecsum = temp1 + temp2;  // two pion vector sum

                if ((Xicandidate->Charge() == -1 && track1->Charge() == -1) ||
                    (Xicandidate->Charge() == +1 && track1->Charge() == +1))
                    continue;  // check only unlike-sign

                if ((vecsum.Rapidity() > fXi1530RapidityCut_high) ||
                    (vecsum.Rapidity() < fXi1530RapidityCut_low))
                    continue;  // rapidity cut

                // Other default cuts
                Double_t fTPCNSigPion = GetTPCnSigma(track1, AliPID::kPion);
                if ((TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut))
                    continue;
                // Xi PID
                if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                    fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                    fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
                } else {  // Xi+ has -proton, +pion
                    fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                    fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
                }
                fTPCNSigBachelorPion = GetTPCnSigma(
                    bTrackXi, AliPID::kPion);  // bachelor is always pion
                if (TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut)
                    continue;
                if (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
                    continue;
                if (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
                    continue;

                // Xi1530Pion DCA zVetex Check
                GetImpactParam(track1, b, bCov);
                pionZ = b[1];
                if (pionZ > fXi1530PionZVertexCut)
                    continue;

                // DCA between daughters Check
                Double_t fDCADist_Lambda =
                    TMath::Abs(Xicandidate->GetDcaV0Daughters());
                Double_t fDCADist_Xi = TMath::Abs(Xicandidate->GetDcaXiDaughters());
                if (fDCADist_Lambda > fDCADist_LambdaDaughtersCut)
                    continue;
                if (fDCADist_Xi > fDCADist_XiDaughtersCut)
                    continue;

                // DCA Lambda to PV Check
                fDCADist_Lambda_PV =
                    TMath::Abs(Xicandidate->GetD(lPosPV[0], lPosPV[1], lPosPV[2]));
                if (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut)
                    continue;

                // CPA Check
                Double_t cX, cY, cZ;
                Xicandidate->GetXYZcascade(cX, cY, cZ);
                Double_t fLambdaCPA;
                if (!fCPAstudy)
                    fLambdaCPA =
                        Xicandidate->GetV0CosineOfPointingAngle(cX, cY, cZ);
                else
                    fLambdaCPA =
                        Xicandidate->GetV0CosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);
                Double_t fXiCPA =
                    Xicandidate->GetCascadeCosineOfPointingAngle(lPosPV[0], lPosPV[1], lPosPV[2]);

                if (fLambdaCPA < fV0CosineOfPointingAngleCut)
                    continue;
                if (fXiCPA < fCascadeCosineOfPointingAngleCut)
                    continue;
                // Xi Mass Window Check
                Double_t fMass_Xi = Xicandidate->GetEffMassXi();
                if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut)
                    continue;

                FillTHnSparse("hInvMass",
                              {(double)kDefaultOption, (double)kMixing,
                               (double)fCent, vecsum.Pt(), vecsum.M()});
                if (IsQAInvMass)
                    fHistos->FillTH1("hTotalInvMass_Mix", vecsum.M());
                if (fExoticFinder2) {
                    for (UInt_t k = 0; k < ntracks_full; k++) {
                        track2 =
                            (AliVTrack*)fEvt->GetTrack(goodtrackfullindices[k]);
                        if (track1->GetID() == track2->GetID())
                            continue;
                        temp3.SetXYZM(track2->Px(), track2->Py(), track2->Pz(),
                                      pionmass);
                        vecsum2 =
                            vecsum + temp3;  // vecsum = Xi1530, temp3=pion
                        // Y cut
                        if ((vecsum2.Rapidity() > fXi1530RapidityCut_high) ||
                            (vecsum2.Rapidity() < fXi1530RapidityCut_low))
                            continue;
                        
                        // Mass window
                        double mXi1530 = vecsum.M();
                        if (TMath::Abs(mXi1530 - massXi1530) > 0.05)
                            continue;

                        int sign2 = kMCTrue;
                        if (track2->Charge() > 0){
                            if (track1->Charge() > 0)
                                sign2 = kMCTrue;
                            else
                                sign2 = kMCTruePS;
                        }
                        else {
                            if (track1->Charge() > 0)
                                sign2 = kINEL10;
                            else
                                sign2 = kINELg010;
                        }
                        FillTHnSparse("hInvMass_hf",
                                      {(double)sign2, (double)fCent,
                                       vecsum2.Pt(), vecsum2.M()});
                        if (IsQAInvMass) {
                            if (track2->Charge() > 0) {
                                if (track1->Charge() > 0)
                                    fHistos->FillTH1("hTotalInvMass_HFppMix",
                                                     vecsum2.M());
                                else
                                    fHistos->FillTH1("hTotalInvMass_HFnpMix",
                                                     vecsum2.M());
                            } else {
                                if (track1->Charge() > 0)
                                    fHistos->FillTH1("hTotalInvMass_HFpnMix",
                                                     vecsum2.M());
                                else
                                    fHistos->FillTH1("hTotalInvMass_HFnnMix",
                                                     vecsum2.M());
                            }
                        }
                    }
                }
            }
        }
    }       // mix loop
}
void AliAnalysisTaskXi1530::FillTracksAOD() {
    AliVTrack* track1;         // charged track, pion
    AliVTrack* track2;         // charged track, pion
    AliAODcascade* Xicandidate;  // Cascade
    tracklist trackpool;

    TLorentzVector temp1, temp2, temp3;
    TLorentzVector vecsum;  // Xi1530 candidate
    TLorentzVector vecsum2; // for
    Double_t fTPCNSigProton, fTPCNSigLambdaPion, fTPCNSigBachelorPion,
        fDCADist_Lambda_PV, fDCADist_Xi_PV;
    Double_t LambdaX, LambdaY, LambdaZ;

    // for DCA value
    Float_t b[2];
    Float_t bCov[3];
    Double_t pionZ = 999;

    const UInt_t ncascade = goodcascadeindices.size();
    const UInt_t ntracks = goodtrackindices.size();
    const UInt_t ntracks_full = goodtrackfullindices.size();

    for (UInt_t sys = 0; sys < (UInt_t)binSystematics.GetNbins(); sys++) {
        // Systematic study loop.
        // sys = 0 -> Default cut option
        // for more details, please check "SysCheck" in header file.
        AliInfo(Form("Sys check! %s", (const char*)SysCheck.at(sys)));
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliAODEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate){
                AliInfo(Form("No Xi! %s", (const char*)SysCheck.at(sys)));
                continue;
            }
            AliAODTrack* pTrackXi = (AliAODTrack*)(Xicandidate->GetDaughter(0));
            AliAODTrack* nTrackXi = (AliAODTrack*)(Xicandidate->GetDaughter(1));
            AliAODTrack* bTrackXi =
                (AliAODTrack*)(Xicandidate->GetDecayVertexXi()->GetDaughter(0));

            if (Xicandidate->ChargeXi() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }
            fTPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion

            temp1.SetXYZM(Xicandidate->MomXiX(), Xicandidate->MomXiY(),
                          Xicandidate->MomXiZ(), Xicandidate->MassXi());

            for (UInt_t j = 0; j < ntracks; j++) {
                track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[j]);
                if (!track1){
                    AliInfo(Form("No track! %s", (const char*)SysCheck.at(sys)));
                    continue;
                }

                if (track1->GetID() == pTrackXi->GetID() ||
                    track1->GetID() == nTrackXi->GetID() ||
                    track1->GetID() == bTrackXi->GetID()){
                        AliInfo(
                            Form("same track! %s", (const char*)SysCheck.at(sys)));
                        continue;
                    }
                    

                // PID Cut Systematic check
                // -------------------------------------------------
                Double_t fTPCNSigPion = GetTPCnSigma(track1, AliPID::kPion);

                // Xi1530Pion PID
                if ((SysCheck.at(sys) != "TPCNsigmaXi1530PionLoose") &&
                    (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut)){
                    AliInfo(Form("pion PID! %f %s", fTPCNSigPion,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }

                if ((SysCheck.at(sys) == "TPCNsigmaXi1530PionTight") &&
                    (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut_tight))
                    continue;
                // Xi PID
                if (SysCheck.at(sys) != "TPCNsigmaXiLoose") {
                    if ((TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut) ||
                        (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut) ||
                        (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut) ){
                        AliInfo(Form("Xi PID! %f %f %f %s", fTPCNSigProton,
                                     fTPCNSigLambdaPion, fTPCNSigBachelorPion,
                                     (const char*)SysCheck.at(sys)));
                        continue;
                    }
                }
                if (SysCheck.at(sys) == "TPCNsigmaXiTight") {
                    if ((TMath::Abs(fTPCNSigProton) > 
                         fTPCNsigLambdaProtonCut_tight) ||
                        (TMath::Abs(fTPCNSigLambdaPion) >
                         fTPCNsigLambdaPionCut_tight) ||
                        (TMath::Abs(fTPCNSigBachelorPion) >
                         fTPCNsigBachelorPionCut_tight))
                            continue;
                }

                // Xi1530Pion DCA zVetex Check
                GetImpactParam(track1, b, bCov);
                pionZ = b[1];
                if ((SysCheck.at(sys) != "Xi1530PionZVertexLoose") &&
                    (pionZ > fXi1530PionZVertexCut)){
                    AliInfo(Form("pionZ! %f %s", pionZ,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "Xi1530PionZVertexTight") &&
                    (pionZ > fXi1530PionZVertexCut_tight))
                    continue;

                // DCA between daughters Check
                Double_t fDCADist_Lambda = TMath::Abs(Xicandidate->DcaV0Daughters());
                Double_t fDCADist_Xi = TMath::Abs(Xicandidate->DcaXiDaughters());
                if ((SysCheck.at(sys) != "DCADistLambdaDaughtersLoose") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut)){
                    AliInfo(Form("DCADistLambdaDaughters! %f %s",
                                 fDCADist_Lambda,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                    
                if ((SysCheck.at(sys) == "DCADistLambdaDaughtersTight") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "DCADistXiDaughtersLoose") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut)){
                    AliInfo(Form("DCADistXiDaughters! %f %s", fDCADist_Xi,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                    
                if ((SysCheck.at(sys) == "DCADistXiDaughtersTight") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut_tight))
                    continue;

                // DCA Lambda to PV Check
                fDCADist_Lambda_PV =
                    TMath::Abs(Xicandidate->DcaV0ToPrimVertex());
                fDCADist_Xi_PV =
                    TMath::Abs(Xicandidate->DcaXiToPrimVertex(lPosPV[0], lPosPV[1], lPosPV[2]));
                if ((SysCheck.at(sys) != "DCADistLambdaPVLoose") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut)){
                    AliInfo(Form("DCADistLambdaPV! %f %s", fDCADist_Lambda_PV,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                    
                if ((SysCheck.at(sys) == "DCADistLambdaPVTight") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_tight))
                    continue;

                // CPA Check
                Double_t lPosXi[3];
                lPosXi[0] = Xicandidate->DecayVertexXiX();
                lPosXi[1] = Xicandidate->DecayVertexXiY();
                lPosXi[2] = Xicandidate->DecayVertexXiZ();
                Double_t fLambdaCPA;
                if (!fCPAstudy)
                    fLambdaCPA = Xicandidate->CosPointingAngle(lPosXi);
                else
                    fLambdaCPA = Xicandidate->CosPointingAngle(lPosPV);
                Double_t fXiCPA =
                    Xicandidate->CosPointingAngleXi(lPosPV[0], lPosPV[1], lPosPV[2]);

                if ((SysCheck.at(sys) != "V0CosineOfPointingAngleLoose") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut)){
                    AliInfo(Form("V0CosineOfPointingAngle! %f %s", fLambdaCPA,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "V0CosineOfPointingAngleTight") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "CascadeCosineOfPointingAngleLoose") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut)){
                    AliInfo(Form("CascadeCosineOfPointingAngle! %f %s", fXiCPA,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                }
                if ((SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut_tight))
                    continue;

                // Xi Mass Window Check
                Double_t fMass_Xi = Xicandidate->MassXi();
                if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_loose) {
                    continue;
                }
                if ((SysCheck.at(sys) != "XiMassWindowLoose") &&
                    (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut)){
                    AliInfo(Form("XiMassWindow! %f %s", fMass_Xi,
                                 (const char*)SysCheck.at(sys)));
                    continue;
                    }
                if ((SysCheck.at(sys) == "XiMassWindowTight") &&
                    (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut_tight))
                    continue;

                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionmass);

                vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
                // Y cut
                if ((vecsum.Rapidity() > fXi1530RapidityCut_high) ||
                    (vecsum.Rapidity() < fXi1530RapidityCut_low))
                    continue;

                // Opening Angle - Not using in normal mode
                if (fExoticFinder) {
                    Double_t angle = temp1.Angle(temp2.Vect());
                    fHistos->FillTH1("hExoOpenAngle", angle);
                    if (TMath::Abs(angle) < 0.0785398)  // 4.5 degree
                        continue;
                }

                auto sign = kAllType;
                if ((Xicandidate->ChargeXi() == -1 && track1->Charge() == +1) ||
                    (Xicandidate->ChargeXi() == +1 && track1->Charge() == -1))
                    sign = kData;  // Unlike sign -> Data
                else
                    sign = kLS;  // like sign bg

                // Phi, Eta
                Double_t xiPx = Xicandidate->Px();
                Double_t xiPy = Xicandidate->Py();
                Double_t xiPz = Xicandidate->Pz();
                Double_t Xi_momsum =
                    TMath::Sqrt(xiPx * xiPx + xiPy * xiPy + xiPz * xiPz);
                Double_t Xi_eta = 0.5 * TMath::Log((Xi_momsum + xiPz) /
                                                   (Xi_momsum - xiPz + 1.e-13));
                Double_t Xi_phi = TMath::Pi() + TMath::ATan2(-xiPy, -xiPx);

                LambdaX = Xicandidate->DecayVertexV0X();
                LambdaY = Xicandidate->DecayVertexV0Y();
                LambdaZ = Xicandidate->DecayVertexV0Z();

                if (IsMC) {
                    if (IsTrueXi1530AOD(Xicandidate,
                                        track1)) {  // MC Association, if it
                                                    // comes from True Xi1530

                        // True Xi1530 signals
                        FillTHnSparse("hInvMass",
                                      {(double)sys, (double)kMCReco,
                                       (double)fCent, vecsum.Pt(), vecsum.M()});
                        if (fQA) {
                            fHistos->FillTH1("hMC_reconstructed_Y",
                                             vecsum.Rapidity());
                            // For cut study
                            fHistos->FillTH1(
                                "hDCADist_Lambda_BTW_Daughters_TrueMC",
                                TMath::Abs(fDCADist_Lambda));
                            fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_TrueMC",
                                             TMath::Abs(fDCADist_Xi));
                            if (Xicandidate->Charge() ==
                                -1) {  // Xi- has +proton, -pion
                                GetImpactParam(pTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC", b[0]);
                                GetImpactParam(nTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC", b[0]);
                            } else {
                                GetImpactParam(pTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC", b[0]);
                                GetImpactParam(nTrackXi, b, bCov);
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC", b[0]);
                            }
                            GetImpactParam(bTrackXi, b, bCov);
                            fHistos->FillTH1(
                                "hDCADist_BachelorPion_to_PV_TrueMC", b[0]);

                            fHistos->FillTH1("hDCADist_lambda_to_PV_TrueMC",
                                             TMath::Abs(fDCADist_Lambda_PV));
                            fHistos->FillTH1("hDCADist_Xi_to_PV_TrueMC",
                                             TMath::Abs(fDCADist_Xi_PV));

                            fHistos->FillTH2("hPhiEta_Xi_TrueMC", Xi_phi,
                                             Xi_eta);
                            fHistos->FillTH2("hLambda_Rxy_TrueMC", LambdaX,
                                             LambdaY);

                            fHistos->FillTH2("hXi_Rxy_TrueMC", lPosXi[0],
                                             lPosXi[1]);
                        }

                    }  // Xi1530 check
                }      // MC
                AliInfo(Form("Passed!! %s", (const char*)SysCheck.at(sys)));
                FillTHnSparse("hInvMass",
                              {(double)sys, (double)sign, (double)fCent,
                               vecsum.Pt(), vecsum.M()});
                if (IsQAInvMass && sys == 0) {
                    if ((int)sign == (int)kData)
                        fHistos->FillTH1("hTotalInvMass_data", vecsum.M());
                    if ((int)sign == (int)kLS)
                        fHistos->FillTH1("hTotalInvMass_LS", vecsum.M());
                }
                if(fExoticFinder2){
                    if (SysCheck.at(sys) == "DefaultOption"){
                        for (UInt_t k = 0; k < ntracks_full; k++) {
                            track2 =
                                (AliVTrack*)fEvt->GetTrack(goodtrackfullindices[k]);
                            if(track1->GetID() == track2->GetID())
                                continue;
                            temp3.SetXYZM(track2->Px(), track2->Py(), track2->Pz(),
                                        pionmass);

                            vecsum2 = vecsum + temp3;  // vecsum = Xi1530, temp3=pion
                            // Y cut
                            if ((vecsum2.Rapidity() > fXi1530RapidityCut_high) ||
                                (vecsum2.Rapidity() < fXi1530RapidityCut_low))
                                continue;

                            // Mass window
                            double mXi1530 = vecsum.M();
                            if (TMath::Abs(mXi1530 - massXi1530) > 0.05)
                                continue;

                            int sign2 = kData;
                            if (track2->Charge() > 0) {
                                if (track1->Charge() > 0)
                                    sign2 = kData;
                                else
                                    sign2 = kLS;
                            }
                            else {
                                if (track1->Charge() > 0)
                                    sign2 = kMixing;
                                else
                                    sign2 = kMCReco;
                            }
                            FillTHnSparse("hInvMass_hf",
                                        {(double)sign2, (double)fCent,
                                        vecsum2.Pt(), vecsum2.M()});
                            if (IsQAInvMass) {
                                if (track2->Charge() > 0) {
                                    if (track1->Charge() > 0)
                                        fHistos->FillTH1("hTotalInvMass_HFpp",
                                                         vecsum2.M());
                                    else
                                        fHistos->FillTH1("hTotalInvMass_HFnp",
                                                         vecsum2.M());
                                } else {
                                    if (track1->Charge() > 0)
                                        fHistos->FillTH1("hTotalInvMass_HFpn",
                                                         vecsum2.M());
                                    else
                                        fHistos->FillTH1("hTotalInvMass_HFnn",
                                                         vecsum2.M());
                                }
                            }
                        }
                    }
                }

                // Fill the QA Histos
                if (fQA) {
                    if (SysCheck.at(sys) == "DefaultOption") {
                        if (IsQAPID)
                            fHistos->FillTH2("hTPCPIDXi1530Pion_cut",
                                             track1->GetTPCmomentum(),
                                             track1->GetTPCsignal());
                        if (Xicandidate->Charge() ==
                            -1) {  // Xi- has +proton, -pion
                            if (IsQAPID) {
                                fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                                 pTrackXi->GetTPCmomentum(),
                                                 pTrackXi->GetTPCsignal());
                                fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                                 nTrackXi->GetTPCmomentum(),
                                                 nTrackXi->GetTPCsignal());
                            }
                            GetImpactParam(pTrackXi, b, bCov);
                            fHistos->FillTH1("hDCADist_LambdaProton_to_PV_cut",
                                             b[0]);
                            GetImpactParam(nTrackXi, b, bCov);
                            fHistos->FillTH1("hDCADist_LambdaPion_to_PV_cut",
                                             b[0]);
                        } else {  // Xi+ has -proton, +pion
                            if (IsQAPID) {
                                fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                                 nTrackXi->GetTPCmomentum(),
                                                 nTrackXi->GetTPCsignal());
                                fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                                 pTrackXi->GetTPCmomentum(),
                                                 pTrackXi->GetTPCsignal());
                            }
                            GetImpactParam(nTrackXi, b, bCov);
                            fHistos->FillTH1("hDCADist_LambdaProton_to_PV_cut",
                                             b[0]);
                            GetImpactParam(pTrackXi, b, bCov);
                            fHistos->FillTH1("hDCADist_LambdaPion_to_PV_cut",
                                             b[0]);
                        }
                        if (IsQAPID) {
                            fHistos->FillTH2("hTPCPIDBachelorPion_cut",
                                             bTrackXi->GetTPCmomentum(),
                                             bTrackXi->GetTPCsignal());

                            // TPC PID Signal
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_cut",
                                             fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_cut",
                                             fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_cut",
                                             fTPCNSigBachelorPion);
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_cut",
                                             fTPCNSigPion);
                        }
                        // DCA QA
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_cut",
                                         fDCADist_Lambda);
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_cut",
                                         fDCADist_Xi);
                        fHistos->FillTH1("hDCADist_lambda_to_PV_cut",
                                         fDCADist_Lambda_PV);
                        fHistos->FillTH1("hDCADist_Xi_to_PV_cut",
                                         TMath::Abs(fDCADist_Xi_PV));
                        fHistos->FillTH1(
                            "hDCADist_BachelorPion_to_PV_cut",
                            TMath::Abs(Xicandidate->DcaBachToPrimVertex()));
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_cut",
                                         pionZ);
                        // CPA QA
                        fHistos->FillTH1("hCosPA_lambda_cut", fLambdaCPA);
                        fHistos->FillTH1("hCosPA_Xi_cut", fXiCPA);

                        // Mass window QA
                        fHistos->FillTH1("hMass_Xi_cut", fMass_Xi);

                        // Eta
                        fHistos->FillTH2("hPhiEta_Xi_cut", Xi_phi, Xi_eta);

                        // XY Radius
                        fHistos->FillTH2("hLambda_Rxy_cut", LambdaX, LambdaY);
                        fHistos->FillTH2("hXi_Rxy_cut", lPosXi[0], lPosXi[1]);
                    }
                    // PID
                    if (IsQAPID) {
                        if (SysCheck.at(sys) == "TPCNsigmaXi1530PionLoose")
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_loose",
                                            fTPCNSigPion);
                        if (SysCheck.at(sys) == "TPCNsigmaXi1530PionTight")
                            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_tight",
                                            fTPCNSigPion);
                        if (SysCheck.at(sys) == "TPCNsigmaXiLoose") {
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_loose",
                                            fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_loose",
                                            fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_loose",
                                            fTPCNSigBachelorPion);
                        }
                        if (SysCheck.at(sys) == "TPCNsigmaXiTight") {
                            fHistos->FillTH1("hTPCPIDsignalLambdaProton_tight",
                                            fTPCNSigProton);
                            fHistos->FillTH1("hTPCPIDsignalLambdaPion_tight",
                                            fTPCNSigLambdaPion);
                            fHistos->FillTH1("hTPCPIDsignalBachelorPion_tight",
                                            fTPCNSigBachelorPion);
                        }
                    }
                    // Xi1530Pion DCA zVetex Check
                    if (SysCheck.at(sys) == "Xi1530PionZVertexLoose")
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_loose",
                                         pionZ);
                    if (SysCheck.at(sys) == "Xi1530PionZVertexTight")
                        fHistos->FillTH1("hDCADist_Xi1530pion_to_PV_tight",
                                         pionZ);

                    // DCA between daughters Check
                    if (SysCheck.at(sys) == "DCADistLambdaDaughtersLoose")
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_loose",
                                         fDCADist_Lambda);
                    if (SysCheck.at(sys) == "DCADistLambdaDaughtersTight")
                        fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_tight",
                                         fDCADist_Lambda);
                    if (SysCheck.at(sys) == "DCADistXiDaughtersLoose")
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_loose",
                                         fDCADist_Xi);
                    if (SysCheck.at(sys) == "DCADistXiDaughtersTight")
                        fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_tight",
                                         fDCADist_Xi);
                    // Lambda DCA zVetex Check
                    if (SysCheck.at(sys) == "DCADistLambdaPVLoose")
                        fHistos->FillTH1("hDCADist_lambda_to_PV_loose",
                                         fDCADist_Lambda_PV);
                    if (SysCheck.at(sys) == "DCADistLambdaPVTight")
                        fHistos->FillTH1("hDCADist_lambda_to_PV_tight",
                                         fDCADist_Lambda_PV);
                    // CPA Check
                    if (SysCheck.at(sys) == "V0CosineOfPointingAngleLoose")
                        fHistos->FillTH1("hCosPA_lambda_loose", fLambdaCPA);
                    if (SysCheck.at(sys) == "V0CosineOfPointingAngleTight")
                        fHistos->FillTH1("hCosPA_lambda_tight", fLambdaCPA);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleLoose")
                        fHistos->FillTH1("hCosPA_Xi_loose", fXiCPA);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight")
                        fHistos->FillTH1("hCosPA_Xi_tight", fXiCPA);

                    // Xi Mass window
                    if (SysCheck.at(sys) == "XiMassWindowLoose")
                        fHistos->FillTH1("hMass_Xi_loose", fMass_Xi);
                    if (SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight")
                        fHistos->FillTH1("hMass_Xi_tight", fMass_Xi);
                }
            }
        }
        AliInfo(Form("Sys check! %.u", sys));
        if ((!fsetsystematics) && (sys == 0))
            break;
    }

    // Event Mixing
    if ((centbin >= 0) && (zbin >= 0) && fsetmixing) {
        eventpool& ep = fEMpool[centbin][zbin];
        Int_t epsize = ep.size();
        if (IsQAEvent) {
            if (epsize < fnMix){
                fHistos->FillTH1("EventQA/hMult_SkippedDataQA", (double)fCent);
                return;
            }
            fHistos->FillTH1("EventQA/hMult_ProcessedDataQA", (double)fCent);
        }

        Int_t nForSkipSameEvent = 0;
        for (auto pool : ep) {
            if (nForSkipSameEvent == (epsize - 1))
                continue;  // same event
            for (auto track : pool)
                trackpool.push_back((AliVTrack*)track);
            nForSkipSameEvent++;
        }
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliAODEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate)
                continue;

            AliAODTrack* pTrackXi = (AliAODTrack*)(Xicandidate->GetDaughter(0));
            AliAODTrack* nTrackXi = (AliAODTrack*)(Xicandidate->GetDaughter(1));
            AliAODTrack* bTrackXi =
                (AliAODTrack*)(Xicandidate->GetDecayVertexXi()->GetDaughter(0));

            if (Xicandidate->ChargeXi() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }
            fTPCNSigBachelorPion = GetTPCnSigma(
                bTrackXi, AliPID::kPion);  // bachelor is always pion

            temp1.SetXYZM(Xicandidate->MomXiX(), Xicandidate->MomXiY(),
                          Xicandidate->MomXiZ(), Xicandidate->MassXi());
            for (UInt_t jt = 0; jt < trackpool.size(); jt++) {
                track1 = (AliVTrack*)trackpool.at(jt);
                if (track1->GetID() == pTrackXi->GetID() ||
                    track1->GetID() == nTrackXi->GetID() ||
                    track1->GetID() == bTrackXi->GetID())
                    continue;

                // PID Cut
                Double_t fTPCNSigPion = GetTPCnSigma(track1, AliPID::kPion);
                if (TMath::Abs(fTPCNSigPion) > fTPCNsigXi1530PionCut)
                    continue;
                if ((TMath::Abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut) ||
                    (TMath::Abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut) ||
                    (TMath::Abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut))
                    continue;

                // Xi1530Pion DCA zVetex Check
                GetImpactParam(track1, b, bCov);
                pionZ = b[1];
                if (pionZ > fXi1530PionZVertexCut)
                    continue;

                // DCA between daughters Check
                Double_t fDCADist_Lambda = TMath::Abs(Xicandidate->DcaV0Daughters());
                Double_t fDCADist_Xi = TMath::Abs(Xicandidate->DcaXiDaughters());
                if (fDCADist_Lambda > fDCADist_LambdaDaughtersCut)
                    continue;
                if (fDCADist_Xi > fDCADist_XiDaughtersCut)
                    continue;
                // DCA Lambda to PV Check
                Double_t fDCADist_Lambda_PV =
                    TMath::Abs(Xicandidate->DcaV0ToPrimVertex());
                if (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut)
                    continue;

                // CPA Check
                Double_t lPosXi[3];
                lPosXi[0] = Xicandidate->DecayVertexXiX();
                lPosXi[1] = Xicandidate->DecayVertexXiY();
                lPosXi[2] = Xicandidate->DecayVertexXiZ();
                Double_t fLambdaCPA;
                if (!fCPAstudy)
                    fLambdaCPA = Xicandidate->CosPointingAngle(lPosXi);
                else
                    fLambdaCPA = Xicandidate->CosPointingAngle(lPosPV);
                Double_t fXiCPA = Xicandidate->CosPointingAngleXi(
                    lPosPV[0], lPosPV[1], lPosPV[2]);

                if (fLambdaCPA < fV0CosineOfPointingAngleCut)
                    continue;
                if (fXiCPA < fCascadeCosineOfPointingAngleCut)
                    continue;

                // Xi Mass Window Check
                Double_t fMass_Xi = Xicandidate->MassXi();
                if (TMath::Abs(fMass_Xi - Ximass) > fXiMassWindowCut)
                    continue;

                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionmass);
                vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
                // Y cut
                if ((vecsum.Rapidity() > fXi1530RapidityCut_high) ||
                    (vecsum.Rapidity() < fXi1530RapidityCut_low))
                    continue;
                // Opening Angle - Not using in normal mode
                if (fExoticFinder) {
                    Double_t angle = temp1.Angle(temp2.Vect());
                    if (TMath::Abs(angle) < 0.0785398)  // 4.5 degree
                        continue;
                }
                // Only use unlike-sign
                if ((Xicandidate->ChargeXi() == -1 && track1->Charge() == -1) ||
                    (Xicandidate->ChargeXi() == +1 && track1->Charge() == +1))
                    continue; 

                FillTHnSparse("hInvMass",
                              {(double)kDefaultOption, (double)kMixing,
                               (double)fCent, vecsum.Pt(), vecsum.M()});
                if (IsQAInvMass)
                    fHistos->FillTH1("hTotalInvMass_Mix", vecsum.M());
                if (fExoticFinder2) {
                    for (UInt_t k = 0; k < ntracks_full; k++) {
                        track2 =
                            (AliVTrack*)fEvt->GetTrack(goodtrackfullindices[k]);
                        if (track1->GetID() == track2->GetID())
                            continue;
                        temp3.SetXYZM(track2->Px(), track2->Py(), track2->Pz(),
                                      pionmass);
                        vecsum2 =
                            vecsum + temp3;  // vecsum = Xi1530, temp3=pion
                        // Y cut
                        if ((vecsum2.Rapidity() > fXi1530RapidityCut_high) ||
                            (vecsum2.Rapidity() < fXi1530RapidityCut_low))
                            continue;
                        
                        // Mass window
                        double mXi1530 = vecsum.M();
                        if (TMath::Abs(mXi1530 - massXi1530) > 0.05)
                            continue;

                        int sign2 = kMCTrue;
                        if (track2->Charge() > 0) {
                            if (track1->Charge() > 0)
                                sign2 = kMCTrue;
                            else
                                sign2 = kMCTruePS;
                        } else {
                            if (track1->Charge() > 0)
                                sign2 = kINEL10;
                            else
                                sign2 = kINELg010;
                        }
                        FillTHnSparse("hInvMass_hf",
                                      {(double)sign2, (double)fCent,
                                       vecsum2.Pt(), vecsum2.M()});
                        if (IsQAInvMass) {
                            if (track2->Charge() > 0) {
                                if (track1->Charge() > 0)
                                    fHistos->FillTH1("hTotalInvMass_HFppMix",
                                                     vecsum2.M());
                                else
                                    fHistos->FillTH1("hTotalInvMass_HFnpMix",
                                                     vecsum2.M());
                            } else {
                                if (track1->Charge() > 0)
                                    fHistos->FillTH1("hTotalInvMass_HFpnMix",
                                                     vecsum2.M());
                                else
                                    fHistos->FillTH1("hTotalInvMass_HFnnMix",
                                                     vecsum2.M());
                            }
                        }
                    }
                }
            }
        }
    }       // mix loop
}
void AliAnalysisTaskXi1530::Terminate(Option_t*) {}

void AliAnalysisTaskXi1530::FillMCinput(AliMCEvent* fMCEvent, Int_t check) {
    // Fill MC input Xi1530 histogram
    // check = 1: INELg0|10
    // check = 2: INEL10
    // check = 3: MB(V0AND)
    // check = 4: After all event cuts
    for (Int_t it = 0; it < fMCEvent->GetNumberOfPrimaries(); it++) {
        TParticle* mcInputTrack =
            (TParticle*)fMCEvent->GetTrack(it)->Particle();
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (TMath::Abs(mcInputTrack->GetPdgCode()) != kXiStarCode)
            continue;
        if (IsPrimaryMC && !mcInputTrack->IsPrimary())
            continue;
        if (fQA)
            fHistos->FillTH1("hMC_generated_Y", mcInputTrack->Y());

        // Y cut
        if ((mcInputTrack->Y() > fXi1530RapidityCut_high) ||
            (mcInputTrack->Y() < fXi1530RapidityCut_low))
            continue;

        if (check == 1)
            FillTHnSparse(
                "hInvMass",
                {(double)kDefaultOption, (double)kINELg010, (double)fCent,
                 mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
        else if (check == 2)
            FillTHnSparse("hInvMass", {(double)kDefaultOption, (double)kINEL10,
                                       (double)fCent, mcInputTrack->Pt(),
                                       mcInputTrack->GetCalcMass()});
        else if (check == 3)
            FillTHnSparse(
                "hInvMass",
                {(double)kDefaultOption, (double)kMCTruePS, (double)fCent,
                 mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
        else if (check == 4)
            FillTHnSparse("hInvMass", {(double)kDefaultOption, (double)kMCTrue,
                                       (double)fCent, mcInputTrack->Pt(),
                                       mcInputTrack->GetCalcMass()});
    }
}
void AliAnalysisTaskXi1530::FillMCinputAOD(AliMCEvent* fMCEvent, Int_t check) {
    // Fill MC input Xi1530 histogram
    // check = 1: INELg0|10
    // check = 2: INEL10
    // check = 3: MB(V0AND)
    // check = 4: After all event cuts

    for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++) {
        AliAODMCParticle* mcInputTrack = (AliAODMCParticle*)fMCArray->At(it);
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (TMath::Abs(mcInputTrack->PdgCode()) != kXiStarCode)
            continue;
        if (IsPrimaryMC && !mcInputTrack->IsPrimary())
            continue;
        if (fQA)
            fHistos->FillTH1("hMC_generated_Y", mcInputTrack->Y());

        // Y cut
        if ((mcInputTrack->Y() > fXi1530RapidityCut_high) ||
            (mcInputTrack->Y() < fXi1530RapidityCut_low))
            continue;

        if (check == 1)
            FillTHnSparse(
                "hInvMass",
                {(double)kDefaultOption, (double)kINELg010, (double)fCent,
                 mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
        else if (check == 2)
            FillTHnSparse("hInvMass", {(double)kDefaultOption, (double)kINEL10,
                                       (double)fCent, mcInputTrack->Pt(),
                                       mcInputTrack->GetCalcMass()});
        else if (check == 3)
            FillTHnSparse(
                "hInvMass",
                {(double)kDefaultOption, (double)kMCTruePS, (double)fCent,
                 mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
        else if (check == 4)
            FillTHnSparse("hInvMass", {(double)kDefaultOption, (double)kMCTrue,
                                       (double)fCent, mcInputTrack->Pt(),
                                       mcInputTrack->GetCalcMass()});
    }
}
void AliAnalysisTaskXi1530::FillMCinputdXi(AliMCEvent* fMCEvent, Int_t check) {
    // Fill MC input Xi1530 histogram
    // check = 1: INEL>0|10
    // check = 2: INEL10
    // check = 3: MB(V0AND)
    // check = 4: After all event cuts
    for (Int_t it = 0; it < fMCEvent->GetNumberOfPrimaries(); it++) {
        TParticle* mcInputTrack =
            (TParticle*)fMCEvent->GetTrack(it)->Particle();
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (!(TMath::Abs(mcInputTrack->GetPdgCode()) == kXiCode))
            continue;

        if (check == 1)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kINELg010, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 2)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kINEL10, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 3)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kMCTruePS, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 4)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kMCTrue, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
    }
}
void AliAnalysisTaskXi1530::FillMCinputdXiAOD(AliMCEvent* fMCEvent,
                                              Int_t check) {
    // Fill MC input Xi1530 histogram
    // check = 1: INEL>0|10
    // check = 2: INEL10
    // check = 3: MB(V0AND)
    // check = 4: After all event cuts

    for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++) {
        AliAODMCParticle* mcInputTrack = (AliAODMCParticle*)fMCArray->At(it);
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (!(TMath::Abs(mcInputTrack->PdgCode()) == kXiCode))
            continue;

        if (check == 1)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kINELg010, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 2)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kINEL10, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 3)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kMCTruePS, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else if (check == 4)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kMCTrue, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
    }
}

THnSparse* AliAnalysisTaskXi1530::CreateTHnSparse(TString name,
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

Long64_t AliAnalysisTaskXi1530::FillTHnSparse(TString name,
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

Long64_t AliAnalysisTaskXi1530::FillTHnSparse(THnSparse* h,
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

TAxis AliAnalysisTaskXi1530::AxisFix(TString name,
                                     int nbin,
                                     Double_t xmin,
                                     Double_t xmax) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(nbin, xmin, xmax);
    axis.SetName(name);
    return axis;
}

TAxis AliAnalysisTaskXi1530::AxisStr(TString name, std::vector<TString> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
    UInt_t i = 1;
    for (auto blabel : bin)
        ax.SetBinLabel(i++, blabel);
    return ax;
}

TAxis AliAnalysisTaskXi1530::AxisVar(TString name, std::vector<Double_t> bin) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(bin.size() - 1, &bin.front());
    axis.SetName(name);
    return axis;
}

TAxis AliAnalysisTaskXi1530::AxisLog(TString name,
                                     int nbin,
                                     Double_t xmin,
                                     Double_t xmax,
                                     Double_t xmin0) {
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    int binoffset = (xmin0 < 0 || (xmin - xmin0) < 1e-9) ? 0 : 1;
    std::vector<Double_t> bin(nbin + 1 + binoffset, 0);
    double logBW3 = (log(xmax) - log(xmin)) / nbin;
    for (int ij = 0; ij <= nbin; ij++)
        bin[ij + binoffset] = xmin * exp(ij * logBW3);
    TAxis axis(nbin, &bin.front());
    axis.SetName(name);
    return axis;
}
Bool_t AliAnalysisTaskXi1530::IsTrueXi1530(AliESDcascade* Xi, AliVTrack* pion) {
    // Check if associated Xi1530 is true Xi1530 in MC set
    if (!Xi)
        return kFALSE;
    if (!pion)
        return kFALSE;

    Bool_t TrueXi1530 = kFALSE;

    AliESDtrack* pTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetPindex()));
    AliESDtrack* nTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetNindex()));
    AliESDtrack* bTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetBindex()));

    TParticle* MCXiD2esd =
        (TParticle*)fMCEvent->GetTrack(TMath::Abs(bTrackXi->GetLabel()))->Particle();
    TParticle* MCLamD1esd;
    TParticle* MCLamD2esd;
    TParticle* MCLamesd;
    TParticle* MCXiesd;
    TParticle* MCXiStaresd;
    TParticle* MCXiStarD2esd;

    if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
        MCLamD1esd = (TParticle*)fMCEvent->GetTrack(TMath::Abs(pTrackXi->GetLabel()))
                         ->Particle();
        MCLamD2esd = (TParticle*)fMCEvent->GetTrack(TMath::Abs(nTrackXi->GetLabel()))
                         ->Particle();
        if (MCLamD1esd->GetMother(0) ==
            MCLamD2esd->GetMother(0)) {  // Same mother(lambda)
            if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                 TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                 TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                     kProtonCode)) {  // Lamda daugthers check #1
                MCLamesd = (TParticle*)fMCEvent
                               ->GetTrack(TMath::Abs(MCLamD1esd->GetMother(0)))
                               ->Particle();
                if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                    kLambdaCode) {  // Lambda check
                    if (MCLamesd->GetMother(0) ==
                        MCXiD2esd->GetMother(
                            0)) {  // Lambda+pion(D2esd) mother check
                        MCXiesd = (TParticle*)fMCEvent
                                      ->GetTrack(TMath::Abs(MCLamesd->GetMother(0)))
                                      ->Particle();
                        if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                            kXiCode) {  // Xi Check
                            MCXiStarD2esd =
                                (TParticle*)fMCEvent
                                    ->GetTrack(TMath::Abs(pion->GetLabel()))
                                    ->Particle();
                            if (MCXiesd->GetMother(0) ==
                                MCXiStarD2esd->GetMother(
                                    0)) {  // Xi+pion mother check
                                MCXiStaresd =
                                    (TParticle*)fMCEvent
                                        ->GetTrack(TMath::Abs(MCXiesd->GetMother(0)))
                                        ->Particle();
                                if (TMath::Abs(MCXiStaresd->GetPdgCode()) ==
                                    kXiStarCode) {  // Xi1530 check
                                    if (IsPrimaryMC) {
                                        if (MCXiStaresd->IsPrimary()) {
                                            TrueXi1530 = kTRUE;
                                        }  // Primary(input) Xi1530 check
                                    } else {
                                        TrueXi1530 = kTRUE;
                                    }
                                }  // Xi1530 check
                            }      // Xi+pion mother check
                        }          // Xi Check
                    }              // Lambda+pion(D2esd) mother check
                }                  // Lambda check
            }                      // Lamda daugthers check
        }                          // Same mother(lambda)
    }                              // D2esd->pion
    return TrueXi1530;
}
Bool_t AliAnalysisTaskXi1530::IsTrueXi1530AOD(AliAODcascade* Xi,
                                              AliVTrack* pion) {
    // Check if associated Xi1530 is true Xi1530 in MC set
    if (!Xi)
        return kFALSE;
    if (!pion)
        return kFALSE;

    Bool_t TrueXi1530 = kFALSE;

    AliAODTrack* pTrackXi = (AliAODTrack*)(Xi->GetDaughter(0));
    AliAODTrack* nTrackXi = (AliAODTrack*)(Xi->GetDaughter(1));
    AliAODTrack* bTrackXi =
        (AliAODTrack*)(Xi->GetDecayVertexXi()->GetDaughter(0));

    AliAODMCParticle* MCXiD2esd =
        (AliAODMCParticle*)fMCArray->At(TMath::Abs(bTrackXi->GetLabel()));
    AliAODMCParticle* MCLamD1esd;
    AliAODMCParticle* MCLamD2esd;
    AliAODMCParticle* MCLamesd;
    AliAODMCParticle* MCXiesd;
    AliAODMCParticle* MCXiStaresd;
    AliAODMCParticle* MCXiStarD2esd;

    if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
        MCLamD1esd = (AliAODMCParticle*)fMCArray->At(TMath::Abs(pTrackXi->GetLabel()));
        MCLamD2esd = (AliAODMCParticle*)fMCArray->At(TMath::Abs(nTrackXi->GetLabel()));
        if (MCLamD1esd->GetMother() ==
            MCLamD2esd->GetMother()) {  // Same mother(lambda)
            if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                 TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                 TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                     kProtonCode)) {  // Lamda daugthers check #1
                MCLamesd = (AliAODMCParticle*)fMCArray->At(
                    TMath::Abs(MCLamD1esd->GetMother()));
                if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                    kLambdaCode) {  // Lambda check
                    if (MCLamesd->GetMother() ==
                        MCXiD2esd
                            ->GetMother()) {  // Lambda+pion(D2esd) mother check
                        MCXiesd = (AliAODMCParticle*)fMCArray->At(
                            TMath::Abs(MCLamesd->GetMother()));
                        if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                            kXiCode) {  // Xi Check
                            MCXiStarD2esd = (AliAODMCParticle*)fMCArray->At(
                                TMath::Abs(pion->GetLabel()));
                            if (MCXiesd->GetMother() ==
                                MCXiStarD2esd
                                    ->GetMother()) {  // Xi+pion mother check
                                MCXiStaresd = (AliAODMCParticle*)fMCArray->At(
                                    TMath::Abs(MCXiesd->GetMother()));
                                if (TMath::Abs(MCXiStaresd->GetPdgCode()) ==
                                    kXiStarCode) {  // Xi1530 check
                                    if (IsPrimaryMC) {
                                        if (MCXiStaresd->IsPrimary()) {
                                            TrueXi1530 = kTRUE;
                                        }  // Primary(input) Xi1530 check
                                    } else {
                                        TrueXi1530 = kTRUE;
                                    }
                                }  // Xi1530 check
                            }      // Xi+pion mother check
                        }          // Xi Check
                    }              // Lambda+pion(D2esd) mother check
                }                  // Lambda check
            }                      // Lamda daugthers check
        }                          // Same mother(lambda)
    }                              // D2esd->pion
    return TrueXi1530;
}
Bool_t AliAnalysisTaskXi1530::IsTrueXi(UInt_t xiIndex) {
    // Check if associated Xi is true Xi in MC set
    AliESDcascade* xiESD;
    AliAODcascade* xiAOD;

    Bool_t TrueXi = kFALSE;
    if (!IsAOD) {
        xiESD = ((AliESDEvent*)fEvt)->GetCascade(xiIndex);
        if (!xiESD)
            return kFALSE;
        AliESDtrack* pTrackXi =
            ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(xiESD->GetPindex()));
        AliESDtrack* nTrackXi =
            ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(xiESD->GetNindex()));
        AliESDtrack* bTrackXi =
            ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(xiESD->GetBindex()));

        TParticle* MCXiD2esd =
            (TParticle*)fMCEvent->GetTrack(TMath::Abs(bTrackXi->GetLabel()))
                ->Particle();
        TParticle* MCLamD1esd;
        TParticle* MCLamD2esd;
        TParticle* MCLamesd;
        TParticle* MCXiesd;

        if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
            MCLamD1esd =
                (TParticle*)fMCEvent->GetTrack(TMath::Abs(pTrackXi->GetLabel()))
                    ->Particle();
            MCLamD2esd =
                (TParticle*)fMCEvent->GetTrack(TMath::Abs(nTrackXi->GetLabel()))
                    ->Particle();
            if (MCLamD1esd->GetMother(0) ==
                MCLamD2esd->GetMother(0)) {  // Same mother(lambda)
                if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                    (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                         kProtonCode)) {  // Lamda daugthers check #1
                    MCLamesd = (TParticle*)fMCEvent
                                   ->GetTrack(TMath::Abs(MCLamD1esd->GetMother(0)))
                                   ->Particle();
                    if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                        kLambdaCode) {  // Lambda check
                        if (MCLamesd->GetMother(0) ==
                            MCXiD2esd->GetMother(
                                0)) {  // Lambda+pion(D2esd) mother check
                            MCXiesd =
                                (TParticle*)fMCEvent
                                    ->GetTrack(TMath::Abs(MCLamesd->GetMother(0)))
                                    ->Particle();
                            if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                                kXiCode) {  // Xi Check
                                TrueXi = kTRUE;
                            }
                        }
                    }
                }
            }
        }
    } else {  // AOD Case
        xiAOD = ((AliAODEvent*)fEvt)->GetCascade(xiIndex);
        if (!xiAOD)
            return kFALSE;
        AliAODTrack* pTrackXi = (AliAODTrack*)(xiAOD->GetDaughter(0));
        AliAODTrack* nTrackXi = (AliAODTrack*)(xiAOD->GetDaughter(1));
        AliAODTrack* bTrackXi =
            (AliAODTrack*)(xiAOD->GetDecayVertexXi()->GetDaughter(0));

        AliAODMCParticle* MCXiD2esd =
            (AliAODMCParticle*)fMCArray->At(TMath::Abs(bTrackXi->GetLabel()));
        AliAODMCParticle* MCLamD1esd;
        AliAODMCParticle* MCLamD2esd;
        AliAODMCParticle* MCLamesd;
        AliAODMCParticle* MCXiesd;

        if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
            MCLamD1esd =
                (AliAODMCParticle*)fMCArray->At(TMath::Abs(pTrackXi->GetLabel()));
            MCLamD2esd =
                (AliAODMCParticle*)fMCArray->At(TMath::Abs(nTrackXi->GetLabel()));
            if (MCLamD1esd->GetMother() ==
                MCLamD2esd->GetMother()) {  // Same mother(lambda)
                if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                    (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                         kProtonCode)) {  // Lamda daugthers check #1
                    MCLamesd = (AliAODMCParticle*)fMCArray->At(
                        TMath::Abs(MCLamD1esd->GetMother()));
                    if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                        kLambdaCode) {  // Lambda check
                        if (MCLamesd->GetMother() ==
                            MCXiD2esd->GetMother()) {  // Lambda+pion(D2esd)
                                                       // mother check
                            MCXiesd = (AliAODMCParticle*)fMCArray->At(
                                TMath::Abs(MCLamesd->GetMother()));
                            if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                                kXiCode) {  // Xi Check
                                TrueXi = kTRUE;
                            }
                        }
                    }
                }
            }
        }
    }
    return TrueXi;
}
double AliAnalysisTaskXi1530::GetTPCnSigma(AliVTrack* track,
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
void AliAnalysisTaskXi1530::GetImpactParam(AliVTrack* track,
                                           Float_t p[2],
                                           Float_t cov[3]) {
    AliNanoAODTrack* nanoT = dynamic_cast<AliNanoAODTrack*>(track);
    if (nanoT)
        nanoT->GetImpactParameters(p[0], p[1]);
    else
        track->GetImpactParameters(p, cov);
}