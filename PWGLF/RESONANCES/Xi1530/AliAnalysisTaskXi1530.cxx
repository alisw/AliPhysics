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
//  Last Modified Date: 2019/03/08
//
////////////////////////////////////////////////////////////////////////////

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliESDcascade.h"
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
#include "TFile.h"
#include "TSystem.h"

// from header
#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "THistManager.h"
//

#include "AliAnalysisTaskXi1530.h"

// Some constants
const Double_t pi = TMath::Pi();
const Double_t pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t Ximass = 1.32171;
enum {
    kData = 1,
    kLS,
    kMixing,
    kMCReco,
    kMCTrue,
    kMCTruePS,
    kAllType
};                                                 // for Physicsl Results
enum { kIsSelected = 1, kPS, kAllNone };           // for V0M signal QA plot
enum { kTrueINELg0 = 1, kTrig, kReco, kGoodVtx };  // for Trigger Efficiency
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
AliAnalysisTaskXi1530RunTable::AliAnalysisTaskXi1530RunTable()
    : fCollisionType(kUnknownCollType) {
    ;
}

AliAnalysisTaskXi1530RunTable::AliAnalysisTaskXi1530RunTable(Int_t runnumber) {
    // Need to be Modified
    if (runnumber >= 256504 && runnumber <= 260014)
        fCollisionType = kPP;  // LHC16kl
    else
        fCollisionType = kUnknownCollType;
}
AliAnalysisTaskXi1530RunTable::~AliAnalysisTaskXi1530RunTable() {
    ;
}

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
    delete fTrackCuts2;
    delete fPIDResponse;
    delete fRunTable;
}
//________________________________________________________________________
void AliAnalysisTaskXi1530::UserCreateOutputObjects() {
    // TrackCuts for Xi1530--------------------------------------------------
    // Primary pion cut(Xi1530pion)
    fTrackCuts = new AliESDtrackCuts();
    fTrackCuts->GetStandardITSTPCTrackCuts2011(kTRUE, kTRUE);
    fTrackCuts->SetEtaRange(-0.8, 0.8);
    fTrackCuts->SetPtRange(0.15, 1e20);

    // secondary particle cut(Xi daugthers)
    fTrackCuts2 = new AliESDtrackCuts();
    // fTrackCuts2 -> GetStandardITSTPCTrackCuts2011(kFALSE,kTRUE); // not
    // primary
    fTrackCuts2->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts2->SetMinNClustersTPC(50);
    fTrackCuts2->SetRequireTPCRefit(kTRUE);
    fTrackCuts2->SetMaxChi2PerClusterTPC(4);
    // fTrackCuts2 -> SetEtaRange(-0.8,0.8);
    fTrackCuts2->SetPtRange(0.15, 1e20);

    // secondary particle cut(Xi daugthers) for systematic study
    fTrackCuts3 = new AliESDtrackCuts();
    fTrackCuts3->GetStandardITSTPCTrackCuts2011(kFALSE, kTRUE);  // not primary
    fTrackCuts3->SetPtRange(0.15, 1e20);
    // ----------------------------------------------------------------------

    fHistos = new THistManager("Xi1530hists");

    auto binType = AxisStr(
        "Type", {"DATA", "LS", "Mixing", "MCReco", "MCTrue", "kMCTruePS"});
    if (!IsMC) {
        if (IsAA && !IsHighMult)
            binCent = AxisFix("Cent", 10, 0, 100);  // for AA study
        else if (!IsHighMult)
            binCent = AxisVar("Cent", {0, 1, 5, 10, 15, 20, 30, 40, 50, 70,
                                       100});  // for kINT7 study
        else
            binCent = AxisVar(
                "Cent", {0, 0.01, 0.03, 0.05, 0.07, 0.1});  // for HM study
    } else
        binCent =
            AxisVar("Cent", {0, 0.01, 0.03, 0.05, 0.07, 0.01, 1, 5, 10, 15, 20,
                             30, 40, 50, 70, 100});  // for kINT7 study

    auto binPt = AxisFix("Pt", 200, 0, 20);
    auto binMass = AxisFix("Mass", 2000, 1.0, 3.0);
    binZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10});
    auto binType_V0M = AxisStr("Type", {"isSelected", "isSelectedPS"});
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
    CreateTHnSparse("hMult", "Multiplicity", 1, {binCent}, "s");
    CreateTHnSparse("hV0MSignal", "V0MSignal", 3,
                    {binType_V0M, binCent, AxisFix("V0MSig", 2000, 0, 2000)},
                    "s");

    auto binTrklet =
        AxisVar("nTrklet", {0, 5, 10, 15, 20, 25, 30, 35, 40, 100});
    if (IsMC) {
        // To get Trigger efficiency in each trk/V0M Multiplicity region
        auto MCType =
            AxisStr("Type", {"TrueINELg0", "Triggered", "Reco", "GoodVtx"});
        // auto binTrklet = AxisVar("nTrklet",{0,5,10,15,20,25,30,35,40,100});
        CreateTHnSparse("htriggered_CINT7", "", 3, {MCType, binCent, binTrklet},
                        "s");  // inv mass distribution of Xi
    }

    std::vector<TString> ent = {
        "All",         "Trigger",     "InCompleteDAQ",
        "No BG",       "No pile up",  "Tracklet>1",
        "Good vertex", "|Zvtx|<10cm", "AliMultSelection",
        "INELg0True"};  // Normal setup
    auto hNofEvt =
        fHistos->CreateTH1("hEventNumbers", "", ent.size(), 0, ent.size());
    for (auto i = 0u; i < ent.size(); i++)
        hNofEvt->GetXaxis()->SetBinLabel(i + 1, ent.at(i).Data());

    // QA Histograms--------------------------------------------------
    //
    if (IsHighMult) {
        fHistos->CreateTH1("hMult_QA", "", 100, 0, 0.1, "s");
        fHistos->CreateTH1("hMult_QA_onlyMult", "", 100, 0, 0.1, "s");
    } else {
        fHistos->CreateTH1("hMult_QA", "", 1000, 0, 100, "s");
        fHistos->CreateTH1("hMult_QA_onlyMult", "", 1000, 0, 100, "s");
    }
    fHistos->CreateTH2("hPhiEta", "", 180, 0, 2 * pi, 40, -2, 2);
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
    fHistos->CreateTH2("hTPCPIDLambdaProton_cut", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("hTPCPIDLambdaPion_cut", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("hTPCPIDBachelorPion_cut", "", 200, 0, 20, 2000, 0, 200);
    fHistos->CreateTH2("hTPCPIDXi1530Pion_cut", "", 200, 0, 20, 2000, 0, 200);
    // Signal
    fHistos->CreateTH1("hTPCPIDsignalLambdaProton_cut", "", 100, -5, 5, "s");
    fHistos->CreateTH1("hTPCPIDsignalLambdaPion_cut", "", 100, -5, 5, "s");
    fHistos->CreateTH1("hTPCPIDsignalBachelorPion_cut", "", 100, -5, 5, "s");
    fHistos->CreateTH1("hTPCPIDsignalXi1530Pion_cut", "", 100, -5, 5, "s");

    // D C A
    // between daughters
    // before
    fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters", "", 300, 0, 3, "s");
    // after
    fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_cut", "", 300, 0, 3, "s");
    fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_cut", "", 300, 0, 3, "s");
    // to PV
    // before
    fHistos->CreateTH1("hDCADist_lambda_to_PV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_Xi_to_PV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_LambdaProton_to_PV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_LambdaPion_to_PV", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_BachelorPion_to_PV", "", 500, 0, 5, "s");
    // after
    fHistos->CreateTH1("hDCADist_lambda_to_PV_cut", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_Xi_to_PV_cut", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_LambdaProton_to_PV_cut", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_LambdaPion_to_PV_cut", "", 500, 0, 5, "s");
    fHistos->CreateTH1("hDCADist_BachelorPion_to_PV_cut", "", 500, 0, 5, "s");

    // C P A
    // before
    fHistos->CreateTH1("hCosPA_lambda", "", 150, 0.85, 1.0, "s");
    fHistos->CreateTH1("hCosPA_Xi", "", 150, 0.85, 1.0, "s");
    // after
    fHistos->CreateTH1("hCosPA_lambda_cut", "", 150, 0.85, 1.0, "s");
    fHistos->CreateTH1("hCosPA_Xi_cut", "", 150, 0.85, 1.0, "s");

    // M a s s   W i n d o w
    fHistos->CreateTH1("hMass_Xi", "", 200, 1.2, 1.4, "s");      // before
    fHistos->CreateTH1("hMass_Xi_cut", "", 200, 1.2, 1.4, "s");  // after

    // E t a
    fHistos->CreateTH2("hPhiEta_Xi", "", 180, 0, 2 * pi, 40, -2, 2);  // before
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

    // Invmass Check
    fHistos->CreateTH1("hTotalInvMass_data", "", 2000, 0.5, 2.5, "s");
    fHistos->CreateTH1("hTotalInvMass_LS", "", 2000, 0.5, 2.5, "s");
    fHistos->CreateTH1("hTotalInvMass_Mix", "", 2000, 0.5, 2.5, "s");

    if (IsMC) {
        // For MC True stduy purpose
        fHistos->CreateTH1("hDCADist_Lambda_BTW_Daughters_TrueMC", "", 300, 0,
                           3, "s");
        fHistos->CreateTH1("hDCADist_Xi_BTW_Daughters_TrueMC", "", 300, 0, 3,
                           "s");
        fHistos->CreateTH1("hDCADist_LambdaProton_to_PV_TrueMC", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("hDCADist_LambdaPion_to_PV_TrueMC", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("hDCADist_BachelorPion_to_PV_TrueMC", "", 500, 0, 5,
                           "s");
        fHistos->CreateTH1("hDCADist_lambda_to_PV_TrueMC", "", 500, 0, 5, "s");
        fHistos->CreateTH1("hDCADist_Xi_to_PV_TrueMC", "", 500, 0, 5, "s");
        fHistos->CreateTH2("hPhiEta_Xi_TrueMC", "", 180, 0, 2 * pi, 40, -2, 2);
        fHistos->CreateTH2("hLambda_Rxy_TrueMC", "", 400, -200, 200, 400, -200,
                           200);
        fHistos->CreateTH2("hXi_Rxy_TrueMC", "", 400, -200, 200, 400, -200,
                           200);

        fHistos->CreateTH1("htriggered_CINT7_true", "", 1000, 0, 1000, "s");
        fHistos->CreateTH1("htriggered_CINT7_trig", "", 1000, 0, 1000, "s");
        fHistos->CreateTH1("htriggered_CINT7_reco", "", 1000, 0, 1000, "s");
        fHistos->CreateTH1("htriggered_CINT7_GoodVtx", "", 1000, 0, 1000, "s");

        fHistos->CreateTH1("hMC_generated_Y", "", 400, -2, 2, "s");
        fHistos->CreateTH1("hMC_reconstructed_Y", "", 400, -2, 2, "s");
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

    // Connect to ESD tree --------------------------------------------------
    event->IsA() == AliESDEvent::Class()
        ? fEvt = dynamic_cast<AliESDEvent*>(event)
        : fEvt = dynamic_cast<AliAODEvent*>(event);
    if (!fEvt)
        return;
    // ----------------------------------------------------------------------

    // Load InputHandler for each event--------------------------------------
    AliInputEventHandler* inputHandler =
        (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();
    // ----------------------------------------------------------------------

    // Preparation for MC ---------------------------------------------------
    if (IsMC) {
        if (fEvt->IsA() == AliESDEvent::Class()) {
            fMCStack = MCEvent()->Stack();
            IsINEL0True = IsMCEventTrueINEL0();
        }  // ESD Case
        else {
            fMCArray = (TClonesArray*)fEvt->FindListObject("mcparticles");
        }  // AOD Case
    }
    // ----------------------------------------------------------------------
    bField = fEvt->GetMagneticField();  // bField for track DCA

    // Event Selection START ************************************************

    // Trigger Check---------------------------------------------------------
    Bool_t IsSelectedTrig;
    if (!IsMC) {
        if (!IsHighMult)
            IsSelectedTrig =
                (inputHandler->IsEventSelected() & AliVEvent::kINT7);
        else
            IsSelectedTrig =
                (inputHandler->IsEventSelected() & AliVEvent::kHighMultV0);
    } else
        IsSelectedTrig = (inputHandler->IsEventSelected() & AliVEvent::kINT7);

    // SPD vs Cluster BG cut-------------------------------------------------
    Bool_t SPDvsClustersBG = kFALSE;
    AliAnalysisUtils* AnalysisUtils = new AliAnalysisUtils();
    if (!AnalysisUtils) {
        AliInfo("No AnalysisUtils Object Found");
        return;
    } else
        SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fEvt);

    // Pile up rejection-----------------------------------------------------
    // Enhanced Pileup rejection method for the higher multiplicity region
    // detail: https://alice-notes.web.cern.ch/node/478
    Bool_t IsNotPileUp = AliPPVsMultUtils::IsNotPileupSPDInMultBins(fEvt);

    // In Complete DAQ Event Cut----------------------------------------------
    Bool_t IncompleteDAQ = fEvt->IsIncompleteDAQ();

    // Multiplicity(centrality) ---------------------------------------------
    ftrackmult = AliESDtrackCuts::GetReferenceMultiplicity(
        ((AliESDEvent*)fEvt), AliESDtrackCuts::kTracklets,
        0.8);                      // tracklet in eta +_0.8
    fCent = GetMultiplicty(fEvt);  // Centrality(AA), Multiplicity(pp)

    // PID response ----------------------------------------------------------
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if (!fPIDResponse)
        AliInfo("No PIDd");
    // -----------------------------------------------------------------------

    // Vertex Check-----------------------------------------------------------
    const AliVVertex* pVtx = fEvt->GetPrimaryVertex();
    // const AliVVertex* trackVtx  = fEvt->GetPrimaryVertexTracks() ;
    const AliVVertex* spdVtx = fEvt->GetPrimaryVertexSPD();
    PVx = pVtx->GetX();
    PVy = pVtx->GetY();
    PVz = pVtx->GetZ();
    fZ = spdVtx->GetZ();

    Bool_t IsGoodVertex =
        spdVtx->GetStatus() &&
        SelectVertex2015pp(((AliESDEvent*)fEvt),
                           kTRUE  // spdVertex->GetDispersion() < 0.04 cm
                                  // && spdVertex->GetZRes()    < 0.25 cm
                           ,
                           kFALSE  // Don't need both trk and spd
                           ,
                           kTRUE);  // z-position difference      < 0.5 cm
    Bool_t IsVtxInZCut = (fabs(fZ) < 10);
    Bool_t IsTrackletinEta1 =
        (AliESDtrackCuts::GetReferenceMultiplicity(
             ((AliESDEvent*)fEvt), AliESDtrackCuts::kTracklets, 1.0) >= 1);

    // Multi Selection--------------------------------------------------------
    // Include:
    //    – INEL>0 selection: At least one SPD tracklet is required within |η| <
    //    1 (IsTrackletinEta1) – Pileup rejection using
    //    AliAnalysisUtils::IsPileUpSPD()                 (IsNotPileUp) – SPD
    //    vertex z resolution < 0.02 cm (IsGoodVertex) -> GetZRes()<0.25 (rough)
    //    – z-position difference between trackand SPD vertex < 0.5 cm
    //    (IsGoodVertex) – vertex z position: |vz| < 10 cm (IsVtxInZCut)
    // Most of them are already choosed in above sections.

    AliMultSelection* MultSelection =
        (AliMultSelection*)fEvt->FindListObject("MultSelection");
    Bool_t IsMultSelcted = MultSelection->IsEventSelected();
    if (IsSelectedTrig && IsMultSelcted)
        fHistos->FillTH1("hMult_QA_onlyMult", (double)fCent);
    // Physics Selection------------------------------------------------------
    IsPS = IsSelectedTrig        // CINT7 Trigger selected
           && !IncompleteDAQ     // No IncompleteDAQ
           && !SPDvsClustersBG   // No SPDvsClusters Background
           && IsNotPileUp        // PileUp rejection
           && IsTrackletinEta1;  // at least 1 tracklet in eta +_1 region.

    IsINEL0Rec =
        IsPS && IsGoodVertex &&
        IsVtxInZCut;  // recontructed INEL > 0 is PS + vtx + Zvtx inside +-10

    // Fill Numver of Events -------------------------------------------------
    fHistos->FillTH1("hEventNumbers", "All", 1);
    if (IsSelectedTrig)
        fHistos->FillTH1("hEventNumbers", "Trigger", 1);
    if (IsSelectedTrig && !IncompleteDAQ)
        fHistos->FillTH1("hEventNumbers", "InCompleteDAQ", 1);
    if (IsSelectedTrig && !IncompleteDAQ && !SPDvsClustersBG)
        fHistos->FillTH1("hEventNumbers", "No BG", 1);
    if (IsSelectedTrig && !IncompleteDAQ && !SPDvsClustersBG && IsNotPileUp)
        fHistos->FillTH1("hEventNumbers", "No pile up", 1);
    if (IsPS)
        fHistos->FillTH1("hEventNumbers", "Tracklet>1", 1);
    if (IsPS && IsGoodVertex)
        fHistos->FillTH1("hEventNumbers", "Good vertex", 1);
    if (IsPS && IsGoodVertex && IsVtxInZCut)
        fHistos->FillTH1("hEventNumbers", "|Zvtx|<10cm", 1);
    if (IsPS && IsGoodVertex && IsVtxInZCut && IsMultSelcted)
        fHistos->FillTH1("hEventNumbers", "AliMultSelection", 1);

    if (IsMC && IsINEL0True)
        fHistos->FillTH1("hEventNumbers", "INELg0True",
                         1);  // For trigger efficiency
    // -----------------------------------------------------------------------

    // ***********************************************************************
    // // Event Selection done

    // Event Mixing pool -----------------------------------------------------
    zbin = binZ.FindBin(fZ) - 1;           // Event mixing z-bin
    centbin = binCent.FindBin(fCent) - 1;  // Event mixing cent bin
    // -----------------------------------------------------------------------

    //  Missing Vetex and Trriger Efficiency ---------------------------------
    if (IsMC) {
        // Int_t trklbin = binTrklet.FindBin(ftrackmult) -1; // # of tracklet
        // bin
        if (IsINEL0True) {
            FillTHnSparse("htriggered_CINT7",
                          {(double)kTrueINELg0, (double)fCent, ftrackmult});
            fHistos->FillTH1("htriggered_CINT7_true", centbin);
        }
        if (IsSelectedTrig) {
            FillTHnSparse("htriggered_CINT7",
                          {(double)kTrig, (double)fCent, ftrackmult});
            fHistos->FillTH1("htriggered_CINT7_trig", centbin);
        }
        if (IsPS) {
            FillTHnSparse("htriggered_CINT7",
                          {(double)kReco, (double)fCent, ftrackmult});
            fHistos->FillTH1("htriggered_CINT7_reco", centbin);
        }
        if (IsINEL0Rec) {
            FillTHnSparse("htriggered_CINT7",
                          {(double)kGoodVtx, (double)fCent, ftrackmult});
            fHistos->FillTH1("htriggered_CINT7_GoodVtx", centbin);
        }
    }
    // -----------------------------------------------------------------------

    // V0M Signal QA --------------------------------------------------------
    // From BeomKyu Kim
    AliVVZERO* lVV0 = fEvt->GetVZEROData();
    Double_t intensity = 0.;
    if (inputHandler->IsEventSelected()) {
        for (int i = 0; i < 64; i++) {
            intensity += lVV0->GetMultiplicity(i);
        }
        FillTHnSparse("hV0MSignal", {kIsSelected, (double)fCent, intensity});
    }
    // ----------------------------------------------------------------------

    // Signal Loss Correction -----------------------------------------------
    if (IsMC && IsSelectedTrig) {
        FillMCinput(fMCStack, kFALSE);
        FillMCinputdXi(fMCStack, kFALSE);
    }
    // ----------------------------------------------------------------------

    // Check tracks and casade, Fill histo************************************
    if ((IsINEL0Rec && IsMultSelcted) ||
        (IsMC &&
         IsINEL0Rec)) {  // In Good Event condition: (IsPS && IsGoodVertex &&
                         // IsVtxInZCut) && IsMultSelcted

        // Draw Multiplicity QA plot in only selected event.
        FillTHnSparse("hMult", {(double)fCent});
        fHistos->FillTH1("hMult_QA", (double)fCent);

        // V0M signal QA
        FillTHnSparse("hV0MSignal", {kPS, (double)fCent, intensity});

        if (IsMC) {
            FillMCinput(
                fMCStack,
                kTRUE);  // Note that MC input Event is filled at this moment.
            FillMCinputdXi(
                fMCStack,
                kTRUE);  // Note that MC input Event is filled at this moment.
        }
        if (this->GoodTracksSelection()  // If Good track
            &&
            this->GoodCascadeSelection())  // and Good cascade is in this event,
            this->FillTracks();            // Fill the histogram
        if (fsetmixing && goodtrackindices.size())
            FillTrackToEventPool();  // use only pion track pool.
    }
    // ***********************************************************************

    PostData(1, fHistos->GetListOfHistograms());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskXi1530::GoodTracksSelection() {
    // Choose Good Tracks from AliESDtracks,
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
    AliVTrack* track;

    for (UInt_t it = 0; it < ntracks; it++) {
        if (fEvt->IsA() == AliESDEvent::Class()) {
            track = (AliESDtrack*)fEvt->GetTrack(it);
            if (!track)
                continue;
            if (!fTrackCuts->AcceptTrack((AliESDtrack*)track))
                continue;
            // if (!track->IsOn(AliVTrack::kITSpureSA)) continue;
            fHistos->FillTH2("hPhiEta", track->Phi(), track->Eta());

            // PID cut for pion
            Double_t fTPCNSigPion =
                fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            fHistos->FillTH2("hTPCPIDXi1530Pion", track->GetTPCmomentum(),
                             track->GetTPCsignal());
            fHistos->FillTH1("hTPCPIDsignalXi1530Pion", fTPCNSigPion);
            if (abs(fTPCNSigPion) > fTPCNsigXi1530PionCut_loose)
                continue;
            fHistos->FillTH2("hTPCPIDXi1530Pion_cut", track->GetTPCmomentum(),
                             track->GetTPCsignal());
            fHistos->FillTH1("hTPCPIDsignalXi1530Pion_cut", fTPCNSigPion);

            // Eta cut
            if (abs(track->Eta()) > fXi1530PionEtaCut)
                continue;

            // Z vertex cut
            if (abs(track->GetZ() - fZ) > fXi1530PionZVertexCut_loose)
                continue;

            // Pion mass window
            // if (fabs(track->M() - pionmass) > 0.007) continue;

            goodtrackindices.push_back(it);
        } else {
            /*
            track = (AliAODTrack*) fEvt ->GetTrack(it);
            if (!track) continue;
            if( ! ((AliAODTrack*) track)->TestFilterBit(fFilterBit)) continue;
            if (track->Pt()<fptcut) continue;
            if (abs(track->Eta())>fetacut) continue;
            fHistos->FillTH2("hPhiEta",track->Phi(),track->Eta());
             */
        }
    }
    return goodtrackindices.size();
}

Bool_t AliAnalysisTaskXi1530::GoodCascadeSelection() {
    // Choose Good Cascade from AliESDcascade
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
    Double_t LambdaX, LambdaY, LambdaZ;
    Double_t fTPCNSigProton, fTPCNSigLambdaPion, fTPCNSigBachelorPion;
    Double_t fDCADist_LambdaProton_PV, fDCADist_LambdaPion_PV;

    fNCascade = 0;
    Bool_t StandardXi = kTRUE;
    for (UInt_t it = 0; it < ncascade; it++) {
        StandardXi = kTRUE;
        if (fEvt->IsA() == AliESDEvent::Class()) {  // ESD case
            Xicandidate = ((AliESDEvent*)fEvt)->GetCascade(it);
            if (!Xicandidate)
                continue;

            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetNindex()))
                continue;
            if (TMath::Abs(Xicandidate->GetPindex()) ==
                TMath::Abs(Xicandidate->GetBindex()))
                continue;
            if (TMath::Abs(Xicandidate->GetNindex()) ==
                TMath::Abs(Xicandidate->GetBindex()))
                continue;

            AliESDtrack* pTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetPindex()));
            AliESDtrack* nTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetNindex()));
            AliESDtrack* bTrackXi =
                ((AliESDEvent*)fEvt)
                    ->GetTrack(TMath::Abs(Xicandidate->GetBindex()));

            // Standard track QA cuts
            if (!fTrackCuts2->AcceptTrack(pTrackXi))
                continue;
            if (!fTrackCuts2->AcceptTrack(nTrackXi))
                continue;
            if (!fTrackCuts2->AcceptTrack(bTrackXi))
                continue;

            // PID cuts for Xi daughters
            if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                fTPCNSigProton =
                    fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion =
                    fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion);
                fHistos->FillTH2("hTPCPIDLambdaProton",
                                 pTrackXi->GetTPCmomentum(),
                                 pTrackXi->GetTPCsignal());
                fHistos->FillTH2("hTPCPIDLambdaPion",
                                 nTrackXi->GetTPCmomentum(),
                                 nTrackXi->GetTPCsignal());

            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton =
                    fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion =
                    fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion);
                fHistos->FillTH2("hTPCPIDLambdaProton",
                                 nTrackXi->GetTPCmomentum(),
                                 nTrackXi->GetTPCsignal());
                fHistos->FillTH2("hTPCPIDLambdaPion",
                                 pTrackXi->GetTPCmomentum(),
                                 pTrackXi->GetTPCsignal());
            }
            fTPCNSigBachelorPion = fPIDResponse->NumberOfSigmasTPC(
                bTrackXi, AliPID::kPion);  // bachelor is always pion
            fHistos->FillTH2("hTPCPIDBachelorPion", bTrackXi->GetTPCmomentum(),
                             bTrackXi->GetTPCsignal());

            fHistos->FillTH1("hTPCPIDsignalLambdaProton", fTPCNSigProton);
            fHistos->FillTH1("hTPCPIDsignalLambdaPion", fTPCNSigLambdaPion);
            fHistos->FillTH1("hTPCPIDsignalBachelorPion", fTPCNSigBachelorPion);

            if (abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut_loose)
                StandardXi = kFALSE;  // PID for proton
            if (abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut_loose)
                StandardXi = kFALSE;  // PID for 1st pion
            if (abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut_loose)
                StandardXi = kFALSE;  // PID for 2nd pion

            // DCA cut
            // DCA between Dautgher particles
            Double_t fDCADist_Lambda = fabs(Xicandidate->GetDcaV0Daughters());
            Double_t fDCADist_Xi = fabs(Xicandidate->GetDcaXiDaughters());
            fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters", fDCADist_Lambda);
            fHistos->FillTH1("hDCADist_Xi_BTW_Daughters", fDCADist_Xi);

            if (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_loose)
                StandardXi = kFALSE;  // DCA proton-pion
            if (fDCADist_Xi > fDCADist_XiDaughtersCut_loose)
                StandardXi = kFALSE;  // DCA Lambda-pion

            // DCA to PV
            Double_t fDCADist_Lambda_PV =
                fabs(Xicandidate->GetD(PVx, PVy, PVz));
            Double_t fDCADist_Xi_PV =
                fabs(Xicandidate->GetDcascade(PVx, PVy, PVz));
            if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                fDCADist_LambdaProton_PV =
                    fabs(pTrackXi->GetD(PVx, PVy, bField));
                fDCADist_LambdaPion_PV = fabs(nTrackXi->GetD(PVx, PVy, bField));
            } else {
                fDCADist_LambdaProton_PV =
                    fabs(nTrackXi->GetD(PVx, PVy, bField));
                fDCADist_LambdaPion_PV = fabs(pTrackXi->GetD(PVx, PVy, bField));
            }
            Double_t fDCADist_BachelorPion_PV =
                fabs(bTrackXi->GetD(PVx, PVy, bField));
            fHistos->FillTH1("hDCADist_lambda_to_PV", fDCADist_Lambda_PV);
            fHistos->FillTH1("hDCADist_Xi_to_PV", fDCADist_Xi_PV);
            fHistos->FillTH1("hDCADist_LambdaProton_to_PV",
                             fDCADist_LambdaProton_PV);
            fHistos->FillTH1("hDCADist_LambdaPion_to_PV",
                             fDCADist_LambdaPion_PV);
            fHistos->FillTH1("hDCADist_BachelorPion_to_PV",
                             fDCADist_BachelorPion_PV);

            if (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_loose)
                StandardXi = kFALSE;  // DCA proton-pion

            // CPA cut
            Double_t fLambdaCPA =
                Xicandidate->GetV0CosineOfPointingAngle(PVx, PVy, PVz);
            Double_t fXiCPA =
                Xicandidate->GetCascadeCosineOfPointingAngle(PVx, PVy, PVz);
            fHistos->FillTH1("hCosPA_lambda", fLambdaCPA);
            fHistos->FillTH1("hCosPA_Xi", fXiCPA);

            if (fLambdaCPA < fV0CosineOfPointingAngleCut_loose)
                StandardXi = kFALSE;
            if (fXiCPA < fCascadeCosineOfPointingAngleCut_loose)
                StandardXi = kFALSE;

            // Mass window cut
            Double_t fMass_Xi = Xicandidate->M();
            fHistos->FillTH1("hMass_Xi", fMass_Xi);
            if (fabs(fMass_Xi - Ximass) > fXiMassWindowCut_loose)
                StandardXi = kFALSE;

            // Eta cut
            if (abs(Xicandidate->Eta()) > fXiEtaCut)
                StandardXi = kFALSE;
            fHistos->FillTH2("hPhiEta_Xi", Xicandidate->Phi(),
                             Xicandidate->Eta());

            // XY Raidus cut(experiemntal)
            Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);
            fHistos->FillTH2("hLambda_Rxy", LambdaX, LambdaY);
            // if(sqrt( pow(LambdaX,2) + pow(LambdaY,2) ) > 100)
            // StandardXi=kFALSE; // NOT USING

            Double_t cX, cY, cZ;
            Xicandidate->GetXYZcascade(cX, cY, cZ);
            fHistos->FillTH2("hXi_Rxy", cX, cY);
            // if(sqrt( pow(cX,2) + pow(cY,2) ) > 100) StandardXi=kFALSE; // NOT
            // USING

            // After selection above
            if (StandardXi) {  // Save only the Xi is good candidate
                FillTHnSparse("hInvMass_dXi",
                              {(int)kData, (double)fCent, Xicandidate->Pt(),
                               Xicandidate->M()});
                if (IsMC)
                    if (IsTrueXi(((AliESDEvent*)fEvt)->GetCascade(it))) {
                        FillTHnSparse("hInvMass_dXi",
                                      {(int)kMCReco, (double)fCent,
                                       Xicandidate->Pt(), Xicandidate->M()});
                    }

                fNCascade++;
                goodcascadeindices.push_back(it);

                // PID QA
                if (Xicandidate->Charge() == -1) {  // Xi- has +proton, -pion
                    fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                } else {  // Xi+ has -proton, +pion
                    fHistos->FillTH2("hTPCPIDLambdaProton_cut",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                    fHistos->FillTH2("hTPCPIDLambdaPion_cut",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                }
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

                // DCA QA
                fHistos->FillTH1("hDCADist_Lambda_BTW_Daughters_cut",
                                 fDCADist_Lambda);
                fHistos->FillTH1("hDCADist_Xi_BTW_Daughters_cut", fDCADist_Xi);
                fHistos->FillTH1("hDCADist_lambda_to_PV_cut",
                                 fDCADist_Lambda_PV);
                fHistos->FillTH1("hDCADist_Xi_to_PV_cut", fDCADist_Xi_PV);
                fHistos->FillTH1("hDCADist_LambdaProton_to_PV_cut",
                                 fDCADist_LambdaProton_PV);
                fHistos->FillTH1("hDCADist_LambdaPion_to_PV_cut",
                                 fDCADist_LambdaPion_PV);
                fHistos->FillTH1("hDCADist_BachelorPion_to_PV_cut",
                                 fDCADist_BachelorPion_PV);

                // CPA QA
                fHistos->FillTH1("hCosPA_lambda_cut", fLambdaCPA);
                fHistos->FillTH1("hCosPA_Xi_cut", fXiCPA);

                // Mass window QA
                fHistos->FillTH1("hMass_Xi_cut", fMass_Xi);

                // Eta
                fHistos->FillTH2("hPhiEta_Xi_cut", Xicandidate->Phi(),
                                 Xicandidate->Eta());

                // XY Radius
                fHistos->FillTH2("hLambda_Rxy_cut", LambdaX, LambdaY);
                fHistos->FillTH2("hXi_Rxy_cut", cX, cY);
            }   // for standard Xi
        }       // ESD case
        else {  // !! NEED TO MODIFY !!
        }       // AOD case
    }           // All Xi loop

    return goodcascadeindices.size();
}

void AliAnalysisTaskXi1530::FillTracks() {
    AliVTrack* track1;           // charged track, pion
    AliESDcascade* Xicandidate;  // Cascade

    TLorentzVector temp1, temp2;
    TLorentzVector vecsum;  // Xi1530 candidate
    Double_t fTPCNSigProton, fTPCNSigLambdaPion, fTPCNSigBachelorPion;
    Double_t fDCADist_LambdaProton_PV, fDCADist_LambdaPion_PV;

    // The following CovMatrix is set so that PropogateToDCA() ignores track
    // errors. Only used to propagate Xi to third pion for XiStar reconstruction
    // Origin: AliPhysics/PWGLF/RESONANCES/extra/AliXiStar.cxx (Dhevan
    // Gangadharan)
    Double_t fCovMatrix[21], xiVtx[3], xiP[3], PiX[3];
    for (Int_t i = 0; i < 21; i++)
        fCovMatrix[i] = 0;
    fCovMatrix[0] = 1, fCovMatrix[2] = 1, fCovMatrix[5] = 1, fCovMatrix[9] = 1,
    fCovMatrix[14] = 1, fCovMatrix[20] = 1;
    AliESDtrack* fXiTrack = new AliESDtrack();  // As a ESD Track

    const UInt_t ncascade = goodcascadeindices.size();
    const UInt_t ntracks = goodtrackindices.size();

    tracklist trackpool;
    if (fsetmixing) {
        eventpool& ep = fEMpool[centbin][zbin];
        if ((int)ep.size() < (int)fnMix)
            return;
        for (auto pool : ep) {
            for (auto track : pool)
                trackpool.push_back((AliVTrack*)track);
        }
    }

    for (UInt_t sys = 0; sys < (UInt_t)binSystematics.GetNbins(); sys++) {
        // Systematic study loop.
        // sys = 0 -> Default cut option
        // for more details, please check "SysCheck" in header file.
        AliInfo(Form("Sys check! %s", (const char*)SysCheck.at(sys)));
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate)
                continue;

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
                fTPCNSigProton =
                    fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion =
                    fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion);
            } else {  // Xi+ has -proton, +pion
                fTPCNSigProton =
                    fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton);
                fTPCNSigLambdaPion =
                    fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion);
            }
            fTPCNSigBachelorPion = fPIDResponse->NumberOfSigmasTPC(
                bTrackXi, AliPID::kPion);  // bachelor is always pion

            temp1.SetXYZM(Xicandidate->Px(), Xicandidate->Py(),
                          Xicandidate->Pz(), Ximass);

            // for PropogateToDCA
            xiVtx[0] = Xicandidate->Xv();
            xiVtx[1] = Xicandidate->Yv();
            xiVtx[2] = Xicandidate->Zv();
            xiP[0] = Xicandidate->Px();
            xiP[1] = Xicandidate->Py();
            xiP[2] = Xicandidate->Pz();
            fXiTrack->Set(xiVtx, xiP, fCovMatrix,
                          Short_t(Xicandidate->Charge()));

            for (UInt_t j = 0; j < ntracks; j++) {
                track1 = (AliVTrack*)fEvt->GetTrack(goodtrackindices[j]);
                if (!track1)
                    continue;

                if (track1->GetID() == pTrackXi->GetID() ||
                    track1->GetID() == nTrackXi->GetID() ||
                    track1->GetID() == bTrackXi->GetID())
                    continue;

                // PID Cut Systematic check
                // -------------------------------------------------
                Double_t fTPCNSigPion =
                    fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);

                // Xi1530Pion PID
                if ((SysCheck.at(sys) != "TPCNsigmaXi1530PionLoose") &&
                    (abs(fTPCNSigPion) > fTPCNsigXi1530PionCut))
                    continue;
                if ((SysCheck.at(sys) == "TPCNsigmaXi1530PionTight") &&
                    (abs(fTPCNSigPion) > fTPCNsigXi1530PionCut_tight))
                    continue;
                // Xi PID
                if (SysCheck.at(sys) != "TPCNsigmaXiLoose") {
                    if (abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut)
                        continue;
                    if (abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
                        continue;
                    if (abs(fTPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
                        continue;
                }
                if (SysCheck.at(sys) != "TPCNsigmaXiTight") {
                    if (abs(fTPCNSigProton) > fTPCNsigLambdaProtonCut_tight)
                        continue;
                    if (abs(fTPCNSigLambdaPion) > fTPCNsigLambdaPionCut_tight)
                        continue;
                    if (abs(fTPCNSigBachelorPion) >
                        fTPCNsigBachelorPionCut_tight)
                        continue;
                }

                // Xi1530Pion DCA zVetex Check
                if ((SysCheck.at(sys) != "Xi1530PionZVertexLoose") &&
                    (abs(track1->GetZ() - fZ) > fXi1530PionZVertexCut))
                    continue;
                if ((SysCheck.at(sys) == "Xi1530PionZVertexTight") &&
                    (abs(track1->GetZ() - fZ) > fXi1530PionZVertexCut_tight))
                    continue;

                // DCA between daughters Check
                Double_t fDCADist_Lambda =
                    fabs(Xicandidate->GetDcaV0Daughters());
                Double_t fDCADist_Xi = fabs(Xicandidate->GetDcaXiDaughters());
                if ((SysCheck.at(sys) != "DCADistLambdaDaughtersLoose") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut))
                    continue;
                if ((SysCheck.at(sys) == "DCADistLambdaDaughtersTight") &&
                    (fDCADist_Lambda > fDCADist_LambdaDaughtersCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "DCADistXiDaughtersLoose") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut))
                    continue;
                if ((SysCheck.at(sys) == "DCADistXiDaughtersTight") &&
                    (fDCADist_Xi > fDCADist_XiDaughtersCut_tight))
                    continue;

                // DCA Lambda to PV Check
                Double_t fDCADist_Lambda_PV =
                    fabs(Xicandidate->GetD(PVx, PVy, PVz));
                if ((SysCheck.at(sys) != "DCADistLambdaPVLoose") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut))
                    continue;
                if ((SysCheck.at(sys) == "DCADistLambdaPVTight") &&
                    (fDCADist_Lambda_PV < fDCADist_Lambda_PVCut_tight))
                    continue;

                // CPA Check
                Double_t fLambdaCPA =
                    Xicandidate->GetV0CosineOfPointingAngle(PVx, PVy, PVz);
                Double_t fXiCPA =
                    Xicandidate->GetCascadeCosineOfPointingAngle(PVx, PVy, PVz);

                if ((SysCheck.at(sys) != "V0CosineOfPointingAngleLoose") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut))
                    continue;
                if ((SysCheck.at(sys) == "V0CosineOfPointingAngleTight") &&
                    (fLambdaCPA < fV0CosineOfPointingAngleCut_tight))
                    continue;
                if ((SysCheck.at(sys) != "CascadeCosineOfPointingAngleLoose") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut))
                    continue;
                if ((SysCheck.at(sys) == "CascadeCosineOfPointingAngleTight") &&
                    (fXiCPA < fCascadeCosineOfPointingAngleCut_tight))
                    continue;

                // Xi Mass Window Check
                Double_t fMass_Xi = Xicandidate->M();
                if ((SysCheck.at(sys) != "XiMassWindowLoose") &&
                    (fabs(fMass_Xi - Ximass) > fXiMassWindowCut))
                    continue;
                if ((SysCheck.at(sys) == "XiMassWindowTight") &&
                    (fabs(fMass_Xi - Ximass) > fXiMassWindowCut_tight))
                    continue;

                // XiTrack Cut Systematic check
                // ---------------------------------------------
                if (SysCheck.at(sys) == "XiTrackCut") {
                    if (!fTrackCuts3->AcceptTrack(pTrackXi))
                        continue;
                    if (!fTrackCuts3->AcceptTrack(nTrackXi))
                        continue;
                    if (!fTrackCuts3->AcceptTrack(bTrackXi))
                        continue;
                }
                temp2.SetXYZM(track1->Px(), track1->Py(), track1->Pz(),
                              pionmass);

                vecsum = temp1 + temp2;  // temp1 = cascade, temp2=pion
                // Y cut
                if (fabs(vecsum.Rapidity()) > fXi1530RapidityCut)
                    continue;

                // PropagateToDCA cut
                track1->GetXYZ(PiX);
                AliVertex* XiStarVtx = new AliVertex(PiX, 0, 0);
                if (!(fXiTrack->PropagateToDCA(XiStarVtx, bField, 3)))
                    continue;

                auto sign = kAllType;
                if ((Xicandidate->Charge() == -1 && track1->Charge() == +1) ||
                    (Xicandidate->Charge() == +1 && track1->Charge() == -1))
                    sign = kData;  // Unlike sign -> Data
                else
                    sign = kLS;  // like sign bg

                if (IsMC) {
                    if (fEvt->IsA() == AliESDEvent::Class()) {  // ESD case
                        if (IsTrueXi1530(Xicandidate,
                                         track1)) {  // MC Association, if it
                                                     // comes from True Xi1530

                            // True Xi1530 signals
                            FillTHnSparse(
                                "hInvMass",
                                {(double)sys, (double)kMCReco, (double)fCent,
                                 vecsum.Pt(), vecsum.M()});
                            fHistos->FillTH1("hMC_reconstructed_Y",
                                             vecsum.Rapidity());
                            // For cut study
                            fHistos->FillTH1(
                                "hDCADist_Lambda_BTW_Daughters_TrueMC",
                                fabs(Xicandidate->GetDcaV0Daughters()));
                            fHistos->FillTH1(
                                "hDCADist_Xi_BTW_Daughters_TrueMC",
                                fabs(Xicandidate->GetDcaXiDaughters()));
                            if (Xicandidate->Charge() ==
                                -1) {  // Xi- has +proton, -pion
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC",
                                    fabs(pTrackXi->GetD(PVx, PVy, bField)));
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC",
                                    fabs(nTrackXi->GetD(PVx, PVy, bField)));
                            } else {
                                fHistos->FillTH1(
                                    "hDCADist_LambdaProton_to_PV_TrueMC",
                                    fabs(nTrackXi->GetD(PVx, PVy, bField)));
                                fHistos->FillTH1(
                                    "hDCADist_LambdaPion_to_PV_TrueMC",
                                    fabs(pTrackXi->GetD(PVx, PVy, bField)));
                            }
                            fHistos->FillTH1(
                                "hDCADist_BachelorPion_to_PV_TrueMC",
                                fabs(bTrackXi->GetD(PVx, PVy, bField)));

                            fHistos->FillTH1(
                                "hDCADist_lambda_to_PV_TrueMC",
                                fabs(Xicandidate->GetD(PVx, PVy, PVz)));
                            fHistos->FillTH1(
                                "hDCADist_Xi_to_PV_TrueMC",
                                fabs(Xicandidate->GetDcascade(PVx, PVy, PVz)));

                            fHistos->FillTH2("hPhiEta_Xi_TrueMC",
                                             Xicandidate->Phi(),
                                             Xicandidate->Eta());

                            Double_t LambdaX, LambdaY, LambdaZ;
                            Xicandidate->GetXYZ(LambdaX, LambdaY, LambdaZ);

                            fHistos->FillTH2("hLambda_Rxy_TrueMC", LambdaX,
                                             LambdaY);

                            Double_t cX, cY, cZ;
                            Xicandidate->GetXYZcascade(cX, cY, cZ);
                            fHistos->FillTH2("hXi_Rxy_TrueMC", cX, cY);
                        }   // Xi1530 check
                    }       // MC ESD
                    else {  // !! NEED TO UPDATE FOR AOD CASE !!
                        //
                    }  // MC AOD
                }      // MC
                FillTHnSparse("hInvMass",
                              {(double)sys, (double)sign, (double)fCent,
                               vecsum.Pt(), vecsum.M()});
                if (sys == 0) {
                    if ((int)sign == (int)kData)
                        fHistos->FillTH1("hTotalInvMass_data", vecsum.M());
                    if ((int)sign == (int)kLS)
                        fHistos->FillTH1("hTotalInvMass_LS", vecsum.M());
                }
            }
        }
        AliInfo(Form("Sys check! %.u", sys));
        if ((!fsetsystematics) && (sys == 0))
            break;
    }

    // Event Mixing
    if (fsetmixing) {
        for (UInt_t i = 0; i < ncascade; i++) {
            Xicandidate =
                ((AliESDEvent*)fEvt)->GetCascade(goodcascadeindices[i]);
            if (!Xicandidate)
                continue;
            temp1.SetXYZM(Xicandidate->Px(), Xicandidate->Py(),
                          Xicandidate->Pz(), Xicandidate->M());

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

                if (fabs(vecsum.Rapidity()) > fXi1530RapidityCut)
                    continue;  // rapidity cut

                FillTHnSparse("hInvMass",
                              {(double)kDefaultOption, (double)kMixing,
                               (double)fCent, vecsum.Pt(), vecsum.M()});
                fHistos->FillTH1("hTotalInvMass_Mix", vecsum.M());
            }
        }
    }  // mix loop
}

void AliAnalysisTaskXi1530::Terminate(Option_t*) {}

Bool_t AliAnalysisTaskXi1530::SelectVertex2015pp(
    AliESDEvent* esd,
    Bool_t checkSPDres,       // enable check on vtx resolution
    Bool_t requireSPDandTrk,  // ask for both trk and SPD vertex
    Bool_t checkProximity)    // apply cut on relative position of spd and trk
                            // verteces
{
    // From
    // AliPhysics/PWGLF/SPECTRA/ChargedHadrons/dNdPtVsMultpp/AliAnalysisTaskPPvsMultINEL0.cxx
    // Original author: Sergio Iga
    if (!esd)
        return kFALSE;

    const AliESDVertex* trkVertex = esd->GetPrimaryVertexTracks();
    const AliESDVertex* spdVertex = esd->GetPrimaryVertexSPD();
    Bool_t hasSPD = spdVertex->GetStatus();
    Bool_t hasTrk = trkVertex->GetStatus();

    // Note that AliVertex::GetStatus checks that N_contributors is > 0
    // reject events if both are explicitly requested and none is available
    if (requireSPDandTrk && !(hasSPD && hasTrk))
        return kFALSE;

    // reject events if none between the SPD or track verteces are available
    // if no trk vertex, try to fall back to SPD vertex;
    if (!hasTrk) {
        if (!hasSPD)
            return kFALSE;
        // on demand check the spd vertex resolution and reject if not satisfied
        if (checkSPDres && !IsGoodSPDvertexRes(spdVertex))
            return kFALSE;
    } else {
        if (hasSPD) {
            // if enabled check the spd vertex resolution and reject if not
            // satisfied if enabled, check the proximity between the spd vertex
            // and trak vertex, and reject if not satisfied
            if (checkSPDres && !IsGoodSPDvertexRes(spdVertex))
                return kFALSE;
            if ((checkProximity &&
                 TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ()) > 0.5))
                return kFALSE;
        }
    }
    return kTRUE;
}
Bool_t AliAnalysisTaskXi1530::IsGoodSPDvertexRes(
    const AliESDVertex* spdVertex) {
    // From
    // AliPhysics/PWGLF/SPECTRA/ChargedHadrons/dNdPtVsMultpp/AliAnalysisTaskPPvsMultINEL0.cxx
    // Original author: Sergio Iga
    if (!spdVertex)
        return kFALSE;
    if (spdVertex->IsFromVertexerZ() &&
        !(spdVertex->GetDispersion() < 0.04 && spdVertex->GetZRes() < 0.25))
        return kFALSE;
    return kTRUE;
}
Double_t AliAnalysisTaskXi1530::GetMultiplicty(AliVEvent* fEvt) {
    // Set multiplicity value
    // fCenttemp:
    //       0-100: Selected, value.
    //       999: Not selected
    //       -999: No MultSection
    //
    Double_t fCenttemp = -999;
    AliMultSelection* MultSelection =
        (AliMultSelection*)fEvt->FindListObject("MultSelection");
    if (MultSelection) {
        fCenttemp = MultSelection->GetMultiplicityPercentile("V0M");
    } else {
        // If this happens, re-check if AliMultSelectionTask ran before your
        // task!
        AliInfo("Didn't find MultSelection!");
        fCenttemp = 999;
    }
    return fCenttemp;
}
void AliAnalysisTaskXi1530::FillMCinput(AliStack* fMCStack, Bool_t PS) {
    // Fill MC input Xi1530 histogram
    for (Int_t it = 0; it < fMCStack->GetNprimary(); it++) {
        TParticle* mcInputTrack = (TParticle*)fMCStack->Particle(it);
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (abs(mcInputTrack->GetPdgCode()) != kXiStarCode)
            continue;
        if (IsPrimaryMC && !mcInputTrack->IsPrimary())
            continue;
        fHistos->FillTH1("hMC_generated_Y", mcInputTrack->Y());

        // Y cut
        if (fabs(mcInputTrack->Y()) > fXi1530RapidityCut)
            continue;

        if (PS)
            FillTHnSparse(
                "hInvMass",
                {(double)kDefaultOption, (double)kMCTruePS, (double)fCent,
                 mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
        else
            FillTHnSparse("hInvMass", {(double)kDefaultOption, (double)kMCTrue,
                                       (double)fCent, mcInputTrack->Pt(),
                                       mcInputTrack->GetCalcMass()});
    }
}
void AliAnalysisTaskXi1530::FillMCinputdXi(AliStack* fMCStack, Bool_t PS) {
    // Fill MC input Xi1530 histogram
    for (Int_t it = 0; it < fMCStack->GetNprimary(); it++) {
        TParticle* mcInputTrack = (TParticle*)fMCStack->Particle(it);
        if (!mcInputTrack) {
            Error("UserExec", "Could not receive MC track %d", it);
            continue;
        }
        if (!(abs(mcInputTrack->GetPdgCode()) == kXiCode))
            continue;
        if (PS)
            FillTHnSparse("hInvMass_dXi",
                          {(double)kMCTruePS, (double)fCent, mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
        else
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

Bool_t AliAnalysisTaskXi1530::IsMCEventTrueINEL0() {
    // From
    // AliPhysics/PWGLF/SPECTRA/ChargedHadrons/dNdPtVsMultpp/AliAnalysisTaskPPvsMultINEL0.cxx
    // Original author: Sergio Iga

    Bool_t isINEL0 = kFALSE;
    for (int iT = 0; iT < fMCStack->GetNtrack(); iT++) {
        TParticle* mcParticle = fMCStack->Particle(iT);
        if (!mcParticle) {
            AliInfo("no mcParticle");
            continue;
        }
        if (!fMCStack->IsPhysicalPrimary(iT))
            continue;
        if (!(mcParticle->Pt() > 0.0))
            continue;
        if (TMath::Abs(mcParticle->Eta()) > 1.0)
            continue;
        if (!(TMath::Abs(mcParticle->GetPDG()->Charge()) == 3))
            continue;
        isINEL0 = kTRUE;
        break;
    }
    return isINEL0;
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
        (TParticle*)fMCStack->Particle(abs(bTrackXi->GetLabel()));
    TParticle* MCLamD1esd;
    TParticle* MCLamD2esd;
    TParticle* MCLamesd;
    TParticle* MCXiesd;
    TParticle* MCXiStaresd;
    TParticle* MCXiStarD2esd;

    if (abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
        MCLamD1esd = (TParticle*)fMCStack->Particle(abs(pTrackXi->GetLabel()));
        MCLamD2esd = (TParticle*)fMCStack->Particle(abs(nTrackXi->GetLabel()));
        if (MCLamD1esd->GetMother(0) ==
            MCLamD2esd->GetMother(0)) {  // Same mother(lambda)
            if ((abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                 abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                (abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                 abs(MCLamD2esd->GetPdgCode()) ==
                     kProtonCode)) {  // Lamda daugthers check #1
                MCLamesd = (TParticle*)fMCStack->Particle(
                    abs(MCLamD1esd->GetMother(0)));
                if (abs(MCLamesd->GetPdgCode()) ==
                    kLambdaCode) {  // Lambda check
                    if (MCLamesd->GetMother(0) ==
                        MCXiD2esd->GetMother(
                            0)) {  // Lambda+pion(D2esd) mother check
                        MCXiesd = (TParticle*)fMCStack->Particle(
                            abs(MCLamesd->GetMother(0)));
                        if (abs(MCXiesd->GetPdgCode()) ==
                            kXiCode) {  // Xi Check
                            MCXiStarD2esd = (TParticle*)fMCStack->Particle(
                                abs(pion->GetLabel()));
                            if (MCXiesd->GetMother(0) ==
                                MCXiStarD2esd->GetMother(
                                    0)) {  // Xi+pion mother check
                                MCXiStaresd = (TParticle*)fMCStack->Particle(
                                    abs(MCXiesd->GetMother(0)));
                                if (abs(MCXiStaresd->GetPdgCode()) ==
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
Bool_t AliAnalysisTaskXi1530::IsTrueXi(AliESDcascade* Xi) {
    // Check if associated Xi1530 is true Xi1530 in MC set
    if (!Xi)
        return kFALSE;

    Bool_t TrueXi = kFALSE;

    AliESDtrack* pTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetPindex()));
    AliESDtrack* nTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetNindex()));
    AliESDtrack* bTrackXi =
        ((AliESDEvent*)fEvt)->GetTrack(TMath::Abs(Xi->GetBindex()));

    TParticle* MCXiD2esd =
        (TParticle*)fMCStack->Particle(abs(bTrackXi->GetLabel()));
    TParticle* MCLamD1esd;
    TParticle* MCLamD2esd;
    TParticle* MCLamesd;
    TParticle* MCXiesd;

    if (abs(MCXiD2esd->GetPdgCode()) == kPionCode) {  // D2esd->pion
        MCLamD1esd = (TParticle*)fMCStack->Particle(abs(pTrackXi->GetLabel()));
        MCLamD2esd = (TParticle*)fMCStack->Particle(abs(nTrackXi->GetLabel()));
        if (MCLamD1esd->GetMother(0) ==
            MCLamD2esd->GetMother(0)) {  // Same mother(lambda)
            if ((abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                 abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                (abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                 abs(MCLamD2esd->GetPdgCode()) ==
                     kProtonCode)) {  // Lamda daugthers check #1
                MCLamesd = (TParticle*)fMCStack->Particle(
                    abs(MCLamD1esd->GetMother(0)));
                if (abs(MCLamesd->GetPdgCode()) ==
                    kLambdaCode) {  // Lambda check
                    if (MCLamesd->GetMother(0) ==
                        MCXiD2esd->GetMother(
                            0)) {  // Lambda+pion(D2esd) mother check
                        MCXiesd = (TParticle*)fMCStack->Particle(
                            abs(MCLamesd->GetMother(0)));
                        if (abs(MCXiesd->GetPdgCode()) ==
                            kXiCode) {  // Xi Check
                            TrueXi = kTRUE;
                        }
                    }
                }
            }
        }
    }
    return TrueXi;
}

void AliAnalysisTaskXi1530::FillTrackToEventPool() {
    // Fill Selected tracks to event mixing pool
    AliVTrack* goodtrack;

    tracklist* etl;
    eventpool* ep;
    // Event mixing pool

    ep = &fEMpool[centbin][zbin];
    ep->push_back(tracklist());
    etl = &(ep->back());
    // Fill selected tracks
    for (UInt_t i = 0; i < goodtrackindices.size(); i++) {
        goodtrack = (AliESDtrack*)fEvt->GetTrack(goodtrackindices[i]);
        if (!goodtrack)
            continue;
        etl->push_back((AliVTrack*)goodtrack->Clone());
    }
    if (!goodtrackindices.size())
        ep->pop_back();
    if ((int)ep->size() > (int)fnMix) {
        for (auto it : ep->front())
            delete it;
        ep->pop_front();
    }
}
