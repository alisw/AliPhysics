/*************************************************************************
 *Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                       *
 *Author: The ALICE Off-line Project.                                    *
 *Contributors are mentioned in the code where appropriate.              *
 *                                                                       *
 *Permission to use, copy, modify and distribute this software and its   *
 *documentation strictly for non-commercial purposes is hereby granted   *
 *without fee, provided that the above copyright notice appears in all   *
 *copies and that both the copyright notice and this permission notice   *
 *appear in the supporting documentation. The authors make no claims     *
 *about the suitability of this software for any purpose. It is          *
 *provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*AliAnalysisTaskXi1530PbPb
 *
 *Test code for the reconstructing Xi(1530)^{0} by using pion + Xi
 *Output could be saved to nTuple by using SetFillnTuple(kTRUE)
 *   -> can be used for TMVA input
 *
 *Author: Bong-Hwi Lim
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
#include "AliKFVertex.h"
#include "AliMultSelectionTask.h"
#include "AliPIDResponse.h"
#include "TChain.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAnalysisTaskTrackMixer.h"
#include "AliAnalysisTaskXi1530PbPb.h"
#include "AliESDcascade.h"
#include "AliESDv0.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t pi = TMath::Pi();
const Double_t pionMass = AliPID::ParticleMass(AliPID::kPion);
const Double_t v0Mass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
const Double_t Ximass = TDatabasePDG::Instance()->GetParticle(kXiMinus)->Mass();

enum
{
    kNormal = 1,
    kAnti = 2
};
enum
{
    kSigmaStarPCode = 3224, // Sigma(1385)+
    kSigmaStarNCode = 3114, // Sigma(1385)-
    kXiCode = 3312,         // Xi-
    kXiZeroCode = 3322,     // Xi0
    kLambdaCode = 3122,     // Lambda
    kProtonCode = 2212,     // Proton+
    kPionCode = 211,        // Pion+
    kXiStarCode = 3324,     // Xi(1530)0
    kXiStarMCode = 3314,    // Xi(1530)-
    kLambdaStarCode = 3124, // Lambda1520
    kSigmaZeroCode = 3212   // Sigma0
};
enum
{
    kXiStar = 1,
    kXiStar_LS,
    kXiStar_MIX,
    kXiStar_NOT,
    kAllType
};
enum
{
    kXiStar_GEN = 1, // 1
    kXiStar_GEN_INEL10,
    kXiStar_GEN_INEL10_IGZ,
    kXiStar_GEN_TRIG,
    kXiStar_REC,
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

class AliAnalysisTaskXi1530PbPb;

ClassImp(AliAnalysisTaskXi1530PbPb)
    AliAnalysisTaskXi1530PbPb::AliAnalysisTaskXi1530PbPb() : AliAnalysisTaskSE(),
                                                             fEventCuts(),
                                                             fCheckTPCGeo(kFALSE),
                                                             fTPCActiveLengthCutDeltaY(3.0),
                                                             fTPCActiveLengthCutDeltaZ(220.0),
                                                             fRequireCutGeoNcrNclLength(130),
                                                             fRequireCutGeoNcrNclGeom1Pt(1.5),
                                                             fCutGeoNcrNclFractionNcr(0.85),
                                                             fCutGeoNcrNclFractionNcl(0.7),
                                                             fTrackCuts(nullptr),
                                                             fPIDResponse(),
                                                             fMixingPool(),
                                                             fEvt(nullptr),
                                                             fMCEvent(nullptr),
                                                             fHistos(nullptr),
                                                             fVertex(nullptr),
                                                             fMCArray(nullptr),
                                                             fIsAOD(kFALSE),
                                                             fIsNano(kFALSE),
                                                             fSetMixing(kFALSE),
                                                             fUseBuiltinMixer(kFALSE),
                                                             fFillQAPlot(kTRUE),
                                                             fIsMC(kFALSE),
                                                             fIsPrimaryMC(kFALSE),
                                                             fFillTree(kFALSE),
                                                             fIsINEL(kFALSE),
                                                             fIsHM(kFALSE),
                                                             fLambdaCPAtoXi(kFALSE),
                                                             fExoticFinder(kFALSE),
                                                             fSkipFillingHisto(kFALSE),
                                                             fEMpool(0),
                                                             fBinCent(),
                                                             fBinZ(),
                                                             fPosPV(),
                                                             fMagField(0),
                                                             fCent(-1),
                                                             fnMix(10),
                                                             fCentBin(-1),
                                                             fZbin(-1),
                                                             fFilterBit(32.0),
                                                             fTPCNsigXi1530PionCut(3.0),
                                                             fXi1530PionEtaCut(0.8),
                                                             fXi1530PionPVzCut(2.5),
                                                             fXi1530PionPVrSigmaCut(7.0),
                                                             fTPCNsigLambdaProtonCut(3.0),
                                                             fTPCNsigLambdaPionCut(3.0),
                                                             fTPCNsigBachelorPionCut(3.0),
                                                             fDCALambdaDaughtersCut(0.8),
                                                             fDCAXiDaughtersCut(0.6),
                                                             fMinDCALambdaPVCut(0.00),
                                                             fMaxDCALambdaPVCut(10.0),
                                                             fMinDCAXiPVCut(0.0),
                                                             fMaxDCAXiPVCut(5.0),
                                                             fMinDCALambdaProtonPVCut(0.0),
                                                             fMinDCALambdaPionPVCut(0.0),
                                                             fMinDCABachelorPionPVCut(0.0),
                                                             fMaxDCALambdaProtonPVCut(10.0),
                                                             fMaxDCALambdaPionPVCut(60.0),
                                                             fMaxDCABachelorPionPVCut(40),
                                                             fV0CPACut(0.995),
                                                             fXiCPACut(0.99),
                                                             fXiEtaCut(0.8),
                                                             fLambdaLowRadius(0.0),
                                                             fLambdaHighRadius(999.0),
                                                             fXiLowRadius(0.0),
                                                             fXiHighRadius(100.0),
                                                             fV0MassWindowCut(0.01),
                                                             fXiMassWindowCut(0.01),
                                                             fExoticMaxOpenAngle(0.0785398),
                                                             fXi1530YCutHigh(0.5),
                                                             fXi1530YCutLow(-0.5),
                                                             fXi1530MassHigh(1.6),
                                                             fXi1530MassLow(1.45),
                                                             fGoodTrackArray(),
                                                             fGoodXiArray()
{
    /// Default constructor
}

//_____________________________________________________________________________
AliAnalysisTaskXi1530PbPb::AliAnalysisTaskXi1530PbPb(const char *name,
                                                     Bool_t MCcase) : AliAnalysisTaskSE(name),
                                                                      fEventCuts(),
                                                                      fCheckTPCGeo(kFALSE),
                                                                      fTPCActiveLengthCutDeltaY(3.0),
                                                                      fTPCActiveLengthCutDeltaZ(220.0),
                                                                      fRequireCutGeoNcrNclLength(130),
                                                                      fRequireCutGeoNcrNclGeom1Pt(1.5),
                                                                      fCutGeoNcrNclFractionNcr(0.85),
                                                                      fCutGeoNcrNclFractionNcl(0.7),
                                                                      fTrackCuts(nullptr),
                                                                      fPIDResponse(),
                                                                      fMixingPool(),
                                                                      fEvt(nullptr),
                                                                      fMCEvent(nullptr),
                                                                      fHistos(nullptr),
                                                                      fVertex(nullptr),
                                                                      fMCArray(nullptr),
                                                                      fIsAOD(kFALSE),
                                                                      fIsNano(kFALSE),
                                                                      fSetMixing(kFALSE),
                                                                      fUseBuiltinMixer(kFALSE),
                                                                      fFillQAPlot(kTRUE),
                                                                      fIsMC(MCcase),
                                                                      fIsPrimaryMC(kFALSE),
                                                                      fFillTree(kFALSE),
                                                                      fIsINEL(kFALSE),
                                                                      fIsHM(kFALSE),
                                                                      fLambdaCPAtoXi(kFALSE),
                                                                      fExoticFinder(kFALSE),
                                                                      fSkipFillingHisto(kFALSE),
                                                                      fEMpool(0),
                                                                      fBinCent(),
                                                                      fBinZ(),
                                                                      fPosPV(),
                                                                      fMagField(0),
                                                                      fCent(-1),
                                                                      fnMix(10),
                                                                      fCentBin(-1),
                                                                      fZbin(-1),
                                                                      fFilterBit(32.0),
                                                                      fTPCNsigXi1530PionCut(3.0),
                                                                      fXi1530PionEtaCut(0.8),
                                                                      fXi1530PionPVzCut(2.5),
                                                                      fXi1530PionPVrSigmaCut(7.0),
                                                                      fTPCNsigLambdaProtonCut(3.0),
                                                                      fTPCNsigLambdaPionCut(3.0),
                                                                      fTPCNsigBachelorPionCut(3.0),
                                                                      fDCALambdaDaughtersCut(0.8),
                                                                      fDCAXiDaughtersCut(0.6),
                                                                      fMinDCALambdaPVCut(0.00),
                                                                      fMaxDCALambdaPVCut(10),
                                                                      fMinDCAXiPVCut(0.0),
                                                                      fMaxDCAXiPVCut(5.0),
                                                                      fMinDCALambdaProtonPVCut(0.0),
                                                                      fMinDCALambdaPionPVCut(0.0),
                                                                      fMinDCABachelorPionPVCut(0.0),
                                                                      fMaxDCALambdaProtonPVCut(10.0),
                                                                      fMaxDCALambdaPionPVCut(60.0),
                                                                      fMaxDCABachelorPionPVCut(40),
                                                                      fV0CPACut(0.995),
                                                                      fXiCPACut(0.99),
                                                                      fXiEtaCut(0.8),
                                                                      fLambdaLowRadius(0.0),
                                                                      fLambdaHighRadius(999.0),
                                                                      fXiLowRadius(0.0),
                                                                      fXiHighRadius(100.0),
                                                                      fV0MassWindowCut(0.01),
                                                                      fXiMassWindowCut(0.01),
                                                                      fExoticMaxOpenAngle(0.0785398),
                                                                      fXi1530YCutHigh(0.5),
                                                                      fXi1530YCutLow(-0.5),
                                                                      fXi1530MassHigh(1.6),
                                                                      fXi1530MassLow(1.45),
                                                                      fGoodTrackArray(),
                                                                      fGoodXiArray()
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TNtupleD::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskXi1530PbPb::~AliAnalysisTaskXi1530PbPb()
{
    if (AliAnalysisManager::GetAnalysisManager()->IsProofMode())
        return;
    if (fHistos)
        delete fHistos;
    if (fTree)
        delete fTree;
    if (!fIsMC)
    {
        delete fTreeXi1530Rec;
    }
}

//___________________________________________________________________
void AliAnalysisTaskXi1530PbPb::SetCutOpen()
{
    // Use for Tree study
    // Pion cuts
    SetFilterbitXi1530Pion(32);
    SetMaxNsigXi1530Pion(5);
    SetMaxEtaXi1530Pion(0.8);
    SetMaxVertexZXi1530Pion(99);
    SetMaxPVrSigmaXi1530Pion(99);

    // Lambda cuts
    SetMaxNsigV0Proton(5);
    SetMaxNsigV0Pion(5);
    SetMaxNsigBachelorPion(5);
    SetMaxDCAV0daughters(999);
    SetMaxDCAXidaughters(999);
    SetMinDCAPVV0(0);
    SetMinDCAPVV0Proton(0);
    SetMinDCAPVV0Pion(0);
    SetMinDCAPVBachelorPion(0);
    SetMinDCAPVXi(0);
    SetMinCPAV0(0.9);
    SetMinCPAXi(0.9);
    SetMaxEtaXi(0.8);
    SetLowRadiusV0(0);
    SetHighRadiusV0(200);
    SetLowRadiusXi(0);
    SetHighRadiusXi(200);
    SetMaxMassWindowV0(0.5);
    SetMaxMassWindowXi(0.5);

    // Xi1530 cut
    SetXi1530RapidityCutHigh(1);
    SetXi1530RapidityCutLow(-1);
}

//_____________________________________________________________________________
void AliAnalysisTaskXi1530PbPb::UserCreateOutputObjects()
{
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

    fHistos = new THistManager("Xi1530hists");
    auto binAnti = AxisStr("AType",
                           {"Normal",
                            "Anti"});
    auto binType = (!fIsMC) ? AxisStr("Type",
                                      {"XiStar",
                                       "XiStar_LS",
                                       "XiStar_mix"})
                            : AxisStr("Type",
                                      {"XiStar",
                                       "XiStar_LS",
                                       "XiStar_mix",
                                       "XiStar_no"});
    auto binTypeMC = AxisStr("Type",
                             {"XiStar_gen",
                              "XiStar_gen_inel10",
                              "XiStar_gen_inel10_igz",
                              "XiStar_gen_trig",
                              "XiStar_rec"});

    std::vector<double> centaxisbin;
    (fIsHM) ? centaxisbin = {0, 0.001, 0.01, 0.05, 0.1}
            : centaxisbin = {-1, 0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    // can be use from pp to PbPb
    if (fIsMC)
        centaxisbin = {-1, 0, 0.001, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    // for general MC
    fBinCent = AxisVar("Cent", centaxisbin);
    auto binPt = AxisFix("Pt", 200, 0, 20);
    auto binMass = AxisFix("Mass", 1800, 1.2, 3.0);
    // fBinZ = AxisFix("Z", 20, -10, 10);	// 1cm diff
    fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10}); // moderate diff

    CreateTHnSparse("Xi1530_data", "Xi1530_data", 5,
                    {binAnti, binType, fBinCent, binPt, binMass}, "s");
    if (fIsMC)
    {
        auto binTypeMCNorm = AxisStr("Type",
                                     {"kAll",
                                      "kINEL10",
                                      "kINEL_trig",
                                      "kINEL_trig_vtx",
                                      "kINEL_trig_vtx10",
                                      "kINELg0",
                                      "kINELg010",
                                      "kINELg0_trig",
                                      "kINELg0_trig_vtx",
                                      "kINELg0_trig_vtx10",
                                      "kSelected"});
        CreateTHnSparse("Xi1530_mc", "Xi1530_mc", 5,
                        {binAnti, binTypeMC, fBinCent, binPt, binMass}, "s");
        CreateTHnSparse("Normalisation", "", 2, {binTypeMCNorm, fBinCent}, "s");
    }

    fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
    if (fIsHM)
        fHistos->CreateTH1("hMultiplicity", "", 100, 0, 0.1, "s");
    else
        fHistos->CreateTH1("hMultiplicity", "", 101, -1, 100, "s");

    fHistos->CreateTH1("hXi1530Data", "", 1800, 1.2, 3.0, "s");
    fHistos->CreateTH1("hXi1530LS", "", 1800, 1.2, 3.0, "s");
    fHistos->CreateTH1("hXi1530Mix", "", 1800, 1.2, 3.0, "s");
    if (fFillQAPlot)
    {
        // BEFORE CUT 	///////////////////////////////////
        // TPC PID dEdX
        fHistos->CreateTH2("QA/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QA/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QA/hTPCPIDBachelorPion", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QA/hTPCPIDXi1530Pion", "", 200, 0, 20, 200, 0, 200);
        // TPC PID Signal
        fHistos->CreateTH1("QA/hTPCPIDsignalLambdaProton", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QA/hTPCPIDsignalLambdaPion", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QA/hTPCPIDsignalBachelorPion", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QA/hTPCPIDsignalXi1530Pion", "", 100, -5, 5, "s");
        // TPC Geometrical Cut
        if (fCheckTPCGeo)
            fHistos->CreateTH1("QA/hTPCGeoCheck", "", 2, -0.5, 1.5, "s");
        // DCA between daughters
        fHistos->CreateTH1("QA/hDCALambdaBTWDaughters", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QA/hDCAXiBTWDaughters", "", 300, 0, 3, "s");
        // DCA to PV
        fHistos->CreateTH1("QA/hDCArPVXi1530Pion", "", 50, 0, 0.5, "s");
        fHistos->CreateTH1("QA/hDCAXi1530PiontoPV", "", 30, 0, 3, "s");
        fHistos->CreateTH1("QA/hDCALambdatoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QA/hDCAXitoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QA/hDCALambdaProtontoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QA/hDCALambdaPiontoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QA/hDCABachelorPiontoPV", "", 50, 0, 5, "s");
        // C P A
        fHistos->CreateTH1("QA/hCPALambda", "", 500, 0.95, 1.0, "s");
        fHistos->CreateTH1("QA/hCPAXi", "", 500, 0.95, 1.0, "s");
        // E t a
        fHistos->CreateTH1("QA/hEtaXi1530Pion", "", 20, -1.0, 1.0);
        fHistos->CreateTH2("QA/hPhiEtaXi", "", 180, 0, 2 * pi, 40, -2, 2);
        // ETC
        fHistos->CreateTH1("QA/hPtXi1530Pion", "", 200, 0, 20);
        fHistos->CreateTH1("QA/hMassV0", "", 80, 1.08, 1.16, "s");
        fHistos->CreateTH1("QA/hMassXi", "", 100, 1.26, 1.36, "s");
        fHistos->CreateTH1("QA/hLambdaRXY", "", 200, 0, 200);
        fHistos->CreateTH1("QA/hXiRXY", "", 200, 0, 200);

        // AFTER CUT	///////////////////////////////////
        // TPC PID dEdX
        fHistos->CreateTH2("QAcut/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QAcut/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QAcut/hTPCPIDBachelorPion", "", 200, 0, 20, 200, 0, 200);
        fHistos->CreateTH2("QAcut/hTPCPIDXi1530Pion", "", 200, 0, 20, 200, 0, 200);
        // TPC PID Signal
        fHistos->CreateTH1("QAcut/hTPCPIDsignalLambdaProton", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QAcut/hTPCPIDsignalLambdaPion", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QAcut/hTPCPIDsignalBachelorPion", "", 100, -5, 5, "s");
        fHistos->CreateTH1("QAcut/hTPCPIDsignalXi1530Pion", "", 100, -5, 5, "s");
        // TPC Geometrical Cut
        if (fCheckTPCGeo)
            fHistos->CreateTH1("QAcut/hTPCGeoCheck", "", 2, -0.5, 1.5, "s");
        // DCA between daughters
        fHistos->CreateTH1("QAcut/hDCALambdaBTWDaughters", "", 300, 0, 3, "s");
        fHistos->CreateTH1("QAcut/hDCAXiBTWDaughters", "", 300, 0, 3, "s");
        // DCA to PV
        fHistos->CreateTH1("QAcut/hDCArPVXi1530Pion", "", 50, 0, 0.5, "s");
        fHistos->CreateTH1("QAcut/hDCAXi1530PiontoPV", "", 30, 0, 3, "s");
        fHistos->CreateTH1("QAcut/hDCALambdatoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QAcut/hDCAXitoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QAcut/hDCALambdaProtontoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QAcut/hDCALambdaPiontoPV", "", 50, 0, 5, "s");
        fHistos->CreateTH1("QAcut/hDCABachelorPiontoPV", "", 50, 0, 5, "s");
        // C P A
        fHistos->CreateTH1("QAcut/hCPALambda", "", 150, 0.85, 1.0, "s");
        fHistos->CreateTH1("QAcut/hCPAXi", "", 150, 0.85, 1.0, "s");
        // E t a
        fHistos->CreateTH1("QAcut/hEtaXi1530Pion", "", 20, -1.0, 1.0);
        fHistos->CreateTH2("QAcut/hPhiEtaXi", "", 180, 0, 2 * pi, 40, -2, 2);
        // ETC
        fHistos->CreateTH1("QAcut/hPtXi1530Pion", "", 200, 0, 20);
        fHistos->CreateTH1("QAcut/hMassV0", "", 80, 1.08, 1.16, "s");
        fHistos->CreateTH1("QAcut/hMassXi", "", 100, 1.26, 1.36, "s");
        fHistos->CreateTH1("QAcut/hLambdaRXY", "", 200, 0, 200);
        fHistos->CreateTH1("QAcut/hXiRXY", "", 200, 0, 200);
        if (fIsMC)
        {
            fHistos->CreateTH1("MCQA/hMCGeneratedXi1530Pt", "", 150, 0, 15, "s");
            fHistos->CreateTH1("MCQA/hMCReconstructedXi1530Pt", "", 150, 0, 15, "s");
        }
    }

    if (fUseBuiltinMixer)
        fEMpool.resize(fBinCent.GetNbins() + 1,
                       std::vector<eventpool>(fBinZ.GetNbins() + 1));

    fTreeXi1530Rec = fIsMC ? &fTreeXi1530Gen : new StructXi1530PbPb;
    OpenFile(2);
    fTree = new TTree("Xi1530Tree", "Xi1530Tree");
    if (fIsMC)
    {
        fTree->Branch("Xi1530TreeMC", &fTreeXi1530Gen);
    }
    else
    {
        fTree->Branch("StructXi1530PbPb", fTreeXi1530Rec);
    }

    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fTree);
}

//_____________________________________________________________________________
void AliAnalysisTaskXi1530PbPb::UserExec(Option_t *)
{
    AliVEvent *event = InputEvent();
    if (!event)
    {
        PostData(1, fHistos->GetListOfHistograms());
        PostData(2, fTree);
        AliInfo("Could not retrieve event");
        return;
    }

    AliNanoAODHeader *nanoHeader =
        dynamic_cast<AliNanoAODHeader *>(fInputEvent->GetHeader());

    event->IsA() == AliESDEvent::Class()
        ? fEvt = dynamic_cast<AliESDEvent *>(event)
        : fEvt = dynamic_cast<AliAODEvent *>(event);
    if (!fIsAOD && (event->IsA() != AliESDEvent::Class()))
        fIsAOD = true;
    if (!fEvt)
    {
        PostData(1, fHistos->GetListOfHistograms());
        PostData(2, fTree);
        return;
    }

    AliInputEventHandler *inputHandler =
        (AliInputEventHandler *)AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler();

    bool IsEvtSelected{false},
        IsINEL0True{false},
        IsSelectedTrig{false},
        IsVtxInZCut{false},
        IsGoodVertex{false};
    if (!nanoHeader)
    {
        IsEvtSelected = fEventCuts.AcceptEvent(event);
        if (fIsMC)
        {
            if (fIsAOD)
                fMCArray =
                    (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = fEventCuts.IsTrueINELgtZero(fEvt, true);
        }

        // fCent = fEventCuts.GetCentrality(0);
        fCent = AliMultSelectionTask::IsINELgtZERO(event)
                    ? fEventCuts.GetCentrality()
                    : -0.5;
        fPIDResponse = (AliPIDResponse *)inputHandler->GetPIDResponse();
        if (!fPIDResponse)
            AliInfo("No PIDd");
        IsVtxInZCut = fEventCuts.PassedCut(AliEventCuts::kVertexPosition);
        IsSelectedTrig = fEventCuts.PassedCut(AliEventCuts::kTrigger);
        IsGoodVertex = fEventCuts.PassedCut(AliEventCuts::kVertexQuality);
    }
    else
    {
        if (!fIsNano)
            fIsNano = kTRUE;
        if (fIsMC)
        {
            if (fIsAOD)
                fMCArray =
                    (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
            fMCEvent = MCEvent();
            IsINEL0True = true;
        }

        IsEvtSelected = true;
        fCent = nanoHeader->GetCentr("V0M");
        static int inel_index = -1;
        if (inel_index < 0)
            inel_index = nanoHeader->GetVarIndex("cstINELgt0");
        if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
            fCent = -0.5;
    }

    if (fIsMC)
    {
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
        PostData(2, fTree);
        return; // event cut
    }

    fHistos->FillTH1("hMultiplicity", (double)fCent);

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

    bool checkPion = GoodTracksSelection();
    bool checkCascade = GoodCascadeSelection();

    if (checkPion && checkCascade)
    {
        if (!fSkipFillingHisto)
            FillTracks(); // Fill the histogram
        if (fFillTree)
            FillTree();
    }

    if (fUseBuiltinMixer && fSetMixing && fGoodTrackArray.size())
    {
        FillTrackToEventPool(); // use only pion track pool.
    }

    PostData(1, fHistos->GetListOfHistograms());
    PostData(2, fTree);
}

//_____________________________________________________________________________
void AliAnalysisTaskXi1530PbPb::Terminate(Option_t *) {}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskXi1530PbPb::GoodTracksSelection()
{
    const UInt_t nTracks = fEvt->GetNumberOfTracks();
    fGoodTrackArray.clear();
    AliVTrack *track;
    Float_t b[2];
    Float_t bCov[3];
    Double_t nTPCNSigPion, pionZ, pionPt, pionSigmaDCA_r, pionDCA_r, lEta;
    Int_t isTPCGeo = 0;

    for (UInt_t it = 0; it < nTracks; it++)
    {
        track = (AliVTrack *)fEvt->GetTrack(it);
        if (!track)
            continue;

        // ---------- Track selection begin ----------
        if (!fIsAOD)
        {
            if (!fTrackCuts->AcceptTrack((AliESDtrack *)track))
                continue;
            if (fCheckTPCGeo)
                isTPCGeo = IsSelectedTPCGeoCut(((AliESDtrack *)track)) ? 1 : 0;
        } // ESD Case
        else
        {
            if (!fIsNano)
            {
                if (!((AliAODTrack *)track)->TestFilterBit(fFilterBit))
                    continue;
                if (fCheckTPCGeo)
                    isTPCGeo = IsSelectedTPCGeoCut(((AliAODTrack *)track)) ? 1 : 0;
            }
            else
            {
                if (!(static_cast<AliNanoAODTrack *>(track)->TestFilterBit(fFilterBit)))
                    continue;
                if (fCheckTPCGeo)
                {
                    static const Int_t tpcGeo_index = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCGeoLength");
                    isTPCGeo = (static_cast<AliNanoAODTrack *>(track)->GetVar(tpcGeo_index) > 0.5) ? 1 : 0;
                }
            }
        }

        GetImpactParam(track, b, bCov);

        pionZ = TMath::Abs(b[1]);
        nTPCNSigPion = GetTPCnSigma(track, AliPID::kPion);
        pionPt = track->Pt();
        pionSigmaDCA_r = (0.0026 + 0.0050 / pionPt);
        pionDCA_r = TMath::Abs(b[0]);
        lEta = track->Eta();
        if (fFillQAPlot)
        {
            fHistos->FillTH1("QA/hDCAXi1530PiontoPV", pionZ);
            fHistos->FillTH1("QA/hDCArPVXi1530Pion", pionDCA_r);
            fHistos->FillTH1("QA/hEtaXi1530Pion", lEta);
            fHistos->FillTH1("QA/hPtXi1530Pion", pionPt);
            fHistos->FillTH1("QA/hTPCPIDsignalXi1530Pion", nTPCNSigPion);
            fHistos->FillTH2("QA/hTPCPIDXi1530Pion", track->GetTPCmomentum(),
                             track->GetTPCsignal());
            if (fCheckTPCGeo)
                fHistos->FillTH1("QA/hTPCGeoCheck", isTPCGeo);
        }

        if (TMath::Abs(nTPCNSigPion) > fTPCNsigXi1530PionCut)
            continue;
        if (TMath::Abs(lEta) > fXi1530PionEtaCut)
            continue;
        if (pionPt < 0.15)
            continue;
        if (pionZ > fXi1530PionPVzCut)
            continue;
        if (pionDCA_r > pionSigmaDCA_r * fXi1530PionPVrSigmaCut)
            continue;
        if (fCheckTPCGeo && isTPCGeo < 0.5)
            continue;

        if (fFillQAPlot)
        {
            fHistos->FillTH1("QAcut/hDCAXi1530PiontoPV", pionZ);
            fHistos->FillTH1("QAcut/hDCArPVXi1530Pion", pionDCA_r);
            fHistos->FillTH1("QAcut/hEtaXi1530Pion", lEta);
            fHistos->FillTH1("QAcut/hPtXi1530Pion", pionPt);
            fHistos->FillTH1("QAcut/hTPCPIDsignalXi1530Pion", nTPCNSigPion);
            fHistos->FillTH2("QAcut/hTPCPIDXi1530Pion", track->GetTPCmomentum(),
                             track->GetTPCsignal());
            if (fCheckTPCGeo)
                fHistos->FillTH1("QAcut/hTPCGeoCheck", isTPCGeo);
        }

        fGoodTrackArray.push_back(it);
    }

    return fGoodTrackArray.size();
}

Bool_t AliAnalysisTaskXi1530PbPb::GoodCascadeSelection()
{
    fGoodXiArray.clear();
    const UInt_t nXi = fEvt->GetNumberOfCascades();

    AliESDcascade *xiESD;
    AliAODcascade *xiAOD;
    Double_t lTPCNSigProton(999), lTPCNSigLambdaPion(999), lTPCNSigBachelorPion(999),
        lDCALambdaProtonPV(999), lDCALambdaPionPV(999), lDCABachelorPionPV(999),
        lDCALambdaDaughters(999), lDCAXiDaughters(999), lDCALambdaPV(999),
        lDCAXiPV(999), lLambdaCPA(999), lXiCPA(999), lMassV0(999), lMassXi(999),
        lXimomsum(999), lXieta(999), lXiphi(999), radiusV0(999), radiusXi(999);
    Double_t lPosXi[3];
    Double_t lPosV0[3];
    Double_t lMomXi[3];
    Float_t b[2];    // Float due to the function input
    Float_t bCov[3]; // Float due to the function input

    Bool_t AcceptedXi = kTRUE;
    if (!fIsAOD)
    {
        // ESD case
        for (UInt_t it = 0; it < nXi; it++)
        {
            AcceptedXi = kTRUE;
            xiESD = ((AliESDEvent *)fEvt)->GetCascade(it);
            if (!xiESD)
            {
                AliInfo("No Cascade!");
                continue;
            }

            if (TMath::Abs(xiESD->GetPindex()) ==
                TMath::Abs(xiESD->GetNindex()))
            {
                AliInfo("Same index p-n");
                continue;
            }

            if (TMath::Abs(xiESD->GetPindex()) ==
                TMath::Abs(xiESD->GetBindex()))
            {
                AliInfo("Same index p-b");
                continue;
            }

            if (TMath::Abs(xiESD->GetNindex()) ==
                TMath::Abs(xiESD->GetBindex()))
            {
                AliInfo("Same index n-b");
                continue;
            }

            AliESDtrack *pTrackXi = (AliESDtrack *)fEvt->GetTrack(TMath::Abs(xiESD->GetPindex()));
            AliESDtrack *nTrackXi = (AliESDtrack *)fEvt->GetTrack(TMath::Abs(xiESD->GetNindex()));
            AliESDtrack *bTrackXi = (AliESDtrack *)fEvt->GetTrack(TMath::Abs(xiESD->GetBindex()));

            // PID cuts for Xi daughters
            if (xiESD->Charge() == -1)
            {
                // Xi- has +proton, -pion
                lTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
                if (fFillQAPlot)
                {
                    fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                    fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                }
            }
            else
            {
                // Xi+ has -proton, +pion
                lTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
                if (fFillQAPlot)
                {
                    fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                    fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                }
            }

            lTPCNSigBachelorPion = GetTPCnSigma(bTrackXi, AliPID::kPion);
            if (fFillQAPlot)
            {
                fHistos->FillTH2("QA/hTPCPIDBachelorPion",
                                 bTrackXi->GetTPCmomentum(),
                                 bTrackXi->GetTPCsignal());
                fHistos->FillTH1("QA/hTPCPIDsignalLambdaProton", lTPCNSigProton);
                fHistos->FillTH1("QA/hTPCPIDsignalLambdaPion", lTPCNSigLambdaPion);
                fHistos->FillTH1("QA/hTPCPIDsignalBachelorPion", lTPCNSigBachelorPion);
            }

            if (TMath::Abs(lTPCNSigProton) > fTPCNsigLambdaProtonCut)
            {
                AliInfo(Form("Lambda proton PID cut failed -  value: %f, cut: %f",
                             lTPCNSigProton, fTPCNsigLambdaProtonCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lTPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
            {
                AliInfo(Form("Lambda pion PID cut failed -  value: %f, cut: %f",
                             lTPCNSigLambdaPion, fTPCNsigLambdaPionCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lTPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
            {
                AliInfo(Form("Bachelor Pion PID cut failed -  value: %f, cut: %f",
                             lTPCNSigBachelorPion,
                             fTPCNsigBachelorPionCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            // DCA cut
            // DCA between Dautgher particles
            lDCALambdaDaughters = TMath::Abs(xiESD->GetDcaV0Daughters());
            lDCAXiDaughters = TMath::Abs(xiESD->GetDcaXiDaughters());
            if (lDCALambdaDaughters > fDCALambdaDaughtersCut)
            {
                AliInfo(Form(
                    "DCA Lambda daughters cut failed -  value: %f, cut: %f",
                    lDCALambdaDaughters, fDCALambdaDaughtersCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA proton-pion
            }

            if (lDCAXiDaughters > fDCAXiDaughtersCut)
            {
                AliInfo(
                    Form("DCA Xi daughters cut failed -  value: %f, cut: %f",
                         lDCAXiDaughters, fDCAXiDaughtersCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA Lambda-pion
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hDCALambdaBTWDaughters", lDCALambdaDaughters);
                fHistos->FillTH1("QA/hDCAXiBTWDaughters", lDCAXiDaughters);
            }

            // DCA to PV cut
            lDCALambdaPV =
                TMath::Abs(xiESD->GetD(fPosPV[0], fPosPV[1], fPosPV[2]));
            lDCAXiPV =
                TMath::Abs(xiESD->GetDcascade(fPosPV[0], fPosPV[1], fPosPV[2]));
            if (xiESD->Charge() == -1)
            {
                // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }
            else
            {
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }

            GetImpactParam(bTrackXi, b, bCov);
            lDCABachelorPionPV = TMath::Abs(b[0]);

            if (lDCALambdaPV < fMinDCALambdaPVCut)
            {
                AliInfo(Form("DCA Lambda PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPV, fMinDCALambdaPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Lambda-vertex
            }
            if (lDCALambdaPV > fMaxDCALambdaPVCut)
            {
                AliInfo(Form("max DCA Lambda PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPV, fMaxDCALambdaPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Lambda-vertex
            }

            if (lDCAXiPV < fMinDCAXiPVCut)
            {
                AliInfo(Form("min DCA Xi PV cut failed -  value: %f, cut: %f",
                             lDCAXiPV, fMinDCAXiPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA Xi-vertex
            }
            if (lDCAXiPV > fMaxDCAXiPVCut)
            {
                AliInfo(Form("max DCA Xi PV cut failed -  value: %f, cut: %f",
                             lDCAXiPV, fMaxDCAXiPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA Xi-vertex
            }

            if (lDCALambdaProtonPV < fMinDCALambdaProtonPVCut)
            {
                AliInfo(Form("min DCA LambdaProton PV cut failed -  value: %f, cut: %f",
                             lDCALambdaProtonPV, fMinDCALambdaProtonPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Proton of Lambda-vertex
            }

            if (lDCALambdaProtonPV > fMaxDCALambdaProtonPVCut)
            {
                AliInfo(Form("max DCA LambdaProton PV cut failed -  value: %f, cut: %f",
                             lDCALambdaProtonPV, fMaxDCALambdaProtonPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Proton of Lambda-vertex
            }

            if (lDCALambdaPionPV < fMinDCALambdaPionPVCut)
            {
                AliInfo(Form("min DCA LambdaPion PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPionPV, fMinDCALambdaPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Pion of Lambda-vertex
            }

            if (lDCALambdaPionPV > fMaxDCALambdaPionPVCut)
            {
                AliInfo(Form("max DCA LambdaPion PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPionPV, fMaxDCALambdaPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Pion of Lambda-vertex
            }

            if (lDCABachelorPionPV < fMinDCABachelorPionPVCut)
            {
                AliInfo(Form("min DCA BachelorPion PV cut failed -  value: %f, cut: %f",
                             lDCABachelorPionPV, fMinDCABachelorPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA BachelorPion of Xi-vertex
            }

            if (lDCABachelorPionPV > fMaxDCABachelorPionPVCut)
            {
                AliInfo(Form("max DCA BachelorPion PV cut failed -  value: %f, cut: %f",
                             lDCABachelorPionPV, fMaxDCABachelorPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA BachelorPion of Xi-vertex
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hDCALambdatoPV", lDCALambdaPV);
                fHistos->FillTH1("QA/hDCAXitoPV", lDCAXiPV);
                fHistos->FillTH1("QA/hDCALambdaProtontoPV", lDCALambdaProtonPV);
                fHistos->FillTH1("QA/hDCALambdaPiontoPV", lDCALambdaPionPV);
                fHistos->FillTH1("QA/hDCABachelorPiontoPV", lDCABachelorPionPV);
            }

            // CPA cut
            lXiCPA = TMath::Abs(xiESD->GetCascadeCosineOfPointingAngle(fPosPV[0], fPosPV[1], fPosPV[2]));
            xiESD->GetXYZcascade(lPosXi[0], lPosXi[1], lPosXi[2]);
            if (fLambdaCPAtoXi)
                lLambdaCPA = TMath::Abs(xiESD->GetV0CosineOfPointingAngle(lPosXi[0], lPosXi[1], lPosXi[2]));
            else
                lLambdaCPA = TMath::Abs(xiESD->GetV0CosineOfPointingAngle(fPosPV[0], fPosPV[1], fPosPV[2]));

            if (lLambdaCPA < fV0CPACut)
            {
                AliInfo(Form("CPA Lambda cut failed -  value: %f, cut: %f",
                             lLambdaCPA, fV0CPACut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (lXiCPA < fXiCPACut)
            {
                AliInfo(Form("CPA Xi cut failed -  value: %f, cut: %f", lXiCPA,
                             fXiCPACut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hCPALambda", lLambdaCPA);
                fHistos->FillTH1("QA/hCPAXi", lXiCPA);
            }

            // Mass window cut
            lMassXi = xiESD->GetEffMassXi();
            lMassV0 = xiESD->GetEffMass();

            if (TMath::Abs(lMassV0 - v0Mass) > fV0MassWindowCut)
            {
                AliInfo(Form("Lambda Mass cut failed -  value: %f, cut: %f", lMassV0 - v0Mass, fV0MassWindowCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lMassXi - Ximass) > fXiMassWindowCut)
            {
                AliInfo(Form("Xi Mass cut failed -  value: %f, cut: %f", TMath::Abs(lMassXi - Ximass), fXiMassWindowCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hMassXi", lMassXi);
                fHistos->FillTH1("QA/hMassV0", lMassV0);
            }

            // Eta cut
            lXieta = xiESD->Eta();
            lXiphi = xiESD->Phi();
            if (TMath::Abs(lXieta) > fXiEtaCut)
            {
                AliInfo(Form("Eta cut failed -  value: %f, cut: %f", xiESD->Eta(), fXiEtaCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
                fHistos->FillTH2("QA/hPhiEtaXi", lXiphi, lXieta);

            // Radius cut
            xiESD->GetXYZ(lPosV0[0], lPosV0[1], lPosV0[2]);
            xiESD->GetXYZcascade(lPosXi[0], lPosXi[1], lPosXi[2]);
            radiusV0 = TMath::Hypot(lPosV0[0], lPosV0[1]);
            radiusXi = TMath::Hypot(lPosXi[0] - fPosPV[0], lPosXi[1] - fPosPV[1]);

            if ((radiusV0 < fLambdaLowRadius) || (radiusV0 > fLambdaHighRadius))
            {
                AliInfo(Form("V0 radius cut failed -  value: %f, cut: %f < value < %f", radiusV0, fLambdaLowRadius, fLambdaHighRadius));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if ((radiusXi < fLambdaLowRadius) || (radiusXi > fLambdaHighRadius))
            {
                AliInfo(Form("Xi radius cut failed -  value: %f, cut: %f < value < %f", radiusXi, fXiLowRadius, fXiHighRadius));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hLambdaRXY", radiusV0);
                fHistos->FillTH1("QA/hXiRXY", radiusXi);
            }

            // After selection above
            if (AcceptedXi)
            {
                fGoodXiArray.push_back(it); // for standard V0
                if (fFillQAPlot)
                {
                    if (xiESD->Charge() == -1)
                    {
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                         pTrackXi->GetTPCmomentum(),
                                         pTrackXi->GetTPCsignal());
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                         nTrackXi->GetTPCmomentum(),
                                         nTrackXi->GetTPCsignal());
                    }
                    else
                    {
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                         nTrackXi->GetTPCmomentum(),
                                         nTrackXi->GetTPCsignal());
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                         pTrackXi->GetTPCmomentum(),
                                         pTrackXi->GetTPCsignal());
                    }

                    fHistos->FillTH2("QAcut/hTPCPIDBachelorPion",
                                     bTrackXi->GetTPCmomentum(),
                                     bTrackXi->GetTPCsignal());
                    fHistos->FillTH1("QAcut/hTPCPIDsignalLambdaProton", lTPCNSigProton);
                    fHistos->FillTH1("QAcut/hTPCPIDsignalLambdaPion", lTPCNSigLambdaPion);
                    fHistos->FillTH1("QAcut/hTPCPIDsignalBachelorPion", lTPCNSigBachelorPion);
                    fHistos->FillTH1("QAcut/hDCALambdaBTWDaughters", lDCALambdaDaughters);
                    fHistos->FillTH1("QAcut/hDCAXiBTWDaughters", lDCAXiDaughters);
                    fHistos->FillTH1("QAcut/hDCALambdatoPV", lDCALambdaPV);
                    fHistos->FillTH1("QAcut/hDCAXitoPV", lDCAXiPV);
                    fHistos->FillTH1("QAcut/hDCALambdaProtontoPV", lDCALambdaProtonPV);
                    fHistos->FillTH1("QAcut/hDCALambdaPiontoPV", lDCALambdaPionPV);
                    fHistos->FillTH1("QAcut/hDCABachelorPiontoPV", lDCABachelorPionPV);
                    fHistos->FillTH1("QAcut/hCPALambda", lLambdaCPA);
                    fHistos->FillTH1("QAcut/hCPAXi", lXiCPA);
                    fHistos->FillTH1("QAcut/hMassXi", lMassXi);
                    fHistos->FillTH1("QAcut/hMassV0", lMassV0);
                    fHistos->FillTH2("QAcut/hPhiEtaXi", lXiphi, lXieta);
                    fHistos->FillTH1("QAcut/hLambdaRXY", radiusV0);
                    fHistos->FillTH1("QAcut/hXiRXY", radiusXi);
                }
            }
        } // All Xi loop
    }     // ESD case
    else
    {
        for (UInt_t it = 0; it < nXi; it++)
        {
            AcceptedXi = kTRUE;
            xiAOD = ((AliAODEvent *)fEvt)->GetCascade(it);
            if (!xiAOD)
            {
                AliInfo("No Cascade!");
                continue;
            }

            if (TMath::Abs(xiAOD->GetPosID()) ==
                TMath::Abs(xiAOD->GetNegID()))
            {
                AliInfo("Same index p-n");
                continue;
            }

            if (TMath::Abs(xiAOD->GetPosID()) ==
                TMath::Abs(xiAOD->GetBachID()))
            {
                AliInfo("Same index p-b");
                continue;
            }

            if (TMath::Abs(xiAOD->GetNegID()) ==
                TMath::Abs(xiAOD->GetBachID()))
            {
                AliInfo("Same index n-b");
                continue;
            }

            AliAODTrack *pTrackXi = (AliAODTrack *)(xiAOD->GetDaughter(0));
            AliAODTrack *nTrackXi = (AliAODTrack *)(xiAOD->GetDaughter(1));
            AliAODTrack *bTrackXi = (AliAODTrack *)(xiAOD->GetDecayVertexXi()->GetDaughter(0));

            // PID cuts for Xi daughters
            if (xiAOD->ChargeXi() == -1)
            {
                // Xi- has +proton, -pion
                lTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
                if (fFillQAPlot)
                {
                    fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                    fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                }
            }
            else
            {
                // Xi+ has -proton, +pion
                lTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
                if (fFillQAPlot)
                {
                    fHistos->FillTH2("QA/hTPCPIDLambdaProton",
                                     nTrackXi->GetTPCmomentum(),
                                     nTrackXi->GetTPCsignal());
                    fHistos->FillTH2("QA/hTPCPIDLambdaPion",
                                     pTrackXi->GetTPCmomentum(),
                                     pTrackXi->GetTPCsignal());
                }
            }

            lTPCNSigBachelorPion = GetTPCnSigma(bTrackXi, AliPID::kPion);
            if (fFillQAPlot)
            {
                fHistos->FillTH2("QA/hTPCPIDBachelorPion",
                                 bTrackXi->GetTPCmomentum(),
                                 bTrackXi->GetTPCsignal());
                fHistos->FillTH1("QA/hTPCPIDsignalLambdaProton", lTPCNSigProton);
                fHistos->FillTH1("QA/hTPCPIDsignalLambdaPion", lTPCNSigLambdaPion);
                fHistos->FillTH1("QA/hTPCPIDsignalBachelorPion", lTPCNSigBachelorPion);
            }

            if (TMath::Abs(lTPCNSigProton) > fTPCNsigLambdaProtonCut)
            {
                AliInfo(Form("Lambda proton PID cut failed -  value: %f, cut: %f",
                             lTPCNSigProton, fTPCNsigLambdaProtonCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lTPCNSigLambdaPion) > fTPCNsigLambdaPionCut)
            {
                AliInfo(Form("Lambda pion PID cut failed -  value: %f, cut: %f",
                             lTPCNSigLambdaPion, fTPCNsigLambdaPionCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lTPCNSigBachelorPion) > fTPCNsigBachelorPionCut)
            {
                AliInfo(Form("Bachelor Pion PID cut failed -  value: %f, cut: %f",
                             lTPCNSigBachelorPion,
                             fTPCNsigBachelorPionCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            // DCA cut
            // DCA between Dautgher particles
            lDCALambdaDaughters = TMath::Abs(xiAOD->DcaV0Daughters());
            lDCAXiDaughters = TMath::Abs(xiAOD->DcaXiDaughters());
            if (lDCALambdaDaughters > fDCALambdaDaughtersCut)
            {
                AliInfo(Form(
                    "DCA Lambda daughters cut failed -  value: %f, cut: %f",
                    lDCALambdaDaughters, fDCALambdaDaughtersCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA proton-pion
            }

            if (lDCAXiDaughters > fDCAXiDaughtersCut)
            {
                AliInfo(
                    Form("DCA Xi daughters cut failed -  value: %f, cut: %f",
                         lDCAXiDaughters, fDCAXiDaughtersCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA Lambda-pion
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hDCALambdaBTWDaughters", lDCALambdaDaughters);
                fHistos->FillTH1("QA/hDCAXiBTWDaughters", lDCAXiDaughters);
            }

            // DCA to PV cut
            lDCALambdaPV = TMath::Abs(xiAOD->DcaV0ToPrimVertex());
            lDCAXiPV = TMath::Abs(xiAOD->DcaXiToPrimVertex(fPosPV[0], fPosPV[1], fPosPV[2]));
            if (xiAOD->Charge() == -1)
            {
                // Xi- has +proton, -pion
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }
            else
            {
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }

            GetImpactParam(bTrackXi, b, bCov);
            lDCABachelorPionPV = TMath::Abs(b[0]);

            if (lDCALambdaPV < fMinDCALambdaPVCut)
            {
                AliInfo(Form("min DCA Lambda PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPV, fMinDCALambdaPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Lambda-vertex
            }

            if (lDCALambdaPV > fMaxDCALambdaPVCut)
            {
                AliInfo(Form("max DCA Lambda PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPV, fMaxDCALambdaPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Lambda-vertex
            }

            if (lDCAXiPV < fMinDCAXiPVCut)
            {
                AliInfo(Form("DCA Xi PV cut failed -  value: %f, cut: %f",
                             lDCAXiPV, fMinDCAXiPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // DCA Xi-vertex
            }

            if (lDCALambdaProtonPV < fMinDCALambdaProtonPVCut)
            {
                AliInfo(Form("min DCA LambdaProton PV cut failed -  value: %f, cut: %f",
                             lDCALambdaProtonPV, fMinDCALambdaProtonPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Proton of Lambda-vertex
            }

            if (lDCALambdaProtonPV > fMaxDCALambdaProtonPVCut)
            {
                AliInfo(Form("max DCA LambdaProton PV cut failed -  value: %f, cut: %f",
                             lDCALambdaProtonPV, fMaxDCALambdaProtonPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Proton of Lambda-vertex
            }

            if (lDCALambdaPionPV < fMinDCALambdaPionPVCut)
            {
                AliInfo(Form("min DCA LambdaPion PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPionPV, fMinDCALambdaPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA Pion of Lambda-vertex
            }

            if (lDCALambdaPionPV > fMaxDCALambdaPionPVCut)
            {
                AliInfo(Form("max DCA LambdaPion PV cut failed -  value: %f, cut: %f",
                             lDCALambdaPionPV, fMaxDCALambdaPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA Pion of Lambda-vertex
            }

            if (lDCABachelorPionPV < fMinDCABachelorPionPVCut)
            {
                AliInfo(Form("min DCA BachelorPion PV cut failed -  value: %f, cut: %f",
                             lDCABachelorPionPV, fMinDCABachelorPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // min DCA BachelorPion of Xi-vertex
            }

            if (lDCABachelorPionPV > fMaxDCABachelorPionPVCut)
            {
                AliInfo(Form("max DCA BachelorPion PV cut failed -  value: %f, cut: %f",
                             lDCABachelorPionPV, fMaxDCABachelorPionPVCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue; // max DCA BachelorPion of Xi-vertex
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hDCALambdatoPV", lDCALambdaPV);
                fHistos->FillTH1("QA/hDCAXitoPV", lDCAXiPV);
                fHistos->FillTH1("QA/hDCALambdaProtontoPV", lDCALambdaProtonPV);
                fHistos->FillTH1("QA/hDCALambdaPiontoPV", lDCALambdaPionPV);
                fHistos->FillTH1("QA/hDCABachelorPiontoPV", lDCABachelorPionPV);
            }

            // CPA cut
            lPosXi[0] = xiAOD->DecayVertexXiX();
            lPosXi[1] = xiAOD->DecayVertexXiY();
            lPosXi[2] = xiAOD->DecayVertexXiZ();
            lXiCPA = TMath::Abs(xiAOD->CosPointingAngleXi(fPosPV[0], fPosPV[1], fPosPV[2]));
            if (fLambdaCPAtoXi)
                lLambdaCPA = TMath::Abs(xiAOD->CosPointingAngle(lPosXi));
            else
                lLambdaCPA = TMath::Abs(xiAOD->CosPointingAngle(fPosPV));

            if (lLambdaCPA < fV0CPACut)
            {
                AliInfo(Form("CPA Lambda cut failed -  value: %f, cut: %f",
                             lLambdaCPA, fV0CPACut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (lXiCPA < fXiCPACut)
            {
                AliInfo(Form("CPA Xi cut failed -  value: %f, cut: %f", lXiCPA,
                             fXiCPACut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hCPALambda", lLambdaCPA);
                fHistos->FillTH1("QA/hCPAXi", lXiCPA);
            }

            // Mass window cut
            lMassXi = xiAOD->MassXi();
            lMassV0 = -999;
            if (xiAOD->ChargeXi() == -1)
                lMassV0 = xiAOD->MassLambda();
            else
                lMassV0 = xiAOD->MassAntiLambda();

            if (TMath::Abs(lMassV0 - v0Mass) > fV0MassWindowCut)
            {
                AliInfo(Form("Lambda Mass cut failed -  value: %f, cut: %f", lMassV0 - v0Mass, fV0MassWindowCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (TMath::Abs(lMassXi - Ximass) > fXiMassWindowCut)
            {
                AliInfo(Form("Xi Mass cut failed -  value: %f, cut: %f", TMath::Abs(lMassXi - Ximass), fXiMassWindowCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hMassXi", lMassXi);
                fHistos->FillTH1("QA/hMassV0", lMassV0);
            }

            // Eta cut
            lMomXi[0] = xiAOD->MomXiX();
            lMomXi[1] = xiAOD->MomXiY();
            lMomXi[2] = xiAOD->MomXiZ();
            lXimomsum =
                TMath::Sqrt(lMomXi[0] * lMomXi[0] + lMomXi[1] * lMomXi[1] +
                            lMomXi[2] * lMomXi[2]);
            lXieta =
                0.5 * TMath::Log((lXimomsum + lMomXi[2]) /
                                 (lXimomsum - lMomXi[2] + 1.e-13));
            lXiphi =
                TMath::Pi() + TMath::ATan2(-lMomXi[1], -lMomXi[0]);
            if (TMath::Abs(lXieta) > fXiEtaCut)
            {
                AliInfo(Form("Eta cut failed -  value: %f, cut: %f", xiAOD->Eta(), fXiEtaCut));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
                fHistos->FillTH2("QA/hPhiEtaXi", lXiphi, lXieta);

            // Radius cut
            lPosV0[0] = xiAOD->DecayVertexV0X();
            lPosV0[1] = xiAOD->DecayVertexV0Y();
            lPosV0[2] = xiAOD->DecayVertexV0Z();
            radiusV0 = TMath::Hypot(lPosV0[0], lPosV0[1]);
            radiusXi = TMath::Hypot(lPosXi[0] - fPosPV[0], lPosXi[1] - fPosPV[1]);

            if ((radiusV0 < fLambdaLowRadius) || (radiusV0 > fLambdaHighRadius))
            {
                AliInfo(Form("V0 radius cut failed -  value: %f, cut: %f < value < %f", radiusV0, fLambdaLowRadius, fLambdaHighRadius));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if ((radiusXi < fLambdaLowRadius) || (radiusXi > fLambdaHighRadius))
            {
                AliInfo(Form("Xi radius cut failed -  value: %f, cut: %f < value < %f", radiusXi, fXiLowRadius, fXiHighRadius));
                if (fFillQAPlot)
                    AcceptedXi = kFALSE;
                else
                    continue;
            }

            if (fFillQAPlot)
            {
                fHistos->FillTH1("QA/hLambdaRXY", radiusV0);
                fHistos->FillTH1("QA/hXiRXY", radiusXi);
            }

            // After selection above
            if (AcceptedXi)
            {
                fGoodXiArray.push_back(it); // for standard V0
                if (fFillQAPlot)
                {
                    if (xiAOD->ChargeXi() == -1)
                    {
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                         pTrackXi->GetTPCmomentum(),
                                         pTrackXi->GetTPCsignal());
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                         nTrackXi->GetTPCmomentum(),
                                         nTrackXi->GetTPCsignal());
                    }
                    else
                    {
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaProton",
                                         nTrackXi->GetTPCmomentum(),
                                         nTrackXi->GetTPCsignal());
                        fHistos->FillTH2("QAcut/hTPCPIDLambdaPion",
                                         pTrackXi->GetTPCmomentum(),
                                         pTrackXi->GetTPCsignal());
                    }

                    fHistos->FillTH2("QAcut/hTPCPIDBachelorPion",
                                     bTrackXi->GetTPCmomentum(),
                                     bTrackXi->GetTPCsignal());
                    fHistos->FillTH1("QAcut/hTPCPIDsignalLambdaProton", lTPCNSigProton);
                    fHistos->FillTH1("QAcut/hTPCPIDsignalLambdaPion", lTPCNSigLambdaPion);
                    fHistos->FillTH1("QAcut/hTPCPIDsignalBachelorPion", lTPCNSigBachelorPion);
                    fHistos->FillTH1("QAcut/hDCALambdaBTWDaughters", lDCALambdaDaughters);
                    fHistos->FillTH1("QAcut/hDCAXiBTWDaughters", lDCAXiDaughters);
                    fHistos->FillTH1("QAcut/hDCALambdatoPV", lDCALambdaPV);
                    fHistos->FillTH1("QAcut/hDCAXitoPV", lDCAXiPV);
                    fHistos->FillTH1("QAcut/hDCALambdaProtontoPV", lDCALambdaProtonPV);
                    fHistos->FillTH1("QAcut/hDCALambdaPiontoPV", lDCALambdaPionPV);
                    fHistos->FillTH1("QAcut/hDCABachelorPiontoPV", lDCABachelorPionPV);
                    fHistos->FillTH1("QAcut/hCPALambda", lLambdaCPA);
                    fHistos->FillTH1("QAcut/hCPAXi", lXiCPA);
                    fHistos->FillTH1("QAcut/hMassXi", lMassXi);
                    fHistos->FillTH1("QAcut/hMassV0", lMassV0);
                    fHistos->FillTH2("QAcut/hPhiEtaXi", lXiphi, lXieta);
                    fHistos->FillTH1("QAcut/hLambdaRXY", radiusV0);
                    fHistos->FillTH1("QAcut/hXiRXY", radiusXi);
                }
            }
        }
    } // All Xi loop
    return fGoodXiArray.size();
}

void AliAnalysisTaskXi1530PbPb::FillTracks()
{
    AliVTrack *track1;
    AliVTrack *track_mix;
    AliESDcascade *xiESD;
    AliAODcascade *xiAOD;
    Bool_t isXiAnti = kFALSE;
    Bool_t SkipMixing = kFALSE;
    Int_t pID, nID, bID;
    auto sign = kAllType;
    auto binAnti = kNormal;

    // Mixing check
    Float_t b[2];
    Float_t bCov[3];
    Double_t nTPCNSigPion, pionZ, pionSigmaDCA_r, pionDCA_r;

    TLorentzVector vecXi, vecPion, vecPionMix;
    TLorentzVector vecXi1530;
    const UInt_t nXi = fGoodXiArray.size();
    const UInt_t nTracks = fGoodTrackArray.size();

    tracklist trackpool;
    if (fSetMixing)
    {
        eventpool &ep = (!fUseBuiltinMixer)
                            ? fMixingPool->GetMixingPool()[fCentBin][fZbin]
                            : fEMpool[fCentBin][fZbin];
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

    for (UInt_t i = 0; i < nXi; i++)
    {
        if (!fIsAOD)
        {
            xiESD = ((AliESDEvent *)fEvt)->GetCascade(fGoodXiArray[i]);
            if (!xiESD)
                continue;
            AliESDtrack *pTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetPindex()));
            AliESDtrack *nTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetNindex()));
            AliESDtrack *bTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetBindex()));

            pID = pTrackXi->GetID();
            nID = nTrackXi->GetID();
            bID = nTrackXi->GetID();

            if (xiESD->Charge() > 0)
                isXiAnti = kTRUE;
            else
                isXiAnti = kFALSE;

            vecXi.SetXYZM(xiESD->Px(), xiESD->Py(), xiESD->Pz(), xiESD->GetEffMassXi());
        }
        else
        {
            xiAOD = ((AliAODEvent *)fEvt)->GetCascade(fGoodXiArray[i]);
            if (!xiAOD)
                continue;
            AliAODTrack *pTrackXi =
                (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack *nTrackXi =
                (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(1));
            AliAODTrack *bTrackXi =
                (AliAODTrack *)(xiAOD->GetDecayVertexXi()->GetDaughter(0));
            pID = pTrackXi->GetID();
            nID = nTrackXi->GetID();
            bID = bTrackXi->GetID();

            if (xiAOD->ChargeXi() > 0)
                isXiAnti = kTRUE;
            else
                isXiAnti = kFALSE;

            vecXi.SetXYZM(xiAOD->MomXiX(), xiAOD->MomXiY(), xiAOD->MomXiZ(), xiAOD->MassXi());
        }

        (isXiAnti) ? binAnti = kAnti : binAnti = kNormal;

        // pion loop start
        for (UInt_t j = 0; j < nTracks; j++)
        {
            track1 = (AliVTrack *)fEvt->GetTrack(fGoodTrackArray[j]);
            if (!track1)
                continue;

            if (track1->GetID() == pID ||
                track1->GetID() == nID ||
                track1->GetID() == bID)
            {
                AliInfo(Form("same track for Xi1530 prion!"));
                continue;
            }

            vecPion.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

            vecXi1530 = vecXi + vecPion; // vecXi = cascade, vecPion=pion
            // Y cut
            if ((vecXi1530.Rapidity() > fXi1530YCutHigh) ||
                (vecXi1530.Rapidity() < fXi1530YCutLow))
                continue;

            // Opening Angle - Not using in normal mode
            if (fExoticFinder)
            {
                Double_t angle = vecXi.Angle(vecPion.Vect());
                if (fFillQAPlot)
                    fHistos->FillTH1("QA/hExoticOpenAngle", angle);
                if (TMath::Abs(angle) < fExoticMaxOpenAngle) // 4.5 degree
                    continue;
            }

            AliInfo("All cut passed!!");
            sign = kAllType;
            if (((isXiAnti) && (track1->Charge() < 0)) ||
                ((!isXiAnti) && (track1->Charge() > 0)))
            {
                sign = kXiStar;
                fHistos->FillTH1("hXi1530Data", vecXi1530.M());
            }
            else
            {
                sign = kXiStar_LS;
                fHistos->FillTH1("hXi1530LS", vecXi1530.M());
            }

            FillTHnSparse("Xi1530_data",
                          {(double)binAnti,
                           (double)sign,
                           (double)fCent,
                           vecXi1530.Pt(),
                           vecXi1530.M()});

            if (fIsMC)
            {
                if (IsTrueXi1530(fGoodXiArray[i], fGoodTrackArray[j]))
                {
                    AliInfo("True Xi1530");
                    FillTHnSparse("Xi1530_mc",
                                  {(double)binAnti,
                                   (double)kXiStar_REC,
                                   (double)fCent,
                                   vecXi1530.Pt(),
                                   vecXi1530.M()});
                    if (fFillQAPlot)
                        fHistos->FillTH1("MCQA/hMCReconstructedXi1530Pt", vecXi1530.Pt());
                }
                else
                {
                    // AliInfo("Not True Xi1530");
                    // MC not true bkg
                    sign = kXiStar_NOT;
                    FillTHnSparse("Xi1530_data",
                                  {(double)binAnti,
                                   (double)sign,
                                   (double)fCent,
                                   vecXi1530.Pt(),
                                   vecXi1530.M()});
                }
            }
        } // pion loop end

        // event mix pion loop start
        if (fSetMixing && !SkipMixing && (fCentBin >= 0) && (fZbin >= 0))
        {
            for (UInt_t jt = 0; jt < trackpool.size(); jt++)
            {
                track_mix = trackpool.at(jt);
                if (track_mix->GetID() == pID ||
                    track_mix->GetID() == nID ||
                    track_mix->GetID() == bID)
                {
                    AliInfo(Form("same track for Xi1530 prion!"));
                    continue;
                }

                if (!fUseBuiltinMixer)
                {
                    // We need to check further.
                    GetImpactParam(track_mix, b, bCov);
                    pionZ = TMath::Abs(b[1]);
                    nTPCNSigPion = GetTPCnSigma(track_mix, AliPID::kPion);
                    pionSigmaDCA_r = (0.0026 + 0.0050 / track_mix->Pt());
                    pionDCA_r = TMath::Abs(b[0]);
                    if (TMath::Abs(nTPCNSigPion) > fTPCNsigXi1530PionCut)
                        continue;
                    if (pionZ > fXi1530PionPVzCut)
                        continue;
                    if (pionDCA_r > pionSigmaDCA_r * fXi1530PionPVrSigmaCut)
                        continue;
                }

                vecPionMix.SetXYZM(track_mix->Px(), track_mix->Py(), track_mix->Pz(), pionMass);
                vecXi1530 = vecXi + vecPionMix;

                // Y cut
                if ((vecXi1530.Rapidity() > fXi1530YCutHigh) ||
                    (vecXi1530.Rapidity() < fXi1530YCutLow))
                    continue;

                // Opening Angle - Not using in normal mode
                if (fExoticFinder)
                {
                    Double_t angle = vecXi.Angle(vecPionMix.Vect());
                    fHistos->FillTH1("QA/hExoticOpenAngle", angle);
                    if (TMath::Abs(angle) < fExoticMaxOpenAngle) // 4.5 degree
                        continue;
                }

                sign = kAllType;
                if (((isXiAnti) && (track_mix->Charge() < 0)) ||
                    ((!isXiAnti) && (track_mix->Charge() > 0)))
                {
                    sign = kXiStar_MIX;

                    fHistos->FillTH1("hXi1530Mix", vecXi1530.M());

                    FillTHnSparse("Xi1530_data",
                                  {(double)binAnti,
                                   (double)sign,
                                   (double)fCent,
                                   vecXi1530.Pt(),
                                   vecXi1530.M()});
                } // Don't conosider LS for mix
            }
        } // event mix pion loop end
    }
}

void AliAnalysisTaskXi1530PbPb::FillTree()
{
    // FillTree will fill the ntuples for Tree Analysis
    // It can be turned on by using SetFillnTuple() Default: Off

    AliVTrack *track1 = nullptr;
    AliESDcascade *xiESD = nullptr;
    AliAODcascade *xiAOD = nullptr;
    Double_t lTPCNSigProton(999), lTPCNSigLambdaPion(999), lTPCNSigBachelorPion(999),
        lDCALambdaProtonPV(999), lDCALambdaPionPV(999), lDCABachelorPionPV(999),
        lDCALambdaDaughters(999), lDCAXiDaughters(999), lDCALambdaPV(999),
        lDCAXiPV(999), lLambdaCPA(999), lXiCPA(999), lMassV0(999), lMassXi(999),
        lXimomsum(999), lXieta(999), radiusV0(999), radiusXi(999);
    Double_t lPosXi[3], lPosV0[3], lMomXi[3], pos[3], cov[6], lPosRsnVtx[3], dztemp[2], covartemp[3], chi2perNDF, dispersion;
    Double_t cv[21], xdummy, ydummy, dca, radiusRsn;
    Float_t b[2];    // Float due to the function input
    Float_t bCov[3]; // Float due to the function input

    Int_t pID, nID, bID;
    Bool_t isXiAnti = kFALSE;

    TLorentzVector vecXi, vecPion;
    TLorentzVector vecXi1530;
    const UInt_t nXi = fGoodXiArray.size();
    const UInt_t nTracks = fGoodTrackArray.size();
    for (UInt_t i = 0; i < nXi; i++)
    {
        if (!fIsAOD)
        {
            xiESD = ((AliESDEvent *)fEvt)->GetCascade(fGoodXiArray[i]);
            if (!xiESD)
                continue;
            AliESDtrack *pTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetPindex()));
            AliESDtrack *nTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetNindex()));
            AliESDtrack *bTrackXi =
                ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetBindex()));

            pID = pTrackXi->GetID();
            nID = nTrackXi->GetID();
            bID = nTrackXi->GetID();

            if (xiESD->Charge() > 0)
                isXiAnti = kTRUE;
            else
                isXiAnti = kFALSE;

            vecXi.SetXYZM(xiESD->Px(), xiESD->Py(), xiESD->Pz(), xiESD->GetEffMassXi());

            if (!isXiAnti)
            {
                lTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            }
            else
            {
                lTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }

            lTPCNSigBachelorPion = GetTPCnSigma(bTrackXi, AliPID::kPion);

            lDCALambdaDaughters = TMath::Abs(xiESD->GetDcaV0Daughters());
            lDCAXiDaughters = TMath::Abs(xiESD->GetDcaXiDaughters());

            lDCALambdaPV =
                TMath::Abs(xiESD->GetD(fPosPV[0], fPosPV[1], fPosPV[2]));
            lDCAXiPV =
                TMath::Abs(xiESD->GetDcascade(fPosPV[0], fPosPV[1], fPosPV[2]));

            if (!isXiAnti)
            {
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }
            else
            {
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }

            GetImpactParam(bTrackXi, b, bCov);
            lDCABachelorPionPV = TMath::Abs(b[0]);

            lXiCPA = TMath::Abs(xiESD->GetCascadeCosineOfPointingAngle(fPosPV[0], fPosPV[1], fPosPV[2]));
            xiESD->GetXYZcascade(lPosXi[0], lPosXi[1], lPosXi[2]);
            if (fLambdaCPAtoXi)
                lLambdaCPA = TMath::Abs(xiESD->GetV0CosineOfPointingAngle(lPosXi[0], lPosXi[1], lPosXi[2]));
            else
                lLambdaCPA = TMath::Abs(xiESD->GetV0CosineOfPointingAngle(fPosPV[0], fPosPV[1], fPosPV[2]));

            lXieta = xiESD->Eta();

            xiESD->GetXYZ(lPosV0[0], lPosV0[1], lPosV0[2]);
            xiESD->GetXYZcascade(lPosXi[0], lPosXi[1], lPosXi[2]);
            radiusV0 = TMath::Hypot(lPosV0[0], lPosV0[1]);
            radiusXi = TMath::Hypot(lPosXi[0] - fPosPV[0], lPosXi[1] - fPosPV[1]);

            lMassXi = xiESD->GetEffMassXi();
            lMassV0 = xiESD->GetEffMass();

            // Tree
            fTreeXi1530Rec->TPCNsigLambdaProton = lTPCNSigProton;
            fTreeXi1530Rec->TPCNsigLambdaPion = lTPCNSigLambdaPion;
            fTreeXi1530Rec->TPCNsigBachelorPion = lTPCNSigBachelorPion;
            fTreeXi1530Rec->DCALambdaDaughters = lDCALambdaDaughters;
            fTreeXi1530Rec->DCAXiDaughters = lDCAXiDaughters;
            fTreeXi1530Rec->DCALambdaPV = lDCALambdaPV;
            fTreeXi1530Rec->DCALambdaProtonPV = lDCALambdaProtonPV;
            fTreeXi1530Rec->DCALambdaPionPV = lDCALambdaPionPV;
            fTreeXi1530Rec->DCABachelorPionPV = lDCABachelorPionPV;
            fTreeXi1530Rec->DCAXiPV = lDCAXiPV;
            fTreeXi1530Rec->V0CPA = lLambdaCPA;
            fTreeXi1530Rec->XiCPA = lXiCPA;
            fTreeXi1530Rec->XiEta = lXieta;
            fTreeXi1530Rec->LambdaRadius = radiusV0;
            fTreeXi1530Rec->XiRadius = radiusXi;
            fTreeXi1530Rec->V0Mass = lMassV0;
            fTreeXi1530Rec->XiMass = lMassXi;
        }
        else
        {
            xiAOD = ((AliAODEvent *)fEvt)->GetCascade(fGoodXiArray[i]);
            if (!xiAOD)
                continue;
            AliAODTrack *pTrackXi =
                (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(0));
            AliAODTrack *nTrackXi =
                (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(1));
            AliAODTrack *bTrackXi =
                (AliAODTrack *)(xiAOD->GetDecayVertexXi()->GetDaughter(0));
            pID = pTrackXi->GetID();
            nID = nTrackXi->GetID();
            bID = nTrackXi->GetID();

            if (xiAOD->ChargeXi() > 0)
                isXiAnti = kTRUE;
            else
                isXiAnti = kFALSE;

            vecXi.SetXYZM(xiAOD->MomXiX(), xiAOD->MomXiY(), xiAOD->MomXiZ(), xiAOD->MassXi());

            if (!isXiAnti)
            {
                lTPCNSigProton = GetTPCnSigma(pTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(nTrackXi, AliPID::kPion);
            }
            else
            {
                lTPCNSigProton = GetTPCnSigma(nTrackXi, AliPID::kProton);
                lTPCNSigLambdaPion = GetTPCnSigma(pTrackXi, AliPID::kPion);
            }

            lTPCNSigBachelorPion = GetTPCnSigma(bTrackXi, AliPID::kPion);

            lDCALambdaDaughters = TMath::Abs(xiAOD->DcaV0Daughters());
            lDCAXiDaughters = TMath::Abs(xiAOD->DcaXiDaughters());

            lDCALambdaPV =
                TMath::Abs(xiAOD->DcaV0ToPrimVertex());
            lDCAXiPV =
                TMath::Abs(xiAOD->DcaXiToPrimVertex(fPosPV[0], fPosPV[1], fPosPV[2]));

            if (!isXiAnti)
            {
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }
            else
            {
                GetImpactParam(nTrackXi, b, bCov);
                lDCALambdaProtonPV = TMath::Abs(b[0]);
                GetImpactParam(pTrackXi, b, bCov);
                lDCALambdaPionPV = TMath::Abs(b[0]);
            }

            GetImpactParam(bTrackXi, b, bCov);
            lDCABachelorPionPV = TMath::Abs(b[0]);

            lPosXi[0] = xiAOD->DecayVertexXiX();
            lPosXi[1] = xiAOD->DecayVertexXiY();
            lPosXi[2] = xiAOD->DecayVertexXiZ();
            lXiCPA = TMath::Abs(xiAOD->CosPointingAngleXi(fPosPV[0], fPosPV[1], fPosPV[2]));
            if (fLambdaCPAtoXi)
                lLambdaCPA = TMath::Abs(xiAOD->CosPointingAngle(lPosXi));
            else
                lLambdaCPA = TMath::Abs(xiAOD->CosPointingAngle(fPosPV));

            lMomXi[0] = xiAOD->MomXiX();
            lMomXi[1] = xiAOD->MomXiY();
            lMomXi[2] = xiAOD->MomXiZ();
            lXimomsum =
                TMath::Sqrt(lMomXi[0] * lMomXi[0] + lMomXi[1] * lMomXi[1] +
                            lMomXi[2] * lMomXi[2]);
            lXieta =
                0.5 * TMath::Log((lXimomsum + lMomXi[2]) /
                                 (lXimomsum - lMomXi[2] + 1.e-13));

            lPosV0[0] = xiAOD->DecayVertexV0X();
            lPosV0[1] = xiAOD->DecayVertexV0Y();
            lPosV0[2] = xiAOD->DecayVertexV0Z();
            radiusV0 = TMath::Hypot(lPosV0[0], lPosV0[1]);
            radiusXi = TMath::Hypot(lPosXi[0] - fPosPV[0], lPosXi[1] - fPosPV[1]);

            lMassXi = xiAOD->MassXi();
            lMassV0 = -999;
            if (!isXiAnti)
                lMassV0 = xiAOD->MassLambda();
            else
                lMassV0 = xiAOD->MassAntiLambda();

            // Tree
            fTreeXi1530Rec->TPCNsigLambdaProton = lTPCNSigProton;
            fTreeXi1530Rec->TPCNsigLambdaPion = lTPCNSigLambdaPion;
            fTreeXi1530Rec->TPCNsigBachelorPion = lTPCNSigBachelorPion;
            fTreeXi1530Rec->DCALambdaDaughters = lDCALambdaDaughters;
            fTreeXi1530Rec->DCAXiDaughters = lDCAXiDaughters;
            fTreeXi1530Rec->DCALambdaPV = lDCALambdaPV;
            fTreeXi1530Rec->DCALambdaProtonPV = lDCALambdaProtonPV;
            fTreeXi1530Rec->DCALambdaPionPV = lDCALambdaPionPV;
            fTreeXi1530Rec->DCABachelorPionPV = lDCABachelorPionPV;
            fTreeXi1530Rec->DCAXiPV = lDCAXiPV;
            fTreeXi1530Rec->V0CPA = lLambdaCPA;
            fTreeXi1530Rec->XiCPA = lXiCPA;
            fTreeXi1530Rec->XiEta = lXieta;
            fTreeXi1530Rec->LambdaRadius = radiusV0;
            fTreeXi1530Rec->XiRadius = radiusXi;
            fTreeXi1530Rec->V0Mass = lMassV0;
            fTreeXi1530Rec->XiMass = lMassXi;
        }

        for (UInt_t j = 0; j < nTracks; j++)
        {
            track1 = (AliVTrack *)fEvt->GetTrack(fGoodTrackArray[j]);
            if (!track1)
                continue;

            if (track1->GetID() == pID ||
                track1->GetID() == nID ||
                track1->GetID() == bID)
                continue;

            vecPion.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), pionMass);

            vecXi1530 = vecXi + vecPion; // vecXi = cascade, vecPion=pion
                                         // Mass cut
            if ((vecXi1530.M() > fXi1530MassHigh) ||
                (vecXi1530.M() < fXi1530MassLow))
                continue;

            // Y cut
            if ((vecXi1530.Rapidity() > fXi1530YCutHigh) ||
                (vecXi1530.Rapidity() < fXi1530YCutLow))
                continue;

            // DCA and CPA cut for resonance
            if (!fIsAOD)
            { // Only for ESD
                for (int i = 0; i < 21; i++)
                    cv[i] = 0;
                const AliESDVertex *vtxT3D = ((AliESDEvent *)fEvt)->GetPrimaryVertex();
                AliExternalTrackParam pionTrk(*(AliESDtrack *)fEvt->GetTrack(fGoodTrackArray[j]));
                AliExternalTrackParam xiTrk(lPosXi, lMomXi, cv, xiESD->Charge());
                AliExternalTrackParam *pPionTrk = &pionTrk, *pXiTrk = &xiTrk;
                pPionTrk->PropagateToDCA(vtxT3D, fMagField, 250, dztemp, covartemp);
                pXiTrk->PropagateToDCA(vtxT3D, fMagField, 250, dztemp, covartemp);

                AliESDVertex *vertexESD = 0;
                AliAODVertex *vertexAOD = 0;

                AliKFParticle::SetField(fMagField);
                AliKFVertex vertexKF;

                AliESDtrack *esdTrackPion = (AliESDtrack *)pPionTrk;
                AliKFParticle daughterKFPion(*esdTrackPion, 211);
                vertexKF.AddDaughter(daughterKFPion);

                AliESDtrack *esdTrackXi = (AliESDtrack *)pXiTrk;
                AliKFParticle daughterKFXi(*esdTrackXi, 211);
                vertexKF.AddDaughter(daughterKFXi);

                vertexESD = new AliESDVertex(vertexKF.Parameters(),
                                             vertexKF.CovarianceMatrix(),
                                             vertexKF.GetChi2(),
                                             vertexKF.GetNContributors());

                vertexESD->GetXYZ(pos);       // position
                vertexESD->GetCovMatrix(cov); // covariance matrix
                chi2perNDF = vertexESD->GetChi2toNDF();
                dispersion = vertexESD->GetDispersion();
                delete vertexESD;
                vertexESD = NULL;
                vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, 2); // Resonance Vertex

                // Values
                vertexAOD->GetXYZ(lPosRsnVtx); // Secondary vertex
                // Double_t deltaPos[3]{lPosRsnVtx[0] - fPosPV[0], lPosRsnVtx[1] - fPosPV[1], lPosRsnVtx[2] - fPosPV[2]};  // Distance from PV
                dca = pPionTrk->GetDCA(pXiTrk, fMagField, xdummy, ydummy); // DCA of resonance daughters
                radiusRsn = TMath::Hypot(lPosRsnVtx[0] - fPosPV[0], lPosRsnVtx[1] - fPosPV[1]);

                // CPA to PV
                TVector3 mom(vecXi1530.Px(), vecXi1530.Py(), vecXi1530.Pz());
                TVector3 fline(lPosRsnVtx[0] - fPosPV[0],
                               lPosRsnVtx[1] - fPosPV[1],
                               lPosRsnVtx[2] - fPosPV[2]);

                Double_t ptot2 = mom.Mag2() * fline.Mag2();
                if (ptot2 > 0)
                {
                    Double_t cos = mom.Dot(fline) / TMath::Sqrt(ptot2);
                    if (cos > 1.0)
                        cos = 1.0;
                    if (cos < -1.0)
                        cos = -1.0;
                    fTreeXi1530Rec->Xi1530CPA = cos;
                }
                else
                    fTreeXi1530Rec->Xi1530CPA = 0.0;

                fTreeXi1530Rec->Xi1530Radius = radiusRsn;
                fTreeXi1530Rec->Xi1530DCA = dca;
            }
            else
            {
                fTreeXi1530Rec->Xi1530Radius = -999;
                fTreeXi1530Rec->Xi1530DCA = -999;
                fTreeXi1530Rec->Xi1530CPA = -999;
            }

            GetImpactParam(track1, b, bCov);

            fTreeXi1530Rec->TPCNsigXi1530Pion = GetTPCnSigma(track1, AliPID::kPion);
            fTreeXi1530Rec->Xi1530PionEta = track1->Eta();
            fTreeXi1530Rec->Xi1530PionPVz = TMath::Abs(b[1]);
            fTreeXi1530Rec->Xi1530PionXYVertexSigma = TMath::Abs(b[0]) / (0.0026 + 0.0050 / track1->Pt());

            fTreeXi1530Rec->Xi1530Pt = vecXi1530.Pt();
            fTreeXi1530Rec->Xi1530Mass = vecXi1530.M();
            fTreeXi1530Rec->Centrality = fCent;
            fTreeXi1530Rec->zVertex = fPosPV[2];

            if (fIsMC)
            {
                if (IsTrueXi1530(fGoodXiArray[i], fGoodTrackArray[j]))
                    fTreeXi1530Gen.MCflag = 1; // MCflag Xi1530
                else
                    fTreeXi1530Gen.MCflag = 2; // MCflag -> not true
                fTreeXi1530Gen.isReconstructed = true;
            }
            // Charge
            auto chargeFlag = 0;
            if (isXiAnti)
                chargeFlag++;
            if (track1->Charge() < 0)
                chargeFlag = chargeFlag + 2;
            fTreeXi1530Rec->Antiflag = chargeFlag; // Antiflag -> Charge flag
                                                   // 0: Xi- + pi+
                                                   // 1: Xi+ + pi+ // LS
                                                   // 2: Xi- + pi- // LS
                                                   // 3: Xi+ + pi-
            fTree->Fill();
        } // pion loop
    }
}

void AliAnalysisTaskXi1530PbPb::FillMCinput(AliMCEvent *fMCEvent, int Fillbin)
{
    int sign = kAllType;
    auto binAnti = kNormal;
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

            if (TMath::Abs(mcInputTrack->GetPdgCode()) != kXiStarCode)
                continue;
            if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
                continue;
            if (fFillQAPlot)
                fHistos->FillTH1("MCQA/hMCGeneratedXi1530Pt", mcInputTrack->Pt());

            // Y cut
            if ((mcInputTrack->Y() > fXi1530YCutHigh) ||
                (mcInputTrack->Y() < fXi1530YCutLow))
                continue;
            if (mcInputTrack->GetPdgCode() > 0)
                binAnti = kNormal;
            else
                binAnti = kAnti;

            sign = kXiStar_GEN + Fillbin;
            FillTHnSparse("Xi1530_mc",
                          {(double)binAnti,
                           (double)sign,
                           (double)fCent,
                           mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
            if (fFillTree)
            {
                fTreeXi1530Gen.isReconstructed = false;
                fTreeXi1530Gen.MCflag = 1;
                fTreeXi1530Gen.ptMC = mcInputTrack->Pt();
                fTreeXi1530Gen.Antiflag = (mcInputTrack->GetPdgCode() > 0) ? 0 : 3;
                fTree->Fill();
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

            if (TMath::Abs(mcInputTrack->GetPdgCode()) != kXiStarCode)
                continue;
            if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
                continue;
            if (fFillQAPlot)
                fHistos->FillTH1("MCQA/hMCGeneratedXi1530Pt", mcInputTrack->Pt());

            // Y cut
            if ((mcInputTrack->Y() > fXi1530YCutHigh) ||
                (mcInputTrack->Y() < fXi1530YCutLow))
                continue;

            if (mcInputTrack->GetPdgCode() > 0)
                binAnti = kNormal;
            else
                binAnti = kAnti;

            sign = kXiStar_GEN + Fillbin;
            FillTHnSparse("Xi1530_mc",
                          {(double)binAnti,
                           (double)sign,
                           (double)fCent,
                           mcInputTrack->Pt(),
                           mcInputTrack->GetCalcMass()});
            if (fFillTree)
            {
                fTreeXi1530Gen.isReconstructed = false;
                fTreeXi1530Gen.MCflag = 1;
                fTreeXi1530Gen.ptMC = mcInputTrack->Pt();
                fTreeXi1530Gen.Antiflag = (mcInputTrack->GetPdgCode() > 0) ? 0 : 3;
                fTree->Fill();
            }
        }
    }
}

Bool_t AliAnalysisTaskXi1530PbPb::IsTrueXi1530(UInt_t xiIndex, UInt_t pionIndex)
{
    Bool_t TrueXi1530 = kFALSE;
    AliVTrack *track1;
    AliESDcascade *xiESD;
    AliAODcascade *xiAOD;

    track1 = (AliVTrack *)fEvt->GetTrack(pionIndex);

    if (!fIsAOD)
    {
        xiESD = ((AliESDEvent *)fEvt)->GetCascade(xiIndex);
        if (!xiESD)
            return kFALSE;
        AliESDtrack *pTrackXi =
            ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetPindex()));
        AliESDtrack *nTrackXi =
            ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetNindex()));
        AliESDtrack *bTrackXi =
            ((AliESDEvent *)fEvt)->GetTrack(TMath::Abs(xiESD->GetBindex()));

        TParticle *MCXiD2esd =
            (TParticle *)fMCEvent->GetTrack(TMath::Abs(bTrackXi->GetLabel()))->Particle();
        TParticle *MCLamD1esd;
        TParticle *MCLamD2esd;
        TParticle *MCLamesd;
        TParticle *MCXiesd;
        TParticle *MCXiStaresd;
        TParticle *MCXiStarD2esd;

        if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode)
        {
            // D2esd->pion
            MCLamD1esd = (TParticle *)fMCEvent->GetTrack(TMath::Abs(pTrackXi->GetLabel()))->Particle();
            MCLamD2esd = (TParticle *)fMCEvent->GetTrack(TMath::Abs(nTrackXi->GetLabel()))->Particle();
            if (MCLamD1esd->GetMother(0) ==
                MCLamD2esd->GetMother(0))
            {
                // Same mother(lambda)
                if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                    (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                         kProtonCode))
                {
                    // Lamda daugthers check #1
                    MCLamesd = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCLamD1esd->GetMother(0)))->Particle();
                    if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                        kLambdaCode)
                    {
                        // Lambda check
                        if (MCLamesd->GetMother(0) ==
                            MCXiD2esd->GetMother(0))
                        {
                            // Lambda+pion(D2esd) mother check
                            MCXiesd = (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCLamesd->GetMother(0)))->Particle();
                            if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                                kXiCode)
                            {
                                // Xi Check
                                MCXiStarD2esd =
                                    (TParticle *)fMCEvent->GetTrack(TMath::Abs(track1->GetLabel()))->Particle();
                                if (MCXiesd->GetMother(0) ==
                                    MCXiStarD2esd->GetMother(0))
                                {
                                    // Xi+pion mother check
                                    MCXiStaresd =
                                        (TParticle *)fMCEvent->GetTrack(TMath::Abs(MCXiesd->GetMother(0)))->Particle();
                                    if (TMath::Abs(MCXiStaresd->GetPdgCode()) ==
                                        kXiStarCode)
                                    {
                                        // Xi1530 check
                                        if (fIsPrimaryMC)
                                        {
                                            if (MCXiStaresd->IsPrimary())
                                            {
                                                TrueXi1530 = kTRUE;
                                            } // Primary(input) Xi1530 check
                                        }
                                        else
                                        {
                                            TrueXi1530 = kTRUE;
                                        }
                                    } // Xi1530 check
                                }     // Xi+pion mother check
                            }         // Xi Check
                        }             // Lambda+pion(D2esd) mother check
                    }                 // Lambda check
                }                     // Lamda daugthers check
            }                         // Same mother(lambda)
        }                             // D2esd->pion
        return TrueXi1530;
    }
    else
    {
        xiAOD = ((AliAODEvent *)fEvt)->GetCascade(xiIndex);
        if (!xiAOD)
            return kFALSE;
        AliAODTrack *pTrackXi =
            (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(0));
        AliAODTrack *nTrackXi =
            (AliAODTrack *)(xiAOD->GetSecondaryVtx()->GetDaughter(1));
        AliAODTrack *bTrackXi =
            (AliAODTrack *)(xiAOD->GetDecayVertexXi()->GetDaughter(0));

        AliAODMCParticle *MCXiD2esd =
            (AliAODMCParticle *)fMCArray->At(TMath::Abs(bTrackXi->GetLabel()));
        AliAODMCParticle *MCLamD1esd;
        AliAODMCParticle *MCLamD2esd;
        AliAODMCParticle *MCLamesd;
        AliAODMCParticle *MCXiesd;
        AliAODMCParticle *MCXiStaresd;
        AliAODMCParticle *MCXiStarD2esd;

        if (TMath::Abs(MCXiD2esd->GetPdgCode()) == kPionCode)
        {
            // D2esd->pion
            MCLamD1esd = (AliAODMCParticle *)fMCArray->At(TMath::Abs(pTrackXi->GetLabel()));
            MCLamD2esd = (AliAODMCParticle *)fMCArray->At(TMath::Abs(nTrackXi->GetLabel()));
            if (MCLamD1esd->GetMother() ==
                MCLamD2esd->GetMother())
            {
                // Same mother(lambda)
                if ((TMath::Abs(MCLamD1esd->GetPdgCode()) == kProtonCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) == kPionCode) ||
                    (TMath::Abs(MCLamD1esd->GetPdgCode()) == kPionCode &&
                     TMath::Abs(MCLamD2esd->GetPdgCode()) ==
                         kProtonCode))
                {
                    // Lamda daugthers check #1
                    MCLamesd = (AliAODMCParticle *)fMCArray->At(
                        TMath::Abs(MCLamD1esd->GetMother()));
                    if (TMath::Abs(MCLamesd->GetPdgCode()) ==
                        kLambdaCode)
                    {
                        // Lambda check
                        if (MCLamesd->GetMother() ==
                            MCXiD2esd->GetMother())
                        {
                            // Lambda+pion(D2esd) mother check
                            MCXiesd = (AliAODMCParticle *)fMCArray->At(
                                TMath::Abs(MCLamesd->GetMother()));
                            if (TMath::Abs(MCXiesd->GetPdgCode()) ==
                                kXiCode)
                            {
                                // Xi Check
                                MCXiStarD2esd = (AliAODMCParticle *)fMCArray->At(
                                    TMath::Abs(track1->GetLabel()));
                                if (MCXiesd->GetMother() ==
                                    MCXiStarD2esd->GetMother())
                                {
                                    // Xi+pion mother check
                                    MCXiStaresd = (AliAODMCParticle *)fMCArray->At(
                                        TMath::Abs(MCXiesd->GetMother()));
                                    if (TMath::Abs(MCXiStaresd->GetPdgCode()) ==
                                        kXiStarCode)
                                    {
                                        // Xi1530 check
                                        if (fIsPrimaryMC)
                                        {
                                            if (MCXiStaresd->IsPrimary())
                                            {
                                                TrueXi1530 = kTRUE;
                                            } // Primary(input) Xi1530 check
                                        }
                                        else
                                        {
                                            TrueXi1530 = kTRUE;
                                        }
                                    } // Xi1530 check
                                }     // Xi+pion mother check
                            }         // Xi Check
                        }             // Lambda+pion(D2esd) mother check
                    }                 // Lambda check
                }                     // Lamda daugthers check
            }                         // Same mother(lambda)
        }                             // D2esd->pion
        return TrueXi1530;
    }
}

THnSparse *AliAnalysisTaskXi1530PbPb::CreateTHnSparse(
    TString name,
    TString title,
    Int_t ndim,
    std::vector<TAxis> bins,
    Option_t *opt)
{
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    const TAxis *axises[bins.size()];
    for (UInt_t i = 0; i < bins.size(); i++)
        axises[i] = &bins[i];
    THnSparse *h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
    return h;
}

Long64_t AliAnalysisTaskXi1530PbPb::FillTHnSparse(TString name,
                                                  std::vector<Double_t> x,
                                                  Double_t w)
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

Long64_t AliAnalysisTaskXi1530PbPb::FillTHnSparse(THnSparse *h,
                                                  std::vector<Double_t> x,
                                                  Double_t w)
{
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    if (int(x.size()) != h->GetNdimensions())
    {
        std::cout << "ERROR : wrong sized of array while Fill " << h->GetName() << std::endl;
        exit(1);
    }

    return h->Fill(&x.front(), w);
}

TAxis AliAnalysisTaskXi1530PbPb::AxisFix(TString name,
                                         int nbin,
                                         Double_t xmin,
                                         Double_t xmax)
{
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(nbin, xmin, xmax);
    axis.SetName(name);
    return axis;
}

TAxis AliAnalysisTaskXi1530PbPb::AxisStr(TString name,
                                         std::vector<TString> bin)
{
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
    UInt_t i = 1;
    for (auto blabel : bin)
        ax.SetBinLabel(i++, blabel);
    return ax;
}

TAxis AliAnalysisTaskXi1530PbPb::AxisVar(TString name,
                                         std::vector<Double_t> bin)
{
    // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
    // Original author: Beomkyu Kim
    TAxis axis(bin.size() - 1, &bin.front());
    axis.SetName(name);
    return axis;
}

double AliAnalysisTaskXi1530PbPb::GetTPCnSigma(AliVTrack *track,
                                               AliPID::EParticleType type)
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

void AliAnalysisTaskXi1530PbPb::FillTrackToEventPool()
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
    for (UInt_t i = 0; i < fGoodTrackArray.size(); i++)
    {
        goodtrack = (AliVTrack *)fEvt->GetTrack(fGoodTrackArray[i]);
        if (!goodtrack)
            continue;
        etl->push_back((AliVTrack *)goodtrack->Clone());
    }

    if (!fGoodTrackArray.size())
        ep->pop_back();
    if ((int)ep->size() > (int)fnMix)
    {
        for (auto it : ep->front())
            delete it;
        ep->pop_front();
    }
}

void AliAnalysisTaskXi1530PbPb::GetImpactParam(AliVTrack *track,
                                               Float_t p[2],
                                               Float_t cov[3])
{
    AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
    if (nanoT)
        nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
    else
        track->GetImpactParameters(p, cov);
}

Bool_t AliAnalysisTaskXi1530PbPb::IsSelectedTPCGeoCut(AliAODTrack *track)
{
    Bool_t checkResult = kTRUE;

    AliESDtrack esdTrack(track);
    esdTrack.SetTPCClusterMap(track->GetTPCClusterMap());
    esdTrack.SetTPCSharedMap(track->GetTPCSharedMap());
    esdTrack.SetTPCPointsF(track->GetTPCNclsF());

    auto nCrossedRowsTPC = esdTrack.GetTPCCrossedRows();
    auto lengthInActiveZoneTPC = esdTrack.GetLengthInActiveZone(
        0, fTPCActiveLengthCutDeltaY, fTPCActiveLengthCutDeltaZ, fMagField);
    auto cutGeoNcrNclLength = fRequireCutGeoNcrNclLength -
                              TMath::Power(TMath::Abs(esdTrack.GetSigned1Pt()),
                                           fRequireCutGeoNcrNclGeom1Pt);

    if (lengthInActiveZoneTPC < cutGeoNcrNclLength)
        checkResult = false;
    if (nCrossedRowsTPC < fCutGeoNcrNclFractionNcr * cutGeoNcrNclLength)
        checkResult = false;
    if (esdTrack.GetTPCncls() < fCutGeoNcrNclFractionNcl * cutGeoNcrNclLength)
        checkResult = false;

    return checkResult;
}
Bool_t AliAnalysisTaskXi1530PbPb::IsSelectedTPCGeoCut(AliESDtrack *track)
{
    Bool_t checkResult = kTRUE;

    auto nCrossedRowsTPC = track->GetTPCCrossedRows();
    auto lengthInActiveZoneTPC = track->GetLengthInActiveZone(
        0, fTPCActiveLengthCutDeltaY, fTPCActiveLengthCutDeltaZ, fMagField);
    auto cutGeoNcrNclLength = fRequireCutGeoNcrNclLength -
                              TMath::Power(TMath::Abs(track->GetSigned1Pt()),
                                           fRequireCutGeoNcrNclGeom1Pt);

    if (lengthInActiveZoneTPC < cutGeoNcrNclLength)
        checkResult = false;
    if (nCrossedRowsTPC < fCutGeoNcrNclFractionNcr * cutGeoNcrNclLength)
        checkResult = false;
    if (track->GetTPCncls() < fCutGeoNcrNclFractionNcl * cutGeoNcrNclLength)
        checkResult = false;

    return checkResult;
}
