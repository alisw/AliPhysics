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

/* AliAnalysisTaskXi1820BH
 *
 *  Test code for the reconstructing Xi(1820)^{0, +-}
 *
 *  Author: Bong-Hwi Lim
 *
 */

#include <TPDGCode.h>
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
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliEventCuts.h"
#include "AliMultSelectionTask.h"
#include "AliPIDResponse.h"
#include "TChain.h"

// for NanoAOD
#include <AliNanoAODHeader.h>
#include <AliNanoAODTrack.h>

#include "AliAODv0.h"
#include "AliAnalysisTaskXi1820BH.h"
#include "AliAnalysisTaskTrackMixer.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "THistManager.h"

const Double_t kaonMass = AliPID::ParticleMass(AliPID::kKaon);
const Double_t k0sMass = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
const Double_t v0Mass = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();

enum
{
  kNormal = 1,
  kAnti = 2
};
enum
{
  kXi1820ZeroCode = 123324, // Xi(1820)0
  kXi1820MCode = 123314,    // Xi(1820)-
  kLambdaCode = 3122,       // Lambda
  kProtonCode = 2212,       // Proton+
  kPionCode = 211,          // Pion+
  kKaonCode = 321,          // Kaon-
  kK0sCode = 310,           // K0s
};
enum // Need to be modified
{
  kXi1820Zero = 1,
  kXi1820M,
  kXi1820Zero_MIX,
  kXi1820M_MIX,
  kXi1820Zero_NOT,
  kXi1820M_NOT,
  kAllType
};
enum // Need to be modified
{
  kXi1820Zero_GEN = 1, // 1
  kXi1820M_GEN,
  kXi1820Zero_GEN_INEL10,
  kXi1820M_GEN_INEL10,
  kXi1820Zero_GEN_INEL10_IGZ,
  kXi1820M_GEN_INEL10_IGZ,
  kXi1820Zero_GEN_TRIG,
  kXi1820M_GEN_TRIG,
  kXi1820Zero_REC,
  kXi1820M_REC
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

class AliAnalysisTaskXi1820BH;

ClassImp(AliAnalysisTaskXi1820BH)
    AliAnalysisTaskXi1820BH::AliAnalysisTaskXi1820BH()
    : AliAnalysisTaskSE(),
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
      fEvt(nullptr),
      fMCEvent(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fMCArray(nullptr),
      fIsNano(kFALSE),
      fSetMixing(kFALSE),
      fSetZero(kFALSE),
      fFillQAPlot(kTRUE),
      fIsMC(kFALSE),
      fIsPrimaryMC(kFALSE),
      fIsINEL(kFALSE),
      fIsHM(kFALSE),
      fUseAsymmCut(kFALSE),
      fOnlyUseOnTheFlyV0(kFALSE),
      fEMpoolKaon(0),
      fEMpoolK0s(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fMagField(0),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fFilterBit(32.0),
      fTPCNsigXi1820ZeroKaonCut(2.0),
      fTOFNsigXi1820ZeroKaonCut(3.0),
      fXi1820ZeroKaonEtaCut(0.8),
      fXi1820ZeroKaonZVertexCut(2.0),
      fXi1820KaonZeroXYVertexSigmaCut(7.0),
      fTPCNsigK0sPionCut(3.0),
      fTPCNsigLambdaProtonCut(3.0),
      fTPCNsigLambdaPionCut(3.0),
      fDCADistK0sDaughtersCut(1.0),
      fDCADistLambdaDaughtersCut(1.0),
      fDCArDistK0sPVCut(0.3),
      fDCArDistLambdaPVCut(0.4),
      fDCAPositiveTrackK0s(0.06),
      fDCANegativeTrackK0s(0.06),
      fDCAPositiveTrackLambda(0.05),
      fDCANegativeTrackLambda(0.05),
      fK0sCosineOfPointingAngleCut(0.97),
      fLambdaCosineOfPointingAngleCut(0.99),
      fMaxK0sRapidity(0.8),
      fMaxLambdaRapidity(0.8),
      fK0sLowRadius(0.5),
      fK0sHighRadius(200.0),
      fLambdaLowRadius(0.5),
      fLambdaHighRadius(200.0),
      fK0sLifetime(30.0),
      fLambdaLifetime(30.0),
      fK0sMassWindowCut(0.01),
      fLambdaMassWindowCut(0.01),
      fXi1820YCutHigh(0.5),
      fXi1820YCutLow(-0.5),
      fXi1820AsymmCutHigh(0.0),
      fXi1820AsymmCutLow(1.0),
      fGoodKaonArray(),
      fGoodK0sArray(),
      fGoodLambdaArray()
{
  /// Default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskXi1820BH::AliAnalysisTaskXi1820BH(const char *name, Bool_t MCcase)
    : AliAnalysisTaskSE(name),
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
      fEvt(nullptr),
      fMCEvent(nullptr),
      fHistos(nullptr),
      fVertex(nullptr),
      fMCArray(nullptr),
      fIsNano(kFALSE),
      fSetMixing(kFALSE),
      fSetZero(kFALSE),
      fFillQAPlot(kTRUE),
      fIsMC(MCcase),
      fIsPrimaryMC(kFALSE),
      fIsINEL(kFALSE),
      fIsHM(kFALSE),
      fUseAsymmCut(kFALSE),
      fOnlyUseOnTheFlyV0(kFALSE),
      fEMpoolKaon(0),
      fEMpoolK0s(0),
      fBinCent(),
      fBinZ(),
      fPosPV(),
      fMagField(0),
      fCent(-1),
      fnMix(10),
      fCentBin(-1),
      fZbin(-1),
      fFilterBit(32),
      fTPCNsigXi1820ZeroKaonCut(2.0),
      fTOFNsigXi1820ZeroKaonCut(3.0),
      fXi1820ZeroKaonEtaCut(0.8),
      fXi1820ZeroKaonZVertexCut(2.0),
      fXi1820KaonZeroXYVertexSigmaCut(7.0),
      fTPCNsigK0sPionCut(3.0),
      fDCADistK0sDaughtersCut(1.0),
      fDCArDistK0sPVCut(0.3),
      fDCAPositiveTrackK0s(0.06),
      fDCANegativeTrackK0s(0.06),
      fK0sCosineOfPointingAngleCut(0.97),
      fMaxK0sRapidity(0.8),
      fK0sLowRadius(0.5),
      fK0sHighRadius(200.0),
      fK0sLifetime(30.0),
      fK0sMassWindowCut(0.01),
      fTPCNsigLambdaProtonCut(3.0),
      fTPCNsigLambdaPionCut(3.0),
      fDCADistLambdaDaughtersCut(1.0),
      fDCArDistLambdaPVCut(0.4),
      fDCAPositiveTrackLambda(0.05),
      fDCANegativeTrackLambda(0.05),
      fLambdaCosineOfPointingAngleCut(0.99),
      fMaxLambdaRapidity(0.8),
      fLambdaLowRadius(0.5),
      fLambdaHighRadius(200.0),
      fLambdaLifetime(30.0),
      fLambdaMassWindowCut(0.01),
      fXi1820YCutHigh(0.5),
      fXi1820YCutLow(-0.5),
      fXi1820AsymmCutHigh(0.0),
      fXi1820AsymmCutLow(1.0),
      fGoodKaonArray(),
      fGoodK0sArray(),
      fGoodLambdaArray()
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskXi1820BH::~AliAnalysisTaskXi1820BH() {}
//_____________________________________________________________________________
void AliAnalysisTaskXi1820BH::UserCreateOutputObjects()
{
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fHistos = new THistManager("Xi1820hists");
  auto binAnti = AxisStr("AType", {"Normal", "Anti"});
  auto binType = (!fIsMC)
                     ? AxisStr("Type", {"Xi1820Zero", "Xi1820M",
                                        "Xi1820Zero_mix", "Xi1820M_mix"})
                     : AxisStr("Type", {"Xi1820Zero", "Xi1820M",
                                        "Xi1820Zero_mix", "Xi1820M_mix",
                                        "Xi1820Zero_no", "Xi1820M_no"});
  auto binTypeMC = AxisStr("Type", {"Xi1820Zero_gen",
                                    "Xi1820M_gen",
                                    "Xi1820Zero_gen_inel10",
                                    "Xi1820M_gen_inel10",
                                    "Xi1820Zero_gen_inel10_igz",
                                    "Xi1820M_gen_inel10_igz",
                                    "Xi1820Zero_gen_trig",
                                    "Xi1820M_gen_trig",
                                    "Xi1820Zero_rec",
                                    "Xi1820M_rec"});

  std::vector<double> centaxisbin;
  (fIsHM)
      ? centaxisbin = {0, 0.001, 0.01, 0.05, 0.1}
      : centaxisbin = {-1, 0, 1, 5, 10, 15, 20, 30, 40,
                       50, 60, 70, 80, 90, 100}; // can be use from pp to PbPb
  if (fIsMC)
    centaxisbin = {-1, 0, 0.001, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20,
                   30, 40, 50, 60, 70, 80, 90, 100}; // for general MC
  fBinCent = AxisVar("Cent", centaxisbin);
  auto binPt = AxisFix("Pt", 200, 0, 20);
  auto binMass = AxisFix("Mass", 1400, 1.6, 3.0);
  auto binMassMC = AxisFix("Mass", 800, 1.6, 2.4);
  fBinZ = AxisVar("Z", {-10, -5, -3, -1, 1, 3, 5, 10}); // moderate diff

  CreateTHnSparse("Xi1820_data", "Xi1820_data", 5,
                  {binAnti, binType, fBinCent, binPt, binMass}, "s");
  if (fIsMC)
  {
    auto binTypeMCNorm = AxisStr(
        "Type", {"kAll", "kINEL10", "kINEL_trig", "kINEL_trig_vtx",
                 "kINEL_trig_vtx10", "kINELg0", "kINELg010", "kINELg0_trig",
                 "kINELg0_trig_vtx", "kINELg0_trig_vtx10", "kSelected"});
    CreateTHnSparse("Xi1820_mc", "Xi1820_mc", 5,
                    {binAnti, binTypeMC, fBinCent, binPt, binMassMC}, "s");
    CreateTHnSparse("Normalisation", "", 2, {binTypeMCNorm, fBinCent}, "s");
  }
  fEventCuts.AddQAplotsToList(fHistos->GetListOfHistograms());
  if (fIsHM)
    fHistos->CreateTH1("hMultiplicity", "", 100, 0, 0.1, "s");
  else
    fHistos->CreateTH1("hMultiplicity", "", 101, -1, 100, "s");
  if (fFillQAPlot)
  {
    fHistos->CreateTH2("QA/hTPCPIDKaon", "", 200, 0, 20, 200, -10, 10);
    fHistos->CreateTH2("QA/hTOFPIDKaon", "", 200, 0, 20, 200, -10, 10);
    fHistos->CreateTH2("QA/hTPCTOFPIDKaon", "", 200, -10, 10, 200, -10, 10);
    fHistos->CreateTH1("QA/hEtaKaon", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QA/hDCAPVKaon", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QA/hDCArPVKaon", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QA/hPtKaon", "", 200, 0, 20);
    if (fCheckTPCGeo)
      fHistos->CreateTH1("QA/hTPCGeoCheckKaon", "", 2, -0.5, 1.5, "s");
    fHistos->CreateTH2("QA/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QA/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 200, 0, 200);
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
    if (fSetZero)
    {
      fHistos->CreateTH2("QA/hTPCPIDK0s", "", 200, 0, 20, 200, 0, 200);
      fHistos->CreateTH1("QA/hEtaK0s", "", 20, -1.0, 1.0);
      fHistos->CreateTH1("QA/hDCAPVK0s", "", 30, 0, 3, "s");
      fHistos->CreateTH1("QA/hDCArPVK0s", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QA/hPtK0s", "", 200, 0, 20);
      fHistos->CreateTH2("QA/hTPCPIDK0sPionP", "", 200, 0, 20, 200, -10, 10);
      fHistos->CreateTH2("QA/hTPCPIDK0sPionN", "", 200, 0, 20, 200, -10, 10);
      fHistos->CreateTH1("QA/hDCA_k0sDaughters", "", 300, 0, 3, "s");
      fHistos->CreateTH1("QA/hDCAk0sPV", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QA/hDCAK0sPVPionN", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QA/hDCAK0sPVPionP", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QA/hCosPAk0s", "", 50, 0.95, 1.0, "s");
      fHistos->CreateTH1("QA/hLifetimek0s", "", 50, 0, 50, "s");
      fHistos->CreateTH1("QA/hYk0s", "", 140, -0.7, 0.7, "s");
      fHistos->CreateTH1("QA/hMassK0s", "", 100, 0.4, 0.6, "s");
      fHistos->CreateTH1("QA/hK0sRxy", "", 200, 0, 200);
    }
    fHistos->CreateTH2("QA/hXi1820Asymm", "", 100, 0, 1, 100, 0, 100);
    if (fIsMC)
    {
      fHistos->CreateTH2("QA/hXi1820Asymm_true", "", 100, 0, 1, 100, 0, 100);
      fHistos->CreateTH2("QA/hXi1820Asymm_true_selected", "", 100, 0, 1, 100,
                         0, 100);
    }
    fHistos->CreateTH1("QA/hLambdaAntiCheck", "", 4, -0.5, 3.5);

    fHistos->CreateTH2("QAcut/hTPCPIDKaon", "", 200, 0, 20, 200, -10, 10);
    fHistos->CreateTH2("QAcut/hTOFPIDKaon", "", 200, 0, 20, 200, -10, 10);
    fHistos->CreateTH2("QAcut/hTPCTOFPIDKaon", "", 200, -10, 10, 200, -10, 10);
    fHistos->CreateTH1("QAcut/hEtaKaon", "", 20, -1.0, 1.0);
    fHistos->CreateTH1("QAcut/hDCAPVKaon", "", 30, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCArPVKaon", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hPtKaon", "", 200, 0, 20);
    if (fCheckTPCGeo)
      fHistos->CreateTH1("QAcut/hTPCGeoCheckKaon", "", 2, -0.5, 1.5, "s");
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaProton", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaProton", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH2("QAcut/hTPCPIDAntiLambdaPion", "", 200, 0, 20, 200, 0, 200);
    fHistos->CreateTH1("QAcut/hDCA_lambdaDaughters", "", 300, 0, 3, "s");
    fHistos->CreateTH1("QAcut/hDCAlambdaPV", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hDCALambdaPVProton", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hDCALambdaPVPion", "", 50, 0, 0.5, "s");
    fHistos->CreateTH1("QAcut/hCosPAlambda", "", 50, 0.95, 1.0, "s");
    fHistos->CreateTH1("QAcut/hLifetimelambda", "", 50, 0, 50, "s");
    fHistos->CreateTH1("QAcut/hYlambda", "", 140, -0.7, 0.7, "s");
    fHistos->CreateTH1("QAcut/hMassLambda", "", 80, 1.08, 1.16, "s");
    fHistos->CreateTH1("QAcut/hLambdaRxy", "", 200, 0, 200);
    if (fSetZero)
    {
      fHistos->CreateTH2("QAcut/hTPCPIDK0s", "", 200, 0, 20, 200, 0, 200);
      fHistos->CreateTH1("QAcut/hEtaK0s", "", 20, -1.0, 1.0);
      fHistos->CreateTH1("QAcut/hDCAPVK0s", "", 30, 0, 3, "s");
      fHistos->CreateTH1("QAcut/hDCArPVK0s", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QAcut/hPtK0s", "", 200, 0, 20);
      fHistos->CreateTH2("QAcut/hTPCPIDK0sPionP", "", 200, 0, 20, 200, -10, 10);
      fHistos->CreateTH2("QAcut/hTPCPIDK0sPionN", "", 200, 0, 20, 200, -10, 10);
      fHistos->CreateTH1("QAcut/hDCA_k0sDaughters", "", 300, 0, 3, "s");
      fHistos->CreateTH1("QAcut/hDCAk0sPV", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QAcut/hDCAK0sPVPionN", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QAcut/hDCAK0sPVPionP", "", 50, 0, 0.5, "s");
      fHistos->CreateTH1("QAcut/hCosPAk0s", "", 50, 0.95, 1.0, "s");
      fHistos->CreateTH1("QAcut/hLifetimek0s", "", 50, 0, 50, "s");
      fHistos->CreateTH1("QAcut/hYk0s", "", 140, -0.7, 0.7, "s");
      fHistos->CreateTH1("QAcut/hMassK0s", "", 100, 0.4, 0.6, "s");
      fHistos->CreateTH1("QAcut/hK0sRxy", "", 200, 0, 200);
    }
    fHistos->CreateTH2("QAcut/hXi1820Asymm", "", 100, 0, 1, 100, 0, 100);
    if (fIsMC)
    {
      fHistos->CreateTH2("QAcut/hXi1820Asymm_true", "", 100, 0, 1, 100, 0,
                         100);
      fHistos->CreateTH2("QAcut/hXi1820Asymm_true_selected", "", 100, 0, 1,
                         100, 0, 100);
    }
    fHistos->CreateTH1("QAcut/hLambdaAntiCheck", "", 4, -0.5, 3.5);
  }

  fEMpoolKaon.resize(fBinCent.GetNbins() + 1, std::vector<eventpool>(fBinZ.GetNbins() + 1));
  fEMpoolK0s.resize(fBinCent.GetNbins() + 1, std::vector<v0pool>(fBinZ.GetNbins() + 1));

  PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskXi1820BH::UserExec(Option_t *)
{
  AliVEvent *event = InputEvent();
  if (!event)
  {
    PostData(1, fHistos->GetListOfHistograms());

    AliInfo("Could not retrieve event");
    return;
  }
  AliNanoAODHeader *nanoHeader = dynamic_cast<AliNanoAODHeader *>(fInputEvent->GetHeader()); // NanoAOD Header
  fEvt = dynamic_cast<AliAODEvent *>(event);
  if (!fEvt)
  {
    PostData(1, fHistos->GetListOfHistograms());

    return;
  }

  AliInputEventHandler *inputHandler =
      (AliInputEventHandler *)AliAnalysisManager::GetAnalysisManager()
          ->GetInputEventHandler();

  if (nanoHeader && !fIsNano)
    fIsNano = kTRUE; // set flag for NanoAOD
  bool IsEvtSelected = (!nanoHeader) ? fEventCuts.AcceptEvent(event) : true;
  fCent = (!nanoHeader) ? AliMultSelectionTask::IsINELgtZERO(event) ? fEventCuts.GetCentrality() : -0.5 : nanoHeader->GetCentr("V0M");
  if (fIsNano) // NanoAOD INEL centrality
  {
    static int inel_index = -1;
    if (inel_index < 0)
      inel_index = nanoHeader->GetVarIndex("cstINELgt0"); // NanoAOD INEL>0
    if ((inel_index > 0) && (nanoHeader->GetVar(inel_index) < 0.5))
      fCent = -0.5; // INEL>0 centrality is saved in the -0.5 bin
  }
  fPIDResponse = (AliPIDResponse *)inputHandler->GetPIDResponse();
  if (!fIsNano && !fPIDResponse)
    AliInfo("No PIDd");
  if (fIsMC)
  {
    fMCArray = (TClonesArray *)fEvt->FindListObject("mcparticles"); // AOD Case
    fMCEvent = MCEvent();
    FillMCEventProperties(fMCEvent);
  }

  if (!IsEvtSelected)
  {
    PostData(1, fHistos->GetListOfHistograms());

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

  bool checkKaon = GoodKaonSelection();
  bool checkLambda = GoodLambdaSelection();

  if (checkKaon && checkLambda)
  {
    FillKaonV0(); // Xi1820^{+-}
  }

  if (fSetZero)
  {
    bool checkK0s = GoodK0sSelection();
    if (checkK0s && checkLambda)
      FillK0sV0(); // Xi1820^{0}?
    if (fSetMixing && fGoodK0sArray.size())
      FillK0sToEventPool();
  }
  if (fSetMixing)
  {
    if (fGoodKaonArray.size())
      FillKaonToEventPool(); // use only pion track pool.
  }
  PostData(1, fHistos->GetListOfHistograms());
}
//_____________________________________________________________________________
void AliAnalysisTaskXi1820BH::Terminate(Option_t *) {}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskXi1820BH::GoodKaonSelection()
{
  // Selection of good kaon tracks
  // pT > 0.15 GeV/c
  // |Eta| < 0.8
  // at least 1 of ITS hit
  // Number of crossed TPC rows > 70
  // Ratio of findable clusters > 0.8
  // DCA_xy to PV < 7 sigma
  // DCA_z to PV < 2 cm
  // TPC PID nSigma < 2
  // TOF PID nSigma < 3

  const UInt_t nTracks = fEvt->GetNumberOfTracks();
  fGoodKaonArray.clear();
  AliVTrack *track;
  Float_t b[2];
  Float_t bCov[3];
  Double_t nTPCNSigKaon, nTOFNSigKaon, kaonZ, kaonpT, kaonSigmaDCA_r, kaonDCA_r, lEta;
  Int_t isTPCGeo;

  for (UInt_t it = 0; it < nTracks; it++)
  {
    track = (AliVTrack *)fEvt->GetTrack(it);
    if (!track)
      continue;

    // ---------- Track selection begin ----------
    if (!fIsNano) // AOD case
    {
      if (!((AliAODTrack *)track)->TestFilterBit(fFilterBit))
        continue;
      if (fCheckTPCGeo)
        isTPCGeo = IsSelectedTPCGeoCut(((AliAODTrack *)track)) ? 1 : 0;
    }
    else // NanoAOD case
    {
      if (!(static_cast<AliNanoAODTrack *>(track)->TestFilterBit(fFilterBit)))
        continue;
      if (fCheckTPCGeo)
      {
        static const Int_t tpcGeo_index = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCGeoLength");
        isTPCGeo = (static_cast<AliNanoAODTrack *>(track)->GetVar(tpcGeo_index) > 0.5) ? 1 : 0;
      }
    }
    // ---------- Track selection end ----------

    GetImpactParam(track, b, bCov);

    kaonZ = b[1];
    nTPCNSigKaon = GetTPCnSigma(track, AliPID::kKaon);
    nTOFNSigKaon = GetTOFnSigma(track, AliPID::kKaon);
    kaonpT = track->Pt();
    kaonSigmaDCA_r = (0.0026 + 0.0050 / kaonpT);
    kaonDCA_r = b[0];
    lEta = track->Eta();
    if (fFillQAPlot)
    {
      fHistos->FillTH1("QA/hDCAPVKaon", kaonZ);
      fHistos->FillTH1("QA/hDCArPVKaon", kaonDCA_r);
      fHistos->FillTH1("QA/hEtaKaon", lEta);
      fHistos->FillTH1("QA/hPtKaon", kaonpT);
      fHistos->FillTH2("QA/hTPCPIDKaon", kaonpT, nTPCNSigKaon);
      fHistos->FillTH2("QA/hTOFPIDKaon", kaonpT, nTOFNSigKaon);
      fHistos->FillTH2("QA/hTPCTOFPIDKaon", nTPCNSigKaon, nTOFNSigKaon);
      if (fCheckTPCGeo)
        fHistos->FillTH1("QA/hTPCGeoCheckKaon", isTPCGeo);
    }

    if (TMath::Abs(nTPCNSigKaon) > fTPCNsigXi1820ZeroKaonCut)
      continue;
    if (TMath::Abs(nTOFNSigKaon) > fTOFNsigXi1820ZeroKaonCut)
      continue;
    if (TMath::Abs(lEta) > fXi1820ZeroKaonEtaCut)
      continue;
    if (kaonpT < 0.15)
      continue;
    if (kaonZ > fXi1820ZeroKaonZVertexCut)
      continue;
    if (kaonDCA_r > kaonSigmaDCA_r * fXi1820KaonZeroXYVertexSigmaCut)
      continue;
    if (fCheckTPCGeo && isTPCGeo < 0.5)
      continue;

    if (fFillQAPlot)
    {
      fHistos->FillTH1("QAcut/hDCAPVKaon", kaonZ);
      fHistos->FillTH1("QAcut/hDCArPVKaon", kaonDCA_r);
      fHistos->FillTH1("QAcut/hEtaKaon", lEta);
      fHistos->FillTH1("QAcut/hPtKaon", kaonpT);
      fHistos->FillTH2("QAcut/hTPCPIDKaon", kaonpT, nTPCNSigKaon);
      fHistos->FillTH2("QAcut/hTOFPIDKaon", kaonpT, nTOFNSigKaon);
      fHistos->FillTH2("QAcut/hTPCTOFPIDKaon", nTPCNSigKaon, nTOFNSigKaon);
      if (fCheckTPCGeo)
        fHistos->FillTH1("QAcut/hTPCGeoCheckKaon", isTPCGeo);
    }

    fGoodKaonArray.push_back(it);
  }
  return fGoodKaonArray.size();
}
Bool_t AliAnalysisTaskXi1820BH::GoodK0sSelection()
{
  fGoodK0sArray.clear();
  const UInt_t nK0s = fEvt->GetNumberOfV0s();

  AliAODv0 *v0AOD;
  AliAODTrack *pTrackV0, *nTrackV0;
  Double_t lDCADist_K0sPiPlus_PV{0}, lDCADist_K0sPiMinus_PV{0}, lDCADistK0s{0},
      lDCADistK0s_PV{0}, lK0sCPA{0}, nTPCNSigKaonPlus{0}, nTPCNSigKaonMinus{0},
      fMassK0s{0}, lV0Radius{0}, lV0TotalMomentum{0}, lLength{0}, lLifetime{0},
      lK0sRap{0};
  Bool_t AcceptedV0 = kTRUE;

  for (UInt_t it = 0; it < nK0s; it++)
  {
    AcceptedV0 = kTRUE;
    v0AOD = ((AliAODEvent *)fEvt)->GetV0(it);
    if (!v0AOD)
      continue;

    if (!fOnlyUseOnTheFlyV0 && v0AOD->GetOnFlyStatus())
      continue;
    if (fOnlyUseOnTheFlyV0 && !v0AOD->GetOnFlyStatus())
      continue;

    if (TMath::Abs(v0AOD->GetPosID()) == TMath::Abs(v0AOD->GetNegID()))
      continue;

    pTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
    nTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

    // filter like-sign V0
    if (TMath::Abs(((pTrackV0->GetSign()) - (nTrackV0->GetSign()))) < 0.1)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // PID cuts
    nTPCNSigKaonMinus = GetTPCnSigma(nTrackV0, AliPID::kPion);
    nTPCNSigKaonPlus = GetTPCnSigma(pTrackV0, AliPID::kPion);

    if ((TMath::Abs(nTPCNSigKaonMinus) > fTPCNsigK0sPionCut) &&
        (TMath::Abs(nTPCNSigKaonPlus) > fTPCNsigK0sPionCut))
    {
      if (fFillQAPlot)
      {
        fHistos->FillTH2("QA/hTPCPIDK0sPionP", pTrackV0->Pt(), nTPCNSigKaonMinus);
        fHistos->FillTH2("QA/hTPCPIDK0sPionN", nTrackV0->Pt(), nTPCNSigKaonPlus);
        AcceptedV0 = kFALSE;
      }
      else
        continue;
    }

    // DCA cut
    // DCA between Dautgher particles
    lDCADistK0s = TMath::Abs(v0AOD->DcaV0Daughters());
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCA_k0sDaughters", lDCADistK0s);
    if (lDCADistK0s > fDCADistK0sDaughtersCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue; // DCA pion-pion

    // DCA to PV V0
    lDCADistK0s_PV = TMath::Abs(v0AOD->DcaV0ToPrimVertex());
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCAk0sPV", lDCADistK0s_PV);
    if (lDCADistK0s_PV > fDCArDistK0sPVCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // DCA to PV daughter
    // Adopted from AliRsnCutV0 fCustomTrackDCACuts
    lDCADist_K0sPiPlus_PV = v0AOD->DcaPosToPrimVertex();
    lDCADist_K0sPiMinus_PV = v0AOD->DcaNegToPrimVertex();
    if ((TMath::Abs(lDCADist_K0sPiPlus_PV) < fDCAPositiveTrackK0s) || (TMath::Abs(lDCADist_K0sPiMinus_PV) < fDCANegativeTrackK0s))
    {
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;
    }
    if (fFillQAPlot)
    {
      fHistos->FillTH1("QA/hDCAK0sPVPionP", lDCADist_K0sPiPlus_PV);
      fHistos->FillTH1("QA/hDCAK0sPVPionN", lDCADist_K0sPiMinus_PV);
    }

    // CPA cut
    lK0sCPA = TMath::Abs(v0AOD->CosPointingAngle(fVertex));
    if (lK0sCPA < fK0sCosineOfPointingAngleCut || lK0sCPA >= 1)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hCosPAk0s", lK0sCPA);

    // Rapidity cut
    lK0sRap = v0AOD->RapK0Short();
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hYk0s", lK0sRap);
    if (TMath::Abs(lK0sRap) > fMaxK0sRapidity)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Radius cut
    lV0Radius = v0AOD->RadiusV0();
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hK0sRxy", lV0Radius);
    if ((lV0Radius < fK0sLowRadius) || (lV0Radius > fK0sHighRadius))
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Life time cut
    lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
    lLength = v0AOD->DecayLength(fPosPV);
    lLifetime = TMath::Abs(k0sMass * lLength / lV0TotalMomentum);
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hLifetimek0s", lLifetime);
    if (lLifetime > fK0sLifetime)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Mass window cut
    fMassK0s = v0AOD->MassK0Short();
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hMassK0s", fMassK0s);
    if (TMath::Abs(fMassK0s - k0sMass) > fK0sMassWindowCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // After selection above
    if (AcceptedV0)
    {
      fGoodK0sArray.push_back(it); // for standard V0
      if (fFillQAPlot)
      {
        fHistos->FillTH2("QAcut/hTPCPIDK0sPionP", pTrackV0->Pt(), nTPCNSigKaonMinus);
        fHistos->FillTH2("QAcut/hTPCPIDK0sPionN", nTrackV0->Pt(), nTPCNSigKaonPlus);
        fHistos->FillTH1("QAcut/hDCA_k0sDaughters", lDCADistK0s);
        fHistos->FillTH1("QAcut/hDCA_k0sDaughters", lDCADistK0s);
        fHistos->FillTH1("QAcut/hDCAk0sPV", lDCADistK0s_PV);
        fHistos->FillTH1("QAcut/hDCAK0sPVPionP", lDCADist_K0sPiPlus_PV);
        fHistos->FillTH1("QAcut/hDCAK0sPVPionN", lDCADist_K0sPiMinus_PV);
        fHistos->FillTH1("QAcut/hDCAk0sPV", lDCADistK0s_PV);
        fHistos->FillTH1("QAcut/hLifetimek0s", lLifetime);
        fHistos->FillTH1("QAcut/hK0sRxy", lV0Radius);
        fHistos->FillTH1("QAcut/hCosPAk0s", lK0sCPA);
        fHistos->FillTH1("QAcut/hYk0s", lK0sRap);
        fHistos->FillTH1("QAcut/hMassK0s", fMassK0s);
      }
    }
  }
  return fGoodK0sArray.size();
}
Bool_t AliAnalysisTaskXi1820BH::GoodLambdaSelection()
{
  fGoodLambdaArray.clear();
  const UInt_t nLamba = fEvt->GetNumberOfV0s();

  AliAODv0 *v0AOD;
  AliAODTrack *pTrackV0, *nTrackV0;
  Bool_t lPIDLambda{kFALSE}, lPIDAntiLambda{kFALSE};
  Double_t lDCADist_LambdaProton_PV{0}, lDCADist_LambdaPion_PV{0},
      lDCADistLambda{0}, lDCADistLambda_PV{0}, lLambdaCPA{0}, nTPCNSigProton{0},
      nTPCNSigAntiProton{0}, nTPCNSigPion{0}, nTPCNSigAntiPion{0}, fMassV0{0}, radius{0},
      lV0TotalMomentum{0}, lLength{0}, lLifetime{0};
  UInt_t isAnti{0}, isAntiCheck{0};

  Bool_t AcceptedV0 = kTRUE;

  for (UInt_t it = 0; it < nLamba; it++)
  {
    lPIDLambda = kFALSE;
    lPIDAntiLambda = kFALSE;
    AcceptedV0 = kTRUE;
    isAntiCheck = 0;
    v0AOD = ((AliAODEvent *)fEvt)->GetV0(it);
    if (!v0AOD)
      continue;

    if (!fOnlyUseOnTheFlyV0 && v0AOD->GetOnFlyStatus())
      continue;
    if (fOnlyUseOnTheFlyV0 && !v0AOD->GetOnFlyStatus())
      continue;

    if (TMath::Abs(v0AOD->GetPosID()) == TMath::Abs(v0AOD->GetNegID()))
      continue;

    pTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
    nTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

    // filter like-sign V0
    if (TMath::Abs(((pTrackV0->GetSign()) - (nTrackV0->GetSign()))) < 0.1)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // PID cuts
    nTPCNSigProton = GetTPCnSigma(pTrackV0, AliPID::kProton);
    nTPCNSigAntiProton = GetTPCnSigma(nTrackV0, AliPID::kProton);
    nTPCNSigPion = GetTPCnSigma(nTrackV0, AliPID::kPion);
    nTPCNSigAntiPion = GetTPCnSigma(pTrackV0, AliPID::kPion);

    if ((TMath::Abs(nTPCNSigProton) > fTPCNsigLambdaProtonCut) &&
        (TMath::Abs(nTPCNSigPion) > fTPCNsigLambdaPionCut) &&
        (TMath::Abs(nTPCNSigAntiProton) > fTPCNsigLambdaProtonCut) &&
        (TMath::Abs(nTPCNSigAntiPion) > fTPCNsigLambdaPionCut))
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    if ((TMath::Abs(nTPCNSigProton) <= fTPCNsigLambdaProtonCut) &&
        (TMath::Abs(nTPCNSigPion) <= fTPCNsigLambdaPionCut))
    {
      lPIDLambda = kTRUE;
      isAnti = 0;
      if (fFillQAPlot)
      {
        fHistos->FillTH2("QA/hTPCPIDLambdaProton", pTrackV0->GetTPCmomentum(),
                         pTrackV0->GetTPCsignal());
        fHistos->FillTH2("QA/hTPCPIDLambdaPion", nTrackV0->GetTPCmomentum(),
                         nTrackV0->GetTPCsignal());
      }
    }
    else if ((TMath::Abs(nTPCNSigAntiProton) <= fTPCNsigLambdaProtonCut) &&
             (TMath::Abs(nTPCNSigAntiPion) <= fTPCNsigLambdaPionCut))
    {
      lPIDAntiLambda = kTRUE;
      isAnti = 1;
      if (fFillQAPlot)
      {
        fHistos->FillTH2("QA/hTPCPIDAntiLambdaProton",
                         nTrackV0->GetTPCmomentum(),
                         nTrackV0->GetTPCsignal());
        fHistos->FillTH2("QA/hTPCPIDAntiLambdaPion",
                         pTrackV0->GetTPCmomentum(),
                         pTrackV0->GetTPCsignal());
      }
    }
    else if (fFillQAPlot)
      AcceptedV0 = kFALSE;
    else
      continue;
    if (lPIDLambda)
      isAntiCheck += 1;
    if (lPIDAntiLambda)
      isAntiCheck += 2;
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hLambdaAntiCheck", isAntiCheck);

    // DCA cut
    lDCADistLambda = TMath::Abs(v0AOD->DcaV0Daughters());
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCA_lambdaDaughters", lDCADistLambda);
    if (lDCADistLambda > fDCADistLambdaDaughtersCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue; // DCA proton-pion

    // DCA to PV daughter
    // Adopted from AliRsnCutV0 fCustomTrackDCACuts
    if (!isAnti)
    {
      lDCADist_LambdaProton_PV = v0AOD->DcaPosToPrimVertex();
      lDCADist_LambdaPion_PV = v0AOD->DcaNegToPrimVertex();
    }
    else
    {
      lDCADist_LambdaPion_PV = v0AOD->DcaPosToPrimVertex();
      lDCADist_LambdaProton_PV = v0AOD->DcaNegToPrimVertex();
    }
    if (TMath::Abs(lDCADist_LambdaProton_PV) < fDCAPositiveTrackLambda)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;
    if (TMath::Abs(lDCADist_LambdaPion_PV) < fDCANegativeTrackLambda)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCALambdaPVProton", lDCADist_LambdaProton_PV);
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCALambdaPVPion", lDCADist_LambdaPion_PV);

    // DCA to PV
    lDCADistLambda_PV = TMath::Abs(v0AOD->DcaV0ToPrimVertex());
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hDCAlambdaPV", lDCADistLambda_PV);
    if (lDCADistLambda_PV > fDCArDistLambdaPVCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // CPA cut
    lLambdaCPA = TMath::Abs(v0AOD->CosPointingAngle(fVertex));
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hCosPAlambda", lLambdaCPA);
    if (lLambdaCPA < fLambdaCosineOfPointingAngleCut || lLambdaCPA >= 1)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Rapidity cut
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hYlambda", v0AOD->RapLambda());
    if (TMath::Abs(v0AOD->RapLambda()) > fMaxLambdaRapidity)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Radius cut
    radius = v0AOD->RadiusV0();
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hLambdaRxy", radius);
    if ((radius < fLambdaLowRadius) || (radius > fLambdaHighRadius))
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Life time cut
    lV0TotalMomentum = TMath::Sqrt(v0AOD->Ptot2V0());
    lLength = v0AOD->DecayLength(fPosPV);
    lLifetime = TMath::Abs(v0Mass * lLength / lV0TotalMomentum);
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hLifetimelambda", lLifetime);
    if (lLifetime > fLambdaLifetime)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // Mass window cut
    fMassV0 = -999;
    if (lPIDLambda)
      fMassV0 = v0AOD->MassLambda();
    else if (lPIDAntiLambda)
      fMassV0 = v0AOD->MassAntiLambda();
    if (fFillQAPlot)
      fHistos->FillTH1("QA/hMassLambda", fMassV0);
    if (TMath::Abs(fMassV0 - v0Mass) > fLambdaMassWindowCut)
      if (fFillQAPlot)
        AcceptedV0 = kFALSE;
      else
        continue;

    // After selection above
    if (AcceptedV0)
    {
      fGoodLambdaArray.push_back({it, isAnti});
      if (fFillQAPlot)
      {
        if (lPIDLambda)
        {
          fHistos->FillTH2("QAcut/hTPCPIDLambdaProton", pTrackV0->GetTPCmomentum(), pTrackV0->GetTPCsignal());
          fHistos->FillTH2("QAcut/hTPCPIDLambdaPion", nTrackV0->GetTPCmomentum(), nTrackV0->GetTPCsignal());
        }
        if (lPIDAntiLambda)
        {
          fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaProton", nTrackV0->GetTPCmomentum(), nTrackV0->GetTPCsignal());
          fHistos->FillTH2("QAcut/hTPCPIDAntiLambdaPion", pTrackV0->GetTPCmomentum(), pTrackV0->GetTPCsignal());
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
  } // All v0 loop

  return fGoodLambdaArray.size();
}
void AliAnalysisTaskXi1820BH::FillKaonV0()
{
  AliVTrack *track1;
  AliVTrack *track_mix;
  AliAODv0 *v0AOD;
  AliAODTrack *pTrackV0, *nTrackV0;
  Bool_t isAnti, isKaonPlus, isTrueXi1820;
  Bool_t SkipMixing = kFALSE;
  Int_t pID, nID;
  int sign = kAllType;
  int binAnti = 0;
  double asym = 0.;

  TLorentzVector vecLambda, vecKaon, vecKaonMix;
  TLorentzVector vecXi1820PM;
  const UInt_t nLamba = fGoodLambdaArray.size();
  const UInt_t nTracks = fGoodKaonArray.size();

  // Prepare mixing pool
  tracklist trackpool;
  if (fSetMixing)
  {
    // Kaon track pool
    eventpool &ep = fEMpoolKaon[fCentBin][fZbin];
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
  // Lambda loop
  for (UInt_t i = 0; i < nLamba; i++)
  {

    v0AOD = ((AliAODEvent *)fEvt)->GetV0(fGoodLambdaArray[i][0]);
    if (!v0AOD)
      continue;
    pTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
    nTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(1));
    pID = pTrackV0->GetID();
    nID = nTrackV0->GetID();

    if (fGoodLambdaArray[i][1] > 0)
      isAnti = true;
    else
      isAnti = false;

    if (!isAnti)
      vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(), v0AOD->MassLambda());
    else
      vecLambda.SetXYZM(v0AOD->MomV0X(), v0AOD->MomV0Y(), v0AOD->MomV0Z(), v0AOD->MassAntiLambda());

    // Kaon loop
    for (UInt_t j = 0; j < nTracks; j++)
    {
      track1 = (AliVTrack *)fEvt->GetTrack(fGoodKaonArray[j]);
      if (!track1)
        continue;

      // Skip duplicated track
      if (track1->GetID() == pID || track1->GetID() == nID)
        continue;
      vecKaon.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), kaonMass);
      vecXi1820PM = vecLambda + vecKaon;
      // Y cut
      if ((vecXi1820PM.Rapidity() > fXi1820YCutHigh) ||
          (vecXi1820PM.Rapidity() < fXi1820YCutLow))
        continue;

      // AsymmCut
      if (fUseAsymmCut)
      {
        auto P1 = (TVector3)vecLambda.Vect();
        auto P2 = (TVector3)vecKaon.Vect();
        asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if (fFillQAPlot)
          fHistos->FillTH2("QA/hXi1820Asymm", asym, fCent);
        if ((asym < fXi1820AsymmCutLow) || (asym > fXi1820AsymmCutHigh))
          continue;
        if (fFillQAPlot)
          fHistos->FillTH2("QAcut/hXi1820Asymm", asym, fCent);
      }

      if (track1->Charge() > 0)
        isKaonPlus = true;
      else
        isKaonPlus = false;

      if (isKaonPlus && !isAnti) // Kaon+ + Lambda -> skip
        continue;
      if (!isKaonPlus && isAnti) // Kaon- + Anti-Lambda -> skip
        continue;

      binAnti = (isAnti) ? kAnti : kNormal; // Anti
      sign = kXi1820M;

      FillTHnSparse("Xi1820_data", {(double)binAnti, (double)sign, (double)fCent, vecXi1820PM.Pt(), vecXi1820PM.M()});

      if (fIsMC)
      {
        isTrueXi1820 = IsTrueXi1820PM(fGoodLambdaArray[i][0], fGoodKaonArray[j], 0);
        if (isTrueXi1820)
        {
          sign = kXi1820M_REC;
          FillTHnSparse("Xi1820_mc", {(double)binAnti, (double)sign, (double)fCent, vecXi1820PM.Pt(), vecXi1820PM.M()});
          if (fUseAsymmCut && fFillQAPlot)
            fHistos->FillTH2("QAcut/hXi1820Asymm_true_selected", asym, fCent);
        }
        else
        {
          // MC not true bkg
          sign = kXi1820M_NOT;
          FillTHnSparse("Xi1820_data", {(double)binAnti, (double)sign, (double)fCent, vecXi1820PM.Pt(), vecXi1820PM.M()});
        }
      }
    } // Kaon loop

    if (fSetMixing && !SkipMixing && (fCentBin >= 0) && (fZbin >= 0))
    {
      for (UInt_t jt = 0; jt < trackpool.size(); jt++)
      {
        track_mix = trackpool.at(jt);
        if (track_mix->GetID() == pID || track_mix->GetID() == nID)
          continue;
        vecKaonMix.SetXYZM(track_mix->Px(), track_mix->Py(), track_mix->Pz(), kaonMass);
        vecXi1820PM = vecLambda + vecKaonMix;
        // Y cut
        if ((vecXi1820PM.Rapidity() > fXi1820YCutHigh) ||
            (vecXi1820PM.Rapidity() < fXi1820YCutLow))
          continue;
        // AsymmCut
        if (fUseAsymmCut)
        {
          auto P1 = (TVector3)vecLambda.Vect();
          auto P2 = (TVector3)vecKaonMix.Vect();
          auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
          if ((asym < fXi1820AsymmCutLow) || (asym > fXi1820AsymmCutHigh))
            continue;
        }

        if (track_mix->Charge() > 0)
          isKaonPlus = true;
        else
          isKaonPlus = false;

        if (fGoodLambdaArray[i][1] > 0)
          isAnti = true;
        else
          isAnti = false;
        if (isKaonPlus && !isAnti)
          continue;
        if (!isKaonPlus && isAnti)
          continue;

        binAnti = (isAnti) ? kAnti : kNormal;
        sign = kXi1820M_MIX;

        FillTHnSparse("Xi1820_data", {(double)binAnti, (double)sign, (double)fCent, vecXi1820PM.Pt(), vecXi1820PM.M()});
      }
    }
  }
}
void AliAnalysisTaskXi1820BH::FillK0sV0()
{
  AliAODv0 *k0s_V0, *k0s_V0_mix, *LambdaV0;
  AliAODTrack *pTrackV0, *nTrackV0, *pTrackK0s, *nTrackK0s, *pTrackK0s_Mix, *nTrackK0s_Mix;
  Bool_t isAnti;
  Bool_t isTrueXi1820;
  Bool_t SkipMixing = kFALSE;
  Int_t pID, nID, pIDK0s, nIDK0s, pIDK0s_Mix, nIDK0s_Mix;
  int binAnti = 0;
  double asym = 0.;

  TLorentzVector vecLambda, vecK0s, vecK0sMix;
  TLorentzVector vecXi1820Zero;
  const UInt_t nLamba = fGoodLambdaArray.size();
  const UInt_t nK0s = fGoodK0sArray.size();
  const UInt_t nTracks = fGoodKaonArray.size();

  // Prepare mixing pool
  v0list k0sMixingPool;
  if (fSetMixing)
  {
    // K0s track pool
    v0pool &ep = fEMpoolK0s[fCentBin][fZbin];
    if ((int)ep.size() < (int)fnMix)
      SkipMixing = kTRUE;
    if (!SkipMixing)
    {
      for (auto pool : ep)
      {
        for (auto v0 : pool)
          k0sMixingPool.push_back((AliAODv0 *)v0);
      }
    }
  }

  for (UInt_t i = 0; i < nLamba; i++) // Lambda loop
  {
    // std::cout << "Lambda loop: " << i << std::endl;
    LambdaV0 = ((AliAODEvent *)fEvt)->GetV0(fGoodLambdaArray[i][0]);
    if (!LambdaV0)
      continue;
    pTrackV0 = (AliAODTrack *)(LambdaV0->GetSecondaryVtx()->GetDaughter(0));
    nTrackV0 = (AliAODTrack *)(LambdaV0->GetSecondaryVtx()->GetDaughter(1));
    pID = pTrackV0->GetID();
    nID = nTrackV0->GetID();

    if (fGoodLambdaArray[i][1] > 0)
      isAnti = true;
    else
      isAnti = false;

    if (!isAnti)
      vecLambda.SetXYZM(LambdaV0->MomV0X(), LambdaV0->MomV0Y(), LambdaV0->MomV0Z(), LambdaV0->MassLambda());
    else
      vecLambda.SetXYZM(LambdaV0->MomV0X(), LambdaV0->MomV0Y(), LambdaV0->MomV0Z(), LambdaV0->MassAntiLambda());

    for (UInt_t j = 0; j < nK0s; j++) // Kaon loop
    {
      // std::cout << "K0s loop: " << j << std::endl;
      k0s_V0 = ((AliAODEvent *)fEvt)->GetV0(fGoodK0sArray[j]);
      if (!k0s_V0)
      {
        continue;
      }
      if (fGoodK0sArray[j] == fGoodLambdaArray[i][0])
      {
        continue;
      }

      pTrackK0s = (AliAODTrack *)(k0s_V0->GetSecondaryVtx()->GetDaughter(0));
      nTrackK0s = (AliAODTrack *)(k0s_V0->GetSecondaryVtx()->GetDaughter(1));
      pIDK0s = pTrackK0s->GetID();
      nIDK0s = nTrackK0s->GetID();

      // Skip duplicated track
      if (pIDK0s == pID || pIDK0s == nID || nIDK0s == pID || nIDK0s == nID)
      {
        continue;
      }

      vecK0s.SetXYZM(k0s_V0->Px(), k0s_V0->Py(), k0s_V0->Pz(), k0sMass);
      vecXi1820Zero = vecLambda + vecK0s;
      // Y cut
      if ((vecXi1820Zero.Rapidity() > fXi1820YCutHigh) ||
          (vecXi1820Zero.Rapidity() < fXi1820YCutLow))
      {
        continue;
      }

      // AsymmCut
      if (fUseAsymmCut)
      {
        auto P1 = (TVector3)vecLambda.Vect();
        auto P2 = (TVector3)vecK0s.Vect();
        asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
        if (fFillQAPlot)
          fHistos->FillTH2("QA/hXi1820Asymm", asym, fCent);
        if ((asym < fXi1820AsymmCutLow) || (asym > fXi1820AsymmCutHigh))
          continue;
        if (fFillQAPlot)
          fHistos->FillTH2("QAcut/hXi1820Asymm", asym, fCent);
      }
      binAnti = (isAnti) ? kAnti : kNormal;
      FillTHnSparse("Xi1820_data", {(double)binAnti, (double)kXi1820Zero, (double)fCent, vecXi1820Zero.Pt(), vecXi1820Zero.M()});

      if (fIsMC)
      {
        isTrueXi1820 = IsTrueXi1820Zero(fGoodLambdaArray[i][0], fGoodK0sArray[j], 0);
        if (isTrueXi1820)
        {
          FillTHnSparse("Xi1820_mc", {(double)binAnti, (double)kXi1820Zero_REC, (double)fCent, vecXi1820Zero.Pt(), vecXi1820Zero.M()});
          if (fUseAsymmCut && fFillQAPlot)
            fHistos->FillTH2("QAcut/hXi1820Asymm_true_selected", asym, fCent);
        }
        else // MC not true bkg
          FillTHnSparse("Xi1820_data", {(double)binAnti, (double)kXi1820Zero_NOT, (double)fCent, vecXi1820Zero.Pt(), vecXi1820Zero.M()});
      }
    } // K0s loop

    if (fSetMixing && !SkipMixing && (fCentBin >= 0) && (fZbin >= 0))
    {
      for (UInt_t jt = 0; jt < k0sMixingPool.size(); jt++)
      {
        k0s_V0_mix = k0sMixingPool.at(jt);
        pTrackK0s_Mix = (AliAODTrack *)(k0s_V0->GetSecondaryVtx()->GetDaughter(0));
        nTrackK0s_Mix = (AliAODTrack *)(k0s_V0->GetSecondaryVtx()->GetDaughter(1));
        pIDK0s_Mix = pTrackK0s_Mix->GetID();
        nIDK0s_Mix = nTrackK0s_Mix->GetID();

        // Skip duplicated track
        if (pIDK0s_Mix == pID || pIDK0s_Mix == nID || nIDK0s_Mix == pID || nIDK0s_Mix == nID)
          continue;
        vecK0sMix.SetXYZM(k0s_V0_mix->Px(), k0s_V0_mix->Py(), k0s_V0_mix->Pz(), k0sMass);
        vecXi1820Zero = vecLambda + vecK0sMix;
        // Y cut
        if ((vecXi1820Zero.Rapidity() > fXi1820YCutHigh) ||
            (vecXi1820Zero.Rapidity() < fXi1820YCutLow))
          continue;
        // AsymmCut
        if (fUseAsymmCut)
        {
          auto P1 = (TVector3)vecLambda.Vect();
          auto P2 = (TVector3)vecK0sMix.Vect();
          auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
          if ((asym < fXi1820AsymmCutLow) || (asym > fXi1820AsymmCutHigh))
            continue;
        }

        if (fGoodLambdaArray[i][1] > 0)
          isAnti = true;
        else
          isAnti = false;

        binAnti = (isAnti) ? kAnti : kNormal;
        FillTHnSparse("Xi1820_data", {(double)binAnti, (double)kXi1820Zero_MIX, (double)fCent, vecXi1820Zero.Pt(), vecXi1820Zero.M()});
      }
    }
  }
}
void AliAnalysisTaskXi1820BH::FillMCEventProperties(AliMCEvent *fMCEvent)
{
  // Fill MC event properties to calculate Trigger efficiency and vertex selection efficiency
  bool IsINEL0True = (fIsNano) ? true : fEventCuts.IsTrueINELgtZero(fEvt, true);
  bool IsVtxInZCut = (fIsNano) ? true : fEventCuts.PassedCut(AliEventCuts::kVertexPosition);
  bool IsSelectedTrig = (fIsNano) ? true : fEventCuts.PassedCut(AliEventCuts::kTrigger);
  bool IsGoodVertex = (fIsNano) ? true : fEventCuts.PassedCut(AliEventCuts::kVertexQuality);
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
        FillTHnSparse("Normalisation",
                      {(int)kINELg0_trig_vtx, (double)fCent});
        if (IsVtxInZCut)
        {
          FillTHnSparse("Normalisation",
                        {(int)kINELg0_trig_vtx10, (double)fCent});
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
        FillTHnSparse("Normalisation",
                      {(int)kINEL_trig_vtx10, (double)fCent});
      }
    }
  }
}
void AliAnalysisTaskXi1820BH::FillMCinput(AliMCEvent *fMCEvent, int Fillbin)
{
  int sign = kAllType;
  int binAnti = 0;
  TLorentzVector vecPart1, vecPart2;

  for (Int_t it = 0; it < fMCArray->GetEntriesFast(); it++)
  {
    AliAODMCParticle *mcInputTrack = (AliAODMCParticle *)fMCArray->At(it);
    if (!mcInputTrack)
    {
      Error("UserExec", "Could not receive MC track %d", it);
      continue;
    }

    Int_t v0PdgCode = mcInputTrack->GetPdgCode();

    if ((TMath::Abs(v0PdgCode) != kXi1820ZeroCode) &&
        (TMath::Abs(v0PdgCode) != kXi1820MCode))
      continue;
    if (fIsPrimaryMC && !mcInputTrack->IsPrimary())
      continue;

    // Y cut
    if ((mcInputTrack->Y() > fXi1820YCutHigh) ||
        (mcInputTrack->Y() < fXi1820YCutLow))
      continue;
    // AsymmCut
    if (fUseAsymmCut)
    {
      auto mcPart1 = (AliAODMCParticle *)fMCArray->At(
          abs(mcInputTrack->GetDaughterFirst()));
      auto mcPart2 = (AliAODMCParticle *)fMCArray->At(
          abs(mcInputTrack->GetDaughterLast()));
      mcPart1->Momentum(vecPart1);
      mcPart2->Momentum(vecPart2);

      auto P1 = (TVector3)vecPart1.Vect();
      auto P2 = (TVector3)vecPart2.Vect();
      auto asym = TMath::Abs(P1.Mag() - P2.Mag()) / (P1.Mag() + P2.Mag());
      if (fFillQAPlot)
        fHistos->FillTH2("QA/hXi1820Asymm_true", asym, fCent);
      if ((asym < fXi1820AsymmCutLow) || (asym > fXi1820AsymmCutHigh))
        continue;
      if (fFillQAPlot)
        fHistos->FillTH2("QAcut/hXi1820Asymm_true", asym, fCent);
    }

    binAnti = (v0PdgCode < 0) ? kAnti : kNormal;
    if (TMath::Abs(v0PdgCode) == kXi1820ZeroCode)
    {
      sign = kXi1820Zero_GEN + (int)Fillbin * 2;
      if (!fSetZero)
        continue;
    }
    if (TMath::Abs(v0PdgCode) == kXi1820MCode)
    {
      sign = kXi1820M_GEN + (int)Fillbin * 2;
    }

    FillTHnSparse("Xi1820_mc",
                  {(double)binAnti, (double)sign, (double)fCent,
                   mcInputTrack->Pt(), mcInputTrack->GetCalcMass()});
  }
}
Bool_t AliAnalysisTaskXi1820BH::IsTrueXi1820PM(UInt_t v0Index, UInt_t pionIndex, UInt_t BkgCheck)
{
  AliVTrack *track1 = (AliVTrack *)fEvt->GetTrack(pionIndex);
  AliAODv0 *v0AOD = ((AliAODEvent *)fEvt)->GetV0(v0Index);
  if (!v0AOD)
    return kFALSE;
  AliAODTrack *pTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
  AliAODTrack *nTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

  AliAODMCParticle *MCLamD1 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(pTrackV0->GetLabel()));
  AliAODMCParticle *MCLamD2 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(nTrackV0->GetLabel()));
  AliAODMCParticle *MCLam;
  AliAODMCParticle *MCPi;
  AliAODMCParticle *MCXi1820PM;

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
  MCLam = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCLamD1->GetMother()));
  if (TMath::Abs(MCLam->GetPdgCode()) != kLambdaCode)
    return kFALSE;
  // Kaon check
  MCPi = (AliAODMCParticle *)fMCArray->At(TMath::Abs(track1->GetLabel()));
  if (TMath::Abs(MCPi->GetPdgCode()) != kKaonCode)
    return kFALSE;

  switch (BkgCheck)
  {
  case 0: // Normal Xi1820PM case
    // pion-Lambda mother check
    if (MCPi->GetMother() != MCLam->GetMother())
      return kFALSE;
    MCXi1820PM = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCLam->GetMother()));
    if (TMath::Abs(MCXi1820PM->GetPdgCode()) != kXi1820MCode)
      return kFALSE;
    if (fIsPrimaryMC)
    {
      if (MCXi1820PM->IsPrimary())
        return kTRUE;
      else
        return kFALSE;
    }
    else
      return kTRUE;
    break;
  }
}
Bool_t AliAnalysisTaskXi1820BH::IsTrueXi1820Zero(UInt_t v0Index, UInt_t k0sIndex, UInt_t BkgCheck)
{
  AliAODv0 *k0sAOD = ((AliAODEvent *)fEvt)->GetV0(k0sIndex);
  AliAODv0 *v0AOD = ((AliAODEvent *)fEvt)->GetV0(v0Index);
  if (!v0AOD)
    return kFALSE;
  if (!k0sAOD)
    return kFALSE;
  AliAODTrack *pTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(0));
  AliAODTrack *nTrackV0 = (AliAODTrack *)(v0AOD->GetSecondaryVtx()->GetDaughter(1));

  AliAODTrack *pTrackK0s = (AliAODTrack *)(k0sAOD->GetSecondaryVtx()->GetDaughter(0));
  AliAODTrack *nTrackK0s = (AliAODTrack *)(k0sAOD->GetSecondaryVtx()->GetDaughter(1));

  AliAODMCParticle *MCLamD1 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(pTrackV0->GetLabel()));
  AliAODMCParticle *MCLamD2 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(nTrackV0->GetLabel()));
  AliAODMCParticle *MCK0sD1 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(pTrackK0s->GetLabel()));
  AliAODMCParticle *MCK0sD2 = (AliAODMCParticle *)fMCArray->At(TMath::Abs(nTrackK0s->GetLabel()));
  AliAODMCParticle *MCLam;
  AliAODMCParticle *MCK0s;
  AliAODMCParticle *MCXi1820Zero;

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
  MCLam = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCLamD1->GetMother()));
  if (TMath::Abs(MCLam->GetPdgCode()) != kLambdaCode)
    return kFALSE;

  // K0s daughter (pion, pion) check
  if ((TMath::Abs(MCK0sD1->GetPdgCode()) != kPionCode ||
       TMath::Abs(MCK0sD2->GetPdgCode()) != kPionCode))
    return kFALSE;
  // K0s duather's mother check
  if (MCK0sD1->GetMother() != MCK0sD2->GetMother())
    return kFALSE;
  MCK0s = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCK0sD1->GetMother()));
  if (TMath::Abs(MCK0s->GetPdgCode()) != kK0sCode)
    return kFALSE;

  switch (BkgCheck)
  {
  case 0: // Normal Xi1820PM case
    // K0s-Lambda mother check
    if (MCK0s->GetMother() != MCLam->GetMother())
      return kFALSE;
    std::cout << "same mother!" << std::endl;
    MCXi1820Zero = (AliAODMCParticle *)fMCArray->At(TMath::Abs(MCLam->GetMother()));
    if (TMath::Abs(MCXi1820Zero->GetPdgCode()) != kXi1820ZeroCode)
      return kFALSE;
    std::cout << "PDG code: " << MCXi1820Zero->GetPdgCode() << std::endl;
    if (fIsPrimaryMC)
    {
      if (MCXi1820Zero->IsPrimary())
        return kTRUE;
      else
        return kFALSE;
    }
    else
      return kTRUE;
    break;
  }
}

THnSparse *AliAnalysisTaskXi1820BH::CreateTHnSparse(
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
Long64_t AliAnalysisTaskXi1820BH::FillTHnSparse(TString name,
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

Long64_t AliAnalysisTaskXi1820BH::FillTHnSparse(THnSparse *h,
                                                std::vector<Double_t> x,
                                                Double_t w)
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
TAxis AliAnalysisTaskXi1820BH::AxisFix(TString name,
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
TAxis AliAnalysisTaskXi1820BH::AxisStr(TString name,
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

TAxis AliAnalysisTaskXi1820BH::AxisVar(TString name,
                                       std::vector<Double_t> bin)
{
  // From AliPhysics/PWGUD/DIFFRACTIVE/Resonance/AliAnalysisTaskf0f2.cxx
  // Original author: Beomkyu Kim
  TAxis axis(bin.size() - 1, &bin.front());
  axis.SetName(name);
  return axis;
}
double AliAnalysisTaskXi1820BH::GetTPCnSigma(AliVTrack *track,
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
double AliAnalysisTaskXi1820BH::GetTOFnSigma(AliVTrack *track,
                                             AliPID::EParticleType type)
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

void AliAnalysisTaskXi1820BH::FillKaonToEventPool()
{
  // Fill Selected tracks to event mixing pool
  if ((fCentBin < 0) || (fZbin < 0))
    return;
  AliVTrack *goodtrack;

  tracklist *etl;
  eventpool *ep;
  // Event mixing pool

  ep = &fEMpoolKaon[fCentBin][fZbin];
  ep->push_back(tracklist());
  etl = &(ep->back());
  // Fill selected tracks
  for (UInt_t i = 0; i < fGoodKaonArray.size(); i++)
  {
    goodtrack = (AliVTrack *)fEvt->GetTrack(fGoodKaonArray[i]);
    if (!goodtrack)
      continue;
    etl->push_back((AliVTrack *)goodtrack->Clone());
  }
  if (!fGoodKaonArray.size())
    ep->pop_back();
  if ((int)ep->size() > (int)fnMix)
  {
    for (auto it : ep->front())
      delete it;
    ep->pop_front();
  }
}
void AliAnalysisTaskXi1820BH::FillK0sToEventPool()
{
  // Fill Selected k0s to event mixing pool
  if ((fCentBin < 0) || (fZbin < 0))
    return;
  AliAODv0 *goodK0s;

  v0list *etl;
  v0pool *ep;
  // Event mixing pool

  ep = &fEMpoolK0s[fCentBin][fZbin];
  ep->push_back(v0list());
  etl = &(ep->back());
  // Fill selected k0s
  for (UInt_t i = 0; i < fGoodK0sArray.size(); i++)
  {
    goodK0s = (AliAODv0 *)fEvt->GetTrack(fGoodK0sArray[i]);
    if (!goodK0s)
      continue;
    etl->push_back((AliAODv0 *)goodK0s->Clone());
  }
  if (!fGoodK0sArray.size())
    ep->pop_back();
  if ((int)ep->size() > (int)fnMix)
  {
    for (auto it : ep->front())
      delete it;
    ep->pop_front();
  }
}
void AliAnalysisTaskXi1820BH::GetImpactParam(AliVTrack *track,
                                             Float_t p[2],
                                             Float_t cov[3])
{
  AliNanoAODTrack *nanoT = dynamic_cast<AliNanoAODTrack *>(track);
  if (nanoT)
    nanoT->AliNanoAODTrack::GetImpactParameters(p[0], p[1]);
  else
    track->GetImpactParameters(p, cov);
}
Bool_t AliAnalysisTaskXi1820BH::IsSelectedTPCGeoCut(AliAODTrack *track)
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
Bool_t AliAnalysisTaskXi1820BH::IsSelectedTPCGeoCut(AliESDtrack *track)
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