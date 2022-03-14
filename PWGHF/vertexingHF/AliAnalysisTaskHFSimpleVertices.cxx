#include <chrono>
#include <ctime>
#include <ratio>
#include <regex>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliNeutralTrackParam.h"
#include "AliESDtrackCuts.h"
#include "AliVertexerTracks.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexingHFUtils.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "DCAFitterN.h"
#include "AliAnalysisTaskHFSimpleVertices.h"
#include "AliLog.h"

#include <TH1F.h>
#include <TSystem.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TObjString.h>


/**************************************************************************
 * Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Implementation of class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//
//*************************************************************************

ClassImp(AliAnalysisTaskHFSimpleVertices)
//______________________________________________________________________________
AliAnalysisTaskHFSimpleVertices::AliAnalysisTaskHFSimpleVertices() :
  AliAnalysisTaskSE("HFSimpleVertices"),
  fOutput{nullptr},
  fHistNEvents{nullptr},
  fHistTrackStatus{nullptr},
  fHistPtAllTracks{nullptr},
  fHistPtSelTracks{nullptr},
  fHistTglAllTracks{nullptr},
  fHistTglSelTracks{nullptr},
  fHistEtaAllTracks{nullptr},
  fHistEtaSelTracks2prong{nullptr},
  fHistEtaSelTracks3prong{nullptr},
  fHistImpParAllTracks{nullptr},
  fHistImpParSelTracks2prong{nullptr},
  fHistImpParSelTracks3prong{nullptr},
  fHistITSmapAllTracks{nullptr},
  fHistITSmapSelTracks{nullptr},
  fHistPrimVertX{nullptr},
  fHistPrimVertY{nullptr},
  fHistPrimVertZ{nullptr},
  fHist2ProngVertX{nullptr},
  fHist2ProngVertY{nullptr},
  fHist2ProngVertZ{nullptr},
  fHist3ProngVertX{nullptr},
  fHist3ProngVertY{nullptr},
  fHist3ProngVertZ{nullptr},
  fHistDist12LcpKpi{nullptr},
  fHistInvMassD0{nullptr},
  fHistPtD0{nullptr},
  fHistYPtD0{nullptr},
  fHistPtD0Dau0{nullptr},
  fHistPtD0Dau1{nullptr},
  fHistImpParD0Dau0{nullptr},
  fHistImpParD0Dau1{nullptr},
  fHistd0Timesd0{nullptr},
  fHistCosPointD0{nullptr},
  fHistDecLenD0{nullptr},
  fHistDecLenXYD0{nullptr},
  fHistImpParErrD0Dau{nullptr},
  fHistDecLenErrD0{nullptr},
  fHistDecLenXYErrD0{nullptr},
  fHistCovMatPrimVXX2Prong{nullptr},
  fHistCovMatSecVXX2Prong{nullptr},
  fHistD0SignalVertX{nullptr},
  fHistD0SignalVertY{nullptr},
  fHistD0SignalVertZ{nullptr},
  fHistInvMassD0Signal{nullptr},
  fHistInvMassD0Refl{nullptr},
  fHistInvMassJpsi{nullptr},
  fHistPtJpsi{nullptr},
  fHistPtJpsiDau0{nullptr},
  fHistPtJpsiDau1{nullptr},
  fHistImpParJpsiDau0{nullptr},
  fHistImpParJpsiDau1{nullptr},
  fHistd0Timesd0Jpsi{nullptr},
  fHistCosPointJpsi{nullptr},
  fHistDecLenJpsi{nullptr},
  fHistDecLenXYJpsi{nullptr},
  fHistDecLenErrJpsi{nullptr},
  fHistDecLenXYErrJpsi{nullptr},
  fHistJpsiSignalVertX{nullptr},
  fHistJpsiSignalVertY{nullptr},
  fHistJpsiSignalVertZ{nullptr},
  fHistInvMassJpsiSignal{nullptr},
  fHistInvMassDplus{nullptr},
  fHistInvMassDplusSignal{nullptr},
  fHistPtDplus{nullptr},
  fHistYPtDplus{nullptr},
  fHistPtDplusDau0{nullptr},
  fHistPtDplusDau1{nullptr},
  fHistPtDplusDau2{nullptr},
  fHistImpParDplusDau0{nullptr},
  fHistImpParDplusDau1{nullptr},
  fHistImpParDplusDau2{nullptr},
  fHistDecLenDplus{nullptr},
  fHistDecLenXYDplus{nullptr},
  fHistNormDecLenXYDplus{nullptr},
  fHistImpParErrDplusDau{nullptr},
  fHistDecLenErrDplus{nullptr},
  fHistDecLenXYErrDplus{nullptr},
  fHistCosPointDplus{nullptr},
  fHistCosPointXYDplus{nullptr},
  fHistImpParXYDplus{nullptr},
  fHistNormIPDplus{nullptr},
  fHistoSumSqImpParDplusDau{nullptr},
  fHistCovMatPrimVXX3Prong{nullptr},
  fHistCovMatSecVXX3Prong{nullptr},
  fHistInvMassDs{nullptr},
  fHistInvMassDsSignal{nullptr},
  fHistInvMassDsRefl{nullptr},
  fHistPtDs{nullptr},
  fHistYPtDs{nullptr},
  fHistDecLenDs{nullptr},
  fHistCosPointDs{nullptr},
  fHistInvMassLc{nullptr},
  fHistPtLc{nullptr},
  fHistYPtLc{nullptr},
  fHistPtLcDau0{nullptr},
  fHistPtLcDau1{nullptr},
  fHistPtLcDau2{nullptr},
  fHistDecLenLc{nullptr},
  fHistCosPointLc{nullptr},
  fHistInvMassK0s{nullptr},
  fHistInvMassLcK0sp{nullptr},
  fHistPtLcK0sp{nullptr},
  fHistCPUTimeTrackVsNTracks{nullptr},
  fHistCPUTimeCandVsNTracks{nullptr},
  fHistWallTimeTrackVsNTracks{nullptr},
  fHistWallTimeCandVsNTracks{nullptr},
  fReadMC(kFALSE),
  fUsePhysSel(kTRUE),
  fTriggerMask(AliVEvent::kAny),
  fSelectOnCentrality(kFALSE),
  fMinCentrality(-1.),
  fMaxCentrality(110.),
  fCentrEstimator("V0M"),
  fCutOnSPDVsTrackVtx(kFALSE),
  fMaxZVert(999.),
  fDo3Prong(kFALSE),
  fMaxDecVertRadius2(8),
  fMassDzero(0.),
  fMassJpsi(0.),
  fMassDplus(0.),
  fMassDs(0.),
  fMassLambdaC(0.),
  fMassK0s(0.),
  fSecVertexerAlgo(0),
  fVertexerTracks{nullptr},
  fO2Vertexer2Prong{},
  fO2Vertexer3Prong{},
  fVertexerPropagateToPCA(true),
  fVertexerMaxR(200.),
  fVertexerMaxDZIni(4.),
  fVertexerMinParamChange(1.e-3),
  fVertexerMinRelChi2Change(0.9),
  fVertexerUseAbsDCA(true),
  fTrackCuts2pr{nullptr},
  fTrackCuts3pr{nullptr},
  fTrackCutsBach{nullptr},
  fTrackCutsV0Dau{nullptr},
  fMaxTracksToProcess(9999999),
  fNPtBinsSingleTrack(6),
  fNPtBinsDzero(25),
  fPtWithoutVtxToll(0.1),
  fMinPtDzero(0.),
  fMaxPtDzero(9999.),
  fMinPtDplus(0.),
  fMaxPtDplus(9999.),
  fMinPtJpsi(0.),
  fMaxPtJpsi(9999.),
  fCandidateCutLevel(1),
  fSelectD0(1),
  fSelectD0bar(1),
  fMinPt3Prong(0.),
  fMaxRapidityCand(-999.),
  fNPtBinsDplus(12),
  fNPtBinsJpsi(9),
  fNPtBinsLc(10),
  fSelectDplus(1),
  fSelectJpsi(1),
  fSelectLcpKpi(1),
  fFindVertexForCascades(kFALSE),
  fMinPtV0(0.),
  fMinCosPointV0(0.),
  fCutOnK0sMass(0.1),
  fEnableCPUTimeCheck(kFALSE),
  fCountTimeInMilliseconds(kFALSE)
{
  //
  InitDefault();

  for(Int_t i=0; i<5; i++){
    fHistPtGenPrompt[i]=0x0;
    fHistPtGenFeeddw[i]=0x0;
    fHistPtGenLimAccPrompt[i]=0x0;
    fHistPtGenLimAccFeeddw[i]=0x0;
    fHistPtRecoGenPtPrompt[i]=0x0;
    fHistPtRecoGenPtFeeddw[i]=0x0;
    fHistPtRecoPrompt[i]=0x0;
    fHistPtRecoFeeddw[i]=0x0;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}


//___________________________________________________________________________
AliAnalysisTaskHFSimpleVertices::~AliAnalysisTaskHFSimpleVertices(){
  //
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;

  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    delete fHistTrackStatus;
    delete fHistPtAllTracks;
    delete fHistPtSelTracks;
    delete fHistTglAllTracks;
    delete fHistTglSelTracks;
    delete fHistEtaAllTracks;
    delete fHistEtaSelTracks2prong;
    delete fHistEtaSelTracks3prong;
    delete fHistImpParAllTracks;
    delete fHistImpParSelTracks2prong;
    delete fHistImpParSelTracks3prong;
    delete fHistITSmapAllTracks;
    delete fHistITSmapSelTracks;
    delete fHistPrimVertX;
    delete fHistPrimVertY;
    delete fHistPrimVertZ;
    delete fHist2ProngVertX;
    delete fHist2ProngVertY;
    delete fHist2ProngVertZ;
    delete fHist3ProngVertX;
    delete fHist3ProngVertY;
    delete fHist3ProngVertZ;
    delete fHistDist12LcpKpi;
    delete fHistInvMassD0;
    delete fHistPtD0;
    delete fHistYPtD0;
    delete fHistPtD0Dau0;
    delete fHistPtD0Dau1;
    delete fHistImpParD0Dau0;
    delete fHistImpParD0Dau1;
    delete fHistd0Timesd0;
    delete fHistCosPointD0;
    delete fHistDecLenD0;
    delete fHistDecLenXYD0;
    delete fHistImpParErrD0Dau;
    delete fHistDecLenErrD0;
    delete fHistDecLenXYErrD0;
    delete fHistCovMatPrimVXX2Prong;
    delete fHistCovMatSecVXX2Prong;
    delete fHistD0SignalVertX;
    delete fHistD0SignalVertY;
    delete fHistD0SignalVertZ;
    delete fHistInvMassD0Signal;
    delete fHistInvMassD0Refl;
    delete fHistInvMassJpsi;
    delete fHistPtJpsi;
    delete fHistPtJpsiDau0;
    delete fHistPtJpsiDau1;
    delete fHistImpParJpsiDau0;
    delete fHistImpParJpsiDau1;
    delete fHistd0Timesd0Jpsi;
    delete fHistCosPointJpsi;
    delete fHistDecLenJpsi;
    delete fHistDecLenXYJpsi;
    delete fHistDecLenErrJpsi;
    delete fHistDecLenXYErrJpsi;
    delete fHistJpsiSignalVertX;
    delete fHistJpsiSignalVertY;
    delete fHistJpsiSignalVertZ;
    delete fHistInvMassJpsiSignal;
    delete fHistInvMassDplus;
    delete fHistInvMassDplusSignal;
    delete fHistPtDplus;
    delete fHistYPtDplus;
    delete fHistPtDplusDau0;
    delete fHistPtDplusDau1;
    delete fHistPtDplusDau2;
    delete fHistImpParDplusDau0;
    delete fHistImpParDplusDau1;
    delete fHistImpParDplusDau2;
    delete fHistDecLenDplus;
    delete fHistDecLenXYDplus;
    delete fHistNormDecLenXYDplus;
    delete fHistImpParErrDplusDau;
    delete fHistDecLenErrDplus;
    delete fHistDecLenXYErrDplus;
    delete fHistCosPointDplus;
    delete fHistCosPointXYDplus;
    delete fHistImpParXYDplus;
    delete fHistNormIPDplus;
    delete fHistoSumSqImpParDplusDau;
    delete fHistCovMatPrimVXX3Prong;
    delete fHistCovMatSecVXX3Prong;
    delete fHistInvMassDs;
    delete fHistInvMassDsSignal;
    delete fHistInvMassDsRefl;
    delete fHistPtDs;
    delete fHistYPtDs;
    delete fHistDecLenDs;
    delete fHistCosPointDs;
    delete fHistInvMassLc;
    delete fHistPtLc;
    delete fHistYPtLc;
    delete fHistPtLcDau0;
    delete fHistPtLcDau1;
    delete fHistPtLcDau2;
    delete fHistDecLenLc;
    delete fHistCosPointLc;
    delete fHistInvMassK0s;
    delete fHistInvMassLcK0sp;
    delete fHistPtLcK0sp;
    for(Int_t i=0; i<5; i++){
      delete fHistPtGenPrompt[i];
      delete fHistPtGenFeeddw[i];
      delete fHistPtGenLimAccPrompt[i];
      delete fHistPtGenLimAccFeeddw[i];
      delete fHistPtRecoGenPtPrompt[i];
      delete fHistPtRecoGenPtFeeddw[i];
      delete fHistPtRecoPrompt[i];
      delete fHistPtRecoFeeddw[i];
    }
    if(fHistCPUTimeTrackVsNTracks) delete fHistCPUTimeTrackVsNTracks;
    if(fHistCPUTimeCandVsNTracks) delete fHistCPUTimeCandVsNTracks;
    if(fHistWallTimeTrackVsNTracks) delete fHistWallTimeTrackVsNTracks;
    if(fHistWallTimeCandVsNTracks) delete fHistWallTimeCandVsNTracks;
  }
  delete fOutput;
  delete fTrackCuts2pr;
  delete fTrackCuts3pr;
  delete fTrackCutsBach;
  delete fTrackCutsV0Dau;
  delete fVertexerTracks;
}

//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::InitDefault(){
  /// initialization with default values

  fMassDzero = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  fMassJpsi = TDatabasePDG::Instance()->GetParticle(443)->Mass();
  fMassDplus = TDatabasePDG::Instance()->GetParticle(411)->Mass();
  fMassDs = TDatabasePDG::Instance()->GetParticle(431)->Mass();
  fMassLambdaC = TDatabasePDG::Instance()->GetParticle(4122)->Mass();
  fMassK0s = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  fTrackCuts2pr = new AliESDtrackCuts("AliESDtrackCuts", "default");
  fTrackCuts2pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts2pr->SetMinNClustersTPC(50);
  fTrackCuts2pr->SetRequireITSRefit(kTRUE);
  fTrackCuts2pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
  // fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // fTrackCuts->SetMaxDCAToVertexZ(3.2);
  // fTrackCuts->SetMaxDCAToVertexXY(2.4);
  // fTrackCuts->SetDCAToVertex2D(kTRUE);

  // no possibility to apply pT-dependent DCA cut in AliESDtrackCuts --> done by hand
  fTrackCuts2pr->SetMaxDCAToVertexXY(1000.);
  fTrackCuts2pr->SetMinDCAToVertexXY(0.);

  fTrackCuts3pr = new AliESDtrackCuts("AliESDtrackCuts", "default3p");
  fTrackCuts3pr->SetPtRange(0., 1.e10);
  // fTrackCuts->SetEtaRange(-0.8, +0.8);
  fTrackCuts3pr->SetMinNClustersTPC(50);
  fTrackCuts3pr->SetRequireITSRefit(kTRUE);
  fTrackCuts3pr->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);

  // no possibility to apply pT-dependent DCA cut in AliESDtrackCuts --> done by hand
  fTrackCuts3pr->SetMaxDCAToVertexXY(1000.);
  fTrackCuts3pr->SetMinDCAToVertexXY(0.);

  fTrackCutsBach = new AliESDtrackCuts("AliESDtrackCuts", "defaultbach");
  fTrackCutsBach->SetPtRange(0.3, 1.e10);
  fTrackCutsBach->SetEtaRange(-0.8, +0.8);
  fTrackCutsBach->SetMinDCAToVertexXYPtDep("0.*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  fTrackCutsBach->SetMaxDCAToVertexXY(1.);
  
  fTrackCutsV0Dau = new AliESDtrackCuts("AliESDtrackCuts", "defaultV0dau");
  fTrackCutsV0Dau->SetPtRange(0.05, 1.e10);
  fTrackCutsV0Dau->SetEtaRange(-1.1, +1.1);
  fTrackCutsV0Dau->SetMinNCrossedRowsTPC(50);
  fTrackCutsV0Dau->SetRequireTPCRefit(kTRUE);
  
  fNPtBinsSingleTrack = 6;
  Double_t defaultPtBinsSingleTrack[7] = {0., 0.5, 1., 1.5, 2., 3., 1000.};
  for (Int_t ib = 0; ib < fNPtBinsSingleTrack + 1; ib++)
      fPtBinLimsSingleTrack[ib] = defaultPtBinsSingleTrack[ib];

  Double_t defaultSingleTrackCuts[6][kNCutVarsSingleTrack] =
    {{0.0025, 1000.},  /* 0   < pt < 0.5    */
     {0.0025, 1000.},  /* 0.5 < pt < 1      */
     {0.0025, 1000.},  /* 1   < pt < 1.5    */
     {0.0025, 1000.},  /* 1.5 < pt < 2      */
     {0.0000, 1000.},  /* 2   < pt < 3      */
     {0.0000, 1000.}}; /* 3   < pt < 1000   */

  for (Int_t ib = 0; ib < fNPtBinsSingleTrack; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsSingleTrack; jc++) {
      fSingleTrackCuts2Prong[ib][jc] = defaultSingleTrackCuts[ib][jc];
      fSingleTrackCuts3Prong[ib][jc] = defaultSingleTrackCuts[ib][jc];
    }
  }

  Double_t defaultPtBins2Prongs[3] = {0., 5., 1000.};
  Double_t defaultPtBins3Prongs[3] = {0., 5., 1000.};
  for (Int_t ib = 0; ib < 2 + 1; ib++) {
    fPtBinLimsDzeroSkims[ib] = defaultPtBins2Prongs[ib];
    fPtBinLimsJpsiSkims[ib] = defaultPtBins2Prongs[ib];
    fPtBinLimsDplusSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsDsSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsLcSkims[ib] = defaultPtBins3Prongs[ib];
    fPtBinLimsXicSkims[ib] = defaultPtBins3Prongs[ib];
  }

  Double_t defaultDzeroSkimCuts[2][kNCutVars2Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultJpsiSkimCuts[2][kNCutVars2Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultDplusSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultDsSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultLcSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  Double_t defaultXicSkimCuts[2][kNCutVars3Prong] =
    {{0., 3., -1.1, 999999.},  /* 0   < pt < 5    */
     {0., 3., -1.1, 999999.}}; /* 5   < pt < 1000   */

  for(Int_t ib=0; ib<fNPtBins; ib++){
    for(Int_t jc=0; jc<kNCutVars2Prong; jc++){
      fDzeroSkimCuts[ib][jc]=defaultDzeroSkimCuts[ib][jc];
      fJpsiSkimCuts[ib][jc]=defaultJpsiSkimCuts[ib][jc];
      fDplusSkimCuts[ib][jc]=defaultDplusSkimCuts[ib][jc];
      fDsSkimCuts[ib][jc]=defaultDsSkimCuts[ib][jc];
      fLcSkimCuts[ib][jc]=defaultLcSkimCuts[ib][jc];
      fXicSkimCuts[ib][jc]=defaultXicSkimCuts[ib][jc];
    }
  }

  fNPtBinsDzero=25;
  Double_t defaultPtBins[26] = {0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0, 16.0, 20.0, 24.0, 36.0, 50.0, 100.0};
  for(Int_t ib=0; ib<fNPtBinsDzero+1; ib++) fPtBinLimsDzero[ib]=defaultPtBins[ib];

  Double_t defaultD0Cuts[25][kNCutVarsDzero] =
    {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* pt<0.5*/
     {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   /* 0.5<pt<1*/
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1<pt<1.5 */
     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  /* 1.5<pt<2 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2<pt<2.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  /* 2.5<pt<3 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3<pt<3.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  /* 3.5<pt<4 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4<pt<4.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 4.5<pt<5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5<pt<5.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 5.5<pt<6 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6<pt<6.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   /* 6.5<pt<7 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7<pt<7.5 */
     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   /* 7.5<pt<8 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 8<pt<9 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 9<pt<10 */
     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   /* 10<pt<12 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0.},   /* 12<pt<16 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 16<pt<20 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 20<pt<24 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 24<pt<36 */
     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  /* 36<pt<50 */
     {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0.}}; /* pt>50 */
  for(Int_t ib=0; ib<fNPtBinsDzero; ib++){
   for(Int_t jc=0; jc<kNCutVarsDzero; jc++){
     fDzeroCuts[ib][jc]=defaultD0Cuts[ib][jc];
   }
  }

  fNPtBinsDplus = 12;
  Double_t defaultPtBinsDplus[13] = {1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36.};
  for (Int_t ib = 0; ib < fNPtBinsDplus + 1; ib++)
    fPtBinLimsDplus[ib] = defaultPtBinsDplus[ib];

  Double_t defaultDplusCuts[12][kNCutVarsDplus] =
    {{0.2, 0.3, 0.3, 0.07, 6., 0.96, 0.985, 2.5},  /* 1<pt<2   */
     {0.2, 0.3, 0.3, 0.07, 5., 0.96, 0.985, 2.5},  /* 2<pt<3   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.980, 2.5},  /* 3<pt<4   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 4<pt<5   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 5<pt<6   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 6<pt<7   */
     {0.2, 0.3, 0.3, 0.10, 5., 0.96, 0.000, 2.5},  /* 7<pt<8   */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 8<pt<10  */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 10<pt<12 */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 12<pt<16 */
     {0.2, 0.3, 0.3, 0.12, 5., 0.96, 0.000, 2.5},  /* 16<pt<24 */
     {0.2, 0.3, 0.3, 0.20, 5., 0.94, 0.000, 2.5}}; /* 24<pt<36 */
  for(Int_t ib=0; ib<fNPtBinsDplus; ib++){
   for(Int_t jc=0; jc<kNCutVarsDplus; jc++){
     fDplusCuts[ib][jc]=defaultDplusCuts[ib][jc];
   }
  }

  fNPtBinsLc = 10;
  Double_t defaultPtBinsLc[11] = {0., 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 24.0, 36.0};
  for (Int_t ib = 0; ib < fNPtBinsLc + 1; ib++)
    fPtBinLimsLc[ib] = defaultPtBinsLc[ib];

  Double_t defaultLcCuts[10][kNCutVarsLc] =
    {{0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* pt<1*/
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 1<pt<2*/
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 2<pt<3 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 3<pt<4 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 4<pt<5 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 5<pt<6 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 6<pt<8 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 8<pt<12 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.},  /* 12<pt<24 */
     {0.400, 0.4, 0.4, 0.4, 0.01, 0.09, 0.005, 0.}}; /* 24<pt<36 */
  for (Int_t ib = 0; ib < fNPtBinsLc; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsLc; jc++) {
      fLcCuts[ib][jc] = defaultLcCuts[ib][jc];
    }
  }

  fNPtBinsJpsi = 9;
  Double_t defaultPtBinsJpsi[10] = {0, 0.5, 1., 2., 3., 4., 5., 7., 10., 15.};
  for (Int_t ib = 0; ib < fNPtBinsJpsi + 1; ib++)
      fPtBinLimsJpsi[ib] = defaultPtBinsJpsi[ib];

  Double_t defaultJpsiCuts[9][kNCutVarsJpsi] =
    {{0.5, 0.2, 0.4, 1},  /* pt<0.5    */
    {0.5, 0.2, 0.4,  1},  /* 0.5<pt<1   */
    {0.5, 0.2, 0.4,  1},  /* 1<pt<2   */
    {0.5, 0.2, 0.4,  1},  /* 2<pt<3   */
    {0.5, 0.2, 0.4,  1},  /* 3<pt<4   */
    {0.5, 0.2, 0.4,  1},  /* 4<pt<5   */
    {0.5, 0.2, 0.4,  1},  /* 5<pt<7   */
    {0.5, 0.2, 0.4,  1},  /* 7<pt<10  */
    {0.5, 0.2, 0.4,  1}}; /* 10<pt<15 */

  for (Int_t ib = 0; ib < fNPtBinsJpsi; ib++) {
    for (Int_t jc = 0; jc < kNCutVarsJpsi; jc++) {
      fJpsiCuts[ib][jc] = defaultJpsiCuts[ib][jc];
    }
  }
}

//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::InitFromJson(TString filename){
  /// read configuration from json file
  if (filename != "" && gSystem->Exec(Form("ls %s > /dev/null", filename.Data())) == 0) {
    printf("------Read configuration from JSON file------\n");

    std::string triggerMaskFromJSON = GetJsonString(filename.Data(), "triggerClassName");

    if (triggerMask.find(triggerMaskFromJSON) != triggerMask.end()) {
      fUsePhysSel = kTRUE;
      fTriggerMask = triggerMask[triggerMaskFromJSON];
    }

    Double_t ptmintrack2 = GetJsonFloat(filename.Data(), "pTMinTrack2Prong");
    printf("Min pt track (2 prong)= %f\n", ptmintrack2);
    if(ptmintrack2>0) fTrackCuts2pr->SetPtRange(ptmintrack2, 1.e10);
    Double_t ptmintrack3 = GetJsonFloat(filename.Data(), "pTMinTrack3Prong");
    printf("Min pt track (3 prong)= %f\n", ptmintrack3);
    if(ptmintrack3>0) fTrackCuts3pr->SetPtRange(ptmintrack3, 1.e10);
    Int_t do3Prongs = GetJsonInteger(filename.Data(), "do3prong");
    printf("do3prong     = %d\n", do3Prongs);
    if(do3Prongs>0) fDo3Prong=kTRUE;
    Int_t selectD0 = GetJsonInteger(filename.Data(), "d_selectionFlagD0");
    printf("d_selectionFlagD0 = %d\n",selectD0);
    if(selectD0>=0) fSelectD0=selectD0;
    Int_t selectD0bar = GetJsonInteger(filename.Data(), "d_selectionFlagD0bar");
    printf("d_selectionFlagD0bar = %d\n",selectD0bar);
    if(selectD0>=0) fSelectD0bar=selectD0bar;
    Int_t selectDplus = GetJsonInteger(filename.Data(), "d_selectionFlagDPlus");
    printf("d_selectionFlagDplus = %d\n",selectDplus);
    if(selectDplus>=0) fSelectDplus=selectDplus;
    Int_t selectJpsi = GetJsonInteger(filename.Data(), "d_selectionFlagJpsi");
    printf("d_selectionFlagJpsi = %d\n",selectJpsi);
    if(selectJpsi>=0) fSelectJpsi=selectJpsi;
    Int_t selectLcpKpi = GetJsonInteger(filename.Data(), "d_selectionFlagLc");
    printf("d_selectionFlagLc = %d\n", selectLcpKpi);
    if (selectLcpKpi >= 0) fSelectLcpKpi=selectLcpKpi;
    Int_t minncluTPC = GetJsonInteger(filename.Data(), "tpcNClsFound");
    if(minncluTPC>0) printf("minncluTPC   = %d\n", minncluTPC);
    fTrackCuts2pr->SetMinNClustersTPC(minncluTPC);
    fTrackCuts3pr->SetMinNClustersTPC(minncluTPC);

    int nptbinlimsSingleTrack = 0;
    float* ptbinsSingleTrack = GetJsonArray(filename.Data(),"pTBinsTrack",nptbinlimsSingleTrack);
    int npt2Prong = 0, npt3Prong = 0, nc2Prong = 0, nc3Prong = 0;
    float** cutsSingleTrack2Prong = GetJsonMatrix(filename.Data(),"cutsTrack2Prong",npt2Prong,nc2Prong);
    float** cutsSingleTrack3Prong = GetJsonMatrix(filename.Data(),"cutsTrack3Prong",npt3Prong,nc3Prong);
    if((nptbinlimsSingleTrack-1 != npt3Prong) || (nptbinlimsSingleTrack-1 != npt2Prong))
      AliFatal("Number of pT bins in JSON for single track cuts of 2-prong and 3-prongs not consistent, please check it");

    for (Int_t ib = 0; ib < nptbinlimsSingleTrack; ib++) {
      fPtBinLimsSingleTrack[ib] = ptbinsSingleTrack[ib];
    }
    for (Int_t ib = 0; ib < nptbinlimsSingleTrack-1; ib++) {
      for (Int_t jc = 0; jc < nc2Prong; jc++) {
        fSingleTrackCuts2Prong[ib][jc] = cutsSingleTrack2Prong[ib][jc];
        AliWarning(Form("2prong %d, %d, %f", ib, jc, cutsSingleTrack2Prong[ib][jc]));
      }
      for (Int_t jc = 0; jc < nc3Prong; jc++) {
        fSingleTrackCuts3Prong[ib][jc] = cutsSingleTrack3Prong[ib][jc];
        AliWarning(Form("3prong %d %d, %d, %f", nc2Prong, ib, jc, cutsSingleTrack2Prong[ib][jc]));
      }
    }
    Double_t etamax2 = GetJsonFloat(filename.Data(), "etaMax2Prong");
    printf("Max eta  (2 prong) = %f\n", etamax2);
    if(etamax2>0) fTrackCuts2pr->SetEtaRange(-etamax2, +etamax2);
    Double_t etamax3 = GetJsonFloat(filename.Data(), "etaMax3Prong");
    printf("Max eta  (3 prong) = %f\n", etamax3);
    if(etamax3>0) fTrackCuts3pr->SetEtaRange(-etamax3, +etamax3);

    // vertexer parameters
    printf("--- DCAFitterN parameters ---\n");
    Int_t b_propdca = GetJsonBool(filename.Data(),"propToDCA");
    if(b_propdca==1){
      fVertexerPropagateToPCA = true;
      printf("propdca = %d\n",fVertexerPropagateToPCA);
    }else if(b_propdca==0){
      fVertexerPropagateToPCA = false;
      printf("propdca = %d\n",fVertexerPropagateToPCA);
    }
    Double_t d_maxr = GetJsonFloat(filename.Data(), "maxRad");
    if(d_maxr>0){
      fMaxDecVertRadius2=d_maxr*d_maxr;
      fVertexerMaxR=d_maxr;
      printf("maxr = %f\n",fVertexerMaxR);
    }
    Double_t d_maxdzini = GetJsonFloat(filename.Data(), "maxDZIni");
    if(d_maxdzini>0){
      fVertexerMaxDZIni=d_maxdzini;
      printf("maxdzini = %f\n",fVertexerMaxDZIni);
    }
    Double_t d_minparamchange = GetJsonFloat(filename.Data(), "minParamChange");
    if(d_minparamchange>0){
      fVertexerMinParamChange=d_minparamchange;
      printf("minparamchange = %f\n",fVertexerMinParamChange);
    }
    Double_t d_minrelchi2change = GetJsonFloat(filename.Data(), "minRelChi2Change");
    if(d_minrelchi2change){
      fVertexerMinRelChi2Change=d_minrelchi2change;
      printf("minrelchi2change = %f\n",fVertexerMinRelChi2Change);
    }
    printf("----------------\n");
    Double_t ptMinCand = GetJsonFloat(filename.Data(), "d_pTCandMin");
    printf("Min pt Dzero cand = %f\n", ptMinCand);
    if(ptMinCand>=0.) fMinPtDzero=ptMinCand;
    Double_t ptMaxCand = GetJsonFloat(filename.Data(), "d_pTCandMax");
    printf("Max pt Dzero cand = %f\n", ptMaxCand);
    if(ptMaxCand>=0. && ptMaxCand>=fMinPtDzero) fMaxPtDzero=ptMaxCand;
    Double_t ptMinCand3 = GetJsonFloat(filename.Data(), "cut3ProngPtCandMin");
    printf("Min pt 3-prong cand = %f\n", ptMinCand3);
    if( ptMinCand3>=0.) fMinPt3Prong=ptMinCand3;

    // Selections used in the skimming
    printf("------- CANDIDATE SELECTIONS FOR SKIMMING -------\n");

    Double_t ptTol = GetJsonFloat(filename.Data(), "pTTolerance");
    if(ptTol > 0)
      fPtWithoutVtxToll = ptTol;

    int nptbinlimsDzeroSkims = 0, ncDzeroSkims = 0, nptDzeroSkims = 0;
    float* ptbinlimsDzeroSkims = GetJsonArray(filename.Data(),"pTBinsD0ToPiK",nptbinlimsDzeroSkims);
    float** cutsDzeroSkims = GetJsonMatrix(filename.Data(),"cutsD0ToPiK",nptDzeroSkims,ncDzeroSkims);
    if(nptbinlimsDzeroSkims-1 != nptDzeroSkims)
      AliFatal("Number of pT bins in JSON for Dzero at skims level not consistent, please check it");

    int nptbinlimsJpsiSkims = 0, ncJpsiSkims = 0, nptJpsiSkims = 0;
    float* ptbinlimsJpsiSkims = GetJsonArray(filename.Data(),"pTBinsJpsiToEE",nptbinlimsJpsiSkims);
    float** cutsJpsiSkims = GetJsonMatrix(filename.Data(),"cutsJpsiToEE",nptJpsiSkims,ncJpsiSkims);
    if(nptbinlimsJpsiSkims-1 != nptJpsiSkims)
      AliFatal("Number of pT bins in JSON for J/psi at skims level not consistent, please check it");

    int nptbinlimsDplusSkims = 0, ncDplusSkims = 0, nptDplusSkims = 0;
    float* ptbinlimsDplusSkims = GetJsonArray(filename.Data(),"pTBinsDPlusToPiKPi",nptbinlimsDplusSkims);
    float** cutsDplusSkims = GetJsonMatrix(filename.Data(),"cutsDPlusToPiKPi",nptDplusSkims,ncDplusSkims);
    if(nptbinlimsDplusSkims-1 != nptDplusSkims)
      AliFatal("Number of pT bins in JSON for Dplus at skims level not consistent, please check it");

    int nptbinlimsDsSkims = 0, ncDsSkims = 0, nptDsSkims = 0;
    float* ptbinlimsDsSkims = GetJsonArray(filename.Data(),"pTBinsDsToPiKK",nptbinlimsDsSkims);
    float** cutsDsSkims = GetJsonMatrix(filename.Data(),"cutsDsToPiKK",nptDsSkims,ncDsSkims);
    if(nptbinlimsDsSkims-1 != nptDsSkims)
      AliFatal("Number of pT bins in JSON for Ds at skims level not consistent, please check it");

    int nptbinlimsLcSkims = 0, ncLcSkims = 0, nptLcSkims = 0;
    float* ptbinlimsLcSkims = GetJsonArray(filename.Data(),"pTBinsLcToPKPi",nptbinlimsLcSkims);
    float** cutsLcSkims = GetJsonMatrix(filename.Data(),"cutsLcToPKPi",nptLcSkims,ncLcSkims);
    if(nptbinlimsLcSkims-1 != nptLcSkims)
      AliFatal("Number of pT bins in JSON for Lc at skims level not consistent, please check it");

    int nptbinlimsXicSkims = 0, ncXicSkims = 0, nptXicSkims = 0;
    float* ptbinlimsXicSkims = GetJsonArray(filename.Data(),"pTBinsXicToPKPi",nptbinlimsXicSkims);
    float** cutsXicSkims = GetJsonMatrix(filename.Data(),"cutsXicToPKPi",nptXicSkims,ncXicSkims);
    if(nptbinlimsXicSkims-1 != nptXicSkims)
      AliFatal("Number of pT bins in JSON for Xic at skims level not consistent, please check it");

    for (Int_t ib = 0; ib < nptbinlimsDzeroSkims; ib++) {
      fPtBinLimsDzeroSkims[ib] = ptbinlimsDzeroSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDzeroSkims; ib++) {
      for (Int_t jc = 0; jc < ncDzeroSkims; jc++) {
        fDzeroSkimCuts[ib][jc] = cutsDzeroSkims[ib][jc];
      }
      AliWarning(Form("Dzero cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; d0xd0  < %f\n", fPtBinLimsDzeroSkims[ib], fPtBinLimsDzeroSkims[ib+1], fDzeroSkimCuts[ib][0], fDzeroSkimCuts[ib][1], fDzeroSkimCuts[ib][2], fDzeroSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsJpsiSkims; ib++) {
      fPtBinLimsJpsiSkims[ib] = ptbinlimsJpsiSkims[ib];
    }
    for (Int_t ib = 0; ib < nptJpsiSkims; ib++) {
      for (Int_t jc = 0; jc < ncJpsiSkims; jc++) {
        fJpsiSkimCuts[ib][jc] = cutsJpsiSkims[ib][jc];
      }
      AliWarning(Form("J/psi cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; d0xd0  < %f\n", fPtBinLimsJpsiSkims[ib], fPtBinLimsJpsiSkims[ib+1], fJpsiSkimCuts[ib][0], fJpsiSkimCuts[ib][1], fJpsiSkimCuts[ib][2], fJpsiSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsDplusSkims; ib++) {
      fPtBinLimsDplusSkims[ib] = ptbinlimsDplusSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDplusSkims; ib++) {
      for (Int_t jc = 0; jc < ncDplusSkims; jc++) {
        fDplusSkimCuts[ib][jc] = cutsDplusSkims[ib][jc];
      }
      AliWarning(Form("Dplus cuts: %f < pt < %f  ; %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsDplusSkims[ib], fPtBinLimsDplusSkims[ib+1], fDplusSkimCuts[ib][0], fDplusSkimCuts[ib][1], fDplusSkimCuts[ib][2], fDplusSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsDsSkims; ib++) {
      fPtBinLimsDsSkims[ib] = ptbinlimsDsSkims[ib];
    }
    for (Int_t ib = 0; ib < nptDsSkims; ib++) {
      for (Int_t jc = 0; jc < ncDsSkims; jc++) {
        fDsSkimCuts[ib][jc] = cutsDsSkims[ib][jc];
      }
      AliWarning(Form("Ds cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsDsSkims[ib], fPtBinLimsDsSkims[ib+1], fDsSkimCuts[ib][0], fDsSkimCuts[ib][1], fDsSkimCuts[ib][2], fDsSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsLcSkims; ib++) {
      fPtBinLimsLcSkims[ib] = ptbinlimsLcSkims[ib];
    }
    for (Int_t ib = 0; ib < nptLcSkims; ib++) {
      for (Int_t jc = 0; jc < ncLcSkims; jc++) {
        fLcSkimCuts[ib][jc] = cutsLcSkims[ib][jc];
      }
      AliWarning(Form("Lc cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsLcSkims[ib], fPtBinLimsLcSkims[ib+1], fLcSkimCuts[ib][0], fLcSkimCuts[ib][1], fLcSkimCuts[ib][2], fLcSkimCuts[ib][3]));
    }

    for (Int_t ib = 0; ib < nptbinlimsXicSkims; ib++) {
      fPtBinLimsXicSkims[ib] = ptbinlimsXicSkims[ib];
    }
    for (Int_t ib = 0; ib < nptXicSkims; ib++) {
      for (Int_t jc = 0; jc < ncXicSkims; jc++) {
        fXicSkimCuts[ib][jc] = cutsXicSkims[ib][jc];
      }
      AliWarning(Form("Xic cuts: %f < pt < %f  ;  %f < mass < %f  ;  cospoint > %f  ; declen > %f\n", fPtBinLimsXicSkims[ib], fPtBinLimsXicSkims[ib+1], fXicSkimCuts[ib][0], fXicSkimCuts[ib][1], fXicSkimCuts[ib][2], fXicSkimCuts[ib][3]));
    }

    
    Double_t cutcpaV0 = GetJsonFloat(filename.Data(), "cosPAV0");
    if(cutcpaV0>0){
      fMinCosPointV0=cutcpaV0;
      printf("cosPAV0 cut = %f\n",fMinCosPointV0);
    }

    Double_t cutinvmV0 = GetJsonFloat(filename.Data(), "cutInvMassV0");
    if(cutinvmV0>0){
      fCutOnK0sMass=cutinvmV0;
      printf("invmass V0 cut = %f\n",fCutOnK0sMass);
    }
    
    
    printf("---------------------------------------------\n");

    printf("---- TEST READOUT OF ARRAYS ----\n");
    int nptbinlims;
    float* ptbins = GetJsonArray(filename.Data(),"pTBins",nptbinlims);
    printf("%d ptbins. Limits: ",nptbinlims-1);
    for(int j=0; j<nptbinlims; j++) printf("%.1f ",ptbins[j]);
    printf("\n");
    int npt,nc;
    float** cuts = GetJsonMatrix(filename.Data(),"D0_to_pi_K_cuts",npt,nc);
    printf("D0 2D cuts: %d pt bins, %d cut variables:\n",npt,nc);
    for(int j=0; j<npt; j++){
      for(int k=0; k<nc; k++){
        printf("%.2f ",cuts[j][k]);
      }
      printf("\n");
    }
    printf("\n");
  }else{
    AliError(Form("Json configuration file %s not found\n",filename.Data()));
  }
}

//___________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::UserCreateOutputObjects() {
  // create output histos

  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");

  fHistNEvents = new TH1F("hNEvents", "Number of processed events",15,-0.5,14.5);
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"All events");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"PhysSel");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"InCentralityClass");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Good vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Pass zSPD-zTrk vert sel");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Pass |zvert|");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Pileup cut");
  fOutput->Add(fHistNEvents);

  fHistTrackStatus = new TH1F("hTrackStatus","",4,-0.5,3.5);
  fHistTrackStatus->GetXaxis()->SetBinLabel(1,"Rejected");
  fHistTrackStatus->GetXaxis()->SetBinLabel(2,"2-prong");
  fHistTrackStatus->GetXaxis()->SetBinLabel(3,"3-prong");
  fHistTrackStatus->GetXaxis()->SetBinLabel(4,"bachelor");
  fOutput->Add(fHistTrackStatus);
  
  // single track histos
  fHistPtAllTracks = new TH1F("hPtAllTracks", " All tracks ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtSelTracks = new TH1F("hPtSelTracks", " Selected tracks ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistTglAllTracks = new TH1F("hTglAllTracks", " All tracks ; tg#lambda", 100, -5., 5.);
  fHistTglSelTracks = new TH1F("hTglSelTracks", " Selected tracks ; tg#lambda", 100, -5., 5.);
  fHistEtaAllTracks = new TH1F("hEtaAllTracks", " All tracks ; #eta", 100, -1, 1.);
  Float_t eta2mincut,eta2maxcut;
  fTrackCuts2pr->GetEtaRange(eta2mincut,eta2maxcut);
  Int_t nBinsETA=100;
  Double_t lowlimETA=-1.;
  Double_t higlimETA=1.;
  if(eta2maxcut<10){
    nBinsETA=static_cast<int>(1.2 * eta2maxcut * 100);
    lowlimETA=1.2*eta2mincut;
    higlimETA=1.2*eta2maxcut;
  }
  fHistEtaSelTracks2prong = new TH1F("hEtaSelTracks2prong", " Selected tracks ; #eta", nBinsETA,lowlimETA,higlimETA);
  Float_t eta3mincut,eta3maxcut;
  fTrackCuts3pr->GetEtaRange(eta3mincut,eta3maxcut);
  nBinsETA=100;
  lowlimETA=-1.;
  higlimETA=1.;
  if(eta3maxcut<10){
    nBinsETA=static_cast<int>(1.2 * eta3maxcut * 100);
    lowlimETA=1.2*eta3mincut;
    higlimETA=1.2*eta3maxcut;
  }
  fHistEtaSelTracks3prong = new TH1F("hEtaSelTracks3prong", " Selected tracks ; #eta", nBinsETA,lowlimETA,higlimETA);
  fHistImpParAllTracks = new TH1F("hImpParAllTracks", " All tracks ; d_{0}^{xy} (cm)", 100, -1, 1.);
  Int_t nBinsDCA=400;
  Double_t lowlimDCA=-2.;
  Double_t higlimDCA=2.;
  fHistImpParSelTracks2prong = new TH1F("hImpParSelTracks2prong", " Selected tracks ; d_{0}^{xy} (cm)", nBinsDCA,lowlimDCA,higlimDCA);
  nBinsDCA=400;
  lowlimDCA=-2.;
  higlimDCA=2.;
  fHistImpParSelTracks3prong = new TH1F("hImpParSelTracks3prong", " Selected tracks ; d_{0}^{xy} (cm)", nBinsDCA,lowlimDCA,higlimDCA);
  fHistITSmapAllTracks = new TH1F("hITSmapAllTracks", " All tracks ; ITS cluster map", 64, -0.5, 63.5);
  fHistITSmapSelTracks = new TH1F("hITSmapSelTracks", " Selected tracks ; ITS cluster map", 64, -0.5, 63.5);
  fOutput->Add(fHistPtAllTracks);
  fOutput->Add(fHistPtSelTracks);
  fOutput->Add(fHistTglAllTracks);
  fOutput->Add(fHistTglSelTracks);
  fOutput->Add(fHistEtaAllTracks);
  fOutput->Add(fHistEtaSelTracks2prong);
  fOutput->Add(fHistEtaSelTracks3prong);
  fOutput->Add(fHistImpParAllTracks);
  fOutput->Add(fHistImpParSelTracks2prong);
  fOutput->Add(fHistImpParSelTracks3prong);
  fOutput->Add(fHistITSmapAllTracks);
  fOutput->Add(fHistITSmapSelTracks);

  // vertex histos
  fHistPrimVertX = new TH1F("hPrimVertX"," Primary Vertex ; x (cm)",100, -0.5, 0.5);
  fHistPrimVertY = new TH1F("hPrimVertY"," Primary Vertex ; y (cm)",100, -0.5, 0.5);
  fHistPrimVertZ = new TH1F("hPrimVertZ"," Primary Vertex ; z (cm)",100, -20.0, 20.0);
  fHist2ProngVertX = new TH1F("h2ProngVertX"," Secondary Vertex ; x (cm)",1000, -2., 2.);
  fHist2ProngVertY = new TH1F("h2ProngVertY"," Secondary Vertex ; y (cm)",1000, -2., 2.);
  fHist2ProngVertZ = new TH1F("h2ProngVertZ"," Secondary Vertex ; z (cm)",1000, -20.0, 20.0);
  fHist3ProngVertX = new TH1F("h3ProngVertX"," Secondary Vertex ; x (cm)",1000, -2.0, 2.0);
  fHist3ProngVertY = new TH1F("h3ProngVertY"," Secondary Vertex ; y (cm)",1000, -2.0, 2.0);
  fHist3ProngVertZ = new TH1F("h3ProngVertZ"," Secondary Vertex ; z (cm)",1000, -20.0, 20.0);
  fHistDist12LcpKpi = new TH1F("hDist12Lc", " ; Dist12 (cm)", 200, 0., 2.0);
  fOutput->Add(fHistPrimVertX);
  fOutput->Add(fHistPrimVertY);
  fOutput->Add(fHistPrimVertZ);
  fOutput->Add(fHist2ProngVertX);
  fOutput->Add(fHist2ProngVertY);
  fOutput->Add(fHist2ProngVertZ);
  fOutput->Add(fHist3ProngVertX);
  fOutput->Add(fHist3ProngVertY);
  fOutput->Add(fHist3ProngVertZ);
  fOutput->Add(fHistDist12LcpKpi);

  // D0 meson candidate histos
  fHistInvMassD0 = new TH1F("hInvMassD0", " ; M_{K#pi} (GeV/c^{2})", 500, 0, 5.0);
  fHistPtD0  = new TH1F("hPtD0", " ; D^{0} p_{T} (GeV/c)", 100, 0, 10.);
  fHistYPtD0  = new TH2F("hYPtD0", " ; D^{0} p_{T} (GeV/c) ; y", 100, 0, 10.,120,-1.2,1.2);
  fHistPtD0Dau0 = new TH1F("hPtD0Dau0", " D^{0} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtD0Dau1 = new TH1F("hPtD0Dau1", " D^{0} prong1 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistImpParD0Dau0 = new TH1F("hImpParD0Dau0", " D^{0} prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistImpParD0Dau1 = new TH1F("hImpParD0Dau1", " D^{0} prong1 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistd0Timesd0 = new TH1F("hd0Timesd0", " d_{0}^{xy}x d_{0}^{xy} (cm^{2})", 500, -1.0, 1.0);
  fHistCosPointD0 = new TH1F("hCosPointD0", " ; cos(#theta_{P})", 110, -1.1, 1.1);
  fHistDecLenD0 = new TH1F("hDecLenD0", " ; Decay Length (cm)", 200, 0., 2.0);
  fHistDecLenXYD0 = new TH1F("hDecLenXYD0", " ; Decay Length xy (cm)", 200, 0., 2.0);
  fHistImpParErrD0Dau = new TH1F("hImpParErrD0Dau", " D^{0} prongs ; #sigma(d_{0}^{xy}) (cm)", 100, -1.0, 1.0);
  fHistDecLenErrD0 = new TH1F("hDecLenErrD0", " ; #sigma(Decay Length) (cm)", 100, 0., 1.0);
  fHistDecLenXYErrD0 = new TH1F("hDecLenXYErrD0", " ; #sigma(Decay Length xy) (cm)", 100, 0., 1.0);
  fHistCovMatPrimVXX2Prong = new TH1F("hCovMatPrimVXX2Prong", " Primary Vertex ; XX element of covariant matrix (cm^2)", 100, 0., 1.0e-4);
  fHistCovMatSecVXX2Prong = new TH1F("hCovMatSecVXX2Prong", " Secondary Vertex 2-prong ; XX element of covariant matrix (cm^2)", 100, 0., 0.2);
  fOutput->Add(fHistPtD0);
  fOutput->Add(fHistYPtD0);
  fOutput->Add(fHistInvMassD0);
  fOutput->Add(fHistPtD0Dau0);
  fOutput->Add(fHistPtD0Dau1);
  fOutput->Add(fHistImpParD0Dau0);
  fOutput->Add(fHistImpParD0Dau1);
  fOutput->Add(fHistd0Timesd0);
  fOutput->Add(fHistCosPointD0);
  fOutput->Add(fHistDecLenD0);
  fOutput->Add(fHistDecLenXYD0);
  fOutput->Add(fHistImpParErrD0Dau);
  fOutput->Add(fHistDecLenErrD0);
  fOutput->Add(fHistDecLenXYErrD0);
  fOutput->Add(fHistCovMatPrimVXX2Prong);
  fOutput->Add(fHistCovMatSecVXX2Prong);

  // MC truth D0 histos
  fHistD0SignalVertX = new TH1F("hD0SignalVertX"," Secondary Vertex ; x (cm)",1000, -2., 2.);
  fHistD0SignalVertY = new TH1F("hD0SignalVertY"," Secondary Vertex ; y (cm)",1000, -2., 2.);
  fHistD0SignalVertZ = new TH1F("hD0SignalVertZ"," Secondary Vertex ; z (cm)",1000, -20.0, 20.0);
  fHistInvMassD0Signal = new TH1F("hInvMassD0Signal", " ; M_{K#pi} (GeV/c^{2})", 500, 0, 5.0);
  fHistInvMassD0Refl = new TH1F("hInvMassD0Refl", " ; M_{K#pi} (GeV/c^{2})", 500, 0, 5.0);
  fOutput->Add(fHistD0SignalVertX);
  fOutput->Add(fHistD0SignalVertY);
  fOutput->Add(fHistD0SignalVertZ);
  fOutput->Add(fHistInvMassD0Signal);
  fOutput->Add(fHistInvMassD0Refl);

  // Jpsi meson candidate histos
  fHistInvMassJpsi = new TH1F("hInvMassJpsi", " ; M_{e^{+}e_{-}} (GeV/c^{2})", 500, 0, 5.0);
  fHistPtJpsi  = new TH1F("hPtJpsi", " ; Jpsi p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtJpsiDau0 = new TH1F("hPtJpsiDau0", " Jpsi prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtJpsiDau1 = new TH1F("hPtJpsiDau1", " Jpsi prong1 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistImpParJpsiDau0 = new TH1F("hImpParJpsiDau0", " Jpsi prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistImpParJpsiDau1 = new TH1F("hImpParJpsiDau1", " Jpsi prong1 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistd0Timesd0Jpsi = new TH1F("hd0Timesd0Jpsi", " d_{0}^{xy}x d_{0}^{xy} (cm^{2})", 500, -1.0, 1.0);
  fHistCosPointJpsi = new TH1F("hCosPointJpsi", " ; cos(#theta_{P})", 110, -1.1, 1.1);
  fHistDecLenJpsi = new TH1F("hDecLenJpsi", " ; Decay Length (cm)", 200, 0., 2.0);
  fHistDecLenXYJpsi = new TH1F("hDecLenXYJpsi", " ; Decay Length xy (cm)", 200, 0., 2.0);
  fHistDecLenErrJpsi = new TH1F("hDecLenErrJpsi", " ; #sigma(Decay Length) (cm)", 100, 0., 1.0);
  fHistDecLenXYErrJpsi = new TH1F("hDecLenXYErrJpsi", " ; #sigma(Decay Length xy) (cm)", 100, 0., 1.0);

  fOutput->Add(fHistPtJpsi);
  fOutput->Add(fHistInvMassJpsi);
  fOutput->Add(fHistPtJpsiDau0);
  fOutput->Add(fHistPtJpsiDau1);
  fOutput->Add(fHistImpParJpsiDau0);
  fOutput->Add(fHistImpParJpsiDau1);
  fOutput->Add(fHistd0Timesd0Jpsi);
  fOutput->Add(fHistCosPointJpsi);
  fOutput->Add(fHistDecLenJpsi);
  fOutput->Add(fHistDecLenXYJpsi);
  fOutput->Add(fHistDecLenErrJpsi);
  fOutput->Add(fHistDecLenXYErrJpsi);

  fHistJpsiSignalVertX = new TH1F("hJpsiSignalVertX"," Secondary Vertex ; x (cm)",1000, -2., 2.);
  fHistJpsiSignalVertY = new TH1F("hJpsiSignalVertY"," Secondary Vertex ; y (cm)",1000, -2., 2.);
  fHistJpsiSignalVertZ = new TH1F("hJpsiSignalVertZ"," Secondary Vertex ; z (cm)",1000, -20.0, 20.0);
  fHistInvMassJpsiSignal = new TH1F("hInvMassJpsiSignal", " ; M_{K#pi} (GeV/c^{2})", 500, 0, 5.0);
  //fHistInvMassD0Refl = new TH1F("hInvMassD0Refl", " ; M_{K#pi} (GeV/c^{2})", 500, 0, 5.0);
  fOutput->Add(fHistJpsiSignalVertX);
  fOutput->Add(fHistJpsiSignalVertY);
  fOutput->Add(fHistJpsiSignalVertZ);
  fOutput->Add(fHistInvMassJpsiSignal);

  // Dplus meson candidate histos
  fHistInvMassDplus = new TH1F("hInvMassDplus", " ; M_{K#pi#pi} (GeV/c^{2})", 350, 1.7, 2.05);
  fHistInvMassDplusSignal = new TH1F("hInvMassDplusSignal", " ; M_{K#pi#pi} (GeV/c^{2})", 350, 1.7, 2.05);
  fHistPtDplus = new TH1F("hPtDplus", " ; D^{+} p_{T} (GeV/c)", 100, 0, 10.);
  fHistYPtDplus  = new TH2F("hYPtDplus", " ; D^{+} p_{T} (GeV/c) ; y", 100, 0, 10.,120,-1.2,1.2);
  fHistPtDplusDau0 = new TH1F("hPtDplusDau0", " D^{+} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtDplusDau1 = new TH1F("hPtDplusDau1", " D^{+} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtDplusDau2 = new TH1F("hPtDplusDau2", " D^{+} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistImpParDplusDau0 = new TH1F("hImpParDplusDau0", " D^{+} prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistImpParDplusDau1 = new TH1F("hImpParDplusDau1", " D^{+} prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistImpParDplusDau2 = new TH1F("hImpParDplusDau2", " D^{+} prong0 ; d_{0}^{xy} (cm)", 100, -1.0, 1.0);
  fHistDecLenDplus = new TH1F("hDecLenDplus", " ; Decay Length (cm)", 200, 0., 2.0);
  fHistDecLenXYDplus = new TH1F("hDecLenXYDplus", " ; Decay Length xy (cm)", 200, 0., 2.0);
  fHistNormDecLenXYDplus = new TH1F("hNormDecLenXYDplus", " ; Norm. Decay Length xy (cm)", 80, 0., 80.);
  fHistImpParErrDplusDau = new TH1F("hImpParErrDplusDau", " D^{+} prongs ; #sigma(d_{0}^{xy}) (cm)", 100, 0., 1.0);
  fHistDecLenErrDplus = new TH1F("hDecLenErrDplus", " ; #sigma(Decay Length) (cm)", 100, 0., 1.0);
  fHistDecLenXYErrDplus = new TH1F("hDecLenXYErrDplus", " ; #sigma(Decay Length xy) (cm)", 100, 0., 1.0);
  fHistCosPointDplus = new TH1F("hCosPointDplus", " ; cos(#theta_{P})", 110, -1.1, 1.1);
  fHistCosPointXYDplus = new TH1F("hCosPointXYDplus", " ; cos(#theta^{xy}_{P})", 110, -1.1, 1.1);
  fHistImpParXYDplus = new TH1F("hImpParXYDplus", " ; d_{0}^{xy} (cm)", 200, -1., 1.);
  fHistNormIPDplus = new TH1F("hNormIPDplus", " ; Norm. IP", 200, -20., 20.);
  fHistoSumSqImpParDplusDau = new TH1F("hSumSqImpParDplusDau", " ; Squared sum of d_{0}^{i} (cm^2)", 100, 0., 1.);
  fHistCovMatPrimVXX3Prong = new TH1F("hCovMatPrimVXX3Prong", " Primary Vertex ; XX element of covariant matrix (cm^2)", 100, 0., 1.0e-4);
  fHistCovMatSecVXX3Prong = new TH1F("hCovMatSecVXX3Prong", " Secondary Vertex 3-prong ; XX element of covariant matrix (cm^2)", 100, 0., 0.2);
  fOutput->Add(fHistInvMassDplus);
  fOutput->Add(fHistInvMassDplusSignal);
  fOutput->Add(fHistPtDplus);
  fOutput->Add(fHistYPtDplus);
  fOutput->Add(fHistPtDplusDau0);
  fOutput->Add(fHistPtDplusDau1);
  fOutput->Add(fHistPtDplusDau2);
  fOutput->Add(fHistImpParDplusDau0);
  fOutput->Add(fHistImpParDplusDau1);
  fOutput->Add(fHistImpParDplusDau2);
  fOutput->Add(fHistDecLenDplus);
  fOutput->Add(fHistDecLenXYDplus);
  fOutput->Add(fHistNormDecLenXYDplus);
  fOutput->Add(fHistImpParErrDplusDau);
  fOutput->Add(fHistDecLenErrDplus);
  fOutput->Add(fHistDecLenXYErrDplus);
  fOutput->Add(fHistCosPointDplus);
  fOutput->Add(fHistCosPointXYDplus);
  fOutput->Add(fHistImpParXYDplus);
  fOutput->Add(fHistNormIPDplus);
  fOutput->Add(fHistoSumSqImpParDplusDau);
  fOutput->Add(fHistCovMatPrimVXX3Prong);
  fOutput->Add(fHistCovMatSecVXX3Prong);

  // Ds->KKpi candidate histos
  fHistInvMassDs = new TH1F("hInvMassDs", " ; M_{KK#pi} (GeV/c^{2})", 500, 1.7, 2.2);
  fHistInvMassDsSignal = new TH1F("hInvMassDsSignal", " ; M_{KK#pi} (GeV/c^{2})", 500, 1.7, 2.2);
  fHistInvMassDsRefl = new TH1F("hInvMassDsRefl", " ; M_{KK#pi} (GeV/c^{2})", 500, 1.7, 2.2);
  fHistPtDs = new TH1F("hPtDs"," ; D_{s} p_{T} (GeV/c)",100, 0, 10.);
  fHistYPtDs  = new TH2F("hYPtDs", " ; D_{s} p_{T} (GeV/c) ; y", 100, 0, 10.,120,-1.2,1.2);
  fHistDecLenDs = new TH1F("hDecLenDs"," ; Decay Length (cm)", 200, 0., 2.0);
  fHistCosPointDs = new TH1F("hCosPointDs", " ; cos(#theta_{P})", 110, -1.1, 1.1);
  fOutput->Add(fHistInvMassDs);
  fOutput->Add(fHistInvMassDsSignal);
  fOutput->Add(fHistInvMassDsRefl);
  fOutput->Add(fHistPtDs);
  fOutput->Add(fHistYPtDs);
  fOutput->Add(fHistDecLenDs);
  fOutput->Add(fHistCosPointDs);

  // Lc pKpi candidate histos
  fHistInvMassLc = new TH1F("hInvMassLc", " ; M_{pK#pi} (GeV/c^{2})", 500, 1.6, 3.1);
  fHistPtLc = new TH1F("hPtLc", " ; #Lambda_{c} p_{T} (GeV/c)", 100, 0, 10.);
  fHistYPtLc  = new TH2F("hYPtLc", " ; #Lambda_{c} p_{T} (GeV/c) ; y", 100, 0, 10.,120,-1.2,1.2);
  fHistPtLcDau0 = new TH1F("hPtLcDau0", " #Lambda_{c} prong0 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtLcDau1 = new TH1F("hPtLcDau1", " #Lambda_{c} prong1 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistPtLcDau2 = new TH1F("hPtLcDau2", " #Lambda_{c} prong2 ; p_{T} (GeV/c)", 100, 0, 10.);
  fHistDecLenLc = new TH1F("hDecLenLc", " ; Decay Length (cm)", 200, 0., 2.0);
  fHistCosPointLc = new TH1F("hCosPointLc", " ; cos(#theta_{P})", 110, -1.1, 1.1);
  fOutput->Add(fHistInvMassLc);
  fOutput->Add(fHistPtLc);
  fOutput->Add(fHistYPtLc);
  fOutput->Add(fHistPtLcDau0);
  fOutput->Add(fHistPtLcDau1);
  fOutput->Add(fHistPtLcDau2);
  fOutput->Add(fHistDecLenLc);
  fOutput->Add(fHistCosPointLc);

  fHistInvMassK0s = new TH1F("hInvMassK0s", " ; M_{#pi#pi} (GeV/c^{2})", 200, 0.4, 0.6);
  fHistInvMassLcK0sp = new TH1F("hInvMassLcK0sp", " ; M_{pK^{0}_{s}} (GeV/c^{2})", 500, 1.6, 3.1);
  fHistPtLcK0sp = new TH1F("hPtLcK0sp", " ; #Lambda_{c} p_{T} (GeV/c)", 100, 0, 10.);
  fOutput->Add(fHistInvMassK0s);
  fOutput->Add(fHistInvMassLcK0sp);
  fOutput->Add(fHistPtLcK0sp);

  // MC gen level histos (for effic)
  fHistPtGenPrompt[0] = new TH1F("hPtGenPromptD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenPrompt[1] = new TH1F("hPtGenPromptDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenPrompt[2] = new TH1F("hPtGenPromptDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenPrompt[3] = new TH1F("hPtGenPromptLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenPrompt[4] = new TH1F("hPtGenPromptLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenFeeddw[0] = new TH1F("hPtGenFeeddwD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenFeeddw[1] = new TH1F("hPtGenFeeddwDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenFeeddw[2] = new TH1F("hPtGenFeeddwDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenFeeddw[3] = new TH1F("hPtGenFeeddwLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenFeeddw[4] = new TH1F("hPtGenFeeddwLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  for(Int_t i=0; i<5; i++){
    fOutput->Add(fHistPtGenPrompt[i]);
    fOutput->Add(fHistPtGenFeeddw[i]);
  }
  fHistPtGenLimAccPrompt[0] = new TH1F("hPtGenLimAccPromptD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccPrompt[1] = new TH1F("hPtGenLimAccPromptDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccPrompt[2] = new TH1F("hPtGenLimAccPromptDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccPrompt[3] = new TH1F("hPtGenLimAccPromptLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccPrompt[4] = new TH1F("hPtGenLimAccPromptLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccFeeddw[0] = new TH1F("hPtGenLimAccFeeddwD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccFeeddw[1] = new TH1F("hPtGenLimAccFeeddwDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccFeeddw[2] = new TH1F("hPtGenLimAccFeeddwDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccFeeddw[3] = new TH1F("hPtGenLimAccFeeddwLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtGenLimAccFeeddw[4] = new TH1F("hPtGenLimAccFeeddwLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  for(Int_t i=0; i<5; i++){
    fOutput->Add(fHistPtGenLimAccPrompt[i]);
    fOutput->Add(fHistPtGenLimAccFeeddw[i]);
  }
  fHistPtRecoGenPtPrompt[0] = new TH1F("hPtRecoGenPtPromptD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtPrompt[1] = new TH1F("hPtRecoGenPtPromptDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtPrompt[2] = new TH1F("hPtRecoGenPtPromptDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtPrompt[3] = new TH1F("hPtRecoGenPtPromptLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtPrompt[4] = new TH1F("hPtRecoGenPtPromptLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtFeeddw[0] = new TH1F("hPtRecoGenPtFeeddwD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtFeeddw[1] = new TH1F("hPtRecoGenPtFeeddwDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtFeeddw[2] = new TH1F("hPtRecoGenPtFeeddwDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtFeeddw[3] = new TH1F("hPtRecoGenPtFeeddwLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoGenPtFeeddw[4] = new TH1F("hPtRecoGenPtFeeddwLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoPrompt[0] = new TH1F("hPtRecoPromptD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoPrompt[1] = new TH1F("hPtRecoPromptDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoPrompt[2] = new TH1F("hPtRecoPromptDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoPrompt[3] = new TH1F("hPtRecoPromptLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoPrompt[4] = new TH1F("hPtRecoPromptLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoFeeddw[0] = new TH1F("hPtRecoFeeddwD0Kpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoFeeddw[1] = new TH1F("hPtRecoFeeddwDplusKpipi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoFeeddw[2] = new TH1F("hPtRecoFeeddwDsKKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoFeeddw[3] = new TH1F("hPtRecoFeeddwLcpKpi"," ; p_{T} (GeV/c)",100,0.,10.);
  fHistPtRecoFeeddw[4] = new TH1F("hPtRecoFeeddwLcpK0s"," ; p_{T} (GeV/c)",100,0.,10.);
  for(Int_t i=0; i<5; i++){
    fOutput->Add(fHistPtRecoGenPtPrompt[i]);
    fOutput->Add(fHistPtRecoGenPtFeeddw[i]);
    fOutput->Add(fHistPtRecoPrompt[i]);
    fOutput->Add(fHistPtRecoFeeddw[i]);
  }

  if(fEnableCPUTimeCheck) {
    TString units = "s";
    if(fCountTimeInMilliseconds)
      units = "ms";
    fHistCPUTimeCandVsNTracks = new TH2F("hCPUTimeCandPerEventVsNTracks", Form("CPU time per event elapsed for candidate selection;# of selected tracks;CPU time / event (%s);entries", units.Data()), 2500, 0., 25000, 5000, 0., 100.);
    fHistWallTimeCandVsNTracks = new TH2F("hWallTimeCandPerEventVsNTracks", Form("Wall time per event elapsed for candidate selection;# of selected tracks;Wall time / event (%s);entries", units.Data()), 2500, 0., 25000, 5000, 0., 100.);
    fHistCPUTimeTrackVsNTracks = new TH2F("hCPUTimeTrackPerEventVsNTracks", Form("CPU time per event elapsed for track selection;# of selected tracks;CPU time / event (%s);entries", units.Data()), 2500, 0., 25000, 5000, 0., 1.);
    fHistWallTimeTrackVsNTracks = new TH2F("hWallTimeTrackPerEventVsNTracks", Form("Wall time per event elapsed for track selection;# of selected tracks;Wall time / event (%s);entries", units.Data()), 2500, 0., 25000, 5000, 0., 1.);
    fOutput->Add(fHistCPUTimeCandVsNTracks);
    fOutput->Add(fHistCPUTimeTrackVsNTracks);
    fOutput->Add(fHistWallTimeCandVsNTracks);
    fOutput->Add(fHistWallTimeTrackVsNTracks);
  }

  PostData(1,fOutput);

}
//______________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::UserExec(Option_t *)
{
  //

  AliESDEvent *esd = (AliESDEvent*) (InputEvent());
  if(!esd) {
    printf("AliAnalysisTaskHFSimpleVertices::UserExec(): bad ESD\n");
    return;
  }

  // AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  // AliPIDResponse *pidResp=inputHandler->GetPIDResponse();

  Double_t centr=-1;
  AliMultSelection* mulSel=0x0;
  if(fSelectOnCentrality){
    mulSel = (AliMultSelection*)esd->FindListObject("MultSelection");
    if(mulSel){
      centr=mulSel->GetMultiplicityPercentile(fCentrEstimator.Data());
    }
  }

  AliMCEvent* mcEvent = nullptr;
  if(fReadMC){
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    // check centrality and z vertex position, which should be selected before counting the generated signals
    if(!fSelectOnCentrality || (mulSel && centr>=fMinCentrality && centr<=fMaxCentrality)){
      const AliVVertex* mcVert=mcEvent->GetPrimaryVertex();
      if(!mcVert){
        Printf("ERROR: generated vertex not available");
        return;
      }
      if(TMath::Abs(mcVert->GetZ()) < fMaxZVert){
        Int_t nParticles=mcEvent->GetNumberOfTracks();
        for (Int_t i=0;i<nParticles;i++){
          AliMCParticle* mcPart = (AliMCParticle*)mcEvent->GetTrack(i);
          TParticle* part = (TParticle*)mcEvent->Particle(i);
          if(!part || !mcPart) continue;
          Int_t absPdg=TMath::Abs(part->GetPdgCode());
          //      Int_t pdg=part->GetPdgCode();
          Int_t iPart=-1;
          Int_t deca=0;
          Bool_t isGoodDecay=kFALSE;
          Int_t dummy[4];
          if(absPdg==421){
            iPart=0;
            deca=AliVertexingHFUtils::CheckD0Decay(mcEvent,i,dummy);
            if(deca==1) isGoodDecay=kTRUE;
          }else if(absPdg==411){
            iPart=1;
            deca=AliVertexingHFUtils::CheckDplusDecay(mcEvent,i,dummy);
            if(deca>0) isGoodDecay=kTRUE;
          }else if(absPdg==431){
            iPart=2;
            deca=AliVertexingHFUtils::CheckDsDecay(mcEvent,i,dummy);
            if(deca==1) isGoodDecay=kTRUE;
          }else if(absPdg==4122){
            iPart=3;
            deca=AliVertexingHFUtils::CheckLcpKpiDecay(mcEvent,i,dummy);
            if(deca>0){
              iPart=3;
              isGoodDecay=kTRUE;
            }else{
              deca=AliVertexingHFUtils::CheckLcV0bachelorDecay(mcEvent,i,dummy);
              if(deca==1){
                iPart=4;
                isGoodDecay=kTRUE;
              }
            }
          }
          if(isGoodDecay){
            Double_t ptgen=part->Pt();
            Double_t ygen=part->Y();
            Int_t isFromB=AliVertexingHFUtils::CheckOrigin(mcEvent,mcPart,kTRUE);
            if(isFromB==4){
              fHistPtGenPrompt[iPart]->Fill(ptgen);
              if(TMath::Abs(ygen)<0.5) fHistPtGenLimAccPrompt[iPart]->Fill(ptgen);
            }else if(isFromB==5){
              fHistPtGenFeeddw[iPart]->Fill(ptgen);
              if(TMath::Abs(ygen)<0.5) fHistPtGenLimAccFeeddw[iPart]->Fill(ptgen);
            }
          }
        }
      }
    }
  }

  Int_t pdgD0dau[2]={321,211};
  Int_t pdgJpsidau[2]={11,11};



  fHistNEvents->Fill(0);
  if(fUsePhysSel){
    Bool_t isPhysSel = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerMask);
    if(!isPhysSel) return;
  }
  fHistNEvents->Fill(1);

  if(fSelectOnCentrality){
    if(centr<fMinCentrality || centr>fMaxCentrality) return;
  }
  fHistNEvents->Fill(2);

  AliESDVertex* primVtxTrk = (AliESDVertex*)esd->GetPrimaryVertex();
  AliESDVertex* primVtxSPD = (AliESDVertex*)esd->GetPrimaryVertexSPD();
  TString titTrc=primVtxTrk->GetTitle();
  if(titTrc.IsNull())return;
  if (primVtxTrk->IsFromVertexer3D() || primVtxTrk->IsFromVertexerZ()) return;
  if (primVtxTrk->GetNContributors() < 2) return;
  if (primVtxSPD->GetNContributors()<1) return;
  fHistNEvents->Fill(3);

  double covTrc[6],covSPD[6];
  primVtxTrk->GetCovarianceMatrix(covTrc);
  primVtxSPD->GetCovarianceMatrix(covSPD);
  double dz = primVtxTrk->GetZ()-primVtxSPD->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (fCutOnSPDVsTrackVtx && (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20)) return; // bad vertexing
  fHistNEvents->Fill(4);

  Float_t zvert=primVtxTrk->GetZ();
  if(TMath::Abs(zvert)>fMaxZVert) return;
  fHistNEvents->Fill(5);

  fHistPrimVertX->Fill(primVtxTrk->GetX());
  fHistPrimVertY->Fill(primVtxTrk->GetY());
  fHistPrimVertZ->Fill(primVtxTrk->GetZ());

  AliAODVertex *vertexAODp = ConvertToAODVertex(primVtxTrk);
  Double_t bzkG = (Double_t)esd->GetMagneticField();
  Int_t totTracks = TMath::Min(fMaxTracksToProcess, esd->GetNumberOfTracks());
  Double_t d0track[2],covd0track[3];


  if(!fVertexerTracks){
    fVertexerTracks = new AliVertexerTracks(bzkG);
  }else{
    Double_t oldField=fVertexerTracks->GetFieldkG();
    if(oldField!=bzkG) fVertexerTracks->SetFieldkG(bzkG);
  }
  fO2Vertexer2Prong.setBz(bzkG);
  fO2Vertexer2Prong.setPropagateToPCA(fVertexerPropagateToPCA);
  fO2Vertexer2Prong.setMaxR(fVertexerMaxR);
  fO2Vertexer2Prong.setMaxDZIni(fVertexerMaxDZIni);
  fO2Vertexer2Prong.setMinParamChange(fVertexerMinParamChange);
  fO2Vertexer2Prong.setMinRelChi2Change(fVertexerMinRelChi2Change);
  fO2Vertexer2Prong.setUseAbsDCA(fVertexerUseAbsDCA);

  fO2Vertexer3Prong.setBz(bzkG);
  fO2Vertexer3Prong.setPropagateToPCA(fVertexerPropagateToPCA);
  fO2Vertexer3Prong.setMaxR(fVertexerMaxR);
  fO2Vertexer3Prong.setMaxDZIni(fVertexerMaxDZIni);
  fO2Vertexer3Prong.setMinParamChange(fVertexerMinParamChange);
  fO2Vertexer3Prong.setMinRelChi2Change(fVertexerMinRelChi2Change);
  fO2Vertexer3Prong.setUseAbsDCA(fVertexerUseAbsDCA);


  std::clock_t clockStartTrack;
  std::chrono::time_point<std::chrono::high_resolution_clock> startTimeTrack;
  if(fEnableCPUTimeCheck) {
    clockStartTrack = std::clock();
    startTimeTrack = std::chrono::high_resolution_clock::now();
  }

  // Apply single track cuts and flag them
  int nTracksSel = 0;
  UChar_t* status = new UChar_t[totTracks];
  for (Int_t iTrack = 0; iTrack < totTracks; iTrack++) {
    status[iTrack] = 0;
    AliESDtrack* track = esd->GetTrack(iTrack);
    track->PropagateToDCA(primVtxTrk, bzkG, 100., d0track, covd0track);
    fHistPtAllTracks->Fill(track->Pt());
    fHistTglAllTracks->Fill(track->GetTgl());
    fHistEtaAllTracks->Fill(track->Eta());
    fHistImpParAllTracks->Fill(d0track[0]);
    fHistITSmapAllTracks->Fill(track->GetITSClusterMap());
    status[iTrack] = SingleTrkCuts(track,primVtxTrk,bzkG,d0track);
    if(status[iTrack] > 0)
      nTracksSel++;
    if(status[iTrack] == 0) fHistTrackStatus->Fill(0);
    if(status[iTrack] & 1) fHistTrackStatus->Fill(1);
    if(status[iTrack] & 2) fHistTrackStatus->Fill(2);
    if(status[iTrack] & 4) fHistTrackStatus->Fill(3);
  }

  std::clock_t clockStartCand;
  std::chrono::time_point<std::chrono::high_resolution_clock> startTimeCand;
  if(fEnableCPUTimeCheck) {
    // count time elapsed for track selection
    auto endTimeTrack = std::chrono::high_resolution_clock::now();
    auto clockEndTrack = std::clock();
    std::chrono::duration<double> wallTimeTrack = std::chrono::duration_cast<std::chrono::duration<double> >(endTimeTrack - startTimeTrack);
    double multFact = 1.;
    if(fCountTimeInMilliseconds)
        multFact = 1000.;
    fHistCPUTimeTrackVsNTracks->Fill(nTracksSel, double(clockEndTrack-clockStartTrack) / CLOCKS_PER_SEC * multFact);
    fHistWallTimeTrackVsNTracks->Fill(nTracksSel, wallTimeTrack.count() * multFact);

    // start time for candidate selection
    clockStartCand = std::clock();
    startTimeCand = std::chrono::high_resolution_clock::now();
  }

  Double_t d02[2]={0.,0.};
  Double_t d03[3] = {0., 0., 0.};
  AliAODRecoDecay* rd4massCalc2 = new AliAODRecoDecay(0x0, 2, 0, d02);
  AliAODRecoDecay* rd4massCalc3 = new AliAODRecoDecay(0x0, 3, 1, d03);
  TObjArray* twoTrackArray = new TObjArray(2);
  TObjArray *twoTrackArrayCasc = new TObjArray(2);
  TObjArray* threeTrackArray = new TObjArray(3);
  Double_t mom0[3], mom1[3], mom2[3];
  Int_t nv0 = esd->GetNumberOfV0s();
  for (Int_t iPosTrack_0 = 0; iPosTrack_0 < totTracks; iPosTrack_0++) {
    AliESDtrack* track_p0 = esd->GetTrack(iPosTrack_0);
    track_p0->GetPxPyPz(mom0);
    if (status[iPosTrack_0] == 0) continue;
    track_p0->PropagateToDCA(primVtxTrk, bzkG, 100., d0track, covd0track);
    fHistPtSelTracks->Fill(track_p0->Pt());
    fHistTglSelTracks->Fill(track_p0->GetTgl());
    if (status[iPosTrack_0] & 1){
      fHistImpParSelTracks2prong->Fill(d0track[0]);
      fHistEtaSelTracks2prong->Fill(track_p0->Eta());
    }
    if (status[iPosTrack_0] & 2){
      fHistImpParSelTracks3prong->Fill(d0track[0]);
      fHistEtaSelTracks3prong->Fill(track_p0->Eta());
    }
    fHistITSmapSelTracks->Fill(track_p0->GetITSClusterMap());
    if (status[iPosTrack_0] & 4){ // good bachelor track
      for(Int_t iv0=0; iv0<nv0; iv0++){
        AliESDv0 *v0 = esd->GetV0(iv0);
        if (!v0) continue;
        Bool_t onFlyStatus=v0->GetOnFlyStatus();
        if(onFlyStatus==kTRUE) continue;
        if(v0->Pt()<fMinPtV0) continue;
        UInt_t labPos = (UInt_t)TMath::Abs(v0->GetPindex());
        UInt_t labNeg = (UInt_t)TMath::Abs(v0->GetNindex());
        AliESDtrack *posVV0track=esd->GetTrack(labPos);
        AliESDtrack *negVV0track=esd->GetTrack(labNeg);
        if( !posVV0track || !negVV0track ) continue;
        // bachelor must not be a v0-track
        if (posVV0track->GetID() == track_p0->GetID() ||
            negVV0track->GetID() == track_p0->GetID()) continue;
        // reject like-sign v0
        if ( posVV0track->Charge() == negVV0track->Charge() ) continue;
        // avoid ghost TPC tracks
        if(!(posVV0track->GetStatus() & AliESDtrack::kTPCrefit) ||
           !(negVV0track->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        //  reject kinks (only necessary on AliESDtracks)
        if (posVV0track->GetKinkIndex(0)>0  || negVV0track->GetKinkIndex(0)>0) continue;
        // track cuts on V0 daughters
        if(!fTrackCutsV0Dau->AcceptTrack(posVV0track) || !fTrackCutsV0Dau->AcceptTrack(negVV0track)) continue;
        // selections on K0s
        if(!SelectV0(v0,primVtxTrk)) continue;
        Double_t xyz[3], pxpypz[3];
        v0->XvYvZv(xyz);
        v0->PxPyPz(pxpypz);
        Double_t cv[21]; for(int i=0; i<21; i++) cv[i]=0;
        AliNeutralTrackParam *trackV0 = new AliNeutralTrackParam(xyz,pxpypz,cv,0);
        twoTrackArrayCasc->AddAt(track_p0,0);
        twoTrackArrayCasc->AddAt(trackV0,1);
        AliESDVertex* vertexCasc = 0x0;
        if(fFindVertexForCascades){
          vertexCasc = ReconstructSecondaryVertex(twoTrackArrayCasc, primVtxTrk);
        }else{
          // assume Cascade decays at the primary vertex
          Double_t pos[3],cov[6],chi2perNDF;
          primVtxTrk->GetXYZ(pos);
          primVtxTrk->GetCovMatrix(cov);
          chi2perNDF = primVtxTrk->GetChi2toNDF();
          vertexCasc = new AliESDVertex(pos,cov,chi2perNDF,2);
        }
        if (vertexCasc == 0x0) {
          delete trackV0;
          twoTrackArrayCasc->Clear();
          continue;
        }
        AliAODVertex* vertexAOD = ConvertToAODVertex(vertexCasc);
        AliAODRecoCascadeHF* theCascade=MakeCascade(twoTrackArrayCasc, vertexAOD, bzkG);
        Double_t ptLcK0sp=theCascade->Pt();
        Double_t invMassLcK0sp=theCascade->InvMassLctoK0sP();
        fHistInvMassLcK0sp->Fill(invMassLcK0sp);
        fHistPtLcK0sp->Fill(ptLcK0sp);
        delete trackV0;
        twoTrackArrayCasc->Clear();
        delete theCascade;
        delete vertexAOD;
        delete vertexCasc;
      }
    }
    if (track_p0->Charge() < 0) continue;
    for (Int_t iNegTrack_0 = 0; iNegTrack_0 < totTracks; iNegTrack_0++) {
      AliESDtrack* track_n0 = esd->GetTrack(iNegTrack_0);
      track_n0->GetPxPyPz(mom1);
      if (track_n0->Charge() > 0) continue;
      if ((status[iPosTrack_0] & 1) == 0) continue;
      if ((status[iNegTrack_0] & 1) == 0) continue;
      twoTrackArray->AddAt(track_p0, 0);
      twoTrackArray->AddAt(track_n0, 1);

      AliESDVertex* trkv = ReconstructSecondaryVertex(twoTrackArray, primVtxTrk);
      if (trkv == 0x0) {
        twoTrackArray->Clear();
        continue;
      }

      double deltax = trkv->GetX() - primVtxTrk->GetX();
      double deltay = trkv->GetY() - primVtxTrk->GetY();
      double deltaz = trkv->GetZ() - primVtxTrk->GetZ();
      double decaylength = TMath::Sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
      double decaylengthxy = TMath::Sqrt(deltax * deltax + deltay * deltay);

      if(SelectInvMassAndPt2prong(twoTrackArray,rd4massCalc2)>0){
        AliAODVertex* vertexAOD = ConvertToAODVertex(trkv);
        AliAODRecoDecayHF2Prong* the2Prong = Make2Prong(twoTrackArray, vertexAOD, bzkG);
        the2Prong->SetOwnPrimaryVtx(vertexAODp);
  Int_t dzeroSel = 3;Int_t jpsiSel = 3;
        if(fCandidateCutLevel == 2){
          if(fSelectD0 + fSelectD0bar > 0){
            dzeroSel = DzeroSelectionCuts(the2Prong);
          }
          if(fSelectJpsi>0){
            jpsiSel=JpsiSelectionCuts(the2Prong,track_p0,track_n0,primVtxTrk,bzkG);
          }
        }else if(fCandidateCutLevel == 1){
          dzeroSel = DzeroSkimCuts(the2Prong);
          jpsiSel=JpsiSkimCuts(the2Prong);
        }else if(fCandidateCutLevel == 3){
          if(fSelectJpsi>0){
            jpsiSel=JpsiSkimCuts(the2Prong)+JpsiSelectionCuts(the2Prong,track_p0,track_n0,primVtxTrk,bzkG);
          }
        }

        Double_t ptD = the2Prong->Pt();
        Double_t rapid = the2Prong->Y(421);
        if(dzeroSel>0 && IsInFiducialAcceptance(ptD,rapid)) {
          Double_t m0 = the2Prong->InvMassD0();
          Double_t m0b = the2Prong->InvMassD0bar();
          Double_t ptDau0 = the2Prong->PtProng(0);
          Double_t ptDau1 = the2Prong->PtProng(1);
          Double_t ipDau0 = the2Prong->Getd0Prong(0);
          Double_t ipDau1 = the2Prong->Getd0Prong(1);
          Double_t d0xd0 = the2Prong->Prodd0d0();
          if (fSelectD0 == 0 || dzeroSel == 1 || dzeroSel == 3) fHistInvMassD0->Fill(m0);
          if (fSelectD0bar == 0 || dzeroSel == 2 || dzeroSel == 3) fHistInvMassD0->Fill(m0b);
          fHistPtD0->Fill(ptD);
          fHistYPtD0->Fill(ptD,rapid);
          fHistPtD0Dau0->Fill(ptDau0);
          fHistPtD0Dau1->Fill(ptDau1);
          fHistImpParD0Dau0->Fill(ipDau0);
          fHistImpParD0Dau1->Fill(ipDau1);
          fHistd0Timesd0->Fill(d0xd0);
          fHistCosPointD0->Fill(the2Prong->CosPointingAngle());
          fHistDecLenD0->Fill(decaylength);
          fHistDecLenXYD0->Fill(decaylengthxy);
          fHistImpParErrD0Dau->Fill(the2Prong->Getd0errProng(0));
          fHistImpParErrD0Dau->Fill(the2Prong->Getd0errProng(1));
          fHistDecLenErrD0->Fill(the2Prong->DecayLengthError());
          fHistDecLenXYErrD0->Fill(the2Prong->DecayLengthXYError());
          Double_t covMatrix[6];
          the2Prong->GetPrimaryVtx()->GetCovMatrix(covMatrix);
          fHistCovMatPrimVXX2Prong->Fill(covMatrix[0]);
          the2Prong->GetSecondaryVtx()->GetCovMatrix(covMatrix);
          fHistCovMatSecVXX2Prong->Fill(covMatrix[0]);
          fHist2ProngVertX->Fill(trkv->GetX());
          fHist2ProngVertY->Fill(trkv->GetY());
          fHist2ProngVertZ->Fill(trkv->GetZ());
          if(fReadMC && mcEvent){
            Int_t labD=MatchToMC(the2Prong,421,mcEvent,2,twoTrackArray,pdgD0dau);
            if(labD>=0){
              fHistD0SignalVertX->Fill(trkv->GetX());
              fHistD0SignalVertY->Fill(trkv->GetY());
              fHistD0SignalVertZ->Fill(trkv->GetZ());
              AliESDtrack* trDau0=(AliESDtrack*)twoTrackArray->UncheckedAt(0);
              Int_t labelDau0=TMath::Abs(trDau0->GetLabel());
              AliMCParticle* partDau0 = (AliMCParticle*)mcEvent->GetTrack(labelDau0);
              Int_t pdgCode = TMath::Abs(partDau0->PdgCode());
              if(pdgCode==211){
                fHistInvMassD0Signal->Fill(m0);
                fHistInvMassD0Refl->Fill(m0b);
              }else if(pdgCode==321){
                fHistInvMassD0Signal->Fill(m0b);
                fHistInvMassD0Refl->Fill(m0);
              }
              AliMCParticle* dmes = (AliMCParticle*)mcEvent->GetTrack(labD);
              if(dmes){
                Int_t orig=AliVertexingHFUtils::CheckOrigin(mcEvent,dmes,kTRUE);
                Double_t ptgen=dmes->Pt();
                if(orig==4){
                  fHistPtRecoPrompt[0]->Fill(ptD);
                  fHistPtRecoGenPtPrompt[0]->Fill(ptgen);
                }else if(orig==5){
                  fHistPtRecoFeeddw[0]->Fill(ptD);
                  fHistPtRecoGenPtFeeddw[0]->Fill(ptgen);
                }
              }
            }
          }
        }

        if((fCandidateCutLevel>2 && jpsiSel>1)||(fCandidateCutLevel<= 2 && jpsiSel>0) ) {
          Double_t m0 = the2Prong->InvMassJPSIee();
          Double_t ptD = the2Prong->Pt();
          Double_t ptDau0 = the2Prong->PtProng(0);
          Double_t ptDau1 = the2Prong->PtProng(1);
          Double_t ipDau0 = the2Prong->Getd0Prong(0);
          Double_t ipDau1 = the2Prong->Getd0Prong(1);
          Double_t d0xd0 = the2Prong->Prodd0d0();
          //  Printf("jpsiSel xxx xxx d0xd0 xxx %f ipDau0 xxx  %f ipDau1 xxx  %f\n",d0xd0,ipDau0,ipDau1);
          if (fSelectJpsi == 0 || jpsiSel == 1 ||jpsiSel == 2) fHistInvMassJpsi->Fill(m0);
          fHistPtJpsi->Fill(ptD);
          fHistPtJpsiDau0->Fill(ptDau0);
          fHistPtJpsiDau1->Fill(ptDau1);
          fHistImpParJpsiDau0->Fill(ipDau0);
          fHistImpParJpsiDau1->Fill(ipDau1);
          fHistd0Timesd0Jpsi->Fill(d0xd0);
          fHistCosPointJpsi->Fill(the2Prong->CosPointingAngle());
          fHistDecLenJpsi->Fill(decaylength);
          fHistDecLenXYJpsi->Fill(decaylengthxy);
          fHistDecLenErrJpsi->Fill(the2Prong->DecayLengthError());
          fHistDecLenXYErrJpsi->Fill(the2Prong->DecayLengthXYError());
          Double_t covMatrixj[6];
          the2Prong->GetPrimaryVtx()->GetCovMatrix(covMatrixj);
          fHistCovMatPrimVXX2Prong->Fill(covMatrixj[0]);
          the2Prong->GetSecondaryVtx()->GetCovMatrix(covMatrixj);
          fHistCovMatSecVXX2Prong->Fill(covMatrixj[0]);
          fHist2ProngVertX->Fill(trkv->GetX());
          fHist2ProngVertY->Fill(trkv->GetY());
          fHist2ProngVertZ->Fill(trkv->GetZ());
          if(fReadMC && mcEvent){
            Int_t labJ=MatchToMC(the2Prong,443,mcEvent,2,twoTrackArray,pdgJpsidau);
            if(labJ>=0){
              fHistJpsiSignalVertX->Fill(trkv->GetX());
              fHistJpsiSignalVertY->Fill(trkv->GetY());
              fHistJpsiSignalVertZ->Fill(trkv->GetZ());
              AliESDtrack* trDau0=(AliESDtrack*)twoTrackArray->UncheckedAt(0);
              Int_t labelDau0=TMath::Abs(trDau0->GetLabel());
              AliMCParticle* partDau0 = (AliMCParticle*)mcEvent->GetTrack(labelDau0);
              Int_t pdgCode = TMath::Abs(partDau0->PdgCode());
              if(pdgCode==11){
                fHistInvMassJpsiSignal->Fill(m0);
              }
            }
          }
        }
        delete the2Prong;
        delete vertexAOD;
      }
      delete trkv;
      if (fDo3Prong) {
        if((status[iPosTrack_0] & 2) == 0) continue;
        if((status[iNegTrack_0] & 2) == 0) continue;
        for (Int_t iPosTrack_1 = iPosTrack_0 + 1; iPosTrack_1 < totTracks; iPosTrack_1++) {
          AliESDtrack* track_p1 = esd->GetTrack(iPosTrack_1);
          if (!track_p1) continue;
          if (track_p1->Charge() < 0) continue;
          track_p1->GetPxPyPz(mom2);
          if ((status[iPosTrack_1] & 2) == 0) continue;
          // order tracks according to charge: +-+
          threeTrackArray->AddAt(track_p0, 0);
          threeTrackArray->AddAt(track_n0, 1);
          threeTrackArray->AddAt(track_p1, 2);
          ProcessTriplet(threeTrackArray, rd4massCalc3,primVtxTrk,vertexAODp,bzkG,decaylength,mcEvent);
          threeTrackArray->Clear();
        }
        for (Int_t iNegTrack_1 = iNegTrack_0 + 1; iNegTrack_1 < totTracks; iNegTrack_1++) {
          AliESDtrack* track_n1 = esd->GetTrack(iNegTrack_1);
          if (!track_n1) continue;
          if (track_n1->Charge() > 0) continue;
          if ((status[iNegTrack_1] & 2) == 0) continue;
          track_n1->GetPxPyPz(mom2);
          // order tracks according to charge: -+-
          threeTrackArray->AddAt(track_n0, 0);
          threeTrackArray->AddAt(track_p0, 1);
          threeTrackArray->AddAt(track_n1, 2);
          ProcessTriplet(threeTrackArray, rd4massCalc3,primVtxTrk,vertexAODp,bzkG,decaylength,mcEvent);
          threeTrackArray->Clear();
        }
      }
      twoTrackArray->Clear();
    }
    //  delete vertexAODp;
  }
  delete[] status;
  delete twoTrackArray;
  delete twoTrackArrayCasc;
  delete threeTrackArray;
  delete rd4massCalc2;
  delete rd4massCalc3;
  delete vertexAODp;

  if(fEnableCPUTimeCheck) {
    auto endTimeCand = std::chrono::high_resolution_clock::now();
    auto clockEndCand = std::clock();
    std::chrono::duration<double> wallTimeCand = std::chrono::duration_cast<std::chrono::duration<double> >(endTimeCand - startTimeCand);
    double multFact = 1.;
    if(fCountTimeInMilliseconds)
        multFact = 1000.;
    fHistCPUTimeCandVsNTracks->Fill(nTracksSel, double(clockEndCand-clockStartCand) / CLOCKS_PER_SEC * multFact);
    fHistWallTimeCandVsNTracks->Fill(nTracksSel, wallTimeCand.count() * multFact);
  }

  PostData(1,fOutput);

}


//______________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::ProcessTriplet(TObjArray* threeTrackArray, AliAODRecoDecay* rd4massCalc3, AliESDVertex* primVtxTrk, AliAODVertex *vertexAODp, float bzkG, double dist12, AliMCEvent* mcEvent){

  Int_t pdgDplusdau[3]={321,211,211};
  Int_t pdgDsdau[3]={321,321,211};

  Int_t massSel = SelectInvMassAndPt3prong(threeTrackArray, rd4massCalc3);
  if (massSel == 0) {
    threeTrackArray->Clear();
    return;
  }
  AliESDVertex* trkv3 = ReconstructSecondaryVertex(threeTrackArray, primVtxTrk);
  if (trkv3 == 0x0) {
    threeTrackArray->Clear();
    return;
  }
  AliAODVertex* vertexAOD3 = ConvertToAODVertex(trkv3);
  AliAODRecoDecayHF3Prong* the3Prong = Make3Prong(threeTrackArray, vertexAOD3, bzkG);
  the3Prong->SetOwnPrimaryVtx(vertexAODp);
  Double_t ptcand_3prong = the3Prong->Pt();
  if (ptcand_3prong < fMinPt3Prong){
    delete the3Prong;
    delete vertexAOD3;
    return;
  }
  Double_t covMatrixPV[6], covMatrixSV[6];
  the3Prong->GetPrimaryVtx()->GetCovMatrix(covMatrixPV);
  the3Prong->GetSecondaryVtx()->GetCovMatrix(covMatrixSV);

  if (massSel & (1 << kbitDplus)) {
    Int_t dplusSel = 1;
    if(fCandidateCutLevel == 2 && fSelectDplus > 0){
      dplusSel = DplusSelectionCuts(the3Prong, bzkG);
    }else if(fCandidateCutLevel == 1){
      dplusSel = DplusSkimCuts(the3Prong);
    }
    Double_t rapid = the3Prong->Y(411);
    if(dplusSel>0 && IsInFiducialAcceptance(ptcand_3prong,rapid)) {
      Double_t sqSumd0Prong = 0;
      for(Int_t iProng = 0; iProng < 3; iProng++)
        sqSumd0Prong += the3Prong->Getd0Prong(iProng) * the3Prong->Getd0Prong(iProng);
      Double_t mplus=the3Prong->InvMassDplus();
      fHistInvMassDplus->Fill(mplus);
      fHistPtDplus->Fill(ptcand_3prong);
      fHistYPtDplus->Fill(ptcand_3prong,rapid);
      fHistPtDplusDau0->Fill(the3Prong->PtProng(0));
      fHistPtDplusDau1->Fill(the3Prong->PtProng(1));
      fHistPtDplusDau2->Fill(the3Prong->PtProng(2));
      fHistImpParDplusDau0->Fill(the3Prong->Getd0Prong(0));
      fHistImpParDplusDau1->Fill(the3Prong->Getd0Prong(1));
      fHistImpParDplusDau2->Fill(the3Prong->Getd0Prong(2));
      fHistDecLenDplus->Fill(the3Prong->DecayLength());
      fHistDecLenXYDplus->Fill(the3Prong->DecayLengthXY());
      fHistNormDecLenXYDplus->Fill(the3Prong->NormalizedDecayLengthXY());
      fHistImpParErrDplusDau->Fill(the3Prong->Getd0errProng(0));
      fHistImpParErrDplusDau->Fill(the3Prong->Getd0errProng(1));
      fHistImpParErrDplusDau->Fill(the3Prong->Getd0errProng(2));
      fHistDecLenErrDplus->Fill(the3Prong->DecayLengthError());
      fHistDecLenXYErrDplus->Fill(the3Prong->DecayLengthXYError());
      fHistCosPointDplus->Fill(the3Prong->CosPointingAngle());
      fHistCosPointXYDplus->Fill(the3Prong->CosPointingAngleXY());
      fHistImpParXYDplus->Fill(the3Prong->ImpParXY());
      fHistNormIPDplus->Fill(AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(the3Prong, bzkG));
      fHistoSumSqImpParDplusDau->Fill(sqSumd0Prong);
      fHistCovMatPrimVXX3Prong->Fill(covMatrixPV[0]);
      fHistCovMatSecVXX3Prong->Fill(covMatrixSV[0]);
      fHist3ProngVertX->Fill(trkv3->GetX());
      fHist3ProngVertY->Fill(trkv3->GetY());
      fHist3ProngVertZ->Fill(trkv3->GetZ());
      if(fReadMC && mcEvent){
        Int_t labD=MatchToMC(the3Prong,411,mcEvent,3,threeTrackArray,pdgDplusdau);
        if(labD>=0){
          fHistInvMassDplusSignal->Fill(mplus);
          AliMCParticle* dmes = (AliMCParticle*)mcEvent->GetTrack(labD);
          if(dmes){
            Int_t orig=AliVertexingHFUtils::CheckOrigin(mcEvent,dmes,kTRUE);
            Double_t ptgen=dmes->Pt();
            if(orig==4){
              fHistPtRecoGenPtPrompt[1]->Fill(ptgen);
              fHistPtRecoPrompt[1]->Fill(ptcand_3prong);
            }else if(orig==5){
              fHistPtRecoGenPtFeeddw[1]->Fill(ptgen);
              fHistPtRecoFeeddw[1]->Fill(ptcand_3prong);
            }
          }
        }
      }
    }
  }
  if (massSel & (1 << kbitDs)) {
    Int_t dsSel = 3;
    if(fCandidateCutLevel >= 1){
      dsSel = DsSkimCuts(the3Prong);
    }
    Double_t rapid = the3Prong->Y(431);
    if(dsSel>0 && IsInFiducialAcceptance(ptcand_3prong,rapid)) {
      Double_t mKKpi=the3Prong->InvMassDsKKpi();
      Double_t mpiKK=the3Prong->InvMassDspiKK();
      if( dsSel==1 || dsSel==3) fHistInvMassDs->Fill(mKKpi);
      if( dsSel==2 || dsSel==3) fHistInvMassDs->Fill(mpiKK);
      fHistPtDs->Fill(ptcand_3prong);
      fHistYPtDs->Fill(ptcand_3prong,rapid);
      fHistDecLenDs->Fill(the3Prong->DecayLength());
      fHistCosPointDs->Fill(the3Prong->CosPointingAngle());
      fHistCovMatPrimVXX3Prong->Fill(covMatrixPV[0]);
      fHistCovMatSecVXX3Prong->Fill(covMatrixSV[0]);
      fHist3ProngVertX->Fill(trkv3->GetX());
      fHist3ProngVertY->Fill(trkv3->GetY());
      fHist3ProngVertZ->Fill(trkv3->GetZ());
      if(fReadMC && mcEvent){
        Int_t labD=MatchToMC(the3Prong,431,mcEvent,3,threeTrackArray,pdgDsdau);
        if(labD>=0){
          AliMCParticle* dmes = (AliMCParticle*)mcEvent->GetTrack(labD);
          if(dmes){
            AliESDtrack* trDau0=(AliESDtrack*)threeTrackArray->UncheckedAt(0);
            Int_t labelDau0=TMath::Abs(trDau0->GetLabel());
            AliMCParticle* partDau0 = (AliMCParticle*)mcEvent->GetTrack(labelDau0);
            if(partDau0){
              Int_t pdgCode = TMath::Abs(partDau0->PdgCode());
              if(pdgCode==211){
                fHistInvMassDsSignal->Fill(mpiKK);
                fHistInvMassDsRefl->Fill(mKKpi);
              }else if(pdgCode==321){
                fHistInvMassDsSignal->Fill(mKKpi);
                fHistInvMassDsRefl->Fill(mpiKK);
              }
            }
            Int_t orig=AliVertexingHFUtils::CheckOrigin(mcEvent,dmes,kTRUE);
            Double_t ptgen=dmes->Pt();
            if(orig==4){
              fHistPtRecoGenPtPrompt[2]->Fill(ptgen);
              fHistPtRecoPrompt[2]->Fill(ptcand_3prong);
            }else if(orig==5){
              fHistPtRecoGenPtFeeddw[2]->Fill(ptgen);
              fHistPtRecoFeeddw[2]->Fill(ptcand_3prong);
            }
          }
        }
      }
    }
  }
  if (massSel & (1 << kbitLc)) {
    fHistDist12LcpKpi->Fill(dist12);
    Int_t lcSel = 3;
    if(fCandidateCutLevel == 2 && fSelectLcpKpi > 0){
      lcSel = LcSelectionCuts(the3Prong);
      if(dist12 <= fLcCuts[0][4]) lcSel=0; //cut on dist12
    }else if(fCandidateCutLevel == 1){
      lcSel = LcSkimCuts(the3Prong);
    }
    Double_t rapid = the3Prong->Y(4122);
    if (lcSel > 0 && IsInFiducialAcceptance(ptcand_3prong,rapid)) {
      Double_t m0LcpKpi = the3Prong->InvMassLcpKpi();
      Double_t m0LcpiKp = the3Prong->InvMassLcpiKp();
      if (fSelectLcpKpi == 0 || lcSel == 1 || lcSel == 3) fHistInvMassLc->Fill(m0LcpKpi);
      if (fSelectLcpKpi == 0 || lcSel == 2 || lcSel == 3) fHistInvMassLc->Fill(m0LcpiKp);
      fHistPtLc->Fill(ptcand_3prong);
      fHistYPtLc->Fill(ptcand_3prong,rapid);
      fHistPtLcDau0->Fill(the3Prong->PtProng(0));
      fHistPtLcDau1->Fill(the3Prong->PtProng(1));
      fHistPtLcDau2->Fill(the3Prong->PtProng(2));
      fHistDecLenLc->Fill(the3Prong->DecayLength());
      fHistCosPointLc->Fill(the3Prong->CosPointingAngle());
      fHistCovMatPrimVXX3Prong->Fill(covMatrixPV[0]);
      fHistCovMatSecVXX3Prong->Fill(covMatrixSV[0]);
      fHist3ProngVertX->Fill(trkv3->GetX());
      fHist3ProngVertY->Fill(trkv3->GetY());
      fHist3ProngVertZ->Fill(trkv3->GetZ());
    }
  }
  delete trkv3;
  delete the3Prong;
  delete vertexAOD3;
  return;
}


//______________________________________________________________________________
Bool_t AliAnalysisTaskHFSimpleVertices::GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float bzkG)
{
  /// fast calculation (no covariance matrix treatment) of track momentum at secondary vertex

  Double_t alpha = tr->GetAlpha();
  Double_t sn = TMath::Sin(alpha), cs = TMath::Cos(alpha);
  Double_t x = tr->GetX(), y = tr->GetParameter()[0], snp = tr->GetParameter()[2];
  Double_t xv = secVert->GetX() * cs + secVert->GetY() * sn;
  Double_t yv = -secVert->GetX() * sn + secVert->GetY() * cs;
  x -= xv;
  y -= yv;
  Double_t crv = tr->GetC(bzkG);
  if (TMath::Abs(bzkG) < 0.000001)
    crv = 0.;
  double csp = TMath::Sqrt((1. - snp) * (1. + snp));

  Double_t tgfv = -(crv * x - snp) / (crv * y + csp);
  cs = 1. / TMath::Sqrt(1 + tgfv * tgfv);
  sn = cs < 1. ? tgfv * cs : 0.;

  x = xv * cs + yv * sn;
  Double_t alpNew = alpha + TMath::ASin(sn);
  Double_t ca = TMath::Cos(alpNew - alpha), sa = TMath::Sin(alpNew - alpha);
  Double_t p2 = tr->GetSnp();
  Double_t xNew = tr->GetX() * ca + tr->GetY() * sa;
  Double_t p2New = p2 * ca - TMath::Sqrt((1. - p2) * (1. + p2)) * sa;
  momentum[0] = tr->GetSigned1Pt();
  momentum[1] = p2New * (x - xNew) * tr->GetC(bzkG);
  momentum[2] = tr->GetTgl();
  Bool_t retCode = tr->Local2GlobalMomentum(momentum, alpNew);
  return retCode;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::SingleTrkCuts(AliESDtrack* trk, AliESDVertex* primVert, Double_t bzkG, Double_t d0track[2])
{
  Double_t covd0track[3];
  if (!trk->PropagateToDCA(primVert, bzkG, kVeryBig, d0track, covd0track)) return kFALSE;
  trk->RelateToVertex(primVert, bzkG, kVeryBig);
  Int_t retCode=0;

  // test first impact parameter
  Int_t iPtTrack = GetPtBinSingleTrack(trk->Pt());
  if(TMath::Abs(d0track[0]) >= fSingleTrackCuts2Prong[iPtTrack][0] && TMath::Abs(d0track[0]) <= fSingleTrackCuts2Prong[iPtTrack][1])
    retCode+=1;
  if(TMath::Abs(d0track[0]) >= fSingleTrackCuts3Prong[iPtTrack][0] && TMath::Abs(d0track[0]) <= fSingleTrackCuts3Prong[iPtTrack][1])
    retCode+=2;

  // test other cuts if impact-parameter tested
  Int_t retCodeCurrent=retCode;
  if((retCodeCurrent == 1 || retCodeCurrent == 3) && !fTrackCuts2pr->AcceptTrack(trk)) retCode-=1;
  if(retCodeCurrent >= 2 && !fTrackCuts3pr->AcceptTrack(trk)) retCode-=2;

  //test cuts for bachelor track
  if(fTrackCutsBach->AcceptTrack(trk)) retCode+=4;
  
  return retCode;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskHFSimpleVertices::SelectV0(AliESDv0 *v0, AliESDVertex* primvtx)
{
  // Apply V0 selections
  Double_t cpa=v0->GetV0CosineOfPointingAngle(primvtx->GetX(),primvtx->GetY(),primvtx->GetZ());
  if(cpa<fMinCosPointV0) return kFALSE;
  Double_t invMassK0s = v0->GetEffMass(2,2);
  if(TMath::Abs(invMassK0s-fMassK0s)>fCutOnK0sMass) return kFALSE;
  fHistInvMassK0s->Fill(invMassK0s);
  return kTRUE;
}
//______________________________________________________________________________
AliESDVertex* AliAnalysisTaskHFSimpleVertices::ReconstructSecondaryVertex(TObjArray* trkArray, AliESDVertex* primvtx)
{

  AliESDVertex* trkv =0x0;
  // printf("------\n");
  // for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
  //   AliExternalTrackParam *tp=(AliExternalTrackParam *)trkArray->At(jt);
  //   printf("%f ",tp->Pt());
  // }
  // printf("\n");
  if(fSecVertexerAlgo==0){
    fVertexerTracks->SetVtxStart(primvtx);
    trkv = (AliESDVertex*)fVertexerTracks->VertexForSelectedESDTracks(trkArray);
    if (trkv->GetNContributors() != trkArray->GetEntriesFast()) return 0x0;
  }else if(fSecVertexerAlgo==1){
    o2::track::TrackParCov* o2Track[3] = {nullptr};
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      o2Track[jt] = static_cast<o2::track::TrackParCov *>((AliExternalTrackParam *)trkArray->At(jt));
    }
    Int_t nVert=0;
    Double_t vertCoord[3];
    Double_t vertCov[6];
    Double_t vertChi2=-999.;
    if(trkArray->GetEntriesFast()==2){
      nVert=fO2Vertexer2Prong.process(*o2Track[0], *o2Track[1]);
      if(nVert){
        fO2Vertexer2Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer2Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer2Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer2Prong.getChi2AtPCACandidate();
      }
    }else if(trkArray->GetEntriesFast()==3){
      nVert=fO2Vertexer3Prong.process(*o2Track[0], *o2Track[1], *o2Track[2]);
      if(nVert){
        fO2Vertexer3Prong.propagateTracksToVertex();
        auto vertPos = fO2Vertexer3Prong.getPCACandidate();
        auto vertCMat = fO2Vertexer3Prong.calcPCACovMatrix().Array();
        for(Int_t ic=0; ic<3; ic++) vertCoord[ic]=vertPos[ic];
        for(Int_t ic=0; ic<6; ic++) vertCov[ic]=vertCMat[ic];
        vertChi2 = fO2Vertexer3Prong.getChi2AtPCACandidate();
      }
    }
    if(nVert) trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
    //   else printf("nVert = %d\n",nVert);
  }else if(fSecVertexerAlgo==2){
    KFParticle  dMesonVert;
    double posmom[6],cov[21];
    for(Int_t jt=0; jt<trkArray->GetEntriesFast(); jt++){
      AliExternalTrackParam* trpar=(AliExternalTrackParam*)trkArray->At(jt);
      trpar->GetXYZ(posmom);
      trpar->GetPxPyPz(posmom+3);
      trpar->GetCovarianceXYZPxPyPz(cov);
      KFParticle trKFpar;
      trKFpar.Create(posmom,cov,trpar->Charge(),0.13957); // mass of the pion for the time being
      dMesonVert.AddDaughter(trKFpar);
    }
    Double_t vertChi2 = dMesonVert.GetChi2() / dMesonVert.GetNDF();
    Double_t vertCoord[3]={dMesonVert.X(),dMesonVert.Y(),dMesonVert.Z()};
    Double_t vertCov[6];
    vertCov[0]=dMesonVert.GetCovariance(0,0);
    vertCov[1]=dMesonVert.GetCovariance(0,1);
    vertCov[2]=dMesonVert.GetCovariance(1,1);
    vertCov[3]=dMesonVert.GetCovariance(0,2);
    vertCov[4]=dMesonVert.GetCovariance(1,2);
    vertCov[5]=dMesonVert.GetCovariance(2,2);
    trkv = new AliESDVertex(vertCoord,vertCov,vertChi2,trkArray->GetEntriesFast());
  }
  if(!trkv) return 0x0;
  Double_t vertRadius2 = trkv->GetX() * trkv->GetX() + trkv->GetY() * trkv->GetY();
  if(vertRadius2>fMaxDecVertRadius2) return 0x0;
  //  trkv->Print("all");
  //  printf("=============\n");
  return trkv;
}
//______________________________________________________________________________
AliAODVertex* AliAnalysisTaskHFSimpleVertices::ConvertToAODVertex(AliESDVertex* trkv)
{
  Double_t pos[3], cov[6], chi2perNDF;
  trkv->GetXYZ(pos);       // position
  trkv->GetCovMatrix(cov); // covariance matrix
  chi2perNDF = trkv->GetChi2toNDF();
  //  double dispersion = trkv->GetDispersion();
  AliAODVertex* vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, trkv->GetNContributors());
  vertexAOD->SetNContributors(trkv->GetNContributors());
  return vertexAOD;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskHFSimpleVertices::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  if(fMaxRapidityCand>-998.){
    if(TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
    else return kTRUE;
  }
  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    if (TMath::Abs(y) > 0.8) return kFALSE;
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5;
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;
  }
  return kTRUE;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DzeroSkimCuts(AliAODRecoDecayHF2Prong* cand)
{
  bool isD0 = true;
  bool isD0bar = true;
  Double_t ptCand = cand->Pt();
  Int_t iPtDzero = GetPtBin(ptCand, fPtBinLimsDzeroSkims, kMaxNPtBins2ProngsSkims);
  if (iPtDzero < 0) return 0;
  Double_t m0=cand->InvMassD0();
  Double_t m0b=cand->InvMassD0bar();
  if(m0 < fDzeroSkimCuts[iPtDzero][0] || m0 > fDzeroSkimCuts[iPtDzero][1]) isD0=false;
  if(m0b < fDzeroSkimCuts[iPtDzero][0] || m0b > fDzeroSkimCuts[iPtDzero][1])  isD0bar=false;
  if (!isD0 && !isD0bar) return 0;
  if (cand->CosPointingAngle() < fDzeroSkimCuts[iPtDzero][2]) return 0;
  if (cand->Prodd0d0() > fDzeroSkimCuts[iPtDzero][3]) return 0;

  Int_t returnValue=0;
  if(isD0) returnValue+=1;
  if(isD0bar) returnValue+=2;
  return returnValue;

}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::JpsiSkimCuts(AliAODRecoDecayHF2Prong* cand)
{
  Double_t ptCand = cand->Pt();
  Int_t iPtJpsi = GetPtBin(ptCand, fPtBinLimsJpsiSkims, kMaxNPtBins2ProngsSkims);
  if (iPtJpsi < 0) return 0;
  Double_t m0=cand->InvMassJPSIee();
  if(m0 < fJpsiSkimCuts[iPtJpsi][0] || m0 > fJpsiSkimCuts[iPtJpsi][1]) return 0;
  if (cand->CosPointingAngle() < fJpsiSkimCuts[iPtJpsi][2]) return 0;
  if (cand->Prodd0d0() > fJpsiSkimCuts[iPtJpsi][3]) return 0;

  return 1;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::JpsiSelectionCuts(AliAODRecoDecayHF2Prong* cand,AliESDtrack* trk_p,AliESDtrack* trk_n,AliESDVertex* primvtx,float bzkG)
{
  Double_t dca[2];
  Double_t dca1[2];
  Double_t covtrack[2]; Double_t covtrack1[2];
  Double_t ptCand = cand->Pt();
  if (ptCand < fMinPtJpsi || ptCand > fMaxPtJpsi) return 0;
  Int_t jPtBin = GetPtBin(ptCand, fPtBinLimsJpsi, fNPtBinsJpsi);
  if (jPtBin==-1) return 0;
  if (TMath::Abs(cand->InvMassJPSIee()-fMassJpsi) > fJpsiCuts[jPtBin][0] ) return 0;
  trk_p->PropagateToDCA(primvtx, bzkG, 100., dca, covtrack);
  trk_n->PropagateToDCA(primvtx, bzkG, 100., dca1, covtrack1);
  if((TMath::Abs(dca[0])>fJpsiCuts[jPtBin][1]||TMath::Abs(dca[0])>fJpsiCuts[jPtBin][2])||(TMath::Abs(dca1[1])>fJpsiCuts[jPtBin][1]||TMath::Abs(dca1[1])>fJpsiCuts[jPtBin][2])) return 0;
  if (cand->Pt2Prong(0) < fJpsiCuts[jPtBin][3] || cand->Pt2Prong(1) < fJpsiCuts[jPtBin][3] ) return 0;

  return 1;

}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DzeroSelectionCuts(AliAODRecoDecayHF2Prong* cand)
{
  bool isD0 = true;
  bool isD0bar = true;
  Double_t ptCand = cand->Pt();
  if (ptCand < fMinPtDzero || ptCand > fMaxPtDzero) return 0;
  Int_t jPtBin = GetPtBin(ptCand, fPtBinLimsDzero, fNPtBinsDzero);
  if (jPtBin==-1) return 0;
  if (cand->Prodd0d0() > fDzeroCuts[jPtBin][7]) return 0;
  if (cand->CosPointingAngle() < fDzeroCuts[jPtBin][8]) return 0;
  if (cand->CosPointingAngleXY() < fDzeroCuts[jPtBin][9]) return 0;
  if (cand->NormalizedDecayLengthXY() < fDzeroCuts[jPtBin][10]) return 0;
  Double_t decayLengthCut = TMath::Min((cand->P() * 0.0066) + 0.01, 0.06);
  if (TMath::Abs(cand->Normalizedd0Prong(0)) < 0.5 || TMath::Abs(cand->Normalizedd0Prong(1)) < 0.5) return 0;
  if (cand->DecayLength() * cand->DecayLength() < decayLengthCut * decayLengthCut) return 0;
  // if (cand->NormalizedDecayLength() * cand->NormalizedDecayLength() < 1.0) return 0;
  if (TMath::Abs(cand->InvMassD0()-fMassDzero) > fDzeroCuts[jPtBin][0] ) isD0=false;
  if (TMath::Abs(cand->InvMassD0bar()-fMassDzero) > fDzeroCuts[jPtBin][0] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (cand->Pt2Prong(0) < fDzeroCuts[jPtBin][4]*fDzeroCuts[jPtBin][4] || cand->Pt2Prong(1) < fDzeroCuts[jPtBin][3]*fDzeroCuts[jPtBin][3] ) isD0=false;
  if (cand->Pt2Prong(0) < fDzeroCuts[jPtBin][3]*fDzeroCuts[jPtBin][3] || cand->Pt2Prong(1) < fDzeroCuts[jPtBin][4]*fDzeroCuts[jPtBin][4] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[jPtBin][6] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[jPtBin][5] ) isD0=false;
  if (TMath::Abs(cand->Getd0Prong(0)) > fDzeroCuts[jPtBin][5] || TMath::Abs(cand->Getd0Prong(1)) > fDzeroCuts[jPtBin][6] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;

  Double_t cosThetaStarD0,cosThetaStarD0bar;
  cand->CosThetaStarD0(cosThetaStarD0,cosThetaStarD0bar);
  if (TMath::Abs(cosThetaStarD0) > fDzeroCuts[jPtBin][2] ) isD0=false;
  if (TMath::Abs(cosThetaStarD0bar) > fDzeroCuts[jPtBin][2] ) isD0bar=false;
  if (!isD0 && !isD0bar) return 0;


  Int_t returnValue=0;
  if(isD0) returnValue+=1;
  if(isD0bar) returnValue+=2;
  return returnValue;

}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DplusSkimCuts(AliAODRecoDecayHF3Prong* cand)
{
  Double_t ptCand = cand->Pt();
  Int_t iPtDplus = GetPtBin(ptCand, fPtBinLimsDplusSkims, kMaxNPtBins3ProngsSkims);
  if (iPtDplus < 0) return 0;
  Double_t m=cand->InvMassDplus();
  if(m < fDplusSkimCuts[iPtDplus][0] || m > fDplusSkimCuts[iPtDplus][1])  return 0;
  if (cand->CosPointingAngle() < fDplusSkimCuts[iPtDplus][2]) return 0;
  if (cand->DecayLength2() < fDplusSkimCuts[iPtDplus][3]*fDplusSkimCuts[iPtDplus][3]) return 0;

  return 1;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DplusSelectionCuts(AliAODRecoDecayHF3Prong* cand, Double_t bzkG)
{
  Double_t ptCand = cand->Pt();
  if (ptCand < fMinPtDplus || ptCand > fMaxPtDplus) return 0;
  Int_t jPtBin = GetPtBin(ptCand, fPtBinLimsDplus, fNPtBinsDplus);
  if (jPtBin==-1) return 0;
  if (cand->Pt2Prong(0) < fDplusCuts[jPtBin][1] * fDplusCuts[jPtBin][1] ||
      cand->Pt2Prong(1) < fDplusCuts[jPtBin][2] * fDplusCuts[jPtBin][2] ||
      cand->Pt2Prong(2) < fDplusCuts[jPtBin][1] * fDplusCuts[jPtBin][1])
      return 0;
  if (TMath::Abs(cand->InvMassDplus() - fMassDplus) > fDplusCuts[jPtBin][0]) return 0;
  if (cand->DecayLength() < fDplusCuts[jPtBin][3]) return 0;
  if (cand->NormalizedDecayLengthXY() < fDplusCuts[jPtBin][4]) return 0;
  if (cand->CosPointingAngle() < fDplusCuts[jPtBin][5]) return 0;
  if (cand->CosPointingAngleXY() < fDplusCuts[jPtBin][6]) return 0;
  Double_t dd0max = AliVertexingHFUtils::ComputeMaxd0MeasMinusExp(cand, bzkG);
  if(TMath::Abs(dd0max) > fDplusCuts[jPtBin][7]) return 0;

  return 1;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::DsSkimCuts(AliAODRecoDecayHF3Prong* cand)
{
  bool isKKpi = true;
  bool ispiKK = true;
  Double_t ptCand = cand->Pt();
  Int_t iPtDs = GetPtBin(ptCand, fPtBinLimsDsSkims, kMaxNPtBins3ProngsSkims);
  if (iPtDs < 0) return 0;
  Double_t mKKpi=cand->InvMassDsKKpi();
  Double_t mpiKK=cand->InvMassDspiKK();
  if(mKKpi < fDsSkimCuts[iPtDs][0] || mKKpi > fDsSkimCuts[iPtDs][1]) isKKpi=false;
  if(mpiKK < fDsSkimCuts[iPtDs][0] || mpiKK > fDsSkimCuts[iPtDs][1]) ispiKK=false;
  if (!isKKpi && !ispiKK) return 0;
  if (cand->CosPointingAngle() < fDsSkimCuts[iPtDs][2]) return 0;
  if (cand->DecayLength2() < fDsSkimCuts[iPtDs][3]*fDsSkimCuts[iPtDs][3]) return 0;

  Int_t returnValue=0;
  if(isKKpi) returnValue+=1;
  if(ispiKK) returnValue+=2;
  return returnValue;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::LcSkimCuts(AliAODRecoDecayHF3Prong* cand)
{
  bool ispKpi = true;
  bool ispiKp = true;
  Double_t ptCand = cand->Pt();
  Int_t iPtLc = GetPtBin(ptCand, fPtBinLimsLcSkims, kMaxNPtBins3ProngsSkims);
  if (iPtLc < 0) return 0;
  Double_t mpKpi=cand->InvMassLcpKpi();
  Double_t mpiKp=cand->InvMassLcpiKp();
  if(mpKpi < fLcSkimCuts[iPtLc][0] || mpKpi > fLcSkimCuts[iPtLc][1]) ispKpi=false;
  if(mpiKp < fLcSkimCuts[iPtLc][0] || mpiKp > fLcSkimCuts[iPtLc][1]) ispiKp=false;
  if (!ispKpi && !ispiKp) return 0;
  if (cand->CosPointingAngle() < fLcSkimCuts[iPtLc][2]) return 0;
  if (cand->DecayLength2() < fLcSkimCuts[iPtLc][3]*fLcSkimCuts[iPtLc][3]) return 0;

  Int_t returnValue=0;
  if(ispKpi) returnValue+=1;
  if(ispiKp) returnValue+=2;
  return returnValue;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::GetPtBin(Double_t ptCand, Double_t* ptBinLims, Double_t nPtBins)
{
  for (Int_t i = 0; i < nPtBins; i++) {
    if (ptCand>=ptBinLims[i] && ptCand<ptBinLims[i+1]){
      return i;
    }
  }
  return -1;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::GetPtBinSingleTrack(Double_t ptTrack)
{
  for (Int_t i = 0; i < fNPtBinsSingleTrack; i++) {
    if (ptTrack>=fPtBinLimsSingleTrack[i] && ptTrack<fPtBinLimsSingleTrack[i+1]){
      return i;
    }
  }
  return -1;
}


//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::SelectInvMassAndPt2prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc2)
{
  Int_t retval = (1 << kbitDzero) + (1 << kbitDzerobar) + (1 << kbitJpsi);
  Double_t momentum[3];
  Double_t px[2], py[2], pz[2];
  for (Int_t iTrack = 0; iTrack < 2; iTrack++) {
    AliESDtrack* track = (AliESDtrack*)trkArray->UncheckedAt(iTrack);
    track->GetPxPyPz(momentum);
    px[iTrack] = momentum[0];
    py[iTrack] = momentum[1];
    pz[iTrack] = momentum[2];
  }
  UInt_t pdg2[2];
  Int_t nprongs = 2;
  rd4massCalc2->SetPxPyPzProngs(nprongs, px, py, pz);
  Double_t ptCand = rd4massCalc2->Pt() + fPtWithoutVtxToll;
  Int_t iPtBinDzero = GetPtBin(ptCand, fPtBinLimsDzeroSkims, kMaxNPtBins2ProngsSkims);
  if(iPtBinDzero < 0) {
    retval &= ~(1 << kbitDzero);
    retval &= ~(1 << kbitDzerobar);
  }
  Int_t iPtBinJpsi = GetPtBin(ptCand, fPtBinLimsJpsiSkims, kMaxNPtBins2ProngsSkims);
  if(iPtBinJpsi < 0) {
    retval &= ~(1 << kbitJpsi);
  }
  Double_t minv2;
  Double_t lolim=fDzeroSkimCuts[iPtBinDzero][0];
  Double_t hilim=fDzeroSkimCuts[iPtBinDzero][1];
  pdg2[0]=211; pdg2[1]=321;   // pi+ K- --> D0
  minv2 = rd4massCalc2->InvMass2(nprongs,pdg2);
  if((retval & (1 << kbitDzero)) && (minv2<lolim*lolim || minv2>hilim*hilim)) retval &= ~(1 << kbitDzero);
  pdg2[0]=321; pdg2[1]=211;   // K+ pi- --> D0bar
  minv2 = rd4massCalc2->InvMass2(nprongs,pdg2);
  if((retval & (1 << kbitDzerobar)) && (minv2<lolim*lolim || minv2>hilim*hilim)) retval &= ~(1 << kbitDzerobar);
  lolim=fJpsiSkimCuts[iPtBinJpsi][0];
  hilim=fJpsiSkimCuts[iPtBinJpsi][1];
  pdg2[0]=11; pdg2[1]=11;  //e+e- -->Jpsi
  minv2 = rd4massCalc2->InvMass2(nprongs,pdg2);
  if((retval & (1 << kbitJpsi)) && (minv2<lolim*lolim || minv2>hilim*hilim)) retval &= ~(1 << kbitJpsi);
  return retval;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFSimpleVertices::SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3)
{

  Int_t retval = (1 << kbitDplus) + (1 << kbitDs) + (1 << kbitLc);
  Double_t momentum[3];
  Double_t px[3], py[3], pz[3];
  for (Int_t iTrack = 0; iTrack < 3; iTrack++) {
    AliESDtrack* track = (AliESDtrack*)trkArray->UncheckedAt(iTrack);
    track->GetPxPyPz(momentum);
    px[iTrack] = momentum[0];
    py[iTrack] = momentum[1];
    pz[iTrack] = momentum[2];
  }
  UInt_t pdg3[3];
  Int_t nprongs = 3;
  rd4massCalc3->SetPxPyPzProngs(nprongs, px, py, pz);
  Double_t ptCand = rd4massCalc3->Pt() + fPtWithoutVtxToll;
  Int_t iPtBinDplus = GetPtBin(ptCand, fPtBinLimsDplusSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinDplus < 0) {
    retval &= ~(1 << kbitDplus);
  }
  Int_t iPtBinDs = GetPtBin(ptCand, fPtBinLimsDsSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinDs < 0) {
    retval &= ~(1 << kbitDs);
  }
  Int_t iPtBinLc = GetPtBin(ptCand, fPtBinLimsLcSkims, kMaxNPtBins3ProngsSkims);
  if(iPtBinLc < 0) {
    retval &= ~(1 << kbitLc);
  }
  Double_t minv2;
  Double_t lolim=fDplusSkimCuts[iPtBinDplus][0];
  Double_t hilim=fDplusSkimCuts[iPtBinDplus][1];
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitDplus)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    retval &= ~(1 << kbitDplus);
  lolim=fDsSkimCuts[iPtBinDs][0];
  hilim=fDsSkimCuts[iPtBinDs][1];
  Bool_t isSelDs[2] = {true, true}; 
  for (Int_t ih = 0; ih < 2; ih++) {
    Int_t k = ih * 2;
    pdg3[k] = 321;
    pdg3[1] = 321;
    pdg3[2 - k] = 211;
    minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
    if ((retval & (1 << kbitDs)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
      isSelDs[ih] = false;
  }
  if(!isSelDs[0] && !isSelDs[1])
   retval &= ~(1 << kbitDs);
  lolim=fLcSkimCuts[iPtBinLc][0];
  hilim=fLcSkimCuts[iPtBinLc][1];
  Bool_t isSelLc[2] = {true, true}; 
  pdg3[0] = 2212;
  pdg3[1] = 321;
  pdg3[2] = 211;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitLc)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    isSelLc[0] = false;
  pdg3[0] = 211;
  pdg3[1] = 321;
  pdg3[2] = 2212;
  minv2 = rd4massCalc3->InvMass2(nprongs, pdg3);
  if ((retval & (1 << kbitLc)) && (minv2 < lolim * lolim || minv2 > hilim * hilim))
    isSelLc[1] = false;
  if(!isSelLc[0] && !isSelLc[1])
   retval &= ~(1 << kbitLc);

  return retval;
}
//______________________________________________________________________________
AliAODRecoDecayHF2Prong* AliAnalysisTaskHFSimpleVertices::Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)twoTrackArray->UncheckedAt(1);

  Double_t px[2], py[2], pz[2], d0[2], d0err[2];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];

  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);

  AliAODRecoDecayHF2Prong* the2Prong = new AliAODRecoDecayHF2Prong(0x0, px, py, pz, d0, d0err, dcap1n1);
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the2Prong->SetOwnSecondaryVtx(ownsecv);
  return the2Prong;
}
//______________________________________________________________________________
AliAODRecoDecayHF3Prong* AliAnalysisTaskHFSimpleVertices::Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t bzkG)
{

  AliESDtrack* track_0 = (AliESDtrack*)threeTrackArray->UncheckedAt(0);
  AliESDtrack* track_1 = (AliESDtrack*)threeTrackArray->UncheckedAt(1);
  AliESDtrack* track_2 = (AliESDtrack*)threeTrackArray->UncheckedAt(2);

  Double_t px[3], py[3], pz[3], d0[3], d0err[3];
  Double_t momentum[3];
  GetTrackMomentumAtSecVert(track_0, secVert, momentum, bzkG);
  px[0] = momentum[0];
  py[0] = momentum[1];
  pz[0] = momentum[2];
  GetTrackMomentumAtSecVert(track_1, secVert, momentum, bzkG);
  px[1] = momentum[0];
  py[1] = momentum[1];
  pz[1] = momentum[2];
  GetTrackMomentumAtSecVert(track_2, secVert, momentum, bzkG);
  px[2] = momentum[0];
  py[2] = momentum[1];
  pz[2] = momentum[2];
  Float_t d0z0f[2], covd0z0f[3];
  track_0->GetImpactParameters(d0z0f, covd0z0f);
  d0[0] = d0z0f[0];
  d0err[0] = TMath::Sqrt(covd0z0f[0]);
  track_1->GetImpactParameters(d0z0f, covd0z0f);
  d0[1] = d0z0f[0];
  d0err[1] = TMath::Sqrt(covd0z0f[0]);
  track_2->GetImpactParameters(d0z0f, covd0z0f);
  d0[2] = d0z0f[0];
  d0err[2] = TMath::Sqrt(covd0z0f[0]);

  Double_t xdummy, ydummy;
  float dcap1n1 = track_0->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap2n1 = track_2->GetDCA(track_1, bzkG, xdummy, ydummy);
  float dcap1p2 = track_0->GetDCA(track_2, bzkG, xdummy, ydummy);
  Double_t dca[3] = {dcap1n1, dcap2n1, dcap1p2};
  Double_t dispersion = 0;
  Double_t dist12 = 0.;
  Double_t dist23 = 0.;
  Short_t charge = (Short_t)(track_0->Charge() + track_1->Charge() + track_2->Charge());

  // construct the candidate passing a NULL pointer for the secondary vertex to avoid creation of TRef
  AliAODRecoDecayHF3Prong* the3Prong = new AliAODRecoDecayHF3Prong(0x0, px, py, pz, d0, d0err, dca, dispersion, dist12, dist23, charge);
  // add a pointer to the secondary vertex via SetOwnSecondaryVtx (no TRef created)
  AliAODVertex* ownsecv=secVert->CloneWithoutRefs();
  the3Prong->SetOwnSecondaryVtx(ownsecv);
  return the3Prong;
}
//______________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskHFSimpleVertices::MakeCascade(TObjArray *twoTrackArray, AliAODVertex* secVert, Double_t bzkG){
  AliAODRecoCascadeHF *theCascade = (AliAODRecoCascadeHF*)Make2Prong(twoTrackArray,secVert,bzkG);
  if(!theCascade) return 0x0;
  AliESDtrack *trackBachelor = (AliESDtrack*)twoTrackArray->UncheckedAt(0);
  theCascade->SetCharge(trackBachelor->Charge());
  theCascade->GetSecondaryVtx()->AddDaughter(trackBachelor);
  return theCascade;
}
//______________________________________________________________________________

Int_t AliAnalysisTaskHFSimpleVertices::LcSelectionCuts(
    AliAODRecoDecayHF3Prong *cand) {
  bool isLcpKpi = true;
  bool isLcpiKp = true;
  Double_t ptCand = cand->Pt();
  if (ptCand < 0. || ptCand >= 36.)
    return 0;
  Int_t jPtBin = GetPtBin(ptCand, fPtBinLimsLc, fNPtBinsLc);
  if (jPtBin == -1)
    return 0;

  if (cand->CosPointingAngle() <= fLcCuts[jPtBin][7])
    return 0;
  if (cand->DecayLength() <= fLcCuts[jPtBin][6])
    return 0;

  if (TMath::Abs(cand->InvMassLcpKpi() - fMassLambdaC) > fLcCuts[jPtBin][0])
    isLcpKpi = false;
  if (TMath::Abs(cand->InvMassLcpiKp() - fMassLambdaC) > fLcCuts[jPtBin][0])
    isLcpiKp = false;
  if (!isLcpKpi && !isLcpiKp)
    return 0;

  if (cand->Pt2Prong(0) < fLcCuts[jPtBin][1] * fLcCuts[jPtBin][1] ||
      cand->Pt2Prong(1) < fLcCuts[jPtBin][2] * fLcCuts[jPtBin][2] ||
      cand->Pt2Prong(2) < fLcCuts[jPtBin][3] * fLcCuts[jPtBin][3])
    isLcpKpi = false;
  if (cand->Pt2Prong(0) < fLcCuts[jPtBin][3] * fLcCuts[jPtBin][3] ||
      cand->Pt2Prong(1) < fLcCuts[jPtBin][2] * fLcCuts[jPtBin][2] ||
      cand->Pt2Prong(2) < fLcCuts[jPtBin][1] * fLcCuts[jPtBin][1])
    isLcpiKp = false;
  if (!isLcpKpi && !isLcpiKp)
    return 0;

  Int_t returnValue = 0;
  if (isLcpKpi)
    returnValue += 1;
  if (isLcpiKp)
    returnValue += 2;
  return returnValue;
}
//----------------------------------------------------------------------------
Int_t AliAnalysisTaskHFSimpleVertices::MatchToMC(AliAODRecoDecay* rd, Int_t pdgabs, AliMCEvent* mcEvent,
                                                 Int_t ndgCk, const TObjArray *trkArray, const Int_t *pdgDg) const {

  Int_t ndg=rd->GetOwnSecondaryVtx()->GetNContributors();
  if(!ndg) {
    AliError("No daughters available");
    return -1;
  }
  if(ndg>10) {
    AliError("Only decays with <10 daughters supported");
    return -1;
  }
  if(ndgCk>0 && ndgCk!=ndg) {
    AliError(Form("Wrong number of daughter PDGs passed: %d requested , AliAODRecoDecay has %d",ndgCk,ndg));
    return -1;
  }
  Int_t dgLabels[10] = {0};

  // loop on daughters and write the labels
  for(Int_t i=0; i<ndg; i++) {
    AliESDtrack* trk = (AliESDtrack*)trkArray->UncheckedAt(i);
    dgLabels[i] = trk->GetLabel();
  }

  Int_t labMom[10]={0,0,0,0,0,0,0,0,0,0};
  Int_t i,j,lab,labMother,pdgMother,pdgPart;
  AliMCParticle *part=0;
  AliMCParticle *mother=0;
  Double_t pxSumDgs=0.,pySumDgs=0.,pzSumDgs=0.;
  Bool_t pdgUsed[10]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
  // loop on daughter labels
  for(i=0; i<ndg; i++) {
    labMom[i]=-1;
    lab = TMath::Abs(dgLabels[i]);
    if(lab<0) {
      printf("daughter with negative label %d\n",lab);
      return -1;
    }
    part = (AliMCParticle*)mcEvent->GetTrack(lab);
    if(!part) {
      printf("no MC particle\n");
      return -1;
    }
    // check the PDG of the daughter, if requested
    if(ndgCk>0) {
      pdgPart=TMath::Abs(part->PdgCode());
      for(j=0; j<ndg; j++) {
        if(!pdgUsed[j] && pdgPart==pdgDg[j]) {
          pdgUsed[j]=kTRUE;
          break;
        }
      }
    }

    mother = part;
    while(mother->GetMother()>=0) {
      labMother=mother->GetMother();
      mother = (AliMCParticle*)mcEvent->GetTrack(labMother);
      if(!mother) {
        printf("no MC mother particle\n");
        break;
      }
      pdgMother = TMath::Abs(mother->PdgCode());
      if(pdgMother==pdgabs) {
        labMom[i]=labMother;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs += part->Px();
        pySumDgs += part->Py();
        pzSumDgs += part->Pz();
        break;
      } else if(pdgMother>pdgabs || pdgMother<10) {
        break;
      }
    }
    if(labMom[i]==-1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters

  // check if the candidate is signal
  labMother=labMom[0];
  // all labels have to be the same and !=-1
  for(i=0; i<ndg; i++) {
    if(labMom[i]==-1)        return -1;
    if(labMom[i]!=labMother) return -1;
  }

  // check that all daughter PDGs are matched
  if(ndgCk>0) {
    for(i=0; i<ndg; i++) {
      if(pdgUsed[i]==kFALSE) return -1;
    }
  }
  mother = (AliMCParticle*)mcEvent->GetTrack(labMother);
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();
  // within 0.1%
  if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13)) > 0.00001 &&
     (TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13)) > 0.00001 &&
     (TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13)) > 0.00001)
    return -1;

  return labMother;
}

//______________________________________________________________________________
std::string AliAnalysisTaskHFSimpleVertices::GetJsonString(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  char* value=0x0;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      value=strtok(line, ":");
      value=strtok(NULL, ":");
      break;
    }
  }
  fclose(fj);
  TString sValue = value;
  sValue.ReplaceAll("\"", "");
  sValue.ReplaceAll("\n", "");
  sValue.ReplaceAll("\t", "");
  sValue.ReplaceAll(" ", "");
  return std::string(sValue.Data());
}
//______________________________________________________________________________
int AliAnalysisTaskHFSimpleVertices::GetJsonBool(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  int value=-1;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      if(temp.Contains("true")) value=1;
      if(temp.Contains("false")) value=0;
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
int AliAnalysisTaskHFSimpleVertices::GetJsonInteger(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  int value=-999;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atoi();
      break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________
float AliAnalysisTaskHFSimpleVertices::GetJsonFloat(const char* jsonFileName, const char* key){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float value=-999.;
  while(!feof(fj)){
    char* rc=fgets(line,500,fj);
    if(rc && strstr(line,key)){
      char* token=strtok(line, ":");
      token=strtok(NULL, ":");
      TString temp=token;
      temp.ReplaceAll("\"","");
      temp.ReplaceAll(",","");
      value=temp.Atof();
      break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
float* AliAnalysisTaskHFSimpleVertices::GetJsonArray(const char* jsonFileName, const char* key, int &size){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float* arrVals=0x0;
  while(!feof(fj)){
    fgets(line,500,fj);
    if(strstr(line,key)){
      TString full="";
      while(!feof(fj)){
        fgets(line,500,fj);
        int len = strlen(line);
        if(line[len-1]=='\n') line[len-1]=0;
        full.Append(line);
        if(strstr(line,"}")) break;
      }
      full.ReplaceAll("\"values\":","");
      full.ReplaceAll(" ","");
      full.ReplaceAll("},","");
      full.ReplaceAll("}","");
      TObjArray* arrStr=full.Tokenize(",");
      size=arrStr->GetEntriesFast();
      arrVals=new float[size];
      for(int j=0; j<size; j++){
        TObjString* sss=(TObjString*)arrStr->At(j);
        TString strval=sss->GetString();
        strval.ReplaceAll("[","");
        strval.ReplaceAll("]","");
        strval.ReplaceAll("\"","");
        arrVals[j]=strval.Atof();
      }
      arrStr->Delete();
      delete arrStr;
    }
  }
  return arrVals;
}

//______________________________________________________________________________
float** AliAnalysisTaskHFSimpleVertices::GetJsonMatrix(const char* jsonFileName, const char* key, int &size1, int &size2){
  FILE* fj=fopen(jsonFileName,"r");
  char line[500];
  float** arrVals=0x0;
  while(!feof(fj)){
    fgets(line,500,fj);
    if(strstr(line,key)){
      TString full="";
      while(!feof(fj)){
        fgets(line,500,fj);
        int len = strlen(line);
        if(line[len-1]=='\n') line[len-1]=0;
        full.Append(line);
        if(strstr(line,"}")) break;
      }
      full.ReplaceAll("\"values\":","");
      full.ReplaceAll(" ","");
      full.ReplaceAll("},","");
      full.ReplaceAll("}","");
      TObjArray* rowArrStr=full.Tokenize("]");
      size1=rowArrStr->GetEntriesFast();
      arrVals=new float*[size1];
      for(int j=0; j<size1; j++){
        TObjString* rowStr=(TObjString*)rowArrStr->At(j);
        TString rowStrVal=rowStr->GetString();
        TObjArray* arrStr=rowStrVal.Tokenize(",");
        size2=arrStr->GetEntriesFast();
        arrVals[j]=new float[size2];
        for(int k=0; k<size2; k++){
          TObjString* sss=(TObjString*)arrStr->At(k);
          TString strval=sss->GetString();
          strval.ReplaceAll(",[","");
          strval.ReplaceAll("[","");
          strval.ReplaceAll("]","");
          strval.ReplaceAll("\"","");
          arrVals[j][k]=strval.Atof();
        }
      }
    }
  }
  return arrVals;
}

//______________________________________________________________________________
void AliAnalysisTaskHFSimpleVertices::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  fHistNEvents= dynamic_cast<TH1F*>(fOutput->FindObject("hNEvents"));
  printf("AliAnalysisTaskHFSimpleVertices::Terminate --- Number of events: read = %.0f\n",fHistNEvents->GetBinContent(1));
  return;
}
