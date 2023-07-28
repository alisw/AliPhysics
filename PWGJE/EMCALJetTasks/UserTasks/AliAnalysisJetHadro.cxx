/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, proviyaded that the above copyright notice appears in all *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purapose. It is         *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//                        Analysis for Jet Hadrochemistry                       //
//                                                                              //
//    This analysis extracts pT-spectra of charged kaons, protons, and pions    //
//                      for the inclusive event and in jets.                    //
//   It is based on particles identification via the dE/dx signal of the TPC.   //
//                                                                              //
// Author: Sierra Weyhmiller <sierra.lisa.weyhmiller@cern.ch>, Yale University  //
//      Author: Mesut Arslandok <mesut.arslandok@cern.ch>, Yale University      //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPDG.h"
#include "AliMathBase.h"
#include "AliESDFMD.h"
#include "AliFMDFloatMap.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliTPCdEdxInfo.h"
#include "AliKFVertex.h"
#include "AliLumiTools.h"
#include "AliKFParticle.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliVParticle.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliTPCParam.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliKFParticle.h"
//#include "AliAnalysisTaskFilteredTree.h"
#include "AliAnalysisJetHadro.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliRunLoader.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
//#include "AliESDtools.h"
#include "AliFJWrapper.h"
#include "AliEmcalJet.h"
#include "AliEmcalJetTask.h"
#include "AliRhoParameter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
#include <TObjString.h>
using namespace std;
using std::cout;
using std::setw;

ClassImp(AliAnalysisJetHadro)

const char* AliAnalysisJetHadro::fEventInfo_centEstStr[] = {"V0M","CL0","CL1"};

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro()
: AliAnalysisTaskEmcalJet("TaskEbyeRatios"), fEventCuts(0), fPIDResponse(0),fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
//fESDtool(nullptr),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsEMCBGconst(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetResonance(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fFilldscaledTree(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fFillFastJet(kTRUE),
fFillEMCJet(kTRUE),
fFillIncTracks(kTRUE),
fFillpTPC(kTRUE),
fFillp(kTRUE),
fFillpT(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBGConst(kTRUE),
fcent_min(0.0),
fcent_max(100.0),
fjetMinPtSub(-1000.0),
fjetMinArea(-1000.0),
fRunFastSimulation(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fIncludeTOF(kFALSE),
fUseCouts(kFALSE),
fRunNumberForExpecteds(0),
fFillExpecteds(kFALSE),
fNSettings(17),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fMomExpec_NBins(2000),
fMomExpec_Low(0.0),
fMomExpec_High(20.0),
fEtaExpec_NBins(9),
fEtaExpec_Low(0.0),
fEtaExpec_High(0.9),
fNdEdxBins(2000),
fDEdxUp(0.0),
fDEdxDown(20.0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fTPCFindable(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(14),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fIsMCPileup(0),
fJetContainer(0),
fbgJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealJets(0),
ftotalJetArea(0),
ftotalNumRealJets(0),
ftotalNumRealJetEvents(0),
ftotalNumIncEvents(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fSystCentEstimatetor(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_dEdx_p(0),
fHistIncTracks_dEdx_pT(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_kin(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_dEdx_p(0),
fHistJetTracks_dEdx_pT(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_kin(0),
fHistIncTracks_mpi(0),
/*
fHistIncTracks_mpi_small(0),
fHistIncTracks_spi_small(0),
fHistIncTracks_mel_small(0),
fHistIncTracks_sel_small(0),
fHistIncTracks_mka_small(0),
fHistIncTracks_ska_small(0),
fHistIncTracks_mpr_small(0),
fHistIncTracks_spr_small(0),
*/
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0),
fEffMatrixNSigmasTOF(0)
{
  // default Constructor
  /* fast compilation test
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  .L AliAnalysisJetHadro.cxx++
  */
}

//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro(const char *name)
: AliAnalysisTaskEmcalJet(name), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
//fESDtool(nullptr),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreejetsEMCBGconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetResonance(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fFilldscaledTree(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fFillFastJet(kTRUE),
fFillEMCJet(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBGConst(kTRUE),
fFillIncTracks(kTRUE),
fFillpTPC(kTRUE),
fFillp(kTRUE),
fFillpT(kTRUE),
fcent_min(0.0),
fcent_max(0.0),
fjetMinPtSub(-1000.0),
fjetMinArea(-1000.0),
fRunFastSimulation(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fIncludeTOF(kFALSE),
fUseCouts(kFALSE),
fRunNumberForExpecteds(0),
fFillExpecteds(kFALSE),
fNSettings(17),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fMomExpec_NBins(2000),
fMomExpec_Low(0.0),
fMomExpec_High(20.0),
fEtaExpec_NBins(9),
fEtaExpec_Low(0.0),
fEtaExpec_High(0.9),
fNdEdxBins(2000),
fDEdxUp(0.0),
fDEdxDown(20.0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fTPCFindable(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(14),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fJetContainer(0),
fbgJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealJets(0),
ftotalJetArea(0),
ftotalNumRealJets(0),
ftotalNumRealJetEvents(0),
ftotalNumIncEvents(0),
fIsMCPileup(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fSystCentEstimatetor(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_dEdx_p(0),
fHistIncTracks_dEdx_pT(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_kin(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_dEdx_p(0),
fHistJetTracks_dEdx_pT(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_kin(0),
fHistIncTracks_mpi(0),
/*
fHistIncTracks_mpi_small(0),
fHistIncTracks_spi_small(0),
fHistIncTracks_mel_small(0),
fHistIncTracks_sel_small(0),
fHistIncTracks_mka_small(0),
fHistIncTracks_ska_small(0),
fHistIncTracks_mpr_small(0),
fHistIncTracks_spr_small(0),
*/
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0),
fEffMatrixNSigmasTOF(0)
{
  //
  //         standard constructur which should be used
  //
  /* fast compilation test
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  .L AliAnalysisJetHadro.cxx++
  */
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:***************** CONSTRUCTOR CALLED: AliAnalysisJetHadro  *****************"<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  // ==========================================
  //
  // ==========================================
  //
  // Define outputs
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());
  DefineOutput(12, TTree::Class());
  DefineOutput(13, TTree::Class());
  // ==========================================

  ftotalJetArea = 0.0;
  ftotalNumRealJets = 0;
  ftotalNumRealJetEvents  = 0;
  ftotalNumIncEvents  = 0;

}
//________________________________________________________________________
AliAnalysisJetHadro::~AliAnalysisJetHadro()
{

  //
  // Destructor
  //
  std::cout << " Info::siweyhmi: ===== In the Destructor ===== " << std::endl;
  if (fHistEmptyEvent)      delete fHistEmptyEvent;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistCentralityImpPar)delete fHistCentralityImpPar;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHistIncTracks_dEdx)          delete fHistIncTracks_dEdx;
  if (fHistIncTracks_dEdx_p)          delete fHistIncTracks_dEdx_p;
  if (fHistIncTracks_dEdx_pT)          delete fHistIncTracks_dEdx_pT;
  if (fHistIncTracks_moms)          delete fHistIncTracks_moms;
  if (fHistIncTracks_moms_p)          delete fHistIncTracks_moms_p;
  if (fHistIncTracks_kin)          delete fHistIncTracks_kin;
  if (fHistJetTracks_dEdx)          delete fHistJetTracks_dEdx;
  if (fHistJetTracks_dEdx_p)          delete fHistJetTracks_dEdx_p;
  if (fHistJetTracks_dEdx_pT)          delete fHistJetTracks_dEdx_pT;
  if (fHistJetTracks_moms)          delete fHistJetTracks_moms;
  if (fHistJetTracks_moms_p)          delete fHistJetTracks_moms_p;
  if (fHistJetTracks_kin)          delete fHistJetTracks_kin;
  /*
  if (fHistIncTracks_mpi_small)          delete fHistIncTracks_mpi_small;
  if (fHistIncTracks_spi_small)          delete fHistIncTracks_spi_small;
  if (fHistIncTracks_mel_small)          delete fHistIncTracks_mel_small;
  if (fHistIncTracks_sel_small)          delete fHistIncTracks_sel_small;
  if (fHistIncTracks_mka_small)          delete fHistIncTracks_mka_small;
  if (fHistIncTracks_ska_small)          delete fHistIncTracks_ska_small;
  if (fHistIncTracks_mpr_small)          delete fHistIncTracks_mpr_small;
  if (fHistIncTracks_spr_small)          delete fHistIncTracks_spr_small;
  */
  if (fHistIncTracks_mpi)          delete fHistIncTracks_mpi;
  if (fHistIncTracks_spi)          delete fHistIncTracks_spi;
  if (fHistIncTracks_mel)          delete fHistIncTracks_mel;
  if (fHistIncTracks_sel)          delete fHistIncTracks_sel;
  if (fHistIncTracks_mka)          delete fHistIncTracks_mka;
  if (fHistIncTracks_ska)          delete fHistIncTracks_ska;
  if (fHistIncTracks_mpr)          delete fHistIncTracks_mpr;
  if (fHistIncTracks_spr)          delete fHistIncTracks_spr;
  if (fHistJet_ptsub_v_area)          delete fHistJet_ptsub_v_area;
  if (fHistJet_kin)          delete fHistJet_kin;
  if (fHistJet_moms)          delete fHistJet_moms;
  if (fHistCent)            delete fHistCent;
  if (fHistPhi)             delete fHistPhi;
  //
  if (fPIDCombined) delete fPIDCombined;
  if (fESDtrackCuts)          delete fESDtrackCuts;
  if (fESDtrackCuts_Bit96)    delete fESDtrackCuts_Bit96;
  if (fESDtrackCuts_Bit128)   delete fESDtrackCuts_Bit128;
  if (fESDtrackCuts_Bit768)   delete fESDtrackCuts_Bit768;
  if (fESDtrackCutsLoose)     delete fESDtrackCutsLoose;
  if (fTreeSRedirector)       delete fTreeSRedirector;
  if (fESDtrackCutsCleanSamp) delete fESDtrackCutsCleanSamp;
  if (fPileUpTightnessCut4) delete fPileUpTightnessCut4;
  if (fPileUpTightnessCut3) delete fPileUpTightnessCut3;
  if (fPileUpTightnessCut2) delete fPileUpTightnessCut2;
  if (fPileUpTightnessCut1) delete fPileUpTightnessCut1;

}
//
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
//
void AliAnalysisJetHadro::Initialize()
{
  //
  // updating parameters in case of changes (standard cuts and the eta window)
  //
  std::cout << " Info::siweyhmi: ===== In the Initialize ===== " << std::endl;
  if (fRunFastSimulation)    { std::cout << " Info::siweyhmi: !!! We are running fast simulation return !!! " << std::endl; return; }
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  //
  // ------------------------------------------------
  //
  // tight DCA cut used by Emil
  fESDtrackCuts_Bit96     = new AliESDtrackCuts("fESDtrackCuts_Bit96","");
  fESDtrackCuts_Bit96->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit96->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit96->SetMinNClustersTPC(70); // ???? should be --> fESDtrackCuts_Bit96->SetMinNCrossedRowsTPC(70);
  fESDtrackCuts_Bit96->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts_Bit96->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit96->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts_Bit96->SetRequireITSRefit(kTRUE);
  fESDtrackCuts_Bit96->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  fESDtrackCuts_Bit96->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexZ(2);
  fESDtrackCuts_Bit96->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts_Bit96->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts_Bit96->SetMaxChi2PerClusterITS(36);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(4);
  }
  //
  // TPC only tracks
  fESDtrackCuts_Bit128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCuts_Bit128->SetName("Bit128");
  fESDtrackCuts_Bit128->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit128->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit128->SetMinNClustersTPC(70);
  fESDtrackCuts_Bit128->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit128->SetDCAToVertex2D(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(4);
  }
  //
  // Hybrid tracks
  fESDtrackCuts_Bit768 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts_Bit768->SetName("Bit768");
  fESDtrackCuts_Bit768->SetEtaRange(-0.9,0.9);
  fESDtrackCuts_Bit768->SetPtRange(0.15,1000.);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit768->SetDCAToVertex2D(kTRUE);
  fESDtrackCuts_Bit768->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit768->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts_Bit768->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  fESDtrackCuts_Bit768->SetRequireITSRefit(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(4);
  }
  //
  // for the systematic check fill all tracks and tag them with cutbit but for MC do not
  //
  fESDtrackCuts = new AliESDtrackCuts("esdTrackCuts","");
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  fESDtrackCuts->SetPtRange(0.15,1000.);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(80);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");
  fESDtrackCuts->SetMaxDCAToVertexXY(2.4);   // hybrid cuts  TODO
  fESDtrackCuts->SetMaxDCAToVertexZ(3.2);    // hybrid cuts  TODO
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetDCAToVertex2D(kTRUE);  // fESDtrackCuts->SetDCAToVertex2D(kFALSE);    TODO
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if (fIncludeITS) {
    // require ITS pixels  -->  Reason for the empty events and structure in phi
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }
  //
  // very loose cuts --> cuts will be tightened using the bipmap
  fESDtrackCutsLoose = new AliESDtrackCuts("esdTrackCutsLoose","");
  fESDtrackCutsLoose->SetEtaRange(-100.,100.);
  fESDtrackCutsLoose->SetPtRange(0.1,100000.);
  fESDtrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsLoose->SetMinNClustersTPC(50);
  fESDtrackCutsLoose->SetMinNCrossedRowsTPC(50);
  fESDtrackCutsLoose->SetMaxDCAToVertexXY(10);   // hybrid cuts  TODO
  fESDtrackCutsLoose->SetMaxDCAToVertexZ(10);    // hybrid cuts  TODO
  //
  // track cuts to be used for v0s
  fESDtrackCutsCleanSamp = new AliESDtrackCuts("AliESDtrackCutsV0","");
  fESDtrackCutsCleanSamp -> SetEtaRange(-1.5,1.5);
  fESDtrackCutsCleanSamp -> SetPtRange(0.,1e10);
  fESDtrackCutsCleanSamp -> SetMinNCrossedRowsTPC(80);
  fESDtrackCutsCleanSamp -> SetRequireTPCRefit(kTRUE);
  fESDtrackCutsCleanSamp -> SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsCleanSamp -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsCleanSamp -> SetMaxChi2PerClusterITS(36);
  fESDtrackCutsCleanSamp -> SetMaxFractionSharedTPCClusters(0.4);
  // ------------------------------------------------

  fPileUpTightnessCut4 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut3 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut2 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut1 = new AliEventCuts(kFALSE);

  fPileUpTightnessCut4->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,4);
  fPileUpTightnessCut3->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,3);
  fPileUpTightnessCut2->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,2);
  fPileUpTightnessCut1->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
  //
  //
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
  std::cout << " Info::siweyhmi: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisJetHadro::UserCreateOutputObjects()
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  if (!fDEdxCheck) Initialize();
  std::cout << " Info::siweyhmi: ===== In the UserCreateOutputObjects ===== " << std::endl;
  // ------------  setup PIDCombined  ---------------
  fPIDCombined = new AliPIDCombined("pidCombined","");
  //
  // **********************   Input handler to get the PID object *********************
  if (!(fRunFastSimulation)) {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
    AliFatal("Input handler needed");
    else {
      fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
      if (!fPIDResponse) std::cout << " Info::siweyhmi: ======= PIDResponse object was not created ====== " << std::endl;
    }
  }
  //
  // ************************************************************************
  //   OpenFile output --> one can open several files
  // ************************************************************************
  //
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  //
  if (fDefaultEventCuts) fEventCuts.AddQAplotsToList(fListHist); /// fList is your output TList
  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  //fHistEmptyEvent        = new TH1F("hEmptyEvent",           "control histogram to count empty events"    , 10,  0., 10.);
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 10,  0., 100.);
  //fHistCentralityImpPar  = new TH1F("hCentralityImpPar",     "control histogram for centrality imppar"    , 10,  0., 100.);
  //fHistImpParam          = new TH1F("hImpParam",             "control histogram for impact parameter"     , 200, 0., 20.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);

  if (fFillpTPC) fHistIncTracks_dEdx    = new TH3F("fHistIncTracks_dEdx",   "dEdx histogram for inclusive tracks"        , fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  if (fFillp) fHistIncTracks_dEdx_p    = new TH3F("fHistIncTracks_dEdx_p",   "dEdx histogram for inclusive tracks"        ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  if (fFillpT) fHistIncTracks_dEdx_pT    = new TH3F("fHistIncTracks_dEdx_pT",   "dEdx histogram for inclusive tracks"        ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  if (fFillpTPC) fHistIncTracks_moms    = new TH2F("fHistIncTracks_moms",   "All mom types for inclusive tracks"         ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High);
  if (fFillp) fHistIncTracks_moms_p    = new TH2F("fHistIncTracks_moms_p",   "All mom types for inclusive tracks"         ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High);

  fHistIncTracks_kin    = new TH3F("fHistIncTracks_kin",     "Kinematics histogram for inclusive tracks"  , 200, 0., 20., 9, 0.0, 0.9, 64, -3.2, 3.2);

  if (fFillpTPC) fHistJetTracks_dEdx    = new TH3F("fHistJetTracks_dEdx",   "dEdx histogram for Jet tracks"        ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  if (fFillp) fHistJetTracks_dEdx_p    = new TH3F("fHistJetTracks_dEdx_p",   "dEdx histogram for Jet tracks"        ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  if (fFillpT) fHistJetTracks_dEdx_pT    = new TH3F("fHistJetTracks_dEdx_pT",   "dEdx histogram for Jet tracks"        ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,   fNdEdxBins,fDEdxUp,fDEdxDown, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  if (fFillpTPC) fHistJetTracks_moms    = new TH2F("fHistJetTracks_moms",   "All mom types for Jet tracks"         ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High);
  if (fFillp) fHistJetTracks_moms_p    = new TH2F("fHistJetTracks_moms_p",   "All mom types for Jet tracks"         ,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High,fMomExpec_NBins,fMomExpec_Low, fMomExpec_High);

  fHistJetTracks_kin    = new TH3F("fHistJetTracks_kin",     "Kinematics histogram for Jet tracks"  , 200, 0., 20., 9, 0.0, 0.9, 64, -3.2, 3.2);

  /*
  fHistIncTracks_mpi_small  = new TH2D("fHistIncTracks_mpi_small",     "Expected mean pion histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_spi_small  = new TH2D("fHistIncTracks_spi_small",     "Expected sigma pion histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mel_small  = new TH2D("fHistIncTracks_mel_small",     "Expected mean electron histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_sel_small  = new TH2D("fHistIncTracks_sel_small",     "Expected sigma electron histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mka_small  = new TH2D("fHistIncTracks_mka_small",     "Expected mean kaon histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_ska_small  = new TH2D("fHistIncTracks_ska_small",     "Expected sigma kaon histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mpr_small  = new TH2D("fHistIncTracks_mpr_small",     "Expected mean proton histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_spr_small  = new TH2D("fHistIncTracks_spr_small",     "Expected sigma proton histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  */

  fHistIncTracks_mpi  = new TH3F("fHistIncTracks_mpi",     "Expected mean pion histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 1000, 0., 1000., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_spi  = new TH3F("fHistIncTracks_spi",     "Expected sigma pion histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 500, 0., 50., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mel  = new TH3F("fHistIncTracks_mel",     "Expected mean electron histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 1000, 0., 1000., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_sel  = new TH3F("fHistIncTracks_sel",     "Expected sigma electron histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 500, 0., 50., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mka  = new TH3F("fHistIncTracks_mka",     "Expected mean kaon histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 1000, 0., 1000., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_ska  = new TH3F("fHistIncTracks_ska",     "Expected sigma kaon histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 500, 0., 50., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistIncTracks_mpr  = new TH3F("fHistIncTracks_mpr",     "Expected mean proton histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 1000, 0., 1000., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);
  fHistIncTracks_spr  = new TH3F("fHistIncTracks_spr",     "Expected sigma proton histogram for inclusive tracks"  , fMomExpec_NBins, fMomExpec_Low, fMomExpec_High, 500, 0., 50., fEtaExpec_NBins, fEtaExpec_Low, fEtaExpec_High);

  fHistJet_ptsub_v_area  = new TH2F("fHistJet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);
  fHistJet_kin  = new TH3F("fHistJet_kin", "Kinematics histogram for Jets"  , 300, 0., 300., 48, -0.6, 0.6, 130, 0., 6.5);
  fHistJet_moms  = new TH2F("fHistJet_moms", "All mom types for jets"  , 300, 0., 300., 600, 0., 600.);

  //fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  if (fMCtrue) {
    //fListHist->Add(fHistCentralityImpPar);
    //fListHist->Add(fHistImpParam);
  }
  fListHist->Add(fHistVertex);

  if (fFillpTPC) fListHist->Add(fHistIncTracks_dEdx);
  if (fFillp) fListHist->Add(fHistIncTracks_dEdx_p);
  if (fFillpT) fListHist->Add(fHistIncTracks_dEdx_pT);
  if (fFillpTPC) fListHist->Add(fHistIncTracks_moms);
  if (fFillp) fListHist->Add(fHistIncTracks_moms_p);
  fListHist->Add(fHistIncTracks_kin);

  if (fFillpTPC) fListHist->Add(fHistJetTracks_dEdx);
  if (fFillp) fListHist->Add(fHistJetTracks_dEdx_p);
  if (fFillpT) fListHist->Add(fHistJetTracks_dEdx_pT);
  if (fFillpTPC) fListHist->Add(fHistJetTracks_moms);
  if (fFillp) fListHist->Add(fHistJetTracks_moms_p);
  fListHist->Add(fHistJetTracks_kin);

/*
  fListHist->Add(fHistIncTracks_mpi_small);
  fListHist->Add(fHistIncTracks_spi_small);

  fListHist->Add(fHistIncTracks_mel_small);
  fListHist->Add(fHistIncTracks_sel_small);

  fListHist->Add(fHistIncTracks_mka_small);
  fListHist->Add(fHistIncTracks_ska_small);

  fListHist->Add(fHistIncTracks_mpr_small);
  fListHist->Add(fHistIncTracks_spr_small);
*/

  fListHist->Add(fHistJet_ptsub_v_area);
  fListHist->Add(fHistJet_kin);
  fListHist->Add(fHistJet_moms);

  fListHist->Add(fHistIncTracks_mpi);
  fListHist->Add(fHistIncTracks_spi);

  fListHist->Add(fHistIncTracks_mel);
  fListHist->Add(fHistIncTracks_sel);

  fListHist->Add(fHistIncTracks_mka);
  fListHist->Add(fHistIncTracks_ska);

  fListHist->Add(fHistIncTracks_mpr);
  fListHist->Add(fHistIncTracks_spr);

  fEventInfo_CentralityEstimates  = new TVectorF(3);
  for (Int_t i=0;i<3;i++) (*fEventInfo_CentralityEstimates)[i]=-10.;
  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fTreejetsEMCconst  = ((*fTreeSRedirector)<<"jetsEMCconst").GetTree();
  fTreejetsEMCBGconst  = ((*fTreeSRedirector)<<"jetsEMCBGconst").GetTree();
  fTreejetsFJ        = ((*fTreeSRedirector)<<"jetsFJ").GetTree();
  fTreejetsFJBG      = ((*fTreeSRedirector)<<"jetsFJBG").GetTree();
  fTreejetsEMC       = ((*fTreeSRedirector)<<"jetsEMC").GetTree();
  fTreejetsEMCBG       = ((*fTreeSRedirector)<<"jetsEMCBG").GetTree();
  fTreejetsFJconst   = ((*fTreeSRedirector)<<"jetsFJconst").GetTree();
  fTreejetsFJBGconst   = ((*fTreeSRedirector)<<"jetsFJBGconst").GetTree();
  fTreejetResonance  = ((*fTreeSRedirector)<<"jetResonance").GetTree();
  fTreejetEvents     = ((*fTreeSRedirector)<<"jeteventInfo").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"tracks").GetTree();
  //
  // ************************************************************************
  //   Send output objects to container
  // ************************************************************************
  //
  PostData(1, fListHist);
  PostData(2, fTreejetsEMCconst);
  PostData(3, fTreejetsEMCBGconst);
  PostData(4, fTreejetsFJ);
  PostData(5, fTreejetsFJBG);
  PostData(6, fTreejetsFJconst);
  PostData(7, fTreejetsFJBGconst);
  PostData(8, fTreejetResonance);
  PostData(9, fTreejetEvents);
  PostData(10, fTreejetsEMC);
  PostData(11, fTreejetsEMCBG);
  PostData(12, fTreeMC);
  PostData(13, fTreeCuts);

  fEventCuts.SetManualMode();

  std::cout << " Info::siweyhmi: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::Run()
{
  //
  // main event loop
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the UserExec ===== " << std::endl;
  //
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  fEventCountInFile++;
  //
  //  get the filename
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf(" Error::siweyhmi: Could not receive input chain"); return kFALSE; }
  TString tmpChunkname = fChunkName;
  TObjString fileName(chain->GetCurrentFile()->GetName());
  fChunkName = (TString)fileName.GetString();
  if (tmpChunkname != fChunkName) std::cout <<  " Info::siweyhmi: ===== Current chunk name is ===== " << fChunkName << std::endl;
  //
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  //
  if (eventHandler) fMCEvent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  TString genheaderName;
  if (fMCEvent){
    genHeader = fMCEvent->GenEventHeader();
    genheaderName = genHeader->GetName();
    if(!genHeader){ Printf(" Error::siweyhmi: Event generator header not available!!!\n"); return kFALSE; }
  }
  //
  // If ESDs exist get some event variables
  //
  fCentrality = -5;
  fCentImpBin =-10.;
  AliCentrality    *esdCentrality = 0x0;
  AliMultSelection *MultSelection = 0x0;
  AliMultSelectionTask *MultSelectionTask = 0x0;
  ULong64_t gid=0;
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (fESD)
  {
    //
    // Init magnetic filed for golden chi2 cut
    fESD->InitMagneticField();
    //
    // event selection
    fPileUpBit=0;
    if (fDefaultEventCuts && fESD){

      if ( (fPassIndex==3 || fPassIndex==2) && fYear>2013){
        //
        // pileup bit: 0bxxxx, where the leftmost bit is the tightest and the rightmost is the loosest cut
        // OOB pileup cut (for Pb-Pb) based on ITS and TPC clusters: 0-> no cut; 1-> default cut (remove all OOB pileup); 2-> looser cut; 3-> even more looser cut; 4-> very loose cut
        if (fPileUpTightnessCut4->AcceptEvent(fESD)) { fPileUpBit |= 1 << 3; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut3->AcceptEvent(fESD)) { fPileUpBit |= 1 << 2; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut2->AcceptEvent(fESD)) { fPileUpBit |= 1 << 1; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        if (fPileUpTightnessCut1->AcceptEvent(fESD)) { fPileUpBit |= 1 << 0; if (fUseCouts) std::cout << "pileupbit: " << std::bitset<4>(fPileUpBit) << std::endl;}
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1); // standard
        if (fUseCouts) std::cout << " aaaaa " << std::endl; //CHANGE
        if (!fEventCuts.AcceptEvent(fESD)) {cout<< "pileup event " << endl; return kFALSE;}
      }
    }
    //
    // calculate TPC mult
    fTPCMult = 0;
    for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack){
      AliESDtrack *track = fESD->GetTrack(itrack);
      if (track->IsOn(AliESDtrack::kTPCin)) fTPCMult++;
    }
    //
    //
    esdCentrality = fESD->GetCentrality();
    MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
    fRunNo = fESD->GetRunNumber();
    //
    // if the run number is specified, fill expecteds only for that run. otherwise fill for any run
    if (fRunNumberForExpecteds > 0 && fRunOnGrid)
    fFillExpecteds = (fRunNo == fRunNumberForExpecteds);
    else if (fRunNumberForExpecteds < 0)
    fFillExpecteds = kFALSE;
    else
    fFillExpecteds = kTRUE;
    //
    static Int_t timeStampCache = -1;
    if (!fMCtrue) {
      fTimeStamp = fESD->GetTimeStampCTPBCCorr();
      const char *ocdb;
      if(!fRunOnGrid) ocdb = Form("local:///cvmfs/alice.cern.ch/calibration/data/%d/OCDB",fYear);
      else ocdb = "raw://";
      //
      // retrieve interaction rate
      if (timeStampCache!=Int_t(fTimeStamp)) timeStampCache=Int_t(fTimeStamp);
      if (!gGrid && timeStampCache>0) {
        AliInfo("Trying to connect to AliEn ...");
        TGrid::Connect("alien://");
      } else {
        fEventInfo_LumiGraph = (TGraph*)AliLumiTools::GetLumiFromCTP(fRunNo,ocdb);
        fIntRate   = fEventInfo_LumiGraph->Eval(fTimeStamp); delete fEventInfo_LumiGraph;
      }
    }
    fEventMult = fESD->GetNumberOfTracks();
    fBField    = fESD->GetMagneticField();
    fEvent     = fESD->GetEventNumberInFile();
    fBeamType  = fESD->GetBeamType();
    //
    // Global id for the event --> which is made with Hashing
    //
    ULong64_t orbitID      = (ULong64_t)fESD->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)fESD->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)fESD->GetPeriodNumber();
    gid = ((periodID << 36) | (orbitID << 12) | bunchCrossID);
    fEventGID = gid;    // uniqe event id for real data
    // fTimeStamp             = fESD->GetTimeStamp();
    // fEventGID              = TMath::Abs(Int_t(TString::Hash(&gid,sizeof(Int_t))));    // uniqe event id for real data
  }

  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  if (!fPIDResponse && !(fRunFastSimulation)) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && fPIDResponse) fPIDResponse->SetCurrentMCEvent(fMCEvent);
  //
  if(fMCtrue)
  {
    //
    // Add different generator particles to PDG Data Base to avoid problems when reading MC generator particles
    AliPDG::AddParticlesToPdgDataBase();
    //
    // ========================== MC =========================
    //
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) { Printf(" Error::siweyhmi: No MC stack available !!!\n"); return kFALSE;}
    //
    if (MultSelection) {
      fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::siweyhmi: Centralitity is taken from MultSelection " << fCentrality << std::endl;
    } else if (esdCentrality) {
      fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::siweyhmi: Centralitity is taken from esdCentrality " << fCentrality << std::endl;
    }
    //
    // impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    AliGenHijingEventHeader *lHIJINGHeader = 0x0;  // event header for HIJING
    AliGenHepMCEventHeader *lHepMCHeader = 0x0;    // event header for EPOS
    //
    // If EPOS based MC
    if ( genheaderName.Contains("EPOSLHC") && fMCEvent ){
      if (genHeader->InheritsFrom(AliGenHepMCEventHeader::Class()))
      lHepMCHeader = (AliGenHepMCEventHeader*)genHeader;
      if (lHepMCHeader ){
        fNHardScatters = lHepMCHeader->Ncoll_hard(); // Number of hard scatterings
        fNProjectileParticipants = lHepMCHeader->Npart_proj(); // Number of projectile participants
        fNTargetParticipants     = lHepMCHeader->Npart_targ(); // Number of target participants
        fNNColl   = lHepMCHeader->Ncoll(); // Number of NN (nucleon-nucleon) collisions
        fNNwColl  = lHepMCHeader->N_Nwounded_collisions(); // Number of N-Nwounded collisions
        fNwNColl  = lHepMCHeader->Nwounded_N_collisions(); // Number of Nwounded-N collisons
        fNwNwColl = lHepMCHeader->Nwounded_Nwounded_collisions();// Number of Nwounded-Nwounded collisions
        fMCImpactParameter = lHepMCHeader->impact_parameter();
        if (fUseCouts)  std::cout << " Info::siweyhmi: EPOS: Centralitity is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenEposEventHeader*) genHeader)->GetName() << std::endl;
        fHistImpParam->Fill(fMCImpactParameter);
        if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
        if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
        if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
        if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
        if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
        if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
        if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
        if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
        if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
        if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
        //fHistCentralityImpPar->Fill(fCentImpBin);
      }
    } else if (!TMath::IsNaN(((AliGenHijingEventHeader*) genHeader)->ImpactParameter()) && fMCEvent){
      //
      // If HIJING based MC
      lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
      fNHardScatters = lHIJINGHeader->HardScatters();
      fNProjectileParticipants = lHIJINGHeader->ProjectileParticipants();
      fNTargetParticipants     = lHIJINGHeader->TargetParticipants();
      fNNColl   = lHIJINGHeader->NN();
      fNNwColl  = lHIJINGHeader->NNw();
      fNwNColl  = lHIJINGHeader->NwN();
      fNwNwColl = lHIJINGHeader->NwNw();
      fMCImpactParameter = lHIJINGHeader->ImpactParameter();
      fHistImpParam->Fill(fMCImpactParameter);
      if (fUseCouts)  std::cout << " Info::siweyhmi: HIJING: Centralitity is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenHijingEventHeader*) genHeader)->GetName() << std::endl;
      if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
      if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
      //fHistCentralityImpPar->Fill(fCentImpBin);
    }
    if (fCentrality<0) fCentrality=fCentImpBin;
    //
    // Use file name in Hashing to create unique event ID
    fEventGIDMC  = TMath::Abs(Int_t(TString::Hash(&fEventCountInFile,sizeof(Int_t))));    // uniqe event id for real data
    // fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile+(1000*TMath::Abs(fMCImpactParameter)));  // ????
    fEventGID    = fEventGIDMC;
    if (fUseCouts) {
      std::cout << " Info::siweyhmi: ========================================================================================== " << std::endl;
      std::cout << " Info::siweyhmi: " << fEventCountInFile << " ----- eventIDMC = " << fEventGIDMC << "   " << fChunkName << std::endl;
      std::cout << " Info::siweyhmi: Centrality = " << fCentrality << " ----- Impact Param = " << fMCImpactParameter << " fCentralityImp = " << fCentImpBin << std::endl;
      std::cout << " Info::siweyhmi: ========================================================================================== " << std::endl;
    }
    //
  }
  //
  if (!(fRunFastSimulation))
  {
    //
    // ========================== Real =========================
    //
    if (fUseCouts) std::cout << " bbbbb " << std::endl; //CHANGE
    if (!fESD)          { Printf(" Error::siweyhmi: fESD not available"); return kFALSE; }
    if (!fESDtrackCuts) { Printf(" Error::siweyhmi: fESDtrackCuts not available"); return kFALSE; }
    //
    // ------------------------------------------------
    // ------- monitor vertex position =---------------
    // ------------------------------------------------
    //
    Bool_t isVertexOk = kTRUE;
    fVertex = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexSPD = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexTPC = fESD->GetPrimaryVertexTracks();
    if( fVertex->GetNContributors()<1) isVertexOk = kFALSE;
    if( fVertex->GetNContributors()>1) {
      vertexSPD = fESD->GetPrimaryVertexSPD();    // SPD vertex
      vertexTPC = fESD->GetPrimaryVertexTPC();    // TPC vertex
      fTPCvZ = vertexTPC->GetZ();
      fSPDvZ = vertexSPD->GetZ();
      fVz    = fVertex->GetZ();
      TString vertexType = fVertex->GetTitle();    // ??? Put condition Abs(vertex-vertexTPC) as a bool_t into ttree
      //if ( vertexType.Contains("vertexer: Z") && (fVertex->GetDispersion() > 0.04 || fVertex->GetZRes() > 0.25) ) isVertexOk = kFALSE; // TODO
    }
    fMultiplicity    = fVertex->GetNContributors();    // fMultiplicity = fESD -> GetNumberOfTracks();
    fNContributors   = fVertex->GetNContributors();
    fMultiplicityMC  = fMultiplicity;
    //
    // ------------------------------------------------
    // ------- event vertex cut along Z ---------------
    // ------------------------------------------------
    //
    // if (fMCtrue && TMath::Abs(fVz) > 10) return;   // For MC put fixed cut
    if (fUseCouts) std::cout << " ccccc " << std::endl; //CHANGE
    if (TMath::Abs(fVz)>10) return kFALSE;
    //
    if (fVertex && isVertexOk) fHistVertex->Fill(fVz);
    else return kFALSE;
    //
    // ------------------------------------------------
    // ---------- Centrality definition ---------------
    // ------------------------------------------------
    //
    Int_t nEst = sizeof(fEventInfo_centEstStr)/sizeof(char*);
    fEventInfo_CentralityEstimates->Zero();  // matchEff counter
    if (fBeamType.CompareTo("A-A") == 0) { // PbPb
      if (MultSelection) {
        if(fSystCentEstimatetor == -1) fCentrality = MultSelection->GetMultiplicityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = MultSelection->GetMultiplicityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=MultSelection->GetMultiplicityPercentile(fEventInfo_centEstStr[i]);
      } else if (esdCentrality) {
        if(fSystCentEstimatetor == -1) fCentrality = esdCentrality->GetCentralityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = esdCentrality->GetCentralityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = esdCentrality->GetCentralityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=esdCentrality->GetCentralityPercentile(fEventInfo_centEstStr[i]);
      } else {
        std::cout << " Info::siweyhmi: Error: There is no cent info " << std::endl;
      }

    }
    //
    // pp
    if (fBeamType.CompareTo("p-p") == 0){
      if (MultSelection) {
        fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      } else if (esdCentrality) {
        fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      } else {
        std::cout << " Info::siweyhmi: Error: There is no cent info " << std::endl;
      }
    }
    //
    if (fUseCouts) {
      std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
      std::cout << " Info::siweyhmi: Event counter = " << fEventCountInFile << " - cent =  " << fCentrality << " = gid = " << gid << " = fEventGID = " << fEventGID << std::endl;
      std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
    }
  }
  fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
  //
  // in case small stat is enough
  if (fUseCouts) std::cout << " dddd " << std::endl; //CHANGE
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)!=0) return kFALSE;

  // if centrality estimation failed
  if (fUseCouts) std::cout << " eeee " << std::endl; //CHANGE
  if (fCentrality > 100.) return kFALSE;
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fESD && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fNumRealJets = 0;

    FillTPCdEdxReal();
    if (fDoEMCJet){
    FindJetsEMC();
    }
    if (fDoFastJet){
      FindJetsFJ();
    }
    FillEventTree();
    if (fisGoodIncEvent==1){
      ftotalNumIncEvents++;
    }
    if (fhasRealFJjet==1){
      ftotalNumRealJetEvents++;
    }
    if (fUseCouts)  std::cout << " Info::siweyhmi: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  //
  // full MC analysis
  //
  if (fMCtrue && fESD  && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fNumRealJets = 0;

    FillTreeMC();
    if (fDoEMCJet){
      FindJetsEMC();
    }
    if (fDoFastJet){
      FindJetsFJ();
    }
    FillEventTree();
    if (fUseCouts)  std::cout << " Info::siweyhmi: (full MC analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisJetHadro::FindJetsEMC()
{

  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsEMC ===== " << std::endl;
  //
  Double_t pT_sub_min = fjetMinPtSub;
  if (fUseCouts) cout << "Minimum jet pt after subtraction is " << fjetMinPtSub << endl;
  // Get the jet container
  fJetContainer = this->GetJetContainer("detJets");
  TString fRhoName = fJetContainer->GetRhoName();
  if (fUseCouts) cout << "Rho Name is " << fRhoName << endl;

  if (fJetContainer->GetRhoParameter()) fjetRhoVal = fJetContainer->GetRhoVal();
  if (fUseCouts) cout << "In FindJetsEMC Rho value is " << fjetRhoVal << endl;

  //
  // Get some jet container properties
  Int_t Njets         = fJetContainer->GetNJets();
  Int_t NAcceptedjets = fJetContainer->GetNAcceptedJets();
  Double_t jetRadius  = fJetContainer->GetJetRadius();
  Double_t jetEtaMin = fJetContainer->GetJetEtaMin();
  Double_t jetEtaMax = fJetContainer->GetJetEtaMax();
  UInt_t jetAcceptanceType = fJetContainer->GetAcceptanceType();
  //
  // loop over jets
  for(auto jet : fJetContainer->accepted())
  {
    //if (TMath::Abs(jet->Eta()) >= 0.9-jetRadius) continue; //this is not needed when using jetcontainer->accepted
    fhasAcceptedEMCjet = 1;
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Double_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Double_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters
    Bool_t IsJetMc = jet->IsMC(); //using >0.0 of jet area inside emcal
    Double_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Double_t jetMCPt = jet->MCPt();
    Double_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE); //bg sub pt

    if (jetptsub > pT_sub_min && JetAreaPt > fjetMinArea)
    {
      fhasRealEMCjet = 1;
    }

    if(!fTreeSRedirector) return;
    if (fFillEMCJet){
      (*fTreeSRedirector)<<"jetsEMC"<<
      "gid="       << fEventGID << //  global event ID
      "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
      "jetRhoVal=" << fjetRhoVal <<
      "jetpt="     << fJetPt    << // jetPt
      "jeteta="    << fJetEta   << // jetEta
      "jetphi="    << fJetPhi   << // jetPhi
      "jetArea="   << JetAreaPt  <<
      "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
      "particleMaxChargedPt="  << particleMaxChargedPt <<
      "jetptsub="  << jetptsub <<
      "cent="      << fCentrality  <<  //  centrality
      "\n";
    }

    if (jetptsub <= pT_sub_min || JetAreaPt <= fjetMinArea) continue;

    Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
    const AliMultiplicity *multObj = fESD->GetMultiplicity();
    Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliESDtrack* esdtrack = (AliESDtrack*)(particle);
      if (!esdtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      fTrackCutBits=0;  // reset the bits for the next track
      //
      if (!esdtrack->GetInnerParam()) continue;               // Ask if track is in the TPC
      //if (!fESDtrackCuts->AcceptTrack(esdtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      if (!(esdtrack->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      Double_t closestPar[3];
      GetExpecteds(esdtrack,closestPar);
      SetCutBitsAndSomeTrackVariables(esdtrack);
      Int_t tpcNcls = esdtrack->GetTPCncls();
      Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
      AliVMultiplicity *multiObj = fESD->GetMultiplicity();
      Int_t nITSClusters = 0;
      for(Int_t j=2;j<6;j++) nITSClusters += multiObj->GetNumberOfITSClusters(j);
      //
      UShort_t tpcFindableCls = esdtrack->GetTPCNclsF();
      UShort_t tpcSharedCls   = esdtrack->GetTPCnclsS();
      Double_t tofSignalTunedOnData = esdtrack->GetTOFsignalTunedOnData();
      Double_t length    = esdtrack->GetIntegratedLength();
      Double_t tofSignal = esdtrack->GetTOFsignal();
      Double_t beta = -.05;
      if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;

      // Fill track constituent information
      if (fFillJetsEMCConst){

        (*fTreeSRedirector)<<"jetsEMCconst"<<
        //
        "gid="       << fEventGID << //  global event ID
        "jetRadius=" << jetRadius << // jet Radius
        "Njets="     << Njets << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetEta
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCconst"<<
          "intrate="   << fIntRate              <<  // interaction rate
          "run="       << fRunNo <<                  // run Number
          "bField="    << fBField <<                 // magnetic filed
          "pileupbit=" << fPileUpBit <<              // flag for pileup selection
          "primMult="  << fNContributors <<          //  #prim tracks
          "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC cluster multiplicity
          "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
          "itsmult="   << itsNumberOfTracklets <<    // ITS multiplicity
          "itsclmult=" << nITSClusters <<    // ITS multiplicity
          "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
          "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
          "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
          "tpcSignalN="    << fTrackTPCSignalN <<  //  number of cl used in dEdx
          "lengthInActiveZone="  << fTrackLengthInActiveZone <<  //  track length in active zone
          "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
          "dcaz="      << fTrackDCAz            <<  // dca cut along z
          "ncltpc="    << fNcl                  <<  // number of clusters
          "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
          "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
          "missCl="    << fMissingCl            <<  // fraction of missing clusters
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "tofSignalTOD=" << tofSignalTunedOnData         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCconst"<<"\n";
      }
    }

  }

  if (fFillJetsEMCBG || fFillJetsEMCBGConst){
  fbgJetContainer = this->GetJetContainer("bgJets");
  //
  // Get some jet container properties
  Int_t Njetsbg         = fbgJetContainer->GetNJets();
  Int_t NAcceptedjetsbg = fbgJetContainer->GetNAcceptedJets();
  Double_t jetRadiusbg  = fbgJetContainer->GetJetRadius();
  Double_t jetEtaMinbg = fbgJetContainer->GetJetEtaMin();
  Double_t jetEtaMaxbg = fbgJetContainer->GetJetEtaMax();
  //
  // loop over jets
  for(auto jet : fbgJetContainer->accepted())
  {
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Double_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Double_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters
    Bool_t IsJetMc = jet->IsMC(); //using >0.0 of jet area inside emcal
    Double_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Double_t jetMCPt = jet->MCPt();
    Double_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE); //bg sub pt

    if (fFillJetsEMCBG){
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"jetsEMCBG"<<
    "gid="       << fEventGID << //  global event ID
    "jetRadius=" << jetRadiusbg << // jet Radius
    "Njets="     << Njetsbg << // Number of jets in event
    "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
    "jetEtaMinbg=" << jetEtaMinbg <<
    "jetEtaMaxbg=" << jetEtaMaxbg <<
    "jetRhoVal=" << fjetRhoVal <<
    "jetpt="     << fJetPt    << // jetPt
    "jeteta="    << fJetEta   << // jetEta
    "jetphi="    << fJetPhi   << // jetPhi
    "jetArea="   << JetAreaPt  <<
    "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
    "particleMaxChargedPt="  << particleMaxChargedPt <<
    "jetptsub="  << jetptsub <<
    "cent="      << fCentrality  <<  //  centrality
    "\n";
    }

    Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
    const AliMultiplicity *multObj = fESD->GetMultiplicity();
    Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliESDtrack* esdtrack = (AliESDtrack*)(particle);
      if (!esdtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      fTrackCutBits=0;  // reset the bits for the next track
      //
      if (!esdtrack->GetInnerParam()) continue;               // Ask if track is in the TPC
      //if (!fESDtrackCuts->AcceptTrack(esdtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      if (!(esdtrack->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(esdtrack);
      Double_t closestPar[3];
      GetExpecteds(esdtrack,closestPar);
      SetCutBitsAndSomeTrackVariables(esdtrack);
      Int_t tpcNcls = esdtrack->GetTPCncls();
      Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
      AliVMultiplicity *multiObj = fESD->GetMultiplicity();
      Int_t nITSClusters = 0;
      for(Int_t j=2;j<6;j++) nITSClusters += multiObj->GetNumberOfITSClusters(j);
      //
      UShort_t tpcFindableCls = esdtrack->GetTPCNclsF();
      UShort_t tpcSharedCls   = esdtrack->GetTPCnclsS();
      Double_t tofSignalTunedOnData = esdtrack->GetTOFsignalTunedOnData();
      Double_t length    = esdtrack->GetIntegratedLength();
      Double_t tofSignal = esdtrack->GetTOFsignal();
      Double_t beta = -.05;
      if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;

      // Fill track constituent information
      if (fFillJetsEMCBGConst){

        (*fTreeSRedirector)<<"jetsEMCBGconst"<<
        //
        "gid="       << fEventGID << //  global event ID
        "jetRadius=" << jetRadiusbg << // jet Radius
        "Njets="     << Njetsbg << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetPhi
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
        "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCBGconst"<<
          "intrate="   << fIntRate              <<  // interaction rate
          "run="       << fRunNo <<                  // run Number
          "bField="    << fBField <<                 // magnetic filed
          "pileupbit=" << fPileUpBit <<              // flag for pileup selection
          "primMult="  << fNContributors <<          //  #prim tracks
          "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC cluster multiplicity
          "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
          "itsmult="   << itsNumberOfTracklets <<    // ITS multiplicity
          "itsclmult=" << nITSClusters <<    // ITS multiplicity
          "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
          "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
          "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
          "tpcSignalN="    << fTrackTPCSignalN <<  //  number of cl used in dEdx
          "lengthInActiveZone="  << fTrackLengthInActiveZone <<  //  track length in active zone
          "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
          "dcaz="      << fTrackDCAz            <<  // dca cut along z
          "ncltpc="    << fNcl                  <<  // number of clusters
          "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
          "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
          "missCl="    << fMissingCl            <<  // fraction of missing clusters
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "tofSignalTOD=" << tofSignalTunedOnData         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCBGconst"<<"\n";
      }
    }

  }
  }

}
//________________________________________________________________________
void AliAnalysisJetHadro::FindJetsFJ()
{
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsFJ ===== " << std::endl;
  //
  // Create jetwrapper with the same settings used in FindJetsEMC
  int nJetRadiusBins = 1;
  int nSettings = 0;
  float fTrackPt = 0.15;
  float fTrackEta = 0.9;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;           // fixed
  float bgJetRadius = 0.2;              // fixed
  //
  float pT_sub_min = fjetMinPtSub;    // can be 40, 60, 80
  std::vector<float>  fJetRadius;
  std::vector<float>  fPtSubMin;
  fJetRadius.resize(nJetRadiusBins);
  fJetRadius[0] = 0.4;
  fJetRadius[1] = 0.2;
  fJetRadius[2] = 0.6;
  fPtSubMin.resize(nJetRadiusBins);

  for (int iset=-1; iset<nSettings; iset++){
    for (int iJetRadius=0; iJetRadius<nJetRadiusBins; iJetRadius++){
      Float_t jetAbsEtaCut = fTrackEta-fJetRadius[iJetRadius];   // fixed

      //
      //SOME CODE FROM NIMA
      AliFJWrapper *fFastJetWrapper;
      fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
      fFastJetWrapper->Clear();
      fFastJetWrapper->SetR(fJetRadius[iJetRadius]);
      fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
      fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme);
      fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
      fFastJetWrapper->SetGhostArea(fGhostArea);
      fFastJetWrapper->SetAreaType(fastjet::AreaType::active_area);
      fFastJetWrapper->SetMaxRap(1);
      std::vector<int> trackTTIndex;
      trackTTIndex.clear();
      //
      // Access and loop over container
      // auto tracks = this->GetTrackContainer("detTracks");
      // for (auto esdtrack : tracks->accepted()){...}
      //
      //std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
      std::vector<fastjet::PseudoJet> particlesEmbedded; //fill this with your event
      double particleEtaCut = 0.9;
      //
      // loop over esd tracks and add their four vector to wrapper --> identical to track container in EMC jet
      for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
        AliESDtrack* track = fESD->GetTrack(iTrack);
        if (!fESDtrackCutsLoose->AcceptTrack(track)) continue;
        if (!track->GetInnerParam()) continue;
        if (!(track->GetTPCsignalN()>0)) continue;
        if (iset==-1) {
          if (!fESDtrackCuts->AcceptTrack(track)) continue; // default cuts which should match EMC jets
        } else {
          Double_t closestPar[3];
          GetExpecteds(track,closestPar);
          SetCutBitsAndSomeTrackVariables(track);
          if (!GetSystematicClassIndex(fTrackCutBits,iset)) continue;
        }
        //
        if (track->Pt() < fTrackPt || TMath::Abs(track->Eta()) >= particleEtaCut) continue;
        fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
        particlesEmbedded.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));// TMath::Sqrt(track->P()*track->P()+0.13957*0.13957) ) );
      }
      //
      // background jet definitions
      fastjet::JetMedianBackgroundEstimator bgE;
      fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2) * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackPt, 1000.0); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
      bgE.set_selector(selectorBG);
      fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, bgJetRadius, fastjet::E_scheme, fastjet::Best); //define the kT jet finding which will do the average background estimation
      fastjet::GhostedAreaSpec ghostSpecBG(particleEtaCut, 1, fGhostArea); //this ghost area might be too small and increase processing time too much
      fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts, ghostSpecBG);
      fastjet::ClusterSequenceArea cluster_seq_BG(particlesEmbedded, jetDefBG, areaDefBG);
      std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
      if (jetsBG.size() > 0) {
        bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
        frhoFJ = bgE.rho();
      }
      if (fUseCouts) std::cout << "frhoFJ is " << frhoFJ << std::endl;
      //
      // start of background jet loop
      for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
        fastjet::PseudoJet jet = jetsBG[ijet];
        //if (jet.pt() < fTrackPt || jet.perp() > 1000.0) || TMath::Abs(jet.eta()) >= bgJetAbsEtaCut) continue; //redundant because of the selector
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - frhoFJ*jetArea;
        Int_t jetNum = jetsBG.size();
        Int_t Njets = jetsBG.size();

        if (fFilljetsFJBGTree)
        {
          (*fTreeSRedirector)<<"jetsFJBG"<<
          "gid="            << fEventGID << //  global event ID
          "syst="           << iset << //  syst setting
          "bjJetRadius="    << bgJetRadius << // jet Radius
          "bgJetAbsEtaCut=" << bgJetAbsEtaCut << //abs eta cut for jet
          "jetNum="         << jetNum <<    //  number of jets
          "jetpt="          << jetpt <<
          "jetphi="         << jetphi <<
          "jeteta="         << jeteta <<
          "jetptsub="       << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="          << frhoFJ << //event rho
          "jetArea="        << jetArea << //jet area
          "cent="           << fCentrality  <<  //  centrality
          "\n";
        }

        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
        Int_t nConstituents = constituents.size();

        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Float_t pt = constituent.pt();
          Float_t phi = constituent.phi();
          Float_t eta = constituent.eta();

          if (fFillJetsFJBGConst){
          (*fTreeSRedirector)<<"jetsFJBGconst"<<
          "pT="        << fPt                   <<  // transverse momentum
          "eta="       << fEta                  <<  //  eta
          "phi="       << fPhi                  <<  //  phi
          "\n";
          }
        }
      } // end of background jet loop
      //
      // background subtraction on the constituent level TODO
      // fastjet::contrib::ConstituentSubtractor subtractorConstituent(&bgE); //add the background estimator to the correct subtractor
      // subtractorConstituent.set_common_bge_for_rho_and_rhom(true); //CHECK : should not be the case since particles have mass
      // subtractorConstituent.set_max_standardDeltaR(0.25); // set the max event wise subtraction distance
      // particlesEmbeddedSubtracted = subtractorConstituent.subtract_event(particlesEmbedded, particleEtaCut); //perform subtraction and fill the subtracted event container
      //

      // run jet finder using wrapper
      fFastJetWrapper->Run();
      std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();
      // start of jet loop
      for (Int_t ijet=0; ijet<Int_t(jets.size()); ijet++){
        fastjet::PseudoJet jet = jets[ijet];
        if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - frhoFJ*jetArea;
        Int_t jetNum = jets.size();
        //
        std::vector<fastjet::PseudoJet> constituents(sorted_by_pt(fFastJetWrapper->GetJetConstituents(ijet)));
        Int_t nConstituents = constituents.size();
        if (nConstituents<1) continue;
        fastjet::PseudoJet highestpt_const = constituents[0];
        if (highestpt_const.pt() > 100.0) continue;
        fhasAcceptedFJjet = 1;

        if (jetptsub > pT_sub_min && jetArea > fjetMinArea)
        {
          fhasRealFJjet = 1;
          fNumRealJets +=1;
          ftotalJetArea += jetArea;
          ftotalNumRealJets +=1;
        }

        if (fFillFastJet){
          (*fTreeSRedirector)<<"jetsFJ"<<
          "gid="          << fEventGID << //  global event ID
          "syst="         << iset << //  syst setting
          "jetNum="       << jetNum <<    //  number of jets
          "jetpt="        << jetpt <<
          "jetphi="       << jetphi <<
          "jeteta="       << jeteta <<
          "nConst="       << nConstituents <<    //  global event ID
          "cent="         << fCentrality           <<  //  centrality
          "jetptsub="     << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="        << frhoFJ << //event rho
          "jetArea="      << jetArea << //jet area
          "\n";
        }

        fHistJet_ptsub_v_area->Fill(jetArea, jetptsub);

        if (jetArea <= fjetMinArea) continue;
        fHistJet_kin->Fill(jetptsub, jeteta, jetphi);
        fHistJet_moms->Fill(jetptsub, jetpt);

        if (jetptsub <= pT_sub_min) continue;

        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Int_t trackIndex = constituent.user_index();
          AliESDtrack* trackConst = fESD->GetTrack(trackIndex);

          //Track cuts start
          fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
          fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
          fTrackCutBits=0;  // reset the bits for the next track
          //
          // --------------------------------------------------------------
          //      Get relevant track info and set cut bits
          // --------------------------------------------------------------
          //
          Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackConst);
          Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(trackConst);
          Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackConst);
          Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackConst);
          Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackConst);
          if (!trackConst->GetInnerParam()) continue;               // Ask if track is in the TPC
          // if (!fESDtrackCutsLoose->AcceptTrack(trackConst))  continue;    // Loose cuts
          if (!(trackConst->GetTPCsignalN()>0)) continue;
          //
          // Get the track variables
          Double_t closestPar[3];
          GetExpecteds(trackConst,closestPar);
          SetCutBitsAndSomeTrackVariables(trackConst);
          Int_t tpcNcls = trackConst->GetTPCncls();
          Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
          Int_t nITSClusters = 0;
          AliVMultiplicity *multiObj = fESD->GetMultiplicity();
          for(Int_t j=2;j<6;j++) nITSClusters += multiObj->GetNumberOfITSClusters(j);

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9){
            if (fFillpTPC) fHistJetTracks_dEdx->Fill(fPtot,fTPCSignal,TMath::Abs(fEta));
            if (fFillp) fHistJetTracks_dEdx_p->Fill(fPVertex,fTPCSignal,TMath::Abs(fEta));
            if (fFillpT) fHistJetTracks_dEdx_pT->Fill(fPt,fTPCSignal,TMath::Abs(fEta));
            if (fFillpTPC) fHistJetTracks_moms->Fill(fPt,fPtot);
            if (fFillp) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
            fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
          }

          if (fFillJetsFJConst){

          (*fTreeSRedirector)<<"jetsFJconst"<<
          "gid="       << fEventGID << //  global event ID
          "syst="      << iset << //  syst setting
          "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
          "jetNum="    << jetNum <<    //  number of jets
          "jetpt="     << jetpt <<     //  global event ID
          "jetphi="    << jetphi <<    //  global event ID
          "jeteta="    << jeteta <<    //  global event ID
          "nConst="    << nConstituents <<    //  global event ID

          "jetptsub="  << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="     << frhoFJ << //event rho
          "jetArea="   << jetArea << //jet area

          "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
          "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts

          "dEdx="      << fTPCSignal            <<  //  dEdx of the track
          "sign="      << fSign                 <<  //  charge
          "ptot="      << fPtot                 <<  //  TPC momentum
          "p="         << fPVertex              <<  //  momentum at vertex
          "pT="        << fPt                   <<  // transverse momentum
          "eta="       << fEta                  <<  //  eta
          "cent="      << fCentrality           <<  //  centrality
          "phi="       << fPhi                  <<  //  phi

          "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
          "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
          "dEdxMeanPi=" << fDEdxPi              <<
          "dEdxSigmaPi=" << fSigmaPi            <<
          "dEdxMeanKa=" << fDEdxKa              <<
          "dEdxSigmaKa=" << fSigmaKa            <<
          "dEdxMeanPr=" << fDEdxPr              <<
          "dEdxSigmaPr=" << fSigmaPr            <<
          "dEdxMeanDe=" << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe;

          if (!fSmallOut){
            (*fTreeSRedirector)<<"jetsFJconst"<<
            "bit96="     << fBit96_base <<    // tight cuts of 2011 tuned data
            "bit128="    << fBit128 <<        // TPC only tracks cuts
            "bit768="    << fBit768 <<        // Hybrid track cuts
            "pixCut="    << ifDCAcutIfNoITSPixel <<    // cut: apply a DCAcut If No ITS Pixel
            "run="       << fRunNo <<                  // run Number
            "bField="    << fBField <<                 // magnetic filed
            "pileupbit=" << fPileUpBit <<              // flag for pileup selection
            "primMult="  << fNContributors <<          //  #prim tracks
            //
            "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
            "itsclmult=" << nITSClusters <<    // ITS multiplicity
            "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
            //
            "eventtime=" << fTimeStamp            <<  // event timeStamp
            "intrate="   << fIntRate              <<  // interaction rate

            "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
            "dcaz="      << fTrackDCAz            <<  // dca cut along z
            "ncltpc="    << fNcl                  <<  // number of clusters
            "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
            "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
            "missCl="    << fMissingCl            <<  // fraction of missing clusters

            "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
            "pitpcpid="  << fNSigmasPiTPC         <<
            "katpcpid="  << fNSigmasKaTPC         <<
            "prtpcpid="  << fNSigmasPrTPC         <<
            "detpcpid="  << fNSigmasDeTPC         <<

            "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
            "pitofpid="  << fNSigmasPiTOF         <<
            "katofpid="  << fNSigmasKaTOF         <<
            "prtofpid="  << fNSigmasPrTOF         <<
            "detofpid="  << fNSigmasDeTOF         <<

            "closestTPCPIDtype=" << closestPar[1]         << //particle type
            "closestTPCPIDmass=" << closestPar[2];
          }
        (*fTreeSRedirector)<<"jetsFJconst"<<"\n";
        }
        }
      } // end of jet loop

    }
  }
}
//________________________________________________________________________
void AliAnalysisJetHadro::FillTPCdEdxReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillTPCdEdxReal ===== " << std::endl;
  // --------------------------------------------------------------
  // Get the event
  AliVEvent *event=InputEvent();
  if (CountEmptyEvents(2)<1) return;
  fisGoodIncEvent = 1;
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
  //
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
    fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track) continue;
    //
    // --------------------------------------------------------------
    //      Get relevant track info and set cut bits
    // --------------------------------------------------------------
    //
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(track);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(track);
    if (!track->GetInnerParam()) continue;               // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    // Get the track variables
    Double_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track);
    Int_t tpcNcls = track->GetTPCncls();
    //
    // --------------------------------------------------------------
    //  Some print out
    // --------------------------------------------------------------
    //
    // Tree for the all cut variables
    if (fUseCouts && fEventCountInFile==5 && !fRunOnGrid) {
      std::cout << " Info::siweyhmi: CutBinMap --> " <<fTrackTPCCrossedRows << " " << fTrackChi2TPC << " " <<  fTrackNewITScut  << std::endl;
      PrintNumInBinary(fTrackCutBits);
    }
    //
    // --------------------------------------------------------------
    //   Fill the trees
    // --------------------------------------------------------------
    //
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    //
    // different dca cuts
    // TMath::Abs(fTrackDCAxy)< 0.3
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);  // 10h tuned loose cut
    //
    UShort_t tpcFindableCls = track->GetTPCNclsF();
    UShort_t tpcSharedCls = track->GetTPCnclsS();
    //
    Double_t tofSignalTunedOnData = track->GetTOFsignalTunedOnData();
    Double_t length = track->GetIntegratedLength();
    Double_t tofSignal = track->GetTOFsignal();
    Double_t beta = -.05;
    if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;

    if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < 0.9){
      if (fFillpTPC) fHistIncTracks_dEdx->Fill(fPtot,fTPCSignal,TMath::Abs(fEta));
      if (fFillp) fHistIncTracks_dEdx_p->Fill(fPVertex,fTPCSignal,TMath::Abs(fEta));
      if (fFillpT) fHistIncTracks_dEdx_pT->Fill(fPt,fTPCSignal,TMath::Abs(fEta));
      if (fFillpTPC) fHistIncTracks_moms->Fill(fPt,fPtot);
      if (fFillp) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
      fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);

      fHistIncTracks_mpi->Fill(fPtot,fDEdxPi,TMath::Abs(fEta));
      fHistIncTracks_spi->Fill(fPtot,fSigmaPi,TMath::Abs(fEta));

      fHistIncTracks_mel->Fill(fPtot,fDEdxEl,TMath::Abs(fEta));
      fHistIncTracks_sel->Fill(fPtot,fSigmaEl,TMath::Abs(fEta));

      fHistIncTracks_mka->Fill(fPtot,fDEdxKa,TMath::Abs(fEta));
      fHistIncTracks_ska->Fill(fPtot,fSigmaKa,TMath::Abs(fEta));

      fHistIncTracks_mpr->Fill(fPtot,fDEdxPr,TMath::Abs(fEta));
      fHistIncTracks_spr->Fill(fPtot,fSigmaPr,TMath::Abs(fEta));
    }


    //

	  if (fFillIncTracks && !fFillOnlyHists)
    {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"tracks"<<
      //
      "gid="       << fEventGID             <<  //  global event ID
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit768="    << fBit768 <<        // Hybrid track cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "phi="       << fPhi                  <<  //  phi
      "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi=" << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa=" << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr=" << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "cent="      << fCentrality;
      if (!fSmallOut){
        (*fTreeSRedirector)<<"tracks"<<
        "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
        "intrate="   << fIntRate              <<  // interaction rate
        "eventtime=" << fTimeStamp            <<  // event timeStamp
        "bit96="     << fBit96_base <<    // tight cuts of 2011 tuned data
        "bit128="    << fBit128 <<        // TPC only tracks cuts
        "pixCut="    << ifDCAcutIfNoITSPixel <<    // cut: apply a DCAcut If No ITS Pixel
        "run="       << fRunNo <<                  // run Number
        "bField="    << fBField <<                 // magnetic filed
        "pileupbit=" << fPileUpBit <<              // flag for pileup selection
        "primMult="  << fNContributors <<          //  #prim tracks
        "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC cluster multiplicity
        "dcabase="   << dcaBaseCut <<  //  TPC multiplicity
        "dca10h="    << dca10h <<  //  TPC multiplicity
        "dca11h="    << dca11h <<  //  TPC multiplicity
        "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
        "itsmult="   << itsNumberOfTracklets <<    // ITS multiplicity
        "itsclmult=" << nITSClusters <<    // ITS multiplicity
        "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
        "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
        "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
        "tpcSignalN="    << fTrackTPCSignalN <<  //  number of cl used in dEdx
        "lengthInActiveZone="  << fTrackLengthInActiveZone <<  //  track length in active zone
        "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
        "dcaz="      << fTrackDCAz            <<  // dca cut along z
        "ncltpc="    << fNcl                  <<  // number of clusters
        "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
        "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
        "missCl="    << fMissingCl            <<  // fraction of missing clusters
        "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
        "pitpcpid="  << fNSigmasPiTPC         <<
        "katpcpid="  << fNSigmasKaTPC         <<
        "prtpcpid="  << fNSigmasPrTPC         <<
        "tofSignal=" << tofSignal         <<
        "tofSignalTOD=" << tofSignalTunedOnData         <<
        "prtofpid="  << fNSigmasPrTOF<<
        "closestTPCPIDtype=" << closestPar[1]         << //particle type
        "closestTPCPIDmass=" << closestPar[2]         << //particle mass
        "dEdxMeanDe=" << fDEdxDe              <<
        "dEdxSigmaDe=" << fSigmaDe            <<
        "beta="      << beta;
      }
      (*fTreeSRedirector)<<"tracks"<<"\n";
    }
    //
    // --------------------------------------------------------------
    //  Fill the THnSparseF for the Expected values form PID response
    // --------------------------------------------------------------
    //
    // define acceptance of interest
    Bool_t etaAcc  = (fEta >=fEtaDown       && fEta<=fEtaUp);
    Bool_t momAcc  = (fPVertex>=fMomDown    && fPVertex<=fMomUp);
    Bool_t dEdxAcc = (fTPCSignal>=fDEdxDown && fTPCSignal<=fDEdxUp);
    Bool_t fAcceptance = (etaAcc && momAcc && dEdxAcc);
    Bool_t nSigmasElTPCCut = (TMath::Abs(fNSigmasElTPC)<2);
    Bool_t nSigmasPiTPCCut = (TMath::Abs(fNSigmasPiTPC)<2);
    Bool_t nSigmasKaTPCCut = (TMath::Abs(fNSigmasKaTPC)<2);
    Bool_t nSigmasPrTPCCut = (TMath::Abs(fNSigmasPrTPC)<2);
    Bool_t nSigmasDeTPCCut = (TMath::Abs(fNSigmasDeTPC)<2);
    Bool_t nSigmaTPCall = (nSigmasElTPCCut || nSigmasPiTPCCut || nSigmasKaTPCCut || nSigmasPrTPCCut || nSigmasDeTPCCut);
    Bool_t ndEdxTPCall  = (fDEdxEl>20 || fDEdxPi>20 || fDEdxKa>20 || fDEdxPr>20 || fDEdxDe>20);

    if(fAcceptance && !fMCtrue){
      if (fFillExpecteds && fEvent < 5 && (fNSigmasPiTPC >= 3 || (fNSigmasPiTPC < 3 && fRandom.Rndm() < 0.001))) {
        Double_t sign = static_cast<Double_t>(fSign);
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"expecteds"<<
        "cent="                 << fCentrality   <<  // centrality
        "sign="                 << sign          <<
        "ptot="                 << fPtot          <<  // momentum
        "eta="                  << fEta           <<  // eta
        "phi="                  << fPhi           <<  // phi
        "dEdxEl="               << fDEdxEl       <<
        "dEdxPi="               << fDEdxPi       <<
        "dEdxKa="               << fDEdxKa       <<
        "dEdxPr="               << fDEdxPr       <<
        "dEdxDe="               << fDEdxDe       <<
        "sigmaEl="               << fSigmaEl     <<
        "sigmaPi="               << fSigmaPi     <<
        "sigmaKa="               << fSigmaKa     <<
        "sigmaPr="               << fSigmaPr     <<
        "sigmaDe="               << fSigmaDe     <<
        "dEdx="                 << fTPCSignal     <<
        //
        "\n";
     }
    }

  }// end of track loop

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillTreeMC()
{
  //
  Int_t trackOrigin = -10;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillTreeMC ===== " << std::endl;
  //
  // ======================================================================
  // ------   reconstructed MC particles with dEdx information-------------
  // ======================================================================
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // Esd track
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    //
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
    Int_t pdg = trackMCgen->Particle()->GetPdgCode();
    //
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    //
    // check the origin of the track
    trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(lab))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(lab))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(lab)) trackOrigin = 2;
    if (trackOrigin<-1) continue; // TODO
    //
    // Track cuts from dtector
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(trackReal);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    //
    if (!trackReal->GetInnerParam()) continue;     // TODO        // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // TODO
    if (!(trackReal->GetTPCsignalN()>0)) continue; // TODO
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
    if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
    if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
    if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
    if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
    if (iPart == -10) continue; // TODO
    //
    Double_t closestPar[3];
    GetExpecteds(trackReal,closestPar);
    SetCutBitsAndSomeTrackVariables(trackReal);
    //
    if (trackReal-> GetInnerParam()){
      fPtotMC       = trackReal->GetInnerParam()->GetP();
      fTPCSignalMC  = trackReal->GetTPCsignal();
    }
    fEtaMC        = trackReal->Eta();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    fMissingCl    = trackReal->GetTPCClusterInfo(3,0,0,159);
    //
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    //
    // different dca cuts
    // TMath::Abs(fTrackDCAxy)< 0.3
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPtMC,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPtMC,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPtMC,1.01);  // 10h tuned loose cut
    //
    UShort_t tpcFindableCls = trackReal->GetTPCNclsF();
    UShort_t tpcSharedCls = trackReal->GetTPCnclsS();
    Float_t dca[2], covar[3];
    trackReal->GetImpactParameters(dca, covar);
    Double_t tofSignalTunedOnData = trackReal->GetTOFsignalTunedOnData();
    Double_t length = trackReal->GetIntegratedLength();
    Double_t tofSignal = trackReal->GetTOFsignal();
    Double_t beta = -.05;
    if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;

    // Bool_t settings[17];
    TVectorF settings(17);
    for (Int_t i = 0; i < 17;i++) {
      settings[i] = (Float_t) GetSystematicClassIndex(fTrackCutBits, i);
    }

    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"fTreeMC"<<
    "isample="   << sampleNo <<                // sample id for subsample method
    "orig="     << trackOrigin <<   // origin of the track
    "part="      << iPart <<
    "gid="       << fEventGID <<  //  global event ID
    "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
    "cutBit="    << fTrackCutBits <<  //  Systematic Cuts
    "settings.="    << &settings <<  //  Systematic settings
    "sign="      << fSignMC <<         // sign
    "ptot="      << fPtotMC <<         // tpc momentum
    "p="         << pMC <<             // vertex momentum
    "pT="        << fPtMC <<           // transverse momentum
    "eta="       << fEtaMC <<          // mc eta
    "phi="       << fPhiMC <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "centimp="   << fCentImpBin <<
    "vZ="        << fVz <<
    "nsigmatofka="  << fNSigmasKaTOF         <<  // interaction rate
    "nsigmatofpr="  << fNSigmasPrTOF         <<  // interaction rate
    //
    "impPar="       << fMCImpactParameter <<      // impact parameter taken from MC event header
    "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
    "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
    "tpcSignalN="           << fTrackTPCSignalN            <<  //  number of cl used in dEdx
    "lengthInActiveZone="           << fTrackLengthInActiveZone            <<  //  track length in active zone
    "tofSignal=" << tofSignal         <<
    "tofSignalTOD=" << tofSignalTunedOnData         <<
    "beta=" << beta         <<
    "dcaxy="     << fTrackDCAxy <<
    "dcaz="      << fTrackDCAz <<
    "cRows="     << fTrackTPCCrossedRows  <<
    "chi2tpc="   << fTrackChi2TPC         <<
    "defCut="    << ifDefaultCuts <<  // default cut
    "bit96="     << fBit96_base <<  // run Number
    "bit128="    << fBit128 <<  // run Number
    "bit768="    << fBit768 <<  // run Number
    "pixCut="    << ifDCAcutIfNoITSPixel <<  // run Number
    "run="       << fRunNo <<  // run Number
    "bField="    << fBField <<  // run Number
    "pileupbit=" << fPileUpBit <<
    "primmult="  << fNContributors <<  //  #prim tracks
    "ncltpc="    << fNcl                  <<  //  centrality
    "ncltpccorr="<< fNclCorr              <<  //  centrality
    "missCl="    << fMissingCl <<
    "chi2tpccorr=" << fTrackChi2TPCcorr         <<
    "dcabase="  << dcaBaseCut <<  //  TPC multiplicity
    "dca10h="   << dca10h <<  //  TPC multiplicity
    "dca11h="   << dca11h <<  //  TPC multiplicity
    "fCdd="      << covar[0] <<
    "fCdz="      << covar[1] <<
    "fCzz="      << covar[2] <<
    "tpcpileup=" << isTPCPileup <<
    "itspileup=" << isITSPileup <<
    "tpcmult="   << fTPCMult <<  //  TPC multiplicity
    "itsmult="   << itsNumberOfTracklets <<
    "itsclmult="   << nITSClusters <<    // ITS multiplicity
    "tpcclmult="   << nTPCClusters <<    // ITS multiplicity
    "\n";

  } // ======= end of track loop for MC dEdx filling =======

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillEventTree()
{
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEventTree ===== " << std::endl;
  (*fTreeSRedirector)<<"jeteventInfo"<<
  //"gid="       << fEventGID << //  global event ID
  "rhoFJ="      << frhoFJ << //event rho
  "rhoEMC="  << fjetRhoVal <<
  "cent="      << fCentrality  <<  //  centrality
  "isGoodIncEvent="   << fisGoodIncEvent <<
  "hasAcceptedFJjet="   << fhasAcceptedFJjet <<
  "hasRealFJjet="   << fhasRealFJjet <<
  "hasAcceptedEMCjet="   << fhasAcceptedEMCjet <<
  "hasRealEMCjet="   << fhasRealEMCjet <<
  "NumRealJets="   << fNumRealJets <<
  "\n";
}
//________________________________________________________________________
void AliAnalysisJetHadro::GetExpecteds(AliESDtrack *track, Double_t closestPar[3])
{

  //
  // bettaGamma is not well deifned below bg=0.01 --> below 200MeV protons and deuterons
  Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  // if (ptotForBetaGamma<ptotForBetaGammaThr) return;
  //
  // --------------------------------------------------------------
  //  Calculates expected sigma and dEdx for a given track and returns colesest expected particle and its index
  // --------------------------------------------------------------
  //
  fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //
  fNSigmasElTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,     fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,     fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  if (ptotForBetaGamma>ptotForBetaGammaThr){
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,   fPIDResponse->GetTOFResponse().GetTimeZero());
    fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron, fPIDResponse->GetTOFResponse().GetTimeZero());
  }

  //
  //
  Int_t nSigmaTmp = (fEventInfo) ? 10000 : 3; //(fEventInfo) ? 10000 : 2
  //
  // Electron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) { //TMath::Abs(fNSigmasElTPC)<nSigmaTmp
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) { //TMath::Abs(fNSigmasElTPC)<nSigmaTmp
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) { ///TMath::Abs(fNSigmasElTPC)<nSigmaTmp
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) { ///TMath::Abs(fNSigmasElTPC)<nSigmaTmp
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) { ///TMath::Abs(fNSigmasElTPC)<nSigmaTmp
    fDEdxDe  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    fSigmaDe = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
  }
  //
  // --------------------------------------------------------------
  //  Find Closest dEdx its corresponding particle and mass
  // --------------------------------------------------------------
  //
  Float_t values[] = {fDEdxEl, fDEdxPi, fDEdxKa, fDEdxPr, fDEdxDe};
  Float_t tpcdEdx = track->GetTPCsignal();
  Float_t smallestDiff = TMath::Abs(tpcdEdx - values[0]);
  Int_t closestIndex = 0;
  for (Int_t i = 0; i < 5; i++) {
    Double_t currentDiff = TMath::Abs(tpcdEdx - values[i]);
    if (currentDiff < smallestDiff) {
      smallestDiff = currentDiff;
      closestIndex = i;
    }
  }
  //
  // TF1 f1("f1","AliExternalTrackParam::BetheBlochAleph(x/0.1)",0.1,5)
  Float_t partMass = 0.;
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  if (closestIndex == 0 ) partMass = pdg->GetParticle(kPDGel)->Mass();  // GetParticle("e+")
  if (closestIndex == 1 ) partMass = pdg->GetParticle(kPDGpi)->Mass();  // GetParticle("pi+")
  if (closestIndex == 2 ) partMass = pdg->GetParticle(kPDGka)->Mass();  // GetParticle("K+")
  if (closestIndex == 3 ) partMass = pdg->GetParticle(kPDGpr)->Mass();  // GetParticle("proton")
  if (closestIndex == 4 ) partMass = 2.01410178;                        // pdg->GetParticle(1000010020)->Mass();
  //
  closestPar[0]=values[closestIndex];
  closestPar[1]=closestIndex;
  closestPar[2]=partMass;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::CheckIfFromResonance(Int_t mcType, AliMCParticle *trackMCgen, Int_t trackIndex, Bool_t parInterest, Double_t ptot, Double_t eta, Double_t cent, Bool_t fillTree)
{

  //
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  //
  TObjString momName="xxx";
  Int_t labMom = trackMCgen->Particle()->GetFirstMother();
  Int_t pdgMom = 0;
  if ( (labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks()) ){
    pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
    momName = fMCStack->Particle(labMom)->GetName();
  }
  //
  //Check if the particle is in the black list of resonances
  for (Int_t ires=0;ires<fNResBins;ires++){

    if ( fResonances[ires]=="xxx" ){
      if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;}
    } else {
      if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
    }
  }
  //
  // dump resonance info
  Int_t pdg = trackMCgen->Particle()->GetPdgCode();
  TObjString parName(trackMCgen->Particle()->GetName());
  if(fEventCountInFile==2 && !fRunOnGrid && !fSmallOut && fillTree) {
    if(!fTreeSRedirector) return kFALSE;
    (*fTreeSRedirector)<<"jetResonance"<<
    "acceptRes="   << acceptRes   <<
    "mcType="      << mcType       <<        // lower edge of momentum bin
    "ptot="        << ptot       <<          // lower edge of momentum bin
    "eta="         << eta     <<             // lower edge of eta bin
    "cent="        << cent        <<         // cent bin
    "parInterest=" << parInterest <<         // only pi, ka, and proton
    "pdg="         << pdg         <<         // pdg of prim particle
    "lab="         << trackIndex  <<         // index of prim particle
    "pdgMom="      << pdgMom      <<         // pdg of mother
    "labMom="      << labMom      <<         // index of mother
    "parName.="    << &parName    <<         //  full path - file name with ESD
    "momName.="    << &momName    <<         //  full path - file name with ESD
    "\n";
  }

  return acceptRes;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp)
{

  Bool_t sisterInAcceptance = kTRUE;
  Bool_t motherInAcceptance = kTRUE;

  Int_t labMom = trackMCgen->GetMother();
  Int_t pdgMom = 0;
  if ((labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks())){
    //
    // Check if the mother is also in the acceptance
    AliMCParticle *momTrack = (AliMCParticle *)fMCEvent->GetTrack(labMom);
    pdgMom = momTrack->Particle()->GetPdgCode();
    if ( pdgMom!=0 ) {
      Int_t    motherSign = momTrack->Charge();
      Double_t motherEta  = (fRapidityType==0) ? momTrack->Eta() :  momTrack->Y();
      Double_t motherMom  = momTrack->P();
      Bool_t etaAcc  = (motherEta>=etaLow && motherEta<etaUp);
      Bool_t momAcc  = (motherMom>=pDown  && motherMom<pUp  );
      if ( !(etaAcc && momAcc) ) motherInAcceptance=kFALSE;
    }
    //
    // Check if the sister is also in the acceptance
    Int_t labSister = momTrack->GetDaughterLast();
    if ((labSister>=0) && (labSister < fMCEvent->GetNumberOfTracks())){
      AliMCParticle *sisterTrack = (AliMCParticle *)fMCEvent->GetTrack(labSister);
      Int_t    sisterSign = sisterTrack->Charge();
      Double_t sisterEta  = (fRapidityType==0) ? sisterTrack->Eta() :  sisterTrack->Y();
      Double_t sisterMom  = sisterTrack->P();
      Bool_t etaAcc  = (sisterEta>=etaLow && sisterEta<etaUp);
      Bool_t momAcc  = (sisterMom>=pDown  && sisterMom<pUp  );
      if ( !(sisterSign!=0 && etaAcc && momAcc) ) sisterInAcceptance=kFALSE;
    }
  }
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  if ( pdgMom!=0 &&  fSisterCheck==0 ) acceptRes = kFALSE;                          // in anycase if mother exist      reject particle
  if ( pdgMom!=0 &&  motherInAcceptance &&  sisterInAcceptance && fSisterCheck==1 ) acceptRes = kFALSE;   // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==2 ) acceptRes = kFALSE;   // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==3 ) acceptRes = kTRUE;    // if sister and mother are in acc accept particle
  //
  if ( pdgMom!=0 && !motherInAcceptance && !sisterInAcceptance && fSisterCheck==4 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  sisterInAcceptance && fSisterCheck==5 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  fSisterCheck==6 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  return acceptRes;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::ApplyDCAcutIfNoITSPixel(AliESDtrack *track)
{

  Float_t p[2],cov[3];
  track->GetImpactParameters(p,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);

  fIsITSpixel01 = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters = track->GetNumberOfITSClusters();

  if (!cov[0] || !cov[2]) {
    return kFALSE;
  } else {
    fPrimRestriction = TMath::Sqrt((p[0]*p[0])/cov[0] + (p[1]*p[1])/cov[2]);
    return (fPrimRestriction<2 && fNITSclusters>2) || (fPrimRestriction<5 && fIsITSpixel01);
  }

}
//________________________________________________________________________
UInt_t AliAnalysisJetHadro::SetCutBitsAndSomeTrackVariables(AliESDtrack *track)
{
  //
  // Set some track variables
  //
  //
  // --------------------------------------------------------------
  //  calculate some variables by hand
  // --------------------------------------------------------------
  //
  // Double_t p[3];
  // track->GetPxPyPz(p);
  // Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  // Double_t pt       = TMath::Sqrt(p[0]*p[0] + p[1]*p[1]);
  // Double_t mass     = track->GetMass();  // assumed to be pion mass
  // Double_t energy   = TMath::Sqrt(mass*mass + momentum*momentum);
  // Float_t eta = -100.;
  // Float_t rap   = -100.;
  // if((momentum != TMath::Abs(p[2]))&&(momentum != 0)) eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  // if((energy != TMath::Abs(p[2]))&&(energy != 0))     rap = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
  //
  // --------------------------------------------------------------
  //  some extra getters
  // --------------------------------------------------------------
  //
  Double_t goldenChi2 = 0.;
  fPVertex = track->P();
  fTheta=track->Theta();
  fSign= track->GetSign();
  fPx  =track->Px();
  fPy  =track->Py();
  fPz  =track->Pz();
  fPt  =track->Pt();
  fY   =track->Y();
  fPhi =track->Phi()-TMath::Pi();
  fEta = track->Eta();
  Float_t pv[2],cov[3];
  track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  fTrackDCAxy = pv[0];
  fTrackDCAz  = pv[1];

  //
  // TPC related quantities
  Bool_t cleanDeTPC = kFALSE;
  if (track->GetInnerParam()){
    fPtot      = track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    fTrackTPCSignalN     = track->GetTPCsignalN();
    fTrackTPCCrossedRows = Float_t(track->GetTPCCrossedRows());
    fTPCShared = track->GetTPCnclsS();
    fTPCFindable = track->GetTPCNclsF();
    fMissingCl = track->GetTPCClusterInfo(3,0,0,159);
    goldenChi2 = track->GetChi2TPCConstrainedVsGlobal(fVertex);
    fTrackLengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    //
    Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t nSigmasDeTPC = (ptotForBetaGamma>ptotForBetaGammaThr) ? fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron) : 0.;
    cleanDeTPC = ((TMath::Abs(nSigmasDeTPC)<=2.));
    fNcl       = track->GetTPCncls();
    // fTrackChi2TPC  = (fNcl>0) ? TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/fNcl)) : -1;  // ???
    fTrackChi2TPC  = (fNcl>0) ? TMath::Abs(track->GetTPCchi2()/fNcl) : -1;  // ???
    //
    // correct for the missing clusters
    fNclCorr = fNcl;
    fTrackChi2TPCcorr = fTrackChi2TPC;
    //
    // --------------------------------------------------------------
    //      Bayesian PID part
    // --------------------------------------------------------------
    //
    if (ptotForBetaGamma>ptotForBetaGammaThr) {
      fPIDCombined->SetDefaultTPCPriors();
      Double_t probTPC[AliPID::kSPECIES]={0.};
      Double_t probTOF[AliPID::kSPECIES]={0.};
      // Get TPC probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      fTrackProbPiTPC = probTPC[AliPID::kPion];
      fTrackProbKaTPC = probTPC[AliPID::kKaon];
      fTrackProbPrTPC = probTPC[AliPID::kProton];
      // fTrackProbDeTPC = probTPC[AliPID::kDeuteron];
      // Get TOF probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
      fTrackProbPiTOF = probTOF[AliPID::kPion];
      fTrackProbKaTOF = probTOF[AliPID::kKaon];
      fTrackProbPrTOF = probTOF[AliPID::kProton];
    }
  }
  //
  // its restrictions
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);
  fIsITSpixel01    = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters    = track->GetNumberOfITSClusters();
  if (cov[0]>0 && cov[1]>0){
    fPrimRestriction = TMath::Sqrt((pv[0]*pv[0])/cov[0] + (pv[1]*pv[1])/cov[2]);
  }
  //
  fTrackRequireITSRefit  = track->IsOn(AliESDtrack::kITSrefit); // track->IsOn(AliESDtrack::kTPCrefit);
  fTrackIsFirstITSlayer  = track->HasPointOnITSLayer(0);
  fTrackIsSecondITSlayer = track->HasPointOnITSLayer(1);
  fTrackNewITScut        = ApplyDCAcutIfNoITSPixel(track);
  //
  Double_t nclsTRD     = (Float_t)track->GetTRDncls();
  Double_t TOFSignalDx = track->GetTOFsignalDx();
  Double_t TOFSignalDz = track->GetTOFsignalDz();
  //
  //
  Float_t nSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,    fPIDResponse->GetTOFResponse().GetTimeZero());
  Float_t nSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,  fPIDResponse->GetTOFResponse().GetTimeZero());
  Float_t nSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron,fPIDResponse->GetTOFResponse().GetTimeZero());
  Bool_t cleanPrTOF = ((TMath::Abs(nSigmasPrTOF)<=fEffMatrixNSigmasTOF));
  Bool_t cleanDeTOF = ((TMath::Abs(nSigmasDeTOF)<=fEffMatrixNSigmasTOF));
  Bool_t cleanKaTOF = ((TMath::Abs(nSigmasKaTOF)<=fEffMatrixNSigmasTOF));
  Bool_t cleanKaTOFTRD = ((TMath::Abs(nSigmasKaTOF)<=1.2) && TOFSignalDz<1. && TOFSignalDx<1. && nclsTRD>100);
  //
  // Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);
  // Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);
  Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);
  Bool_t dcaLoose   = TMath::Abs(fTrackDCAxy)<0.4;  // 10h tuned loose cut

  //
  // Systematic settings
  fTrackCutBits=0;
  //
  // Crossed rows
  if (fCorrectForMissCl==1){
    if (fNcl>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fNcl>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNcl>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
  } else if (fCorrectForMissCl==2){
    if (fNclCorr>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fNclCorr>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNclCorr>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
    // cout <<  "ncls  = " << fNclCorr <<  " --- " << fNcl << endl;
  } else {
    if (fTrackTPCCrossedRows>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fTrackTPCCrossedRows>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fTrackTPCCrossedRows>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
  }
  //
  // Special treatment of the 2018 pass3 and 2015 pass2 data
  // Chi2 TPC
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    if (fTrackChi2TPC<2.2) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);
    if (fTrackChi2TPC<2.5) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);
    if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);
  } else {
    //
    // correction for missing clusters
    if (fCorrectForMissCl==2){
      if (fTrackChi2TPCcorr<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPCcorr<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPCcorr<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    } else {
      if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPC<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPC<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    }
  }
  //
  // Shared TPC clusters
  Bool_t sharedCls = kFALSE;
  Bool_t sharedClsLoose = kFALSE;
  if (fTrackTPCCrossedRows > 0 && fNcl > 0) {
    sharedCls = (fTPCShared / fTrackTPCCrossedRows < 0.25) && (fTPCShared / static_cast<Float_t>(fNcl) < 0.3);
    sharedClsLoose = fTPCShared / fTrackTPCCrossedRows < 0.25;
  }
  if (sharedCls) (fTrackCutBits |= 1 << kSharedCls);
  if (sharedClsLoose) (fTrackCutBits |= 1 << kSharedClsLoose);
  //
  // Found TPC clusters
  if (fTPCFindable > 0) {
    if (fTrackTPCCrossedRows / fTPCFindable > 0.80) (fTrackCutBits |= 1 << kFindableCls);
    if (fTrackTPCCrossedRows / fTPCFindable > 0.85) (fTrackCutBits |= 1 << kFindableClsTight);
    if (fTrackTPCCrossedRows / fTPCFindable > 0.75) (fTrackCutBits |= 1 << kFindableClsLoose);
  }
  //
  // DCAxy
  if (dcaBaseCut) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDep);
  if (dcaLoose)   (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepLarge);
  //
  // DCAz
  if (TMath::Abs(fTrackDCAz)<0.15) (fTrackCutBits |= 1 << kVertexZSmall);
  if (TMath::Abs(fTrackDCAz)<1.00) (fTrackCutBits |= 1 << kVertexZ);
  //
  // Event vertex z
  if (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15) (fTrackCutBits |= 1 << kEventVertexZ);
  if (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1 ) (fTrackCutBits |= 1 << kEventVertexZLarge);
  //
  // track length cut --> dangerous cuts because it creates momentum dependent efficiency
  if (fTrackLengthInActiveZone>=90)  (fTrackCutBits |= 1 << kActiveZone);
  //
  // NCl in dEdx calculation
  if (fTrackTPCSignalN>=60) (fTrackCutBits |= 1 << kTPCSignalNSmall);
  if (fTrackTPCSignalN>=70) (fTrackCutBits |= 1 << kTPCSignalN);
  if (fTrackTPCSignalN>=80) (fTrackCutBits |= 1 << kTPCSignalNLarge);
  //
  // pile-up
  if (!fMCtrue) { // real data
    if (fPileUpBit & 1 << 0) (fTrackCutBits |= 1 << kPileup);
    if (fPileUpBit & 1 << 1) (fTrackCutBits |= 1 << kPileupLoose);
  } else {
    if (!fIsMCPileup) (fTrackCutBits |= 1 << kPileup);
    fTrackCutBits |= 1 << kPileupLoose; // fill for all events if no pileup rejection
  }
  //
  // B field polarity
  if (fBField > 0) (fTrackCutBits |= 1 << kBFieldPos);
  if (fBField < 0) (fTrackCutBits |= 1 << kBFieldNeg);
  //
  // --------------------------------------------------------------------
  //                    Clean sample selections
  // --------------------------------------------------------------------
  //
  // variable nsigma TOF protons and kaons for amplitude estimation
  if (cleanPrTOF) (fTrackCutBits |= 1 << kCleanPrTOF);
  if (cleanKaTOF) (fTrackCutBits |= 1 << kCleanKaTOF);
  //
  // Clean Kaons protons and deuterons
  if (cleanKaTOFTRD)        (fTrackCutBits |= 1 << kCleanKaTOFTRD);
  if (fTrackProbKaTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbKaTOF);
  if (fTrackProbPrTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbPrTOF);
  if (cleanDeTOF && cleanDeTPC) (fTrackCutBits |= 1 << kCleanDeTOF);
  //
  return fTrackCutBits;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::GetSystematicClassIndex(UInt_t cut,Int_t syst)
{
  /*
  syst:
  0 -->  Reference
  1 -->  CRows70
  2 -->  CRows90
  3 -->  ActiveZone
  4 -->  Chi2TPCSmall
  5 -->  Chi2TPCLarge
  6 -->  kMaxDCAToVertexXYPtDepLarge
  7 -->  kVertexZSmall
  8 -->  kEventVertexZLarge
  9 -->  kSharedCls
  10 -->  kFindableClsTight
  11 -->  kFindableClsLoose
  12 -->  kPileupLoose
  13 -->  kBFieldPos
  14 -->  kBFieldNeg
  15 -->  kTPCSignalNSmall
  16 -->  kTPCSignalNLarge
  */

  std::vector<Int_t> fCutArr;

  switch(syst) {

    case kCutReference:   // 0 -->  Reference
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutCrossedRowsTPC70:  // 1 -->  kNCrossedRowsTPC70
    {
      fCutArr = {kNCrossedRowsTPC70,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutCrossedRowsTPC90:  // 2 -->  kNCrossedRowsTPC90
    {
      fCutArr = {kNCrossedRowsTPC90,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutActiveZone:  // 3 -->  kActiveZone
    {
      fCutArr = {kActiveZone,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCSmall:   // 4 -->  kMaxChi2PerClusterTPCSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCSmall, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCLarge:   // 5 -->  kMaxChi2PerClusterTPCLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCLarge, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepLarge:   // 6 -->  kMaxDCAToVertexXYPtDepLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutVertexZSmall:   // 7 -->  kVertexZSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutEventVertexZLarge:  // 8 -->  kEventVertexZLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutSharedCls:   // 9 -->  kSharedClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedClsLoose, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsTight:   // 10 -->  kFindableClsTight
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsTight,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsLoose:   // 11 -->  kFindableClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsLoose,kTPCSignalN};
    }
    break;
    //
    case kCutPileupLoose:   // 12 -->  kPileupLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileupLoose, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutBFieldPos:   // 13 -->  kBFieldPos
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldPos};
    }
    break;
    //
    case kCutBFieldNeg:   // 14 --> kBFieldNeg
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldNeg};
    }
    break;
    //
    case kCutTPCSignalNSmall:   // 15 --> kTPCSignalNSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNSmall};
    }
    break;
    //
    case kCutTPCSignalNLarge:   // 16 --> kTPCSignalNLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNLarge};
    }
    break;
    //
    default:
    {
      fCutArr = {};
    }

  }
  //
  //  Apply conditions
  for (UInt_t i=0;i<fCutArr.size();i++){
    if( ((cut >> fCutArr[i]) & 1) == 0 )
    {
      return kFALSE;
    }
  }
  return kTRUE;

}
//________________________________________________________________________
Int_t AliAnalysisJetHadro::CountEmptyEvents(Int_t counterBin)
{

  //
  // count Empty Events
  //
  Int_t emptyCount=0;
  for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack) {   // Track loop
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track->GetInnerParam()) continue;
    Float_t momtrack = track->GetInnerParam()->GetP();
    if (momtrack<fMomDown || momtrack>fMomUp) continue;
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    if (track->GetTPCsignalN()<60) continue;
    if (track->GetTPCsignal()>0) emptyCount++;
  }
  //
  // check if the event is empty
  if (emptyCount<1) {
    //fHistEmptyEvent->Fill(counterBin);
    std::cout << " Info::siweyhmi: Empty event in " << fChunkName << std::endl;
  }
  if (fUseCouts) std::cout << " Info::siweyhmi: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  return emptyCount;

}
//________________________________________________________________________
void AliAnalysisJetHadro::PrintNumInBinary(UInt_t num)
{
  TString bin="";
  Int_t numberOfBits = sizeof(UInt_t)*8;
  for (Int_t i=numberOfBits-1; i>=0; i--) {
    Bool_t isBitSet = (num & (1<<i));
    if (isBitSet) {
      bin+="1";
    } else {
      bin+="0";
    }
  }
  std::cout << "Info::siweyhmi: fTrackCutBits = " << bin << std::endl;
}
//________________________________________________________________________
void AliAnalysisJetHadro::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  std::cout << " Info::siweyhmi: ===== In the Terminate ===== " << std::endl;
/*
  for (int i=0; i<fMomExpec_NBins; i++){
    for (int j=0; j<fEtaExpec_NBins; j++){
    TH1D *test_mpi = fHistIncTracks_mpi->ProjectionY("test_mpi",i+1,i+1,j+1,j+1,"e");
    Double_t binx = fHistIncTracks_mpi_small->GetXaxis()->GetBinCenter(i+1);
    Double_t biny = fHistIncTracks_mpi_small->GetYaxis()->GetBinCenter(j+1);
    fHistIncTracks_mpi_small->Fill(binx,biny,test_mpi->GetMean());
    TH1D *test_spi = fHistIncTracks_spi->ProjectionY("test_spi",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_spi_small->Fill(binx,biny,test_spi->GetMean());

    TH1D *test_mel = fHistIncTracks_mel->ProjectionY("test_mel",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_mel_small->Fill(binx,biny,test_mel->GetMean());
    TH1D *test_sel = fHistIncTracks_sel->ProjectionY("test_sel",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_sel_small->Fill(binx,biny,test_sel->GetMean());

    TH1D *test_mka = fHistIncTracks_mka->ProjectionY("test_mka",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_mka_small->Fill(binx,biny,test_mka->GetMean());
    TH1D *test_ska = fHistIncTracks_ska->ProjectionY("test_ska",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_ska_small->Fill(binx,biny,test_ska->GetMean());

    TH1D *test_mpr = fHistIncTracks_mpr->ProjectionY("test_mpr",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_mpr_small->Fill(binx,biny,test_mpr->GetMean());
    TH1D *test_spr = fHistIncTracks_spr->ProjectionY("test_spr",i+1,i+1,j+1,j+1,"e");
    fHistIncTracks_spr_small->Fill(binx,biny,test_spr->GetMean());
    }
  }
*/

  std::cout << "The totalJetArea of all jets that pass our cuts is "   << ftotalJetArea << std::endl;
  std::cout << "The totalNumRealJets that pass our cuts is "   << ftotalNumRealJets << std::endl;
  std::cout << "The totalNumIncEvents is "   << ftotalNumIncEvents << std::endl;
  std::cout << "The totalNumRealJetEvents is "   << ftotalNumRealJetEvents << std::endl;

}
