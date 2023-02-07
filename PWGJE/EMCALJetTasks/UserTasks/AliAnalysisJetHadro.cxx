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
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TRandom.h"
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
#include "AliESDv0KineCuts.h"
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
#include "AliESDv0Cuts.h"
#include "AliESDv0.h"
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
#include "AliAnalysisTaskFilteredTree.h"
#include "AliAnalysisJetHadro.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliRunLoader.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliESDtools.h"
#include "AliFJWrapper.h"
#include "AliEmcalJet.h"
#include "AliEmcalJetTask.h"
#include "AliRhoParameter.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
using namespace std;
using std::cout;
using std::setw;

ClassImp(AliAnalysisJetHadro)

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro()
: AliAnalysisTaskEmcalJet("TaskEbyeRatios"), fEventCuts(0), fPIDResponse(0),fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit96_spd(0),
fESDtrackCuts_Bit96_sdd(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreeMC(0x0),
fTreejetsFJ(0x0),
fTreeBGjetsFJ(0x0),
fTreeCuts(0x0),
fTreejetsFJconst(0x0),
fTreeResonance(0x0),
fTreeEvents(0x0),
fTreeDScaled(0x0),
fTreejetsFJGen(0x0),
fTreejetEMC(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHndEdx(),
fHnExpected(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracks(kFALSE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fFillBGJetsFJTree(kTRUE),
fFilldscaledTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fNSettings(22),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
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
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(10),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
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
fHasV0FirstITSlayer(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0),
fEventInfo_LumiGraph(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  // default Constructor
  /* fast compilation test
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisJetHadro.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisJetHadro.cxx++
  */
}

//________________________________________________________________________
AliAnalysisJetHadro::AliAnalysisJetHadro(const char *name)
: AliAnalysisTaskEmcalJet(name), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit96_spd(0),
fESDtrackCuts_Bit96_sdd(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreeMC(0x0),
fTreejetsFJ(0x0),
fTreeBGjetsFJ(0x0),
fTreeCuts(0x0),
fTreejetsFJconst(0x0),
fTreeResonance(0x0),
fTreeEvents(0x0),
fTreeDScaled(0x0),
fTreejetsFJGen(0x0),
fTreejetEMC(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHndEdx(),
fHnExpected(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaDown(0),
fEtaUp(0),
fNEtaBins(0),
fPercentageOfEvents(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fWeakAndMaterial(kFALSE),
fEffMatrix(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracks(kFALSE),
fFillOnlyHists(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fFillBGJetsFJTree(kTRUE),
fFilldscaledTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fFillDnchDeta(kFALSE),
fIncludeTOF(kFALSE),
fUseThnSparse(kFALSE),
fUseCouts(kFALSE),
fNSettings(22),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
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
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fNResModeMC(2),
fNCentbinsData(10),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
frhoFJ(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
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
fHasV0FirstITSlayer(0),
fSystCrossedRows(0),
fSystDCAxy(0),
fSystChi2(0),
fSystVz(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistdEdxTPC(0),
fHistArmPod(0),
fEventInfo_LumiGraph(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  //
  //         standard constructur which should be used
  //
  /* fast compilation test
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/sw/ubuntu1404_x86-64/AliPhysics/master-1/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/AliPhysics/PWGPP");
  .L /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AliAnalysisJetHadro.cxx++
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisJetHadro.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisJetHadro.cxx++
  */
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:***************** CONSTRUCTOR CALLED: AliAnalysisJetHadro  *****************"<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
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

}
//________________________________________________________________________
AliAnalysisJetHadro::~AliAnalysisJetHadro()
{

  //
  // Destructor
  //
  std::cout << " Info::siweyhmi: ===== In the Destructor ===== " << std::endl;
  if (fHistPosEffMatrixRec) delete fHistPosEffMatrixRec;
  if (fHistNegEffMatrixRec) delete fHistNegEffMatrixRec;
  if (fHistPosEffMatrixGen) delete fHistPosEffMatrixGen;
  if (fHistNegEffMatrixGen) delete fHistNegEffMatrixGen;
  if (fHistPosEffMatrixScanRec) delete fHistPosEffMatrixScanRec;
  if (fHistNegEffMatrixScanRec) delete fHistNegEffMatrixScanRec;
  if (fHistPosEffMatrixScanGen) delete fHistPosEffMatrixScanGen;
  if (fHistNegEffMatrixScanGen) delete fHistNegEffMatrixScanGen;
  if (fHistdEdxTPC)         delete fHistdEdxTPC;
  if (fHistEmptyEvent)      delete fHistEmptyEvent;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistCentralityImpPar)delete fHistCentralityImpPar;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHistArmPod)          delete fHistArmPod;
  if (fHistCent)            delete fHistCent;
  if (fHistPhi)             delete fHistPhi;
  if (fHistGenMult)         delete fHistGenMult;
  if (fHistRapDistFullAccPr)   delete fHistRapDistFullAccPr;
  if (fHistRapDistFullAccAPr)  delete fHistRapDistFullAccAPr;
  if (fHistInvK0s)          delete fHistInvK0s;
  if (fHistInvLambda)       delete fHistInvLambda;
  if (fHistInvAntiLambda)   delete fHistInvAntiLambda;
  if (fHistInvPhoton)       delete fHistInvPhoton;
  if (fHndEdx)              delete fHndEdx;
  if (fHnExpected)          delete fHnExpected;
  if (fPIDCombined)         delete fPIDCombined;
  if (fESDtrackCuts)          delete fESDtrackCuts;
  if (fESDtrackCuts_Bit96)    delete fESDtrackCuts_Bit96;
  if (fESDtrackCuts_Bit96_spd)    delete fESDtrackCuts_Bit96_spd;
  if (fESDtrackCuts_Bit96_sdd)    delete fESDtrackCuts_Bit96_sdd;
  if (fESDtrackCuts_Bit128)   delete fESDtrackCuts_Bit128;
  if (fESDtrackCuts_Bit768)   delete fESDtrackCuts_Bit768;
  if (fESDtrackCutsLoose)     delete fESDtrackCutsLoose;
  if (fESDtrackCutsV0)        delete fESDtrackCutsV0;
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
  if (fRunFastHighMomentCal) { std::cout << " Info::siweyhmi: !!! We are running fast high moment calculation return !!! " << std::endl; return; }
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  //
  // ------------------------------------------------
  //
  // tight DCA cut used by Emil
  fESDtrackCuts_Bit96_spd = new AliESDtrackCuts("Bit96_spd",""); fESDtrackCuts_Bit96_spd->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  fESDtrackCuts_Bit96_sdd = new AliESDtrackCuts("Bit96_sdd",""); fESDtrackCuts_Bit96_sdd->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
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
  fESDtrackCuts_Bit768->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit768->SetPtRange(0.,100000.);
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
  fESDtrackCuts->SetEtaRange(-100.,100.);
  fESDtrackCuts->SetPtRange(0.,100000.);
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
  fESDtrackCutsLoose->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsLoose->SetMaxFractionSharedTPCClusters(0.4);
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
  // V0 selection
  // ------------------------------------------------
  fESDtrackCutsV0 = new AliESDv0Cuts("AliESDCutsV0","");
  fESDtrackCutsV0 ->SetMaxDcaV0Daughters(1.0);
  // ------------------------------------------------
  //
  // Special selection for clean samples
  fV0OpenCuts   = new AliESDv0KineCuts();
  fV0StrongCuts = new AliESDv0KineCuts();
  SetSpecialV0Cuts(fV0OpenCuts);
  SetSpecialV0Cuts(fV0StrongCuts);

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
  if (!(fRunFastSimulation || fRunFastHighMomentCal)) {
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
  //   histogram of splines
  // ************************************************************************
  //
  if(!fEffMatrix && !fMCtrue)
  {
    // 0 --> assumed particle: 0. electron, 1. pion, 2. kaon, 3. proton, 5. deuteron
    // 1 --> sign
    // 2 --> Centrality
    // 3 --> eta
    // 4 --> momentum
    // 5 --> Expected Sigma of a given track
    // 6 --> Expected Mean of a given track
    const Int_t nExpectedbins = 7;
    //                                        0    1,    2,                 3,           4,        5,       6
    Int_t   binsExpected[nExpectedbins]  = {  5,   3,  fNCentbinsData,   fNEtaBins,   fNMomBins,   240,   4000 };  // ????
    Double_t xminExpected[nExpectedbins] = {  0., -2.,   0.,             fEtaDown,    fMomDown,     1.,   20.  };
    Double_t xmaxExpected[nExpectedbins] = {  5.,  2.,  80.,             fEtaUp,      fMomUp,      61.,   1020.};
    TString axisNameExpected[nExpectedbins]   = {"particleType","sign","Centrality"    ,"eta" ,"momentum" ,"ExSigma","ExMean"};
    TString axisTitleExpected[nExpectedbins]  = {"particleType","sign","Centrality [%]","#eta","#it{p} (GeV/#it{c})", "#sigma","#mu"};
    fHnExpected = new THnSparseF("hExpected","hExpected",nExpectedbins,binsExpected,xminExpected,xmaxExpected);
    std::cout << "Number of centrality bins: " << fxCentBins.size() << std::endl;
    fHnExpected->GetAxis(2)->Set(fNCentbinsData-1,fxCentBins.data());
    for (Int_t iaxis=0; iaxis<nExpectedbins;iaxis++){
      fHnExpected->GetAxis(iaxis)->SetName(axisNameExpected[iaxis]);
      fHnExpected->GetAxis(iaxis)->SetTitle(axisTitleExpected[iaxis]);
    }
    fListHist->Add(fHnExpected);

  }
  //
  // ************************************************************************
  //   dEdx histograms
  // ************************************************************************
  //
  // 0 --> sign
  // 1 --> Centrality
  // 2 --> eta
  // 3 --> momentum
  // 4 --> TPC dEdx
  //                                  0,         1,           2,           3,            4
  const Int_t nhistbins = 5;
  Int_t dEdxnBins = Int_t((fDEdxUp-fDEdxDown)/fDEdxBinWidth);
  TString axisNamedEdx[nhistbins]   = {"Sign" ,"Centrality"     ,"eta"  ,"momentum"  ,"TPC dEdx"};
  TString axisTitledEdx[nhistbins]  = {"Sign" ,"Centrality [%]" ,"#eta" ,"#it{p} (GeV/#it{c})" ,"TPC d#it{E}/d#it{x} Signal (a.u.)"};
  if (fUseThnSparse)
  {
    //
    // inclusive spectra
    Int_t   binsdEdx[nhistbins]  = { 2,  fNCentbinsData, fNEtaBins,   fNMomBins,    dEdxnBins};
    Double_t xmindEdx[nhistbins] = {-2,  0.,             fEtaDown,    fMomDown,     fDEdxDown};
    Double_t xmaxdEdx[nhistbins] = { 2,  80.,            fEtaUp,      fMomUp,       fDEdxUp};
    fHndEdx= new THnSparseF("hdEdx","Inclusive dEdx Spectrum"  ,nhistbins,binsdEdx,xmindEdx,xmaxdEdx);
    fHndEdx->GetAxis(1)->Set(fNCentbinsData-1,fxCentBins.data());
    // Set the branch names
    for (Int_t iaxis=0; iaxis<nhistbins;iaxis++){
      fHndEdx   ->GetAxis(iaxis)->SetName(axisNamedEdx[iaxis]);
      fHndEdx   ->GetAxis(iaxis)->SetTitle(axisTitledEdx[iaxis]);
    }
    fListHist->Add(fHndEdx);
  }
  //
  // ************************************************************************
  //   Efficiency matrix histograms
  // ************************************************************************
  //
  if(fEffMatrix && !fRunOnGrid)
  {
    const Int_t ndim=5;
    Int_t nbins0[ndim]  ={3,8, 40       ,16        ,50  };
    Double_t xmin0[ndim]={0,0, fMomDown ,fEtaDown  ,0.  };
    Double_t xmax0[ndim]={3,80,fMomUp   ,fEtaUp    ,6.25};
    fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec","fHistPosEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec","fHistNegEffMatrixRec",ndim, nbins0,xmin0,xmax0);
    fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen","fHistPosEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen","fHistNegEffMatrixGen",ndim, nbins0,xmin0,xmax0);
    TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
    TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
    for (Int_t iEff=0; iEff<ndim;iEff++){
      fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
      fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    }
    fListHist->Add(fHistPosEffMatrixRec);
    fListHist->Add(fHistNegEffMatrixRec);
    fListHist->Add(fHistPosEffMatrixGen);
    fListHist->Add(fHistNegEffMatrixGen);
    //
    //
    const Int_t ndimScan=6;
    Int_t nbinsScan[ndimScan]   = {2, fNSettings,           3,8, 40       ,16   };
    Double_t xminScan[ndimScan] = {0, 0,                    0,0, fMomDown ,fEtaDown };
    Double_t xmaxScan[ndimScan] = {2, Double_t(fNSettings), 3,80,fMomUp   ,fEtaUp };
    fHistPosEffMatrixScanRec  =new THnF("fHistPosEffMatrixScanRec","fHistPosEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistNegEffMatrixScanRec  =new THnF("fHistNegEffMatrixScanRec","fHistNegEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistPosEffMatrixScanGen  =new THnF("fHistPosEffMatrixScanGen","fHistPosEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
    fHistNegEffMatrixScanGen  =new THnF("fHistNegEffMatrixScanGen","fHistNegEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
    TString axisNameEffScan[ndimScan]  = {"detector", "systematic", "particle"      ,"Centrality"     ,"momentum"      ,"eta"  };
    TString axisTitleEffScan[ndimScan] = {"detector", "setting",    "particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta"};
    for (Int_t iEff=0; iEff<ndimScan;iEff++){
      fHistPosEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistNegEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistPosEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
      fHistNegEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    }
    fListHist->Add(fHistPosEffMatrixScanRec);
    fListHist->Add(fHistNegEffMatrixScanRec);
    fListHist->Add(fHistPosEffMatrixScanGen);
    fListHist->Add(fHistNegEffMatrixScanGen);
  }
  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  Int_t nSafeBins = (fRunFastSimulation) ? 1200 : 5;
  fHistEmptyEvent        = new TH1F("hEmptyEvent",           "control histogram to count empty events"    , 10,  0., 10.);
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 10,  0., 100.);
  fHistCentralityImpPar  = new TH1F("hCentralityImpPar",     "control histogram for centrality imppar"    , 10,  0., 100.);
  fHistImpParam          = new TH1F("hImpParam",             "control histogram for impact parameter"     , 200, 0., 20.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);
  fHistGenMult           = new TH1F("hGenPrMult",            "generated protons"                          , fGenprotonBins,0., 200.);
  fHistRapDistFullAccPr  = new TH2F("hRapDistFullAccPr",     "rapidity dist of protons in full acceptance"    , nSafeBins, -15, 15., 10, 0., 100.);
  fHistRapDistFullAccAPr = new TH2F("hRapDistFullAccApr",    "rapidity dist of antiprotons in full acceptance", nSafeBins, -15, 15., 10, 0., 100.);
  fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  fListHist->Add(fHistCentralityImpPar);
  fListHist->Add(fHistImpParam);
  fListHist->Add(fHistVertex);
  fListHist->Add(fHistGenMult);
  fListHist->Add(fHistRapDistFullAccPr);
  fListHist->Add(fHistRapDistFullAccAPr);
  //
  // ************************************************************************
  //   Clean sample helper histograms
  // ************************************************************************
  //
  if (fUseCouts && fFillArmPodTree)
  {
    fHistArmPod            = new TH2F("hArmPod",           "Armenteros-Podolanski plot"                     , 100,-1.,1., 110,0.,0.22);
    fHistInvK0s            = new TH1F("fHistInvK0s",       "control histogram for K0s invariant mass"       , 1000, 0.3,  0.70);
    fHistInvLambda         = new TH1F("fHistInvLambda",    "control histogram for lambda invariant mass"    , 1000, 1.07, 1.16);
    fHistInvAntiLambda     = new TH1F("fHistInvAntiLambda","control histogram for antilambda invariant mass", 1000, 1.07, 1.16);
    fHistInvPhoton         = new TH1F("fHistInvPhoton",    "control histogram for photon invariant mass"    , 1000, 0.,   0.05);
    fListHist->Add(fHistInvK0s);
    fListHist->Add(fHistInvLambda);
    fListHist->Add(fHistInvAntiLambda);
    fListHist->Add(fHistInvPhoton);
    fListHist->Add(fHistArmPod);
  }
  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fArmPodTree    = ((*fTreeSRedirector)<<"fArmPodTree").GetTree();
  fTreejetsEMCconst    = ((*fTreeSRedirector)<<"jetsEMCconst").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreejetsFJ    = ((*fTreeSRedirector)<<"jetsFJ").GetTree();
  fTreeBGjetsFJ  = ((*fTreeSRedirector)<<"BGjetsFJ").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"tracks").GetTree();
  fTreejetEMC       = ((*fTreeSRedirector)<<"jetEMC").GetTree();
  fTreejetsFJconst = ((*fTreeSRedirector)<<"jetsFJconst").GetTree();
  fTreeResonance = ((*fTreeSRedirector)<<"resonance").GetTree();
  fTreeEvents    = ((*fTreeSRedirector)<<"eventInfo").GetTree();
  fTreeDScaled   = ((*fTreeSRedirector)<<"dscaled").GetTree();
  fTreejetsFJGen = ((*fTreeSRedirector)<<"jetsFJGen").GetTree();
  //
  // ************************************************************************
  //   Send output objects to container
  // ************************************************************************
  //
  PostData(1, fListHist);
  PostData(2, fArmPodTree);
  PostData(3, fTreejetsEMCconst);
  PostData(4, fTreeMC);
  PostData(5, fTreejetsFJ);
  PostData(6, fTreeBGjetsFJ);
  PostData(7, fTreeCuts);
  PostData(8, fTreejetsFJconst);
  PostData(9, fTreeResonance);
  PostData(10, fTreeEvents);
  PostData(11, fTreeDScaled);
  PostData(12, fTreejetsFJGen);
  PostData(13, fTreejetEMC);

  std::cout << " Info::siweyhmi: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
Bool_t AliAnalysisJetHadro::Run()
{
  //
  // main event loop
  //
  if (fRunOnGrid) fUseCouts=kFALSE; // for security
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
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,0); // do not apply any pile cut
        // if (!fEventCuts.AcceptEvent(fESD)) {cout<< "pileup event " << endl; return;}
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
  if (!fPIDResponse && !(fRunFastSimulation || fRunFastHighMomentCal)) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
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
        fHistCentralityImpPar->Fill(fCentImpBin);
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
      fHistCentralityImpPar->Fill(fCentImpBin);
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
  if (!(fRunFastSimulation || fRunFastHighMomentCal))
  {
    //
    // ========================== Real =========================
    //
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
      if ( vertexType.Contains("vertexer: Z") && (fVertex->GetDispersion() > 0.04 || fVertex->GetZRes() > 0.25) ) isVertexOk = kFALSE; // TODO
    }
    fMultiplicity    = fVertex->GetNContributors();    // fMultiplicity = fESD -> GetNumberOfTracks();
    fNContributors   = fVertex->GetNContributors();
    fMultiplicityMC  = fMultiplicity;
    //
    // ------------------------------------------------
    // ------- event vertex cut along Z ---------------
    // ------------------------------------------------
    //
    // if (fMCtrue && TMath::Abs(fVz) > 15) return;   // For MC put fixed cut
    if (fDefaultTrackCuts && (TMath::Abs(fVz)>7 || TMath::Abs(fVz)<0.15) ) return kFALSE;
    else if (TMath::Abs(fVz)>15) return kFALSE;
    //
    if (fVertex && isVertexOk) fHistVertex->Fill(fVz);
    else return kFALSE;
    //
    // ------------------------------------------------
    // ---------- Centrality definition ---------------
    // ------------------------------------------------
    //
    if (fBeamType.CompareTo("A-A") == 0) { // PbPb
      if (MultSelection) {
        fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      } else if (esdCentrality) {
        fCentrality = esdCentrality->GetCentralityPercentile("V0M");
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
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)!=0) return kFALSE;
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fFillTracks && fESD){
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = -2.0;
    fjetRhoVal = -2.0;

    FillTPCdEdxReal();
    FindJetsEMC();
    FindJetsFJ();
    if (fFillArmPodTree) FillCleanSamples();
    FillEventTree(); //rhoEMC, rhoFJ, cent, yesjet, yesgoodjet
    if (fUseCouts)  std::cout << " Info::siweyhmi: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  //
  // full MC analysis
  //
  if (fMCtrue && fESD){
    FillMCFull();
    FindJetsFJGen();
    FindJetsEMC();
    FillEffMatrix();
    FindJetsFJ();
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
  Double_t pT_sub_min = 10.0;
  // Get the jet container
  fJetContainer = this->GetJetContainer("detJets");
  TString fRhoName = fJetContainer->GetRhoName();
  if (fUseCouts) cout << "Rho Name is " << fRhoName << endl;

  fjetRhoVal = -2.0;
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
    if (TMath::Abs(jet->Eta()) >= 0.5) continue;
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
    Double_t jetPtSub = jet->PtSub(fjetRhoVal, kFALSE); //bg sub pt

    if (jetPtSub > pT_sub_min)
    {
      fhasRealEMCjet = 1;
    }

    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"jetEMC"<<
    "gid="       << fEventGID << //  global event ID
    "jetRadius=" << jetRadius << // jet Radius
    "Njets=" << Njets << // Number of jets in event
    "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
    "jetRhoVal=" << fjetRhoVal <<
    "jetEtaMin=" << jetEtaMin << //min eta cut for jet
    "jetEtaMax=" << jetEtaMax << //max eta cut for jet
    "jetAcceptanceType=" << jetAcceptanceType << //what detector (and fidicuial = account for R) acceptance used
    "jetpt="     << fJetPt    << // jetPt
    "jeteta="    << fJetEta   << // jetEta
    "jetphi="    << fJetPhi   << // jetEta
    "jetM="    << JetM   <<
    "jetLabel="    << JetLabel   <<
    "jetArea="    << JetAreaPt  <<
    "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
    "IsJetMc="    << IsJetMc  <<
    "particleMaxChargedPt="  << particleMaxChargedPt <<
    "jetMCPt=" << jetMCPt <<
    "jetPtSub=" << jetPtSub <<
    "cent="      << fCentrality  <<  //  centrality
    "\n";

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      if (TMath::Abs(jet->Eta()) >= 0.5) continue;
      const AliVParticle* particle = jet->Track(i);
      AliESDtrack* esdtrack = (AliESDtrack*)(particle);
      if (!esdtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;

      fTrackCutBits=0;  // reset the bits for the next track

      //
      // --------------------------------------------------------------
      //      Get relevant track info and set cut bits
      // --------------------------------------------------------------
      //
      Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(esdtrack);
      Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(esdtrack);
      if (!esdtrack->GetInnerParam()) continue;               // Ask if track is in the TPC
      if (!fESDtrackCutsLoose->AcceptTrack(esdtrack))  continue;    // Loose cuts
      if (!(esdtrack->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      Float_t closestPar[3];
      GetExpecteds(esdtrack,closestPar);
      SetCutBitsAndSomeTrackVariables(esdtrack,0);
      Int_t tpcNcls = esdtrack->GetTPCncls();
      Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
      AliVMultiplicity *multiObj = fESD->GetMultiplicity();
      Int_t nITSClusters = 0;
      for(Int_t j=2;j<6;j++) nITSClusters += multiObj->GetNumberOfITSClusters(j);
      //
      // Fill track constituent information
      (*fTreeSRedirector)<<"jetsEMCconst"<<
      "gid="       << fEventGID << //  global event ID
      "jetRadius=" << jetRadius << // jet Radius
      "Njets=" << Njets << // Number of jets in event
      "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
      "jetRhoVal=" << fjetRhoVal <<
      "jetEtaMin=" << jetEtaMin << //min eta cut for jet
      "jetEtaMax=" << jetEtaMax << //max eta cut for jet
      "jetAcceptanceType=" << jetAcceptanceType << //what detector (and fidicuial = account for R) acceptance used
      "jetpt="     << fJetPt    << // jetPt
      "jeteta="    << fJetEta   << // jetEta
      "jetphi="    << fJetPhi   << // jetEta
      "jetM="    << JetM   <<
      "jetLabel="    << JetLabel   <<
      "jetArea="    << JetAreaPt  <<
      "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
      "IsJetMc="    << IsJetMc  <<
      "particleMaxChargedPt="  << particleMaxChargedPt <<
      "jetMCPt=" << jetMCPt <<
      "jetPtSub=" << jetPtSub <<

      //
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit768="    << fBit768 <<        // Hybrid track cuts
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
      "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "cent="      << fCentrality           <<  //  centrality
      //
      "phi="       << fPhi                  <<  //  phi
      "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
      "dcaz="      << fTrackDCAz            <<  // dca cut along z
      "ncltpc="    << fNcl                  <<  // number of clusters
      "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
      "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
      "missCl="    << fMissingCl            <<  // fraction of missing clusters
      //
      "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
      "pitpcpid="  << fNSigmasPiTPC         <<
      "katpcpid="  << fNSigmasKaTPC         <<
      "prtpcpid="  << fNSigmasPrTPC         <<
      "detpcpid="  << fNSigmasDeTPC         <<
      //
      "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
      "pitofpid="  << fNSigmasPiTOF         <<
      "katofpid="  << fNSigmasKaTOF         <<
      "prtofpid="  << fNSigmasPrTOF         <<
      "detofpid="  << fNSigmasDeTOF         <<
      //
      "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi="  << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa="  << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr="  << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "dEdxMeanDe="  << fDEdxDe              <<
      "dEdxSigmaDe=" << fSigmaDe            <<
      "closestTPCPIDtype=" << closestPar[1]         << //particle type
      "closestTPCPIDmass=" << closestPar[2]         << //particle mass
      "\n";
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
  Float_t jetRadius = 0.4;
  Float_t bgJetRadius = 0.2;
  Float_t jetAbsEtaCut = 0.5;
  Float_t bgJetAbsEtaCut = 0.7;
  Float_t pT_sub_min = 10.0;

  //SOME CODE FROM NIMA
  AliFJWrapper *fFastJetWrapper;
  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
  fFastJetWrapper->Clear();
  fFastJetWrapper->SetR(jetRadius);
  fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
  fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme); // fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
  fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
  fFastJetWrapper->SetGhostArea(0.005);
  fFastJetWrapper->SetAreaType(fastjet::AreaType::passive_area);
  std::vector<int> trackTTIndex;
  trackTTIndex.clear();
  //
  // Access and loop over container
  // auto tracks = this->GetTrackContainer("detTracks");
  // for (auto esdtrack : tracks->accepted()){...}
  //
  std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
  std::vector<fastjet::PseudoJet> particlesEmbedded; //fill this with your event
  double particleEtaCut = 0.9;
  //
  // loop over esd tracks and add their four vector to wrapper
  for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = fESD->GetTrack(iTrack);
    if (!track->GetInnerParam()) continue;
    if (!fESDtrackCuts->AcceptTrack(track)) continue;
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    if (track->Pt() < 0.15 || TMath::Abs(track->Eta()) >= 0.9) continue;
    fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
    particlesEmbedded.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), TMath::Sqrt(track->P()*track->P()+0.13957*0.13957) ) ); //filling with event
  }

  fastjet::JetMedianBackgroundEstimator bgE;
  fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2) * fastjet::SelectorAbsEtaMax(particleEtaCut); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
  //bgE.set_selector(selectorBG);
  fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, bgJetRadius, fastjet::E_scheme, fastjet::Best); //define the kT jet finding which will do the average background estimation
  fastjet::GhostedAreaSpec ghostSpecBG(particleEtaCut, 1, 0.005); //this ghost area might be too small and increase processing time too much
  fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts, ghostSpecBG);
  fastjet::ClusterSequenceArea cluster_seq_BG(particlesEmbedded, jetDefBG, areaDefBG);
  std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
  bgE.set_jets(jetsBG);  // give the kT jets to the background estimator

  frhoFJ = bgE.rho();

  for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
    fastjet::PseudoJet jet = jetsBG[ijet];
    if (jet.pt() < 0.15 || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= bgJetAbsEtaCut) continue;
    Float_t jetpt = jet.pt();
    Float_t jetphi = jet.phi();
    Float_t jeteta = jet.eta();
    Float_t jetArea = jet.area();
    Float_t jetptsub = jetpt - frhoFJ*jetArea;
    Int_t jetNum = jetsBG.size();
    Int_t Njets = jetsBG.size();

    if (fFillBGJetsFJTree)
    {
      (*fTreeSRedirector)<<"BGjetsFJ"<<
      "gid="       << fEventGID << //  global event ID
      "bjJetRadius=" << bgJetRadius << // jet Radius
      "bgJetAbsEtaCut=" << bgJetAbsEtaCut << //abs eta cut for jet
      "jetNum="    << jetNum <<    //  number of jets
      "jetpt="     << jetpt <<     //  global event ID
      "jetphi="    << jetphi <<    //  global event ID
      "jeteta="    << jeteta <<    //  global event ID
      "jetptsub=" << jetptsub << //bg sub jet pt (pt - rho*Area)
      "rhoFJ="      << frhoFJ << //event rho
      "jetArea="    << jetArea << //jet area
      "cent="      << fCentrality  <<  //  centrality
      "\n";
    }
  }

  fastjet::contrib::ConstituentSubtractor subtractorConstituent(&bgE); //add the background estimator to the correct subtractor
  subtractorConstituent.set_common_bge_for_rho_and_rhom(true); //CHECK : should not be the case since particles have mass
  subtractorConstituent.set_max_standardDeltaR(0.25); // set the max event wise subtraction distance
  particlesEmbeddedSubtracted = subtractorConstituent.subtract_event(particlesEmbedded, particleEtaCut); //perform subtraction and fill the subtracted event container

  fFastJetWrapper->Run();
  std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();

  for (Int_t ijet=0; ijet<Int_t(jets.size()); ijet++){
    fastjet::PseudoJet jet = jets[ijet];
    if (jet.pt() < 0.15 || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
    fhasAcceptedFJjet = 1;
    Float_t jetpt = jet.pt();
    Float_t jetphi = jet.phi();
    Float_t jeteta = jet.eta();
    Float_t jetArea = jet.area();
    Float_t jetptsub = jetpt - frhoFJ*jetArea;
    Int_t jetNum = jets.size();

    if (jetptsub > pT_sub_min)
    {
      fhasRealFJjet = 1;
    }

    //
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
    Int_t nConstituents = constituents.size();

    (*fTreeSRedirector)<<"jetsFJ"<<
    "gid="       << fEventGID << //  global event ID
    "jetRadius=" << jetRadius << // jet Radius
    "jetAbsEtaCut=" << jetAbsEtaCut << //abs eta cut for jet
    "jetNum="    << jetNum <<    //  number of jets
    "jetpt="     << jetpt <<     //  global event ID
    "jetphi="    << jetphi <<    //  global event ID
    "jeteta="    << jeteta <<    //  global event ID
    "nConst="    << nConstituents <<    //  global event ID
    "cent="      << fCentrality           <<  //  centrality

    "jetptsub="  << jetptsub << //bg sub jet pt (pt - rho*Area)
    "rhoFJ="       << frhoFJ << //event rho
    "jetArea="   << jetArea << //jet area
    "\n";

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
      Bool_t fBit96_spd    = fESDtrackCuts_Bit96_spd->AcceptTrack(trackConst);
      Bool_t fBit96_sdd    = fESDtrackCuts_Bit96_sdd->AcceptTrack(trackConst);
      Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(trackConst);
      Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackConst);
      Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackConst);
      Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackConst);
      if (!trackConst->GetInnerParam()) continue;               // Ask if track is in the TPC
      if (!fESDtrackCutsLoose->AcceptTrack(trackConst))  continue;    // Loose cuts
      if (!(trackConst->GetTPCsignalN()>0)) continue;
      //
      // Get the track variables
      Float_t closestPar[3];
      GetExpecteds(trackConst,closestPar);
      SetCutBitsAndSomeTrackVariables(trackConst,0);
      Int_t tpcNcls = trackConst->GetTPCncls();
      Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
      Int_t nITSClusters = 0;
      AliVMultiplicity *multiObj = fESD->GetMultiplicity();
      for(Int_t j=2;j<6;j++) nITSClusters += multiObj->GetNumberOfITSClusters(j);
      //
      // different dca cuts
      // TMath::Abs(fTrackDCAxy)< 0.3
      Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
      Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);   // 10h tuned loose cut
      Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);     // 10h tuned loose cut


      (*fTreeSRedirector)<<"jetsFJconst"<<
      "gid="       << fEventGID << //  global event ID
      "jetRadius=" << jetRadius << // jet Radius
      "jetAbsEtaCut=" << jetAbsEtaCut << //abs eta cut for jet
      "jetNum="    << jetNum <<    //  number of jets
      "jetpt="     << jetpt <<     //  global event ID
      "jetphi="    << jetphi <<    //  global event ID
      "jeteta="    << jeteta <<    //  global event ID
      "nConst="    << nConstituents <<    //  global event ID

      "jetptsub="  << jetptsub << //bg sub jet pt (pt - rho*Area)
      "rhoFJ="       << frhoFJ << //event rho
      "jetArea="   << jetArea << //jet area

      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit96spd="  << fBit96_spd <<     // cut account for the regions where no ITS-TPC match
      "bit96sdd="  << fBit96_sdd <<     // cut account for the regions where no ITS-TPC match
      "bit96="     << fBit96_base <<    // tight cuts of 2011 tuned data
      "bit128="    << fBit128 <<        // TPC only tracks cuts
      "bit768="    << fBit768 <<        // Hybrid track cuts
      "pixCut="    << ifDCAcutIfNoITSPixel <<    // cut: apply a DCAcut If No ITS Pixel
      "run="       << fRunNo <<                  // run Number
      "bField="    << fBField <<                 // magnetic filed
      "pileupbit=" << fPileUpBit <<              // flag for pileup selection
      "primMult="  << fNContributors <<          //  #prim tracks
      //
      "dcabase="  << dcaBaseCut <<  //  TPC multiplicity
      "dca10h="   << dca10h <<  //  TPC multiplicity
      "dca11h="   << dca11h <<  //  TPC multiplicity
      //
      "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
      "itsclmult=" << nITSClusters <<    // ITS multiplicity
      "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
      //
      "eventtime=" << fTimeStamp            <<  // event timeStamp
      "intrate="   << fIntRate              <<  // interaction rate
      "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "cent="      << fCentrality           <<  //  centrality
      //
      "phi="       << fPhi                  <<  //  phi
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

      "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi=" << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa=" << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr=" << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "dEdxMeanDe=" << fDEdxDe              <<
      "dEdxSigmaDe=" << fSigmaDe            <<

      "closestTPCPIDtype=" << closestPar[1]         << //particle type
      "closestTPCPIDmass=" << closestPar[2]         << //particle mass
      "\n";
    }

  }

}

//________________________________________________________________________
void AliAnalysisJetHadro::FillEventTree()
{
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEventTree ===== " << std::endl;
  (*fTreeSRedirector)<<"eventInfo"<<
  "gid="       << fEventGID << //  global event ID
  "rhoFJ="      << frhoFJ << //event rho
  "rhoEMC="  << fjetRhoVal <<
  "cent="      << fCentrality  <<  //  centrality
  "hasAcceptedFJjet="   << fhasAcceptedFJjet <<
  "hasRealFJjet="   << fhasRealFJjet <<
  "hasAcceptedEMCjet="   << fhasAcceptedEMCjet <<
  "hasRealEMCjet="   << fhasRealEMCjet <<
  "\n";
}

//________________________________________________________________________
void AliAnalysisJetHadro::FindJetsFJGen()
{

  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsFJGen ===== " << std::endl;
  //
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  { // track loop
    //
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    if (isTPCPileup || isITSPileup) continue;
    //
    // initialize the dummy particle id
    fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    //
    // check the origin of the track
    Int_t trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(iTrack))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(iTrack))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(iTrack)) trackOrigin = 2;
    if (trackOrigin<-1) continue;
    //
    Float_t ptotMCgen = trackMCgen->P();
    Float_t pTMCgen   = trackMCgen->Pt();
    Float_t phiMCGen  = trackMCgen->Phi();
    Float_t etaMCgen  = trackMCgen->Eta();
    Float_t rapMCgen  = trackMCgen->Y();
    //
    // select particle of interest
    Int_t iPart = -10;
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Int_t sign = (pdg>0) ? 1 : -1;
    if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; fPiMCgen = iPart;} // select pi+
    if (TMath::Abs(pdg) == kPDGka) {iPart = 2; fKaMCgen = iPart;} // select ka+
    if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; fPrMCgen = iPart;} // select pr+
    if (iPart == -10) continue;
    //
    // Resonance control
    Bool_t parInterest = (fPiMCgen>-1||fKaMCgen>-1||fPrMCgen>-1) ? kTRUE : kFALSE;
    Bool_t acceptRes = CheckIfFromResonance(1,trackMCgen,iTrack,parInterest,ptotMCgen,etaMCgen,fCentrality,kTRUE);
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;

    (*fTreeSRedirector)<<"jetsFJGen"<<
    "gid="       << fEventGID <<
    "acceptRes=" << acceptRes <<                // sample id for subsample method
    "part="      << iPart <<                // sample id for subsample method
    "origin="    << trackOrigin <<
    "sign="      << sign <<         // sign
    "p="         << ptotMCgen <<             // vertex momentum
    "pT="        << pTMCgen <<           // transverse momentum
    "eta="       << etaMCgen <<          // mc eta
    "rap="       << rapMCgen <<          // mc eta
    "phi="       << phiMCGen <<          // mc eta
    "cent="      << fCentrality <<     // Centrality
    "vZ="        << fVz <<
    "\n";


  } // ======= end of track loop for generated particles to see distributions =======

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
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
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
    Bool_t fBit96_spd    = fESDtrackCuts_Bit96_spd->AcceptTrack(track);
    Bool_t fBit96_sdd    = fESDtrackCuts_Bit96_sdd->AcceptTrack(track);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(track);
    if (!track->GetInnerParam()) continue;               // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    // Get the track variables
    Float_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track,0);
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
    if (fFillTracks && !fFillOnlyHists)
    {
      if(!fTreeSRedirector) return;
      if (fFillDistributions) FillTrackVariables(track);
      (*fTreeSRedirector)<<"tracks"<<
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "bit96spd="  << fBit96_spd <<     // cut account for the regions where no ITS-TPC match
      "bit96sdd="  << fBit96_sdd <<     // cut account for the regions where no ITS-TPC match
      "bit96="     << fBit96_base <<    // tight cuts of 2011 tuned data
      "bit128="    << fBit128 <<        // TPC only tracks cuts
      "bit768="    << fBit768 <<        // Hybrid track cuts
      "pixCut="    << ifDCAcutIfNoITSPixel <<    // cut: apply a DCAcut If No ITS Pixel
      "run="       << fRunNo <<                  // run Number
      "bField="    << fBField <<                 // magnetic filed
      "pileupbit=" << fPileUpBit <<              // flag for pileup selection
      "primMult="  << fNContributors <<          //  #prim tracks
      "tpcClMult=" << tpcClusterMultiplicity <<  //  TPC cluster multiplicity
      //
      "dcabase="  << dcaBaseCut <<  //  TPC multiplicity
      "dca10h="   << dca10h <<  //  TPC multiplicity
      "dca11h="   << dca11h <<  //  TPC multiplicity
      //
      "tpcmult="   << fTPCMult <<                //  TPC track multiplicity
      "itsmult="   << itsNumberOfTracklets <<    // ITS multiplicity
      "itsclmult=" << nITSClusters <<    // ITS multiplicity
      "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
      //
      "gid="       << fEventGID             <<  //  global event ID
      "eventtime=" << fTimeStamp            <<  // event timeStamp
      "intrate="   << fIntRate              <<  // interaction rate
      "cutBit="    << fTrackCutBits         <<  //  Systematic Cuts
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "cent="      << fCentrality           <<  //  centrality
      //
      "phi="       << fPhi                  <<  //  phi
      "dcaxy="     << fTrackDCAxy           <<  // dca cut on xy plane
      "dcaz="      << fTrackDCAz            <<  // dca cut along z
      "ncltpc="    << fNcl                  <<  // number of clusters
      "cRows="     << fTrackTPCCrossedRows  <<  // crossed Rows in TPC
      "chi2tpc="   << fTrackChi2TPC         <<  // TPC chi2
      "missCl="    << fMissingCl            <<  // fraction of missing clusters
      //
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

      "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi=" << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa=" << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr=" << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "dEdxMeanDe=" << fDEdxDe              <<
      "dEdxSigmaDe=" << fSigmaDe            <<

      "closestTPCPIDtype=" << closestPar[1]         << //particle type
      "closestTPCPIDmass=" << closestPar[2]         << //particle mass
      "\n";
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

    if(!fEffMatrix && fAcceptance && !fMCtrue){
      Double_t exMean[5]  = {fDEdxEl,  fDEdxPi,  fDEdxKa,  fDEdxPr,  fDEdxDe};
      Double_t exSigma[5] = {fSigmaEl, fSigmaPi, fSigmaKa, fSigmaPr, fSigmaDe};
      if ( ndEdxTPCall && nSigmaTPCall && ifDefaultCuts )
      {
        for (Int_t iPart = 0; iPart< 5; iPart++){
          Double_t weightExpected[7] = {Double_t(iPart),Double_t(fSign),fCentrality,fEta,fPtot, exSigma[iPart], exMean[iPart]};
          fHnExpected->Fill(weightExpected);
        }
      }
    }
    //
    // --------------------------------------------------------------
    //  Fill thnsparse for inclusive data
    // --------------------------------------------------------------
    //
    if (fUseThnSparse){
      Double_t trackdEdx[5] = {Double_t(fSign),fCentrality, fEta,fPtot, fTPCSignal};
      if(fUseThnSparse) fHndEdx->Fill(trackdEdx);
    }

  }// end of track loop

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillTrackVariables(AliESDtrack *track)
{

  /*

  // ITSTPC standard cuts
  // TPC
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  // ITS
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

  // Additional cuts
  esdTrackCuts->SetDCAToVertex2D(kTRUE);
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);

  */


  //
  // Apply downclaing condition
  Double_t qP      = track->Charge()/track->P();
  if ( !(fRandom.Rndm()*(qP*qP) < 0.05) ) return;
  AliTPCdEdxInfo tpcdEdxInfo;
  //
  // Calculate further variables
  Int_t tpcCrossedRows=0, tpcSignalN=0;
  Double_t eta   = track->Eta();
  Double_t tgl   = track->Pz()/track->Pt();
  Double_t phi   = track->Phi()-TMath::Pi(); // ????
  Int_t sign     = track->GetSign();
  Double_t phi85 = track->GetParameterAtRadius(85,5,7);
  ULong64_t flag = track->GetStatus();
  Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
  Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
  Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
  Int_t nclTPC   = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
  Int_t nclITS   = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
  Int_t nclTRD   = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
  Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
  Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
  Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
  Double_t itsdEdx = track->GetITSsignal();
  Double_t trddEdx = track->GetTRDsignal();
  Double_t tpcdEdx = track->GetTPCsignal();
  Double_t ptot0   = track->GetP();
  Double_t pt      = track->Pt();
  UChar_t itsclmap = track->GetITSClusterMap();
  Float_t pv[2],cov[3];
  track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  //
  // --------------------------------------------------------------
  //      TPC related observables
  // --------------------------------------------------------------
  //
  Double_t sharedTPCClusters=0.;
  Double_t phiTPC=-100.;
  Double_t ptotTPC=0.;
  Double_t lengthInActiveZone=0.;
  Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
  if (track->GetInnerParam()){
    track->GetTPCdEdxInfo(tpcdEdxInfo);
    phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
    ptotTPC = track->GetInnerParam()->GetP();
    lengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    track->GetImpactParametersTPC(pTPC,covTPC);
    tpcCrossedRows = track->GetTPCCrossedRows();
    tpcSignalN = track->GetTPCsignalN();
    Float_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track,0);
    if(nclTPC) sharedTPCClusters = static_cast<Double_t>(track->GetTPCnclsS())/static_cast<Double_t>(nclTPC);

  }
  //
  // pt dependent nclusters cut
  TFormula *f1NClustersTPCLinearPtDepSmall = new TFormula("f1NClustersTPCLinearPtDepSmall","60.+30./20.*x");
  TFormula *f1NClustersTPCLinearPtDep      = new TFormula("f1NClustersTPCLinearPtDep"     ,"70.+30./20.*x");
  TFormula *f1NClustersTPCLinearPtDepLArge = new TFormula("f1NClustersTPCLinearPtDepLArge","80.+30./20.*x");
  Int_t nClsPtDepSmall = (Int_t)(f1NClustersTPCLinearPtDepSmall->Eval(track->Pt()));
  Int_t nClsPtDep      = (Int_t)(f1NClustersTPCLinearPtDep     ->Eval(track->Pt()));
  Int_t nClsPtDepLarge = (Int_t)(f1NClustersTPCLinearPtDepLArge->Eval(track->Pt()));
  delete f1NClustersTPCLinearPtDepSmall;
  delete f1NClustersTPCLinearPtDep;
  delete f1NClustersTPCLinearPtDepLArge;
  //
  // --------------------------------------------------------------
  //   Fill downscaled tree
  // --------------------------------------------------------------
  //
  if(!fTreeSRedirector) return;
  if (fFilldscaledTree)
  {
    (*fTreeSRedirector)<<"dscaled"<<
    "gid="                  << fEventGID             <<  //  global event ID
    "eventtime="            << fTimeStamp            <<
    "intrate="              << fIntRate              <<  // interaction rate
    "cutBit="               << fTrackCutBits         <<  //  Systematic Cuts
    "dEdx="                 << fTPCSignal            <<  //  dEdx of the track
    "sign="                 << fSign                 <<  //  charge
    "ptot="                 << fPtot                 <<  //  TPC momentum
    "p="                    << fPVertex              <<  //  TPC momentum
    "pT="                   << fPt                   <<
    "eta="                  << fEta                  <<  //  eta
    "cent="                 << fCentrality           <<  //  centrality
    "nclTPC="               << nclTPC                <<  //  #ITS clusters
    "dcaxy="                << fTrackDCAxy           <<  // fD pv[0]
    "dcaz="                 << fTrackDCAz            <<  // fZ pv[1]
    //
    "isOnITS="              << isOnITS         <<
    "isOnTRD="              << isOnTRD         <<
    //
    "itsclmap="             << itsclmap              <<  //  vertex Z
    "flag="                 << flag                  <<
    "sharedTPCClusters="    << sharedTPCClusters     <<
    "tpcSignalN="           << tpcSignalN            <<  //  number of cl used in dEdx
    "cRows="                << tpcCrossedRows        <<  //  crossed rows
    "lengthInActiveZone="   << lengthInActiveZone    <<  //  fTrackLengthInActiveZone in TPC
    "phi="                  << phi                   <<  //  ph
    "phiTPC="               << phiTPC                <<
    "phi85="                << phi85                 <<
    "qP="                   << qP                    <<  //  charge/momentu,
    "tgl="                  << tgl                   <<  //  tangent
    "pt="                   << pt                    <<  //  pT
    //
    "vz="                   << fVz                   <<  //  vertex Z
    "fdTPC="                << pTPC[0]               <<
    "fzTPC="                << pTPC[1]               <<
    "fCddTPC="              << covTPC[0]                <<  //  DCAxy
    "fCdzTPC="              << covTPC[1]                <<  //  DCAz
    "fCzzTPC="              << covTPC[2]                <<  //  DCAz
    //
    "fd="                   << pv[0]               <<
    "fz="                   << pv[1]               <<
    "fCdd="                 << cov[0]                <<  //  DCAxy
    "fCdz="                 << cov[1]                <<  //  DCAz
    "fCzz="                 << cov[2]                <<  //  DCAz
    //
    "itsdEdx="              << itsdEdx               <<
    "trddEdx="              << trddEdx               <<
    "nclits="               << nclITS                <<  //  #ITS clusters
    "ncltrd="               << nclTRD                <<  //  #TRD clusters
    "chi2tpc="              << chi2TPC               <<  //  TPC chi2
    "chi2its="              << chi2ITS               <<  //  ITS chi2
    "chi2trd="              << chi2TRD               <<  //  TRD chi2
    //
    "itspixel01="           << fIsITSpixel01         <<
    "primRes="              << fPrimRestriction      <<
    "\n";
  }

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillMCFull()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  Int_t runNumber = fESD->GetRunNumber();
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillMCFull ===== " << std::endl;
  //
  // ======================================================================
  // For Marian
  AliTPCdEdxInfo tpcdEdxInfo;
  // ======================================================================
  //
  // Fill dEdx tree for MC closure
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // initialize the dummy particle id
    fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.; fDeMC =-100.; fMuMC =-100.;
    // Esd track
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
    TParticle *trackMC  = fMCStack->Particle(lab);
    Int_t pdg = trackMC->GetPdgCode();
    //
    Float_t dcaprim[2], covprim[3];
    Float_t dcaweak[2], covweak[3];
    Float_t dcamaterial[2], covmaterial[3];
    if(fMCStack->IsSecondaryFromMaterial(lab))  trackReal->GetImpactParameters(dcamaterial, covmaterial);
    if(fMCStack->IsSecondaryFromWeakDecay(lab)) trackReal->GetImpactParameters(dcaweak,     covweak);
    if(fMCStack->IsPhysicalPrimary(lab))        trackReal->GetImpactParameters(dcaprim,     covprim);
    //
    // track cuts
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t fBit96_spd = fESDtrackCuts_Bit96_spd->AcceptTrack(trackReal);
    Bool_t fBit96_sdd = fESDtrackCuts_Bit96_sdd->AcceptTrack(trackReal);
    Bool_t fBit96_base = fESDtrackCuts_Bit96->AcceptTrack(trackReal);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!fWeakAndMaterial) {
      if (!fMCStack -> IsPhysicalPrimary(lab)) continue;
      if (!trackReal -> GetInnerParam()) continue;
      if (!fESDtrackCutsLoose -> AcceptTrack(trackReal)) continue; // real track cuts
      if (!(trackReal->GetTPCsignalN()>0)) continue;
    }
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGel) iPart = 0; // select el
    if (TMath::Abs(pdg) == kPDGpi) iPart = 1; // select pi
    if (TMath::Abs(pdg) == kPDGka) iPart = 2; // select ka
    if (TMath::Abs(pdg) == kPDGpr) iPart = 3; // select pr
    if (TMath::Abs(pdg) == kPDGde) iPart = 4; // select de
    if (TMath::Abs(pdg) == kPDGmu) iPart = 5; // select mu
    //
    if (iPart == -10) continue;
    //
    if (iPart == 0 ) fElMC = trackReal->GetTPCsignal();
    if (iPart == 1 ) fPiMC = trackReal->GetTPCsignal();
    if (iPart == 2 ) fKaMC = trackReal->GetTPCsignal();
    if (iPart == 3 ) fPrMC = trackReal->GetTPCsignal();
    if (iPart == 4 ) fDeMC = trackReal->GetTPCsignal();
    if (iPart == 5 ) fMuMC = trackReal->GetTPCsignal();
    //
    fEtaMC        = trackReal->Eta();
    fPtotMC       = trackReal->GetInnerParam()->GetP();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    fTPCSignalMC  = trackReal->GetTPCsignal();
    fPxMC         = trackReal->Px();
    fPyMC         = trackReal->Py();
    fPzMC         = trackReal->Pz();
    fMissingCl    = trackReal->GetTPCClusterInfo(3,0,0,159);
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    trackReal->GetTPCdEdxInfo(tpcdEdxInfo);
    Float_t pv[2],cov[3];
    trackReal->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    Int_t itsNcls = trackReal->GetITSNcls();
    //
    if (trackReal -> GetInnerParam())  {
      Float_t closestPar[3];
      GetExpecteds(trackReal,closestPar);
      SetCutBitsAndSomeTrackVariables(trackReal,iPart);
    }
    //
    // --------------------------------------------------------------
    //                        Fill the trees
    // --------------------------------------------------------------
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    if (fWeakAndMaterial){
      (*fTreeSRedirector)<<"fTreeMC"<<
      "defCut="    << ifDefaultCuts <<  // default cuts
      "bit96spd="  << fBit96_spd <<  // run Number
      "bit96sdd="  << fBit96_sdd <<  // run Number
      "bit96="     << fBit96_base <<  // run Number
      "bit128="    << fBit128 <<  // run Number
      "bit768="    << fBit768 <<  // run Number
      "pixCut="    << ifDCAcutIfNoITSPixel <<  // run Number
      "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
      "gid="       << fEventGID             <<  //  global event ID
      "part="      << iPart <<
      "ptot="      << fPtotMC <<         // mc momentum
      "p="         << pMC <<
      "pT="        << fPtMC <<         // mc momentum
      "eta="       << fEtaMC <<          // mc eta
      "phi="       << fPhiMC <<          // mc eta
      "cent="      << fCentrality <<     // Centrality
      "centimp="   << fCentImpBin <<
      "sign="      << fSignMC <<         // sign
      "fZ="        << fTrackDCAz <<
      "fD="        << fTrackDCAxy <<
      "fCdd="      << cov[0] <<
      "fCdz="      << cov[1] <<
      "fCzz="      << cov[2] <<
      "missCl="    << fMissingCl <<
      "dEdxInfo.=" << &tpcdEdxInfo <<
      "nMult="     << fNContributors <<
      "tpcMult="   << fTPCMult <<
      "itsNcls="   << itsNcls <<
      "dcaxyprim="      << dcaprim[0] <<           // electron dEdx
      "dcaxyweak="      << dcaweak[0] <<           // electron dEdx
      "dcaxymaterial="  << dcamaterial[0] <<           // electron dEdx
      "dcazprim="       << dcaprim[1] <<           // electron dEdx
      "dcazweak="       << dcaweak[1] <<           // electron dEdx
      "dcazmaterial="   << dcamaterial[1] <<           // electron dEdx
      "\n";
    } else {
      (*fTreeSRedirector)<<"fTreeMC"<<
      "isample="   << sampleNo <<        // sample id for subsample method
      "gid="       << fEventGID <<       //  global event ID
      "part="      << iPart <<
      "dEdx="      << fTPCSignalMC <<    // dEdx of mc track
      "cutBit="    << fTrackCutBits <<   //  Systematic Cuts
      "sign="      << fSignMC <<         // sign
      "ptot="      << fPtotMC <<         // TPC momentum
      "p="         << pMC <<             // momentum at vertex
      "pT="        << fPtMC <<           // mc momentum
      "eta="       << fEtaMC <<          // mc eta
      "phi="       << fPhiMC <<          // mc eta
      "cent="      << fCentrality <<     // Centrality
      "centimp="   << fCentImpBin <<
      //
      "dcaxy="     << fTrackDCAxy           <<
      "dcaz="      << fTrackDCAz            <<
      "ncltpc="    << fNcl                  <<  //  centrality
      "ncltpccorr="<< fNclCorr                  <<  //  centrality
      "cRows="     << fTrackTPCCrossedRows  <<
      "chi2tpc="   << fTrackChi2TPC         <<
      "chi2tpccorr=" << fTrackChi2TPCcorr         <<
      // "intrate="   << fIntRate              <<  // interaction rate
      "\n";
    }

  } // ======= end of track loop for MC dEdx filling =======

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillEffMatrix()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  Int_t trackOrigin = -10;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEffMatrix ===== " << std::endl;
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------   reconstructed MC particles  ------------------------------
  // -----------------------------------------------------------------------------------------
  //
  //
  Bool_t bEventVertexZSmall      = (TMath::Abs(fVz)<6 && TMath::Abs(fVz)>0.2);
  Bool_t bCutReference           = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge      = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  Bool_t bEventVertexZALICE      = (TMath::Abs(fVz)<10);
  Bool_t bEventVertexZALICETight = (TMath::Abs(fVz)<7);
  //
  // loop over tracks
  //
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  { // track loop
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    Int_t lab = TMath::Abs(trackReal->GetLabel()); // avoid from negatif labels, they include some garbage

    //
    Float_t etaMCrec = (fRapidityType==0) ? trackReal->Eta() :  trackReal->Y();
    if ((etaMCrec<fEtaDown) || (etaMCrec>fEtaUp)) continue;   // eta [-0.8,0.8]
    if (!fMCStack->IsPhysicalPrimary(lab)) continue;   // MC primary track check
    //
    Bool_t ifDefaultCuts = fESDtrackCuts->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!trackReal->GetInnerParam()) continue;   // If track in TPC
    if (!ifDefaultCuts) continue;   // apply esdtrack cuts
    if (!ifDCAcutIfNoITSPixel) continue;   // apply esdtrack cuts
    //
    // access generated track info
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
    Int_t pdg = trackMCgen->Particle()->GetPdgCode();  //  Int_t pdg = trackMC->GetPdgCode();   TODO
    //
    // get track info
    Float_t fPtRec = 0.;
    if(fUsePtCut==1) fPtRec = trackReal->P();
    if(fUsePtCut==2) fPtRec = trackReal->Pt();
    if ((fPtRec<fMomDown) || (fPtRec>fMomUp)) continue;
    //
    Int_t fPartID  = -10;
    if (TMath::Abs(pdg) == kPDGpi) fPartID=0; // select pi
    if (TMath::Abs(pdg) == kPDGka) fPartID=1; // select ka
    if (TMath::Abs(pdg) == kPDGpr) fPartID=2; // select pr
    if (fPartID == -10) continue;
    //
    // Loop over all track settings
    Float_t closestPar[3];
    GetExpecteds(trackReal,closestPar);
    SetCutBitsAndSomeTrackVariables(trackReal,fPartID);
    //
    // Fill TPC eff matrix
    Double_t xxxRec[6]={0.,0.,Float_t(fPartID),fCentrality,fPtRec,etaMCrec};
    if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRec);
    if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRec);
    //
    // Fill TOF eff matrix
    Bool_t piTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kPion))  <=2.5);
    Bool_t kaTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kKaon))  <=2.5);
    Bool_t prTOF = (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(trackReal, AliPID::kProton))<=2.5);
    if ( (piTOF && fPartID==0) ||  (kaTOF && fPartID==1) || (prTOF && fPartID==2) ) {
      Double_t xxxRecTOF[6]={1.,0.,Float_t(fPartID),fCentrality,fPtRec,etaMCrec};
      if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecTOF);
      if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecTOF);
    }


  } // ======= end of rec track loop =======
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------   MC generated pure MC particles  --------------------------
  // -----------------------------------------------------------------------------------------
  //
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  { // track loop
    //
    // Select real trigger event and reject other pile up vertices
    Bool_t isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iTrack,fMCEvent);
    Bool_t isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    if (isTPCPileup || isITSPileup) continue;
    //
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    if (!trackMCgen) continue;
    //
    // apply primary track and acceptance cuts
    if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
    //
    // Aplly eta acceptance
    Float_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
    if ((etaMCgen<fEtaDown) || (etaMCgen>fEtaUp))  continue;
    //
    // Aplly momentum acceptance
    Float_t ptotMCgen =0.;
    if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
    if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
    if ((ptotMCgen<fMomDown) || (ptotMCgen>fMomUp)) continue;
    //
    // Efficiency matices for individual particles
    Int_t fPartID  = -10;
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    if (TMath::Abs(pdg) == kPDGpi) fPartID=0; // select pi
    if (TMath::Abs(pdg) == kPDGka) fPartID=1; // select ka
    if (TMath::Abs(pdg) == kPDGpr) fPartID=2; // select pr
    if (fPartID == -10) continue;
    //
    Double_t xxxGen[6]={0.,0.,Float_t(fPartID),fCentrality,ptotMCgen,etaMCgen};
    if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGen);
    if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGen);
    // generated does not know about TOF
    Double_t xxxGenTOF[6]={1.,0.,Float_t(fPartID),fCentrality,ptotMCgen,etaMCgen};
    if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenTOF);
    if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenTOF);

  } // ======= end of gen track loop =======

}
//________________________________________________________________________
void AliAnalysisJetHadro::FillCleanSamples()
{

  // Fill Clean Pions from K0s
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillCleanSamples ===== " << std::endl;
  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  AliKFParticle::SetField(fESD->GetMagneticField());
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  Double_t mm[3] = {0,0,0};
  const Double_t cProtonMass  =TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass    =TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass=TDatabasePDG::Instance()->GetParticle(11)->Mass();
  //
  // Selection From Ionut
  const AliESDVertex *primaryVertex = fESD->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
  if(fV0OpenCuts) {
    fV0OpenCuts->SetEvent(fESD);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  if(fV0StrongCuts) {
    fV0StrongCuts->SetEvent(fESD);
    fV0StrongCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  //
  //
  TObjArray* listCrossV0 = fESDtrackCutsV0->GetAcceptedV0s(fESD);
  Int_t nGoodV0s         = listCrossV0->GetEntries();
  delete listCrossV0;
  //
  // Loop over V0s
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s; iV0MI++) {
    //
    Int_t v0purity = 0;
    AliESDv0 * fV0s = fESD->GetV0(iV0MI);
    Int_t lOnFlyStatus = 0;
    lOnFlyStatus = fV0s->GetOnFlyStatus();
    if (!lOnFlyStatus) {fTrackCutBits=0; continue;}
    //
    AliESDtrack* trackPosTest = fESD->GetTrack(fV0s->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0s->GetNindex());
    //
    // ----------------------------------------------------------------------------------------------------------
    //  Selections from ionuts
    // ----------------------------------------------------------------------------------------------------------
    //
    if(trackPosTest->GetSign() == trackNegTest->GetSign()) {fTrackCutBits=0; continue;}
    Bool_t v0ChargesAreCorrect = (trackPosTest->GetSign()==+1 ? kTRUE : kFALSE);
    trackPosTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetNindex()) : trackPosTest);
    trackNegTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetPindex()) : trackNegTest);
    //
    Bool_t goodK0s = kTRUE, goodLambda = kTRUE, goodALambda = kTRUE, goodGamma = kTRUE;
    if(fV0OpenCuts) {
      goodK0s = kFALSE, goodLambda = kFALSE, goodALambda = kFALSE, goodGamma = kFALSE;
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0OpenCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) {
        goodK0s = kTRUE;
        if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(trackPosTest) || !fK0sPionCuts->IsSelected(trackNegTest))) goodK0s = kFALSE;
      }
      if (processV0 && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodLambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackPosTest)) goodLambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackNegTest)) goodLambda = kFALSE;
      }
      if (processV0 && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodALambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackNegTest)) goodALambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackPosTest)) goodALambda = kFALSE;
      }
      if (processV0 && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) {
        goodGamma = kTRUE;
        if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(trackPosTest) || !fGammaElectronCuts->IsSelected(trackNegTest))) goodGamma = kFALSE;
      }
    }
    //
    Bool_t veryGoodK0s = kFALSE, veryGoodLambda = kFALSE, veryGoodALambda = kFALSE, veryGoodGamma = kFALSE;
    if(fV0StrongCuts && (goodK0s || goodLambda || goodALambda || goodGamma)) {
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0StrongCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && goodK0s     && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) veryGoodK0s = kTRUE;
      if (processV0 && goodLambda  && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodLambda = kTRUE;
      if (processV0 && goodALambda && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodALambda = kTRUE;
      if (processV0 && goodGamma   && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) veryGoodGamma = kTRUE;
    }
    //
    //
    if( goodK0s || goodLambda || goodALambda || goodGamma ) v0purity = 1;
    if( veryGoodK0s || veryGoodLambda || veryGoodALambda || veryGoodGamma ) v0purity = 2;
    //
    // ----------------------------------------------------------------------------------------------------------
    //  My cuts
    // ----------------------------------------------------------------------------------------------------------
    //
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackPosTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackNegTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!trackPosTest->GetInnerParam()) {fTrackCutBits=0; continue;}
    if (!trackNegTest->GetInnerParam()) {fTrackCutBits=0; continue;}

    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0s->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0s->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    fV0s->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);

    if ((vecP.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}
    if ((vecN.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}
    Double_t thetaP  = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN  = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    if ( ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) <0.00001) {fTrackCutBits=0; continue;}
    fAlfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    fQt   = vecP.Mag()*sin(thetaP);
    if (fUseCouts) fHistArmPod->Fill(fAlfa,fQt);
    // fV0s->ChangeMassHypothesis(22);   // ????
    // fV0s->ChangeMassHypothesis(310); // ????
    //
    // main armentoros podolanki cuts
    if (TMath::Abs(fAlfa)>0.9) {fTrackCutBits=0; continue;}
    if (fQt >0.22) {fTrackCutBits=0; continue;}
    if (fQt >0.02 && fQt<0.12 && TMath::Abs(fAlfa)<0.4) {fTrackCutBits=0; continue;}
    SelectCleanSamplesFromV0s(fV0s,trackPosTest,trackNegTest);
    //
    TLorentzVector posE, negE, photon, posP, negP, posPi, negPi, lambda, antiLambda, kaon, posProton, k0sProton;
    negE.SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
    posE.SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    posPi.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    negP.SetXYZM(mn[0],mn[1],mn[2],cProtonMass);
    posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
    kaon=posPi+negPi;
    photon=posE+negE;
    lambda=posP+negPi;
    antiLambda=posPi+negP;
    //
    //
    Float_t ptotForBetaGamma0 = trackPosTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGamma1 = trackNegTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t posNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kPion));
    Double_t negNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kPion));
    Double_t posNTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kProton)) : 0.;
    Double_t negNTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kProton)) : 0.;
    Double_t posNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kElectron));
    Double_t negNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kElectron));
    //
    // --------------------------------------------------------------
    //  Invariant mass cuts
    // --------------------------------------------------------------
    //
    Bool_t isK0sMass        = (kaon.M()>0.485 && kaon.M()<0.51); // (kaon.M()>0.490 && kaon.M()<0.504);
    Bool_t isLambdaMass     = (lambda.M()>1.112 && lambda.M()<1.119); //(lambda.M()>1.113 && lambda.M()<1.118);
    Bool_t isAntiLambdaMass = (antiLambda.M()>1.112 && antiLambda.M()<1.119); // (antiLambda.M()>1.113 && antiLambda.M()<1.118);
    Bool_t isPhotonMass     = (photon.M()<0.005); // (photon.M()<0.005);
    Double_t oneLegSigma = 3.5;
    if (fQt<0.02 && TMath::Abs(fAlfa)<0.5){
      // Cuts concerns Electrons
      //
      // Apply one leg cut for electrons
      if (!(negNTPCSigmaEl<2. || posNTPCSigmaEl<2.)) {fTrackCutBits=0; continue;}
      //
      if (fUseCouts) fHistInvPhoton->Fill(photon.M());
      if (isK0sMass) {fTrackCutBits=0; continue;}
      if (isLambdaMass) {fTrackCutBits=0; continue;}
      if (isAntiLambdaMass) {fTrackCutBits=0; continue;}
      if (!isPhotonMass) {fTrackCutBits=0; continue;}
    } else {
      if (fUseCouts) {
        fHistInvK0s->Fill(kaon.M());
        fHistInvLambda->Fill(lambda.M());
        fHistInvAntiLambda->Fill(antiLambda.M());
      }
      if ( !(isK0sMass || isLambdaMass || isAntiLambdaMass) ) {fTrackCutBits=0; continue;}
      //
      // Apply one leg cut for K0s
      if (fQt>0.11 && (!(negNTPCSigmaPi < oneLegSigma || posNTPCSigmaPi < oneLegSigma))) {fTrackCutBits=0; continue;}
      // //
      // // Apply one leg cut for antilambda
      // if (fQt<1.1 && fAlfa<0) {
      //   if (!(negNTPCSigmaPr < oneLegSigma || posNTPCSigmaPi < oneLegSigma)) {fTrackCutBits=0; continue;}
      //   if (negNTPCSigmaPr > 4) {fTrackCutBits=0; continue;}
      // }
      // //
      // // Apply one leg cut for lambda
      // if (fQt<1.1 && fAlfa>0) {
      //   if (!(negNTPCSigmaPi < oneLegSigma || posNTPCSigmaPr < oneLegSigma)) {fTrackCutBits=0; continue;}
      //   if (posNTPCSigmaPr > 4) {fTrackCutBits=0; continue;}
      // }
    }
    //
    // Set the variables to be filled in the tree
    Bool_t selectPosLegs = (posNTPCSigmaPi<3 || posNTPCSigmaPr<3 || posNTPCSigmaEl<3);
    Bool_t selectNegLegs = (negNTPCSigmaPi<3 || negNTPCSigmaPr<3 || negNTPCSigmaEl<3);
    //
    // positive leg
    UInt_t cutBit0 = 0, cutBit1=0, fTrackCutBits=0;
    cutBit0 = SetCutBitsAndSomeTrackVariables(trackPosTest,0);
    Float_t dEdx0    = trackPosTest->GetTPCsignal();
    Float_t itsdEdx0 = trackPosTest->GetITSsignal();
    Float_t ptot0  = trackPosTest->GetInnerParam()->GetP();
    Float_t p0     = trackPosTest->P();
    Float_t pT0    = trackPosTest->Pt();
    Float_t eta0   = trackPosTest->Eta();
    Int_t sign0    = trackPosTest->Charge();
    Float_t phi0   = trackPosTest->Phi();
    Float_t nSigmasPiTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // negative leg
    cutBit1 = SetCutBitsAndSomeTrackVariables(trackNegTest,0);
    Float_t dEdx1  = trackNegTest->GetTPCsignal();
    Float_t itsdEdx1 = trackNegTest->GetITSsignal();
    Float_t ptot1  = trackNegTest->GetInnerParam()->GetP();
    Float_t p1     = trackNegTest->P();
    Float_t pT1    = trackNegTest->Pt();
    Float_t eta1   = trackNegTest->Eta();
    Int_t sign1    = trackNegTest->Charge();
    Float_t phi1   = trackNegTest->Phi();
    Float_t nSigmasPiTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // --------------------------------------------------------------
    //  Fill Clean Samples tree
    // --------------------------------------------------------------
    //
    if (selectNegLegs || selectPosLegs)
    {
      if (fFillArmPodTree)
      {
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"fArmPodTree"<<
        "gid="                  << fEventGID             <<  //  global event ID
        "eventtime="            << fTimeStamp            <<
        "intrate="              << fIntRate              <<  // interaction rate
        "piFromK0="             << fCleanPionsFromK0    <<  // K0s cut for pions
        "v0haspixel="           << fHasV0FirstITSlayer  <<  // ITS pixel cout
        "purity="               << v0purity             <<
        "qt="                   << fQt                  <<  // qT
        "alfa="                 << fAlfa                <<  // alpha
        "cent="                 << fCentrality          <<  // centrality
        //
        "cutBit0="               << cutBit0        <<  // cut bits
        "itsdEdx0="              << itsdEdx0       <<  // TPC dEdx
        "dEdx0="                 << dEdx0          <<  // TPC dEdx
        "sign0="                 << sign0          <<
        "ptot0="                 << ptot0          <<  // momentum
        "p0="                    << p0             <<
        "pT0="                   << pT0            <<
        "eta0="                  << eta0           <<  // eta
        "phi0="                  << phi0           <<  // eta
        "nSigmasPiTOF0="         << nSigmasPiTOF0  <<  // TOF nsigma cut for pions
        "nSigmasPrTOF0="         << nSigmasPrTOF0  <<  // TOF nsigma cut for protons
        //
        "cutBit1="               << cutBit1        <<  // cut bits
        "dEdx1="                 << dEdx1          <<  // TPC dEdx
        "itsdEdx1="              << itsdEdx1       <<  // TPC dEdx
        "sign1="                 << sign1          <<
        "ptot1="                 << ptot1          <<  // momentum
        "p1="                    << p1             <<
        "pT1="                   << pT1            <<
        "eta1="                  << eta1           <<  // eta
        "phi1="                  << phi1           <<  // eta
        "nSigmasPiTOF1="         << nSigmasPiTOF1  <<  // TOF nsigma cut for pions
        "nSigmasPrTOF1="         << nSigmasPrTOF1  <<  // TOF nsigma cut for protons
        //
        "\n";

      }
    }
    cutBit0 = 0, cutBit1=0; fTrackCutBits=0;

  } // end of V0 loop

}
//________________________________________________________________________
void AliAnalysisJetHadro::GetExpecteds(AliESDtrack *track, Float_t closestPar[3])
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
  Int_t nSigmaTmp = 2;
  //
  // Electron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
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
  if(fEventCountInFile==2 && !fRunOnGrid && fillTree) {
    if(!fTreeSRedirector) return kFALSE;
    (*fTreeSRedirector)<<"resonance"<<
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
void AliAnalysisJetHadro::SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1)
{
  //
  // SetAliases and Metadata for the V0 trees
  //
  AliKFParticle kfparticle; //
  AliAnalysisTaskFilteredTree filteredV0;
  // Int_t type=filteredV0.GetKFParticle(v0,fESD,kfparticle);
  filteredV0.GetKFParticle(v0,fESD,kfparticle);
  //
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  Double_t massLambda = pdg->GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg->GetParticle("K0")->Mass();
  // const Double_t massProton  =pdg->GetParticle(kPDGpr)->Mass();
  // const Double_t massPion    =pdg->GetParticle(kPDGpi)->Mass();

  const Double_t livetimeK0=2.684341668932;  // livetime in cm (surpisely missing info in PDG - see root forum)
  const Double_t livetimeLambda=7.8875395;  // livetime in cm (missing info in PDG - see root forum)
  fHasTrack0FirstITSlayer = track0->HasPointOnITSLayer(0);
  fHasTrack1FirstITSlayer = track1->HasPointOnITSLayer(0);
  fHasV0FirstITSlayer = (fHasTrack0FirstITSlayer||fHasTrack1FirstITSlayer);
  //
  //
  Double_t v0Rr = v0->GetRr();  // rec position of the vertex CKBrev
  Double_t v0P  = v0->P();      // TMath::Sqrt(Px()*Px()+Py()*Py()+Pz()*Pz())
  Double_t livetimeLikeK0 = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massK0)*(v0P/massK0)+1)*livetimeK0));
  Double_t livetimeLikeLambda = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massLambda)*(v0P/massLambda)+1)*livetimeLambda));
  Double_t livetimeLikeGamma = v0Rr/80.;
  // Double_t livetimeLikeBkg   = v0Rr/80.;

  // delta of mass
  Double_t K0Delta = v0->GetEffMass(2,2)-massK0;        //   tree->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  Double_t LDelta  = v0->GetEffMass(4,2)-massLambda;    //   tree->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  Double_t ALDelta = v0->GetEffMass(2,4)-massLambda;    //   tree->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  Double_t EDelta  = v0->GetEffMass(0,0);               //   tree->SetAlias("EDelta","(v0.GetEffMass(0,0))");

  // pull of the mass
  if (v0->GetKFInfo(2,2,1)==0. || v0->GetKFInfo(4,2,1)==0. || v0->GetKFInfo(2,4,1)==0. || v0->GetKFInfo(0,0,1)==0.) return;
  Double_t K0Pull = (v0->GetEffMass(2,2)-massK0)/v0->GetKFInfo(2,2,1);        //   tree->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  Double_t LPull  = (v0->GetEffMass(4,2)-massLambda)/v0->GetKFInfo(4,2,1);    //   tree->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  Double_t ALPull = (v0->GetEffMass(2,4)-massLambda)/v0->GetKFInfo(2,4,1);    //   tree->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  Double_t EPull  = EDelta/v0->GetKFInfo(0,0,1);                              //   tree->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");
  Double_t K0Like0 = TMath::Exp(-K0Pull*K0Pull)*livetimeLikeK0;
  Double_t LLike0  = TMath::Exp(-LPull*LPull)*livetimeLikeLambda;
  Double_t ALLike0 = TMath::Exp(-ALPull*ALPull)*livetimeLikeLambda;
  Double_t ELike0  = TMath::Exp(-abs(EPull)*0.2)*livetimeLikeGamma;
  Double_t V0Like  = TMath::Exp(-TMath::ACos(v0->GetV0CosineOfPointingAngle())*v0Rr/0.36)*TMath::Exp(-TMath::Sqrt(kfparticle.GetChi2())/0.5);

  //   tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts
  Int_t ntracks = fESD->GetNumberOfTracks();
  Double_t BkgLike = 0.000005*ntracks;    // backround coeefecint  to be fitted - depends on other cuts
  Double_t LLike = LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  Double_t ALLike = ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  Double_t tr0NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kPion));
  Double_t tr1NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion));
  Float_t ptotForBetaGamma0 = track0->GetInnerParam()->GetP();
  Float_t ptotForBetaGamma1 = track1->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  Double_t tr0NTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kProton)) : 0.;
  Double_t tr1NTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kProton)) : 0.;
  //   treeV0->SetAlias("cleanPion0FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   treeV0->SetAlias("cleanPion1FromK0","K0Like0>0.05&&K0Like0>(LLike0+ALLike0+ELike0)*3&&abs(K0Delta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8&&v0.PtArmV0()>0.06");
  fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  //   treeV0->SetAlias("cleanPion0FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.006&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DProton-50)<8");
  fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPr<2;
  //   treeV0->SetAlias("cleanPion1FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.006&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DProton-50)<8");
  fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPr<2;
  //   treeV0->SetAlias("cleanProton0FromLambda","LLike>0.05&&LLike0>(K0Like0+ALLike0+ELike0)*3&&abs(LDelta)<0.001&&V0Like>0.1&&abs(track1.fTPCsignal/dEdx1DPion-50)<8");
  fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2;
  //   treeV0->SetAlias("cleanProton1FromLambda","ALLike>0.05&&ALLike0>(K0Like0+LLike0+ELike0)*3&&abs(ALDelta)<0.001&&V0Like>0.1&&abs(track0.fTPCsignal/dEdx0DPion-50)<8");
  fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2;
  fCleanPionsFromK0 =  (fCleanPion0FromK0 || fCleanPion1FromK0);

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
UInt_t AliAnalysisJetHadro::SetCutBitsAndSomeTrackVariables(AliESDtrack *track, Int_t particleType)
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
    //this block is important
    fPtot      = track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    fTrackTPCSignalN     = track->GetTPCsignalN();
    fTrackTPCCrossedRows = Float_t(track->GetTPCCrossedRows());
    fTPCShared = track->GetTPCnclsS();
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
  Bool_t cleanPrTOF = ((TMath::Abs(nSigmasPrTOF)<=2.5));
  Bool_t cleanDeTOF = ((TMath::Abs(nSigmasDeTOF)<=3.0));
  Bool_t cleanKaTOF = ((TMath::Abs(nSigmasKaTOF)<=2.5));
  Bool_t cleanKaTOFTRD = ((TMath::Abs(nSigmasKaTOF)<=1.2) && TOFSignalDz<1. && TOFSignalDx<1. && nclsTRD>100);
  //
  // Systematic settings
  fTrackCutBits=0;
  //
  if (fTrackTPCCrossedRows>=60)  (fTrackCutBits |= 1 << kNCrossedRowsTPC60);
  if (fTrackTPCCrossedRows>=80)  (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
  if (fTrackTPCCrossedRows>=100) (fTrackCutBits |= 1 << kNCrossedRowsTPC100);
  //
  // Special treatment of the 2018 pass3 and 2015 pass2 data
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    if (fTrackChi2TPC<2.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);
    if (fTrackChi2TPC<2.5) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);
    if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);
  } else {
    if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
    if (fTrackChi2TPC<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
    if (fTrackChi2TPC<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
  }
  Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
  Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);    // 10h tuned loose cut
  Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);  // 10h tuned loose cut
  Bool_t dcaLoose   = TMath::Abs(fTrackDCAxy)<0.4;  // 10h tuned loose cut
  if (dca10h)     (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepSmall);
  if (dcaBaseCut) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDep);
  if (dcaLoose)   (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepLarge);
  //
  if (TMath::Abs(fTrackDCAz)<2.0) (fTrackCutBits |= 1 << kVertexZSmall);
  if (TMath::Abs(fTrackDCAz)<3.0) (fTrackCutBits |= 1 << kVertexZ);
  if (TMath::Abs(fTrackDCAz)<4.0) (fTrackCutBits |= 1 << kVertexZLarge);
  //
  if (TMath::Abs(fVz)<6 && TMath::Abs(fVz)>0.2 ) (fTrackCutBits |= 1 << kEventVertexZSmall);
  if (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15) (fTrackCutBits |= 1 << kEventVertexZ);
  if (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1 ) (fTrackCutBits |= 1 << kEventVertexZLarge);
  if (TMath::Abs(fVz)<10 ) (fTrackCutBits |= 1 << kEventVertexZALICE);
  if (TMath::Abs(fVz)<7  ) (fTrackCutBits |= 1 << kEventVertexZALICETight);
  //
  if (fTrackRequireITSRefit)                           (fTrackCutBits |= 1 << kRequireITSRefit);
  if (fTrackIsFirstITSlayer || fTrackIsSecondITSlayer) (fTrackCutBits |= 1 << kPixelRequirementITS);
  if (fTrackNewITScut)                                 (fTrackCutBits |= 1 << kNewITSCut);
  //
  // dangerous cuts
  if (fTrackLengthInActiveZone>=90)  (fTrackCutBits |= 1 << kActiveZoneSmall);
  if (fTrackLengthInActiveZone>=100) (fTrackCutBits |= 1 << kActiveZone);
  if (fTrackLengthInActiveZone>=110) (fTrackCutBits |= 1 << kActiveZoneLarge);
  //
  if (fTrackTPCSignalN>=60) (fTrackCutBits |= 1 << kTPCSignalNSmall);
  if (fTrackTPCSignalN>=70) (fTrackCutBits |= 1 << kTPCSignalN);
  if (fTrackTPCSignalN>=80) (fTrackCutBits |= 1 << kTPCSignalNLarge);
  //
  // --------------------------------------------------------------------
  //                    Clean sample selections
  // --------------------------------------------------------------------
  //
  // 2.5 nsigma TOF protons and kaons for amplitude estimation
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
void AliAnalysisJetHadro::SetSpecialV0Cuts(AliESDv0KineCuts* cuts)
{

  cuts->SetMode(0, 0);    // cuts->SetMode(mode, type); mode 0: purely kinematical selection   type 0: pp 1:PbPb not yet ready
  //
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  //
  // gamma cuts
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.05);
  // K0s cuts
  cuts->SetK0CutChi2NDF(10.0);
  Float_t cosPointK0s[2] = {0.0, 0.02};
  cuts->SetK0CutCosPoint(cosPointK0s);
  Float_t cutDCAK0s[2] = {0.0, 0.2};
  cuts->SetK0CutDCA(cutDCAK0s);
  Float_t vtxRK0s[2] = {2.0, 30.0};
  cuts->SetK0CutVertexR(vtxRK0s);
  Float_t k0sInvMass[2] = {0.486, 0.508};
  cuts->SetK0CutInvMass(k0sInvMass);
  // Lambda and anti-Lambda cuts
  cuts->SetLambdaCutChi2NDF(10.0);
  Float_t cosPointLambda[2] = {0.0, 0.02};
  cuts->SetLambdaCutCosPoint(cosPointLambda);
  Float_t cutDCALambda[2] = {0.0, 0.2};
  cuts->SetLambdaCutDCA(cutDCALambda);
  Float_t vtxRLambda[2] = {2.0, 40.0};
  cuts->SetLambdaCutVertexR(vtxRLambda);
  Float_t lambdaInvMass[2] = {1.11, 1.12};
  cuts->SetLambdaCutInvMass(lambdaInvMass);

}
//________________________________________________________________________
void AliAnalysisJetHadro::BinLogAxis(TH1 *h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the BinLogAxis ===== " << std::endl;
  TAxis *axis       = h->GetXaxis();
  Int_t bins        = axis->GetNbins();

  Double_t from     = axis->GetXmin();
  Double_t to       = axis->GetXmax();
  std::vector<double>  newBins;
  newBins.resize(bins + 1);

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins.data());

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
    fHistEmptyEvent->Fill(counterBin);
    std::cout << " Info::siweyhmi: Empty event in " << fChunkName << std::endl;
  }
  if (fUseCouts) std::cout << " Info::siweyhmi: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  return emptyCount;

}
//
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

}
