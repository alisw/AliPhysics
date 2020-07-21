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
//full code, 18r.pbpb_legotrain_01
// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
//#include "AliTriggerAnalysis.h"
#include "AliAODMCHeader.h"

// my headers
#include "AliAnalysisTaskUpcFourPi.h"

ClassImp(AliAnalysisTaskUpcFourPi);

using std::cout;
using std::endl;

//trees for UPC EtaC analysis,fList1KStar
// christopher.anson@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcFourPi::AliAnalysisTaskUpcFourPi() //first constructor, contains all histograms
	: AliAnalysisTaskSE(), fType(0), fTracking(0), isMC(kFALSE), fRunTree(kFALSE), fRunHist(kTRUE), fRunSystematics(kFALSE), fPIDResponse(0), fEtaCK0sChannelTree(0), fEtaCTree(0), fMeritCutChoice(0),
	fRunNum(0), fPerNum(0), fOrbNum(0), fL0inputs(0), fL1inputs(0),
	fTOFmask(0), fIsPhysicsSelected(kFALSE),
	fVtxContrib(0), fVtxChi2(0), fVtxNDF(0), fSpdVtxContrib(0),
	fBCrossNum(0), fNtracklets(0), fNLooseTracks(0),
	fZNAenergy(0), fZNCenergy(0), fZPAenergy(0), fZPCenergy(0), fZDCAtime(0), fZDCCtime(0), fV0Adecision(0), fV0Cdecision(0), fADAdecision(0), fADCdecision(0),
	fDataFilnam(0), fRecoPass(0), fEvtNum(0),
	fJPsiAODTracks(0), fJPsiESDTracks(0), fEtaCAODTracks(0), fEtaCESDTracks(0), fGenPart(0),
	fEtaCCandidatesPerChannel(0), fEtaCLowPtCandidatesPerChannel(0), fAllPtVsMinvEtaC(0), fAllMinvEtaCLowPt(0), fChannelVsMinvEtaC(0),
	//trigger hists and lists
	fListTrig(0), fHistCcup29TriggersPerRun(0), fHistCcup30TriggersPerRun(0), fHistCcup31TriggersPerRun(0),
	fHistZedTriggersPerRun(0), fHistMBTriggersPerRun(0), fHistCentralTriggersPerRun(0), fHistSemiCentralTriggersPerRun(0), fHistCtrueTriggersPerRun(0),
	//PId/ZDC hists and lists
	fListHist(0), fListHistPID(0), fListHistUnusedPID(0), fListHistPIDcuts(0), fHistNTracks(0), fHistPIDCheck(0),
	fListHistZDC(0), fHistZDCCuts(0),
	fPionTPCvsPionTOFLowP(0), fPionTPCvsPionTOFMidP(0), fPionTPCvsPionTOFHighP(0),
	fKaonTPCvsKaonTOFLowP(0), fKaonTPCvsKaonTOFMidP(0), fKaonTPCvsKaonTOFHighP(0),
	fPionTPCvsKaonTPCLowP(0), fPionTPCvsKaonTPCMidP(0), fPionTPCvsKaonTPCHighP(0), fPionTOFvsKaonTOFMidP(0), fPionTOFvsKaonTOFPionMidP(0), fPionTOFvsKaonTOFKaonMidP(0),
	fDedxVsPAll(0), fDedxVsPPion(0), fDedxVsPKaon(0), fDedxVsPBoth(0),
	fTPCdEdxVsTOFbetaAll(0), fTPCdEdxVsTOFbetaPionsWithPID(0), fTPCdEdxVsTOFbetaKaonsWithPID(0),
	fFourTrackMissing(0), fSixTrackMissing(0),
	fTOFbetaVsPAll(0), fTOFbetaVsPPion(0), fTOFbetaVsPKaon(0), fTOFbetaVsPBoth(0),
	fHistZDCAenergy(0), fHistZDCCenergy(0), fHistZDCAtime(0), fHistZDCCtime(0), fHistZDCImpactParameter(0), fHistZDCAImpactParameter(0), fHistZDCCImpactParameter(0),
	fHistFourTracksNboth(0), fHistFourTracksNneither(0), fHistSixTracksNboth(0), fHistSixTracksNneither(0),
	fListSystematics(0), fListJPsiLoose(0), fListJPsiTight(0), fListEtaCLoose(0), fListEtaCTight(0),
	fListNKaonChange(0), fListNPionChange(0), f4KnKaonVsnewKaon(0),
	fList4TrackPID(0), fList6TrackPID(0), fTOFIntegratedLength(0),
	//Kstar hists and lists
	fListHistKstar(0), fList2KStar(0), fList1KStar(0), fList0KStar(0), fListHelicityCuts(0), fList2KStarDiagnostic(0), fList1KStarDiagnostic(0), fList0KStarDiagnostic(0),
	fList2KStarEtaC(0), fList1KStarEtaC(0), fList0KStarEtaC(0), fHistNeventsEtaCKstarChannel(0), fKstarEventCandidates(0),
	fHistNeventsEtaC(0), fMPiKvsMPiK(0), fKstarMPiKvsMPiK(0),
	f2KstarEtaVsMinvEtaC(0), f2KstarEtaVsMinvEtaC400MeVPtMax(0), f2KstarEtaVsMinvEtaC100MeVPtMax(0), f2KstarSumPzVsMinvEtaC(0), f2KstarScalarSumP(0), f2KstarVectorSumPt(0), 
	f1KstarEtaVsMinvEtaC(0), f1KstarEtaVsMinvEtaC400MeVPtMax(0), f1KstarEtaVsMinvEtaC100MeVPtMax(0), f1KstarSumPzVsMinvEtaC(0), f1KstarScalarSumP(0), f1KstarVectorSumPt(0), 
	f0KstarEtaVsMinvEtaC(0), f0KstarEtaVsMinvEtaC400MeVPtMax(0), f0KstarEtaVsMinvEtaC100MeVPtMax(0), f0KstarSumPzVsMinvEtaC(0), f0KstarScalarSumP(0), f0KstarVectorSumPt(0),
	f2KstarPtPiPlus(0), f2KstarPtPiMinus(0), f2KstarPtKPlus(0), f2KstarPtKMinus(0), f2KstarTPCsignalPion(0), f2KstarTPCsignalKaon(0), f2KstarDedxVsPtPion(0), f2KstarDedxVsPtKaon(0), f2KstarTPCsignalVsQPtPion(0), f2KstarTPCsignalVsQPtKaon(0), f2KstarPtVsMinvFirstKstar(0), f2KstarPtVsMinvSecondKstar(0), f2KstarPtVsMinvEtaC(0),
	f1KstarPtPiPlus(0), f1KstarPtPiMinus(0), f1KstarPtKPlus(0), f1KstarPtKMinus(0), f1KstarTPCsignalPion(0), f1KstarTPCsignalKaon(0), f1KstarDedxVsPtPion(0), f1KstarDedxVsPtKaon(0), f1KstarTPCsignalVsQPtPion(0), f1KstarTPCsignalVsQPtKaon(0), f1KstarPtVsMinvKstar(0), f1KstarPtVsMinvOtherPiKcombo(0), f1KstarPtVsMinvEtaC(0),
	f0KstarPtPiPlus(0), f0KstarPtPiMinus(0), f0KstarPtKPlus(0), f0KstarPtKMinus(0), f0KstarTPCsignalPion(0), f0KstarTPCsignalKaon(0), f0KstarDedxVsPtPion(0), f0KstarDedxVsPtKaon(0), f0KstarTPCsignalVsQPtPion(0), f0KstarTPCsignalVsQPtKaon(0), f0KstarPtVsMinvFirstPiKcombo(0), f0KstarPtVsMinvSecondPiKcombo(0), f0KstarPtVsMinvEtaC(0),
	fKstarParentPx(0), fKstarParentPy(0), fKstarParentPz(0), fKstarDaughterParentAngle(0), fKstarDaughterParentCosAngle(0), fKstarDaughterDaughterAngle(0), fKstarDaughterDaughterCosAngle(0), fKstarDaughterPtotal(0), fKstarDaughterPtotalNorm(0),
	fKstarParentPxCheck(0), fKstarParentPyCheck(0), fKstarParentPzCheck(0), fKstarDaughterParentAngleCheck(0), fKstarDaughterParentCosAngleCheck(0), fKstarDaughterDaughterAngleCheck(0), fKstarDaughterDaughterCosAngleCheck(0), fKstarDaughterPtotalCheck(0), fKstarDaughterPtotalNormCheck(0),
	//2rho/4pi hists and lists
	fListHist2Rho4Pion(0), fList2RhoCandidates(0), fList2Rho(0), fList1Rho(0), fList0Rho(0), fList2RhoDiagnostic(0), fList1RhoDiagnostic(0), fList0RhoDiagnostic(0),
	fList2RhoEtaC(0), fList1RhoEtaC(0), fList0RhoEtaC(0), f4PiEventCandidates(0), f4PinPionVsnewPion(0),
	fHistNeventsEtaCRhoChannel(0), fHistNeventsFourPicharge(0), f2PairRhoM2VsPiPiM2(0), f1PairRhoM2VsPiPiM2(0), f2RhoM2VsPiPiM2(0),
	f2Rho2PairPtVsMinvEtaC(0), f2Rho2PairPtVsMinvFirstPair(0), f2Rho2PairPtVsMinvSecondPair(0), f2Rho2PairEtaVsMinvEtaC(0), f2Rho2PairEtaVsMinvEtaC400MeVPtMax(0), f2Rho2PairEtaVsMinvEtaC100MeVPtMax(0),
	f2Rho1PairPtVsMinvEtaC(0), f2Rho1PairPtVsMinvRhoPair(0), f2Rho1PairPtVsMinvNonRhoPair(0), f2Rho1PairEtaVsMinvEtaC(0), f2Rho1PairEtaVsMinvEtaC400MeVPtMax(0), f2Rho1PairEtaVsMinvEtaC100MeVPtMax(0),
	f4PionPtVsMinvEtaC(0), f4PionPtVsMinvRho(0), f4PionEtaVsMinvEtaC(0), f4PionEtaVsMinvEtaC400MeVPtMax(0), f4PionEtaVsMinvEtaC100MeVPtMax(0),
	f2Rho2PairSumPzVsMinvEtaC(0), f2Rho2PairScalarSumP(0), f2Rho2PairVectorSumPt(0),
	f2Rho1PairSumPzVsMinvEtaC(0), f2Rho1PairScalarSumP(0), f2Rho1PairVectorSumPt(0),
	f4PionSumPzVsMinvEtaC(0), f4PionScalarSumP(0), f4PionVectorSumPt(0),
	f2RhoParentPx(0), f2RhoParentPy(0), f2RhoParentPz(0), f2RhoDaughterParentAngle(0), f2RhoDaughterParentCosAngle(0), f2RhoDaughterDaughterAngle(0), f2RhoDaughterDaughterCosAngle(0), f2RhoDaughterPtotal(0),
	f2RhoParentPxCheck(0), f2RhoParentPyCheck(0), f2RhoParentPzCheck(0), f2RhoDaughterParentAngleCheck(0), f2RhoDaughterParentCosAngleCheck(0), f2RhoDaughterDaughterAngleCheck(0), f2RhoDaughterDaughterCosAngleCheck(0), f2RhoDaughterPtotalCheck(0),
	//3PiPi/4K/k0s hists and lists
	fListHistK0s3PiPi4K(0), fList3PiPi(0), fList4K(0), fListK0Short(0), fList3PiPiEtaC(0), fList4KEtaC(0), fListK0ShortEtaC(0), fList3PiPiDiagnostic(0), fList4KDiagnostic(0),
	fListK0ShortDiagnostic(0), fListK0ShortPID(0), f6PiEventCandidates(0),
	fHistNeventsEtaCK0sChannel(0), fHistK0sCandidatesPerEvent(0), fK0sPtVsMinvEtaC(0), fK0sEtaVsMinvEtaC(0), fK0sEtaVsMinvEtaC400MeVPtMax(0), fK0sEtaVsMinvEtaC100MeVPtMax(0),
	fK0sPosDaughterPt(0), fK0sNegDaughterPt(0), fK0sPosVsNegDaughterPt(0), fK0sPionPt(0), fK0sKaonPt(0), fK0sPtVsMinvK0s(0),
	fK0sPtVsMinvKPi(0), fK0sM2K0sVsM2KPi(0), fK0sM2K0sPiVsM2KPi(0), fK0sM2K0sKVsM2KPi(0), fK0sDecayLength(0),
	fHistFourTracksNpion(0), fHistSixTracksNpion(0), fHistNK0sPion(0), fHistNProngFound(0), fHistFourTracksNkaon(0), fHistSixTracksNkaon(0), fHistPostFourTracksNkaon(0), fHistPostFourTracksNpion(0), fK0sEventCandidates(0),
	fK0DaughterDca(0), fK0sDcaToPrimVertex(0), fK0sDaughterDcaToPrimVertex(0), fK0sMassDistribution(0), fV0DecayLength(0), fV0Eta(0), fV0CosPointingAngle(0), fV0sMassK0s(0), fK0GoodTracks(0),
	fK0sTOFbetaVsPAll(0), fK0sTOFbetaVsPPion(0), fK0sDedxVsPAll(0), fK0sDedxVsPPion(0),

	fHistNeventsEtaC3PiPiChannel(0), f3PiPiPtVsMinvEtaC(0), f3PiPiEtaVsMinvEtaC(0), f3PiPiEtaVsMinvEtaC400MeVPtMax(0), f3PiPiEtaVsMinvEtaC100MeVPtMax(0),
	fHistNeventsEtaC4KaonChannel(0), f4KEventCandidates(0), f4KaonPtVsMinvEtaC(0), f4KaonEtaVsMinvEtaC(0), f4KaonEtaVsMinvEtaC400MeVPtMax(0), f4KaonEtaVsMinvEtaC100MeVPtMax(0),
	f4KaonPtVsMinvKK(0), f4KVs2KMinv(0), f4KVs2KMinvSquared(0), fM2KKVsM2KK(0),
	f3PiPiSumPzVsMinvEtaC(0), f3PiPiScalarSumP(0), f3PiPiVectorSumPt(0), 
	f4KaonSumPzVsMinvEtaC(0), f4KaonScalarSumP(0), f4KaonVectorSumPt(0),
	fK0sSumPzVsMinvEtaC(0), fK0sScalarSumP(0), fK0sVectorSumPt(0),
	
	fList2K4Pi(0), fList2K4PiEtaC(0), fList2K4PiDiagnostic(0), fHistNeventsEtaC2K4PiChannel(0), f2K4PiEventCandidates(0), f2K4PiPtVsMinvEtaC(0), f2K4PiEtaVsMinvEtaC(0), f2K4PiEtaVsMinvEtaC400MeVPtMax(0),
	f2K4PiEtaVsMinvEtaC100MeVPtMax(0), f2K4PiSumPzVsMinvEtaC(0), f2K4PiScalarSumP(0), f2K4PiVectorSumPt(0)
{

//Dummy constructor

}//AliAnalysisTaskUpcFourPi

//_____________________________________________________________________________
AliAnalysisTaskUpcFourPi::AliAnalysisTaskUpcFourPi(const char *name) //second constructor, contains all histograms
	: AliAnalysisTaskSE(name), fType(0), fTracking(0), isMC(kFALSE), fRunTree(kFALSE), fRunHist(kTRUE), fRunSystematics(kFALSE), fPIDResponse(0), fEtaCK0sChannelTree(0), fEtaCTree(0), fMeritCutChoice(0),
	fRunNum(0), fPerNum(0), fOrbNum(0), fL0inputs(0), fL1inputs(0),
	fTOFmask(0), fIsPhysicsSelected(kFALSE),
	fVtxContrib(0), fVtxChi2(0), fVtxNDF(0), fSpdVtxContrib(0),
	fBCrossNum(0), fNtracklets(0), fNLooseTracks(0),
	fZNAenergy(0), fZNCenergy(0), fZPAenergy(0), fZPCenergy(0), fZDCAtime(0), fZDCCtime(0), fV0Adecision(0), fV0Cdecision(0), fADAdecision(0), fADCdecision(0),
	fDataFilnam(0), fRecoPass(0), fEvtNum(0),
	fJPsiAODTracks(0), fJPsiESDTracks(0), fEtaCAODTracks(0), fEtaCESDTracks(0), fGenPart(0),
	fEtaCCandidatesPerChannel(0), fEtaCLowPtCandidatesPerChannel(0), fAllPtVsMinvEtaC(0), fAllMinvEtaCLowPt(0), fChannelVsMinvEtaC(0),

	//trigger hists and lists
	fListTrig(0), fHistCcup29TriggersPerRun(0), fHistCcup30TriggersPerRun(0), fHistCcup31TriggersPerRun(0),
	fHistZedTriggersPerRun(0), fHistMBTriggersPerRun(0), fHistCentralTriggersPerRun(0), fHistSemiCentralTriggersPerRun(0), fHistCtrueTriggersPerRun(0),

	//PId/ZDC hists and lists
	fListHist(0), fListHistPID(0), fListHistUnusedPID(0), fListHistPIDcuts(0), fHistNTracks(0),
	fListHistZDC(0), fHistZDCCuts(0),
	fPionTPCvsPionTOFLowP(0), fPionTPCvsPionTOFMidP(0), fPionTPCvsPionTOFHighP(0),
	fKaonTPCvsKaonTOFLowP(0), fKaonTPCvsKaonTOFMidP(0), fKaonTPCvsKaonTOFHighP(0),
	fPionTPCvsKaonTPCLowP(0), fPionTPCvsKaonTPCMidP(0), fPionTPCvsKaonTPCHighP(0), fPionTOFvsKaonTOFMidP(0), fPionTOFvsKaonTOFPionMidP(0), fPionTOFvsKaonTOFKaonMidP(0),
	fDedxVsPAll(0), fDedxVsPPion(0), fDedxVsPKaon(0), fDedxVsPBoth(0), fHistPIDCheck(0),
	fTPCdEdxVsTOFbetaAll(0), fTPCdEdxVsTOFbetaPionsWithPID(0), fTPCdEdxVsTOFbetaKaonsWithPID(0),
	fFourTrackMissing(0), fSixTrackMissing(0),
	fTOFbetaVsPAll(0), fTOFbetaVsPPion(0), fTOFbetaVsPKaon(0), fTOFbetaVsPBoth(0),
	fHistZDCAenergy(0), fHistZDCCenergy(0), fHistZDCAtime(0), fHistZDCCtime(0), fHistZDCImpactParameter(0), fHistZDCAImpactParameter(0), fHistZDCCImpactParameter(0),
	fListSystematics(0), fListJPsiLoose(0), fListJPsiTight(0), fListEtaCLoose(0), fListEtaCTight(0),
	fHistFourTracksNboth(0), fHistFourTracksNneither(0), fHistSixTracksNboth(0), fHistSixTracksNneither(0),
	fListNKaonChange(0), fListNPionChange(0), f4KnKaonVsnewKaon(0),
	fList4TrackPID(0), fList6TrackPID(0), fTOFIntegratedLength(0),
	//Kstar hists and lists
	fListHistKstar(0), fList2KStar(0), fList1KStar(0), fList0KStar(0), fListHelicityCuts(0), fList2KStarDiagnostic(0), fList1KStarDiagnostic(0), fList0KStarDiagnostic(0),
	fList2KStarEtaC(0), fList1KStarEtaC(0), fList0KStarEtaC(0), fHistNeventsEtaCKstarChannel(0), fKstarEventCandidates(0),
	fHistNeventsEtaC(0), fMPiKvsMPiK(0), fKstarMPiKvsMPiK(0),
	f2KstarEtaVsMinvEtaC(0), f2KstarEtaVsMinvEtaC400MeVPtMax(0), f2KstarEtaVsMinvEtaC100MeVPtMax(0), f2KstarSumPzVsMinvEtaC(0), f2KstarScalarSumP(0), f2KstarVectorSumPt(0),
	f1KstarEtaVsMinvEtaC(0), f1KstarEtaVsMinvEtaC400MeVPtMax(0), f1KstarEtaVsMinvEtaC100MeVPtMax(0), f1KstarSumPzVsMinvEtaC(0), f1KstarScalarSumP(0), f1KstarVectorSumPt(0),
	f0KstarEtaVsMinvEtaC(0), f0KstarEtaVsMinvEtaC400MeVPtMax(0), f0KstarEtaVsMinvEtaC100MeVPtMax(0), f0KstarSumPzVsMinvEtaC(0), f0KstarScalarSumP(0), f0KstarVectorSumPt(0),
	f2KstarPtPiPlus(0), f2KstarPtPiMinus(0), f2KstarPtKPlus(0), f2KstarPtKMinus(0), f2KstarTPCsignalPion(0), f2KstarTPCsignalKaon(0), f2KstarDedxVsPtPion(0), f2KstarDedxVsPtKaon(0), f2KstarTPCsignalVsQPtPion(0), f2KstarTPCsignalVsQPtKaon(0), f2KstarPtVsMinvFirstKstar(0), f2KstarPtVsMinvSecondKstar(0), f2KstarPtVsMinvEtaC(0),
	f1KstarPtPiPlus(0), f1KstarPtPiMinus(0), f1KstarPtKPlus(0), f1KstarPtKMinus(0), f1KstarTPCsignalPion(0), f1KstarTPCsignalKaon(0), f1KstarDedxVsPtPion(0), f1KstarDedxVsPtKaon(0), f1KstarTPCsignalVsQPtPion(0), f1KstarTPCsignalVsQPtKaon(0), f1KstarPtVsMinvKstar(0), f1KstarPtVsMinvOtherPiKcombo(0), f1KstarPtVsMinvEtaC(0),
	f0KstarPtPiPlus(0), f0KstarPtPiMinus(0), f0KstarPtKPlus(0), f0KstarPtKMinus(0), f0KstarTPCsignalPion(0), f0KstarTPCsignalKaon(0), f0KstarDedxVsPtPion(0), f0KstarDedxVsPtKaon(0), f0KstarTPCsignalVsQPtPion(0), f0KstarTPCsignalVsQPtKaon(0), f0KstarPtVsMinvFirstPiKcombo(0), f0KstarPtVsMinvSecondPiKcombo(0), f0KstarPtVsMinvEtaC(0),
	fKstarParentPx(0), fKstarParentPy(0), fKstarParentPz(0), fKstarDaughterParentAngle(0), fKstarDaughterParentCosAngle(0), fKstarDaughterDaughterAngle(0), fKstarDaughterDaughterCosAngle(0), fKstarDaughterPtotal(0), fKstarDaughterPtotalNorm(0),
	fKstarParentPxCheck(0), fKstarParentPyCheck(0), fKstarParentPzCheck(0), fKstarDaughterParentAngleCheck(0), fKstarDaughterParentCosAngleCheck(0), fKstarDaughterDaughterAngleCheck(0), fKstarDaughterDaughterCosAngleCheck(0), fKstarDaughterPtotalCheck(0), fKstarDaughterPtotalNormCheck(0),
	//2rho/4pi hists and lists
	fListHist2Rho4Pion(0), fList2RhoCandidates(0), fList2Rho(0), fList1Rho(0), fList0Rho(0), fList2RhoDiagnostic(0), fList1RhoDiagnostic(0), fList0RhoDiagnostic(0),
	fList2RhoEtaC(0), fList1RhoEtaC(0), fList0RhoEtaC(0), f4PiEventCandidates(0), f4PinPionVsnewPion(0),
	fHistNeventsEtaCRhoChannel(0), fHistNeventsFourPicharge(0), f2PairRhoM2VsPiPiM2(0), f1PairRhoM2VsPiPiM2(0), f2RhoM2VsPiPiM2(0),
	f2Rho2PairPtVsMinvEtaC(0), f2Rho2PairPtVsMinvFirstPair(0), f2Rho2PairPtVsMinvSecondPair(0), f2Rho2PairEtaVsMinvEtaC(0), f2Rho2PairEtaVsMinvEtaC400MeVPtMax(0), f2Rho2PairEtaVsMinvEtaC100MeVPtMax(0),
	f2Rho1PairPtVsMinvEtaC(0), f2Rho1PairPtVsMinvRhoPair(0), f2Rho1PairPtVsMinvNonRhoPair(0), f2Rho1PairEtaVsMinvEtaC(0), f2Rho1PairEtaVsMinvEtaC400MeVPtMax(0), f2Rho1PairEtaVsMinvEtaC100MeVPtMax(0),
	f4PionPtVsMinvEtaC(0), f4PionPtVsMinvRho(0), f4PionEtaVsMinvEtaC(0), f4PionEtaVsMinvEtaC400MeVPtMax(0), f4PionEtaVsMinvEtaC100MeVPtMax(0),
	f2Rho2PairSumPzVsMinvEtaC(0), f2Rho2PairScalarSumP(0), f2Rho2PairVectorSumPt(0),
	f2Rho1PairSumPzVsMinvEtaC(0), f2Rho1PairScalarSumP(0), f2Rho1PairVectorSumPt(0),
	f4PionSumPzVsMinvEtaC(0), f4PionScalarSumP(0), f4PionVectorSumPt(0),
	f2RhoParentPx(0), f2RhoParentPy(0), f2RhoParentPz(0), f2RhoDaughterParentAngle(0), f2RhoDaughterParentCosAngle(0), f2RhoDaughterDaughterAngle(0), f2RhoDaughterDaughterCosAngle(0), f2RhoDaughterPtotal(0),
	f2RhoParentPxCheck(0), f2RhoParentPyCheck(0), f2RhoParentPzCheck(0), f2RhoDaughterParentAngleCheck(0), f2RhoDaughterParentCosAngleCheck(0), f2RhoDaughterDaughterAngleCheck(0), f2RhoDaughterDaughterCosAngleCheck(0), f2RhoDaughterPtotalCheck(0),
	//3PiPi/4K/k0s hists and lists
	fListHistK0s3PiPi4K(0), fList3PiPi(0), fList4K(0), fListK0Short(0), fList3PiPiEtaC(0), fList4KEtaC(0), fListK0ShortEtaC(0), fList3PiPiDiagnostic(0), fList4KDiagnostic(0),
	fListK0ShortDiagnostic(0), fListK0ShortPID(0), f6PiEventCandidates(0),
	fHistNeventsEtaCK0sChannel(0), fHistK0sCandidatesPerEvent(0), fK0sPtVsMinvEtaC(0), fK0sEtaVsMinvEtaC(0), fK0sEtaVsMinvEtaC400MeVPtMax(0), fK0sEtaVsMinvEtaC100MeVPtMax(0),
	fK0sPosDaughterPt(0), fK0sNegDaughterPt(0), fK0sPosVsNegDaughterPt(0), fK0sPionPt(0), fK0sKaonPt(0), fK0sPtVsMinvK0s(0),
	fK0sPtVsMinvKPi(0), fK0sM2K0sVsM2KPi(0), fK0sM2K0sPiVsM2KPi(0), fK0sM2K0sKVsM2KPi(0), fK0sDecayLength(0),
	fHistFourTracksNpion(0), fHistSixTracksNpion(0), fHistNK0sPion(0), fHistNProngFound(0), fHistFourTracksNkaon(0), fHistSixTracksNkaon(0), fHistPostFourTracksNkaon(0), fHistPostFourTracksNpion(0), fK0sEventCandidates(0),
	fK0DaughterDca(0), fK0sDcaToPrimVertex(0), fK0sDaughterDcaToPrimVertex(0), fK0sMassDistribution(0), fV0DecayLength(0), fV0Eta(0), fV0CosPointingAngle(0), fV0sMassK0s(0), fK0GoodTracks(0),
	fK0sTOFbetaVsPAll(0), fK0sTOFbetaVsPPion(0), fK0sDedxVsPAll(0), fK0sDedxVsPPion(0),

	fHistNeventsEtaC3PiPiChannel(0), f3PiPiPtVsMinvEtaC(0), f3PiPiEtaVsMinvEtaC(0), f3PiPiEtaVsMinvEtaC400MeVPtMax(0), f3PiPiEtaVsMinvEtaC100MeVPtMax(0),
	fHistNeventsEtaC4KaonChannel(0), f4KEventCandidates(0), f4KaonPtVsMinvEtaC(0), f4KaonEtaVsMinvEtaC(0), f4KaonEtaVsMinvEtaC400MeVPtMax(0), f4KaonEtaVsMinvEtaC100MeVPtMax(0),
	f4KaonPtVsMinvKK(0), f4KVs2KMinv(0), f4KVs2KMinvSquared(0), fM2KKVsM2KK(0),
	f3PiPiSumPzVsMinvEtaC(0), f3PiPiScalarSumP(0), f3PiPiVectorSumPt(0),
	f4KaonSumPzVsMinvEtaC(0), f4KaonScalarSumP(0), f4KaonVectorSumPt(0),
	fK0sSumPzVsMinvEtaC(0), fK0sScalarSumP(0), fK0sVectorSumPt(0),

	fList2K4Pi(0), fList2K4PiEtaC(0), fList2K4PiDiagnostic(0), fHistNeventsEtaC2K4PiChannel(0), f2K4PiEventCandidates(0), f2K4PiPtVsMinvEtaC(0), f2K4PiEtaVsMinvEtaC(0), f2K4PiEtaVsMinvEtaC400MeVPtMax(0),
	f2K4PiEtaVsMinvEtaC100MeVPtMax(0), f2K4PiSumPzVsMinvEtaC(0), f2K4PiScalarSumP(0), f2K4PiVectorSumPt(0)
{

  // Constructor
  if( strstr(name,"ESD") ) fType = 0;
  if( strstr(name,"AOD") ) fType = 1;

  //cout << "##### fType = " << fType << " so AOD selected." << endl;

  fMeritCutChoice = 4;  //case for selecting best K0s candidate
  
  Init();

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}//AliAnalysisTaskUpcFourPi

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::Init()
{
  
  for(Int_t i=0; i<ntrg; i++) { //records which triggers were triggered on each event
  	fTrigger[i] = kFALSE;
	fTriggerInputsMC[i] = kFALSE;
	}
  for(Int_t i=0; i<7; i++) { //arrays storing up to 7 values for up to 7 tracks of the TPC's PID sigma values for each particle type
	fPIDTPCMuon[i] = -666;
	fPIDTPCElectron[i] = -666;
	fPIDTPCPion[i] = -666;
	fPIDTPCKaon[i] = -666;
	fPIDTPCProton[i] = -666;
	
	fPIDTOFMuon[i] = -666;
	fPIDTOFElectron[i] = -666;
	fPIDTOFPion[i] = -666;
	fPIDTOFKaon[i] = -666;
	fPIDTOFProton[i] = -666;
	
	fIsVtxContributor[i] = kFALSE;
	}
  for(Int_t i=0; i<7; i++) { //These are being used for K0 channel
	fPIDTPCMuonPos[i] = -666;
	fPIDTPCElectronPos[i] = -666;
	fPIDTPCPionPos[i] = -666;
	fPIDTPCKaonPos[i] = -666;
	fPIDTPCProtonPos[i] = -666;
	
	fPIDTOFMuonPos[i] = -666;
	fPIDTOFElectronPos[i] = -666;
	fPIDTOFPionPos[i] = -666;
	fPIDTOFKaonPos[i] = -666;
	fPIDTOFProtonPos[i] = -666;

	fPIDTPCMuonNeg[i] = -666;
	fPIDTPCElectronNeg[i] = -666;
	fPIDTPCPionNeg[i] = -666;
	fPIDTPCKaonNeg[i] = -666;
	fPIDTPCProtonNeg[i] = -666;
	
	fPIDTOFMuonNeg[i] = -666;
	fPIDTOFElectronNeg[i] = -666;
	fPIDTOFPionNeg[i] = -666;
	fPIDTOFKaonNeg[i] = -666;
	fPIDTOFProtonNeg[i] = -666;
  }
  for(Int_t i=0; i<3; i++){ //vertex positions and errors, used in trees only
  	fVtxPos[i] = -666; 
	fMCVtxPos[i] = -666;
	fVtxErr[i] = -666;
	fKfVtxPos[i] = -666;
	fSpdVtxPos[i] = -666;
	}

  //cout << "##### end of Init()" << endl;
}//Init

//_____________________________________________________________________________
AliAnalysisTaskUpcFourPi::~AliAnalysisTaskUpcFourPi() 
{
  // Destructor
	if (fListTrig) {
		delete fListTrig;
		fListTrig = 0x0;
	}
	if (fListHist) {
		delete fListHist;
		fListHist = 0x0;
	}
	if (fListHistKstar) {
		delete fListHistKstar;
		fListHistKstar = 0x0;
	}
	if (fListHist2Rho4Pion) {
		delete fListHist2Rho4Pion;
		fListHist2Rho4Pion = 0x0;
	}
	if (fListHistK0s3PiPi4K) {
		delete fListHistK0s3PiPi4K;
		fListHistK0s3PiPi4K = 0x0;
	}
}//~AliAnalysisTaskUpcFourPi

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::UserCreateOutputObjects() //use the names defined in the constructor to makes specific data structures
{
	//cout << "##### Start of UserCreateOutputObjects()" << endl;

	//PID response
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();

	//input file
	fDataFilnam = new TObjString();
	fDataFilnam->SetString("");

	//##### Create Track Array objects. #####
	fEtaCAODTracks = new TClonesArray("AliAODTrack", 1000);
	fEtaCESDTracks = new TClonesArray("AliESDtrack", 1000);
	fGenPart = new TClonesArray("TParticle", 1000);

	//##### Create an Output Container and trigger class histograms. #####
	//NOTE: These histograms count the actual number of events analyzed from each trigger class in each run.
	//      I think this information must be compared to the number recorded/produced per run and luminosity
	//      per run to get the actual luminosity analyzed in your analysis. Something like 
	//                    Lum_analyzed = (N_analyzed/N_recorded) * Lum_recorded.

	//##### Create some Output Containers. #####
	fListTrig = new TList();
	fListTrig->SetOwner();

	fListHist = new TList();
	fListHist->SetOwner();

	fListHistKstar = new TList();
	fListHistKstar->SetOwner();

	fListHist2Rho4Pion = new TList();
	fListHist2Rho4Pion->SetOwner();

	fListHistK0s3PiPi4K = new TList();
	fListHistK0s3PiPi4K->SetOwner();

	//ListHist sub-lists
	fListHistZDC = new TList();
	fListHistZDC->SetOwner();
	fListHistZDC->SetName("ZDC Histograms");
	fListHist->Add(fListHistZDC);

	fListHistPID = new TList();
	fListHistPID->SetOwner();
	fListHistPID->SetName("PID Histograms");
	fListHist->Add(fListHistPID);

	fListHistPIDcuts = new TList();
	fListHistPIDcuts->SetOwner();
	fListHistPIDcuts->SetName("PID cuts");
	fListHist->Add(fListHistPIDcuts);

	fListHistUnusedPID = new TList();
	fListHistUnusedPID->SetOwner();
	fListHistUnusedPID->SetName("TPC and TOF Signals");
	fListHistPIDcuts->Add(fListHistUnusedPID);

	fListSystematics = new TList();
	fListSystematics->SetOwner();
	fListSystematics->SetName("Systematics");
	fListHist->Add(fListSystematics);

	//Kstar sub-lists
	fList2KStar = new TList();
	fList2KStar->SetOwner();
	fList2KStar->SetName("Two KStar");
	fListHistKstar->Add(fList2KStar);

	fList1KStar = new TList();
	fList1KStar->SetOwner();
	fList1KStar->SetName("One KStar");
	fListHistKstar->Add(fList1KStar);

	fList0KStar = new TList();
	fList0KStar->SetOwner();
	fList0KStar->SetName("Zero KStar");
	fListHistKstar->Add(fList0KStar);

	fListHelicityCuts = new TList();
	fListHelicityCuts->SetOwner();
	fListHelicityCuts->SetName("Helicity Cuts");
	fListHistKstar->Add(fListHelicityCuts);

	fList2KStarDiagnostic = new TList();
	fList2KStarDiagnostic->SetOwner();
	fList2KStarDiagnostic->SetName("Diagnostic");
	fList2KStar->Add(fList2KStarDiagnostic);

	fList1KStarDiagnostic = new TList();
	fList1KStarDiagnostic->SetOwner();
	fList1KStarDiagnostic->SetName("Diagnostic");
	fList1KStar->Add(fList1KStarDiagnostic);

	fList0KStarDiagnostic = new TList();
	fList0KStarDiagnostic->SetOwner();
	fList0KStarDiagnostic->SetName("Diagnostic");
	fList0KStar->Add(fList0KStarDiagnostic);

	fList2KStarEtaC = new TList();
	fList2KStarEtaC->SetOwner();
	fList2KStarEtaC->SetName("EtaC");
	fList2KStar->Add(fList2KStarEtaC);

	fList1KStarEtaC = new TList();
	fList1KStarEtaC->SetOwner();
	fList1KStarEtaC->SetName("EtaC");
	fList1KStar->Add(fList1KStarEtaC);

	fList0KStarEtaC = new TList();
	fList0KStarEtaC->SetOwner();
	fList0KStarEtaC->SetName("EtaC");
	fList0KStar->Add(fList0KStarEtaC);

	//4pi sub-lists
	fList2Rho = new TList();
	fList2Rho->SetOwner();
	fList2Rho->SetName("Two Rho");
	fListHist2Rho4Pion->Add(fList2Rho);

	fList1Rho = new TList();
	fList1Rho->SetOwner();
	fList1Rho->SetName("One Rho");
	fListHist2Rho4Pion->Add(fList1Rho);

	fList0Rho = new TList();
	fList0Rho->SetOwner();
	fList0Rho->SetName("Zero Rho");
	fListHist2Rho4Pion->Add(fList0Rho);

	fList2RhoCandidates = new TList();
	fList2RhoCandidates->SetOwner();
	fList2RhoCandidates->SetName("Rho Candidates");
	fListHist2Rho4Pion->Add(fList2RhoCandidates);

	fList2RhoDiagnostic = new TList();
	fList2RhoDiagnostic->SetOwner();
	fList2RhoDiagnostic->SetName("Diagnostic");
	fList2Rho->Add(fList2RhoDiagnostic);

	fList1RhoDiagnostic = new TList();
	fList1RhoDiagnostic->SetOwner();
	fList1RhoDiagnostic->SetName("Diagnostic");
	fList1Rho->Add(fList1RhoDiagnostic);

	fList0RhoDiagnostic = new TList();
	fList0RhoDiagnostic->SetOwner();
	fList0RhoDiagnostic->SetName("Diagnostic");
	fList0Rho->Add(fList0RhoDiagnostic);

	fList2RhoEtaC = new TList();
	fList2RhoEtaC->SetOwner();
	fList2RhoEtaC->SetName("EtaC");
	fList2Rho->Add(fList2RhoEtaC);

	fList1RhoEtaC = new TList();
	fList1RhoEtaC->SetOwner();
	fList1RhoEtaC->SetName("EtaC");
	fList1Rho->Add(fList1RhoEtaC);

	fList0RhoEtaC = new TList();
	fList0RhoEtaC->SetOwner();
	fList0RhoEtaC->SetName("EtaC");
	fList0Rho->Add(fList0RhoEtaC);

	//3pipi sub-lists
	fList3PiPi = new TList();
	fList3PiPi->SetOwner();
	fList3PiPi->SetName("Three Pi Pi");
	fListHistK0s3PiPi4K->Add(fList3PiPi);

	fList4K = new TList();
	fList4K->SetOwner();
	fList4K->SetName("Four Kaon");
	fListHistK0s3PiPi4K->Add(fList4K);

	fListK0Short = new TList();
	fListK0Short->SetOwner();
	fListK0Short->SetName("K0 Short");
	fListHistK0s3PiPi4K->Add(fListK0Short);

	fList2K4Pi = new TList();
	fList2K4Pi->SetOwner();
	fList2K4Pi->SetName("Two K Four Pi");
	fListHistK0s3PiPi4K->Add(fList2K4Pi);

	fList3PiPiEtaC = new TList();
	fList3PiPiEtaC->SetOwner();
	fList3PiPiEtaC->SetName("EtaC");
	fList3PiPi->Add(fList3PiPiEtaC);

	fList4KEtaC = new TList();
	fList4KEtaC->SetOwner();
	fList4KEtaC->SetName("EtaC");
	fList4K->Add(fList4KEtaC);

	fListK0ShortEtaC = new TList();
	fListK0ShortEtaC->SetOwner();
	fListK0ShortEtaC->SetName("EtaC");
	fListK0Short->Add(fListK0ShortEtaC);

	fList2K4PiEtaC = new TList();
	fList2K4PiEtaC->SetOwner();
	fList2K4PiEtaC->SetName("EtaC");
	fList2K4Pi->Add(fList2K4PiEtaC);

	fList3PiPiDiagnostic = new TList();
	fList3PiPiDiagnostic->SetOwner();
	fList3PiPiDiagnostic->SetName("Diagnostic");
	fList3PiPi->Add(fList3PiPiDiagnostic);

	fList4KDiagnostic = new TList();
	fList4KDiagnostic->SetOwner();
	fList4KDiagnostic->SetName("Diagnostic");
	fList4K->Add(fList4KDiagnostic);

	fListK0ShortDiagnostic = new TList();
	fListK0ShortDiagnostic->SetOwner();
	fListK0ShortDiagnostic->SetName("Diagnostic");
	fListK0Short->Add(fListK0ShortDiagnostic);

	fList2K4PiDiagnostic = new TList();
	fList2K4PiDiagnostic->SetOwner();
	fList2K4PiDiagnostic->SetName("Diagnostic");
	fList2K4Pi->Add(fList2K4PiDiagnostic);

	fListK0ShortPID = new TList();
	fListK0ShortPID->SetOwner();
	fListK0ShortPID->SetName("PID Checks");
	fListK0Short->Add(fListK0ShortPID);

	fList4TrackPID = new TList();
	fList4TrackPID->SetOwner();
	fList4TrackPID->SetName("Four Track PID");
	fListHistPID->Add(fList4TrackPID);

	fList6TrackPID = new TList();
	fList6TrackPID->SetOwner();
	fList6TrackPID->SetName("Six Track PID");
	fListHistPID->Add(fList6TrackPID);

	//##### Define trigger histograms. #####
	fHistCcup29TriggersPerRun = new TH1D("fHistCcup29TriggersPerRun", "fHistCcup29TriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistCcup29TriggersPerRun); 
	fHistCcup30TriggersPerRun = new TH1D("fHistCcup30TriggersPerRun", "fHistCcup30TriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistCcup30TriggersPerRun);
	fHistCcup31TriggersPerRun = new TH1D("fHistCcup31TriggersPerRun", "fHistCcup31TriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistCcup31TriggersPerRun);

	fHistZedTriggersPerRun = new TH1D("fHistZedTriggersPerRun", "fHistZedTriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistZedTriggersPerRun);
	fHistMBTriggersPerRun = new TH1D("fHistMBTriggersPerRun", "fHistMBTriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistMBTriggersPerRun);
	fHistCentralTriggersPerRun = new TH1D("fHistCentralTriggersPerRun", "fHistCentralTriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistCentralTriggersPerRun);
	fHistSemiCentralTriggersPerRun = new TH1D("fHistSemiCentralTriggersPerRun", "fHistSemiCentralTriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistSemiCentralTriggersPerRun);
	fHistCtrueTriggersPerRun = new TH1D("fHistCtrueTriggersPerRun", "fHistCtrueTriggersPerRun", 935, 296689.5, 297624.5);
	fListTrig->Add(fHistCtrueTriggersPerRun);

	//##### Define histograms for ZDC information. #####
	TString CutNameZDC[4] = { "CCUP4","< 8 neutrons","0 netrons","No timing" };
	fHistZDCCuts = new TH1D("fHistZDCCuts", "fHistZDCCuts", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) fHistZDCCuts->GetXaxis()->SetBinLabel(i + 1, CutNameZDC[i].Data());
	fListHistZDC->Add(fHistZDCCuts);

	fHistZDCAenergy = new TH1D("fHistZDCAenergy", "fHistZDCAenergy", 1500, -3000., 12000.);
	fListHistZDC->Add(fHistZDCAenergy);
	fHistZDCCenergy = new TH1D("fHistZDCCenergy", "fHistZDCCenergy", 1500, -3000., 12000.);
	fListHistZDC->Add(fHistZDCCenergy);
	fHistZDCAtime = new TH1D("fHistZDCAtime", "fHistZDCAtime", 1600, -8., 8.);
	fListHistZDC->Add(fHistZDCAtime);
	fHistZDCCtime = new TH1D("fHistZDCCtime", "fHistZDCCtime", 1600, -8., 8.);
	fListHistZDC->Add(fHistZDCCtime);
	fHistZDCImpactParameter = new TH1D("fHistZDCImpactParameter", "fHistZDCImpactParameter", 100., 0., 10000.); //540,0.,1000); //Test large |b| values for UPC.
	fListHistZDC->Add(fHistZDCImpactParameter);                                                  //Test |b| in femtometers (units in fm).
	fHistZDCAImpactParameter = new TH1D("fHistZDCAImpactParameter", "fHistZDCCImpactParameter", 100., 0., 100.*pow(10, -15)); //540,-TMath::Pi(),2.*TMath::Pi());
	fListHistZDC->Add(fHistZDCAImpactParameter);                                                 //Test |b| in hundreds of femtometers (units in m).
	fHistZDCCImpactParameter = new TH1D("fHistZDCCImpactParameter", "fHistZDCCImpactparameter", 100., 0., 1.*pow(10, -15));//540,-TMath::Pi(),2.*TMath::Pi());
	fListHistZDC->Add(fHistZDCCImpactParameter);                                                 //Test |b| in femtometers (units in m).

	//##### Define Histograms for PID #####	
	fPionTPCvsPionTOFLowP = new TH2D("fPionTPCvsPionTOFLowP", "fPionTPCvsPionTOFLowP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsPionTOFLowP);
	fKaonTPCvsKaonTOFLowP = new TH2D("fKaonTPCvsKaonTOFLowP", "fKaonTPCvsKaonTOFLowP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fKaonTPCvsKaonTOFLowP);
	fPionTPCvsKaonTPCLowP = new TH2D("fPionTPCvsKaonTPCLowP", "fPionTPCvsKaonTPCLowP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsKaonTPCLowP);

	fPionTPCvsPionTOFMidP = new TH2D("fPionTPCvsPionTOFMidP", "fPionTPCvsPionTOFMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsPionTOFMidP);
	fKaonTPCvsKaonTOFMidP = new TH2D("fKaonTPCvsKaonTOFMidP", "fKaonTPCvsKaonTOFMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fKaonTPCvsKaonTOFMidP);
	fPionTPCvsKaonTPCMidP = new TH2D("fPionTPCvsKaonTPCMidP", "fPionTPCvsKaonTPCMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsKaonTPCMidP);
	fPionTOFvsKaonTOFMidP = new TH2D("fPionTOFvsKaonTOFMidP", "fPionTOFvsKaonTOFMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTOFvsKaonTOFMidP);
	fPionTOFvsKaonTOFPionMidP = new TH2D("fPionTOFvsKaonTOFPionMidP", "fPionTOFvsKaonTOFPionMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTOFvsKaonTOFPionMidP);
	fPionTOFvsKaonTOFKaonMidP = new TH2D("fPionTOFvsKaonTOFKaonMidP", "fPionTOFvsKaonTOFKaonMidP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTOFvsKaonTOFKaonMidP);

	fPionTPCvsPionTOFHighP = new TH2D("fPionTPCvsPionTOFHighP", "fPionTPCvsPionTOFHighP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsPionTOFHighP);
	fKaonTPCvsKaonTOFHighP = new TH2D("fKaonTPCvsKaonTOFHighP", "fKaonTPCvsKaonTOFHighP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fKaonTPCvsKaonTOFHighP);
	fPionTPCvsKaonTPCHighP = new TH2D("fPionTPCvsKaonTPCHighP", "fPionTPCvsKaonTPCHighP", 900, -30, 60, 900, -30, 60);
	fListHistPIDcuts->Add(fPionTPCvsKaonTPCHighP);

	fHistNTracks = new TH1D("nTracks", "Number of Tracks per Event", 12, -0.5, 11.5);
	fListHistPID->Add(fHistNTracks);
	fHistFourTracksNpion = new TH1D("Four Track nPion", "nPion Per Four Track Event", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistFourTracksNpion);
	fHistFourTracksNkaon = new TH1D("Four Track nKaon", "nKaon Per Four Track Event", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistFourTracksNkaon);
	fHistFourTracksNboth = new TH1D("Four Track nBoth", "nBoth Per Four Track Event", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistFourTracksNboth);
	fFourTrackMissing = new TH1D("Four Track nMissingTOF", "nMissingTOF Per Four Track Event", 9, -0.5, 8.5);
	fList4TrackPID->Add(fFourTrackMissing);
	fHistFourTracksNneither = new TH1D("Four Track nNeither", "nNeither Per Four Track Event", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistFourTracksNneither);
	fHistPostFourTracksNpion = new TH1D("Four Track nPion Check", "Change in nPion (4 tracks) after analysis (should be 0)", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistPostFourTracksNpion);
	fHistPostFourTracksNkaon = new TH1D("Four Track nKaon Check", "Change in nKaon (4 tracks) after analysis (should be 0)", 9, -0.5, 8.5);
	fList4TrackPID->Add(fHistPostFourTracksNkaon);
	fHistSixTracksNpion = new TH1D("Four Track nPion", "nPion Per Four Track Event", 9, -0.5, 8.5);
	fList6TrackPID->Add(fHistSixTracksNpion);
	fHistSixTracksNkaon = new TH1D("Four Track nKaon", "nKaon Per Four Track Event", 9, -0.5, 8.5);
	fList6TrackPID->Add(fHistSixTracksNkaon);
	fHistSixTracksNboth = new TH1D("Four Track nBoth", "nBoth Per Four Track Event", 9, -0.5, 8.5);
	fList6TrackPID->Add(fHistSixTracksNboth);
	fSixTrackMissing = new TH1D("Four Track nMissingTOF", "nMissingTOF Per Four Track Event", 9, -0.5, 8.5);
	fList6TrackPID->Add(fSixTrackMissing);
	fHistSixTracksNneither = new TH1D("Four Track nNeither", "nNeither Per Four Track Event", 9, -0.5, 8.5);
	fList6TrackPID->Add(fHistSixTracksNneither);
	f4KnKaonVsnewKaon = new TH2D("Two KK Channel Missing TOF check", "4KK nKaon V. newKaon", 15, -0.5, 4.5, 15, -0.5, 4.5);
	fListHistPID->Add(f4KnKaonVsnewKaon);
	f4PinPionVsnewPion = new TH2D("Two PiPi Channel Missing TOF check", "4PiPi nPion V. newPion", 15, -0.5, 4.5, 15, -0.5, 4.5);
	fListHistPID->Add(f4PinPionVsnewPion);

	fHistPIDCheck = new TH1D("fHistPIDCheck", "fHistPIDCheck", 9, -0.5, 8.5);
	fListHistPID->Add(fHistPIDCheck);
	fListNPionChange = new TH1D("fListNPionChange", "fListNPionChange", 9, -4.5, 4.5);
	fListHistPID->Add(fListNPionChange);
	fListNKaonChange = new TH1D("fListNKaonChange", "fListNKaonChange", 9, -4.5, 4.5);
	fListHistPID->Add(fListNKaonChange);
	
	fDedxVsPAll = new TH2D("dEdx V P All", "TPC dE/dx v. Momentum, no PID", 1000, 0., 10., 3000, 0., 300.);
	fListHistUnusedPID->Add(fDedxVsPAll);
	fDedxVsPPion = new TH2D("dEdx V P Pion", "TPC dE/dx v. Momentum, Pion PID", 1000, 0., 10., 3000, 0., 300.);
	fListHistUnusedPID->Add(fDedxVsPPion);
	fDedxVsPKaon = new TH2D("dEdx V P Kaon", "TPC dE/dx v. Momentum, Kaon PID", 1000, 0., 10., 3000, 0., 300.);
	fListHistUnusedPID->Add(fDedxVsPKaon);
	fDedxVsPBoth = new TH2D("new dEdx V P Both", "TPC dE/dx v. Momentum, Both new PID", 1000, 0., 10., 3000, 0., 300.);
	fListHistUnusedPID->Add(fDedxVsPBoth);
	fTOFIntegratedLength = new TH1D("Integrated Track Length", "Integrated Track Length for non-missing tracks", 1000, 0, 1000);
	fListHistUnusedPID->Add(fTOFIntegratedLength);
	fTOFbetaVsPAll = new TH2D("beta V P All", "TOF #beta .v Momentum, No PID", 1000, 0., 10., 200, 0., 2.);
	fListHistUnusedPID->Add(fTOFbetaVsPAll);
	fTOFbetaVsPPion = new TH2D("beta v P Pions", "TOF #beta v. Momentum, Pion PID", 1000, 0., 10., 200, 0., 2.);
	fListHistUnusedPID->Add(fTOFbetaVsPPion);
	fTOFbetaVsPKaon = new TH2D("beta V P Kaons", "TOF #beta V. Momentum, Kaon PID", 1000, 0., 10., 200, 0., 2.);
	fListHistUnusedPID->Add(fTOFbetaVsPKaon);
	fTOFbetaVsPBoth = new TH2D("new beta V P Both", "TOF #beta v. Momentum, Both new PID", 1000, 0., 10., 200, 0., 2.);
	fListHistUnusedPID->Add(fTOFbetaVsPBoth);
	fTPCdEdxVsTOFbetaAll = new TH2D("fTPCdEdxVsTOFbetaAll", "fTPCdEdxVsTOFbetaAll", 200, 0., 2., 3000, 0., 300.);
	fListHistUnusedPID->Add(fTPCdEdxVsTOFbetaAll);
	fTPCdEdxVsTOFbetaPionsWithPID = new TH2D("fTPCdEdxVsTOFbetaPionsWithPID", "fTPCdEdxVsTOFbetaPionsWithPID", 200, 0., 2., 3000, 0., 300.);
	fListHistUnusedPID->Add(fTPCdEdxVsTOFbetaPionsWithPID);
	fTPCdEdxVsTOFbetaKaonsWithPID = new TH2D("fTPCdEdxVsTOFbetaKaonsWithPID", "fTPCdEdxVsTOFbetaKaonsWithPID", 200, 0., 2., 3000, 0., 300.);
	fListHistUnusedPID->Add(fTPCdEdxVsTOFbetaKaonsWithPID);

	TString CutNameEtaC[8] = { "Analyzed","Triggered","Vertex cut","V0 decision","Neutron ZDC cut","Four Tracks, >1 SPD hit","#eta_{C} Channel Match","#eta_{C} candidates" };
	fHistNeventsEtaC = new TH1D("EtaC Events Histogram", "#eta_{C} Events at Each Analysis Stage", 8, 0.5, 8.5);
	for (Int_t i = 0; i < 8; i++) fHistNeventsEtaC->GetXaxis()->SetBinLabel(i + 1, CutNameEtaC[i].Data());
	fListHist->Add(fHistNeventsEtaC);

	TString ChannelNames[9] = { "K* K*","K* K^{+} #pi^{-} + c.c.","K^{+} #pi^{-} K^{-} #pi^{+} + c.c.","#rho^{0} #rho^{0}","2(#pi^{+} #pi^{-})","3(#pi^{+} #pi^{-})",
		"2(K^{+} K^{-})","K^{0} K^{-} #pi^{+} #pi^{-} #pi^{+} + c.c.","K^{+} K^{-} 2(#pi^{+} #pi^{-})" };
	fEtaCCandidatesPerChannel = new TH1D("Candidates Per Channel", "#eta_{C} Candidates Per Decay Channel", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fEtaCCandidatesPerChannel->GetXaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fEtaCCandidatesPerChannel);

	fEtaCLowPtCandidatesPerChannel = new TH1D("Candidates Per Channel Low Pt", "#eta_{C} Candidates Per Decay Channel, Low Pt", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fEtaCLowPtCandidatesPerChannel->GetXaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fEtaCLowPtCandidatesPerChannel);

	fAllPtVsMinvEtaC = new TH2D("Pt V Minv EtaC Candidates", "Transverse Momentum V. Invariant Mass, all #eta_{C} candidates", 1200, 0., 6., 500, 0., 5.);
	fListHist->Add(fAllPtVsMinvEtaC);

	fAllMinvEtaCLowPt = new TH1D("Minv EtaC Candidates, Low Pt", "Invariant Mass, Low Pt #eta_{C} candidates", 300, 0., 6.);
	fListHist->Add(fAllMinvEtaCLowPt);

	//fAllMinvEtaCLowPt->Fit("gaus");

	fChannelVsMinvEtaC = new TH2D("Channel V Minv EtaC, Low Pt", "Decay Channel V. Invariant Mass, Low Pt #eta_{C} candidates", 300, 0., 6., 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fChannelVsMinvEtaC->GetYaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fChannelVsMinvEtaC);




	//##### Define Decay Channel histograms #####

	//Helicity Cut histograms
	fKstarParentPx = new TH1D("fKstarParentPx", "Parent p_{x} in rest frame (k*(892) channels", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPx);
	fKstarParentPy = new TH1D("fKstarParentPy", "Parent p_{y} in rest frame (k*(892) channels", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPy);
	fKstarParentPz = new TH1D("fKstarParentPz", "Parent p_{z} in rest frame (k*(892) channels", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPz);
	fKstarDaughterParentAngle = new TH1D("fKstarDaughterParentAngle", "Angle between daughter and parent", 3800, -190., 190.);
	fListHelicityCuts->Add(fKstarDaughterParentAngle);
	fKstarDaughterParentCosAngle = new TH1D("fKstarDaughterParentCosAngle", "Cosine of angle between daughter and parent", 2200, -1.1, 1.1);
	fListHelicityCuts->Add(fKstarDaughterParentCosAngle);
	fKstarDaughterDaughterAngle = new TH1D("fKstarDaughterDaughterAngle", "Angle between the two daughters", 2000, -10., 190.);
	fListHelicityCuts->Add(fKstarDaughterDaughterAngle);
	fKstarDaughterDaughterCosAngle = new TH1D("fKstarDaughterDaughterCosAngle", "Cos(Angle) between two daughters", 220, -1.1, 1.1);
	fListHelicityCuts->Add(fKstarDaughterDaughterCosAngle);
	fKstarDaughterPtotal = new TH1D("fKstarDaughterPtotal", "Momentum sum of two daughters", 1000, -5., 5.);
	fListHelicityCuts->Add(fKstarDaughterPtotal);
	fKstarDaughterPtotalNorm = new TH1D("fKstarDaughterPtotalNorm", "Normalized momentum sum of two daughters", 1000, -5., 5.);
	fListHelicityCuts->Add(fKstarDaughterPtotalNorm);

	//Helicity Cut histograms - Check histos
	fKstarParentPxCheck = new TH1D("fKstarParentPxCheck", "Parent p_{x} in rest frame (k*(892) channels (cuts passed)", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPxCheck);
	fKstarParentPyCheck = new TH1D("fKstarParentPyCheck", "Parent p_{y} in rest frame (k*(892) channels (cuts passed)", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPyCheck);
	fKstarParentPzCheck = new TH1D("fKstarParentPzCheck", "Parent p_{z} in rest frame (k*(892) channels (cuts passed)", 200, -1., 1.);
	fListHelicityCuts->Add(fKstarParentPzCheck);
	fKstarDaughterParentAngleCheck = new TH1D("fKstarDaughterParentAngleCheck", "Angle between daughter and parent (cuts passed)", 3800, -190., 190.);
	fListHelicityCuts->Add(fKstarDaughterParentAngleCheck);
	fKstarDaughterParentCosAngleCheck = new TH1D("fKstarDaughterParentCosAngleCheck", "Cosine of angle between daughter and parent (cuts passed)", 2200, -1.1, 1.1);
	fListHelicityCuts->Add(fKstarDaughterParentCosAngleCheck);
	fKstarDaughterDaughterAngleCheck = new TH1D("fKstarDaughterDaughterAngleCheck", "Angle between the two daughters (cuts passed)", 2000, -10., 190.);
	fListHelicityCuts->Add(fKstarDaughterDaughterAngleCheck);
	fKstarDaughterDaughterCosAngleCheck = new TH1D("fKstarDaughterDaughterCosAngleCheck", "Cos(Angle) between two daughters (cuts passed)", 220, -1.1, 1.1);
	fListHelicityCuts->Add(fKstarDaughterDaughterCosAngleCheck);
	fKstarDaughterPtotalCheck = new TH1D("fKstarDaughterPtotalCheck", "Momentum sum of two daughters (cuts passed)", 1000, -5., 5.);
	fListHelicityCuts->Add(fKstarDaughterPtotalCheck);
	fKstarDaughterPtotalNormCheck = new TH1D("fKstarDaughterPtotalNormCheck", "Normalized momentum sum of two daughters (cuts passed)", 1000, -5., 5.);
	fListHelicityCuts->Add(fKstarDaughterPtotalNormCheck);

	//##### Define 2 rho(0) channel histograms.
	TString CutNameEtaCRhoChannel[7] = { "Four Pions","non-zero net charge","zero net charge (#eta_{C})","two sets of #rho's",
		"one set of #rho's","no sets of #rho's","2+ tracks with pT>0.4" }; 

	TString CutNameFourPiCharge[2] = { "Four Pions","non-zero net charge"};
	
	//Create a histotram to count events after five event cuts and with different event/track combinations. 
	fHistNeventsEtaCRhoChannel = new TH1D("EtaC Events Rho Channel", "#eta_{C} Events at Each 2(#pi^{+} #pi^{-}) Channel Analysis Stage", 7, 0.5, 7.5);
	for (Int_t i = 0; i < 7; i++) fHistNeventsEtaCRhoChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaCRhoChannel[i].Data());
	fListHist2Rho4Pion->Add(fHistNeventsEtaCRhoChannel);

	//Create a histotram to count four pion events after the event cut for charge 
	fHistNeventsFourPicharge = new TH1D("Four Pi Events non zero charge", "Four Pi Events Non zero charge", 2, 0.5, 2.5);
	for (Int_t i = 0; i < 2; i++) fHistNeventsFourPicharge->GetXaxis()->SetBinLabel(i + 1, CutNameFourPiCharge[i].Data());
	fListHist2Rho4Pion->Add(fHistNeventsFourPicharge);
	//fListHist ->Add(fHistNeventsFourPicharge);

	TString FourPiEventCandidatesNames[4] = { "Four good tracks","4-#pi 0-M/B","3-#pi 1-M/B","2-#pi 2-M/B" };

	f4PiEventCandidates = new TH1D("Two PiPi Event Candidates", "Candidates for 2(#pi^{+} #pi^{-}) Channel Analysis", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) f4PiEventCandidates->GetXaxis()->SetBinLabel(i + 1, FourPiEventCandidatesNames[i].Data());
	fListHist2Rho4Pion->Add(f4PiEventCandidates);

	f2PairRhoM2VsPiPiM2 = new TH2D("f2PairRhoM2VsPiPiM2", "f2PairRhoM2VsPiPiM2", 500, 0., 0.5, 500, 0., 0.05);
	fListHist2Rho4Pion->Add(f2PairRhoM2VsPiPiM2);
	f1PairRhoM2VsPiPiM2 = new TH2D("f1PairRhoM2VsPiPiM2", "f1PairRhoM2VsPiPiM2", 500, 0., 0.5, 500, 0., 0.05);
	fListHist2Rho4Pion->Add(f1PairRhoM2VsPiPiM2);
	f2RhoM2VsPiPiM2 = new TH2D("f2RhoM2VsPiPiM2", "f2RhoM2VsPiPiM2", 500, 0., 0.5, 500, 0., 0.05);
	fListHist2Rho4Pion->Add(f2RhoM2VsPiPiM2);

	f2Rho2PairPtVsMinvEtaC = new TH2D("f2Rho2PairPtVsMinvEtaC", "f2Rho2PairPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList2Rho->Add(f2Rho2PairPtVsMinvEtaC);
	f2Rho2PairPtVsMinvFirstPair = new TH2D("f2Rho2PairPtVsMinvFirstPair", "f2Rho2PairPtVsMinvFirstPair", 600, 0., 6., 500, 0., 5.);
	fList2Rho->Add(f2Rho2PairPtVsMinvFirstPair);
	f2Rho2PairPtVsMinvSecondPair = new TH2D("f2Rho2PairPtVsMinvSecondPair", "f2Rho2PairPtVsMinvSecondPair", 600, 0., 6., 500, 0., 5.);
	fList2Rho->Add(f2Rho2PairPtVsMinvSecondPair);
	f2Rho2PairEtaVsMinvEtaC = new TH2D("f2Rho2PairEtaVsMinvEtaC", "f2Rho2PairEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList2RhoEtaC->Add(f2Rho2PairEtaVsMinvEtaC);
	f2Rho2PairEtaVsMinvEtaC400MeVPtMax = new TH2D("f2Rho2PairEtaVsMinvEtaC400MeVPtMax", "f2Rho2PairEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList2RhoEtaC->Add(f2Rho2PairEtaVsMinvEtaC400MeVPtMax);
	f2Rho2PairEtaVsMinvEtaC100MeVPtMax = new TH2D("f2Rho2PairEtaVsMinvEtaC100MeVPtMax", "f2Rho2PairEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList2RhoEtaC->Add(f2Rho2PairEtaVsMinvEtaC100MeVPtMax);
	f2Rho2PairSumPzVsMinvEtaC = new TH2D("f2Rho2PairSumPzVsMinvEtaC", "f2Rho2PairSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList2RhoEtaC->Add(f2Rho2PairSumPzVsMinvEtaC);
	f2Rho2PairScalarSumP = new TH1D("f2Rho2PairScalarSumP", "f2Rho2PairScalarSumP", 1000, 0., 10.);
	fList2RhoEtaC->Add(f2Rho2PairScalarSumP);
	f2Rho2PairVectorSumPt = new TH1D("f2Rho2PairVectorSumPt", "f2Rho2PairVectorSumPt", 1000, 0., 1.);
	fList2RhoEtaC->Add(f2Rho2PairVectorSumPt);

	f2Rho1PairPtVsMinvEtaC = new TH2D("f2Rho1PairPtVsMinvEtaC", "f2Rho1PairPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList1Rho->Add(f2Rho1PairPtVsMinvEtaC);
	f2Rho1PairPtVsMinvRhoPair = new TH2D("f2Rho1PairPtVsMinvRhoPair", "f2Rho1PairPtVsMinvRhoPair", 600, 0., 6., 500, 0., 5.);
	fList1Rho->Add(f2Rho1PairPtVsMinvRhoPair);
	f2Rho1PairPtVsMinvNonRhoPair = new TH2D("f2Rho1PairPtVsMinvNonRhoPair", "f2Rho1PairPtVsMinvNonRhoPair", 600, 0., 6., 500, 0., 5.);
	fList1Rho->Add(f2Rho1PairPtVsMinvNonRhoPair);
	f2Rho1PairEtaVsMinvEtaC = new TH2D("f2Rho1PairEtaVsMinvEtaC", "f2Rho1PairEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList1RhoEtaC->Add(f2Rho1PairEtaVsMinvEtaC);
	f2Rho1PairEtaVsMinvEtaC400MeVPtMax = new TH2D("f2Rho1PairEtaVsMinvEtaC400MeVPtMax", "f2Rho1PairEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList1RhoEtaC->Add(f2Rho1PairEtaVsMinvEtaC400MeVPtMax);
	f2Rho1PairEtaVsMinvEtaC100MeVPtMax = new TH2D("f2Rho1PairEtaVsMinvEtaC100MeVPtMax", "f2Rho1PairEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList1RhoEtaC->Add(f2Rho1PairEtaVsMinvEtaC100MeVPtMax);
	f2Rho1PairSumPzVsMinvEtaC = new TH2D("f2Rho1PairSumPzVsMinvEtaC", "f2Rho1PairSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList1RhoEtaC->Add(f2Rho1PairSumPzVsMinvEtaC);
	f2Rho1PairScalarSumP = new TH1D("f2Rho1PairScalarSumP", "f2Rho1PairScalarSumP", 1000, 0., 10.);
	fList1RhoEtaC->Add(f2Rho1PairScalarSumP);
	f2Rho1PairVectorSumPt = new TH1D("f2Rho1PairVectorSumPt", "f2Rho1PairVectorSumPt", 1000, 0., 1.);
	fList1RhoEtaC->Add(f2Rho1PairVectorSumPt);

	f4PionPtVsMinvEtaC = new TH2D("f4PionPtVsMinvEtaC", "f4PionPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList0Rho->Add(f4PionPtVsMinvEtaC);
	f4PionPtVsMinvRho = new TH2D("f4PionPtVsMinvRho", "f4PionPtVsMinvRho", 600, 0., 6., 500, 0., 5.);
	fList0Rho->Add(f4PionPtVsMinvRho);
	f4PionEtaVsMinvEtaC = new TH2D("f4PionEtaVsMinvEtaC", "f4PionEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList0Rho->Add(f4PionEtaVsMinvEtaC);
	f4PionEtaVsMinvEtaC400MeVPtMax = new TH2D("f4PionEtaVsMinvEtaC400MeVPtMax", "f4PionEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList0RhoEtaC->Add(f4PionEtaVsMinvEtaC400MeVPtMax);
	f4PionEtaVsMinvEtaC100MeVPtMax = new TH2D("f4PionEtaVsMinvEtaC100MeVPtMax", "f4PionEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList0RhoEtaC->Add(f4PionEtaVsMinvEtaC100MeVPtMax);
	f4PionSumPzVsMinvEtaC = new TH2D("f4PionSumPzVsMinvEtaC", "f4PionSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList0RhoEtaC->Add(f4PionSumPzVsMinvEtaC);
	f4PionScalarSumP = new TH1D("f4PionScalarSumP", "f4PionScalarSumP", 1000, 0., 10.);
	fList0RhoEtaC->Add(f4PionScalarSumP);
	f4PionVectorSumPt = new TH1D("f4PionVectorSumPt", "f4PionVectorSumPt", 1000, 0., 1.);
	fList0RhoEtaC->Add(f4PionVectorSumPt);


	//2Rho0 channel Helicity Cut histograms
	f2RhoParentPx = new TH1D("f2RhoParentPx", "Parent p_{x} in rest frame (2#rho^{0} channel", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPx);
	f2RhoParentPy = new TH1D("f2RhoParentPy", "Parent p_{y} in rest frame (2#rho^{0} channel", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPy);
	f2RhoParentPz = new TH1D("f2RhoParentPz", "Parent p_{z} in rest frame (2#rho^{0} channel", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPz);
	f2RhoDaughterParentAngle = new TH1D("f2RhoDaughterParentAngle", "Angle between daughter and parent", 3800, -190., 190.);
	fList2RhoCandidates->Add(f2RhoDaughterParentAngle);
	f2RhoDaughterParentCosAngle = new TH1D("f2RhoDaughterParentCosAngle", "Cosine of angle between daughter and parent", 2200, -1.1, 1.1);
	fList2RhoCandidates->Add(f2RhoDaughterParentCosAngle);
	f2RhoDaughterDaughterAngle = new TH1D("f2RhoDaughterDaughterAngle", "Angle between the two daughters", 2000, -10., 190.);
	fList2RhoCandidates->Add(f2RhoDaughterDaughterAngle);
	f2RhoDaughterDaughterCosAngle = new TH1D("f2RhoDaughterDaughterCosAngle", "Cos(Angle) between two daughters", 220, -1.1, 1.1);
	fList2RhoCandidates->Add(f2RhoDaughterDaughterCosAngle);
	f2RhoDaughterPtotal = new TH1D("f2RhoDaughterPtotal", "Momentum sum of two daughters", 1000, -5., 5.);
	fList2RhoCandidates->Add(f2RhoDaughterPtotal);

	//2Rho0 channel Helicity Cut histograms - Check histos
	f2RhoParentPxCheck = new TH1D("f2RhoParentPxCheck", "Parent p_{x} in rest frame (2#rho^{0} channel (cuts passed)", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPxCheck);
	f2RhoParentPyCheck = new TH1D("f2RhoParentPyCheck", "Parent p_{y} in rest frame (2#rho^{0} channel (cuts passed)", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPyCheck);
	f2RhoParentPzCheck = new TH1D("f2RhoParentPzCheck", "Parent p_{z} in rest frame (2#rho^{0} channel (cuts passed)", 200, -1., 1.);
	fList2RhoCandidates->Add(f2RhoParentPzCheck);
	f2RhoDaughterParentAngleCheck = new TH1D("f2RhoDaughterParentAngleCheck", "Angle between daughter and parent (cuts passed)", 3800, -190., 190.);
	fList2RhoCandidates->Add(f2RhoDaughterParentAngleCheck);
	f2RhoDaughterParentCosAngleCheck = new TH1D("f2RhoDaughterParentCosAngleCheck", "Cosine of angle between daughter and parent (cuts passed)", 2200, -1.1, 1.1);
	fList2RhoCandidates->Add(f2RhoDaughterParentCosAngleCheck);
	f2RhoDaughterDaughterAngleCheck = new TH1D("f2RhoDaughterDaughterAngleCheck", "Angle between the two daughters (cuts passed)", 2000, -10., 190.);
	fList2RhoCandidates->Add(f2RhoDaughterDaughterAngleCheck);
	f2RhoDaughterDaughterCosAngleCheck = new TH1D("f2RhoDaughterDaughterCosAngleCheck", "Cos(Angle) between two daughters (cuts passed)", 220, -1.1, 1.1);
	fList2RhoCandidates->Add(f2RhoDaughterDaughterCosAngleCheck);
	f2RhoDaughterPtotalCheck = new TH1D("f2RhoDaughterPtotalCheck", "Momentum sum of two daughters (cuts passed)", 1000, -5., 5.);
	fList2RhoCandidates->Add(f2RhoDaughterPtotalCheck);

	InitSystematics();

	//cout << "##### End of UserCreateOutputObjects()" << endl;

	PostData(1, fListTrig);
	PostData(2, fListHist);
	PostData(3, fListHistKstar);
	PostData(4, fListHist2Rho4Pion);
	PostData(5, fListHistK0s3PiPi4K);
}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::InitSystematics()
{
	//Initialize variables needed for RunAODsystematics (not currently in use)
	fListEtaCLoose = new TList();
	fListEtaCLoose->SetOwner();
	fListEtaCLoose->SetName("EtaCLoose");
	fListSystematics->Add(fListEtaCLoose);

	TH1D *fHistEtaCNClusLoose = new TH1D("EtaCNClusLoose", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCLoose->Add(fHistEtaCNClusLoose);

	TH1D *fHistEtaCChi2Loose = new TH1D("EtaCChi2Loose", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCLoose->Add(fHistEtaCChi2Loose);

	TH1D *fHistEtaCDCAzLoose = new TH1D("EtaCDCAzLoose", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCLoose->Add(fHistEtaCDCAzLoose);

	TH1D *fHistEtaCDCAxyLoose = new TH1D("EtaCDCAxyLoose", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCLoose->Add(fHistEtaCDCAxyLoose);

	TH1D *fHistEtaCITShitsLoose = new TH1D("EtaCITShitsLoose", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCLoose->Add(fHistEtaCITShitsLoose);


	fListEtaCTight = new TList();
	fListEtaCTight->SetOwner();
	fListEtaCTight->SetName("EtaCTight");
	fListSystematics->Add(fListEtaCTight);

	TH1D *fHistEtaCNClusTight = new TH1D("EtaCNClusTight", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCTight->Add(fHistEtaCNClusTight);

	TH1D *fHistEtaCChi2Tight = new TH1D("EtaCChi2Tight", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCTight->Add(fHistEtaCChi2Tight);

	TH1D *fHistEtaCDCAzTight = new TH1D("EtaCDCAzTight", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCTight->Add(fHistEtaCDCAzTight);

	TH1D *fHistEtaCDCAxyTight = new TH1D("EtaCDCAxyTight", "Invariant mass of #psi(2S) candidates", 50, 2.5, 5.5);
	fListEtaCTight->Add(fHistEtaCDCAxyTight);
}

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::UserExec(Option_t *) 
{

	//cout<<"#################### Beginning Next event... ##################"<<endl;
	//fType is 1 for AODs (in constructor)
	if (fType == 0) {
		RunESDtrig();
		if (fRunHist) RunESDhist();
		if (fRunTree) RunESDtree();
	}
	//fRunHist set to TRUE, fRunTree set to FALSE in constructor
	if (fType == 1) {
		RunAODtrig();
		if (fRunHist) RunAODhist();
		if (fRunTree) RunAODtree();
	}

}//UserExec

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunAODtrig()
{
  //  cout << "########## Beginning of RunAODtrig()" << endl;
  //##### This function fills the event, triggerClass, and run histograms. #####
  //get event from the input root file
  AliAODEvent *aod = (AliAODEvent*) InputEvent();
  if(!aod) return;

  fRunNum = aod ->GetRunNumber();
  //Trigger
  TString trigger = aod->GetFiredTriggerClasses();
  
  if(trigger.Contains("CCUP29")) fHistCcup29TriggersPerRun->Fill(fRunNum); //CCUP29 triggers
  if(trigger.Contains("CCUP30")) fHistCcup30TriggersPerRun->Fill(fRunNum); //CCUP7 triggers
  if(trigger.Contains("CCUP31")) fHistCcup31TriggersPerRun->Fill(fRunNum); //CCUP2 triggers
  
  if(trigger.Contains("CTRUE-B")) fHistCtrueTriggersPerRun->Fill(fRunNum); //CTRUE triggers
  
  fL1inputs = aod->GetHeader()->GetL1TriggerInputs();
  if(fL1inputs & (1 << 18)) fHistZedTriggersPerRun->Fill(fRunNum); //1ZED trigger inputs
  
  //MB, Central and SemiCentral triggers
  AliCentrality *centrality = aod->GetCentrality();
  UInt_t selectionMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  Double_t percentile = centrality->GetCentralityPercentileUnchecked("V0M");
  //Double_t percentile = centrality->GetCentralityPercentile("V0M");
  
  //if(((selectionMask & AliVEvent::kMB) == AliVEvent::kMB) && percentile<=80 && percentile>=0) fHistMBTriggersPerRun->Fill(fRunNum);
  
  //if(((selectionMask & AliVEvent::kCentral) == AliVEvent::kCentral) && percentile<=6 && percentile>=0 && (trigger.Contains("CCUP"))) fHistCentralTriggersPerRun->Fill(fRunNum);

  //if(((selectionMask & AliVEvent::kSemiCentral) == AliVEvent::kSemiCentral) && percentile<=50 && percentile>=15) fHistSemiCentralTriggersPerRun->Fill(fRunNum);

  //  cout << "########## End of RunAODtrig()" << endl;

PostData(1, fListTrig);

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunAODhist()
{
	//#####  Get particle data from the PDG database  #####
	TDatabasePDG *pdgdat = TDatabasePDG::Instance();

	TParticlePDG *partKaon = pdgdat->GetParticle(321);
	Double_t kaonMass = partKaon->Mass();

	TParticlePDG *partPion = pdgdat->GetParticle(211);
	Double_t pionMass = partPion->Mass();

	TParticlePDG *partKstar = pdgdat->GetParticle(313);
	Double_t kStarMass = partKstar->Mass();
	Double_t kStarWidth = partKstar->Width();

	TParticlePDG *partK0short = pdgdat->GetParticle(310);
	Double_t k0ShortMass = partK0short->Mass();
	Double_t k0ShortWidth = partK0short->Width();

	TParticlePDG *partRho = pdgdat->GetParticle(113);
	Double_t rhoMass = partRho->Mass();
	Double_t rhoWidth = partRho->Width();

	Double_t etaCMass = 2.98;
	Double_t etaCWidth = 0.032;

	//##### Get Event #####
	AliAODEvent *aod = (AliAODEvent*)InputEvent();
	if (!aod) return; //if return condition is met, next event is analyzed from the top of RunAODhist()
	//cout<<"Event number: "<<((TTree*) GetInputData(0))->GetTree()->GetReadEntry()<<endl;

	//place all AOD events in first bin
	fHistNeventsEtaC->Fill(1);

	//##### Apply Event Cuts. #####
	
	TString trigger = aod->GetFiredTriggerClasses();
	//Trigger condition -- must not be MC event, must contain and CCUP trigger to continue
	if (!isMC && !trigger.Contains("CCUP")) return;

	//place all triggered events in the second bin
	fHistNeventsEtaC->Fill(2);

	//vertex cut -- must have 2 or more tracks contributing to primary vertex
	AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
	fVtxContrib = fAODVertex->GetNContributors();

	if (fVtxContrib < 4) return;
	fHistNeventsEtaC->Fill(3);

	//v0 decision cut -- V0A and V0C must be empty to continue
	AliAODVZERO *fV0data = aod->GetVZEROData();
	fV0Adecision = fV0data->GetV0ADecision();
	fV0Cdecision = fV0data->GetV0CDecision();

	if (fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;
	fHistNeventsEtaC->Fill(4);

	//neutron ZDC cut -- must have <8 neutrons
	AliAODZDC *fZDCdata = aod->GetZDCData();
	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
	fZDCAtime = fZDCdata->GetZNATime();
	fZDCCtime = fZDCdata->GetZNCTime();

	fHistZDCAenergy->Fill(fZNAenergy);
	fHistZDCCenergy->Fill(fZNCenergy);
	fHistZDCAtime->Fill(fZDCAtime);

	//fHistZDCAtime->Fit("gaus");
	fHistZDCCtime->Fill(fZDCCtime);
	fHistZDCImpactParameter->Fill(fZDCdata->GetImpactParameter());
	fHistZDCAImpactParameter->Fill(fZDCdata->GetImpactParamSideA());
	fHistZDCCImpactParameter->Fill(fZDCdata->GetImpactParamSideC());

	if (trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1); //used to compair the number of neutrons per event, binned by number of neutrons
	if (fZNAenergy < 11500 && fZNCenergy < 11500) fHistZDCCuts->Fill(2); //was 8200, fills events with low ZNA and ZNC energy, <8 neutrons
	if (fZNAenergy < 1500 && fZNCenergy < 1500) fHistZDCCuts->Fill(3); //was 683, fills events with very low ZNA and ZNC energy, 0 neutrons
	if (fZDCAtime == 0 && fZDCCtime == 0) fHistZDCCuts->Fill(4);

	if (fZNAenergy > 11500 || fZNCenergy > 11500) return; // was >8200 in 2011 code, >1500 for 0n0n or <= 1500 for XnXn; one ZDC <= 1500 AND other ZDC > 1500 for 0nXn or Xn0n
	fHistNeventsEtaC->Fill(5);

	//##### Systematics - cut variation (Skipping this function right now.) #####
	if (fRunSystematics) RunAODsystematics(aod);

	//##### Begin Event Analysis #####
	Int_t nGoodTracks = 0;
	Int_t trackIndex[7] = { -1,-1,-1,-1,-1,-1,-1 };
	UInt_t nSpdHits = 0;

	//##### Track Cuts - Identify and count "good quality" tracks. (Same, standard cuts from Psi2S code.) #####
	for (Int_t itr = 0; itr < aod->GetNumberOfTracks(); itr++) { //check all tracks for quality
		AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
		if (!trk) continue; //continue restarts the for loop, starting with the next iteration
		if (!(trk->TestFilterBit(1 << 0))) continue;
		if (!(trk->GetStatus() & AliESDtrack::kTPCrefit)) continue;
		if (!(trk->GetStatus() & AliESDtrack::kITSrefit)) continue;
		if (trk->GetTPCNcls() < 50) continue; //requires at least 50 good TPC clusters
		if (trk->Chi2perNDF() > 4) continue; //requires (at least one?) chi^2 per degree of freedom systematvalue to be <4

		Double_t dca[2] = { 0.0,0.0 }, cov[3] = { 0.0,0.0,0.0 }; //distance of closest approach
		AliAODTrack* trk_clone = (AliAODTrack*)trk->Clone("trk_clone");
		if (!trk_clone->PropagateToDCA(fAODVertex, aod->GetMagneticField(), 300., dca, cov)) continue;
		delete trk_clone;
		if (TMath::Abs(dca[1]) > 2) continue;
		Double_t cut_DCAxy = 4 * (0.0182 + 0.0350 / TMath::Power(trk->Pt(), 1.01));
		if (TMath::Abs(dca[0]) > cut_DCAxy) continue;
		if ((trk->HasPointOnITSLayer(0)) || (trk->HasPointOnITSLayer(1))) nSpdHits++;
		trackIndex[nGoodTracks] = itr;
		nGoodTracks++;

		if (nGoodTracks > 10) break;
	}
	if (nSpdHits > 1) fHistNTracks->Fill(nGoodTracks);

	//##### Particle Identification #####
	TLorentzVector vPion[6], vKaon[6];
	Short_t qKaon[6], qPion[6];
	UInt_t nKaon = 0, nPion = 0, nBoth = 0, nNeither = 0;
	UInt_t nPreKaon = 0, nPrePion = 0; //added by Alec
	Double_t fRecTPCsignalPion[6], fRecTPCsignalKaon[6];
	Float_t eventStartTime = 0.;
	Double_t beta = 0;
	Int_t nTracksWithoutTOFinfo = 0;
	Double_t PCutLow = 0.6, PCutHigh = 3.;

	if (( nGoodTracks == 4) && nSpdHits > 1) {
		for (Int_t i = 0; i < nGoodTracks; i++) {
			AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
			if (!trk) AliFatal("Not a standard AOD");
			//Diagnostic figures regarding nSigma TOF vs TPC, etc.

			//Get nsigma info for PID
			fPIDTPCMuon[i] = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kMuon);
			fPIDTPCElectron[i] = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kElectron);
			fPIDTPCPion[i] = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kPion);
			fPIDTPCKaon[i] = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kKaon);
			fPIDTPCProton[i] = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton);

			if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trk) == AliPIDResponse::kDetPidOk) { //3 = kTOF
				fPIDTOFMuon[i] = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kMuon);
				fPIDTOFElectron[i] = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kElectron);
				fPIDTOFPion[i] = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kPion);
				fPIDTOFKaon[i] = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kKaon);
				fPIDTOFProton[i] = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton);
			}
			else {
				fPIDTOFMuon[i] = -999.;
				fPIDTOFElectron[i] = -999.;
				fPIDTOFPion[i] = -999.;
				fPIDTOFKaon[i] = -999.;
				fPIDTOFProton[i] = -999.;
			}

			//nSigmaTPC/TOF vs nSigmaTPC/TOF for all particles with Pion mass hypothesis and Kaon mass hypothesis. Histograms divided by pT range where the TPC and TOF signals do/don't overlap
			if (trk->P() < PCutLow) {
				fPionTPCvsPionTOFLowP->Fill(fPIDTOFPion[i], fPIDTPCPion[i]);
				fKaonTPCvsKaonTOFLowP->Fill(fPIDTOFKaon[i], fPIDTPCKaon[i]);
				fPionTPCvsKaonTPCLowP->Fill(fPIDTPCKaon[i], fPIDTPCPion[i]);
			}
			else if (trk->P() >= PCutLow && trk->P() < PCutHigh) {
				fPionTPCvsPionTOFMidP->Fill(fPIDTOFPion[i], fPIDTPCPion[i]);
				fKaonTPCvsKaonTOFMidP->Fill(fPIDTOFKaon[i], fPIDTPCKaon[i]);
				fPionTPCvsKaonTPCMidP->Fill(fPIDTPCKaon[i], fPIDTPCPion[i]);
				fPionTOFvsKaonTOFMidP->Fill(fPIDTOFKaon[i], fPIDTOFPion[i]);
				if (fabs(fPIDTPCPion[i]) < 3.) fPionTOFvsKaonTOFPionMidP->Fill(fPIDTOFKaon[i], fPIDTOFPion[i]);
				if (fabs(fPIDTPCKaon[i]) < 3.) fPionTOFvsKaonTOFKaonMidP->Fill(fPIDTOFKaon[i], fPIDTOFPion[i]);
			}
			else if (trk->P() >= PCutHigh) {
				fPionTPCvsPionTOFHighP->Fill(fPIDTOFPion[i], fPIDTPCPion[i]);
				fKaonTPCvsKaonTOFHighP->Fill(fPIDTOFKaon[i], fPIDTPCKaon[i]);
				fPionTPCvsKaonTPCHighP->Fill(fPIDTPCKaon[i], fPIDTPCPion[i]);
			}
			//TPC signal (dE/dx) or TOF signal vs TOF beta
			eventStartTime = aod->GetTOFHeader()->GetDefaultEventTimeVal(); //may have to define this differently. An array of start times is returned
			if (fPIDTOFPion[i] != -999.) fTOFIntegratedLength->Fill(trk->GetIntegratedLength()); //dagnostic histogram for intergrated track length
			if (trk->GetIntegratedLength() >= 360. && trk->GetIntegratedLength() <= 800. && trk->GetTOFsignal() > 0. && eventStartTime < 999990.0) {
				beta = (trk->GetIntegratedLength()*0.01) / ((trk->GetTOFsignal() - eventStartTime)*(1.e-12)*TMath::C()); //0.01 converts cm to m, 1e-12 converts ps to sec
				//NOTE: sometimes one can calculate beta without being able to calculate TOF nSigma!

				fTPCdEdxVsTOFbetaAll->Fill(beta, trk->GetTPCsignal());
				fTOFbetaVsPAll->Fill(trk->P(), beta);
				if (((fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && ((fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns pion and TOF returns pion, both, or missing
					|| (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns pion
					&& (beta > 0.55 || beta == -999.)) {
					fTPCdEdxVsTOFbetaPionsWithPID->Fill(beta, trk->GetTPCsignal());
					fTOFbetaVsPPion->Fill(trk->P(), beta);
				}
				else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && ((fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns kaon and TOF returns kaon, both, or missing
					|| (fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns kaon
					&& (beta > 0.55 || beta == -999.)) {
					fTPCdEdxVsTOFbetaKaonsWithPID->Fill(beta, trk->GetTPCsignal());
					fTOFbetaVsPKaon->Fill(trk->P(), beta);
				}
				else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)  //if TPC returns kaon and TOF returns pion
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)  //if TPC returns pion and TOF returns kaon
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.))) //if TPC and TOF both return both or missing
					&& (beta > 0.55 || beta == -999.)) {
					fTOFbetaVsPBoth->Fill(trk->P(), beta);
				}
			}
			fDedxVsPAll->Fill(trk->P(), trk->GetTPCsignal());
			
			//##### Particle Identification #####
			//Identify pions
			if (((fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && ((fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns pion and TOF returns pion, both, or missing
				|| (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns pion
				&& (beta > 0.55 || beta == -999.)) { //apparently being able to calculate a beta value is not perfectly correlated with having a TOF signal (only missing TOF info cases should get a beta value with no TOF value though)
				fDedxVsPPion->Fill(trk->P(), trk->GetTPCsignal());
				fRecTPCsignalPion[nPion] = trk->GetTPCsignal(); //set the first slot, corresponding to the nPion number, in each array to have the info from this track
				qPion[nPion] = trk->Charge();
				vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
				nPion++; //move on to pion #2
			}
			//Identify tracks which are likely to be a kaon
			else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && ((fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns kaon and TOF returns kaon, both, or missing
				|| (fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns kaon
				&& (beta > 0.55 || beta == -999.)) {
				fDedxVsPKaon->Fill(trk->P(), trk->GetTPCsignal());
				fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
				qKaon[nKaon] = trk->Charge();
				vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
				nKaon++;
			}
			//Identify tracks which could be either a kaon or pion
			else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)  //if TPC returns kaon and TOF returns pion
				|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)  //if TPC returns pion and TOF returns kaon
				|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.))) //if TPC and TOF both return both or missing
				&& (beta > 0.55 || beta == -999.)) {
				fDedxVsPBoth->Fill(trk->P(), trk->GetTPCsignal());
				nBoth++;
			}
			//Identify tracks which are likely neither a kaon nor pion
			else nNeither++;

			if (fPIDTOFPion[i] == -999.) nTracksWithoutTOFinfo++; //count number ofr tracks with missing TOF signal
		}
		if (nGoodTracks == 4) {
			fHistFourTracksNpion->Fill(nPion);
			fHistFourTracksNkaon->Fill(nKaon);
			fFourTrackMissing->Fill(nTracksWithoutTOFinfo);
			fHistFourTracksNboth->Fill(nBoth);
			fHistFourTracksNneither->Fill(nNeither);

			//diagnostic histograms -- events in these bins should match with number of events analyzed
			f4PiEventCandidates->Fill(1);
			if (nPion == 4) f4PiEventCandidates->Fill(2);
			if (nPion == 3 && nBoth == 1) f4PiEventCandidates->Fill(3);
			if (nPion == 2 && nBoth == 2) f4PiEventCandidates->Fill(4);

		}
		//PID check counts total identified tracks -- should always be four
		fHistPIDCheck->Fill(nPion + nKaon + nBoth + nNeither);
		nPrePion = nPion;
	}
	//if (nSpdHits > 1) cout << "PID block Completed. Good Tracks: " << nGoodTracks << endl;
	if (nGoodTracks != 4 || nSpdHits < 2) return;
	//if ((qpos - qneg) != 0) return;
	if (nGoodTracks == 4 || nSpdHits > 1) fHistNeventsEtaC->Fill(6);
	
	//##### Analysis by Decay Channel #####
	Int_t nHighPtTracks = 0;
	Double_t SumPz = 0, VectorSumPt = 0, ScalarSumP = 0;
	TLorentzVector vCandidate;
	TVector3 sumPtVector;
	UInt_t newPion = 0;  //Added by Alec

	//EtaC->RhoRho, 2(Pi+ Pi-) Channels
	newPion = 0;
	nHighPtTracks = 0;
	Int_t nPiMinus = 0;
	Int_t nPiPlus = 0;
	Int_t nRhoPairs = 0;
	Int_t caseOne = 0;
	Int_t caseTwo = 0;
	Bool_t goodRho[4] = { kFALSE, kFALSE, kFALSE, kFALSE };
	TLorentzVector vPionMinus[4], vPionPlus[4], vRho[4];
	Double_t boostInfoRhoZero[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };
	Double_t boostInfoRhoOne[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };
	Double_t boostInfoRhoTwo[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };
	Double_t boostInfoRhoThree[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };

	if (nGoodTracks == 4 && nSpdHits > 1) {
		//cout << "##### 4Pi channel Started." << endl;

		//##### Deal with missing tracks due to TOF
		if ((nPion == 2 && nBoth == 2) || (nPion == 3 && nBoth == 1)) {
			for (int i = 0; i < 4; i++) { //loop over the four tracks
				AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
				if ((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)  //if both
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.)
					|| (fPIDTOFPion[i] == -999. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.)) { //or if missing
					fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
					qPion[nPion] = trk->Charge();
					vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass); 
					nPion++;
					newPion++;
				}
			}
		}

		if (nKaon+nNeither == 0) f4PinPionVsnewPion->Fill(newPion, nPion-newPion); //events with only pions and missings

		//Analyze good events, fill the histos.
		if ((nPion == 4)) {
			//cout << "##### Four Pi channel analyzing..." << endl;
			fHistNeventsEtaCRhoChannel->Fill(1);
			fHistNeventsFourPicharge->Fill(1);
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3]) != 0) {
				fHistNeventsEtaCRhoChannel->Fill(2); //non-zero net charge
				fHistNeventsFourPicharge->Fill(2);
			}
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3]) == 0) {
				fHistNeventsEtaCRhoChannel->Fill(3); //zero net charge
				
				vCandidate = vPion[0] + vPion[1] + vPion[2] + vPion[3];
				//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

				SumPz = vPion[0].Pz() + vPion[1].Pz() + vPion[2].Pz() + vPion[3].Pz();
				ScalarSumP = vPion[0].P() + vPion[1].P() + vPion[2].P() + vPion[3].P();
				sumPtVector = vPion[0].Vect() + vPion[1].Vect() + vPion[2].Vect() + vPion[3].Vect();
				VectorSumPt = sumPtVector.Pt();

				//Number of tracks with pT>0.4 GeV/c
				for (Int_t aa = 0; aa < 4; aa++) if (vPion[aa].Pt() > 0.1) nHighPtTracks++;

				//Get masses of potential intermediate Rho's
				if (qPion[0] < 0) {
					vPionMinus[nPiMinus] = vPion[0];
					nPiMinus++;
				}
				else {
					vPionPlus[nPiPlus] = vPion[0];
					nPiPlus++;
				}
				if (qPion[1] < 0) {
					vPionMinus[nPiMinus] = vPion[1];
					nPiMinus++;
				}
				else {
					vPionPlus[nPiPlus] = vPion[1];
					nPiPlus++;
				}
				if (qPion[2] < 0) {
					vPionMinus[nPiMinus] = vPion[2];
					nPiMinus++;
				}
				else {
					vPionPlus[nPiPlus] = vPion[2];
					nPiPlus++;
				}
				if (qPion[3] < 0) {
					vPionMinus[nPiMinus] = vPion[3];
					nPiMinus++;
				}
				else {
					vPionPlus[nPiPlus] = vPion[3];
					nPiPlus++;
				}
				//Either 0 and 1 are rho's or 2 and 3. If both sets are rho's choose best set.
				vRho[0] = vPionMinus[0] + vPionPlus[0];
				vRho[1] = vPionMinus[1] + vPionPlus[1];
				vRho[2] = vPionMinus[1] + vPionPlus[0];
				vRho[3] = vPionMinus[0] + vPionPlus[1];

				//Extract boostInfo for the two potential rho0's in the two sets.
				BoostCut(vPionMinus[0], vPionPlus[0], vRho[0], boostInfoRhoZero);
				BoostCut(vPionMinus[1], vPionPlus[1], vRho[1], boostInfoRhoOne);
				BoostCut(vPionMinus[1], vPionPlus[0], vRho[2], boostInfoRhoTwo);
				BoostCut(vPionMinus[0], vPionPlus[1], vRho[3], boostInfoRhoThree);

				//Fill histos with first set of rho0's
				f2RhoParentPx->Fill(boostInfoRhoZero[9]);                 f2RhoParentPx->Fill(boostInfoRhoOne[9]);
				f2RhoParentPy->Fill(boostInfoRhoZero[10]);                f2RhoParentPy->Fill(boostInfoRhoOne[10]);
				f2RhoParentPz->Fill(boostInfoRhoZero[11]);                f2RhoParentPz->Fill(boostInfoRhoOne[11]);
				f2RhoDaughterParentAngle->Fill(boostInfoRhoZero[4]);      f2RhoDaughterParentAngle->Fill(boostInfoRhoOne[4]);
				f2RhoDaughterParentAngle->Fill(boostInfoRhoZero[5]);      f2RhoDaughterParentAngle->Fill(boostInfoRhoOne[5]);
				f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoZero[4] * TMath::Pi() / 180.));      f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoOne[4] * TMath::Pi() / 180.));
				f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoZero[5] * TMath::Pi() / 180.));      f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoOne[5] * TMath::Pi() / 180.));
				f2RhoDaughterDaughterAngle->Fill(boostInfoRhoZero[6]);    f2RhoDaughterDaughterAngle->Fill(boostInfoRhoOne[6]);
				f2RhoDaughterDaughterCosAngle->Fill(boostInfoRhoZero[7]); f2RhoDaughterDaughterCosAngle->Fill(boostInfoRhoOne[7]);
				f2RhoDaughterPtotal->Fill(boostInfoRhoZero[8]);           f2RhoDaughterPtotal->Fill(boostInfoRhoOne[8]);

				//Fill histos with second set of rho0's
				f2RhoParentPx->Fill(boostInfoRhoTwo[9]);                 f2RhoParentPx->Fill(boostInfoRhoThree[9]);
				f2RhoParentPy->Fill(boostInfoRhoTwo[10]);                f2RhoParentPy->Fill(boostInfoRhoThree[10]);
				f2RhoParentPz->Fill(boostInfoRhoTwo[11]);                f2RhoParentPz->Fill(boostInfoRhoThree[11]);
				f2RhoDaughterParentAngle->Fill(boostInfoRhoTwo[4]);      f2RhoDaughterParentAngle->Fill(boostInfoRhoThree[4]);
				f2RhoDaughterParentAngle->Fill(boostInfoRhoTwo[5]);      f2RhoDaughterParentAngle->Fill(boostInfoRhoThree[5]);
				f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoTwo[4] * TMath::Pi() / 180.));      f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoThree[4] * TMath::Pi() / 180.));
				f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoTwo[5] * TMath::Pi() / 180.));      f2RhoDaughterParentCosAngle->Fill(TMath::Cos(boostInfoRhoThree[5] * TMath::Pi() / 180.));
				f2RhoDaughterDaughterAngle->Fill(boostInfoRhoTwo[6]);    f2RhoDaughterDaughterAngle->Fill(boostInfoRhoThree[6]);
				f2RhoDaughterDaughterCosAngle->Fill(boostInfoRhoTwo[7]); f2RhoDaughterDaughterCosAngle->Fill(boostInfoRhoThree[7]);
				f2RhoDaughterPtotal->Fill(boostInfoRhoTwo[8]);           f2RhoDaughterPtotal->Fill(boostInfoRhoThree[8]);

				//Apply helicity cut
				if (fabs(boostInfoRhoZero[4]) > 45. && fabs(boostInfoRhoZero[5]) > 45. && fabs(boostInfoRhoZero[6]) > 160. && fabs(boostInfoRhoZero[8]) < 0.2) goodRho[0] = kTRUE;
				else goodRho[0] = kFALSE;
				if (fabs(boostInfoRhoOne[4]) > 45. && fabs(boostInfoRhoOne[5]) > 45. && fabs(boostInfoRhoOne[6]) > 160. && fabs(boostInfoRhoOne[8]) < 0.2) goodRho[1] = kTRUE;
				else goodRho[1] = kFALSE;
				if (fabs(boostInfoRhoTwo[4]) > 45. && fabs(boostInfoRhoTwo[5]) > 45. && fabs(boostInfoRhoTwo[6]) > 160. && fabs(boostInfoRhoTwo[8]) < 0.2) goodRho[2] = kTRUE;
				else goodRho[2] = kFALSE;
				if (fabs(boostInfoRhoThree[4]) > 45. && fabs(boostInfoRhoThree[5]) > 45. && fabs(boostInfoRhoThree[6]) > 160. && fabs(boostInfoRhoThree[8]) < 0.2) goodRho[3] = kTRUE;
				else goodRho[3] = kFALSE;

				//Fill check histos for good rho0's
				if (goodRho[0]) {
					f2RhoParentPxCheck->Fill(boostInfoRhoZero[9]);
					f2RhoParentPyCheck->Fill(boostInfoRhoZero[10]);
					f2RhoParentPzCheck->Fill(boostInfoRhoZero[11]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoZero[4]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoZero[5]);
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoZero[4] * TMath::Pi() / 180.));
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoZero[5] * TMath::Pi() / 180.));
					f2RhoDaughterDaughterAngleCheck->Fill(boostInfoRhoZero[6]);
					f2RhoDaughterDaughterCosAngleCheck->Fill(boostInfoRhoZero[7]);
					f2RhoDaughterPtotalCheck->Fill(boostInfoRhoZero[8]);
				}
				if (goodRho[1]) {
					f2RhoParentPxCheck->Fill(boostInfoRhoOne[9]);
					f2RhoParentPyCheck->Fill(boostInfoRhoOne[10]);
					f2RhoParentPzCheck->Fill(boostInfoRhoOne[11]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoOne[4]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoOne[5]);
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoOne[4] * TMath::Pi() / 180.));
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoOne[5] * TMath::Pi() / 180.));
					f2RhoDaughterDaughterAngleCheck->Fill(boostInfoRhoOne[6]);
					f2RhoDaughterDaughterCosAngleCheck->Fill(boostInfoRhoOne[7]);
					f2RhoDaughterPtotalCheck->Fill(boostInfoRhoOne[8]);
				}
				if (goodRho[2]) {
					f2RhoParentPxCheck->Fill(boostInfoRhoTwo[9]);
					f2RhoParentPyCheck->Fill(boostInfoRhoTwo[10]);
					f2RhoParentPzCheck->Fill(boostInfoRhoTwo[11]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoTwo[4]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoTwo[5]);
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoTwo[4] * TMath::Pi() / 180.));
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoTwo[5] * TMath::Pi() / 180.));
					f2RhoDaughterDaughterAngleCheck->Fill(boostInfoRhoTwo[6]);
					f2RhoDaughterDaughterCosAngleCheck->Fill(boostInfoRhoTwo[7]);
					f2RhoDaughterPtotalCheck->Fill(boostInfoRhoTwo[8]);
				}
				if (goodRho[3]) {
					f2RhoParentPxCheck->Fill(boostInfoRhoThree[9]);
					f2RhoParentPyCheck->Fill(boostInfoRhoThree[10]);
					f2RhoParentPzCheck->Fill(boostInfoRhoThree[11]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoThree[4]);
					f2RhoDaughterParentAngleCheck->Fill(boostInfoRhoThree[5]);
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoThree[4] * TMath::Pi() / 180.));
					f2RhoDaughterParentCosAngleCheck->Fill(TMath::Cos(boostInfoRhoThree[5] * TMath::Pi() / 180.));
					f2RhoDaughterDaughterAngleCheck->Fill(boostInfoRhoThree[6]);
					f2RhoDaughterDaughterCosAngleCheck->Fill(boostInfoRhoThree[7]);
					f2RhoDaughterPtotalCheck->Fill(boostInfoRhoThree[8]);
				}

				//##### Turn off Helicity Cut #####
				goodRho[0] = kTRUE;    goodRho[1] = kTRUE;    goodRho[2] = kTRUE;    goodRho[3] = kTRUE;

				//Identify sets with 2 rho0's
				if (vRho[0].M() < (rhoMass + rhoWidth) && vRho[0].M() > (rhoMass - rhoWidth) &&
					vRho[1].M() < (rhoMass + rhoWidth) && vRho[1].M() > (rhoMass - rhoWidth) &&
					goodRho[0] && goodRho[1]) {
					nRhoPairs++; caseOne = 1;
				}
				if (vRho[2].M() < (rhoMass + rhoWidth) && vRho[2].M() > (rhoMass - rhoWidth) &&
					vRho[3].M() < (rhoMass + rhoWidth) && vRho[3].M() > (rhoMass - rhoWidth) &&
					goodRho[2] && goodRho[3]) {
					nRhoPairs++; caseTwo = 1;
				}

				//Dalitz plots -- look at rho-pi+ vs pi--pi+ minv
				f2RhoM2VsPiPiM2->Fill((vRho[0].M()*vPionPlus[1].M()), (vPionMinus[1].M()*vPionPlus[1].M()));
				f2RhoM2VsPiPiM2->Fill((vRho[1].M()*vPionPlus[0].M()), (vPionMinus[0].M()*vPionPlus[0].M()));
				f2RhoM2VsPiPiM2->Fill((vRho[2].M()*vPionPlus[1].M()), (vPionMinus[0].M()*vPionPlus[1].M()));
				f2RhoM2VsPiPiM2->Fill((vRho[3].M()*vPionPlus[0].M()), (vPionMinus[1].M()*vPionPlus[0].M()));

				if (nRhoPairs == 2) {
					fHistNeventsEtaCRhoChannel->Fill(4);
					fHistNeventsEtaC->Fill(7);
					fEtaCCandidatesPerChannel->Fill(4);
					if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(4);
					if (sqrt(pow(fabs(vRho[0].M() - rhoMass), 2.) + pow(fabs(vRho[1].M() - rhoMass), 2.)) < sqrt(pow(fabs(vRho[2].M() - rhoMass), 2.) + pow(fabs(vRho[3].M() - rhoMass), 2.))) {
						f2Rho2PairPtVsMinvFirstPair->Fill(vRho[0].M(), vRho[0].Pt());
						f2Rho2PairPtVsMinvFirstPair->Fill(vRho[1].M(), vRho[1].Pt());
						f2Rho2PairPtVsMinvSecondPair->Fill(vRho[2].M(), vRho[2].Pt());
						f2Rho2PairPtVsMinvSecondPair->Fill(vRho[3].M(), vRho[3].Pt());

						//Dalitz plots -- look at rho-pi+ vs pi--pi+ minv
						f2PairRhoM2VsPiPiM2->Fill((vRho[0].M()*vPionPlus[1].M()), (vPionMinus[1].M()*vPionPlus[1].M()));
						f2PairRhoM2VsPiPiM2->Fill((vRho[1].M()*vPionPlus[0].M()), (vPionMinus[0].M()*vPionPlus[0].M()));
					}
					else {
						f2Rho2PairPtVsMinvFirstPair->Fill(vRho[2].M(), vRho[2].Pt());
						f2Rho2PairPtVsMinvFirstPair->Fill(vRho[3].M(), vRho[3].Pt());

						f2Rho2PairPtVsMinvSecondPair->Fill(vRho[0].M(), vRho[0].Pt());
						f2Rho2PairPtVsMinvSecondPair->Fill(vRho[1].M(), vRho[1].Pt());

						f2PairRhoM2VsPiPiM2->Fill((vRho[2].M()*vPionPlus[1].M()), (vPionMinus[0].M()*vPionPlus[1].M()));
						f2PairRhoM2VsPiPiM2->Fill((vRho[3].M()*vPionPlus[0].M()), (vPionMinus[1].M()*vPionPlus[0].M()));
					}
					f2Rho2PairPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt()); //4Pi final states with 2 intermediate rho's.
					fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					if (vCandidate.Pt() < 0.11) {
						fChannelVsMinvEtaC->Fill(vCandidate.M(), 4);
						fAllMinvEtaCLowPt->Fill(vCandidate.M());
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(8);
					}
					f2Rho2PairEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.4) f2Rho2PairEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.1) f2Rho2PairEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					f2Rho2PairSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
					f2Rho2PairScalarSumP->Fill(ScalarSumP);
					f2Rho2PairVectorSumPt->Fill(VectorSumPt);
					if (nHighPtTracks > 1) fHistNeventsEtaCRhoChannel->Fill(7);
				}

				else if (nRhoPairs == 1) {
					fHistNeventsEtaCRhoChannel->Fill(5);
					fHistNeventsEtaC->Fill(7);
					fEtaCCandidatesPerChannel->Fill(4);
					if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(4);
					if (caseOne == 1) {
						f2Rho1PairPtVsMinvRhoPair->Fill(vRho[0].M(), vRho[0].Pt());
						f2Rho1PairPtVsMinvRhoPair->Fill(vRho[1].M(), vRho[1].Pt());
						f2Rho1PairPtVsMinvNonRhoPair->Fill(vRho[2].M(), vRho[2].Pt());
						f2Rho1PairPtVsMinvNonRhoPair->Fill(vRho[3].M(), vRho[3].Pt());

						//Dalitz plots -- look at rho-pi+ vs pi--pi+ minv
						f1PairRhoM2VsPiPiM2->Fill((vRho[0].M()*vPionPlus[1].M()), (vPionMinus[1].M()*vPionPlus[1].M()));
						f1PairRhoM2VsPiPiM2->Fill((vRho[1].M()*vPionPlus[0].M()), (vPionMinus[0].M()*vPionPlus[0].M()));
					}
					else if (caseTwo == 1) {
						f2Rho1PairPtVsMinvRhoPair->Fill(vRho[2].M(), vRho[2].Pt());
						f2Rho1PairPtVsMinvRhoPair->Fill(vRho[3].M(), vRho[3].Pt());
						f2Rho1PairPtVsMinvNonRhoPair->Fill(vRho[0].M(), vRho[0].Pt());
						f2Rho1PairPtVsMinvNonRhoPair->Fill(vRho[1].M(), vRho[1].Pt());

						f1PairRhoM2VsPiPiM2->Fill((vRho[2].M()*vPionPlus[1].M()), (vPionMinus[0].M()*vPionPlus[1].M()));
						f1PairRhoM2VsPiPiM2->Fill((vRho[3].M()*vPionPlus[0].M()), (vPionMinus[1].M()*vPionPlus[0].M()));
					}
					f2Rho1PairPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt()); //4Pi final states with 2 intermediate rho's.
					fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					if (vCandidate.Pt() < 0.11) {
						fChannelVsMinvEtaC->Fill(vCandidate.M(), 4);
						fAllMinvEtaCLowPt->Fill(vCandidate.M());
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(8);
					}
					f2Rho1PairEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.4) f2Rho1PairEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.1) f2Rho1PairEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					f2Rho1PairSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
					f2Rho1PairScalarSumP->Fill(ScalarSumP);
					f2Rho1PairVectorSumPt->Fill(VectorSumPt);
					if (nHighPtTracks > 1) fHistNeventsEtaCRhoChannel->Fill(7);
				}

				else if (nRhoPairs == 0) {
					fHistNeventsEtaCRhoChannel->Fill(6);
					fHistNeventsEtaC->Fill(7);
					fEtaCCandidatesPerChannel->Fill(5);
					if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(5);
					f4PionPtVsMinvRho->Fill(vRho[0].M(), vRho[0].Pt());
					f4PionPtVsMinvRho->Fill(vRho[1].M(), vRho[1].Pt());
					f4PionPtVsMinvRho->Fill(vRho[2].M(), vRho[2].Pt());
					f4PionPtVsMinvRho->Fill(vRho[3].M(), vRho[3].Pt());

					f4PionPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt()); //All 4Pi final states.
					fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					if (vCandidate.Pt() < 0.11) {
						fChannelVsMinvEtaC->Fill(vCandidate.M(), 5);
						fAllMinvEtaCLowPt->Fill(vCandidate.M());
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(8);
					}
					f4PionEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.4) f4PionEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.1) f4PionEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					f4PionSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
					f4PionScalarSumP->Fill(ScalarSumP);
					f4PionVectorSumPt->Fill(VectorSumPt);
					if (nHighPtTracks > 1) fHistNeventsEtaCRhoChannel->Fill(7);
				}


			}
		}
		//cout << "##### 4Pi channel Complete." << endl;
		//Remove the partiles assumed to match because of missing TOF for the next analysis
		nPion -= newPion;
	}
	//End EtaC->RhoRho Channel

	if (nGoodTracks == 4 && nSpdHits > 1) {
		fHistPostFourTracksNpion->Fill(nPion);
		fListNPionChange->Fill(nPion - nPrePion);
	}
	//  cout << "##### End of RunAODHist()" << endl;
	
	PostData(2, fListHist);
	PostData(3, fListHistKstar);
	PostData(4, fListHist2Rho4Pion);
	PostData(5, fListHistK0s3PiPi4K);
}

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunAODtree()
{

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunAODMC(AliAODEvent *aod)
{

}//RunAODMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunESDtrig()
{

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunESDhist()
{

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunESDtree()
{

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::RunESDMC(AliESDEvent* esd)
{

}//RunESDMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcFourPi::Terminate(Option_t *) 
{

  //cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcFourPi::GetMedian(Double_t *daArray) {
    // Allocate an array of the same size and sort it.
    Double_t dpSorted[4];
    for (Int_t i = 0; i < 4; ++i) {
        dpSorted[i] = daArray[i];
    }
    for (Int_t i = 3; i > 0; --i) {
        for (Int_t j = 0; j < i; ++j) {
            if (dpSorted[j] > dpSorted[j+1]) {
                Double_t dTemp = dpSorted[j];
                dpSorted[j] = dpSorted[j+1];
                dpSorted[j+1] = dTemp;
            }
        }
    }

    // Middle or average of middle values in the sorted array.
    Double_t dMedian = 0.0;
    dMedian = (dpSorted[2] + dpSorted[1])/2.0;
    
    return dMedian;
}

bool AliAnalysisTaskUpcFourPi::CheckMeritCutWinner(int cutChoice, double oldPars[3], double newPars[3]){
	//performs "merit cut" judgement check on v0s with shared daughters, using one of three criteria.
	//if cutChoice = 4 (it does by default), it uses all three criteria, needed 2 of 3 'points'

	//cout << "#################### Beginning of CheckMeritCutWinner()" << endl;

	bool newV0Wins = kFALSE;
	double pardiff[3] = { newPars[0] - oldPars[0],
						 newPars[1] - oldPars[1],
						 newPars[2] - oldPars[2] };
	if (cutChoice > 0 && cutChoice < 4) {
		if (pardiff[cutChoice] <= 0.) newV0Wins = kTRUE;
	}
	else if (cutChoice == 4) {
		int newWinCount = 0;
		for (int i = 0; i < 3; i++) { if (pardiff[i] <= 0) newWinCount++; } //I think that original pardiff[i+1] was an error.
		if (newWinCount > 1) newV0Wins = kTRUE;
	}
	else {};

	//cout << "#################### End of CheckMeritCutWinner()" << endl;

	return newV0Wins;
}

void AliAnalysisTaskUpcFourPi::BoostCut(TLorentzVector d1, TLorentzVector d2, TLorentzVector parent, Double_t *boostInfo) {

   //boost to parent rest frame
   TVector3 boostVector( -1.*parent.Px()/parent.E(), -1.*parent.Py()/parent.E(), -1.*parent.Pz()/parent.E());
   TLorentzVector d1_newSys(d1.Px(),d1.Py(),d1.Pz(),d1.E());
   d1_newSys.Boost(boostVector);
   TLorentzVector d2_newSys(d2.Px(),d2.Py(),d2.Pz(),d2.E());
   d2_newSys.Boost(boostVector);
   TLorentzVector parent_newSys(parent.Px(),parent.Py(),parent.Pz(),parent.E());
   parent_newSys.Boost(boostVector);

   //Get daughter components perpendicular to parent momentum
   Double_t d1perpMom = d1_newSys.Perp(parent.Vect());
   Double_t d2perpMom = d2_newSys.Perp(parent.Vect());

   //Get daughter components parallel to parent momentum
   Double_t d1parMom = sqrt(d1_newSys.Vect().Mag2() - d1_newSys.Perp2(parent.Vect()));
   Double_t d2parMom = sqrt(d2_newSys.Vect().Mag2() - d2_newSys.Perp2(parent.Vect()));

   //Get angle relative to parent momentum (in degrees)
   Double_t d1perpAngle = 180.*atan2(d1perpMom, d1parMom)/TMath::Pi();
   Double_t d2perpAngle = 180.*atan2(d2perpMom, d2parMom)/TMath::Pi();

   //Get angle between two daughters in the parent rest frame
   Double_t d1d2angle = d1_newSys.Angle(d2_newSys.Vect());
   Double_t d1d2cosAngle = cos(d1d2angle);
   d1d2angle = 180.*d1d2angle/TMath::Pi(); //convert to degrees

   //Get relative momentum of two daughters in parent rest frame
   Double_t d1d2Ptotal = sqrt(pow((d1_newSys.Px() + d2_newSys.Px()),2.) + pow((d1_newSys.Py() + d2_newSys.Py()),2.) + pow((d1_newSys.Pz() + d2_newSys.Pz()),2.));
   Double_t d1Mag = sqrt( pow(d1_newSys.Px(),2.) + pow(d1_newSys.Py(),2.) + pow(d1_newSys.Pz(),2.) );
   Double_t d2Mag = sqrt( pow(d2_newSys.Px(),2.) + pow(d2_newSys.Py(),2.) + pow(d2_newSys.Pz(),2.) );
   Double_t d1d2PtotalNorm = sqrt( pow((d1_newSys.Px()/d1Mag + d2_newSys.Px()/d2Mag),2.) + pow((d1_newSys.Py()/d1Mag + d2_newSys.Py()/d2Mag),2.) + pow((d1_newSys.Pz()/d1Mag + d2_newSys.Pz()/d2Mag),2.) );

   //   boostInfo[0] = { d1perpMom, d2perpMom, d1parMom, d2parMom, d1perpAngle, d2perpAngle, d1d2angle, d1d2cosAngle, d1d2RelativeMom, parent_newSys.Px(), parent_newSys.Py(), parent_newSys.Pz() };

   boostInfo[0] = d1perpMom;
   boostInfo[1] = d2perpMom;
   boostInfo[2] = d1parMom;
   boostInfo[3] = d2parMom;
   boostInfo[4] = d1perpAngle;
   boostInfo[5] = d2perpAngle;
   boostInfo[6] = d1d2angle;
   boostInfo[7] = d1d2cosAngle;
   boostInfo[8] = d1d2Ptotal;
   boostInfo[9] = parent_newSys.Px();
   boostInfo[10] = parent_newSys.Py();
   boostInfo[11] = parent_newSys.Pz();
   boostInfo[12] = d1d2PtotalNorm;

   return; // boostInfo;
 }

void AliAnalysisTaskUpcFourPi::RunAODsystematics(AliAODEvent* aod)
{

}
