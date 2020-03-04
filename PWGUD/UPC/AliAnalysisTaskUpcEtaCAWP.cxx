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
//full code, Alec's Local copy
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
#include "AliAnalysisTaskUpcEtaCAWP.h"

ClassImp(AliAnalysisTaskUpcEtaCAWP);

using std::cout;
using std::endl;

//trees for UPC EtaC analysis,fList1KStar
// christopher.anson@cern.ch

//_____________________________________________________________________________
AliAnalysisTaskUpcEtaCAWP::AliAnalysisTaskUpcEtaCAWP() //first constructor, contains all histograms
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
	fHistNeventsEtaCRhoChannel(0), f2PairRhoM2VsPiPiM2(0), f1PairRhoM2VsPiPiM2(0), f2RhoM2VsPiPiM2(0),
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

}//AliAnalysisTaskUpcEtaCAWP

//_____________________________________________________________________________
AliAnalysisTaskUpcEtaCAWP::AliAnalysisTaskUpcEtaCAWP(const char *name) //second constructor, contains all histograms
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
	fHistNeventsEtaCRhoChannel(0), f2PairRhoM2VsPiPiM2(0), f1PairRhoM2VsPiPiM2(0), f2RhoM2VsPiPiM2(0),
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
}//AliAnalysisTaskUpcEtaCAWP

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::Init()
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
AliAnalysisTaskUpcEtaCAWP::~AliAnalysisTaskUpcEtaCAWP() 
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
}//~AliAnalysisTaskUpcEtaCAWP

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::UserCreateOutputObjects() //use the names defined in the constructor to makes specific data structures
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

	TString CutNameEtaC[5] = { "Analyzed","Triggered","#eta_{C} Charge/nTracks Match","#eta_{C} Channel Match","#eta_{C} Candidates"};
	fHistNeventsEtaC = new TH1D("Events Histogram", "Collision Events at Each Analysis Stage", 5, 0.5, 5.5);
	for (Int_t i = 0; i < 5; i++) fHistNeventsEtaC->GetXaxis()->SetBinLabel(i + 1, CutNameEtaC[i].Data());
	fListHist->Add(fHistNeventsEtaC);

	TString ChannelNames[9] = { "K* K*","K* K^{+} #pi^{-} + c.c.","K^{+} #pi^{-} K^{-} #pi^{+} + c.c.","#rho^{0} #rho^{0}","2(#pi^{+} #pi^{-})","3(#pi^{+} #pi^{-})",
		"2(K^{+} K^{-})","K^{0} K^{-} #pi^{+} #pi^{-} #pi^{+} + c.c.","K^{+} K^{-} 2(#pi^{+} #pi^{-})" };
	fEtaCCandidatesPerChannel = new TH1D("Candidates Per Channel", "Events Per Decay Channel", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fEtaCCandidatesPerChannel->GetXaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fEtaCCandidatesPerChannel);

	fEtaCLowPtCandidatesPerChannel = new TH1D("Candidates Per Channel Low Pt", "Events Per Decay Channel, Low Pt", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fEtaCLowPtCandidatesPerChannel->GetXaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fEtaCLowPtCandidatesPerChannel);
	
	fAllPtVsMinvEtaC = new TH2D("Pt V Minv EtaC Candidates", "Transverse Momentum V. Invariant Mass, all Channel Matches", 300, 0., 6., 500, 0., 5.);
	fListHist->Add(fAllPtVsMinvEtaC);

	fAllMinvEtaCLowPt = new TH1D("Minv EtaC Candidates, Low Pt", "Invariant Mass, Low Pt Events", 300, 0., 6.);
	fListHist->Add(fAllMinvEtaCLowPt);

	fChannelVsMinvEtaC = new TH2D("Channel V Minv EtaC, Low Pt", "Decay Channel V. Invariant Mass, Low Pt Events", 300, 0., 6., 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fChannelVsMinvEtaC->GetYaxis()->SetBinLabel(i + 1, ChannelNames[i].Data());
	fListHist->Add(fChannelVsMinvEtaC);


	//##### Define Decay Channel histograms #####
	//##### Define k*(892) channel histograms.
	TString CutNameEtaCKstarChannel[9] = { "Two-K, Two-#pi","Like sign kaons","Like sign pions","Like sign both","Opposite sign (#eta_{C})",
		"2Kstar events","1Kstar events","0Kstar events","2+ tracks with pT>0.4" };

	fHistNeventsEtaCKstarChannel = new TH1D("EtaC Events Kstar Channel", "#eta_{C} Events at Each K* Channel Analysis Stage", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fHistNeventsEtaCKstarChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaCKstarChannel[i].Data());
	fListHistKstar->Add(fHistNeventsEtaCKstarChannel);

	TString KstarEventCandidatesNames[7] = { "Four good tracks","2-#pi 2-K 0-M/B","2-#pi 1-K 1-M/B","2-#pi 0-K 2-M/B",
		"1-#pi 2-K 1-M/B","0-#pi 2-K 2-M/B","1-#pi 1-K 2-M/B" };

	fKstarEventCandidates = new TH1D("Kstar Event Candidates", "Candidates for K* Channel Analysis", 7, 0.5, 7.5);
	for (Int_t i = 0; i < 7; i++) fKstarEventCandidates->GetXaxis()->SetBinLabel(i + 1, KstarEventCandidatesNames[i].Data());
	fListHistKstar->Add(fKstarEventCandidates);

	//A Dalitz plot. These are sometimes useful to look at in these types of studies. plots PiK pair masses
	fMPiKvsMPiK = new TH2D("fMPiKvsMPiK", "fMPiKvsMPiK", 300, 0.0, 0.3, 500, 0.0, 0.05);
	fListHistKstar->Add(fMPiKvsMPiK);
	fKstarMPiKvsMPiK = new TH2D("fKstarMPiKvsMPiK", "fKstarMPiKvsMPiK", 300, 0.0, 0.3, 500, 0.0, 0.05);
	fListHistKstar->Add(fKstarMPiKvsMPiK);

	//Define 2 k*(892) histograms.
	f2KstarPtVsMinvEtaC = new TH2D("f2KstarPtVsMinvEtaC", "f2KstarPtVsMinvEtaC", 1000, 0., 6., 1000, 0., 10.);
	fList2KStar->Add(f2KstarPtVsMinvEtaC);
	f2KstarPtVsMinvFirstKstar = new TH2D("f2KstarPtVsMinvFirstKstar", "f2KstarPtVsMinvFirstKstar", 600, 0., 6., 500, 0., 5.);
	fList2KStar->Add(f2KstarPtVsMinvFirstKstar);
	f2KstarPtVsMinvSecondKstar = new TH2D("f2KstarPtVsMinvSecondKstar", "f2KstarPtVsMinvSecondKstar", 600, 0., 6., 500, 0., 5.);
	fList2KStar->Add(f2KstarPtVsMinvSecondKstar);

	f2KstarPtPiPlus = new TH1D("f2KstarPtPiPlus", "f2KstarPtPiPlus", 300, 0., 3.);
	fList2KStarDiagnostic->Add(f2KstarPtPiPlus);
	f2KstarPtPiMinus = new TH1D("f2KstarPtPiMinus", "f2KstarPtPiMinus", 300, 0., 3.);
	fList2KStarDiagnostic->Add(f2KstarPtPiMinus);
	f2KstarPtKPlus = new TH1D("f2KstarPtKPlus", "f2KstarPtKPlus", 300, 0., 3.);
	fList2KStarDiagnostic->Add(f2KstarPtKPlus);
	f2KstarPtKMinus = new TH1D("f2KstarPtKMinus", "f2KstarPtKminus", 300, 0., 3.);
	fList2KStarDiagnostic->Add(f2KstarPtKMinus);
	f2KstarTPCsignalPion = new TH2D("f2KstarTPCsignalPion", "f2KstarTPCsignalPion", 600, 0., 300., 600, 0., 300.);
	fList2KStarDiagnostic->Add(f2KstarTPCsignalPion);
	f2KstarTPCsignalKaon = new TH2D("f2KstarTPCsignalKaon", "f2KstarTPCsignalKaon", 600, 0., 300., 600, 0., 300.);
	fList2KStarDiagnostic->Add(f2KstarTPCsignalKaon);
	f2KstarDedxVsPtPion = new TH2D("f2KstarDedxVsPtPion", "f2KstarDedxVsPtPion", 1000, 0., 10., 3000, 0., 300.);
	fList2KStarDiagnostic->Add(f2KstarDedxVsPtPion);
	f2KstarDedxVsPtKaon = new TH2D("f2KstarDedxVsPtKaon", "f2KstarDedxVsPtKaon", 1000, 0., 10., 3000, 0., 300.);
	fList2KStarDiagnostic->Add(f2KstarDedxVsPtKaon);


	f2KstarEtaVsMinvEtaC = new TH2D("f2KstarEtaVsMinvEtaC", "f2KstarEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList2KStarEtaC->Add(f2KstarEtaVsMinvEtaC);
	f2KstarEtaVsMinvEtaC400MeVPtMax = new TH2D("f2KstarEtaVsMinvEtaC400MeVPtMax", "f2KstarEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList2KStarEtaC->Add(f2KstarEtaVsMinvEtaC400MeVPtMax);
	f2KstarEtaVsMinvEtaC100MeVPtMax = new TH2D("f2KstarEtaVsMinvEtaC100MeVPtMax", "f2KstarEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList2KStarEtaC->Add(f2KstarEtaVsMinvEtaC100MeVPtMax);
	f2KstarSumPzVsMinvEtaC = new TH2D("f2KstarSumPzVsMinvEtaC", "f2KstarSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList2KStarEtaC->Add(f2KstarSumPzVsMinvEtaC);
	f2KstarScalarSumP = new TH1D("f2KstarScalarSumP", "f2KstarScalarSumP", 1000, 0., 10.);
	fList2KStarEtaC->Add(f2KstarScalarSumP);
	f2KstarVectorSumPt = new TH1D("f2KstarVectorSumPt", "f2KstarVectorSumPt", 1000, 0., 1.);
	fList2KStarEtaC->Add(f2KstarVectorSumPt);

	//Define 1 k*(892) channel histograms
	f1KstarPtVsMinvEtaC = new TH2D("f1KstarPtVsMinvEtaC", "f1KstarPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList1KStar->Add(f1KstarPtVsMinvEtaC);
	f1KstarPtVsMinvKstar = new TH2D("f1KstarPtVsMinvKstar", "f1KstarPtVsMinvKstar", 600, 0., 10., 500, 0., 5.);
	fList1KStar->Add(f1KstarPtVsMinvKstar);
	f1KstarPtVsMinvOtherPiKcombo = new TH2D("f1KstarPtVsMinvOtherPiKcombo", "f1KstarPtVsMinvOtherPiKcombo", 600, 0., 6., 500, 0., 5.);
	fList1KStar->Add(f1KstarPtVsMinvOtherPiKcombo);

	f1KstarPtPiPlus = new TH1D("f1KstarPtPiPlus", "f1KstarPtPiPlus", 300, 0., 3.);
	fList1KStarDiagnostic->Add(f1KstarPtPiPlus);
	f1KstarPtPiMinus = new TH1D("f1KstarPtPiMinus", "f1KstarPtPiMinus", 300, 0., 3.);
	fList1KStarDiagnostic->Add(f1KstarPtPiMinus);
	f1KstarPtKPlus = new TH1D("f1KstarPtKPlus", "f1KstarPtKPlus", 300, 0., 3.);
	fList1KStarDiagnostic->Add(f1KstarPtKPlus);
	f1KstarPtKMinus = new TH1D("f1KstarPtKMinus", "f1KstarPtKMinus", 300, 0., 3.);
	fList1KStarDiagnostic->Add(f1KstarPtKMinus);
	f1KstarTPCsignalPion = new TH2D("f1KstarTPCsignalPion", "f1KstarTPCsignalPion", 600, 0., 300., 600, 0., 300.);
	fList1KStarDiagnostic->Add(f1KstarTPCsignalPion);
	f1KstarTPCsignalKaon = new TH2D("f1KstarTPCsignalKaon", "f1KstarTPCsignalKaon", 600, 0., 300., 600, 0., 300.);
	fList1KStarDiagnostic->Add(f1KstarTPCsignalKaon);
	f1KstarDedxVsPtPion = new TH2D("f1KstarDedxVsPtPion", "f1KstarDedxVsPtPion", 1000, 0., 10., 3000, 0., 300.);
	fList1KStarDiagnostic->Add(f1KstarDedxVsPtPion);
	f1KstarDedxVsPtKaon = new TH2D("f1KstarDedxVsPtKaon", "f1KstarDedxVsPtKaon", 1000, 0., 10., 3000, 0., 300.);
	fList1KStarDiagnostic->Add(f1KstarDedxVsPtKaon);

	f1KstarEtaVsMinvEtaC = new TH2D("f1KstarEtaVsMinvEtaC", "f1KstarEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList1KStarEtaC->Add(f1KstarEtaVsMinvEtaC);
	f1KstarEtaVsMinvEtaC400MeVPtMax = new TH2D("f1KstarEtaVsMinvEtaC400MeVPtMax", "f1KstarEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList1KStarEtaC->Add(f1KstarEtaVsMinvEtaC400MeVPtMax);
	f1KstarEtaVsMinvEtaC100MeVPtMax = new TH2D("f1KstarEtaVsMinvEtaC100MeVPtMax", "f1KstarEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList1KStarEtaC->Add(f1KstarEtaVsMinvEtaC100MeVPtMax);
	f1KstarSumPzVsMinvEtaC = new TH2D("f1KstarSumPzVsMinvEtaC", "f1KstarSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList1KStarEtaC->Add(f1KstarSumPzVsMinvEtaC);
	f1KstarScalarSumP = new TH1D("f1KstarScalarSumP", "f1KstarScalarSumP", 1000, 0., 10.);
	fList1KStarEtaC->Add(f1KstarScalarSumP);
	f1KstarVectorSumPt = new TH1D("f1KstarVectorSumPt", "f1KstarVectorSumPt", 1000, 0., 1.);
	fList1KStarEtaC->Add(f1KstarVectorSumPt);

	//Define 0 k*(892) channel histograms
	f0KstarPtVsMinvEtaC = new TH2D("f0KstarPtVsMinvEtaC", "f0KstarPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList0KStar->Add(f0KstarPtVsMinvEtaC);
	f0KstarPtVsMinvFirstPiKcombo = new TH2D("f0KstarPtVsMinvFirstPiKcombo", "f0KstarPtVsMinvFirstPiKcombo", 600, 0., 6., 500, 0., 5.);
	fList0KStar->Add(f0KstarPtVsMinvFirstPiKcombo);
	f0KstarPtVsMinvSecondPiKcombo = new TH2D("f0KstarPtVsMinvSecondPiKcombo", "f0KstarPtVsMinvSecondPiKcombo", 600, 0., 6., 500, 0., 5.);
	fList0KStar->Add(f0KstarPtVsMinvSecondPiKcombo);

	f0KstarPtPiPlus = new TH1D("f0KstarPtPiPlus", "f0KstarPtPiPlus", 300, 0., 3.);
	fList0KStarDiagnostic->Add(f0KstarPtPiPlus);
	f0KstarPtPiMinus = new TH1D("f0KstarPtPiMinus", "fK0starPtPiMinus", 300, 0., 3.);
	fList0KStarDiagnostic->Add(f0KstarPtPiMinus);
	f0KstarPtKPlus = new TH1D("f0KstarPtKPlus", "f0KstarPtKPlus", 300, 0., 3.);
	fList0KStarDiagnostic->Add(f0KstarPtKPlus);
	f0KstarPtKMinus = new TH1D("f0KstarPtKMinus", "f0KstarPtKMinus", 300, 0., 3.);
	fList0KStarDiagnostic->Add(f0KstarPtKMinus);
	f0KstarTPCsignalPion = new TH2D("f0KstarTPCsignalPion", "f0KstarTPCsignalPion", 600, 0., 300., 600, 0., 300.);
	fList0KStarDiagnostic->Add(f0KstarTPCsignalPion);
	f0KstarTPCsignalKaon = new TH2D("f0KstarTPCsignalKaon", "f0KstarTPCsignalKaon", 600, 0., 300., 600, 0., 300.);
	fList0KStarDiagnostic->Add(f0KstarTPCsignalKaon);
	f0KstarDedxVsPtPion = new TH2D("f0KstarDedxVsPtPion", "f0KstarDedxVsPtPion", 1000, 0., 10., 3000, 0., 300.);
	fList0KStarDiagnostic->Add(f0KstarDedxVsPtPion);
	f0KstarDedxVsPtKaon = new TH2D("f0KstarDedxVsPtKaon", "f0KstarDedxVsPtKaon", 1000, 0., 10., 3000, 0., 300.);
	fList0KStarDiagnostic->Add(f0KstarDedxVsPtKaon);


	f0KstarEtaVsMinvEtaC = new TH2D("f0KstarEtaVsMinvEtaC", "f0KstarEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList0KStarEtaC->Add(f0KstarEtaVsMinvEtaC);
	f0KstarEtaVsMinvEtaC400MeVPtMax = new TH2D("f0KstarEtaVsMinvEtaC400MeVPtMax", "f0KstarEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList0KStarEtaC->Add(f0KstarEtaVsMinvEtaC400MeVPtMax);
	f0KstarEtaVsMinvEtaC100MeVPtMax = new TH2D("f0KstarEtaVsMinvEtaC100MeVPtMax", "f0KstarEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList0KStarEtaC->Add(f0KstarEtaVsMinvEtaC100MeVPtMax);
	f0KstarSumPzVsMinvEtaC = new TH2D("f0KstarSumPzVsMinvEtaC", "f0KstarSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList0KStarEtaC->Add(f0KstarSumPzVsMinvEtaC);
	f0KstarScalarSumP = new TH1D("f0KstarScalarSumP", "f0KstarScalarSumP", 1000, 0., 10.);
	fList0KStarEtaC->Add(f0KstarScalarSumP);
	f0KstarVectorSumPt = new TH1D("f0KstarVectorSumPt", "f0KstarVectorSumPt", 1000, 0., 1.);
	fList0KStarEtaC->Add(f0KstarVectorSumPt);

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

	//##### Define K+,K-,2(Pi+,Pi-) channel histograms.
	TString CutNameEtaC2K4PiChannel[7] = { "Four Pions Two Kaons","non-zero net charge","zero net charge (#eta_{C})","2+ tracks with pT>0.4" };
	
	fHistNeventsEtaC2K4PiChannel = new TH1D("EtaC Events 2K4Pi Channel", "#eta_{C} Events at Each K^{+} K^{-} 2(Pi^{+} Pi^{-}) Channel Analysis Stage", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) fHistNeventsEtaC2K4PiChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaC2K4PiChannel[i].Data());
	fListHistK0s3PiPi4K->Add(fHistNeventsEtaC2K4PiChannel);

	TString TwoKFourPiEventCandidatesNames[8] = { "Six good tracks","2-K 4-Pi 0-M/B","1-K 4-Pi 1-M/B","0-K 4-Pi 2-M/B","2-K 3-Pi 1-M/B","2-K 2-Pi 2-M/B","2-K 1-Pi 3-M/B","1-K 3-Pi 2-M/B" };

	f2K4PiEventCandidates = new TH1D("Two K Four Pi Event Candidates", "Candidates for K^{+} K^{-} 2(Pi^{+} Pi^{-}) Channel Analysis", 8, 0.5, 8.5);
	for (Int_t i = 0; i < 8; i++) f2K4PiEventCandidates->GetXaxis()->SetBinLabel(i + 1, TwoKFourPiEventCandidatesNames[i].Data());
	fList2K4Pi->Add(f2K4PiEventCandidates);
	
	f2K4PiPtVsMinvEtaC = new TH2D("Pt V Minv EtaC", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Pt V. Minv, #eta_{C}", 600, 0., 6., 500, 0., 5.);
	fList2K4Pi->Add(f2K4PiPtVsMinvEtaC);

	f2K4PiEtaVsMinvEtaC = new TH2D("Eta V Minv EtaC", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Eta V. Minv, #eta_{C}", 1000, 0., 10., 2000, -10., 10.);
	fList2K4PiEtaC->Add(f2K4PiEtaVsMinvEtaC);
	f2K4PiEtaVsMinvEtaC400MeVPtMax = new TH2D("Eta V Minv EtaC 0.4 Pt Max", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Eta V. Minv, #eta_{C}, 0.4 GeV/c Pt_{max}", 1000, 0., 10., 2000, -10., 10.);
	fList2K4PiEtaC->Add(f2K4PiEtaVsMinvEtaC400MeVPtMax);
	f2K4PiEtaVsMinvEtaC100MeVPtMax = new TH2D("Eta V Minv EtaC 0.1 Pt Max", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Eta V. Minv, #eta_{C}, 0.1 GeV/c Pt_{max}", 1000, 0., 10., 2000, -10., 10.);
	fList2K4PiEtaC->Add(f2K4PiEtaVsMinvEtaC100MeVPtMax);
	f2K4PiSumPzVsMinvEtaC = new TH2D("Sum Pz V Minv EtaC", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Sum Pz V. Minv, #eta_{C}", 1000, 0., 10., 2000, -10., 10.);
	fList2K4PiEtaC->Add(f2K4PiSumPzVsMinvEtaC);
	f2K4PiScalarSumP = new TH1D("Scalar Sum P", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Scalar Sum P, #eta_{C}", 1000, 0., 10.);
	fList2K4PiEtaC->Add(f2K4PiScalarSumP);
	f2K4PiVectorSumPt = new TH1D("Vector Sum Pt", "K^{+} K^{-} 2(Pi^{+} Pi^{-}) Vector Sum Pt, #eta_{C}", 1000, 0., 1.);
	fList2K4PiEtaC->Add(f2K4PiVectorSumPt);

	//##### Define 2 rho(0) channel histograms.
	TString CutNameEtaCRhoChannel[7] = { "Four Pions","non-zero net charge","zero net charge (#eta_{C})","two sets of #rho's",
		"one set of #rho's","no sets of #rho's","2+ tracks with pT>0.4" }; 
	
	//Create a histotram to count events after five event cuts and with different event/track combinations. 
	fHistNeventsEtaCRhoChannel = new TH1D("EtaC Events Rho Channel", "#eta_{C} Events at Each 2(#pi^{+} #pi^{-}) Channel Analysis Stage", 7, 0.5, 7.5);
	for (Int_t i = 0; i < 7; i++) fHistNeventsEtaCRhoChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaCRhoChannel[i].Data());
	fListHist2Rho4Pion->Add(fHistNeventsEtaCRhoChannel);

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

	

	//##### Define 3PiPi channel histograms
	TString CutNameEtaC3PiPiChannel[4] = { "Six Pions","non-zero net charge","zero net charge (#eta_{C})","2+ tracks with pT>0.4" };

	fHistNeventsEtaC3PiPiChannel = new TH1D("EtaC Events 3PiPi Channel", "#eta_{C} Events at Each 3(#pi^{+} #pi^{-}) Channel Analysis Stage", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) fHistNeventsEtaC3PiPiChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaC3PiPiChannel[i].Data());
	fListHistK0s3PiPi4K->Add(fHistNeventsEtaC3PiPiChannel);

	TString SixPiEventCandidatesNames[5] = { "Six good tracks","6-#pi 0-M/B","5-#pi 1-M/B","4-#pi 2-M/B","3-#pi 3-M/B" };

	f6PiEventCandidates = new TH1D("Three PiPi Event Candidates", "Candidates for 3(#pi^{+} #pi^{-}) Channel Analysis", 5, 0.5, 5.5);
	for (Int_t i = 0; i < 5; i++) f6PiEventCandidates->GetXaxis()->SetBinLabel(i + 1, SixPiEventCandidatesNames[i].Data());
	fList3PiPi->Add(f6PiEventCandidates);

	f3PiPiPtVsMinvEtaC = new TH2D("f3PiPiPtVsMinvEtaC", "f3PiPiPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList3PiPi->Add(f3PiPiPtVsMinvEtaC);
	f3PiPiEtaVsMinvEtaC = new TH2D("f3PiPiEtaVsMinvEtaC", "f3PiPiEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList3PiPiEtaC->Add(f3PiPiEtaVsMinvEtaC);
	f3PiPiEtaVsMinvEtaC400MeVPtMax = new TH2D("f3PiPiEtaVsMinvEtaC400MeVPtMax", "f3PiPiEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList3PiPiEtaC->Add(f3PiPiEtaVsMinvEtaC400MeVPtMax);
	f3PiPiEtaVsMinvEtaC100MeVPtMax = new TH2D("f3PiPiEtaVsMinvEtaC100MeVPtMax", "f3PiPiEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList3PiPiEtaC->Add(f3PiPiEtaVsMinvEtaC100MeVPtMax);
	f3PiPiSumPzVsMinvEtaC = new TH2D("f3PiPiSumPzVsMinvEtaC", "f3PiPiSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList3PiPiEtaC->Add(f3PiPiSumPzVsMinvEtaC);
	f3PiPiScalarSumP = new TH1D("f3PiPiScalarSumP", "f3PiPiScalarSumP", 1000, 0., 10.);
	fList3PiPiEtaC->Add(f3PiPiScalarSumP);
	f3PiPiVectorSumPt = new TH1D("f3PiPiVectorSumPt", "f3PiPiVectorSumP", 1000, 0., 1.);
	fList3PiPiEtaC->Add(f3PiPiVectorSumPt);

	//##### Define 4 kaon channel histograms
	TString CutNameEtaC4KaonChannel[4] = { "Four Kaons","non-zero net charge","zero net charge (#eta_{C})","2+ tracks with pT>0.4" };

	fHistNeventsEtaC4KaonChannel = new TH1D("EtaC Events 2KK Channel", "#eta_{C} Events at Each 2(K^{+} K^{-}) Channel Analysis Stage", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) fHistNeventsEtaC4KaonChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaC4KaonChannel[i].Data());
	fListHistK0s3PiPi4K->Add(fHistNeventsEtaC4KaonChannel);

	TString FourKEventCandidatesNames[4] = { "four good tracks","4-K 0-M/B","3-K 1-M/B","2-K 2-M/B" };

	f4KEventCandidates = new TH1D("Two KK Event Candidates", "Candidates for 2(K^{+} K^{-}) Channel Analysis", 4, 0.5, 4.5);
	for (Int_t i = 0; i < 4; i++) f4KEventCandidates->GetXaxis()->SetBinLabel(i + 1, FourKEventCandidatesNames[i].Data());
	fList4K->Add(f4KEventCandidates);
	
	f4KaonPtVsMinvEtaC = new TH2D("f4KaonPtVsMinvEtaC", "f4KaonPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fList4K->Add(f4KaonPtVsMinvEtaC);
	f4KaonPtVsMinvKK = new TH2D("f4KaonPtVsMinvKK", "f4KaonPtVsMinvKK", 600, 0., 6., 500, 0., 5.);
	fList4KDiagnostic->Add(f4KaonPtVsMinvKK);
	f4KVs2KMinv = new TH2D("f4KVs2KMinv", "f4KVs2KMinv", 600, 0., 6., 500, 0., 5.);
	fList4KDiagnostic->Add(f4KVs2KMinv);
	fM2KKVsM2KK = new TH2D("fM2KKVsM2KK", "fM2KKVsM2KK", 600, 0., 6., 500, 0., 5.);
	fList4KDiagnostic->Add(fM2KKVsM2KK);
	f4KVs2KMinvSquared = new TH2D("f4KVs2KMinvSquared", "f4KVs2KMinvSquared", 2500, 0., 25., 2500, 0., 25.);
	fList4KDiagnostic->Add(f4KVs2KMinvSquared);
	f4KaonEtaVsMinvEtaC = new TH2D("f4KaonEtaVsMinvEtaC", "f4KaonEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList4KEtaC->Add(f4KaonEtaVsMinvEtaC);
	f4KaonEtaVsMinvEtaC400MeVPtMax = new TH2D("f4KaonEtaVsMinvEtaC400MeVPtMax", "f4KaonEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList4KEtaC->Add(f4KaonEtaVsMinvEtaC400MeVPtMax);
	f4KaonEtaVsMinvEtaC100MeVPtMax = new TH2D("f4KaonEtaVsMinvEtaC100MeVPtMax", "f4KaonEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fList4KEtaC->Add(f4KaonEtaVsMinvEtaC100MeVPtMax);
	f4KaonSumPzVsMinvEtaC = new TH2D("f4KaonSumPzVsMinvEtaC", "f4KaonSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fList4KEtaC->Add(f4KaonSumPzVsMinvEtaC);
	f4KaonScalarSumP = new TH1D("f4KaonScalarSumP", "f4KaonScalarSumP", 1000, 0., 10.);
	fList4KEtaC->Add(f4KaonScalarSumP);
	f4KaonVectorSumPt = new TH1D("f4KaonVectorSumPt", "f4KaonVectorSumPt", 1000, 0., 1.);
	fList4KEtaC->Add(f4KaonVectorSumPt);



	//##### Define K0 channel Histograms
	TString CutNameEtaCK0sChannel[9] = { "4-6 Good Tracks & >0 K^{0}s","Six Good Enough Tracks","1-K 3-#pi 1-K^{0}s","negative total charge","positive total charge",
		"K^{-} and 0 total charge","K^{+}and 0 total charge","Total with correct charges (#eta_{C})","2+ tracks with pT>0.4" };

	fHistNeventsEtaCK0sChannel = new TH1D("EtaC Events K0short Channel", "#eta_{C} Events at Each K^{0} short Channel Analysis Stage", 9, 0.5, 9.5);
	for (Int_t i = 0; i < 9; i++) fHistNeventsEtaCK0sChannel->GetXaxis()->SetBinLabel(i + 1, CutNameEtaCK0sChannel[i].Data());
	fListHistK0s3PiPi4K->Add(fHistNeventsEtaCK0sChannel);

	TString K0sEventCandidateNames[8] = { "Six Good Enough Tracks","Opp sign K^{0}s daughters","1-K 3-#pi 1-K^{0}s 0-M/B","1-K 2-#pi 1-K^{0}s 1-M/B",
		"0-K 3-#pi 1-K^{0}s 1-M/B","1-K 1-#pi 1-K^{0}s 2-M/B","0-K 2-#pi 1-K^{0}s 2-M/B","1-K 1-#pi 1-K^{0}s 3-M/B" };

	fK0sEventCandidates = new TH1D("K0 short Event Candidates", "Candidates for K0 Short Channel Analysis ", 8, 0.5, 8.5);
	for (Int_t i = 0; i < 8; i++) fK0sEventCandidates->GetXaxis()->SetBinLabel(i + 1, K0sEventCandidateNames[i].Data());
	fListK0Short->Add(fK0sEventCandidates);

	fHistK0sCandidatesPerEvent = new TH1D("K0 short Per Event", "Number of K0 short Candidates per Event", 10, -0.5, 9.5);
	fListK0Short->Add(fHistK0sCandidatesPerEvent);
	fK0sPtVsMinvEtaC = new TH2D("fK0sPtVsMinvEtaC", "fK0sPtVsMinvEtaC", 600, 0., 6., 500, 0., 5.);
	fListK0Short->Add(fK0sPtVsMinvEtaC);
	fK0sPtVsMinvK0s = new TH2D("fK0sPtVsMinvK0s", "fK0sPtVsMinvK0s", 600, 0., 6., 500, 0., 5.);
	fListK0Short->Add(fK0sPtVsMinvK0s);
	fK0sPtVsMinvKPi = new TH2D("fK0sPtVsMinvKPi", "fK0sPtVsMinvKPi", 600, 0., 6., 500, 0., 5.);
	fListK0Short->Add(fK0sPtVsMinvKPi);

	fV0DecayLength = new TH1D("V0 Decay Length", "V0 Decay Length", 200, 0., 200.);
	fListK0ShortDiagnostic->Add(fV0DecayLength);
	fV0Eta = new TH1D("V0 Eta", "V0 Eta", 160, -0.8, 0.8);
	fListK0ShortDiagnostic->Add(fV0Eta);
	fV0CosPointingAngle = new TH1D("V0 Cos Pointing Angle", "V0 Cos Pointing Angle", 500, 0., 1.);
	fListK0ShortDiagnostic->Add(fV0CosPointingAngle);
	fV0sMassK0s = new TH1D("V0 Mass K0 short Hypothesis", "V0 Mass K0 short Hypothesis", 600, 0.2, 0.8);
	fListK0ShortDiagnostic->Add(fV0sMassK0s);
	fK0DaughterDca = new TH1D("K0s Daughter Dca", "K0s Daughter Dca", 200, 0., 20.);
	fListK0ShortDiagnostic->Add(fK0DaughterDca);
	fK0sDcaToPrimVertex = new TH1D("K0s Dca To Prim Vertex", "K0s Dca To Prim Vertex", 200, 0., 20.);
	fListK0ShortDiagnostic->Add(fK0sDcaToPrimVertex);
	fK0sDaughterDcaToPrimVertex = new TH1D("K0s Daughter Dca To PrimVertex", "K0s Daughter Dca To Prim Vertex", 200, 0., 20.);
	fListK0ShortDiagnostic->Add(fK0sDaughterDcaToPrimVertex);
	fK0sMassDistribution = new TH1D("K0s Mass Distribution", "Mass Distribution, K^{0} Short Hypothesis", 600, 0.2, 0.8);
	fListK0ShortDiagnostic->Add(fK0sMassDistribution);
	fHistNProngFound = new TH1D("nProngFound Per Event", "nProngFound Per Event", 10, -0.5, 9.5);
	fListK0ShortDiagnostic->Add(fHistNProngFound);
	fK0GoodTracks = new TH1D("K0 Good tracks", "number of K0 tracks which are 'Good'", 10, -0.5, 9.5);
	fListK0ShortDiagnostic->Add(fK0GoodTracks);
	fHistNK0sPion = new TH1D("nK0sPion Per Event", "nK0sPion Per Event (should always be 2)", 10, -0.5, 9.5);
	fListK0ShortDiagnostic->Add(fHistNK0sPion);
	fK0sEtaVsMinvEtaC = new TH2D("fK0sEtaVsMinvEtaC", "fK0sEtaVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fListK0ShortEtaC->Add(fK0sEtaVsMinvEtaC);
	fK0sEtaVsMinvEtaC400MeVPtMax = new TH2D("fK0sEtaVsMinvEtaC400MeVPtMax", "fK0sEtaVsMinvEtaC400MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fListK0ShortEtaC->Add(fK0sEtaVsMinvEtaC400MeVPtMax);
	fK0sEtaVsMinvEtaC100MeVPtMax = new TH2D("fK0sEtaVsMinvEtaC100MeVPtMax", "fK0sEtaVsMinvEtaC100MeVPtMax", 1000, 0., 10., 2000, -10., 10.);
	fListK0ShortEtaC->Add(fK0sEtaVsMinvEtaC100MeVPtMax);
	fK0sSumPzVsMinvEtaC = new TH2D("fK0sSumPzVsMinvEtaC", "fK0sSumPzVsMinvEtaC", 1000, 0., 10., 2000, -10., 10.);
	fListK0ShortEtaC->Add(fK0sSumPzVsMinvEtaC);
	fK0sScalarSumP = new TH1D("fK0sScalarSumP", "fK0sScalarSumP", 1000, 0., 10.);
	fListK0ShortEtaC->Add(fK0sScalarSumP);
	fK0sVectorSumPt = new TH1D("fK0sVectorSumPt", "fK0sVectorSumPt", 1000, 0., 1.);
	fListK0ShortEtaC->Add(fK0sVectorSumPt);
	fK0sDecayLength = new TH1D("fK0sDecayLength", "fK0sDecayLength", 300, 0., 30.);
	fListK0ShortDiagnostic->Add(fK0sDecayLength);
	fK0sPosDaughterPt = new TH1D("fK0sPosDaughterPt", "fK0sPosDaughterPt", 1000, 0., 10.);
	fListK0ShortDiagnostic->Add(fK0sPosDaughterPt);
	fK0sNegDaughterPt = new TH1D("fK0sNegDaughterPt", "fK0sNegDaughterPt", 1000, 0., 10.);
	fListK0ShortDiagnostic->Add(fK0sNegDaughterPt);
	fK0sPosVsNegDaughterPt = new TH2D("fK0sPosVsNegDaughterPt", "fK0sPosVsNegDaughterPt", 1000, 0., 10., 100, 0., 10.);
	fListK0ShortDiagnostic->Add(fK0sPosVsNegDaughterPt);
	fK0sPionPt = new TH1D("fK0sPionPt", "fK0sPionPt", 1000, 0., 10.);
	fListK0ShortDiagnostic->Add(fK0sPionPt);
	fK0sKaonPt = new TH1D("fK0sKaonPt", "fK0sKaonPt", 1000, 0., 10.);
	fListK0ShortDiagnostic->Add(fK0sKaonPt);
	fK0sM2K0sVsM2KPi = new TH2D("fK0sM2K0sVsM2KPi", "fK0sM2K0sVsM2KPi", 600, 0., 6., 500, 0., 5.);
	fListK0ShortDiagnostic->Add(fK0sM2K0sVsM2KPi);
	fK0sM2K0sPiVsM2KPi = new TH2D("fK0sM2K0sPiVsM2KPi", "fK0sM2K0sPiVsM2KPi", 600, 0., 6., 500, 0., 5.);
	fListK0ShortDiagnostic->Add(fK0sM2K0sPiVsM2KPi);
	fK0sM2K0sKVsM2KPi = new TH2D("fK0sM2K0sKVsM2KPi", "fK0sM2K0sKVsM2KPi", 600, 0., 6., 500, 0., 5.);
	fListK0ShortDiagnostic->Add(fK0sM2K0sKVsM2KPi);

	fK0sTOFbetaVsPAll = new TH2D("beta V P All", "#beta V P All K0s Channel", 1000, 0., 10., 200, 0., 2.);
	fListK0ShortPID->Add(fK0sTOFbetaVsPAll);
	fK0sTOFbetaVsPPion = new TH2D("beta V P Pion", "#beta V P Pion K0s Channel", 1000, 0., 10., 200, 0., 2.);
	fListK0ShortPID->Add(fK0sTOFbetaVsPPion);
	fK0sDedxVsPAll = new TH2D("Dedx V P All", "Dedx V P All K0s Channel", 1000, 0., 10., 3000, 0., 300.);
	fListK0ShortPID->Add(fK0sDedxVsPAll);
	fK0sDedxVsPPion = new TH2D("Dedx V P Pion", "Dedx V P Pion K0s Channel", 1000, 0., 10., 3000, 0., 300.);
	fListK0ShortPID->Add(fK0sDedxVsPPion);

	InitSystematics();

	//cout << "##### End of UserCreateOutputObjects()" << endl;

	PostData(1, fListTrig);
	PostData(2, fListHist);
	PostData(3, fListHistKstar);
	PostData(4, fListHist2Rho4Pion);
	PostData(5, fListHistK0s3PiPi4K);
}//UserCreateOutputObjects

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::InitSystematics()
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
void AliAnalysisTaskUpcEtaCAWP::UserExec(Option_t *) 
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
void AliAnalysisTaskUpcEtaCAWP::RunAODtrig()
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
void AliAnalysisTaskUpcEtaCAWP::RunAODhist()
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

	Double_t etaCMass = 2.984;
	Double_t etaCWidth = 0.032;

	//##### Get Event #####
	AliAODEvent *aod = (AliAODEvent*)InputEvent();
	if (!aod) return; //if return condition is met, stop analyzing this event, and next event is analyzed from the top of RunAODhist()
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

	if (fVtxContrib < 2) return;

	//#####v0 decision cut -- V0A and V0C must be empty to continue
	AliAODVZERO *fV0data = aod->GetVZEROData();
	fV0Adecision = fV0data->GetV0ADecision();
	fV0Cdecision = fV0data->GetV0CDecision();

	if (fV0Adecision != AliAODVZERO::kV0Empty || fV0Cdecision != AliAODVZERO::kV0Empty) return;

	//#####neutron ZDC cut -- must have <8 neutrons
	AliAODZDC *fZDCdata = aod->GetZDCData();
	fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
	fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
	fZDCAtime = fZDCdata->GetZNATime();
	fZDCCtime = fZDCdata->GetZNCTime();

	fHistZDCAenergy->Fill(fZNAenergy);
	fHistZDCCenergy->Fill(fZNCenergy);
	fHistZDCAtime->Fill(fZDCAtime);
	fHistZDCCtime->Fill(fZDCCtime);
	fHistZDCImpactParameter->Fill(fZDCdata->GetImpactParameter());
	fHistZDCAImpactParameter->Fill(fZDCdata->GetImpactParamSideA());
	fHistZDCCImpactParameter->Fill(fZDCdata->GetImpactParamSideC());

	if (trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1); //used to compair the number of neutrons per event, binned by number of neutrons
	if (fZNAenergy < 11500 && fZNCenergy < 11500) fHistZDCCuts->Fill(2); //was 8200, fills events with low ZNA and ZNC energy, <8 neutrons
	if (fZNAenergy < 1500 && fZNCenergy < 1500) fHistZDCCuts->Fill(3); //was 683, fills events with very low ZNA and ZNC energy, 0 neutrons
	if (fZDCAtime == 0 && fZDCCtime == 0) fHistZDCCuts->Fill(4);

	if (fZNAenergy > 11500 || fZNCenergy > 11500) return; // was >8200 in 2011 code, >1500 for 0n0n or <= 1500 for XnXn; one ZDC <= 1500 AND other ZDC > 1500 for 0nXn or Xn0n

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

	//##### Diagnostic figures regarding nSigma TOF vs TPC, etc. #####
	TLorentzVector vPion[6], vKaon[6];
	Short_t qKaon[6], qPion[6];
	UInt_t nKaon = 0, nPion = 0, nBoth = 0, nNeither = 0;
	UInt_t nPreKaon = 0, nPrePion = 0; //added by Alec
	Double_t fRecTPCsignalPion[6], fRecTPCsignalKaon[6];
	Float_t eventStartTime = 0.;
	Double_t beta = 0;
	Int_t nTracksWithoutTOFinfo = 0;
	Double_t PCutLow = 0.6, PCutHigh = 3.;
	Int_t qpos = 0, qneg = 0;

	if ((nGoodTracks == 6 || nGoodTracks == 4) && nSpdHits > 1) {
		for (Int_t i = 0; i < nGoodTracks; i++) {
			beta = -999.;
			AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
			if (!trk) AliFatal("Not a standard AOD");
			//calculate net charge of good tracks -- count number of positive and number of negative as integers so we don't need to add shorts
			if (trk->Charge() > 0) {
				qpos++;
			} 
			else if (trk->Charge() < 0) {
				qneg++;
			}

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

			//nSigmaTPC/TOF vs nSigmaTPC/TOF for all particles with Pion mass hypothesis and Kaon mass hypothesis.      
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
			//TPC signal (dE/dx) vs TOF signal (time in ps) or TOF beta
			eventStartTime = aod->GetTOFHeader()->GetDefaultEventTimeVal(); //may have to define this differently. An array of start times is returned
			if (fPIDTOFPion[i] != -999.) fTOFIntegratedLength->Fill(trk->GetIntegratedLength());
			if (trk->GetIntegratedLength() >= 360. && trk->GetIntegratedLength() <= 800. && trk->GetTOFsignal() > 0. && eventStartTime < 999990.0) {
				beta = (trk->GetIntegratedLength()*0.01) / ((trk->GetTOFsignal() - eventStartTime)*(1.e-12)*TMath::C()); //0.01 converts cm to m, 1e-12 converts ps to sec

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
			//Identify tracks which are likely to be a pion
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

			if (fPIDTOFPion[i] == -999.) nTracksWithoutTOFinfo++; //Track if missing TOF signal
		}

		if (nGoodTracks == 4) {
			fHistFourTracksNpion->Fill(nPion);
			fHistFourTracksNkaon->Fill(nKaon);
			fFourTrackMissing->Fill(nTracksWithoutTOFinfo);
			fHistFourTracksNboth->Fill(nBoth);
			fHistFourTracksNneither->Fill(nNeither);

			//diagnostic histograms -- events in these bins should match with number of events analyzed
			f4KEventCandidates->Fill(1);
			if (nKaon == 4) f4KEventCandidates->Fill(2);
			if (nKaon == 3 && nBoth == 1) f4KEventCandidates->Fill(3);
			if (nKaon == 2 && nBoth == 2) f4KEventCandidates->Fill(4);

			f4PiEventCandidates->Fill(1);
			if (nPion == 4) f4PiEventCandidates->Fill(2);
			if (nPion == 3 && nBoth == 1) f4PiEventCandidates->Fill(3);
			if (nPion == 2 && nBoth == 2) f4PiEventCandidates->Fill(4);

			fKstarEventCandidates->Fill(1);
			if (nPion == 2 && nKaon == 2) fKstarEventCandidates->Fill(2);
			if (nPion == 2 && nKaon == 1 && nBoth == 1) fKstarEventCandidates->Fill(3);
			if (nPion == 2 && nKaon == 0 && nBoth == 2) fKstarEventCandidates->Fill(4);
			if (nPion == 1 && nKaon == 2 && nBoth == 1) fKstarEventCandidates->Fill(5);
			if (nPion == 0 && nKaon == 2 && nBoth == 2) fKstarEventCandidates->Fill(6);
			if (nPion == 1 && nKaon == 1 && nBoth == 2) fKstarEventCandidates->Fill(7);
		}
		else if (nGoodTracks == 6) {
			fHistSixTracksNpion->Fill(nPion);
			fHistSixTracksNkaon->Fill(nKaon);
			fSixTrackMissing->Fill(nTracksWithoutTOFinfo);
			fHistSixTracksNboth->Fill(nBoth);
			fHistSixTracksNneither->Fill(nNeither);

			f6PiEventCandidates->Fill(1);
			if (nPion == 6) f6PiEventCandidates->Fill(2);
			if (nPion == 5 && nBoth == 1) f6PiEventCandidates->Fill(3);
			if (nPion == 4 && nBoth == 2) f6PiEventCandidates->Fill(4);
			if (nPion == 3 && nBoth == 3) f6PiEventCandidates->Fill(5);

			f2K4PiEventCandidates->Fill(1);
			if (nPion == 4 && nKaon == 2) f2K4PiEventCandidates->Fill(2);
			if (nPion == 4 && nKaon == 1 && nBoth == 1) f2K4PiEventCandidates->Fill(3);
			if (nPion == 4 && nKaon == 0 && nBoth == 2) f2K4PiEventCandidates->Fill(4);
			if (nPion == 3 && nKaon == 2 && nBoth == 1) f2K4PiEventCandidates->Fill(5);
			if (nPion == 2 && nKaon == 2 && nBoth == 2) f2K4PiEventCandidates->Fill(6);
			if (nPion == 1 && nKaon == 2 && nBoth == 3) f2K4PiEventCandidates->Fill(7);
			if (nPion == 3 && nKaon == 1 && nBoth == 2) f2K4PiEventCandidates->Fill(8);
		}
		//PID check counts total identified tracks -- should always be four or six
		fHistPIDCheck->Fill(nPion + nKaon + nBoth + nNeither);
		nPreKaon = nKaon; //diagnostic, makes sure nKaon and nPion don't change because of assuming missing/both tracks
		nPrePion = nPion;
	}
	//if (nSpdHits > 1) cout << "PID block Completed. Good Tracks: " << nGoodTracks << endl;
	if (nGoodTracks < 4 || nGoodTracks > 6 || nSpdHits < 2) return;
	if ((qpos - qneg) != 0) return;
	if (nGoodTracks == 4 || nGoodTracks == 6 || nSpdHits > 1) fHistNeventsEtaC->Fill(3);
	
	//##### Analysis by Decay Channel #####
	Int_t nHighPtTracks = 0;
	Double_t SumPz = 0, VectorSumPt = 0, ScalarSumP = 0;
	TLorentzVector vCandidate;
	TVector3 sumPtVector;
	UInt_t newKaon = 0, newPion = 0;  //Added by Alec

	//##### 2Kstar, Kstar K+ Pi-, and K+ K- Pi+ Pi- channels #####
	TLorentzVector vKstar[2];
	Bool_t goodPairA = kFALSE;
	Bool_t goodPairB = kFALSE;
	Double_t boostInfoA[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };
	Double_t boostInfoB[13] = { -999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999.,-999. };
	Bool_t firstGuess = kTRUE;

	if (nGoodTracks == 4 && nSpdHits > 1) {
		//cout << "##### Kstar channel Started" << endl;
		for (Int_t pidLoopCounter = 0; pidLoopCounter < 1; pidLoopCounter++) { //loop over twice in case of 1 missing pion and 1 missing kaon, will double count the same event with both assumptions
			//Reset all variables for the next loop (in case of missing PID)
			newKaon = 0; newPion = 0; //added by Alec
			nHighPtTracks = 0;
			for (int aa = 0; aa < 13; aa++) {
				boostInfoA[aa] = -999.;
				boostInfoB[aa] = -999.;
			}
			//cout << "Before assumptions with " << nKaon << " Kaons and " << nPion << " Pions and " << nBoth << " M/B." << endl;
			//##### Deal with tracks missing TOF information. #####
			//Any events missing TOF info which can be completed to have the products for this decay channel will be assumed to belong to this channel.
			if ((nPion == 2 && nBoth == 2) || (nPion == 2 && nKaon == 1 && nBoth == 1)) {
				for (int i = 0; i < 4; i++) { //loop over the four tracks
					AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
					//if PID returns "both", identify the particle as a kaon
					if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.) 
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)  
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.))) 
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
						qKaon[nKaon] = trk->Charge();
						vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
						nKaon++;
						newKaon++;
					}
				}
			}
			if (((nKaon-newKaon) == 2 && nBoth == 2) || ((nKaon - newKaon) == 2 && nPion == 1 && nBoth == 1)) {
				for (int i = 0; i < 4; i++) { //loop over the four tracks
					AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
					//if PID returns "both", identify the particle as a pion
					if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
						qPion[nPion] = trk->Charge();
						vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
						nPion++;
						newPion++;
					}
				}
			}
			if ((nKaon - newKaon) == 1 && (nPion - newPion) == 1 && nBoth == 2) {
				//cout << "One Kaon, One Pion with two Missing tracks" << endl;
				if (firstGuess) {
					//cout << "First guess" << endl;
					for (int i = 0; i < 4; i++) { //loop over the four tracks
						if (nPion == 1) {
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							//identify the first "both" as a pion
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
								qPion[nPion] = trk->Charge();
								vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
								nPion++;
								newPion++;
							}
						}
						else {
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							//identify the second "both" as a kaon
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
								qKaon[nKaon] = trk->Charge();
								vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
								nKaon++;
								newKaon++;
							}
						}
						//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
					}
					pidLoopCounter--;
					firstGuess = kFALSE;
				}
				else {
					//cout << "Second Guess" << endl;
					for (int i = 0; i < 4; i++) { //loop over the four tracks
						if (nKaon == 1) {
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							//identify the first "both" as a kaon
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
								qKaon[nKaon] = trk->Charge();
								vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
								nKaon++;
								newKaon++;
							}
						}
						else {
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							//identify the second "both" as a pion
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
								qPion[nPion] = trk->Charge();
								vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
								nPion++;
								newPion++;
							}
						}
						//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
					}
				}
			}
			//cout << "Track assumptions complete with " << newKaon << " new Kaons and " << newPion << " new Pions." << endl;

			//##### Analyze events, fill the histos. #####
			//Select events with products matching this decay channel
			if ((nPion == 2) && (nKaon == 2)) {
				//cout << "##### Kstar channel Analyzing..." << endl;
				fHistNeventsEtaCKstarChannel->Fill(1);
				if (qKaon[0] * qKaon[1] > 0) fHistNeventsEtaCKstarChannel->Fill(2);
				if (qPion[0] * qPion[1] > 0) fHistNeventsEtaCKstarChannel->Fill(3);
				if ((qKaon[0] * qKaon[1] > 0) && (qPion[0] * qPion[1] > 0)) fHistNeventsEtaCKstarChannel->Fill(4); //events where both kaon and pion pairs have the same sign
				if ((qKaon[0] * qKaon[1] < 0) && (qPion[0] * qPion[1] < 0)) { //events with same sign, and total charge 0
					fHistNeventsEtaCKstarChannel->Fill(5);
					vCandidate = vPion[0] + vPion[1] + vKaon[0] + vKaon[1]; //EtaC candidate mass.
					//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

					//Belle Cuts - Some possible cuts used by the Belle Collaboration. Skip for now, histograms filled to inspect for cuts #####
					SumPz = vPion[0].Pz() + vPion[1].Pz() + vKaon[0].Pz() + vKaon[1].Pz();
					ScalarSumP = vPion[0].P() + vPion[1].P() + vKaon[0].P() + vKaon[1].P();
					sumPtVector = vPion[0].Vect() + vPion[1].Vect() + vKaon[0].Vect() + vKaon[1].Vect();
					VectorSumPt = sumPtVector.Pt();

					//Number of tracks with pT>0.4 GeV/c
					for (Int_t aa = 0; aa < 2; aa++) {
						if (vPion[aa].Pt() > 0.4) nHighPtTracks++;
						if (vKaon[aa].Pt() > 0.4) nHighPtTracks++;
					}
					if (nHighPtTracks > 1) fHistNeventsEtaCKstarChannel->Fill(9); //two or more tracks with higher momentum -- why is this a bin?

					//##### Define the Kstar(892) candidates, each event has two candidates for each K/Pi pair #####
					if ((qPion[0] * qKaon[0] < 0) && (qPion[1] * qKaon[1] < 0)) {
						vKstar[0] = vPion[0] + vKaon[0];
						vKstar[1] = vPion[1] + vKaon[1];
						BoostCut(vPion[0], vKaon[0], vKstar[0], boostInfoA);
						BoostCut(vPion[1], vKaon[1], vKstar[1], boostInfoB);
					}
					else if ((qPion[0] * qKaon[1] < 0) && (qPion[1] * qKaon[0] < 0)) {
						vKstar[0] = vPion[0] + vKaon[1];
						vKstar[1] = vPion[1] + vKaon[0];
						BoostCut(vPion[0], vKaon[1], vKstar[0], boostInfoA);
						BoostCut(vPion[1], vKaon[0], vKstar[1], boostInfoB);
					}

					//Boost/Helicity cut
					fKstarParentPx->Fill(boostInfoA[9]);
					fKstarParentPy->Fill(boostInfoA[10]);
					fKstarParentPz->Fill(boostInfoA[11]);
					fKstarDaughterParentAngle->Fill(boostInfoA[4]);
					fKstarDaughterParentAngle->Fill(boostInfoA[5]);
					fKstarDaughterDaughterAngle->Fill(boostInfoA[6]);
					fKstarDaughterDaughterCosAngle->Fill(boostInfoA[7]);
					fKstarDaughterPtotal->Fill(boostInfoA[8]);
					fKstarDaughterPtotalNorm->Fill(boostInfoA[12]);
					//*****Changed to <45. as a check.*****
					if (fabs(boostInfoA[4]) > 45. && fabs(boostInfoA[5]) > 45. && fabs(boostInfoA[6]) > 160. && fabs(boostInfoA[8]) < 0.2) goodPairA = kTRUE; //&& fabs(boostInfoA[7]) < 0.95
					else goodPairA = kFALSE;

					fKstarParentPx->Fill(boostInfoB[9]);
					fKstarParentPy->Fill(boostInfoB[10]);
					fKstarParentPz->Fill(boostInfoB[11]);
					fKstarDaughterParentAngle->Fill(boostInfoB[4]);
					fKstarDaughterParentAngle->Fill(boostInfoB[5]);
					fKstarDaughterDaughterAngle->Fill(boostInfoB[6]);
					fKstarDaughterDaughterCosAngle->Fill(boostInfoB[7]);
					fKstarDaughterPtotal->Fill(boostInfoB[8]);
					fKstarDaughterPtotalNorm->Fill(boostInfoB[12]);
					//*****Changed to < 45. as a check.*****
					if (fabs(boostInfoB[4]) > 45. && fabs(boostInfoB[5]) > 45. && fabs(boostInfoB[6]) > 160. && fabs(boostInfoB[8]) < 0.2) goodPairB = kTRUE; //&& fabs(boostInfoB[7]) < 0.95
					else goodPairB = kFALSE;

					//Now fill check histos with only pairs that pass
					if (goodPairA) {
						fKstarParentPxCheck->Fill(boostInfoA[9]);
						fKstarParentPyCheck->Fill(boostInfoA[10]);
						fKstarParentPzCheck->Fill(boostInfoA[11]);
						fKstarDaughterParentAngleCheck->Fill(boostInfoA[4]);
						fKstarDaughterParentAngleCheck->Fill(boostInfoA[5]);
						fKstarDaughterDaughterAngleCheck->Fill(boostInfoA[6]);
						fKstarDaughterDaughterCosAngleCheck->Fill(boostInfoA[7]);
						fKstarDaughterPtotalCheck->Fill(boostInfoA[8]);
						fKstarDaughterPtotalNormCheck->Fill(boostInfoA[12]);
					}
					if (goodPairB) {
						fKstarParentPxCheck->Fill(boostInfoB[9]);
						fKstarParentPyCheck->Fill(boostInfoB[10]);
						fKstarParentPzCheck->Fill(boostInfoB[11]);
						fKstarDaughterParentAngleCheck->Fill(boostInfoB[4]);
						fKstarDaughterParentAngleCheck->Fill(boostInfoB[5]);
						fKstarDaughterDaughterAngleCheck->Fill(boostInfoB[6]);
						fKstarDaughterDaughterCosAngleCheck->Fill(boostInfoB[7]);
						fKstarDaughterPtotalCheck->Fill(boostInfoB[8]);
						fKstarDaughterPtotalNormCheck->Fill(boostInfoB[12]);
					}

					//##### Turn off Helicity Cut #####
					goodPairA = kTRUE;    goodPairB = kTRUE;

					//Fill Dalitz plot with PiK masses (PiK) vs Pi Pi
					if ((qPion[0] * qKaon[0] < 0) && (qPion[1] * qKaon[1] < 0)) {
						fMPiKvsMPiK->Fill(vPion[0].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
						fMPiKvsMPiK->Fill(vPion[1].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
					}
					else if ((qPion[0] * qKaon[1] < 0) && (qPion[1] * qKaon[0] < 0)) {
						fMPiKvsMPiK->Fill(vPion[0].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
						fMPiKvsMPiK->Fill(vPion[1].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
					}

					//##### Fill histos - Determine number of intermediate k*(892)s. #####
					//2 Kstar case
					if ((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
						(vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth)) && 
						goodPairA && goodPairB) { //if both candidates are within the Kstar width of the Kstar mass and pass the Boost/Helicity cut
						fHistNeventsEtaCKstarChannel->Fill(6); //2 kstar events bin
						fHistNeventsEtaC->Fill(4);
						fEtaCCandidatesPerChannel->Fill(1);
						if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(1);
						
						if (qPion[0] > 0 && qPion[1] < 0) { //the and statement is not necessary, since we've already cut all non opposide charges since it has a Kstar mass
							f2KstarPtPiPlus->Fill(vPion[0].Pt());
							f2KstarPtPiMinus->Fill(vPion[1].Pt());
						}
						else {
							f2KstarPtPiPlus->Fill(vPion[1].Pt());
							f2KstarPtPiMinus->Fill(vPion[0].Pt());
						}
						if (qKaon[0] > 0 && qKaon[1] < 0) {
							f2KstarPtKPlus->Fill(vKaon[0].Pt());
							f2KstarPtKMinus->Fill(vKaon[1].Pt());
						}
						else {
							f2KstarPtKPlus->Fill(vKaon[1].Pt());
							f2KstarPtKMinus->Fill(vKaon[0].Pt());
						}

						if ((qPion[0] * qKaon[0] < 0) && (qPion[1] * qKaon[1] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[0].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
							fKstarMPiKvsMPiK->Fill(vPion[1].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
						}
						else if ((qPion[0] * qKaon[1] < 0) && (qPion[1] * qKaon[0] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[0].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
							fKstarMPiKvsMPiK->Fill(vPion[1].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
						}

						//Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
						f2KstarTPCsignalPion->Fill(fRecTPCsignalPion[0], fRecTPCsignalPion[1]);
						f2KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0], fRecTPCsignalKaon[1]);
						f2KstarDedxVsPtPion->Fill(vPion[0].Pt(), fRecTPCsignalPion[0]);		f2KstarDedxVsPtPion->Fill(vPion[1].Pt(), fRecTPCsignalPion[1]); //vs P???
						f2KstarDedxVsPtKaon->Fill(vKaon[0].Pt(), fRecTPCsignalKaon[0]);		f2KstarDedxVsPtKaon->Fill(vKaon[1].Pt(), fRecTPCsignalKaon[1]);
						
						//Fill intermediate Kstar histos, one for each candidate
						f2KstarPtVsMinvFirstKstar->Fill(vKstar[0].M(), vKstar[0].Pt());
						f2KstarPtVsMinvSecondKstar->Fill(vKstar[1].M(), vKstar[1].Pt());

						//Fill EtaC histos
						f2KstarPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						if (vCandidate.Pt() < 0.11) {
							fChannelVsMinvEtaC->Fill(vCandidate.M(), 1);
							fAllMinvEtaCLowPt->Fill(vCandidate.M());
							if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
						}
						f2KstarEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.4) f2KstarEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.1) f2KstarEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						//Belle cut histos
						f2KstarSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
						f2KstarScalarSumP->Fill(ScalarSumP);
						f2KstarVectorSumPt->Fill(VectorSumPt);
					}

					//1 Kstar case
					else if ((vKstar[0].M() < (kStarMass + kStarWidth)) && (vKstar[0].M() > (kStarMass - kStarWidth)) &&
						((vKstar[1].M() > (kStarMass + kStarWidth)) || (vKstar[1].M() < (kStarMass - kStarWidth))) &&
						goodPairA) {
						//Fill using first Kstar candidate
						fHistNeventsEtaCKstarChannel->Fill(7);
						fHistNeventsEtaC->Fill(4);
						fEtaCCandidatesPerChannel->Fill(2);
						if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(2);
						if (qPion[0] > 0 && qPion[1] < 0) {
							f1KstarPtPiPlus->Fill(vPion[0].Pt());
							f1KstarPtPiMinus->Fill(vPion[1].Pt());
						}
						else {
							f1KstarPtPiPlus->Fill(vPion[1].Pt());
							f1KstarPtPiMinus->Fill(vPion[0].Pt());
						}
						if (qKaon[0] > 0 && qKaon[1] < 0) {
							f1KstarPtKPlus->Fill(vKaon[0].Pt());
							f1KstarPtKMinus->Fill(vKaon[1].Pt());
						}
						else {
							f1KstarPtKPlus->Fill(vKaon[1].Pt());
							f1KstarPtKMinus->Fill(vKaon[0].Pt());
						}

						if ((qPion[0] * qKaon[0] < 0) && (qPion[1] * qKaon[1] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[0].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
						}
						else if ((qPion[0] * qKaon[1] < 0) && (qPion[1] * qKaon[0] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[0].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
						}

						//Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
						f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0], fRecTPCsignalPion[1]);
						f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0], fRecTPCsignalKaon[1]);
						f1KstarDedxVsPtPion->Fill(vPion[0].Pt(), fRecTPCsignalPion[0]);		f1KstarDedxVsPtPion->Fill(vPion[1].Pt(), fRecTPCsignalPion[1]);
						f1KstarDedxVsPtKaon->Fill(vKaon[0].Pt(), fRecTPCsignalKaon[0]);		f1KstarDedxVsPtKaon->Fill(vKaon[1].Pt(), fRecTPCsignalKaon[1]);

						//Fill intermediate Kstar histos, one for each candidate
						f1KstarPtVsMinvKstar->Fill(vKstar[0].M(), vKstar[0].Pt());
						f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[1].M(), vKstar[1].Pt());

						//Fill EtaC histos
						f1KstarPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						if (vCandidate.Pt() < 0.11) {
							fChannelVsMinvEtaC->Fill(vCandidate.M(), 2);
							fAllMinvEtaCLowPt->Fill(vCandidate.M());
							if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
						}
						f1KstarEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.4) f1KstarEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.1) f1KstarEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						//Belle cut histos
						f1KstarSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
						f1KstarScalarSumP->Fill(ScalarSumP);
						f1KstarVectorSumPt->Fill(VectorSumPt);
					}
					else if ((vKstar[1].M() < (kStarMass + kStarWidth)) && (vKstar[1].M() > (kStarMass - kStarWidth)) &&
						((vKstar[0].M() > (kStarMass + kStarWidth)) || (vKstar[0].M() < (kStarMass - kStarWidth))) &&
						goodPairB) {
						//Fill using second Kstar candidate
						fHistNeventsEtaCKstarChannel->Fill(7);
						fHistNeventsEtaC->Fill(4);
						fEtaCCandidatesPerChannel->Fill(2);
						if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(2);

						if (qPion[0] > 0 && qPion[1] < 0) {
							f1KstarPtPiPlus->Fill(vPion[0].Pt());
							f1KstarPtPiMinus->Fill(vPion[1].Pt());
						}
						else {
							f1KstarPtPiPlus->Fill(vPion[1].Pt());
							f1KstarPtPiMinus->Fill(vPion[0].Pt());
						}
						if (qKaon[0] > 0 && qKaon[1] < 0) {
							f1KstarPtKPlus->Fill(vKaon[0].Pt());
							f1KstarPtKMinus->Fill(vKaon[1].Pt());
						}
						else {
							f1KstarPtKPlus->Fill(vKaon[1].Pt());
							f1KstarPtKMinus->Fill(vKaon[0].Pt());
						}

						if ((qPion[0] * qKaon[0] < 0) && (qPion[1] * qKaon[1] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[1].M()*vKaon[1].M(), vPion[0].M()*vPion[1].M());
						}
						else if ((qPion[0] * qKaon[1] < 0) && (qPion[1] * qKaon[0] < 0)) {
							fKstarMPiKvsMPiK->Fill(vPion[1].M()*vKaon[0].M(), vPion[0].M()*vPion[1].M());
						}

						//Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
						f1KstarTPCsignalPion->Fill(fRecTPCsignalPion[0], fRecTPCsignalPion[1]);
						f1KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0], fRecTPCsignalKaon[1]);
						f1KstarDedxVsPtPion->Fill(vPion[0].Pt(), fRecTPCsignalPion[0]);		f1KstarDedxVsPtPion->Fill(vPion[1].Pt(), fRecTPCsignalPion[1]);
						f1KstarDedxVsPtKaon->Fill(vKaon[0].Pt(), fRecTPCsignalKaon[0]);		f1KstarDedxVsPtKaon->Fill(vKaon[1].Pt(), fRecTPCsignalKaon[1]);

						//Fill intermediate Kstar histos, one for each candidate
						f1KstarPtVsMinvKstar->Fill(vKstar[1].M(), vKstar[1].Pt());
						f1KstarPtVsMinvOtherPiKcombo->Fill(vKstar[0].M(), vKstar[0].Pt());

						//Fill EtaC histos
						f1KstarPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						if (vCandidate.Pt() < 0.11) {
							fChannelVsMinvEtaC->Fill(vCandidate.M(), 2);
							fAllMinvEtaCLowPt->Fill(vCandidate.M());
							if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
						}
						f1KstarEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.4) f1KstarEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.1) f1KstarEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						//Belle cut histos
						f1KstarSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
						f1KstarScalarSumP->Fill(ScalarSumP);
						f1KstarVectorSumPt->Fill(VectorSumPt);
					}

					//0 Kstar case
					else {
						fHistNeventsEtaCKstarChannel->Fill(8);
						fHistNeventsEtaC->Fill(4);
						fEtaCCandidatesPerChannel->Fill(3);
						if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(3);

						if (qPion[0] > 0 && qPion[1] < 0) {
							f0KstarPtPiPlus->Fill(vPion[0].Pt());
							f0KstarPtPiMinus->Fill(vPion[1].Pt());
						}
						else {
							f0KstarPtPiPlus->Fill(vPion[1].Pt());
							f0KstarPtPiMinus->Fill(vPion[0].Pt());
						}
						if (qKaon[0] > 0 && qKaon[1] < 0) {
							f0KstarPtKPlus->Fill(vKaon[0].Pt());
							f0KstarPtKMinus->Fill(vKaon[1].Pt());
						}
						else {
							f0KstarPtKPlus->Fill(vKaon[1].Pt());
							f0KstarPtKMinus->Fill(vKaon[0].Pt());
						}

						//Fill TPC signal correlation for Pions-Pion correlation and Kaon-Kaon correlation
						f0KstarTPCsignalPion->Fill(fRecTPCsignalPion[0], fRecTPCsignalPion[1]);
						f0KstarTPCsignalKaon->Fill(fRecTPCsignalKaon[0], fRecTPCsignalKaon[1]);
						f0KstarDedxVsPtPion->Fill(vPion[0].Pt(), fRecTPCsignalPion[0]);		f0KstarDedxVsPtPion->Fill(vPion[1].Pt(), fRecTPCsignalPion[1]);
						f0KstarDedxVsPtKaon->Fill(vKaon[0].Pt(), fRecTPCsignalKaon[0]);		f0KstarDedxVsPtKaon->Fill(vKaon[1].Pt(), fRecTPCsignalKaon[1]);
						
						//Fill intermediate Kstar histos, one for each candidate
						f0KstarPtVsMinvFirstPiKcombo->Fill(vKstar[0].M(), vKstar[0].Pt());
						f0KstarPtVsMinvSecondPiKcombo->Fill(vKstar[1].M(), vKstar[1].Pt());
						
						//Fill EtaC histos
						f0KstarPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
						if (vCandidate.Pt() < 0.11) {
							fChannelVsMinvEtaC->Fill(vCandidate.M(), 3);
							fAllMinvEtaCLowPt->Fill(vCandidate.M());
							if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
						}
						f0KstarEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.4) f0KstarEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						if (vCandidate.Pt() < 0.1) f0KstarEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
						//Belle cut histos
						f0KstarSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
						f0KstarScalarSumP->Fill(ScalarSumP);
						f0KstarVectorSumPt->Fill(VectorSumPt);
					}
				}
			}
			//Remove the partiles assumed to match because of missing TOF for the next analysis
			nKaon -= newKaon;
			nPion -= newPion;

			//cout << "Kstar channel Completed"<< endl;
		}//end pidLoopCounter
	}
	//End Kstar Channel

	//EtaC -> K+ K- 2(pi+ pi-) Channel
	firstGuess = kTRUE;
	if (nGoodTracks == 6 && nSpdHits > 1) {
		//cout << "##### 2K4Pi channel Started" << endl;
		for (Int_t pidLoopCounter = 0; pidLoopCounter < 1; pidLoopCounter++) { //loop over twice in case of 1 missing pion and 1 missing kaon, will double count the same event with both assumptions
			//Reset all variables for the next loop (in case of missing PID)
			newKaon = 0; newPion = 0; //added by Alec
			nHighPtTracks = 0;
			//cout << "Before assumptions with " << nKaon << " Kaons and " << nPion << " Pions and " << nBoth << " M/B." << endl;
			//##### Deal with tracks missing TOF information. #####
			//Any events missing TOF info which can be completed to have the products for this decay channel will be assumed to belong to this channel.
			if ((nPion == 4 && nBoth == 2) || (nPion == 4 && nKaon == 1 && nBoth == 1)) {
				for (int i = 0; i < 6; i++) { //loop over the four tracks
					AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
					if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
						qKaon[nKaon] = trk->Charge();
						vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
						nKaon++;
						newKaon++;
					}
				}
			}
			if (((nKaon - newKaon) == 2 && nPion == 1 && nBoth == 3) || ((nKaon - newKaon) == 2 && nPion == 2 && nBoth == 2) || ((nKaon - newKaon) == 2 && nPion == 3 && nBoth == 1)) {
				for (int i = 0; i < 6; i++) { //loop over the four tracks
					AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
					if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
						qPion[nPion] = trk->Charge();
						vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
						nPion++;
						newPion++;
					}
				}
			}
			if ((nKaon - newKaon) == 1 && (nPion - newPion) == 3 && nBoth == 2) {
				//cout << "One Kaon, One Pion with two Missing tracks" << endl;
				if (firstGuess) {
					//cout << "First guess" << endl;
					for (int i = 0; i < 6; i++) { //loop over the four tracks
						if (nPion == 3) {
							AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
								qPion[nPion] = trk->Charge();
								vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
								nPion++;
								newPion++; //added by Alec
							}
						}
						else {
							AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
								qKaon[nKaon] = trk->Charge();
								vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
								nKaon++;
								newKaon++; //added by Alec
							}
						}
						//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
					}
					pidLoopCounter--;
					firstGuess = kFALSE;
				}
				else {
					//cout << "Second Guess" << endl;
					for (int i = 0; i < 6; i++) { //loop over the four tracks
						if (nKaon == 1) {
							AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
								qKaon[nKaon] = trk->Charge();
								vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
								nKaon++;
								newKaon++; //added by Alec
							}
						}
						else {
							AliAODTrack* trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
								qPion[nPion] = trk->Charge();
								vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
								nPion++;
								newPion++; //added by Alec
							}
						}
						//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
					}
				}
			}
			//cout << "Track assumptions complete with " << newKaon << " new Kaons and " << newPion << " new Pions." << endl;
			
			//##### Analyze events, fill the histos. #####
			if ((nPion == 4) && (nKaon == 2)) {
				//cout << "##### 2K4Pi channel Analyzing..." << endl;
				fHistNeventsEtaC2K4PiChannel->Fill(1);
				if ((qPion[0] + qPion[1] + qPion[2] + qPion[3] + qKaon[0] + qKaon[1]) != 0) fHistNeventsEtaC2K4PiChannel->Fill(2); //events with wrong total charge
				if ((qPion[0] + qPion[1] + qPion[2] + qPion[3] + qKaon[0] + qKaon[1]) == 0) { //events with same sign, and total charge 0
					fHistNeventsEtaC2K4PiChannel->Fill(3);
					fHistNeventsEtaC->Fill(4);
					fEtaCCandidatesPerChannel->Fill(9);

					vCandidate = vPion[0] + vPion[1] + vPion[2] + vPion[3] + vKaon[0] + vKaon[1]; //EtaC candidate mass.
					//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;
					if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(9);

					//Belle Cuts - Some possible cuts used by the Belle Collaboration. Skip for now, histograms filled to inspect for cuts #####
					SumPz = vPion[0].Pz() + vPion[1].Pz() + vPion[2].Pz() + vPion[3].Pz() + vKaon[0].Pz() + vKaon[1].Pz();
					ScalarSumP = vPion[0].P() + vPion[1].P() + vPion[2].P() + vPion[3].P() + vKaon[0].P() + vKaon[1].P();
					sumPtVector = vPion[0].Vect() + vPion[1].Vect() + vPion[2].Vect() + vPion[3].Vect() + vKaon[0].Vect() + vKaon[1].Vect();
					VectorSumPt = sumPtVector.Pt();

					//Number of tracks with pT>0.4 GeV/c
					for (Int_t aa = 0; aa < 2; aa++) {
						if (vKaon[aa].Pt() > 0.4) nHighPtTracks++;
					}
					for (Int_t aa = 0; aa < 4; aa++) {
						if (vPion[aa].Pt() > 0.4) nHighPtTracks++;
					}
					if (nHighPtTracks > 1) fHistNeventsEtaC2K4PiChannel->Fill(4); //two or more tracks with higher momentum -- why is this a bin?

					//Fill EtaC histos
					f2K4PiPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					if (vCandidate.Pt() < 0.11) {
						fChannelVsMinvEtaC->Fill(vCandidate.M(), 9);
						fAllMinvEtaCLowPt->Fill(vCandidate.M());
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
					}
					f2K4PiEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.4) f2K4PiEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					if (vCandidate.Pt() < 0.1) f2K4PiEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
					//Belle cut histos
					f2K4PiSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
					f2K4PiScalarSumP->Fill(ScalarSumP);
					f2K4PiVectorSumPt->Fill(VectorSumPt);
				}
			}
			
			//Remove the partiles assumed to match because of missing TOF for the next analysis
			nKaon -= newKaon;
			nPion -= newPion;
			//cout << "##### 2K4Pi channel Complete." << endl;
		}//end pidLoopCounter
	}
	//end K+ K- 2(Pi+ Pi-) channel

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
				if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
					&& (beta > 0.55 || beta == -999.)) {
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
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3]) != 0) fHistNeventsEtaCRhoChannel->Fill(2); //non-zero net charge
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3]) == 0) {
				fHistNeventsEtaCRhoChannel->Fill(3); //zero net charge
				vCandidate = vPion[0] + vPion[1] + vPion[2] + vPion[3];
				//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

				SumPz = vPion[0].Pz() + vPion[1].Pz() + vPion[2].Pz() + vPion[3].Pz();
				ScalarSumP = vPion[0].P() + vPion[1].P() + vPion[2].P() + vPion[3].P();
				sumPtVector = vPion[0].Vect() + vPion[1].Vect() + vPion[2].Vect() + vPion[3].Vect();
				VectorSumPt = sumPtVector.Pt();

				//Number of tracks with pT>0.4 GeV/c
				for (Int_t aa = 0; aa < 4; aa++) if (vPion[aa].Pt() > 0.4) nHighPtTracks++;

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
					fHistNeventsEtaC->Fill(4);
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
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
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
					fHistNeventsEtaC->Fill(4);
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
					f2Rho1PairPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt()); //4Pi final states with 1 intermediate rho
					fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
					if (vCandidate.Pt() < 0.11) {
						fChannelVsMinvEtaC->Fill(vCandidate.M(), 4);
						fAllMinvEtaCLowPt->Fill(vCandidate.M());
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
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
					fHistNeventsEtaC->Fill(4);
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
						if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
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



  //EtaC->3(pi+pi-) Channel
	newPion = 0;
	nHighPtTracks = 0;
	
	if (nGoodTracks == 6 && nSpdHits > 1) {
		//cout << "##### 3PiPi channel Started." << endl;
		//##### missing TOF
		if ((nPion == 3 && nBoth == 3) || (nPion == 4 && nBoth == 2) || (nPion == 5 && nBoth == 1)) {
			for (int i = 0; i < 6; i++) { //loop over the four tracks
				AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
				if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
					&& (beta > 0.55 || beta == -999.)) {
					fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
					qPion[nPion] = trk->Charge();
					vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
					nPion++;
					newPion++;
				}
			}
		}

		//##### Analyze good events, fill the histos.
		if ((nPion == 6)) {
			//cout << "##### 3PiPi Channel Analyzing..." << endl;
			fHistNeventsEtaC3PiPiChannel->Fill(1); //6 pions
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3] + qPion[4] + qPion[5]) != 0) fHistNeventsEtaC3PiPiChannel->Fill(2); //non-zero net charge
			if ((qPion[0] + qPion[1] + qPion[2] + qPion[3] + qPion[4] + qPion[5]) == 0) {
				fHistNeventsEtaC3PiPiChannel->Fill(3); //zero net charge, candidate
				fHistNeventsEtaC->Fill(4);
				fEtaCCandidatesPerChannel->Fill(6);

				vCandidate = vPion[0] + vPion[1] + vPion[2] + vPion[3] + vPion[4] + vPion[5];
				
				if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(6);

				SumPz = vPion[0].Pz() + vPion[1].Pz() + vPion[2].Pz() + vPion[3].Pz() + vPion[4].Pz() + vPion[5].Pz();
				ScalarSumP = vPion[0].P() + vPion[1].P() + vPion[2].P() + vPion[3].P() + vPion[4].P() + vPion[5].P();
				sumPtVector = vPion[0].Vect() + vPion[1].Vect() + vPion[2].Vect() + vPion[3].Vect() + vPion[4].Vect() + vPion[5].Vect();
				VectorSumPt = sumPtVector.Pt();
				//Number of tracks with pT>0.4 GeV/c
				for (Int_t aa = 0; aa < 6; aa++) if (vPion[aa].Pt() > 0.4) nHighPtTracks++;

				f3PiPiPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
				fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
				if (vCandidate.Pt() < 0.11) {
					fChannelVsMinvEtaC->Fill(vCandidate.M(), 6);
					fAllMinvEtaCLowPt->Fill(vCandidate.M());
					if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
				}
				f3PiPiEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
				if (vCandidate.Pt() < 0.4) f3PiPiEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
				if (vCandidate.Pt() < 0.1) f3PiPiEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
				f3PiPiSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
				f3PiPiScalarSumP->Fill(ScalarSumP);
				f3PiPiVectorSumPt->Fill(VectorSumPt);
				if (nHighPtTracks > 1) fHistNeventsEtaC3PiPiChannel->Fill(4);
				//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

				//Fill Dalitz plots?
			}
		}
		//cout << "##### 3PiPi channel Complete." << endl;
		//Remove the partiles assumed to match because of missing TOF for the next analysis
		nPion -= newPion;
	}
	//End EtaC->3(pi+pi-) Channel

  //EtaC->2(k+k-) Channel
	newKaon = 0;
	nHighPtTracks = 0;
	Int_t nKMinus = 0;
	Int_t nKPlus = 0;
	TLorentzVector vKK[4], vKaonMinus[4], vKaonPlus[4];

	if (nGoodTracks == 4 && nSpdHits > 1) {
		//cout << "##### 2KK channel Started." << endl;
		//##### Missing TOF
		if ((nKaon == 2 && nBoth == 2) || (nKaon == 3 && nBoth == 1)) {
			for (int i = 0; i < 4; i++) { //loop over the four tracks
				AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
				if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
					|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
					&& (beta > 0.55 || beta == -999.)) {
					fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
					qKaon[nKaon] = trk->Charge();
					vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
					nKaon++;
					newKaon++;
				}
			}
		}
		if (nPion + nNeither == 0) f4KnKaonVsnewKaon->Fill(newKaon,nKaon-newKaon);

		//Analyze good events, fill the histos.
		if ((nKaon == 4)) {
			//cout << "##### 2KK channel Analyzing..." << endl;
			fHistNeventsEtaC4KaonChannel->Fill(1); //4 kaons
			if ((qKaon[0] + qKaon[1] + qKaon[2] + qKaon[3]) != 0) fHistNeventsEtaC4KaonChannel->Fill(2); //non-zero net charge
			if ((qKaon[0] + qKaon[1] + qKaon[2] + qKaon[3]) == 0) {
				fHistNeventsEtaC4KaonChannel->Fill(3); //zero net charge
				fHistNeventsEtaC->Fill(4);
				fEtaCCandidatesPerChannel->Fill(7);

				vCandidate = vKaon[0] + vKaon[1] + vKaon[2] + vKaon[3];
				if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(7);

				SumPz = vKaon[0].Pz() + vKaon[1].Pz() + vKaon[2].Pz() + vKaon[3].Pz();
				ScalarSumP = vKaon[0].P() + vKaon[1].P() + vKaon[2].P() + vKaon[3].P();
				sumPtVector = vKaon[0].Vect() + vKaon[1].Vect() + vKaon[2].Vect() + vKaon[3].Vect();
				VectorSumPt = sumPtVector.Pt();
				//Number of tracks with pT>0.4 GeV/c
				for (Int_t aa = 0; aa < 4; aa++) if (vKaon[aa].Pt() > 0.4) nHighPtTracks++;

				f4KaonPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
				fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
				if (vCandidate.Pt() < 0.11) {
					fChannelVsMinvEtaC->Fill(vCandidate.M(), 7);
					fAllMinvEtaCLowPt->Fill(vCandidate.M());
					if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
				}
				f4KaonEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
				if (vCandidate.Pt() < 0.4) f4KaonEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
				if (vCandidate.Pt() < 0.1) f4KaonEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
				f4KaonSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
				f4KaonScalarSumP->Fill(ScalarSumP);
				f4KaonVectorSumPt->Fill(VectorSumPt);
				if (nHighPtTracks > 1) fHistNeventsEtaC4KaonChannel->Fill(4);
				//cout << "mEtaC " << vCandidate.M() << ", ptEtaC " << vCandidate.Pt() << endl;

				//Now the intermediate K+K- pairs (all combinations)
				for (Int_t numKK = 0; numKK < 4; numKK++) vKK[numKK].SetPtEtaPhiM(0., 0., 0., 0.);

				//Get masses of potential intermediate Rho's
				if (qKaon[0] < 0) {
					vKaonMinus[nKMinus] = vKaon[0];
					nKMinus++;
				}
				else {
					vKaonPlus[nKPlus] = vKaon[0];
					nKPlus++;
				}
				if (qKaon[1] < 0) {
					vKaonMinus[nKMinus] = vKaon[1];
					nKMinus++;
				}
				else {
					vKaonPlus[nKPlus] = vKaon[1];
					nKPlus++;
				}
				if (qKaon[2] < 0) {
					vKaonMinus[nKMinus] = vKaon[2];
					nKMinus++;
				}
				else {
					vKaonPlus[nKPlus] = vKaon[2];
					nKPlus++;
				}
				if (qKaon[3] < 0) {
					vKaonMinus[nKMinus] = vKaon[3];
					nKMinus++;
				}
				else {
					vKaonPlus[nKPlus] = vKaon[3];
					nKPlus++;
				}
				//Either 0 and 1 are rho's or 2 and 3. If both sets are rho's choose best set.
				vKK[0] = vKaonMinus[0] + vKaonPlus[0];
				vKK[1] = vKaonMinus[1] + vKaonPlus[1];
				vKK[2] = vKaonMinus[1] + vKaonPlus[0];
				vKK[3] = vKaonMinus[0] + vKaonPlus[1];

				for (Int_t numKK = 0; numKK < 4; numKK++) {
					f4KaonPtVsMinvKK->Fill(vKK[numKK].M(), vKK[numKK].Pt());
				}

				//Dalitz plots -- For all cases look at 4K vs 2K minv
				for (Int_t KKIndex = 0; KKIndex < 4; KKIndex++) {
					f4KVs2KMinv->Fill(vCandidate.M(), vKK[KKIndex].M());
					f4KVs2KMinvSquared->Fill((vCandidate.M()*vCandidate.M()), (vRho[KKIndex].M()*vKK[KKIndex].M()));
				}
				fM2KKVsM2KK->Fill((vKK[0].M()*vKK[0].M()), (vKK[1].M()*vKK[1].M()));
				fM2KKVsM2KK->Fill((vKK[2].M()*vKK[2].M()), (vKK[3].M()*vKK[3].M()));
			}
		}
		//Remove the partiles assumed to match because of missing TOF for the next analysis
		//cout << "##### 2KK channel Completed." << endl;
		nKaon -= newKaon;
	}
	//End EtaC->2(k+k-) Channel

	if (nGoodTracks == 4 && nSpdHits > 1) {

		fHistPostFourTracksNkaon->Fill(nKaon);
		fHistPostFourTracksNpion->Fill(nPion);
		fListNPionChange->Fill(nPion - nPrePion);
		fListNKaonChange->Fill(nKaon - nPreKaon);
	}

	
    //EtaC->K0short case (using V0s)
	//make it harder to get into here
	if ((nGoodTracks == 4 && ((nPion == 3 && (nBoth + nKaon) == 1) || (nPion == 2 && nBoth == 2) || (nPion == 2 && nKaon == 1 && nBoth == 1))) 
		|| (nGoodTracks == 5 && ((nPion == 4 && (nBoth + nKaon) == 1) || (nPion == 3 && nBoth == 1 && nKaon == 1) || (nPion == 3 && nBoth == 2) || (nPion == 2 && nBoth == 3) || (nPion == 2 && nKaon == 1 && nBoth == 2)))
		|| (nGoodTracks == 6 && ((nPion == 5 && (nBoth + nKaon) == 1) || (nPion == 4 && nBoth == 1 && nKaon == 1) || (nPion == 4 && nBoth == 2) || (nPion == 3 && nBoth == 3) || (nPion == 3 && nKaon == 1 && nBoth == 2) || (nPion == 2 && nKaon == 1 && nBoth == 3) || (nPion == 2 && nBoth == 4)))
		&& nSpdHits > 1) {
		//K0s cut definitions
		Bool_t fCutCheck = kFALSE; //change if cut checks are implimented.

		const float kMinDCAPrimaryPion = 0.4; //Loosen or remove
		const float kMaxDCADaughtersK0 = 1.0; //originally was 0.3;
		const float kMaxDCAK0 = 0.3;
		const float kMaxDLK0 = 30.0;
		const float kMinDLK0 = 0.2; //fMinDecayLength;
		const float kEtaCut = 0.8;
		const float kMinCosAngle = 0.98; //originally was0.99;

										 //for cut checks
		double kCheckMassLow[3] = { 0.49, 0.48, 0.45 };
		double kCheckMassHigh[3] = { 0.505, 0.515, 0.550 };
		double kCheckDCAK0[3] = { 0.1, 0.3, 1.0 };
		double kCheckDCAPi[3] = { 1.0, 0.4, 0.1 };
		double kCheckDCAPiPi[3] = { 0.1, 0.3, 1.0 };

		int bestK0sCandidateIndex = 0;
		int k0ShortCount = 0;
		int UsedNMissingDaughters = 0;
		ULong_t statusPos = 0;
		ULong_t statusNeg = 0;
		int nMissingDaughters = 0;

		//First identify any K0 candidates, keep track of the best candidate, and mark its daughter pions
		for (int i = 0; i < aod->GetNumberOfV0s(); i++) {
			bool goodPiPlus = kFALSE;
			bool goodPiMinus = kFALSE;
			short int pos0or1;
			short int neg0or1;

			//    cout << " Num V0s " << aod->GetNumberOfV0s() << endl;
			if (aod->GetNumberOfV0s() > 10) break;

			//load v0 track
			AliAODv0* v0 = aod->GetV0(i);
			if (!v0) continue;
			Bool_t fOnlineCase = kFALSE; //Change this code if implementing online case.
			if (fOnlineCase) { //need to define this somewhere (in input to analysis perhaps).
				if (!(v0->GetOnFlyStatus())) continue;
			} //for online
			else { //offline case, currently used
				if (v0->GetOnFlyStatus()) continue; //GetOnFlyStatus is true if this v0 is reconstructed
			}

			//for on-the-fly ordering
			AliAODTrack* tempTrack = (AliAODTrack*)v0->GetDaughter(0);
			if (tempTrack->Charge() > 0) { pos0or1 = 0; neg0or1 = 1; }
			else { pos0or1 = 1; neg0or1 = 0; }

			//load daughter tracks
			AliAODTrack* prongTrackPos = (AliAODTrack*)v0->GetDaughter(pos0or1);
			AliAODTrack* prongTrackNeg = (AliAODTrack*)v0->GetDaughter(neg0or1);
			if (!prongTrackPos) continue;
			if (!prongTrackNeg) continue;

			//daughter cuts
			if (v0->PtProng(pos0or1) < 0.15) continue;
			if (v0->PtProng(neg0or1) < 0.15) continue;
			if (fabs(v0->EtaProng(pos0or1)) > 0.8) continue;
			if (fabs(v0->EtaProng(neg0or1)) > 0.8) continue;
			//both daughters must be good enough
			if (!(prongTrackPos->TestFilterBit(1 << 0)) || !(prongTrackNeg->TestFilterBit(1 << 0))) continue;
			if (prongTrackPos->GetTPCNcls() < 50 || prongTrackNeg->GetTPCNcls() < 50) continue;
			if (prongTrackPos->Chi2perNDF() > 4 || prongTrackNeg->Chi2perNDF() > 4) continue;

			//load status for PID
			statusPos = prongTrackPos->GetStatus();
			prongTrackPos->SetAODEvent(aod);
			statusNeg = prongTrackNeg->GetStatus();
			prongTrackNeg->SetAODEvent(aod);

			//cout << "#################### Before PID block 2" << endl;

			//Get nsigma info for PID
			fPIDTPCMuonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos, AliPID::kMuon);
			fPIDTPCElectronPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos, AliPID::kElectron);
			fPIDTPCPionPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos, AliPID::kPion);
			fPIDTPCKaonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos, AliPID::kKaon);
			fPIDTPCProtonPos[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackPos, AliPID::kProton);

			if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, prongTrackPos) == AliPIDResponse::kDetPidOk) { //3 = kTOF
				fPIDTOFMuonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos, AliPID::kMuon);
				fPIDTOFElectronPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos, AliPID::kElectron);
				fPIDTOFPionPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos, AliPID::kPion);
				fPIDTOFKaonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos, AliPID::kKaon);
				fPIDTOFProtonPos[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackPos, AliPID::kProton);
			}
			else {
				fPIDTOFMuonPos[i] = -999.;
				fPIDTOFElectronPos[i] = -999.;
				fPIDTOFPionPos[i] = -999.;
				fPIDTOFKaonPos[i] = -999.;
				fPIDTOFProtonPos[i] = -999.;
			}

			fPIDTPCMuonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg, AliPID::kMuon);
			fPIDTPCElectronNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg, AliPID::kElectron);
			fPIDTPCPionNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg, AliPID::kPion);
			fPIDTPCKaonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg, AliPID::kKaon);
			fPIDTPCProtonNeg[i] = fPIDResponse->NumberOfSigmasTPC(prongTrackNeg, AliPID::kProton);

			if (fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, prongTrackNeg) == AliPIDResponse::kDetPidOk) { //3 = kTOF
				fPIDTOFMuonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg, AliPID::kMuon);
				fPIDTOFElectronNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg, AliPID::kElectron);
				fPIDTOFPionNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg, AliPID::kPion);
				fPIDTOFKaonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg, AliPID::kKaon);
				fPIDTOFProtonNeg[i] = fPIDResponse->NumberOfSigmasTOF(prongTrackNeg, AliPID::kProton);
			}
			else {
				fPIDTOFMuonNeg[i] = -999.;
				fPIDTOFElectronNeg[i] = -999.;
				fPIDTOFPionNeg[i] = -999.;
				fPIDTOFKaonNeg[i] = -999.;
				fPIDTOFProtonNeg[i] = -999.;
			}

			//    cout << "#################### After PID block 2" << endl;

			eventStartTime = aod->GetTOFHeader()->GetDefaultEventTimeVal(); //may have to define this differently. An array of start times is returned
			if (prongTrackPos->GetIntegratedLength() >= 360. && prongTrackPos->GetIntegratedLength() <= 800. && prongTrackPos->GetTOFsignal() > 0. && eventStartTime < 999990.0) {
				beta = (prongTrackPos->GetIntegratedLength()*0.01) / ((prongTrackPos->GetTOFsignal() - eventStartTime)*(1.e-12)*TMath::C()); //0.01 converts cm to m, 1e-12 converts ps to sec

				fK0sTOFbetaVsPAll->Fill(prongTrackPos->P(), beta);
				if (fabs(fPIDTPCPionPos[i]) < 3. || fabs(fPIDTOFPionPos[i]) < 3.) {//if TPC or TOF returns pion
					fK0sTOFbetaVsPPion->Fill(prongTrackPos->P(), beta);
				}
			}
			if (prongTrackNeg->GetIntegratedLength() >= 360. && prongTrackNeg->GetIntegratedLength() <= 800. && prongTrackNeg->GetTOFsignal() > 0. && eventStartTime < 999990.0) {
				beta = (prongTrackNeg->GetIntegratedLength()*0.01) / ((prongTrackNeg->GetTOFsignal() - eventStartTime)*(1.e-12)*TMath::C()); //0.01 converts cm to m, 1e-12 converts ps to sec

				fK0sTOFbetaVsPAll->Fill(prongTrackNeg->P(), beta);
				if (fabs(fPIDTPCPionNeg[i]) < 3. || fabs(fPIDTOFPionNeg[i]) < 3.) {
					fK0sTOFbetaVsPPion->Fill(prongTrackNeg->P(), beta);
				}
			}
			fK0sDedxVsPAll->Fill(prongTrackPos->P(), prongTrackPos->GetTPCsignal());
			fK0sDedxVsPAll->Fill(prongTrackNeg->P(), prongTrackNeg->GetTPCsignal());
			if (fabs(fPIDTPCPionPos[i]) < 3. || fabs(fPIDTOFPionPos[i]) < 3.) {//if TPC or TOF returns pion
				goodPiPlus = kTRUE;
				fK0sDedxVsPPion->Fill(prongTrackPos->P(), prongTrackPos->GetTPCsignal());
			}
			if (fabs(fPIDTPCPionNeg[i]) < 3. || fabs(fPIDTOFPionNeg[i]) < 3.) {
				goodPiMinus = kTRUE;
				fK0sDedxVsPPion->Fill(prongTrackNeg->P(), prongTrackNeg->GetTPCsignal());
			}
			//cout << "End v0 PID" << endl;
			//Skip PID using TOF for now

			//cout << "V0 mass " << v0->MassK0Short() << ", Decay Length " << v0->DecayLength(fAODVertex) << ", Eta " << v0->Eta() << ", CosPA " << v0->CosPointingAngle(fAODVertex) << endl;
			fV0DecayLength->Fill(v0->DecayLength(fAODVertex));
			fV0Eta->Fill(v0->Eta());
			fV0CosPointingAngle->Fill(v0->CosPointingAngle(fAODVertex));
			fV0sMassK0s->Fill(v0->MassK0Short());
			if (!goodPiMinus || !goodPiPlus) continue;
			if (v0->Eta() > kEtaCut) continue;
			if (v0->CosPointingAngle(fAODVertex) < kMinCosAngle) continue;
			if (v0->DecayLength(fAODVertex) > kMaxDLK0) continue;
			if (v0->DecayLength(fAODVertex) < kMinDLK0) continue;


			//PRINT    cout << "v0 candidate selected" << endl;
			//We now have selected a V0 that has decayed to two pions
			double v0Dca = v0->DcaV0ToPrimVertex();
			fK0sDcaToPrimVertex->Fill(v0Dca);
			fK0sDaughterDcaToPrimVertex->Fill(v0->DcaNegToPrimVertex());
			fK0sDaughterDcaToPrimVertex->Fill(v0->DcaPosToPrimVertex());
			fK0DaughterDca->Fill(v0->DcaV0Daughters());
			fK0sMassDistribution->Fill(v0->MassK0Short());
			if (!fCutCheck) {
				//PRINT      cout << "cut check" << endl;
				if (v0->DcaNegToPrimVertex() < kMinDCAPrimaryPion) continue; //kMinDCAPrimaryPion = 0.4
				if (v0->DcaPosToPrimVertex() < kMinDCAPrimaryPion) continue;
				if (v0->DcaV0Daughters() > kMaxDCADaughtersK0) continue; //kMaxDCADaughtersK0 = 1.0
				if (v0Dca > kMaxDCAK0) continue;
				if (v0->MassK0Short() < 0.48 && v0->MassK0Short() > 0.515) continue;
				//cout << "cut check passed" << endl;
			}
			else {
				//cout << "No cut check" << endl;
				if (v0->DcaNegToPrimVertex() < kCheckDCAPi[2]) continue;
				if (v0->DcaPosToPrimVertex() < kCheckDCAPi[2]) continue;
				if (v0->DcaV0Daughters() > kCheckDCAPiPi[2]) continue;
				if (v0Dca > kCheckDCAK0[2]) continue;
				if (v0->MassK0Short() < kCheckMassLow[2] && v0->MassK0Short() > kCheckMassHigh[2]) continue;
				//cout << "no cut check passed" << endl;
			}

			//Skipping some MC stuff that is commented out

			//##### Check for the best K0s.
			//We will have one or two candidates only (unless there are v0s that might be from pile up.) We shall require:
			//(1) closest mass to K0s mass
			//(2) DCA K0s closest to primary vertex, and
			//(3) DCA daughters closest to each other (at point of decay).

			//I will add a histo fNumberV0sPerEvent in events, expecting it to be low and another with fNumber K0sPerEvent.
			//Will skip shared daughters check. I just want the best one.

			bool v0JudgeNew; //true if new v0 beats old
			double newV0Pars[3] = { 0., 0., 0. };
			double oldV0Pars[3] = { 0., 0., 0. };
			newV0Pars[0] = fabs(v0->MassK0Short() - k0ShortMass);
			newV0Pars[1] = v0Dca;
			newV0Pars[2] = v0->DcaV0Daughters(); //parameters used in merit cut

			if (i == 0) {
				oldV0Pars[0] = fabs(v0->MassK0Short() - k0ShortMass);
				oldV0Pars[1] = v0Dca;
				oldV0Pars[2] = v0->DcaV0Daughters(); //first time through, skip comparison
				bestK0sCandidateIndex = i;
				UsedNMissingDaughters = nMissingDaughters;
				//      cout << "v0 compared:" << i << endl;
			}
			else {
				v0JudgeNew = CheckMeritCutWinner(fMeritCutChoice, oldV0Pars, newV0Pars); //true if new wins
				if (v0JudgeNew) {
					oldV0Pars[0] = fabs(v0->MassK0Short() - k0ShortMass);
					oldV0Pars[1] = v0Dca;
					oldV0Pars[2] = v0->DcaV0Daughters(); //new beats old, set oldV0Par values for next pass
					bestK0sCandidateIndex = i; //update index of best K0s candidate
					UsedNMissingDaughters = nMissingDaughters; //Update number of missing daughters for the current, best K0s candidate
															   //cout << "v0 compared:" << i << endl;
				}
			}
			k0ShortCount++;

			//    cout << "new k0ShortCount: " << k0ShortCount << endl;
		}	//Just found the index for the best k0s candidate
		fHistK0sCandidatesPerEvent->Fill(k0ShortCount); //Diagnostic histo.
		if (k0ShortCount > 0) {
			fHistNeventsEtaCK0sChannel->Fill(1);

			//Now create data structures for Pions from best K0s candidate.
			AliAODv0* v0 = aod->GetV0(bestK0sCandidateIndex);
			AliAODTrack* tempTrack = (AliAODTrack*)v0->GetDaughter(0);
			short int pos0or1;
			short int neg0or1;
			bool orderswitch = kFALSE;
			if (tempTrack->Charge() > 0) { pos0or1 = 0; neg0or1 = 1; }
			else { pos0or1 = 1; neg0or1 = 0; orderswitch = kTRUE; }
			AliAODTrack* prongTrackPos = (AliAODTrack*)v0->GetDaughter(pos0or1);
			AliAODTrack* prongTrackNeg = (AliAODTrack*)v0->GetDaughter(neg0or1);

			//PRINT    cout << "D" << endl;

			//Six track loop
			UInt_t nK0sPion = 0, nK0GoodTracks = 0;
			nPion = 0; nKaon = 0;
			Int_t currentID = -666;
			Int_t posTrackIndex = -666;
			Int_t negTrackIndex = -666;
			Bool_t skipTrack = kFALSE;
			Int_t nProngFound = 0;
			Int_t trackIndex[7] = { -1,-1,-1,-1,-1,-1,-1 };
			Short_t qK0sPion[7];
			nGoodTracks = 0; 
			Double_t fRecTPCsignalK0sPion[7];
			Double_t dca[2] = { 0.0,0.0 }, cov[3] = { 0.0,0.0,0.0 };

			for (Int_t itr = 0; itr < aod->GetNumberOfTracks(); itr++) {
				AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
				if (!trk) continue;
				if (!(trk->TestFilterBit(1 << 0))) continue;
				if (trk->GetTPCNcls() < 50) continue;
				if (trk->Chi2perNDF() > 4) continue;
				currentID = trk->GetID();
				skipTrack = kFALSE;

				if (currentID == prongTrackPos->GetID()) { posTrackIndex = itr; skipTrack = kTRUE; nProngFound++; }
				if (currentID == prongTrackNeg->GetID()) { negTrackIndex = itr; skipTrack = kTRUE; nProngFound++; }
				if (!skipTrack) {
					AliAODTrack* trk_clone = (AliAODTrack*)trk->Clone("trk_clone");
					if (!(trk->GetStatus() & AliESDtrack::kTPCrefit)) continue;
					if (!(trk->GetStatus() & AliESDtrack::kITSrefit)) continue;
					if (!trk_clone->PropagateToDCA(fAODVertex, aod->GetMagneticField(), 300., dca, cov)) continue;
					delete trk_clone;
					if (TMath::Abs(dca[1]) > 2) continue;
					Double_t cut_DCAxy = 4 * (0.0182 + 0.0350 / TMath::Power(trk->Pt(), 1.01));
					if (TMath::Abs(dca[0]) > cut_DCAxy) continue;
					if ((trk->HasPointOnITSLayer(0)) || (trk->HasPointOnITSLayer(1))) nSpdHits++;
				}
				else {
					AliAODTrack* trk_clone = (AliAODTrack*)trk->Clone("trk_clone");
					Double_t cut_DCAxy = 4 * (0.0182 + 0.0350 / TMath::Power(trk->Pt(), 1.01));
					if ((trk->GetStatus() & AliESDtrack::kTPCrefit) && (trk->GetStatus() & AliESDtrack::kITSrefit) 
						&& trk_clone->PropagateToDCA(fAODVertex, aod->GetMagneticField(), 300., dca, cov) && TMath::Abs(dca[1]) < 2 && TMath::Abs(dca[0]) < cut_DCAxy ) {
						nK0GoodTracks++;
					}
					delete trk_clone;
				}
				trackIndex[nGoodTracks] = itr;
				nGoodTracks++;

				if (nGoodTracks > 6) break; //Alec's change -- was 4
			}//Track loop
			fK0GoodTracks->Fill(nK0GoodTracks);
			fHistNeventsEtaCK0sChannel->Fill(2);
			fHistNProngFound->Fill(nProngFound); //possible to find 6 good tracks before getting to the not good track which is the V0 ID'd earlier
			//PRINT    cout << "E" << endl;
			//now identify the additional kaon and pion, excluding the pions from K0s daughters by fID.
			firstGuess = kTRUE;
			TLorentzVector vK0sPion[6], vKPiK0sChannel, vK0s, vK0sPi, vK0sK;
			nKaon = 0; nPion = 0; nK0sPion = 0; nBoth = 0;
			qK0sPion[0] = 0; qK0sPion[1] = 0;
			//PID
			if (nProngFound == 2 && nGoodTracks == 6) { //Alec's change -- was 4 good tracks
				for (Int_t i = 0; i < 6; i++) { //changed from 4
					AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
					if (!trk) AliFatal("Not a standard AOD");

					//Here I need to identify Pions and Kaons.
					if (trackIndex[i] == posTrackIndex) {
						fRecTPCsignalK0sPion[0] = trk->GetTPCsignal();
						qK0sPion[nK0sPion] = trk->Charge();
						vK0sPion[nK0sPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
						nK0sPion++;
					}
					else if (trackIndex[i] == negTrackIndex) {
						fRecTPCsignalK0sPion[1] = trk->GetTPCsignal();
						qK0sPion[nK0sPion] = trk->Charge();
						vK0sPion[nK0sPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
						nK0sPion++;
					}
					else if (((fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && ((fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns pion and TOF returns pion, both, or missing
						|| (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns pion
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
						qPion[nPion] = trk->Charge();
						vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
						nPion++;
					}
					else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && ((fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.) || (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3.) || fPIDTOFPion[i] == -999.))  //if TPC returns kaon and TOF returns kaon, both, or missing
						|| (fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3. && fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3.))  //if TPC returns both and TOF returns kaon
						&& (beta > 0.55 || beta == -999.)) {
						fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
						qKaon[nKaon] = trk->Charge();
						vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
						nKaon++;
					}
					else if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)  //if TPC returns kaon and TOF returns pion
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)  //if TPC returns pion and TOF returns kaon
						|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.))) //if TPC and TOF both return both or missing
						&& (beta > 0.55 || beta == -999.)) {
						nBoth++;
					}

					if (nPion > 3 || nKaon > 1) break; //Alec's change -- was nPion > 1
				}
				fK0sEventCandidates->Fill(1);
				if ((qK0sPion[0] + qK0sPion[1]) == 0) fK0sEventCandidates->Fill(2);
				if (nPion == 3 && nKaon == 1 && nBoth == 0) fK0sEventCandidates->Fill(3);
				if (nPion == 2 && nKaon == 1 && nBoth == 1) fK0sEventCandidates->Fill(4);
				if (nPion == 3 && nKaon == 0 && nBoth == 1) fK0sEventCandidates->Fill(5);
				if (nPion == 1 && nKaon == 1 && nBoth == 2) fK0sEventCandidates->Fill(6);
				if (nPion == 2 && nKaon == 0 && nBoth == 2) fK0sEventCandidates->Fill(7);
				if (nPion == 0 && nKaon == 1 && nBoth == 3) fK0sEventCandidates->Fill(8);
				fHistNK0sPion->Fill(nK0sPion);
				for (Int_t pidLoopCounter = 0; pidLoopCounter < 1; pidLoopCounter++) {
					newKaon = 0; newPion = 0; nHighPtTracks = 0;

					//##### Deal with tracks missing TOF information. #####
					//Any events missing TOF info which can be completed to have the products for this decay channel will be assumed to belong to this channel.
					if (nPion == 3 && nBoth == 1) {
						for (int i = 0; i < 6; i++) { //loop over the four tracks
							if ((trackIndex[i] == posTrackIndex) || (trackIndex[i] == negTrackIndex)) continue;
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
								qKaon[nKaon] = trk->Charge();
								vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
								nKaon++;
								newKaon++;
							}
						}
					}
					if (((nKaon - newKaon) == 1 && nBoth == 3) || ((nKaon - newKaon) == 1 && nPion == 1 && nBoth == 2) || ((nKaon - newKaon) == 1 && nPion == 2 && nBoth == 1)) {
						for (int i = 0; i < 6; i++) { //loop over the four tracks
							if ((trackIndex[i] == posTrackIndex) || (trackIndex[i] == negTrackIndex)) continue;
							AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
							if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
								|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
								&& (beta > 0.55 || beta == -999.)) {
								fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
								qPion[nPion] = trk->Charge();
								vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
								nPion++;
								newPion++;
							}
						}
					}
					if ((nKaon - newKaon) == 0 && (nPion - newPion) == 2 && nBoth == 2) {
						//cout << "One Kaon, One Pion with two Missing tracks" << endl;
						if (firstGuess) {
							//cout << "First guess" << endl;
							for (int i = 0; i < 6; i++) { //loop over the four tracks
								if ((trackIndex[i] == posTrackIndex) || (trackIndex[i] == negTrackIndex)) continue;
								if (nPion == 2) {
									AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
									if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
										&& (beta > 0.55 || beta == -999.)) {
										fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
										qPion[nPion] = trk->Charge();
										vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
										nPion++;
										newPion++; //added by Alec
									}
								}
								else {
									AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
									if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
										&& (beta > 0.55 || beta == -999.)) {
										fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
										qKaon[nKaon] = trk->Charge();
										vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
										nKaon++;
										newKaon++; //added by Alec
									}
								}
								//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
							}
							pidLoopCounter--;
							firstGuess = kFALSE;
						}
						else {
							//cout << "Second Guess" << endl;
							for (int i = 0; i < 6; i++) { //loop over the four tracks
								if ((trackIndex[i] == posTrackIndex) || (trackIndex[i] == negTrackIndex)) continue;
								if (nKaon == 0) {
									AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
									if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
										&& (beta > 0.55 || beta == -999.)) {
										fRecTPCsignalKaon[nKaon] = trk->GetTPCsignal();
										qKaon[nKaon] = trk->Charge();
										vKaon[nKaon].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), kaonMass);
										nKaon++;
										newKaon++; //added by Alec
									}
								}
								else {
									AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[i]));
									if (((fabs(fPIDTPCPion[i]) > 3. && fabs(fPIDTPCKaon[i]) < 3. && fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) > 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) > 3. && fabs(fPIDTOFPion[i]) > 3. && fabs(fPIDTOFKaon[i]) < 3.)
										|| (fabs(fPIDTPCPion[i]) < 3. && fabs(fPIDTPCKaon[i]) < 3. && (fabs(fPIDTOFPion[i]) < 3. && fabs(fPIDTOFKaon[i]) < 3. || fPIDTOFPion[i] == -999.)))
										&& (beta > 0.55 || beta == -999.)) {
										fRecTPCsignalPion[nPion] = trk->GetTPCsignal();
										qPion[nPion] = trk->Charge();
										vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
										nPion++;
										newPion++; //added by Alec
									}
								}
								//cout << "Good Track #" << i + 1 << ", " << nPion << " Pions and " << nKaon << " Kaons ID'd" << endl;
							}
						}
					}

					if (nPion == 3 && nK0sPion == 2 && nKaon == 1) {
						fHistNeventsEtaCK0sChannel->Fill(3);
						if (qKaon[0] < 0 && (qPion[0] + qPion[1] + qPion[2]) < 0) fHistNeventsEtaCK0sChannel->Fill(4); //K- and sumQPions < 0 K-pi+pi-pi- or K-pi-pi-pi-
						if (qKaon[0] > 0 && (qPion[0] + qPion[1] + qPion[2]) > 0) fHistNeventsEtaCK0sChannel->Fill(5); //K+ and sumQPions > 0 K+pi+pi+pi- or K+pi+pi+pi+
						if (qKaon[0] < 0 && (qPion[0] + qPion[1] + qPion[2]) > 0 && (qK0sPion[0] + qK0sPion[1]) == 0) fHistNeventsEtaCK0sChannel->Fill(6); //Requires K-Pi+Pi+Pi- (K- events)
						if (qKaon[0] > 0 && (qPion[0] + qPion[1] + qPion[2]) < 0 && (qK0sPion[0] + qK0sPion[1]) == 0) fHistNeventsEtaCK0sChannel->Fill(7); //Requires K+Pi-Pi+Pi- (K+ events)
																																							   //Continue for good cases
						if (((qKaon[0] < 0) && ((qPion[0] + qPion[1] + qPion[2]) > 0) && ((qK0sPion[0] + qK0sPion[1]) == 0)) || ((qKaon[0] > 0) && ((qPion[0] + qPion[1] + qPion[2]) < 0) && ((qK0sPion[0] + qK0sPion[1]) == 0))) {

							//PRINT	  cout << "Filling K0s histos" << endl;

							fHistNeventsEtaCK0sChannel->Fill(8);
							fEtaCCandidatesPerChannel->Fill(8);
							//Now we need to use the V0s to identify a K0s and select best K0s candidate.
							//Individual Track Ptś
							fK0sPosDaughterPt->Fill(vK0sPion[0].Pt());
							fK0sNegDaughterPt->Fill(vK0sPion[1].Pt());
							fK0sPosVsNegDaughterPt->Fill(vK0sPion[1].Pt(), vK0sPion[0].Pt());
							fK0sPionPt->Fill(vPion[0].Pt());
							fK0sKaonPt->Fill(vKaon[0].Pt());
							//Compute K0s info
							vK0s = vK0sPion[0] + vK0sPion[1];
							fK0sPtVsMinvK0s->Fill(vK0s.M(), vK0s.Pt());
							//fK0sMinv->Fill(vK0s.M());
							//Compute PiK info
							vKPiK0sChannel = vPion[0] + vPion[1] + vPion[2] + vKaon[0];
							fK0sPtVsMinvKPi->Fill(vKPiK0sChannel.M(), vKPiK0sChannel.Pt());
							//fKPiMinvK0sChannel->Fill(vKPiK0sChannel.M());
							//Dalitz plot K0s vs PiK combo
							vK0sPi = vK0s + vPion[0];
							vK0sK = vK0s + vKaon[0];
							fK0sM2K0sVsM2KPi->Fill(vKPiK0sChannel.M()*vKPiK0sChannel.M(), vK0s.M()*vK0s.M());
							fK0sM2K0sPiVsM2KPi->Fill(vKPiK0sChannel.M()*vKPiK0sChannel.M(), vK0sPi.M()*vK0sPi.M());
							fK0sM2K0sKVsM2KPi->Fill(vKPiK0sChannel.M()*vKPiK0sChannel.M(), vK0sK.M()*vK0sK.M());
							//Compute EtaC info
							fHistNeventsEtaC->Fill(4);

							vCandidate = vK0sPion[0] + vK0sPion[1] + vPion[0] + vPion[1] + vPion[2] + vKaon[0];
							if (vCandidate.Pt() < 0.11) fEtaCLowPtCandidatesPerChannel->Fill(8);

							//Sum Pz
							SumPz = vK0sPion[0].Pz() + vK0sPion[1].Pz() + vPion[0].Pz() + vPion[1].Pz() + vPion[2].Pz() + vKaon[0].Pz();
							//Scalar Sum P
							ScalarSumP = vK0sPion[0].P() + vK0sPion[1].P() + vPion[0].P() + vPion[1].P() + vPion[2].P() + vKaon[0].P();
							//Vector Sum Pt
							sumPtVector = vK0sPion[0].Vect() + vK0sPion[1].Vect() + vPion[0].Vect() + vPion[1].Vect() + vPion[2].Vect() + vKaon[0].Vect();
							VectorSumPt = sumPtVector.Pt();
							//Number of tracks with pT>0.4 GeV/c
							for (Int_t aa = 0; aa < 2; aa++) if (vK0sPion[aa].Pt() > 0.4) nHighPtTracks++;
							for (Int_t aa = 0; aa < 3; aa++) if (vPion[aa].Pt() > 0.4) nHighPtTracks++;
							if (vKaon[0].Pt() > 0.4) nHighPtTracks++;

							fK0sPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
							fAllPtVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Pt());
							if (vCandidate.Pt() < 0.11) {
								fChannelVsMinvEtaC->Fill(vCandidate.M(), 8);
								fAllMinvEtaCLowPt->Fill(vCandidate.M());
								if (vCandidate.M() > (etaCMass - etaCWidth) && vCandidate.M() < (etaCMass + etaCWidth)) fHistNeventsEtaC->Fill(5);
							}
							fK0sEtaVsMinvEtaC->Fill(vCandidate.M(), vCandidate.Eta());
							if (vCandidate.Pt() < 0.4) fK0sEtaVsMinvEtaC400MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
							if (vCandidate.Pt() < 0.1) fK0sEtaVsMinvEtaC100MeVPtMax->Fill(vCandidate.M(), vCandidate.Eta());
							fK0sSumPzVsMinvEtaC->Fill(vCandidate.M(), SumPz);
							fK0sScalarSumP->Fill(ScalarSumP);
							fK0sVectorSumPt->Fill(VectorSumPt);
							if (nHighPtTracks > 1) fHistNeventsEtaCK0sChannel->Fill(9);
							//fEtaCMinvK0sChannel->Fill(vCandidate.M());
							fK0sDecayLength->Fill(v0->DecayLength(fAODVertex));
						}
					}
				}
			}
		}
	}
	//  cout << "##### End of RunAODHist()" << endl;
	
	PostData(2, fListHist);
	PostData(3, fListHistKstar);
	PostData(4, fListHist2Rho4Pion);
	PostData(5, fListHistK0s3PiPi4K);
}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunAODtree()
{

}//RunAOD

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunAODMC(AliAODEvent *aod)
{

}//RunAODMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunESDtrig()
{

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunESDhist()
{

}

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunESDtree()
{

}//RunESD

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::RunESDMC(AliESDEvent* esd)
{

}//RunESDMC

//_____________________________________________________________________________
void AliAnalysisTaskUpcEtaCAWP::Terminate(Option_t *) 
{

  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
Double_t AliAnalysisTaskUpcEtaCAWP::GetMedian(Double_t *daArray) {
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

bool AliAnalysisTaskUpcEtaCAWP::CheckMeritCutWinner(int cutChoice, double oldPars[3], double newPars[3]){
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

void AliAnalysisTaskUpcEtaCAWP::BoostCut(TLorentzVector d1, TLorentzVector d2, TLorentzVector parent, Double_t *boostInfo) {

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

void AliAnalysisTaskUpcEtaCAWP::RunAODsystematics(AliAODEvent* aod)
{

  /*
  Double_t fJPsiSels[4];

  fJPsiSels[0] =   70; //min number of TPC clusters
  fJPsiSels[1] =   4; //chi2
  fJPsiSels[2] =   2; //DCAz
  fJPsiSels[3] =   1; // DCAxy 1x 

  Double_t fJPsiSelsMid[4];

  fJPsiSelsMid[0] =   70; //min number of TPC clusters
  fJPsiSelsMid[1] =   4; //chi2
  fJPsiSelsMid[2] =   2; //DCAz
  fJPsiSelsMid[3] =   1; // DCAxy 1x 
  
  Double_t fJPsiSelsLoose[4];

  fJPsiSelsLoose[0] =   60; //min number of TPC clusters
  fJPsiSelsLoose[1] =   5; //chi2
  fJPsiSelsLoose[2] =   3; //DCAz
  fJPsiSelsLoose[3] =   2; // DCAxy 2x 

  Double_t fJPsiSelsTight[4];

  fJPsiSelsTight[0] =   80; //min number of TPC clusters
  fJPsiSelsTight[1] =   3.5; //chi2
  fJPsiSelsTight[2] =   1; //DCAz
  fJPsiSelsTight[3] =   0.5; // DCAxy 0.5x 

  Int_t nGoodTracks = 0;
  Int_t trackIndex[5] = {-1,-1,-1,-1,-1};
  
  TLorentzVector vLepton[4], vPion[4], vCandidate, vDilepton;
  Short_t qLepton[4],qPion[4];
  UInt_t nLepton=0, nPion=0, nHighPt=0;
  Double_t fRecTPCsignal[5], fRecTPCsignalDist;
  Int_t fChannel = 0;

  AliAODVertex *fAODVertex = aod->GetPrimaryVertex();
  
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  
  TParticlePDG *partKaon = pdgdat->GetParticle( 321 );
  Double_t kaonMass = partKaon->Mass();
    
  TParticlePDG *partPion = pdgdat->GetParticle( 211 );
  Double_t pionMass = partPion->Mass();

  TParticlePDG *partKstar = pdgdat->GetParticle( 313 );
  Double_t kStarMass = partKstar->Mass();
  Double_t kStarWidth = partKstar->Width();

  TParticlePDG *partK0short = pdgdat->GetParticle( 310 );
  Double_t k0ShortMass = partK0short->Mass();
  Double_t k0ShortWidth = partK0short->Width();

  
for(Int_t i=0; i<5; i++){
	  //cout<<"Loose sytematics, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsLoose[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
  //Two track loop
  nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if(i!=4){ if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;}
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
     
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > 1) nHighPt++;     
      		fRecTPCsignal[nLepton] = trk->GetTPCsignal();     
      		qLepton[nLepton] = trk->Charge();
      		if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				mass[nLepton] = 0;
				}
      		if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				mass[nLepton] = 1;
				}
       		nLepton++;		
  		}		
  	if(nLepton == 2){
		if(qLepton[0]*qLepton[1] < 0 && nHighPt > 0 && (mass[0]!=-1 || mass[1]!=-1)){
			vCandidate = vLepton[0]+vLepton[1];		  
  			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  			if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  			else { 
				fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  				if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
				}
			if(fChannel == -1 && vCandidate.Pt()<0.15) ((TH1D*)(fListJPsiLoose->At(i)))->Fill(vCandidate.M()); 
  			if(fChannel == 1 && vCandidate.Pt()<0.3) ((TH1D*)(fListJPsiLoose->At(i)))->Fill(vCandidate.M()); 	
			}
		}
  }
}//loose cuts

for(Int_t i=0; i<4; i++){
	  //cout<<"Tight sytematics, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsTight[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
  //Two track loop
  nGoodTracks = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((!trk->HasPointOnITSLayer(0))&&(!trk->HasPointOnITSLayer(1))) continue;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
     
      trackIndex[nGoodTracks] = itr;
      nGoodTracks++;
				  
      if(nGoodTracks > 2) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nHighPt=0;
  
  if(nGoodTracks == 2){
  	  for(Int_t k=0; k<2; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");
    
      		if(trk->Pt() > 1) nHighPt++;     
      		fRecTPCsignal[nLepton] = trk->GetTPCsignal();     
      		qLepton[nLepton] = trk->Charge();
      		if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
				mass[nLepton] = 0;
				}
      		if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      				vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
				mass[nLepton] = 1;
				}
       		nLepton++;		
  		}		
  	if(nLepton == 2){
		if(qLepton[0]*qLepton[1] < 0 && nHighPt > 0 && (mass[0]!=-1 || mass[1]!=-1)){
			vCandidate = vLepton[0]+vLepton[1];		  
  			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  			if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  			else { 
				fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  				if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
				}
			if(fChannel == -1 && vCandidate.Pt()<0.15) ((TH1D*)(fListJPsiTight->At(i)))->Fill(vCandidate.M()); 
  			if(fChannel == 1 && vCandidate.Pt()<0.3) ((TH1D*)(fListJPsiTight->At(i)))->Fill(vCandidate.M()); 	
			}
		}
  }
}//tight cuts

//---------------------------------------------EtaC------------------------------------------------------------------------

  Double_t fEtaCSels[4];

  fEtaCSels[0] =   50; //min number of TPC clusters
  fEtaCSels[1] =   4; //chi2
  fEtaCSels[2] =   2; //DCAz
  fEtaCSels[3] =   4; // DCAxy 1x 

  Double_t fEtaCSelsMid[4];

  fEtaCSelsMid[0] =   50; //min number of TPC clusters
  fEtaCSelsMid[1] =   4; //chi2
  fEtaCSelsMid[2] =   2; //DCAz
  fEtaCSelsMid[3] =   4; // DCAxy 1x 
  
  Double_t fEtaCSelsLoose[4];

  fEtaCSelsLoose[0] =   60; //min number of TPC clusters
  fEtaCSelsLoose[1] =   5; //chi2
  fEtaCSelsLoose[2] =   3; //DCAz
  fEtaCSelsLoose[3] =   6; // DCAxy 2x 

  Double_t fEtaCSelsTight[4];

  fEtaCSelsTight[0] =   70; //min number of TPC clusters
  fEtaCSelsTight[1] =   3.5; //chi2
  fEtaCSelsTight[2] =   1; //DCAz
  fEtaCSelsTight[3] =   2; // DCAxy 0.5x 

  nGoodTracks = 0; nLepton=0; nHighPt=0; fChannel = 0;
  Int_t nSpdHits = 0;
  Double_t trackPt[5]={0,0,0,0,0};

for(Int_t i=0; i<5; i++){
	  //cout<<"Loose systematics psi2s, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsLoose[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
 
  //Four track loop
  nGoodTracks = 0; nSpdHits = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
  nLepton=0; nPion=0; nHighPt=0;
  
  if(nGoodTracks == 4){
  	  if(i!=4){ if(nSpdHits<2) continue;} 
    	  MeanPt = GetMedian(trackPt);
  	  for(Int_t k=0; k<4; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > MeanPt){   
      			fRecTPCsignal[nLepton] = trk->GetTPCsignal();      
      			qLepton[nLepton] = trk->Charge();
      			if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
					mass[nLepton] = 0;
					}
      			if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
					mass[nLepton] = 1;
					}
			nLepton++;
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}	      
    		}
	if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0) && mass[0] != -1 && mass[1] != -1){
  		vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  		if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  		else { 
			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  			if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
			}			
		if(fChannel == -1) if(vDilepton.M() > 3.0 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.15) ((TH1D*)(fListEtaCLoose->At(i)))->Fill(vCandidate.M());		
  		if(fChannel == 1) if(vDilepton.M() > 2.6 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.3) ((TH1D*)(fListEtaCLoose->At(i)))->Fill(vCandidate.M());
	}
  }   
}//loose cuts

for(Int_t i=0; i<4; i++){
	  //cout<<"Tight systematics psi2s, cut"<<i<<endl;
	  for(Int_t j=0; j<4; j++){
		  if(i==j) fJPsiSels[j] = fJPsiSelsTight[i];
		  else fJPsiSels[j] = fJPsiSelsMid[j];
	  }
 
  //Four track loop
  nGoodTracks = 0; nSpdHits = 0;
  for(Int_t itr=0; itr<aod ->GetNumberOfTracks(); itr++) {
    AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(itr));
    if( !trk ) continue;
    if(!(trk->TestFilterBit(1<<0))) continue;

      if(!(trk->GetStatus() & AliESDtrack::kTPCrefit) ) continue;
      if(!(trk->GetStatus() & AliESDtrack::kITSrefit) ) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* trk_clone=(AliAODTrack*)trk->Clone("trk_clone");
      if(!trk_clone->PropagateToDCA(fAODVertex,aod->GetMagneticField(),300.,dca,cov)) continue;
      delete trk_clone;
      Double_t cut_DCAxy = (0.0182 + 0.0350/TMath::Power(trk->Pt(),1.01));
      
      if(trk->GetTPCNcls() < fJPsiSels[0])continue;
      if(trk->Chi2perNDF() > fJPsiSels[1])continue;
      if(TMath::Abs(dca[1]) > fJPsiSels[2]) continue;      
      if(TMath::Abs(dca[0]) > fJPsiSels[3]*cut_DCAxy) continue;
      if((trk->HasPointOnITSLayer(0))||(trk->HasPointOnITSLayer(1))) nSpdHits++;
     
      trackIndex[nGoodTracks] = itr;
      trackPt[nGoodTracks] = trk->Pt();
      nGoodTracks++;
				  
      if(nGoodTracks > 4) break;  
  }//Track loop
    
  Int_t mass[3]={-1,-1,-1};
  fChannel = 0;
    nLepton=0; nPion=0; nHighPt=0;
  
  if(nGoodTracks == 4){
  	  if(nSpdHits<2) continue; 
    	  MeanPt = GetMedian(trackPt);
  	  for(Int_t k=0; k<4; k++){
                AliAODTrack *trk = dynamic_cast<AliAODTrack*>(aod->GetTrack(trackIndex[k]));
                if(!trk) AliFatal("Not a standard AOD");

      		if(trk->Pt() > MeanPt){   
      			fRecTPCsignal[nLepton] = trk->GetTPCsignal();      
      			qLepton[nLepton] = trk->Charge();
      			if(fRecTPCsignal[nLepton] > 40 && fRecTPCsignal[nLepton] < 70){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
					mass[nLepton] = 0;
					}
      			if(fRecTPCsignal[nLepton] > 70 && fRecTPCsignal[nLepton] < 100){
      					vLepton[nLepton].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
					mass[nLepton] = 1;
					}
			nLepton++;
			}
		else{
			qPion[nPion] = trk->Charge();
			vPion[nPion].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);
			nPion++;
			}	      
    		}
	if((qLepton[0]*qLepton[1] < 0) && (qPion[0]*qPion[1] < 0) && mass[0] != -1 && mass[1] != -1){
  		vCandidate = vLepton[0]+vLepton[1]+vPion[0]+vPion[1];
  		vDilepton = vLepton[0]+vLepton[1];
		fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-56,2)+TMath::Power(fRecTPCsignal[1]-56,2));
  		if (fRecTPCsignalDist < 3.6*4.0) fChannel = -1;
  		else { 
			fRecTPCsignalDist = TMath::Sqrt(TMath::Power(fRecTPCsignal[0]-78,2)+TMath::Power(fRecTPCsignal[1]-78,2));
  			if (fRecTPCsignalDist < 4.1*4.0) fChannel = 1; 
			}			
		if(fChannel == -1) if(vDilepton.M() > 3.0 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.15) ((TH1D*)(fListEtaCTight->At(i)))->Fill(vCandidate.M());		
  		if(fChannel == 1) if(vDilepton.M() > 2.6 && vDilepton.M() < 3.2 && vCandidate.Pt()<0.3) ((TH1D*)(fListEtaCTight->At(i)))->Fill(vCandidate.M());
	}
  }   
}//Tight cuts
  */
}
