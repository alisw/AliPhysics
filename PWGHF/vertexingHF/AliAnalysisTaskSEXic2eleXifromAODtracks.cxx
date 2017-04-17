/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//               Xic->eXi analysis code
//
//  Input: AOD
//  Output: TTree or THnSparse (mass vs pT vs Centrality)
//
//-------------------------------------------------------------------------
//
//                 Authors: Y.S Watanabe(a)
//  (a) CNS, the University of Tokyo
//  Contatcs: wyosuke@cns.s.u-tokyo.ac.jp
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEXic2eleXifromAODtracks.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTOFPIDResponse.h"
#include "AliAODPidHF.h"
#include "AliInputEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliNeutralTrackParam.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliVertexerTracks.h"
#include "AliEventPoolManager.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEXic2eleXifromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEXic2eleXifromAODtracks::AliAnalysisTaskSEXic2eleXifromAODtracks() : 
  AliAnalysisTaskSE(),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHEventPlane(0),
  fHNTrackletvsZ(0),
  fHNTrackletCorrvsZ(0),
  fAnalCuts(0),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(kFALSE),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fEleVariablesTree(0),
  fCascVariablesTree(0),
  fSingleVariablesTree(0),
  fMCVariablesTree(0),
  fMCEleVariablesTree(0),
  fMCCascVariablesTree(0),
  fMCGenPairVariablesTree(0),
  fCorrelationVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateCascVariables(),
  fCandidateSingleVariables(),
  fCandidateMCVariables(),
  fCandidateMCEleVariables(),
  fCandidateMCCascVariables(),
  fCandidateMCGenPairVariables(),
  fCorrelationVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fEventPlane(0),
  fQ(0),
  fQSub1(0),
  fQSub2(0),
  fBzkG(0),
  fCentrality(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fUseCentralitySPDTracklet(kFALSE),
  fUseEventPlane(0),
  fEvNumberCounter(0),
  fMCEventType(-9999),
  fMCDoPairAnalysis(kFALSE),
  fHistoEleXiMass(0),
  fHistoEleXiMassRS(0),
  fHistoEleXiMassWS(0),
  fHistoEleXiMassRSMix(0),
  fHistoEleXiMassWSMix(0),
  fHistoEleXiMassRSSide(0),
  fHistoEleXiMassWSSide(0),
  fHistoEleXiMassRS1(0),
  fHistoEleXiMassWS1(0),
  fHistoEleXiMassRSMix1(0),
  fHistoEleXiMassWSMix1(0),
  fHistoEleXiMassRSSide1(0),
  fHistoEleXiMassWSSide1(0),
  fHistoEleXiMassRS2(0),
  fHistoEleXiMassWS2(0),
  fHistoEleXiMassRSMix2(0),
  fHistoEleXiMassWSMix2(0),
  fHistoEleXiMassRSSide2(0),
  fHistoEleXiMassWSSide2(0),
  fHistoEleXiMassAway(0),
  fHistoEleXiMassRSAway(0),
  fHistoEleXiMassWSAway(0),
  fHistoEleXiMassRSMixAway(0),
  fHistoEleXiMassWSMixAway(0),
  fHistoEleXiMassRSSideAway(0),
  fHistoEleXiMassWSSideAway(0),
  fHistoEleXiMassRS1Away(0),
  fHistoEleXiMassWS1Away(0),
  fHistoEleXiMassRSMix1Away(0),
  fHistoEleXiMassWSMix1Away(0),
  fHistoEleXiMassRSSide1Away(0),
  fHistoEleXiMassWSSide1Away(0),
  fHistoEleXiMassRS2Away(0),
  fHistoEleXiMassWS2Away(0),
  fHistoEleXiMassRSMix2Away(0),
  fHistoEleXiMassWSMix2Away(0),
  fHistoEleXiMassRSSide2Away(0),
  fHistoEleXiMassWSSide2Away(0),
  fHistoEleXiMassvsElePtRS(0),
  fHistoEleXiMassvsElePtWS(0),
  fHistoEleXiMassvsElePtRSMix(0),
  fHistoEleXiMassvsElePtWSMix(0),
  fHistoEleXiMassvsElePtRSSide(0),
  fHistoEleXiMassvsElePtWSSide(0),
  fHistoEleXiMassvsElePtRS1(0),
  fHistoEleXiMassvsElePtWS1(0),
  fHistoEleXiMassvsElePtRSMix1(0),
  fHistoEleXiMassvsElePtWSMix1(0),
  fHistoEleXiMassvsElePtRSSide1(0),
  fHistoEleXiMassvsElePtWSSide1(0),
  fHistoEleXiMassvsElePtRS2(0),
  fHistoEleXiMassvsElePtWS2(0),
  fHistoEleXiMassvsElePtRSMix2(0),
  fHistoEleXiMassvsElePtWSMix2(0),
  fHistoEleXiMassvsElePtRSSide2(0),
  fHistoEleXiMassvsElePtWSSide2(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleXiMassMCS(0),
  fHistoEleXiMassMCS1(0),
  fHistoEleXiMassMCS2(0),
  fHistoEleXiMassXibMCS(0),
  fHistoEleXiMassXibMCS1(0),
  fHistoEleXiMassXibMCS2(0),
  fHistoEleXiMassPromptMCS(0),
  fHistoEleXiMassPromptMCS1(0),
  fHistoEleXiMassPromptMCS2(0),
  fHistoEleXiMassBFeeddownMCS(0),
  fHistoEleXiMassBFeeddownMCS1(0),
  fHistoEleXiMassBFeeddownMCS2(0),
  fHistoEleXiMassMCGen(0),
  fHistoEleXiMassvsElePtMCS(0),
  fHistoEleXiMassvsElePtMCGen(0),
  fHistoEleXiMassvsElePtMCS1(0),
  fHistoEleXiMassvsElePtMCGen1(0),
  fHistoEleXiMassvsElePtMCS2(0),
  fHistoEleXiMassvsElePtMCGen2(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsXiPtRS(0),
  fHistoElePtvsXiPtWS(0),
  fHistoElePtvsXiPtRSMix(0),
  fHistoElePtvsXiPtWSMix(0),
  fHistoElePtvsXiPtMCS(0),
  fHistoElePtvsXiPtvsXicPtMCS(0),
  fHistoElePtvsXiPtMCGen(0),
  fHistoElePtvsXiPtvsXicPtMCGen(0),
  fHistoElePtvsXiPtMCXicGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoXiMassvsPtMCS(0),
  fHistoXiMassvsPtMCGen(0),
  fHistoOmegaMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTOFPIDSelTPC(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
  fHistoMassConversionsMin(0),
  fHistoMassConversionsSameSignMin(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoXiQovPtvsPhi(0),
	fHistoXicMCGen(0),
	fHistoXicMCGen1(0),
	fHistoXicMCGen2(0),
	fHistoXicMCS(0),
	fHistoXicMCS1(0),
	fHistoXicMCS2(0),
	fHistoXibMCGen(0),
	fHistoXibMCGenWithXic(0),
	fHistoXibMCS(0),
	fHistoXicElectronMCGen(0),
	fHistoXicElectronMCGen1(0),
	fHistoXicElectronMCGen2(0),
	fHistoXicElectronMCS(0),
	fHistoXicElectronMCS1(0),
	fHistoXicElectronMCS2(0),
	fHistoElectronMCGen(0),
	fHistoBottomElectronMCGen(0),
	fHistoCharmElectronMCGen(0),
	fHistoXiMCGen(0),
	fHistoLambdaPtvsDl(0),
	fHistoLambdaPtvsDlSide(0),
	fHistoLambdaPtvsDlMCS(0),
	fHistoLambdaPtvsDR(0),
	fHistoLambdaPtvsDRSide(0),
	fHistoLambdaPtvsDRMCS(0),
	fHistoEleXiPtvsRapidityRS(0),
	fHistoEleXiPtvsRapidityWS(0),
	fHistoEleXiPtvsRapidityMCS(0),
	fHistoCorrelationVariablesvsEleXiPt(0),
	fHistoCorrelationVariablesvsEleXiPtMix(0),
	fHistoCorrelationVariablesvsEleXiPtMC(0),
	fHistoCorrelationVariablesvsElePt(0),
	fHistoCorrelationVariablesvsElePtMix(0),
	fHistoCorrelationVariablesvsElePtMC(0),
	fHistoCorrelationVariablesvsXiPt(0),
	fHistoCorrelationVariablesvsXiPtMix(0),
	fHistoCorrelationVariablesvsXiPtMC(0),
	fHistoMassVariablesvsEleXiPt(0),
	fHistoMassVariablesvsEleXiPtMix(0),
	fHistoMassVariablesvsEleXiPtMC(0),
	fHistoMassVariablesvsElePt(0),
	fHistoMassVariablesvsElePtMix(0),
	fHistoMassVariablesvsElePtMC(0),
	fHistoMassVariablesvsXiPt(0),
	fHistoMassVariablesvsXiPtMix(0),
	fHistoMassVariablesvsXiPtMC(0),
	fHistoResponseElePt(0),
	fHistoResponseXiPt(0),
	fHistoResponseEleXiPt(0),
	fHistoResponseXiPtvsEleXiPt(0),
	fHistoResponseXiPtXib(0),
	fHistoResponseEleXiPtXib(0),
	fHistoResponseMcGenXibPtvsXicPt(0),
  fHistoPi0MCGen(0),
  fHistoElectronPi0Total(0),
  fHistoElectronPi0Tag(0),
  fHistoEtaMCGen(0),
  fHistoElectronEtaTotal(0),
  fHistoElectronEtaTag(0),
  fHistoKaonMCGen(0),
  fHistoD0MCGen(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonXivsRunNumber(0),
	fHistoMCEventType(0),
	fHistoMCXic0Decays(0),
	fHistoMCDeltaPhiccbar(0),
	fHistoMCNccbar(0),
	fRefMult(9.26),
  fGTI(0),fGTIndex(0), fTrackBuffSize(19000),
  fHistodPhiSdEtaSElectronProtonR125RS(0),
  fHistodPhiSdEtaSElectronProtonR125WS(0),
  fHistodPhiSdEtaSElectronProtonR125RSMix(0),
  fHistodPhiSdEtaSElectronProtonR125WSMix(0),
  fHistodPhiSdEtaSElectronPionR125RS(0),
  fHistodPhiSdEtaSElectronPionR125WS(0),
  fHistodPhiSdEtaSElectronPionR125RSMix(0),
  fHistodPhiSdEtaSElectronPionR125WSMix(0),
  fHistodPhiSdEtaSElectronBachelorR125RS(0),
  fHistodPhiSdEtaSElectronBachelorR125WS(0),
  fHistodPhiSdEtaSElectronBachelorR125RSMix(0),
  fHistodPhiSdEtaSElectronBachelorR125WSMix(0),
  fDoEventMixing(0),
  fMixWithoutConversionFlag(kFALSE),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
  fNRPBins					(0), 
	fNOfPools(1),
  fPoolIndex(-9999),
  nextResVec(),
  reservoirsReady(),
  m_ReservoirE(),
  m_ReservoirL1(),
  m_ReservoirL2(),
  m_ReservoirVarsE(),
  m_ReservoirVarsL1(),
  m_ReservoirVarsL2()
{
  //
  // Default Constructor. 
  //
	for(Int_t i=0;i<23;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
	}
	for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}

//___________________________________________________________________________
AliAnalysisTaskSEXic2eleXifromAODtracks::AliAnalysisTaskSEXic2eleXifromAODtracks(const Char_t* name,
									     AliRDHFCutsXictoeleXifromAODtracks* analCuts, 
									     Bool_t writeVariableTree) :
  AliAnalysisTaskSE(name),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputAll(0),
  fListCuts(0),
  fCEvents(0),
  fHTrigger(0),
  fHCentrality(0),
  fHEventPlane(0),
  fHNTrackletvsZ(0),
  fHNTrackletCorrvsZ(0),
  fAnalCuts(analCuts),
  fIsEventSelected(kFALSE),
  fWriteVariableTree(writeVariableTree),
  fWriteEachVariableTree(kFALSE),
  fWriteMCVariableTree(kFALSE),
  fVariablesTree(0),
  fEleVariablesTree(0),
  fCascVariablesTree(0),
  fSingleVariablesTree(0),
  fMCVariablesTree(0),
  fMCEleVariablesTree(0),
  fMCCascVariablesTree(0),
  fMCGenPairVariablesTree(0),
  fCorrelationVariablesTree(0),
  fReconstructPrimVert(kFALSE),
  fIsMB(kFALSE),
  fIsSemi(kFALSE),
  fIsCent(kFALSE),
  fIsINT7(kFALSE),
  fIsEMC7(kFALSE),
  fCandidateVariables(),
  fCandidateEleVariables(),
  fCandidateCascVariables(),
  fCandidateSingleVariables(),
  fCandidateMCVariables(),
  fCandidateMCEleVariables(),
  fCandidateMCCascVariables(),
  fCandidateMCGenPairVariables(),
  fCorrelationVariables(),
  fVtx1(0),
  fV1(0),
  fVtxZ(0),
  fEventPlane(0),
  fQ(0),
  fQSub1(0),
  fQSub2(0),
  fBzkG(0),
  fCentrality(0),
  fRunNumber(0),
  fTriggerCheck(0),
  fUseCentralityV0M(kFALSE),
  fUseCentralitySPDTracklet(kFALSE),
  fUseEventPlane(0),
  fEvNumberCounter(0),
  fMCEventType(-9999),
  fMCDoPairAnalysis(kFALSE),
  fHistoEleXiMass(0),
  fHistoEleXiMassRS(0),
  fHistoEleXiMassWS(0),
  fHistoEleXiMassRSMix(0),
  fHistoEleXiMassWSMix(0),
  fHistoEleXiMassRSSide(0),
  fHistoEleXiMassWSSide(0),
  fHistoEleXiMassRS1(0),
  fHistoEleXiMassWS1(0),
  fHistoEleXiMassRSMix1(0),
  fHistoEleXiMassWSMix1(0),
  fHistoEleXiMassRSSide1(0),
  fHistoEleXiMassWSSide1(0),
  fHistoEleXiMassRS2(0),
  fHistoEleXiMassWS2(0),
  fHistoEleXiMassRSMix2(0),
  fHistoEleXiMassWSMix2(0),
  fHistoEleXiMassRSSide2(0),
  fHistoEleXiMassWSSide2(0),
  fHistoEleXiMassAway(0),
  fHistoEleXiMassRSAway(0),
  fHistoEleXiMassWSAway(0),
  fHistoEleXiMassRSMixAway(0),
  fHistoEleXiMassWSMixAway(0),
  fHistoEleXiMassRSSideAway(0),
  fHistoEleXiMassWSSideAway(0),
  fHistoEleXiMassRS1Away(0),
  fHistoEleXiMassWS1Away(0),
  fHistoEleXiMassRSMix1Away(0),
  fHistoEleXiMassWSMix1Away(0),
  fHistoEleXiMassRSSide1Away(0),
  fHistoEleXiMassWSSide1Away(0),
  fHistoEleXiMassRS2Away(0),
  fHistoEleXiMassWS2Away(0),
  fHistoEleXiMassRSMix2Away(0),
  fHistoEleXiMassWSMix2Away(0),
  fHistoEleXiMassRSSide2Away(0),
  fHistoEleXiMassWSSide2Away(0),
  fHistoEleXiMassvsElePtRS(0),
  fHistoEleXiMassvsElePtWS(0),
  fHistoEleXiMassvsElePtRSMix(0),
  fHistoEleXiMassvsElePtWSMix(0),
  fHistoEleXiMassvsElePtRSSide(0),
  fHistoEleXiMassvsElePtWSSide(0),
  fHistoEleXiMassvsElePtRS1(0),
  fHistoEleXiMassvsElePtWS1(0),
  fHistoEleXiMassvsElePtRSMix1(0),
  fHistoEleXiMassvsElePtWSMix1(0),
  fHistoEleXiMassvsElePtRSSide1(0),
  fHistoEleXiMassvsElePtWSSide1(0),
  fHistoEleXiMassvsElePtRS2(0),
  fHistoEleXiMassvsElePtWS2(0),
  fHistoEleXiMassvsElePtRSMix2(0),
  fHistoEleXiMassvsElePtWSMix2(0),
  fHistoEleXiMassvsElePtRSSide2(0),
  fHistoEleXiMassvsElePtWSSide2(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleXiMassMCS(0),
  fHistoEleXiMassMCS1(0),
  fHistoEleXiMassMCS2(0),
  fHistoEleXiMassXibMCS(0),
  fHistoEleXiMassXibMCS1(0),
  fHistoEleXiMassXibMCS2(0),
  fHistoEleXiMassPromptMCS(0),
  fHistoEleXiMassPromptMCS1(0),
  fHistoEleXiMassPromptMCS2(0),
  fHistoEleXiMassBFeeddownMCS(0),
  fHistoEleXiMassBFeeddownMCS1(0),
  fHistoEleXiMassBFeeddownMCS2(0),
  fHistoEleXiMassMCGen(0),
  fHistoEleXiMassvsElePtMCS(0),
  fHistoEleXiMassvsElePtMCGen(0),
  fHistoEleXiMassvsElePtMCS1(0),
  fHistoEleXiMassvsElePtMCGen1(0),
  fHistoEleXiMassvsElePtMCS2(0),
  fHistoEleXiMassvsElePtMCGen2(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsXiPtRS(0),
  fHistoElePtvsXiPtWS(0),
  fHistoElePtvsXiPtRSMix(0),
  fHistoElePtvsXiPtWSMix(0),
  fHistoElePtvsXiPtMCS(0),
  fHistoElePtvsXiPtvsXicPtMCS(0),
  fHistoElePtvsXiPtMCGen(0),
  fHistoElePtvsXiPtvsXicPtMCGen(0),
  fHistoElePtvsXiPtMCXicGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistod0Bach(0),
  fHistoXiMassvsPt(0),
  fHistoXiMassvsPtMCS(0),
  fHistoXiMassvsPtMCGen(0),
  fHistoOmegaMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTOFPIDSelTPC(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
  fHistoMassConversionsMin(0),
  fHistoMassConversionsSameSignMin(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoXiQovPtvsPhi(0),
	fHistoXicMCGen(0),
	fHistoXicMCGen1(0),
	fHistoXicMCGen2(0),
	fHistoXicMCS(0),
	fHistoXicMCS1(0),
	fHistoXicMCS2(0),
	fHistoXibMCGen(0),
	fHistoXibMCGenWithXic(0),
	fHistoXibMCS(0),
	fHistoXicElectronMCGen(0),
	fHistoXicElectronMCGen1(0),
	fHistoXicElectronMCGen2(0),
	fHistoXicElectronMCS(0),
	fHistoXicElectronMCS1(0),
	fHistoXicElectronMCS2(0),
	fHistoElectronMCGen(0),
	fHistoBottomElectronMCGen(0),
	fHistoCharmElectronMCGen(0),
	fHistoXiMCGen(0),
	fHistoLambdaPtvsDl(0),
	fHistoLambdaPtvsDlSide(0),
	fHistoLambdaPtvsDlMCS(0),
	fHistoLambdaPtvsDR(0),
	fHistoLambdaPtvsDRSide(0),
	fHistoLambdaPtvsDRMCS(0),
	fHistoEleXiPtvsRapidityRS(0),
	fHistoEleXiPtvsRapidityWS(0),
	fHistoEleXiPtvsRapidityMCS(0),
	fHistoCorrelationVariablesvsEleXiPt(0),
	fHistoCorrelationVariablesvsEleXiPtMix(0),
	fHistoCorrelationVariablesvsEleXiPtMC(0),
	fHistoCorrelationVariablesvsElePt(0),
	fHistoCorrelationVariablesvsElePtMix(0),
	fHistoCorrelationVariablesvsElePtMC(0),
	fHistoCorrelationVariablesvsXiPt(0),
	fHistoCorrelationVariablesvsXiPtMix(0),
	fHistoCorrelationVariablesvsXiPtMC(0),
	fHistoMassVariablesvsEleXiPt(0),
	fHistoMassVariablesvsEleXiPtMix(0),
	fHistoMassVariablesvsEleXiPtMC(0),
	fHistoMassVariablesvsElePt(0),
	fHistoMassVariablesvsElePtMix(0),
	fHistoMassVariablesvsElePtMC(0),
	fHistoMassVariablesvsXiPt(0),
	fHistoMassVariablesvsXiPtMix(0),
	fHistoMassVariablesvsXiPtMC(0),
	fHistoResponseElePt(0),
	fHistoResponseXiPt(0),
	fHistoResponseEleXiPt(0),
	fHistoResponseXiPtvsEleXiPt(0),
	fHistoResponseXiPtXib(0),
	fHistoResponseEleXiPtXib(0),
	fHistoResponseMcGenXibPtvsXicPt(0),
  fHistoPi0MCGen(0),
  fHistoElectronPi0Total(0),
  fHistoElectronPi0Tag(0),
  fHistoEtaMCGen(0),
  fHistoElectronEtaTotal(0),
  fHistoElectronEtaTag(0),
  fHistoKaonMCGen(0),
  fHistoD0MCGen(0),
	fCounter(0),
	fHistonEvtvsRunNumber(0),
	fHistonElevsRunNumber(0),
	fHistonXivsRunNumber(0),
	fHistoMCEventType(0),
	fHistoMCXic0Decays(0),
	fHistoMCDeltaPhiccbar(0),
	fHistoMCNccbar(0),
	fRefMult(9.26),
  fGTI(0),fGTIndex(0), fTrackBuffSize(19000),
  fHistodPhiSdEtaSElectronProtonR125RS(0),
  fHistodPhiSdEtaSElectronProtonR125WS(0),
  fHistodPhiSdEtaSElectronProtonR125RSMix(0),
  fHistodPhiSdEtaSElectronProtonR125WSMix(0),
  fHistodPhiSdEtaSElectronPionR125RS(0),
  fHistodPhiSdEtaSElectronPionR125WS(0),
  fHistodPhiSdEtaSElectronPionR125RSMix(0),
  fHistodPhiSdEtaSElectronPionR125WSMix(0),
  fHistodPhiSdEtaSElectronBachelorR125RS(0),
  fHistodPhiSdEtaSElectronBachelorR125WS(0),
  fHistodPhiSdEtaSElectronBachelorR125RSMix(0),
  fHistodPhiSdEtaSElectronBachelorR125WSMix(0),
  fDoEventMixing(0),
  fMixWithoutConversionFlag(kFALSE),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
  fNRPBins					(0), 
  fNOfPools(1),
  fPoolIndex(-9999),
  nextResVec(),
  reservoirsReady(),
  m_ReservoirE(),
  m_ReservoirL1(),
  m_ReservoirL2(),
  m_ReservoirVarsE(),
  m_ReservoirVarsL1(),
  m_ReservoirVarsL2()
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEXic2eleXifromAODtracks","Calling Constructor");

	for(Int_t i=0;i<23;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
	}
	for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;

  DefineOutput(1,TList::Class());  //conters
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());  //conters
  DefineOutput(4,TTree::Class());  //My private output
  DefineOutput(5,TTree::Class());  //My private output
  DefineOutput(6,TTree::Class());  //My private output
  DefineOutput(7,TTree::Class());  //My private output
  DefineOutput(8,AliNormalizationCounter::Class());
  DefineOutput(9,TTree::Class());  //My private output
  DefineOutput(10,TTree::Class());  //My private output
  DefineOutput(11,TTree::Class());  //My private output
  DefineOutput(12,TTree::Class());  //My private output
}

//___________________________________________________________________________
AliAnalysisTaskSEXic2eleXifromAODtracks::~AliAnalysisTaskSEXic2eleXifromAODtracks() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEXic2eleXifromAODtracks","Calling Destructor");

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fOutputAll) {
    delete fOutputAll;
    fOutputAll = 0;
  }

  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }


  if (fAnalCuts) {
    delete fAnalCuts;
    fAnalCuts = 0;
  }

  if (fVariablesTree) {
    delete fVariablesTree;
    fVariablesTree = 0;
  }
  if (fEleVariablesTree) {
    delete fEleVariablesTree;
    fEleVariablesTree = 0;
  }
  if (fCascVariablesTree) {
    delete fCascVariablesTree;
    fCascVariablesTree = 0;
  }
  if (fSingleVariablesTree) {
    delete fSingleVariablesTree;
    fSingleVariablesTree = 0;
  }
  if (fMCVariablesTree) {
    delete fMCVariablesTree;
    fMCVariablesTree = 0;
  }
  if (fMCEleVariablesTree) {
    delete fMCEleVariablesTree;
    fMCEleVariablesTree = 0;
  }
  if (fMCCascVariablesTree) {
    delete fMCCascVariablesTree;
    fMCCascVariablesTree = 0;
  }
  if (fMCGenPairVariablesTree) {
    delete fMCGenPairVariablesTree;
    fMCGenPairVariablesTree = 0;
  }
  if (fCorrelationVariablesTree) {
    delete fCorrelationVariablesTree;
    fCorrelationVariablesTree = 0;
  }
	if(fCounter){
		delete fCounter;
		fCounter = 0;
	}

  for(Int_t i = 0;i<fNOfPools;i++){
    for(Int_t j=0;j<fNumberOfEventsForMixing;j++){
      while(!m_ReservoirE[i][j].empty()){
        delete m_ReservoirE[i][j].back();
        m_ReservoirE[i][j].pop_back();
      }
      while(!m_ReservoirL1[i][j].empty()){
        delete m_ReservoirL1[i][j].back();
        m_ReservoirL1[i][j].pop_back();
      }
      while(!m_ReservoirL2[i][j].empty()){
        delete m_ReservoirL2[i][j].back();
        m_ReservoirL2[i][j].pop_back();
      }
      while(!m_ReservoirVarsE[i][j].empty()){
        delete m_ReservoirVarsE[i][j].back();
        m_ReservoirVarsE[i][j].pop_back();
      }
      while(!m_ReservoirVarsL1[i][j].empty()){
        delete m_ReservoirVarsL1[i][j].back();
        m_ReservoirVarsL1[i][j].pop_back();
      }
      while(!m_ReservoirVarsL2[i][j].empty()){
        delete m_ReservoirVarsL2[i][j].back();
        m_ReservoirVarsL2[i][j].pop_back();
      }
    }
  }

  if (fGTI)
    delete[] fGTI;
  fGTI=0;
  if (fGTIndex)
    delete[] fGTIndex;
  fGTIndex=0;
}

//_________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::Init() {
  //
  // Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsXictoeleXifromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec(Option_t *)
{
  //
  // UserExec
  //

  if (!fInputEvent) {
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fCEvents->Fill(1);
	fEvNumberCounter++;

  //------------------------------------------------
  // First check if the event has proper B
  //------------------------------------------------
  fBzkG = (Double_t)aodEvent->GetMagneticField(); 
  AliKFParticle::SetField(fBzkG);
  if (TMath::Abs(fBzkG)<0.001) {
    return;
  }
  fCEvents->Fill(2);

	Int_t countTr=0;
	Double_t countCorr=0;
  if(fUseCentralitySPDTracklet)
  {
    AliAODTracklets* tracklets=aodEvent->GetTracklets();
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>-1.0 && eta<1.0) countTr++;
    }
    AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
    Bool_t isVtxOk=kFALSE;
    if(vtx1){
      if(vtx1->GetNContributors()>0){
        fCEvents->Fill(8);
        isVtxOk=kTRUE;
      }
    }

    countCorr=countTr;
    if(isVtxOk){
      TProfile* estimatorAvg = GetEstimatorHistogram(aodEvent);
      countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countTr,vtx1->GetZ(),fRefMult));
    }
  }


  if(fUseCentralitySPDTracklet){
    fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo,countCorr);
  }else{
    fCounter->StoreEvent(aodEvent,fAnalCuts,fUseMCInfo);
  }

  fIsEventSelected = fAnalCuts->IsEventSelected(aodEvent);

  //------------------------------------------------
  // MC analysis setting
  //------------------------------------------------
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader=0;
  if (fUseMCInfo) {
    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fCEvents->Fill(6); // in case of MC events
  
    // load MC header
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      AliError("AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec: MC header branch not found!\n");
      return;
    }
    fCEvents->Fill(7); // in case of MC events
  
    Double_t zMCVertex = mcHeader->GetVtxZ();
    if (TMath::Abs(zMCVertex) > fAnalCuts->GetMaxVtxZ()) {
      AliDebug(2,Form("Event rejected: abs(zVtxMC)=%f > fAnalCuts->GetMaxVtxZ()=%f",zMCVertex,fAnalCuts->GetMaxVtxZ()));
      return;
    } else {
      fCEvents->Fill(17); // in case of MC events
    }
    if ((TMath::Abs(zMCVertex) < fAnalCuts->GetMaxVtxZ()) && (!fAnalCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnalCuts->IsEventRejectedDueToTrigger())) {
			Bool_t selevt = MakeMCAnalysis(mcArray);
			if(!selevt) return;
		}
  }

  //------------------------------------------------
  // Event selection 
  //------------------------------------------------
  fVtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if (!fVtx1) return;

  Double_t pos[3],cov[6];
  fVtx1->GetXYZ(pos);
  fVtx1->GetCovarianceMatrix(cov);
  fV1 = new AliESDVertex(pos,cov,100.,100,fVtx1->GetName());
	fVtxZ = pos[2];

  Bool_t fIsTriggerNotOK = fAnalCuts->IsEventRejectedDueToTrigger();
  Bool_t fIsPhysSelNotOK = fAnalCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t fIsNoVertex = fAnalCuts->IsEventRejectedDueToNotRecoVertex();
  if(!fIsTriggerNotOK && !fIsPhysSelNotOK && !fIsNoVertex && fabs(fVtx1->GetZ())<fAnalCuts->GetMaxVtxZ()) fCEvents->Fill(3);
  if(!fIsEventSelected) {
    delete fV1;
    return;
  }
  fCEvents->Fill(4);

	fHNTrackletvsZ->Fill(fVtxZ,countTr);
	fHNTrackletCorrvsZ->Fill(fVtxZ,countCorr);

  fIsMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  fIsSemi=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  fIsCent=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral); 
  fIsINT7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);  
  fIsEMC7=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);   
  fTriggerCheck = fIsMB+2*fIsSemi+4*fIsCent+8*fIsINT7+16*fIsEMC7;
  if(fIsMB) fHTrigger->Fill(1);
  if(fIsSemi) fHTrigger->Fill(2);
  if(fIsCent) fHTrigger->Fill(3);
  if(fIsINT7) fHTrigger->Fill(4);
  if(fIsEMC7) fHTrigger->Fill(5);
  if(fIsMB|fIsSemi|fIsCent) fHTrigger->Fill(7);
  if(fIsINT7|fIsEMC7) fHTrigger->Fill(8);
  if(fIsMB&fIsSemi) fHTrigger->Fill(10);
  if(fIsMB&fIsCent) fHTrigger->Fill(11);
  if(fIsINT7&fIsEMC7) fHTrigger->Fill(12);

	if(fUseCentralityV0M){
		AliCentrality *cent = aodEvent->GetCentrality();
		fCentrality = cent->GetCentralityPercentile("V0M");
	}else if(fUseCentralitySPDTracklet){
		if(countCorr>=0 && countCorr<=0) fCentrality = 5.;
    else if(countCorr>=1 && countCorr<=8) fCentrality = 15.;
		else if(countCorr>=9 && countCorr<=13) fCentrality = 25.;
		else if(countCorr>=14 && countCorr<=19) fCentrality = 35.;
		else if(countCorr>=20 && countCorr<=30) fCentrality = 45.;
		else if(countCorr>=31 && countCorr<=49) fCentrality = 55.;
		else fCentrality = 65.;
	}else{
		fCentrality = 1.;
	}
  if(fCentrality<0.||fCentrality>100.-0.0000001) {
    delete fV1;
    return;
  }
  fHCentrality->Fill(fCentrality);

  if(fUseEventPlane>0){
    AliEventplane *pl=aodEvent->GetEventplane();
    if(!pl){
      AliError("AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
      fCEvents->Fill(18);
      return;
    }
    Double_t ep_v0m = GetPhi0Pi(pl->GetEventplane("V0",aodEvent,2));
    Double_t ep_v0a = GetPhi0Pi(pl->GetEventplane("V0A",aodEvent,2));
    Double_t ep_v0c = GetPhi0Pi(pl->GetEventplane("V0C",aodEvent,2));
    Double_t ep_tpc = GetPhi0Pi(pl->GetEventplane("Q"));
    if(fUseEventPlane==1)
      fEventPlane = ep_v0m;
    if(fUseEventPlane==2)
      fEventPlane = ep_v0a;
    if(fUseEventPlane==3)
      fEventPlane = ep_v0c;
    if(fUseEventPlane==4){
      fEventPlane = ep_tpc;
      fQSub1 = pl->GetQsub1();
      fQSub2 = pl->GetQsub2();
      fQ = pl->GetQVector();
			if(!fQSub1 || !fQSub2){
				AliError("AliAnalysisTaskSEXic2eleXifromAODtracks::UserExec: no Q vectors");
				fCEvents->Fill(19);
				return;
			}
		}
  }

  fHEventPlane->Fill(fEventPlane);
	fRunNumber = aodEvent->GetRunNumber();

	Int_t runnumber_offset = 0;
	Int_t runnumber = aodEvent->GetRunNumber();
	if(runnumber<=131000&&runnumber>=114000){
		runnumber_offset = 114000;//lhc10bcde
	}else if(runnumber<=196000&&runnumber>=195000){
		runnumber_offset = 195000;//lhc13bc
	}else if(runnumber<=170593&&runnumber>=167902){
		runnumber_offset = 167902;//lhc11h
	}else if(runnumber<=246994&&runnumber>=244824){
		runnumber_offset = 244824;//lhc15o
	}
	fHistonEvtvsRunNumber->Fill(runnumber-runnumber_offset,1.);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  //Int_t nv0 = aodEvent->GetNumberOfV0s();
  fCEvents->Fill(5);


  //------------------------------------------------
  // Main analysis done in this function
  //------------------------------------------------
  fAnalCuts->SetMagneticField(fBzkG);
  fAnalCuts->SetPrimaryVertex(pos);
  MakeAnalysis(aodEvent,mcArray);

  PostData(1,fOutput);
  PostData(3,fOutputAll);
  PostData(4,fVariablesTree);
  PostData(5,fEleVariablesTree);
  PostData(6,fCascVariablesTree);
  PostData(7,fMCVariablesTree);
  PostData(8,fCounter);    
  PostData(9,fMCEleVariablesTree);
  PostData(10,fMCCascVariablesTree);
  //PostData(11,fMCGenPairVariablesTree);
  PostData(11,fSingleVariablesTree);
  PostData(12,fCorrelationVariablesTree);

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  //AliInfo("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    AliError("fOutput not available");
    return;
  }

  fOutputAll = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputAll) {     
    AliError("fOutputAll not available");
    return;
  }

  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::UserCreateOutputObjects() 
{ 
  //
  // UserCreateOutputObject
  //
  //AliInfo(Form("CreateOutputObjects of task %s\n", GetName()));

  //------------------------------------------------
  // output object setting
  //------------------------------------------------
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");
  DefineGeneralHistograms(); // define general histograms
  PostData(1,fOutput);

  fOutputAll = new TList();
  fOutputAll->SetOwner();
  fOutputAll->SetName("anahisto");
  DefineAnalysisHistograms(); // define general histograms
  PostData(3,fOutputAll);

  DefineTreeVariables();
  PostData(4,fVariablesTree);

  DefineEleTreeVariables();
  PostData(5,fEleVariablesTree);

  DefineCascTreeVariables();
  PostData(6,fCascVariablesTree);

  DefineMCTreeVariables();
  PostData(7,fMCVariablesTree);

  DefineMCEleTreeVariables();
  PostData(9,fMCEleVariablesTree);

  DefineMCCascTreeVariables();
  PostData(10,fMCCascVariablesTree);

  //DefineMCGenPairTreeVariables();
  //PostData(11,fMCGenPairVariablesTree);
  DefineSingleTreeVariables();
  PostData(11,fSingleVariablesTree);

  DefineCorrelationTreeVariables();
  PostData(12,fCorrelationVariablesTree);


  //Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(8)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  if(fUseCentralitySPDTracklet){
    fCounter->SetStudyMultiplicity(kTRUE,1.);
  }
  fCounter->Init();
  PostData(8,fCounter);

	if(fDoEventMixing){
		fNOfPools=fNCentBins*fNzVtxBins*fNRPBins;
    m_ReservoirE.resize(fNOfPools,std::vector<std::vector<TLorentzVector *> > (fNumberOfEventsForMixing));
    m_ReservoirL1.resize(fNOfPools,std::vector<std::vector<TLorentzVector *> > (fNumberOfEventsForMixing));
    m_ReservoirL2.resize(fNOfPools,std::vector<std::vector<TLorentzVector *> > (fNumberOfEventsForMixing));
    m_ReservoirVarsE.resize(fNOfPools,std::vector<std::vector<TVector *>  > (fNumberOfEventsForMixing));
    m_ReservoirVarsL1.resize(fNOfPools,std::vector<std::vector<TVector *>  > (fNumberOfEventsForMixing));
    m_ReservoirVarsL2.resize(fNOfPools,std::vector<std::vector<TVector *>  > (fNumberOfEventsForMixing));
    nextResVec.resize(fNOfPools,0);
    reservoirsReady.resize(fNOfPools,kFALSE);

    for(Int_t s=0; s<fNOfPools; s++) {
      for(Int_t k=0;k<fNumberOfEventsForMixing;k++){
        m_ReservoirE[s][k].clear();
        m_ReservoirL1[s][k].clear();
        m_ReservoirL2[s][k].clear();
        m_ReservoirVarsE[s][k].clear();
        m_ReservoirVarsL1[s][k].clear();
        m_ReservoirVarsL2[s][k].clear();
      }
    }
	}

  fGTI = new AliAODTrack *[fTrackBuffSize]; // Array of pointers 
  fGTIndex = new Int_t [fTrackBuffSize]; // Array of index 

  return;
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  // Main Analysis part
  //
  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------
  if(fDoEventMixing){
    fPoolIndex=GetPoolIndex(fVtxZ,fCentrality,fEventPlane);
    Int_t nextRes( nextResVec[fPoolIndex] );
    while(!m_ReservoirE[fPoolIndex][nextRes].empty()){
      delete m_ReservoirE[fPoolIndex][nextRes].back();
      m_ReservoirE[fPoolIndex][nextRes].pop_back();
    }
    while(!m_ReservoirL1[fPoolIndex][nextRes].empty()){
      delete m_ReservoirL1[fPoolIndex][nextRes].back();
      m_ReservoirL1[fPoolIndex][nextRes].pop_back();
    }
    while(!m_ReservoirL2[fPoolIndex][nextRes].empty()){
      delete m_ReservoirL2[fPoolIndex][nextRes].back();
      m_ReservoirL2[fPoolIndex][nextRes].pop_back();
    }
    while(!m_ReservoirVarsE[fPoolIndex][nextRes].empty()){
      delete m_ReservoirVarsE[fPoolIndex][nextRes].back();
      m_ReservoirVarsE[fPoolIndex][nextRes].pop_back();
    }
    while(!m_ReservoirVarsL1[fPoolIndex][nextRes].empty()){
      delete m_ReservoirVarsL1[fPoolIndex][nextRes].back();
      m_ReservoirVarsL1[fPoolIndex][nextRes].pop_back();
    }
    while(!m_ReservoirVarsL2[fPoolIndex][nextRes].empty()){
      delete m_ReservoirVarsL2[fPoolIndex][nextRes].back();
      m_ReservoirVarsL2[fPoolIndex][nextRes].pop_back();
    }
  }

  if(fAnalCuts->GetProdAODFilterBit()==7){
    ResetGlobalTrackReference();
    // ..and set it
    for (Int_t iTrack=0;iTrack<aodEvent->GetNumberOfTracks();iTrack++){
      // cast needed since the event now returns AliVTrack instead of AliAODTrack
      AliAODTrack *track = dynamic_cast<AliAODTrack *>(aodEvent->GetTrack(iTrack));
      if (!track) continue;

      // Store the reference of the global tracks
      StoreGlobalTrackReference(track,iTrack);
    }
  }


  Int_t nCascs= aodEvent->GetNumberOfCascades();
  Int_t nTracks= aodEvent->GetNumberOfTracks();

  Bool_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Bool_t  seleCascFlags[nCascs];
  Int_t     nSeleCasc=0;
  SelectCascade(aodEvent,nCascs,nSeleCasc,seleCascFlags,mcArray);

	Int_t runnumber_offset = 0;
	Int_t runnumber = aodEvent->GetRunNumber();
	if(runnumber<=131000&&runnumber>=114000){
		runnumber_offset = 114000;//lhc10bcde
	}else if(runnumber<=196000&&runnumber>=195000){
		runnumber_offset = 195000;//lhc13bc
	}else if(runnumber<=170593&&runnumber>=167902){
		runnumber_offset = 167902;//lhc11h
	}else if(runnumber<=246994&&runnumber>=244824){
		runnumber_offset = 244824;//lhc15o
	}
	fHistonElevsRunNumber->Fill(runnumber-runnumber_offset,nSeleTrks);
	fHistonXivsRunNumber->Fill(runnumber-runnumber_offset,nSeleCasc);

  if(nSeleTrks==0 || nSeleCasc==0) return;

  //------------------------------------------------
  // Fill pool and single tree
  //------------------------------------------------
  for (Int_t itrk = 0; itrk<nTracks; itrk++) {
    if(!seleTrkFlags[itrk]) continue;
    AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);

    FillElectronROOTObjects(trk,aodEvent,mcArray);
  }

  for (Int_t icasc = 0; icasc<nCascs; icasc++) {
    if(!seleCascFlags[icasc]) continue;
    AliAODcascade *casc = aodEvent->GetCascade(icasc);
    if(!casc) continue;

		FillCascROOTObjects(casc,aodEvent,mcArray);
  }

  if(fWriteEachVariableTree)
    return;

  //------------------------------------------------
  // Cascade loop 
  //------------------------------------------------
  for (Int_t icasc = 0; icasc<nCascs; icasc++) {
    if(!seleCascFlags[icasc]) continue;
    AliAODcascade *casc = aodEvent->GetCascade(icasc);
    if(!casc) continue;

    AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
    AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));

    //------------------------------------------------
    // track loop 
    //------------------------------------------------
    for (Int_t itrk = 0; itrk<nTracks; itrk++) {
      if(!seleTrkFlags[itrk]) continue;
      AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);
      //if(trk->GetID()<0) continue;

      //if(!fAnalCuts->SelectWithRoughCuts(v0,trk)) continue;
      
      //TPC only track (BIT 7) does not have PID information 
      //In addition to that, TPC only tracks does not have good DCA resolution
      //(according to femtoscopy code)
      AliAODTrack *trkpid = 0;
      if(fAnalCuts->GetProdAODFilterBit()==7){
        trkpid = fGTI[-trk->GetID()-1];
      }else{
        trkpid = trk;
      }

      Int_t cpid = cptrack->GetID();
      Int_t cnid = cntrack->GetID();
      Int_t cbid = cbtrack->GetID();
      Int_t lpid = trkpid->GetID();
      if((cpid==lpid)||(cnid==lpid)||(cbid==lpid)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(casc,trk,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
      if(!secVert) continue;

      AliAODRecoCascadeHF *exobj = MakeCascadeHF(casc,trk,trkpid,aodEvent,secVert);
      if(!exobj) {
	continue;
      }

      FillROOTObjects(exobj, casc,trk,trkpid,aodEvent,mcArray);

      exobj->GetSecondaryVtx()->RemoveDaughters();
      exobj->UnsetOwnPrimaryVtx();
      delete exobj;exobj=NULL;
      delete secVert;
    }
  }

  if(fDoEventMixing){
    DoEventMixingWithPools(fPoolIndex);

    Int_t nextRes( nextResVec[fPoolIndex] );
    nextRes++;
    if( nextRes>=fNumberOfEventsForMixing ){
      nextRes = 0;
      reservoirsReady[fPoolIndex] = kTRUE;
    }
    nextResVec[fPoolIndex] = nextRes;
  }
}


////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineTreeVariables() 
{
  //
  // Define tree variables
  //

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 92;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassEleXi";
  fCandidateVariableNames[ 2]="EleXiPt";
  fCandidateVariableNames[ 3]="EleXiPx";
  fCandidateVariableNames[ 4]="EleXiPy";
  fCandidateVariableNames[ 5]="EleXiPz";
  fCandidateVariableNames[ 6]="ElePx";
  fCandidateVariableNames[ 7]="ElePy";
  fCandidateVariableNames[ 8]="ElePz";
  fCandidateVariableNames[ 9]="XiPx";
  fCandidateVariableNames[10]="XiPy";
  fCandidateVariableNames[11]="XiPz";
  fCandidateVariableNames[12]="XiCharge";
  fCandidateVariableNames[13]="MassXi";
  fCandidateVariableNames[14]="MassLambda";
  fCandidateVariableNames[15]="Eled0";
  fCandidateVariableNames[16]="Xid0";
  fCandidateVariableNames[17]="nSigmaTPCele";
  fCandidateVariableNames[18]="nSigmaTOFele";
  fCandidateVariableNames[19]="nSigmaTPCpr_etrk";
  fCandidateVariableNames[20]="nSigmaTOFpr_etrk";
  fCandidateVariableNames[21]="nSigmaTPCka_etrk";
  fCandidateVariableNames[22]="nSigmaTOFka_etrk";
  fCandidateVariableNames[23]="nSigmaTPCv0pr";
  fCandidateVariableNames[24]="nSigmaTOFv0pr";
  fCandidateVariableNames[25]="nSigmaTPCv0pi";
  fCandidateVariableNames[26]="nSigmaTOFv0pi";
  fCandidateVariableNames[27]="nSigmaTPCbachpi";
  fCandidateVariableNames[28]="nSigmaTOFbachpi";
  fCandidateVariableNames[29]="EleCharge";
  fCandidateVariableNames[30]="Mixing";
  fCandidateVariableNames[31]="DcaXiDaughters";
  fCandidateVariableNames[32]="DcaV0Daughters";
  fCandidateVariableNames[33]="DecayLengthXi";
  fCandidateVariableNames[34]="CosPointingAngleXi";
  fCandidateVariableNames[35]="DcaV0toPrimVertex";
  fCandidateVariableNames[36]="DcaPostoPrimVertex";
  fCandidateVariableNames[37]="DcaNegtoPrimVertex";
  fCandidateVariableNames[38]="DcaBachtoPrimVertex";
  fCandidateVariableNames[39]="DecayLengthV0";
  fCandidateVariableNames[40]="CosPointingAngleV0";

  fCandidateVariableNames[41]="mcpdgxic";
  fCandidateVariableNames[42]="mclabxic";
  fCandidateVariableNames[43]="mcxicpx";
  fCandidateVariableNames[44]="mcxicpy";
  fCandidateVariableNames[45]="mcxicpz";
  fCandidateVariableNames[46]="mcelepx";
  fCandidateVariableNames[47]="mcelepy";
  fCandidateVariableNames[48]="mcelepz";
  fCandidateVariableNames[49]="mccascpx";
  fCandidateVariableNames[50]="mccascpy";
  fCandidateVariableNames[51]="mccascpz";

  fCandidateVariableNames[52]="mcpdgele";
  fCandidateVariableNames[53]="mcpdgcasc";
  fCandidateVariableNames[54]="mcpdgmomele";
  fCandidateVariableNames[55]="mcpdgmomcasc";
  fCandidateVariableNames[56]="mcpdggrmomele";
  fCandidateVariableNames[57]="mcpdggrmomcasc";
  fCandidateVariableNames[58]="mcngenele";
  fCandidateVariableNames[59]="mcngencasc";

  fCandidateVariableNames[60]="nSigmaTPCpi_etrk";
  fCandidateVariableNames[61]="nSigmaTOFpi_etrk";

  fCandidateVariableNames[62]="V0PosPx";
  fCandidateVariableNames[63]="V0PosPy";
  fCandidateVariableNames[64]="V0PosPz";
  fCandidateVariableNames[65]="V0NegPx";
  fCandidateVariableNames[66]="V0NegPy";
  fCandidateVariableNames[67]="V0NegPz";
  fCandidateVariableNames[68]="V0VertX";
  fCandidateVariableNames[69]="V0VertY";
  fCandidateVariableNames[70]="V0VertZ";
  fCandidateVariableNames[71]="BachPx";
  fCandidateVariableNames[72]="BachPy";
  fCandidateVariableNames[73]="BachPz";
  fCandidateVariableNames[74]="XiVertX";
  fCandidateVariableNames[75]="XiVertY";
  fCandidateVariableNames[76]="XiVertZ";
  fCandidateVariableNames[77]="PrimVertX";
  fCandidateVariableNames[78]="PrimVertY";
  fCandidateVariableNames[79]="PrimVertZ";

  fCandidateVariableNames[80]="MassOmega";

	fCandidateVariableNames[81]= "EleITSMatch";
	fCandidateVariableNames[82]= "BachITSMatch";
	fCandidateVariableNames[83]= "V0PosITSMatch";
	fCandidateVariableNames[84]= "V0NegITSMatch";

	fCandidateVariableNames[85]= "TPCNclsF";
	fCandidateVariableNames[86]= "TPCNcls";
	fCandidateVariableNames[87]= "TPCNclsS";
	fCandidateVariableNames[88]= "IsXiPeakReagion";


  fCandidateVariableNames[89]="MagneticField";
  fCandidateVariableNames[90]="EvNumber";
  fCandidateVariableNames[91]="RunNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *exobj, AliAODcascade *casc, AliAODTrack *trk, AliAODTrack *trkpid, AliAODEvent *event, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!trk) return;
	if(!trkpid) return;
	if(!casc) return;

	for(Int_t i=0;i<92;i++){
		fCandidateVariables[i] = -9999.;
	}


  AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
  AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
//	if(cptrack->Charge()<0 && cntrack->Charge()>0){
//		cptrack =   (AliAODTrack*)(casc->GetDaughter(1));
//		cntrack =   (AliAODTrack*)(casc->GetDaughter(0));
//	}
//  Double_t d0z0[2],covd0z0[3];
//  cptrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging
//  cntrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging
//  cbtrack->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);//propagate to primary vertex for debugging

	Double_t pxe = exobj->PxProng(0);
	Double_t pye = exobj->PyProng(0);
	Double_t pze = exobj->PzProng(0);
	Double_t mome = sqrt(pxe*pxe+pye*pye+pze*pze);
	Double_t Ee = sqrt(mome*mome+0.000510998928*0.000510998928);
	Double_t pxv = exobj->PxProng(1);
	Double_t pyv = exobj->PyProng(1);
	Double_t pzv = exobj->PzProng(1);
	Double_t momv = sqrt(pxv*pxv+pyv*pyv+pzv*pzv);
	Double_t Ev = sqrt(momv*momv+1.32171*1.32171);
	Double_t cosoa = (pxe*pxv+pye*pyv+pze*pzv)/mome/momv;
	Double_t Esum = Ee + Ev;

	Double_t uxe = pxe/mome;
	Double_t uye = pye/mome;
	Double_t uze = pze/mome;
	Double_t lf = -2.*(pxv*uxe+pyv*uye+pzv*uze);
	Double_t pxv_flip = pxv + lf * uxe;
	Double_t pyv_flip = pyv + lf * uye;
	Double_t pzv_flip = pzv + lf * uze;
	Double_t pxsum_flip = pxe + pxv_flip;
	Double_t pysum_flip = pye + pyv_flip;
	Double_t pzsum_flip = pze + pzv_flip;
	Double_t mexi_flip = sqrt(Esum*Esum-pxsum_flip*pxsum_flip-pysum_flip*pysum_flip-pzsum_flip*pzsum_flip);
	Double_t ptexi_flip = sqrt(pxsum_flip*pxsum_flip+pysum_flip*pysum_flip);
 
  Double_t minmass_ee = 9999.;
  Bool_t isconv = fAnalCuts->TagConversions(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmass_ee);
  Double_t minmasslike_ee = 9999.;
  Bool_t isconv_like = fAnalCuts->TagConversionsSameSign(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmasslike_ee);


  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);

  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3312};
  fCandidateVariables[ 1] = exobj->InvMass(2,pdgdg);
  fCandidateVariables[ 2] = exobj->Pt();
  fCandidateVariables[ 3] = exobj->Px();
  fCandidateVariables[ 4] = exobj->Py();
  fCandidateVariables[ 5] = exobj->Pz();
  fCandidateVariables[ 6] = exobj->PxProng(0);
  fCandidateVariables[ 7] = exobj->PyProng(0);
  fCandidateVariables[ 8] = exobj->PzProng(0);
  fCandidateVariables[ 9] = exobj->PxProng(1);
  fCandidateVariables[10] = exobj->PyProng(1);
  fCandidateVariables[11] = exobj->PzProng(1);
  fCandidateVariables[12] = casc->ChargeXi();
  fCandidateVariables[13] = casc->MassXi();
	if(casc->ChargeXi()<0)
		fCandidateVariables[14] = casc->MassLambda();
	else
		fCandidateVariables[14] = casc->MassAntiLambda();
  fCandidateVariables[15] = exobj->Getd0Prong(0);
  fCandidateVariables[16] = exobj->Getd0Prong(1);

	Double_t nSigmaTPCele = -9999.;
	Double_t nSigmaTOFele = -9999.;
  if(fAnalCuts->GetIsUsePID())
  {
		nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
		nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);
    fCandidateVariables[17] = nSigmaTPCele;
    fCandidateVariables[18] = nSigmaTOFele;

		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kProton);
		Double_t nSigmaTOFpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kProton);
		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kKaon);
		Double_t nSigmaTOFka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kKaon);
		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kPion);
		Double_t nSigmaTOFpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kPion);
    fCandidateVariables[19] = nSigmaTPCpr_etrk;
    fCandidateVariables[20] = nSigmaTOFpr_etrk;
    fCandidateVariables[21] = nSigmaTPCka_etrk;
    fCandidateVariables[22] = nSigmaTOFka_etrk;
    fCandidateVariables[60] = nSigmaTPCpi_etrk;
    fCandidateVariables[61] = nSigmaTOFpi_etrk;
  }

	Double_t nSigmaTPCv0pr=-9999.;
	Double_t nSigmaTOFv0pr=-9999.;
	Double_t nSigmaTPCv0pi=-9999.;
	Double_t nSigmaTOFv0pi=-9999.;
	Double_t nSigmaTPCbachpi=-9999.;
	Double_t nSigmaTOFbachpi=-9999.;
	if(fAnalCuts->GetUseCascadePID())
	{
		if(casc->ChargeXi()>0){
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
			nSigmaTPCbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kPion);
			nSigmaTOFbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kPion);
		}else{
			nSigmaTPCv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidCascPr()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
			nSigmaTPCbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTPC(cbtrack,AliPID::kPion);
			nSigmaTOFbachpi = fAnalCuts->GetPidCascPi()->GetPidResponse()->NumberOfSigmasTOF(cbtrack,AliPID::kPion);
		}
      fCandidateVariables[23] = nSigmaTPCv0pr;
      fCandidateVariables[24] = nSigmaTOFv0pr;
      fCandidateVariables[25] = nSigmaTPCv0pi;
      fCandidateVariables[26] = nSigmaTOFv0pi;
      fCandidateVariables[27] = nSigmaTPCbachpi;
      fCandidateVariables[28] = nSigmaTOFbachpi;
  }
  fCandidateVariables[29] = trk->Charge();
  fCandidateVariables[30] = 0;
  fCandidateVariables[31] = casc->DcaXiDaughters();
  fCandidateVariables[32] = casc->DcaV0Daughters();
  fCandidateVariables[33] = casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateVariables[34] = casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateVariables[35] = casc->DcaV0ToPrimVertex();
  fCandidateVariables[36] = casc->DcaPosToPrimVertex();
  fCandidateVariables[37] = casc->DcaNegToPrimVertex();
  fCandidateVariables[38] = casc->DcaBachToPrimVertex();
  fCandidateVariables[39] = casc->DecayLengthV0();
  fCandidateVariables[40] = casc->CosPointingAngle(casc->GetDecayVertexXi());

  AliAODMCParticle *mcxic = 0;
  AliAODMCParticle *mcele = 0;
  AliAODMCParticle *mccasc = 0;
  Int_t mclabxic = 0;
	Int_t mcpdgele_array[100];
	Int_t mcpdgcasc_array[100];
	Int_t mclabelele_array[100];
	Int_t mclabelcasc_array[100];
	Int_t mcngen_ele = -9999;
	Int_t mcngen_casc = -9999;

	if(fUseMCInfo && mcArray){
		mclabxic =  MatchToMC(exobj,mcArray,mcpdgele_array, mcpdgcasc_array,mclabelele_array,mclabelcasc_array,mcngen_ele,mcngen_casc);

    if(mclabxic>-1){
      mcxic = (AliAODMCParticle*) mcArray->At(mclabxic);
			if(mclabelele_array[0]>=0)
				mcele = (AliAODMCParticle*) mcArray->At(mclabelele_array[0]);
			if(mclabelcasc_array[0]>=0)
				mccasc = (AliAODMCParticle*) mcArray->At(mclabelcasc_array[0]);
			if(mcxic){
				fCandidateVariables[41] = mcxic->GetPdgCode();
				fCandidateVariables[42] = mcxic->Label();
				fCandidateVariables[43] = mcxic->Px();
				fCandidateVariables[44] = mcxic->Py();
				fCandidateVariables[45] = mcxic->Pz();
			}
			if(mcele){
				fCandidateVariables[46] = mcele->Px();
				fCandidateVariables[47] = mcele->Py();
				fCandidateVariables[48] = mcele->Pz();
			}
			if(mccasc){
				fCandidateVariables[49] = mccasc->Px();
				fCandidateVariables[50] = mccasc->Py();
				fCandidateVariables[51] = mccasc->Pz();
			}
		}
		fCandidateVariables[52] = mcpdgele_array[0];
		fCandidateVariables[53] = mcpdgcasc_array[0];
		fCandidateVariables[54] = mcpdgele_array[1];
		fCandidateVariables[55] = mcpdgcasc_array[1];
		fCandidateVariables[56] = mcpdgele_array[2];
		fCandidateVariables[57] = mcpdgcasc_array[2];
		fCandidateVariables[58] = mcngen_ele;
		fCandidateVariables[59] = mcngen_casc;
	}
	fCandidateVariables[62] = casc->MomPosX();
	fCandidateVariables[63] = casc->MomPosY();
	fCandidateVariables[64] = casc->MomPosZ();
	fCandidateVariables[65] = casc->MomNegX();
	fCandidateVariables[66] = casc->MomNegY();
	fCandidateVariables[67] = casc->MomNegZ();
	fCandidateVariables[68] = casc->DecayVertexV0X();
	fCandidateVariables[69] = casc->DecayVertexV0Y();
	fCandidateVariables[70] = casc->DecayVertexV0Z();
	fCandidateVariables[71] = casc->MomBachX();
	fCandidateVariables[72] = casc->MomBachY();
	fCandidateVariables[73] = casc->MomBachZ();
	fCandidateVariables[74] = casc->DecayVertexXiX();
	fCandidateVariables[75] = casc->DecayVertexXiY();
	fCandidateVariables[76] = casc->DecayVertexXiZ();
	fCandidateVariables[77] = fVtx1->GetX();
	fCandidateVariables[78] = fVtx1->GetY();
	fCandidateVariables[79] = fVtx1->GetZ();

	fCandidateVariables[80] = casc->MassOmega();

	if(trk) fCandidateVariables[81] = trk->GetITSClusterMap();
	if(cbtrack) fCandidateVariables[82] = cbtrack->GetITSClusterMap();
	if(cptrack) fCandidateVariables[83] = cptrack->GetITSClusterMap();
	if(cntrack) fCandidateVariables[84] = cntrack->GetITSClusterMap();

  fCandidateVariables[85] = trk->GetTPCNclsF();
  fCandidateVariables[86] = trk->GetTPCNcls();
  fCandidateVariables[87] = trk->GetTPCnclsS();
  fCandidateVariables[88] = fAnalCuts->IsPeakRegion(casc);

  fCandidateVariables[89] = fBzkG;
  fCandidateVariables[90] = fEvNumberCounter;
  fCandidateVariables[91] = fRunNumber;


//  if(fWriteVariableTree)
//    fVariablesTree->Fill();

  Double_t dphis_ele_pr, detas_ele_pr,dphis_ele_pi, detas_ele_pi, dphis_ele_bach, detas_ele_bach;
  dphis_ele_pr = 9999.;detas_ele_pr = 9999.;dphis_ele_pi = 9999.;detas_ele_pi = 9999.;dphis_ele_bach=9999.;detas_ele_bach=9999.;
  //fAnalCuts->GetdPhiSdEtaSR125(trk,cptrack,cntrack,cbtrack,fBzkG,posVtx, dphis_ele_pr,detas_ele_pr,dphis_ele_pi,detas_ele_pi, dphis_ele_bach, detas_ele_bach);


	Double_t cont[4];
	cont[0] = exobj->InvMass(2,pdgdg);
	cont[1] = exobj->Pt();
	cont[2] = exobj->Getd0Prong(0)*trk->Charge();
	cont[3] = fCentrality;
	fHistoEleXiMass->Fill(cont);

	Double_t cont_flip[4];
	cont_flip[0] = mexi_flip;
	cont_flip[1] = ptexi_flip;
	cont_flip[2] = 0.0;
	cont_flip[3] = fCentrality;

	Double_t cont2[3];
	cont2[0] = exobj->InvMass(2,pdgdg);
	cont2[1] = trk->Pt();
	cont2[2] = fCentrality;

	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = trk->Pt();
	cont_eleptvseta[1] = trk->Eta();
	cont_eleptvseta[2] = fCentrality;

	Double_t cont_eleptvsxipt[3];
	cont_eleptvsxipt[0] = trk->Pt();
	cont_eleptvsxipt[1] = sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
	cont_eleptvsxipt[2] = fCentrality;

	Double_t cont_eleptvsd0[3];
	cont_eleptvsd0[0] = trk->Pt();
	cont_eleptvsd0[1] = exobj->Getd0Prong(0)*trk->Charge();
	cont_eleptvsd0[2] = fCentrality;

	Double_t exobj_mass = exobj->InvMass(2,pdgdg);
	Double_t exobj_px = exobj->Px();
	Double_t exobj_py = exobj->Py();
	Double_t exobj_pz = exobj->Pz();
	Double_t exobj_E = sqrt(exobj_mass*exobj_mass+exobj_px*exobj_px+exobj_py*exobj_py+exobj_pz*exobj_pz);
	Double_t exobj_rap = 0.5*log((exobj_E+exobj_pz)/(exobj_E-exobj_pz));

  //
  // Old strategy only look at mass and pairpt
  //
	if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsPeakRegion(casc))
	{
		if(trk->Charge()*casc->ChargeXi()<0){
			fHistoEleXiMassRS->Fill(cont);
			if(trk->Charge()>0) fHistoEleXiMassRS1->Fill(cont);
			else  fHistoEleXiMassRS2->Fill(cont);
			fHistoEleXiMassvsElePtRS->Fill(cont2);
			if(trk->Charge()>0) fHistoEleXiMassvsElePtRS1->Fill(cont2);
			else  fHistoEleXiMassvsElePtRS2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleXiMassMax()){
				fHistoEleXiPtvsRapidityRS->Fill(exobj->Pt(),exobj_rap);
				fHistoElePtRS->Fill(trk->Pt(),fCentrality);
				fHistoElePtvsEtaRS->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtRS->Fill(cont_eleptvsxipt);
				fHistoElePtvsd0RS->Fill(cont_eleptvsd0);
				for(Int_t ih=0;ih<23;ih++){
					Double_t cont_eleptvscutvars[3];
					cont_eleptvscutvars[0] = exobj->Pt();
					cont_eleptvscutvars[2] = fCentrality;

					if(ih==0){
						cont_eleptvscutvars[1] = trk->GetTPCNcls();
					}else if(ih==1){
						cont_eleptvscutvars[1] = trk->GetTPCsignalN();
					}else if(ih==2){
						cont_eleptvscutvars[1] = nSigmaTPCele;
					}else if(ih==3){
						cont_eleptvscutvars[1] = nSigmaTOFele;
					}else if(ih==4){
						cont_eleptvscutvars[1] = trk->Eta();
					}else if(ih==5){
						cont_eleptvscutvars[1] = trk->GetITSNcls();
					}else if(ih==6){
						if(casc->ChargeXi()<0)
							cont_eleptvscutvars[1] = casc->MassLambda();
						else
							cont_eleptvscutvars[1] = casc->MassAntiLambda();
					}else if(ih==7){
						cont_eleptvscutvars[1] = casc->MassXi();
					}else if(ih==8){
						Double_t lPosV0[3];
						lPosV0[0] = casc->DecayVertexV0X();
						lPosV0[1] = casc->DecayVertexV0Y();
						lPosV0[2] = casc->DecayVertexV0Z();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
					}else if(ih==9){
						Double_t lPosXi[3];
						lPosXi[0] = casc->DecayVertexXiX();
						lPosXi[1] = casc->DecayVertexXiY();
						lPosXi[2] = casc->DecayVertexXiZ();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
					}else if(ih==10){
						cont_eleptvscutvars[1] = casc->DcaV0Daughters();
					}else if(ih==11){
						cont_eleptvscutvars[1] = casc->DcaXiDaughters();
					}else if(ih==12){
						cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
					}else if(ih==13){
						if(casc->ChargeXi()<0.)
							cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
					}else if(ih==14){
						if(casc->ChargeXi()>0.)
							cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
					}else if(ih==15){
						cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
					}else if(ih==16){
						cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
					}else if(ih==17){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
					}else if(ih==18){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
					}else if(ih==19){
						cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
					}else if(ih==20){
						cont_eleptvscutvars[1] =  casc->Eta();
					}else if(ih==21){
						cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
					}else if(ih==22){
						Double_t xipx = exobj->PxProng(1);
						Double_t xipy = exobj->PyProng(1);
						Double_t xipz = exobj->PzProng(1);
						Double_t epx = exobj->PxProng(0);
						Double_t epy = exobj->PyProng(0);
						Double_t epz = exobj->PzProng(0);
						cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
					}else{
						cont_eleptvscutvars[1] = -9999.;
					}

					fHistoElePtvsCutVarsRS[ih]->Fill(cont_eleptvscutvars);
				}
			}
			if(fUseMCInfo){
				if(mcxic){
					Int_t pdgcode = mcxic->GetPdgCode();
					cont2[1] = mcele->Pt();
					if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgcasc_array[1])==4132){
						Int_t labmotherxic = mcxic->GetMother();
						Bool_t isbottomfd = kFALSE;
						if(labmotherxic>=0)
						{
							AliAODMCParticle *motherxic = (AliAODMCParticle*)mcArray->At(labmotherxic);
							Int_t pdgmotherxic = motherxic->GetPdgCode();
							if(TMath::Abs(pdgmotherxic)==511||TMath::Abs(pdgmotherxic)==521||TMath::Abs(pdgmotherxic)==5122||TMath::Abs(pdgmotherxic)==5132||TMath::Abs(pdgmotherxic)==5232||TMath::Abs(pdgmotherxic)==5332){
								isbottomfd = kTRUE;
							}
						}

						fHistoEleXiMassMCS->Fill(cont);
						if(trk->Charge()>0) fHistoEleXiMassMCS1->Fill(cont);
						else  fHistoEleXiMassMCS2->Fill(cont);

						if(isbottomfd){
							fHistoEleXiMassBFeeddownMCS->Fill(cont);
							if(trk->Charge()>0) fHistoEleXiMassBFeeddownMCS1->Fill(cont);
							else  fHistoEleXiMassBFeeddownMCS2->Fill(cont);
						}else{
							fHistoEleXiMassPromptMCS->Fill(cont);
							if(trk->Charge()>0) fHistoEleXiMassPromptMCS1->Fill(cont);
							else  fHistoEleXiMassPromptMCS2->Fill(cont);
						}

						fHistoEleXiMassvsElePtMCS->Fill(cont2);
						if(trk->Charge()>0) fHistoEleXiMassvsElePtMCS1->Fill(cont2);
						else fHistoEleXiMassvsElePtMCS2->Fill(cont2);
						if(cont[0]<fAnalCuts->GetEleXiMassMax()){
							fHistoEleXiPtvsRapidityMCS->Fill(exobj->Pt(),exobj_rap);
							fHistoElePtMCS->Fill(mcele->Pt(),fCentrality);
							fHistoElePtvsEtaMCS->Fill(cont_eleptvseta);
							fHistoElePtvsXiPtMCS->Fill(cont_eleptvsxipt);
							fHistoElePtvsd0MCS->Fill(cont_eleptvsd0);

							Double_t cont_xic[3];
							cont_xic[0] = mcxic->Pt();
							cont_xic[1] = mcxic->Y();
							cont_xic[2] = fCentrality;
							fHistoXicMCS->Fill(cont_xic);
							if(trk->Charge()>0) fHistoXicMCS1->Fill(cont_xic);
							else fHistoXicMCS2->Fill(cont_xic);

							Double_t cont_mcele[3];
							cont_mcele[0] = mcele->Pt();
							cont_mcele[1] = mcele->Eta();
							cont_mcele[2] = fCentrality;
							fHistoXicElectronMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoXicElectronMCS1->Fill(cont_mcele);
							else fHistoXicElectronMCS2->Fill(cont_mcele);

							fHistoResponseElePt->Fill(mcxic->Pt(),trk->Pt());
							fHistoResponseXiPt->Fill(mcxic->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
							fHistoResponseEleXiPt->Fill(mcxic->Pt(),exobj->Pt());
							fHistoResponseXiPtvsEleXiPt->Fill(exobj->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));

							Double_t cont_eleptvsxiptvsxicpt[4];
							cont_eleptvsxiptvsxicpt[0] = cont_eleptvsxipt[0];
							cont_eleptvsxiptvsxicpt[1] = cont_eleptvsxipt[1];
							cont_eleptvsxiptvsxicpt[2] = mcxic->Pt();
							cont_eleptvsxiptvsxicpt[3] = cont_eleptvsxipt[2];
							fHistoElePtvsXiPtvsXicPtMCS->Fill(cont_eleptvsxiptvsxicpt);

							if(isbottomfd){
								fHistoElePtvsd0BFeeddownMCS->Fill(cont_eleptvsd0);
							}else{
								fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
							}

							for(Int_t ih=0;ih<23;ih++){
								Double_t cont_eleptvscutvars[3];
								cont_eleptvscutvars[0] = exobj->Pt();
								cont_eleptvscutvars[2] = fCentrality;

								if(ih==0){
									cont_eleptvscutvars[1] = trk->GetTPCNcls();
								}else if(ih==1){
									cont_eleptvscutvars[1] = trk->GetTPCsignalN();
								}else if(ih==2){
									cont_eleptvscutvars[1] = nSigmaTPCele;
								}else if(ih==3){
									cont_eleptvscutvars[1] = nSigmaTOFele;
								}else if(ih==4){
									cont_eleptvscutvars[1] = trk->Eta();
								}else if(ih==5){
									cont_eleptvscutvars[1] = trk->GetITSNcls();
								}else if(ih==6){
									if(casc->ChargeXi()<0)
										cont_eleptvscutvars[1] = casc->MassLambda();
									else
										cont_eleptvscutvars[1] = casc->MassAntiLambda();
								}else if(ih==7){
									cont_eleptvscutvars[1] = casc->MassXi();
								}else if(ih==8){
									Double_t lPosV0[3];
									lPosV0[0] = casc->DecayVertexV0X();
									lPosV0[1] = casc->DecayVertexV0Y();
									lPosV0[2] = casc->DecayVertexV0Z();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
								}else if(ih==9){
									Double_t lPosXi[3];
									lPosXi[0] = casc->DecayVertexXiX();
									lPosXi[1] = casc->DecayVertexXiY();
									lPosXi[2] = casc->DecayVertexXiZ();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
								}else if(ih==10){
									cont_eleptvscutvars[1] = casc->DcaV0Daughters();
								}else if(ih==11){
									cont_eleptvscutvars[1] = casc->DcaXiDaughters();
								}else if(ih==12){
									cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
								}else if(ih==13){
									if(casc->ChargeXi()<0.)
										cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
								}else if(ih==14){
									if(casc->ChargeXi()>0.)
										cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
								}else if(ih==15){
									cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
								}else if(ih==16){
									cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
								}else if(ih==17){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
								}else if(ih==18){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
								}else if(ih==19){
									cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
								}else if(ih==20){
									cont_eleptvscutvars[1] =  casc->Eta();
								}else if(ih==21){
									cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
								}else if(ih==22){
									Double_t xipx = exobj->PxProng(1);
									Double_t xipy = exobj->PyProng(1);
									Double_t xipz = exobj->PzProng(1);
									Double_t epx = exobj->PxProng(0);
									Double_t epy = exobj->PyProng(0);
									Double_t epz = exobj->PzProng(0);
									cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
								}else{
									cont_eleptvscutvars[1] = -9999.;
								}

								fHistoElePtvsCutVarsMCS[ih]->Fill(cont_eleptvscutvars);
							}
						}
					}
				}
			}
      fHistodPhiSdEtaSElectronProtonR125RS->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125RS->Fill(dphis_ele_pi,detas_ele_pi);
      fHistodPhiSdEtaSElectronBachelorR125RS->Fill(dphis_ele_bach,detas_ele_bach);
		}else{
			fHistoEleXiMassWS->Fill(cont);
			if(trk->Charge()>0) fHistoEleXiMassWS1->Fill(cont);
			else  fHistoEleXiMassWS2->Fill(cont);
			fHistoEleXiMassvsElePtWS->Fill(cont2);
			if(trk->Charge()>0) fHistoEleXiMassvsElePtWS1->Fill(cont2);
			else  fHistoEleXiMassvsElePtWS2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleXiMassMax()){
				fHistoEleXiPtvsRapidityWS->Fill(exobj->Pt(),exobj_rap);
				fHistoElePtWS->Fill(trk->Pt(),fCentrality);
				fHistoElePtvsEtaWS->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtWS->Fill(cont_eleptvsxipt);
				fHistoElePtvsd0WS->Fill(cont_eleptvsd0);

				for(Int_t ih=0;ih<23;ih++){
					Double_t cont_eleptvscutvars[3];
					cont_eleptvscutvars[0] = exobj->Pt();
					cont_eleptvscutvars[2] = fCentrality;

					if(ih==0){
						cont_eleptvscutvars[1] = trk->GetTPCNcls();
					}else if(ih==1){
						cont_eleptvscutvars[1] = trk->GetTPCsignalN();
					}else if(ih==2){
						cont_eleptvscutvars[1] = nSigmaTPCele;
					}else if(ih==3){
						cont_eleptvscutvars[1] = nSigmaTOFele;
					}else if(ih==4){
						cont_eleptvscutvars[1] = trk->Eta();
					}else if(ih==5){
						cont_eleptvscutvars[1] = trk->GetITSNcls();
					}else if(ih==6){
						if(casc->ChargeXi()<0)
							cont_eleptvscutvars[1] = casc->MassLambda();
						else
							cont_eleptvscutvars[1] = casc->MassAntiLambda();
					}else if(ih==7){
						cont_eleptvscutvars[1] = casc->MassXi();
					}else if(ih==8){
						Double_t lPosV0[3];
						lPosV0[0] = casc->DecayVertexV0X();
						lPosV0[1] = casc->DecayVertexV0Y();
						lPosV0[2] = casc->DecayVertexV0Z();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
					}else if(ih==9){
						Double_t lPosXi[3];
						lPosXi[0] = casc->DecayVertexXiX();
						lPosXi[1] = casc->DecayVertexXiY();
						lPosXi[2] = casc->DecayVertexXiZ();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
					}else if(ih==10){
						cont_eleptvscutvars[1] = casc->DcaV0Daughters();
					}else if(ih==11){
						cont_eleptvscutvars[1] = casc->DcaXiDaughters();
					}else if(ih==12){
						cont_eleptvscutvars[1] = casc->DcaBachToPrimVertex();
					}else if(ih==13){
						if(casc->ChargeXi()<0.)
							cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
					}else if(ih==14){
						if(casc->ChargeXi()>0.)
							cont_eleptvscutvars[1] = casc->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = casc->DcaNegToPrimVertex();
					}else if(ih==15){
						cont_eleptvscutvars[1] =  casc->CosPointingAngle(casc->GetDecayVertexXi());
					}else if(ih==16){
						cont_eleptvscutvars[1] =  casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
					}else if(ih==17){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
					}else if(ih==18){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
					}else if(ih==19){
						cont_eleptvscutvars[1] =  nSigmaTPCbachpi;
					}else if(ih==20){
						cont_eleptvscutvars[1] =  casc->Eta();
					}else if(ih==21){
						cont_eleptvscutvars[1] =  0.5*TMath::Log((sqrt(casc->Ptot2Xi())+casc->MomXiZ())/(sqrt(casc->Ptot2Xi())-casc->MomXiZ()));
					}else if(ih==22){
						Double_t xipx = exobj->PxProng(1);
						Double_t xipy = exobj->PyProng(1);
						Double_t xipz = exobj->PzProng(1);
						Double_t epx = exobj->PxProng(0);
						Double_t epy = exobj->PyProng(0);
						Double_t epz = exobj->PzProng(0);
						cont_eleptvscutvars[1] = acos((xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz));
					}else{
						cont_eleptvscutvars[1] = -9999.;
					}

					fHistoElePtvsCutVarsWS[ih]->Fill(cont_eleptvscutvars);
				}
			}
      if(fUseMCInfo){
				if(mcxic){
					Int_t pdgcode = mcxic->GetPdgCode();
					Double_t cont_xib[3];
					cont_xib[0] = mcxic->Pt();
					cont_xib[1] = mcxic->Y();
					cont_xib[2] = fCentrality;

					if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgcasc_array[1])==4132 && abs(mcpdgcasc_array[2])==5132){
						fHistoEleXiMassXibMCS->Fill(cont);
						if(trk->Charge()>0) fHistoEleXiMassXibMCS1->Fill(cont);
						else  fHistoEleXiMassXibMCS2->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleXiMassMax()){
							fHistoXibMCS->Fill(cont_xib);
							fHistoResponseXiPtXib->Fill(mcxic->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
							fHistoResponseEleXiPtXib->Fill(mcxic->Pt(),exobj->Pt());
						}
          }
					if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgcasc_array[1])==4232 && abs(mcpdgcasc_array[2])==5132){
						fHistoEleXiMassXibMCS->Fill(cont);
						if(trk->Charge()>0) fHistoEleXiMassXibMCS1->Fill(cont);
						else  fHistoEleXiMassXibMCS2->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleXiMassMax()){
							fHistoXibMCS->Fill(cont_xib);
							fHistoResponseXiPtXib->Fill(mcxic->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
							fHistoResponseEleXiPtXib->Fill(mcxic->Pt(),exobj->Pt());
						}
          }
					if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgcasc_array[1])==4132 && abs(mcpdgcasc_array[2])==5232){
						fHistoEleXiMassXibMCS->Fill(cont);
						if(trk->Charge()>0) fHistoEleXiMassXibMCS1->Fill(cont);
						else  fHistoEleXiMassXibMCS2->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleXiMassMax()){
							fHistoXibMCS->Fill(cont_xib);
							fHistoResponseXiPtXib->Fill(mcxic->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
							fHistoResponseEleXiPtXib->Fill(mcxic->Pt(),exobj->Pt());
						}
          }
					if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgcasc_array[1])==4232 && abs(mcpdgcasc_array[2])==5232){
						fHistoEleXiMassXibMCS->Fill(cont);
						if(trk->Charge()>0) fHistoEleXiMassXibMCS1->Fill(cont);
						else  fHistoEleXiMassXibMCS2->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleXiMassMax()){
							fHistoXibMCS->Fill(cont_xib);
							fHistoResponseXiPtXib->Fill(mcxic->Pt(),sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)));
							fHistoResponseEleXiPtXib->Fill(mcxic->Pt(),exobj->Pt());
						}
          }
        }
      }
      fHistodPhiSdEtaSElectronProtonR125WS->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125WS->Fill(dphis_ele_pi,detas_ele_pi);
      fHistodPhiSdEtaSElectronBachelorR125WS->Fill(dphis_ele_bach,detas_ele_bach);
		}

	}

	//if( exobj->InvMass(2,pdgdg)<10. && cosoa < 0. && fAnalCuts->IsPeakRegion(casc))
	if( mexi_flip <10. && cosoa < 0. && fAnalCuts->IsPeakRegion(casc))
	{
		if(trk->Charge()*casc->ChargeXi()<0){
			fHistoEleXiMassRSAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleXiMassRS1Away->Fill(cont_flip);
			else  fHistoEleXiMassRS2Away->Fill(cont_flip);
		}else{
			fHistoEleXiMassWSAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleXiMassWS1Away->Fill(cont_flip);
			else  fHistoEleXiMassWS2Away->Fill(cont_flip);
		}
	}

	if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsSideBand(casc))
	{
		if(trk->Charge()*casc->ChargeXi()<0){
			fHistoEleXiMassRSSide->Fill(cont);
			if(trk->Charge()>0) fHistoEleXiMassRSSide1->Fill(cont);
			else  fHistoEleXiMassRSSide2->Fill(cont);
			fHistoEleXiMassvsElePtRSSide->Fill(cont2);
			if(trk->Charge()>0) fHistoEleXiMassvsElePtRSSide1->Fill(cont2);
			else  fHistoEleXiMassvsElePtRSSide2->Fill(cont2);
		}else{
			fHistoEleXiMassWSSide->Fill(cont);
			if(trk->Charge()>0) fHistoEleXiMassWSSide1->Fill(cont);
			else  fHistoEleXiMassWSSide2->Fill(cont);
			fHistoEleXiMassvsElePtWSSide->Fill(cont2);
			if(trk->Charge()>0) fHistoEleXiMassvsElePtWSSide1->Fill(cont2);
			else  fHistoEleXiMassvsElePtWSSide2->Fill(cont2);
		}
	}

	//if( exobj->InvMass(2,pdgdg)<10. && cosoa < 0. && fAnalCuts->IsSideBand(casc))
	if( mexi_flip< 10. && cosoa < 0. && fAnalCuts->IsSideBand(casc))
	{
		if(trk->Charge()*casc->ChargeXi()<0){
			fHistoEleXiMassRSSideAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleXiMassRSSide1Away->Fill(cont_flip);
			else  fHistoEleXiMassRSSide2Away->Fill(cont_flip);
		}else{
			fHistoEleXiMassWSSideAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleXiMassWSSide1Away->Fill(cont_flip);
			else  fHistoEleXiMassWSSide2Away->Fill(cont_flip);
		}
	}

  //
  // New strategy: Fully analyze correlation
  //
  for(Int_t iv=0;iv<15;iv++){
    fCorrelationVariables[iv] = -9999.;
  }
	Double_t cont_cor_nd[7];
  for(Int_t iv=0;iv<7;iv++){
    cont_cor_nd[iv] = -9999.;
  }
	Double_t cont_mass_nd[8];
  for(Int_t iv=0;iv<8;iv++){
    cont_mass_nd[iv] = -9999.;
  }

  fCorrelationVariables[0] = sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2));
  fCorrelationVariables[1] = trk->Pt();
  fCorrelationVariables[2] = fAnalCuts->DeltaPhi(casc,trk);
  fCorrelationVariables[3] = fAnalCuts->DeltaEta(casc,trk);
  fCorrelationVariables[5] = exobj->Getd0Prong(0);
  if(!fUseMCInfo) fCorrelationVariables[6] = 0;
  else fCorrelationVariables[6] = 1;
  if(trk->Charge()>0){
    if(casc->ChargeXi()<0) fCorrelationVariables[7] = 0;
    else fCorrelationVariables[7] = 2;
  }else if(trk->Charge()<0){
    if(casc->ChargeXi()<0) fCorrelationVariables[7] = 1;
    else fCorrelationVariables[7] = 3;
  }
  fCorrelationVariables[8] = (Int_t)isconv + 2 * (Int_t)isconv_like;
  fCorrelationVariables[10] = fCentrality;
  fCorrelationVariables[11] = exobj->Pt();
  fCorrelationVariables[12] = exobj->InvMass(2,pdgdg);
  fCorrelationVariables[13] = fAnalCuts->CosOpeningAngle(casc,trk);

	cont_cor_nd[0] = exobj->Pt();
	cont_cor_nd[1] = fAnalCuts->DeltaPhi(casc,trk);
	cont_cor_nd[2] = 1.;//not used
	if(trk->Charge()>0){
		if(casc->ChargeXi()<0) cont_cor_nd[3] = 0;
		else cont_cor_nd[3] = 2;
	}else if(trk->Charge()<0){
		if(casc->ChargeXi()<0) cont_cor_nd[3] = 3;
		else cont_cor_nd[3] = 1;
	}
	cont_cor_nd[4] = fCorrelationVariables[8];
	cont_cor_nd[5] = 0;
	cont_cor_nd[6] = fCentrality;

  if(fUseMCInfo && FromSemileptonicDecays(mcpdgele_array)>0){
    if(mcxic){
      Int_t pdgcode = mcxic->GetPdgCode();
      if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgcasc_array[1])==4132){
        fCorrelationVariables[9] = 1;
				cont_cor_nd[5] = 1;
      }
      if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgcasc_array[1])==4132 && abs(mcpdgcasc_array[2])==5132){
        fCorrelationVariables[9] = 11;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgcasc_array[1])==4232 && abs(mcpdgcasc_array[2])==5132){
        fCorrelationVariables[9] = 12;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgcasc_array[1])==4132 && abs(mcpdgcasc_array[2])==5232){
        fCorrelationVariables[9] = 13;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgcasc_array[1])==4232 && abs(mcpdgcasc_array[2])==5232){
        fCorrelationVariables[9] = 14;
				cont_cor_nd[5] = 6;
      }
      if(FromSemileptonicDecays(mcpdgele_array)==3){
        fCorrelationVariables[9] = 16;
      }
    }
    if(fCorrelationVariables[9]<0){
      Bool_t lam_from_bottom = HaveBottomInHistory(mcpdgcasc_array);
      Bool_t lam_from_charm = HaveCharmInHistory(mcpdgcasc_array);
      if(FromSemileptonicDecays(mcpdgele_array)==1){
        if(lam_from_bottom) fCorrelationVariables[9] = 1011;
        else if(lam_from_charm) fCorrelationVariables[9] = 1012;
        else  fCorrelationVariables[9] = 1013;
				cont_cor_nd[5] = 7;
      }
      if(FromSemileptonicDecays(mcpdgele_array)==2){
        if(lam_from_bottom) fCorrelationVariables[9] = 1014;
        else if(lam_from_charm) fCorrelationVariables[9] = 1015;
        else  fCorrelationVariables[9] = 1016;
				cont_cor_nd[5] = 8;
      }
      if(FromSemileptonicDecays(mcpdgele_array)==1 && HaveBottomInHistory(mcpdgele_array)){
        if(lam_from_bottom) fCorrelationVariables[9] = 1017;
        else if(lam_from_charm) fCorrelationVariables[9] = 1018;
        else  fCorrelationVariables[9] = 1019;
				cont_cor_nd[5] = 9;
      }
      if(FromSemileptonicDecays(mcpdgele_array)==3){
        if(HaveBottomInHistory(mcpdgele_array)) fCorrelationVariables[9] = 1021;
        else if(HaveCharmInHistory(mcpdgele_array)) fCorrelationVariables[9] = 1022;
        else  fCorrelationVariables[9] = 1023;
      }
    }
  }
  if(fUseMCInfo) fCorrelationVariables[14] = mcpdgele_array[1];

  if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsPeakRegion(casc))
  {
		if(fWriteVariableTree)
			fCorrelationVariablesTree->Fill();

		if(fUseMCInfo){
			if(exobj->InvMass(2,pdgdg)<fAnalCuts->GetEleXiMassMax())
				fHistoCorrelationVariablesvsEleXiPtMC->Fill(cont_cor_nd);
			cont_cor_nd[0] = trk->Pt();
			fHistoCorrelationVariablesvsElePtMC->Fill(cont_cor_nd);
			cont_cor_nd[0] = sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2));
			fHistoCorrelationVariablesvsXiPtMC->Fill(cont_cor_nd);
		}else{
			if(exobj->InvMass(2,pdgdg)<fAnalCuts->GetEleXiMassMax())
				fHistoCorrelationVariablesvsEleXiPt->Fill(cont_cor_nd);
			cont_cor_nd[0] = trk->Pt();
			fHistoCorrelationVariablesvsElePt->Fill(cont_cor_nd);
			cont_cor_nd[0] = sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2));
			fHistoCorrelationVariablesvsXiPt->Fill(cont_cor_nd);
		}
  }

	cont_mass_nd[0] =  exobj->InvMass(2,pdgdg);
	cont_mass_nd[1] =  cont_cor_nd[0];
	cont_mass_nd[4] =  cont_cor_nd[3];
	cont_mass_nd[5] =  cont_cor_nd[4];
	cont_mass_nd[6] =  cont_cor_nd[5];
	cont_mass_nd[7] =  cont_cor_nd[6];
	if(fAnalCuts->IsPeakRegion(casc)) cont_mass_nd[3] = 1;
	if(fAnalCuts->IsSideBand(casc)) cont_mass_nd[3] = 0;
	if(fAnalCuts->IsSelected(exobj,AliRDHFCuts::kCandidate)) cont_mass_nd[2]=1;
	if(mexi_flip < 10.&& cosoa < 0.) cont_mass_nd[2]=0;
	if(fUseMCInfo){
		fHistoMassVariablesvsEleXiPtMC->Fill(cont_cor_nd);
		cont_cor_nd[0] = trk->Pt();
		fHistoMassVariablesvsElePtMC->Fill(cont_cor_nd);
		cont_cor_nd[0] = sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2));
		fHistoMassVariablesvsXiPtMC->Fill(cont_cor_nd);
	}else{
		fHistoMassVariablesvsEleXiPt->Fill(cont_cor_nd);
		cont_cor_nd[0] = trk->Pt();
		fHistoMassVariablesvsElePt->Fill(cont_cor_nd);
		cont_cor_nd[0] = sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2));
		fHistoMassVariablesvsXiPt->Fill(cont_cor_nd);
	}


  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMixROOTObjects(TLorentzVector *trke, TLorentzVector *casc, TVector *elevars, TVector *cascvars, Int_t chargexi) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //
	if(!trke) return;
	if(!casc) return;


	for(Int_t i=0;i<90;i++){
		fCandidateVariables[i] = -9999.;
	}

	Double_t pxe = trke->Px();
	Double_t pye = trke->Py();
	Double_t pze = trke->Pz();
	Double_t mome = sqrt(pxe*pxe+pye*pye+pze*pze);
	Double_t Ee = sqrt(mome*mome+0.000510998928*0.000510998928);

	Double_t pxv = casc->Px();
	Double_t pyv = casc->Py();
	Double_t pzv = casc->Pz();
	Double_t momv = sqrt(pxv*pxv+pyv*pyv+pzv*pzv);
	Double_t Ev = sqrt(momv*momv+1.32171*1.32171);

	Double_t cosoa = (pxe*pxv+pye*pyv+pze*pzv)/mome/momv;

	Double_t pxsum = pxe + pxv;
	Double_t pysum = pye + pyv;
	Double_t pzsum = pze + pzv;
	Double_t Esum = Ee + Ev;
	Double_t mexi = sqrt(Esum*Esum-pxsum*pxsum-pysum*pysum-pzsum*pzsum);

	Double_t uxe = pxe/mome;
	Double_t uye = pye/mome;
	Double_t uze = pze/mome;
	Double_t lf = -2.*(pxv*uxe+pyv*uye+pzv*uze);
	Double_t pxv_flip = pxv + lf * uxe;
	Double_t pyv_flip = pyv + lf * uye;
	Double_t pzv_flip = pzv + lf * uze;
	Double_t pxsum_flip = pxe + pxv_flip;
	Double_t pysum_flip = pye + pyv_flip;
	Double_t pzsum_flip = pze + pzv_flip;
	Double_t mexi_flip = sqrt(Esum*Esum-pxsum_flip*pxsum_flip-pysum_flip*pysum_flip-pzsum_flip*pzsum_flip);
	Double_t ptexi_flip = sqrt(pxsum_flip*pxsum_flip+pysum_flip*pysum_flip);

  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);

  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3312};
  fCandidateVariables[ 1] = mexi;
  fCandidateVariables[ 2] = sqrt(pxsum*pxsum+pysum*pysum);
  fCandidateVariables[ 3] = pxsum;
  fCandidateVariables[ 4] = pysum;
  fCandidateVariables[ 5] = pzsum;
  fCandidateVariables[ 6] = pxe;
  fCandidateVariables[ 7] = pye;
  fCandidateVariables[ 8] = pze;
  fCandidateVariables[ 9] = pxv;
  fCandidateVariables[10] = pyv;
  fCandidateVariables[11] = pzv;
  fCandidateVariables[12] = chargexi;
  fCandidateVariables[13] = casc->M();
  fCandidateVariables[29] = trke->T();
  fCandidateVariables[30] = 1;//mixing
	fCandidateVariables[77] = fVtx1->GetX();
	fCandidateVariables[78] = fVtx1->GetY();
	fCandidateVariables[79] = fVtx1->GetZ();
  fCandidateVariables[88] = fEvNumberCounter;

//  if(fWriteVariableTree)
//    fVariablesTree->Fill();

	Double_t cont[4];
	cont[0] = mexi;
	cont[1] = sqrt(pxsum*pxsum+pysum*pysum);
	cont[2] = 0.0;
	cont[3] = fCentrality;

	Double_t cont_flip[4];
	cont_flip[0] = mexi_flip;
	cont_flip[1] = ptexi_flip;
	cont_flip[2] = 0.0;
	cont_flip[3] = fCentrality;

	Double_t cont2[3];
	cont2[0] = mexi;
	cont2[1] = sqrt(pxe*pxe+pye*pye);
	cont2[2] = fCentrality;

	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = trke->Pt();
	cont_eleptvseta[1] = trke->Eta();
	cont_eleptvseta[2] = fCentrality;

	Double_t cont_eleptvsxipt[3];
	cont_eleptvsxipt[0] = trke->Pt();
	cont_eleptvsxipt[1] = casc->Pt();
	cont_eleptvsxipt[2] = fCentrality;

	Double_t cont_eleptvsd0[3];
	cont_eleptvsd0[0] = trke->Pt();
	cont_eleptvsd0[1] = 0.;
	cont_eleptvsd0[2] = fCentrality;

  Double_t xyzR125_ele[3], xyzR125_pr[3], xyzR125_pi[3], xyzR125_bach[3];
  xyzR125_ele[0] = (*elevars)[0];
  xyzR125_ele[1] = (*elevars)[1];
  xyzR125_ele[2] = (*elevars)[2];
  xyzR125_pr[0] = (*cascvars)[0];
  xyzR125_pr[1] = (*cascvars)[1];
  xyzR125_pr[2] = (*cascvars)[2];
  xyzR125_pi[0] = (*cascvars)[3];
  xyzR125_pi[1] = (*cascvars)[4];
  xyzR125_pi[2] = (*cascvars)[5];
  xyzR125_bach[0] = (*cascvars)[6];
  xyzR125_bach[1] = (*cascvars)[7];
  xyzR125_bach[2] = (*cascvars)[8];

  Double_t rdhfcutvars[12];
  rdhfcutvars[0] = xyzR125_ele[0];
  rdhfcutvars[1] = xyzR125_ele[1];
  rdhfcutvars[2] = xyzR125_ele[2];
  rdhfcutvars[3] = xyzR125_pr[0];
  rdhfcutvars[4] = xyzR125_pr[1];
  rdhfcutvars[5] = xyzR125_pr[2];
  rdhfcutvars[6] = xyzR125_pi[0];
  rdhfcutvars[7] = xyzR125_pi[1];
  rdhfcutvars[8] = xyzR125_pi[2];
  rdhfcutvars[9] = xyzR125_bach[0];
  rdhfcutvars[10] = xyzR125_bach[1];
  rdhfcutvars[11] = xyzR125_bach[2];

  Double_t dphis_ele_pr = fAnalCuts->dPhiSR125(xyzR125_ele,xyzR125_pr);
  Double_t detas_ele_pr = fAnalCuts->dEtaSR125(xyzR125_ele,xyzR125_pr);
  Double_t dphis_ele_pi = fAnalCuts->dPhiSR125(xyzR125_ele,xyzR125_pi);
  Double_t detas_ele_pi = fAnalCuts->dEtaSR125(xyzR125_ele,xyzR125_pi);
  Double_t dphis_ele_bach = fAnalCuts->dPhiSR125(xyzR125_ele,xyzR125_bach);
  Double_t detas_ele_bach = fAnalCuts->dEtaSR125(xyzR125_ele,xyzR125_bach);

	//if(mexi<10. && cosoa>0. && fAnalCuts->IsPeakRegion(casc))
  if(fAnalCuts->IsSelected(trke,casc,rdhfcutvars,AliRDHFCuts::kCandidate) &&  fAnalCuts->IsPeakRegion(casc))
	{
		if(((int)trke->T())*chargexi<0){
			fHistoEleXiMassRSMix->Fill(cont);
			if(trke->T()>0) fHistoEleXiMassRSMix1->Fill(cont);
			else  fHistoEleXiMassRSMix2->Fill(cont);
			fHistoEleXiMassvsElePtRSMix->Fill(cont2);
			if(trke->T()>0) fHistoEleXiMassvsElePtRSMix1->Fill(cont2);
			else  fHistoEleXiMassvsElePtRSMix2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleXiMassMax()){
				fHistoElePtRSMix->Fill(trke->Pt(),fCentrality);
				fHistoElePtvsEtaRSMix->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtRSMix->Fill(cont_eleptvsxipt);
				fHistoElePtvsd0RSMix->Fill(cont_eleptvsd0);
			}
      fHistodPhiSdEtaSElectronProtonR125RSMix->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125RSMix->Fill(dphis_ele_pi,detas_ele_pi);
      fHistodPhiSdEtaSElectronBachelorR125RSMix->Fill(dphis_ele_bach,detas_ele_bach);
		}else{
			fHistoEleXiMassWSMix->Fill(cont);
			if(trke->T()>0) fHistoEleXiMassWSMix1->Fill(cont);
			else  fHistoEleXiMassWSMix2->Fill(cont);
			fHistoEleXiMassvsElePtWSMix->Fill(cont2);
			if(trke->T()>0) fHistoEleXiMassvsElePtWSMix1->Fill(cont2);
			else  fHistoEleXiMassvsElePtWSMix2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleXiMassMax()){
				fHistoElePtWSMix->Fill(trke->Pt(),fCentrality);
				fHistoElePtvsEtaWSMix->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtWSMix->Fill(cont_eleptvsxipt);
				fHistoElePtvsd0WSMix->Fill(cont_eleptvsd0);
			}
      fHistodPhiSdEtaSElectronProtonR125WSMix->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125WSMix->Fill(dphis_ele_pi,detas_ele_pi);
      fHistodPhiSdEtaSElectronBachelorR125WSMix->Fill(dphis_ele_bach,detas_ele_bach);
		}
	}

	//if(mexi < 10. && cosoa<0. && fAnalCuts->IsPeakRegion(casc))
	if(mexi_flip < 10. && cosoa<0. && fAnalCuts->IsPeakRegion(casc))
	{
		if(((int)trke->T())*chargexi<0){
			fHistoEleXiMassRSMixAway->Fill(cont_flip);
			if(trke->T()>0) fHistoEleXiMassRSMix1Away->Fill(cont_flip);
			else  fHistoEleXiMassRSMix2Away->Fill(cont_flip);
		}else{
			fHistoEleXiMassWSMixAway->Fill(cont_flip);
			if(trke->T()>0) fHistoEleXiMassWSMix1Away->Fill(cont_flip);
			else  fHistoEleXiMassWSMix2Away->Fill(cont_flip);
		}
	}

  //
  // New strategy: Fully analyze correlation
  //
  for(Int_t iv=0;iv<15;iv++){
    fCorrelationVariables[iv] = -9999.;
  }

	Double_t cont_cor_nd[7];
  for(Int_t iv=0;iv<7;iv++){
    cont_cor_nd[iv] = -9999.;
  }
	Double_t cont_mass_nd[8];
  for(Int_t iv=0;iv<8;iv++){
    cont_mass_nd[iv] = -9999.;
  }

  fCorrelationVariables[0] = casc->Pt();
  fCorrelationVariables[1] = trke->Pt();
  fCorrelationVariables[2] = TVector2::Phi_mpi_pi(casc->Phi()-trke->Phi());
  fCorrelationVariables[3] = casc->Eta()-trke->Eta();
  fCorrelationVariables[5] = (*elevars)[5];
  fCorrelationVariables[6] = 2;
  if(trke->T()>0){
    if(chargexi<0) fCorrelationVariables[7] = 0;
    else fCorrelationVariables[7] = 2;
  }else if(trke->T()<0){
    if(chargexi<0) fCorrelationVariables[7] = 1;
    else fCorrelationVariables[7] = 3;
  }
  fCorrelationVariables[8] = (*elevars)[6];
  fCorrelationVariables[9] = (*elevars)[7];
  fCorrelationVariables[10] = fCentrality;
  fCorrelationVariables[11] = sqrt(pxsum*pxsum+pysum*pysum);
  fCorrelationVariables[12] = mexi;
  fCorrelationVariables[13] = cosoa;

	cont_cor_nd[0] =  sqrt(pxsum*pxsum+pysum*pysum);
	cont_cor_nd[1] =  TVector2::Phi_mpi_pi(casc->Phi()-trke->Phi());
	cont_cor_nd[2] = 1.;//not used
  if(trke->T()>0){
    if(chargexi<0) cont_cor_nd[3] = 0;
    else cont_cor_nd[3] = 2;
  }else if(trke->T()<0){
    if(chargexi<0) cont_cor_nd[3] = 3;
    else cont_cor_nd[3] = 1;
  }
	cont_cor_nd[4] = fCorrelationVariables[8];
	cont_cor_nd[5] = 0;
	if(fabs(fCorrelationVariables[9]-1013)<0.001) cont_cor_nd[5] = 7;
	if(fabs(fCorrelationVariables[9]-1016)<0.001) cont_cor_nd[5] = 8;
	if(fabs(fCorrelationVariables[9]-1019)<0.001) cont_cor_nd[5] = 9;
	cont_cor_nd[6] = fCentrality;

  if(fAnalCuts->IsSelected(trke,casc,rdhfcutvars,AliRDHFCuts::kCandidate) &&  fAnalCuts->IsPeakRegion(casc))
  {
		if(fWriteVariableTree)
			fCorrelationVariablesTree->Fill();

		if(mexi<fAnalCuts->GetEleXiMassMax())
			fHistoCorrelationVariablesvsEleXiPtMix->Fill(cont_cor_nd);
		cont_cor_nd[0] = trke->Pt();
		fHistoCorrelationVariablesvsElePtMix->Fill(cont_cor_nd);
		cont_cor_nd[0] = casc->Pt();
		fHistoCorrelationVariablesvsXiPtMix->Fill(cont_cor_nd);
  }

	cont_mass_nd[0] =  mexi;
	cont_mass_nd[1] =  cont_cor_nd[0];
	cont_mass_nd[4] =  cont_cor_nd[3];
	cont_mass_nd[5] =  cont_cor_nd[4];
	cont_mass_nd[6] =  cont_cor_nd[5];
	cont_mass_nd[7] =  cont_cor_nd[6];
	if(fAnalCuts->IsPeakRegion(casc)) cont_mass_nd[3] = 1;
	if(fAnalCuts->IsSideBand(casc)) cont_mass_nd[3] = 0;
	if(fAnalCuts->IsSelected(trke,casc,rdhfcutvars,AliRDHFCuts::kCandidate) ) cont_mass_nd[2]=1;
	if(mexi_flip < 10.&& cosoa < 0.) cont_mass_nd[2]=0;
	fHistoMassVariablesvsEleXiPtMix->Fill(cont_mass_nd);
	cont_mass_nd[0] = trke->Pt();
	fHistoMassVariablesvsElePtMix->Fill(cont_mass_nd);
	cont_mass_nd[0] = casc->Pt();
	fHistoMassVariablesvsXiPtMix->Fill(cont_mass_nd);

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineSingleTreeVariables() 
{
  //
  /// Define single tree variables
  //

  const char* nameoutput = GetOutputSlot(11)->GetContainer()->GetName();
  fSingleVariablesTree = new TTree(nameoutput,"single variables tree");
  Int_t nVar = 16;
  fCandidateSingleVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Px";
  fCandidateVariableNames[ 1]="Py";
  fCandidateVariableNames[ 2]="Pz";
  fCandidateVariableNames[ 3]="Charge";
  fCandidateVariableNames[ 4]="XiMass";
  fCandidateVariableNames[ 5]="Bz";
  fCandidateVariableNames[ 6]="Centrality";
  fCandidateVariableNames[ 7]="PrimVertZ";
  fCandidateVariableNames[ 8]="EvNumber";
  fCandidateVariableNames[ 9]="RunNumber";
  fCandidateVariableNames[10]="EventPlane";
  fCandidateVariableNames[11]="Vars0";//e: dca, X: 
  fCandidateVariableNames[12]="Vars1";//e: trk ID, X: trk ID (bach)
  fCandidateVariableNames[13]="Vars2";//e: nSigma TPC, X: trk ID (pos)
  fCandidateVariableNames[14]="Vars3";//e: nSigma TOF X: trk ID (neg)
  fCandidateVariableNames[15]="Vars4";//e: nSigma ITS  X: 

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fSingleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateSingleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineEleTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fEleVariablesTree = new TTree(nameoutput,"electron variables tree");
  Int_t nVar = 26;
  fCandidateEleVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="ElePx";
  fCandidateVariableNames[ 1]="ElePy";
  fCandidateVariableNames[ 2]="ElePz";
  fCandidateVariableNames[ 3]="TPCChi2overNDF";
  fCandidateVariableNames[ 4]="ITSNcls";
  fCandidateVariableNames[ 5]="TPCNcls";
  fCandidateVariableNames[ 6]="TPCNclsPID";
  fCandidateVariableNames[ 7]="TPCNclsRatio";
  fCandidateVariableNames[ 8]="d0R";
  fCandidateVariableNames[ 9]="d0Z";
  fCandidateVariableNames[10]="ITSClusterMap";
  fCandidateVariableNames[11]="nSigmaTPCele";
  fCandidateVariableNames[12]="nSigmaTOFele";
  fCandidateVariableNames[13]="nSigmaTPCpi";
  fCandidateVariableNames[14]="nSigmaTPCka";
  fCandidateVariableNames[15]="nSigmaTPCpr";
  fCandidateVariableNames[16]="EvNumber";
  fCandidateVariableNames[17]="EleCharge";
  fCandidateVariableNames[18]="ElePdgCode";
  fCandidateVariableNames[19]="EleMotherPdgCode";
  fCandidateVariableNames[20]="mcelepx";
  fCandidateVariableNames[21]="mcelepy";
  fCandidateVariableNames[22]="mcelepz";
  fCandidateVariableNames[23]="Centrality";
  fCandidateVariableNames[24]="PrimVertZ";
  fCandidateVariableNames[25]="RunNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fEleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateEleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillElectronROOTObjects(AliAODTrack *trk, AliAODEvent *event, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree 
  //

	if(!trk) return;

  AliAODTrack *trkpid = 0;
  if(fAnalCuts->GetProdAODFilterBit()==7){
    trkpid = fGTI[-trk->GetID()-1];
  }else{
    trkpid = trk;
  }

	fHistoBachPt->Fill(trk->Pt());
	fHistoElectronQovPtvsPhi->Fill(trk->Phi(),(Double_t)trk->Charge()/trk->Pt());

  Double_t d0z0[2],covd0z0[3];
  trk->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);
  fHistod0Bach->Fill(d0z0[0]);

  Double_t minmass_ee = 9999.;
  Double_t minmasslike_ee = 9999.;
  Bool_t isconv = fAnalCuts->TagConversions(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmass_ee);
  Bool_t isconv_like = fAnalCuts->TagConversionsSameSign(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmasslike_ee);
  fHistoMassConversionsMin->Fill(minmass_ee);
  fHistoMassConversionsSameSignMin->Fill(minmasslike_ee);

  Int_t mcetype = -9999;
	Int_t pdgEle = -9999;
	Int_t pdgEleMother = -9999;
	Float_t mcepx = -9999;
	Float_t mcepy = -9999;
	Float_t mcepz = -9999;
  Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;
  for(Int_t i=0;i<100;i++){
    pdgarray_ele[i]=-9999;
    labelarray_ele[i]=-9999;
  }
  ngen_ele = -9999;
  if(fUseMCInfo)
  {
    Int_t labEle = trk->GetLabel();
    if(labEle<0) return;
    AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labEle);
    if(!mcetrk) return;
    pdgEle = mcetrk->GetPdgCode();
    if(abs(pdgEle)!=11) return;
    mcepx = mcetrk->Px();
    mcepy = mcetrk->Py();
    mcepz = mcetrk->Pz();

    Int_t labEleMother = mcetrk->GetMother();
    if(labEleMother>-1){
      AliAODMCParticle *mcemothertrk = (AliAODMCParticle*)mcArray->At(labEleMother);
      if(mcemothertrk){
        pdgEleMother = mcemothertrk->GetPdgCode();
      }
    }

    GetMCDecayHistory(mcetrk,mcArray,pdgarray_ele,labelarray_ele,ngen_ele);

    Int_t hfe_flag = 0;
    Bool_t gamma_flag = kFALSE;
    Bool_t pi0_flag = kFALSE;
    Bool_t eta_flag = kFALSE;
    Double_t pt_pi0 = -9999.;
    Double_t pt_eta = -9999.;
    if(abs(pdgarray_ele[0])>400&&abs(pdgarray_ele[0])<440){
      hfe_flag = 1;
    }
    if(abs(pdgarray_ele[0])>4000&&abs(pdgarray_ele[0])<4400){
      hfe_flag = 1;
    }
    if(abs(pdgarray_ele[0])>500&&abs(pdgarray_ele[0])<540){
      hfe_flag = 2;
    }
    if(abs(pdgarray_ele[0])>5000&&abs(pdgarray_ele[0])<5400){
      hfe_flag = 2;
    }
    if(abs(pdgarray_ele[0])==22){
      gamma_flag = kTRUE;
    }
    if(!gamma_flag){
      fHistoBachPtMCS->Fill(trk->Pt());
    }
    if((abs(pdgarray_ele[0])==22) && (abs(pdgarray_ele[1])==111)){
      pi0_flag = kTRUE;
      AliAODMCParticle *mctrkm = (AliAODMCParticle*)mcArray->At(labelarray_ele[1]);
      pt_pi0 = mctrkm->Pt();
    }
    if(abs(pdgarray_ele[0])==111){
      pi0_flag = kTRUE;
      AliAODMCParticle *mctrkm = (AliAODMCParticle*)mcArray->At(labelarray_ele[0]);
      pt_pi0 = mctrkm->Pt();
    }
    if((abs(pdgarray_ele[0])==22) && (abs(pdgarray_ele[1])==221)){
      eta_flag = kTRUE;
      AliAODMCParticle *mctrkm = (AliAODMCParticle*)mcArray->At(labelarray_ele[1]);
      pt_eta = mctrkm->Pt();
    }
    if(abs(pdgarray_ele[0])==221){
      eta_flag = kTRUE;
      AliAODMCParticle *mctrkm = (AliAODMCParticle*)mcArray->At(labelarray_ele[0]);
      pt_eta = mctrkm->Pt();
    }

    if(pi0_flag){
      Double_t cont_pi0[3];
      cont_pi0[0] = trk->Pt();
      cont_pi0[1] = pt_pi0;
      cont_pi0[2] = fCentrality;
      fHistoElectronPi0Total->Fill(cont_pi0);
      if(isconv) fHistoElectronPi0Tag->Fill(cont_pi0);
    }

    if(eta_flag){
      Double_t cont_eta[3];
      cont_eta[0] = trk->Pt();
      cont_eta[1] = pt_eta;
      cont_eta[2] = fCentrality;
      fHistoElectronEtaTotal->Fill(cont_eta);
      if(isconv) fHistoElectronEtaTag->Fill(cont_eta);
    }

    if(hfe_flag==0) return;

    if(hfe_flag==1){
      mcetype = 1013;
    }
    if(hfe_flag==2){
      mcetype = 1016;
    }
    if(hfe_flag==1 && HaveBottomInHistory(pdgarray_ele)){
      mcetype = 1019;
    }
  }


	if(fDoEventMixing && !(fMixWithoutConversionFlag && isconv)){
    Double_t pv[3];
    pv[0] = fVtx1->GetX();
    pv[1] = fVtx1->GetY();
    pv[2] = fVtx1->GetZ();
    Double_t xyzR125[3] = {9999.,9999.,9999.};
    if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(trk,fBzkG,pv,xyzR125);
    TVector *varvec = new TVector(8);
    (*varvec)[0] = xyzR125[0];
    (*varvec)[1] = xyzR125[1];
    (*varvec)[2] = xyzR125[2];
    (*varvec)[3] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
    (*varvec)[4] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);
    (*varvec)[5] = d0z0[0];
    (*varvec)[6] = (Int_t)isconv + 2 * (Int_t)isconv_like;
    (*varvec)[7] = mcetype;

    Int_t nextRes( nextResVec[fPoolIndex] );
    m_ReservoirE[fPoolIndex][nextRes].push_back(new TLorentzVector(trk->Px(),trk->Py(),trk->Pz(),trk->Charge()));
    m_ReservoirVarsE[fPoolIndex][nextRes].push_back(varvec);
	}

	if(!fWriteEachVariableTree) return;


	for(Int_t i=0;i<26;i++){
		fCandidateEleVariables[i] = -9999.;
	}
	for(Int_t i=0;i<16;i++){
		fCandidateSingleVariables[i] = -9999.;
	}

  fCandidateSingleVariables[ 0] = trk->Px();
  fCandidateSingleVariables[ 1] = trk->Py();
  fCandidateSingleVariables[ 2] = trk->Pz();
  fCandidateSingleVariables[ 3] = trk->Charge();
  fCandidateSingleVariables[ 4] = -9999.;//not lambda
  fCandidateSingleVariables[ 5] = fBzkG;
  fCandidateSingleVariables[ 6] = fCentrality;
  fCandidateSingleVariables[ 7] = fVtxZ;
  fCandidateSingleVariables[ 8] = fEvNumberCounter;
  fCandidateSingleVariables[ 9] = fRunNumber;
  if(fUseEventPlane==4){
    //Double_t evplane = GetEventPlaneForCandidate(trkpid,0,event->GetEventplane(),fQ,fQSub1,fQSub2); // remove autocorrelations
    //fCandidateSingleVariables[10] = evplane;
    fCandidateSingleVariables[10] = fEventPlane;
  }else{
    fCandidateSingleVariables[10] = fEventPlane;
  }
  fCandidateSingleVariables[11] = d0z0[0];
  fCandidateSingleVariables[12] = trkpid->GetID();
  if(fAnalCuts->GetIsUsePID())
  {
    Double_t nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
    Double_t nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
    Double_t nSigmaITSele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasITS(trk,AliPID::kElectron);
    fCandidateSingleVariables[13] = nSigmaTPCele;
    fCandidateSingleVariables[14] = nSigmaTOFele;
    fCandidateSingleVariables[15] = nSigmaITSele;
  }

	fSingleVariablesTree->Fill();

//  fCandidateEleVariables[ 0] = trk->Px();
//  fCandidateEleVariables[ 1] = trk->Py();
//  fCandidateEleVariables[ 2] = trk->Pz();
//  fCandidateEleVariables[ 3] = trk->Chi2perNDF();
//  fCandidateEleVariables[ 4] = trk->GetITSNcls();
//  fCandidateEleVariables[ 5] = trk->GetTPCncls();
//  fCandidateEleVariables[ 6] = trk->GetTPCsignalN();
//	if(trk->GetTPCNclsF()>0) 
//		fCandidateEleVariables[ 7] = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();
//
//  fCandidateEleVariables[ 8] = d0z0[0];
//  fCandidateEleVariables[ 9] = d0z0[1];
//	Int_t itsmap = trk->GetITSClusterMap();
//	Int_t bit1 = 1;
//	Int_t bit2 = 2;
//	Bool_t spdfirst = (itsmap & bit1) == bit1;
//	Bool_t spdsecond = (itsmap & bit2) == bit2;
//  fCandidateEleVariables[10] = ((Int_t)spdfirst) + 2 * ((Int_t)spdsecond);
//
//  if(fAnalCuts->GetIsUsePID())
//  {
//		Double_t nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
//		Double_t nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
//		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
//		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
//		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
//    fCandidateEleVariables[11] = nSigmaTPCele;
//    fCandidateEleVariables[12] = nSigmaTOFele;
//    fCandidateEleVariables[13] = nSigmaTPCpi_etrk;
//    fCandidateEleVariables[14] = nSigmaTPCka_etrk;
//    fCandidateEleVariables[15] = nSigmaTPCpr_etrk;
//  }
//  fCandidateEleVariables[16] = fEvNumberCounter;
//  fCandidateEleVariables[17] = trk->Charge();
//  fCandidateEleVariables[18] = pdgEle;
//  fCandidateEleVariables[19] = pdgEleMother;
//  fCandidateEleVariables[20] = mcepx;
//  fCandidateEleVariables[21] = mcepy;
//  fCandidateEleVariables[22] = mcepz;
//  fCandidateEleVariables[23] = fCentrality;
//  fCandidateEleVariables[24] = fVtxZ;
//  fCandidateEleVariables[25] = fRunNumber;
//
//
//	fEleVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineCascTreeVariables() 
{
  //
  // Define V0 tree variables
  //

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fCascVariablesTree = new TTree(nameoutput,"cascade variables tree");
  Int_t nVar = 26;
  fCandidateCascVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassXi";
  fCandidateVariableNames[ 2]="XiPx";
  fCandidateVariableNames[ 3]="XiPy";
  fCandidateVariableNames[ 4]="XiPz";
  fCandidateVariableNames[ 5]="InvMassLambda";
  fCandidateVariableNames[ 6]="DcaXiDaughters";
  fCandidateVariableNames[ 7]="DcaV0Daughters";
  fCandidateVariableNames[ 8]="DecayLengthXi";
  fCandidateVariableNames[ 9]="CosPointingAngleXi";
  fCandidateVariableNames[10]="DcaV0toPrimVertex";
  fCandidateVariableNames[11]="DcaPostoPrimVertex";
  fCandidateVariableNames[12]="DcaNegtoPrimVertex";
  fCandidateVariableNames[13]="DcaBachtoPrimVertex";
  fCandidateVariableNames[14]="DecayLengthV0";
  fCandidateVariableNames[15]="CosPointingAngleV0";
  fCandidateVariableNames[16]="XiCharge";
  fCandidateVariableNames[17]="XiPdgCode";
  fCandidateVariableNames[18]="XiMotherPdgCode";
  fCandidateVariableNames[19]="mcxipx";
  fCandidateVariableNames[20]="mcxipy";
  fCandidateVariableNames[21]="mcxipz";
  fCandidateVariableNames[22]="labcasc";
  fCandidateVariableNames[23]="RunNumber";
  fCandidateVariableNames[24]="PrimVertZ";
  fCandidateVariableNames[25]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fCascVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateCascVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillCascROOTObjects(AliAODcascade *casc, AliAODEvent *event, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteVariableTree (tree not implemented yet)
  //
	if(!casc) return;
  AliAODTrack *cptrack =  (AliAODTrack*)(casc->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(casc->GetDaughter(1));
  AliAODTrack *cbtrack =  (AliAODTrack*)(casc->GetDecayVertexXi()->GetDaughter(0));
  if(!cptrack) return;
  if(!cntrack) return;
  if(!cbtrack) return;

	fHistoXiMassvsPt->Fill(casc->MassXi(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
	fHistoOmegaMassvsPt->Fill(casc->MassOmega(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));
	Double_t momxix = casc->MomXiX();
	Double_t momxiy = casc->MomXiY();
	Double_t phi_alice = atan2(momxiy,momxix);
	if(phi_alice<0.) phi_alice += 2 * M_PI;
	fHistoXiQovPtvsPhi->Fill(phi_alice,(Double_t)casc->ChargeXi()/sqrt(momxix*momxix+momxiy*momxiy));

  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);
  Double_t ptotlam = TMath::Sqrt(pow(casc->Px(),2)+pow(casc->Py(),2)+pow(casc->Pz(),2));
  Double_t dl = sqrt(pow(casc->DecayVertexV0X()-posVtx[0],2)+pow(casc->DecayVertexV0Y()-posVtx[1],2)+pow(casc->DecayVertexV0Z()-posVtx[2],2));
  Double_t v0propdl = dl*mlamPDG/ptotlam;
	if(fAnalCuts->IsPeakRegion(casc)){
		fHistoLambdaPtvsDl->Fill(casc->Pt(),v0propdl);
		fHistoLambdaPtvsDR->Fill(casc->Pt(),dl);
	}
	if(fAnalCuts->IsSideBand(casc)){
		fHistoLambdaPtvsDlSide->Fill(casc->Pt(),v0propdl);
		fHistoLambdaPtvsDRSide->Fill(casc->Pt(),dl);
	}

	Int_t xipdgcode = -9999;
	Int_t ximotherpdgcode = -9999;
	Float_t mcxipx = -9999.;
	Float_t mcxipy = -9999.;
	Float_t mcxipz = -9999.;
	Int_t labcasc = -9999.;
	if(fUseMCInfo){
		Int_t pdgDgcasc[2]={211,3122};
		Int_t pdgDgv0[2]={2212,211};
		labcasc = MatchToMCCascade(casc,3312,pdgDgcasc,pdgDgv0,mcArray); // the cascade
		if(labcasc<0) return;

		fHistoXiMassvsPtMCS->Fill(casc->MassXi(),sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY()));

		AliAODMCParticle *mccasctrk = (AliAODMCParticle*)mcArray->At(labcasc);
		if(!mccasctrk) return;

		fHistoLambdaPtvsDlMCS->Fill(casc->Pt(),v0propdl);
		fHistoLambdaPtvsDRMCS->Fill(casc->Pt(),dl);

		Bool_t hfxi_flag = kFALSE;
		xipdgcode = mccasctrk->GetPdgCode();
		Int_t labcascmother = mccasctrk->GetMother();
		if(labcascmother>=0){
			AliAODMCParticle *mothercasc = (AliAODMCParticle*)mcArray->At(labcascmother);
			if(mothercasc){
				ximotherpdgcode = mothercasc->GetPdgCode();
				if(abs(ximotherpdgcode)>4000&&abs(ximotherpdgcode)<4400){
					hfxi_flag = kTRUE;
				}
			}
		}
		if(!hfxi_flag) return;
		mcxipx = mccasctrk->Px();
		mcxipy = mccasctrk->Py();
		mcxipz = mccasctrk->Pz();
	}


  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
  xyz[0]=casc->DecayVertexXiX();
  xyz[1]=casc->DecayVertexXiY();
  xyz[2]=casc->DecayVertexXiZ();
  pxpypz[0]=casc->MomXiX();
  pxpypz[1]=casc->MomXiY();
  pxpypz[2]=casc->MomXiZ();
  casc->GetCovarianceXYZPxPyPz(cv);
  sign=casc->ChargeXi();
  AliExternalTrackParam	*trackCasc = new AliExternalTrackParam(xyz,pxpypz,cv,sign);
  trackCasc->PropagateToDCA(fVtx1,fBzkG,kVeryBig);
  Double_t momcasc_new[3]={-9999,-9999,-9999.};
  trackCasc->GetPxPyPz(momcasc_new);
  delete trackCasc;

	if(fDoEventMixing){
    Int_t nextRes( nextResVec[fPoolIndex] );
		TLorentzVector *lv = new TLorentzVector();
		lv->SetXYZM(momcasc_new[0],momcasc_new[1],momcasc_new[2],casc->MassXi());
    Double_t xyzR125pr[3]={9999.,9999.,9999.};
    Double_t xyzR125pi[3]={9999.,9999.,9999.};
    Double_t xyzR125bach[3]={9999.,9999.,9999.};
		if(casc->ChargeXi()>0){
      m_ReservoirL1[fPoolIndex][nextRes].push_back(lv);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.){
        fAnalCuts->SetSftPosR125(cptrack,fBzkG,posVtx,xyzR125pr);
        fAnalCuts->SetSftPosR125(cntrack,fBzkG,posVtx,xyzR125pi);
        fAnalCuts->SetSftPosR125(cbtrack,fBzkG,posVtx,xyzR125bach);
      }
      TVector *varvec = new TVector(9);
      (*varvec)[0] = xyzR125pr[0];
      (*varvec)[1] = xyzR125pr[1];
      (*varvec)[2] = xyzR125pr[2];
      (*varvec)[3] = xyzR125pi[0];
      (*varvec)[4] = xyzR125pi[1];
      (*varvec)[5] = xyzR125pi[2];
      (*varvec)[6] = xyzR125bach[0];
      (*varvec)[7] = xyzR125bach[1];
      (*varvec)[8] = xyzR125bach[2];
      m_ReservoirVarsL1[fPoolIndex][nextRes].push_back(varvec);
    }else{
      m_ReservoirL2[fPoolIndex][nextRes].push_back(lv);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.){
        fAnalCuts->SetSftPosR125(cntrack,fBzkG,posVtx,xyzR125pr);
        fAnalCuts->SetSftPosR125(cptrack,fBzkG,posVtx,xyzR125pi);
        fAnalCuts->SetSftPosR125(cbtrack,fBzkG,posVtx,xyzR125bach);
      }
      TVector *varvec = new TVector(9);
      (*varvec)[0] = xyzR125pr[0];
      (*varvec)[1] = xyzR125pr[1];
      (*varvec)[2] = xyzR125pr[2];
      (*varvec)[3] = xyzR125pi[0];
      (*varvec)[4] = xyzR125pi[1];
      (*varvec)[5] = xyzR125pi[2];
      (*varvec)[6] = xyzR125bach[0];
      (*varvec)[7] = xyzR125bach[1];
      (*varvec)[8] = xyzR125bach[2];
      m_ReservoirVarsL2[fPoolIndex][nextRes].push_back(varvec);
    }
	}

	if(!fWriteEachVariableTree) return;

	for(Int_t i=0;i<26;i++){
		fCandidateCascVariables[i] = -9999.;
	}
	for(Int_t i=0;i<16;i++){
		fCandidateSingleVariables[i] = -9999.;
	}
  fCandidateSingleVariables[ 0] = momcasc_new[0];
  fCandidateSingleVariables[ 1] = momcasc_new[1];
  fCandidateSingleVariables[ 2] = momcasc_new[2];
  fCandidateSingleVariables[ 3] = casc->ChargeXi();
  fCandidateSingleVariables[ 4] = casc->MassXi();
  fCandidateSingleVariables[ 5] = fBzkG;
  fCandidateSingleVariables[ 6] = fCentrality;
  fCandidateSingleVariables[ 7] = fVtxZ;
  fCandidateSingleVariables[ 8] = fEvNumberCounter;
  fCandidateSingleVariables[ 9] = fRunNumber;
  if(fUseEventPlane==4){
    //Double_t evplane = GetEventPlaneForCandidate(0,v0,event->GetEventplane(),fQ,fQSub1,fQSub2); // remove autocorrelations
    //fCandidateSingleVariables[10] = evplane;
    fCandidateSingleVariables[10] = fEventPlane;

  }else{
    fCandidateSingleVariables[10] = fEventPlane;
  }
  fCandidateSingleVariables[11] = casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);
  fCandidateSingleVariables[12] = cbtrack->GetID();
  fCandidateSingleVariables[13] = cptrack->GetID();
  fCandidateSingleVariables[14] = cntrack->GetID();
  fCandidateSingleVariables[15] = casc->DecayLengthV0();

  fSingleVariablesTree->Fill();

//  fCandidateCascVariables[ 0] = fCentrality;
//  fCandidateCascVariables[ 1] = casc->MassXi();
//  fCandidateCascVariables[ 2] = momcasc_new[0];//casc->MomXiX();
//  fCandidateCascVariables[ 3] = momcasc_new[1];//casc->MomXiY();
//  fCandidateCascVariables[ 4] = momcasc_new[2];//casc->MomXiZ();
//	if(casc->ChargeXi()<0)
//		fCandidateCascVariables[ 5] = casc->MassLambda();
//	else
//		fCandidateCascVariables[ 5] = casc->MassAntiLambda();
//
//  fCandidateCascVariables[ 6] = casc->DcaXiDaughters();
//  fCandidateCascVariables[ 7] = casc->DcaV0Daughters();
//  fCandidateCascVariables[ 8] = casc->DecayLengthXi(posVtx[0],posVtx[1],posVtx[2]);
//  fCandidateCascVariables[ 9] = casc->CosPointingAngleXi(posVtx[0],posVtx[1],posVtx[2]);
//  fCandidateCascVariables[10] = casc->DcaV0ToPrimVertex();
//  fCandidateCascVariables[11] = casc->DcaPosToPrimVertex();
//  fCandidateCascVariables[12] = casc->DcaNegToPrimVertex();
//  fCandidateCascVariables[13] = casc->DcaBachToPrimVertex();
//  fCandidateCascVariables[14] = casc->DecayLengthV0();
//  fCandidateCascVariables[15] = casc->CosPointingAngle(casc->GetDecayVertexXi());
//  fCandidateCascVariables[16] = casc->ChargeXi();
//  fCandidateCascVariables[17] = xipdgcode;
//  fCandidateCascVariables[18] = ximotherpdgcode;
//  fCandidateCascVariables[19] = mcxipx;
//  fCandidateCascVariables[20] = mcxipy;
//  fCandidateCascVariables[21] = mcxipz;
//  fCandidateCascVariables[22] = labcasc;
//  fCandidateCascVariables[23] = fRunNumber;
//  fCandidateCascVariables[24] = fVtxZ;
//  fCandidateCascVariables[25] = fEvNumberCounter;
//
//
//	fCascVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineMCTreeVariables() 
{
  //
  // Define electron tree variables
  //

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fMCVariablesTree = new TTree(nameoutput,"MC variables tree");
  Int_t nVar = 16;
  fCandidateMCVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="DecayType";
  fCandidateVariableNames[ 2]="XicPx";
  fCandidateVariableNames[ 3]="XicPy";
  fCandidateVariableNames[ 4]="XicPz";
  fCandidateVariableNames[ 5]="ElePx";
  fCandidateVariableNames[ 6]="ElePy";
  fCandidateVariableNames[ 7]="ElePz";
  fCandidateVariableNames[ 8]="CascPx";
  fCandidateVariableNames[ 9]="CascPy";
  fCandidateVariableNames[10]="CascPz";
  fCandidateVariableNames[11]="PdgCode";
  fCandidateVariableNames[12]="ElePdgCode";
  fCandidateVariableNames[13]="CascPdgCode";
  fCandidateVariableNames[14]="RunNumber";
  fCandidateVariableNames[15]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mccascpart, Int_t decaytype) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mcpart) return;
	if(!mcepart) return;
	if(!mccascpart) return;

	for(Int_t i=0;i<14;i++){
		fCandidateMCVariables[i] = -9999.;
	}

	fCandidateMCVariables[ 0] = fCentrality;
	fCandidateMCVariables[ 1] = decaytype;
	fCandidateMCVariables[ 2] = mcpart->Px();
	fCandidateMCVariables[ 3] = mcpart->Py();
	fCandidateMCVariables[ 4] = mcpart->Pz();
	fCandidateMCVariables[ 5] = mcepart->Px();
	fCandidateMCVariables[ 6] = mcepart->Py();
	fCandidateMCVariables[ 7] = mcepart->Pz();
	fCandidateMCVariables[ 8] = mccascpart->Px();
	fCandidateMCVariables[ 9] = mccascpart->Py();
	fCandidateMCVariables[10] = mccascpart->Pz();
	fCandidateMCVariables[11] = mcpart->GetPdgCode();
	fCandidateMCVariables[12] = mcepart->GetPdgCode();
	fCandidateMCVariables[13] = mccascpart->GetPdgCode();
	fCandidateMCVariables[14] = fRunNumber;
	fCandidateMCVariables[15] = fEvNumberCounter;

	Double_t epx = mcepart->Px();
	Double_t epy = mcepart->Py();
	Double_t epz = mcepart->Pz();
	Double_t eE = sqrt(epx*epx+epy*epy+epz*epz+0.000511*0.000511);
	Double_t cascpx = mccascpart->Px();
	Double_t cascpy = mccascpart->Py();
	Double_t cascpz = mccascpart->Pz();
	Double_t cascE = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+1.32171*1.32171);

	Double_t InvMassEleXi = sqrt(pow(eE+cascE,2)-pow(epx+cascpx,2)-pow(epy+cascpy,2)-pow(epz+cascpz,2));

	Double_t cont[4];
	cont[0] = InvMassEleXi;
	cont[1] = mcpart->Pt();
	cont[2] = 0.0;
	cont[3] = fCentrality;
	Double_t cont2[3];
	cont2[0] = InvMassEleXi;
	cont2[1] = mcepart->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = mcepart->Pt();
	cont_eleptvseta[1] = mcepart->Eta();
	cont_eleptvseta[2] = fCentrality;
	Double_t cont_eleptvsxipt[3];
	cont_eleptvsxipt[0] = mcepart->Pt();
	cont_eleptvsxipt[1] = mccascpart->Pt();
	cont_eleptvsxipt[2] = fCentrality;
	Double_t cont_eleptvsxiptvsxicpt[4];
	cont_eleptvsxiptvsxicpt[0] = mcepart->Pt();
	cont_eleptvsxiptvsxicpt[1] = mccascpart->Pt();
	cont_eleptvsxiptvsxicpt[2] = mcpart->Pt();
	cont_eleptvsxiptvsxicpt[3] = fCentrality;

	Double_t contmc[3];
	contmc[0] = mcpart->Pt();
	contmc[1] = mcpart->Y();
	contmc[2] = fCentrality;
	Double_t contmcele[3];
	contmcele[0] = mcepart->Pt();
	contmcele[1] = mcepart->Eta();
	contmcele[2] = fCentrality;


	AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
	Float_t etamin, etamax;
	esdcuts->GetEtaRange(etamin,etamax);

	if(decaytype==0){
		fHistoXicMCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoXicMCGen1->Fill(contmc);//4132 is particle
		if(mcpart->GetPdgCode()<0) fHistoXicMCGen2->Fill(contmc);//-4132 is anti-particle
		fHistoXicElectronMCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoXicElectronMCGen1->Fill(contmcele);//-11 is positron
		if(mcepart->GetPdgCode()>0) fHistoXicElectronMCGen2->Fill(contmcele);//11 is electron
		fHistoEleXiMassMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleXiMassvsElePtMCGen->Fill(cont2);
			if(mcepart->GetPdgCode()<0) fHistoEleXiMassvsElePtMCGen1->Fill(cont2);//-11 is positron
			if(mcepart->GetPdgCode()>0) fHistoEleXiMassvsElePtMCGen2->Fill(cont2);//11 is electron
			if(InvMassEleXi<fAnalCuts->GetEleXiMassMax()){
				fHistoElePtMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsXiPtMCGen->Fill(cont_eleptvsxipt);
			}
		}
		if(fabs(mcpart->Y())<0.7){
			if(InvMassEleXi<fAnalCuts->GetEleXiMassMax()){
				fHistoElePtvsXiPtMCXicGen->Fill(cont_eleptvsxipt);
				fHistoElePtvsXiPtvsXicPtMCGen->Fill(cont_eleptvsxiptvsxicpt);
			}
		}
	}else if(decaytype==10){
		fHistoXibMCGen->Fill(contmc);
	}

	if(fWriteMCVariableTree)
		fMCVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineMCEleTreeVariables() 
{
  //
  // Define mc gen electron tree variables
  //

  const char* nameoutput = GetOutputSlot(9)->GetContainer()->GetName();
  fMCEleVariablesTree = new TTree(nameoutput,"MC Ele variables tree");
  Int_t nVar = 8;
  fCandidateMCEleVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="ElePx";
  fCandidateVariableNames[ 2]="ElePy";
  fCandidateVariableNames[ 3]="ElePz";
  fCandidateVariableNames[ 4]="ElePdgCode";
  fCandidateVariableNames[ 5]="EleMotherPdgCode";
  fCandidateVariableNames[ 6]="RunNumber";
  fCandidateVariableNames[ 7]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCEleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCEleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMCEleROOTObjects(AliAODMCParticle *mcepart, TClonesArray *mcArray) 
{
  //
  // Fill tree depending on fWriteMCVariableTree 
  //
  if(!mcepart) return;

  Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;
  GetMCDecayHistory(mcepart,mcArray,pdgarray_ele,labelarray_ele,ngen_ele);
  Bool_t ele_from_bottom = HaveBottomInHistory(pdgarray_ele);
  Bool_t ele_from_charm = HaveCharmInHistory(pdgarray_ele);
  Int_t semi_flag = FromSemileptonicDecays(pdgarray_ele);

  Double_t contmc[3];
  contmc[0] = mcepart->Pt();
  contmc[1] = mcepart->Eta();
  contmc[2] = fCentrality;

  if(semi_flag==1 && !ele_from_bottom){
    fHistoCharmElectronMCGen->Fill(contmc);
  }
  if(semi_flag==1 && ele_from_bottom){
    fHistoBottomElectronMCGen->Fill(contmc);
  }
  if(semi_flag==2 ){
    fHistoBottomElectronMCGen->Fill(contmc);
  }


  Bool_t hfe_flag = kFALSE;
  Int_t labemother = mcepart->GetMother();
  Int_t pdgmotherele = -9999;
  if(labemother>=0){
    AliAODMCParticle *motherele = (AliAODMCParticle*)mcArray->At(labemother);
    pdgmotherele = motherele->GetPdgCode();
    if(abs(pdgmotherele)>4000&&abs(pdgmotherele)<4400){
      hfe_flag = kTRUE;
    }
  }
  if(!hfe_flag) return;

  fHistoElectronMCGen->Fill(contmc);

  for(Int_t i=0;i<8;i++){
    fCandidateMCEleVariables[i] = -9999.;
  }

	fCandidateMCEleVariables[ 0] = fCentrality;
	fCandidateMCEleVariables[ 1] = mcepart->Px();
	fCandidateMCEleVariables[ 2] = mcepart->Py();
	fCandidateMCEleVariables[ 3] = mcepart->Pz();
	fCandidateMCEleVariables[ 4] = mcepart->GetPdgCode();
	fCandidateMCEleVariables[ 5] = pdgmotherele;
	fCandidateMCEleVariables[ 6] = fRunNumber;
	fCandidateMCEleVariables[ 7] = fEvNumberCounter;

	// This function makes output too heavy  (do not use if you have output size limitation)
	//if(fWriteMCVariableTree && fWriteEachVariableTree && mcepart->Pt()>0.4 && fabs(mcepart->Eta())<1.0)
		//fMCEleVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineMCCascTreeVariables() 
{
  //
  // Define Mc cascade tree variables
  //

  const char* nameoutput = GetOutputSlot(10)->GetContainer()->GetName();
  fMCCascVariablesTree = new TTree(nameoutput,"MC cascade variables tree");
  Int_t nVar = 8;
  fCandidateMCCascVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="CascPx";
  fCandidateVariableNames[ 2]="CascPy";
  fCandidateVariableNames[ 3]="CascPz";
  fCandidateVariableNames[ 4]="CascPdgCode";
  fCandidateVariableNames[ 5]="CascMotherPdgCode";
  fCandidateVariableNames[ 6]="RunNumber";
  fCandidateVariableNames[ 7]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCCascVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCCascVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMCCascROOTObjects(AliAODMCParticle *mccascpart, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mccascpart) return;

	for(Int_t i=0;i<8;i++){
		fCandidateMCCascVariables[i] = -9999.;
	}

	Bool_t hfxi_flag = kFALSE;
	Int_t labcascmother = mccascpart->GetMother();
	Int_t pdgmothercasc = -9999;
	if(labcascmother>=0){
		AliAODMCParticle *mothercasc = (AliAODMCParticle*)mcArray->At(labcascmother);
		if(mothercasc){
			pdgmothercasc = mothercasc->GetPdgCode();
			if(abs(pdgmothercasc)>4000&&abs(pdgmothercasc)<4400){
				hfxi_flag = kTRUE;
			}
		}
	}
	if(!hfxi_flag) return;

	Double_t contmc[3];
	contmc[0] = mccascpart->Pt();
	contmc[1] = mccascpart->Eta();
	contmc[2] = fCentrality;
	fHistoXiMCGen->Fill(contmc);

	fCandidateMCCascVariables[ 0] = fCentrality;
	fCandidateMCCascVariables[ 1] = mccascpart->Px();
	fCandidateMCCascVariables[ 2] = mccascpart->Py();
	fCandidateMCCascVariables[ 3] = mccascpart->Pz();
	fCandidateMCCascVariables[ 4] = mccascpart->GetPdgCode();
	fCandidateMCCascVariables[ 5] = pdgmothercasc;
	fCandidateMCCascVariables[ 6] = fRunNumber;
	fCandidateMCCascVariables[ 7] = fEvNumberCounter;

	// This function makes output too heavy  (do not use if you have output size limitation)
	//if(fWriteMCVariableTree && fWriteEachVariableTree && mccascpart->Pt()>0.4 && fabs(mccascpart->Eta())<1.0)
		//fMCCascVariablesTree->Fill();
}

////__________________________________________________________________________
void  AliAnalysisTaskSEXic2eleXifromAODtracks::DefineGeneralHistograms() {
  //
  // This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",18,-0.5,17.5);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetBinLabel(1,"X1");
  fCEvents->GetXaxis()->SetBinLabel(2,"Analyzed events");
  fCEvents->GetXaxis()->SetBinLabel(3,"AliAODVertex exists");
  fCEvents->GetXaxis()->SetBinLabel(4,"TriggerOK");
  fCEvents->GetXaxis()->SetBinLabel(5,"IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(6,"CascadesHF exists");
  fCEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fCEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fCEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fCEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fCEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fCEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fCEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fCEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  fCEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fCEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fCEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnalCuts->GetMaxVtxZ()));
  //fCEvents->GetXaxis()->SetTitle("");
  fCEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger","counter",18,-0.5,17.5);
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fHCentrality = new TH1F("fHCentrality","conter",100,0.,100.);
  fHEventPlane = new TH1F("fHEventPlane","conter",100,-3.14,3.14);
	fHNTrackletvsZ = new TH2F("fHNTrackletvsZ","N_{tracklet} vs z",30,-15.,15.,120,-0.5,119.5);
	fHNTrackletCorrvsZ = new TH2F("fHNTrackletCorrvsZ","N_{tracklet} vs z",30,-15.,15.,120,-0.5,119.5);

  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);
  fOutput->Add(fHEventPlane);
	fOutput->Add(fHNTrackletvsZ);
	fOutput->Add(fHNTrackletCorrvsZ);

  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSEXic2eleXifromAODtracks::DefineAnalysisHistograms() 
{
  //
  // Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[4]=		{22	,100	, 100   ,10};
  Double_t xmin_base[4]={1.3,0		, -0.5 ,0.00};
  Double_t xmax_base[4]={5.7,20.	, 0.5  ,100};

  fHistoEleXiMass = new THnSparseF("fHistoEleXiMass","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMass);
  fHistoEleXiMassRS = new THnSparseF("fHistoEleXiMassRS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS);
  fHistoEleXiMassWS = new THnSparseF("fHistoEleXiMassWS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS);
  fHistoEleXiMassRSMix = new THnSparseF("fHistoEleXiMassRSMix","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix);
  fHistoEleXiMassWSMix = new THnSparseF("fHistoEleXiMassWSMix","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix);
  fHistoEleXiMassRSSide = new THnSparseF("fHistoEleXiMassRSSide","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSide);
  fHistoEleXiMassWSSide = new THnSparseF("fHistoEleXiMassWSSide","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSide);
  fHistoEleXiMassRS1 = new THnSparseF("fHistoEleXiMassRS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS1);
  fHistoEleXiMassWS1 = new THnSparseF("fHistoEleXiMassWS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS1);
  fHistoEleXiMassRSMix1 = new THnSparseF("fHistoEleXiMassRSMix1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix1);
  fHistoEleXiMassWSMix1 = new THnSparseF("fHistoEleXiMassWSMix1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix1);
  fHistoEleXiMassRSSide1 = new THnSparseF("fHistoEleXiMassRSSide1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSide1);
  fHistoEleXiMassWSSide1 = new THnSparseF("fHistoEleXiMassWSSide1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSide1);
  fHistoEleXiMassRS2 = new THnSparseF("fHistoEleXiMassRS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS2);
  fHistoEleXiMassWS2 = new THnSparseF("fHistoEleXiMassWS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS2);
  fHistoEleXiMassRSMix2 = new THnSparseF("fHistoEleXiMassRSMix2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix2);
  fHistoEleXiMassWSMix2 = new THnSparseF("fHistoEleXiMassWSMix2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix2);
  fHistoEleXiMassRSSide2 = new THnSparseF("fHistoEleXiMassRSSide2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSide2);
  fHistoEleXiMassWSSide2 = new THnSparseF("fHistoEleXiMassWSSide2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSide2);

  fHistoEleXiMassRSAway = new THnSparseF("fHistoEleXiMassRSAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSAway);
  fHistoEleXiMassWSAway = new THnSparseF("fHistoEleXiMassWSAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSAway);
  fHistoEleXiMassRSMixAway = new THnSparseF("fHistoEleXiMassRSMixAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMixAway);
  fHistoEleXiMassWSMixAway = new THnSparseF("fHistoEleXiMassWSMixAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMixAway);
  fHistoEleXiMassRSSideAway = new THnSparseF("fHistoEleXiMassRSSideAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSideAway);
  fHistoEleXiMassWSSideAway = new THnSparseF("fHistoEleXiMassWSSideAway","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSideAway);
  fHistoEleXiMassRS1Away = new THnSparseF("fHistoEleXiMassRS1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS1Away);
  fHistoEleXiMassWS1Away = new THnSparseF("fHistoEleXiMassWS1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS1Away);
  fHistoEleXiMassRSMix1Away = new THnSparseF("fHistoEleXiMassRSMix1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix1Away);
  fHistoEleXiMassWSMix1Away = new THnSparseF("fHistoEleXiMassWSMix1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix1Away);
  fHistoEleXiMassRSSide1Away = new THnSparseF("fHistoEleXiMassRSSide1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSide1Away);
  fHistoEleXiMassWSSide1Away = new THnSparseF("fHistoEleXiMassWSSide1Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSide1Away);
  fHistoEleXiMassRS2Away = new THnSparseF("fHistoEleXiMassRS2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRS2Away);
  fHistoEleXiMassWS2Away = new THnSparseF("fHistoEleXiMassWS2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWS2Away);
  fHistoEleXiMassRSMix2Away = new THnSparseF("fHistoEleXiMassRSMix2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSMix2Away);
  fHistoEleXiMassWSMix2Away = new THnSparseF("fHistoEleXiMassWSMix2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSMix2Away);
  fHistoEleXiMassRSSide2Away = new THnSparseF("fHistoEleXiMassRSSide2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassRSSide2Away);
  fHistoEleXiMassWSSide2Away = new THnSparseF("fHistoEleXiMassWSSide2Away","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassWSSide2Away);

  Int_t bins_base_elept[3]=		{10	,100		,10};
  Double_t xmin_base_elept[3]={1.3,0		,0.00};
  Double_t xmax_base_elept[3]={3.3,10.	,100};

  fHistoEleXiMassvsElePtRS = new THnSparseF("fHistoEleXiMassvsElePtRS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRS);
  fHistoEleXiMassvsElePtWS = new THnSparseF("fHistoEleXiMassvsElePtWS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWS);
  fHistoEleXiMassvsElePtRSMix = new THnSparseF("fHistoEleXiMassvsElePtRSMix","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSMix);
  fHistoEleXiMassvsElePtWSMix = new THnSparseF("fHistoEleXiMassvsElePtWSMix","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSMix);
  fHistoEleXiMassvsElePtRSSide = new THnSparseF("fHistoEleXiMassvsElePtRSSide","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSSide);
  fHistoEleXiMassvsElePtWSSide = new THnSparseF("fHistoEleXiMassvsElePtWSSide","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSSide);
  fHistoEleXiMassvsElePtRS1 = new THnSparseF("fHistoEleXiMassvsElePtRS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRS1);
  fHistoEleXiMassvsElePtWS1 = new THnSparseF("fHistoEleXiMassvsElePtWS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWS1);
  fHistoEleXiMassvsElePtRSMix1 = new THnSparseF("fHistoEleXiMassvsElePtRSMix1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSMix1);
  fHistoEleXiMassvsElePtWSMix1 = new THnSparseF("fHistoEleXiMassvsElePtWSMix1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSMix1);
  fHistoEleXiMassvsElePtRSSide1 = new THnSparseF("fHistoEleXiMassvsElePtRSSide1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSSide1);
  fHistoEleXiMassvsElePtWSSide1 = new THnSparseF("fHistoEleXiMassvsElePtWSSide1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSSide1);
  fHistoEleXiMassvsElePtRS2 = new THnSparseF("fHistoEleXiMassvsElePtRS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRS2);
  fHistoEleXiMassvsElePtWS2 = new THnSparseF("fHistoEleXiMassvsElePtWS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWS2);
  fHistoEleXiMassvsElePtRSMix2 = new THnSparseF("fHistoEleXiMassvsElePtRSMix2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSMix2);
  fHistoEleXiMassvsElePtWSMix2 = new THnSparseF("fHistoEleXiMassvsElePtWSMix2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSMix2);
  fHistoEleXiMassvsElePtRSSide2 = new THnSparseF("fHistoEleXiMassvsElePtRSSide2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtRSSide2);
  fHistoEleXiMassvsElePtWSSide2 = new THnSparseF("fHistoEleXiMassvsElePtWSSide2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtWSSide2);

  fHistoElePtRS=new TH2F("fHistoElePtRS","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRS);
  fHistoElePtWS=new TH2F("fHistoElePtWS","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWS);
  fHistoElePtRSMix=new TH2F("fHistoElePtRSMix","Right-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtRSMix);
  fHistoElePtWSMix=new TH2F("fHistoElePtWSMix","Wrong-sign e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtWSMix);

  fHistoEleXiMassMCS = new THnSparseF("fHistoEleXiMassMCS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCS);
  fHistoEleXiMassMCS1 = new THnSparseF("fHistoEleXiMassMCS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCS1);
  fHistoEleXiMassMCS2 = new THnSparseF("fHistoEleXiMassMCS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCS2);
  fHistoEleXiMassXibMCS = new THnSparseF("fHistoEleXiMassXibMCS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassXibMCS);
  fHistoEleXiMassXibMCS1 = new THnSparseF("fHistoEleXiMassXibMCS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassXibMCS1);
  fHistoEleXiMassXibMCS2 = new THnSparseF("fHistoEleXiMassXibMCS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassXibMCS2);
  fHistoEleXiMassPromptMCS = new THnSparseF("fHistoEleXiMassPromptMCS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassPromptMCS);
  fHistoEleXiMassPromptMCS1 = new THnSparseF("fHistoEleXiMassPromptMCS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassPromptMCS1);
  fHistoEleXiMassPromptMCS2 = new THnSparseF("fHistoEleXiMassPromptMCS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassPromptMCS2);
  fHistoEleXiMassBFeeddownMCS = new THnSparseF("fHistoEleXiMassBFeeddownMCS","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassBFeeddownMCS);
  fHistoEleXiMassBFeeddownMCS1 = new THnSparseF("fHistoEleXiMassBFeeddownMCS1","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassBFeeddownMCS1);
  fHistoEleXiMassBFeeddownMCS2 = new THnSparseF("fHistoEleXiMassBFeeddownMCS2","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassBFeeddownMCS2);
  fHistoEleXiMassMCGen = new THnSparseF("fHistoEleXiMassMCGen","",4,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleXiMassMCGen);
  fHistoEleXiMassvsElePtMCS = new THnSparseF("fHistoEleXiMassvsElePtMCS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCS);
  fHistoEleXiMassvsElePtMCGen = new THnSparseF("fHistoEleXiMassvsElePtMCGen","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCGen);
  fHistoEleXiMassvsElePtMCS1 = new THnSparseF("fHistoEleXiMassvsElePtMCS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCS1);
  fHistoEleXiMassvsElePtMCGen1 = new THnSparseF("fHistoEleXiMassvsElePtMCGen1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCGen1);
  fHistoEleXiMassvsElePtMCS2 = new THnSparseF("fHistoEleXiMassvsElePtMCS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCS2);
  fHistoEleXiMassvsElePtMCGen2 = new THnSparseF("fHistoEleXiMassvsElePtMCGen2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleXiMassvsElePtMCGen2);
  fHistoElePtMCS=new TH2F("fHistoElePtMCS","MC S e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtMCS);
  fHistoElePtMCGen=new TH2F("fHistoElePtMCGen","MC Gen e pt",100,0.,10.,10,0.,100.);
  fOutputAll->Add(fHistoElePtMCGen);

  Int_t bins_eleptvseta[3]=		{50,20	,10};
  Double_t xmin_eleptvseta[3]={0.,-1.	,0.0};
  Double_t xmax_eleptvseta[3]={5.,1.	,100};

  fHistoElePtvsEtaRS = new THnSparseF("fHistoElePtvsEtaRS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaRS);
  fHistoElePtvsEtaWS = new THnSparseF("fHistoElePtvsEtaWS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaWS);
  fHistoElePtvsEtaRSMix = new THnSparseF("fHistoElePtvsEtaRSMix","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaRSMix);
  fHistoElePtvsEtaWSMix = new THnSparseF("fHistoElePtvsEtaWSMix","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaWSMix);
  fHistoElePtvsEtaMCS = new THnSparseF("fHistoElePtvsEtaMCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaMCS);
  fHistoElePtvsEtaMCGen = new THnSparseF("fHistoElePtvsEtaMCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaMCGen);

  Int_t bins_eleptvsxipt[3]=	{50,20	,10};
  Double_t xmin_eleptvsxipt[3]={0.,0.	,0.0};
  Double_t xmax_eleptvsxipt[3]={5.,5.	,100};

  fHistoElePtvsXiPtRS = new THnSparseF("fHistoElePtvsXiPtRS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtRS);
  fHistoElePtvsXiPtWS = new THnSparseF("fHistoElePtvsXiPtWS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtWS);
  fHistoElePtvsXiPtRSMix = new THnSparseF("fHistoElePtvsXiPtRSMix","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtRSMix);
  fHistoElePtvsXiPtWSMix = new THnSparseF("fHistoElePtvsXiPtWSMix","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtWSMix);
  fHistoElePtvsXiPtMCS = new THnSparseF("fHistoElePtvsXiPtMCS","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCS);
  fHistoElePtvsXiPtMCGen = new THnSparseF("fHistoElePtvsXiPtMCGen","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCGen);
  fHistoElePtvsXiPtMCXicGen = new THnSparseF("fHistoElePtvsXiPtMCXicGen","",3,bins_eleptvsxipt,xmin_eleptvsxipt,xmax_eleptvsxipt);
  fOutputAll->Add(fHistoElePtvsXiPtMCXicGen);

  Int_t bins_eleptvsxiptvsxicpt[4]=	{50,20,10,10};
  Double_t xmin_eleptvsxiptvsxicpt[4]={0.,0.,0.,0.0};
  Double_t xmax_eleptvsxiptvsxicpt[4]={5.,5.,10.,100};
  fHistoElePtvsXiPtvsXicPtMCS = new THnSparseF("fHistoElePtvsXiPtvsXicPtMCS","",4,bins_eleptvsxiptvsxicpt,xmin_eleptvsxiptvsxicpt,xmax_eleptvsxiptvsxicpt);
  fOutputAll->Add(fHistoElePtvsXiPtvsXicPtMCS);
  fHistoElePtvsXiPtvsXicPtMCGen = new THnSparseF("fHistoElePtvsXiPtvsXicPtMCGen","",4,bins_eleptvsxiptvsxicpt,xmin_eleptvsxiptvsxicpt,xmax_eleptvsxiptvsxicpt);
  fOutputAll->Add(fHistoElePtvsXiPtvsXicPtMCGen);

  Int_t bins_eleptvsd0[3]=	{50 ,50	,10};
  Double_t xmin_eleptvsd0[3]={0.,-0.2	,0.0};
  Double_t xmax_eleptvsd0[3]={5.,0.2	,100};

  fHistoElePtvsd0RS = new THnSparseF("fHistoElePtvsd0RS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0RS);
  fHistoElePtvsd0WS = new THnSparseF("fHistoElePtvsd0WS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0WS);
  fHistoElePtvsd0RSMix = new THnSparseF("fHistoElePtvsd0RSMix","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0RSMix);
  fHistoElePtvsd0WSMix = new THnSparseF("fHistoElePtvsd0WSMix","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0WSMix);
  fHistoElePtvsd0MCS = new THnSparseF("fHistoElePtvsd0MCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0MCS);
  fHistoElePtvsd0PromptMCS = new THnSparseF("fHistoElePtvsd0PromptMCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0PromptMCS);
  fHistoElePtvsd0BFeeddownMCS = new THnSparseF("fHistoElePtvsd0BFeeddownMCS","",3,bins_eleptvsd0,xmin_eleptvsd0,xmax_eleptvsd0);
  fOutputAll->Add(fHistoElePtvsd0BFeeddownMCS);

  //------------------------------------------------
  // checking histograms
  //------------------------------------------------
  fHistoBachPt = new TH1F("fHistoBachPt","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPt);
  fHistoBachPtMCS = new TH1F("fHistoBachPtMCS","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCS);
  fHistoBachPtMCGen = new TH1F("fHistoBachPtMCGen","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCGen);
  fHistod0Bach = new TH1F("fHistod0Bach","Bachelor d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0Bach);
  fHistoXiMassvsPt=new TH2F("fHistoXiMassvsPt","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPt);
  fHistoXiMassvsPtMCS=new TH2F("fHistoXiMassvsPtMCS","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPtMCS);
  fHistoXiMassvsPtMCGen=new TH2F("fHistoXiMassvsPtMCGen","Xi mass",100,1.32-0.05,1.32+0.05,20,0.,10.);
  fOutputAll->Add(fHistoXiMassvsPtMCGen);
  fHistoOmegaMassvsPt=new TH2F("fHistoOmegaMassvsPt","Omega mass",100,1.67-0.05,1.67+0.05,20,0.,10.);
  fOutputAll->Add(fHistoOmegaMassvsPt);

  fHistoElectronTPCPID=new TH2F("fHistoElectronTPCPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCPID);
  fHistoElectronTOFPID=new TH2F("fHistoElectronTOFPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFPID);
  fHistoElectronTPCSelPID=new TH2F("fHistoElectronTPCSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCSelPID);
  fHistoElectronTOFSelPID=new TH2F("fHistoElectronTOFSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFSelPID);
  fHistoElectronTPCPIDSelTOF=new TH2F("fHistoElectronTPCPIDSelTOF","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOF);
  fHistoElectronTOFPIDSelTPC=new TH2F("fHistoElectronTOFPIDSelTPC","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTOFPIDSelTPC);
  fHistoElectronTPCPIDSelTOFSmallEta=new TH2F("fHistoElectronTPCPIDSelTOFSmallEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFSmallEta);
  fHistoElectronTPCPIDSelTOFLargeEta=new TH2F("fHistoElectronTPCPIDSelTOFLargeEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFLargeEta);
  fHistoMassConversionsMin=new TH1F("fHistoMassConversionsMin","",500,0,0.5);
  fOutputAll->Add(fHistoMassConversionsMin);
  fHistoMassConversionsSameSignMin=new TH1F("fHistoMassConversionsSameSignMin","",500,0,0.5);
  fOutputAll->Add(fHistoMassConversionsSameSignMin);

	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i]=new TH2F(Form("fHistoElectronTPCPIDSelTOFEtaDep[%d]",i),"",10,0.,5.,500,-10.,10.);
		fOutputAll->Add(fHistoElectronTPCPIDSelTOFEtaDep[i]);
	}

  fHistoElectronQovPtvsPhi=new TH2F("fHistoElectronQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoElectronQovPtvsPhi);
  fHistoXiQovPtvsPhi=new TH2F("fHistoXiQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoXiQovPtvsPhi);

  Int_t bins_xicmcgen[3]=	{100 ,20	,10};
  Double_t xmin_xicmcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_xicmcgen[3]={20.,1.0	,100};
  fHistoXicMCGen = new THnSparseF("fHistoXicMCGen","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCGen);
  fHistoXicMCGen1 = new THnSparseF("fHistoXicMCGen1","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCGen1);
  fHistoXicMCGen2 = new THnSparseF("fHistoXicMCGen2","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCGen2);
  fHistoXicMCS = new THnSparseF("fHistoXicMCS","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCS);
  fHistoXicMCS1 = new THnSparseF("fHistoXicMCS1","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCS1);
  fHistoXicMCS2 = new THnSparseF("fHistoXicMCS2","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fOutputAll->Add(fHistoXicMCS2);

  Int_t bins_xibmcgen[3]=	{100 ,20	,10};
  Double_t xmin_xibmcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_xibmcgen[3]={50.,1.0	,100};
  fHistoXibMCGen = new THnSparseF("fHistoXibMCGen","",3,bins_xibmcgen,xmin_xibmcgen,xmax_xibmcgen);
  fOutputAll->Add(fHistoXibMCGen);
  fHistoXibMCS = new THnSparseF("fHistoXibMCS","",3,bins_xibmcgen,xmin_xibmcgen,xmax_xibmcgen);
  fOutputAll->Add(fHistoXibMCS);

  Int_t bins_xibmcgen_withxic[3]=	{50 ,100	,100};
  Double_t xmin_xibmcgen_withxic[3]={0.,-5.,-5.};
  Double_t xmax_xibmcgen_withxic[3]={50.,5.,5.};
  fHistoXibMCGenWithXic = new THnSparseF("fHistoXibMCGenWithXic","",3,bins_xibmcgen_withxic,xmin_xibmcgen_withxic,xmax_xibmcgen_withxic);
  fOutputAll->Add(fHistoXibMCGenWithXic);

  Int_t bins_elemcgen[3]=	{100 ,20	,10};
  Double_t xmin_elemcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_elemcgen[3]={10.,1.0	,100};
  fHistoElectronMCGen = new THnSparseF("fHistoElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronMCGen);
  fHistoBottomElectronMCGen = new THnSparseF("fHistoBottomElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoBottomElectronMCGen);
  fHistoCharmElectronMCGen = new THnSparseF("fHistoCharmElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoCharmElectronMCGen);
  fHistoXicElectronMCGen = new THnSparseF("fHistoXicElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCGen);
  fHistoXicElectronMCGen1 = new THnSparseF("fHistoXicElectronMCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCGen1);
  fHistoXicElectronMCGen2 = new THnSparseF("fHistoXicElectronMCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCGen2);
  fHistoXicElectronMCS = new THnSparseF("fHistoXicElectronMCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCS);
  fHistoXicElectronMCS1 = new THnSparseF("fHistoXicElectronMCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCS1);
  fHistoXicElectronMCS2 = new THnSparseF("fHistoXicElectronMCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoXicElectronMCS2);

  Int_t bins_ximcgen[3]=	{50 ,20	,10};
  Double_t xmin_ximcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_ximcgen[3]={10.,1.0	,100};
  fHistoXiMCGen = new THnSparseF("fHistoXiMCGen","",3,bins_ximcgen,xmin_ximcgen,xmax_ximcgen);
  fOutputAll->Add(fHistoXiMCGen);

  fHistoLambdaPtvsDl=new TH2F("fHistoLambdaPtvsDl","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDl);
  fHistoLambdaPtvsDlSide=new TH2F("fHistoLambdaPtvsDlSide","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlSide);
  fHistoLambdaPtvsDlMCS=new TH2F("fHistoLambdaPtvsDlMCS","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlMCS);
  fHistoLambdaPtvsDR=new TH2F("fHistoLambdaPtvsDR","Lambda pt vs dl",20,0.,10.,80,0.,160.);
  fOutputAll->Add(fHistoLambdaPtvsDR);
  fHistoLambdaPtvsDRSide=new TH2F("fHistoLambdaPtvsDRSide","Lambda pt vs dl",20,0.,10.,80,0.,160.);
  fOutputAll->Add(fHistoLambdaPtvsDRSide);
  fHistoLambdaPtvsDRMCS=new TH2F("fHistoLambdaPtvsDRMCS","Lambda pt vs dl",20,0.,10.,80,0.,160.);
  fOutputAll->Add(fHistoLambdaPtvsDRMCS);

  fHistoEleXiPtvsRapidityRS=new TH2F("fHistoEleXiPtvsRapidityRS","EleXi pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleXiPtvsRapidityRS);
  fHistoEleXiPtvsRapidityWS=new TH2F("fHistoEleXiPtvsRapidityWS","EleXi pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleXiPtvsRapidityWS);
  fHistoEleXiPtvsRapidityMCS=new TH2F("fHistoEleXiPtvsRapidityMCS","EleXi pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleXiPtvsRapidityMCS);

  fHistoResponseElePt = new TH2D("fHistoResponseElePt","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoResponseElePt);
  fHistoResponseXiPt = new TH2D("fHistoResponseXiPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseXiPt);
  fHistoResponseEleXiPt = new TH2D("fHistoResponseEleXiPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleXiPt);
  fHistoResponseXiPtvsEleXiPt = new TH2D("fHistoResponseXiPtvsEleXiPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseXiPtvsEleXiPt);
  fHistoResponseXiPtXib = new TH2D("fHistoResponseXiPtXib","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseXiPtXib);
  fHistoResponseEleXiPtXib = new TH2D("fHistoResponseEleXiPtXib","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleXiPtXib);
  fHistoResponseMcGenXibPtvsXicPt = new TH2D("fHistoResponseMcGenXibPtvsXicPt","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseMcGenXibPtvsXicPt);


  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonElevsRunNumber=new TH1F("fHistonElevsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonElevsRunNumber);
  fHistonXivsRunNumber=new TH1F("fHistonXivsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonXivsRunNumber);
  fHistoMCEventType=new TH1F("fHistoMCEventType","",4,-0.5,3.5);
  fOutputAll->Add(fHistoMCEventType);
  fHistoMCXic0Decays=new TH1F("fHistoMCXic0Decays","",4,0.5,4.5);
  fOutputAll->Add(fHistoMCXic0Decays);
  fHistoMCDeltaPhiccbar=new TH1F("fHistoMCDeltaPhiccbar","",100,0.,3.2);
  fOutputAll->Add(fHistoMCDeltaPhiccbar);
  fHistoMCNccbar=new TH1F("fHistoMCNccbar","",100,-0.5,99.5);
  fOutputAll->Add(fHistoMCNccbar);
 
	fHistodPhiSdEtaSElectronProtonR125RS = new TH2D("fHistodPhiSdEtaSElectronProtonR125RS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronProtonR125RS);
	fHistodPhiSdEtaSElectronProtonR125WS = new TH2D("fHistodPhiSdEtaSElectronProtonR125WS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronProtonR125WS);
	fHistodPhiSdEtaSElectronProtonR125RSMix = new TH2D("fHistodPhiSdEtaSElectronProtonR125RSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronProtonR125RSMix);
	fHistodPhiSdEtaSElectronProtonR125WSMix = new TH2D("fHistodPhiSdEtaSElectronProtonR125WSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronProtonR125WSMix);
	fHistodPhiSdEtaSElectronPionR125RS = new TH2D("fHistodPhiSdEtaSElectronPionR125RS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronPionR125RS);
	fHistodPhiSdEtaSElectronPionR125WS = new TH2D("fHistodPhiSdEtaSElectronPionR125WS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronPionR125WS);
	fHistodPhiSdEtaSElectronPionR125RSMix = new TH2D("fHistodPhiSdEtaSElectronPionR125RSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronPionR125RSMix);
	fHistodPhiSdEtaSElectronPionR125WSMix = new TH2D("fHistodPhiSdEtaSElectronPionR125WSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronPionR125WSMix);
	fHistodPhiSdEtaSElectronBachelorR125RS = new TH2D("fHistodPhiSdEtaSElectronBachelorR125RS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronBachelorR125RS);
	fHistodPhiSdEtaSElectronBachelorR125WS = new TH2D("fHistodPhiSdEtaSElectronBachelorR125WS","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronBachelorR125WS);
	fHistodPhiSdEtaSElectronBachelorR125RSMix = new TH2D("fHistodPhiSdEtaSElectronBachelorR125RSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronBachelorR125RSMix);
	fHistodPhiSdEtaSElectronBachelorR125WSMix = new TH2D("fHistodPhiSdEtaSElectronBachelorR125WSMix","",50,0.,0.2,50,0.,0.2);
  fOutputAll->Add(fHistodPhiSdEtaSElectronBachelorR125WSMix);

	for(Int_t ih=0;ih<23;ih++){
		Int_t bins_eleptvscutvars[3];
		Double_t xmin_eleptvscutvars[3];
		Double_t xmax_eleptvscutvars[3];

		bins_eleptvscutvars[0] = 20;//electron pT bin
		xmin_eleptvscutvars[0] = 0.;
		xmax_eleptvscutvars[0] = 20.;
		bins_eleptvscutvars[2] = 10;//centrality bin
		xmin_eleptvscutvars[2] = 0.;
		xmax_eleptvscutvars[2] = 100.;

		if(ih==0 || ih==1){
			//0: TPC Ncluster 1: TPC ncluster PID
			bins_eleptvscutvars[1] = 40;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 160.;
		}else if(ih==2 || ih==3){
			//2: nSigma(TPC,e) 3: nSigma(TOF,e)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==4){
			//4: eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==5){
			//5: nITS cluster
			bins_eleptvscutvars[1] = 7;
			xmin_eleptvscutvars[1] = -0.5;
			xmax_eleptvscutvars[1] = 6.5;
		}else if(ih==6){
			//6: Lambda mass
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 1.1156-0.03;
			xmax_eleptvscutvars[1] = 1.1156+0.03;
		}else if(ih==7){
			//7: Xi mass
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 1.32-0.03;
			xmax_eleptvscutvars[1] = 1.32+0.03;
		}else if(ih==8 || ih==9){
			//8: Rfid Lambda, 9: Rfid Xi
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==10 || ih==11){
			//11: DCA Xi, 10: Dca V0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 2.;
		}else if(ih==12 || ih==13 || ih==14){
			//12: DCA Bachto prim, 13: DCA V0pr to prim 14: DCA V0pi to prim
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 0.5;
		}else if(ih==15 || ih==16){
			//16: CosPAXi, 15: CosPAv0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.95;
			xmax_eleptvscutvars[1] = 1.0;
		}else if(ih==17 || ih==18 || ih==19){
			//17: nSigma(TPC, bach)  18: nSigma(TPC, v0pr), 19: nSigma(TPC,v0pi)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==20 || ih==21){
			//20: V0 eta 21:  Xi eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==22){
			//20: Opening angle
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 3.141592/2;
		}

		fHistoElePtvsCutVarsRS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsRS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsRS[ih]);
		fHistoElePtvsCutVarsWS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsWS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsWS[ih]);
		fHistoElePtvsCutVarsMCS[ih] = new THnSparseF(Form("fHistoElePtvsCutVarsMCS[%d]",ih),"",3,bins_eleptvscutvars,xmin_eleptvscutvars,xmax_eleptvscutvars);
		fOutputAll->Add(fHistoElePtvsCutVarsMCS[ih]);
	}

  Int_t bins_eletag[3]=	{20 ,40	,10};
  Double_t xmin_eletag[3]={0.,0.	,0.0};
  Double_t xmax_eletag[3]={10.,20	,100};
  fHistoElectronPi0Total = new THnSparseF("fHistoElectronPi0Total","",3,bins_eletag,xmin_eletag,xmax_eletag);
  fOutputAll->Add(fHistoElectronPi0Total);
  fHistoElectronPi0Tag = new THnSparseF("fHistoElectronPi0Tag","",3,bins_eletag,xmin_eletag,xmax_eletag);
  fOutputAll->Add(fHistoElectronPi0Tag);
  fHistoElectronEtaTotal = new THnSparseF("fHistoElectronEtaTotal","",3,bins_eletag,xmin_eletag,xmax_eletag);
  fOutputAll->Add(fHistoElectronEtaTotal);
  fHistoElectronEtaTag = new THnSparseF("fHistoElectronEtaTag","",3,bins_eletag,xmin_eletag,xmax_eletag);
  fOutputAll->Add(fHistoElectronEtaTag);

  fHistoPi0MCGen = new TH1F("fHistoPi0MCGen","",100,0.,20.);
  fOutputAll->Add(fHistoPi0MCGen);
  fHistoEtaMCGen = new TH1F("fHistoEtaMCGen","",100,0.,20.);
  fOutputAll->Add(fHistoEtaMCGen);
  fHistoKaonMCGen = new TH1F("fHistoKaonMCGen","",100,0.,20.);
  fOutputAll->Add(fHistoKaonMCGen);
  fHistoD0MCGen = new TH1F("fHistoD0MCGen","",100,0.,20.);
  fOutputAll->Add(fHistoD0MCGen);

	//Axis 0: Pt
	//Axis 1: Dphi
	//Axis 2: proper dl
	//Axis 3: Sign Type
	//Axis 4: Conv Type
	//Axis 5: MC Type
	//Axis 6: Centrality
  Int_t bins_cor_nd[7]=	{100 , 20, 20, 4, 3, 10, 10};
  Double_t xmin_cor_nd[7]={0.,-M_PI,0.,-0.5,-0.5,-0.5,0.};
  Double_t xmax_cor_nd[7]={20.,M_PI,40.,3.5,2.5,9.5,100.};
  Double_t xmax_cor_nd2[7]={10.,M_PI,40.,3.5,2.5,9.5,100.};
  fHistoCorrelationVariablesvsEleXiPt = new THnSparseF("fHistoCorrelationVariablesvsEleXiPt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsEleXiPtMix = new THnSparseF("fHistoCorrelationVariablesvsEleXiPtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsEleXiPtMC = new THnSparseF("fHistoCorrelationVariablesvsEleXiPtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsElePt = new THnSparseF("fHistoCorrelationVariablesvsElePt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsElePtMix = new THnSparseF("fHistoCorrelationVariablesvsElePtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsElePtMC = new THnSparseF("fHistoCorrelationVariablesvsElePtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsXiPt = new THnSparseF("fHistoCorrelationVariablesvsXiPt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsXiPtMix = new THnSparseF("fHistoCorrelationVariablesvsXiPtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsXiPtMC = new THnSparseF("fHistoCorrelationVariablesvsXiPtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);

	//Axis 0: Mass
	//Axis 1: Pt
	//Axis 2: Near or Away
	//Axis 3: peak or Sideband
	//Axis 4: Sign Type
	//Axis 5: Conv Type
	//Axis 6: MC Type
	//Axis 7: Centrality
  Int_t bins_mass_nd[8]=	{22,100 , 2, 2, 4, 3, 10, 10};
  Double_t xmin_mass_nd[8]={1.3,0.,-0.5,-0.5,-0.5,-0.5,-0.5,0.};
  Double_t xmax_mass_nd[8]={5.7,20.,1.5,1.5,3.5,2.5,9.5,100.};
  Double_t xmax_mass_nd2[8]={5.7,10.,1.5,1.5,3.5,2.5,9.5,100.};
  fHistoMassVariablesvsEleXiPt = new THnSparseF("fHistoMassVariablesvsEleXiPt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsEleXiPtMix = new THnSparseF("fHistoMassVariablesvsEleXiPtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsEleXiPtMC = new THnSparseF("fHistoMassVariablesvsEleXiPtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsElePt = new THnSparseF("fHistoMassVariablesvsElePt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsElePtMix = new THnSparseF("fHistoMassVariablesvsElePtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsElePtMC = new THnSparseF("fHistoMassVariablesvsElePtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsXiPt = new THnSparseF("fHistoMassVariablesvsXiPt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsXiPtMix = new THnSparseF("fHistoMassVariablesvsXiPtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsXiPtMC = new THnSparseF("fHistoMassVariablesvsXiPtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);


  fOutputAll->Add(fHistoCorrelationVariablesvsEleXiPt);
  fOutputAll->Add(fHistoCorrelationVariablesvsEleXiPtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsEleXiPtMC);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePt);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePtMC);
  fOutputAll->Add(fHistoCorrelationVariablesvsXiPt);
  fOutputAll->Add(fHistoCorrelationVariablesvsXiPtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsXiPtMC);

  fOutputAll->Add(fHistoMassVariablesvsEleXiPt);
  fOutputAll->Add(fHistoMassVariablesvsEleXiPtMix);
  fOutputAll->Add(fHistoMassVariablesvsEleXiPtMC);
  fOutputAll->Add(fHistoMassVariablesvsElePt);
  fOutputAll->Add(fHistoMassVariablesvsElePtMix);
  fOutputAll->Add(fHistoMassVariablesvsElePtMC);
  fOutputAll->Add(fHistoMassVariablesvsXiPt);
  fOutputAll->Add(fHistoMassVariablesvsXiPtMix);
  fOutputAll->Add(fHistoMassVariablesvsXiPtMC);

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSEXic2eleXifromAODtracks::MakeCascadeHF(AliAODcascade *casc, AliAODTrack *part, AliAODTrack *partpid, AliAODEvent * aod, AliAODVertex *secVert) 
{
  //
  // Create AliAODRecoCascadeHF object from the argument
  //

  if(!casc) return 0x0;
  if(!part) return 0x0;
  if(!aod) return 0x0;

  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(casc,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;
  Double_t posprim[3]; primVertexAOD->GetXYZ(posprim);

  //------------------------------------------------
  // DCA between tracks
  //------------------------------------------------
  AliESDtrack *esdtrack = new AliESDtrack((AliVTrack*)partpid);

  Double_t xyz[3], pxpypz[3], cv[21]; Short_t sign;
  xyz[0]=casc->DecayVertexXiX();
  xyz[1]=casc->DecayVertexXiY();
  xyz[2]=casc->DecayVertexXiZ();
  pxpypz[0]=casc->MomXiX();
  pxpypz[1]=casc->MomXiY();
  pxpypz[2]=casc->MomXiZ();
  casc->GetCovarianceXYZPxPyPz(cv);
  sign=casc->ChargeXi();
  AliExternalTrackParam	*trackCasc = new AliExternalTrackParam(xyz,pxpypz,cv,sign);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackCasc,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackCasc->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }
  Double_t momcasc_new[3]={-9999,-9999,-9999.};
  trackCasc->GetPxPyPz(momcasc_new);

  Double_t px[2],py[2],pz[2];
  px[0] = part->Px(); py[0] = part->Py(); pz[0] = part->Pz(); 
  px[1] = momcasc_new[0]; py[1] = momcasc_new[1]; pz[1] = momcasc_new[2]; 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];

  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  d0[0]= d0z0bach[0];
  d0err[0] = TMath::Sqrt(covd0z0bach[0]);

  Double_t d0z0casc[2],covd0z0casc[3];
  trackCasc->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0casc,covd0z0casc);
  d0[1]= d0z0casc[0];
  d0err[1] = TMath::Sqrt(covd0z0casc[0]);

  //------------------------------------------------
  // Create AliAODRecoCascadeHF
  //------------------------------------------------
  Short_t charge = part->Charge();
  AliAODRecoCascadeHF *theCascade = new AliAODRecoCascadeHF(secVert,charge,px,py,pz,d0,d0err,dca);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack) delete esdtrack;
      if(trackCasc) delete trackCasc;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)part->GetID(),(UShort_t)trackCasc->GetID()};
  theCascade->SetProngIDs(2,id);

	theCascade->GetSecondaryVtx()->AddDaughter(part);
	theCascade->GetSecondaryVtx()->AddDaughter(casc);
  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackCasc) delete trackCasc;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent* aod)
{
  //
  // Make an array of tracks which should not be used in primary vertex calculation and 
  // Call PrimaryVertex function
  //

  TObjArray *TrackArray = new TObjArray(3);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk);
  TrackArray->AddAt(cptrk1,0);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(0));
  TrackArray->AddAt(cascptrack,1);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)casc->GetDaughter(1));
  TrackArray->AddAt(cascntrack,2);
  AliESDtrack *cascbtrack = new AliESDtrack((AliVTrack*)casc->GetDecayVertexXi()->GetDaughter(0));
  TrackArray->AddAt(cascbtrack,3);
  
  AliAODVertex *newvert  = PrimaryVertex(TrackArray,aod);
  
  for(Int_t i=0;i<4;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)TrackArray->UncheckedAt(i);
      delete tesd;
    }
  TrackArray->Clear();
  delete TrackArray;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::PrimaryVertex(const TObjArray *trkArray,
								   AliVEvent *event)
{
  //
  //Used only for pp
  //copied from AliAnalysisVertexingHF (except for the following 3 lines)
  //

  Bool_t fRecoPrimVtxSkippingTrks = kTRUE;
  Bool_t fRmTrksFromPrimVtx = kFALSE;

  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;
  
  //vertexESD = new AliESDVertex(*fV1);
  

  if(!fRecoPrimVtxSkippingTrks && !fRmTrksFromPrimVtx) { 
    // primary vertex from the input event
    
    vertexESD = new AliESDVertex(*fV1);
    
  } else {
    // primary vertex specific to this candidate
    
    Int_t nTrks = trkArray->GetEntriesFast();
    AliVertexerTracks *vertexer = new AliVertexerTracks(event->GetMagneticField());
    
    if(fRecoPrimVtxSkippingTrks) { 
      // recalculating the vertex
      
      if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraint")) {
	Float_t diamondcovxy[3];
	event->GetDiamondCovXY(diamondcovxy);
	Double_t pos[3]={event->GetDiamondX(),event->GetDiamondY(),0.};
	Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
	AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
	vertexer->SetVtxStart(diamond);
	delete diamond; diamond=NULL;
	if(strstr(fV1->GetTitle(),"VertexerTracksWithConstraintOnlyFitter")) 
	  vertexer->SetOnlyFitter();
      }
      Int_t skipped[1000];
      Int_t nTrksToSkip=0,id;
      AliExternalTrackParam *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliExternalTrackParam*)trkArray->UncheckedAt(i);
	id = (Int_t)t->GetID();
	if(id<0) continue;
	skipped[nTrksToSkip++] = id;
      }
      // TEMPORARY FIX
      // For AOD, skip also tracks without covariance matrix
      Double_t covtest[21];
      for(Int_t j=0; j<event->GetNumberOfTracks(); j++) {
	AliVTrack *vtrack = (AliVTrack*)event->GetTrack(j);
	if(!vtrack->GetCovarianceXYZPxPyPz(covtest)) {
	  id = (Int_t)vtrack->GetID();
	  if(id<0) continue;
	  skipped[nTrksToSkip++] = id;
	}
      }
      for(Int_t ijk=nTrksToSkip; ijk<1000; ijk++) skipped[ijk]=-1;
      //
      vertexer->SetSkipTracks(nTrksToSkip,skipped);
      vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(event); 
      
    } else if(fRmTrksFromPrimVtx && nTrks>0) { 
      // removing the prongs tracks
      
      TObjArray rmArray(nTrks);
      UShort_t *rmId = new UShort_t[nTrks];
      AliESDtrack *esdTrack = 0;
      AliESDtrack *t = 0;
      for(Int_t i=0; i<nTrks; i++) {
	t = (AliESDtrack*)trkArray->UncheckedAt(i);
	esdTrack = new AliESDtrack(*t);
	rmArray.AddLast(esdTrack);
	if(esdTrack->GetID()>=0) {
	  rmId[i]=(UShort_t)esdTrack->GetID();
	} else {
	  rmId[i]=9999;
	}
      }
      Float_t diamondxy[2]={static_cast<Float_t>(event->GetDiamondX()),static_cast<Float_t>(event->GetDiamondY())};
      vertexESD = vertexer->RemoveTracksFromVertex(fV1,&rmArray,rmId,diamondxy);
      delete [] rmId; rmId=NULL;
      rmArray.Delete();
      
    }
    
    delete vertexer; vertexer=NULL;
    if(!vertexESD) return vertexAOD;
    if(vertexESD->GetNContributors()<=0) { 
      //AliDebug(2,"vertexing failed"); 
      delete vertexESD; vertexESD=NULL;
      return vertexAOD;
    }
    
    
  }
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  
  return vertexAOD;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSEXic2eleXifromAODtracks::ReconstructSecondaryVertex(AliAODcascade *casc, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
  //
	
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(casc,part,aod);
    if(!primVertexAOD){
      primVertexAOD = fVtx1;
    }else{
      unsetvtx = kTRUE;
    }
  }else{
    primVertexAOD = fVtx1;
  }
  if(!primVertexAOD) return 0x0;

  AliESDVertex * vertexESD = new AliESDVertex(*fV1);

  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  AliAODVertex *secVert = new AliAODVertex(pos,cov,chi2perNDF);

  return secVert;
}
//________________________________________________________________________
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::MatchToMC(AliAODRecoCascadeHF *exobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_casc, Int_t *labelarray_ele, Int_t *labelarray_casc,  Int_t &ngen_ele, Int_t &ngen_casc) 
{
  //
  // Match to MC
  //
	for(Int_t i=0;i<100;i++){
		pdgarray_ele[i] = -9999;
		labelarray_ele[i] = -9999;
		pdgarray_casc[i] = -9999;
		labelarray_casc[i] = -9999;
	}
	ngen_ele = 0;
	ngen_casc = 0;

  AliVTrack *trk = dynamic_cast<AliVTrack*>(exobj->GetBachelor());
  if(!trk) return -1;
  Int_t labEle = trk->GetLabel();
	if(labEle<0) return -1;
	AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labEle);
	if(!mcetrk) return -1;
	labelarray_ele[0] = labEle;
	pdgarray_ele[0] = mcetrk->GetPdgCode();
	ngen_ele ++;

  AliAODMCParticle *mcprimele=0;
  mcprimele = mcetrk;
  while(mcprimele->GetMother()>=0) {
    Int_t labprim_ele=mcprimele->GetMother();
    AliAODMCParticle *tmcprimele = (AliAODMCParticle*)mcArray->At(labprim_ele);
    if(!tmcprimele) {
			break;
    }

    mcprimele = tmcprimele;
		pdgarray_ele[ngen_ele] = mcprimele->GetPdgCode();
		labelarray_ele[ngen_ele] = labprim_ele;
		ngen_ele ++;
		if(ngen_ele==100) break;
  }

  AliAODcascade *theCascade = dynamic_cast<AliAODcascade*>(exobj->GetCascade());
	if(!theCascade) return -1;

	Int_t pdgDgcasc[2]={211,3122};
	Int_t pdgDgv0[2]={2212,211};
  Int_t labcasc = MatchToMCCascade(theCascade,3312,pdgDgcasc,pdgDgv0,mcArray); // the cascade
  if(labcasc<0) return -1;

	AliAODMCParticle *mccasc = (AliAODMCParticle*)mcArray->At(labcasc);
	if(!mccasc) return -1;
	labelarray_casc[0] = labcasc;
	pdgarray_casc[0] = mccasc->GetPdgCode();
	ngen_casc ++;

  AliAODMCParticle *mcprimcasc=0;
  mcprimcasc = mccasc;
  while(mcprimcasc->GetMother()>=0) {
    Int_t labprim_casc=mcprimcasc->GetMother();
    AliAODMCParticle *tmcprimcasc = (AliAODMCParticle*)mcArray->At(labprim_casc);
    if(!tmcprimcasc) {
			break;
    }

    mcprimcasc = tmcprimcasc;
		pdgarray_casc[ngen_casc] = mcprimcasc->GetPdgCode();
		labelarray_casc[ngen_casc] = labprim_casc;
		ngen_casc ++;
		if(ngen_casc==100) break;
  }

	Bool_t same_flag = kFALSE;
	Int_t matchedlabel=-9999;
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_casc;ivmc++){
			if(labelarray_ele[iemc]==labelarray_casc[ivmc]){
				same_flag = kTRUE;
				matchedlabel = labelarray_ele[iemc];
				break;
			}
		}
		if(same_flag) break;
	}

	return matchedlabel;

}
//________________________________________________________________________
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::MatchToMCCascade(AliAODcascade *theCascade, Int_t pdgabscasc, Int_t *pdgDgcasc, Int_t *pdgDgv0, TClonesArray *mcArray) const // the cascade
{
	//
	// Matching to MC of cascade
	//

	AliAODTrack *cptrack = (AliAODTrack*) theCascade->GetDaughter(0);
	if(!cptrack) return -1;
	Int_t label_p = TMath::Abs(cptrack->GetLabel());
	if(label_p<0) return -1;
	AliAODTrack *cntrack = (AliAODTrack*) theCascade->GetDaughter(1);
	if(!cntrack) return -1;
	Int_t label_n = TMath::Abs(cntrack->GetLabel());
	if(label_n<0) return -1;
	Int_t labv0 = theCascade->MatchToMC(pdgDgcasc[1],mcArray,2,pdgDgv0);
	if(labv0<0) return -1;
	AliAODMCParticle *mcpartv0= (AliAODMCParticle*) mcArray->At(labv0);

	AliAODTrack *cbtrack = (AliAODTrack*) theCascade->GetDecayVertexXi()->GetDaughter(0);
	if(!cbtrack) return -1;

	Int_t label_b = TMath::Abs(cbtrack->GetLabel());
	if(label_b<0) return -1;

	AliAODMCParticle *mcpartb= (AliAODMCParticle*) mcArray->At(label_b);
	Int_t pdgb = TMath::Abs(mcpartb->GetPdgCode());
	if(pdgb!=pdgDgcasc[0]) return -1;

	AliAODMCParticle *mcmotherv0=mcpartv0;
	Bool_t isFromXiv0 = kFALSE;
	Int_t labxiv0 = mcmotherv0->GetMother();
	if(labxiv0<0) return -1;
	mcmotherv0 =  (AliAODMCParticle*) mcArray->At(labxiv0);
	if(mcmotherv0){
		Int_t pdg = TMath::Abs(mcmotherv0 ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXiv0 = kTRUE;
		}
	}
	if(!isFromXiv0) return -1;

	AliAODMCParticle *mcmotherb=mcpartb;
	Bool_t isFromXib = kFALSE;
	Int_t labxib = mcmotherb->GetMother();
	if(labxib<0) return -1;
	mcmotherb =  (AliAODMCParticle*) mcArray->At(labxib);
	if(mcmotherb){
		Int_t pdg = TMath::Abs(mcmotherb ->GetPdgCode());
		if(pdg==pdgabscasc){
			isFromXib = kTRUE;
		}
	}
	if(!isFromXib) return -1;

	if(labxiv0!=labxib) return -1;//Bachelor and V0 should come from the same Xi

	return labxib;
}
//________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
  //const Int_t entries = event->GetNumberOfTracks();
  if(trkEntries==0) return;
  
  nSeleTrks=0;
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
    //if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;

    if(!fAnalCuts) continue;
    
    AliAODTrack *aodt = (AliAODTrack*)track;

    if(fAnalCuts->GetProdUseAODFilterBit()){
      Int_t filterbit = fAnalCuts->GetProdAODFilterBit();
      if(filterbit==7){
        if(!aodt->TestFilterBit(BIT(filterbit))) continue;
      }else{
        if(!aodt->TestFilterMask(BIT(filterbit))) continue;
      }
    }

    AliAODTrack *aodtpid = 0;
    if(fAnalCuts->GetProdAODFilterBit()==7){
      aodtpid = fGTI[-aodt->GetID()-1];
    }else{
      aodtpid = aodt;
    }

		Double_t nsigma_tpcele = -9999;
		Double_t nsigma_tofele = -9999;
		if(fAnalCuts->GetIsUsePID()){
			nsigma_tpcele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(aodtpid,AliPID::kElectron);
			nsigma_tofele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(aodtpid,AliPID::kElectron);
		}

    if(fAnalCuts->SingleTrkCutsNoPID(aodt,aodtpid,fVtx1)){
			fHistoElectronTPCPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFPID->Fill(aodt->Pt(),nsigma_tofele);
			if(fabs(nsigma_tofele)<3.){
				fHistoElectronTPCPIDSelTOF->Fill(aodt->Pt(),nsigma_tpcele);
				Double_t eleeta = aodt->Eta();
				if(fabs(eleeta)<0.6)
					fHistoElectronTPCPIDSelTOFSmallEta->Fill(aodt->Pt(),nsigma_tpcele);
				if(fabs(eleeta)>0.6 && fabs(eleeta)<0.8)
					fHistoElectronTPCPIDSelTOFLargeEta->Fill(aodt->Pt(),nsigma_tpcele);
				if(eleeta>-0.8 && eleeta<-0.6){
					fHistoElectronTPCPIDSelTOFEtaDep[0]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.6&&eleeta<-0.4){
					fHistoElectronTPCPIDSelTOFEtaDep[1]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.4&&eleeta<-0.2){
					fHistoElectronTPCPIDSelTOFEtaDep[2]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>-0.2&&eleeta<0.0){
					fHistoElectronTPCPIDSelTOFEtaDep[3]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.0&&eleeta<0.2){
					fHistoElectronTPCPIDSelTOFEtaDep[4]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.2&&eleeta<0.4){
					fHistoElectronTPCPIDSelTOFEtaDep[5]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.4&&eleeta<0.6){
					fHistoElectronTPCPIDSelTOFEtaDep[6]->Fill(aodt->Pt(),nsigma_tpcele);
				}else if(eleeta>0.6&&eleeta<0.8){
					fHistoElectronTPCPIDSelTOFEtaDep[7]->Fill(aodt->Pt(),nsigma_tpcele);
				}
			}
			if(nsigma_tpcele>-0.5&&nsigma_tpcele<3.){
				fHistoElectronTOFPIDSelTPC->Fill(aodt->Pt(),nsigma_tofele);
			}
		}
    if(fAnalCuts->SingleTrkCuts(aodt,aodtpid,fVtx1)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
			fHistoElectronTPCSelPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFSelPID->Fill(aodt->Pt(),nsigma_tofele);
    }
  } // end loop on tracks
}
//________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::SelectCascade( const AliVEvent *event,Int_t nCascs,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray)
{
  //
  // Select good Casc using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

	Double_t primVtx[3];
	fVtx1->GetXYZ(primVtx);

  nSeleCasc = 0;
  for(Int_t icasc=0;icasc<nCascs;icasc++)
    {
      seleCascFlags[icasc] = kFALSE;
      AliAODcascade *casc = ((AliAODEvent*)event)->GetCascade(icasc);

      if(!fAnalCuts) continue;
      if(fAnalCuts->SingleCascadeCuts(casc,primVtx)){
				seleCascFlags[icasc] = kTRUE;
				nSeleCasc++;
      }
    }
}
//_________________________________________________________________
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult, Double_t rp){
	//
  // check in which of the pools the current event falls
	//

  Int_t theBinZ=TMath::BinarySearch(fNzVtxBins,fZvtxBins,zvert);
  if(theBinZ<0 || theBinZ>=fNzVtxBins) return -1;
  Int_t theBinM=TMath::BinarySearch(fNCentBins,fCentBins,mult);
  if(theBinM<0 || theBinM>=fNCentBins) return -1;
  //return fNCentBins*theBinZ+theBinM;
  Int_t theBinR=TMath::BinarySearch(fNRPBins,fRPBins,rp);
  if(theBinR<0 || theBinR>=fNRPBins) return -1;
  return fNzVtxBins*fNCentBins*theBinR+fNCentBins*theBinZ+theBinM;
}
//_________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::DoEventMixingWithPools(Int_t poolIndex)
{
  //
  // perform mixed event analysis
  //
  Int_t nextRes( nextResVec[poolIndex] );
  Int_t KiddiePool = m_ReservoirE[poolIndex].size();
  if( !reservoirsReady[poolIndex] )  KiddiePool = nextRes;

  if( KiddiePool>0 )
  {
    for(Int_t j=0;j<KiddiePool;j++){
      if( j!=nextRes ){
        FillBackground(m_ReservoirE[poolIndex][nextRes],m_ReservoirVarsE[poolIndex][nextRes],m_ReservoirL1[poolIndex][j],m_ReservoirVarsL1[poolIndex][j],1);
        FillBackground(m_ReservoirE[poolIndex][j],m_ReservoirVarsE[poolIndex][j],m_ReservoirL1[poolIndex][nextRes],m_ReservoirVarsL1[poolIndex][nextRes],1);
        FillBackground(m_ReservoirE[poolIndex][nextRes],m_ReservoirVarsE[poolIndex][nextRes],m_ReservoirL2[poolIndex][j],m_ReservoirVarsL2[poolIndex][j],-1);
        FillBackground(m_ReservoirE[poolIndex][j],m_ReservoirVarsE[poolIndex][j],m_ReservoirL2[poolIndex][nextRes],m_ReservoirVarsL2[poolIndex][nextRes],-1);
      }
    }
  }
}
//_________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillBackground(std::vector<TLorentzVector * > mixTypeE,std::vector<TVector * > mixTypeEVars, std::vector<TLorentzVector * > mixTypeL, std::vector<TVector * > mixTypeLVars, Int_t charge_xi)
{     
  //
  // Fill background
  //
  int nEle = mixTypeE.size();
  int nCasc = mixTypeL.size();
  for(Int_t ie=0;ie<nEle;ie++){
    TLorentzVector* trke=mixTypeE[ie];
    if(!trke) continue;
    TVector *elevars = mixTypeEVars[ie];
    for(Int_t iv=0;iv<nCasc;iv++){
      TLorentzVector* casc=mixTypeL[iv];
      TVector *cascvars = mixTypeLVars[iv];
      if(!casc) continue;
      FillMixROOTObjects(trke,casc,elevars,cascvars,charge_xi);
    }
  }
  return;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskSEXic2eleXifromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
	//
  // Analyze AliAODmcparticle
	//

	Int_t nmcpart = mcArray->GetEntriesFast();

	Int_t mcevttype = 0;
	Int_t nccbar = 0;
	Bool_t sigmaevent = kFALSE;
	if(fMCEventType==1 || fMCEventType==2 || fMCEventType==11 || fMCEventType==12){
		//1: c quark event
		//2: b quark event
		//11: near side c-cbar event
		//12: away side c-cbar event
		Int_t ncquark = 0;
		Int_t ncbarquark = 0;
		Double_t phi_c = -9999.;
		Double_t phi_cbar = -9999.;
		for(Int_t i=0;i<nmcpart;i++)
		{
			AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
			if(TMath::Abs(mcpart->GetPdgCode())==4){
				if(fabs(mcpart->Y())<1.5){
					if(mcpart->GetPdgCode()==4){
						phi_c = mcpart->Phi();
						ncquark++;
					}
					if(mcpart->GetPdgCode()==-4){
						phi_cbar = mcpart->Phi();
						ncbarquark++;
					}
					if(mcevttype==0){
						mcevttype = 1;
					}else if(mcevttype==1){
						mcevttype = 1;
					}else if(mcevttype==2){
						mcevttype = 3;
					}else if(mcevttype==3){
						mcevttype = 3;
					}
					nccbar++;
				}
			}
			if(TMath::Abs(mcpart->GetPdgCode())==5){
				if(!mcpart->IsPhysicalPrimary()) continue;
				if(fabs(mcpart->Y())<1.5){
					if(mcevttype==0){
						mcevttype = 2;
					}else if(mcevttype==1){
						mcevttype = 3;
					}else if(mcevttype==2){
						mcevttype = 2;
					}else if(mcevttype==3){
						mcevttype = 3;
					}
				}
			}
		}

		if(fMCEventType==1||fMCEventType==11||fMCEventType==12){
			if((mcevttype==2)||(mcevttype==0)||(mcevttype==3)) return kFALSE;
		}else if(fMCEventType==2){
			if((mcevttype==1)||(mcevttype==0)||(mcevttype==3)) return kFALSE;
		}

		if(fMCEventType>10){
			if(ncquark!=1) return kFALSE;
			if(ncbarquark!=1) return kFALSE;
			Double_t dphi = fabs(phi_c - phi_cbar);
			if(dphi>2*M_PI) dphi -= 2*M_PI;
			if(dphi>M_PI) dphi = 2*M_PI-dphi;
			if(fMCEventType==11 && dphi>M_PI/3.) return kFALSE;
			if(fMCEventType==12 && dphi<2*M_PI/3.) return kFALSE;
			fHistoMCDeltaPhiccbar->Fill(dphi);
		}

		fHistoMCEventType->Fill(mcevttype);
		fHistoMCNccbar->Fill(nccbar);
	}

	for(Int_t i=0;i<nmcpart;i++)
	{
		AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
		if(TMath::Abs(mcpart->GetPdgCode())==4132){
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3312){
					xi_flag = kTRUE;
					mccascpart = mcdau;
				}
			}

			Int_t decaytype = -9999;
			if(e_flag && xi_flag) decaytype = 0;

			if(e_flag&&xi_flag)
				fHistoMCXic0Decays->Fill(1);
			if(!e_flag&&xi_flag)
				fHistoMCXic0Decays->Fill(2);
			if(e_flag&&!xi_flag)
				fHistoMCXic0Decays->Fill(3);
			if(!e_flag&&!xi_flag)
				fHistoMCXic0Decays->Fill(4);

			FillMCROOTObjects(mcpart,mcepart,mccascpart,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==5132 || TMath::Abs(mcpart->GetPdgCode())==5232){
			Bool_t e_flag = kFALSE;
			Bool_t xic_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mcxicpart = 0;
			AliAODMCParticle *mccascpart = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==4132 || TMath::Abs(mcdau->GetPdgCode())==4232 ){
					xic_flag = kTRUE;
					mcxicpart = mcdau;
				}
			}

			Bool_t xi_flag = kFALSE;
      if(e_flag && xic_flag){
        for(Int_t idau=mcxicpart->GetFirstDaughter();idau<mcxicpart->GetLastDaughter()+1;idau++)
        {
          if(idau<0) break;
          AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
          if(!mcdau) continue;
          if(TMath::Abs(mcdau->GetPdgCode())==3312){
            xi_flag = kTRUE;
            mccascpart = mcdau;
          }
        }
      }

			if(xic_flag){
				Double_t contmc_withxic[3];
				contmc_withxic[0] = mcpart->Pt();
				contmc_withxic[1] = mcpart->Y();
				contmc_withxic[2] = mcxicpart->Y();
				if(fabs(mcxicpart->Y())<1.){
					fHistoResponseMcGenXibPtvsXicPt->Fill(mcpart->Pt(),mcxicpart->Pt());
				}
				fHistoXibMCGenWithXic->Fill(contmc_withxic);
			}

			Int_t decaytype = -9999;
			if(e_flag && xic_flag && xi_flag) decaytype = 10;
			FillMCROOTObjects(mcpart,mcepart,mccascpart,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==11 && mcpart->GetStatus()==1){
			AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
			Float_t etamin, etamax;
			esdcuts->GetEtaRange(etamin,etamax);
			if(fabs(mcpart->Eta())<etamax){
				Bool_t gamma_flag = kFALSE;
				Int_t labmother = mcpart->GetMother();
				if(labmother>=0){
					AliAODMCParticle *mcmother = (AliAODMCParticle*) mcArray->At(labmother);
					Int_t pdgmother = mcmother->GetPdgCode();
					if(TMath::Abs(pdgmother)==22) gamma_flag = kTRUE;
				}
				if(!gamma_flag) fHistoBachPtMCGen->Fill(mcpart->Pt());
			}
			FillMCEleROOTObjects(mcpart, mcArray);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==3312){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdCascEtaRange(etamin,etamax);
			fAnalCuts->GetProdCascRapRange(rapmin,rapmax);

			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax)){
				fHistoXiMassvsPtMCGen->Fill(1.32171, mcpart->Pt());
			}
			FillMCCascROOTObjects(mcpart, mcArray);
		}

		if(TMath::Abs(mcpart->GetPdgCode())==111){
      if(fabs(mcpart->Y())<0.5){
        fHistoPi0MCGen->Fill(mcpart->Pt());
      }
    }
		if(TMath::Abs(mcpart->GetPdgCode())==221){
      if(fabs(mcpart->Y())<0.5){
        fHistoEtaMCGen->Fill(mcpart->Pt());
      }
    }
		if(TMath::Abs(mcpart->GetPdgCode())==321){
      if(fabs(mcpart->Y())<0.5){
        fHistoKaonMCGen->Fill(mcpart->Pt());
      }
    }
		if(TMath::Abs(mcpart->GetPdgCode())==421){
      if(fabs(mcpart->Y())<0.5){
        fHistoD0MCGen->Fill(mcpart->Pt());
      }
    }
	}

	if(fMCDoPairAnalysis)
	{
		for(Int_t i=0;i<nmcpart;i++)
		{
			AliAODMCParticle *mcparte = (AliAODMCParticle*) mcArray->At(i);
			if(!mcparte) continue;
			if(TMath::Abs(mcparte->GetPdgCode())!=11) continue;
			if(mcparte->GetStatus()!=1) continue;
			if(mcparte->Pt()<0.4) continue;//Apply rough cuts
			if(fabs(mcparte->Eta())>0.8) continue;//Apply rough cuts
			for(Int_t j=0;j<nmcpart;j++)
			{
				AliAODMCParticle *mcpartv = (AliAODMCParticle*) mcArray->At(j);
				if(!mcpartv) continue;
				if(TMath::Abs(mcpartv->GetPdgCode())!=3312) continue;
				if(mcpartv->Pt()<0.4) continue;//Apply rough cuts
				if(fabs(mcpartv->Eta())>0.8) continue;//Apply rough cuts
				if(mcpartv->GetNDaughters()!=2) continue;//Apply rough cuts

				FillMCGenPairROOTObjects(mcparte,mcpartv,mcArray);
			}
		}
		return kFALSE;
	}

	return kTRUE;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineMCGenPairTreeVariables() 
{
  //
  // Define mc pair tree variables
  //

  const char* nameoutput = GetOutputSlot(11)->GetContainer()->GetName();
  fMCGenPairVariablesTree = new TTree(nameoutput,"MC pair variables tree");
  Int_t nVar = 38;
  fCandidateMCGenPairVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

	fCandidateVariableNames[ 0] = "InvMassEleXi";
	fCandidateVariableNames[ 1] = "EleXiPx";
	fCandidateVariableNames[ 2] = "EleXiPy";
	fCandidateVariableNames[ 3] = "EleXiPz";
	fCandidateVariableNames[ 4] = "ElePdgCode";
	fCandidateVariableNames[ 5] = "ElePx";
	fCandidateVariableNames[ 6] = "ElePy";
	fCandidateVariableNames[ 7] = "ElePz";
	fCandidateVariableNames[ 8] = "XiPdgCode";
	fCandidateVariableNames[ 9] = "XiPx";
	fCandidateVariableNames[10] = "XiPy";
	fCandidateVariableNames[11] = "XiPz";
	fCandidateVariableNames[12] = "SameFlag";
	fCandidateVariableNames[13] = "EleNGeneration";
	fCandidateVariableNames[14] = "EleGen1PDG";
	fCandidateVariableNames[15] = "EleGen2PDG";
	fCandidateVariableNames[16] = "EleGen3PDG";
	fCandidateVariableNames[17] = "EleGen4PDG";
	fCandidateVariableNames[18] = "EleGen5PDG";
	fCandidateVariableNames[19] = "EleGen6PDG";
	fCandidateVariableNames[20] = "EleGen7PDG";
	fCandidateVariableNames[21] = "EleGen8PDG";
	fCandidateVariableNames[22] = "EleGen9PDG";
	fCandidateVariableNames[23] = "EleGen10PDG";
	fCandidateVariableNames[24] = "ElePrimPDG";
	fCandidateVariableNames[25] = "XiNGeneration";
	fCandidateVariableNames[26] = "XiGen1PDG";
	fCandidateVariableNames[27] = "XiGen2PDG";
	fCandidateVariableNames[28] = "XiGen3PDG";
	fCandidateVariableNames[29] = "XiGen4PDG";
	fCandidateVariableNames[30] = "XiGen5PDG";
	fCandidateVariableNames[31] = "XiGen6PDG";
	fCandidateVariableNames[32] = "XiGen7PDG";
	fCandidateVariableNames[33] = "XiGen8PDG";
	fCandidateVariableNames[34] = "XiGen9PDG";
	fCandidateVariableNames[35] = "XiGen10PDG";
	fCandidateVariableNames[36] = "XiPrimPDG";
	fCandidateVariableNames[37] = "MatchedPDG";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCGenPairVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCGenPairVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::FillMCGenPairROOTObjects(AliAODMCParticle *mcparte, AliAODMCParticle *mcpartv, TClonesArray *mcArray) 
{
  //
  // Fill histograms or mc pair analysis tree 
  //
	for(Int_t i=0;i<38;i++){
		fCandidateMCGenPairVariables[i] = -9999.;
	}

	TLorentzVector vele, vxi, vexi;
	vele.SetXYZM(mcparte->Px(),mcparte->Py(),mcparte->Pz(),0.000510998928);
	vxi.SetXYZM(mcpartv->Px(),mcpartv->Py(),mcpartv->Pz(),1.32171);
	vexi = vele + vxi;

	Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;
	Int_t pdgarray_xi[100], labelarray_xi[100], ngen_xi;
	GetMCDecayHistory(mcparte,mcArray,pdgarray_ele,labelarray_ele,ngen_ele);
	GetMCDecayHistory(mcpartv,mcArray,pdgarray_xi,labelarray_xi,ngen_xi);

	Bool_t same_flag = kFALSE;
	Int_t matched_pdg = -999999;
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_xi;ivmc++){
			if(labelarray_ele[iemc]==labelarray_xi[ivmc]){
				same_flag = kTRUE;
				matched_pdg = pdgarray_ele[iemc];
				break;
			}
		}
		if(same_flag) break;
	}
	Int_t pdgprim_ele = pdgarray_ele[ngen_ele-1];
	Int_t pdgprim_xi = pdgarray_xi[ngen_xi-1];

	fCandidateMCGenPairVariables[ 0] = vexi.M();
	fCandidateMCGenPairVariables[ 1] = vexi.Px();
	fCandidateMCGenPairVariables[ 2] = vexi.Py();
	fCandidateMCGenPairVariables[ 3] = vexi.Pz();
	fCandidateMCGenPairVariables[ 4] = mcparte->GetPdgCode();
	fCandidateMCGenPairVariables[ 5] = vele.Px();
	fCandidateMCGenPairVariables[ 6] = vele.Py();
	fCandidateMCGenPairVariables[ 7] = vele.Pz();
	fCandidateMCGenPairVariables[ 8] = mcpartv->GetPdgCode();
	fCandidateMCGenPairVariables[ 9] = vxi.Px();
	fCandidateMCGenPairVariables[10] = vxi.Py();
	fCandidateMCGenPairVariables[11] = vxi.Pz();
	fCandidateMCGenPairVariables[12] = (Float_t)same_flag;
	fCandidateMCGenPairVariables[13] = (Float_t)ngen_ele;
	fCandidateMCGenPairVariables[14] = (Float_t)pdgarray_ele[0];
	fCandidateMCGenPairVariables[15] = (Float_t)pdgarray_ele[1];
	fCandidateMCGenPairVariables[16] = (Float_t)pdgarray_ele[2];
	fCandidateMCGenPairVariables[17] = (Float_t)pdgarray_ele[3];
	fCandidateMCGenPairVariables[18] = (Float_t)pdgarray_ele[4];
	fCandidateMCGenPairVariables[19] = (Float_t)pdgarray_ele[5];
	fCandidateMCGenPairVariables[20] = (Float_t)pdgarray_ele[6];
	fCandidateMCGenPairVariables[21] = (Float_t)pdgarray_ele[7];
	fCandidateMCGenPairVariables[22] = (Float_t)pdgarray_ele[8];
	fCandidateMCGenPairVariables[23] = (Float_t)pdgarray_ele[9];
	fCandidateMCGenPairVariables[24] = (Float_t)pdgarray_ele[ngen_ele-1];
	fCandidateMCGenPairVariables[25] = (Float_t)ngen_xi;
	fCandidateMCGenPairVariables[26] = (Float_t)pdgarray_xi[0];
	fCandidateMCGenPairVariables[27] = (Float_t)pdgarray_xi[1];
	fCandidateMCGenPairVariables[28] = (Float_t)pdgarray_xi[2];
	fCandidateMCGenPairVariables[29] = (Float_t)pdgarray_xi[3];
	fCandidateMCGenPairVariables[30] = (Float_t)pdgarray_xi[4];
	fCandidateMCGenPairVariables[31] = (Float_t)pdgarray_xi[5];
	fCandidateMCGenPairVariables[32] = (Float_t)pdgarray_xi[6];
	fCandidateMCGenPairVariables[33] = (Float_t)pdgarray_xi[7];
	fCandidateMCGenPairVariables[34] = (Float_t)pdgarray_xi[8];
	fCandidateMCGenPairVariables[35] = (Float_t)pdgarray_xi[9];
	fCandidateMCGenPairVariables[36] = (Float_t)pdgarray_xi[ngen_xi-1];
	fCandidateMCGenPairVariables[37] = (Float_t) matched_pdg;

	fMCGenPairVariablesTree->Fill();
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSEXic2eleXifromAODtracks::DefineCorrelationTreeVariables() 
{
  //
  // Define mc pair tree variables
  //

  const char* nameoutput = GetOutputSlot(12)->GetContainer()->GetName();
  fCorrelationVariablesTree = new TTree(nameoutput,"Correlation variables tree");
  Int_t nVar = 15;
  fCorrelationVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[0] = "XiPt";
  fCandidateVariableNames[1] = "ElePt";
  fCandidateVariableNames[2] = "DeltaPhi";
  fCandidateVariableNames[3] = "DeltaEta";
  fCandidateVariableNames[4] = "V0ProperDecayLength";
  fCandidateVariableNames[5] = "Eled0";
  fCandidateVariableNames[6] = "FGMixMC";
  fCandidateVariableNames[7] = "SignType";
  fCandidateVariableNames[8] = "Convtype";
  fCandidateVariableNames[9] = "MCType";
  fCandidateVariableNames[10] = "Centrality";
  fCandidateVariableNames[11] = "EleXiPt";
  fCandidateVariableNames[12] = "EleXiMass";
  fCandidateVariableNames[13] = "EleXiCosOA";
  fCandidateVariableNames[14] = "MCEleMother";


  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fCorrelationVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCorrelationVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}

////-------------------------------------------------------------------------------
void	AliAnalysisTaskSEXic2eleXifromAODtracks::GetMCDecayHistory(AliAODMCParticle *mcpart, TClonesArray *mcArray, Int_t *pdgarray, Int_t *labelarray, Int_t &ngen)
{
//
// MC decay history
//

	for(Int_t i=0;i<100;i++){
		pdgarray[i] = -9999;
		labelarray[i] = -9999;
	}
	ngen = 0;

	AliAODMCParticle *mcprim = mcpart;
  while(mcprim->GetMother()>=0) {
    Int_t lab_prim=mcprim->GetMother();

    AliAODMCParticle *tmcprim = (AliAODMCParticle*)mcArray->At(lab_prim);
    if(!tmcprim) {
			break;
    }
		if((TMath::Abs(tmcprim->GetPdgCode())<10) || (TMath::Abs(tmcprim->GetPdgCode())==21)) break;

		mcprim = tmcprim;

		pdgarray[ngen] = mcprim->GetPdgCode();
		labelarray[ngen] = lab_prim;

		ngen ++;
		if(ngen == 100) break;
	}
}
//________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::StoreGlobalTrackReference(AliAODTrack *track, Int_t index){
  //
  // Stores the pointer to the global track
  // copied from femtoscopy/k0analysis/plamanalysis
  //

  // Check that the id is positive
  if(track->GetID()<0){
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    return;
  }

  // Check id is not too big for buffer
  if(track->GetID()>=fTrackBuffSize){
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
        ,track->GetID(),fTrackBuffSize);
    return;
  }

  // Warn if we overwrite a track
  if(fGTI[track->GetID()]){
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if( (!track->GetFilterMap()) &&
        (!track->GetTPCNcls())   )
      return;

    // Imagine the other way around,
    // the zero map zero clusters track
    // is stored and the good one wants 
    // to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if( fGTI[track->GetID()]->GetFilterMap() ||
        fGTI[track->GetID()]->GetTPCNcls()   ){
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
          (fGTI[track->GetID()])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
          (fGTI[track->GetID()])->GetFilterMap(),track->GetFilterMap());
    }
  } // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  // 	   ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[track->GetID()]) = track;
  (fGTIndex[track->GetID()]) = index;
}
//________________________________________________________________________
void AliAnalysisTaskSEXic2eleXifromAODtracks::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
    fGTIndex[i]=-9999;
  }
}
//________________________________________________________________________
Int_t AliAnalysisTaskSEXic2eleXifromAODtracks::FromSemileptonicDecays(Int_t *history){
  // 
  // Check the mother 
  // 
  if(abs(history[1])==411) return 1;
  if(abs(history[1])==421) return 1;
  if(abs(history[1])==431) return 1;
  if(abs(history[1])==4122) return 1;
  if(abs(history[1])==4132) return 1;
  if(abs(history[1])==4232) return 1;
  if(abs(history[1])==4332) return 1;

  if(abs(history[1])==511) return 2;
  if(abs(history[1])==521) return 2;
  if(abs(history[1])==531) return 2;
  if(abs(history[1])==5122) return 2;
  if(abs(history[1])==5132) return 2;
  if(abs(history[1])==5232) return 2;
  if(abs(history[1])==5332) return 2;

  if(abs(history[1])==130) return 3;
  if(abs(history[1])==310) return 3;
  if(abs(history[1])==311) return 3;
  if(abs(history[1])==321) return 3;
  if(abs(history[1])==3122) return 3;
  if(abs(history[1])==3312) return 3;
  if(abs(history[1])==3322) return 3;
  if(abs(history[1])==3334) return 3;

  return 0;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSEXic2eleXifromAODtracks::HaveCharmInHistory(Int_t *history){
  // 
  // Check if the particle is from charm hadrons
  // 
  for(Int_t ih=0;ih<10;ih++){
    if(abs(history[ih])==411) return kTRUE;
    if(abs(history[ih])==421) return kTRUE;
    if(abs(history[ih])==431) return kTRUE;
    if(abs(history[ih])==4122) return kTRUE;
    if(abs(history[ih])==4132) return kTRUE;
    if(abs(history[ih])==4232) return kTRUE;
    if(abs(history[ih])==4332) return kTRUE;
  }
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSEXic2eleXifromAODtracks::HaveBottomInHistory(Int_t *history){
  // 
  // Check if the particle is from Bottom hadrons
  // 
  for(Int_t ih=0;ih<10;ih++){
    if(abs(history[ih])==511) return kTRUE;
    if(abs(history[ih])==521) return kTRUE;
    if(abs(history[ih])==531) return kTRUE;
    if(abs(history[ih])==5122) return kTRUE;
    if(abs(history[ih])==5132) return kTRUE;
    if(abs(history[ih])==5232) return kTRUE;
    if(abs(history[ih])==5332) return kTRUE;
  }
  return kFALSE;
}
//____________________________________________________________________________
TProfile* AliAnalysisTaskSEXic2eleXifromAODtracks::GetEstimatorHistogram(const AliVEvent* event){
	/// Get Estimator Histogram from period event->GetRunNumber();
	///
	/// If you select SPD tracklets in |eta|<1 you should use type == 1
	///

	Int_t runNo  = event->GetRunNumber();
	Int_t period = -1;   // pp: 0-LHC10b, 1-LHC10c, 2-LHC10d, 3-LHC10e

	if(runNo>114930 && runNo<117223) period = 0;
	if(runNo>119158 && runNo<120830) period = 1;
	if(runNo>122373 && runNo<126438) period = 2;
	if(runNo>127711 && runNo<130851) period = 3;
	if(period<0 || period>3) return 0;


	return fMultEstimatorAvg[period];
}

//________________________________________________________________________
Float_t AliAnalysisTaskSEXic2eleXifromAODtracks::GetEventPlaneForCandidate(AliAODTrack* trk, AliAODcascade *casc,AliEventplane *pl, TVector2* q,TVector2* qsub1,TVector2* qsub2){
	//
  // remove autocorrelations (not implemented yet)
	//
 
  return fEventPlane;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSEXic2eleXifromAODtracks::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
  return result;
}
