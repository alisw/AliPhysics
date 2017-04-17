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
//               Lc->e Lambda  analysis code
//
//  Input: AOD
//  Output: TTree and/or THnSparse (mass vs pT vs Centrality)
//
//  Cuts:
//  TTree: SingleCuts on V0 and electron
//  THnSparse: In addition to that, IsSelected(obj, kCandidate) applied
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
#include "AliAnalysisTaskSELc2eleLambdafromAODtracks.h"
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
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSELc2eleLambdafromAODtracks);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSELc2eleLambdafromAODtracks::AliAnalysisTaskSELc2eleLambdafromAODtracks() : 
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
  fV0VariablesTree(0),
  fSingleVariablesTree(0),
  fMCVariablesTree(0),
  fMCEleVariablesTree(0),
  fMCV0VariablesTree(0),
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
  fCandidateV0Variables(),
  fCandidateSingleVariables(),
  fCandidateMCVariables(),
  fCandidateMCEleVariables(),
  fCandidateMCV0Variables(),
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
  fHistoEleLambdaMass(0),
  fHistoEleLambdaMassRS(0),
  fHistoEleLambdaMassRS1(0),
  fHistoEleLambdaMassRS2(0),
  fHistoEleLambdaMassWS(0),
  fHistoEleLambdaMassWS1(0),
  fHistoEleLambdaMassWS2(0),
  fHistoEleLambdaMassRSMix(0),
  fHistoEleLambdaMassRSMix1(0),
  fHistoEleLambdaMassRSMix2(0),
  fHistoEleLambdaMassWSMix(0),
  fHistoEleLambdaMassWSMix1(0),
  fHistoEleLambdaMassWSMix2(0),
  fHistoEleLambdaMassRSSide(0),
  fHistoEleLambdaMassRSSide1(0),
  fHistoEleLambdaMassRSSide2(0),
  fHistoEleLambdaMassWSSide(0),
  fHistoEleLambdaMassWSSide1(0),
  fHistoEleLambdaMassWSSide2(0),
  fHistoEleLambdaMassRSAway(0),
  fHistoEleLambdaMassRS1Away(0),
  fHistoEleLambdaMassRS2Away(0),
  fHistoEleLambdaMassWSAway(0),
  fHistoEleLambdaMassWS1Away(0),
  fHistoEleLambdaMassWS2Away(0),
  fHistoEleLambdaMassRSMixAway(0),
  fHistoEleLambdaMassRSMix1Away(0),
  fHistoEleLambdaMassRSMix2Away(0),
  fHistoEleLambdaMassWSMixAway(0),
  fHistoEleLambdaMassWSMix1Away(0),
  fHistoEleLambdaMassWSMix2Away(0),
  fHistoEleLambdaMassRSSideAway(0),
  fHistoEleLambdaMassRSSide1Away(0),
  fHistoEleLambdaMassRSSide2Away(0),
  fHistoEleLambdaMassWSSideAway(0),
  fHistoEleLambdaMassWSSide1Away(0),
  fHistoEleLambdaMassWSSide2Away(0),
  fHistoEleLambdaMassvsElePtRS(0),
  fHistoEleLambdaMassvsElePtWS(0),
  fHistoEleLambdaMassvsElePtRSMix(0),
  fHistoEleLambdaMassvsElePtWSMix(0),
  fHistoEleLambdaMassvsElePtRSSide(0),
  fHistoEleLambdaMassvsElePtWSSide(0),
  fHistoEleLambdaMassvsElePtRS1(0),
  fHistoEleLambdaMassvsElePtWS1(0),
  fHistoEleLambdaMassvsElePtRSMix1(0),
  fHistoEleLambdaMassvsElePtWSMix1(0),
  fHistoEleLambdaMassvsElePtRSSide1(0),
  fHistoEleLambdaMassvsElePtWSSide1(0),
  fHistoEleLambdaMassvsElePtRS2(0),
  fHistoEleLambdaMassvsElePtWS2(0),
  fHistoEleLambdaMassvsElePtRSMix2(0),
  fHistoEleLambdaMassvsElePtWSMix2(0),
  fHistoEleLambdaMassvsElePtRSSide2(0),
  fHistoEleLambdaMassvsElePtWSSide2(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleLambdaMassMCS(0),
  fHistoEleLambdaMassLambdabMCS(0),
  fHistoEleLambdaMassXibMCS(0),
  fHistoEleLambdaMassMCGen(0),
  fHistoEleLambdaMassvsElePtMCS(0),
  fHistoEleLambdaMassvsElePtMCGen(0),
  fHistoEleLambdaMassvsElePtMCS1(0),
  fHistoEleLambdaMassvsElePtMCGen1(0),
  fHistoEleLambdaMassvsElePtMCS2(0),
  fHistoEleLambdaMassvsElePtMCGen2(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsLambdaPtRS(0),
  fHistoElePtvsLambdaPtWS(0),
  fHistoElePtvsLambdaPtRSMix(0),
  fHistoElePtvsLambdaPtWSMix(0),
  fHistoElePtvsLambdaPtMCS(0),
  fHistoElePtvsLambdaPtvsLcPtMCS(0),
  fHistoElePtvsLambdaPtMCGen(0),
  fHistoElePtvsLambdaPtvsLcPtMCGen(0),
  fHistoElePtvsLambdaPtMCLcGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoEleLambdaMassFeeddownXic0MCS(0),
  fHistoEleLambdaMassFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS1(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS2(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2(0),
  fHistoElePtFeeddownXic0MCS(0),
  fHistoElePtFeeddownXic0MCGen(0),
  fHistoElePtvsEtaFeeddownXic0MCS(0),
  fHistoElePtvsEtaFeeddownXic0MCGen(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCS(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2(0),
  fHistoElePtFeeddownXicPlusMCS(0),
  fHistoElePtFeeddownXicPlusMCGen(0),
  fHistoElePtvsEtaFeeddownXicPlusMCS(0),
  fHistoElePtvsEtaFeeddownXicPlusMCGen(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistoBachFracSharedITS(0),
  fHistoBachChi2NclsITS(0),
  fHistod0Bach(0),
  fHistoLambdaMassvsPt(0),
  fHistoLambdaMassvsPtMCS(0),
  fHistoLambdaMassvsPtMCGen(0),
  fHistoLambdaPtvsEta(0),
  fHistoLambdaPtvsEtaMCS(0),
  fHistoLambdaPtvsEtaMCGen(0),
  fHistoLambdaArmenteros(0),
  fHistoLambdaArmenterosMCS(0),
  fHistoLambdaPtvsDl(0),
  fHistoLambdaPtvsDlSide(0),
  fHistoLambdaPtvsDlMCS(0),
  fHistoLambdaPtvsDlFeeddownXi0MCS(0),
  fHistoLambdaPtvsDlFeeddownXiMinusMCS(0),
  fHistoLambdaPtvsDlFeeddownOmegaMCS(0),
  fHistoK0sMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronITSPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronITSSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
  fHistoElectronITSPIDSelTPCTOF(0),
  fHistoElectronITSPIDSelTPCTOFSmallEta(0),
  fHistoElectronITSPIDSelTPCTOFLargeEta(0),
  fHistoMassConversionsMin(0),
  fHistoMassConversionsSameSignMin(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoLambdaQovPtvsPhi(0),
	fHistoLcMCGen(0),
	fHistoLcMCGen1(0),
	fHistoLcMCGen2(0),
	fHistoLbMCGen(0),
	fHistoLbMCGen1(0),
	fHistoLbMCGen2(0),
	fHistoLcMCGenWithEleA(0),
	fHistoLcMCGenWithEleB(0),
	fHistoLbMCGenWithEleA(0),
	fHistoLbMCGenWithEleB(0),
	fHistoLbMCGenWithLambdac(0),
	fHistoLcMCS(0),
	fHistoLcMCS1(0),
	fHistoLcMCS2(0),
	fHistoLbMCS(0),
	fHistoLbMCS1(0),
	fHistoLbMCS2(0),
	fHistoFeedDownXic0MCGen(0),
	fHistoFeedDownXic0MCGen1(0),
	fHistoFeedDownXic0MCGen2(0),
	fHistoFeedDownXic0MCS(0),
	fHistoFeedDownXic0MCS1(0),
	fHistoFeedDownXic0MCS2(0),
	fHistoFeedDownXicPlusMCGen(0),
	fHistoFeedDownXicPlusMCGen1(0),
	fHistoFeedDownXicPlusMCGen2(0),
	fHistoFeedDownXicPlusMCS(0),
	fHistoFeedDownXicPlusMCS1(0),
	fHistoFeedDownXicPlusMCS2(0),
	fHistoFeedDownSigmaMCGen(0),
	fHistoFeedDownSigmaMCGen1(0),
	fHistoFeedDownSigmaMCGen2(0),
	fHistoFeedDownSigmaMCS(0),
	fHistoFeedDownSigmaMCS1(0),
	fHistoFeedDownSigmaMCS2(0),
	fHistoLcElectronMCGen(0),
	fHistoLcElectronMCGen1(0),
	fHistoLcElectronMCGen2(0),
	fHistoLcElectronMCS(0),
	fHistoLcElectronMCS1(0),
	fHistoLcElectronMCS2(0),
	fHistoLbElectronMCGen(0),
	fHistoLbElectronMCGen1(0),
	fHistoLbElectronMCGen2(0),
	fHistoLbElectronMCS(0),
	fHistoLbElectronMCS1(0),
	fHistoLbElectronMCS2(0),
	fHistoElectronFeedDownXic0MCGen(0),
	fHistoElectronFeedDownXic0MCGen1(0),
	fHistoElectronFeedDownXic0MCGen2(0),
	fHistoElectronFeedDownXic0MCS(0),
	fHistoElectronFeedDownXic0MCS1(0),
	fHistoElectronFeedDownXic0MCS2(0),
	fHistoElectronFeedDownXicPlusMCGen(0),
	fHistoElectronFeedDownXicPlusMCGen1(0),
	fHistoElectronFeedDownXicPlusMCGen2(0),
	fHistoElectronFeedDownXicPlusMCS(0),
	fHistoElectronFeedDownXicPlusMCS1(0),
	fHistoElectronFeedDownXicPlusMCS2(0),
	fHistoElectronFeedDownXibMCGen(0),
	fHistoElectronFeedDownXibMCGen1(0),
	fHistoElectronFeedDownXibMCGen2(0),
	fHistoElectronFeedDownXibMCS(0),
	fHistoElectronFeedDownXibMCS1(0),
	fHistoElectronFeedDownXibMCS2(0),
	fHistoElectronMCGen(0),
	fHistoBottomElectronMCGen(0),
	fHistoCharmElectronMCGen(0),
	fHistoLambdaMCGen(0),
	fHistoElePtvsV0dlRS(0),
	fHistoElePtvsV0dlRS1(0),
	fHistoElePtvsV0dlRS2(0),
	fHistoElePtvsV0dlRSSide(0),
	fHistoElePtvsV0dlRSSide1(0),
	fHistoElePtvsV0dlRSSide2(0),
	fHistoElePtvsV0dlRSMix(0),
	fHistoElePtvsV0dlRSMix1(0),
	fHistoElePtvsV0dlRSMix2(0),
	fHistoElePtvsV0dlWS(0),
	fHistoElePtvsV0dlWS1(0),
	fHistoElePtvsV0dlWS2(0),
	fHistoElePtvsV0dlWSSide(0),
	fHistoElePtvsV0dlWSSide1(0),
	fHistoElePtvsV0dlWSSide2(0),
	fHistoElePtvsV0dlWSMix(0),
	fHistoElePtvsV0dlWSMix1(0),
	fHistoElePtvsV0dlWSMix2(0),
	fHistoElePtvsV0dlMCS(0),
	fHistoElePtvsV0dlMCS1(0),
	fHistoElePtvsV0dlMCS2(0),
	fHistoElePtvsV0dlFeedDownXic0MCS(0),
	fHistoElePtvsV0dlFeedDownXic0MCS1(0),
	fHistoElePtvsV0dlFeedDownXic0MCS2(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS1(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS2(0),
	fHistoElePtvsV0dcaRS(0),
	fHistoElePtvsV0dcaRS1(0),
	fHistoElePtvsV0dcaRS2(0),
	fHistoElePtvsV0dcaRSSide(0),
	fHistoElePtvsV0dcaRSSide1(0),
	fHistoElePtvsV0dcaRSSide2(0),
	fHistoElePtvsV0dcaRSMix(0),
	fHistoElePtvsV0dcaRSMix1(0),
	fHistoElePtvsV0dcaRSMix2(0),
	fHistoElePtvsV0dcaWS(0),
	fHistoElePtvsV0dcaWS1(0),
	fHistoElePtvsV0dcaWS2(0),
	fHistoElePtvsV0dcaWSSide(0),
	fHistoElePtvsV0dcaWSSide1(0),
	fHistoElePtvsV0dcaWSSide2(0),
	fHistoElePtvsV0dcaWSMix(0),
	fHistoElePtvsV0dcaWSMix1(0),
	fHistoElePtvsV0dcaWSMix2(0),
	fHistoElePtvsV0dcaMCS(0),
	fHistoElePtvsV0dcaMCS1(0),
	fHistoElePtvsV0dcaMCS2(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS1(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS2(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS1(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS2(0),
	fHistoEleLambdaPtvsV0dlRS(0),
	fHistoEleLambdaPtvsV0dlRS1(0),
	fHistoEleLambdaPtvsV0dlRS2(0),
	fHistoEleLambdaPtvsV0dlRSSide(0),
	fHistoEleLambdaPtvsV0dlRSSide1(0),
	fHistoEleLambdaPtvsV0dlRSSide2(0),
	fHistoEleLambdaPtvsV0dlRSMix(0),
	fHistoEleLambdaPtvsV0dlRSMix1(0),
	fHistoEleLambdaPtvsV0dlRSMix2(0),
	fHistoEleLambdaPtvsV0dlWS(0),
	fHistoEleLambdaPtvsV0dlWS1(0),
	fHistoEleLambdaPtvsV0dlWS2(0),
	fHistoEleLambdaPtvsV0dlWSSide(0),
	fHistoEleLambdaPtvsV0dlWSSide1(0),
	fHistoEleLambdaPtvsV0dlWSSide2(0),
	fHistoEleLambdaPtvsV0dlWSMix(0),
	fHistoEleLambdaPtvsV0dlWSMix1(0),
	fHistoEleLambdaPtvsV0dlWSMix2(0),
	fHistoEleLambdaPtvsV0dlMCS(0),
	fHistoEleLambdaPtvsV0dlMCS1(0),
	fHistoEleLambdaPtvsV0dlMCS2(0),
	fHistoEleLambdaPtvsV0dlPromptMCS(0),
	fHistoEleLambdaPtvsV0dlPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS1(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS2(0),
	fHistoEleLambdaPtvsV0dlXibMCS(0),
	fHistoEleLambdaPtvsV0dlXibMCS1(0),
	fHistoEleLambdaPtvsV0dlXibMCS2(0),
	fHistoEleLambdaPtvsV0dlRSAway(0),
	fHistoEleLambdaPtvsV0dlRS1Away(0),
	fHistoEleLambdaPtvsV0dlRS2Away(0),
	fHistoEleLambdaPtvsV0dlRSSideAway(0),
	fHistoEleLambdaPtvsV0dlRSSide1Away(0),
	fHistoEleLambdaPtvsV0dlRSSide2Away(0),
	fHistoEleLambdaPtvsV0dlRSMixAway(0),
	fHistoEleLambdaPtvsV0dlRSMix1Away(0),
	fHistoEleLambdaPtvsV0dlRSMix2Away(0),
	fHistoEleLambdaPtvsV0dlWSAway(0),
	fHistoEleLambdaPtvsV0dlWS1Away(0),
	fHistoEleLambdaPtvsV0dlWS2Away(0),
	fHistoEleLambdaPtvsV0dlWSSideAway(0),
	fHistoEleLambdaPtvsV0dlWSSide1Away(0),
	fHistoEleLambdaPtvsV0dlWSSide2Away(0),
	fHistoEleLambdaPtvsV0dlWSMixAway(0),
	fHistoEleLambdaPtvsV0dlWSMix1Away(0),
	fHistoEleLambdaPtvsV0dlWSMix2Away(0),
	fHistoCorrelationVariablesvsEleLambdaPt(0),
	fHistoCorrelationVariablesvsEleLambdaPtMix(0),
	fHistoCorrelationVariablesvsEleLambdaPtMC(0),
	fHistoCorrelationVariablesvsElePt(0),
	fHistoCorrelationVariablesvsElePtMix(0),
	fHistoCorrelationVariablesvsElePtMC(0),
	fHistoCorrelationVariablesvsLambdaPt(0),
	fHistoCorrelationVariablesvsLambdaPtMix(0),
	fHistoCorrelationVariablesvsLambdaPtMC(0),
	fHistoMassVariablesvsEleLambdaPt(0),
	fHistoMassVariablesvsEleLambdaPtMix(0),
	fHistoMassVariablesvsEleLambdaPtMC(0),
	fHistoMassVariablesvsElePt(0),
	fHistoMassVariablesvsElePtMix(0),
	fHistoMassVariablesvsElePtMC(0),
	fHistoMassVariablesvsLambdaPt(0),
	fHistoMassVariablesvsLambdaPtMix(0),
	fHistoMassVariablesvsLambdaPtMC(0),
	fHistoResponseElePt(0),
	fHistoResponseElePt1(0),
	fHistoResponseElePt2(0),
	fHistoResponseEleLambdaPt(0),
	fHistoResponseEleLambdaPt1(0),
	fHistoResponseEleLambdaPt2(0),
	fHistoResponseEleLambdaPtFeeddownXic0(0),
	fHistoResponseEleLambdaPtFeeddownXic01(0),
	fHistoResponseEleLambdaPtFeeddownXic02(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus1(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus2(0),
	fHistoResponseEleLambdaPtFeeddownSigma(0),
	fHistoResponseEleLambdaPtFeeddownSigma1(0),
	fHistoResponseEleLambdaPtFeeddownSigma2(0),
	fHistoResponseEleLambdaPtLambdab(0),
	fHistoResponseElePtLambdab(0),
	fHistoResponseLambdaPtLambdab(0),
	fHistoResponseLambdaPt(0),
	fHistoResponseLambdaPtFeeddownXic0(0),
	fHistoResponseLambdaPtFeeddownXicPlus(0),
	fHistoResponseLambdaPtFeeddownSigma(0),
	fHistoResponseLambdaPtvsEleLambdaPt(0),
	fHistoResponseLambdaPtvsEleLambdaPtLambdab(0),
	fHistoResponseMcGenLambdabPtvsLambdacPt(0),
	fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt(0),
	fHistoEleLambdaPtvsRapidityRS(0),
	fHistoEleLambdaPtvsRapidityWS(0),
	fHistoEleLambdaPtvsRapidityMCS(0),
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
  fHistonLambdavsRunNumber(0),
  fHistoMCEventType(0),
  fHistoMCDeltaPhiccbar(0),
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
  fDoEventMixing(0),
  fMixWithoutConversionFlag(kFALSE),
  fNumberOfEventsForMixing		(5),
  fNzVtxBins					(0), 
  fNCentBins					(0),
  fNRPBins					(0), 
  fNOfPools(1),
  fPoolIndex(-9999),
	fHistoPoolMonitor(0),
	fHistoPoolIDCounter(0),
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
  /// Default Constructor.
  //
	for(Int_t i=0;i<17;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
		fHistoElectronITSPIDSelTPCTOFEtaDep[i] = 0;
	}
	for(Int_t i=0; i<4; i++) fMultEstimatorAvg[i]=0;
}

//___________________________________________________________________________
AliAnalysisTaskSELc2eleLambdafromAODtracks::AliAnalysisTaskSELc2eleLambdafromAODtracks(const Char_t* name,
									     AliRDHFCutsLctoeleLambdafromAODtracks* analCuts, 
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
  fV0VariablesTree(0),
  fSingleVariablesTree(0),
  fMCVariablesTree(0),
  fMCEleVariablesTree(0),
  fMCV0VariablesTree(0),
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
  fCandidateV0Variables(),
  fCandidateSingleVariables(),
  fCandidateMCVariables(),
  fCandidateMCEleVariables(),
  fCandidateMCV0Variables(),
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
  fHistoEleLambdaMass(0),
  fHistoEleLambdaMassRS(0),
  fHistoEleLambdaMassRS1(0),
  fHistoEleLambdaMassRS2(0),
  fHistoEleLambdaMassWS(0),
  fHistoEleLambdaMassWS1(0),
  fHistoEleLambdaMassWS2(0),
  fHistoEleLambdaMassRSMix(0),
  fHistoEleLambdaMassRSMix1(0),
  fHistoEleLambdaMassRSMix2(0),
  fHistoEleLambdaMassWSMix(0),
  fHistoEleLambdaMassWSMix1(0),
  fHistoEleLambdaMassWSMix2(0),
  fHistoEleLambdaMassRSSide(0),
  fHistoEleLambdaMassRSSide1(0),
  fHistoEleLambdaMassRSSide2(0),
  fHistoEleLambdaMassWSSide(0),
  fHistoEleLambdaMassWSSide1(0),
  fHistoEleLambdaMassWSSide2(0),
  fHistoEleLambdaMassRSAway(0),
  fHistoEleLambdaMassRS1Away(0),
  fHistoEleLambdaMassRS2Away(0),
  fHistoEleLambdaMassWSAway(0),
  fHistoEleLambdaMassWS1Away(0),
  fHistoEleLambdaMassWS2Away(0),
  fHistoEleLambdaMassRSMixAway(0),
  fHistoEleLambdaMassRSMix1Away(0),
  fHistoEleLambdaMassRSMix2Away(0),
  fHistoEleLambdaMassWSMixAway(0),
  fHistoEleLambdaMassWSMix1Away(0),
  fHistoEleLambdaMassWSMix2Away(0),
  fHistoEleLambdaMassRSSideAway(0),
  fHistoEleLambdaMassRSSide1Away(0),
  fHistoEleLambdaMassRSSide2Away(0),
  fHistoEleLambdaMassWSSideAway(0),
  fHistoEleLambdaMassWSSide1Away(0),
  fHistoEleLambdaMassWSSide2Away(0),
  fHistoEleLambdaMassvsElePtRS(0),
  fHistoEleLambdaMassvsElePtWS(0),
  fHistoEleLambdaMassvsElePtRSMix(0),
  fHistoEleLambdaMassvsElePtWSMix(0),
  fHistoEleLambdaMassvsElePtRSSide(0),
  fHistoEleLambdaMassvsElePtWSSide(0),
  fHistoEleLambdaMassvsElePtRS1(0),
  fHistoEleLambdaMassvsElePtWS1(0),
  fHistoEleLambdaMassvsElePtRSMix1(0),
  fHistoEleLambdaMassvsElePtWSMix1(0),
  fHistoEleLambdaMassvsElePtRSSide1(0),
  fHistoEleLambdaMassvsElePtWSSide1(0),
  fHistoEleLambdaMassvsElePtRS2(0),
  fHistoEleLambdaMassvsElePtWS2(0),
  fHistoEleLambdaMassvsElePtRSMix2(0),
  fHistoEleLambdaMassvsElePtWSMix2(0),
  fHistoEleLambdaMassvsElePtRSSide2(0),
  fHistoEleLambdaMassvsElePtWSSide2(0),
  fHistoElePtRS(0),
  fHistoElePtWS(0),
  fHistoElePtRSMix(0),
  fHistoElePtWSMix(0),
  fHistoEleLambdaMassMCS(0),
  fHistoEleLambdaMassLambdabMCS(0),
  fHistoEleLambdaMassXibMCS(0),
  fHistoEleLambdaMassMCGen(0),
  fHistoEleLambdaMassvsElePtMCS(0),
  fHistoEleLambdaMassvsElePtMCGen(0),
  fHistoEleLambdaMassvsElePtMCS1(0),
  fHistoEleLambdaMassvsElePtMCGen1(0),
  fHistoEleLambdaMassvsElePtMCS2(0),
  fHistoEleLambdaMassvsElePtMCGen2(0),
  fHistoElePtMCS(0),
  fHistoElePtMCGen(0),
  fHistoElePtvsEtaRS(0),
  fHistoElePtvsEtaWS(0),
  fHistoElePtvsEtaRSMix(0),
  fHistoElePtvsEtaWSMix(0),
  fHistoElePtvsEtaMCS(0),
  fHistoElePtvsEtaMCGen(0),
  fHistoElePtvsLambdaPtRS(0),
  fHistoElePtvsLambdaPtWS(0),
  fHistoElePtvsLambdaPtRSMix(0),
  fHistoElePtvsLambdaPtWSMix(0),
  fHistoElePtvsLambdaPtMCS(0),
  fHistoElePtvsLambdaPtvsLcPtMCS(0),
  fHistoElePtvsLambdaPtMCGen(0),
  fHistoElePtvsLambdaPtvsLcPtMCGen(0),
  fHistoElePtvsLambdaPtMCLcGen(0),
  fHistoElePtvsd0RS(0),
  fHistoElePtvsd0WS(0),
  fHistoElePtvsd0RSMix(0),
  fHistoElePtvsd0WSMix(0),
  fHistoElePtvsd0MCS(0),
  fHistoElePtvsd0PromptMCS(0),
  fHistoElePtvsd0BFeeddownMCS(0),
  fHistoEleLambdaMassFeeddownXic0MCS(0),
  fHistoEleLambdaMassFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS1(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS2(0),
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2(0),
  fHistoElePtFeeddownXic0MCS(0),
  fHistoElePtFeeddownXic0MCGen(0),
  fHistoElePtvsEtaFeeddownXic0MCS(0),
  fHistoElePtvsEtaFeeddownXic0MCGen(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCS(0),
  fHistoElePtvsLambdaPtFeeddownXic0MCGen(0),
  fHistoEleLambdaMassFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2(0),
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2(0),
  fHistoElePtFeeddownXicPlusMCS(0),
  fHistoElePtFeeddownXicPlusMCGen(0),
  fHistoElePtvsEtaFeeddownXicPlusMCS(0),
  fHistoElePtvsEtaFeeddownXicPlusMCGen(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS(0),
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen(0),
  fHistoBachPt(0),
  fHistoBachPtMCS(0),
  fHistoBachPtMCGen(0),
  fHistoBachFracSharedITS(0),
  fHistoBachChi2NclsITS(0),
  fHistod0Bach(0),
  fHistoLambdaMassvsPt(0),
  fHistoLambdaMassvsPtMCS(0),
  fHistoLambdaMassvsPtMCGen(0),
  fHistoLambdaPtvsEta(0),
  fHistoLambdaPtvsEtaMCS(0),
  fHistoLambdaPtvsEtaMCGen(0),
  fHistoLambdaArmenteros(0),
  fHistoLambdaArmenterosMCS(0),
  fHistoLambdaPtvsDl(0),
  fHistoLambdaPtvsDlSide(0),
  fHistoLambdaPtvsDlMCS(0),
  fHistoLambdaPtvsDlFeeddownXi0MCS(0),
  fHistoLambdaPtvsDlFeeddownXiMinusMCS(0),
  fHistoLambdaPtvsDlFeeddownOmegaMCS(0),
  fHistoK0sMassvsPt(0),
  fHistoElectronTPCPID(0),
  fHistoElectronTOFPID(0),
  fHistoElectronITSPID(0),
  fHistoElectronTPCSelPID(0),
  fHistoElectronTOFSelPID(0),
  fHistoElectronITSSelPID(0),
  fHistoElectronTPCPIDSelTOF(0),
  fHistoElectronTPCPIDSelTOFSmallEta(0),
  fHistoElectronTPCPIDSelTOFLargeEta(0),
  fHistoElectronITSPIDSelTPCTOF(0),
  fHistoElectronITSPIDSelTPCTOFSmallEta(0),
  fHistoElectronITSPIDSelTPCTOFLargeEta(0),
  fHistoMassConversionsMin(0),
  fHistoMassConversionsSameSignMin(0),
	fHistoElectronQovPtvsPhi(0),
	fHistoLambdaQovPtvsPhi(0),
	fHistoLcMCGen(0),
	fHistoLcMCGen1(0),
	fHistoLcMCGen2(0),
	fHistoLbMCGen(0),
	fHistoLbMCGen1(0),
	fHistoLbMCGen2(0),
	fHistoLcMCGenWithEleA(0),
	fHistoLcMCGenWithEleB(0),
	fHistoLbMCGenWithEleA(0),
	fHistoLbMCGenWithEleB(0),
	fHistoLbMCGenWithLambdac(0),
	fHistoLcMCS(0),
	fHistoLcMCS1(0),
	fHistoLcMCS2(0),
	fHistoLbMCS(0),
	fHistoLbMCS1(0),
	fHistoLbMCS2(0),
	fHistoFeedDownXic0MCGen(0),
	fHistoFeedDownXic0MCGen1(0),
	fHistoFeedDownXic0MCGen2(0),
	fHistoFeedDownXic0MCS(0),
	fHistoFeedDownXic0MCS1(0),
	fHistoFeedDownXic0MCS2(0),
	fHistoFeedDownXicPlusMCGen(0),
	fHistoFeedDownXicPlusMCGen1(0),
	fHistoFeedDownXicPlusMCGen2(0),
	fHistoFeedDownXicPlusMCS(0),
	fHistoFeedDownXicPlusMCS1(0),
	fHistoFeedDownXicPlusMCS2(0),
	fHistoFeedDownSigmaMCGen(0),
	fHistoFeedDownSigmaMCGen1(0),
	fHistoFeedDownSigmaMCGen2(0),
	fHistoFeedDownSigmaMCS(0),
	fHistoFeedDownSigmaMCS1(0),
	fHistoFeedDownSigmaMCS2(0),
	fHistoLcElectronMCGen(0),
	fHistoLcElectronMCGen1(0),
	fHistoLcElectronMCGen2(0),
	fHistoLcElectronMCS(0),
	fHistoLcElectronMCS1(0),
	fHistoLcElectronMCS2(0),
	fHistoLbElectronMCGen(0),
	fHistoLbElectronMCGen1(0),
	fHistoLbElectronMCGen2(0),
	fHistoLbElectronMCS(0),
	fHistoLbElectronMCS1(0),
	fHistoLbElectronMCS2(0),
	fHistoElectronFeedDownXic0MCGen(0),
	fHistoElectronFeedDownXic0MCGen1(0),
	fHistoElectronFeedDownXic0MCGen2(0),
	fHistoElectronFeedDownXic0MCS(0),
	fHistoElectronFeedDownXic0MCS1(0),
	fHistoElectronFeedDownXic0MCS2(0),
	fHistoElectronFeedDownXicPlusMCGen(0),
	fHistoElectronFeedDownXicPlusMCGen1(0),
	fHistoElectronFeedDownXicPlusMCGen2(0),
	fHistoElectronFeedDownXicPlusMCS(0),
	fHistoElectronFeedDownXicPlusMCS1(0),
	fHistoElectronFeedDownXicPlusMCS2(0),
	fHistoElectronFeedDownXibMCGen(0),
	fHistoElectronFeedDownXibMCGen1(0),
	fHistoElectronFeedDownXibMCGen2(0),
	fHistoElectronFeedDownXibMCS(0),
	fHistoElectronFeedDownXibMCS1(0),
	fHistoElectronFeedDownXibMCS2(0),
	fHistoElectronMCGen(0),
	fHistoBottomElectronMCGen(0),
	fHistoCharmElectronMCGen(0),
	fHistoLambdaMCGen(0),
	fHistoElePtvsV0dlRS(0),
	fHistoElePtvsV0dlRS1(0),
	fHistoElePtvsV0dlRS2(0),
	fHistoElePtvsV0dlRSSide(0),
	fHistoElePtvsV0dlRSSide1(0),
	fHistoElePtvsV0dlRSSide2(0),
	fHistoElePtvsV0dlRSMix(0),
	fHistoElePtvsV0dlRSMix1(0),
	fHistoElePtvsV0dlRSMix2(0),
	fHistoElePtvsV0dlWS(0),
	fHistoElePtvsV0dlWS1(0),
	fHistoElePtvsV0dlWS2(0),
	fHistoElePtvsV0dlWSSide(0),
	fHistoElePtvsV0dlWSSide1(0),
	fHistoElePtvsV0dlWSSide2(0),
	fHistoElePtvsV0dlWSMix(0),
	fHistoElePtvsV0dlWSMix1(0),
	fHistoElePtvsV0dlWSMix2(0),
	fHistoElePtvsV0dlMCS(0),
	fHistoElePtvsV0dlMCS1(0),
	fHistoElePtvsV0dlMCS2(0),
	fHistoElePtvsV0dlFeedDownXic0MCS(0),
	fHistoElePtvsV0dlFeedDownXic0MCS1(0),
	fHistoElePtvsV0dlFeedDownXic0MCS2(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS1(0),
	fHistoElePtvsV0dlFeedDownXicPlusMCS2(0),
	fHistoElePtvsV0dcaRS(0),
	fHistoElePtvsV0dcaRS1(0),
	fHistoElePtvsV0dcaRS2(0),
	fHistoElePtvsV0dcaRSSide(0),
	fHistoElePtvsV0dcaRSSide1(0),
	fHistoElePtvsV0dcaRSSide2(0),
	fHistoElePtvsV0dcaRSMix(0),
	fHistoElePtvsV0dcaRSMix1(0),
	fHistoElePtvsV0dcaRSMix2(0),
	fHistoElePtvsV0dcaWS(0),
	fHistoElePtvsV0dcaWS1(0),
	fHistoElePtvsV0dcaWS2(0),
	fHistoElePtvsV0dcaWSSide(0),
	fHistoElePtvsV0dcaWSSide1(0),
	fHistoElePtvsV0dcaWSSide2(0),
	fHistoElePtvsV0dcaWSMix(0),
	fHistoElePtvsV0dcaWSMix1(0),
	fHistoElePtvsV0dcaWSMix2(0),
	fHistoElePtvsV0dcaMCS(0),
	fHistoElePtvsV0dcaMCS1(0),
	fHistoElePtvsV0dcaMCS2(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS1(0),
	fHistoElePtvsV0dcaFeedDownXic0MCS2(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS1(0),
	fHistoElePtvsV0dcaFeedDownXicPlusMCS2(0),
	fHistoEleLambdaPtvsV0dlRS(0),
	fHistoEleLambdaPtvsV0dlRS1(0),
	fHistoEleLambdaPtvsV0dlRS2(0),
	fHistoEleLambdaPtvsV0dlRSSide(0),
	fHistoEleLambdaPtvsV0dlRSSide1(0),
	fHistoEleLambdaPtvsV0dlRSSide2(0),
	fHistoEleLambdaPtvsV0dlRSMix(0),
	fHistoEleLambdaPtvsV0dlRSMix1(0),
	fHistoEleLambdaPtvsV0dlRSMix2(0),
	fHistoEleLambdaPtvsV0dlWS(0),
	fHistoEleLambdaPtvsV0dlWS1(0),
	fHistoEleLambdaPtvsV0dlWS2(0),
	fHistoEleLambdaPtvsV0dlWSSide(0),
	fHistoEleLambdaPtvsV0dlWSSide1(0),
	fHistoEleLambdaPtvsV0dlWSSide2(0),
	fHistoEleLambdaPtvsV0dlWSMix(0),
	fHistoEleLambdaPtvsV0dlWSMix1(0),
	fHistoEleLambdaPtvsV0dlWSMix2(0),
	fHistoEleLambdaPtvsV0dlMCS(0),
	fHistoEleLambdaPtvsV0dlMCS1(0),
	fHistoEleLambdaPtvsV0dlMCS2(0),
	fHistoEleLambdaPtvsV0dlPromptMCS(0),
	fHistoEleLambdaPtvsV0dlPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1(0),
	fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS1(0),
	fHistoEleLambdaPtvsV0dlLambdabMCS2(0),
	fHistoEleLambdaPtvsV0dlXibMCS(0),
	fHistoEleLambdaPtvsV0dlXibMCS1(0),
	fHistoEleLambdaPtvsV0dlXibMCS2(0),
	fHistoEleLambdaPtvsV0dlRSAway(0),
	fHistoEleLambdaPtvsV0dlRS1Away(0),
	fHistoEleLambdaPtvsV0dlRS2Away(0),
	fHistoEleLambdaPtvsV0dlRSSideAway(0),
	fHistoEleLambdaPtvsV0dlRSSide1Away(0),
	fHistoEleLambdaPtvsV0dlRSSide2Away(0),
	fHistoEleLambdaPtvsV0dlRSMixAway(0),
	fHistoEleLambdaPtvsV0dlRSMix1Away(0),
	fHistoEleLambdaPtvsV0dlRSMix2Away(0),
	fHistoEleLambdaPtvsV0dlWSAway(0),
	fHistoEleLambdaPtvsV0dlWS1Away(0),
	fHistoEleLambdaPtvsV0dlWS2Away(0),
	fHistoEleLambdaPtvsV0dlWSSideAway(0),
	fHistoEleLambdaPtvsV0dlWSSide1Away(0),
	fHistoEleLambdaPtvsV0dlWSSide2Away(0),
	fHistoEleLambdaPtvsV0dlWSMixAway(0),
	fHistoEleLambdaPtvsV0dlWSMix1Away(0),
	fHistoEleLambdaPtvsV0dlWSMix2Away(0),
	fHistoCorrelationVariablesvsEleLambdaPt(0),
	fHistoCorrelationVariablesvsEleLambdaPtMix(0),
	fHistoCorrelationVariablesvsEleLambdaPtMC(0),
	fHistoCorrelationVariablesvsElePt(0),
	fHistoCorrelationVariablesvsElePtMix(0),
	fHistoCorrelationVariablesvsElePtMC(0),
	fHistoCorrelationVariablesvsLambdaPt(0),
	fHistoCorrelationVariablesvsLambdaPtMix(0),
	fHistoCorrelationVariablesvsLambdaPtMC(0),
	fHistoMassVariablesvsEleLambdaPt(0),
	fHistoMassVariablesvsEleLambdaPtMix(0),
	fHistoMassVariablesvsEleLambdaPtMC(0),
	fHistoMassVariablesvsElePt(0),
	fHistoMassVariablesvsElePtMix(0),
	fHistoMassVariablesvsElePtMC(0),
	fHistoMassVariablesvsLambdaPt(0),
	fHistoMassVariablesvsLambdaPtMix(0),
	fHistoMassVariablesvsLambdaPtMC(0),
	fHistoResponseElePt(0),
	fHistoResponseElePt1(0),
	fHistoResponseElePt2(0),
	fHistoResponseEleLambdaPt(0),
	fHistoResponseEleLambdaPt1(0),
	fHistoResponseEleLambdaPt2(0),
	fHistoResponseEleLambdaPtFeeddownXic0(0),
	fHistoResponseEleLambdaPtFeeddownXic01(0),
	fHistoResponseEleLambdaPtFeeddownXic02(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus1(0),
	fHistoResponseEleLambdaPtFeeddownXicPlus2(0),
	fHistoResponseEleLambdaPtFeeddownSigma(0),
	fHistoResponseEleLambdaPtFeeddownSigma1(0),
	fHistoResponseEleLambdaPtFeeddownSigma2(0),
	fHistoResponseEleLambdaPtLambdab(0),
	fHistoResponseElePtLambdab(0),
	fHistoResponseLambdaPtLambdab(0),
	fHistoResponseLambdaPt(0),
	fHistoResponseLambdaPtFeeddownXic0(0),
	fHistoResponseLambdaPtFeeddownXicPlus(0),
	fHistoResponseLambdaPtFeeddownSigma(0),
	fHistoResponseLambdaPtvsEleLambdaPt(0),
	fHistoResponseLambdaPtvsEleLambdaPtLambdab(0),
	fHistoResponseMcGenLambdabPtvsLambdacPt(0),
	fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt(0),
	fHistoEleLambdaPtvsRapidityRS(0),
	fHistoEleLambdaPtvsRapidityWS(0),
	fHistoEleLambdaPtvsRapidityMCS(0),
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
	fHistonLambdavsRunNumber(0),
	fHistoMCEventType(0),
	fHistoMCDeltaPhiccbar(0),
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
  fDoEventMixing(0),
  fMixWithoutConversionFlag(kFALSE),
	fNumberOfEventsForMixing		(5),
	fNzVtxBins					(0), 
	fNCentBins					(0),
  fNRPBins					(0), 
	fNOfPools(1),
  fPoolIndex(-9999),
	fHistoPoolMonitor(0),
	fHistoPoolIDCounter(0),
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
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSELc2eleLambdafromAODtracks","Calling Constructor");

	for(Int_t i=0;i<17;i++){
		fHistoElePtvsCutVarsRS[i] = 0;
		fHistoElePtvsCutVarsWS[i] = 0;
		fHistoElePtvsCutVarsMCS[i] = 0;
	}
	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i] = 0;
		fHistoElectronITSPIDSelTPCTOFEtaDep[i] = 0;
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
AliAnalysisTaskSELc2eleLambdafromAODtracks::~AliAnalysisTaskSELc2eleLambdafromAODtracks() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSELc2eleLambdafromAODtracks","Calling Destructor");

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
  if (fV0VariablesTree) {
    delete fV0VariablesTree;
    fV0VariablesTree = 0;
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
  if (fMCV0VariablesTree) {
    delete fMCV0VariablesTree;
    fMCV0VariablesTree = 0;
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::Init() {
  //
  /// Initialization
  //
  //

  fIsEventSelected=kFALSE;

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsLctoeleLambdafromAODtracks(*fAnalCuts));
  PostData(2,fListCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec(Option_t *)
{
  //
  /// UserExec
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
      AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: MC header branch not found!\n");
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
  if(!fIsTriggerNotOK) fCEvents->Fill(3);
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
		fCentrality = fAnalCuts->GetCentrality(aodEvent);
  }else if(fUseCentralitySPDTracklet){
    if(countCorr>=0 && countCorr<=0) fCentrality = 5.;
    else if(countCorr>=1 && countCorr<=8) fCentrality = 15.;
    else if(countCorr>= 9 && countCorr<=13) fCentrality = 25.;
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
		if(fUseEventPlane<10){
			//OLD framework

			AliEventplane *pl=aodEvent->GetEventplane();
			if(!pl){
				AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
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
					AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: no Q vectors");
					fCEvents->Fill(19);
					return;
				}
			}
		}else{
			//New framework
			AliQnCorrectionsManager *flowQnVectorMgr;
			AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
				dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
			if (flowQnVectorTask != NULL) {
				flowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
			}
			else {
				AliFatal("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!");
			}

			if(fUseEventPlane==11){
				AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: No such event plane method");
				return;
			}else if(fUseEventPlane==12){
				AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: No such event plane method");
				return;
			}else if(fUseEventPlane==13){
				AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: No such event plane method");
				return;
			}else if(fUseEventPlane==14){
				const AliQnCorrectionsQnVector *myQnVector;
				myQnVector = flowQnVectorMgr->GetDetectorQnVector("TPC");
				if(!myQnVector) {
					AliError("AliAnalysisTaskSELc2eleLambdafromAODtracks::UserExec: no Q vectors");
					fCEvents->Fill(19);
					return;
				}
				Double_t ep_tpc = myQnVector->EventPlane(2);
				fEventPlane = GetPhi0Pi((Float_t)ep_tpc);
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
  PostData(6,fV0VariablesTree);
  PostData(7,fMCVariablesTree);
  PostData(8,fCounter);    
  PostData(9,fMCEleVariablesTree);
  PostData(10,fMCV0VariablesTree);
  //PostData(11,fMCGenPairVariablesTree);
  PostData(11,fSingleVariablesTree);
  PostData(12,fCorrelationVariablesTree);

  fIsEventSelected=kFALSE;

  delete fV1;
  return;
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::Terminate(Option_t*)
{    
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.
  
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::UserCreateOutputObjects() 
{ 
  ///
  /// UserCreateOutputObject
  ///
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

  DefineV0TreeVariables();
  PostData(6,fV0VariablesTree);

  DefineMCTreeVariables();
  PostData(7,fMCVariablesTree);

  DefineMCEleTreeVariables();
  PostData(9,fMCEleVariablesTree);

  DefineMCV0TreeVariables();
  PostData(10,fMCV0VariablesTree);

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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeAnalysis
(
 AliAODEvent *aodEvent, TClonesArray *mcArray
 )
{
  //
  /// Main Analysis part
  //
  if(fDoEventMixing){
    fPoolIndex=GetPoolIndex(fVtxZ,fCentrality,fEventPlane);
		fHistoPoolIDCounter->Fill(fPoolIndex);
		if(fPoolIndex<0) return;

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

  //------------------------------------------------
  // Select good track before hand to save time
  //------------------------------------------------

  Int_t nV0s= aodEvent->GetNumberOfV0s();
  Int_t nTracks= aodEvent->GetNumberOfTracks();

  Bool_t  seleTrkFlags[nTracks];
  Int_t nSeleTrks=0;
  SelectTrack(aodEvent,nTracks,nSeleTrks,seleTrkFlags,mcArray);

  Bool_t  seleV0Flags[nV0s];
  Int_t     nSeleV0=0;
  SelectV0(aodEvent,nV0s,nSeleV0,seleV0Flags,mcArray);

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
	fHistonLambdavsRunNumber->Fill(runnumber-runnumber_offset,nSeleV0);

  if(nSeleTrks==0 || nSeleV0 ==0) return;

  //------------------------------------------------
  // Fill Pools and single Root Objects
  //------------------------------------------------
  for (Int_t itrk = 0; itrk<nTracks; itrk++) {
    if(!seleTrkFlags[itrk]) continue;
    AliAODTrack *trk = (AliAODTrack*)aodEvent->GetTrack(itrk);

    FillElectronROOTObjects(trk,aodEvent,mcArray);
  }

  for (Int_t iv0 = 0; iv0<nV0s; iv0++) {
    if(!seleV0Flags[iv0]) continue;
    AliAODv0 *v0 = aodEvent->GetV0(iv0);
    if(!v0) continue;

    FillV0ROOTObjects(v0, aodEvent, mcArray);
  }

  if(fWriteEachVariableTree)
    return;

  //------------------------------------------------
  // V0 loop 
  //------------------------------------------------
  for (Int_t iv0 = 0; iv0<nV0s; iv0++) {
    if(!seleV0Flags[iv0]) continue;
    AliAODv0 *v0 = aodEvent->GetV0(iv0);
    if(!v0) continue;

    AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
    AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));

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
      Int_t lpid = trkpid->GetID();
      if((cpid==lpid)||(cnid==lpid)) continue;

      AliAODVertex *secVert = ReconstructSecondaryVertex(v0,trk,aodEvent);//Fake, prim vertex is just used as secondary vertex. place holder for future
      if(!secVert) continue;

      AliAODRecoCascadeHF *elobj = MakeCascadeHF(v0,trk,trkpid,aodEvent,secVert);
      if(!elobj) {
	continue;
      }

      FillROOTObjects(elobj,v0,trk,trkpid,aodEvent,mcArray);

      elobj->GetSecondaryVtx()->RemoveDaughters();
      elobj->UnsetOwnPrimaryVtx();
      delete elobj;elobj=NULL;
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineTreeVariables() 
{
  ///
  /// Define tree variables
  ///

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fVariablesTree = new TTree(nameoutput,"Candidates variables tree");
  Int_t nVar = 94;
  fCandidateVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="InvMassEleLambda";
  fCandidateVariableNames[ 2]="EleLambdaPt";
  fCandidateVariableNames[ 3]="EleLambdaPx";
  fCandidateVariableNames[ 4]="EleLambdaPy";
  fCandidateVariableNames[ 5]="EleLambdaPz";
  fCandidateVariableNames[ 6]="ElePx";
  fCandidateVariableNames[ 7]="ElePy";
  fCandidateVariableNames[ 8]="ElePz";
  fCandidateVariableNames[ 9]="V0Px";
  fCandidateVariableNames[10]="V0Py";
  fCandidateVariableNames[11]="V0Pz";
  fCandidateVariableNames[12]="AntiLambdaFlag";
  fCandidateVariableNames[13]="MassLambda";
  fCandidateVariableNames[14]="MassAntiLambda";
  fCandidateVariableNames[15]="Eled0";
  fCandidateVariableNames[16]="V0d0";
  fCandidateVariableNames[17]="nSigmaTPCele";
  fCandidateVariableNames[18]="nSigmaTOFele";
  fCandidateVariableNames[19]="nSigmaTPCv0pr";
  fCandidateVariableNames[20]="nSigmaTOFv0pr";
  fCandidateVariableNames[21]="EleCharge";
  fCandidateVariableNames[22]="ProtonPx";
  fCandidateVariableNames[23]="ProtonPy";
  fCandidateVariableNames[24]="ProtonPz";
  fCandidateVariableNames[25]="PiPx";
  fCandidateVariableNames[26]="PiPy";
  fCandidateVariableNames[27]="PiPz";
  fCandidateVariableNames[28]="mcpdglc";
  fCandidateVariableNames[29]="mclablc";
  fCandidateVariableNames[30]="mcpdgmomele";
  fCandidateVariableNames[31]="mcpdgmomv0";
  fCandidateVariableNames[32]="Mixing";
  fCandidateVariableNames[33]="mcpdgele";
  fCandidateVariableNames[34]="nSigmaTPCpr_etrk";
  fCandidateVariableNames[35]="nSigmaTOFpr_etrk";
  fCandidateVariableNames[36]="nSigmaTPCka_etrk";
  fCandidateVariableNames[37]="nSigmaTOFka_etrk";
  fCandidateVariableNames[38]="MassK0Short";
  fCandidateVariableNames[39]="mcpdggrmomele";
  fCandidateVariableNames[40]="mcpdggrmomv0";
  fCandidateVariableNames[41]="mcngenele";
  fCandidateVariableNames[42]="mcngenv0";
  fCandidateVariableNames[43]="mclcpx";
  fCandidateVariableNames[44]="mclcpy";
  fCandidateVariableNames[45]="mclcpz";
  fCandidateVariableNames[46]="mcelepx";
  fCandidateVariableNames[47]="mcelepy";
  fCandidateVariableNames[48]="mcelepz";
  fCandidateVariableNames[49]="mcv0px";
  fCandidateVariableNames[50]="mcv0py";
  fCandidateVariableNames[51]="mcv0pz";
  fCandidateVariableNames[52]="nSigmaTPCpi_etrk";
  fCandidateVariableNames[53]="nSigmaTOFpi_etrk";
  fCandidateVariableNames[54]="PrimVertx";
  fCandidateVariableNames[55]="PrimVerty";
  fCandidateVariableNames[56]="PrimVertz";
  fCandidateVariableNames[57]="V0Vertx";
  fCandidateVariableNames[58]="V0Verty";
  fCandidateVariableNames[59]="V0Vertz";

  fCandidateVariableNames[60]="DcaV0PrToPrimVertex";
  fCandidateVariableNames[61]="DcaV0PiToPrimVertex";
  fCandidateVariableNames[62]="DcaV0daughters";
  fCandidateVariableNames[63]="V0CosPointingAngle";
  fCandidateVariableNames[64]="V0ProperDecayLength";
  fCandidateVariableNames[65]="MassK0Short2";

  fCandidateVariableNames[66]="nSigmaTPCv0pi";
  fCandidateVariableNames[67]="nSigmaTOFv0pi";

	fCandidateVariableNames[68]= "EleITSMatch";
	fCandidateVariableNames[69]= "V0PosITSMatch";
	fCandidateVariableNames[70]= "V0NegITSMatch";
	fCandidateVariableNames[71]= "IsV0PeakRegion";
	fCandidateVariableNames[72]= "mcpdgv0";

	fCandidateVariableNames[73]= "mcv0pospx";
	fCandidateVariableNames[74]= "mcv0pospy";
	fCandidateVariableNames[75]= "mcv0pospz";
	fCandidateVariableNames[76]= "mcpdgv0pos";
	fCandidateVariableNames[77]= "mcv0negpx";
	fCandidateVariableNames[78]= "mcv0negpy";
	fCandidateVariableNames[79]= "mcv0negpz";
	fCandidateVariableNames[80]= "mcpdgv0neg";

	fCandidateVariableNames[81]= "ProtonPxPV";
	fCandidateVariableNames[82]= "ProtonPyPV";
	fCandidateVariableNames[83]= "ProtonPzPV";
	fCandidateVariableNames[84]= "PiPxPV";
	fCandidateVariableNames[85]= "PiPyPV";
	fCandidateVariableNames[86]= "PiPzPV";

	fCandidateVariableNames[87]= "mcv0vertx";
	fCandidateVariableNames[88]= "mcv0verty";
	fCandidateVariableNames[89]= "mcv0vertz";

  fCandidateVariableNames[90]="OnFlyStatus";
  fCandidateVariableNames[91]="MagneticField";
  fCandidateVariableNames[92]="EvNumber";
  fCandidateVariableNames[93]="RunNumber";


  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillROOTObjects(AliAODRecoCascadeHF *elobj, AliAODv0 *v0, AliAODTrack *trk, AliAODTrack *trkpid, AliAODEvent *event, TClonesArray *mcArray) 
{
  ///
  /// Fill histograms or tree depending on fWriteVariableTree
  ///
	if(!trk) return;
	if(!trkpid) return;
	if(!v0) return;

	for(Int_t i=0;i<94;i++)
  {
		fCandidateVariables[i] = -9999.;
	}

	Bool_t anti_lambda_flag = kFALSE;
	if(fabs(v0->MassAntiLambda()-1.115683)<fAnalCuts->GetProdV0MassTolLambdaRough()) anti_lambda_flag = kTRUE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
//	if(cptrack->Charge()<0 && cntrack->Charge()>0){
//		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
//		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
//	}

	Double_t v0px = elobj->PxProng(1);
	Double_t v0py = elobj->PyProng(1);
	Double_t v0pz = elobj->PzProng(1);
	Double_t momv0 = sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz);
	Double_t Ev0 = sqrt(momv0*momv0+1.115683*1.115683);
	Double_t epx = elobj->PxProng(0);
	Double_t epy = elobj->PyProng(0);
	Double_t epz = elobj->PzProng(0);
	Double_t mome = sqrt(epx*epx+epy*epy+epz*epz);
	Double_t Ee = sqrt(mome*mome+0.000510998928*0.000510998928);
  Double_t cosoa =  (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);
	Double_t Esum = Ee + Ev0;

	Double_t uxe = epx/mome;
	Double_t uye = epy/mome;
	Double_t uze = epz/mome;
	Double_t lf = -2.*(v0px*uxe+v0py*uye+v0pz*uze);
	Double_t pxv_flip = v0px + lf * uxe;
	Double_t pyv_flip = v0py + lf * uye;
	Double_t pzv_flip = v0pz + lf * uze;
	Double_t pxsum_flip = epx + pxv_flip;
	Double_t pysum_flip = epy + pyv_flip;
	Double_t pzsum_flip = epz + pzv_flip;
	Double_t melam_flip = sqrt(Esum*Esum-pxsum_flip*pxsum_flip-pysum_flip*pysum_flip-pzsum_flip*pzsum_flip);
	Double_t ptelam_flip = sqrt(pxsum_flip*pxsum_flip+pysum_flip*pysum_flip);

  Double_t minmass_ee = 9999.;
  Bool_t isconv = fAnalCuts->TagConversions(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmass_ee);
  Double_t minmasslike_ee = 9999.;
  Bool_t isconv_like = fAnalCuts->TagConversionsSameSign(trk,fGTIndex,(AliAODEvent*)event,event->GetNumberOfTracks(),minmasslike_ee);

  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3122};
  fCandidateVariables[ 1] = elobj->InvMass(2,pdgdg);
  fCandidateVariables[ 2] = elobj->Pt();
  fCandidateVariables[ 3] = elobj->Px();
  fCandidateVariables[ 4] = elobj->Py();
  fCandidateVariables[ 5] = elobj->Pz();
  fCandidateVariables[ 6] = elobj->PxProng(0);
  fCandidateVariables[ 7] = elobj->PyProng(0);
  fCandidateVariables[ 8] = elobj->PzProng(0);
  fCandidateVariables[ 9] = elobj->PxProng(1);
  fCandidateVariables[10] = elobj->PyProng(1);
  fCandidateVariables[11] = elobj->PzProng(1);
  fCandidateVariables[12] = anti_lambda_flag;
  fCandidateVariables[13] = v0->MassLambda();
  fCandidateVariables[14] = v0->MassAntiLambda();
  fCandidateVariables[15] = elobj->Getd0Prong(0);
  fCandidateVariables[16] = elobj->Getd0Prong(1);

  Double_t nSigmaTPCele=-9999.;
  Double_t nSigmaTOFele=-9999.;
  Double_t nSigmaTPCv0pr=-9999.;
  Double_t nSigmaTOFv0pr=-9999.;
  Double_t nSigmaTPCv0pi=-9999.;
  Double_t nSigmaTOFv0pi=-9999.;
  if(fAnalCuts->GetIsUsePID())
  {
		nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
		nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);
    fCandidateVariables[17] = nSigmaTPCele;
    fCandidateVariables[18] = nSigmaTOFele;
  }

	if(fAnalCuts->GetUseLambdaPID())
	{
		if(anti_lambda_flag){
			nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
		}else{
			nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
			nSigmaTOFv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
			nSigmaTPCv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
			nSigmaTOFv0pi = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
		}
      fCandidateVariables[19] = nSigmaTPCv0pr;
      fCandidateVariables[20] = nSigmaTOFv0pr;
      fCandidateVariables[66] = nSigmaTPCv0pi;
      fCandidateVariables[67] = nSigmaTOFv0pi;
  }
  fCandidateVariables[21] = trk->Charge();

	if(anti_lambda_flag){
		fCandidateVariables[22] = v0->MomNegX(); //cntrack->Px();
		fCandidateVariables[23] = v0->MomNegY(); //cntrack->Py();
		fCandidateVariables[24] = v0->MomNegZ(); //cntrack->Pz();
		fCandidateVariables[25] = v0->MomPosX(); //cptrack->Px();
		fCandidateVariables[26] = v0->MomPosY(); //cptrack->Py();
		fCandidateVariables[27] = v0->MomPosZ(); //cptrack->Pz();
	}else{
		fCandidateVariables[22] =  v0->MomPosX(); //cptrack->Px();
		fCandidateVariables[23] =  v0->MomPosY(); //cptrack->Py();
		fCandidateVariables[24] =  v0->MomPosZ(); //cptrack->Pz();
		fCandidateVariables[25] =  v0->MomNegX(); //cntrack->Px();
		fCandidateVariables[26] =  v0->MomNegY(); //cntrack->Py();
		fCandidateVariables[27] =  v0->MomNegZ(); //cntrack->Pz();
	}

  AliAODMCParticle *mclc = 0;
  AliAODMCParticle *mcele = 0;
  AliAODMCParticle *mcv0 = 0;
  Int_t mclablc = 0;
	Int_t mcpdgele_array[100];
	Int_t mcpdgv0_array[100];
	Int_t mclabelele_array[100];
	Int_t mclabelv0_array[100];
	Int_t mcngen_ele=-9999;
	Int_t mcngen_v0=-9999;

	if(fUseMCInfo && mcArray){
    mclablc =  MatchToMC(elobj,mcArray,mcpdgele_array, mcpdgv0_array,mclabelele_array,mclabelv0_array,mcngen_ele,mcngen_v0);
    if(mclabelele_array[0]>=0){
      mcele = (AliAODMCParticle*) mcArray->At(mclabelele_array[0]);
			if(mcele){
				fCandidateVariables[46] = mcele->Px();
				fCandidateVariables[47] = mcele->Py();
				fCandidateVariables[48] = mcele->Pz();
			}
    }
    if(mclabelv0_array[0]>=0){
      mcv0 = (AliAODMCParticle*) mcArray->At(mclabelv0_array[0]);
			if(mcv0){
				fCandidateVariables[49] = mcv0->Px();
				fCandidateVariables[50] = mcv0->Py();
				fCandidateVariables[51] = mcv0->Pz();
				fCandidateVariables[72] = mcv0->GetPdgCode();
			}
    }else{
      Int_t pdgdgv0[2]={211,211};
      Int_t labV0 = v0->MatchToMC(310,mcArray,2,pdgdgv0); // the V0
      if(labV0>=0){
        mcv0 = (AliAODMCParticle*) mcArray->At(labV0);
        if(mcv0){
          fCandidateVariables[49] = mcv0->Px();
          fCandidateVariables[50] = mcv0->Py();
          fCandidateVariables[51] = mcv0->Pz();
          fCandidateVariables[72] = mcv0->GetPdgCode();
        }
      }
    }

    Int_t labPos = cptrack->GetLabel();
    if(labPos>=0){
      AliAODMCParticle *mcv0pos = (AliAODMCParticle*)mcArray->At(labPos);
      fCandidateVariables[73] = mcv0pos->Px();
      fCandidateVariables[74] = mcv0pos->Py();
      fCandidateVariables[75] = mcv0pos->Pz();
      fCandidateVariables[76] = mcv0pos->GetPdgCode();
      fCandidateVariables[87] = mcv0pos->Xv();
      fCandidateVariables[88] = mcv0pos->Yv();
      fCandidateVariables[89] = mcv0pos->Zv();
    }
    Int_t labNeg = cntrack->GetLabel();
    if(labNeg>=0){
      AliAODMCParticle *mcv0neg = (AliAODMCParticle*)mcArray->At(labNeg);
      fCandidateVariables[77] = mcv0neg->Px();
      fCandidateVariables[78] = mcv0neg->Py();
      fCandidateVariables[79] = mcv0neg->Pz();
      fCandidateVariables[80] = mcv0neg->GetPdgCode();
    }

    if(mclablc>-1){
      mclc = (AliAODMCParticle*) mcArray->At(mclablc);
			if(mclc){
				fCandidateVariables[28] = mclc->GetPdgCode();
				fCandidateVariables[29] = mclc->Label();
				fCandidateVariables[43] = mclc->Px();
				fCandidateVariables[44] = mclc->Py();
				fCandidateVariables[45] = mclc->Pz();
			}
		}
		fCandidateVariables[30] = mcpdgele_array[1];
		fCandidateVariables[31] = mcpdgv0_array[1];
		fCandidateVariables[33] = mcpdgele_array[0];
		fCandidateVariables[39] = mcpdgele_array[2];
		fCandidateVariables[40] = mcpdgv0_array[2];
		fCandidateVariables[41] = mcngen_ele;
		fCandidateVariables[42] = mcngen_v0;
	}
	fCandidateVariables[32] = 0;

  if(fAnalCuts->GetIsUsePID())
  {
		Double_t nSigmaTPCpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kProton);
		Double_t nSigmaTOFpr_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kProton);
		Double_t nSigmaTPCka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kKaon);
		Double_t nSigmaTOFka_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kKaon);
		Double_t nSigmaTPCpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kPion);
		Double_t nSigmaTOFpi_etrk = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kPion);
    fCandidateVariables[34] = nSigmaTPCpr_etrk;
    fCandidateVariables[35] = nSigmaTOFpr_etrk;
    fCandidateVariables[36] = nSigmaTPCka_etrk;
    fCandidateVariables[37] = nSigmaTOFka_etrk;
    fCandidateVariables[52] = nSigmaTPCpi_etrk;
    fCandidateVariables[53] = nSigmaTOFpi_etrk;
  }
  fCandidateVariables[38] = v0->MassK0Short();

  fCandidateVariables[54] = fVtx1->GetX();
  fCandidateVariables[55] = fVtx1->GetY();
  fCandidateVariables[56] = fVtx1->GetZ();
  fCandidateVariables[57] = v0->DecayVertexV0X();
  fCandidateVariables[58] = v0->DecayVertexV0Y();
  fCandidateVariables[59] = v0->DecayVertexV0Z();

	Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
	Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
  if(!anti_lambda_flag){
		fCandidateVariables[60] = lDcaPosToPrimVertex;
		fCandidateVariables[61] = lDcaNegToPrimVertex;
  }else{
		fCandidateVariables[60] = lDcaNegToPrimVertex;
		fCandidateVariables[61] = lDcaPosToPrimVertex;
  }
	fCandidateVariables[62] = v0->DcaV0Daughters();
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);
  fCandidateVariables[63] = v0->CosPointingAngle(posVtx); 
  Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
  fCandidateVariables[64] = v0->DecayLengthV0(posVtx)*1.1157/ptotlam;
  fCandidateVariables[65] = v0->MassK0Short();

	if(trk) fCandidateVariables[68] = trk->GetITSClusterMap();
	if(cptrack) fCandidateVariables[69] = cptrack->GetITSClusterMap();
	if(cntrack) fCandidateVariables[70] = cntrack->GetITSClusterMap();

  fCandidateVariables[71] = fAnalCuts->IsPeakRegion(v0);

  if(!anti_lambda_flag){
    fCandidateVariables[81] = cptrack->Px();
    fCandidateVariables[82] = cptrack->Py();
    fCandidateVariables[83] = cptrack->Pz();
    fCandidateVariables[84] = cntrack->Px();
    fCandidateVariables[85] = cntrack->Py();
    fCandidateVariables[86] = cntrack->Pz();
  }else{
    fCandidateVariables[81] = cntrack->Px();
    fCandidateVariables[82] = cntrack->Py();
    fCandidateVariables[83] = cntrack->Pz();
    fCandidateVariables[84] = cptrack->Px();
    fCandidateVariables[85] = cptrack->Py();
    fCandidateVariables[86] = cptrack->Pz();
  }
  fCandidateVariables[90] = v0->GetOnFlyStatus();

  fCandidateVariables[91] = fBzkG;
  fCandidateVariables[92] = fEvNumberCounter;
  fCandidateVariables[93] = fRunNumber;


//  if(fWriteVariableTree)
//    fVariablesTree->Fill();

  Double_t dphis_ele_pr, detas_ele_pr,dphis_ele_pi, detas_ele_pi;
  dphis_ele_pr = 9999.;detas_ele_pr = 9999.;dphis_ele_pi = 9999.;detas_ele_pi = 9999.;
  //fAnalCuts->GetdPhiSdEtaSR125(trk,cptrack,cntrack,fBzkG,posVtx, dphis_ele_pr,detas_ele_pr,dphis_ele_pi,detas_ele_pi);

	Double_t cont[3];
	cont[0] = elobj->InvMass(2,pdgdg);
	cont[1] = elobj->Pt();
	cont[2] = fCentrality;
	fHistoEleLambdaMass->Fill(cont);

	Double_t cont_flip[3];
	cont_flip[0] = melam_flip;
	cont_flip[1] = ptelam_flip;
	cont_flip[2] = fCentrality;

	Double_t cont2[3];
	cont2[0] = elobj->InvMass(2,pdgdg);
	cont2[1] = trk->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = trk->Pt();
	cont_eleptvseta[1] = trk->Eta();
	cont_eleptvseta[2] = fCentrality;

	Double_t cont_eleptvslambdapt[3];
	cont_eleptvslambdapt[0] = trk->Pt();
	cont_eleptvslambdapt[1] = v0->Pt();
	cont_eleptvslambdapt[2] = fCentrality;

	Double_t cont_eleptvsd0[3];
	cont_eleptvsd0[0] = trk->Pt();
	cont_eleptvsd0[1] = elobj->Getd0Prong(0)*trk->Charge();
	cont_eleptvsd0[2] = fCentrality;

	Double_t cont_eleptvsv0dl[3];
	cont_eleptvsv0dl[0] = trk->Pt();
	cont_eleptvsv0dl[1] = v0->DecayLengthV0(posVtx)*1.115683/ptotlam;
	cont_eleptvsv0dl[2] = fCentrality;

	Double_t cont_elelamptvsv0dl[4];
	cont_elelamptvsv0dl[0] = elobj->Pt();
	cont_elelamptvsv0dl[1] = v0->DecayLengthV0(posVtx)*1.115683/ptotlam;
	cont_elelamptvsv0dl[2] = elobj->Getd0Prong(0)*trk->Charge();
	cont_elelamptvsv0dl[3] = fCentrality;

	Double_t cont_elelamptvsv0dl_flip[4];
	cont_elelamptvsv0dl_flip[0] = ptelam_flip;
	cont_elelamptvsv0dl_flip[1] = v0->DecayLengthV0(posVtx)*1.115683/ptotlam;
	cont_elelamptvsv0dl_flip[2] = 0.0;
	cont_elelamptvsv0dl_flip[3] = fCentrality;


	Double_t cont_eleptvsv0dca[3];
	cont_eleptvsv0dca[0] = trk->Pt();
	cont_eleptvsv0dca[1] = v0->DcaV0ToPrimVertex();
	cont_eleptvsv0dca[2] = fCentrality;

	Double_t elobj_mass = elobj->InvMass(2,pdgdg);
	Double_t elobj_px = elobj->Px();
	Double_t elobj_py = elobj->Py();
	Double_t elobj_pz = elobj->Pz();
	Double_t elobj_E = sqrt(elobj_mass*elobj_mass+elobj_px*elobj_px+elobj_py*elobj_py+elobj_pz*elobj_pz);
	Double_t elobj_rap = 0.5*log((elobj_E+elobj_pz)/(elobj_E-elobj_pz));

	//
	// Old strategy only look at mass and pairpt
	//
	if(fAnalCuts->IsSelected(elobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsPeakRegion(v0))
	{
		if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
			fHistoEleLambdaMassRS->Fill(cont);
			if(trk->Charge()>0) fHistoEleLambdaMassRS1->Fill(cont);
			else fHistoEleLambdaMassRS2->Fill(cont);

			fHistoEleLambdaMassvsElePtRS->Fill(cont2);
			if(trk->Charge()>0) fHistoEleLambdaMassvsElePtRS1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtRS2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsRapidityRS->Fill(elobj->Pt(),elobj_rap);
				fHistoElePtRS->Fill(trk->Pt(),fCentrality);
				fHistoElePtvsEtaRS->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtRS->Fill(cont_eleptvslambdapt);
				fHistoElePtvsd0RS->Fill(cont_eleptvsd0);

				fHistoElePtvsV0dlRS->Fill(cont_eleptvsv0dl);
				if(trk->Charge()>0) fHistoElePtvsV0dlRS1->Fill(cont_eleptvsv0dl);
				else fHistoElePtvsV0dlRS2->Fill(cont_eleptvsv0dl);

				fHistoEleLambdaPtvsV0dlRS->Fill(cont_elelamptvsv0dl);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlRS1->Fill(cont_elelamptvsv0dl);
				else fHistoEleLambdaPtvsV0dlRS2->Fill(cont_elelamptvsv0dl);

				fHistoElePtvsV0dcaRS->Fill(cont_eleptvsv0dca);
				if(trk->Charge()>0) fHistoElePtvsV0dcaRS1->Fill(cont_eleptvsv0dca);
				else fHistoElePtvsV0dcaRS2->Fill(cont_eleptvsv0dca);

				for(Int_t ih=0;ih<17;ih++){
					Double_t cont_eleptvscutvars[3];
					cont_eleptvscutvars[0] = elobj->Pt();
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
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->MassLambda();
						else
							cont_eleptvscutvars[1] = v0->MassAntiLambda();
					}else if(ih==7){
						Double_t lPosV0[3];
						lPosV0[0] = v0->DecayVertexV0X();
						lPosV0[1] = v0->DecayVertexV0Y();
						lPosV0[2] = v0->DecayVertexV0Z();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
					}else if(ih==8){
						cont_eleptvscutvars[1] = v0->DcaV0Daughters();
					}else if(ih==9){
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
					}else if(ih==10){
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
						else
							cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
					}else if(ih==11){
						cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
					}else if(ih==12){
						cont_eleptvscutvars[1] =  v0->MassK0Short();
					}else if(ih==13){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
					}else if(ih==14){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
					}else if(ih==15){
						cont_eleptvscutvars[1] =  v0->Eta();
					}else if(ih==16){
						Double_t v0px = elobj->PxProng(1);
						Double_t v0py = elobj->PyProng(1);
						Double_t v0pz = elobj->PzProng(1);
						Double_t epx = elobj->PxProng(0);
						Double_t epy = elobj->PyProng(0);
						Double_t epz = elobj->PzProng(0);
						cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
					}else{
						cont_eleptvscutvars[1] = -9999.;
					}

					fHistoElePtvsCutVarsRS[ih]->Fill(cont_eleptvscutvars);
				}
			}
			if(fUseMCInfo){
				if(mclc){
					Int_t pdgcode = mclc->GetPdgCode();
					Double_t cont_mclc[3];
					cont_mclc[0] = mclc->Pt();
					cont_mclc[1] = mclc->Y();
					cont_mclc[2] = fCentrality;
					Double_t cont_mcele[3];
					cont_mcele[0] = mcele->Pt();
					cont_mcele[1] = mcele->Eta();
					cont_mcele[2] = fCentrality;

					Int_t labmotherlc = mclc->GetMother();
					Bool_t isbottomfd = kFALSE;
					if(labmotherlc>=0){
						AliAODMCParticle *motherlc = (AliAODMCParticle*)mcArray->At(labmotherlc);
						Int_t pdgmotherlc = motherlc->GetPdgCode();
						if(TMath::Abs(pdgmotherlc)==511||TMath::Abs(pdgmotherlc)==521||TMath::Abs(pdgmotherlc)==5122||TMath::Abs(pdgmotherlc)==5132||TMath::Abs(pdgmotherlc)==5232||TMath::Abs(pdgmotherlc)==5332){
							isbottomfd = kTRUE;
						}
					}
					cont2[1] = mcele->Pt();

					if(abs(pdgcode)==4122 && abs(mcpdgele_array[1])==4122 && abs(mcpdgv0_array[1])==4122){
						fHistoEleLambdaMassMCS->Fill(cont);
						fHistoEleLambdaMassvsElePtMCS->Fill(cont2);
						if(trk->Charge()>0) fHistoEleLambdaMassvsElePtMCS1->Fill(cont2);
						else fHistoEleLambdaMassvsElePtMCS2->Fill(cont2);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoEleLambdaPtvsRapidityMCS->Fill(elobj->Pt(),elobj_rap);
							fHistoElePtMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaMCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtMCS->Fill(cont_eleptvslambdapt);
							fHistoElePtvsd0MCS->Fill(cont_eleptvsd0);

							fHistoElePtvsV0dlMCS->Fill(cont_eleptvsv0dl);
							if(trk->Charge()>0) fHistoElePtvsV0dlMCS1->Fill(cont_eleptvsv0dl);
							else fHistoElePtvsV0dlMCS2->Fill(cont_eleptvsv0dl);

							fHistoElePtvsV0dcaMCS->Fill(cont_eleptvsv0dca);
							if(trk->Charge()>0) fHistoElePtvsV0dcaMCS1->Fill(cont_eleptvsv0dca);
							else fHistoElePtvsV0dcaMCS2->Fill(cont_eleptvsv0dca);

							fHistoEleLambdaPtvsV0dlMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlMCS2->Fill(cont_elelamptvsv0dl);

							if(isbottomfd){
								fHistoEleLambdaPtvsV0dlBFeeddownMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlBFeeddownMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlBFeeddownMCS2->Fill(cont_elelamptvsv0dl);
							}else{
								fHistoEleLambdaPtvsV0dlPromptMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlPromptMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlPromptMCS2->Fill(cont_elelamptvsv0dl);
							}

							fHistoLcMCS->Fill(cont_mclc);
							if(trk->Charge()>0) fHistoLcMCS1->Fill(cont_mclc);
							else fHistoLcMCS2->Fill(cont_mclc);

							fHistoLcElectronMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoLcElectronMCS1->Fill(cont_mcele);
							else fHistoLcElectronMCS2->Fill(cont_mcele);

							fHistoResponseElePt->Fill(mclc->Pt(),trk->Pt());
							if(trk->Charge()>0) fHistoResponseElePt1->Fill(mclc->Pt(),trk->Pt());
							else fHistoResponseElePt2->Fill(mclc->Pt(),trk->Pt());
							fHistoResponseEleLambdaPt->Fill(mclc->Pt(),elobj->Pt());
							if(trk->Charge()>0) fHistoResponseEleLambdaPt1->Fill(mclc->Pt(),trk->Pt());
							else fHistoResponseEleLambdaPt2->Fill(mclc->Pt(),trk->Pt());

							fHistoResponseLambdaPt->Fill(mclc->Pt(),v0->Pt());
							fHistoResponseLambdaPtvsEleLambdaPt->Fill(elobj->Pt(),v0->Pt());

							Double_t cont_eleptvslambdaptvslcpt[4];
							cont_eleptvslambdaptvslcpt[0] = cont_eleptvslambdapt[0];
							cont_eleptvslambdaptvslcpt[1] = cont_eleptvslambdapt[1];
							cont_eleptvslambdaptvslcpt[2] = mclc->Pt();
							cont_eleptvslambdaptvslcpt[3] = cont_eleptvslambdapt[2];
							fHistoElePtvsLambdaPtvsLcPtMCS->Fill(cont_eleptvslambdaptvslcpt);

							Double_t cont_allpt[4];
							cont_allpt[0] = mclc->Pt();
							cont_allpt[1] = elobj->Pt();
							cont_allpt[2] = trk->Pt();
							cont_allpt[3] = v0->Pt();
							fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt->Fill(cont_allpt);

							if(isbottomfd){
								fHistoElePtvsd0BFeeddownMCS->Fill(cont_eleptvsd0);
							}else{
								fHistoElePtvsd0PromptMCS->Fill(cont_eleptvsd0);
							}

							for(Int_t ih=0;ih<17;ih++){
								Double_t cont_eleptvscutvars[3];
								cont_eleptvscutvars[0] = elobj->Pt();
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
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->MassLambda();
									else
										cont_eleptvscutvars[1] = v0->MassAntiLambda();
								}else if(ih==7){
									Double_t lPosV0[3];
									lPosV0[0] = v0->DecayVertexV0X();
									lPosV0[1] = v0->DecayVertexV0Y();
									lPosV0[2] = v0->DecayVertexV0Z();
									cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
								}else if(ih==8){
									cont_eleptvscutvars[1] = v0->DcaV0Daughters();
								}else if(ih==9){
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
									else
										cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
								}else if(ih==10){
									if(!anti_lambda_flag)
										cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
									else
										cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
								}else if(ih==11){
									cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
								}else if(ih==12){
									cont_eleptvscutvars[1] =  v0->MassK0Short();
								}else if(ih==13){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
								}else if(ih==14){
									cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
								}else if(ih==15){
									cont_eleptvscutvars[1] =  v0->Eta();
								}else if(ih==16){
									Double_t v0px = elobj->PxProng(1);
									Double_t v0py = elobj->PyProng(1);
									Double_t v0pz = elobj->PzProng(1);
									Double_t epx = elobj->PxProng(0);
									Double_t epy = elobj->PyProng(0);
									Double_t epz = elobj->PzProng(0);
									cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
								}else{
									cont_eleptvscutvars[1] = -9999.;
								}

								fHistoElePtvsCutVarsMCS[ih]->Fill(cont_eleptvscutvars);
							}
						}
					}
					if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgv0_array[1])==3312){
						fHistoEleLambdaMassFeeddownXic0MCS->Fill(cont);
						fHistoEleLambdaMassvsElePtFeeddownXic0MCS->Fill(cont2);
						if(trk->Charge()>0) fHistoEleLambdaMassvsElePtFeeddownXic0MCS1->Fill(cont2);
						else fHistoEleLambdaMassvsElePtFeeddownXic0MCS2->Fill(cont2);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoFeedDownXic0MCS->Fill(cont_mclc);
							if(trk->Charge()>0) fHistoFeedDownXic0MCS1->Fill(cont_mclc);
							else fHistoFeedDownXic0MCS2->Fill(cont_mclc);

							fHistoElectronFeedDownXic0MCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXic0MCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXic0MCS2->Fill(cont_mcele);

							fHistoElePtvsV0dlFeedDownXic0MCS->Fill(cont_eleptvsv0dl);
							if(trk->Charge()>0) fHistoElePtvsV0dlFeedDownXic0MCS1->Fill(cont_eleptvsv0dl);
							else fHistoElePtvsV0dlFeedDownXic0MCS2->Fill(cont_eleptvsv0dl);

							fHistoEleLambdaPtvsV0dlFeedDownXic0MCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2->Fill(cont_elelamptvsv0dl);

							if(isbottomfd){
								fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2->Fill(cont_elelamptvsv0dl);
							}else{
								fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2->Fill(cont_elelamptvsv0dl);
							}

							fHistoElePtvsV0dcaFeedDownXic0MCS->Fill(cont_eleptvsv0dca);
							if(trk->Charge()>0) fHistoElePtvsV0dcaFeedDownXic0MCS->Fill(cont_eleptvsv0dca);
							else fHistoElePtvsV0dcaFeedDownXic0MCS->Fill(cont_eleptvsv0dca);

							fHistoResponseEleLambdaPtFeeddownXic0->Fill(mclc->Pt(),elobj->Pt());
							if(trk->Charge()>0) fHistoResponseEleLambdaPtFeeddownXic01->Fill(mclc->Pt(),trk->Pt());
							else fHistoResponseEleLambdaPtFeeddownXic02->Fill(mclc->Pt(),trk->Pt());

							fHistoResponseLambdaPtFeeddownXic0->Fill(mclc->Pt(),v0->Pt());

							fHistoElePtFeeddownXic0MCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaFeeddownXic0MCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtFeeddownXic0MCS->Fill(cont_eleptvslambdapt);
						}
					}
					if(abs(pdgcode)==4232 && abs(mcpdgele_array[1])==4232 && abs(mcpdgv0_array[1])==3322){
						fHistoEleLambdaMassFeeddownXicPlusMCS->Fill(cont);
						fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS->Fill(cont2);
						if(trk->Charge()>0) fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1->Fill(cont2);
						else fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2->Fill(cont2);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoFeedDownXicPlusMCS->Fill(cont_mclc);
							if(trk->Charge()>0) fHistoFeedDownXicPlusMCS1->Fill(cont_mclc);
							else fHistoFeedDownXicPlusMCS2->Fill(cont_mclc);

							fHistoElectronFeedDownXicPlusMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXicPlusMCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXicPlusMCS2->Fill(cont_mcele);

							fHistoElePtvsV0dlFeedDownXicPlusMCS->Fill(cont_eleptvsv0dl);
							if(trk->Charge()>0) fHistoElePtvsV0dlFeedDownXicPlusMCS1->Fill(cont_eleptvsv0dl);
							else fHistoElePtvsV0dlFeedDownXicPlusMCS2->Fill(cont_eleptvsv0dl);

							fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2->Fill(cont_elelamptvsv0dl);

							if(isbottomfd){
								fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2->Fill(cont_elelamptvsv0dl);
							}else{
								fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2->Fill(cont_elelamptvsv0dl);
							}

							fHistoElePtvsV0dcaFeedDownXicPlusMCS->Fill(cont_eleptvsv0dca);
							if(trk->Charge()>0) fHistoElePtvsV0dcaFeedDownXicPlusMCS->Fill(cont_eleptvsv0dca);
							else fHistoElePtvsV0dcaFeedDownXicPlusMCS->Fill(cont_eleptvsv0dca);

							fHistoResponseEleLambdaPtFeeddownXicPlus->Fill(mclc->Pt(),elobj->Pt());
							if(trk->Charge()>0) fHistoResponseEleLambdaPtFeeddownXicPlus1->Fill(mclc->Pt(),trk->Pt());
							else fHistoResponseEleLambdaPtFeeddownXicPlus2->Fill(mclc->Pt(),trk->Pt());

							fHistoResponseLambdaPtFeeddownXicPlus->Fill(mclc->Pt(),v0->Pt());

							fHistoElePtFeeddownXicPlusMCS->Fill(trk->Pt(),fCentrality);
							fHistoElePtvsEtaFeeddownXicPlusMCS->Fill(cont_eleptvseta);
							fHistoElePtvsLambdaPtFeeddownXicPlusMCS->Fill(cont_eleptvslambdapt);
						}
					}
					if(abs(pdgcode)==4122 && abs(mcpdgele_array[1])==4122 && ( abs(mcpdgv0_array[1])==3212 || abs(mcpdgv0_array[1])==3214 || abs(mcpdgv0_array[1])==3224) ){
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoFeedDownSigmaMCS->Fill(cont_mclc);
							if(trk->Charge()>0) fHistoFeedDownSigmaMCS1->Fill(cont_mclc);
							else fHistoFeedDownSigmaMCS2->Fill(cont_mclc);

							fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2->Fill(cont_elelamptvsv0dl);

							if(isbottomfd){
								fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2->Fill(cont_elelamptvsv0dl);
							}else{
								fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS->Fill(cont_elelamptvsv0dl);
								if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1->Fill(cont_elelamptvsv0dl);
								else fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2->Fill(cont_elelamptvsv0dl);
							}
							fHistoResponseEleLambdaPtFeeddownSigma->Fill(mclc->Pt(),elobj->Pt());
							if(trk->Charge()>0) fHistoResponseEleLambdaPtFeeddownSigma1->Fill(mclc->Pt(),trk->Pt());
							else fHistoResponseEleLambdaPtFeeddownSigma2->Fill(mclc->Pt(),trk->Pt());

							fHistoResponseLambdaPtFeeddownSigma->Fill(mclc->Pt(),v0->Pt());
						}
					}
				}
			}
      fHistodPhiSdEtaSElectronProtonR125RS->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125RS->Fill(dphis_ele_pi,detas_ele_pi);
		}else if((trk->Charge()<0 && !anti_lambda_flag) || (trk->Charge()>0 && anti_lambda_flag)){
			fHistoEleLambdaMassWS->Fill(cont);
			if(trk->Charge()>0) fHistoEleLambdaMassWS1->Fill(cont);
			else fHistoEleLambdaMassWS2->Fill(cont);

			fHistoEleLambdaMassvsElePtWS->Fill(cont2);
			if(trk->Charge()>0) fHistoEleLambdaMassvsElePtWS1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtWS2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsRapidityWS->Fill(elobj->Pt(),elobj_rap);
				fHistoElePtWS->Fill(trk->Pt(),fCentrality);
				fHistoElePtvsEtaWS->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtWS->Fill(cont_eleptvslambdapt);
				fHistoElePtvsd0WS->Fill(cont_eleptvsd0);

				fHistoElePtvsV0dlWS->Fill(cont_eleptvsv0dl);
				if(trk->Charge()>0) fHistoElePtvsV0dlWS1->Fill(cont_eleptvsv0dl);
				else fHistoElePtvsV0dlWS2->Fill(cont_eleptvsv0dl);

				fHistoEleLambdaPtvsV0dlWS->Fill(cont_elelamptvsv0dl);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlWS1->Fill(cont_elelamptvsv0dl);
				else fHistoEleLambdaPtvsV0dlWS2->Fill(cont_elelamptvsv0dl);

				fHistoElePtvsV0dcaWS->Fill(cont_eleptvsv0dca);
				if(trk->Charge()>0) fHistoElePtvsV0dcaWS1->Fill(cont_eleptvsv0dca);
				else fHistoElePtvsV0dcaWS2->Fill(cont_eleptvsv0dca);

				for(Int_t ih=0;ih<17;ih++){
					Double_t cont_eleptvscutvars[3];
					cont_eleptvscutvars[0] = elobj->Pt();
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
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->MassLambda();
						else
							cont_eleptvscutvars[1] = v0->MassAntiLambda();
					}else if(ih==7){
						Double_t lPosV0[3];
						lPosV0[0] = v0->DecayVertexV0X();
						lPosV0[1] = v0->DecayVertexV0Y();
						lPosV0[2] = v0->DecayVertexV0Z();
						cont_eleptvscutvars[1] = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
					}else if(ih==8){
						cont_eleptvscutvars[1] = v0->DcaV0Daughters();
					}else if(ih==9){
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
						else
							cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
					}else if(ih==10){
						if(!anti_lambda_flag)
							cont_eleptvscutvars[1] = v0->DcaNegToPrimVertex();
						else
							cont_eleptvscutvars[1] = v0->DcaPosToPrimVertex();
					}else if(ih==11){
						cont_eleptvscutvars[1] =  v0->CosPointingAngle(posVtx);
					}else if(ih==12){
						cont_eleptvscutvars[1] =  v0->MassK0Short();
					}else if(ih==13){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pr;
					}else if(ih==14){
						cont_eleptvscutvars[1] =  nSigmaTPCv0pi;
					}else if(ih==15){
						cont_eleptvscutvars[1] =  v0->Eta();
					}else if(ih==16){
						Double_t v0px = elobj->PxProng(1);
						Double_t v0py = elobj->PyProng(1);
						Double_t v0pz = elobj->PzProng(1);
						Double_t epx = elobj->PxProng(0);
						Double_t epy = elobj->PyProng(0);
						Double_t epz = elobj->PzProng(0);
						cont_eleptvscutvars[1] = acos((v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz));
					}else{
						cont_eleptvscutvars[1] = -9999.;
					}

					fHistoElePtvsCutVarsWS[ih]->Fill(cont_eleptvscutvars);
				}
			}
			if(fUseMCInfo){
				if(mclc){
					Int_t pdgcode = mclc->GetPdgCode();
					Double_t cont_mclb[3];
					cont_mclb[0] = mclc->Pt();
					cont_mclb[1] = mclc->Y();
					cont_mclb[2] = fCentrality;
					Double_t cont_mcele[3];
					cont_mcele[0] = mcele->Pt();
					cont_mcele[1] = mcele->Eta();
					cont_mcele[2] = fCentrality;

					if(abs(pdgcode)==5122 && abs(mcpdgele_array[1])==5122 && abs(mcpdgv0_array[1])==4122 && abs(mcpdgv0_array[2])==5122){
						fHistoEleLambdaMassLambdabMCS->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoResponseEleLambdaPtLambdab->Fill(mclc->Pt(),elobj->Pt());
							fHistoResponseElePtLambdab->Fill(mclc->Pt(),trk->Pt());
							fHistoResponseLambdaPtLambdab->Fill(mclc->Pt(),v0->Pt());
							fHistoResponseLambdaPtvsEleLambdaPtLambdab->Fill(elobj->Pt(),v0->Pt());

							fHistoLbMCS->Fill(cont_mclb);
							if(trk->Charge()>0) fHistoLbMCS1->Fill(cont_mclb);
							else fHistoLbMCS2->Fill(cont_mclb);

							fHistoEleLambdaPtvsV0dlLambdabMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlLambdabMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlLambdabMCS2->Fill(cont_elelamptvsv0dl);

							fHistoLbElectronMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoLbElectronMCS1->Fill(cont_mcele);
							else fHistoLbElectronMCS2->Fill(cont_mcele);
						}
          }
					if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgv0_array[1])==4132 && abs(mcpdgv0_array[2])==5132){
						fHistoEleLambdaMassXibMCS->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoEleLambdaPtvsV0dlXibMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlXibMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlXibMCS2->Fill(cont_elelamptvsv0dl);

							fHistoElectronFeedDownXibMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXibMCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXibMCS2->Fill(cont_mcele);
						}
          }
					if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgv0_array[1])==4232 && abs(mcpdgv0_array[2])==5132){
						fHistoEleLambdaMassXibMCS->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoEleLambdaPtvsV0dlXibMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlXibMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlXibMCS2->Fill(cont_elelamptvsv0dl);

							fHistoElectronFeedDownXibMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXibMCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXibMCS2->Fill(cont_mcele);
						}
          }
					if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgv0_array[1])==4132 && abs(mcpdgv0_array[2])==5232){
						fHistoEleLambdaMassXibMCS->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoEleLambdaPtvsV0dlXibMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlXibMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlXibMCS2->Fill(cont_elelamptvsv0dl);

							fHistoElectronFeedDownXibMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXibMCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXibMCS2->Fill(cont_mcele);
						}
          }
					if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgv0_array[1])==4232 && abs(mcpdgv0_array[2])==5232){
						fHistoEleLambdaMassXibMCS->Fill(cont);
						if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
							fHistoEleLambdaPtvsV0dlXibMCS->Fill(cont_elelamptvsv0dl);
							if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlXibMCS1->Fill(cont_elelamptvsv0dl);
							else fHistoEleLambdaPtvsV0dlXibMCS2->Fill(cont_elelamptvsv0dl);

							fHistoElectronFeedDownXibMCS->Fill(cont_mcele);
							if(trk->Charge()>0) fHistoElectronFeedDownXibMCS1->Fill(cont_mcele);
							else fHistoElectronFeedDownXibMCS2->Fill(cont_mcele);
						}
          }
				}
			}
      fHistodPhiSdEtaSElectronProtonR125WS->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125WS->Fill(dphis_ele_pi,detas_ele_pi);
		}
	}

	//if( elobj->InvMass(2,pdgdg) < 10. && cosoa < 0. && fAnalCuts->IsPeakRegion(v0))
	if( melam_flip < 10. && cosoa < 0. && fAnalCuts->IsPeakRegion(v0))
	{
		if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
			fHistoEleLambdaMassRSAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleLambdaMassRS1Away->Fill(cont_flip);
			else fHistoEleLambdaMassRS2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlRSAway->Fill(cont_elelamptvsv0dl_flip);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlRS1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlRS2Away->Fill(cont_elelamptvsv0dl_flip);
			}
		}else{
			fHistoEleLambdaMassWSAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleLambdaMassWS1Away->Fill(cont_flip);
			else fHistoEleLambdaMassWS2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlWSAway->Fill(cont_elelamptvsv0dl_flip);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlWS1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlWS2Away->Fill(cont_elelamptvsv0dl_flip);
			}
		}
	}

	if(fAnalCuts->IsSelected(elobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsSideBand(v0))
	{
		if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
			fHistoEleLambdaMassRSSide->Fill(cont);
			if(trk->Charge()>0) fHistoEleLambdaMassRSSide1->Fill(cont);
			else fHistoEleLambdaMassRSSide2->Fill(cont);

			fHistoEleLambdaMassvsElePtRSSide->Fill(cont2);
			if(trk->Charge()>0) fHistoEleLambdaMassvsElePtRSSide1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtRSSide2->Fill(cont2);

			fHistoElePtvsV0dlRSSide->Fill(cont_eleptvsv0dl);
			if(trk->Charge()>0) fHistoElePtvsV0dlRSSide1->Fill(cont_eleptvsv0dl);
			else fHistoElePtvsV0dlRSSide2->Fill(cont_eleptvsv0dl);

			fHistoEleLambdaPtvsV0dlRSSide->Fill(cont_elelamptvsv0dl);
			if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlRSSide1->Fill(cont_elelamptvsv0dl);
			else fHistoEleLambdaPtvsV0dlRSSide2->Fill(cont_elelamptvsv0dl);

			fHistoElePtvsV0dcaRSSide->Fill(cont_eleptvsv0dca);
			if(trk->Charge()>0) fHistoElePtvsV0dcaRSSide1->Fill(cont_eleptvsv0dca);
			else fHistoElePtvsV0dcaRSSide2->Fill(cont_eleptvsv0dca);

		}else if((trk->Charge()<0 && !anti_lambda_flag) || (trk->Charge()>0 && anti_lambda_flag)){
			fHistoEleLambdaMassWSSide->Fill(cont);
			if(trk->Charge()>0) fHistoEleLambdaMassWSSide1->Fill(cont);
			else fHistoEleLambdaMassWSSide2->Fill(cont);

			fHistoEleLambdaMassvsElePtWSSide->Fill(cont2);
			if(trk->Charge()>0) fHistoEleLambdaMassvsElePtWSSide1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtWSSide2->Fill(cont2);

			fHistoElePtvsV0dlWSSide->Fill(cont_eleptvsv0dl);
			if(trk->Charge()>0) fHistoElePtvsV0dlWSSide1->Fill(cont_eleptvsv0dl);
			else fHistoElePtvsV0dlWSSide2->Fill(cont_eleptvsv0dl);

			fHistoEleLambdaPtvsV0dlWSSide->Fill(cont_elelamptvsv0dl);
			if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlWSSide1->Fill(cont_elelamptvsv0dl);
			else fHistoEleLambdaPtvsV0dlWSSide2->Fill(cont_elelamptvsv0dl);

			fHistoElePtvsV0dcaWSSide->Fill(cont_eleptvsv0dca);
			if(trk->Charge()>0) fHistoElePtvsV0dcaWSSide1->Fill(cont_eleptvsv0dca);
			else fHistoElePtvsV0dcaWSSide2->Fill(cont_eleptvsv0dca);
		}
	}

	//if( elobj->InvMass(2,pdgdg) < 10. && cosoa < 0. && fAnalCuts->IsSideBand(v0))
	if( melam_flip < 10. && cosoa < 0. && fAnalCuts->IsSideBand(v0))
	{
		if((trk->Charge()>0 && !anti_lambda_flag) || (trk->Charge()<0 && anti_lambda_flag)){
			fHistoEleLambdaMassRSSideAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleLambdaMassRSSide1Away->Fill(cont_flip);
			else fHistoEleLambdaMassRSSide2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlRSSideAway->Fill(cont_elelamptvsv0dl_flip);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlRSSide1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlRSSide2Away->Fill(cont_elelamptvsv0dl_flip);
			}
		}else{
			fHistoEleLambdaMassWSSideAway->Fill(cont_flip);
			if(trk->Charge()>0) fHistoEleLambdaMassWSSide1Away->Fill(cont_flip);
			else fHistoEleLambdaMassWSSide2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlWSSideAway->Fill(cont_elelamptvsv0dl_flip);
				if(trk->Charge()>0) fHistoEleLambdaPtvsV0dlWSSide1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlWSSide2Away->Fill(cont_elelamptvsv0dl_flip);
			}
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

  fCorrelationVariables[0] = v0->Pt();
  fCorrelationVariables[1] = trk->Pt();
  fCorrelationVariables[2] = fAnalCuts->DeltaPhi(v0,trk);
  fCorrelationVariables[3] = fAnalCuts->DeltaEta(v0,trk);
  fCorrelationVariables[4] = v0->DecayLengthV0(posVtx)*1.115683/ptotlam;;
  fCorrelationVariables[5] = elobj->Getd0Prong(0);
  if(!fUseMCInfo) fCorrelationVariables[6] = 0;
  else fCorrelationVariables[6] = 1;
  if(trk->Charge()>0){
    if(!anti_lambda_flag) fCorrelationVariables[7] = 0;
    else fCorrelationVariables[7] = 2;
  }else if(trk->Charge()<0){
    if(!anti_lambda_flag) fCorrelationVariables[7] = 1;
    else fCorrelationVariables[7] = 3;
  }
  fCorrelationVariables[8] = (Int_t)isconv + 2 * (Int_t)isconv_like;
  fCorrelationVariables[10] = fCentrality;
  fCorrelationVariables[11] = elobj->Pt();
	fCorrelationVariables[12] = elobj->InvMass(2,pdgdg);
	fCorrelationVariables[13] = fAnalCuts->CosOpeningAngle(v0,trk);

	cont_cor_nd[0] = elobj->Pt();
	cont_cor_nd[1] = fAnalCuts->DeltaPhi(v0,trk);
	cont_cor_nd[2] = v0->DecayLengthV0(posVtx)*1.115683/ptotlam;;
  if(trk->Charge()>0){
    if(!anti_lambda_flag) cont_cor_nd[3] = 0;
    else cont_cor_nd[3] = 2;
  }else if(trk->Charge()<0){
    if(!anti_lambda_flag) cont_cor_nd[3] = 3;
    else cont_cor_nd[3] = 1;
  }
	cont_cor_nd[4] = fCorrelationVariables[8];
	cont_cor_nd[5] = 0;
	cont_cor_nd[6] = fCentrality;

  if(fUseMCInfo && FromSemileptonicDecays(mcpdgele_array)>0){
    if(mclc){
      Int_t pdgcode = mclc->GetPdgCode();
      if(abs(pdgcode)==4122 && abs(mcpdgele_array[1])==4122 && abs(mcpdgv0_array[1])==4122){
        fCorrelationVariables[9] = 1;
				cont_cor_nd[5] = 1;
      }
      if(abs(pdgcode)==4122 && abs(mcpdgele_array[1])==4122 && ( abs(mcpdgv0_array[1])==3212 || abs(mcpdgv0_array[1])==3214 || abs(mcpdgv0_array[1])==3224) ){
        fCorrelationVariables[9] = 2;
				cont_cor_nd[5] = 2;
      }
      if(abs(pdgcode)==4132 && abs(mcpdgele_array[1])==4132 && abs(mcpdgv0_array[1])==3312){
        fCorrelationVariables[9] = 3;
				cont_cor_nd[5] = 3;
      }
      if(abs(pdgcode)==4232 && abs(mcpdgele_array[1])==4232 && abs(mcpdgv0_array[1])==3322){
        fCorrelationVariables[9] = 4;
				cont_cor_nd[5] = 4;
      }
      if(abs(pdgcode)==5122 && abs(mcpdgele_array[1])==5122 && abs(mcpdgv0_array[1])==4122 && abs(mcpdgv0_array[2])==5122){
        fCorrelationVariables[9] = 11;
				cont_cor_nd[5] = 5;
      }
      if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgv0_array[1])==4132 && abs(mcpdgv0_array[2])==5132){
        fCorrelationVariables[9] = 12;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5132 && abs(mcpdgele_array[1])==5132 && abs(mcpdgv0_array[1])==4232 && abs(mcpdgv0_array[2])==5132){
        fCorrelationVariables[9] = 13;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgv0_array[1])==4132 && abs(mcpdgv0_array[2])==5232){
        fCorrelationVariables[9] = 14;
				cont_cor_nd[5] = 6;
      }
      if(abs(pdgcode)==5232 && abs(mcpdgele_array[1])==5232 && abs(mcpdgv0_array[1])==4232 && abs(mcpdgv0_array[2])==5232){
        fCorrelationVariables[9] = 15;
				cont_cor_nd[5] = 6;
      }
      if(FromSemileptonicDecays(mcpdgele_array)==3){
        fCorrelationVariables[9] = 16;
      }
    }
    if(fCorrelationVariables[9]<0){
      Bool_t lam_from_bottom = HaveBottomInHistory(mcpdgv0_array);
      Bool_t lam_from_charm = HaveCharmInHistory(mcpdgv0_array);
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

	if(fAnalCuts->IsSelected(elobj,AliRDHFCuts::kCandidate) && fAnalCuts->IsPeakRegion(v0))
	{
		if(fWriteVariableTree)
			fCorrelationVariablesTree->Fill();

		if(fUseMCInfo){
			if(elobj->InvMass(2,pdgdg)<fAnalCuts->GetEleLambdaMassMax())
				fHistoCorrelationVariablesvsEleLambdaPtMC->Fill(cont_cor_nd);
			cont_cor_nd[0] = trk->Pt();
			fHistoCorrelationVariablesvsElePtMC->Fill(cont_cor_nd);
			cont_cor_nd[0] = v0->Pt();
			fHistoCorrelationVariablesvsLambdaPtMC->Fill(cont_cor_nd);
		}else{
			if(elobj->InvMass(2,pdgdg)<fAnalCuts->GetEleLambdaMassMax())
				fHistoCorrelationVariablesvsEleLambdaPt->Fill(cont_cor_nd);
			cont_cor_nd[0] = trk->Pt();
			fHistoCorrelationVariablesvsElePt->Fill(cont_cor_nd);
			cont_cor_nd[0] = v0->Pt();
			fHistoCorrelationVariablesvsLambdaPt->Fill(cont_cor_nd);
		}
  }

	cont_mass_nd[0] =  elobj->InvMass(2,pdgdg);
	cont_mass_nd[1] =  cont_cor_nd[0];
	cont_mass_nd[4] =  cont_cor_nd[3];
	cont_mass_nd[5] =  cont_cor_nd[4];
	cont_mass_nd[6] =  cont_cor_nd[5];
	cont_mass_nd[7] =  cont_cor_nd[6];

	if(fAnalCuts->IsPeakRegion(v0)) cont_mass_nd[3] = 1;
	if(fAnalCuts->IsSideBand(v0)) cont_mass_nd[3] = 0;
	if(fAnalCuts->IsSelected(elobj,AliRDHFCuts::kCandidate)) cont_mass_nd[2]=1;
	if(melam_flip < 10.&& cosoa < 0.) cont_mass_nd[2]=0;
	if(fUseMCInfo){
		fHistoMassVariablesvsEleLambdaPtMC->Fill(cont_mass_nd);
		cont_mass_nd[1] = trk->Pt();
		fHistoMassVariablesvsElePtMC->Fill(cont_mass_nd);
		cont_mass_nd[1] = v0->Pt();
		fHistoMassVariablesvsLambdaPtMC->Fill(cont_mass_nd);
	}else{
		fHistoMassVariablesvsEleLambdaPt->Fill(cont_mass_nd);
		cont_mass_nd[1] = trk->Pt();
		fHistoMassVariablesvsElePt->Fill(cont_mass_nd);
		cont_mass_nd[1] = v0->Pt();
		fHistoMassVariablesvsLambdaPt->Fill(cont_mass_nd);
	}

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMixROOTObjects(TLorentzVector *trke, TLorentzVector *v0, TVector *elevars, TVector *v0vars, Int_t chargepr) 
{
  ///
  /// Fill histograms or tree depending on fWriteVariableTree
  ///
	if(!trke) return;
	if(!v0) return;

	for(Int_t i=0;i<94;i++)
  {
		fCandidateVariables[i] = -9999.;
	}


	Double_t pxe = trke->Px();
	Double_t pye = trke->Py();
	Double_t pze = trke->Pz();
	Double_t mome = sqrt(pxe*pxe+pye*pye+pze*pze);
	Double_t Ee = sqrt(mome*mome+0.000510998928*0.000510998928);

	Double_t pxv = v0->Px();
	Double_t pyv = v0->Py();
	Double_t pzv = v0->Pz();
	Double_t momv = sqrt(pxv*pxv+pyv*pyv+pzv*pzv);
	Double_t Ev = sqrt(momv*momv+1.115683*1.115683);

	Double_t cosoa = (pxe*pxv+pye*pyv+pze*pzv)/mome/momv;

	Double_t pxsum = pxe + pxv;
	Double_t pysum = pye + pyv;
	Double_t pzsum = pze + pzv;
	Double_t Esum = Ee + Ev;
	Double_t mel = sqrt(Esum*Esum-pxsum*pxsum-pysum*pysum-pzsum*pzsum);

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
	Double_t mel_flip = sqrt(Esum*Esum-pxsum_flip*pxsum_flip-pysum_flip*pysum_flip-pzsum_flip*pzsum_flip);
	Double_t ptel_flip = sqrt(pxsum_flip*pxsum_flip+pysum_flip*pysum_flip);

  fCandidateVariables[ 0] = fCentrality;
	UInt_t pdgdg[2]={11,3122};
  fCandidateVariables[ 1] = mel;
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
	if(chargepr>0){
		fCandidateVariables[12] = 0;
		fCandidateVariables[13] = v0->M();
		fCandidateVariables[14] = 0.;
	}else{
		fCandidateVariables[12] = 1;
		fCandidateVariables[13] = 0;
		fCandidateVariables[14] = v0->M();
	}

  fCandidateVariables[17] = (*elevars)[3];
  fCandidateVariables[18] = (*elevars)[4];

  fCandidateVariables[21] = trke->T();
  fCandidateVariables[22] = (*v0vars)[9];
  fCandidateVariables[23] = (*v0vars)[10];
  fCandidateVariables[24] = (*v0vars)[11];
  fCandidateVariables[25] = (*v0vars)[12];
  fCandidateVariables[26] = (*v0vars)[13];
  fCandidateVariables[27] = (*v0vars)[14];

	fCandidateVariables[32] = 1;

  fCandidateVariables[54] = fVtx1->GetX();
  fCandidateVariables[55] = fVtx1->GetY();
  fCandidateVariables[56] = fVtx1->GetZ();
  fCandidateVariables[57] = (*v0vars)[6];
  fCandidateVariables[58] = (*v0vars)[7];
  fCandidateVariables[59] = (*v0vars)[8];
  fCandidateVariables[64] = (*v0vars)[15];

  fCandidateVariables[91] = fBzkG;
  fCandidateVariables[92] = fEvNumberCounter;
  fCandidateVariables[93] = fRunNumber;


//  if(fWriteVariableTree)
//    fVariablesTree->Fill();

	Double_t cont[3];
	cont[0] = mel;
	cont[1] = sqrt(pxsum*pxsum+pysum*pysum);
	cont[2] = fCentrality;
	fHistoEleLambdaMass->Fill(cont);

	Double_t cont_flip[3];
	cont_flip[0] = mel_flip;
	cont_flip[1] = ptel_flip;
	cont_flip[2] = fCentrality;

	Double_t cont2[3];
	cont2[0] = mel;
	cont2[1] = sqrt(pxe*pxe+pye*pye);
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = trke->Pt();
	cont_eleptvseta[1] = trke->Eta();
	cont_eleptvseta[2] = fCentrality;

	Double_t cont_eleptvslambdapt[3];
	cont_eleptvslambdapt[0] = trke->Pt();
	cont_eleptvslambdapt[1] = v0->Pt();
	cont_eleptvslambdapt[2] = fCentrality;

	Double_t cont_eleptvsd0[3];
	cont_eleptvsd0[0] = trke->Pt();
	cont_eleptvsd0[1] = 0.;
	cont_eleptvsd0[2] = fCentrality;

	Double_t cont_eleptvsv0dl[3];
	cont_eleptvsv0dl[0] = trke->Pt();
	cont_eleptvsv0dl[1] = (*v0vars)[15];
	cont_eleptvsv0dl[2] = fCentrality;

	Double_t cont_eleptvsv0dca[3];
	cont_eleptvsv0dca[0] = trke->Pt();
	cont_eleptvsv0dca[1] = (*v0vars)[16];
	cont_eleptvsv0dca[2] = fCentrality;

	Double_t cont_elelamptvsv0dl[4];
	cont_elelamptvsv0dl[0] = sqrt(pxsum*pxsum+pysum*pysum);
	cont_elelamptvsv0dl[1] = (*v0vars)[15];
	cont_elelamptvsv0dl[2] = 0.0;
	cont_elelamptvsv0dl[3] = fCentrality;

	Double_t cont_elelamptvsv0dl_flip[4];
	cont_elelamptvsv0dl_flip[0] = ptel_flip;
	cont_elelamptvsv0dl_flip[1] = (*v0vars)[15];
	cont_elelamptvsv0dl_flip[2] = 0.0;
	cont_elelamptvsv0dl_flip[3] = fCentrality;

  Double_t xyzR125_ele[3], xyzR125_pr[3], xyzR125_pi[3];
  xyzR125_ele[0] = (*elevars)[0];
  xyzR125_ele[1] = (*elevars)[1];
  xyzR125_ele[2] = (*elevars)[2];
  xyzR125_pr[0] = (*v0vars)[0];
  xyzR125_pr[1] = (*v0vars)[1];
  xyzR125_pr[2] = (*v0vars)[2];
  xyzR125_pi[0] = (*v0vars)[3];
  xyzR125_pi[1] = (*v0vars)[4];
  xyzR125_pi[2] = (*v0vars)[5];

  Double_t rdhfcutvars[9];
  rdhfcutvars[0] = xyzR125_ele[0];
  rdhfcutvars[1] = xyzR125_ele[1];
  rdhfcutvars[2] = xyzR125_ele[2];
  rdhfcutvars[3] = xyzR125_pr[0];
  rdhfcutvars[4] = xyzR125_pr[1];
  rdhfcutvars[5] = xyzR125_pr[2];
  rdhfcutvars[6] = xyzR125_pi[0];
  rdhfcutvars[7] = xyzR125_pi[1];
  rdhfcutvars[8] = xyzR125_pi[2];

  Double_t dphis_ele_pr = fAnalCuts->dPhiSR125(xyzR125_ele,xyzR125_pr);
  Double_t detas_ele_pr = fAnalCuts->dEtaSR125(xyzR125_ele,xyzR125_pr);
  Double_t dphis_ele_pi = fAnalCuts->dPhiSR125(xyzR125_ele,xyzR125_pi);
  Double_t detas_ele_pi = fAnalCuts->dEtaSR125(xyzR125_ele,xyzR125_pi);

	//if(mel < 10. && cosoa>0. &&  fAnalCuts->IsPeakRegion(v0))
	if(fAnalCuts->IsSelected(trke,v0,rdhfcutvars,AliRDHFCuts::kCandidate) &&  fAnalCuts->IsPeakRegion(v0))
	{

		if(((int)trke->T())*chargepr>0){
			fHistoEleLambdaMassRSMix->Fill(cont);
			if(trke->T()>0) fHistoEleLambdaMassRSMix1->Fill(cont);
			else fHistoEleLambdaMassRSMix2->Fill(cont);

			fHistoEleLambdaMassvsElePtRSMix->Fill(cont2);
			if(trke->T()>0) fHistoEleLambdaMassvsElePtRSMix1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtRSMix2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtRSMix->Fill(trke->Pt(),fCentrality);
				fHistoElePtvsEtaRSMix->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtRSMix->Fill(cont_eleptvslambdapt);
				fHistoElePtvsd0RSMix->Fill(cont_eleptvsd0);

				fHistoElePtvsV0dlRSMix->Fill(cont_eleptvsv0dl);
				if(trke->T()>0) fHistoElePtvsV0dlRSMix1->Fill(cont_eleptvsv0dl);
				else fHistoElePtvsV0dlRSMix2->Fill(cont_eleptvsv0dl);

				fHistoEleLambdaPtvsV0dlRSMix->Fill(cont_elelamptvsv0dl);
				if(trke->T()>0) fHistoEleLambdaPtvsV0dlRSMix1->Fill(cont_elelamptvsv0dl);
				else fHistoEleLambdaPtvsV0dlRSMix2->Fill(cont_elelamptvsv0dl);

				fHistoElePtvsV0dcaRSMix->Fill(cont_eleptvsv0dca);
				if(trke->T()>0) fHistoElePtvsV0dcaRSMix1->Fill(cont_eleptvsv0dca);
				else fHistoElePtvsV0dcaRSMix2->Fill(cont_eleptvsv0dca);
			}
      fHistodPhiSdEtaSElectronProtonR125RSMix->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125RSMix->Fill(dphis_ele_pi,detas_ele_pi);
		}else{
			fHistoEleLambdaMassWSMix->Fill(cont);
			if(trke->T()>0) fHistoEleLambdaMassWSMix1->Fill(cont);
			else fHistoEleLambdaMassWSMix2->Fill(cont);

			fHistoEleLambdaMassvsElePtWSMix->Fill(cont2);
			if(trke->T()>0) fHistoEleLambdaMassvsElePtWSMix1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtWSMix2->Fill(cont2);
			if(cont[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtWSMix->Fill(trke->Pt(),fCentrality);
				fHistoElePtvsEtaWSMix->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtWSMix->Fill(cont_eleptvslambdapt);
				fHistoElePtvsd0WSMix->Fill(cont_eleptvsd0);

				fHistoElePtvsV0dlWSMix->Fill(cont_eleptvsv0dl);
				if(trke->T()>0) fHistoElePtvsV0dlWSMix1->Fill(cont_eleptvsv0dl);
				else fHistoElePtvsV0dlWSMix2->Fill(cont_eleptvsv0dl);

				fHistoEleLambdaPtvsV0dlWSMix->Fill(cont_elelamptvsv0dl);
				if(trke->T()>0) fHistoEleLambdaPtvsV0dlWSMix1->Fill(cont_elelamptvsv0dl);
				else fHistoEleLambdaPtvsV0dlWSMix2->Fill(cont_elelamptvsv0dl);

				fHistoElePtvsV0dcaWSMix->Fill(cont_eleptvsv0dca);
				if(trke->T()>0) fHistoElePtvsV0dcaWSMix1->Fill(cont_eleptvsv0dca);
				else fHistoElePtvsV0dcaWSMix2->Fill(cont_eleptvsv0dca);
			}
      fHistodPhiSdEtaSElectronProtonR125WSMix->Fill(dphis_ele_pr,detas_ele_pr);
      fHistodPhiSdEtaSElectronPionR125WSMix->Fill(dphis_ele_pi,detas_ele_pi);
		}
	}

	//if(mel < 10. && cosoa<0. &&  fAnalCuts->IsPeakRegion(v0))
	if(mel_flip < 10. && cosoa<0. &&  fAnalCuts->IsPeakRegion(v0))
	{
		if(((int)trke->T())*chargepr>0){
			fHistoEleLambdaMassRSMixAway->Fill(cont_flip);
			if(trke->T()>0) fHistoEleLambdaMassRSMix1Away->Fill(cont_flip);
			else fHistoEleLambdaMassRSMix2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlRSMixAway->Fill(cont_elelamptvsv0dl_flip);
				if(trke->T()>0) fHistoEleLambdaPtvsV0dlRSMix1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlRSMix2Away->Fill(cont_elelamptvsv0dl_flip);
			}
		}else{
			fHistoEleLambdaMassWSMixAway->Fill(cont_flip);
			if(trke->T()>0) fHistoEleLambdaMassWSMix1Away->Fill(cont_flip);
			else fHistoEleLambdaMassWSMix2Away->Fill(cont_flip);
			if(cont_flip[0]<fAnalCuts->GetEleLambdaMassMax()){
				fHistoEleLambdaPtvsV0dlWSMixAway->Fill(cont_elelamptvsv0dl_flip);
				if(trke->T()>0) fHistoEleLambdaPtvsV0dlWSMix1Away->Fill(cont_elelamptvsv0dl_flip);
				else fHistoEleLambdaPtvsV0dlWSMix2Away->Fill(cont_elelamptvsv0dl_flip);
			}
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

  fCorrelationVariables[0] = v0->Pt();
  fCorrelationVariables[1] = trke->Pt();
  fCorrelationVariables[2] = TVector2::Phi_mpi_pi(v0->Phi()-trke->Phi());
  fCorrelationVariables[3] = v0->Eta()-trke->Eta();
  fCorrelationVariables[4] = (*v0vars)[15];
  fCorrelationVariables[5] = (*elevars)[5];
  fCorrelationVariables[6] = 2;
  if(trke->T()>0){
    if(chargepr>0) fCorrelationVariables[7] = 0;
    else fCorrelationVariables[7] = 2;
  }else if(trke->T()<0){
    if(chargepr>0) fCorrelationVariables[7] = 1;
    else fCorrelationVariables[7] = 3;
  }
  fCorrelationVariables[8] = (*elevars)[6];
  fCorrelationVariables[9] = (*elevars)[7];
  fCorrelationVariables[10] = fCentrality;
  fCorrelationVariables[11] = sqrt(pxsum*pxsum+pysum*pysum);
  fCorrelationVariables[12] = mel;
  fCorrelationVariables[13] = cosoa;

	cont_cor_nd[0] =  sqrt(pxsum*pxsum+pysum*pysum);
	cont_cor_nd[1] =  TVector2::Phi_mpi_pi(v0->Phi()-trke->Phi());
	cont_cor_nd[2] = (*v0vars)[15];
  if(trke->T()>0){
    if(chargepr>0) cont_cor_nd[3] = 0;
    else cont_cor_nd[3] = 2;
  }else if(trke->T()<0){
    if(chargepr>0) cont_cor_nd[3] = 3;
    else cont_cor_nd[3] = 1;
  }
	cont_cor_nd[4] = fCorrelationVariables[8];
	cont_cor_nd[5] = 0;
	if(fabs(fCorrelationVariables[9]-1013)<0.001) cont_cor_nd[5] = 7;
	if(fabs(fCorrelationVariables[9]-1016)<0.001) cont_cor_nd[5] = 8;
	if(fabs(fCorrelationVariables[9]-1019)<0.001) cont_cor_nd[5] = 9;
	cont_cor_nd[6] = fCentrality;

	if(fAnalCuts->IsSelected(trke,v0,rdhfcutvars,AliRDHFCuts::kCandidate) &&  fAnalCuts->IsPeakRegion(v0))
	{
		if(fWriteVariableTree)
			fCorrelationVariablesTree->Fill();

		if(mel<fAnalCuts->GetEleLambdaMassMax())
			fHistoCorrelationVariablesvsEleLambdaPtMix->Fill(cont_cor_nd);
		cont_cor_nd[0] = trke->Pt();
		fHistoCorrelationVariablesvsElePtMix->Fill(cont_cor_nd);
		cont_cor_nd[0] = v0->Pt();
		fHistoCorrelationVariablesvsLambdaPtMix->Fill(cont_cor_nd);
  }

	cont_mass_nd[0] =  mel;
	cont_mass_nd[1] =  cont_cor_nd[0];
	cont_mass_nd[4] =  cont_cor_nd[3];
	cont_mass_nd[5] =  cont_cor_nd[4];
	cont_mass_nd[6] =  cont_cor_nd[5];
	cont_mass_nd[7] =  cont_cor_nd[6];

	if(fAnalCuts->IsPeakRegion(v0)) cont_mass_nd[3] = 1;
	if(fAnalCuts->IsSideBand(v0)) cont_mass_nd[3] = 0;
	if(fAnalCuts->IsSelected(trke,v0,rdhfcutvars,AliRDHFCuts::kCandidate)) cont_mass_nd[2]=1;
	if(mel_flip < 10.&& cosoa < 0.) cont_mass_nd[2]=0;
	fHistoMassVariablesvsEleLambdaPtMix->Fill(cont_mass_nd);
	cont_mass_nd[0] = trke->Pt();
	fHistoMassVariablesvsElePtMix->Fill(cont_mass_nd);
	cont_mass_nd[0] = v0->Pt();
	fHistoMassVariablesvsLambdaPtMix->Fill(cont_mass_nd);

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineSingleTreeVariables() 
{
  //
  /// Define single tree variables
  //

  const char* nameoutput = GetOutputSlot(11)->GetContainer()->GetName();
  fSingleVariablesTree = new TTree(nameoutput,"single variables tree");
  Int_t nVar = 30;
  fCandidateSingleVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Px";
  fCandidateVariableNames[ 1]="Py";
  fCandidateVariableNames[ 2]="Pz";
  fCandidateVariableNames[ 3]="Charge";
  fCandidateVariableNames[ 4]="LambdaMass";
  fCandidateVariableNames[ 5]="Bz";
  fCandidateVariableNames[ 6]="Centrality";
  fCandidateVariableNames[ 7]="PrimVertZ";
  fCandidateVariableNames[ 8]="EvNumber";
  fCandidateVariableNames[ 9]="RunNumber";
  fCandidateVariableNames[10]="EventPlane";
  fCandidateVariableNames[11]="Vars0";//e: dca, L: proper dl
  fCandidateVariableNames[12]="Vars1";//e: trk ID, L: trk ID (pos)
  fCandidateVariableNames[13]="Vars2";//e: nSigma TPC, L: trk ID (neg)
  fCandidateVariableNames[14]="Vars3";//e: nSigma TOF L: DCA V0 to PV
  fCandidateVariableNames[15]="Vars4";//e: nSigma ITS  L: DCA V0 daughters
  fCandidateVariableNames[16]="MCVars0";// PDG 
  fCandidateVariableNames[17]="MCVars1";// Mother PDG
  fCandidateVariableNames[18]="MCVars2";// Grand Mother PDG
  fCandidateVariableNames[19]="MCVars3";// Matched Label
  fCandidateVariableNames[20]="MCVars4";// Matched PDG
  fCandidateVariableNames[21]="Vars5";//e: Golden chi2  L: Pos dca to PV
  fCandidateVariableNames[22]="Vars6";//e: nCls ITS  L: Neg dca to PV
  fCandidateVariableNames[23]="Vars7";//e: Chi2 ITS  L: V0 cosine of pointing angle
  fCandidateVariableNames[24]="Vars8";//e: nCls shared ITS  L: pT arm
  fCandidateVariableNames[25]="Vars9";//e: Not used  L: Alpha arm
  fCandidateVariableNames[26]="Vars10";//e: Not used  L: nSigmaTPC proton
  fCandidateVariableNames[27]="Vars11";//e: Not used  L: nSigmaTPC pion
  fCandidateVariableNames[28]="Vars12";//e: Not used  L: nSigmaTOF proton
  fCandidateVariableNames[29]="Vars13";//e: Not used  L: nSigmaTOF pion
  //fCandidateVariableNames[10]="PdgMother";
  //fCandidateVariableNames[14]="LabMother";
  //fCandidateVariableNames[15]="PdgGrMother";
  //fCandidateVariableNames[16]="LabGrMother";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fSingleVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateSingleVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineEleTreeVariables() 
{
  //
  /// Define electron tree variables
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillElectronROOTObjects(AliAODTrack *trk, AliAODEvent *event, TClonesArray *mcArray) 
{
  //
  /// Fill histograms or tree depending on fWriteVariableTree
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

  //Int_t ncls_its = trk->GetITSclusters(0);
  Int_t ncls_its = 0;
  for(Int_t i=0;i<6;i++){
    if(trk->HasPointOnITSLayer(i)) ncls_its++;
  }
  Int_t ncls_its_shared = 0;
  for(Int_t i=0;i<6;i++){
    if(trk->HasPointOnITSLayer(i) && trk->HasSharedPointOnITSLayer(i)) ncls_its_shared++;
  }

	if(ncls_its>0){
		fHistoBachFracSharedITS->Fill((Float_t)ncls_its_shared/(Float_t)ncls_its);
		fHistoBachChi2NclsITS->Fill(trk->GetITSchi2()/(Float_t)ncls_its);
	}else{
		fHistoBachFracSharedITS->Fill(-9999);
		fHistoBachChi2NclsITS->Fill(-9999);
	}

  Double_t d0z0[2],covd0z0[3];
  trk->PropagateToDCA(fVtx1,fBzkG,kVeryBig,d0z0,covd0z0);

  Double_t minmass_ee = 9999.;
  Double_t minmasslike_ee = 9999.;
  Bool_t isconv = fAnalCuts->TagConversions(trk,fGTIndex,event,event->GetNumberOfTracks(),minmass_ee);
  Bool_t isconv_like = fAnalCuts->TagConversionsSameSign(trk,fGTIndex,event,event->GetNumberOfTracks(),minmasslike_ee);
  fHistoMassConversionsMin->Fill(minmass_ee);
  fHistoMassConversionsSameSignMin->Fill(minmasslike_ee);

  Int_t pdgEle = -9999;
  Int_t pdgEleMother = -9999;
  Float_t mcelepx = -9999;
  Float_t mcelepy = -9999;
  Float_t mcelepz = -9999;
  Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;//for write eavh variables tree
  for(Int_t i=0;i<100;i++){
    pdgarray_ele[i]=-9999;
    labelarray_ele[i]=-9999;
  }
  ngen_ele=-9999;
  Int_t mcetype = -9999;
  if(fUseMCInfo)
  {
    Int_t labEle = trk->GetLabel();
    if(labEle<0) return;
    AliAODMCParticle *mcetrk = (AliAODMCParticle*)mcArray->At(labEle);
    if(!mcetrk) return;
    if(abs(mcetrk->GetPdgCode())!=11) return;

    pdgEle = mcetrk->GetPdgCode();
    mcelepx = mcetrk->Px();
    mcelepy = mcetrk->Py();
    mcelepz = mcetrk->Pz();

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
    Double_t xyzR125[3]={9999.,9999.,9999.};
    if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(trk,fBzkG,pv,xyzR125);
    TVector *varvec = new TVector(8);
    (*varvec)[0] = xyzR125[0];
    (*varvec)[1] = xyzR125[1];
    (*varvec)[2] = xyzR125[2];
    (*varvec)[3] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
    (*varvec)[4] = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);
    (*varvec)[5] = d0z0[0];
    (*varvec)[6] = (Int_t)isconv + 2*(Int_t)isconv_like;
    (*varvec)[7] = mcetype;

    Int_t nextRes( nextResVec[fPoolIndex] );
    m_ReservoirE[fPoolIndex][nextRes].push_back(new TLorentzVector(trk->Px(),trk->Py(),trk->Pz(),trk->Charge()));
    m_ReservoirVarsE[fPoolIndex][nextRes].push_back(varvec);
	}

	if(!fWriteEachVariableTree) return;


  //test propagation
//  cout<<"B: "<<fBzkG<<endl;
//  cout<<"trk: "<<trk->Charge()<<" "<<trk->Px()<<" "<<trk->Py()<<" "<<trk->Pz()<<endl;
//  AliExternalTrackParam etp;
//  etp.CopyFromVTrack(trk);
//  Double_t xyz_test[3]={-9999.,-9999.,-9999.};  
//  etp.GetXYZ(xyz_test);
//  cout<<"Before: "<<xyz_test[0]<<" "<<xyz_test[1]<<" "<<xyz_test[2]<<endl;
//  etp.PropagateTo(80.,(Float_t)fBzkG);
//  etp.GetXYZ(xyz_test); // GetXYZ returns global coordinates
//  cout<<"After: "<<xyz_test[0]<<" "<<xyz_test[1]<<" "<<xyz_test[2]<<endl;

	for(Int_t i=0;i<26;i++){
		fCandidateEleVariables[i] = -9999.;
	}
	for(Int_t i=0;i<30;i++){
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
  fCandidateSingleVariables[12] = trk->GetID();
  if(fAnalCuts->GetIsUsePID())
  {
		Double_t nSigmaTPCele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
		Double_t nSigmaTOFele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);
		Double_t nSigmaITSele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasITS(trk,AliPID::kElectron);
    fCandidateSingleVariables[13] = nSigmaTPCele;
    fCandidateSingleVariables[14] = nSigmaTOFele;
    fCandidateSingleVariables[15] = nSigmaITSele;
  }
  fCandidateSingleVariables[16] = pdgEle;
  fCandidateSingleVariables[17] = pdgarray_ele[0];
  fCandidateSingleVariables[18] = pdgarray_ele[1];
  fCandidateSingleVariables[19] = -9999;//mcinfo[0];
  fCandidateSingleVariables[20] = -9999;//mcinfo[1];

  fCandidateSingleVariables[21] = -9999;//trk->GetChi2TPCConstrainedVsGlobal(fVtx1);
  fCandidateSingleVariables[22] = (Float_t) ncls_its;
  fCandidateSingleVariables[23] = trk->GetITSchi2();
  fCandidateSingleVariables[24] = (Float_t)ncls_its_shared;
  fCandidateSingleVariables[25] = -9999.;//not used
  //26-30 not used

//  if(fUseMCInfo){
//    fCandidateSingleVariables[10] = pdgarray_ele[0];
//    fCandidateSingleVariables[14] = labelarray_ele[0];
//    fCandidateSingleVariables[15] = pdgarray_ele[1];
//    fCandidateSingleVariables[16] = labelarray_ele[1];
//  }

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
//  fCandidateEleVariables[20] = mcelepx;
//  fCandidateEleVariables[21] = mcelepy;
//  fCandidateEleVariables[22] = mcelepz;
//  fCandidateEleVariables[23] = fCentrality;
//  fCandidateEleVariables[24] = fVtxZ;
//  fCandidateEleVariables[25] = fRunNumber;
//
//	fHistod0Bach->Fill(d0z0[0]);
//
//	fEleVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineV0TreeVariables() 
{
  //
  /// Define V0 tree variables
  //

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fV0VariablesTree = new TTree(nameoutput,"v0 variables tree");
  Int_t nVar = 36;
  fCandidateV0Variables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="V0Px";
  fCandidateVariableNames[ 1]="V0Py";
  fCandidateVariableNames[ 2]="V0Pz";
  fCandidateVariableNames[ 3]="MassLambda";
  fCandidateVariableNames[ 4]="MassAntiLambda";
  fCandidateVariableNames[ 5]="ProtonPx";
  fCandidateVariableNames[ 6]="ProtonPy";
  fCandidateVariableNames[ 7]="ProtonPz";
  fCandidateVariableNames[ 8]="PionPx";
  fCandidateVariableNames[ 9]="PionPy";
  fCandidateVariableNames[10]="PionPz";
  fCandidateVariableNames[11]="RfidV0";
  fCandidateVariableNames[12]="DcaV0PrToPrimVertex";
  fCandidateVariableNames[13]="DcaV0PiToPrimVertex";
  fCandidateVariableNames[14]="DcaV0daughters";
  fCandidateVariableNames[15]="V0CosPointingAngle";
  fCandidateVariableNames[16]="V0ProperDecayLength";
  fCandidateVariableNames[17]="MassK0Short";
  fCandidateVariableNames[18]="nSigmaTPCpr";
  fCandidateVariableNames[19]="nSigmaTPCpi";
  fCandidateVariableNames[20]="TPCNCrossV0Pr";
  fCandidateVariableNames[21]="TPCNCrossV0Pi";
  fCandidateVariableNames[22]="TPCNCrossRatioV0Pr";
  fCandidateVariableNames[23]="TPCNCrossRatioV0Pi";
  fCandidateVariableNames[24]="V0PdgCode";
  fCandidateVariableNames[25]="V0MotherPdgCode";
  fCandidateVariableNames[26]="mcv0px";
  fCandidateVariableNames[27]="mcv0py";
  fCandidateVariableNames[28]="mcv0pz";
  fCandidateVariableNames[29]="mcv0vertx";
  fCandidateVariableNames[30]="mcv0verty";
  fCandidateVariableNames[31]="mcv0vertz";
  fCandidateVariableNames[32]="EvNumber";
  fCandidateVariableNames[33]="Centrality";
  fCandidateVariableNames[34]="PrimVertZ";
  fCandidateVariableNames[35]="RunNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fV0VariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateV0Variables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }

  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillV0ROOTObjects(AliAODv0 *v0, AliAODEvent *event, TClonesArray *mcArray) 
{
  //
  /// Fill histograms or tree depending on fWriteVariableTree
  //
	if(!v0) return;
  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
	if(!cptrack) return;
	if(!cntrack) return;
//	if(cptrack->Charge()<0 && cntrack->Charge()>0){
//		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
//		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
//	}
  //test propagation
//  cout<<"B: "<<fBzkG<<endl;
//  cout<<"trk: "<<cptrack->Charge()<<" "<<v0->MomPosX()<<" "<<v0->MomPosY()<<" "<<v0->MomPosZ()<<endl;
//  cout<<"vert: "<<v0->DecayVertexV0X()<<" "<<v0->DecayVertexV0Y()<<" "<<v0->DecayVertexV0Z()<<endl;
//  AliExternalTrackParam etp;
//  etp.CopyFromVTrack(cptrack);
//  Double_t xyz_test[3]={-9999.,-9999.,-9999.};  
//  etp.PropagateTo(10.,(Float_t)fBzkG);
//  etp.GetXYZ(xyz_test); // GetXYZ returns global coordinates
//  cout<<"10: "<<xyz_test[0]<<" "<<xyz_test[1]<<" "<<xyz_test[2]<<endl;
//  etp.PropagateTo(120.,(Float_t)fBzkG);
//  etp.GetXYZ(xyz_test); // GetXYZ returns global coordinates
//  cout<<"120: "<<xyz_test[0]<<" "<<xyz_test[1]<<" "<<xyz_test[2]<<endl;

  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t posVtx[3] = {0.,0.,0.};
  fVtx1->GetXYZ(posVtx);
  Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
  Double_t v0propdl = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
  Double_t cont_arm[3];
  cont_arm[0] = v0->Pt();
  cont_arm[1] = v0->PtArmV0();
  cont_arm[2] = v0->AlphaV0();

	if(TMath::Abs(v0->MassLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
		fHistoLambdaMassvsPt->Fill(v0->MassLambda(),v0->Pt());
	}
	if(TMath::Abs(v0->MassAntiLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
		fHistoLambdaMassvsPt->Fill(v0->MassAntiLambda(),v0->Pt());
	}
	fHistoK0sMassvsPt->Fill(v0->MassK0Short(),v0->Pt());

	if(fAnalCuts->IsPeakRegion(v0)){
		fHistoLambdaPtvsDl->Fill(v0->Pt(),v0propdl);
		fHistoLambdaPtvsEta->Fill(v0->Pt(),v0->Eta());
		fHistoLambdaArmenteros->Fill(cont_arm);
	}
	if(fAnalCuts->IsSideBand(v0)){
		fHistoLambdaPtvsDlSide->Fill(v0->Pt(),v0propdl);
	}

	Double_t momv0x = v0->MomV0X();
	Double_t momv0y = v0->MomV0Y();
	Double_t phi_alice = atan2(momv0y,momv0x);
	if(phi_alice<0.) phi_alice += 2 * M_PI;
	fHistoLambdaQovPtvsPhi->Fill(phi_alice,1./sqrt(momv0x*momv0x+momv0y*momv0y));

  if(fDoEventMixing){
    Int_t nextRes( nextResVec[fPoolIndex] );
    Double_t posVtx[3] = {0.,0.,0.};
    fVtx1->GetXYZ(posVtx);
    TLorentzVector *lv = new TLorentzVector();
    Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
    Double_t xyzR125pr[3] = {9999.,9999.,9999.};
    Double_t xyzR125pi[3] = {9999.,9999.,9999.};
    if(TMath::Abs(v0->MassLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
      lv->SetXYZM(v0->Px(),v0->Py(),v0->Pz(),v0->MassLambda());
      m_ReservoirL1[fPoolIndex][nextRes].push_back(lv);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(cptrack,fBzkG,posVtx,xyzR125pr);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(cntrack,fBzkG,posVtx,xyzR125pi);
      TVector *varvec = new TVector(17);
      (*varvec)[0] = xyzR125pr[0];
      (*varvec)[1] = xyzR125pr[1];
      (*varvec)[2] = xyzR125pr[2];
      (*varvec)[3] = xyzR125pi[0];
      (*varvec)[4] = xyzR125pi[1];
      (*varvec)[5] = xyzR125pi[2];
      (*varvec)[6] = v0->DecayVertexV0X();
      (*varvec)[7] = v0->DecayVertexV0Y();
      (*varvec)[8] = v0->DecayVertexV0Z();
      (*varvec)[9] = v0->MomPosX();
      (*varvec)[10] = v0->MomPosY();
      (*varvec)[11] = v0->MomPosZ();
      (*varvec)[12] = v0->MomNegX();
      (*varvec)[13] = v0->MomNegY();
      (*varvec)[14] = v0->MomNegZ();
      (*varvec)[15] = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
      (*varvec)[16] = v0->DcaV0ToPrimVertex();
      m_ReservoirVarsL1[fPoolIndex][nextRes].push_back(varvec);
    }else{
      lv->SetXYZM(v0->Px(),v0->Py(),v0->Pz(),v0->MassAntiLambda());
      m_ReservoirL2[fPoolIndex][nextRes].push_back(lv);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(cntrack,fBzkG,posVtx,xyzR125pr);
      if(fAnalCuts->GetCuts()[2]>0. || fAnalCuts->GetCuts()[3]>0.) fAnalCuts->SetSftPosR125(cptrack,fBzkG,posVtx,xyzR125pi);
      TVector *varvec = new TVector(17);
      (*varvec)[0] = xyzR125pr[0];
      (*varvec)[1] = xyzR125pr[1];
      (*varvec)[2] = xyzR125pr[2];
      (*varvec)[3] = xyzR125pi[0];
      (*varvec)[4] = xyzR125pi[1];
      (*varvec)[5] = xyzR125pi[2];
      (*varvec)[6] = v0->DecayVertexV0X();
      (*varvec)[7] = v0->DecayVertexV0Y();
      (*varvec)[8] = v0->DecayVertexV0Z();
      (*varvec)[9] = v0->MomNegX();
      (*varvec)[10] = v0->MomNegY();
      (*varvec)[11] = v0->MomNegZ();
      (*varvec)[12] = v0->MomPosX();
      (*varvec)[13] = v0->MomPosY();
      (*varvec)[14] = v0->MomPosZ();
      (*varvec)[15] = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
      (*varvec)[16] = v0->DcaV0ToPrimVertex();
      m_ReservoirVarsL2[fPoolIndex][nextRes].push_back(varvec);
    }
  }

	if(!fWriteEachVariableTree) return;

	Int_t v0pdgcode = -9999;
	Int_t v0motherpdgcode = -9999;
	Float_t mcv0px = -9999;
	Float_t mcv0py = -9999;
	Float_t mcv0pz = -9999;
	Float_t mcv0vertx = -9999;
	Float_t mcv0verty = -9999;
	Float_t mcv0vertz = -9999;
  Int_t pdgarray_v0[100], labelarray_v0[100], ngen_v0;//for write eavh variables tree
	if(fUseMCInfo)
	{
		Int_t pdgdgv0[2]={2212,211};
		Int_t labV0 = v0->MatchToMC(3122,mcArray,2,pdgdgv0); // the V0
		if(labV0>=0){
			if(TMath::Abs(v0->MassLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
				fHistoLambdaMassvsPtMCS->Fill(v0->MassLambda(),v0->Pt());
			}
			if(TMath::Abs(v0->MassAntiLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
				fHistoLambdaMassvsPtMCS->Fill(v0->MassAntiLambda(),v0->Pt());
			}
      if(fAnalCuts->IsPeakRegion(v0)){
        fHistoLambdaPtvsEtaMCS->Fill(v0->Pt(),v0->Eta());
        fHistoLambdaArmenterosMCS->Fill(cont_arm);
      }

			AliAODMCParticle *mcv0trk = (AliAODMCParticle*)mcArray->At(labV0);
			if(mcv0trk){
				mcv0px = mcv0trk->Px();
				mcv0py = mcv0trk->Py();
				mcv0pz = mcv0trk->Pz();

				Int_t labptrk = cptrack->GetLabel();
				if(labptrk>=0){
					AliAODMCParticle *mcv0ptrk = (AliAODMCParticle*)mcArray->At(labptrk);
					if(mcv0ptrk){
						mcv0vertx = mcv0ptrk->Xv();
						mcv0verty = mcv0ptrk->Yv();
						mcv0vertz = mcv0ptrk->Zv();
					}
				}

				GetMCDecayHistory(mcv0trk,mcArray,pdgarray_v0,labelarray_v0,ngen_v0);

				Bool_t hfv0_flag = kFALSE;
				Bool_t promptv0_flag = kFALSE;
				v0pdgcode = mcv0trk->GetPdgCode();
				Int_t labv0mother = mcv0trk->GetMother();
				if(labv0mother>=0){
					AliAODMCParticle *motherv0 = (AliAODMCParticle*)mcArray->At(labv0mother);
					if(motherv0){
						v0motherpdgcode = motherv0->GetPdgCode();
						if(abs(v0motherpdgcode)>4000&&abs(v0motherpdgcode)<4400){
							hfv0_flag = kTRUE;
						}
						if(abs(v0motherpdgcode)==3322){
							fHistoLambdaPtvsDlFeeddownXi0MCS->Fill(v0->Pt(),v0propdl);
						}else if(abs(v0motherpdgcode)==3312){
							fHistoLambdaPtvsDlFeeddownXiMinusMCS->Fill(v0->Pt(),v0propdl);
						}else if(abs(v0motherpdgcode)==3334){
							fHistoLambdaPtvsDlFeeddownOmegaMCS->Fill(v0->Pt(),v0propdl);
						}else{
							fHistoLambdaPtvsDlMCS->Fill(v0->Pt(),v0propdl);
							promptv0_flag = kTRUE;
						}
					}
					//if(!hfv0_flag) return;
					//if(!promptv0_flag) return;
				}

			}else{
				for(Int_t i=0;i<100;i++){
					pdgarray_v0[i] = -9999;
					labelarray_v0[i] = -9999;
				}
				ngen_v0 = -9999;
			}
		}
	}

	for(Int_t i=0;i<36;i++){
		fCandidateV0Variables[i] = -9999.;
	}

	for(Int_t i=0;i<30;i++){
		fCandidateSingleVariables[i] = -9999.;
	}

  Bool_t  isparticle = kFALSE;
  if(TMath::Abs(v0->MassLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()){
    isparticle = kTRUE;
  }

  fCandidateSingleVariables[ 0] = v0->Px();
  fCandidateSingleVariables[ 1] = v0->Py();
  fCandidateSingleVariables[ 2] = v0->Pz();
  if(isparticle){
    fCandidateSingleVariables[ 3] = 1.;
    fCandidateSingleVariables[ 4] = v0->MassLambda();
  }else{
    fCandidateSingleVariables[ 3] = -1.;
    fCandidateSingleVariables[ 4] = v0->MassAntiLambda();
  }
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
  fCandidateSingleVariables[11] = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
  fCandidateSingleVariables[12] = cptrack->GetID();
  fCandidateSingleVariables[13] = cntrack->GetID();
  fCandidateSingleVariables[14] = v0->DcaV0ToPrimVertex();
  fCandidateSingleVariables[15] = v0->DcaV0Daughters();

  fCandidateSingleVariables[16] = v0pdgcode;
  fCandidateSingleVariables[17] = pdgarray_v0[0];
  fCandidateSingleVariables[18] = pdgarray_v0[1];
  fCandidateSingleVariables[19] = -9999;//mcinfo[0];
  fCandidateSingleVariables[20] = -9999;//mcinfo[1];

  fCandidateSingleVariables[21] = v0->DcaPosToPrimVertex();
  fCandidateSingleVariables[22] = v0->DcaNegToPrimVertex();
  fCandidateSingleVariables[23] = v0->CosPointingAngle(posVtx);
//  Float_t fracsharedpostpc = -9999;
//  if(cptrack->GetTPCNcls()>0) fracsharedpostpc = (Float_t)cptrack->GetTPCnclsS()/(Float_t)cptrack->GetTPCncls();
//  Float_t fracsharednegtpc = -9999;
//  if(cntrack->GetTPCNcls()>0) fracsharednegtpc = (Float_t)cntrack->GetTPCnclsS()/(Float_t)cntrack->GetTPCncls();
//  fCandidateSingleVariables[24] = fracsharedpostpc;
//  fCandidateSingleVariables[25] = fracsharednegtpc;

  fCandidateSingleVariables[24] = v0->PtArmV0();
  fCandidateSingleVariables[25] = v0->AlphaV0();

  if(fAnalCuts->GetUseLambdaPID())
  {
    if(isparticle){
      fCandidateSingleVariables[26] = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
      fCandidateSingleVariables[27] = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
      fCandidateSingleVariables[28] = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
      fCandidateSingleVariables[29] = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
    }else{
      fCandidateSingleVariables[26] = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
      fCandidateSingleVariables[27] = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
      fCandidateSingleVariables[28] = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
      fCandidateSingleVariables[29] = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
    }
  }

//  if(fUseMCInfo){
//    fCandidateSingleVariables[10] = pdgarray_v0[0];
//    fCandidateSingleVariables[14] = labelarray_v0[0];
//    fCandidateSingleVariables[15] = pdgarray_v0[1];
//    fCandidateSingleVariables[16] = labelarray_v0[1];
//  }

	fSingleVariablesTree->Fill();
//  fCandidateV0Variables[ 0] = v0->Px();
//  fCandidateV0Variables[ 1] = v0->Py();
//  fCandidateV0Variables[ 2] = v0->Pz();
//  fCandidateV0Variables[ 3] = v0->MassLambda();
//  fCandidateV0Variables[ 4] = v0->MassAntiLambda();
//
//	Bool_t isparticle = kTRUE;
//	if(fabs(v0->MassAntiLambda()-mlamPDG)<fAnalCuts->GetProdV0MassTolLambdaRough()) isparticle=kFALSE;
//
//	if(isparticle){
//		fCandidateV0Variables[ 5] = v0->MomPosX(); //cptrack->Px();
//		fCandidateV0Variables[ 6] = v0->MomPosY(); //cptrack->Py();
//		fCandidateV0Variables[ 7] = v0->MomPosZ(); //cptrack->Pz();
//		fCandidateV0Variables[ 8] = v0->MomNegX(); //cntrack->Px();
//		fCandidateV0Variables[ 9] = v0->MomNegY(); //cntrack->Py();
//		fCandidateV0Variables[10] = v0->MomNegZ(); //cntrack->Pz();
//	}else{
//		fCandidateV0Variables[ 5] =  v0->MomNegX(); //cntrack->Px();
//		fCandidateV0Variables[ 6] =  v0->MomNegY(); //cntrack->Py();
//		fCandidateV0Variables[ 7] =  v0->MomNegZ(); //cntrack->Pz();
//		fCandidateV0Variables[ 8] =  v0->MomPosX(); //cptrack->Px();
//		fCandidateV0Variables[ 9] =  v0->MomPosY(); //cptrack->Py();
//		fCandidateV0Variables[10] =  v0->MomPosZ(); //cptrack->Pz();
//	}
//
//  Double_t lPosV0[3];
//  lPosV0[0] = v0->DecayVertexV0X();
//  lPosV0[1] = v0->DecayVertexV0Y();
//  lPosV0[2] = v0->DecayVertexV0Z();
//  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
//	fCandidateV0Variables[11] = decayvertV0;
//
//	Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
//	Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
//  if(isparticle){
//		fCandidateV0Variables[12] = lDcaPosToPrimVertex;
//		fCandidateV0Variables[13] = lDcaNegToPrimVertex;
//  }else{
//		fCandidateV0Variables[12] = lDcaNegToPrimVertex;
//		fCandidateV0Variables[13] = lDcaPosToPrimVertex;
//  }
//	fCandidateV0Variables[14] = v0->DcaV0Daughters();
//  fCandidateV0Variables[15] = v0->CosPointingAngle(posVtx); 
//  fCandidateV0Variables[16] = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
//  fCandidateV0Variables[17] = v0->MassK0Short();
//
//  if(fAnalCuts->GetUseLambdaPID())
//  {
//		if(isparticle){
//			Double_t nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
//			Double_t nSigmaTPCv0pi = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
//			fCandidateV0Variables[18] = nSigmaTPCv0pr;
//			fCandidateV0Variables[19] = nSigmaTPCv0pi;
//		}else{
//			Double_t nSigmaTPCv0pr = fAnalCuts->GetPidProton()->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
//			Double_t nSigmaTPCv0pi = fAnalCuts->GetPidPion()->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
//			fCandidateV0Variables[18] = nSigmaTPCv0pr;
//			fCandidateV0Variables[19] = nSigmaTPCv0pi;
//		}
//  }
//	if(isparticle){
//		fCandidateV0Variables[20] = cptrack->GetTPCClusterInfo(2,1);
//		fCandidateV0Variables[21] = cntrack->GetTPCClusterInfo(2,1);
//		if(cptrack->GetTPCNclsF()>0)
//			fCandidateV0Variables[22] = (Float_t) cptrack->GetTPCClusterInfo(2,1)/(Float_t)cptrack->GetTPCNclsF();
//		if(cntrack->GetTPCNclsF()>0)
//			fCandidateV0Variables[23] =(Float_t)  cntrack->GetTPCClusterInfo(2,1)/(Float_t)cntrack->GetTPCNclsF();
//	}else{
//		fCandidateV0Variables[20] = cntrack->GetTPCClusterInfo(2,1);
//		fCandidateV0Variables[21] = cptrack->GetTPCClusterInfo(2,1);
//		if(cntrack->GetTPCNclsF()>0)
//			fCandidateV0Variables[22] = (Float_t) cntrack->GetTPCClusterInfo(2,1)/(Float_t)cntrack->GetTPCNclsF();
//		if(cptrack->GetTPCNclsF()>0)
//			fCandidateV0Variables[23] = (Float_t) cptrack->GetTPCClusterInfo(2,1)/(Float_t)cptrack->GetTPCNclsF();
//	}
//	fCandidateV0Variables[24] = v0pdgcode;
//	fCandidateV0Variables[25] = v0motherpdgcode;
//	fCandidateV0Variables[26] = mcv0px;
//	fCandidateV0Variables[27] = mcv0py;
//	fCandidateV0Variables[28] = mcv0pz;
//	fCandidateV0Variables[29] = mcv0vertx;
//	fCandidateV0Variables[30] = mcv0verty;
//	fCandidateV0Variables[31] = mcv0vertz;
//	fCandidateV0Variables[32] = fEvNumberCounter;
//	fCandidateV0Variables[33] = fCentrality;
//	fCandidateV0Variables[34] = fVtxZ;
//	fCandidateV0Variables[35] = fRunNumber;
//
//
//		fV0VariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineMCTreeVariables() 
{
  ///
  /// Define electron tree variables
  ///

  const char* nameoutput = GetOutputSlot(7)->GetContainer()->GetName();
  fMCVariablesTree = new TTree(nameoutput,"MC variables tree");
  Int_t nVar = 11;
  fCandidateMCVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="DecayType";
  fCandidateVariableNames[ 2]="LcPx";
  fCandidateVariableNames[ 3]="LcPy";
  fCandidateVariableNames[ 4]="LcPz";
  fCandidateVariableNames[ 5]="ElePx";
  fCandidateVariableNames[ 6]="ElePy";
  fCandidateVariableNames[ 7]="ElePz";
  fCandidateVariableNames[ 8]="V0Px";
  fCandidateVariableNames[ 9]="V0Py";
  fCandidateVariableNames[10]="V0Pz";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMCROOTObjects(AliAODMCParticle *mcpart, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype) 
{
  //
  /// Fill histograms or tree depending on fWriteMCVariableTree
  //
	if(!mcpart) return;
	if(!mcepart) return;
	if(!mcv0part) return;

	for(Int_t i=0;i<11;i++){
		fCandidateMCVariables[i] = -9999.;
	}

	fCandidateMCVariables[ 0] = fCentrality;
	fCandidateMCVariables[ 1] = (Float_t) decaytype;
	fCandidateMCVariables[ 2] = mcpart->Px();
	fCandidateMCVariables[ 3] = mcpart->Py();
	fCandidateMCVariables[ 4] = mcpart->Pz();
	fCandidateMCVariables[ 5] = mcepart->Px();
	fCandidateMCVariables[ 6] = mcepart->Py();
	fCandidateMCVariables[ 7] = mcepart->Pz();
	fCandidateMCVariables[ 8] = mcv0part->Px();
	fCandidateMCVariables[ 9] = mcv0part->Py();
	fCandidateMCVariables[10] = mcv0part->Pz();

	Double_t epx = mcepart->Px();
	Double_t epy = mcepart->Py();
	Double_t epz = mcepart->Pz();
	Double_t eE = sqrt(epx*epx+epy*epy+epz*epz+0.000511*0.000511);
	Double_t v0px = mcv0part->Px();
	Double_t v0py = mcv0part->Py();
	Double_t v0pz = mcv0part->Pz();
	Double_t v0E = sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz+1.1157*1.1157);

	Double_t InvMassEleLambda = sqrt(pow(eE+v0E,2)-pow(epx+v0px,2)-pow(epy+v0py,2)-pow(epz+v0pz,2));

	Double_t cont[3];
	cont[0] = InvMassEleLambda;
	cont[1] = mcpart->Pt();
	cont[2] = fCentrality;
	Double_t cont2[3];
	cont2[0] = InvMassEleLambda;
	cont2[1] = mcepart->Pt();
	cont2[2] = fCentrality;
	Double_t cont_eleptvseta[3];
	cont_eleptvseta[0] = mcepart->Pt();
	cont_eleptvseta[1] = mcepart->Eta();
	cont_eleptvseta[2] = fCentrality;
	Double_t cont_eleptvslambdapt[3];
	cont_eleptvslambdapt[0] = mcepart->Pt();
	cont_eleptvslambdapt[1] = mcv0part->Pt();
	cont_eleptvslambdapt[2] = fCentrality;
	Double_t cont_eleptvslambdaptvslcpt[4];
	cont_eleptvslambdaptvslcpt[0] = mcepart->Pt();
	cont_eleptvslambdaptvslcpt[1] = mcv0part->Pt();
	cont_eleptvslambdaptvslcpt[2] = mcpart->Pt();
	cont_eleptvslambdaptvslcpt[3] = fCentrality;
	Double_t contmc[3];
	contmc[0] = mcpart->Pt();
	contmc[1] = mcpart->Y();
	contmc[2] = fCentrality;
	Double_t contmcele[3];
	contmcele[0] = mcepart->Pt();
	contmcele[1] = mcepart->Eta();
	contmcele[2] = fCentrality;
	Double_t contmc_withele[3];
	contmc_withele[0] = mcpart->Pt();
	contmc_withele[1] = mcpart->Y();
	contmc_withele[2] = mcepart->Y();

	AliESDtrackCuts *esdcuts = fAnalCuts->GetTrackCuts();
	Float_t etamin, etamax;
	esdcuts->GetEtaRange(etamin,etamax);

	if(decaytype==0){
		fHistoLcMCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoLcMCGen1->Fill(contmc);
		if(mcpart->GetPdgCode()<0) fHistoLcMCGen2->Fill(contmc);
		fHistoLcMCGenWithEleA->Fill(contmc_withele);
    if(fabs(mcepart->Eta())<0.8){
      fHistoLcMCGenWithEleB->Fill(mcpart->Pt(),mcepart->Pt());
    }
		fHistoLcElectronMCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoLcElectronMCGen1->Fill(contmcele);
		if(mcepart->GetPdgCode()>0) fHistoLcElectronMCGen2->Fill(contmcele);
		fHistoEleLambdaMassMCGen->Fill(cont);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtMCGen->Fill(cont2);
			if(mcepart->GetPdgCode()<0) fHistoEleLambdaMassvsElePtMCGen1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtMCGen2->Fill(cont2);
			if(InvMassEleLambda<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtMCGen->Fill(cont_eleptvslambdapt);
			}
		}
		if(fabs(mcpart->Y())<0.7){
			if(InvMassEleLambda<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtvsLambdaPtMCLcGen->Fill(cont_eleptvslambdapt);
				fHistoElePtvsLambdaPtvsLcPtMCGen->Fill(cont_eleptvslambdaptvslcpt);
			}
		}
	}else if(decaytype==1){
		fHistoFeedDownXic0MCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoFeedDownXic0MCGen1->Fill(contmc);
		if(mcpart->GetPdgCode()<0) fHistoFeedDownXic0MCGen2->Fill(contmc);
		fHistoEleLambdaMassFeeddownXic0MCGen->Fill(cont);
		fHistoElectronFeedDownXic0MCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoElectronFeedDownXic0MCGen1->Fill(contmcele);
		if(mcepart->GetPdgCode()>0) fHistoElectronFeedDownXic0MCGen2->Fill(contmcele);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtFeeddownXic0MCGen->Fill(cont2);
			if(mcepart->GetPdgCode()<0) fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2->Fill(cont2);
			if(InvMassEleLambda<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtFeeddownXic0MCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaFeeddownXic0MCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtFeeddownXic0MCGen->Fill(cont_eleptvslambdapt);
			}
		}
	}else if(decaytype==2){
		fHistoFeedDownXicPlusMCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoFeedDownXicPlusMCGen1->Fill(contmc);
		if(mcpart->GetPdgCode()<0) fHistoFeedDownXicPlusMCGen2->Fill(contmc);
		fHistoEleLambdaMassFeeddownXicPlusMCGen->Fill(cont);
		fHistoElectronFeedDownXicPlusMCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoElectronFeedDownXicPlusMCGen1->Fill(contmcele);
		if(mcepart->GetPdgCode()>0) fHistoElectronFeedDownXicPlusMCGen2->Fill(contmcele);
		if(fabs(mcepart->Eta())<etamax){
			fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen->Fill(cont2);
			if(mcepart->GetPdgCode()<0) fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1->Fill(cont2);
			else fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2->Fill(cont2);
			if(InvMassEleLambda<fAnalCuts->GetEleLambdaMassMax()){
				fHistoElePtFeeddownXicPlusMCGen->Fill(mcepart->Pt(),fCentrality);
				fHistoElePtvsEtaFeeddownXicPlusMCGen->Fill(cont_eleptvseta);
				fHistoElePtvsLambdaPtFeeddownXicPlusMCGen->Fill(cont_eleptvslambdapt);
			}
		}
	}else if(decaytype==3){
		fHistoFeedDownSigmaMCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoFeedDownSigmaMCGen1->Fill(contmc);
		if(mcpart->GetPdgCode()<0) fHistoFeedDownSigmaMCGen2->Fill(contmc);
	}else if(decaytype==10){
		fHistoLbElectronMCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoLbElectronMCGen1->Fill(contmcele);
		if(mcepart->GetPdgCode()>0) fHistoLbElectronMCGen2->Fill(contmcele);
		fHistoLbMCGen->Fill(contmc);
		if(mcpart->GetPdgCode()>0) fHistoLbMCGen1->Fill(contmc);
		if(mcpart->GetPdgCode()<0) fHistoLbMCGen2->Fill(contmc);

		fHistoLcMCGenWithEleA->Fill(contmc_withele);
    if(fabs(mcepart->Eta())<0.8){
      fHistoLcMCGenWithEleB->Fill(mcpart->Pt(),mcepart->Pt());
    }
	}else if(decaytype==11){
		fHistoElectronFeedDownXibMCGen->Fill(contmcele);
		if(mcepart->GetPdgCode()<0) fHistoElectronFeedDownXibMCGen1->Fill(contmcele);
		if(mcepart->GetPdgCode()>0) fHistoElectronFeedDownXibMCGen2->Fill(contmcele);
	}

	if(fWriteMCVariableTree)
		fMCVariablesTree->Fill();
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineMCEleTreeVariables() 
{
  //
  // Define electron tree variables
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMCEleROOTObjects(AliAODMCParticle *mcepart, TClonesArray *mcArray) 
{
  //
  // Fill tree depending on fWriteMCVariableTree 
  //
	if(!mcepart) return;

  Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;
  GetMCDecayHistory(mcepart,mcArray,pdgarray_ele,labelarray_ele,ngen_ele);
  Bool_t ele_from_bottom = HaveBottomInHistory(pdgarray_ele);
  Bool_t ele_from_charm = HaveCharmInHistory(pdgarray_ele);
  //Int_t semi_flag = FromSemileptonicDecays(pdgarray_ele);
  Int_t semi_flag = 0;
	if(pdgarray_ele[0]>400&&pdgarray_ele[0]<450) semi_flag = 1;
	if(pdgarray_ele[0]>4000&&pdgarray_ele[0]<4500) semi_flag = 1;
	if(pdgarray_ele[0]>500&&pdgarray_ele[0]<550) semi_flag = 2;
	if(pdgarray_ele[0]>5000&&pdgarray_ele[0]<5500) semi_flag = 2;

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

//	if(fWriteMCVariableTree && fWriteEachVariableTree && mcepart->Pt()>0.4 && fabs(mcepart->Eta())<1.0 )
//		fMCEleVariablesTree->Fill();

}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineMCV0TreeVariables() 
{
  //
  // Define Mc v0 tree variables
  //

  const char* nameoutput = GetOutputSlot(10)->GetContainer()->GetName();
  fMCV0VariablesTree = new TTree(nameoutput,"MC v0 variables tree");
  Int_t nVar = 8;
  fCandidateMCV0Variables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[ 0]="Centrality";
  fCandidateVariableNames[ 1]="V0Px";
  fCandidateVariableNames[ 2]="V0Py";
  fCandidateVariableNames[ 3]="V0Pz";
  fCandidateVariableNames[ 4]="V0PdgCode";
  fCandidateVariableNames[ 5]="V0MotherPdgCode";
  fCandidateVariableNames[ 6]="RunNumber";
  fCandidateVariableNames[ 7]="EvNumber";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCV0VariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCV0Variables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}
////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMCV0ROOTObjects(AliAODMCParticle *mcv0part, TClonesArray *mcArray) 
{
  //
  // Fill histograms or tree depending on fWriteMCVariableTree 
  //
	if(!mcv0part) return;

	for(Int_t i=0;i<8;i++){
		fCandidateMCV0Variables[i] = -9999.;
	}

	Bool_t hfv0_flag = kFALSE;
	Int_t labv0mother = mcv0part->GetMother();
	Int_t pdgmotherv0 = -9999;
	if(labv0mother>=0){
		AliAODMCParticle *motherv0 = (AliAODMCParticle*)mcArray->At(labv0mother);
		if(motherv0){
			pdgmotherv0 = motherv0->GetPdgCode();
			if(abs(pdgmotherv0)>4000&&abs(pdgmotherv0)<4400){
				hfv0_flag = kTRUE;
			}
		}
	}
	if(!hfv0_flag) return;

	Double_t contmc[3];
	contmc[0] = mcv0part->Pt();
	contmc[1] = mcv0part->Eta();
	contmc[2] = fCentrality;
	fHistoLambdaMCGen->Fill(contmc);

	fCandidateMCV0Variables[ 0] = fCentrality;
	fCandidateMCV0Variables[ 1] = mcv0part->Px();
	fCandidateMCV0Variables[ 2] = mcv0part->Py();
	fCandidateMCV0Variables[ 3] = mcv0part->Pz();
	fCandidateMCV0Variables[ 4] = mcv0part->GetPdgCode();
	fCandidateMCV0Variables[ 5] = pdgmotherv0;
	fCandidateMCV0Variables[ 6] = fRunNumber;
	fCandidateMCV0Variables[ 7] = fEvNumberCounter;

	if(fWriteMCVariableTree && fWriteEachVariableTree && mcv0part->Pt()>0.4 && fabs(mcv0part->Eta())<1.0 )
		fMCV0VariablesTree->Fill();
}


////__________________________________________________________________________
void  AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineGeneralHistograms() {
  //
  /// This is to define general histograms
  //

  fCEvents = new TH1F("fCEvents","conter",20,-0.5,19.5);
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
  fCEvents->GetXaxis()->SetBinLabel(19,"No RP");
  fCEvents->GetXaxis()->SetBinLabel(20,"No Qvector (only for TPC event plane)");
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
	fHistoPoolMonitor =  new TH2F("fHistoPoolMonitor","Pool ID vs depth",4000,-0.5,3999.5,50,-0.5,49.5);
	fHistoPoolIDCounter =  new TH1F("fHistoPoolIDCounter","Pool ID",4003,-3.5,3999.5);

  fOutput->Add(fCEvents);
  fOutput->Add(fHTrigger);
  fOutput->Add(fHCentrality);
  fOutput->Add(fHEventPlane);
  fOutput->Add(fHNTrackletvsZ);
  fOutput->Add(fHNTrackletCorrvsZ);
	fOutput->Add(fHistoPoolMonitor);
	fOutput->Add(fHistoPoolIDCounter);

  return;
}
//__________________________________________________________________________
void  AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineAnalysisHistograms() 
{
  //
  /// Define analyis histograms
  //
	
  //------------------------------------------------
  // Basic histogram
  //------------------------------------------------
  Int_t bins_base[3]=		{23,100		,10};
  Double_t xmin_base[3]={1.1,0		,0.00};
  Double_t xmax_base[3]={5.7,20.	,100};
  fHistoEleLambdaMass = new THnSparseF("fHistoEleLambdaMass","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMass);
  fHistoEleLambdaMassRS = new THnSparseF("fHistoEleLambdaMassRS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS);
  fHistoEleLambdaMassRS1 = new THnSparseF("fHistoEleLambdaMassRS1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS1);
  fHistoEleLambdaMassRS2 = new THnSparseF("fHistoEleLambdaMassRS2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS2);
  fHistoEleLambdaMassWS = new THnSparseF("fHistoEleLambdaMassWS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS);
  fHistoEleLambdaMassWS1 = new THnSparseF("fHistoEleLambdaMassWS1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS1);
  fHistoEleLambdaMassWS2 = new THnSparseF("fHistoEleLambdaMassWS2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS2);
  fHistoEleLambdaMassRSMix = new THnSparseF("fHistoEleLambdaMassRSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix);
  fHistoEleLambdaMassRSMix1 = new THnSparseF("fHistoEleLambdaMassRSMix1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix1);
  fHistoEleLambdaMassRSMix2 = new THnSparseF("fHistoEleLambdaMassRSMix2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix2);
  fHistoEleLambdaMassWSMix = new THnSparseF("fHistoEleLambdaMassWSMix","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix);
  fHistoEleLambdaMassWSMix1 = new THnSparseF("fHistoEleLambdaMassWSMix1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix1);
  fHistoEleLambdaMassWSMix2 = new THnSparseF("fHistoEleLambdaMassWSMix2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix2);
  fHistoEleLambdaMassRSSide = new THnSparseF("fHistoEleLambdaMassRSSide","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSide);
  fHistoEleLambdaMassRSSide1 = new THnSparseF("fHistoEleLambdaMassRSSide1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSide1);
  fHistoEleLambdaMassRSSide2 = new THnSparseF("fHistoEleLambdaMassRSSide2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSide2);
  fHistoEleLambdaMassWSSide = new THnSparseF("fHistoEleLambdaMassWSSide","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSide);
  fHistoEleLambdaMassWSSide1 = new THnSparseF("fHistoEleLambdaMassWSSide1","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSide1);
  fHistoEleLambdaMassWSSide2 = new THnSparseF("fHistoEleLambdaMassWSSide2","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSide2);

  fHistoEleLambdaMassRSAway = new THnSparseF("fHistoEleLambdaMassRSAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSAway);
  fHistoEleLambdaMassRS1Away = new THnSparseF("fHistoEleLambdaMassRS1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS1Away);
  fHistoEleLambdaMassRS2Away = new THnSparseF("fHistoEleLambdaMassRS2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRS2Away);
  fHistoEleLambdaMassWSAway = new THnSparseF("fHistoEleLambdaMassWSAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSAway);
  fHistoEleLambdaMassWS1Away = new THnSparseF("fHistoEleLambdaMassWS1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS1Away);
  fHistoEleLambdaMassWS2Away = new THnSparseF("fHistoEleLambdaMassWS2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWS2Away);
  fHistoEleLambdaMassRSMixAway = new THnSparseF("fHistoEleLambdaMassRSMixAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMixAway);
  fHistoEleLambdaMassRSMix1Away = new THnSparseF("fHistoEleLambdaMassRSMix1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix1Away);
  fHistoEleLambdaMassRSMix2Away = new THnSparseF("fHistoEleLambdaMassRSMix2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSMix2Away);
  fHistoEleLambdaMassWSMixAway = new THnSparseF("fHistoEleLambdaMassWSMixAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMixAway);
  fHistoEleLambdaMassWSMix1Away = new THnSparseF("fHistoEleLambdaMassWSMix1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix1Away);
  fHistoEleLambdaMassWSMix2Away = new THnSparseF("fHistoEleLambdaMassWSMix2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSMix2Away);
  fHistoEleLambdaMassRSSideAway = new THnSparseF("fHistoEleLambdaMassRSSideAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSideAway);
  fHistoEleLambdaMassRSSide1Away = new THnSparseF("fHistoEleLambdaMassRSSide1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSide1Away);
  fHistoEleLambdaMassRSSide2Away = new THnSparseF("fHistoEleLambdaMassRSSide2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassRSSide2Away);
  fHistoEleLambdaMassWSSideAway = new THnSparseF("fHistoEleLambdaMassWSSideAway","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSideAway);
  fHistoEleLambdaMassWSSide1Away = new THnSparseF("fHistoEleLambdaMassWSSide1Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSide1Away);
  fHistoEleLambdaMassWSSide2Away = new THnSparseF("fHistoEleLambdaMassWSSide2Away","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassWSSide2Away);

  Int_t bins_base_elept[3]=		{10,100		,10};
  Double_t xmin_base_elept[3]={1.1,0		,0.00};
  Double_t xmax_base_elept[3]={3.1,10.	,100};
  fHistoEleLambdaMassvsElePtRS = new THnSparseF("fHistoEleLambdaMassvsElePtRS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRS);
  fHistoEleLambdaMassvsElePtWS = new THnSparseF("fHistoEleLambdaMassvsElePtWS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWS);
  fHistoEleLambdaMassvsElePtRSMix = new THnSparseF("fHistoEleLambdaMassvsElePtRSMix","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSMix);
  fHistoEleLambdaMassvsElePtWSMix = new THnSparseF("fHistoEleLambdaMassvsElePtWSMix","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSMix);
  fHistoEleLambdaMassvsElePtRSSide = new THnSparseF("fHistoEleLambdaMassvsElePtRSSide","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSSide);
  fHistoEleLambdaMassvsElePtWSSide = new THnSparseF("fHistoEleLambdaMassvsElePtWSSide","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSSide);
  fHistoEleLambdaMassvsElePtRS1 = new THnSparseF("fHistoEleLambdaMassvsElePtRS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRS1);
  fHistoEleLambdaMassvsElePtWS1 = new THnSparseF("fHistoEleLambdaMassvsElePtWS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWS1);
  fHistoEleLambdaMassvsElePtRSMix1 = new THnSparseF("fHistoEleLambdaMassvsElePtRSMix1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSMix1);
  fHistoEleLambdaMassvsElePtWSMix1 = new THnSparseF("fHistoEleLambdaMassvsElePtWSMix1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSMix1);
  fHistoEleLambdaMassvsElePtRSSide1 = new THnSparseF("fHistoEleLambdaMassvsElePtRSSide1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSSide1);
  fHistoEleLambdaMassvsElePtWSSide1 = new THnSparseF("fHistoEleLambdaMassvsElePtWSSide1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSSide1);
  fHistoEleLambdaMassvsElePtRS2 = new THnSparseF("fHistoEleLambdaMassvsElePtRS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRS2);
  fHistoEleLambdaMassvsElePtWS2 = new THnSparseF("fHistoEleLambdaMassvsElePtWS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWS2);
  fHistoEleLambdaMassvsElePtRSMix2 = new THnSparseF("fHistoEleLambdaMassvsElePtRSMix2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSMix2);
  fHistoEleLambdaMassvsElePtWSMix2 = new THnSparseF("fHistoEleLambdaMassvsElePtWSMix2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSMix2);
  fHistoEleLambdaMassvsElePtRSSide2 = new THnSparseF("fHistoEleLambdaMassvsElePtRSSide2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtRSSide2);
  fHistoEleLambdaMassvsElePtWSSide2 = new THnSparseF("fHistoEleLambdaMassvsElePtWSSide2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtWSSide2);

  fHistoEleLambdaMassMCS = new THnSparseF("fHistoEleLambdaMassMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassMCS);
  fHistoEleLambdaMassLambdabMCS = new THnSparseF("fHistoEleLambdaMassLambdabMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassLambdabMCS);
  fHistoEleLambdaMassXibMCS = new THnSparseF("fHistoEleLambdaMassXibMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassXibMCS);
  fHistoEleLambdaMassMCGen = new THnSparseF("fHistoEleLambdaMassMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassMCGen);
  fHistoEleLambdaMassvsElePtMCS = new THnSparseF("fHistoEleLambdaMassvsElePtMCS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCS);
  fHistoEleLambdaMassvsElePtMCGen = new THnSparseF("fHistoEleLambdaMassvsElePtMCGen","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCGen);
  fHistoEleLambdaMassvsElePtMCS1 = new THnSparseF("fHistoEleLambdaMassvsElePtMCS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCS1);
  fHistoEleLambdaMassvsElePtMCGen1 = new THnSparseF("fHistoEleLambdaMassvsElePtMCGen1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCGen1);
  fHistoEleLambdaMassvsElePtMCS2 = new THnSparseF("fHistoEleLambdaMassvsElePtMCS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCS2);
  fHistoEleLambdaMassvsElePtMCGen2 = new THnSparseF("fHistoEleLambdaMassvsElePtMCGen2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtMCGen2);

  fHistoElePtRS = new TH2F("fHistoElePtRS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtRS);
  fHistoElePtWS = new TH2F("fHistoElePtWS","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtWS);
  fHistoElePtRSMix = new TH2F("fHistoElePtRSMix","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtRSMix);
  fHistoElePtWSMix = new TH2F("fHistoElePtWSMix","",100,0.,10.,10,0,100);
  fOutputAll->Add(fHistoElePtWSMix);
  fHistoElePtMCS = new TH2F("fHistoElePtMCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtMCS);
  fHistoElePtMCGen = new TH2F("fHistoElePtMCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtMCGen);

  Int_t bins_eleptvseta[3]=		{10,10	,10};
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

  Int_t bins_eleptvslambdapt[3]=	{10,10	,10};
  Double_t xmin_eleptvslambdapt[3]={0.,0.	,0.0};
  Double_t xmax_eleptvslambdapt[3]={5.,5.	,100};

  fHistoElePtvsLambdaPtRS = new THnSparseF("fHistoElePtvsLambdaPtRS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtRS);
  fHistoElePtvsLambdaPtWS = new THnSparseF("fHistoElePtvsLambdaPtWS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtWS);
  fHistoElePtvsLambdaPtRSMix = new THnSparseF("fHistoElePtvsLambdaPtRSMix","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtRSMix);
  fHistoElePtvsLambdaPtWSMix = new THnSparseF("fHistoElePtvsLambdaPtWSMix","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtWSMix);
  fHistoElePtvsLambdaPtMCS = new THnSparseF("fHistoElePtvsLambdaPtMCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCS);
  fHistoElePtvsLambdaPtMCGen = new THnSparseF("fHistoElePtvsLambdaPtMCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCGen);
  fHistoElePtvsLambdaPtMCLcGen = new THnSparseF("fHistoElePtvsLambdaPtMCLcGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtMCLcGen);

  Int_t bins_eleptvslambdaptvslcpt[4]=	{10,10,10,10};
  Double_t xmin_eleptvslambdaptvslcpt[4]={0.,0.,0.,0.0};
  Double_t xmax_eleptvslambdaptvslcpt[4]={5.,5.,10.,100};
  fHistoElePtvsLambdaPtvsLcPtMCS = new THnSparseF("fHistoElePtvsLambdaPtvsLcPtMCS","",4,bins_eleptvslambdaptvslcpt,xmin_eleptvslambdaptvslcpt,xmax_eleptvslambdaptvslcpt);
  fOutputAll->Add(fHistoElePtvsLambdaPtvsLcPtMCS);
  fHistoElePtvsLambdaPtvsLcPtMCGen = new THnSparseF("fHistoElePtvsLambdaPtvsLcPtMCGen","",4,bins_eleptvslambdaptvslcpt,xmin_eleptvslambdaptvslcpt,xmax_eleptvslambdaptvslcpt);
  fOutputAll->Add(fHistoElePtvsLambdaPtvsLcPtMCGen);

  Int_t bins_allpt[4]=	{10,10,20,20};
  Double_t xmin_allpt[4]={0.,0.,0.,0.0};
  Double_t xmax_allpt[4]={20.,20.,10.,10};
  fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt = new THnSparseF("fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt","",4,bins_allpt,xmin_allpt,xmax_allpt);
  fOutputAll->Add(fHistoLcPtvseleLambdaPtvsElePtvsLambdaPt);

  fHistoEleLambdaPtvsRapidityRS=new TH2F("fHistoEleLambdaPtvsRapidityRS","EleLambda pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleLambdaPtvsRapidityRS);
  fHistoEleLambdaPtvsRapidityWS=new TH2F("fHistoEleLambdaPtvsRapidityWS","EleLambda pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleLambdaPtvsRapidityWS);
  fHistoEleLambdaPtvsRapidityMCS=new TH2F("fHistoEleLambdaPtvsRapidityMCS","EleLambda pt vs rap",20,0.,20.,40,-2.,2.);
  fOutputAll->Add(fHistoEleLambdaPtvsRapidityMCS);

  Int_t bins_eleptvsd0[3]=	{10 ,10	,10};
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


	//Feeddown from Xic0
  fHistoEleLambdaMassFeeddownXic0MCS = new THnSparseF("fHistoEleLambdaMassFeeddownXic0MCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXic0MCS);
  fHistoEleLambdaMassFeeddownXic0MCGen = new THnSparseF("fHistoEleLambdaMassFeeddownXic0MCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXic0MCGen);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCS);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCGen","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCGen);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS1 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCS1);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCGen1);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCS2 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCS2);
  fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXic0MCGen2);
  fHistoElePtFeeddownXic0MCS = new TH2F("fHistoElePtFeeddownXic0MCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXic0MCS);
  fHistoElePtFeeddownXic0MCGen = new TH2F("fHistoElePtFeeddownXic0MCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXic0MCGen);
  fHistoElePtvsEtaFeeddownXic0MCS = new THnSparseF("fHistoElePtvsEtaFeeddownXic0MCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXic0MCS);
  fHistoElePtvsEtaFeeddownXic0MCGen = new THnSparseF("fHistoElePtvsEtaFeeddownXic0MCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXic0MCGen);
  fHistoElePtvsLambdaPtFeeddownXic0MCS = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXic0MCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXic0MCS);
  fHistoElePtvsLambdaPtFeeddownXic0MCGen = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXic0MCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXic0MCGen);

	//Feeddown from XicPlus
  fHistoEleLambdaMassFeeddownXicPlusMCS = new THnSparseF("fHistoEleLambdaMassFeeddownXicPlusMCS","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXicPlusMCS);
  fHistoEleLambdaMassFeeddownXicPlusMCGen = new THnSparseF("fHistoEleLambdaMassFeeddownXicPlusMCGen","",3,bins_base,xmin_base,xmax_base);
  fOutputAll->Add(fHistoEleLambdaMassFeeddownXicPlusMCGen);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS1);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen1);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCS2);
  fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2 = new THnSparseF("fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2","",3,bins_base_elept,xmin_base_elept,xmax_base_elept);
  fOutputAll->Add(fHistoEleLambdaMassvsElePtFeeddownXicPlusMCGen2);
  fHistoElePtFeeddownXicPlusMCS = new TH2F("fHistoElePtFeeddownXicPlusMCS","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXicPlusMCS);
  fHistoElePtFeeddownXicPlusMCGen = new TH2F("fHistoElePtFeeddownXicPlusMCGen","",100,0,10,10,0,100);
  fOutputAll->Add(fHistoElePtFeeddownXicPlusMCGen);
  fHistoElePtvsEtaFeeddownXicPlusMCS = new THnSparseF("fHistoElePtvsEtaFeeddownXicPlusMCS","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXicPlusMCS);
  fHistoElePtvsEtaFeeddownXicPlusMCGen = new THnSparseF("fHistoElePtvsEtaFeeddownXicPlusMCGen","",3,bins_eleptvseta,xmin_eleptvseta,xmax_eleptvseta);
  fOutputAll->Add(fHistoElePtvsEtaFeeddownXicPlusMCGen);
  fHistoElePtvsLambdaPtFeeddownXicPlusMCS = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXicPlusMCS","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXicPlusMCS);
  fHistoElePtvsLambdaPtFeeddownXicPlusMCGen = new THnSparseF("fHistoElePtvsLambdaPtFeeddownXicPlusMCGen","",3,bins_eleptvslambdapt,xmin_eleptvslambdapt,xmax_eleptvslambdapt);
  fOutputAll->Add(fHistoElePtvsLambdaPtFeeddownXicPlusMCGen);

  //------------------------------------------------
  // checking histograms
  //------------------------------------------------
  fHistoBachPt = new TH1F("fHistoBachPt","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPt);
  fHistoBachPtMCS = new TH1F("fHistoBachPtMCS","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCS);
  fHistoBachPtMCGen = new TH1F("fHistoBachPtMCGen","Bachelor p_{T}",100,0.0,5.0);
  fOutputAll->Add(fHistoBachPtMCGen);
  fHistoBachFracSharedITS = new TH1F("fHistoBachFracSharedITS","Fraction shared",100,0.0,1.0);
  fOutputAll->Add(fHistoBachFracSharedITS);
  fHistoBachChi2NclsITS = new TH1F("fHistoBachChi2NclsITS","#chi^{2}/Ncls",100,0.0,50.0);
  fOutputAll->Add(fHistoBachChi2NclsITS);
  fHistod0Bach = new TH1F("fHistod0Bach","Bachelor d_{0}",100,-0.5,0.5);
  fOutputAll->Add(fHistod0Bach);
  fHistoLambdaMassvsPt=new TH2F("fHistoLambdaMassvsPt","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPt);
  fHistoLambdaMassvsPtMCS=new TH2F("fHistoLambdaMassvsPtMCS","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPtMCS);
  fHistoLambdaMassvsPtMCGen=new TH2F("fHistoLambdaMassvsPtMCGen","Lambda mass",100,1.116-0.05,1.116+0.05,20,0.,10.);
  fOutputAll->Add(fHistoLambdaMassvsPtMCGen);
  fHistoLambdaPtvsEta=new TH2F("fHistoLambdaPtvsEta","Lambda mass",20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoLambdaPtvsEta);
  fHistoLambdaPtvsEtaMCS=new TH2F("fHistoLambdaPtvsEtaMCS","Lambda mass",20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoLambdaPtvsEtaMCS);
  fHistoLambdaPtvsEtaMCGen=new TH2F("fHistoLambdaPtvsEtaMCGen","Lambda mass",20,0.,10.,20,-1.,1.);
  fOutputAll->Add(fHistoLambdaPtvsEtaMCGen);
  fHistoK0sMassvsPt=new TH2F("fHistoK0sMassvsPt","K0s mass",100,0.497-0.05,0.497+0.05,20,0.,10.);
  fOutputAll->Add(fHistoK0sMassvsPt);
  fHistoLambdaPtvsDl=new TH2F("fHistoLambdaPtvsDl","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDl);
  fHistoLambdaPtvsDlSide=new TH2F("fHistoLambdaPtvsDlSide","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlSide);
  fHistoLambdaPtvsDlMCS=new TH2F("fHistoLambdaPtvsDlMCS","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlMCS);
  fHistoLambdaPtvsDlFeeddownXi0MCS=new TH2F("fHistoLambdaPtvsDlFeeddownXi0MCS","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlFeeddownXi0MCS);
  fHistoLambdaPtvsDlFeeddownXiMinusMCS=new TH2F("fHistoLambdaPtvsDlFeeddownXiMinusMCS","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlFeeddownXiMinusMCS);
  fHistoLambdaPtvsDlFeeddownOmegaMCS=new TH2F("fHistoLambdaPtvsDlFeeddownOmegaMCS","Lambda pt vs dl",20,0.,10.,20,0.,40.);
  fOutputAll->Add(fHistoLambdaPtvsDlFeeddownOmegaMCS);

  Int_t bins_armenteros[3]=	{20 ,100	,200};
  Double_t xmin_armenteros[3]={0.,0.	,-1.};
  Double_t xmax_armenteros[3]={20.,1.	,1.};
  fHistoLambdaArmenteros = new THnSparseF("fHistoLambdaArmenteros","",3,bins_armenteros,xmin_armenteros,xmax_armenteros);
  fOutputAll->Add(fHistoLambdaArmenteros);
  fHistoLambdaArmenterosMCS = new THnSparseF("fHistoLambdaArmenterosMCS","",3,bins_armenteros,xmin_armenteros,xmax_armenteros);
  fOutputAll->Add(fHistoLambdaArmenterosMCS);

  fHistoElectronTPCPID=new TH2F("fHistoElectronTPCPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCPID);
  fHistoElectronTOFPID=new TH2F("fHistoElectronTOFPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFPID);
  fHistoElectronITSPID=new TH2F("fHistoElectronITSPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronITSPID);
  fHistoElectronTPCSelPID=new TH2F("fHistoElectronTPCSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTPCSelPID);
  fHistoElectronTOFSelPID=new TH2F("fHistoElectronTOFSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronTOFSelPID);
  fHistoElectronITSSelPID=new TH2F("fHistoElectronITSSelPID","",50,0.,5.,50,-20.,20.);
  fOutputAll->Add(fHistoElectronITSSelPID);
  fHistoElectronTPCPIDSelTOF=new TH2F("fHistoElectronTPCPIDSelTOF","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOF);
  fHistoElectronTPCPIDSelTOFSmallEta=new TH2F("fHistoElectronTPCPIDSelTOFSmallEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFSmallEta);
  fHistoElectronTPCPIDSelTOFLargeEta=new TH2F("fHistoElectronTPCPIDSelTOFLargeEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronTPCPIDSelTOFLargeEta);
  fHistoElectronITSPIDSelTPCTOF=new TH2F("fHistoElectronITSPIDSelTPCTOF","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronITSPIDSelTPCTOF);
  fHistoElectronITSPIDSelTPCTOFSmallEta=new TH2F("fHistoElectronITSPIDSelTPCTOFSmallEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronITSPIDSelTPCTOFSmallEta);
  fHistoElectronITSPIDSelTPCTOFLargeEta=new TH2F("fHistoElectronITSPIDSelTPCTOFLargeEta","",10,0.,5.,500,-10.,10.);
  fOutputAll->Add(fHistoElectronITSPIDSelTPCTOFLargeEta);

  fHistoMassConversionsMin=new TH1F("fHistoMassConversionsMin","",500,0,0.5);
  fOutputAll->Add(fHistoMassConversionsMin);
  fHistoMassConversionsSameSignMin=new TH1F("fHistoMassConversionsSameSignMin","",500,0,0.5);
  fOutputAll->Add(fHistoMassConversionsSameSignMin);

	for(Int_t i=0;i<8;i++){
		fHistoElectronTPCPIDSelTOFEtaDep[i]=new TH2F(Form("fHistoElectronTPCPIDSelTOFEtaDep[%d]",i),"",10,0.,5.,500,-10.,10.);
		fOutputAll->Add(fHistoElectronTPCPIDSelTOFEtaDep[i]);
		fHistoElectronITSPIDSelTPCTOFEtaDep[i]=new TH2F(Form("fHistoElectronITSPIDSelTPCTOFEtaDep[%d]",i),"",10,0.,5.,500,-10.,10.);
		fOutputAll->Add(fHistoElectronITSPIDSelTPCTOFEtaDep[i]);
	}
  fHistoElectronQovPtvsPhi=new TH2F("fHistoElectronQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoElectronQovPtvsPhi);
  fHistoLambdaQovPtvsPhi=new TH2F("fHistoLambdaQovPtvsPhi","",70,0.,7.,50,-2.,2.);
  fOutputAll->Add(fHistoLambdaQovPtvsPhi);

  Int_t bins_lcmcgen_withelea[3]=	{20 ,100	,100};
  Double_t xmin_lcmcgen_withelea[3]={0.,-5.,-5.};
  Double_t xmax_lcmcgen_withelea[3]={20.,5.,5.};
  fHistoLcMCGenWithEleA = new THnSparseF("fHistoLcMCGenWithEleA","",3,bins_lcmcgen_withelea,xmin_lcmcgen_withelea,xmax_lcmcgen_withelea);
  fOutputAll->Add(fHistoLcMCGenWithEleA);
  fHistoLcMCGenWithEleB = new TH2D("fHistoLcMCGenWithEleB","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoLcMCGenWithEleB);

  fHistoLbMCGenWithEleA = new THnSparseF("fHistoLbMCGenWithEleA","",3,bins_lcmcgen_withelea,xmin_lcmcgen_withelea,xmax_lcmcgen_withelea);
  fOutputAll->Add(fHistoLbMCGenWithEleA);
  fHistoLbMCGenWithEleB = new TH2D("fHistoLbMCGenWithEleB","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoLbMCGenWithEleB);

  Int_t bins_lbmcgen_withlc[3]=	{50 ,100	,100};
  Double_t xmin_lbmcgen_withlc[3]={0.,-5.,-5.};
  Double_t xmax_lbmcgen_withlc[3]={50.,5.,5.};
  fHistoLbMCGenWithLambdac = new THnSparseF("fHistoLbMCGenWithLambdac","",3,bins_lbmcgen_withlc,xmin_lbmcgen_withlc,xmax_lbmcgen_withlc);
  fOutputAll->Add(fHistoLbMCGenWithLambdac);

  Int_t bins_lcmcgen[3]=	{100 ,20	,10};
  Double_t xmin_lcmcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_lcmcgen[3]={20.,1.0	,100};
  Int_t bins_lbmcgen[3]=	{100 ,20	,10};
  Double_t xmin_lbmcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_lbmcgen[3]={50.,1.0	,100};
  fHistoLcMCGen = new THnSparseF("fHistoLcMCGen","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCGen);
  fHistoLcMCGen1 = new THnSparseF("fHistoLcMCGen1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCGen1);
  fHistoLcMCGen2 = new THnSparseF("fHistoLcMCGen2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCGen2);
  fHistoLbMCGen = new THnSparseF("fHistoLbMCGen","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCGen);
  fHistoLbMCGen1 = new THnSparseF("fHistoLbMCGen1","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCGen1);
  fHistoLbMCGen2 = new THnSparseF("fHistoLbMCGen2","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCGen2);
  fHistoFeedDownXic0MCGen = new THnSparseF("fHistoFeedDownXic0MCGen","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCGen);
  fHistoFeedDownXic0MCGen1 = new THnSparseF("fHistoFeedDownXic0MCGen1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCGen1);
  fHistoFeedDownXic0MCGen2 = new THnSparseF("fHistoFeedDownXic0MCGen2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCGen2);
  fHistoFeedDownXicPlusMCGen = new THnSparseF("fHistoFeedDownXicPlusMCGen","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCGen);
  fHistoFeedDownXicPlusMCGen1 = new THnSparseF("fHistoFeedDownXicPlusMCGen1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCGen1);
  fHistoFeedDownXicPlusMCGen2 = new THnSparseF("fHistoFeedDownXicPlusMCGen2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCGen2);
  fHistoFeedDownSigmaMCGen = new THnSparseF("fHistoFeedDownSigmaMCGen","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCGen);
  fHistoFeedDownSigmaMCGen1 = new THnSparseF("fHistoFeedDownSigmaMCGen1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCGen1);
  fHistoFeedDownSigmaMCGen2 = new THnSparseF("fHistoFeedDownSigmaMCGen2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCGen2);

  fHistoLcMCS = new THnSparseF("fHistoLcMCS","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCS);
  fHistoLcMCS1 = new THnSparseF("fHistoLcMCS1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCS1);
  fHistoLcMCS2 = new THnSparseF("fHistoLcMCS2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoLcMCS2);
  fHistoLbMCS = new THnSparseF("fHistoLbMCS","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCS);
  fHistoLbMCS1 = new THnSparseF("fHistoLbMCS1","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCS1);
  fHistoLbMCS2 = new THnSparseF("fHistoLbMCS2","",3,bins_lbmcgen,xmin_lbmcgen,xmax_lbmcgen);
  fOutputAll->Add(fHistoLbMCS2);
  fHistoFeedDownXic0MCS = new THnSparseF("fHistoFeedDownXic0MCS","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCS);
  fHistoFeedDownXic0MCS1 = new THnSparseF("fHistoFeedDownXic0MCS1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCS1);
  fHistoFeedDownXic0MCS2 = new THnSparseF("fHistoFeedDownXic0MCS2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXic0MCS2);
  fHistoFeedDownXicPlusMCS = new THnSparseF("fHistoFeedDownXicPlusMCS","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCS);
  fHistoFeedDownXicPlusMCS1 = new THnSparseF("fHistoFeedDownXicPlusMCS1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCS1);
  fHistoFeedDownXicPlusMCS2 = new THnSparseF("fHistoFeedDownXicPlusMCS2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownXicPlusMCS2);
  fHistoFeedDownSigmaMCS = new THnSparseF("fHistoFeedDownSigmaMCS","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCS);
  fHistoFeedDownSigmaMCS1 = new THnSparseF("fHistoFeedDownSigmaMCS1","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCS1);
  fHistoFeedDownSigmaMCS2 = new THnSparseF("fHistoFeedDownSigmaMCS2","",3,bins_lcmcgen,xmin_lcmcgen,xmax_lcmcgen);
  fOutputAll->Add(fHistoFeedDownSigmaMCS2);

  Int_t bins_elemcgen[3]=	{100 ,20	,10};
  Double_t xmin_elemcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_elemcgen[3]={10.,1.0	,100};
  fHistoElectronMCGen = new THnSparseF("fHistoElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronMCGen);
  fHistoBottomElectronMCGen = new THnSparseF("fHistoBottomElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoBottomElectronMCGen);
  fHistoCharmElectronMCGen = new THnSparseF("fHistoCharmElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoCharmElectronMCGen);
  fHistoLcElectronMCGen = new THnSparseF("fHistoLcElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCGen);
  fHistoLcElectronMCGen1 = new THnSparseF("fHistoLcElectronMCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCGen1);
  fHistoLcElectronMCGen2 = new THnSparseF("fHistoLcElectronMCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCGen2);
  fHistoLbElectronMCGen = new THnSparseF("fHistoLbElectronMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCGen);
  fHistoLbElectronMCGen1 = new THnSparseF("fHistoLbElectronMCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCGen1);
  fHistoLbElectronMCGen2 = new THnSparseF("fHistoLbElectronMCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCGen2);
  fHistoElectronFeedDownXic0MCGen = new THnSparseF("fHistoElectronFeedDownXic0MCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCGen);
  fHistoElectronFeedDownXic0MCGen1 = new THnSparseF("fHistoElectronFeedDownXic0MCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCGen1);
  fHistoElectronFeedDownXic0MCGen2 = new THnSparseF("fHistoElectronFeedDownXic0MCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCGen2);
  fHistoElectronFeedDownXicPlusMCGen = new THnSparseF("fHistoElectronFeedDownXicPlusMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCGen);
  fHistoElectronFeedDownXicPlusMCGen1 = new THnSparseF("fHistoElectronFeedDownXicPlusMCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCGen1);
  fHistoElectronFeedDownXicPlusMCGen2 = new THnSparseF("fHistoElectronFeedDownXicPlusMCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCGen2);
  fHistoElectronFeedDownXibMCGen = new THnSparseF("fHistoElectronFeedDownXibMCGen","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCGen);
  fHistoElectronFeedDownXibMCGen1 = new THnSparseF("fHistoElectronFeedDownXibMCGen1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCGen1);
  fHistoElectronFeedDownXibMCGen2 = new THnSparseF("fHistoElectronFeedDownXibMCGen2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCGen2);
  fHistoLcElectronMCS = new THnSparseF("fHistoLcElectronMCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCS);
  fHistoLcElectronMCS1 = new THnSparseF("fHistoLcElectronMCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCS1);
  fHistoLcElectronMCS2 = new THnSparseF("fHistoLcElectronMCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLcElectronMCS2);
  fHistoLbElectronMCS = new THnSparseF("fHistoLbElectronMCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCS);
  fHistoLbElectronMCS1 = new THnSparseF("fHistoLbElectronMCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCS1);
  fHistoLbElectronMCS2 = new THnSparseF("fHistoLbElectronMCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoLbElectronMCS2);
  fHistoElectronFeedDownXic0MCS = new THnSparseF("fHistoElectronFeedDownXic0MCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCS);
  fHistoElectronFeedDownXic0MCS1 = new THnSparseF("fHistoElectronFeedDownXic0MCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCS1);
  fHistoElectronFeedDownXic0MCS2 = new THnSparseF("fHistoElectronFeedDownXic0MCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXic0MCS2);
  fHistoElectronFeedDownXicPlusMCS = new THnSparseF("fHistoElectronFeedDownXicPlusMCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCS);
  fHistoElectronFeedDownXicPlusMCS1 = new THnSparseF("fHistoElectronFeedDownXicPlusMCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCS1);
  fHistoElectronFeedDownXicPlusMCS2 = new THnSparseF("fHistoElectronFeedDownXicPlusMCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXicPlusMCS2);
  fHistoElectronFeedDownXibMCS = new THnSparseF("fHistoElectronFeedDownXibMCS","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCS);
  fHistoElectronFeedDownXibMCS1 = new THnSparseF("fHistoElectronFeedDownXibMCS1","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCS1);
  fHistoElectronFeedDownXibMCS2 = new THnSparseF("fHistoElectronFeedDownXibMCS2","",3,bins_elemcgen,xmin_elemcgen,xmax_elemcgen);
  fOutputAll->Add(fHistoElectronFeedDownXibMCS2);

  Int_t bins_lambdamcgen[3]=	{50 ,20	,10};
  Double_t xmin_lambdamcgen[3]={0.,-1.0	,0.0};
  Double_t xmax_lambdamcgen[3]={10.,1.0	,100};
  fHistoLambdaMCGen = new THnSparseF("fHistoLambdaMCGen","",3,bins_lambdamcgen,xmin_lambdamcgen,xmax_lambdamcgen);
  fOutputAll->Add(fHistoLambdaMCGen);

  Int_t bins_eleptvsv0dl[3]=	{100 ,20	,10};
  Double_t xmin_eleptvsv0dl[3]={0.,0.	,0.0};
  Double_t xmax_eleptvsv0dl[3]={10.,40.	,100};
  fHistoElePtvsV0dlRS = new THnSparseF("fHistoElePtvsV0dlRS","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRS);
  fHistoElePtvsV0dlRS1 = new THnSparseF("fHistoElePtvsV0dlRS1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRS1);
  fHistoElePtvsV0dlRS2 = new THnSparseF("fHistoElePtvsV0dlRS2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRS2);
  fHistoElePtvsV0dlRSSide = new THnSparseF("fHistoElePtvsV0dlRSSide","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSSide);
  fHistoElePtvsV0dlRSSide1 = new THnSparseF("fHistoElePtvsV0dlRSSide1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSSide1);
  fHistoElePtvsV0dlRSSide2 = new THnSparseF("fHistoElePtvsV0dlRSSide2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSSide2);
  fHistoElePtvsV0dlRSMix = new THnSparseF("fHistoElePtvsV0dlRSMix","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSMix);
  fHistoElePtvsV0dlRSMix1 = new THnSparseF("fHistoElePtvsV0dlRSMix1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSMix1);
  fHistoElePtvsV0dlRSMix2 = new THnSparseF("fHistoElePtvsV0dlRSMix2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlRSMix2);
  fHistoElePtvsV0dlWS = new THnSparseF("fHistoElePtvsV0dlWS","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWS);
  fHistoElePtvsV0dlWS1 = new THnSparseF("fHistoElePtvsV0dlWS1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWS1);
  fHistoElePtvsV0dlWS2 = new THnSparseF("fHistoElePtvsV0dlWS2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWS2);
  fHistoElePtvsV0dlWSSide = new THnSparseF("fHistoElePtvsV0dlWSSide","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSSide);
  fHistoElePtvsV0dlWSSide1 = new THnSparseF("fHistoElePtvsV0dlWSSide1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSSide1);
  fHistoElePtvsV0dlWSSide2 = new THnSparseF("fHistoElePtvsV0dlWSSide2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSSide2);
  fHistoElePtvsV0dlWSMix = new THnSparseF("fHistoElePtvsV0dlWSMix","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSMix);
  fHistoElePtvsV0dlWSMix1 = new THnSparseF("fHistoElePtvsV0dlWSMix1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSMix1);
  fHistoElePtvsV0dlWSMix2 = new THnSparseF("fHistoElePtvsV0dlWSMix2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlWSMix2);
  fHistoElePtvsV0dlMCS = new THnSparseF("fHistoElePtvsV0dlMCS","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlMCS);
  fHistoElePtvsV0dlMCS1 = new THnSparseF("fHistoElePtvsV0dlMCS1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlMCS1);
  fHistoElePtvsV0dlMCS2 = new THnSparseF("fHistoElePtvsV0dlMCS2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlMCS2);
  fHistoElePtvsV0dlFeedDownXic0MCS = new THnSparseF("fHistoElePtvsV0dlFeedDownXic0MCS","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXic0MCS);
  fHistoElePtvsV0dlFeedDownXic0MCS1 = new THnSparseF("fHistoElePtvsV0dlFeedDownXic0MCS1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXic0MCS1);
  fHistoElePtvsV0dlFeedDownXic0MCS2 = new THnSparseF("fHistoElePtvsV0dlFeedDownXic0MCS2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXic0MCS2);
  fHistoElePtvsV0dlFeedDownXicPlusMCS = new THnSparseF("fHistoElePtvsV0dlFeedDownXicPlusMCS","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXicPlusMCS);
  fHistoElePtvsV0dlFeedDownXicPlusMCS1 = new THnSparseF("fHistoElePtvsV0dlFeedDownXicPlusMCS1","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXicPlusMCS1);
  fHistoElePtvsV0dlFeedDownXicPlusMCS2 = new THnSparseF("fHistoElePtvsV0dlFeedDownXicPlusMCS2","",3,bins_eleptvsv0dl,xmin_eleptvsv0dl,xmax_eleptvsv0dl);
  fOutputAll->Add(fHistoElePtvsV0dlFeedDownXicPlusMCS2);

  Int_t bins_eleptvsv0dca[3]=	{100 ,20	,10};
  Double_t xmin_eleptvsv0dca[3]={0.,0.	,0.0};
  Double_t xmax_eleptvsv0dca[3]={10.,1.	,100};
  fHistoElePtvsV0dcaRS = new THnSparseF("fHistoElePtvsV0dcaRS","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRS);
  fHistoElePtvsV0dcaRS1 = new THnSparseF("fHistoElePtvsV0dcaRS1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRS1);
  fHistoElePtvsV0dcaRS2 = new THnSparseF("fHistoElePtvsV0dcaRS2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRS2);
  fHistoElePtvsV0dcaRSSide = new THnSparseF("fHistoElePtvsV0dcaRSSide","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSSide);
  fHistoElePtvsV0dcaRSSide1 = new THnSparseF("fHistoElePtvsV0dcaRSSide1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSSide1);
  fHistoElePtvsV0dcaRSSide2 = new THnSparseF("fHistoElePtvsV0dcaRSSide2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSSide2);
  fHistoElePtvsV0dcaRSMix = new THnSparseF("fHistoElePtvsV0dcaRSMix","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSMix);
  fHistoElePtvsV0dcaRSMix1 = new THnSparseF("fHistoElePtvsV0dcaRSMix1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSMix1);
  fHistoElePtvsV0dcaRSMix2 = new THnSparseF("fHistoElePtvsV0dcaRSMix2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaRSMix2);
  fHistoElePtvsV0dcaWS = new THnSparseF("fHistoElePtvsV0dcaWS","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWS);
  fHistoElePtvsV0dcaWS1 = new THnSparseF("fHistoElePtvsV0dcaWS1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWS1);
  fHistoElePtvsV0dcaWS2 = new THnSparseF("fHistoElePtvsV0dcaWS2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWS2);
  fHistoElePtvsV0dcaWSSide = new THnSparseF("fHistoElePtvsV0dcaWSSide","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSSide);
  fHistoElePtvsV0dcaWSSide1 = new THnSparseF("fHistoElePtvsV0dcaWSSide1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSSide1);
  fHistoElePtvsV0dcaWSSide2 = new THnSparseF("fHistoElePtvsV0dcaWSSide2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSSide2);
  fHistoElePtvsV0dcaWSMix = new THnSparseF("fHistoElePtvsV0dcaWSMix","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSMix);
  fHistoElePtvsV0dcaWSMix1 = new THnSparseF("fHistoElePtvsV0dcaWSMix1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSMix1);
  fHistoElePtvsV0dcaWSMix2 = new THnSparseF("fHistoElePtvsV0dcaWSMix2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaWSMix2);
  fHistoElePtvsV0dcaMCS = new THnSparseF("fHistoElePtvsV0dcaMCS","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaMCS);
  fHistoElePtvsV0dcaMCS1 = new THnSparseF("fHistoElePtvsV0dcaMCS1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaMCS1);
  fHistoElePtvsV0dcaMCS2 = new THnSparseF("fHistoElePtvsV0dcaMCS2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaMCS2);
  fHistoElePtvsV0dcaFeedDownXic0MCS = new THnSparseF("fHistoElePtvsV0dcaFeedDownXic0MCS","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXic0MCS);
  fHistoElePtvsV0dcaFeedDownXic0MCS1 = new THnSparseF("fHistoElePtvsV0dcaFeedDownXic0MCS1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXic0MCS1);
  fHistoElePtvsV0dcaFeedDownXic0MCS2 = new THnSparseF("fHistoElePtvsV0dcaFeedDownXic0MCS2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXic0MCS2);
  fHistoElePtvsV0dcaFeedDownXicPlusMCS = new THnSparseF("fHistoElePtvsV0dcaFeedDownXicPlusMCS","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXicPlusMCS);
  fHistoElePtvsV0dcaFeedDownXicPlusMCS1 = new THnSparseF("fHistoElePtvsV0dcaFeedDownXicPlusMCS1","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXicPlusMCS1);
  fHistoElePtvsV0dcaFeedDownXicPlusMCS2 = new THnSparseF("fHistoElePtvsV0dcaFeedDownXicPlusMCS2","",3,bins_eleptvsv0dca,xmin_eleptvsv0dca,xmax_eleptvsv0dca);
  fOutputAll->Add(fHistoElePtvsV0dcaFeedDownXicPlusMCS2);

  Int_t bins_elelamptvsv0dl[4]=	{100 ,20	,100, 10};
  Double_t xmin_elelamptvsv0dl[4]={0.,0.	,-0.5, 0.0};
  Double_t xmax_elelamptvsv0dl[4]={20.,40.,0.5	,100};
  fHistoEleLambdaPtvsV0dlRS = new THnSparseF("fHistoEleLambdaPtvsV0dlRS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRS);
  fHistoEleLambdaPtvsV0dlRS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlRS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRS1);
  fHistoEleLambdaPtvsV0dlRS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlRS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRS2);
  fHistoEleLambdaPtvsV0dlRSSide = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSide","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSide);
  fHistoEleLambdaPtvsV0dlRSSide1 = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSide1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSide1);
  fHistoEleLambdaPtvsV0dlRSSide2 = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSide2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSide2);
  fHistoEleLambdaPtvsV0dlRSMix = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMix","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMix);
  fHistoEleLambdaPtvsV0dlRSMix1 = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMix1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMix1);
  fHistoEleLambdaPtvsV0dlRSMix2 = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMix2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMix2);
  fHistoEleLambdaPtvsV0dlWS = new THnSparseF("fHistoEleLambdaPtvsV0dlWS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWS);
  fHistoEleLambdaPtvsV0dlWS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlWS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWS1);
  fHistoEleLambdaPtvsV0dlWS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlWS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWS2);
  fHistoEleLambdaPtvsV0dlWSSide = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSide","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSide);
  fHistoEleLambdaPtvsV0dlWSSide1 = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSide1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSide1);
  fHistoEleLambdaPtvsV0dlWSSide2 = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSide2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSide2);
  fHistoEleLambdaPtvsV0dlWSMix = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMix","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMix);
  fHistoEleLambdaPtvsV0dlWSMix1 = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMix1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMix1);
  fHistoEleLambdaPtvsV0dlWSMix2 = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMix2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMix2);
  fHistoEleLambdaPtvsV0dlMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlMCS);
  fHistoEleLambdaPtvsV0dlMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlMCS1);
  fHistoEleLambdaPtvsV0dlMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlMCS2);
  fHistoEleLambdaPtvsV0dlPromptMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlPromptMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlPromptMCS);
  fHistoEleLambdaPtvsV0dlPromptMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlPromptMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlPromptMCS1);
  fHistoEleLambdaPtvsV0dlPromptMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlPromptMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlPromptMCS2);
  fHistoEleLambdaPtvsV0dlBFeeddownMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlBFeeddownMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlBFeeddownMCS);
  fHistoEleLambdaPtvsV0dlBFeeddownMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlBFeeddownMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlBFeeddownMCS1);
  fHistoEleLambdaPtvsV0dlBFeeddownMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlBFeeddownMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlBFeeddownMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownXic0MCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0MCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0MCS);
  fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0MCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0MCS2);
  fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS);
  fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0PromptMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS);
  fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXic0BFeeddownMCS2);

  fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusPromptMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownXicPlusBFeeddownMCS2);

  fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaPromptMCS2);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS1);
  fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlFeedDownSigmaBFeeddownMCS2);
  fHistoEleLambdaPtvsV0dlLambdabMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlLambdabMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlLambdabMCS);
  fHistoEleLambdaPtvsV0dlLambdabMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlLambdabMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlLambdabMCS1);
  fHistoEleLambdaPtvsV0dlLambdabMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlLambdabMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlLambdabMCS2);
  fHistoEleLambdaPtvsV0dlXibMCS = new THnSparseF("fHistoEleLambdaPtvsV0dlXibMCS","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlXibMCS);
  fHistoEleLambdaPtvsV0dlXibMCS1 = new THnSparseF("fHistoEleLambdaPtvsV0dlXibMCS1","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlXibMCS1);
  fHistoEleLambdaPtvsV0dlXibMCS2 = new THnSparseF("fHistoEleLambdaPtvsV0dlXibMCS2","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlXibMCS2);

  fHistoEleLambdaPtvsV0dlRSAway = new THnSparseF("fHistoEleLambdaPtvsV0dlRSAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSAway);
  fHistoEleLambdaPtvsV0dlRS1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRS1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRS1Away);
  fHistoEleLambdaPtvsV0dlRS2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRS2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRS2Away);
  fHistoEleLambdaPtvsV0dlRSSideAway = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSideAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSideAway);
  fHistoEleLambdaPtvsV0dlRSSide1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSide1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSide1Away);
  fHistoEleLambdaPtvsV0dlRSSide2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRSSide2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSSide2Away);
  fHistoEleLambdaPtvsV0dlRSMixAway = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMixAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMixAway);
  fHistoEleLambdaPtvsV0dlRSMix1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMix1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMix1Away);
  fHistoEleLambdaPtvsV0dlRSMix2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlRSMix2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlRSMix2Away);
  fHistoEleLambdaPtvsV0dlWSAway = new THnSparseF("fHistoEleLambdaPtvsV0dlWSAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSAway);
  fHistoEleLambdaPtvsV0dlWS1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWS1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWS1Away);
  fHistoEleLambdaPtvsV0dlWS2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWS2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWS2Away);
  fHistoEleLambdaPtvsV0dlWSSideAway = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSideAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSideAway);
  fHistoEleLambdaPtvsV0dlWSSide1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSide1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSide1Away);
  fHistoEleLambdaPtvsV0dlWSSide2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWSSide2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSSide2Away);
  fHistoEleLambdaPtvsV0dlWSMixAway = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMixAway","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMixAway);
  fHistoEleLambdaPtvsV0dlWSMix1Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMix1Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMix1Away);
  fHistoEleLambdaPtvsV0dlWSMix2Away = new THnSparseF("fHistoEleLambdaPtvsV0dlWSMix2Away","",4,bins_elelamptvsv0dl,xmin_elelamptvsv0dl,xmax_elelamptvsv0dl);
  fOutputAll->Add(fHistoEleLambdaPtvsV0dlWSMix2Away);

  fHistoResponseElePt = new TH2D("fHistoResponseElePt","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoResponseElePt);
  fHistoResponseElePt1 = new TH2D("fHistoResponseElePt1","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoResponseElePt1);
  fHistoResponseElePt2 = new TH2D("fHistoResponseElePt2","",100,0.,20.,100,0.,10.);
  fOutputAll->Add(fHistoResponseElePt2);
  fHistoResponseEleLambdaPt = new TH2D("fHistoResponseEleLambdaPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPt);
  fHistoResponseEleLambdaPt1 = new TH2D("fHistoResponseEleLambdaPt1","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPt1);
  fHistoResponseEleLambdaPt2 = new TH2D("fHistoResponseEleLambdaPt2","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPt2);
  fHistoResponseEleLambdaPtFeeddownXic0 = new TH2D("fHistoResponseEleLambdaPtFeeddownXic0","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXic0);
  fHistoResponseEleLambdaPtFeeddownXic01 = new TH2D("fHistoResponseEleLambdaPtFeeddownXic01","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXic01);
  fHistoResponseEleLambdaPtFeeddownXic02 = new TH2D("fHistoResponseEleLambdaPtFeeddownXic02","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXic02);
  fHistoResponseEleLambdaPtFeeddownXicPlus = new TH2D("fHistoResponseEleLambdaPtFeeddownXicPlus","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXicPlus);
  fHistoResponseEleLambdaPtFeeddownXicPlus1 = new TH2D("fHistoResponseEleLambdaPtFeeddownXicPlus1","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXicPlus1);
  fHistoResponseEleLambdaPtFeeddownXicPlus2 = new TH2D("fHistoResponseEleLambdaPtFeeddownXicPlus2","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownXicPlus2);
  fHistoResponseEleLambdaPtFeeddownSigma = new TH2D("fHistoResponseEleLambdaPtFeeddownSigma","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownSigma);
  fHistoResponseEleLambdaPtFeeddownSigma1 = new TH2D("fHistoResponseEleLambdaPtFeeddownSigma1","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownSigma1);
  fHistoResponseEleLambdaPtFeeddownSigma2 = new TH2D("fHistoResponseEleLambdaPtFeeddownSigma2","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtFeeddownSigma2);
  fHistoResponseEleLambdaPtLambdab = new TH2D("fHistoResponseEleLambdaPtLambdab","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseEleLambdaPtLambdab);
  fHistoResponseElePtLambdab = new TH2D("fHistoResponseElePtLambdab","",100,0.,50.,100,0.,10.);
  fOutputAll->Add(fHistoResponseElePtLambdab);
  fHistoResponseLambdaPtLambdab = new TH2D("fHistoResponseLambdaPtLambdab","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtLambdab);
  fHistoResponseLambdaPt = new TH2D("fHistoResponseLambdaPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPt);
  fHistoResponseLambdaPtFeeddownXic0 = new TH2D("fHistoResponseLambdaPtFeeddownXic0","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtFeeddownXic0);
  fHistoResponseLambdaPtFeeddownXicPlus = new TH2D("fHistoResponseLambdaPtFeeddownXicPlus","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtFeeddownXicPlus);
  fHistoResponseLambdaPtFeeddownSigma = new TH2D("fHistoResponseLambdaPtFeeddownSigma","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtFeeddownSigma);
  fHistoResponseLambdaPtvsEleLambdaPt = new TH2D("fHistoResponseLambdaPtvsEleLambdaPt","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtvsEleLambdaPt);
  fHistoResponseLambdaPtvsEleLambdaPtLambdab = new TH2D("fHistoResponseLambdaPtvsEleLambdaPtLambdab","",100,0.,20.,100,0.,20.);
  fOutputAll->Add(fHistoResponseLambdaPtvsEleLambdaPtLambdab);
  fHistoResponseMcGenLambdabPtvsLambdacPt = new TH2D("fHistoResponseMcGenLambdabPtvsLambdacPt","",100,0.,50.,100,0.,20.);
  fOutputAll->Add(fHistoResponseMcGenLambdabPtvsLambdacPt);

  fHistonEvtvsRunNumber=new TH1F("fHistonEvtvsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonEvtvsRunNumber);
  fHistonElevsRunNumber=new TH1F("fHistonElevsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonElevsRunNumber);
  fHistonLambdavsRunNumber=new TH1F("fHistonLambdavsRunNumber","",20000,-0.5,19999.5);
  fOutputAll->Add(fHistonLambdavsRunNumber);
  fHistoMCEventType=new TH1F("fHistoMCEventType","",4,-0.5,3.5);
  fOutputAll->Add(fHistoMCEventType);
  fHistoMCDeltaPhiccbar=new TH1F("fHistoMCDeltaPhiccbar","",100,0.,3.2);
  fOutputAll->Add(fHistoMCDeltaPhiccbar);

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

	for(Int_t ih=0;ih<17;ih++){
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
			//7: Rfid Lambda
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 5.;
		}else if(ih==8){
			//10: Dca V0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 2.;
		}else if(ih==9 || ih==10 ){
			//9: DCA V0pr to prim 10: DCA V0pi to prim
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.;
			xmax_eleptvscutvars[1] = 0.5;
		}else if(ih==11){
			//11: CosPAv0
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = 0.95;
			xmax_eleptvscutvars[1] = 1.0;
		}else if(ih==12){
			//12:K0s masss
			bins_eleptvscutvars[1] = 50;
			xmin_eleptvscutvars[1] = 0.497-0.03;
			xmax_eleptvscutvars[1] = 0.497+0.03;
		}else if(ih==13 || ih==14){
			//13: nSigmaTPC(pr), nSigma(pi)
			bins_eleptvscutvars[1] = 20;
			xmin_eleptvscutvars[1] = -5;
			xmax_eleptvscutvars[1] = 5;
		}else if(ih==15){
			//15: eta
			bins_eleptvscutvars[1] = 30;
			xmin_eleptvscutvars[1] = -1.5;
			xmax_eleptvscutvars[1] = 1.5;
		}else if(ih==16){
			//16: Opening angle
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
  fHistoCorrelationVariablesvsEleLambdaPt = new THnSparseF("fHistoCorrelationVariablesvsEleLambdaPt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsEleLambdaPtMix = new THnSparseF("fHistoCorrelationVariablesvsEleLambdaPtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsEleLambdaPtMC = new THnSparseF("fHistoCorrelationVariablesvsEleLambdaPtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsElePt = new THnSparseF("fHistoCorrelationVariablesvsElePt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsElePtMix = new THnSparseF("fHistoCorrelationVariablesvsElePtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsElePtMC = new THnSparseF("fHistoCorrelationVariablesvsElePtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd2);
  fHistoCorrelationVariablesvsLambdaPt = new THnSparseF("fHistoCorrelationVariablesvsLambdaPt","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsLambdaPtMix = new THnSparseF("fHistoCorrelationVariablesvsLambdaPtMix","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);
  fHistoCorrelationVariablesvsLambdaPtMC = new THnSparseF("fHistoCorrelationVariablesvsLambdaPtMC","",7,bins_cor_nd,xmin_cor_nd,xmax_cor_nd);

	//Axis 0: Mass
	//Axis 1: Pt
	//Axis 2: Near or Away
	//Axis 3: peak or Sideband
	//Axis 4: Sign Type
	//Axis 5: Conv Type
	//Axis 6: MC Type
	//Axis 7: Centrality
  Int_t bins_mass_nd[8]=	{23,100 , 2, 2, 4, 3, 10, 10};
  Double_t xmin_mass_nd[8]={1.1,0.,-0.5,-0.5,-0.5,-0.5,-0.5,0.};
  Double_t xmax_mass_nd[8]={5.7,20.,1.5,1.5,3.5,2.5,9.5,100.};
  Double_t xmax_mass_nd2[8]={5.7,10.,1.5,1.5,3.5,2.5,9.5,100.};
  fHistoMassVariablesvsEleLambdaPt = new THnSparseF("fHistoMassVariablesvsEleLambdaPt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsEleLambdaPtMix = new THnSparseF("fHistoMassVariablesvsEleLambdaPtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsEleLambdaPtMC = new THnSparseF("fHistoMassVariablesvsEleLambdaPtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsElePt = new THnSparseF("fHistoMassVariablesvsElePt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsElePtMix = new THnSparseF("fHistoMassVariablesvsElePtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsElePtMC = new THnSparseF("fHistoMassVariablesvsElePtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd2);
  fHistoMassVariablesvsLambdaPt = new THnSparseF("fHistoMassVariablesvsLambdaPt","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsLambdaPtMix = new THnSparseF("fHistoMassVariablesvsLambdaPtMix","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);
  fHistoMassVariablesvsLambdaPtMC = new THnSparseF("fHistoMassVariablesvsLambdaPtMC","",8,bins_mass_nd,xmin_mass_nd,xmax_mass_nd);

  fOutputAll->Add(fHistoCorrelationVariablesvsEleLambdaPt);
  fOutputAll->Add(fHistoCorrelationVariablesvsEleLambdaPtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsEleLambdaPtMC);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePt);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsElePtMC);
  fOutputAll->Add(fHistoCorrelationVariablesvsLambdaPt);
  fOutputAll->Add(fHistoCorrelationVariablesvsLambdaPtMix);
  fOutputAll->Add(fHistoCorrelationVariablesvsLambdaPtMC);

  fOutputAll->Add(fHistoMassVariablesvsEleLambdaPt);
  fOutputAll->Add(fHistoMassVariablesvsEleLambdaPtMix);
  fOutputAll->Add(fHistoMassVariablesvsEleLambdaPtMC);
  fOutputAll->Add(fHistoMassVariablesvsElePt);
  fOutputAll->Add(fHistoMassVariablesvsElePtMix);
  fOutputAll->Add(fHistoMassVariablesvsElePtMC);
  fOutputAll->Add(fHistoMassVariablesvsLambdaPt);
  fOutputAll->Add(fHistoMassVariablesvsLambdaPtMix);
  fOutputAll->Add(fHistoMassVariablesvsLambdaPtMC);

  return;
}

//________________________________________________________________________
AliAODRecoCascadeHF* AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeCascadeHF(AliAODv0 *v0, AliAODTrack *part, AliAODTrack *partpid, AliAODEvent * aod, AliAODVertex *secVert) 
{
  ///
  /// Create AliAODRecoCascadeHF object from the argument
  ///

  if(!v0) return 0x0;
  if(!part) return 0x0;
  if(!aod) return 0x0;

  //------------------------------------------------
  // PrimaryVertex
  //------------------------------------------------
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
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

  AliNeutralTrackParam *trackV0=NULL;
  const AliVTrack *trackVV0 = dynamic_cast<const AliVTrack*>(v0);
  if(trackVV0)  trackV0 = new AliNeutralTrackParam(trackVV0);

  Double_t xdummy, ydummy;
  Double_t dca = esdtrack->GetDCA(trackV0,fBzkG,xdummy,ydummy);


  //------------------------------------------------
  // Propagate all tracks to the secondary vertex and calculate momentum there
  //------------------------------------------------
	
  Double_t d0z0bach[2],covd0z0bach[3];
  if(sqrt(pow(secVert->GetX(),2)+pow(secVert->GetY(),2))<1.){
    part->PropagateToDCA(secVert,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(secVert,fBzkG,kVeryBig);
  }else{
    part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
    trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig);
  }
  Double_t momv0_new[3]={-9999,-9999,-9999.};
  trackV0->GetPxPyPz(momv0_new);

  Double_t px[2],py[2],pz[2];
  px[0] = part->Px(); py[0] = part->Py(); pz[0] = part->Pz(); 
  px[1] = momv0_new[0]; py[1] = momv0_new[1]; pz[1] = momv0_new[2]; 

  //------------------------------------------------
  // d0
  //------------------------------------------------
  Double_t d0[3],d0err[3];

  part->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0bach,covd0z0bach);
  d0[0]= d0z0bach[0];
  d0err[0] = TMath::Sqrt(covd0z0bach[0]);

  Double_t d0z0v0[2],covd0z0v0[3];
  trackV0->PropagateToDCA(primVertexAOD,fBzkG,kVeryBig,d0z0v0,covd0z0v0);
  d0[1]= d0z0v0[0];
  d0err[1] = TMath::Sqrt(covd0z0v0[0]);

  //------------------------------------------------
  // Create AliAODRecoCascadeHF
  //------------------------------------------------
  Short_t charge = part->Charge();
  AliAODRecoCascadeHF *theCascade = new AliAODRecoCascadeHF(secVert,charge,px,py,pz,d0,d0err,dca);
  if(!theCascade)  
    {
      if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
      if(esdtrack) delete esdtrack;
      if(trackV0) delete trackV0;
      return 0x0;
    }
  theCascade->SetOwnPrimaryVtx(primVertexAOD);
  UShort_t id[2]={(UShort_t)part->GetID(),(UShort_t)trackV0->GetID()};
  theCascade->SetProngIDs(2,id);

	theCascade->GetSecondaryVtx()->AddDaughter(part);
	theCascade->GetSecondaryVtx()->AddDaughter(v0);

  if(unsetvtx) delete primVertexAOD; primVertexAOD=NULL;
  if(esdtrack) delete esdtrack;
  if(trackV0) delete trackV0;

  return theCascade;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::CallPrimaryVertex(AliAODv0 *v0, AliAODTrack *trk, AliAODEvent* aod)
{
  //
  /// Make an array of tracks which should not be used in primary vertex calculation and
  /// Call PrimaryVertex function
  //

  TObjArray *TrackArray = new TObjArray(3);
  
  AliESDtrack *cptrk1 = new AliESDtrack((AliVTrack*)trk);
  TrackArray->AddAt(cptrk1,0);
  
  AliESDtrack *cascptrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(0));
  TrackArray->AddAt(cascptrack,1);
  AliESDtrack *cascntrack = new AliESDtrack((AliVTrack*)v0->GetDaughter(1));
  TrackArray->AddAt(cascntrack,2);
  
  AliAODVertex *newvert  = PrimaryVertex(TrackArray,aod);
  
  for(Int_t i=0;i<3;i++)
    {
      AliESDtrack *tesd = (AliESDtrack*)TrackArray->UncheckedAt(i);
      delete tesd;
    }
  TrackArray->Clear();
  delete TrackArray;
  
  return newvert;
}

//________________________________________________________________________
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::PrimaryVertex(const TObjArray *trkArray,
								   AliVEvent *event)
{
  //
  /// Used only for pp
  /// copied from AliAnalysisVertexingHF (except for the following 3 lines)
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
AliAODVertex* AliAnalysisTaskSELc2eleLambdafromAODtracks::ReconstructSecondaryVertex(AliAODv0 *v0, AliAODTrack *part, AliAODEvent * aod) 
{
  //
  // Reconstruct secondary vertex from trkArray (Copied from AliAnalysisVertexingHF)
	// Currently only returns Primary vertex (can we reconstruct secondary vertex from e - v0??)
  //
	
  AliAODVertex *primVertexAOD;
  Bool_t unsetvtx = kFALSE;
  if(fReconstructPrimVert){
    primVertexAOD = CallPrimaryVertex(v0,part,aod);
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
Int_t AliAnalysisTaskSELc2eleLambdafromAODtracks::MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_v0, Int_t *labelarray_ele, Int_t *labelarray_v0,  Int_t &ngen_ele, Int_t &ngen_v0) 
{
  //
  // Match to MC
  //
	for(Int_t i=0;i<100;i++){
		pdgarray_ele[i] = -9999;
		labelarray_ele[i] = -9999;
		pdgarray_v0[i] = -9999;
		labelarray_v0[i] = -9999;
	}
	ngen_ele = 0;
	ngen_v0 = 0;

  AliVTrack *trk = dynamic_cast<AliVTrack*>(elobj->GetBachelor());
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

  AliAODv0 *theV0 = dynamic_cast<AliAODv0*>(elobj->Getv0());
	if(!theV0) return -1;
	Int_t pdgdgv0[2]={2212,211};
  Int_t labV0 = theV0->MatchToMC(3122,mcArray,2,pdgdgv0); // the V0
	if(labV0<0) return -1;
	AliAODMCParticle *mcv0 = (AliAODMCParticle*)mcArray->At(labV0);
	if(!mcv0) return -1;
	labelarray_v0[0] = labV0;
	pdgarray_v0[0] = mcv0->GetPdgCode();
	ngen_v0 ++;

  AliAODMCParticle *mcprimv0=0;
  mcprimv0 = mcv0;
  while(mcprimv0->GetMother()>=0) {
    Int_t labprim_v0=mcprimv0->GetMother();
    AliAODMCParticle *tmcprimv0 = (AliAODMCParticle*)mcArray->At(labprim_v0);
    if(!tmcprimv0) {
			break;
    }

    mcprimv0 = tmcprimv0;
		pdgarray_v0[ngen_v0] = mcprimv0->GetPdgCode();
		labelarray_v0[ngen_v0] = labprim_v0;
		ngen_v0 ++;
		if(ngen_v0==100) break;
  }

	Bool_t same_flag = kFALSE;
	Int_t matchedlabel=-9999;
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_v0;ivmc++){
			if(labelarray_ele[iemc]==labelarray_v0[ivmc]){
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray)
{
  //
  // Select good tracks using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //
  
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
		Double_t nsigma_itsele = -9999;
		if(fAnalCuts->GetIsUsePID()){
			nsigma_tpcele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(aodtpid,AliPID::kElectron);
			nsigma_tofele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(aodtpid,AliPID::kElectron);
			nsigma_itsele = fAnalCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasITS(aodtpid,AliPID::kElectron);
		}

    if(fAnalCuts->SingleTrkCutsNoPID(aodt,aodtpid,fVtx1)){
			fHistoElectronTPCPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFPID->Fill(aodt->Pt(),nsigma_tofele);
			fHistoElectronITSPID->Fill(aodt->Pt(),nsigma_itsele);
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

				if(nsigma_tpcele>-0.5 && nsigma_tpcele<3.){
					fHistoElectronITSPIDSelTPCTOF->Fill(aodt->Pt(),nsigma_itsele);
					if(fabs(eleeta)<0.6)
						fHistoElectronITSPIDSelTPCTOFSmallEta->Fill(aodt->Pt(),nsigma_itsele);
					
					if(fabs(eleeta)>0.6 && fabs(eleeta)<0.8)
						fHistoElectronITSPIDSelTPCTOFLargeEta->Fill(aodt->Pt(),nsigma_itsele);
					
					if(eleeta>-0.8 && eleeta<-0.6){
						fHistoElectronITSPIDSelTPCTOFEtaDep[0]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>-0.6&&eleeta<-0.4){
						fHistoElectronITSPIDSelTPCTOFEtaDep[1]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>-0.4&&eleeta<-0.2){
						fHistoElectronITSPIDSelTPCTOFEtaDep[2]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>-0.2&&eleeta<0.0){
						fHistoElectronITSPIDSelTPCTOFEtaDep[3]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>0.0&&eleeta<0.2){
						fHistoElectronITSPIDSelTPCTOFEtaDep[4]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>0.2&&eleeta<0.4){
						fHistoElectronITSPIDSelTPCTOFEtaDep[5]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>0.4&&eleeta<0.6){
						fHistoElectronITSPIDSelTPCTOFEtaDep[6]->Fill(aodt->Pt(),nsigma_itsele);
					}else if(eleeta>0.6&&eleeta<0.8){
						fHistoElectronITSPIDSelTPCTOFEtaDep[7]->Fill(aodt->Pt(),nsigma_itsele);
					}
				}
			}
		}

    if(fAnalCuts->SingleTrkCuts(aodt,aodtpid,fVtx1)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
			fHistoElectronTPCSelPID->Fill(aodt->Pt(),nsigma_tpcele);
			fHistoElectronTOFSelPID->Fill(aodt->Pt(),nsigma_tofele);
			fHistoElectronITSSelPID->Fill(aodt->Pt(),nsigma_itsele);
    }
  } // end loop on tracks
}
//________________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::SelectV0( const AliVEvent *event,Int_t nV0s,Int_t &nSeleV0, Bool_t *seleV0Flags, TClonesArray *mcArray)
{
  //
  // Select good V0 using fAnalCuts (AliRDHFCuts object) and return the array of their ids
  //

  nSeleV0 = 0;
  for(Int_t iv0=0;iv0<nV0s;iv0++)
    {
      seleV0Flags[iv0] = kFALSE;
      AliAODv0 *v0 = ((AliAODEvent*)event)->GetV0(iv0);

      if(!fAnalCuts) continue;
      if(fAnalCuts->SingleV0Cuts(v0,fVtx1)){
				seleV0Flags[iv0] = kTRUE;
				nSeleV0++;
      }
    }
}
//_________________________________________________________________
Int_t AliAnalysisTaskSELc2eleLambdafromAODtracks::GetPoolIndex(Double_t zvert, Double_t mult, Double_t rp){
	//
  // check in which of the pools the current event falls
	//

  Int_t theBinZ=TMath::BinarySearch(fNzVtxBins,fZvtxBins,zvert);
  if(theBinZ<0 || theBinZ>=fNzVtxBins) return -1;
  Int_t theBinM=TMath::BinarySearch(fNCentBins,fCentBins,mult);
  if(theBinM<0 || theBinM>=fNCentBins) return -2;
  //return fNCentBins*theBinZ+theBinM;
  Int_t theBinR=TMath::BinarySearch(fNRPBins,fRPBins,rp);
  if(theBinR<0 || theBinR>=fNRPBins) return -3;
  return fNzVtxBins*fNCentBins*theBinR+fNCentBins*theBinZ+theBinM;
}
//_________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DoEventMixingWithPools(Int_t poolIndex)
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
				fHistoPoolMonitor->Fill(poolIndex,nextRes);
      }
    }
  }
}

//_________________________________________________________________
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillBackground(std::vector<TLorentzVector * > mixTypeE,std::vector<TVector * > mixTypeEVars, std::vector<TLorentzVector * > mixTypeL, std::vector<TVector * > mixTypeLVars, Int_t charge_v0pr)
{     
  //
  // Fill background
  //
  int nEle = mixTypeE.size();
  int nV0 = mixTypeL.size();
  for(Int_t ie=0;ie<nEle;ie++){
    TLorentzVector* trke=mixTypeE[ie];
    if(!trke) continue;
    TVector *elevars = mixTypeEVars[ie];
    for(Int_t iv=0;iv<nV0;iv++){
      TLorentzVector* v0=mixTypeL[iv];
      TVector *v0vars = mixTypeLVars[iv];
      if(!v0) continue;
      FillMixROOTObjects(trke,v0,elevars,v0vars,charge_v0pr);
    }
  }
  return;
}
//_________________________________________________________________
Bool_t AliAnalysisTaskSELc2eleLambdafromAODtracks::MakeMCAnalysis(TClonesArray *mcArray)
{
	//
  // Analyze AliAODmcparticle
	//

	Int_t nmcpart = mcArray->GetEntriesFast();

	Int_t mcevttype = 0;
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
				}
			}
			if(TMath::Abs(mcpart->GetPdgCode())==5){
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
	}

	for(Int_t i=0;i<nmcpart;i++)
	{
		AliAODMCParticle *mcpart = (AliAODMCParticle*) mcArray->At(i);
		if(TMath::Abs(mcpart->GetPdgCode())==4122){
			//cout<<"Lambdac"<<endl;
			Bool_t e_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			Bool_t sigma_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3122){
					lam_flag = kTRUE;
					mcv0part = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3212){
					sigma_flag = kTRUE;
					mcv0part = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3214){
					sigma_flag = kTRUE;
					mcv0part = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3224){
					sigma_flag = kTRUE;
					mcv0part = mcdau;
				}
			}

			Int_t decaytype = -9999;
			if(e_flag && lam_flag) decaytype = 0;
			if(e_flag && sigma_flag) decaytype = 3;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==4132){
			//cout<<"Lambdac"<<endl;
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
			AliAODMCParticle *mcv0part = 0;
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
					for(Int_t idauxi=mccascpart->GetFirstDaughter();idauxi<mccascpart->GetLastDaughter()+1;idauxi++)
					{
						if(idauxi<0) break;
						AliAODMCParticle *mcdauxi = (AliAODMCParticle*) mcArray->At(idauxi);
						if(!mcdauxi) continue;
						if(TMath::Abs(mcdauxi->GetPdgCode())==3122){
							lam_flag = kTRUE;
							mcv0part = mcdauxi;
						}
					}
				}
			}
			Int_t decaytype = -9999;
			if(e_flag && xi_flag && lam_flag) decaytype = 1;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}
		if(TMath::Abs(mcpart->GetPdgCode())==4232){
			//cout<<"Lambdac"<<endl;
			Bool_t e_flag = kFALSE;
			Bool_t xi_flag = kFALSE;
			Bool_t lam_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mccascpart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==3322){
					xi_flag = kTRUE;
					mccascpart = mcdau;
					for(Int_t idauxi=mccascpart->GetFirstDaughter();idauxi<mccascpart->GetLastDaughter()+1;idauxi++)
					{
						if(idauxi<0) break;
						AliAODMCParticle *mcdauxi = (AliAODMCParticle*) mcArray->At(idauxi);
						if(!mcdauxi) continue;
						if(TMath::Abs(mcdauxi->GetPdgCode())==3122){
							lam_flag = kTRUE;
							mcv0part = mcdauxi;
						}
					}
				}
			}
			Int_t decaytype = -9999;
			if(e_flag && xi_flag && lam_flag) decaytype = 2;

			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}

		if(TMath::Abs(mcpart->GetPdgCode())==5122){
			Bool_t e_flag = kFALSE;
			Bool_t lc_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mclcpart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==4122){
					lc_flag = kTRUE;
					mclcpart = mcdau;
				}
			}

			Bool_t lam_flag = kFALSE;
      if(e_flag && lc_flag){
        for(Int_t idau=mclcpart->GetFirstDaughter();idau<mclcpart->GetLastDaughter()+1;idau++)
        {
          if(idau<0) break;
          AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
          if(!mcdau) continue;
          if(TMath::Abs(mcdau->GetPdgCode())==3122){
            lam_flag = kTRUE;
            mcv0part = mcdau;
          }
        }
      }

			if(lc_flag){
				Double_t contmc_withlc[3];
				contmc_withlc[0] = mcpart->Pt();
				contmc_withlc[1] = mcpart->Y();
				contmc_withlc[2] = mclcpart->Y();
				if(fabs(mclcpart->Y())<1.){
					fHistoResponseMcGenLambdabPtvsLambdacPt->Fill(mcpart->Pt(),mclcpart->Pt());
				}
				fHistoLbMCGenWithLambdac->Fill(contmc_withlc);
			}

			Int_t decaytype = -9999;
			if(e_flag && lc_flag && lam_flag) decaytype = 10;
			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
		}

		if(TMath::Abs(mcpart->GetPdgCode())==5132 || TMath::Abs(mcpart->GetPdgCode())==5232){
			Bool_t e_flag = kFALSE;
			Bool_t xc_flag = kFALSE;
			AliAODMCParticle *mcepart = 0;
			AliAODMCParticle *mclcpart = 0;
			AliAODMCParticle *mcv0part = 0;
			for(Int_t idau=mcpart->GetFirstDaughter();idau<mcpart->GetLastDaughter()+1;idau++)
			{
				if(idau<0) break;
				AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
				if(!mcdau) continue;
				if(TMath::Abs(mcdau->GetPdgCode())==11){
					e_flag = kTRUE;
					mcepart = mcdau;
				}
				if(TMath::Abs(mcdau->GetPdgCode())==4132 || TMath::Abs(mcdau->GetPdgCode())==4232){
					xc_flag = kTRUE;
					mclcpart = mcdau;
				}
			}

			Bool_t lam_flag = kFALSE;
      if(e_flag && xc_flag){
        for(Int_t idau=mclcpart->GetFirstDaughter();idau<mclcpart->GetLastDaughter()+1;idau++)
        {
          if(idau<0) break;
          AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
          if(!mcdau) continue;
          if(TMath::Abs(mcdau->GetPdgCode())==3312 || TMath::Abs(mcdau->GetPdgCode())==3322){
            lam_flag = kTRUE;
            mcv0part = mcdau;
          }
        }
      }
			Int_t decaytype = -9999;
			if(e_flag && xc_flag && lam_flag) decaytype = 11;
			FillMCROOTObjects(mcpart,mcepart,mcv0part,decaytype);
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
		if(TMath::Abs(mcpart->GetPdgCode())==3122){
			Double_t etamin, etamax, rapmin, rapmax;
			fAnalCuts->GetProdV0EtaRange(etamin,etamax);
			fAnalCuts->GetProdV0RapRange(rapmin,rapmax);

			if((fabs(mcpart->Y())<rapmax) && (fabs(mcpart->Eta())<etamax)){
				fHistoLambdaMassvsPtMCGen->Fill(1.115683, mcpart->Pt());
				fHistoLambdaPtvsEtaMCGen->Fill(mcpart->Pt(),mcpart->Eta());
			}
			FillMCV0ROOTObjects(mcpart, mcArray);
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
				if(TMath::Abs(mcpartv->GetPdgCode())!=3122) continue;
				if(mcpartv->Pt()<0.4) continue;//Apply rough cuts
				if(fabs(mcpartv->Eta())>0.8) continue;//Apply rough cuts
				if(mcpartv->GetNDaughters()!=2) continue;

				FillMCGenPairROOTObjects(mcparte,mcpartv,mcArray);
			}
		}
		return kFALSE;
	}

	return kTRUE;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineMCGenPairTreeVariables() 
{
  //
  // Define mc pair tree variables
  //

  const char* nameoutput = GetOutputSlot(11)->GetContainer()->GetName();
  fMCGenPairVariablesTree = new TTree(nameoutput,"MC pair variables tree");
  Int_t nVar = 38;
  fCandidateMCGenPairVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

	fCandidateVariableNames[ 0] = "InvMassEleLambda";
	fCandidateVariableNames[ 1] = "EleLambdaPx";
	fCandidateVariableNames[ 2] = "EleLambdaPy";
	fCandidateVariableNames[ 3] = "EleLambdaPz";
	fCandidateVariableNames[ 4] = "ElePdgCode";
	fCandidateVariableNames[ 5] = "ElePx";
	fCandidateVariableNames[ 6] = "ElePy";
	fCandidateVariableNames[ 7] = "ElePz";
	fCandidateVariableNames[ 8] = "LambdaPdgCode";
	fCandidateVariableNames[ 9] = "LambdaPx";
	fCandidateVariableNames[10] = "LambdaPy";
	fCandidateVariableNames[11] = "LambdaPz";
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
	fCandidateVariableNames[25] = "LamNGeneration";
	fCandidateVariableNames[26] = "LambdaGen1PDG";
	fCandidateVariableNames[27] = "LambdaGen2PDG";
	fCandidateVariableNames[28] = "LambdaGen3PDG";
	fCandidateVariableNames[29] = "LambdaGen4PDG";
	fCandidateVariableNames[30] = "LambdaGen5PDG";
	fCandidateVariableNames[31] = "LambdaGen6PDG";
	fCandidateVariableNames[32] = "LambdaGen7PDG";
	fCandidateVariableNames[33] = "LambdaGen8PDG";
	fCandidateVariableNames[34] = "LambdaGen9PDG";
	fCandidateVariableNames[35] = "LambdaGen10PDG";
	fCandidateVariableNames[36] = "LambdaPrimPDG";
	fCandidateVariableNames[37] = "MatchedPDG";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fMCGenPairVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCandidateMCGenPairVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::FillMCGenPairROOTObjects(AliAODMCParticle *mcparte, AliAODMCParticle *mcpartv, TClonesArray *mcArray) 
{
  //
  // Fill histograms or mc pair analysis tree 
  //
	for(Int_t i=0;i<38;i++){
		fCandidateMCGenPairVariables[i] = -9999.;
	}

	TLorentzVector vele, vlam, velam;
	vele.SetXYZM(mcparte->Px(),mcparte->Py(),mcparte->Pz(),0.000510998928);
	vlam.SetXYZM(mcpartv->Px(),mcpartv->Py(),mcpartv->Pz(),1.115683);
	velam = vele + vlam;

	Int_t pdgarray_ele[100], labelarray_ele[100], ngen_ele;
	Int_t pdgarray_lam[100], labelarray_lam[100], ngen_lam;
	GetMCDecayHistory(mcparte,mcArray,pdgarray_ele,labelarray_ele,ngen_ele);
	GetMCDecayHistory(mcpartv,mcArray,pdgarray_lam,labelarray_lam,ngen_lam);

	Bool_t same_flag = kFALSE;
	Int_t matched_pdg = -999999;
	for(Int_t iemc=0;iemc<ngen_ele;iemc++){
		for(Int_t ivmc=0;ivmc<ngen_lam;ivmc++){
			if(labelarray_ele[iemc]==labelarray_lam[ivmc]){
				same_flag = kTRUE;
				matched_pdg = pdgarray_ele[iemc];
				break;
			}
		}
		if(same_flag) break;
	}
	Int_t pdgprim_ele = pdgarray_ele[ngen_ele-1];
	Int_t pdgprim_lam = pdgarray_lam[ngen_lam-1];

	fCandidateMCGenPairVariables[ 0] = velam.M();
	fCandidateMCGenPairVariables[ 1] = velam.Px();
	fCandidateMCGenPairVariables[ 2] = velam.Py();
	fCandidateMCGenPairVariables[ 3] = velam.Pz();
	fCandidateMCGenPairVariables[ 4] = mcparte->GetPdgCode();
	fCandidateMCGenPairVariables[ 5] = vele.Px();
	fCandidateMCGenPairVariables[ 6] = vele.Py();
	fCandidateMCGenPairVariables[ 7] = vele.Pz();
	fCandidateMCGenPairVariables[ 8] = mcpartv->GetPdgCode();
	fCandidateMCGenPairVariables[ 9] = vlam.Px();
	fCandidateMCGenPairVariables[10] = vlam.Py();
	fCandidateMCGenPairVariables[11] = vlam.Pz();
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
	fCandidateMCGenPairVariables[25] = (Float_t)ngen_lam;
	fCandidateMCGenPairVariables[26] = (Float_t)pdgarray_lam[0];
	fCandidateMCGenPairVariables[27] = (Float_t)pdgarray_lam[1];
	fCandidateMCGenPairVariables[28] = (Float_t)pdgarray_lam[2];
	fCandidateMCGenPairVariables[29] = (Float_t)pdgarray_lam[3];
	fCandidateMCGenPairVariables[30] = (Float_t)pdgarray_lam[4];
	fCandidateMCGenPairVariables[31] = (Float_t)pdgarray_lam[5];
	fCandidateMCGenPairVariables[32] = (Float_t)pdgarray_lam[6];
	fCandidateMCGenPairVariables[33] = (Float_t)pdgarray_lam[7];
	fCandidateMCGenPairVariables[34] = (Float_t)pdgarray_lam[8];
	fCandidateMCGenPairVariables[35] = (Float_t)pdgarray_lam[9];
	fCandidateMCGenPairVariables[36] = (Float_t)pdgarray_lam[ngen_lam-1];
	fCandidateMCGenPairVariables[37] = (Float_t) matched_pdg;

	fMCGenPairVariablesTree->Fill();
}

////-------------------------------------------------------------------------------
void AliAnalysisTaskSELc2eleLambdafromAODtracks::DefineCorrelationTreeVariables() 
{
  //
  // Define mc pair tree variables
  //

  const char* nameoutput = GetOutputSlot(12)->GetContainer()->GetName();
  fCorrelationVariablesTree = new TTree(nameoutput,"Correlation variables tree");
  Int_t nVar = 15;
  fCorrelationVariables = new Float_t [nVar];
  TString * fCandidateVariableNames = new TString[nVar];

  fCandidateVariableNames[0] = "V0Pt";
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
  fCandidateVariableNames[11] = "EleLambdaPt";
  fCandidateVariableNames[12] = "EleLambdaMass";
  fCandidateVariableNames[13] = "EleLambdaCosOA";
  fCandidateVariableNames[14] = "MCEleMother";


  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fCorrelationVariablesTree->Branch(fCandidateVariableNames[ivar].Data(),&fCorrelationVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
  }
  return;
}


////-------------------------------------------------------------------------------
void	AliAnalysisTaskSELc2eleLambdafromAODtracks::GetMCDecayHistory(AliAODMCParticle *mcpart, TClonesArray *mcArray, Int_t *pdgarray, Int_t *labelarray, Int_t &ngen)
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::StoreGlobalTrackReference(AliAODTrack *track, Int_t index){
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
void AliAnalysisTaskSELc2eleLambdafromAODtracks::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
    fGTIndex[i]=-9999;
  }
}
//________________________________________________________________________
Int_t AliAnalysisTaskSELc2eleLambdafromAODtracks::FromSemileptonicDecays(Int_t *history){
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
Bool_t AliAnalysisTaskSELc2eleLambdafromAODtracks::HaveCharmInHistory(Int_t *history){
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
Bool_t AliAnalysisTaskSELc2eleLambdafromAODtracks::HaveBottomInHistory(Int_t *history){
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
TProfile* AliAnalysisTaskSELc2eleLambdafromAODtracks::GetEstimatorHistogram(const AliVEvent* event){
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
Float_t AliAnalysisTaskSELc2eleLambdafromAODtracks::GetEventPlaneForCandidate(AliAODTrack* trk, AliAODv0 *v0,AliEventplane *pl, TVector2* q,TVector2* qsub1,TVector2* qsub2){
	//
  // remove autocorrelations 
	//
 
  TArrayF* qx = pl->GetQContributionXArray();
  TArrayF* qy = pl->GetQContributionYArray();
  TVector2 qcopy = *q;

	if(trk && v0){
    AliAODTrack *track0 = (AliAODTrack*)v0->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)v0->GetDaughter(1);  

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}

    TVector2 q2;
    if((trk->GetID()) < qx->fN){
      q2.Set(qx->At(trk->GetID()),qy->At(trk->GetID()));}
      
    qcopy = qcopy -(q0+q1+q2);
	}else if(trk){
    TVector2 q0;
    if((trk->GetID()) < qx->fN){
      q0.Set(qx->At(trk->GetID()),qy->At(trk->GetID()));}
      
    qcopy = qcopy -q0;
	}else if(v0){
    AliAODTrack *track0 = (AliAODTrack*)v0->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)v0->GetDaughter(1);  

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
      
    qcopy = qcopy -(q0+q1);
	}
  
  return qcopy.Phi()/2.;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSELc2eleLambdafromAODtracks::GetPhi0Pi(Float_t phi){
	//
	// Sets the phi angle in the range 0-pi
	//
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
  return result;
}
