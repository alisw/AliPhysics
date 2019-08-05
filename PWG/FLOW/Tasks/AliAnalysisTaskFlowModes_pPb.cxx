/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// AliAnalysisTaskFlowModes_pPb - ALICE Flow framework
//
// ALICE analysis task for universal study of flow.
// Note: So far implemented only for AOD analysis!
// Author: Naghmeh Mohammadi
// Modified by: Barnabas Porfy, CERN Summer Student, 2019
// Generic framework by: You Zhou

#include <vector>

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TList.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAnalysisTaskFlowModes_pPb.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVTrack.h"
#include "AliESDpid.h"
#include "AliFlowBayesianPID.h"


class AliAnalysisTaskFlowModes_pPb;

ClassImp(AliAnalysisTaskFlowModes_pPb); // classimp: necessary for root

Int_t AliAnalysisTaskFlowModes_pPb::fHarmonics[] = {2,3,4,5,6};
Int_t AliAnalysisTaskFlowModes_pPb::fMixedHarmonics[] = {422,633,523,6222};// 422: v4{psi2}, 523: v5{psi2, psi3}, 633: v6{psi3} and 6222: v6{psi2}

Double_t AliAnalysisTaskFlowModes_pPb::fEtaGap[] = {0.,0.4,0.8};



AliAnalysisTaskFlowModes_pPb::AliAnalysisTaskFlowModes_pPb() : AliAnalysisTaskSE(),
fCheckPileUp(kFALSE),
fUsePileUpSPD(kFALSE),
fModifySPDDefaultParams(kFALSE),
fMinVtxPileUpContrSPD(),
fMinPileUpZdistSPD(),
fRejectOutOfBunchPileUp(kFALSE),

fEventAOD(0x0),
fPIDResponse(0x0),
fPIDCombined(0x0),
fFlowNUAWeightsFile(0x0),
fFlowNUEWeightsFile(0x0),
fInit(kFALSE),
fIndexSampling(0),
fIndexCentrality(-1),
fEventCounter(0),
fNumEventsAnalyse(50),
fRunNumber(-1),
fExtraPileUp(kFALSE),
fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),

// FlowPart containers
fVectorCharged(0x0),
fVectorPion(0x0),
fVectorKaon(0x0),
fVectorProton(0x0),

// analysis selection
fRunMode(kFull),
fAnalType(kAOD),
fSampling(kFALSE),
fFillQA(kTRUE),
fProcessCharged(kFALSE),
fProcessPID(kFALSE),
fBayesianResponse(NULL),

// flow related
fCutFlowRFPsPtMin(0),
fCutFlowRFPsPtMax(0),
fFlowPOIsPtMin(0),
fFlowPOIsPtMax(10.),
fFlowCentMin(0),
fFlowCentMax(150),
fFlowCentNumBins(150),
fCutFlowDoFourCorrelations(kTRUE),
fDoOnlyMixedCorrelations(kFALSE),
fFlowFillWeights(kFALSE),
fFlowUseNUAWeights(kFALSE),
fFlowUseNUEWeights(kFALSE),
fFlowNUAWeightsPath(),
fFlowNUEWeightsPath(),
fPositivelyChargedRef(kFALSE),
fNegativelyChargedRef(kFALSE),
fPositivelyChargedPOI(kFALSE),
fNegativelyChargedPOI(kFALSE),



// events selection
fPVtxCutZ(0.),
fColSystem(kPPb),
//fColSystem(kPbPb),
fMultEstimator(),
fTrigger(0),
fFullCentralityRange(kTRUE),
// charged tracks selection
fCutChargedEtaMax(0),
fCutChargedPtMax(0),
fCutChargedPtMin(0),
fCutChargedDCAzMax(0),
fCutChargedDCAxyMax(0),
fCutChargedTrackFilterBit(0),
fCutChargedNumTPCclsMin(0),
fMaxChi2perTPCcls(0),

// PID tracks selection
fESDpid(),
fCutPIDUseAntiProtonOnly(kFALSE),
fPIDnsigma(kTRUE),
fPIDnsigmaCombination(2),
fPIDbayesian(kFALSE),
fCutPIDnSigmaPionMax(3),
fCutPIDnSigmaKaonMax(3),
fCutPIDnSigmaProtonMax(3),
fCutPIDnSigmaTPCRejectElectron(3),
fCutPIDnSigmaCombinedNoTOFrejection(kFALSE),
fCurrCentr(0.0),
fParticleProbability(0.9),

// output lists
fQAEvents(0x0),
fQACharged(0x0),
fQAPID(0x0),
fFlowWeights(0x0),
fFlowRefs(0x0),
fFlowCharged(0x0),
fFlowPID(0x0),

// flow histograms & profiles
fh3BeforeNUAWeightsRefs(0x0),
fh3BeforeNUAWeightsCharged(0x0),
fh3BeforeNUAWeightsPion(0x0),
fh3BeforeNUAWeightsKaon(0x0),
fh3BeforeNUAWeightsProton(0x0),

fh3AfterNUAWeightsRefs(0x0),
fh3AfterNUAWeightsCharged(0x0),
fh3AfterNUAWeightsPion(0x0),
fh3AfterNUAWeightsKaon(0x0),
fh3AfterNUAWeightsProton(0x0),

fhBeforeNUEWeightsRefs(0x0),
fhBeforeNUEWeightsCharged(0x0),
fhBeforeNUEWeightsPion(0x0),
fhBeforeNUEWeightsKaon(0x0),
fhBeforeNUEWeightsProton(0x0),

fhAfterNUEWeightsRefs(0x0),
fhAfterNUEWeightsCharged(0x0),
fhAfterNUEWeightsPion(0x0),
fhAfterNUEWeightsKaon(0x0),
fhAfterNUEWeightsProton(0x0),

fh3NUAWeightRefsPlus(0x0),
fh3NUAWeightChargedPlus(0x0),
fh3NUAWeightPionPlus(0x0),
fh3NUAWeightKaonPlus(0x0),
fh3NUAWeightProtonPlus(0x0),

fh3NUAWeightRefsMinus(0x0),
fh3NUAWeightChargedMinus(0x0),
fh3NUAWeightPionMinus(0x0),
fh3NUAWeightKaonMinus(0x0),
fh3NUAWeightProtonMinus(0x0),

// event histograms
fhEventSampling(0x0),
fhEventCentrality(0x0),
fh2EventCentralityNumSelCharged(0x0),
fhEventCounter(0x0),

// charged histogram
fh2RefsMult(0x0),
fh2RefsPt(0x0),
fh2RefsEta(0x0),
fh2RefsPhi(0x0),
fhChargedCounter(0x0),

// PID histogram
fh2PIDPionMult(0x0),
fh2PIDPionPt(0x0),
fh2PIDPionPhi(0x0),
fh2PIDPionEta(0x0),
fhPIDPionCharge(0x0),
fh2PIDKaonMult(0x0),
fh2PIDKaonPt(0x0),
fh2PIDKaonPhi(0x0),
fh2PIDKaonEta(0x0),
fhPIDKaonCharge(0x0),
fh2PIDProtonMult(0x0),
fh2PIDProtonPt(0x0),
fh2PIDProtonPhi(0x0),
fh2PIDProtonEta(0x0),
fhPIDProtonCharge(0x0),
fh2PIDPionTPCdEdx(0x0),
fh2PIDPionTOFbeta(0x0),
fh2PIDKaonTPCdEdx(0x0),
fh2PIDKaonTOFbeta(0x0),
fh2PIDProtonTPCdEdx(0x0),
fh2PIDProtonTOFbeta(0x0),
fh2PIDPionTPCnSigmaPion(0x0),
fh2PIDPionTOFnSigmaPion(0x0),
fh2PIDPionTPCnSigmaKaon(0x0),
fh2PIDPionTOFnSigmaKaon(0x0),
fh2PIDPionTPCnSigmaProton(0x0),
fh2PIDPionTOFnSigmaProton(0x0),
fh2PIDKaonTPCnSigmaPion(0x0),
fh2PIDKaonTOFnSigmaPion(0x0),
fh2PIDKaonTPCnSigmaKaon(0x0),
fh2PIDKaonTOFnSigmaKaon(0x0),
fh2PIDKaonTPCnSigmaProton(0x0),
fh2PIDKaonTOFnSigmaProton(0x0),
fh2PIDProtonTPCnSigmaPion(0x0),
fh2PIDProtonTOFnSigmaPion(0x0),
fh2PIDProtonTPCnSigmaKaon(0x0),
fh2PIDProtonTOFnSigmaKaon(0x0),
fh2PIDProtonTPCnSigmaProton(0x0),
fh2PIDProtonTOFnSigmaProton(0x0),
fhNUEWeightRefsPlus(0x0),
fhNUEWeightRefsMinus(0x0),
fhNUEWeightChargedPlus(0x0),
fhNUEWeightChargedMinus(0x0),
fhNUEWeightPionPlus(0x0),
fhNUEWeightPionMinus(0x0),
fhNUEWeightKaonPlus(0x0),
fhNUEWeightKaonMinus(0x0),
fhNUEWeightProtonPlus(0x0),
fhNUEWeightProtonMinus(0x0)
{
    SetPriors(); //init arrays
    // New PID procedure (Bayesian Combined PID)
    // allocating here is necessary because we don't
    // stream this member
    fBayesianResponse = new AliFlowBayesianPID();
    fBayesianResponse->SetNewTrackParam();
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskFlowModes_pPb::AliAnalysisTaskFlowModes_pPb(const char* name) : AliAnalysisTaskSE(name),
    //pPb
fCheckPileUp(kFALSE),
fUsePileUpSPD(kFALSE),
fModifySPDDefaultParams(kFALSE),
fMinVtxPileUpContrSPD(),
fMinPileUpZdistSPD(),
fRejectOutOfBunchPileUp(kFALSE),

fEventAOD(0x0),
fPIDResponse(0x0),
fPIDCombined(0x0),
fFlowNUAWeightsFile(0x0),
fFlowNUEWeightsFile(0x0),
fInit(kFALSE),
fIndexSampling(0),
fIndexCentrality(-1),
fEventCounter(0),
fNumEventsAnalyse(50),
fRunNumber(-1),
fExtraPileUp(kFALSE),
fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),

// FlowPart containers
fVectorCharged(0x0),
fVectorPion(0x0),
fVectorKaon(0x0),
fVectorProton(0x0),

// analysis selection
fRunMode(kFull),
fAnalType(kAOD),
fSampling(kFALSE),
fFillQA(kTRUE),
fProcessCharged(kFALSE),
fProcessPID(kFALSE),
fBayesianResponse(NULL),

// flow related
fCutFlowRFPsPtMin(0),
fCutFlowRFPsPtMax(0),
fFlowPOIsPtMin(0),
fFlowPOIsPtMax(10.),
fFlowCentMin(0),
fFlowCentMax(150),
fFlowCentNumBins(150),
fCutFlowDoFourCorrelations(kTRUE),
fDoOnlyMixedCorrelations(kFALSE),
fFlowFillWeights(kFALSE),
fFlowUseNUAWeights(kFALSE),
fFlowUseNUEWeights(kFALSE),
fFlowNUAWeightsPath(),
fFlowNUEWeightsPath(),
fPositivelyChargedRef(kFALSE),
fNegativelyChargedRef(kFALSE),
fPositivelyChargedPOI(kFALSE),
fNegativelyChargedPOI(kFALSE),

// events selection
fPVtxCutZ(0.),
fColSystem(kPPb),
//fColSystem(kPbPb),
fMultEstimator(),
fTrigger(0),
//TRIGGER ERRORS TODO
fFullCentralityRange(kTRUE),
// charged tracks selection
fCutChargedEtaMax(0),
fCutChargedPtMax(0),
fCutChargedPtMin(0),
fCutChargedDCAzMax(0),
fCutChargedDCAxyMax(0),
fCutChargedTrackFilterBit(0),
fCutChargedNumTPCclsMin(0),
fMaxChi2perTPCcls(0),
// PID tracks selection
fESDpid(),
fCutPIDUseAntiProtonOnly(kFALSE),
fPIDnsigma(kTRUE),
fPIDnsigmaCombination(2.),
fPIDbayesian(kFALSE),
fCutPIDnSigmaPionMax(3),
fCutPIDnSigmaKaonMax(3),
fCutPIDnSigmaProtonMax(3),
fCutPIDnSigmaTPCRejectElectron(3),
fCutPIDnSigmaCombinedNoTOFrejection(kFALSE),
fCurrCentr(0.0),
fParticleProbability(0.9),


// output lists
fQAEvents(0x0),
fQACharged(0x0),
fQAPID(0x0),
fFlowWeights(0x0),
fFlowRefs(0x0),
fFlowCharged(0x0),
fFlowPID(0x0),

// flow histograms & profiles
fh3BeforeNUAWeightsRefs(0x0),
fh3BeforeNUAWeightsCharged(0x0),
fh3BeforeNUAWeightsPion(0x0),
fh3BeforeNUAWeightsKaon(0x0),
fh3BeforeNUAWeightsProton(0x0),

fh3AfterNUAWeightsRefs(0x0),
fh3AfterNUAWeightsCharged(0x0),
fh3AfterNUAWeightsPion(0x0),
fh3AfterNUAWeightsKaon(0x0),
fh3AfterNUAWeightsProton(0x0),

fhBeforeNUEWeightsRefs(0x0),
fhBeforeNUEWeightsCharged(0x0),
fhBeforeNUEWeightsPion(0x0),
fhBeforeNUEWeightsKaon(0x0),
fhBeforeNUEWeightsProton(0x0),

fhAfterNUEWeightsRefs(0x0),
fhAfterNUEWeightsCharged(0x0),
fhAfterNUEWeightsPion(0x0),
fhAfterNUEWeightsKaon(0x0),
fhAfterNUEWeightsProton(0x0),

fh3NUAWeightRefsPlus(0x0),
fh3NUAWeightChargedPlus(0x0),
fh3NUAWeightPionPlus(0x0),
fh3NUAWeightKaonPlus(0x0),
fh3NUAWeightProtonPlus(0x0),

fh3NUAWeightRefsMinus(0x0),
fh3NUAWeightChargedMinus(0x0),
fh3NUAWeightPionMinus(0x0),
fh3NUAWeightKaonMinus(0x0),
fh3NUAWeightProtonMinus(0x0),

// event histograms
fhEventSampling(0x0),
fhEventCentrality(0x0),
fh2EventCentralityNumSelCharged(0x0),
fhEventCounter(0x0),

// charged histogram
fh2RefsMult(0x0),
fh2RefsPt(0x0),
fh2RefsEta(0x0),
fh2RefsPhi(0x0),
fhChargedCounter(0x0),

// PID histogram
fh2PIDPionMult(0x0),
fh2PIDPionPt(0x0),
fh2PIDPionPhi(0x0),
fh2PIDPionEta(0x0),
fhPIDPionCharge(0x0),
fh2PIDKaonMult(0x0),
fh2PIDKaonPt(0x0),
fh2PIDKaonPhi(0x0),
fh2PIDKaonEta(0x0),
fhPIDKaonCharge(0x0),
fh2PIDProtonMult(0x0),
fh2PIDProtonPt(0x0),
fh2PIDProtonPhi(0x0),
fh2PIDProtonEta(0x0),
fhPIDProtonCharge(0x0),
fh2PIDPionTPCdEdx(0x0),
fh2PIDPionTOFbeta(0x0),
fh2PIDKaonTPCdEdx(0x0),
fh2PIDKaonTOFbeta(0x0),
fh2PIDProtonTPCdEdx(0x0),
fh2PIDProtonTOFbeta(0x0),
fh2PIDPionTPCnSigmaPion(0x0),
fh2PIDPionTOFnSigmaPion(0x0),
fh2PIDPionTPCnSigmaKaon(0x0),
fh2PIDPionTOFnSigmaKaon(0x0),
fh2PIDPionTPCnSigmaProton(0x0),
fh2PIDPionTOFnSigmaProton(0x0),
fh2PIDKaonTPCnSigmaPion(0x0),
fh2PIDKaonTOFnSigmaPion(0x0),
fh2PIDKaonTPCnSigmaKaon(0x0),
fh2PIDKaonTOFnSigmaKaon(0x0),
fh2PIDKaonTPCnSigmaProton(0x0),
fh2PIDKaonTOFnSigmaProton(0x0),
fh2PIDProtonTPCnSigmaPion(0x0),
fh2PIDProtonTOFnSigmaPion(0x0),
fh2PIDProtonTPCnSigmaKaon(0x0),
fh2PIDProtonTOFnSigmaKaon(0x0),
fh2PIDProtonTPCnSigmaProton(0x0),
fh2PIDProtonTOFnSigmaProton(0x0),
fhNUEWeightRefsPlus(0x0),
fhNUEWeightRefsMinus(0x0),
fhNUEWeightChargedPlus(0x0),
fhNUEWeightChargedMinus(0x0),
fhNUEWeightPionPlus(0x0),
fhNUEWeightPionMinus(0x0),
fhNUEWeightKaonPlus(0x0),
fhNUEWeightKaonMinus(0x0),
fhNUEWeightProtonPlus(0x0),
fhNUEWeightProtonMinus(0x0)


{
    
    SetPriors(); //init arrays
    // New PID procedure (Bayesian Combined PID)
    fBayesianResponse = new AliFlowBayesianPID();
    fBayesianResponse->SetNewTrackParam();
    
    // Flow vectors
    for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    {
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
            fFlowVecQpos[iHarm][iPower] = TComplex(0,0,kFALSE);
            fFlowVecQneg[iHarm][iPower] = TComplex(0,0,kFALSE);
            
            for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
            {
                fFlowVecPpos[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
                fFlowVecPneg[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
                fFlowVecS[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
            }//endfor(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
        }//endfor(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    }//endfor(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    
    // Flow profiles & histograms
    if(fDoOnlyMixedCorrelations){
        for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics; iMixedHarm++)
        {
            for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
            {
                for(Short_t iSample(0); iSample < fNumSamples; iSample++)
                {
                    fpMixedRefsCor4[iSample][iGap][iMixedHarm] = 0x0;
                    if(iMixedHarm==3)fpMixedRefsCor6[iSample][iGap] = 0x0;
                    
                    fpMixedChargedCor3Pos[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedChargedCor3Neg[iSample][iGap][iMixedHarm] = 0x0;
                    
                    fpMixedPionCor3Pos[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedPionCor3Neg[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedKaonCor3Pos[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedKaonCor3Neg[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedProtonCor3Pos[iSample][iGap][iMixedHarm] = 0x0;
                    fpMixedProtonCor3Neg[iSample][iGap][iMixedHarm] = 0x0;
                    
                    if(iMixedHarm==3){
                        fpMixedChargedCor4Pos[iSample][iGap] = 0x0;
                        fpMixedChargedCor4Neg[iSample][iGap] = 0x0;
                        
                        fpMixedPionCor4Pos[iSample][iGap] = 0x0;
                        fpMixedPionCor4Neg[iSample][iGap] = 0x0;
                        fpMixedKaonCor4Pos[iSample][iGap] = 0x0;
                        fpMixedKaonCor4Neg[iSample][iGap] = 0x0;
                        fpMixedProtonCor4Pos[iSample][iGap] = 0x0;
                        fpMixedProtonCor4Neg[iSample][iGap] = 0x0;
                    }
                }
                
            }//endfor(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        }//endfor(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics; iMixedHarm++)
    }//endif(fDoOnlyMixedCorrelations)
    
    if(!fDoOnlyMixedCorrelations){
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++){
            //fpRefsCor4[iHarm] = 0x0;
            //fp2ChargedCor4[iHarm] = 0x0;
            //fp2PionCor4[iHarm] = 0x0;
            //fp2KaonCor4[iHarm] = 0x0;
            //fp2ProtonCor4[iHarm] = 0x0;
            for(Short_t iGap(0); iGap < fNumEtaGap; iGap++){
                // mean Qx,Qy
                fpMeanQxRefsPos[iGap][iHarm] = 0x0;
                fpMeanQxRefsNeg[iGap][iHarm] = 0x0;
                fpMeanQyRefsPos[iGap][iHarm] = 0x0;
                fpMeanQyRefsNeg[iGap][iHarm] = 0x0;
                for(Short_t iSample(0); iSample < fNumSamples; iSample++)
                {
                    
                    fpRefsCor2[iSample][iGap][iHarm] = 0x0;
                    fp2ChargedCor2Pos[iSample][iGap][iHarm] = 0x0;
                    fp2ChargedCor2Neg[iSample][iGap][iHarm] = 0x0;
                    fp2PionCor2Pos[iSample][iGap][iHarm] = 0x0;
                    fp2PionCor2Neg[iSample][iGap][iHarm] = 0x0;
                    fp2KaonCor2Pos[iSample][iGap][iHarm] = 0x0;
                    fp2KaonCor2Neg[iSample][iGap][iHarm] = 0x0;
                    fp2ProtonCor2Pos[iSample][iGap][iHarm] = 0x0;
                    fp2ProtonCor2Neg[iSample][iGap][iHarm] = 0x0;
                }
            }//endfor(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        }//endfor(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
    }//endif(!fDoOnlyMixedCorrelations)
    
    // QA histograms
    for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
    {
        // Event histograms
        fhQAEventsPVz[iQA] = 0x0;
        fhQAEventsNumContrPV[iQA] = 0x0;
        fhQAEventsNumSPDContrPV[iQA] = 0x0;
        fhQAEventsDistPVSPD[iQA] = 0x0;
        fhQAEventsSPDresol[iQA] = 0x0;
        fhQAEventsPileUp[iQA] = 0x0;
        fhQAEventsCentralityOutliers[iQA] = 0x0;
        fhEventsMultTOFFilterbit32[iQA] = 0x0;
        // charged
        fhQAChargedMult[iQA] = 0x0;
        fhQAChargedPt[iQA] = 0x0;
        fhQAChargedEta[iQA] = 0x0;
        fhQAChargedPhi[iQA] = 0x0;
        fhQAChargedCharge[iQA] = 0x0;
        fhQAChargedFilterBit[iQA] = 0x0;
        fhQAChargedNumTPCcls[iQA] = 0x0;
        fhQAChargedDCAxy[iQA] = 0x0;
        fhQAChargedDCAz[iQA] = 0x0;
        
        // PID
        fhQAPIDTPCstatus[iQA] = 0x0;
        fhQAPIDTOFstatus[iQA] = 0x0;
        fhQAPIDTPCdEdx[iQA] = 0x0;
        fhQAPIDTOFbeta[iQA] = 0x0;
        
        fh3PIDPionTPCTOFnSigmaPion[iQA] = 0x0;
        fh3PIDPionTPCTOFnSigmaKaon[iQA] = 0x0;
        fh3PIDPionTPCTOFnSigmaProton[iQA] = 0x0;
        
        fh3PIDKaonTPCTOFnSigmaPion[iQA] = 0x0;
        fh3PIDKaonTPCTOFnSigmaKaon[iQA] = 0x0;
        fh3PIDKaonTPCTOFnSigmaProton[iQA] = 0x0;
        
        fh3PIDProtonTPCTOFnSigmaPion[iQA] = 0x0;
        fh3PIDProtonTPCTOFnSigmaKaon[iQA] = 0x0;
        fh3PIDProtonTPCTOFnSigmaProton[iQA] = 0x0;
        
        
    }//endfor(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
    
    // defining input/output
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());
    DefineOutput(4, TList::Class());
    DefineOutput(5, TList::Class());
    DefineOutput(6, TList::Class());
    DefineOutput(7, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskFlowModes_pPb::~AliAnalysisTaskFlowModes_pPb()
{
    // destructor
    // if(fPIDCombined)
    // {
    //   delete fPIDCombined;
    // }
    
    // deleting FlowPart vectors (containers)
    if(fVectorCharged) delete fVectorCharged;
    if(fVectorPion) delete fVectorPion;
    if(fVectorKaon) delete fVectorKaon;
    if(fVectorProton) delete fVectorProton;
    
    // deleting output lists
    if(fFlowWeights) delete fFlowWeights;
    if(fFlowRefs) delete fFlowRefs;
    if(fFlowCharged) delete fFlowCharged;
    if(fFlowPID) delete fFlowPID;
    
    if(fQAEvents) delete fQAEvents;
    if(fQACharged) delete fQACharged;
    if(fQAPID) delete fQAPID;
    //deleting histograms
    if(fhNUEWeightRefsPlus)     delete fhNUEWeightRefsPlus;
    if(fhNUEWeightRefsMinus)    delete fhNUEWeightRefsMinus;
    if(fhNUEWeightChargedPlus)  delete fhNUEWeightChargedPlus;
    if(fhNUEWeightChargedMinus) delete fhNUEWeightChargedMinus;
    if(fhNUEWeightPionPlus)     delete fhNUEWeightPionPlus;
    if(fhNUEWeightPionMinus)    delete fhNUEWeightPionMinus;
    if(fhNUEWeightKaonPlus)     delete fhNUEWeightKaonPlus;
    if(fhNUEWeightKaonMinus)    delete fhNUEWeightKaonMinus;
    if(fhNUEWeightProtonPlus)   delete fhNUEWeightProtonPlus;
    if(fhNUEWeightProtonMinus)  delete fhNUEWeightProtonMinus;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::UserCreateOutputObjects()
{
    // create output objects
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // *************************************************************
    
    // list all parameters used in this analysis
    ListParameters();
    
    // task initialization
    fInit = InitializeTask();
    if(!fInit) return;
    
    // creating output lists
    fFlowRefs = new TList();
    fFlowRefs->SetOwner(kTRUE);
    fFlowRefs->SetName("fFlowRefs");
    fFlowCharged = new TList();
    fFlowCharged->SetOwner(kTRUE);
    fFlowCharged->SetName("fFlowCharged");
    fFlowPID = new TList();
    fFlowPID->SetOwner(kTRUE);
    fFlowPID->SetName("fFlowPID");
    fFlowWeights = new TList();
    fFlowWeights->SetOwner(kTRUE);
    fFlowWeights->SetName("fFlowWeights");
    
    fQAEvents = new TList();
    fQAEvents->SetOwner(kTRUE);
    fQACharged = new TList();
    fQACharged->SetOwner(kTRUE);
    fQAPID = new TList();
    fQAPID->SetOwner(kTRUE);
    
    // creating particle vectors
    fVectorCharged = new std::vector<FlowPart>;
    fVectorPion = new std::vector<FlowPart>;
    fVectorKaon = new std::vector<FlowPart>;
    fVectorProton = new std::vector<FlowPart>;
    
    fVectorCharged->reserve(300);
    if(fProcessPID) { fVectorPion->reserve(200); fVectorKaon->reserve(100); fVectorProton->reserve(100); }
    
    // creating histograms
    // event histogram
    
    fhEventSampling = new TH2D("fhEventSampling","Event sampling; centrality/multiplicity; sample index", fFlowCentNumBins,0,fFlowCentNumBins, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);
    fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); centrality/multiplicity",fMultEstimator.Data()), fFlowCentNumBins,0,fFlowCentNumBins);
    fQAEvents->Add(fhEventCentrality);
    fh2EventCentralityNumSelCharged = new TH2D("fh2EventCentralityNumSelCharged",Form("Event centrality (%s) vs. N^{sel}_{ch}; N^{sel}_{ch}; centrality/multiplicity",fMultEstimator.Data()), 3000,0,3000, fFlowCentNumBins,0,fFlowCentNumBins);
    fQAEvents->Add(fh2EventCentralityNumSelCharged);
    
    const Short_t iEventCounterBins = 11;//10 it was
    TString sEventCounterLabel[iEventCounterBins] = {"Input","Physics selection OK","Centr. Est. Consis. OK","PV OK","SPD Vtx OK","Pileup MV OK","Out-of-bunch Pileup OK","Vtx Consis. OK","PV #it{z} OK","ESD TPC Mult. Diff. OK","Selected"};
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Short_t i(0); i < iEventCounterBins; i++) fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fQAEvents->Add(fhEventCounter);
    
    // flow histograms & profiles
    // weights
    if(fFlowFillWeights || fRunMode == kFillWeights)
    {
        fh3BeforeNUAWeightsRefs = new TH3D("fh3BeforeNUAWeightsRefs","Weights: Refs; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5, 20, -10,10);
        fh3BeforeNUAWeightsRefs->Sumw2();
        fFlowWeights->Add(fh3BeforeNUAWeightsRefs);
        fh3BeforeNUAWeightsCharged = new TH3D("fh3BeforeNUAWeightsCharged","Weights: Charged; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3BeforeNUAWeightsCharged->Sumw2();
        fFlowWeights->Add(fh3BeforeNUAWeightsCharged);
        fh3BeforeNUAWeightsPion = new TH3D("fh3BeforeNUAWeightsPion","Weights: #pi; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3BeforeNUAWeightsPion->Sumw2();
        fFlowWeights->Add(fh3BeforeNUAWeightsPion);
        fh3BeforeNUAWeightsKaon = new TH3D("fh3BeforeNUAWeightsKaon","Weights: K; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3BeforeNUAWeightsKaon->Sumw2();
        fFlowWeights->Add(fh3BeforeNUAWeightsKaon);
        fh3BeforeNUAWeightsProton = new TH3D("fh3BeforeNUAWeightsProton","Weights: p; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3BeforeNUAWeightsProton->Sumw2();
        fFlowWeights->Add(fh3BeforeNUAWeightsProton);
        
        fhBeforeNUEWeightsRefs = new TH1D("fhBeforeNUEWeightsRefs","Weights: Refs; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhBeforeNUEWeightsRefs->Sumw2();
        fFlowWeights->Add(fhBeforeNUEWeightsRefs);
        fhBeforeNUEWeightsCharged = new TH1D("fhBeforeNUEWeightsCharged","Weights: Charged; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhBeforeNUEWeightsCharged->Sumw2();
        fFlowWeights->Add(fhBeforeNUEWeightsCharged);
        fhBeforeNUEWeightsPion = new TH1D("fhBeforeNUEWeightsPion","Weights: #pi; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhBeforeNUEWeightsPion->Sumw2();
        fFlowWeights->Add(fhBeforeNUEWeightsPion);
        fhBeforeNUEWeightsKaon = new TH1D("fhBeforeNUEWeightsKaon","Weights: K; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhBeforeNUEWeightsKaon->Sumw2();
        fFlowWeights->Add(fhBeforeNUEWeightsKaon);
        fhBeforeNUEWeightsProton = new TH1D("fhBeforeNUEWeightsProton","Weights: p; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhBeforeNUEWeightsProton->Sumw2();
        fFlowWeights->Add(fhBeforeNUEWeightsProton);
    }//if(fFlowFillWeights || fRunMode == kFillWeights)
    
    if(fFlowUseNUAWeights)
    {
        fh3AfterNUAWeightsRefs = new TH3D("fh3AfterNUAWeightsRefs","Weights: Refs; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5, 20, -10,10);
        fh3AfterNUAWeightsRefs->Sumw2();
        fFlowWeights->Add(fh3AfterNUAWeightsRefs);
        fh3AfterNUAWeightsCharged = new TH3D("fh3AfterNUAWeightsCharged","Weights: Charged; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3AfterNUAWeightsCharged->Sumw2();
        fFlowWeights->Add(fh3AfterNUAWeightsCharged);
        fh3AfterNUAWeightsPion = new TH3D("fh3AfterNUAWeightsPion","Weights: #pi; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3AfterNUAWeightsPion->Sumw2();
        fFlowWeights->Add(fh3AfterNUAWeightsPion);
        fh3AfterNUAWeightsKaon = new TH3D("fh3AfterNUAWeightsKaon","Weights: K; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3AfterNUAWeightsKaon->Sumw2();
        fFlowWeights->Add(fh3AfterNUAWeightsKaon);
        fh3AfterNUAWeightsProton = new TH3D("fh3AfterNUAWeightsProton","Weights: p; #varphi; #eta; Prim. vtx_{z} (cm)", 100,0,TMath::TwoPi(), 151,-1.5,1.5,20, -10,10);
        fh3AfterNUAWeightsProton->Sumw2();
        fFlowWeights->Add(fh3AfterNUAWeightsProton);
    }//if(fFlowUseNUAWeights)
    if(fFlowUseNUEWeights){
        fhAfterNUEWeightsRefs = new TH1D("fhAfterNUEWeightsRefs","Weights: Refs; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhAfterNUEWeightsRefs->Sumw2();
        fFlowWeights->Add(fhAfterNUEWeightsRefs);
        fhAfterNUEWeightsCharged = new TH1D("fhAfterNUEWeightsCharged","Weights: Charged; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhAfterNUEWeightsCharged->Sumw2();
        fFlowWeights->Add(fhAfterNUEWeightsCharged);
        fhAfterNUEWeightsPion = new TH1D("fhAfterNUEWeightsPion","Weights: #pi; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhAfterNUEWeightsPion->Sumw2();
        fFlowWeights->Add(fhAfterNUEWeightsPion);
        fhAfterNUEWeightsKaon = new TH1D("fhAfterNUEWeightsKaon","Weights: K; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhAfterNUEWeightsKaon->Sumw2();
        fFlowWeights->Add(fhAfterNUEWeightsKaon);
        fhAfterNUEWeightsProton = new TH1D("fhAfterNUEWeightsProton","Weights: p; #it{p}_{T} (GeV/#it{c})",fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fhAfterNUEWeightsProton->Sumw2();
        fFlowWeights->Add(fhAfterNUEWeightsProton);
    }//if(fFlowUseNUEWeights)
    
    
    //mixed correlations
    //
    if(fDoOnlyMixedCorrelations){
        for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics; iMixedHarm++)
        {
            for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)//For now only for nonoverlapping subevents...
            {
                for(Short_t iSample(0); iSample < fNumSamples; iSample++)
                {
                    //reference flow for mixed harmonics 422,633,523
                    if(iMixedHarm!=3){
                        fpMixedRefsCor4[iSample][iGap][iMixedHarm] = new TProfile(Form("fpRefs_<4>_MixedHarm%d_gap%02.2g_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<4>> | Gap %g | v%d | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax);
                        fpMixedRefsCor4[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                        fFlowRefs->Add(fpMixedRefsCor4[iSample][iGap][iMixedHarm]);
                    }
                    if(iMixedHarm==3){
                        //reference flow for mixed harmonics 6222
                        fpMixedRefsCor6[iSample][iGap] = new TProfile(Form("fpRefs_<6>_MixedHarm%d_gap%02.2g_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<6>> | Gap %g | v%d  | sample %d ; centrality/multiplicity;",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax);
                        fpMixedRefsCor6[iSample][iGap]->Sumw2(kTRUE);
                        fFlowRefs->Add(fpMixedRefsCor6[iSample][iGap]);
                    }
                    
                    if(fProcessCharged){
                        if(iMixedHarm!=3){
                            fpMixedChargedCor3Pos[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Charged_<3>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<3'>> | Gap %g | v%d  | POIs pos | sample %d ; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedChargedCor3Pos[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowCharged->Add(fpMixedChargedCor3Pos[iSample][iGap][iMixedHarm]);
                            
                            fpMixedChargedCor3Neg[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Charged_<3>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<3'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedChargedCor3Neg[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowCharged->Add(fpMixedChargedCor3Neg[iSample][iGap][iMixedHarm]);
                        }
                        if(iMixedHarm==3){
                            fpMixedChargedCor4Pos[iSample][iGap] = new TProfile2D(Form("fp2Charged_<4>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<4'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedChargedCor4Pos[iSample][iGap]->Sumw2(kTRUE);
                            fFlowCharged->Add(fpMixedChargedCor4Pos[iSample][iGap]);
                            
                            fpMixedChargedCor4Neg[iSample][iGap] = new TProfile2D(Form("fp2Charged_<4>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<4'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedChargedCor4Neg[iSample][iGap]->Sumw2(kTRUE);
                            fFlowCharged->Add(fpMixedChargedCor4Neg[iSample][iGap]);
                        }
                    }//if(fProcessCharged)
                    if(fProcessPID){
                        if(iMixedHarm!=3){
                            fpMixedPionCor3Pos[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Pion_<3>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Pion: <<3'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedPionCor3Pos[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedPionCor3Pos[iSample][iGap][iMixedHarm]);
                            
                            fpMixedPionCor3Neg[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Pion_<3>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Pion: <<3'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedPionCor3Neg[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedPionCor3Neg[iSample][iGap][iMixedHarm]);
                            
                            fpMixedKaonCor3Pos[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Kaon_<3>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Pion: <<3'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedKaonCor3Pos[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedKaonCor3Pos[iSample][iGap][iMixedHarm]);
                            
                            fpMixedKaonCor3Neg[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Kaon_<3>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Kaon: <<3'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedKaonCor3Neg[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedKaonCor3Neg[iSample][iGap][iMixedHarm]);
                            
                            fpMixedProtonCor3Pos[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Proton_<3>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Proton: <<3'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedProtonCor3Pos[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedProtonCor3Pos[iSample][iGap][iMixedHarm]);
                            
                            fpMixedProtonCor3Neg[iSample][iGap][iMixedHarm] = new TProfile2D(Form("fp2Proton_<3>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Proton: <<3'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedProtonCor3Neg[iSample][iGap][iMixedHarm]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedProtonCor3Neg[iSample][iGap][iMixedHarm]);
                        }
                        if(iMixedHarm==3){
                            fpMixedPionCor4Pos[iSample][iGap] = new TProfile2D(Form("fp2Pion_<4>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Pion: <<4'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedPionCor4Pos[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedPionCor4Pos[iSample][iGap]);
                            
                            fpMixedPionCor4Neg[iSample][iGap] = new TProfile2D(Form("fp2Pion_<4>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Pion: <<4'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedPionCor4Neg[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedPionCor4Neg[iSample][iGap]);
                            
                            fpMixedKaonCor4Pos[iSample][iGap] = new TProfile2D(Form("fp2Kaon_<4>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Kaon: <<4'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedKaonCor4Pos[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedKaonCor4Pos[iSample][iGap]);
                            
                            fpMixedKaonCor4Neg[iSample][iGap] = new TProfile2D(Form("fp2Kaon_<4>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Kaon: <<4'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedKaonCor4Neg[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedKaonCor4Neg[iSample][iGap]);
                            
                            fpMixedProtonCor4Pos[iSample][iGap] = new TProfile2D(Form("fp2Proton_<4>_MixedHarm%d_gap%02.2g_Pos_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Proton: <<4'>> | Gap %g | v%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedProtonCor4Pos[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedProtonCor4Pos[iSample][iGap]);
                            
                            fpMixedProtonCor4Neg[iSample][iGap] = new TProfile2D(Form("fp2Proton_<4>_MixedHarm%d_gap%02.2g_Neg_sample%d",fMixedHarmonics[iMixedHarm],10*fEtaGap[iGap],iSample),Form("Proton: <<4'>> | Gap %g | v%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fMixedHarmonics[iMixedHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                            fpMixedProtonCor4Neg[iSample][iGap]->Sumw2(kTRUE);
                            fFlowPID->Add(fpMixedProtonCor4Neg[iSample][iGap]);
                        }
                    }//if(fProcessPID)
                }//for(Short_t iSample(0); iSample < fNumSamples; iSample++)
            }//for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        }//for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics; iMixedHarm++)
    }//if(fDoOnlyMixedCorrelations)
    
    if(!fDoOnlyMixedCorrelations){
        // correlations
        for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
        {
            for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
            {
                for(Short_t iSample(0); iSample < fNumSamples; iSample++)
                {
                    fpRefsCor2[iSample][iGap][iHarm] = new TProfile(Form("fpRefs_<2>_harm%d_gap%02.2g_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Ref: <<2>> | Gap %g | n=%d | sample %d; centrality/multiplicity;",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax);
                    fpRefsCor2[iSample][iGap][iHarm]->Sumw2(kTRUE);
                    fFlowRefs->Add(fpRefsCor2[iSample][iGap][iHarm]);
                    
                    if(fProcessCharged)
                    {
                        fp2ChargedCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2ChargedCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowCharged->Add(fp2ChargedCor2Pos[iSample][iGap][iHarm]);
                        
                        
                        fp2ChargedCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Charged_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("Charged: <<2'>> | Gap %g | n=%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2ChargedCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowCharged->Add(fp2ChargedCor2Neg[iSample][iGap][iHarm]);
                        
                    }//if(fProcessCharged)
                    
                    if(fProcessPID)
                    {
                        fp2PionCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d   | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2PionCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2PionCor2Pos[iSample][iGap][iHarm]);
                        
                        fp2KaonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2KaonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2KaonCor2Pos[iSample][iGap][iHarm]);
                        
                        fp2ProtonCor2Pos[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Pos_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d  | POIs pos | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2ProtonCor2Pos[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2ProtonCor2Pos[iSample][iGap][iHarm]);
                        
                        fp2PionCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Pion_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID #pi: <<2'>> | Gap %g | n=%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2PionCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2PionCor2Neg[iSample][iGap][iHarm]);
                        
                        fp2KaonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Kaon_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID K: <<2'>> | Gap %g | n=%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2KaonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2KaonCor2Neg[iSample][iGap][iHarm]);
                        
                        fp2ProtonCor2Neg[iSample][iGap][iHarm] = new TProfile2D(Form("fp2Proton_<2>_harm%d_gap%02.2g_Neg_sample%d",fHarmonics[iHarm],10*fEtaGap[iGap],iSample),Form("PID p: <<2'>> | Gap %g | n=%d  | POIs neg | sample %d; centrality/multiplicity; #it{p}_{T} (GeV/c)",fEtaGap[iGap],fHarmonics[iHarm],iSample), fFlowCentNumBins,fFlowCentMin,fFlowCentMax, fFlowPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
                        fp2ProtonCor2Neg[iSample][iGap][iHarm]->Sumw2(kTRUE);
                        fFlowPID->Add(fp2ProtonCor2Neg[iSample][iGap][iHarm]);
                    }//endif(fProcessPID)
                }//for(Short_t iSample(0); iSample < fNumSamples; iSample++)
            }//for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
        }//for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
    }//endif(!fDoOnlyMixedCorrelations)
    
    // charged (tracks) histograms
    fh2RefsMult = new TH2D("fh2RefsMult","RFPs: Centrality: Multiplicity; centrality; multiplicity",fFlowCentNumBins,0,fFlowCentNumBins,1000,0,1000);
    fQACharged->Add(fh2RefsMult);
    fh2RefsPt = new TH2D("fh2RefsPt","RFPs: #it{p}_{T};  centrality; #it{p}_{T} (GeV/#it{c})",fFlowCentNumBins,0,fFlowCentNumBins,100,0,10);
    fQACharged->Add(fh2RefsPt);
    fh2RefsEta = new TH2D("fh2RefsEta","RFPs: #eta; centrality; #eta",fFlowCentNumBins,0,fFlowCentNumBins, 151,-1.5,1.5);
    fQACharged->Add(fh2RefsEta);
    fh2RefsPhi = new TH2D("fh2RefsPhi","RFPs: #varphi; centrality; #varphi",fFlowCentNumBins,0,fFlowCentNumBins,100,0,TMath::TwoPi());
    fQACharged->Add(fh2RefsPhi);
    
    if(fProcessCharged)
    {
        TString sChargedCounterLabel[] = {"Input","Out-of-bunch PLP","FB","#TPC-Cls","TPC-Chi2","DCA-z","DCA-xy","Eta","Selected"};
        const Short_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
        fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
        for(Short_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
        fQACharged->Add(fhChargedCounter);
    } // endif {fProcessCharged}
    
    // PID tracks histograms
    if(fProcessPID)
    {
        fh2PIDPionMult = new TH2D("fh2PIDPionMult","PID: #pi: Multiplicity; centrality; multiplicity", fFlowCentNumBins,0,fFlowCentNumBins,200,0,200);
        fQAPID->Add(fh2PIDPionMult);
        fh2PIDKaonMult = new TH2D("fh2PIDKaonMult","PID: K: Multiplicity; centrality; multiplicity", fFlowCentNumBins,0,fFlowCentNumBins,100,0,100);
        fQAPID->Add(fh2PIDKaonMult);
        fh2PIDProtonMult = new TH2D("fh2PIDProtonMult","PID: p: Multiplicity; centrality; multiplicity", fFlowCentNumBins,0,fFlowCentNumBins,100,0,100);
        fQAPID->Add(fh2PIDProtonMult);
        
        if(fFillQA)
        {
            fh2PIDPionPt = new TH2D("fh2PIDPionPt","PID: #pi: centrality vs. #it{p}_{T}; centrality; #it{p}_{T}", fFlowCentNumBins,0,fFlowCentNumBins,100,0.,10.);
            fQAPID->Add(fh2PIDPionPt);
            fh2PIDPionPhi = new TH2D("fh2PIDPionPhi","PID: #pi: centrality vs. #varphi; centrality; #varphi", fFlowCentNumBins,0,fFlowCentNumBins,100,0,TMath::TwoPi());
            fQAPID->Add(fh2PIDPionPhi);
            fh2PIDPionEta = new TH2D("fh2PIDPionEta","PID: #pi: centrality vs. #eta; centrality; #eta", fFlowCentNumBins,0,fFlowCentNumBins,151,-1.5,1.5);
            fQAPID->Add(fh2PIDPionEta);
            fhPIDPionCharge = new TH1D("fhPIDPionCharge","PID: #pi: charge; charge", 3,-1.5,1.5);
            fQAPID->Add(fhPIDPionCharge);
            fh2PIDKaonPt = new TH2D("fh2PIDKaonPt","PID: K: centrality vs. #it{p}_{T}; centrality; #it{p}_{T}", fFlowCentNumBins,0,fFlowCentNumBins,100,0.,10.);
            fQAPID->Add(fh2PIDKaonPt);
            fh2PIDKaonPhi = new TH2D("fh2PIDKaonPhi","PID: K: centrality vs. #varphi; centrality; #varphi", fFlowCentNumBins,0,fFlowCentNumBins,100,0,TMath::TwoPi());
            fQAPID->Add(fh2PIDKaonPhi);
            fh2PIDKaonEta = new TH2D("fh2PIDKaonEta","PID: K: centrality vs. #eta; centrality; #eta", fFlowCentNumBins,0,fFlowCentNumBins,151,-1.5,1.5);
            fQAPID->Add(fh2PIDKaonEta);
            fhPIDKaonCharge = new TH1D("fhPIDKaonCharge","PID: K: charge; charge", 3,-1.5,1.5);
            fQAPID->Add(fhPIDKaonCharge);
            fh2PIDProtonPt = new TH2D("fh2PIDProtonPt","PID: p: centrality vs. #it{p}_{T}; centrality; #it{p}_{T}", fFlowCentNumBins,0,fFlowCentNumBins,100,0.,10.);
            fQAPID->Add(fh2PIDProtonPt);
            fh2PIDProtonPhi = new TH2D("fh2PIDProtonPhi","PID: p: centrality vs. #varphi; centrality; #varphi", fFlowCentNumBins,0,fFlowCentNumBins,100,0,TMath::TwoPi());
            fQAPID->Add(fh2PIDProtonPhi);
            fh2PIDProtonEta = new TH2D("fh2PIDProtonEta","PID: p: centrality vs. #eta; centrality; #eta", fFlowCentNumBins,0,fFlowCentNumBins,151,-1.5,1.5);
            fQAPID->Add(fh2PIDProtonEta);
            fhPIDProtonCharge = new TH1D("fhPIDProtonCharge","PID: p: charge; charge", 3,-1.5,1.5);
            fQAPID->Add(fhPIDProtonCharge);
            fh2PIDPionTPCdEdx = new TH2D("fh2PIDPionTPCdEdx","PID: #pi: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
            fQAPID->Add(fh2PIDPionTPCdEdx);
            fh2PIDPionTOFbeta = new TH2D("fh2PIDPionTOFbeta","PID: #pi: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
            fQAPID->Add(fh2PIDPionTOFbeta);
            fh2PIDKaonTPCdEdx = new TH2D("fh2PIDKaonTPCdEdx","PID: K: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
            fQAPID->Add(fh2PIDKaonTPCdEdx);
            fh2PIDKaonTOFbeta = new TH2D("fh2PIDKaonTOFbeta","PID: K: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
            fQAPID->Add(fh2PIDKaonTOFbeta);
            fh2PIDProtonTPCdEdx = new TH2D("fh2PIDProtonTPCdEdx","PID: p: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
            fQAPID->Add(fh2PIDProtonTPCdEdx);
            fh2PIDProtonTOFbeta = new TH2D("fh2PIDProtonTOFbeta","PID: p: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
            fQAPID->Add(fh2PIDProtonTOFbeta);
            fh2PIDPionTPCnSigmaPion = new TH2D("fh2PIDPionTPCnSigmaPion","PID: #pi: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTPCnSigmaPion);
            fh2PIDPionTOFnSigmaPion = new TH2D("fh2PIDPionTOFnSigmaPion","PID: #pi: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTOFnSigmaPion);
            fh2PIDPionTPCnSigmaKaon = new TH2D("fh2PIDPionTPCnSigmaKaon","PID: #pi: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTPCnSigmaKaon);
            fh2PIDPionTOFnSigmaKaon = new TH2D("fh2PIDPionTOFnSigmaKaon","PID: #pi: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTOFnSigmaKaon);
            fh2PIDPionTPCnSigmaProton = new TH2D("fh2PIDPionTPCnSigmaProton","PID: #pi: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTPCnSigmaProton);
            fh2PIDPionTOFnSigmaProton = new TH2D("fh2PIDPionTOFnSigmaProton","PID: #pi: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDPionTOFnSigmaProton);
            fh2PIDKaonTPCnSigmaPion = new TH2D("fh2PIDKaonTPCnSigmaPion","PID: K: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTPCnSigmaPion);
            fh2PIDKaonTOFnSigmaPion = new TH2D("fh2PIDKaonTOFnSigmaPion","PID: K: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTOFnSigmaPion);
            fh2PIDKaonTPCnSigmaKaon = new TH2D("fh2PIDKaonTPCnSigmaKaon","PID: K: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTPCnSigmaKaon);
            fh2PIDKaonTOFnSigmaKaon = new TH2D("fh2PIDKaonTOFnSigmaKaon","PID: K: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTOFnSigmaKaon);
            fh2PIDKaonTPCnSigmaProton = new TH2D("fh2PIDKaonTPCnSigmaProton","PID: K: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTPCnSigmaProton);
            fh2PIDKaonTOFnSigmaProton = new TH2D("fh2PIDKaonTOFnSigmaProton","PID: K: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDKaonTOFnSigmaProton);
            fh2PIDProtonTPCnSigmaPion = new TH2D("fh2PIDProtonTPCnSigmaPion","PID: p: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTPCnSigmaPion);
            fh2PIDProtonTOFnSigmaPion = new TH2D("fh2PIDProtonTOFnSigmaPion","PID: p: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTOFnSigmaPion);
            fh2PIDProtonTPCnSigmaKaon = new TH2D("fh2PIDProtonTPCnSigmaKaon","PID: p: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTPCnSigmaKaon);
            fh2PIDProtonTOFnSigmaKaon = new TH2D("fh2PIDProtonTOFnSigmaKaon","PID: p: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTOFnSigmaKaon);
            fh2PIDProtonTPCnSigmaProton = new TH2D("fh2PIDProtonTPCnSigmaProton","PID: p: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTPCnSigmaProton);
            fh2PIDProtonTOFnSigmaProton = new TH2D("fh2PIDProtonTOFnSigmaProton","PID: p: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 42,-11,10);
            fQAPID->Add(fh2PIDProtonTOFnSigmaProton);
        }
        
    } //endif {fProcessPID}
    
    const Short_t iNBinsPIDstatus = 4;
    TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};
    const Short_t iNFilterMapBinBins = 32;
    
    // QA histograms
    if(fFillQA)
    {
        TString sQAindex[fiNumIndexQA] = {"Before", "After"};
        for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
        {
            // EVENTs QA histograms
            fhQAEventsPVz[iQA] = new TH1D(Form("fhQAEventsPVz_%s",sQAindex[iQA].Data()), "QA Events: PV-#it{z}", 101,-50,50);
            fQAEvents->Add(fhQAEventsPVz[iQA]);
            fhQAEventsNumContrPV[iQA] = new TH1D(Form("fhQAEventsNumContrPV_%s",sQAindex[iQA].Data()), "QA Events: Number of contributors to AOD PV", 20,0,20);
            fQAEvents->Add(fhQAEventsNumContrPV[iQA]);
            fhQAEventsNumSPDContrPV[iQA] = new TH1D(Form("fhQAEventsNumSPDContrPV_%s",sQAindex[iQA].Data()), "QA Events: SPD contributors to PV", 20,0,20);
            fQAEvents->Add(fhQAEventsNumSPDContrPV[iQA]);
            fhQAEventsDistPVSPD[iQA] = new TH1D(Form("fhQAEventsDistPVSPD_%s",sQAindex[iQA].Data()), "QA Events: PV SPD vertex", 50,0,5);
            fQAEvents->Add(fhQAEventsDistPVSPD[iQA]);
            fhQAEventsSPDresol[iQA] = new TH1D(Form("fhQAEventsSPDresol_%s",sQAindex[iQA].Data()), "QA Events: SPD resolution", 150,0,15);
            fQAEvents->Add(fhQAEventsSPDresol[iQA]);
            
            fhQAEventsCentralityOutliers[iQA] = new TH2D(Form("fhQAEventsCentralityOutliers_%s",sQAindex[iQA].Data()),"QA Events: Centrality distribution; centrality percentile (V0M);centrality percentile (CL1)",100,0,100,100,0,100);
            fQAEvents->Add(fhQAEventsCentralityOutliers[iQA]);
            
            fhQAEventsPileUp[iQA] = new TH2D(Form("fhQAEventsPileUp_%s",sQAindex[iQA].Data()),"QA Events: TPC vs. ESD multiplicity; TPC multiplicity; ESD multiplicity",500,0,6000,500,0,6000);
            fQAEvents->Add(fhQAEventsPileUp[iQA]);
            fhEventsMultTOFFilterbit32[iQA] = new TH2D(Form("fhEventsMultTOFFilterbit32_%s",sQAindex[iQA].Data()),"filterbit32 vs. TOF multiplicity; multiplicity(fb32);multiplicity(fb32+TOF)", 4000,0,4000,2000,0,2000);
            fQAEvents->Add(fhEventsMultTOFFilterbit32[iQA]);
            // Charged tracks QA
            if(fProcessCharged)
            {
                fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 1500,0,1500);
                fQACharged->Add(fhQAChargedMult[iQA]);
                fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
                fQACharged->Add(fhQAChargedCharge[iQA]);
                fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 100,0.,10.);
                fQACharged->Add(fhQAChargedPt[iQA]);
                fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 151,-1.5,1.5);
                fQACharged->Add(fhQAChargedEta[iQA]);
                fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
                fQACharged->Add(fhQAChargedPhi[iQA]);
                fhQAChargedFilterBit[iQA] = new TH1D(Form("fhQAChargedFilterBit_%s",sQAindex[iQA].Data()), "QA Charged: Filter bit",iNFilterMapBinBins,0,iNFilterMapBinBins);
                for(Int_t j = 0x0; j < iNFilterMapBinBins; j++) fhQAChargedFilterBit[iQA]->GetXaxis()->SetBinLabel(j+1, Form("%g",TMath::Power(2,j)));
                fQACharged->Add(fhQAChargedFilterBit[iQA]);
                fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
                fQACharged->Add(fhQAChargedNumTPCcls[iQA]);
                fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
                fQACharged->Add(fhQAChargedDCAxy[iQA]);
                fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
                fQACharged->Add(fhQAChargedDCAz[iQA]);
            } // endif {fProcessCharged}
            
            // PID tracks QA
            if(fProcessPID)
            {
                fhQAPIDTPCstatus[iQA] = new TH1D(Form("fhQAPIDTPCstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
                fQAPID->Add(fhQAPIDTPCstatus[iQA]);
                fhQAPIDTPCdEdx[iQA] = new TH2D(Form("fhQAPIDTPCdEdx_%s",sQAindex[iQA].Data()),"QA PID: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 100,0,10, 131,-10,1000);
                fQAPID->Add(fhQAPIDTPCdEdx[iQA]);
                fhQAPIDTOFstatus[iQA] = new TH1D(Form("fhQAPIDTOFstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
                fQAPID->Add(fhQAPIDTOFstatus[iQA]);
                fhQAPIDTOFbeta[iQA] = new TH2D(Form("fhQAPIDTOFbeta_%s",sQAindex[iQA].Data()),"QA PID: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 100,0,10, 101,-0.1,1.5);
                fQAPID->Add(fhQAPIDTOFbeta[iQA]);
                
                fh3PIDPionTPCTOFnSigmaPion[iQA] = new TH3D(Form("fh3PIDPionTPCTOFnSigmaPion_%s",sQAindex[iQA].Data()),"QA PID: #pi: TPC-TOF n#sigma (#pi hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDPionTPCTOFnSigmaPion[iQA]);
                fh3PIDPionTPCTOFnSigmaKaon[iQA] = new TH3D(Form("fh3PIDPionTPCTOFnSigmaKaon_%s",sQAindex[iQA].Data()),"QA PID: #pi: TPC-TOF n#sigma (K hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDPionTPCTOFnSigmaKaon[iQA]);
                fh3PIDPionTPCTOFnSigmaProton[iQA] = new TH3D(Form("fh3PIDPionTPCTOFnSigmaProton_%s",sQAindex[iQA].Data()),"QA PID: #pi: TPC-TOF n#sigma (p hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDPionTPCTOFnSigmaProton[iQA]);
                
                fh3PIDKaonTPCTOFnSigmaPion[iQA] = new TH3D(Form("fh3PIDKaonTPCTOFnSigmaPion_%s",sQAindex[iQA].Data()),"QA PID: K: TPC-TOF n#sigma (#pi hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDKaonTPCTOFnSigmaPion[iQA]);
                fh3PIDKaonTPCTOFnSigmaKaon[iQA] = new TH3D(Form("fh3PIDKaonTPCTOFnSigmaKaon_%s",sQAindex[iQA].Data()),"QA PID: K: TPC-TOF n#sigma (K hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDKaonTPCTOFnSigmaKaon[iQA]);
                fh3PIDKaonTPCTOFnSigmaProton[iQA] = new TH3D(Form("fh3PIDKaonTPCTOFnSigmaProton_%s",sQAindex[iQA].Data()),"QA PID: K: TPC-TOF n#sigma (p hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDKaonTPCTOFnSigmaProton[iQA]);
                
                fh3PIDProtonTPCTOFnSigmaPion[iQA] = new TH3D(Form("fh3PIDProtonTPCTOFnSigmaPion_%s",sQAindex[iQA].Data()),"QA PID: p: TPC-TOF n#sigma (#pi hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDProtonTPCTOFnSigmaPion[iQA]);
                fh3PIDProtonTPCTOFnSigmaKaon[iQA] = new TH3D(Form("fh3PIDProtonTPCTOFnSigmaKaon_%s",sQAindex[iQA].Data()),"QA PID: p: TPC-TOF n#sigma (K hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDProtonTPCTOFnSigmaKaon[iQA]);
                fh3PIDProtonTPCTOFnSigmaProton[iQA] = new TH3D(Form("fh3PIDProtonTPCTOFnSigmaProton_%s",sQAindex[iQA].Data()),"QA PID: p: TPC-TOF n#sigma (p hyp.); n#sigma^{TPC}; n#sigma^{TOF}; #it{p}_{T}", 22,-11,10, 22,-11,10, 100,0,10);
                fQAPID->Add(fh3PIDProtonTPCTOFnSigmaProton[iQA]);
                
                
                
                for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
                {
                    fhQAPIDTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
                    fhQAPIDTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
                }
            } // endif {fProcessPID}
        }
    }
    
    // posting data (mandatory)
    PostData(1, fFlowRefs);
    PostData(2, fFlowCharged);
    PostData(3, fFlowPID);
    PostData(4, fQAEvents);
    PostData(5, fQACharged);
    PostData(6, fQAPID);
    PostData(7, fFlowWeights);
    
    return;
}//AliAnalysisTaskFlowModes_pPb::UserCreateOutputObjects()
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::ListParameters()
{
    // lists all task parameters
    // *************************************************************
    printf("\n======= List of parameters ========================================\n");
    printf("   -------- Analysis task ---------------------------------------\n");
    printf("      fRunMode: (RunMode) %d\n",    fRunMode);
    printf("      fAnalType: (AnalType) %d\n",    fAnalType);
    printf("      fFillQA: (Bool_t) %s\n",    fFillQA ? "kTRUE" : "kFALSE");
    printf("      fProcessCharged: (Bool_t) %s\n",    fProcessCharged ? "kTRUE" : "kFALSE");
    printf("      fProcessPID: (Bool_t) %s\n",    fProcessPID ? "kTRUE" : "kFALSE");
    printf("   -------- Flow related ----------------------------------------\n");
    printf("      fCutFlowDoFourCorrelations: (Bool_t) %s\n",    fCutFlowDoFourCorrelations ? "kTRUE" : "kFALSE");
    printf("      fCutFlowRFPsPtMin: (Float_t) %g (GeV/c)\n",    fCutFlowRFPsPtMin);
    printf("      fCutFlowRFPsPtMax: (Float_t) %g (GeV/c)\n",    fCutFlowRFPsPtMax);
    printf("      fFlowPOIsPtMin: (Float_t) %g (GeV/c)\n",    fFlowPOIsPtMin);
    printf("      fFlowPOIsPtMax: (Float_t) %g (GeV/c)\n",    fFlowPOIsPtMax);
    printf("      fFlowCentNumBins: (Int_t) %d (GeV/c)\n",    fFlowCentNumBins);
    printf("      fFlowCentMin: (Int_t) %d (GeV/c)\n",    fFlowCentMin);
    printf("      fFlowCentMax: (Int_t) %d (GeV/c)\n",    fFlowCentMax);
    printf("      fFlowUseNUAWeights: (Bool_t) %s\n",    fFlowUseNUAWeights ? "kTRUE" : "kFALSE");
    printf("      fFlowUseNUEWeights: (Bool_t) %s\n",    fFlowUseNUEWeights ? "kTRUE" : "kFALSE");
    printf("      fPositivelyChargedRef: (Bool_t) %s\n", fPositivelyChargedRef ? "kTRUE" : "kFALSE");
    printf("      fNegativelyChargedRef: (Bool_t) %s\n", fNegativelyChargedRef ? "kTRUE" : "kFALSE");
    printf("      fPositivelyChargedPOI: (Bool_t) %s\n", fPositivelyChargedPOI ? "kTRUE" : "kFALSE");
    printf("      fNegativelyChargedPOI: (Bool_t) %s\n", fNegativelyChargedPOI ? "kTRUE" : "kFALSE");
    printf("      fFlowNUAWeightsPath: (TString) '%s' \n",    fFlowNUAWeightsPath.Data());
    printf("      fFlowNUEWeightsPath: (TString) '%s' \n",    fFlowNUEWeightsPath.Data());
    printf("   -------- Events ----------------------------------------------\n");
    printf("      fTrigger: (Short_t) %d\n",    fTrigger);
    printf("      fMultEstimator: (TString) '%s'\n",    fMultEstimator.Data());
    printf("      fPVtxCutZ: (Double_t) %g (cm)\n",    fPVtxCutZ);
    printf("      fFullCentralityRange: (Bool_t) runs over %s centrality range \n",    fFullCentralityRange ? "0-100%" : "50-100%");
    printf("   -------- Charge tracks ---------------------------------------\n");
    printf("      fCutChargedTrackFilterBit: (UInt) %d\n",    fCutChargedTrackFilterBit);
    printf("      fCutChargedNumTPCclsMin: (UShort_t) %d\n",    fCutChargedNumTPCclsMin);
    printf("      fCutChargedEtaMax: (Float_t) %g\n",    fCutChargedEtaMax);
    printf("      fCutChargedPtMin: (Float_t) %g (GeV/c)\n",    fCutChargedPtMin);
    printf("      fCutChargedPtMax: (Float_t) %g (GeV/c)\n",    fCutChargedPtMax);
    printf("      fCutChargedDCAzMax: (Float_t) %g (cm)\n",    fCutChargedDCAzMax);
    printf("      fCutChargedDCAxyMax: (Float_t) %g (cm)\n",    fCutChargedDCAxyMax);
    printf("   -------- PID (pi,K,p) tracks ---------------------------------\n");
    printf("      fCutPIDUseAntiProtonOnly: (Bool_t) %s\n",  fCutPIDUseAntiProtonOnly ? "kTRUE" : "kFALSE");
    printf("      fCutPIDnSigmaCombinedNoTOFrejection: (Bool_t) %s\n",  fCutPIDnSigmaCombinedNoTOFrejection ? "kTRUE" : "kFALSE");
    printf("      fCutPIDnSigmaPionMax: (Double_t) %g\n",    fCutPIDnSigmaPionMax);
    printf("      fCutPIDnSigmaKaonMax: (Double_t) %g\n",    fCutPIDnSigmaKaonMax);
    printf("=====================================================================\n\n"); 
    
    return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::InitializeTask()
{
    // called once on beginning of task (within UserCreateOutputObjects method)
    // check if task parameters are specified and valid
    // returns kTRUE if succesfull
    // *************************************************************
    
    printf("====== InitializeTask AliAnalysisTaskFlowModes_pPb =========================\n");
    
    if(fAnalType != kESD && fAnalType != kAOD)
    {
        ::Error("InitializeTask","Analysis type not specified! Terminating!");
        return kFALSE;
    }
    
    if(fAnalType == kESD)
    {
        ::Error("InitializeTask","Analysis type: ESD not implemented! Terminating!");
        return kFALSE;
    }
    
    if(fColSystem != kPP && fColSystem != kPbPb && fColSystem != kPPb)
    {
        ::Error("InitializeTask","Collisional system not specified! Terminating!");
        return kFALSE;
    }
    
    // checking PID response
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse)
    {
        ::Error("InitializeTask","AliPIDResponse object not found! Terminating!");
        return kFALSE;
    }
    
    fPIDCombined = new AliPIDCombined();
    if(!fPIDCombined)
    {
        ::Error("InitializeTask","AliPIDCombined object not found! Terminating!");
        return kFALSE;
    }
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetSelectedSpecies(5); // all particle species
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask
    
    // checking cut setting
    ::Info("InitializeTask","Checking task parameters setting conflicts (ranges, etc)");
    if(fCutFlowRFPsPtMin > 0. && fCutFlowRFPsPtMax > 0. && fCutFlowRFPsPtMin > fCutFlowRFPsPtMax)
    {
        ::Error("InitializeTask","Cut: RFPs Pt range wrong!");
        return kFALSE;
    }
    
    // upper-case for multiplicity estimator
    fMultEstimator.ToUpper();
    
    // checking for weights source file
    if(fFlowUseNUAWeights && !fFlowNUAWeightsPath.EqualTo(""))
    {
        fFlowNUAWeightsFile = TFile::Open(Form("alien:///%s",fFlowNUAWeightsPath.Data()));
        if(!fFlowNUAWeightsFile)
        {
            ::Error("InitializeTask","NUA flow weights file not found");
            return kFALSE;
        }
    }
    if(fFlowUseNUEWeights && !fFlowNUEWeightsPath.EqualTo(""))
    {
        fFlowNUEWeightsFile = TFile::Open(Form("alien:///%s",fFlowNUEWeightsPath.Data()));
        if(!fFlowNUEWeightsFile)
        {
            ::Error("InitializeTask","NUE flow weights file not found");
            return kFALSE;
        }
    }
    
    ::Info("InitializeTask","Initialization succesfull!");
    printf("======================================================================\n\n");
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::UserExec(Option_t *)
{
    // main method called for each event (event loop)
    // *************************************************************
    
    if(!fInit) return; // check if initialization succesfull
    
    // local event counter check: if running in test mode, it runs until the 50 events are succesfully processed
    if(fRunMode == kTest && fEventCounter >= fNumEventsAnalyse) return;
    
    // event selection
    fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!EventSelection()) return;
    
    // processing of selected event
    if(!ProcessEvent()) return;
    
    // posting data (mandatory)
    PostData(1, fFlowRefs);
    PostData(2, fFlowCharged);
    PostData(3, fFlowPID);
    PostData(4, fQAEvents);
    PostData(5, fQACharged);
    PostData(6, fQAPID);
    PostData(7, fFlowWeights);
    
    return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::EventSelection()
{
    // main (envelope) method for event selection
    // Specific event selection methods are called from here
    // returns kTRUE if event pass all selection criteria
    // *************************************************************
    
    Bool_t eventSelected = kFALSE;
    
    if(!fEventAOD) return kFALSE;
    
    // Fill event QA BEFORE cuts
    if(fFillQA) FillEventsQA(0);
    
    // event selection for small systems pp in Run2
    if(fColSystem == kPP) eventSelected = IsEventSelected_pp();
    
    // event selection for PbPb in Run2
    if(fColSystem == kPbPb) eventSelected = IsEventSelected_PbPb();
    
    // eventSelected = IsEventSelected_PbPb();
    
    // event selection for pPb in Run2
    if(fColSystem == kPPb) eventSelected = IsEventSelected_pPb();
    return eventSelected;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::IsEventSelected_PbPb()
{
    // Event selection for PbPb collisions recorded in Run 2 year 2015
    // PbPb (LHC15o)
    // return kTRUE if event passes all criteria, kFALSE otherwise
    // *************************************************************
    fhEventCounter->Fill("Input",1);
    
    // Physics selection (trigger)
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
    UInt_t fSelectMask = inputHandler->IsEventSelected();
    
    Bool_t isTriggerSelected = kFALSE;
    switch(fTrigger) // check for high multiplicity trigger
    {
        case 0:
            isTriggerSelected = fSelectMask& AliVEvent::kINT7;
            break;
            
        case 1:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
            break;
            
        case 2:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD;
            break;
            
        default: isTriggerSelected = kFALSE;
    }
    if(!isTriggerSelected){  return kFALSE;}
    
    // events passing physics selection
    fhEventCounter->Fill("Physics selection OK",1);
    
    // get centrality from AliMultSelection
    
    AliMultSelection* fMultSelection = 0x0;
    fMultSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    Double_t centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    
    // cut on consistency between centrality estimators: VOM vs CL1
    
    double fEstimatorsCorrelationCoef[2];
    double fEstimatorsSigmaPars[4];
    double fDeltaEstimatorNsigma[2];
    
    fEstimatorsCorrelationCoef[0] = 0.0157497;
    fEstimatorsCorrelationCoef[1] = 0.973488;
    fEstimatorsSigmaPars[0] = 0.673612;
    fEstimatorsSigmaPars[1] = 0.0290718;
    fEstimatorsSigmaPars[2] = -0.000546728;
    fEstimatorsSigmaPars[3] = 5.82749e-06;
    fDeltaEstimatorNsigma[0] = 5.;
    fDeltaEstimatorNsigma[1] = 5.5;
    
    const double center = centrCL1 * fEstimatorsCorrelationCoef[1] + fEstimatorsCorrelationCoef[0];
    const double sigma = fEstimatorsSigmaPars[0] + fEstimatorsSigmaPars[1] * centrCL1 + fEstimatorsSigmaPars[2] * centrCL1 * centrCL1 + fEstimatorsSigmaPars[3] * centrCL1 * centrCL1 * centrCL1;
    if (centrV0M < center - fDeltaEstimatorNsigma[0] * sigma && centrV0M > center + fDeltaEstimatorNsigma[1] * sigma) {return kFALSE;}
    //if(fabs(centrV0M-centrCL1)>7.5) return kFALSE;
    fhEventCounter->Fill("Centr. Est. Consis. OK",1);
    
    
    // primary vertex selection
    const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
    if(!vtx || vtx->GetNContributors() < 1){ return kFALSE;}
    fhEventCounter->Fill("PV OK",1);
    
    // SPD vertex selection
    const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertexSPD());
    
    
    if (((AliAODHeader*)fEventAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return kFALSE;
    
    if (fEventAOD->IsIncompleteDAQ()) return kFALSE;
    
    if (vtx->GetNContributors() < 2 || vtxSPD->GetNContributors()<1) return kFALSE;
    double cov[6]={0}; double covSPD[6]={0};
    vtx->GetCovarianceMatrix(cov);
    vtxSPD->GetCovarianceMatrix(covSPD);
    
    //Special selections for SPD vertex
    double zRes = TMath::Sqrt(covSPD[5]);
    double fMaxResol=0.25;
    if ( vtxSPD->IsFromVertexerZ() && (zRes>fMaxResol)) return kFALSE;
    fhEventCounter->Fill("SPD Vtx OK",1);
    
    
    // check for multi-vertexer pile-up
    /*
     const int    kMinPileUpContrib = 5;
     const double kMaxPileUpChi2 = 5.0;
     const double kMinWDist = 15;
     const AliAODVertex* vtxPileUp = 0;
     int nPileUp = 0;
     
     nPileUp=fEventAOD->GetNumberOfPileupVerticesTracks();
     
     if (nPileUp) {
     if (vtx == vtxSPD) return kTRUE; // there are pile-up vertices but no primary
     Int_t bcPrim = vtPrm->GetBC();
     for (int iPileUp=0;iPileUp<nPileUp;iPileUp++) {
     vtxPileUp = (const AliAODVertex*)fEventAOD->GetPileupVertexTracks(iPileUp);
     if (vtxPileUp->GetNContributors() < kMinPileUpContrib) continue;
     if (vtxPileUp->GetChi2perNDF() > kMaxPileUpChi2) continue;
     int bcPlp = vtxPileUp->GetBC(); ///newly added
     if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other Bunch crossing (BC)
     double wDst = GetWDist(vtx,vtxPileUp);
     if (wDst<kMinWDist) continue;
     return kTRUE; // pile-up: well separated vertices
     }
     }
     */
    //  fhEventCounter->Fill("Pileup MV OK",1);
    /////////////////////////////
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5);
    utils.SetMinWDistMV(15);
    utils.SetCheckPlpFromDifferentBCMV(kTRUE);
    
    Bool_t isPileupFromMV = utils.IsPileUpMV(fEventAOD);
    
    if(isPileupFromMV) return kFALSE;
    fhEventCounter->Fill("Pileup MV OK",1);
    
    //Bool_t fRejectOutOfBunchPileUp = kFALSE;
    //if(fRejectOutOfBunchPileUp) // out-of-bunch rejection (provided by Christian)
    //{
    //out-of-bunch
    //  if (utils.IsOutOfBunchPileUp(fEventAOD)){ return kFALSE; }
    //}
    fhEventCounter->Fill("Out-of-bunch Pileup OK",1);
    /////////////////////////////
    // check vertex consistency
    double dz = vtx->GetZ() - vtxSPD->GetZ();
    double errTot = TMath::Sqrt(cov[5]+covSPD[5]);
    double err = TMath::Sqrt(cov[5]);
    double nsigTot = dz/errTot;
    double nsig = dz/err;
    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsig)>20) return kFALSE;
    fhEventCounter->Fill("Vtx Consis. OK",1);
    
    
    const Double_t aodVtxZ = vtx->GetZ();
    if( TMath::Abs(aodVtxZ) > fPVtxCutZ ){return kFALSE;}
    fhEventCounter->Fill("PV #it{z} OK",1);
    
    // cut on # ESD tracks vs # TPC only tracks
    const Int_t nTracks = fEventAOD->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)fEventAOD->GetHeader())->GetNumberOfESDTracks();
    Int_t multTrk = 0;
    //Int_t multTrkBefC = 0;
    Int_t multTrkTOF = 0;
    Int_t multTPC = 0;
    for (Int_t it = 0; it < nTracks; it++) {
        AliAODTrack* AODTrk = (AliAODTrack*)fEventAOD->GetTrack(it);
        if (!AODTrk){ delete AODTrk; continue; }
        if (AODTrk->TestFilterBit(128)) {multTPC++;}
        if (AODTrk->TestFilterBit(32)){
            multTrk++;
            if ( TMath::Abs(AODTrk->GetTOFsignalDz()) <= 10 && AODTrk->GetTOFsignal() >= 12000 && AODTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        }
    } // end of for (Int_t it = 0; it < nTracks; it++)
    Double_t multTPCn = multTPC;
    Double_t multEsdn = multEsd;
    Double_t fESDvsTPConlyLinearCut[2];
    
    fESDvsTPConlyLinearCut[0] = 700.;
    fESDvsTPConlyLinearCut[1] = 3.38;
    
    if(fExtraPileUp && (multEsdn > fESDvsTPConlyLinearCut[0] + fESDvsTPConlyLinearCut[1] * multTPCn)) return kFALSE;
    if(!fExtraPileUp && (multEsdn > 15000 + fESDvsTPConlyLinearCut[1] * multTPCn)) return kFALSE;
    fhEventCounter->Fill("ESD TPC Mult. Diff. OK",1);
    
    Int_t fTOFvsFB32nSigmaCut[2];
    fTOFvsFB32nSigmaCut[0] = 4.;
    fTOFvsFB32nSigmaCut[1] = 4.;
    
    Double_t multTrkn = multTrk;
    Double_t multTrkTOFn = multTrkTOF;
    
    Double_t  fTOFvsFB32correlationPars[4];
    Double_t  fTOFvsFB32sigmaPars[6];
    
    fTOFvsFB32correlationPars[0] = -1.0178;
    fTOFvsFB32correlationPars[1] = 0.333132;
    fTOFvsFB32correlationPars[2] = 9.10282e-05;
    fTOFvsFB32correlationPars[3] = -1.61861e-08;
    
    fTOFvsFB32sigmaPars[0] = 1.47848;
    fTOFvsFB32sigmaPars[1] = 0.0385923;
    fTOFvsFB32sigmaPars[2] = -5.06153e-05;
    fTOFvsFB32sigmaPars[3] = 4.37641e-08;
    fTOFvsFB32sigmaPars[4] = -1.69082e-11;
    fTOFvsFB32sigmaPars[5] = 2.35085e-15;
    
    //Double_t mu32tof = PolN(multTrkn,fTOFvsFB32correlationPars,3);
    //Double_t sigma32tof = PolN(multTrkn,fTOFvsFB32sigmaPars, 5);
    
    Double_t mu32tof = fTOFvsFB32correlationPars[0] + fTOFvsFB32correlationPars[1]* multTrkn + fTOFvsFB32correlationPars[2]* pow(multTrkn,2) + fTOFvsFB32correlationPars[3]* pow(multTrkn,3);
    Double_t sigma32tof = fTOFvsFB32sigmaPars[0] + fTOFvsFB32sigmaPars[1]* multTrkn + fTOFvsFB32sigmaPars[2]* pow(multTrkn,2) + fTOFvsFB32sigmaPars[3]* pow(multTrkn,3) + fTOFvsFB32sigmaPars[4]* pow(multTrkn,4) + fTOFvsFB32sigmaPars[5]* pow(multTrkn,5);
    
    
    if (fExtraPileUp && (multTrkTOFn > mu32tof + fTOFvsFB32nSigmaCut[0] * sigma32tof || multTrkTOFn < mu32tof - fTOFvsFB32nSigmaCut[1] * sigma32tof)) return kFALSE;
    //if(fExtraPileUp && multTrkTOFn< (-32+ 0.32*multTrkn+0.000037*multTrkn*multTrkn)) return kFALSE;
    //if(fExtraPileUp && multTrkTOFn> (13+0.46*multTrkn+0.000018*multTrkn*multTrkn)) return kFALSE;
    fhEventCounter->Fill("TOF fb32 Mult. correlation OK",1);
    
    fhEventCounter->Fill("Selected",1);
    
    // Fill event QA AFTER cuts
    if(fFillQA) FillEventsQA(1);
    
    return kTRUE;
    
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::IsEventSelected_pp()
{
    // Event selection for small system collision recorder in Run 2 year 2016
    // pp (LHC16kl), pPb (LHC16rqts)
    // return kTRUE if event passes all criteria, kFALSE otherwise
    // *************************************************************
    
    fhEventCounter->Fill("Input",1);
    
    // Physics selection (trigger)
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
    UInt_t fSelectMask = inputHandler->IsEventSelected();
    
    Bool_t isTriggerSelected = kFALSE;
    switch(fTrigger) // check for high multiplicity trigger
    {
        case 0:
            isTriggerSelected = fSelectMask& AliVEvent::kINT7;
            break;
            
        case 1:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
            break;
            
        case 2:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD;
            break;
            
        default: isTriggerSelected = kFALSE;
    }
    
    if(!isTriggerSelected)
        return kFALSE;
    
    // events passing physics selection
    fhEventCounter->Fill("Physics selection OK",1);
    
    // primary vertex selection
    const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
    if(!vtx || vtx->GetNContributors() < 1)
        return kFALSE;
    fhEventCounter->Fill("PV OK",1);
    
    // SPD vertex selection
    const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertexSPD());
    
    Double_t dMaxResol = 0.25; // suggested from DPG
    Double_t cov[6] = {0};
    vtxSPD->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
    fhEventCounter->Fill("SPD Vtx OK",1);
    
    // PileUp rejection included in Physics selection
    // but with values for high mult pp (> 5 contrib) => for low ones: do manually (> 3 contrib)
    
    /*
     if(fTrigger == 0 && fAOD->IsPileupFromSPD(3,0.8) )
     {
     return kFALSE;
     }
     */
    
    //fhEventCounter->Fill("Pileup SPD OK",1);
    
    // pileup rejection from multivertexer
    
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(5);
    utils.SetMaxPlpChi2MV(5);
    utils.SetMinWDistMV(15);
    utils.SetCheckPlpFromDifferentBCMV(kFALSE);
    Bool_t isPileupFromMV = utils.IsPileUpMV(fEventAOD);
    utils.SetCheckPlpFromDifferentBCMV(kTRUE);
    
    if(isPileupFromMV) return kFALSE;
    fhEventCounter->Fill("Pileup MV OK",1);
    
    //Bool_t fRejectOutOfBunchPileUp = kTRUE; //you have to add this to the task and set from header
    //if(fRejectOutOfBunchPileUp) // out-of-bunch rejection (provided by Christian)
    //{
    //out-of-bunch
    //  if (utils.IsOutOfBunchPileUp(fEventAOD)){ return kFALSE; }
    //}
    fhEventCounter->Fill("Out-of-bunch Pileup OK",1);
    
    //   if (utils.IsSPDClusterVsTrackletBG(fEventAOD))
    //   {
    //     return kFALSE;
    //   }
    //
    //   fhEventCounter->Fill("SPDClTrBG OK",1);
    //
    //   // SPD pileup
    //   if (utils.IsPileUpSPD(fEventAOD))
    //   {
    //     return kFALSE;
    //   }
    //
    //   fhEventCounter->Fill("SPDPU OK",1);
    // }
    
    //fhEventCounter->Fill("Utils OK",1);
    
    // cutting on PV z-distance
    const Double_t aodVtxZ = vtx->GetZ();
    if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
    {
        return kFALSE;
    }
    fhEventCounter->Fill("PV #it{z} OK",1);
    
    fhEventCounter->Fill("Selected",1);
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::IsEventSelected_pPb()
{
    // Event selection for pPb collisions recorded in Run 2 year 2016
    // pPb (LHC16q)
    // return kTRUE if event passes all criteria, kFALSE otherwise
    // *************************************************************
    fhEventCounter->Fill("Input",1);
    
    // Physics selection (trigger)
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
    UInt_t fSelectMask = inputHandler->IsEventSelected();
    
    Bool_t isTriggerSelected = kFALSE;
    switch(fTrigger) // check for high multiplicity trigger
    {
        case 0:
            isTriggerSelected = fSelectMask& AliVEvent::kINT7;
            break;
            
        case 1:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultV0;
            break;
            
        case 2:
            isTriggerSelected = fSelectMask& AliVEvent::kHighMultSPD;
            break;
        case 3:
            isTriggerSelected = fSelectMask& AliVEvent::kINT7;
            break;
            
        default: isTriggerSelected = kFALSE;
    }
    if(!isTriggerSelected){  return kFALSE;}
     
    // events passing physics selection
    fhEventCounter->Fill("Physics selection OK",1);


    // get centrality from AliMultSelection
    
    AliMultSelection* fMultSelection = 0x0;
    fMultSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    Double_t centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    
    // cut on consistency between centrality estimators: VOM vs CL1
    fhEventCounter->Fill("Centr. Est. Consis. OK",1);
    
    
    // primary vertex selection
    const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
    if(!vtx || vtx->GetNContributors() < 1){ return kFALSE;}
    fhEventCounter->Fill("PV OK",1);
    
    // SPD vertex selection
    
    // check for pile up

//    fUtils = new AliAnalysisUtils();
    AliAnalysisUtils utils;
    
    if(fCheckPileUp)
    {
        utils.SetUseMVPlpSelection(kTRUE);
        if(fUsePileUpSPD)
        {
            utils.SetUseMVPlpSelection(kFALSE);
            if(fModifySPDDefaultParams)
            {
                utils.SetMinPlpContribSPD(fMinVtxPileUpContrSPD);
                utils.SetMinPlpZdistSPD(fMinPileUpZdistSPD);
            }
        }
        //fUtils->SetUseOutOfBunchPileUp(kTRUE);
        if(utils.IsPileUpEvent(fEventAOD)) return kFALSE;
    }
    fhEventCounter->Fill("SPD Vtx OK",1);
    
    // check for multi-vertexer pile-up
    /////////////////////////////
    utils.SetMinPlpContribMV(5);    // min. multiplicity of the pile-up vertex to consider
    utils.SetMaxPlpChi2MV(5);   // max chi2 per contributor of the pile-up vertex to consider.
                                // Vertices with too large chi2 are likely to be fake combinatoruals
                                // or debris of other vertices.
    utils.SetMinWDistMV(15);    // minimum weighted distance in Z between 2 vertices
                                // (i.e. (zv1-zv2)/sqrt(sigZv1^2+sigZv2^2) ) to consider for pile-up candidate
    utils.SetCheckPlpFromDifferentBCMV(kFALSE); // if true, the vertex with |BCID|>2 will trigger pile-up tag regarless
                                                // its weighted distance to other vertex (previous setting)
    fhEventCounter->Fill("Pileup MV OK",1);

    if(fRejectOutOfBunchPileUp) // out-of-bunch rejection (provided by Christian)
    {
    //out-of-bunch
      if (utils.IsOutOfBunchPileUp(fEventAOD)){ return kFALSE; }
    }

    fhEventCounter->Fill("Out-of-bunch Pileup OK",1);
    
    fhEventCounter->Fill("Vtx Consis. OK",1);
    
    const Double_t aodVtxZ = vtx->GetZ();
    if( TMath::Abs(aodVtxZ) > fPVtxCutZ ){return kFALSE;}
    fhEventCounter->Fill("PV #it{z} OK",1);
    
    fhEventCounter->Fill("ESD TPC Mult. Diff. OK",1);
    fhEventCounter->Fill("TOF fb32 Mult. correlation OK",1);
    
    fhEventCounter->Fill("Selected",1);
    
    // Fill event QA AFTER cuts
    if(fFillQA) FillEventsQA(1);
    
    return kTRUE;
    
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillEventsQA(const Short_t iQAindex)
{
    // Filling various QA plots related with event selection
    // *************************************************************
    
    const AliAODVertex* aodVtx = fEventAOD->GetPrimaryVertex();
    const Double_t dVtxZ = aodVtx->GetZ();
    const Int_t iNumContr = aodVtx->GetNContributors();
    const AliAODVertex* spdVtx = fEventAOD->GetPrimaryVertexSPD();
    const Int_t iNumContrSPD = spdVtx->GetNContributors();
    const Double_t spdVtxZ = spdVtx->GetZ();
    
    fhQAEventsPVz[iQAindex]->Fill(dVtxZ);
    fhQAEventsNumContrPV[iQAindex]->Fill(iNumContr);
    fhQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
    fhQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));
    
    // SPD vertexer resolution
    Double_t cov[6] = {0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    fhQAEventsSPDresol[iQAindex]->Fill(zRes);
    
    AliMultSelection* fMultSelection = 0x0;
    fMultSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    Double_t centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
/*    
    if(iQAindex==1){
        // cut on consistency between centrality estimators: VOM vs CL1
        double fEstimatorsCorrelationCoef[2];
        double fEstimatorsSigmaPars[4];
        double fDeltaEstimatorNsigma[2];
        
        fEstimatorsCorrelationCoef[0] = 0.0157497;
        fEstimatorsCorrelationCoef[1] = 0.973488;
        fEstimatorsSigmaPars[0] = 0.673612;
        fEstimatorsSigmaPars[1] = 0.0290718;
        fEstimatorsSigmaPars[2] = -0.000546728;
        fEstimatorsSigmaPars[3] = 5.82749e-06;
        fDeltaEstimatorNsigma[0] = 5.;
        fDeltaEstimatorNsigma[1] = 5.5;
        
        const double center = centrCL1 * fEstimatorsCorrelationCoef[1] + fEstimatorsCorrelationCoef[0];
        const double sigma = fEstimatorsSigmaPars[0] + fEstimatorsSigmaPars[1] * centrCL1 + fEstimatorsSigmaPars[2] * centrCL1 * centrCL1 + fEstimatorsSigmaPars[3] * centrCL1 * centrCL1 * centrCL1;
        if (centrV0M < center - fDeltaEstimatorNsigma[0] * sigma && centrV0M > center + fDeltaEstimatorNsigma[1] * sigma) return;
    }
*/
    fhQAEventsCentralityOutliers[iQAindex]->Fill(centrV0M,centrCL1);
    
    const Int_t nTracks = fEventAOD->GetNumberOfTracks();
/*
    Int_t multEsd = ((AliAODHeader*)fEventAOD->GetHeader())->GetNumberOfESDTracks();
    Int_t multTPC = 0;
    Int_t multTrk = 0;
    Int_t multTrkTOF = 0;
    for (Int_t it = 0; it < nTracks; it++) {
        AliAODTrack* AODTrk = (AliAODTrack*)fEventAOD->GetTrack(it);
        if (!AODTrk){ delete AODTrk; continue; }
        if (AODTrk->TestFilterBit(128)) {multTPC++;}
        if (AODTrk->TestFilterBit(32)){
            multTrk++;
            if ( TMath::Abs(AODTrk->GetTOFsignalDz()) <= 10 && AODTrk->GetTOFsignal() >= 12000 && AODTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        }
        
    }// end of for (Int_t it = 0; it < nTracks; it++)
    if(iQAindex==1){
        Double_t multTPCn = multTPC;
        Double_t multEsdn = multEsd;
        
        //Double_t multESDTPCDif = multEsdn - multTPCn*3.38;
        //if (multESDTPCDif > 700.) return;//15000
        
        Double_t fESDvsTPConlyLinearCut[2];
        
        fESDvsTPConlyLinearCut[0] = 700.;
        fESDvsTPConlyLinearCut[1] = 3.38;
        
        if(fExtraPileUp && (multEsdn > fESDvsTPConlyLinearCut[0] + fESDvsTPConlyLinearCut[1] * multTPCn)) return;
        if(!fExtraPileUp && (multEsdn > 15000 + fESDvsTPConlyLinearCut[1] * multTPCn)) return;
        
    }
    fhQAEventsPileUp[iQAindex]->Fill(multTPC,multEsd);
    
    if(iQAindex==1){
        Double_t multTrkn = multTrk;
        Double_t multTrkTOFn = multTrkTOF;
        
        Int_t fTOFvsFB32nSigmaCut[2];
        fTOFvsFB32nSigmaCut[0] = 4.;
        fTOFvsFB32nSigmaCut[1] = 4.;
        
        Double_t  fTOFvsFB32correlationPars[4];
        Double_t  fTOFvsFB32sigmaPars[6];
        
        fTOFvsFB32correlationPars[0] = -1.0178;
        fTOFvsFB32correlationPars[1] = 0.333132;
        fTOFvsFB32correlationPars[2] = 9.10282e-05;
        fTOFvsFB32correlationPars[3] = -1.61861e-08;
        
        fTOFvsFB32sigmaPars[0] = 1.47848;
        fTOFvsFB32sigmaPars[1] = 0.0385923;
        fTOFvsFB32sigmaPars[2] = -5.06153e-05;
        fTOFvsFB32sigmaPars[3] = 4.37641e-08;
        fTOFvsFB32sigmaPars[4] = -1.69082e-11;
        fTOFvsFB32sigmaPars[5] = 2.35085e-15;
        
        //Double_t mu32tof = PolN(multTrkn,fTOFvsFB32correlationPars,3);
        //Double_t sigma32tof = PolN(multTrkn,fTOFvsFB32sigmaPars, 5);
        
        Double_t mu32tof = fTOFvsFB32correlationPars[0] + fTOFvsFB32correlationPars[1]* multTrkn + fTOFvsFB32correlationPars[2]* pow(multTrkn,2) + fTOFvsFB32correlationPars[3]* pow(multTrkn,3);
        Double_t sigma32tof = fTOFvsFB32sigmaPars[0] + fTOFvsFB32sigmaPars[1]* multTrkn + fTOFvsFB32sigmaPars[2]* pow(multTrkn,2) + fTOFvsFB32sigmaPars[3]* pow(multTrkn,3) + fTOFvsFB32sigmaPars[4]* pow(multTrkn,4) + fTOFvsFB32sigmaPars[5]* pow(multTrkn,5);
        
        if (fExtraPileUp && (multTrkTOFn > mu32tof + fTOFvsFB32nSigmaCut[0] * sigma32tof || multTrkTOFn < mu32tof - fTOFvsFB32nSigmaCut[1] * sigma32tof)) return;
        
        //if(fExtraPileUp && multTrkTOFn< (-32+ 0.32*multTrkn+0.000037*multTrkn*multTrkn)) return;
        //if(fExtraPileUp && multTrkTOFn> (13+0.46*multTrkn+0.000018*multTrkn*multTrkn)) return;
    }
    fhEventsMultTOFFilterbit32[iQAindex]->Fill(multTrk,multTrkTOF);
    
*/
    return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::Filtering() //void
{
    
    // main (envelope) method for filtering all particles of interests (POIs) in selected events
    // All POIs passing selection criteria are saved to relevant TClonesArray for further processing
    // return kTRUE if succesfull (no errors in process)
    // *************************************************************
    if(!fProcessCharged && !fProcessPID) // if neither is ON, filtering is skipped
        return kFALSE; //return;
    fVectorCharged->clear();
    FilterCharged();
    
    // estimate centrality & assign indexes (centrality/percentile, ...)
    if(fColSystem == kPbPb  || fColSystem == kPPb){
        
        fIndexCentrality = GetCentralityIndex();
        if(fIndexCentrality < 0) return kFALSE; // return; not succesfull estimation
/*
        Double_t Mult = fVectorCharged->size();
        if(fExtraPileUp && fIndexCentrality< (-1.5*TMath::Power(Mult,0.46)-0.6*TMath::Log(Mult)*TMath::Log(Mult)+81)) return kFALSE;
        if(fExtraPileUp && fIndexCentrality>(-2.3*TMath::Power(Mult,0.39)-0.9*TMath::Log(Mult)*TMath::Log(Mult)+110)) return kFALSE;
*/
        
    }
    if(fColSystem == kPP){fIndexCentrality = 1;}

    
    fh2EventCentralityNumSelCharged->Fill(fVectorCharged->size(),fIndexCentrality);
    
    fhEventCentrality->Fill(fIndexCentrality);
    
    fIndexSampling = GetSamplingIndex();
    
    fhEventSampling->Fill(fIndexCentrality,fIndexSampling);
    
    if(fProcessPID)
    {
        
        fVectorPion->clear();
        fVectorKaon->clear();
        fVectorProton->clear();
        FilterPID();
    }
    
    return kTRUE; //return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FilterCharged()
{
    // Filtering input charged tracks
    // If track passes all requirements as defined in IsChargedSelected(),
    // the relevant properties (pT, eta, phi) are stored in FlowPart struct
    // and pushed to relevant vector container.
    // return kFALSE if any complications occurs
    // *************************************************************
    const Short_t iNumTracks = fEventAOD->GetNumberOfTracks();
    if(iNumTracks < 1) return;
    
    AliAODTrack* track = 0x0;
    Int_t iNumRefs = 0;
    Double_t NUAweight = 0;
    Double_t NUEweight = 0;
    
    for(Short_t iTrack(0); iTrack < iNumTracks; iTrack++)
    {
        track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
        if(!track) continue;
        
        if(fFillQA) FillQACharged(0,track); // QA before selection
        
        if(fNegativelyChargedRef==kTRUE && track->Charge()>0) continue;
        if(fPositivelyChargedRef==kTRUE && track->Charge()<0) continue;
        
        if(IsChargedSelected(track))
        {
            fVectorCharged->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kCharged) );
            
            if(fRunMode == kFillWeights || fFlowFillWeights){
                fh3BeforeNUAWeightsCharged->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ());
                fhBeforeNUEWeightsCharged->Fill(track->Pt());
            }
            if(fFlowUseNUAWeights)
            {
                if(track->Charge()>0){ NUAweight = fh3NUAWeightChargedPlus->GetBinContent( fh3NUAWeightChargedPlus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                if(track->Charge()<0){ NUAweight = fh3NUAWeightChargedMinus->GetBinContent( fh3NUAWeightChargedMinus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                
                fh3AfterNUAWeightsCharged->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ(),NUAweight);
            }
            if(fFlowUseNUEWeights){
                if(track->Charge()>0){NUEweight = fhNUEWeightChargedPlus->GetBinContent( fhNUEWeightChargedPlus->FindBin(track->Pt()));}
                if(track->Charge()<0){NUEweight = fhNUEWeightChargedMinus->GetBinContent( fhNUEWeightChargedMinus->FindBin(track->Pt()));}
                
                fhAfterNUEWeightsCharged->Fill(track->Pt(),NUEweight);
            }
            
            if(fFillQA) FillQACharged(1,track); // QA after selection
            
            // Filling refs QA plots
            if(fCutFlowRFPsPtMin > 0. && track->Pt() >= fCutFlowRFPsPtMin && fCutFlowRFPsPtMax > 0. && track->Pt() <= fCutFlowRFPsPtMax)
            {
                iNumRefs++;
                if(fRunMode == kFillWeights || fFlowFillWeights){
                    fh3BeforeNUAWeightsRefs->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ());
                    fhBeforeNUEWeightsRefs->Fill(track->Pt());
                }
                if(fFlowUseNUAWeights)
                {
                    if(track->Charge()>0){ NUAweight = fh3NUAWeightRefsPlus->GetBinContent( fh3NUAWeightRefsPlus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    if(track->Charge()<0){ NUAweight = fh3NUAWeightRefsMinus->GetBinContent( fh3NUAWeightRefsMinus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    
                    fh3AfterNUAWeightsRefs->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ(),NUAweight);
                }
                if(fFlowUseNUEWeights)
                {
                    if(track->Charge()>0) {NUEweight = fhNUEWeightRefsPlus->GetBinContent( fhNUEWeightRefsPlus->FindBin(track->Pt()) );}
                    if(track->Charge()<0) {NUEweight = fhNUEWeightRefsMinus->GetBinContent( fhNUEWeightRefsMinus->FindBin(track->Pt()) );}
                    
                    //if(track->Charge()>0 && fIndexCentrality<10){NUEweight = fhNUEWeightRefsPlus[fIndexCentrality/5]->GetBinContent( fhNUEWeightRefsPlus[fIndexCentrality/5]->FindBin(track->Pt()) );}
                    //if(track->Charge()>0 && fIndexCentrality>=10 && fIndexCentrality<60){NUEweight = fhNUEWeightRefsPlus[1+fIndexCentrality/10]->GetBinContent( fhNUEWeightRefsPlus[1+fIndexCentrality/10]->FindBin(track->Pt()) );}
                    //if(track->Charge()<0 && fIndexCentrality<10){NUEweight = fhNUEWeightRefsMinus[fIndexCentrality/5]->GetBinContent( fhNUEWeightRefsMinus[fIndexCentrality/5]->FindBin(track->Pt()) );}
                    //if(track->Charge()<0 && fIndexCentrality>=10 && fIndexCentrality<60){NUEweight = fhNUEWeightRefsMinus[1+fIndexCentrality/10]->GetBinContent( fhNUEWeightRefsMinus[1+fIndexCentrality/10]->FindBin(track->Pt()) );}
                    
                    fhAfterNUEWeightsRefs->Fill(track->Pt(),NUEweight);
                }
                FillQARefs(1,track);
            }
        }
    }
    
    // fill QA charged multiplicity
    fh2RefsMult->Fill(fIndexCentrality,iNumRefs);
    if(fFillQA)
    {
        fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
        fhQAChargedMult[1]->Fill(fVectorCharged->size());
    }
    
    return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::IsChargedSelected(const AliAODTrack* track)
{
    // Selection of charged track
    // returns kTRUE if track pass all requirements, kFALSE otherwise
    // *************************************************************
    if(!track) return kFALSE;
    fhChargedCounter->Fill("Input",1);

    //out-of-bunch pileup
    if( !(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) ) return kFALSE;
    if( !track->GetTOFBunchCrossing()==0) return kFALSE;
    fhChargedCounter->Fill("Out-of-bunch PLP",1);

    // filter bit
    if( !track->TestFilterBit(fCutChargedTrackFilterBit) ) return kFALSE;
    fhChargedCounter->Fill("FB",1);
    
    // number of TPC clusters (additional check for not ITS-standalone tracks)
    if( track->GetTPCNcls() < fCutChargedNumTPCclsMin && fCutChargedTrackFilterBit != 2) return kFALSE;
    fhChargedCounter->Fill("#TPC-Cls",1);
    
    // track chi2 per space points
    Double_t chi2TPC =0.;
    chi2TPC = track->Chi2perNDF();
    if (fMaxChi2perTPCcls > 0. && chi2TPC > fMaxChi2perTPCcls) return kFALSE;
    fhChargedCounter->Fill("TPC-Chi2",1);
    
    // track DCA coordinates
    // note AliAODTrack::XYZAtDCA() works only for constrained tracks
    Double_t dTrackXYZ[3] = {0};
    Double_t dVertexXYZ[3] = {0.};
    Double_t dDCAXYZ[3] = {0.};
    if( fCutChargedDCAzMax > 0. || fCutChargedDCAxyMax > 0.)
    {
        const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
        if(!vertex) return kFALSE; // event does not have a PV
        
        track->GetXYZ(dTrackXYZ);
        vertex->GetXYZ(dVertexXYZ);
        
        for(Short_t i(0); i < 3; i++)
            dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
    }
    
    if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fCutChargedDCAzMax) return kFALSE;
    fhChargedCounter->Fill("DCA-z",1);
    
    if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
    fhChargedCounter->Fill("DCA-xy",1);
    
    // pseudorapidity (eta)
    if(fCutChargedEtaMax > 0. && TMath::Abs(track->Eta()) > fCutChargedEtaMax) return kFALSE;
    fhChargedCounter->Fill("Eta",1);
    
    // track passing all criteria
    fhChargedCounter->Fill("Selected",1);
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillQARefs(const Short_t iQAindex, const AliAODTrack* track)
{
    // Filling various QA plots related to RFPs subset of charged track selection
    // *************************************************************
    
    if(!track) return;
   if(iQAindex == 0) return; // NOTE implemented only for selected RFPs
    
    fh2RefsPt->Fill(fIndexCentrality,track->Pt());
    fh2RefsEta->Fill(fIndexCentrality,track->Eta());
    fh2RefsPhi->Fill(fIndexCentrality,track->Phi());
    
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillQACharged(const Short_t iQAindex, const AliAODTrack* track)
{
    // Filling various QA plots related to charged track selection
    // *************************************************************
    if(!track) return;
    
    // filter bit testing
    for(Short_t i(0); i < 32; i++)
    {
        if(track->TestFilterBit(TMath::Power(2.,i)))
            fhQAChargedFilterBit[iQAindex]->Fill(i);
    }
    
    // track charge
    fhQAChargedCharge[iQAindex]->Fill(track->Charge());
    
    // number of TPC clusters
    fhQAChargedNumTPCcls[iQAindex]->Fill(track->GetTPCNcls());
    
    // track DCA
    Double_t dDCAXYZ[3] = {-999., -999., -999.};
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(vertex)
    {
        Double_t dTrackXYZ[3] = {-999., -999., -999.};
        Double_t dVertexXYZ[3] = {-999., -999., -999.};
        
        track->GetXYZ(dTrackXYZ);
        vertex->GetXYZ(dVertexXYZ);
        
        for(Short_t i(0); i < 3; i++)
            dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
    }
    fhQAChargedDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]));
    fhQAChargedDCAz[iQAindex]->Fill(dDCAXYZ[2]);
    
    // kinematics
    fhQAChargedPt[iQAindex]->Fill(track->Pt());
    fhQAChargedPhi[iQAindex]->Fill(track->Phi());
    fhQAChargedEta[iQAindex]->Fill(track->Eta());
    
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FilterPID()
{
    // If track passes all requirements as defined in IsPIDSelected() (and species dependent),
    // the relevant properties (pT, eta, phi) are stored in FlowPart struct
    // and pushed to relevant vector container.
    // return kFALSE if any complications occurs
    // *************************************************************
    
    const Short_t iNumTracks = fEventAOD->GetNumberOfTracks();
    if(iNumTracks < 1) return;
    
    PartSpecies species = kUnknown;
    AliAODTrack* track = 0x0;
    Double_t NUAweight = 0;
    Double_t NUEweight = 0;
    
    for(Short_t iTrack(0); iTrack < iNumTracks; iTrack++)
    {
        track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
        if(!track) continue;
        
        // PID tracks are subset of selected charged tracks (same quality requirements)
        if(!IsChargedSelected(track)) continue;
        
        if(fFillQA) FillPIDQA(0,track,kUnknown);   // filling QA for tracks before selection (but after charged criteria applied)
        
        if(fNegativelyChargedPOI==kTRUE && track->Charge()>0) continue;
        if(fPositivelyChargedPOI==kTRUE && track->Charge()<0) continue;
        
        // PID track selection (return most favourable species)
        PartSpecies species = IsPIDSelected(track);
        // check if only protons should be used
        if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) species = kUnknown;
        
        // selection of PID tracks
        switch (species)
        {
            case kPion:
                fVectorPion->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kPion, fPDGMassPion, track->Px(), track->Py(), track->Pz()) );
                if(fRunMode == kFillWeights || fFlowFillWeights){
                    fh3BeforeNUAWeightsPion->Fill(track->Phi(), track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ());
                    fhBeforeNUEWeightsPion->Fill(track->Pt());
                }
                if(fFlowUseNUAWeights)
                {
                    if(track->Charge() > 0){ NUAweight = fh3NUAWeightPionPlus->GetBinContent( fh3NUAWeightPionPlus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    if(track->Charge() < 0){ NUAweight = fh3NUAWeightPionMinus->GetBinContent( fh3NUAWeightPionMinus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    
                    fh3AfterNUAWeightsPion->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ(),NUAweight);
                }
                if(fFlowUseNUEWeights)
                {
                    if(track->Charge()>0){NUEweight = fhNUEWeightPionPlus->GetBinContent( fhNUEWeightPionPlus->FindBin(track->Pt()));}
                    if(track->Charge() < 0){ NUEweight = fhNUEWeightPionMinus->GetBinContent( fhNUEWeightPionMinus->FindBin(track->Pt()) );}
                    
                    fhAfterNUEWeightsPion->Fill(track->Pt(),NUEweight);
                }
                break;
            case kKaon:
                fVectorKaon->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kKaon, fPDGMassKaon, track->Px(), track->Py(), track->Pz()) );
                if(fRunMode == kFillWeights || fFlowFillWeights){
                    fh3BeforeNUAWeightsKaon->Fill(track->Phi(), track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ());
                    fhBeforeNUEWeightsKaon->Fill(track->Pt());
                }
                if(fFlowUseNUAWeights)
                {
                    if(track->Charge() > 0){NUAweight = fh3NUAWeightKaonPlus->GetBinContent( fh3NUAWeightKaonPlus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    if(track->Charge() < 0){NUAweight = fh3NUAWeightKaonMinus->GetBinContent( fh3NUAWeightKaonMinus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    
                    fh3AfterNUAWeightsKaon->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ(),NUAweight);
                }
                if(fFlowUseNUEWeights)
                {
                    if(track->Charge() > 0){ NUEweight = fhNUEWeightKaonPlus->GetBinContent( fhNUEWeightKaonPlus->FindBin(track->Pt()) );}
                    if(track->Charge() < 0){ NUEweight = fhNUEWeightKaonMinus->GetBinContent( fhNUEWeightKaonMinus->FindBin(track->Pt()) );}
                    
                    fhAfterNUEWeightsKaon->Fill(track->Pt(),NUEweight);
                }
                break;
            case kProton:
                fVectorProton->emplace_back( FlowPart(track->Pt(),track->Phi(),track->Eta(), track->Charge(), kProton, fPDGMassProton, track->Px(), track->Py(), track->Pz()) );
                if(fRunMode == kFillWeights || fFlowFillWeights){
                    fh3BeforeNUAWeightsProton->Fill(track->Phi(), track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ());
                    fhBeforeNUEWeightsProton->Fill(track->Pt());
                }
                if(fFlowUseNUAWeights)
                {
                    if(track->Charge() > 0){NUAweight = fh3NUAWeightProtonPlus->GetBinContent( fh3NUAWeightProtonPlus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    if(track->Charge() < 0){NUAweight = fh3NUAWeightProtonMinus->GetBinContent( fh3NUAWeightProtonMinus->FindBin(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ()) );}
                    
                    fh3AfterNUAWeightsProton->Fill(track->Phi(),track->Eta(),fEventAOD->GetPrimaryVertex()->GetZ(),NUAweight);
                }
                if(fFlowUseNUEWeights)
                {
                    if(track->Charge() > 0){ NUEweight = fhNUEWeightProtonPlus->GetBinContent( fhNUEWeightProtonPlus->FindBin(track->Pt()) );}
                    if(track->Charge() < 0 ){ NUEweight = fhNUEWeightProtonMinus->GetBinContent( fhNUEWeightProtonMinus->FindBin(track->Pt()));}
                    
                    fhAfterNUEWeightsProton->Fill(track->Pt(),NUEweight);
                }
                break;
            default:
                break;
        }
        
        //if(fFillQA) FillPIDQA(1,track,species); // filling QA for tracks AFTER selection
    }
    
    fh2PIDPionMult->Fill(fIndexCentrality,fVectorPion->size());
    fh2PIDKaonMult->Fill(fIndexCentrality,fVectorKaon->size());
    fh2PIDProtonMult->Fill(fIndexCentrality,fVectorProton->size());
    
    return;
}
//_____________________________________________________________________________
AliAnalysisTaskFlowModes_pPb::PartSpecies AliAnalysisTaskFlowModes_pPb::IsPIDSelected(const AliAODTrack* track)
{
    // Selection of PID tracks (pi,K,p) - track identification
    // nSigma cutting is used
    // returns AliAnalysisTaskFlowModes_pPb::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
    // *************************************************************
    
    // checking detector states
    AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
    AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
    
    Bool_t bIsTPCok = (pidStatusTPC == AliPIDResponse::kDetPidOk);
    Bool_t bIsTOFok = ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME) && (track->GetTOFsignal() > 12000) && (track->GetTOFsignal() < 100000)); // checking TOF
    //if(!bIsTPCok) return kUnknown;
    
    const Double_t dPt = track->Pt();
    const Double_t dP = track->P();
    
    // use nSigma cuts (based on combination of TPC / TOF nSigma cuts)
    
    Double_t dNumSigmaTPC[5] = {-99,-99,-99,-99,-99}; // TPC nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    Double_t dNumSigmaTOF[5] = {-99,-99,-99,-99,-99}; // TOF nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    
    Float_t *probabilities;
    Float_t mismProb;
    Float_t ProbBayes[5] = {0,0,0,0,0}; //0=el, 1=mu, 2=pi, 3=ka, 4=pr, 5=deuteron, 6=triton, 7=He3
    
    Double_t probTPC[AliPID::kSPECIES]={0.};
    Double_t probTOF[AliPID::kSPECIES]={0.};
    Double_t probTPCTOF[AliPID::kSPECIES]={0.};
    
    
    // filling nSigma arrays
    if(bIsTPCok) // should be anyway
    {
        dNumSigmaTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
        dNumSigmaTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
        dNumSigmaTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
        dNumSigmaTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
        dNumSigmaTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    }
    
    if(bIsTOFok) // should be anyway
    {
        dNumSigmaTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
        dNumSigmaTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
        dNumSigmaTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
        dNumSigmaTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
        dNumSigmaTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }
    
    
    if(fPIDbayesian){
        if(!TPCTOFagree(track)){return kUnknown;}
        
        fBayesianResponse->SetDetResponse(fEventAOD, fCurrCentr,AliESDpid::kTOF_T0); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
        if(fEventAOD->GetTOFHeader()){
            fESDpid.SetTOFResponse(fEventAOD,AliESDpid::kTOF_T0);
        }
        fBayesianResponse->SetDetAND(1);
    }
    
    Int_t ParticleFlag[]={0,0,0,0};//Unknown (electron,muon), pion, kaon, proton
    
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
    UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
    if (detUsed  == (UInt_t)fPIDCombined->GetDetectorMask() ) {  // TPC is available
        
        // TPC nSigma cuts
        if(dP>0.2 && dP <= 0.5)
        {
            if(fPIDnsigma){
                Double_t dMinSigmasTPC = TMath::MinElement(5,dNumSigmaTPC);
                // electron rejection
                if(dMinSigmasTPC == dNumSigmaTPC[0] && TMath::Abs(dNumSigmaTPC[0]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1; //return kUnknown;
                if(dMinSigmasTPC == dNumSigmaTPC[0] && TMath::Abs(dNumSigmaTPC[1]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1; //return kUnknown;
                if(dMinSigmasTPC == dNumSigmaTPC[2] && TMath::Abs(dNumSigmaTPC[2]) <= fCutPIDnSigmaPionMax) ParticleFlag[1] = 1;//return kPion;
                if(dMinSigmasTPC == dNumSigmaTPC[3] && TMath::Abs(dNumSigmaTPC[3]) <= fCutPIDnSigmaKaonMax) ParticleFlag[2] = 1;//return kKaon;
                if(dMinSigmasTPC == dNumSigmaTPC[4] && TMath::Abs(dNumSigmaTPC[4]) <= fCutPIDnSigmaProtonMax) ParticleFlag[3] = 1;//return kProton;
            }
            if(fPIDbayesian){
                ProbBayes[0] = probTPC[0];
                ProbBayes[1] = probTPC[1];
                ProbBayes[2] = probTPC[2];
                ProbBayes[3] = probTPC[3];
                ProbBayes[4] = probTPC[4];
                
                Double_t dMaxBayesianProb = TMath::MaxElement(5,ProbBayes);
                if(dMaxBayesianProb > fParticleProbability){
                    if(dMaxBayesianProb == ProbBayes[0] && TMath::Abs(dNumSigmaTPC[0]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1;// return kUnknown;
                    if(dMaxBayesianProb == ProbBayes[1] && TMath::Abs(dNumSigmaTPC[1]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1;// return kUnknown;
                    if(dMaxBayesianProb == ProbBayes[2] && TMath::Abs(dNumSigmaTPC[2]) <= fCutPIDnSigmaPionMax) ParticleFlag[1] = 1;//return kPion;
                    if(dMaxBayesianProb == ProbBayes[3] && TMath::Abs(dNumSigmaTPC[3]) <= fCutPIDnSigmaKaonMax) ParticleFlag[2] = 1;//return kKaon;
                    if(dMaxBayesianProb == ProbBayes[4] && TMath::Abs(dNumSigmaTPC[4]) <= fCutPIDnSigmaProtonMax) ParticleFlag[3] = 1;//return kProton;
                }
            }
        }
        
        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
        detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
        if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()){
            AliAODPid* pidObj = track->GetDetPid();
            if ((detUsed >= AliPIDResponse::kDetTOF) && (pidObj && pidObj->GetTOFsignal() < 99999)){
                // combined TPC + TOF nSigma cuts
                if(dP > 0.5) // && < 4 GeV TODO once TPC dEdx parametrisation is available
                {
                    Double_t dNumSigmaCombined[5] = {-99,-99,-99,-99,-99};
                    
                    // discard candidates if no TOF is available if cut is on
                    if(fCutPIDnSigmaCombinedNoTOFrejection && !bIsTOFok) ParticleFlag[0] = 1;// return kUnknown;
                    
                    // calculating combined nSigmas
                    for(Short_t i(0); i < 5; i++)
                    {
                        if(bIsTOFok) { dNumSigmaCombined[i] = TMath::Sqrt(dNumSigmaTPC[i]*dNumSigmaTPC[i] + dNumSigmaTOF[i]*dNumSigmaTOF[i]); }
                        else { dNumSigmaCombined[i] = dNumSigmaTPC[i]; }
                    }
                    
                    if(fPIDnsigma){
                        Double_t dMinSigmasCombined = TMath::MinElement(5,dNumSigmaCombined);
                        
                        // electron rejection
                        if(dMinSigmasCombined == dNumSigmaCombined[0] && TMath::Abs(dNumSigmaCombined[0]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1;// return kUnknown; electron
                        //muon rejection
                        if(dMinSigmasCombined == dNumSigmaCombined[1] && TMath::Abs(dNumSigmaCombined[1]) <= fCutPIDnSigmaTPCRejectElectron) ParticleFlag[0] = 1;// return kUnknown; muon
                        
                        switch (fPIDnsigmaCombination) {
                            case 1:
                                //combination 1
                                if(dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 3.) ParticleFlag[1] = 1;//return kPion;
                                if(dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 2.5) ParticleFlag[2] = 1; //return kKaon;
                                if(dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 3.) ParticleFlag[3] = 1;//return kProton;
                                break;
                            case 2:
                                //combination 2
                                if(dP < 2.5 && dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 3.) ParticleFlag[1] = 1; //return kPion;
                                if(dP > 2.5 && dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 2.) ParticleFlag[1] = 1; //return kPion;
                                
                                if(dP < 2. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 2.5) ParticleFlag[2] = 1; //return kKaon;
                                if(dP > 2. && dP < 3. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 2.) ParticleFlag[2] = 1; //return kKaon;
                                if(dP > 3. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 1.5) ParticleFlag[2] = 1; //return kKaon;
                                
                                if(dP < 3. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 3.) ParticleFlag[3] = 1; //return kProton;
                                if(dP > 3. && dP < 5. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 2.) ParticleFlag[3] = 1; //return kProton;
                                if(dP > 5. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 1.5) ParticleFlag[3] = 1; //return kProton;
                                break;
                            case 3:
                                //combination 3
                                if(dP < 2.5 && dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 3.) ParticleFlag[1] = 1; //return kPion;
                                if(dP > 2.5 && dP < 4. && dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 1.5) ParticleFlag[1] = 1; //return kPion;
                                if(dP > 4. && dMinSigmasCombined == dNumSigmaCombined[2] && TMath::Abs(dNumSigmaCombined[2]) <= 1.) ParticleFlag[1] = 1; //return kPion;
                                
                                if(dP < 2. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 2.5) ParticleFlag[2] = 1; //return kKaon;
                                if(dP > 2. && dP < 3. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 1.5)ParticleFlag[2] = 1; //return kKaon;
                                if(dP > 3. && dMinSigmasCombined == dNumSigmaCombined[3] && TMath::Abs(dNumSigmaCombined[3]) <= 1.) ParticleFlag[2] = 1; //return kKaon;
                                
                                if(dP < 3. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 3.) ParticleFlag[3] = 1; //return kProton;
                                if(dP > 3. && dP < 5. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 2.) ParticleFlag[3] = 1; //return kProton;
                                if(dP > 5. && dMinSigmasCombined == dNumSigmaCombined[4] && TMath::Abs(dNumSigmaCombined[4]) <= 1.) ParticleFlag[3] = 1; //return kProton;
                                break;
                            default:
                                break;
                        }
                    }
                    if(fPIDbayesian){
                        ProbBayes[0] = probTPCTOF[0];
                        ProbBayes[1] = probTPCTOF[1];
                        ProbBayes[2] = probTPCTOF[2];
                        ProbBayes[3] = probTPCTOF[3];
                        ProbBayes[4] = probTPCTOF[4];
                        Double_t dMaxBayesianProb = TMath::MaxElement(5,ProbBayes);
                        if(dMaxBayesianProb > fParticleProbability){
                            if(dMaxBayesianProb == ProbBayes[0] && TMath::Abs(dNumSigmaCombined[0]) <= fCutPIDnSigmaPionMax) ParticleFlag[0] = 1; //return kUnknown;
                            if(dMaxBayesianProb == ProbBayes[1] && TMath::Abs(dNumSigmaCombined[1]) <= fCutPIDnSigmaPionMax) ParticleFlag[0] = 1; //return kUnknown;
                            if(dMaxBayesianProb == ProbBayes[2] && TMath::Abs(dNumSigmaCombined[2]) <= fCutPIDnSigmaPionMax) ParticleFlag[1] = 1; //return kPion;
                            if(dMaxBayesianProb == ProbBayes[3] && TMath::Abs(dNumSigmaCombined[3]) <= fCutPIDnSigmaKaonMax) ParticleFlag[2] = 1; //return kKaon;
                            if(dMaxBayesianProb == ProbBayes[4] && TMath::Abs(dNumSigmaCombined[4]) <= fCutPIDnSigmaProtonMax) ParticleFlag[3] = 1; //return kProton;
                        }else{ParticleFlag[0] = 1;}
                    }
                }//if(dP > 0.5)
            }//if ((detUsed >= AliPIDResponse::kDetTOF) && (pidObj && pidObj->GetTOFsignal() < 99999))
        }//TPC or TOF if (detUsed)
    }//TPC if (detUsed)
    
    PartSpecies species;
    if(!ParticleFlag[0] && ParticleFlag[1] && !ParticleFlag[2] && !ParticleFlag[3]){
        species = kPion;
    }else if(!ParticleFlag[0] && !ParticleFlag[1] && ParticleFlag[2] && !ParticleFlag[3]){
        species = kKaon;
    }else if(!ParticleFlag[0] && !ParticleFlag[1] && !ParticleFlag[2] && ParticleFlag[3]){
        species = kProton;
    }else{species = kUnknown;}
    
    if(fFillQA) FillPIDQA(1,track,species); // filling QA for tracks AFTER selection
    
    return species;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillPIDQA(const Short_t iQAindex, const AliAODTrack* track, const PartSpecies species)
{
    // Filling various QA plots related to PID (pi,K,p) track selection
    // *************************************************************
    if(!track) return;
    if(!fPIDResponse || !fPIDCombined)
    {
        ::Error("FillPIDQA","AliPIDResponse or AliPIDCombined object not found!");
        return;
    }
    
    // TPC & TOF statuses & measures
    AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
    AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
    
    fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
    fhQAPIDTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );
    
    Bool_t bIsTPCok = (pidStatusTPC == AliPIDResponse::kDetPidOk);
    //Bool_t bIsTOFok = ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
    
    Bool_t bIsTOFok = ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME) && (track->GetTOFsignal() > 12000) && (track->GetTOFsignal() < 100000)); // checking TOF
    
    Double_t dNumSigmaTPC[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    Double_t dNumSigmaTOF[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    
    Double_t dTPCdEdx = -5; // TPC dEdx for selected particle
    Double_t dTOFbeta = -0.05; //TOF beta for selected particle
    
    Double_t dP = track->P();
    Double_t dPt = track->Pt();
    
    // detector status dependent
    if(bIsTPCok)
    {
        dNumSigmaTPC[0] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        dNumSigmaTPC[1] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
        dNumSigmaTPC[2] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        dNumSigmaTPC[3] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        dNumSigmaTPC[4] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        
        dTPCdEdx = track->GetTPCsignal();
        fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), dTPCdEdx);
    }
    else // TPC status not OK
    {
        dNumSigmaTPC[0] = -11.;
        dNumSigmaTPC[1] = -11.;
        dNumSigmaTPC[2] = -11.;
        dNumSigmaTPC[3] = -11.;
        dNumSigmaTPC[4] = -11.;
        
        fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), -5.);
    }
    
    if(bIsTOFok)
    {
        dNumSigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        dNumSigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
        dNumSigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        dNumSigmaTOF[3] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        dNumSigmaTOF[4] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        
        Double_t dTOF[5];
        track->GetIntegratedTimes(dTOF);
        dTOFbeta = dTOF[0] / track->GetTOFsignal();
        fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);
    }
    else // TOF status not OK
    {
        dNumSigmaTOF[0] = -11.;
        dNumSigmaTOF[1] = -11.;
        dNumSigmaTOF[2] = -11.;
        dNumSigmaTOF[3] = -11.;
        dNumSigmaTOF[4] = -11.;
        
        fhQAPIDTOFbeta[iQAindex]->Fill(track->P(),-0.05);
    }
    
    
    // species dependent QA
    switch (species)
    {
        case kPion:
            fh2PIDPionPt->Fill(fIndexCentrality,track->Pt());
            fh2PIDPionPhi->Fill(fIndexCentrality,track->Phi());
            fh2PIDPionEta->Fill(fIndexCentrality,track->Eta());
            fhPIDPionCharge->Fill(track->Charge());
            fh2PIDPionTPCdEdx->Fill(dPt,dTPCdEdx);
            fh2PIDPionTOFbeta->Fill(dPt,dTOFbeta);
            fh2PIDPionTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
            fh2PIDPionTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
            fh2PIDPionTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
            fh2PIDPionTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
            fh2PIDPionTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
            fh2PIDPionTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
            
            fh3PIDPionTPCTOFnSigmaPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dPt);
            fh3PIDPionTPCTOFnSigmaKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dPt);
            fh3PIDPionTPCTOFnSigmaProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dPt);
            
            break;
            
        case kKaon:
            fh2PIDKaonPt->Fill(fIndexCentrality,track->Pt());
            fh2PIDKaonPhi->Fill(fIndexCentrality,track->Phi());
            fh2PIDKaonEta->Fill(fIndexCentrality,track->Eta());
            fhPIDKaonCharge->Fill(track->Charge());
            fh2PIDKaonTPCdEdx->Fill(dP,dTPCdEdx);
            fh2PIDKaonTOFbeta->Fill(dP,dTOFbeta);
            fh2PIDKaonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
            fh2PIDKaonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
            fh2PIDKaonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
            fh2PIDKaonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
            fh2PIDKaonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
            fh2PIDKaonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
            
            fh3PIDKaonTPCTOFnSigmaPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dPt);
            fh3PIDKaonTPCTOFnSigmaKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dPt);
            fh3PIDKaonTPCTOFnSigmaProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dPt);
            
            break;
            
        case kProton:
            fh2PIDProtonPt->Fill(fIndexCentrality,track->Pt());
            fh2PIDProtonPhi->Fill(fIndexCentrality,track->Phi());
            fh2PIDProtonEta->Fill(fIndexCentrality,track->Eta());
            fhPIDProtonCharge->Fill(track->Charge());
            fh2PIDProtonTPCdEdx->Fill(dP,dTPCdEdx);
            fh2PIDProtonTOFbeta->Fill(dP,dTOFbeta);
            fh2PIDProtonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
            fh2PIDProtonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
            fh2PIDProtonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
            fh2PIDProtonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
            fh2PIDProtonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
            fh2PIDProtonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
            
            fh3PIDProtonTPCTOFnSigmaPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dPt);
            fh3PIDProtonTPCTOFnSigmaKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dPt);
            fh3PIDProtonTPCTOFnSigmaProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dPt);
            
            break;
            
        default:
            break;
    }
    
    
    //
    
    
    return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskFlowModes_pPb::ProcessEvent()
{
    
    // main method for processing of (selected) event:
    // - Filtering of tracks / particles for flow calculations
    // - Phi,eta,pt weights for generic framework are calculated if specified
    // - Flow calculations
    // returns kTRUE if succesfull
    // *************************************************************
    
    // printf("======= EVENT ================\n");
    
    
    // checking the run number for applying weights & loading TList with weights
    if(fRunNumber < 0 || fRunNumber != fEventAOD->GetRunNumber() )
    {
        fRunNumber = fEventAOD->GetRunNumber();
        
        if(fFlowUseNUAWeights && fFlowNUAWeightsFile)
        {
            TDirectory* dirFlowNUAWeights = (TDirectory*) fFlowNUAWeightsFile->Get(Form("000%d",fRunNumber));
            if(!dirFlowNUAWeights) {::Error("ProcessEvent","TList from flow weights not found."); return kFALSE; }
            fh3NUAWeightRefsPlus = (TH3D*) dirFlowNUAWeights->Get("ChargedPlus"); if(!fh3NUAWeightRefsPlus) { ::Error("ProcessEvent","Positive Refs NUA weights not found"); return kFALSE; }
            fh3NUAWeightRefsMinus = (TH3D*) dirFlowNUAWeights->Get("ChargedMinus"); if(!fh3NUAWeightRefsMinus) { ::Error("ProcessEvent","Negative Refs NUA weights not found"); return kFALSE; }
            
            fh3NUAWeightChargedPlus = (TH3D*) dirFlowNUAWeights->Get("ChargedPlus"); if(!fh3NUAWeightChargedPlus) { ::Error("ProcessEvent","Positive Charged NUA weights not found"); return kFALSE; }
            fh3NUAWeightChargedMinus = (TH3D*) dirFlowNUAWeights->Get("ChargedMinus"); if(!fh3NUAWeightChargedMinus) { ::Error("ProcessEvent","Nagative Charged NUA weights not found"); return kFALSE; }
            
            fh3NUAWeightPionPlus = (TH3D*) dirFlowNUAWeights->Get("PionPlus"); if(!fh3NUAWeightPionPlus) { ::Error("ProcessEvent","Positive Pion NUA weights not found"); return kFALSE; }
            fh3NUAWeightPionMinus = (TH3D*) dirFlowNUAWeights->Get("PionMinus"); if(!fh3NUAWeightPionMinus) { ::Error("ProcessEvent","Negative Pion NUA weights not found"); return kFALSE; }
            
            
            fh3NUAWeightKaonPlus = (TH3D*) dirFlowNUAWeights->Get("KaonPlus"); if(!fh3NUAWeightKaonPlus) { ::Error("ProcessEvent","Positive Kaon NUA weights not found"); return kFALSE; }
            fh3NUAWeightKaonMinus = (TH3D*) dirFlowNUAWeights->Get("KaonMinus"); if(!fh3NUAWeightKaonMinus) { ::Error("ProcessEvent","Negative Kaon NUA weights not found"); return kFALSE; }
            
            
            fh3NUAWeightProtonPlus = (TH3D*) dirFlowNUAWeights->Get("ProtonPlus"); if(!fh3NUAWeightProtonPlus) { ::Error("ProcessEvent","Positive Proton NUA weights not found"); return kFALSE; }
            fh3NUAWeightProtonMinus = (TH3D*) dirFlowNUAWeights->Get("ProtonMinus"); if(!fh3NUAWeightProtonMinus) { ::Error("ProcessEvent","Negative Proton NUA weights not found"); return kFALSE; }
            
            
        }
    }
    
    if(fColSystem == kPbPb || fColSystem == kPPb){
        fIndexCentrality = GetCentralityIndex();
        if(fIndexCentrality < 0) { return kFALSE;}
    }
    if(fFlowUseNUEWeights && fFlowNUEWeightsFile)
    {
        //TDirectory* dirFlowNUEWeights = 0x0;
        const char* gCentrality[] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60"};
        TString indexCentrality;// = 0x0;
        if(fIndexCentrality>=0. && fIndexCentrality<=5.){indexCentrality = "0-5";}//gCentrality[0];}
        if(fIndexCentrality>5. && fIndexCentrality<=10.){indexCentrality = "5-10";}//gCentrality[1];}
        if(fIndexCentrality>10. && fIndexCentrality<=20.){indexCentrality = "10-20";}//gCentrality[2];}
        if(fIndexCentrality>20. && fIndexCentrality<=30.){indexCentrality = "20-30";}//gCentrality[3];}
        if(fIndexCentrality>30. && fIndexCentrality<=40.){indexCentrality = "30-40";}//gCentrality[4];}
        if(fIndexCentrality>40. && fIndexCentrality<=50.){indexCentrality = "40-50";}//gCentrality[5];}
        if(fIndexCentrality>50. && fIndexCentrality<=60.){indexCentrality = "50-60";}//gCentrality[6];}
        if(fIndexCentrality<0. || fIndexCentrality>60.){return kFALSE;}
        
        TDirectory* dirFlowNUEWeights = (TDirectory*) fFlowNUEWeightsFile->Get(indexCentrality);
        if(!dirFlowNUEWeights) { ::Error("ProcessEvent","TDirectoy from NUE weights not found."); return kFALSE; }
        fhNUEWeightRefsPlus = dynamic_cast<TH1F*>( dirFlowNUEWeights->Get("ChargedPlus"));
        if(!fhNUEWeightRefsPlus) { ::Error("ProcessEvent","Positive Refs NUE weights not found"); return kFALSE; }
        fhNUEWeightRefsMinus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("ChargedMinus"));
        if(!fhNUEWeightRefsMinus) { ::Error("ProcessEvent","Negative Refs NUE weights not found"); return kFALSE; }
        fhNUEWeightChargedPlus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("ChargedPlus"));
        if(!fhNUEWeightChargedPlus) { ::Error("ProcessEvent","Positive Charged NUE weights not found"); return kFALSE; }
        fhNUEWeightChargedMinus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("ChargedMinus"));
        if(!fhNUEWeightChargedMinus) { ::Error("ProcessEvent","Negative Charged NUE weights not found"); return kFALSE; }
        fhNUEWeightPionPlus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("PionPlus"));
        if(!fhNUEWeightPionPlus) { ::Error("ProcessEvent","Positive Pion NUE weights not found"); return kFALSE; }
        fhNUEWeightPionMinus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("PionMinus"));
        if(!fhNUEWeightPionMinus) { ::Error("ProcessEvent","Negative Pion NUE weights not found"); return kFALSE; }
        fhNUEWeightKaonPlus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("KaonPlus"));
        if(!fhNUEWeightKaonPlus) { ::Error("ProcessEvent","Positive Kaon NUE weights not found"); return kFALSE; }
        fhNUEWeightKaonMinus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("KaonMinus"));
        if(!fhNUEWeightKaonMinus) { ::Error("ProcessEvent","Negative Kaon NUE weights not found"); return kFALSE; }
        fhNUEWeightProtonPlus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("ProtonPlus"));
        if(!fhNUEWeightProtonPlus) { ::Error("ProcessEvent","Positive Proton NUE weights not found"); return kFALSE; }
        fhNUEWeightProtonMinus = dynamic_cast<TH1F*>(dirFlowNUEWeights->Get("ProtonMinus"));
        if(!fhNUEWeightProtonMinus) { ::Error("ProcessEvent","Negative Proton NUE weights not found"); return kFALSE; }
    }
    
    // filtering particles
    if(!Filtering()) return kFALSE;
    // at this point, centrality index (percentile) should be properly estimated, if not, skip event
    
    // if running in kFillWeights mode, skip the remaining part
    if(fRunMode == kFillWeights) { fEventCounter++; return kTRUE; }
    
    // checking if there is at least one charged track selected;
    // if not, event is skipped: unable to compute Reference flow (and thus any differential flow)
    if(fVectorCharged->size() < 1)
        return kFALSE;
    // at this point, all particles fullfiling relevant POIs (and REFs) criteria are filled in TClonesArrays
    
    // >>>> flow starts here <<<<
    // >>>> Flow a la General Framework <<<<
    for(Short_t iGap(0); iGap < fNumEtaGap; iGap++)
    {
        
        // Reference (pT integrated) flow
        DoFlowRefs(iGap);
        
        // pT differential
        if(fProcessCharged)
        {
            // charged track flow
            DoFlowCharged(iGap);
        }
        if(fProcessPID)
        {
            const Int_t iSizePion = fVectorPion->size();
            const Int_t iSizeKaon = fVectorKaon->size();
            const Int_t iSizeProton = fVectorProton->size();
            
            // pi,K,p flow
            if(iSizePion > 0) DoFlowPID(iGap,kPion);
            if(iSizeKaon > 0) DoFlowPID(iGap,kKaon);
            if(iSizeProton > 0) DoFlowPID(iGap,kProton);
        }
    } // endfor {iGap} eta gaps
    
    fEventCounter++; // counter of processed events
    
    return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::DoFlowRefs(const Short_t iEtaGapIndex)
{
    
    // Estimate <2>, <4> and <6> for reference flow for all harmonics based on relevant flow vectors
    // *************************************************************
    
    //Float_t dEtaGap = fEtaGap[iEtaGapIndex];
    Short_t iHarmonics = 0;
    Double_t Cn2 = 0;
    Double_t Dn4GapP = 0;
    Double_t Dn6GapP = 0;
    TComplex vector = TComplex(0,0,kFALSE);
    Double_t dValue = 999;
    
    FillRefsVectors(iEtaGapIndex); // filling RFPs (Q) flow vectors
    if(!fDoOnlyMixedCorrelations){
        // estimating <2>
        Cn2 = TwoGap(0,0).Re();
        if(Cn2 != 0)
        {
            for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
            {
                iHarmonics = fHarmonics[iHarm];
                vector = TwoGap(iHarmonics,-iHarmonics);
                dValue = vector.Re()/Cn2;
                // printf("Gap (RFPs): %g Harm %d | Dn2: %g | fFlowVecQpos[0][0]: %g | fFlowVecQneg[0][0]: %g | fIndexCentrality %d\n\n", dEtaGap,iHarmonics,Cn2,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re(),fIndexCentrality);
                if( TMath::Abs(dValue < 1) )
                    fpRefsCor2[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, dValue, Cn2);
                
            }
        }
    }
    //for mixed harmonics
    if(fDoOnlyMixedCorrelations){
        //estimating <4>
        Dn4GapP = FourGapPos(0,0,0,0).Re();
        if(Dn4GapP != 0)
        {
            // (2,2 | 2,2)_gap , referece flow for v4/psi2
            TComplex Four_2222_GapP = FourGapPos(2, 2, -2, -2);
            double c4_2222_GapP = Four_2222_GapP.Re()/Dn4GapP;
            fpMixedRefsCor4[fIndexSampling][iEtaGapIndex][0]->Fill(fIndexCentrality, c4_2222_GapP, Dn4GapP);
            // (3,3 | 3,3)_gap, referece flow for v6/psi3
            TComplex Four_3333_GapP = FourGapPos(3, 3, -3, -3);
            double c4_3333_GapP = Four_3333_GapP.Re()/Dn4GapP;
            fpMixedRefsCor4[fIndexSampling][iEtaGapIndex][1]->Fill(fIndexCentrality, c4_3333_GapP, Dn4GapP);
            // (3,2 | 3,2)_gap, reference flow for v5/psi23
            TComplex Four_3232_GapP = FourGapPos(3, 2, -3, -2);
            double c4_3232_GapP = Four_3232_GapP.Re()/Dn4GapP;
            fpMixedRefsCor4[fIndexSampling][iEtaGapIndex][2]->Fill(fIndexCentrality, c4_3232_GapP, Dn4GapP);
        }
        //estimating <6>
        Dn6GapP = SixGapPos(0,0,0,0,0,0).Re();
        if(Dn6GapP != 0)
        {
            // (2,2,2 | 2,2,2)_gap, reference flow for v6/psi2
            TComplex Six_222222_GapP = SixGapPos(2, 2, 2, -2, -2, -2);
            double c6_222222_GapP = Six_222222_GapP.Re()/Dn6GapP;
            fpMixedRefsCor6[fIndexSampling][iEtaGapIndex]->Fill(fIndexCentrality, c6_222222_GapP, Dn6GapP);
        }
    }
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::DoFlowCharged(const Short_t iEtaGapIndex)
{
    // Estimate <2'>, <3'> and <4'> for pT diff flow of charged tracks for all harmonics based on relevant flow vectors
    // *************************************************************
    
    FillPOIsVectors(iEtaGapIndex,kCharged);  // filling POIs (P,S) flow vectors
    
    const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
    
    //Float_t dEtaGap = fEtaGap[iEtaGapIndex];
    Short_t iHarmonics = 0;
    Double_t Dn2 = 0;
    Double_t DDn3GapP=0;
    Double_t DDn4GapP=0;
    TComplex vector = TComplex(0,0,kFALSE);
    Double_t dValue = 999;
    
    
    for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
    {
        if(!fDoOnlyMixedCorrelations){
            // estimating <2'>
            // POIs in positive eta
            Dn2 = TwoDiffGapPos(0,0,iPt).Re();
            if(Dn2 != 0)
            {
                for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
                {
                    iHarmonics = fHarmonics[iHarm];
                    vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
                    dValue = vector.Re()/Dn2;
                    if( TMath::Abs(dValue < 1) )
                        fp2ChargedCor2Pos[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
                }
            }
            
            // POIs in negative eta
            Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
            if(Dn2 != 0)
            {
                for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
                {
                    iHarmonics = fHarmonics[iHarm];
                    vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
                    dValue = vector.Re()/Dn2;
                    if( TMath::Abs(dValue < 1) )
                        fp2ChargedCor2Neg[fIndexSampling][iEtaGapIndex][iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
                }
            }
        }
        if(fDoOnlyMixedCorrelations){
            // estimating <3'>
            // POIs in positive eta
            DDn3GapP = ThreeDiffGapPos(0,0,0,iPt).Re();
            if(DDn3GapP!=0)
            {
                for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics-1; iMixedHarm++)
                {
                    if(iMixedHarm==0){ vector = ThreeDiffGapPos(4,-2,-2,iPt);}
                    if(iMixedHarm==1){ vector = ThreeDiffGapPos(6,-3,-3,iPt);}
                    if(iMixedHarm==2){ vector = ThreeDiffGapPos(5,-3,-2,iPt);}
                    dValue = vector.Re()/DDn3GapP;
                    if( TMath::Abs(dValue < 1) ){
                        fpMixedChargedCor3Pos[fIndexSampling][iEtaGapIndex][iMixedHarm]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn3GapP);
                    }
                }
            }
            // POIs in negative eta
            DDn3GapP = ThreeDiffGapNeg(0,0,0,iPt).Re();
            if(DDn3GapP!=0)
            {
                for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics-1; iMixedHarm++)
                {
                    if(iMixedHarm==0){ vector = ThreeDiffGapNeg(4,-2,-2,iPt);}
                    if(iMixedHarm==1){ vector = ThreeDiffGapNeg(6,-3,-3,iPt);}
                    if(iMixedHarm==2){ vector = ThreeDiffGapNeg(5,-2,-3,iPt);}
                    dValue = vector.Re()/DDn3GapP;
                    if( TMath::Abs(dValue < 1) ){
                        fpMixedChargedCor3Neg[fIndexSampling][iEtaGapIndex][iMixedHarm]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn3GapP);
                    }
                }
            }
            // estimating <4'>
            // POIs in positive eta
            DDn4GapP = Four13DiffGapPos(0,0,0,0,iPt).Re();
            if(DDn4GapP!=0)
            {
                vector = Four13DiffGapPos(6,-2,-2,-2,iPt);
                dValue = vector.Re()/DDn4GapP;
                if( TMath::Abs(dValue < 1) ){
                    fpMixedChargedCor4Pos[fIndexSampling][iEtaGapIndex]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn4GapP);
                }
            }
            // POIs in negative eta
            DDn4GapP = Four13DiffGapNeg(0,0,0,0,iPt).Re();
            if(DDn4GapP!=0)
            {
                vector = Four13DiffGapNeg(6,-2,-2,-2,iPt);
                dValue = vector.Re()/DDn4GapP;
                if( TMath::Abs(dValue < 1) ){
                    fpMixedChargedCor4Neg[fIndexSampling][iEtaGapIndex]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn4GapP);
                }
            }
        }
    } // endfor {iPt}
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::DoFlowPID(const Short_t iEtaGapIndex, const PartSpecies species)
{
    // Estimate <2'>, <3'> and <4'> for pT diff flow of pi/K/p tracks for all harmonics based on relevant flow vectors
    // *************************************************************
    
    TProfile2D** profile2Pos = 0x0;
    TProfile2D** profile2Neg = 0x0;
    //TProfile2D** profile4 = 0x0;
    TProfile2D** profile3Pos = 0x0;
    TProfile2D** profile3Neg = 0x0;
    TProfile2D* profile4Pos = 0x0;
    TProfile2D* profile4Neg = 0x0;
    
    
    if(!fDoOnlyMixedCorrelations){
        
        switch (species)
        {
            case kPion:
                profile2Pos = fp2PionCor2Pos[fIndexSampling][iEtaGapIndex];
                profile2Neg = fp2PionCor2Neg[fIndexSampling][iEtaGapIndex];
                //profile4 = fp2PionCor4;
                break;
                
            case kKaon:
                profile2Pos = fp2KaonCor2Pos[fIndexSampling][iEtaGapIndex];
                profile2Neg = fp2KaonCor2Neg[fIndexSampling][iEtaGapIndex];
                //profile4 = fp2KaonCor4;
                break;
                
            case kProton:
                profile2Pos = fp2ProtonCor2Pos[fIndexSampling][iEtaGapIndex];
                profile2Neg = fp2ProtonCor2Neg[fIndexSampling][iEtaGapIndex];
                //profile4 = fp2ProtonCor4;
                break;
                
            default:
                ::Error("DoFlowPID","Unexpected species! Terminating!");
                return;
        }
    }
    
    if(fDoOnlyMixedCorrelations){
        
        switch (species)
        {
            case kPion:
                profile3Pos = fpMixedPionCor3Pos[fIndexSampling][iEtaGapIndex];
                profile3Neg = fpMixedPionCor3Neg[fIndexSampling][iEtaGapIndex];
                
                profile4Pos = fpMixedPionCor4Pos[fIndexSampling][iEtaGapIndex];
                profile4Neg = fpMixedPionCor4Neg[fIndexSampling][iEtaGapIndex];
                break;
                
            case kKaon:
                profile3Pos = fpMixedKaonCor3Pos[fIndexSampling][iEtaGapIndex];
                profile3Neg = fpMixedKaonCor3Neg[fIndexSampling][iEtaGapIndex];
                
                profile4Pos = fpMixedKaonCor4Pos[fIndexSampling][iEtaGapIndex];
                profile4Neg = fpMixedKaonCor4Neg[fIndexSampling][iEtaGapIndex];
                break;
                
            case kProton:
                profile3Pos = fpMixedProtonCor3Pos[fIndexSampling][iEtaGapIndex];
                profile3Neg = fpMixedProtonCor3Neg[fIndexSampling][iEtaGapIndex];
                
                profile4Pos = fpMixedProtonCor4Pos[fIndexSampling][iEtaGapIndex];
                profile4Neg = fpMixedProtonCor4Neg[fIndexSampling][iEtaGapIndex];
                break;
                
            default:
                ::Error("DoFlowPID with mixed harmonics","Unexpected species! Terminating!");
                return;
        }
    }
    
    FillPOIsVectors(iEtaGapIndex,species); // Filling POIs vectors
    
    const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
    
    //Float_t dEtaGap = fEtaGap[iEtaGapIndex];
    Short_t iHarmonics = 0;
    Double_t Dn2 = 0;
    TComplex vector = TComplex(0,0,kFALSE);
    Double_t dValue = 999;
    Double_t DDn3GapP=0;
    Double_t DDn4GapP=0;
    
    
    // filling POIs (P,S) flow vectors
    
    for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
    {
        if(!fDoOnlyMixedCorrelations){
            // estimating <2'>
            // POIs in positive eta
            Dn2 = TwoDiffGapPos(0,0,iPt).Re();
            if(Dn2 != 0)
            {
                for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
                {
                    iHarmonics = fHarmonics[iHarm];
                    vector = TwoDiffGapPos(iHarmonics,-iHarmonics,iPt);
                    dValue = vector.Re()/Dn2;
                    if( TMath::Abs(dValue < 1) )
                        profile2Pos[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
                }
            }
            // POIs in negative eta
            Dn2 = TwoDiffGapNeg(0,0,iPt).Re();
            if(Dn2 != 0)
            {
                for(Short_t iHarm(0); iHarm < fNumHarmonics; iHarm++)
                {
                    iHarmonics = fHarmonics[iHarm];
                    vector = TwoDiffGapNeg(iHarmonics,-iHarmonics,iPt);
                    dValue = vector.Re()/Dn2;
                    if( TMath::Abs(dValue < 1) )
                        profile2Neg[iHarm]->Fill(fIndexCentrality, iPt*dPtBinWidth, dValue, Dn2);
                }
            }
        }
        if(fDoOnlyMixedCorrelations){
            // estimating <3'>
            // POIs in positive eta
            DDn3GapP = ThreeDiffGapPos(0,0,0,iPt).Re();
            if(DDn3GapP!=0)
            {
                for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics-1; iMixedHarm++)
                {
                    if(iMixedHarm==0){ vector = ThreeDiffGapPos(4,-2,-2,iPt);}
                    if(iMixedHarm==1){ vector = ThreeDiffGapPos(6,-3,-3,iPt);}
                    if(iMixedHarm==2){ vector = ThreeDiffGapPos(5,-3,-2,iPt);}
                    dValue = vector.Re()/DDn3GapP;
                    if( TMath::Abs(dValue < 1) ){
                        profile3Pos[iMixedHarm]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn3GapP);
                    }
                }
            }
            // POIs in negative eta
            DDn3GapP = ThreeDiffGapNeg(0,0,0,iPt).Re();
            if(DDn3GapP!=0)
            {
                for(Short_t iMixedHarm(0); iMixedHarm < fNumMixedHarmonics-1; iMixedHarm++)
                {
                    if(iMixedHarm==0){ vector = ThreeDiffGapNeg(4,-2,-2,iPt);}
                    if(iMixedHarm==1){ vector = ThreeDiffGapNeg(6,-3,-3,iPt);}
                    if(iMixedHarm==2){ vector = ThreeDiffGapNeg(5,-2,-3,iPt);}
                    dValue = vector.Re()/DDn3GapP;
                    if( TMath::Abs(dValue < 1) ){
                        profile3Neg[iMixedHarm]->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn3GapP);
                    }
                }
            }
            
            //estimating <4'>
            //POIs in positive eta
            DDn4GapP = Four13DiffGapPos(0,0,0,0,iPt).Re();
            if(DDn4GapP!=0)
            {
                vector = Four13DiffGapPos(6,-2,-2,-2,iPt);
                dValue = vector.Re()/DDn4GapP;
                if( TMath::Abs(dValue < 1) ){
                    profile4Pos->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn4GapP);
                }
            }
            //POIs in negative eta
            DDn4GapP = Four13DiffGapNeg(0,0,0,0,iPt).Re();
            if(DDn4GapP!=0)
            {
                vector = Four13DiffGapNeg(6,-2,-2,-2,iPt);
                dValue = vector.Re()/DDn4GapP;
                if( TMath::Abs(dValue < 1) ){
                    profile4Neg->Fill(fIndexCentrality,iPt*dPtBinWidth, dValue, DDn4GapP);
                }
            }
        }
    } // endfor {iPt}
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillRefsVectors(const Short_t iEtaGapIndex)
{
    // Filling Q flow vector with RFPs
    // return kTRUE if succesfull (i.e. no error occurs), kFALSE otherwise
    // *************************************************************
    const Float_t dEtaGap = fEtaGap[iEtaGapIndex];
    TH3D* h3NUAWeights = 0x0;
    TH1F* hNUEWeights = 0x0;
    Double_t dNUAWeight = 1.;
    Double_t dNUEWeight = 1.;
    
    // clearing output (global) flow vectors
    ResetRFPsVector(fFlowVecQpos);
    ResetRFPsVector(fFlowVecQneg);
    
    Double_t dQcosPos, dQcosNeg, dQsinPos, dQsinNeg;
    
    // Double_t dQcosPos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
    // Double_t dQcosNeg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
    // Double_t dQsinPos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
    // Double_t dQsinNeg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax] = {0};
    
    for (auto part = fVectorCharged->begin(); part != fVectorCharged->end(); part++)
    {
        // checking species of used particles (just for double checking purpose)
        if( part->species != kCharged)
        {
            ::Warning("FillRefsVectors","Unexpected part. species (%d) in selected sample (expected %d)",part->species,kCharged);
            continue;
        }
        if(fFlowUseNUAWeights)
        {
            if(part->charge >0) h3NUAWeights = fh3NUAWeightRefsPlus;
            if(part->charge <0) h3NUAWeights = fh3NUAWeightRefsMinus;
            
            if(!h3NUAWeights) { ::Error("FillRefsVectors","Histogram with NUA weights not found."); continue; }
        }
        
        if(fFlowUseNUEWeights)
        {
            if(part->charge >0) hNUEWeights = fhNUEWeightRefsPlus;
            if(part->charge <0) hNUEWeights = fhNUEWeightRefsMinus;
            
            if(!hNUEWeights) { ::Error("FillRefsVectors","Histogram with NUE weights not found."); return; }
        }
        // RFPs pT check
        if(fCutFlowRFPsPtMin > 0. && part->pt < fCutFlowRFPsPtMin)
            continue;
        
        if(fCutFlowRFPsPtMax > 0. && part->pt > fCutFlowRFPsPtMax)
            continue;
        
        // 0-ing variables
        dQcosPos = 0;
        dQcosNeg = 0;
        dQsinPos = 0;
        dQsinNeg = 0;
        
        // loading weights if needed
        if(fFlowUseNUAWeights && h3NUAWeights)
        {
            dNUAWeight = h3NUAWeights->GetBinContent(h3NUAWeights->FindBin(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ()));
            if(dNUAWeight <= 0) dNUAWeight = 1.;
            //if(iEtaGapIndex == 0)
            fh3AfterNUAWeightsRefs->Fill(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ(), dNUAWeight);
        }
        if(fFlowUseNUEWeights && hNUEWeights)
        {
            dNUEWeight = hNUEWeights->GetBinContent(hNUEWeights->FindBin(part->pt));
            if(dNUEWeight <= 0) dNUEWeight = 1.;
            //if(iEtaGapIndex == 0)
            fhAfterNUEWeightsRefs->Fill(part->pt, dNUEWeight);
        }
        
        // RPF candidate passing all criteria: start filling flow vectors
        
        if(part->eta > dEtaGap / 2)
        {
            // RFP in positive eta acceptance
            for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++){
                for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
                {
                    dQcosPos = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Cos(iHarm * part->phi);
                    dQsinPos = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Sin(iHarm * part->phi);
                    fFlowVecQpos[iHarm][iPower] += TComplex(dQcosPos,dQsinPos,kFALSE);
                }
            }
        }
        if(part->eta < -dEtaGap / 2 )
        {
            // RFP in negative eta acceptance
            for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++){
                for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
                {
                    dQcosNeg = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Cos(iHarm * part->phi);
                    dQsinNeg = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Sin(iHarm * part->phi);
                    fFlowVecQneg[iHarm][iPower] += TComplex(dQcosNeg,dQsinNeg,kFALSE);
                }
            }
        }
    } // endfor {tracks} particle loop
    
    // // filling local flow vectors to global flow vector arrays
    // for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    //   for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    //   {
    //     fFlowVecQpos[iHarm][iPower] = TComplex(dQcosPos[iHarm][iPower],dQsinPos[iHarm][iPower],kFALSE);
    //     if(dEtaGap > -1)
    //       fFlowVecQneg[iHarm][iPower] = TComplex(dQcosNeg[iHarm][iPower],dQsinNeg[iHarm][iPower],kFALSE);
    //   }
    
    // printf("RFPs EtaGap %g : number %g (pos) %g (neg) \n", dEtaGap,fFlowVecQpos[0][0].Re(),fFlowVecQneg[0][0].Re());
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::FillPOIsVectors(const Short_t iEtaGapIndex, const PartSpecies species, const Short_t iMassIndex)
{
    // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
    // *************************************************************
    
    if(species == kUnknown) return;
    
    Double_t dNUAWeight = 1.;  // for generic framework != 1
    Double_t dNUEWeight = 1.;  // for generic framework != 1
    
    // clearing output (global) flow vectors
    ResetPOIsVector(fFlowVecPpos);
    ResetPOIsVector(fFlowVecPneg);
    ResetPOIsVector(fFlowVecS);
    
    std::vector<FlowPart>* vector = 0x0;
    //TH3D* hist = 0x0;
    //Double_t dMassLow = 0, dMassHigh = 0;
    TH3D* h3NUAWeights = 0x0;
    TH1F* hNUEWeights = 0x0;
    
    // swich based on species
    switch (species)
    {
        case kCharged:
            vector = fVectorCharged;
            break;
            
        case kPion:
            vector = fVectorPion;
            break;
            
        case kKaon:
            vector = fVectorKaon;
            break;
            
        case kProton:
            vector = fVectorProton;
            break;
            
        default:
            ::Error("FillPOIsVectors","Selected species unknown.");
            return;
    }//switch(species)
    
    const Double_t dEtaGap = fEtaGap[iEtaGapIndex];
    //const Double_t dMass = (dMassLow+dMassHigh)/2;
    
    Short_t iPtBin = 0;
    Double_t dCos = 0, dSin = 0;
    
    for (auto part = vector->begin(); part != vector->end(); part++)
    {
        
        // swich based on species
        switch (species)
        {
            case kCharged:
                if(fFlowUseNUAWeights) {
                    if(part->charge >0) h3NUAWeights = fh3NUAWeightChargedPlus;
                    if(part->charge <0) h3NUAWeights = fh3NUAWeightChargedMinus;
                }
                if(fFlowUseNUEWeights) {
                    if(part->charge >0) hNUEWeights = fhNUEWeightChargedPlus;
                    if(part->charge <0) hNUEWeights = fhNUEWeightChargedMinus;
                }
                break;
                
            case kPion:
                if(fFlowUseNUAWeights) {
                    if(part->charge >0) h3NUAWeights = fh3NUAWeightPionPlus;
                    if(part->charge <0) h3NUAWeights = fh3NUAWeightPionMinus;
                }
                if(fFlowUseNUEWeights) {
                    if(part->charge >0) hNUEWeights = fhNUEWeightPionPlus;
                    if(part->charge <0) hNUEWeights = fhNUEWeightPionMinus;
                }
                break;
                
            case kKaon:
                if(fFlowUseNUAWeights) {
                    if(part->charge >0) h3NUAWeights = fh3NUAWeightKaonPlus;
                    if(part->charge <0) h3NUAWeights = fh3NUAWeightKaonMinus;
                }
                if(fFlowUseNUEWeights) {
                    if(part->charge >0) hNUEWeights = fhNUEWeightKaonPlus;
                    if(part->charge <0) hNUEWeights = fhNUEWeightKaonMinus;
                }
                break;
                
            case kProton:
                if(fFlowUseNUAWeights) {
                    if(part->charge >0) h3NUAWeights = fh3NUAWeightProtonPlus;
                    if(part->charge <0) h3NUAWeights = fh3NUAWeightProtonMinus;
                }
                if(fFlowUseNUEWeights) {
                    if(part->charge >0) hNUEWeights = fhNUEWeightProtonPlus;
                    if(part->charge <0) hNUEWeights = fhNUEWeightProtonMinus;
                }
                break;
                
            default:
                ::Error("FillPOIsVectors","Selected species unknown.");
                return;
        }//switch(species)
        
        if(fFlowUseNUAWeights && !h3NUAWeights) { ::Error("FillPOIsVectors","Histogram with NUA weights not found."); continue; }
        if(fFlowUseNUEWeights && !hNUEWeights) { ::Error("FillPOIsVectors","Histogram with NUE weights not found."); return; }
        
        // checking species of used particles (just for double checking purpose)
        if( part->species != species)
        {
            ::Warning("FillPOIsVectors","Unexpected part. species (%d) in selected sample (expected %d)",part->species,species);
            continue;
        }//endif{part->species != species}
        
        // assign iPtBin based on particle momenta
        iPtBin = GetPOIsPtBinIndex(part->pt);
        
        // 0-ing variables
        dCos = 0;
        dSin = 0;
        
        // POIs candidate passing all criteria: start filling flow vectors
        
        // loading weights if needed
        if(fFlowUseNUAWeights && h3NUAWeights)
        {
            dNUAWeight = h3NUAWeights->GetBinContent(h3NUAWeights->FindBin(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ()));
            if(dNUAWeight <= 0){ dNUAWeight = 1.;}
            //if(iEtaGapIndex == 0){
            switch (species)
            {
                case kCharged:
                    fh3AfterNUAWeightsCharged->Fill(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ(), dNUAWeight);
                    break;
                case kPion:
                    fh3AfterNUAWeightsPion->Fill(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ(), dNUAWeight);
                    break;
                case kKaon:
                    fh3AfterNUAWeightsKaon->Fill(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ(), dNUAWeight);
                    break;
                case kProton:
                    fh3AfterNUAWeightsProton->Fill(part->phi,part->eta,fEventAOD->GetPrimaryVertex()->GetZ(), dNUAWeight);
                    break;
                default:
                    ::Error("Fill fh3AfterNUAWeights","Selected species unknown.");
                    return;
            }
            //}
        }//endif{fFlowUseNUAWeights && h3NUAWeights}
        
        if(fFlowUseNUEWeights && hNUEWeights)
        {
            dNUEWeight = hNUEWeights->GetBinContent(hNUEWeights->FindBin(part->pt));
            if(dNUEWeight <= 0) dNUEWeight = 1.;
            //if(iEtaGapIndex == 0){
            switch (species)
            {
                case kCharged:
                    fhAfterNUEWeightsCharged->Fill(part->pt, dNUEWeight);
                    break;
                case kPion:
                    fhAfterNUEWeightsPion->Fill(part->pt, dNUEWeight);
                    break;
                case kKaon:
                    fhAfterNUEWeightsKaon->Fill(part->pt, dNUEWeight);
                    break;
                case kProton:
                    fhAfterNUEWeightsProton->Fill(part->pt, dNUEWeight);
                    break;
                default:
                    ::Error("Fill fhAfterNUEWeights","Selected species unknown.");
                    return;
            }
            //}
        }//endif{fFlowUseNUEWeights && hNUEWeights}
        
        if(part->eta > dEtaGap / 2 )
        {
            // particle in positive eta acceptance
            for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++){
                for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
                {
                    dCos = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Cos(iHarm * part->phi);
                    dSin = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Sin(iHarm * part->phi);
                    fFlowVecPpos[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);
                }//endfor{Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++}
            }//endfor{Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++}
        }//endif{part->eta > dEtaGap / 2 }
        if(part->eta < -dEtaGap / 2 ){
            // particle in negative eta acceptance
            for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++){
                for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
                {
                    dCos = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Cos(iHarm * part->phi);
                    dSin = TMath::Power(dNUAWeight*dNUEWeight,iPower) * TMath::Sin(iHarm * part->phi);
                    fFlowVecPneg[iHarm][iPower][iPtBin] += TComplex(dCos,dSin,kFALSE);
                }//endfor{Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++}
            }//endfor{Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++}
        }//endif{part->eta < -dEtaGap / 2 }
    } // endfor {tracks}
    return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowModes_pPb::GetPOIsPtBinIndex(const Double_t pt)
{
    // Return POIs pT bin index based on pT value
    // *************************************************************
    const Double_t dPtBinWidth = (fFlowPOIsPtMax - fFlowPOIsPtMin) / fFlowPOIsPtNumBins;
    // printf("Pt %g | index %d\n",pt,(Short_t) (pt / dPtBinWidth) );
    return (Short_t) (pt / dPtBinWidth);
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::ResetRFPsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
    // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
    // *************************************************************
    for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
            array[iHarm][iPower] = TComplex(0,0,kFALSE);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::ResetPOIsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins])
{
    for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
            for(Short_t iPt(0); iPt < fFlowPOIsPtNumBins; iPt++)
                array[iHarm][iPower][iPt] = TComplex(0,0,kFALSE);
    return;
}
//_____________________________________________________________________________
void AliAnalysisTaskFlowModes_pPb::ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
    // List all values of given flow vector TComplex array
    // *************************************************************
    printf(" ### Listing (TComplex) flow vector array ###########################\n");
    for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    {
        printf("Harm %d (power):",iHarm);
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
            printf("|(%d) %g+%g(i)",iPower, array[iHarm][iPower].Re(), array[iHarm][iPower].Im());
        }
        printf("\n");
    }
    return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowModes_pPb::GetSamplingIndex()
{
    // Assessing sampling index based on generated random number
    // returns centrality index
    // *************************************************************
    
    Short_t index = 0x0;
    
    if(fSampling && fNumSamples > 1)
    {
        TRandom3 rr(0);
        Double_t ranNum = rr.Rndm(); // getting random number in (0,1)
        Double_t generated = ranNum * fNumSamples; // getting random number in range (0, fNumSamples)
        // finding right index for sampling based on generated number and total number of samples
        for(Short_t i(0); i < fNumSamples; i++)
        {
            if(generated < (i+1) )
            {
                index = i;
                break;
            }
        }
    }
    
    return index;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskFlowModes_pPb::GetCentralityIndex()
{
    // Estimating centrality percentile based on selected estimator.
    // (Default) If no multiplicity estimator is specified (fMultEstimator == '' || Charged), percentile is estimated as number of selected / filtered charged tracks.
    // If a valid multiplicity estimator is specified, centrality percentile is estimated via AliMultSelection
    // otherwise -1 is returned (and event is skipped)
    // *************************************************************
    fMultEstimator.ToUpper();
    if(
       fMultEstimator.EqualTo("V0A") || fMultEstimator.EqualTo("V0C") || fMultEstimator.EqualTo("V0M") ||
       fMultEstimator.EqualTo("CL0") || fMultEstimator.EqualTo("CL1") || fMultEstimator.EqualTo("ZNA") ||
       fMultEstimator.EqualTo("ZNC") || fMultEstimator.EqualTo("TRK") ){
        
        // some of supported AliMultSelection estimators (listed above)
        Float_t dPercentile = 300;
        
        // checking AliMultSelection
        AliMultSelection* multSelection = 0x0;
        multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
        if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1;}
        
        dPercentile = multSelection->GetMultiplicityPercentile(fMultEstimator.Data());
/*
        if(fFullCentralityRange){
            if(dPercentile > 100 || dPercentile < 0) 
            { 
                AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); 
                return -1;
            }
            else 
            { 
                return dPercentile;
            }
        }
        if(!fFullCentralityRange){
            if(dPercentile > 100 || dPercentile <50) { AliWarning("Centrality percentile estimated not within 50-100 range. Returning -1"); return -1;}
            else { return dPercentile;}
        }
*/
    if(fFullCentralityRange)
    {
        return dPercentile;
    }
    }
    else if(fMultEstimator.EqualTo("") || fMultEstimator.EqualTo("CHARGED"))
    {
        // assigning centrality based on number of selected charged tracks
        return fVectorCharged->size();
    }
    else
    {
        AliWarning(Form("Multiplicity estimator '%s' not supported. Returning -1\n",fMultEstimator.Data()));
        return -1;
    }
    
    
    return -1;
}
//____________________________
Double_t AliAnalysisTaskFlowModes_pPb::GetWDist(const AliAODVertex* v0, const AliAODVertex* v1) {
    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
        return kFALSE;
    }
    static TMatrixDSym vVb(3);
    double dist = -1;
    double dx = v0->GetX()-v1->GetX();
    double dy = v0->GetY()-v1->GetY();
    double dz = v0->GetZ()-v1->GetZ();
    double cov0[6],cov1[6];
    v0->GetCovarianceMatrix(cov0);
    v1->GetCovarianceMatrix(cov1);
    vVb(0,0) = cov0[0]+cov1[0];
    vVb(1,1) = cov0[2]+cov1[2];
    vVb(2,2) = cov0[5]+cov1[5];
    vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
    vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
    vVb.InvertFast();
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz+2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;
}
//_________________________________________________
void AliAnalysisTaskFlowModes_pPb::Terminate(Option_t* option)
{
    // called on end of task, after all events are processed
    // *************************************************************
    
    return;
}
//_____________________________________________________________________________
// Set of methods returning given complex flow vector based on flow harmonics (n) and weight power indexes (p)
// a la General Framework implementation.
// Q: flow vector of RFPs (with/out eta gap)
// P: flow vector of POIs (with/out eta gap) (in usual notation p)
// S: flow vector of overlaping RFPs and POIs (in usual notation q)

TComplex AliAnalysisTaskFlowModes_pPb::Q(const Short_t n, const Short_t p)
{
    if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
    else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::QGapPos(const Short_t n, const Short_t p)
{
    if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
    else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::QGapNeg(const Short_t n, const Short_t p)
{
    if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
    else return fFlowVecQneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::P(const Short_t n, const Short_t p, const Short_t pt)
{
    if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p][pt]);
    else return fFlowVecPpos[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::PGapPos(const Short_t n, const Short_t p, const Short_t pt)
{
    if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p][pt]);
    else return fFlowVecPpos[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::PGapNeg(const Short_t n, const Short_t p, const Short_t pt)
{
    if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p][pt]);
    else return fFlowVecPneg[n][p][pt];
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::S(const Short_t n, const Short_t p, const Short_t pt)
{
    if(n < 0) return TComplex::Conjugate(fFlowVecS[-n][p][pt]);
    else return fFlowVecS[n][p][pt];
}
//____________________________________________________________________

// Set of flow calculation methods for cumulants of different orders with/out eta gap

TComplex AliAnalysisTaskFlowModes_pPb::Two(const Short_t n1, const Short_t n2)
{
    TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::TwoGap(const Short_t n1, const Short_t n2)
{
    TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::TwoDiff(const Short_t n1, const Short_t n2, const Short_t pt)
{
    TComplex formula = P(n1,1,pt)*Q(n2,1) - S(n1+n2,1,pt);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::TwoDiffGapPos(const Short_t n1, const Short_t n2, const Short_t pt)
{
    TComplex formula = PGapPos(n1,1,pt)*QGapNeg(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::TwoDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t pt)
{
    TComplex formula = PGapNeg(n1,1,pt)*QGapPos(n2,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::ThreeDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3,const Short_t pt)
{
    TComplex formula = PGapPos(n1,1,pt)*QGapNeg(n2,1)*QGapNeg(n3,1)- PGapPos(n1,1,pt)*QGapNeg(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::ThreeDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3,const Short_t pt)
{
    TComplex formula = PGapNeg(n1,1,pt)*QGapPos(n2,1)*QGapPos(n3,1)- PGapNeg(n1,1,pt)*QGapPos(n2+n3,2);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
    + 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
    return formula;
}
// //____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt)
{
    TComplex formula = P(n1,1,pt)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2,pt)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2,pt)*Q(n4,1)
    - P(n1,1,pt)*Q(n2+n3,2)*Q(n4,1)+2.*S(n1+n2+n3,3,pt)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2,pt)
    + Q(n2+n3,2)*S(n1+n4,2,pt)-P(n1,1,pt)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2,pt)*Q(n2+n4,2)
    + 2.*Q(n3,1)*S(n1+n2+n4,3,pt)-P(n1,1,pt)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2,pt)*Q(n3+n4,2)
    + 2.*Q(n2,1)*S(n1+n3+n4,3,pt)+2.*P(n1,1,pt)*Q(n2+n3+n4,3)-6.*S(n1+n2+n3+n4,4,pt);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::FourGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4) //  n1+n2 = n3+n4;  n1, n2 from P & n3, n4 from M
{
    TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapPos(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3+n4,2)+QGapPos(n1+n2,2)*QGapNeg(n3+n4,2);
    //TComplex *out = (TComplex*) &formula;
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::FourGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4) //  n1+n2 = n3+n4;  n1, n2 from M & n3, n4 from P
{
    TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)-QGapNeg(n1+n2,2)*QGapPos(n3,1)*QGapPos(n4,1)-QGapNeg(n1,1)*QGapNeg(n2,1)*QGapPos(n3+n4,2)+QGapNeg(n1+n2,2)*QGapPos(n3+n4,2);
    //TComplex *out = (TComplex*) &formula;
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::Four13DiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4,const Short_t pt) // n1 = n2 + n3 + n4
{
    TComplex formula = PGapPos(n1,1,pt)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)- PGapPos(n1,1,pt)*QGapNeg(n2+n3,2)*QGapNeg(n4,1) - PGapPos(n1,1,pt)*QGapNeg(n2+n4,2)*QGapNeg(n3,1) - PGapPos(n1,1,pt)*QGapNeg(n3+n4,2)*QGapNeg(n2,1) + 2.*PGapPos(n1,1,pt)*QGapNeg(n2+n3+n4,3);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::Four13DiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4,const Short_t pt) // n1 = n2 + n3 + n4
{
    TComplex formula = PGapNeg(n1,1,pt)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)- PGapNeg(n1,1,pt)*QGapPos(n2+n3,2)*QGapPos(n4,1) - PGapNeg(n1,1,pt)*QGapPos(n2+n4,2)*QGapPos(n3,1) - PGapNeg(n1,1,pt)*QGapPos(n3+n4,2)*QGapPos(n2,1) + 2.*PGapNeg(n1,1,pt)*QGapPos(n2+n3+n4,3);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::SixGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t n5, const Short_t n6) // n1 + n2 + n3 = n4 + n5 + n6
{
    TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1) - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1) - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapNeg(n4+n6,2)*QGapNeg(n5,1) - QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapNeg(n5+n6,2)*QGapNeg(n4,1) + 2.*QGapPos(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapNeg(n4+n5+n6,3)- QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1) + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1) + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapNeg(n4+n6,2)*QGapNeg(n5,1) + QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapNeg(n5+n6,2)*QGapNeg(n4,1) - 2.*QGapPos(n1+n2,2)*QGapPos(n3,1)*QGapNeg(n4+n5+n6,3) - QGapPos(n1+n3,2)*QGapPos(n2,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1) + QGapPos(n1+n3,2)*QGapPos(n2,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1) + QGapPos(n1+n3,2)*QGapPos(n2,1)*QGapNeg(n4+n6,2)*QGapNeg(n5,1) + QGapPos(n1+n3,2)*QGapPos(n2,1)*QGapNeg(n5+n6,2)*QGapNeg(n4,1) - 2.*QGapPos(n1+n3,2)*QGapPos(n2,1)*QGapNeg(n4+n5+n6,3) - QGapPos(n2+n3,2)*QGapPos(n1,1)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1) + QGapPos(n2+n3,2)*QGapPos(n1,1)*QGapNeg(n4+n5,2)*QGapNeg(n6,1) + QGapPos(n2+n3,2)*QGapPos(n1,1)*QGapNeg(n4+n6,2)*QGapNeg(n5,1) + QGapPos(n2+n3,2)*QGapPos(n1,1)*QGapNeg(n5+n6,2)*QGapNeg(n4,1) - 2.*QGapPos(n2+n3,2)*QGapPos(n1,1)*QGapNeg(n4+n5+n6,3)+ 2.*QGapPos(n1+n2+n3,3)*QGapNeg(n4,1)*QGapNeg(n5,1)*QGapNeg(n6,1) - 2.*QGapPos(n1+n2+n3,3)*QGapNeg(n4+n5,2)*QGapNeg(n6,1) - 2.*QGapPos(n1+n2+n3,3)*QGapNeg(n4+n6,2)*QGapNeg(n5,1) - 2.*QGapPos(n1+n2+n3,3)*QGapNeg(n5+n6,2)*QGapNeg(n4,1) + 4.*QGapPos(n1+n2+n3,3)*QGapNeg(n4+n5+n6,3);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskFlowModes_pPb::SixGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t n5, const Short_t n6) // n1 + n2 + n3 = n4 + n5 + n6
{
    
    TComplex formula = QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1) - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapPos(n4+n5,2)*QGapPos(n6,1) - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapPos(n4+n6,2)*QGapPos(n5,1) - QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapPos(n5+n6,2)*QGapPos(n4,1) + 2.*QGapNeg(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapPos(n4+n5+n6,3)- QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1) + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4+n5,2)*QGapPos(n6,1) + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4+n6,2)*QGapPos(n5,1) + QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n5+n6,2)*QGapPos(n4,1) - 2.*QGapNeg(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4+n5+n6,3) - QGapNeg(n1+n3,2)*QGapNeg(n2,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1) + QGapNeg(n1+n3,2)*QGapNeg(n2,1)*QGapPos(n4+n5,2)*QGapPos(n6,1) + QGapNeg(n1+n3,2)*QGapNeg(n2,1)*QGapPos(n4+n6,2)*QGapPos(n5,1) + QGapNeg(n1+n3,2)*QGapNeg(n2,1)*QGapPos(n5+n6,2)*QGapPos(n4,1) - 2.*QGapNeg(n1+n3,2)*QGapNeg(n2,1)*QGapPos(n4+n5+n6,3) - QGapNeg(n2+n3,2)*QGapNeg(n1,1)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1) + QGapNeg(n2+n3,2)*QGapNeg(n1,1)*QGapPos(n4+n5,2)*QGapPos(n6,1) + QGapNeg(n2+n3,2)*QGapNeg(n1,1)*QGapPos(n4+n6,2)*QGapPos(n5,1) + QGapNeg(n2+n3,2)*QGapNeg(n1,1)*QGapPos(n5+n6,2)*QGapPos(n4,1) - 2.*QGapNeg(n2+n3,2)*QGapNeg(n1,1)*QGapPos(n4+n5+n6,3)+ 2.*QGapNeg(n1+n2+n3,3)*QGapPos(n4,1)*QGapPos(n5,1)*QGapPos(n6,1) -2.*QGapNeg(n1+n2+n3,3)*QGapPos(n4+n5,2)*QGapPos(n6,1) - 2.*QGapNeg(n1+n2+n3,3)*QGapPos(n4+n6,2)*QGapPos(n5,1) - 2.*QGapNeg(n1+n2+n3,3)*QGapPos(n5+n6,2)*QGapPos(n4,1) + 4.*QGapNeg(n1+n2+n3,3)*QGapPos(n4+n5+n6,3);
    return formula;
}

//-----------------------------------------------------------------------
void AliAnalysisTaskFlowModes_pPb::SetPriors(Float_t centrCur){
    //set priors for the bayesian pid selection
    fCurrCentr = centrCur;
    
    fBinLimitPID[0] = 0.300000;
    fBinLimitPID[1] = 0.400000;
    fBinLimitPID[2] = 0.500000;
    fBinLimitPID[3] = 0.600000;
    fBinLimitPID[4] = 0.700000;
    fBinLimitPID[5] = 0.800000;
    fBinLimitPID[6] = 0.900000;
    fBinLimitPID[7] = 1.000000;
    fBinLimitPID[8] = 1.200000;
    fBinLimitPID[9] = 1.400000;
    fBinLimitPID[10] = 1.600000;
    fBinLimitPID[11] = 1.800000;
    fBinLimitPID[12] = 2.000000;
    fBinLimitPID[13] = 2.200000;
    fBinLimitPID[14] = 2.400000;
    fBinLimitPID[15] = 2.600000;
    fBinLimitPID[16] = 2.800000;
    fBinLimitPID[17] = 3.000000;
    
    for(Int_t i=18;i<fgkPIDptBin;i++){
        fBinLimitPID[i] = 3.0 + 0.2 * (i-17);
    }
    
    // 0-10%
    if(centrCur < 10){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0168;
        fC[1][4] = 0.0122;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0272;
        fC[2][4] = 0.0070;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0562;
        fC[3][4] = 0.0258;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0861;
        fC[4][4] = 0.0496;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1168;
        fC[5][4] = 0.0740;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1476;
        fC[6][4] = 0.0998;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1810;
        fC[7][4] = 0.1296;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2240;
        fC[8][4] = 0.1827;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2812;
        fC[9][4] = 0.2699;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3328;
        fC[10][4] = 0.3714;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3780;
        fC[11][4] = 0.4810;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.4125;
        fC[12][4] = 0.5771;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.4486;
        fC[13][4] = 0.6799;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.4840;
        fC[14][4] = 0.7668;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4971;
        fC[15][4] = 0.8288;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4956;
        fC[16][4] = 0.8653;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.5173;
        fC[17][4] = 0.9059;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 10-20%
    else if(centrCur < 20){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0132;
        fC[1][4] = 0.0088;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0283;
        fC[2][4] = 0.0068;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0577;
        fC[3][4] = 0.0279;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0884;
        fC[4][4] = 0.0534;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1179;
        fC[5][4] = 0.0794;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1480;
        fC[6][4] = 0.1058;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1807;
        fC[7][4] = 0.1366;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2219;
        fC[8][4] = 0.1891;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2804;
        fC[9][4] = 0.2730;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3283;
        fC[10][4] = 0.3660;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3710;
        fC[11][4] = 0.4647;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.4093;
        fC[12][4] = 0.5566;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.4302;
        fC[13][4] = 0.6410;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.4649;
        fC[14][4] = 0.7055;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4523;
        fC[15][4] = 0.7440;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4591;
        fC[16][4] = 0.7799;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.4804;
        fC[17][4] = 0.8218;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 20-30%
    else if(centrCur < 30){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0102;
        fC[1][4] = 0.0064;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0292;
        fC[2][4] = 0.0066;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0597;
        fC[3][4] = 0.0296;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0900;
        fC[4][4] = 0.0589;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1199;
        fC[5][4] = 0.0859;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1505;
        fC[6][4] = 0.1141;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1805;
        fC[7][4] = 0.1454;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2221;
        fC[8][4] = 0.2004;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2796;
        fC[9][4] = 0.2838;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3271;
        fC[10][4] = 0.3682;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3648;
        fC[11][4] = 0.4509;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3988;
        fC[12][4] = 0.5339;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.4315;
        fC[13][4] = 0.5995;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.4548;
        fC[14][4] = 0.6612;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4744;
        fC[15][4] = 0.7060;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4899;
        fC[16][4] = 0.7388;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.4411;
        fC[17][4] = 0.7293;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 30-40%
    else if(centrCur < 40){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0102;
        fC[1][4] = 0.0048;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0306;
        fC[2][4] = 0.0079;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0617;
        fC[3][4] = 0.0338;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0920;
        fC[4][4] = 0.0652;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1211;
        fC[5][4] = 0.0955;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1496;
        fC[6][4] = 0.1242;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1807;
        fC[7][4] = 0.1576;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2195;
        fC[8][4] = 0.2097;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2732;
        fC[9][4] = 0.2884;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3204;
        fC[10][4] = 0.3679;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3564;
        fC[11][4] = 0.4449;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3791;
        fC[12][4] = 0.5052;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.4062;
        fC[13][4] = 0.5647;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.4234;
        fC[14][4] = 0.6203;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4441;
        fC[15][4] = 0.6381;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4629;
        fC[16][4] = 0.6496;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.4293;
        fC[17][4] = 0.6491;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 40-50%
    else if(centrCur < 50){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0093;
        fC[1][4] = 0.0057;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0319;
        fC[2][4] = 0.0075;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0639;
        fC[3][4] = 0.0371;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0939;
        fC[4][4] = 0.0725;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1224;
        fC[5][4] = 0.1045;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1520;
        fC[6][4] = 0.1387;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1783;
        fC[7][4] = 0.1711;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2202;
        fC[8][4] = 0.2269;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2672;
        fC[9][4] = 0.2955;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3191;
        fC[10][4] = 0.3676;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3434;
        fC[11][4] = 0.4321;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3692;
        fC[12][4] = 0.4879;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.3993;
        fC[13][4] = 0.5377;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.3818;
        fC[14][4] = 0.5547;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4003;
        fC[15][4] = 0.5484;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4281;
        fC[16][4] = 0.5383;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.3960;
        fC[17][4] = 0.5374;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 50-60%
    else if(centrCur < 60){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0076;
        fC[1][4] = 0.0032;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0329;
        fC[2][4] = 0.0085;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0653;
        fC[3][4] = 0.0423;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0923;
        fC[4][4] = 0.0813;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1219;
        fC[5][4] = 0.1161;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1519;
        fC[6][4] = 0.1520;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1763;
        fC[7][4] = 0.1858;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2178;
        fC[8][4] = 0.2385;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2618;
        fC[9][4] = 0.3070;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.3067;
        fC[10][4] = 0.3625;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3336;
        fC[11][4] = 0.4188;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3706;
        fC[12][4] = 0.4511;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.3765;
        fC[13][4] = 0.4729;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.3942;
        fC[14][4] = 0.4855;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.4051;
        fC[15][4] = 0.4762;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.3843;
        fC[16][4] = 0.4763;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.4237;
        fC[17][4] = 0.4773;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 60-70%
    else if(centrCur < 70){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0071;
        fC[1][4] = 0.0012;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0336;
        fC[2][4] = 0.0097;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0662;
        fC[3][4] = 0.0460;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0954;
        fC[4][4] = 0.0902;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1181;
        fC[5][4] = 0.1306;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1481;
        fC[6][4] = 0.1662;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1765;
        fC[7][4] = 0.1963;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2155;
        fC[8][4] = 0.2433;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2580;
        fC[9][4] = 0.3022;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.2872;
        fC[10][4] = 0.3481;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3170;
        fC[11][4] = 0.3847;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3454;
        fC[12][4] = 0.4258;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.3580;
        fC[13][4] = 0.4299;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.3903;
        fC[14][4] = 0.4326;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.3690;
        fC[15][4] = 0.4491;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.4716;
        fC[16][4] = 0.4298;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.3875;
        fC[17][4] = 0.4083;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 70-80%
    else if(centrCur < 80){
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0075;
        fC[1][4] = 0.0007;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0313;
        fC[2][4] = 0.0124;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0640;
        fC[3][4] = 0.0539;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0923;
        fC[4][4] = 0.0992;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1202;
        fC[5][4] = 0.1417;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1413;
        fC[6][4] = 0.1729;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1705;
        fC[7][4] = 0.1999;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.2103;
        fC[8][4] = 0.2472;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2373;
        fC[9][4] = 0.2916;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.2824;
        fC[10][4] = 0.3323;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.3046;
        fC[11][4] = 0.3576;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3585;
        fC[12][4] = 0.4003;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.3461;
        fC[13][4] = 0.3982;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.3362;
        fC[14][4] = 0.3776;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.3071;
        fC[15][4] = 0.3500;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.2914;
        fC[16][4] = 0.3937;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.3727;
        fC[17][4] = 0.3877;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    // 80-100%
    else{
        fC[0][0] = 0.005;
        fC[0][1] = 0.005;
        fC[0][2] = 1.0000;
        fC[0][3] = 0.0010;
        fC[0][4] = 0.0010;
        
        fC[1][0] = 0.005;
        fC[1][1] = 0.005;
        fC[1][2] = 1.0000;
        fC[1][3] = 0.0060;
        fC[1][4] = 0.0035;
        
        fC[2][0] = 0.005;
        fC[2][1] = 0.005;
        fC[2][2] = 1.0000;
        fC[2][3] = 0.0323;
        fC[2][4] = 0.0113;
        
        fC[3][0] = 0.005;
        fC[3][1] = 0.005;
        fC[3][2] = 1.0000;
        fC[3][3] = 0.0609;
        fC[3][4] = 0.0653;
        
        fC[4][0] = 0.005;
        fC[4][1] = 0.005;
        fC[4][2] = 1.0000;
        fC[4][3] = 0.0922;
        fC[4][4] = 0.1076;
        
        fC[5][0] = 0.005;
        fC[5][1] = 0.005;
        fC[5][2] = 1.0000;
        fC[5][3] = 0.1096;
        fC[5][4] = 0.1328;
        
        fC[6][0] = 0.005;
        fC[6][1] = 0.005;
        fC[6][2] = 1.0000;
        fC[6][3] = 0.1495;
        fC[6][4] = 0.1779;
        
        fC[7][0] = 0.005;
        fC[7][1] = 0.005;
        fC[7][2] = 1.0000;
        fC[7][3] = 0.1519;
        fC[7][4] = 0.1989;
        
        fC[8][0] = 0.005;
        fC[8][1] = 0.005;
        fC[8][2] = 1.0000;
        fC[8][3] = 0.1817;
        fC[8][4] = 0.2472;
        
        fC[9][0] = 0.005;
        fC[9][1] = 0.005;
        fC[9][2] = 1.0000;
        fC[9][3] = 0.2429;
        fC[9][4] = 0.2684;
        
        fC[10][0] = 0.005;
        fC[10][1] = 0.005;
        fC[10][2] = 1.0000;
        fC[10][3] = 0.2760;
        fC[10][4] = 0.3098;
        
        fC[11][0] = 0.005;
        fC[11][1] = 0.005;
        fC[11][2] = 1.0000;
        fC[11][3] = 0.2673;
        fC[11][4] = 0.3198;
        
        fC[12][0] = 0.005;
        fC[12][1] = 0.005;
        fC[12][2] = 1.0000;
        fC[12][3] = 0.3165;
        fC[12][4] = 0.3564;
        
        fC[13][0] = 0.005;
        fC[13][1] = 0.005;
        fC[13][2] = 1.0000;
        fC[13][3] = 0.3526;
        fC[13][4] = 0.3011;
        
        fC[14][0] = 0.005;
        fC[14][1] = 0.005;
        fC[14][2] = 1.0000;
        fC[14][3] = 0.3788;
        fC[14][4] = 0.3011;
        
        fC[15][0] = 0.005;
        fC[15][1] = 0.005;
        fC[15][2] = 1.0000;
        fC[15][3] = 0.3788;
        fC[15][4] = 0.3011;
        
        fC[16][0] = 0.005;
        fC[16][1] = 0.005;
        fC[16][2] = 1.0000;
        fC[16][3] = 0.3788;
        fC[16][4] = 0.3011;
        
        fC[17][0] = 0.005;
        fC[17][1] = 0.005;
        fC[17][2] = 1.0000;
        fC[17][3] = 0.3788;
        fC[17][4] = 0.3011;
        
        for(Int_t i=18;i<fgkPIDptBin;i++){
            fC[i][0] = fC[17][0];
            fC[i][1] = fC[17][1];
            fC[i][2] = fC[17][2];
            fC[i][3] = fC[17][3];
            fC[i][4] = fC[17][4];
        }
    }
    
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t AliAnalysisTaskFlowModes_pPb::TPCTOFagree(const AliVTrack *track)
{
    //check pid agreement between TPC and TOF
    Bool_t status = kFALSE;
    
    const Float_t c = 2.99792457999999984e-02;
    
    Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};
    
    
    Double_t exptimes[9];
    track->GetIntegratedTimes(exptimes);
    
    Float_t dedx = track->GetTPCsignal();
    
    Float_t p = track->P();
    Float_t time = track->GetTOFsignal()- fESDpid.GetTOFResponse().GetStartTime(p);
    Float_t tl = exptimes[0]*c; // to work both for ESD and AOD
    
    Float_t betagammares =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);
    
    Float_t betagamma1 = tl/(time-5 *betagammares) * 33.3564095198152043;
    
    //  printf("betagamma1 = %f\n",betagamma1);
    
    if(betagamma1 < 0.1) betagamma1 = 0.1;
    
    if(betagamma1 < 0.99999) betagamma1 /= TMath::Sqrt(1-betagamma1*betagamma1);
    else betagamma1 = 100;
    
    Float_t betagamma2 = tl/(time+5 *betagammares) * 33.3564095198152043;
    //  printf("betagamma2 = %f\n",betagamma2);
    
    if(betagamma2 < 0.1) betagamma2 = 0.1;
    
    if(betagamma2 < 0.99999) betagamma2 /= TMath::Sqrt(1-betagamma2*betagamma2);
    else betagamma2 = 100;
    
    
    Float_t momtpc=track->GetTPCmomentum();
    
    for(Int_t i=0;i < 5;i++){
        Float_t resolutionTOF =  fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[i], mass[i]);
        if(TMath::Abs(exptimes[i] - time) < 5 * resolutionTOF){
            Float_t dedxExp = 0;
            if(i==0) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kElectron);
            else if(i==1) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kMuon);
            else if(i==2) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kPion);
            else if(i==3) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kKaon);
            else if(i==4) dedxExp =  fESDpid.GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kProton);
            
            Float_t resolutionTPC = 2;
            if(i==0) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kElectron);
            else if(i==1) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kMuon);
            else if(i==2) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kPion);
            else if(i==3) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kKaon);
            else if(i==4) resolutionTPC =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);
            
            if(TMath::Abs(dedx - dedxExp) < 3 * resolutionTPC){
                status = kTRUE;
            }
        }
    }
    
    Float_t bb1 =  fESDpid.GetTPCResponse().Bethe(betagamma1);
    Float_t bb2 =  fESDpid.GetTPCResponse().Bethe(betagamma2);
    Float_t bbM =  fESDpid.GetTPCResponse().Bethe((betagamma1+betagamma2)*0.5);
    
    
    //  status = kFALSE;
    // for nuclei
    Float_t resolutionTOFpr =   fESDpid.GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);
    Float_t resolutionTPCpr =   fESDpid.GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);
    if(TMath::Abs(dedx-bb1) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
        status = kTRUE;
    }
    else if(TMath::Abs(dedx-bb2) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
        status = kTRUE;
    }
    else if(TMath::Abs(dedx-bbM) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
        status = kTRUE;
    }
    
    return status;
}


