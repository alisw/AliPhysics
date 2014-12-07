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
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <TAxis.h>
// #include "TObjArray.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>

#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODv0.h>
#include <AliAODInputHandler.h>

#include "AliAnalysisTaskProtonLambda.h"
#include <AliCentrality.h>
//#include "AliAODpid.h"
#include <AliPID.h>
#include <AliPIDResponse.h>
// #include <../STEER/STEER/AliV0.h>
#include <AliExternalTrackParam.h>
//#include <AliAODTrack.h>
//#include <AliESDtrack.h>

//#include "EventCollection.h"

// Task to study femtoscopic proton-lambda correlations
// Author: Hans Beck

ClassImp(AliAnalysisTaskProtonLambda)
//ClassImp(AliAnalysisTaskProtonLambda::GlobalTrackInfo)
//ClassImp(AliAnalysisTaskProtonLambda::GTIContainer)

#endif  
//________________________________________________________________________
AliAnalysisTaskProtonLambda::AliAnalysisTaskProtonLambda() 
  : AliAnalysisTaskSE(),
  fkUseOnTheFly(kTRUE),
  fkAbsZvertexCut(10.0),
  fkCentCut(20.0),
  fkLamMass(1.115683),
  fkProMass(0.9382720),
  fkPioMass(0.13957018),
  fPIDResponse(0), 
  fTpcResponse(0),
  fFemtoBuffer(0),
  fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
  fOutput2Part(0),
  fGTI(0),fTrackBuffSize(18000),
  fHistGoodEvent(0),
  // fHistPrimaryVertexPosXY(0), fHistPrimaryVertexPosZ(0),        
  // fHistTrackMultiplicity(0),    
  // fHistShareV0pos(0),fHistShareV0neg(0),
  // fHistPosTofBeforeCut(0), fHistPosTofAfterCut(0),           
  // fHistNegTofBeforeCut(0), fHistNegTofAfterCut(0),           
  // fHistPosTpcBeforeCut(0), fHistPosTpcAfterCut(0),            
  // fHistNegTpcBeforeCut(0), fHistNegTpcAfterCut(0),            
  // fHistGoodV0(0), fHistCorrectSigns(0),              
  // fHistDcaPosToPrimVertex(0), fHistDcaNegToPrimVertex(0),        
  // fHistDcaPosToPrimVertexZoom(0), fHistDcaNegToPrimVertexZoom(0),  
  // fHistRadiusV0(0), fHistDecayLengthV0(0), fHistDcaV0Daughters(0),          
  // fHistChi2(0), fHistCosPointAngle(0), fHistCosPointAngleZoom(0),
  fHistSideBandOffLam(0), fHistSideBandOffALam(0), fHistTPCNclsPosOffLam(0),       
  fHistTPCNclsNegOffLam(0), fHistTPCNclsPosOffALam(0), fHistTPCNclsNegOffALam(0),      
  // fHistPosNsigmaTpcOffLam(0), fHistPosNsigmaTpcOffALam(0), fHistNegNsigmaTpcOffLam(0),
  // fHistNegNsigmaTpcOffALam(0), fHistUseTofOffLam(0), fHistUseTofOffALam(0),
  // fHistDcaPosOffLam(0), fHistDcaPosOffALam(0), fHistDcaNegOffLam(0),           
  // fHistDcaNegOffALam(0), fHistDcaV0DaughtersOffLam(0), fHistDcaV0DaughtersOffALam(0),  
  // fHistCosPointLamOff(0), fHistCosPointALamOff(0), fHistCosPointLamZoomOff(0),     
  // fHistCosPointALamZoomOff(0), fHistV0RadiusLamOff(0), fHistV0RadiusALamOff(0),        
  // fHistV0DecayLengthLamOff(0), fHistV0DecayLengthALamOff(0), fHistDcaV0PriVertexLamOff(0),     
  // fHistDcaV0PriVertexALamOff(0),
  fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
  // fHistPtVsMassLambdaOff(0), fHistPtVsMassAntiLambdaOff(0),
  fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
  // fHistPtVsYLambdaOff(0), fHistPtVsYAntiLambdaOff(0),       
  fHistSideBandOnLam(0), fHistSideBandOnALam(0),
  // fHistLikeSignOnLam(0), fHistLikeSignOnALam(0),         
  fHistTPCNclsPosOnLam(0), fHistTPCNclsNegOnLam(0), fHistTPCNclsPosOnALam(0),fHistTPCNclsNegOnALam(0),     
  // fHistPosNsigmaTpcOnLam(0), fHistPosNsigmaTpcOnALam(0), fHistNegNsigmaTpcOnLam(0), fHistNegNsigmaTpcOnALam(0),        
  // fHistUseTofOnLam(0),fHistUseTofOnALam(0),fHistDcaPosOnLam(0),fHistDcaPosOnALam(0),fHistDcaNegOnLam(0),               
  // fHistDcaNegOnALam(0),fHistDcaV0DaughtersOnLam(0),fHistDcaV0DaughtersOnALam(0),fHistCosPointLamOn(0),             
  // fHistCosPointALamOn(0),fHistCosPointLamZoomOn(0),fHistCosPointALamZoomOn(0),fHistV0RadiusLamOn(0),             
  // fHistV0RadiusALamOn(0),fHistV0DecayLengthLamOn(0),fHistV0DecayLengthALamOn(0),fHistDcaV0PriVertexLamOn(0),       
  // fHistDcaV0PriVertexALamOn(0),
  // fHistChi2TPCPosLamOn(0),  fHistChi2TPCPosALamOn(0),  fHistChi2TPCNegLamOn(0),  fHistChi2TPCNegALamOn(0),
  // fHistMinvTPConlyLamOn(0),  fHistMinvTPConlyALamOn(0),
  fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
  // fHistPtVsMassLambdaOn(0), fHistPtVsMassAntiLambdaOn(0),
  fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
  // fHistPtVsYLambdaOn(0), fHistPtVsYAntiLambdaOn(0),
  // fHistMomDiffLam(0),fHistMomDiffALam(0),fHistMomDiffBgLam(0),fHistMomDiffBgALam(0),
  // fHistMomDiffWoSPDLam(0),fHistMomDiffWoSPDALam(0),fHistMomDiffWoSPDBgLam(0),fHistMomDiffWoSPDBgALam(0),
  fPriHistShare(0),
  // fPriHistPosNsigmaTof(0),
  fPriHistPosNsigmaTofVsP(0),fPriHistPosNsigmaTofVsPt(0),     
  // fPriHistNegNsigmaTof(0),
  fPriHistNegNsigmaTofVsP(0),fPriHistNegNsigmaTofVsPt(0),fPriHistTOFsignalPosVsP(0),      
  fPriHistTOFsignalPosVsPt(0),fPriHistTOFsignalNegVsP(0),fPriHistTOFsignalNegVsPt(0),fPriHistHybridTOFsigPosWoTPC(0), 
  fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
  // fPriHistHasTofPos(0),          
  fPriHistTPCsignalPos(0),
  // fPriHistNsigmaTPCPos(0), fPriHistTPCsignalTOFcutPos(0),fPriHistNsigmaTPCTOFcutPos(0),   
  fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
  // fPriHistHasTofNeg(0),         
  fPriHistTPCsignalNeg(0),
  // fPriHistNsigmaTPCNeg(0),fPriHistTPCsignalTOFcutNeg(0),fPriHistNsigmaTPCTOFcutNeg(0),    
  fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
  fPriHistDCAxyYPtPro(0),fPriHistDCAxyYPtAPro(0),
  // f2HistLamLamMeanMinDistProReal(0),
  // f2HistLamLamMeanMinDistPioReal(0),f2HistLamProMeanMinDistProReal(0),f2HistALamALamMeanMinDistAProReal(0), 
  // f2HistALamALamMeanMinDistPioReal(0),f2HistALamAProMeanMinDistAProReal(0),
  // f2HistSftLamLamMeanMinDistProReal(0),
  // f2HistSftLamLamMeanMinDistPioReal(0),f2HistSftLamProMeanMinDistProReal(0),f2HistSftALamALamMeanMinDistAProReal(0), 
  // f2HistSftALamALamMeanMinDistPioReal(0),f2HistSftALamAProMeanMinDistAProReal(0),
  // f2HistSftIrocLamLamMeanMinDistProReal(0),
  // f2HistSftIrocLamLamMeanMinDistPioReal(0),f2HistSftIrocLamProMeanMinDistProReal(0),f2HistSftIrocALamALamMeanMinDistAProReal(0), 
  // f2HistSftIrocALamALamMeanMinDistPioReal(0),f2HistSftIrocALamAProMeanMinDistAProReal(0),
  // f2HistSftOrocLamLamMeanMinDistProReal(0),
  // f2HistSftOrocLamLamMeanMinDistPioReal(0),f2HistSftOrocLamProMeanMinDistProReal(0),f2HistSftOrocALamALamMeanMinDistAProReal(0), 
  // f2HistSftOrocALamALamMeanMinDistPioReal(0),f2HistSftOrocALamAProMeanMinDistAProReal(0),
  // f2HistMtLamLamReal(0), 
  f2HistMtLamProReal(0), 
  // f2HistMtALamALamReal(0), 
  f2HistMtALamAProReal(0),
  // f2HistMtLowQLamLamReal(0), 
  f2HistMtLowQLamProReal(0), 
  // f2HistMtLowQALamALamReal(0), 
  f2HistMtLowQALamAProReal(0),
  LamProReal(0),ALamAProReal(0),
  // f3HistLamLamQinvReal(0),               
  // f3HistALamALamQinvReal(0),f3HistLamLamMinvReal(0),               
  // f3HistLamProMinvReal(0),f3HistALamALamMinvReal(0),f3HistALamAProMinvReal(0),
  // f2HistBgLamBgLamMeanMinDistProReal(0),f2HistBgLamBgLamMeanMinDistPioReal(0),
  // f2HistBgLamProMeanMinDistProReal(0),f2HistBgALamBgALamMeanMinDistAProReal(0),
  // f2HistBgALamBgALamMeanMinDistPioReal(0),f2HistBgALamAProMeanMinDistAProReal(0),
  // f2HistSftBgLamBgLamMeanMinDistProReal(0),f2HistSftBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftBgLamProMeanMinDistProReal(0),f2HistSftBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftBgALamBgALamMeanMinDistPioReal(0),f2HistSftBgALamAProMeanMinDistAProReal(0),
  // f2HistSftIrocBgLamBgLamMeanMinDistProReal(0),f2HistSftIrocBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftIrocBgLamProMeanMinDistProReal(0),f2HistSftIrocBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftIrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftIrocBgALamAProMeanMinDistAProReal(0),
  // f2HistSftOrocBgLamBgLamMeanMinDistProReal(0),f2HistSftOrocBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftOrocBgLamProMeanMinDistProReal(0),f2HistSftOrocBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftOrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftOrocBgALamAProMeanMinDistAProReal(0),
  BgLamProReal(0),BgALamAProReal(0),
  // f3HistBgLamBgLamQinvReal(0),             
  // f3HistBgALamBgALamQinvReal(0),
  // f2HistLamLamMeanMinDistProMixed(0),f2HistLamLamMeanMinDistPioMixed(0),
  // f2HistLamProMeanMinDistProMixed(0),f2HistALamALamMeanMinDistAProMixed(0),   
  // f2HistALamALamMeanMinDistPioMixed(0),f2HistALamAProMeanMinDistAProMixed(0),
  // f2HistSftLamLamMeanMinDistProMixed(0),f2HistSftLamLamMeanMinDistPioMixed(0),
  // f2HistSftLamProMeanMinDistProMixed(0),f2HistSftALamALamMeanMinDistAProMixed(0),   
  // f2HistSftALamALamMeanMinDistPioMixed(0),f2HistSftALamAProMeanMinDistAProMixed(0),
  // f2HistSftIrocLamLamMeanMinDistProMixed(0),f2HistSftIrocLamLamMeanMinDistPioMixed(0),
  // f2HistSftIrocLamProMeanMinDistProMixed(0),f2HistSftIrocALamALamMeanMinDistAProMixed(0),   
  // f2HistSftIrocALamALamMeanMinDistPioMixed(0),f2HistSftIrocALamAProMeanMinDistAProMixed(0),
  // f2HistSftOrocLamLamMeanMinDistProMixed(0),f2HistSftOrocLamLamMeanMinDistPioMixed(0),
  // f2HistSftOrocLamProMeanMinDistProMixed(0),f2HistSftOrocALamALamMeanMinDistAProMixed(0),   
  // f2HistSftOrocALamALamMeanMinDistPioMixed(0),f2HistSftOrocALamAProMeanMinDistAProMixed(0),
  LamProMixed(0),ALamAProMixed(0),
  // f3HistLamLamQinvMixed(0),                
  // f3HistALamALamQinvMixed(0),f3HistLamLamMinvMixed(0),                
  // f3HistLamProMinvMixed(0),f3HistALamALamMinvMixed(0),f3HistALamAProMinvMixed(0),
  // f2HistBgLamBgLamMeanMinDistProMixed(0),f2HistBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistBgLamProMeanMinDistProMixed(0),f2HistBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistBgALamBgALamMeanMinDistPioMixed(0),f2HistBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftBgLamBgLamMeanMinDistProMixed(0),f2HistSftBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftBgLamProMeanMinDistProMixed(0),f2HistSftBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftBgALamBgALamMeanMinDistPioMixed(0),f2HistSftBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftIrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftIrocBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftIrocBgLamProMeanMinDistProMixed(0),f2HistSftIrocBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftIrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftIrocBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftOrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftOrocBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftOrocBgLamProMeanMinDistProMixed(0),f2HistSftOrocBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftOrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftOrocBgALamAProMeanMinDistAProMixed(0),
  BgLamProMixed(0),BgALamAProMixed(0)
  // f3HistBgLamBgLamQinvMixed(0),             
  // f3HistBgALamBgALamQinvMixed(0)
{
  // Dummy constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::AliAnalysisTaskProtonLambda(const char *name) 
  : AliAnalysisTaskSE(name),
  fkUseOnTheFly(kTRUE),
  fkAbsZvertexCut(10.0),
  fkCentCut(20.0),
  fkLamMass(1.115683),
  fkProMass(0.9382720),
  fkPioMass(0.13957018),
  fPIDResponse(0), 
  fTpcResponse(0),
  fFemtoBuffer(0),
  fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
  fOutput2Part(0),
  fGTI(0),fTrackBuffSize(18000),
  fHistGoodEvent(0),
  // fHistPrimaryVertexPosXY(0), fHistPrimaryVertexPosZ(0),        
  // fHistTrackMultiplicity(0),    
  // fHistShareV0pos(0),fHistShareV0neg(0),
  // fHistPosTofBeforeCut(0), fHistPosTofAfterCut(0),           
  // fHistNegTofBeforeCut(0), fHistNegTofAfterCut(0),           
  // fHistPosTpcBeforeCut(0), fHistPosTpcAfterCut(0),            
  // fHistNegTpcBeforeCut(0), fHistNegTpcAfterCut(0),            
  // fHistGoodV0(0), fHistCorrectSigns(0),              
  // fHistDcaPosToPrimVertex(0), fHistDcaNegToPrimVertex(0),        
  // fHistDcaPosToPrimVertexZoom(0), fHistDcaNegToPrimVertexZoom(0),  
  // fHistRadiusV0(0), fHistDecayLengthV0(0), fHistDcaV0Daughters(0),          
  // fHistChi2(0), fHistCosPointAngle(0), fHistCosPointAngleZoom(0),
  fHistSideBandOffLam(0), fHistSideBandOffALam(0), fHistTPCNclsPosOffLam(0),       
  fHistTPCNclsNegOffLam(0), fHistTPCNclsPosOffALam(0), fHistTPCNclsNegOffALam(0),      
  // fHistPosNsigmaTpcOffLam(0), fHistPosNsigmaTpcOffALam(0), fHistNegNsigmaTpcOffLam(0),
  // fHistNegNsigmaTpcOffALam(0), fHistUseTofOffLam(0), fHistUseTofOffALam(0),
  // fHistDcaPosOffLam(0), fHistDcaPosOffALam(0), fHistDcaNegOffLam(0),           
  // fHistDcaNegOffALam(0), fHistDcaV0DaughtersOffLam(0), fHistDcaV0DaughtersOffALam(0),  
  // fHistCosPointLamOff(0), fHistCosPointALamOff(0), fHistCosPointLamZoomOff(0),     
  // fHistCosPointALamZoomOff(0), fHistV0RadiusLamOff(0), fHistV0RadiusALamOff(0),        
  // fHistV0DecayLengthLamOff(0), fHistV0DecayLengthALamOff(0), fHistDcaV0PriVertexLamOff(0),     
  // fHistDcaV0PriVertexALamOff(0),
  fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
  // fHistPtVsMassLambdaOff(0), fHistPtVsMassAntiLambdaOff(0),
  fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
  // fHistPtVsYLambdaOff(0), fHistPtVsYAntiLambdaOff(0),       
  fHistSideBandOnLam(0), fHistSideBandOnALam(0),
  // fHistLikeSignOnLam(0), fHistLikeSignOnALam(0),         
  fHistTPCNclsPosOnLam(0), fHistTPCNclsNegOnLam(0), fHistTPCNclsPosOnALam(0),fHistTPCNclsNegOnALam(0),     
  // fHistPosNsigmaTpcOnLam(0), fHistPosNsigmaTpcOnALam(0), fHistNegNsigmaTpcOnLam(0), fHistNegNsigmaTpcOnALam(0),        
  // fHistUseTofOnLam(0),fHistUseTofOnALam(0),fHistDcaPosOnLam(0),fHistDcaPosOnALam(0),fHistDcaNegOnLam(0),               
  // fHistDcaNegOnALam(0),fHistDcaV0DaughtersOnLam(0),fHistDcaV0DaughtersOnALam(0),fHistCosPointLamOn(0),             
  // fHistCosPointALamOn(0),fHistCosPointLamZoomOn(0),fHistCosPointALamZoomOn(0),fHistV0RadiusLamOn(0),             
  // fHistV0RadiusALamOn(0),fHistV0DecayLengthLamOn(0),fHistV0DecayLengthALamOn(0),fHistDcaV0PriVertexLamOn(0),       
  // fHistDcaV0PriVertexALamOn(0),
  // fHistChi2TPCPosLamOn(0),  fHistChi2TPCPosALamOn(0),  fHistChi2TPCNegLamOn(0),  fHistChi2TPCNegALamOn(0),
  // fHistMinvTPConlyLamOn(0),  fHistMinvTPConlyALamOn(0),
  fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
  // fHistPtVsMassLambdaOn(0), fHistPtVsMassAntiLambdaOn(0),
  fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
  // fHistPtVsYLambdaOn(0), fHistPtVsYAntiLambdaOn(0),
  // fHistMomDiffLam(0),fHistMomDiffALam(0),fHistMomDiffBgLam(0),fHistMomDiffBgALam(0),
  // fHistMomDiffWoSPDLam(0),fHistMomDiffWoSPDALam(0),fHistMomDiffWoSPDBgLam(0),fHistMomDiffWoSPDBgALam(0),
  fPriHistShare(0),
  // fPriHistPosNsigmaTof(0),
  fPriHistPosNsigmaTofVsP(0),fPriHistPosNsigmaTofVsPt(0),     
  // fPriHistNegNsigmaTof(0),
  fPriHistNegNsigmaTofVsP(0),fPriHistNegNsigmaTofVsPt(0),fPriHistTOFsignalPosVsP(0),      
  fPriHistTOFsignalPosVsPt(0),fPriHistTOFsignalNegVsP(0),fPriHistTOFsignalNegVsPt(0),fPriHistHybridTOFsigPosWoTPC(0), 
  fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
  // fPriHistHasTofPos(0),          
  fPriHistTPCsignalPos(0),
  // fPriHistNsigmaTPCPos(0), fPriHistTPCsignalTOFcutPos(0),fPriHistNsigmaTPCTOFcutPos(0),   
  fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
  // fPriHistHasTofNeg(0),         
  fPriHistTPCsignalNeg(0),
  // fPriHistNsigmaTPCNeg(0),fPriHistTPCsignalTOFcutNeg(0),fPriHistNsigmaTPCTOFcutNeg(0),    
  fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
  fPriHistDCAxyYPtPro(0),fPriHistDCAxyYPtAPro(0),
  // f2HistLamLamMeanMinDistProReal(0),
  // f2HistLamLamMeanMinDistPioReal(0),f2HistLamProMeanMinDistProReal(0),f2HistALamALamMeanMinDistAProReal(0), 
  // f2HistALamALamMeanMinDistPioReal(0),f2HistALamAProMeanMinDistAProReal(0),
  // f2HistSftLamLamMeanMinDistProReal(0),
  // f2HistSftLamLamMeanMinDistPioReal(0),f2HistSftLamProMeanMinDistProReal(0),f2HistSftALamALamMeanMinDistAProReal(0), 
  // f2HistSftALamALamMeanMinDistPioReal(0),f2HistSftALamAProMeanMinDistAProReal(0),
  // f2HistSftIrocLamLamMeanMinDistProReal(0),
  // f2HistSftIrocLamLamMeanMinDistPioReal(0),f2HistSftIrocLamProMeanMinDistProReal(0),f2HistSftIrocALamALamMeanMinDistAProReal(0), 
  // f2HistSftIrocALamALamMeanMinDistPioReal(0),f2HistSftIrocALamAProMeanMinDistAProReal(0),
  // f2HistSftOrocLamLamMeanMinDistProReal(0),
  // f2HistSftOrocLamLamMeanMinDistPioReal(0),f2HistSftOrocLamProMeanMinDistProReal(0),f2HistSftOrocALamALamMeanMinDistAProReal(0), 
  // f2HistSftOrocALamALamMeanMinDistPioReal(0),f2HistSftOrocALamAProMeanMinDistAProReal(0),
  // f2HistMtLamLamReal(0), 
  f2HistMtLamProReal(0), 
  // f2HistMtALamALamReal(0), 
  f2HistMtALamAProReal(0),
  // f2HistMtLowQLamLamReal(0), 
  f2HistMtLowQLamProReal(0), 
  // f2HistMtLowQALamALamReal(0), 
  f2HistMtLowQALamAProReal(0),
  LamProReal(0),ALamAProReal(0),
  // f3HistLamLamQinvReal(0),               
  // f3HistALamALamQinvReal(0),f3HistLamLamMinvReal(0),               
  // f3HistLamProMinvReal(0),f3HistALamALamMinvReal(0),f3HistALamAProMinvReal(0),
  // f2HistBgLamBgLamMeanMinDistProReal(0),f2HistBgLamBgLamMeanMinDistPioReal(0),
  // f2HistBgLamProMeanMinDistProReal(0),f2HistBgALamBgALamMeanMinDistAProReal(0),
  // f2HistBgALamBgALamMeanMinDistPioReal(0),f2HistBgALamAProMeanMinDistAProReal(0),
  // f2HistSftBgLamBgLamMeanMinDistProReal(0),f2HistSftBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftBgLamProMeanMinDistProReal(0),f2HistSftBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftBgALamBgALamMeanMinDistPioReal(0),f2HistSftBgALamAProMeanMinDistAProReal(0),
  // f2HistSftIrocBgLamBgLamMeanMinDistProReal(0),f2HistSftIrocBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftIrocBgLamProMeanMinDistProReal(0),f2HistSftIrocBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftIrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftIrocBgALamAProMeanMinDistAProReal(0),
  // f2HistSftOrocBgLamBgLamMeanMinDistProReal(0),f2HistSftOrocBgLamBgLamMeanMinDistPioReal(0),
  // f2HistSftOrocBgLamProMeanMinDistProReal(0),f2HistSftOrocBgALamBgALamMeanMinDistAProReal(0),
  // f2HistSftOrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftOrocBgALamAProMeanMinDistAProReal(0),
  BgLamProReal(0),BgALamAProReal(0),
  // f3HistBgLamBgLamQinvReal(0),             
  // f3HistBgALamBgALamQinvReal(0),
  // f2HistLamLamMeanMinDistProMixed(0),f2HistLamLamMeanMinDistPioMixed(0),
  // f2HistLamProMeanMinDistProMixed(0),f2HistALamALamMeanMinDistAProMixed(0),   
  // f2HistALamALamMeanMinDistPioMixed(0),f2HistALamAProMeanMinDistAProMixed(0),
  // f2HistSftLamLamMeanMinDistProMixed(0),f2HistSftLamLamMeanMinDistPioMixed(0),
  // f2HistSftLamProMeanMinDistProMixed(0),f2HistSftALamALamMeanMinDistAProMixed(0),   
  // f2HistSftALamALamMeanMinDistPioMixed(0),f2HistSftALamAProMeanMinDistAProMixed(0),
  // f2HistSftIrocLamLamMeanMinDistProMixed(0),f2HistSftIrocLamLamMeanMinDistPioMixed(0),
  // f2HistSftIrocLamProMeanMinDistProMixed(0),f2HistSftIrocALamALamMeanMinDistAProMixed(0),   
  // f2HistSftIrocALamALamMeanMinDistPioMixed(0),f2HistSftIrocALamAProMeanMinDistAProMixed(0),
  // f2HistSftOrocLamLamMeanMinDistProMixed(0),f2HistSftOrocLamLamMeanMinDistPioMixed(0),
  // f2HistSftOrocLamProMeanMinDistProMixed(0),f2HistSftOrocALamALamMeanMinDistAProMixed(0),   
  // f2HistSftOrocALamALamMeanMinDistPioMixed(0),f2HistSftOrocALamAProMeanMinDistAProMixed(0),
  LamProMixed(0),ALamAProMixed(0),
  // f3HistLamLamQinvMixed(0),                
  // f3HistALamALamQinvMixed(0),f3HistLamLamMinvMixed(0),                
  // f3HistLamProMinvMixed(0),f3HistALamALamMinvMixed(0),f3HistALamAProMinvMixed(0),
  // f2HistBgLamBgLamMeanMinDistProMixed(0),f2HistBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistBgLamProMeanMinDistProMixed(0),f2HistBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistBgALamBgALamMeanMinDistPioMixed(0),f2HistBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftBgLamBgLamMeanMinDistProMixed(0),f2HistSftBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftBgLamProMeanMinDistProMixed(0),f2HistSftBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftBgALamBgALamMeanMinDistPioMixed(0),f2HistSftBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftIrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftIrocBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftIrocBgLamProMeanMinDistProMixed(0),f2HistSftIrocBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftIrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftIrocBgALamAProMeanMinDistAProMixed(0),
  // f2HistSftOrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftOrocBgLamBgLamMeanMinDistPioMixed(0),
  // f2HistSftOrocBgLamProMeanMinDistProMixed(0),f2HistSftOrocBgALamBgALamMeanMinDistAProMixed(0),
  // f2HistSftOrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftOrocBgALamAProMeanMinDistAProMixed(0),
  BgLamProMixed(0),BgALamAProMixed(0)
  // f3HistBgLamBgLamQinvMixed(0),             
  // f3HistBgALamBgALamQinvMixed(0)
{
  // Constructor
  fPrimaryVtxPosition[0]=0;
  fPrimaryVtxPosition[1]=0;
  fPrimaryVtxPosition[2]=0;

  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::~AliAnalysisTaskProtonLambda() {
  // Destructor, go through the data member and delete them

  // fPIDResponse is just a pointer to the pid response task,
  // we don't create it so we don't delete it. It comes from 
  // the AliInputEventHandler

  if (fTpcResponse){
    delete fTpcResponse;
    fTpcResponse=0;
  }
  if(fFemtoBuffer){
    delete fFemtoBuffer;
    fFemtoBuffer=0;
  }
  // fAOD also just comes from a function of the AliAnalysisTaskSE
  // fPrimaryVtx comes from the fAOD

  // The lists containing the histograms
  if (fOutputList){
    fOutputList->Delete();
    delete fOutputList;
    fOutputList=0;
  }
  if (fOutputPrimaries){
    fOutputPrimaries->Delete();
    delete fOutputPrimaries;
    fOutputPrimaries=0;
  }
  if (fOutput2Part){
    fOutput2Part->Delete();
    delete fOutput2Part;
    fOutput2Part=0;
  }

  // Array, note the [] with the delete
  if (fGTI)
    delete[] fGTI;
  fGTI=0;

}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::UserCreateOutputObjects()
{
  // Create histograms and other objects and variables
  // Called once

  // Get the PID response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man){AliError("Couldn't get the analysis manager!");}
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler){AliError("Couldn't get the input handler!");}
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse){AliError("Couldn't get the PID response task!");}

  // Create dE/dx spectra cut. use it by calling
  // fTpcResponse->GetExpectedSignal(mom, AliPID::kProton)
  fTpcResponse = new AliTPCPIDResponse();
  Double_t AlephParameters[5];
  // They are only valid for data, see $ALICE_ROOT/PWG2/SPECTRA/AliProtonAnalysisBase.cxx
  // for monte carlo parameters
  AlephParameters[0] = 0.0283086;
  AlephParameters[1] = 2.63394e+01;
  AlephParameters[2] = 5.04114e-11;
  AlephParameters[3] = 2.12543e+00;
  AlephParameters[4] = 4.88663e+00;
  fTpcResponse->SetBetheBlochParameters(AlephParameters[0],AlephParameters[1],AlephParameters[2],AlephParameters[3],AlephParameters[4]);
  
  // Create the buffer for event mixing
  // Standard values are
  //  fkZvertexBins(10),
  //  fkCentBins(10),
  //  fkMixBuff(5),
  //  fkPriTrackLim(100),
  //  fkV0Lim(50),
  //  fFemtoBuffer = new FemtoBuffer(10,10,5,100,50,fkAbsZvertexCut,fkCentCut);
  fFemtoBuffer = new FemtoBuffer(4,2,7,100,50,fkAbsZvertexCut,fkCentCut);

  // In AODs, TPC only tracks don't have the pid information stored.
  // Also, the TPC only tracks don't have any resolution in the DCAxy
  // to distinguish between primaries and secondaries so we need the
  // corresponding global track for every TPC only track. The way to do 
  // this is to just store the pointer to the global track for every id.
  fGTI = new AliAODTrack *[fTrackBuffSize]; // Array of pointers 

  // Create the output list
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputPrimaries = new TList();
  fOutputPrimaries->SetOwner();
  fOutput2Part = new TList();
  fOutput2Part->SetOwner();

  // Invariant mass binning for lambdas
  const Int_t nMinvBins = 140;
  const Float_t minvLowEdge=1.060683, minvHiEdge=1.200683;

  // Control hist for event cuts
  fHistGoodEvent = new TH1F("h1GoodEvent","No of events passing the cuts.",10,-.5,9.5);
  fOutputList->Add(fHistGoodEvent);

  // Primary Vertex:
  // fHistPrimaryVertexPosXY       = new TH2F("h2PrimaryVertexPosXY", "Primary Vertex Position XY;Primary Vertex Position X (cm);Primary Vertex Position Y (cm)",100,-0.5,0.5,100,-0.5,0.5);
  // fOutputList->Add(fHistPrimaryVertexPosXY);
  // fHistPrimaryVertexPosZ       = new TH1F("h1PrimaryVertexPosZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-12.0,12.0);
  // fOutputList->Add(fHistPrimaryVertexPosZ);

  // // Multiplicity
  // fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000);
  // fOutputList->Add(fHistTrackMultiplicity);

  // //
  // //  V0 histograms
  // //
  // // Shared clusters
  // fHistShareV0pos = new TH1F("h1ShareV0pos","Shared clusters pos V0 daughters;#shared clusters;counts"
  // 			     ,160,0,160);
  // fOutputList->Add(fHistShareV0pos);
  // fHistShareV0neg = new TH1F("h1ShareV0neg","Shared clusters neg V0 daughters;#shared clusters;counts"
  // 			     ,160,0,160);
  // fOutputList->Add(fHistShareV0neg);
  
  // // PID via TPC dE/dx
  // fHistPosTpcBeforeCut = new TH2F ("h2PosTpcBeforeCut","TPC signal (pos daughter) before cut;momentum (GeV/c);TPC signal",40,0,4,100,0,400);
  // fOutputList->Add(fHistPosTpcBeforeCut);
  // fHistPosTpcAfterCut  = new TH2F ("h2PosTpcAfterCut","TPC signal (pos daughter) after cut;momentum (GeV/c);TPC signal",40,0,4,100,0,400);
  // fOutputList->Add(fHistPosTpcAfterCut);
  // fHistNegTpcBeforeCut = new TH2F ("h2NegTpcBeforeCut","TPC signal (neg daughter) before cut;momentum (GeV/c);TPC signal",40,0,4,100,0,400);
  // fOutputList->Add(fHistNegTpcBeforeCut);
  // fHistNegTpcAfterCut  = new TH2F ("h2NegTpcAfterCut","TPC signal (neg daughter) after cut;momentum (GeV/c);TPC signal",40,0,4,100,0,400);
  // fOutputList->Add(fHistNegTpcAfterCut);

  // // Histograms comparing offline and on-the-fly
  // fHistGoodV0                  = new TH2F("h2GoodV0","0: all, 1: two daughters, 2: like-sign, 3: 80 clusters4: tpcrefit;id;Status",10,-.5,9.5,2,-.5,1.5);
  // fOutputList->Add(fHistGoodV0);
  // fHistCorrectSigns            = new TH2F ("h2CorrectSigns","0: correct, 1: swapped, 2: like-sign;sign;Status",3,-.5,2.5,2,-.5,1.5);
  // fOutputList->Add(fHistCorrectSigns);
  // fHistDcaPosToPrimVertex      = new TH2F("h2DcaPosToPrimVertex", "Positive V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  // fOutputList->Add(fHistDcaPosToPrimVertex);
  // fHistDcaNegToPrimVertex      = new TH2F("h2DcaNegToPrimVertex", "Negative V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  // fOutputList->Add(fHistDcaNegToPrimVertex);
  // fHistDcaPosToPrimVertexZoom  = new TH2F("h2DcaPosToPrimVertexZoom", "Positive V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  // fOutputList->Add(fHistDcaPosToPrimVertexZoom);
  // fHistDcaNegToPrimVertexZoom  = new TH2F("h2DcaNegToPrimVertexZoom", "Negative V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  // fOutputList->Add(fHistDcaNegToPrimVertexZoom);
  // fHistRadiusV0                = new TH2F("h2RadiusV0", "Radius;Radius(cm);Status",1000,0,100,2,-0.5,1.5);
  // fOutputList->Add(fHistRadiusV0);
  // fHistDecayLengthV0           = new TH2F("h2DecayLengthV0", "V0s decay Length;decay length(cm);Status", 200, 0, 100,2,-0.5,1.5);
  // fOutputList->Add(fHistDecayLengthV0);
  // fHistDcaV0Daughters          = new TH2F("h2DcaV0Daughters", "DCA between daughters;dca(cm);Status", 160, 0, 4,2,-0.5,1.5);
  // fOutputList->Add(fHistDcaV0Daughters);
  // fHistChi2                    = new TH2F("h2Chi2", "V0s chi2;chi2;Status", 12, 0, 1.2,2,-0.5,1.5);
  // fOutputList->Add(fHistChi2);
  // fHistCosPointAngle           = new TH2F("h2CosPointAngle", "Cosine of V0's pointing angle", 100,0,1,2,-0.5,1.5);
  // fOutputList->Add(fHistCosPointAngle);
  // fHistCosPointAngleZoom       = new TH2F("h2CosPointAngleZoom", "Cosine of V0's pointing angle", 100,0.9,1,2,-0.5,1.5);
  // fOutputList->Add(fHistCosPointAngleZoom);

  //
  // V0 offline distributons
  //
  
  // Invariant mass distribution for the side band background
  fHistSideBandOffLam = new TH1F ("h1SideBandOffLam","m_{inv}(#Lambda) w/o any cuts;m_{inv}(#Lambda)",nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistSideBandOffLam);
  fHistSideBandOffALam  = new TH1F ("h1SideBandOffALam","m_{inv}(#bar{#Lambda}) w/o any cuts;m_{inv}(#bar{#Lambda})",nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistSideBandOffALam);

  // Number of TPC clusters
  fHistTPCNclsPosOffLam = new TH2F ("h2MassLamOffTPCNclsPos","m_{inv}(#Lambda) vs NTPCcls(pos);NTPCcls(pos);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsPosOffLam);
  fHistTPCNclsNegOffLam = new TH2F ("h2MassLamOffTPCNclsNeg","m_{inv}(#Lambda) vs NTPCcls(neg);NTPCcls(neg);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsNegOffLam);
  fHistTPCNclsPosOffALam = new TH2F ("h2MassALamOffTPCNclsPos","m_{inv}(#bar{#Lambda}) vs NTPCcls(pos);NTPCcls(pos);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsPosOffALam);
  fHistTPCNclsNegOffALam = new TH2F ("h2MassALamOffTPCNclsNeg","m_{inv}(#bar{#Lambda}) vs NTPCcls(neg);NTPCcls(neg);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsNegOffALam);

  // fHistPosNsigmaTpcOffLam  = new TH2F ("h2PosNsigmaTpcOffLam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistPosNsigmaTpcOffLam);
  // fHistPosNsigmaTpcOffALam = new TH2F ("h2PosNsigmaTpcOffALam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
  // fOutputList->Add(fHistPosNsigmaTpcOffALam);
  // fHistNegNsigmaTpcOffLam  = new TH2F ("h2NegNsigmaTpcOffLam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
  // fOutputList->Add(fHistNegNsigmaTpcOffLam);
  // fHistNegNsigmaTpcOffALam = new TH2F ("h2NegNsigmaTpcOffALam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
  // fOutputList->Add(fHistNegNsigmaTpcOffALam);
  // fHistUseTofOffLam   = new TH2F ("h2UseTofOffLam","0: no tof or within 5sigma, 1: tof more than 5sigma;m_{inv};TOF",2,-.5,1.5,nMinvBins,minvLowEdge,minvHiEdge);        
  // fOutputList->Add(fHistUseTofOffLam);
  // fHistUseTofOffALam  = new TH2F ("h2UseTofOffALam","0: no tof or within 5sigma, 1: tof more than 5sigma;m_{inv};TOF",2,-.5,1.5,nMinvBins,minvLowEdge,minvHiEdge);        
  // fOutputList->Add(fHistUseTofOffALam);

  // // DCA of daughters to primary vertex
  // fHistDcaPosOffLam = new TH2F ("h2DcaPosOffLam","m_{inv}(#Lambda) vs dca pos daughter;dca (cm);m_{inv}(p#pi^{-})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaPosOffLam);
  // fHistDcaPosOffALam = new TH2F ("h2DcaPosOffALam","m_{inv}(#bar{#Lambda}) vs dca pos daughter;dca (cm);m_{inv}(#bar{p}#pi^{+})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaPosOffALam);
  // fHistDcaNegOffLam = new TH2F ("h2DcaNegOffLam","m_{inv}(#Lambda) vs dca neg daughter;dca (cm);m_{inv}(p#pi^{-})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaNegOffLam);
  // fHistDcaNegOffALam = new TH2F ("h2DcaNegOffALam","m_{inv}(#bar{#Lambda}) vs dca neg daughter;dca (cm);m_{inv}(#bar{p}#pi^{+})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaNegOffALam);

  // // DCA of V0 daughters
  // fHistDcaV0DaughtersOffLam = new TH2F ("h2DcaLamDaughtersOff","DCA of #Lambda daughters vs minv;dca(cm);minv",20,0,2,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaV0DaughtersOffLam);
  // fHistDcaV0DaughtersOffALam = new TH2F ("h2DcaALamDaughtersOff","DCA of #bar{#Lambda} daughters vs minv;dca(cm);minv",20,0,2,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaV0DaughtersOffALam);

  // // Cosine of pointing angle
  // fHistCosPointLamOff = new TH2F ("h2CosPointLamOff","m_{inv}(#Lambda) vs cos(pointing angle);cos(pointing angle);m_{inv}(#Lambda)",10,0.99,1.0,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistCosPointLamOff);
  // fHistCosPointALamOff = new TH2F ("h2CosPointALamOff","m_{inv}(#bar{#Lambda}) vs cos(pointing angle);cos(pointing angle);m_{inv}(#bar{#Lambda})",10,0.99,1.0,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistCosPointALamOff);
  // fHistCosPointLamZoomOff = new TH2F ("h2CosPointLamZoomOff","m_{inv}(#Lambda) vs cos(pointing angle);cos(pointing angle);m_{inv}(#Lambda)",10,0.999,1.0,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistCosPointLamZoomOff);
  // fHistCosPointALamZoomOff = new TH2F ("h2CosPointALamZoomOff","m_{inv}(#bar{#Lambda}) vs cos(pointing angle);cos(pointing angle);m_{inv}(#bar{#Lambda})",10,0.999,1.0,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistCosPointALamZoomOff);

  // // Radius of V0 vertex position
  // fHistV0RadiusLamOff = new TH2F ("h2V0RadiusLamOff","m_{inv}(#Lambda) vs V0 radius of V0 vertex;radius(cm);m_{inv}",20,0,10,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistV0RadiusLamOff);
  // fHistV0RadiusALamOff = new TH2F ("h2V0RadiusALamOff","m_{inv}(#bar{#Lambda}) vs V0 radius of V0 vertex;radius(cm);m_{inv}",20,0,10,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistV0RadiusALamOff);

  // // Decay length of V0
  // fHistV0DecayLengthLamOff = new TH2F ("h2V0DecayLengthLamOff","m_{inv}(#Lambda) vs decay length of V0;decay length (cm);m_{inv}(#Lambda)",100,0,20,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistV0DecayLengthLamOff);
  // fHistV0DecayLengthALamOff = new TH2F ("h2V0DecayLengthALamOff","m_{inv}(#bar{#Lambda}) vs decay length of V0;decay length (cm);m_{inv}(#bar{#Lambda})",100,0,20,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistV0DecayLengthALamOff);

  // // DCA of primary vertex and V0
  // fHistDcaV0PriVertexLamOff = new TH2F ("h2DcaV0PriVertexLamOff","m_{inv}(#Lambda) vs dca (V0 - prim. vertex);dca(cm);m_{inv}(#Lambda)",200,0,20,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaV0PriVertexLamOff);
  // fHistDcaV0PriVertexALamOff = new TH2F ("h2DcaV0PriVertexALamOff","m_{inv}(#bar{#Lambda}) vs dca (V0 - prim. vertex);dca(cm);m_{inv}(#bar{#Lambda})",200,0,20,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistDcaV0PriVertexALamOff);

  // Invariant mass, invariant mass vs pt and y-pt
  fHistMassLambdaOff            = new TH1F("h1MassLambdaOff", "#Lambda^{0} Offline candidates;M(p#pi^{-}) (GeV/c^{2});Counts", nMinvBins, minvLowEdge, minvHiEdge);
  fOutputList->Add(fHistMassLambdaOff);
  fHistMassAntiLambdaOff          = new TH1F("h1MassAntiLambdaOff", "#bar{#Lambda}^{0} Offline candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", nMinvBins, minvLowEdge, minvHiEdge);
  fOutputList->Add(fHistMassAntiLambdaOff);
  // fHistPtVsMassLambdaOff         = new TH2F("h2PtVsMassLambdaOff","#Lambda^{0} Offline candidates;p_{t} (GeV/c);M(p#pi^{-}) (GeV/c^{2})",100,0,10,nMinvBins, minvLowEdge, minvHiEdge);
  // fOutputList->Add(fHistPtVsMassLambdaOff);
  // fHistPtVsMassAntiLambdaOff     = new TH2F("h2PtVsMassAntiLambdaOff","#bar{#Lambda}^{0} Offline candidates;p_{t} (GeV/c);M(#bar{p}#pi^{+}) (GeV/c^{2})",100,0,10,nMinvBins, minvLowEdge, minvHiEdge);
  // fOutputList->Add(fHistPtVsMassAntiLambdaOff);
  // fHistPtVsYLambdaOff          = new TH2F("h2PtVsYLambdaOff", "#Lambda^{0} Offline candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  // fOutputList->Add(fHistPtVsYLambdaOff);
  // fHistPtVsYAntiLambdaOff      = new TH2F("h2PtVsYAntiLambdaOff", "#bar{#Lambda}^{0} Offline candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  // fOutputList->Add(fHistPtVsYAntiLambdaOff);

    // 3d y pt mass
  fHistYPtMassLamOff = new TH3F ("h3YPtMassLamOff","m_{inv}(#Lambda) vs y and pt;y;pt;mass",30,-1.5,1.5,30,0,15,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistYPtMassLamOff);
  fHistYPtMassALamOff = new TH3F ("h3YPtMassALamOff","m_{inv}(#bar{#Lambda}) vs y and pt;y;pt;mass",30,-1.5,1.5,30,0,15,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistYPtMassALamOff);

  //
  // V0 on-the-fly distributons
  //
  
  // Invariant mass distribution for the side band background
  fHistSideBandOnLam = new TH1F ("h1SideBandOnLam","m_{inv}(#Lambda) w/o any cuts;m_{inv}(#Lambda)",nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistSideBandOnLam);
  fHistSideBandOnALam  = new TH1F ("h1SideBandOnALam","m_{inv}(#bar{#Lambda}) w/o any cuts;m_{inv}(#bar{#Lambda})",nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistSideBandOnALam);
  
  // // Like-sign
  // fHistLikeSignOnLam = new TH2F ("h2MassLamOnLikeSign"," 0 = ok, 1 = swapped, 2 = like sign;sign;m_{inv} p#pi^{-}",3,-.5,2.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistLikeSignOnLam);
  // fHistLikeSignOnALam = new TH2F ("h2MassALamOnLikeSign"," 0 = ok, 1 = swapped, 2= like sign;sign;m_{inv} #bar{p}#pi^{+}",3,-.5,2.5,nMinvBins,minvLowEdge,minvHiEdge);
  // fOutputList->Add(fHistLikeSignOnALam);

  // Number of TPC clusters
  fHistTPCNclsPosOnLam = new TH2F ("h2MassLamOnTPCNclsPos","m_{inv}(#Lambda) vs NTPCcls(pos);NTPCcls(pos);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsPosOnLam);
  fHistTPCNclsNegOnLam = new TH2F ("h2MassLamOnTPCNclsNeg","m_{inv}(#Lambda) vs NTPCcls(neg);NTPCcls(neg);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsNegOnLam);
  fHistTPCNclsPosOnALam = new TH2F ("h2MassALamOnTPCNclsPos","m_{inv}(#bar{#Lambda}) vs NTPCcls(pos);NTPCcls(pos);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsPosOnALam);
  fHistTPCNclsNegOnALam = new TH2F ("h2MassALamOnTPCNclsNeg","m_{inv}(#bar{#Lambda}) vs NTPCcls(neg);NTPCcls(neg);minv",18,0,180,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistTPCNclsNegOnALam);

//  fHistPosNsigmaTpcOnLam  = new TH2F ("h2PosNsigmaTpcOnLam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistPosNsigmaTpcOnLam);
//   fHistPosNsigmaTpcOnALam = new TH2F ("h2PosNsigmaTpcOnALam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
//   fOutputList->Add(fHistPosNsigmaTpcOnALam);
//   fHistNegNsigmaTpcOnLam  = new TH2F ("h2NegNsigmaTpcOnLam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
//   fOutputList->Add(fHistNegNsigmaTpcOnLam);
//   fHistNegNsigmaTpcOnALam = new TH2F ("h2NegNsigmaTpcOnALam","minv vs nsigma tpc pos;Nsigma TPC;m_{inv}",50,0,5,nMinvBins,minvLowEdge,minvHiEdge);        
//   fOutputList->Add(fHistNegNsigmaTpcOnALam);
//   fHistUseTofOnLam   = new TH2F ("h2UseTofOnLam","0: no tof or within 5sigma, 1: tof more than 5sigma;m_{inv};TOF",2,-.5,1.5,nMinvBins,minvLowEdge,minvHiEdge);        
//   fOutputList->Add(fHistUseTofOnLam);
//   fHistUseTofOnALam  = new TH2F ("h2UseTofOnALam","0: no tof or within 5sigma, 1: tof more than 5sigma;m_{inv};TOF",2,-.5,1.5,nMinvBins,minvLowEdge,minvHiEdge);        
//   fOutputList->Add(fHistUseTofOnALam);

//   // DCA of daughters to primary vertex
//   fHistDcaPosOnLam = new TH2F ("h2DcaPosOnLam","m_{inv}(#Lambda) vs dca pos daughter;dca (cm);m_{inv}(p#pi^{-})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaPosOnLam);
//   fHistDcaPosOnALam = new TH2F ("h2DcaPosOnALam","m_{inv}(#bar{#Lambda}) vs dca pos daughter;dca (cm);m_{inv}(#bar{p}#pi^{+})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaPosOnALam);
//   fHistDcaNegOnLam = new TH2F ("h2DcaNegOnLam","m_{inv}(#Lambda) vs dca neg daughter;dca (cm);m_{inv}(p#pi^{-})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaNegOnLam);
//   fHistDcaNegOnALam = new TH2F ("h2DcaNegOnALam","m_{inv}(#bar{#Lambda}) vs dca neg daughter;dca (cm);m_{inv}(#bar{p}#pi^{+})",50,0,0.5,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaNegOnALam);

//   // DCA of V0 daughters
//   fHistDcaV0DaughtersOnLam = new TH2F ("h2DcaLamDaughtersOn","DCA of #Lambda daughters vs minv;dca(cm);minv",20,0,2,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaV0DaughtersOnLam);
//   fHistDcaV0DaughtersOnALam = new TH2F ("h2DcaALamDaughtersOn","DCA of #bar{#Lambda} daughters vs minv;dca(cm);minv",20,0,2,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaV0DaughtersOnALam);

// // Cosine of pointing angle
//   fHistCosPointLamOn = new TH2F ("h2CosPointLamOn","m_{inv}(#Lambda) vs cos(pointing angle);cos(pointing angle);m_{inv}(#Lambda)",10,0.99,1.0,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistCosPointLamOn);
//   fHistCosPointALamOn = new TH2F ("h2CosPointALamOn","m_{inv}(#bar{#Lambda}) vs cos(pointing angle);cos(pointing angle);m_{inv}(#bar{#Lambda})",10,0.99,1.0,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistCosPointALamOn);
//   fHistCosPointLamZoomOn = new TH2F ("h2CosPointLamZoomOn","m_{inv}(#Lambda) vs cos(pointing angle);cos(pointing angle);m_{inv}(#Lambda)",10,0.999,1.0,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistCosPointLamZoomOn);
//   fHistCosPointALamZoomOn = new TH2F ("h2CosPointALamZoomOn","m_{inv}(#bar{#Lambda}) vs cos(pointing angle);cos(pointing angle);m_{inv}(#bar{#Lambda})",10,0.999,1.0,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistCosPointALamZoomOn);

//   // Radius of V0 vertex position
//   fHistV0RadiusLamOn = new TH2F ("h2V0RadiusLamOn","m_{inv}(#Lambda) vs V0 radius of V0 vertex;radius(cm);m_{inv}",20,0,10,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistV0RadiusLamOn);
//   fHistV0RadiusALamOn = new TH2F ("h2V0RadiusALamOn","m_{inv}(#bar{#Lambda}) vs V0 radius of V0 vertex;radius(cm);m_{inv}",20,0,10,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistV0RadiusALamOn);

//   // Decay length of V0
//   fHistV0DecayLengthLamOn = new TH2F ("h2V0DecayLengthLamOn","m_{inv}(#Lambda) vs decay length of V0;decay length (cm);m_{inv}(#Lambda)",100,0,20,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistV0DecayLengthLamOn);
//   fHistV0DecayLengthALamOn = new TH2F ("h2V0DecayLengthALamOn","m_{inv}(#bar{#Lambda}) vs decay length of V0;decay length (cm);m_{inv}(#bar{#Lambda})",100,0,20,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistV0DecayLengthALamOn);

//   // DCA of primary vertex and V0
//   fHistDcaV0PriVertexLamOn = new TH2F ("h2DcaV0PriVertexLamOn","m_{inv}(#Lambda) vs dca (V0 - prim. vertex);dca(cm);m_{inv}(#Lambda)",200,0,20,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaV0PriVertexLamOn);
//   fHistDcaV0PriVertexALamOn = new TH2F ("h2DcaV0PriVertexALamOn","m_{inv}(#bar{#Lambda}) vs dca (V0 - prim. vertex);dca(cm);m_{inv}(#bar{#Lambda})",200,0,20,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistDcaV0PriVertexALamOn);

//   // Chi2 of TPC clusters
//   fHistChi2TPCPosLamOn = new TH2F ("h2Chi2TPCPosLamOn","m_{inv}(#Lambda) vs TPC #Chi^{2} / ndf pos daughter;TPC #Chi^{2}/ndf pos daughter;m_{inv}(#Lambda)",100,0.,10.,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistChi2TPCPosLamOn);
//   fHistChi2TPCPosALamOn = new TH2F ("h2Chi2TPCPosALamOn","m_{inv}(#bar{#Lambda}) vs TPC #Chi^{2} / ndf pos daughter;TPC #Chi^{2}/ndf pos daughter;m_{inv}(#bar{#Lambda})",100,0.,10.,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistChi2TPCPosALamOn);
//   fHistChi2TPCNegLamOn = new TH2F ("h2Chi2TPCNegLamOn","m_{inv}(#Lambda) vs TPC #Chi^{2} / ndf neg daughter;TPC #Chi^{2}/ndf neg daughter;m_{inv}(#Lambda)",100,0.,10.,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistChi2TPCNegLamOn);
//   fHistChi2TPCNegALamOn = new TH2F ("h2Chi2TPCNegALamOn","m_{inv}(#bar{#Lambda}) vs TPC #Chi^{2} / ndf neg daughter;TPC #Chi^{2}/ndf neg daughter;m_{inv}(#bar{#Lambda})",100,0.,10.,nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistChi2TPCNegALamOn);
//   // Invariant mass with TPC only tracks
//   fHistMinvTPConlyLamOn = new TH1F ("h1MinvTPConlyLamOn","m_{inv}(#Lambda) using TPC only daughters;m_{inv}(p#pi^{-})[GeV/c^{2}]",nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistMinvTPConlyLamOn);
//   fHistMinvTPConlyALamOn = new TH1F ("h1MinvTPConlyALamOn","m_{inv}(#bar{#Lambda}) using TPC only daughters;m_{inv}(#bar{p}#pi^{+})[GeV/c^{2}]",nMinvBins,minvLowEdge,minvHiEdge);
//   fOutputList->Add(fHistMinvTPConlyALamOn);

  // Invariant mass, invariant mass vs pt and y-pt
  fHistMassLambdaOn            = new TH1F("h1MassLambdaOn", "#Lambda^{0} Online candidates;M(p#pi^{-}) (GeV/c^{2});Counts", nMinvBins, minvLowEdge, minvHiEdge);
  fOutputList->Add(fHistMassLambdaOn);
  fHistMassAntiLambdaOn          = new TH1F("h1MassAntiLambdaOn", "#bar{#Lambda}^{0} Online candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", nMinvBins, minvLowEdge, minvHiEdge);
  fOutputList->Add(fHistMassAntiLambdaOn);
  // fHistPtVsMassLambdaOn         = new TH2F("h2PtVsMassLambdaOn","#Lambda^{0} Online candidates;p_{t} (GeV/c);M(p#pi^{-}) (GeV/c^{2})",100,0,10,nMinvBins, minvLowEdge, minvHiEdge);
  // fOutputList->Add(fHistPtVsMassLambdaOn);
  // fHistPtVsMassAntiLambdaOn     = new TH2F("h2PtVsMassAntiLambdaOn","#bar{#Lambda}^{0} Online candidates;p_{t} (GeV/c);M(#bar{p}#pi^{+}) (GeV/c^{2})",100,0,10,nMinvBins, minvLowEdge, minvHiEdge);
  // fOutputList->Add(fHistPtVsMassAntiLambdaOn);
  // fHistPtVsYLambdaOn          = new TH2F("h2PtVsYLambdaOn", "#Lambda^{0} Online candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  // fOutputList->Add(fHistPtVsYLambdaOn);
  // fHistPtVsYAntiLambdaOn      = new TH2F("h2PtVsYAntiLambdaOn", "#bar{#Lambda}^{0} Online candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  // fOutputList->Add(fHistPtVsYAntiLambdaOn);

  // 3d y pt mass
  fHistYPtMassLamOn = new TH3F ("h3YPtMassLamOn","m_{inv}(#Lambda) vs y and pt;y;pt;mass",30,-1.5,1.5,30,0,15,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistYPtMassLamOn);
  fHistYPtMassALamOn = new TH3F ("h3YPtMassALamOn","m_{inv}(#bar{#Lambda}) vs y and pt;y;pt;mass",30,-1.5,1.5,30,0,15,nMinvBins,minvLowEdge,minvHiEdge);
  fOutputList->Add(fHistYPtMassALamOn);

  // // Momentum difference of standard (on-the-fly/offline) V0 and TPC only V0
  // Int_t nPBins=200; Float_t AbsPRange=1.;
  // fHistMomDiffLam = new TH3F ("h3MomDiffLam","momentum difference #DeltaP standard V0 / TPConly V0 #Lambda;#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 			      ,nPBins,-AbsPRange,AbsPRange
  // 			      ,nPBins,-AbsPRange,AbsPRange
  // 			      ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffLam);
  // fHistMomDiffALam = new TH3F ("h3MomDiffALam","momentum difference #DeltaP standard V0 / TPConly V0 #bar{#Lamdba};#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 			       ,nPBins,-AbsPRange,AbsPRange
  // 			       ,nPBins,-AbsPRange,AbsPRange
  // 			       ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffALam);
  // fHistMomDiffBgLam = new TH3F ("h3MomDiffBgLam","momentum difference #DeltaP standard V0 / TPConly V0 Bg#Lambda;#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 				,nPBins,-AbsPRange,AbsPRange
  // 				,nPBins,-AbsPRange,AbsPRange
  // 				,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffBgLam);
  // fHistMomDiffBgALam = new TH3F ("h3MomDiffBgALam","momentum difference #DeltaP standard V0 / TPConly V0 Bg#bar{#Lambda};#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 				 ,nPBins,-AbsPRange,AbsPRange
  // 				 ,nPBins,-AbsPRange,AbsPRange
  // 				 ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffBgALam);

  // // Same momentum difference with rejecting tracks with SPD hits
  // fHistMomDiffWoSPDLam = new TH3F ("h3MomDiffWoSPDLam","momentum difference #DeltaP standard V0 / TPConly V0 #Lambda;#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 			      ,nPBins,-AbsPRange,AbsPRange
  // 			      ,nPBins,-AbsPRange,AbsPRange
  // 			      ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffWoSPDLam);
  // fHistMomDiffWoSPDALam = new TH3F ("h3MomDiffWoSPDALam","momentum difference #DeltaP standard V0 / TPConly V0 #bar{#Lamdba};#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 			       ,nPBins,-AbsPRange,AbsPRange
  // 			       ,nPBins,-AbsPRange,AbsPRange
  // 			       ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffWoSPDALam);
  // fHistMomDiffWoSPDBgLam = new TH3F ("h3MomDiffWoSPDBgLam","momentum difference #DeltaP standard V0 / TPConly V0 Bg#Lambda;#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 				,nPBins,-AbsPRange,AbsPRange
  // 				,nPBins,-AbsPRange,AbsPRange
  // 				,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffWoSPDBgLam);
  // fHistMomDiffWoSPDBgALam = new TH3F ("h3MomDiffWoSPDBgALam","momentum difference #DeltaP standard V0 / TPConly V0 Bg#bar{#Lambda};#DeltaP_{x}[GeV/c];#DeltaP_{y}[GeV/c];#DeltaP_{z}[GeV/c]"
  // 				 ,nPBins,-AbsPRange,AbsPRange
  // 				 ,nPBins,-AbsPRange,AbsPRange
  // 				 ,nPBins,-AbsPRange,AbsPRange);
  // fOutputList->Add(fHistMomDiffWoSPDBgALam);

  //
  // Distributions for the primaries 
  //
  // Shared clusters
  fPriHistShare = new TH1F ("h1PriShare","Shared clusters, primaries;#shared clusters;counts",
			    160,0,160);
  fOutputPrimaries->Add(fPriHistShare);

  // Nsigma TOF distribution when forcing TOF
  // fPriHistPosNsigmaTof = new TH1F ("h1PosNsigmaTof","Nsigma distribution for positives;n_{#sigma,TOF}(proton);counts",200,-50,50);
  // fOutputPrimaries->Add(fPriHistPosNsigmaTof);
  fPriHistPosNsigmaTofVsP = new TH2F ("h2PosNsigmaTofVsP","Nsigma distribution for positives;total momentum [GeV/c];n_{#sigma,TOF}(proton);counts",20,0,5,200,-50,50);
  fOutputPrimaries->Add(fPriHistPosNsigmaTofVsP);
  fPriHistPosNsigmaTofVsPt   = new TH2F ("h2PosNsigmaTofVsPt","Nsigma distribution for positives;transverse momentum [GeV/c];n_{#sigma,TOF}(proton);counts",20,0,5,200,-50,50);
  fOutputPrimaries->Add(fPriHistPosNsigmaTofVsPt);

  // fPriHistNegNsigmaTof = new TH1F ("h1NegNsigmaTof","Nsigma distribution for negatives;n_{#sigma,TOF}(anti-proton);counts",200,-50,50);
  // fOutputPrimaries->Add(fPriHistNegNsigmaTof);
  fPriHistNegNsigmaTofVsP = new TH2F ("h2NegNsigmaTofVsP","Nsigma distribution for negatives;total momentum [GeV/c];n_{#sigma,TOF}(anti-proton);counts",20,0,5,200,-50,50);
  fOutputPrimaries->Add(fPriHistNegNsigmaTofVsP);
  fPriHistNegNsigmaTofVsPt   = new TH2F ("h2NegNsigmaTofVsPt","Nsigma distribution for negatives;transverse momentum [GeV/c];n_{#sigma,TOF}(anti-proton);counts",20,0,5,200,-50,50);
  fOutputPrimaries->Add(fPriHistNegNsigmaTofVsPt);
  fPriHistTOFsignalPosVsP = new TH2F ("h2TOFsignalPosVsP","tof signal vs p (positives);p [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]",20,0.0,5.0,120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistTOFsignalPosVsP);
  fPriHistTOFsignalPosVsPt = new TH2F ("h2TOFsignalPosVsPt","tof signal vs pt (positives);pt [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]",20,0.0,5.0,120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistTOFsignalPosVsPt);
  fPriHistTOFsignalNegVsP = new TH2F ("h2TOFsignalNegVsP","tof signal vs p (negatives);p [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]",20,0.0,5.0,120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistTOFsignalNegVsP);
  fPriHistTOFsignalNegVsPt = new TH2F ("h2TOFsignalNegVsPt","tof signal vs pt (negatives);pt [GeV/c];t_{meas} - t_{0} - t_{expected} [ps]",20,0.0,5.0,120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistTOFsignalNegVsPt);
  // Hybrid analysis
  fPriHistHybridTOFsigPosWoTPC = new TH1F ("h1HybridTOFsigPosWoTPC","tof signal pos (p=.75-1.0GeV) w/o dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]",120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistHybridTOFsigPosWoTPC);
  fPriHistHybridTOFsigPosTPCok = new TH1F ("h1HybridTOFsigPosTPCok","tof signal pos (p=.75-1.0GeV) with dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]",120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistHybridTOFsigPosTPCok);
  fPriHistHybridTOFsigNegWoTPC = new TH1F ("h1HybridTOFsigNegWoTPC","tof signal neg (p=.75-1.0GeV) w/o dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]",120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistHybridTOFsigNegWoTPC);
  fPriHistHybridTOFsigNegTPCok = new TH1F ("h1HybridTOFsigNegTPCok","tof signal neg (p=.75-1.0GeV) with dedx sel.;t_{meas} - t_{0} - t_{expected} [ps]",120,-10000.0,5000.0);
  fOutputPrimaries->Add(fPriHistHybridTOFsigNegTPCok);
  // dEdx analysis
  // fPriHistHasTofPos = new TH1F ("h1HasTofPos","Positives: 0 = no TOF, 1 = TOFpid bit there",2,-.5,1.5);
  // fOutputPrimaries->Add(fPriHistHasTofPos);
  fPriHistTPCsignalPos = new TH2F ("h2TPCsignalPos","TPC signal for positives;p_{tot};dEdx",40,0,4,100,0,400);
  fOutputPrimaries->Add(fPriHistTPCsignalPos);
  // fPriHistNsigmaTPCPos = new TH2F ("h2NsigmaTPCPos","Nsigma TPC for positives;p_{tot};N_{#sigma}",40,0,4,100,-5.0,5.0);
  // fOutputPrimaries->Add(fPriHistNsigmaTPCPos);
  // fPriHistTPCsignalTOFcutPos = new TH2F ("h2TPCsignalTOFcutPos","TPC signal for positives using a +/- 10sigma TOF cut;p_{tot};dEdx",40,0.0,4.0,100,0.0,400.0);
  // fOutputPrimaries->Add(fPriHistTPCsignalTOFcutPos);
  // fPriHistNsigmaTPCTOFcutPos = new TH2F ("h2NsigmaTPCTOFcutPos","Nsigma TPC for positives using a +/- 10sigma TOF cut;p_{tot};N_{#sigma}",40,0.0,4.0,100,-5.0,5.0);
  // fOutputPrimaries->Add(fPriHistNsigmaTPCTOFcutPos);

  // fPriHistHasTofNeg = new TH1F ("h1HasTofNeg","Negatives: 0 = no TOF, 1 = TOFpid bit there",2,-.5,1.5);
  // fOutputPrimaries->Add(fPriHistHasTofNeg);
  fPriHistTPCsignalNeg = new TH2F ("h2TPCsignalNeg","TPC signal for negatives;p_{tot};dEdx",40,0.0,4.0,100,0.0,400.0);
  fOutputPrimaries->Add(fPriHistTPCsignalNeg);
  // fPriHistNsigmaTPCNeg = new TH2F ("h2NsigmaTPCNeg","Nsigma TPC for negatives;p_{tot};N_{#sigma}",40,0.0,4.0,100,-5.0,5.0);
  // fOutputPrimaries->Add(fPriHistNsigmaTPCNeg);
  // fPriHistTPCsignalTOFcutNeg = new TH2F ("h2TPCsignalTOFcutNeg","TPC signal for negatives using a +/- 10sigma TOF cut;p_{tot};dEdx",40,0.0,4.0,100,0.0,400.0);
  // fOutputPrimaries->Add(fPriHistTPCsignalTOFcutNeg);
  // fPriHistNsigmaTPCTOFcutNeg = new TH2F ("h2NsigmaTPCTOFcutNeg","Nsigma TPC for negatives using a +/- 10sigma TOF cut;p_{tot};N_{#sigma}",40,0.0,4.0,100,-5.0,5.0);
  // fOutputPrimaries->Add(fPriHistNsigmaTPCTOFcutNeg);  

  fPriHistTPCsignalLowPPos = new TH2F ("h2TPCsignalLowPPos","dEdx for low momenta, positives",20,0.1,0.3,3000,0,3000);
  fOutputPrimaries->Add(fPriHistTPCsignalLowPPos);
  fPriHistTPCsignalMedPPos = new TH2F ("h2TPCsignalMedPPos","dEdx for medium momenta, positives",60,0.3,0.9,500,0,500);
  fOutputPrimaries->Add(fPriHistTPCsignalMedPPos);
  fPriHistTPCsignalHigPPos = new TH2F ("h2TPCsignalHigPPos","dEdx for high momenta, positives",100,0.9,1.9,120,0,120);
  fOutputPrimaries->Add(fPriHistTPCsignalHigPPos);
  fPriHistTPCsignalLowPNeg = new TH2F ("h2TPCsignalLowPNeg","dEdx for low momenta, negatives",20,0.1,0.3,3000,0,3000);
  fOutputPrimaries->Add(fPriHistTPCsignalLowPNeg);
  fPriHistTPCsignalMedPNeg = new TH2F ("h2TPCsignalMedPNeg","dEdx for medium momenta, negatives",60,0.3,0.9,500,0,500);
  fOutputPrimaries->Add(fPriHistTPCsignalMedPNeg);
  fPriHistTPCsignalHigPNeg = new TH2F ("h2TPCsignalHigPNeg","dEdx for high momenta, negatives",100,0.9,1.9,120,0,120);
  fOutputPrimaries->Add(fPriHistTPCsignalHigPNeg);
  
  //  Common for all protons

  // DCA xy distribution to determine primaries, secondaries from weak decay and secondaries from material
  fPriHistDCAxyYPtPro = new TH3F ("h3DCAxyYPtPro","DCAxy vs (y,pt) protons",100,-3.,3.,30,-1.5,1.5,14,0.,3.5);
  fOutputPrimaries->Add(fPriHistDCAxyYPtPro);
  fPriHistDCAxyYPtAPro = new TH3F ("h3DCAxyYPtAPro","DCAxy vs (y,pt) anti-protons",100,-3.,3.,30,-1.5,1.5,14,0.,3.5);
  fOutputPrimaries->Add(fPriHistDCAxyYPtAPro);

  //  2 particle histograms fOutput2Part
  // Common binning for TTR
  Int_t nDistBins=200;
  Float_t distLow=0.,distHig=20.;
  //  Two-track resolution: real events
  // f2HistLamLamMeanMinDistProReal = new TH2F ("h2LamLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistLamLamMeanMinDistProReal);
  // f2HistLamLamMeanMinDistPioReal = new TH2F ("h2LamLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistLamLamMeanMinDistPioReal);
  // f2HistLamProMeanMinDistProReal = new TH2F ("h2LamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistLamProMeanMinDistProReal);
  // f2HistALamALamMeanMinDistAProReal = new TH2F ("h2ALamALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistALamALamMeanMinDistAProReal);
  // f2HistALamALamMeanMinDistPioReal = new TH2F ("h2ALamALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistALamALamMeanMinDistPioReal);
  // f2HistALamAProMeanMinDistAProReal = new TH2F ("h2ALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistALamAProMeanMinDistAProReal);

  // f2HistSftLamLamMeanMinDistProReal = new TH2F ("h2SftLamLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftLamLamMeanMinDistProReal);
  // f2HistSftLamLamMeanMinDistPioReal = new TH2F ("h2SftLamLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftLamLamMeanMinDistPioReal);
  // f2HistSftLamProMeanMinDistProReal = new TH2F ("h2SftLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftLamProMeanMinDistProReal);
  // f2HistSftALamALamMeanMinDistAProReal = new TH2F ("h2SftALamALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftALamALamMeanMinDistAProReal);
  // f2HistSftALamALamMeanMinDistPioReal = new TH2F ("h2SftALamALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftALamALamMeanMinDistPioReal);
  // f2HistSftALamAProMeanMinDistAProReal = new TH2F ("h2SftALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftALamAProMeanMinDistAProReal);

  // f2HistSftIrocLamLamMeanMinDistProReal = new TH2F ("h2SftIrocLamLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocLamLamMeanMinDistProReal);
  // f2HistSftIrocLamLamMeanMinDistPioReal = new TH2F ("h2SftIrocLamLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocLamLamMeanMinDistPioReal);
  // f2HistSftIrocLamProMeanMinDistProReal = new TH2F ("h2SftIrocLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocLamProMeanMinDistProReal);
  // f2HistSftIrocALamALamMeanMinDistAProReal = new TH2F ("h2SftIrocALamALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocALamALamMeanMinDistAProReal);
  // f2HistSftIrocALamALamMeanMinDistPioReal = new TH2F ("h2SftIrocALamALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocALamALamMeanMinDistPioReal);
  // f2HistSftIrocALamAProMeanMinDistAProReal = new TH2F ("h2SftIrocALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocALamAProMeanMinDistAProReal);

  // f2HistSftOrocLamLamMeanMinDistProReal = new TH2F ("h2SftOrocLamLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocLamLamMeanMinDistProReal);
  // f2HistSftOrocLamLamMeanMinDistPioReal = new TH2F ("h2SftOrocLamLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocLamLamMeanMinDistPioReal);
  // f2HistSftOrocLamProMeanMinDistProReal = new TH2F ("h2SftOrocLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocLamProMeanMinDistProReal);
  // f2HistSftOrocALamALamMeanMinDistAProReal = new TH2F ("h2SftOrocALamALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocALamALamMeanMinDistAProReal);
  // f2HistSftOrocALamALamMeanMinDistPioReal = new TH2F ("h2SftOrocALamALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocALamALamMeanMinDistPioReal);
  // f2HistSftOrocALamAProMeanMinDistAProReal = new TH2F ("h2SftOrocALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocALamAProMeanMinDistAProReal);
  
  // Mt of the pairs
  Int_t nMtBins=25;
  Float_t mtLow=1.0,mtHig=3.5;
  // f2HistMtLamLamReal = new TH1F("h1MtLamLamReal"
  // 				,"m_{t}(#Lambda #Lambda);m{t} [GeV];counts"
  // 				,nMtBins,mtLow,mtHig);
  // fOutput2Part->Add(f2HistMtLamLamReal);
  f2HistMtLamProReal = new TH1F("h1MtLamProReal"
				,"m_{t}(p #Lambda);m{t} [GeV];counts"
				,nMtBins,mtLow,mtHig);
  fOutput2Part->Add(f2HistMtLamProReal);
  // f2HistMtALamALamReal = new TH1F("h1MtALamALamReal"
  // 				  ,"m_{t}(#bar{#Lambda} #bar{#Lambda});m{t} [GeV];counts"
  // 				  ,nMtBins,mtLow,mtHig);
  // fOutput2Part->Add(f2HistMtALamALamReal);
  f2HistMtALamAProReal  = new TH1F("h1MtALamAProReal"
				   ,"m_{t}(#bar{#Lambda} #bar{#Lambda});m{t} [GeV];counts"
				   ,nMtBins,mtLow,mtHig);
  fOutput2Part->Add(f2HistMtALamAProReal);
  // The same only filling for low q pairs
  // f2HistMtLowQLamLamReal = new TH1F("h1MtLowQLamLamReal"
  // 				,"m_{t}(#Lambda #Lambda);m{t} [GeV];counts"
  // 				,nMtBins,mtLow,mtHig);
  // fOutput2Part->Add(f2HistMtLowQLamLamReal);
  f2HistMtLowQLamProReal = new TH1F("h1MtLowQLamProReal"
				,"m_{t}(p #Lambda);m{t} [GeV];counts"
				,nMtBins,mtLow,mtHig);
  fOutput2Part->Add(f2HistMtLowQLamProReal);
  // f2HistMtLowQALamALamReal = new TH1F("h1MtLowQALamALamReal"
  // 				  ,"m_{t}(#bar{#Lambda} #bar{#Lambda});m{t} [GeV];counts"
  // 				  ,nMtBins,mtLow,mtHig);
  // fOutput2Part->Add(f2HistMtLowQALamALamReal);
  f2HistMtLowQALamAProReal  = new TH1F("h1MtLowQALamAProReal"
				   ,"m_{t}(#bar{#Lambda} #bar{#Lambda});m{t} [GeV];counts"
				   ,nMtBins,mtLow,mtHig);
  fOutput2Part->Add(f2HistMtLowQALamAProReal);

  // Common qinv binning
  Int_t nQinvBins = 400; // also for minv
  Float_t QinvLow = 0.0;
  Float_t QinvHig = 2.5;

  // Sept'12 Use a THnSparse for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  Int_t HnSpBins[4]={nQinvBins,nDistBins,nDistBins,nQinvBins};
  Double_t HnSpMin[4]={QinvLow,distLow,distLow,QinvLow};
  Double_t HnSpMax[4]={QinvHig,distHig,distHig,QinvHig};
  LamProReal = new THnSparseF("HnSp4LamProReal","lamProRealQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(LamProReal);
  ALamAProReal = new THnSparseF("HnSp4ALamAProReal","alamAProRealQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(ALamAProReal);

  // Qinv: real events
  // Since March 4th 2012 do corr. fcts vs distances
  // f3HistLamLamQinvReal = new TH3F ("h3LamLamQinvReal", "Qinv LamLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 				   ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamLamQinvReal);
  // f3HistALamALamQinvReal = new TH3F ("h3ALamALamQinvReal", "Qinv ALamALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min (d) #pi [cm]"
  // 				     ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamALamQinvReal);
  // // minv (H2 dibaryon??)
  // f3HistLamLamMinvReal = new TH3F ("h3LamLamMinvReal", "Minv LamLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 				   ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamLamMinvReal);
  // f3HistLamProMinvReal = new TH3F ("h3LamProMinvReal", "Minv LamPro;q_{inv} [GeV/c];<d> p [cm];min(d) p [cm]"
  // 				   ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamProMinvReal);
  // f3HistALamALamMinvReal = new TH3F ("h3ALamALamMinvReal", "Minv ALamALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min(d) #pi [cm]"
  // 				     ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamALamMinvReal);
  // f3HistALamAProMinvReal = new TH3F ("h3ALamAProMinvReal", "Minv ALamAPro;q_{inv} [GeV/c];<d> #bar{p} [cm];min(d) #bar{p} [cm]"
  // 				     ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamAProMinvReal);
  
  
  // // Two-track resolution: mixed events
  // f2HistLamLamMeanMinDistProMixed = new TH2F ("h2LamLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistLamLamMeanMinDistProMixed);
  // f2HistLamLamMeanMinDistPioMixed = new TH2F ("h2LamLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistLamLamMeanMinDistPioMixed);
  // f2HistLamProMeanMinDistProMixed = new TH2F ("h2LamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistLamProMeanMinDistProMixed);
  // f2HistALamALamMeanMinDistAProMixed = new TH2F ("h2ALamALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistALamALamMeanMinDistAProMixed);
  // f2HistALamALamMeanMinDistPioMixed = new TH2F ("h2ALamALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistALamALamMeanMinDistPioMixed);
  // f2HistALamAProMeanMinDistAProMixed = new TH2F ("h2ALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistALamAProMeanMinDistAProMixed);

  // f2HistSftLamLamMeanMinDistProMixed = new TH2F ("h2SftLamLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftLamLamMeanMinDistProMixed);
  // f2HistSftLamLamMeanMinDistPioMixed = new TH2F ("h2SftLamLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftLamLamMeanMinDistPioMixed);
  // f2HistSftLamProMeanMinDistProMixed = new TH2F ("h2SftLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftLamProMeanMinDistProMixed);
  // f2HistSftALamALamMeanMinDistAProMixed = new TH2F ("h2SftALamALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftALamALamMeanMinDistAProMixed);
  // f2HistSftALamALamMeanMinDistPioMixed = new TH2F ("h2SftALamALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftALamALamMeanMinDistPioMixed);
  // f2HistSftALamAProMeanMinDistAProMixed = new TH2F ("h2SftALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftALamAProMeanMinDistAProMixed);

  // f2HistSftIrocLamLamMeanMinDistProMixed = new TH2F ("h2SftIrocLamLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocLamLamMeanMinDistProMixed);
  // f2HistSftIrocLamLamMeanMinDistPioMixed = new TH2F ("h2SftIrocLamLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocLamLamMeanMinDistPioMixed);
  // f2HistSftIrocLamProMeanMinDistProMixed = new TH2F ("h2SftIrocLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocLamProMeanMinDistProMixed);
  // f2HistSftIrocALamALamMeanMinDistAProMixed = new TH2F ("h2SftIrocALamALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocALamALamMeanMinDistAProMixed);
  // f2HistSftIrocALamALamMeanMinDistPioMixed = new TH2F ("h2SftIrocALamALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocALamALamMeanMinDistPioMixed);
  // f2HistSftIrocALamAProMeanMinDistAProMixed = new TH2F ("h2SftIrocALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocALamAProMeanMinDistAProMixed);

  // f2HistSftOrocLamLamMeanMinDistProMixed = new TH2F ("h2SftOrocLamLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocLamLamMeanMinDistProMixed);
  // f2HistSftOrocLamLamMeanMinDistPioMixed = new TH2F ("h2SftOrocLamLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocLamLamMeanMinDistPioMixed);
  // f2HistSftOrocLamProMeanMinDistProMixed = new TH2F ("h2SftOrocLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocLamProMeanMinDistProMixed);
  // f2HistSftOrocALamALamMeanMinDistAProMixed = new TH2F ("h2SftOrocALamALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocALamALamMeanMinDistAProMixed);
  // f2HistSftOrocALamALamMeanMinDistPioMixed = new TH2F ("h2SftOrocALamALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocALamALamMeanMinDistPioMixed);
  // f2HistSftOrocALamAProMeanMinDistAProMixed = new TH2F ("h2SftOrocALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocALamAProMeanMinDistAProMixed);

  // Sept'12 Use a THnSparse for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  LamProMixed = new THnSparseF("HnSp4LamProMixed","lamProMixedQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(LamProMixed);
  ALamAProMixed = new THnSparseF("HnSp4ALamAProMixed","alamAProMixedQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(ALamAProMixed);

  // // Qinv: mixed events
  // f3HistLamLamQinvMixed = new TH3F ("h3LamLamQinvMixed", "Qinv LamLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 				    ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamLamQinvMixed);
  // f3HistALamALamQinvMixed = new TH3F ("h3ALamALamQinvMixed", "Qinv ALamALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min(d) #pi [cm]"
  // 				      ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamALamQinvMixed);
  // // minv (H2 di-baryon??)
  // f3HistLamLamMinvMixed = new TH3F ("h3LamLamMinvMixed", "Minv LamLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 				    ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamLamMinvMixed);
  // f3HistLamProMinvMixed = new TH3F ("h3LamProMinvMixed", "Minv LamPro;q_{inv} [GeV/c];<d> p [cm];min(d) p [cm]"
  // 				    ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistLamProMinvMixed);
  // f3HistALamALamMinvMixed = new TH3F ("h3ALamALamMinvMixed", "Minv ALamALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min(d) #pi [cm]"
  // 				      ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamALamMinvMixed);
  // f3HistALamAProMinvMixed = new TH3F ("h3ALamAProMinvMixed", "Minv ALamAPro;q_{inv} [GeV/c];<d> #bar{p} [cm];min(d) #bar{p} [cm]"
  // 				      ,nQinvBins,2.0,3.0,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistALamAProMinvMixed);

  // Same for Background (anti-)lambdas

  // // Two-track resolution: real events
  // f2HistBgLamBgLamMeanMinDistProReal = new TH2F ("h2BgLamBgLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgLamBgLamMeanMinDistProReal);
  // f2HistBgLamBgLamMeanMinDistPioReal = new TH2F ("h2BgLamBgLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgLamBgLamMeanMinDistPioReal);
  // f2HistBgLamProMeanMinDistProReal = new TH2F ("h2BgLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistBgLamProMeanMinDistProReal);
  // f2HistBgALamBgALamMeanMinDistAProReal = new TH2F ("h2BgALamBgALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistBgALamBgALamMeanMinDistAProReal);
  // f2HistBgALamBgALamMeanMinDistPioReal = new TH2F ("h2BgALamBgALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgALamBgALamMeanMinDistPioReal);
  // f2HistBgALamAProMeanMinDistAProReal = new TH2F ("h2BgALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgALamAProMeanMinDistAProReal);

  // f2HistSftBgLamBgLamMeanMinDistProReal = new TH2F ("h2SftBgLamBgLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgLamBgLamMeanMinDistProReal);
  // f2HistSftBgLamBgLamMeanMinDistPioReal = new TH2F ("h2SftBgLamBgLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgLamBgLamMeanMinDistPioReal);
  // f2HistSftBgLamProMeanMinDistProReal = new TH2F ("h2SftBgLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftBgLamProMeanMinDistProReal);
  // f2HistSftBgALamBgALamMeanMinDistAProReal = new TH2F ("h2SftBgALamBgALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftBgALamBgALamMeanMinDistAProReal);
  // f2HistSftBgALamBgALamMeanMinDistPioReal = new TH2F ("h2SftBgALamBgALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgALamBgALamMeanMinDistPioReal);
  // f2HistSftBgALamAProMeanMinDistAProReal = new TH2F ("h2SftBgALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgALamAProMeanMinDistAProReal);

  // f2HistSftIrocBgLamBgLamMeanMinDistProReal = new TH2F ("h2SftIrocBgLamBgLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgLamBgLamMeanMinDistProReal);
  // f2HistSftIrocBgLamBgLamMeanMinDistPioReal = new TH2F ("h2SftIrocBgLamBgLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgLamBgLamMeanMinDistPioReal);
  // f2HistSftIrocBgLamProMeanMinDistProReal = new TH2F ("h2SftIrocBgLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocBgLamProMeanMinDistProReal);
  // f2HistSftIrocBgALamBgALamMeanMinDistAProReal = new TH2F ("h2SftIrocBgALamBgALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocBgALamBgALamMeanMinDistAProReal);
  // f2HistSftIrocBgALamBgALamMeanMinDistPioReal = new TH2F ("h2SftIrocBgALamBgALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgALamBgALamMeanMinDistPioReal);
  // f2HistSftIrocBgALamAProMeanMinDistAProReal = new TH2F ("h2SftIrocBgALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgALamAProMeanMinDistAProReal);

  // f2HistSftOrocBgLamBgLamMeanMinDistProReal = new TH2F ("h2SftOrocBgLamBgLamMeanMinDistProReal","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgLamBgLamMeanMinDistProReal);
  // f2HistSftOrocBgLamBgLamMeanMinDistPioReal = new TH2F ("h2SftOrocBgLamBgLamMeanMinDistPioReal","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgLamBgLamMeanMinDistPioReal);
  // f2HistSftOrocBgLamProMeanMinDistProReal = new TH2F ("h2SftOrocBgLamProMeanMinDistProReal","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocBgLamProMeanMinDistProReal);
  // f2HistSftOrocBgALamBgALamMeanMinDistAProReal = new TH2F ("h2SftOrocBgALamBgALamMeanMinDistAProReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocBgALamBgALamMeanMinDistAProReal);
  // f2HistSftOrocBgALamBgALamMeanMinDistPioReal = new TH2F ("h2SftOrocBgALamBgALamMeanMinDistPioReal","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgALamBgALamMeanMinDistPioReal);
  // f2HistSftOrocBgALamAProMeanMinDistAProReal = new TH2F ("h2SftOrocBgALamAProMeanMinDistProReal","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgALamAProMeanMinDistAProReal);

  // Sept'12 Use a THnSparse for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  BgLamProReal = new THnSparseF("HnSp4BgLamProReal","lamProRealQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(BgLamProReal);
  BgALamAProReal = new THnSparseF("HnSp4BgALamAProReal","alamAProRealQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(BgALamAProReal);

  // Qinv: real events
  // f3HistBgLamBgLamQinvReal = new TH3F ("h3BgLamBgLamQinvReal", "Qinv BgLamBgLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 				       ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistBgLamBgLamQinvReal);
  // f3HistBgALamBgALamQinvReal = new TH3F ("h3BgALamBgALamQinvReal", "Qinv BgALamBgALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min(d) #pi [cm]"
  // 					 ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistBgALamBgALamQinvReal);
  
  // // Two-track resolution: mixed events
  // f2HistBgLamBgLamMeanMinDistProMixed = new TH2F ("h2BgLamBgLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgLamBgLamMeanMinDistProMixed);
  // f2HistBgLamBgLamMeanMinDistPioMixed = new TH2F ("h2BgLamBgLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgLamBgLamMeanMinDistPioMixed);
  // f2HistBgLamProMeanMinDistProMixed = new TH2F ("h2BgLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistBgLamProMeanMinDistProMixed);
  // f2HistBgALamBgALamMeanMinDistAProMixed = new TH2F ("h2BgALamBgALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistBgALamBgALamMeanMinDistAProMixed);
  // f2HistBgALamBgALamMeanMinDistPioMixed = new TH2F ("h2BgALamBgALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgALamBgALamMeanMinDistPioMixed);
  // f2HistBgALamAProMeanMinDistAProMixed = new TH2F ("h2BgALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistBgALamAProMeanMinDistAProMixed);

  // f2HistSftBgLamBgLamMeanMinDistProMixed = new TH2F ("h2SftBgLamBgLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgLamBgLamMeanMinDistProMixed);
  // f2HistSftBgLamBgLamMeanMinDistPioMixed = new TH2F ("h2SftBgLamBgLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgLamBgLamMeanMinDistPioMixed);
  // f2HistSftBgLamProMeanMinDistProMixed = new TH2F ("h2SftBgLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftBgLamProMeanMinDistProMixed);
  // f2HistSftBgALamBgALamMeanMinDistAProMixed = new TH2F ("h2SftBgALamBgALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftBgALamBgALamMeanMinDistAProMixed);
  // f2HistSftBgALamBgALamMeanMinDistPioMixed = new TH2F ("h2SftBgALamBgALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgALamBgALamMeanMinDistPioMixed);
  // f2HistSftBgALamAProMeanMinDistAProMixed = new TH2F ("h2SftBgALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftBgALamAProMeanMinDistAProMixed);

  // f2HistSftIrocBgLamBgLamMeanMinDistProMixed = new TH2F ("h2SftIrocBgLamBgLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgLamBgLamMeanMinDistProMixed);
  // f2HistSftIrocBgLamBgLamMeanMinDistPioMixed = new TH2F ("h2SftIrocBgLamBgLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgLamBgLamMeanMinDistPioMixed);
  // f2HistSftIrocBgLamProMeanMinDistProMixed = new TH2F ("h2SftIrocBgLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocBgLamProMeanMinDistProMixed);
  // f2HistSftIrocBgALamBgALamMeanMinDistAProMixed = new TH2F ("h2SftIrocBgALamBgALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftIrocBgALamBgALamMeanMinDistAProMixed);
  // f2HistSftIrocBgALamBgALamMeanMinDistPioMixed = new TH2F ("h2SftIrocBgALamBgALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgALamBgALamMeanMinDistPioMixed);
  // f2HistSftIrocBgALamAProMeanMinDistAProMixed = new TH2F ("h2SftIrocBgALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftIrocBgALamAProMeanMinDistAProMixed);

  // f2HistSftOrocBgLamBgLamMeanMinDistProMixed = new TH2F ("h2SftOrocBgLamBgLamMeanMinDistProMixed","#Lambda#Lambda Mean vs min dist of decay protons;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgLamBgLamMeanMinDistProMixed);
  // f2HistSftOrocBgLamBgLamMeanMinDistPioMixed = new TH2F ("h2SftOrocBgLamBgLamMeanMinDistPioMixed","#Lambda#Lambda Mean vs min dist of decay pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgLamBgLamMeanMinDistPioMixed);
  // f2HistSftOrocBgLamProMeanMinDistProMixed = new TH2F ("h2SftOrocBgLamProMeanMinDistProMixed","p#Lambda Mean vs min dist of pri. p - dec p;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocBgLamProMeanMinDistProMixed);
  // f2HistSftOrocBgALamBgALamMeanMinDistAProMixed = new TH2F ("h2SftOrocBgALamBgALamMeanMinDistAProMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist decay #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig); 
  // fOutput2Part->Add(f2HistSftOrocBgALamBgALamMeanMinDistAProMixed);
  // f2HistSftOrocBgALamBgALamMeanMinDistPioMixed = new TH2F ("h2SftOrocBgALamBgALamMeanMinDistPioMixed","#bar{#Lambda}#bar{#Lambda} Mean vs min dist of dec. pions;mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgALamBgALamMeanMinDistPioMixed);
  // f2HistSftOrocBgALamAProMeanMinDistAProMixed = new TH2F ("h2SftOrocBgALamAProMeanMinDistProMixed","#bar{p}#bar{#Lambda} Mean vs min dist of pri. #bar{p} - dec #bar{p};mean dist. [cm];min dist [cm]",nDistBins,distLow,distHig,nDistBins,distLow,distHig);
  // fOutput2Part->Add(f2HistSftOrocBgALamAProMeanMinDistAProMixed);


  // Sept'12 Use a THnSparse for (Bg)(A)Lam(A)Pro with 4 dimensions:
  // qinv (ppri,ppri), mean dist (ppri,ppri), min dist(ppri,ppri)
  // qinv (lam pro)
  BgLamProMixed = new THnSparseF("HnSp4BgLamProMixed","lamProMixedQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(BgLamProMixed);
  BgALamAProMixed = new THnSparseF("HnSp4BgALamAProMixed","alamAProMixedQppMeanMinQlamp"
			      ,4,HnSpBins,HnSpMin,HnSpMax);
  fOutput2Part->Add(BgALamAProMixed);

  // Qinv: mixed events
  // f3HistBgLamBgLamQinvMixed = new TH3F ("h3BgLamBgLamQinvMixed", "Qinv BgLamBgLam;q_{inv} [GeV/c];min(d) p [cm];min(d) #pi [cm]"
  // 					,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistBgLamBgLamQinvMixed);
  // f3HistBgALamBgALamQinvMixed = new TH3F ("h3BgALamBgALamQinvMixed", "Qinv BgALamBgALam;q_{inv} [GeV/c];min(d) #bar{p} [cm];min(d) #pi [cm]"
  // 					  ,nQinvBins,QinvLow,QinvHig,100,0.,10.,100,0.,10.);
  // fOutput2Part->Add(f3HistBgALamBgALamQinvMixed);
  
  // Post the data
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);
  PostData(3, fOutput2Part);

}

//________________________________________________________________________
void AliAnalysisTaskProtonLambda::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  // Fill a control histogram
  fHistGoodEvent->Fill(0.0);

  // Get the event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }

  // Fill a control histogram
  fHistGoodEvent->Fill(1.0);  

  // Get the centrality selection
  AliCentrality *centrality=NULL;
  centrality = fAOD->GetCentrality();
  if (!centrality) {
    printf ("ERROR: couldn't get the AliCentrality\n");
    return;
  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(2.0);  

  // Check the fQuality flag of the centrality task
  // for details see
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies#How_we_determine_centrality
  if (centrality->GetQuality()){
    return;
  }

  // Fill a control histogram
  fHistGoodEvent->Fill(3.0);  

  // Analyze only 20% most central events using multiplicity in V0 detector (standard)
  Float_t centralityPercentile = centrality->GetCentralityPercentileUnchecked("V0M");
  if ( centralityPercentile > fkCentCut){
    return;
  }

  // Fill a control histogram
  fHistGoodEvent->Fill(4.0);  

  // Primary vertex, GetPrimaryVertex() returns the "best" reconstructed vertex
  fPrimaryVtx = fAOD->GetPrimaryVertex();
  if (!fPrimaryVtx){
    printf ("ERROR: no primary vertex\n");
    return;
  }

  // Fill a control histogram
  fHistGoodEvent->Fill(5.0);  
  fPrimaryVtx->GetXYZ(fPrimaryVtxPosition);
  // fHistPrimaryVertexPosXY->Fill(fPrimaryVtxPosition[0],fPrimaryVtxPosition[1]);
  // fHistPrimaryVertexPosZ->Fill(fPrimaryVtxPosition[2]);
  
  // Zvertex cut, probably done anyhow in centrality task
  if (TMath::Abs(fPrimaryVtxPosition[2]) > fkAbsZvertexCut)
    return;
  
  // Fill a control histogram
  fHistGoodEvent->Fill(6.0);

  // Multiplicity
  if (!(fAOD->GetNumberOfTracks())) {
    return;
  }
  
  // Fill a control histogram
  fHistGoodEvent->Fill(7.0);

  // fHistTrackMultiplicity->Fill(fAOD->GetNumberOfTracks());

  // Set up the event buffer to store this event
  fFemtoBuffer->ShiftAndAdd(fAOD);

  // // Debugging: print number of stored tracks in the event
  // for(UChar_t i=0;i<fFemtoBuffer->GetMixBuffSize();i++)
  //   printf("iMix: %u, NPro %u, NAPro %u, NLam %u, NALam %u"
  // 	   "NBgLam %u, NBgALam %u\n"
  // 	   ,i
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNPro()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNAPro()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNLam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNALam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNBgLam()
  // 	   ,fFemtoBuffer->GetEvt(i)->GetNBgALam()
  // 	   );
  // printf("\n");


  // Reset the reference array to the global tracks..
  ResetGlobalTrackReference();
  // ..and set it
  AliAODTrack *track=NULL;
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
    track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) AliFatal("Not a standard AOD");
    if (!track) continue;
    
    // Store the reference of the global tracks
    StoreGlobalTrackReference(track);
  }
  
  // V0 loop
  const Int_t nV0s = fAOD->GetNumberOfV0s();
  AliAODv0 *v0=NULL;
  AliAODTrack *pTrack=NULL;
  AliAODTrack *nTrack=NULL;
  for (Int_t iV0 = 0; iV0 < nV0s; iV0++) {
    v0 = fAOD->GetV0(iV0);

    // Skip if V0 is not there
    if((!v0))
      continue;

    // Check that the array fGTI isn't too small
    // for the track ids
    if(v0->GetPosID() >= fTrackBuffSize||
       v0->GetNegID() >= fTrackBuffSize)
      continue;

    // This is AODs: find the track for given id:
    pTrack=fGTI[v0->GetPosID()];
    nTrack=fGTI[v0->GetNegID()];
	
    // Skip if one of the daughter is not there
    if ((!pTrack) || (!nTrack)) continue;

    // Famous crossed rows / findable clusters cut,
    // rejects split tracks very well
    // (Don't do it for the V0s as we require 80 clusters 
    // and reject shared clusters)
    //    if( (!acceptTrack(pTrack)) || (!acceptTrack(nTrack)) )
    //      continue;

    // Reject tracks with shared clusters
    if(!GoodTPCFitMapSharedMap(pTrack,nTrack))
      continue;

    // Analysis done seperately for offline and on-the-fly
    if (!(v0->GetOnFlyStatus()))
      ProcessOffline(v0, pTrack, nTrack);
    else
      ProcessOnTheFly(v0, pTrack, nTrack);

    // V0s get added to the mixed events in the 'Process..' fcts
    
  } // End of V0 loop
  

  // Loop over primary tracks
  for (Int_t iTrack=0;iTrack<fAOD->GetNumberOfTracks();iTrack++){
    track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) AliFatal("Not a standard AOD");
    if (!track) continue;
    
    if(!track->TestFilterBit(128))
      continue;

    // Famous crossed rows / findable clusters cut,
    // rejects split tracks very well
    if(!acceptTrack(track))
      continue;

    // Reject tracks with shared clusters
    if(!GoodTPCFitMapSharedMap(track))
      continue;

    // Check that the array fGTI isn't too small
    // for the track id
    if(-track->GetID()-1 >= fTrackBuffSize)
      continue;

    // Without a corresponding global track it's useless
    if(!fGTI[-track->GetID()-1]){
      printf ("No global info! iTrack %d, ID %d\n",iTrack,track->GetID());
      continue;
    }

    // Visualization of TPC dE/dx
    FillDedxHist(track);

    // Depending on momentum choose pid method
    if (track->P() < 0.75){
       ProcessTPC(track);
    }
    else if (track->P() < 1.0){
       ProcessHybrid(track);
    }
    else if (track->P() < 3.25){
      ProcessTOF(track);
    }

    
    // Tracks get added to the mixed events in the 'Process..' fcts

  } // End of loop over primary tracks

  // Track cuts do not allow for split tracks

  //
  // TODO: Use Adam's shared cluster cut!
  //


  // Cleaning procedure for lambdas & lambdas, lambdas & protons,
  // anti-lambdas & anti-lambdas, anti-lambdas & protons + (anti-)lambda background
  CleaningProcedure();

  // Process real events
  ProcessReal();
  ProcessRealBackground();
  
  // Process mixed events
  ProcessMixed();
  ProcessMixedBackground();

  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fOutputPrimaries);
  PostData(3, fOutput2Part);

}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessOffline(AliAODv0 *v0, AliAODTrack *pTrack, AliAODTrack *nTrack) 
{

  // For clarity in code: Fill some hists with on-the-fly status
  //  const Float_t kOnTheFlyStat = 0.0;

  // All cuts are checked with invariant mass histograms
  //  v0->ChangeMassHypothesis(3122);
  Float_t minvLam = v0->MassLambda();
  //  v0->ChangeMassHypothesis(-3122);
  Float_t minvALam = v0->MassAntiLambda();
  // Cosine as local variable as this is some computation
  const Float_t lCosPoint = v0->CosPointingAngle(fPrimaryVtxPosition);

  // Also calculate a V0 momentum with TPC only daughters
  //  Double_t TPConlyV0Mom[3], TPConlyV0MinvLam=0, TPConlyV0MinvALam=0;
  //  getTPConlyV0Info(pTrack, nTrack,
  //		   TPConlyV0Mom, TPConlyV0MinvLam, TPConlyV0MinvALam);

  // Fill a minv hist w/o any cuts. Select background from the sideband
  fHistSideBandOffLam->Fill(minvLam);
  fHistSideBandOffALam->Fill(minvALam);
  // Fill the event buffer w/ background
  if (!fkUseOnTheFly){
    if ( TMath::Abs(minvLam - fkLamMass) > 0.015 &&
 	 TMath::Abs(minvLam - fkLamMass) < 0.035 ){
      fFemtoBuffer->GetEvt(0)->AddBgLam(v0, pTrack, nTrack);
    }
    if ( TMath::Abs(minvALam - fkLamMass) > 0.015 &&
 	 TMath::Abs(minvALam - fkLamMass) < 0.035 ){
      fFemtoBuffer->GetEvt(0)->AddBgALam(v0, pTrack, nTrack);
    }
  }

  // Control histogram: fill all v0s
  // fHistGoodV0->Fill(0.0,kOnTheFlyStat);
  // fHistGoodV0->Fill(1.0,kOnTheFlyStat);

  // Require 80 TPC clusters for both pos and neg daughter
  fHistTPCNclsPosOffLam->Fill(pTrack->GetTPCNcls(),minvLam);
  fHistTPCNclsNegOffLam->Fill(nTrack->GetTPCNcls(),minvLam);
  fHistTPCNclsPosOffALam->Fill(pTrack->GetTPCNcls(),minvALam);
  fHistTPCNclsNegOffALam->Fill(nTrack->GetTPCNcls(),minvALam);

  if ( ( (pTrack->GetTPCNcls()) < 80 ) || ( (nTrack->GetTPCNcls()) < 80 ) ) 
    return;
  //  fHistGoodV0->Fill(2.0,kOnTheFlyStat);

  // Require a maximum dca of the daughters of 0.6cm
  // fHistDcaV0DaughtersOffLam->Fill(v0->DcaV0Daughters(),minvLam);
  // fHistDcaV0DaughtersOffALam->Fill(v0->DcaV0Daughters(),minvALam);
  // fHistDcaV0Daughters->Fill(v0->DcaV0Daughters(),kOnTheFlyStat);
  if (v0->DcaV0Daughters() > 0.6)
    return;
  //  fHistGoodV0->Fill(3.0,kOnTheFlyStat);
  
  // Force TPC PID to be present
  if (!(pTrack->GetStatus() & AliVTrack::kTPCpid) ||
      !(nTrack->GetStatus() & AliVTrack::kTPCpid))
    return;
  //  fHistGoodV0->Fill(4.0,kOnTheFlyStat);

  // Visualize TPC signal before performing selection
  // fHistPosTpcBeforeCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  // fHistNegTpcBeforeCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  // The Nsigma distribution for TPC dE/dx
  // fHistPosNsigmaTpcOffLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)),minvLam);
  // fHistPosNsigmaTpcOffALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)),minvALam);
  // fHistNegNsigmaTpcOffLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)),minvLam);
  // fHistNegNsigmaTpcOffALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)),minvALam);
  // Perform cut on TPC dE/dx
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)) > 3.4)
    minvLam=0.0;
  // else 
  //   fHistPosTpcAfterCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > 4.4)
    minvLam=0.0;
  // else
  //   fHistNegTpcAfterCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > 4.2)
    minvALam=0.0;
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)) > 3.4)
    minvALam=0.0;

  // Don't use a tof cut for pions

  // Check whether to use a 5sigma tof cut or none for protons
  // if (pTrack->GetStatus() & AliVTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOffLam->Fill(1.0,minvLam);
  //   else
  //     fHistUseTofOffLam->Fill(0.0,minvLam);
  // }
  // else
  //   fHistUseTofOffLam->Fill(0.0,minvLam);
  // Check whether to use a 5sigma tof cut or none for anti-protons
  // if (nTrack->GetStatus() & AliVTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOffALam->Fill(1.0,minvALam);
  //   else
  //     fHistUseTofOffALam->Fill(0.0,minvALam);
  // }
  // else
  //   fHistUseTofOffALam->Fill(0.0,minvALam);

  // Don't use a TOF cut for offline
  
  // Don't need to check for sign of pairs as this is always
  // correct for offline finder

  // Don't need to check for TPC refit as it is required
  // by the offline finder itself

  //
  // Require a minimum distance between daughters and primary vertex
  //
  // Fill histograms with the distributions before cutting
  // fHistDcaPosOffLam->Fill(v0->DcaPosToPrimVertex(),minvLam);
  // fHistDcaPosOffALam->Fill(v0->DcaPosToPrimVertex(),minvALam);
  // fHistDcaNegOffLam->Fill(v0->DcaNegToPrimVertex(),minvLam);
  // fHistDcaNegOffALam->Fill(v0->DcaNegToPrimVertex(),minvALam);
  
  // fHistDcaPosToPrimVertex->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertex->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosToPrimVertexZoom->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertexZoom->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  
  // Do the cut
  if (v0->DcaPosToPrimVertex() < 0.1)
    minvLam=0.0;
  if (v0->DcaPosToPrimVertex() < 0.3)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.1)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.3)
    minvLam=0.0;

  // Cosine of pointing angle. Computed at the beginning.
  // Fill historgrams before cutting
  // fHistCosPointLamOff->Fill(lCosPoint,minvLam);
  // fHistCosPointALamOff->Fill(lCosPoint,minvALam);
  // fHistCosPointLamZoomOff->Fill(lCosPoint,minvLam);
  // fHistCosPointALamZoomOff->Fill(lCosPoint,minvALam);
  
  // fHistCosPointAngle->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointAngleZoom->Fill(lCosPoint,kOnTheFlyStat);
  
  // Do the cut in cos (pointing angle) 
  // (note the difference 0.9996 for offline and 0.9999 for on-the-fly)
  if (lCosPoint < 0.9996)
    return;
  
  // fHistGoodV0->Fill(7.0,kOnTheFlyStat);
  
  // Fill some histograms with cut variables
  // fHistChi2->Fill(v0->Chi2V0(),kOnTheFlyStat);
  
  // Idea to cut on the radius
  // fHistRadiusV0->Fill(v0->RadiusV0(),kOnTheFlyStat);
  // fHistV0RadiusLamOff->Fill(v0->RadiusV0(),minvLam);
  // fHistV0RadiusALamOff->Fill(v0->RadiusV0(),minvALam);

  // Idea to cut on the decay length
  // fHistDecayLengthV0->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),kOnTheFlyStat);
  // fHistV0DecayLengthLamOff->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvLam);
  // fHistV0DecayLengthALamOff->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvALam);
 
  // Idea to cut on DCA of V0 and primay vertex
  // fHistDcaV0PriVertexLamOff->Fill(v0->DcaV0ToPrimVertex(),minvLam);
  // fHistDcaV0PriVertexALamOff->Fill(v0->DcaV0ToPrimVertex(),minvALam);
     
  // Fill some invariant mass distributions
  fHistMassLambdaOff->Fill(minvLam);
  fHistMassAntiLambdaOff->Fill(minvALam);
  // fHistPtVsMassLambdaOff->Fill(v0->Pt(),minvLam);
  // fHistPtVsMassAntiLambdaOff->Fill(v0->Pt(),minvALam);

  // 3d histogram: rapidity, pt and mass
  fHistYPtMassLamOff->Fill(v0->Y(3122),v0->Pt(),minvLam);
  fHistYPtMassALamOff->Fill(v0->Y(-3122),v0->Pt(),minvALam);
  
  // Invariant mass cut lambda :: fill a y-pt hist
  // if ( TMath::Abs(minvLam - fkLamMass) < 0.01 ){
  //   fHistPtVsYLambdaOff->Fill(v0->Pt(),v0->Y(3122));
  // }
  // // Invariant mass cut anti-lambda :: fill a y-pt hist
  // if ( TMath::Abs(minvALam - fkLamMass) < 0.01 ){
  //   fHistPtVsYAntiLambdaOff->Fill(v0->Pt(),v0->Y(-3122));
  // }

  // Fill the mixed events when offline V0 finder is used
  if (!fkUseOnTheFly){
    // Highest significance for minv +/- 4 MeV
    if ( TMath::Abs(minvLam - fkLamMass) < 0.004 ){
      fFemtoBuffer->GetEvt(0)->AddLam(v0, pTrack, nTrack);
    }
    if ( TMath::Abs(minvALam - fkLamMass) < 0.004 ){
      fFemtoBuffer->GetEvt(0)->AddALam(v0, pTrack, nTrack);
    }
  }
}   
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessOnTheFly(AliAODv0 *v0, AliAODTrack *pTrack, AliAODTrack *nTrack) 
{
  // For clarity in code: Fill some hists with on-the-fly status
  //  const Float_t kOnTheFlyStat = 1.0;

  // All cuts are checked with invariant mass histograms
  Float_t minvLam = v0->MassLambda();
  Float_t minvALam = v0->MassAntiLambda();
  const Float_t lCosPoint = v0->CosPointingAngle(fPrimaryVtxPosition);

  // Control histogram: fill all v0s
  //  fHistGoodV0->Fill(0.0,kOnTheFlyStat);
  // Control hist: after require two daughter tracks
  //  fHistGoodV0->Fill(1.0,kOnTheFlyStat);
  // Check the right sign of the tracks (mainly on-the-fly)
  if (pTrack->Charge() > 0 && nTrack->Charge() < 0){
    // Correct assignment
    // fHistCorrectSigns->Fill(0.0,kOnTheFlyStat);

    // fHistLikeSignOnLam->Fill(0.0,minvLam);
    // fHistLikeSignOnALam->Fill(0.0,minvALam);    
  }
  else if (pTrack->Charge() < 0 && nTrack->Charge() > 0){
    // Swapped sings
    //    fHistCorrectSigns->Fill(1.0,kOnTheFlyStat);

    pTrack = fGTI[v0->GetNegID()];
    nTrack = fGTI[v0->GetPosID()];
    

    // See http://savannah.cern.ch/bugs/?90749
    // For AODs it depends on with which root version 
    // the AODs got produced.

    // See above: swapping mass assignment
    minvLam = v0->MassAntiLambda();
    minvALam = v0->MassLambda();

    //    fHistLikeSignOnLam->Fill(1.0,minvLam);
    //    fHistLikeSignOnALam->Fill(1.0,minvALam);    
  }
  else {
    // Like sign pairs
    //    fHistCorrectSigns->Fill(2.0,kOnTheFlyStat);
    
    //    fHistLikeSignOnLam->Fill(2.0,minvLam);
    //    fHistLikeSignOnALam->Fill(2.0,minvALam);    

    // Don't use like sign-pairs
    return;
  }
  //  fHistGoodV0->Fill(2.0,kOnTheFlyStat);

  // V0 momentum
  Double_t V0Mom[3];
  v0->PxPyPz(V0Mom);
  // Also calculate a V0 momentum with TPC only daughters
  //  Double_t TPConlyV0Mom[3], TPConlyV0MinvLam=0, TPConlyV0MinvALam=0;
  //  getTPConlyV0Info(pTrack, nTrack,
  //		   TPConlyV0Mom, TPConlyV0MinvLam, TPConlyV0MinvALam);

   // Fill a minv hist w/o any cuts. Select background from the sideband
  fHistSideBandOnLam->Fill(minvLam);
  fHistSideBandOnALam->Fill(minvALam);
  // Fill the event buffer w/ background
  if (fkUseOnTheFly){
    // Select side band aka background lambdas
    if (TMath::Abs(minvLam - fkLamMass) > 0.015 &&
	TMath::Abs(minvLam - fkLamMass) < 0.035 ){
      
      fFemtoBuffer->GetEvt(0)->AddBgLam(v0, pTrack, nTrack);
      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffBgLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //			    V0Mom[1] - TPConlyV0Mom[1],
      //			    V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
      // No SPD hits
      //	fHistMomDiffWoSPDBgLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //				   V0Mom[1] - TPConlyV0Mom[1],
      //				   V0Mom[2] - TPConlyV0Mom[2]);

      //    }
    } // End of background lambdas
    // Select side band aka background anti-lambdas
    if ( TMath::Abs(minvALam - fkLamMass) > 0.015 &&
  	 TMath::Abs(minvALam - fkLamMass) < 0.035 ){

      fFemtoBuffer->GetEvt(0)->AddBgALam(v0, pTrack, nTrack);
      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffBgALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //			       V0Mom[1] - TPConlyV0Mom[1],
      //			       V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
      // No SPD hits
      //	fHistMomDiffWoSPDBgALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //			      V0Mom[1] - TPConlyV0Mom[1],
      //				      V0Mom[2] - TPConlyV0Mom[2]);
      //      } // No SPD hits
    } // End of background anti-lambda
  } // End of if use on-the-fly finder

  //
  // Require 80 TPC clusters for both daughters
  //
  // There's a lambda signal for 0-9 clusters of the proton 
  // as it's for 110-120?!
  // There was a bug in the finding of the global track, since 
  // fixing it, offline is fine (and the problem looks less 
  // severe for on-the-fly). Still there is a problem here. 
  // There are tracks with 0 clusters. This is not the case
  // for the offline finder. The speculation would be that 
  // 1-9 clusters are treated correctly also here, it's just
  // the 0 cluster tracks. Should be a filter issue: on-the-fly
  // finds a V0, stores the daughter but info doesn't get written.
  if(pTrack->GetTPCNcls()){
    // More than zero clusters
    fHistTPCNclsPosOnLam->Fill(pTrack->GetTPCNcls(),minvLam);
    fHistTPCNclsPosOnALam->Fill(pTrack->GetTPCNcls(),minvALam);
  }
  else {
    // Zero clusters, fill the underflow to distinguish
    fHistTPCNclsPosOnLam->Fill(-1,minvLam);
    fHistTPCNclsPosOnALam->Fill(-1,minvALam);
  }
  if(nTrack->GetTPCNcls()){
    // More than zero clusters
    fHistTPCNclsNegOnLam->Fill(nTrack->GetTPCNcls(),minvLam);
    fHistTPCNclsNegOnALam->Fill(nTrack->GetTPCNcls(),minvALam);
  }
  else {
    // Zero clusters, fill the underflow to distinguish
    fHistTPCNclsNegOnLam->Fill(-1,minvLam);
    fHistTPCNclsNegOnALam->Fill(-1,minvALam);
  }
  
  // Do the cut on the TPC clusters, 0 OR at least 80
  if ( ( pTrack->GetTPCNcls() < 80 && pTrack->GetTPCNcls() ) ||
       ( nTrack->GetTPCNcls() < 80 && nTrack->GetTPCNcls() ) ) 
    return;
  //  fHistGoodV0->Fill(3.0,kOnTheFlyStat);

  // Require a maximum dca of the daughters of 0.2cm
  // fHistDcaV0DaughtersOnLam->Fill(v0->DcaV0Daughters(),minvLam);
  // fHistDcaV0DaughtersOnALam->Fill(v0->DcaV0Daughters(),minvALam);
  // fHistDcaV0Daughters->Fill(v0->DcaV0Daughters(),kOnTheFlyStat);
  if (v0->DcaV0Daughters() > 0.2)
    return;
  //  fHistGoodV0->Fill(4.0,kOnTheFlyStat);
  
  // Require cosine of pointing angle bigger than 0.9999
  // fHistCosPointAngle->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointAngleZoom->Fill(lCosPoint,kOnTheFlyStat);
  // fHistCosPointLamOn->Fill(lCosPoint,minvLam);
  // fHistCosPointALamOn->Fill(lCosPoint,minvALam);
  // fHistCosPointLamZoomOn->Fill(lCosPoint,minvLam);
  // fHistCosPointALamZoomOn->Fill(lCosPoint,minvALam);
  if (lCosPoint<0.9999)
    return;
  //  fHistGoodV0->Fill(5.0,kOnTheFlyStat);
  // Force TPC PID to be present
  if (!(pTrack->GetStatus() & AliVTrack::kTPCpid) ||
      !(nTrack->GetStatus() & AliVTrack::kTPCpid)) {
    // No TPC pid present for this track
    return;
  }
  //  fHistGoodV0->Fill(6.0,kOnTheFlyStat);
  // Visualize TPC signal before performing selection
  // fHistPosTpcBeforeCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  // fHistNegTpcBeforeCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  // // The Nsigma distribution for TPC dE/dx
  // fHistPosNsigmaTpcOnLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)),minvLam);
  // fHistPosNsigmaTpcOnALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)),minvALam);
  // fHistNegNsigmaTpcOnLam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)),minvLam);
  // fHistNegNsigmaTpcOnALam->Fill(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)),minvALam);

  // Perform cut on TPC dE/dx
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)) > 3.7)
    minvLam=0.0;
  // else 
  //   fHistPosTpcAfterCut->Fill(pTrack->P(),pTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > 3.8)
    minvLam=0.0;
  // else
  //   fHistNegTpcAfterCut->Fill(nTrack->P(),nTrack->GetTPCsignal());
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > 4.2)
    minvALam=0.0;
  if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)) > 3.9)
    minvALam=0.0;

  // Don't use a tof cut for pions

  // Check whether to use a 5sigma tof cut or none for protons
  // if (pTrack->GetStatus() & AliVTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOnLam->Fill(1.0,minvLam);
  //   else
  //     fHistUseTofOnLam->Fill(0.0,minvLam);
  // }
  // else
  //   fHistUseTofOnLam->Fill(0.0,minvLam);
  // // Check whether to use a 5sigma tof cut or none for anti-protons
  // if (nTrack->GetStatus() & AliVTrack::kTOFpid){
  //   if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
  //     fHistUseTofOnALam->Fill(1.0,minvALam);
  //   else
  //     fHistUseTofOnALam->Fill(0.0,minvALam);
  // }
  // else
  //   fHistUseTofOnALam->Fill(0.0,minvALam);

  // Reject (anti-)protons with more than 5sigma TOF
  if (nTrack->GetStatus() & AliVTrack::kTOFpid){
    if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton)) > 5.0)
      minvALam=0.0;
  }
  if (pTrack->GetStatus() & AliVTrack::kTOFpid){
    if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton)) > 5.0)
      minvLam=0.0;
  }
  
  // Don't require TPC refit. You would kill nearly your whole signal  

  // Distance between daughters and primary vertex
  // fHistDcaPosToPrimVertex->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertex->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosToPrimVertexZoom->Fill(v0->DcaPosToPrimVertex(),kOnTheFlyStat);
  // fHistDcaNegToPrimVertexZoom->Fill(v0->DcaNegToPrimVertex(),kOnTheFlyStat);
  // fHistDcaPosOnLam->Fill(v0->DcaPosToPrimVertex(),minvLam);
  // fHistDcaPosOnALam->Fill(v0->DcaPosToPrimVertex(),minvALam);
  // fHistDcaNegOnLam->Fill(v0->DcaNegToPrimVertex(),minvLam);
  // fHistDcaNegOnALam->Fill(v0->DcaNegToPrimVertex(),minvALam);
  // Require at least 0.02 cm distance from the primary vertex for the (anti-)protons
  if (v0->DcaPosToPrimVertex() < 0.02)
    minvLam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.02)
    minvALam=0.0;
  // Require at least 0.05 cm distance from the primary vertex for the pions
  if (v0->DcaPosToPrimVertex() < 0.05)
    minvALam=0.0;
  if (v0->DcaNegToPrimVertex() < 0.05)
    minvLam=0.0;
  
  // Fill some histograms with cut variables
  //  fHistChi2->Fill(v0->Chi2V0(),kOnTheFlyStat);
  
  
  // Idea to cut on the radius
  // fHistRadiusV0->Fill(v0->RadiusV0(),kOnTheFlyStat);
  // fHistV0RadiusLamOn->Fill(v0->RadiusV0(),minvLam);
  // fHistV0RadiusALamOn->Fill(v0->RadiusV0(),minvALam);
  
  // Idea to cut on the decay length
  // fHistDecayLengthV0->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),kOnTheFlyStat);
  // fHistV0DecayLengthLamOn->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvLam);
  // fHistV0DecayLengthALamOn->Fill(v0->DecayLengthV0(fPrimaryVtxPosition),minvALam);

  // Idea to cut on DCA of V0 and primay vertex
  // fHistDcaV0PriVertexLamOn->Fill(v0->DcaV0ToPrimVertex(),minvLam);
  // fHistDcaV0PriVertexALamOn->Fill(v0->DcaV0ToPrimVertex(),minvALam);

  // TPC Chi2 / number of degrees of freedom
  // A cut on at least 80 clusters is already done before,
  // no concern to divide by zero
  // fHistChi2TPCPosLamOn->Fill(pTrack->Chi2perNDF(),minvLam);
  // fHistChi2TPCPosALamOn->Fill(pTrack->Chi2perNDF(),minvALam);
  // fHistChi2TPCNegLamOn->Fill(nTrack->Chi2perNDF(),minvLam);
  // fHistChi2TPCNegALamOn->Fill(nTrack->Chi2perNDF(),minvALam);
  // Don't cut like Chi2/ndf < 4! One might throw away the tracks
  // with Chi2/ndf roughly one as they are good primaries

  // Fill some invariant mass distributions
  fHistMassLambdaOn->Fill(minvLam);
  fHistMassAntiLambdaOn->Fill(minvALam);
  // fHistPtVsMassLambdaOn->Fill(v0->Pt(),minvLam);
  // fHistPtVsMassAntiLambdaOn->Fill(v0->Pt(),minvALam);

  // TPC only invariant mass distributions
  //  if(minvLam > .1){
    // Lambda is good
    //    fHistMinvTPConlyLamOn->Fill(TPConlyV0MinvLam);
  //  }
  //  if (minvALam > .1){
    // Anti-lambda is good
    //    fHistMinvTPConlyALamOn->Fill(TPConlyV0MinvALam);
  //  }
  
  // 3d histogram: rapidity, pt and mass
  fHistYPtMassLamOn->Fill(v0->Y(3122),v0->Pt(),minvLam);
  fHistYPtMassALamOn->Fill(v0->Y(-3122),v0->Pt(),minvALam);
  
  // // Invariant mass cut lambda :: fill a y-pt hists
  // if ( TMath::Abs(minvLam - fkLamMass) < 0.01 ){
  //   fHistPtVsYLambdaOn->Fill(v0->Pt(),v0->Y(3122));
  // }
  // // Invariant mass cut anti-lambda :: fill a y-pt hists
  // if ( TMath::Abs(minvALam - fkLamMass) < 0.01 ){
  //   fHistPtVsYAntiLambdaOn->Fill(v0->Pt(),v0->Y(-3122));
  // }
  
  // Fill the mixed events when on-the-fly V0 finder is used
  if (fkUseOnTheFly){

    // Highest significance for minv +/- 4 MeV
    if ( TMath::Abs(minvLam - fkLamMass) < 0.004 ){
      fFemtoBuffer->GetEvt(0)->AddLam(v0, pTrack, nTrack);
      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //		    V0Mom[1] - TPConlyV0Mom[1],
      //		    V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
      //	// No SPD hits
      //	fHistMomDiffWoSPDLam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //				   V0Mom[1] - TPConlyV0Mom[1],
      //				   V0Mom[2] - TPConlyV0Mom[2]);
      // } // No SPD hits
    } // Good lambda
    if ( TMath::Abs(minvALam - fkLamMass) < 0.004 ) {
      fFemtoBuffer->GetEvt(0)->AddALam(v0, pTrack, nTrack);
      // Momentum difference of standard V0 / TPC only V0
      //      fHistMomDiffALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //			    V0Mom[1] - TPConlyV0Mom[1],
      //			    V0Mom[2] - TPConlyV0Mom[2]);
      // Same excluding V0s with daughters with SPD hits
      //      if( !(pTrack->HasPointOnITSLayer(0) || pTrack->HasPointOnITSLayer(1) ||
      //	    nTrack->HasPointOnITSLayer(0) || nTrack->HasPointOnITSLayer(1) )){
	// No SPD hits
      //	fHistMomDiffWoSPDALam->Fill(V0Mom[0] - TPConlyV0Mom[0], 
      //				   V0Mom[1] - TPConlyV0Mom[1],
      //				   V0Mom[2] - TPConlyV0Mom[2]);
      //    } // No SPD hits
    } // Good anti-lambda
  } // Use on-the-fly finder for Femto analysis
} // ProcessOnTheFly
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessTOF(AliAODTrack* track) 
{
  // Request the kTOFpid bit. There are tracks with kTOFout and wihthout kTOFpid,
  // but these tracks have a bad TOF signal.
  if(!((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid))
    return;

  // TOF signal corrected for expected time and (if neccessary) for start time
  Float_t corrTOFsig = GetCorrectedTOFSignal(track);

  // Distinguish between charges
  if (track->Charge() > 0){
    // Simple Nsigma TOF distribution
    //    fPriHistPosNsigmaTof->Fill(fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of total momentum
    fPriHistPosNsigmaTofVsP->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of transverse momentum
    fPriHistPosNsigmaTofVsPt->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    
    // Try the tof signal instead of nsigma
    fPriHistTOFsignalPosVsP->Fill(track->P(), corrTOFsig);
    fPriHistTOFsignalPosVsPt->Fill(track->Pt(), corrTOFsig);
    
  }
  else if (track->Charge() < 0){
    // Simple Nsigma TOF distribution
    //    fPriHistNegNsigmaTof->Fill(fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of total momentum
    fPriHistNegNsigmaTofVsP->Fill(track->P(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    // Nsigma TOF in bins of transverse momentum
    fPriHistNegNsigmaTofVsPt->Fill(track->Pt(),fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton));
    
    // Try the tof signal instead of nsigma
    fPriHistTOFsignalNegVsP->Fill(track->P(), corrTOFsig);
    fPriHistTOFsignalNegVsPt->Fill(track->Pt(), corrTOFsig);
  }
  
  // Final judging: simple first idea. 
  // min -800 up to 2 GeV and 0 up to 3.25GeV
  if (track->P() < 2.0){
    //    if (corrTOFsig > -800.0){
    // In AODs, the resolution is better, do -500 (AODs) 
    // instead of -800 (ESDs)
    if (corrTOFsig > -500.0) {
      // Create additional TPC only constrained tp pri. vtx track parameters
      //	constrainTrack(track);
      if (track->Charge()>0){
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddPro(track);
	}
      }
      else{
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddAPro(track);
	}
      }
    }
  }
  else if (track->P() < 3.25){
    if (corrTOFsig > 0){
      // Create additional TPC only constrained tp pri. vtx track parameters
      //	constrainTrack(track);
      if (track->Charge()>0){
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddPro(track);
	}
      }
      else{
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddAPro(track);
	}
      }
    }
  }
} // End of void ProcessTOF
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessTPC(AliAODTrack* track){

  // Require the TPCpid bit
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTPCpid))
    return;
    
  // In contrast to ESDs one doesn't check for AliESDtrack::kTOFpid
  // but for AliVTrack::kTOFout?? 
  // Check how many particles have TOFout bit
  // if (track->Charge() > 0){
  //   if ((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid)
  //     fPriHistHasTofPos->Fill(1.0);
  //   else 
  //     fPriHistHasTofPos->Fill(0.0);
  // }
  // else{
  //   if ((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid)
  //     fPriHistHasTofNeg->Fill(1.0);
  //   else 
  //     fPriHistHasTofNeg->Fill(0.0);
  // }

  // For all plots <dE/dx> vs p one should use
  // the momentum at the inner wall of the TPC.

  // Use a TOF cut and fill the same dE/dx histograms
  // Bool_t acceptedTOF=kFALSE;
  // if ((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid){
  //   if (fPIDResponse->NumberOfSigmasTOF((fGTI[-track->GetID()-1]), AliPID::kProton) > -10.0)
  //     acceptedTOF=kTRUE;
  // }
  // if (acceptedTOF){
  //     if (track->Charge() > 0){
  // 	fPriHistTPCsignalTOFcutPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 (fGTI[-track->GetID()-1])->GetTPCsignal());
  // 	fPriHistNsigmaTPCTOFcutPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 fPIDResponse->NumberOfSigmasTPC((fGTI[-track->GetID()-1]), AliPID::kProton));
  //     }
  //     else{
  // 	fPriHistTPCsignalTOFcutNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 (fGTI[-track->GetID()-1])->GetTPCsignal());
  // 	fPriHistNsigmaTPCTOFcutNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
  // 					 fPIDResponse->NumberOfSigmasTPC((fGTI[-track->GetID()-1]), AliPID::kProton));
  //     }
  // }	
    
  // A first idea of a cut: use the spectra cut.
  // (should perhaps change for momenta ~ 0.75 GeV)
  if ( ((fGTI[-track->GetID()-1])->GetTPCsignal() > 
	fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
					AliPID::kProton))
      // New since Sept 10th 2012: Also use a cut to reject deuterons.
      // I checked: The cut is good!
       && ((fGTI[-track->GetID()-1])->GetTPCsignal() <
	   2.0*fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
				      AliPID::kProton))
                        ) {
    // Distinguish between charges
    if (track->Charge()>0){
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the femto event
	fFemtoBuffer->GetEvt(0)->AddPro(track);
      }
    }
    else{
      // Cut .1 cm on DCAxy and fill a histogram
      if(goodDCA(track)){
	// Add to the femto event
	fFemtoBuffer->GetEvt(0)->AddAPro(track);
      }
    }
  }
} // End of void ProcessTPC
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessHybrid(AliAODTrack *track){
  
  // Intermediate momentum: use dEdx for a pre-selection
  // and do the pid with tof
  
  // Boolean for extra! tpc pid cuts
  Bool_t acceptTPC = kTRUE;

  // Require the TPCpid bit
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTPCpid))
    acceptTPC = kFALSE;
 
  // Pre-selection cut with TPC, don't return immediately to be able
  // to visualize the effect
  if (acceptTPC){
    // Do a mild dEdx cut
    if ((fGTI[-track->GetID()-1])->GetTPCsignal() < 
	fTpcResponse->GetExpectedSignal((fGTI[-track->GetID()-1])->GetTPCmomentum(),
					AliPID::kElectron))
      acceptTPC = kFALSE;
  }
    
  // Ask for TOF pid flag and fill
  if (!((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid))
    return;
  
  // The corrected TOF signal
  Double_t corrTOFsig = GetCorrectedTOFSignal(track);
  
  // Distinguish between charges
  if (track->Charge() > 0) {
    // Fill the tof signal w/o dedx pre-selection
    fPriHistHybridTOFsigPosWoTPC->Fill(corrTOFsig);
    // Do the pre-selection
    if (acceptTPC){
      fPriHistHybridTOFsigPosTPCok->Fill(corrTOFsig);

      // Do the tof cut
      // Sept '12: also include an upper cut
      if ( (corrTOFsig > -1000.0) && (corrTOFsig < 1250.) ){
	// Create additional TPC only constrained to pri. vtx track parameters
	//	constrainTrack(track);
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // Add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddPro(track);
	}
      }
    }
  }
  else {
    // Fill the tof signal w/o dedx pre-selection
    fPriHistHybridTOFsigNegWoTPC->Fill(corrTOFsig);
    // Do the pre-selection
    if (acceptTPC){
      fPriHistHybridTOFsigNegTPCok->Fill(corrTOFsig);
      
      // Do the tof cut
      // Sept '12: also include an upper cut
      if ( (corrTOFsig > -1000.0) && (corrTOFsig < 1250.) ){
	// Create additional TPC only constrained to pri. vtx track parameters
	//	constrainTrack(track);
	// Cut .1 cm on DCAxy and fill a histogram
	if(goodDCA(track)){
	  // add to the femto event
	  fFemtoBuffer->GetEvt(0)->AddAPro(track);
	}
      }
    }
  }
} // End of ProcessHybrid
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::CleaningProcedure() {
  // fFemtoBuffer->GetEvt(0) pointer must be set
  // Checks that no tracks are shared between Lam & Lam, Lam & Pro, ALam & ALam, ALam & APro

  // printf ("Cleaning procedure. Lam: %d, ALam: %d, Pro: %d, APro:%d\n"
  // 	  ,fFemtoBuffer->GetEvt(0)->GetNLam(),fFemtoBuffer->GetEvt(0)->GetNALam(),fFemtoBuffer->GetEvt(0)->GetNPro(),fFemtoBuffer->GetEvt(0)->GetNAPro());

  //
  // Check for lambdas..
  //
  for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNLam();i++) {
    if (!fFemtoBuffer->GetEvt(0)->fLamTracks[i].UseIt())
      continue;
    // Unique track id's for first V0
    Int_t posId1 = fFemtoBuffer->GetEvt(0)->fLamTracks[i].fPosDaughter.fID;
    Int_t negId1 = fFemtoBuffer->GetEvt(0)->fLamTracks[i].fNegDaughter.fID;

    // .. & lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNLam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fLamTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fLamTracks[j].fPosDaughter.fID;
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fLamTracks[j].fPosDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){

	// printf ("shared track lamlam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
  	// 	posId1, posId2, negId1, negId2);
	
  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fLamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fLamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fLamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fLamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fLamTracks[i].UseIt())
      continue;

    // .. & protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNPro();j++){
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fProTracks[j].fID;
      
      // If V0 and proton share a track
      if (posId1 == posId2){
	//  	printf ("shared track lam p! id:%d\n",posId1);
  	
	// Remove the proton
  	fFemtoBuffer->GetEvt(0)->fProTracks[j].SetBadFlag();
      }
      
    } // Proton loop

  } // First V0 loop

  //
  // Check for anti-lambdas..
  //
  for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNALam();i++){
    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[i].UseIt())
      continue;
    // Unique track id's for first V0
    Int_t posId1 = fFemtoBuffer->GetEvt(0)->fALamTracks[i].fPosDaughter.fID;
    Int_t negId1 = fFemtoBuffer->GetEvt(0)->fALamTracks[i].fNegDaughter.fID;

    // .. & anti-lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNALam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fALamTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fALamTracks[j].fPosDaughter.fID;
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fALamTracks[j].fNegDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){
	
	// printf ("shared track ALamALam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
  	// 	posId1, posId2, negId1, negId2);

  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fALamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fALamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fALamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fALamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd anti-V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[i].UseIt())
      continue;
    
    // .. & anti-protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNAPro();j++){
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fAProTracks[j].fID;
      
      // If V0 and proton share a track
      if (negId1 == negId2){
	//  	printf ("shared track alam ap! id:%d\n",posId1);

  	// Remove the proton
  	fFemtoBuffer->GetEvt(0)->fAProTracks[j].SetBadFlag();
      }
      
    } // Anti-proton loop

  } // First anti-V0 loop

  //
  // Do the same with the side band background.
  // Discard background when sharing track with primary proton.
  //

   for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNBgLam();i++){
    if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
      continue;
    // Unique track id's for first V0
    Int_t posId1 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fPosDaughter.fID;
    Int_t negId1 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fNegDaughter.fID;

    // .. & lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNBgLam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fPosDaughter.fID;
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fNegDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){

	// printf ("shared track bglambglam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
	// 	posId1, posId2, negId1, negId2);
	
  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fBgLamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
      continue;

    // .. & protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNPro();j++) {
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[j].UseIt())
  	continue;
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].UseIt())
  	continue;

      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fProTracks[j].fID;
      
      // If V0 and proton share a track
      if (posId1 == posId2){
	//  	printf ("shared track bglam p! id:%d\n",posId1);
  	// Remove the background lambda
  	fFemtoBuffer->GetEvt(0)->fBgLamTracks[i].SetBadFlag();
      }
      
    } // Proton loop

  } // First V0 loop

  //
  // Check for anti-lambdas..
  //
  for (Int_t i=0;i<fFemtoBuffer->GetEvt(0)->GetNBgALam();i++){
    if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].UseIt())
      continue;
    // Unique track id's for first V0
    Int_t posId1 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fPosDaughter.fID;
    Int_t negId1 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fNegDaughter.fID;

    // .. & anti-lambdas
    for (Int_t j=i+1;j<fFemtoBuffer->GetEvt(0)->GetNBgALam();j++){
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].UseIt())
  	continue;
      // Unique track id's for second V0
      Int_t posId2 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fPosDaughter.fID;
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fNegDaughter.fID;
      
      // If V0s share a track remove one
      if (posId1 == posId2 || negId1 == negId2){
	
	// printf ("shared track BgALamBgALam! posId1: %d, posId2: %d, negId1: %d, negId2: %d\n",
  	// 	posId1, posId2, negId1, negId2);

  	// Use a criterion to select best V0
  	if (fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].fCosPoint > fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].fCosPoint){
  	  fFemtoBuffer->GetEvt(0)->fBgALamTracks[j].SetBadFlag();
  	}
  	else{
  	  fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].SetBadFlag();
  	}
      }
      
    } // Scnd anti-V0 loop

    if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].UseIt())
      continue;
    
    // .. & anti-protons
    for (Int_t j=0;j<fFemtoBuffer->GetEvt(0)->GetNAPro();j++){
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[j].UseIt())
  	continue;
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].UseIt())
  	continue;
      
      // Unique track id's for second V0
      Int_t negId2 = fFemtoBuffer->GetEvt(0)->fAProTracks[j].fID;
      
      // If V0 and proton share a track
      if (negId1 == negId2){
	//  	printf ("shared track bgalam ap! id:%d\n",posId1);
  	// Remove the background anti-lambda
  	fFemtoBuffer->GetEvt(0)->fBgALamTracks[i].SetBadFlag();
      }
      
    } // Anti-proton loop

  } // First anti-V0 loop

  
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessReal() {
  // Process real events
  
  // // Count the number of pairs before TTR cut
  // Int_t nLamLamPairsWoTTR = 0,nLamProPairsWoTTR=0,nALamALamPairsWoTTR=0,nALamAProPairsWoTTR=0;
  // // and with TTR cut
  // Int_t nLamLamPairs = 0,nLamProPairs=0,nALamALamPairs=0,nALamAProPairs=0;

  // Declare numbers to speed up the code
  Int_t iLam,//iLam2,
    iRadius,iPro,iALam,
    //iALam2,
    iAPro,nMeasPro,
    //nMeasPio,
    nMeasAPro;
  //  Int_t nMeasProIroc,nMeasPioIroc,nMeasAProIroc,nMeasProOroc,nMeasPioOroc,nMeasAProOroc;
  // Float_t distPro,distPio,minDistPro,meanDistPro,minDistPio,meanDistPio,
  //   distAPro,minDistAPro,meanDistAPro;
  Float_t distSftPro,//distSftPio,
    minDistSftPro,meanDistSftPro,//minDistSftPio,meanDistSftPio,
    distSftAPro,minDistSftAPro,meanDistSftAPro;
  // Float_t minDistSftIrocPro,meanDistSftIrocPro,minDistSftIrocPio,meanDistSftIrocPio,
  //   minDistSftIrocAPro,meanDistSftIrocAPro;
  // Float_t minDistSftOrocPro,meanDistSftOrocPro,minDistSftOrocPio,meanDistSftOrocPio,
  //   minDistSftOrocAPro,meanDistSftOrocAPro;
  
  // printf("Real event, NLam: %d, NPro %d, NALam %d, NAPro %d\n",
  // 	 fFemtoBuffer->GetEvt(0)->GetNLam(),
  // 	 fFemtoBuffer->GetEvt(0)->GetNPro(),
  // 	 fFemtoBuffer->GetEvt(0)->GetNALam(),
  // 	 fFemtoBuffer->GetEvt(0)->GetNAPro()
  // 	 );
  // Lambda loop
  for (iLam = 0; iLam < fFemtoBuffer->GetEvt(0)->GetNLam(); iLam++){

    // Skip if unUseIt() entry
    if (!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].UseIt())
      continue;
    
    // // Second lambda loop
    // for (iLam2 = iLam+1; iLam2 < fFemtoBuffer->GetEvt(0)->GetNLam(); iLam2++){

    //   // Skip if unUseIt() entry
    //   if (!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2].UseIt())
    // 	continue;

    //   // Count the number of lam-lam pairs
    //   nLamLamPairsWoTTR++;

    //   // Reset the distances for each pair
    //   minDistPro=999.0;meanDistPro=0.0;minDistPio=999.0;meanDistPio=0.0;
    //   minDistSftPro=999.0;meanDistSftPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
    //   minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
    //   minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
    //   // Reset the number of measurements for the mean
    //   nMeasPro=0;nMeasPio=0;nMeasProIroc=0;nMeasPioIroc=0;
    //   nMeasProOroc=0;nMeasPioOroc=0;

    //   // Check for two-track resolution
    //   for (iRadius=0;iRadius<9;iRadius++){
    // 	// Get the spatial distance at each radius
    // 	distPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2].fPosDaughter.fXglobal[iRadius]);
    // 	distPio = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2].fNegDaughter.fXglobal[iRadius]);
    // 	// Shifted distances
    // 	distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));
    // 	distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

    // 	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
    // 	if (distPro > -1.0) {
    // 	  // Minimum distance
    // 	  if (distPro < minDistPro)
    // 	    minDistPro = distPro;
    // 	  if (distSftPro < minDistSftPro)
    // 	    minDistSftPro = distSftPro;
    // 	  // Mean distance
    // 	  meanDistPro+=distPro;
    // 	  meanDistSftPro+=distSftPro;
    // 	  nMeasPro++;
	
    // 	  // IROC
    // 	  if (iRadius<3){
    // 	    if (distSftPro < minDistSftIrocPro)
    // 	      minDistSftIrocPro = distSftPro;
    // 	    meanDistSftIrocPro+=distSftPro;
    // 	    nMeasProIroc++;
    // 	  }
    // 	  // OROC
    // 	  else {
    // 	    if (distSftPro < minDistSftOrocPro)
    // 	      minDistSftOrocPro = distSftPro;
    // 	    meanDistSftOrocPro+=distSftPro;
    // 	    nMeasProOroc++;	    
    // 	  }
    // 	}
    // 	if (distPio > -1.0){
    // 	  // Minimum distance
    // 	  if (distPio < minDistPio)
    // 	    minDistPio = distPio;
    // 	  if (distSftPio < minDistSftPio)
    // 	    minDistSftPio = distSftPio;
    // 	  // Mean distance
    // 	  meanDistPio+=distPio;
    // 	  meanDistSftPio+=distSftPio;
    // 	  nMeasPio++; 
	  
    // 	  // IROC
    // 	  if (iRadius<3){
    // 	    if (distSftPio < minDistSftIrocPio)
    // 	      minDistSftIrocPio = distSftPio;
    // 	    meanDistSftIrocPio+=distSftPio;
    // 	    nMeasPioIroc++;
    // 	  }
    // 	  // OROC
    // 	  else {
    // 	    if (distSftPio < minDistSftOrocPio)
    // 	      minDistSftOrocPio = distSftPio;
    // 	    meanDistSftOrocPio+=distSftPio;
    // 	    nMeasPioOroc++;	    
    // 	  }

    // 	}

    //   } // Loop over iRadius
      
    //   // Require at least one measurement
    //   if ( (!nMeasPio) || (!nMeasPro) )
    // 	continue;

    //   // Divide by the number of measurements to get the mean
    //   meanDistPro /= (Float_t)nMeasPro;
    //   meanDistPio /= (Float_t)nMeasPio;
    //   meanDistSftPro /= (Float_t)nMeasPro;
    //   meanDistSftPio /= (Float_t)nMeasPio;
      
    //   // Fill the two track resolution histograms
    //   f2HistLamLamMeanMinDistProReal->Fill(meanDistPro,minDistPro);
    //   f2HistLamLamMeanMinDistPioReal->Fill(meanDistPio,minDistPio);

    //   f2HistSftLamLamMeanMinDistProReal->Fill(meanDistSftPro,minDistSftPro);
    //   f2HistSftLamLamMeanMinDistPioReal->Fill(meanDistSftPio,minDistSftPio);

    //   // Fill IROC / OROC histograms only with at least one measurement
    //   if (nMeasProIroc){
    // 	meanDistSftIrocPro /= (Float_t)nMeasProIroc;
    // 	f2HistSftIrocLamLamMeanMinDistProReal->Fill(meanDistSftIrocPro,minDistSftIrocPro);
    //   }
    //   if (nMeasPioIroc){
    // 	meanDistSftIrocPio /= (Float_t)nMeasPioIroc;
    // 	f2HistSftIrocLamLamMeanMinDistPioReal->Fill(meanDistSftIrocPio,minDistSftIrocPio);
    //   }
    //   if (nMeasProOroc){
    // 	meanDistSftOrocPro /= (Float_t)nMeasProOroc;
    // 	f2HistSftOrocLamLamMeanMinDistProReal->Fill(meanDistSftOrocPro,minDistSftOrocPro);
    //   }
    //   if (nMeasPioOroc){
    // 	meanDistSftOrocPio /= (Float_t)nMeasPioOroc;
    // 	f2HistSftOrocLamLamMeanMinDistPioReal->Fill(meanDistSftOrocPio,minDistSftOrocPio);
    //   }

    //   // Do a cut (value needs to be refined)
    //   //      if ( meanDistSftPro < 2.0 || meanDistSftPio < 2.0 )
    //   //	continue;
      
    //   // Count the number of pairs
    //   nLamLamPairs++;

    //   // Mt of the pair
    //   f2HistMtLamLamReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
    //    				  fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2]));

    //   // Fill the qinv, minv histogram
    //   f3HistLamLamQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2]),minDistSftPro,minDistSftPio);
    //   f3HistLamLamMinvReal->Fill(Minv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2]),minDistSftPro,minDistSftPio);
      
    //   // Mt of the pair fr low q pairs only
    //   if(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
    // 	      fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2]) < .2)
    // 	f2HistMtLowQLamLamReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
    // 					fFemtoBuffer->GetEvt(0)->fLamTracks[iLam2]));

    // } // Second lambda loop

    // Proton loop
    for (iPro=0;iPro<fFemtoBuffer->GetEvt(0)->GetNPro();iPro++){

      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
  	continue;

      //      printf(".");
      //      nLamProPairsWoTTR++;

      // Reset the distances for each pair
      // minDistPro=999.0;meanDistPro=0.0;
      minDistSftPro=999.0;meanDistSftPro=0.0;
      //minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;
      // Reset the number of measurements for the mean
      nMeasPro=0;//nMeasProIroc=0;nMeasProOroc=0;
      
      // Check for two-track resolution
      for (iRadius=0;iRadius<9;iRadius++){
  	// Get the spatial distance at each radius
	//  	distPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fProTracks[iPro].fXglobal[iRadius]);
  	distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fProTracks[iPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
	if (distSftPro > -1.0){
	  // Minimum distance
	  // if (distPro < minDistPro)
	  //   minDistPro = distPro;
	  if (distSftPro < minDistSftPro)
	    minDistSftPro = distSftPro;
	  // Mean distance
	  //	  meanDistPro+=distPro;
	  meanDistSftPro+=distSftPro;
	  nMeasPro++;
	  
	  // // IROC                                   
          // if (iRadius<3){
          //   if (distSftPro < minDistSftIrocPro)
          //     minDistSftIrocPro = distSftPro;
          //   meanDistSftIrocPro+=distSftPro;
          //   nMeasProIroc++;
          // }
          // // OROC                                   
          // else {
          //   if (distSftPro < minDistSftOrocPro)
          //     minDistSftOrocPro = distSftPro;
          //   meanDistSftOrocPro+=distSftPro;
          //   nMeasProOroc++;
          // }


	}
      } // Loop over iRadius

      // Require at least one measurement
      if ( !nMeasPro )
        continue;

      // Divide by the number of measurements to get the mean
      //      meanDistPro /= (Float_t)nMeasPro;
      meanDistSftPro /= (Float_t)nMeasPro;
      
      // Fill the two track resolution histogram
      // f2HistLamProMeanMinDistProReal->Fill(meanDistPro,minDistPro);
      // f2HistSftLamProMeanMinDistProReal->Fill(meanDistSftPro,minDistSftPro);

      // // Fill IROC / OROC histograms only with at least one measurement                                
      // if (nMeasProIroc){
      //   meanDistSftIrocPro /= (Float_t)nMeasProIroc;
      //   f2HistSftIrocLamProMeanMinDistProReal->Fill(meanDistSftIrocPro,minDistSftIrocPro);
      // }
      // if (nMeasProOroc){
      //   meanDistSftOrocPro /= (Float_t)nMeasProOroc;
      //   f2HistSftOrocLamProMeanMinDistProReal->Fill(meanDistSftOrocPro,minDistSftOrocPro);
      // }
  	
      // Do a cut (value needs to be refined)
      //      if ( meanDistSftPro < 2.0 )
      //  	continue;

      // Look at possible residual correlations of the daughters
      // vs mean dist of protons
      //      f2HistLamPosDProQinvReal->Fill(QinvProPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter, fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro);
      //      f2HistLamNegDProQinvReal->Fill(QinvPioPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter, fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro);
  

      // Count the number of pairs
      //      nLamProPairs++;

      // Mt of the pair
      f2HistMtLamProReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
				  fFemtoBuffer->GetEvt(0)->fProTracks[iPro]));

      // THnSparse with qinvpropro, mean dist propro, min dist propro, qinv lampro
      Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter,
		   fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),
	meanDistSftPro,
	minDistSftPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
	     fFemtoBuffer->GetEvt(0)->fProTracks[iPro])
      };
      LamProReal->Fill(x);

      // Fill the qinv histogram, using TPC only momentum for primary protons
      //      f3HistLamProQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro,minDistSftPro); // Using THnSparse since sept '12
      //      f3HistLamProMinvReal->Fill(Minv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro,minDistSftPro);
      // and using TPC only plus primary vertex constraint momentum for primary protons
      //      f3HistLamProQinvConstrReal->Fill(QinvConstr(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro,minDistSftPro);

      // Mt of the pair for low q pairs only
      if(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
	      fFemtoBuffer->GetEvt(0)->fProTracks[iPro]) < .2)
	f2HistMtLowQLamProReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
					fFemtoBuffer->GetEvt(0)->fProTracks[iPro]));
      
    }// Proton loop

  }// Lambda loop


  
  // Anti-lambda loop
  for (iALam = 0; iALam < fFemtoBuffer->GetEvt(0)->GetNALam(); iALam++){

    // Skip if unUseIt() entry
    if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].UseIt())
      continue;

    // // Second anti-lambda loop
    // for (iALam2 = iALam+1; iALam2 < fFemtoBuffer->GetEvt(0)->GetNALam(); iALam2++){

    //   // Skip if unUseIt() entry
    //   if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2].UseIt())
    // 	continue;

    //   // Count the number of antilam-antilam pairs
    //   nALamALamPairsWoTTR++;

    //   // Reset the distances for each pair
    //   minDistAPro=999.0;meanDistAPro=0.0;minDistPio=999.0;meanDistPio=0.0;
    //   minDistSftAPro=999.0;meanDistSftAPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
    //   minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
    //   minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
    //   // Reset the number of measurements for the mean
    //   nMeasAPro=0;nMeasPio=0;nMeasAProIroc=0;nMeasPioIroc=0;nMeasAProOroc=0;nMeasPioOroc=0;

    //   // Check for two-track resolution
    //   for (iRadius=0;iRadius<9;iRadius++){
    // 	// Get the spatial distance at each radius
    // 	distPio = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2].fPosDaughter.fXglobal[iRadius]);
    // 	distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2].fNegDaughter.fXglobal[iRadius]);
    // 	// Shifted distances
    // 	distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));
    // 	distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

    // 	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
    //     if (distAPro > -1.0){
    //       // Minimum distance
    //       if (distAPro < minDistAPro)
    //         minDistAPro = distAPro;
    // 	  if (distSftAPro < minDistSftAPro)
    //         minDistSftAPro = distSftAPro;
    //       // Mean distance
    //       meanDistAPro+=distAPro;
    //       meanDistSftAPro+=distSftAPro;
    //       nMeasAPro++;

    // 	  // IROC                                                      
    //       if (iRadius<3){
    //         if (distSftAPro < minDistSftIrocAPro)
    //           minDistSftIrocAPro = distSftAPro;
    //         meanDistSftIrocAPro+=distSftAPro;
    //         nMeasAProIroc++;
    //       }
    //       // OROC                                                         
    //       else {
    //         if (distSftAPro < minDistSftOrocAPro)
    //           minDistSftOrocAPro = distSftAPro;
    //         meanDistSftOrocAPro+=distSftAPro;
    //         nMeasAProOroc++;
    //       }

    //     }
    //     if (distPio > -1.0){
    //       // Minimum distance
    //       if (distPio < minDistPio)
    //         minDistPio = distPio;
    //       if (distSftPio < minDistSftPio)
    //         minDistSftPio = distSftPio;
    //       // Mean distance
    //       meanDistPio+=distPio;
    //       meanDistSftPio+=distSftPio;
    //       nMeasPio++; 

    // 	  // IROC                                                         
    //       if (iRadius<3){
    //         if (distSftPio < minDistSftIrocPio)
    //           minDistSftIrocPio = distSftPio;
    //         meanDistSftIrocPio+=distSftPio;
    //         nMeasPioIroc++;
    //       }
    //       // OROC                                                         
    //       else {
    //         if (distSftPio < minDistSftOrocPio)
    //           minDistSftOrocPio = distSftPio;
    //         meanDistSftOrocPio+=distSftPio;
    //         nMeasPioOroc++;
    //       }
    //     }
    //   } // Loop over iRadius
      
    //   // Require at least one measurement
    //   if ( (!nMeasPio) || (!nMeasAPro) )
    //     continue;

    //   // Divide by the number of measurements to get the mean
    //   meanDistAPro /= (Float_t)nMeasAPro;
    //   meanDistPio /= (Float_t)nMeasPio;
    //   meanDistSftAPro /= (Float_t)nMeasAPro;
    //   meanDistSftPio /= (Float_t)nMeasPio;

    //   // Fill the two track resolution histograms
    //   f2HistALamALamMeanMinDistAProReal->Fill(meanDistAPro,minDistAPro);
    //   f2HistALamALamMeanMinDistPioReal->Fill(meanDistPio,minDistPio);

    //   f2HistSftALamALamMeanMinDistAProReal->Fill(meanDistSftAPro,minDistSftAPro);
    //   f2HistSftALamALamMeanMinDistPioReal->Fill(meanDistSftPio,minDistSftPio);

    //   // Fill IROC / OROC histograms only with at least one measurement   
    //   if (nMeasAProIroc){
    //     meanDistSftAPro /= (Float_t)nMeasAProIroc;
    //     f2HistSftIrocALamALamMeanMinDistAProReal->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
    //   }
    //   if (nMeasPioIroc){
    //     meanDistSftPio /= (Float_t)nMeasPioIroc;
    //     f2HistSftIrocALamALamMeanMinDistPioReal->Fill(meanDistSftIrocPio,minDistSftIrocPio);
    //   }
    //   if (nMeasAProOroc){
    //     meanDistSftAPro /= (Float_t)nMeasAProOroc;
    //     f2HistSftOrocALamALamMeanMinDistAProReal->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
    //   }
    //   if (nMeasPioOroc){
    //     meanDistSftPio /= (Float_t)nMeasPioOroc;
    //     f2HistSftOrocALamALamMeanMinDistPioReal->Fill(meanDistSftOrocPio,minDistSftOrocPio);
    //   }


    //   //      // Do a cut (value needs to be refined)
    //   //      if ( meanDistSftAPro < 2.0 || meanDistSftPio < 2.0 )
    //   //  	continue;
      
    //   // Count the number of pairs
    //   nALamALamPairs++;

    //   // Mt of the pair
    //   f2HistMtALamALamReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
    // 				    fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2]));

    //   // Fill the qinv, minv histogram 
    //   f3HistALamALamQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2]),minDistSftAPro,minDistSftPio);
    //   f3HistALamALamMinvReal->Fill(Minv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2]),minDistSftAPro,minDistSftPio);

    //   // Mt of the pair for low q pairs only
    //   if(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
    // 	      fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2]) < .2)
    // 	f2HistMtLowQALamALamReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
    // 				      fFemtoBuffer->GetEvt(0)->fALamTracks[iALam2]));

    // } // Second lambda loop

    // AProton loop
    for (iAPro=0;iAPro<fFemtoBuffer->GetEvt(0)->GetNAPro();iAPro++){

      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
  	continue;

      //      nALamAProPairsWoTTR++;
      
      // Reset the distances for each pair
      //      minDistAPro=999.0;meanDistAPro=0.0;
      minDistSftAPro=999.0;meanDistSftAPro=0.0;
      // minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;
      // minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;
      // Reset the number of measurements for the mean
      nMeasAPro=0;//nMeasAProIroc=0;nMeasAProOroc=0;

      // Check for two-track resolution
      for (iRadius=0;iRadius<9;iRadius++){
  	// Get the spatial distance at each radius
	//  	distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].fXglobal[iRadius]);
  	distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
        if (distSftAPro > -1.0){
          // Minimum distance
          // if (distAPro < minDistAPro)
          //   minDistAPro = distAPro;
	  if (distSftAPro < minDistSftAPro)
            minDistSftAPro = distSftAPro;
	  // Mean distance
	  //          meanDistAPro+=distAPro;
	  meanDistSftAPro+=distSftAPro;
	  nMeasAPro++;

	//   // IROC                                                             
        //   if (iRadius<3){
        //     if (distSftAPro < minDistSftIrocAPro)
        //       minDistSftIrocAPro = distSftAPro;
        //     meanDistSftIrocAPro+=distSftAPro;
        //     nMeasAProIroc++;
        //   }
        //   // OROC                                                             
        //   else {
        //     if (distSftAPro < minDistSftOrocAPro)
        //       minDistSftOrocAPro = distSftAPro;
        //     meanDistSftOrocAPro+=distSftAPro;
        //     nMeasAProOroc++;
        //   }

	}
      } // Loop over iRadius
      
      // Require at least one measurement
      if ( !nMeasAPro )
        continue;

      // Divide by the number of measurements to get the mean
      //      meanDistAPro /= (Float_t)nMeasAPro;
      meanDistSftAPro /= (Float_t)nMeasAPro;

      // Fill the two track resolution histogram
      // f2HistALamAProMeanMinDistAProReal->Fill(meanDistAPro,minDistAPro);
      // f2HistSftALamAProMeanMinDistAProReal->Fill(meanDistSftAPro,minDistSftAPro);

      // // Fill IROC / OROC histograms only with at least one measurement
      // if (nMeasAProIroc){
      //   meanDistSftIrocAPro /= (Float_t)nMeasAProIroc;
      //   f2HistSftIrocALamAProMeanMinDistAProReal->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
      // }
      // if (nMeasAProOroc){
      //   meanDistSftOrocAPro /= (Float_t)nMeasAProOroc;
      //   f2HistSftOrocALamAProMeanMinDistAProReal->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
      // }
  	
      // Look at possible residual correlations
      //      f2HistALamPosDAProQinvReal->Fill(QinvPioPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter, fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro);
      //      f2HistALamNegDAProQinvReal->Fill(QinvProPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter, fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro);

      
      //      // Do a cut (value needs to be refined)
      //      if ( meanDistSftAPro < 2.0 )
      //  	continue;

      // Count the number of pairs
      //      nALamAProPairs++;

      // Mt of the pair
      f2HistMtALamAProReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
				    fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]));
      // Use THnSparse since sept '12
      Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter, 
		   fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),
	meanDistSftAPro,
	minDistSftAPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
	     fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro])
      };
      ALamAProReal->Fill(x);

      // Fill the qinv histogram using TPC only momentum for the primary protons
      //      f3HistALamAProQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro); // Using THnSparse since Sept '12, see above
      //      f3HistALamAProMinvReal->Fill(Minv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);
      // and using TPC only plus primary vertex constraint momentum for primary protons
      //      f3HistALamAProQinvConstrReal->Fill(QinvConstr(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);

      // Mt of the pair for low q pairs only
      if(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
	      fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]) <.2)
	f2HistMtLowQALamAProReal->Fill(mt(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
					  fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]));

    }// AProton loop

  }// ALambda loop


  // Print the number of found pairs in this event
  // printf ("w/o  TTR cut: LamLamPairs: %d, LamProPairs: %d, ALamALamPairs: %d, ALamAProPairs: %d\n",
  // 	  nLamLamPairsWoTTR,nLamProPairsWoTTR,nALamALamPairsWoTTR,nALamAProPairsWoTTR);
  // printf ("with TTR cut: LamLamPairs: %d, LamProPairs: %d, ALamALamPairs: %d, ALamAProPairs: %d\n",
  // 	  nLamLamPairs,nLamProPairs,nALamALamPairs,nALamAProPairs);
  

}

//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessMixed() {
  // Process mixed events

  // Declare numbers to speed up the code
  Int_t iLam,//iLam2,
    iRadius,iPro,iALam,
    //iALam2,
    iAPro,nMeasPro,
    //nMeasPio,
    nMeasAPro;
  //  Int_t nMeasProIroc,nMeasPioIroc,nMeasAProIroc,nMeasProOroc,nMeasPioOroc,nMeasAProOroc;
  // Float_t distPro,distPio,minDistPro,meanDistPro,minDistPio,meanDistPio,
  //   distAPro,minDistAPro,meanDistAPro;
  Float_t distSftPro,
    //distSftPio,
    minDistSftPro,meanDistSftPro,//minDistSftPio,meanDistSftPio,
    distSftAPro,minDistSftAPro,meanDistSftAPro;
  // Float_t minDistSftIrocPro,meanDistSftIrocPro,minDistSftIrocPio,meanDistSftIrocPio,
  //   minDistSftIrocAPro,meanDistSftIrocAPro;
  // Float_t minDistSftOrocPro,meanDistSftOrocPro,minDistSftOrocPio,meanDistSftOrocPio,
  //   minDistSftOrocAPro,meanDistSftOrocAPro;
  
  // Loop over the event buffer
  for (UChar_t iMix = 1;iMix<fFemtoBuffer->GetMixBuffSize();iMix++){
    // Lambda loop
    for (iLam = 0; iLam < fFemtoBuffer->GetEvt(0)->GetNLam(); iLam++){
      
      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].UseIt())
  	continue;
      
      // // Second lambda loop
      // for (iLam2 = 0; iLam2 < (fFemtoBuffer->GetEvt(iMix))->GetNLam(); iLam2++){
	
      // 	// Skip if unUseIt() entry
      // 	if (!(fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2].UseIt())
      // 	  continue;
	
      // 	// Reset the distances for each pair
      // 	minDistPro=999.0;meanDistPro=0.0;minDistPio=999.0;meanDistPio=0.0;
      // 	minDistSftPro=999.0;meanDistSftPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
      // 	minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
      // 	minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
      // 	// Reset the number of measurements for the mean
      // 	nMeasPro=0;nMeasPio=0;nMeasProIroc=0;nMeasPioIroc=0;nMeasProOroc=0;nMeasPioOroc=0;
	
      // 	// Check for two-track resolution
      // 	for (iRadius=0;iRadius<9;iRadius++){
      // 	  // Get the spatial distance at each radius
      // 	  distPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2].fPosDaughter.fXglobal[iRadius]);
      // 	  distPio = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2].fNegDaughter.fXglobal[iRadius]);
      // 	  // Shifted distances
      // 	  distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));
      // 	  distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

      // 	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
      // 	  if (distPro > -1.0){
      // 	    // Minimum distance
      // 	    if (distPro < minDistPro)
      // 	      minDistPro = distPro;
      // 	    if (distSftPro < minDistSftPro)
      // 	      minDistSftPro = distSftPro;
      // 	    // Mean distance
      // 	    meanDistPro+=distPro;
      // 	    meanDistSftPro+=distSftPro;
      // 	    nMeasPro++;
	    
      // 	    // IROC                                                         
      // 	    if (iRadius<3){
      // 	      if (distSftPro < minDistSftIrocPro)
      // 		minDistSftIrocPro = distSftPro;
      // 	      meanDistSftIrocPro+=distSftPro;
      // 	      nMeasProIroc++;
      // 	    }
      // 	    // OROC                                                         
      // 	    else {
      // 	      if (distSftPro < minDistSftOrocPro)
      // 		minDistSftOrocPro = distSftPro;
      // 	      meanDistSftOrocPro+=distSftPro;
      // 	      nMeasProOroc++;
      // 	    }

      // 	  }
      // 	  if (distPio > -1.0){
      // 	    // Minimum distance
      // 	    if (distPio < minDistPio)
      // 	      minDistPio = distPio;
      // 	    if (distSftPio < minDistSftPio)
      // 	      minDistSftPio = distSftPio;
      // 	    // Mean distance
      // 	    meanDistPio+=distPio;
      // 	    meanDistSftPio+=distSftPio;
      // 	    nMeasPio++; 

      // 	    // IROC                                                         
      // 	    if (iRadius<3){
      // 	      if (distSftPio < minDistSftIrocPio)
      // 		minDistSftIrocPio = distSftPio;
      // 	      meanDistSftIrocPio+=distSftPio;
      // 	      nMeasPioIroc++;
      // 	    }
      // 	    // OROC                                                         
      // 	    else {
      // 	      if (distSftPio < minDistSftOrocPio)
      // 		minDistSftOrocPio = distSftPio;
      // 	      meanDistSftOrocPio+=distSftPio;
      // 	      nMeasPioOroc++;
      // 	    }
	    
      // 	  }
	  
      // 	} // Loop over iRadius
	
      // 	// Require at least one measurement
      // 	if ( (!nMeasPio) || (!nMeasPro) )
      // 	  continue;
	
      // 	// Divide by the number of measurements to get the mean
      // 	meanDistPro /= (Float_t)nMeasPro;
      // 	meanDistPio /= (Float_t)nMeasPio;
      // 	meanDistSftPro /= (Float_t)nMeasPro;
      // 	meanDistSftPio /= (Float_t)nMeasPio;
	
      // 	// Fill the two track resolution histograms
      // 	f2HistLamLamMeanMinDistProMixed->Fill(meanDistPro,minDistPro);
      // 	f2HistLamLamMeanMinDistPioMixed->Fill(meanDistPio,minDistPio);

      // 	f2HistSftLamLamMeanMinDistProMixed->Fill(meanDistSftPro,minDistSftPro);
      // 	f2HistSftLamLamMeanMinDistPioMixed->Fill(meanDistSftPio,minDistSftPio);

      // 	// Fill IROC / OROC histograms only with at least one measurement   
      // 	if (nMeasProIroc){
      // 	  meanDistSftIrocPro /= (Float_t)nMeasProIroc;
      // 	  f2HistSftIrocLamLamMeanMinDistProMixed->Fill(meanDistSftIrocPro,minDistSftIrocPro);
      // 	}
      // 	if (nMeasPioIroc){
      // 	  meanDistSftIrocPio /= (Float_t)nMeasPioIroc;
      // 	  f2HistSftIrocLamLamMeanMinDistPioMixed->Fill(meanDistSftIrocPio,minDistSftIrocPio);
      // 	}
      // 	if (nMeasProOroc){
      // 	  meanDistSftOrocPro /= (Float_t)nMeasProOroc;
      // 	  f2HistSftOrocLamLamMeanMinDistProMixed->Fill(meanDistSftOrocPro,minDistSftOrocPro);
      // 	}
      // 	if (nMeasPioOroc){
      // 	  meanDistSftOrocPio /= (Float_t)nMeasPioOroc;
      // 	  f2HistSftOrocLamLamMeanMinDistPioMixed->Fill(meanDistSftOrocPio,minDistSftOrocPio);
      // 	}
     	
      // 	//  	// Do a cut (value needs to be refined)
      // 	//  	if ( meanDistSftPro < 2.0 || meanDistSftPio < 2.0 )
      // 	//  	  continue;
	
      // 	// Fill the qinv, minv histogram
      // 	f3HistLamLamQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], (fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2]),minDistSftPro,minDistSftPio);
      // 	f3HistLamLamMinvMixed->Fill(Minv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], (fFemtoBuffer->GetEvt(iMix))->fLamTracks[iLam2]),minDistSftPro,minDistSftPio);
      
      // } // Second lambda loop
      
      // Proton loop
      for (iPro=0;iPro<(fFemtoBuffer->GetEvt(iMix))->GetNPro();iPro++){
	
  	// Skip if unUseIt() entry
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].UseIt())
  	  continue;
	
	// Reset the distances for each pair
	//	minDistPro=999.0;meanDistPro=0.0;	
	minDistSftPro=999.0;meanDistSftPro=0.0;
	// minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;
	// minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;
	// Reset the number of measurements for the mean
	nMeasPro=0;//nMeasProIroc=0;nMeasProOroc=0;

  	// Check for two-track resolution
  	for (iRadius=0;iRadius<9;iRadius++){
  	  // Get the spatial distance at each radius
	  //  	  distPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].fXglobal[iRadius]);
  	  distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
	  if (distSftPro > -1.0){
	    // Minimum distance
	    // if (distPro < minDistPro)
	    //   minDistPro = distPro;
	    if (distSftPro < minDistSftPro)
            minDistSftPro = distSftPro;
	    // Mean distance
	    //	    meanDistPro+=distPro;
	    meanDistSftPro+=distSftPro;
	    nMeasPro++;

	  //   // IROC                                                         
	  //   if (iRadius<3){
	  //     if (distSftPro < minDistSftIrocPro)
	  // 	minDistSftIrocPro = distSftPro;
	  //     meanDistSftIrocPro+=distSftPro;
	  //     nMeasProIroc++;
	  //   }
	  //   // OROC                                                         
	  //   else {
	  //     if (distSftPro < minDistSftOrocPro)
	  // 	minDistSftOrocPro = distSftPro;
	  //     meanDistSftOrocPro+=distSftPro;
	  //     nMeasProOroc++;
	  //   }

	  }
	} // Loop over iRadius
	
	// Require at least one measurement
	if ( !nMeasPro )
	  continue;

	// Divide by the number of measurements to get the mean
	//	meanDistPro /= (Float_t)nMeasPro;
	meanDistSftPro /= (Float_t)nMeasPro;
	
  	// // Fill the two track resolution histogram
  	// f2HistLamProMeanMinDistProMixed->Fill(meanDistPro,minDistPro);
	// f2HistSftLamProMeanMinDistProMixed->Fill(meanDistSftPro,minDistSftPro);

	// // Fill IROC / OROC histograms only with at least one measurement
	// if (nMeasProIroc){
	//   meanDistSftIrocPro /= (Float_t)nMeasProIroc;
	//   f2HistSftIrocLamProMeanMinDistProMixed->Fill(meanDistSftIrocPro,minDistSftIrocPro);
	// }
	// if (nMeasProOroc){
	//   meanDistSftOrocPro /= (Float_t)nMeasProOroc;
	//   f2HistSftOrocLamProMeanMinDistProMixed->Fill(meanDistSftOrocPro,minDistSftOrocPro);
	// }

	// Look at possible residual correlations with the daughters 
	//	f2HistLamPosDProQinvMixed->Fill(QinvProPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter, fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]),meanDistSftPro);
	//	f2HistLamNegDProQinvMixed->Fill(QinvPioPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fNegDaughter, fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]),meanDistSftPro);

  	
	//  	// Do a cut (value needs to be refined)
	//  	if ( meanDistSftPro < 2.0 )
	//  	  continue;

	// Use THnSparse since sept '12
	Double_t x[4]={
	  QinvProPro(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam].fPosDaughter,
		     fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]),
	  meanDistSftPro,
	  minDistSftPro,
	  Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam],
	       fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro])
	};
	LamProMixed->Fill(x);

  	// Fill the qinv histogram using TPC only momentum for the primary protons
      //  	f3HistLamProQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], (fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro]),meanDistSftPro,minDistSftPro);
	//  	f3HistLamProMinvMixed->Fill(Minv(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], (fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro]),meanDistSftPro,minDistSftPro);
	// and using TPC only plus primary vertex constraint for the primary proton
	//	f3HistLamProQinvConstrMixed->Fill(QinvConstr(fFemtoBuffer->GetEvt(0)->fLamTracks[iLam], (fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro]),meanDistSftPro,minDistSftPro);

      }// Proton loop
    }// Lambda loop
    
    
    // Anti-lambda loop
    for (iALam = 0; iALam < fFemtoBuffer->GetEvt(0)->GetNALam(); iALam++){
      
      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].UseIt())
  	continue;
      
      // // Second anti-lambda loop
      // for (iALam2 = 0; iALam2 < (fFemtoBuffer->GetEvt(iMix))->GetNALam(); iALam2++){
	
      // 	// Skip if unUseIt() entry
      // 	if (!(fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2].UseIt())
      // 	  continue;
	
      // 	// Reset the distances for each pair
      // 	minDistAPro=999.0;meanDistAPro=0.0;minDistPio=999.0;meanDistPio=0.0;
      // 	minDistSftAPro=999.0;meanDistSftAPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
      // 	minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
      // 	minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
      // 	// Reset the number of measurements for the mean
      // 	nMeasAPro=0;nMeasPio=0;nMeasAProIroc=0;nMeasPioIroc=0;nMeasAProOroc=0;nMeasPioOroc=0;

      // 	// Check for two-track resolution
      // 	for (iRadius=0;iRadius<9;iRadius++){
      // 	  // Get the spatial distance at each radius
      // 	  distPio = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2].fPosDaughter.fXglobal[iRadius]);
      // 	  distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2].fNegDaughter.fXglobal[iRadius]);
      // 	  // Shifted distances                                              
      // 	  distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));
      // 	  distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));


      // 	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
      // 	  if (distAPro > -1.0){
      // 	    // Minimum distance
      // 	    if (distAPro < minDistAPro)
      // 	      minDistAPro = distAPro;
      // 	    if (distSftAPro < minDistSftAPro)
      // 	      minDistSftAPro = distSftAPro;
      // 	    // Mean distance
      // 	    meanDistAPro+=distAPro;
      // 	    meanDistSftAPro+=distSftAPro;
      // 	    nMeasAPro++;
      // 	    // IROC                                                         
      // 	    if (iRadius<3){
      // 	      if (distSftAPro < minDistSftIrocAPro)
      // 		minDistSftIrocAPro = distSftAPro;
      // 	      meanDistSftIrocAPro+=distSftAPro;
      // 	      nMeasAProIroc++;
      // 	    }
      // 	    // OROC	    
      // 	    else {
      // 	      if (distSftAPro < minDistSftOrocAPro)
      // 		minDistSftOrocAPro = distSftAPro;
      // 	      meanDistSftOrocAPro+=distSftAPro;
      // 	      nMeasAProOroc++;
      // 	    }
      // 	  }
      // 	  if (distPio > -1.0){
      //     // Minimum distance
      // 	    if (distPio < minDistPio)
      // 	      minDistPio = distPio;
      // 	    if (distSftPio < minDistSftPio)
      // 	      minDistSftPio = distSftPio;
      // 	    // Mean distance
      // 	    meanDistPio+=distPio;
      // 	    meanDistSftPio+=distSftPio;
      // 	    nMeasPio++; 
      // 	    // IROC
      // 	    if (iRadius<3){
      // 	      if (distSftPio < minDistSftIrocPio)
      // 		minDistSftIrocPio = distSftPio;
      // 	      meanDistSftIrocPio+=distSftPio;
      // 	      nMeasPioIroc++;
      // 	    }
      // 	    // OROC
      // 	    else {
      // 	      if (distSftPio < minDistSftOrocPio)
      // 		minDistSftOrocPio = distSftPio;
      // 	      meanDistSftOrocPio+=distSftPio;
      // 	      nMeasPioOroc++;
      // 	    }
      // 	  }
      // 	} // Loop over iRadius
	
      // 	// Require at least one measurement
      // 	if ( (!nMeasPio) || (!nMeasAPro) )
      // 	  continue;
	
      // 	// Divide by the number of measurements to get the mean
      // 	meanDistAPro /= (Float_t)nMeasAPro;
      // 	meanDistPio /= (Float_t)nMeasPio;
      // 	meanDistSftAPro /= (Float_t)nMeasAPro;
      // 	meanDistSftPio /= (Float_t)nMeasPio;
	
      // 	// Fill the two track resolution histograms
      // 	f2HistALamALamMeanMinDistAProMixed->Fill(meanDistAPro,minDistAPro);
      // 	f2HistALamALamMeanMinDistPioMixed->Fill(meanDistPio,minDistPio);
      // 	f2HistSftALamALamMeanMinDistAProMixed->Fill(meanDistSftAPro,minDistSftAPro);
      // 	f2HistSftALamALamMeanMinDistPioMixed->Fill(meanDistSftPio,minDistSftPio);

      // 	// Fill IROC / OROC histograms only with at least one measurement
      // 	if (nMeasAProIroc){
      // 	  meanDistSftIrocAPro /= (Float_t)nMeasAProIroc;
      // 	  f2HistSftIrocALamALamMeanMinDistAProMixed->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
      // 	}
      // 	if (nMeasPioIroc){
      // 	  meanDistSftIrocPio /= (Float_t)nMeasPioIroc;
      // 	  f2HistSftIrocALamALamMeanMinDistPioMixed->Fill(meanDistSftIrocPio,minDistSftIrocPio);
      // 	}
      // 	if (nMeasAProOroc){
      // 	  meanDistSftOrocAPro /= (Float_t)nMeasAProOroc;
      // 	  f2HistSftOrocALamALamMeanMinDistAProMixed->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
      // 	}
      // 	if (nMeasPioOroc){
      // 	  meanDistSftOrocPio /= (Float_t)nMeasPioOroc;
      // 	  f2HistSftOrocALamALamMeanMinDistPioMixed->Fill(meanDistSftOrocPio,minDistSftOrocPio);
      // 	}
      // 	//  	// Do a cut (value needs to be refined)
      // 	//  	if ( meanDistSftAPro < 2.0 || meanDistSftPio < 2.0 )
      // 	//  	  continue;

      // 	// Fill the qinv, minv histogram
      // 	f3HistALamALamQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], (fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2]),minDistSftAPro,minDistSftPio);
      // 	f3HistALamALamMinvMixed->Fill(Minv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], (fFemtoBuffer->GetEvt(iMix))->fALamTracks[iALam2]),minDistSftAPro,minDistSftPio);
	
      // } // Second lambda loop

      // AProton loop
      for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(iMix))->GetNAPro();iAPro++){
	
  	// Skip if unUseIt() entry
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].UseIt())
  	  continue;
	
	// Reset the distances for each pair
	//	minDistAPro=999.0;meanDistAPro=0.0;	
	minDistSftAPro=999.0;meanDistSftAPro=0.0;
	// minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;
	// minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;
	// Reset the number of measurements for the mean
	nMeasAPro=0;//nMeasAProIroc=0;nMeasAProOroc=0;

  	// Check for two-track resolution
  	for (iRadius=0;iRadius<9;iRadius++){
  	  // Get the spatial distance at each radius
	  //  	  distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].fXglobal[iRadius]);
  	  distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
	  if (distSftAPro > -1.0){
	    // Minimum distance
	    // if (distAPro < minDistAPro)
	    //   minDistAPro = distAPro;
	    if (distSftAPro < minDistSftAPro)
	      minDistSftAPro = distSftAPro;
	    // Mean distance
	    //	    meanDistAPro+=distAPro;
	    meanDistSftAPro+=distSftAPro;
	    nMeasAPro++;
	  //   // IROC                                                                    
          //   if (iRadius<3){
          //     if (distSftAPro < minDistSftIrocAPro)
          //       minDistSftIrocAPro = distSftAPro;
          //     meanDistSftIrocAPro+=distSftAPro;
          //     nMeasAProIroc++;
          //   }
          //   // OROC                                                                    
          //   else {
          //     if (distSftAPro < minDistSftOrocAPro)
          //       minDistSftOrocAPro = distSftAPro;
          //     meanDistSftOrocAPro+=distSftAPro;
          //     nMeasAProOroc++;
          //   }

	  }
	} // Loop over iRadius
      
	// Require at least one measurement
	if ( !nMeasAPro )
	  continue;
	
	// Divide by the number of measurements to get the mean
	//	meanDistAPro /= (Float_t)nMeasAPro;
	meanDistSftAPro /= (Float_t)nMeasAPro;
	
  	// // Fill the two track resolution histogram
  	// f2HistALamAProMeanMinDistAProMixed->Fill(meanDistAPro,minDistAPro);
	// f2HistSftALamAProMeanMinDistAProMixed->Fill(meanDistSftAPro,minDistSftAPro);

	// // Fill IROC / OROC histograms only with at least one measurement   
	// if (nMeasAProIroc){
	//   meanDistSftAPro /= (Float_t)nMeasAProIroc;
	//   f2HistSftIrocALamAProMeanMinDistAProMixed->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
	// }
	// if (nMeasAProOroc){
	//   meanDistSftAPro /= (Float_t)nMeasAProOroc;
	//   f2HistSftOrocALamAProMeanMinDistAProMixed->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
	// }

	// Look at possible residual correlations for the daughters
	//	f2HistALamPosDAProQinvMixed->Fill(QinvPioPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fPosDaughter, fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]),meanDistSftAPro);
	//	f2HistALamNegDAProQinvMixed->Fill(QinvProPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter, fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]),meanDistSftAPro);

  	// Do a cut (value needs to be refined)
	//  	if ( meanDistSftAPro < 2.0 )
	//  	  continue;
	
	// Use THnSparse since sept '12
	Double_t x[4]={
	  QinvProPro(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam].fNegDaughter,
		     fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]),
	  meanDistSftAPro,
	  minDistSftAPro,
	  Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam],
	       fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro])
	};
	ALamAProMixed->Fill(x);


  	// Fill the qinv histogram using TPC only momentum for the primary proton
	//  	f3HistALamAProQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], (fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);
	//  	f3HistALamAProMinvMixed->Fill(Minv(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], (fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);
	// and using TPC only plus primary vertex constraint for the primary proton
	//  	f3HistALamAProQinvConstrMixed->Fill(QinvConstr(fFemtoBuffer->GetEvt(0)->fALamTracks[iALam], (fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);
	
      }// AProton loop
    }// ALambda loop
    
  }// Event buffer loop

}// End of void ProcessMixed 
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessRealBackground() {
  // Process real events with background lambdas
  
  // Declare numbers to speed up the code
  Int_t iBgLam,
    //iBgLam2,
    iRadius,iPro,iBgALam,
    //iBgALam2,
    iAPro,nMeasPro,
    //nMeasPio,
    nMeasAPro;
  //  Int_t nMeasProIroc,nMeasPioIroc,nMeasAProIroc,nMeasProOroc,nMeasPioOroc,nMeasAProOroc;
  // Float_t distPro,distPio,minDistPro=999.0,meanDistPro=0.0,minDistPio=999.0,meanDistPio=0.0,
  //   distAPro,minDistAPro=999.0,meanDistAPro=0.0;
  Float_t distSftPro,
    //distSftPio,
    minDistSftPro=999.0,meanDistSftPro=0.0,
    //minDistSftPio=999.0,meanDistSftPio=0.0,
    distSftAPro,minDistSftAPro=999.0,meanDistSftAPro=0.0;
  // Float_t minDistSftIrocPro=999.0,meanDistSftIrocPro=0.0,minDistSftIrocPio=999.0,meanDistSftIrocPio=0.0,
  //   minDistSftIrocAPro=999.0,meanDistSftIrocAPro=0.0;
  // Float_t minDistSftOrocPro=999.0,meanDistSftOrocPro=0.0,minDistSftOrocPio=999.0,meanDistSftOrocPio=0.0,
  //   minDistSftOrocAPro=999.0,meanDistSftOrocAPro=0.0;
  
      
  // BgLambda loop
  for (iBgLam = 0; iBgLam < fFemtoBuffer->GetEvt(0)->GetNBgLam(); iBgLam++){

    // Skip if unUseIt() entry
    if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].UseIt())
      continue;

    // // Second lambda loop
    // for (iBgLam2 = iBgLam+1; iBgLam2 < fFemtoBuffer->GetEvt(0)->GetNBgLam(); iBgLam2++){

    //   // Skip if unUseIt() entry
    //   if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2].UseIt())
    // 	continue;

    //   // Reset the distances for each pair
    //   minDistPro=999.0;meanDistPro=0.0;minDistPio=999.0;meanDistPio=0.0;
    //   minDistSftPro=999.0;meanDistSftPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
    //   minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
    //   minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
    //   // Reset the number of measurements for the mean
    //   nMeasPro=0;nMeasPio=0;nMeasProIroc=0;nMeasPioIroc=0;nMeasProOroc=0;nMeasPioOroc=0;

    //   // Check for two-track resolution
    //   for (iRadius=0;iRadius<9;iRadius++){
    // 	// Get the spatial distance at each radius
    // 	distPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2].fPosDaughter.fXglobal[iRadius]);
    // 	distPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2].fNegDaughter.fXglobal[iRadius]);
    // 	// Shifted distances
    // 	distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));
    // 	distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));
	
    // 	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
    //     if (distPro > -1.0){
    //       // Minimum distance
    //       if (distPro < minDistPro)
    //         minDistPro = distPro;
    // 	  if (distSftPro < minDistSftPro)
    //         minDistSftPro = distSftPro;
    //       // Mean distance
    //       meanDistPro+=distPro;
    // 	  meanDistSftPro+=distSftPro;
    //       nMeasPro++;

    // 	  // IROC                                                         
    //       if (iRadius<3){
    //         if (distSftPro < minDistSftIrocPro)
    //           minDistSftIrocPro = distSftPro;
    //         meanDistSftIrocPro+=distSftPro;
    //         nMeasProIroc++;
    //       }
    //       // OROC                                                         
    //       else {
    //         if (distSftPro < minDistSftOrocPro)
    //           minDistSftOrocPro = distSftPro;
    //         meanDistSftOrocPro+=distSftPro;
    //         nMeasProOroc++;
    //       }
    //     }
	
    //     if (distPio > -1.0){
    //       // Minimum distance
    //       if (distPio < minDistPio)
    //         minDistPio = distPio;
    // 	  if (distSftPio < minDistSftPio)
    //         minDistSftPio = distSftPio;
    //       // Mean distance
    //       meanDistPio+=distPio;
    // 	  meanDistSftPio+=distSftPio;
    //       nMeasPio++; 
    // 	  // IROC                                                         
    //       if (iRadius<3){
    //         if (distSftPio < minDistSftIrocPio)
    //           minDistSftIrocPio = distSftPio;
    //         meanDistSftIrocPio+=distSftPio;
    //         nMeasPioIroc++;
    //       }
    //       // OROC                                                         
    //       else {
    //         if (distSftPio < minDistSftOrocPio)
    //           minDistSftOrocPio = distSftPio;
    //         meanDistSftOrocPio+=distSftPio;
    //         nMeasPioOroc++;
    //       }
    //     }

    //   } // Loop over iRadius


      
    //   // Require at least one measurement
    //   if ( (!nMeasPio) || (!nMeasPro) )
    //     continue;
      
    //   // Divide by the number of measurements to get the mean
    //   meanDistPro /= (Float_t)nMeasPro;
    //   meanDistPio /= (Float_t)nMeasPio;
    //   meanDistSftPro /= (Float_t)nMeasPro;
    //   meanDistSftPio /= (Float_t)nMeasPio;

    //   // Fill the two track resolution histograms
    //   f2HistBgLamBgLamMeanMinDistProReal->Fill(meanDistPro,minDistPro);
    //   f2HistBgLamBgLamMeanMinDistPioReal->Fill(meanDistPio,minDistPio);
    //   f2HistSftBgLamBgLamMeanMinDistProReal->Fill(meanDistSftPro,minDistSftPro);
    //   f2HistSftBgLamBgLamMeanMinDistPioReal->Fill(meanDistSftPio,minDistSftPio);

    //   // Fill IROC / OROC histograms only with at least one measurement   
    //   if (nMeasProIroc){
    //     meanDistSftIrocPro /= (Float_t)nMeasProIroc;
    //     f2HistSftIrocBgLamBgLamMeanMinDistProReal->Fill(meanDistSftIrocPro,minDistSftIrocPro);
    //   }
    //   if (nMeasPioIroc){
    //     meanDistSftIrocPio /= (Float_t)nMeasPioIroc;
    //     f2HistSftIrocBgLamBgLamMeanMinDistPioReal->Fill(meanDistSftIrocPio,minDistSftIrocPio);
    //   }
    //   if (nMeasProOroc){
    //     meanDistSftOrocPro /= (Float_t)nMeasProOroc;
    //     f2HistSftOrocBgLamBgLamMeanMinDistProReal->Fill(meanDistSftOrocPro,minDistSftOrocPro);
    //   }
    //   if (nMeasPioOroc){
    //     meanDistSftOrocPio /= (Float_t)nMeasPioOroc;
    //     f2HistSftOrocBgLamBgLamMeanMinDistPioReal->Fill(meanDistSftOrocPio,minDistSftOrocPio);
    //   }


     	
    //   //      // Do a cut (value needs to be refined)
    // 	//      if ( meanDistSftPro < 2.0 || meanDistSftPio < 2.0 )
    //   //	continue;
      
    //   // Fill the qinv, minv histogram
    //   f3HistBgLamBgLamQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam], fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam2]),minDistSftPro,minDistSftPio);

    // } // Second lambda loop


    // Proton loop
    for (iPro=0;iPro<fFemtoBuffer->GetEvt(0)->GetNPro();iPro++){

      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fProTracks[iPro].UseIt())
  	continue;

      // Reset the distances for each pair
      //      minDistPro=999.0;meanDistPro=0.0;
      minDistSftPro=999.0;meanDistSftPro=0.0;
      //      minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;
      //      minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;
      // Reset the number of measurements for the mean
      nMeasPro=0;//nMeasProIroc=0;nMeasProOroc=0;

      // Check for two-track resolution
      for (iRadius=0;iRadius<9;iRadius++){
  	// Get the spatial distance at each radius
	//  	distPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fProTracks[iPro].fXglobal[iRadius]);
  	distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fProTracks[iPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
        if (distSftPro > -1.0){
          // Minimum distance
          // if (distPro < minDistPro)
          //   minDistPro = distPro;
          if (distSftPro < minDistSftPro)
            minDistSftPro = distSftPro;
          // Mean distance
	  //          meanDistPro+=distPro;
	  meanDistSftPro+=distSftPro;
          nMeasPro++;
	  // // IROC                                                         
          // if (iRadius<3){
          //   if (distSftPro < minDistSftIrocPro)
          //     minDistSftIrocPro = distSftPro;
          //   meanDistSftIrocPro+=distSftPro;
          //   nMeasProIroc++;
          // }
          // // OROC                                                         
          // else {
          //   if (distSftPro < minDistSftOrocPro)
          //     minDistSftOrocPro = distSftPro;
          //   meanDistSftOrocPro+=distSftPro;
          //   nMeasProOroc++;
          // }

        }

      } // Loop over iRadius
      
      // Require at least one measurement
      if ( !nMeasPro )
        continue;

      // Divide by the number of measurements to get the mean
      //      meanDistPro /= (Float_t)nMeasPro;
      meanDistSftPro /= (Float_t)nMeasPro;

      // Fill the two track resolution histogram
      // f2HistBgLamProMeanMinDistProReal->Fill(meanDistPro,minDistPro);
      // f2HistSftBgLamProMeanMinDistProReal->Fill(meanDistSftPro,minDistSftPro);

      // // Fill IROC / OROC histograms only with at least one measurement
      // if (nMeasProIroc){
      //   meanDistSftIrocPro /= (Float_t)nMeasProIroc;
      //   f2HistSftIrocBgLamProMeanMinDistProReal->Fill(meanDistSftIrocPro,minDistSftIrocPro);
      // }
      // if (nMeasProOroc){
      //   meanDistSftOrocPro /= (Float_t)nMeasProOroc;
      //   f2HistSftOrocBgLamProMeanMinDistProReal->Fill(meanDistSftOrocPro,minDistSftOrocPro);
      // }

  	
      // Do a cut (value needs to be refined)
      //      if ( meanDistSftPro < 2.0 )
      //  	continue;

      // Since sept '12 do THnSparse with qinvpropro, 
      // mean dist propro, min dist propro, qinv lampro
      Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter,
		   fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),
	meanDistSftPro,
	minDistSftPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
	     fFemtoBuffer->GetEvt(0)->fProTracks[iPro])
      };
      BgLamProReal->Fill(x);

      // Fill the qinv histogram
      //      f3HistBgLamProQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam], fFemtoBuffer->GetEvt(0)->fProTracks[iPro]),meanDistSftPro,minDistSftPro);

    }// Proton loop
  }// BgLambda loop


  // Anti-lambda loop
  for (iBgALam = 0; iBgALam < fFemtoBuffer->GetEvt(0)->GetNBgALam(); iBgALam++){

    // Skip if unUseIt() entry
    if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].UseIt())
      continue;
    
    // // Second anti-lambda loop
    // for (iBgALam2 = iBgALam+1; iBgALam2 < fFemtoBuffer->GetEvt(0)->GetNBgALam(); iBgALam2++){

    //   // Skip if unUseIt() entry
    //   if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].UseIt())
    // 	continue;

    //   // Reset the distances for each pair
    //   minDistAPro=999.0;meanDistAPro=0.0;minDistPio=999.0;meanDistPio=0.0;
    //   minDistSftAPro=999.0;meanDistSftAPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
    //   minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
    //   minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
    //   // Reset the number of measurements for the mean
    //   nMeasAPro=0;nMeasPio=0;
    //   nMeasAProIroc=0;nMeasPioIroc=0;
    //   nMeasAProOroc=0;nMeasPioOroc=0;

    //   // Check for two-track resolution
    //   for (iRadius=0;iRadius<9;iRadius++){
    // 	// Get the spatial distance at each radius
    // 	distPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fPosDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fPosDaughter.fXglobal[iRadius]);
    // 	distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fNegDaughter.fXglobal[iRadius]);
    // 	// Shifted distances                                              
    //     distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fPosDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));
    //     distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

    // 	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
    //     if (distAPro > -1.0){
    //       // Minimum distance
    //       if (distAPro < minDistAPro)
    //         minDistAPro = distAPro;
    // 	  if (distSftAPro < minDistSftAPro)
    //         minDistSftAPro = distSftAPro;
    //       // Mean distance
    //       meanDistAPro+=distAPro;
    // 	  meanDistSftAPro+=distSftAPro;
    //       nMeasAPro++;

    // 	  // IROC                                                         
    //       if (iRadius<3){
    //         if (distSftAPro < minDistSftIrocAPro)
    //           minDistSftIrocAPro = distSftAPro;
    //         meanDistSftIrocAPro+=distSftAPro;
    //         nMeasAProIroc++;
    //       }
    //       // OROC                                                                          
    //       else {
    //         if (distSftAPro < minDistSftOrocAPro)
    //           minDistSftOrocAPro = distSftAPro;
    //         meanDistSftOrocAPro+=distSftAPro;
    //         nMeasAProOroc++;
    //       }

    //     }
    //     if (distPio > -1.0){
    //       // Minimum distance
    //       if (distPio < minDistPio)
    //         minDistPio = distPio;
    // 	  if (distSftPio < minDistSftPio)
    //         minDistSftPio = distSftPio;
    //       // Mean distance
    //       meanDistPio+=distPio;
    // 	  meanDistSftPio+=distSftPio;
    //       nMeasPio++; 
    // 	  // IROC                                                         
    //       if (iRadius<3){
    //         if (distSftPio < minDistSftIrocPio)
    //           minDistSftIrocPio = distSftPio;
    //         meanDistSftIrocPio+=distSftPio;
    //         nMeasPioIroc++;
    //       }
    //       // OROC                                                         
    //       else {
    //         if (distSftPio < minDistSftOrocPio)
    //           minDistSftOrocPio = distSftPio;
    //         meanDistSftOrocPio+=distSftPio;
    //         nMeasPioOroc++;
    //       }

    //     }

    //   } // Loop over iRadius
      
    //   // Require at least one measurement
    //   if ( (!nMeasPio) || (!nMeasAPro) )
    //     continue;

    //   // Divide by the number of measurements to get the mean
    //   meanDistAPro /= (Float_t)nMeasAPro;
    //   meanDistPio /= (Float_t)nMeasPio;
    //   meanDistSftAPro /= (Float_t)nMeasAPro;
    //   meanDistSftPio /= (Float_t)nMeasPio;

    //   // Fill the two track resolution histograms
    //   f2HistBgALamBgALamMeanMinDistAProReal->Fill(meanDistAPro,minDistAPro);
    //   f2HistBgALamBgALamMeanMinDistPioReal->Fill(meanDistPio,minDistPio);
    //   f2HistSftBgALamBgALamMeanMinDistAProReal->Fill(meanDistSftAPro,minDistSftAPro);
    //   f2HistSftBgALamBgALamMeanMinDistPioReal->Fill(meanDistSftPio,minDistSftPio);
    //   // Fill IROC / OROC histograms only with at least one measurement
    //   if (nMeasAProIroc){
    //     meanDistSftAPro /= (Float_t)nMeasAProIroc;
    //     f2HistSftIrocBgALamBgALamMeanMinDistAProReal->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
    //   }
    //   if (nMeasPioIroc){
    //     meanDistSftPio /= (Float_t)nMeasPioIroc;
    //     f2HistSftIrocBgALamBgALamMeanMinDistPioReal->Fill(meanDistSftIrocPio,minDistSftIrocPio);
    //   }
    //   if (nMeasAProOroc){
    //     meanDistSftAPro /= (Float_t)nMeasAProOroc;
    //     f2HistSftOrocBgALamBgALamMeanMinDistAProReal->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
    //   }
    //   if (nMeasPioOroc){
    //     meanDistSftPio /= (Float_t)nMeasPioOroc;
    //     f2HistSftOrocBgALamBgALamMeanMinDistPioReal->Fill(meanDistSftOrocPio,minDistSftOrocPio);
    //   }

     	
    //   //      // Do a cut (value needs to be refined)
    //   //      if ( meanDistSftAPro < 2.0 || meanDistSftPio < 2.0 )
    //   //  	continue;
      
    //   // Fill the qinv histogram
    //   f3HistBgALamBgALamQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2]),minDistSftAPro,minDistSftPio);
      
    // } // Second lambda loop


    // AProton loop
    for (iAPro=0;iAPro<fFemtoBuffer->GetEvt(0)->GetNAPro();iAPro++){

      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].UseIt())
  	continue;

      // Reset the distances for each pair
      //      minDistAPro=999.0;meanDistAPro=0.0;
      minDistSftAPro=999.0;meanDistSftAPro=0.0;
      // minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;
      // minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;
      // Reset the number of measurements for the mean
      nMeasAPro=0;//nMeasAProIroc=0;nMeasAProOroc=0;
      
      // Check for two-track resolution
      for (iRadius=0;iRadius<9;iRadius++){
  	// Get the spatial distance at each radius
	//  	distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXglobal[iRadius],fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].fXglobal[iRadius]);
	distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXshifted[iRadius],fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),fFemtoBuffer->GetEvt(0));

	// calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
        if (distSftAPro > -1.0){
          // Minimum distance
          // if (distAPro < minDistAPro)
          //   minDistAPro = distAPro;
	  if (distSftAPro < minDistSftAPro)
            minDistSftAPro = distSftAPro;
          // Mean distance
	  //          meanDistAPro+=distAPro;
	  meanDistSftAPro+=distSftAPro;
          nMeasAPro++;

	  // // IROC                                                         
          // if (iRadius<3){
          //   if (distSftAPro < minDistSftIrocAPro)
          //     minDistSftIrocAPro = distSftAPro;
          //   meanDistSftIrocAPro+=distSftAPro;
          //   nMeasAProIroc++;
          // }
          // // OROC                                                                          
          // else {
          //   if (distSftAPro < minDistSftOrocAPro)
          //     minDistSftOrocAPro = distSftAPro;
          //   meanDistSftOrocAPro+=distSftAPro;
          //   nMeasAProOroc++;
          // }
        }
      } // Loop over iRadius
      
      // Require at least one measurement
      if ( !nMeasAPro )
        continue;

      // Divide by the number of measurements to get the mean
      //      meanDistAPro /= (Float_t)nMeasAPro;
      meanDistSftAPro /= (Float_t)nMeasAPro;

      // // Fill the two track resolution histogram
      // f2HistBgALamAProMeanMinDistAProReal->Fill(meanDistAPro,minDistAPro);
      // f2HistSftBgALamAProMeanMinDistAProReal->Fill(meanDistSftAPro,minDistSftAPro);
      // // Fill IROC / OROC histograms only with at least one measurement   
      // if (nMeasAProIroc){
      //   meanDistSftAPro /= (Float_t)nMeasAProIroc;
      //   f2HistSftIrocBgALamAProMeanMinDistAProReal->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
      // }
      // if (nMeasAProOroc){
      //   meanDistSftAPro /= (Float_t)nMeasAProOroc;
      //   f2HistSftOrocBgALamAProMeanMinDistAProReal->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
      // }

      
      // Do a cut (value needs to be refined)
      //      if ( meanDistSftAPro < 2.0 )
      //  	continue;

      // Since sept '12 do THnSparse with qinvpropro, 
      // mean dist propro, min dist propro, qinv lampro
      Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter,
		   fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),
	meanDistSftAPro,
	minDistSftAPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], 
	     fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro])
      };
      BgALamAProReal->Fill(x);

      // Fill the qinv histogram
      //      f3HistBgALamAProQinvReal->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], fFemtoBuffer->GetEvt(0)->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);

    }// AProton loop
  }// BgALambda loop
} // End of void ProcessRealBackground

//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ProcessMixedBackground() {
  // Process mixed events

  // Declare numbers to speed up the code
  Int_t iBgLam,
    //iBgLam2,
    iRadius,iPro,iBgALam,
    //iBgALam2,
    iAPro,nMeasPro,
    //nMeasPio,
    nMeasAPro;
  //  Int_t nMeasProIroc,nMeasPioIroc,nMeasAProIroc,nMeasProOroc,nMeasPioOroc,nMeasAProOroc;
  // Float_t distPro,distPio,minDistPro,meanDistPro,minDistPio,meanDistPio,
  //   distAPro,minDistAPro,meanDistAPro;
  Float_t distSftPro,//distSftPio,
    minDistSftPro,meanDistSftPro,//minDistSftPio,meanDistSftPio,
    distSftAPro,minDistSftAPro,meanDistSftAPro;
  // Float_t minDistSftIrocPro,meanDistSftIrocPro,minDistSftIrocPio,meanDistSftIrocPio,
  //   minDistSftIrocAPro,meanDistSftIrocAPro;
  // Float_t minDistSftOrocPro,meanDistSftOrocPro,minDistSftOrocPio,meanDistSftOrocPio,
  //   minDistSftOrocAPro,meanDistSftOrocAPro;

  
  // Loop over the event buffer
  for (UChar_t iMix = 1;iMix<fFemtoBuffer->GetMixBuffSize();iMix++){
    // BgLambda loop
    for (iBgLam = 0; iBgLam < fFemtoBuffer->GetEvt(0)->GetNBgLam(); iBgLam++){

      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].UseIt())
  	continue;
      
      // // Second lambda loop
      // for (iBgLam2 = 0; iBgLam2 < (fFemtoBuffer->GetEvt(iMix))->GetNBgLam(); iBgLam2++){

      // 	// Skip if unUseIt() entry
      // 	if (!(fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2].UseIt())
      // 	  continue;
	
      // 	// Reset the distances for each pair
      // 	minDistPro=999.0;meanDistPro=0.0;minDistPio=999.0;meanDistPio=0.0;
      // 	minDistSftPro=999.0;meanDistSftPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
      //   minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
      //   minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;
      // 	// Reset the number of measurements for the mean
      // 	nMeasPro=0;nMeasPio=0;nMeasProIroc=0;nMeasPioIroc=0;nMeasProOroc=0;nMeasPioOroc=0;

      // 	// Check for two-track resolution
      // 	for (iRadius=0;iRadius<9;iRadius++){
      // 	  // Get the spatial distance at each radius
      // 	  distPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2].fPosDaughter.fXglobal[iRadius]);
      // 	  distPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fNegDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2].fNegDaughter.fXglobal[iRadius]);
      // 	  // Shifted distances                                            
      //     distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));
      //     distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fNegDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

      // 	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
      // 	  if (distPro > -1.0){
      // 	    // Minimum distance
      // 	    if (distPro < minDistPro)
      // 	      minDistPro = distPro;
      // 	    if (distSftPro < minDistSftPro)
      //         minDistSftPro = distSftPro;
      // 	    // Mean distance
      // 	    meanDistPro+=distPro;
      // 	    meanDistSftPro+=distSftPro;
      // 	    nMeasPro++;
      // 	    // IROC                                                                          
      //       if (iRadius<3){
      //         if (distSftPro < minDistSftIrocPro)
      //           minDistSftIrocPro = distSftPro;
      //         meanDistSftIrocPro+=distSftPro;
      //         nMeasProIroc++;
      //       }
      //       // OROC                                                                          
      //       else {
      //         if (distSftPro < minDistSftOrocPro)
      //           minDistSftOrocPro = distSftPro;
      //         meanDistSftOrocPro+=distSftPro;
      //         nMeasProOroc++;
      //       }
      // 	  }
      // 	  if (distPio > -1.0){
      // 	    // Minimum distance
      // 	    if (distPio < minDistPio)
      // 	      minDistPio = distPio;
      // 	    if (distSftPio < minDistSftPio)
      //         minDistSftPio = distSftPio;
      // 	    // Mean distance
      // 	    meanDistPio+=distPio;
      // 	    meanDistSftPio+=distSftPio;
      // 	    nMeasPio++; 
      // 	    // IROC                                                                          
      //       if (iRadius<3){
      //         if (distSftPio < minDistSftIrocPio)
      //           minDistSftIrocPio = distSftPio;
      //         meanDistSftIrocPio+=distSftPio;
      //         nMeasPioIroc++;
      //       }
      //       // OROC                                                                          
      //       else {
      //         if (distSftPio < minDistSftOrocPio)
      //           minDistSftOrocPio = distSftPio;
      //         meanDistSftOrocPio+=distSftPio;
      //         nMeasPioOroc++;
      //       }
      // 	  }	  
      // 	} // Loop over iRadius
	
      // 	// Require at least one measurement
      // 	if ( (!nMeasPio) || (!nMeasPro) )
      // 	  continue;
	
      // 	// Divide by the number of measurements to get the mean
      // 	meanDistPro /= (Float_t)nMeasPro;
      // 	meanDistPio /= (Float_t)nMeasPio;
      // 	meanDistSftPro /= (Float_t)nMeasPro;
      //   meanDistSftPio /= (Float_t)nMeasPio;
	
      // 	// Fill the two track resolution histograms
      // 	f2HistBgLamBgLamMeanMinDistProMixed->Fill(meanDistPro,minDistPro);
      // 	f2HistBgLamBgLamMeanMinDistPioMixed->Fill(meanDistPio,minDistPio);
      // 	f2HistSftBgLamBgLamMeanMinDistProMixed->Fill(meanDistSftPro,minDistSftPro);
      //   f2HistSftBgLamBgLamMeanMinDistPioMixed->Fill(meanDistSftPio,minDistSftPio);
      // 	// Fill IROC / OROC histograms only with at least one measurement
      //   if (nMeasProIroc){
      //     meanDistSftPro /= (Float_t)nMeasProIroc;
      //     f2HistSftIrocBgLamBgLamMeanMinDistProMixed->Fill(meanDistSftIrocPro,minDistSftIrocPro);
      //   }
      //   if (nMeasPioIroc){
      //     meanDistSftPio /= (Float_t)nMeasPioIroc;
      //     f2HistSftIrocBgLamBgLamMeanMinDistPioMixed->Fill(meanDistSftIrocPio,minDistSftIrocPio);
      //   }
      //   if (nMeasProOroc){
      //     meanDistSftPro /= (Float_t)nMeasProOroc;
      //     f2HistSftOrocBgLamBgLamMeanMinDistProMixed->Fill(meanDistSftOrocPro,minDistSftOrocPro);
      //   }
      //   if (nMeasPioOroc){
      //     meanDistSftPio /= (Float_t)nMeasPioOroc;
      //     f2HistSftOrocBgLamBgLamMeanMinDistPioMixed->Fill(meanDistSftOrocPio,minDistSftOrocPio);
      //   }

      // 	//  	// Do a cut (value needs to be refined)
      // 	//  	if ( meanDistSftPro < 2.0 || meanDistSftPio < 2.0 )
      // 	//  	  continue;
	
      // 	// Fill the qinv, minv histogram
      // 	f3HistBgLamBgLamQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam], (fFemtoBuffer->GetEvt(iMix))->fBgLamTracks[iBgLam2]),minDistSftPro,minDistSftPio);
      
      // } // Second lambda loop

      
      // Proton loop
      for (iPro=0;iPro<(fFemtoBuffer->GetEvt(iMix))->GetNPro();iPro++){
	
  	// Skip if unUseIt() entry
  	if (!(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].UseIt())
  	  continue;
	
	// Reset the distances for each pair
	//	minDistPro=999.0;meanDistPro=0.0;
	minDistSftPro=999.0;meanDistSftPro=0.0;
	// minDistSftIrocPro=999.0;meanDistSftIrocPro=0.0;
        // minDistSftOrocPro=999.0;meanDistSftOrocPro=0.0;
	// Reset the number of measurements for the mean
	nMeasPro=0;//nMeasProIroc=0;nMeasProOroc=0;

  	// Check for two-track resolution
  	for (iRadius=0;iRadius<9;iRadius++){
  	  // Get the spatial distance at each radius
	  //  	  distPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].fXglobal[iRadius]);
	  distSftPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter.fXshifted[iRadius],
				(fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

	  
	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
	  if (distSftPro > -1.0){
	    // Minimum distance
	    // if (distPro < minDistPro)
	    //   minDistPro = distPro;
	    if (distSftPro < minDistSftPro)
	      minDistSftPro = distSftPro;
	    // Mean distance
	    //	    meanDistPro+=distPro;
	    meanDistSftPro+=distSftPro;
	    nMeasPro++;
	    // // IROC                                                                    
            // if (iRadius<3){
            //   if (distSftPro < minDistSftIrocPro)
            //     minDistSftIrocPro = distSftPro;
            //   meanDistSftIrocPro+=distSftPro;
            //   nMeasProIroc++;
            // }
            // // OROC                                                                    
            // else {
            //   if (distSftPro < minDistSftOrocPro)
            //     minDistSftOrocPro = distSftPro;
            //   meanDistSftOrocPro+=distSftPro;
            //   nMeasProOroc++;
            // }
	  }
	} // Loop over iRadius
	
	// Require at least one measurement
	if ( !nMeasPro )
	  continue;

	// Divide by the number of measurements to get the mean
	//	meanDistPro /= (Float_t)nMeasPro;
	meanDistSftPro /= (Float_t)nMeasPro;
	
  	// // Fill the two track resolution histogram
  	// f2HistBgLamProMeanMinDistProMixed->Fill(meanDistPro,minDistPro);
	// f2HistSftBgLamProMeanMinDistProMixed->Fill(meanDistSftPro,minDistSftPro);
	// // Fill IROC / OROC histograms only with at least one measurement              
        // if (nMeasProIroc){
        //   meanDistSftIrocPro /= (Float_t)nMeasProIroc;
        //   f2HistSftIrocBgLamProMeanMinDistProMixed->Fill(meanDistSftIrocPro,minDistSftIrocPro);
        // }
        // if (nMeasProOroc){
        //   meanDistSftOrocPro /= (Float_t)nMeasProOroc;
        //   f2HistSftOrocBgLamProMeanMinDistProMixed->Fill(meanDistSftOrocPro,minDistSftOrocPro);
        // }

	//  	// Do a cut (value needs to be refined)
	//  	if ( meanDistSftPro < 2.0 )
	//  	  continue;
	
	      // Since sept '12 do THnSparse with qinvpropro, 
      // mean dist propro, min dist propro, qinv lampro
      Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam].fPosDaughter,
		   fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro]),
	meanDistSftPro,
	minDistSftPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam],
	     fFemtoBuffer->GetEvt(iMix)->fProTracks[iPro])
      };
      BgLamProMixed->Fill(x);


  	// Fill the qinv histogram
	//  	f3HistBgLamProQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgLamTracks[iBgLam], (fFemtoBuffer->GetEvt(iMix))->fProTracks[iPro]),meanDistSftPro,minDistSftPro);
	
      }// Proton loop
    }// BgLambda loop
    
    
    // Anti-lambda loop
    for (iBgALam = 0; iBgALam < fFemtoBuffer->GetEvt(0)->GetNBgALam(); iBgALam++){
      
      // Skip if unUseIt() entry
      if (!fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].UseIt())
  	continue;
      
      // // Second anti-lambda loop
      // for (iBgALam2 = 0; iBgALam2 < (fFemtoBuffer->GetEvt(iMix))->GetNBgALam(); iBgALam2++){
	
      // 	// Skip if unUseIt() entry
      // 	if (!(fFemtoBuffer->GetEvt(iMix))->fBgALamTracks[iBgALam2].UseIt())
      // 	  continue;

      // 	// Reset the distances for each pair
      // 	minDistAPro=999.0;meanDistAPro=0.0;minDistPio=999.0;meanDistPio=0.0;
      // 	minDistSftAPro=999.0;meanDistSftAPro=0.0;minDistSftPio=999.0;meanDistSftPio=0.0;
      //   minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;minDistSftIrocPio=999.0;meanDistSftIrocPio=0.0;
      //   minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;minDistSftOrocPio=999.0;meanDistSftOrocPio=0.0;

      // 	// Reset the number of measurements for the mean
      // 	nMeasAPro=0;nMeasPio=0;nMeasAProIroc=0;nMeasPioIroc=0;nMeasAProOroc=0;nMeasPioOroc=0;

      // 	// Check for two-track resolution
      // 	for (iRadius=0;iRadius<9;iRadius++){
      // 	  // Get the spatial distance at each radius
      // 	  distPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fPosDaughter.fXglobal[iRadius],
      // 			     (fFemtoBuffer->GetEvt(iMix))->fBgALamTracks[iBgALam2].fPosDaughter.fXglobal[iRadius]);
      // 	  distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXglobal[iRadius],
      // 			      (fFemtoBuffer->GetEvt(iMix))->fBgALamTracks[iBgALam2].fNegDaughter.fXglobal[iRadius]);
      // 	  // Shifted distances                                                         
      //     distSftPio = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fPosDaughter.fXshifted[iRadius],
      // 				fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fPosDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));
      //     distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXshifted[iRadius],
      // 				 fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam2].fNegDaughter.fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));

      // 	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
      // 	  if (distAPro > -1.0){
      // 	    // Minimum distance
      // 	    if (distAPro < minDistAPro)
      // 	      minDistAPro = distAPro;
      // 	    if (distSftAPro < minDistSftAPro)
      //         minDistSftAPro = distSftAPro;
      // 	    // Mean distance
      // 	    meanDistAPro+=distAPro;
      // 	    meanDistSftAPro+=distSftAPro;
      // 	    nMeasAPro++;
      // 	    // IROC                                                                    
      //       if (iRadius<3){
      //         if (distSftAPro < minDistSftIrocAPro)
      //           minDistSftIrocAPro = distSftAPro;
      //         meanDistSftIrocAPro+=distSftAPro;
      //         nMeasAProIroc++;
      //       }
      //       // OROC                                                                    
      //       else {
      //         if (distSftAPro < minDistSftOrocAPro)
      //           minDistSftOrocAPro = distSftAPro;
      //         meanDistSftOrocAPro+=distSftAPro;
      //         nMeasAProOroc++;
      //       }
      // 	  }

      // 	  if (distPio > -1.0){
      // 	    // Minimum distance
      // 	    if (distPio < minDistPio)
      // 	      minDistPio = distPio;
      // 	    if (distSftPio < minDistSftPio)
      //         minDistSftPio = distSftPio;
      // 	    // Mean distance
      // 	    meanDistPio+=distPio;
      // 	    meanDistSftPio+=distSftPio;
      // 	    nMeasPio++; 
      // 	    // IROC                                                                    
      //       if (iRadius<3){
      //         if (distSftPio < minDistSftIrocPio)
      //           minDistSftIrocPio = distSftPio;
      //         meanDistSftIrocPio+=distSftPio;
      //         nMeasPioIroc++;
      //       }
      //       // OROC                                                                    
      //       else {
      //         if (distSftPio < minDistSftOrocPio)
      //           minDistSftOrocPio = distSftPio;
      //         meanDistSftOrocPio+=distSftPio;
      //         nMeasPioOroc++;
      //       }
      // 	  }
	  
      // 	} // Loop over iRadius
      
      // 	// Require at least one measurement
      // 	if ( (!nMeasPio) || (!nMeasAPro) )
      // 	  continue;
	
      // 	// Divide by the number of measurements to get the mean
      // 	meanDistAPro /= (Float_t)nMeasAPro;
      // 	meanDistPio /= (Float_t)nMeasPio;
      // 	meanDistSftAPro /= (Float_t)nMeasAPro;
      //   meanDistSftPio /= (Float_t)nMeasPio;
	
      // 	// Fill the two track resolution histograms
      // 	f2HistBgALamBgALamMeanMinDistAProMixed->Fill(meanDistAPro,minDistAPro);
      // 	f2HistBgALamBgALamMeanMinDistPioMixed->Fill(meanDistPio,minDistPio);
      // 	f2HistSftBgALamBgALamMeanMinDistAProMixed->Fill(meanDistSftAPro,minDistSftAPro);
      //   f2HistSftBgALamBgALamMeanMinDistPioMixed->Fill(meanDistSftPio,minDistSftPio);

      // 	// Fill IROC / OROC histograms only with at least one measurement              
      //   if (nMeasAProIroc){
      //     meanDistSftIrocAPro /= (Float_t)nMeasAProIroc;
      //     f2HistSftIrocBgALamBgALamMeanMinDistAProMixed->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
      //   }
      //   if (nMeasPioIroc){
      //     meanDistSftPio /= (Float_t)nMeasPioIroc;
      //     f2HistSftIrocBgALamBgALamMeanMinDistPioMixed->Fill(meanDistSftIrocPio,minDistSftIrocPio);
      //   }
      //   if (nMeasAProOroc){
      // 	  meanDistSftAPro /= (Float_t)nMeasAProOroc;
      //     f2HistSftOrocBgALamBgALamMeanMinDistAProMixed->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
      //   }
      //   if (nMeasPioOroc){
      //     meanDistSftPio /= (Float_t)nMeasPioOroc;
      //     f2HistSftOrocBgALamBgALamMeanMinDistPioMixed->Fill(meanDistSftOrocPio,minDistSftOrocPio);
      //   }     	
      // 	//  	// Do a cut (value needs to be refined)
      // 	//  	if ( meanDistSftAPro < 2.0 || meanDistSftPio < 2.0 )
      // 	//  	  continue;
	
      // 	// Fill the qinv, minv histogram
      // 	f3HistBgALamBgALamQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], (fFemtoBuffer->GetEvt(iMix))->fBgALamTracks[iBgALam2]),minDistSftAPro,minDistSftPio);

      // } // Second lambda loop

      // AProton loop
      for (iAPro=0;iAPro<(fFemtoBuffer->GetEvt(iMix))->GetNAPro();iAPro++){
	
  	// Skip if unUseIt() entry
  	if (!(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].UseIt())
  	  continue;
	
	// Reset the distances for each pair
	//	minDistAPro=999.0;meanDistAPro=0.0;
	minDistSftAPro=999.0;meanDistSftAPro=0.0;
        // minDistSftIrocAPro=999.0;meanDistSftIrocAPro=0.0;
        // minDistSftOrocAPro=999.0;meanDistSftOrocAPro=0.0;
	// Reset the number of measurements for the mean
	nMeasAPro=0;//nMeasAProIroc=0;nMeasAProOroc=0;

  	// Check for two-track resolution
  	for (iRadius=0;iRadius<9;iRadius++){
  	  // Get the spatial distance at each radius
	  //  	  distAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXglobal[iRadius],(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].fXglobal[iRadius]);
	  distSftAPro = calcDist(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter.fXshifted[iRadius],(fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro].fXshifted[iRadius]);//,fFemtoBuffer->GetEvt(0),(fFemtoBuffer->GetEvt(iMix)));
	  
	  // calcDist returns -2.0 if one shouldn't use this distance, see calcDist function
	  if (distSftAPro > -1.0){
	    // Minimum distance
	    //	    if (distAPro < minDistAPro)
	    //	      minDistAPro = distAPro;
	    if (distSftAPro < minDistSftAPro)
              minDistSftAPro = distSftAPro;
	    // Mean distance
	    //	    meanDistAPro+=distAPro;
	    meanDistSftAPro+=distSftAPro;
	    nMeasAPro++;
	    // // IROC                                                                    
            // if (iRadius<3){
            //   if (distSftAPro < minDistSftIrocAPro)
            //     minDistSftIrocAPro = distSftAPro;
            //   meanDistSftIrocAPro+=distSftAPro;
            //   nMeasAProIroc++;
            // }
            // // OROC                                                                    
            // else {
            //   if (distSftAPro < minDistSftOrocAPro)
            //     minDistSftOrocAPro = distSftAPro;
            //   meanDistSftOrocAPro+=distSftAPro;
            //   nMeasAProOroc++;
            // }
	  }
	} // Loop over iRadius
	
	// Require at least one measurement
	if ( !nMeasAPro )
	  continue;

	// Divide by the number of measurements to get the mean
	//	meanDistAPro /= (Float_t)nMeasAPro;
	meanDistSftAPro /= (Float_t)nMeasAPro;

  	// // Fill the two track resolution histogram
  	// f2HistBgALamAProMeanMinDistAProMixed->Fill(meanDistAPro,minDistAPro);
	// f2HistSftBgALamAProMeanMinDistAProMixed->Fill(meanDistSftAPro,minDistSftAPro);
  	
	// // Fill IROC / OROC histograms only with at least one measurement              
        // if (nMeasAProIroc){
        //   meanDistSftAPro /= (Float_t)nMeasAProIroc;
        //   f2HistSftIrocBgALamAProMeanMinDistAProMixed->Fill(meanDistSftIrocAPro,minDistSftIrocAPro);
        // }
        // if (nMeasAProOroc){
        //   meanDistSftAPro /= (Float_t)nMeasAProOroc;
        //   f2HistSftOrocBgALamAProMeanMinDistAProMixed->Fill(meanDistSftOrocAPro,minDistSftOrocAPro);
	// }
	//  	// Do a cut (value needs to be refined)
	//	if ( meanDistSftAPro < 2.0 )
	//  	  continue;

	// Use THnSparse since Sept '12
	Double_t x[4]={
	QinvProPro(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam].fNegDaughter,
		   fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro]),
	meanDistSftAPro,
	minDistSftAPro,
	Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], 
	     fFemtoBuffer->GetEvt(iMix)->fAProTracks[iAPro])
      };
      BgALamAProMixed->Fill(x);

  	// Fill the qinv histogram
	//  	f3HistBgALamAProQinvMixed->Fill(Qinv(fFemtoBuffer->GetEvt(0)->fBgALamTracks[iBgALam], (fFemtoBuffer->GetEvt(iMix))->fAProTracks[iAPro]),meanDistSftAPro,minDistSftAPro);
	
	}// AProton loop
    }// BgALambda loop
    
  }// Event buffer loop

}// End of void ProcessMixedBackground

//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Qinv(FemtoBufferV0 v01, FemtoBufferV0 v02){
  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  // Always using lambda mass (no mass difference found yet for lam <-> alam (see PDG))

  // printf("v01 px %3.2f py %3.2f pz %3.2f"
  // 	 "v02 px %3.2f py %3.2f pz %3.2f"
  // 	 "\n"
  // 	 ,v01.fP[0],v01.fP[1],v01.fP[2]
  // 	 ,v02.fP[0],v02.fP[1],v02.fP[2]
  // 	 );

  //Double_t e1 = t1->GetE(mPart1);
  Double_t e1 = sqrt(fkLamMass*fkLamMass + v01.fP[0]*v01.fP[0]+v01.fP[1]*v01.fP[1]+v01.fP[2]*v01.fP[2]);
  //Double_t e2 = t2->GetE(mPart2);
  Double_t e2 = sqrt(fkLamMass*fkLamMass + v02.fP[0]*v02.fP[0]+v02.fP[1]*v02.fP[1]+v02.fP[2]*v02.fP[2]);
  Double_t qinvL;
  Double_t qP;
  Double_t pinv;
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  qinvL = (e1-e2) * (e1-e2) - ( (v01.fP[0]-v02.fP[0])*(v01.fP[0]-v02.fP[0]) + (v01.fP[1]-v02.fP[1])*(v01.fP[1]-v02.fP[1]) + (v01.fP[2]-v02.fP[2])*(v01.fP[2]-v02.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  qP    = (e1-e2)   * (e1+e2)
        - (v01.fP[0]-v02.fP[0]) * (v01.fP[0]+v02.fP[0])
  	- (v01.fP[1]-v02.fP[1]) * (v01.fP[1]+v02.fP[1])
  	- (v01.fP[2]-v02.fP[2]) * (v01.fP[2]+v02.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  pinv  = (e1+e2) * (e1+e2) - ( (v01.fP[0]+v02.fP[0])*(v01.fP[0]+v02.fP[0])
  			       +(v01.fP[1]+v02.fP[1])*(v01.fP[1]+v02.fP[1])
  			       +(v01.fP[2]+v02.fP[2])*(v01.fP[2]+v02.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Qinv(FemtoBufferV0 v0, FemtoBufferTrack track) {
  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  //  Always using lambda mass (no mass difference found yet for lam <-> alam (see PDG))
  
  //  Double_t e1 = t1->GetE(mPart1);
  Double_t e1 = sqrt(fkLamMass*fkLamMass + v0.fP[0]*v0.fP[0]+v0.fP[1]*v0.fP[1]+v0.fP[2]*v0.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  Double_t e2 = sqrt(fkProMass*fkProMass + track.fP[0]*track.fP[0]+track.fP[1]*track.fP[1]+track.fP[2]*track.fP[2]);
  Double_t qinvL;
  Double_t qP;
  Double_t pinv;
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  qinvL = (e1-e2) * (e1-e2) - ( (v0.fP[0]-track.fP[0])*(v0.fP[0]-track.fP[0]) + (v0.fP[1]-track.fP[1])*(v0.fP[1]-track.fP[1]) + (v0.fP[2]-track.fP[2])*(v0.fP[2]-track.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  qP    = (e1-e2)   * (e1+e2)
        - (v0.fP[0]-track.fP[0]) * (v0.fP[0]+track.fP[0])
  	- (v0.fP[1]-track.fP[1]) * (v0.fP[1]+track.fP[1])
  	- (v0.fP[2]-track.fP[2]) * (v0.fP[2]+track.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  pinv  = (e1+e2) * (e1+e2) - ( (v0.fP[0]+track.fP[0])*(v0.fP[0]+track.fP[0])
  			       +(v0.fP[1]+track.fP[1])*(v0.fP[1]+track.fP[1])
  			       +(v0.fP[2]+track.fP[2])*(v0.fP[2]+track.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Qinv(FemtoBufferTrack track, FemtoBufferV0 v0){
  return Qinv(v0, track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::QinvProPro(FemtoBufferTrack proTrack1, FemtoBufferTrack proTrack2) {
  // Same as above, with different masses for the tracks,
  // here both tracks are protons

  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  
  //  PDG_t e1 = t1->GetE(mPart1);
  Double_t e1 = sqrt(fkProMass*fkProMass + proTrack1.fP[0]*proTrack1.fP[0]+proTrack1.fP[1]*proTrack1.fP[1]+proTrack1.fP[2]*proTrack1.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  Double_t e2 = sqrt(fkProMass*fkProMass + proTrack2.fP[0]*proTrack2.fP[0]+proTrack2.fP[1]*proTrack2.fP[1]+proTrack2.fP[2]*proTrack2.fP[2]);
  Double_t qinvL;
  Double_t qP;
  Double_t pinv;
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  qinvL = (e1-e2) * (e1-e2) - ( (proTrack1.fP[0]-proTrack2.fP[0])*(proTrack1.fP[0]-proTrack2.fP[0]) + (proTrack1.fP[1]-proTrack2.fP[1])*(proTrack1.fP[1]-proTrack2.fP[1]) + (proTrack1.fP[2]-proTrack2.fP[2])*(proTrack1.fP[2]-proTrack2.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  qP    = (e1-e2)   * (e1+e2)
        - (proTrack1.fP[0]-proTrack2.fP[0]) * (proTrack1.fP[0]+proTrack2.fP[0])
  	- (proTrack1.fP[1]-proTrack2.fP[1]) * (proTrack1.fP[1]+proTrack2.fP[1])
  	- (proTrack1.fP[2]-proTrack2.fP[2]) * (proTrack1.fP[2]+proTrack2.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  pinv  = (e1+e2) * (e1+e2) - ( (proTrack1.fP[0]+proTrack2.fP[0])*(proTrack1.fP[0]+proTrack2.fP[0])
  			       +(proTrack1.fP[1]+proTrack2.fP[1])*(proTrack1.fP[1]+proTrack2.fP[1])
  			       +(proTrack1.fP[2]+proTrack2.fP[2])*(proTrack1.fP[2]+proTrack2.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::QinvPioPro(FemtoBufferTrack pioTrack, FemtoBufferTrack proTrack) {
  // Same as above, with different masses for the tracks,
  // here both tracks are protons

  // Copied from NA49. See http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Qinv
  
  //  PDG_t e1 = t1->GetE(mPart1);
  Double_t e1 = sqrt(fkPioMass*fkPioMass + pioTrack.fP[0]*pioTrack.fP[0]+pioTrack.fP[1]*pioTrack.fP[1]+pioTrack.fP[2]*pioTrack.fP[2]);
  //  Double_t e2 = t2->GetE(mPart2);
  Double_t e2 = sqrt(fkProMass*fkProMass + proTrack.fP[0]*proTrack.fP[0]+proTrack.fP[1]*proTrack.fP[1]+proTrack.fP[2]*proTrack.fP[2]);
  Double_t qinvL;
  Double_t qP;
  Double_t pinv;
  
  // First calculate -Qinv^2  as usual : 
  //qinvL = (e1-e2) * (e1-e2) - Q(t1,t2) * Q(t1,t2);
  qinvL = (e1-e2) * (e1-e2) - ( (pioTrack.fP[0]-proTrack.fP[0])*(pioTrack.fP[0]-proTrack.fP[0]) + (pioTrack.fP[1]-proTrack.fP[1])*(pioTrack.fP[1]-proTrack.fP[1]) + (pioTrack.fP[2]-proTrack.fP[2])*(pioTrack.fP[2]-proTrack.fP[2]) );  

  //Qx(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()-t2->GetPx()); };
  //static Double_t Qy(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPy()-t2->GetPy()); };
  //static Double_t Qz(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPz()-t2->GetPz()); };
  //static Double_t  Q(T49ParticleRoot* t1,T49ParticleRoot* t2)
  //               { return  TMath::Sqrt(Qx(t1,t2)*Qx(t1,t2)+Qy(t1,t2)*Qy(t1,t2)+Qz(t1,t2)*Qz(t1,t2)); };
  


  //static Double_t Px(T49ParticleRoot* t1,T49ParticleRoot* t2) { return (t1->GetPx()+t2->GetPx()); };
  //qP    = (e1-e2)   * (e1+e2)
  //      - Qx(t1,t2) * Px(t1,t2)
  //      - Qy(t1,t2) * Py(t1,t2)
  //      - Qz(t1,t2) * Pz(t1,t2);
  qP    = (e1-e2)   * (e1+e2)
        - (pioTrack.fP[0]-proTrack.fP[0]) * (pioTrack.fP[0]+proTrack.fP[0])
  	- (pioTrack.fP[1]-proTrack.fP[1]) * (pioTrack.fP[1]+proTrack.fP[1])
  	- (pioTrack.fP[2]-proTrack.fP[2]) * (pioTrack.fP[2]+proTrack.fP[2]);

  //pinv  = (e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2);
  pinv  = (e1+e2) * (e1+e2) - ( (pioTrack.fP[0]+proTrack.fP[0])*(pioTrack.fP[0]+proTrack.fP[0])
  			       +(pioTrack.fP[1]+proTrack.fP[1])*(pioTrack.fP[1]+proTrack.fP[1])
  			       +(pioTrack.fP[2]+proTrack.fP[2])*(pioTrack.fP[2]+proTrack.fP[2]));

  return TMath::Sqrt(qP*qP/pinv - qinvL);
}
//________________________________________________________________________
// Float_t AliAnalysisTaskProtonLambda::QinvConstr(FemtoBufferV0 v0, FemtoBufferTrack track) {
//   // Same as Qinv(v0,track) but with constrained momentum for the track

//   // Check whether constrained momentum is there
//   if ((track.fPconstr[0]<0.00001)&&(track.fPconstr[1]<0.00001)&&(track.fPconstr[2]<0.00001))
//     return Qinv(v0,track);

//   // Standard Qinv(v0, track), just with constrained momentum instead of TPC only momentum
//   Double_t e1 = sqrt(fkLamMass*fkLamMass + v0.fP[0]*v0.fP[0]+v0.fP[1]*v0.fP[1]+v0.fP[2]*v0.fP[2]);
//   Double_t e2 = sqrt(fkProMass*fkProMass + track.fPconstr[0]*track.fPconstr[0]+track.fPconstr[1]*track.fPconstr[1]+track.fPconstr[2]*track.fPconstr[2]);
//   Double_t qinvL;
//   Double_t qP;
//   Double_t pinv;
//   qinvL = (e1-e2) * (e1-e2) - ( (v0.fP[0]-track.fPconstr[0])*(v0.fP[0]-track.fPconstr[0]) + (v0.fP[1]-track.fPconstr[1])*(v0.fP[1]-track.fPconstr[1]) + (v0.fP[2]-track.fPconstr[2])*(v0.fP[2]-track.fPconstr[2]) );  
//   qP    = (e1-e2)   * (e1+e2)
//         - (v0.fP[0]-track.fPconstr[0]) * (v0.fP[0]+track.fPconstr[0])
//   	- (v0.fP[1]-track.fPconstr[1]) * (v0.fP[1]+track.fPconstr[1])
//   	- (v0.fP[2]-track.fPconstr[2]) * (v0.fP[2]+track.fPconstr[2]);
//   pinv  = (e1+e2) * (e1+e2) - ( (v0.fP[0]+track.fPconstr[0])*(v0.fP[0]+track.fPconstr[0])
//   			       +(v0.fP[1]+track.fPconstr[1])*(v0.fP[1]+track.fPconstr[1])
//   			       +(v0.fP[2]+track.fPconstr[2])*(v0.fP[2]+track.fPconstr[2]));

//   return TMath::Sqrt(qP*qP/pinv - qinvL);
// }
// //________________________________________________________________________
// Float_t AliAnalysisTaskProtonLambda::QinvConstr(FemtoBufferTrack track, FemtoBufferV0 v0){
//   return QinvConstr(v0, track);
// }
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Minv(FemtoBufferV0 v01, FemtoBufferV0 v02){
  // Taken from NA49. See 
  // http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Minv
  
  //  Double_t e1 = t1->GetE(mPart1);
  //  Double_t e2 = t2->GetE(mPart2);
  //  GetE(Float_t mass)  { return sqrt(GetP()*GetP()+mass*mass); }  
  Float_t e1 = TMath::Sqrt(v01.fP[0]*v01.fP[0] + v01.fP[1]*v01.fP[1] + v01.fP[2]*v01.fP[2] 
  			   + fkLamMass*fkLamMass);
  Float_t e2 = TMath::Sqrt(v02.fP[0]*v02.fP[0] + v02.fP[1]*v02.fP[1] + v02.fP[2]*v02.fP[2] 
  			   + fkLamMass*fkLamMass);
  
  // return TMath::Sqrt((e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2));
  return TMath::Sqrt((e1+e2) * (e1+e2) - (  (v01.fP[0]+v02.fP[0])*(v01.fP[0]+v02.fP[0])
  					   +(v01.fP[1]+v02.fP[1])*(v01.fP[1]+v02.fP[1])
  					   +(v01.fP[2]+v02.fP[2])*(v01.fP[2]+v02.fP[2])));

}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Minv(FemtoBufferV0 v0, FemtoBufferTrack track){
  // Taken from NA49. See 
  // http://na49info.web.cern.ch/na49info/na49/Software/minidst/ana/html/src/T49Tool.cxx.html#T49Tool:Minv
  
  //  Double_t e1 = t1->GetE(mPart1);
  //  Double_t e2 = t2->GetE(mPart2);
  //  GetE(Float_t mass)  { return sqrt(GetP()*GetP()+mass*mass); }  
  Float_t e1 = TMath::Sqrt(v0.fP[0]*v0.fP[0] + v0.fP[1]*v0.fP[1] + v0.fP[2]*v0.fP[2] 
  			   + fkLamMass*fkLamMass);
  Float_t e2 = TMath::Sqrt(track.fP[0]*track.fP[0] + track.fP[1]*track.fP[1] + track.fP[2]*track.fP[2] 
  			   + fkProMass*fkProMass);
  
  // return TMath::Sqrt((e1+e2) * (e1+e2) - P(t1,t2) * P(t1,t2));
  return TMath::Sqrt((e1+e2) * (e1+e2) - (  (v0.fP[0]+track.fP[0])*(v0.fP[0]+track.fP[0])
  					   +(v0.fP[1]+track.fP[1])*(v0.fP[1]+track.fP[1])
  					   +(v0.fP[2]+track.fP[2])*(v0.fP[2]+track.fP[2])));

}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::Minv(FemtoBufferTrack track, FemtoBufferV0 v0){
    return Minv(v0, track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::calcDist(const Float_t r1[3], const Float_t r2[3]){
  // Return the spatial distance of two space vectors r1 and r2
  
  // Return 'error' when no position is given.
  // // When a particle doesn't manage to get to a certain radius,
  // // the function GetXYZAt returns position 0,0,0
  // if ( (r1[0] < 0.0001 && r1[1] < 0.0001 && r1[1] < 0.0001) ||
  //      (r2[0] < 0.0001 && r2[1] < 0.0001 && r2[1] < 0.0001) )
  //   return -2.0;
  
  // The above stuff is stupid, this discards every 
  // track just going in the negative direction.
  // Also we don't use the GetXYZAt anymore and our
  // 'bad position' value is -9999.,-9999.-9999.
  if ( (r1[0] < -9998. && r1[1] < -9998. && r1[1] < -9998.) ||
       (r2[0] < -9998. && r2[1] < -9998. && r2[1] < -9998.) )
    return -2.0;
  

  return TMath::Sqrt((r1[0]-r2[0])*(r1[0]-r2[0]) +
   		     (r1[1]-r2[1])*(r1[1]-r2[1]) +
   		     (r1[2]-r2[2])*(r1[2]-r2[2]) );
}
//________________________________________________________________________
// This function is no longer needed
// Float_t AliAnalysisTaskProtonLambda::calcDistSft(const Float_t r1[3], const Float_t r2[3], const FemtoBufferEvent *evt1, const FemtoBufferEvent *evt2){
//   // Return the spatial distance of two space vectors r1 and r2
//   // With each event shifted to (0,0,0)
  
//   // // Return 'error' when no position is given.
//   // // When a particle doesn't manage to get to a certain radius,
//   // // the function GetXYZAt returns position 0,0,0
//   // if ( (r1[0] < 0.0001 && r1[1] < 0.0001 && r1[2] < 0.0001) ||
//   //      (r2[0] < 0.0001 && r2[1] < 0.0001 && r2[2] < 0.0001) )
//   //   return -2.0;
//   // The above stuff is stupid, this discards every 
//   // track just going in the negative direction.
//   // Also we don't use the GetXYZAt anymore and our
//   // 'bad position' value is -9999.,-9999.-9999.
//   if ( (r1[0] < -9998. && r1[1] < -9998. && r1[1] < -9998.) ||
//        (r2[0] < -9998. && r2[1] < -9998. && r2[1] < -9998.) )
//     return -2.0;

//   // Get the vertex postions 
//   Double_t vtx1[3],vtx2[3];
//   evt1->GetVtxPos(vtx1);
//   evt2->GetVtxPos(vtx2);
  
//   // Calculate shifted positions
//   Double_t r1Sft[3],r2Sft[3];
//   for (Int_t i=0;i<3;i++){
//     r1Sft[i]=r1[i] - vtx1[i];
//     r2Sft[i]=r2[i] - vtx2[i];
//   }

//   // Return shifted distances
//   return TMath::Sqrt((r1Sft[0]-r2Sft[0])*(r1Sft[0]-r2Sft[0]) +
//    		     (r1Sft[1]-r2Sft[1])*(r1Sft[1]-r2Sft[1]) +
//    		     (r1Sft[2]-r2Sft[2])*(r1Sft[2]-r2Sft[2]) );
// }
//________________________________________________________________________
// void AliAnalysisTaskProtonLambda::constrainTrack(AliAODTrack *track) {
//   // Abuses data members of the AliAODTrack to store a set of track
//   // parameters for TPC only constrained to the primary vtx
//   // plus a bool whether constraining was successful
//   if (!track->GetConstrainedParam())
//     return;
  
//   // Constrain track to pri vtx, set the bool to successful / unsuccessful
//   track->SetTOFcluster(track->RelateToVertexTPC(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(),
// 						5.0, const_cast<AliExternalTrackParam *> (track->GetConstrainedParam())));
// }
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda::goodDCA(AliAODTrack *track) {
  // Get the DCAxy and DCAz. There also exists a TPC only 
  // impact parameter, but this has not enough resolution 
  // to discriminate between primaries, secondaries and material
  Float_t xy=0.,rap=RapidityProton(track),pt=track->Pt();
  xy = DCAxy(fGTI[-track->GetID()-1], fAOD);
  // Fill the DCAxy histograms
  if (track->Charge() > 0){
    fPriHistDCAxyYPtPro->Fill(xy,rap,pt);
  }
  else{
    fPriHistDCAxyYPtAPro->Fill(xy,rap,pt);
  }
  // Do a cut. 0.1 cm shows highest significance for primaries
  if (xy>.1)
    return kFALSE;
  return kTRUE;
}
//_______________________________________________________________
Float_t AliAnalysisTaskProtonLambda::RapidityProton(AliAODTrack *track){
  // Can't find how to set the assumed mass for the AliAODTrack.
  // Same stuff as in AliAODTrack::Y() just with proton mass
  Double_t e = TMath::Sqrt(track->P()*track->P() + fkProMass*fkProMass);
  Double_t pz = track->Pz();
  if (e != TMath::Abs(pz)) { // energy was not equal to pz
    return 0.5*TMath::Log((e+pz)/(e-pz));
  } else { // energy was equal to pz
    return -999.;
  }
}
//________________________________________________________________________
// void AliAnalysisTaskProtonLambda::getTPConlyV0Info(const AliAODTrack *posDaughter,const AliAODTrack *negDaughter, Double_t tpcV0Mom[3], Double_t TPConlyV0MinvLam, Double_t TPConlyV0MinvALam){
//   //
//   // Calculates a V0 with the TPC only parameters
//   //

//   // Duplicate the V0
//   AliV0 tpcV0;
//   // Get the TPC only track parameters from the daughters
//   const AliExternalTrackParam *pParam = 0, *nParam = 0;
//   pParam = posDaughter->GetTPCInnerParam();
//   nParam = negDaughter->GetTPCInnerParam();
//   // Protection if there's no TPC only track parameters
//   if(!pParam||!nParam)
//     return;
//   // Set the tpcV0 daughters to the TPC only ones
//   tpcV0.SetParamP(*pParam);
//   tpcV0.SetParamN(*nParam);
//   // Calculate the new properties of the V0
//   Double_t vertex[3];
//   fFemtoBuffer->GetEvt(0)->GetVtxPos(vertex);
//   tpcV0.Update((Float_t *) vertex);
//   // Get the updated momentum
//   tpcV0.fPxPyPz(tpcV0Mom);
//   // New TPC only mass, lambda..
//   tpcV0.ChangeMassHypothesis(3122);
//   TPConlyV0MinvLam = tpcV0.GetEffMass();
//   // ..anti-lambda.
//   tpcV0.ChangeMassHypothesis(-3122);
//   TPConlyV0MinvALam = tpcV0.GetEffMass();
// }
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::DCAxy(const AliAODTrack *track, const AliVEvent *evt){
  // Note that AliAODTrack::PropagateToDCA() changes the track. 
  // Don't know whether this is what one wants?
  if(!track){
    printf("Pointer to track is zero!\n");
    return -9999.;
  }

  // Create an external parameter from the AODtrack
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  // Propagation through the beam pipe would need a correction 
  // for material, I guess.
  if(etp.GetX()>3.) {
    printf("This method can be used only for propagation inside the beam pipe\n");
    printf("  id: %d, filtermap: %d\n",track->GetID(),track->GetFilterMap());
    return -9999.; 
  }
  // Do the propagation
  Double_t dca[2]={-9999.,-9999.},covar[3]={0.,0.,0.};
  if(!etp.PropagateToDCA(evt->GetPrimaryVertex(),evt->GetMagneticField(),10.,dca,covar)) return -9999.;
  // return the DCAxy
  return dca[0];
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FillDedxHist(const AliVTrack *track){
  // This is for visualization. Fill the the dE/dx histograms
  // for all tracks, not only for those, where only the TPC
  // is used for PID. Thus avoiding the sharp cut off at a 
  // momentum of 0.75 GeV/c.

  if(!(fGTI[-track->GetID()-1])){
    printf("Warning: No global track info there!\n");
    return;
  }

  // TPC signal and Nsigma. See STEER/STEERBase/AliPIDResponse.h for how 
  // NSigmaTPC works (and refrain from banging your head against the wall
  // when you see it).
  // Positive tracks
  if (track->Charge() > 0){
    fPriHistTPCsignalPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    // fPriHistNsigmaTPCPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
    // 			       fPIDResponse->GetTPCResponse().GetNumberOfSigmas((fGTI[-track->GetID()-1])->GetTPCmomentum()
    // 										,(fGTI[-track->GetID()-1])->GetTPCsignal()
    // 										,(fGTI[-track->GetID()-1])->GetTPCsignalN()
    // 										,AliPID::kProton));
    // Fill histograms in three momentum ranges
    fPriHistTPCsignalLowPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalMedPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalHigPPos->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());  
    
  }
  // Negative tracks
  else{ 
    fPriHistTPCsignalNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    // fPriHistNsigmaTPCNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),
    // 			       fPIDResponse->GetTPCResponse().GetNumberOfSigmas((fGTI[-track->GetID()-1])->GetTPCmomentum()
    // 										,(fGTI[-track->GetID()-1])->GetTPCsignal()
    // 										,(fGTI[-track->GetID()-1])->GetTPCsignalN()
    // 										,AliPID::kProton));
    // Fill histograms in three momentum ranges
    fPriHistTPCsignalLowPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalMedPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());
    fPriHistTPCsignalHigPNeg->Fill((fGTI[-track->GetID()-1])->GetTPCmomentum(),(fGTI[-track->GetID()-1])->GetTPCsignal());  
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::StoreGlobalTrackReference(AliAODTrack *track){
  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see they don't have
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no 
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  // 	   //	   "not adding it's info! "
  // 	   "ID: %d FilterMap: %d\n"
  // 	   ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }
  
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
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::ResetGlobalTrackReference(){
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for(UShort_t i=0;i<fTrackBuffSize;i++){
    fGTI[i]=0;
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda::acceptTrack(const AliAODTrack *track){
  // Apply additional track cuts

  // In the documents
  // https://alisoft.cern.ch/AliRoot/trunk/TPC/doc/Definitions/Definitions.pdf
  // TPC people describe the cut strategy for the TPC. It is explicitly
  // stated that a cut on the number of crossed rows and a cut on the
  // number of crossed rows over findable clusters is recommended to 
  // remove fakes. In the pdf a cut value of .83 on the ratio 
  // is stated, no value for the number of crossed rows. Looking at the 
  // AliESDtrackCuts.cxx one sees that exactly this cut is used with
  // 0.8 on the ratio and 70 on the crossed rows.

  // Checked the filter task and AliAODTrack and AliESDtrack and
  // AliESDtrackCuts and the Definitions.pdf:
  // The function to get the findable clusters is GetTPCNclsF()
  
  // For the number fo crossed rows for ESD tracks, the function
  // GetTPCCrossedRows() usually is used. Looking at the AliESDtrack.cxx
  // one sees that it's just an alias (with additional caching) for
  // GetTPCClusterInfo(2, 1); The identical function exists in the
  // AliAODTrack.cxx

  // I checked: for AOD073 both, the number of crossed rows and the
  // number of findable clusters, are there.

  // WARNING: in LHC10h pass2 the cluster map is wrong for 
  // sector 0 / 18. It's used in the calculation of
  // the number of crossed rows!

  Float_t nCrossed = track->GetTPCClusterInfo(2, 1);
  if(nCrossed<70)
    return kFALSE;
  if(!track->GetTPCNclsF())
    return kFALSE; // Note that the AliESDtrackCuts would here return kTRUE
  if((nCrossed/track->GetTPCNclsF()) < .8)
    return kFALSE;
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda::GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,
							   const AliAODTrack *nTrack){
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for positive and negative daughters from V0s

  // Get the shared maps
  const TBits posSharedMap = pTrack->GetTPCSharedMap();
  const TBits negSharedMap = nTrack->GetTPCSharedMap();
  // Fill a control histogram
  //  fHistShareV0pos->Fill(posSharedMap.CountBits());
  //  fHistShareV0neg->Fill(negSharedMap.CountBits());
  // Reject shared clusters
  if( ((posSharedMap.CountBits()) >= 1) ||
      ((negSharedMap.CountBits()) >= 1)){
    // Bad tracks, have too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonLambda::GoodTPCFitMapSharedMap(const AliAODTrack *track){
  // Rejects tracks with shared clusters after filling a control histogram
  // This overload is used for primaries

  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  // Fill a control histogram
  fPriHistShare->Fill(sharedMap.CountBits());
  // Reject shared clusters
  if((sharedMap.CountBits()) >= 1){
    // Bad track, has too many shared clusters!
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::GetCorrectedTOFSignal(const AliVTrack *track){
  // Return the corrected TOF signal, see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

  // Check for the global track
  if(!(fGTI[-track->GetID()-1])){
    printf("Warning: no corresponding global track found!\n");
    return -9999.;
  }

  // Request the TOFpid bit
  if(!((fGTI[-track->GetID()-1])->GetStatus() & AliVTrack::kTOFpid))
    return -9999.;

  // The expected time
  Double_t expectedTimes[AliPID::kSPECIES];
  (fGTI[-track->GetID()-1])->GetIntegratedTimes(expectedTimes);

  // Check for TOF header
  if(fAOD->GetTOFHeader()){
    // New AODs without start time subtraction
    return ((fGTI[-track->GetID()-1])->GetTOFsignal()
	    - expectedTimes[AliPID::kProton]
	    - fPIDResponse->GetTOFResponse().GetStartTime(track->P()));
  }

  // Old AODs with start time already subtracted
  return ((fGTI[-track->GetID()-1])->GetTOFsignal()
	   - expectedTimes[AliPID::kProton]);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::mt(FemtoBufferTrack track, FemtoBufferV0 v0) {
  // Overloaded function
  return mt(v0,track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::mt(FemtoBufferV0 v0, FemtoBufferTrack track){
  // Returns the transverse mass of the pair assuming 
  // proton mass for track and lambda mass for v0

  // Following Phys Rev C 83, 054906
  return TMath::Sqrt(ktSquared(v0,track) +
		     TMath::Power((0.5*(fkLamMass + fkProMass)),2));
}

//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::mt(FemtoBufferV0 v01, FemtoBufferV0 v02){
  // Returns the transverse mass of the pair assuming 
  // lambda mass for both v0

  // Following Phys Rev C 83, 054906
  return TMath::Sqrt(ktSquared(v01,v02) +
		     TMath::Power(fkLamMass,2));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::ktSquared(FemtoBufferV0 v01, FemtoBufferV0 v02){
  // Returns the kt squared
  // kt = 1/2 * | (vector{pt1} + vector{pt2}) |
  // kt = 1/2 * | ({px1+px2}, {py1+py2}) |
  // kt2 = 1/2*1/2 * ( (px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) )
  return .5*.5*(  (v01.fP[0] + v02.fP[0])*(v01.fP[0] + v02.fP[0])
		+ (v01.fP[1] + v02.fP[1])*(v01.fP[1] + v02.fP[1]));
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::ktSquared(FemtoBufferTrack track, FemtoBufferV0 v0){
  // Overloaded function
  return ktSquared(v0,track);
}
//________________________________________________________________________
Float_t AliAnalysisTaskProtonLambda::ktSquared(FemtoBufferV0 v0, FemtoBufferTrack track){
  // Returns the kt squared
  // kt = 1/2 * | (vector{pt1} + vector{pt2}) |
  // kt = 1/2 * | ({px1+px2}, {py1+py2}) |
  // kt2 = 1/2*1/2 * ( (px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) )
  return .5*.5*(  (v0.fP[0] + track.fP[0])*(v0.fP[0] + track.fP[0])
		+ (v0.fP[1] + track.fP[1])*(v0.fP[1] + track.fP[1]));
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::AliAnalysisTaskProtonLambda(const AliAnalysisTaskProtonLambda& atpl)
  // Not implemented, only initializing the const data member as the compiler complains.
  // Implementation is straight forward, though.
  : AliAnalysisTaskSE(atpl),
    fkUseOnTheFly(atpl.fkUseOnTheFly),
    fkAbsZvertexCut(atpl.fkAbsZvertexCut),
    fkCentCut(atpl.fkCentCut),
    fkLamMass(atpl.fkLamMass),
    fkProMass(atpl.fkProMass),
    fkPioMass(atpl.fkPioMass),
    
    fPIDResponse(0), 
    fTpcResponse(0),
    fFemtoBuffer(0),
    fAOD(0), fPrimaryVtx(0), fOutputList(0), fOutputPrimaries(0),
    fOutput2Part(0),
    fGTI(0),    
    fTrackBuffSize(atpl.fTrackBuffSize),
    fHistGoodEvent(0),
    // fHistPrimaryVertexPosXY(0), fHistPrimaryVertexPosZ(0),        
    // fHistTrackMultiplicity(0),    
    // fHistShareV0pos(0),fHistShareV0neg(0),
    // fHistPosTofBeforeCut(0), fHistPosTofAfterCut(0),           
    // fHistNegTofBeforeCut(0), fHistNegTofAfterCut(0),           
    // fHistPosTpcBeforeCut(0), fHistPosTpcAfterCut(0),            
    // fHistNegTpcBeforeCut(0), fHistNegTpcAfterCut(0),            
    // fHistGoodV0(0), fHistCorrectSigns(0),              
    // fHistDcaPosToPrimVertex(0), fHistDcaNegToPrimVertex(0),        
    // fHistDcaPosToPrimVertexZoom(0), fHistDcaNegToPrimVertexZoom(0),  
    // fHistRadiusV0(0), fHistDecayLengthV0(0), fHistDcaV0Daughters(0),          
    // fHistChi2(0), fHistCosPointAngle(0), fHistCosPointAngleZoom(0),
    fHistSideBandOffLam(0), fHistSideBandOffALam(0), fHistTPCNclsPosOffLam(0),       
    fHistTPCNclsNegOffLam(0), fHistTPCNclsPosOffALam(0), fHistTPCNclsNegOffALam(0),      
    // fHistPosNsigmaTpcOffLam(0), fHistPosNsigmaTpcOffALam(0), fHistNegNsigmaTpcOffLam(0),
    // fHistNegNsigmaTpcOffALam(0), fHistUseTofOffLam(0), fHistUseTofOffALam(0),
    // fHistDcaPosOffLam(0), fHistDcaPosOffALam(0), fHistDcaNegOffLam(0),           
    // fHistDcaNegOffALam(0), fHistDcaV0DaughtersOffLam(0), fHistDcaV0DaughtersOffALam(0),  
    // fHistCosPointLamOff(0), fHistCosPointALamOff(0), fHistCosPointLamZoomOff(0),     
    // fHistCosPointALamZoomOff(0), fHistV0RadiusLamOff(0), fHistV0RadiusALamOff(0),        
    // fHistV0DecayLengthLamOff(0), fHistV0DecayLengthALamOff(0), fHistDcaV0PriVertexLamOff(0),     
    // fHistDcaV0PriVertexALamOff(0),
    fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),        
    // fHistPtVsMassLambdaOff(0), fHistPtVsMassAntiLambdaOff(0),
    fHistYPtMassLamOff(0), fHistYPtMassALamOff(0),
    // fHistPtVsYLambdaOff(0), fHistPtVsYAntiLambdaOff(0),       
    fHistSideBandOnLam(0), fHistSideBandOnALam(0),
    // fHistLikeSignOnLam(0), fHistLikeSignOnALam(0),         
    fHistTPCNclsPosOnLam(0), fHistTPCNclsNegOnLam(0), fHistTPCNclsPosOnALam(0),fHistTPCNclsNegOnALam(0),     
    // fHistPosNsigmaTpcOnLam(0), fHistPosNsigmaTpcOnALam(0), fHistNegNsigmaTpcOnLam(0), fHistNegNsigmaTpcOnALam(0),        
    // fHistUseTofOnLam(0),fHistUseTofOnALam(0),fHistDcaPosOnLam(0),fHistDcaPosOnALam(0),fHistDcaNegOnLam(0),               
    // fHistDcaNegOnALam(0),fHistDcaV0DaughtersOnLam(0),fHistDcaV0DaughtersOnALam(0),fHistCosPointLamOn(0),             
    // fHistCosPointALamOn(0),fHistCosPointLamZoomOn(0),fHistCosPointALamZoomOn(0),fHistV0RadiusLamOn(0),             
    // fHistV0RadiusALamOn(0),fHistV0DecayLengthLamOn(0),fHistV0DecayLengthALamOn(0),fHistDcaV0PriVertexLamOn(0),       
    // fHistDcaV0PriVertexALamOn(0),
    // fHistChi2TPCPosLamOn(0),  fHistChi2TPCPosALamOn(0),  fHistChi2TPCNegLamOn(0),  fHistChi2TPCNegALamOn(0),
    // fHistMinvTPConlyLamOn(0),  fHistMinvTPConlyALamOn(0),
    fHistMassLambdaOn(0),fHistMassAntiLambdaOn(0),
    // fHistPtVsMassLambdaOn(0), fHistPtVsMassAntiLambdaOn(0),
    fHistYPtMassLamOn(0),fHistYPtMassALamOn(0),
    // fHistPtVsYLambdaOn(0), fHistPtVsYAntiLambdaOn(0),
    // fHistMomDiffLam(0),fHistMomDiffALam(0),fHistMomDiffBgLam(0),fHistMomDiffBgALam(0),
    // fHistMomDiffWoSPDLam(0),fHistMomDiffWoSPDALam(0),fHistMomDiffWoSPDBgLam(0),fHistMomDiffWoSPDBgALam(0),
    fPriHistShare(0),
    // fPriHistPosNsigmaTof(0),
    fPriHistPosNsigmaTofVsP(0),fPriHistPosNsigmaTofVsPt(0),     
    // fPriHistNegNsigmaTof(0),
    fPriHistNegNsigmaTofVsP(0),fPriHistNegNsigmaTofVsPt(0),fPriHistTOFsignalPosVsP(0),      
    fPriHistTOFsignalPosVsPt(0),fPriHistTOFsignalNegVsP(0),fPriHistTOFsignalNegVsPt(0),fPriHistHybridTOFsigPosWoTPC(0), 
    fPriHistHybridTOFsigPosTPCok(0),fPriHistHybridTOFsigNegWoTPC(0),fPriHistHybridTOFsigNegTPCok(0),
    // fPriHistHasTofPos(0),          
    fPriHistTPCsignalPos(0),
    // fPriHistNsigmaTPCPos(0), fPriHistTPCsignalTOFcutPos(0),fPriHistNsigmaTPCTOFcutPos(0),   
    fPriHistTPCsignalLowPPos(0),fPriHistTPCsignalMedPPos(0),fPriHistTPCsignalHigPPos(0),
    // fPriHistHasTofNeg(0),         
    fPriHistTPCsignalNeg(0),
    // fPriHistNsigmaTPCNeg(0),fPriHistTPCsignalTOFcutNeg(0),fPriHistNsigmaTPCTOFcutNeg(0),    
    fPriHistTPCsignalLowPNeg(0),fPriHistTPCsignalMedPNeg(0),fPriHistTPCsignalHigPNeg(0),
    fPriHistDCAxyYPtPro(0),fPriHistDCAxyYPtAPro(0),
    // f2HistLamLamMeanMinDistProReal(0),
    // f2HistLamLamMeanMinDistPioReal(0),f2HistLamProMeanMinDistProReal(0),f2HistALamALamMeanMinDistAProReal(0), 
    // f2HistALamALamMeanMinDistPioReal(0),f2HistALamAProMeanMinDistAProReal(0),
    // f2HistSftLamLamMeanMinDistProReal(0),
    // f2HistSftLamLamMeanMinDistPioReal(0),f2HistSftLamProMeanMinDistProReal(0),f2HistSftALamALamMeanMinDistAProReal(0), 
    // f2HistSftALamALamMeanMinDistPioReal(0),f2HistSftALamAProMeanMinDistAProReal(0),
    // f2HistSftIrocLamLamMeanMinDistProReal(0),
    // f2HistSftIrocLamLamMeanMinDistPioReal(0),f2HistSftIrocLamProMeanMinDistProReal(0),f2HistSftIrocALamALamMeanMinDistAProReal(0), 
    // f2HistSftIrocALamALamMeanMinDistPioReal(0),f2HistSftIrocALamAProMeanMinDistAProReal(0),
    // f2HistSftOrocLamLamMeanMinDistProReal(0),
    // f2HistSftOrocLamLamMeanMinDistPioReal(0),f2HistSftOrocLamProMeanMinDistProReal(0),f2HistSftOrocALamALamMeanMinDistAProReal(0), 
    // f2HistSftOrocALamALamMeanMinDistPioReal(0),f2HistSftOrocALamAProMeanMinDistAProReal(0),
    // f2HistMtLamLamReal(0), 
    f2HistMtLamProReal(0), 
    // f2HistMtALamALamReal(0), 
    f2HistMtALamAProReal(0),
    // f2HistMtLowQLamLamReal(0), 
    f2HistMtLowQLamProReal(0), 
    // f2HistMtLowQALamALamReal(0), 
    f2HistMtLowQALamAProReal(0),
    LamProReal(0),ALamAProReal(0),
    // f3HistLamLamQinvReal(0),               
    // f3HistALamALamQinvReal(0),f3HistLamLamMinvReal(0),               
    // f3HistLamProMinvReal(0),f3HistALamALamMinvReal(0),f3HistALamAProMinvReal(0),
    // f2HistBgLamBgLamMeanMinDistProReal(0),f2HistBgLamBgLamMeanMinDistPioReal(0),
    // f2HistBgLamProMeanMinDistProReal(0),f2HistBgALamBgALamMeanMinDistAProReal(0),
    // f2HistBgALamBgALamMeanMinDistPioReal(0),f2HistBgALamAProMeanMinDistAProReal(0),
    // f2HistSftBgLamBgLamMeanMinDistProReal(0),f2HistSftBgLamBgLamMeanMinDistPioReal(0),
    // f2HistSftBgLamProMeanMinDistProReal(0),f2HistSftBgALamBgALamMeanMinDistAProReal(0),
    // f2HistSftBgALamBgALamMeanMinDistPioReal(0),f2HistSftBgALamAProMeanMinDistAProReal(0),
    // f2HistSftIrocBgLamBgLamMeanMinDistProReal(0),f2HistSftIrocBgLamBgLamMeanMinDistPioReal(0),
    // f2HistSftIrocBgLamProMeanMinDistProReal(0),f2HistSftIrocBgALamBgALamMeanMinDistAProReal(0),
    // f2HistSftIrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftIrocBgALamAProMeanMinDistAProReal(0),
    // f2HistSftOrocBgLamBgLamMeanMinDistProReal(0),f2HistSftOrocBgLamBgLamMeanMinDistPioReal(0),
    // f2HistSftOrocBgLamProMeanMinDistProReal(0),f2HistSftOrocBgALamBgALamMeanMinDistAProReal(0),
    // f2HistSftOrocBgALamBgALamMeanMinDistPioReal(0),f2HistSftOrocBgALamAProMeanMinDistAProReal(0),
    BgLamProReal(0),BgALamAProReal(0),
    // f3HistBgLamBgLamQinvReal(0),             
    // f3HistBgALamBgALamQinvReal(0),
    // f2HistLamLamMeanMinDistProMixed(0),f2HistLamLamMeanMinDistPioMixed(0),
    // f2HistLamProMeanMinDistProMixed(0),f2HistALamALamMeanMinDistAProMixed(0),   
    // f2HistALamALamMeanMinDistPioMixed(0),f2HistALamAProMeanMinDistAProMixed(0),
    // f2HistSftLamLamMeanMinDistProMixed(0),f2HistSftLamLamMeanMinDistPioMixed(0),
    // f2HistSftLamProMeanMinDistProMixed(0),f2HistSftALamALamMeanMinDistAProMixed(0),   
    // f2HistSftALamALamMeanMinDistPioMixed(0),f2HistSftALamAProMeanMinDistAProMixed(0),
    // f2HistSftIrocLamLamMeanMinDistProMixed(0),f2HistSftIrocLamLamMeanMinDistPioMixed(0),
    // f2HistSftIrocLamProMeanMinDistProMixed(0),f2HistSftIrocALamALamMeanMinDistAProMixed(0),   
    // f2HistSftIrocALamALamMeanMinDistPioMixed(0),f2HistSftIrocALamAProMeanMinDistAProMixed(0),
    // f2HistSftOrocLamLamMeanMinDistProMixed(0),f2HistSftOrocLamLamMeanMinDistPioMixed(0),
    // f2HistSftOrocLamProMeanMinDistProMixed(0),f2HistSftOrocALamALamMeanMinDistAProMixed(0),   
    // f2HistSftOrocALamALamMeanMinDistPioMixed(0),f2HistSftOrocALamAProMeanMinDistAProMixed(0),
    LamProMixed(0),ALamAProMixed(0),
    // f3HistLamLamQinvMixed(0),                
    // f3HistALamALamQinvMixed(0),f3HistLamLamMinvMixed(0),                
    // f3HistLamProMinvMixed(0),f3HistALamALamMinvMixed(0),f3HistALamAProMinvMixed(0),
    // f2HistBgLamBgLamMeanMinDistProMixed(0),f2HistBgLamBgLamMeanMinDistPioMixed(0),
    // f2HistBgLamProMeanMinDistProMixed(0),f2HistBgALamBgALamMeanMinDistAProMixed(0),
    // f2HistBgALamBgALamMeanMinDistPioMixed(0),f2HistBgALamAProMeanMinDistAProMixed(0),
    // f2HistSftBgLamBgLamMeanMinDistProMixed(0),f2HistSftBgLamBgLamMeanMinDistPioMixed(0),
    // f2HistSftBgLamProMeanMinDistProMixed(0),f2HistSftBgALamBgALamMeanMinDistAProMixed(0),
    // f2HistSftBgALamBgALamMeanMinDistPioMixed(0),f2HistSftBgALamAProMeanMinDistAProMixed(0),
    // f2HistSftIrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftIrocBgLamBgLamMeanMinDistPioMixed(0),
    // f2HistSftIrocBgLamProMeanMinDistProMixed(0),f2HistSftIrocBgALamBgALamMeanMinDistAProMixed(0),
    // f2HistSftIrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftIrocBgALamAProMeanMinDistAProMixed(0),
    // f2HistSftOrocBgLamBgLamMeanMinDistProMixed(0),f2HistSftOrocBgLamBgLamMeanMinDistPioMixed(0),
    // f2HistSftOrocBgLamProMeanMinDistProMixed(0),f2HistSftOrocBgALamBgALamMeanMinDistAProMixed(0),
    // f2HistSftOrocBgALamBgALamMeanMinDistPioMixed(0),f2HistSftOrocBgALamAProMeanMinDistAProMixed(0),
    BgLamProMixed(0),BgALamAProMixed(0)
    // f3HistBgLamBgLamQinvMixed(0),             
    // f3HistBgALamBgALamQinvMixed(0)

{
  // Copy constructor
  printf("Copy constructor not implemented\n");
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda& AliAnalysisTaskProtonLambda::operator=(const AliAnalysisTaskProtonLambda& atpl)
{
  if(this!=&atpl){
  // One operation with the atpl to get rid of the warning unused parameter
  fPrimaryVtxPosition[0]=atpl.fPrimaryVtxPosition[0];
  printf("Assignment operator not implemented\n");
  }
  return *this;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}
//________________________________________________________________________
//
//
//     Classes in the class AliAnalysisTaskProtonLambda
//         FemtoBuffer, FemtoBufferEvent, FemtoBufferV0 and FemtoBufferTrack
//
//________________________________________________________________________
//
//                        FemtoBufferTrack
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferTrack::FemtoBufferTrack():
  fID(65535)
{
  // Standard constructor, initialize everything with values indicating 
  // a track that should not be used
  
  // No idea how to initialize the arrays nicely like the fID(65535)..
  for (UChar_t i=0;i<3;i++){
    fP[i]=-9999.;
    for (UChar_t j=0;j<9;j++){
      //      fXglobal[j][i]=-9999.;
      fXshifted[j][i]=-9999.;
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferTrack::FemtoBufferTrack(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]):
  fID(65535)  
{
  // Constructor

  // Use the function to have the code in one place
  Set(track,bfield,priVtx);
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferTrack::GetGlobalPositionAtGlobalRadii(const AliAODTrack *track, const Float_t bfield){
  // Function not used, do dummy operations to get rid of warnings
  Float_t a=bfield;
  a=track->P();
  
  // // Gets the global position of the track at nine different radii in the TPC
  // // track is the track you want to propagate
  // // bfield is the magnetic field of your event
  // // globalPositionsAtRadii is the array of global positions in the radii and xyz
  
  // // We have two versions of the two track resolution plots in our proton-lambda task:
  // // a) with all events shifted to (0,0,0), b) without shift.
  // // For a) we should compare the tracks at shifted radii,
  // // for b) we should still use the global radii. This function here is for b).

  // // Initialize the array to something indicating there was no propagation
  // for(Int_t i=0;i<9;i++){
  //   for(Int_t j=0;j<3;j++){
  //     fXglobal[i][j]=-9999.;
  //   }
  // }

  //  // Make a copy of the track to not change parameters of the track
  // AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  // //  printf("\nAfter CopyFromVTrack\n");
  // //  etp.Print();
 
  // // The global position of the the track
  // Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // // Counter for which radius we want
  // Int_t iR=0; 
  // // The radii at which we get the global positions
  // // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  // // Compare squared radii for faster code
  // Float_t RSquaredWanted[9]={85.*85.,105.*105.,125.*125.,145.*145.,165.*165.,
  // 			     185.*185.,205.*205.,225.*225.,245.*245.}; 
  // // The global radius we are at, squared. Compare squared radii for faster code
  // Float_t globalRadiusSquared=0;

  // // Propagation is done in local x of the track
  // for (Float_t x = 58.;x<247.;x+=1.){
  //   // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
  //   // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
  //   // the track is straight, i.e. has inifinite pt and doesn't get bent. 
  //   // If the track's momentum is smaller than infinite, it will develop a y-component,
  //   // which adds to the global radius

  //   // Stop if the propagation was not succesful. This can happen for low pt tracks
  //   // that don't reach outer radii
  //   if(!etp.PropagateTo(x,bfield))break;
  //   etp.GetXYZ(xyz); // GetXYZ returns global coordinates

  //   // No shifting for global radii
  //   globalRadiusSquared = (xyz[0])*(xyz[0])
  //                       + (xyz[1])*(xyz[1]);

  //   // Roughly reached the radius we want
  //   if(globalRadiusSquared > RSquaredWanted[iR]){
      
  //     // Bigger loop has bad precision, we're nearly one centimeter too far, 
  //     // go back in small steps.
  //     while (globalRadiusSquared>RSquaredWanted[iR]){
  // 	x-=.1;
  // 	//	printf("propagating to x %5.2f\n",x);
  // 	if(!etp.PropagateTo(x,bfield))break;
  // 	etp.GetXYZ(xyz); // GetXYZ returns global coordinates

  // 	// No shifting for global radii
  // 	globalRadiusSquared = (xyz[0])*(xyz[0])
  // 	                    + (xyz[1])*(xyz[1]);
  //     }
  //     //      printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",TMath::Sqrt(globalRadiusSquared),x,xyz[0],xyz[1],xyz[2]);
  //     fXglobal[iR][0]=xyz[0];
  //     fXglobal[iR][1]=xyz[1];
  //     fXglobal[iR][2]=xyz[2];
  //     // Indicate we want the next radius    
  //     iR+=1;
  //   }
  //   if(iR>=8){
  //     // TPC edge reached
  //     return;
  //   }
  // }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferTrack::GetShiftedPositionAtShiftedRadii(const AliAODTrack *track, const Float_t bfield, const Float_t priVtx[3]){
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  // globalPositionsAtRadii is the array of global positions in the radii and xyz
  
  // Initialize the array to something indicating there was no propagation
  for(Int_t i=0;i<9;i++){
    for(Int_t j=0;j<3;j++){
      fXshifted[i][j]=-9999.;
    }
  }

   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp; etp.CopyFromVTrack(track);
  //  printf("\nAfter CopyFromVTrack\n");
  //  etp.Print();
 
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // Counter for which radius we want
  Int_t iR=0; 
  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  // Compare squared radii for faster code
  Float_t RSquaredWanted[9]={85.*85.,105.*105.,125.*125.,145.*145.,165.*165.,
			     185.*185.,205.*205.,225.*225.,245.*245.}; 
  // The shifted radius we are at, squared. Compare squared radii for faster code
  Float_t shiftedRadiusSquared=0;

  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Without shifting the primary vertex to (0.,0.,0.) the next line would just be
    // WRONG: globalRadiusSquared = xyz[0]*xyz[0]+xyz[1]*xyz[1];
    // but as we shift the primary vertex we want to compare positions at shifted radii.
    // I can't draw in ASCII but please take a piece of paper and just visualize it once.

    // Changing plus to minus on July10th2012
    shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                         + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted[iR]){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted[iR]){
	x-=.1;
	//	printf("propagating to x %5.2f\n",x);
	if(!etp.PropagateTo(x,bfield))break;
	etp.GetXYZ(xyz); // GetXYZ returns global coordinates
	// Added the shifting also here on July11th2012
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      //      printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",TMath::Sqrt(globalRadiusSquared),x,xyz[0],xyz[1],xyz[2]);
      fXshifted[iR][0]=xyz[0]-priVtx[0];
      fXshifted[iR][1]=xyz[1]-priVtx[1];
      fXshifted[iR][2]=xyz[2]-priVtx[2];
      // Indicate we want the next radius    
      iR+=1;
    }
    if(iR>=8){
      // TPC edge reached
      return;
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Double_t priVtx[3]){
  // Overloaded function
  Float_t priVtxPos[3]={static_cast<Float_t>(priVtx[0]),static_cast<Float_t>(priVtx[1]),static_cast<Float_t>(priVtx[2])};
  Set(track,bfield,priVtxPos);
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferTrack::Set(const AliAODTrack *track,const Float_t bfield,const Float_t priVtx[3]){
  // Set the properties of this to the AliAODtrack
  //
  //    UShort_t fID;               //! Unique track id (->AliAODTrack.h), UShort_t goes to 65000
  //    Double_t fP[3];             //! Momentum of track
  //    Float_t  fXglobal[9][3];    //! Global positions at different global radii
  //    Float_t  fXshifted[9][3];   //! Shifted positions at different shifted radii


  // Set the ID, a good ID also indicates to use the track
  if(track->GetID() >=0){
    // global tracks, i.e. v0 daughters
    fID = track->GetID();
  }
  else {
    // e.g. tpc only tracks, i.e. primary protons
    fID = -track->GetID()-1;

  }
  // Set the momentum
  track->PxPyPz(fP);  
  //  GetGlobalPositionAtGlobalRadii(track,bfield);
  GetShiftedPositionAtShiftedRadii(track,bfield,priVtx);

}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferTrack::FemtoBufferTrack(const FemtoBufferTrack& fbt):
  fID(fbt.fID)
 {
  // Copy constructor

  for (UChar_t i=0;i<3;i++){
    fP[i]=fbt.fP[i];
    for (UChar_t j=0;j<9;j++){
      //      fXglobal[j][i]=fbt.fXglobal[j][i];
      fXshifted[j][i]=fbt.fXshifted[j][i];
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferTrack& AliAnalysisTaskProtonLambda::FemtoBufferTrack::operator=(const FemtoBufferTrack& fbt){
  // Assignment operator, from wikipedia :)
  
  // Protect against self-assignment
  if(this != &fbt){
    fID = fbt.fID;
    for (UChar_t i=0;i<3;i++){
      fP[i]=fbt.fP[i];
      for (UChar_t j=0;j<9;j++){
	//	fXglobal[j][i]=fbt.fXglobal[j][i];
	fXshifted[j][i]=fbt.fXshifted[j][i];
      }
    }
  }
  // By convention, always return *this (Could it be the convention is called c++?)
  return *this;
}
//________________________________________________________________________
//
//                        FemtoBufferV0
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferV0::FemtoBufferV0():
  fCosPoint(-9999.),
  fPosDaughter(),
  fNegDaughter()
{
  // Dummy constructor, set everything so it
  // indicates a V0 which should not be used
  fP[0]=-9999.;
  fP[1]=-9999.;
  fP[2]=-9999.;
  // C++11 provides initializer lists, it should work like
  //class C
  //{
  //int x[4];
  //public:
  //C(): x{0,1,2,3} {}
  //};
  // and http://clang.llvm.org/cxx_status.html says, they have it in clang 3.1,
  // but it doesn't seem to work! :/

}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferV0::FemtoBufferV0(const AliAODv0 *v0, const AliAODTrack *posDaughter, const AliAODTrack *negDaughter, const Double_t bfield, Double_t priVtxPos[3]):
  fCosPoint(-9999.),
  fPosDaughter(),
  fNegDaughter()
{
  // Constructor, set the properties of this to these of the AliAODv0

  // Use Set function to keep code in one place. Only constant data member
  // would require the FemtoBuff() : fbla(), fblup() {} method
  Set(v0,posDaughter,negDaughter,bfield,priVtxPos);
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferV0::Set(const AliAODv0 *v0, const AliAODTrack *posDaughter, const AliAODTrack *negDaughter, const Double_t bfield, Double_t priVtxPos[3])
{
  // Set the properties of this to these of the AliAODv0
  fCosPoint=v0->CosPointingAngle(priVtxPos);
  v0->PxPyPz(fP);
  // printf("Set px %3.2f, py %3.2f, pz %3.2f\n"
  // 	 ,fP[0],fP[1],fP[2]
  // 	 );
  // The daughters
  fPosDaughter.Set(posDaughter,bfield,priVtxPos);
  fNegDaughter.Set(negDaughter,bfield,priVtxPos);
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferV0::FemtoBufferV0(const FemtoBufferV0 &fbv):
  fCosPoint(fbv.fCosPoint),
  fPosDaughter(fbv.fPosDaughter),
  fNegDaughter(fbv.fNegDaughter)
  //,fP{fbv.fP[0],fbv.fP[1],fbv.fP[2]} // C++11
{
  // Copy constructor
  fP[0] = fbv.fP[0]; // C++03
  fP[1] = fbv.fP[1];
  fP[2] = fbv.fP[2];
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferV0& AliAnalysisTaskProtonLambda::FemtoBufferV0::operator=(const FemtoBufferV0 &fbv){
  // Assignment operator

  // Protect against self-assignment
  if(this != &fbv){
    fCosPoint=fbv.fCosPoint;
    fP[0]=fbv.fP[0];
    fP[1]=fbv.fP[1];
    fP[2]=fbv.fP[2];
    fPosDaughter=fbv.fPosDaughter;
    fNegDaughter=fbv.fNegDaughter;
  }
  return *this;
}
//________________________________________________________________________
//
//                        FemtoBufferEvent
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent::FemtoBufferEvent():
  fPriTrackLim(0),fV0Lim(0)
  ,fProTracks(0),fAProTracks(0)
  ,fLamTracks(0),fALamTracks(0)
  ,fBgLamTracks(0),fBgALamTracks(0)
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-9999.)
{
  // Standard constructor, all pointer to zero
  fPriVtxPos[0]=-9999.;
  fPriVtxPos[1]=-9999.;
  fPriVtxPos[2]=-9999.;

  printf("This constructor has zero size in the arrays!\n");
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent::FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff,const Double_t bfield,const Double_t priVtxPos[3]):
  fPriTrackLim(priTrackBuff),fV0Lim(V0Buff)
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgALamTracks(new FemtoBufferV0[fV0Lim])
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-bfield)
  //  ,fPriVtxPos{priVtxPos[0],priVtxPos[1],priVtxPos[2]} // This is C++11
{
  // Constructor.
  fPriVtxPos[0] = priVtxPos[0]; // This is some old C++
  fPriVtxPos[1] = priVtxPos[1];
  fPriVtxPos[2] = priVtxPos[2];
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent::FemtoBufferEvent(const UShort_t priTrackBuff,const UShort_t V0Buff):
  fPriTrackLim(priTrackBuff),fV0Lim(V0Buff)
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgALamTracks(new FemtoBufferV0[fV0Lim])
  ,fNProTracks(0),fNAProTracks(0),fNLamTracks(0),fNALamTracks(0)
  ,fNBgLamTracks(0),fNBgALamTracks(0)
  ,fBfield(-9999.)
  //  ,fPriVtxPos{-9999.,-9999.,-9999.} // This is C++11
{  
  // Constructor. fBfield and fPriVtxPos not needed yet, can be set later.
  fPriVtxPos[0] = -9999.; // This is C++03
  fPriVtxPos[1] = -9999.;
  fPriVtxPos[2] = -9999.;

  //  printf("constructed eventwith NBgLam: %u\n",fNBgLamTracks);
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent::FemtoBufferEvent(const FemtoBufferEvent &fbe):
  fPriTrackLim(fbe.GetPriTrackLim())
  ,fV0Lim(fbe.GetV0Lim())
  ,fProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fAProTracks(new FemtoBufferTrack[fPriTrackLim])
  ,fLamTracks (new FemtoBufferV0[fV0Lim])
  ,fALamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgLamTracks(new FemtoBufferV0[fV0Lim])
  ,fBgALamTracks(new FemtoBufferV0[fV0Lim])
  ,fNProTracks(fbe.GetNPro()),fNAProTracks(fbe.GetNAPro())
  ,fNLamTracks(fbe.GetNLam()),fNALamTracks(fbe.GetNALam())
  ,fNBgLamTracks(fbe.GetNBgLam()),fNBgALamTracks(fbe.GetNBgALam())
  ,fBfield(fbe.GetBfield())
{
  // Copy constructor
  fbe.GetVtxPos(fPriVtxPos);
  // Avoid to much creation and deletion of objects
  UShort_t i;
  // Copy the primary tracks
  for (i=0;i<fPriTrackLim;i++){
    fProTracks[i]=fbe.fProTracks[i];
    fAProTracks[i]=fbe.fAProTracks[i];
  }
  // Copy the V0s
  for (i=0;i<fV0Lim;i++){
    fLamTracks[i]=fbe.fLamTracks[i];
    fALamTracks[i]=fbe.fALamTracks[i];
    fBgLamTracks[i]=fbe.fBgLamTracks[i];
    fBgALamTracks[i]=fbe.fBgALamTracks[i];
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent& AliAnalysisTaskProtonLambda::FemtoBufferEvent::operator=(const FemtoBufferEvent &fbe){
  // Assignment operator

  // Protect against self-assignment
  if(this!=&fbe){
    // Well, we use arrays of a constant size to avoid
    // excessive memory allocation and won't give this up.
    // So we'll only copy as much as fits on the left side
    // from the right side.
    // DON'T COPY THE ARRAY SIZES fV0Lim AND fPriTrackLim !!!
    if(fPriTrackLim < fbe.GetPriTrackLim() 
       || fV0Lim < fbe.GetV0Lim()){
      // AliWarning(Form("Trying to assign too big event (buffer %d/%d) to"
      // 		    " this (buffer %d/%d). Only partially copying.",
      // 		    fbe.GetPriTrackLim(),fbe.GetV0Lim(),
      // 		    fPriTrackLim,fV0Lim));
      printf("Trying to assign too big event (buffer %d/%d) to"
    		    " this (buffer %d/%d). Only partially copying.\n",
	     fbe.GetPriTrackLim(),fbe.GetV0Lim(),
	     fPriTrackLim,fV0Lim);
    }
    // Always start with the easy stuff :)
    fbe.GetVtxPos(fPriVtxPos);
    fBfield = fbe.GetBfield();
    // Number of tracks is minimum of array size of 'this'
    // and the number of tracks from the right side
    fNProTracks = TMath::Min(fPriTrackLim,fbe.GetNPro());
    fNAProTracks = TMath::Min(fPriTrackLim,fbe.GetNAPro());
    fNLamTracks = TMath::Min(fV0Lim,fbe.GetNLam());
    fNALamTracks = TMath::Min(fV0Lim,fbe.GetNALam());
    fNBgLamTracks = TMath::Min(fV0Lim,fbe.GetNBgLam());
    fNBgALamTracks = TMath::Min(fV0Lim,fbe.GetNBgALam());
    
    // Avoid creation and deletion of 'i' for every loop
    UShort_t i;
    // Copy primary tracks. No need to set a 'bad track'
    // flag for the entries above GetNPro() (...) as
    // above everything is bad by definition.
    // Protons
    for (i=0;i<GetNPro();i++)
      fProTracks[i]=fbe.fProTracks[i];
    // Anti-protons
    for (i=0;i<GetNAPro();i++)
      fAProTracks[i]=fbe.fAProTracks[i];
    // Copy the V0s 
    // Lambdas
    for (i=0;i<GetNLam();i++){
      fLamTracks[i]=fbe.fLamTracks[i];
    }
    // Anti-lambdas
    for (i=0;i<GetNALam();i++){
      fALamTracks[i]=fbe.fALamTracks[i];
    }
    // Background lambdas
    for (i=0;i<GetNBgLam();i++){
      fBgLamTracks[i]=fbe.fBgLamTracks[i];
    }
    // Background anti-lambdas
    for (i=0;i<GetNBgALam();i++){
      fBgALamTracks[i]=fbe.fBgALamTracks[i];
    }  
  }
  return *this;
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBufferEvent::~FemtoBufferEvent(){
  // Destructor

  // Delete the arrays of tracks,
  // note the [] with the delete
  if(fProTracks){
    delete[] fProTracks;
    fProTracks=0;
  }
  if(fAProTracks){
    delete[] fAProTracks;
    fAProTracks=0;
  }
  if(fLamTracks){
    delete[] fLamTracks;
    fLamTracks=0;
  }
  if(fALamTracks){
    delete[] fALamTracks;
    fALamTracks=0;
  }
  if(fBgLamTracks){
    delete[] fBgLamTracks;
    fBgLamTracks=0;
  }
  if(fBgALamTracks){
    delete[] fBgALamTracks;
    fBgALamTracks=0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::Reset(const Double_t bfield, const Double_t priVtxPos[3]){
  // Reset the old event, i.e., make clear 'here is no info'
  // by setting the 'number of stored ...' to zero
  fNProTracks=0;
  fNAProTracks=0;
  fNLamTracks=0;
  fNALamTracks=0;
  fNBgLamTracks=0;
  fNBgALamTracks=0;
  
  // And set the new event properties 
  fBfield = bfield;
  fPriVtxPos[0]=priVtxPos[0];
  fPriVtxPos[1]=priVtxPos[1];
  fPriVtxPos[2]=priVtxPos[2];
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddPro(const AliAODTrack *track){
  // Add a proton to this event

  // Check whether there is still space in the array
  if(fNProTracks > fPriTrackLim-1){
    // AliWarning(Form("Cannot add proton, array size (%d) too small"
    // 		    ,fPriTrackLim));
    printf("Cannot add proton, array size (%d) too small\n"
    		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array
  fProTracks[fNProTracks].Set(track,fBfield,fPriVtxPos);
  fNProTracks++;
  //  printf("Added proton %d/%d\n",fNProTracks,fPriTrackLim);

}  
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddAPro(const AliAODTrack *track){
  // Add a anti-proton to this event

  // Check whether there is still space in the array
  if(fNAProTracks > fPriTrackLim-1){
    // AliWarning(Form("Cannot add anti-proton, array size (%d) too small"
    // 		    ,fPriTrackLim));
    printf("Cannot add anti-proton, array size (%d) too small\n"
		    ,fPriTrackLim);
    return;
  }
  // Add the V0 at the end of the array
  fAProTracks[fNAProTracks].Set(track,fBfield,fPriVtxPos);
  fNAProTracks++;
}  
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNLamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add lambda, array size (%d) too small"
    // 		    ,fV0Lim));
    printf("Cannot add lambda, array size (%d) too small"
		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fLamTracks[fNLamTracks].Set(v0,posDaughter,negDaughter,
			      fBfield,fPriVtxPos);
  fNLamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNALamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add anti-lambda, array size (%d) too small"
    // 		    ,fV0Lim));
    printf("Cannot add anti-lambda, array size (%d) too small\n"
    		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fALamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNALamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddBgLam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNBgLamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add background lambda,"
    // 		    " array size (%d) too small"
    // 		    ,fV0Lim));
    // printf("Cannot add background lambda,"
    // 	   "already stored %d" 
    // 	   " array size (%d) too small\n"
    // 	   ,fNBgLamTracks
    // 	   ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fBgLamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNBgLamTracks++;
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBufferEvent::AddBgALam(const AliAODv0 *v0,const AliAODTrack *posDaughter,const AliAODTrack *negDaughter){
  // Adds a lambda with it's daughters to the event

  // Check whether there is still space in the array
  if(fNBgALamTracks > fV0Lim-1){
    // AliWarning(Form("Cannot add background anti-lambda,"
    // 		    " array size (%d) too small"
    // 		    ,fV0Lim));
    //    printf("Cannot add background anti-lambda,"
    //		    " array size (%d) too small\n"
    //		    ,fV0Lim);
    return;
  }
 
  // Add the V0 at the end of the array
  fBgALamTracks[fNALamTracks].Set(v0,posDaughter,negDaughter,
				fBfield,fPriVtxPos);
  fNBgALamTracks++;
}
//________________________________________________________________________
//
//                        FemtoBuffer
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBuffer::FemtoBuffer() :
  fkZvertexBins(0),
  fkCentBins(0),
  fkMixBuffSize(0),
  fkPriTrackLim(0),
  fkV0Lim(0),
  fZvertexAxis(0),
  fCentAxis(0),
  fCurEvt(0),
  fEC(0)
{
  // Dummy constructor, create arrays with zero size
  // Note that some data member are constant, you
  // won't be able to create the FemtoBuffer first with this
  // constructor and then set the appropiate size.
  
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBuffer::FemtoBuffer(const UChar_t ZvertexBins,const UChar_t CentBins,const UChar_t MixBuff,const UShort_t PriTrackLim,const UShort_t V0Lim, const Float_t AbsZvertexCut,const Float_t CentCut) :
  fkZvertexBins(ZvertexBins),
  fkCentBins(CentBins),
  fkMixBuffSize(MixBuff),
  fkPriTrackLim(PriTrackLim),
  fkV0Lim(V0Lim),
  fZvertexAxis(new TAxis(fkZvertexBins,-AbsZvertexCut,AbsZvertexCut)),
  fCentAxis(new TAxis (fkCentBins,0.0,CentCut)),
  fCurEvt(new FemtoBufferEvent *[fkMixBuffSize]),
  fEC(new FemtoBufferEvent ***[fkZvertexBins])
{
  // Constructor, creates at once all events with all tracks
  //  printf ("Creating with pritracklim %d and v0lim %d\n",fkPriTrackLim,fkV0Lim);

  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new FemtoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new FemtoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new FemtoBufferEvent(fkPriTrackLim,fkV0Lim);
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBuffer::FemtoBuffer(const AliAnalysisTaskProtonLambda::FemtoBuffer &fb) :
  fkZvertexBins(fb.fkZvertexBins),
  fkCentBins(fb.fkCentBins),
  fkMixBuffSize(fb.fkMixBuffSize),
  fkPriTrackLim(fb.fkPriTrackLim),
  fkV0Lim(fb.fkV0Lim),
  fZvertexAxis(new TAxis(*(fb.fZvertexAxis))),
  fCentAxis(new TAxis (*(fb.fCentAxis))),
  fCurEvt(new FemtoBufferEvent *[fkMixBuffSize]),
  fEC(new FemtoBufferEvent ***[fkZvertexBins])
{
  // Copy constructor. Linux complains not having this and 
  // compiling this task with aliroot

  printf("FemtoBuffer ctor not tested yet, be cautious\n");
  
  // Create the array step by step
  // Bins in z of the primary vertex position. Do this as
  // the detector looks different from a different z coordinate
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    fEC[iZBin] = new FemtoBufferEvent **[fkCentBins];
    // Bins in centrality
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      fEC[iZBin][iCentBin] = new FemtoBufferEvent *[fkMixBuffSize];
      // The number of events to keep for one mixing class
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	// Create an event to hold the info for mixing
	fEC[iZBin][iCentBin][iMixBuff] = new FemtoBufferEvent(*(fb.fEC[iZBin][iCentBin][iMixBuff]));
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBuffer& AliAnalysisTaskProtonLambda::FemtoBuffer::operator=(const AliAnalysisTaskProtonLambda::FemtoBuffer& fb){
  //Assignment operator
  if(this!=&fb){
    printf("FemtoBuffer assignment operator not implemented\n");
  }
  return *this;
  
}
//________________________________________________________________________
AliAnalysisTaskProtonLambda::FemtoBuffer::~FemtoBuffer(){
  // Destructor
  // The axes to fin the correct bins
  if(fZvertexAxis){
    delete fZvertexAxis;
    fZvertexAxis=0;
  }
  if(fCentAxis){
    delete fCentAxis;
    fCentAxis=0;
  }
  // fCurEvt is an array of pointer
  if(fCurEvt){
    delete[] fCurEvt;
    fCurEvt=0;
  }
  // Delete all the events and the pointer to them
  for (UChar_t iZBin=0;iZBin<fkZvertexBins;iZBin++){
    for (UChar_t iCentBin=0;iCentBin<fkCentBins;iCentBin++){
      for(UChar_t iMixBuff=0;iMixBuff<fkMixBuffSize;iMixBuff++){
	if(fEC[iZBin][iCentBin][iMixBuff]){
	  delete fEC[iZBin][iCentBin][iMixBuff];
	  fEC[iZBin][iCentBin][iMixBuff]=0;
	}
      }
      if(fEC[iZBin][iCentBin]){
	delete fEC[iZBin][iCentBin];
	fEC[iZBin][iCentBin]=0;
      }
    }
    if(fEC[iZBin]){
      delete fEC[iZBin];
      fEC[iZBin]=0;
    }
  }
  if(fEC){
    delete fEC;
    fEC=0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBuffer::ShiftAndAdd(AliAODEvent *evt){
  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly
  Double_t priVtxPos[3];
  evt->GetPrimaryVertex()->GetXYZ(priVtxPos);
  //  printf("Mag field: %f\n",evt->GetMagneticField());
  ShiftAndAdd(evt->GetMagneticField(),
	      priVtxPos,
	      evt->GetCentrality()->GetCentralityPercentileUnchecked("V0M"));
}
//________________________________________________________________________
void AliAnalysisTaskProtonLambda::FemtoBuffer::ShiftAndAdd(const Double_t bfield,const Double_t priVtxPos[3],const Float_t centrality){
  // Shift the events in the appropiate centrality / zvertex bin and set the 
  // current event pointer correctly

  // Find the correct centrality/zvertex bin 
  const UChar_t ZvertexBin = fZvertexAxis->FindFixBin(priVtxPos[2]) - 1; // -1 for array starting at 0
  const UChar_t CentBin = fCentAxis->FindFixBin(centrality) - 1;// -1 for array starting at 0

  // The new current event is the old last event
  fCurEvt[0] = fEC[ZvertexBin][CentBin][fkMixBuffSize-1];

  // Shift the pointer, starting from the back
  UChar_t iMix;
  for(iMix=fkMixBuffSize-1;iMix>0;iMix--){
    fEC[ZvertexBin][CentBin][iMix] = fEC[ZvertexBin][CentBin][iMix-1];
  }
  // And reset the zero'th one
  fEC[ZvertexBin][CentBin][0] = fCurEvt[0];
  fEC[ZvertexBin][CentBin][0]->Reset(bfield,priVtxPos);
  // Also set the pointer to the other events..
  for (iMix=1;iMix<fkMixBuffSize;iMix++){
    fCurEvt[iMix] = fEC[ZvertexBin][CentBin][iMix];
  }
}
