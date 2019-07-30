/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
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

//*****************************************************
//   Class AliEventShape
//   author: Jaroslav Priskin
//   
//
//*****************************************************

#include "AliAnalysisSphericityTask.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TDatabasePDG.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <AliMultSelection.h> 
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include <AliVTrack.h>
#include <TTreeStream.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <TParticle.h>
#include "AliPPVsMultUtils.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliCentrality.h" 
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 
#include "AliESDcascade.h"
#include "AliTPCPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliTransverseEventShape.h"

// STL includes
#include <iostream>
using namespace std;

ClassImp(AliAnalysisSphericityTask)
//_____________________________________________________________________________
AliAnalysisSphericityTask::AliAnalysisSphericityTask():
AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fPPVsMultUtils(0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
  fCent(1),            
  fUseHybrid(0x0),
  fTrackFilterHybrid1(0x0),
  fTrackFilterHybrid2(0x0),
  fTrackFilterESA(0x0),
  fMinMultESA(0x0),
  fSizeStepESA(0x0),
  fIsAbsEtaESA(0x0),
  fEtaMaxCutESA(0x0),
  fEtaMinCutESA(0x0),
  fPtMaxCutESA(0x0),
  fPtMinCutESA(0x0),
  fNrec(100), 
  fCentEst("V0M"),
  fPIDMode(kSigma),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  fkExtraSelections(0),
  fkExtraSelectionsCut(0),
  fInvMassCutLambda(0),
  fInvMassCutKaon(0),
  fInvMassCutXi(0),
  fInvMassCutOmega(0),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),       //def: kFALSE
  fStrangeness(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinCent(0.0),
  fMaxCent(100.0),
  fStoreMcIn(kFALSE),       //def: kFALSE takes from AnalysisMC
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fOADBPath(),
  fPIDResponse(0x0),
  fOldRun(0),
  fRecoPass(0),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0), 
  fMinDaughterTpcClusters(80),
  fkQualityCutTPCrefit(0),
  hVtxBeforeCuts(0x0), 
  hVtxAfterCuts(0x0),
  hn1(0x0),
  hn2(0x0),
  hso(0x0),
  hcent(0x0),
  hst(0x0),
  fhphiSt(0x0),     
  fhptSt(0x0),    
  fhetaSt(0x0),        
  hsphericityVSmulti(0x0),     
  hsphericityVSpT(0x0),  
  hmultiplicityVSpT(0x0),
  hmultiplicityVSpTbin(0x0),     
  hsphericityVScent(0x0),     
  hsphericityVSMEANpT(0x0),      
  hmultiplicityVSMEANpT(0x0),  
  hmultiplicityVSMEANpTbin(0x0),          
  hsphericity(0x0),    
  hetaVSphi(0x0),       
  hetaVSphiJET(0x0),     
  hetaVSphiISO(0x0),        
  hetaVSphiMID(0x0),        
  HMultRef(0x0),
  fHistTrackMultiplicity(0),
  fHistMassKaon(0), 
  fHistMassXiMinus(0),
  fHistMassXiPlus(0),
  fHistMassOmegaMinus(0),
  fHistMassOmegaPlus(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fNSigma1(0), 
  fNRatio1(0), 
  fNSigma2(0), 
  fNRatio2(0), 
  fNBoundP(0),
  fV0Cuts(0),
  fCascadeCuts(0),
  fRerunVertexers(0),
  fMaxV0Rapidity(1.),
  fHistXiArmenteros(0),
  fHistV0Armenteros(0),
  fStrangePart(kFALSE),
  fFindKaon(kFALSE),
  fFindLambda(kFALSE),
  fFindAntiLambda(kFALSE),
  fFindXiPlus(kFALSE),
  fFindXiMinus(kFALSE),
  fFindOmegaPlus(kFALSE),
  fFindOmegaMinus(kFALSE),
  hsphericityVSmultiL(0x0),
  hsphericityL(0x0),
  hsphericityVScentL(0x0),     
  hsphericityVSMEANpTL(0x0),
  hsphericityVSpTL(0x0),
  hmultiplicityVSMEANpTL(0x0),  
  hmultiplicityVSMEANpTbinL(0x0),
  hmultiplicityVSpTL(0x0),
  hmultiplicityVSpTbinL(0x0), 
  hsphericityVSmultiAL(0x0),
  hsphericityAL(0x0),
  hsphericityVScentAL(0x0),     
  hsphericityVSMEANpTAL(0x0),
  hsphericityVSpTAL(0x0),
  hmultiplicityVSMEANpTAL(0x0),  
  hmultiplicityVSMEANpTbinAL(0x0),
  hmultiplicityVSpTAL(0x0),
  hmultiplicityVSpTbinAL(0x0),
  hsphericityVSmultiLO(0x0),
  hsphericityLO(0x0),
  hsphericityVScentLO(0x0),     
  hsphericityVSMEANpTLO(0x0),
  hsphericityVSpTLO(0x0),
  hmultiplicityVSMEANpTLO(0x0),  
  hmultiplicityVSMEANpTbinLO(0x0),
  hmultiplicityVSpTLO(0x0),
  hmultiplicityVSpTbinLO(0x0), 
  hsphericityVSmultiALLO(0x0),
  hsphericityALLO(0x0),
  hsphericityVScentALLO(0x0),     
  hsphericityVSMEANpTALLO(0x0),
  hsphericityVSpTALLO(0x0),
  hmultiplicityVSMEANpTALLO(0x0),  
  hmultiplicityVSMEANpTbinALLO(0x0),
  hmultiplicityVSpTALLO(0x0),
  hmultiplicityVSpTbinALLO(0x0),
  hsphericityVSmultiX(0x0),
  hsphericityX(0x0),
  hsphericityVScentX(0x0),     
  hsphericityVSMEANpTX(0x0),
  hsphericityVSpTX(0x0),
  hmultiplicityVSMEANpTX(0x0),  
  hmultiplicityVSMEANpTbinX(0x0),
  hmultiplicityVSpTX(0x0),
  hmultiplicityVSpTbinX(0x0), 
  hsphericityVSmultiXP(0x0),
  hsphericityXP(0x0),
  hsphericityVScentXP(0x0),     
  hsphericityVSMEANpTXP(0x0),
  hsphericityVSpTXP(0x0),
  hmultiplicityVSMEANpTXP(0x0),  
  hmultiplicityVSMEANpTbinXP(0x0),
  hmultiplicityVSpTXP(0x0),
  hmultiplicityVSpTbinXP(0x0),
  hsphericityVSmultiXN(0x0),
  hsphericityXN(0x0),
  hsphericityVScentXN(0x0),     
  hsphericityVSMEANpTXN(0x0),
  hsphericityVSpTXN(0x0),
  hmultiplicityVSMEANpTXN(0x0),  
  hmultiplicityVSMEANpTbinXN(0x0),
  hmultiplicityVSpTXN(0x0),
  hmultiplicityVSpTbinXN(0x0), 
  hsphericityVSmultiXNXP(0x0),
  hsphericityXNXP(0x0),
  hsphericityVScentXNXP(0x0),     
  hsphericityVSMEANpTXNXP(0x0),
  hsphericityVSpTXNXP(0x0),
  hmultiplicityVSMEANpTXNXP(0x0),  
  hmultiplicityVSMEANpTbinXNXP(0x0),
  hmultiplicityVSpTXNXP(0x0),
  hmultiplicityVSpTbinXNXP(0x0),
  hsphericityVSmultiO(0x0),
  hsphericityO(0x0),
  hsphericityVScentO(0x0),     
  hsphericityVSMEANpTO(0x0),
  hsphericityVSpTO(0x0),
  hmultiplicityVSMEANpTO(0x0),  
  hmultiplicityVSMEANpTbinO(0x0),
  hmultiplicityVSpTO(0x0),
  hmultiplicityVSpTbinO(0x0), 
  hsphericityVSmultiOP(0x0),
  hsphericityOP(0x0),
  hsphericityVScentOP(0x0),     
  hsphericityVSMEANpTOP(0x0),
  hsphericityVSpTOP(0x0),
  hmultiplicityVSMEANpTOP(0x0),  
  hmultiplicityVSMEANpTbinOP(0x0),
  hmultiplicityVSpTOP(0x0),
  hmultiplicityVSpTbinOP(0x0),
  hsphericityVSmultiON(0x0),
  hsphericityON(0x0),
  hsphericityVScentON(0x0),     
  hsphericityVSMEANpTON(0x0),
  hsphericityVSpTON(0x0),
  hmultiplicityVSMEANpTON(0x0),  
  hmultiplicityVSMEANpTbinON(0x0),
  hmultiplicityVSpTON(0x0),
  hmultiplicityVSpTbinON(0x0), 
  hsphericityVSmultiONOP(0x0),
  hsphericityONOP(0x0),
  hsphericityVScentONOP(0x0),     
  hsphericityVSMEANpTONOP(0x0),
  hsphericityVSpTONOP(0x0),
  hmultiplicityVSMEANpTONOP(0x0),  
  hmultiplicityVSMEANpTbinONOP(0x0),
  hmultiplicityVSpTONOP(0x0),
  hmultiplicityVSpTbinONOP(0x0),
  hsphericityVSmultiKA(0x0),
  hsphericityKA(0x0),
  hsphericityVScentKA(0x0),     
  hsphericityVSMEANpTKA(0x0),
  hsphericityVSpTKA(0x0),
  hmultiplicityVSMEANpTKA(0x0),  
  hmultiplicityVSMEANpTbinKA(0x0),
  hmultiplicityVSpTKA(0x0),
  hmultiplicityVSpTbinKA(0x0),

  hsphericityVSmultiL2(0x0),
  hsphericityL2(0x0),
  hsphericityVScentL2(0x0),     
  hsphericityVSMEANpTL2(0x0),
  hsphericityVSpTL2(0x0),
  hmultiplicityVSMEANpTL2(0x0),  
  hmultiplicityVSMEANpTbinL2(0x0),
  hmultiplicityVSpTL2(0x0),
  hmultiplicityVSpTbinL2(0x0),
  hsphericityVSmultiALO(0x0),
  hsphericityALO(0x0),
  hsphericityVScentALO(0x0),     
  hsphericityVSMEANpTALO(0x0),
  hsphericityVSpTALO(0x0),
  hmultiplicityVSMEANpTALO(0x0),  
  hmultiplicityVSMEANpTbinALO(0x0),
  hmultiplicityVSpTALO(0x0),
  hmultiplicityVSpTbinALO(0x0),
  hsphericityVSmultiLOO(0x0),
  hsphericityLOO(0x0),
  hsphericityVScentLOO(0x0),     
  hsphericityVSMEANpTLOO(0x0),
  hsphericityVSpTLOO(0x0),
  hmultiplicityVSMEANpTLOO(0x0),  
  hmultiplicityVSMEANpTbinLOO(0x0),
  hmultiplicityVSpTLOO(0x0),
  hmultiplicityVSpTbinLOO(0x0), 
  hsphericityVSmultiALLOO(0x0),
  hsphericityALLOO(0x0),
  hsphericityVScentALLOO(0x0),     
  hsphericityVSMEANpTALLOO(0x0),
  hsphericityVSpTALLOO(0x0),
  hmultiplicityVSMEANpTALLOO(0x0),  
  hmultiplicityVSMEANpTbinALLOO(0x0),
  hmultiplicityVSpTALLOO(0x0),
  hmultiplicityVSpTbinALLOO(0x0), 
  hsphericityVSmultiXO(0x0),
  hsphericityXO(0x0),
  hsphericityVScentXO(0x0),     
  hsphericityVSMEANpTXO(0x0),
  hsphericityVSpTXO(0x0),
  hmultiplicityVSMEANpTXO(0x0),  
  hmultiplicityVSMEANpTbinXO(0x0),
  hmultiplicityVSpTXO(0x0),
  hmultiplicityVSpTbinXO(0x0),
  hsphericityVSmultiXPO(0x0),
  hsphericityXPO(0x0),
  hsphericityVScentXPO(0x0),     
  hsphericityVSMEANpTXPO(0x0),
  hsphericityVSpTXPO(0x0),
  hmultiplicityVSMEANpTXPO(0x0),  
  hmultiplicityVSMEANpTbinXPO(0x0),
  hmultiplicityVSpTXPO(0x0),
  hmultiplicityVSpTbinXPO(0x0),
  hsphericityVSmultiXNO(0x0),
  hsphericityXNO(0x0),
  hsphericityVScentXNO(0x0),     
  hsphericityVSMEANpTXNO(0x0),
  hsphericityVSpTXNO(0x0),
  hmultiplicityVSMEANpTXNO(0x0),  
  hmultiplicityVSMEANpTbinXNO(0x0),
  hmultiplicityVSpTXNO(0x0),
  hmultiplicityVSpTbinXNO(0x0), 
  hsphericityVSmultiXNXPO(0x0),
  hsphericityXNXPO(0x0),
  hsphericityVScentXNXPO(0x0),     
  hsphericityVSMEANpTXNXPO(0x0),
  hsphericityVSpTXNXPO(0x0),
  hmultiplicityVSMEANpTXNXPO(0x0),  
  hmultiplicityVSMEANpTbinXNXPO(0x0),
  hmultiplicityVSpTXNXPO(0x0),
  hmultiplicityVSpTbinXNXPO(0x0),
  hsphericityVSmultiOO(0x0),
  hsphericityOO(0x0),
  hsphericityVScentOO(0x0),     
  hsphericityVSMEANpTOO(0x0),
  hsphericityVSpTOO(0x0),
  hmultiplicityVSMEANpTOO(0x0),  
  hmultiplicityVSMEANpTbinOO(0x0),
  hmultiplicityVSpTOO(0x0),
  hmultiplicityVSpTbinOO(0x0),
  hsphericityVSmultiOPO(0x0),
  hsphericityOPO(0x0),
  hsphericityVScentOPO(0x0),     
  hsphericityVSMEANpTOPO(0x0),
  hsphericityVSpTOPO(0x0),
  hmultiplicityVSMEANpTOPO(0x0),  
  hmultiplicityVSMEANpTbinOPO(0x0),
  hmultiplicityVSpTOPO(0x0),
  hmultiplicityVSpTbinOPO(0x0),
  hsphericityVSmultiONO(0x0),
  hsphericityONO(0x0),
  hsphericityVScentONO(0x0),     
  hsphericityVSMEANpTONO(0x0),
  hsphericityVSpTONO(0x0),
  hmultiplicityVSMEANpTONO(0x0),  
  hmultiplicityVSMEANpTbinONO(0x0),
  hmultiplicityVSpTONO(0x0),
  hmultiplicityVSpTbinONO(0x0), 
  hsphericityVSmultiONOPO(0x0),
  hsphericityONOPO(0x0),
  hsphericityVScentONOPO(0x0),     
  hsphericityVSMEANpTONOPO(0x0),
  hsphericityVSpTONOPO(0x0),
  hmultiplicityVSMEANpTONOPO(0x0),  
  hmultiplicityVSMEANpTbinONOPO(0x0),
  hmultiplicityVSpTONOPO(0x0),
  hmultiplicityVSpTbinONOPO(0x0),
  hsphericityVSmultiKAO(0x0),
  hsphericityKAO(0x0),
  hsphericityVScentKAO(0x0),     
  hsphericityVSMEANpTKAO(0x0),
  hsphericityVSpTKAO(0x0),
  hmultiplicityVSMEANpTKAO(0x0),  
  hmultiplicityVSMEANpTbinKAO(0x0),
  hmultiplicityVSpTKAO(0x0),
  hmultiplicityVSpTbinKAO(0x0),
  fESASelection(0x0)
  
{
  //default constructor
}

AliAnalysisSphericityTask::AliAnalysisSphericityTask(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fPPVsMultUtils(0),
  fMC(0x0),
  fMCStack(0x0),
  fMCArray(0x0),
  fCent(1),                   
  fUseHybrid(0x0),
  fTrackFilterHybrid1(0x0),
  fTrackFilterHybrid2(0x0),
  fTrackFilterESA(0x0),
  fMinMultESA(0x0),
  fSizeStepESA(0x0),
  fIsAbsEtaESA(0x0),
  fEtaMaxCutESA(0x0),
  fEtaMinCutESA(0x0),
  fPtMaxCutESA(0x0),
  fPtMinCutESA(0x0),
  fNrec(100),             
  fCentEst("V0M"),
  fPIDMode(kSigma),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  fAnalysisPbPb(kFALSE),
  fkExtraSelections(0),
  fkExtraSelectionsCut(0),
  fInvMassCutLambda(0),
  fInvMassCutKaon(0),
  fInvMassCutXi(0),
  fInvMassCutOmega(0),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),       //def: kFALSE
  fStrangeness(kFALSE),   
  fVtxCut(10.0),  
  fEtaCut(0.9),  
  fMinCent(0.0),
  fMaxCent(100.0),
  fStoreMcIn(kFALSE),       //def: kFALSE takes from AnalysisMC
  fMcProcessType(-999),
  fTriggeredEventMB(-999),
  fVtxStatus(-999),
  fZvtx(-999),
  fZvtxMC(-999),
  fRun(-999),
  fEventId(-999),
  fListOfObjects(0),
  fOADBPath(),
  fPIDResponse(0x0),
  fOldRun(0),
  fRecoPass(0),
  fMinDaughterTpcClusters(80),
  fkQualityCutTPCrefit(0), 
  hVtxBeforeCuts(0x0),
  hVtxAfterCuts(0x0),
  hn1(0x0),
  hn2(0x0),
  hso(0x0),
  hcent(0x0),
  hst(0x0),
  fhphiSt(0x0),      
  fhptSt(0x0),     
  fhetaSt(0x0),        
  hsphericityVSmulti(0x0),        
  hsphericityVSpT(0x0),      
  hmultiplicityVSpT(0x0),
  hmultiplicityVSpTbin(0x0),      
  hsphericityVScent(0x0),     
  hsphericityVSMEANpT(0x0), 
  hmultiplicityVSMEANpT(0x0),  
  hmultiplicityVSMEANpTbin(0x0),        
  hsphericity(0x0),      
  hetaVSphi(0x0),     
  hetaVSphiJET(0x0),     
  hetaVSphiISO(0x0),    
  hetaVSphiMID(0x0),    
  HMultRef(0x0),
  fHistTrackMultiplicity(0),
  fHistMassKaon(0),
  fHistMassXiMinus(0),
  fHistMassXiPlus(0),
  fHistMassOmegaMinus(0),
  fHistMassOmegaPlus(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fNSigma1(0), 
  fNRatio1(0), 
  fNSigma2(0), 
  fNRatio2(0),  
  fNBoundP(0),
  fV0Cuts(0),
  fCascadeCuts(0),
  fRerunVertexers(0),
  fMaxV0Rapidity(1.),
  fHistXiArmenteros(0),
  fHistV0Armenteros(0),
  fStrangePart(kFALSE),
  fFindKaon(kFALSE),
  fFindLambda(kFALSE),
  fFindAntiLambda(kFALSE),
  fFindXiPlus(kFALSE),
  fFindXiMinus(kFALSE),
  fFindOmegaPlus(kFALSE),
  fFindOmegaMinus(kFALSE),
  hsphericityVSmultiL(0x0),
  hsphericityL(0x0),
  hsphericityVScentL(0x0),     
  hsphericityVSMEANpTL(0x0),
  hsphericityVSpTL(0x0),
  hmultiplicityVSMEANpTL(0x0),  
  hmultiplicityVSMEANpTbinL(0x0),
  hmultiplicityVSpTL(0x0),
  hmultiplicityVSpTbinL(0x0),
  hsphericityVSmultiAL(0x0),
  hsphericityAL(0x0),
  hsphericityVScentAL(0x0),     
  hsphericityVSMEANpTAL(0x0),
  hsphericityVSpTAL(0x0),
  hmultiplicityVSMEANpTAL(0x0),  
  hmultiplicityVSMEANpTbinAL(0x0),
  hmultiplicityVSpTAL(0x0),
  hmultiplicityVSpTbinAL(0x0),
  hsphericityVSmultiLO(0x0),
  hsphericityLO(0x0),
  hsphericityVScentLO(0x0),     
  hsphericityVSMEANpTLO(0x0),
  hsphericityVSpTLO(0x0),
  hmultiplicityVSMEANpTLO(0x0),  
  hmultiplicityVSMEANpTbinLO(0x0),
  hmultiplicityVSpTLO(0x0),
  hmultiplicityVSpTbinLO(0x0), 
  hsphericityVSmultiALLO(0x0),
  hsphericityALLO(0x0),
  hsphericityVScentALLO(0x0),     
  hsphericityVSMEANpTALLO(0x0),
  hsphericityVSpTALLO(0x0),
  hmultiplicityVSMEANpTALLO(0x0),  
  hmultiplicityVSMEANpTbinALLO(0x0),
  hmultiplicityVSpTALLO(0x0),
  hmultiplicityVSpTbinALLO(0x0), 
  hsphericityVSmultiX(0x0),
  hsphericityX(0x0),
  hsphericityVScentX(0x0),     
  hsphericityVSMEANpTX(0x0),
  hsphericityVSpTX(0x0),
  hmultiplicityVSMEANpTX(0x0),  
  hmultiplicityVSMEANpTbinX(0x0),
  hmultiplicityVSpTX(0x0),
  hmultiplicityVSpTbinX(0x0),
  hsphericityVSmultiXP(0x0),
  hsphericityXP(0x0),
  hsphericityVScentXP(0x0),     
  hsphericityVSMEANpTXP(0x0),
  hsphericityVSpTXP(0x0),
  hmultiplicityVSMEANpTXP(0x0),  
  hmultiplicityVSMEANpTbinXP(0x0),
  hmultiplicityVSpTXP(0x0),
  hmultiplicityVSpTbinXP(0x0),
  hsphericityVSmultiXN(0x0),
  hsphericityXN(0x0),
  hsphericityVScentXN(0x0),     
  hsphericityVSMEANpTXN(0x0),
  hsphericityVSpTXN(0x0),
  hmultiplicityVSMEANpTXN(0x0),  
  hmultiplicityVSMEANpTbinXN(0x0),
  hmultiplicityVSpTXN(0x0),
  hmultiplicityVSpTbinXN(0x0), 
  hsphericityVSmultiXNXP(0x0),
  hsphericityXNXP(0x0),
  hsphericityVScentXNXP(0x0),     
  hsphericityVSMEANpTXNXP(0x0),
  hsphericityVSpTXNXP(0x0),
  hmultiplicityVSMEANpTXNXP(0x0),  
  hmultiplicityVSMEANpTbinXNXP(0x0),
  hmultiplicityVSpTXNXP(0x0),
  hmultiplicityVSpTbinXNXP(0x0),
  hsphericityVSmultiO(0x0),
  hsphericityO(0x0),
  hsphericityVScentO(0x0),     
  hsphericityVSMEANpTO(0x0),
  hsphericityVSpTO(0x0),
  hmultiplicityVSMEANpTO(0x0),  
  hmultiplicityVSMEANpTbinO(0x0),
  hmultiplicityVSpTO(0x0),
  hmultiplicityVSpTbinO(0x0),
  hsphericityVSmultiOP(0x0),
  hsphericityOP(0x0),
  hsphericityVScentOP(0x0),     
  hsphericityVSMEANpTOP(0x0),
  hsphericityVSpTOP(0x0),
  hmultiplicityVSMEANpTOP(0x0),  
  hmultiplicityVSMEANpTbinOP(0x0),
  hmultiplicityVSpTOP(0x0),
  hmultiplicityVSpTbinOP(0x0),
  hsphericityVSmultiON(0x0),
  hsphericityON(0x0),
  hsphericityVScentON(0x0),     
  hsphericityVSMEANpTON(0x0),
  hsphericityVSpTON(0x0),
  hmultiplicityVSMEANpTON(0x0),  
  hmultiplicityVSMEANpTbinON(0x0),
  hmultiplicityVSpTON(0x0),
  hmultiplicityVSpTbinON(0x0), 
  hsphericityVSmultiONOP(0x0),
  hsphericityONOP(0x0),
  hsphericityVScentONOP(0x0),     
  hsphericityVSMEANpTONOP(0x0),
  hsphericityVSpTONOP(0x0),
  hmultiplicityVSMEANpTONOP(0x0),  
  hmultiplicityVSMEANpTbinONOP(0x0),
  hmultiplicityVSpTONOP(0x0),
  hmultiplicityVSpTbinONOP(0x0),
  hsphericityVSmultiKA(0x0),
  hsphericityKA(0x0),
  hsphericityVScentKA(0x0),     
  hsphericityVSMEANpTKA(0x0),
  hsphericityVSpTKA(0x0),
  hmultiplicityVSMEANpTKA(0x0),  
  hmultiplicityVSMEANpTbinKA(0x0),
  hmultiplicityVSpTKA(0x0),
  hmultiplicityVSpTbinKA(0x0),

  hsphericityVSmultiL2(0x0),
  hsphericityL2(0x0),
  hsphericityVScentL2(0x0),     
  hsphericityVSMEANpTL2(0x0),
  hsphericityVSpTL2(0x0),
  hmultiplicityVSMEANpTL2(0x0),  
  hmultiplicityVSMEANpTbinL2(0x0),
  hmultiplicityVSpTL2(0x0),
  hmultiplicityVSpTbinL2(0x0),
  hsphericityVSmultiALO(0x0),
  hsphericityALO(0x0),
  hsphericityVScentALO(0x0),     
  hsphericityVSMEANpTALO(0x0),
  hsphericityVSpTALO(0x0),
  hmultiplicityVSMEANpTALO(0x0),  
  hmultiplicityVSMEANpTbinALO(0x0),
  hmultiplicityVSpTALO(0x0),
  hmultiplicityVSpTbinALO(0x0),
  hsphericityVSmultiLOO(0x0),
  hsphericityLOO(0x0),
  hsphericityVScentLOO(0x0),     
  hsphericityVSMEANpTLOO(0x0),
  hsphericityVSpTLOO(0x0),
  hmultiplicityVSMEANpTLOO(0x0),  
  hmultiplicityVSMEANpTbinLOO(0x0),
  hmultiplicityVSpTLOO(0x0),
  hmultiplicityVSpTbinLOO(0x0), 
  hsphericityVSmultiALLOO(0x0),
  hsphericityALLOO(0x0),
  hsphericityVScentALLOO(0x0),     
  hsphericityVSMEANpTALLOO(0x0),
  hsphericityVSpTALLOO(0x0),
  hmultiplicityVSMEANpTALLOO(0x0),  
  hmultiplicityVSMEANpTbinALLOO(0x0),
  hmultiplicityVSpTALLOO(0x0),
  hmultiplicityVSpTbinALLOO(0x0), 
  hsphericityVSmultiXO(0x0),
  hsphericityXO(0x0),
  hsphericityVScentXO(0x0),     
  hsphericityVSMEANpTXO(0x0),
  hsphericityVSpTXO(0x0),
  hmultiplicityVSMEANpTXO(0x0),  
  hmultiplicityVSMEANpTbinXO(0x0),
  hmultiplicityVSpTXO(0x0),
  hmultiplicityVSpTbinXO(0x0),
  hsphericityVSmultiXPO(0x0),
  hsphericityXPO(0x0),
  hsphericityVScentXPO(0x0),     
  hsphericityVSMEANpTXPO(0x0),
  hsphericityVSpTXPO(0x0),
  hmultiplicityVSMEANpTXPO(0x0),  
  hmultiplicityVSMEANpTbinXPO(0x0),
  hmultiplicityVSpTXPO(0x0),
  hmultiplicityVSpTbinXPO(0x0),
  hsphericityVSmultiXNO(0x0),
  hsphericityXNO(0x0),
  hsphericityVScentXNO(0x0),     
  hsphericityVSMEANpTXNO(0x0),
  hsphericityVSpTXNO(0x0),
  hmultiplicityVSMEANpTXNO(0x0),  
  hmultiplicityVSMEANpTbinXNO(0x0),
  hmultiplicityVSpTXNO(0x0),
  hmultiplicityVSpTbinXNO(0x0), 
  hsphericityVSmultiXNXPO(0x0),
  hsphericityXNXPO(0x0),
  hsphericityVScentXNXPO(0x0),     
  hsphericityVSMEANpTXNXPO(0x0),
  hsphericityVSpTXNXPO(0x0),
  hmultiplicityVSMEANpTXNXPO(0x0),  
  hmultiplicityVSMEANpTbinXNXPO(0x0),
  hmultiplicityVSpTXNXPO(0x0),
  hmultiplicityVSpTbinXNXPO(0x0),
  hsphericityVSmultiOO(0x0),
  hsphericityOO(0x0),
  hsphericityVScentOO(0x0),     
  hsphericityVSMEANpTOO(0x0),
  hsphericityVSpTOO(0x0),
  hmultiplicityVSMEANpTOO(0x0),  
  hmultiplicityVSMEANpTbinOO(0x0),
  hmultiplicityVSpTOO(0x0),
  hmultiplicityVSpTbinOO(0x0),
  hsphericityVSmultiOPO(0x0),
  hsphericityOPO(0x0),
  hsphericityVScentOPO(0x0),     
  hsphericityVSMEANpTOPO(0x0),
  hsphericityVSpTOPO(0x0),
  hmultiplicityVSMEANpTOPO(0x0),  
  hmultiplicityVSMEANpTbinOPO(0x0),
  hmultiplicityVSpTOPO(0x0),
  hmultiplicityVSpTbinOPO(0x0),
  hsphericityVSmultiONO(0x0),
  hsphericityONO(0x0),
  hsphericityVScentONO(0x0),     
  hsphericityVSMEANpTONO(0x0),
  hsphericityVSpTONO(0x0),
  hmultiplicityVSMEANpTONO(0x0),  
  hmultiplicityVSMEANpTbinONO(0x0),
  hmultiplicityVSpTONO(0x0),
  hmultiplicityVSpTbinONO(0x0), 
  hsphericityVSmultiONOPO(0x0),
  hsphericityONOPO(0x0),
  hsphericityVScentONOPO(0x0),     
  hsphericityVSMEANpTONOPO(0x0),
  hsphericityVSpTONOPO(0x0),
  hmultiplicityVSMEANpTONOPO(0x0),  
  hmultiplicityVSMEANpTbinONOPO(0x0),
  hmultiplicityVSpTONOPO(0x0),
  hmultiplicityVSpTbinONOPO(0x0),
  hsphericityVSmultiKAO(0x0),
  hsphericityKAO(0x0),
  hsphericityVScentKAO(0x0),     
  hsphericityVSMEANpTKAO(0x0),
  hsphericityVSpTKAO(0x0),
  hmultiplicityVSMEANpTKAO(0x0),  
  hmultiplicityVSMEANpTbinKAO(0x0),
  hmultiplicityVSpTKAO(0x0),
  hmultiplicityVSpTbinKAO(0x0),
  fESASelection(0x0)
{
  // Default constructor (should not be used)
  DefineOutput(1, TList::Class());
}

AliAnalysisSphericityTask::~AliAnalysisSphericityTask() {
  //
  // Destructor
  //
  if (fESASelection) {
    delete fESASelection;
    fESASelection = 0x0;
  }  
}

//______________________________________________________________________________
void AliAnalysisSphericityTask::UserCreateOutputObjects()
{ 
  
  fRandom = new TRandom(0); // 0 means random seed

  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");
  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(lmcEvtHandler) inputHandler->CreatePIDResponse(kTRUE);
  else
  inputHandler->CreatePIDResponse(kFALSE);
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");
  fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

  fTPCpid = dynamic_cast<AliTPCPIDResponse&>(fPIDResponse->GetTPCResponse());

  //Helper
  if(! fESASelection ) {
    fESASelection = new AliTransverseEventShape();
    fESASelection->SetUseHybridESA(fUseHybrid);
    fESASelection->SetTrackFilterESAHyb1(fTrackFilterHybrid1);
    fESASelection->SetTrackFilterESAHyb2(fTrackFilterHybrid2);
    fESASelection->SetTrackFilterESA(fTrackFilterESA);
    fESASelection->SetMinMultForESA(fMinMultESA);
    fESASelection->SetStepSizeESA(fSizeStepESA);
    fESASelection->SetIsEtaAbsESA(fIsAbsEtaESA);
    fESASelection->SetTrackEtaMinESA(fEtaMinCutESA);
    fESASelection->SetTrackEtaMaxESA(fEtaMaxCutESA);
    fESASelection->SetTrackPtMinESA(fPtMinCutESA);
    fESASelection->SetTrackPtMaxESA(fPtMaxCutESA);
    fESASelection->Init();
  }
  
  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //bining
  //used recommended bins from analysis notes 
    const Int_t nMultbins = 9;
    Double_t Multbins[nMultbins+1]={0.0,35.0,76.0,149.0,261.0,426.0,649.0,966.0,1294.0,1601.0};     //multiplicity bins
  
    const Int_t nBinsSt = 5;
    Double_t BinSt[nBinsSt+1]={0.0, 0.1, 0.3, 0.5, 0.7, 1.0};                                       //recomended sphericity bins
    
  //
  // Histograms
  if(fStrangeness){
  hn1=new TH1D("ESD - only strange","Statistics of Events; ;Events",11,-1,10);
  hn1->GetXaxis()->SetBinLabel(2, "Processed");
	hn1->GetXaxis()->SetBinLabel(3, "PileUpTrigger cut");;
	hn1->GetXaxis()->SetBinLabel(4, "Trigger cut");
  hn1->GetXaxis()->SetBinLabel(5, "Vertex position cut");
  hn1->GetXaxis()->SetBinLabel(6, "Centrality cut");
  hn1->GetXaxis()->SetBinLabel(7, "Sphericity");
  fListOfObjects->Add(hn1);
  }
  else{
  hn1=new TH1D("ESD","Statistics of Events; ;Events",11,-1,10);
  hn1->GetXaxis()->SetBinLabel(2, "Processed");
	hn1->GetXaxis()->SetBinLabel(3, "PileUpTrigger cut");;
	hn1->GetXaxis()->SetBinLabel(4, "Trigger cut");
  hn1->GetXaxis()->SetBinLabel(5, "Vertex position cut");
  hn1->GetXaxis()->SetBinLabel(6, "Centrality cut");
  hn1->GetXaxis()->SetBinLabel(7, "Sphericity");
  fListOfObjects->Add(hn1); 
  }
  hn2=new TH1D("EVENTS (after cut)","Statistics of Tracks; ;Tracks",11,-1,10);
  hn2->GetXaxis()->SetBinLabel(2, "Processed");
	hn2->GetXaxis()->SetBinLabel(3, "Cuts in pseudorapidity");;
	hn2->GetXaxis()->SetBinLabel(4, "Cuts in #it{p}_{T}");
	hn2->GetXaxis()->SetBinLabel(5, "Quality cuts");
  fListOfObjects->Add(hn2);
  
  hVtxBeforeCuts = new TH1D("hVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(hVtxBeforeCuts);  
  hVtxAfterCuts = new TH1D("hVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fListOfObjects->Add(hVtxAfterCuts);

  hcent = new TH1D("hcent",";centrality; events",110,0.0,110.0);
  fListOfObjects->Add(hcent);
  hso = new TH1D("hso",";spherocity; entries",200,0.0,1.0);
  fListOfObjects->Add(hso);
  fhetaSt = new TH1D("fhetaSt","ESD;#eta; entries",30,-1.,1.);
  fListOfObjects->Add(fhetaSt);
	fhphiSt = new TH1D("fhphiSt","ESD;#phi (rad); entries", 64, 0, 2*TMath::Pi());
  fListOfObjects->Add(fhphiSt);
	fhptSt  = new TH1D("fhptSt","ESD;#it{p}_{T} (GeV/#it{c}); entries", 500, 0, 5);
  fListOfObjects->Add(fhptSt);

  hst = new TH1D("hst",";Sphericity #it{S}_{T}; entries",200,0.0,1.0);
  fListOfObjects->Add(hst);

  hetaVSphi = new TH2D(" Profile"," #eta VS #phi;#eta; #phi (rad)", 30, -1., 1., 64, 0, 2*TMath::Pi());                  
  fListOfObjects->Add(hetaVSphi);
  hetaVSphiISO = new TH2D(" Isotropic profile"," #eta VS #phi;#eta; #phi (rad)", 30, -1., 1., 64, 0, 2*TMath::Pi());               
  fListOfObjects->Add(hetaVSphiISO);
  hetaVSphiJET = new TH2D(" Jetty-like profile"," #eta VS #phi;#eta; #phi (rad)", 30, -1., 1., 64, 0, 2*TMath::Pi());                   
  fListOfObjects->Add(hetaVSphiJET);
  hetaVSphiMID = new TH2D(" Mid-Spherocity profile"," #eta VS #phi;#eta; #phi (rad)", 30, -1., 1., 64, 0, 2*TMath::Pi());                     
  fListOfObjects->Add(hetaVSphiMID);
  HMultRef= new TH1D("HMultRef",";Multiplicity; entries",1600,0.0,1600);
  fListOfObjects->Add(HMultRef);

  /////////////////////////////////////////////////////////////// MAIN /////////////////////////////////////////////////////////////////////////////

  hsphericityVSmulti = new TH2D("Sphericity w Bins"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmulti);
  hsphericity = new TH2D("Sphericity wo Bins"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericity); 

  hsphericityVSpT = new TH2D("Sphericity wo Bins to pT"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpT);
  hsphericityVSMEANpT = new TH2D("Sphericity wo Bins to <pT>"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpT);
  
  hsphericityVScent = new TH2D("Sphericity wo Bins to Centrality"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScent);

  hmultiplicityVSpTbin = new TH2D("Multiplicity w Bins to pT"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbin);  
  hmultiplicityVSMEANpTbin = new TH2D("Multiplicity w Bins to <pT>"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbin);
  hmultiplicityVSpT = new TH2D("Multiplicity wo Bins to pT"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpT);  
  hmultiplicityVSMEANpT = new TH2D("Multiplicity wo Bins to <pT>"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpT);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if (!fHistTrackMultiplicity) {
      fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000);
    fListOfObjects->Add(fHistTrackMultiplicity);
  }

  if (!fHistMassKaon) {
    fHistMassKaon = new TH1F("fHistMassKaon","#Kappa^{0} candidates;M(#Kappa^{0}) (GeV/c^{2});Counts",150,0.485,0.515);
    fListOfObjects->Add(fHistMassKaon);
  }
  if (!fHistMassLambda) {
    fHistMassLambda = new TH1F("fHistMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts",150,1.05,1.2);
    fListOfObjects->Add(fHistMassLambda);
  }
  if (!fHistMassAntiLambda) {
    fHistMassAntiLambda = new TH1F("fHistMassAntiLambda","#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts",150,1.05,1.2);
    fListOfObjects->Add(fHistMassAntiLambda);
  }
  if (! fHistMassXiMinus) {
    fHistMassXiMinus = new TH1F("fHistMassXiMinus","#Xi^{-} candidates;M( #Lambda , #pi^{-} ) (GeV/c^{2});Counts", 200,1.25,1.4);
    fListOfObjects->Add(fHistMassXiMinus);
  }
  if (! fHistMassXiPlus) {
    fHistMassXiPlus = new TH1F("fHistMassXiPlus","#Xi^{+} candidates;M( #bar{#Lambda}^{0} , #pi^{+} )(GeV/c^{2});Counts",200,1.25,1.4);
    fListOfObjects->Add(fHistMassXiPlus);
  }
  if (! fHistMassOmegaMinus) {
    fHistMassOmegaMinus = new TH1F("fHistMassOmegaMinus","#Omega^{-} candidates;M( #Lambda , K^{-} ) (GeV/c^{2});Counts", 120, 1.62, 1.74);
    fListOfObjects->Add(fHistMassOmegaMinus);
  }
  if (! fHistMassOmegaPlus) {
    fHistMassOmegaPlus = new TH1F("fHistMassOmegaPlus","#Omega^{+} candidates;M( #bar{#Lambda}^{0} , K^{+} )(GeV/c^{2});Counts",120, 1.62, 1.74);
    fListOfObjects->Add(fHistMassOmegaPlus);
  }
  if(! fHistXiArmenteros) {
	  fHistXiArmenteros = new TH2F( "fHistXiArmenteros", "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
	  fListOfObjects->Add(fHistXiArmenteros);
  }
  if(! fHistV0Armenteros) {
	  fHistV0Armenteros = new TH2F( "fHistV0Armenteros", "#alpha_{Arm}(V0 cand.) Vs Pt_{Arm}(V0 cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)", 140, -1.2, 1.2, 300, 0., 0.3);
	  fListOfObjects->Add(fHistV0Armenteros);
  }
   
   /////////////////////////////////////////////////////////////// Lambda /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiL = new TH2D("Sphericity w Bins #Lambda"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiL);
  hsphericityL = new TH2D("Sphericity wo Bins #Lambda"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityL); 

  hsphericityVSpTL = new TH2D("Sphericity wo Bins to pT #Lambda"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTL);
  hsphericityVSMEANpTL = new TH2D("Sphericity wo Bins to <pT> #Lambda"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTL);
  
  hsphericityVScentL = new TH2D("Sphericity wo Bins to Centrality #Lambda"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentL);

  hmultiplicityVSpTbinL = new TH2D("Multiplicity w Bins to pT #Lambda"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinL);  
  hmultiplicityVSMEANpTbinL = new TH2D("Multiplicity w Bins to <pT> #Lambda"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinL);
  hmultiplicityVSpTL = new TH2D("Multiplicity wo Bins to pT #Lambda"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTL);  
  hmultiplicityVSMEANpTL = new TH2D("Multiplicity wo Bins to <pT> #Lambda"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTL);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Anti Lambda /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiAL = new TH2D("Sphericity w Bins #bar{#Lambda}^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiAL);
  hsphericityAL = new TH2D("Sphericity wo Bins #bar{#Lambda}^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityAL); 

  hsphericityVSpTAL = new TH2D("Sphericity wo Bins to pT #bar{#Lambda}^{0}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTAL);
  hsphericityVSMEANpTAL = new TH2D("Sphericity wo Bins to <pT> #bar{#Lambda}^{0}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTAL);
  
  hsphericityVScentAL = new TH2D("Sphericity wo Bins to Centrality #bar{#Lambda}^{0}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentAL);

  hmultiplicityVSpTbinAL = new TH2D("Multiplicity w Bins to pT #bar{#Lambda}^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinAL);  
  hmultiplicityVSMEANpTbinAL = new TH2D("Multiplicity w Bins to <pT> #bar{#Lambda}^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinAL);
  hmultiplicityVSpTAL = new TH2D("Multiplicity wo Bins to pT #bar{#Lambda}^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTAL);  
  hmultiplicityVSMEANpTAL = new TH2D("Multiplicity wo Bins to <pT> #bar{#Lambda}^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTAL);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Lambda Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiLO = new TH2D("Sphericity w Bins #Lambda^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiLO);
  hsphericityLO = new TH2D("Sphericity wo Bins #Lambda^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityLO); 

  hsphericityVSpTLO = new TH2D("Sphericity wo Bins to pT #Lambda^{0}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTLO);
  hsphericityVSMEANpTLO = new TH2D("Sphericity wo Bins to <pT> #Lambda^{0}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTLO);
  
  hsphericityVScentLO = new TH2D("Sphericity wo Bins to Centrality #Lambda^{0}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentLO);

  hmultiplicityVSpTbinLO = new TH2D("Multiplicity w Bins to pT #Lambda^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinLO);  
  hmultiplicityVSMEANpTbinLO = new TH2D("Multiplicity w Bins to <pT> #Lambda^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinLO);
  hmultiplicityVSpTLO = new TH2D("Multiplicity wo Bins to pT #Lambda^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTLO);  
  hmultiplicityVSMEANpTLO = new TH2D("Multiplicity wo Bins to <pT> #Lambda^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTLO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Lambda or Anti Lambda Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiALLO = new TH2D("Sphericity w Bins #Lambda^{0} or #bar{#Lambda}^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiALLO);
  hsphericityALLO = new TH2D("Sphericity wo Bins #Lambda^{0} or #bar{#Lambda}^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityALLO); 

  hsphericityVSpTALLO = new TH2D("Sphericity wo Bins to pT #Lambda^{0} or #bar{#Lambda}^{0}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTALLO);
  hsphericityVSMEANpTALLO = new TH2D("Sphericity wo Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTALLO);
  
  hsphericityVScentALLO = new TH2D("Sphericity wo Bins to Centrality #Lambda^{0} or #bar{#Lambda}^{0}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentALLO);

  hmultiplicityVSpTbinALLO = new TH2D("Multiplicity w Bins to pT #Lambda^{0} or #bar{#Lambda}^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinALLO);  
  hmultiplicityVSMEANpTbinALLO = new TH2D("Multiplicity w Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinALLO);
  hmultiplicityVSpTALLO = new TH2D("Multiplicity wo Bins to pT #Lambda^{0} or #bar{#Lambda}^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTALLO);  
  hmultiplicityVSMEANpTALLO = new TH2D("Multiplicity wo Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTALLO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////// Xi /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiX = new TH2D("Sphericity w Bins #Xi"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiX);
  hsphericityX = new TH2D("Sphericity wo Bins #Xi"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityX); 

  hsphericityVSpTX = new TH2D("Sphericity wo Bins to pT #Xi"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTX);
  hsphericityVSMEANpTX = new TH2D("Sphericity wo Bins to <pT> #Xi"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTX);
  
  hsphericityVScentX = new TH2D("Sphericity wo Bins to Centrality #Xi"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentX);

  hmultiplicityVSpTbinX = new TH2D("Multiplicity w Bins to pT #Xi"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinX);  
  hmultiplicityVSMEANpTbinX = new TH2D("Multiplicity w Bins to <pT> #Xi"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinX);
  hmultiplicityVSpTX = new TH2D("Multiplicity wo Bins to pT #Xi"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTX);  
  hmultiplicityVSMEANpTX = new TH2D("Multiplicity wo Bins to <pT> #Xi"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTX);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Xi Plus /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXP = new TH2D("Sphericity w Bins #bar{#Xi}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXP);
  hsphericityXP = new TH2D("Sphericity wo Bins #bar{#Xi}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXP); 

  hsphericityVSpTXP = new TH2D("Sphericity wo Bins to pT #bar{#Xi}^{+}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXP);
  hsphericityVSMEANpTXP = new TH2D("Sphericity wo Bins to <pT> #bar{#Xi}^{+}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXP);
  
  hsphericityVScentXP = new TH2D("Sphericity wo Bins to Centrality #bar{#Xi}^{+}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXP);

  hmultiplicityVSpTbinXP = new TH2D("Multiplicity w Bins to pT #bar{#Xi}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXP);  
  hmultiplicityVSMEANpTbinXP = new TH2D("Multiplicity w Bins to <pT> #bar{#Xi}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXP);
  hmultiplicityVSpTXP = new TH2D("Multiplicity wo Bins to pT #bar{#Xi}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXP);  
  hmultiplicityVSMEANpTXP = new TH2D("Multiplicity wo Bins to <pT> #bar{#Xi}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXP);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Xi Minus /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXN = new TH2D("Sphericity w Bins #Xi^{-}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXN);
  hsphericityXN = new TH2D("Sphericity wo Bins #Xi^{-}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXN); 

  hsphericityVSpTXN = new TH2D("Sphericity wo Bins to pT #Xi^{-}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXN);
  hsphericityVSMEANpTXN = new TH2D("Sphericity wo Bins to <pT> #Xi^{-}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXN);
  
  hsphericityVScentXN = new TH2D("Sphericity wo Bins to Centrality #Xi^{-}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXN);

  hmultiplicityVSpTbinXN = new TH2D("Multiplicity w Bins to pT #Xi^{-}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXN);  
  hmultiplicityVSMEANpTbinXN = new TH2D("Multiplicity w Bins to <pT> #Xi^{-}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXN);
  hmultiplicityVSpTXN = new TH2D("Multiplicity wo Bins to pT #Xi^{-}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXN);  
  hmultiplicityVSMEANpTXN = new TH2D("Multiplicity wo Bins to <pT> #Xi^{-}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXN);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Xi Plus or Xi Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXNXP = new TH2D("Sphericity w Bins #Xi^{-} or #bar{#Xi}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXNXP);
  hsphericityXNXP = new TH2D("Sphericity wo Bins #Xi^{-} or #bar{#Xi}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXNXP); 

  hsphericityVSpTXNXP = new TH2D("Sphericity wo Bins to pT #Xi^{-} or #bar{#Xi}^{+}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXNXP);
  hsphericityVSMEANpTXNXP = new TH2D("Sphericity wo Bins to <pT> #Xi^{-} or #bar{#Xi}^{+}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXNXP);
  
  hsphericityVScentXNXP = new TH2D("Sphericity wo Bins to Centrality #Xi^{-} or #bar{#Xi}^{+}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXNXP);

  hmultiplicityVSpTbinXNXP = new TH2D("Multiplicity w Bins to pT #Xi^{-} or #bar{#Xi}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXNXP);  
  hmultiplicityVSMEANpTbinXNXP = new TH2D("Multiplicity w Bins to <pT> #Xi^{-} or #bar{#Xi}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXNXP);
  hmultiplicityVSpTXNXP = new TH2D("Multiplicity wo Bins to pT #Xi^{-} or #bar{#Xi}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXNXP);  
  hmultiplicityVSMEANpTXNXP = new TH2D("Multiplicity wo Bins to <pT> #Xi^{-} or #bar{#Xi}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXNXP);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Omega /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiO = new TH2D("Sphericity w Bins #Omega"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiO);
  hsphericityO = new TH2D("Sphericity wo Bins #Omega"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityO); 

  hsphericityVSpTO = new TH2D("Sphericity wo Bins to pT #Omega"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTO);
  hsphericityVSMEANpTO = new TH2D("Sphericity wo Bins to <pT> #Omega"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTO);
  
  hsphericityVScentO = new TH2D("Sphericity wo Bins to Centrality #Omega"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentO);

  hmultiplicityVSpTbinO = new TH2D("Multiplicity w Bins to pT #Omega"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinO);  
  hmultiplicityVSMEANpTbinO = new TH2D("Multiplicity w Bins to <pT> #Omega"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinO);
  hmultiplicityVSpTO = new TH2D("Multiplicity wo Bins to pT #Omega"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTO);  
  hmultiplicityVSMEANpTO = new TH2D("Multiplicity wo Bins to <pT> #Omega"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Omega Plus /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiOP = new TH2D("Sphericity w Bins #bar{#Omega}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiOP);
  hsphericityOP = new TH2D("Sphericity wo Bins #bar{#Omega}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityOP); 

  hsphericityVSpTOP = new TH2D("Sphericity wo Bins to pT #bar{#Omega}^{+}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTOP);
  hsphericityVSMEANpTOP = new TH2D("Sphericity wo Bins to <pT> #bar{#Omega}^{+}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTOP);
  
  hsphericityVScentOP = new TH2D("Sphericity wo Bins to Centrality #bar{#Omega}^{+}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentOP);

  hmultiplicityVSpTbinOP = new TH2D("Multiplicity w Bins to pT #bar{#Omega}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinOP);  
  hmultiplicityVSMEANpTbinOP = new TH2D("Multiplicity w Bins to <pT> #bar{#Omega}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinOP);
  hmultiplicityVSpTOP = new TH2D("Multiplicity wo Bins to pT #bar{#Omega}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTOP);  
  hmultiplicityVSMEANpTOP = new TH2D("Multiplicity wo Bins to <pT> #bar{#Omega}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTOP);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Omega Minus /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiON = new TH2D("Sphericity w Bins #Omega^{-}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiON);
  hsphericityON = new TH2D("Sphericity wo Bins #Omega^{-}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityON); 

  hsphericityVSpTON = new TH2D("Sphericity wo Bins to pT #Omega^{-}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTON);
  hsphericityVSMEANpTON = new TH2D("Sphericity wo Bins to <pT> #Omega^{-}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTON);
  
  hsphericityVScentON = new TH2D("Sphericity wo Bins to Centrality #Omega^{-}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentON);

  hmultiplicityVSpTbinON = new TH2D("Multiplicity w Bins to pT #Omega^{-}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinON);  
  hmultiplicityVSMEANpTbinON = new TH2D("Multiplicity w Bins to <pT> #Omega^{-}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinON);
  hmultiplicityVSpTON = new TH2D("Multiplicity wo Bins to pT #Omega^{-}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTON);  
  hmultiplicityVSMEANpTON = new TH2D("Multiplicity wo Bins to <pT> #Omega^{-}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTON);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Omega Plus or Omega Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiONOP = new TH2D("Sphericity w Bins #Omega^{-} or #bar{#Omega}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiONOP);
  hsphericityONOP = new TH2D("Sphericity wo Bins #Omega^{-} or #bar{#Omega}^{+}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityONOP); 

  hsphericityVSpTONOP = new TH2D("Sphericity wo Bins to pT #Omega^{-} or #bar{#Omega}^{+}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTONOP);
  hsphericityVSMEANpTONOP = new TH2D("Sphericity wo Bins to <pT> #Omega^{-} or #bar{#Omega}^{+}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTONOP);
  
  hsphericityVScentONOP = new TH2D("Sphericity wo Bins to Centrality #Omega^{-} or #bar{#Omega}^{+}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentONOP);

  hmultiplicityVSpTbinONOP = new TH2D("Multiplicity w Bins to pT #Omega^{-} or #bar{#Omega}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinONOP);  
  hmultiplicityVSMEANpTbinONOP = new TH2D("Multiplicity w Bins to <pT> #Omega^{-} or #bar{#Omega}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinONOP);
  hmultiplicityVSpTONOP = new TH2D("Multiplicity wo Bins to pT #Omega^{-} or #bar{#Omega}^{+}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTONOP);  
  hmultiplicityVSMEANpTONOP = new TH2D("Multiplicity wo Bins to <pT> #Omega^{-} or #bar{#Omega}^{+}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTONOP);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Kaon /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiKA = new TH2D("Sphericity w Bins #Kappa^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiKA);
  hsphericityKA = new TH2D("Sphericity wo Bins #Kappa^{0}"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityKA); 

  hsphericityVSpTKA = new TH2D("Sphericity wo Bins to pT #Kappa^{0}"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTKA);
  hsphericityVSMEANpTKA = new TH2D("Sphericity wo Bins to <pT> #Kappa^{0}"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTKA);
  
  hsphericityVScentKA = new TH2D("Sphericity wo Bins to Centrality #Kappa^{0}"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentKA);

  hmultiplicityVSpTbinKA = new TH2D("Multiplicity w Bins to pT #Kappa^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinKA);  
  hmultiplicityVSMEANpTbinKA = new TH2D("Multiplicity w Bins to <pT> #Kappa^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinKA);
  hmultiplicityVSpTKA = new TH2D("Multiplicity wo Bins to pT #Kappa^{0}"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTKA);  
  hmultiplicityVSMEANpTKA = new TH2D("Multiplicity wo Bins to <pT> #Kappa^{0}"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTKA);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////// ONLY WO OTHER STRANGE //////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////// Lambda /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiL2 = new TH2D("Sphericity w Bins #Lambda Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiL2);
  hsphericityL2 = new TH2D("Sphericity wo Bins #Lambda Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityL2); 

  hsphericityVSpTL2 = new TH2D("Sphericity wo Bins to pT #Lambda Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTL2);
  hsphericityVSMEANpTL2 = new TH2D("Sphericity wo Bins to <pT> #Lambda Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTL2);
  
  hsphericityVScentL2 = new TH2D("Sphericity wo Bins to Centrality #Lambda Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentL2);

  hmultiplicityVSpTbinL2 = new TH2D("Multiplicity w Bins to pT #Lambda Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinL2);  
  hmultiplicityVSMEANpTbinL2 = new TH2D("Multiplicity w Bins to <pT> #Lambda Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinL2);
  hmultiplicityVSpTL2 = new TH2D("Multiplicity wo Bins to pT #Lambda Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTL2);  
  hmultiplicityVSMEANpTL2 = new TH2D("Multiplicity wo Bins to <pT> #Lambda Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTL2);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Anti Lambda Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiALO = new TH2D("Sphericity w Bins #bar{#Lambda}^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiALO);
  hsphericityALO = new TH2D("Sphericity wo Bins #bar{#Lambda}^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityALO); 

  hsphericityVSpTALO = new TH2D("Sphericity wo Bins to pT #bar{#Lambda}^{0} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTALO);
  hsphericityVSMEANpTALO = new TH2D("Sphericity wo Bins to <pT> #bar{#Lambda}^{0} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTALO);
  
  hsphericityVScentALO = new TH2D("Sphericity wo Bins to Centrality #bar{#Lambda}^{0} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentALO);

  hmultiplicityVSpTbinALO = new TH2D("Multiplicity w Bins to pT #bar{#Lambda}^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinALO);  
  hmultiplicityVSMEANpTbinALO = new TH2D("Multiplicity w Bins to <pT> #bar{#Lambda}^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinALO);
  hmultiplicityVSpTALO = new TH2D("Multiplicity wo Bins to pT #bar{#Lambda}^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTALO);  
  hmultiplicityVSMEANpTALO = new TH2D("Multiplicity wo Bins to <pT> #bar{#Lambda}^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTALO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Lambda Only 2 /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiLOO = new TH2D("Sphericity w Bins #Lambda^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiLOO);
  hsphericityLOO = new TH2D("Sphericity wo Bins #Lambda^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityLOO); 

  hsphericityVSpTLOO = new TH2D("Sphericity wo Bins to pT #Lambda^{0} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTLOO);
  hsphericityVSMEANpTLOO = new TH2D("Sphericity wo Bins to <pT> #Lambda^{0} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTLOO);
  
  hsphericityVScentLOO = new TH2D("Sphericity wo Bins to Centrality #Lambda^{0} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentLOO);

  hmultiplicityVSpTbinLOO = new TH2D("Multiplicity w Bins to pT #Lambda^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinLOO);  
  hmultiplicityVSMEANpTbinLOO = new TH2D("Multiplicity w Bins to <pT> #Lambda^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinLOO);
  hmultiplicityVSpTLOO = new TH2D("Multiplicity wo Bins to pT #Lambda^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTLOO);  
  hmultiplicityVSMEANpTLOO = new TH2D("Multiplicity wo Bins to <pT> #Lambda^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTLOO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Only Lambda or Anti Lambda Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiALLOO = new TH2D("Sphericity w Bins #Lambda^{0} or #bar{#Lambda}^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiALLOO);
  hsphericityALLOO = new TH2D("Sphericity wo Bins #Lambda^{0} or #bar{#Lambda}^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityALLOO); 

  hsphericityVSpTALLOO = new TH2D("Sphericity wo Bins to pT #Lambda^{0} or #bar{#Lambda}^{0} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTALLOO);
  hsphericityVSMEANpTALLOO = new TH2D("Sphericity wo Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTALLOO);
  
  hsphericityVScentALLOO = new TH2D("Sphericity wo Bins to Centrality #Lambda^{0} or #bar{#Lambda}^{0} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentALLOO);

  hmultiplicityVSpTbinALLOO = new TH2D("Multiplicity w Bins to pT #Lambda^{0} or #bar{#Lambda}^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinALLOO);  
  hmultiplicityVSMEANpTbinALLOO = new TH2D("Multiplicity w Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinALLOO);
  hmultiplicityVSpTALLOO = new TH2D("Multiplicity wo Bins to pT #Lambda^{0} or #bar{#Lambda}^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTALLOO);  
  hmultiplicityVSMEANpTALLOO = new TH2D("Multiplicity wo Bins to <pT> #Lambda^{0} or #bar{#Lambda}^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTALLOO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////// Xi Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXO = new TH2D("Sphericity w Bins #Xi Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXO);
  hsphericityXO = new TH2D("Sphericity wo Bins #Xi Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXO); 

  hsphericityVSpTXO = new TH2D("Sphericity wo Bins to pT #Xi Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXO);
  hsphericityVSMEANpTXO = new TH2D("Sphericity wo Bins to <pT> #Xi Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXO);
  
  hsphericityVScentXO = new TH2D("Sphericity wo Bins to Centrality #Xi Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXO);

  hmultiplicityVSpTbinXO = new TH2D("Multiplicity w Bins to pT #Xi Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXO);  
  hmultiplicityVSMEANpTbinXO = new TH2D("Multiplicity w Bins to <pT> #Xi Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXO);
  hmultiplicityVSpTXO = new TH2D("Multiplicity wo Bins to pT #Xi Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXO);  
  hmultiplicityVSMEANpTXO = new TH2D("Multiplicity wo Bins to <pT> #Xi Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Xi Plus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXPO = new TH2D("Sphericity w Bins #bar{#Xi}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXPO);
  hsphericityXPO = new TH2D("Sphericity wo Bins #bar{#Xi}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXPO); 

  hsphericityVSpTXPO = new TH2D("Sphericity wo Bins to pT #bar{#Xi}^{+} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXPO);
  hsphericityVSMEANpTXPO = new TH2D("Sphericity wo Bins to <pT> #bar{#Xi}^{+} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXPO);
  
  hsphericityVScentXPO = new TH2D("Sphericity wo Bins to Centrality #bar{#Xi}^{+} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXPO);

  hmultiplicityVSpTbinXPO = new TH2D("Multiplicity w Bins to pT #bar{#Xi}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXPO);  
  hmultiplicityVSMEANpTbinXPO = new TH2D("Multiplicity w Bins to <pT> #bar{#Xi}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXPO);
  hmultiplicityVSpTXPO = new TH2D("Multiplicity wo Bins to pT #bar{#Xi}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXPO);  
  hmultiplicityVSMEANpTXPO = new TH2D("Multiplicity wo Bins to <pT> #bar{#Xi}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXPO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Xi Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXNO = new TH2D("Sphericity w Bins #Xi^{-} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXNO);
  hsphericityXNO = new TH2D("Sphericity wo Bins #Xi^{-} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXNO); 

  hsphericityVSpTXNO = new TH2D("Sphericity wo Bins to pT #Xi^{-} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXNO);
  hsphericityVSMEANpTXNO = new TH2D("Sphericity wo Bins to <pT> #Xi^{-} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXNO);
  
  hsphericityVScentXNO = new TH2D("Sphericity wo Bins to Centrality #Xi^{-} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXNO);

  hmultiplicityVSpTbinXNO = new TH2D("Multiplicity w Bins to pT #Xi^{-} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXNO);  
  hmultiplicityVSMEANpTbinXNO = new TH2D("Multiplicity w Bins to <pT> #Xi^{-} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXNO);
  hmultiplicityVSpTXNO = new TH2D("Multiplicity wo Bins to pT #Xi^{-} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXNO);  
  hmultiplicityVSMEANpTXNO = new TH2D("Multiplicity wo Bins to <pT> #Xi^{-} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXNO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Xi Plus Only or Xi Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiXNXPO = new TH2D("Sphericity w Bins #Xi^{-} or #bar{#Xi}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiXNXPO);
  hsphericityXNXPO = new TH2D("Sphericity wo Bins #Xi^{-} or #bar{#Xi}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityXNXPO); 

  hsphericityVSpTXNXPO = new TH2D("Sphericity wo Bins to pT #Xi^{-} or #bar{#Xi}^{+} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTXNXPO);
  hsphericityVSMEANpTXNXPO = new TH2D("Sphericity wo Bins to <pT> #Xi^{-} or #bar{#Xi}^{+} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTXNXPO);
  
  hsphericityVScentXNXPO = new TH2D("Sphericity wo Bins to Centrality #Xi^{-} or #bar{#Xi}^{+} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentXNXPO);

  hmultiplicityVSpTbinXNXPO = new TH2D("Multiplicity w Bins to pT #Xi^{-} or #bar{#Xi}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinXNXPO);  
  hmultiplicityVSMEANpTbinXNXPO = new TH2D("Multiplicity w Bins to <pT> #Xi^{-} or #bar{#Xi}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinXNXPO);
  hmultiplicityVSpTXNXPO = new TH2D("Multiplicity wo Bins to pT #Xi^{-} or #bar{#Xi}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTXNXPO);  
  hmultiplicityVSMEANpTXNXPO = new TH2D("Multiplicity wo Bins to <pT> #Xi^{-} or #bar{#Xi}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTXNXPO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Omega Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiOO = new TH2D("Sphericity w Bins #Omega Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiOO);
  hsphericityOO = new TH2D("Sphericity wo Bins #Omega Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityOO); 

  hsphericityVSpTOO = new TH2D("Sphericity wo Bins to pT #Omega Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTOO);
  hsphericityVSMEANpTOO = new TH2D("Sphericity wo Bins to <pT> #Omega Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTOO);
  
  hsphericityVScentOO = new TH2D("Sphericity wo Bins to Centrality #Omega Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentOO);

  hmultiplicityVSpTbinOO = new TH2D("Multiplicity w Bins to pT #Omega Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinOO);  
  hmultiplicityVSMEANpTbinOO = new TH2D("Multiplicity w Bins to <pT> #Omega Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinOO);
  hmultiplicityVSpTOO = new TH2D("Multiplicity wo Bins to pT #Omega Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTOO);  
  hmultiplicityVSMEANpTOO = new TH2D("Multiplicity wo Bins to <pT> #Omega Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTOO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////////// Omega Plus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiOPO = new TH2D("Sphericity w Bins #bar{#Omega}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiOPO);
  hsphericityOPO = new TH2D("Sphericity wo Bins #bar{#Omega}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityOPO); 

  hsphericityVSpTOPO = new TH2D("Sphericity wo Bins to pT #bar{#Omega}^{+} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTOPO);
  hsphericityVSMEANpTOPO = new TH2D("Sphericity wo Bins to <pT> #bar{#Omega}^{+} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTOPO);
  
  hsphericityVScentOPO = new TH2D("Sphericity wo Bins to Centrality #bar{#Omega}^{+} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentOPO);

  hmultiplicityVSpTbinOPO = new TH2D("Multiplicity w Bins to pT #bar{#Omega}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinOPO);  
  hmultiplicityVSMEANpTbinOPO = new TH2D("Multiplicity w Bins to <pT> #bar{#Omega}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinOPO);
  hmultiplicityVSpTOPO = new TH2D("Multiplicity wo Bins to pT #bar{#Omega}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTOPO);  
  hmultiplicityVSMEANpTOPO = new TH2D("Multiplicity wo Bins to <pT> #bar{#Omega}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTOPO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Omega Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiONO = new TH2D("Sphericity w Bins #Omega^{-} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiONO);
  hsphericityONO = new TH2D("Sphericity wo Bins #Omega^{-} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityONO); 

  hsphericityVSpTONO = new TH2D("Sphericity wo Bins to pT #Omega^{-} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTONO);
  hsphericityVSMEANpTONO = new TH2D("Sphericity wo Bins to <pT> #Omega^{-} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTONO);
  
  hsphericityVScentONO = new TH2D("Sphericity wo Bins to Centrality #Omega^{-} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentONO);

  hmultiplicityVSpTbinONO = new TH2D("Multiplicity w Bins to pT #Omega^{-} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinONO);  
  hmultiplicityVSMEANpTbinONO = new TH2D("Multiplicity w Bins to <pT> #Omega^{-} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinONO);
  hmultiplicityVSpTONO = new TH2D("Multiplicity wo Bins to pT #Omega^{-} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTONO);  
  hmultiplicityVSMEANpTONO = new TH2D("Multiplicity wo Bins to <pT> #Omega^{-} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTONO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////// Omega Plus Only or Omega Minus Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiONOPO = new TH2D("Sphericity w Bins #Omega^{-} or #bar{#Omega}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiONOPO);
  hsphericityONOPO = new TH2D("Sphericity wo Bins #Omega^{-} or #bar{#Omega}^{+} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityONOPO); 

  hsphericityVSpTONOPO = new TH2D("Sphericity wo Bins to pT #Omega^{-} or #bar{#Omega}^{+} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTONOPO);
  hsphericityVSMEANpTONOPO = new TH2D("Sphericity wo Bins to <pT> #Omega^{-} or #bar{#Omega}^{+} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTONOPO);
  
  hsphericityVScentONOPO = new TH2D("Sphericity wo Bins to Centrality #Omega^{-} or #bar{#Omega}^{+} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentONOPO);

  hmultiplicityVSpTbinONOPO = new TH2D("Multiplicity w Bins to pT #Omega^{-} or #bar{#Omega}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinONOPO);  
  hmultiplicityVSMEANpTbinONOPO = new TH2D("Multiplicity w Bins to <pT> #Omega^{-} or #bar{#Omega}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinONOPO);
  hmultiplicityVSpTONOPO = new TH2D("Multiplicity wo Bins to pT #Omega^{-} or #bar{#Omega}^{+} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTONOPO);  
  hmultiplicityVSMEANpTONOPO = new TH2D("Multiplicity wo Bins to <pT> #Omega^{-} or #bar{#Omega}^{+} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTONOPO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////// Kaon Only /////////////////////////////////////////////////////////////////////////////
  
  hsphericityVSmultiKAO = new TH2D("Sphericity w Bins #Kappa^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}", nMultbins, Multbins, nBinsSt, BinSt);              
  fListOfObjects->Add(hsphericityVSmultiKAO);
  hsphericityKAO = new TH2D("Sphericity wo Bins #Kappa^{0} Only"," Sphericity vs Multiplicity; Multiplicity #it{N}_{CH}; Sphericity #it{S}_{T}",1600, 0., 1600., 100, 0., 1.);     
  fListOfObjects->Add(hsphericityKAO); 

  hsphericityVSpTKAO = new TH2D("Sphericity wo Bins to pT #Kappa^{0} Only"," Sphericity vs #it{p}_{T};Sphericity #it{S}_{T}; #it{p}_{T} (GeV/#it{c})",100, 0.0, 1., 500, 0., 5.);                   
  fListOfObjects->Add(hsphericityVSpTKAO);
  hsphericityVSMEANpTKAO = new TH2D("Sphericity wo Bins to <pT> #Kappa^{0} Only"," Sphericity vs <#it{p}_{T}>;Sphericity #it{S}_{T}; <#it{p}_{T}> (GeV/#it{c})",100, 0.0, 1., 300, 0.3, 0.9);                   
  fListOfObjects->Add(hsphericityVSMEANpTKAO);
  
  hsphericityVScentKAO = new TH2D("Sphericity wo Bins to Centrality #Kappa^{0} Only"," Sphericity vs Centrality;Sphericity #it{S}_{T}; Centrality [#it{%}]",100, 0.0, 1., 100, 0., 90.);                 
  fListOfObjects->Add(hsphericityVScentKAO);

  hmultiplicityVSpTbinKAO = new TH2D("Multiplicity w Bins to pT #Kappa^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",nMultbins, Multbins, 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTbinKAO);  
  hmultiplicityVSMEANpTbinKAO = new TH2D("Multiplicity w Bins to <pT> #Kappa^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",nMultbins, Multbins, 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTbinKAO);
  hmultiplicityVSpTKAO = new TH2D("Multiplicity wo Bins to pT #Kappa^{0} Only"," Multiplicity vs #it{p}_{T};Multiplicity #it{N}_{CH}; #it{p}_{T} (GeV/#it{c})",1600, 0.0, 1600., 500, 0., 5.);              
  fListOfObjects->Add(hmultiplicityVSpTKAO);  
  hmultiplicityVSMEANpTKAO = new TH2D("Multiplicity wo Bins to <pT> #Kappa^{0} Only"," Multiplicity vs <#it{p}_{T}>;Multiplicity #it{N}_{CH}; <#it{p}_{T}> (GeV/#it{c})",1600, 0.0, 1600., 300, 0.3, 0.9);   
  fListOfObjects->Add(hmultiplicityVSMEANpTKAO);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Post output data.
  PostData(1, fListOfObjects);
  
}

//______________________________________________________________________________
void AliAnalysisSphericityTask::UserExec(Option_t *) 
{
  // Main loop
    
  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }
  AliESDEvent* lESDEvent =(AliESDEvent*)event;

  hn1->Fill(0);
  Bool_t isPileup = kFALSE;

  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);
    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    isPileup = fESD->IsPileupFromSPD();
    if(fPileUpRej)
      if(isPileup)
	      return;
  } else {
    fAOD = dynamic_cast<AliAODEvent*>(event);
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
    isPileup = fAOD->IsPileupFromSPD();
    if(fPileUpRej)
      if(isPileup)
	      return;    
  }
  
  if (fAnalysisMC) {
    if (fAnalysisType == "ESD"){
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(!fMC){
	      Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	      this->Dump();
	      return;
      }    
    } else { // AOD
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if(fMC)
	      fMC->Dump();
      fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
      if(!fMCArray){
	      Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
	      this->Dump();
	      return;
      }    
    }
  }
  
  hn1->Fill(1);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Get trigger decision
  fTriggeredEventMB = 0; //init
    
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigBit){
    fTriggeredEventMB = 1;  //event triggered as minimum bias
  }
  
  // real data that are not triggered we skip
  if(!fAnalysisMC && !fTriggeredEventMB)
    return; 
  hn1->Fill(2);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Double_t tPrimaryVtxPosition[3];
  Double_t lMagneticField = event->GetMagneticField();
  Double_t lPt = 0;
  Double_t lRapLambda = 0;
  Double_t lInvMassLambda = 0, lInvMassAntiLambda = 0, lInvMassK0s = 0;
  Double_t lAlphaV0 = 0;
  Double_t lPtArmV0 = 0;
  Bool_t   lOnFlyStatus;
  Double_t lPz = 0;
  Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  Double_t lDcaPosToPrimVertexV0 = 0, lDcaNegToPrimVertexV0 = 0;
  Double_t lV0CosineOfPointingAngle = 0;
  Double_t lV0Radius = 0;
  Double_t tDecayVertexV0[3];
  Int_t    ncascades = 0;
           ncascades = event->GetNumberOfCascades();
  Double_t lInvMassXiMinus = 0.;
  Double_t lInvMassXiPlus = 0.;
  Double_t lInvMassOmegaMinus = 0.;
  Double_t lInvMassOmegaPlus = 0.;
  Short_t  lChargeXi = -2;
  Double_t lRapXi   = -20.0, lRapOmega = -20.0;
  Double_t lXiTransvMom  = 0.;
  Double_t lAlphaXi = 0;
  Double_t lPtArmXi = 0;
  Double_t lXiMomX = 0. , lXiMomY = 0., lXiMomZ = 0.;
  Double_t lV0quality = 0.;
  Double_t lInvMassLambdaAsCascDghter = -1;
  Double_t lDcaXiDaughters = -1. ;
  Double_t lXiCosineOfPointingAngle = -1. ;
  Double_t lXiRadius2D = -1000. ;
  Double_t lDcaV0DaughtersXi = -1.;		
  Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
  Double_t lDcaPosToPrimVertexXi  = -1.;
  Double_t lDcaNegToPrimVertexXi  = -1.;
  Double_t lV0CosineOfPointingAngleXi = -1. ;
  Double_t lV0RadiusXi = -1000.0;	
  Double_t lPosXi[3] = {-1000.0, -1000.0, -1000.0};
  Double_t lPosV0Xi[3] = {-1000.0, -1000.0, -1000.0};
  
  fStrangePart=kFALSE;
  fFindKaon=kFALSE;
  fFindLambda=kFALSE;
  fFindAntiLambda=kFALSE;
  fFindXiPlus=kFALSE;
  fFindXiMinus=kFALSE;
  fFindOmegaPlus=kFALSE;
  fFindOmegaMinus=kFALSE;
  	
  if (fAnalysisType == "ESD"){
    const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
    if(vtxESD->GetNContributors()<1) {
      // SPD vertex
      vtxESD = fESD->GetPrimaryVertexSPD();
      /* quality checks on SPD-vertex */
      if (vtxESD->IsFromVertexerZ() && (vtxESD->GetDispersion() > 0.04 || vtxESD->GetZRes() > 0.25))  
	      fZvtx  = -1599; //vertex = 0x0; //
      else if (vtxESD->GetNContributors()<1) 
	      fZvtx  = -999; //vertex = 0x0; //
      else {
        const AliESDVertex *trkVtx = fESD->GetPrimaryVertex();
        if (!trkVtx || trkVtx->GetNContributors()<=0) return;
        TString vtxTtl = trkVtx->GetTitle();
        if (!vtxTtl.Contains("VertexerTracks")) return;
        Float_t zvtx = trkVtx->GetZ();
        const AliESDVertex* spdVtx = fESD->GetPrimaryVertexSPD();
        if (spdVtx->GetNContributors()<=0) return;
        TString vtxTyp = spdVtx->GetTitle();  
        Double_t cov[6]={0};
        spdVtx->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
        if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
	      fZvtx = vtxESD->GetZ();
        }
    }  
    else{  
      const AliESDVertex *trkVtx = fESD->GetPrimaryVertex();
      if (!trkVtx || trkVtx->GetNContributors()<=0) return;
      TString vtxTtl = trkVtx->GetTitle();
      if (!vtxTtl.Contains("VertexerTracks")) return;
      Float_t zvtx = trkVtx->GetZ();
      const AliESDVertex* spdVtx = fESD->GetPrimaryVertexSPD();
      if (spdVtx->GetNContributors()<=0) return;
      TString vtxTyp = spdVtx->GetTitle();  
      Double_t cov[6]={0};
      spdVtx->GetCovarianceMatrix(cov);
      Double_t zRes = TMath::Sqrt(cov[5]);
      if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
      if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
      fZvtx = vtxESD->GetZ();
    }
  }
  else // AOD
    fZvtx = GetVertex(fAOD);
  
  hVtxBeforeCuts->Fill(fZvtx);
  
  //cut on the z position of vertex
  if (TMath::Abs(fZvtx) > fVtxCut) {	
    return;
  }

  const AliVVertex *primaryVtx = event->GetPrimaryVertex();
  tPrimaryVtxPosition[0] = primaryVtx->GetX();
  tPrimaryVtxPosition[1] = primaryVtx->GetY();
  tPrimaryVtxPosition[2] = primaryVtx->GetZ();

  hn1->Fill(3);

  //cut on centrality
  Float_t CentPercentile = 0;
  Int_t ntracklet = 0;
  if(fAnalysisType == "ESD"){
    AliCentrality *centrality = fESD->GetCentrality();
    CentPercentile = centrality->GetCentralityPercentile(fCentEst);
    if (!(fMinCent < CentPercentile && CentPercentile < fMaxCent))	
      return;
  }
  hn1->Fill(4);                                                  //fills to centrality cut

  Int_t     IndxTrksMult = -1;
	Int_t     IndxV0MMult = -1;
  Double_t  totaleta = 0;
  Double_t  totalpt =  0;
  Double_t  totalphi = 0;
  Int_t     TrackMult03=0;   
  Double_t  spherocity=0;
  Double_t  sphericity=0;
  Int_t     Event_Type=0;
  Int_t     nTracks = fESD->GetNumberOfTracks();
  
  if(fTriggeredEventMB) {    // only analyze triggered events

    TrackMult03=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);    
    HMultRef->Fill(TrackMult03);

    ntracklet = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8); 
    fHistTrackMultiplicity->Fill(ntracklet);        

    hcent->Fill(CentPercentile);

    if(fStrangeness)    //analyze only event containing stragne partitions
    {  
      for (Int_t iXi = 0; iXi < ncascades; iXi++){            // This is the begining of the Cascade loop
      	AliESDcascade *xi = lESDEvent->GetCascade(iXi);
	      if (!xi) continue;

	    lChargeXi = xi->Charge();
	    lRapXi    = xi->RapXi();
	    lRapOmega = xi->RapOmega();
      lAlphaXi  = xi->AlphaXi();
	    lPtArmXi  = xi->PtArmXi();
	    xi->GetPxPyPz(lXiMomX, lXiMomY, lXiMomZ);
	    lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX + lXiMomY*lXiMomY);
      lDcaXiDaughters	= xi->GetDcaXiDaughters();
      lDcaV0DaughtersXi = xi->GetDcaV0Daughters(); 
	    lXiCosineOfPointingAngle = xi->GetCascadeCosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	    lV0CosineOfPointingAngleXi = xi->GetV0CosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	    lDcaV0ToPrimVertexXi = xi->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	 						       						     
	    UInt_t lIdxPosXi = (UInt_t) TMath::Abs( xi->GetPindex());
      UInt_t lIdxNegXi = (UInt_t) TMath::Abs( xi->GetNindex());
      UInt_t lBachIdx = (UInt_t) TMath::Abs( xi->GetBindex());
      if(lBachIdx == lIdxNegXi) continue; 
      if(lBachIdx == lIdxPosXi) continue;

	    AliESDtrack *pTrackXi = ((AliESDEvent*)event)->GetTrack( lIdxPosXi );
	    AliESDtrack *nTrackXi = ((AliESDEvent*)event)->GetTrack( lIdxNegXi );
	    AliESDtrack *bachTrackXi = ((AliESDEvent*)event)->GetTrack( lBachIdx );
                if (!pTrackXi || !nTrackXi || !bachTrackXi )  continue;

	    const AliExternalTrackParam *pExtTrack = pTrackXi->GetTPCInnerParam();
      const AliExternalTrackParam *nExtTrack = nTrackXi->GetTPCInnerParam();
	    const AliExternalTrackParam *bachExtTrack = bachTrackXi->GetTPCInnerParam();
		  if (!pExtTrack || !nExtTrack || !bachExtTrack )  continue;

      // Daughter track quality cuts
        if (!IsAccepted(pTrackXi))    continue;
	      if (!IsAccepted(nTrackXi))    continue;
	      if (!IsAccepted(bachTrackXi)) continue;

      if( bachTrackXi->Charge() < 0 )	{
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,3312); 	
		  lInvMassXiMinus = xi->GetEffMassXi();                                  
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,3334);	
		  lInvMassOmegaMinus = xi->GetEffMassXi();	
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,3312); 
	    }
	    if( bachTrackXi->Charge() >  0 ){
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,-3312);	
		  lInvMassXiPlus = xi->GetEffMassXi();
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,-3334);	
		  lInvMassOmegaPlus = xi->GetEffMassXi();
		  lV0quality = 0.;
		  xi->ChangeMassHypothesis(lV0quality,-3312);
	    }
      
	    lInvMassLambdaAsCascDghter = xi->GetEffMass();
	    lDcaBachToPrimVertexXi = TMath::Abs(bachTrackXi->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField));
	    lDcaPosToPrimVertexXi = TMath::Abs(pTrackXi->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField));
	    lDcaNegToPrimVertexXi = TMath::Abs(nTrackXi->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField));
	    xi->GetXYZcascade(lPosXi[0],lPosXi[1],lPosXi[2]); 
	    lXiRadius2D = TMath::Sqrt(lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1]);
	    xi->GetXYZ(lPosV0Xi[0],lPosV0Xi[1],lPosV0Xi[2]); 
	    lV0RadiusXi = TMath::Sqrt(lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1]);
 	
	    Double_t lXiTotalMomentum = TMath::Sqrt(lXiMomX*lXiMomX+lXiMomY*lXiMomY+lXiMomZ*lXiMomZ);
	    Double_t lCtauXi = 1.321 * TMath::Sqrt(
			      TMath::Power( lPosXi[0] - tPrimaryVtxPosition[0] ,2) +
			      TMath::Power( lPosXi[1] - tPrimaryVtxPosition[1] ,2) +
			      TMath::Power( lPosXi[2] - tPrimaryVtxPosition[2] ,2)
			      );
      	  lCtauXi /= (lXiTotalMomentum+1e-10);    //avoid division by zero, to be sure

	    Double_t lCtauOmega = 1.672 * TMath::Sqrt(
			      TMath::Power( lPosXi[0] - tPrimaryVtxPosition[0] ,2) +
			      TMath::Power( lPosXi[1] - tPrimaryVtxPosition[1] ,2) +
			      TMath::Power( lPosXi[2] - tPrimaryVtxPosition[2] ,2)
			      );
      	  lCtauOmega /= (lXiTotalMomentum+1e-10); //avoid division by zero, to be sure
      	       
      if((IsProton(pTrackXi) && IsPion(nTrackXi)) ||( IsProton(nTrackXi) && IsPion(pTrackXi))){
		
	      if(fkExtraSelections){
	        if(lDcaPosToPrimVertexXi < 0.1)continue;
	        if(lDcaNegToPrimVertexXi < 0.1)continue;
	        if(lDcaV0DaughtersXi > 0.8)continue;
	        if(lV0CosineOfPointingAngleXi < 0.998)continue;
	        if(lV0RadiusXi < 3.0)continue;
	        if(lDcaV0ToPrimVertexXi < 0.1)continue;
	        if(TMath::Abs(lInvMassLambdaAsCascDghter - 1.11568) > 0.005)continue;
	        if(lDcaBachToPrimVertexXi < 0.03)continue;
	        if(lDcaXiDaughters > 0.3)continue;
	        if(lXiCosineOfPointingAngle < 0.9992)continue;
	        if(lXiRadius2D < 1.0)continue;
	      }
             
	      if( lChargeXi < 0 ){
		      if(IsPion(bachTrackXi) && TMath::Abs(lRapXi) < 0.6 && lCtauXi < 15)
          {
            if(fkExtraSelectionsCut){
				      if (TMath::Abs(lInvMassXiMinus -1.3217) < fInvMassCutXi) 
                {
				        fHistMassXiMinus->Fill(lInvMassXiMinus);	
				        fFindXiMinus=kTRUE;
                fStrangePart=kTRUE;
				        fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
                }
            }
            else
            {
              fHistMassXiMinus->Fill(lInvMassXiMinus);	
				      fFindXiMinus=kTRUE;
              fStrangePart=kTRUE;
				      fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
            }
		      }
		      if(IsKaon(bachTrackXi) && TMath::Abs(lRapOmega) < 0.6 && lCtauOmega < 8 && (TMath::Abs( lInvMassXiMinus -1.3217 ) > 0.010)){
            if(fkExtraSelectionsCut){
		  		    if (TMath::Abs(lInvMassOmegaMinus -1.67245) < fInvMassCutOmega) 
                {
				        fHistMassOmegaMinus->Fill(lInvMassOmegaMinus);	
				        fFindOmegaMinus=kTRUE;
                fStrangePart=kTRUE;
				        fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
                }
            }
            else
            {
              fHistMassOmegaMinus->Fill(lInvMassOmegaMinus);	
				      fFindOmegaMinus=kTRUE;
              fStrangePart=kTRUE;
				      fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
            }
		      }
	      }
	      if( lChargeXi > 0 ){
		      if(IsPion(bachTrackXi) && TMath::Abs(lRapXi) < 0.6 && lCtauXi < 15){
            if(fkExtraSelectionsCut){
		  		    if (TMath::Abs(lInvMassXiPlus -1.3217) < fInvMassCutXi) 
                {
			  	      fHistMassXiPlus->Fill(lInvMassXiPlus);	
			  	      fFindXiPlus=kTRUE;
			  	      fStrangePart=kTRUE;
                fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
                }
            }
            else{
              fHistMassXiPlus->Fill(lInvMassXiPlus);	
			  	    fFindXiPlus=kTRUE;
			  	    fStrangePart=kTRUE;
              fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
            }
          
          }
		      if(IsKaon(bachTrackXi) && TMath::Abs(lRapOmega) < 0.6 && lCtauOmega < 8 && (TMath::Abs( lInvMassXiPlus -1.3217 ) > 0.010)){
            if(fkExtraSelectionsCut){
		  		    if (TMath::Abs(lInvMassOmegaPlus -1.67245) < fInvMassCutOmega) 
                {
			  	      fHistMassOmegaPlus->Fill(lInvMassOmegaPlus);	
                fFindOmegaPlus=kTRUE;
			  	      fStrangePart=kTRUE;
			  	      fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
                }
            }
            else{
              fHistMassOmegaPlus->Fill(lInvMassOmegaPlus);	
              fFindOmegaPlus=kTRUE;
			  	    fStrangePart=kTRUE;
			  	    fHistXiArmenteros->Fill(lChargeXi*lAlphaXi,lPtArmXi);
            }
		      }
	      }
 	    }
    }//end of Cascade loop
    
    if(fRerunVertexers){	
	        lESDEvent->ResetV0s();
		      AliV0vertexer lV0vtxer;
		      lV0vtxer.SetCuts(fV0Cuts);
	        lV0vtxer.Tracks2V0vertices(lESDEvent);
    }

    Int_t nv0s = 0;
    nv0s = event->GetNumberOfV0s();

    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {// This is the begining of the V0 loop
	    AliESDv0 *v0 = lESDEvent->GetV0(iV0);
	    if (!v0) continue;
          
	    lOnFlyStatus = v0->GetOnFlyStatus();
	    if (lOnFlyStatus) continue;
	    lPz = v0->Pz();
	    lPt = v0->Pt();	  
	    lRapLambda  = v0->RapLambda();
	    lAlphaV0 = v0->AlphaV0();
	    lPtArmV0 = v0->PtArmV0();

	    if (TMath::Abs(lRapLambda) > fMaxV0Rapidity) continue;
	  
	    UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
	    UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

	    AliESDtrack *pTrack=((AliESDEvent*)event)->GetTrack(lKeyPos);
	    AliESDtrack *nTrack=((AliESDEvent*)event)->GetTrack(lKeyNeg);
	    if (!pTrack || !nTrack) {
	      Printf("ERROR: Could not retreive one of the daughter track");
	    continue;
	    }
	    // Filter like-sign V0 (next: add counter and distribution)
	    if ( pTrack->GetSign() == nTrack->GetSign()){
	      continue;
	    }
	  
	    AliExternalTrackParam *pExtTrack = (AliExternalTrackParam *)pTrack->GetTPCInnerParam();
	    AliExternalTrackParam *nExtTrack = (AliExternalTrackParam *)nTrack->GetTPCInnerParam();
		  if (!pExtTrack || !nExtTrack)  continue;

      // Tracks quality cuts 
	    if (!IsAccepted(pTrack)) continue;
	    if (!IsAccepted(nTrack)) continue;

      // Getting invariant mass infos directly from ESD
	    v0->ChangeMassHypothesis(310);
	    lInvMassK0s = v0->GetEffMass();
	    v0->ChangeMassHypothesis(3122);
	    lInvMassLambda = v0->GetEffMass();
	    v0->ChangeMassHypothesis(-3122);
	    lInvMassAntiLambda = v0->GetEffMass();
      
	    v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
	    lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
	    lDcaPosToPrimVertexV0 = TMath::Abs(pTrack->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField));
	    lDcaNegToPrimVertexV0 = TMath::Abs(nTrack->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField));
	    lDcaV0Daughters = v0->GetDcaV0Daughters();	  
	    lDcaV0ToPrimVertex = v0->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	    lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	  
	    Double_t tV0mom[3];
          v0->GetPxPyPz(tV0mom[0],tV0mom[1],tV0mom[2]);
	    Double_t lV0TotalMomentum = TMath::Sqrt(tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2]);
	    Double_t lCtau = 1.115683 * TMath::Sqrt(
			      TMath::Power( tDecayVertexV0[0] - tPrimaryVtxPosition[0] ,2) +
			      TMath::Power( tDecayVertexV0[1] - tPrimaryVtxPosition[1] ,2) +
			      TMath::Power( tDecayVertexV0[2] - tPrimaryVtxPosition[2] ,2)
			      );
      	  lCtau /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

	    if(fkExtraSelections){
	  	  if(lDcaPosToPrimVertexV0 < 0.1) continue;
	  	  if(lDcaNegToPrimVertexV0 < 0.1) continue;
		    if(lDcaV0Daughters > 1.0) continue;
		    if(lV0CosineOfPointingAngle < 0.998) continue;
		    if(lV0Radius < 0.9) continue;
	    }
                
      if(IsKaon (pTrack))
        {
           if(fkExtraSelectionsCut)
            {
			      if (TMath::Abs(lInvMassK0s -0.497611) < fInvMassCutKaon ) {
              fStrangePart=kTRUE;
			        fHistMassKaon->Fill(lInvMassK0s);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindKaon=kTRUE;
              }
            }
            else 
            {
              fStrangePart=kTRUE;
			        fHistMassKaon->Fill(lInvMassK0s);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindKaon=kTRUE;
            }
			  }
      if(IsKaon (nTrack))
      {
           if(fkExtraSelectionsCut)
            {
			      if (TMath::Abs(lInvMassK0s -0.497611) < fInvMassCutKaon) {
              fStrangePart=kTRUE;
			        fHistMassKaon->Fill(lInvMassK0s);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindKaon=kTRUE;
				      }
            }
            else
            {
              fStrangePart=kTRUE;
			        fHistMassKaon->Fill(lInvMassK0s);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindKaon=kTRUE;
            }
      }
        
	    if((IsProton (pTrack) && IsPion (nTrack)) || (IsProton (nTrack) && IsPion (pTrack))) 
	      {
			    if ( lAlphaV0 > 0){
           if(fkExtraSelectionsCut)
            {
			      if (TMath::Abs(lInvMassLambda -1.11568) < fInvMassCutLambda ) {
              fStrangePart=kTRUE;
			        fHistMassLambda->Fill(lInvMassLambda);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindLambda=kTRUE;
              }
            }
            else 
            {
              fStrangePart=kTRUE;
			        fHistMassLambda->Fill(lInvMassLambda);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindLambda=kTRUE;
            }
			    }
	        if ( lAlphaV0 < 0){
           if(fkExtraSelectionsCut)
            {
			      if (TMath::Abs(lInvMassAntiLambda -1.11568) < fInvMassCutLambda) {
              fStrangePart=kTRUE;
			        fHistMassAntiLambda->Fill(lInvMassAntiLambda);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindAntiLambda=kTRUE;
				      }
            }
            else
            {
              fStrangePart=kTRUE;
			        fHistMassAntiLambda->Fill(lInvMassAntiLambda);	
			        fHistV0Armenteros->Fill(lAlphaV0,lPtArmV0);
              fFindAntiLambda=kTRUE;
            }
			    }
        }
      } // This is the end of the V0 loop
    if(!fStrangePart) return;
    } // end of strangeness analysis

    ///////////////////////////////////////// Sphericity part //////////////////////////////

	  spherocity = fESASelection->GetEventShape( event, "SO", kTRUE );
    sphericity = fESASelection->GetEventShape( event, "ST", kTRUE );
        
    cout<<"not cuts: So="<<spherocity<<"\t St = "<< sphericity <<endl;
    
    // Events with non-measured spherocity
    if( spherocity < 0 )
          return;
      
    // Jetty-like Events
    if( spherocity < 0.3 )
      {
        Event_Type=0;
      }

    // Mid-Sphericity Events
    if( 0.3 < spherocity && spherocity < 0.7 )
      {
        Event_Type=2;
      }  

    // Isotropic Events
    if( spherocity > 0.7 )
      {
        Event_Type=1;
      } 

    ///////////////////////////////////////////// event statistics ///////////////////////////////////
    
    Int_t   nRec = 0;

	  for(Int_t i1 = 0; i1 < nTracks; ++i1) // event statistics loop
    {
      AliESDtrack* Track = 0;
		  Track = fESD->GetTrack(i1);
      
      if(!Track)
			continue;
      hn2->Fill(0);
      
      Float_t eta  = Track->Eta();
		  Float_t pt   = Track->Pt();
		  Float_t phi  = Track->Phi();

      if(fIsAbsEtaESA){  //cuts in pseudorapidity
			if( TMath::Abs(eta) > fEtaMaxCutESA || TMath::Abs(eta) < fEtaMinCutESA )          
				continue;
		  }
		  else{
			if( eta > fEtaMaxCutESA || eta < fEtaMinCutESA )
				continue;
		  }
      hn2->Fill(1);
		  //cuts in pt
		  if( pt > fPtMaxCutESA || pt <  fPtMinCutESA )
			continue;
      hn2->Fill(2);
		  //quality cuts
		  if(!fUseHybrid){//golden track cuts
			if(!fTrackFilterESA->IsSelected(Track))
				continue;
		  }
		  else{//hybrid track cuts
			Bool_t cutset1 = kFALSE;
			Bool_t cutset2 = kFALSE;
			cutset1 = fTrackFilterHybrid1->IsSelected(Track);
			cutset2 = (!fTrackFilterESA->IsSelected(Track)) && fTrackFilterHybrid2->IsSelected(Track);
			if(!(cutset1 || cutset2))
				continue;
		  }

      nRec ++;
      totaleta += eta;
      totalpt += pt;
      totalphi += phi;
      fhetaSt->Fill(eta);
			fhphiSt->Fill(phi);
			fhptSt->Fill(pt);
      hetaVSphi->Fill(eta,phi);
      hn2->Fill(3);                              
      
      if(Event_Type == 1)
        hetaVSphiISO->Fill(eta,phi);
      
      if(Event_Type == 0)
        hetaVSphiJET->Fill(eta,phi);

      if(Event_Type == 2)
        hetaVSphiMID->Fill(eta,phi);

      ///////////////////////////////////////////
      hsphericityVSpT->Fill(sphericity,pt);            
      hmultiplicityVSpT->Fill(TrackMult03,pt); 
      hmultiplicityVSpTbin->Fill(TrackMult03,pt);
      ///////////////////////////////////////////

      if(fFindAntiLambda && fFindLambda)  //Events with Lambdas
        {
          hsphericityVSpTL->Fill(sphericity,pt);
          hmultiplicityVSpTL->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinL->Fill(TrackMult03,pt);
        }  
      if(fFindAntiLambda)  //Events with Anti Lambdas
        {
          hsphericityVSpTAL->Fill(sphericity,pt);
          hmultiplicityVSpTAL->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinAL->Fill(TrackMult03,pt);
        }
      if(fFindLambda)  //Events with Lambdas Only
        {
          hsphericityVSpTLO->Fill(sphericity,pt);
          hmultiplicityVSpTLO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinLO->Fill(TrackMult03,pt);
        }
      if(fFindLambda || fFindAntiLambda)  //Events with Lambdas or Anti Lambdas Only
        {
          hsphericityVSpTALLO->Fill(sphericity,pt);
          hmultiplicityVSpTALLO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinALLO->Fill(TrackMult03,pt);
        }
      if(fFindXiPlus && fFindXiMinus)  //Events with Xis
        {
          hsphericityVSpTX->Fill(sphericity,pt);
          hmultiplicityVSpTX->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinX->Fill(TrackMult03,pt);
        }  
      if(fFindXiPlus)  //Events with Xis+
        {
          hsphericityVSpTXP->Fill(sphericity,pt);
          hmultiplicityVSpTXP->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXP->Fill(TrackMult03,pt);
        }
      if(fFindXiMinus)  //Events with Xis- Only
        {
          hsphericityVSpTXN->Fill(sphericity,pt);
          hmultiplicityVSpTXN->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXN->Fill(TrackMult03,pt);
        }
      if(fFindXiPlus || fFindXiMinus)  //Events with Xis+ or Xis- Only
        {
          hsphericityVSpTXNXP->Fill(sphericity,pt);
          hmultiplicityVSpTXNXP->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXNXP->Fill(TrackMult03,pt);
        }
        if(fFindOmegaPlus && fFindOmegaMinus)  //Events with Omegas
        {
          hsphericityVSpTO->Fill(sphericity,pt);
          hmultiplicityVSpTO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinO->Fill(TrackMult03,pt);
        }  
      if(fFindOmegaPlus)  //Events with Omegas+
        {
          hsphericityVSpTOP->Fill(sphericity,pt);
          hmultiplicityVSpTOP->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinOP->Fill(TrackMult03,pt);
        }
      if(fFindOmegaMinus)  //Events with Omegas- Only
        {
          hsphericityVSpTON->Fill(sphericity,pt);
          hmultiplicityVSpTON->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinON->Fill(TrackMult03,pt);
        }
      if(fFindOmegaPlus || fFindOmegaMinus)  //Events with Omegas+ or Omegas- Only
        {
          hsphericityVSpTONOP->Fill(sphericity,pt);
          hmultiplicityVSpTONOP->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinONOP->Fill(TrackMult03,pt);
        }
      if(fFindKaon)  //Events with Kaons
        {
          hsphericityVSpTKA->Fill(sphericity,pt);
          hmultiplicityVSpTKA->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinKA->Fill(TrackMult03,pt);
        }
      //////////////////////////////////////////////////////////////////////////////////
      
      if((fFindAntiLambda && fFindLambda) && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Lambdas
        {
          hsphericityVSpTL2->Fill(sphericity,pt);
          hmultiplicityVSpTL2->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinL2->Fill(TrackMult03,pt);
        }  
      if(fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Anti Lambdas
        {
          hsphericityVSpTALO->Fill(sphericity,pt);
          hmultiplicityVSpTALO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinALO->Fill(TrackMult03,pt);
        }
      if(fFindLambda && !fFindAntiLambda && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Lambdas Only
        {
          hsphericityVSpTLOO->Fill(sphericity,pt);
          hmultiplicityVSpTLOO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinLOO->Fill(TrackMult03,pt);
        }
      if((fFindLambda || fFindAntiLambda) && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Lambdas or Anti Lambdas Only
        {
          hsphericityVSpTALLOO->Fill(sphericity,pt);
          hmultiplicityVSpTALLOO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinALLOO->Fill(TrackMult03,pt);
        }
      if((fFindXiPlus && fFindXiMinus) && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindAntiLambda && !fFindLambda)  //Events with Xis
        {
          hsphericityVSpTXO->Fill(sphericity,pt);
          hmultiplicityVSpTXO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXO->Fill(TrackMult03,pt);
        }  
      if(fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindAntiLambda && !fFindLambda)  //Events with Xis+
        {
          hsphericityVSpTXPO->Fill(sphericity,pt);
          hmultiplicityVSpTXPO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXPO->Fill(TrackMult03,pt);
        }
      if(fFindXiMinus && !fFindXiPlus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindAntiLambda && !fFindLambda)  //Events with Xis- Only
        {
          hsphericityVSpTXNO->Fill(sphericity,pt);
          hmultiplicityVSpTXNO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXNO->Fill(TrackMult03,pt);
        }
      if((fFindXiPlus || fFindXiMinus) && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindAntiLambda && !fFindLambda)  //Events with Xis+ or Xis- Only
        {
          hsphericityVSpTXNXPO->Fill(sphericity,pt);
          hmultiplicityVSpTXNXPO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinXNXPO->Fill(TrackMult03,pt);
        }
        if((fFindOmegaPlus && fFindOmegaMinus) && !fFindKaon && !fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus)  //Events with Omegas
        {
          hsphericityVSpTOO->Fill(sphericity,pt);
          hmultiplicityVSpTOO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinOO->Fill(TrackMult03,pt);
        }  
      if(fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus)  //Events with Omegas+
        {
          hsphericityVSpTOPO->Fill(sphericity,pt);
          hmultiplicityVSpTOPO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinOPO->Fill(TrackMult03,pt);
        }
      if(fFindOmegaMinus && !fFindOmegaPlus && !fFindKaon && !fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus)  //Events with Omegas- Only
        {
          hsphericityVSpTONO->Fill(sphericity,pt);
          hmultiplicityVSpTONO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinONO->Fill(TrackMult03,pt);
        }
      if((fFindOmegaPlus || fFindOmegaMinus) && !fFindKaon && !fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus)  //Events with Omegas+ or Omegas- Only
        {
          hsphericityVSpTONOPO->Fill(sphericity,pt);
          hmultiplicityVSpTONOPO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinONOPO->Fill(TrackMult03,pt);
        }
      if(fFindKaon && !fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus && !fFindOmegaMinus && !fFindOmegaPlus)  //Events with Kaons
        {
          hsphericityVSpTKAO->Fill(sphericity,pt);
          hmultiplicityVSpTKAO->Fill(TrackMult03,pt); 
          hmultiplicityVSpTbinKAO->Fill(TrackMult03,pt);
        }
  	}   // end of event statistics loop
    
    totalpt = totalpt/nRec;     //mean pT
    
    /////////////////////////////////////////////////////////////////////////////////
    
    hn1->Fill(5);
    hso->Fill(spherocity);
    hsphericityVSmulti->Fill(TrackMult03,sphericity);         
    hsphericityVSMEANpT->Fill(sphericity,totalpt);     
    hmultiplicityVSMEANpT->Fill(TrackMult03,totalpt); 
    hmultiplicityVSMEANpTbin->Fill(TrackMult03,totalpt);     
    hsphericity->Fill(TrackMult03,sphericity);             
    hst->Fill(sphericity);                        
    hsphericityVScent->Fill(sphericity,CentPercentile);      
    hVtxAfterCuts->Fill(fZvtx);
    
    ////////////////////////////////////////////////////////////////////////////////

    if(fFindAntiLambda && fFindLambda)  //Events with Lambdas
    {
      hsphericityVSmultiL->Fill(TrackMult03,sphericity);
      hsphericityL->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTL->Fill(sphericity,totalpt);
      hsphericityVScentL->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTL->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinL->Fill(TrackMult03,totalpt);
    }
    if(fFindAntiLambda)  //Events with Anti Lambdas
    {
      hsphericityVSmultiAL->Fill(TrackMult03,sphericity);
      hsphericityAL->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTAL->Fill(sphericity,totalpt);
      hsphericityVScentAL->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTAL->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinAL->Fill(TrackMult03,totalpt);
    }
    if(fFindLambda)  //Events with Lambdas Only
    {
      hsphericityVSmultiLO->Fill(TrackMult03,sphericity);
      hsphericityLO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTLO->Fill(sphericity,totalpt);
      hsphericityVScentLO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTLO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinLO->Fill(TrackMult03,totalpt);
    }
    if(fFindLambda || fFindAntiLambda)  //Events with Lambdas or Anti Lambdas Only
    {
      hsphericityVSmultiALLO->Fill(TrackMult03,sphericity);
      hsphericityALLO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTALLO->Fill(sphericity,totalpt);
      hsphericityVScentALLO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTALLO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinALLO->Fill(TrackMult03,totalpt);
    }
    if(fFindXiPlus && fFindXiMinus)  //Events with Xis
    {
      hsphericityVSmultiX->Fill(TrackMult03,sphericity);
      hsphericityX->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTX->Fill(sphericity,totalpt);
      hsphericityVScentX->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTX->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinX->Fill(TrackMult03,totalpt);
    }
    if(fFindXiPlus)  //Events with Xis+
    {
      hsphericityVSmultiXP->Fill(TrackMult03,sphericity);
      hsphericityXP->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXP->Fill(sphericity,totalpt);
      hsphericityVScentXP->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXP->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXP->Fill(TrackMult03,totalpt);
    }
    if(fFindXiMinus)  //Events with Xis- Only
    {
      hsphericityVSmultiXN->Fill(TrackMult03,sphericity);
      hsphericityXN->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXN->Fill(sphericity,totalpt);
      hsphericityVScentXN->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXN->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXN->Fill(TrackMult03,totalpt);
    }
    if(fFindXiPlus || fFindXiMinus)  //Events with Xis+ or Xis- Only
    {
      hsphericityVSmultiXNXP->Fill(TrackMult03,sphericity);
      hsphericityXNXP->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXNXP->Fill(sphericity,totalpt);
      hsphericityVScentXNXP->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXNXP->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXNXP->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaPlus && fFindOmegaMinus)  //Events with Omegas
    {
      hsphericityVSmultiO->Fill(TrackMult03,sphericity);
      hsphericityO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTO->Fill(sphericity,totalpt);
      hsphericityVScentO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinO->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaPlus)  //Events with Omegas+
    {
      hsphericityVSmultiOP->Fill(TrackMult03,sphericity);
      hsphericityOP->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTOP->Fill(sphericity,totalpt);
      hsphericityVScentOP->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTOP->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinOP->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaMinus)  //Events with Omegas- Only
    {
      hsphericityVSmultiON->Fill(TrackMult03,sphericity);
      hsphericityON->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTON->Fill(sphericity,totalpt);
      hsphericityVScentON->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTON->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinON->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaPlus || fFindOmegaMinus)  //Events with Omegas+ or Omegas- Only
    {
      hsphericityVSmultiONOP->Fill(TrackMult03,sphericity);
      hsphericityONOP->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTONOP->Fill(sphericity,totalpt);
      hsphericityVScentONOP->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTONOP->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinONOP->Fill(TrackMult03,totalpt);
    }
    if(fFindKaon)  //Events with Kaons
    {
      hsphericityVSmultiKA->Fill(TrackMult03,sphericity);
      hsphericityKA->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTKA->Fill(sphericity,totalpt);
      hsphericityVScentKA->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTKA->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinKA->Fill(TrackMult03,totalpt);
    }
    ////////////////////////////////////////////////////////////////////////////////

    if((fFindAntiLambda && fFindLambda) && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Only Lambdas
    {
      hsphericityVSmultiL2->Fill(TrackMult03,sphericity);
      hsphericityL2->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTL2->Fill(sphericity,totalpt);
      hsphericityVScentL2->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTL2->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinL2->Fill(TrackMult03,totalpt);
    }
    if(fFindAntiLambda && !fFindLambda && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Anti Lambdas
    {
      hsphericityVSmultiALO->Fill(TrackMult03,sphericity);
      hsphericityALO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTALO->Fill(sphericity,totalpt);
      hsphericityVScentALO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTALO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinALO->Fill(TrackMult03,totalpt);
    }
    if(fFindLambda && !fFindAntiLambda && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Lambdas Only
    {
      hsphericityVSmultiLOO->Fill(TrackMult03,sphericity);
      hsphericityLOO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTLOO->Fill(sphericity,totalpt);
      hsphericityVScentLOO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTLOO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinLOO->Fill(TrackMult03,totalpt);
    }
    if((fFindLambda || fFindAntiLambda) && !fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon)  //Events with Lambdas or Anti Lambdas Only
    {
      hsphericityVSmultiALLOO->Fill(TrackMult03,sphericity);
      hsphericityALLOO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTALLOO->Fill(sphericity,totalpt);
      hsphericityVScentALLOO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTALLOO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinALLOO->Fill(TrackMult03,totalpt);
    }
    if((fFindXiPlus && fFindXiMinus) && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindLambda && !fFindAntiLambda)  //Events with Xis
    {
      hsphericityVSmultiXO->Fill(TrackMult03,sphericity);
      hsphericityXO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXO->Fill(sphericity,totalpt);
      hsphericityVScentXO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXO->Fill(TrackMult03,totalpt);
    }
    if(fFindXiPlus && !fFindXiMinus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindLambda && !fFindAntiLambda)  //Events with Xis+
    {
      hsphericityVSmultiXPO->Fill(TrackMult03,sphericity);
      hsphericityXPO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXPO->Fill(sphericity,totalpt);
      hsphericityVScentXPO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXPO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXPO->Fill(TrackMult03,totalpt);
    }
    if(fFindXiMinus && !fFindXiPlus && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindLambda && !fFindAntiLambda)  //Events with Xis- Only
    {
      hsphericityVSmultiXNO->Fill(TrackMult03,sphericity);
      hsphericityXNO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXNO->Fill(sphericity,totalpt);
      hsphericityVScentXNO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXNO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXNO->Fill(TrackMult03,totalpt);
    }
    if((fFindXiPlus || fFindXiMinus) && !fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindLambda && !fFindAntiLambda)  //Events with Xis+ or Xis- Only
    {
      hsphericityVSmultiXNXPO->Fill(TrackMult03,sphericity);
      hsphericityXNXPO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTXNXPO->Fill(sphericity,totalpt);
      hsphericityVScentXNXPO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTXNXPO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinXNXPO->Fill(TrackMult03,totalpt);
    }
    if((fFindOmegaPlus && fFindOmegaMinus) && !fFindKaon && !fFindLambda && !fFindAntiLambda && !fFindXiMinus && !fFindXiPlus)  //Events with Omegas
    {
      hsphericityVSmultiOO->Fill(TrackMult03,sphericity);
      hsphericityOO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTOO->Fill(sphericity,totalpt);
      hsphericityVScentOO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTOO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinOO->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaPlus && !fFindOmegaMinus && !fFindKaon && !fFindLambda && !fFindAntiLambda && !fFindXiMinus && !fFindXiPlus)  //Events with Omegas+
    {
      hsphericityVSmultiOPO->Fill(TrackMult03,sphericity);
      hsphericityOPO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTOPO->Fill(sphericity,totalpt);
      hsphericityVScentOPO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTOPO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinOPO->Fill(TrackMult03,totalpt);
    }
    if(fFindOmegaMinus && !fFindOmegaPlus && !fFindKaon && !fFindLambda && !fFindAntiLambda && !fFindXiMinus && !fFindXiPlus)  //Events with Omegas- Only
    {
      hsphericityVSmultiONO->Fill(TrackMult03,sphericity);
      hsphericityONO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTONO->Fill(sphericity,totalpt);
      hsphericityVScentONO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTONO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinONO->Fill(TrackMult03,totalpt);
    }
    if((fFindOmegaPlus || fFindOmegaMinus) && !fFindKaon && !fFindLambda && !fFindAntiLambda && !fFindXiMinus && !fFindXiPlus)  //Events with Omegas+ or Omegas- Only
    {
      hsphericityVSmultiONOPO->Fill(TrackMult03,sphericity);
      hsphericityONOPO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTONOPO->Fill(sphericity,totalpt);
      hsphericityVScentONOPO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTONOPO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinONOPO->Fill(TrackMult03,totalpt);
    }
    if(fFindKaon && !fFindLambda && !fFindAntiLambda && !fFindXiMinus && !fFindXiPlus && !fFindOmegaMinus && !fFindOmegaPlus)  //Events with Kaons
    {
      hsphericityVSmultiKAO->Fill(TrackMult03,sphericity);
      hsphericityKAO->Fill(TrackMult03,sphericity);         
      hsphericityVSMEANpTKAO->Fill(sphericity,totalpt);
      hsphericityVScentKAO->Fill(sphericity,CentPercentile);
      hmultiplicityVSMEANpTKAO->Fill(TrackMult03,totalpt); 
      hmultiplicityVSMEANpTbinKAO->Fill(TrackMult03,totalpt);
    }
  }   //end of analysis of triggered events
   
  // Post output data.
  PostData(1, fListOfObjects);
  
}     //end of main loop

//_____________________________________________________________________________
Bool_t AliAnalysisSphericityTask::IsAccepted(AliESDtrack* track) {
  
   Double_t gP = 0.0, gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
   AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
     gP = 0.0; gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
    }
    else {
      gP = tpcTrack->P();
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
    }

  // if(gPt < 0.25) return kFALSE; 
  // Checks if the track is excluded from the cuts
  if(fkQualityCutTPCrefit){
  	ULong_t Status = track->GetStatus();
     	if((Status&AliESDtrack::kTPCrefit) == 0) 
       	return kFALSE;
  }

  Int_t nClustersTPC = track->GetTPCclusters(0x0);

  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

    if(nClustersTPC < fMinDaughterTpcClusters) 
      return kFALSE;

    if(track->GetTPCsignalN() < fMinDaughterTpcClusters) 
      return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Double_t AliAnalysisSphericityTask::Rapidity(Double_t gPx, 
					      Double_t gPy, 
					      Double_t gPz, Int_t fType) const {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  if (fType == 0) fMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (fType == 1) fMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  
  Double_t gP = TMath::Sqrt(TMath::Power(gPx,2) + 
                           TMath::Power(gPy,2) + 
			   TMath::Power(gPz,2));
  Double_t energy = TMath::Sqrt(gP*gP + fMass*fMass);
  Double_t y = -999;
  if(energy != gPz) 
    y = 0.5*TMath::Log((energy + gPz)/(energy - gPz));

  return y;
}

//_____________________________________________________________________________
Bool_t AliAnalysisSphericityTask::IsProton(AliESDtrack *track) {
  //Function that checks if a track is a proton
  
   Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;
  
  //Ratio of the measured over the theoretical dE/dx a la STAR
  if(fPIDMode == kRatio) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      gEta = tpcTrack->Eta();
    }
    
    Double_t normalizeddEdx = -10.;
    if((track->GetTPCsignal() > 0.0) && (fTPCpid.GetExpectedSignal(gP,AliPID::kProton) > 0.0))
      normalizeddEdx = TMath::Log(track->GetTPCsignal()/fTPCpid.GetExpectedSignal(gP,AliPID::kProton));

    if (gP <= fNBoundP) if(normalizeddEdx >= fNRatio1) return kTRUE;
    if (gP >  fNBoundP) if(normalizeddEdx >= fNRatio2) return kTRUE;
  }//kRatio PID mode

  //Definition of an N-sigma area around the dE/dx vs P band
 else if(fPIDMode == kSigma) {
   
    Double_t nsigma = 100.0;
    
    Double_t mom = track->GetP();
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in)
      mom = in->GetP();

    nsigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
  
    if (mom <= fNBoundP) if(nsigma <= fNSigma1) return kTRUE;
    if (mom >  fNBoundP) if(nsigma <= fNSigma2) return kTRUE;
 }//kSigma PID method 

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisSphericityTask::IsPion(AliESDtrack *track) {
  //Function that checks if a track is a proton
   Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;

  //Ratio of the measured over the theoretical dE/dx a la STAR
  if(fPIDMode == kRatio) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      gEta = tpcTrack->Eta();
    }
   
    Double_t normalizeddEdx = -10.;
    if((track->GetTPCsignal() > 0.0) && (fTPCpid.GetExpectedSignal(gP,AliPID::kPion) > 0.0))
      normalizeddEdx = TMath::Log(track->GetTPCsignal()/fTPCpid.GetExpectedSignal(gP,AliPID::kPion));

    if (gP <= fNBoundP) if(normalizeddEdx >= fNRatio1) return kTRUE;
    if (gP >  fNBoundP) if(normalizeddEdx >= fNRatio2) return kTRUE;
      return kTRUE;
  }//kRatio PID mode

  //Definition of an N-sigma area around the dE/dx vs P band
  else if(fPIDMode == kSigma) {
   
    Double_t nsigma = 100.0;
    
    Double_t mom = track->GetP();
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in)
      mom = in->GetP();

    nsigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
  
    if (mom <= fNBoundP) if(nsigma <= fNSigma1) return kTRUE;
    if (mom >  fNBoundP) if(nsigma <= fNSigma2) return kTRUE;
 }//kSigma PID method 

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisSphericityTask::IsKaon(AliESDtrack *track) {
  //Function that checks if a track is a proton
   Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;

  //Ratio of the measured over the theoretical dE/dx a la STAR
  if(fPIDMode == kRatio) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      gEta = tpcTrack->Eta();
    }
   
    Double_t normalizeddEdx = -10.;
    if((track->GetTPCsignal() > 0.0) && (fTPCpid.GetExpectedSignal(gP,AliPID::kKaon) > 0.0))
      normalizeddEdx = TMath::Log(track->GetTPCsignal()/fTPCpid.GetExpectedSignal(gP,AliPID::kKaon));

    if (gP <= fNBoundP) if(normalizeddEdx >= fNRatio1) return kTRUE;
    if (gP >  fNBoundP) if(normalizeddEdx >= fNRatio2) return kTRUE;
      return kTRUE;
  }//kRatio PID mode

  //Definition of an N-sigma area around the dE/dx vs P band
  else if(fPIDMode == kSigma) {
   
    Double_t nsigma = 100.0;
    
    Double_t mom = track->GetP();
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in)
      mom = in->GetP();

    nsigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
  
    if (mom <= fNBoundP) if(nsigma <= fNSigma1) return kTRUE;
    if (mom >  fNBoundP) if(nsigma <= fNSigma2) return kTRUE;
 }//kSigma PID method 

  return kFALSE;
}

//_____________________________________________________________________________
Float_t AliAnalysisSphericityTask::GetVertex(const AliVEvent* event) const
{
  Float_t zvtx = -999;
  
  const AliVVertex* primaryVertex = event->GetPrimaryVertex();   
  if(primaryVertex->GetNContributors()>0)
    zvtx = primaryVertex->GetZ();
  
  return zvtx;
}


//_____________________________________________________________________________
Float_t AliAnalysisSphericityTask::GetTest(){
  
  return 10.0;
  
}


//_____________________________________________________________________________
void AliAnalysisSphericityTask::Terminate(Option_t *)
{  
  TFile* fout = new TFile("Sphericity_qa.root", "RECREATE"); 
  if (fESASelection)
    {
      fESASelection->SaveHistosSt();
    }
  
  fout->Write();
  fout->Close();
  
  Printf("Writing result to Sphericity_qa.root");
}