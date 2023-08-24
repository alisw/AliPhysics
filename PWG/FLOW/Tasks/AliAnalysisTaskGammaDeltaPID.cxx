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
/* $Id: AliAnalysisTaskGammaDeltaPID.cxx ver: 2.0                     $   */
/* Simple Task to fill V0 and ZDC Energies for Gain Calibration           */
/* Works with 15o and 18q/r. Support to be added for LHC10h               */
/* Developer: Md Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)        */
/* Last Modified: Aug 23, 2021,  First version committed                  */
/* Last Modified: Oct 08, 2021,  Second version committed                 */
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskGammaDeltaPID.h"
#include "AliInputEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliVEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliTimeRangeCut.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliCentrality.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliVParticle.h"
#include "AliAODVZERO.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliAODZDC.h"
#include "AliStack.h"
#include "AliAODv0.h"

#include "TMatrixDSym.h"
#include "TParticle.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TChain.h"
#include "TMath.h"
#include "stdio.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include <iostream>
#include <vector>


using namespace std;

ClassImp(AliAnalysisTaskGammaDeltaPID)

AliAnalysisTaskGammaDeltaPID::AliAnalysisTaskGammaDeltaPID(const char *name):
  AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fListZDCCorr(NULL),
  fListTemp(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  
  gHarmonic(0),
  gParticleID(0),
  fFilterBit(1),
  fTPCclustMin(70),
  gOldRunNumber(1),  
  fMinVzCut(-10.0),
  fMaxVzCut(+10.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzzMax(3.2),
  fMinEtaCut(-0.8),
  fMaxEtaCut(+0.8),
  fEtaGapNeg(-0.1),
  fEtaGapPos(+0.1),
  fTrkChi2Min(0.1),
  fTrkChi2Max(4.0),    
  fTPCdEdxMin(10.),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fVertexZEvent(-99.),
  sDetectorForEP("bla1"), 
  sCentrEstimator("V0M"),  
  bUseKinkTracks(kFALSE),
  bSkipPileUpCut(kFALSE),
  bSkipNestedLoop(kFALSE),
  bUseV0EventPlane(kFALSE),
  bAnalysLambdaPairs(kFALSE),
  bUseZDCSpectatorPlane(kFALSE),
  gTypeOfRecentering(-1),
  bRecenterFailOrNot(kTRUE),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersNsigma(3.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  fV0DaughterUseTOF(kFALSE),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fMassMean(1.115683),
  fLambdaMassCut(0.005),
  vecPosEPTrkID(0),
  vecNegEPTrkID(0),
  vecPosEPTrkPhi(0),
  vecNegEPTrkPhi(0),
  vecPosEPTrkNUAWgt(0),
  vecNegEPTrkNUAWgt(0),
  
  fHistVertexZcm(NULL),
  fHistAnalysisInfo(NULL),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fDebugwEventCount(NULL),
  fHCorrectQNxV0C(NULL),
  fHCorrectQNyV0C(NULL),    
  fHCorrectQNxV0A(NULL),
  fHCorrectQNyV0A(NULL),  
  fHCorrectQ3xV0C(NULL),
  fHCorrectQ3yV0C(NULL),    
  fHCorrectQ3xV0A(NULL),
  fHCorrectQ3yV0A(NULL),  
  fHCorrectMCPosChrg(NULL),
  fHCorrectMCPosPion(NULL),
  fHCorrectMCPosKaon(NULL),
  fHCorrectMCPosProt(NULL),
  fHCorrectMCNegChrg(NULL),
  fHCorrectMCNegPion(NULL),
  fHCorrectMCNegKaon(NULL),
  fHCorrectMCNegProt(NULL),
  fHCorrectTPCQnxEtaPos(NULL),
  fHCorrectTPCQnyEtaPos(NULL),
  fHCorrectTPCQnxEtaNeg(NULL),
  fHCorrectTPCQnyEtaNeg(NULL),
  fHCorrectTPCQ3xEtaPos(NULL),
  fHCorrectTPCQ3yEtaPos(NULL),
  fHCorrectTPCQ3xEtaNeg(NULL),
  fHCorrectTPCQ3yEtaNeg(NULL),   
  
  fHistTPCPsiNPosPlane(NULL),
  fHistTPCPsiNNegPlane(NULL),
  fHistTPCPsi3PosPlane(NULL),
  fHistTPCPsi3NegPlane(NULL),
  fHistTPCPsi4PosPlane(NULL),
  fHistTPCPsi4NegPlane(NULL),
  
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  fHistV0CPsi3EventPlane(NULL),
  fHistV0APsi3EventPlane(NULL),
  fHistTPCPosqVectorvsCent(NULL),
  fHistTPCNegqVectorvsCent(NULL),
  fHistV0CDetqVectorvsCent(NULL),
  fHistV0ADetqVectorvsCent(NULL),

  fHCorrectV0ChWeghts(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsCL1After(NULL), 
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsV0MAfter(NULL), 
  fHistCentCL0VsV0MBefore(NULL),
  fHistCentCL0VsV0MAfter(NULL), 
  fHistTPCVsESDTrkBefore(NULL), 
  fHistTPCVsESDTrkAfter(NULL),
  fHCorrectNUAChrgPos(NULL),
  fHCorrectNUAChrgNeg(NULL),
  fHCorrectNUAkPIDPos(NULL),
  fHCorrectNUAkPIDNeg(NULL),
  
  fHZDCCparameters(NULL),
  fHZDCAparameters(NULL),
  
  hZDCCPsiSpectatorPlane(NULL),
  hZDCAPsiSpectatorPlane(NULL),
  hZDCCAPsiSpectatorPlane(NULL),
	
  hZDCCPsivsTPCPsi2(NULL),
  hZDCCPsivsTPCPosPsi2(NULL),
  hZDCCPsivsTPCNegPsi2(NULL),
  hZDCCPsivsV0CAPsi2(NULL),
  hZDCCPsivsV0CPsi2(NULL),
  hZDCCPsivsV0APsi2(NULL),
	
  hZDCAPsivsTPCPsi2(NULL),
  hZDCAPsivsTPCPosPsi2(NULL),
  hZDCAPsivsTPCNegPsi2(NULL),
  hZDCAPsivsV0CAPsi2(NULL),
  hZDCAPsivsV0CPsi2(NULL),
  hZDCAPsivsV0APsi2(NULL),
	
  hZDCCAPsivsTPCPsi2(NULL),
  hZDCCAPsivsTPCPosPsi2(NULL),
  hZDCCAPsivsTPCNegPsi2(NULL),
  hZDCCAPsivsV0CAPsi2(NULL),
  hZDCCAPsivsV0CPsi2(NULL),
  hZDCCAPsivsV0APsi2(NULL),
  
  hV2TPCvsCent(NULL),
  hV2V0CvsCent(NULL),
  hV2V0AvsCent(NULL),
  hV2V0CAvsCent(NULL),
  hV2ZDCCvsCent(NULL),
  hV2ZDCAvsCent(NULL),
  hV2ZDCCACombinevsCent(NULL),
  hV2ZDCCAvsCent(NULL),
  
  hPOIQRe(NULL),
  hPOIQIm(NULL),
  hPOIMul(NULL),
  hFlowQC2Pro(NULL),
  hFlowQC4Pro(NULL),
  hFlowwCov24Pro(NULL),
  
  hAvg3pC112vsCentPP(NULL),
  hAvg3pC112vsCentNN(NULL),    
  hAvg3pC112vsCentOS(NULL),
  hAvg3pC123vsCentPP(NULL),
  hAvg3pC123vsCentNN(NULL),    
  hAvg3pC123vsCentOS(NULL),
  hAvgDelta1vsCentPP(NULL),
  hAvgDelta1vsCentNN(NULL),    
  hAvgDelta1vsCentOS(NULL),
  hAvgDelta2vsCentPP(NULL),
  hAvgDelta2vsCentNN(NULL),    
  hAvgDelta2vsCentOS(NULL),
  hAvgDelta3vsCentPP(NULL),
  hAvgDelta3vsCentNN(NULL),    
  hAvgDelta3vsCentOS(NULL),
  hAvgDelta4vsCentPP(NULL),
  hAvgDelta4vsCentNN(NULL),    
  hAvgDelta4vsCentOS(NULL),
  
  hAvg3pC112vsCentPPUsingZDCC(NULL),
  hAvg3pC112vsCentNNUsingZDCC(NULL),
  hAvg3pC112vsCentOSUsingZDCC(NULL),
  hAvg3pC112vsCentPPUsingZDCA(NULL),
  hAvg3pC112vsCentNNUsingZDCA(NULL),
  hAvg3pC112vsCentOSUsingZDCA(NULL),
  hAvg3pC112vsCentPPUsingZDCCACombine(NULL),
  hAvg3pC112vsCentNNUsingZDCCACombine(NULL),
  hAvg3pC112vsCentOSUsingZDCCACombine(NULL),
  hAvg3pC112vsCentPPUsingZDCCA(NULL),
  hAvg3pC112vsCentNNUsingZDCCA(NULL),
  hAvg3pC112vsCentOSUsingZDCCA(NULL),
  
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQ3XvsCentV0C(NULL),
  hAvgQ3YvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hAvgQ3XvsCentV0A(NULL),
  hAvgQ3YvsCentV0A(NULL),
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hTPCPsi4Correlation(NULL),  
  hV0CV0APsiNCorrelation(NULL),
  hV0CTPCPsiNCorrelation(NULL),
  hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  hV0CTPCPsi3Correlation(NULL),
  hV0ATPCPsi3Correlation(NULL),  
  
  fAvgCos2PsivsCentEtaPos(NULL),
  fAvgSin2PsivsCentEtaPos(NULL),
  fAvgCos2PsivsCentEtaNeg(NULL),
  fAvgSin2PsivsCentEtaNeg(NULL),
  fAvgCos3PsivsCentEtaPos(NULL),
  fAvgSin3PsivsCentEtaPos(NULL),
  fAvgCos3PsivsCentEtaNeg(NULL),
  fAvgSin3PsivsCentEtaNeg(NULL),
  fAvgCos4PsivsCentEtaPos(NULL),
  fAvgSin4PsivsCentEtaPos(NULL),
  fAvgCos4PsivsCentEtaNeg(NULL),
  fAvgSin4PsivsCentEtaNeg(NULL),
  hAvgV0ChannelsvsVz(NULL),

  fHistV0Pt(NULL),              
  fHistV0Eta(NULL),              
  fHistV0DcatoPrimVertex(NULL), 
  fHistV0CPA(NULL),             
  fHistV0DecayLength(NULL),     
  fProfileDelta_Lambda_hPos(NULL),
  fProfileDelta_Lambda_hNeg(NULL),
  fProfileDelta_Lambda_Proton(NULL),   
  fProfileDelta_Lambda_AntiProton(NULL),
  fProfileDelta_AntiLambda_hPos(NULL),
  fProfileDelta_AntiLambda_hNeg(NULL),
  fProfileDelta_AntiLambda_Proton(NULL),
  fProfileDelta_AntiLambda_AntiProton(NULL),    
  fProfileGammaTPC_Lambda_hPos(NULL),
  fProfileGammaTPC_Lambda_hNeg(NULL),
  fProfileGammaTPC_Lambda_Proton(NULL),
  fProfileGammaTPC_Lambda_AntiProton(NULL),
  fProfileGammaTPC_AntiLambda_hPos(NULL),
  fProfileGammaTPC_AntiLambda_hNeg(NULL),
  fProfileGammaTPC_AntiLambda_Proton(NULL),
  fProfileGammaTPC_AntiLambda_AntiProton(NULL),  
  fProfileGammaV0C_Lambda_hPos(NULL),
  fProfileGammaV0C_Lambda_hNeg(NULL),
  fProfileGammaV0C_Lambda_Proton(NULL),
  fProfileGammaV0C_Lambda_AntiProton(NULL),
  fProfileGammaV0C_AntiLambda_hPos(NULL),
  fProfileGammaV0C_AntiLambda_hNeg(NULL),
  fProfileGammaV0C_AntiLambda_Proton(NULL),
  fProfileGammaV0C_AntiLambda_AntiProton(NULL),
  fProfileGammaV0A_Lambda_hPos(NULL),
  fProfileGammaV0A_Lambda_hNeg(NULL),
  fProfileGammaV0A_Lambda_Proton(NULL),
  fProfileGammaV0A_Lambda_AntiProton(NULL),
  fProfileGammaV0A_AntiLambda_hPos(NULL),
  fProfileGammaV0A_AntiLambda_hNeg(NULL),
  fProfileGammaV0A_AntiLambda_Proton(NULL),
  fProfileGammaV0A_AntiLambda_AntiProton(NULL),
  fProfileGammaZNC_Lambda_hPos(NULL),
  fProfileGammaZNC_Lambda_hNeg(NULL),
  fProfileGammaZNC_Lambda_Proton(NULL),
  fProfileGammaZNC_Lambda_AntiProton(NULL),
  fProfileGammaZNC_AntiLambda_hPos(NULL),
  fProfileGammaZNC_AntiLambda_hNeg(NULL),
  fProfileGammaZNC_AntiLambda_Proton(NULL),
  fProfileGammaZNC_AntiLambda_AntiProton(NULL),
  fProfileGammaZNA_Lambda_hPos(NULL),
  fProfileGammaZNA_Lambda_hNeg(NULL),
  fProfileGammaZNA_Lambda_Proton(NULL),
  fProfileGammaZNA_Lambda_AntiProton(NULL),
  fProfileGammaZNA_AntiLambda_hPos(NULL),
  fProfileGammaZNA_AntiLambda_hNeg(NULL),
  fProfileGammaZNA_AntiLambda_Proton(NULL),
  fProfileGammaZNA_AntiLambda_AntiProton(NULL),
  hEmptyPointerFortheList(NULL),
  bCheckPIDFlow(kFALSE)  
{
  
  //std::vector<Int_t> vecPosEPTrkID = {0};
  //std::vector<Int_t> vecNegEPTrkID = {0};
  
  for(int i=0;i<3;i++){
    fCurrentVtx[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    fHistLambdaPt[i]              = NULL;            
    fHistLambdaEta[i]             = NULL;
    fHistLambdaPhi[i]             = NULL;
    fHistLambdaDcaToPrimVertex[i] = NULL;
    fHistLambdaCPA[i]             = NULL;
    fHistLambdaDecayLength[i]     = NULL;
    fHistLambdaMass[i]            = NULL;
    fHistAntiLambdaPt[i]          = NULL;
    fHistAntiLambdaEta[i]         = NULL;
    fHistAntiLambdaPhi[i]         = NULL;
    fHistAntiLambdaDcaToPrimVertex[i] = NULL;
    fHistAntiLambdaCPA[i]         = NULL;
    fHistAntiLambdaDecayLength[i] = NULL;
    fHistAntiLambdaMass[i]        = NULL;
    fProfileLambdaMassVsPt[i]     = NULL;
    fProfileAntiLambdaMassVsPt[i] = NULL;
  }  

  for (int i = 0; i < 5; i++) {
    fProfile2DRawFlowCentPthPos[i] = NULL;
    fProfile2DRawFlowCentPthNeg[i] = NULL;  
    fProfile2DRawFlowCentPtProton[i] = NULL;
    fProfile2DRawFlowCentPtAntiProton[i] = NULL;
  }
  for (int i = 0; i < 6; i++) {
    fProfile2DRawFlowCentPtLambda[i] = NULL;
    fProfile2DRawFlowCentPtAntiLambda[i] = NULL;
  }  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaDeltaPID::AliAnalysisTaskGammaDeltaPID():
  AliAnalysisTaskSE(),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fListZDCCorr(NULL),
  fListTemp(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),

  gHarmonic(0),
  gParticleID(0),
  fFilterBit(1),
  fTPCclustMin(70),
  gOldRunNumber(1),  
  fMinVzCut(-10.0),
  fMaxVzCut(+10.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzzMax(3.2),
  fMinEtaCut(-0.8),
  fMaxEtaCut(+0.8),
  fEtaGapNeg(-0.1),
  fEtaGapPos(+0.1),
  fTrkChi2Min(0.1),
  fTrkChi2Max(4.0),    
  fTPCdEdxMin(10.),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fVertexZEvent(-99.),
  sDetectorForEP("bla2"), 
  sCentrEstimator("V0M"),  
  bUseKinkTracks(kFALSE),
  bSkipPileUpCut(kFALSE),
  bSkipNestedLoop(kFALSE),
  bUseV0EventPlane(kFALSE),
  bAnalysLambdaPairs(kFALSE),
  bUseZDCSpectatorPlane(kFALSE),
  gTypeOfRecentering(-1),
  bRecenterFailOrNot(kTRUE),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersNsigma(3.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  fV0DaughterUseTOF(kFALSE),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fMassMean(1.115683),
  fLambdaMassCut(0.005),
  vecPosEPTrkID(0),
  vecNegEPTrkID(0),
  vecPosEPTrkPhi(0),
  vecNegEPTrkPhi(0),
  vecPosEPTrkNUAWgt(0),
  vecNegEPTrkNUAWgt(0),
  
  fHistVertexZcm(NULL),
  fHistAnalysisInfo(NULL),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fDebugwEventCount(NULL),
  fHCorrectQNxV0C(NULL),
  fHCorrectQNyV0C(NULL),    
  fHCorrectQNxV0A(NULL),
  fHCorrectQNyV0A(NULL),  
  fHCorrectQ3xV0C(NULL),
  fHCorrectQ3yV0C(NULL),    
  fHCorrectQ3xV0A(NULL),
  fHCorrectQ3yV0A(NULL),   
  fHCorrectMCPosChrg(NULL),
  fHCorrectMCPosPion(NULL),
  fHCorrectMCPosKaon(NULL),
  fHCorrectMCPosProt(NULL),
  fHCorrectMCNegChrg(NULL),
  fHCorrectMCNegPion(NULL),
  fHCorrectMCNegKaon(NULL),
  fHCorrectMCNegProt(NULL), 
  fHCorrectTPCQnxEtaPos(NULL),
  fHCorrectTPCQnyEtaPos(NULL),
  fHCorrectTPCQnxEtaNeg(NULL),
  fHCorrectTPCQnyEtaNeg(NULL),
  fHCorrectTPCQ3xEtaPos(NULL),
  fHCorrectTPCQ3yEtaPos(NULL),
  fHCorrectTPCQ3xEtaNeg(NULL),
  fHCorrectTPCQ3yEtaNeg(NULL),   

  
  fHistTPCPsiNPosPlane(NULL),
  fHistTPCPsiNNegPlane(NULL),
  fHistTPCPsi3PosPlane(NULL),
  fHistTPCPsi3NegPlane(NULL),
  fHistTPCPsi4PosPlane(NULL),
  fHistTPCPsi4NegPlane(NULL),
  
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  fHistV0CPsi3EventPlane(NULL),
  fHistV0APsi3EventPlane(NULL),
  fHistTPCPosqVectorvsCent(NULL),
  fHistTPCNegqVectorvsCent(NULL),
  fHistV0CDetqVectorvsCent(NULL),
  fHistV0ADetqVectorvsCent(NULL),

  fHCorrectV0ChWeghts(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsCL1After(NULL), 
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsV0MAfter(NULL), 
  fHistCentCL0VsV0MBefore(NULL),
  fHistCentCL0VsV0MAfter(NULL), 
  fHistTPCVsESDTrkBefore(NULL), 
  fHistTPCVsESDTrkAfter(NULL),
  fHCorrectNUAChrgPos(NULL),
  fHCorrectNUAChrgNeg(NULL),
  fHCorrectNUAkPIDPos(NULL),
  fHCorrectNUAkPIDNeg(NULL),
  
  fHZDCCparameters(NULL),
  fHZDCAparameters(NULL),
  
  hZDCCPsiSpectatorPlane(NULL),
  hZDCAPsiSpectatorPlane(NULL),
  hZDCCAPsiSpectatorPlane(NULL),
	
  hZDCCPsivsTPCPsi2(NULL),
  hZDCCPsivsTPCPosPsi2(NULL),
  hZDCCPsivsTPCNegPsi2(NULL),
  hZDCCPsivsV0CAPsi2(NULL),
  hZDCCPsivsV0CPsi2(NULL),
  hZDCCPsivsV0APsi2(NULL),
	
  hZDCAPsivsTPCPsi2(NULL),
  hZDCAPsivsTPCPosPsi2(NULL),
  hZDCAPsivsTPCNegPsi2(NULL),
  hZDCAPsivsV0CAPsi2(NULL),
  hZDCAPsivsV0CPsi2(NULL),
  hZDCAPsivsV0APsi2(NULL),
	
  hZDCCAPsivsTPCPsi2(NULL),
  hZDCCAPsivsTPCPosPsi2(NULL),
  hZDCCAPsivsTPCNegPsi2(NULL),
  hZDCCAPsivsV0CAPsi2(NULL),
  hZDCCAPsivsV0CPsi2(NULL),
  hZDCCAPsivsV0APsi2(NULL),
  
  hV2TPCvsCent(NULL),
  hV2V0CvsCent(NULL),
  hV2V0AvsCent(NULL),
  hV2V0CAvsCent(NULL),
  hV2ZDCCvsCent(NULL),
  hV2ZDCAvsCent(NULL),
  hV2ZDCCACombinevsCent(NULL),
  hV2ZDCCAvsCent(NULL),
  
  hPOIQRe(NULL),
  hPOIQIm(NULL),
  hPOIMul(NULL),
  hFlowQC2Pro(NULL),
  hFlowQC4Pro(NULL),
  hFlowwCov24Pro(NULL),
  
  hAvg3pC112vsCentPP(NULL),
  hAvg3pC112vsCentNN(NULL),    
  hAvg3pC112vsCentOS(NULL),
  hAvg3pC123vsCentPP(NULL),
  hAvg3pC123vsCentNN(NULL),    
  hAvg3pC123vsCentOS(NULL),
  hAvgDelta1vsCentPP(NULL),
  hAvgDelta1vsCentNN(NULL),    
  hAvgDelta1vsCentOS(NULL),
  hAvgDelta2vsCentPP(NULL),
  hAvgDelta2vsCentNN(NULL),    
  hAvgDelta2vsCentOS(NULL),
  hAvgDelta3vsCentPP(NULL),
  hAvgDelta3vsCentNN(NULL),    
  hAvgDelta3vsCentOS(NULL),
  hAvgDelta4vsCentPP(NULL),
  hAvgDelta4vsCentNN(NULL),    
  hAvgDelta4vsCentOS(NULL),
  
  hAvg3pC112vsCentPPUsingZDCC(NULL),
  hAvg3pC112vsCentNNUsingZDCC(NULL),
  hAvg3pC112vsCentOSUsingZDCC(NULL),
  hAvg3pC112vsCentPPUsingZDCA(NULL),
  hAvg3pC112vsCentNNUsingZDCA(NULL),
  hAvg3pC112vsCentOSUsingZDCA(NULL),
  hAvg3pC112vsCentPPUsingZDCCACombine(NULL),
  hAvg3pC112vsCentNNUsingZDCCACombine(NULL),
  hAvg3pC112vsCentOSUsingZDCCACombine(NULL),
  hAvg3pC112vsCentPPUsingZDCCA(NULL),
  hAvg3pC112vsCentNNUsingZDCCA(NULL),
  hAvg3pC112vsCentOSUsingZDCCA(NULL),
  
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQ3XvsCentV0C(NULL),
  hAvgQ3YvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hAvgQ3XvsCentV0A(NULL),
  hAvgQ3YvsCentV0A(NULL),
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hTPCPsi4Correlation(NULL),    
  hV0CV0APsiNCorrelation(NULL),
  hV0CTPCPsiNCorrelation(NULL),
  hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  hV0CTPCPsi3Correlation(NULL),
  hV0ATPCPsi3Correlation(NULL),  
  
  fAvgCos2PsivsCentEtaPos(NULL),
  fAvgSin2PsivsCentEtaPos(NULL),
  fAvgCos2PsivsCentEtaNeg(NULL),
  fAvgSin2PsivsCentEtaNeg(NULL),
  fAvgCos3PsivsCentEtaPos(NULL),
  fAvgSin3PsivsCentEtaPos(NULL),
  fAvgCos3PsivsCentEtaNeg(NULL),
  fAvgSin3PsivsCentEtaNeg(NULL),
  fAvgCos4PsivsCentEtaPos(NULL),
  fAvgSin4PsivsCentEtaPos(NULL),
  fAvgCos4PsivsCentEtaNeg(NULL),
  fAvgSin4PsivsCentEtaNeg(NULL),  
  hAvgV0ChannelsvsVz(NULL),
  
  fHistV0Pt(NULL),              
  fHistV0Eta(NULL),              
  fHistV0DcatoPrimVertex(NULL), 
  fHistV0CPA(NULL),             
  fHistV0DecayLength(NULL),     
  fProfileDelta_Lambda_hPos(NULL),
  fProfileDelta_Lambda_hNeg(NULL),
  fProfileDelta_Lambda_Proton(NULL),   
  fProfileDelta_Lambda_AntiProton(NULL),
  fProfileDelta_AntiLambda_hPos(NULL),
  fProfileDelta_AntiLambda_hNeg(NULL),
  fProfileDelta_AntiLambda_Proton(NULL),
  fProfileDelta_AntiLambda_AntiProton(NULL),
  fProfileGammaTPC_Lambda_hPos(NULL),
  fProfileGammaTPC_Lambda_hNeg(NULL),
  fProfileGammaTPC_Lambda_Proton(NULL),
  fProfileGammaTPC_Lambda_AntiProton(NULL),
  fProfileGammaTPC_AntiLambda_hPos(NULL),
  fProfileGammaTPC_AntiLambda_hNeg(NULL),
  fProfileGammaTPC_AntiLambda_Proton(NULL),
  fProfileGammaTPC_AntiLambda_AntiProton(NULL),    
  fProfileGammaV0C_Lambda_hPos(NULL),
  fProfileGammaV0C_Lambda_hNeg(NULL),
  fProfileGammaV0C_Lambda_Proton(NULL),
  fProfileGammaV0C_Lambda_AntiProton(NULL),
  fProfileGammaV0C_AntiLambda_hPos(NULL),
  fProfileGammaV0C_AntiLambda_hNeg(NULL),
  fProfileGammaV0C_AntiLambda_Proton(NULL),
  fProfileGammaV0C_AntiLambda_AntiProton(NULL),
  fProfileGammaV0A_Lambda_hPos(NULL),
  fProfileGammaV0A_Lambda_hNeg(NULL),
  fProfileGammaV0A_Lambda_Proton(NULL),
  fProfileGammaV0A_Lambda_AntiProton(NULL),
  fProfileGammaV0A_AntiLambda_hPos(NULL),
  fProfileGammaV0A_AntiLambda_hNeg(NULL),
  fProfileGammaV0A_AntiLambda_Proton(NULL),
  fProfileGammaV0A_AntiLambda_AntiProton(NULL),
  fProfileGammaZNC_Lambda_hPos(NULL),
  fProfileGammaZNC_Lambda_hNeg(NULL),
  fProfileGammaZNC_Lambda_Proton(NULL),
  fProfileGammaZNC_Lambda_AntiProton(NULL),
  fProfileGammaZNC_AntiLambda_hPos(NULL),
  fProfileGammaZNC_AntiLambda_hNeg(NULL),
  fProfileGammaZNC_AntiLambda_Proton(NULL),
  fProfileGammaZNC_AntiLambda_AntiProton(NULL),
  fProfileGammaZNA_Lambda_hPos(NULL),
  fProfileGammaZNA_Lambda_hNeg(NULL),
  fProfileGammaZNA_Lambda_Proton(NULL),
  fProfileGammaZNA_Lambda_AntiProton(NULL),
  fProfileGammaZNA_AntiLambda_hPos(NULL),
  fProfileGammaZNA_AntiLambda_hNeg(NULL),
  fProfileGammaZNA_AntiLambda_Proton(NULL),
  fProfileGammaZNA_AntiLambda_AntiProton(NULL),
  hEmptyPointerFortheList(NULL),
  bCheckPIDFlow(kFALSE)
{
  //std::vector<Int_t> vecPosEPTrkID = {0};
  //std::vector<Int_t> vecNegEPTrkID = {0};
  
  for(int i=0;i<3;i++){
    fCurrentVtx[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    fHistLambdaPt[i]              = NULL;            
    fHistLambdaEta[i]             = NULL;
    fHistLambdaPhi[i]             = NULL;
    fHistLambdaDcaToPrimVertex[i] = NULL;
    fHistLambdaCPA[i]             = NULL;
    fHistLambdaDecayLength[i]     = NULL;
    fHistLambdaMass[i]            = NULL;
    fHistAntiLambdaPt[i]          = NULL;
    fHistAntiLambdaEta[i]         = NULL;
    fHistAntiLambdaPhi[i]         = NULL;
    fHistAntiLambdaDcaToPrimVertex[i] = NULL;
    fHistAntiLambdaCPA[i]         = NULL;
    fHistAntiLambdaDecayLength[i] = NULL;
    fHistAntiLambdaMass[i]        = NULL;
    fProfileLambdaMassVsPt[i]     = NULL;
    fProfileAntiLambdaMassVsPt[i] = NULL;
  }  
  for (int i = 0; i < 5; i++) {
    fProfile2DRawFlowCentPthPos[i] = NULL;
    fProfile2DRawFlowCentPthNeg[i] = NULL;  
    fProfile2DRawFlowCentPtProton[i] = NULL;
    fProfile2DRawFlowCentPtAntiProton[i] = NULL;

  }
  for (int i = 0; i < 6; i++) {
    fProfile2DRawFlowCentPtLambda[i] = NULL;
    fProfile2DRawFlowCentPtAntiLambda[i] = NULL;
  }
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskGammaDeltaPID::~AliAnalysisTaskGammaDeltaPID()
{
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!
  if(fListHist)      delete fListHist;  
  if(fListTRKCorr)   delete fListTRKCorr;
  if(fListNUACorr)   delete fListNUACorr;
  if(fListV0MCorr)   delete fListV0MCorr;
  if(fListZDCCorr)   delete fListZDCCorr;
  if(fListTemp)      delete fListTemp;

  if(fV0CutPU)      delete fV0CutPU;
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU;   
}










//________________ Define Histograms _______________
void AliAnalysisTaskGammaDeltaPID::UserCreateOutputObjects()
{
  //Get The Input Hander:
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  if (!inputHandler) {  printf("\n***** ERROR *****\n Input handler missing, Status:QUIT!\n");    exit(1);}
  
  //// obtain the PID response object if needed:
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) {  printf("\n***** ERROR *****\n fPIDResponse missing, Status:QUIT!\n");    exit(1);}    
    
  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);
  

  //----------- User's Analysis histograms Defined here --------------
  Char_t  name[100];
  Char_t title[100]; 
  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables
  

  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  fListTemp = new TList();
  fListTemp->SetOwner(kTRUE);
  
  SetupQAHistograms();
  SetupAnalysisHistograms();
  SetupPileUpRemovalFunctions();
  SetupEventAndTaskConfigInfo();
  if (bCheckPIDFlow) SetupLambdaPIDFlowHistograms();



 
    
  cout<<"Ncls: "<<fTPCclustMin<<" Harm: "<<gHarmonic<<" POI: "<<gParticleID<<" nsigTPC: "<<fNSigmaTPCCut<<" nsigCirc: "<<fNSigmaTOFCut<<endl;
  cout<<"FB: "<<fFilterBit<<" chi2min: "<<fTrkChi2Min<<" chi2max: "<<fTrkChi2Max<<" etaMin: "<<fMinEtaCut<<" etaMax: "<<fMaxEtaCut<<endl;
  cout<<"dEdxMin: "<<fTPCdEdxMin<<" dcaXY: "<<fDCAxyMax<<" dcaZ: "<<fDCAzzMax<<"  VzLow: "<<fMinVzCut<<" VzHigh: "<<fMaxVzCut<<endl;
  cout<<"minPt: "<<fMinPtCut<<" maxPt: "<<fMaxPtCut<<" etaGapNeg: "<<fEtaGapNeg<<" etaGapPos: "<<fEtaGapPos<<" oldRun: "<<gOldRunNumber<<endl;
  cout<<"Centrality Estimator: "<<sCentrEstimator.Data()<<" EventPlane Detector:"<<sDetectorForEP.Data()<<endl;
  
  //Set Up MC Correction file
  if(fListTRKCorr && !fHCorrectMCPosChrg){
    GetMCCorrectionHist();  
  }
  //cout<<"\n ******  Bug Testing Mode.. So we exit here ****** \n"<<endl; exit(111);
  
     
  PostData(1,fListHist);
}










//____________________________ Call Event by Event ___________________________________
void AliAnalysisTaskGammaDeltaPID::UserExec(Option_t*) {
 
  //std::cout<<" Info:UserExec() called ..!!!\n";
  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  fESD = dynamic_cast <AliESDEvent*> (InputEvent());
  if(!(fESD || fAOD)) {
    printf("ERROR: fESD & fAOD not available\n");
    return;
  }

  fDebugwEventCount->Fill(0.1);


  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent){
    printf("ERROR: fVevent not available\n");
    return;
  }



  
  //-------------- Vtx cuts ---------------
  const AliVVertex *pointVtx = fVevent->GetPrimaryVertex();
  Double_t pVtxZ = -999, pVtxX = -999, pVtxY=-999;
  pVtxZ  = pointVtx->GetZ();
  pVtxX  = pointVtx->GetX();
  pVtxY  = pointVtx->GetY();
  
  if(pVtxZ < fMinVzCut || pVtxZ > fMaxVzCut){
    return;
  }

  fVertexZEvent = pVtxZ;   /// Global variable as needed by many functions.
  
  fDebugwEventCount->Fill(1.1);

  
  /// Print the Cuts or Variables Ranges As set in AddTask: if this Doesn't print what we want then we are in Trouble!!!
  //cout<<"Ncls: "<<fTPCclustMin<<" Harm: "<<gHarmonic<<" POI: "<<gParticleID<<" nsigTPC: "<<fNSigmaTPCCut<<" nsigCirc: "<<fNSigmaTOFCut<<endl;
  //cout<<"FB: "<<fFilterBit<<" chi2min: "<<fTrkChi2Min<<" chi2max: "<<fTrkChi2Max<<" etaMin: "<<fMinEtaCut<<" etaMax: "<<fMaxEtaCut<<endl;
  //cout<<"dEdxMin: "<<fTPCdEdxMin<<" dcaXY: "<<fDCAxyMax<<" dcaZ: "<<fDCAzzMax<<"  VzLow: "<<fMinVzCut<<" VzHigh: "<<fMaxVzCut<<endl;
  //cout<<"minPt: "<<fMinPtCut<<" maxPt: "<<fMaxPtCut<<" etaGapNeg: "<<fEtaGapNeg<<" etaGapPos: "<<fEtaGapPos<<" oldRun: "<<gOldRunNumber<<endl;
  //cout<<" Cent: "<<sCentrEstimator.Data()<<" EventPlane Det:"<<sDetectorForEP.Data()<<endl;
  

  

  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  Float_t centrCL1   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    exit(110);
  }

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 
 
  centrality = centrV0M;  // This Is Default, Flag for Estimator is set in AddTask.
  
  if(sCentrEstimator=="CL0"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL0");;
  }
  else if(sCentrEstimator=="CL1"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL1");;
  }
  else if(sCentrEstimator=="V0C"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0C");
  }
  else if(sCentrEstimator=="V0A"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0A");
  }
  else if(sCentrEstimator=="TRK"){
    centrality = fMultSelection->GetMultiplicityPercentile("TRK");
  }
  else{
    centrality = centrV0M;
  }


  fCentDistBeforCut->Fill(centrality);
  
  if(centrality<0 || centrality>90){ // Minimum,Maximum Centrality Cut
    return;
  }

  fDebugwEventCount->Fill(2.1);

  
  Int_t ntracks = fAOD->GetNumberOfTracks();
  if(ntracks < 4) return;    // Minimum 4 tracks per event for the Analysis!!

  fDebugwEventCount->Fill(3.1);
  
  
  /// ----> Get Magnetic field and RunNo.---------
  // Float_t fMagField = fAOD->GetMagneticField();
  Int_t runNumber = fAOD->GetRunNumber();

  //Set Up Correction Map for this run:
  if(runNumber!=gOldRunNumber) {
    if(fListNUACorr){
      GetNUACorrectionHist(runNumber,gParticleID); 
    }
    if(fListV0MCorr){
      GetV0MCorrectionHist(runNumber);
    }
    if(fListZDCCorr){
	  GetZDCCorrectionHist(runNumber);
	}
    gOldRunNumber = runNumber;
  }
  //----------------------------------------------

  UInt_t period = fAOD->GetPeriodNumber();
  UInt_t orbit24 = fAOD->GetOrbitNumber();
  
  if (period > 255) { // 8 bits
	  cout<<"invalid period number"<<endl;
	  period = 255;
	  orbit24 = (1<<24)-1;
  }
  
  
  if (orbit24 >= (1<<24)) { // 24 bits
	  cout<<"invalid orbit number"<<endl;
	  period = 255;
	  orbit24 = (1<<24)-1;
  }
  
  UInt_t orbit = period * (1<<24) + orbit24;

  Double_t fOrbitNumber = static_cast<double>(orbit)/1000000.; // scale down by 10^6 to fit to the scale used in least square fit. In least square fit, orbit number is scaled down by 10^6 to weight down the contribution to it. Maybe not necessary to scale down, but it should make the fit more stable in principle
  
  Bool_t kPileupEvent = kFALSE;

  //Double_t fSumQnxNeg = 0, fSumQnyNeg = 0, fSumQnxPos = 0, fSumQnyPos = 0;
  Double_t fSumQnxNeg[3] = {0,}; // Array: Q2x, Q3x, Q4x,...
  Double_t fSumQnyNeg[3] = {0,}; // Array: Q2y, Q3y, Q4y,... 
  Double_t fSumQnxPos[3] = {0,};
  Double_t fSumQnyPos[3] = {0,};  
  
  Double_t fMultNeg = 0, fMultPos = 0;

  vecPosEPTrkID.clear();
  vecNegEPTrkID.clear();
  vecPosEPTrkPhi.clear();
  vecNegEPTrkPhi.clear();
  vecPosEPTrkNUAWgt.clear();
  vecNegEPTrkNUAWgt.clear();

  vector<Int_t>().swap(vecPosEPTrkID);
  vector<Int_t>().swap(vecNegEPTrkID);   
  vector<Double_t>().swap(vecPosEPTrkPhi);
  vector<Double_t>().swap(vecNegEPTrkPhi);
  vector<Double_t>().swap(vecPosEPTrkNUAWgt);
  vector<Double_t>().swap(vecNegEPTrkNUAWgt);
  
  kPileupEvent = GetTPCQvectAndRemovePileUp2018(fAOD, fSumQnxNeg, fSumQnyNeg, fSumQnxPos, fSumQnyPos, fMultNeg, fMultPos);
  
  if(kPileupEvent) return;  // If not a PileUp event, then We have TPC q vectors for EP.
  fDebugwEventCount->Fill(4.1);

  
  if(fMultNeg<0.1 || fMultPos<0.1) return;  //// This means there is not enough track in the the event. 
  fDebugwEventCount->Fill(5.1);
  
  Double_t fQ2xNeg = 0., fQ2yNeg = 0., fQ2xPos = 0., fQ2yPos = 0.;  ///q = Q/Mult Vector from TPC;
  Double_t fQ3xNeg = 0., fQ3yNeg = 0., fQ3xPos = 0., fQ3yPos = 0.;  ///q = Q/Mult Vector from TPC;
  Double_t fQ4xNeg = 0., fQ4yNeg = 0., fQ4xPos = 0., fQ4yPos = 0.;  ///q = Q/Mult Vector from TPC;

  fQ2xNeg = fSumQnxNeg[0]/fMultNeg;
  fQ2yNeg = fSumQnyNeg[0]/fMultNeg;
  fQ3xNeg = fSumQnxNeg[1]/fMultNeg;
  fQ3yNeg = fSumQnyNeg[1]/fMultNeg;
  fQ4xNeg = fSumQnxNeg[2]/fMultNeg;
  fQ4yNeg = fSumQnyNeg[2]/fMultNeg;
  
  fQ2xPos = fSumQnxPos[0]/fMultPos;
  fQ2yPos = fSumQnyPos[0]/fMultPos;
  fQ3xPos = fSumQnxPos[1]/fMultPos;
  fQ3yPos = fSumQnyPos[1]/fMultPos;
  fQ4xPos = fSumQnxPos[2]/fMultPos;
  fQ4yPos = fSumQnyPos[2]/fMultPos;

  //cout<<"Before rec: q2xn "<<fQ2xNeg<<"\t q2yn "<<fQ2yNeg<<"\t q2xp "<<fQ2xPos<<"\t q2yp "<<fQ2yPos<<endl;

  // *** Rihan: Temporarily Turned off. Uncomment after end of Test.!
  ApplyTPCqVectRecenter(centrV0M, 2, fQ2xNeg, fQ2yNeg, fQ2xPos, fQ2yPos);
  ApplyTPCqVectRecenter(centrV0M, 3, fQ3xNeg, fQ3yNeg, fQ3xPos, fQ3yPos);  

  //cout<<"After  rec: q2xn "<<fQ2xNeg<<"\t q2yn "<<fQ2yNeg<<"\t q2xp "<<fQ2xPos<<"\t q2yp "<<fQ2yPos<<endl;
  //cout<<"------- Bug Testing Mode... we exit here....... "<<endl; return;



  
  /// Fill <Q> vector for TPC event plane :
  fAvgCos2PsivsCentEtaPos->Fill(centrV0M,fQ2xPos); 
  fAvgSin2PsivsCentEtaPos->Fill(centrV0M,fQ2yPos); 
  fAvgCos2PsivsCentEtaNeg->Fill(centrV0M,fQ2xNeg);
  fAvgSin2PsivsCentEtaNeg->Fill(centrV0M,fQ2yNeg);

  fAvgCos3PsivsCentEtaPos->Fill(centrV0M,fQ3xPos); 
  fAvgSin3PsivsCentEtaPos->Fill(centrV0M,fQ3yPos); 
  fAvgCos3PsivsCentEtaNeg->Fill(centrV0M,fQ3xNeg);
  fAvgSin3PsivsCentEtaNeg->Fill(centrV0M,fQ3yNeg);

  fAvgCos4PsivsCentEtaPos->Fill(centrV0M,fQ4xPos); 
  fAvgSin4PsivsCentEtaPos->Fill(centrV0M,fQ4yPos); 
  fAvgCos4PsivsCentEtaNeg->Fill(centrV0M,fQ4xNeg);
  fAvgSin4PsivsCentEtaNeg->Fill(centrV0M,fQ4yNeg);  

  
  
  /// TPC q2-vectors for ESE:
  fHistTPCPosqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQ2xPos*fQ2xPos + fQ2yPos*fQ2yPos)); // centrV0M hardcoded to prevent mismatch!
  fHistTPCNegqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQ2xNeg*fQ2xNeg + fQ2yNeg*fQ2yNeg));
  


  /// TPC Event Planes:
  Double_t fPsiNTPCPos = 0., fPsiNTPCNeg = 0.;
  Double_t fPsi3TPCPos = 0., fPsi3TPCNeg = 0.;
  Double_t fPsi4TPCPos = 0., fPsi4TPCNeg = 0.;
  
  Double_t fPsiNTPC = 0.;
  
  if(fQ2xPos != 0 && fQ2yPos != 0){
    fPsiNTPCPos = (1./2)*TMath::ATan2(fQ2yPos,fQ2xPos); // @Shi NTPC really is just 2TPC
    if(fPsiNTPCPos < 0) fPsiNTPCPos += TMath::TwoPi()/2;
    fPsi3TPCPos = (1./3)*TMath::ATan2(fQ3yPos,fQ3xPos);
    if(fPsi3TPCPos < 0) fPsi3TPCPos += TMath::TwoPi()/3;
    fPsi4TPCPos = (1./4)*TMath::ATan2(fQ4yPos,fQ4xPos);
    if(fPsi4TPCPos < 0) fPsi4TPCPos += TMath::TwoPi()/4;      
  }
  if(fQ2xNeg != 0 && fQ2yNeg != 0){
    fPsiNTPCNeg = (1./2)*TMath::ATan2(fQ2yNeg,fQ2xNeg);
    if(fPsiNTPCNeg < 0) fPsiNTPCNeg += TMath::TwoPi()/2;
    fPsi3TPCNeg = (1./3)*TMath::ATan2(fQ3yNeg,fQ3xNeg);
    if(fPsi3TPCNeg < 0) fPsi3TPCNeg += TMath::TwoPi()/3;
    fPsi4TPCNeg = (1./4)*TMath::ATan2(fQ4yNeg,fQ4xNeg);
    if(fPsi4TPCNeg < 0) fPsi4TPCNeg += TMath::TwoPi()/4;     
  }

  if(fQ2xPos != 0 && fQ2yPos != 0 && fQ2xNeg != 0 && fQ2yNeg != 0) {
	fPsiNTPC = (1./2)*TMath::ATan2(fQ2yPos+fQ2yNeg, fQ2xPos+fQ2xNeg);
	if(fPsiNTPC < 0) fPsiNTPC += TMath::TwoPi()/2;
  }
  
  fHistTPCPsiNPosPlane->Fill(centrality,fPsiNTPCPos);
  fHistTPCPsiNNegPlane->Fill(centrality,fPsiNTPCNeg);
  fHistTPCPsi3PosPlane->Fill(centrality,fPsi3TPCPos);
  fHistTPCPsi3NegPlane->Fill(centrality,fPsi3TPCNeg); 
  fHistTPCPsi4PosPlane->Fill(centrality,fPsi4TPCPos);
  fHistTPCPsi4NegPlane->Fill(centrality,fPsi4TPCNeg); 
  
  hTPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsiNTPCNeg));    /// TPC Psi2 Resolution
  hTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsi3TPCPos - 3*fPsi3TPCNeg));    /// TPC Psi3 Resolution
  hTPCPsi4Correlation->Fill(centrality,TMath::Cos(4*fPsi3TPCPos - 4*fPsi3TPCNeg));    /// TPC Psi3 Resolution



  
  
  ///==========>  Get V0 Event Planes <===============
     
  Double_t fQnxV0C=0, fQnyV0C=0, fQnxV0A=0, fQnyV0A=0; 
  Bool_t kPassV0 = GetGainCorrectedV0Qvector(fAOD, fVertexZEvent, 2, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A); // @Shi hard coded 2, fQnxV0 is really fQ2xV0

  Double_t fQ3xV0C=0, fQ3yV0C=0, fQ3xV0A=0, fQ3yV0A=0;
  kPassV0 = GetGainCorrectedV0Qvector(fAOD, fVertexZEvent, 3, fQ3xV0C, fQ3yV0C, fQ3xV0A, fQ3yV0A);

  
  if(!kPassV0) return;           /// V0 does not have signal for this event.  
  fDebugwEventCount->Fill(6.1);


  ApplyV0XqVectRecenter(centrCL1, 2, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A);

  
  ////---- fill the <Q> vector from V0A/C vs Cent:-----------
  hAvgQNXvsCentV0C->Fill(centrCL1,fQnxV0C);   /// to Avoid self correlation, V0 <Q> is filled with CL1 centrality! 
  hAvgQNYvsCentV0C->Fill(centrCL1,fQnyV0C);  
  hAvgQNXvsCentV0A->Fill(centrCL1,fQnxV0A);
  hAvgQNYvsCentV0A->Fill(centrCL1,fQnyV0A);

  //cout<<"------- Bug Testing Mode... we exit here....... "<<endl; return;

  
  ApplyV0XqVectRecenter(centrCL1, 3, fQ3xV0C, fQ3yV0C, fQ3xV0A, fQ3yV0A);
  hAvgQ3XvsCentV0C->Fill(centrCL1,fQ3xV0C);
  hAvgQ3YvsCentV0C->Fill(centrCL1,fQ3yV0C);  
  hAvgQ3XvsCentV0A->Fill(centrCL1,fQ3xV0A);
  hAvgQ3YvsCentV0A->Fill(centrCL1,fQ3yV0A);  


  /// V0 q2-vectors for ESE:
  fHistV0CDetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0C*fQnxV0C + fQnyV0C*fQnyV0C)); // centrCL1 hardcoded to prevent mismatch!
  fHistV0ADetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0A*fQnxV0A + fQnyV0A*fQnyV0A));
  

  
  ///--------> Get V0A and V0C Event Planes 

  Double_t fPsiNV0C = 0., fPsiNV0A = 0.;
  Double_t fPsi3V0C = 0., fPsi3V0A = 0.;
  Double_t fPsiNV0 = 0.;
  
  Bool_t kPassZNC = kFALSE;
    
  fPsiNV0C = (1./gHarmonic)*TMath::ATan2(fQnyV0C,fQnxV0C); // @Shi actually gHarmonic does not really make sense here, fQnxV0 is hard coded to be fQ2xV0
  if(fPsiNV0C < 0) fPsiNV0C += TMath::TwoPi()/gHarmonic;
  fPsiNV0A = (1./gHarmonic)*TMath::ATan2(fQnyV0A,fQnxV0A);
  if(fPsiNV0A < 0) fPsiNV0A += TMath::TwoPi()/gHarmonic;
  fHistV0CPsiNEventPlane->Fill(centrality,fPsiNV0C);
  fHistV0APsiNEventPlane->Fill(centrality,fPsiNV0A);

  fPsiNV0 = (1./gHarmonic)*TMath::ATan2(fQnyV0C+fQnyV0A,fQnxV0C+fQnxV0A);
  if(fPsiNV0 < 0) fPsiNV0 += TMath::TwoPi()/gHarmonic;
  
  fPsi3V0C = (1./3)*TMath::ATan2(fQ3yV0C,fQ3xV0C);
  if(fPsi3V0C < 0) fPsi3V0C += TMath::TwoPi()/3;      
  fPsi3V0A = (1./3)*TMath::ATan2(fQ3yV0A,fQ3xV0A);
  if(fPsi3V0A < 0) fPsi3V0A += TMath::TwoPi()/3;      
  fHistV0CPsi3EventPlane->Fill(centrality,fPsi3V0C);
  fHistV0APsi3EventPlane->Fill(centrality,fPsi3V0A);    

  /// V0A, V0C Resolutions:
  hV0CV0APsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNV0A    - 2*fPsiNV0C));
  hV0CTPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsiNV0C));
  hV0ATPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsiNV0A));

  hV0CV0APsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNV0A    - 3*fPsiNV0C));
  hV0CTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsiNV0C));
  hV0ATPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsiNV0A));

  ///--------> Get ZDC Event Planes
  Double_t fQxZNCC=0, fQyZNCC=0, fQxZNCA=0, fQyZNCA=0; 
  Double_t fPsiZNCC = 0., fPsiZNCA = 0., fPsiZNCCA = 0; // fPsiZNCCA combine ZNCC and ZNCA
  
  AliAODZDC *aodZDC = fAOD->GetZDCData();
	
  if(!aodZDC) {
    printf("\n ********* Error: could not find ZDC data ************ \n ");
  }
  else if(aodZDC && bUseZDCSpectatorPlane){
    const Double_t *fZNATowerRawAOD = aodZDC->GetZNATowerEnergy();
    const Double_t *fZNCTowerRawAOD = aodZDC->GetZNCTowerEnergy();
	
	if((fZNATowerRawAOD[0]<0) || (fZNATowerRawAOD[1]<0) || (fZNATowerRawAOD[2]<0) || (fZNATowerRawAOD[3]<0) || (fZNATowerRawAOD[4] < 0)) {
		bRecenterFailOrNot = kFALSE; 
		return;
	}
	
	if((fZNCTowerRawAOD[0]<0) || (fZNCTowerRawAOD[1]<0) || (fZNCTowerRawAOD[2]<0) || (fZNCTowerRawAOD[3]<0) || (fZNCTowerRawAOD[4] < 0)) {
		bRecenterFailOrNot = kFALSE; 
		return;
	}
	
	Double_t towZNCraw1GainEq = 0, towZNCraw2GainEq = 0, towZNCraw3GainEq = 0, towZNCraw4GainEq = 0;
	towZNCraw1GainEq = fZNCTowerRawAOD[1]*fHZDCCparameters->GetBinContent(1);
	towZNCraw2GainEq = fZNCTowerRawAOD[2]*fHZDCCparameters->GetBinContent(2);
	towZNCraw3GainEq = fZNCTowerRawAOD[3]*fHZDCCparameters->GetBinContent(3);
	towZNCraw4GainEq = fZNCTowerRawAOD[4]*fHZDCCparameters->GetBinContent(4);

	Double_t towZNAraw1GainEq = 0, towZNAraw2GainEq = 0, towZNAraw3GainEq = 0, towZNAraw4GainEq = 0;
	towZNAraw1GainEq = fZNATowerRawAOD[1]*fHZDCAparameters->GetBinContent(1);
	towZNAraw2GainEq = fZNATowerRawAOD[2]*fHZDCAparameters->GetBinContent(2);
	towZNAraw3GainEq = fZNATowerRawAOD[3]*fHZDCAparameters->GetBinContent(3);
	towZNAraw4GainEq = fZNATowerRawAOD[4]*fHZDCAparameters->GetBinContent(4);
	
	const Double_t xZDCC[4] = {-1, 1, -1, 1}; // directional vector
    const Double_t yZDCC[4] = {-1, -1, 1, 1};
    const Double_t xZDCA[4] = {1, -1, 1, -1};
    const Double_t yZDCA[4] = {-1, -1, 1, 1};
    
    Double_t towZNC[5] = {fZNCTowerRawAOD[0], towZNCraw1GainEq, towZNCraw2GainEq, towZNCraw3GainEq, towZNCraw4GainEq};
    Double_t towZNA[5] = {fZNATowerRawAOD[0], towZNAraw1GainEq, towZNAraw2GainEq, towZNAraw3GainEq, towZNAraw4GainEq};
    
	
    Double_t EZNC = 0, wZNC = 0, denZNC = 0, numXZNC = 0, numYZNC = 0;
    Double_t EZNA = 0, wZNA = 0, denZNA = 0, numXZNA = 0, numYZNA = 0; 

    for(Int_t i=0; i<4; i++){
		// ZNC part
        // get energy
        EZNC = towZNC[i+1];
        
        // build ZDCC centroid
        wZNC = TMath::Max(0., 4.0 + TMath::Log(towZNC[i+1]/fZNCTowerRawAOD[0]));
        numXZNC += xZDCC[i]*wZNC;
        numYZNC += yZDCC[i]*wZNC;
        denZNC += wZNC;
        
        // ZNA part
        // get energy
        EZNA = towZNA[i+1];

        // build ZDCA centroid
        wZNA = TMath::Max(0., 4.0 + TMath::Log(towZNA[i+1]/fZNATowerRawAOD[0]));
        numXZNA += xZDCA[i]*wZNA;
        numYZNA += yZDCA[i]*wZNA;
        denZNA += wZNA;
        
	}
	
	bRecenterFailOrNot = kTRUE;
	
	if (denZNC==0) {bRecenterFailOrNot = kFALSE; return;}
	if (denZNA==0) {bRecenterFailOrNot = kFALSE; return;}

	Double_t ZDCCxPosFromLogWeight = numXZNC/denZNC;
	Double_t ZDCCyPosFromLogWeight = numYZNC/denZNC;
	Double_t ZDCAxPosFromLogWeight = numXZNA/denZNA;
	Double_t ZDCAyPosFromLogWeight = numYZNA/denZNA;
	
	Double_t ZDCCAvgxPosFromVtxFit = 0;
	Double_t ZDCCAvgyPosFromVtxFit = 0;
	
	Double_t ZDCAAvgxPosFromVtxFit = 0;
	Double_t ZDCAAvgyPosFromVtxFit = 0;

	if (gTypeOfRecentering == 0) { // 1st order centrality+vtxPos, 
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pVtxX + fHZDCCparameters->GetBinContent(8)*pVtxY + fHZDCCparameters->GetBinContent(9)*pVtxZ + fHZDCCparameters->GetBinContent(10);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(11)*centrality + fHZDCCparameters->GetBinContent(12)*pVtxX + fHZDCCparameters->GetBinContent(13)*pVtxY + fHZDCCparameters->GetBinContent(14)*pVtxZ + fHZDCCparameters->GetBinContent(15);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pVtxX + fHZDCAparameters->GetBinContent(8)*pVtxY + fHZDCAparameters->GetBinContent(9)*pVtxZ + fHZDCAparameters->GetBinContent(10);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(11)*centrality + fHZDCAparameters->GetBinContent(12)*pVtxX + fHZDCAparameters->GetBinContent(13)*pVtxY + fHZDCAparameters->GetBinContent(14)*pVtxZ + fHZDCAparameters->GetBinContent(15);
	} else if (gTypeOfRecentering == 1) { // 3rd order centrality+vtxpos
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pow(centrality,2) + fHZDCCparameters->GetBinContent(8)*pow(centrality,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(13)*centrality + fHZDCCparameters->GetBinContent(14)*pow(centrality,2) + fHZDCCparameters->GetBinContent(15)*pow(centrality,3) + fHZDCCparameters->GetBinContent(16)*pVtxX + fHZDCCparameters->GetBinContent(17)*pVtxY + fHZDCCparameters->GetBinContent(18)*pVtxZ + fHZDCCparameters->GetBinContent(19);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pow(centrality,2) + fHZDCAparameters->GetBinContent(8)*pow(centrality,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(13)*centrality + fHZDCAparameters->GetBinContent(14)*pow(centrality,2) + fHZDCAparameters->GetBinContent(15)*pow(centrality,3) + fHZDCAparameters->GetBinContent(16)*pVtxX + fHZDCAparameters->GetBinContent(17)*pVtxY + fHZDCAparameters->GetBinContent(18)*pVtxZ + fHZDCAparameters->GetBinContent(19);
	} else if (gTypeOfRecentering == 2) { // 3rd order centrality+vtxpos+orbitNum
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pow(centrality,2) + fHZDCCparameters->GetBinContent(8)*pow(centrality,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12)*fOrbitNumber + fHZDCCparameters->GetBinContent(13);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(14)*centrality + fHZDCCparameters->GetBinContent(15)*pow(centrality,2) + fHZDCCparameters->GetBinContent(16)*pow(centrality,3) + fHZDCCparameters->GetBinContent(17)*pVtxX + fHZDCCparameters->GetBinContent(18)*pVtxY + fHZDCCparameters->GetBinContent(19)*pVtxZ + fHZDCCparameters->GetBinContent(20)*fOrbitNumber + fHZDCCparameters->GetBinContent(21);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pow(centrality,2) + fHZDCAparameters->GetBinContent(8)*pow(centrality,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12)*fOrbitNumber + fHZDCAparameters->GetBinContent(13);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(14)*centrality + fHZDCAparameters->GetBinContent(15)*pow(centrality,2) + fHZDCAparameters->GetBinContent(16)*pow(centrality,3) + fHZDCAparameters->GetBinContent(17)*pVtxX + fHZDCAparameters->GetBinContent(18)*pVtxY + fHZDCAparameters->GetBinContent(19)*pVtxZ + fHZDCAparameters->GetBinContent(20)*fOrbitNumber + fHZDCAparameters->GetBinContent(21);
	} else if (gTypeOfRecentering == 3) { // 1st order tow0 + vtxpos	
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*pVtxX + fHZDCCparameters->GetBinContent(7)*pVtxY + fHZDCCparameters->GetBinContent(8)*pVtxZ + fHZDCCparameters->GetBinContent(9)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(10);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(11)*pVtxX + fHZDCCparameters->GetBinContent(12)*pVtxY + fHZDCCparameters->GetBinContent(13)*pVtxZ + fHZDCCparameters->GetBinContent(14)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(15);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*pVtxX + fHZDCAparameters->GetBinContent(7)*pVtxY + fHZDCAparameters->GetBinContent(8)*pVtxZ + fHZDCAparameters->GetBinContent(9)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(10);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(11)*pVtxX + fHZDCAparameters->GetBinContent(12)*pVtxY + fHZDCAparameters->GetBinContent(13)*pVtxZ + fHZDCAparameters->GetBinContent(14)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(15);
	} else if (gTypeOfRecentering == 4) { // 3rd order tow0 + vtxpos
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*pVtxX + fHZDCCparameters->GetBinContent(7)*pVtxY + fHZDCCparameters->GetBinContent(8)*pVtxZ + fHZDCCparameters->GetBinContent(9)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(10)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(11)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(12);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(13)*pVtxX + fHZDCCparameters->GetBinContent(14)*pVtxY + fHZDCCparameters->GetBinContent(15)*pVtxZ + fHZDCCparameters->GetBinContent(16)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(17)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(18)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(19);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*pVtxX + fHZDCAparameters->GetBinContent(7)*pVtxY + fHZDCAparameters->GetBinContent(8)*pVtxZ + fHZDCAparameters->GetBinContent(9)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(10)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(11)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(12);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(13)*pVtxX + fHZDCAparameters->GetBinContent(14)*pVtxY + fHZDCAparameters->GetBinContent(15)*pVtxZ + fHZDCAparameters->GetBinContent(16)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(17)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(18)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(19);
	} else if (gTypeOfRecentering == 5) { // 5th order tow0 + vtxpos
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*pVtxX + fHZDCCparameters->GetBinContent(7)*pVtxY + fHZDCCparameters->GetBinContent(8)*pVtxZ + fHZDCCparameters->GetBinContent(9)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(10)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(11)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(12)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(13)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(14);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(15)*pVtxX + fHZDCCparameters->GetBinContent(16)*pVtxY + fHZDCCparameters->GetBinContent(17)*pVtxZ + fHZDCCparameters->GetBinContent(18)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(19)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(20)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(21)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(22)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(23);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*pVtxX + fHZDCAparameters->GetBinContent(7)*pVtxY + fHZDCAparameters->GetBinContent(8)*pVtxZ + fHZDCAparameters->GetBinContent(9)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(10)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(11)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(12)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(13)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(14);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(15)*pVtxX + fHZDCAparameters->GetBinContent(16)*pVtxY + fHZDCAparameters->GetBinContent(17)*pVtxZ + fHZDCAparameters->GetBinContent(18)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(19)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(20)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(21)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(22)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(23);
	} else if (gTypeOfRecentering == 6) { // 5th order tow0 + vtxpos + orbitNum
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*pVtxX + fHZDCCparameters->GetBinContent(7)*pVtxY + fHZDCCparameters->GetBinContent(8)*pVtxZ + fHZDCCparameters->GetBinContent(9)*fOrbitNumber + fHZDCCparameters->GetBinContent(10)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(11)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(12)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(13)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(14)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(15);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(16)*pVtxX + fHZDCCparameters->GetBinContent(17)*pVtxY + fHZDCCparameters->GetBinContent(18)*pVtxZ + fHZDCCparameters->GetBinContent(19)*fOrbitNumber + fHZDCCparameters->GetBinContent(20)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(21)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(22)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(23)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(24)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(25);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*pVtxX + fHZDCAparameters->GetBinContent(7)*pVtxY + fHZDCAparameters->GetBinContent(8)*pVtxZ + fHZDCAparameters->GetBinContent(9)*fOrbitNumber + fHZDCAparameters->GetBinContent(10)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(11)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(12)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(13)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(14)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(15);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(16)*pVtxX + fHZDCAparameters->GetBinContent(17)*pVtxY + fHZDCAparameters->GetBinContent(18)*pVtxZ + fHZDCAparameters->GetBinContent(19)*fOrbitNumber + fHZDCAparameters->GetBinContent(20)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(21)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(22)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(23)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(24)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(25);
	} else if (gTypeOfRecentering == 7) { // 5th order tow0 + 3rd order centrality + vtxpos + orbitNum
		ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pow(centrality,2) + fHZDCCparameters->GetBinContent(8)*pow(centrality,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12)*fOrbitNumber + fHZDCCparameters->GetBinContent(13)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(14)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(15)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(16)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(17)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(18);
		ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(19)*centrality + fHZDCCparameters->GetBinContent(20)*pow(centrality,2) + fHZDCCparameters->GetBinContent(21)*pow(centrality,3) + fHZDCCparameters->GetBinContent(22)*pVtxX + fHZDCCparameters->GetBinContent(23)*pVtxY + fHZDCCparameters->GetBinContent(24)*pVtxZ + fHZDCCparameters->GetBinContent(25)*fOrbitNumber + fHZDCCparameters->GetBinContent(26)*fZNCTowerRawAOD[0]/100000 + fHZDCCparameters->GetBinContent(27)*pow(fZNCTowerRawAOD[0]/100000,2) + fHZDCCparameters->GetBinContent(28)*pow(fZNCTowerRawAOD[0]/100000,3) + fHZDCCparameters->GetBinContent(29)*pow(fZNCTowerRawAOD[0]/100000,4) + fHZDCCparameters->GetBinContent(30)*pow(fZNCTowerRawAOD[0]/100000,5) + fHZDCCparameters->GetBinContent(31);
		
		ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pow(centrality,2) + fHZDCAparameters->GetBinContent(8)*pow(centrality,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12)*fOrbitNumber + fHZDCAparameters->GetBinContent(13)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(14)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(15)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(16)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(17)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(18);
		ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(19)*centrality + fHZDCAparameters->GetBinContent(20)*pow(centrality,2) + fHZDCAparameters->GetBinContent(21)*pow(centrality,3) + fHZDCAparameters->GetBinContent(22)*pVtxX + fHZDCAparameters->GetBinContent(23)*pVtxY + fHZDCAparameters->GetBinContent(24)*pVtxZ + fHZDCAparameters->GetBinContent(25)*fOrbitNumber + fHZDCAparameters->GetBinContent(26)*fZNATowerRawAOD[0]/100000 + fHZDCAparameters->GetBinContent(27)*pow(fZNATowerRawAOD[0]/100000,2) + fHZDCAparameters->GetBinContent(28)*pow(fZNATowerRawAOD[0]/100000,3) + fHZDCAparameters->GetBinContent(29)*pow(fZNATowerRawAOD[0]/100000,4) + fHZDCAparameters->GetBinContent(30)*pow(fZNATowerRawAOD[0]/100000,5) + fHZDCAparameters->GetBinContent(31);
	}
	
	fQxZNCC = ZDCCxPosFromLogWeight - ZDCCAvgxPosFromVtxFit;
	fQyZNCC = ZDCCyPosFromLogWeight - ZDCCAvgyPosFromVtxFit;
	
	fQxZNCA = ZDCAxPosFromLogWeight - ZDCAAvgxPosFromVtxFit;
	fQyZNCA = ZDCAyPosFromLogWeight - ZDCAAvgyPosFromVtxFit;
	
	// Event plane

	fPsiZNCA = TMath::ATan2(fQyZNCA,fQxZNCA); // Psi_{1,A} spectator plane -pi to pi
	if (fPsiZNCA < 0) { // Psi_{1,A} should be differ to Psi_{1,C} by pi. 
	  fPsiZNCA = fPsiZNCA + TMath::Pi();
	} else if (fPsiZNCA >= 0) {
	  fPsiZNCA = fPsiZNCA - TMath::Pi();
	}

	fPsiZNCC = TMath::ATan2(fQyZNCC,fQxZNCC); // Psi_{1,C} spectator plane 

	fPsiZNCCA = TMath::ATan2((fQyZNCC-fQyZNCA),(fQxZNCC-fQxZNCA));

	// QAs for ZDC 
	hZDCCPsiSpectatorPlane->Fill(fPsiZNCC);
	hZDCAPsiSpectatorPlane->Fill(fPsiZNCA);
	hZDCCAPsiSpectatorPlane->Fill(fPsiZNCCA);
	
	Double_t fPsiZNCCProj = fPsiZNCC;
	Double_t fPsiZNCAProj = fPsiZNCA;
	Double_t fPsiZNCCAProj = fPsiZNCCA;
	if (fPsiZNCCProj<0) fPsiZNCCProj += TMath::TwoPi()/2;
	if (fPsiZNCAProj<0) fPsiZNCAProj += TMath::TwoPi()/2;
	if (fPsiZNCCAProj<0) fPsiZNCCAProj += TMath::TwoPi()/2;
	
	hZDCCPsivsZDCAPsi->Fill(fPsiZNCC, fPsiZNCA);
	hZDCCPsivsTPCPsi2->Fill(fPsiZNCC, fPsiNTPC);
	hZDCCPsivsTPCPosPsi2->Fill(fPsiZNCC, fPsiNTPCPos);
	hZDCCPsivsTPCNegPsi2->Fill(fPsiZNCC, fPsiNTPCNeg);
	hZDCCPsivsV0CAPsi2->Fill(fPsiZNCC, fPsiNV0);
	hZDCCPsivsV0CPsi2->Fill(fPsiZNCC, fPsiNV0C);
	hZDCCPsivsV0APsi2->Fill(fPsiZNCC, fPsiNV0A);
	
	hZDCAPsivsTPCPsi2->Fill(fPsiZNCA, fPsiNTPC);
	hZDCAPsivsTPCPosPsi2->Fill(fPsiZNCA, fPsiNTPCPos);
	hZDCAPsivsTPCNegPsi2->Fill(fPsiZNCA, fPsiNTPCNeg);
	hZDCAPsivsV0CAPsi2->Fill(fPsiZNCA, fPsiNV0);
	hZDCAPsivsV0CPsi2->Fill(fPsiZNCA, fPsiNV0C);
	hZDCAPsivsV0APsi2->Fill(fPsiZNCA, fPsiNV0A);
	
	hZDCCAPsivsTPCPsi2->Fill(fPsiZNCCA, fPsiNTPC);
	hZDCCAPsivsTPCPosPsi2->Fill(fPsiZNCCA, fPsiNTPCPos);
	hZDCCAPsivsTPCNegPsi2->Fill(fPsiZNCCA, fPsiNTPCNeg);
	hZDCCAPsivsV0CAPsi2->Fill(fPsiZNCCA, fPsiNV0);
	hZDCCAPsivsV0CPsi2->Fill(fPsiZNCCA, fPsiNV0C);
	hZDCCAPsivsV0APsi2->Fill(fPsiZNCCA, fPsiNV0A);
  }
  
  

  //If we skip the nested loops and main Analysis:
  if(bSkipNestedLoop){
 
    fHistVertexZcm->Fill(pVtxZ);
    fCentDistAfterCut->Fill(centrality);
    
    PostData(1,fListHist);
    return;     //Just fill QAs and get out..
  }

  Double_t fSelectedV0PsiN = 0;
  Double_t fSelectedV0Psi3 = 0;
  Double_t fSelectedTPCPsiN = 0;
  Double_t fSelectedTPCPsi3 = 0;
  Int_t gEPeta = 0;

  if(sDetectorForEP.Contains("TPCPos")||sDetectorForEP.Contains("TPCpos")){
    gEPeta = 1;
    bUseV0EventPlane = kFALSE;
    fSelectedTPCPsiN  = fPsiNTPCPos;
    fSelectedTPCPsi3  = fPsi3TPCPos;    
  }
  else if(sDetectorForEP.Contains("TPCNeg")||sDetectorForEP.Contains("TPCneg")){
    gEPeta = -1;
    bUseV0EventPlane = kFALSE;
    fSelectedTPCPsiN  = fPsiNTPCNeg;
    fSelectedTPCPsi3  = fPsi3TPCNeg;
  }
  else if(sDetectorForEP.Contains("TPC")){
    gEPeta = 1;
    bUseV0EventPlane = kFALSE;
    fSelectedTPCPsiN  = fPsiNTPCPos;
    fSelectedTPCPsi3  = fPsi3TPCPos;
  }
  else if(sDetectorForEP.Contains("V0C")){
    bUseV0EventPlane = kTRUE;
    fSelectedV0PsiN  = fPsiNV0C;
    fSelectedV0Psi3  = fPsi3V0C;
  }
  else if(sDetectorForEP.Contains("V0A")){
    bUseV0EventPlane = kTRUE;
    fSelectedV0PsiN  = fPsiNV0A;
    fSelectedV0Psi3  = fPsi3V0A;
  }
  else if(sDetectorForEP.Contains("V0")){
    bUseV0EventPlane = kTRUE;
    fSelectedV0PsiN  = fPsiNV0A;
    fSelectedV0Psi3  = fPsi3V0A;
  }
  else{ 
    bUseV0EventPlane = kFALSE;
  }
  







  

  /// We Are About to Start Main Analysis Below: //
  
 
  Int_t kPIDtrk1=gParticleID;   /// gParticleID is Set From AddTask. Both Identified..
  Int_t kPIDtrk2=gParticleID;   /// 0 = hadron (h-h), 1 = Pi-Pi, 2 = K-K, 3 = Prot-Prot, 

  ///For single Identified cases:
  
  if(gParticleID==10){
    kPIDtrk1 = 0; //Ch
    kPIDtrk2 = 1; //Pion
  }
  else if(gParticleID==20){
    kPIDtrk1 = 0; //Ch
    kPIDtrk2 = 2; //Kaon
  }
  else if(gParticleID==30){
    kPIDtrk1 = 0; //Ch
    kPIDtrk2 = 3; //P
  }


///Track variables:
  Int_t   trk1Chrg=0, trk1TpcNC=0;
  Int_t   trk2Chrg=0, trk2TpcNC=0;
  Bool_t  bPIDoktrk1=kFALSE, bPIDoktrk2=kFALSE;

  
  Double_t trk1Pt=0,trk1Phi=0,trk1Eta=0,trk1DCAxy=0.0, trk1DCAz=0.0,trk1Chi2=0,trk1dEdx=0,trk1Wgt=1.0;
  Double_t trk2Pt=0,trk2Phi=0,trk2Eta=0,trk2DCAxy=0.0, trk2DCAz=0.0,trk2Chi2=0,trk2dEdx=0,trk2Wgt=1.0;  
  Double_t wgtComb1Ch = 1.0, wgtComb1PID = 1.0;
  Double_t wgtComb2part = 1.0; 
  Double_t ptWgtMCChtrk1  = 1.0, WgtNUAChtrk1   = 1.0;
  Double_t ptWgtMCChtrk2  = 1.0, WgtNUAChtrk2   = 1.0;
  Double_t ptWgtMCPIDtrk1 = 1.0, WgtNUAPIDtrk1  = 1.0;
  Double_t ptWgtMCPIDtrk2 = 1.0, WgtNUAPIDtrk2  = 1.0;
  Double_t wgt1PIDparticle= 1.0, wgt2PIDparticle= 1.0; 

  Double_t localSumQ2x =0,localSumQ2y=0;
  Double_t localSumQ3x =0,localSumQ3y=0;
  Double_t localSumQ2xs =0,localSumQ2ys=0;
  Double_t localSumQ3xs =0,localSumQ3ys=0;  
  
  Double_t localMultTPC=0,localMultTPCs=0; 


  // Chun Zheng: vectors to contain information for Lambda-x pairing..
  vector<Double_t> vecPt;
  vector<Int_t>    vecID;
  vector<Double_t> vecEta;
  vector<Double_t> vecPhi;
  vector<Int_t>    vecPDGCode;
  vector<Double_t> vecNUAWeight;     //Charge
  vector<Double_t> vecNUEWeight;     //Charge
  vector<Double_t> vecNUAWeightPID;  //Charge
  vector<Double_t> vecNUEWeightPID;  //Charge
  Bool_t isItPiontrk1 = kFALSE, isItKaontrk1 = kFALSE, isItProttrk1 = kFALSE;

  
  ///----------> Starting Analysis track Loop -----------
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { 

    AliAODTrack* AODtrack1 = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(iTrack));
    if(!AODtrack1) continue;
    
    if(AODtrack1->TestFilterBit(fFilterBit)) {  //// Only use FB tracks. 

      trk1Pt    = AODtrack1->Pt();
      trk1Phi   = AODtrack1->Phi();
      trk1Eta   = AODtrack1->Eta();
      trk1Chrg  = AODtrack1->Charge();
      trk1Chi2  = AODtrack1->Chi2perNDF();
      trk1TpcNC = AODtrack1->GetTPCNcls();
      trk1DCAxy = AODtrack1->DCA();
      trk1DCAz  = AODtrack1->ZAtDCA(); 
      trk1dEdx  = AODtrack1->GetDetPid()->GetTPCsignal();  
     
      //Apply track cuts here:
      if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && (trk1dEdx >= fTPCdEdxMin) && (trk1TpcNC >= fTPCclustMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fTrkChi2Max) && TMath::Abs(trk1Chrg)) {

	WgtNUAChtrk1  = 1.0;   
	WgtNUAPIDtrk1 = 1.0;
	ptWgtMCChtrk1 = 1.0;  
	ptWgtMCPIDtrk1 = 1.0;


	WgtNUAChtrk1  = GetNUAWeightForTrack(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg);
	WgtNUAPIDtrk1 = GetNUAWeightForTrackPID(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg);
	ptWgtMCChtrk1 = GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,0); // this is eff
	ptWgtMCChtrk1 = 1./ptWgtMCChtrk1; // weight = 1/eff
	ptWgtMCPIDtrk1= GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,kPIDtrk1);
	ptWgtMCPIDtrk1 = 1./ptWgtMCPIDtrk1;

	wgtComb1Ch  = WgtNUAChtrk1*ptWgtMCChtrk1;    /// Charge
	wgtComb1PID = WgtNUAPIDtrk1*ptWgtMCPIDtrk1;  /// PID


	
	//----- For Chunzheng ---------->
	/// fill the information in vector array for Lambda-X pair
       
	if(bAnalysLambdaPairs){

	  Int_t code = 0;

	  isItProttrk1 = kFALSE; 

	  isItProttrk1 = CheckPIDofParticle(AODtrack1,3); // 3=proton
	
	  if(trk1Chrg > 0) {
      code = 999;
	    if(isItProttrk1) code = 2212;
	  } else{  /// 
      code = -999;
	    if(isItProttrk1) code = -2212;
	  }
	  
	  Int_t trk1ID = AODtrack1->GetID();//unique in a event
	  vecPDGCode.push_back(code);
	  vecPhi.push_back(trk1Phi);
	  vecEta.push_back(trk1Eta);
	  vecPt.push_back(trk1Pt);
	  vecID.push_back(trk1ID);
	  vecNUAWeight.push_back(WgtNUAChtrk1);  //if Particle is considered Unidentified
	  vecNUEWeight.push_back(ptWgtMCChtrk1); //if Particle is considered Unidentified
	  vecNUAWeightPID.push_back(WgtNUAPIDtrk1);  //if Particle is considered Identified (good purity upto pT<2 GeV)
	  vecNUEWeightPID.push_back(ptWgtMCPIDtrk1); //if Particle is considered Identified (good purity upto pT<2 GeV)

	  continue; /// Skip Analysing CME PID if Lambda-X study is running. Otherwise Jobs may reach TTL and killed.	  
	}

	
	//<---- Chunzheng Filled his Vectors for single Particle id and weights.
	//<---- Rest of the Lambda Analysis is done outside this track loop.




	
	
	/// Rihan: The part below is only relevant for CME Analysis (no Lambda):
	
	bPIDoktrk1=kFALSE;
	bPIDoktrk2=kFALSE;
	
	bPIDoktrk1 = CheckPIDofParticle(AODtrack1,kPIDtrk1);   // check if track1 is of desired PID request #1,

	if(!bPIDoktrk1)
	  bPIDoktrk2 = CheckPIDofParticle(AODtrack1,kPIDtrk2);  // if not request #1, then check if track satisfies request #2

	if(!bPIDoktrk1 && !bPIDoktrk2)    /// the track-1 is neither of the desired PIDs, then skip this track. 
	  continue;

       
	Double_t fPsiNEvent = 0, fPsi3Event = 0;

	///Choose whether to use TPC or V0EP
	if(bUseV0EventPlane){
	  fPsiNEvent = fSelectedV0PsiN;
	  fPsi3Event = fSelectedV0Psi3;
	}
	else{ 
	  fPsiNEvent = fSelectedTPCPsiN;
	  fPsi3Event = fSelectedTPCPsi3;	    
	}

	if(gEPeta < 0){
	  localSumQ2x = fSumQnxNeg[0];            /// We need the full Q-sum. Then remove only qx,qy for current track
	  localSumQ2y = fSumQnyNeg[0];
	  localSumQ3x = fSumQnxNeg[1];
	  localSumQ3y = fSumQnyNeg[1];
	  localMultTPC= fMultNeg;
	}
	else{
	  localSumQ2x = fSumQnxPos[0];
	  localSumQ2y = fSumQnyPos[0];
	  localSumQ3x = fSumQnxPos[1];
	  localSumQ3y = fSumQnyPos[1];
	  localMultTPC= fMultPos;
	}
	
	///remove Autocorrelation for track-1 only if both EP and track1 are on the same eta side.
	if(!bUseV0EventPlane && trk1Pt < 2.0 && (gEPeta*trk1Eta) > 0){   /// we used pT<2.0 tracks for EP.
	  if(gEPeta > 0 && trk1Eta >= 0.1) { // exclude tracks falling into eta gap
	    localSumQ2x -= WgtNUAChtrk1*TMath::Cos(2*trk1Phi);   /// wgts and phi of track1
	    localSumQ2y -= WgtNUAChtrk1*TMath::Sin(2*trk1Phi);
	    localSumQ3x -= WgtNUAChtrk1*TMath::Cos(3*trk1Phi);
	    localSumQ3y -= WgtNUAChtrk1*TMath::Sin(3*trk1Phi);
	    localMultTPC-= WgtNUAChtrk1; 	  
      }
      else if(gEPeta < 0 && trk1Eta <= -0.1) {
	    localSumQ2x -= WgtNUAChtrk1*TMath::Cos(2*trk1Phi);   /// wgts and phi of track1
	    localSumQ2y -= WgtNUAChtrk1*TMath::Sin(2*trk1Phi);
	    localSumQ3x -= WgtNUAChtrk1*TMath::Cos(3*trk1Phi);
	    localSumQ3y -= WgtNUAChtrk1*TMath::Sin(3*trk1Phi);
	    localMultTPC-= WgtNUAChtrk1; 	  
      }
	}
	
    
	// save v2 for each plane
	if(localMultTPC>0){
		Double_t localSumQ2xTemp = localSumQ2x/localMultTPC;
		Double_t localSumQ2yTemp = localSumQ2y/localMultTPC;
	
		Double_t fPsi2TPC = 0;
		
		fPsi2TPC = (1./2)*TMath::ATan2(localSumQ2y,localSumQ2xTemp);
		if(fPsi2TPC < 0) fPsi2TPC += TMath::TwoPi()/2;
		
		if(gParticleID==0){
		  hV2TPCvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsi2TPC), wgtComb1Ch);
		  hV2V0CvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNV0C),wgtComb1Ch);
		  hV2V0AvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNV0A),wgtComb1Ch);
		  hV2V0CAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - (fPsiNV0C+fPsiNV0A)),wgtComb1Ch);
		  hV2ZDCCvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCC),wgtComb1Ch);
		  hV2ZDCAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCA),wgtComb1Ch);
		  hV2ZDCCACombinevsCent->Fill(centrality,TMath::Cos(2*trk1Phi - (fPsiZNCC+fPsiZNCA)),wgtComb1Ch);
		  hV2ZDCCAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCCA),wgtComb1Ch);
		  
		  // for four particle cumulant 
		  hPOIQRe->Fill(0.5,wgtComb1Ch*TMath::Cos(2*trk1Phi)); // w*cos(2phi)
		  hPOIQRe->Fill(1.5,pow(wgtComb1Ch,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
		  hPOIQRe->Fill(2.5,pow(wgtComb1Ch,3)*TMath::Cos(2*trk1Phi)); // w^3*cos(2phi)
		  hPOIQIm->Fill(0.5,wgtComb1Ch*TMath::Sin(2*trk1Phi)); // w*sin(2phi)
		  hPOIQIm->Fill(1.5,pow(wgtComb1Ch,2)*TMath::Sin(4*trk1Phi)); // w^2*sin(4phi)
		  hPOIQIm->Fill(2.5,pow(wgtComb1Ch,3)*TMath::Sin(2*trk1Phi)); // w^3*sin(2phi)
		  hPOIMul->Fill(0.5,pow(wgtComb1Ch,0)); // 1
		  hPOIMul->Fill(1.5,pow(wgtComb1Ch,1)); // w
		  hPOIMul->Fill(2.5,pow(wgtComb1Ch,2)); // w^2
		  hPOIMul->Fill(3.5,pow(wgtComb1Ch,3)); // w^3
		  hPOIMul->Fill(4.5,pow(wgtComb1Ch,4)); // w^4
		}
		else if(gParticleID < 10){/// Both Particles are Identified. Not considering gParticleID>10
		  hV2TPCvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsi2TPC), wgtComb1PID);
		  hV2V0CvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNV0C),wgtComb1PID);
		  hV2V0AvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNV0A),wgtComb1PID);
		  hV2V0CAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - (fPsiNV0C+fPsiNV0A)),wgtComb1PID);
		  hV2ZDCCvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCC),wgtComb1PID);
		  hV2ZDCAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCA),wgtComb1PID);
		  hV2ZDCCACombinevsCent->Fill(centrality,TMath::Cos(2*trk1Phi - (fPsiZNCC+fPsiZNCA)),wgtComb1PID);
		  hV2ZDCCAvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiZNCCA),wgtComb1PID);
		  
		  // for four particle cumulant 
		  hPOIQRe->Fill(0.5,wgtComb1PID*TMath::Cos(2*trk1Phi)); // w*cos(2phi)
		  hPOIQRe->Fill(1.5,pow(wgtComb1PID,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
		  hPOIQRe->Fill(2.5,pow(wgtComb1PID,3)*TMath::Cos(2*trk1Phi)); // w^3*cos(2phi)
		  hPOIQIm->Fill(0.5,wgtComb1PID*TMath::Sin(2*trk1Phi)); // w*sin(2phi)
		  hPOIQIm->Fill(1.5,pow(wgtComb1PID,2)*TMath::Sin(4*trk1Phi)); // w^2*sin(4phi)
		  hPOIQIm->Fill(2.5,pow(wgtComb1PID,3)*TMath::Sin(2*trk1Phi)); // w^3*sin(2phi)
		  hPOIMul->Fill(0.5,pow(wgtComb1PID,0)); // 1
		  hPOIMul->Fill(1.5,pow(wgtComb1PID,1)); // w
		  hPOIMul->Fill(2.5,pow(wgtComb1PID,2)); // w^2
		  hPOIMul->Fill(3.5,pow(wgtComb1PID,3)); // w^3
		  hPOIMul->Fill(4.5,pow(wgtComb1PID,4)); // w^4
		  
		}
	}
	
	///---> 2nd track Loop   
	for(Int_t jTrack = 0; jTrack < ntracks; jTrack++) { 

	  /// skip autocorrelation:
	  if(jTrack==iTrack) continue;

	  
	  AliAODTrack* AODtrack2 = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(jTrack));
	  if(!AODtrack2) continue;
    
	  if(AODtrack2->TestFilterBit(fFilterBit)) {  //// Only use FB tracks. 

	    trk2Pt    = AODtrack2->Pt();
	    trk2Phi   = AODtrack2->Phi();
	    trk2Eta   = AODtrack2->Eta();
	    trk2Chrg  = AODtrack2->Charge();
	    trk2Chi2  = AODtrack2->Chi2perNDF();
	    trk2TpcNC = AODtrack2->GetTPCNcls();
	    trk2DCAxy = AODtrack2->DCA();
	    trk2DCAz  = AODtrack2->ZAtDCA();            
	    trk2dEdx  = AODtrack2->GetDetPid()->GetTPCsignal();       

	    //Apply track cuts for second track
	    if((trk2Pt <= fMaxPtCut) && (trk2Pt >= fMinPtCut) && (trk2Eta <= fMaxEtaCut) && (trk2Eta >= fMinEtaCut) && (trk2dEdx >= fTPCdEdxMin) && (trk2TpcNC >= fTPCclustMin) && (trk2Chi2 >= fTrkChi2Min) && (trk2Chi2 <= fTrkChi2Max) && TMath::Abs(trk2Chrg)) {


	      if(bPIDoktrk1){
		bPIDoktrk2 = CheckPIDofParticle(AODtrack2,kPIDtrk2); //
	      }
	      else if(bPIDoktrk2){
		bPIDoktrk1 = CheckPIDofParticle(AODtrack2,kPIDtrk1); 
	      }

	      if(!bPIDoktrk1 && !bPIDoktrk2)    /// the track1 and track2 is None of the desired PID, skip this pair..
		continue;

	      //// If I am here then I have both partners. Lets fill correlator Histograms:
	      

	      wgtComb2part = 1.0;
	      
	      if(gParticleID==0){ /// Both Particles are UN-Identified
		WgtNUAChtrk2  = GetNUAWeightForTrack(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg);
		ptWgtMCChtrk2 = GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,0);		
		ptWgtMCChtrk2 = 1./ptWgtMCChtrk2;
		wgtComb2part  = WgtNUAChtrk1*ptWgtMCChtrk1 * WgtNUAChtrk2*ptWgtMCChtrk2;      /// Combined weight for Ch
	      }
	      else if(gParticleID < 10){/// Both Particles are Identified.			
		WgtNUAPIDtrk2  = GetNUAWeightForTrackPID(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg);      
		ptWgtMCPIDtrk2 = GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,kPIDtrk2);
		ptWgtMCPIDtrk2 = 1./ptWgtMCPIDtrk2;
		wgtComb2part  = WgtNUAPIDtrk1*ptWgtMCPIDtrk1 * WgtNUAPIDtrk2*ptWgtMCPIDtrk2;  /// Combined weight for PID
	      }
	      else if(gParticleID > 10){/// Single Particle is Identified.
		WgtNUAChtrk2   = GetNUAWeightForTrack(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg);    // first track is Charged
		ptWgtMCChtrk2  = GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,0);                // first track is Charged
		ptWgtMCChtrk2 = 1./ptWgtMCChtrk2;
		WgtNUAPIDtrk2  = GetNUAWeightForTrackPID(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg); // 2nd is PID
		ptWgtMCPIDtrk2 = GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,kPIDtrk2);         // 2nd is PID
		ptWgtMCPIDtrk2 = 1./ptWgtMCPIDtrk2;
		wgtComb2part  = WgtNUAChtrk1*ptWgtMCChtrk1 * WgtNUAPIDtrk2*ptWgtMCPIDtrk2;        // Wgt for Hybrid Ch+PID Analysis.
	      }
	   


	      ///remove Autocorrelation for track-2 only if both EP and track2 are on the same eta side.
	      
	    localSumQ2xs = localSumQ2x;                          /// We need the the Qsum (with first track q removed), For Each track-2 
		localSumQ2ys = localSumQ2y;                          /// Otherwise we reduce the sum to zero in this 2nd loop.!!
		localSumQ3xs = localSumQ3x;
		localSumQ3ys = localSumQ3y;
		localMultTPCs = localMultTPC;
		
	    if(!bUseV0EventPlane && trk2Pt < 2.0 && (gEPeta*trk2Eta) > 0){   /// We used pT<2.0 tracks for EP.
		  if(gEPeta > 0 && trk2Eta >= 0.1) {
		    localSumQ2xs -= WgtNUAChtrk2*TMath::Cos(2*trk2Phi);  /// wgts and phi of track2
		    localSumQ2ys -= WgtNUAChtrk2*TMath::Sin(2*trk2Phi);
		    localSumQ3xs -= WgtNUAChtrk2*TMath::Cos(3*trk2Phi);
		    localSumQ3ys -= WgtNUAChtrk2*TMath::Sin(3*trk2Phi);
		    localMultTPCs-= WgtNUAChtrk2;
		  }
		  else if(gEPeta < 0 && trk2Eta <= -0.1) {
		    localSumQ2xs -= WgtNUAChtrk2*TMath::Cos(2*trk2Phi);  /// wgts and phi of track2
		    localSumQ2ys -= WgtNUAChtrk2*TMath::Sin(2*trk2Phi);
		    localSumQ3xs -= WgtNUAChtrk2*TMath::Cos(3*trk2Phi);
		    localSumQ3ys -= WgtNUAChtrk2*TMath::Sin(3*trk2Phi);
		    localMultTPCs-= WgtNUAChtrk2;
		  }
	    }

	    if(!bUseV0EventPlane){
 		
		  if(localMultTPCs>0){
		    localSumQ2xs = localSumQ2xs/localMultTPCs;
		    localSumQ2ys = localSumQ2ys/localMultTPCs;
		    localSumQ3xs = localSumQ3xs/localMultTPCs;
		    localSumQ3ys = localSumQ3ys/localMultTPCs;
		  }
		
		  fPsiNEvent = (1./2)*TMath::ATan2(localSumQ2ys,localSumQ2xs);
		  if(fPsiNEvent < 0) fPsiNEvent += TMath::TwoPi()/2;		
	  	  fPsi3Event = (1./3)*TMath::ATan2(localSumQ3ys,localSumQ3xs);
		  if(fPsi3Event < 0) fPsi3Event += TMath::TwoPi()/3;
		  //fPsi4TPCPos = (1./4)*TMath::ATan2(fQ4yPos,fQ4xPos);
		  //if(fPsi4TPCPos < 0) fPsi4TPCPos += TMath::TwoPi()/4;  		
	    }
	      
	      if(trk1Chrg*trk2Chrg < 0){ //Opposite sign	
		hAvg3pC112vsCentOS->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb2part);
		hAvg3pC123vsCentOS->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb2part);
				
		hAvgDelta1vsCentOS->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb2part);
		hAvgDelta2vsCentOS->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta3vsCentOS->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta4vsCentOS->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb2part);		
	      }		
	      else if(trk1Chrg > 0 && trk2Chrg > 0){ ///pos-pos       
		hAvg3pC112vsCentPP->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb2part);
		hAvg3pC123vsCentPP->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb2part);
				
		hAvgDelta1vsCentPP->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb2part);
		hAvgDelta2vsCentPP->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta3vsCentPP->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta4vsCentPP->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb2part);		
	      }
	      //else if(trk1Chrg < 0 && trk2Chrg < 0){  ///this is obviously last option (and every microseconds counts!)
	      else{
		hAvg3pC112vsCentNN->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb2part);
		hAvg3pC123vsCentNN->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb2part);
		
		hAvgDelta1vsCentNN->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb2part);
		hAvgDelta2vsCentNN->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta3vsCentNN->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb2part);
		hAvgDelta4vsCentNN->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb2part);		
	      }
	      
	      
	      // ZDC spectator plane part
	      if(bUseZDCSpectatorPlane) {
			  if(trk1Chrg*trk2Chrg < 0){ //Opposite sign	
			hAvg3pC112vsCentOSUsingZDCC->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCC),wgtComb2part);
			hAvg3pC112vsCentOSUsingZDCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCA),wgtComb2part);
			hAvg3pC112vsCentOSUsingZDCCACombine->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - (fPsiZNCC+fPsiZNCA)),wgtComb2part);
			hAvg3pC112vsCentOSUsingZDCCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCCA),wgtComb2part);
					
			  }		
			  else if(trk1Chrg > 0 && trk2Chrg > 0){ ///pos-pos       
			hAvg3pC112vsCentPPUsingZDCC->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCC),wgtComb2part);
			hAvg3pC112vsCentPPUsingZDCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCA),wgtComb2part);
			hAvg3pC112vsCentPPUsingZDCCACombine->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - (fPsiZNCC+fPsiZNCA)),wgtComb2part);
			hAvg3pC112vsCentPPUsingZDCCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCCA),wgtComb2part);
					
			  }
			  //else if(trk1Chrg < 0 && trk2Chrg < 0){  ///this is obviously last option (and every microseconds counts!)
			  else{
			hAvg3pC112vsCentNNUsingZDCC->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCC),wgtComb2part);
			hAvg3pC112vsCentNNUsingZDCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCA),wgtComb2part);
			hAvg3pC112vsCentNNUsingZDCCACombine->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - (fPsiZNCC+fPsiZNCA)),wgtComb2part);
			hAvg3pC112vsCentNNUsingZDCCA->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiZNCCA),wgtComb2part);
			
			  }
			  
		  }
	      
	      
	    }//j-track trackCuts applied
	  }//j-track FB validated
	}///j-track loop ends
	

      }//----> i-track => All trackCuts applied.     
    }//-----> i-track => FB is validated.    
  }///-----> i-track loop Ends <--------
 


  

  if(bAnalysLambdaPairs) {

    ///------------------------------------------
    ///   Lambda PID:  @Chunzheng
    ///
    /// This part Below runs the loop over v0
    ///  and fills only lambda-xxx pairs:
    ///  
    ///------------------------------------------

    /// TPC Eta-sub Event Planes:
    // SumQnxTPCPos and SumQnyTPCPos are from +ve Eta side of TPC
    // SumQnxTPCNeg and SumQnyTPCNeg are from -Ve Eta side of TPC
    
    
 
    ////////-----> Starting V0 Loop -----------
    vector<Int_t>    vecLambdaCode;
    vector<Double_t> vecLambdaPhi;
    vector<Double_t> vecLambdaPt;
    vector<Int_t>    vecDaughterPosID; // Pos Daughter ID
    vector<Int_t>    vecDaughterNegID; // Neg Daughter ID

    fCurrentVtx[0] = -999;
    fCurrentVtx[1] = -999;
    fCurrentVtx[2] = -999;
  
    AliAODVertex* fVtx = fAOD->GetPrimaryVertex(); 
    fVtx->GetXYZ(fCurrentVtx);

    Int_t nV0s = fAOD->GetNumberOfV0s();
  
    for (Int_t iV0 = 0; iV0 < nV0s; iV0++) {
      AliAODv0 *v0 = fAOD->GetV0(iV0);
      if(!v0) continue;
      //Basic kinematic variable
      Double_t pt  = v0->Pt();
      Double_t eta = v0->PseudoRapV0();
      Double_t dcaToPV = v0->DcaV0ToPrimVertex();//DCA to Primary Vertex
      Double_t CPA = v0->CosPointingAngle(fCurrentVtx);//cosine pointing angle
      Double_t dl  = v0->DecayLengthV0(fCurrentVtx);
      fHistV0Pt->Fill(pt);
      fHistV0Eta->Fill(eta);
      fHistV0DcatoPrimVertex->Fill(dcaToPV);
      fHistV0CPA->Fill(CPA);
      fHistV0DecayLength->Fill(dl);

      //V0 cut
      if(!IsGoodV0(v0)) continue;
    
      //V0 daughters cut
      AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
      AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
    
      if(!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
      Float_t nDcaPV = v0->DcaNegToPrimVertex();
      Float_t pDcaPV = v0->DcaPosToPrimVertex();
      if( nDcaPV<fDaughtersDCAToPrimVtxMin || pDcaPV<fDaughtersDCAToPrimVtxMin) continue;

      Int_t code = GetLambdaCode(pTrack,nTrack);
      if(fabs(code) != 3122) continue;
      TVector2 Vt(v0->MomV0X(), v0->MomV0Y());
      Double_t phi = Vt.Phi() > 0 ? Vt.Phi() : Vt.Phi() + TMath::TwoPi();
      Int_t id_posDaughter = v0->GetPosID();
      Int_t id_negDaughter = v0->GetNegID();

      if (code == 3122) {
        Double_t massLambda  = v0->MassLambda();
        fHistLambdaPt[0]              -> Fill(pt);
        fHistLambdaEta[0]             -> Fill(eta);
        fHistLambdaPhi[0]             -> Fill(phi);
        fHistLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
        fHistLambdaCPA[0]             -> Fill(CPA);
        fHistLambdaDecayLength[0]     -> Fill(dl);
        fHistLambdaMass[0]            -> Fill(massLambda);
        fProfileLambdaMassVsPt[0]     -> Fill(pt, massLambda);

        if (TMath::Abs(massLambda - fMassMean) < fLambdaMassCut) {
          //if a particle has been used as daughter particle before(It happends), we have to refuse a new one.
          if (find(vecDaughterPosID.begin(), vecDaughterPosID.end(), id_posDaughter) != vecDaughterPosID.end()) continue;
          if (find(vecDaughterNegID.begin(), vecDaughterNegID.end(), id_negDaughter) != vecDaughterNegID.end()) continue;
          fHistLambdaPt[1]              -> Fill(pt);
          fHistLambdaEta[1]             -> Fill(eta);
          fHistLambdaPhi[1]             -> Fill(phi);
          fHistLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
          fHistLambdaCPA[1]             -> Fill(CPA);
          fHistLambdaDecayLength[1]     -> Fill(dl);
          fHistLambdaMass[1]            -> Fill(massLambda);
          fProfileLambdaMassVsPt[1]     -> Fill(pt, massLambda);

          vecLambdaCode.push_back(code);
          vecLambdaPhi.push_back(phi);
          vecLambdaPt.push_back(pt);
          vecDaughterPosID.push_back(id_posDaughter);
          vecDaughterNegID.push_back(id_negDaughter);
        }
      } 

      if (code == -3122) {
        Double_t massAntiLambda  = v0->MassAntiLambda();
        fHistAntiLambdaPt[0]              -> Fill(pt);
        fHistAntiLambdaEta[0]             -> Fill(eta);
        fHistAntiLambdaPhi[0]             -> Fill(phi);
        fHistAntiLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
        fHistAntiLambdaCPA[0]             -> Fill(CPA);
        fHistAntiLambdaDecayLength[0]     -> Fill(dl);
        fHistAntiLambdaMass[0]            -> Fill(massAntiLambda);
        fProfileAntiLambdaMassVsPt[0]     -> Fill(pt, massAntiLambda);

        if (TMath::Abs(massAntiLambda - fMassMean) < fLambdaMassCut) {
          if (find(vecDaughterPosID.begin(), vecDaughterPosID.end(), id_posDaughter) != vecDaughterPosID.end()) continue;
          if (find(vecDaughterNegID.begin(), vecDaughterNegID.end(), id_negDaughter) != vecDaughterNegID.end()) continue;
          fHistAntiLambdaPt[1]              -> Fill(pt);
          fHistAntiLambdaEta[1]             -> Fill(eta);
          fHistAntiLambdaPhi[1]             -> Fill(phi);
          fHistAntiLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
          fHistAntiLambdaCPA[1]             -> Fill(CPA);
          fHistAntiLambdaDecayLength[1]     -> Fill(dl);
          fHistAntiLambdaMass[1]            -> Fill(massAntiLambda);
          fProfileAntiLambdaMassVsPt[1]     -> Fill(pt, massAntiLambda);
          
          vecLambdaCode.push_back(code);
          vecLambdaPhi.push_back(phi);
          vecLambdaPt.push_back(pt);
          vecDaughterPosID.push_back(id_posDaughter);
          vecDaughterNegID.push_back(id_negDaughter);
        }
      } 
    }

    //// Now fill Lambda-X correlations:
    //Lambda - X
    for (vector<Double_t>::size_type jLambda = 0; jLambda < vecLambdaPhi.size(); jLambda++) {
      Int_t    code_lambda    = vecLambdaCode[jLambda];
      Double_t phi_lambda     = vecLambdaPhi[jLambda];
      Double_t pt_lambda      = vecLambdaPt[jLambda];
      Int_t    id_posDaughter = vecDaughterPosID[jLambda];
      Int_t    id_negDaughter = vecDaughterNegID[jLambda];

      //Remove AutoCorrelation:
      Double_t fPosTPCQxTemp = 0., fPosTPCQyTemp = 0.; ///Get the Total sum of Pos TPC Qx,Qy  locally, then Remove AutoCorr if needed.
      Double_t fNegTPCQxTemp = 0., fNegTPCQyTemp = 0.; ///Get the Total sum of Neg TPC Qx,Qy  locally, then Remove AutoCorr if needed.
      Double_t fPosTPCMult = 0.,   fNegTPCMult = 0.;
      Double_t qxPosTPC = 0., qyPosTPC = 0.;
      Double_t qxNegTPC = 0., qyNegTPC = 0.;

      // use EP from opposite eta than the charged track! One way to remove AutoCorrelation.
      fNegTPCQxTemp = fSumQnxNeg[0];
      fNegTPCQyTemp = fSumQnyNeg[0];
      fNegTPCMult   = fMultNeg;

      fPosTPCQxTemp = fSumQnxPos[0];
      fPosTPCQyTemp = fSumQnyPos[0];
      fPosTPCMult   = fMultPos;
      
      //for NegEP
      vector<Int_t>::iterator iterPosDaughterNegTPC = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_posDaughter);
      vector<Int_t>::iterator iterNegDaughterNegTPC = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_negDaughter);
      
      //for PosEP
      vector<Int_t>::iterator iterPosDaughterPosTPC = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_posDaughter);
      vector<Int_t>::iterator iterNegDaughterPosTPC = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_negDaughter);
	
      if (iterPosDaughterNegTPC != vecNegEPTrkID.end()) {
        Int_t iDaughter = distance(vecNegEPTrkID.begin(), iterPosDaughterNegTPC);
        qxNegTPC +=  vecNegEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecNegEPTrkPhi[iDaughter]);
        qyNegTPC +=  vecNegEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecNegEPTrkPhi[iDaughter]);
        fNegTPCMult -= vecNegEPTrkNUAWgt[iDaughter];
      }

      if (iterNegDaughterNegTPC != vecNegEPTrkID.end()) {
        Int_t iDaughter = distance(vecNegEPTrkID.begin(), iterNegDaughterNegTPC);
        qxNegTPC += vecNegEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecNegEPTrkPhi[iDaughter]);
        qyNegTPC += vecNegEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecNegEPTrkPhi[iDaughter]);
        fNegTPCMult -= vecNegEPTrkNUAWgt[iDaughter];
      }

      if (iterPosDaughterPosTPC != vecPosEPTrkID.end()) {
        Int_t iDaughter = distance(vecPosEPTrkID.begin(), iterPosDaughterPosTPC);
        qxPosTPC += vecPosEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecPosEPTrkPhi[iDaughter]);
        qyPosTPC += vecPosEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecPosEPTrkPhi[iDaughter]);
        fPosTPCMult -= vecPosEPTrkNUAWgt[iDaughter];
      }
    
      if (iterNegDaughterPosTPC != vecPosEPTrkID.end()) {
        Int_t iDaughter = distance(vecPosEPTrkID.begin(), iterNegDaughterPosTPC);
        qxPosTPC +=  vecPosEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecPosEPTrkPhi[iDaughter]);
        qyPosTPC +=  vecPosEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecPosEPTrkPhi[iDaughter]);
        fPosTPCMult -= vecPosEPTrkNUAWgt[iDaughter];
      }

      fPosTPCQxTemp -= qxPosTPC;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used.. 
      fPosTPCQyTemp -= qyPosTPC;

      fNegTPCQxTemp -= qxNegTPC;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used.. 
      fNegTPCQyTemp -= qyNegTPC;
    
      if(fPosTPCMult > 0. && fNegTPCMult > 0.) {
        fPosTPCQxTemp = fPosTPCQxTemp/fPosTPCMult;
        fPosTPCQyTemp = fPosTPCQyTemp/fPosTPCMult;
        fNegTPCQxTemp = fNegTPCQxTemp/fNegTPCMult;
        fNegTPCQyTemp = fNegTPCQyTemp/fNegTPCMult;
      } else continue;

      Double_t fPsiNPosTPCNoAuto = (1./2)*TMath::ATan2(fPosTPCQyTemp,fPosTPCQxTemp);   //AutoCorrelation Removed Pos EP.
      if(fPsiNPosTPCNoAuto < 0) fPsiNPosTPCNoAuto += TMath::TwoPi()/2.0;
      Double_t fPsiNNegTPCNoAuto = (1./2)*TMath::ATan2(fNegTPCQyTemp,fNegTPCQxTemp);   //AutoCorrelation Removed Neg EP.
      if(fPsiNNegTPCNoAuto < 0) fPsiNNegTPCNoAuto += TMath::TwoPi()/2.0;


      //Check the PID Flow
      if(bCheckPIDFlow) {

        Double_t v2TmpTPCPos = TMath::Cos(2*(phi_lambda - fPsiNPosTPCNoAuto));
        Double_t v2TmpTPCNeg = TMath::Cos(2*(phi_lambda - fPsiNNegTPCNoAuto));
        Double_t v2TmpV0C = TMath::Cos(2*(phi_lambda - fPsiNV0C));
        Double_t v2TmpV0A = TMath::Cos(2*(phi_lambda - fPsiNV0A));
        Double_t v2TmpZNC = TMath::Cos(2*(phi_lambda - fPsiZNCC));
        Double_t v2TmpZNA = TMath::Cos(2*(phi_lambda - fPsiZNCA));

        if(code_lambda == 3122) {
          fProfile2DRawFlowCentPtLambda[0]->Fill(centrality, pt_lambda, v2TmpTPCPos);
          fProfile2DRawFlowCentPtLambda[1]->Fill(centrality, pt_lambda, v2TmpTPCNeg);
          fProfile2DRawFlowCentPtLambda[2]->Fill(centrality, pt_lambda, v2TmpV0C);
          fProfile2DRawFlowCentPtLambda[3]->Fill(centrality, pt_lambda, v2TmpV0A);
          fProfile2DRawFlowCentPtLambda[4]->Fill(centrality, pt_lambda, v2TmpZNC);
          fProfile2DRawFlowCentPtLambda[5]->Fill(centrality, pt_lambda, v2TmpZNA);
        }

        if(code_lambda == -3122) {
          fProfile2DRawFlowCentPtAntiLambda[0]->Fill(centrality, pt_lambda, v2TmpTPCPos);
          fProfile2DRawFlowCentPtAntiLambda[1]->Fill(centrality, pt_lambda, v2TmpTPCNeg);
          fProfile2DRawFlowCentPtAntiLambda[2]->Fill(centrality, pt_lambda, v2TmpV0C);
          fProfile2DRawFlowCentPtAntiLambda[3]->Fill(centrality, pt_lambda, v2TmpV0A);
          fProfile2DRawFlowCentPtAntiLambda[4]->Fill(centrality, pt_lambda, v2TmpZNC);
          fProfile2DRawFlowCentPtAntiLambda[5]->Fill(centrality, pt_lambda, v2TmpZNA);
        }
      }

      for (vector<Double_t>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
        Int_t    id   = vecID[iTrk];
        Int_t    code = vecPDGCode[iTrk];
        Double_t pt   = vecPt[iTrk];
        Double_t eta  = vecEta[iTrk];
        Double_t phi  = vecPhi[iTrk];
        if(id == id_posDaughter || id == id_negDaughter) continue;

        Double_t fPsiNTPCNoAuto = -999;
        if (eta > 0.) fPsiNTPCNoAuto = fPsiNNegTPCNoAuto;
        else fPsiNTPCNoAuto = fPsiNPosTPCNoAuto;
  
        Double_t delta     = TMath::Cos(phi_lambda - phi);
        Double_t gammaTPC  = TMath::Cos(phi_lambda + phi - 2 *fPsiNTPCNoAuto);
        Double_t gammaV0C  = TMath::Cos(phi_lambda + phi - 2 *fPsiNV0C);
        Double_t gammaV0A  = TMath::Cos(phi_lambda + phi - 2 *fPsiNV0A);
        Double_t gammaZNC  = TMath::Cos(phi_lambda + phi - 2 *fPsiZNCC);
        Double_t gammaZNA  = TMath::Cos(phi_lambda + phi - 2 *fPsiZNCA);

        if (code > 0 && code_lambda ==  3122) {
          fProfileDelta_Lambda_hPos        -> Fill(centrality, delta);
          fProfileGammaTPC_Lambda_hPos     -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_Lambda_hPos     -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_Lambda_hPos     -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_Lambda_hPos     -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_Lambda_hPos     -> Fill(centrality, gammaZNA);
        }
        if (code < 0 && code_lambda ==  3122) {
          fProfileDelta_Lambda_hNeg        -> Fill(centrality, delta);
          fProfileGammaTPC_Lambda_hNeg     -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_Lambda_hNeg     -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_Lambda_hNeg     -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_Lambda_hNeg     -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_Lambda_hNeg     -> Fill(centrality, gammaZNA);
        } 
        if (code > 0 && code_lambda == -3122) {
          fProfileDelta_AntiLambda_hPos    -> Fill(centrality, delta);
          fProfileGammaTPC_AntiLambda_hPos -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_AntiLambda_hPos -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_AntiLambda_hPos -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_AntiLambda_hPos -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_AntiLambda_hPos -> Fill(centrality, gammaZNA);
        }
        if (code < 0 && code_lambda == -3122) {
          fProfileDelta_AntiLambda_hNeg    -> Fill(centrality, delta);
          fProfileGammaTPC_AntiLambda_hNeg -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_AntiLambda_hNeg -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_AntiLambda_hNeg -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_AntiLambda_hNeg -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_AntiLambda_hNeg -> Fill(centrality, gammaZNA);
        }
        if (code ==  2212 && code_lambda ==  3122) {
          fProfileDelta_Lambda_Proton      -> Fill(centrality, delta);
          fProfileGammaTPC_Lambda_Proton   -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_Lambda_Proton   -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_Lambda_Proton   -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_Lambda_Proton   -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_Lambda_Proton   -> Fill(centrality, gammaZNA);
        }
        if (code == -2212 && code_lambda ==  3122) {
          fProfileDelta_Lambda_AntiProton    -> Fill(centrality, delta);
          fProfileGammaTPC_Lambda_AntiProton -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_Lambda_AntiProton -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_Lambda_AntiProton -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_Lambda_AntiProton -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_Lambda_AntiProton -> Fill(centrality, gammaZNA);
        }
        if (code ==  2212 && code_lambda == -3122) {
          fProfileDelta_AntiLambda_Proton    -> Fill(centrality, delta);
          fProfileGammaTPC_AntiLambda_Proton -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_AntiLambda_Proton -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_AntiLambda_Proton -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_AntiLambda_Proton -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_AntiLambda_Proton -> Fill(centrality, gammaZNA);
        }
        if (code == -2212 && code_lambda == -3122) {
          fProfileDelta_AntiLambda_AntiProton    -> Fill(centrality, delta);
          fProfileGammaTPC_AntiLambda_AntiProton -> Fill(centrality, gammaTPC);
          fProfileGammaV0C_AntiLambda_AntiProton -> Fill(centrality, gammaV0C);
          fProfileGammaV0A_AntiLambda_AntiProton -> Fill(centrality, gammaV0A);
          fProfileGammaZNC_AntiLambda_AntiProton -> Fill(centrality, gammaZNC);
          fProfileGammaZNA_AntiLambda_AntiProton -> Fill(centrality, gammaZNA);
        }
      }// (Anti)Lambda-X pair done 
    }/// loop over charge particle array

    if(bCheckPIDFlow) {
      for (vector<Double_t>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
        Int_t    code = vecPDGCode[iTrk];
        Double_t pt   = vecPt[iTrk];
        Double_t eta  = vecEta[iTrk];
        Double_t phi  = vecPhi[iTrk];

        Double_t v2TmpTPCPos = TMath::Cos(2*(phi - fPsiNTPCPos));
        Double_t v2TmpTPCNeg = TMath::Cos(2*(phi - fPsiNTPCNeg));
        Double_t v2TmpV0C    = TMath::Cos(2*(phi - fPsiNV0C));
        Double_t v2TmpV0A    = TMath::Cos(2*(phi - fPsiNV0A));
        Double_t v2TmpZNC    = TMath::Cos(2*(phi - fPsiZNCC));
        Double_t v2TmpZNA    = TMath::Cos(2*(phi - fPsiZNCA));
        
        if (code > 0) {
          if (eta > 0) fProfile2DRawFlowCentPthPos[0]->Fill(centrality, pt, v2TmpTPCNeg);
          else         fProfile2DRawFlowCentPthPos[0]->Fill(centrality, pt, v2TmpTPCPos);
          fProfile2DRawFlowCentPthPos[1]->Fill(centrality, pt, v2TmpV0C);
          fProfile2DRawFlowCentPthPos[2]->Fill(centrality, pt, v2TmpV0A);
          fProfile2DRawFlowCentPthPos[3]->Fill(centrality, pt, v2TmpZNC);
          fProfile2DRawFlowCentPthPos[4]->Fill(centrality, pt, v2TmpZNA);
        }
        if (code < 0) {
          if (eta > 0) fProfile2DRawFlowCentPthNeg[0]->Fill(centrality, pt, v2TmpTPCNeg);
          else         fProfile2DRawFlowCentPthNeg[0]->Fill(centrality, pt, v2TmpTPCPos);
          fProfile2DRawFlowCentPthNeg[1]->Fill(centrality, pt, v2TmpV0C);
          fProfile2DRawFlowCentPthNeg[2]->Fill(centrality, pt, v2TmpV0A);
          fProfile2DRawFlowCentPthNeg[3]->Fill(centrality, pt, v2TmpZNC);
          fProfile2DRawFlowCentPthNeg[4]->Fill(centrality, pt, v2TmpZNA);
        }
        if (code == 2212) {
          if (eta > 0) fProfile2DRawFlowCentPtProton[0]->Fill(centrality, pt, v2TmpTPCNeg);
          else         fProfile2DRawFlowCentPtProton[0]->Fill(centrality, pt, v2TmpTPCPos);
          fProfile2DRawFlowCentPtProton[1]->Fill(centrality, pt, v2TmpV0C);
          fProfile2DRawFlowCentPtProton[2]->Fill(centrality, pt, v2TmpV0A);
          fProfile2DRawFlowCentPtProton[3]->Fill(centrality, pt, v2TmpZNC);
          fProfile2DRawFlowCentPtProton[4]->Fill(centrality, pt, v2TmpZNA);
        }
        if (code == -2212) {
          if (eta > 0) fProfile2DRawFlowCentPtAntiProton[0]->Fill(centrality, pt, v2TmpTPCNeg);
          else         fProfile2DRawFlowCentPtAntiProton[0]->Fill(centrality, pt, v2TmpTPCPos);
          fProfile2DRawFlowCentPtAntiProton[1]->Fill(centrality, pt, v2TmpV0C);
          fProfile2DRawFlowCentPtAntiProton[2]->Fill(centrality, pt, v2TmpV0A);
          fProfile2DRawFlowCentPtAntiProton[3]->Fill(centrality, pt, v2TmpZNC);
          fProfile2DRawFlowCentPtAntiProton[4]->Fill(centrality, pt, v2TmpZNA);
        }
      }////--------- PID Flow hist are Filled ----------
    }
  }////--------- Lambda Histograms are Filled ----------
 
  CalculateVnfromQC(centrality); //@ Shi calculate v2{4} for ZDC resolution
  ResetEventByEventQuantities();

  fDebugwEventCount->Fill(9.1); ///Left for Analysis
  
  fHistVertexZcm->Fill(pVtxZ);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

}//---------------- UserExec ----------------------

void AliAnalysisTaskGammaDeltaPID::CalculateVnfromQC(Double_t cent){

	Double_t QRe = hPOIQRe->GetBinContent(1);
	Double_t QIm = hPOIQIm->GetBinContent(1);
	Double_t Q2Re2 = hPOIQRe->GetBinContent(2);
	Double_t Q2Im2 = hPOIQIm->GetBinContent(2);
	Double_t QRe3 = hPOIQRe->GetBinContent(3);
	Double_t QIm3 = hPOIQIm->GetBinContent(3);

	Double_t QM0 = hPOIMul->GetBinContent(1);
	Double_t QM  = hPOIMul->GetBinContent(2);
	Double_t QM2 = hPOIMul->GetBinContent(3);
	Double_t QM3 = hPOIMul->GetBinContent(4);
	Double_t QM4 = hPOIMul->GetBinContent(5);
	
	Double_t IQM2 = QM*QM-QM2;
    Double_t WQM2 = IQM2;
    Double_t IQC2 = -99.;
    
    Bool_t Q2f = kFALSE, Q4f = kFALSE;
    
	if(QM0>1) {
	  IQC2 = (QRe*QRe+QIm*QIm-QM2)/IQM2;
	  hFlowQC2Pro->Fill(cent,IQC2,WQM2);
	  Q2f = kTRUE;
	}
      
	Double_t IQM4 = QM*QM*QM*QM - 6.*QM2*QM*QM + 8.*QM3*QM + 3.*QM2*QM2 - 6.*QM4;
    Double_t WQM4 = IQM4;
    Double_t IQC4 = -99.;
    if(QM0>3) {
      IQC4 = ((QRe*QRe+QIm*QIm)*(QRe*QRe+QIm*QIm)
                  - 2.*(QRe*QRe*Q2Re2+2.*QRe*QIm*Q2Im2-QIm*QIm*Q2Re2)
                  + 8.*(QRe3*QRe+QIm3*QIm)
                  + (Q2Re2*Q2Re2+Q2Im2*Q2Im2)
                  - 4.*QM2*(QRe*QRe+QIm*QIm)
                  - 6.*QM4 + 2.*QM2*QM2) / IQM4;

      hFlowQC4Pro->Fill(cent,IQC4,WQM4);
      Q4f = kTRUE;
    }
    
    // product of correlations or covariances
	if(Q2f && Q4f) {
	  hFlowwCov24Pro->Fill(cent,IQC2*IQC4,WQM2*WQM4);
	}

}


void AliAnalysisTaskGammaDeltaPID::ResetEventByEventQuantities(){
	if(hPOIQRe) hPOIQRe->Reset();    // 
	if(hPOIQIm) hPOIQIm->Reset();    // 
	if(hPOIMul) hPOIMul->Reset();    // 
	
}

void AliAnalysisTaskGammaDeltaPID::SetupAnalysisHistograms(){
  
  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables
  Char_t  name[100];
  Char_t title[100];

  //// CME PID ANALYSIS Histograms/Profiles:

  //// different types of v2
  hV2TPCvsCent = new TProfile("hV2TPCvsCent"," v2_{TPC} vs Cent; Cent; v2_{TPC}",18,0,90);
  hV2TPCvsCent->Sumw2();
  fListHist->Add(hV2TPCvsCent);
  hV2V0CvsCent = new TProfile("hV2V0CvsCent"," v2_{V0C} vs Cent; Cent; v2_{V0C}",18,0,90);
  hV2V0CvsCent->Sumw2();
  fListHist->Add(hV2V0CvsCent);
  hV2V0AvsCent = new TProfile("hV2V0AvsCent"," v2_{V0A} vs Cent; Cent; v2_{V0A}",18,0,90);
  hV2V0AvsCent->Sumw2();
  fListHist->Add(hV2V0AvsCent);
  hV2V0CAvsCent = new TProfile("hV2V0CAvsCent"," v2_{V0CA} vs Cent; Cent; v2_{V0CA}",18,0,90);
  hV2V0CAvsCent->Sumw2();
  fListHist->Add(hV2V0CAvsCent);
  hV2ZDCCvsCent = new TProfile("hV2ZDCCvsCent"," v2_{ZDCC} vs Cent; Cent; v2_{ZDCC}",18,0,90);
  hV2ZDCCvsCent->Sumw2();
  fListHist->Add(hV2ZDCCvsCent);
  hV2ZDCAvsCent = new TProfile("hV2ZDCAvsCent"," v2_{ZDCA} vs Cent; Cent; v2_{ZDCA}",18,0,90);
  hV2ZDCAvsCent->Sumw2();
  fListHist->Add(hV2ZDCAvsCent);
  hV2ZDCCACombinevsCent = new TProfile("hV2ZDCCACombinevsCent"," v2_{ZDCC+ZDCA} vs Cent; Cent; v2_{ZDCC+ZDCA}",18,0,90);
  hV2ZDCCACombinevsCent->Sumw2();
  fListHist->Add(hV2ZDCCACombinevsCent);
  hV2ZDCCAvsCent = new TProfile("hV2ZDCCAvsCent"," v2_{ZDCCA} vs Cent; Cent; v2_{ZDCCA}",18,0,90);
  hV2ZDCCAvsCent->Sumw2();
  fListHist->Add(hV2ZDCCAvsCent);
  
  hPOIQRe = new TH1D("hPOIQRe"," hPOIQRe",3,0,3);
  hPOIQRe->Sumw2();
  fListTemp->Add(hPOIQRe);
  hPOIQIm = new TH1D("hPOIQIm"," hPOIQIm",3,0,3);
  hPOIQIm->Sumw2();
  fListTemp->Add(hPOIQIm);
  hPOIMul = new TH1D("hPOIMul"," hPOIMul",5,0,5);
  hPOIMul->Sumw2();
  fListTemp->Add(hPOIMul);
  
  hFlowQC2Pro = new TProfile("hFlowQC2Pro"," QC2 vs Cent; Cent; QC2",18,0,90);
  hFlowQC2Pro->Sumw2();
  fListHist->Add(hFlowQC2Pro);
  hFlowQC4Pro = new TProfile("hFlowQC4Pro"," QC4 vs Cent; Cent; QC4",18,0,90);
  hFlowQC4Pro->Sumw2();
  fListHist->Add(hFlowQC4Pro);
  hFlowwCov24Pro = new TProfile("hFlowwCov24Pro","Cov24Pro vs Cent; Cent; v2_{ZDCCA}",18,0,90);
  hFlowwCov24Pro->Sumw2();
  fListHist->Add(hFlowwCov24Pro);
  
  //// 3p correlators PID:
  hAvg3pC112vsCentPP = new TProfile("hAvg3pC112vsCentPP"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC112vsCentPP->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPP);
  hAvg3pC112vsCentNN = new TProfile("hAvg3pC112vsCentNN"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);  
  hAvg3pC112vsCentNN->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNN);
  hAvg3pC112vsCentOS = new TProfile("hAvg3pC112vsCentOS"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC112vsCentOS->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOS);

  hAvg3pC123vsCentPP = new TProfile("hAvg3pC123vsCentPP"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvg3pC123vsCentPP->Sumw2();
  fListHist->Add(hAvg3pC123vsCentPP);
  hAvg3pC123vsCentNN = new TProfile("hAvg3pC123vsCentNN"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);  
  hAvg3pC123vsCentNN->Sumw2();
  fListHist->Add(hAvg3pC123vsCentNN);
  hAvg3pC123vsCentOS = new TProfile("hAvg3pC123vsCentOS"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvg3pC123vsCentOS->Sumw2();
  fListHist->Add(hAvg3pC123vsCentOS);  

  /// 2p Correlator PID:  
  hAvgDelta1vsCentPP = new TProfile("hAvgDelta1vsCentPP"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta1vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta1vsCentPP);
  hAvgDelta1vsCentNN = new TProfile("hAvgDelta1vsCentNN"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);    
  hAvgDelta1vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta1vsCentNN);  
  hAvgDelta1vsCentOS = new TProfile("hAvgDelta1vsCentOS"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta1vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta1vsCentOS);

  hAvgDelta2vsCentPP = new TProfile("hAvgDelta2vsCentPP"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta2vsCentPP);
  hAvgDelta2vsCentNN = new TProfile("hAvgDelta2vsCentNN"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);    
  hAvgDelta2vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta2vsCentNN);  
  hAvgDelta2vsCentOS = new TProfile("hAvgDelta2vsCentOS"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta2vsCentOS);

  hAvgDelta3vsCentPP = new TProfile("hAvgDelta3vsCentPP"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta3vsCentPP);
  hAvgDelta3vsCentNN = new TProfile("hAvgDelta3vsCentNN"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);    
  hAvgDelta3vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta3vsCentNN);  
  hAvgDelta3vsCentOS = new TProfile("hAvgDelta3vsCentOS"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta3vsCentOS);

  hAvgDelta4vsCentPP = new TProfile("hAvgDelta4vsCentPP"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta4vsCentPP);
  hAvgDelta4vsCentNN = new TProfile("hAvgDelta4vsCentNN"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);    
  hAvgDelta4vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta4vsCentNN);  
  hAvgDelta4vsCentOS = new TProfile("hAvgDelta4vsCentOS"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta4vsCentOS);
  
  /// 3p correlators PID using SP from ZDC
  hAvg3pC112vsCentPPUsingZDCC = new TProfile("hAvg3pC112vsCentPPUsingZDCC"," <3p> vs Cent; Cent; <#gamma112> Using ZDCC",18,0,90);
  hAvg3pC112vsCentPPUsingZDCC->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPPUsingZDCC);
  hAvg3pC112vsCentNNUsingZDCC = new TProfile("hAvg3pC112vsCentNNUsingZDCC"," <3p> vs Cent; Cent; <#gamma112> Using ZDCC",18,0,90);  
  hAvg3pC112vsCentNNUsingZDCC->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNNUsingZDCC);
  hAvg3pC112vsCentOSUsingZDCC = new TProfile("hAvg3pC112vsCentOSUsingZDCC"," <3p> vs Cent; Cent; <#gamma112> Using ZDCC",18,0,90);
  hAvg3pC112vsCentOSUsingZDCC->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOSUsingZDCC);
  
  hAvg3pC112vsCentPPUsingZDCA = new TProfile("hAvg3pC112vsCentPPUsingZDCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCA",18,0,90);
  hAvg3pC112vsCentPPUsingZDCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPPUsingZDCA);
  hAvg3pC112vsCentNNUsingZDCA = new TProfile("hAvg3pC112vsCentNNUsingZDCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCA",18,0,90);  
  hAvg3pC112vsCentNNUsingZDCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNNUsingZDCA);
  hAvg3pC112vsCentOSUsingZDCA = new TProfile("hAvg3pC112vsCentOSUsingZDCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCA",18,0,90);
  hAvg3pC112vsCentOSUsingZDCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOSUsingZDCA);
  
  hAvg3pC112vsCentPPUsingZDCCACombine = new TProfile("hAvg3pC112vsCentPPUsingZDCCACombine"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA Combine",18,0,90);
  hAvg3pC112vsCentPPUsingZDCCACombine->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPPUsingZDCCACombine);
  hAvg3pC112vsCentNNUsingZDCCACombine = new TProfile("hAvg3pC112vsCentNNUsingZDCCACombine"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA Combine",18,0,90);  
  hAvg3pC112vsCentNNUsingZDCCACombine->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNNUsingZDCCACombine);
  hAvg3pC112vsCentOSUsingZDCCACombine = new TProfile("hAvg3pC112vsCentOSUsingZDCCACombine"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA Combine",18,0,90);
  hAvg3pC112vsCentOSUsingZDCCACombine->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOSUsingZDCCACombine);
  
  hAvg3pC112vsCentPPUsingZDCCA = new TProfile("hAvg3pC112vsCentPPUsingZDCCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA",18,0,90);
  hAvg3pC112vsCentPPUsingZDCCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPPUsingZDCCA);
  hAvg3pC112vsCentNNUsingZDCCA = new TProfile("hAvg3pC112vsCentNNUsingZDCCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA",18,0,90);  
  hAvg3pC112vsCentNNUsingZDCCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNNUsingZDCCA);
  hAvg3pC112vsCentOSUsingZDCCA = new TProfile("hAvg3pC112vsCentOSUsingZDCCA"," <3p> vs Cent; Cent; <#gamma112> Using ZDCCA",18,0,90);
  hAvg3pC112vsCentOSUsingZDCCA->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOSUsingZDCCA);

  ///TPC Event Plane Correlations for Resolution Estimation
  hTPCPsiNCorrelation = new TProfile("hTPCPsiNCorrelation",Form("Psi_{%d} Res. vs Cent; Cent; Resolution",gHarmonic),90,0,90);
  fListHist->Add(hTPCPsiNCorrelation);
  hTPCPsi3Correlation = new TProfile("hTPCPsi3Correlation",Form("Psi_{%d} Res. vs Cent; Cent; Resolution",3),90,0,90);
  fListHist->Add(hTPCPsi3Correlation);
  hTPCPsi4Correlation = new TProfile("hTPCPsi4Correlation",Form("Psi_{%d} Res. vs Cent; Cent; Resolution",4),90,0,90);
  fListHist->Add(hTPCPsi4Correlation);    

  ///V0X-V0X/TPC  Event Plane Correlations for Resolution:
  hV0CV0APsiNCorrelation = new TProfile("hV0CV0APsiNCorrelation",Form("V0C-V0A Psi%d; Cent; Resolution",gHarmonic),90,0,90);
  fListHist->Add(hV0CV0APsiNCorrelation);
  hV0CTPCPsiNCorrelation = new TProfile("hV0CTPCPsiNCorrelation",Form("V0C-TPC Psi%d; Cent; Resolution",gHarmonic),90,0,90);
  fListHist->Add(hV0CTPCPsiNCorrelation);
  hV0ATPCPsiNCorrelation = new TProfile("hV0ATPCPsiNCorrelation",Form("V0A-TPC Psi%d; Cent; Resolution",gHarmonic),90,0,90);
  fListHist->Add(hV0ATPCPsiNCorrelation);

  hV0CV0APsi3Correlation = new TProfile("hV0CV0APsi3Correlation",Form("V0C-V0A Psi%d; Cent; Resolution",3),90,0,90);
  fListHist->Add(hV0CV0APsi3Correlation);
  hV0CTPCPsi3Correlation = new TProfile("hV0CTPCPsi3Correlation",Form("V0C-TPC Psi%d; Cent; Resolution",3),90,0,90);
  fListHist->Add(hV0CTPCPsi3Correlation);
  hV0ATPCPsi3Correlation = new TProfile("hV0ATPCPsi3Correlation",Form("V0A-TPC Psi%d; Cent; Resolution",3),90,0,90);
  fListHist->Add(hV0ATPCPsi3Correlation);
   
}







void AliAnalysisTaskGammaDeltaPID::SetupQAHistograms(){

  Double_t fVzBinsZDC[21] = {-10,-8,-7,-6,-5,-4,-3,-2,-1,-0.5, 0, 0.5, 1.,2.,3.,4.,5.,6.,7.,8.,10.};
  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables

  Int_t gMaxGlobalmult  = 4000;
  Int_t gMaxTPCcorrmult = 5000;
  Int_t gMaxESDtracks   = 20000;
  
  /// Bunch of QA histograms:  
  fHistVertexZcm = new TH1F("fHistVertexZcm"," V_{z}; V_{z} cm; events ",100,-10,10);
  fListHist->Add(fHistVertexZcm);
  //fHistVxvsVzMinBias = new TProfile("fHistVxvsVzMinBias"," <Vx> vs Vz; V_{z}cm; <V_{x}> cm", 20,fVzBinsZDC);
  //fListHist->Add(fHistVxvsVzMinBias);
  //fHistVyvsVzMinBias = new TProfile("fHistVyvsVzMinBias"," <Vy> vs Vz; V_{z}cm; <V_{y}> cm", 20,fVzBinsZDC);
  //fListHist->Add(fHistVyvsVzMinBias);

  fCentDistBeforCut = new TH1F("fCentDistBeforCut","Cent w/o any Cuts; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistBeforCut);
  fCentDistAfterCut = new TH1F("fCentDistAfterCut","Cent with all Cut; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistAfterCut);
  fDebugwEventCount = new TH1F("fDebugwEventCount","Event no. count; ; no.Events ",10,0,10);
  fListHist->Add(fDebugwEventCount);

  fDebugwEventCount->GetXaxis()->SetBinLabel(1,"AOD exist");
  fDebugwEventCount->GetXaxis()->SetBinLabel(2,"Vz cuts");
  fDebugwEventCount->GetXaxis()->SetBinLabel(3,"Centrality cut");
  fDebugwEventCount->GetXaxis()->SetBinLabel(4,"nTrkAOD>4");
  fDebugwEventCount->GetXaxis()->SetBinLabel(5,"PileUp cut");
  fDebugwEventCount->GetXaxis()->SetBinLabel(6,"TPCmult>0.1");
  fDebugwEventCount->GetXaxis()->SetBinLabel(7,"V0Mult>0.0001");
  fDebugwEventCount->GetXaxis()->SetBinLabel(10,"Analyzed");
  
  /// 2D QA:
  
  ///TPC Event Planes:
  fHistTPCPsiNPosPlane = new TH2F("fHistTPCPsiNPosPlane",Form("#Psi_{n}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsiNPosPlane);
  fHistTPCPsiNNegPlane = new TH2F("fHistTPCPsiNNegPlane",Form("#Psi_{n}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsiNNegPlane);
  fHistTPCPsi3PosPlane = new TH2F("fHistTPCPsi3PosPlane",Form("#Psi_{3}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi3PosPlane);
  fHistTPCPsi3NegPlane = new TH2F("fHistTPCPsi3NegPlane",Form("#Psi_{3}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi3NegPlane);
  fHistTPCPsi4PosPlane = new TH2F("fHistTPCPsi4PosPlane",Form("#Psi_{3}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi4PosPlane);
  fHistTPCPsi4NegPlane = new TH2F("fHistTPCPsi4NegPlane",Form("#Psi_{3}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi4NegPlane);  
  
  /// V0 Event Planes:
  fHistV0CPsiNEventPlane = new TH2F("fHistV0CPsiNEventPlane",Form("#Psi_{n}(V0C); centrality; #Psi_{%d,V0C}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0CPsiNEventPlane);
  fHistV0APsiNEventPlane = new TH2F("fHistV0APsiNEventPlane",Form("#Psi_{n}(V0A); centrality; #Psi_{%d,V0A}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0APsiNEventPlane);  
  fHistV0CPsi3EventPlane = new TH2F("fHistV0CPsi3EventPlane",Form("#Psi_{3}(V0C); centrality; #Psi_{%d,V0C}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0CPsi3EventPlane);
  fHistV0APsi3EventPlane = new TH2F("fHistV0APsi3EventPlane",Form("#Psi_{3}(V0A); centrality; #Psi_{%d,V0A}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0APsi3EventPlane);    
  ///q-Vector for ESE:
  fHistTPCPosqVectorvsCent = new TH2F("fHistTPCPosqVectorvsCent","q TPC pos vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistTPCPosqVectorvsCent);
  fHistTPCNegqVectorvsCent = new TH2F("fHistTPCNegqVectorvsCent","q TPC neg vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistTPCNegqVectorvsCent);
  fHistV0CDetqVectorvsCent = new TH2F("fHistV0CDetqVectorvsCent","q V0Cdet  vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistV0CDetqVectorvsCent);
  fHistV0ADetqVectorvsCent = new TH2F("fHistV0ADetqVectorvsCent","q V0Adet  vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistV0ADetqVectorvsCent);  

  /// ZDC Spectator Planes:
  hZDCCPsiSpectatorPlane = new TH1F("hZDCCPsiSpectatorPlane","hZDCCPsiSpectatorPlane",100,-TMath::Pi(),TMath::Pi());
  hZDCCPsiSpectatorPlane->Sumw2();
  hZDCCPsiSpectatorPlane->SetXTitle("#Psi_{ZDCC}");
  fListHist->Add(hZDCCPsiSpectatorPlane);
  
  hZDCAPsiSpectatorPlane = new TH1F("hZDCAPsiSpectatorPlane","hZDCAPsiSpectatorPlane",100,-TMath::Pi(),TMath::Pi());
  hZDCAPsiSpectatorPlane->Sumw2();
  hZDCAPsiSpectatorPlane->SetXTitle("#Psi_{ZDCA}");
  fListHist->Add(hZDCAPsiSpectatorPlane);
  
  hZDCCAPsiSpectatorPlane = new TH1F("hZDCCAPsiSpectatorPlane","hZDCCAPsiSpectatorPlane",100,-TMath::Pi(),TMath::Pi());
  hZDCCAPsiSpectatorPlane->Sumw2();
  hZDCCAPsiSpectatorPlane->SetXTitle("#Psi_{ZDCCA}");
  fListHist->Add(hZDCCAPsiSpectatorPlane);
  
  hZDCCPsivsZDCAPsi = new TH2F("hZDCCPsivsZDCAPsi","hZDCCPsivsZDCAPsi",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  hZDCCPsivsZDCAPsi->Sumw2();
  hZDCCPsivsZDCAPsi->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsZDCAPsi->SetYTitle("#Psi_{ZDCA}");
  fListHist->Add(hZDCCPsivsZDCAPsi);
  
  hZDCCPsivsTPCPsi2 = new TH2F("hZDCCPsivsTPCPsi2","hZDCCPsivsTPCPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsTPCPsi2->Sumw2();
  hZDCCPsivsTPCPsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsTPCPsi2->SetYTitle("#Psi_{TPC}");
  fListHist->Add(hZDCCPsivsTPCPsi2);
  
  hZDCCPsivsTPCPosPsi2 = new TH2F("hZDCCPsivsTPCPosPsi2","hZDCCPsivsTPCPosPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsTPCPosPsi2->Sumw2();
  hZDCCPsivsTPCPosPsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsTPCPosPsi2->SetYTitle("#Psi_{TPC+}");
  fListHist->Add(hZDCCPsivsTPCPosPsi2);
  
  hZDCCPsivsTPCNegPsi2 = new TH2F("hZDCCPsivsTPCNegPsi2","hZDCCPsivsTPCNegPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsTPCNegPsi2->Sumw2();
  hZDCCPsivsTPCNegPsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsTPCNegPsi2->SetYTitle("#Psi_{TPC-}");
  fListHist->Add(hZDCCPsivsTPCNegPsi2);
  
  hZDCCPsivsV0CAPsi2 = new TH2F("hZDCCPsivsV0CAPsi2","hZDCCPsivsV0CAPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsV0CAPsi2->Sumw2();
  hZDCCPsivsV0CAPsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsV0CAPsi2->SetYTitle("#Psi_{V0CA}");
  fListHist->Add(hZDCCPsivsV0CAPsi2);
  
  hZDCCPsivsV0CPsi2 = new TH2F("hZDCCPsivsV0CPsi2","hZDCCPsivsV0CPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsV0CPsi2->Sumw2();
  hZDCCPsivsV0CPsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsV0CPsi2->SetYTitle("#Psi_{V0C}");
  fListHist->Add(hZDCCPsivsV0CPsi2);
  
  hZDCCPsivsV0APsi2 = new TH2F("hZDCCPsivsV0APsi2","hZDCCPsivsV0APsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCPsivsV0APsi2->Sumw2();
  hZDCCPsivsV0APsi2->SetXTitle("#Psi_{ZDCC}");
  hZDCCPsivsV0APsi2->SetYTitle("#Psi_{V0A}");
  fListHist->Add(hZDCCPsivsV0APsi2);
  
  
  hZDCAPsivsTPCPsi2 = new TH2F("hZDCAPsivsTPCPsi2","hZDCAPsivsTPCPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsTPCPsi2->Sumw2();
  hZDCAPsivsTPCPsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsTPCPsi2->SetYTitle("#Psi_{TPC}");
  fListHist->Add(hZDCAPsivsTPCPsi2);
  
  hZDCAPsivsTPCPosPsi2 = new TH2F("hZDCAPsivsTPCPosPsi2","hZDCAPsivsTPCPosPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsTPCPosPsi2->Sumw2();
  hZDCAPsivsTPCPosPsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsTPCPosPsi2->SetYTitle("#Psi_{TPC+}");
  fListHist->Add(hZDCAPsivsTPCPosPsi2);
  
  hZDCAPsivsTPCNegPsi2 = new TH2F("hZDCAPsivsTPCNegPsi2","hZDCAPsivsTPCNegPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsTPCNegPsi2->Sumw2();
  hZDCAPsivsTPCNegPsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsTPCNegPsi2->SetYTitle("#Psi_{TPC-}");
  fListHist->Add(hZDCAPsivsTPCNegPsi2);
  
  hZDCAPsivsV0CAPsi2 = new TH2F("hZDCAPsivsV0CAPsi2","hZDCAPsivsV0CAPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsV0CAPsi2->Sumw2();
  hZDCAPsivsV0CAPsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsV0CAPsi2->SetYTitle("#Psi_{V0CA}");
  fListHist->Add(hZDCAPsivsV0CAPsi2);
  
  hZDCAPsivsV0CPsi2 = new TH2F("hZDCAPsivsV0CPsi2","hZDCAPsivsV0CPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsV0CPsi2->Sumw2();
  hZDCAPsivsV0CPsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsV0CPsi2->SetYTitle("#Psi_{V0C}");
  fListHist->Add(hZDCAPsivsV0CPsi2);
  
  hZDCAPsivsV0APsi2 = new TH2F("hZDCAPsivsV0APsi2","hZDCAPsivsV0APsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCAPsivsV0APsi2->Sumw2();
  hZDCAPsivsV0APsi2->SetXTitle("#Psi_{ZDCA}");
  hZDCAPsivsV0APsi2->SetYTitle("#Psi_{V0A}");
  fListHist->Add(hZDCAPsivsV0APsi2);
  
  
  hZDCCAPsivsTPCPsi2 = new TH2F("hZDCCAPsivsTPCPsi2","hZDCCAPsivsTPCPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsTPCPsi2->Sumw2();
  hZDCCAPsivsTPCPsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsTPCPsi2->SetYTitle("#Psi_{TPC}");
  fListHist->Add(hZDCCAPsivsTPCPsi2);
  
  hZDCCAPsivsTPCPosPsi2 = new TH2F("hZDCCAPsivsTPCPosPsi2","hZDCCAPsivsTPCPosPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsTPCPosPsi2->Sumw2();
  hZDCCAPsivsTPCPosPsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsTPCPosPsi2->SetYTitle("#Psi_{TPC+}");
  fListHist->Add(hZDCCAPsivsTPCPosPsi2);
  
  hZDCCAPsivsTPCNegPsi2 = new TH2F("hZDCCAPsivsTPCNegPsi2","hZDCCAPsivsTPCNegPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsTPCNegPsi2->Sumw2();
  hZDCCAPsivsTPCNegPsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsTPCNegPsi2->SetYTitle("#Psi_{TPC-}");
  fListHist->Add(hZDCCAPsivsTPCNegPsi2);
  
  hZDCCAPsivsV0CAPsi2 = new TH2F("hZDCCAPsivsV0CAPsi2","hZDCCAPsivsV0CAPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsV0CAPsi2->Sumw2();
  hZDCCAPsivsV0CAPsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsV0CAPsi2->SetYTitle("#Psi_{V0CA}");
  fListHist->Add(hZDCCAPsivsV0CAPsi2);
  
  hZDCCAPsivsV0CPsi2 = new TH2F("hZDCCAPsivsV0CPsi2","hZDCCAPsivsV0CPsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsV0CPsi2->Sumw2();
  hZDCCAPsivsV0CPsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsV0CPsi2->SetYTitle("#Psi_{V0C}");
  fListHist->Add(hZDCCAPsivsV0CPsi2);
  
  hZDCCAPsivsV0APsi2 = new TH2F("hZDCCAPsivsV0APsi2","hZDCCAPsivsV0APsi2",100,-TMath::Pi(),TMath::Pi(),100,0,TMath::Pi());
  hZDCCAPsivsV0APsi2->Sumw2();
  hZDCCAPsivsV0APsi2->SetXTitle("#Psi_{ZDCCA}");
  hZDCCAPsivsV0APsi2->SetYTitle("#Psi_{V0A}");
  fListHist->Add(hZDCCAPsivsV0APsi2);

  
  /// Centrality Correlations and PileUp QA
  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MAfter);

  fHistCentCL0VsV0MBefore = new TH2F("fHistCentCL0VsV0MBefore","Before;Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHist->Add(fHistCentCL0VsV0MBefore);
  fHistCentCL0VsV0MAfter  = new TH2F("fHistCentCL0VsV0MAfter"," After; Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHist->Add(fHistCentCL0VsV0MAfter);

  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHist->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHist->Add(fHistTPCVsESDTrkAfter);



  /// V0 Mult vs Channel,Vz:
  hAvgV0ChannelsvsVz = new TProfile2D("hAvgV0ChannelsvsVz"," <V0 Ch(n)> ; Channel;  V_{z}(cm);",64,0,64,20,fVzBinsZDC);
  fListHist->Add(hAvgV0ChannelsvsVz);

  ///  <Q> for the V0C
  hAvgQNXvsCentV0C = new TProfile("hAvgQNXvsCentV0C"," V0C <Q2x> vs Cent; Cent; <Q_{2,x}>",90,0,90);
  fListHist->Add(hAvgQNXvsCentV0C);
  hAvgQNYvsCentV0C = new TProfile("hAvgQNYvsCentV0C"," V0C <Q2y> vs Cent; Cent; <Q_{2,y}>",90,0,90);
  fListHist->Add(hAvgQNYvsCentV0C);
  hAvgQ3XvsCentV0C = new TProfile("hAvgQ3XvsCentV0C"," V0C <Q3x> vs Cent; Cent; <Q_{3,x}>",90,0,90);
  fListHist->Add(hAvgQ3XvsCentV0C);
  hAvgQ3YvsCentV0C = new TProfile("hAvgQ3YvsCentV0C"," V0C <Q3y> vs Cent; Cent; <Q_{3,y}>",90,0,90);
  fListHist->Add(hAvgQ3YvsCentV0C);
  ///  <Q> for the V0A
  hAvgQNXvsCentV0A = new TProfile("hAvgQNXvsCentV0A"," V0A <Q2x> vs Cent; Cent; <Q_{2,x}>",90,0,90);
  fListHist->Add(hAvgQNXvsCentV0A);
  hAvgQNYvsCentV0A = new TProfile("hAvgQNYvsCentV0A"," V0A <Q2y> vs Cent; Cent; <Q_{2,y}>",90,0,90);
  fListHist->Add(hAvgQNYvsCentV0A);
  hAvgQ3XvsCentV0A = new TProfile("hAvgQ3XvsCentV0A"," V0A <Q3x> vs Cent; Cent; <Q_{3,x}>",90,0,90);
  fListHist->Add(hAvgQ3XvsCentV0A);
  hAvgQ3YvsCentV0A = new TProfile("hAvgQ3YvsCentV0A"," V0A <Q3y> vs Cent; Cent; <Q_{3,y}>",90,0,90);  
  fListHist->Add(hAvgQ3YvsCentV0A);

  /// Now <Q> for the TPC..  Note:gHarmonic is set as global variable!
  fAvgCos2PsivsCentEtaPos = new TProfile("fAvgCos2PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgCos2PsivsCentEtaPos);
  fAvgSin2PsivsCentEtaPos = new TProfile("fAvgSin2PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgSin2PsivsCentEtaPos);
  fAvgCos2PsivsCentEtaNeg = new TProfile("fAvgCos2PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgCos2PsivsCentEtaNeg);
  fAvgSin2PsivsCentEtaNeg = new TProfile("fAvgSin2PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgSin2PsivsCentEtaNeg);

  fAvgCos3PsivsCentEtaPos = new TProfile("fAvgCos3PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgCos3PsivsCentEtaPos);
  fAvgSin3PsivsCentEtaPos = new TProfile("fAvgSin3PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgSin3PsivsCentEtaPos);
  fAvgCos3PsivsCentEtaNeg = new TProfile("fAvgCos3PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgCos3PsivsCentEtaNeg);
  fAvgSin3PsivsCentEtaNeg = new TProfile("fAvgSin3PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgSin3PsivsCentEtaNeg);  

  fAvgCos4PsivsCentEtaPos = new TProfile("fAvgCos4PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",4),90,0,90);
  fListHist->Add(fAvgCos4PsivsCentEtaPos);
  fAvgSin4PsivsCentEtaPos = new TProfile("fAvgSin4PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",4),90,0,90);
  fListHist->Add(fAvgSin4PsivsCentEtaPos);
  fAvgCos4PsivsCentEtaNeg = new TProfile("fAvgCos4PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",4),90,0,90);
  fListHist->Add(fAvgCos4PsivsCentEtaNeg);
  fAvgSin4PsivsCentEtaNeg = new TProfile("fAvgSin4PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",4),90,0,90);
  fListHist->Add(fAvgSin4PsivsCentEtaNeg);





  /////////   LAMBDA Histograms   //////////


  /////--------> chunzheng: V0 QA and other histograms:
  //V0-mother informations:
  fHistV0Pt = new TH1D("hV0Pt","", 200, 0., 20.);
  fListHist->Add(fHistV0Pt);
  fHistV0Eta = new TH1D("hV0Eta","", 200, -10., 10.);
  fListHist->Add(fHistV0Eta);
  fHistV0DcatoPrimVertex = new TH1D("hV0DcaToPrimVertex","",200, 0., 20.);
  fListHist->Add(fHistV0DcatoPrimVertex);
  fHistV0CPA = new TH1D("hV0CPA","", 1000, 0.9, 1.);
  fListHist->Add(fHistV0CPA);
  fHistV0DecayLength = new TH1D("hV0DecayLength","",500,0,500.);
  fListHist->Add(fHistV0DecayLength);


  Char_t name[10];

  for(int i=0; i<2; i++){
    ///// Case 0 = before cut, case 1 = afterCut.
    if(i==0)
      sprintf(name,"Before");
    else
       sprintf(name,"After");
    /// Lambdas:
    fHistLambdaPt[i] = new TH1D(Form("hLambdaPt_%sMassCut",name),"", 200, 0., 20.);
    fListHist->Add(fHistLambdaPt[i]);
    fHistLambdaEta[i] = new TH1D(Form("hLambdaEta_%sMassCut",name),"",200, -10., 10.);
    fListHist->Add(fHistLambdaEta[i]);
    fHistLambdaPhi[i] = new TH1D(Form("hLambdaPhi_%sMassCut",name),"", 360, 0., TMath::TwoPi());
    fListHist->Add(fHistLambdaPhi[i]);
    fHistLambdaDcaToPrimVertex[i] = new TH1D(Form("hLambdaDcaToPrimVertex_%sMassCut",name),"",200, 0., 20.);
    fListHist->Add(fHistLambdaDcaToPrimVertex[i]);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",name),"",200, 0.9, 1.);
    fListHist->Add(fHistLambdaCPA[i]);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",name),"", 250, 0., 500.);
    fListHist->Add(fHistLambdaDecayLength[i]);
    fHistLambdaMass[i] = new TH1D(Form("hLambdaMass_%sMassCut",name),"",1000,1.,1.25); //  Current bin size = 0.00025
    fListHist->Add(fHistLambdaMass[i]);
    fProfileLambdaMassVsPt[i] = new TProfile(Form("pLambdaMassVsPt_%sMassCut",name),"",200,0,20);
    fListHist->Add(fProfileLambdaMassVsPt[i]);

    // AntiLambdas
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",name),"", 200, 0., 20.);
    fListHist->Add(fHistAntiLambdaPt[i]);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",name),"",200, -10., 10.);
    fListHist->Add(fHistAntiLambdaEta[i]);
    fHistAntiLambdaPhi[i] = new TH1D(Form("hAntiLambdaPhi_%sMassCut",name),"", 360, 0, TMath::TwoPi());
    fListHist->Add(fHistAntiLambdaPhi[i]);
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",name),"",200, 0., 20.);
    fListHist->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",name),"",200, 0.9, 1.);
    fListHist->Add(fHistAntiLambdaCPA[i]);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",name),"", 250, 0., 500.);
    fListHist->Add(fHistAntiLambdaDecayLength[i]);
    fHistAntiLambdaMass[i] = new TH1D(Form("hAntiLambdaMass_%sMassCut",name),"",1000,1.,1.25); // Current bin size = 0.00025
    fListHist->Add(fHistAntiLambdaMass[i]);
    fProfileAntiLambdaMassVsPt[i] = new TProfile(Form("pAntiLambdaMassVsPt_%sMassCut",name),"",200,0,20);
    fListHist->Add(fProfileAntiLambdaMassVsPt[i]);    
  }


  ///Gamma Correlators:
  //TPC Plane
  ///Lambda - X
  fProfileGammaTPC_Lambda_hPos = new TProfile("fProfileGammaTPC_Lambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_Lambda_hPos);
  fProfileGammaTPC_Lambda_hNeg = new TProfile("fProfileGammaTPC_Lambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_Lambda_hNeg);
  fProfileGammaTPC_Lambda_Proton = new TProfile("fProfileGammaTPC_Lambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_Lambda_Proton);
  fProfileGammaTPC_Lambda_AntiProton = new TProfile("fProfileGammaTPC_Lambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaTPC_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaTPC_AntiLambda_hPos = new TProfile("fProfileGammaTPC_AntiLambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_AntiLambda_hPos);
  fProfileGammaTPC_AntiLambda_hNeg = new TProfile("fProfileGammaTPC_AntiLambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_AntiLambda_hNeg);
  fProfileGammaTPC_AntiLambda_Proton = new TProfile("fProfileGammaTPC_AntiLambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaTPC_AntiLambda_Proton);
  fProfileGammaTPC_AntiLambda_AntiProton = new TProfile("fProfileGammaTPC_AntiLambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaTPC_AntiLambda_AntiProton);

  //V0C Plane
  ///Lambda - X
  fProfileGammaV0C_Lambda_hPos = new TProfile("fProfileGammaV0C_Lambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_Lambda_hPos);
  fProfileGammaV0C_Lambda_hNeg = new TProfile("fProfileGammaV0C_Lambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_Lambda_hNeg);
  fProfileGammaV0C_Lambda_Proton = new TProfile("fProfileGammaV0C_Lambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_Lambda_Proton);
  fProfileGammaV0C_Lambda_AntiProton = new TProfile("fProfileGammaV0C_Lambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaV0C_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaV0C_AntiLambda_hPos = new TProfile("fProfileGammaV0C_AntiLambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_AntiLambda_hPos);
  fProfileGammaV0C_AntiLambda_hNeg = new TProfile("fProfileGammaV0C_AntiLambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_AntiLambda_hNeg);
  fProfileGammaV0C_AntiLambda_Proton = new TProfile("fProfileGammaV0C_AntiLambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0C_AntiLambda_Proton);
  fProfileGammaV0C_AntiLambda_AntiProton = new TProfile("fProfileGammaV0C_AntiLambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaV0C_AntiLambda_AntiProton);

  //V0A Plane
  ///Lambda - X
  fProfileGammaV0A_Lambda_hPos = new TProfile("fProfileGammaV0A_Lambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_Lambda_hPos);
  fProfileGammaV0A_Lambda_hNeg = new TProfile("fProfileGammaV0A_Lambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_Lambda_hNeg);
  fProfileGammaV0A_Lambda_Proton = new TProfile("fProfileGammaV0A_Lambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_Lambda_Proton);
  fProfileGammaV0A_Lambda_AntiProton = new TProfile("fProfileGammaV0A_Lambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaV0A_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaV0A_AntiLambda_hPos = new TProfile("fProfileGammaV0A_AntiLambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_AntiLambda_hPos);
  fProfileGammaV0A_AntiLambda_hNeg = new TProfile("fProfileGammaV0A_AntiLambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_AntiLambda_hNeg);
  fProfileGammaV0A_AntiLambda_Proton = new TProfile("fProfileGammaV0A_AntiLambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileGammaV0A_AntiLambda_Proton);
  fProfileGammaV0A_AntiLambda_AntiProton = new TProfile("fProfileGammaV0A_AntiLambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileGammaV0A_AntiLambda_AntiProton);

  //ZNC Plane
  ///Lambda - X
  fProfileGammaZNC_Lambda_hPos = new TProfile("fProfileGammaZNC_Lambda_hPos","",18,0,90);
  fListHist->Add(fProfileGammaZNC_Lambda_hPos);
  fProfileGammaZNC_Lambda_hNeg = new TProfile("fProfileGammaZNC_Lambda_hNeg","",18,0,90);
  fListHist->Add(fProfileGammaZNC_Lambda_hNeg);
  fProfileGammaZNC_Lambda_Proton = new TProfile("fProfileGammaZNC_Lambda_Proton","",18,0,90);
  fListHist->Add(fProfileGammaZNC_Lambda_Proton);
  fProfileGammaZNC_Lambda_AntiProton = new TProfile("fProfileGammaZNC_Lambda_AntiProton","",18,0,90);
  fListHist->Add(fProfileGammaZNC_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaZNC_AntiLambda_hPos = new TProfile("fProfileGammaZNC_AntiLambda_hPos","",18,0,90);
  fListHist->Add(fProfileGammaZNC_AntiLambda_hPos);
  fProfileGammaZNC_AntiLambda_hNeg = new TProfile("fProfileGammaZNC_AntiLambda_hNeg","",18,0,90);
  fListHist->Add(fProfileGammaZNC_AntiLambda_hNeg);
  fProfileGammaZNC_AntiLambda_Proton = new TProfile("fProfileGammaZNC_AntiLambda_Proton","",18,0,90);
  fListHist->Add(fProfileGammaZNC_AntiLambda_Proton);
  fProfileGammaZNC_AntiLambda_AntiProton = new TProfile("fProfileGammaZNC_AntiLambda_AntiProton","",18,0,90);
  fListHist->Add(fProfileGammaZNC_AntiLambda_AntiProton);

  //ZNA Plane
  ///Lambda - X
  fProfileGammaZNA_Lambda_hPos = new TProfile("fProfileGammaZNA_Lambda_hPos","",18,0,90);
  fListHist->Add(fProfileGammaZNA_Lambda_hPos);
  fProfileGammaZNA_Lambda_hNeg = new TProfile("fProfileGammaZNA_Lambda_hNeg","",18,0,90);
  fListHist->Add(fProfileGammaZNA_Lambda_hNeg);
  fProfileGammaZNA_Lambda_Proton = new TProfile("fProfileGammaZNA_Lambda_Proton","",18,0,90);
  fListHist->Add(fProfileGammaZNA_Lambda_Proton);
  fProfileGammaZNA_Lambda_AntiProton = new TProfile("fProfileGammaZNA_Lambda_AntiProton","",18,0,90);
  fListHist->Add(fProfileGammaZNA_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaZNA_AntiLambda_hPos = new TProfile("fProfileGammaZNA_AntiLambda_hPos","",18,0,90);
  fListHist->Add(fProfileGammaZNA_AntiLambda_hPos);
  fProfileGammaZNA_AntiLambda_hNeg = new TProfile("fProfileGammaZNA_AntiLambda_hNeg","",18,0,90);
  fListHist->Add(fProfileGammaZNA_AntiLambda_hNeg);
  fProfileGammaZNA_AntiLambda_Proton = new TProfile("fProfileGammaZNA_AntiLambda_Proton","",18,0,90);
  fListHist->Add(fProfileGammaZNA_AntiLambda_Proton);
  fProfileGammaZNA_AntiLambda_AntiProton = new TProfile("fProfileGammaZNA_AntiLambda_AntiProton","",18,0,90);
  fListHist->Add(fProfileGammaZNA_AntiLambda_AntiProton);
  
  ///Delta Correlators:
   ///Lambda - X
  fProfileDelta_Lambda_hPos = new TProfile("fProfileDelta_Lambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_Lambda_hPos);
  fProfileDelta_Lambda_hNeg = new TProfile("fProfileDelta_Lambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_Lambda_hNeg);
  fProfileDelta_Lambda_Proton = new TProfile("fProfileDelta_Lambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_Lambda_Proton);
  fProfileDelta_Lambda_AntiProton = new TProfile("fProfileDelta_Lambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileDelta_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileDelta_AntiLambda_hPos = new TProfile("fProfileDelta_AntiLambda_hPos", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_AntiLambda_hPos);
  fProfileDelta_AntiLambda_hNeg = new TProfile("fProfileDelta_AntiLambda_hNeg", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_AntiLambda_hNeg);
  fProfileDelta_AntiLambda_Proton = new TProfile("fProfileDelta_AntiLambda_Proton", "", 18, 0., 90.);
  fListHist->Add(fProfileDelta_AntiLambda_Proton);
  fProfileDelta_AntiLambda_AntiProton = new TProfile("fProfileDelta_AntiLambda_AntiProton", "", 18, 0., 90.);  
  fListHist->Add(fProfileDelta_AntiLambda_AntiProton);

  ///// <-------- chunzheng's Histograms....
  
  
}





Bool_t AliAnalysisTaskGammaDeltaPID::IsGoodV0(AliAODv0 *aodV0) 
{
  if (!aodV0) {
    AliError(Form("ERROR: Could not retrieve aodV0"));
    return kFALSE;
  }
  // Offline reconstructed V0 only
  if ( aodV0->GetOnFlyStatus() ) return kFALSE;
  // Get daughters and check them     
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return kFALSE;
  }
  // Unlike signs of daughters
  if ( myTrackNegTest->Charge() == myTrackPosTest->Charge() ) return kFALSE;
  // Cosinus of pointing angle      
  double dCPA = aodV0->CosPointingAngle(fCurrentVtx);
  // cut on Cosinus of pointing angle
  if ( dCPA < fV0CPAMin ) return kFALSE;
  // DCA of V0
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if ( TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax ) return kFALSE;
  // V0 path length before decay
  double dDecayLength = aodV0->DecayLengthV0(fCurrentVtx);
  if ( dDecayLength > fV0DecayLengthMax ) return kFALSE;
  if ( dDecayLength < fV0DecayLengthMin ) return kFALSE;
  // DCA between daughters
  double dDCA = aodV0->DcaV0Daughters();
  if ( dDCA > fV0DcaBetweenDaughtersMax ) return kFALSE;
  double dPt = aodV0->Pt();
  if( dPt < fV0PtMin ) return kFALSE;
  double dRapidity = aodV0->RapLambda();
  if( TMath::Abs(dRapidity) > fV0RapidityMax ) return kFALSE;
  return kTRUE;
}



Bool_t AliAnalysisTaskGammaDeltaPID::IsGoodDaughterTrack(const AliAODTrack *track)
{
  // TPC refit
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  //No kinks
  if(Int_t(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return kFALSE;
  // Maximum value of transverse momentum 
  double dPt = track->Pt();
  if (dPt > fDaughtersPtMax) return kFALSE;
  // Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (TMath::Abs(dEta) > fDaughtersEtaMax) return kFALSE;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return kFALSE;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  // [number of crossed rows]>0.8  [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  return kTRUE;
}

Int_t AliAnalysisTaskGammaDeltaPID::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack) 
{

  Bool_t IsLambda     = kFALSE;
  Bool_t IsAntiLambda = kFALSE;
  Int_t  code = 0;

  //Λ-->(p+)+(π-)
  Float_t nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton));//TPC p+
  Float_t nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion));//TPC π-
  //(Λ-)-->(p-)+(π+)
  Float_t nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion));//TPC π+
  Float_t nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton));//TPC p-

  IsLambda = (nSigTPCPosProton < fV0PosProtonTPCNsigma) && (nSigTPCNegPion < fV0NegPionTPCNsigma);
  IsAntiLambda = (nSigTPCNegProton < fV0NegProtonTPCNsigma) && (nSigTPCPosPion < fV0PosPionTPCNsigma);

  if (fV0DaughterUseTOF)
  {
    Float_t nSigTOFPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack,AliPID::kProton));//TOF p+
    Float_t nSigTOFNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack,AliPID::kPion));//TOF π-
    Float_t nSigTOFPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack,AliPID::kPion));//TOF π+
    Float_t nSigTOFNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack,AliPID::kProton));//TOF p-

    IsLambda *= (nSigTOFPosProton < fV0PosProtonTOFNsigma && nSigTOFNegPion < fV0NegPionTOFNsigma);
    IsAntiLambda *= (nSigTOFNegProton < fV0NegProtonTOFNsigma && nSigTOFPosPion < fV0PosPionTOFNsigma);
  }
  
  if (IsLambda)     code =  3122;
  if (IsAntiLambda) code = -3122;
  if (IsLambda && IsAntiLambda) code = 0;

  return code;
}





void AliAnalysisTaskGammaDeltaPID::SetupPileUpRemovalFunctions(){
  
  ////==========> LHC18q/r PileUp Removal Functions: ---- Do not Remove them !!! -----
  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
  
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
  
  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",  0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
  //--------------------------------------------------------------------------------------

}













void AliAnalysisTaskGammaDeltaPID::SetupEventAndTaskConfigInfo(){

 
  ///--  Analysis cuts and settings: -------
  fHistAnalysisInfo = new TH1F("fHistAnalysisInfo","Analysis Settings and cuts",40, 0., 40);  
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(1,"VzCutLow");
  fHistAnalysisInfo->SetBinContent(1,fMinVzCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(2,"VzCutHigh");
  fHistAnalysisInfo->SetBinContent(2,fMaxVzCut);
  if(sCentrEstimator=="V0M"){
    fHistAnalysisInfo->GetXaxis()->SetBinLabel(3,"CentV0M");
    fHistAnalysisInfo->SetBinContent(3,1);
  }
  else if(sCentrEstimator=="CL0"){
    fHistAnalysisInfo->GetXaxis()->SetBinLabel(4,"CentCL0");
    fHistAnalysisInfo->SetBinContent(4,1);
  }
  else if(sCentrEstimator=="CL1"){
    fHistAnalysisInfo->GetXaxis()->SetBinLabel(4,"CentCL1");
    fHistAnalysisInfo->SetBinContent(4,1);
  }  
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(5,"PileUpCut");
  if(!bSkipPileUpCut){
    fHistAnalysisInfo->SetBinContent(5,1);
  }
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(6,"#Psi harmonic");
  fHistAnalysisInfo->SetBinContent(6,gHarmonic);  

  fHistAnalysisInfo->GetXaxis()->SetBinLabel(7,"NUA Applied");
  if(fListNUACorr){
    fHistAnalysisInfo->SetBinContent(7,1);
  }  
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(8,"NUE Applied");
  if(fListTRKCorr){
    fHistAnalysisInfo->SetBinContent(8,1);
  }
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(9,"V0 Gain/Q ");
  if(fListV0MCorr){
    fHistAnalysisInfo->SetBinContent(9,1);
  }
  

  fHistAnalysisInfo->GetXaxis()->SetBinLabel(10,"ZDC Gain/Q "); 
  //if(fListV0MCorr){
  //fHistAnalysisInfo->SetBinContent(10,1);
  //}
  
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(11,"Reserved");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(12,"Reserved");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(13,"Reserved");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(14,"Reserved");


  fHistAnalysisInfo->GetXaxis()->SetBinLabel(15,"FilterBit");
  fHistAnalysisInfo->SetBinContent(15,fFilterBit);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(16,"N_{cls,TPC}");
  fHistAnalysisInfo->SetBinContent(16,fTPCclustMin);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(17,"Chi^{2}max");
  fHistAnalysisInfo->SetBinContent(17,fTrkChi2Max);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(18,"PtLowAna");
  fHistAnalysisInfo->SetBinContent(18,fMinPtCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(19,"PtHighAna");
  fHistAnalysisInfo->SetBinContent(19,fMaxPtCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(20,"EtaLowAna");    
  fHistAnalysisInfo->SetBinContent(20,fMinEtaCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(21,"EtaHighAna");
  fHistAnalysisInfo->SetBinContent(21,fMaxEtaCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(22,"ParticlePairID");
  fHistAnalysisInfo->SetBinContent(22,gParticleID);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(23,"nSigTPC");
  fHistAnalysisInfo->SetBinContent(23,fNSigmaTPCCut);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(24,"nSigCircular");
  fHistAnalysisInfo->SetBinContent(24,fNSigmaTOFCut);
    
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(25,"Reserved");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(26,"Reserved");
  
  fListHist->Add(fHistAnalysisInfo);
}

void AliAnalysisTaskGammaDeltaPID::SetupLambdaPIDFlowHistograms() {
  // Double_t centbins[10] = {0., 10., 20., 30., 40., 50., 60., 70, 80, 90};
  // Double_t pTbins[13]   = {0.0, 0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 5.0};
  // Double_t phibins[361] = {0.};
  // for (Int_t i = 0; i <= 360; i++) phibins[i] = i * 4.*TMath::Pi()/360.;
  Char_t name1[20];
  for (Int_t i = 0; i < 5; i++) {
    if(i==0) sprintf(name1,"TPC");
    if(i==1) sprintf(name1,"V0C");
    if(i==2) sprintf(name1,"V0A");
    if(i==3) sprintf(name1,"ZNC");
    if(i==4) sprintf(name1,"ZNA");
    fProfile2DRawFlowCentPthPos[i] = new TProfile2D(Form("fProfile2DRawFlowCentPthPos_%s",name1),Form("fProfile2DRawFlowCentPthPos_%s",name1), 10, 0., 100., 25, 0., 5.);
    fProfile2DRawFlowCentPthNeg[i] = new TProfile2D(Form("fProfile2DRawFlowCentPthNeg_%s",name1),Form("fProfile2DRawFlowCentPthNeg_%s",name1), 10, 0., 100., 25, 0., 5.);
    fProfile2DRawFlowCentPtProton[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtProton_%s",name1),Form("fProfile2DRawFlowCentPtProton_%s",name1), 10, 0., 100., 25, 0., 5.);
    fProfile2DRawFlowCentPtAntiProton[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtAntiProton_%s",name1),Form("fProfile2DRawFlowCentPtAntiProtong_%s",name1), 10, 0., 100., 25, 0., 5.);
    fListHist->Add(fProfile2DRawFlowCentPthPos[i]);
    fListHist->Add(fProfile2DRawFlowCentPthNeg[i]);
    fListHist->Add(fProfile2DRawFlowCentPtProton[i]);
    fListHist->Add(fProfile2DRawFlowCentPtAntiProton[i]);
  }
  Char_t name2[20];
  for (Int_t i = 0; i < 6; i++) {
    if(i==0) sprintf(name2,"TPCPos");
    if(i==1) sprintf(name2,"TPCNeg");
    if(i==2) sprintf(name2,"V0C");
    if(i==3) sprintf(name2,"V0A");
    if(i==4) sprintf(name2,"ZNC");
    if(i==5) sprintf(name2,"ZNA");
    fProfile2DRawFlowCentPtLambda[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtLambda_%s",name2),Form("fProfile2DRawFlowCentPtLambda_%s",name2), 10, 0., 100., 25, 0., 5.);
    fProfile2DRawFlowCentPtAntiLambda[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtAntiLambda_%s",name2),Form("fProfile2DRawFlowCentPtAntiLambda_%s",name2), 10, 0., 100., 25, 0., 5.);
    fListHist->Add(fProfile2DRawFlowCentPtLambda[i]);
    fListHist->Add(fProfile2DRawFlowCentPtAntiLambda[i]);
  }
}



void AliAnalysisTaskGammaDeltaPID::GetNUACorrectionHist(Int_t run, Int_t kParticleID) {
  if(fListNUACorr){
    //charge:
    fHCorrectNUAChrgPos = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,run));  // 0 = Charge
    fHCorrectNUAChrgNeg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,run));  // 0 = Charge
    ///PID pion=1,Kaon=2,Prot=3
    if(kParticleID==10) kParticleID = 1;
    if(kParticleID==20) kParticleID = 2;
    if(kParticleID==30) kParticleID = 3;
    
    fHCorrectNUAkPIDPos = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",kParticleID,run)); 
    fHCorrectNUAkPIDNeg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",kParticleID,run)); 
    
    if(fHCorrectNUAChrgPos && fHCorrectNUAChrgNeg){
      cout<<" ===========> Info:: Setting up NUA corrections for Charge, run = "<<run<<endl;   
    }
    if(fHCorrectNUAkPIDPos && fHCorrectNUAkPIDNeg){
      cout<<" ===========> Info:: Setting up NUA corrections for PID, run = "<<run<<endl;   
    }    
    /// Now Get Average Q vector:
    fHCorrectTPCQnxEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCosNPsivsCentEtaPosRun%d",run));
    fHCorrectTPCQnyEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSinNPsivsCentEtaPosRun%d",run));
    fHCorrectTPCQnxEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCosNPsivsCentEtaNegRun%d",run));
    fHCorrectTPCQnyEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSinNPsivsCentEtaNegRun%d",run));

    fHCorrectTPCQ3xEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCos3PsivsCentEtaPosRun%d",run));
    fHCorrectTPCQ3yEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSin3PsivsCentEtaPosRun%d",run));
    fHCorrectTPCQ3xEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCos3PsivsCentEtaNegRun%d",run));
    fHCorrectTPCQ3yEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSin3PsivsCentEtaNegRun%d",run));


    ///For debug purpose:
    if(fHCorrectTPCQnxEtaPos && fHCorrectTPCQnyEtaPos){
      cout<<" ===========> Info:: Found TPC <Q2> vectors for EtaPos run "<<run<<endl;   
    }
    if(fHCorrectTPCQnxEtaNeg && fHCorrectTPCQnyEtaNeg){
      cout<<" ===========> Info:: Found TPC <Q2> vectors for EtaNeg run "<<run<<endl;   
    }
    if(fHCorrectTPCQ3xEtaPos && fHCorrectTPCQ3yEtaPos){
      cout<<" ===========> Info:: Found TPC <Q3> vectors for EtaPos run "<<run<<endl;   
    }
    if(fHCorrectTPCQ3xEtaNeg && fHCorrectTPCQ3yEtaNeg){
      cout<<" ===========> Info:: Found TPC <Q3> vectors for EtaNeg run "<<run<<endl;   
    }
    
    
  }
  else {
    printf("\n ******** ::Warning:: ********* \n No NUA Correction Histograms found for run %d, Using Wgts = 1.0 \n",run);
  }
}

void AliAnalysisTaskGammaDeltaPID::GetMCCorrectionHist(){

    if(fListTRKCorr) {
    /// Default: Centrality Independent MC efficiency:    
      fHCorrectMCPosChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
      fHCorrectMCPosPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionPos");
      fHCorrectMCPosKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonPos");
      fHCorrectMCPosProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtPos");

      fHCorrectMCNegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");
      fHCorrectMCNegPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionNeg");
      fHCorrectMCNegKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonNeg");
      fHCorrectMCNegProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtNeg");
    }
    // Debug:
    if(fHCorrectMCPosChrg && fHCorrectMCNegChrg){
      printf(" ===========> Info:: Tracking Efficiency Found for Charge \n ");
    }
    if(fHCorrectMCPosPion && fHCorrectMCNegPion){
      printf(" ===========> Info:: Tracking Efficiency Found for Pions \n ");
    }
    if(fHCorrectMCPosKaon && fHCorrectMCNegKaon){
      printf(" ===========> Info:: Tracking Efficiency Found for Kaons \n ");
    }
    if(fHCorrectMCPosProt && fHCorrectMCNegProt){
      printf(" ===========> Info:: Tracking Efficiency Found for Proton \n ");
    }
}

void AliAnalysisTaskGammaDeltaPID::GetV0MCorrectionHist(Int_t run){ 

  if(fListV0MCorr){
    //V0 Channel Gains:
    fHCorrectV0ChWeghts = (TH2F *) fListV0MCorr->FindObject(Form("hWgtV0ChannelsvsVzRun%d",run));
    if(fHCorrectV0ChWeghts){
      printf("\n ===========> Info:: V0 Channel Weights Found for Run %d \n ",run);
    }
    //Get V0A, V0C <Q> Vectors:
    fHCorrectQNxV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0CRun%d",run));
    fHCorrectQNyV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0CRun%d",run));    
    fHCorrectQNxV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0ARun%d",run));
    fHCorrectQNyV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0ARun%d",run));
	
    fHCorrectQ3xV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0CRun%d",run));
    fHCorrectQ3yV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0CRun%d",run));    
    fHCorrectQ3xV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0ARun%d",run));
    fHCorrectQ3yV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0ARun%d",run));    
    if(fHCorrectQNxV0C && fHCorrectQNyV0C && fHCorrectQNxV0A && fHCorrectQNyV0A){
      printf(" ===========> Info:: V0A,V0C <Q> Found for Run %d \n ",run);
    }    
  }
  else{
    fHCorrectV0ChWeghts=NULL;
  } 
}

void AliAnalysisTaskGammaDeltaPID::GetZDCCorrectionHist(Int_t run){ 
  if(fListZDCCorr){
	fHZDCCparameters = (TH1D*)(fListZDCCorr->FindObject(Form("Run %d", run))->FindObject(Form("fZDCCparameters[%d]",run)));
	fHZDCAparameters = (TH1D*)(fListZDCCorr->FindObject(Form("Run %d", run))->FindObject(Form("fZDCAparameters[%d]",run)));
	if(fHZDCCparameters && fHZDCAparameters){
      printf("\n ===========> Info:: ZDC Channel Weights Found for Run %d \n ",run);
    }
  }
  else{
	fHZDCCparameters=NULL;
	fHZDCAparameters=NULL;
	printf("\n ===========> Info:: ZDC Channel Weights NOT Found for Run %d \n ",run);
  }
}

Bool_t AliAnalysisTaskGammaDeltaPID::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){
  
  if(pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check
  
  Bool_t bPIDokay = kFALSE;

  if(!fPIDResponse){
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..
  
  Double_t nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  Double_t trkPtPID  = ftrack->Pt();
  Int_t trkChargePID = ftrack->Charge();

  ///Pion => 
  if(pidToCheck==1){ 
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);//Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if(trkPtPID<=0.5 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut){
      bPIDokay = kTRUE;
    }
    // Using TPCTOF RMS cut for higher pt: 
    else if(trkPtPID>0.5 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut){ 
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  ///Kaon => 
  else if(pidToCheck==2){
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);	
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
	
    if(trkPtPID<=0.45 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut){
      bPIDokay = kTRUE;
    }
    else if(trkPtPID>0.45 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut){
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  ///proton => 
  else if(pidToCheck==3){///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);    
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
	
    if(trkPtPID<=0.6 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut){
      bPIDokay = kTRUE;
      if(trkChargePID>0 && trkPtPID<0.4) bPIDokay = kFALSE;  
    }
    else if(trkPtPID>0.6 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut){
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  ///Lambda => 
  else if(pidToCheck==4){///not defined yet. 
    return kFALSE;
    // Rihan: call Lamda PID function from Chunseng..
  }
  else{
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3,4  (Charge Pion, Kaon, Proton, Lambda)\n return with kFALSE \n");
    return kFALSE;
  }
  
  return kFALSE;
}



Bool_t AliAnalysisTaskGammaDeltaPID::CheckEventIsPileUp2018(AliAODEvent *faod) {

  /// Todo Rihan: I can check for PileUp and get TPC event Plane in Same Function
  /// Utilizing same track loop. This method would save time..
  
  Bool_t BisPileup=kFALSE;

  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = faod->GetNumberOfTracks();

  Int_t multTrk = 0;

  for (Int_t it = 0; it < nTracks; it++) {
    
    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);

    if (!aodTrk){
      delete aodTrk;
      continue;
    }

    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
    }
  }

  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);

  if (centrCL0 < fCenCutLowPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (centrCL0 > fCenCutHighPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) {
    BisPileup=kTRUE;
  }     
  if (multV0On < fV0CutPU->Eval(multV0Tot)) {
    BisPileup=kTRUE;
  }
  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    BisPileup=kTRUE;
  }
  if (faod->IsIncompleteDAQ()) {
    BisPileup=kTRUE;
  }    
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  fHistCentCL0VsV0MBefore->Fill(centrV0M,centrCL0);
  fHistTPCVsESDTrkBefore->Fill(multTrk,multEsd);  
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTrk);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);

  if (!BisPileup)
    {      
      fHistCentCL0VsV0MAfter->Fill(centrV0M,centrCL0);
      fHistTPCVsESDTrkAfter->Fill(multTrk,multEsd);  
      fHistTPConlyVsCL1After->Fill(centrCL1,multTrk);
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);      
    }

  return BisPileup; 
}



Bool_t AliAnalysisTaskGammaDeltaPID::GetTPCQvectAndRemovePileUp2018(AliAODEvent *faod,Double_t *qnxEtaNeg,Double_t *qnyEtaNeg,Double_t *qnxEtaPos,Double_t *qnyEtaPos,Double_t& multNeg,Double_t& multPos)
{
  Bool_t BisPileup=kFALSE;
  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = faod->GetNumberOfTracks();
  Int_t multTrk = 0;


  
  Int_t    fTPCclustMinforEP= 70;      // Fixed for EP calculation
  Int_t    trk1Chrg=0,trk1NClus=0;     //
  Double_t fMinPtCutforEP   = 0.2;    // Fixed for EP calculation
  Double_t fMaxPtCutforEP   = 2.0;    // Fixed for EP calculation
  Double_t fEtaGapPosforEP  = 0.1;    // could be made variable in AddTask Macro.
  Double_t fEtaGapNegforEP  =-0.1;    // could be made variable in AddTask Macro.

  Double_t trk1Pt=0, trk1Phi=0, trk1Eta=0, trk1Chi2=0, trk1dEdx=0, trk1Wgt=1.0;  
  Double_t SumQ2xTPCPos = 0., SumQ2yTPCPos = 0., SumQ2xTPCNeg = 0., SumQ2yTPCNeg = 0;
  Double_t SumQ3xTPCPos = 0., SumQ3yTPCPos = 0., SumQ3xTPCNeg = 0., SumQ3yTPCNeg = 0;
  Double_t SumQ4xTPCPos = 0., SumQ4yTPCPos = 0., SumQ4xTPCNeg = 0., SumQ4yTPCNeg = 0;
  
  Double_t fWgtMultTPCPos=0., fWgtMultTPCNeg=0;

  
    
  for (Int_t it = 0; it < nTracks; it++) {
    
    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);

    if (!aodTrk){
      //delete aodTrk;
      continue;
    }
    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
    }
    
    /// Now get TPC Q-Vectors:
    if(!aodTrk->TestFilterBit(fFilterBit))    continue;  //// Only use FB tracks. 

    trk1Pt    = aodTrk->Pt();
    trk1Phi   = aodTrk->Phi();
    trk1Eta   = aodTrk->Eta();
    trk1Chrg  = aodTrk->Charge();
    trk1Chi2  = aodTrk->Chi2perNDF();
    trk1NClus = aodTrk->GetTPCNcls();           
    trk1dEdx  = aodTrk->GetDetPid()->GetTPCsignal();  
  
    //trk1DCAxy = aodTrk->DCA();
    //trk1DCAz  = aodTrk->ZAtDCA();
      
    //Apply track cuts for EP  here:
    if(trk1NClus< fTPCclustMinforEP)        continue;
    if((trk1Eta < -0.8) || (trk1Eta > 0.8)) continue;    
    if(trk1Pt < fMinPtCutforEP) continue;
    if(trk1Pt > fMaxPtCutforEP) continue;
    if(trk1Chi2 < 0.1)          continue;
    if(trk1Chi2 > 4.0)          continue;
    if(trk1dEdx < 10)           continue;
    if(!TMath::Abs(trk1Chrg))   continue;

    Int_t trk1ID = aodTrk->GetID();         
      
    trk1Wgt = GetNUAWeightForTrack(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg);
    
    ///Used Pt as weight for Better resolution:
    if(trk1Eta >= fEtaGapPosforEP){
      SumQ2xTPCPos   += trk1Wgt*TMath::Cos(2*trk1Phi);
      SumQ2yTPCPos   += trk1Wgt*TMath::Sin(2*trk1Phi);
      SumQ3xTPCPos   += trk1Wgt*TMath::Cos(3*trk1Phi);
      SumQ3yTPCPos   += trk1Wgt*TMath::Sin(3*trk1Phi);
      SumQ4xTPCPos   += trk1Wgt*TMath::Cos(4*trk1Phi);
      SumQ4yTPCPos   += trk1Wgt*TMath::Sin(4*trk1Phi);
      
      fWgtMultTPCPos += trk1Wgt; 
      vecPosEPTrkID.push_back(trk1ID);
      vecPosEPTrkPhi.push_back(trk1Phi);
      vecPosEPTrkNUAWgt.push_back(trk1Wgt);
    }
    else if(trk1Eta <= fEtaGapNegforEP){
      SumQ2xTPCNeg   += trk1Wgt*TMath::Cos(2*trk1Phi);
      SumQ2yTPCNeg   += trk1Wgt*TMath::Sin(2*trk1Phi);
      SumQ3xTPCNeg   += trk1Wgt*TMath::Cos(3*trk1Phi);
      SumQ3yTPCNeg   += trk1Wgt*TMath::Sin(3*trk1Phi);
      SumQ4xTPCNeg   += trk1Wgt*TMath::Cos(4*trk1Phi);
      SumQ4yTPCNeg   += trk1Wgt*TMath::Sin(4*trk1Phi);      
      
      fWgtMultTPCNeg += trk1Wgt;
      vecNegEPTrkID.push_back(trk1ID);
      vecNegEPTrkPhi.push_back(trk1Phi);
      vecNegEPTrkNUAWgt.push_back(trk1Wgt);
    }
        
  }//AOD track loop



  
  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t  multV0a   = aodV0->GetMTotV0A(); // return total multiplicity in V0A 0-31
  Float_t  multV0c   = aodV0->GetMTotV0C(); // return total multiplicity in V0C 32-63
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA(); // Sum of the trigger (clock=10) charge on A side.
  UShort_t multV0cOn = aodV0->GetTriggerChargeC(); // Sum of the trigger (clock=10) charge on C side.
  UShort_t multV0On  = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot); // not needed?

  if (centrCL0 < fCenCutLowPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (centrCL0 > fCenCutHighPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) {
    BisPileup=kTRUE;
  }     
  if (multV0On < fV0CutPU->Eval(multV0Tot)) {
    BisPileup=kTRUE;
  }
  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    BisPileup=kTRUE;
  }
  if (faod->IsIncompleteDAQ()) {
    BisPileup=kTRUE;
  }    
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks(); // should we have a pileup cut on multESD and multTPC

  fHistCentCL0VsV0MBefore->Fill(centrV0M,centrCL0);
  fHistTPCVsESDTrkBefore->Fill(multTrk,multEsd);  
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTrk);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);

  if (!BisPileup) {      
    fHistCentCL0VsV0MAfter->Fill(centrV0M,centrCL0);
    fHistTPCVsESDTrkAfter->Fill(multTrk,multEsd);  
    fHistTPConlyVsCL1After->Fill(centrCL1,multTrk);
    fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);      
  }

  multNeg = fWgtMultTPCNeg;
  multPos = fWgtMultTPCPos;
  
  /// Set The q vector values for Event Plane: 
  if(fWgtMultTPCPos<0.1 || fWgtMultTPCNeg<0.1){        /// this means there is not enough tracks in this event!!

    qnxEtaNeg[0] = 0;
    qnxEtaNeg[1] = 0;
    qnxEtaNeg[2] = 0;    
    qnyEtaNeg[0] = 0;
    qnyEtaNeg[1] = 0;
    qnyEtaNeg[2] = 0;

    qnxEtaPos[0] = 0;
    qnxEtaPos[1] = 0;
    qnxEtaPos[2] = 0;
    qnyEtaPos[0] = 0;
    qnyEtaPos[1] = 0;
    qnyEtaPos[2] = 0;     
    
  }
  else{
    qnxEtaNeg[0] = SumQ2xTPCNeg;
    qnxEtaNeg[1] = SumQ3xTPCNeg;
    qnxEtaNeg[2] = SumQ4xTPCNeg;    
    qnyEtaNeg[0] = SumQ2yTPCNeg;
    qnyEtaNeg[1] = SumQ3yTPCNeg;
    qnyEtaNeg[2] = SumQ4yTPCNeg;

    qnxEtaPos[0] = SumQ2xTPCPos;
    qnxEtaPos[1] = SumQ3xTPCPos;
    qnxEtaPos[2] = SumQ4xTPCPos;    
    qnyEtaPos[0] = SumQ2yTPCPos;
    qnyEtaPos[1] = SumQ3yTPCPos;
    qnyEtaPos[2] = SumQ4yTPCPos;    
  }
  
  return BisPileup;  

}


Double_t AliAnalysisTaskGammaDeltaPID::GetNUAWeightForTrack(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

  Double_t WgtNUAtrk1 = 1.0;
  Int_t    iBinforNUA = 1;
  
  if(gChrg>0){
    if(fHCorrectNUAChrgPos){ /// safety measures for breaks!
      iBinforNUA = fHCorrectNUAChrgPos->FindBin(fVtxZ,fPhi,fEta);
      WgtNUAtrk1 = fHCorrectNUAChrgPos->GetBinContent(iBinforNUA);
    }
  }
  else{
    if(fHCorrectNUAChrgNeg){ /// safety measures for breaks!
      iBinforNUA = fHCorrectNUAChrgNeg->FindBin(fVtxZ,fPhi,fEta);
      WgtNUAtrk1 = fHCorrectNUAChrgNeg->GetBinContent(iBinforNUA);
    }  
  }
  return WgtNUAtrk1;
}

Double_t AliAnalysisTaskGammaDeltaPID::GetNUAWeightForTrackPID(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

  Double_t WgtNUAtrk1 = 1.0;
  Int_t    iBinforNUA = 0;

  /// It is responsibility of User to Set the Correct NUA file for the PID needed. 
  if(gChrg>0){
    if(fHCorrectNUAkPIDPos){ /// safety measures for breaks!
      iBinforNUA = fHCorrectNUAkPIDPos->FindBin(fVtxZ,fPhi,fEta);
      WgtNUAtrk1 = fHCorrectNUAkPIDPos->GetBinContent(iBinforNUA);
    }
  }
  else{
    if(fHCorrectNUAkPIDNeg){ /// safety measures for breaks!
      iBinforNUA = fHCorrectNUAkPIDNeg->FindBin(fVtxZ,fPhi,fEta);
      WgtNUAtrk1 = fHCorrectNUAkPIDNeg->GetBinContent(iBinforNUA);
    }  
  }
  
  return WgtNUAtrk1;
}



void AliAnalysisTaskGammaDeltaPID::ApplyTPCqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t& qxEtaNeg, Double_t& qyEtaNeg,Double_t& qxEtaPos,Double_t& qyEtaPos){

  Int_t icentbin = 0;

  if(gPsiN==2){ /// <Q2> correction for Psi2:  NOTE: Proper File should be set in AddTask!! 
    if(fHCorrectTPCQnxEtaPos && fHCorrectTPCQnyEtaPos){ ///prevents Code break!
      icentbin  = fHCorrectTPCQnxEtaPos->FindBin(fCent);
      qxEtaPos -= fHCorrectTPCQnxEtaPos->GetBinContent(icentbin); 
      qyEtaPos -= fHCorrectTPCQnyEtaPos->GetBinContent(icentbin);   /// Careful with Names => Qx and Qy DONT MISS-match!!!
      //cout<<" =======> Info::  Applied recenter for TPC Eta Pos Psi2 "<<endl;
    }
    if(fHCorrectTPCQnxEtaNeg && fHCorrectTPCQnyEtaNeg){ 
      icentbin  = fHCorrectTPCQnxEtaNeg->FindBin(fCent);
      qxEtaNeg -= fHCorrectTPCQnxEtaNeg->GetBinContent(icentbin);
      qyEtaNeg -= fHCorrectTPCQnyEtaNeg->GetBinContent(icentbin);   
      //cout<<" =======> Info::  Applied recenter for TPC Eta Neg Psi2 "<<endl;
    }
    
  }
  else if(gPsiN==3){  ///<Q3> correction for Psi3:  
    if(fHCorrectTPCQ3xEtaPos && fHCorrectTPCQ3yEtaPos){ 
      icentbin  = fHCorrectTPCQ3xEtaPos->FindBin(fCent);
      qxEtaPos -= fHCorrectTPCQ3xEtaPos->GetBinContent(icentbin); 
      qyEtaPos -= fHCorrectTPCQ3yEtaPos->GetBinContent(icentbin);
      //cout<<" =======> Info::  Applied recenter for TPC Eta Pos Psi3 "<<endl;
    }
    if(fHCorrectTPCQ3xEtaNeg && fHCorrectTPCQ3yEtaNeg){
      icentbin  = fHCorrectTPCQ3xEtaNeg->FindBin(fCent);
      qxEtaNeg -= fHCorrectTPCQ3xEtaNeg->GetBinContent(icentbin);
      qyEtaNeg -= fHCorrectTPCQ3yEtaNeg->GetBinContent(icentbin);
      //cout<<" =======> Info::  Applied recenter for TPC Eta Neg Psi3 "<<endl;
    }
  }
  else if(gPsiN==4){  // to be implemented!
    qxEtaPos -= 0.;
    qyEtaPos -= 0.;
    qxEtaNeg -= 0.;
    qyEtaNeg -= 0.;    
  }
  else{
    qxEtaPos -= 0.;
    qyEtaPos -= 0.;
    qxEtaNeg -= 0.;
    qyEtaNeg -= 0.;    
  }
  
  return;
}

void  AliAnalysisTaskGammaDeltaPID::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

  Int_t icentbin = 0;
  Double_t avgqx=0,avgqy=0; 
  //cout<<" => Before qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<"\tqnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  if(gPsiN==3){  ///<Q> correction for Psi3:  
    if(fHCorrectQ3xV0C && fHCorrectQ3yV0C){
      icentbin = fHCorrectQ3xV0C->FindBin(fCent);
      avgqx = fHCorrectQ3xV0C->GetBinContent(icentbin);
      avgqy = fHCorrectQ3yV0C->GetBinContent(icentbin);
      qnxV0C -= avgqx;
      qnyV0C -= avgqy;	
      //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    if(fHCorrectQ3xV0A && fHCorrectQ3yV0A){
      icentbin = fHCorrectQ3xV0A->FindBin(fCent);
      avgqx = fHCorrectQ3xV0A->GetBinContent(icentbin);
      avgqy = fHCorrectQ3yV0A->GetBinContent(icentbin);
      qnxV0A -= avgqx;
      qnyV0A -= avgqy;
      //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A "<<qnxV0A<<"\tqnyV0A"<<qnyV0A<<endl;
  }
  else{ /// Proper File which contain <q> for harmonic 'N' should be set in AddTask!! 
    if(fHCorrectQNxV0C && fHCorrectQNyV0C){
      icentbin = fHCorrectQNxV0C->FindBin(fCent);
      avgqx = fHCorrectQNxV0C->GetBinContent(icentbin);
      avgqy = fHCorrectQNyV0C->GetBinContent(icentbin);
      qnxV0C -= avgqx;
      qnyV0C -= avgqy;      
      //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    if(fHCorrectQNxV0A && fHCorrectQNyV0A){
      icentbin = fHCorrectQNxV0A->FindBin(fCent);
      avgqx = fHCorrectQNxV0A->GetBinContent(icentbin);
      avgqy = fHCorrectQNyV0A->GetBinContent(icentbin);
      qnxV0A -= avgqx;
      qnyV0A -= avgqy;           
      //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  }

 
  return;

}





Double_t AliAnalysisTaskGammaDeltaPID::GetMCEfficiencyWeightForTrack(Double_t fPt,Int_t gChrg,Int_t kPID){
  Int_t    iBinforNUE = 1;
  Double_t ptWgtMCNUE = 1.0;
	  
  if(gChrg>0){
    if(kPID==0){
      if(fHCorrectMCPosChrg){
	iBinforNUE = fHCorrectMCPosChrg->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCPosChrg->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==1){
      if(fHCorrectMCPosPion){
	iBinforNUE = fHCorrectMCPosPion->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCPosPion->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==2){
      if(fHCorrectMCPosKaon){
	iBinforNUE = fHCorrectMCPosKaon->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCPosKaon->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==3){
      if(fHCorrectMCPosProt){
	iBinforNUE = fHCorrectMCPosProt->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCPosProt->GetBinContent(iBinforNUE);
      }
    }    
  }/// +ve charge
  else{
    if(kPID==0){
      if(fHCorrectMCNegChrg){
	iBinforNUE = fHCorrectMCNegChrg->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCNegChrg->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==1){
      if(fHCorrectMCNegPion){
	iBinforNUE = fHCorrectMCNegPion->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCNegPion->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==2){
      if(fHCorrectMCNegKaon){
	iBinforNUE = fHCorrectMCNegKaon->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCNegKaon->GetBinContent(iBinforNUE);
      }
    }
    else if(kPID==3){
      if(fHCorrectMCNegProt){
	iBinforNUE = fHCorrectMCNegProt->FindBin(fPt);  
	ptWgtMCNUE = fHCorrectMCNegProt->GetBinContent(iBinforNUE);
      }
    }
  }/// -ve charge

  if(ptWgtMCNUE==0 || ptWgtMCNUE>100){ /// that means something went wrong or empty bin?
    ptWgtMCNUE = 1.0;  
  }
  return ptWgtMCNUE;
}//function ends
	


Bool_t AliAnalysisTaskGammaDeltaPID::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

  const AliAODVZERO *fAODV0 = (AliAODVZERO *) faod->GetVZEROData();
  Float_t fMultV0 = 0.;
  Float_t fPhiV0  = 0.;
  Float_t fV0chGain = 1.0;

  Double_t fQxV0CHarmN=0,fQyV0CHarmN=0,fQxV0AHarmN=0,fQyV0AHarmN=0;

  Double_t fSumMV0A = 0;
  Double_t fSumMV0C = 0;
  Int_t ibinV0=0;

  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

    fMultV0 = fAODV0->GetMultiplicity(iV0);

    /// V0 Channel Gain Correction:
    if(fHCorrectV0ChWeghts){ 
      ibinV0    = fHCorrectV0ChWeghts->FindBin(fVtxZ,iV0);
      fV0chGain = fHCorrectV0ChWeghts->GetBinContent(ibinV0); 
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    
    hAvgV0ChannelsvsVz->Fill(iV0+0.5, fVtxZ, fMultV0);

    fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);

    if(iV0 < 32){
      qnxV0C   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0C   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0C += fMultV0;
    }
    else if(iV0 >= 32){
      qnxV0A   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0A   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0A += fMultV0;
    } 
  }///V0 Channel loop

  /// Now the q vectors:
  if(fSumMV0A<=1e-4 || fSumMV0C<=1e-4){
    qnxV0C = 0;
    qnyV0C = 0;
    qnxV0A = 0;
    qnyV0A = 0;    
    return kFALSE;       
  }
  else{
    qnxV0C = qnxV0C/fSumMV0C;
    qnyV0C = qnyV0C/fSumMV0C;
    qnxV0A = qnxV0A/fSumMV0A;
    qnyV0A = qnyV0A/fSumMV0A;    
    return kTRUE;  
  }
  
}


void AliAnalysisTaskGammaDeltaPID::Terminate(Option_t *)  {
  //fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  //if (!fOutputList) return;
}





