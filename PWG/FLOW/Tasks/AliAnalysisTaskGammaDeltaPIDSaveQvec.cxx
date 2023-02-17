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
/* $Id: AliAnalysisTaskGammaDeltaPIDSaveQvec.cxx ver: 2.0                     $   */
/* Simple Task to fill V0 and ZDC Energies for Gain Calibration           */
/* Works with 15o and 18q/r. Support to be added for LHC10h               */
/* Developer: Md Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)        */
/* Last Modified: Aug 23, 2021,  First version committed                  */
/* Last Modified: Oct 08, 2021,  Second version committed                 */
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskGammaDeltaPIDSaveQvec.h"
#include "AliAnalysisTaskGammaDeltaPIDSaveQvecEvent.h"
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

ClassImp(AliAnalysisTaskGammaDeltaPIDSaveQvec)

AliAnalysisTaskGammaDeltaPIDSaveQvec::AliAnalysisTaskGammaDeltaPIDSaveQvec(const char *name):
  AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fEventList(NULL),
  fTempList(NULL),
  treeEvent(NULL),
  fpQvecEvent(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
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
  fHCorrectNUAkPIDPosPion(NULL),
  fHCorrectNUAkPIDNegPion(NULL),
  fHCorrectNUAkPIDPosKaon(NULL),
  fHCorrectNUAkPIDNegKaon(NULL),
  fHCorrectNUAkPIDPosProton(NULL),
  fHCorrectNUAkPIDNegProton(NULL),
  
  hAvgV2TPCvsCent(NULL),
  hAvgV2TPCvsCentPion(NULL),
  hAvgV2TPCvsCentKaon(NULL),
  hAvgV2TPCvsCentProton(NULL),
  
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
  
  hAvg3pC112vsCentOSPionPion(NULL),
  hAvg3pC123vsCentOSPionPion(NULL),
  hAvgDelta1vsCentOSPionPion(NULL),
  hAvgDelta2vsCentOSPionPion(NULL),
  hAvgDelta3vsCentOSPionPion(NULL),
  hAvgDelta4vsCentOSPionPion(NULL),
  hAvg3pC112vsCentOSKaonKaon(NULL),
  hAvg3pC123vsCentOSKaonKaon(NULL),
  hAvgDelta1vsCentOSKaonKaon(NULL),
  hAvgDelta2vsCentOSKaonKaon(NULL),
  hAvgDelta3vsCentOSKaonKaon(NULL),
  hAvgDelta4vsCentOSKaonKaon(NULL),
  hAvg3pC112vsCentOSProtonProton(NULL),
  hAvg3pC123vsCentOSProtonProton(NULL),
  hAvgDelta1vsCentOSProtonProton(NULL),
  hAvgDelta2vsCentOSProtonProton(NULL),
  hAvgDelta3vsCentOSProtonProton(NULL),
  hAvgDelta4vsCentOSProtonProton(NULL),

  hAvg3pC112vsCentOSPionCharge(NULL),
  hAvg3pC123vsCentOSPionCharge(NULL),
  hAvgDelta1vsCentOSPionCharge(NULL),
  hAvgDelta2vsCentOSPionCharge(NULL),
  hAvgDelta3vsCentOSPionCharge(NULL),
  hAvgDelta4vsCentOSPionCharge(NULL),		
  hAvg3pC112vsCentOSKaonCharge(NULL),
  hAvg3pC123vsCentOSKaonCharge(NULL),
  hAvgDelta1vsCentOSKaonCharge(NULL),
  hAvgDelta2vsCentOSKaonCharge(NULL),
  hAvgDelta3vsCentOSKaonCharge(NULL),
  hAvgDelta4vsCentOSKaonCharge(NULL),
  hAvg3pC112vsCentOSProtonCharge(NULL),
  hAvg3pC123vsCentOSProtonCharge(NULL),
  hAvgDelta1vsCentOSProtonCharge(NULL),
  hAvgDelta2vsCentOSProtonCharge(NULL),
  hAvgDelta3vsCentOSProtonCharge(NULL),
  hAvgDelta4vsCentOSProtonCharge(NULL),	
  

  hAvg3pC112vsCentPPPionPion(NULL),
  hAvg3pC123vsCentPPPionPion(NULL),
  hAvgDelta1vsCentPPPionPion(NULL),
  hAvgDelta2vsCentPPPionPion(NULL),
  hAvgDelta3vsCentPPPionPion(NULL),
  hAvgDelta4vsCentPPPionPion(NULL),
  hAvg3pC112vsCentPPKaonKaon(NULL),
  hAvg3pC123vsCentPPKaonKaon(NULL),
  hAvgDelta1vsCentPPKaonKaon(NULL),
  hAvgDelta2vsCentPPKaonKaon(NULL),
  hAvgDelta3vsCentPPKaonKaon(NULL),
  hAvgDelta4vsCentPPKaonKaon(NULL),
  hAvg3pC112vsCentPPProtonProton(NULL),
  hAvg3pC123vsCentPPProtonProton(NULL),
  hAvgDelta1vsCentPPProtonProton(NULL),
  hAvgDelta2vsCentPPProtonProton(NULL),
  hAvgDelta3vsCentPPProtonProton(NULL),
  hAvgDelta4vsCentPPProtonProton(NULL),

  hAvg3pC112vsCentPPPionCharge(NULL),
  hAvg3pC123vsCentPPPionCharge(NULL),
  hAvgDelta1vsCentPPPionCharge(NULL),
  hAvgDelta2vsCentPPPionCharge(NULL),
  hAvgDelta3vsCentPPPionCharge(NULL),
  hAvgDelta4vsCentPPPionCharge(NULL),		
  hAvg3pC112vsCentPPKaonCharge(NULL),
  hAvg3pC123vsCentPPKaonCharge(NULL),
  hAvgDelta1vsCentPPKaonCharge(NULL),
  hAvgDelta2vsCentPPKaonCharge(NULL),
  hAvgDelta3vsCentPPKaonCharge(NULL),
  hAvgDelta4vsCentPPKaonCharge(NULL),
  hAvg3pC112vsCentPPProtonCharge(NULL),
  hAvg3pC123vsCentPPProtonCharge(NULL),
  hAvgDelta1vsCentPPProtonCharge(NULL),
  hAvgDelta2vsCentPPProtonCharge(NULL),
  hAvgDelta3vsCentPPProtonCharge(NULL),
  hAvgDelta4vsCentPPProtonCharge(NULL),	

  
  hAvg3pC112vsCentNNPionPion(NULL),
  hAvg3pC123vsCentNNPionPion(NULL),
  hAvgDelta1vsCentNNPionPion(NULL),
  hAvgDelta2vsCentNNPionPion(NULL),
  hAvgDelta3vsCentNNPionPion(NULL),
  hAvgDelta4vsCentNNPionPion(NULL),
  hAvg3pC112vsCentNNKaonKaon(NULL),
  hAvg3pC123vsCentNNKaonKaon(NULL),
  hAvgDelta1vsCentNNKaonKaon(NULL),
  hAvgDelta2vsCentNNKaonKaon(NULL),
  hAvgDelta3vsCentNNKaonKaon(NULL),
  hAvgDelta4vsCentNNKaonKaon(NULL),
  hAvg3pC112vsCentNNProtonProton(NULL),
  hAvg3pC123vsCentNNProtonProton(NULL),
  hAvgDelta1vsCentNNProtonProton(NULL),
  hAvgDelta2vsCentNNProtonProton(NULL),
  hAvgDelta3vsCentNNProtonProton(NULL),
  hAvgDelta4vsCentNNProtonProton(NULL),

  hAvg3pC112vsCentNNPionCharge(NULL),
  hAvg3pC123vsCentNNPionCharge(NULL),
  hAvgDelta1vsCentNNPionCharge(NULL),
  hAvgDelta2vsCentNNPionCharge(NULL),
  hAvgDelta3vsCentNNPionCharge(NULL),
  hAvgDelta4vsCentNNPionCharge(NULL),		
  hAvg3pC112vsCentNNKaonCharge(NULL),
  hAvg3pC123vsCentNNKaonCharge(NULL),
  hAvgDelta1vsCentNNKaonCharge(NULL),
  hAvgDelta2vsCentNNKaonCharge(NULL),
  hAvgDelta3vsCentNNKaonCharge(NULL),
  hAvgDelta4vsCentNNKaonCharge(NULL),
  hAvg3pC112vsCentNNProtonCharge(NULL),
  hAvg3pC123vsCentNNProtonCharge(NULL),
  hAvgDelta1vsCentNNProtonCharge(NULL),
  hAvgDelta2vsCentNNProtonCharge(NULL),
  hAvgDelta3vsCentNNProtonCharge(NULL),
  hAvgDelta4vsCentNNProtonCharge(NULL),	
  
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
  hAvgV0ChannelsvsVz(NULL)
{
  
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQRe[c][h] = NULL;
      fCMEQIm[c][h] = NULL;
      fCMEMult[c][h] = NULL;
      fCMEQRePion[c][h] = NULL;
      fCMEQImPion[c][h] = NULL;
      fCMEMultPion[c][h] = NULL;
      fCMEQReKaon[c][h] = NULL;
      fCMEQImKaon[c][h] = NULL;
      fCMEMultKaon[c][h] = NULL;
      fCMEQReProton[c][h] = NULL;
      fCMEQImProton[c][h] = NULL;
      fCMEMultProton[c][h] = NULL;
    }
  }
  //@shi reset CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQReBothCharge[c][h] = NULL;
      fCMEQImBothCharge[c][h] = NULL;
      fCMEMultBothCharge[c][h] = NULL;
      fCMEQRePionBothCharge[c][h] = NULL;
      fCMEQImPionBothCharge[c][h] = NULL;
      fCMEMultPionBothCharge[c][h] = NULL;
      fCMEQReKaonBothCharge[c][h] = NULL;
      fCMEQImKaonBothCharge[c][h] = NULL;
      fCMEMultKaonBothCharge[c][h] = NULL;
      fCMEQReProtonBothCharge[c][h] = NULL;
      fCMEQImProtonBothCharge[c][h] = NULL;
      fCMEMultProtonBothCharge[c][h] = NULL;
    }
  }

  for(Int_t c=0;c<2;c++) {
	fCMEQ2Re4[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2[c] = NULL; // w^3*sin(2phi)
    fCMEw0[c] = NULL;    // w^0
    fCMEw1[c] = NULL;    // w^1
    fCMEw2[c] = NULL;    // w^2
    fCMEw3[c] = NULL;    // w^3
    fCMEw4[c] = NULL;    // w^4
  
	fCMEQ2Re4Pion[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Pion[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Pion[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Pion[c] = NULL; // w^3*sin(2phi)
    fCMEw0Pion[c] = NULL;    // w^0
    fCMEw1Pion[c] = NULL;    // w^1
    fCMEw2Pion[c] = NULL;    // w^2
    fCMEw3Pion[c] = NULL;    // w^3
    fCMEw4Pion[c] = NULL;    // w^4
    
    fCMEQ2Re4Kaon[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Kaon[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Kaon[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Kaon[c] = NULL; // w^3*sin(2phi)
    fCMEw0Kaon[c] = NULL;    // w^0
    fCMEw1Kaon[c] = NULL;    // w^1
    fCMEw2Kaon[c] = NULL;    // w^2
    fCMEw3Kaon[c] = NULL;    // w^3
    fCMEw4Kaon[c] = NULL;    // w^4
    
    fCMEQ2Re4Proton[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Proton[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Proton[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Proton[c] = NULL; // w^3*sin(2phi)
    fCMEw0Proton[c] = NULL;    // w^0
    fCMEw1Proton[c] = NULL;    // w^1
    fCMEw2Proton[c] = NULL;    // w^2
    fCMEw3Proton[c] = NULL;    // w^3
    fCMEw4Proton[c] = NULL;    // w^4
  }
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaDeltaPIDSaveQvec::AliAnalysisTaskGammaDeltaPIDSaveQvec():
  AliAnalysisTaskSE(),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fEventList(NULL),
  fTempList(NULL),
  treeEvent(NULL),
  fpQvecEvent(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),   
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
  fHCorrectNUAkPIDPosPion(NULL),
  fHCorrectNUAkPIDNegPion(NULL),
  fHCorrectNUAkPIDPosKaon(NULL),
  fHCorrectNUAkPIDNegKaon(NULL),
  fHCorrectNUAkPIDPosProton(NULL),
  fHCorrectNUAkPIDNegProton(NULL),
  
  hAvgV2TPCvsCent(NULL),
  hAvgV2TPCvsCentPion(NULL),
  hAvgV2TPCvsCentKaon(NULL),
  hAvgV2TPCvsCentProton(NULL),
  
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
  
  hAvg3pC112vsCentOSPionPion(NULL),
  hAvg3pC123vsCentOSPionPion(NULL),
  hAvgDelta1vsCentOSPionPion(NULL),
  hAvgDelta2vsCentOSPionPion(NULL),
  hAvgDelta3vsCentOSPionPion(NULL),
  hAvgDelta4vsCentOSPionPion(NULL),
  hAvg3pC112vsCentOSKaonKaon(NULL),
  hAvg3pC123vsCentOSKaonKaon(NULL),
  hAvgDelta1vsCentOSKaonKaon(NULL),
  hAvgDelta2vsCentOSKaonKaon(NULL),
  hAvgDelta3vsCentOSKaonKaon(NULL),
  hAvgDelta4vsCentOSKaonKaon(NULL),
  hAvg3pC112vsCentOSProtonProton(NULL),
  hAvg3pC123vsCentOSProtonProton(NULL),
  hAvgDelta1vsCentOSProtonProton(NULL),
  hAvgDelta2vsCentOSProtonProton(NULL),
  hAvgDelta3vsCentOSProtonProton(NULL),
  hAvgDelta4vsCentOSProtonProton(NULL),

  hAvg3pC112vsCentOSPionCharge(NULL),
  hAvg3pC123vsCentOSPionCharge(NULL),
  hAvgDelta1vsCentOSPionCharge(NULL),
  hAvgDelta2vsCentOSPionCharge(NULL),
  hAvgDelta3vsCentOSPionCharge(NULL),
  hAvgDelta4vsCentOSPionCharge(NULL),		
  hAvg3pC112vsCentOSKaonCharge(NULL),
  hAvg3pC123vsCentOSKaonCharge(NULL),
  hAvgDelta1vsCentOSKaonCharge(NULL),
  hAvgDelta2vsCentOSKaonCharge(NULL),
  hAvgDelta3vsCentOSKaonCharge(NULL),
  hAvgDelta4vsCentOSKaonCharge(NULL),
  hAvg3pC112vsCentOSProtonCharge(NULL),
  hAvg3pC123vsCentOSProtonCharge(NULL),
  hAvgDelta1vsCentOSProtonCharge(NULL),
  hAvgDelta2vsCentOSProtonCharge(NULL),
  hAvgDelta3vsCentOSProtonCharge(NULL),
  hAvgDelta4vsCentOSProtonCharge(NULL),	
  

  hAvg3pC112vsCentPPPionPion(NULL),
  hAvg3pC123vsCentPPPionPion(NULL),
  hAvgDelta1vsCentPPPionPion(NULL),
  hAvgDelta2vsCentPPPionPion(NULL),
  hAvgDelta3vsCentPPPionPion(NULL),
  hAvgDelta4vsCentPPPionPion(NULL),
  hAvg3pC112vsCentPPKaonKaon(NULL),
  hAvg3pC123vsCentPPKaonKaon(NULL),
  hAvgDelta1vsCentPPKaonKaon(NULL),
  hAvgDelta2vsCentPPKaonKaon(NULL),
  hAvgDelta3vsCentPPKaonKaon(NULL),
  hAvgDelta4vsCentPPKaonKaon(NULL),
  hAvg3pC112vsCentPPProtonProton(NULL),
  hAvg3pC123vsCentPPProtonProton(NULL),
  hAvgDelta1vsCentPPProtonProton(NULL),
  hAvgDelta2vsCentPPProtonProton(NULL),
  hAvgDelta3vsCentPPProtonProton(NULL),
  hAvgDelta4vsCentPPProtonProton(NULL),

  hAvg3pC112vsCentPPPionCharge(NULL),
  hAvg3pC123vsCentPPPionCharge(NULL),
  hAvgDelta1vsCentPPPionCharge(NULL),
  hAvgDelta2vsCentPPPionCharge(NULL),
  hAvgDelta3vsCentPPPionCharge(NULL),
  hAvgDelta4vsCentPPPionCharge(NULL),		
  hAvg3pC112vsCentPPKaonCharge(NULL),
  hAvg3pC123vsCentPPKaonCharge(NULL),
  hAvgDelta1vsCentPPKaonCharge(NULL),
  hAvgDelta2vsCentPPKaonCharge(NULL),
  hAvgDelta3vsCentPPKaonCharge(NULL),
  hAvgDelta4vsCentPPKaonCharge(NULL),
  hAvg3pC112vsCentPPProtonCharge(NULL),
  hAvg3pC123vsCentPPProtonCharge(NULL),
  hAvgDelta1vsCentPPProtonCharge(NULL),
  hAvgDelta2vsCentPPProtonCharge(NULL),
  hAvgDelta3vsCentPPProtonCharge(NULL),
  hAvgDelta4vsCentPPProtonCharge(NULL),	

  
  hAvg3pC112vsCentNNPionPion(NULL),
  hAvg3pC123vsCentNNPionPion(NULL),
  hAvgDelta1vsCentNNPionPion(NULL),
  hAvgDelta2vsCentNNPionPion(NULL),
  hAvgDelta3vsCentNNPionPion(NULL),
  hAvgDelta4vsCentNNPionPion(NULL),
  hAvg3pC112vsCentNNKaonKaon(NULL),
  hAvg3pC123vsCentNNKaonKaon(NULL),
  hAvgDelta1vsCentNNKaonKaon(NULL),
  hAvgDelta2vsCentNNKaonKaon(NULL),
  hAvgDelta3vsCentNNKaonKaon(NULL),
  hAvgDelta4vsCentNNKaonKaon(NULL),
  hAvg3pC112vsCentNNProtonProton(NULL),
  hAvg3pC123vsCentNNProtonProton(NULL),
  hAvgDelta1vsCentNNProtonProton(NULL),
  hAvgDelta2vsCentNNProtonProton(NULL),
  hAvgDelta3vsCentNNProtonProton(NULL),
  hAvgDelta4vsCentNNProtonProton(NULL),

  hAvg3pC112vsCentNNPionCharge(NULL),
  hAvg3pC123vsCentNNPionCharge(NULL),
  hAvgDelta1vsCentNNPionCharge(NULL),
  hAvgDelta2vsCentNNPionCharge(NULL),
  hAvgDelta3vsCentNNPionCharge(NULL),
  hAvgDelta4vsCentNNPionCharge(NULL),		
  hAvg3pC112vsCentNNKaonCharge(NULL),
  hAvg3pC123vsCentNNKaonCharge(NULL),
  hAvgDelta1vsCentNNKaonCharge(NULL),
  hAvgDelta2vsCentNNKaonCharge(NULL),
  hAvgDelta3vsCentNNKaonCharge(NULL),
  hAvgDelta4vsCentNNKaonCharge(NULL),
  hAvg3pC112vsCentNNProtonCharge(NULL),
  hAvg3pC123vsCentNNProtonCharge(NULL),
  hAvgDelta1vsCentNNProtonCharge(NULL),
  hAvgDelta2vsCentNNProtonCharge(NULL),
  hAvgDelta3vsCentNNProtonCharge(NULL),
  hAvgDelta4vsCentNNProtonCharge(NULL),	
  
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
  hAvgV0ChannelsvsVz(NULL)
{
  
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQRe[c][h] = NULL;
      fCMEQIm[c][h] = NULL;
      fCMEMult[c][h] = NULL;
      fCMEQRePion[c][h] = NULL;
      fCMEQImPion[c][h] = NULL;
      fCMEMultPion[c][h] = NULL;
      fCMEQReKaon[c][h] = NULL;
      fCMEQImKaon[c][h] = NULL;
      fCMEMultKaon[c][h] = NULL;
      fCMEQReProton[c][h] = NULL;
      fCMEQImProton[c][h] = NULL;
      fCMEMultProton[c][h] = NULL;
    }
  }
  //@shi reset CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQReBothCharge[c][h] = NULL;
      fCMEQImBothCharge[c][h] = NULL;
      fCMEMultBothCharge[c][h] = NULL;
      fCMEQRePionBothCharge[c][h] = NULL;
      fCMEQImPionBothCharge[c][h] = NULL;
      fCMEMultPionBothCharge[c][h] = NULL;
      fCMEQReKaonBothCharge[c][h] = NULL;
      fCMEQImKaonBothCharge[c][h] = NULL;
      fCMEMultKaonBothCharge[c][h] = NULL;
      fCMEQReProtonBothCharge[c][h] = NULL;
      fCMEQImProtonBothCharge[c][h] = NULL;
      fCMEMultProtonBothCharge[c][h] = NULL;
    }
  }
  
  for(Int_t c=0;c<2;c++) {
	fCMEQ2Re4[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2[c] = NULL; // w^3*sin(2phi)
    fCMEw0[c] = NULL;    // w^0
    fCMEw1[c] = NULL;    // w^1
    fCMEw2[c] = NULL;    // w^2
    fCMEw3[c] = NULL;    // w^3
    fCMEw4[c] = NULL;    // w^4
  
	fCMEQ2Re4Pion[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Pion[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Pion[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Pion[c] = NULL; // w^3*sin(2phi)
    fCMEw0Pion[c] = NULL;    // w^0
    fCMEw1Pion[c] = NULL;    // w^1
    fCMEw2Pion[c] = NULL;    // w^2
    fCMEw3Pion[c] = NULL;    // w^3
    fCMEw4Pion[c] = NULL;    // w^4
    
    fCMEQ2Re4Kaon[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Kaon[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Kaon[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Kaon[c] = NULL; // w^3*sin(2phi)
    fCMEw0Kaon[c] = NULL;    // w^0
    fCMEw1Kaon[c] = NULL;    // w^1
    fCMEw2Kaon[c] = NULL;    // w^2
    fCMEw3Kaon[c] = NULL;    // w^3
    fCMEw4Kaon[c] = NULL;    // w^4
    
    fCMEQ2Re4Proton[c] = NULL; // w^2*cos(4phi)
    fCMEQ3Re2Proton[c] = NULL; // w^3*cos(2phi)
    fCMEQ2Im4Proton[c] = NULL; // w^2*sin(4phi)
    fCMEQ3Im2Proton[c] = NULL; // w^3*sin(2phi)
    fCMEw0Proton[c] = NULL;    // w^0
    fCMEw1Proton[c] = NULL;    // w^1
    fCMEw2Proton[c] = NULL;    // w^2
    fCMEw3Proton[c] = NULL;    // w^3
    fCMEw4Proton[c] = NULL;    // w^4
  }
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskGammaDeltaPIDSaveQvec::~AliAnalysisTaskGammaDeltaPIDSaveQvec()
{
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!
  if(fTempList)      delete fTempList;
  if(fListHist)      delete fListHist;  
  if(fListTRKCorr)   delete fListTRKCorr;
  if(fListNUACorr)   delete fListNUACorr;
  if(fListV0MCorr)   delete fListV0MCorr;

  if(fV0CutPU)      delete fV0CutPU;
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU;   
}










//________________ Define Histograms _______________
void AliAnalysisTaskGammaDeltaPIDSaveQvec::UserCreateOutputObjects()
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

  // base list to hold all temp objects:
  fTempList = new TList();
  fTempList->SetName("temp");
  fTempList->SetOwner(kTRUE);
  
  SetupQAHistograms();
  SetupAnalysisHistograms();
  SetupPileUpRemovalFunctions();
  SetupEventAndTaskConfigInfo();

  SetupQvecSavingObjects();

 
    
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
void AliAnalysisTaskGammaDeltaPIDSaveQvec::UserExec(Option_t*) {
 
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
  
  UShort_t bunchCrossNumber = fAOD->GetBunchCrossNumber();
  UInt_t timeStamp = fAOD->GetTimeStamp();
  
  fpQvecEvent->setRawPeriod(period);
  fpQvecEvent->setRawOrbitNumber24(orbit24);
  fpQvecEvent->setOrbitNumber(orbit);
  fpQvecEvent->setBunchCrossNumber(bunchCrossNumber);
  fpQvecEvent->setTimeStamp(timeStamp);
  
  Bool_t kPileupEvent = kFALSE;

  //Double_t fSumQnxNeg = 0, fSumQnyNeg = 0, fSumQnxPos = 0, fSumQnyPos = 0;
  Double_t fSumQnxNeg[3] = {0,}; // Array: Q2x, Q3x, Q4x,...
  Double_t fSumQnyNeg[3] = {0,}; // Array: Q2y, Q3y, Q4y,... 
  Double_t fSumQnxPos[3] = {0,};
  Double_t fSumQnyPos[3] = {0,};  
  
  Double_t fMultNeg = 0, fMultPos = 0;

  
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
  
  if(fQ2xPos != 0 && fQ2yPos != 0){
    fPsiNTPCPos = (1./2)*TMath::ATan2(fQ2yPos,fQ2xPos);
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
     
  Double_t fQnxV0C=0, fQnyV0C=0, fSumMV0C = 0, fQnxV0A=0, fQnyV0A=0, fSumMV0A = 0; 
  Bool_t kPassV0 = GetGainCorrectedV0Qvector(fAOD, fVertexZEvent, 2, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A, fSumMV0C, fSumMV0A);

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
    
  fPsiNV0C = (1./gHarmonic)*TMath::ATan2(fQnyV0C,fQnxV0C);
  if(fPsiNV0C < 0) fPsiNV0C += TMath::TwoPi()/gHarmonic;
  fPsiNV0A = (1./gHarmonic)*TMath::ATan2(fQnyV0A,fQnxV0A);
  if(fPsiNV0A < 0) fPsiNV0A += TMath::TwoPi()/gHarmonic;
  fHistV0CPsiNEventPlane->Fill(centrality,fPsiNV0C);
  fHistV0APsiNEventPlane->Fill(centrality,fPsiNV0A);

 
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



  

  //If we skip the nested loops and main Analysis:
  if(bSkipNestedLoop){
 
    fHistVertexZcm->Fill(pVtxZ);
    fCentDistAfterCut->Fill(centrality);
    
    PostData(1,fListHist);
    return;     //Just fill QAs and get out..
  }


  Double_t fSelectedV0PsiN = 0;
  Double_t fSelectedV0Psi3 = 0;
  

  if(sDetectorForEP.Contains("V0C")){
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
  
  //--------> Get V0A and V0C tower data
  const AliAODVZERO *fAODV0 = (AliAODVZERO *) fAOD->GetVZEROData();
  fpQvecEvent->setTowV0Craw0(fAODV0->GetMultiplicity(0));
  fpQvecEvent->setTowV0Craw1(fAODV0->GetMultiplicity(1));
  fpQvecEvent->setTowV0Craw2(fAODV0->GetMultiplicity(2));
  fpQvecEvent->setTowV0Craw3(fAODV0->GetMultiplicity(3));
  fpQvecEvent->setTowV0Craw4(fAODV0->GetMultiplicity(4));
  fpQvecEvent->setTowV0Craw5(fAODV0->GetMultiplicity(5));
  fpQvecEvent->setTowV0Craw6(fAODV0->GetMultiplicity(6));
  fpQvecEvent->setTowV0Craw7(fAODV0->GetMultiplicity(7));
  fpQvecEvent->setTowV0Craw8(fAODV0->GetMultiplicity(8));
  fpQvecEvent->setTowV0Craw9(fAODV0->GetMultiplicity(9));
  fpQvecEvent->setTowV0Craw10(fAODV0->GetMultiplicity(10));
  fpQvecEvent->setTowV0Craw11(fAODV0->GetMultiplicity(11));
  fpQvecEvent->setTowV0Craw12(fAODV0->GetMultiplicity(12));
  fpQvecEvent->setTowV0Craw13(fAODV0->GetMultiplicity(13));
  fpQvecEvent->setTowV0Craw14(fAODV0->GetMultiplicity(14));
  fpQvecEvent->setTowV0Craw15(fAODV0->GetMultiplicity(15));
  fpQvecEvent->setTowV0Craw16(fAODV0->GetMultiplicity(16));
  fpQvecEvent->setTowV0Craw17(fAODV0->GetMultiplicity(17));
  fpQvecEvent->setTowV0Craw18(fAODV0->GetMultiplicity(18));
  fpQvecEvent->setTowV0Craw19(fAODV0->GetMultiplicity(19));
  fpQvecEvent->setTowV0Craw20(fAODV0->GetMultiplicity(20));
  fpQvecEvent->setTowV0Craw21(fAODV0->GetMultiplicity(21));
  fpQvecEvent->setTowV0Craw22(fAODV0->GetMultiplicity(22));
  fpQvecEvent->setTowV0Craw23(fAODV0->GetMultiplicity(23));
  fpQvecEvent->setTowV0Craw24(fAODV0->GetMultiplicity(24));
  fpQvecEvent->setTowV0Craw25(fAODV0->GetMultiplicity(25));
  fpQvecEvent->setTowV0Craw26(fAODV0->GetMultiplicity(26));
  fpQvecEvent->setTowV0Craw27(fAODV0->GetMultiplicity(27));
  fpQvecEvent->setTowV0Craw28(fAODV0->GetMultiplicity(28));
  fpQvecEvent->setTowV0Craw29(fAODV0->GetMultiplicity(29));
  fpQvecEvent->setTowV0Craw30(fAODV0->GetMultiplicity(30));
  fpQvecEvent->setTowV0Craw31(fAODV0->GetMultiplicity(31));
  
  fpQvecEvent->setTowV0Araw0(fAODV0->GetMultiplicity(32));
  fpQvecEvent->setTowV0Araw1(fAODV0->GetMultiplicity(33));
  fpQvecEvent->setTowV0Araw2(fAODV0->GetMultiplicity(34));
  fpQvecEvent->setTowV0Araw3(fAODV0->GetMultiplicity(35));
  fpQvecEvent->setTowV0Araw4(fAODV0->GetMultiplicity(36));
  fpQvecEvent->setTowV0Araw5(fAODV0->GetMultiplicity(37));
  fpQvecEvent->setTowV0Araw6(fAODV0->GetMultiplicity(38));
  fpQvecEvent->setTowV0Araw7(fAODV0->GetMultiplicity(39));
  fpQvecEvent->setTowV0Araw8(fAODV0->GetMultiplicity(40));
  fpQvecEvent->setTowV0Araw9(fAODV0->GetMultiplicity(41));
  fpQvecEvent->setTowV0Araw10(fAODV0->GetMultiplicity(42));
  fpQvecEvent->setTowV0Araw11(fAODV0->GetMultiplicity(43));
  fpQvecEvent->setTowV0Araw12(fAODV0->GetMultiplicity(44));
  fpQvecEvent->setTowV0Araw13(fAODV0->GetMultiplicity(45));
  fpQvecEvent->setTowV0Araw14(fAODV0->GetMultiplicity(46));
  fpQvecEvent->setTowV0Araw15(fAODV0->GetMultiplicity(47));
  fpQvecEvent->setTowV0Araw16(fAODV0->GetMultiplicity(48));
  fpQvecEvent->setTowV0Araw17(fAODV0->GetMultiplicity(49));
  fpQvecEvent->setTowV0Araw18(fAODV0->GetMultiplicity(50));
  fpQvecEvent->setTowV0Araw19(fAODV0->GetMultiplicity(51));
  fpQvecEvent->setTowV0Araw20(fAODV0->GetMultiplicity(52));
  fpQvecEvent->setTowV0Araw21(fAODV0->GetMultiplicity(53));
  fpQvecEvent->setTowV0Araw22(fAODV0->GetMultiplicity(54));
  fpQvecEvent->setTowV0Araw23(fAODV0->GetMultiplicity(55));
  fpQvecEvent->setTowV0Araw24(fAODV0->GetMultiplicity(56));
  fpQvecEvent->setTowV0Araw25(fAODV0->GetMultiplicity(57));
  fpQvecEvent->setTowV0Araw26(fAODV0->GetMultiplicity(58));
  fpQvecEvent->setTowV0Araw27(fAODV0->GetMultiplicity(59));
  fpQvecEvent->setTowV0Araw28(fAODV0->GetMultiplicity(60));
  fpQvecEvent->setTowV0Araw29(fAODV0->GetMultiplicity(61));
  fpQvecEvent->setTowV0Araw30(fAODV0->GetMultiplicity(62));
  fpQvecEvent->setTowV0Araw31(fAODV0->GetMultiplicity(63));
  
  //=============== Get the ZDC data ==================== @Shi
 
  AliAODZDC *aodZDC = fAOD->GetZDCData();

  if(!aodZDC) {
    printf("\n ********* Error: could not find ZDC data ************ \n ");
  }
  else if(aodZDC){
    const Double_t *fZNATowerRawAOD = aodZDC->GetZNATowerEnergy();
    const Double_t *fZNCTowerRawAOD = aodZDC->GetZNCTowerEnergy();
    
    //const Double_t *fZPATowerRawAOD = aodZDC->GetZPATowerEnergy();
    //const Double_t *fZPCTowerRawAOD = aodZDC->GetZPCTowerEnergy();  
    //const Int_t nZDCChannel = 5;
    
	fpQvecEvent->setTowZNCraw0(fZNCTowerRawAOD[0]);
	fpQvecEvent->setTowZNCraw1(fZNCTowerRawAOD[1]);
	fpQvecEvent->setTowZNCraw2(fZNCTowerRawAOD[2]);
	fpQvecEvent->setTowZNCraw3(fZNCTowerRawAOD[3]);
	fpQvecEvent->setTowZNCraw4(fZNCTowerRawAOD[4]);

	fpQvecEvent->setTowZNAraw0(fZNATowerRawAOD[0]);
	fpQvecEvent->setTowZNAraw1(fZNATowerRawAOD[1]);
	fpQvecEvent->setTowZNAraw2(fZNATowerRawAOD[2]);
	fpQvecEvent->setTowZNAraw3(fZNATowerRawAOD[3]);
	fpQvecEvent->setTowZNAraw4(fZNATowerRawAOD[4]);
	

  }

  // ================================================================================> Set event info 
    
  fpQvecEvent->setRunNum(runNumber);
  fpQvecEvent->setCentrality(centrality);
  fpQvecEvent->setVtxPosX(pVtxX);
  fpQvecEvent->setVtxPosY(pVtxY);
  fpQvecEvent->setVtxPosZ(pVtxZ);
  
  fpQvecEvent->setVZCRe(fQnxV0C);
  fpQvecEvent->setVZCIm(fQnyV0C);
  fpQvecEvent->setVZCM(fSumMV0C);
  fpQvecEvent->setVZARe(fQnxV0A);
  fpQvecEvent->setVZAIm(fQnyV0A);
  fpQvecEvent->setVZAM(fSumMV0A);

  

  /// We Are About to Start Main Analysis Below: //
  
 
  Int_t kPIDtrk1=gParticleID;   /// gParticleID is Set From AddTask. Both Identified..
  Int_t kPIDtrk2=gParticleID;   /// 0 = hadron (h-h), 1 = Pi-Pi, 2 = K-K, 3 = Prot-Prot, 

  ///For single Identified cases:
  
  //if(gParticleID==10){
    //kPIDtrk1 = 0; //Ch
    //kPIDtrk2 = 1; //Pion
  //}
  //else if(gParticleID==20){
    //kPIDtrk1 = 0; //Ch
    //kPIDtrk2 = 2; //Kaon
  //}
  //else if(gParticleID==30){
    //kPIDtrk1 = 0; //Ch
    //kPIDtrk2 = 3; //P
  //}


///Track variables:
  Int_t   trk1Chrg=0, trk1TpcNC=0;
  Int_t   trk2Chrg=0, trk2TpcNC=0;
  //Bool_t  bPIDoktrk1=kFALSE, bPIDoktrk2=kFALSE;

  
  Double_t trk1Pt=0,trk1Phi=0,trk1Eta=0,trk1DCAxy=0.0, trk1DCAz=0.0,trk1Chi2=0,trk1dEdx=0,trk1Wgt=1.0;
  Double_t trk2Pt=0,trk2Phi=0,trk2Eta=0,trk2DCAxy=0.0, trk2DCAz=0.0,trk2Chi2=0,trk2dEdx=0,trk2Wgt=1.0;  
  Double_t wgtComb1Ch = 1.0;
  Double_t wgtComb2Ch = 1.0; 
  Double_t ptWgtMCChtrk1  = 1.0, WgtNUAChtrk1   = 1.0;
  Double_t ptWgtMCChtrk2  = 1.0, WgtNUAChtrk2   = 1.0;
  Double_t ptWgtMCPIDtrk1 = 1.0, WgtNUAPIDtrk1  = 1.0;
  Double_t ptWgtMCPIDtrk2 = 1.0, WgtNUAPIDtrk2  = 1.0;
  Double_t wgt1PIDparticle= 1.0, wgt2PIDparticle= 1.0; 

  Double_t localSumQ2x =0,localSumQ2y=0;
  Double_t localSumQ3x =0,localSumQ3y=0;
  Double_t localMultTPC=0; 


  Bool_t isTrk1Pion = kFALSE;
  Bool_t isTrk1Kaon = kFALSE;
  Bool_t isTrk1Proton = kFALSE;
  Double_t ptWgtMCPiontrk1 = 1.0;
  Double_t ptWgtMCKaontrk1 = 1.0;
  Double_t ptWgtMCProtontrk1 = 1.0;
  Bool_t isTrk2Pion = kFALSE;
  Bool_t isTrk2Kaon = kFALSE;
  Bool_t isTrk2Proton = kFALSE;
  Double_t ptWgtMCPionTrk2 = 1.0;
  Double_t ptWgtMCKaonTrk2 = 1.0;
  Double_t ptWgtMCProtonTrk2 = 1.0;
  Double_t wgtComb1PIDPion = 1.0;
  Double_t wgtComb1PIDKaon = 1.0;
  Double_t wgtComb1PIDProton = 1.0;
  Double_t wgtComb2PIDPion = 1.0;
  Double_t wgtComb2PIDKaon = 1.0;
  Double_t wgtComb2PIDProton = 1.0;
  Int_t chargeIndex = 0;
  
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
     
      
      //Apply track cuts for TPC EP here:
      //if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && !((trk1Eta >= fEtaGapNeg) && (trk1Eta <= fEtaGapPos)) && (trk1dEdx >= fTPCdEdxMin) && (trk1TpcNC >= fTPCclustMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fTrkChi2Max) && TMath::Abs(trk1Chrg)) {
	  if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && (trk1dEdx >= fTPCdEdxMin) && (trk1TpcNC >= fTPCclustMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fTrkChi2Max) && TMath::Abs(trk1Chrg)) {

	WgtNUAChtrk1  = 1.0;   
	WgtNUAPIDtrk1 = 1.0;
	ptWgtMCChtrk1 = 1.0;  
	ptWgtMCPIDtrk1 = 1.0;


	WgtNUAChtrk1  = GetNUAWeightForTrack(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg);          
	ptWgtMCChtrk1 = GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,0);

	wgtComb1Ch  = WgtNUAChtrk1*ptWgtMCChtrk1;    /// Charge

	
	ptWgtMCPiontrk1 = 1.0;
	ptWgtMCKaontrk1 = 1.0;
	ptWgtMCProtontrk1 = 1.0;
	wgtComb1PIDPion = 1.0;
	wgtComb1PIDKaon = 1.0;
	wgtComb1PIDProton = 1.0;
	
	// pion 
	isTrk1Pion = kFALSE;
	isTrk1Pion = CheckPIDofParticle(AODtrack1,1); // 1 for pion
	if (isTrk1Pion) {
		ptWgtMCPiontrk1= GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,1);
		WgtNUAPIDtrk1 = GetNUAWeightForTrackPID(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg,1); 
		wgtComb1PIDPion = WgtNUAPIDtrk1*ptWgtMCPiontrk1;
	}
	
	// kaon
	isTrk1Kaon = kFALSE;
	isTrk1Kaon = CheckPIDofParticle(AODtrack1,2); // 2 for kaon
	if (isTrk1Kaon) {
		ptWgtMCKaontrk1= GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,2);
		WgtNUAPIDtrk1 = GetNUAWeightForTrackPID(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg,2); 
		wgtComb1PIDKaon = WgtNUAPIDtrk1*ptWgtMCKaontrk1;
	}
	
	// proton
	isTrk1Proton = kFALSE;
	isTrk1Proton = CheckPIDofParticle(AODtrack1,3); // 3 for proton
	if (isTrk1Proton) {
		ptWgtMCProtontrk1= GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,3);
		WgtNUAPIDtrk1 = GetNUAWeightForTrackPID(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg,3); 
		wgtComb1PIDProton = WgtNUAPIDtrk1*ptWgtMCProtontrk1;
	}
	
	
	/// Rihan: The part below is only relevant for CME Analysis Only (no Lambda):
	
	//bPIDoktrk1=kFALSE;
	//bPIDoktrk2=kFALSE;
	
	//bPIDoktrk1 = CheckPIDofParticle(AODtrack1,kPIDtrk1);   // check if track1 is of desired PID request #1,

	//if(!bPIDoktrk1)
	  //bPIDoktrk2 = CheckPIDofParticle(AODtrack1,kPIDtrk2);  // if not request #1, then check if track satisfies request #2

	//if(!bPIDoktrk1 && !bPIDoktrk2)    /// the track-1 is neither of the desired PIDs, then skip this track. 
	  //continue;

       
	Double_t fPsiNEvent = 0, fPsi3Event = 0;

	///Choose whether to use TPC or V0EP
	if(bUseV0EventPlane){
	  fPsiNEvent = fSelectedV0PsiN;
	  fPsi3Event = fSelectedV0Psi3;
	}
	else{ 
	  if(trk1Eta >= 0){
	    fPsiNEvent = fPsiNTPCNeg;
	    fPsi3Event = fPsi3TPCNeg;	    
	  }
	  else{
	    fPsiNEvent = fPsiNTPCPos;
	    fPsi3Event = fPsi3TPCPos;
	  }
	}

	// save V2
	if(trk1Eta >= 0){
	  localSumQ2x = fSumQnxNeg[0]; /// We need the full Q-sum. Then remove only current track-2
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
	
	
	if(localMultTPC>0){
	  localSumQ2x = localSumQ2x/localMultTPC;
	  localSumQ2y = localSumQ2y/localMultTPC;
	
	
	  fPsiNTPCPos = (1./2)*TMath::ATan2(localSumQ2y,localSumQ2x);
	  if(fPsiNTPCPos < 0) fPsiNTPCPos += TMath::TwoPi()/2;	

	  // force using TPC to save the results
	  fPsiNEvent = fPsiNTPCPos;
	
	  hAvgV2TPCvsCent->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNEvent),wgtComb1Ch);
	  if (isTrk1Pion) {
		hAvgV2TPCvsCentPion->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNEvent),wgtComb1PIDPion);
	  } else if (isTrk1Kaon) {
		hAvgV2TPCvsCentKaon->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNEvent),wgtComb1PIDKaon);
	  } else if (isTrk1Proton) {
		hAvgV2TPCvsCentProton->Fill(centrality,TMath::Cos(2*trk1Phi - 2*fPsiNEvent),wgtComb1PIDProton);
	  }
	}
	
	// save v2 for each plane
	Double_t fPsi2TPC = 0;
	
	fPsi2TPC = (1./2)*TMath::ATan2(localSumQ2y,localSumQ2x);
	if(fPsi2TPC < 0) fPsi2TPC += TMath::TwoPi()/2;
	
	

	


	// ================================ save Qvec ===================================

	chargeIndex = (trk1Chrg > 0. ? 0 : 1);
	
	for (Int_t h=0;h<fCRCnHar;h++) { // h = 0, 1
		// all charged (Particle is UN-Identified)
		fCMEQRe[chargeIndex][h]->Fill(trk1Eta,wgtComb1Ch*TMath::Cos((h+1.)*trk1Phi)); // w*cos(phi) and w*cos(2phi)
		fCMEQIm[chargeIndex][h]->Fill(trk1Eta,wgtComb1Ch*TMath::Sin((h+1.)*trk1Phi));
		fCMEMult[chargeIndex][h]->Fill(trk1Eta,wgtComb1Ch);
		fCMEQRe[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.)*TMath::Cos((h+1.)*trk1Phi));
		fCMEQIm[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.)*TMath::Sin((h+1.)*trk1Phi));
		fCMEMult[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.));
		if (h == 0) {
			fCMEQ2Re4[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
			fCMEQ3Re2[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,3)*TMath::Cos(2*trk1Phi)); // w^3*cos(2phi)
			fCMEQ2Im4[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,2)*TMath::Sin(4*trk1Phi)); // w^2*sin(4phi)
			fCMEQ3Im2[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,3)*TMath::Sin(2*trk1Phi)); // w^3*sin(2phi)
			fCMEw0[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,0));
			fCMEw1[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,1));
			fCMEw2[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,2));
			fCMEw3[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,3));
			fCMEw4[chargeIndex]->Fill(trk1Eta,pow(wgtComb1Ch,4));
		}
		
		fCMEQReBothCharge[0][h]->Fill(trk1Eta,wgtComb1Ch*TMath::Cos((h+1.)*trk1Phi));
		fCMEQImBothCharge[0][h]->Fill(trk1Eta,wgtComb1Ch*TMath::Sin((h+1.)*trk1Phi));
		fCMEMultBothCharge[0][h]->Fill(trk1Eta,wgtComb1Ch);
		fCMEQReBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.)*TMath::Cos((h+1.)*trk1Phi));
		fCMEQImBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.)*TMath::Sin((h+1.)*trk1Phi));
		fCMEMultBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1Ch,2.));
		
		if (trk1Pt <= 2.0) { // beyond 2 GeV the purity drops
			if (isTrk1Pion) {
				// pion 
				fCMEQRePion[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDPion*TMath::Cos((h+1.)*trk1Phi)); // w*cos(phi) and w*cos(2phi)
				fCMEQImPion[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDPion*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultPion[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDPion);
				fCMEQRePion[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImPion[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultPion[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.));
				if (h == 0) {
					fCMEQ2Re4Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Re2Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,3)*TMath::Cos(2*trk1Phi)); // w^2*cos(4phi)
					fCMEQ2Im4Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,2)*TMath::Sin(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Im2Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,3)*TMath::Sin(2*trk1Phi)); // w^2*cos(4phi)
					fCMEw0Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,0));
					fCMEw1Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,1));
					fCMEw2Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,2));
					fCMEw3Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,3));
					fCMEw4Pion[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDPion,4));
				}
				
				fCMEQRePionBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDPion*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImPionBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDPion*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultPionBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDPion);
				fCMEQRePionBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImPionBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultPionBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDPion,2.));
			} else if (isTrk1Kaon) {
				// kaon
				fCMEQReKaon[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDKaon*TMath::Cos((h+1.)*trk1Phi)); // w*cos(phi) and w*cos(2phi)
				fCMEQImKaon[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDKaon*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultKaon[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDKaon);
				fCMEQReKaon[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImKaon[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultKaon[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.));
				if (h == 0) {
					fCMEQ2Re4Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Re2Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,3)*TMath::Cos(2*trk1Phi)); // w^2*cos(4phi)
					fCMEQ2Im4Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2)*TMath::Sin(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Im2Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,3)*TMath::Sin(2*trk1Phi)); // w^2*cos(4phi)
					fCMEw0Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,0));
					fCMEw1Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,1));
					fCMEw2Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2));
					fCMEw3Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,3));
					fCMEw4Kaon[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDKaon,4));
				}
				
				fCMEQReKaonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDKaon*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImKaonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDKaon*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultKaonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDKaon);
				fCMEQReKaonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImKaonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultKaonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDKaon,2.));
			} else if (isTrk1Proton) {
				// proton
				fCMEQReProton[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDProton*TMath::Cos((h+1.)*trk1Phi)); // w*cos(phi) and w*cos(2phi)
				fCMEQImProton[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDProton*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultProton[chargeIndex][h]->Fill(trk1Eta,wgtComb1PIDProton);
				fCMEQReProton[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImProton[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultProton[2+chargeIndex][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.));
				if (h == 0) {
					fCMEQ2Re4Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,2)*TMath::Cos(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Re2Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,3)*TMath::Cos(2*trk1Phi)); // w^2*cos(4phi)
					fCMEQ2Im4Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,2)*TMath::Sin(4*trk1Phi)); // w^2*cos(4phi)
					fCMEQ3Im2Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,3)*TMath::Sin(2*trk1Phi)); // w^2*cos(4phi)
					fCMEw0Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,0));
					fCMEw1Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,1));
					fCMEw2Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,2));
					fCMEw3Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,3));
					fCMEw4Proton[chargeIndex]->Fill(trk1Eta,pow(wgtComb1PIDProton,4));
				}
				
				fCMEQReProtonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDProton*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImProtonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDProton*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultProtonBothCharge[0][h]->Fill(trk1Eta,wgtComb1PIDProton);
				fCMEQReProtonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.)*TMath::Cos((h+1.)*trk1Phi));
				fCMEQImProtonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.)*TMath::Sin((h+1.)*trk1Phi));
				fCMEMultProtonBothCharge[1][h]->Fill(trk1Eta,pow(wgtComb1PIDProton,2.));
			}
		}
		
	}
	

	///---> 2nd track Loop   
	for(Int_t jTrack = 0; jTrack < iTrack; jTrack++) { 

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
	    //if((trk2Pt <= fMaxPtCut) && (trk2Pt >= fMinPtCut) && (trk2Eta <= fMaxEtaCut) && (trk2Eta >= fMinEtaCut) && !((trk2Eta >= fEtaGapNeg) && (trk2Eta <= fEtaGapPos)) && (trk2dEdx >= fTPCdEdxMin) && (trk2TpcNC >= fTPCclustMin) && (trk2Chi2 >= fTrkChi2Min) && (trk2Chi2 <= fTrkChi2Max) && TMath::Abs(trk2Chrg)) {
	    if((trk2Pt <= fMaxPtCut) && (trk2Pt >= fMinPtCut) && (trk2Eta <= fMaxEtaCut) && (trk2Eta >= fMinEtaCut) && (trk2dEdx >= fTPCdEdxMin) && (trk2TpcNC >= fTPCclustMin) && (trk2Chi2 >= fTrkChi2Min) && (trk2Chi2 <= fTrkChi2Max) && TMath::Abs(trk2Chrg)) {


		WgtNUAChtrk2  = 1.0;   
		WgtNUAPIDtrk2 = 1.0;
		ptWgtMCChtrk2 = 1.0;  
		ptWgtMCPIDtrk2 = 1.0;


		WgtNUAChtrk2  = GetNUAWeightForTrack(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg);    
		ptWgtMCChtrk2 = GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,0);
		
		wgtComb2Ch  = WgtNUAChtrk2*ptWgtMCChtrk2;      /// Combined weight for trk2 Ch
		
		ptWgtMCPionTrk2 = 1.0;
		ptWgtMCKaonTrk2 = 1.0;
		ptWgtMCProtonTrk2 = 1.0;
		wgtComb2PIDPion = 1.0;
		wgtComb2PIDKaon = 1.0;
		wgtComb2PIDProton = 1.0;
		
		
		
		// pion 
		isTrk2Pion = kFALSE;
		isTrk2Pion = CheckPIDofParticle(AODtrack2,1); // 1 for pion
		if (isTrk2Pion) {
			ptWgtMCPionTrk2= GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,1);
			WgtNUAPIDtrk2 = GetNUAWeightForTrackPID(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg,1);
			wgtComb2PIDPion = WgtNUAPIDtrk2*ptWgtMCPionTrk2;
		}
		
		// kaon
		isTrk2Kaon = kFALSE;
		isTrk2Kaon = CheckPIDofParticle(AODtrack2,2); // 2 for kaon
		if (isTrk2Kaon) {
			ptWgtMCKaonTrk2= GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,2);
			WgtNUAPIDtrk2 = GetNUAWeightForTrackPID(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg,2);
			wgtComb2PIDKaon = WgtNUAPIDtrk2*ptWgtMCKaonTrk2;
		}
		
		// proton
		isTrk2Proton = kFALSE;
		isTrk2Proton = CheckPIDofParticle(AODtrack2,3); // 2 for proton
		if (isTrk2Proton) {
			ptWgtMCProtonTrk2= GetMCEfficiencyWeightForTrack(trk2Pt,trk2Chrg,3);
			WgtNUAPIDtrk2 = GetNUAWeightForTrackPID(fVertexZEvent,trk2Phi,trk2Eta,trk2Chrg,3);
			wgtComb2PIDProton = WgtNUAPIDtrk2*ptWgtMCProtonTrk2;
		}
	

	      //// If I am here then I have both partners. Lets fill correlator Histograms:
	   


	      //Remove EP-POIs AutoCorrelation: only for TPC EP
	      if(trk2Eta*trk1Eta < 0){
		
		if(trk1Eta >= 0){
		  localSumQ2x = fSumQnxNeg[0]; /// We need the full Q-sum. Then remove only current track-2
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
		
		if(trk2Pt < 2.0 && !((trk2Eta >= fEtaGapNeg) && (trk2Eta <= fEtaGapPos))){ // 
		  localSumQ2x -= WgtNUAChtrk2*trk2Pt*TMath::Cos(2*trk2Phi);
		  localSumQ2y -= WgtNUAChtrk2*trk2Pt*TMath::Sin(2*trk2Phi);
		  localSumQ3x -= WgtNUAChtrk2*trk2Pt*TMath::Cos(3*trk2Phi);
		  localSumQ3y -= WgtNUAChtrk2*trk2Pt*TMath::Sin(3*trk2Phi);
		  localMultTPC-= WgtNUAChtrk2*trk2Pt;                        /// Rihan Todo: Remove PtWeights from Q sum and MultSum!
		}
		
		if(localMultTPC>0){
		  localSumQ2x = localSumQ2x/localMultTPC;
		  localSumQ2y = localSumQ2y/localMultTPC;
		}
		
		fPsiNTPCPos = (1./2)*TMath::ATan2(localSumQ2y,localSumQ2x);
		if(fPsiNTPCPos < 0) fPsiNTPCPos += TMath::TwoPi()/2;
		fPsi3TPCPos = (1./3)*TMath::ATan2(localSumQ3y,localSumQ3x);
		if(fPsi3TPCPos < 0) fPsi3TPCPos += TMath::TwoPi()/3;
		//fPsi4TPCPos = (1./4)*TMath::ATan2(fQ4yPos,fQ4xPos);
		//if(fPsi4TPCPos < 0) fPsi4TPCPos += TMath::TwoPi()/4;  		
	      }

		// force using TPC to save the results
		fPsiNEvent = fPsiNTPCPos;
		fPsi3Event = fPsi3TPCPos;
	      
	      if(trk1Chrg*trk2Chrg < 0){ //Opposite sign	
		    hAvg3pC112vsCentOS->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2Ch);
		    hAvg3pC123vsCentOS->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2Ch);
				
		    hAvgDelta1vsCentOS->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2Ch);
		    hAvgDelta2vsCentOS->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2Ch);
		    hAvgDelta3vsCentOS->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2Ch);
		    hAvgDelta4vsCentOS->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2Ch);		
		
			if(trk1Pt <= 2.0 && trk2Pt <= 2.0) {
				if(isTrk1Pion && isTrk2Pion) {
					hAvg3pC112vsCentOSPionPion->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvg3pC123vsCentOSPionPion->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2PIDPion);
							
					hAvgDelta1vsCentOSPionPion->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta2vsCentOSPionPion->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta3vsCentOSPionPion->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta4vsCentOSPionPion->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);		
				} else if(isTrk1Kaon && isTrk2Kaon) {
					hAvg3pC112vsCentOSKaonKaon->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvg3pC123vsCentOSKaonKaon->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2PIDKaon);
							
					hAvgDelta1vsCentOSKaonKaon->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta2vsCentOSKaonKaon->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta3vsCentOSKaonKaon->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta4vsCentOSKaonKaon->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);	
				} else if(isTrk1Proton && isTrk2Proton) {
					hAvg3pC112vsCentOSProtonProton->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvg3pC123vsCentOSProtonProton->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2PIDProton);
							
					hAvgDelta1vsCentOSProtonProton->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta2vsCentOSProtonProton->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta3vsCentOSProtonProton->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta4vsCentOSProtonProton->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);	
				}
					
				if(isTrk1Pion || isTrk2Pion) {
					if(isTrk1Pion) {
						hAvg3pC112vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2Ch);
						hAvg3pC123vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2Ch);
								
						hAvgDelta1vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta2vsCentOSPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta3vsCentOSPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta4vsCentOSPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);		
					} else if(isTrk2Pion) {
						hAvg3pC112vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDPion);
						hAvg3pC123vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDPion);
								
						hAvgDelta1vsCentOSPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta2vsCentOSPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta3vsCentOSPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta4vsCentOSPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);	
					}
				} else if(isTrk1Kaon || isTrk2Kaon) {
					if(isTrk1Kaon) {
						hAvg3pC112vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2Ch);
						hAvg3pC123vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2Ch);
								
						hAvgDelta1vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta2vsCentOSKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta3vsCentOSKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta4vsCentOSKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);	
					} else if(isTrk2Kaon) {
						hAvg3pC112vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDKaon);
						hAvg3pC123vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDKaon);
								
						hAvgDelta1vsCentOSKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta2vsCentOSKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta3vsCentOSKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta4vsCentOSKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);	
					}
				} else if(isTrk1Proton || isTrk2Proton) {
					if(isTrk1Proton) {
						hAvg3pC112vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2Ch);
						hAvg3pC123vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2Ch);
								
						hAvgDelta1vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta2vsCentOSProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta3vsCentOSProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta4vsCentOSProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);	
					} else if(isTrk2Proton) {
						hAvg3pC112vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDProton);
						hAvg3pC123vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDProton);
								
						hAvgDelta1vsCentOSProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta2vsCentOSProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta3vsCentOSProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta4vsCentOSProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);	
					}
					
				}
		    }
		
	      }		
	      else if(trk1Chrg > 0 && trk2Chrg > 0){		      
		    hAvg3pC112vsCentPP->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb2Ch);
		    hAvg3pC123vsCentPP->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb2Ch);
				
		    hAvgDelta1vsCentPP->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb2Ch);
		    hAvgDelta2vsCentPP->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb2Ch);
		    hAvgDelta3vsCentPP->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb2Ch);
		    hAvgDelta4vsCentPP->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb2Ch);		
		    
		    if(trk1Pt <= 2.0 && trk2Pt <= 2.0) {
				if(isTrk1Pion && isTrk2Pion) {
					hAvg3pC112vsCentPPPionPion->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvg3pC123vsCentPPPionPion->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2PIDPion);
							
					hAvgDelta1vsCentPPPionPion->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta2vsCentPPPionPion->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta3vsCentPPPionPion->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta4vsCentPPPionPion->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);		
				} else if(isTrk1Kaon && isTrk2Kaon) {
					hAvg3pC112vsCentPPKaonKaon->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvg3pC123vsCentPPKaonKaon->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2PIDKaon);
							
					hAvgDelta1vsCentPPKaonKaon->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta2vsCentPPKaonKaon->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta3vsCentPPKaonKaon->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta4vsCentPPKaonKaon->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);	
				} else if(isTrk1Proton && isTrk2Proton) {
					hAvg3pC112vsCentPPProtonProton->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvg3pC123vsCentPPProtonProton->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2PIDProton);
							
					hAvgDelta1vsCentPPProtonProton->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta2vsCentPPProtonProton->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta3vsCentPPProtonProton->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta4vsCentPPProtonProton->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);	
				}
					
				if(isTrk1Pion || isTrk2Pion) {
					if(isTrk1Pion) {
						hAvg3pC112vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2Ch);
						hAvg3pC123vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2Ch);
								
						hAvgDelta1vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta2vsCentPPPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta3vsCentPPPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta4vsCentPPPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);		
					} else if(isTrk2Pion) {
						hAvg3pC112vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDPion);
						hAvg3pC123vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDPion);
								
						hAvgDelta1vsCentPPPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta2vsCentPPPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta3vsCentPPPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta4vsCentPPPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);	
					}
				} else if(isTrk1Kaon || isTrk2Kaon) {
					if(isTrk1Kaon) {
						hAvg3pC112vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2Ch);
						hAvg3pC123vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2Ch);
								
						hAvgDelta1vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta2vsCentPPKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta3vsCentPPKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta4vsCentPPKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);	
					} else if(isTrk2Kaon) {
						hAvg3pC112vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDKaon);
						hAvg3pC123vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDKaon);
								
						hAvgDelta1vsCentPPKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta2vsCentPPKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta3vsCentPPKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta4vsCentPPKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);	
					}
				} else if(isTrk1Proton || isTrk2Proton) {
					if(isTrk1Proton) {
						hAvg3pC112vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2Ch);
						hAvg3pC123vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2Ch);
								
						hAvgDelta1vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta2vsCentPPProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta3vsCentPPProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta4vsCentPPProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);	
					} else if(isTrk2Proton) {
						hAvg3pC112vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDProton);
						hAvg3pC123vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDProton);
								
						hAvgDelta1vsCentPPProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta2vsCentPPProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta3vsCentPPProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta4vsCentPPProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);	
					}
				}
			}
	      }
	      //else if(trk1Chrg < 0 && trk2Chrg < 0){  ///this is obvious!
	      else{
		    hAvg3pC112vsCentNN->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb2Ch);
		    hAvg3pC123vsCentNN->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb2Ch);
		
		    hAvgDelta1vsCentNN->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb2Ch);
		    hAvgDelta2vsCentNN->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb2Ch);
		    hAvgDelta3vsCentNN->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb2Ch);
		    hAvgDelta4vsCentNN->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb2Ch);
		    
		    if(trk1Pt <= 2.0 && trk2Pt <= 2.0) {
				if(isTrk1Pion && isTrk2Pion) {
					hAvg3pC112vsCentNNPionPion->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvg3pC123vsCentNNPionPion->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2PIDPion);
							
					hAvgDelta1vsCentNNPionPion->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta2vsCentNNPionPion->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta3vsCentNNPionPion->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);
					hAvgDelta4vsCentNNPionPion->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2PIDPion);		
				} else if(isTrk1Kaon && isTrk2Kaon) {
					hAvg3pC112vsCentNNKaonKaon->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvg3pC123vsCentNNKaonKaon->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2PIDKaon);
							
					hAvgDelta1vsCentNNKaonKaon->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta2vsCentNNKaonKaon->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta3vsCentNNKaonKaon->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);
					hAvgDelta4vsCentNNKaonKaon->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2PIDKaon);	
				} else if(isTrk1Proton && isTrk2Proton) {
					hAvg3pC112vsCentNNProtonProton->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvg3pC123vsCentNNProtonProton->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2PIDProton);
							
					hAvgDelta1vsCentNNProtonProton->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta2vsCentNNProtonProton->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta3vsCentNNProtonProton->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);
					hAvgDelta4vsCentNNProtonProton->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2PIDProton);	
				}
					
				if(isTrk1Pion || isTrk2Pion) {
					if(isTrk1Pion) {
						hAvg3pC112vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDPion*wgtComb2Ch);
						hAvg3pC123vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDPion*wgtComb2Ch);
								
						hAvgDelta1vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta2vsCentNNPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta3vsCentNNPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);
						hAvgDelta4vsCentNNPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDPion*wgtComb2Ch);		
					} else if(isTrk2Pion) {
						hAvg3pC112vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDPion);
						hAvg3pC123vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDPion);
								
						hAvgDelta1vsCentNNPionCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta2vsCentNNPionCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta3vsCentNNPionCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);
						hAvgDelta4vsCentNNPionCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDPion);	
					}
				} else if(isTrk1Kaon || isTrk2Kaon) {
					if(isTrk1Kaon) {
						hAvg3pC112vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDKaon*wgtComb2Ch);
						hAvg3pC123vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDKaon*wgtComb2Ch);
								
						hAvgDelta1vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta2vsCentNNKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta3vsCentNNKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);
						hAvgDelta4vsCentNNKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDKaon*wgtComb2Ch);	
					} else if(isTrk2Kaon) {
						hAvg3pC112vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDKaon);
						hAvg3pC123vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDKaon);
								
						hAvgDelta1vsCentNNKaonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta2vsCentNNKaonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta3vsCentNNKaonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);
						hAvgDelta4vsCentNNKaonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDKaon);	
					}
				} else if(isTrk1Proton || isTrk2Proton) {
					if(isTrk1Proton) {
						hAvg3pC112vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1PIDProton*wgtComb2Ch);
						hAvg3pC123vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1PIDProton*wgtComb2Ch);
								
						hAvgDelta1vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta2vsCentNNProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta3vsCentNNProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);
						hAvgDelta4vsCentNNProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1PIDProton*wgtComb2Ch);	
					} else if(isTrk2Proton) {
						hAvg3pC112vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgtComb1Ch*wgtComb2PIDProton);
						hAvg3pC123vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgtComb1Ch*wgtComb2PIDProton);
								
						hAvgDelta1vsCentNNProtonCharge->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta2vsCentNNProtonCharge->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta3vsCentNNProtonCharge->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);
						hAvgDelta4vsCentNNProtonCharge->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgtComb1Ch*wgtComb2PIDProton);	
					}
				}
			}
	      }
	      
	      
	    }//j-track cuts
	  }//j-track FB validated
	}///j-track loop
	

      }//----> i-track loop => All trackCuts applied.     
    }//-----> i-track loop => FB is validated.    
  }///-----> i-track loop Ends here.<--------

  CalculateCMESPPP();
  
  // fill event information
  treeEvent->Fill();
  
  // Reset event-by-event variables
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      if(fCMEQRe[c][h]) fCMEQRe[c][h]->Reset();
      if(fCMEQIm[c][h]) fCMEQIm[c][h]->Reset();
      if(fCMEMult[c][h]) fCMEMult[c][h]->Reset();
      if(fCMEQRePion[c][h]) fCMEQRePion[c][h]->Reset();
      if(fCMEQImPion[c][h]) fCMEQImPion[c][h]->Reset();
      if(fCMEMultPion[c][h]) fCMEMultPion[c][h]->Reset();
      if(fCMEQReKaon[c][h]) fCMEQReKaon[c][h]->Reset();
      if(fCMEQImKaon[c][h]) fCMEQImKaon[c][h]->Reset();
      if(fCMEMultKaon[c][h]) fCMEMultKaon[c][h]->Reset();
      if(fCMEQReProton[c][h]) fCMEQReProton[c][h]->Reset();
      if(fCMEQImProton[c][h]) fCMEQImProton[c][h]->Reset();
      if(fCMEMultProton[c][h]) fCMEMultProton[c][h]->Reset();
    }
  }
  //@shi reset CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      if(fCMEQReBothCharge[c][h]) fCMEQReBothCharge[c][h]->Reset();
      if(fCMEQImBothCharge[c][h]) fCMEQImBothCharge[c][h]->Reset();
      if(fCMEMultBothCharge[c][h]) fCMEMultBothCharge[c][h]->Reset();
      if(fCMEQRePionBothCharge[c][h]) fCMEQRePionBothCharge[c][h]->Reset();
      if(fCMEQImPionBothCharge[c][h]) fCMEQImPionBothCharge[c][h]->Reset();
      if(fCMEMultPionBothCharge[c][h]) fCMEMultPionBothCharge[c][h]->Reset();
      if(fCMEQReKaonBothCharge[c][h]) fCMEQReKaonBothCharge[c][h]->Reset();
      if(fCMEQImKaonBothCharge[c][h]) fCMEQImKaonBothCharge[c][h]->Reset();
      if(fCMEMultKaonBothCharge[c][h]) fCMEMultKaonBothCharge[c][h]->Reset();
      if(fCMEQReProtonBothCharge[c][h]) fCMEQReProtonBothCharge[c][h]->Reset();
      if(fCMEQImProtonBothCharge[c][h]) fCMEQImProtonBothCharge[c][h]->Reset();
      if(fCMEMultProtonBothCharge[c][h]) fCMEMultProtonBothCharge[c][h]->Reset();
    }
  }
  
  for(Int_t c=0;c<2;c++) {
	if(fCMEQ2Re4[c]) fCMEQ2Re4[c]->Reset(); // w^2*cos(4phi)
    if(fCMEQ3Re2[c]) fCMEQ3Re2[c]->Reset(); // w^3*cos(2phi)
    if(fCMEQ2Im4[c]) fCMEQ2Im4[c]->Reset(); // w^2*sin(4phi)
    if(fCMEQ3Im2[c]) fCMEQ3Im2[c]->Reset(); // w^3*sin(2phi)
    if(fCMEw0[c]) fCMEw0[c]->Reset();    // w^0
    if(fCMEw1[c]) fCMEw1[c]->Reset();    // w^1
    if(fCMEw2[c]) fCMEw2[c]->Reset();    // w^2
    if(fCMEw3[c]) fCMEw3[c]->Reset();    // w^3
    if(fCMEw4[c]) fCMEw4[c]->Reset();    // w^4
  
	if(fCMEQ2Re4Pion[c]) fCMEQ2Re4Pion[c]->Reset(); // w^2*cos(4phi)
    if(fCMEQ3Re2Pion[c]) fCMEQ3Re2Pion[c]->Reset(); // w^3*cos(2phi)
    if(fCMEQ2Im4Pion[c]) fCMEQ2Im4Pion[c]->Reset(); // w^2*sin(4phi)
    if(fCMEQ3Im2Pion[c]) fCMEQ3Im2Pion[c]->Reset(); // w^3*sin(2phi)
    if(fCMEw0Pion[c]) fCMEw0Pion[c]->Reset();    // w^0
    if(fCMEw1Pion[c]) fCMEw1Pion[c]->Reset();    // w^1
    if(fCMEw2Pion[c]) fCMEw2Pion[c]->Reset();    // w^2
    if(fCMEw3Pion[c]) fCMEw3Pion[c]->Reset();    // w^3
    if(fCMEw4Pion[c]) fCMEw4Pion[c]->Reset();    // w^4
    
    if(fCMEQ2Re4Kaon[c]) fCMEQ2Re4Kaon[c]->Reset(); // w^2*cos(4phi)
    if(fCMEQ3Re2Kaon[c]) fCMEQ3Re2Kaon[c]->Reset(); // w^3*cos(2phi)
    if(fCMEQ2Im4Kaon[c]) fCMEQ2Im4Kaon[c]->Reset(); // w^2*sin(4phi)
    if(fCMEQ3Im2Kaon[c]) fCMEQ3Im2Kaon[c]->Reset(); // w^3*sin(2phi)
    if(fCMEw0Kaon[c]) fCMEw0Kaon[c]->Reset();    // w^0
    if(fCMEw1Kaon[c]) fCMEw1Kaon[c]->Reset();    // w^1
    if(fCMEw2Kaon[c]) fCMEw2Kaon[c]->Reset();    // w^2
    if(fCMEw3Kaon[c]) fCMEw3Kaon[c]->Reset();    // w^3
    if(fCMEw4Kaon[c]) fCMEw4Kaon[c]->Reset();    // w^4
    
    if(fCMEQ2Re4Proton[c]) fCMEQ2Re4Proton[c]->Reset(); // w^2*cos(4phi)
    if(fCMEQ3Re2Proton[c]) fCMEQ3Re2Proton[c]->Reset(); // w^3*cos(2phi)
    if(fCMEQ2Im4Proton[c]) fCMEQ2Im4Proton[c]->Reset(); // w^2*sin(4phi)
    if(fCMEQ3Im2Proton[c]) fCMEQ3Im2Proton[c]->Reset(); // w^3*sin(2phi)
    if(fCMEw0Proton[c]) fCMEw0Proton[c]->Reset();    // w^0
    if(fCMEw1Proton[c]) fCMEw1Proton[c]->Reset();    // w^1
    if(fCMEw2Proton[c]) fCMEw2Proton[c]->Reset();    // w^2
    if(fCMEw3Proton[c]) fCMEw3Proton[c]->Reset();    // w^3
    if(fCMEw4Proton[c]) fCMEw4Proton[c]->Reset();    // w^4
  }
  
  
  fDebugwEventCount->Fill(9.1); ///Left for Analysis
  
  fHistVertexZcm->Fill(pVtxZ);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

}//---------------- UserExec ----------------------



void AliAnalysisTaskGammaDeltaPIDSaveQvec::SetupAnalysisHistograms(){
  
  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables
  Char_t  name[100];
  Char_t title[100];

  //// V2
  hAvgV2TPCvsCent = new TProfile("hAvgV2TPCvsCent","v2",18,0,90);
  hAvgV2TPCvsCent->Sumw2();
  fListHist->Add(hAvgV2TPCvsCent);
  hAvgV2TPCvsCentPion = new TProfile("hAvgV2TPCvsCentPion","v2 pion",18,0,90);
  hAvgV2TPCvsCentPion->Sumw2();
  fListHist->Add(hAvgV2TPCvsCentPion);
  hAvgV2TPCvsCentKaon = new TProfile("hAvgV2TPCvsCentKaon","v2 kaon",18,0,90);
  hAvgV2TPCvsCentKaon->Sumw2();
  fListHist->Add(hAvgV2TPCvsCentKaon);
  hAvgV2TPCvsCentProton = new TProfile("hAvgV2TPCvsCentProton","v2 proton",18,0,90);
  hAvgV2TPCvsCentProton->Sumw2();
  fListHist->Add(hAvgV2TPCvsCentProton);

  //// CME PID ANALYSIS Histograms/Profiles:

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
  
   
  /// Pion, Kaon, Proton PID correlators
  hAvg3pC112vsCentOSPionPion = new TProfile("hAvg3pC112vsCentOSPionPion"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSPionPion = new TProfile("hAvg3pC123vsCentOSPionPion"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSPionPion = new TProfile("hAvgDelta1vsCentOSPionPion"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSPionPion = new TProfile("hAvgDelta2vsCentOSPionPion"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSPionPion = new TProfile("hAvgDelta3vsCentOSPionPion"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSPionPion = new TProfile("hAvgDelta4vsCentOSPionPion"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentOSKaonKaon = new TProfile("hAvg3pC112vsCentOSKaonKaon"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSKaonKaon = new TProfile("hAvg3pC123vsCentOSKaonKaon"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSKaonKaon = new TProfile("hAvgDelta1vsCentOSKaonKaon"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSKaonKaon = new TProfile("hAvgDelta2vsCentOSKaonKaon"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSKaonKaon = new TProfile("hAvgDelta3vsCentOSKaonKaon"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSKaonKaon = new TProfile("hAvgDelta4vsCentOSKaonKaon"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentOSProtonProton = new TProfile("hAvg3pC112vsCentOSProtonProton"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSProtonProton = new TProfile("hAvg3pC123vsCentOSProtonProton"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSProtonProton = new TProfile("hAvgDelta1vsCentOSProtonProton"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSProtonProton = new TProfile("hAvgDelta2vsCentOSProtonProton"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSProtonProton = new TProfile("hAvgDelta3vsCentOSProtonProton"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSProtonProton = new TProfile("hAvgDelta4vsCentOSProtonProton"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);

  hAvg3pC112vsCentOSPionPion->Sumw2();
  hAvg3pC123vsCentOSPionPion->Sumw2();
  hAvgDelta1vsCentOSPionPion->Sumw2();
  hAvgDelta2vsCentOSPionPion->Sumw2();
  hAvgDelta3vsCentOSPionPion->Sumw2();
  hAvgDelta4vsCentOSPionPion->Sumw2();
  hAvg3pC112vsCentOSKaonKaon->Sumw2();
  hAvg3pC123vsCentOSKaonKaon->Sumw2();
  hAvgDelta1vsCentOSKaonKaon->Sumw2();
  hAvgDelta2vsCentOSKaonKaon->Sumw2();
  hAvgDelta3vsCentOSKaonKaon->Sumw2();
  hAvgDelta4vsCentOSKaonKaon->Sumw2();
  hAvg3pC112vsCentOSProtonProton->Sumw2();
  hAvg3pC123vsCentOSProtonProton->Sumw2();
  hAvgDelta1vsCentOSProtonProton->Sumw2();
  hAvgDelta2vsCentOSProtonProton->Sumw2();
  hAvgDelta3vsCentOSProtonProton->Sumw2();
  hAvgDelta4vsCentOSProtonProton->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentOSPionPion);
  fListHist->Add(hAvg3pC123vsCentOSPionPion);
  fListHist->Add(hAvgDelta1vsCentOSPionPion);
  fListHist->Add(hAvgDelta2vsCentOSPionPion);
  fListHist->Add(hAvgDelta3vsCentOSPionPion);
  fListHist->Add(hAvgDelta4vsCentOSPionPion);
  fListHist->Add(hAvg3pC112vsCentOSKaonKaon);
  fListHist->Add(hAvg3pC123vsCentOSKaonKaon);
  fListHist->Add(hAvgDelta1vsCentOSKaonKaon);
  fListHist->Add(hAvgDelta2vsCentOSKaonKaon);
  fListHist->Add(hAvgDelta3vsCentOSKaonKaon);
  fListHist->Add(hAvgDelta4vsCentOSKaonKaon);
  fListHist->Add(hAvg3pC112vsCentOSProtonProton);
  fListHist->Add(hAvg3pC123vsCentOSProtonProton);
  fListHist->Add(hAvgDelta1vsCentOSProtonProton);
  fListHist->Add(hAvgDelta2vsCentOSProtonProton);
  fListHist->Add(hAvgDelta3vsCentOSProtonProton);
  fListHist->Add(hAvgDelta4vsCentOSProtonProton);
  
  hAvg3pC112vsCentOSPionCharge = new TProfile("hAvg3pC112vsCentOSPionCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSPionCharge = new TProfile("hAvg3pC123vsCentOSPionCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSPionCharge = new TProfile("hAvgDelta1vsCentOSPionCharge"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSPionCharge = new TProfile("hAvgDelta2vsCentOSPionCharge"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSPionCharge = new TProfile("hAvgDelta3vsCentOSPionCharge"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSPionCharge = new TProfile("hAvgDelta4vsCentOSPionCharge"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentOSKaonCharge = new TProfile("hAvg3pC112vsCentOSKaonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSKaonCharge = new TProfile("hAvg3pC123vsCentOSKaonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSKaonCharge = new TProfile("hAvgDelta1vsCentOSKaonCharge"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSKaonCharge = new TProfile("hAvgDelta2vsCentOSKaonCharge"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSKaonCharge = new TProfile("hAvgDelta3vsCentOSKaonCharge"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSKaonCharge = new TProfile("hAvgDelta4vsCentOSKaonCharge"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentOSProtonCharge = new TProfile("hAvg3pC112vsCentOSProtonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentOSProtonCharge = new TProfile("hAvg3pC123vsCentOSProtonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentOSProtonCharge = new TProfile("hAvgDelta1vsCentOSProtonCharge"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentOSProtonCharge = new TProfile("hAvgDelta2vsCentOSProtonCharge"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentOSProtonCharge = new TProfile("hAvgDelta3vsCentOSProtonCharge"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentOSProtonCharge = new TProfile("hAvgDelta4vsCentOSProtonCharge"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  
  hAvg3pC112vsCentOSPionCharge->Sumw2();
  hAvg3pC123vsCentOSPionCharge->Sumw2();
  hAvgDelta1vsCentOSPionCharge->Sumw2();
  hAvgDelta2vsCentOSPionCharge->Sumw2();
  hAvgDelta3vsCentOSPionCharge->Sumw2();
  hAvgDelta4vsCentOSPionCharge->Sumw2();
  hAvg3pC112vsCentOSKaonCharge->Sumw2();
  hAvg3pC123vsCentOSKaonCharge->Sumw2();
  hAvgDelta1vsCentOSKaonCharge->Sumw2();
  hAvgDelta2vsCentOSKaonCharge->Sumw2();
  hAvgDelta3vsCentOSKaonCharge->Sumw2();
  hAvgDelta4vsCentOSKaonCharge->Sumw2();
  hAvg3pC112vsCentOSProtonCharge->Sumw2();
  hAvg3pC123vsCentOSProtonCharge->Sumw2();
  hAvgDelta1vsCentOSProtonCharge->Sumw2();
  hAvgDelta2vsCentOSProtonCharge->Sumw2();
  hAvgDelta3vsCentOSProtonCharge->Sumw2();
  hAvgDelta4vsCentOSProtonCharge->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentOSPionCharge);
  fListHist->Add(hAvg3pC123vsCentOSPionCharge);
  fListHist->Add(hAvgDelta1vsCentOSPionCharge);
  fListHist->Add(hAvgDelta2vsCentOSPionCharge);
  fListHist->Add(hAvgDelta3vsCentOSPionCharge);
  fListHist->Add(hAvgDelta4vsCentOSPionCharge);
  fListHist->Add(hAvg3pC112vsCentOSKaonCharge);
  fListHist->Add(hAvg3pC123vsCentOSKaonCharge);
  fListHist->Add(hAvgDelta1vsCentOSKaonCharge);
  fListHist->Add(hAvgDelta2vsCentOSKaonCharge);
  fListHist->Add(hAvgDelta3vsCentOSKaonCharge);
  fListHist->Add(hAvgDelta4vsCentOSKaonCharge);
  fListHist->Add(hAvg3pC112vsCentOSProtonCharge);
  fListHist->Add(hAvg3pC123vsCentOSProtonCharge);
  fListHist->Add(hAvgDelta1vsCentOSProtonCharge);
  fListHist->Add(hAvgDelta2vsCentOSProtonCharge);
  fListHist->Add(hAvgDelta3vsCentOSProtonCharge);
  fListHist->Add(hAvgDelta4vsCentOSProtonCharge);

  hAvg3pC112vsCentPPPionPion = new TProfile("hAvg3pC112vsCentPPPionPion"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPPionPion = new TProfile("hAvg3pC123vsCentPPPionPion"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPPionPion = new TProfile("hAvgDelta1vsCentPPPionPion"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPPionPion = new TProfile("hAvgDelta2vsCentPPPionPion"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPPionPion = new TProfile("hAvgDelta3vsCentPPPionPion"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPPionPion = new TProfile("hAvgDelta4vsCentPPPionPion"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentPPKaonKaon = new TProfile("hAvg3pC112vsCentPPKaonKaon"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPKaonKaon = new TProfile("hAvg3pC123vsCentPPKaonKaon"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPKaonKaon = new TProfile("hAvgDelta1vsCentPPKaonKaon"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPKaonKaon = new TProfile("hAvgDelta2vsCentPPKaonKaon"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPKaonKaon = new TProfile("hAvgDelta3vsCentPPKaonKaon"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPKaonKaon = new TProfile("hAvgDelta4vsCentPPKaonKaon"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentPPProtonProton = new TProfile("hAvg3pC112vsCentPPProtonProton"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPProtonProton = new TProfile("hAvg3pC123vsCentPPProtonProton"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPProtonProton = new TProfile("hAvgDelta1vsCentPPProtonProton"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPProtonProton = new TProfile("hAvgDelta2vsCentPPProtonProton"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPProtonProton = new TProfile("hAvgDelta3vsCentPPProtonProton"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPProtonProton = new TProfile("hAvgDelta4vsCentPPProtonProton"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);

  hAvg3pC112vsCentPPPionPion->Sumw2();
  hAvg3pC123vsCentPPPionPion->Sumw2();
  hAvgDelta1vsCentPPPionPion->Sumw2();
  hAvgDelta2vsCentPPPionPion->Sumw2();
  hAvgDelta3vsCentPPPionPion->Sumw2();
  hAvgDelta4vsCentPPPionPion->Sumw2();
  hAvg3pC112vsCentPPKaonKaon->Sumw2();
  hAvg3pC123vsCentPPKaonKaon->Sumw2();
  hAvgDelta1vsCentPPKaonKaon->Sumw2();
  hAvgDelta2vsCentPPKaonKaon->Sumw2();
  hAvgDelta3vsCentPPKaonKaon->Sumw2();
  hAvgDelta4vsCentPPKaonKaon->Sumw2();
  hAvg3pC112vsCentPPProtonProton->Sumw2();
  hAvg3pC123vsCentPPProtonProton->Sumw2();
  hAvgDelta1vsCentPPProtonProton->Sumw2();
  hAvgDelta2vsCentPPProtonProton->Sumw2();
  hAvgDelta3vsCentPPProtonProton->Sumw2();
  hAvgDelta4vsCentPPProtonProton->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentPPPionPion);
  fListHist->Add(hAvg3pC123vsCentPPPionPion);
  fListHist->Add(hAvgDelta1vsCentPPPionPion);
  fListHist->Add(hAvgDelta2vsCentPPPionPion);
  fListHist->Add(hAvgDelta3vsCentPPPionPion);
  fListHist->Add(hAvgDelta4vsCentPPPionPion);
  fListHist->Add(hAvg3pC112vsCentPPKaonKaon);
  fListHist->Add(hAvg3pC123vsCentPPKaonKaon);
  fListHist->Add(hAvgDelta1vsCentPPKaonKaon);
  fListHist->Add(hAvgDelta2vsCentPPKaonKaon);
  fListHist->Add(hAvgDelta3vsCentPPKaonKaon);
  fListHist->Add(hAvgDelta4vsCentPPKaonKaon);
  fListHist->Add(hAvg3pC112vsCentPPProtonProton);
  fListHist->Add(hAvg3pC123vsCentPPProtonProton);
  fListHist->Add(hAvgDelta1vsCentPPProtonProton);
  fListHist->Add(hAvgDelta2vsCentPPProtonProton);
  fListHist->Add(hAvgDelta3vsCentPPProtonProton);
  fListHist->Add(hAvgDelta4vsCentPPProtonProton);
  
  hAvg3pC112vsCentPPPionCharge = new TProfile("hAvg3pC112vsCentPPPionCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPPionCharge = new TProfile("hAvg3pC123vsCentPPPionCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPPionCharge = new TProfile("hAvgDelta1vsCentPPPionCharge"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPPionCharge = new TProfile("hAvgDelta2vsCentPPPionCharge"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPPionCharge = new TProfile("hAvgDelta3vsCentPPPionCharge"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPPionCharge = new TProfile("hAvgDelta4vsCentPPPionCharge"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentPPKaonCharge = new TProfile("hAvg3pC112vsCentPPKaonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPKaonCharge = new TProfile("hAvg3pC123vsCentPPKaonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPKaonCharge = new TProfile("hAvgDelta1vsCentPPKaonCharge"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPKaonCharge = new TProfile("hAvgDelta2vsCentPPKaonCharge"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPKaonCharge = new TProfile("hAvgDelta3vsCentPPKaonCharge"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPKaonCharge = new TProfile("hAvgDelta4vsCentPPKaonCharge"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentPPProtonCharge = new TProfile("hAvg3pC112vsCentPPProtonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentPPProtonCharge = new TProfile("hAvg3pC123vsCentPPProtonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentPPProtonCharge = new TProfile("hAvgDelta1vsCentPPProtonCharge"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentPPProtonCharge = new TProfile("hAvgDelta2vsCentPPProtonCharge"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentPPProtonCharge = new TProfile("hAvgDelta3vsCentPPProtonCharge"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentPPProtonCharge = new TProfile("hAvgDelta4vsCentPPProtonCharge"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  
  hAvg3pC112vsCentPPPionCharge->Sumw2();
  hAvg3pC123vsCentPPPionCharge->Sumw2();
  hAvgDelta1vsCentPPPionCharge->Sumw2();
  hAvgDelta2vsCentPPPionCharge->Sumw2();
  hAvgDelta3vsCentPPPionCharge->Sumw2();
  hAvgDelta4vsCentPPPionCharge->Sumw2();
  hAvg3pC112vsCentPPKaonCharge->Sumw2();
  hAvg3pC123vsCentPPKaonCharge->Sumw2();
  hAvgDelta1vsCentPPKaonCharge->Sumw2();
  hAvgDelta2vsCentPPKaonCharge->Sumw2();
  hAvgDelta3vsCentPPKaonCharge->Sumw2();
  hAvgDelta4vsCentPPKaonCharge->Sumw2();
  hAvg3pC112vsCentPPProtonCharge->Sumw2();
  hAvg3pC123vsCentPPProtonCharge->Sumw2();
  hAvgDelta1vsCentPPProtonCharge->Sumw2();
  hAvgDelta2vsCentPPProtonCharge->Sumw2();
  hAvgDelta3vsCentPPProtonCharge->Sumw2();
  hAvgDelta4vsCentPPProtonCharge->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentPPPionCharge);
  fListHist->Add(hAvg3pC123vsCentPPPionCharge);
  fListHist->Add(hAvgDelta1vsCentPPPionCharge);
  fListHist->Add(hAvgDelta2vsCentPPPionCharge);
  fListHist->Add(hAvgDelta3vsCentPPPionCharge);
  fListHist->Add(hAvgDelta4vsCentPPPionCharge);
  fListHist->Add(hAvg3pC112vsCentPPKaonCharge);
  fListHist->Add(hAvg3pC123vsCentPPKaonCharge);
  fListHist->Add(hAvgDelta1vsCentPPKaonCharge);
  fListHist->Add(hAvgDelta2vsCentPPKaonCharge);
  fListHist->Add(hAvgDelta3vsCentPPKaonCharge);
  fListHist->Add(hAvgDelta4vsCentPPKaonCharge);
  fListHist->Add(hAvg3pC112vsCentPPProtonCharge);
  fListHist->Add(hAvg3pC123vsCentPPProtonCharge);
  fListHist->Add(hAvgDelta1vsCentPPProtonCharge);
  fListHist->Add(hAvgDelta2vsCentPPProtonCharge);
  fListHist->Add(hAvgDelta3vsCentPPProtonCharge);
  fListHist->Add(hAvgDelta4vsCentPPProtonCharge);
  
  hAvg3pC112vsCentNNPionPion = new TProfile("hAvg3pC112vsCentNNPionPion"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNPionPion = new TProfile("hAvg3pC123vsCentNNPionPion"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNPionPion = new TProfile("hAvgDelta1vsCentNNPionPion"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNPionPion = new TProfile("hAvgDelta2vsCentNNPionPion"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNPionPion = new TProfile("hAvgDelta3vsCentNNPionPion"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNPionPion = new TProfile("hAvgDelta4vsCentNNPionPion"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentNNKaonKaon = new TProfile("hAvg3pC112vsCentNNKaonKaon"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNKaonKaon = new TProfile("hAvg3pC123vsCentNNKaonKaon"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNKaonKaon = new TProfile("hAvgDelta1vsCentNNKaonKaon"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNKaonKaon = new TProfile("hAvgDelta2vsCentNNKaonKaon"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNKaonKaon = new TProfile("hAvgDelta3vsCentNNKaonKaon"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNKaonKaon = new TProfile("hAvgDelta4vsCentNNKaonKaon"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentNNProtonProton = new TProfile("hAvg3pC112vsCentNNProtonProton"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNProtonProton = new TProfile("hAvg3pC123vsCentNNProtonProton"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNProtonProton = new TProfile("hAvgDelta1vsCentNNProtonProton"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNProtonProton = new TProfile("hAvgDelta2vsCentNNProtonProton"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNProtonProton = new TProfile("hAvgDelta3vsCentNNProtonProton"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNProtonProton = new TProfile("hAvgDelta4vsCentNNProtonProton"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);

  hAvg3pC112vsCentNNPionPion->Sumw2();
  hAvg3pC123vsCentNNPionPion->Sumw2();
  hAvgDelta1vsCentNNPionPion->Sumw2();
  hAvgDelta2vsCentNNPionPion->Sumw2();
  hAvgDelta3vsCentNNPionPion->Sumw2();
  hAvgDelta4vsCentNNPionPion->Sumw2();
  hAvg3pC112vsCentNNKaonKaon->Sumw2();
  hAvg3pC123vsCentNNKaonKaon->Sumw2();
  hAvgDelta1vsCentNNKaonKaon->Sumw2();
  hAvgDelta2vsCentNNKaonKaon->Sumw2();
  hAvgDelta3vsCentNNKaonKaon->Sumw2();
  hAvgDelta4vsCentNNKaonKaon->Sumw2();
  hAvg3pC112vsCentNNProtonProton->Sumw2();
  hAvg3pC123vsCentNNProtonProton->Sumw2();
  hAvgDelta1vsCentNNProtonProton->Sumw2();
  hAvgDelta2vsCentNNProtonProton->Sumw2();
  hAvgDelta3vsCentNNProtonProton->Sumw2();
  hAvgDelta4vsCentNNProtonProton->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentNNPionPion);
  fListHist->Add(hAvg3pC123vsCentNNPionPion);
  fListHist->Add(hAvgDelta1vsCentNNPionPion);
  fListHist->Add(hAvgDelta2vsCentNNPionPion);
  fListHist->Add(hAvgDelta3vsCentNNPionPion);
  fListHist->Add(hAvgDelta4vsCentNNPionPion);
  fListHist->Add(hAvg3pC112vsCentNNKaonKaon);
  fListHist->Add(hAvg3pC123vsCentNNKaonKaon);
  fListHist->Add(hAvgDelta1vsCentNNKaonKaon);
  fListHist->Add(hAvgDelta2vsCentNNKaonKaon);
  fListHist->Add(hAvgDelta3vsCentNNKaonKaon);
  fListHist->Add(hAvgDelta4vsCentNNKaonKaon);
  fListHist->Add(hAvg3pC112vsCentNNProtonProton);
  fListHist->Add(hAvg3pC123vsCentNNProtonProton);
  fListHist->Add(hAvgDelta1vsCentNNProtonProton);
  fListHist->Add(hAvgDelta2vsCentNNProtonProton);
  fListHist->Add(hAvgDelta3vsCentNNProtonProton);
  fListHist->Add(hAvgDelta4vsCentNNProtonProton);
  
  hAvg3pC112vsCentNNPionCharge = new TProfile("hAvg3pC112vsCentNNPionCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNPionCharge = new TProfile("hAvg3pC123vsCentNNPionCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNPionCharge = new TProfile("hAvgDelta1vsCentNNPionCharge"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNPionCharge = new TProfile("hAvgDelta2vsCentNNPionCharge"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNPionCharge = new TProfile("hAvgDelta3vsCentNNPionCharge"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNPionCharge = new TProfile("hAvgDelta4vsCentNNPionCharge"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentNNKaonCharge = new TProfile("hAvg3pC112vsCentNNKaonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNKaonCharge = new TProfile("hAvg3pC123vsCentNNKaonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNKaonCharge = new TProfile("hAvgDelta1vsCentNNKaonCharge"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNKaonCharge = new TProfile("hAvgDelta2vsCentNNKaonCharge"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNKaonCharge = new TProfile("hAvgDelta3vsCentNNKaonCharge"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNKaonCharge = new TProfile("hAvgDelta4vsCentNNKaonCharge"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvg3pC112vsCentNNProtonCharge = new TProfile("hAvg3pC112vsCentNNProtonCharge"," <3p> vs Cent; Cent; <#gamma112>",18,0,90);
  hAvg3pC123vsCentNNProtonCharge = new TProfile("hAvg3pC123vsCentNNProtonCharge"," <3p> vs Cent; Cent; <#gamma123>",18,0,90);
  hAvgDelta1vsCentNNProtonCharge = new TProfile("hAvgDelta1vsCentNNProtonCharge"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta2vsCentNNProtonCharge = new TProfile("hAvgDelta2vsCentNNProtonCharge"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta3vsCentNNProtonCharge = new TProfile("hAvgDelta3vsCentNNProtonCharge"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);
  hAvgDelta4vsCentNNProtonCharge = new TProfile("hAvgDelta4vsCentNNProtonCharge"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",18,0,90);

  hAvg3pC112vsCentNNPionCharge->Sumw2();
  hAvg3pC123vsCentNNPionCharge->Sumw2();
  hAvgDelta1vsCentNNPionCharge->Sumw2();
  hAvgDelta2vsCentNNPionCharge->Sumw2();
  hAvgDelta3vsCentNNPionCharge->Sumw2();
  hAvgDelta4vsCentNNPionCharge->Sumw2();
  hAvg3pC112vsCentNNKaonCharge->Sumw2();
  hAvg3pC123vsCentNNKaonCharge->Sumw2();
  hAvgDelta1vsCentNNKaonCharge->Sumw2();
  hAvgDelta2vsCentNNKaonCharge->Sumw2();
  hAvgDelta3vsCentNNKaonCharge->Sumw2();
  hAvgDelta4vsCentNNKaonCharge->Sumw2();
  hAvg3pC112vsCentNNProtonCharge->Sumw2();
  hAvg3pC123vsCentNNProtonCharge->Sumw2();
  hAvgDelta1vsCentNNProtonCharge->Sumw2();
  hAvgDelta2vsCentNNProtonCharge->Sumw2();
  hAvgDelta3vsCentNNProtonCharge->Sumw2();
  hAvgDelta4vsCentNNProtonCharge->Sumw2();
  
  fListHist->Add(hAvg3pC112vsCentNNPionCharge);
  fListHist->Add(hAvg3pC123vsCentNNPionCharge);
  fListHist->Add(hAvgDelta1vsCentNNPionCharge);
  fListHist->Add(hAvgDelta2vsCentNNPionCharge);
  fListHist->Add(hAvgDelta3vsCentNNPionCharge);
  fListHist->Add(hAvgDelta4vsCentNNPionCharge);
  fListHist->Add(hAvg3pC112vsCentNNKaonCharge);
  fListHist->Add(hAvg3pC123vsCentNNKaonCharge);
  fListHist->Add(hAvgDelta1vsCentNNKaonCharge);
  fListHist->Add(hAvgDelta2vsCentNNKaonCharge);
  fListHist->Add(hAvgDelta3vsCentNNKaonCharge);
  fListHist->Add(hAvgDelta4vsCentNNKaonCharge);
  fListHist->Add(hAvg3pC112vsCentNNProtonCharge);
  fListHist->Add(hAvg3pC123vsCentNNProtonCharge);
  fListHist->Add(hAvgDelta1vsCentNNProtonCharge);
  fListHist->Add(hAvgDelta2vsCentNNProtonCharge);
  fListHist->Add(hAvgDelta3vsCentNNProtonCharge);
  fListHist->Add(hAvgDelta4vsCentNNProtonCharge);

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
  
  /// Q vec related histograms
  Double_t fCRCEtaBinEdges[fCMEnEtaBin+1] = {-0.8, -0.1, 0, 0.1, 0.8};
  
  for(Int_t c=0;c<4;c++) {
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQRe[c][h] = new TH1D(Form("fCMEQRe[%d][%d]",c,h),Form("fCMEQRe[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQRe[c][h]);
      fCMEQIm[c][h] = new TH1D(Form("fCMEQIm[%d][%d]",c,h),Form("fCMEQIm[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQIm[c][h]);
      fCMEMult[c][h] = new TH1D(Form("fCMEMult[%d][%d]",c,h),Form("fCMEMult[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMult[c][h]);
      
      fCMEQRePion[c][h] = new TH1D(Form("fCMEQRePion[%d][%d]",c,h),Form("fCMEQRePion[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQRePion[c][h]);
      fCMEQImPion[c][h] = new TH1D(Form("fCMEQImPion[%d][%d]",c,h),Form("fCMEQImPion[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImPion[c][h]);
      fCMEMultPion[c][h] = new TH1D(Form("fCMEMultPion[%d][%d]",c,h),Form("fCMEMultPion[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultPion[c][h]);
      
      fCMEQReKaon[c][h] = new TH1D(Form("fCMEQReKaon[%d][%d]",c,h),Form("fCMEQReKaon[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQReKaon[c][h]);
      fCMEQImKaon[c][h] = new TH1D(Form("fCMEQImKaon[%d][%d]",c,h),Form("fCMEQImKaon[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImKaon[c][h]);
      fCMEMultKaon[c][h] = new TH1D(Form("fCMEMultKaon[%d][%d]",c,h),Form("fCMEMultKaon[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultKaon[c][h]);
      
      fCMEQReProton[c][h] = new TH1D(Form("fCMEQReProton[%d][%d]",c,h),Form("fCMEQReProton[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQReProton[c][h]);
      fCMEQImProton[c][h] = new TH1D(Form("fCMEQImProton[%d][%d]",c,h),Form("fCMEQImProton[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImProton[c][h]);
      fCMEMultProton[c][h] = new TH1D(Form("fCMEMultProton[%d][%d]",c,h),Form("fCMEMultProton[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultProton[c][h]);
    }
  }
  

  //@shi book CME Qvector for spectator plane participant plane method
  for(Int_t c=0;c<2;c++) { // c is the index for the power of weight: weight^c*cos((h+1)*phi)
    for (Int_t h=0;h<fCRCnHar;h++) {
      fCMEQReBothCharge[c][h] = new TH1D(Form("fCMEQReBothCharge[%d][%d]",c,h),Form("fCMEQReBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQReBothCharge[c][h]);
      fCMEQImBothCharge[c][h] = new TH1D(Form("fCMEQImBothCharge[%d][%d]",c,h),Form("fCMEQImBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImBothCharge[c][h]);
      fCMEMultBothCharge[c][h] = new TH1D(Form("fCMEMultBothCharge[%d][%d]",c,h),Form("fCMEMultBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultBothCharge[c][h]);
      
      fCMEQRePionBothCharge[c][h] = new TH1D(Form("fCMEQRePionBothCharge[%d][%d]",c,h),Form("fCMEQRePionBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQRePionBothCharge[c][h]);
      fCMEQImPionBothCharge[c][h] = new TH1D(Form("fCMEQImPionBothCharge[%d][%d]",c,h),Form("fCMEQImPionBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImPionBothCharge[c][h]);
      fCMEMultPionBothCharge[c][h] = new TH1D(Form("fCMEMultPionBothCharge[%d][%d]",c,h),Form("fCMEMultPionBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultPionBothCharge[c][h]);
      
      fCMEQReKaonBothCharge[c][h] = new TH1D(Form("fCMEQReKaonBothCharge[%d][%d]",c,h),Form("fCMEQReKaonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQReKaonBothCharge[c][h]);
      fCMEQImKaonBothCharge[c][h] = new TH1D(Form("fCMEQImKaonBothCharge[%d][%d]",c,h),Form("fCMEQImKaonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImKaonBothCharge[c][h]);
      fCMEMultKaonBothCharge[c][h] = new TH1D(Form("fCMEMultKaonBothCharge[%d][%d]",c,h),Form("fCMEMultKaonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultKaonBothCharge[c][h]);
      
      fCMEQReProtonBothCharge[c][h] = new TH1D(Form("fCMEQReProtonBothCharge[%d][%d]",c,h),Form("fCMEQReProtonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQReProtonBothCharge[c][h]);
      fCMEQImProtonBothCharge[c][h] = new TH1D(Form("fCMEQImProtonBothCharge[%d][%d]",c,h),Form("fCMEQImProtonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEQImProtonBothCharge[c][h]);
      fCMEMultProtonBothCharge[c][h] = new TH1D(Form("fCMEMultProtonBothCharge[%d][%d]",c,h),Form("fCMEMultProtonBothCharge[%d][%d]",c,h),fCMEnEtaBin,fCRCEtaBinEdges);
      fTempList->Add(fCMEMultProtonBothCharge[c][h]);
    }
  }
  
  for(Int_t c=0;c<2;c++) {
	fCMEQ2Re4[c] = new TH1D(Form("fCMEQ2Re4[%d]",c),Form("fCMEQ2Re4[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Re4[c]);
    fCMEQ3Re2[c] = new TH1D(Form("fCMEQ3Re2[%d]",c),Form("fCMEQ3Re2[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Re2[c]);
    fCMEQ2Im4[c] = new TH1D(Form("fCMEQ2Im4[%d]",c),Form("fCMEQ2Im4[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Im4[c]);
    fCMEQ3Im2[c] = new TH1D(Form("fCMEQ3Im2[%d]",c),Form("fCMEQ3Im2[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Im2[c]);
    fCMEw0[c] = new TH1D(Form("fCMEw0[%d]",c),Form("fCMEw0[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw0[c]);
    fCMEw1[c] = new TH1D(Form("fCMEw1[%d]",c),Form("fCMEw1[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw1[c]);
    fCMEw2[c] = new TH1D(Form("fCMEw2[%d]",c),Form("fCMEw2[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw2[c]);
    fCMEw3[c] = new TH1D(Form("fCMEw3[%d]",c),Form("fCMEw3[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw3[c]);
    fCMEw4[c] = new TH1D(Form("fCMEw4[%d]",c),Form("fCMEw4[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw4[c]);
    
    fCMEQ2Re4Pion[c] = new TH1D(Form("fCMEQ2Re4Pion[%d]",c),Form("fCMEQ2Re4Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Re4Pion[c]);
    fCMEQ3Re2Pion[c] = new TH1D(Form("fCMEQ3Re2Pion[%d]",c),Form("fCMEQ3Re2Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Re2Pion[c]);
    fCMEQ2Im4Pion[c] = new TH1D(Form("fCMEQ2Im4Pion[%d]",c),Form("fCMEQ2Im4Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Im4Pion[c]);
    fCMEQ3Im2Pion[c] = new TH1D(Form("fCMEQ3Im2Pion[%d]",c),Form("fCMEQ3Im2Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Im2Pion[c]);
    fCMEw0Pion[c] = new TH1D(Form("fCMEw0Pion[%d]",c),Form("fCMEw0Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw0Pion[c]);
    fCMEw1Pion[c] = new TH1D(Form("fCMEw1Pion[%d]",c),Form("fCMEw1Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw1Pion[c]);
    fCMEw2Pion[c] = new TH1D(Form("fCMEw2Pion[%d]",c),Form("fCMEw2Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw2Pion[c]);
    fCMEw3Pion[c] = new TH1D(Form("fCMEw3Pion[%d]",c),Form("fCMEw3Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw3Pion[c]);
    fCMEw4Pion[c] = new TH1D(Form("fCMEw4Pion[%d]",c),Form("fCMEw4Pion[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw4Pion[c]);
    
    fCMEQ2Re4Kaon[c] = new TH1D(Form("fCMEQ2Re4Kaon[%d]",c),Form("fCMEQ2Re4Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Re4Kaon[c]);
    fCMEQ3Re2Kaon[c] = new TH1D(Form("fCMEQ3Re2Kaon[%d]",c),Form("fCMEQ3Re2Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Re2Kaon[c]);
    fCMEQ2Im4Kaon[c] = new TH1D(Form("fCMEQ2Im4Kaon[%d]",c),Form("fCMEQ2Im4Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Im4Kaon[c]);
    fCMEQ3Im2Kaon[c] = new TH1D(Form("fCMEQ3Im2Kaon[%d]",c),Form("fCMEQ3Im2Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Im2Kaon[c]);
    fCMEw0Kaon[c] = new TH1D(Form("fCMEw0Kaon[%d]",c),Form("fCMEw0Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw0Kaon[c]);
    fCMEw1Kaon[c] = new TH1D(Form("fCMEw1Kaon[%d]",c),Form("fCMEw1Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw1Kaon[c]);
    fCMEw2Kaon[c] = new TH1D(Form("fCMEw2Kaon[%d]",c),Form("fCMEw2Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw2Kaon[c]);
    fCMEw3Kaon[c] = new TH1D(Form("fCMEw3Kaon[%d]",c),Form("fCMEw3Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw3Kaon[c]);
    fCMEw4Kaon[c] = new TH1D(Form("fCMEw4Kaon[%d]",c),Form("fCMEw4Kaon[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw4Kaon[c]);
    
    fCMEQ2Re4Proton[c] = new TH1D(Form("fCMEQ2Re4Proton[%d]",c),Form("fCMEQ2Re4Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Re4Proton[c]);
    fCMEQ3Re2Proton[c] = new TH1D(Form("fCMEQ3Re2Proton[%d]",c),Form("fCMEQ3Re2Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Re2Proton[c]);
    fCMEQ2Im4Proton[c] = new TH1D(Form("fCMEQ2Im4Proton[%d]",c),Form("fCMEQ2Im4Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ2Im4Proton[c]);
    fCMEQ3Im2Proton[c] = new TH1D(Form("fCMEQ3Im2Proton[%d]",c),Form("fCMEQ3Im2Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEQ3Im2Proton[c]);
    fCMEw0Proton[c] = new TH1D(Form("fCMEw0Proton[%d]",c),Form("fCMEw0Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw0Proton[c]);
    fCMEw1Proton[c] = new TH1D(Form("fCMEw1Proton[%d]",c),Form("fCMEw1Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw1Proton[c]);
    fCMEw2Proton[c] = new TH1D(Form("fCMEw2Proton[%d]",c),Form("fCMEw2Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw2Proton[c]);
    fCMEw3Proton[c] = new TH1D(Form("fCMEw3Proton[%d]",c),Form("fCMEw3Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw3Proton[c]);
    fCMEw4Proton[c] = new TH1D(Form("fCMEw4Proton[%d]",c),Form("fCMEw4Proton[%d]",c),fCMEnEtaBin,fCRCEtaBinEdges);
    fTempList->Add(fCMEw4Proton[c]);
  }
  
}



void AliAnalysisTaskGammaDeltaPIDSaveQvec::CalculateCMESPPP()
{
  //************************************************ Get all variables ****************************
  //*********************************************** TPC part *************************************
  Int_t h = 0; //@Shi used for TPC and v0 part. For ZDCpart, it is set to 1
  Double_t e = 1E-5;
  
  Double_t uPNReTPC=0., uPNImTPC=0., uPN2ReTPC=0., uPN2ImTPC=0., uPN2Re2TPC=0., uPN2Im2TPC=0., uPNMTPC=0., uPN2MTPC=0.;
  Double_t uPNReTPCPosEta=0., uPNImTPCPosEta=0., uPN2ReTPCPosEta=0., uPN2ImTPCPosEta=0., uPN2Re2TPCPosEta=0., uPN2Im2TPCPosEta=0., uPNMTPCPosEta=0., uPN2MTPCPosEta=0.;
  Double_t uPNReTPCNegEta=0., uPNImTPCNegEta=0., uPN2ReTPCNegEta=0., uPN2ImTPCNegEta=0., uPN2Re2TPCNegEta=0., uPN2Im2TPCNegEta=0., uPNMTPCNegEta=0., uPN2MTPCNegEta=0.;
  
  Double_t uPReTPCPosEta=0., uPImTPCPosEta=0., uP2ReTPCPosEta=0., uP2ImTPCPosEta=0., uP2Re2TPCPosEta=0., uP2Im2TPCPosEta=0., uPMTPCPosEta=0., uP2MTPCPosEta=0.;
  Double_t uNReTPCPosEta=0., uNImTPCPosEta=0., uN2ReTPCPosEta=0., uN2ImTPCPosEta=0., uN2Re2TPCPosEta=0., uN2Im2TPCPosEta=0., uNMTPCPosEta=0., uN2MTPCPosEta=0.;
  
  Double_t uPReTPCNegEta=0., uPImTPCNegEta=0., uP2ReTPCNegEta=0., uP2ImTPCNegEta=0., uP2Re2TPCNegEta=0., uP2Im2TPCNegEta=0., uPMTPCNegEta=0., uP2MTPCNegEta=0.;
  Double_t uNReTPCNegEta=0., uNImTPCNegEta=0., uN2ReTPCNegEta=0., uN2ImTPCNegEta=0., uN2Re2TPCNegEta=0., uN2Im2TPCNegEta=0., uNMTPCNegEta=0., uN2MTPCNegEta=0.;

  Double_t uP4Re2TPCPosEta=0., uP4Im2TPCPosEta=0., uP2Re3TPCPosEta=0., uP2Im3TPCPosEta=0.;
  Double_t uN4Re2TPCPosEta=0., uN4Im2TPCPosEta=0., uN2Re3TPCPosEta=0., uN2Im3TPCPosEta=0.;
  
  Double_t uP4Re2TPCNegEta=0., uP4Im2TPCNegEta=0., uP2Re3TPCNegEta=0., uP2Im3TPCNegEta=0.;
  Double_t uN4Re2TPCNegEta=0., uN4Im2TPCNegEta=0., uN2Re3TPCNegEta=0., uN2Im3TPCNegEta=0.;

  Double_t uP0MTPCPosEta=0., uP3MTPCPosEta=0., uP4MTPCPosEta=0.;
  Double_t uN0MTPCPosEta=0., uN3MTPCPosEta=0., uN4MTPCPosEta=0.;
  Double_t uP0MTPCNegEta=0., uP3MTPCNegEta=0., uP4MTPCNegEta=0.;
  Double_t uN0MTPCNegEta=0., uN3MTPCNegEta=0., uN4MTPCNegEta=0.;
  
  Double_t uP4Re2TPCSubPosEta=0., uP4Im2TPCSubPosEta=0., uP2Re3TPCSubPosEta=0., uP2Im3TPCSubPosEta=0.;
  Double_t uN4Re2TPCSubPosEta=0., uN4Im2TPCSubPosEta=0., uN2Re3TPCSubPosEta=0., uN2Im3TPCSubPosEta=0.;
  
  Double_t uP4Re2TPCSubNegEta=0., uP4Im2TPCSubNegEta=0., uP2Re3TPCSubNegEta=0., uP2Im3TPCSubNegEta=0.;
  Double_t uN4Re2TPCSubNegEta=0., uN4Im2TPCSubNegEta=0., uN2Re3TPCSubNegEta=0., uN2Im3TPCSubNegEta=0.;

  Double_t uP0MTPCSubPosEta=0., uP3MTPCSubPosEta=0., uP4MTPCSubPosEta=0.;
  Double_t uN0MTPCSubPosEta=0., uN3MTPCSubPosEta=0., uN4MTPCSubPosEta=0.;
  Double_t uP0MTPCSubNegEta=0., uP3MTPCSubNegEta=0., uP4MTPCSubNegEta=0.;
  Double_t uN0MTPCSubNegEta=0., uN3MTPCSubNegEta=0., uN4MTPCSubNegEta=0.;

  Double_t uPNReTPCPion=0., uPNImTPCPion=0., uPN2ReTPCPion=0., uPN2ImTPCPion=0., uPN2Re2TPCPion=0., uPN2Im2TPCPion=0., uPNMTPCPion=0., uPN2MTPCPion=0.;
  Double_t uPNReTPCPionPosEta=0., uPNImTPCPionPosEta=0., uPN2ReTPCPionPosEta=0., uPN2ImTPCPionPosEta=0., uPN2Re2TPCPionPosEta=0., uPN2Im2TPCPionPosEta=0., uPNMTPCPionPosEta=0., uPN2MTPCPionPosEta=0.;
  Double_t uPNReTPCPionNegEta=0., uPNImTPCPionNegEta=0., uPN2ReTPCPionNegEta=0., uPN2ImTPCPionNegEta=0., uPN2Re2TPCPionNegEta=0., uPN2Im2TPCPionNegEta=0., uPNMTPCPionNegEta=0., uPN2MTPCPionNegEta=0.;
  
  Double_t uPReTPCPionPosEta=0., uPImTPCPionPosEta=0., uP2ReTPCPionPosEta=0., uP2ImTPCPionPosEta=0., uP2Re2TPCPionPosEta=0., uP2Im2TPCPionPosEta=0., uPMTPCPionPosEta=0., uP2MTPCPionPosEta=0.;
  Double_t uNReTPCPionPosEta=0., uNImTPCPionPosEta=0., uN2ReTPCPionPosEta=0., uN2ImTPCPionPosEta=0., uN2Re2TPCPionPosEta=0., uN2Im2TPCPionPosEta=0., uNMTPCPionPosEta=0., uN2MTPCPionPosEta=0.;
  
  Double_t uPReTPCPionNegEta=0., uPImTPCPionNegEta=0., uP2ReTPCPionNegEta=0., uP2ImTPCPionNegEta=0., uP2Re2TPCPionNegEta=0., uP2Im2TPCPionNegEta=0., uPMTPCPionNegEta=0., uP2MTPCPionNegEta=0.;
  Double_t uNReTPCPionNegEta=0., uNImTPCPionNegEta=0., uN2ReTPCPionNegEta=0., uN2ImTPCPionNegEta=0., uN2Re2TPCPionNegEta=0., uN2Im2TPCPionNegEta=0., uNMTPCPionNegEta=0., uN2MTPCPionNegEta=0.;
  
  Double_t uP4Re2TPCPionPosEta=0., uP4Im2TPCPionPosEta=0., uP2Re3TPCPionPosEta=0., uP2Im3TPCPionPosEta=0.;
  Double_t uN4Re2TPCPionPosEta=0., uN4Im2TPCPionPosEta=0., uN2Re3TPCPionPosEta=0., uN2Im3TPCPionPosEta=0.;
  
  Double_t uP4Re2TPCPionNegEta=0., uP4Im2TPCPionNegEta=0., uP2Re3TPCPionNegEta=0., uP2Im3TPCPionNegEta=0.;
  Double_t uN4Re2TPCPionNegEta=0., uN4Im2TPCPionNegEta=0., uN2Re3TPCPionNegEta=0., uN2Im3TPCPionNegEta=0.;

  Double_t uP0MTPCPionPosEta=0., uP3MTPCPionPosEta=0., uP4MTPCPionPosEta=0.;
  Double_t uN0MTPCPionPosEta=0., uN3MTPCPionPosEta=0., uN4MTPCPionPosEta=0.;
  Double_t uP0MTPCPionNegEta=0., uP3MTPCPionNegEta=0., uP4MTPCPionNegEta=0.;
  Double_t uN0MTPCPionNegEta=0., uN3MTPCPionNegEta=0., uN4MTPCPionNegEta=0.;
  
  Double_t uP4Re2TPCPionSubPosEta=0., uP4Im2TPCPionSubPosEta=0., uP2Re3TPCPionSubPosEta=0., uP2Im3TPCPionSubPosEta=0.;
  Double_t uN4Re2TPCPionSubPosEta=0., uN4Im2TPCPionSubPosEta=0., uN2Re3TPCPionSubPosEta=0., uN2Im3TPCPionSubPosEta=0.;
  
  Double_t uP4Re2TPCPionSubNegEta=0., uP4Im2TPCPionSubNegEta=0., uP2Re3TPCPionSubNegEta=0., uP2Im3TPCPionSubNegEta=0.;
  Double_t uN4Re2TPCPionSubNegEta=0., uN4Im2TPCPionSubNegEta=0., uN2Re3TPCPionSubNegEta=0., uN2Im3TPCPionSubNegEta=0.;

  Double_t uP0MTPCPionSubPosEta=0., uP3MTPCPionSubPosEta=0., uP4MTPCPionSubPosEta=0.;
  Double_t uN0MTPCPionSubPosEta=0., uN3MTPCPionSubPosEta=0., uN4MTPCPionSubPosEta=0.;
  Double_t uP0MTPCPionSubNegEta=0., uP3MTPCPionSubNegEta=0., uP4MTPCPionSubNegEta=0.;
  Double_t uN0MTPCPionSubNegEta=0., uN3MTPCPionSubNegEta=0., uN4MTPCPionSubNegEta=0.;
  
  Double_t uPNReTPCKaon=0., uPNImTPCKaon=0., uPN2ReTPCKaon=0., uPN2ImTPCKaon=0., uPN2Re2TPCKaon=0., uPN2Im2TPCKaon=0., uPNMTPCKaon=0., uPN2MTPCKaon=0.;
  Double_t uPNReTPCKaonPosEta=0., uPNImTPCKaonPosEta=0., uPN2ReTPCKaonPosEta=0., uPN2ImTPCKaonPosEta=0., uPN2Re2TPCKaonPosEta=0., uPN2Im2TPCKaonPosEta=0., uPNMTPCKaonPosEta=0., uPN2MTPCKaonPosEta=0.;
  Double_t uPNReTPCKaonNegEta=0., uPNImTPCKaonNegEta=0., uPN2ReTPCKaonNegEta=0., uPN2ImTPCKaonNegEta=0., uPN2Re2TPCKaonNegEta=0., uPN2Im2TPCKaonNegEta=0., uPNMTPCKaonNegEta=0., uPN2MTPCKaonNegEta=0.;
  
  Double_t uPReTPCKaonPosEta=0., uPImTPCKaonPosEta=0., uP2ReTPCKaonPosEta=0., uP2ImTPCKaonPosEta=0., uP2Re2TPCKaonPosEta=0., uP2Im2TPCKaonPosEta=0., uPMTPCKaonPosEta=0., uP2MTPCKaonPosEta=0.;
  Double_t uNReTPCKaonPosEta=0., uNImTPCKaonPosEta=0., uN2ReTPCKaonPosEta=0., uN2ImTPCKaonPosEta=0., uN2Re2TPCKaonPosEta=0., uN2Im2TPCKaonPosEta=0., uNMTPCKaonPosEta=0., uN2MTPCKaonPosEta=0.;
  
  Double_t uPReTPCKaonNegEta=0., uPImTPCKaonNegEta=0., uP2ReTPCKaonNegEta=0., uP2ImTPCKaonNegEta=0., uP2Re2TPCKaonNegEta=0., uP2Im2TPCKaonNegEta=0., uPMTPCKaonNegEta=0., uP2MTPCKaonNegEta=0.;
  Double_t uNReTPCKaonNegEta=0., uNImTPCKaonNegEta=0., uN2ReTPCKaonNegEta=0., uN2ImTPCKaonNegEta=0., uN2Re2TPCKaonNegEta=0., uN2Im2TPCKaonNegEta=0., uNMTPCKaonNegEta=0., uN2MTPCKaonNegEta=0.;
  
  Double_t uP4Re2TPCKaonPosEta=0., uP4Im2TPCKaonPosEta=0., uP2Re3TPCKaonPosEta=0., uP2Im3TPCKaonPosEta=0.;
  Double_t uN4Re2TPCKaonPosEta=0., uN4Im2TPCKaonPosEta=0., uN2Re3TPCKaonPosEta=0., uN2Im3TPCKaonPosEta=0.;
  
  Double_t uP4Re2TPCKaonNegEta=0., uP4Im2TPCKaonNegEta=0., uP2Re3TPCKaonNegEta=0., uP2Im3TPCKaonNegEta=0.;
  Double_t uN4Re2TPCKaonNegEta=0., uN4Im2TPCKaonNegEta=0., uN2Re3TPCKaonNegEta=0., uN2Im3TPCKaonNegEta=0.;

  Double_t uP0MTPCKaonPosEta=0., uP3MTPCKaonPosEta=0., uP4MTPCKaonPosEta=0.;
  Double_t uN0MTPCKaonPosEta=0., uN3MTPCKaonPosEta=0., uN4MTPCKaonPosEta=0.;
  Double_t uP0MTPCKaonNegEta=0., uP3MTPCKaonNegEta=0., uP4MTPCKaonNegEta=0.;
  Double_t uN0MTPCKaonNegEta=0., uN3MTPCKaonNegEta=0., uN4MTPCKaonNegEta=0.;
  
  Double_t uP4Re2TPCKaonSubPosEta=0., uP4Im2TPCKaonSubPosEta=0., uP2Re3TPCKaonSubPosEta=0., uP2Im3TPCKaonSubPosEta=0.;
  Double_t uN4Re2TPCKaonSubPosEta=0., uN4Im2TPCKaonSubPosEta=0., uN2Re3TPCKaonSubPosEta=0., uN2Im3TPCKaonSubPosEta=0.;
  
  Double_t uP4Re2TPCKaonSubNegEta=0., uP4Im2TPCKaonSubNegEta=0., uP2Re3TPCKaonSubNegEta=0., uP2Im3TPCKaonSubNegEta=0.;
  Double_t uN4Re2TPCKaonSubNegEta=0., uN4Im2TPCKaonSubNegEta=0., uN2Re3TPCKaonSubNegEta=0., uN2Im3TPCKaonSubNegEta=0.;

  Double_t uP0MTPCKaonSubPosEta=0., uP3MTPCKaonSubPosEta=0., uP4MTPCKaonSubPosEta=0.;
  Double_t uN0MTPCKaonSubPosEta=0., uN3MTPCKaonSubPosEta=0., uN4MTPCKaonSubPosEta=0.;
  Double_t uP0MTPCKaonSubNegEta=0., uP3MTPCKaonSubNegEta=0., uP4MTPCKaonSubNegEta=0.;
  Double_t uN0MTPCKaonSubNegEta=0., uN3MTPCKaonSubNegEta=0., uN4MTPCKaonSubNegEta=0.;
  
  Double_t uPNReTPCProton=0., uPNImTPCProton=0., uPN2ReTPCProton=0., uPN2ImTPCProton=0., uPN2Re2TPCProton=0., uPN2Im2TPCProton=0., uPNMTPCProton=0., uPN2MTPCProton=0.;
  Double_t uPNReTPCProtonPosEta=0., uPNImTPCProtonPosEta=0., uPN2ReTPCProtonPosEta=0., uPN2ImTPCProtonPosEta=0., uPN2Re2TPCProtonPosEta=0., uPN2Im2TPCProtonPosEta=0., uPNMTPCProtonPosEta=0., uPN2MTPCProtonPosEta=0.;
  Double_t uPNReTPCProtonNegEta=0., uPNImTPCProtonNegEta=0., uPN2ReTPCProtonNegEta=0., uPN2ImTPCProtonNegEta=0., uPN2Re2TPCProtonNegEta=0., uPN2Im2TPCProtonNegEta=0., uPNMTPCProtonNegEta=0., uPN2MTPCProtonNegEta=0.;
  
  Double_t uPReTPCProtonPosEta=0., uPImTPCProtonPosEta=0., uP2ReTPCProtonPosEta=0., uP2ImTPCProtonPosEta=0., uP2Re2TPCProtonPosEta=0., uP2Im2TPCProtonPosEta=0., uPMTPCProtonPosEta=0., uP2MTPCProtonPosEta=0.;
  Double_t uNReTPCProtonPosEta=0., uNImTPCProtonPosEta=0., uN2ReTPCProtonPosEta=0., uN2ImTPCProtonPosEta=0., uN2Re2TPCProtonPosEta=0., uN2Im2TPCProtonPosEta=0., uNMTPCProtonPosEta=0., uN2MTPCProtonPosEta=0.;
  
  Double_t uPReTPCProtonNegEta=0., uPImTPCProtonNegEta=0., uP2ReTPCProtonNegEta=0., uP2ImTPCProtonNegEta=0., uP2Re2TPCProtonNegEta=0., uP2Im2TPCProtonNegEta=0., uPMTPCProtonNegEta=0., uP2MTPCProtonNegEta=0.;
  Double_t uNReTPCProtonNegEta=0., uNImTPCProtonNegEta=0., uN2ReTPCProtonNegEta=0., uN2ImTPCProtonNegEta=0., uN2Re2TPCProtonNegEta=0., uN2Im2TPCProtonNegEta=0., uNMTPCProtonNegEta=0., uN2MTPCProtonNegEta=0.; 
  
  Double_t uP4Re2TPCProtonPosEta=0., uP4Im2TPCProtonPosEta=0., uP2Re3TPCProtonPosEta=0., uP2Im3TPCProtonPosEta=0.;
  Double_t uN4Re2TPCProtonPosEta=0., uN4Im2TPCProtonPosEta=0., uN2Re3TPCProtonPosEta=0., uN2Im3TPCProtonPosEta=0.;
  
  Double_t uP4Re2TPCProtonNegEta=0., uP4Im2TPCProtonNegEta=0., uP2Re3TPCProtonNegEta=0., uP2Im3TPCProtonNegEta=0.;
  Double_t uN4Re2TPCProtonNegEta=0., uN4Im2TPCProtonNegEta=0., uN2Re3TPCProtonNegEta=0., uN2Im3TPCProtonNegEta=0.;

  Double_t uP0MTPCProtonPosEta=0., uP3MTPCProtonPosEta=0., uP4MTPCProtonPosEta=0.;
  Double_t uN0MTPCProtonPosEta=0., uN3MTPCProtonPosEta=0., uN4MTPCProtonPosEta=0.;
  Double_t uP0MTPCProtonNegEta=0., uP3MTPCProtonNegEta=0., uP4MTPCProtonNegEta=0.;
  Double_t uN0MTPCProtonNegEta=0., uN3MTPCProtonNegEta=0., uN4MTPCProtonNegEta=0.;
  
  Double_t uP4Re2TPCProtonSubPosEta=0., uP4Im2TPCProtonSubPosEta=0., uP2Re3TPCProtonSubPosEta=0., uP2Im3TPCProtonSubPosEta=0.;
  Double_t uN4Re2TPCProtonSubPosEta=0., uN4Im2TPCProtonSubPosEta=0., uN2Re3TPCProtonSubPosEta=0., uN2Im3TPCProtonSubPosEta=0.;
  
  Double_t uP4Re2TPCProtonSubNegEta=0., uP4Im2TPCProtonSubNegEta=0., uP2Re3TPCProtonSubNegEta=0., uP2Im3TPCProtonSubNegEta=0.;
  Double_t uN4Re2TPCProtonSubNegEta=0., uN4Im2TPCProtonSubNegEta=0., uN2Re3TPCProtonSubNegEta=0., uN2Im3TPCProtonSubNegEta=0.;

  Double_t uP0MTPCProtonSubPosEta=0., uP3MTPCProtonSubPosEta=0., uP4MTPCProtonSubPosEta=0.;
  Double_t uN0MTPCProtonSubPosEta=0., uN3MTPCProtonSubPosEta=0., uN4MTPCProtonSubPosEta=0.;
  Double_t uP0MTPCProtonSubNegEta=0., uP3MTPCProtonSubNegEta=0., uP4MTPCProtonSubNegEta=0.;
  Double_t uN0MTPCProtonSubNegEta=0., uN3MTPCProtonSubNegEta=0., uN4MTPCProtonSubNegEta=0.;
  
  Double_t uPNReTPCSubPosEta = 0, uPNImTPCSubPosEta = 0, uPN2ReTPCSubPosEta = 0, uPN2ImTPCSubPosEta = 0, uPN2Re2TPCSubPosEta = 0, uPN2Im2TPCSubPosEta = 0, uPNMTPCSubPosEta = 0, uPN2MTPCSubPosEta = 0, uPReTPCSubPosEta = 0, uPImTPCSubPosEta = 0, uP2ReTPCSubPosEta = 0, uP2ImTPCSubPosEta = 0, uP2Re2TPCSubPosEta = 0, uP2Im2TPCSubPosEta = 0, uPMTPCSubPosEta = 0, uP2MTPCSubPosEta = 0, uNReTPCSubPosEta = 0, uNImTPCSubPosEta = 0, uN2ReTPCSubPosEta = 0, uN2ImTPCSubPosEta = 0, uN2Re2TPCSubPosEta = 0, uN2Im2TPCSubPosEta = 0, uNMTPCSubPosEta = 0, uN2MTPCSubPosEta = 0, uPNReTPCPionSubPosEta = 0, uPNImTPCPionSubPosEta = 0, uPN2ReTPCPionSubPosEta = 0, uPN2ImTPCPionSubPosEta = 0, uPN2Re2TPCPionSubPosEta = 0, uPN2Im2TPCPionSubPosEta = 0, uPNMTPCPionSubPosEta = 0, uPN2MTPCPionSubPosEta = 0, uPReTPCPionSubPosEta = 0, uPImTPCPionSubPosEta = 0, uP2ReTPCPionSubPosEta = 0, uP2ImTPCPionSubPosEta = 0, uP2Re2TPCPionSubPosEta = 0, uP2Im2TPCPionSubPosEta = 0, uPMTPCPionSubPosEta = 0, uP2MTPCPionSubPosEta = 0, uNReTPCPionSubPosEta = 0, uNImTPCPionSubPosEta = 0, uN2ReTPCPionSubPosEta = 0, uN2ImTPCPionSubPosEta = 0, uN2Re2TPCPionSubPosEta = 0, uN2Im2TPCPionSubPosEta = 0, uNMTPCPionSubPosEta = 0, uN2MTPCPionSubPosEta = 0, uPNReTPCKaonSubPosEta = 0, uPNImTPCKaonSubPosEta = 0, uPN2ReTPCKaonSubPosEta = 0, uPN2ImTPCKaonSubPosEta = 0, uPN2Re2TPCKaonSubPosEta = 0, uPN2Im2TPCKaonSubPosEta = 0, uPNMTPCKaonSubPosEta = 0, uPN2MTPCKaonSubPosEta = 0, uPReTPCKaonSubPosEta = 0, uPImTPCKaonSubPosEta = 0, uP2ReTPCKaonSubPosEta = 0, uP2ImTPCKaonSubPosEta = 0, uP2Re2TPCKaonSubPosEta = 0, uP2Im2TPCKaonSubPosEta = 0, uPMTPCKaonSubPosEta = 0, uP2MTPCKaonSubPosEta = 0, uNReTPCKaonSubPosEta = 0, uNImTPCKaonSubPosEta = 0, uN2ReTPCKaonSubPosEta = 0, uN2ImTPCKaonSubPosEta = 0, uN2Re2TPCKaonSubPosEta = 0, uN2Im2TPCKaonSubPosEta = 0, uNMTPCKaonSubPosEta = 0, uN2MTPCKaonSubPosEta = 0, uPNReTPCProtonSubPosEta = 0, uPNImTPCProtonSubPosEta = 0, uPN2ReTPCProtonSubPosEta = 0, uPN2ImTPCProtonSubPosEta = 0, uPN2Re2TPCProtonSubPosEta = 0, uPN2Im2TPCProtonSubPosEta = 0, uPNMTPCProtonSubPosEta = 0, uPN2MTPCProtonSubPosEta = 0, uPReTPCProtonSubPosEta = 0, uPImTPCProtonSubPosEta = 0, uP2ReTPCProtonSubPosEta = 0, uP2ImTPCProtonSubPosEta = 0, uP2Re2TPCProtonSubPosEta = 0, uP2Im2TPCProtonSubPosEta = 0, uPMTPCProtonSubPosEta = 0, uP2MTPCProtonSubPosEta = 0, uNReTPCProtonSubPosEta = 0, uNImTPCProtonSubPosEta = 0, uN2ReTPCProtonSubPosEta = 0, uN2ImTPCProtonSubPosEta = 0, uN2Re2TPCProtonSubPosEta = 0, uN2Im2TPCProtonSubPosEta = 0, uNMTPCProtonSubPosEta = 0, uN2MTPCProtonSubPosEta = 0;
  
  Double_t uPNReTPCSubNegEta = 0, uPNImTPCSubNegEta = 0, uPN2ReTPCSubNegEta = 0, uPN2ImTPCSubNegEta = 0, uPN2Re2TPCSubNegEta = 0, uPN2Im2TPCSubNegEta = 0, uPNMTPCSubNegEta = 0, uPN2MTPCSubNegEta = 0, uPReTPCSubNegEta = 0, uPImTPCSubNegEta = 0, uP2ReTPCSubNegEta = 0, uP2ImTPCSubNegEta = 0, uP2Re2TPCSubNegEta = 0, uP2Im2TPCSubNegEta = 0, uPMTPCSubNegEta = 0, uP2MTPCSubNegEta = 0, uNReTPCSubNegEta = 0, uNImTPCSubNegEta = 0, uN2ReTPCSubNegEta = 0, uN2ImTPCSubNegEta = 0, uN2Re2TPCSubNegEta = 0, uN2Im2TPCSubNegEta = 0, uNMTPCSubNegEta = 0, uN2MTPCSubNegEta = 0, uPNReTPCPionSubNegEta = 0, uPNImTPCPionSubNegEta = 0, uPN2ReTPCPionSubNegEta = 0, uPN2ImTPCPionSubNegEta = 0, uPN2Re2TPCPionSubNegEta = 0, uPN2Im2TPCPionSubNegEta = 0, uPNMTPCPionSubNegEta = 0, uPN2MTPCPionSubNegEta = 0, uPReTPCPionSubNegEta = 0, uPImTPCPionSubNegEta = 0, uP2ReTPCPionSubNegEta = 0, uP2ImTPCPionSubNegEta = 0, uP2Re2TPCPionSubNegEta = 0, uP2Im2TPCPionSubNegEta = 0, uPMTPCPionSubNegEta = 0, uP2MTPCPionSubNegEta = 0, uNReTPCPionSubNegEta = 0, uNImTPCPionSubNegEta = 0, uN2ReTPCPionSubNegEta = 0, uN2ImTPCPionSubNegEta = 0, uN2Re2TPCPionSubNegEta = 0, uN2Im2TPCPionSubNegEta = 0, uNMTPCPionSubNegEta = 0, uN2MTPCPionSubNegEta = 0, uPNReTPCKaonSubNegEta = 0, uPNImTPCKaonSubNegEta = 0, uPN2ReTPCKaonSubNegEta = 0, uPN2ImTPCKaonSubNegEta = 0, uPN2Re2TPCKaonSubNegEta = 0, uPN2Im2TPCKaonSubNegEta = 0, uPNMTPCKaonSubNegEta = 0, uPN2MTPCKaonSubNegEta = 0, uPReTPCKaonSubNegEta = 0, uPImTPCKaonSubNegEta = 0, uP2ReTPCKaonSubNegEta = 0, uP2ImTPCKaonSubNegEta = 0, uP2Re2TPCKaonSubNegEta = 0, uP2Im2TPCKaonSubNegEta = 0, uPMTPCKaonSubNegEta = 0, uP2MTPCKaonSubNegEta = 0, uNReTPCKaonSubNegEta = 0, uNImTPCKaonSubNegEta = 0, uN2ReTPCKaonSubNegEta = 0, uN2ImTPCKaonSubNegEta = 0, uN2Re2TPCKaonSubNegEta = 0, uN2Im2TPCKaonSubNegEta = 0, uNMTPCKaonSubNegEta = 0, uN2MTPCKaonSubNegEta = 0, uPNReTPCProtonSubNegEta = 0, uPNImTPCProtonSubNegEta = 0, uPN2ReTPCProtonSubNegEta = 0, uPN2ImTPCProtonSubNegEta = 0, uPN2Re2TPCProtonSubNegEta = 0, uPN2Im2TPCProtonSubNegEta = 0, uPNMTPCProtonSubNegEta = 0, uPN2MTPCProtonSubNegEta = 0, uPReTPCProtonSubNegEta = 0, uPImTPCProtonSubNegEta = 0, uP2ReTPCProtonSubNegEta = 0, uP2ImTPCProtonSubNegEta = 0, uP2Re2TPCProtonSubNegEta = 0, uP2Im2TPCProtonSubNegEta = 0, uPMTPCProtonSubNegEta = 0, uP2MTPCProtonSubNegEta = 0, uNReTPCProtonSubNegEta = 0, uNImTPCProtonSubNegEta = 0, uN2ReTPCProtonSubNegEta = 0, uN2ImTPCProtonSubNegEta = 0, uN2Re2TPCProtonSubNegEta = 0, uN2Im2TPCProtonSubNegEta = 0, uNMTPCProtonSubNegEta = 0, uN2MTPCProtonSubNegEta = 0;
  
  
  
  for(Int_t EBin=1; EBin<=fCMEQRe[0][0]->GetNbinsX(); EBin++) {
    // both charge all region [0][h]: first index is the power of weight, h is cos((h+1)phi)
    uPNReTPC += fCMEQReBothCharge[0][h]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(dPhi)
    uPNImTPC += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
    uPN2ReTPC += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(2dPhi)
    uPN2ImTPC += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
    uPN2Re2TPC += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin); // pow(SpecWeig*wPhiEta,2)*TMath::Cos(2dPhi)
    uPN2Im2TPC += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
    uPNMTPC += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
    uPN2MTPC += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
    
    // PID all region [0][h]: first index is the power of weight, h is cos((h+1)phi)
    uPNReTPCPion += fCMEQRePionBothCharge[0][h]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(dPhi)
    uPNImTPCPion += fCMEQImPionBothCharge[0][h]->GetBinContent(EBin);
    uPN2ReTPCPion += fCMEQRePionBothCharge[0][h+1]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(2dPhi)
    uPN2ImTPCPion += fCMEQImPionBothCharge[0][h+1]->GetBinContent(EBin);
    uPN2Re2TPCPion += fCMEQRePionBothCharge[1][h+1]->GetBinContent(EBin); // pow(SpecWeig*wPhiEta,2)*TMath::Cos(2dPhi)
    uPN2Im2TPCPion += fCMEQImPionBothCharge[1][h+1]->GetBinContent(EBin);
    uPNMTPCPion += fCMEMultPionBothCharge[0][h]->GetBinContent(EBin);
    uPN2MTPCPion += fCMEMultPionBothCharge[1][h]->GetBinContent(EBin);
    
    uPNReTPCKaon += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(dPhi)
    uPNImTPCKaon += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
    uPN2ReTPCKaon += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(2dPhi)
    uPN2ImTPCKaon += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
    uPN2Re2TPCKaon += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin); // pow(SpecWeig*wPhiEta,2)*TMath::Cos(2dPhi)
    uPN2Im2TPCKaon += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
    uPNMTPCKaon += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
    uPN2MTPCKaon += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
    
    uPNReTPCKaon += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(dPhi)
    uPNImTPCKaon += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
    uPN2ReTPCKaon += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin); // SpecWeig*wPhiEta*TMath::Cos(2dPhi)
    uPN2ImTPCKaon += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
    uPN2Re2TPCKaon += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin); // pow(SpecWeig*wPhiEta,2)*TMath::Cos(2dPhi)
    uPN2Im2TPCKaon += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
    uPNMTPCKaon += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
    uPN2MTPCKaon += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
    
    if (EBin == fCMEQRe[0][h]->FindBin(0.4)) { // positive eta region
		// both charge pos eta region
		uPNReTPCPosEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPosEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPosEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPosEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPosEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPosEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPosEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPosEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCPosEta += fCMEQRe[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPosEta += fCMEQIm[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPosEta += fCMEQRe[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPosEta += fCMEQIm[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPosEta += fCMEQRe[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPosEta += fCMEQIm[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPosEta += fCMEMult[0][h]->GetBinContent(EBin); // w
		uP2MTPCPosEta += fCMEMult[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCPosEta += fCMEQ2Re4[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCPosEta += fCMEQ2Im4[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCPosEta += fCMEQ3Re2[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCPosEta += fCMEQ3Im2[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCPosEta += fCMEw0[0]->GetBinContent(EBin); // w^0
		uP3MTPCPosEta += fCMEw3[0]->GetBinContent(EBin); // w^3
		uP4MTPCPosEta += fCMEw4[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCPosEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCPosEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCPosEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPosEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPosEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPosEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCPosEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCPosEta += fCMEMult[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCPosEta += fCMEQ2Re4[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCPosEta += fCMEQ2Im4[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCPosEta += fCMEQ3Re2[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCPosEta += fCMEQ3Im2[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCPosEta += fCMEw0[1]->GetBinContent(EBin); // w^0
		uN3MTPCPosEta += fCMEw3[1]->GetBinContent(EBin); // w^3
		uN4MTPCPosEta += fCMEw4[1]->GetBinContent(EBin); // w^4
		
		// ========================== Pion =======================
		// both charge pos eta region
		uPNReTPCPionPosEta += fCMEQRePionBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPionPosEta += fCMEQImPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPionPosEta += fCMEQRePionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPionPosEta += fCMEQImPionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPionPosEta += fCMEQRePionBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPionPosEta += fCMEQImPionBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPionPosEta += fCMEMultPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPionPosEta += fCMEMultPionBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCPionPosEta += fCMEQRePion[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPionPosEta += fCMEQImPion[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPionPosEta += fCMEQRePion[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPionPosEta += fCMEQImPion[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPionPosEta += fCMEQRePion[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPionPosEta += fCMEQImPion[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPionPosEta += fCMEMultPion[0][h]->GetBinContent(EBin); // w
		uP2MTPCPionPosEta += fCMEMultPion[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCPionPosEta += fCMEQ2Re4Pion[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCPionPosEta += fCMEQ2Im4Pion[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCPionPosEta += fCMEQ3Re2Pion[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCPionPosEta += fCMEQ3Im2Pion[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCPionPosEta += fCMEw0Pion[0]->GetBinContent(EBin); // w^0
		uP3MTPCPionPosEta += fCMEw3Pion[0]->GetBinContent(EBin); // w^3
		uP4MTPCPionPosEta += fCMEw4Pion[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCPionPosEta += fCMEQRePion[1][h]->GetBinContent(EBin);
		uNImTPCPionPosEta += fCMEQImPion[1][h]->GetBinContent(EBin);
		uN2ReTPCPionPosEta += fCMEQRePion[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPionPosEta += fCMEQImPion[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPionPosEta += fCMEQRePion[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPionPosEta += fCMEQImPion[3][h+1]->GetBinContent(EBin);
		uNMTPCPionPosEta += fCMEMultPion[1][h]->GetBinContent(EBin);
		uN2MTPCPionPosEta += fCMEMultPion[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCPionPosEta += fCMEQ2Re4Pion[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCPionPosEta += fCMEQ2Im4Pion[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCPionPosEta += fCMEQ3Re2Pion[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCPionPosEta += fCMEQ3Im2Pion[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCPionPosEta += fCMEw0Pion[1]->GetBinContent(EBin); // w^0
		uN3MTPCPionPosEta += fCMEw3Pion[1]->GetBinContent(EBin); // w^3
		uN4MTPCPionPosEta += fCMEw4Pion[1]->GetBinContent(EBin); // w^4
		
		// ========================== Kaon =======================
		// both charge pos eta region
		uPNReTPCKaonPosEta += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCKaonPosEta += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCKaonPosEta += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCKaonPosEta += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCKaonPosEta += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCKaonPosEta += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCKaonPosEta += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCKaonPosEta += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCKaonPosEta += fCMEQReKaon[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCKaonPosEta += fCMEQImKaon[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCKaonPosEta += fCMEQReKaon[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCKaonPosEta += fCMEQImKaon[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCKaonPosEta += fCMEQReKaon[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCKaonPosEta += fCMEQImKaon[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCKaonPosEta += fCMEMultKaon[0][h]->GetBinContent(EBin); // w
		uP2MTPCKaonPosEta += fCMEMultKaon[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCKaonPosEta += fCMEQ2Re4Kaon[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCKaonPosEta += fCMEQ2Im4Kaon[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCKaonPosEta += fCMEQ3Re2Kaon[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCKaonPosEta += fCMEQ3Im2Kaon[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCKaonPosEta += fCMEw0Kaon[0]->GetBinContent(EBin); // w^0
		uP3MTPCKaonPosEta += fCMEw3Kaon[0]->GetBinContent(EBin); // w^3
		uP4MTPCKaonPosEta += fCMEw4Kaon[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCKaonPosEta += fCMEQReKaon[1][h]->GetBinContent(EBin);
		uNImTPCKaonPosEta += fCMEQImKaon[1][h]->GetBinContent(EBin);
		uN2ReTPCKaonPosEta += fCMEQReKaon[1][h+1]->GetBinContent(EBin);
		uN2ImTPCKaonPosEta += fCMEQImKaon[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCKaonPosEta += fCMEQReKaon[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCKaonPosEta += fCMEQImKaon[3][h+1]->GetBinContent(EBin);
		uNMTPCKaonPosEta += fCMEMultKaon[1][h]->GetBinContent(EBin);
		uN2MTPCKaonPosEta += fCMEMultKaon[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCKaonPosEta += fCMEQ2Re4Kaon[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCKaonPosEta += fCMEQ2Im4Kaon[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCKaonPosEta += fCMEQ3Re2Kaon[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCKaonPosEta += fCMEQ3Im2Kaon[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCKaonPosEta += fCMEw0Kaon[1]->GetBinContent(EBin); // w^0
		uN3MTPCKaonPosEta += fCMEw3Kaon[1]->GetBinContent(EBin); // w^3
		uN4MTPCKaonPosEta += fCMEw4Kaon[1]->GetBinContent(EBin); // w^4
		
		// ========================== Proton =======================
		// both charge pos eta region
		uPNReTPCProtonPosEta += fCMEQReProtonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCProtonPosEta += fCMEQImProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCProtonPosEta += fCMEQReProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCProtonPosEta += fCMEQImProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCProtonPosEta += fCMEQReProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCProtonPosEta += fCMEQImProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCProtonPosEta += fCMEMultProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCProtonPosEta += fCMEMultProtonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCProtonPosEta += fCMEQReProton[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCProtonPosEta += fCMEQImProton[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCProtonPosEta += fCMEQReProton[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCProtonPosEta += fCMEQImProton[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCProtonPosEta += fCMEQReProton[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCProtonPosEta += fCMEQImProton[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCProtonPosEta += fCMEMultProton[0][h]->GetBinContent(EBin); // w
		uP2MTPCProtonPosEta += fCMEMultProton[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCProtonPosEta += fCMEQ2Re4Proton[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCProtonPosEta += fCMEQ2Im4Proton[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCProtonPosEta += fCMEQ3Re2Proton[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCProtonPosEta += fCMEQ3Im2Proton[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCProtonPosEta += fCMEw0Proton[0]->GetBinContent(EBin); // w^0
		uP3MTPCProtonPosEta += fCMEw3Proton[0]->GetBinContent(EBin); // w^3
		uP4MTPCProtonPosEta += fCMEw4Proton[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCProtonPosEta += fCMEQReProton[1][h]->GetBinContent(EBin);
		uNImTPCProtonPosEta += fCMEQImProton[1][h]->GetBinContent(EBin);
		uN2ReTPCProtonPosEta += fCMEQReProton[1][h+1]->GetBinContent(EBin);
		uN2ImTPCProtonPosEta += fCMEQImProton[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCProtonPosEta += fCMEQReProton[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCProtonPosEta += fCMEQImProton[3][h+1]->GetBinContent(EBin);
		uNMTPCProtonPosEta += fCMEMultProton[1][h]->GetBinContent(EBin);
		uN2MTPCProtonPosEta += fCMEMultProton[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCProtonPosEta += fCMEQ2Re4Proton[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCProtonPosEta += fCMEQ2Im4Proton[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCProtonPosEta += fCMEQ3Re2Proton[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCProtonPosEta += fCMEQ3Im2Proton[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCProtonPosEta += fCMEw0Proton[1]->GetBinContent(EBin); // w^0
		uN3MTPCProtonPosEta += fCMEw3Proton[1]->GetBinContent(EBin); // w^3
		uN4MTPCProtonPosEta += fCMEw4Proton[1]->GetBinContent(EBin); // w^4
	} else if (EBin == fCMEQRe[0][h]->FindBin(-0.4)) { // negative eta region
		// both charge neg eta region
		uPNReTPCNegEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCNegEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCNegEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCNegEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCNegEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCNegEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCNegEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCNegEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// pegitive charge negitive eta region
		uPReTPCNegEta += fCMEQRe[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCNegEta += fCMEQIm[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCNegEta += fCMEQRe[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCNegEta += fCMEQIm[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCNegEta += fCMEQRe[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCNegEta += fCMEQIm[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCNegEta += fCMEMult[0][h]->GetBinContent(EBin); // w
		uP2MTPCNegEta += fCMEMult[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCNegEta += fCMEQ2Re4[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCNegEta += fCMEQ2Im4[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCNegEta += fCMEQ3Re2[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCNegEta += fCMEQ3Im2[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCNegEta += fCMEw0[0]->GetBinContent(EBin); // w^0
		uP3MTPCNegEta += fCMEw3[0]->GetBinContent(EBin); // w^3
		uP4MTPCNegEta += fCMEw4[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCNegEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCNegEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCNegEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCNegEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCNegEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCNegEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCNegEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCNegEta += fCMEMult[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCNegEta += fCMEQ2Re4[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCNegEta += fCMEQ2Im4[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCNegEta += fCMEQ3Re2[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCNegEta += fCMEQ3Im2[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCNegEta += fCMEw0[1]->GetBinContent(EBin); // w^0
		uN3MTPCNegEta += fCMEw3[1]->GetBinContent(EBin); // w^3
		uN4MTPCNegEta += fCMEw4[1]->GetBinContent(EBin); // w^4
		
		// ========================== Pion =======================
		// both charge neg eta region
		uPNReTPCPionNegEta += fCMEQRePionBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPionNegEta += fCMEQImPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPionNegEta += fCMEQRePionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPionNegEta += fCMEQImPionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPionNegEta += fCMEQRePionBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPionNegEta += fCMEQImPionBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPionNegEta += fCMEMultPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPionNegEta += fCMEMultPionBothCharge[1][h]->GetBinContent(EBin);
		
		// pegitive charge negitive eta region
		uPReTPCPionNegEta += fCMEQRePion[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPionNegEta += fCMEQImPion[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPionNegEta += fCMEQRePion[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPionNegEta += fCMEQImPion[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPionNegEta += fCMEQRePion[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPionNegEta += fCMEQImPion[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPionNegEta += fCMEMultPion[0][h]->GetBinContent(EBin); // w
		uP2MTPCPionNegEta += fCMEMultPion[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCPionNegEta += fCMEQ2Re4Pion[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCPionNegEta += fCMEQ2Im4Pion[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCPionNegEta += fCMEQ3Re2Pion[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCPionNegEta += fCMEQ3Im2Pion[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCPionNegEta += fCMEw0Pion[0]->GetBinContent(EBin); // w^0
		uP3MTPCPionNegEta += fCMEw3Pion[0]->GetBinContent(EBin); // w^3
		uP4MTPCPionNegEta += fCMEw4Pion[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCPionNegEta += fCMEQRePion[1][h]->GetBinContent(EBin);
		uNImTPCPionNegEta += fCMEQImPion[1][h]->GetBinContent(EBin);
		uN2ReTPCPionNegEta += fCMEQRePion[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPionNegEta += fCMEQImPion[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPionNegEta += fCMEQRePion[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPionNegEta += fCMEQImPion[3][h+1]->GetBinContent(EBin);
		uNMTPCPionNegEta += fCMEMultPion[1][h]->GetBinContent(EBin);
		uN2MTPCPionNegEta += fCMEMultPion[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCPionNegEta += fCMEQ2Re4Pion[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCPionNegEta += fCMEQ2Im4Pion[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCPionNegEta += fCMEQ3Re2Pion[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCPionNegEta += fCMEQ3Im2Pion[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCPionNegEta += fCMEw0Pion[1]->GetBinContent(EBin); // w^0
		uN3MTPCPionNegEta += fCMEw3Pion[1]->GetBinContent(EBin); // w^3
		uN4MTPCPionNegEta += fCMEw4Pion[1]->GetBinContent(EBin); // w^4
		
		// ========================== Kaon =======================
		// both charge Neg eta region
		uPNReTPCKaonNegEta += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCKaonNegEta += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCKaonNegEta += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCKaonNegEta += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCKaonNegEta += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCKaonNegEta += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCKaonNegEta += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCKaonNegEta += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge negitive eta region
		uPReTPCKaonNegEta += fCMEQReKaon[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCKaonNegEta += fCMEQImKaon[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCKaonNegEta += fCMEQReKaon[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCKaonNegEta += fCMEQImKaon[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCKaonNegEta += fCMEQReKaon[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCKaonNegEta += fCMEQImKaon[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCKaonNegEta += fCMEMultKaon[0][h]->GetBinContent(EBin); // w
		uP2MTPCKaonNegEta += fCMEMultKaon[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCKaonNegEta += fCMEQ2Re4Kaon[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCKaonNegEta += fCMEQ2Im4Kaon[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCKaonNegEta += fCMEQ3Re2Kaon[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCKaonNegEta += fCMEQ3Im2Kaon[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCKaonNegEta += fCMEw0Kaon[0]->GetBinContent(EBin); // w^0
		uP3MTPCKaonNegEta += fCMEw3Kaon[0]->GetBinContent(EBin); // w^3
		uP4MTPCKaonNegEta += fCMEw4Kaon[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCKaonNegEta += fCMEQReKaon[1][h]->GetBinContent(EBin);
		uNImTPCKaonNegEta += fCMEQImKaon[1][h]->GetBinContent(EBin);
		uN2ReTPCKaonNegEta += fCMEQReKaon[1][h+1]->GetBinContent(EBin);
		uN2ImTPCKaonNegEta += fCMEQImKaon[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCKaonNegEta += fCMEQReKaon[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCKaonNegEta += fCMEQImKaon[3][h+1]->GetBinContent(EBin);
		uNMTPCKaonNegEta += fCMEMultKaon[1][h]->GetBinContent(EBin);
		uN2MTPCKaonNegEta += fCMEMultKaon[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCKaonNegEta += fCMEQ2Re4Kaon[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCKaonNegEta += fCMEQ2Im4Kaon[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCKaonNegEta += fCMEQ3Re2Kaon[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCKaonNegEta += fCMEQ3Im2Kaon[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCKaonNegEta += fCMEw0Kaon[1]->GetBinContent(EBin); // w^0
		uN3MTPCKaonNegEta += fCMEw3Kaon[1]->GetBinContent(EBin); // w^3
		uN4MTPCKaonNegEta += fCMEw4Kaon[1]->GetBinContent(EBin); // w^4
		
		// ========================== Proton =======================
		// both charge neg eta region
		uPNReTPCProtonNegEta += fCMEQReProtonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCProtonNegEta += fCMEQImProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCProtonNegEta += fCMEQReProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCProtonNegEta += fCMEQImProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCProtonNegEta += fCMEQReProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCProtonNegEta += fCMEQImProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCProtonNegEta += fCMEMultProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCProtonNegEta += fCMEMultProtonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge negitive eta region
		uPReTPCProtonNegEta += fCMEQReProton[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCProtonNegEta += fCMEQImProton[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCProtonNegEta += fCMEQReProton[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCProtonNegEta += fCMEQImProton[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCProtonNegEta += fCMEQReProton[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCProtonNegEta += fCMEQImProton[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCProtonNegEta += fCMEMultProton[0][h]->GetBinContent(EBin); // w
		uP2MTPCProtonNegEta += fCMEMultProton[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCProtonNegEta += fCMEQ2Re4Proton[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCProtonNegEta += fCMEQ2Im4Proton[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCProtonNegEta += fCMEQ3Re2Proton[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCProtonNegEta += fCMEQ3Im2Proton[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCProtonNegEta += fCMEw0Proton[0]->GetBinContent(EBin); // w^0
		uP3MTPCProtonNegEta += fCMEw3Proton[0]->GetBinContent(EBin); // w^3
		uP4MTPCProtonNegEta += fCMEw4Proton[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCProtonNegEta += fCMEQReProton[1][h]->GetBinContent(EBin);
		uNImTPCProtonNegEta += fCMEQImProton[1][h]->GetBinContent(EBin);
		uN2ReTPCProtonNegEta += fCMEQReProton[1][h+1]->GetBinContent(EBin);
		uN2ImTPCProtonNegEta += fCMEQImProton[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCProtonNegEta += fCMEQReProton[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCProtonNegEta += fCMEQImProton[3][h+1]->GetBinContent(EBin);
		uNMTPCProtonNegEta += fCMEMultProton[1][h]->GetBinContent(EBin);
		uN2MTPCProtonNegEta += fCMEMultProton[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCProtonNegEta += fCMEQ2Re4Proton[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCProtonNegEta += fCMEQ2Im4Proton[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCProtonNegEta += fCMEQ3Re2Proton[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCProtonNegEta += fCMEQ3Im2Proton[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCProtonNegEta += fCMEw0Proton[1]->GetBinContent(EBin); // w^0
		uN3MTPCProtonNegEta += fCMEw3Proton[1]->GetBinContent(EBin); // w^3
		uN4MTPCProtonNegEta += fCMEw4Proton[1]->GetBinContent(EBin); // w^4
	} else if (EBin == fCMEQRe[0][h]->FindBin(0.05)) { // sub pos eta region
		// both charge pos eta region
		uPNReTPCSubPosEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCSubPosEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCSubPosEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCSubPosEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCSubPosEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCSubPosEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCSubPosEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCSubPosEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCSubPosEta += fCMEQRe[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCSubPosEta += fCMEQIm[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCSubPosEta += fCMEQRe[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCSubPosEta += fCMEQIm[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCSubPosEta += fCMEQRe[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCSubPosEta += fCMEQIm[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCSubPosEta += fCMEMult[0][h]->GetBinContent(EBin); // w
		uP2MTPCSubPosEta += fCMEMult[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCSubPosEta += fCMEQ2Re4[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCSubPosEta += fCMEQ2Im4[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCSubPosEta += fCMEQ3Re2[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCSubPosEta += fCMEQ3Im2[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCSubPosEta += fCMEw0[0]->GetBinContent(EBin); // w^0
		uP3MTPCSubPosEta += fCMEw3[0]->GetBinContent(EBin); // w^3
		uP4MTPCSubPosEta += fCMEw4[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCSubPosEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCSubPosEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCSubPosEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCSubPosEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCSubPosEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCSubPosEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCSubPosEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCSubPosEta += fCMEMult[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCSubPosEta += fCMEQ2Re4[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCSubPosEta += fCMEQ2Im4[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCSubPosEta += fCMEQ3Re2[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCSubPosEta += fCMEQ3Im2[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCSubPosEta += fCMEw0[1]->GetBinContent(EBin); // w^0
		uN3MTPCSubPosEta += fCMEw3[1]->GetBinContent(EBin); // w^3
		uN4MTPCSubPosEta += fCMEw4[1]->GetBinContent(EBin); // w^4
		
		// ========================== Pion =======================
		// both charge pos eta region
		uPNReTPCPionSubPosEta += fCMEQRePionBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPionSubPosEta += fCMEQImPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPionSubPosEta += fCMEQRePionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPionSubPosEta += fCMEQImPionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPionSubPosEta += fCMEQRePionBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPionSubPosEta += fCMEQImPionBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPionSubPosEta += fCMEMultPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPionSubPosEta += fCMEMultPionBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCPionSubPosEta += fCMEQRePion[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPionSubPosEta += fCMEQImPion[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPionSubPosEta += fCMEQRePion[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPionSubPosEta += fCMEQImPion[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPionSubPosEta += fCMEQRePion[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPionSubPosEta += fCMEQImPion[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPionSubPosEta += fCMEMultPion[0][h]->GetBinContent(EBin); // w
		uP2MTPCPionSubPosEta += fCMEMultPion[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCPionSubPosEta += fCMEQ2Re4Pion[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCPionSubPosEta += fCMEQ2Im4Pion[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCPionSubPosEta += fCMEQ3Re2Pion[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCPionSubPosEta += fCMEQ3Im2Pion[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCPionSubPosEta += fCMEw0Pion[0]->GetBinContent(EBin); // w^0
		uP3MTPCPionSubPosEta += fCMEw3Pion[0]->GetBinContent(EBin); // w^3
		uP4MTPCPionSubPosEta += fCMEw4Pion[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCPionSubPosEta += fCMEQRePion[1][h]->GetBinContent(EBin);
		uNImTPCPionSubPosEta += fCMEQImPion[1][h]->GetBinContent(EBin);
		uN2ReTPCPionSubPosEta += fCMEQRePion[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPionSubPosEta += fCMEQImPion[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPionSubPosEta += fCMEQRePion[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPionSubPosEta += fCMEQImPion[3][h+1]->GetBinContent(EBin);
		uNMTPCPionSubPosEta += fCMEMultPion[1][h]->GetBinContent(EBin);
		uN2MTPCPionSubPosEta += fCMEMultPion[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCPionSubPosEta += fCMEQ2Re4Pion[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCPionSubPosEta += fCMEQ2Im4Pion[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCPionSubPosEta += fCMEQ3Re2Pion[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCPionSubPosEta += fCMEQ3Im2Pion[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCPionSubPosEta += fCMEw0Pion[1]->GetBinContent(EBin); // w^0
		uN3MTPCPionSubPosEta += fCMEw3Pion[1]->GetBinContent(EBin); // w^3
		uN4MTPCPionSubPosEta += fCMEw4Pion[1]->GetBinContent(EBin); // w^4
		
		// ========================== Kaon =======================
		// both charge pos eta region
		uPNReTPCKaonSubPosEta += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCKaonSubPosEta += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCKaonSubPosEta += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCKaonSubPosEta += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCKaonSubPosEta += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCKaonSubPosEta += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCKaonSubPosEta += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCKaonSubPosEta += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCKaonSubPosEta += fCMEQReKaon[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCKaonSubPosEta += fCMEQImKaon[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCKaonSubPosEta += fCMEQReKaon[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCKaonSubPosEta += fCMEQImKaon[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCKaonSubPosEta += fCMEQReKaon[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCKaonSubPosEta += fCMEQImKaon[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCKaonSubPosEta += fCMEMultKaon[0][h]->GetBinContent(EBin); // w
		uP2MTPCKaonSubPosEta += fCMEMultKaon[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCKaonSubPosEta += fCMEQ2Re4Kaon[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCKaonSubPosEta += fCMEQ2Im4Kaon[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCKaonSubPosEta += fCMEQ3Re2Kaon[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCKaonSubPosEta += fCMEQ3Im2Kaon[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCKaonSubPosEta += fCMEw0Kaon[0]->GetBinContent(EBin); // w^0
		uP3MTPCKaonSubPosEta += fCMEw3Kaon[0]->GetBinContent(EBin); // w^3
		uP4MTPCKaonSubPosEta += fCMEw4Kaon[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCKaonSubPosEta += fCMEQReKaon[1][h]->GetBinContent(EBin);
		uNImTPCKaonSubPosEta += fCMEQImKaon[1][h]->GetBinContent(EBin);
		uN2ReTPCKaonSubPosEta += fCMEQReKaon[1][h+1]->GetBinContent(EBin);
		uN2ImTPCKaonSubPosEta += fCMEQImKaon[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCKaonSubPosEta += fCMEQReKaon[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCKaonSubPosEta += fCMEQImKaon[3][h+1]->GetBinContent(EBin);
		uNMTPCKaonSubPosEta += fCMEMultKaon[1][h]->GetBinContent(EBin);
		uN2MTPCKaonSubPosEta += fCMEMultKaon[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCKaonSubPosEta += fCMEQ2Re4Kaon[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCKaonSubPosEta += fCMEQ2Im4Kaon[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCKaonSubPosEta += fCMEQ3Re2Kaon[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCKaonSubPosEta += fCMEQ3Im2Kaon[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCKaonSubPosEta += fCMEw0Kaon[1]->GetBinContent(EBin); // w^0
		uN3MTPCKaonSubPosEta += fCMEw3Kaon[1]->GetBinContent(EBin); // w^3
		uN4MTPCKaonSubPosEta += fCMEw4Kaon[1]->GetBinContent(EBin); // w^4
		
		// ========================== Proton =======================
		// both charge pos eta region
		uPNReTPCProtonSubPosEta += fCMEQReProtonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCProtonSubPosEta += fCMEQImProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCProtonSubPosEta += fCMEQReProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCProtonSubPosEta += fCMEQImProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCProtonSubPosEta += fCMEQReProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCProtonSubPosEta += fCMEQImProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCProtonSubPosEta += fCMEMultProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCProtonSubPosEta += fCMEMultProtonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge positive eta region
		uPReTPCProtonSubPosEta += fCMEQReProton[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCProtonSubPosEta += fCMEQImProton[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCProtonSubPosEta += fCMEQReProton[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCProtonSubPosEta += fCMEQImProton[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCProtonSubPosEta += fCMEQReProton[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCProtonSubPosEta += fCMEQImProton[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCProtonSubPosEta += fCMEMultProton[0][h]->GetBinContent(EBin); // w
		uP2MTPCProtonSubPosEta += fCMEMultProton[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCProtonSubPosEta += fCMEQ2Re4Proton[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCProtonSubPosEta += fCMEQ2Im4Proton[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCProtonSubPosEta += fCMEQ3Re2Proton[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCProtonSubPosEta += fCMEQ3Im2Proton[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCProtonSubPosEta += fCMEw0Proton[0]->GetBinContent(EBin); // w^0
		uP3MTPCProtonSubPosEta += fCMEw3Proton[0]->GetBinContent(EBin); // w^3
		uP4MTPCProtonSubPosEta += fCMEw4Proton[0]->GetBinContent(EBin); // w^4
		
		// negative charge positive eta region
		uNReTPCProtonSubPosEta += fCMEQReProton[1][h]->GetBinContent(EBin);
		uNImTPCProtonSubPosEta += fCMEQImProton[1][h]->GetBinContent(EBin);
		uN2ReTPCProtonSubPosEta += fCMEQReProton[1][h+1]->GetBinContent(EBin);
		uN2ImTPCProtonSubPosEta += fCMEQImProton[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCProtonSubPosEta += fCMEQReProton[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCProtonSubPosEta += fCMEQImProton[3][h+1]->GetBinContent(EBin);
		uNMTPCProtonSubPosEta += fCMEMultProton[1][h]->GetBinContent(EBin);
		uN2MTPCProtonSubPosEta += fCMEMultProton[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCProtonSubPosEta += fCMEQ2Re4Proton[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCProtonSubPosEta += fCMEQ2Im4Proton[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCProtonSubPosEta += fCMEQ3Re2Proton[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCProtonSubPosEta += fCMEQ3Im2Proton[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCProtonSubPosEta += fCMEw0Proton[1]->GetBinContent(EBin); // w^0
		uN3MTPCProtonSubPosEta += fCMEw3Proton[1]->GetBinContent(EBin); // w^3
		uN4MTPCProtonSubPosEta += fCMEw4Proton[1]->GetBinContent(EBin); // w^4
	} else if (EBin == fCMEQRe[0][h]->FindBin(-0.05)) { // negative eta region
		// both charge neg eta region
		uPNReTPCSubNegEta += fCMEQReBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCSubNegEta += fCMEQImBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCSubNegEta += fCMEQReBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCSubNegEta += fCMEQImBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCSubNegEta += fCMEQReBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCSubNegEta += fCMEQImBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCSubNegEta += fCMEMultBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCSubNegEta += fCMEMultBothCharge[1][h]->GetBinContent(EBin);
		
		// pegitive charge negitive eta region
		uPReTPCSubNegEta += fCMEQRe[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCSubNegEta += fCMEQIm[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCSubNegEta += fCMEQRe[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCSubNegEta += fCMEQIm[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCSubNegEta += fCMEQRe[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCSubNegEta += fCMEQIm[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCSubNegEta += fCMEMult[0][h]->GetBinContent(EBin); // w
		uP2MTPCSubNegEta += fCMEMult[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCSubNegEta += fCMEQ2Re4[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCSubNegEta += fCMEQ2Im4[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCSubNegEta += fCMEQ3Re2[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCSubNegEta += fCMEQ3Im2[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCSubNegEta += fCMEw0[0]->GetBinContent(EBin); // w^0
		uP3MTPCSubNegEta += fCMEw3[0]->GetBinContent(EBin); // w^3
		uP4MTPCSubNegEta += fCMEw4[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCSubNegEta += fCMEQRe[1][h]->GetBinContent(EBin);
		uNImTPCSubNegEta += fCMEQIm[1][h]->GetBinContent(EBin);
		uN2ReTPCSubNegEta += fCMEQRe[1][h+1]->GetBinContent(EBin);
		uN2ImTPCSubNegEta += fCMEQIm[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCSubNegEta += fCMEQRe[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCSubNegEta += fCMEQIm[3][h+1]->GetBinContent(EBin);
		uNMTPCSubNegEta += fCMEMult[1][h]->GetBinContent(EBin);
		uN2MTPCSubNegEta += fCMEMult[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCSubNegEta += fCMEQ2Re4[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCSubNegEta += fCMEQ2Im4[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCSubNegEta += fCMEQ3Re2[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCSubNegEta += fCMEQ3Im2[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCSubNegEta += fCMEw0[1]->GetBinContent(EBin); // w^0
		uN3MTPCSubNegEta += fCMEw3[1]->GetBinContent(EBin); // w^3
		uN4MTPCSubNegEta += fCMEw4[1]->GetBinContent(EBin); // w^4
		
		// ========================== Pion =======================
		// both charge neg eta region
		uPNReTPCPionSubNegEta += fCMEQRePionBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCPionSubNegEta += fCMEQImPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCPionSubNegEta += fCMEQRePionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCPionSubNegEta += fCMEQImPionBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCPionSubNegEta += fCMEQRePionBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCPionSubNegEta += fCMEQImPionBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCPionSubNegEta += fCMEMultPionBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCPionSubNegEta += fCMEMultPionBothCharge[1][h]->GetBinContent(EBin);
		
		// pegitive charge negitive eta region
		uPReTPCPionSubNegEta += fCMEQRePion[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCPionSubNegEta += fCMEQImPion[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCPionSubNegEta += fCMEQRePion[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCPionSubNegEta += fCMEQImPion[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCPionSubNegEta += fCMEQRePion[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCPionSubNegEta += fCMEQImPion[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCPionSubNegEta += fCMEMultPion[0][h]->GetBinContent(EBin); // w
		uP2MTPCPionSubNegEta += fCMEMultPion[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCPionSubNegEta += fCMEQ2Re4Pion[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCPionSubNegEta += fCMEQ2Im4Pion[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCPionSubNegEta += fCMEQ3Re2Pion[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCPionSubNegEta += fCMEQ3Im2Pion[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCPionSubNegEta += fCMEw0Pion[0]->GetBinContent(EBin); // w^0
		uP3MTPCPionSubNegEta += fCMEw3Pion[0]->GetBinContent(EBin); // w^3
		uP4MTPCPionSubNegEta += fCMEw4Pion[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCPionSubNegEta += fCMEQRePion[1][h]->GetBinContent(EBin);
		uNImTPCPionSubNegEta += fCMEQImPion[1][h]->GetBinContent(EBin);
		uN2ReTPCPionSubNegEta += fCMEQRePion[1][h+1]->GetBinContent(EBin);
		uN2ImTPCPionSubNegEta += fCMEQImPion[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCPionSubNegEta += fCMEQRePion[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCPionSubNegEta += fCMEQImPion[3][h+1]->GetBinContent(EBin);
		uNMTPCPionSubNegEta += fCMEMultPion[1][h]->GetBinContent(EBin);
		uN2MTPCPionSubNegEta += fCMEMultPion[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCPionSubNegEta += fCMEQ2Re4Pion[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCPionSubNegEta += fCMEQ2Im4Pion[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCPionSubNegEta += fCMEQ3Re2Pion[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCPionSubNegEta += fCMEQ3Im2Pion[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCPionSubNegEta += fCMEw0Pion[1]->GetBinContent(EBin); // w^0
		uN3MTPCPionSubNegEta += fCMEw3Pion[1]->GetBinContent(EBin); // w^3
		uN4MTPCPionSubNegEta += fCMEw4Pion[1]->GetBinContent(EBin); // w^4
		
		// ========================== Kaon =======================
		// both charge Neg eta region
		uPNReTPCKaonSubNegEta += fCMEQReKaonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCKaonSubNegEta += fCMEQImKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCKaonSubNegEta += fCMEQReKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCKaonSubNegEta += fCMEQImKaonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCKaonSubNegEta += fCMEQReKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCKaonSubNegEta += fCMEQImKaonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCKaonSubNegEta += fCMEMultKaonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCKaonSubNegEta += fCMEMultKaonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge negitive eta region
		uPReTPCKaonSubNegEta += fCMEQReKaon[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCKaonSubNegEta += fCMEQImKaon[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCKaonSubNegEta += fCMEQReKaon[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCKaonSubNegEta += fCMEQImKaon[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCKaonSubNegEta += fCMEQReKaon[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCKaonSubNegEta += fCMEQImKaon[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCKaonSubNegEta += fCMEMultKaon[0][h]->GetBinContent(EBin); // w
		uP2MTPCKaonSubNegEta += fCMEMultKaon[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCKaonSubNegEta += fCMEQ2Re4Kaon[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCKaonSubNegEta += fCMEQ2Im4Kaon[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCKaonSubNegEta += fCMEQ3Re2Kaon[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCKaonSubNegEta += fCMEQ3Im2Kaon[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCKaonSubNegEta += fCMEw0Kaon[0]->GetBinContent(EBin); // w^0
		uP3MTPCKaonSubNegEta += fCMEw3Kaon[0]->GetBinContent(EBin); // w^3
		uP4MTPCKaonSubNegEta += fCMEw4Kaon[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCKaonSubNegEta += fCMEQReKaon[1][h]->GetBinContent(EBin);
		uNImTPCKaonSubNegEta += fCMEQImKaon[1][h]->GetBinContent(EBin);
		uN2ReTPCKaonSubNegEta += fCMEQReKaon[1][h+1]->GetBinContent(EBin);
		uN2ImTPCKaonSubNegEta += fCMEQImKaon[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCKaonSubNegEta += fCMEQReKaon[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCKaonSubNegEta += fCMEQImKaon[3][h+1]->GetBinContent(EBin);
		uNMTPCKaonSubNegEta += fCMEMultKaon[1][h]->GetBinContent(EBin);
		uN2MTPCKaonSubNegEta += fCMEMultKaon[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCKaonSubNegEta += fCMEQ2Re4Kaon[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCKaonSubNegEta += fCMEQ2Im4Kaon[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCKaonSubNegEta += fCMEQ3Re2Kaon[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCKaonSubNegEta += fCMEQ3Im2Kaon[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCKaonSubNegEta += fCMEw0Kaon[1]->GetBinContent(EBin); // w^0
		uN3MTPCKaonSubNegEta += fCMEw3Kaon[1]->GetBinContent(EBin); // w^3
		uN4MTPCKaonSubNegEta += fCMEw4Kaon[1]->GetBinContent(EBin); // w^4
		
		// ========================== Proton =======================
		// both charge neg eta region
		uPNReTPCProtonSubNegEta += fCMEQReProtonBothCharge[0][h]->GetBinContent(EBin);
		uPNImTPCProtonSubNegEta += fCMEQImProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2ReTPCProtonSubNegEta += fCMEQReProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2ImTPCProtonSubNegEta += fCMEQImProtonBothCharge[0][h+1]->GetBinContent(EBin);
		uPN2Re2TPCProtonSubNegEta += fCMEQReProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPN2Im2TPCProtonSubNegEta += fCMEQImProtonBothCharge[1][h+1]->GetBinContent(EBin);
		uPNMTPCProtonSubNegEta += fCMEMultProtonBothCharge[0][h]->GetBinContent(EBin);
		uPN2MTPCProtonSubNegEta += fCMEMultProtonBothCharge[1][h]->GetBinContent(EBin);
		
		// positive charge negitive eta region
		uPReTPCProtonSubNegEta += fCMEQReProton[0][h]->GetBinContent(EBin); // w*cos(phi)
		uPImTPCProtonSubNegEta += fCMEQImProton[0][h]->GetBinContent(EBin); // w*sin(phi)
		uP2ReTPCProtonSubNegEta += fCMEQReProton[0][h+1]->GetBinContent(EBin); // w*cos(2phi)
		uP2ImTPCProtonSubNegEta += fCMEQImProton[0][h+1]->GetBinContent(EBin); // w*sin(2phi)
		uP2Re2TPCProtonSubNegEta += fCMEQReProton[2][h+1]->GetBinContent(EBin); // w^2*cos(2phi)
		uP2Im2TPCProtonSubNegEta += fCMEQImProton[2][h+1]->GetBinContent(EBin); // w^2*sin(2phi)
		uPMTPCProtonSubNegEta += fCMEMultProton[0][h]->GetBinContent(EBin); // w
		uP2MTPCProtonSubNegEta += fCMEMultProton[2][h]->GetBinContent(EBin); // w^2
		
		uP4Re2TPCProtonSubNegEta += fCMEQ2Re4Proton[0]->GetBinContent(EBin); // w^2*cos(4phi)
		uP4Im2TPCProtonSubNegEta += fCMEQ2Im4Proton[0]->GetBinContent(EBin); // w^2*sin(4phi)
		uP2Re3TPCProtonSubNegEta += fCMEQ3Re2Proton[0]->GetBinContent(EBin); // w^3*cos(2phi)
		uP2Im3TPCProtonSubNegEta += fCMEQ3Im2Proton[0]->GetBinContent(EBin); // w^3*sin(2phi)
		uP0MTPCProtonSubNegEta += fCMEw0Proton[0]->GetBinContent(EBin); // w^0
		uP3MTPCProtonSubNegEta += fCMEw3Proton[0]->GetBinContent(EBin); // w^3
		uP4MTPCProtonSubNegEta += fCMEw4Proton[0]->GetBinContent(EBin); // w^4
		
		// negative charge negitive eta region
		uNReTPCProtonSubNegEta += fCMEQReProton[1][h]->GetBinContent(EBin);
		uNImTPCProtonSubNegEta += fCMEQImProton[1][h]->GetBinContent(EBin);
		uN2ReTPCProtonSubNegEta += fCMEQReProton[1][h+1]->GetBinContent(EBin);
		uN2ImTPCProtonSubNegEta += fCMEQImProton[1][h+1]->GetBinContent(EBin);
		uN2Re2TPCProtonSubNegEta += fCMEQReProton[3][h+1]->GetBinContent(EBin);
		uN2Im2TPCProtonSubNegEta += fCMEQImProton[3][h+1]->GetBinContent(EBin);
		uNMTPCProtonSubNegEta += fCMEMultProton[1][h]->GetBinContent(EBin);
		uN2MTPCProtonSubNegEta += fCMEMultProton[3][h]->GetBinContent(EBin);
		
		uN4Re2TPCProtonSubNegEta += fCMEQ2Re4Proton[1]->GetBinContent(EBin); // w^2*cos(4phi)
		uN4Im2TPCProtonSubNegEta += fCMEQ2Im4Proton[1]->GetBinContent(EBin); // w^2*sin(4phi)
		uN2Re3TPCProtonSubNegEta += fCMEQ3Re2Proton[1]->GetBinContent(EBin); // w^3*cos(2phi)
		uN2Im3TPCProtonSubNegEta += fCMEQ3Im2Proton[1]->GetBinContent(EBin); // w^3*sin(2phi)
		uN0MTPCProtonSubNegEta += fCMEw0Proton[1]->GetBinContent(EBin); // w^0
		uN3MTPCProtonSubNegEta += fCMEw3Proton[1]->GetBinContent(EBin); // w^3
		uN4MTPCProtonSubNegEta += fCMEw4Proton[1]->GetBinContent(EBin); // w^4
		
	}
  }
  
  // set fpQvecEvent
  // Pos Eta 0.1<|eta|<0.8, Neg Eta -0.8<|eta|<-0.1
  //cout<<"==> uPReTPCPosEta + uNReTPCPosEta = "<<uPReTPCPosEta + uNReTPCPosEta<<", uPNReTPCPosEta = "<<uPNReTPCPosEta<<endl;
  //cout<<"==> uPNReTPC = "<<uPNReTPC<<", uPReTPCSubNegEta + uPReTPCNegEta + uPReTPCPosEta + uPReTPCSubPosEta + uNReTPCSubNegEta + uNReTPCNegEta + uNReTPCPosEta + uNReTPCSubPosEta = "<<uPReTPCSubNegEta + uPReTPCNegEta + uPReTPCPosEta + uPReTPCSubPosEta + uNReTPCSubNegEta + uNReTPCNegEta + uNReTPCPosEta + uNReTPCSubPosEta<<endl;
  fpQvecEvent->setTPCRePosChPosEta( uPReTPCPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCImPosChPosEta( uPImTPCPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPC2RePosChPosEta( uP2ReTPCPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPC2ImPosChPosEta( uP2ImTPCPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPC2Re2PosChPosEta( uP2Re2TPCPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPC2Im2PosChPosEta( uP2Im2TPCPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCMPosChPosEta( uPMTPCPosEta );   // w ch+ eta+
  fpQvecEvent->setTPC2MPosChPosEta( uP2MTPCPosEta );   // w^2 ch+ eta+
  fpQvecEvent->setTPC4Re2PosChPosEta( uP4Re2TPCPosEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2PosChPosEta( uP4Im2TPCPosEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3PosChPosEta( uP2Re3TPCPosEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3PosChPosEta( uP2Im3TPCPosEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MPosChPosEta( uP0MTPCPosEta );   // w^0 ch+ eta+
  fpQvecEvent->setTPC3MPosChPosEta( uP3MTPCPosEta );   // w^3 ch+ eta+
  fpQvecEvent->setTPC4MPosChPosEta( uP4MTPCPosEta );   // w^4 ch+ eta+
  
  fpQvecEvent->setTPCRePosChNegEta( uPReTPCNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCImPosChNegEta( uPImTPCNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPC2RePosChNegEta( uP2ReTPCNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPC2ImPosChNegEta( uP2ImTPCNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPC2Re2PosChNegEta( uP2Re2TPCNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPC2Im2PosChNegEta( uP2Im2TPCNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCMPosChNegEta( uPMTPCNegEta );   // w ch+ eta-
  fpQvecEvent->setTPC2MPosChNegEta( uP2MTPCNegEta );   // w^2 ch+ eta-
  fpQvecEvent->setTPC4Re2PosChNegEta( uP4Re2TPCNegEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2PosChNegEta( uP4Im2TPCNegEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3PosChNegEta( uP2Re3TPCNegEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3PosChNegEta( uP2Im3TPCNegEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MPosChNegEta( uP0MTPCNegEta );   // w^0 ch+ eta-
  fpQvecEvent->setTPC3MPosChNegEta( uP3MTPCNegEta );   // w^3 ch+ eta-
  fpQvecEvent->setTPC4MPosChNegEta( uP4MTPCNegEta );   // w^4 ch+ eta-
  
  fpQvecEvent->setTPCReNegChPosEta( uNReTPCPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCImNegChPosEta( uNImTPCPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPC2ReNegChPosEta( uN2ReTPCPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPC2ImNegChPosEta( uN2ImTPCPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPC2Re2NegChPosEta( uN2Re2TPCPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPC2Im2NegChPosEta( uN2Im2TPCPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCMNegChPosEta( uNMTPCPosEta );   // w ch- eta+
  fpQvecEvent->setTPC2MNegChPosEta( uN2MTPCPosEta );   // w^2  h- eta+
  fpQvecEvent->setTPC4Re2NegChPosEta( uN4Re2TPCPosEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2NegChPosEta( uN4Im2TPCPosEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3NegChPosEta( uN2Re3TPCPosEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3NegChPosEta( uN2Im3TPCPosEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MNegChPosEta( uN0MTPCPosEta );   // w^0 ch- eta+
  fpQvecEvent->setTPC3MNegChPosEta( uN3MTPCPosEta );   // w^3 ch- eta+
  fpQvecEvent->setTPC4MNegChPosEta( uN4MTPCPosEta );   // w^4 ch- eta+
  
  fpQvecEvent->setTPCReNegChNegEta( uNReTPCNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCImNegChNegEta( uNImTPCNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPC2ReNegChNegEta( uN2ReTPCNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPC2ImNegChNegEta( uN2ImTPCNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPC2Re2NegChNegEta( uN2Re2TPCNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPC2Im2NegChNegEta( uN2Im2TPCNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCMNegChNegEta( uNMTPCNegEta );   // w ch- eta-
  fpQvecEvent->setTPC2MNegChNegEta( uN2MTPCNegEta );   // w^2 ch- eta-
  fpQvecEvent->setTPC4Re2NegChNegEta( uN4Re2TPCNegEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2NegChNegEta( uN4Im2TPCNegEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3NegChNegEta( uN2Re3TPCNegEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3NegChNegEta( uN2Im3TPCNegEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MNegChNegEta( uN0MTPCNegEta );   // w^0 ch- eta-
  fpQvecEvent->setTPC3MNegChNegEta( uN3MTPCNegEta );   // w^3 ch- eta-
  fpQvecEvent->setTPC4MNegChNegEta( uN4MTPCNegEta );   // w^4 ch- eta-
  
  fpQvecEvent->setTPCPionRePosChPosEta( uPReTPCPionPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCPionImPosChPosEta( uPImTPCPionPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCPion2RePosChPosEta( uP2ReTPCPionPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCPion2ImPosChPosEta( uP2ImTPCPionPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCPion2Re2PosChPosEta( uP2Re2TPCPionPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCPion2Im2PosChPosEta( uP2Im2TPCPionPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCPionMPosChPosEta( uPMTPCPionPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCPion2MPosChPosEta( uP2MTPCPionPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCPionRePosChNegEta( uPReTPCPionNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCPionImPosChNegEta( uPImTPCPionNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCPion2RePosChNegEta( uP2ReTPCPionNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCPion2ImPosChNegEta( uP2ImTPCPionNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCPion2Re2PosChNegEta( uP2Re2TPCPionNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCPion2Im2PosChNegEta( uP2Im2TPCPionNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCPionMPosChNegEta( uPMTPCPionNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCPion2MPosChNegEta( uP2MTPCPionNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCPionReNegChPosEta( uNReTPCPionPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCPionImNegChPosEta( uNImTPCPionPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCPion2ReNegChPosEta( uN2ReTPCPionPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCPion2ImNegChPosEta( uN2ImTPCPionPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCPion2Re2NegChPosEta( uN2Re2TPCPionPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCPion2Im2NegChPosEta( uN2Im2TPCPionPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCPionMNegChPosEta( uNMTPCPionPosEta );   // w ch- eta+
  fpQvecEvent->setTPCPion2MNegChPosEta( uN2MTPCPionPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCPionReNegChNegEta( uNReTPCPionNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCPionImNegChNegEta( uNImTPCPionNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCPion2ReNegChNegEta( uN2ReTPCPionNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCPion2ImNegChNegEta( uN2ImTPCPionNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCPion2Re2NegChNegEta( uN2Re2TPCPionNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCPion2Im2NegChNegEta( uN2Im2TPCPionNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCPionMNegChNegEta( uNMTPCPionNegEta );   // w ch- eta-
  fpQvecEvent->setTPCPion2MNegChNegEta( uN2MTPCPionNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCPion4Re2PosChPosEta( uP4Re2TPCPionPosEta );
  fpQvecEvent->setTPCPion4Im2PosChPosEta( uP4Im2TPCPionPosEta );
  fpQvecEvent->setTPCPion2Re3PosChPosEta( uP2Re3TPCPionPosEta );
  fpQvecEvent->setTPCPion2Im3PosChPosEta( uP2Im3TPCPionPosEta );
  fpQvecEvent->setTPCPion0MPosChPosEta( uP0MTPCPionPosEta );
  fpQvecEvent->setTPCPion3MPosChPosEta( uP3MTPCPionPosEta );
  fpQvecEvent->setTPCPion4MPosChPosEta( uP4MTPCPionPosEta );
  fpQvecEvent->setTPCPion4Re2PosChNegEta( uP4Re2TPCPionNegEta );
  fpQvecEvent->setTPCPion4Im2PosChNegEta( uP4Im2TPCPionNegEta );
  fpQvecEvent->setTPCPion2Re3PosChNegEta( uP2Re3TPCPionNegEta );
  fpQvecEvent->setTPCPion2Im3PosChNegEta( uP2Im3TPCPionNegEta );
  fpQvecEvent->setTPCPion0MPosChNegEta( uP0MTPCPionNegEta );
  fpQvecEvent->setTPCPion3MPosChNegEta( uP3MTPCPionNegEta );
  fpQvecEvent->setTPCPion4MPosChNegEta( uP4MTPCPionNegEta );
  fpQvecEvent->setTPCPion4Re2NegChPosEta( uN4Re2TPCPionPosEta );
  fpQvecEvent->setTPCPion4Im2NegChPosEta( uN4Im2TPCPionPosEta );
  fpQvecEvent->setTPCPion2Re3NegChPosEta( uN2Re3TPCPionPosEta );
  fpQvecEvent->setTPCPion2Im3NegChPosEta( uN2Im3TPCPionPosEta );
  fpQvecEvent->setTPCPion0MNegChPosEta( uN0MTPCPionPosEta );
  fpQvecEvent->setTPCPion3MNegChPosEta( uN3MTPCPionPosEta );
  fpQvecEvent->setTPCPion4MNegChPosEta( uN4MTPCPionPosEta );
  fpQvecEvent->setTPCPion4Re2NegChNegEta( uN4Re2TPCPionNegEta );
  fpQvecEvent->setTPCPion4Im2NegChNegEta( uN4Im2TPCPionNegEta );
  fpQvecEvent->setTPCPion2Re3NegChNegEta( uN2Re3TPCPionNegEta );
  fpQvecEvent->setTPCPion2Im3NegChNegEta( uN2Im3TPCPionNegEta );
  fpQvecEvent->setTPCPion0MNegChNegEta( uN0MTPCPionNegEta );
  fpQvecEvent->setTPCPion3MNegChNegEta( uN3MTPCPionNegEta );
  fpQvecEvent->setTPCPion4MNegChNegEta( uN4MTPCPionNegEta );
  
  fpQvecEvent->setTPCKaonRePosChPosEta( uPReTPCKaonPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCKaonImPosChPosEta( uPImTPCKaonPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCKaon2RePosChPosEta( uP2ReTPCKaonPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCKaon2ImPosChPosEta( uP2ImTPCKaonPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCKaon2Re2PosChPosEta( uP2Re2TPCKaonPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCKaon2Im2PosChPosEta( uP2Im2TPCKaonPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCKaonMPosChPosEta( uPMTPCKaonPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCKaon2MPosChPosEta( uP2MTPCKaonPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCKaonRePosChNegEta( uPReTPCKaonNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCKaonImPosChNegEta( uPImTPCKaonNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCKaon2RePosChNegEta( uP2ReTPCKaonNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCKaon2ImPosChNegEta( uP2ImTPCKaonNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCKaon2Re2PosChNegEta( uP2Re2TPCKaonNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCKaon2Im2PosChNegEta( uP2Im2TPCKaonNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCKaonMPosChNegEta( uPMTPCKaonNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCKaon2MPosChNegEta( uP2MTPCKaonNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCKaonReNegChPosEta( uNReTPCKaonPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCKaonImNegChPosEta( uNImTPCKaonPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCKaon2ReNegChPosEta( uN2ReTPCKaonPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCKaon2ImNegChPosEta( uN2ImTPCKaonPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCKaon2Re2NegChPosEta( uN2Re2TPCKaonPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCKaon2Im2NegChPosEta( uN2Im2TPCKaonPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCKaonMNegChPosEta( uNMTPCKaonPosEta );   // w ch- eta+
  fpQvecEvent->setTPCKaon2MNegChPosEta( uN2MTPCKaonPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCKaonReNegChNegEta( uNReTPCKaonNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCKaonImNegChNegEta( uNImTPCKaonNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCKaon2ReNegChNegEta( uN2ReTPCKaonNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCKaon2ImNegChNegEta( uN2ImTPCKaonNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCKaon2Re2NegChNegEta( uN2Re2TPCKaonNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCKaon2Im2NegChNegEta( uN2Im2TPCKaonNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCKaonMNegChNegEta( uNMTPCKaonNegEta );   // w ch- eta-
  fpQvecEvent->setTPCKaon2MNegChNegEta( uN2MTPCKaonNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCKaon4Re2PosChPosEta( uP4Re2TPCKaonPosEta );
  fpQvecEvent->setTPCKaon4Im2PosChPosEta( uP4Im2TPCKaonPosEta );
  fpQvecEvent->setTPCKaon2Re3PosChPosEta( uP2Re3TPCKaonPosEta );
  fpQvecEvent->setTPCKaon2Im3PosChPosEta( uP2Im3TPCKaonPosEta );
  fpQvecEvent->setTPCKaon0MPosChPosEta( uP0MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon3MPosChPosEta( uP3MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon4MPosChPosEta( uP4MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon4Re2PosChNegEta( uP4Re2TPCKaonNegEta );
  fpQvecEvent->setTPCKaon4Im2PosChNegEta( uP4Im2TPCKaonNegEta );
  fpQvecEvent->setTPCKaon2Re3PosChNegEta( uP2Re3TPCKaonNegEta );
  fpQvecEvent->setTPCKaon2Im3PosChNegEta( uP2Im3TPCKaonNegEta );
  fpQvecEvent->setTPCKaon0MPosChNegEta( uP0MTPCKaonNegEta );
  fpQvecEvent->setTPCKaon3MPosChNegEta( uP3MTPCKaonNegEta );
  fpQvecEvent->setTPCKaon4MPosChNegEta( uP4MTPCKaonNegEta );
  fpQvecEvent->setTPCKaon4Re2NegChPosEta( uN4Re2TPCKaonPosEta );
  fpQvecEvent->setTPCKaon4Im2NegChPosEta( uN4Im2TPCKaonPosEta );
  fpQvecEvent->setTPCKaon2Re3NegChPosEta( uN2Re3TPCKaonPosEta );
  fpQvecEvent->setTPCKaon2Im3NegChPosEta( uN2Im3TPCKaonPosEta );
  fpQvecEvent->setTPCKaon0MNegChPosEta( uN0MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon3MNegChPosEta( uN3MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon4MNegChPosEta( uN4MTPCKaonPosEta );
  fpQvecEvent->setTPCKaon4Re2NegChNegEta( uN4Re2TPCKaonNegEta );
  fpQvecEvent->setTPCKaon4Im2NegChNegEta( uN4Im2TPCKaonNegEta );
  fpQvecEvent->setTPCKaon2Re3NegChNegEta( uN2Re3TPCKaonNegEta );
  fpQvecEvent->setTPCKaon2Im3NegChNegEta( uN2Im3TPCKaonNegEta );
  fpQvecEvent->setTPCKaon0MNegChNegEta( uN0MTPCKaonNegEta );
  fpQvecEvent->setTPCKaon3MNegChNegEta( uN3MTPCKaonNegEta );
  fpQvecEvent->setTPCKaon4MNegChNegEta( uN4MTPCKaonNegEta );
  
  fpQvecEvent->setTPCProtonRePosChPosEta( uPReTPCProtonPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCProtonImPosChPosEta( uPImTPCProtonPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCProton2RePosChPosEta( uP2ReTPCProtonPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCProton2ImPosChPosEta( uP2ImTPCProtonPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCProton2Re2PosChPosEta( uP2Re2TPCProtonPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCProton2Im2PosChPosEta( uP2Im2TPCProtonPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCProtonMPosChPosEta( uPMTPCProtonPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCProton2MPosChPosEta( uP2MTPCProtonPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCProtonRePosChNegEta( uPReTPCProtonNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCProtonImPosChNegEta( uPImTPCProtonNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCProton2RePosChNegEta( uP2ReTPCProtonNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCProton2ImPosChNegEta( uP2ImTPCProtonNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCProton2Re2PosChNegEta( uP2Re2TPCProtonNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCProton2Im2PosChNegEta( uP2Im2TPCProtonNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCProtonMPosChNegEta( uPMTPCProtonNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCProton2MPosChNegEta( uP2MTPCProtonNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCProtonReNegChPosEta( uNReTPCProtonPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCProtonImNegChPosEta( uNImTPCProtonPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCProton2ReNegChPosEta( uN2ReTPCProtonPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCProton2ImNegChPosEta( uN2ImTPCProtonPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCProton2Re2NegChPosEta( uN2Re2TPCProtonPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCProton2Im2NegChPosEta( uN2Im2TPCProtonPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCProtonMNegChPosEta( uNMTPCProtonPosEta );   // w ch- eta+
  fpQvecEvent->setTPCProton2MNegChPosEta( uN2MTPCProtonPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCProtonReNegChNegEta( uNReTPCProtonNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCProtonImNegChNegEta( uNImTPCProtonNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCProton2ReNegChNegEta( uN2ReTPCProtonNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCProton2ImNegChNegEta( uN2ImTPCProtonNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCProton2Re2NegChNegEta( uN2Re2TPCProtonNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCProton2Im2NegChNegEta( uN2Im2TPCProtonNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCProtonMNegChNegEta( uNMTPCProtonNegEta );   // w ch- eta-
  fpQvecEvent->setTPCProton2MNegChNegEta( uN2MTPCProtonNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCProton4Re2PosChPosEta( uP4Re2TPCProtonPosEta );
  fpQvecEvent->setTPCProton4Im2PosChPosEta( uP4Im2TPCProtonPosEta );
  fpQvecEvent->setTPCProton2Re3PosChPosEta( uP2Re3TPCProtonPosEta );
  fpQvecEvent->setTPCProton2Im3PosChPosEta( uP2Im3TPCProtonPosEta );
  fpQvecEvent->setTPCProton0MPosChPosEta( uP0MTPCProtonPosEta );
  fpQvecEvent->setTPCProton3MPosChPosEta( uP3MTPCProtonPosEta );
  fpQvecEvent->setTPCProton4MPosChPosEta( uP4MTPCProtonPosEta );
  fpQvecEvent->setTPCProton4Re2PosChNegEta( uP4Re2TPCProtonNegEta );
  fpQvecEvent->setTPCProton4Im2PosChNegEta( uP4Im2TPCProtonNegEta );
  fpQvecEvent->setTPCProton2Re3PosChNegEta( uP2Re3TPCProtonNegEta );
  fpQvecEvent->setTPCProton2Im3PosChNegEta( uP2Im3TPCProtonNegEta );
  fpQvecEvent->setTPCProton0MPosChNegEta( uP0MTPCProtonNegEta );
  fpQvecEvent->setTPCProton3MPosChNegEta( uP3MTPCProtonNegEta );
  fpQvecEvent->setTPCProton4MPosChNegEta( uP4MTPCProtonNegEta );
  fpQvecEvent->setTPCProton4Re2NegChPosEta( uN4Re2TPCProtonPosEta );
  fpQvecEvent->setTPCProton4Im2NegChPosEta( uN4Im2TPCProtonPosEta );
  fpQvecEvent->setTPCProton2Re3NegChPosEta( uN2Re3TPCProtonPosEta );
  fpQvecEvent->setTPCProton2Im3NegChPosEta( uN2Im3TPCProtonPosEta );
  fpQvecEvent->setTPCProton0MNegChPosEta( uN0MTPCProtonPosEta );
  fpQvecEvent->setTPCProton3MNegChPosEta( uN3MTPCProtonPosEta );
  fpQvecEvent->setTPCProton4MNegChPosEta( uN4MTPCProtonPosEta );
  fpQvecEvent->setTPCProton4Re2NegChNegEta( uN4Re2TPCProtonNegEta );
  fpQvecEvent->setTPCProton4Im2NegChNegEta( uN4Im2TPCProtonNegEta );
  fpQvecEvent->setTPCProton2Re3NegChNegEta( uN2Re3TPCProtonNegEta );
  fpQvecEvent->setTPCProton2Im3NegChNegEta( uN2Im3TPCProtonNegEta );
  fpQvecEvent->setTPCProton0MNegChNegEta( uN0MTPCProtonNegEta );
  fpQvecEvent->setTPCProton3MNegChNegEta( uN3MTPCProtonNegEta );
  fpQvecEvent->setTPCProton4MNegChNegEta( uN4MTPCProtonNegEta );
  
  // SubPos Eta 0<|eta|<0.1, SubNeg Eta -0.1<|eta|<0
  fpQvecEvent->setTPCRePosChSubPosEta( uPReTPCSubPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCImPosChSubPosEta( uPImTPCSubPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPC2RePosChSubPosEta( uP2ReTPCSubPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPC2ImPosChSubPosEta( uP2ImTPCSubPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPC2Re2PosChSubPosEta( uP2Re2TPCSubPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPC2Im2PosChSubPosEta( uP2Im2TPCSubPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCMPosChSubPosEta( uPMTPCSubPosEta );   // w ch+ eta+
  fpQvecEvent->setTPC2MPosChSubPosEta( uP2MTPCSubPosEta );   // w^2 ch+ eta+
  fpQvecEvent->setTPC4Re2PosChSubPosEta( uP4Re2TPCSubPosEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2PosChSubPosEta( uP4Im2TPCSubPosEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3PosChSubPosEta( uP2Re3TPCSubPosEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3PosChSubPosEta( uP2Im3TPCSubPosEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MPosChSubPosEta( uP0MTPCSubPosEta );   // w^0 ch+ eta+
  fpQvecEvent->setTPC3MPosChSubPosEta( uP3MTPCSubPosEta );   // w^3 ch+ eta+
  fpQvecEvent->setTPC4MPosChSubPosEta( uP4MTPCSubPosEta );   // w^4 ch+ eta+
  
  fpQvecEvent->setTPCRePosChSubNegEta( uPReTPCSubNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCImPosChSubNegEta( uPImTPCSubNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPC2RePosChSubNegEta( uP2ReTPCSubNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPC2ImPosChSubNegEta( uP2ImTPCSubNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPC2Re2PosChSubNegEta( uP2Re2TPCSubNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPC2Im2PosChSubNegEta( uP2Im2TPCSubNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCMPosChSubNegEta( uPMTPCSubNegEta );   // w ch+ eta-
  fpQvecEvent->setTPC2MPosChSubNegEta( uP2MTPCSubNegEta );   // w^2 ch+ eta-
  fpQvecEvent->setTPC4Re2PosChSubNegEta( uP4Re2TPCSubNegEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2PosChSubNegEta( uP4Im2TPCSubNegEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3PosChSubNegEta( uP2Re3TPCSubNegEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3PosChSubNegEta( uP2Im3TPCSubNegEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MPosChSubNegEta( uP0MTPCSubNegEta );   // w^0 ch+ eta-
  fpQvecEvent->setTPC3MPosChSubNegEta( uP3MTPCSubNegEta );   // w^3 ch+ eta-
  fpQvecEvent->setTPC4MPosChSubNegEta( uP4MTPCSubNegEta );   // w^4 ch+ eta-
  
  fpQvecEvent->setTPCReNegChSubPosEta( uNReTPCSubPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCImNegChSubPosEta( uNImTPCSubPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPC2ReNegChSubPosEta( uN2ReTPCSubPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPC2ImNegChSubPosEta( uN2ImTPCSubPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPC2Re2NegChSubPosEta( uN2Re2TPCSubPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPC2Im2NegChSubPosEta( uN2Im2TPCSubPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCMNegChSubPosEta( uNMTPCSubPosEta );   // w ch- eta+
  fpQvecEvent->setTPC2MNegChSubPosEta( uN2MTPCSubPosEta );   // w^2  h- eta+
  fpQvecEvent->setTPC4Re2NegChSubPosEta( uN4Re2TPCSubPosEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2NegChSubPosEta( uN4Im2TPCSubPosEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3NegChSubPosEta( uN2Re3TPCSubPosEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3NegChSubPosEta( uN2Im3TPCSubPosEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MNegChSubPosEta( uN0MTPCSubPosEta );   // w^0 ch- eta+
  fpQvecEvent->setTPC3MNegChSubPosEta( uN3MTPCSubPosEta );   // w^3 ch- eta+
  fpQvecEvent->setTPC4MNegChSubPosEta( uN4MTPCSubPosEta );   // w^4 ch- eta+
  
  fpQvecEvent->setTPCReNegChSubNegEta( uNReTPCSubNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCImNegChSubNegEta( uNImTPCSubNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPC2ReNegChSubNegEta( uN2ReTPCSubNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPC2ImNegChSubNegEta( uN2ImTPCSubNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPC2Re2NegChSubNegEta( uN2Re2TPCSubNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPC2Im2NegChSubNegEta( uN2Im2TPCSubNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCMNegChSubNegEta( uNMTPCSubNegEta );   // w ch- eta-
  fpQvecEvent->setTPC2MNegChSubNegEta( uN2MTPCSubNegEta );   // w^2 ch- eta-
  fpQvecEvent->setTPC4Re2NegChSubNegEta( uN4Re2TPCSubNegEta); // w^2*cos(4phi)
  fpQvecEvent->setTPC4Im2NegChSubNegEta( uN4Im2TPCSubNegEta); // w^2*sin(4phi)
  fpQvecEvent->setTPC2Re3NegChSubNegEta( uN2Re3TPCSubNegEta); // w^3*cos(2phi)
  fpQvecEvent->setTPC2Im3NegChSubNegEta( uN2Im3TPCSubNegEta); // w^3*sin(2phi)
  fpQvecEvent->setTPC0MNegChSubNegEta( uN0MTPCSubNegEta );   // w^0 ch- eta-
  fpQvecEvent->setTPC3MNegChSubNegEta( uN3MTPCSubNegEta );   // w^3 ch- eta-
  fpQvecEvent->setTPC4MNegChSubNegEta( uN4MTPCSubNegEta );   // w^4 ch- eta-
  
  fpQvecEvent->setTPCPionRePosChSubPosEta( uPReTPCPionSubPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCPionImPosChSubPosEta( uPImTPCPionSubPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCPion2RePosChSubPosEta( uP2ReTPCPionSubPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCPion2ImPosChSubPosEta( uP2ImTPCPionSubPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCPion2Re2PosChSubPosEta( uP2Re2TPCPionSubPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCPion2Im2PosChSubPosEta( uP2Im2TPCPionSubPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCPionMPosChSubPosEta( uPMTPCPionSubPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCPion2MPosChSubPosEta( uP2MTPCPionSubPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCPionRePosChSubNegEta( uPReTPCPionSubNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCPionImPosChSubNegEta( uPImTPCPionSubNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCPion2RePosChSubNegEta( uP2ReTPCPionSubNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCPion2ImPosChSubNegEta( uP2ImTPCPionSubNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCPion2Re2PosChSubNegEta( uP2Re2TPCPionSubNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCPion2Im2PosChSubNegEta( uP2Im2TPCPionSubNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCPionMPosChSubNegEta( uPMTPCPionSubNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCPion2MPosChSubNegEta( uP2MTPCPionSubNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCPionReNegChSubPosEta( uNReTPCPionSubPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCPionImNegChSubPosEta( uNImTPCPionSubPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCPion2ReNegChSubPosEta( uN2ReTPCPionSubPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCPion2ImNegChSubPosEta( uN2ImTPCPionSubPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCPion2Re2NegChSubPosEta( uN2Re2TPCPionSubPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCPion2Im2NegChSubPosEta( uN2Im2TPCPionSubPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCPionMNegChSubPosEta( uNMTPCPionSubPosEta );   // w ch- eta+
  fpQvecEvent->setTPCPion2MNegChSubPosEta( uN2MTPCPionSubPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCPionReNegChSubNegEta( uNReTPCPionSubNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCPionImNegChSubNegEta( uNImTPCPionSubNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCPion2ReNegChSubNegEta( uN2ReTPCPionSubNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCPion2ImNegChSubNegEta( uN2ImTPCPionSubNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCPion2Re2NegChSubNegEta( uN2Re2TPCPionSubNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCPion2Im2NegChSubNegEta( uN2Im2TPCPionSubNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCPionMNegChSubNegEta( uNMTPCPionSubNegEta );   // w ch- eta-
  fpQvecEvent->setTPCPion2MNegChSubNegEta( uN2MTPCPionSubNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCPion4Re2PosChSubPosEta( uP4Re2TPCPionSubPosEta );
  fpQvecEvent->setTPCPion4Im2PosChSubPosEta( uP4Im2TPCPionSubPosEta );
  fpQvecEvent->setTPCPion2Re3PosChSubPosEta( uP2Re3TPCPionSubPosEta );
  fpQvecEvent->setTPCPion2Im3PosChSubPosEta( uP2Im3TPCPionSubPosEta );
  fpQvecEvent->setTPCPion0MPosChSubPosEta( uP0MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion3MPosChSubPosEta( uP3MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion4MPosChSubPosEta( uP4MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion4Re2PosChSubNegEta( uP4Re2TPCPionSubNegEta );
  fpQvecEvent->setTPCPion4Im2PosChSubNegEta( uP4Im2TPCPionSubNegEta );
  fpQvecEvent->setTPCPion2Re3PosChSubNegEta( uP2Re3TPCPionSubNegEta );
  fpQvecEvent->setTPCPion2Im3PosChSubNegEta( uP2Im3TPCPionSubNegEta );
  fpQvecEvent->setTPCPion0MPosChSubNegEta( uP0MTPCPionSubNegEta );
  fpQvecEvent->setTPCPion3MPosChSubNegEta( uP3MTPCPionSubNegEta );
  fpQvecEvent->setTPCPion4MPosChSubNegEta( uP4MTPCPionSubNegEta );
  fpQvecEvent->setTPCPion4Re2NegChSubPosEta( uN4Re2TPCPionSubPosEta );
  fpQvecEvent->setTPCPion4Im2NegChSubPosEta( uN4Im2TPCPionSubPosEta );
  fpQvecEvent->setTPCPion2Re3NegChSubPosEta( uN2Re3TPCPionSubPosEta );
  fpQvecEvent->setTPCPion2Im3NegChSubPosEta( uN2Im3TPCPionSubPosEta );
  fpQvecEvent->setTPCPion0MNegChSubPosEta( uN0MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion3MNegChSubPosEta( uN3MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion4MNegChSubPosEta( uN4MTPCPionSubPosEta );
  fpQvecEvent->setTPCPion4Re2NegChSubNegEta( uN4Re2TPCPionSubNegEta );
  fpQvecEvent->setTPCPion4Im2NegChSubNegEta( uN4Im2TPCPionSubNegEta );
  fpQvecEvent->setTPCPion2Re3NegChSubNegEta( uN2Re3TPCPionSubNegEta );
  fpQvecEvent->setTPCPion2Im3NegChSubNegEta( uN2Im3TPCPionSubNegEta );
  fpQvecEvent->setTPCPion0MNegChSubNegEta( uN0MTPCPionSubNegEta );
  fpQvecEvent->setTPCPion3MNegChSubNegEta( uN3MTPCPionSubNegEta );
  fpQvecEvent->setTPCPion4MNegChSubNegEta( uN4MTPCPionSubNegEta );
  
  fpQvecEvent->setTPCKaonRePosChSubPosEta( uPReTPCKaonSubPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCKaonImPosChSubPosEta( uPImTPCKaonSubPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCKaon2RePosChSubPosEta( uP2ReTPCKaonSubPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCKaon2ImPosChSubPosEta( uP2ImTPCKaonSubPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCKaon2Re2PosChSubPosEta( uP2Re2TPCKaonSubPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCKaon2Im2PosChSubPosEta( uP2Im2TPCKaonSubPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCKaonMPosChSubPosEta( uPMTPCKaonSubPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCKaon2MPosChSubPosEta( uP2MTPCKaonSubPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCKaonRePosChSubNegEta( uPReTPCKaonSubNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCKaonImPosChSubNegEta( uPImTPCKaonSubNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCKaon2RePosChSubNegEta( uP2ReTPCKaonSubNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCKaon2ImPosChSubNegEta( uP2ImTPCKaonSubNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCKaon2Re2PosChSubNegEta( uP2Re2TPCKaonSubNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCKaon2Im2PosChSubNegEta( uP2Im2TPCKaonSubNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCKaonMPosChSubNegEta( uPMTPCKaonSubNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCKaon2MPosChSubNegEta( uP2MTPCKaonSubNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCKaonReNegChSubPosEta( uNReTPCKaonSubPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCKaonImNegChSubPosEta( uNImTPCKaonSubPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCKaon2ReNegChSubPosEta( uN2ReTPCKaonSubPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCKaon2ImNegChSubPosEta( uN2ImTPCKaonSubPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCKaon2Re2NegChSubPosEta( uN2Re2TPCKaonSubPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCKaon2Im2NegChSubPosEta( uN2Im2TPCKaonSubPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCKaonMNegChSubPosEta( uNMTPCKaonSubPosEta );   // w ch- eta+
  fpQvecEvent->setTPCKaon2MNegChSubPosEta( uN2MTPCKaonSubPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCKaonReNegChSubNegEta( uNReTPCKaonSubNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCKaonImNegChSubNegEta( uNImTPCKaonSubNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCKaon2ReNegChSubNegEta( uN2ReTPCKaonSubNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCKaon2ImNegChSubNegEta( uN2ImTPCKaonSubNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCKaon2Re2NegChSubNegEta( uN2Re2TPCKaonSubNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCKaon2Im2NegChSubNegEta( uN2Im2TPCKaonSubNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCKaonMNegChSubNegEta( uNMTPCKaonSubNegEta );   // w ch- eta-
  fpQvecEvent->setTPCKaon2MNegChSubNegEta( uN2MTPCKaonSubNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCKaon4Re2PosChSubPosEta( uP4Re2TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4Im2PosChSubPosEta( uP4Im2TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon2Re3PosChSubPosEta( uP2Re3TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon2Im3PosChSubPosEta( uP2Im3TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon0MPosChSubPosEta( uP0MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon3MPosChSubPosEta( uP3MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4MPosChSubPosEta( uP4MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4Re2PosChSubNegEta( uP4Re2TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon4Im2PosChSubNegEta( uP4Im2TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon2Re3PosChSubNegEta( uP2Re3TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon2Im3PosChSubNegEta( uP2Im3TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon0MPosChSubNegEta( uP0MTPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon3MPosChSubNegEta( uP3MTPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon4MPosChSubNegEta( uP4MTPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon4Re2NegChSubPosEta( uN4Re2TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4Im2NegChSubPosEta( uN4Im2TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon2Re3NegChSubPosEta( uN2Re3TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon2Im3NegChSubPosEta( uN2Im3TPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon0MNegChSubPosEta( uN0MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon3MNegChSubPosEta( uN3MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4MNegChSubPosEta( uN4MTPCKaonSubPosEta );
  fpQvecEvent->setTPCKaon4Re2NegChSubNegEta( uN4Re2TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon4Im2NegChSubNegEta( uN4Im2TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon2Re3NegChSubNegEta( uN2Re3TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon2Im3NegChSubNegEta( uN2Im3TPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon0MNegChSubNegEta( uN0MTPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon3MNegChSubNegEta( uN3MTPCKaonSubNegEta );
  fpQvecEvent->setTPCKaon4MNegChSubNegEta( uN4MTPCKaonSubNegEta );
  
  fpQvecEvent->setTPCProtonRePosChSubPosEta( uPReTPCProtonSubPosEta ); // w * cos(theta+) eta+
  fpQvecEvent->setTPCProtonImPosChSubPosEta( uPImTPCProtonSubPosEta ); // w * sin(theta+) eta+
  fpQvecEvent->setTPCProton2RePosChSubPosEta( uP2ReTPCProtonSubPosEta ); // w * cos(2theta+) eta+
  fpQvecEvent->setTPCProton2ImPosChSubPosEta( uP2ImTPCProtonSubPosEta ); // w * sin(2theta+) eta+
  fpQvecEvent->setTPCProton2Re2PosChSubPosEta( uP2Re2TPCProtonSubPosEta ); // w^2 * cos(2theta+) eta+
  fpQvecEvent->setTPCProton2Im2PosChSubPosEta( uP2Im2TPCProtonSubPosEta ); // w^2 * sin(2theta+) eta+
  fpQvecEvent->setTPCProtonMPosChSubPosEta( uPMTPCProtonSubPosEta );   // w ch+ eta+
  fpQvecEvent->setTPCProton2MPosChSubPosEta( uP2MTPCProtonSubPosEta );   // w^2 ch+ eta+
  
  fpQvecEvent->setTPCProtonRePosChSubNegEta( uPReTPCProtonSubNegEta ); // w * cos(theta+) eta-
  fpQvecEvent->setTPCProtonImPosChSubNegEta( uPImTPCProtonSubNegEta ); // w * sin(theta+) eta-
  fpQvecEvent->setTPCProton2RePosChSubNegEta( uP2ReTPCProtonSubNegEta ); // w * cos(2theta+) eta-
  fpQvecEvent->setTPCProton2ImPosChSubNegEta( uP2ImTPCProtonSubNegEta ); // w * sin(2theta+) eta-
  fpQvecEvent->setTPCProton2Re2PosChSubNegEta( uP2Re2TPCProtonSubNegEta ); // w^2 * cos(2theta+) eta-
  fpQvecEvent->setTPCProton2Im2PosChSubNegEta( uP2Im2TPCProtonSubNegEta ); // w^2 * sin(2theta+) eta-
  fpQvecEvent->setTPCProtonMPosChSubNegEta( uPMTPCProtonSubNegEta );   // w ch+ eta-
  fpQvecEvent->setTPCProton2MPosChSubNegEta( uP2MTPCProtonSubNegEta );   // w^2 ch+ eta-
  
  fpQvecEvent->setTPCProtonReNegChSubPosEta( uNReTPCProtonSubPosEta ); // w * cos(theta-) eta+
  fpQvecEvent->setTPCProtonImNegChSubPosEta( uNImTPCProtonSubPosEta ); // w * sin(theta-) eta+
  fpQvecEvent->setTPCProton2ReNegChSubPosEta( uN2ReTPCProtonSubPosEta ); // w * cos(2theta-) eta+
  fpQvecEvent->setTPCProton2ImNegChSubPosEta( uN2ImTPCProtonSubPosEta ); // w * sin(2theta-) eta+
  fpQvecEvent->setTPCProton2Re2NegChSubPosEta( uN2Re2TPCProtonSubPosEta ); // w^2 * cos(2theta-) eta+
  fpQvecEvent->setTPCProton2Im2NegChSubPosEta( uN2Im2TPCProtonSubPosEta ); // w^2 * sin(2theta-) eta+
  fpQvecEvent->setTPCProtonMNegChSubPosEta( uNMTPCProtonSubPosEta );   // w ch- eta+
  fpQvecEvent->setTPCProton2MNegChSubPosEta( uN2MTPCProtonSubPosEta );   // w^2  h- eta+
  
  fpQvecEvent->setTPCProtonReNegChSubNegEta( uNReTPCProtonSubNegEta ); // w * cos(theta-) eta-
  fpQvecEvent->setTPCProtonImNegChSubNegEta( uNImTPCProtonSubNegEta ); // w * sin(theta-) eta-
  fpQvecEvent->setTPCProton2ReNegChSubNegEta( uN2ReTPCProtonSubNegEta ); // w * cos(2theta-) eta-
  fpQvecEvent->setTPCProton2ImNegChSubNegEta( uN2ImTPCProtonSubNegEta ); // w * sin(2theta-) eta-
  fpQvecEvent->setTPCProton2Re2NegChSubNegEta( uN2Re2TPCProtonSubNegEta ); // w^2 * cos(2theta-) eta-
  fpQvecEvent->setTPCProton2Im2NegChSubNegEta( uN2Im2TPCProtonSubNegEta ); // w^2 * sin(2theta-) eta-
  fpQvecEvent->setTPCProtonMNegChSubNegEta( uNMTPCProtonSubNegEta );   // w ch- eta-
  fpQvecEvent->setTPCProton2MNegChSubNegEta( uN2MTPCProtonSubNegEta );   // w^2 ch- eta-
  
  fpQvecEvent->setTPCProton4Re2PosChSubPosEta( uP4Re2TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4Im2PosChSubPosEta( uP4Im2TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton2Re3PosChSubPosEta( uP2Re3TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton2Im3PosChSubPosEta( uP2Im3TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton0MPosChSubPosEta( uP0MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton3MPosChSubPosEta( uP3MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4MPosChSubPosEta( uP4MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4Re2PosChSubNegEta( uP4Re2TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton4Im2PosChSubNegEta( uP4Im2TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton2Re3PosChSubNegEta( uP2Re3TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton2Im3PosChSubNegEta( uP2Im3TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton0MPosChSubNegEta( uP0MTPCProtonSubNegEta );
  fpQvecEvent->setTPCProton3MPosChSubNegEta( uP3MTPCProtonSubNegEta );
  fpQvecEvent->setTPCProton4MPosChSubNegEta( uP4MTPCProtonSubNegEta );
  fpQvecEvent->setTPCProton4Re2NegChSubPosEta( uN4Re2TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4Im2NegChSubPosEta( uN4Im2TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton2Re3NegChSubPosEta( uN2Re3TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton2Im3NegChSubPosEta( uN2Im3TPCProtonSubPosEta );
  fpQvecEvent->setTPCProton0MNegChSubPosEta( uN0MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton3MNegChSubPosEta( uN3MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4MNegChSubPosEta( uN4MTPCProtonSubPosEta );
  fpQvecEvent->setTPCProton4Re2NegChSubNegEta( uN4Re2TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton4Im2NegChSubNegEta( uN4Im2TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton2Re3NegChSubNegEta( uN2Re3TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton2Im3NegChSubNegEta( uN2Im3TPCProtonSubNegEta );
  fpQvecEvent->setTPCProton0MNegChSubNegEta( uN0MTPCProtonSubNegEta );
  fpQvecEvent->setTPCProton3MNegChSubNegEta( uN3MTPCProtonSubNegEta );
  fpQvecEvent->setTPCProton4MNegChSubNegEta( uN4MTPCProtonSubNegEta );
  
}




void AliAnalysisTaskGammaDeltaPIDSaveQvec::SetupQAHistograms(){

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



  
  
}







void AliAnalysisTaskGammaDeltaPIDSaveQvec::SetupPileUpRemovalFunctions(){
  
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













void AliAnalysisTaskGammaDeltaPIDSaveQvec::SetupEventAndTaskConfigInfo(){

 
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


void AliAnalysisTaskGammaDeltaPIDSaveQvec::SetupQvecSavingObjects() {
  
  fpQvecEvent = new AliAnalysisTaskGammaDeltaPIDSaveQvecEvent();
  
  fEventList = new TList();
  fEventList->SetName("Event List");
  fEventList->SetOwner(kTRUE);
  fListHist->Add(fEventList);
  
  treeEvent = new TTree("events", "event");
  treeEvent->Branch("event", &fpQvecEvent);
  
  fEventList->Add(treeEvent);
	
}

void AliAnalysisTaskGammaDeltaPIDSaveQvec::GetNUACorrectionHist(Int_t run, Int_t kParticleID) {
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
    
    fHCorrectNUAkPIDPosPion = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID10Pos_Run%d",run)); 
    fHCorrectNUAkPIDNegPion = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID10Neg_Run%d",run));
    fHCorrectNUAkPIDPosKaon = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID20Pos_Run%d",run)); 
    fHCorrectNUAkPIDNegKaon = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID20Neg_Run%d",run));
    fHCorrectNUAkPIDPosProton = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID30Pos_Run%d",run)); 
    fHCorrectNUAkPIDNegProton = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID30Neg_Run%d",run));
    
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

void AliAnalysisTaskGammaDeltaPIDSaveQvec::GetMCCorrectionHist(){

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

void AliAnalysisTaskGammaDeltaPIDSaveQvec::GetV0MCorrectionHist(Int_t run){ 

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


Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvec::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){
  
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
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvec::CheckEventIsPileUp2018(AliAODEvent *faod) {

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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetTPCQvectAndRemovePileUp2018(AliAODEvent *faod,Double_t *qnxEtaNeg,Double_t *qnyEtaNeg,Double_t *qnxEtaPos,Double_t *qnyEtaPos,Double_t& multNeg,Double_t& multPos)
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
      SumQ2xTPCPos   += trk1Wgt*trk1Pt*TMath::Cos(2*trk1Phi); // trk1Wgt no MC ??????
      SumQ2yTPCPos   += trk1Wgt*trk1Pt*TMath::Sin(2*trk1Phi);
      SumQ3xTPCPos   += trk1Wgt*trk1Pt*TMath::Cos(3*trk1Phi);
      SumQ3yTPCPos   += trk1Wgt*trk1Pt*TMath::Sin(3*trk1Phi);
      SumQ4xTPCPos   += trk1Wgt*trk1Pt*TMath::Cos(4*trk1Phi);
      SumQ4yTPCPos   += trk1Wgt*trk1Pt*TMath::Sin(4*trk1Phi);
      
      fWgtMultTPCPos += trk1Wgt; 
    }
    else if(trk1Eta <= fEtaGapNegforEP){
      SumQ2xTPCNeg   += trk1Wgt*trk1Pt*TMath::Cos(2*trk1Phi);
      SumQ2yTPCNeg   += trk1Wgt*trk1Pt*TMath::Sin(2*trk1Phi);
      SumQ3xTPCNeg   += trk1Wgt*trk1Pt*TMath::Cos(3*trk1Phi);
      SumQ3yTPCNeg   += trk1Wgt*trk1Pt*TMath::Sin(3*trk1Phi);
      SumQ4xTPCNeg   += trk1Wgt*trk1Pt*TMath::Cos(4*trk1Phi);
      SumQ4yTPCNeg   += trk1Wgt*trk1Pt*TMath::Sin(4*trk1Phi);      
      
      fWgtMultTPCNeg += trk1Wgt;
    }
        
  }//AOD track loop



  
  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t  multV0a   = aodV0->GetMTotV0A();
  Float_t  multV0c   = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On  = multV0aOn + multV0cOn;

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


Double_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetNUAWeightForTrack(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

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

Double_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetNUAWeightForTrackPID(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

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


Double_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetNUAWeightForTrackPID(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg,Int_t gPID){

  Double_t WgtNUAtrk1 = 1.0;
  Int_t    iBinforNUA = 0;

  /// It is responsibility of User to Set the Correct NUA file for the PID needed. 
  if(gChrg>0){
	if(gPID==1) {// pion
      if(fHCorrectNUAkPIDPosPion){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDPosPion->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDPosPion->GetBinContent(iBinforNUA);
      }
	} else if(gPID==2) {//kaon
	  if(fHCorrectNUAkPIDPosKaon){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDPosKaon->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDPosKaon->GetBinContent(iBinforNUA);
      }
	} else if (gPID==3) {//proton
	  if(fHCorrectNUAkPIDPosProton){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDPosProton->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDPosProton->GetBinContent(iBinforNUA);
      }
	} else {
	  printf("\n\n **WARNING** ::GetNUAWeightForTrackPID() Unknown PID.\n\n");
	}
  }
  else{
    if(gPID==1) {// pion
      if(fHCorrectNUAkPIDNegPion){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDNegPion->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDNegPion->GetBinContent(iBinforNUA);
      }
	} else if(gPID==2) {//kaon
	  if(fHCorrectNUAkPIDNegKaon){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDNegKaon->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDNegKaon->GetBinContent(iBinforNUA);
      }
	} else if (gPID==3) {//proton
	  if(fHCorrectNUAkPIDNegProton){ /// safety measures for breaks!
        iBinforNUA = fHCorrectNUAkPIDNegProton->FindBin(fVtxZ,fPhi,fEta);
        WgtNUAtrk1 = fHCorrectNUAkPIDNegProton->GetBinContent(iBinforNUA);
      }
	} else {
	  printf("\n\n **WARNING** ::GetNUAWeightForTrackPID() Unknown PID.\n\n");
	}
  }
  
  return WgtNUAtrk1;
}


void AliAnalysisTaskGammaDeltaPIDSaveQvec::ApplyTPCqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t& qxEtaNeg, Double_t& qyEtaNeg,Double_t& qxEtaPos,Double_t& qyEtaPos){

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

void  AliAnalysisTaskGammaDeltaPIDSaveQvec::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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





Double_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetMCEfficiencyWeightForTrack(Double_t fPt,Int_t gChrg,Int_t kPID){
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
	


Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvec::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A, Double_t &sumMultV0C, Double_t &sumMultV0A){

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
    sumMultV0C = fSumMV0C;
    sumMultV0A = fSumMV0A;
    return kTRUE;  
  }
  
}

void AliAnalysisTaskGammaDeltaPIDSaveQvec::Terminate(Option_t *)  {
  //fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  //if (!fOutputList) return;
}





