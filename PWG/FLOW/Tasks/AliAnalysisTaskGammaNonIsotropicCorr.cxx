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
/* $Id: AliAnalysisTaskGammaNonIsotropicCorr.cxx ver: 2.0                     $   */
/* Simple Task to fill V0 and ZDC Energies for Gain Calibration           */
/* Works with 15o and 18q/r. Support to be added for LHC10h               */
/* Developer: Md Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)        */
/* Last Modified: Aug 23, 2021,  First version committed                  */
/* Last Modified: Oct 08, 2021,  Second version committed                 */
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskGammaNonIsotropicCorr.h"
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

ClassImp(AliAnalysisTaskGammaNonIsotropicCorr)

AliAnalysisTaskGammaNonIsotropicCorr::AliAnalysisTaskGammaNonIsotropicCorr(const char *name):
  AliAnalysisTaskSE(name),
  whichData(0),
  period("0"),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fTempList(NULL),
  fListV0MCorr(NULL),   
  fListV0MCorrRunByRun(NULL),
  fListZDCCorr(NULL),
  fListZDCCorrRunByRun(NULL),
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
  fTPCsharedCut(56), // 70*0.8
  bUseTPCCrossedRows(kFALSE),
  fMinVzCut(-10.0),
  fMaxVzCut(+10.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(-1), // 2.4
  fDCAzzMax(-1), // 3.2
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
  
  //fHistTPCPsiNPosPlane(NULL),
  //fHistTPCPsiNNegPlane(NULL),
  //fHistTPCPsi3PosPlane(NULL),
  //fHistTPCPsi3NegPlane(NULL),
  //fHistTPCPsi4PosPlane(NULL),
  //fHistTPCPsi4NegPlane(NULL),
  
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  //fHistTPCPosqVectorvsCent(NULL),
  //fHistTPCNegqVectorvsCent(NULL),
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
  fHV0Cparameters(NULL),
  fHV0Aparameters(NULL),
  fHZDCCparameters(NULL),
  fHZDCAparameters(NULL),
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hTPCPsi4Correlation(NULL),  
  hV0CV0APsiNCorrelation(NULL),
  //hV0CTPCPsiNCorrelation(NULL),
  //hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  //hV0CTPCPsi3Correlation(NULL),
  //hV0ATPCPsi3Correlation(NULL),  
  
  //fAvgCos2PsivsCentEtaPos(NULL),
  //fAvgSin2PsivsCentEtaPos(NULL),
  //fAvgCos2PsivsCentEtaNeg(NULL),
  //fAvgSin2PsivsCentEtaNeg(NULL),
  //fAvgCos3PsivsCentEtaPos(NULL),
  //fAvgSin3PsivsCentEtaPos(NULL),
  //fAvgCos3PsivsCentEtaNeg(NULL),
  //fAvgSin3PsivsCentEtaNeg(NULL),
  //fAvgCos4PsivsCentEtaPos(NULL),
  //fAvgSin4PsivsCentEtaPos(NULL),
  //fAvgCos4PsivsCentEtaNeg(NULL),
  //fAvgSin4PsivsCentEtaNeg(NULL),
  hAvgV0ChannelsvsVz(NULL)
{
  centrality = -99.0;
  fPsi2V0C = 0.;
  fPsi2V0A = 0.;
  fPsiZDCC = 0.;
  fPsiZDCA = 0.;
  fPsiZDCCA = 0; 
  fCMEQReRP = NULL;
  fCMEQImRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent = NULL;
  f2pCorrelatorCosPsiDiffPerEvent = NULL;
  f2pCorrelatorCos2PsiDiffPerEvent = NULL;
  fCMEQRePOIPos = NULL;
  fCMEQImPOIPos = NULL;
  fCMEQRePOINeg = NULL;
  fCMEQImPOINeg = NULL;
  fNITV0OS = NULL;
  fNITZDCOS = NULL;
  fNITV0POIPos = NULL;
  fNITZDCPOIPos = NULL;
  fNITV0POINeg = NULL;
  fNITZDCPOINeg = NULL;

  fRePEBEOS = NULL;
  fImPEBEOS = NULL;
  f2pCorrelatorCosPsiDiffOS = NULL;
  f2pCorrelatorCos2PsiDiffOS = NULL;
  fRePEBEPP = NULL;
  fImPEBEPP = NULL;
  f2pCorrelatorCosPsiDiffPP = NULL;
  f2pCorrelatorCos2PsiDiffPP = NULL;
  fRePEBENN = NULL;
  fImPEBENN = NULL;
  f2pCorrelatorCosPsiDiffNN = NULL;
  f2pCorrelatorCos2PsiDiffNN = NULL;
  
  // 2-p correlator
  fCMESPPPCenBin = 18;
  f2pCorrelatorCos2PsiDiff2PsiV0CRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0ARP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCARP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCARP = NULL;
  
  f2pCorrelatorCosPsiDiff = NULL;
  f2pCorrelatorCos2PsiDiff = NULL;
  
  f2pCorrelatorCosPsiDiffEBE = NULL;
  f2pCorrelatorCos2PsiDiffEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE = NULL;
  
  // 3-p correlator
  f3pCorrelatorTPCOSPro = NULL; // cos[n(psi1+psi2-2phi3)] OS
  f3pCorrelatorV0COSPro = NULL;
  f3pCorrelatorV0AOSPro = NULL;
  f3pCorrelatorZDCCOSPro = NULL;
  f3pCorrelatorZDCAOSPro = NULL;
  f3pCorrelatorZDCCAOSPro = NULL;
  
  f3pCorrelatorTPCPPPro = NULL; // cos[n(psi1+psi2-2phi3)] PP
  f3pCorrelatorV0CPPPro = NULL;
  f3pCorrelatorV0APPPro = NULL;
  f3pCorrelatorZDCCPPPro = NULL;
  f3pCorrelatorZDCAPPPro = NULL;
  f3pCorrelatorZDCCAPPPro = NULL;
  
  f3pCorrelatorTPCNNPro = NULL; // cos[n(psi1+psi2-2phi3)] NN
  f3pCorrelatorV0CNNPro = NULL;
  f3pCorrelatorV0ANNPro = NULL;
  f3pCorrelatorZDCCNNPro = NULL;
  f3pCorrelatorZDCANNPro = NULL;
  f3pCorrelatorZDCCANNPro = NULL;
  
  f3pCorrelatorTPCSSPro = NULL; // cos[n(psi1+psi2-2phi3)] SS
  f3pCorrelatorV0CSSPro = NULL;
  f3pCorrelatorV0ASSPro = NULL;
  f3pCorrelatorZDCCSSPro = NULL;
  f3pCorrelatorZDCASSPro = NULL;
  f3pCorrelatorZDCCASSPro = NULL;
  
  // non-isotropic terms 3-p 
  for(Int_t i = 0; i < 10; i++) {
    fNonIsotropicTermsPro[i] = NULL;
    fNonIsotropicTermsOSPro[i] = NULL;
    fNonIsotropicTermsPPPro[i] = NULL;
    fNonIsotropicTermsNNPro[i] = NULL;
    fNonIsotropicTermsSSPro[i] = NULL;
  }
	
  for(Int_t i = 0; i < 8; i++) {
    fNonIsotropicTermsV0Pro[i] = NULL;
    fNonIsotropicTermsV0OSPro[i] = NULL;
    fNonIsotropicTermsV0PPPro[i] = NULL;
    fNonIsotropicTermsV0NNPro[i] = NULL;
    fNonIsotropicTermsV0SSPro[i] = NULL;
  }
	
  for(Int_t i = 0; i < 12; i++) {
    fNonIsotropicTermsZDCPro[i] = NULL;
    fNonIsotropicTermsZDCOSPro[i] = NULL;
    fNonIsotropicTermsZDCPPPro[i] = NULL;
    fNonIsotropicTermsZDCNNPro[i] = NULL;
    fNonIsotropicTermsZDCSSPro[i] = NULL;
  }
	
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaNonIsotropicCorr::AliAnalysisTaskGammaNonIsotropicCorr():
  AliAnalysisTaskSE(),
  whichData(0),
  period("0"),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fTempList(NULL),
  fListV0MCorr(NULL),   
  fListV0MCorrRunByRun(NULL),
  fListZDCCorr(NULL),
  fListZDCCorrRunByRun(NULL),
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
  fTPCsharedCut(56), // 70*0.8
  bUseTPCCrossedRows(kFALSE),
  fMinVzCut(-10.0),
  fMaxVzCut(+10.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(-1),
  fDCAzzMax(-1),
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

  
  //fHistTPCPsiNPosPlane(NULL),
  //fHistTPCPsiNNegPlane(NULL),
  //fHistTPCPsi3PosPlane(NULL),
  //fHistTPCPsi3NegPlane(NULL),
  //fHistTPCPsi4PosPlane(NULL),
  //fHistTPCPsi4NegPlane(NULL),
  
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  //fHistTPCPosqVectorvsCent(NULL),
  //fHistTPCNegqVectorvsCent(NULL),
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
  fHV0Cparameters(NULL),
  fHV0Aparameters(NULL),
  fHZDCCparameters(NULL),
  fHZDCAparameters(NULL),
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hTPCPsi4Correlation(NULL),    
  hV0CV0APsiNCorrelation(NULL),
  //hV0CTPCPsiNCorrelation(NULL),
  //hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  //hV0CTPCPsi3Correlation(NULL),
  //hV0ATPCPsi3Correlation(NULL),  
  
  //fAvgCos2PsivsCentEtaPos(NULL),
  //fAvgSin2PsivsCentEtaPos(NULL),
  //fAvgCos2PsivsCentEtaNeg(NULL),
  //fAvgSin2PsivsCentEtaNeg(NULL),
  //fAvgCos3PsivsCentEtaPos(NULL),
  //fAvgSin3PsivsCentEtaPos(NULL),
  //fAvgCos3PsivsCentEtaNeg(NULL),
  //fAvgSin3PsivsCentEtaNeg(NULL),
  //fAvgCos4PsivsCentEtaPos(NULL),
  //fAvgSin4PsivsCentEtaPos(NULL),
  //fAvgCos4PsivsCentEtaNeg(NULL),
  //fAvgSin4PsivsCentEtaNeg(NULL),  
  hAvgV0ChannelsvsVz(NULL)
{
  centrality = -99.0;
  fPsi2V0C = 0.;
  fPsi2V0A = 0.;
  fPsiZDCC = 0.;
  fPsiZDCA = 0.;
  fPsiZDCCA = 0; 
  fCMEQReRP = NULL;
  fCMEQImRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent = NULL;
  f2pCorrelatorCosPsiDiffPerEvent = NULL;
  f2pCorrelatorCos2PsiDiffPerEvent = NULL;
  fCMEQRePOIPos = NULL;
  fCMEQImPOIPos = NULL;
  fCMEQRePOINeg = NULL;
  fCMEQImPOINeg = NULL;
  fNITV0OS = NULL;
  fNITZDCOS = NULL;
  fNITV0POIPos = NULL;
  fNITZDCPOIPos = NULL;
  fNITV0POINeg = NULL;
  fNITZDCPOINeg = NULL;

  fRePEBEOS = NULL;
  fImPEBEOS = NULL;
  f2pCorrelatorCosPsiDiffOS = NULL;
  f2pCorrelatorCos2PsiDiffOS = NULL;
  fRePEBEPP = NULL;
  fImPEBEPP = NULL;
  f2pCorrelatorCosPsiDiffPP = NULL;
  f2pCorrelatorCos2PsiDiffPP = NULL;
  fRePEBENN = NULL;
  fImPEBENN = NULL;
  f2pCorrelatorCosPsiDiffNN = NULL;
  f2pCorrelatorCos2PsiDiffNN = NULL;
  
  // 2-p correlator
  fCMESPPPCenBin = 18;
  f2pCorrelatorCos2PsiDiff2PsiV0CRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0ARP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCRP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCARP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCARP = NULL;
  
  f2pCorrelatorCosPsiDiff = NULL;
  f2pCorrelatorCos2PsiDiff = NULL;
  
  f2pCorrelatorCosPsiDiffEBE = NULL;
  f2pCorrelatorCos2PsiDiffEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE = NULL;
  
  // 3-p correlator
  f3pCorrelatorTPCOSPro = NULL; // cos[n(psi1+psi2-2phi3)] OS
  f3pCorrelatorV0COSPro = NULL;
  f3pCorrelatorV0AOSPro = NULL;
  f3pCorrelatorZDCCOSPro = NULL;
  f3pCorrelatorZDCAOSPro = NULL;
  f3pCorrelatorZDCCAOSPro = NULL;
  
  f3pCorrelatorTPCPPPro = NULL; // cos[n(psi1+psi2-2phi3)] PP
  f3pCorrelatorV0CPPPro = NULL;
  f3pCorrelatorV0APPPro = NULL;
  f3pCorrelatorZDCCPPPro = NULL;
  f3pCorrelatorZDCAPPPro = NULL;
  f3pCorrelatorZDCCAPPPro = NULL;
  
  f3pCorrelatorTPCNNPro = NULL; // cos[n(psi1+psi2-2phi3)] NN
  f3pCorrelatorV0CNNPro = NULL;
  f3pCorrelatorV0ANNPro = NULL;
  f3pCorrelatorZDCCNNPro = NULL;
  f3pCorrelatorZDCANNPro = NULL;
  f3pCorrelatorZDCCANNPro = NULL;
  
  f3pCorrelatorTPCSSPro = NULL; // cos[n(psi1+psi2-2phi3)] SS
  f3pCorrelatorV0CSSPro = NULL;
  f3pCorrelatorV0ASSPro = NULL;
  f3pCorrelatorZDCCSSPro = NULL;
  f3pCorrelatorZDCASSPro = NULL;
  f3pCorrelatorZDCCASSPro = NULL;
    
  // non-isotropic terms 3-p 
  for(Int_t i = 0; i < 10; i++) {
    fNonIsotropicTermsPro[i] = NULL;
    fNonIsotropicTermsOSPro[i] = NULL;
    fNonIsotropicTermsPPPro[i] = NULL;
    fNonIsotropicTermsNNPro[i] = NULL;
    fNonIsotropicTermsSSPro[i] = NULL;
  }
	
  for(Int_t i = 0; i < 8; i++) {
    fNonIsotropicTermsV0Pro[i] = NULL;
    fNonIsotropicTermsV0OSPro[i] = NULL;
    fNonIsotropicTermsV0PPPro[i] = NULL;
    fNonIsotropicTermsV0NNPro[i] = NULL;
    fNonIsotropicTermsV0SSPro[i] = NULL;
  }
	
  for(Int_t i = 0; i < 12; i++) {
    fNonIsotropicTermsZDCPro[i] = NULL;
    fNonIsotropicTermsZDCOSPro[i] = NULL;
    fNonIsotropicTermsZDCPPPro[i] = NULL;
    fNonIsotropicTermsZDCNNPro[i] = NULL;
    fNonIsotropicTermsZDCSSPro[i] = NULL;
  }
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskGammaNonIsotropicCorr::~AliAnalysisTaskGammaNonIsotropicCorr()
{
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!
  if(fTempList)      delete fTempList;
  if(fListHist)      delete fListHist;  
  if(fListV0MCorr)   delete fListV0MCorr;
  if(fListV0MCorrRunByRun)   delete fListV0MCorrRunByRun;
  if(fListZDCCorr)   delete fListZDCCorr;
  if(fListZDCCorrRunByRun)   delete fListZDCCorrRunByRun;

  if(fV0CutPU)      delete fV0CutPU;
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU;   
}


//________________ Define Histograms _______________
void AliAnalysisTaskGammaNonIsotropicCorr::UserCreateOutputObjects()
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
  if (whichData == 2018 && period == 'q')
	SetupPileUpRemovalFunctions18qPass3();
  else if (whichData == 2018 && period == 'r')
    SetupPileUpRemovalFunctions18rPass3();
  //SetupPileUpRemovalFunctions();
  SetupEventAndTaskConfigInfo();


 
    
  cout<<"Ncls: "<<fTPCclustMin<<" TPCsharedCut: "<<fTPCsharedCut<<" UseTPCCrossedRows: "<<bUseTPCCrossedRows<<" Harm: "<<gHarmonic<<" POI: "<<gParticleID<<" nsigTPC: "<<fNSigmaTPCCut<<" nsigCirc: "<<fNSigmaTOFCut<<endl;
  cout<<"FB: "<<fFilterBit<<" chi2min: "<<fTrkChi2Min<<" chi2max: "<<fTrkChi2Max<<" etaMin: "<<fMinEtaCut<<" etaMax: "<<fMaxEtaCut<<endl;
  cout<<"dEdxMin: "<<fTPCdEdxMin<<" dcaXY: "<<fDCAxyMax<<" dcaZ: "<<fDCAzzMax<<"  VzLow: "<<fMinVzCut<<" VzHigh: "<<fMaxVzCut<<endl;
  cout<<"minPt: "<<fMinPtCut<<" maxPt: "<<fMaxPtCut<<" etaGapNeg: "<<fEtaGapNeg<<" etaGapPos: "<<fEtaGapPos<<" oldRun: "<<gOldRunNumber<<endl;
  cout<<"Centrality Estimator: "<<sCentrEstimator.Data()<<" EventPlane Detector:"<<sDetectorForEP.Data()<<endl;
  
  //cout<<"\n ******  Bug Testing Mode.. So we exit here ****** \n"<<endl; exit(111);
  
     
  PostData(1,fListHist);
}










//____________________________ Call Event by Event ___________________________________
void AliAnalysisTaskGammaNonIsotropicCorr::UserExec(Option_t*) {
 
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
  

  

  centrality = -99.0;
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
    
  Double_t fOrbitNumber = static_cast<double>(orbit)/1000000.;
  
  Bool_t kPileupEvent = kFALSE;

  //Double_t fSumQnxNeg = 0, fSumQnyNeg = 0, fSumQnxPos = 0, fSumQnyPos = 0;
  //Double_t fSumQnxNeg[3] = {0,}; // Array: Q2x, Q3x, Q4x,...
  //Double_t fSumQnyNeg[3] = {0,}; // Array: Q2y, Q3y, Q4y,... 
  //Double_t fSumQnxPos[3] = {0,};
  //Double_t fSumQnyPos[3] = {0,};  
  
  //Double_t fMultNeg = 0, fMultPos = 0;
  if (!bSkipPileUpCut) {
  kPileupEvent = CheckEventIsPileUp2018(fAOD);

  if(kPileupEvent) return;  // If not a PileUp event, then We have TPC q vectors for EP.
  }
  fDebugwEventCount->Fill(4.1);

  
  //if(fMultNeg<0.1 || fMultPos<0.1) return;  //// This means there is not enough track in the the event. 
  
  
  //Double_t fQ2xNeg = 0., fQ2yNeg = 0., fQ2xPos = 0., fQ2yPos = 0.;  ///q = Q/Mult Vector from TPC;
  //Double_t fQ3xNeg = 0., fQ3yNeg = 0., fQ3xPos = 0., fQ3yPos = 0.;  ///q = Q/Mult Vector from TPC;
  //Double_t fQ4xNeg = 0., fQ4yNeg = 0., fQ4xPos = 0., fQ4yPos = 0.;  ///q = Q/Mult Vector from TPC;

  //fQ2xNeg = fSumQnxNeg[0]/fMultNeg;
  //fQ2yNeg = fSumQnyNeg[0]/fMultNeg;
  //fQ3xNeg = fSumQnxNeg[1]/fMultNeg;
  //fQ3yNeg = fSumQnyNeg[1]/fMultNeg;
  //fQ4xNeg = fSumQnxNeg[2]/fMultNeg;
  //fQ4yNeg = fSumQnyNeg[2]/fMultNeg;
  
  //fQ2xPos = fSumQnxPos[0]/fMultPos;
  //fQ2yPos = fSumQnyPos[0]/fMultPos;
  //fQ3xPos = fSumQnxPos[1]/fMultPos;
  //fQ3yPos = fSumQnyPos[1]/fMultPos;
  //fQ4xPos = fSumQnxPos[2]/fMultPos;
  //fQ4yPos = fSumQnyPos[2]/fMultPos;


  //cout<<"Before rec: q2xn "<<fQ2xNeg<<"\t q2yn "<<fQ2yNeg<<"\t q2xp "<<fQ2xPos<<"\t q2yp "<<fQ2yPos<<endl;

  // *** Rihan: Temporarily Turned off. Uncomment after end of Test.!
  //ApplyTPCqVectRecenter(centrV0M, 2, fQ2xNeg, fQ2yNeg, fQ2xPos, fQ2yPos);
  //ApplyTPCqVectRecenter(centrV0M, 3, fQ3xNeg, fQ3yNeg, fQ3xPos, fQ3yPos);  

  //cout<<"After  rec: q2xn "<<fQ2xNeg<<"\t q2yn "<<fQ2yNeg<<"\t q2xp "<<fQ2xPos<<"\t q2yp "<<fQ2yPos<<endl;
  //cout<<"------- Bug Testing Mode... we exit here....... "<<endl; return;



  
  /// Fill <Q> vector for TPC event plane :
  //fAvgCos2PsivsCentEtaPos->Fill(centrV0M,fQ2xPos); 
  //fAvgSin2PsivsCentEtaPos->Fill(centrV0M,fQ2yPos); 
  //fAvgCos2PsivsCentEtaNeg->Fill(centrV0M,fQ2xNeg);
  //fAvgSin2PsivsCentEtaNeg->Fill(centrV0M,fQ2yNeg);

  //fAvgCos3PsivsCentEtaPos->Fill(centrV0M,fQ3xPos); 
  //fAvgSin3PsivsCentEtaPos->Fill(centrV0M,fQ3yPos); 
  //fAvgCos3PsivsCentEtaNeg->Fill(centrV0M,fQ3xNeg);
  //fAvgSin3PsivsCentEtaNeg->Fill(centrV0M,fQ3yNeg);

  //fAvgCos4PsivsCentEtaPos->Fill(centrV0M,fQ4xPos); 
  //fAvgSin4PsivsCentEtaPos->Fill(centrV0M,fQ4yPos); 
  //fAvgCos4PsivsCentEtaNeg->Fill(centrV0M,fQ4xNeg);
  //fAvgSin4PsivsCentEtaNeg->Fill(centrV0M,fQ4yNeg);  

  
  
  /// TPC q2-vectors for ESE:
  //fHistTPCPosqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQ2xPos*fQ2xPos + fQ2yPos*fQ2yPos)); // centrV0M hardcoded to prevent mismatch!
  //fHistTPCNegqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQ2xNeg*fQ2xNeg + fQ2yNeg*fQ2yNeg));
  


  /// TPC Event Planes:
  //Double_t fPsiNTPCPos = 0., fPsiNTPCNeg = 0.;
  //Double_t fPsi3TPCPos = 0., fPsi3TPCNeg = 0.;
  //Double_t fPsi4TPCPos = 0., fPsi4TPCNeg = 0.;
  
  //if(fQ2xPos != 0 && fQ2yPos != 0){
    //fPsiNTPCPos = (1./2)*TMath::ATan2(fQ2yPos,fQ2xPos);
    //if(fPsiNTPCPos < 0) fPsiNTPCPos += TMath::TwoPi()/2;
    //fPsi3TPCPos = (1./3)*TMath::ATan2(fQ3yPos,fQ3xPos);
    //if(fPsi3TPCPos < 0) fPsi3TPCPos += TMath::TwoPi()/3;
    //fPsi4TPCPos = (1./4)*TMath::ATan2(fQ4yPos,fQ4xPos);
    //if(fPsi4TPCPos < 0) fPsi4TPCPos += TMath::TwoPi()/4;      
  //}
  //if(fQ2xNeg != 0 && fQ2yNeg != 0){
    //fPsiNTPCNeg = (1./2)*TMath::ATan2(fQ2yNeg,fQ2xNeg);
    //if(fPsiNTPCNeg < 0) fPsiNTPCNeg += TMath::TwoPi()/2;
    //fPsi3TPCNeg = (1./3)*TMath::ATan2(fQ3yNeg,fQ3xNeg);
    //if(fPsi3TPCNeg < 0) fPsi3TPCNeg += TMath::TwoPi()/3;
    //fPsi4TPCNeg = (1./4)*TMath::ATan2(fQ4yNeg,fQ4xNeg);
    //if(fPsi4TPCNeg < 0) fPsi4TPCNeg += TMath::TwoPi()/4;     
  //}

  //fHistTPCPsiNPosPlane->Fill(centrality,fPsiNTPCPos);
  //fHistTPCPsiNNegPlane->Fill(centrality,fPsiNTPCNeg);
  //fHistTPCPsi3PosPlane->Fill(centrality,fPsi3TPCPos);
  //fHistTPCPsi3NegPlane->Fill(centrality,fPsi3TPCNeg); 
  //fHistTPCPsi4PosPlane->Fill(centrality,fPsi4TPCPos);
  //fHistTPCPsi4NegPlane->Fill(centrality,fPsi4TPCNeg); 



  //hTPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsiNTPCNeg));    /// TPC Psi2 Resolution
  //hTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsi3TPCPos - 3*fPsi3TPCNeg));    /// TPC Psi3 Resolution
  //hTPCPsi4Correlation->Fill(centrality,TMath::Cos(4*fPsi3TPCPos - 4*fPsi3TPCNeg));    /// TPC Psi3 Resolution



  
  
  ///==========>  Get V0 Event Planes <===============
     
  Double_t fQnxV0C=0, fQnyV0C=0, fSumMV0C = 0, fQnxV0A=0, fQnyV0A=0, fSumMV0A = 0; 
  Bool_t kPassV0 = GetGainCorrectedV0Qvector(fAOD, fVertexZEvent, 2, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A, fSumMV0C, fSumMV0A);


  if(!kPassV0) return;           /// V0 does not have signal for this event.  
  fDebugwEventCount->Fill(5.1);


  ApplyV0XqVectRecenter(centrCL1, 2, fQnxV0C, fQnyV0C, fQnxV0A, fQnyV0A);

  
  ////---- fill the <Q> vector from V0A/C vs Cent:-----------
  hAvgQNXvsCentV0C->Fill(centrCL1,fQnxV0C);   /// to Avoid self correlation, V0 <Q> is filled with CL1 centrality! 
  hAvgQNYvsCentV0C->Fill(centrCL1,fQnyV0C);  
  hAvgQNXvsCentV0A->Fill(centrCL1,fQnxV0A);
  hAvgQNYvsCentV0A->Fill(centrCL1,fQnyV0A);
  
  //cout<<"------- Bug Testing Mode... we exit here....... "<<endl; return;

  


  /// V0 q2-vectors for ESE:
  fHistV0CDetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0C*fQnxV0C + fQnyV0C*fQnyV0C)); // centrCL1 hardcoded to prevent mismatch!
  fHistV0ADetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0A*fQnxV0A + fQnyV0A*fQnyV0A));
  
  
  ///--------> Get V0A and V0C Event Planes 
    
  fPsi2V0C = (1./gHarmonic)*TMath::ATan2(fQnyV0C,fQnxV0C);
  if(fPsi2V0C < 0) fPsi2V0C += TMath::TwoPi()/gHarmonic;
  fPsi2V0A = (1./gHarmonic)*TMath::ATan2(fQnyV0A,fQnxV0A);
  if(fPsi2V0A < 0) fPsi2V0A += TMath::TwoPi()/gHarmonic;
  fHistV0CPsiNEventPlane->Fill(centrality,fPsi2V0C);
  fHistV0APsiNEventPlane->Fill(centrality,fPsi2V0A);  
 


  /// V0A, V0C Resolutions:
  hV0CV0APsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsi2V0A    - 2*fPsi2V0C));
  //hV0CTPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsi2V0C));
  //hV0ATPCPsiNCorrelation->Fill(centrality,TMath::Cos(2*fPsiNTPCPos - 2*fPsi2V0A));

  hV0CV0APsi3Correlation->Fill(centrality,TMath::Cos(3*fPsi2V0A    - 3*fPsi2V0C));
  //hV0CTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsi2V0C));
  //hV0ATPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsi2V0A));


  
  //=============== Get the ZDC data ==================== @Shi

  Double_t fQxZNCC=0, fQyZNCC=0, fQxZNCA=0, fQyZNCA=0; 
  
  
  AliAODZDC *aodZDC = fAOD->GetZDCData();
	
  if(!aodZDC) {
    printf("\n ********* Error: could not find ZDC data ************ \n ");
  }
  else if(aodZDC){
    const Double_t *fZNATowerRawAOD = aodZDC->GetZNATowerEnergy();
    const Double_t *fZNCTowerRawAOD = aodZDC->GetZNCTowerEnergy();
	
	if((fZNATowerRawAOD[0]<0) || (fZNATowerRawAOD[1]<0) || (fZNATowerRawAOD[2]<0) || (fZNATowerRawAOD[3]<0) || (fZNATowerRawAOD[4] < 0)) {
		return;
	}
	
	if((fZNCTowerRawAOD[0]<0) || (fZNCTowerRawAOD[1]<0) || (fZNCTowerRawAOD[2]<0) || (fZNCTowerRawAOD[3]<0) || (fZNCTowerRawAOD[4] < 0)) {
		return;
	}
	fDebugwEventCount->Fill(6.1);
	
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
	
	
	if (denZNC==0) {return;}
	if (denZNA==0) {return;}
	
	fDebugwEventCount->Fill(7.1);
	
	Double_t ZDCCxPosFromLogWeight = numXZNC/denZNC;
	Double_t ZDCCyPosFromLogWeight = numYZNC/denZNC;
	Double_t ZDCAxPosFromLogWeight = numXZNA/denZNA;
	Double_t ZDCAyPosFromLogWeight = numYZNA/denZNA;
	
	Double_t ZDCCAvgxPosFromVtxFit = 0;
	Double_t ZDCCAvgyPosFromVtxFit = 0;
	
	Double_t ZDCAAvgxPosFromVtxFit = 0;
	Double_t ZDCAAvgyPosFromVtxFit = 0;

	// 3rd order centrality+vtxpos+orbitNum
	// hard code centrality to centrV0M since it is generated based on centrV0M
	ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrV0M + fHZDCCparameters->GetBinContent(7)*pow(centrV0M,2) + fHZDCCparameters->GetBinContent(8)*pow(centrV0M,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12)*fOrbitNumber + fHZDCCparameters->GetBinContent(13);
	ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(14)*centrV0M + fHZDCCparameters->GetBinContent(15)*pow(centrV0M,2) + fHZDCCparameters->GetBinContent(16)*pow(centrV0M,3) + fHZDCCparameters->GetBinContent(17)*pVtxX + fHZDCCparameters->GetBinContent(18)*pVtxY + fHZDCCparameters->GetBinContent(19)*pVtxZ + fHZDCCparameters->GetBinContent(20)*fOrbitNumber + fHZDCCparameters->GetBinContent(21);
	
	ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrV0M + fHZDCAparameters->GetBinContent(7)*pow(centrV0M,2) + fHZDCAparameters->GetBinContent(8)*pow(centrV0M,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12)*fOrbitNumber + fHZDCAparameters->GetBinContent(13);
	ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(14)*centrV0M + fHZDCAparameters->GetBinContent(15)*pow(centrV0M,2) + fHZDCAparameters->GetBinContent(16)*pow(centrV0M,3) + fHZDCAparameters->GetBinContent(17)*pVtxX + fHZDCAparameters->GetBinContent(18)*pVtxY + fHZDCAparameters->GetBinContent(19)*pVtxZ + fHZDCAparameters->GetBinContent(20)*fOrbitNumber + fHZDCAparameters->GetBinContent(21);
	
	
	fQxZNCC = ZDCCxPosFromLogWeight - ZDCCAvgxPosFromVtxFit;
	fQyZNCC = ZDCCyPosFromLogWeight - ZDCCAvgyPosFromVtxFit;
	
	fQxZNCA = ZDCAxPosFromLogWeight - ZDCAAvgxPosFromVtxFit;
	fQyZNCA = ZDCAyPosFromLogWeight - ZDCAAvgyPosFromVtxFit;
	
	// Event plane

	fPsiZDCA = TMath::ATan2(fQyZNCA,fQxZNCA); // Psi_{1,A} spectator plane -pi to pi
	if (fPsiZDCA < 0) { // Psi_{1,A} should be differ to Psi_{1,C} by pi. 
	  fPsiZDCA = fPsiZDCA + TMath::Pi();
	} else if (fPsiZDCA >= 0) {
	  fPsiZDCA = fPsiZDCA - TMath::Pi();
	}

	fPsiZDCC = TMath::ATan2(fQyZNCC,fQxZNCC); // Psi_{1,C} spectator plane 

	fPsiZDCCA = TMath::ATan2((fQyZNCC-fQyZNCA),(fQxZNCC-fQxZNCA));
  }
  
  /// We Are About to Start Main Analysis Below: //
  //Int_t kPIDtrk1=gParticleID;   /// gParticleID is Set From AddTask. Both Identified..
  //Int_t kPIDtrk2=gParticleID;   /// 0 = hadron (h-h), 1 = Pi-Pi, 2 = K-K, 3 = Prot-Prot, 

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
 
  
  Double_t trk1Pt=0,trk1Phi=0,trk1Eta=0,trk1DCAxy=0.0, trk1DCAz=0.0,trk1Chi2=0,trk1dEdx=0,trk1Wgt=1.0;
  Double_t trk2Pt=0,trk2Phi=0,trk2Eta=0,trk2DCAxy=0.0, trk2DCAz=0.0,trk2Chi2=0,trk2dEdx=0,trk2Wgt=1.0;

  Double_t posTrk1[3] = {0 , 0, 0};
  Double_t posTrk2[3] = {0 , 0, 0};
  //Double_t vTrk1[3] = {0 , 0, 0};
  //Double_t vTrk2[3] = {0 , 0, 0};
  

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
	  //trk1DCAxy = AODtrack1->DCA(); // do not use DCA for FB 768
	  //trk1DCAz  = AODtrack1->ZAtDCA();            
	  trk1dEdx  = AODtrack1->GetDetPid()->GetTPCsignal();  
	  
	  
	  // for constrained TPConly tracks
      if(fFilterBit == 128){
	  	trk1DCAxy = AODtrack1->DCA();      // this is the DCA from global track (not exactly what is cut on)
	    trk1DCAz  = AODtrack1->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
      }
      else{
	    //const AliVVertex *vertex = fVevent->GetPrimaryVertex();
	    //vertex->GetXYZ(vTrk1);
	    AODtrack1->GetXYZ(posTrk1);
	    
		
	    trk1DCAxy  = TMath::Sqrt((posTrk1[0] - pVtxX)*(posTrk1[0] - pVtxX) + (posTrk1[1] - pVtxY)*(posTrk1[1] - pVtxY));
	    trk1DCAz   = posTrk1[2] - pVtxZ;
      }
      
      // DCA cut
      if(fDCAxyMax>0 && fDCAzzMax>0){
		if(TMath::Sqrt((trk1DCAxy*trk1DCAxy)/(fDCAxyMax*fDCAxyMax)+(trk1DCAz*trk1DCAz)/(fDCAzzMax*fDCAzzMax)) > 1 ){
		  continue;  // 2D cut
		}
      }

      //if (fDCAxyMax>0 && fDCAzzMax>0 && trk1DCAxy>=fDCAxyMax && trk1DCAz>=fDCAzzMax) { // if fDCAxyMax, fDCAzzMax is set to be less than 0, no cut applied
		//continue;
	  //}
	  
	  // Crossed Row cut
	  if (bUseTPCCrossedRows){ // either crossed row or ncluster Min should be applied
		if ((Float_t)AODtrack1->GetTPCNCrossedRows() < (120 - (5/(Float_t)AODtrack1->Pt())) ){
		  continue;
		}
	  }
	  
	  // n cluster Min cut, do not use together with crossed row
	  if (fTPCclustMin > 0) { 
		if (trk1TpcNC < fTPCclustMin) {
		  continue;
	    }
	  }
	  
	  // shared cluster 
	  if( fTPCsharedCut > 0 && AODtrack1->GetTPCnclsS() > fTPCsharedCut){
		continue;
      }
      
	  //Apply track cuts for TPC EP here:
	  //if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && !((trk1Eta >= fEtaGapNeg) && (trk1Eta <= fEtaGapPos)) && (trk1dEdx >= fTPCdEdxMin) && (trk1TpcNC >= fTPCclustMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fTrkChi2Max) && TMath::Abs(trk1Chrg)) {
	  if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && (trk1dEdx >= fTPCdEdxMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fTrkChi2Max) && TMath::Abs(trk1Chrg)) {
		// ================================ save Qvec ===================================
		
		// Calculate Re[Q_{m,k}] and Im[Q_{m,k}], (m = 1,2,3,4,5,6 and k = 0,1,2,3) for this event:
		for(Int_t h=0;h<2;h++) 
		{
		  // RP All, OS
		  // S_{p,k} is in it
		  fCMEQReRP->Fill(h+0.5, TMath::Cos((h+1.)*trk1Phi)); 
		  fCMEQImRP->Fill(h+0.5, TMath::Sin((h+1.)*trk1Phi)); 
		
		  if (trk1Chrg>0) {
			fCMEQRePOIPos->Fill(h+0.5, TMath::Cos((h+1.)*trk1Phi)); 
			fCMEQImPOIPos->Fill(h+0.5, TMath::Sin((h+1.)*trk1Phi)); 
		  } 
		
		  if (trk1Chrg<0) {
			fCMEQRePOINeg->Fill(h+0.5, TMath::Cos((h+1.)*trk1Phi)); 
			fCMEQImPOINeg->Fill(h+0.5, TMath::Sin((h+1.)*trk1Phi)); 
		  }
		
		}
	
	
		// Calculate <<cos(2a-2Psi_V0)>> for RP, POI OS
		f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C)));
		f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent->Fill(1.5, TMath::Cos(2*(trk1Phi-fPsi2V0A)));
		f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCC)));
		f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->Fill(1.5, TMath::Cos(2*(trk1Phi-fPsiZDCA)));
		f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->Fill(2.5, TMath::Cos(2*(trk1Phi-fPsiZDCCA)));
		
		f2pCorrelatorCos2PsiDiff2PsiV0CRP->Fill(centrality, TMath::Cos(2*(trk1Phi-fPsi2V0C))); //<cos(2psi1-2phi_V0C)> // no event-by-event
		f2pCorrelatorCos2PsiDiff2PsiV0ARP->Fill(centrality, TMath::Cos(2*(trk1Phi-fPsi2V0A))); //<cos(2psi1-2phi_V0A)> // no event-by-event
		// Calculate <<cos(2a-2Psi_ZDC)>> for RP, POI OS
		f2pCorrelatorCos2PsiDiff2PsiZDCCRP->Fill(centrality, TMath::Cos(2*(trk1Phi-fPsiZDCC))); //<cos(2psi1-2phi_ZDCC)> // no event-by-event
		f2pCorrelatorCos2PsiDiff2PsiZDCARP->Fill(centrality, TMath::Cos(2*(trk1Phi-fPsiZDCA))); //<cos(2psi1-2phi_ZDCA)> // no event-by-event
		f2pCorrelatorCos2PsiDiff2PsiZDCCARP->Fill(centrality, TMath::Cos(2*(trk1Phi-fPsiZDCCA))); //<cos(2psi1-2phi_ZDCCA)> // no event-by-event
	
		fNITV0OS->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		fNITV0OS->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		fNITV0OS->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		fNITV0OS->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		fNITZDCOS->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZDCC)); // <<cos(psi1-2phi_ZDCC)>>
		fNITZDCOS->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZDCC)); // <<sin(psi1-2phi_ZDCC)>>
		fNITZDCOS->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZDCA)); // <<cos(psi1-2phi_ZDCA)>>
		fNITZDCOS->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZDCA)); // <<sin(psi1-2phi_ZDCA)>>
		fNITZDCOS->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZDCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		fNITZDCOS->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZDCCA)); // <<sin(psi1-2phi_ZDCCA)>>
	
		if (trk1Chrg>0) {
		  // Calculate <<cos(2a-2Psi_V0)>> for POI Pos
		  //f2pCorrelatorCos2PsiDiff2PsiV0POIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C))); 
		  //f2pCorrelatorCos2PsiDiff2PsiV0POIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0A)));
		  // Calculate <<cos(2a-2Psi_ZDC)>> for POI Pos
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCC)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCAPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCA)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCAPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCCA)));
		
		  fNITV0POIPos->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		  fNITV0POIPos->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		  fNITV0POIPos->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		  fNITV0POIPos->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		  fNITZDCPOIPos->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZDCC)); // <<cos(psi1-2phi_ZDCC)>>
		  fNITZDCPOIPos->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZDCC)); // <<sin(psi1-2phi_ZDCC)>>
		  fNITZDCPOIPos->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZDCA)); // <<cos(psi1-2phi_ZDCA)>>
		  fNITZDCPOIPos->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZDCA)); // <<sin(psi1-2phi_ZDCA)>>
		  fNITZDCPOIPos->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZDCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		  fNITZDCPOIPos->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZDCCA)); // <<sin(psi1-2phi_ZDCCA)>>
		}
	
		if (trk1Chrg<0) {
		  // Calculate <<cos(2a-2Psi_V0)>> for POI Neg
		  //f2pCorrelatorCos2PsiDiff2PsiV0POINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C))); 
		  //f2pCorrelatorCos2PsiDiff2PsiV0POINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0A)));
		  // Calculate <<cos(2a-2Psi_ZDC)>> for POI Neg
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCC)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCAPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCA)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCAPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZDCCA)));
		
		  fNITV0POINeg->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		  fNITV0POINeg->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		  fNITV0POINeg->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		  fNITV0POINeg->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		  fNITZDCPOINeg->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZDCC)); // <<cos(psi1-2phi_ZDCC)>>
		  fNITZDCPOINeg->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZDCC)); // <<sin(psi1-2phi_ZDCC)>>
		  fNITZDCPOINeg->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZDCA)); // <<cos(psi1-2phi_ZDCA)>>
		  fNITZDCPOINeg->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZDCA)); // <<sin(psi1-2phi_ZDCA)>>
		  fNITZDCPOINeg->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZDCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		  fNITZDCPOINeg->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZDCCA)); // <<sin(psi1-2phi_ZDCCA)>>
		}

		///---> 2nd track Loop   
		for(Int_t jTrack = iTrack+1; jTrack < ntracks; jTrack++) {  // garanteed no overlap

		  /// skip autocorrelation:
		  //if(jTrack==iTrack) continue; ///////////////////// Delete later

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

			
			// for constrained TPConly tracks
		    if(fFilterBit == 128){
			  trk2DCAxy = AODtrack2->DCA();      // this is the DCA from global track (not exactly what is cut on)
			  trk2DCAz  = AODtrack2->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
		    }
		    else{
			  //const AliVVertex *vertex = fVevent->GetPrimaryVertex();
			  //vertex->GetXYZ(vTrk2);
			  AODtrack2->GetXYZ(posTrk2);
			  trk2DCAxy  = TMath::Sqrt((posTrk2[0] - pVtxX)*(posTrk2[0] - pVtxX) + (posTrk2[1] - pVtxY)*(posTrk2[1] - pVtxY));
			  trk2DCAz   = posTrk2[2] - pVtxZ;
		    }
		    
		    if(fDCAxyMax>0 && fDCAzzMax>0){
		      if(TMath::Sqrt((trk2DCAxy*trk2DCAxy)/(fDCAxyMax*fDCAxyMax)+(trk2DCAz*trk2DCAz)/(fDCAzzMax*fDCAzzMax)) > 1 ){
		        continue;  // 2D cut
		      }
      		}
      		
		    //if (fDCAxyMax>0 && fDCAzzMax>0 && trk2DCAxy>=fDCAxyMax && trk2DCAz>=fDCAzzMax) { // if fDCAxyMax, fDCAzzMax is set to be less than 0, no cut applied
			  //continue;
			//}
			
			if (bUseTPCCrossedRows){
			  if ((Float_t)AODtrack2->GetTPCNCrossedRows() < (120 - (5/(Float_t)AODtrack2->Pt())) ){
				continue;
			  }
			}
			
			if (fTPCclustMin > 0) { 
			  if (trk2TpcNC < fTPCclustMin) {
			    continue;
			  }
		    }
		    
		    // shared cluster 
			if( fTPCsharedCut > 0 && AODtrack2->GetTPCnclsS() > fTPCsharedCut){
			  continue;
			}
  
			//Apply track cuts for second track
			if((trk2Pt <= fMaxPtCut) && (trk2Pt >= fMinPtCut) && (trk2Eta <= fMaxEtaCut) && (trk2Eta >= fMinEtaCut) && (trk2dEdx >= fTPCdEdxMin) && (trk2Chi2 >= fTrkChi2Min) && (trk2Chi2 <= fTrkChi2Max) && TMath::Abs(trk2Chrg)) {
				// v2 for TPC
				f2pCorrelatorCosPsiDiff->Fill(centrality, TMath::Cos(trk1Phi-trk2Phi)); // no event-by-event
				f2pCorrelatorCos2PsiDiff->Fill(centrality, TMath::Cos(2*(trk1Phi-trk2Phi))); // no event-by-event
				
				f2pCorrelatorCosPsiDiffPerEvent->Fill(0.5, TMath::Cos(trk1Phi-trk2Phi));
				f2pCorrelatorCos2PsiDiffPerEvent->Fill(0.5, TMath::Cos(2*(trk1Phi-trk2Phi)));
				
				if(trk1Chrg*trk2Chrg < 0){ //Opposite sign	
					fRePEBEOS->Fill(0.5, TMath::Cos(trk1Phi+trk2Phi));
					fImPEBEOS->Fill(0.5, TMath::Sin(trk1Phi+trk2Phi));
					
					f2pCorrelatorCosPsiDiffOS->Fill(0.5, TMath::Cos(trk1Phi-trk2Phi)); 
					f2pCorrelatorCos2PsiDiffOS->Fill(0.5, TMath::Cos(2*(trk1Phi-trk2Phi)));
				} else if(trk1Chrg > 0 && trk2Chrg > 0){		      
					fRePEBEPP->Fill(0.5, TMath::Cos(trk1Phi+trk2Phi));
					fImPEBEPP->Fill(0.5, TMath::Sin(trk1Phi+trk2Phi));
					
					f2pCorrelatorCosPsiDiffPP->Fill(0.5, TMath::Cos(trk1Phi-trk2Phi));
					f2pCorrelatorCos2PsiDiffPP->Fill(0.5, TMath::Cos(2*(trk1Phi-trk2Phi)));
					
				} else{ //else if(trk1Chrg < 0 && trk2Chrg < 0){  ///this is obvious!
					fRePEBENN->Fill(0.5, TMath::Cos(trk1Phi+trk2Phi));
					fImPEBENN->Fill(0.5, TMath::Sin(trk1Phi+trk2Phi));
					
					f2pCorrelatorCosPsiDiffNN->Fill(0.5, TMath::Cos(trk1Phi-trk2Phi));
					f2pCorrelatorCos2PsiDiffNN->Fill(0.5, TMath::Cos(2*(trk1Phi-trk2Phi)));
				}
			}//j-track cuts
		  }//j-track FB validated
		}///j-track loop
      }//----> i-track loop => All trackCuts applied.     
    }//-----> i-track loop => FB is validated.    
  }///-----> i-track loop Ends here.<--------
  

  Calculate2pCorrelator();
  CalculateNonIsotropicTerms();
  CalculateDifferential3pCorrelator();
  
  // Reset event-by-event variables
  ResetEventByEventQuantities();
  
  fDebugwEventCount->Fill(8.1); ///Left for Analysis
  
  fHistVertexZcm->Fill(pVtxZ);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

}//---------------- UserExec ----------------------

void AliAnalysisTaskGammaNonIsotropicCorr::ResetEventByEventQuantities(){
  centrality = -99.0;
  fPsi2V0C = 0.;
  fPsi2V0A = 0.;
  fPsiZDCC = 0.;
  fPsiZDCA = 0.;
  fPsiZDCCA = 0; 
  fCMEQReRP->Reset();
  fCMEQImRP->Reset();
  f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent->Reset();
  f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->Reset();
  f2pCorrelatorCosPsiDiffPerEvent->Reset();
  f2pCorrelatorCos2PsiDiffPerEvent->Reset();
  fCMEQRePOIPos->Reset();
  fCMEQImPOIPos->Reset();
  fCMEQRePOINeg->Reset();
  fCMEQImPOINeg->Reset();
  fNITV0OS->Reset();
  fNITZDCOS->Reset();
  fNITV0POIPos->Reset();
  fNITZDCPOIPos->Reset();
  fNITV0POINeg->Reset();
  fNITZDCPOINeg->Reset();

  fRePEBEOS->Reset();
  fImPEBEOS->Reset();
  f2pCorrelatorCosPsiDiffOS->Reset();
  f2pCorrelatorCos2PsiDiffOS->Reset();
  fRePEBEPP->Reset();
  fImPEBEPP->Reset();
  f2pCorrelatorCosPsiDiffPP->Reset();
  f2pCorrelatorCos2PsiDiffPP->Reset();
  fRePEBENN->Reset();
  fImPEBENN->Reset();
  f2pCorrelatorCosPsiDiffNN->Reset();
  f2pCorrelatorCos2PsiDiffNN->Reset();
  
}

void AliAnalysisTaskGammaNonIsotropicCorr::SetupAnalysisHistograms(){
  
  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables
  Char_t  name[100];
  Char_t title[100];

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
  //hV0CTPCPsiNCorrelation = new TProfile("hV0CTPCPsiNCorrelation",Form("V0C-TPC Psi%d; Cent; Resolution",gHarmonic),90,0,90);
  //fListHist->Add(hV0CTPCPsiNCorrelation);
  //hV0ATPCPsiNCorrelation = new TProfile("hV0ATPCPsiNCorrelation",Form("V0A-TPC Psi%d; Cent; Resolution",gHarmonic),90,0,90);
  //fListHist->Add(hV0ATPCPsiNCorrelation);

  hV0CV0APsi3Correlation = new TProfile("hV0CV0APsi3Correlation",Form("V0C-V0A Psi%d; Cent; Resolution",3),90,0,90);
  fListHist->Add(hV0CV0APsi3Correlation);
  //hV0CTPCPsi3Correlation = new TProfile("hV0CTPCPsi3Correlation",Form("V0C-TPC Psi%d; Cent; Resolution",3),90,0,90);
  //fListHist->Add(hV0CTPCPsi3Correlation);
  //hV0ATPCPsi3Correlation = new TProfile("hV0ATPCPsi3Correlation",Form("V0A-TPC Psi%d; Cent; Resolution",3),90,0,90);
  //fListHist->Add(hV0ATPCPsi3Correlation);
  
  /// Q vec related histograms
  Double_t fCRCEtaBinEdges[fCMEnEtaBin+1] = {-0.8, -0.1, 0, 0.1, 0.8};
  
  /// 2-p correlator 
  f2pCorrelatorCos2PsiDiff2PsiV0CRP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0CRP", "f2pCorrelatorCos2PsiDiff2PsiV0CRP; Centrality; #LTcos(2#phi_{RP}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiV0CRP->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiV0CRP);
  
  f2pCorrelatorCos2PsiDiff2PsiV0ARP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0ARP", "f2pCorrelatorCos2PsiDiff2PsiV0ARP; Centrality; #LTcos(2#phi_{RP}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiV0ARP->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiV0ARP);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCCRP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCCRP", "f2pCorrelatorCos2PsiDiff2PsiZDCCRP; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCCRP->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCCRP);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCARP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCARP", "f2pCorrelatorCos2PsiDiff2PsiZDCARP; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCARP->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCARP);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCCARP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCCARP", "f2pCorrelatorCos2PsiDiff2PsiZDCCARP; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCCARP->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCCARP);
  
  f2pCorrelatorCosPsiDiff = new TProfile("f2pCorrelatorCosPsiDiff", "f2pCorrelatorCosPsiDiff; Centrality; #LTcos(#phi_{1}-#phi_{2})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCosPsiDiff->Sumw2();
  fListHist->Add(f2pCorrelatorCosPsiDiff);
  
  f2pCorrelatorCos2PsiDiff = new TProfile("f2pCorrelatorCos2PsiDiff", "f2pCorrelatorCos2PsiDiff; Centrality; #LTcos(2#phi_{1}-2#phi_{2})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff);
  
  f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE", "f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE; Centrality; #LTcos(2#phi_{RP}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE);
  
  f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE", "f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE; Centrality; #LTcos(2#phi_{RP}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE", "f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE", "f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE", "f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE; Centrality; #LTcos(2#phi_{RP}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE);
  
  f2pCorrelatorCosPsiDiffEBE = new TProfile("f2pCorrelatorCosPsiDiffEBE", "f2pCorrelatorCosPsiDiffEBE; Centrality; #LTcos(#phi_{1}-#phi_{2})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCosPsiDiffEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCosPsiDiffEBE);
  
  f2pCorrelatorCos2PsiDiffEBE = new TProfile("f2pCorrelatorCos2PsiDiffEBE", "f2pCorrelatorCos2PsiDiffEBE; Centrality; #LTcos(2#phi_{1}-2#phi_{2})#GT",fCMESPPPCenBin,0,90);
  f2pCorrelatorCos2PsiDiffEBE->Sumw2();
  fListHist->Add(f2pCorrelatorCos2PsiDiffEBE);

  // 3-p correlator
  f3pCorrelatorTPCOSPro = new TProfile("f3pCorrelatorTPCOSPro","f3pCorrelatorTPCOSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT",fCMESPPPCenBin,0.,90.); // cos[n(psi1+psi2-2phi3)] OS
  f3pCorrelatorTPCOSPro->Sumw2();
  fListHist->Add(f3pCorrelatorTPCOSPro);
  
  f3pCorrelatorV0COSPro = new TProfile("f3pCorrelatorV0COSPro","f3pCorrelatorV0COSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0COSPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0COSPro);
  
  f3pCorrelatorV0AOSPro = new TProfile("f3pCorrelatorV0AOSPro","f3pCorrelatorV0AOSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0AOSPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0AOSPro);
  
  f3pCorrelatorZDCCOSPro = new TProfile("f3pCorrelatorZDCCOSPro","f3pCorrelatorZDCCOSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCOSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCOSPro);
  
  f3pCorrelatorZDCAOSPro = new TProfile("f3pCorrelatorZDCAOSPro","f3pCorrelatorZDCAOSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCAOSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCAOSPro);
  
  f3pCorrelatorZDCCAOSPro = new TProfile("f3pCorrelatorZDCCAOSPro","f3pCorrelatorZDCCAOSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCAOSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCAOSPro);
    
    
  f3pCorrelatorTPCPPPro = new TProfile("f3pCorrelatorTPCPPPro","f3pCorrelatorTPCPPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT",fCMESPPPCenBin,0.,90.); // cos[n(psi1+psi2-2phi3)] PP
  f3pCorrelatorTPCPPPro->Sumw2();
  fListHist->Add(f3pCorrelatorTPCPPPro);
  
  f3pCorrelatorV0CPPPro = new TProfile("f3pCorrelatorV0CPPPro","f3pCorrelatorV0CPPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0CPPPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0CPPPro);
  
  f3pCorrelatorV0APPPro = new TProfile("f3pCorrelatorV0APPPro","f3pCorrelatorV0APPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0APPPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0APPPro);
  
  f3pCorrelatorZDCCPPPro = new TProfile("f3pCorrelatorZDCCPPPro","f3pCorrelatorZDCCPPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCPPPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCPPPro);
  
  f3pCorrelatorZDCAPPPro = new TProfile("f3pCorrelatorZDCAPPPro","f3pCorrelatorZDCAPPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCAPPPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCAPPPro);
  
  f3pCorrelatorZDCCAPPPro = new TProfile("f3pCorrelatorZDCCAPPPro","f3pCorrelatorZDCCAPPPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCAPPPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCAPPPro);
  
  
  f3pCorrelatorTPCNNPro = new TProfile("f3pCorrelatorTPCNNPro","f3pCorrelatorTPCNNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT",fCMESPPPCenBin,0.,90.); // cos[n(psi1+psi2-2phi3)] NN
  f3pCorrelatorTPCNNPro->Sumw2();
  fListHist->Add(f3pCorrelatorTPCNNPro);
  
  f3pCorrelatorV0CNNPro = new TProfile("f3pCorrelatorV0CNNPro","f3pCorrelatorV0CNNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0CNNPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0CNNPro);
  
  f3pCorrelatorV0ANNPro = new TProfile("f3pCorrelatorV0ANNPro","f3pCorrelatorV0ANNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0ANNPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0ANNPro);
  
  f3pCorrelatorZDCCNNPro = new TProfile("f3pCorrelatorZDCCNNPro","f3pCorrelatorZDCCNNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCNNPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCNNPro);
  
  f3pCorrelatorZDCANNPro = new TProfile("f3pCorrelatorZDCANNPro","f3pCorrelatorZDCANNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCANNPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCANNPro);
  
  f3pCorrelatorZDCCANNPro = new TProfile("f3pCorrelatorZDCCANNPro","f3pCorrelatorZDCCANNPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCANNPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCANNPro);
  
  
  f3pCorrelatorTPCSSPro = new TProfile("f3pCorrelatorTPCSSPro","f3pCorrelatorTPCSSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#phi_{3})#GT",fCMESPPPCenBin,0.,90.); // cos[n(psi1+psi2-2phi3)] SS
  f3pCorrelatorTPCSSPro->Sumw2();
  fListHist->Add(f3pCorrelatorTPCSSPro);
  
  f3pCorrelatorV0CSSPro = new TProfile("f3pCorrelatorV0CSSPro","f3pCorrelatorV0CSSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0C})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0CSSPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0CSSPro);
  
  f3pCorrelatorV0ASSPro = new TProfile("f3pCorrelatorV0ASSPro","f3pCorrelatorV0ASSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{V0A})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorV0ASSPro->Sumw2();
  fListHist->Add(f3pCorrelatorV0ASSPro);
  
  f3pCorrelatorZDCCSSPro = new TProfile("f3pCorrelatorZDCCSSPro","f3pCorrelatorZDCCSSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCC})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCSSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCSSPro);
  
  f3pCorrelatorZDCASSPro = new TProfile("f3pCorrelatorZDCASSPro","f3pCorrelatorZDCASSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCASSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCASSPro);
  
  f3pCorrelatorZDCCASSPro = new TProfile("f3pCorrelatorZDCCASSPro","f3pCorrelatorZDCCASSPro; Centrality; #LTcos(#phi_{1}+#phi_{2}-2#Psi_{ZDCCA})#GT",fCMESPPPCenBin,0.,90.);
  f3pCorrelatorZDCCASSPro->Sumw2();
  fListHist->Add(f3pCorrelatorZDCCASSPro);
	
  /// non-isotropic terms
  TString NonIsotropicTermsTPCRPTitle[10] = {"#LT#LTcos(#phi_{RP})#GT#GT", "#LT#LTsin(#phi_{RP})#GT#GT", "#LT#LTcos(2#phi_{RP})#GT#GT", "#LT#LTsin(2#phi_{RP})#GT#GT", "#LT#LTcos(#phi_{1}+#phi_{2})#GT#GT", "#LT#LTsin(#phi_{1}+#phi_{2})#GT#GT", "#LT#LTcos(2#phi_{1}-#phi_{2})#GT#GT", "#LT#LTsin(2#phi_{1}-#phi_{2})#GT#GT", "#LT#LTcos(#phi_{1}-#phi_{2}-#phi_{3})#GT#GT", "#LT#LTsin(#phi_{1}-#phi_{2}-#phi_{3})#GT#GT"};
    
  TString NonIsotropicTermsTPCTitle[10] = {"#LT#LTcos(#phi_{POI1})#GT#GT", "#LT#LTsin(#phi_{POI1})#GT#GT", "#LT#LTcos(#phi_{POI2})#GT#GT", "#LT#LTsin(#phi_{POI2})#GT#GT", "#LT#LTcos(#phi_{POI1}-2#phi_{RP})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#phi_{RP})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#phi_{RP})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#phi_{RP})#GT#GT", "#LT#LTcos(#phi_{POI1}+#phi_{POI2})#GT#GT", "#LT#LTsin(#phi_{POI1}+#phi_{POI2})#GT#GT"};
    
  // non-isotropic terms 3-p 
  for(Int_t i = 0; i < 10; i++) {
    fNonIsotropicTermsPro[i] = new TProfile(Form("fNonIsotropicTermsPro[%d]",i),Form("fNonIsotropicTermsPro[%d]; Centrality; %s",i,NonIsotropicTermsTPCRPTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsPro[i]);
    
    fNonIsotropicTermsOSPro[i] = new TProfile(Form("fNonIsotropicTermsOSPro[%d]",i),Form("fNonIsotropicTermsOSPro[%d]; Centrality; %s",i,NonIsotropicTermsTPCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsOSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsOSPro[i]);
    
    fNonIsotropicTermsPPPro[i] = new TProfile(Form("fNonIsotropicTermsPPPro[%d]",i),Form("fNonIsotropicTermsPPPro[%d]; Centrality; %s",i,NonIsotropicTermsTPCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsPPPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsPPPro[i]);
    
    fNonIsotropicTermsNNPro[i] = new TProfile(Form("fNonIsotropicTermsNNPro[%d]",i),Form("fNonIsotropicTermsNNPro[%d]; Centrality; %s",i,NonIsotropicTermsTPCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsNNPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsNNPro[i]);
    
    fNonIsotropicTermsSSPro[i] = new TProfile(Form("fNonIsotropicTermsSSPro[%d]",i),Form("fNonIsotropicTermsSSPro[%d]; Centrality; %s",i,NonIsotropicTermsTPCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsSSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsSSPro[i]);
  }
  
  // Cos(fPsi2V0C), Sin(fPsi2V0C), Cos(fPsi2V0A), Sin(fPsi2V0A)
                                        // Cos(2*fPsi2V0C), Sin(2*fPsi2V0C), Cos(2*fPsi2V0A), Sin(2*fPsi2V0A)
  TString NonIsotropicTermsV0RPTitle[8] = {"#LT#LTcos(#Psi_{V0C})#GT#GT", "#LT#LTsin(#Psi_{V0C})#GT#GT", "#LT#LTcos(#Psi_{V0A})#GT#GT", "#LT#LTsin(#Psi_{V0A})#GT#GT", "#LT#LTcos(2#Psi_{V0C})#GT#GT", "#LT#LTsin(2#Psi_{V0C})#GT#GT", "#LT#LTcos(2#Psi_{V0A})#GT#GT", "#LT#LTsin(2#Psi_{V0A})#GT#GT"};
  
  TString NonIsotropicTermsV0Title[8] = {"#LT#LTcos(#phi_{POI1}-2#Psi_{V0C})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#Psi_{V0C})#GT#GT", "#LT#LTcos(#phi_{POI1}-2#Psi_{V0A})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#Psi_{V0A})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#Psi_{V0C})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#Psi_{V0C})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#Psi_{V0A})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#Psi_{V0A})#GT#GT"};
  
  
  for(Int_t i = 0; i < 8; i++) {
    fNonIsotropicTermsV0Pro[i] = new TProfile(Form("fNonIsotropicTermsV0Pro[%d]",i),Form("fNonIsotropicTermsV0Pro[%d]; Centrality; %s",i,NonIsotropicTermsV0RPTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsV0Pro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsV0Pro[i]);
    
    fNonIsotropicTermsV0OSPro[i] = new TProfile(Form("fNonIsotropicTermsV0OSPro[%d]",i),Form("fNonIsotropicTermsV0OSPro[%d]; Centrality; %s",i,NonIsotropicTermsV0Title[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsV0OSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsV0OSPro[i]);
    
    fNonIsotropicTermsV0PPPro[i] = new TProfile(Form("fNonIsotropicTermsV0PPPro[%d]",i),Form("fNonIsotropicTermsV0PPPro[%d]; Centrality; %s",i,NonIsotropicTermsV0Title[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsV0PPPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsV0PPPro[i]);
    
    fNonIsotropicTermsV0NNPro[i] = new TProfile(Form("fNonIsotropicTermsV0NNPro[%d]",i),Form("fNonIsotropicTermsV0NNPro[%d]; Centrality; %s",i,NonIsotropicTermsV0Title[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsV0NNPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsV0NNPro[i]);
    
    fNonIsotropicTermsV0SSPro[i] = new TProfile(Form("fNonIsotropicTermsV0SSPro[%d]",i),Form("fNonIsotropicTermsV0SSPro[%d]; Centrality; %s",i,NonIsotropicTermsV0Title[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsV0SSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsV0SSPro[i]);
  }
  
  TString NonIsotropicTermsZDCRPTitle[12] = {"#LT#LTcos(#Psi_{ZDCC})#GT#GT", "#LT#LTsin(#Psi_{ZDCC})#GT#GT", "#LT#LTcos(#Psi_{ZDCA})#GT#GT", "#LT#LTsin(#Psi_{ZDCA})#GT#GT", "#LT#LTcos(#Psi_{ZDCCA})#GT#GT", "#LT#LTsin(#Psi_{ZDCCA})#GT#GT", "#LT#LTcos(2#Psi_{ZDCC})#GT#GT", "#LT#LTsin(2#Psi_{ZDCC})#GT#GT", "#LT#LTcos(2#Psi_{ZDCA})#GT#GT", "#LT#LTsin(2#Psi_{ZDCA})#GT#GT", "#LT#LTcos(2#Psi_{ZDCCA})#GT#GT", "#LT#LTsin(2#Psi_{ZDCCA})#GT#GT"};
  
  TString NonIsotropicTermsZDCTitle[12] = {"#LT#LTcos(#phi_{POI1}-2#Psi_{ZDCC})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#Psi_{ZDCC})#GT#GT", "#LT#LTcos(#phi_{POI1}-2#Psi_{ZDCA})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#Psi_{ZDCA})#GT#GT", "#LT#LTcos(#phi_{POI1}-2#Psi_{ZDCCA})#GT#GT", "#LT#LTsin(#phi_{POI1}-2#Psi_{ZDCCA})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#Psi_{ZDCC})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#Psi_{ZDCC})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#Psi_{ZDCA})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#Psi_{ZDCA})#GT#GT", "#LT#LTcos(#phi_{POI2}-2#Psi_{ZDCCA})#GT#GT", "#LT#LTsin(#phi_{POI2}-2#Psi_{ZDCCA})#GT#GT"};
  
  for(Int_t i = 0; i < 12; i++) {
    fNonIsotropicTermsZDCPro[i] = new TProfile(Form("fNonIsotropicTermsZDCPro[%d]",i),Form("fNonIsotropicTermsZDCPro[%d]; Centrality; %s",i,NonIsotropicTermsZDCRPTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsZDCPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsZDCPro[i]);
    
    fNonIsotropicTermsZDCOSPro[i] = new TProfile(Form("fNonIsotropicTermsZDCOSPro[%d]",i),Form("fNonIsotropicTermsZDCOSPro[%d]; Centrality; %s",i,NonIsotropicTermsZDCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsZDCOSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsZDCOSPro[i]);
    
    fNonIsotropicTermsZDCPPPro[i] = new TProfile(Form("fNonIsotropicTermsZDCPPPro[%d]",i),Form("fNonIsotropicTermsZDCPPPro[%d]; Centrality; %s",i,NonIsotropicTermsZDCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsZDCPPPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsZDCPPPro[i]);
    
    fNonIsotropicTermsZDCNNPro[i] = new TProfile(Form("fNonIsotropicTermsZDCNNPro[%d]",i),Form("fNonIsotropicTermsZDCNNPro[%d]; Centrality; %s",i,NonIsotropicTermsZDCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsZDCNNPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsZDCNNPro[i]);
    
    fNonIsotropicTermsZDCSSPro[i] = new TProfile(Form("fNonIsotropicTermsZDCSSPro[%d]",i),Form("fNonIsotropicTermsZDCSSPro[%d]; Centrality; %s",i,NonIsotropicTermsZDCTitle[i].Data()),fCMESPPPCenBin,0.,90.); 
    fNonIsotropicTermsZDCSSPro[i]->Sumw2();
    fListHist->Add(fNonIsotropicTermsZDCSSPro[i]);
  }
  
}

//=======================================================================================================================

void AliAnalysisTaskGammaNonIsotropicCorr::Calculate2pCorrelator()
{
	f2pCorrelatorCosPsiDiffEBE->Fill(centrality, f2pCorrelatorCosPsiDiffPerEvent->GetBinContent(1));
	f2pCorrelatorCos2PsiDiffEBE->Fill(centrality, f2pCorrelatorCos2PsiDiffPerEvent->GetBinContent(1));
	f2pCorrelatorCos2PsiDiff2PsiV0CRPEBE->Fill(centrality, f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent->GetBinContent(1)); //<cos(2psi1-2phi_V0C)> // no event-by-event
	f2pCorrelatorCos2PsiDiff2PsiV0ARPEBE->Fill(centrality, f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent->GetBinContent(2)); //<cos(2psi1-2phi_V0A)> // no event-by-event
	// Calculate <<cos(2a-2Psi_ZDC)>> for RP, POI OS
	f2pCorrelatorCos2PsiDiff2PsiZDCCRPEBE->Fill(centrality, f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->GetBinContent(1)); //<cos(2psi1-2phi_ZDCC)> // no event-by-event
	f2pCorrelatorCos2PsiDiff2PsiZDCARPEBE->Fill(centrality, f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->GetBinContent(2)); //<cos(2psi1-2phi_ZDCA)> // no event-by-event
	f2pCorrelatorCos2PsiDiff2PsiZDCCARPEBE->Fill(centrality, f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent->GetBinContent(3)); //<cos(2psi1-2phi_ZDCCA)> // no event-by-event		
}

//=======================================================================================================================

void AliAnalysisTaskGammaNonIsotropicCorr::CalculateNonIsotropicTerms()
{
 // Calculate non-isotropic terms which appear in the decomposition of 3-p correlator <cos[n(phi1+phi2-2phi3)]>.
 
 // a) Calculate without using particle weights;
 
 // For detector with uniform acceptance all these terms vanish. These non-isotropic terms are stored in fNonIsotropicTermsPro.
 // Binning of fNonIsotropicTermsPro is organized as follows:
 //  1st bin: <<cos(n*phi1)>>
 //  2nd bin: <<sin(n*phi1)>>
 //  3rd bin: <<cos(2n*phi1)>>
 //  4th bin: <<sin(2n*phi1)>>
 //  5th bin: <<cos(n*(phi1+phi2)>>
 //  6th bin: <<sin(n*(phi1+phi2)>>
 //  7th bin: <<cos(n*(2phi1-phi2)>>
 //  8th bin: <<sin(n*(2phi1-phi2)>>
 //  9th bin: <<cos(n*(phi1-phi2-phi3)>> // not needed
 // 10th bin: <<sin(n*(phi1-phi2-phi3)>> // not needed 
	
 // a) Calculate without using particle weights:
  // Multiplicity (number of RPs):
  Double_t dMult = fCMEQReRP->GetBinEntries(1); // (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonics n and 2n: 
  Double_t dReQ1n = fCMEQReRP->GetBinContent(1)*dMult; // (*fReQnk)(0,0) Qvec TPC cos(phi) 
  Double_t dReQ2n = fCMEQReRP->GetBinContent(2)*dMult; // (*fReQnk)(1,0); 
  Double_t dImQ1n = fCMEQImRP->GetBinContent(1)*dMult; // (*fImQnk)(0,0);
  Double_t dImQ2n = fCMEQImRP->GetBinContent(2)*dMult; // (*fImQnk)(1,0);
  // 1-particle terms:
  Double_t cosP1n = 0.; // <cos(n*(phi1))>  
  Double_t sinP1n = 0.; // <sin(n*(phi1))>
  Double_t cosP2n = 0.; // <cos(2n*(phi1))>  
  Double_t sinP2n = 0.; // <sin(2n*(phi1))>
  if(dMult>0) 
  { 
   cosP1n = dReQ1n/dMult; 
   sinP1n = dImQ1n/dMult;
   cosP2n = dReQ2n/dMult; 
   sinP2n = dImQ2n/dMult;   
   // All-event avarages:
   fNonIsotropicTermsPro[0]->Fill(centrality,cosP1n,dMult); // <<cos(n*(phi1))>> 
   fNonIsotropicTermsPro[1]->Fill(centrality,sinP1n,dMult); // <<sin(n*(phi1))>>   
   fNonIsotropicTermsPro[2]->Fill(centrality,cosP2n,dMult); // <<cos(2n*(phi1))>> 
   fNonIsotropicTermsPro[3]->Fill(centrality,sinP2n,dMult); // <<sin(2n*(phi1))>>    
  } // end of if(dMult>0) 
  // 2-particle terms:
  Double_t cosP1nP1n = 0.; // <cos(n*(phi1+phi2))>
  Double_t sinP1nP1n = 0.; // <sin(n*(phi1+phi2))>
  Double_t cosP2nM1n = 0.; // <cos(n*(2phi1-phi2))>
  Double_t sinP2nM1n = 0.; // <sin(n*(2phi1-phi2))>
  if(dMult>1)
  {
   cosP1nP1n = (pow(dReQ1n,2)-pow(dImQ1n,2)-dReQ2n)/(dMult*(dMult-1)); 
   sinP1nP1n = (2.*dReQ1n*dImQ1n-dImQ2n)/(dMult*(dMult-1)); 
   cosP2nM1n = (dReQ2n*dReQ1n+dImQ2n*dImQ1n-dReQ1n)/(dMult*(dMult-1)); 
   sinP2nM1n = (dImQ2n*dReQ1n-dReQ2n*dImQ1n-dImQ1n)/(dMult*(dMult-1)); 
   // All-event avarages:
   fNonIsotropicTermsPro[4]->Fill(centrality,cosP1nP1n,dMult*(dMult-1.)); // <<cos(n*(phi1+phi2))>> 
   fNonIsotropicTermsPro[5]->Fill(centrality,sinP1nP1n,dMult*(dMult-1.)); // <<sin(n*(phi1+phi2))>>   
   fNonIsotropicTermsPro[6]->Fill(centrality,cosP2nM1n,dMult*(dMult-1.)); // <<cos(n*(2phi1-phi2))>> 
   fNonIsotropicTermsPro[7]->Fill(centrality,sinP2nM1n,dMult*(dMult-1.)); // <<sin(n*(2phi1-phi2))>>   
  } // end of if(dMult>1) 
  // 3-particle: correct and ready but not needed, hence commented out.
  Double_t cosP1nM1nM1n = 0.; // <cos(n*(phi1-phi2-phi3))>
  Double_t sinP1nM1nM1n = 0.; // <sin(n*(phi1-phi2-phi3))>
  if(dMult>2)
  {
   cosP1nM1nM1n = (dReQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))-dReQ1n*dReQ2n-dImQ1n*dImQ2n-2.*(dMult-1)*dReQ1n)
                / (dMult*(dMult-1)*(dMult-2)); 
   sinP1nM1nM1n = (-dImQ1n*(pow(dReQ1n,2)+pow(dImQ1n,2))+dReQ1n*dImQ2n-dImQ1n*dReQ2n+2.*(dMult-1)*dImQ1n)
                / (dMult*(dMult-1)*(dMult-2));              
   // All-events avarages:
   fNonIsotropicTermsPro[8]->Fill(centrality,cosP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<cos(n*(phi1-phi2-phi3))>> 
   fNonIsotropicTermsPro[9]->Fill(centrality,sinP1nM1nM1n,dMult*(dMult-1.)*(dMult-2.)); // <<sin(n*(phi1-phi2-phi3))>>    
  } // end of if(dMult>2)
  
  // V0 and ZDC
  fNonIsotropicTermsV0Pro[0]->Fill(centrality, TMath::Cos(fPsi2V0C));
  fNonIsotropicTermsV0Pro[1]->Fill(centrality, TMath::Sin(fPsi2V0C));
  fNonIsotropicTermsV0Pro[2]->Fill(centrality, TMath::Cos(fPsi2V0A));
  fNonIsotropicTermsV0Pro[3]->Fill(centrality, TMath::Sin(fPsi2V0A));
  fNonIsotropicTermsV0Pro[4]->Fill(centrality, TMath::Cos(2*fPsi2V0C));
  fNonIsotropicTermsV0Pro[5]->Fill(centrality, TMath::Sin(2*fPsi2V0C));
  fNonIsotropicTermsV0Pro[6]->Fill(centrality, TMath::Cos(2*fPsi2V0A));
  fNonIsotropicTermsV0Pro[7]->Fill(centrality, TMath::Sin(2*fPsi2V0A));
  
  fNonIsotropicTermsZDCPro[0]->Fill(centrality, TMath::Cos(fPsiZDCC));
  fNonIsotropicTermsZDCPro[1]->Fill(centrality, TMath::Sin(fPsiZDCC));
  fNonIsotropicTermsZDCPro[2]->Fill(centrality, TMath::Cos(fPsiZDCA));
  fNonIsotropicTermsZDCPro[3]->Fill(centrality, TMath::Sin(fPsiZDCA));
  fNonIsotropicTermsZDCPro[4]->Fill(centrality, TMath::Cos(fPsiZDCCA));
  fNonIsotropicTermsZDCPro[5]->Fill(centrality, TMath::Sin(fPsiZDCCA));
  fNonIsotropicTermsZDCPro[6]->Fill(centrality, TMath::Cos(2*fPsiZDCC));
  fNonIsotropicTermsZDCPro[7]->Fill(centrality, TMath::Sin(2*fPsiZDCC));
  fNonIsotropicTermsZDCPro[8]->Fill(centrality, TMath::Cos(2*fPsiZDCA));
  fNonIsotropicTermsZDCPro[9]->Fill(centrality, TMath::Sin(2*fPsiZDCA));
  fNonIsotropicTermsZDCPro[10]->Fill(centrality, TMath::Cos(2*fPsiZDCCA));
  fNonIsotropicTermsZDCPro[11]->Fill(centrality, TMath::Sin(2*fPsiZDCCA));
  

} // end of void AliAnalysisTaskGammaNonIsotropicCorr::CalculateNonIsotropicTerms()

//=======================================================================================================================

void AliAnalysisTaskGammaNonIsotropicCorr::CalculateDifferential3pCorrelator()
{
 // Calculate differential 3-p azimuthal correlator cos[n(psi1+psi2-2phi3)] in terms of Q_{2n}, p_{n}, q1_{n} and q2_{n}.
 
 // a) Calculate differential 3-p correlator without using particle weights;
 // b) Calculate non-isotropic terms for 3-p correlator.

 // a) Calculate differential 3-p correlator without using particle weights: 
   Int_t iBinCounter = 0;
   
  // Multiplicity (number of RPs):
  Double_t dMult = fCMEQReRP->GetBinEntries(1); // (*fSpk)(0,0);
  // Real and imaginary parts of non-weighted Q-vectors (Q_{n,0}) evaluated in harmonic 2n: 
  Double_t dReQ2n = fCMEQReRP->GetBinContent(2)*dMult; // (*fReQnk)(1,0);
  Double_t dImQ2n = fCMEQImRP->GetBinContent(2)*dMult; // (*fImQnk)(1,0);
  
  Double_t dReQ2nV0C = cos(2*fPsi2V0C);
  Double_t dImQ2nV0C = sin(2*fPsi2V0C);
  Double_t dReQ2nV0A = cos(2*fPsi2V0A);
  Double_t dImQ2nV0A = sin(2*fPsi2V0A);
  Double_t dReQ2nZDCC = cos(2*fPsiZDCC);
  Double_t dImQ2nZDCC = sin(2*fPsiZDCC);
  Double_t dReQ2nZDCA = cos(2*fPsiZDCA);
  Double_t dImQ2nZDCA = sin(2*fPsiZDCA);
  Double_t dReQ2nZDCCA = cos(2*fPsiZDCCA);
  Double_t dImQ2nZDCCA = sin(2*fPsiZDCCA);
  
  Double_t Intgp1nRe = 0, Intgp1nIm = 0, Intgmp = 0, Intgoverlap1 = 0, Intgoverlap2 = 0, Intgweight = 0;
  
  // real and imaginary parts of p_{n} OS, PP and NN: 
  Double_t p1nReOS = fRePEBEOS->GetBinContent(1)*fRePEBEOS->GetBinEntries(1); // <cos(dPsi1+dPsi2)>
  Double_t p1nImOS = fImPEBEOS->GetBinContent(1)*fImPEBEOS->GetBinEntries(1); // <sin(dPsi1+dPsi2)>
  
  Double_t p1nRePP = fRePEBEPP->GetBinContent(1)*fRePEBEPP->GetBinEntries(1); // <cos(dPsi1+dPsi2)>
  Double_t p1nImPP = fImPEBEPP->GetBinContent(1)*fImPEBEPP->GetBinEntries(1); // <sin(dPsi1+dPsi2)>
  Double_t p1nReNN = fRePEBENN->GetBinContent(1)*fRePEBENN->GetBinEntries(1); // <cos(dPsi1+dPsi2)>
  Double_t p1nImNN = fImPEBENN->GetBinContent(1)*fImPEBENN->GetBinEntries(1); // <sin(dPsi1+dPsi2)>
  
  // overlap 1: to be improved (terminology)
  Double_t overlap1OS = f2pCorrelatorCosPsiDiffOS->GetBinContent(1)*f2pCorrelatorCosPsiDiffOS->GetBinEntries(1); // POI 1 <cos(dPsi1-dPsi2)>
  Double_t overlap1PP = f2pCorrelatorCosPsiDiffPP->GetBinContent(1)*f2pCorrelatorCosPsiDiffPP->GetBinEntries(1); // POI 1 <cos(dPsi1-dPsi2)>
  Double_t overlap1NN = f2pCorrelatorCosPsiDiffNN->GetBinContent(1)*f2pCorrelatorCosPsiDiffNN->GetBinEntries(1); // POI 1 <cos(dPsi1-dPsi2)>
  // overlap 2: to be improved (terminology)
  Double_t overlap2OS = f2pCorrelatorCosPsiDiffOS->GetBinContent(1)*f2pCorrelatorCosPsiDiffOS->GetBinEntries(1); // POI 2 <cos(dPsi1-dPsi2)> same as POI1
  Double_t overlap2PP = f2pCorrelatorCosPsiDiffPP->GetBinContent(1)*f2pCorrelatorCosPsiDiffPP->GetBinEntries(1); // POI 2 <cos(dPsi1-dPsi2)> same as POI1
  Double_t overlap2NN = f2pCorrelatorCosPsiDiffNN->GetBinContent(1)*f2pCorrelatorCosPsiDiffNN->GetBinEntries(1); // POI 2 <cos(dPsi1-dPsi2)> same as POI1
  // number of pairs of POIs in particular (p1+p2)/2 or |p1-p2| bin:
  Double_t mpOS = f2pCorrelatorCosPsiDiffOS->GetBinEntries(1);
  Double_t mpPP = f2pCorrelatorCosPsiDiffPP->GetBinEntries(1);
  Double_t mpNN = f2pCorrelatorCosPsiDiffNN->GetBinEntries(1);
  // number of pairs of POI1/RP and POI2 in particular (p1+p2)/2 or |p1-p2| bin:
  Double_t mOverlap1OS = f2pCorrelatorCosPsiDiffOS->GetBinEntries(1); // same number as POI1 and POI2 are all RP except for charges
  Double_t mOverlap1PP = f2pCorrelatorCosPsiDiffPP->GetBinEntries(1); // same number as POI1 and POI2 are all RP except for charges
  Double_t mOverlap1NN = f2pCorrelatorCosPsiDiffNN->GetBinEntries(1); // same number as POI1 and POI2 are all RP except for charges
  // number of pairs of POI2/RP and POI1 in particular (p1+p2)/2 or |p1-p2| bin:
  Double_t mOverlap2OS = f2pCorrelatorCosPsiDiffOS->GetBinEntries(1);
  Double_t mOverlap2PP = f2pCorrelatorCosPsiDiffPP->GetBinEntries(1);
  Double_t mOverlap2NN = f2pCorrelatorCosPsiDiffNN->GetBinEntries(1);
  // e-b-e weight for cos[n(psi1+psi2-2phi3)]:
  Double_t weightOS = mpOS*dMult-mOverlap1OS-mOverlap2OS;  // if POI1, POI2 and RP selection is same except for charge, mp=mOverlap1=mOverlap2
  Double_t weightForV0andZDCOS = mpOS; // no overlap between psi1, psi2 and V0,ZDC. dMult for V0 and ZDC is just 1
  Double_t weightPP = mpPP*dMult-mOverlap1PP-mOverlap2PP;  // if POI1, POI2 and RP selection is same except for charge, mp=mOverlap1=mOverlap2
  Double_t weightForV0andZDCPP = mpPP; // no overlap between psi1, psi2 and V0,ZDC. dMult for V0 and ZDC is just 1
  Double_t weightNN = mpNN*dMult-mOverlap1NN-mOverlap2NN;  // if POI1, POI2 and RP selection is same except for charge, mp=mOverlap1=mOverlap2
  Double_t weightForV0andZDCNN = mpNN; // no overlap between psi1, psi2 and V0,ZDC. dMult for V0 and ZDC is just 1
    
  Double_t cosP2nphi1M1npsi2M1npsi2OS = 0; // cos[n(psi1+psi2-2phi3)]
  Double_t cosP2nphi1M1npsi2M1npsiV0COS = 0, cosP2nphi1M1npsi2M1npsiV0AOS = 0; // cos[n(psi1+psi2-2phi_V0)]
  Double_t cosP2nphi1M1npsi2M1npsiZDCCOS = 0, cosP2nphi1M1npsi2M1npsiZDCAOS = 0, cosP2nphi1M1npsi2M1npsiZDCCAOS = 0; // cos[n(psi1+psi2-2phi_ZDC)]
  Double_t cosP2nphi1M1npsi2M1npsi2PP = 0; // cos[n(psi1+psi2-2phi3)]
  Double_t cosP2nphi1M1npsi2M1npsiV0CPP = 0, cosP2nphi1M1npsi2M1npsiV0APP = 0; // cos[n(psi1+psi2-2phi_V0)]
  Double_t cosP2nphi1M1npsi2M1npsiZDCCPP = 0, cosP2nphi1M1npsi2M1npsiZDCAPP = 0, cosP2nphi1M1npsi2M1npsiZDCCAPP = 0; // cos[n(psi1+psi2-2phi_ZDC)]
  Double_t cosP2nphi1M1npsi2M1npsi2NN = 0; // cos[n(psi1+psi2-2phi3)]
  Double_t cosP2nphi1M1npsi2M1npsiV0CNN = 0, cosP2nphi1M1npsi2M1npsiV0ANN = 0; // cos[n(psi1+psi2-2phi_V0)]
  Double_t cosP2nphi1M1npsi2M1npsiZDCCNN = 0, cosP2nphi1M1npsi2M1npsiZDCANN = 0, cosP2nphi1M1npsi2M1npsiZDCCANN = 0; // cos[n(psi1+psi2-2phi_ZDC)]
    
  if(weightOS>0.)
  {
     cosP2nphi1M1npsi2M1npsi2OS = (p1nReOS*dReQ2n+p1nImOS*dImQ2n-overlap1OS-overlap2OS)/(weightOS);
     f3pCorrelatorTPCOSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsi2OS,weightOS);
     
     // Calculate differential 3-p correlator for V0
     cosP2nphi1M1npsi2M1npsiV0COS = (p1nReOS*dReQ2nV0C+p1nImOS*dImQ2nV0C)/(weightForV0andZDCOS);
     cosP2nphi1M1npsi2M1npsiV0AOS = (p1nReOS*dReQ2nV0A+p1nImOS*dImQ2nV0A)/(weightForV0andZDCOS);
     // Calculate differential 3-p correlator for ZDC
     cosP2nphi1M1npsi2M1npsiZDCCOS = (p1nReOS*dReQ2nZDCC+p1nImOS*dImQ2nZDCC)/(weightForV0andZDCOS);
     cosP2nphi1M1npsi2M1npsiZDCAOS = (p1nReOS*dReQ2nZDCA+p1nImOS*dImQ2nZDCA)/(weightForV0andZDCOS);
     cosP2nphi1M1npsi2M1npsiZDCCAOS = (p1nReOS*dReQ2nZDCCA+p1nImOS*dImQ2nZDCCA)/(weightForV0andZDCOS);
     
     f3pCorrelatorV0COSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0COS,weightForV0andZDCOS);
     f3pCorrelatorV0AOSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0AOS,weightForV0andZDCOS);
     f3pCorrelatorZDCCOSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCOS,weightForV0andZDCOS);
     f3pCorrelatorZDCAOSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCAOS,weightForV0andZDCOS);
     f3pCorrelatorZDCCAOSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCAOS,weightForV0andZDCOS);
  }
  if(weightPP>0.)
  {
     cosP2nphi1M1npsi2M1npsi2PP = (p1nRePP*dReQ2n+p1nImPP*dImQ2n-overlap1PP-overlap2PP)/(weightPP);
     f3pCorrelatorTPCPPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsi2PP,weightPP);
     f3pCorrelatorTPCSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsi2PP,weightPP);
     
     // Calculate differential 3-p correlator for V0
     cosP2nphi1M1npsi2M1npsiV0CPP = (p1nRePP*dReQ2nV0C+p1nImPP*dImQ2nV0C)/(weightForV0andZDCPP);
     cosP2nphi1M1npsi2M1npsiV0APP = (p1nRePP*dReQ2nV0A+p1nImPP*dImQ2nV0A)/(weightForV0andZDCPP);
     // Calculate differential 3-p correlator for ZDC
     cosP2nphi1M1npsi2M1npsiZDCCPP = (p1nRePP*dReQ2nZDCC+p1nImPP*dImQ2nZDCC)/(weightForV0andZDCPP);
     cosP2nphi1M1npsi2M1npsiZDCAPP = (p1nRePP*dReQ2nZDCA+p1nImPP*dImQ2nZDCA)/(weightForV0andZDCPP);
     cosP2nphi1M1npsi2M1npsiZDCCAPP = (p1nRePP*dReQ2nZDCCA+p1nImPP*dImQ2nZDCCA)/(weightForV0andZDCPP);
     
     f3pCorrelatorV0CPPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0CPP,weightForV0andZDCPP);
     f3pCorrelatorV0APPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0APP,weightForV0andZDCPP);
     f3pCorrelatorZDCCPPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCPP,weightForV0andZDCPP);
     f3pCorrelatorZDCAPPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCAPP,weightForV0andZDCPP);
     f3pCorrelatorZDCCAPPPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCAPP,weightForV0andZDCPP);
     
     f3pCorrelatorV0CSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0CPP,weightForV0andZDCPP);
     f3pCorrelatorV0ASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0APP,weightForV0andZDCPP);
     f3pCorrelatorZDCCSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCPP,weightForV0andZDCPP);
     f3pCorrelatorZDCASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCAPP,weightForV0andZDCPP);
     f3pCorrelatorZDCCASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCAPP,weightForV0andZDCPP);
  }
  if(weightNN>0.)
  {
     cosP2nphi1M1npsi2M1npsi2NN = (p1nReNN*dReQ2n+p1nImNN*dImQ2n-overlap1NN-overlap2NN)/(weightNN);
     f3pCorrelatorTPCNNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsi2NN,weightNN);
     f3pCorrelatorTPCSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsi2NN,weightNN);
     
     // Calculate differential 3-p correlator for V0
     cosP2nphi1M1npsi2M1npsiV0CNN = (p1nReNN*dReQ2nV0C+p1nImNN*dImQ2nV0C)/(weightForV0andZDCNN);
     cosP2nphi1M1npsi2M1npsiV0ANN = (p1nReNN*dReQ2nV0A+p1nImNN*dImQ2nV0A)/(weightForV0andZDCNN);
     // Calculate differential 3-p correlator for ZDC
     cosP2nphi1M1npsi2M1npsiZDCCNN = (p1nReNN*dReQ2nZDCC+p1nImNN*dImQ2nZDCC)/(weightForV0andZDCNN);
     cosP2nphi1M1npsi2M1npsiZDCANN = (p1nReNN*dReQ2nZDCA+p1nImNN*dImQ2nZDCA)/(weightForV0andZDCNN);
     cosP2nphi1M1npsi2M1npsiZDCCANN = (p1nReNN*dReQ2nZDCCA+p1nImNN*dImQ2nZDCCA)/(weightForV0andZDCNN);
     
     f3pCorrelatorV0CNNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0CNN,weightForV0andZDCNN);
     f3pCorrelatorV0ANNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0ANN,weightForV0andZDCNN);
     f3pCorrelatorZDCCNNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCNN,weightForV0andZDCNN);
     f3pCorrelatorZDCANNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCANN,weightForV0andZDCNN);
     f3pCorrelatorZDCCANNPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCANN,weightForV0andZDCNN);
     
     f3pCorrelatorV0CSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0CNN,weightForV0andZDCNN);
     f3pCorrelatorV0ASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiV0ANN,weightForV0andZDCNN);
     f3pCorrelatorZDCCSSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCNN,weightForV0andZDCNN);
     f3pCorrelatorZDCASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCANN,weightForV0andZDCNN);
     f3pCorrelatorZDCCASSPro->Fill(centrality,cosP2nphi1M1npsi2M1npsiZDCCANN,weightForV0andZDCNN);
  }
     
  // b) Calculate non-isotropic terms for 3-p correlator.
  // non-isotropic terms, 1st POI:
  Double_t p1nRePOI1OS = fCMEQReRP->GetBinContent(1)*fCMEQReRP->GetBinEntries(1); //Cos(dPsi1) POI1, since POI OS is same as RP
  Double_t p1nImPOI1OS = fCMEQImRP->GetBinContent(1)*fCMEQImRP->GetBinEntries(1); //Sin(dPsi1) POI1
  Double_t mpPOI1OS = fCMEQReRP->GetBinEntries(1);
  Double_t q1nRePOI1OS = fCMEQReRP->GetBinContent(1)*fCMEQReRP->GetBinEntries(1); //Cos(dPsi1) POI1&RP, same as above as POI and RP overlap completely
  Double_t q1nImPOI1OS = fCMEQImRP->GetBinContent(1)*fCMEQImRP->GetBinEntries(1); //Sin(dPsi1) POI1&RP
  Double_t mqPOI1OS = fCMEQReRP->GetBinEntries(1); 
  
  Double_t p1nRePOI1PP = fCMEQRePOIPos->GetBinContent(1)*fCMEQRePOIPos->GetBinEntries(1); //Cos(dPsi1) POI1
  Double_t p1nImPOI1PP = fCMEQImPOIPos->GetBinContent(1)*fCMEQImPOIPos->GetBinEntries(1); //Sin(dPsi1) POI1
  Double_t mpPOI1PP = fCMEQRePOIPos->GetBinEntries(1);
  Double_t q1nRePOI1PP = fCMEQRePOIPos->GetBinContent(1)*fCMEQRePOIPos->GetBinEntries(1); //Cos(dPsi1) POI1&RP, same as above as POI and RP overlap completely
  Double_t q1nImPOI1PP = fCMEQImPOIPos->GetBinContent(1)*fCMEQImPOIPos->GetBinEntries(1); //Sin(dPsi1) POI1&RP
  Double_t mqPOI1PP = fCMEQRePOIPos->GetBinEntries(1); 
  
  Double_t p1nRePOI1NN = fCMEQRePOINeg->GetBinContent(1)*fCMEQRePOINeg->GetBinEntries(1); //Cos(dPsi1) POI1
  Double_t p1nImPOI1NN = fCMEQImPOINeg->GetBinContent(1)*fCMEQImPOINeg->GetBinEntries(1); //Sin(dPsi1) POI1
  Double_t mpPOI1NN = fCMEQRePOINeg->GetBinEntries(1);
  Double_t q1nRePOI1NN = fCMEQRePOINeg->GetBinContent(1)*fCMEQRePOINeg->GetBinEntries(1); //Cos(dPsi1) POI1&RP, same as above as POI and RP overlap completely
  Double_t q1nImPOI1NN = fCMEQImPOINeg->GetBinContent(1)*fCMEQImPOINeg->GetBinEntries(1); //Sin(dPsi1) POI1&RP
  Double_t mqPOI1NN = fCMEQRePOINeg->GetBinEntries(1); 
  
  // non-isotropic terms, 2nd POI:
  Double_t p1nRePOI2OS = fCMEQReRP->GetBinContent(1)*fCMEQReRP->GetBinEntries(1); //Cos(dPsi1) POI2, same as POI1 if POI1 and POI2 has exactly same selection
  Double_t p1nImPOI2OS = fCMEQImRP->GetBinContent(1)*fCMEQImRP->GetBinEntries(1); //Sin(dPsi1) POI2
  Double_t mpPOI2OS = fCMEQReRP->GetBinEntries(1);
  Double_t q1nRePOI2OS = fCMEQReRP->GetBinContent(1)*fCMEQReRP->GetBinEntries(1); //Cos(dPsi1) POI2&RP
  Double_t q1nImPOI2OS = fCMEQImRP->GetBinContent(1)*fCMEQImRP->GetBinEntries(1); //Sin(dPsi1) POI2&RP
  Double_t mqPOI2OS = fCMEQReRP->GetBinEntries(1);
  
  Double_t p1nRePOI2PP = fCMEQRePOIPos->GetBinContent(1)*fCMEQRePOIPos->GetBinEntries(1); //Cos(dPsi1) POI2
  Double_t p1nImPOI2PP = fCMEQImPOIPos->GetBinContent(1)*fCMEQImPOIPos->GetBinEntries(1); //Sin(dPsi1) POI2
  Double_t mpPOI2PP = fCMEQRePOIPos->GetBinEntries(1);
  Double_t q1nRePOI2PP = fCMEQRePOIPos->GetBinContent(1)*fCMEQRePOIPos->GetBinEntries(1); //Cos(dPsi1) POI2&RP, same as above as POI and RP overlap completely
  Double_t q1nImPOI2PP = fCMEQImPOIPos->GetBinContent(1)*fCMEQImPOIPos->GetBinEntries(1); //Sin(dPsi1) POI2&RP
  Double_t mqPOI2PP = fCMEQRePOIPos->GetBinEntries(1); 
  
  Double_t p1nRePOI2NN = fCMEQRePOINeg->GetBinContent(1)*fCMEQRePOINeg->GetBinEntries(1); //Cos(dPsi1) POI2
  Double_t p1nImPOI2NN = fCMEQImPOINeg->GetBinContent(1)*fCMEQImPOINeg->GetBinEntries(1); //Sin(dPsi1) POI2
  Double_t mpPOI2NN = fCMEQRePOINeg->GetBinEntries(1);
  Double_t q1nRePOI2NN = fCMEQRePOINeg->GetBinContent(1)*fCMEQRePOINeg->GetBinEntries(1); //Cos(dPsi1) POI2&RP, same as above as POI and RP overlap completely
  Double_t q1nImPOI2NN = fCMEQImPOINeg->GetBinContent(1)*fCMEQImPOINeg->GetBinEntries(1); //Sin(dPsi1) POI2&RP
  Double_t mqPOI2NN = fCMEQRePOINeg->GetBinEntries(1); 
  
  // Fill all-event profiles:   
  if(weightOS>0. && mpPOI1OS>0.)
  { 
   fNonIsotropicTermsOSPro[0]->Fill(centrality,p1nRePOI1OS/mpPOI1OS,mpPOI1OS); // <<cos(#psi_{POI_1})>>
   fNonIsotropicTermsOSPro[1]->Fill(centrality,p1nImPOI1OS/mpPOI1OS,mpPOI1OS); // <<sin(#psi_{POI_1})>>
  }
  if(weightOS>0. && mpPOI2OS>0.)  
  {
   fNonIsotropicTermsOSPro[2]->Fill(centrality,p1nRePOI2OS/mpPOI2OS,mpPOI2OS); // <<cos(#psi_{POI_2})>>
   fNonIsotropicTermsOSPro[3]->Fill(centrality,p1nImPOI2OS/mpPOI2OS,mpPOI2OS); // <<sin(#psi_{POI_2})>>
  }
  if(weightOS>0. && mpPOI1OS*dMult-mqPOI1OS>0.)
  { 
   fNonIsotropicTermsOSPro[4]->Fill(centrality,(p1nRePOI1OS*dReQ2n+p1nImPOI1OS*dImQ2n-q1nRePOI1OS)/(mpPOI1OS*dMult-mqPOI1OS),mpPOI1OS*dMult-mqPOI1OS); // <<cos(#psi_{POI_1}-2*phi)>>
   fNonIsotropicTermsOSPro[5]->Fill(centrality, (p1nImPOI1OS*dReQ2n-p1nRePOI1OS*dImQ2n+q1nImPOI1OS)/(mpPOI1OS*dMult-mqPOI1OS),mpPOI1OS*dMult-mqPOI1OS); // <<sin(#psi_{POI_1}-2*phi)>>
   
   fNonIsotropicTermsV0OSPro[0]->Fill(centrality,(p1nRePOI1OS*dReQ2nV0C+p1nImPOI1OS*dImQ2nV0C)/(mpPOI1OS),mpPOI1OS); // <<cos(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0OSPro[1]->Fill(centrality, (p1nImPOI1OS*dReQ2nV0C-p1nRePOI1OS*dImQ2nV0C)/(mpPOI1OS),mpPOI1OS); // <<sin(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0OSPro[2]->Fill(centrality,(p1nRePOI1OS*dReQ2nV0A+p1nImPOI1OS*dImQ2nV0A)/(mpPOI1OS),mpPOI1OS); // <<cos(#psi_{POI_1}-2*phi_{V0A})>>
   fNonIsotropicTermsV0OSPro[3]->Fill(centrality, (p1nImPOI1OS*dReQ2nV0A-p1nRePOI1OS*dImQ2nV0A)/(mpPOI1OS),mpPOI1OS); // <<sin(#psi_{POI_1}-2*phi_{V0A})>>
   
   fNonIsotropicTermsZDCOSPro[0]->Fill(centrality,(p1nRePOI1OS*dReQ2nZDCC+p1nImPOI1OS*dImQ2nZDCC)/(mpPOI1OS),mpPOI1OS); // <<cos(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCOSPro[1]->Fill(centrality, (p1nImPOI1OS*dReQ2nZDCC-p1nRePOI1OS*dImQ2nZDCC)/(mpPOI1OS),mpPOI1OS); // <<sin(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCOSPro[2]->Fill(centrality,(p1nRePOI1OS*dReQ2nZDCA+p1nImPOI1OS*dImQ2nZDCA)/(mpPOI1OS),mpPOI1OS); // <<cos(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCOSPro[3]->Fill(centrality, (p1nImPOI1OS*dReQ2nZDCA-p1nRePOI1OS*dImQ2nZDCA)/(mpPOI1OS),mpPOI1OS); // <<sin(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCOSPro[4]->Fill(centrality,(p1nRePOI1OS*dReQ2nZDCCA+p1nImPOI1OS*dImQ2nZDCCA)/(mpPOI1OS),mpPOI1OS); // <<cos(#psi_{POI_1}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCOSPro[5]->Fill(centrality, (p1nImPOI1OS*dReQ2nZDCCA-p1nRePOI1OS*dImQ2nZDCCA)/(mpPOI1OS),mpPOI1OS); // <<sin(#psi_{POI_1}-2*phi_{ZDCCA})>>
  }
  if(weightOS>0. && mpPOI2OS*dMult-mqPOI2OS>0.)
  { 
   fNonIsotropicTermsOSPro[6]->Fill(centrality,(p1nRePOI2OS*dReQ2n+p1nImPOI2OS*dImQ2n-q1nRePOI2OS)/(mpPOI2OS*dMult-mqPOI2OS),mpPOI2OS*dMult-mqPOI2OS); // <<cos(#psi_{POI_2}-2*phi)>>
   fNonIsotropicTermsOSPro[7]->Fill(centrality,(p1nImPOI2OS*dReQ2n-p1nRePOI2OS*dImQ2n+q1nImPOI2OS)/(mpPOI2OS*dMult-mqPOI2OS),mpPOI2OS*dMult-mqPOI2OS); // <<sin(#psi_{POI_2}-2*phi)>>
   
   fNonIsotropicTermsV0OSPro[4]->Fill(centrality,(p1nRePOI2OS*dReQ2nV0C+p1nImPOI2OS*dImQ2nV0C)/(mpPOI2OS),mpPOI2OS); // <<cos(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0OSPro[5]->Fill(centrality, (p1nImPOI2OS*dReQ2nV0C-p1nRePOI2OS*dImQ2nV0C)/(mpPOI2OS),mpPOI2OS); // <<sin(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0OSPro[6]->Fill(centrality,(p1nRePOI2OS*dReQ2nV0A+p1nImPOI2OS*dImQ2nV0A)/(mpPOI2OS),mpPOI2OS); // <<cos(#psi_{POI_2}-2*phi_{V0A})>>
   fNonIsotropicTermsV0OSPro[7]->Fill(centrality, (p1nImPOI2OS*dReQ2nV0A-p1nRePOI2OS*dImQ2nV0A)/(mpPOI2OS),mpPOI2OS); // <<sin(#psi_{POI_2}-2*phi_{V0A})>>
   
   fNonIsotropicTermsZDCOSPro[6]->Fill(centrality,(p1nRePOI2OS*dReQ2nZDCC+p1nImPOI2OS*dImQ2nZDCC)/(mpPOI2OS),mpPOI2OS); // <<cos(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCOSPro[7]->Fill(centrality, (p1nImPOI2OS*dReQ2nZDCC-p1nRePOI2OS*dImQ2nZDCC)/(mpPOI2OS),mpPOI2OS); // <<sin(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCOSPro[8]->Fill(centrality,(p1nRePOI2OS*dReQ2nZDCA+p1nImPOI2OS*dImQ2nZDCA)/(mpPOI2OS),mpPOI2OS); // <<cos(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCOSPro[9]->Fill(centrality, (p1nImPOI2OS*dReQ2nZDCA-p1nRePOI2OS*dImQ2nZDCA)/(mpPOI2OS),mpPOI2OS); // <<sin(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCOSPro[10]->Fill(centrality,(p1nRePOI2OS*dReQ2nZDCCA+p1nImPOI2OS*dImQ2nZDCCA)/(mpPOI2OS),mpPOI2OS); // <<cos(#psi_{POI_2}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCOSPro[11]->Fill(centrality, (p1nImPOI2OS*dReQ2nZDCCA-p1nRePOI2OS*dImQ2nZDCCA)/(mpPOI2OS),mpPOI2OS); // <<sin(#psi_{POI_2}-2*phi_{ZDCCA})>>
  }
  if(weightOS>0. && mpOS>0.)
  {
   fNonIsotropicTermsOSPro[8]->Fill(centrality,p1nReOS/mpOS,mpOS); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   fNonIsotropicTermsOSPro[9]->Fill(centrality,p1nImOS/mpOS,mpOS); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
  }
  
  if(weightPP>0. && mpPOI1PP>0.)
  { 
   fNonIsotropicTermsPPPro[0]->Fill(centrality,p1nRePOI1PP/mpPOI1PP,mpPOI1PP); // <<cos(#psi_{POI_1})>>
   fNonIsotropicTermsPPPro[1]->Fill(centrality,p1nImPOI1PP/mpPOI1PP,mpPOI1PP); // <<sin(#psi_{POI_1})>>
   
   fNonIsotropicTermsSSPro[0]->Fill(centrality,p1nRePOI1PP/mpPOI1PP,mpPOI1PP); // <<cos(#psi_{POI_1})>>
   fNonIsotropicTermsSSPro[1]->Fill(centrality,p1nImPOI1PP/mpPOI1PP,mpPOI1PP); // <<sin(#psi_{POI_1})>>
  }
  if(weightPP>0. && mpPOI2PP>0.)  
  {
   fNonIsotropicTermsPPPro[2]->Fill(centrality,p1nRePOI2PP/mpPOI2PP,mpPOI2PP); // <<cos(#psi_{POI_2})>>
   fNonIsotropicTermsPPPro[3]->Fill(centrality,p1nImPOI2PP/mpPOI2PP,mpPOI2PP); // <<sin(#psi_{POI_2})>>
   
   fNonIsotropicTermsSSPro[2]->Fill(centrality,p1nRePOI2PP/mpPOI2PP,mpPOI2PP); // <<cos(#psi_{POI_2})>>
   fNonIsotropicTermsSSPro[3]->Fill(centrality,p1nImPOI2PP/mpPOI2PP,mpPOI2PP); // <<sin(#psi_{POI_2})>>
  }
  if(weightPP>0. && mpPOI1PP*dMult-mqPOI1PP>0.)
  { 
   // TPC
   fNonIsotropicTermsPPPro[4]->Fill(centrality,(p1nRePOI1PP*dReQ2n+p1nImPOI1PP*dImQ2n-q1nRePOI1PP)/(mpPOI1PP*dMult-mqPOI1PP),mpPOI1PP*dMult-mqPOI1PP); // <<cos(#psi_{POI_1}-2*phi)>>
   fNonIsotropicTermsPPPro[5]->Fill(centrality, (p1nImPOI1PP*dReQ2n-p1nRePOI1PP*dImQ2n+q1nImPOI1PP)/(mpPOI1PP*dMult-mqPOI1PP),mpPOI1PP*dMult-mqPOI1PP); // <<sin(#psi_{POI_1}-2*phi)>>
   
   fNonIsotropicTermsSSPro[4]->Fill(centrality,(p1nRePOI1PP*dReQ2n+p1nImPOI1PP*dImQ2n-q1nRePOI1PP)/(mpPOI1PP*dMult-mqPOI1PP),mpPOI1PP*dMult-mqPOI1PP); // <<cos(#psi_{POI_1}-2*phi)>>
   fNonIsotropicTermsSSPro[5]->Fill(centrality, (p1nImPOI1PP*dReQ2n-p1nRePOI1PP*dImQ2n+q1nImPOI1PP)/(mpPOI1PP*dMult-mqPOI1PP),mpPOI1PP*dMult-mqPOI1PP); // <<sin(#psi_{POI_1}-2*phi)>>
   
   // V0
   fNonIsotropicTermsV0PPPro[0]->Fill(centrality,(p1nRePOI1PP*dReQ2nV0C+p1nImPOI1PP*dImQ2nV0C)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0PPPro[1]->Fill(centrality, (p1nImPOI1PP*dReQ2nV0C-p1nRePOI1PP*dImQ2nV0C)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0PPPro[2]->Fill(centrality,(p1nRePOI1PP*dReQ2nV0A+p1nImPOI1PP*dImQ2nV0A)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{V0A})>>
   fNonIsotropicTermsV0PPPro[3]->Fill(centrality, (p1nImPOI1PP*dReQ2nV0A-p1nRePOI1PP*dImQ2nV0A)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{V0A})>>
   
   fNonIsotropicTermsV0SSPro[0]->Fill(centrality,(p1nRePOI1PP*dReQ2nV0C+p1nImPOI1PP*dImQ2nV0C)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[1]->Fill(centrality, (p1nImPOI1PP*dReQ2nV0C-p1nRePOI1PP*dImQ2nV0C)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[2]->Fill(centrality,(p1nRePOI1PP*dReQ2nV0A+p1nImPOI1PP*dImQ2nV0A)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{V0A})>>
   fNonIsotropicTermsV0SSPro[3]->Fill(centrality, (p1nImPOI1PP*dReQ2nV0A-p1nRePOI1PP*dImQ2nV0A)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{V0A})>>
   
   // ZDC
   fNonIsotropicTermsZDCPPPro[0]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCC+p1nImPOI1PP*dImQ2nZDCC)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCPPPro[1]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCC-p1nRePOI1PP*dImQ2nZDCC)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCPPPro[2]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCA+p1nImPOI1PP*dImQ2nZDCA)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCPPPro[3]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCA-p1nRePOI1PP*dImQ2nZDCA)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCPPPro[4]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCCA+p1nImPOI1PP*dImQ2nZDCCA)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCPPPro[5]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCCA-p1nRePOI1PP*dImQ2nZDCCA)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCCA})>>
   
   fNonIsotropicTermsZDCSSPro[0]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCC+p1nImPOI1PP*dImQ2nZDCC)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[1]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCC-p1nRePOI1PP*dImQ2nZDCC)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[2]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCA+p1nImPOI1PP*dImQ2nZDCA)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[3]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCA-p1nRePOI1PP*dImQ2nZDCA)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[4]->Fill(centrality,(p1nRePOI1PP*dReQ2nZDCCA+p1nImPOI1PP*dImQ2nZDCCA)/(mpPOI1PP),mpPOI1PP); // <<cos(#psi_{POI_1}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCSSPro[5]->Fill(centrality, (p1nImPOI1PP*dReQ2nZDCCA-p1nRePOI1PP*dImQ2nZDCCA)/(mpPOI1PP),mpPOI1PP); // <<sin(#psi_{POI_1}-2*phi_{ZDCCA})>>
  }
  if(weightPP>0. && mpPOI2PP*dMult-mqPOI2PP>0.)
  { 
   // TPC
   fNonIsotropicTermsPPPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2n+p1nImPOI2PP*dImQ2n-q1nRePOI2PP)/(mpPOI2PP*dMult-mqPOI2PP),mpPOI2PP*dMult-mqPOI2PP); // <<cos(#psi_{POI_2}-2*phi)>>
   fNonIsotropicTermsPPPro[7]->Fill(centrality,(p1nImPOI2PP*dReQ2n-p1nRePOI2PP*dImQ2n+q1nImPOI2PP)/(mpPOI2PP*dMult-mqPOI2PP),mpPOI2PP*dMult-mqPOI2PP); // <<sin(#psi_{POI_2}-2*phi)>>
   
   fNonIsotropicTermsSSPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2n+p1nImPOI2PP*dImQ2n-q1nRePOI2PP)/(mpPOI2PP*dMult-mqPOI2PP),mpPOI2PP*dMult-mqPOI2PP); // <<cos(#psi_{POI_2}-2*phi)>>
   fNonIsotropicTermsSSPro[7]->Fill(centrality,(p1nImPOI2PP*dReQ2n-p1nRePOI2PP*dImQ2n+q1nImPOI2PP)/(mpPOI2PP*dMult-mqPOI2PP),mpPOI2PP*dMult-mqPOI2PP); // <<sin(#psi_{POI_2}-2*phi)>>
   
   // V0
   fNonIsotropicTermsV0PPPro[4]->Fill(centrality,(p1nRePOI2PP*dReQ2nV0C+p1nImPOI2PP*dImQ2nV0C)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0PPPro[5]->Fill(centrality, (p1nImPOI2PP*dReQ2nV0C-p1nRePOI2PP*dImQ2nV0C)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0PPPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2nV0A+p1nImPOI2PP*dImQ2nV0A)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{V0A})>>
   fNonIsotropicTermsV0PPPro[7]->Fill(centrality, (p1nImPOI2PP*dReQ2nV0A-p1nRePOI2PP*dImQ2nV0A)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{V0A})>>
   
   fNonIsotropicTermsV0SSPro[4]->Fill(centrality,(p1nRePOI2PP*dReQ2nV0C+p1nImPOI2PP*dImQ2nV0C)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[5]->Fill(centrality, (p1nImPOI2PP*dReQ2nV0C-p1nRePOI2PP*dImQ2nV0C)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2nV0A+p1nImPOI2PP*dImQ2nV0A)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{V0A})>>
   fNonIsotropicTermsV0SSPro[7]->Fill(centrality, (p1nImPOI2PP*dReQ2nV0A-p1nRePOI2PP*dImQ2nV0A)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{V0A})>>
   
   // ZDC
   fNonIsotropicTermsZDCPPPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCC+p1nImPOI2PP*dImQ2nZDCC)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCPPPro[7]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCC-p1nRePOI2PP*dImQ2nZDCC)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCPPPro[8]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCA+p1nImPOI2PP*dImQ2nZDCA)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCPPPro[9]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCA-p1nRePOI2PP*dImQ2nZDCA)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCPPPro[10]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCCA+p1nImPOI2PP*dImQ2nZDCCA)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCPPPro[11]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCCA-p1nRePOI2PP*dImQ2nZDCCA)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCCA})>>
   
   fNonIsotropicTermsZDCSSPro[6]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCC+p1nImPOI2PP*dImQ2nZDCC)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[7]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCC-p1nRePOI2PP*dImQ2nZDCC)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[8]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCA+p1nImPOI2PP*dImQ2nZDCA)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[9]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCA-p1nRePOI2PP*dImQ2nZDCA)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[10]->Fill(centrality,(p1nRePOI2PP*dReQ2nZDCCA+p1nImPOI2PP*dImQ2nZDCCA)/(mpPOI2PP),mpPOI2PP); // <<cos(#psi_{POI_2}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCSSPro[11]->Fill(centrality, (p1nImPOI2PP*dReQ2nZDCCA-p1nRePOI2PP*dImQ2nZDCCA)/(mpPOI2PP),mpPOI2PP); // <<sin(#psi_{POI_2}-2*phi_{ZDCCA})>>
  }
  if(weightPP>0. && mpPP>0.)
  {
   fNonIsotropicTermsPPPro[8]->Fill(centrality,p1nRePP/mpPP,mpPP); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   fNonIsotropicTermsPPPro[9]->Fill(centrality,p1nImPP/mpPP,mpPP); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
   
   fNonIsotropicTermsSSPro[8]->Fill(centrality,p1nRePP/mpPP,mpPP); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   fNonIsotropicTermsSSPro[9]->Fill(centrality,p1nImPP/mpPP,mpPP); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
  }
  
  if(weightNN>0. && mpPOI1NN>0.)
  { 
   fNonIsotropicTermsNNPro[0]->Fill(centrality,p1nRePOI1NN/mpPOI1NN,mpPOI1NN); // <<cos(#psi_{POI_1})>>
   fNonIsotropicTermsNNPro[1]->Fill(centrality,p1nImPOI1NN/mpPOI1NN,mpPOI1NN); // <<sin(#psi_{POI_1})>>
   
   fNonIsotropicTermsSSPro[0]->Fill(centrality,p1nRePOI1NN/mpPOI1NN,mpPOI1NN); // <<cos(#psi_{POI_1})>>
   fNonIsotropicTermsSSPro[1]->Fill(centrality,p1nImPOI1NN/mpPOI1NN,mpPOI1NN); // <<sin(#psi_{POI_1})>>
  }
  if(weightNN>0. && mpPOI2NN>0.)  
  {
   fNonIsotropicTermsNNPro[2]->Fill(centrality,p1nRePOI2NN/mpPOI2NN,mpPOI2NN); // <<cos(#psi_{POI_2})>>
   fNonIsotropicTermsNNPro[3]->Fill(centrality,p1nImPOI2NN/mpPOI2NN,mpPOI2NN); // <<sin(#psi_{POI_2})>>
   
   fNonIsotropicTermsSSPro[2]->Fill(centrality,p1nRePOI2NN/mpPOI2NN,mpPOI2NN); // <<cos(#psi_{POI_2})>>
   fNonIsotropicTermsSSPro[3]->Fill(centrality,p1nImPOI2NN/mpPOI2NN,mpPOI2NN); // <<sin(#psi_{POI_2})>>
  }
  if(weightNN>0. && mpPOI1NN*dMult-mqPOI1NN>0.)
  { 
   // TPC
   fNonIsotropicTermsNNPro[4]->Fill(centrality,(p1nRePOI1NN*dReQ2n+p1nImPOI1NN*dImQ2n-q1nRePOI1NN)/(mpPOI1NN*dMult-mqPOI1NN),mpPOI1NN*dMult-mqPOI1NN); // <<cos(#psi_{POI_1}-2*phi)>>
   fNonIsotropicTermsNNPro[5]->Fill(centrality, (p1nImPOI1NN*dReQ2n-p1nRePOI1NN*dImQ2n+q1nImPOI1NN)/(mpPOI1NN*dMult-mqPOI1NN),mpPOI1NN*dMult-mqPOI1NN); // <<sin(#psi_{POI_1}-2*phi)>>
   
   fNonIsotropicTermsNNPro[4]->Fill(centrality,(p1nRePOI1NN*dReQ2n+p1nImPOI1NN*dImQ2n-q1nRePOI1NN)/(mpPOI1NN*dMult-mqPOI1NN),mpPOI1NN*dMult-mqPOI1NN); // <<cos(#psi_{POI_1}-2*phi)>>
   fNonIsotropicTermsNNPro[5]->Fill(centrality, (p1nImPOI1NN*dReQ2n-p1nRePOI1NN*dImQ2n+q1nImPOI1NN)/(mpPOI1NN*dMult-mqPOI1NN),mpPOI1NN*dMult-mqPOI1NN); // <<sin(#psi_{POI_1}-2*phi)>>
   
   // V0
   fNonIsotropicTermsV0NNPro[0]->Fill(centrality,(p1nRePOI1NN*dReQ2nV0C+p1nImPOI1NN*dImQ2nV0C)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[1]->Fill(centrality, (p1nImPOI1NN*dReQ2nV0C-p1nRePOI1NN*dImQ2nV0C)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[2]->Fill(centrality,(p1nRePOI1NN*dReQ2nV0A+p1nImPOI1NN*dImQ2nV0A)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{V0A})>>
   fNonIsotropicTermsV0NNPro[3]->Fill(centrality, (p1nImPOI1NN*dReQ2nV0A-p1nRePOI1NN*dImQ2nV0A)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{V0A})>>
   
   fNonIsotropicTermsV0NNPro[0]->Fill(centrality,(p1nRePOI1NN*dReQ2nV0C+p1nImPOI1NN*dImQ2nV0C)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[1]->Fill(centrality, (p1nImPOI1NN*dReQ2nV0C-p1nRePOI1NN*dImQ2nV0C)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[2]->Fill(centrality,(p1nRePOI1NN*dReQ2nV0A+p1nImPOI1NN*dImQ2nV0A)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{V0A})>>
   fNonIsotropicTermsV0NNPro[3]->Fill(centrality, (p1nImPOI1NN*dReQ2nV0A-p1nRePOI1NN*dImQ2nV0A)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{V0A})>>
   
   // ZDC
   fNonIsotropicTermsZDCNNPro[0]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCC+p1nImPOI1NN*dImQ2nZDCC)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[1]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCC-p1nRePOI1NN*dImQ2nZDCC)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[2]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCA+p1nImPOI1NN*dImQ2nZDCA)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[3]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCA-p1nRePOI1NN*dImQ2nZDCA)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[4]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCCA+p1nImPOI1NN*dImQ2nZDCCA)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCNNPro[5]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCCA-p1nRePOI1NN*dImQ2nZDCCA)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCCA})>>
   
   fNonIsotropicTermsZDCNNPro[0]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCC+p1nImPOI1NN*dImQ2nZDCC)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[1]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCC-p1nRePOI1NN*dImQ2nZDCC)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[2]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCA+p1nImPOI1NN*dImQ2nZDCA)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[3]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCA-p1nRePOI1NN*dImQ2nZDCA)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[4]->Fill(centrality,(p1nRePOI1NN*dReQ2nZDCCA+p1nImPOI1NN*dImQ2nZDCCA)/(mpPOI1NN),mpPOI1NN); // <<cos(#psi_{POI_1}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCNNPro[5]->Fill(centrality, (p1nImPOI1NN*dReQ2nZDCCA-p1nRePOI1NN*dImQ2nZDCCA)/(mpPOI1NN),mpPOI1NN); // <<sin(#psi_{POI_1}-2*phi_{ZDCCA})>>
  }
  if(weightNN>0. && mpPOI2NN*dMult-mqPOI2NN>0.)
  { 
   // TPC
   fNonIsotropicTermsNNPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2n+p1nImPOI2NN*dImQ2n-q1nRePOI2NN)/(mpPOI2NN*dMult-mqPOI2NN),mpPOI2NN*dMult-mqPOI2NN); // <<cos(#psi_{POI_2}-2*phi)>>
   fNonIsotropicTermsNNPro[7]->Fill(centrality,(p1nImPOI2NN*dReQ2n-p1nRePOI2NN*dImQ2n+q1nImPOI2NN)/(mpPOI2NN*dMult-mqPOI2NN),mpPOI2NN*dMult-mqPOI2NN); // <<sin(#psi_{POI_2}-2*phi)>>
   
   fNonIsotropicTermsSSPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2n+p1nImPOI2NN*dImQ2n-q1nRePOI2NN)/(mpPOI2NN*dMult-mqPOI2NN),mpPOI2NN*dMult-mqPOI2NN); // <<cos(#psi_{POI_2}-2*phi)>>
   fNonIsotropicTermsSSPro[7]->Fill(centrality,(p1nImPOI2NN*dReQ2n-p1nRePOI2NN*dImQ2n+q1nImPOI2NN)/(mpPOI2NN*dMult-mqPOI2NN),mpPOI2NN*dMult-mqPOI2NN); // <<sin(#psi_{POI_2}-2*phi)>>
   
   // V0
   fNonIsotropicTermsV0NNPro[4]->Fill(centrality,(p1nRePOI2NN*dReQ2nV0C+p1nImPOI2NN*dImQ2nV0C)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[5]->Fill(centrality, (p1nImPOI2NN*dReQ2nV0C-p1nRePOI2NN*dImQ2nV0C)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0NNPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2nV0A+p1nImPOI2NN*dImQ2nV0A)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{V0A})>>
   fNonIsotropicTermsV0NNPro[7]->Fill(centrality, (p1nImPOI2NN*dReQ2nV0A-p1nRePOI2NN*dImQ2nV0A)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{V0A})>>
   
   fNonIsotropicTermsV0SSPro[4]->Fill(centrality,(p1nRePOI2NN*dReQ2nV0C+p1nImPOI2NN*dImQ2nV0C)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[5]->Fill(centrality, (p1nImPOI2NN*dReQ2nV0C-p1nRePOI2NN*dImQ2nV0C)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{V0C})>>
   fNonIsotropicTermsV0SSPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2nV0A+p1nImPOI2NN*dImQ2nV0A)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{V0A})>>
   fNonIsotropicTermsV0SSPro[7]->Fill(centrality, (p1nImPOI2NN*dReQ2nV0A-p1nRePOI2NN*dImQ2nV0A)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{V0A})>>
   
   // ZDC
   fNonIsotropicTermsZDCNNPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCC+p1nImPOI2NN*dImQ2nZDCC)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[7]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCC-p1nRePOI2NN*dImQ2nZDCC)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCNNPro[8]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCA+p1nImPOI2NN*dImQ2nZDCA)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[9]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCA-p1nRePOI2NN*dImQ2nZDCA)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCNNPro[10]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCCA+p1nImPOI2NN*dImQ2nZDCCA)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCNNPro[11]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCCA-p1nRePOI2NN*dImQ2nZDCCA)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCCA})>>
   
   fNonIsotropicTermsZDCSSPro[6]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCC+p1nImPOI2NN*dImQ2nZDCC)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[7]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCC-p1nRePOI2NN*dImQ2nZDCC)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCC})>>
   fNonIsotropicTermsZDCSSPro[8]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCA+p1nImPOI2NN*dImQ2nZDCA)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[9]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCA-p1nRePOI2NN*dImQ2nZDCA)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCA})>>
   fNonIsotropicTermsZDCSSPro[10]->Fill(centrality,(p1nRePOI2NN*dReQ2nZDCCA+p1nImPOI2NN*dImQ2nZDCCA)/(mpPOI2NN),mpPOI2NN); // <<cos(#psi_{POI_2}-2*phi_{ZDCCA})>>
   fNonIsotropicTermsZDCSSPro[11]->Fill(centrality, (p1nImPOI2NN*dReQ2nZDCCA-p1nRePOI2NN*dImQ2nZDCCA)/(mpPOI2NN),mpPOI2NN); // <<sin(#psi_{POI_2}-2*phi_{ZDCCA})>>
  }
  if(weightNN>0. && mpNN>0.)
  {
   fNonIsotropicTermsNNPro[8]->Fill(centrality,p1nReNN/mpNN,mpNN); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   fNonIsotropicTermsNNPro[9]->Fill(centrality,p1nImNN/mpNN,mpNN); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
   
   fNonIsotropicTermsSSPro[8]->Fill(centrality,p1nReNN/mpNN,mpNN); // <<cos(#psi_{POI_1}+#psi_{POI_2})>>
   fNonIsotropicTermsSSPro[9]->Fill(centrality,p1nImNN/mpNN,mpNN); // <<sin(#psi_{POI_1}+#psi_{POI_2})>>   
  }

} // end of void AliAnalysisTaskGammaNonIsotropicCorr::CalculateDifferential3pCorrelator() 


// ==========================================================================================================

void AliAnalysisTaskGammaNonIsotropicCorr::SetupQAHistograms(){

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
  fDebugwEventCount->GetXaxis()->SetBinLabel(6,"V0Mult>0.0001");
  fDebugwEventCount->GetXaxis()->SetBinLabel(7,"ZDC towE>0");
  fDebugwEventCount->GetXaxis()->SetBinLabel(8,"denZNCA>0");
  fDebugwEventCount->GetXaxis()->SetBinLabel(9,"Analyzed");
  
  /// 2D QA:
  
  ///TPC Event Planes:
  //fHistTPCPsiNPosPlane = new TH2F("fHistTPCPsiNPosPlane",Form("#Psi_{n}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsiNPosPlane);
  //fHistTPCPsiNNegPlane = new TH2F("fHistTPCPsiNNegPlane",Form("#Psi_{n}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsiNNegPlane);
  //fHistTPCPsi3PosPlane = new TH2F("fHistTPCPsi3PosPlane",Form("#Psi_{3}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",3), 18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsi3PosPlane);
  //fHistTPCPsi3NegPlane = new TH2F("fHistTPCPsi3NegPlane",Form("#Psi_{3}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",3), 18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsi3NegPlane);
  //fHistTPCPsi4PosPlane = new TH2F("fHistTPCPsi4PosPlane",Form("#Psi_{3}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",3), 18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsi4PosPlane);
  //fHistTPCPsi4NegPlane = new TH2F("fHistTPCPsi4NegPlane",Form("#Psi_{3}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",3), 18,0,90,50,0,3.14159);
  //fListHist->Add(fHistTPCPsi4NegPlane);  
  
  /// V0 Event Planes:
  fHistV0CPsiNEventPlane = new TH2F("fHistV0CPsiNEventPlane",Form("#Psi_{n}(V0C); centrality; #Psi_{%d,V0C}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0CPsiNEventPlane);
  fHistV0APsiNEventPlane = new TH2F("fHistV0APsiNEventPlane",Form("#Psi_{n}(V0A); centrality; #Psi_{%d,V0A}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0APsiNEventPlane);  
 
  ///q-Vector for ESE:
  //fHistTPCPosqVectorvsCent = new TH2F("fHistTPCPosqVectorvsCent","q TPC pos vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  //fListHist->Add(fHistTPCPosqVectorvsCent);
  //fHistTPCNegqVectorvsCent = new TH2F("fHistTPCNegqVectorvsCent","q TPC neg vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  //fListHist->Add(fHistTPCNegqVectorvsCent);
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
  ///  <Q> for the V0A
  hAvgQNXvsCentV0A = new TProfile("hAvgQNXvsCentV0A"," V0A <Q2x> vs Cent; Cent; <Q_{2,x}>",90,0,90);
  fListHist->Add(hAvgQNXvsCentV0A);
  hAvgQNYvsCentV0A = new TProfile("hAvgQNYvsCentV0A"," V0A <Q2y> vs Cent; Cent; <Q_{2,y}>",90,0,90);
  fListHist->Add(hAvgQNYvsCentV0A);
  
  /// Now <Q> for the TPC..  Note:gHarmonic is set as global variable!
  //fAvgCos2PsivsCentEtaPos = new TProfile("fAvgCos2PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  //fListHist->Add(fAvgCos2PsivsCentEtaPos);
  //fAvgSin2PsivsCentEtaPos = new TProfile("fAvgSin2PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  //fListHist->Add(fAvgSin2PsivsCentEtaPos);
  //fAvgCos2PsivsCentEtaNeg = new TProfile("fAvgCos2PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  //fListHist->Add(fAvgCos2PsivsCentEtaNeg);
  //fAvgSin2PsivsCentEtaNeg = new TProfile("fAvgSin2PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  //fListHist->Add(fAvgSin2PsivsCentEtaNeg);

  //fAvgCos3PsivsCentEtaPos = new TProfile("fAvgCos3PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  //fListHist->Add(fAvgCos3PsivsCentEtaPos);
  //fAvgSin3PsivsCentEtaPos = new TProfile("fAvgSin3PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  //fListHist->Add(fAvgSin3PsivsCentEtaPos);
  //fAvgCos3PsivsCentEtaNeg = new TProfile("fAvgCos3PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  //fListHist->Add(fAvgCos3PsivsCentEtaNeg);
  //fAvgSin3PsivsCentEtaNeg = new TProfile("fAvgSin3PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  //fListHist->Add(fAvgSin3PsivsCentEtaNeg);  

  //fAvgCos4PsivsCentEtaPos = new TProfile("fAvgCos4PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",4),90,0,90);
  //fListHist->Add(fAvgCos4PsivsCentEtaPos);
  //fAvgSin4PsivsCentEtaPos = new TProfile("fAvgSin4PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",4),90,0,90);
  //fListHist->Add(fAvgSin4PsivsCentEtaPos);
  //fAvgCos4PsivsCentEtaNeg = new TProfile("fAvgCos4PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",4),90,0,90);
  //fListHist->Add(fAvgCos4PsivsCentEtaNeg);
  //fAvgSin4PsivsCentEtaNeg = new TProfile("fAvgSin4PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",4),90,0,90);
  //fListHist->Add(fAvgSin4PsivsCentEtaNeg);

  /// Event by event quantities
  fCMEQReRP = new TProfile("fCMEQReRP", "fCMEQReRP",2,0,2);
  fTempList->Add(fCMEQReRP);
  
  fCMEQImRP = new TProfile("fCMEQImRP", "fCMEQImRP",2,0,2);
  fTempList->Add(fCMEQImRP);
  
  f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent", "f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent",2,0,2);
  fTempList->Add(f2pCorrelatorCos2PsiDiff2PsiV0RPPerEvent);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent", "f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent",3,0,3);
  fTempList->Add(f2pCorrelatorCos2PsiDiff2PsiZDCRPPerEvent);
  
  f2pCorrelatorCosPsiDiffPerEvent = new TProfile("f2pCorrelatorCosPsiDiffPerEvent", "f2pCorrelatorCosPsiDiffPerEvent",1,0,1);
  fTempList->Add(f2pCorrelatorCosPsiDiffPerEvent);
  
  f2pCorrelatorCos2PsiDiffPerEvent = new TProfile("f2pCorrelatorCos2PsiDiffPerEvent", "f2pCorrelatorCos2PsiDiffPerEvent",1,0,1);
  fTempList->Add(f2pCorrelatorCos2PsiDiffPerEvent);
  
  fCMEQRePOIPos = new TProfile("fCMEQRePOIPos", "fCMEQRePOIPos",2,0,2);
  fTempList->Add(fCMEQRePOIPos);
  
  fCMEQImPOIPos = new TProfile("fCMEQImPOIPos", "fCMEQImPOIPos",2,0,2);
  fTempList->Add(fCMEQImPOIPos);
  
  fCMEQRePOINeg = new TProfile("fCMEQRePOINeg", "fCMEQRePOINeg",2,0,2);
  fTempList->Add(fCMEQRePOINeg);
  
  fCMEQImPOINeg = new TProfile("fCMEQImPOINeg", "fCMEQImPOINeg",2,0,2);
  fTempList->Add(fCMEQImPOINeg);
  
  fNITV0OS = new TProfile("fNITV0OS", "fNITV0OS",4,0,4);
  fTempList->Add(fNITV0OS);
  
  fNITZDCOS = new TProfile("fNITZDCOS", "fNITZDCOS",6,0,6);
  fTempList->Add(fNITZDCOS);
  
  fNITV0POIPos = new TProfile("fNITV0POIPos", "fNITV0POIPos",4,0,4);
  fTempList->Add(fNITV0POIPos);
  
  fNITZDCPOIPos = new TProfile("fNITZDCPOIPos", "fNITZDCPOIPos",6,0,6);
  fTempList->Add(fNITZDCPOIPos);

  fNITV0POINeg = new TProfile("fNITV0POINeg", "fNITV0POINeg",4,0,4);
  fTempList->Add(fNITV0POINeg);
  
  fNITZDCPOINeg = new TProfile("fNITZDCPOINeg", "fNITZDCPOINeg",6,0,6);
  fTempList->Add(fNITZDCPOINeg);

  
  
  
  fRePEBEOS = new TProfile("fRePEBEOS", "fRePEBEOS",1,0,1);
  fTempList->Add(fRePEBEOS);
  
  fImPEBEOS = new TProfile("fImPEBEOS", "fImPEBEOS",1,0,1);
  fTempList->Add(fImPEBEOS);
  
  f2pCorrelatorCosPsiDiffOS = new TProfile("f2pCorrelatorCosPsiDiffOS", "f2pCorrelatorCosPsiDiffOS",1,0,1);
  fTempList->Add(f2pCorrelatorCosPsiDiffOS);
  
  f2pCorrelatorCos2PsiDiffOS = new TProfile("f2pCorrelatorCos2PsiDiffOS", "f2pCorrelatorCos2PsiDiffOS",1,0,1);
  fTempList->Add(f2pCorrelatorCos2PsiDiffOS);
  
  fRePEBEPP = new TProfile("fRePEBEPP", "fRePEBEPP",1,0,1);
  fTempList->Add(fRePEBEPP);
  
  fImPEBEPP = new TProfile("fImPEBEPP", "fImPEBEPP",1,0,1);
  fTempList->Add(fImPEBEPP);
  
  f2pCorrelatorCosPsiDiffPP = new TProfile("f2pCorrelatorCosPsiDiffPP", "f2pCorrelatorCosPsiDiffPP",1,0,1);
  fTempList->Add(f2pCorrelatorCosPsiDiffPP);
  
  f2pCorrelatorCos2PsiDiffPP = new TProfile("f2pCorrelatorCos2PsiDiffPP", "f2pCorrelatorCos2PsiDiffPP",1,0,1);
  fTempList->Add(f2pCorrelatorCos2PsiDiffPP);
  
  fRePEBENN = new TProfile("fRePEBENN", "fRePEBENN",1,0,1);
  fTempList->Add(fRePEBENN);
  
  fImPEBENN = new TProfile("fImPEBENN", "fImPEBENN",1,0,1);
  fTempList->Add(fImPEBENN);
  
  f2pCorrelatorCosPsiDiffNN = new TProfile("f2pCorrelatorCosPsiDiffNN", "f2pCorrelatorCosPsiDiffNN",1,0,1);
  fTempList->Add(f2pCorrelatorCosPsiDiffNN);
  
  f2pCorrelatorCos2PsiDiffNN = new TProfile("f2pCorrelatorCos2PsiDiffNN", "f2pCorrelatorCos2PsiDiffNN",1,0,1);
  fTempList->Add(f2pCorrelatorCos2PsiDiffNN);
 
}





void AliAnalysisTaskGammaNonIsotropicCorr::SetupPileUpRemovalFunctions18qPass3() { //@Shi for 2018 period Pass3 data
	// 18q pass3
	fSPDCutPU = new TF1("fSPDCutPU", "480. + 3.95*x", 0, 50000);
   
    Double_t parV0[8] = {41.3226, 0.822835, 0.0880984, 206.961, 3.56337, 0.0965816, -0.00076483, 2.11591e-06};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);
   
    Double_t parV0CL0[6] = {0.362458, 0.962768, 0.995134, 0.0331353, -0.000692428, 6.59962e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);
   
    Double_t parFB32[9] = {-812.555, 6.38397, 5379.01, -0.394814, 0.0296228, -26.1633, 317.365, -0.842175, 0.0165651};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
	
}

void AliAnalysisTaskGammaNonIsotropicCorr::SetupPileUpRemovalFunctions18rPass3() { //@Shi for 2018 period Pass3 data
    // 18r pass3

    fSPDCutPU = new TF1("fSPDCutPU", "480. + 3.95*x", 0, 50000);
   
    Double_t parV0[8] = {42.4921, 0.823255, 0.0824939, 139.826, 7.27032, 0.0488425, -0.00045769, 1.40891e-06};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);
   
    Double_t parV0CL0[6] = {0.317973, 0.961823, 1.02383, 0.0330231, -0.000721551, 6.92564e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);
   
    Double_t parFB32[9] = {-817.169, 6.40836, 5380.3, -0.394358, 0.0295209, -25.9573, 316.586, -0.843951, 0.0165442};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
	
}

void AliAnalysisTaskGammaNonIsotropicCorr::SetupEventAndTaskConfigInfo(){

 
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
  if(0){
    fHistAnalysisInfo->SetBinContent(7,1);
  }  
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(8,"NUE Applied");
  if(0){
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



void AliAnalysisTaskGammaNonIsotropicCorr::GetV0MCorrectionHist(Int_t run){ 

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

void AliAnalysisTaskGammaNonIsotropicCorr::GetZDCCorrectionHist(Int_t run){ 
  if(fListZDCCorr){
	fListZDCCorrRunByRun = (TList*) fListZDCCorr->FindObject(Form("Run %d", run));
	fHZDCCparameters = (TH1D*) fListZDCCorrRunByRun->FindObject(Form("fZDCCparameters[%d]",run));
	fHZDCAparameters = (TH1D*) fListZDCCorrRunByRun->FindObject(Form("fZDCAparameters[%d]",run));
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

Bool_t AliAnalysisTaskGammaNonIsotropicCorr::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){
  
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



Bool_t AliAnalysisTaskGammaNonIsotropicCorr::CheckEventIsPileUp2018(AliAODEvent *faod) {

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




void  AliAnalysisTaskGammaNonIsotropicCorr::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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



Bool_t AliAnalysisTaskGammaNonIsotropicCorr::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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



Bool_t AliAnalysisTaskGammaNonIsotropicCorr::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A, Double_t &sumMultV0C, Double_t &sumMultV0A){

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

void AliAnalysisTaskGammaNonIsotropicCorr::Terminate(Option_t *)  {
  //fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  //if (!fOutputList) return;
}




