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
/* $Id: AliAnalysisTaskGammaDeltaPIDSaveQvecSimple.cxx ver: 2.0                     $   */
/* Simple Task to fill V0 and ZDC Energies for Gain Calibration           */
/* Works with 15o and 18q/r. Support to be added for LHC10h               */
/* Developer: Md Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)        */
/* Last Modified: Aug 23, 2021,  First version committed                  */
/* Last Modified: Oct 08, 2021,  Second version committed                 */
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskGammaDeltaPIDSaveQvecSimple.h"
#include "AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple.h"
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

ClassImp(AliAnalysisTaskGammaDeltaPIDSaveQvecSimple)

AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::AliAnalysisTaskGammaDeltaPIDSaveQvecSimple(const char *name):
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
  fEventList(NULL),
  fTempList(NULL),
  treeEvent(NULL),
  fpQvecEvent(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fListZDCCorr(NULL),
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
  fCMEQReRP = NULL;
  fCMEQImRP = NULL;
  fCMEQRePOIPos = NULL;
  fCMEQImPOIPos = NULL;
  fCMEQRePOINeg = NULL;
  fCMEQImPOINeg = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0RP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCRP = NULL;
  fNITV0OS = NULL;
  fNITZDCOS = NULL;
  fNITV0POIPos = NULL;
  fNITZDCPOIPos = NULL;
  fNITV0POINeg = NULL;
  fNITZDCPOINeg = NULL;

  f2pCorrelatorCosPsiDiff = NULL;
  f2pCorrelatorCos2PsiDiff = NULL;
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
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::AliAnalysisTaskGammaDeltaPIDSaveQvecSimple():
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
  fEventList(NULL),
  fTempList(NULL),
  treeEvent(NULL),
  fpQvecEvent(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),   
  fListZDCCorr(NULL),
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
  fCMEQReRP = NULL;
  fCMEQImRP = NULL;
  fCMEQRePOIPos = NULL;
  fCMEQImPOIPos = NULL;
  fCMEQRePOINeg = NULL;
  fCMEQImPOINeg = NULL;
  f2pCorrelatorCos2PsiDiff2PsiV0RP = NULL;
  f2pCorrelatorCos2PsiDiff2PsiZDCRP = NULL;
  fNITV0OS = NULL;
  fNITZDCOS = NULL;
  fNITV0POIPos = NULL;
  fNITZDCPOIPos = NULL;
  fNITV0POINeg = NULL;
  fNITZDCPOINeg = NULL;

  f2pCorrelatorCosPsiDiff = NULL;
  f2pCorrelatorCos2PsiDiff = NULL;
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
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::~AliAnalysisTaskGammaDeltaPIDSaveQvecSimple()
{
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!
  if(fTempList)      delete fTempList;
  if(fListHist)      delete fListHist;  
  if(fListTRKCorr)   delete fListTRKCorr;
  if(fListNUACorr)   delete fListNUACorr;
  if(fListV0MCorr)   delete fListV0MCorr;
  if(fListZDCCorr)   delete fListZDCCorr;

  if(fV0CutPU)      delete fV0CutPU;
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU;   
}


//________________ Define Histograms _______________
void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::UserCreateOutputObjects()
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

  SetupQvecSavingObjects();

 
    
  cout<<"Ncls: "<<fTPCclustMin<<" TPCsharedCut: "<<fTPCsharedCut<<" UseTPCCrossedRows: "<<bUseTPCCrossedRows<<" Harm: "<<gHarmonic<<" POI: "<<gParticleID<<" nsigTPC: "<<fNSigmaTPCCut<<" nsigCirc: "<<fNSigmaTOFCut<<endl;
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
void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::UserExec(Option_t*) {
 
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
  
  fpQvecEvent->setOrbitNumber(orbit);
  
  Double_t fOrbitNumber = static_cast<double>(orbit)/1000000.;
  
  Bool_t kPileupEvent = kFALSE;

  //Double_t fSumQnxNeg = 0, fSumQnyNeg = 0, fSumQnxPos = 0, fSumQnyPos = 0;
  //Double_t fSumQnxNeg[3] = {0,}; // Array: Q2x, Q3x, Q4x,...
  //Double_t fSumQnyNeg[3] = {0,}; // Array: Q2y, Q3y, Q4y,... 
  //Double_t fSumQnxPos[3] = {0,};
  //Double_t fSumQnyPos[3] = {0,};  
  
  //Double_t fMultNeg = 0, fMultPos = 0;

  kPileupEvent = CheckEventIsPileUp2018(fAOD);

  if(kPileupEvent) return;  // If not a PileUp event, then We have TPC q vectors for EP.
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

  Double_t fPsi2V0C = 0., fPsi2V0A = 0.;
    
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
  
  //=============== Get the ZDC data ==================== @Shi

  Double_t fQxZNCC=0, fQyZNCC=0, fQxZNCA=0, fQyZNCA=0; 
  Double_t fPsiZNCC = 0., fPsiZNCA = 0., fPsiZNCCA = 0; // fPsiZNCCA combine ZNCC and ZNCA
  
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
	
	Double_t ZDCCxPosFromLogWeight = numXZNC/denZNC;
	Double_t ZDCCyPosFromLogWeight = numYZNC/denZNC;
	Double_t ZDCAxPosFromLogWeight = numXZNA/denZNA;
	Double_t ZDCAyPosFromLogWeight = numYZNA/denZNA;
	
	Double_t ZDCCAvgxPosFromVtxFit = 0;
	Double_t ZDCCAvgyPosFromVtxFit = 0;
	
	Double_t ZDCAAvgxPosFromVtxFit = 0;
	Double_t ZDCAAvgyPosFromVtxFit = 0;

	// 3rd order centrality+vtxpos+orbitNum
	ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pow(centrality,2) + fHZDCCparameters->GetBinContent(8)*pow(centrality,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12)*fOrbitNumber + fHZDCCparameters->GetBinContent(13);
	ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(14)*centrality + fHZDCCparameters->GetBinContent(15)*pow(centrality,2) + fHZDCCparameters->GetBinContent(16)*pow(centrality,3) + fHZDCCparameters->GetBinContent(17)*pVtxX + fHZDCCparameters->GetBinContent(18)*pVtxY + fHZDCCparameters->GetBinContent(19)*pVtxZ + fHZDCCparameters->GetBinContent(20)*fOrbitNumber + fHZDCCparameters->GetBinContent(21);
	
	ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pow(centrality,2) + fHZDCAparameters->GetBinContent(8)*pow(centrality,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12)*fOrbitNumber + fHZDCAparameters->GetBinContent(13);
	ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(14)*centrality + fHZDCAparameters->GetBinContent(15)*pow(centrality,2) + fHZDCAparameters->GetBinContent(16)*pow(centrality,3) + fHZDCAparameters->GetBinContent(17)*pVtxX + fHZDCAparameters->GetBinContent(18)*pVtxY + fHZDCAparameters->GetBinContent(19)*pVtxZ + fHZDCAparameters->GetBinContent(20)*fOrbitNumber + fHZDCAparameters->GetBinContent(21);
	
	
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
      if (fDCAxyMax>0 && fDCAzzMax>0 && trk1DCAxy>=fDCAxyMax && trk1DCAz>=fDCAzzMax) { // if fDCAxyMax, fDCAzzMax is set to be less than 0, no cut applied
		continue;
	  }
	  
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
		f2pCorrelatorCos2PsiDiff2PsiV0RP->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C))); //<cos(2psi1-2phi_V0C)>
		f2pCorrelatorCos2PsiDiff2PsiV0RP->Fill(1.5, TMath::Cos(2*(trk1Phi-fPsi2V0A))); //<cos(2psi1-2phi_V0A)>
		// Calculate <<cos(2a-2Psi_ZDC)>> for RP, POI OS
		f2pCorrelatorCos2PsiDiff2PsiZDCRP->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCC))); //<cos(2psi1-2phi_ZDCC)>
		f2pCorrelatorCos2PsiDiff2PsiZDCRP->Fill(1.5, TMath::Cos(2*(trk1Phi-fPsiZNCA))); //<cos(2psi1-2phi_ZDCA)>
		f2pCorrelatorCos2PsiDiff2PsiZDCRP->Fill(2.5, TMath::Cos(2*(trk1Phi-fPsiZNCCA))); //<cos(2psi1-2phi_ZDCCA)>
	
		fNITV0OS->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		fNITV0OS->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		fNITV0OS->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		fNITV0OS->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		fNITZDCOS->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZNCC)); // <<cos(psi1-2phi_ZDCC)>>
		fNITZDCOS->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZNCC)); // <<sin(psi1-2phi_ZDCC)>>
		fNITZDCOS->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZNCA)); // <<cos(psi1-2phi_ZDCA)>>
		fNITZDCOS->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZNCA)); // <<sin(psi1-2phi_ZDCA)>>
		fNITZDCOS->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZNCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		fNITZDCOS->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZNCCA)); // <<sin(psi1-2phi_ZDCCA)>>
	
		if (trk1Chrg>0) {
		  // Calculate <<cos(2a-2Psi_V0)>> for POI Pos
		  //f2pCorrelatorCos2PsiDiff2PsiV0POIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C))); 
		  //f2pCorrelatorCos2PsiDiff2PsiV0POIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0A)));
		  // Calculate <<cos(2a-2Psi_ZDC)>> for POI Pos
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCC)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCAPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCA)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCAPOIPos->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCCA)));
		
		  fNITV0POIPos->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		  fNITV0POIPos->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		  fNITV0POIPos->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		  fNITV0POIPos->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		  fNITZDCPOIPos->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZNCC)); // <<cos(psi1-2phi_ZDCC)>>
		  fNITZDCPOIPos->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZNCC)); // <<sin(psi1-2phi_ZDCC)>>
		  fNITZDCPOIPos->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZNCA)); // <<cos(psi1-2phi_ZDCA)>>
		  fNITZDCPOIPos->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZNCA)); // <<sin(psi1-2phi_ZDCA)>>
		  fNITZDCPOIPos->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZNCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		  fNITZDCPOIPos->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZNCCA)); // <<sin(psi1-2phi_ZDCCA)>>
		}
	
		if (trk1Chrg<0) {
		  // Calculate <<cos(2a-2Psi_V0)>> for POI Neg
		  //f2pCorrelatorCos2PsiDiff2PsiV0POINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0C))); 
		  //f2pCorrelatorCos2PsiDiff2PsiV0POINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsi2V0A)));
		  // Calculate <<cos(2a-2Psi_ZDC)>> for POI Neg
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCC)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCAPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCA)));
		  //f2pCorrelatorCos2PsiDiff2PsiZDCCAPOINeg->Fill(0.5, TMath::Cos(2*(trk1Phi-fPsiZNCCA)));
		
		  fNITV0POINeg->Fill(0.5, TMath::Cos(trk1Phi-2*fPsi2V0C)); // <<cos(psi1-2phi_V0C)>>
		  fNITV0POINeg->Fill(1.5, TMath::Sin(trk1Phi-2*fPsi2V0C)); // <<sin(psi1-2phi_V0C)>>
		  fNITV0POINeg->Fill(2.5, TMath::Cos(trk1Phi-2*fPsi2V0A)); // <<cos(psi1-2phi_V0A)>>
		  fNITV0POINeg->Fill(3.5, TMath::Sin(trk1Phi-2*fPsi2V0A)); // <<sin(psi1-2phi_V0A)>>

		  fNITZDCPOINeg->Fill(0.5, TMath::Cos(trk1Phi-2*fPsiZNCC)); // <<cos(psi1-2phi_ZDCC)>>
		  fNITZDCPOINeg->Fill(1.5, TMath::Sin(trk1Phi-2*fPsiZNCC)); // <<sin(psi1-2phi_ZDCC)>>
		  fNITZDCPOINeg->Fill(2.5, TMath::Cos(trk1Phi-2*fPsiZNCA)); // <<cos(psi1-2phi_ZDCA)>>
		  fNITZDCPOINeg->Fill(3.5, TMath::Sin(trk1Phi-2*fPsiZNCA)); // <<sin(psi1-2phi_ZDCA)>>
		  fNITZDCPOINeg->Fill(4.5, TMath::Cos(trk1Phi-2*fPsiZNCCA)); // <<cos(psi1-2phi_ZDCCA)>>
		  fNITZDCPOINeg->Fill(5.5, TMath::Sin(trk1Phi-2*fPsiZNCCA)); // <<sin(psi1-2phi_ZDCCA)>>
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
		    
		    if (fDCAxyMax>0 && fDCAzzMax>0 && trk2DCAxy>=fDCAxyMax && trk2DCAz>=fDCAzzMax) { // if fDCAxyMax, fDCAzzMax is set to be less than 0, no cut applied
			  continue;
			}
			
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
				f2pCorrelatorCosPsiDiff->Fill(0.5, TMath::Cos(trk1Phi-trk2Phi)); 
				f2pCorrelatorCos2PsiDiff->Fill(0.5, TMath::Cos(2*(trk1Phi-trk2Phi))); 
				
				
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
  

  CalculateCMESPPP();

  // fill event information
  treeEvent->Fill();
  
  // Reset event-by-event variables
  ResetEventByEventQuantities();
  
  fDebugwEventCount->Fill(8.1); ///Left for Analysis
  
  fHistVertexZcm->Fill(pVtxZ);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

}//---------------- UserExec ----------------------

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::ResetEventByEventQuantities(){

  fCMEQReRP->Reset();
  fCMEQImRP->Reset();
  fCMEQRePOIPos->Reset();
  fCMEQImPOIPos->Reset();
  fCMEQRePOINeg->Reset();
  fCMEQImPOINeg->Reset();
  f2pCorrelatorCos2PsiDiff2PsiV0RP->Reset();
  f2pCorrelatorCos2PsiDiff2PsiZDCRP->Reset();
  fNITV0OS->Reset();
  fNITZDCOS->Reset();
  fNITV0POIPos->Reset();
  fNITZDCPOIPos->Reset();
  fNITV0POINeg->Reset();
  fNITZDCPOINeg->Reset();

  f2pCorrelatorCosPsiDiff->Reset();
  f2pCorrelatorCos2PsiDiff->Reset();
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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupAnalysisHistograms(){
  
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
  
  
  
}

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::CalculateCMESPPP()
{
  //************************************************ Get all variables ****************************
  //*********************************************** TPC part *************************************
  Int_t h = 0; //@Shi used for TPC and v0 part. For ZDCpart, it is set to 1
  Double_t e = 1E-5;
  
  // =========> begin same multiplicity #1 All number of tracks ===========
  Double_t uRPReTPC = fCMEQReRP->GetBinContent(1); // Qvec TPC cos(phi) 
  Double_t uRPImTPC = fCMEQImRP->GetBinContent(1); // Qvec TPC sin(phi) 
  Double_t uRPMultTPC = fCMEQReRP->GetBinEntries(1); 
  
  Double_t uRP2ReTPC = fCMEQReRP->GetBinContent(2); // Qvec TPC cos(2phi)
  Double_t uRP2ImTPC = fCMEQImRP->GetBinContent(2); // Qvec TPC sin(2phi)
  // end  same multiplicity #1 All number of tracks ===========
  
  // =========> begin same multiplicity #2 All Positive particles ===========
  Double_t uPOIPosReTPC = fCMEQRePOIPos->GetBinContent(1); // Cos(phi_POIPos)
  Double_t uPOIPosMult = fCMEQRePOIPos->GetBinEntries(1);
  Double_t uPOIPosImTPC = fCMEQImPOIPos->GetBinContent(1);
  // end same multiplicity #2 All Positive particles ===========

  // =========> begin same multiplicity #3 All Negative particles ===========
  Double_t uPOINegReTPC = fCMEQRePOINeg->GetBinContent(1); // Cos(phi_POINeg)
  Double_t uPOINegMult = fCMEQRePOINeg->GetBinEntries(1);
  Double_t uPOINegImTPC = fCMEQImPOINeg->GetBinContent(1);
  // end same multiplicity #3 All Negative particles ===========

  if(uPOIPosMult<0.1 || uPOINegMult<0.1) return;  //// This means there is not enough track in the the event. 
  
  // =========> begin same multiplicity #2 All Positive particles =========== (appear before no need to save)
  Double_t uPOIPos2ReTPC = fCMEQRePOIPos->GetBinContent(2); // Cos(2phi_POIPos)
  Double_t uPOIPos2ImTPC = fCMEQImPOIPos->GetBinContent(2);
  // end same multiplicity #2 All Positive particles ===========

  // =========> begin same multiplicity #3 All Negative particles =========== (appear before no need to save)
  Double_t uPOINeg2ReTPC = fCMEQRePOINeg->GetBinContent(2); // Cos(2phi_POINeg)
  Double_t uPOINeg2ImTPC = fCMEQImPOINeg->GetBinContent(2);
  // end same multiplicity #3 All Negative particles ===========

  // =========> begin same multiplicity #1 All number of tracks ===========
  Double_t u2pCorrelatorCos2PsiDiff2PsiV0CRP = f2pCorrelatorCos2PsiDiff2PsiV0RP->GetBinContent(1); //<cos(2psi1-2phi_V0C)>
  Double_t u2pCorrelatorCos2PsiDiff2PsiV0ARP = f2pCorrelatorCos2PsiDiff2PsiV0RP->GetBinContent(2); //<cos(2psi1-2phi_V0A)>

  Double_t u2pCorrelatorCos2PsiDiff2PsiZDCCRP = f2pCorrelatorCos2PsiDiff2PsiZDCRP->GetBinContent(1); //<cos(2psi1-2phi_ZDCC)>
  Double_t u2pCorrelatorCos2PsiDiff2PsiZDCARP = f2pCorrelatorCos2PsiDiff2PsiZDCRP->GetBinContent(2); //<cos(2psi1-2phi_ZDCA)>
  Double_t u2pCorrelatorCos2PsiDiff2PsiZDCCARP = f2pCorrelatorCos2PsiDiff2PsiZDCRP->GetBinContent(3); //<cos(2psi1-2phi_ZDCCA)>

	    
  Double_t uNITCosPsidiff2PsiV0COS = fNITV0OS->GetBinContent(1); // <<cos(psi1-2phi_V0C)>> 
  Double_t uNITSinPsidiff2PsiV0COS = fNITV0OS->GetBinContent(2); // <<sin(psi1-2phi_V0C)>>
  Double_t uNITCosPsidiff2PsiV0AOS = fNITV0OS->GetBinContent(3); // <<cos(psi1-2phi_V0A)>>
  Double_t uNITSinPsidiff2PsiV0AOS = fNITV0OS->GetBinContent(4); // <<sin(psi1-2phi_V0A)>>

  Double_t uNITCosPsidiff2PsiZDCCOS = fNITZDCOS->GetBinContent(1); // <<cos(psi1-2phi_ZDCC)>>
  Double_t uNITSinPsidiff2PsiZDCCOS = fNITZDCOS->GetBinContent(2); // <<sin(psi1-2phi_ZDCC)>>
  Double_t uNITCosPsidiff2PsiZDCAOS = fNITZDCOS->GetBinContent(3); // <<cos(psi1-2phi_ZDCA)>>
  Double_t uNITSinPsidiff2PsiZDCAOS = fNITZDCOS->GetBinContent(4); // <<sin(psi1-2phi_ZDCA)>>
  Double_t uNITCosPsidiff2PsiZDCCAOS = fNITZDCOS->GetBinContent(5); // <<cos(psi1-2phi_ZDCCA)>>
  Double_t uNITSinPsidiff2PsiZDCCAOS = fNITZDCOS->GetBinContent(6); // <<sin(psi1-2phi_ZDCCA)>>
  // end  same multiplicity #1 All number of tracks ===========


  // =========> begin same multiplicity #2 All Positive particles ===========
  Double_t uNITCosPsidiff2PsiV0CPOIPos = fNITV0POIPos->GetBinContent(1); // <<cos(psi1-2phi_V0C)>> 
  Double_t uNITSinPsidiff2PsiV0CPOIPos = fNITV0POIPos->GetBinContent(2); // <<sin(psi1-2phi_V0C)>>
  Double_t uNITCosPsidiff2PsiV0APOIPos = fNITV0POIPos->GetBinContent(3); // <<cos(psi1-2phi_V0A)>>
  Double_t uNITSinPsidiff2PsiV0APOIPos = fNITV0POIPos->GetBinContent(4); // <<sin(psi1-2phi_V0A)>>

  Double_t uNITCosPsidiff2PsiZDCCPOIPos = fNITZDCPOIPos->GetBinContent(1); // <<cos(psi1-2phi_ZDCC)>>
  Double_t uNITSinPsidiff2PsiZDCCPOIPos = fNITZDCPOIPos->GetBinContent(2); // <<sin(psi1-2phi_ZDCC)>>
  Double_t uNITCosPsidiff2PsiZDCAPOIPos = fNITZDCPOIPos->GetBinContent(3); // <<cos(psi1-2phi_ZDCA)>>
  Double_t uNITSinPsidiff2PsiZDCAPOIPos = fNITZDCPOIPos->GetBinContent(4); // <<sin(psi1-2phi_ZDCA)>>
  Double_t uNITCosPsidiff2PsiZDCCAPOIPos = fNITZDCPOIPos->GetBinContent(5); // <<cos(psi1-2phi_ZDCCA)>>
  Double_t uNITSinPsidiff2PsiZDCCAPOIPos = fNITZDCPOIPos->GetBinContent(6); // <<sin(psi1-2phi_ZDCCA)>>
  // end same multiplicity #2 All Positive particles ===========

  // =========> begin same multiplicity #3 All Negative particles ===========
  Double_t uNITCosPsidiff2PsiV0CPOINeg = fNITV0POINeg->GetBinContent(1); // <<cos(psi1-2phi_V0C)>> 
  Double_t uNITSinPsidiff2PsiV0CPOINeg = fNITV0POINeg->GetBinContent(2); // <<sin(psi1-2phi_V0C)>>
  Double_t uNITCosPsidiff2PsiV0APOINeg = fNITV0POINeg->GetBinContent(3); // <<cos(psi1-2phi_V0A)>>
  Double_t uNITSinPsidiff2PsiV0APOINeg = fNITV0POINeg->GetBinContent(4); // <<sin(psi1-2phi_V0A)>>

  Double_t uNITCosPsidiff2PsiZDCCPOINeg = fNITZDCPOINeg->GetBinContent(1); // <<cos(psi1-2phi_ZDCC)>>
  Double_t uNITSinPsidiff2PsiZDCCPOINeg = fNITZDCPOINeg->GetBinContent(2); // <<sin(psi1-2phi_ZDCC)>>
  Double_t uNITCosPsidiff2PsiZDCAPOINeg = fNITZDCPOINeg->GetBinContent(3); // <<cos(psi1-2phi_ZDCA)>>
  Double_t uNITSinPsidiff2PsiZDCAPOINeg = fNITZDCPOINeg->GetBinContent(4); // <<sin(psi1-2phi_ZDCA)>>
  Double_t uNITCosPsidiff2PsiZDCCAPOINeg = fNITZDCPOINeg->GetBinContent(5); // <<cos(psi1-2phi_ZDCCA)>>
  Double_t uNITSinPsidiff2PsiZDCCAPOINeg = fNITZDCPOINeg->GetBinContent(6); // <<sin(psi1-2phi_ZDCCA)>>
  // end same multiplicity #3 All Negative particles ===========


  // =========> begin same multiplicity #4 All RP pairs ===========
  Double_t u2pCorrelatorCosPsiDiff = f2pCorrelatorCosPsiDiff->GetBinContent(1); // <cos(dPsi1-dPsi2)>
  Double_t u2pCorrelatorCos2PsiDiff = f2pCorrelatorCos2PsiDiff->GetBinContent(1); // <cos(2(dPsi1-dPsi2))>
  Double_t u2pCorrelatorRPMult = f2pCorrelatorCosPsiDiff->GetBinEntries(1);
  // end same multiplicity #4 All RP pairs ===========

  // =========> begin same multiplicity #5 All POI OS pairs ===========
  Double_t u2pCorrelatorCosPsiSumPOIOS = fRePEBEOS->GetBinContent(1); // <cos(dPsi1+dPsi2)>
  Double_t u2pCorrelatorPOIOSMult = fRePEBEOS->GetBinEntries(1);
  Double_t u2pCorrelatorSinPsiSumPOIOS = fImPEBEOS->GetBinContent(1); // <sin(dPsi1+dPsi2)>
  Double_t u2pCorrelatorCosPsiDiffPOIOS = f2pCorrelatorCosPsiDiffOS->GetBinContent(1); // <cos(dPsi1-dPsi2)>
  Double_t u2pCorrelatorCos2PsiDiffPOIOS = f2pCorrelatorCos2PsiDiffOS->GetBinContent(1); // <cos(2(dPsi1-dPsi2))>
  // end same multiplicity #5 All POI OS pairs ===========

  // =========> begin same multiplicity #6 All POI PP pairs ===========
  Double_t u2pCorrelatorCosPsiSumPOIPP = fRePEBEPP->GetBinContent(1); // <cos(dPsi1+dPsi2)>
  Double_t u2pCorrelatorPOIPPMult = fRePEBEPP->GetBinEntries(1);
  Double_t u2pCorrelatorSinPsiSumPOIPP = fImPEBEPP->GetBinContent(1); // <sin(dPsi1+dPsi2)>
  Double_t u2pCorrelatorCosPsiDiffPOIPP = f2pCorrelatorCosPsiDiffPP->GetBinContent(1); // <cos(dPsi1-dPsi2)>
  Double_t u2pCorrelatorCos2PsiDiffPOIPP = f2pCorrelatorCos2PsiDiffPP->GetBinContent(1); // <cos(2(dPsi1-dPsi2))>
  // end same multiplicity #6 All POI PP pairs ===========

  // =========> begin same multiplicity #7 All POI NN pairs ===========
  Double_t u2pCorrelatorCosPsiSumPOINN = fRePEBENN->GetBinContent(1); // <cos(dPsi1+dPsi2)>
  Double_t u2pCorrelatorPOINNMult = fRePEBENN->GetBinEntries(1);
  Double_t u2pCorrelatorSinPsiSumPOINN = fImPEBENN->GetBinContent(1); // <sin(dPsi1+dPsi2)>
  Double_t u2pCorrelatorCosPsiDiffPOINN = f2pCorrelatorCosPsiDiffNN->GetBinContent(1); // <cos(dPsi1-dPsi2)>
  Double_t u2pCorrelatorCos2PsiDiffPOINN = f2pCorrelatorCos2PsiDiffNN->GetBinContent(1); // <cos(2(dPsi1-dPsi2))>
  
  // end same multiplicity #7 All POI NN pairs ===========

    
  // set fpQvecEvent
  fpQvecEvent->setRPReTPC( uRPReTPC ); // w * cos(theta+) eta+
  fpQvecEvent->setRPImTPC( uRPImTPC );
  fpQvecEvent->setRPMultTPC( uRPMultTPC );
  
  fpQvecEvent->setRP2ReTPC( uRP2ReTPC );
  fpQvecEvent->setRP2ImTPC( uRP2ImTPC );
  
  fpQvecEvent->setPOIPosReTPC( uPOIPosReTPC );
  fpQvecEvent->setPOIPosMult( uPOIPosMult );
  fpQvecEvent->setPOIPosImTPC( uPOIPosImTPC );
  
  fpQvecEvent->setPOINegReTPC( uPOINegReTPC );
  fpQvecEvent->setPOINegMult( uPOINegMult );
  fpQvecEvent->setPOINegImTPC( uPOINegImTPC );
  
  fpQvecEvent->setPOIPos2ReTPC( uPOIPos2ReTPC );
  fpQvecEvent->setPOIPos2ImTPC( uPOIPos2ImTPC );
  
  fpQvecEvent->setPOINeg2ReTPC( uPOINeg2ReTPC );
  fpQvecEvent->setPOINeg2ImTPC( uPOINeg2ImTPC );

  fpQvecEvent->set2pCorrelatorCos2PsiDiff2PsiV0CRP( u2pCorrelatorCos2PsiDiff2PsiV0CRP );
  fpQvecEvent->set2pCorrelatorCos2PsiDiff2PsiV0ARP( u2pCorrelatorCos2PsiDiff2PsiV0ARP );
  
  fpQvecEvent->set2pCorrelatorCos2PsiDiff2PsiZDCCRP( u2pCorrelatorCos2PsiDiff2PsiZDCCRP );
  fpQvecEvent->set2pCorrelatorCos2PsiDiff2PsiZDCARP( u2pCorrelatorCos2PsiDiff2PsiZDCARP );
  fpQvecEvent->set2pCorrelatorCos2PsiDiff2PsiZDCCARP( u2pCorrelatorCos2PsiDiff2PsiZDCCARP );

  fpQvecEvent->setNITCosPsidiff2PsiV0COS( uNITCosPsidiff2PsiV0COS );
  fpQvecEvent->setNITSinPsidiff2PsiV0COS( uNITSinPsidiff2PsiV0COS );
  fpQvecEvent->setNITCosPsidiff2PsiV0AOS( uNITCosPsidiff2PsiV0AOS );
  fpQvecEvent->setNITSinPsidiff2PsiV0AOS( uNITSinPsidiff2PsiV0AOS );

  fpQvecEvent->setNITCosPsidiff2PsiZDCCOS( uNITCosPsidiff2PsiZDCCOS );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCOS( uNITSinPsidiff2PsiZDCCOS );
  fpQvecEvent->setNITCosPsidiff2PsiZDCAOS( uNITCosPsidiff2PsiZDCAOS );
  fpQvecEvent->setNITSinPsidiff2PsiZDCAOS( uNITSinPsidiff2PsiZDCAOS );
  fpQvecEvent->setNITCosPsidiff2PsiZDCCAOS( uNITCosPsidiff2PsiZDCCAOS );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCAOS( uNITSinPsidiff2PsiZDCCAOS );

  fpQvecEvent->setNITCosPsidiff2PsiV0CPOIPos( uNITCosPsidiff2PsiV0CPOIPos );
  fpQvecEvent->setNITSinPsidiff2PsiV0CPOIPos( uNITSinPsidiff2PsiV0CPOIPos );
  fpQvecEvent->setNITCosPsidiff2PsiV0APOIPos( uNITCosPsidiff2PsiV0APOIPos );
  fpQvecEvent->setNITSinPsidiff2PsiV0APOIPos( uNITSinPsidiff2PsiV0APOIPos );

  fpQvecEvent->setNITCosPsidiff2PsiZDCCPOIPos( uNITCosPsidiff2PsiZDCCPOIPos );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCPOIPos( uNITSinPsidiff2PsiZDCCPOIPos );
  fpQvecEvent->setNITCosPsidiff2PsiZDCAPOIPos( uNITCosPsidiff2PsiZDCAPOIPos );
  fpQvecEvent->setNITSinPsidiff2PsiZDCAPOIPos( uNITSinPsidiff2PsiZDCAPOIPos );
  fpQvecEvent->setNITCosPsidiff2PsiZDCCAPOIPos( uNITCosPsidiff2PsiZDCCAPOIPos );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCAPOIPos( uNITSinPsidiff2PsiZDCCAPOIPos );

  fpQvecEvent->setNITCosPsidiff2PsiV0CPOINeg( uNITCosPsidiff2PsiV0CPOINeg );
  fpQvecEvent->setNITSinPsidiff2PsiV0CPOINeg( uNITSinPsidiff2PsiV0CPOINeg );
  fpQvecEvent->setNITCosPsidiff2PsiV0APOINeg( uNITCosPsidiff2PsiV0APOINeg );
  fpQvecEvent->setNITSinPsidiff2PsiV0APOINeg( uNITSinPsidiff2PsiV0APOINeg );

  fpQvecEvent->setNITCosPsidiff2PsiZDCCPOINeg( uNITCosPsidiff2PsiZDCCPOINeg );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCPOINeg( uNITSinPsidiff2PsiZDCCPOINeg );
  fpQvecEvent->setNITCosPsidiff2PsiZDCAPOINeg( uNITCosPsidiff2PsiZDCAPOINeg );
  fpQvecEvent->setNITSinPsidiff2PsiZDCAPOINeg( uNITSinPsidiff2PsiZDCAPOINeg );
  fpQvecEvent->setNITCosPsidiff2PsiZDCCAPOINeg( uNITCosPsidiff2PsiZDCCAPOINeg );
  fpQvecEvent->setNITSinPsidiff2PsiZDCCAPOINeg( uNITSinPsidiff2PsiZDCCAPOINeg );
  
  fpQvecEvent->set2pCorrelatorCosPsiDiff( u2pCorrelatorCosPsiDiff );
  fpQvecEvent->set2pCorrelatorCos2PsiDiff( u2pCorrelatorCos2PsiDiff );
  fpQvecEvent->set2pCorrelatorRPMult( u2pCorrelatorRPMult );

  fpQvecEvent->set2pCorrelatorCosPsiSumPOIOS( u2pCorrelatorCosPsiSumPOIOS );
  fpQvecEvent->set2pCorrelatorPOIOSMult( u2pCorrelatorPOIOSMult );
  fpQvecEvent->set2pCorrelatorSinPsiSumPOIOS( u2pCorrelatorSinPsiSumPOIOS );
  fpQvecEvent->set2pCorrelatorCosPsiDiffPOIOS( u2pCorrelatorCosPsiDiffPOIOS );
  fpQvecEvent->set2pCorrelatorCos2PsiDiffPOIOS( u2pCorrelatorCos2PsiDiffPOIOS );

  fpQvecEvent->set2pCorrelatorCosPsiSumPOIPP( u2pCorrelatorCosPsiSumPOIPP );
  fpQvecEvent->set2pCorrelatorPOIPPMult( u2pCorrelatorPOIPPMult );
  fpQvecEvent->set2pCorrelatorSinPsiSumPOIPP( u2pCorrelatorSinPsiSumPOIPP );
  fpQvecEvent->set2pCorrelatorCosPsiDiffPOIPP( u2pCorrelatorCosPsiDiffPOIPP );
  fpQvecEvent->set2pCorrelatorCos2PsiDiffPOIPP( u2pCorrelatorCos2PsiDiffPOIPP );
  
  fpQvecEvent->set2pCorrelatorCosPsiSumPOINN( u2pCorrelatorCosPsiSumPOINN );
  fpQvecEvent->set2pCorrelatorPOINNMult( u2pCorrelatorPOINNMult );
  fpQvecEvent->set2pCorrelatorSinPsiSumPOINN( u2pCorrelatorSinPsiSumPOINN );
  fpQvecEvent->set2pCorrelatorCosPsiDiffPOINN( u2pCorrelatorCosPsiDiffPOINN );
  fpQvecEvent->set2pCorrelatorCos2PsiDiffPOINN( u2pCorrelatorCos2PsiDiffPOINN );
}


void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupQAHistograms(){

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
  
  fCMEQRePOIPos = new TProfile("fCMEQRePOIPos", "fCMEQRePOIPos",2,0,2);
  fTempList->Add(fCMEQRePOIPos);
  
  fCMEQImPOIPos = new TProfile("fCMEQImPOIPos", "fCMEQImPOIPos",2,0,2);
  fTempList->Add(fCMEQImPOIPos);
  
  fCMEQRePOINeg = new TProfile("fCMEQRePOINeg", "fCMEQRePOINeg",2,0,2);
  fTempList->Add(fCMEQRePOINeg);
  
  fCMEQImPOINeg = new TProfile("fCMEQImPOINeg", "fCMEQImPOINeg",2,0,2);
  fTempList->Add(fCMEQImPOINeg);
  
  f2pCorrelatorCos2PsiDiff2PsiV0RP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiV0RP", "f2pCorrelatorCos2PsiDiff2PsiV0RP",2,0,2);
  fTempList->Add(f2pCorrelatorCos2PsiDiff2PsiV0RP);
  
  f2pCorrelatorCos2PsiDiff2PsiZDCRP = new TProfile("f2pCorrelatorCos2PsiDiff2PsiZDCRP", "f2pCorrelatorCos2PsiDiff2PsiZDCRP",3,0,3);
  fTempList->Add(f2pCorrelatorCos2PsiDiff2PsiZDCRP);

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

  
  f2pCorrelatorCosPsiDiff = new TProfile("f2pCorrelatorCosPsiDiff", "f2pCorrelatorCosPsiDiff",1,0,1);
  fTempList->Add(f2pCorrelatorCosPsiDiff);
  
  f2pCorrelatorCos2PsiDiff = new TProfile("f2pCorrelatorCos2PsiDiff", "f2pCorrelatorCos2PsiDiff",1,0,1);
  fTempList->Add(f2pCorrelatorCos2PsiDiff);
  
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





void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupPileUpRemovalFunctions18qPass3() { //@Shi for 2018 period Pass3 data
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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupPileUpRemovalFunctions18rPass3() { //@Shi for 2018 period Pass3 data
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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupEventAndTaskConfigInfo(){

 
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


void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupQvecSavingObjects() {
  
  fpQvecEvent = new AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple();
  
  fEventList = new TList();
  fEventList->SetName("Event List");
  fEventList->SetOwner(kTRUE);
  fListHist->Add(fEventList);
  
  treeEvent = new TTree("events", "event");
  treeEvent->Branch("event", &fpQvecEvent);
  
  fEventList->Add(treeEvent);
	
}

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetNUACorrectionHist(Int_t run, Int_t kParticleID) {
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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetMCCorrectionHist(){

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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetV0MCorrectionHist(Int_t run){ 

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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetZDCCorrectionHist(Int_t run){ 
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

Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){
  
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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::CheckEventIsPileUp2018(AliAODEvent *faod) {

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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetTPCQvectAndRemovePileUp2018(AliAODEvent *faod,Double_t *qnxEtaNeg,Double_t *qnyEtaNeg,Double_t *qnxEtaPos,Double_t *qnyEtaPos,Double_t& multNeg,Double_t& multPos)
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


Double_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetNUAWeightForTrack(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

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

Double_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetNUAWeightForTrackPID(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg){

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


Double_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetNUAWeightForTrackPID(Double_t fVtxZ,Double_t fPhi,Double_t fEta,Int_t gChrg,Int_t gPID){

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


void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::ApplyTPCqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t& qxEtaNeg, Double_t& qyEtaNeg,Double_t& qxEtaPos,Double_t& qyEtaPos){

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

void  AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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





Double_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetMCEfficiencyWeightForTrack(Double_t fPt,Int_t gChrg,Int_t kPID){
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
	


Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

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



Bool_t AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A, Double_t &sumMultV0C, Double_t &sumMultV0A){

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

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::Terminate(Option_t *)  {
  //fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  //if (!fOutputList) return;
}





