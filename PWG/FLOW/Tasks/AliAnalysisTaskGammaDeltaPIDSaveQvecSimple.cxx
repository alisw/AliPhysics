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
  }
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::AliAnalysisTaskGammaDeltaPIDSaveQvecSimple():
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
  }
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
  Double_t wgtComb1Ch = 1.0;
  Double_t ptWgtMCChtrk1  = 1.0, WgtNUAChtrk1   = 1.0;

  Double_t localSumQ2x =0,localSumQ2y=0;
  Double_t localSumQ3x =0,localSumQ3y=0;
  Double_t localMultTPC=0; 

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
		ptWgtMCChtrk1 = 1.0;  


		WgtNUAChtrk1  = GetNUAWeightForTrack(fVertexZEvent,trk1Phi,trk1Eta,trk1Chrg);          
		ptWgtMCChtrk1 = GetMCEfficiencyWeightForTrack(trk1Pt,trk1Chrg,0);

		wgtComb1Ch  = WgtNUAChtrk1*ptWgtMCChtrk1;    /// Charge

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
		}
	

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
  }
  
  fDebugwEventCount->Fill(9.1); ///Left for Analysis
  
  fHistVertexZcm->Fill(pVtxZ);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

}//---------------- UserExec ----------------------



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
  }
  
}

void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::CalculateCMESPPP()
{
  //************************************************ Get all variables ****************************
  //*********************************************** TPC part *************************************
  Int_t h = 0; //@Shi used for TPC and v0 part. For ZDCpart, it is set to 1
  Double_t e = 1E-5;
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
  
  Double_t uPReTPCSubPosEta = 0, uPImTPCSubPosEta = 0, uP2ReTPCSubPosEta = 0, uP2ImTPCSubPosEta = 0, uP2Re2TPCSubPosEta = 0, uP2Im2TPCSubPosEta = 0, uPMTPCSubPosEta = 0, uP2MTPCSubPosEta = 0, uNReTPCSubPosEta = 0, uNImTPCSubPosEta = 0, uN2ReTPCSubPosEta = 0, uN2ImTPCSubPosEta = 0, uN2Re2TPCSubPosEta = 0, uN2Im2TPCSubPosEta = 0, uNMTPCSubPosEta = 0, uN2MTPCSubPosEta = 0;
  
  Double_t uPReTPCSubNegEta = 0, uPImTPCSubNegEta = 0, uP2ReTPCSubNegEta = 0, uP2ImTPCSubNegEta = 0, uP2Re2TPCSubNegEta = 0, uP2Im2TPCSubNegEta = 0, uPMTPCSubNegEta = 0, uP2MTPCSubNegEta = 0, uNReTPCSubNegEta = 0, uNImTPCSubNegEta = 0, uN2ReTPCSubNegEta = 0, uN2ImTPCSubNegEta = 0, uN2Re2TPCSubNegEta = 0, uN2Im2TPCSubNegEta = 0, uNMTPCSubNegEta = 0, uN2MTPCSubNegEta = 0;
  
  
  
  for(Int_t EBin=1; EBin<=fCMEQRe[0][0]->GetNbinsX(); EBin++) {
    if (EBin == fCMEQRe[0][h]->FindBin(0.4)) { // positive eta region
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
	} else if (EBin == fCMEQRe[0][h]->FindBin(-0.4)) { // negative eta region
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
	} else if (EBin == fCMEQRe[0][h]->FindBin(0.05)) { // sub pos eta region
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
	} else if (EBin == fCMEQRe[0][h]->FindBin(-0.05)) { // negative eta region
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







void AliAnalysisTaskGammaDeltaPIDSaveQvecSimple::SetupPileUpRemovalFunctions(){
  
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





