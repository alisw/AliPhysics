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
/* $Id: AliAnalysisTaskGammaDeltaPID.cxx  by Md Rihan Haque           $   */
/* Simple Task to fill V0 and ZDC Energies for Gain Calibration           */
/* Works with 15o and 18q/r. Support to be added for LHC10h               */
/* Developer: Md Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)        */
/* Last Modified: Aug 23, 2021,  first version committed                  */
////////////////////////////////////////////////////////////////////////////


//-- general include---
#include "TChain.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"
#include "TParticle.h"
#include "TMath.h"
#include "stdio.h"
#include "Riostream.h"

//---- manager and handler---
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//---V0 and ZDC info---
#include "AliAODZDC.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"
//----V0 particles:
#include "AliAODv0.h"

//---AOD,ESD event--
#include "AliESDEvent.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliTimeRangeCut.h"

//----- For MC event------
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

//----for PID-----
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

//----- Vevent and tracks
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliCentrality.h"


//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaDeltaPID.h"

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskGammaDeltaPID)

AliAnalysisTaskGammaDeltaPID::AliAnalysisTaskGammaDeltaPID(const char *name): AliAnalysisTaskSE(name),
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
  fCentralityMin(0),
  fCentralityMax(90),
  gHarmonic(2),
  gParticleID(1),
  fFilterBit(1),
  fTPCclustMin(70),
  gOldRunNumber(11),
  bUseKinkTracks(kFALSE),
  fNSigmaTPCCut(2.0),
  fNSigmaTOFCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fChi2Max(4.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bSkipPileUpCut(kFALSE), 
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  sDetectorforEP("kNone"),
  bUseV0EventPlane(kFALSE),  
  bSkipAnalysis(kFALSE),
  
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
  fMassMean(1.115683),
  fLambdaMassCut(0.005),
  bAnalysLambdaPairs(kFALSE),

  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistCL0VsV0MBefore(NULL),  
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistCL0VsV0MAfter(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistPileUpCount(NULL),
  fHistAnalysisInfo(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),   
  
  fAvgCosNPsivsCentEtaPos(NULL),
  fAvgSinNPsivsCentEtaPos(NULL),
  fAvgCosNPsivsCentEtaNeg(NULL),
  fAvgSinNPsivsCentEtaNeg(NULL),
  fAvgCos3PsivsCentEtaPos(NULL),
  fAvgSin3PsivsCentEtaPos(NULL),
  fAvgCos3PsivsCentEtaNeg(NULL),
  fAvgSin3PsivsCentEtaNeg(NULL),  
  fHistV0CDetqVectorvsCent(NULL),
  fHistV0ADetqVectorvsCent(NULL),
  fHistTPCPosqVectorvsCent(NULL),
  fHistTPCNegqVectorvsCent(NULL),
  
  fHistVertexZcm(NULL),
  fHistVxvsVzMinBias(NULL),
  fHistVyvsVzMinBias(NULL),  
  hAvgZNACh0vsCentVz(NULL),
  hAvgZNCCh0vsCentVz(NULL),
  hAvgZNACh1vsCentVz(NULL),
  hAvgZNCCh1vsCentVz(NULL),
  hAvgZNACh2vsCentVz(NULL),
  hAvgZNCCh2vsCentVz(NULL),
  hAvgZNACh3vsCentVz(NULL),
  hAvgZNCCh3vsCentVz(NULL),
  hAvgZNACh4vsCentVz(NULL),
  hAvgZNCCh4vsCentVz(NULL),
  hAvgV0ChannelsvsVz(NULL),
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQ3XvsCentV0C(NULL),
  hAvgQ3YvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hAvgQ3XvsCentV0A(NULL),
  hAvgQ3YvsCentV0A(NULL),  

  fHistTPCPsiNPosPlane(NULL),
  fHistTPCPsiNNegPlane(NULL),
  fHistTPCPsi3PosPlane(NULL),
  fHistTPCPsi3NegPlane(NULL),
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  fHistV0CPsi3EventPlane(NULL),
  fHistV0APsi3EventPlane(NULL),  
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hV0CV0APsiNCorrelation(NULL),
  hV0CTPCPsiNCorrelation(NULL),
  hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  hV0CTPCPsi3Correlation(NULL),
  hV0ATPCPsi3Correlation(NULL),  
  
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

  fHCorrectV0ChWeghts(NULL),
  fHCorrectNUAChrgPos(NULL),
  fHCorrectNUAChrgNeg(NULL),
  fHCorrectNUAkPIDPos(NULL),
  fHCorrectNUAkPIDNeg(NULL),
  fHCorrectMCPosChrg(NULL),
  fHCorrectMCPosPion(NULL),
  fHCorrectMCPosKaon(NULL),
  fHCorrectMCPosProt(NULL),
  fHCorrectMCNegChrg(NULL),
  fHCorrectMCNegPion(NULL),
  fHCorrectMCNegKaon(NULL),
  fHCorrectMCNegProt(NULL),  
  fHCorrectTPCQNxEtaPos(NULL),
  fHCorrectTPCQNyEtaPos(NULL),
  fHCorrectTPCQNxEtaNeg(NULL),
  fHCorrectTPCQNyEtaNeg(NULL),
  fHCorrectTPCQ3xEtaPos(NULL),
  fHCorrectTPCQ3yEtaPos(NULL),
  fHCorrectTPCQ3xEtaNeg(NULL),
  fHCorrectTPCQ3yEtaNeg(NULL),    
  fHCorrectQNxV0C(NULL),
  fHCorrectQNyV0C(NULL),    
  fHCorrectQNxV0A(NULL),
  fHCorrectQNyV0A(NULL),  
  fHCorrectQ3xV0C(NULL),
  fHCorrectQ3yV0C(NULL),    
  fHCorrectQ3xV0A(NULL),
  fHCorrectQ3yV0A(NULL),
  
  fHistV0Pt(NULL),              
  fHistV0Eta(NULL),              
  fHistV0DcatoPrimVertex(NULL), 
  fHistV0CPA(NULL),             
  fHistV0DecayLength(NULL),     

  fProfileGammaTPC_Lambda_hPos(NULL),
  fProfileGammaTPC_Lambda_hNeg(NULL),
  fProfileGammaTPC_Lambda_Proton(NULL),
  fProfileGammaTPC_Lambda_AntiProton(NULL),
  fProfileGammaTPC_AntiLambda_hPos(NULL),
  fProfileGammaTPC_AntiLambda_hNeg(NULL),
  fProfileGammaTPC_AntiLambda_Proton(NULL),
  fProfileGammaTPC_AntiLambda_AntiProton(NULL)
{

  for(int i=0;i<3;i++){
    fCurrentVtx[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    fHistLambdaPt[i]  = NULL;            
    fHistLambdaEta[i] = NULL;
    fHistLambdaDcaToPrimVertex[i] = NULL;
    fHistLambdaCPA[i] = NULL;
    fHistLambdaDecayLength[i] = NULL;
    fHistLambdaMass[i] = NULL;
    fProfileLambdaMassVsPt[i] = NULL;
    fHistAntiLambdaPt[i] = NULL;
    fHistAntiLambdaEta[i] = NULL;
    fHistAntiLambdaDcaToPrimVertex[i] = NULL;
    fHistAntiLambdaCPA[i] = NULL;
    fHistAntiLambdaDecayLength[i] = NULL;
    fHistAntiLambdaMass[i] = NULL;
    fProfileAntiLambdaMassVsPt[i] = NULL;
  }  
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskGammaDeltaPID::AliAnalysisTaskGammaDeltaPID():
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
  fCentralityMin(0),
  fCentralityMax(90),
  gHarmonic(2),
  gParticleID(1),
  fFilterBit(1),
  fTPCclustMin(70),
  gOldRunNumber(11),
  bUseKinkTracks(kFALSE),
  fNSigmaTPCCut(2.0),
  fNSigmaTOFCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fChi2Max(4.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bSkipPileUpCut(kFALSE), 
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  sDetectorforEP("kNone"),
  bUseV0EventPlane(kFALSE),   
  bSkipAnalysis(kFALSE),
  
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
  fMassMean(1.115683),
  fLambdaMassCut(0.005),
  bAnalysLambdaPairs(kFALSE),
  
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),  
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistCL0VsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistCL0VsV0MAfter(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistPileUpCount(NULL),
  fHistAnalysisInfo(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  
  fAvgCosNPsivsCentEtaPos(NULL),
  fAvgSinNPsivsCentEtaPos(NULL),
  fAvgCosNPsivsCentEtaNeg(NULL),
  fAvgSinNPsivsCentEtaNeg(NULL),
  fAvgCos3PsivsCentEtaPos(NULL),
  fAvgSin3PsivsCentEtaPos(NULL),
  fAvgCos3PsivsCentEtaNeg(NULL),
  fAvgSin3PsivsCentEtaNeg(NULL),  
  fHistV0CDetqVectorvsCent(NULL),
  fHistV0ADetqVectorvsCent(NULL),
  fHistTPCPosqVectorvsCent(NULL),
  fHistTPCNegqVectorvsCent(NULL),
  
  fHistVertexZcm(NULL),
  fHistVxvsVzMinBias(NULL),
  fHistVyvsVzMinBias(NULL),  
  hAvgZNACh0vsCentVz(NULL),
  hAvgZNCCh0vsCentVz(NULL),
  hAvgZNACh1vsCentVz(NULL),
  hAvgZNCCh1vsCentVz(NULL),
  hAvgZNACh2vsCentVz(NULL),
  hAvgZNCCh2vsCentVz(NULL),
  hAvgZNACh3vsCentVz(NULL),
  hAvgZNCCh3vsCentVz(NULL),
  hAvgZNACh4vsCentVz(NULL),
  hAvgZNCCh4vsCentVz(NULL),
  hAvgV0ChannelsvsVz(NULL),
  hAvgQNXvsCentV0C(NULL),
  hAvgQNYvsCentV0C(NULL),
  hAvgQ3XvsCentV0C(NULL),
  hAvgQ3YvsCentV0C(NULL),
  hAvgQNXvsCentV0A(NULL),
  hAvgQNYvsCentV0A(NULL),
  hAvgQ3XvsCentV0A(NULL),
  hAvgQ3YvsCentV0A(NULL),

  fHistTPCPsiNPosPlane(NULL),
  fHistTPCPsiNNegPlane(NULL),
  fHistTPCPsi3PosPlane(NULL),
  fHistTPCPsi3NegPlane(NULL),  
  fHistV0CPsiNEventPlane(NULL),
  fHistV0APsiNEventPlane(NULL),
  fHistV0CPsi3EventPlane(NULL),
  fHistV0APsi3EventPlane(NULL),
  hTPCPsiNCorrelation(NULL),
  hTPCPsi3Correlation(NULL),
  hV0CV0APsiNCorrelation(NULL),
  hV0CTPCPsiNCorrelation(NULL),
  hV0ATPCPsiNCorrelation(NULL),
  hV0CV0APsi3Correlation(NULL),
  hV0CTPCPsi3Correlation(NULL),
  hV0ATPCPsi3Correlation(NULL), 
  
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

  fHCorrectV0ChWeghts(NULL),
  fHCorrectNUAChrgPos(NULL),
  fHCorrectNUAChrgNeg(NULL),
  fHCorrectNUAkPIDPos(NULL),
  fHCorrectNUAkPIDNeg(NULL),
  fHCorrectMCPosChrg(NULL),
  fHCorrectMCPosPion(NULL),
  fHCorrectMCPosKaon(NULL),
  fHCorrectMCPosProt(NULL),
  fHCorrectMCNegChrg(NULL),
  fHCorrectMCNegPion(NULL),
  fHCorrectMCNegKaon(NULL),
  fHCorrectMCNegProt(NULL),
  fHCorrectTPCQNxEtaPos(NULL),
  fHCorrectTPCQNyEtaPos(NULL),
  fHCorrectTPCQNxEtaNeg(NULL),
  fHCorrectTPCQNyEtaNeg(NULL),
  fHCorrectTPCQ3xEtaPos(NULL),
  fHCorrectTPCQ3yEtaPos(NULL),
  fHCorrectTPCQ3xEtaNeg(NULL),
  fHCorrectTPCQ3yEtaNeg(NULL),  
  fHCorrectQNxV0C(NULL),
  fHCorrectQNyV0C(NULL),    
  fHCorrectQNxV0A(NULL),
  fHCorrectQNyV0A(NULL),  
  fHCorrectQ3xV0C(NULL),
  fHCorrectQ3yV0C(NULL),    
  fHCorrectQ3xV0A(NULL),
  fHCorrectQ3yV0A(NULL),

  fHistV0Pt(NULL),              
  fHistV0Eta(NULL),              
  fHistV0DcatoPrimVertex(NULL), 
  fHistV0CPA(NULL),             
  fHistV0DecayLength(NULL), 

  fProfileGammaTPC_Lambda_hPos(NULL),
  fProfileGammaTPC_Lambda_hNeg(NULL),
  fProfileGammaTPC_Lambda_Proton(NULL),
  fProfileGammaTPC_Lambda_AntiProton(NULL),
  fProfileGammaTPC_AntiLambda_hPos(NULL),
  fProfileGammaTPC_AntiLambda_hNeg(NULL),
  fProfileGammaTPC_AntiLambda_Proton(NULL),
  fProfileGammaTPC_AntiLambda_AntiProton(NULL)   
{

  for(int i=0;i<3;i++){
    fCurrentVtx[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    fHistLambdaPt[i]  = NULL;            
    fHistLambdaEta[i] = NULL;
    fHistLambdaDcaToPrimVertex[i] = NULL;
    fHistLambdaCPA[i] = NULL;
    fHistLambdaDecayLength[i] = NULL;
    fHistLambdaMass[i] = NULL;
    fProfileLambdaMassVsPt[i] = NULL;
    fHistAntiLambdaPt[i] = NULL;
    fHistAntiLambdaEta[i] = NULL;
    fHistAntiLambdaDcaToPrimVertex[i] = NULL;
    fHistAntiLambdaCPA[i] = NULL;
    fHistAntiLambdaDecayLength[i] = NULL;
    fHistAntiLambdaMass[i] = NULL;
    fProfileAntiLambdaMassVsPt[i] = NULL;
  }  
  
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskGammaDeltaPID::~AliAnalysisTaskGammaDeltaPID()
{
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!

  //Delete the clones /// TList clones
  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;

  //delete the functions: //they are 'new'
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fV0CutPU)      delete fV0CutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU; 
      
}










//________________ Define Histograms _______________
void AliAnalysisTaskGammaDeltaPID::UserCreateOutputObjects()
{
  //std::cout<<"\n UserCreateOutputObject: function begins...\n"<<endl; 
  //Get The Input Hander:
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  
  if (!inputHandler) {  printf("\n***** ERROR *****\n Input handler missing, Status:QUIT!\n");    exit(1);}




  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);
  
  //// obtain the PID response object if needed:
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) {  printf("\n***** ERROR *****\n fPIDResponse missing, Status:QUIT!\n");    exit(1);}    
    
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  SetupEventAndTaskConfigInfo();


  /// Bunch of QA histograms:  
  fCentDistBeforCut = new TH1F("fCentDistBeforCut","Cent w/o any Cuts; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistBeforCut);
  fCentDistAfterCut = new TH1F("fCentDistAfterCut","Cent with all Cut; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistAfterCut);

  Int_t gMaxGlobalmult  = 4000;
  Int_t gMaxTPCcorrmult = 5000;
  Int_t gMaxESDtracks   = 20000;


  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MAfter);

  fHistCL0VsV0MBefore = new TH2F("fHistCL0VsV0MBefore","Before;Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHist->Add(fHistCL0VsV0MBefore);
  fHistCL0VsV0MAfter  = new TH2F("fHistCL0VsV0MAfter","After; Cent(V0M); Cent(CL0) ",100,0,100,100,0,100);
  fListHist->Add(fHistCL0VsV0MAfter);

  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHist->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHist->Add(fHistTPCVsESDTrkAfter);




  
  //----------- User's Analysis histograms Here: --------------
  Char_t  name[100];
  Char_t title[100];
 

  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90}; // Usual Bins for Observables
  
  Double_t fVzBinsZDC[21] = {-10,-8,-7,-6,-5,-4,-3,-2,-1,-0.5, 0, 0.5, 1.,2.,3.,4.,5.,6.,7.,8.,10.};
  Double_t fCentBinsZDC[41] = {0.};
  for(int i=1;i<41;i++)  fCentBinsZDC[i] = i*2.0;
  
  ///
  fAvgCosNPsivsCentEtaPos = new TProfile("fAvgCosNPsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgCosNPsivsCentEtaPos);
  fAvgSinNPsivsCentEtaPos = new TProfile("fAvgSinNPsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgSinNPsivsCentEtaPos);
  fAvgCosNPsivsCentEtaNeg = new TProfile("fAvgCosNPsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgCosNPsivsCentEtaNeg);
  fAvgSinNPsivsCentEtaNeg = new TProfile("fAvgSinNPsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",gHarmonic),90,0,90);
  fListHist->Add(fAvgSinNPsivsCentEtaNeg);

  fAvgCos3PsivsCentEtaPos = new TProfile("fAvgCos3PsivsCentEtaPos",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgCos3PsivsCentEtaPos);
  fAvgSin3PsivsCentEtaPos = new TProfile("fAvgSin3PsivsCentEtaPos",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgSin3PsivsCentEtaPos);
  fAvgCos3PsivsCentEtaNeg = new TProfile("fAvgCos3PsivsCentEtaNeg",Form("<cos(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgCos3PsivsCentEtaNeg);
  fAvgSin3PsivsCentEtaNeg = new TProfile("fAvgSin3PsivsCentEtaNeg",Form("<sin(%d#Psi)> vs cent",3),90,0,90);
  fListHist->Add(fAvgSin3PsivsCentEtaNeg);
  
  ///q-Vector for ESE:
  fHistV0CDetqVectorvsCent = new TH2F("fHistV0CDetqVectorvsCent","q V0Cdet  vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistV0CDetqVectorvsCent);
  fHistV0ADetqVectorvsCent = new TH2F("fHistV0ADetqVectorvsCent","q V0Adet  vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistV0ADetqVectorvsCent);  
  fHistTPCPosqVectorvsCent = new TH2F("fHistTPCPosqVectorvsCent","q TPC pos vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistTPCPosqVectorvsCent);
  fHistTPCNegqVectorvsCent = new TH2F("fHistTPCNegqVectorvsCent","q TPC neg vs Cent; Cent; q vect;",10,centRange,200,0,1.0);
  fListHist->Add(fHistTPCNegqVectorvsCent);


  
  fHistVertexZcm     = new TH1F("fHistVertexZcm"," V_{z}; V_{z} cm; events ",100,-10,10);
  fListHist->Add(fHistVertexZcm);
  fHistVxvsVzMinBias = new TProfile("fHistVxvsVzMinBias"," <Vx> vs Vz; V_{z}cm; <V_{x}> cm", 20,fVzBinsZDC);
  fListHist->Add(fHistVxvsVzMinBias);
  fHistVyvsVzMinBias = new TProfile("fHistVyvsVzMinBias"," <Vy> vs Vz; V_{z}cm; <V_{y}> cm", 20,fVzBinsZDC);
  fListHist->Add(fHistVyvsVzMinBias);


  ///Neutron Calorimeters:
  hAvgZNACh0vsCentVz = new TProfile2D("hAvgZNACh0vsCentVz"," <ZNA> Ch0; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNACh0vsCentVz);
  hAvgZNCCh0vsCentVz = new TProfile2D("hAvgZNCCh0vsCentVz"," <ZNC> Ch0; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNCCh0vsCentVz);
  hAvgZNACh1vsCentVz = new TProfile2D("hAvgZNACh1vsCentVz"," <ZNA> Ch1; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNACh1vsCentVz);
  hAvgZNCCh1vsCentVz = new TProfile2D("hAvgZNCCh1vsCentVz"," <ZNC> Ch1; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNCCh1vsCentVz);
  hAvgZNACh2vsCentVz = new TProfile2D("hAvgZNACh2vsCentVz"," <ZNA> Ch2; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNACh2vsCentVz);
  hAvgZNCCh2vsCentVz = new TProfile2D("hAvgZNCCh2vsCentVz"," <ZNC> Ch2; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNCCh2vsCentVz);
  hAvgZNACh3vsCentVz = new TProfile2D("hAvgZNACh3vsCentVz"," <ZNA> Ch3; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNACh3vsCentVz);
  hAvgZNCCh3vsCentVz = new TProfile2D("hAvgZNCCh3vsCentVz"," <ZNC> Ch3; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNCCh3vsCentVz);
  hAvgZNACh4vsCentVz = new TProfile2D("hAvgZNACh4vsCentVz"," <ZNA> Ch4; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNACh4vsCentVz);
  hAvgZNCCh4vsCentVz = new TProfile2D("hAvgZNCCh4vsCentVz"," <ZNC> Ch4; Cent; V_{z}(cm);",40,fCentBinsZDC,20,fVzBinsZDC);
  fListHist->Add(hAvgZNCCh4vsCentVz);


  ///Neutron Calorimeters:
  hAvgV0ChannelsvsVz = new TProfile2D("hAvgV0ChannelsvsVz"," <V0 Ch(n)> ; Channel;  V_{z}(cm);",64,0,64,20,fVzBinsZDC);
  fListHist->Add(hAvgV0ChannelsvsVz);


  hAvgQNXvsCentV0C = new TProfile("hAvgQNXvsCentV0C"," V0C <Q2x> vs Cent; Cent; <Q_{2,x}>",90,0,90);
  fListHist->Add(hAvgQNXvsCentV0C);
  hAvgQNYvsCentV0C = new TProfile("hAvgQNYvsCentV0C"," V0C <Q2y> vs Cent; Cent; <Q_{2,y}>",90,0,90);
  fListHist->Add(hAvgQNYvsCentV0C);
  hAvgQ3XvsCentV0C = new TProfile("hAvgQ3XvsCentV0C"," V0C <Q3x> vs Cent; Cent; <Q_{3,x}>",90,0,90);
  fListHist->Add(hAvgQ3XvsCentV0C);
  hAvgQ3YvsCentV0C = new TProfile("hAvgQ3YvsCentV0C"," V0C <Q3y> vs Cent; Cent; <Q_{3,y}>",90,0,90);
  fListHist->Add(hAvgQ3YvsCentV0C);
  hAvgQNXvsCentV0A = new TProfile("hAvgQNXvsCentV0A"," V0A <Q2x> vs Cent; Cent; <Q_{2,x}>",90,0,90);
  fListHist->Add(hAvgQNXvsCentV0A);
  hAvgQNYvsCentV0A = new TProfile("hAvgQNYvsCentV0A"," V0A <Q2y> vs Cent; Cent; <Q_{2,y}>",90,0,90);
  fListHist->Add(hAvgQNYvsCentV0A);
  hAvgQ3XvsCentV0A = new TProfile("hAvgQ3XvsCentV0A"," V0A <Q3x> vs Cent; Cent; <Q_{3,x}>",90,0,90);
  fListHist->Add(hAvgQ3XvsCentV0A);
  hAvgQ3YvsCentV0A = new TProfile("hAvgQ3YvsCentV0A"," V0A <Q3y> vs Cent; Cent; <Q_{3,y}>",90,0,90);  
  fListHist->Add(hAvgQ3YvsCentV0A);



  //// CME PID ANALYSIS Histograms/Profiles:

  ///TPC Event Planes:
  fHistTPCPsiNPosPlane = new TH2F("fHistTPCPsiNPosPlane",Form("#Psi_{n}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsiNPosPlane);
  fHistTPCPsiNNegPlane = new TH2F("fHistTPCPsiNNegPlane",Form("#Psi_{n}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsiNNegPlane);
  fHistTPCPsi3PosPlane = new TH2F("fHistTPCPsi3PosPlane",Form("#Psi_{3}(#eta+); centrality; #Psi_{%d,#eta+}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi3PosPlane);
  fHistTPCPsi3NegPlane = new TH2F("fHistTPCPsi3NegPlane",Form("#Psi_{3}(#eta-); centrality; #Psi_{%d,#eta-}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistTPCPsi3NegPlane);

  /// V0 Event Planes:
  fHistV0CPsiNEventPlane = new TH2F("fHistV0CPsiNEventPlane",Form("#Psi_{n}(V0C); centrality; #Psi_{%d,V0C}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0CPsiNEventPlane);
  fHistV0APsiNEventPlane = new TH2F("fHistV0APsiNEventPlane",Form("#Psi_{n}(V0A); centrality; #Psi_{%d,V0A}(rad); events",gHarmonic),18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0APsiNEventPlane);  
  fHistV0CPsi3EventPlane = new TH2F("fHistV0CPsi3EventPlane",Form("#Psi_{3}(V0C); centrality; #Psi_{%d,V0C}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0CPsi3EventPlane);
  fHistV0APsi3EventPlane = new TH2F("fHistV0APsi3EventPlane",Form("#Psi_{3}(V0A); centrality; #Psi_{%d,V0A}(rad); events",3), 18,0,90,50,0,3.14159);
  fListHist->Add(fHistV0APsi3EventPlane);  



  
  ///TPC Event Plane Correlations:
  hTPCPsiNCorrelation = new TProfile("hTPCPsiNCorrelation",Form("Psi_{%d} Res. vs Cent; Cent; Resolution",gHarmonic),90,0,90);
  //hTPCEPCorrelation->Sumw2();
  fListHist->Add(hTPCPsiNCorrelation);
  hTPCPsi3Correlation = new TProfile("hTPCPsi3Correlation",Form("Psi_{%d} Res. vs Cent; Cent; Resolution",gHarmonic),90,0,90);
  //hTPCEPCorrelation->Sumw2();
  fListHist->Add(hTPCPsi3Correlation);  

  ///V0X-V0X/TPC  Event Plane Correlations:
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
  


  

  //// 3p correlators PID:
  hAvg3pC112vsCentPP = new TProfile("hAvg3pC112vsCentPP"," <3p> vs Cent; Cent; <#gamma112>",90,0,90);
  hAvg3pC112vsCentPP->Sumw2();
  fListHist->Add(hAvg3pC112vsCentPP);
  hAvg3pC112vsCentNN = new TProfile("hAvg3pC112vsCentNN"," <3p> vs Cent; Cent; <#gamma112>",90,0,90);  
  hAvg3pC112vsCentNN->Sumw2();
  fListHist->Add(hAvg3pC112vsCentNN);
  hAvg3pC112vsCentOS = new TProfile("hAvg3pC112vsCentOS"," <3p> vs Cent; Cent; <#gamma112>",90,0,90);
  hAvg3pC112vsCentOS->Sumw2();
  fListHist->Add(hAvg3pC112vsCentOS);

  hAvg3pC123vsCentPP = new TProfile("hAvg3pC123vsCentPP"," <3p> vs Cent; Cent; <#gamma123>",90,0,90);
  hAvg3pC123vsCentPP->Sumw2();
  fListHist->Add(hAvg3pC123vsCentPP);
  hAvg3pC123vsCentNN = new TProfile("hAvg3pC123vsCentNN"," <3p> vs Cent; Cent; <#gamma123>",90,0,90);  
  hAvg3pC123vsCentNN->Sumw2();
  fListHist->Add(hAvg3pC123vsCentNN);
  hAvg3pC123vsCentOS = new TProfile("hAvg3pC123vsCentOS"," <3p> vs Cent; Cent; <#gamma123>",90,0,90);
  hAvg3pC123vsCentOS->Sumw2();
  fListHist->Add(hAvg3pC123vsCentOS);  

  /// 2p Correlator PID:  
  hAvgDelta1vsCentPP = new TProfile("hAvgDelta1vsCentPP"," Delta1 PP vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta1vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta1vsCentPP);
  hAvgDelta1vsCentNN = new TProfile("hAvgDelta1vsCentNN"," Delta1 NN vs Cent; Cent; <#delta_{m,n}>",90,0,90);    
  hAvgDelta1vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta1vsCentNN);  
  hAvgDelta1vsCentOS = new TProfile("hAvgDelta1vsCentOS"," Delta1 OS vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta1vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta1vsCentOS);

  hAvgDelta2vsCentPP = new TProfile("hAvgDelta2vsCentPP"," Delta2 PP vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta2vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta2vsCentPP);
  hAvgDelta2vsCentNN = new TProfile("hAvgDelta2vsCentNN"," Delta2 NN vs Cent; Cent; <#delta_{m,n}>",90,0,90);    
  hAvgDelta2vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta2vsCentNN);  
  hAvgDelta2vsCentOS = new TProfile("hAvgDelta2vsCentOS"," Delta2 OS vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta2vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta2vsCentOS);

  hAvgDelta3vsCentPP = new TProfile("hAvgDelta3vsCentPP"," Delta3 PP vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta3vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta3vsCentPP);
  hAvgDelta3vsCentNN = new TProfile("hAvgDelta3vsCentNN"," Delta3 NN vs Cent; Cent; <#delta_{m,n}>",90,0,90);    
  hAvgDelta3vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta3vsCentNN);  
  hAvgDelta3vsCentOS = new TProfile("hAvgDelta3vsCentOS"," Delta3 OS vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta3vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta3vsCentOS);

  hAvgDelta4vsCentPP = new TProfile("hAvgDelta4vsCentPP"," Delta4 PP vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta4vsCentPP->Sumw2();
  fListHist->Add(hAvgDelta4vsCentPP);
  hAvgDelta4vsCentNN = new TProfile("hAvgDelta4vsCentNN"," Delta4 NN vs Cent; Cent; <#delta_{m,n}>",90,0,90);    
  hAvgDelta4vsCentNN->Sumw2();
  fListHist->Add(hAvgDelta4vsCentNN);  
  hAvgDelta4vsCentOS = new TProfile("hAvgDelta4vsCentOS"," Delta4 OS vs Cent; Cent; <#delta_{m,n}>",90,0,90);
  hAvgDelta4vsCentOS->Sumw2();
  fListHist->Add(hAvgDelta4vsCentOS);
  
  
  

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


  ///// Case 0 = before cut, case 1 = afterCut.
  //Char_t name2[10];
  
  for(int i=0; i<2; i++){
    if(i==0)
      sprintf(name,"Before");
    else
       sprintf(name,"After");
    /// Lambdas:
    fHistLambdaPt[i] = new TH1D(Form("hLambdaPt_%sMassCut",name),"", 200, 0., 20.);
    fListHist->Add(fHistLambdaPt[i]);
    fHistLambdaEta[i] = new TH1D(Form("hLambdaEta_%sMassCut",name),"",200, -10., 10.);
    fListHist->Add(fHistLambdaEta[i]);
    fHistLambdaDcaToPrimVertex[i] = new TH1D(Form("hLambdaDcaToPrimVertex_%sMassCut",name),"",200, 0., 20.);
    fListHist->Add(fHistLambdaDcaToPrimVertex[i]);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",name),"",200, 0.9, 1.);
    fListHist->Add(fHistLambdaCPA[i]);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",name),"", 250, 0., 500.);
    fListHist->Add(fHistLambdaDecayLength[i]);
    fHistLambdaMass[i] = new TH1D(Form("hLambdaMass_%sMassCut",name),"",250,1.,1.25); //  Current bin size = 0.001
    fListHist->Add(fHistLambdaMass[i]);
    fProfileLambdaMassVsPt[i] = new TProfile(Form("pLambdaMassVsPt_%sMassCut",name),"",200,0,20);
    fListHist->Add(fProfileLambdaMassVsPt[i]);

    // AntiLambdas
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",name),"", 200, 0., 20.);
    fListHist->Add(fHistAntiLambdaPt[i]);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",name),"",200, -10., 10.);
    fListHist->Add(fHistAntiLambdaEta[i]);
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",name),"",200, 0., 20.);
    fListHist->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",name),"",200, 0.9, 1.);
    fListHist->Add(fHistAntiLambdaCPA[i]);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",name),"", 250, 0., 500.);
    fListHist->Add(fHistAntiLambdaDecayLength[i]);
    fHistAntiLambdaMass[i] = new TH1D(Form("hAntiLambdaMass_%sMassCut",name),"",250,1.,1.25); // Current bin size = 0.001
    fListHist->Add(fHistAntiLambdaMass[i]);
    fProfileAntiLambdaMassVsPt[i] = new TProfile(Form("pAntiLambdaMassVsPt_%sMassCut",name),"",200,0,20);
    fListHist->Add(fProfileAntiLambdaMassVsPt[i]);    
  }


  ///Gamma Correlators:
   ///Lambda - X
  fProfileGammaTPC_Lambda_hPos = new TProfile("fProfileGammaTPC_Lambda_hPos","",20,0,100);
  fListHist->Add(fProfileGammaTPC_Lambda_hPos);
  fProfileGammaTPC_Lambda_hNeg = new TProfile("fProfileGammaTPC_Lambda_hNeg","",20,0,100);
  fListHist->Add(fProfileGammaTPC_Lambda_hNeg);
  fProfileGammaTPC_Lambda_Proton = new TProfile("fProfileGammaTPC_Lambda_Proton","",20,0,100);
  fListHist->Add(fProfileGammaTPC_Lambda_Proton);
  fProfileGammaTPC_Lambda_AntiProton = new TProfile("fProfileGammaTPC_Lambda_AntiProton","",20,0,100);  
  fListHist->Add(fProfileGammaTPC_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaTPC_AntiLambda_hPos = new TProfile("fProfileGammaTPC_AntiLambda_hPos","",20,0,100);
  fListHist->Add(fProfileGammaTPC_AntiLambda_hPos);
  fProfileGammaTPC_AntiLambda_hNeg = new TProfile("fProfileGammaTPC_AntiLambda_hNeg","",20,0,100);
  fListHist->Add(fProfileGammaTPC_AntiLambda_hNeg);
  fProfileGammaTPC_AntiLambda_Proton = new TProfile("fProfileGammaTPC_AntiLambda_Proton","",20,0,100);
  fListHist->Add(fProfileGammaTPC_AntiLambda_Proton);
  fProfileGammaTPC_AntiLambda_AntiProton = new TProfile("fProfileGammaTPC_AntiLambda_AntiProton","",20,0,100);  
  fListHist->Add(fProfileGammaTPC_AntiLambda_AntiProton);

  ///// <-------- chunzheng's Histograms....
  
  
    
  //// ====> Correction maps: To be implemented later:

  if(fListTRKCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist found for MC tracking Efficiency.\n"<<std::endl;
  }
  else{
    std::cout<<"\n ******* WARNING *****\n No TList for MC Efficiency Correction..!!\n using TrkWgt = 1.0 \n "<<std::endl;
  }

  if(fListNUACorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist found for NUA Correction.\n"<<std::endl;
    //fListNUACorr->ls(); 
  }
  else{
    std::cout<<"\n ******* WARNING *****\n No TList NUA Correction..!!\n using NUAWgt = 1.0 \n "<<std::endl;
  }

  if(fListV0MCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist fount for V0 gain Correction.\n"<<std::endl;
    //fListV0MCorr->ls();
  }
  else{
    std::cout<<"\n ******* WARNING *****\n No TList V0 gain Correction..!! using V0 Ch.gains = 1.0 \n "<<std::endl;
  }


  std::cout<<"\n ==========> UserCreateOutputObject::Info() The Analysis Settings: <============ "<<std::endl;  

  char pairname[10];//

  if(gParticleID==0){
    sprintf(pairname,"h-h");
  }
  else if(gParticleID==1){
    sprintf(pairname,"pi-pi");
  }
  else if(gParticleID==2){
    sprintf(pairname,"K-K");
  }
  else if(gParticleID==3){
    sprintf(pairname,"p-p");
  }
  else if(gParticleID==10){
    sprintf(pairname,"pi-Ch");
  }
  else if(gParticleID==20){
    sprintf(pairname,"K-Ch");
  }
  else if(gParticleID==30){
    sprintf(pairname,"P-Ch");
  }
  else if(bAnalysLambdaPairs){
    sprintf(pairname,"Lambda-x"); 
    std::cout<<"Lmbda CPAcut: "<<fV0CPAMin<<"\tMassCut: "<<fLambdaMassCut<<"\tDCAtoVtx: "<<fV0DCAToPrimVtxMax<<"\tNsigmaDaut: "<<fDaughtersNsigma<<"\tDautDCA: "<<fDaughtersDCAToPrimVtxMin<<std::endl;
  }
  else{
    sprintf(pairname,"ch-ch");
  }


  if(sDetectorforEP.Contains("V0")){
    if(sDetectorforEP.Contains("V0C")){
      fHistAnalysisInfo->SetBinContent(11,1);
    }
    else{
      fHistAnalysisInfo->SetBinContent(12,1);
    }
  }
  else{
    fHistAnalysisInfo->SetBinContent(13,1);
  }
  
  
  std::cout<<"\n FB: "<<fFilterBit<<", Ncls: "<<fTPCclustMin<<", Harmonic: "<<gHarmonic<<", pair: "<<pairname<<" detforEP = "<<sDetectorforEP<<",\n\n"<<std::endl;



  
  ////==========> LHC18q/r PileUp Removal Functions: ---- Do not Remove them !!! -----
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
  
  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
  //--------------------------------------------------------------------------------------





  
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




  
  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  Float_t centrCL1   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 
  centrality = centrV0M;  // This Is Default, changes below for other options set in AddTask Macro:
  
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
  




  //------ Pile up cut function called------

  Bool_t kPileupEvent = kFALSE;

  kPileupEvent = CheckEventIsPileUp2018(fAOD);
  if(kPileupEvent)  return;
    

  if(centrality<fCentralityMin || centrality>fCentralityMax){ 
    return;
  }





  

  Int_t ntracks = fAOD->GetNumberOfTracks();
  if(ntracks < 4) return;                        // Minimum 4 tracks per event!!

 // 1 = Pion, 2 = Kaon, 3 = Proton, 4 - Lambda, To be Set From AddTask
  Int_t gRequiredPID = 1;
  //Int_t gRequiredPID2 = 1;  // For now we consider like particles. 

  

  //////----> Get Magnetic field and RunNo.---------
  // Float_t fMagField = fAOD->GetMagneticField();
  Int_t runNumber = fAOD->GetRunNumber();
  //------------------------------------------------


  //Load NUA and V0M correction map run by run: And only if Runnumber changes.

  if(runNumber!=gOldRunNumber) {
   
    GetNUACorrectionHist(runNumber,0); // kPID = 0 = charge 

    if(fListTRKCorr) GetMCCorrectionHist();  //gRequiredPID is redundant. We read in all PID and Charge Efficiency.
    
    if(fListV0MCorr) {
     GetV0MCorrectionHist(runNumber);
    }

    gOldRunNumber = runNumber;
  }
  ///gOldRunNumber is global. it remembers which run was used in last event, So that
  ///we dont waste time re-setting calibration histogram for each event. 
  //------------------------------------------
 












  

  Float_t fMultTPCFull = 0;  // TPC mult estimate
  Float_t fMultGlobal  = 0;  // global track multiplicity
  

  


  
 


  Double_t gPsiN = gHarmonic;        /// gHarmonic is  EP harmonic. set in AddTask


  //=============== Get the V0 data ====================
  
  const AliAODVZERO *fAODV0 = (AliAODVZERO *) fAOD->GetVZEROData();
  Float_t fMultV0 = 0;
  Float_t fPhiV0;
  Float_t fV0chGain = 1.0;

  Double_t fQxV0CHarmN=0,fQyV0CHarmN=0,fQxV0AHarmN=0,fQyV0AHarmN=0;
  Double_t fQxV0CHarm3=0,fQyV0CHarm3=0,fQxV0AHarm3=0,fQyV0AHarm3=0;  
  Double_t fSumMV0A = 0,fSumMV0C = 0;
  Double_t fSelectedV0PsiN=0.,fSelectedV0Psi3=0.;
  Int_t    ibinV0, icentbin=1;

  
  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

    fMultV0 = fAODV0->GetMultiplicity(iV0);

    /// V0 Channel Gain Correction:
    if(fHCorrectV0ChWeghts){ 
      ibinV0 = fHCorrectV0ChWeghts->FindBin(pVtxZ,iV0);
      fV0chGain = fHCorrectV0ChWeghts->GetBinContent(ibinV0); 
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    
    hAvgV0ChannelsvsVz->Fill(iV0+0.5, pVtxZ, fMultV0);

    fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);

    if(iV0 < 32){
      fQxV0CHarmN  += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      fQyV0CHarmN  += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fQxV0CHarm3  += TMath::Cos(3*fPhiV0) * fMultV0;       /// TO be removed after full pass with Gain correction.
      fQyV0CHarm3  += TMath::Sin(3*fPhiV0) * fMultV0;
      fSumMV0C     += fMultV0;
    }
    else if(iV0 >= 32){
      fQxV0AHarmN  += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      fQyV0AHarmN  += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fQxV0AHarm3  += TMath::Cos(3*fPhiV0) * fMultV0;       /// TO be removed after full pass with Gain correction.
      fQyV0AHarm3  += TMath::Sin(3*fPhiV0) * fMultV0;
      fSumMV0A     += fMultV0;
    } 
  }


  if(fSumMV0A<=1e-4 || fSumMV0C<=1e-4)   return;         /// this means there is not enough tracks in V0 Detectors in this event!!


  //Scaled by Multiplicity to remove Mult-fluctuations:
  Double_t fQnxV0C = fQxV0CHarmN/fSumMV0C;
  Double_t fQnyV0C = fQyV0CHarmN/fSumMV0C;  
  Double_t fQnxV0A = fQxV0AHarmN/fSumMV0A;
  Double_t fQnyV0A = fQyV0AHarmN/fSumMV0A;

  Double_t fQ3xV0C = fQxV0CHarm3/fSumMV0C;
  Double_t fQ3yV0C = fQyV0CHarm3/fSumMV0C;  
  Double_t fQ3xV0A = fQxV0AHarm3/fSumMV0A;
  Double_t fQ3yV0A = fQyV0AHarm3/fSumMV0A;  


  /// Correct the Q vector from V0 detectors:

  ///<Q> correction for PsiN:
  if(fHCorrectQNxV0C && fHCorrectQNyV0C){ ///prevents Code break!
    icentbin = fHCorrectQNxV0C->FindBin(centrCL1);   ///Hardcode CL1 centrality to avoid mismatched Centrality estimator! 
    fQnxV0C -= fHCorrectQNxV0C->GetBinContent(icentbin); 
    fQnyV0C -= fHCorrectQNyV0C->GetBinContent(icentbin);   /// ***Careful with Names => Qx and Qy DONT MISS-match!!!
  }
  if(fHCorrectQNxV0A && fHCorrectQNyV0A){ 
    icentbin = fHCorrectQNxV0A->FindBin(centrCL1);
    fQnxV0A -= fHCorrectQNxV0A->GetBinContent(icentbin); 
    fQnyV0A -= fHCorrectQNyV0A->GetBinContent(icentbin);  
  }

  ///<Q> correction for Psi3:
  if(fHCorrectQ3xV0C && fHCorrectQ3yV0C){
    icentbin = fHCorrectQ3xV0C->FindBin(centrCL1);
    fQ3xV0C -= fHCorrectQ3xV0C->GetBinContent(icentbin); 
    fQ3yV0C -= fHCorrectQ3yV0C->GetBinContent(icentbin);   
  }
  if(fHCorrectQ3xV0A && fHCorrectQ3yV0A){
    icentbin = fHCorrectQ3xV0A->FindBin(centrCL1);
    fQ3xV0A -= fHCorrectQ3xV0A->GetBinContent(icentbin); 
    fQ3yV0A -= fHCorrectQ3yV0A->GetBinContent(icentbin);  
  }
  
  
  ////---- fill the <Q> vector from V0A/C vs Cent:-----------
  hAvgQNXvsCentV0C->Fill(centrCL1,fQnxV0C); /// to Avoid self correlation, V0 <Q> is filled with CL1 centrality! 
  hAvgQNYvsCentV0C->Fill(centrCL1,fQnyV0C);  
  hAvgQNXvsCentV0A->Fill(centrCL1,fQnxV0A);
  hAvgQNYvsCentV0A->Fill(centrCL1,fQnyV0A);

  hAvgQ3XvsCentV0C->Fill(centrCL1,fQ3xV0C);
  hAvgQ3YvsCentV0C->Fill(centrCL1,fQ3yV0C);  
  hAvgQ3XvsCentV0A->Fill(centrCL1,fQ3xV0A);
  hAvgQ3YvsCentV0A->Fill(centrCL1,fQ3yV0A);  



  ///---- Get the <Q> vector Correction for V0------
  ///. To be implemented when <Q> is available after 2nd Pass..


  
  /// Get V0A and V0C Event Planes:
  Double_t fPsiNV0C = 0., fPsiNV0A = 0.;
  Double_t fPsi3V0C = 0., fPsi3V0A = 0.;

    
  if(fSumMV0C > 0){
    fPsiNV0C = (1./gPsiN)*TMath::ATan2(fQnyV0C,fQnxV0C);
    if(fPsiNV0C < 0) fPsiNV0C += TMath::TwoPi()/gPsiN;

    fPsi3V0C = (1./3)*TMath::ATan2(fQ3yV0C,fQ3xV0C);
    if(fPsi3V0C < 0) fPsi3V0C += TMath::TwoPi()/3;      
  }
  if(fSumMV0A > 0){
    fPsiNV0A = (1./gPsiN)*TMath::ATan2(fQnyV0A,fQnxV0A);
    if(fPsiNV0A < 0) fPsiNV0A += TMath::TwoPi()/gPsiN;
    
    fPsi3V0A = (1./3)*TMath::ATan2(fQ3yV0A,fQ3xV0A);
    if(fPsi3V0A < 0) fPsi3V0A += TMath::TwoPi()/3;      
  }  

  
  fHistV0CPsiNEventPlane->Fill(centrality,fPsiNV0C);
  fHistV0APsiNEventPlane->Fill(centrality,fPsiNV0A);
  fHistV0CPsi3EventPlane->Fill(centrality,fPsi3V0C);
  fHistV0APsi3EventPlane->Fill(centrality,fPsi3V0A);    


  /// if V0 Event Plane is not selected then Default is TPC EP
  
  if(sDetectorforEP.Contains("V0"))
    bUseV0EventPlane = kTRUE;
  
  if(bUseV0EventPlane){
    if(sDetectorforEP.Contains("V0C")){
      fSelectedV0PsiN = fPsiNV0C;
      fSelectedV0Psi3 = fPsi3V0C;
    }
    else if(sDetectorforEP.Contains("V0A")){
      fSelectedV0PsiN = fPsiNV0A;
      fSelectedV0Psi3 = fPsi3V0A;
    }
    else{ /// Default is V0A EP
      fSelectedV0PsiN = fPsiNV0A;
      fSelectedV0Psi3 = fPsi3V0A;
    }
  }


  

  
  //=============== Get the ZDC data ====================
 
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
  
    hAvgZNACh0vsCentVz->Fill(centrality,pVtxZ,fZNATowerRawAOD[0]);
    hAvgZNCCh0vsCentVz->Fill(centrality,pVtxZ,fZNCTowerRawAOD[0]);
    hAvgZNACh1vsCentVz->Fill(centrality,pVtxZ,fZNATowerRawAOD[1]);
    hAvgZNCCh1vsCentVz->Fill(centrality,pVtxZ,fZNCTowerRawAOD[1]);    
    hAvgZNACh2vsCentVz->Fill(centrality,pVtxZ,fZNATowerRawAOD[2]);
    hAvgZNCCh2vsCentVz->Fill(centrality,pVtxZ,fZNCTowerRawAOD[2]);
    hAvgZNACh3vsCentVz->Fill(centrality,pVtxZ,fZNATowerRawAOD[3]);
    hAvgZNCCh3vsCentVz->Fill(centrality,pVtxZ,fZNCTowerRawAOD[3]);    
    hAvgZNACh4vsCentVz->Fill(centrality,pVtxZ,fZNATowerRawAOD[4]);
    hAvgZNCCh4vsCentVz->Fill(centrality,pVtxZ,fZNCTowerRawAOD[4]);
        
  }
  
  

  //// ZDC Event plane: To be Implemented Later....!!



  ///------------ Get the TPC Event Planes........

  
  ///Track variables:
  Float_t trk1Pt=0,trk1Phi=0,trk1Eta=0,trk1DCAxy=0.0, trk1DCAz=0.0,trk1Chi2=0,trk1dEdx=0,trk1Wgt=1.0;
  Float_t trk2Pt=0,trk2Phi=0,trk2Eta=0,trk2DCAxy=0.0, trk2DCAz=0.0,trk2Chi2=0,trk2dEdx=0,trk2Wgt=1.0;  
  Int_t   trk1Chrg=0, trk1TpcNC=0;
  Int_t   trk2Chrg=0, trk2TpcNC=0;
  
  
  ////PID variables:
  Double_t nSigTOFpiontrk1=-99, nSigTPCpiontrk1=-99;
  Double_t nSigTOFkaontrk1=-99, nSigTPCkaontrk1=-99;
  Double_t nSigTOFprottrk1=-99, nSigTPCprottrk1=-99;
  
  Double_t nSigTOFpiontrk2=-99, nSigTPCpiontrk2=-99;
  Double_t nSigTOFkaontrk2=-99, nSigTPCkaontrk2=-99;
  Double_t nSigTOFprottrk2=-99, nSigTPCprottrk2=-99;
  
  Bool_t   isItPiontrk1 = kFALSE, isItKaontrk1 = kFALSE, isItProttrk1 = kFALSE;
  Bool_t   isItPiontrk2 = kFALSE, isItKaontrk2 = kFALSE, isItProttrk2 = kFALSE;
  
  Double_t ptWgtMCtrk1 = 1.0, WgtNUAtrk1 = 1.0, ptWgtMCtrk2 = 1.0, WgtNUAtrk2 = 1.0;
  Double_t wgt1particle=1.0, wgt2particle = 1.0; 

  
  Int_t    ptBinMC = 1, iBinNUA = 1;
  Double_t binCont = 1.0;
  
  Float_t fWgtEvent = 1.0; //Event Weight if Any.



 //// Variables for TPC EP calculation:
  Double_t SumQnxTPCPos = 0.,SumQnyTPCPos = 0;
  Double_t SumQnxTPCNeg = 0.,SumQnyTPCNeg = 0;
  ///For Psi3 :
  Double_t SumQ3xTPCPos = 0.,SumQ3yTPCPos = 0;
  Double_t SumQ3xTPCNeg = 0.,SumQ3yTPCNeg = 0;
  
  Double_t fWgtMultTPCPos=0.,fWgtMultTPCNeg=0;
  
  Int_t   fTPCclustMinforEP= 70;     //Fixed for EP calculation
  Float_t fMinPtCutforEP   = 0.2;    //Fixed for EP calculation
  Float_t fMaxPtCutforEP   = 2.0;    //Fixed for EP calculation
  Float_t fEtaGapPosforEP  = 0.1;    //could be made variable in AddTask Macro..
  Float_t fEtaGapNegforEP  =-0.1;    //could be made variable in AddTask Macro..
  Float_t fTrkChi2MinforEP = 0.1;    //Fixed for EP calculation
  Float_t fTrkChi2MaxforEP = 4.0;    //Fixed for EP calculation


  vector<Int_t>  vecPosEPTrkID;
  vector<Int_t>  vecNegEPTrkID;
  ///----------> Starting track Loop for TPC EVENT Plane -----------
  
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

      
      //Apply track cuts for EP  here:
      if((trk1Pt <= fMaxPtCutforEP) && (trk1Pt >= fMinPtCutforEP) && (trk1Eta >= -0.8) && (trk1Eta <= 0.8) && (trk1dEdx >= 10) && (trk1TpcNC >= fTPCclustMinforEP) && (trk1Chi2 >= fTrkChi2MinforEP) && (trk1Chi2 <= fTrkChi2MaxforEP) && TMath::Abs(trk1Chrg)) {


	//ptWgtMCtrk1 = 1.0;  
	WgtNUAtrk1 = 1.0;   

	if(trk1Chrg>0){
	  if(fHCorrectNUAChrgPos){ /// safety measures for breaks!
	    iBinNUA = fHCorrectNUAChrgPos->FindBin(pVtxZ,trk1Phi,trk1Eta);
	    WgtNUAtrk1 = fHCorrectNUAChrgPos->GetBinContent(iBinNUA);
	  }
	  //if(fHCorrectMCPosChrg){
	  //iBinNUA = fHCorrectMCPosChrg->FindBin(trk1Pt);
	  //ptWgtMCtrk1 = fHCorrectMCPosChrg->GetBinContent(iBinNUA);
	  //}
	}
	else{
	  if(fHCorrectNUAChrgNeg){ /// safety measures for breaks!
	    iBinNUA = fHCorrectNUAChrgNeg->FindBin(pVtxZ,trk1Phi,trk1Eta);
	    WgtNUAtrk1 = fHCorrectNUAChrgNeg->GetBinContent(iBinNUA);
	  }
	  //if(fHCorrectMCNegChrg){
	  //iBinNUA = fHCorrectMCNegChrg->FindBin(trk1Pt);
	  //ptWgtMCtrk1 = fHCorrectMCNegChrg->GetBinContent(iBinNUA);
	  //}	  
	}
	  
	///wgt1particle = WgtNUAtrk1*ptWgtMCtrk1;
	
	wgt1particle = WgtNUAtrk1;   /// No MC wgt for event plane, As we used charged particles. 

	
	Int_t trk1ID = AODtrack1->GetID(); //unique in a event
	
	
	///Used Pt as weight for Better resolution:
	if(trk1Eta >= fEtaGapPosforEP){
	  SumQnxTPCPos   += trk1Pt*wgt1particle*TMath::Cos(gPsiN*trk1Phi);
	  SumQnyTPCPos   += trk1Pt*wgt1particle*TMath::Sin(gPsiN*trk1Phi);
	  SumQ3xTPCPos   += trk1Pt*wgt1particle*TMath::Cos(3*trk1Phi);
	  SumQ3yTPCPos   += trk1Pt*wgt1particle*TMath::Sin(3*trk1Phi);
	  
	  fWgtMultTPCPos += trk1Pt*wgt1particle;
	  vecPosEPTrkID.push_back(trk1ID);
	}
	else if(trk1Eta <= fEtaGapNegforEP){
	  SumQnxTPCNeg   += trk1Pt*wgt1particle*TMath::Cos(gPsiN*trk1Phi);
	  SumQnyTPCNeg   += trk1Pt*wgt1particle*TMath::Sin(gPsiN*trk1Phi);
	  SumQ3xTPCNeg   += trk1Pt*wgt1particle*TMath::Cos(3*trk1Phi);
	  SumQ3yTPCNeg   += trk1Pt*wgt1particle*TMath::Sin(3*trk1Phi);
	  
	  fWgtMultTPCNeg += trk1Pt*wgt1particle;
	  vecNegEPTrkID.push_back(trk1ID);
	}
	 	
     }//----> track loop => All trackCuts applied.     
    }//-----> track loop => FB is validated.    
  }///------> track loop Ends here.<--------





  

  if(fWgtMultTPCPos<0.1 || fWgtMultTPCNeg<0.1) return;         /// this means there is not enough tracks in this event!!


  //Scaled by Multiplicity to remove Event by event fluctuation.
  Double_t fQnxPos = SumQnxTPCPos/fWgtMultTPCPos;
  Double_t fQnyPos = SumQnyTPCPos/fWgtMultTPCPos;
  Double_t fQnxNeg = SumQnxTPCNeg/fWgtMultTPCNeg;
  Double_t fQnyNeg = SumQnyTPCNeg/fWgtMultTPCNeg;

  Double_t fQ3xPos = SumQ3xTPCPos/fWgtMultTPCPos;
  Double_t fQ3yPos = SumQ3yTPCPos/fWgtMultTPCPos;
  Double_t fQ3xNeg = SumQ3xTPCNeg/fWgtMultTPCNeg;
  Double_t fQ3yNeg = SumQ3yTPCNeg/fWgtMultTPCNeg;  
		    

  ///------ Apply TPC <Q> vector corrections:--------
  
  ///<Q> correction for Psi2:
  if(fHCorrectTPCQNxEtaPos && fHCorrectTPCQNyEtaPos){ ///prevents Code break!
    icentbin = fHCorrectTPCQNxEtaPos->FindBin(centrV0M); //centrV0M is hardcoded to avoid mismatch in Estimator.
    fQnxPos -= fHCorrectTPCQNxEtaPos->GetBinContent(icentbin); 
    fQnyPos -= fHCorrectTPCQNyEtaPos->GetBinContent(icentbin);   /// ***Careful with Names => Qx and Qy DONT MISS-match!!!
  }
  if(fHCorrectTPCQNxEtaNeg && fHCorrectTPCQNyEtaNeg){ 
    icentbin = fHCorrectTPCQNxEtaNeg->FindBin(centrV0M);
    fQnxNeg -= fHCorrectTPCQNxEtaNeg->GetBinContent(icentbin);
    fQnyNeg -= fHCorrectTPCQNyEtaNeg->GetBinContent(icentbin); 
  }

  ///<Q> correction for Psi3:
  if(fHCorrectTPCQ3xEtaPos && fHCorrectTPCQ3yEtaPos){ 
    icentbin = fHCorrectTPCQ3xEtaPos->FindBin(centrV0M);
    fQ3xPos -= fHCorrectTPCQ3xEtaPos->GetBinContent(icentbin); 
    fQ3yPos -= fHCorrectTPCQ3yEtaPos->GetBinContent(icentbin);   
  }
  if(fHCorrectTPCQ3xEtaNeg && fHCorrectTPCQ3yEtaNeg){ 
    icentbin = fHCorrectTPCQ3xEtaNeg->FindBin(centrV0M);
    fQ3xNeg -= fHCorrectTPCQ3xEtaNeg->GetBinContent(icentbin);
    fQ3yNeg -= fHCorrectTPCQ3yEtaNeg->GetBinContent(icentbin);   
  }
  //-------------------------------------------------


  
  /// <Q> vector for Psi2 TPC event plane:
  fAvgCosNPsivsCentEtaPos->Fill(centrV0M,fQnxPos);   //centrV0M is hardcoded to avoid mismatch in Estimator.
  fAvgSinNPsivsCentEtaPos->Fill(centrV0M,fQnyPos); 
  fAvgCosNPsivsCentEtaNeg->Fill(centrV0M,fQnxNeg);
  fAvgSinNPsivsCentEtaNeg->Fill(centrV0M,fQnyNeg);

  /// <Q> vector for Psi3:
  fAvgCos3PsivsCentEtaPos->Fill(centrV0M,fQ3xPos); 
  fAvgSin3PsivsCentEtaPos->Fill(centrV0M,fQ3yPos); 
  fAvgCos3PsivsCentEtaNeg->Fill(centrV0M,fQ3xNeg);
  fAvgSin3PsivsCentEtaNeg->Fill(centrV0M,fQ3yNeg);


  

  //// Get TPC Event Planes:
  Double_t fPsiNTPCPos = 0., fPsiNTPCNeg = 0.;
  Double_t fPsi3TPCPos = 0., fPsi3TPCNeg = 0.;

    
  
  if(fWgtMultTPCPos > 0){
    fPsiNTPCPos = (1./gPsiN)*TMath::ATan2(fQnyPos,fQnxPos);
    if(fPsiNTPCPos < 0) fPsiNTPCPos += TMath::TwoPi()/gPsiN;

    fPsi3TPCPos = (1./3)*TMath::ATan2(fQ3yPos,fQ3xPos);
    if(fPsi3TPCPos < 0) fPsi3TPCPos += TMath::TwoPi()/3;      
  }
  if(fWgtMultTPCNeg > 0){
    fPsiNTPCNeg = (1./gPsiN)*TMath::ATan2(fQnyNeg,fQnxNeg);
    if(fPsiNTPCNeg < 0) fPsiNTPCNeg += TMath::TwoPi()/gPsiN;

    fPsi3TPCNeg = (1./3)*TMath::ATan2(fQ3yNeg,fQ3xNeg);
    if(fPsi3TPCNeg < 0) fPsi3TPCNeg += TMath::TwoPi()/3;     
  }  


  
  fHistTPCPsiNPosPlane->Fill(centrality,fPsiNTPCPos);
  fHistTPCPsiNNegPlane->Fill(centrality,fPsiNTPCNeg);
  fHistTPCPsi3PosPlane->Fill(centrality,fPsi3TPCPos);
  fHistTPCPsi3NegPlane->Fill(centrality,fPsi3TPCNeg);  
  
  hTPCPsiNCorrelation->Fill(centrality,TMath::Cos(gPsiN*fPsiNTPCPos - gPsiN*fPsiNTPCNeg)); ///TPC Psi_N Resolution
  hTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsi3TPCPos - 3*fPsi3TPCNeg));         ///TPC Psi3 Resolution


  /// Get Full TPC EP for V0-Correlation:
  Double_t fQnxTPCFull = fQnxPos + fQnxNeg;   //By corrections applied above they are already Recentered.
  Double_t fQnyTPCFull = fQnyPos + fQnyNeg;
  
  
  hV0CV0APsiNCorrelation->Fill(centrality,TMath::Cos(gPsiN*fPsiNV0A - gPsiN*fPsiNV0C));
  hV0CTPCPsiNCorrelation->Fill(centrality,TMath::Cos(gPsiN*fPsiNTPCPos - gPsiN*fPsiNV0C));
  hV0ATPCPsiNCorrelation->Fill(centrality,TMath::Cos(gPsiN*fPsiNTPCPos - gPsiN*fPsiNV0A));

  hV0CV0APsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNV0A - 3*fPsiNV0C));
  hV0CTPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsiNV0C));
  hV0ATPCPsi3Correlation->Fill(centrality,TMath::Cos(3*fPsiNTPCPos - 3*fPsiNV0A));

  //// Store q-vectors for ESE: (Make sure the qx,qy vectors are recentered in the above..!!)

  
  /// For V0X, we use CL1 Centrality removes autocorrelation Biases  // ***hardcoded to prevent accidentally using another Estimator.!
  fHistV0CDetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0C*fQnxV0C + fQnyV0C*fQnyV0C));
  fHistV0ADetqVectorvsCent->Fill(centrCL1,TMath::Sqrt(fQnxV0A*fQnxV0A + fQnyV0A*fQnyV0A));
  
  /// For TPC, we use V0 Centrality   // ***hardcoded to prevent accidentally using another Estimator.!
  fHistTPCPosqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQnxPos*fQnxPos + fQnyPos*fQnyPos));
  fHistTPCNegqVectorvsCent->Fill(centrV0M,TMath::Sqrt(fQnxNeg*fQnxNeg + fQnyNeg*fQnyNeg));
  

  






  
  if(bSkipAnalysis){
    //Okay we are skipping the nested loops and other Analyses.
    //So just fill QAs here and get out..
    fHistVertexZcm->Fill(pVtxZ);
    fHistVxvsVzMinBias->Fill(pVtxZ,pVtxX);
    fHistVyvsVzMinBias->Fill(pVtxZ,pVtxY);
    fCentDistAfterCut->Fill(centrality);  
  
    PostData(1,fListHist);

    /// Note: this Option is added so that we can run Faster a full pass over data
    /// to save any Any Q-vector etc, which is event-by-event quantity.
    return; 
  }













  

  ///I can Store PID in arrays to save time. But this works fast as well:

  
  //@chunzheng vector to contain Charge Particle information for Lambda-x pairing..
  vector<Double_t> vecPt;
  vector<Double_t> vecEta;
  vector<Double_t> vecPhi;
  vector<Int_t>    vecID;
  vector<Int_t>  vecPDGCode;
  vector<Double_t> vecNUAWeight;  //Charge
  vector<Double_t> vecNUEWeight;   //Charge
  vector<Double_t> vecNUAWeightPID;  //Charge
  vector<Double_t> vecNUEWeightPID;   //Charge

 
  Int_t kPIDtrk1=gParticleID;   /// gParticleID is Set From AddTask
  Int_t kPIDtrk2=gParticleID;   /// 0 = hadron (h-h), 1 = Pi-Pi, 2 = K-K, 3 = Prot-Prot, 

  /// If we want single Identified:    
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

  
  
  Bool_t bPIDoktrk1=kFALSE, bPIDoktrk2=kFALSE;
  Double_t ptWgtMCPIDtrk1 = 1.0, WgtNUAPIDtrk1 = 1.0;
  Double_t ptWgtMCPIDtrk2 = 1.0, WgtNUAPIDtrk2 = 1.0;
  Double_t wgt1PIDparticle=1.0, wgt2PIDparticle = 1.0; 

  
  
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
      if((trk1Pt <= fMaxPtCut) && (trk1Pt >= fMinPtCut) && (trk1Eta <= fMaxEtaCut) && (trk1Eta >= fMinEtaCut) && (trk1dEdx >= fdEdxMin) && (trk1TpcNC >= fTPCclustMin) && (trk1Chi2 >= fTrkChi2Min) && (trk1Chi2 <= fChi2Max) && TMath::Abs(trk1Chrg)) {





	ptWgtMCtrk1 = 1.0;   
	WgtNUAtrk1  = 1.0;   

	ptWgtMCPIDtrk1 = 1.0;
	WgtNUAPIDtrk1 = 1.0;
	
	/// Rihan: I should move the following part to a function (todo)
	
	//------------ Get the NUA and MC weight for Track1-----------------
	if(trk1Chrg>0){
	  if(fHCorrectNUAChrgPos){ /// safety measures for breaks!
	    iBinNUA = fHCorrectNUAChrgPos->FindBin(pVtxZ,trk1Phi,trk1Eta);
	    WgtNUAtrk1 = fHCorrectNUAChrgPos->GetBinContent(iBinNUA);
	  }
	  if(fHCorrectNUAkPIDPos){
	    //iBinNUA = fHCorrectNUAChrgPos->FindBin(pVtxZ,trk1Phi,trk1Eta); //same bin 
	    WgtNUAPIDtrk1 = fHCorrectNUAkPIDPos->GetBinContent(iBinNUA);
	  }

	  ///--------- Charge & PID NUE: -------------
	  if(fHCorrectMCPosChrg){
	    iBinNUA = fHCorrectMCPosChrg->FindBin(trk1Pt);  // re-use iBinNUA
	    ptWgtMCtrk1 = fHCorrectMCPosChrg->GetBinContent(iBinNUA);
	  }
	  //now PID NUE:
	  if(kPIDtrk1==1 && fHCorrectMCPosPion){
	    iBinNUA = fHCorrectMCPosPion->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCPosPion->GetBinContent(iBinNUA);
	  }
	  else if(kPIDtrk1==2 && fHCorrectMCPosKaon){
	    iBinNUA = fHCorrectMCPosKaon->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCPosKaon->GetBinContent(iBinNUA);
	  }
	  else if(kPIDtrk1==3 && fHCorrectMCPosProt){
	    iBinNUA = fHCorrectMCPosProt->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCPosProt->GetBinContent(iBinNUA);
	  }//----------------------------------	  
	}
	else{
	  if(fHCorrectNUAChrgNeg){ /// safety measures for breaks!
	    iBinNUA = fHCorrectNUAChrgNeg->FindBin(pVtxZ,trk1Phi,trk1Eta);
	    WgtNUAtrk1 = fHCorrectNUAChrgNeg->GetBinContent(iBinNUA);
	  }
	  if(fHCorrectNUAkPIDNeg){
	    //iBinNUA = fHCorrectNUAChrgNeg->FindBin(pVtxZ,trk1Phi,trk1Eta); //same bin 
	    WgtNUAPIDtrk1 = fHCorrectNUAkPIDNeg->GetBinContent(iBinNUA);
	  }

	  ///--------- Charge & PID NUE: -------------
	  if(fHCorrectMCNegChrg){
	    iBinNUA = fHCorrectMCNegChrg->FindBin(trk1Pt);  // re-use iBinNUA
	    ptWgtMCtrk1 = fHCorrectMCNegChrg->GetBinContent(iBinNUA);
	  }
	  //now PID NUE:
	  if(kPIDtrk1==1 && fHCorrectMCNegPion){
	    iBinNUA = fHCorrectMCNegPion->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCNegPion->GetBinContent(iBinNUA);
	  }
	  else if(kPIDtrk1==2 && fHCorrectMCNegKaon){
	    iBinNUA = fHCorrectMCNegKaon->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCNegKaon->GetBinContent(iBinNUA);
	  }
	  else if(kPIDtrk1==3 && fHCorrectMCNegProt){
	    iBinNUA = fHCorrectMCNegProt->FindBin(trk1Pt);
	    ptWgtMCPIDtrk1 = fHCorrectMCNegProt->GetBinContent(iBinNUA);
	  }//----------------------------------	  	  
	}
	//------------------------------------------------------------
	

	wgt1particle = WgtNUAtrk1*ptWgtMCtrk1;        ///Charge
	wgt1PIDparticle = WgtNUAPIDtrk1*ptWgtMCPIDtrk1;  ///PID


	
	//----- For Chunzheng ---------->
	/// fill the information in vector array for Lambda-X pair

	if(bAnalysLambdaPairs){

	  Int_t code = 0;
	  isItPiontrk1 = kFALSE;
	  isItKaontrk1 = kFALSE;
	  isItProttrk1 = kFALSE; 

	  isItPiontrk1 = CheckPIDofParticle(AODtrack1,1); // 1=pion
	  isItKaontrk1 = CheckPIDofParticle(AODtrack1,2); // 2=Kaon
	  isItProttrk1 = CheckPIDofParticle(AODtrack1,3); // 3=proton
	
	  if(trk1Chrg > 0) {
	    if(isItPiontrk1) code = 211;
	    if(isItKaontrk1) code = 321;
	    if(isItProttrk1) code = 2212;
	  }
	  else{  /// 
	    if(isItPiontrk1) code = -211;
	    if(isItKaontrk1) code = -321;
	    if(isItProttrk1) code = -2212;
	  }

	  
	  Int_t trk1ID = AODtrack1->GetID();//unique in a event
	  vecPDGCode.push_back(code);
	  vecPhi.push_back(trk1Phi);
	  vecEta.push_back(trk1Eta);
	  vecPt.push_back(trk1Pt);
	  vecID.push_back(trk1ID);
	  vecNUAWeight.push_back(WgtNUAtrk1);  //if Particle is considered Unidentified
	  vecNUEWeight.push_back(ptWgtMCtrk1); //if Particle is considered Unidentified
	  vecNUAWeightPID.push_back(WgtNUAPIDtrk1);  //if Particle is considered Identified (good purity upto pT<2 GeV)
	  vecNUEWeightPID.push_back(ptWgtMCPIDtrk1); //if Particle is considered Identified (good purity upto pT<2 GeV)

	  continue; /// Skip Analysing CME PID if Lambda-X study is running. Otherwise Jobs may reach TTL and killed.
	  
	}
	
	//<----- Chunzheng Filled his Vectors for single Particle id and weights.
	// Rest of the Lambda Analysis is done outside this track loop.



	
	/// Rihan: The part below is only relevant for PID CME Analysis:




	
	bPIDoktrk1=kFALSE;
	bPIDoktrk2=kFALSE;



	/// Check if track1 is any of the desired Pair:
	
	bPIDoktrk1 = CheckPIDofParticle(AODtrack1,kPIDtrk1);  //check if track1 is of PID Req#1

	if(!bPIDoktrk1)
	  bPIDoktrk2 = CheckPIDofParticle(AODtrack1,kPIDtrk2); //if track1 != PID Req#1, then check if it is PID Req#2

	if(!bPIDoktrk1 && !bPIDoktrk2)    /// the track-1 is None of the desired PID, hence skip this track. 
	  continue;


	///Decide which TPC EP to use based on eta of 1st track.
	
	Double_t fPsiNEvent = 0, fPsi3Event = 0;

	///Choose whether to use TPC or V0EP
	if(bUseV0EventPlane){
	  fPsiNEvent = fSelectedV0PsiN;
	  fPsi3Event = fSelectedV0Psi3;
	}
	else{ 
	  if(trk1Eta > 0){
	    fPsiNEvent = fPsiNTPCNeg;
	    fPsi3Event = fPsi3TPCNeg;
	  }
	  else{
	    fPsiNEvent = fPsiNTPCPos;
	    fPsi3Event = fPsi3TPCPos;
	  }
	}
	//// If I am here then I have one partner. Lets look for the other one: ////
	

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
	    if((trk2Pt <= fMaxPtCut) && (trk2Pt >= fMinPtCut) && (trk2Eta <= fMaxEtaCut) && (trk2Eta >= fMinEtaCut) && (trk2dEdx >= fdEdxMin) && (trk2TpcNC >= fTPCclustMin) && (trk2Chi2 >= fTrkChi2Min) && (trk2Chi2 <= fChi2Max) && TMath::Abs(trk2Chrg)) {


	      if(bPIDoktrk1){
		bPIDoktrk2 = CheckPIDofParticle(AODtrack2,kPIDtrk2); //
	      }
	      else if(bPIDoktrk2){
		bPIDoktrk1 = CheckPIDofParticle(AODtrack2,kPIDtrk1);
	      }

	      if(!bPIDoktrk1 && !bPIDoktrk2)    /// the track1 and track2 is None of the desired PID.
		continue;



	      //// If I am here then I have both partners. Lets fill correlator Histograms:
	
	      ptWgtMCtrk2 = 1.0;   
	      WgtNUAtrk2  = 1.0;
	      ptWgtMCPIDtrk2 = 1.0;
	      WgtNUAPIDtrk2 = 1.0;

	      /// Rihan: I should move the following part to a function =>
	      //------------ Get the NUA and MC weight for Track2-----------------
	      if(trk2Chrg>0){
		if(fHCorrectNUAChrgPos){ /// safety measures for breaks!
		  iBinNUA = fHCorrectNUAChrgPos->FindBin(pVtxZ,trk2Phi,trk2Eta);
		  WgtNUAtrk2 = fHCorrectNUAChrgPos->GetBinContent(iBinNUA);
		}
		if(fHCorrectNUAkPIDPos){
		  //iBinNUA = fHCorrectNUAChrgPos->FindBin(pVtxZ,trk2Phi,trk2Eta); //same bin 
		  WgtNUAPIDtrk2 = fHCorrectNUAkPIDPos->GetBinContent(iBinNUA);
		}

		///--------- Charge & PID NUE: -------------
		if(fHCorrectMCPosChrg){
		  iBinNUA = fHCorrectMCPosChrg->FindBin(trk2Pt);  // re-use iBinNUA
		  ptWgtMCtrk2 = fHCorrectMCPosChrg->GetBinContent(iBinNUA);
		}
		//now PID NUE:
		if(kPIDtrk2==1 && fHCorrectMCPosPion){
		  iBinNUA = fHCorrectMCPosPion->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCPosPion->GetBinContent(iBinNUA);
		}
		else if(kPIDtrk2==2 && fHCorrectMCPosKaon){
		  iBinNUA = fHCorrectMCPosKaon->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCPosKaon->GetBinContent(iBinNUA);
		}
		else if(kPIDtrk2==3 && fHCorrectMCPosProt){
		  iBinNUA = fHCorrectMCPosProt->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCPosProt->GetBinContent(iBinNUA);
		}//----------------------------------	  
		
	      }
	      else{
		if(fHCorrectNUAChrgNeg){ /// safety measures for breaks!
		  iBinNUA = fHCorrectNUAChrgNeg->FindBin(pVtxZ,trk2Phi,trk2Eta);
		  WgtNUAtrk2 = fHCorrectNUAChrgNeg->GetBinContent(iBinNUA);
		}
		if(fHCorrectNUAkPIDNeg){
		  //iBinNUA = fHCorrectNUAChrgNeg->FindBin(pVtxZ,trk2Phi,trk2Eta); //same bin 
		  WgtNUAPIDtrk2 = fHCorrectNUAkPIDNeg->GetBinContent(iBinNUA);
		}

		///--------- Charge & PID NUE: -------------
		if(fHCorrectMCNegChrg){
		  iBinNUA = fHCorrectMCNegChrg->FindBin(trk2Pt);  // re-use iBinNUA
		  ptWgtMCtrk2 = fHCorrectMCNegChrg->GetBinContent(iBinNUA);
		}
		//now PID NUE:
		if(kPIDtrk2==1 && fHCorrectMCNegPion){
		  iBinNUA = fHCorrectMCNegPion->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCNegPion->GetBinContent(iBinNUA);
		}
		else if(kPIDtrk2==2 && fHCorrectMCNegKaon){
		  iBinNUA = fHCorrectMCNegKaon->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCNegKaon->GetBinContent(iBinNUA);
		}
		else if(kPIDtrk2==3 && fHCorrectMCNegProt){
		  iBinNUA = fHCorrectMCNegProt->FindBin(trk2Pt);
		  ptWgtMCPIDtrk2 = fHCorrectMCNegProt->GetBinContent(iBinNUA);
		}//----------------------------------	 		
	      }
	      //------------------------------------------------------------
	


	      
	      wgt2particle   = WgtNUAtrk1*ptWgtMCtrk1 * WgtNUAtrk2*ptWgtMCtrk2;     //Combined weight.. 
	      wgt2PIDparticle = WgtNUAPIDtrk1*ptWgtMCPIDtrk1 * WgtNUAPIDtrk2*ptWgtMCPIDtrk2;  ///PID



	      
	      if(trk1Chrg*trk2Chrg < 0){ //Opposite sign	
		hAvg3pC112vsCentOS->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgt2particle);
		hAvg3pC123vsCentOS->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgt2particle);
				
		hAvgDelta1vsCentOS->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgt2particle);
		hAvgDelta2vsCentOS->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta3vsCentOS->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta4vsCentOS->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgt2particle);		
	      }		
	      else if(trk1Chrg > 0 && trk2Chrg > 0){		      
		hAvg3pC112vsCentPP->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgt2particle);
		hAvg3pC123vsCentPP->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgt2particle);
				
		hAvgDelta1vsCentPP->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgt2particle);
		hAvgDelta2vsCentPP->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta3vsCentPP->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta4vsCentPP->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgt2particle);		
	      }
	      //else if(trk1Chrg < 0 && trk2Chrg < 0){  ///this is obvious!
	      else{
		hAvg3pC112vsCentNN->Fill(centrality,TMath::Cos(trk1Phi +  trk2Phi  - 2*fPsiNEvent),wgt2particle);
		hAvg3pC123vsCentNN->Fill(centrality,TMath::Cos(trk1Phi + 2*trk2Phi - 3*fPsi3Event),wgt2particle);
		
		hAvgDelta1vsCentNN->Fill(centrality,TMath::Cos(trk1Phi - trk2Phi),wgt2particle);
		hAvgDelta2vsCentNN->Fill(centrality,TMath::Cos(2.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta3vsCentNN->Fill(centrality,TMath::Cos(3.*(trk1Phi - trk2Phi)),wgt2particle);
		hAvgDelta4vsCentNN->Fill(centrality,TMath::Cos(4.*(trk1Phi - trk2Phi)),wgt2particle);		
	      }
	      
	      
	    }//j-track cuts
	  }//j-track FB validated
	}///j-track loop



      	

      }//----> i-track loop => All trackCuts applied.     
    }//-----> i-track loop => FB is validated.    
  }///-----> i-track loop Ends here.<--------
 









  

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
    vector<Double_t> vecLambdaPhi;
    vector<Int_t>  vecLambdaPosID; // Daughter ID
    vector<Int_t>  vecLambdaNegID; // Daughter ID

    vector<Double_t> vecAntiLambdaPhi;
    vector<Int_t>  vecAntiLambdaPosID; // Daughter ID
    vector<Int_t>  vecAntiLambdaNegID; // Daughter ID

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

      //daughers Nsigma: Rihan-> There could be improvement: which case to choose!!
      Float_t nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton));//TPC p+
      Float_t nSigTOFPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack,AliPID::kProton));//TOF p+
      Float_t nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion));//TPC π-
      Float_t nSigTOFNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack,AliPID::kPion));//TOF π-

      Float_t nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion));//TPC π+
      Float_t nSigTOFPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack,AliPID::kPion));//TOF π+
      Float_t nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton));//TPC p-
      Float_t nSigTOFNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack,AliPID::kProton));//TOF p-

      TVector2 Vt(v0->MomV0X(), v0->MomV0Y());
      Double_t phi = Vt.Phi();
      Int_t id_posDaughter = v0->GetPosID();
      Int_t id_negDaughter = v0->GetNegID();

      Double_t massLambda     = v0->MassLambda();
      Double_t massAntiLambda = v0->MassAntiLambda();

      //Float_t fDautPosEta = pTrack->Eta();
      //Float_t fDautNegEta = nTrack->Eta();
      
      //Lambda PID 
      //Λ-->(p+)+(π-)
      if(nSigTPCPosProton < fDaughtersNsigma && nSigTOFPosProton < fDaughtersNsigma && nSigTPCNegPion < fDaughtersNsigma && nSigTOFNegPion < fDaughtersNsigma)  
	{
	  fHistLambdaPt[0]->Fill(pt);
	  fHistLambdaEta[0]->Fill(eta);
	  fHistLambdaDcaToPrimVertex[0]->Fill(dcaToPV);
	  fHistLambdaCPA[0]->Fill(CPA);
	  fHistLambdaDecayLength[0]->Fill(dl);
	  fHistLambdaMass[0]->Fill(massLambda);
	  fProfileLambdaMassVsPt[0]->Fill(pt,massLambda);

	  if(TMath::Abs(massLambda - fMassMean) < fLambdaMassCut)
	    {
	      fHistLambdaPt[1]->Fill(pt);
	      fHistLambdaEta[1]->Fill(eta);
	      fHistLambdaDcaToPrimVertex[1]->Fill(dcaToPV);
	      fHistLambdaCPA[1]->Fill(CPA);
	      fHistLambdaDecayLength[1]->Fill(dl);
	      fHistLambdaMass[1]->Fill(massLambda);
	      fProfileLambdaMassVsPt[1]->Fill(pt,massLambda);

	      vecLambdaPhi.push_back(phi);
	      vecLambdaPosID.push_back(id_posDaughter);
	      vecLambdaNegID.push_back(id_negDaughter);
	    }
	}

      //AntiLambda PID
      //(-Λ)-->(p-)+(π+)
      if(nSigTPCNegProton < fDaughtersNsigma && nSigTOFNegProton < fDaughtersNsigma && nSigTPCPosPion < fDaughtersNsigma && nSigTOFPosPion < fDaughtersNsigma) 
	{
	  fHistAntiLambdaPt[0]->Fill(pt);
	  fHistAntiLambdaEta[0]->Fill(eta);
	  fHistAntiLambdaDcaToPrimVertex[0]->Fill(dcaToPV);
	  fHistAntiLambdaCPA[0]->Fill(CPA);
	  fHistAntiLambdaDecayLength[0]->Fill(dl);
	  fHistAntiLambdaMass[0]->Fill(massAntiLambda);
	  fProfileAntiLambdaMassVsPt[0]->Fill(pt,massAntiLambda);

	  if(TMath::Abs(massAntiLambda - fMassMean) < fLambdaMassCut)
	    {
	      fHistAntiLambdaEta[1]->Fill(eta);
	      fHistAntiLambdaPt[1]->Fill(pt);
	      fHistAntiLambdaDcaToPrimVertex[1]->Fill(dcaToPV);
	      fHistAntiLambdaCPA[1]->Fill(CPA);
	      fHistAntiLambdaDecayLength[1]->Fill(dl);
	      fHistAntiLambdaMass[1]->Fill(massAntiLambda);
	      fProfileAntiLambdaMassVsPt[1]->Fill(pt,massAntiLambda);

	      vecAntiLambdaPhi.push_back(phi);
	      vecAntiLambdaPosID.push_back(id_posDaughter);
	      vecAntiLambdaNegID.push_back(id_negDaughter);
	    }
	}
    
    }

    //// Now fill Lambda-X correlations:
    for (vector<Double_t>::size_type iTrk = 0; iTrk < vecPt.size(); iTrk++) {
      Int_t    id_1   = vecID[iTrk];
      Int_t    code_1 = vecPDGCode[iTrk];
      Double_t pt_1   = vecPt[iTrk];
      Double_t eta_1  = vecEta[iTrk];
      Double_t phi_1  = vecPhi[iTrk];



      //Lambda - X
      for (vector<double>::size_type jLambda = 0; jLambda < vecLambdaPhi.size(); jLambda++)
	{
	  Double_t phi_lambda    =  vecLambdaPhi[jLambda];
	  Short_t  id_posDaughter = vecLambdaPosID[jLambda];
	  Short_t  id_negDaughter = vecLambdaNegID[jLambda];
	  if(id_1 == id_posDaughter || id_1 == id_negDaughter) continue;  // checking if charged particle is daughter of Lambda itself..


	 
	  //double qx = TMath::Cos(2*phi_1);   //<--- Not needed! as we use EP from opposite eta of track!
	  //double qy = TMath::Sin(2*phi_1);
	  
	 
	  Double_t fPsiNNoAuto = 0.;
	  
	  Double_t fTPCQxTemp = 0, fTPCQyTemp = 0; /// Get the Total sum of Qx,Qy  locally, then Remove AutoCorr if needed.
	  Double_t qx=0, qy=0;

	  
	  if(bUseV0EventPlane){ /// If we want to use V0 Event Plane, no Auto,short-range Correlation.
	    fPsiNNoAuto = fSelectedV0PsiN;	 
	  }
	  else{ //Use TPC EP and Remove AutoCorrelation:
	    
	    if(eta_1 > 0) { // use EP from opposite eta than the charged track! One way to remove AutoCorrelation.
	      
	      fTPCQxTemp = SumQnxTPCNeg;   
	      fTPCQyTemp = SumQnyTPCNeg;
	    
	      if(find(vecNegEPTrkID.begin(),vecNegEPTrkID.end(), id_posDaughter) != vecNegEPTrkID.end()){

		vector<int>::iterator iter = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_posDaughter);
		if (iter != vecNegEPTrkID.end()){
		  int iPosDaughter = distance(vecNegEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iPosDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iPosDaughter]);
		}
	      }
	    
	      if(find(vecNegEPTrkID.begin(),vecNegEPTrkID.end(), id_negDaughter) != vecNegEPTrkID.end()){
	      
		vector<int>::iterator iter = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_negDaughter);
		if (iter != vecNegEPTrkID.end()){
		  int iNegDaughter = distance(vecNegEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iNegDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iNegDaughter]);
		}
	      }
	    }///for -ve EP
	    else{
	    
	      fTPCQxTemp = SumQnxTPCPos;   
	      fTPCQyTemp = SumQnyTPCPos;

	      if(find(vecPosEPTrkID.begin(),vecPosEPTrkID.end(), id_posDaughter) != vecPosEPTrkID.end()){

		vector<int>::iterator iter = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_posDaughter);
		if (iter != vecPosEPTrkID.end()){
		  int iPosDaughter = distance(vecPosEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iPosDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iPosDaughter]);
		}
	      }
	    
	      if(find(vecPosEPTrkID.begin(),vecPosEPTrkID.end(), id_negDaughter) != vecPosEPTrkID.end()){
	      
		vector<int>::iterator iter = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_negDaughter);
		if (iter != vecPosEPTrkID.end()){
		  int iNegDaughter = distance(vecPosEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iNegDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iNegDaughter]);
		}
	      }
	    }/// for PosEP 

	  
	    fTPCQxTemp -= qx;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used.. 
	    fTPCQyTemp -= qy;   	      
	    fPsiNNoAuto = (1./gPsiN)*TMath::ATan2(fTPCQyTemp,fTPCQxTemp);   //AutoCorrelation Removed EP.
	    if(fPsiNNoAuto < 0) fPsiNNoAuto += TMath::TwoPi()/gPsiN;  	  

	  }///case if TPC Event Plane is used...

	  

	        
	  Double_t delta = TMath::Cos(phi_lambda - phi_1);
	  Double_t gammaTPC  = TMath::Cos(phi_lambda + phi_1 - 2 *fPsiNNoAuto);


	  if(code_1 > 0){
	    fProfileGammaTPC_Lambda_hPos->Fill(centrality,gammaTPC);

	    if(code_1==2212) fProfileGammaTPC_Lambda_Proton->Fill(centrality,gammaTPC);
	  }
	  else{
	    fProfileGammaTPC_Lambda_hNeg->Fill(centrality,gammaTPC);
	    if(code_1==-2212) fProfileGammaTPC_Lambda_AntiProton->Fill(centrality,gammaTPC);
	  }
      
	}// Lambda-X pair done



      for (vector<double>::size_type jAntiLambda = 0; jAntiLambda < vecAntiLambdaPhi.size(); jAntiLambda++)
	{
	  Double_t phi_antiLambda = vecAntiLambdaPhi[jAntiLambda];
	  Short_t  id_posDaughter = vecAntiLambdaPosID[jAntiLambda];
	  Short_t  id_negDaughter = vecAntiLambdaNegID[jAntiLambda];
	  if(id_1 == id_posDaughter || id_1 == id_negDaughter) continue;


	 
	  Double_t fPsiNNoAuto = 0.;
	  Double_t fTPCQxTemp = 0, fTPCQyTemp = 0; ///Get the Total sum of Qx,Qy  locally, then Remove AutoCorr if needed.
	  Double_t qx=0, qy=0;
	  
	  if(bUseV0EventPlane){ /// If we want to use V0 Event Plane, no Auto,short-range Correlation.
	    fPsiNNoAuto = fSelectedV0PsiN;	 
	  }
	  else{ //Use TPC EP and Remove AutoCorrelation:
	    
	    if(eta_1 > 0){ // use EP from opposite eta than the charged track! One way to remove AutoCorrelation.

	      fTPCQxTemp = SumQnxTPCNeg;   
	      fTPCQyTemp = SumQnyTPCNeg;
	    
	      if(find(vecNegEPTrkID.begin(),vecNegEPTrkID.end(), id_posDaughter) != vecNegEPTrkID.end()){

		vector<int>::iterator iter = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_posDaughter);
		if (iter != vecNegEPTrkID.end()){
		  int iPosDaughter = distance(vecNegEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iPosDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iPosDaughter]);
		}
	      }
	    
	      if(find(vecNegEPTrkID.begin(),vecNegEPTrkID.end(), id_negDaughter) != vecNegEPTrkID.end()){
	      
		vector<int>::iterator iter = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_negDaughter);
		if (iter != vecNegEPTrkID.end()){
		  int iNegDaughter = distance(vecNegEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iNegDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iNegDaughter]);
		}
	      }
	    }///for -ve EP
	    else{
	    
	      fTPCQxTemp = SumQnxTPCPos;   
	      fTPCQyTemp = SumQnyTPCPos;

	      if(find(vecPosEPTrkID.begin(),vecPosEPTrkID.end(), id_posDaughter) != vecPosEPTrkID.end()){

		vector<int>::iterator iter = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_posDaughter);
		if (iter != vecPosEPTrkID.end()){
		  int iPosDaughter = distance(vecPosEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iPosDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iPosDaughter]);
		}
	      }
	    
	      if(find(vecPosEPTrkID.begin(),vecPosEPTrkID.end(), id_negDaughter) != vecPosEPTrkID.end()){
	      
		vector<int>::iterator iter = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_negDaughter);
		if (iter != vecPosEPTrkID.end()){
		  int iNegDaughter = distance(vecPosEPTrkID.begin(), iter);
		  qx += TMath::Cos(2*vecPhi[iNegDaughter]);   
		  qy += TMath::Sin(2*vecPhi[iNegDaughter]);
		}
	      }
	    }/// for PosEP 

	  
	    fTPCQxTemp -= qx;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used.. 
	    fTPCQyTemp -= qy;

	    fPsiNNoAuto = (1./gPsiN)*TMath::ATan2(fTPCQyTemp,fTPCQxTemp);
	    if(fPsiNNoAuto < 0) fPsiNNoAuto += TMath::TwoPi()/gPsiN; 
	  }


	  Double_t delta = TMath::Cos(phi_antiLambda - phi_1);
	  Double_t gammaTPC = TMath::Cos(phi_antiLambda + phi_1 - 2*fPsiNNoAuto);

	  //antiLambda - h+
	  if(code_1 > 0){
	    fProfileGammaTPC_AntiLambda_hPos->Fill(centrality,gammaTPC);
	    if(code_1==2212) fProfileGammaTPC_AntiLambda_Proton->Fill(centrality,gammaTPC);	
	  }
	  else{       
	    fProfileGammaTPC_AntiLambda_hNeg->Fill(centrality,gammaTPC);
	    if(code_1==-2212) fProfileGammaTPC_AntiLambda_AntiProton->Fill(centrality,gammaTPC);
	  }
     
	}//anti-Lambda-X done
    
  
    }/// loop over charge particle array

  }////--------- Lambda Histograms are Filled ----------
  


  



  
  //Last lines in Event loop
  fHistVertexZcm->Fill(pVtxZ);
  fHistVxvsVzMinBias->Fill(pVtxZ,pVtxX);
  fHistVyvsVzMinBias->Fill(pVtxZ,pVtxY);
  fCentDistAfterCut->Fill(centrality);  
  //Post the Histograms:  
  PostData(1,fListHist);

  // std::cout<<" Info:UserExec()  Call Finished ..!!!\n";
}//---------------- UserExec ----------------------



Bool_t AliAnalysisTaskGammaDeltaPID::IsGoodV0(AliAODv0 *aodV0) 
{
  if (!aodV0) {
    AliError(Form("ERROR: Could not retrieve aodV0"));
    return kFALSE;
  }
  //* Offline reconstructed V0 only
  if ( aodV0->GetOnFlyStatus() ) return kFALSE;
  //* Get daughters and check them     
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return kFALSE;
  }
  //* Unlike signs of daughters
  if ( myTrackNegTest->Charge() == myTrackPosTest->Charge() ) return kFALSE;
  //* Cosinus of pointing angle      
  double dCPA = aodV0->CosPointingAngle(fCurrentVtx);
  //* cut on Cosinus of pointing angle
  if ( dCPA < fV0CPAMin ) return kFALSE;
  //* DCA of V0
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if ( TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax ) return kFALSE;
  //* V0 parh length before decay
  double dDecayLength = aodV0->DecayLengthV0(fCurrentVtx);
  if ( dDecayLength > fV0DecayLengthMax ) return kFALSE;
  if ( dDecayLength < fV0DecayLengthMin ) return kFALSE;
  //* DCA between daughters
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
  //* TPC refit
  if ( !track->IsOn(AliAODTrack::kTPCrefit) ) return kFALSE;
  //* Maximum value of transverse momentum 
  double dPt = track->Pt();
  if (dPt > fDaughtersPtMax) return kFALSE;
  //* Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (TMath::Abs(dEta) > fDaughtersEtaMax) return kFALSE;
  //* Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return kFALSE;
  //* Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  //* [number of crossed rows]>0.8 * [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  return kTRUE;
}














Bool_t AliAnalysisTaskGammaDeltaPID::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){

  if(pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check
  
  Bool_t bPIDokay = kFALSE;

  if(!fPIDResponse){
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cut from AddTaskMacro for TPC only pid! Although someone barely needs to change it.
  
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
  else if(pidToCheck==4){///

  }
  else{
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3,4  (Charge Pion, Kaon, Proton, Lambda)\n return with kFALSE \n");
    return kFALSE;
  }

  return kFALSE;
}



Bool_t AliAnalysisTaskGammaDeltaPID::CheckEventIsPileUp2018(AliAODEvent *faod) {


  Bool_t BisPileup=kFALSE;

  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
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


  if (centrCL0 < fCenCutLowPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 1st**************************:"<<endl;
      BisPileup=kTRUE;
    }

  if (centrCL0 > fCenCutHighPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 2nd**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls))
    {
      //cout<<"*****************hi i am in 3rd**************************:"<<endl;
      BisPileup=kTRUE;
    }
        

  if (multV0On < fV0CutPU->Eval(multV0Tot))
    {
      //cout<<"*****************hi i am in 4th**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 5th**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0)
    BisPileup=kTRUE;


  if (faod->IsIncompleteDAQ())
    BisPileup=kTRUE;
        
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);
  fHistTPCVsESDTrkBefore->Fill(multTrk,multEsd);  
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTrk);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);

  if (!BisPileup)
    {      
      fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);
      fHistTPCVsESDTrkAfter->Fill(multTrk,multEsd);  
      fHistTPConlyVsCL1After->Fill(centrCL1,multTrk);
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);      
    }

 return BisPileup; 

}














Bool_t AliAnalysisTaskGammaDeltaPID::CheckEventIsPileUp(AliAODEvent *faod) {

  //cout<<"---PileUp cut for LHC15o dataset--"<<endl;

  Bool_t BisPileup=kFALSE;


  Double_t centrV0M=300;
  Double_t centrCL1=300;
  Double_t centrCL0=300;
  Double_t centrTRK=300;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **Error** ::UserExec() AliMultSelection object not found. Quit Job!! \n\n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
  centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");


  //-- pile-up a la Dobrin for LHC15o -----
  if(PileUpMultiVertex(faod)) {
    fHistPileUpCount->Fill(0.5);    
    BisPileup=kTRUE;
  }
  Int_t isPileup = faod->IsPileupFromSPD(3);
  if(isPileup != 0) {
    fHistPileUpCount->Fill(1.5);
    BisPileup=kTRUE;          
  }
  if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    fHistPileUpCount->Fill(2.5);
    BisPileup=kTRUE;
  }
  if(faod->IsIncompleteDAQ())  {
    fHistPileUpCount->Fill(3.5);
    BisPileup=kTRUE;
  }
  if(fabs(centrV0M-centrCL1)> 5.0)  {//
    fHistPileUpCount->Fill(4.5);
    BisPileup=kTRUE;
  }



  // check vertex consistency
  const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
  const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

  if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
    fHistPileUpCount->Fill(5.5);
    BisPileup=kTRUE;
  }

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);

  double dz = vtTrc->GetZ() - vtSPD->GetZ();

  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;

  if(TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
    fHistPileUpCount->Fill(6.5);
    BisPileup=kTRUE;
  }




  Int_t multTPC = 0;
  Int_t multTPCAll = 0;
  Int_t multITSfb96 = 0;
  Int_t multITSfb32 = 0;

  Int_t multTPCFE = 0;
  Int_t multGlobal = 0;
  Int_t multTPCuncut = 0;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  const Int_t nTracks = faod->GetNumberOfTracks();

  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    
    AliAODTrack* track = (AliAODTrack*)faod->GetTrack(iTracks);
    if(!track)  continue;
  
    //---------- old method -----------
    if(track->TestFilterBit(128))
      multTPCAll++;

    //----------------------------------
    if(track->TestFilterBit(1))  multTPCuncut++;
    if(track->TestFilterBit(32)) multITSfb32++;


    if(track->Pt()<0.2 || track->Pt()>10.0 || TMath::Abs(track->Eta())>0.8 || track->GetTPCNcls()<fTPCclustMin || track->GetTPCsignal()<10.0)
      continue;
    if(track->GetDetPid() && track->Chi2perNDF() > 0.2) multTPC++;
    if(track->TestFilterBit(1) && track->Chi2perNDF()>0.2)  multTPCFE++;
    if(!track->TestFilterBit(16) || track->Chi2perNDF()<0.1)   continue;
                
    Double_t b[2]    = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
                
    AliAODTrack copy(*track);
    Double_t magField = faod->GetMagneticField();
                
    if(magField!=0){     
      if(track->PropagateToDCA(faod->GetPrimaryVertex(), magField, 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) 
	multGlobal++;    
    }
  }


  Double_t multESDTPCDif  = multEsd  - fPileUpSlopeParm*multTPCAll;
  
  Bool_t  bIsOutLier=kFALSE;
  if(multTPC < (-20.0+1.15*multGlobal) || multTPC > (200.+1.45*multGlobal)) { bIsOutLier = kTRUE;}

  
  fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);
  fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);  
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTPCAll);
  

  if(multESDTPCDif > fPileUpConstParm) { 
    fHistPileUpCount->Fill(7.5);
    BisPileup=kTRUE;
  }
  if(BisPileup==kFALSE) {
    if(!fMultSelection->GetThisEventIsNotPileup())  BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
    if(BisPileup)     fHistPileUpCount->Fill(9.5);
  }  
  if(!BisPileup) {
    fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);  
    fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  
    fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
    fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
  }
  


  return BisPileup; 
} //-------pile up function ------









Bool_t AliAnalysisTaskGammaDeltaPID::PileUpMultiVertex(const AliAODEvent* faod)
{  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=faod->GetNumberOfPileupVerticesTracks()))
    return kFALSE;

  vtPrm = faod->GetPrimaryVertex();
  if(vtPrm == faod->GetPrimaryVertexSPD())
    return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)faod->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return kTRUE; // pile-up: well separated vertices
  }
  return kFALSE;
}




double AliAnalysisTaskGammaDeltaPID::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    AliDebug(2,"\n\n ::GetWDist => One of vertices is not valid\n\n");
    return 0;
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
  if (!vVb.IsValid()) {
    AliDebug(2,"Singular Matrix\n");
    return dist;
  }
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
    +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;

}











void AliAnalysisTaskGammaDeltaPID::SetupEventAndTaskConfigInfo(){
  
  fHistPileUpCount = new TH1F("fHistPileUpCount","PileUp Counts (LHC15o only)", 15, 0., 15.);
  fHistPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fHistPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fHistPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultComb08");
  fHistPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fHistPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>5.0");
  fHistPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fHistPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  Int_t puConst = fPileUpConstParm;
  fHistPileUpCount->GetXaxis()->SetBinLabel(8,Form("multESDTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(9,Form("multGlobTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(10,"PileUpMultSelTask");
  fListHist->Add(fHistPileUpCount);



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
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(9,"V0 Gain/Q App");
  if(fListV0MCorr){
    fHistAnalysisInfo->SetBinContent(9,1);
  }
  

  fHistAnalysisInfo->GetXaxis()->SetBinLabel(10,"ZDC Gain/Q App"); 
  //if(fListV0MCorr){
  //fHistAnalysisInfo->SetBinContent(10,1);
  //}
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(11,"V0CEP");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(12,"V0AEP");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(13,"TPCEP");
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(14,"Reserved");


  fHistAnalysisInfo->GetXaxis()->SetBinLabel(15,"FilterBit");
  fHistAnalysisInfo->SetBinContent(15,fFilterBit);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(16,"N_{cls,TPC}");
  fHistAnalysisInfo->SetBinContent(16,fTPCclustMin);
  fHistAnalysisInfo->GetXaxis()->SetBinLabel(17,"Chi^{2}max");
  fHistAnalysisInfo->SetBinContent(17,fChi2Max);
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
}//----------SetupEventAndTaskConfigInfo-----------




Int_t AliAnalysisTaskGammaDeltaPID::GetCentralityScaled0to10(Double_t fCent){

 Int_t cIndex = 0;

 if(fCent<5.0) {
   cIndex  = 0; 
 }
 else if(fCent>=5.0 && fCent<10){
   cIndex  = 1;
 }
 else if(fCent>=10.0) {
   cIndex = abs(fCent/10.0)+1;
 }
 return cIndex;
 
}//------------GetCentralityScaled0to10------------







void AliAnalysisTaskGammaDeltaPID::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{
  if(fListNUACorr){

    fHCorrectNUAChrgPos = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,run));  // 0 = Charge
    fHCorrectNUAChrgNeg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,run));  // 0 = Charge
    ///PID pion=1,Kaon=2,Prot=3
    fHCorrectNUAkPIDPos = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",kParticleID,run)); 
    fHCorrectNUAkPIDNeg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",kParticleID,run)); 
    
    //if(fHCorrectNUAChrgPos && fHCorrectNUAChrgNeg){
    //cout<<"\n=========== Info:: Setting up NUA corrections for run "<<run<<"============"<<endl;   
    //}
    /// Now Get Average Q vector:
    fHCorrectTPCQNxEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCosNPsivsCentEtaPosRun%d",run));
    fHCorrectTPCQNyEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSinNPsivsCentEtaPosRun%d",run));
    fHCorrectTPCQNxEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCosNPsivsCentEtaNegRun%d",run));
    fHCorrectTPCQNyEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSinNPsivsCentEtaNegRun%d",run));

    fHCorrectTPCQ3xEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCos3PsivsCentEtaPosRun%d",run));
    fHCorrectTPCQ3yEtaPos = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSin3PsivsCentEtaPosRun%d",run));
    fHCorrectTPCQ3xEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgCos3PsivsCentEtaNegRun%d",run));
    fHCorrectTPCQ3yEtaNeg = (TH1D *) fListNUACorr->FindObject(Form("fHisAvgSin3PsivsCentEtaNegRun%d",run));
             
    //if(fHCorrectNUAChrgPos && fHCorrectNUAChrgNeg){
    //cout<<"\n=========== Info:: Found TPC <Q> vectors for run "<<run<<"============"<<endl;   
    //}

  }
  else {
    //printf("\n ******** Error/Warning: No NUA Correction found for run %d, Use PID Wgt = 1.0 ********* \n",run);
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
}

void AliAnalysisTaskGammaDeltaPID::GetV0MCorrectionHist(Int_t run){ 

  if(fListV0MCorr){
    
    fHCorrectV0ChWeghts = (TH2F *) fListV0MCorr->FindObject(Form("hWgtV0ChannelsvsVzRun%d",run));

    if(fHCorrectV0ChWeghts){ // Load <Q> vector if Gain Correction is available, otherwise no Point.
      //printf("\n ::Info() V0 Channel Weights Found for Run %d ",run);
      fHCorrectQNxV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0CRun%d",run));
      fHCorrectQNyV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0CRun%d",run));    
      fHCorrectQNxV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0ARun%d",run));
      fHCorrectQNyV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0ARun%d",run));
	
      fHCorrectQ3xV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0CRun%d",run));
      fHCorrectQ3yV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0CRun%d",run));    
      fHCorrectQ3xV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0ARun%d",run));
      fHCorrectQ3yV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0ARun%d",run));
    }
  }
  else{
    fHCorrectV0ChWeghts=NULL;
  } 
}











