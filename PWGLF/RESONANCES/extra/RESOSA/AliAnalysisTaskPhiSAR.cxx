/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
( * Author: The ALICE Off-line Project.                                    *
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
/* AliAnalysisTaskPhiSAR.cxx*/
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TList.h"
#include "TRandom3.h"
#include "TVectorD.h"
#include "TObjArray.h"
#include "TH2.h"
#include "TFile.h"
#include "TPRegexp.h"
#include "TChain.h"
#include "TF1.h"
#include "TSpline.h"
#include "TRandom2.h"

#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
//#include "AliFlowTrackSimple.h"
//#include "AliFlowTrackCuts.h"

//#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include <AliVParticle.h>
#include <AliVTrack.h>
#include "AliVVertex.h"
#include <AliLog.h>
#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliITSPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <AliTRDPIDResponse.h>
#include <AliTOFPIDResponse.h>
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"
#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPhiSAR.h"

//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
//#include "AliTimeRangeCut.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPhiSAR)

//____________________________________________________//
AliAnalysisTaskPhiSAR::AliAnalysisTaskPhiSAR() // All data members should be initialised here
:AliAnalysisTaskSE(0),
  fOutput(0),
  fAliESDtrackCuts(0x0),
  fESDtrackCutsForRP(0x0),
  fPIDResponse(0x0),
  fMultSelection(NULL),
  fFlowQnVectorMgr(0x0),
  fTriggerAna(0x0),
  fAnalysisLevel("ESD"),
  fTrackType("GLOBAL"),
  fRPTrackType("TPC"),
  fSytCutType("PtDCAXY7s"),
  fSystematicCuts(kFALSE),
  fShiftList(NULL),
  fUseShift(kFALSE),
  fAcceptTPC(kFALSE),
  fCheckPileUp(kFALSE),
  fIgnoreTPCzRange(kFALSE),
  fIgnoreTPCzRangeMax(FLT_MAX),
  fIgnoreTPCzRangeMin(-FLT_MAX),
  fCutMinimalTPCdedx(kFALSE),
  fMinimalTPCdedx(0.),
  fUsercuts(kFALSE),
  kaonMass(0.49368),
  pionMass(0.13957),
  kstarMass(0.89594),
  phiMass(1.01946),
  fCurrentEventCentralityBin(-1),
  fCurrentEventDoughterOne(0),
  fCurrentEventDoughterTwo(0),
  fCurrentEventVx(0),
  fCurrentEventVy(0),
  fCurrentEventVz(0),
  fBufferPointer(0),
  fQVector(0),
  fQVZeroA(0),
  fQVZeroC(0),
  fPairRapidity(1.0),
  fDeltaAngle(0),
  fHistZVertex(0), 
  fHistCentralityEvtCount(0),
  fHistEventCount(0),
  fHistCentDist(0),
  fHistPionVsKaonMult(0),  
  fHistTPCClusterRP(0),
  fHistITSClusterRP(0),
  fHistChiSqrPerNdfTPCRP(0),
  fHistChiSqrPerNdfITSRP(0),
  fHistdEdxKaonTpc(0),
  fHistdEdxKaonTof(0),
  fHistdEdxPionTpc(0),
  fHistdEdxPionTof(0),
  fHistTPCClusterPOI(0),
  fHistITSClusterPOI(0),
  fHistChiSqrPerNdfTPCPOI(0),
  fHistChiSqrPerNdfITSPOI(0),
  fHistnCrossRowsPOI(0),
  fHistratioCrossedRowsOverFindableClustersTPC(0),
  fHistDCAxyPOI(0),
  fHistDCAzPOI(0),
  fProCosResThreeSubEventPsiAB(0),
  fProCosResThreeSubEventPsiAC(0),
  fProCosResThreeSubEventPsiBC(0),
  fShiftCosTerm_V0A(0),
  fShiftSinTerm_V0A(0),
  fShiftCosTerm_V0C(0),
  fShiftSinTerm_V0C(0),
  fFullCosTerm_V0A(0),  
  fFullSinTerm_V0A(0),
  fFullCosTerm_V0C(0),  
  fFullSinTerm_V0C(0)
{
  // Dummy constructor ALWAYS needed for I/O.
  kaonMass = 0.49368;
  pionMass = 0.13957;
  kstarMass = 0.89594;
  phiMass  = 1.01946;

  fhInvMassSAEPvzeroA= 0;
  fhInvMassLikePPSAEPvzeroA = 0;
  fhInvMassLikeMMSAEPvzeroA = 0;
  fhInvMassMixSAEPvzeroA = 0;
  fhInvMassSAEPvzeroC = 0;
  fhInvMassMixSAEPvzeroC = 0;



  for(Int_t i = 0; i < kDimBuf; i++){
    if(i < kCenBin)
      {
	/*     	fhInvMassSAEPvzeroA[i] = 0;
	fhInvMassLikePPSAEPvzeroA[i] = 0;
	fhInvMassLikeMMSAEPvzeroA[i] = 0;
	fhInvMassMixSAEPvzeroA[i] = 0;
	fhInvMassSAEPvzeroC[i] = 0;
	//	fhInvMassLikePPSAEPvzeroC[i] = 0;
	//	fhInvMassLikeMMSAEPvzeroC[i] = 0;
	fhInvMassMixSAEPvzeroC[i] = 0;
*/       

	fhInvMassv2[i]=0;
	fhInvMassSinv2[i]=0;
	fhNumInvMassvsPtPhi[i]=0;
	fhDenInvMassvsPtPhi[i]=0;
	fhDenInvMassvsPtPhi_likepp[i]=0;
	fhDenInvMassvsPtPhi_likemm[i]=0;

	fHistEPTPC[i] = 0;
	fHistEPV0A[i] = 0;
	fHistEPV0C[i] = 0;
	fHistEPFMDA[i] = 0;
	fHistEPFMDC[i] = 0;
	fHistEPZDCA[i] = 0;
	fHistEPZDCC[i] = 0;
	fHistEPT0A[i] = 0;
	fHistEPT0C[i] = 0;

	fHistEPV0A_rw[i] = 0;
	fHistEPV0C_rw[i] = 0;
	fHistEPV0A_plain[i] = 0;
	fHistEPV0C_plain[i] = 0;
	fHistEPV0A_rec[i] = 0;
	fHistEPV0C_rec[i] = 0;
	fHistEPV0A_recShift[i] = 0;
	fHistEPV0C_recShift[i] = 0;    
	fHistEPV0A_al[i] = 0;
	fHistEPV0C_al[i] = 0;


	fHistEPV0A_tw[i] = 0;
	fHistEPV0C_tw[i] = 0;
	fHistEPV0A_sc[i] = 0;
	fHistEPV0C_sc[i] = 0;

      } 
    if(i < 100){
      fBufferEventNEvents[i] = 0.0;
      fBufferEventFull[i] = 0.0;
    }
    
    fBufferEventDoughterOne[i] = 0;
    fBufferEventDoughterTwo[i] = 0;
    fBufferEventPsi[i] = 0;
        
    for(Int_t j = 0; j < kDim; j++){
      
      fBufferEventDoughterOneCharge[i][j] = 0;
      fBufferEventDoughterOneId[i][j] = 0;
      fBufferEventDoughterOneIn[i][j] = 0;
      fBufferEventDoughterOnePx[i][j] = 0;
      fBufferEventDoughterOnePy[i][j] = 0;
      fBufferEventDoughterOnePz[i][j] = 0;
      
      fBufferEventDoughterTwoCharge[i][j] = 0;
      fBufferEventDoughterTwoId[i][j] = 0;
      fBufferEventDoughterTwoIn[i][j] = 0;
      fBufferEventDoughterTwoPx[i][j] = 0;
      fBufferEventDoughterTwoPy[i][j] = 0;
      fBufferEventDoughterTwoPz[i][j] = 0;
      
      if(i == 0){
	fCurrentEventDoughterOneCharge[j] = 0;
	fCurrentEventDoughterOneId[j] = 0;
	fCurrentEventDoughterOneTrackNumber[j] = 0;
	fCurrentEventDoughterOneIn[j] = 0;
	fCurrentEventDoughterOnePx[j] = 0;
	fCurrentEventDoughterOnePy[j] = 0;
	fCurrentEventDoughterOnePz[j] = 0;
	fCurrentEventDoughterOneQx[j] = 0;
	fCurrentEventDoughterOneQy[j] = 0;
	
	fCurrentEventDoughterTwoCharge[j] = 0;
	fCurrentEventDoughterTwoIn[j] = 0;
	fCurrentEventDoughterTwoId[j] = 0;
	fCurrentEventDoughterTwoTrackNumber[j] = 0;
	fCurrentEventDoughterTwoPx[j] = 0;
	fCurrentEventDoughterTwoPy[j] = 0;
	fCurrentEventDoughterTwoPz[j] = 0;
	fCurrentEventDoughterTwoQx[j] = 0;
	fCurrentEventDoughterTwoQy[j] = 0;
      }
    }
  }
}
//________________________________________________________________________
AliAnalysisTaskPhiSAR::AliAnalysisTaskPhiSAR(const char *name, const Bool_t useshift)//All data members should be initialised here
  :AliAnalysisTaskSE(name),
   fOutput(0),
   fAliESDtrackCuts(0x0),
   fESDtrackCutsForRP(0x0),
   fPIDResponse(0x0),
   fMultSelection(NULL),
   fFlowQnVectorMgr(0x0),
   fTriggerAna(0x0),
   fAnalysisLevel("ESD"),
   fTrackType("GLOBAL"),
   fRPTrackType("TPC"),
   fSytCutType("PtDCAXY7s"),
   fSystematicCuts(kFALSE),
   fShiftList(NULL),
   fUseShift(useshift),
   fAcceptTPC(kFALSE),
   fCheckPileUp(kFALSE),
   fIgnoreTPCzRange(kFALSE),
   fIgnoreTPCzRangeMax(FLT_MAX),
   fIgnoreTPCzRangeMin(-FLT_MAX),
   fCutMinimalTPCdedx(kFALSE),
   fMinimalTPCdedx(0.),
   fUsercuts(kFALSE),
   kaonMass(0.49368),
   pionMass(0.13957),
   kstarMass(0.89594),
   phiMass(1.01946),
   fCurrentEventCentralityBin(-1),
   fCurrentEventDoughterOne(0),
   fCurrentEventDoughterTwo(0),
   fCurrentEventVx(0),
   fCurrentEventVy(0),
   fCurrentEventVz(0),
   fBufferPointer(0),
   fQVector(0),
   fQVZeroA(0),
   fQVZeroC(0),
   fPairRapidity(1.0),
   fDeltaAngle(0),
   fHistZVertex(0), 
   fHistCentralityEvtCount(0),
   fHistEventCount(0),
   fHistCentDist(0),   
   fHistPionVsKaonMult(0),
   fHistTPCClusterRP(0),
   fHistITSClusterRP(0),
   fHistChiSqrPerNdfTPCRP(0),
   fHistChiSqrPerNdfITSRP(0),
   fHistdEdxKaonTpc(0),
   fHistdEdxKaonTof(0),
   fHistdEdxPionTpc(0),
   fHistdEdxPionTof(0),
   fHistTPCClusterPOI(0),
   fHistITSClusterPOI(0),
   fHistChiSqrPerNdfTPCPOI(0),
   fHistChiSqrPerNdfITSPOI(0),
   fHistnCrossRowsPOI(0),
   fHistratioCrossedRowsOverFindableClustersTPC(0),
   fHistDCAxyPOI(0),
   fHistDCAzPOI(0),
   fProCosResThreeSubEventPsiAB(0),
   fProCosResThreeSubEventPsiAC(0),
   fProCosResThreeSubEventPsiBC(0),
   fShiftCosTerm_V0A(0),
   fShiftSinTerm_V0A(0),
   fShiftCosTerm_V0C(0),
   fShiftSinTerm_V0C(0),
   fFullCosTerm_V0A(0),  
   fFullSinTerm_V0A(0),
   fFullCosTerm_V0C(0),  
   fFullSinTerm_V0C(0)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  kaonMass = 0.49368;
  pionMass = 0.13957;
  kstarMass = 0.89594;
  phiMass = 1.01946;

  printf("===================================================================================\n");
  printf("\n                HI I AM INSIDE .cxx after calling from ADDTASK                            \n");
  printf("===================================================================================\n");
  std::cout<<"**********************TASK NAME************************:"<<name<<std::endl;

  fhInvMassSAEPvzeroA= 0;
  fhInvMassLikePPSAEPvzeroA = 0;
  fhInvMassLikeMMSAEPvzeroA = 0;
  fhInvMassMixSAEPvzeroA= 0;
  fhInvMassSAEPvzeroC = 0;
  fhInvMassMixSAEPvzeroC = 0;
  for(Int_t i = 0; i < kDimBuf; i++){
    
    if(i < kCenBin)
      {
	/*
	fhInvMassLikePPSAEPvzeroA[i] = 0;
	fhInvMassLikeMMSAEPvzeroA[i] = 0;
	fhInvMassMixSAEPvzeroA[i] = 0;
	fhInvMassSAEPvzeroC[i] = 0;
	//fhInvMassLikePPSAEPvzeroC[i] = 0;
	//	fhInvMassLikeMMSAEPvzeroC[i] = 0;
	fhInvMassMixSAEPvzeroC[i] = 0;
	
*/
	fhInvMassv2[i]=0;
	fhInvMassSinv2[i]=0;
	fhNumInvMassvsPtPhi[i]=0;
	fhDenInvMassvsPtPhi[i]=0;
	fhDenInvMassvsPtPhi_likepp[i]=0;
	fhDenInvMassvsPtPhi_likemm[i]=0;
	
	fHistEPTPC[i] = 0;
	fHistEPV0A[i] = 0;
	fHistEPV0C[i] = 0;
	fHistEPFMDA[i] = 0;
	fHistEPFMDC[i] = 0;
	fHistEPZDCA[i] = 0;
	fHistEPZDCC[i] = 0;
	fHistEPT0A[i] = 0;
	fHistEPT0C[i] = 0;

	fHistEPV0A_rw[i] = 0;
	fHistEPV0C_rw[i] = 0;
	fHistEPV0A_plain[i] = 0;
	fHistEPV0C_plain[i] = 0;
	fHistEPV0A_rec[i] = 0;
	fHistEPV0C_rec[i] = 0;
	fHistEPV0A_recShift[i] = 0;
	fHistEPV0C_recShift[i] = 0;   
	fHistEPV0A_al[i] = 0;
	fHistEPV0C_al[i] = 0;

	fHistEPV0A_tw[i] = 0;
	fHistEPV0C_tw[i] = 0;
	fHistEPV0A_sc[i] = 0;
	fHistEPV0C_sc[i] = 0;
      } 

  if(i < 100){
    fBufferEventNEvents[i] = 0.0;
    fBufferEventFull[i] = 0.0;
  }
  
  fBufferEventDoughterOne[i] = 0;
  fBufferEventDoughterTwo[i] = 0;
  fBufferEventPsi[i] = 0;
  
  for(Int_t j = 0; j < kDim; j++){
    
    fBufferEventDoughterOneCharge[i][j] = 0;
    fBufferEventDoughterOneId[i][j] = 0;
    fBufferEventDoughterOneIn[i][j] = 0;
    fBufferEventDoughterOnePx[i][j] = 0;
    fBufferEventDoughterOnePy[i][j] = 0;
    fBufferEventDoughterOnePz[i][j] = 0;
    
    fBufferEventDoughterTwoCharge[i][j] = 0;
    fBufferEventDoughterTwoId[i][j] = 0;
      fBufferEventDoughterTwoIn[i][j] = 0;
      fBufferEventDoughterTwoPx[i][j] = 0;
      fBufferEventDoughterTwoPy[i][j] = 0;
      fBufferEventDoughterTwoPz[i][j] = 0;
      
      if( i == 0){
	fCurrentEventDoughterOneCharge[j] = 0;
	fCurrentEventDoughterOneId[j] = 0;
	fCurrentEventDoughterOneTrackNumber[j] = 0;
	fCurrentEventDoughterOneIn[j] = 0;
	fCurrentEventDoughterOnePx[j] = 0;
	fCurrentEventDoughterOnePy[j] = 0;
	fCurrentEventDoughterOnePz[j] = 0;
	fCurrentEventDoughterOneQx[j] = 0;
	fCurrentEventDoughterOneQy[j] = 0;
	
	fCurrentEventDoughterTwoCharge[j] = 0;
	fCurrentEventDoughterTwoIn[j] = 0;
	fCurrentEventDoughterTwoId[j] = 0;
	fCurrentEventDoughterTwoTrackNumber[j] = 0;
	fCurrentEventDoughterTwoPx[j] = 0;
	fCurrentEventDoughterTwoPy[j] = 0;
	fCurrentEventDoughterTwoPz[j] = 0;
	fCurrentEventDoughterTwoQx[j] = 0;
	fCurrentEventDoughterTwoQy[j] = 0;
      }
  }
  }
  
  DefineInput(0, TChain::Class());
  // Input slot #1 is needed for the weights input file
  if(fUseShift) {
    DefineInput(1, TList::Class());
  }
  DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskPhiSAR::~AliAnalysisTaskPhiSAR()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if(fAliESDtrackCuts) delete fAliESDtrackCuts;
  if(fESDtrackCutsForRP) delete fESDtrackCutsForRP;
  if(fPIDResponse)     delete fPIDResponse;
  if(fTriggerAna)      delete fTriggerAna;
  if(fUseShift)        delete fShiftList;
}

//________________________________________________________________________
void AliAnalysisTaskPhiSAR::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
  AliLog::SetClassDebugLevel("AliAnalysisTaskPhiSAR",10);
  
  
  //----------------------------------FlowQnVectorCorrections----------------------------------------------- 
  AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
    dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
  if (flowQnVectorTask != NULL) {
    fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
  else {
    AliFatal("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!");
  }
  

  //input hander
  if (fTriggerAna) delete fTriggerAna;
  fTriggerAna = new AliTriggerAnalysis;
  
  if(fAliESDtrackCuts) delete fAliESDtrackCuts;
  fAliESDtrackCuts = new AliESDtrackCuts();
  //if(fTrackType.CompareTo("GLOBAL")==0)fAliESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
  if(fTrackType.CompareTo("GLOBAL")==0)fAliESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  if(fTrackType.CompareTo("TPC")==0)fAliESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fAliESDtrackCuts->SetPtRange(0.1,20.0);
  fAliESDtrackCuts->SetEtaRange(-0.8,0.8);

  ///////////   For Systematic study   ///////////////////
  //Default parameters
  
  TString   PtDcaFormula              = "0.0182+0.0350/pt^1.01";//7sigma // 7*(0.0026+0.0050/pt^1.01)
  Double_t  dcazmax                   = 2.0;
  Double_t  maxX2TPCcls               = 4.0;
  Double_t  maxX2ITScls               = 36.0;
  Double_t  minCrossedRows            = 70.0;
  Double_t  minRatioClsCrRowsOverFCls = 0.8;

  

  if (fSystematicCuts){
    if(fSytCutType.Contains("PtDCAXY4s")) {PtDcaFormula = "0.0104+0.020/pt^1.01";}
    if(fSytCutType.Contains("PtDCAXY5s")) {PtDcaFormula = "0.013+0.025/pt^1.01";}
    if(fSytCutType.Contains("PtDCAXY6s")) {PtDcaFormula = "0.0156+0.030/pt^1.01";}
    if(fSytCutType.Contains("PtDCAXY7s")) {PtDcaFormula = "0.0182+0.035/pt^1.01";}//D
    if(fSytCutType.Contains("PtDCAXY8s")) {PtDcaFormula = "0.0208+0.040/pt^1.01";}
    if(fSytCutType.Contains("PtDCAXY9s")) {PtDcaFormula = "0.0234+0.045/pt^1.01";}
    
    if(fSytCutType.Contains("FixDCAZp1")) {dcazmax = 0.1;}
    if(fSytCutType.Contains("FixDCAZ1")) {dcazmax = 1.;}
    if(fSytCutType.Contains("FixDCAZ2")) {dcazmax = 2.;}//D
    if(fSytCutType.Contains("FixDCAZ3")) {dcazmax = 3.;}
    
    if(fSytCutType.Contains("NCrRows60")){minCrossedRows = 60;}
    if(fSytCutType.Contains("NCrRows70")){minCrossedRows = 70;}
    if(fSytCutType.Contains("NCrRows80")){minCrossedRows = 80;}//D
    if(fSytCutType.Contains("NCrRows90")){minCrossedRows = 90;}
    if(fSytCutType.Contains("NCrRows100")){minCrossedRows = 100;}
    
    if(fSytCutType.Contains("RClsCrRowsOvFCls0.7")){minRatioClsCrRowsOverFCls = 0.7;}
    if(fSytCutType.Contains("RClsCrRowsOvFCls0.8")){minRatioClsCrRowsOverFCls = 0.8;}//D
    if(fSytCutType.Contains("RClsCrRowsOvFCls0.9")){minRatioClsCrRowsOverFCls = 0.9;}
    
    if(fSytCutType.Contains("ChiSqrPerTPCCls3")) {maxX2TPCcls = 3.0;}
    if(fSytCutType.Contains("ChiSqrPerTPCCls4")) {maxX2TPCcls = 4.0;}//D
    if(fSytCutType.Contains("ChiSqrPerTPCCls5")) {maxX2TPCcls = 5.0;}
    
    if(fSytCutType.Contains("ChiSqrPerITSCls25")) {maxX2ITScls = 25.0;}
    if(fSytCutType.Contains("ChiSqrPerITSCls36")) {maxX2ITScls = 36.0;}//D
    if(fSytCutType.Contains("ChiSqrPerITSCls49")) {maxX2ITScls = 49.0;}

    
    fAliESDtrackCuts->SetMinNCrossedRowsTPC(minCrossedRows);
    fAliESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioClsCrRowsOverFCls);
    fAliESDtrackCuts->SetMaxChi2PerClusterTPC(maxX2TPCcls);
    fAliESDtrackCuts->SetMaxDCAToVertexXYPtDep(PtDcaFormula.Data());
    fAliESDtrackCuts->SetMaxDCAToVertexZ(dcazmax);
    fAliESDtrackCuts->SetMaxChi2PerClusterITS(maxX2ITScls);
  }
  
  if(fESDtrackCutsForRP) delete fESDtrackCutsForRP;
  fESDtrackCutsForRP = new AliESDtrackCuts();
  //if(fRPTrackType.CompareTo("GLOBAL")==0)fESDtrackCutsForRP = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
  if(fRPTrackType.CompareTo("GLOBAL")==0)fESDtrackCutsForRP = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  if(fRPTrackType.CompareTo("TPC")==0)fESDtrackCutsForRP = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCutsForRP->SetPtRange(0.15,20.0);
  fESDtrackCutsForRP->SetEtaRange(-0.8,0.8);
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");
  
  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");
 
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  
  Float_t pi = TMath::Pi();
  printf("========================okkkk <<");

  // Create histograms
  /*******************************************************
   *                    Event related Plots              *
   *******************************************************/
  //Event QA Histogram
  fHistZVertex            = new TH1F("fHistZVertex", "Z vertex distribution", 200,-20,20);       
  fHistCentralityEvtCount = new TH1F("fHistCentralityEvtCount", "Centrality Event Count", 11,-1,10);
  fHistEventCount         = new TH1F("fHistEventCount","fHistEventCount",10,0,10);
  fHistCentDist           = new TH1F("fHistCentDist","fHistCentDist",111,-1,110);
  //Q Vector
  //Eta Sub-Event Mult QA RP
  //POI Kaon Correlation
  fHistPionVsKaonMult           = new TH2F("fHistPionVsKaonMult","PionMult Vs KaonMult",2000,0,2000,500,0,500);

  /*******************************************************
   *                    Track related Plots              *
   *******************************************************/
  //Track QA Plot for RP
  fHistTPCClusterRP       = new TH1F("fHistTPCClusterRP","fHistTPCClusterRP",160,0,160);
  fHistITSClusterRP       = new TH1F("fHistITSClusterRP","fHistITSClusterRP",10,0,10);
  fHistChiSqrPerNdfTPCRP  = new TH1F("fHistChiSqrPerNdfTPCRP","fHistChiSqrPerNdfTPCRP",120,-1,5);
  fHistChiSqrPerNdfITSRP  = new TH1F("fHistChiSqrPerNdfITSRP","fHistChiSqrPerNdfITSRP",120,-1,5);
  
  //PID Selection for POI
  fHistdEdxKaonTpc        = new TH2F("fHistdEdxKaonTpc","Tpc  dEdx vs Mom kaon",1000,-10.,10.,500,0.,500);
  fHistdEdxKaonTof        = new TH2F("fHistdEdxKaonTof","Tof #beta vs Mom kaon",1000,-10.,10.,150,0.,1.5);
  fHistdEdxPionTpc        = new TH2F("fHistdEdxPionTpc","Tpc  dEdx vs Mom pion",1000,-10.,10.,500,0.,500);
  fHistdEdxPionTof        = new TH2F("fHistdEdxPionTof","Tof #beta vs Mom pion",1000,-10.,10.,150,0.,1.5);
  fHistTPCClusterPOI      = new TH1F("fHistTPCClusterPOI","fHistTPCClusterPOI",160,0,160);
  fHistITSClusterPOI      = new TH1F("fHistITSClusterPOI","fHistITSClusterPOI",100,0,50);
  fHistChiSqrPerNdfTPCPOI = new TH1F("fHistChiSqrPerNdfTPCPOI","fHistChiSqrPerNdfTPCPOI",120,-1,5);
  fHistChiSqrPerNdfITSPOI = new TH1F("fHistChiSqrPerNdfITSPOI","fHistChiSqrPerNdfITSPOI",120,-1,5);
  fHistnCrossRowsPOI      = new TH1F("fHistnCrossRowsPOI","fHistnCrossRowsPOI",160,0,160);
  fHistratioCrossedRowsOverFindableClustersTPC      = new TH1F("fHistratioCrossedRowsOverFindableClustersTPC","fHistratioCrossedRowsOverFindableClustersTPC",100,0,2);
  fHistDCAxyPOI      = new TH1F("fHistDCAxyPOI","fHistDCAxyPOI",200,-5,5);
  fHistDCAzPOI      = new TH1F("fHistDCAzPOI","fHistDCAzPOI",200,-5,5);
  
  //Event Plane Resolution Calculation (RP)

  fProCosResThreeSubEventPsiAB = new TProfile("fProCosResThreeSubEventPsiAB","fProCosResThreeSubEventPsiAB",11,-1,10,-1.,1.);
  fProCosResThreeSubEventPsiAC = new TProfile("fProCosResThreeSubEventPsiAC","fProCosResThreeSubEventPsiAC",11,-1,10,-1.,1.);
  fProCosResThreeSubEventPsiBC = new TProfile("fProCosResThreeSubEventPsiBC","fProCosResThreeSubEventPsiBC",11,-1,10,-1.,1.);

  //shift correction
  fFullCosTerm_V0A = new TProfile2D("fFullCosTerm_V0A","fFullCosTerm_V0A",10,0,10,20,1,21,-2.,2.);
  fFullSinTerm_V0A = new TProfile2D("fFullSinTerm_V0A","fFullSinTerm_V0A",10,0,10,20,1,21,-2.,2.);
  fFullCosTerm_V0C = new TProfile2D("fFullCosTerm_V0C","fFullCosTerm_V0C",10,0,10,20,1,21,-2.,2.);
  fFullSinTerm_V0C = new TProfile2D("fFullSinTerm_V0C","fFullSinTerm_V0C",10,0,10,20,1,21,-2.,2.);
  
  //Just QA
  

  /**************************************************************
   *     Invariant mass Histograms for different centrality     *
   *      Phi-Distribution of RP for Phi-weight calculation     *
   **************************************************************/

  //  for(Int_t i=0;i<kCenBin;i++)
  // {
      /*
      //This 3D is only for inveriant mass method
      fhInvMassSAEPvzeroA[i] = new TH3F(Form("fhInvMassSAEPvzeroA_Cen_%d",i),Form("fhInvMassSAEPvzeroA_SA_pt_Cen_%d",i),
					100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassLikePPSAEPvzeroA[i] = new TH3F(Form("fhInvMassLikePPSAEPvzeroA_Cen_%d",i),Form("fhInvMassLikePPSAEPvzeroA_SA_pt_Cen_%d",i),
					      100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassLikeMMSAEPvzeroA[i] = new TH3F(Form("fhInvMassLikeMMSAEPvzeroA_Cen_%d",i),Form("fhInvMassLikeMMSAEPvzeroA_SA_pt_Cen_%d",i),
					      100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassMixSAEPvzeroA[i] = new TH3F(Form("fhInvMassMixSAEPvzeroA_Cen_%d",i),Form("fhInvMassMixSAEPvzeroA_SA_pt_Cen_%d",i),
					   100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      
      fhInvMassSAEPvzeroC[i] = new TH3F(Form("fhInvMassSAEPvzeroC_Cen_%d",i),Form("fhInvMassSAEPvzeroC_SA_pt_Cen_%d",i),
					100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassLikePPSAEPvzeroC[i] = new TH3F(Form("fhInvMassLikePPSAEPvzeroC_Cen_%d",i),Form("fhInvMassLikePPSAEPvzeroC_SA_pt_Cen_%d",i),
					      100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassLikeMMSAEPvzeroC[i] = new TH3F(Form("fhInvMassLikeMMSAEPvzeroC_Cen_%d",i),Form("fhInvMassLikeMMSAEPvzeroC_SA_pt_Cen_%d",i),
					      100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);
      fhInvMassMixSAEPvzeroC[i] = new TH3F(Form("fhInvMassMixSAEPvzeroC_Cen_%d",i),Form("fhInvMassMixSAEPvzeroC_SA_pt_Cen_%d",i),
					   100,0.0,10.0,10,0.0,1.0, 250, 0.95, 1.2);

      */

  //This 3D is only for inveriant mass method
      
      Int_t bins1[4] = {100, 200, 90,18};
      Double_t xmin1[4] = {0.0, -1.0, 0.6,0.0};
      Double_t xmax1[4] = {10.0, 1.0, 1.5,90.0};
   
      fhInvMassSAEPvzeroA = new THnSparseD(Form("fhInvMassSAEPvzeroA"),Form("fhInvMassSAEPvzeroA_SA_pt"),4,bins1,xmin1,xmax1);  
      fhInvMassSAEPvzeroC = new THnSparseD(Form("fhInvMassSAEPvzeroC"),Form("fhInvMassSAEPvzeroC_SA_pt"),4,bins1,xmin1,xmax1);   
      //      fhInvMassSAEPvzeroA[i] = new THnSparseD(Form("fhInvMassSAEPvzeroA_Cen_%d",i),Form("fhInvMassSAEPvzeroA_SA_pt_Cen_%d",i),3,bins1,xmin1,xmax1);
      //  fhInvMassSAEPvzeroC[i] = new THnSparseD(Form("fhInvMassSAEPvzeroC_Cen_%d",i),Form("fhInvMassSAEPvzeroC_SA_pt_Cen_%d",i),3,bins1,xmin1,xmax1);
      
      Int_t bins[4] = {100,12, 90,18};
      Double_t xmin[4] = {0.0, 0.0, 0.6,0.0};
      Double_t xmax[4] = {10.0,pi, 1.5,90.0};



      fhInvMassMixSAEPvzeroA = new THnSparseD(Form("fhInvMassMixSAEPvzeroA"),Form("fhInvMassMixSAEPvzeroA_SA_pt"),4,bins,xmin,xmax);
      fhInvMassMixSAEPvzeroC = new THnSparseD(Form("fhInvMassMixSAEPvzeroC"),Form("fhInvMassMixSAEPvzeroC_SA_pt"),4,bins,xmin,xmax);     
      fhInvMassLikePPSAEPvzeroA = new THnSparseD(Form("fhInvMassLikePPSAEPvzeroA"),Form("fhInvMassLikePPSAEPvzeroA_SA_pt"),4,bins,xmin,xmax);
      fhInvMassLikeMMSAEPvzeroA = new THnSparseD(Form("fhInvMassLikeMMSAEPvzeroA"),Form("fhInvMassLikeMMSAEPvzeroA_SA_pt"),4,bins,xmin,xmax);

      /*
      fhInvMassMixSAEPvzeroA[i] = new THnSparseD(Form("fhInvMassMixSAEPvzeroA_Cen_%d",i),Form("fhInvMassMixSAEPvzeroA_SA_pt_Cen_%d",i),3,bins,xmin,xmax);
      fhInvMassMixSAEPvzeroC[i] = new THnSparseD(Form("fhInvMassMixSAEPvzeroC_Cen_%d",i),Form("fhInvMassMixSAEPvzeroC_SA_pt_Cen_%d",i),3,bins,xmin,xmax);     
      fhInvMassLikePPSAEPvzeroA[i] = new THnSparseD(Form("fhInvMassLikePPSAEPvzeroA_Cen_%d",i),Form("fhInvMassLikePPSAEPvzeroA_SA_pt_Cen_%d",i),3,bins,xmin,xmax);
      fhInvMassLikeMMSAEPvzeroA[i] = new THnSparseD(Form("fhInvMassLikeMMSAEPvzeroA_Cen_%d",i),Form("fhInvMassLikeMMSAEPvzeroA_SA_pt_Cen_%d",i),3,bins,xmin,xmax);
*/





      /*
      //This 3D is only for inveriant mass method
      fhInvMassv2[i]    = new TH3F(Form("InvMassv2_Cen_%d",i),Form("InvMass_v2_pt_Cen_%d",i),
				   100,0.0,10.0,200,-1.0,1.0, 90, 0.6, 1.5);
      fhInvMassSinv2[i] = new TH3F(Form("InvMassv2_Cen_Sin%d",i),Form("InvMass_v2_pt_Cen_Sin%d",i),
				   100,0.0,10.0,200,-1.0,1.0, 90, 0.6, 1.5);
      
      //For phi-psi bin method Also needed for invariant method
      fhNumInvMassvsPtPhi[i] = new TH3F(Form("NumInvMassvsPtPhi_Cen_%d",i), Form("cen %d same event spectrum",i), 
					100,0.,10.,12,0.,pi, 90, 0.6, 1.5);
      fhDenInvMassvsPtPhi[i] = new TH3F(Form("DenInvMassvsPtPhi_Cen_%d",i), Form("cen %d mixed event spectrum",i),
					100,0.,10.,12,0.,pi, 90, 0.6, 1.5);

      fhDenInvMassvsPtPhi_likepp[i] = new TH3F(Form("NumInvMassvsPtPhi_likepp_Cen_%d",i), Form("cen %d like pp event spectrum",i), 
					100,0.,10.,12,0.,pi, 90, 0.6, 1.5);
      fhDenInvMassvsPtPhi_likemm[i] = new TH3F(Form("DenInvMassvsPtPhi_likemm_Cen_%d",i), Form("cen %d like mm event spectrum",i),
					100,0.,10.,12,0.,pi, 90, 0.6, 1.5);
      */
      //Phi-distribution for RP

      //Multiplicity Distribution POI
    
      for(Int_t i=0;i<kCenBin;i++)
	{
      fHistEPTPC[i] = new TH1F(Form("fHistEPTPC%d",i),Form("fHistEPTPC%d",i),180,0.,pi);
      fHistEPV0A[i] =  new TH1F(Form("fHistEPV0AA%d",i),Form("fHistEPV0AA%d",i),180,0.,pi);
      fHistEPV0C[i] =  new TH1F(Form("fHistEPV0CC%d",i),Form("fHistEPV0CC%d",i),180,0.,pi);
      fHistEPFMDA[i] =  new TH1F(Form("fHistEPFMDA%d",i),Form("fHistEPFMDA%d",i),180,0.,pi);
      fHistEPFMDC[i] =  new TH1F(Form("fHistEPFMDC%d",i),Form("fHistEPFMDC%d",i),180,0.,pi);
      fHistEPZDCA[i] =  new TH1F(Form("fHistEPZDCA%d",i),Form("fHistEPZDCA%d",i),180,0.,pi);
      fHistEPZDCC[i] =  new TH1F(Form("fHistEPZDCC%d",i),Form("fHistEPZDCC%d",i),180,0.,pi);
      fHistEPT0A[i] =  new TH1F(Form("fHistEPT0A%d",i),Form("fHistEPT0A%d",i),180,0.,pi);
      fHistEPT0C[i] =  new TH1F(Form("fHistEPT0C%d",i),Form("fHistEPT0C%d",i),180,0.,pi);

      fHistEPV0A_rw[i] = new TH1F(Form("fHistEPV0A_rw%d",i),Form("fHistEPV0A_rw%d",i),180,0.,pi);
      fHistEPV0C_rw[i] = new TH1F(Form("fHistEPV0C_rw%d",i),Form("fHistEPV0C_rw%d",i),180,0.,pi);
      fHistEPV0A_plain[i] = new TH1F(Form("fHistEPV0A_plain%d",i),Form("fHistEPV0A_plain%d",i),180,0.,pi);
      fHistEPV0C_plain[i] = new TH1F(Form("fHistEPV0C_plain%d",i),Form("fHistEPV0C_plain%d",i),180,0.,pi);
      fHistEPV0A_rec[i] = new TH1F(Form("fHistEPV0A_rec%d",i),Form("fHistEPV0A_rec%d",i),180,0.,pi);
      fHistEPV0C_rec[i] = new TH1F(Form("fHistEPV0C_rec%d",i),Form("fHistEPV0C_rec%d",i),180,0.,pi);
      fHistEPV0A_recShift[i] = new TH1F(Form("fHistEPV0A_recShift%d",i),Form("fHistEPV0A_recShift%d",i),180,0.,pi);
      fHistEPV0C_recShift[i] = new TH1F(Form("fHistEPV0C_recShift%d",i),Form("fHistEPV0C_recShift%d",i),180,0.,pi);
      fHistEPV0A_al[i] = new TH1F(Form("fHistEPV0A_al%d",i),Form("fHistEPV0A_al%d",i),180,0.,pi);
      fHistEPV0C_al[i] = new TH1F(Form("fHistEPV0C_al%d",i),Form("fHistEPV0C_al%d",i),180,0.,pi);

      fHistEPV0A_tw[i] = new TH1F(Form("fHistEPV0A_tw%d",i),Form("fHistEPV0A_tw%d",i),180,0.,pi);
      fHistEPV0C_tw[i] = new TH1F(Form("fHistEPV0C_tw%d",i),Form("fHistEPV0C_tw%d",i),180,0.,pi);
      fHistEPV0A_sc[i] = new TH1F(Form("fHistEPV0A_sc%d",i),Form("fHistEPV0A_sc%d",i),180,0.,pi);
      fHistEPV0C_sc[i] = new TH1F(Form("fHistEPV0C_sc%d",i),Form("fHistEPV0C_sc%d",i),180,0.,pi);

      /*
      fOutput->Add(fhInvMassv2[i]);
      fOutput->Add(fhInvMassSinv2[i]);
      fOutput->Add(fhNumInvMassvsPtPhi[i]);
      fOutput->Add(fhDenInvMassvsPtPhi[i]);
      */
      //      fOutput->Add( fhDenInvMassvsPtPhi_likepp[i]);
      //  fOutput->Add( fhDenInvMassvsPtPhi_likemm[i]);

      
      fOutput->Add(fhInvMassSAEPvzeroA);
      fOutput->Add(fhInvMassSAEPvzeroC);     
      fOutput->Add(fhInvMassLikePPSAEPvzeroA);
      fOutput->Add(fhInvMassLikeMMSAEPvzeroA);
      fOutput->Add(fhInvMassMixSAEPvzeroA);
      fOutput->Add(fhInvMassMixSAEPvzeroC);      

      /*     fOutput->Add(fhInvMassSAEPvzeroA[i]);
      fOutput->Add(fhInvMassSAEPvzeroC[i]);     
      fOutput->Add(fhInvMassLikePPSAEPvzeroA[i]);
      fOutput->Add(fhInvMassLikeMMSAEPvzeroA[i]);
      fOutput->Add(fhInvMassMixSAEPvzeroA[i]);
      fOutput->Add(fhInvMassMixSAEPvzeroC[i]);
      */
      
      
      /*    //
      fOutput->Add(fHistEPTPC[i]);
      fOutput->Add(fHistEPV0A[i]);
      fOutput->Add(fHistEPV0C[i]);
      fOutput->Add(fHistEPFMDA[i]);
      fOutput->Add(fHistEPFMDC[i]);
      fOutput->Add(fHistEPZDCA[i]);
      fOutput->Add(fHistEPZDCC[i]);
      fOutput->Add(fHistEPT0A[i]);
      fOutput->Add(fHistEPT0C[i]);
      fOutput->Add(fHistEPV0A_rw[i]);
      fOutput->Add(fHistEPV0C_rw[i]);
      fOutput->Add(fHistEPV0A_plain[i]);
      fOutput->Add(fHistEPV0C_plain[i]);
      fOutput->Add(fHistEPV0A_rec[i]);
      fOutput->Add(fHistEPV0C_rec[i]);
      fOutput->Add(fHistEPV0A_recShift[i]);
      fOutput->Add(fHistEPV0C_recShift[i]);
      fOutput->Add(fHistEPV0A_al[i]);
      fOutput->Add(fHistEPV0C_al[i]);
      fOutput->Add(fHistEPV0A_tw[i]);
      fOutput->Add(fHistEPV0C_tw[i]);
      fOutput->Add(fHistEPV0A_sc[i]);
      fOutput->Add(fHistEPV0C_sc[i]);
      */
    }
  // NEW HISTO should be defined here, with a sensible name,
  fOutput->Add(fHistZVertex);
  fOutput->Add(fHistCentralityEvtCount);
  /*
  fOutput->Add(fHistEventCount);
  fOutput->Add(fHistCentDist);
 
  fOutput->Add(fHistPionVsKaonMult);
  
  fOutput->Add(fHistTPCClusterRP);
  fOutput->Add(fHistITSClusterRP);
  fOutput->Add(fHistChiSqrPerNdfTPCRP);
  fOutput->Add(fHistChiSqrPerNdfITSRP);

  fOutput->Add(fHistdEdxKaonTpc);
  fOutput->Add(fHistdEdxKaonTof);
  fOutput->Add(fHistdEdxPionTpc);
  fOutput->Add(fHistdEdxPionTof);
  fOutput->Add(fHistTPCClusterPOI);
  fOutput->Add(fHistITSClusterPOI);
  fOutput->Add(fHistChiSqrPerNdfTPCPOI);
  fOutput->Add(fHistChiSqrPerNdfITSPOI);
  fOutput->Add(fHistnCrossRowsPOI);
  fOutput->Add(fHistratioCrossedRowsOverFindableClustersTPC);
  fOutput->Add(fHistDCAxyPOI);
  fOutput->Add(fHistDCAzPOI);
  */
  fOutput->Add(fProCosResThreeSubEventPsiAB);
  fOutput->Add(fProCosResThreeSubEventPsiAC);
  fOutput->Add(fProCosResThreeSubEventPsiBC);

  fOutput->Add(fFullCosTerm_V0A);
  fOutput->Add(fFullSinTerm_V0A);
  fOutput->Add(fFullCosTerm_V0C);
  fOutput->Add(fFullSinTerm_V0C);

  
  PostData(1, fOutput);

  
  // For shift correction ***************************
  
  if(fUseShift) {
    if(!fShiftList) {
      std::cout<<"WARNING: fWeightsList or fShiftList is NULL."<<std::endl;
      exit(0);  
    }
    fShiftList->SetOwner();
    fShiftList->SetName("ShiftCorr");
    //fShiftList->ls();
    
    //For v2 ++++++++++++++++
    if(fShiftList->FindObject("full_cos_Term_V0A"))  {
      fShiftCosTerm_V0A = dynamic_cast<TProfile2D*>
	(fShiftList->FindObject("full_cos_Term_V0A"));
      //fOutput->Add(fShiftCosTerm_v2);
    } else {
      std::cout<<"WARNING: profile with cosine av.in shift corr. for V0A is not accessible"<<std::endl;
      exit(0);
    }
    
    if(fShiftList->FindObject("full_sin_Term_V0A"))  {
      fShiftSinTerm_V0A = dynamic_cast<TProfile2D*>
	(fShiftList->FindObject("full_sin_Term_V0A"));
      //fOutput->Add(fShiftSinTerm_v2);
    } else {
      std::cout<<"WARNING: profile with sine av. in shift corr. for V0A is not accessible"<<std::endl;
      exit(0);
    }
    
    if(fShiftList->FindObject("full_cos_Term_V0C"))  {
      fShiftCosTerm_V0C = dynamic_cast<TProfile2D*>
	(fShiftList->FindObject("full_cos_Term_V0C"));
      //fOutput->Add(fShiftCosTerm_v2);
    } else {
      std::cout<<"WARNING: profile with cosine av.in shift corr. for V0C is not accessible"<<std::endl;
      exit(0);
    }
    
    if(fShiftList->FindObject("full_sin_Term_V0C"))  {
      fShiftSinTerm_V0C = dynamic_cast<TProfile2D*>
	(fShiftList->FindObject("full_sin_Term_V0C"));
      //fOutput->Add(fShiftSinTerm_v2);
    } else {
      std::cout<<"WARNING: profile with sine av. in shift corr. for V0C is not accessible"<<std::endl;
      exit(0);
    }
  } // end of if(fUseShift)
  
  // Post data for ALL output slots >0 here, to get at least an empty histogram
  //PostData(21,fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskPhiSAR::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  // Create pointer to reconstructed event

  
  Int_t myHarmonic = 2;
  const AliQnCorrectionsQnVector *myQnVector;
  Double_t myEventPlaneTPC = 0.0;
  Double_t myEventPlaneVZEROA = 0.0;
  Double_t myEventPlaneVZEROC = 0.0;
  Double_t myEventPlaneZDCA = 0.0;
  Double_t myEventPlaneZDCC = 0.0;
  Double_t myEventPlaneFMDA = 0.0;
  Double_t myEventPlaneFMDC = 0.0;
  Double_t myEventPlaneT0A = 0.0;
  Double_t myEventPlaneT0C = 0.0;
  
  double myEventPlaneV0A_rw=0, myEventPlaneV0C_rw=0; 
  double myEventPlaneV0A_plain=0, myEventPlaneV0C_plain=0; 
  double myEventPlaneV0A_rec=0, myEventPlaneV0C_rec=0; 
  double myEventPlaneV0A_al=0, myEventPlaneV0C_al=0; 
  double myEventPlaneV0A_tw=0, myEventPlaneV0C_tw=0; 
  double myEventPlaneV0A_sc=0, myEventPlaneV0C_sc=0; 

  Double_t psi_Shift_V0A =0.,psi_Shift_V0C =0.;

  AliEventplane *esdEP;
  TVector2 qq1,qq2,rm1,rm2,vzeroa,vzeroc,vzeroa_new,vzeroc_new;
  Double_t psiV0A=0., psiV0C=0.;
  TObjArray tracklistusedinEP;

  if (fAnalysisLevel.CompareTo("ESD")==0){
    AliVEvent *event = InputEvent();
    if (!event) {Printf("ERROR: Could not retrieve event\n"); PostData(1, fOutput);return; }
    if (!fPIDResponse){Printf("ERROR: Could not retrieve PID Reasponse\n"); PostData(1, fOutput);return; }
    ReSet();
    
    /*
    AliTimeRangeCut  *fTimeRangeCut;
    fTimeRangeCut = new AliTimeRangeCut;
    fTimeRangeCut->InitFromEvent(event);
    if (fTimeRangeCut->CutEvent(event))
      return;
    */

    if(PassEvent(event)){
      
      Int_t centBin = GetCentrality(event);  
      
      Int_t nRP = 0, passtrack = 0;
      TObjArray* tracklist = new TObjArray;
      AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
      //AliCentrality *centrality = esd->GetCentrality();
      if (esd){    
	//if (!(fRunNumber == esd->GetRunNumber())) {
	//fRunNumber = esd->GetRunNumber();
	//}
	esdEP = esd->GetEventplane();
	if(fRPTrackType.CompareTo("GLOBAL")==0) tracklist = fESDtrackCutsForRP->GetAcceptedTracks(esd,kFALSE);
	if(fRPTrackType.CompareTo("TPC")==0)    tracklist = fESDtrackCutsForRP->GetAcceptedTracks(esd,kTRUE);
	//Total number of tracks after the track cuts
	const Int_t nt = tracklist->GetEntries();
	
	// get the fully corrected Qn vector from VZEROA sub-detector 
	fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);

	myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPC");
      if (myQnVector != NULL) myEventPlaneTPC = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneTPC < 0.)myEventPlaneTPC = myEventPlaneTPC + TMath::Pi();
      fHistEPTPC[centBin]->Fill(myEventPlaneTPC);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA");
      if (myQnVector != NULL) myEventPlaneVZEROA = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneVZEROA < 0.)myEventPlaneVZEROA = myEventPlaneVZEROA + TMath::Pi();
      fHistEPV0A[centBin]->Fill(myEventPlaneVZEROA);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC");
      if (myQnVector != NULL) myEventPlaneVZEROC = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneVZEROC < 0.)myEventPlaneVZEROC = myEventPlaneVZEROC + TMath::Pi();
      fHistEPV0C[centBin]->Fill(myEventPlaneVZEROC);      

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("ZDCA");
      if (myQnVector != NULL) myEventPlaneZDCA = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneZDCA < 0.)myEventPlaneZDCA = myEventPlaneZDCA + TMath::Pi();
      fHistEPZDCA[centBin]->Fill(myEventPlaneZDCA);


      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("ZDCC");
      if (myQnVector != NULL) myEventPlaneZDCC = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneZDCC < 0.)myEventPlaneZDCC = myEventPlaneZDCC + TMath::Pi();
      fHistEPZDCC[centBin]->Fill(myEventPlaneZDCC);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("FMDA");
      if (myQnVector != NULL) myEventPlaneFMDA = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneFMDA < 0.)myEventPlaneFMDA = myEventPlaneFMDA + TMath::Pi();
      fHistEPFMDA[centBin]->Fill(myEventPlaneFMDA);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("FMDC");
      if (myQnVector != NULL) myEventPlaneFMDC = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneFMDC < 0.)myEventPlaneFMDC = myEventPlaneFMDC + TMath::Pi();
      fHistEPFMDC[centBin]->Fill(myEventPlaneFMDC);      

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TZEROA");
      if (myQnVector != NULL) myEventPlaneT0A = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneT0A < 0.)myEventPlaneT0A = myEventPlaneT0A + TMath::Pi();
      fHistEPT0A[centBin]->Fill(myEventPlaneT0A);      

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TZEROC");
      if (myQnVector != NULL) myEventPlaneT0C = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneT0C < 0.)myEventPlaneT0C = myEventPlaneT0C + TMath::Pi();
      fHistEPT0C[centBin]->Fill(myEventPlaneT0C);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","raw","raw");
      if (myQnVector != NULL) myEventPlaneV0A_rw = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_rw < 0.)myEventPlaneV0A_rw = myEventPlaneV0A_rw + TMath::Pi();
      fHistEPV0A_rw[centBin]->Fill(myEventPlaneV0A_rw);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","raw","raw");
      if (myQnVector != NULL) myEventPlaneV0C_rw = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_rw < 0.)myEventPlaneV0C_rw = myEventPlaneV0C_rw + TMath::Pi();
      fHistEPV0C_rw[centBin]->Fill(myEventPlaneV0C_rw);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","plain","plain");
      if (myQnVector != NULL) myEventPlaneV0A_plain = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_plain < 0.)myEventPlaneV0A_plain = myEventPlaneV0A_plain + TMath::Pi();
      fHistEPV0A_plain[centBin]->Fill(myEventPlaneV0A_plain);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","plain","plain");
      if (myQnVector != NULL) myEventPlaneV0C_plain = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_plain < 0.)myEventPlaneV0C_plain = myEventPlaneV0C_plain + TMath::Pi();
      fHistEPV0C_plain[centBin]->Fill(myEventPlaneV0C_plain);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rec","rec");
      if (myQnVector != NULL) myEventPlaneV0A_rec = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_rec < 0.)myEventPlaneV0A_rec = myEventPlaneV0A_rec + TMath::Pi();
      fHistEPV0A_rec[centBin]->Fill(myEventPlaneV0A_rec);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rec","rec");
      if (myQnVector != NULL) myEventPlaneV0C_rec = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_rec < 0.)myEventPlaneV0C_rec = myEventPlaneV0C_rec + TMath::Pi();
      fHistEPV0C_rec[centBin]->Fill(myEventPlaneV0C_rec);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","align","align");
      if (myQnVector != NULL) myEventPlaneV0A_al = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_al < 0.)myEventPlaneV0A_al = myEventPlaneV0A_al + TMath::Pi();
      fHistEPV0A_al[centBin]->Fill(myEventPlaneV0A_al);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","align","align");
      if (myQnVector != NULL) myEventPlaneV0C_al = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_al < 0.)myEventPlaneV0C_al = myEventPlaneV0C_al + TMath::Pi();
      fHistEPV0C_al[centBin]->Fill(myEventPlaneV0C_al);


      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","twist","twist");
      if (myQnVector != NULL) myEventPlaneV0A_tw = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_tw < 0.)myEventPlaneV0A_tw = myEventPlaneV0A_tw + TMath::Pi();
      fHistEPV0A_tw[centBin]->Fill(myEventPlaneV0A_tw);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","twist","twist");
      if (myQnVector != NULL) myEventPlaneV0C_tw = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_tw < 0.)myEventPlaneV0C_tw = myEventPlaneV0C_tw + TMath::Pi();
      fHistEPV0C_tw[centBin]->Fill(myEventPlaneV0C_tw);


      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","scale","scale");
      if (myQnVector != NULL) myEventPlaneV0A_sc = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0A_sc < 0.)myEventPlaneV0A_sc = myEventPlaneV0A_sc + TMath::Pi();
      fHistEPV0A_sc[centBin]->Fill(myEventPlaneV0A_sc);

      myQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","scale","scale");
      if (myQnVector != NULL) myEventPlaneV0C_sc = myQnVector->EventPlane(myHarmonic);
      if(myEventPlaneV0C_sc < 0.)myEventPlaneV0C_sc = myEventPlaneV0C_sc + TMath::Pi();
      fHistEPV0C_sc[centBin]->Fill(myEventPlaneV0C_sc);


      //Calculate shift correction factor
      if(!fUseShift){
	for(Int_t kk=1;kk<21;kk++){
	  //Shift correction for Psi_2 ==========
	  fFullCosTerm_V0A->Fill(centBin,kk, TMath::Cos(2*kk*myEventPlaneV0A_rec));
	  fFullSinTerm_V0A->Fill(centBin,kk, TMath::Sin(2*kk*myEventPlaneV0A_rec));
	  
	  fFullCosTerm_V0C->Fill(centBin,kk, TMath::Cos(2*kk*myEventPlaneV0C_rec));
	  fFullSinTerm_V0C->Fill(centBin,kk, TMath::Sin(2*kk*myEventPlaneV0C_rec));
	}
      }///===========================
      
      //When shift correction is kTRUE
      if(fUseShift){
	for(Int_t kk=1;kk<21;kk++){
	  //Shift Correction for v2 +++++++++++++++
	  psi_Shift_V0A+= (1.0/(1.0*kk))*(-(fShiftSinTerm_V0A->GetBinContent(centBin+1,kk))*cos(2.*kk*myEventPlaneV0A_rec) + (fShiftCosTerm_V0A->GetBinContent(centBin+1,kk))*sin(2.*kk*myEventPlaneV0A_rec));
	  psi_Shift_V0C+= (1.0/(1.0*kk))*(-(fShiftSinTerm_V0C->GetBinContent(centBin+1,kk))*cos(2.*kk*myEventPlaneV0C_rec) + (fShiftCosTerm_V0C->GetBinContent(centBin+1,kk))*sin(2.*kk*myEventPlaneV0C_rec));
	}
	
	if(psi_Shift_V0A < 0) psi_Shift_V0A +=  TMath::Pi();
	if(psi_Shift_V0A > TMath::Pi()) psi_Shift_V0A -= TMath::Pi();
	
	if(psi_Shift_V0C < 0) psi_Shift_V0C += TMath::Pi();
	if(psi_Shift_V0C > TMath::Pi()) psi_Shift_V0C -= TMath::Pi();
	
	myEventPlaneV0A_rec = myEventPlaneV0A_rec + psi_Shift_V0A;
	if(myEventPlaneV0A_rec < 0) myEventPlaneV0A_rec += TMath::Pi();
	if(myEventPlaneV0A_rec > TMath::Pi()) myEventPlaneV0A_rec -= TMath::Pi();
	
	myEventPlaneV0C_rec = myEventPlaneV0C_rec + psi_Shift_V0C;
	if(myEventPlaneV0C_rec < 0) myEventPlaneV0C_rec += TMath::Pi();
	if(myEventPlaneV0C_rec > TMath::Pi()) myEventPlaneV0C_rec -= TMath::Pi();
      }
      
      fHistEPV0A_recShift[centBin]->Fill(myEventPlaneV0A_rec);
      fHistEPV0C_recShift[centBin]->Fill(myEventPlaneV0C_rec);

      //raw, plain, rec, align, twist, scale
      // if(myEventPlaneTPC!=0)std::cout<<" My event plane angle ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "<<myEventPlaneTPC<<"\t"<<myEventPlaneVZEROA<<"\t"<<myEventPlaneVZEROC<<"\t"<<myEventPlaneZDCA<<"\t"<<myEventPlaneZDCC<<"\t"<<myEventPlaneFMDA<<"\t"<<myEventPlaneFMDC<<std::endl;

       TVector2 QfullTPC;
       TVector2 QfullV0A;
       TVector2 QfullV0C;
       //myEventPlaneV0A_tw       
       /*       QfullTPC.Set(TMath::Cos(2.0*myEventPlaneTPC), TMath::Sin(2.0*myEventPlaneTPC));
       QfullV0A.Set(TMath::Cos(2.0*myEventPlaneV0A_rec), TMath::Sin(2.0*myEventPlaneV0A_rec));
       QfullV0C.Set(TMath::Cos(2.0*myEventPlaneV0C_rec), TMath::Sin(2.0*myEventPlaneV0C_rec));
       */
       QfullTPC.Set(TMath::Cos(2.0*myEventPlaneTPC), TMath::Sin(2.0*myEventPlaneTPC));
       QfullV0A.Set(TMath::Cos(2.0*myEventPlaneV0A_tw), TMath::Sin(2.0*myEventPlaneV0A_tw));
       QfullV0C.Set(TMath::Cos(2.0*myEventPlaneV0C_tw), TMath::Sin(2.0*myEventPlaneV0C_tw));



       fQVZeroA = new TVector2(QfullV0A);
       fQVZeroC = new TVector2(QfullV0C);
       fQVector = new TVector2(QfullTPC);

       //if(fQVZeroA->Mod()!=0)fHistEPV0A[centBin]->Fill(fQVZeroA->Phi()/2.);//Centrality wise V0A EP
       //if(fQVZeroC->Mod()!=0)fHistEPV0C[centBin]->Fill(fQVZeroC->Phi()/2.);//Centrality wise V0A EP
       
       //std::cout<<"In user e  ==== === "<<fQVZeroA->Phi()/2.<<"\t"<<fQVZeroA_new->Phi()/2.<<"\t"<<fQVZeroC->Phi()/2.<<"\t"<<fQVZeroC_new->Phi()/2.<<std::endl;
       
	if(centBin >=0 && centBin < kCenBin){
	  
	  //for resolution correction
	  fProCosResThreeSubEventPsiAB->Fill(centBin, TMath::Cos(2.0 *(myEventPlaneV0A_tw - myEventPlaneTPC)));
	  //	  fProCosResThreeSubEventPsiAC->Fill(centBin, TMath::Cos(2.0 *(myEventPlaneV0A_rec - myEventPlaneV0C_rec)))
	    fProCosResThreeSubEventPsiAC->Fill(centBin,myEventPlaneV0A_tw);
	  //	  fProCosResThreeSubEventPsiBC->Fill(centBin, TMath::Cos(2.0 *(myEventPlaneTPC - myEventPlaneV0C_rec)));
	  fProCosResThreeSubEventPsiBC->Fill(centBin, TMath::Cos(2.0 *centBin));
	}
	
	//for resolution correction
	fProCosResThreeSubEventPsiAB->Fill(-1, TMath::Cos(2.0 *(myEventPlaneV0A_tw - myEventPlaneTPC)));
	fProCosResThreeSubEventPsiAC->Fill(-1, TMath::Cos(2.0 *(myEventPlaneV0A_tw - myEventPlaneV0C_tw)));
	fProCosResThreeSubEventPsiBC->Fill(-1, TMath::Cos(2.0 *(myEventPlaneTPC - myEventPlaneV0C_tw)));
	  
	
	fHistZVertex->Fill(fCurrentEventVz);
	fHistCentralityEvtCount->Fill(centBin);
	


	
	
	Int_t ntracks = event->GetNumberOfTracks();
	//	std::cout<<" no of track  "<<ntracks<<endl;
	for(Int_t i = 0; i < ntracks; i++) {
	  
	  AliVTrack   *track = (AliVTrack*)event->GetTrack(i);
	  if(!track) { 
	    AliError(Form("ERROR: Could not retrieve esdtrack %d",i)); 
	    continue; 
	  }
	  
	  if(PassTrack(track))
	    {

	      fHistDCAzPOI->Fill(1);
	      passtrack++;
	      
	      Float_t pTrackTpc = track->GetTPCmomentum();
	      Float_t pTrack    = track->P();
	      Float_t pTTrack   = track->Pt();
	      Float_t phi       = track->Phi();
	      Float_t eta       = track->Eta();
	      Short_t charge    = track->Charge();
	      Int_t   trackId   = track->GetID();
	      Float_t TpcSignal = track->GetTPCsignal();
	      
	      //Double_t nSigma   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	      //needs for auto correlation removeal. But presently the kaons are not added to the event plane calculation.
	      Float_t Qxi       = TMath::Cos(2.*phi);//cos(2phi)
	      Float_t Qyi       = TMath::Sin(2.*phi);//sin(2phi)
	      //Double_t nSigCut  = GetNSigmaCut(pTrackTpc);
	      AliESDtrack *esdt  = dynamic_cast<AliESDtrack*>(track);
	      Float_t dcaxy = 0.0;
	      Float_t dcaz  = 0.0;
	      esdt->GetImpactParameters(dcaxy,dcaz);
	      Int_t nClustersITS = esdt->GetITSclusters(0);
	      Int_t nClustersTPC = -1;
	      if(fAliESDtrackCuts->GetRequireTPCStandAlone()) {
		nClustersTPC = esdt->GetTPCNclsIter1();
	      }
	      else {
		nClustersTPC = esdt->GetTPCclusters(0);
	      }
	      Float_t chi2PerClusterITS = -1;
	      Float_t chi2PerClusterTPC = -1;
	      if (nClustersITS!=0)
		chi2PerClusterITS = esdt->GetITSchi2()/Float_t(nClustersITS);
	      if (nClustersTPC!=0) {
		if(fAliESDtrackCuts->GetRequireTPCStandAlone()) {
		  chi2PerClusterTPC = esdt->GetTPCchi2Iter1()/Float_t(nClustersTPC);
		} else {
		  chi2PerClusterTPC = esdt->GetTPCchi2()/Float_t(nClustersTPC);
		}
	      }
	      Float_t nCrossedRowsTPC = esdt->GetTPCCrossedRows();
	      Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
	      if (esdt->GetTPCNclsF()>0) {
		ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / esdt->GetTPCNclsF();
	      }
	      	      
	      //one for pion other for kaon sub
	      Bool_t   isTOF     = MatchTOF(track);
	      Double_t tofbeta;
	      if(isTOF){
		tofbeta = GetTOFBeta(track);
		//		std::cout<<" tof  ============== loop"<< tofbeta<<endl;
	      }
	      //	      fHistTPCClusterPOI->Fill(1);
	      fHistTPCClusterPOI->Fill(nClustersTPC);
	      fHistITSClusterPOI->Fill(nClustersITS);
	      fHistChiSqrPerNdfTPCPOI->Fill(chi2PerClusterTPC);
	      fHistChiSqrPerNdfITSPOI->Fill(chi2PerClusterITS);
	      fHistnCrossRowsPOI->Fill(nCrossedRowsTPC);
	      fHistratioCrossedRowsOverFindableClustersTPC->Fill(ratioCrossedRowsOverFindableClustersTPC);
	      fHistDCAxyPOI->Fill(dcaxy);
	      //	      fHistDCAzPOI->Fill(dcaz);
	      
	      //select Pions POI
	      if(IsSelectedPION(track))// for pi-/pi+
		{
		  //		  std::cout<<" pion selected  ======--"<<endl; 
		  fCurrentEventDoughterTwoCharge[fCurrentEventDoughterTwo] = charge;
		  fCurrentEventDoughterTwoPx[fCurrentEventDoughterTwo] = track->Px();
		  fCurrentEventDoughterTwoPy[fCurrentEventDoughterTwo] = track->Py();
		  fCurrentEventDoughterTwoPz[fCurrentEventDoughterTwo] = track->Pz();
		  fCurrentEventDoughterTwoId[fCurrentEventDoughterTwo] = trackId;
		  fCurrentEventDoughterTwoTrackNumber[fCurrentEventDoughterTwo] = passtrack;
		  fCurrentEventDoughterTwoQx[fCurrentEventDoughterTwo] = Qxi;
		  fCurrentEventDoughterTwoQy[fCurrentEventDoughterTwo] = Qyi;
		 
		  if(tracklistusedinEP.Contains(track))
		    fCurrentEventDoughterTwoIn[fCurrentEventDoughterTwo] = 1;//dummy
		  else fCurrentEventDoughterTwoIn[fCurrentEventDoughterTwo] = 0;//dummy
		  fCurrentEventDoughterTwo++;
		  fHistdEdxPionTpc->Fill(pTrackTpc*charge, TpcSignal);
		  if(isTOF)fHistdEdxPionTof->Fill(pTrack*charge, tofbeta);
		}
	      //select Kaons POI
	      //if(!IsSelectedPION(track) && IsSelectedKAON(track))//K+/K-
	      if(IsSelectedKAON(track))//K+/K-
		{
		  //		  std::cout<<" kaon selected  ======--"<<endl; 
		  fCurrentEventDoughterOneCharge[fCurrentEventDoughterOne] = charge;
		  fCurrentEventDoughterOnePx[fCurrentEventDoughterOne] = track->Px();
		  fCurrentEventDoughterOnePy[fCurrentEventDoughterOne] = track->Py();
		  fCurrentEventDoughterOnePz[fCurrentEventDoughterOne] = track->Pz();
		  fCurrentEventDoughterOneId[fCurrentEventDoughterOne] = trackId;
		  fCurrentEventDoughterOneTrackNumber[fCurrentEventDoughterOne] = passtrack;
		  fCurrentEventDoughterOneQx[fCurrentEventDoughterOne] = Qxi;
		  fCurrentEventDoughterOneQy[fCurrentEventDoughterOne] = Qyi;
		  if(tracklistusedinEP.Contains(track))
		    fCurrentEventDoughterOneIn[fCurrentEventDoughterOne] = 1;//dummy
		  else fCurrentEventDoughterOneIn[fCurrentEventDoughterOne] = 0;//dummy
		  fCurrentEventDoughterOne++;
		  fHistdEdxKaonTpc->Fill(pTrackTpc*charge, TpcSignal);
		  if(isTOF)fHistdEdxKaonTof->Fill(pTrack*charge, tofbeta);
		}
	      //}//IsPIDold
	    }//pass track
	  //fHistDCAxyPOI->Fill(1);


	}//track loop
	//	  fHistDCAzPOI->Fill(1);
	
      }//isESD
      
      //std::cout<<" Total particle after track cuts   "<<tracklist->GetEntries()<<"  "<<passtrack<<"  "<<nRP<<std::endl;
      
      fHistPionVsKaonMult->Fill(fCurrentEventDoughterTwo,fCurrentEventDoughterOne);
      //fBufferPointer = (Int_t)(fCurrentEventCentralityBin)*10 + (Int_t)((fCurrentEventVz+10.0)/2.0);
      // std::cout<<centBin<<"   inside ============= eventloop "<<fCurrentEventCentralityBin<<endl;
      //       std::cout<<"x  "<<fBufferPointer<<"    "<<endl;
      //      std::cout<< "y   "<<fCurrentEventDoughterTwo<<endl;
      fCurrentEventCentralityBin=centBin;

      //      cout<" bufffferrrr check  "<< fBufferEventNEvents[fBufferPointer] <<endl;
      if(fBufferPointer < 100 || fBufferPointer >= 0){
	//Same event invariant mass
	MakeRealPair(fQVZeroA, fQVZeroC);
       	MakeLikePair(fQVZeroA, fQVZeroC);
	MakeMixedPair(fBufferPointer,fQVZeroA, fQVZeroC);
	//Copy current event to buffer
       	CopyCurrentToBuffer(fBufferPointer,fQVector);
      }
      //tracklistusedinEP.Clear();
      tracklist->Clear();
      delete tracklist;
      tracklist = 0;
      
    }//passed event
  }
  
  PostData(1, fOutput);
}
  //________________________________________________________________________
void AliAnalysisTaskPhiSAR::Terminate(Option_t *) 
{
  // Draw result to screen, or perform fitting, normalizations
  // Called once at the end of the query
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPhiSAR::PassEvent(AliVEvent *evt){
  //AliVEvent *evt=InputEvent();
  // string to sum messages
  TString msg("");

  //Everything 
  fHistEventCount->Fill(0);
  //Physics selection 
  Bool_t isSelected=0;
  if (evt->InheritsFrom(AliESDEvent::Class())) {
    // type ESD
    // ESD specific check: Physics Selection
    // --> if this is failed, the event is rejected
    //isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB ||AliVEvent::kCentral || AliVEvent::kSemiCentral);
    //isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    //    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    UInt_t  maskisSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
    isSelected = (maskisSelected & AliVEvent::kINT7)  == AliVEvent::kINT7;
    fHistEventCount->Fill(1);
    if (!isSelected) {
      //  x++;
      AliDebugClass(1, "Event does not pass physics selections");
      return kFALSE;
    }
    //    y++;
    fHistEventCount->Fill(2);
  }
  else {
    AliError(Form("Bad input event class: %s", evt->ClassName()));
    return kFALSE;
  }
  //  cout<<x<<"    events    "<<y<<endl;
  fHistEventCount->Fill(3);
  Int_t          nContributors;
  Double_t       zVertex;
  // retrieve ESD event
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(evt);
  if (esd) {
    // pile-up check
    if (fCheckPileUp) 
      {
	if (esd->IsPileupFromSPD()) return kFALSE;
      }
    // get the best primary vertex:
    // first try the one with tracks
    const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
    const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
    const AliESDVertex *vTPC  = esd->GetPrimaryVertexTPC();
    Int_t               ncTrk = -1;
    Int_t               ncSPD = -1;
    Int_t               ncTPC = -1;
    Double_t            vzTrk = 1000000.0;
    Double_t            vzSPD = 1000000.0;
    Double_t            vzTPC = 1000000.0;
    
    //if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
    //if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
    //if (vTPC) vzTPC = TMath::Abs(vTPC->GetZv());
    
    if (vTrk) vzTrk = (vTrk->GetZ());
    if (vSPD) vzSPD = (vSPD->GetZ());
    if (vTPC) vzTPC = (vTPC->GetZ());
    
    if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
    if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
    if (vTPC) ncTPC = (Int_t)vTPC->GetNContributors();
    if (vTrk && ncTrk > 0) {
      nContributors = ncTrk;
      zVertex = vzTrk;
    } else if (vSPD && ncSPD > 0) {
      nContributors = ncSPD;
      zVertex = vzSPD;
    } else if (vTPC && ncTPC > 0) {
      if (!fAcceptTPC)
	return kFALSE;
      else {
	nContributors = ncTPC;
	zVertex = vzTPC;
      }
    } else
      return kFALSE;
  }
  if((TMath::Abs(zVertex) > 10.0 ) || (nContributors < 0)) return kFALSE;
  fHistEventCount->Fill(4);
  
  //Centrality
  Int_t centBin = GetCentrality(evt);
  if(centBin < 0 || centBin > 8) return kFALSE; 
  fHistEventCount->Fill(5);
  
  fCurrentEventVz = zVertex;
  
  return kTRUE;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPhiSAR::GetCentrality(AliVEvent *evt )
{
  //AliVEvent *event = InputEvent();
  evt = InputEvent();
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(evt);

  /*
  AliCentrality *centrality = esd->GetCentrality();
  if (!centrality) {
    AliError("Cannot compute centrality!");
    return -1.0;
  }
  //if(centrality->IsEventInCentralityClass(post(80.,90.,"V0M"))       fCurrentEventCentralityBin = 0;
  if(centrality->IsEventInCentralityClass(70.,80.,"V0M"))       fCurrentEventCentralityBin = 8;
  else if(centrality->IsEventInCentralityClass(60.,70.,"V0M"))  fCurrentEventCentralityBin = 7;
  else if(centrality->IsEventInCentralityClass(50.,60.,"V0M"))  fCurrentEventCentralityBin = 6;
  else if(centrality->IsEventInCentralityClass(40.,50.,"V0M"))  fCurrentEventCentralityBin = 5;
  else if(centrality->IsEventInCentralityClass(30.,40.,"V0M"))  fCurrentEventCentralityBin = 4;
  else if(centrality->IsEventInCentralityClass(20.,30.,"V0M"))  fCurrentEventCentralityBin = 3;
  else if(centrality->IsEventInCentralityClass(10.,20.,"V0M"))  fCurrentEventCentralityBin = 2;
  else if(centrality->IsEventInCentralityClass( 5.,10.,"V0M"))  fCurrentEventCentralityBin = 1;
  else if(centrality->IsEventInCentralityClass( 0., 5.,"V0M"))  fCurrentEventCentralityBin = 0;
  else fCurrentEventCentralityBin = -1;
  */
  Float_t centralitydet = -99.0;
  Float_t centrV0M   = -99.0;

  
  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    return -1;
  }
  else
    {
  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    }
 
  centralitydet = centrV0M;  // This Is Always Default, changes below for other options:
  
  
  ///----> Get 0-10 index of centralitydet:
  Int_t cent10bin = -99, fCurrentEventCentralityBin = -1;

  if(centralitydet<5.0) {
    cent10bin  = 0; 
  }
  else if(centralitydet>=5.0 && centralitydet<10){
    cent10bin  = 1;
  }
  else if(centralitydet>=10.0) {
    cent10bin = abs(centralitydet/10.0)+1;
  }

  fCurrentEventCentralityBin = cent10bin;
  

  
  return fCurrentEventCentralityBin;

  
}


////////// new centrality 
Float_t AliAnalysisTaskPhiSAR::GetCentralityValue(AliVEvent *evt )
{
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(evt);

  Float_t centralitydet = -99.0;
  Float_t centrV0M   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this                                
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    return -1;
  }


 else
   {
     centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
   }

  centralitydet = centrV0M;  // This Is Always Default, changes below for other options:                                                         
 
return centralitydet;


}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPhiSAR::PassTrack(AliVTrack *trackv)
{
  Bool_t  pass=kTRUE;
  AliESDtrack *esdt  = dynamic_cast<AliESDtrack*>(trackv);
  if(!esdt) pass=kFALSE;
  //check cuts on ESD tracks
  Float_t dcaxy = 0.0;
  Float_t dcaz  = 0.0;
  if(esdt){
    if(fTrackType.CompareTo("GLOBAL") || fTrackType.CompareTo("TPC")){
      esdt->GetImpactParameters(dcaxy,dcaz);
      const AliExternalTrackParam* pout = esdt->GetOuterParam();
      const AliExternalTrackParam* pin = esdt->GetInnerParam();
      if (fIgnoreTPCzRange)
	{
	  if(pin && pout)
	    {
	      Double_t zin = pin->GetZ();
	      Double_t zout = pout->GetZ();
	      if (zin*zout<0) pass=kFALSE;   //reject if cross the membrane
	      if (zin < fIgnoreTPCzRangeMin || zin > fIgnoreTPCzRangeMax) pass=kFALSE;
	      if (zout < fIgnoreTPCzRangeMin || zout > fIgnoreTPCzRangeMax) pass=kFALSE;
	    }
	}
      if(fCutMinimalTPCdedx) 
	{
	  if(esdt->GetTPCsignal() < fMinimalTPCdedx) pass=kFALSE;
	}
      //some stuff is still handled by AliESDtrackCuts class - delegate
      if (fAliESDtrackCuts)
	{
	  if (!fAliESDtrackCuts->AcceptTrack(esdt)) pass=kFALSE;
	}
      //      pass =kTRUE;
      }

  }//if esd track

  return pass;
}

//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPhiSAR::MakeRealPair(TVector2 * QvA, TVector2 * QvC)
{
  
  //  std::cout<<"  + "<<endl;
  Double_t EPAngle = QvC->Phi()/2.;
  Double_t weight = 1;
  TLorentzVector DgtTwo(0.,0.,0.,0.);
  TLorentzVector DgtOne(0.,0.,0.,0.);
  TLorentzVector Mother(0.,0.,0.,0.);
  //subhash
  Double_t QvA_X = QvA->X();
  Double_t QvA_Y = QvA->Y();

  Double_t QvC_X = QvC->X();
  Double_t QvC_Y = QvC->Y();
  TVector2 QvA_new, QvC_new;
  AliVEvent *event = InputEvent();
  Float_t centralitydet = GetCentralityValue(event);
  
  Double_t DgtOneCos2phi, DgtOneSin2phi, DgtTwoCos2phi,DgtTwoSin2phi;
  for(Int_t i = 0; i < fCurrentEventDoughterTwo; i++)
    {
      DgtTwo.SetXYZM(fCurrentEventDoughterTwoPx[i],fCurrentEventDoughterTwoPy[i],fCurrentEventDoughterTwoPz[i],pionMass);//sub for K+K-
      //  std::cout<<"  ++ "<<endl;
      for(Int_t j = 0; j < fCurrentEventDoughterOne; j++)
	{
	  DgtOne.SetXYZM(fCurrentEventDoughterOnePx[j],fCurrentEventDoughterOnePy[j],fCurrentEventDoughterOnePz[j],kaonMass); //sub for K+K-
	  //  std::cout<<"  +++ "<<endl;
	  // Mother =(DgtOne)+(DgtTwo)
	  Mother = DgtOne + DgtTwo; 
	  Double_t rapidity = Mother.Rapidity();
	  if( TMath::Abs(rapidity) > fPairRapidity) continue;
	  if(fCurrentEventDoughterOneCharge[j] + fCurrentEventDoughterTwoCharge[i] == 0 )
	    {	      
	      //  std::cout<<"  ++++ "<<endl;
	      Double_t phiAngle = Mother.Phi(); //gives distribution from (-pi to pi)
	      if(phiAngle < 0.) phiAngle += TMath::TwoPi();//gives distribution from (0 to 2pi)
	      
	      QvA_new.Set(QvA_X, QvA_Y);
	      QvC_new.Set(QvC_X, QvC_Y);
	      //std::cout<<"The values ===== "<<Qv_X<<"\t"<< TMath::Cos(2.*DgtOne.Phi())<<"\t"<< TMath::Cos(2.*DgtTwo.Phi()) <<"\t"<<Qv_Y <<"\t"<< TMath::Sin(2.*DgtOne.Phi()) <<"\t"<< TMath::Sin(2.*DgtTwo.Phi())<<std::endl;
	      
	      Double_t EPAngle_new = QvC_new.Phi()/2.;
	      //-------------
	      //phiAngle -= EPAngle;
	      phiAngle -= EPAngle_new;
	      Double_t cosv2 = TMath::Cos(2.0*phiAngle);
	      Double_t sinv2 = TMath::Sin(2.0*phiAngle);
	      Double_t dphi = phiAngle;
	      if(dphi <  0) dphi +=  TMath::Pi();
	      if(dphi > TMath::Pi()) dphi -=  TMath::Pi();
	      // Now dphi distribution is from (0 to pi)

	      Double_t costhetastarA = CosThetaStar(Mother,DgtOne,QvA_new);
	      Double_t costhetastarC = CosThetaStar(Mother,DgtOne,QvC_new);
	      Double_t cosphipsiA = CosPhiPsi(Mother,DgtOne,QvA_new);
	      Double_t cosphipsiC = CosPhiPsi(Mother,DgtOne,QvC_new);
	      
	      Double_t arrA[4]={Mother.Pt(),cosv2,Mother.M(),centralitydet};
	      Double_t arrC[4]={Mother.Pt(),sinv2,Mother.M(),centralitydet};
	      Double_t arrmix[4]={Mother.Pt(),dphi,Mother.M(),centralitydet};

	      if(fCurrentEventCentralityBin >=0 && fCurrentEventCentralityBin < kCenBin) 
		{
		  //if(Mother.Pt()< 0.1 || Mother.Pt()>= 10) continue;
		  // Double_t y[0]=Mother.Pt();
		  // Double_t z[0]=Mother.Pt();
		  //    Double_t arrA[5]={Mother.Pt(),costhetastarA,cosphipsiA,Mother.M(),centralitydet};
		  // Double_t arrC[5]={Mother.Pt(),costhetastarC,cosphipsiC,Mother.M(),centralitydet};
		         
		  // for(Int_t k=0;k<3;k++)
		  // {
		  //		  cout<<Mother.Pt()<<"    "<<cosv2<<"    "<<Mother.M()<<"   "<<centralitydet<<endl;
		  fhInvMassSAEPvzeroA->Fill(arrA,weight);
		  fhInvMassSAEPvzeroC->Fill(arrC,weight);
		  fhInvMassMixSAEPvzeroA->Fill(arrmix,weight);
		      //}
		  // fhInvMassv2[fCurrentEventCentralityBin]->Fill(Mother.Pt(),cosv2,Mother.M(),weight);
		  // fhInvMassSinv2[fCurrentEventCentralityBin]->Fill(Mother.Pt(),sinv2,Mother.M(),weight);
		  // fhNumInvMassvsPtPhi[fCurrentEventCentralityBin]->Fill(Mother.Pt(),dphi,Mother.M(),weight);
		  
		}
	    }//charge conservation loop	  
	  Mother.Clear();
	  DgtOne.Clear();
	}//K+ loop
      DgtTwo.Clear();
    }//K- loop

  return 0;
}

//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPhiSAR::MakeLikePair(TVector2 * QvA, TVector2 * QvC)
{
  
  Double_t EPAngle = QvC->Phi()/2.;
  Double_t weight = 1;
  TLorentzVector DgtTwo(0.,0.,0.,0.);
  TLorentzVector DgtOne(0.,0.,0.,0.);
  TLorentzVector Mother(0.,0.,0.,0.);
  //subhash
  Double_t QvA_X = QvA->X();
  Double_t QvA_Y = QvA->Y();

  Double_t QvC_X = QvC->X();
  Double_t QvC_Y = QvC->Y();
  TVector2 QvA_new, QvC_new;

  AliVEvent *event = InputEvent();
  Float_t centralitydet = GetCentralityValue(event);

  Double_t DgtOneCos2phi, DgtOneSin2phi, DgtTwoCos2phi,DgtTwoSin2phi;
  for(Int_t i = 0; i < fCurrentEventDoughterTwo; i++)
    {
      DgtTwo.SetXYZM(fCurrentEventDoughterTwoPx[i],fCurrentEventDoughterTwoPy[i],fCurrentEventDoughterTwoPz[i],pionMass); //sub for K+K-

      for(Int_t j = i+1; j < fCurrentEventDoughterOne; j++)
	{
	  DgtOne.SetXYZM(fCurrentEventDoughterOnePx[j],fCurrentEventDoughterOnePy[j],fCurrentEventDoughterOnePz[j],kaonMass); //sub for K+K-

	  // Mother =(DgtOne)+(DgtTwo)
	  Mother = DgtOne + DgtTwo; 
	  Double_t rapidity = Mother.Rapidity();
	  if( TMath::Abs(rapidity) > fPairRapidity) continue;
	  if(fCurrentEventDoughterOneCharge[j] + fCurrentEventDoughterTwoCharge[i] == 0 )continue;
	  if(fCurrentEventDoughterOneId[j] == fCurrentEventDoughterTwoId[i])continue;
	  
	  if(fCurrentEventDoughterOneCharge[j] + fCurrentEventDoughterTwoCharge[i] != 0 )
	    {	      
	      Double_t phiAngle = Mother.Phi(); //gives distribution from (-pi to pi)
	      if(phiAngle < 0.) phiAngle += TMath::TwoPi();//gives distribution from (0 to 2pi)
	      
	      QvA_new.Set(QvA_X, QvA_Y);
	      QvC_new.Set(QvC_X, QvC_Y);
	      //std::cout<<"The values ===== "<<Qv_X<<"\t"<< TMath::Cos(2.*DgtOne.Phi())<<"\t"<< TMath::Cos(2.*DgtTwo.Phi()) <<"\t"<<Qv_Y <<"\t"<< TMath::Sin(2.*DgtOne.Phi()) <<"\t"<< TMath::Sin(2.*DgtTwo.Phi())<<std::endl;
	      
	      Double_t EPAngle_new = QvC_new.Phi()/2.;
	      //-------------
	      //phiAngle -= EPAngle;
	      phiAngle -= EPAngle_new;
	      Double_t cosv2 = TMath::Cos(2.0*phiAngle);
	      Double_t sinv2 = TMath::Sin(2.0*phiAngle);
	      Double_t dphi = phiAngle;
	      if(dphi <  0) dphi +=  TMath::Pi();
	      if(dphi > TMath::Pi()) dphi -=  TMath::Pi();
	      // Now dphi distribution is from (0 to pi)

	      Double_t costhetastarA = CosThetaStar(Mother,DgtOne,QvA_new);
	      Double_t costhetastarC = CosThetaStar(Mother,DgtOne,QvC_new);
	      Double_t cosphipsiA = CosPhiPsi(Mother,DgtOne,QvA_new);
	      Double_t cosphipsiC = CosPhiPsi(Mother,DgtOne,QvC_new);
	      
	      Double_t arrA[4]={Mother.Pt(),dphi,Mother.M(),centralitydet};


	      if(fCurrentEventCentralityBin >=0 && fCurrentEventCentralityBin < kCenBin && (fCurrentEventDoughterOneCharge[i] + fCurrentEventDoughterTwoCharge[j] ==2)) 
		{
		  fhInvMassLikePPSAEPvzeroA->Fill(arrA);
		}
	      if(fCurrentEventCentralityBin >=0 && fCurrentEventCentralityBin < kCenBin && (fCurrentEventDoughterOneCharge[i] + fCurrentEventDoughterTwoCharge[j] ==-2)) 
		{
		  fhInvMassLikeMMSAEPvzeroA->Fill(arrA);
		}
	    }//charge conservation loop	  
	  Mother.Clear();
	  DgtOne.Clear();
	}//K+ loop
      DgtTwo.Clear();
    }//K- loop
  
  return 0;
}

//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPhiSAR::MakeMixedPair(Int_t bufferPointer, TVector2 * QvA, TVector2 * QvC)
{
  
  TLorentzVector DgtTwo(0,0,0,0);
  TLorentzVector DgtOne(0,0,0,0);
  TLorentzVector Mother(0,0,0,0);
  Double_t weight = 1;
  Double_t EPAngle = QvC->Phi()/2.;
  //Float_t pi = TMath::Pi();

  Double_t QvA_Xm = QvA->X();
  Double_t QvA_Ym = QvA->Y();
  Double_t QvC_Xm = QvC->X();
  Double_t QvC_Ym = QvC->Y();
  TVector2 QvA_newm, QvC_newm;
  AliVEvent *event = InputEvent();
  Float_t centralitydet = GetCentralityValue(event);

  Double_t DgtOneCos2phi, DgtOneSin2phi, DgtTwoCos2phi,DgtTwoSin2phi;
  Int_t MaxBufferPointer = 100;
  //  cout<<"-----------mixed "<<fBufferEventNEvents[bufferPointer]<<endl;
  //  for(Int_t k = 0; k < fBufferEventNEvents[bufferPointer]; k++)
  for(Int_t k = 0; k <5; k++)
    {
      for(Int_t i = 0; i < fCurrentEventDoughterTwo; i++){
	DgtTwo.SetXYZM(fCurrentEventDoughterTwoPx[i],fCurrentEventDoughterTwoPy[i],fCurrentEventDoughterTwoPz[i],pionMass);//sub for K+K-
	for(Int_t j = 0; j < fBufferEventDoughterOne[k*MaxBufferPointer+bufferPointer]; j++){
	  if(fCurrentEventDoughterTwoCharge[i] + fBufferEventDoughterOneCharge[k*MaxBufferPointer+bufferPointer][j] != 0) continue;
	  DgtOne.SetXYZM(fBufferEventDoughterOnePx[k*MaxBufferPointer+bufferPointer][j],  fBufferEventDoughterOnePy[k*MaxBufferPointer+bufferPointer][j], fBufferEventDoughterOnePz[k*MaxBufferPointer+bufferPointer][j], kaonMass);//sub for K+K-

	  Mother = DgtTwo + DgtOne;
	  
	  Double_t rapidity = Mother.Rapidity();
	  if( TMath::Abs(rapidity) > fPairRapidity) continue;
	  
	  Double_t phiAngle = Mother.Phi();//gives distribution from (-pi to pi)
	  if(phiAngle < 0.) phiAngle += TMath::TwoPi();//gives distribution from (0 to 2pi)
	  
	  QvA_newm.Set(QvA_Xm, QvA_Ym);
	  QvC_newm.Set(QvC_Xm, QvC_Ym);

	  Double_t EPAngle_newm = QvC_newm.Phi()/2.;
	  //-------------
	  //phiAngle -= EPAngle;
	  phiAngle -= EPAngle_newm;
	  Double_t dphi = phiAngle;//gives distribution from (-pi to 2pi)
	  if(dphi <  0) dphi +=  TMath::Pi();
	  if(dphi > TMath::Pi()) dphi -=  TMath::Pi();
	  // now dphi is from 0 to pi

	  Double_t costhetastarA = CosThetaStar(Mother,DgtTwo,QvA_newm);
	  Double_t costhetastarC = CosThetaStar(Mother,DgtTwo,QvC_newm);
	  Double_t cosphipsiA = CosPhiPsi(Mother,DgtTwo,QvA_newm);
	  Double_t cosphipsiC = CosPhiPsi(Mother,DgtTwo,QvC_newm);

	  Double_t arrA[4]={Mother.Pt(),dphi,Mother.M(),centralitydet};

	  if(fCurrentEventCentralityBin >=0 && fCurrentEventCentralityBin < kCenBin) 
	    {  
	      fhInvMassMixSAEPvzeroC->Fill(arrA);
		
	    }
	  
	}  //buffer event DgtOne loop
      }   //current event DgtTwo loop

      Mother.Clear();
      //QvA_newm.Clear();
      //QvC_newm.Clear();
      /*
      for(Int_t i = 0; i < fBufferEventDoughterTwo[k*MaxBufferPointer+bufferPointer]; i++)
	{
	  
	  DgtTwo.SetXYZM(fBufferEventDoughterTwoPx[k*MaxBufferPointer+bufferPointer][i],fBufferEventDoughterTwoPy[k*MaxBufferPointer+bufferPointer][i],fBufferEventDoughterTwoPz[k*MaxBufferPointer+bufferPointer][i],pionMass);//sub for pi- now     

	  for(Int_t j = 0; j < fCurrentEventDoughterOne; j++)
	    {
	      if(fBufferEventDoughterTwoCharge[k*MaxBufferPointer+bufferPointer][i] + fCurrentEventDoughterOneCharge[j] != 0) continue;
	      DgtOne.SetXYZM(fCurrentEventDoughterOnePx[j],fCurrentEventDoughterOnePy[j],fCurrentEventDoughterOnePz[j],kaonMass);

	      Mother = DgtTwo + DgtOne;

	      Double_t rapidity = Mother.Rapidity();
	      if( TMath::Abs(rapidity) > fPairRapidity) continue;
	      
	      Double_t phiAngle = Mother.Phi(); //gives distribution from (-pi to pi)
	      if(phiAngle < 0.) phiAngle += TMath::TwoPi();//gives distribution from (0 to 2pi)

	      //QvA_newm.Set(QvA_Xm, QvA_Ym);
	      //QvC_newm.Set(QvC_Xm, QvC_Ym);
	      
	      Double_t EPAngle_newm = QvC_newm.Phi()/2.;
	      phiAngle -= EPAngle_newm;
	      //phiAngle -= EPAngle;	 
	      Double_t dphi = phiAngle;
	      if(dphi <  0) dphi +=  TMath::Pi();
	      if(dphi > TMath::Pi()) dphi -=  TMath::Pi();	      
	      //now dphi is from 0 to pi

	      Double_t costhetastarA = CosThetaStar(Mother,DgtOne,QvA_newm);
	      Double_t costhetastarC = CosThetaStar(Mother,DgtOne,QvC_newm);
	      
	      if(fCurrentEventCentralityBin >= 0 && fCurrentEventCentralityBin < kCenBin)
		{
		  //if(Mother.Pt() < 0.1 || Mother.Pt()>= 10.) continue;
		  fhInvMassMixSAEPvzeroA[fCurrentEventCentralityBin]->Fill(Mother.Pt(), costhetastarA, Mother.M(),weight);
		  fhInvMassMixSAEPvzeroC[fCurrentEventCentralityBin]->Fill(Mother.Pt(), costhetastarC, Mother.M(),weight);
		}
	    }  //current event DgtOne loop
	}  //buffer event DgtTwo loop
      Mother.Clear();
      */
    }
  
  return 0;
}

//--------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------
void AliAnalysisTaskPhiSAR::CopyCurrentToBuffer(Int_t bufferPointer,TVector2 *Qv)
{
  Double_t EPAngle = Qv->Phi()/2.;
  Int_t Nmix = 5;
  Int_t MaxBufferPointer = 100;
  Int_t iran = 0;
  if(fBufferEventNEvents[bufferPointer] >= Nmix) fBufferEventFull[bufferPointer] = 1;
  //  cout<<"  bffffer  "<<fBufferEventNEvents[bufferPointer]<<endl;
  TRandom3 *Rand = new TRandom3(iran++);
  Int_t eventPointer = -1;
  if(fBufferEventFull[bufferPointer]) 
    { 
      // full - random rewrite one
      do 
	{ Double_t rrr = Rand->Rndm();
	  eventPointer = (Int_t)(Nmix*(1.0 - rrr));
	  //std::cout<<"***********************value of event pointer***********************:"<<eventPointer<<std::endl;
	} while(eventPointer < 0 || eventPointer >= Nmix);
      //  std::cout<<" do loop   "<<eventPointer<<endl;    
} 
  else 
    { // not full
      //      eventPointer = fBufferEventNEvents[bufferPointer];
      eventPointer = 1;
      //      eventPointer      
      //      std::cout<<" outside randam   "<<eventPointer<<endl;
    }
  delete gRandom;
  //  std::cout<<"xy   "<<fCurrentEventDoughterTwo<<endl;
  //  fBufferEventPsi[bufferPointer+MaxBufferPointer*eventPointer] = EPAngle;
  //  std::cout<<bufferPointer<<"   "<<MaxBufferPointer<<" bffffffffff   "<<eventPointer<<"    "<< fCurrentEventDoughterTwo << "   "<<fCurrentEventDoughterOne<<endl;
  fBufferEventDoughterTwo[bufferPointer+MaxBufferPointer*eventPointer]= fCurrentEventDoughterTwo;
  for(Int_t i = 0; i < fCurrentEventDoughterTwo; i++){
    fBufferEventDoughterTwoCharge[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoCharge[i];
    fBufferEventDoughterTwoPx[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoPx[i];
    fBufferEventDoughterTwoPy[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoPy[i];
    fBufferEventDoughterTwoPz[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoPz[i];
    fBufferEventDoughterTwoId[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoId[i];
    fBufferEventDoughterTwoIn[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterTwoIn[i];
  }
  fBufferEventDoughterOne[bufferPointer+MaxBufferPointer*eventPointer]=fCurrentEventDoughterOne;
  for(Int_t i=0; i< fCurrentEventDoughterOne; i++){
    fBufferEventDoughterOneCharge[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOneCharge[i];
    fBufferEventDoughterOnePx[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOnePx[i];
    fBufferEventDoughterOnePy[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOnePy[i];
    fBufferEventDoughterOnePz[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOnePz[i];
    fBufferEventDoughterOneId[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOneId[i];
    fBufferEventDoughterOneIn[bufferPointer+MaxBufferPointer*eventPointer][i]=fCurrentEventDoughterOneIn[i];
  }
  if(fBufferEventNEvents[bufferPointer] < Nmix) 
    {
      fBufferEventNEvents[bufferPointer] += 1; 
    }

}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskPhiSAR::GetDipAngle(Double_t aX,Double_t aY,Double_t aZ,Double_t bX,Double_t bY,Double_t bZ)
{ 
  Double_t p,adotb;
  adotb = ((TMath::Sqrt(aX*aX+aY*aY))*(TMath::Sqrt(bX*bX+bY*bY))) + aZ*bZ;
  p     = (TMath::Sqrt(aX*aX+aY*aY+aZ*aZ))*(TMath::Sqrt(bX*bX+bY*bY+bZ*bZ));
  Double_t theta;
  if(p != 0)
    theta = TMath::ACos(adotb/p);  
  else theta = -9999.;

  return theta;
}

//__________________________________________________________________________________________________

Double_t AliAnalysisTaskPhiSAR::CosThetaStar(TLorentzVector mother, TLorentzVector daughter0,TVector2& Qvect)
{
  //
  // Return cosine of angle of one daughter to the resonance momentum in its rest frame
  //
  //
  TVector3 EPVect;
  EPVect.SetXYZ(Qvect.X(),Qvect.Y(),0.);
  TVector3 BeamVect;
  BeamVect.SetXYZ(0,0,1);

  TVector3 normal = EPVect.Cross(BeamVect);

  TVector3 UnitNormal = normal.Unit();

  //std::cout<<"Magnitude of normal vector  ==  "<<normal.Mag()<<"\t"<<UnitNormal.Mag()<<"\t"<<Qvect.Mod()<<"\t"<<EPVect.Mag()<<"\t"<<BeamVect.Mag()<<std::endl;
  
  //TVector3 momentumM(mother.Vect());
  //TVector3 normal(mother.Y() / momentumM.Mag(), -mother.X() / momentumM.Mag(), 0.0);
  //TVector3 normal(mother.Y() / mother.Pt(), -mother.X() / mother.Pt(), 0.0);
  
  // Computes components of beta
  Double_t betaX = -mother.X() / mother.E();
  Double_t betaY = -mother.Y() / mother.E();
  Double_t betaZ = -mother.Z() / mother.E();
  
  // Computes Lorentz transformation of the momentum of the first daughter
  // into the rest frame of the mother and theta*
  daughter0.Boost(betaX, betaY, betaZ);
  TVector3 momentumD = daughter0.Vect();
  
  Double_t cosThetaStar = TMath::Abs(UnitNormal.Dot(momentumD) / momentumD.Mag());
  
  return cosThetaStar;
}
//----------------------------------------------------------------------------------------------------------------



Double_t AliAnalysisTaskPhiSAR::CosPhiPsi(TLorentzVector mother, TLorentzVector daughter0,TVector2& Qvect)
{
  //
  // Return cosine of angle of one daughter to the resonance momentum in its rest frame
  //
  //
  TVector3 EPVect;
  EPVect.SetXYZ(Qvect.X(),Qvect.Y(),0.);
  TVector3 BeamVect;
  BeamVect.SetXYZ(0,0,1);

  Double_t psi=0.5*TMath::ATan2(Qvect.Y(),Qvect.X());
  if (psi < 0.) psi += TMath::Pi();
  if (psi < TMath::Pi()) psi -= TMath::Pi();


  Double_t phiAngle = mother.Phi();//gives distribution from (-pi to pi)                                                               
  if(phiAngle < 0.) phiAngle += TMath::TwoPi();//gives distribution from (0 to 2pi)                                                    
 
  
  Double_t phipsi = phiAngle-psi;
  
  return phipsi;
}
//----------------------------------------------------------------------------------------------------------------



//________________________________________________________________________
void AliAnalysisTaskPhiSAR::SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts){
  if(fAliESDtrackCuts){ 
    delete fAliESDtrackCuts;
    fAliESDtrackCuts = 0;
  }
  fUsercuts = kTRUE;
  fAliESDtrackCuts = trackcuts;
}
//________________________________________________________________________
void AliAnalysisTaskPhiSAR::SetPOIAndRPTrackType(TString tracktype, TString rptracktype){
  fTrackType = tracktype;
  fRPTrackType = rptracktype;

  std::cout<<"Setting track type for POI and RP : "<<fTrackType.Data()<<"  "<<fRPTrackType.Data()<<std::endl;
  /*
  if (!fUsercuts) {
    if (fTrackType.CompareTo("GLOBAL")==0){ 
      fAliESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
      std::cout<<"GetStandardITSTPCTrackCuts2010  Track cuts will be used to select the track =============== "<<std::endl;
    }	
    if (fTrackType.CompareTo("TPC")==0){  
      fAliESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    }
    fAliESDtrackCuts->SetPtRange(0.15,20.);
    fAliESDtrackCuts->SetEtaRange(-0.8,0.8);
  }
  */
}
//--------------------------------------------------------------------
Bool_t AliAnalysisTaskPhiSAR::IsSelectedPION(AliVTrack *track)
{
  //
  // Check flags
  //

  
  if(!GetStatus(track)) return kFALSE;
  
  //this part from ALICE_ROOT/PWGLF/RESONANCE/AliRsnCutDaughterKStar2010.cxx
  Bool_t   isTOF  = MatchTOF(track);
  Double_t pTPC   = track->GetTPCmomentum();
  Double_t p      = track->P();
  Double_t nsTPC  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
  Double_t nsTOF  = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
  Double_t maxTPC = 1E20;
  Double_t maxTOF = 1E20;
  Bool_t accept   = kFALSE;
  //  std::cout<<pTPC<<"   "<<p<<endl;
  // applies the cut differently depending on the PID and the momentum
  if (isTOF)
    {
      //AliDebugClass(21, Form("Checking Default PID , TOF is Present : nsigma = %f", nsTOF));
      // TPC: 3sigma cut for all
      if (nsTPC > 3.0) return kFALSE;
      // TOF: 3sigma below 1.5 GeV, 2sigma above
      //if (p < 1.5) 
      maxTOF = 3.0; 
      //else maxTOF = 2.0;
      //     if (nsTOF <= maxTOF) return kTRUE;
      
      accept = (nsTOF <= maxTOF);      
//std::cout << " Inside (isTOF)  nsTPC    " <<nsTPC << " nsTOF    "<<nsTOF<< "       accecptTOF "<<accept<<std::endl;

    } 
  else 
    {
      // TPC: 
      // below 350 MeV: 5sigma
      // between 350 and 500 MeV: 3sigma
      // pions above 500 MeV: 2sigma
      //AliDebugClass(21, Form("Checking Default PID , No TOF Only, TPC Present : nsigma = %f", nsTPC));
      if (pTPC <= 0.35) 
	maxTPC = 5.0;
      else if (pTPC <= 0.5)
	maxTPC = 3.0;
      //else if (pTPC > 0.5 &&  AliPID::kPion == 2)
      else if (pTPC > 0.5)
	maxTPC = 2.0;
      //      if (nsTPC <= maxTPC) return kTRUE;
      accept = (nsTPC <= maxTPC) ;

    
  //      std::cout << " Inside (isTPC)  nsTPC    " <<pTPC<<"   "<<nsTPC << " nsTOF    "<<nsTOF<< "       accecptTPC "<<accept<<std::endl; 

    }//
  //  std::cout << " Inside "<<accept<<endl;
  return accept;


  //  return kFALSE;
}
//______________________________________________________________________________

//-----------------------------------------------------
Bool_t AliAnalysisTaskPhiSAR::IsSelectedKAON(AliVTrack *track)
{
  //
  // Check flags
  //
 
  if(!GetStatus(track)) return kFALSE;
  //this part from ALICE_ROOT/PWGLF/RESONANCE/AliRsnCutDaughterKStar2010.cxx
  Bool_t   isTOF  = MatchTOF(track);
  Double_t pTPC   = track->GetTPCmomentum();
  Double_t p      = track->P();
  Double_t nsTPC  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
  Double_t nsTOF  = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
  Double_t maxTPC = 1E20;
  Double_t maxTOF = 1E20;
  Bool_t accept   = kFALSE;
  //  cout<<nsTOF
  // applies the cut differently depending on the PID and the momentum
  if (isTOF)
    {
      //AliDebugClass(21, Form("Checking Default PID , TOF is Present : nsigma = %f", nsTOF));
      
      // TPC: 3sigma cut for all
      if (nsTPC > 2.0) return kFALSE;
      // TOF: 3sigma below 1.5 GeV, 2sigma above
      //if (p < 1.5) 
      maxTOF = 3.0;
      //else maxTOF = 2.0;
      //      if (nsTOF <= maxTOF)  return kTRUE;
     accept = (nsTOF <= maxTOF);

    } 
  else 
    {
      // TPC: 
      // below 350 MeV: 5sigma
      // between 350 and 500 MeV: 3sigma
      // pions above 500 MeV: 2sigma
      // kaons above 500 MeV: 2sigma
      //AliDebugClass(21, Form("Checking Default PID , No TOF Only, TPC Present : nsigma = %f", nsTPC));
      if (pTPC <= 0.35) 
	maxTPC = 5.0;
      else if (pTPC <= 0.5)
	maxTPC = 3.0;
      //else if (pTPC > 0.5 && AliPID::kKaon == 3)
      else if (pTPC > 0.5)
      maxTPC = 2.0;
      //      if (nsTPC <= maxTPC)  return kTRUE;
      accept = (nsTPC <= maxTPC);
      //      if (nsTPC <= maxTPC)  return kTRUE;
      //std::cout << " Inside (isTPC)  nsTPC    " <<pTPC<<"   "<<nsTPC << " nsTOF    "<<nsTOF<< "       accecptTPC "<<accept<<std::endl; 
    }
 
  //  accept=kTRUE;
  return accept; 
  //  return kTRUE; 

}
//______________________________________________________________________________

Double_t AliAnalysisTaskPhiSAR::GetTOFBeta(AliVTrack *track)
{
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(track);
  if(!esdtrack) return 0.;
  const Double_t c = 2.99792457999999984e-02;  
  Double_t p = esdtrack->GetP();
  Double_t l = esdtrack->GetIntegratedLength();  
  Double_t trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t timeTOF = esdtrack->GetTOFsignal()- trackT0; 
  return l/timeTOF/c;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskPhiSAR::GetStatus(const AliVTrack *vtrack)
{
  //
  // Checks if the track has matched the TPC and ITS detector
  //
  if (!vtrack) {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  if ((vtrack->GetStatus() & AliESDtrack::kTPCin   ) == 0) return kFALSE;
  if ((vtrack->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
  if ((vtrack->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
  
  return kTRUE;
}
//______________________________________________________________________________
//Bool_t AliAnalysisTaskPhiSAR::MatchTOF(const AliVTrack *vtrack)
Bool_t AliAnalysisTaskPhiSAR::MatchTOF(AliVTrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//
  
  if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      //      std::cout<<"1"<<endl;
      return kFALSE;
   }
  //  cout<<vtrack->GetStatus()<<" tof  "<<AliESDtrack::kTOFout<<"  "<<AliESDtrack::kTOFpid<<"  "<<AliESDtrack::kTIME<<endl;
    if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
  //  if (!(vtrack->GetStatus() & AliESDtrack::kTOFpid)) return kFALSE;
    if (!(vtrack->GetStatus() & AliESDtrack::kTIME  )) return kFALSE;
    Float_t probmis=fPIDResponse->GetTOFMismatchProbability(vtrack);  
    if(probmis>0.01) return kFALSE;
    //AliESDtrack *esdtrack;//  = dynamic_cast<AliESDtrack*>(vtrack);
    AliESDtrack *esdtrack = dynamic_cast < AliESDtrack * > (vtrack);
    Double_t l=esdtrack->GetIntegratedLength();
    if(l<350) return kFALSE;
//  if(vtrack->IsOn(AliESDtrack::kTOFout) && vtrack->IsOn(AliESDtrack::kTOFpid) && vtrack->IsOn(AliESDtrack::kTIME) && !vtrack->IsOn(AliESDtrack::kTOFmismatch))  return kTRUE;
 
  return kTRUE;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskPhiSAR::GetNSigmaCut(Double_t pTrack)
{
  Double_t nSigCut = 0;
  if(pTrack <  0.35) nSigCut = 3.0;
  if(pTrack >= 0.35 && pTrack < 0.6) nSigCut = 2.0;
  if(pTrack >= 0.6) nSigCut = 1.0;
  return nSigCut;
}
//______________________________________________________________________________
void AliAnalysisTaskPhiSAR::ReSet()
{
  for(Int_t i = 0; i < kDim; i++){
    
    fCurrentEventDoughterOneCharge[i] = 0;
    fCurrentEventDoughterOneId[i] = 0;
    fCurrentEventDoughterOneIn[i] = 0;
    fCurrentEventDoughterOnePx[i] = 0;
    fCurrentEventDoughterOnePy[i] = 0;
    fCurrentEventDoughterOnePz[i] = 0;
    fCurrentEventDoughterOneQx[i] = 0;
    fCurrentEventDoughterOneQy[i] = 0;
    
    fCurrentEventDoughterTwoCharge[i] = 0;
    fCurrentEventDoughterTwoIn[i] = 0;
    fCurrentEventDoughterTwoId[i] = 0;
    fCurrentEventDoughterTwoPx[i] = 0;
    fCurrentEventDoughterTwoPy[i] = 0;
    fCurrentEventDoughterTwoPz[i] = 0;
    fCurrentEventDoughterTwoQx[i] = 0;
    fCurrentEventDoughterTwoQy[i] = 0;
  }
  fCurrentEventCentralityBin = -1;
  fCurrentEventDoughterOne   = 0;
  fCurrentEventDoughterTwo  = 0;
  fCurrentEventVx      = 0;
  fCurrentEventVy      = 0;
  fCurrentEventVz      = 0;
}
//======================================================================


